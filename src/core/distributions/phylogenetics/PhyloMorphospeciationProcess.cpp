#include "PhyloMorphospeciationProcess.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>

#include "RbBoolean.h"
#include "DiscreteTaxonData.h"
#include "DistributionPoisson.h"
#include "FossilizedBirthDeathSpeciationProcess.h"
#include "HomologousDiscreteCharacterData.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RbMathLogic.h"
#include "RlUserInterface.h"
#include "StandardState.h"
#include "StochasticNode.h"
#include "StringUtilities.h"
#include "TopologyNode.h"
#include "TreeDistributionProperties.h"

using namespace RevBayesCore;


PhyloMorphospeciationProcess::PhyloMorphospeciationProcess(const TypedDagNode<Tree> *t, const std::vector<size_t> &k, const TypedDagNode<double> *eps, const TypedDagNode<double> *eps_a, Coding cod, const TypedDagNode<double> *ret, size_t shared) :
    TypedDistribution<AbstractHomologousDiscreteCharacterData>( new HomologousDiscreteCharacterData<StandardState>() ),
    tree( t ),
    epsilon( eps ),
    epsilon_a( eps_a ),
    process( NULL ),
    num_states( k ),
    states_given( k.empty() == false ),
    shared_states( shared ),
    coding( cod ),
    retention( ret )
{
    resolveProcess();

    if ( epsilon_a != NULL && epsilon == NULL )
    {
        throw RbException("dnPhyloMorphospeciation: epsilon_a needs epsilon.");
    }

    // with one change per speciation the characters share their speciations, so dropping some characters
    // is not a condition on each character alone
    if ( coding != ALL && epsilon == NULL )
    {
        throw RbException("dnPhyloMorphospeciation: a coding other than \"all\" needs epsilon.");
    }
    if ( retention != NULL && coding == VARIABLE )
    {
        throw RbException("dnPhyloMorphospeciation: retention needs coding=\"informative\". \"variable\" keeps every autapomorphy, as retention 1 would.");
    }
    if ( retention != NULL && coding == ALL )
    {
        throw RbException("dnPhyloMorphospeciation: retention needs coding=\"informative\", since under \"all\" no character was dropped.");
    }

    addParameter( tree );
    addParameter( epsilon );
    addParameter( epsilon_a );
    addParameter( retention );
    for (const DagNode *n : process->getRateParameters()) addParameter( n );

    for (size_t K : num_states)
    {
        if ( K < 2 ) throw RbException("dnPhyloMorphospeciation: every character needs at least two states.");
    }
    if ( shared_states == 1 ) throw RbException("dnPhyloMorphospeciation: every character needs at least two states.");

    changed = std::vector<size_t>( tree->getValue().getNumberOfNodes(), 0 );

    if ( states_given ) redrawValue();
}


PhyloMorphospeciationProcess* PhyloMorphospeciationProcess::clone( void ) const
{
    return new PhyloMorphospeciationProcess( *this );
}


void PhyloMorphospeciationProcess::resolveProcess( void )
{
    const StochasticNode<Tree> *n = dynamic_cast<const StochasticNode<Tree> *>( tree );
    process = ( n == NULL ) ? NULL : const_cast<FossilizedBirthDeathSpeciationProcess *>(
                  dynamic_cast<const FossilizedBirthDeathSpeciationProcess *>( &n->getDistribution() ) );
    if ( process == NULL )
    {
        throw RbException("dnPhyloMorphospeciation requires a tree drawn from dnFBDSP.");
    }
}


/** The species whose lineage runs through this node: follow the continuing children to a tip. */
long PhyloMorphospeciationProcess::speciesOf( const TopologyNode &node ) const
{
    const TopologyNode *n = &node;
    while ( n->isTip() == false )
    {
        const TopologyNode *next = NULL;
        for (size_t c = 0; c < n->getNumberOfChildren(); c++)
        {
            if ( n->getChild(c).continuesParentSpecies() ) next = &n->getChild(c);
        }
        if ( next == NULL ) return -1;
        n = next;
    }
    return long( n->getIndex() );
}


/** The species of each node, from its continuing child. -1 unless every node below keeps the parent
    species on exactly one of two children, since anything else has probability zero under the model. */
long PhyloMorphospeciationProcess::assignSpecies( const TopologyNode &node )
{
    size_t idx = node.getIndex();
    if ( node.isTip() ) return species[idx] = long( idx );
    if ( node.getNumberOfChildren() != 2 ) return -1;

    size_t continuing = 0;
    for (size_t c = 0; c < 2; c++)
    {
        if ( assignSpecies( node.getChild(c) ) < 0 ) return -1;
        if ( node.getChild(c).continuesParentSpecies() )
        {
            continuing++;
            kept_child[idx] = c;
        }
    }
    if ( continuing != 1 ) return -1;
    budded[idx] = species[ node.getChild( 1 - kept_child[idx] ).getIndex() ];
    return species[idx] = species[ node.getChild( kept_child[idx] ).getIndex() ];
}


/** Each character's chance to change at a budding and at an anagenesis: epsilon, and epsilon_a if given. */
bool PhyloMorphospeciationProcess::changeProbabilities( double &e_budding, double &e_anagenesis ) const
{
    e_budding = e_anagenesis = epsilon->getValue();
    if ( epsilon_a != NULL ) e_anagenesis = epsilon_a->getValue();
    return e_budding >= 0.0 && e_budding <= 1.0 && e_anagenesis >= 0.0 && e_anagenesis <= 1.0;
}


/** int lambda p dt and int lambda_a dt over the part of each node's branch above its species' tau_1, by node index. */
void PhyloMorphospeciationProcess::hiddenIntegrals( const std::vector<long> &sp, const std::vector<double> &tau, std::vector<double> &budding, std::vector<double> &anagenetic ) const
{
    const Tree &t = tree->getValue();
    budding.assign( t.getNumberOfNodes(), 0.0 );
    anagenetic.assign( t.getNumberOfNodes(), 0.0 );
    std::vector<size_t> at;
    std::vector<std::pair<double,double> > spans;
    for (const TopologyNode *node : t.getNodes())
    {
        if ( node->isRoot() ) continue;
        long s = sp[ node->getIndex() ];
        if ( s < 0 ) continue;

        // the species' sampled range carries no change, so only the part above tau_1 counts
        double lo = std::max( node->getAge(), tau[s] );
        double hi = node->getParent().getAge();
        if ( hi > lo )
        {
            at.push_back( node->getIndex() );
            spans.push_back( std::make_pair( lo, hi ) );
        }
    }
    std::vector<double> b = process->hiddenSpeciationIntegrals( spans );
    std::vector<double> a = process->anageneticIntegrals( spans );
    for (size_t k = 0; k < at.size(); k++)
    {
        budding[ at[k] ]    = b[k] - a[k];
        anagenetic[ at[k] ] = a[k];
    }
}


/** tau_1 of every taxon, the species of every node, the Markov chain length of every branch, and with
    epsilon the Poisson means of its hidden events. False unless every node keeps the parent species on
    exactly one child, since anything else has probability zero under the model. */
bool PhyloMorphospeciationProcess::prepare( void )
{
    if ( process->hasSymmetricSpeciation() )
    {
        throw RbException("dnPhyloMorphospeciation: symmetric speciation is not part of the model.");
    }

    RbVector<double> f;
    process->executeMethod( "getFirstAppearances", std::vector<const DagNode*>(), f );
    first.assign( f.size(), 0.0 );
    for (size_t i = 0; i < f.size(); i++) first[i] = f[i];

    const Tree &t = tree->getValue();
    size_t nodes = t.getNumberOfNodes();
    species.assign( nodes, -1 );
    kept_child.assign( nodes, 0 );
    budded.assign( nodes, -1 );
    if ( assignSpecies( t.getRoot() ) < 0 ) return false;
    if ( changed.size() != nodes )
    {
        changed.assign( nodes, 0 );
        latents_ready = false;
    }

    // in an extended tree, a parent ending as a sampled ancestor at the node was replaced by anagenesis
    anagenetic.assign( nodes, 0 );
    bool extended = static_cast<const TreeDistributionProperties *>( process )->isExtended();
    for (const TopologyNode *node : t.getNodes())
    {
        size_t idx = node->getIndex();
        if ( extended && node->isTip() == false && node->getChild( kept_child[idx] ).isSampledAncestorTip() ) anagenetic[idx] = 1;
    }

    double n = double( num_states.size() );
    hiddenIntegrals( species, first, mean_events[BUDDING], mean_events[ANAGENESIS] );
    stem.assign( nodes, 0.0 );
    for (size_t i = 0; i < nodes; i++)
    {
        if ( n > 0 ) stem[i] = ( mean_events[BUDDING][i] + mean_events[ANAGENESIS][i] ) / n;
    }
    return true;
}


/** A start from which every morphology is reachable: one hidden event of each kind wherever it can happen. */
void PhyloMorphospeciationProcess::initializeLatents( void )
{
    const Tree &t = tree->getValue();
    size_t nodes = t.getNumberOfNodes();

    for (int m = 0; m < 2; m++)
    {
        events[m].assign( nodes, 0 );
        for (size_t i = 0; i < nodes; i++)
        {
            double mean = mean_events[m][i];
            if ( t.getNode(i).isRoot() == false && mean > 0.0 ) events[m][i] = std::max( 1L, long( std::llround( mean ) ) );
        }
    }
    latents_ready = true;
}


/** Whether some resolution of a character's ambiguous scores varies, and whether one is parsimony informative. */
static void possiblePatterns( const std::vector<RbBitSet> &scored, size_t K, bool &varies, bool &informative )
{
    std::vector<size_t> counts( K, 0 );
    std::vector<const RbBitSet*> open;
    double resolutions = 1.0;
    for (const RbBitSet &b : scored)
    {
        if ( b.count() == 1 ) counts[ b.find_first() ]++;
        else
        {
            open.push_back( &b );
            resolutions *= double( b.count() );
        }
    }
    varies = informative = false;

    // too many to try: keep the character
    if ( resolutions > 1E5 )
    {
        varies = informative = true;
        return;
    }
    std::function<void(size_t)> walk = [&]( size_t i ) {
        if ( varies && informative ) return;
        if ( i == open.size() )
        {
            size_t present = 0, repeated = 0;
            for (size_t x : counts)
            {
                present += ( x > 0 );
                repeated += ( x > 1 );
            }
            varies = varies || present > 1;
            informative = informative || repeated > 1;
            return;
        }
        for (size_t k = 0; k < K; k++)
        {
            if ( open[i]->test(k) == false ) continue;
            counts[k]++;
            walk( i + 1 );
            counts[k]--;
        }
    };
    walk( 0 );
}


/** The characters the likelihood uses: all but those the coding gives probability zero, with a warning. */
void PhyloMorphospeciationProcess::selectCharacters( void )
{
    const Tree &t = tree->getValue();
    chars.clear();
    autapomorphic.clear();
    size_t dropped = 0;
    for (size_t c = 0; c < num_states.size(); c++)
    {
        std::vector<RbBitSet> scored;
        for (size_t i = 0; i < t.getNumberOfTips(); i++)
        {
            const TopologyNode &tip = t.getTipNode( i );
            const DiscreteCharacterState &x = this->value->getTaxonData( tip.getName() ).getCharacter( c );
            if ( x.isGapState() || x.isMissingState() ) continue;
            RbBitSet bits = x.getState();
            for (size_t k = num_states[c]; k < bits.size(); k++)
            {
                if ( bits.test(k) ) throw RbException() << "dnPhyloMorphospeciation: " << tip.getName() << " has state " << k
                                                        << " at character " << c + 1 << ", which has " << num_states[c] << " states.";
            }
            bits.resize( num_states[c] );
            scored.push_back( bits );
        }

        bool varies = true, informative = true;
        if ( coding != ALL ) possiblePatterns( scored, num_states[c], varies, informative );
        if ( varies == false || ( coding == INFORMATIVE && retention == NULL && informative == false ) )
        {
            dropped++;
            continue;
        }
        chars.push_back( c );
        autapomorphic.push_back( retention != NULL && informative == false );
    }

    if ( dropped > 0 && warned == false )
    {
        static const char *names[] = { "all", "variable", "informative" };
        RBOUT( "WARNING: " + StringUtilities::toString( dropped ) + ( dropped == 1 ? " character is" : " characters are" ) + " impossible under coding=\"" + names[coding] + "\" and will be excluded." );
        warned = true;
    }
}


/** The observed states of every tip, laid out flat, rebuilt only when the data change. */
void PhyloMorphospeciationProcess::buildTipPartials( void )
{
    selectCharacters();
    size_t N = chars.size();
    const Tree &t = tree->getValue();
    size_t T = t.getNumberOfTips();

    columns.assign( N, 0 );
    offset.assign( N, 0 );
    total_states = 0;
    for (size_t c = 0; c < N; c++)
    {
        columns[c] = num_states[ chars[c] ];
        offset[c] = total_states;
        total_states += columns[c];
    }
    max_states = 2;
    for (size_t K : columns) max_states = std::max( max_states, K );

    tip_partials.assign( T, std::vector<double>( total_states, 1.0 ) );
    tip_coded.assign( T, 0 );
    std::vector< std::vector<char> > scored( N, std::vector<char>( T, 0 ) );
    for (size_t i = 0; i < T; i++)
    {
        const TopologyNode &tip = t.getTipNode( i );
        const AbstractDiscreteTaxonData &td = this->value->getTaxonData( tip.getName() );
        std::vector<double> &v = tip_partials[ tip.getIndex() ];
        for (size_t c = 0; c < N; c++)
        {
            const DiscreteCharacterState &x = td.getCharacter( chars[c] );
            if ( x.isGapState() || x.isMissingState() ) continue;
            tip_coded[ tip.getIndex() ] = 1;
            scored[c][ tip.getIndex() ] = 1;

            RbBitSet bits = x.getState();
            for (size_t k = 0; k < columns[c]; k++) v[ offset[c] + k ] = ( k < bits.size() && bits.test(k) ) ? 1.0 : 0.0;
        }
    }

    // characters with the same K scored on the same tips share their probability of being dropped
    group.assign( N, 0 );
    group_states.clear();
    group_scored.clear();
    group_top.clear();
    group_offset.clear();
    group_total = 0;
    std::vector<size_t> first_member;
    if ( coding != ALL )
    {
        std::map<std::string, size_t> seen;
        for (size_t c = 0; c < N; c++)
        {
            std::string key = StringUtilities::toString( columns[c] ) + ":" + std::string( scored[c].begin(), scored[c].end() );
            auto it = seen.find( key );
            if ( it == seen.end() )
            {
                size_t n = 0;
                for (char x : scored[c]) n += x;
                it = seen.insert( std::make_pair( key, group_states.size() ) ).first;
                group_states.push_back( columns[c] );
                group_scored.push_back( n );
                // a variable character is dropped only with no labels, and with j labels at least j + 1 tips are scored
                group_top.push_back( coding == VARIABLE ? 0 : std::min( columns[c] - 1, n - 1 ) );
                group_offset.push_back( group_total );
                group_total += 3 * ( group_top.back() + 1 );
                first_member.push_back( c );
            }
            group[c] = it->second;
        }
    }

    // at a scored tip its state is the node's: 0 with no labels, or the one label. Elsewhere any state, with no labels
    tip_groups.assign( T, std::vector<double>( group_total, 0.0 ) );
    for (size_t i = 0; i < T; i++)
    {
        for (size_t g = 0; g < group_states.size(); g++)
        {
            double *x = &tip_groups[i][ group_offset[g] ];
            if ( scored[ first_member[g] ][i] )
            {
                x[0] = 1.0;
                if ( group_top[g] > 0 ) x[4] = 1.0;
            }
            else
            {
                x[0] = 1.0;
                x[2] = 1.0;
            }
        }
    }

    binomial.assign( max_states + 1, std::vector<double>( max_states + 1, 0.0 ) );
    for (size_t n = 0; n <= max_states; n++)
    {
        binomial[n][0] = 1.0;
        for (size_t k = 1; k <= n; k++) binomial[n][k] = binomial[n-1][k-1] + ( k < n ? binomial[n-1][k] : 0.0 );
    }

    tips_stale = false;
    epoch++;

    ones.v.assign( total_states, 1.0 );
    ones.ln_scale.assign( N, 0.0 );
    ones.u.assign( group_total, 0.0 );
    for (size_t g = 0; g < group_states.size(); g++)
    {
        ones.u[ group_offset[g] ] = 1.0;
        ones.u[ group_offset[g] + 2 ] = 1.0;
    }
    ones.u_scale.assign( group_states.size(), 0.0 );
    ones.epoch = epoch;
    ones.version = next_version++;
}


/** Apply P = e I + (1 - e) 11'/K to character c. Every Mk transition matrix has this form. */
void PhyloMorphospeciationProcess::mix( std::vector<double> &v, size_t c, double e ) const
{
    double K = double( columns[c] );
    double *x = &v[ offset[c] ];
    double mean = 0.0;
    for (size_t k = 0; k < columns[c]; k++) mean += x[k];
    mean *= (1.0 - e) / K;
    for (size_t k = 0; k < columns[c]; k++) x[k] = e * x[k] + mean;
}


/** mix for group g: every state's partial in a class is the same, so the kernel acts on the classes. */
void PhyloMorphospeciationProcess::mixGroup( std::vector<double> &u, size_t g, double e ) const
{
    double K = double( group_states[g] );
    for (size_t j = 0; j <= group_top[g]; j++)
    {
        double *x = &u[ group_offset[g] + 3*j ];
        double others = K - 1.0 - double( j );
        double mean = ( x[0] + double( j ) * x[1] + others * x[2] ) * (1.0 - e) / K;
        for (int s = 0; s < 3; s++) x[s] = e * x[s] + mean;
        if ( j == 0 ) x[1] = 0.0;
        if ( others == 0.0 ) x[2] = 0.0;
    }
}


double PhyloMorphospeciationProcess::decay( size_t c, double e ) const
{
    double K = double( columns[c] );
    return 1.0 - e * K / (K - 1.0);
}


/**
 * log P(a pattern the coding drops) for group g. By symmetry, with background b and a set of j labels, every
 * choice is as likely as background 0 and labels 1..j: K choices of b and C(K-1, j) of the labels. With j up to
 * n - 2 the background is the one state at two tips or more, but with j = n - 1 every tip differs and any of the
 * n could be the background.
 */
double PhyloMorphospeciationProcess::lnDropped( const Slot &root, size_t g ) const
{
    double K = double( group_states[g] );
    size_t n = group_scored[g];
    const double *x = &root.u[ group_offset[g] ];

    // the root starts from the uniform distribution
    auto U = [&]( size_t j ) { return ( x[3*j] + double( j ) * x[3*j+1] + ( K - 1.0 - double( j ) ) * x[3*j+2] ) / K; };
    double invariant = K * U(0), uninformative = invariant;
    for (size_t j = 1; j <= group_top[g]; j++)
    {
        double ways = K * binomial[ group_states[g] - 1 ][j];
        uninformative += ( j + 1 < n ) ? ways * U(j) : ways * U(j) / double( n );
    }

    double p = invariant;
    if ( coding == INFORMATIVE ) p = ( retention == NULL ) ? uninformative : invariant + (1.0 - retention->getValue()) * ( uninformative - invariant );
    return ( p > 0.0 ) ? std::log( p ) + root.u_scale[g] : RbConstants::Double::neginf;
}


/** log(1 - P(dropped)). */
static double lnKept( double ln_dropped )
{
    return ( ln_dropped < 0.0 ) ? std::log1p( -std::exp( ln_dropped ) ) : RbConstants::Double::neginf;
}


/** Carry each character along a node's branch, or force it to the identity or full randomization. */
void PhyloMorphospeciationProcess::carry( Slot &s, size_t node, BranchForcing branch ) const
{
    if ( branch == IDENTITY ) return;

    // the decay depends on a character only through its number of states
    std::vector<double> by_K( max_states + 1, 1.0 );
    if ( branch == RANDOMIZED )
    {
        std::fill( by_K.begin(), by_K.end(), 0.0 );
    }
    else if ( epsilon == NULL )
    {
        double length = stem[node];
        if ( length <= 0.0 ) return;
        for (size_t K = 2; K <= max_states; K++) by_K[K] = std::exp( -length * double(K) / (double(K) - 1.0) );
    }
    else
    {
        if ( node >= events[BUDDING].size() ) return;
        long nb = events[BUDDING][node], na = events[ANAGENESIS][node];
        if ( nb == 0 && na == 0 ) return;
        for (size_t K = 2; K <= max_states; K++)
        {
            double k = double( K );
            by_K[K] = std::pow( 1.0 - change_prob[BUDDING] * k / (k - 1.0), double( nb ) ) * std::pow( 1.0 - change_prob[ANAGENESIS] * k / (k - 1.0), double( na ) );
        }
    }
    for (size_t c = 0; c < columns.size(); c++) mix( s.v, c, by_K[ columns[c] ] );
    for (size_t g = 0; g < group_states.size(); g++) mixGroup( s.u, g, by_K[ group_states[g] ] );
}


bool PhyloMorphospeciationProcess::changesAt( size_t node, size_t c ) const
{
    return changed[node] == c;
}


void PhyloMorphospeciationProcess::prepareCache( void )
{
    if ( tips_stale ) buildTipPartials();
    size_t nodes = tree->getValue().getNumberOfNodes();
    if ( cache.size() != nodes )
    {
        cache.assign( nodes, std::vector<Slot>( 2 ) );
        current.assign( nodes, 0 );
    }
    markCoded();
}


/** Each node has a coded tip below it, by node index, from which tips are coded. */
std::vector<char> PhyloMorphospeciationProcess::codedNodes( const std::vector<char> &tips ) const
{
    const Tree &t = tree->getValue();
    std::vector<char> out( t.getNumberOfNodes(), 0 );
    std::vector<const TopologyNode*> order;
    order.push_back( &t.getRoot() );
    for (size_t k = 0; k < order.size(); k++)
    {
        for (size_t c = 0; c < order[k]->getNumberOfChildren(); c++) order.push_back( &order[k]->getChild(c) );
    }
    // children come after their parents, so the reverse order fills each node from its children
    for (size_t k = order.size(); k-- > 0; )
    {
        const TopologyNode *n = order[k];
        size_t idx = n->getIndex();
        if ( n->isTip() )
        {
            out[idx] = ( idx < tips.size() ) ? tips[idx] : 0;
            continue;
        }
        for (size_t c = 0; c < n->getNumberOfChildren(); c++) out[idx] |= out[ n->getChild(c).getIndex() ];
    }
    return out;
}


void PhyloMorphospeciationProcess::markCoded( void )
{
    coded = codedNodes( tip_coded );
}


/** Every kernel here has rows summing to one and every change sums to one over its states, so a subtree with
    nothing coded has partial one, whatever its tree, ages or hidden counts, and is never computed. */
/** Whether a node is coded can change with the tree, which would bring back a summed-out count from nothing, so
    only a tip coded for no character, which never becomes coded, has its counts summed out. */
bool PhyloMorphospeciationProcess::summedOut( const TopologyNode &node ) const
{
    return node.isTip() && node.getIndex() < tip_coded.size() && tip_coded[ node.getIndex() ] == 0;
}


const PhyloMorphospeciationProcess::Slot& PhyloMorphospeciationProcess::slotOf( const TopologyNode &node ) const
{
    size_t idx = node.getIndex();
    if ( coded[idx] == 0 ) return ones;
    return cache[idx][ current[idx] ];
}


/** The partial at interior node j from its continuing child a and budding child b, scaled per character. */
void PhyloMorphospeciationProcess::combine( const Slot &a, const Slot &b, long j, Forcing forcing, Slot &out ) const
{
    out.v.resize( total_states );
    out.ln_scale.resize( columns.size() );
    double e = ( epsilon != NULL ) ? change_prob[ anagenetic[j] ? ANAGENESIS : BUDDING ] : 0.0;
    for (size_t c = 0; c < columns.size(); c++)
    {
        // the weight on a change at this character
        double w = ( forcing != NONE ) ? double( forcing == CHANGED ) : ( epsilon != NULL ? e : double( changesAt( j, c ) ) );
        size_t o = offset[c], K = columns[c];
        double total = 0.0;
        if ( w > 0.0 ) for (size_t k = 0; k < K; k++) total += b.v[o+k];

        double largest = 0.0;
        for (size_t k = 0; k < K; k++)
        {
            // a change goes by the Mk kernel, uniform over the other states
            out.v[o+k] = a.v[o+k] * ( (1.0 - w) * b.v[o+k] + ( w > 0.0 ? w * ( total - b.v[o+k] ) / double(K - 1) : 0.0 ) );
            largest = std::max( largest, out.v[o+k] );
        }
        out.ln_scale[c] = a.ln_scale[c] + b.ln_scale[c];
        if ( largest > 0.0 )
        {
            for (size_t k = 0; k < K; k++) out.v[o+k] /= largest;
            out.ln_scale[c] += std::log( largest );
        }
    }

    // each group: the change to b, then its labels split between the children. Where the node's state is a label,
    // that label is at a tip below exactly one child, and to the other it is just another state
    out.u.assign( group_total, 0.0 );
    out.u_scale.resize( group_states.size() );
    std::vector<double> changed_b = b.u;
    for (size_t g = 0; g < group_states.size(); g++)
    {
        double w = ( forcing != NONE ) ? double( forcing == CHANGED ) : e;
        mixGroup( changed_b, g, 1.0 - w * double( group_states[g] ) / ( double( group_states[g] ) - 1.0 ) );
        const double *x = &a.u[ group_offset[g] ], *y = &changed_b[ group_offset[g] ];
        double *z = &out.u[ group_offset[g] ];
        double largest = 0.0;
        for (size_t j = 0; j <= group_top[g]; j++)
        {
            for (size_t i = 0; i <= j; i++)
            {
                const double *p = x + 3*i, *q = y + 3*(j - i);
                double ways = binomial[j][i];
                z[3*j]   += ways * p[0] * q[0];
                z[3*j+2] += ways * p[2] * q[2];
                if ( i > 0 )     z[3*j+1] += binomial[j-1][i-1] * p[1] * q[2];
                if ( i < j )     z[3*j+1] += binomial[j-1][i] * p[2] * q[1];
            }
            for (int s = 0; s < 3; s++) largest = std::max( largest, z[3*j+s] );
        }
        out.u_scale[g] = a.u_scale[g] + b.u_scale[g];
        if ( largest > 0.0 )
        {
            for (size_t k = 0; k < 3 * ( group_top[g] + 1 ); k++) z[k] /= largest;
            out.u_scale[g] += std::log( largest );
        }
    }
}


/** Everything a node's partial depends on apart from its children's partials. */
void PhyloMorphospeciationProcess::signature( const TopologyNode &node, Slot &sig ) const
{
    size_t idx = node.getIndex();
    sig.epoch  = epoch;
    sig.length = stem[idx];
    if ( epsilon != NULL && idx < events[BUDDING].size() )
    {
        for (int m = 0; m < 2; m++)
        {
            sig.events[m] = events[m][idx];
            sig.change_prob[m] = change_prob[m];
        }
    }
    if ( node.isTip() == false )
    {
        sig.kept   = kept_child[idx];
        sig.budded = budded[idx];
        sig.xi = ( epsilon == NULL ) ? changed[idx] : size_t( anagenetic[idx] );
    }
}


/** A node's partial at the top of its branch, recomputed only when something it depends on differs from both cached slots. */
const PhyloMorphospeciationProcess::Slot& PhyloMorphospeciationProcess::evaluate( const TopologyNode &node )
{
    size_t idx = node.getIndex();
    if ( coded[idx] == 0 ) return ones;
    Slot sig;
    signature( node, sig );
    if ( node.isTip() == false )
    {
        for (size_t c = 0; c < 2; c++) sig.child_version[c] = evaluate( node.getChild(c) ).version;
    }

    std::vector<Slot> &slots = cache[idx];
    for (size_t s = 0; s < 2; s++)
    {
        const Slot &x = slots[s];
        if ( x.version != 0 && x.epoch == sig.epoch && x.length == sig.length && x.kept == sig.kept && x.budded == sig.budded &&
             x.xi == sig.xi && x.events[0] == sig.events[0] && x.events[1] == sig.events[1] &&
             x.change_prob[0] == sig.change_prob[0] && x.change_prob[1] == sig.change_prob[1] &&
             x.child_version[0] == sig.child_version[0] && x.child_version[1] == sig.child_version[1] )
        {
            current[idx] = s;
            return x;
        }
    }

    // keep the last state for a rejection to return to
    size_t s = 1 - current[idx];
    Slot &out = slots[s];
    if ( node.isTip() )
    {
        out.v = tip_partials[idx];
        out.ln_scale.assign( columns.size(), 0.0 );
        out.u = tip_groups[idx];
        out.u_scale.assign( group_states.size(), 0.0 );
    }
    else
    {
        const Slot &a = slotOf( node.getChild( sig.kept ) );
        const Slot &b = slotOf( node.getChild( 1 - sig.kept ) );
        combine( a, b, long( idx ), NONE, out );
    }
    carry( out, idx, AS_IS );

    out.epoch = sig.epoch;
    out.length = sig.length;
    out.events[0] = sig.events[0];
    out.events[1] = sig.events[1];
    out.change_prob[0] = sig.change_prob[0];
    out.change_prob[1] = sig.change_prob[1];
    out.kept = sig.kept;
    out.budded = sig.budded;
    out.xi = sig.xi;
    out.child_version[0] = sig.child_version[0];
    out.child_version[1] = sig.child_version[1];
    out.version = next_version++;
    current[idx] = s;
    return out;
}


std::vector<double> PhyloMorphospeciationProcess::rootLnLikelihoods( const Slot &root ) const
{
    // the root starts from the stationary distribution, uniform for Mk, so its stem drops out
    std::vector<double> out( columns.size() );
    for (size_t c = 0; c < columns.size(); c++)
    {
        double total = 0.0;
        for (size_t k = 0; k < columns[c]; k++) total += root.v[ offset[c] + k ];
        out[c] = std::log( total / double( columns[c] ) ) + root.ln_scale[c];
    }
    // then log P(dropped) for each group
    for (size_t g = 0; g < group_states.size(); g++) out.push_back( lnDropped( root, g ) );
    return out;
}


/** log L per used character, then log P(dropped) per group, with one node's speciation or branch forced; the cache must be current. */
std::vector<double> PhyloMorphospeciationProcess::forcedLnLikelihoods( const Forced &f )
{
    const Tree &t = tree->getValue();
    if ( f.node < 0 ) return rootLnLikelihoods( slotOf( t.getRoot() ) );
    const TopologyNode *node = &t.getNode( size_t( f.node ) );

    Slot below, here;
    const Slot *from = NULL;                                // the path's partial at the child just below
    size_t from_idx = 0;
    for ( ; node != NULL; node = node->isRoot() ? NULL : &node->getParent() )
    {
        size_t idx = node->getIndex();
        bool start = ( from == NULL );
        if ( node->isTip() )
        {
            here.v = tip_partials[idx];
            here.ln_scale.assign( columns.size(), 0.0 );
            here.u = tip_groups[idx];
            here.u_scale.assign( group_states.size(), 0.0 );
        }
        else
        {
            size_t k = kept_child[idx];
            const Slot *a = &slotOf( node->getChild(k) );
            const Slot *b = &slotOf( node->getChild(1-k) );
            if ( from != NULL && node->getChild(k).getIndex() == from_idx ) a = from;
            if ( from != NULL && node->getChild(1-k).getIndex() == from_idx ) b = from;
            combine( *a, *b, long( idx ), start ? f.forcing : NONE, here );
        }
        carry( here, idx, start ? f.branch : AS_IS );
        std::swap( below, here );
        from = &below;
        from_idx = idx;
    }
    return rootLnLikelihoods( below );
}


double PhyloMorphospeciationProcess::computeLnProbability( void )
{
    if ( num_states.empty() ) return 0.0;

    // a Gibbs draw just prepared everything, and only its own latent variable has changed since
    if ( fresh == false )
    {
        if ( epsilon != NULL && changeProbabilities( change_prob[BUDDING], change_prob[ANAGENESIS] ) == false ) return RbConstants::Double::neginf;
        if ( prepare() == false ) return RbConstants::Double::neginf;
        if ( epsilon != NULL && latents_ready == false ) initializeLatents();
        prepareCache();
    }
    fresh = false;

    const Tree &t = tree->getValue();
    double ln_retention = 0.0;
    if ( retention != NULL )
    {
        double r = retention->getValue();
        if ( r < 0.0 || r > 1.0 ) return RbConstants::Double::neginf;
        ln_retention = std::log( r );
    }

    // each character's likelihood over its chance to be kept, and with retention times that chance
    double lnProb = 0.0;
    std::vector<double> ln_l = rootLnLikelihoods( evaluate( t.getRoot() ) );
    size_t N = chars.size();
    for (size_t c = 0; c < N; c++)
    {
        lnProb += ln_l[c];
        if ( coding != ALL ) lnProb -= lnKept( ln_l[ N + group[c] ] );
        if ( autapomorphic[c] ) lnProb += ln_retention;
    }

    // every interior node is one speciation
    if ( epsilon == NULL )
    {
        // each speciation's character is drawn uniformly
        size_t m = t.getNumberOfNodes() - t.getNumberOfTips();
        return lnProb - double( m ) * std::log( double( num_states.size() ) );
    }

    // the hidden events on each branch are Poisson, and a branch with no time above tau_1 holds none. Above a tip
    // coded for no character they never touch the likelihood, so they are summed out and their density with them
    for (int m = 0; m < 2; m++)
    {
        for (const TopologyNode *node : t.getNodes())
        {
            if ( node->isRoot() || summedOut( *node ) ) continue;
            size_t k = node->getIndex();
            double mean = mean_events[m][k];
            if ( mean <= 0.0 )
            {
                if ( events[m][k] > 0 ) return RbConstants::Double::neginf;
                continue;
            }
            lnProb += RbStatistics::Poisson::lnPdf( mean, int( events[m][k] ) );
        }
    }

    return lnProb;
}


/**
 * Gibbs draw of the character that changed at one speciation, chosen uniformly.
 *
 * \return The hastings ratio, which cancels the posterior ratio.
 */
double PhyloMorphospeciationProcess::resampleChangedCharacter( void )
{
    if ( prepare() == false ) return RbConstants::Double::neginf;
    prepareCache();
    const Tree &t = tree->getValue();
    evaluate( t.getRoot() );

    // every interior node is one speciation
    std::vector<size_t> eligible;
    for (const TopologyNode *n : t.getNodes()) if ( n->isTip() == false ) eligible.push_back( n->getIndex() );
    if ( eligible.empty() ) return 0.0;
    RandomNumberGenerator *rng = GLOBAL_RNG;
    size_t node = eligible[ std::min( eligible.size() - 1, size_t( rng->uniform01() * eligible.size() ) ) ];

    Forced yes_f, no_f;
    yes_f.node = no_f.node = long( node );
    yes_f.forcing = CHANGED;
    no_f.forcing  = UNCHANGED;
    std::vector<double> yes = forcedLnLikelihoods( yes_f );
    std::vector<double> no  = forcedLnLikelihoods( no_f );

    size_t N = chars.size();
    std::vector<double> ln_probs( N );
    for (size_t c = 0; c < N; c++)
    {
        ln_probs[c] = yes[c] - no[c];
        if ( std::isnan( ln_probs[c] ) ) ln_probs[c] = RbConstants::Double::neginf;
    }
    double top = *std::max_element( ln_probs.begin(), ln_probs.end() );
    if ( RbMath::isFinite( top ) == false ) return RbConstants::Double::neginf;

    std::vector<double> probs( N );
    double total = 0.0;
    for (size_t c = 0; c < N; c++)
    {
        probs[c] = std::exp( ln_probs[c] - top );
        total += probs[c];
    }
    double u = rng->uniform01() * total;
    size_t c_new = 0;
    while ( c_new + 1 < N && u >= probs[c_new] )
    {
        u -= probs[c_new];
        c_new++;
    }

    size_t c_old = changed[node];
    stored_node      = node;
    stored_character = c_old;
    reallocated      = true;
    changed[node]    = c_new;

    return ln_probs[c_old] - ln_probs[c_new];
}


/** With epsilon: rates, tree, hidden counts and cache brought up to date for a Gibbs draw. */
bool PhyloMorphospeciationProcess::ready( void )
{
    if ( changeProbabilities( change_prob[BUDDING], change_prob[ANAGENESIS] ) == false ) return false;
    if ( prepare() == false ) return false;
    if ( latents_ready == false ) initializeLatents();
    prepareCache();
    evaluate( tree->getValue().getRoot() );
    return true;
}


/** Gibbs draw of one kind of hidden event on one branch that can hold events, chosen uniformly. */
double PhyloMorphospeciationProcess::resampleEvents( void )
{
    if ( ready() == false ) return RbConstants::Double::neginf;

    std::vector<size_t> eligible;
    for (const TopologyNode *node : tree->getValue().getNodes())
    {
        size_t i = node->getIndex();
        if ( node->isRoot() == false && summedOut( *node ) == false && ( mean_events[BUDDING][i] > 0.0 || mean_events[ANAGENESIS][i] > 0.0 ) ) eligible.push_back( i );
    }
    if ( eligible.empty() ) return 0.0;

    RandomNumberGenerator *rng = GLOBAL_RNG;
    size_t node = eligible[ std::min( eligible.size() - 1, size_t( rng->uniform01() * eligible.size() ) ) ];
    return drawEvents( node );
}


/**
 * Gibbs draw of the count of one kind of hidden event on one node's branch.
 * L(N) = (1 - x) L0 + x L1 per character, with x the branch's decay.
 *
 * \return The hastings ratio, which cancels the posterior ratio.
 */
double PhyloMorphospeciationProcess::drawEvents( size_t node )
{
    const TopologyNode &nd = tree->getValue().getNode( node );
    if ( nd.isRoot() ) return 0.0;

    std::vector<Mode> kinds;
    if ( mean_events[BUDDING][node] > 0.0 )    kinds.push_back( BUDDING );
    if ( mean_events[ANAGENESIS][node] > 0.0 ) kinds.push_back( ANAGENESIS );
    if ( kinds.empty() ) return 0.0;
    RandomNumberGenerator *rng = GLOBAL_RNG;
    Mode m = kinds[ std::min( kinds.size() - 1, size_t( rng->uniform01() * kinds.size() ) ) ];
    Mode other = ( m == BUDDING ) ? ANAGENESIS : BUDDING;

    double mean = mean_events[m][node];
    long N_old = events[m][node];

    // nothing below is coded, so the likelihood is one whatever the count and the conditional is its prior
    if ( coded[node] == 0 )
    {
        long N_new = RbStatistics::Poisson::rv( mean, *rng );
        undo        = UNDO_EVENTS;
        undo_mode   = m;
        undo_index  = node;
        undo_events = N_old;
        events[m][node] = N_new;
        fresh = true;
        return RbStatistics::Poisson::lnPdf( mean, int( N_old ) ) - RbStatistics::Poisson::lnPdf( mean, int( N_new ) );
    }

    Forced one, zero;
    one.node = zero.node = long( node );
    one.branch  = IDENTITY;
    zero.branch = RANDOMIZED;
    std::vector<double> l1 = forcedLnLikelihoods( one );
    std::vector<double> l0 = forcedLnLikelihoods( zero );

    size_t n = chars.size();
    std::vector<double> base( n ), x( n ), d0( n, 0.0 ), d1( n, 0.0 );
    double ln_bound = 0.0;
    for (size_t c = 0; c < n; c++)
    {
        double K = double( columns[c] );
        base[c] = decay( c, change_prob[m] );
        x[c] = std::pow( decay( c, change_prob[other] ), double( events[other][node] ) );
        // x lies in [-1/(K-1), 1], so L(N) is at most L1 or K/(K-1) L0
        ln_bound += std::max( l1[c], l0[c] + std::log( K / (K - 1.0) ) );

        // and P(dropped), linear in x too, is largest at an end of that interval
        if ( coding != ALL )
        {
            d0[c] = std::exp( l0[ n + group[c] ] );
            d1[c] = std::exp( l1[ n + group[c] ] );
            double worst = std::max( d1[c], ( K * d0[c] - d1[c] ) / (K - 1.0) );
            double ln_kept = lnKept( std::log( std::max( worst, 0.0 ) ) );
            if ( RbMath::isFinite( ln_kept ) ) ln_bound -= ln_kept;
        }
    }

    double ln_mean = std::log( mean );
    std::vector<double> w;
    double total = RbConstants::Double::neginf;
    for (long N = 0; ; N++)
    {
        double lw = double( N ) * ln_mean - mean - std::lgamma( double( N ) + 1.0 );
        for (size_t c = 0; c < n && RbMath::isFinite( lw ); c++)
        {
            double top = std::max( l0[c], l1[c] );
            double v = RbMath::isFinite( top ) ? ( 1.0 - x[c] ) * std::exp( l0[c] - top ) + x[c] * std::exp( l1[c] - top ) : 0.0;
            lw += ( v > 0.0 ) ? top + std::log( v ) : RbConstants::Double::neginf;
            if ( coding != ALL && RbMath::isFinite( lw ) )
            {
                double z = ( 1.0 - x[c] ) * d0[c] + x[c] * d1[c];
                lw -= lnKept( std::log( std::max( z, 0.0 ) ) );
            }
        }
        if ( std::isnan( lw ) ) lw = RbConstants::Double::neginf;
        w.push_back( lw );
        if ( RbMath::isFinite( lw ) ) total = ( RbMath::isFinite( total ) ? std::max( total, lw ) + std::log1p( std::exp( -std::fabs( total - lw ) ) ) : lw );
        for (size_t c = 0; c < n; c++) x[c] *= base[c];

        // P(X >= k) <= P(X = k) (k + 1)/(k + 1 - mean) for k > mean
        double k1 = double( N + 1 );
        if ( N >= N_old && k1 > mean && RbMath::isFinite( total ) )
        {
            double tail = k1 * ln_mean - mean - std::lgamma( k1 + 1.0 ) + std::log( (k1 + 1.0) / (k1 + 1.0 - mean) );
            if ( ln_bound + tail < total - 30.0 ) break;
        }
        if ( N > 10000000 ) break;
    }
    if ( RbMath::isFinite( total ) == false ) return RbConstants::Double::neginf;

    double u = rng->uniform01();
    long N_new = long( w.size() ) - 1;
    for (size_t N = 0; N < w.size(); N++)
    {
        double p = std::exp( w[N] - total );
        if ( u < p )
        {
            N_new = long( N );
            break;
        }
        u -= p;
    }

    undo        = UNDO_EVENTS;
    undo_mode   = m;
    undo_index  = node;
    undo_events = N_old;
    events[m][node] = N_new;
    fresh = true;

    // leaving a state of probability zero is always accepted
    return RbMath::isFinite( w[N_old] ) ? w[N_old] - w[N_new] : 0.0;
}


void PhyloMorphospeciationProcess::keepSpecialization( const DagNode *affecter )
{
    reallocated = false;
    undo = UNDO_NONE;
    fresh = false;
}


/** Anything but this node's own touch, after a Gibbs draw, can change what the draw prepared. */
void PhyloMorphospeciationProcess::touchSpecialization( const DagNode *toucher, bool touchAll )
{
    if ( toucher != this->dag_node ) fresh = false;
}


void PhyloMorphospeciationProcess::restoreSpecialization( const DagNode *restorer )
{
    if ( reallocated ) changed[stored_node] = stored_character;
    reallocated = false;
    fresh = false;

    if ( undo == UNDO_EVENTS ) events[undo_mode][undo_index] = undo_events;
    undo = UNDO_NONE;
}


void PhyloMorphospeciationProcess::setValue( AbstractHomologousDiscreteCharacterData *v, bool force )
{
    TypedDistribution<AbstractHomologousDiscreteCharacterData>::setValue( v, force );
    tips_stale = true;

    // new data may be impossible under the old latent variables
    latents_ready = false;

    if ( states_given )
    {
        if ( this->value->getNumberOfCharacters() != num_states.size() )
        {
            throw RbException() << "dnPhyloMorphospeciation: nStates has " << num_states.size() << " entries for "
                                << this->value->getNumberOfCharacters() << " characters.";
        }
        return;
    }

    // K for each character: the one shared, or its largest observed state and at least two
    size_t N = this->value->getNumberOfCharacters();
    if ( shared_states > 0 )
    {
        num_states.assign( N, shared_states );
        return;
    }
    num_states.assign( N, 2 );
    for (size_t i = 0; i < this->value->getNumberOfTaxa(); i++)
    {
        const AbstractDiscreteTaxonData &td = this->value->getTaxonData( i );
        for (size_t c = 0; c < N; c++)
        {
            const DiscreteCharacterState &x = td.getCharacter( c );
            if ( x.isGapState() || x.isMissingState() ) continue;
            RbBitSet bits = x.getState();
            for (size_t k = 0; k < bits.size(); k++)
            {
                if ( bits.test(k) ) num_states[c] = std::max( num_states[c], k + 1 );
            }
        }
    }
}


void PhyloMorphospeciationProcess::simulateDown( const TopologyNode &node, std::vector<size_t> &state, std::vector< std::vector<size_t> > &tips )
{
    RandomNumberGenerator *rng = GLOBAL_RNG;
    if ( node.isTip() )
    {
        tips[node.getIndex()] = state;
        return;
    }

    // a change draws uniformly from the other states
    auto change = [&]( std::vector<size_t> &s, size_t c ) {
        size_t k = size_t( rng->uniform01() * (num_states[c] - 1) );
        s[c] = ( k >= s[c] ) ? k + 1 : k;
    };

    for (size_t i = 0; i < node.getNumberOfChildren(); i++)
    {
        const TopologyNode &child = node.getChild(i);
        std::vector<size_t> s = state;
        if ( child.continuesParentSpecies() == false )
        {
            size_t j = node.getIndex();
            if ( epsilon == NULL ) change( s, changed[j] );
            else for (size_t c = 0; c < s.size(); c++) if ( rng->uniform01() < change_prob[ anagenetic[j] ? ANAGENESIS : BUDDING ] ) change( s, c );
        }

        if ( epsilon == NULL )
        {
            // along the branch: stay with probability e, or redraw uniformly over all K states
            double length = stem[child.getIndex()];
            if ( length > 0.0 )
            {
                for (size_t c = 0; c < s.size(); c++)
                {
                    double K = double( num_states[c] );
                    if ( rng->uniform01() >= std::exp( -length * K / (K - 1.0) ) ) s[c] = size_t( rng->uniform01() * K );
                }
            }
        }
        else
        {
            // one hidden event at a time, each changing every character with its kind's probability
            for (int m = 0; m < 2; m++)
            {
                for (long e = 0; e < events[m][ child.getIndex() ]; e++)
                {
                    for (size_t c = 0; c < s.size(); c++) if ( rng->uniform01() < change_prob[m] ) change( s, c );
                }
            }
        }
        simulateDown( child, s, tips );
    }
}


void PhyloMorphospeciationProcess::redrawValue( void )
{
    if ( num_states.empty() ) return;
    if ( prepare() == false ) throw RbException("dnPhyloMorphospeciation: the tree has symmetric speciation, which is not part of the model.");

    RandomNumberGenerator *rng = GLOBAL_RNG;
    size_t N = num_states.size();
    const Tree &t = tree->getValue();
    if ( epsilon == NULL )
    {
        changed.assign( t.getNumberOfNodes(), 0 );
        for (size_t &c : changed) c = size_t( rng->uniform01() * N );
    }
    else
    {
        if ( changeProbabilities( change_prob[BUDDING], change_prob[ANAGENESIS] ) == false )
        {
            throw RbException("dnPhyloMorphospeciation: epsilon and epsilon_a must lie in [0, 1].");
        }
        for (int m = 0; m < 2; m++)
        {
            events[m].assign( t.getNumberOfNodes(), 0 );
            for (const TopologyNode *node : t.getNodes())
            {
                if ( node->isRoot() ) continue;
                double mean = mean_events[m][ node->getIndex() ];
                if ( mean > 0.0 ) events[m][ node->getIndex() ] = RbStatistics::Poisson::rv( mean, *rng );
            }
        }
        latents_ready = true;
    }

    std::vector<size_t> root( N );
    for (size_t c = 0; c < N; c++) root[c] = size_t( rng->uniform01() * num_states[c] );
    std::vector< std::vector<size_t> > tips( t.getNumberOfNodes() );
    simulateDown( t.getRoot(), root, tips );

    // a character the coding would drop is drawn again, along the same history, until it is kept. Under
    // retention one varying only by autapomorphies is kept with that probability
    double r = ( retention != NULL ) ? retention->getValue() : 1.0;
    if ( r < 0.0 || r > 1.0 ) throw RbException("dnPhyloMorphospeciation: retention must lie in [0, 1].");
    std::vector<size_t> pending;
    for (size_t c = 0; c < N && coding != ALL; c++) pending.push_back( c );
    for (size_t tries = 0; pending.empty() == false; tries++)
    {
        if ( tries == 100000 ) throw RbException("dnPhyloMorphospeciation: a character was dropped by its coding in 100000 draws, so it cannot be simulated here.");
        std::vector<size_t> again_c;
        for (size_t c : pending)
        {
            std::vector<size_t> counts( num_states[c], 0 );
            for (const TopologyNode *tip : t.getNodes()) if ( tip->isTip() ) counts[ tips[tip->getIndex()][c] ]++;
            size_t present = 0, repeated = 0;
            for (size_t k : counts)
            {
                present += ( k > 0 );
                repeated += ( k > 1 );
            }
            bool kept = present > 1 && ( coding == VARIABLE || repeated > 1 || ( retention != NULL && rng->uniform01() < r ) );
            if ( kept == false ) again_c.push_back( c );
        }
        pending = again_c;
        if ( pending.empty() ) break;

        for (size_t c = 0; c < N; c++) root[c] = size_t( rng->uniform01() * num_states[c] );
        std::vector< std::vector<size_t> > again( t.getNumberOfNodes() );
        simulateDown( t.getRoot(), root, again );
        for (size_t c : pending)
        {
            for (const TopologyNode *tip : t.getNodes()) if ( tip->isTip() ) tips[tip->getIndex()][c] = again[tip->getIndex()][c];
        }
    }

    static const std::string labels = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    HomologousDiscreteCharacterData<StandardState> *data = new HomologousDiscreteCharacterData<StandardState>();
    for (const TopologyNode *tip : t.getNodes())
    {
        if ( tip->isTip() == false ) continue;
        DiscreteTaxonData<StandardState> td( tip->getName() );
        for (size_t c = 0; c < N; c++)
        {
            td.addCharacter( StandardState( labels.substr( tips[tip->getIndex()][c], 1 ), labels.substr( 0, num_states[c] ) ) );
        }
        data->addTaxonData( td );
    }

    delete this->value;
    this->value = data;
    tips_stale = true;
}


void PhyloMorphospeciationProcess::executeMethod( const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv ) const
{
    rv.clear();
    const Tree &t = tree->getValue();
    size_t T = t.getNumberOfTips();

    // the interior node at which each species began, -1 for the species at the origin
    std::vector<long> began( T, -1 );
    for (const TopologyNode *node : t.getNodes())
    {
        if ( node->isTip() ) continue;
        for (size_t c = 0; c < node->getNumberOfChildren(); c++)
        {
            if ( node->getChild(c).continuesParentSpecies() ) continue;
            long j = speciesOf( node->getChild(c) );
            if ( j >= 0 ) began[j] = long( node->getIndex() );
        }
    }

    if ( n == "getChangedCharacters" )
    {
        if ( epsilon != NULL )
        {
            throw RbException("dnPhyloMorphospeciation: with epsilon, any number of characters can change at a speciation, and which ones is summed over.");
        }
        for (size_t i = 0; i < T; i++)
        {
            bool ok = began[i] >= 0 && size_t( began[i] ) < changed.size();
            rv.push_back( ok ? double( changed[ began[i] ] + 1 ) : 0.0 );
        }
    }
    else if ( n == "getSpeciationModes" )
    {
        // 1 for budding, 2 for anagenesis, 0 for the species at the origin
        bool extended = static_cast<const TreeDistributionProperties *>( process )->isExtended();
        std::vector<double> mode( T, 1.0 );
        for (const TopologyNode *node : t.getNodes())
        {
            if ( node->isTip() || node->getNumberOfChildren() != 2 ) continue;
            for (size_t c = 0; c < 2; c++)
            {
                const TopologyNode &kept = node->getChild(c);
                if ( kept.continuesParentSpecies() == false || extended == false || kept.isSampledAncestorTip() == false ) continue;
                long j = speciesOf( node->getChild( 1 - c ) );
                if ( j >= 0 ) mode[j] = 2.0;
            }
        }
        for (size_t i = 0; i < T; i++) rv.push_back( began[i] >= 0 ? mode[i] : 0.0 );
    }
    else if ( n == "getChangeProbabilities" )
    {
        if ( epsilon == NULL ) throw RbException() << "dnPhyloMorphospeciation: " << n << "() needs epsilon.";
        double e_b = 0.0, e_a = 0.0;
        changeProbabilities( e_b, e_a );
        rv.push_back( e_b );
        rv.push_back( e_a );
    }
    else if ( n == "getHiddenEventCounts" || n == "getHiddenEventMeans" )
    {
        if ( epsilon == NULL ) throw RbException() << "dnPhyloMorphospeciation: " << n << "() needs epsilon.";

        // by node index, so entry i is the branch above node i, both kinds of event together
        std::vector<double> b, a;
        std::vector<long> sp( t.getNumberOfNodes() );
        for (size_t i = 0; i < sp.size(); i++) sp[i] = speciesOf( t.getNode( i ) );
        RbVector<double> f;
        process->executeMethod( "getFirstAppearances", std::vector<const DagNode*>(), f );
        std::vector<double> tau( f.size() );
        for (size_t i = 0; i < f.size(); i++) tau[i] = f[i];
        hiddenIntegrals( sp, tau, b, a );

        // where the counts are summed out, draw them from their Poisson marginal, or give its mean
        bool mean = n == "getHiddenEventCounts" && args.size() > 0 && static_cast<const TypedDagNode<Boolean>* >( args[0] )->getValue();
        // a character the coding excludes does not count
        std::vector<size_t> used = chars;
        if ( tips_stale )
        {
            used.clear();
            for (size_t c = 0; c < this->value->getNumberOfCharacters(); c++) used.push_back( c );
        }
        std::vector<char> tip_has( T, 0 );
        for (size_t i = 0; i < T; i++)
        {
            const AbstractDiscreteTaxonData &td = this->value->getTaxonData( t.getTipNode( i ).getName() );
            for (size_t k = 0; k < used.size() && tip_has[ t.getTipNode( i ).getIndex() ] == 0; k++)
            {
                const DiscreteCharacterState &x = td.getCharacter( used[k] );
                if ( x.isGapState() == false && x.isMissingState() == false ) tip_has[ t.getTipNode( i ).getIndex() ] = 1;
            }
        }

        for (size_t i = 0; i < t.getNumberOfNodes(); i++)
        {
            const TopologyNode &node = t.getNode( i );
            if ( node.isRoot() )
            {
                rv.push_back( 0.0 );
            }
            else if ( n == "getHiddenEventMeans" || ( node.isTip() && tip_has[i] == 0 && mean ) )
            {
                rv.push_back( b[i] + a[i] );
            }
            else if ( node.isTip() && tip_has[i] == 0 )
            {
                RandomNumberGenerator *rng = GLOBAL_RNG;
                double draw = 0.0;
                if ( b[i] > 0.0 ) draw += RbStatistics::Poisson::rv( b[i], *rng );
                if ( a[i] > 0.0 ) draw += RbStatistics::Poisson::rv( a[i], *rng );
                rv.push_back( draw );
            }
            else
            {
                bool ok = latents_ready && i < events[BUDDING].size();
                rv.push_back( ok ? double( events[BUDDING][i] + events[ANAGENESIS][i] ) : 0.0 );
            }
        }
    }
    else
    {
        throw RbException() << "dnPhyloMorphospeciation does not have a member method called '" << n << "'.";
    }
}


/** The latent variables by node index, '|'-separated, for clamped characters only. */
std::string PhyloMorphospeciationProcess::getHiddenStateString( void ) const
{
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false ) return "";

    std::stringstream ss;
    if ( epsilon == NULL )
    {
        ss << "single";
        for (size_t i = 0; i < changed.size(); i++) ss << ( i == 0 ? "|" : "," ) << changed[i];
        return ss.str();
    }

    if ( latents_ready == false ) return "";
    ss << "multi";
    for (int m = 0; m < 2; m++)
    {
        ss << "|";
        for (size_t k = 0; k < events[m].size(); k++) ss << ( k == 0 ? "" : "," ) << events[m][k];
    }
    return ss.str();
}


void PhyloMorphospeciationProcess::setHiddenStateFromString( const std::string &s )
{
    std::vector<std::string> fields;
    StringUtilities::stringSplit( s, "|", fields );
    size_t nodes = tree->getValue().getNumberOfNodes();

    auto split = [&]( const std::string &x ) {
        std::vector<std::string> out;
        StringUtilities::stringSplit( x, ",", out );
        return out;
    };

    if ( epsilon == NULL )
    {
        std::vector<std::string> xi = ( fields.size() == 2 && fields[0] == "single" ) ? split( fields[1] ) : std::vector<std::string>();
        if ( xi.size() != nodes ) throw RbException("dnPhyloMorphospeciation: the checkpoint does not hold one changed character per node.");
        changed.assign( nodes, 0 );
        for (size_t i = 0; i < nodes; i++) changed[i] = std::stoul( xi[i] );
        return;
    }

    if ( fields.size() != 3 || fields[0] != "multi" ) throw RbException("dnPhyloMorphospeciation: the checkpoint does not hold hidden counts.");
    for (int m = 0; m < 2; m++)
    {
        std::vector<std::string> counts = split( fields[1 + m] );
        if ( counts.size() != nodes ) throw RbException("dnPhyloMorphospeciation: the checkpoint does not hold hidden counts for every branch.");
        events[m].assign( nodes, 0 );
        for (size_t k = 0; k < nodes; k++) events[m][k] = std::stol( counts[k] );
    }
    latents_ready = true;
    fresh = false;
}


void PhyloMorphospeciationProcess::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    if ( oldP == tree )
    {
        tree = static_cast<const TypedDagNode<Tree>* >( newP );
        resolveProcess();
    }
    if ( oldP == epsilon )
    {
        epsilon = static_cast<const TypedDagNode<double>* >( newP );
    }
    if ( oldP == epsilon_a )
    {
        epsilon_a = static_cast<const TypedDagNode<double>* >( newP );
    }
    if ( oldP == retention )
    {
        retention = static_cast<const TypedDagNode<double>* >( newP );
    }
}
