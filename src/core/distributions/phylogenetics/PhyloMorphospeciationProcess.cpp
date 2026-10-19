#include "PhyloMorphospeciationProcess.h"

#include <algorithm>
#include <cmath>

#include "DiscreteTaxonData.h"
#include "FossilizedBirthDeathSpeciationProcess.h"
#include "HomologousDiscreteCharacterData.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "StandardState.h"
#include "StochasticNode.h"
#include "StringUtilities.h"
#include "TopologyNode.h"

using namespace RevBayesCore;


PhyloMorphospeciationProcess::PhyloMorphospeciationProcess(const TypedDagNode<Tree> *t, const std::vector<size_t> &k) :
    TypedDistribution<AbstractHomologousDiscreteCharacterData>( new HomologousDiscreteCharacterData<StandardState>() ),
    tree( t ),
    process( NULL ),
    num_states( k ),
    states_given( k.empty() == false )
{
    resolveProcess();

    addParameter( tree );
    for (const DagNode *n : process->getRateParameters()) addParameter( n );

    for (size_t K : num_states)
    {
        if ( K < 2 ) throw RbException("dnPhyloMorphospeciation: every character needs at least two states.");
    }

    changed = std::vector<size_t>( tree->getValue().getNumberOfTips(), 0 );

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


bool PhyloMorphospeciationProcess::hasSpeciationEvent( size_t taxon ) const
{
    return long( taxon ) != speciesOf( tree->getValue().getRoot() );
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


/** tau_1 of every taxon, the species of every node, and the Markov chain length of every branch. */
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

    double n = double( num_states.size() );
    stem.assign( nodes, 0.0 );
    std::vector<size_t> at;
    std::vector<std::pair<double,double> > spans;
    for (const TopologyNode *node : t.getNodes())
    {
        if ( node->isRoot() ) continue;

        // the species' sampled range carries no change, so only the part above tau_1 counts
        double lo = std::max( node->getAge(), first[ species[ node->getIndex() ] ] );
        double hi = node->getParent().getAge();
        if ( hi > lo && n > 0 )
        {
            at.push_back( node->getIndex() );
            spans.push_back( std::make_pair( lo, hi ) );
        }
    }
    std::vector<double> b = process->hiddenSpeciationIntegrals( spans );
    for (size_t k = 0; k < at.size(); k++) stem[ at[k] ] = b[k] / n;
    return true;
}


/** The observed states of every tip, laid out flat, rebuilt only when the data change. */
void PhyloMorphospeciationProcess::buildTipPartials( void )
{
    size_t N = num_states.size();
    offset.assign( N, 0 );
    total_states = 0;
    for (size_t c = 0; c < N; c++)
    {
        offset[c] = total_states;
        total_states += num_states[c];
    }

    const Tree &t = tree->getValue();
    tip_partials.assign( t.getNumberOfTips(), std::vector<double>( total_states, 1.0 ) );
    for (size_t i = 0; i < t.getNumberOfTips(); i++)
    {
        const TopologyNode &tip = t.getTipNode( i );
        const AbstractDiscreteTaxonData &td = this->value->getTaxonData( tip.getName() );
        std::vector<double> &v = tip_partials[ tip.getIndex() ];
        for (size_t c = 0; c < N; c++)
        {
            const DiscreteCharacterState &x = td.getCharacter( c );
            if ( x.isGapState() || x.isMissingState() ) continue;
            RbBitSet bits = x.getState();
            for (size_t k = 0; k < num_states[c]; k++) v[ offset[c] + k ] = 0.0;
            for (size_t k = 0; k < bits.size(); k++)
            {
                if ( bits.test(k) == false ) continue;
                if ( k >= num_states[c] ) throw RbException() << "dnPhyloMorphospeciation: " << tip.getName() << " has state " << k
                                                              << " at character " << c + 1 << ", which has " << num_states[c] << " states.";
                v[ offset[c] + k ] = 1.0;
            }
        }
    }
    tips_stale = false;
    epoch++;
}


/** Carry each character along a Mk chain of the given length: P = e I + (1 - e) 11'/K. */
void PhyloMorphospeciationProcess::carry( std::vector<double> &v, double length ) const
{
    if ( length <= 0.0 ) return;
    for (size_t c = 0; c < num_states.size(); c++)
    {
        double K = double( num_states[c] );
        double e = std::exp( -length * K / (K - 1.0) );
        double *x = &v[ offset[c] ];
        double mean = 0.0;
        for (size_t k = 0; k < num_states[c]; k++) mean += x[k];
        mean *= (1.0 - e) / K;
        for (size_t k = 0; k < num_states[c]; k++) x[k] = e * x[k] + mean;
    }
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
}


/** The partial at node j's parent from its continuing child a and budding child b, scaled per character. */
void PhyloMorphospeciationProcess::combine( const Slot &a, const Slot &b, long j, Forcing forcing, Slot &out ) const
{
    out.v.resize( total_states );
    out.ln_scale.resize( num_states.size() );
    for (size_t c = 0; c < num_states.size(); c++)
    {
        bool here = ( forcing == NONE ) ? changed[j] == c : forcing == CHANGED;
        size_t o = offset[c], K = num_states[c];
        double total = 0.0;
        if ( here ) for (size_t k = 0; k < K; k++) total += b.v[o+k];

        double largest = 0.0;
        for (size_t k = 0; k < K; k++)
        {
            // the new species differs at this character: the Mk kernel, uniform over the other states
            out.v[o+k] = a.v[o+k] * ( here ? ( total - b.v[o+k] ) / double(K - 1) : b.v[o+k] );
            largest = std::max( largest, out.v[o+k] );
        }
        out.ln_scale[c] = a.ln_scale[c] + b.ln_scale[c];
        if ( largest > 0.0 )
        {
            for (size_t k = 0; k < K; k++) out.v[o+k] /= largest;
            out.ln_scale[c] += std::log( largest );
        }
    }
}


/** A node's partial at the top of its branch, recomputed only when something it depends on differs from both cached slots. */
const PhyloMorphospeciationProcess::Slot& PhyloMorphospeciationProcess::evaluate( const TopologyNode &node )
{
    size_t idx = node.getIndex();
    Slot sig;
    sig.epoch  = epoch;
    sig.length = stem[idx];
    if ( node.isTip() == false )
    {
        sig.kept   = kept_child[idx];
        sig.budded = budded[idx];
        sig.xi     = changed[ budded[idx] ];
        for (size_t c = 0; c < 2; c++) sig.child_version[c] = evaluate( node.getChild(c) ).version;
    }

    std::vector<Slot> &slots = cache[idx];
    for (size_t s = 0; s < 2; s++)
    {
        const Slot &x = slots[s];
        if ( x.version != 0 && x.epoch == sig.epoch && x.length == sig.length && x.kept == sig.kept && x.budded == sig.budded &&
             x.xi == sig.xi && x.child_version[0] == sig.child_version[0] && x.child_version[1] == sig.child_version[1] )
        {
            current[idx] = s;
            return x;
        }
    }

    // overwrite the slot not used last, which keeps the last state for a rejection to return to
    size_t s = 1 - current[idx];
    Slot &out = slots[s];
    if ( node.isTip() )
    {
        out.v = tip_partials[idx];
        out.ln_scale.assign( num_states.size(), 0.0 );
    }
    else
    {
        const Slot &a = cache[ node.getChild( sig.kept ).getIndex() ][ current[ node.getChild( sig.kept ).getIndex() ] ];
        const Slot &b = cache[ node.getChild( 1 - sig.kept ).getIndex() ][ current[ node.getChild( 1 - sig.kept ).getIndex() ] ];
        combine( a, b, sig.budded, NONE, out );
    }
    carry( out.v, sig.length );

    out.epoch = sig.epoch;
    out.length = sig.length;
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
    std::vector<double> out( num_states.size() );
    for (size_t c = 0; c < num_states.size(); c++)
    {
        double total = 0.0;
        for (size_t k = 0; k < num_states[c]; k++) total += root.v[ offset[c] + k ];
        out[c] = std::log( total / double( num_states[c] ) ) + root.ln_scale[c];
    }
    return out;
}


/** log L per character with one species' change forced. Only the path from its speciation to the root differs from the cache. */
std::vector<double> PhyloMorphospeciationProcess::forcedLnLikelihoods( size_t taxon, Forcing forcing )
{
    const Tree &t = tree->getValue();
    const TopologyNode *node = NULL;
    for (const TopologyNode *n : t.getNodes())
    {
        if ( n->isTip() == false && budded[ n->getIndex() ] == long( taxon ) ) node = n;
    }
    if ( node == NULL ) return rootLnLikelihoods( cache[ t.getRoot().getIndex() ][ current[ t.getRoot().getIndex() ] ] );

    Slot below, here;
    const Slot *from = NULL;                                // the path's partial at the child just below
    size_t from_idx = 0;
    for ( ; node != NULL; node = node->isRoot() ? NULL : &node->getParent() )
    {
        size_t idx = node->getIndex();
        size_t k = kept_child[idx];
        const Slot *a = &cache[ node->getChild(k).getIndex() ][ current[ node->getChild(k).getIndex() ] ];
        const Slot *b = &cache[ node->getChild(1-k).getIndex() ][ current[ node->getChild(1-k).getIndex() ] ];
        if ( from != NULL && node->getChild(k).getIndex() == from_idx ) a = from;
        if ( from != NULL && node->getChild(1-k).getIndex() == from_idx ) b = from;

        combine( *a, *b, budded[idx], from == NULL ? forcing : NONE, here );
        carry( here.v, stem[idx] );
        std::swap( below, here );
        from = &below;
        from_idx = idx;
    }
    return rootLnLikelihoods( below );
}


double PhyloMorphospeciationProcess::computeLnProbability( void )
{
    if ( num_states.empty() ) return 0.0;
    if ( prepare() == false ) return RbConstants::Double::neginf;
    prepareCache();

    const Tree &t = tree->getValue();
    double lnProb = 0.0;
    for (double x : rootLnLikelihoods( evaluate( t.getRoot() ) )) lnProb += x;

    // each speciation's character is drawn uniformly, and every interior node is one speciation
    size_t m = t.getNumberOfNodes() - t.getNumberOfTips();
    return lnProb - double( m ) * std::log( double( num_states.size() ) );
}


/** The full conditional of one species' changed character, from one pass with it changed everywhere and one with it nowhere. */
std::vector<double> PhyloMorphospeciationProcess::changedCharacterLnProbs( size_t taxon )
{
    if ( prepare() == false ) return std::vector<double>( num_states.size(), RbConstants::Double::neginf );
    prepareCache();
    evaluate( tree->getValue().getRoot() );

    std::vector<double> yes = forcedLnLikelihoods( taxon, CHANGED );
    std::vector<double> no  = forcedLnLikelihoods( taxon, UNCHANGED );

    double none = 0.0;
    for (double x : no) none += x;

    std::vector<double> out( yes.size() );
    for (size_t c = 0; c < yes.size(); c++) out[c] = none + yes[c] - no[c];
    return out;
}


void PhyloMorphospeciationProcess::setChangedCharacter( size_t taxon, size_t c )
{
    stored_taxon     = taxon;
    stored_character = changed[taxon];
    reallocated      = true;
    changed[taxon]   = c;
}


void PhyloMorphospeciationProcess::keepSpecialization( const DagNode *affecter )
{
    reallocated = false;
}


void PhyloMorphospeciationProcess::restoreSpecialization( const DagNode *restorer )
{
    if ( reallocated ) changed[stored_taxon] = stored_character;
    reallocated = false;
}


void PhyloMorphospeciationProcess::setValue( AbstractHomologousDiscreteCharacterData *v, bool force )
{
    TypedDistribution<AbstractHomologousDiscreteCharacterData>::setValue( v, force );
    tips_stale = true;

    if ( states_given )
    {
        if ( this->value->getNumberOfCharacters() != num_states.size() )
        {
            throw RbException() << "dnPhyloMorphospeciation: nStates has " << num_states.size() << " entries for "
                                << this->value->getNumberOfCharacters() << " characters.";
        }
        return;
    }

    // K for each character: its largest observed state, and at least two
    size_t N = this->value->getNumberOfCharacters();
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
    for (size_t i = 0; i < node.getNumberOfChildren(); i++)
    {
        const TopologyNode &child = node.getChild(i);
        std::vector<size_t> s = state;
        if ( child.continuesParentSpecies() == false )
        {
            size_t c = changed[ speciesOf( child ) ];
            size_t K = num_states[c];
            size_t k = size_t( rng->uniform01() * (K - 1) );
            s[c] = ( k >= s[c] ) ? k + 1 : k;
        }
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
        simulateDown( child, s, tips );
    }
}


void PhyloMorphospeciationProcess::redrawValue( void )
{
    if ( num_states.empty() ) return;
    if ( prepare() == false ) throw RbException("dnPhyloMorphospeciation: the tree has symmetric speciation, which is not part of the model.");

    RandomNumberGenerator *rng = GLOBAL_RNG;
    size_t N = num_states.size();
    for (size_t &c : changed) c = size_t( rng->uniform01() * N );

    const Tree &t = tree->getValue();
    std::vector<size_t> root( N );
    for (size_t c = 0; c < N; c++) root[c] = size_t( rng->uniform01() * num_states[c] );
    std::vector< std::vector<size_t> > tips( t.getNumberOfNodes() );
    simulateDown( t.getRoot(), root, tips );

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
    if ( n != "getChangedCharacters" )
    {
        throw RbException() << "dnPhyloMorphospeciation does not have a member method called '" << n << "'.";
    }
    rv.clear();
    for (size_t i = 0; i < changed.size(); i++)
    {
        rv.push_back( hasSpeciationEvent( i ) ? double( changed[i] + 1 ) : 0.0 );
    }
}


/**
 * The changed character of each species, which a checkpoint cannot recover from the data. Only for
 * clamped characters, whose value comes back with the data.
 */
std::string PhyloMorphospeciationProcess::getHiddenStateString( void ) const
{
    if ( this->dag_node == NULL || this->dag_node->isClamped() == false ) return "";

    std::stringstream ss;
    ss << "single";
    for (size_t i = 0; i < changed.size(); i++) ss << ( i == 0 ? "|" : "," ) << changed[i];
    return ss.str();
}


void PhyloMorphospeciationProcess::setHiddenStateFromString( const std::string &s )
{
    std::vector<std::string> fields, xi;
    StringUtilities::stringSplit( s, "|", fields );
    if ( fields.size() == 2 && fields[0] == "single" ) StringUtilities::stringSplit( fields[1], ",", xi );

    size_t T = tree->getValue().getNumberOfTips();
    if ( xi.size() != T ) throw RbException("dnPhyloMorphospeciation: the checkpoint does not hold one changed character per species.");
    for (size_t i = 0; i < T; i++) changed[i] = std::stoul( xi[i] );
}


void PhyloMorphospeciationProcess::swapParameterInternal( const DagNode *oldP, const DagNode *newP )
{
    if ( oldP == tree )
    {
        tree = static_cast<const TypedDagNode<Tree>* >( newP );
        resolveProcess();
    }
}
