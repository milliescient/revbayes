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


/** tau_1 of every taxon, and the Markov chain length of every branch. False unless every node keeps
    the parent species on exactly one child, since anything else has probability zero under the model. */
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
    double n = double( num_states.size() );
    stem = std::vector<double>( t.getNumberOfNodes(), 0.0 );
    std::vector<size_t> at;
    std::vector<std::pair<double,double> > spans;
    kept_child.assign( t.getNumberOfNodes(), 0 );
    budded.assign( t.getNumberOfNodes(), -1 );
    for (const TopologyNode *node : t.getNodes())
    {
        if ( node->isTip() == false )
        {
            if ( node->getNumberOfChildren() != 2 ) return false;
            size_t continuing = 0;
            for (size_t c = 0; c < 2; c++)
            {
                if ( node->getChild(c).continuesParentSpecies() )
                {
                    continuing++;
                    kept_child[ node->getIndex() ] = c;
                }
            }
            if ( continuing != 1 ) return false;
            budded[ node->getIndex() ] = speciesOf( node->getChild( 1 - kept_child[ node->getIndex() ] ) );
        }
        if ( node->isRoot() ) continue;
        long s = speciesOf( *node );
        if ( s < 0 ) return false;

        // the species' sampled range carries no change, so only the part above tau_1 counts
        double lo = std::max( node->getAge(), first[s] );
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


void PhyloMorphospeciationProcess::partial( const TopologyNode &node, std::vector<double> &ln_scale, long forced_taxon, Forcing forcing )
{
    size_t idx = node.getIndex();
    std::vector<double> &v = buffer[idx];
    if ( node.isTip() )
    {
        v = tip_partials[idx];
        return;
    }

    const TopologyNode &kc = node.getChild( kept_child[idx] );
    const TopologyNode &bc = node.getChild( 1 - kept_child[idx] );
    partial( kc, ln_scale, forced_taxon, forcing );
    partial( bc, ln_scale, forced_taxon, forcing );
    carry( buffer[ kc.getIndex() ], stem[ kc.getIndex() ] );
    carry( buffer[ bc.getIndex() ], stem[ bc.getIndex() ] );
    const std::vector<double> &a = buffer[ kc.getIndex() ];
    const std::vector<double> &b = buffer[ bc.getIndex() ];

    long j = budded[idx];
    for (size_t c = 0; c < num_states.size(); c++)
    {
        bool here = ( j == forced_taxon ) ? forcing == CHANGED : changed[j] == c;
        size_t o = offset[c], K = num_states[c];
        double total = 0.0;
        if ( here ) for (size_t k = 0; k < K; k++) total += b[o+k];

        double largest = 0.0;
        for (size_t k = 0; k < K; k++)
        {
            // the new species differs at this character: the Mk kernel, uniform over the other states
            v[o+k] = a[o+k] * ( here ? ( total - b[o+k] ) / double(K - 1) : b[o+k] );
            largest = std::max( largest, v[o+k] );
        }
        if ( largest > 0.0 )
        {
            for (size_t k = 0; k < K; k++) v[o+k] /= largest;
            ln_scale[c] += std::log( largest );
        }
    }
}


std::vector<double> PhyloMorphospeciationProcess::lnLikelihoodByCharacter( long forced_taxon, Forcing forcing )
{
    size_t N = num_states.size();
    const Tree &t = tree->getValue();
    if ( tips_stale ) buildTipPartials();
    if ( buffer.size() != t.getNumberOfNodes() || buffer[0].size() != total_states )
    {
        buffer.assign( t.getNumberOfNodes(), std::vector<double>( total_states, 0.0 ) );
    }

    std::vector<double> ln_scale( N, 0.0 );
    const TopologyNode &root = t.getRoot();
    partial( root, ln_scale, forced_taxon, forcing );
    const std::vector<double> &v = buffer[ root.getIndex() ];

    // the root starts from the stationary distribution, uniform for Mk, so its stem drops out
    std::vector<double> out( N );
    for (size_t c = 0; c < N; c++)
    {
        double total = 0.0;
        for (size_t k = 0; k < num_states[c]; k++) total += v[ offset[c] + k ];
        out[c] = std::log( total / double( num_states[c] ) ) + ln_scale[c];
    }
    return out;
}


double PhyloMorphospeciationProcess::computeLnProbability( void )
{
    if ( num_states.empty() ) return 0.0;
    if ( prepare() == false ) return RbConstants::Double::neginf;

    double lnProb = 0.0;
    for (double x : lnLikelihoodByCharacter( -1, NONE )) lnProb += x;

    // each speciation's character is drawn uniformly
    size_t m = 0;
    for (size_t i = 0; i < changed.size(); i++) m += hasSpeciationEvent( i );
    return lnProb - double( m ) * std::log( double( num_states.size() ) );
}


/** The full conditional of one species' changed character, from one pass with it changed everywhere and one with it nowhere. */
std::vector<double> PhyloMorphospeciationProcess::changedCharacterLnProbs( size_t taxon )
{
    if ( prepare() == false ) return std::vector<double>( num_states.size(), RbConstants::Double::neginf );
    std::vector<double> yes = lnLikelihoodByCharacter( long( taxon ), CHANGED );
    std::vector<double> no  = lnLikelihoodByCharacter( long( taxon ), UNCHANGED );

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
