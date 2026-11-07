#include <optional>
#include <array>
#include <functional>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <ostream>
#include <string>
#include <type_traits>
#include <vector>
#include <map>
#include <limits>

#include "AbstractBirthDeathProcess.h"
#include "AbstractFossilizedBirthDeathRangeProcess.h"
#include "DistributionExponential.h"
#include "FossilizedBirthDeathSpeciationProcess.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbMathCombinatorialFunctions.h"
#include "RbException.h"
#include "RbMathFunctions.h"
#include "RbMathLogic.h"
#include "RlUserInterface.h"
#include "StringUtilities.h"
#include "StochasticNode.h"
#include "Taxon.h"
#include "TimeInterval.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TypedDagNode.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { template <class valueType> class RbVector; }

using namespace RevBayesCore;

/**
 * Constructor. 
 * We delegate most parameters to the base class and initialize the members.
 *
 * \param[in]    s              Speciation rates.
 * \param[in]    e              Extinction rates.
 * \param[in]    p              Fossil sampling rates.
 * \param[in]    c              Fossil observation counts.
 * \param[in]    r              Instantaneous sampling probabilities.
 * \param[in]    a              Anagenetic speciation rates.
 * \param[in]    b              Symmetric speciation probability.
 * \param[in]    t              Rate change times.
 * \param[in]    cdt            Condition of the process (time/sampling/survival).
 * \param[in]    tn             Taxa.
 * \param[in]    c              Complete sampling?
 */
FossilizedBirthDeathSpeciationProcess::FossilizedBirthDeathSpeciationProcess(const TypedDagNode<double> *ra,
                                                           const DagNode *inspeciation,
                                                           const DagNode *inextinction,
                                                           const DagNode *inpsi,
                                                           const TypedDagNode<double> *inrho,
                                                           const DagNode *inlambda_a,
                                                           const DagNode *inbeta,
                                                           const TypedDagNode< RbVector<double> > *intimes,
                                                           const std::string &incondition,
                                                           const std::vector<Taxon> &intaxa,
                                                           bool comp,
                                                           double inpresent,
                                                           const TypedDagNode<RbVector<double> >* inmixture_rates,
                                                           const TypedDagNode<Simplex>* inmixture_weights) :
    AbstractBirthDeathProcess(ra, incondition, intaxa, true, NULL),
    AbstractFossilizedBirthDeathRangeProcess(inspeciation, inextinction, inpsi, inrho, intimes, incondition, intaxa, comp, NULL, inpresent, inmixture_rates, inmixture_weights)
{
    for(std::vector<const DagNode*>::iterator it = range_parameters.begin(); it != range_parameters.end(); it++)
    {
        addParameter(*it);
    }

    homogeneous_lambda_a             = NULL;
    homogeneous_beta                 = NULL;
    heterogeneous_lambda_a           = NULL;
    heterogeneous_beta               = NULL;

    heterogeneous_lambda_a = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inlambda_a);
    homogeneous_lambda_a   = dynamic_cast<const TypedDagNode<double >*>(inlambda_a);
    heterogeneous_beta     = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inbeta);
    homogeneous_beta       = dynamic_cast<const TypedDagNode<double >*>(inbeta);

    addParameter( homogeneous_lambda_a );
    addParameter( heterogeneous_lambda_a );
    addParameter( homogeneous_beta );
    addParameter( heterogeneous_beta );

    RbException no_timeline_err = RbException("No time intervals provided for heterogeneous fossilized birth death process");

    RbException inconsistent_rates = RbException("Inconsistent number of rates in fossilized birth death process.");

    if ( heterogeneous_lambda_a != NULL )
    {
        if ( timeline == NULL ) throw no_timeline_err;

        if ( heterogeneous_lambda_a->getValue().size() != num_intervals ) throw inconsistent_rates;
    }
    if ( heterogeneous_beta != NULL )
    {
        if ( timeline == NULL ) throw no_timeline_err;

        if ( heterogeneous_beta->getValue().size() != num_intervals ) throw inconsistent_rates;
    }

    species        = std::vector<SpeciesEntry>(taxa.size());
    next_species   = std::vector<SpeciesEntry>(taxa.size());
    stored_species = std::vector<SpeciesEntry>(taxa.size());

    anagenetic    = std::vector<double>(num_intervals, 0.0);
    symmetric     = std::vector<double>(num_intervals, 0.0);

    
    redrawValue();

    updateRanges();
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself 
 */
FossilizedBirthDeathSpeciationProcess* FossilizedBirthDeathSpeciationProcess::clone( void ) const
{
    return new FossilizedBirthDeathSpeciationProcess( *this );
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double FossilizedBirthDeathSpeciationProcess::computeLnProbabilityDivergenceTimes( void ) const
{
    double lnProb = const_cast<FossilizedBirthDeathSpeciationProcess*>(this)->computeLnProbabilityRanges();

    // a tree this process cannot produce
    if ( invalid_continuation == true || hasTiedEvents() == true ) return RbConstants::Double::neginf;

    lnProb += computeLnProbabilityTimes();

    lnProb += speciationModeLnProb();

    return lnProb;
}


/**
 * Two events at one age have probability zero, except a symmetric speciation at its sampled ancestor's age, which the
 * range carries into, and the tips at the present.
 */
bool FossilizedBirthDeathSpeciationProcess::hasTiedEvents( void ) const
{
    double present = times.front();
    std::vector<double> ages;
    for (const TopologyNode *n : value->getNodes())
    {
        if ( n->isSampledAncestorTip() == true || ( n->isTip() == true && n->getAge() == present ) ) continue;
        if ( n->isRoot() == false && n->isTip() == false && n->getAge() == n->getParent().getAge()
             && n->getParent().isSampledAncestorParent() == true && n->continuesParentSpecies() == true ) continue;
        ages.push_back( n->getAge() );
    }
    std::sort( ages.begin(), ages.end() );

    return std::adjacent_find( ages.begin(), ages.end() ) != ages.end();
}


/**
 * A budding pays 1-beta; a symmetric speciation pays 2 beta lambda for both daughters, who each already carry
 * a lambda at their birth, since its daughters have no order while a budding's two orientations are separate
 * states. A range no taxon samples pays its birth and q down to the symmetric node that ends it.
 */
double FossilizedBirthDeathSpeciationProcess::speciationModeLnProb( void ) const
{
    double lnProb = 0.0;

    for (double a : budding_ages)
    {
        lnProb += log( 1.0 - symmetric[ findIndex(a) ] );
    }

    for (double a : symmetric_ages)
    {
        size_t i = findIndex(a);
        lnProb += log( 2.0 * symmetric[i] ) - log( birth[i] );
    }

    for (const UnsampledRange &r : unsampled)
    {
        size_t ti = findIndex(r.top);
        size_t ei = findIndex(r.end);

        if ( num_components > 1 )
        {
            lnProb += unsampledMixtureLnProb( r, mixture_weights );
        }
        else
        {
            lnProb += lnQFromPresent( ti, r.top ) - lnQFromPresent( ei, r.end );
        }

        if ( r.start == RangeStart::SPECIATION ) lnProb += log( birth[ti] );

        // the ancestor's end was charged mu, as computeLnProbabilityTimes cancels for a taxon
        if ( r.start == RangeStart::SAMPLED_ANCESTOR )
        {
            double r_ana, r_sym;
            successorRates( r.top, r_ana, r_sym );
            lnProb += log( r.anagenetic ? r_ana : r_sym ) - log( death[ti] );
        }
    }

    return lnProb;
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double FossilizedBirthDeathSpeciationProcess::computeLnProbabilityTimes( void ) const
{
    double lnProb = 0.0;

    // a successor's range began where its sampled ancestor's ended, by anagenesis or by a symmetric speciation
    // whose other daughter left no sample, as its flag says; this cancels the ancestor's mu and the successor's lambda
    for (size_t i = 0; i < taxa.size(); i++)
    {
        if ( species[i].parent_is_sa == false ) continue;

        double y = ranges[i].divergence;
        size_t yi = findIndex( y );
        double r_ana, r_sym;
        successorRates( y, r_ana, r_sym );
        lnProb += log( species[i].sa_anagenetic ? r_ana : r_sym ) - log( birth[yi] ) - log( death[yi] );
    }

    return lnProb;
}


/**
 * Close range i at d. The symmetric node pays for its own; a successor's term cancels the mu at a sampled ancestor;
 * otherwise extinction, or an anagenetic or symmetric speciation that leaves no sampled descendant, each successor
 * drawing its own psi component.
 */
double FossilizedBirthDeathSpeciationProcess::rangeEndTerm( size_t i, size_t di, double d, double p_d, double pbar_d ) const
{
    if ( species[i].ends_symmetric == true ) return 0.0;

    if ( species[i].ends_at_sa == true ) return log( death[di] );

    double pd = ( p_d >= 0.0 ? p_d : p( di, d ) );
    double pb = ( pbar_d >= 0.0 ? pbar_d : pd );
    return log( death[di] + anagenetic[di] * pb + birth[di] * symmetric[di] * pb * pb );
}


/** The range terms also read the anagenetic rate and the chance of symmetric speciation. */
void FossilizedBirthDeathSpeciationProcess::appendIntervalRates( size_t j, std::vector<double> &v ) const
{
    AbstractFossilizedBirthDeathRangeProcess::appendIntervalRates( j, v );
    v.push_back( anagenetic[j] );
    v.push_back( symmetric[j] );
}


/** At a sampled ancestor of age y, the successor's start: anagenesis, or a symmetric speciation whose other daughter left no sample. */
void FossilizedBirthDeathSpeciationProcess::successorRates( double y, double &r_ana, double &r_sym ) const
{
    size_t yi = findIndex( y );
    double pbar = 0.0;
    if ( num_components > 1 )
    {
        std::vector<double> ph;
        double lp;
        mixtureAt( y, ph, lp );
        for (size_t h = 0; h < num_components; h++) pbar += mixture_weights[h] * ph[h];
    }
    else
    {
        pbar = p( yi, y );
    }
    r_ana = anagenetic[yi];
    r_sym = 2.0 * birth[yi] * symmetric[yi] * pbar;
}


double FossilizedBirthDeathSpeciationProcess::getMaxTaxonAge( const TopologyNode& node ) const
{
    if( node.isTip() )
    {
        return ranges[node.getIndex()].first;
    }
    else
    {
        double max = 0;

        for( size_t i = 0; i < node.getNumberOfChildren(); i++)
        {
            max = std::max( getMaxTaxonAge( node.getChild(i) ), max );
        }

        return max;
    }
}


double FossilizedBirthDeathSpeciationProcess::lnProbTreeShape(void) const
{
    // from an unlabeled oriented tree to a labeled one

    return - RbMath::lnFactorial( value->getNumberOfExtantTips() );
}


/**
 * Compute the probability of survival if the process starts with one species at time start and ends at time end.
 *
 * \param[in]    start      Start time of the process.
 * \param[in]    end        End/stopping time of the process.
 *
 * \return Probability of survival.
 */
double FossilizedBirthDeathSpeciationProcess::pSurvival(double start, double end) const
{
    return AbstractFossilizedBirthDeathRangeProcess::p(findIndex(start), start, true);
}


/**
 * q_i(t)
 */
double FossilizedBirthDeathSpeciationProcess::q( size_t i, double t, bool tilde ) const
{

    if ( t == times[i] ) return 0.0;

    // get the parameters
    double b = birth[i];
    double d = death[i];
    double f = fossil[i];
    double r = (i == 0 ? homogeneous_rho->getValue() : 0.0);
    double ti = times[i];

    double diff = b - d - f;
    double dt   = t - ti;

    double A = sqrt( diff*diff + 4.0*b*f);
    double B = ( (1.0 - 2.0*(1.0-r)*p_i[i] )*b + d + f ) / A;

    double ln_e = -A*dt;

    double tmp = (1.0 + B) + exp(ln_e)*(1.0 - B);

    double q = log(4.0) + ln_e - 2.0*log(tmp);

    if (tilde)
    {
        q = 0.5 * (q - (b+d+f)*dt);

        double a = anagenetic[i];
        double s = symmetric[i];

        q = - a * dt - s * (b + d + f) * dt + (1.0 - s) * q;
    }

    return q;
}


/**
 *
 */
void FossilizedBirthDeathSpeciationProcess::redrawValue(SimulationCondition c)
{
    redrawValue();
}


/** Set the tree, adopting its annotations where it carries them and drawing the rest. */
void FossilizedBirthDeathSpeciationProcess::setValue(Tree *v, bool force)
{
    AbstractBirthDeathProcess::setValue(v, force);
    cached_events_level = 0;

    // the per-taxon terms are read by tip index, and the tip moves read each tip's status, which a Newick drops
    const std::vector<TopologyNode*> nodes = this->getValue().getNodes();
    for (size_t i = 0; i < this->getValue().getNumberOfTips(); i++)
    {
        size_t j = find(taxa.begin(), taxa.end(), nodes[i]->getTaxon()) - taxa.begin();
        nodes[i]->setIndex(j);
        if ( j < taxa.size() ) nodes[i]->setTaxon( taxa[j] );
    }
    this->getValue().orderNodesByIndex();

    // a sampled ancestor at a symmetric speciation sits beside a zero-length branch, which a Newick cannot tell
    for (size_t i = 0; i < this->getValue().getNumberOfTips(); i++)
    {
        TopologyNode &tip = this->getValue().getNode(i);
        std::optional<std::string> sa = tip.eraseNodeParameter( "sa" );
        if ( sa.has_value() == true && *sa == "1" ) tip.setSampledAncestor( true );
    }

    // flags before ages: the flags name each origination time
    if ( adoptSpeciesLabels() == false )
    {
        normalizeContinuationFlags();
    }

    prepareProbComputation();
    updateRanges();

    // a tree from outside carries no status
    for (size_t i = 0; i < taxa.size(); ++i)
    {
        ranges[i].survived = ( ranges[i].death == getPresent() );
    }

    // redraw against this tree
    if ( adoptAppearances() == false )
    {
        for (size_t i = 0; i < taxa.size(); ++i)
        {
            initializeFirstLast(i);
        }
    }
}


/** The tree annotated with the species labels and appearances a plain Newick drops. */
std::string FossilizedBirthDeathSpeciationProcess::getHiddenStateString( void ) const
{
    Tree t = *this->value;

    size_t unsampled_count = 0;
    labelSpecies( t.getRoot(), rangeLabel( t.getRoot(), unsampled_count ), unsampled_count );

    // a successor's species label is its own either way, so its mode needs saying
    for (TopologyNode *node : t.getNodes())
    {
        if ( isModeSuccessor( *node ) == true ) node->setNodeParameter( "mode", node->continuesParentSpecies() ? "A" : "S" );
    }

    for (size_t i = 0; i < taxa.size(); ++i)
    {
        TopologyNode &tip = t.getNode(i);

        std::stringstream first, last;
        first.precision( 17 );
        last.precision( 17 );
        first << ranges[i].first;
        last  << ranges[i].last;

        tip.setNodeParameter( "FAD", first.str() );
        tip.setNodeParameter( "LAD", last.str() );
        // singleton status is latent only under presence reporting; otherwise the occurrence count sets it
        if ( ranges[i].singleton != countSingleton( i ) ) tip.setNodeParameter( "singleton", ranges[i].singleton ? "1" : "0" );
        tip.setNodeParameter( "survived", ranges[i].survived ? "1" : "0" );
        // a Newick reads a zero-length tip as a sampled ancestor only beside a positive branch
        if ( tip.isSampledAncestorTip() == true ) tip.setNodeParameter( "sa", "1" );
        if ( num_components > 1 ) tip.setNodeParameter( "component", StringUtilities::to_string( ranges[i].component ) );
    }

    // exact ages, since rounded branch lengths move a survivor off the present
    for (TopologyNode *node : t.getNodes())
    {
        std::stringstream age;
        age.precision( 17 );
        age << node->getAge();
        node->setNodeParameter( "age", age.str() );
    }

    return t.getNewickRepresentation( false );
}


void FossilizedBirthDeathSpeciationProcess::setHiddenStateFromString( const std::string &s )
{
    // read back as a time tree
    Tree *t = this->value->clone();
    t->initFromString( s );

    // exact ages, since a collapsed symmetric node must sit at its sampled ancestor's age
    for (TopologyNode *node : t->getNodes())
    {
        std::optional<std::string> age = node->eraseNodeParameter( "age" );
        if ( age.has_value() == false ) throw RbException("A dnFBDSP checkpoint must give every node's age.");
        node->setAge( atof( age->c_str() ), false );
    }

    setValue( t );
}


/** The range each node belongs to, stamped top down. */
void FossilizedBirthDeathSpeciationProcess::labelSpecies( TopologyNode &node, const std::string &s, size_t &unsampled_count ) const
{
    node.setNodeParameter( "sp", s );

    std::vector<TopologyNode*> children = node.getChildren();
    for (size_t c = 0; c < children.size(); c++)
    {
        // a child that does not continue begins its own range at this node
        std::string cs = continuesRange( *children[c] ) ? s : rangeLabel( *children[c], unsampled_count );

        labelSpecies( *children[c], cs, unsampled_count );
    }
}


/** A range no taxon samples still needs a label of its own, so the flags along it survive a checkpoint. */
std::string FossilizedBirthDeathSpeciationProcess::rangeLabel( const TopologyNode &node, size_t &unsampled_count ) const
{
    int k = continuingSpecies( node );
    if ( k >= 0 ) return StringUtilities::to_string( k );

    return "u" + StringUtilities::to_string( unsampled_count++ );
}


/** A child continues its parent's range exactly when the two labels match. */
bool FossilizedBirthDeathSpeciationProcess::adoptSpeciesLabels( void )
{
    Tree &t = this->getValue();

    size_t n = t.getNumberOfNodes();

    std::vector<std::string> label( n );
    bool annotated = false;

    for (size_t i = 0; i < n; ++i)
    {
        const TopologyNode &node = t.getNode(i);

        std::optional<std::string> s = node.getNodeParameter( "sp" );
        if ( s.has_value() == false ) continue;

        label[ node.getIndex() ] = *s;
        annotated = true;
    }

    if ( annotated == false ) return false;

    for (size_t i = 0; i < n; ++i)
    {
        TopologyNode &node = t.getNode(i);

        node.eraseNodeParameter( "sp" );

        if ( node.isRoot() == true ) continue;

        const std::string &mine   = label[ node.getIndex() ];
        const std::string &parent = label[ node.getParent().getIndex() ];

        node.setContinuesParentSpecies( mine.empty() == false && mine == parent );
    }

    // a successor's mode, anagenesis where none is given and the rates allow it
    for (size_t i = 0; i < n; ++i)
    {
        TopologyNode &node = t.getNode(i);
        std::optional<std::string> m = node.eraseNodeParameter( "mode" );
        if ( isModeSuccessor( node ) == false ) continue;
        node.setContinuesParentSpecies( m.has_value() ? *m == "A" : defaultSuccessorAnagenetic() );
    }

    return true;
}


/** Whether range i is a singleton by its occurrence count, as it is unless presence reporting samples it. */
bool FossilizedBirthDeathSpeciationProcess::countSingleton( size_t i ) const
{
    size_t count = 0;
    for (const std::pair<const TimeInterval, size_t> &o : taxa[i].getOccurrences()) count += o.second;
    return count < 2;
}


bool FossilizedBirthDeathSpeciationProcess::adoptAppearances( void )
{
    Tree &t = this->getValue();

    bool annotated = false;

    for (size_t i = 0; i < taxa.size(); ++i)
    {
        TopologyNode &tip = t.getNode(i);

        std::optional<std::string> f = tip.eraseNodeParameter( "FAD" );
        std::optional<std::string> l = tip.eraseNodeParameter( "LAD" );

        if ( f.has_value() == true )
        {
            ranges[i].first = atof( f->c_str() );
            annotated = true;
        }
        if ( l.has_value() == true )
        {
            ranges[i].last = atof( l->c_str() );
            annotated = true;
        }

        std::optional<std::string> single   = tip.eraseNodeParameter( "singleton" );
        std::optional<std::string> survived = tip.eraseNodeParameter( "survived" );

        ranges[i].singleton = ( single.has_value() == true ? *single == "1" : countSingleton( i ) );
        if ( survived.has_value() == true ) ranges[i].survived  = ( *survived == "1" );

        std::optional<std::string> component = tip.eraseNodeParameter( "component" );
        if ( component.has_value() == true && size_t( atoi( component->c_str() ) ) < num_components ) ranges[i].component = atoi( component->c_str() );

        if ( annotated == true ) dirty_taxa[i] = true;
    }

    return annotated;
}


void FossilizedBirthDeathSpeciationProcess::setOccurrences( const std::vector<Taxon> &t )
{
    AbstractFossilizedBirthDeathRangeProcess::setOccurrences( t );

    // keep the tree's own taxon copy in step
    AbstractRootedTreeDistribution::taxa = taxa;
}


void FossilizedBirthDeathSpeciationProcess::setMcmcMode(bool tf)
{
    AbstractBirthDeathProcess::setMcmcMode(tf);
    if ( tf == true ) warnIfNoResampleMove();
    if ( tf == true ) warnIfNoReportingNode();
    if ( tf == true ) warnIfInertNTaxaCondition();
    if ( tf == true ) warnIfNoSurvivalMove();
    if ( tf == true ) warnIfNoSingletonMove();
}


void FossilizedBirthDeathSpeciationProcess::redrawValue(void)
{
    // draw the ranges, then a uniform budding topology on them
    drawRanges();

    RandomNumberGenerator* rng = GLOBAL_RNG;
    size_t n = taxa.size();

    // the extinction must sit below the oldest appearance
    double present = times.front();
    for (size_t i = 0; i < n; ++i)
    {
        if ( ranges[i].death > ranges[i].first ) ranges[i].death = rng->uniform01() * (ranges[i].first - present) + present;
    }

    // one tip per taxon at its drawn end age; top[i] tracks the current top of lineage i's subtree
    std::vector<TopologyNode*> top(n);
    for (size_t i = 0; i < n; ++i)
    {
        TopologyNode* tip = new TopologyNode( taxa[i], i );
        tip->setTipAgeUnconstrained( true );
        tip->setAge( ranges[i].death );
        top[i] = tip;
    }

    // births oldest first, each inside a living lineage already placed
    double origin = getOriginAge();

    std::vector<size_t> order(n);
    for (size_t i = 0; i < n; ++i) order[i] = i;
    std::sort( order.begin(), order.end(), [this](size_t a, size_t b){ return ranges[a].first > ranges[b].first; } );

    size_t root_lineage = order[0];
    std::vector<size_t> parent(n, n);   // parent[k] = lineage k buds off; n marks the origin
    ranges[root_lineage].divergence = origin;
    for (size_t idx = 1; idx < n; ++idx)
    {
        size_t k = order[idx];

        // already-placed lineages that are old enough to bud k
        std::vector<size_t> cand;
        for (size_t p = 0; p < idx; ++p)
        {
            if ( ranges[order[p]].divergence > ranges[k].first ) cand.push_back( order[p] );
        }

        // none old enough: bud from the origin and let the density reject it
        size_t j = root_lineage;
        if ( cand.empty() == false )
        {
            size_t pick = size_t( rng->uniform01() * cand.size() );
            if ( pick >= cand.size() ) pick = cand.size() - 1;

            j = cand[pick];
        }

        // birth in (max(first_k, d_j), b_j)
        double lo = std::max( ranges[k].first, ranges[j].death );
        if ( lo >= ranges[j].divergence ) lo = 0.0;

        ranges[k].divergence = rng->uniform01()*(ranges[j].divergence - lo) + lo;
        parent[k] = j;
    }

    // build the tree youngest birth first, so each lineage's subtree is complete before it attaches
    std::vector<size_t> byBirth(n);
    for (size_t i = 0; i < n; ++i) byBirth[i] = i;
    std::sort( byBirth.begin(), byBirth.end(), [this](size_t a, size_t b){ return ranges[a].divergence < ranges[b].divergence; } );

    for (size_t idx = 0; idx < n; ++idx)
    {
        size_t k = byBirth[idx];
        if ( parent[k] == n ) continue;   // the origin lineage stays the root
        size_t j = parent[k];

        // budding node at b_k: j continues and k is new, recorded on the children
        TopologyNode* node = new TopologyNode();
        node->setAge( ranges[k].divergence );
        node->addChild( top[j] );
        node->addChild( top[k] );
        top[j]->setParent( node );
        top[k]->setParent( node );
        top[j]->setContinuesParentSpecies( true );
        top[k]->setContinuesParentSpecies( false );
        top[j] = node;
    }

    Tree *psi = new Tree();
    psi->setRooted( true );
    psi->setRoot( top[root_lineage], true );

    delete this->value;
    this->value = psi;

    // index tips by taxon order, as the likelihood expects
    const std::vector<TopologyNode*> nodes = this->getValue().getNodes();
    for( size_t i = 0; i < this->getValue().getNumberOfTips(); i++)
    {
        size_t j = find(taxa.begin(), taxa.end(), nodes[i]->getTaxon()) - taxa.begin();
        nodes[i]->setIndex(j);
    }
    this->getValue().orderNodesByIndex();

    // exactly one continuation per node
    normalizeContinuationFlags();
}



/** beta in the interval containing age, read from the parameter. */
double FossilizedBirthDeathSpeciationProcess::symmetricAt( double age ) const
{
    if ( homogeneous_beta != NULL )
    {
        return homogeneous_beta->getValue();
    }

    const RbVector<double>& probs = heterogeneous_beta->getValue();
    size_t i = findIndex( age );

    return i < probs.size() ? probs[i] : 0.0;
}


bool FossilizedBirthDeathSpeciationProcess::hasSymmetricSpeciation( void ) const
{
    if ( homogeneous_beta != NULL )
    {
        return homogeneous_beta->getValue() > 0.0;
    }
    if ( heterogeneous_beta != NULL )
    {
        const RbVector<double>& probs = heterogeneous_beta->getValue();
        for (size_t i = 0; i < probs.size(); ++i)
        {
            if ( probs[i] > 0.0 ) return true;
        }
    }

    return false;
}


/** int_{t0}^{t1} lambda_a dt for each span of ages t0 <= t1. */
std::vector<double> FossilizedBirthDeathSpeciationProcess::anageneticIntegrals(const std::vector<std::pair<double,double> > &spans) const
{
    prepareProbComputation();

    std::vector<double> out;
    out.reserve( spans.size() );
    for (const std::pair<double,double> &span : spans)
    {
        double t0 = span.first, t1 = span.second;
        double total = 0.0;
        for (size_t i = findIndex(t0); t0 < t1; i++)
        {
            double end = ( i + 1 < num_intervals ) ? std::min( t1, times[i+1] ) : t1;
            total += anagenetic[i] * (end - t0);
            t0 = end;
        }
        out.push_back( total );
    }

    return out;
}


/**
 * Over each span of ages t0 <= t1, the Poisson means of the hidden events that change an observed lineage's
 * species: buddings whose parent side leaves no sample, lambda (1-beta) p; symmetric speciations whose other
 * daughter leaves none, 2 lambda beta p; and anageneses, lambda_a.
 */
void FossilizedBirthDeathSpeciationProcess::hiddenEventIntegrals(const std::vector<std::pair<double,double> > &spans, std::vector<double> &budding, std::vector<double> &bifurcation, std::vector<double> &anagenesis) const
{
    anagenesis = anageneticIntegrals( spans );
    budding.assign( spans.size(), 0.0 );
    bifurcation.assign( spans.size(), 0.0 );

    for (size_t k = 0; k < spans.size(); k++)
    {
        double t0 = spans[k].first, t1 = spans[k].second;

        // int lambda p dt, split by each interval's beta
        for (size_t i = findIndex(t0); t0 < t1; i++)
        {
            double end = ( i + 1 < num_intervals ) ? std::min( t1, times[i+1] ) : t1;
            double lp  = 0.5 * ( q(i, end) - q(i, t0) + (birth[i] + death[i] + fossil[i]) * (end - t0) );
            budding[k]     += ( 1.0 - symmetric[i] ) * lp;
            bifurcation[k] += 2.0 * symmetric[i] * lp;
            t0 = end;
        }
    }
}


/** The stem a branch lies on: from its range's top down to its first appearance, or to the symmetric node ending it. */
FossilizedBirthDeathSpeciationProcess::Stem FossilizedBirthDeathSpeciationProcess::stemOf( const TopologyNode &node, long s ) const
{
    Stem st;
    if ( s >= 0 )
    {
        st.top = ranges[s].divergence;
        st.end = ranges[s].first;
        st.end_component = long( ranges[s].component );
        mixtureStemSource( size_t( s ), st.source );
        return st;
    }

    // a range no taxon samples starts fresh and is summed over the component it ends in
    const TopologyNode *n = &node;
    while ( n->isRoot() == false && continuesRange( *n ) == true ) n = &n->getParent();
    st.top = ( n->isRoot() == true ? getOriginAge() : n->getParent().getAge() );

    n = &node;
    for (bool down = true; down == true; )
    {
        down = false;
        for (size_t c = 0; c < n->getNumberOfChildren(); c++)
        {
            if ( continuesRange( n->getChild(c) ) == true ) { n = &n->getChild(c); down = true; break; }
        }
    }
    st.end = n->getAge();
    st.end_component = -1;
    st.source = mixture_weights;
    return st;
}


/**
 * Under the psi mixture, the hidden events that change a lineage's species on each span of its stem, by forward-backward
 * over the stem's components: a(t) from the stem's top, g(t) to its end, and E[N] = int a(t) r(t) (w . g(t)) dt / (a . g).
 * A budding hands the lineage on from a parent left in its own component k, at lambda (1-beta) p_k.
 */
void FossilizedBirthDeathSpeciationProcess::mixtureHiddenEvents(const std::vector<std::pair<double,double> > &spans, const std::vector<Stem> &stems, std::vector<double> &budding, std::vector<double> &bifurcation, std::vector<double> &anagenesis) const
{
    prepareProbComputation();

    size_t H = num_components;
    std::vector<size_t> all( H );
    for (size_t h = 0; h < H; h++) all[h] = h;

    budding.assign( spans.size(), 0.0 );
    bifurcation.assign( spans.size(), 0.0 );
    anagenesis.assign( spans.size(), 0.0 );

    double nu_max = 0.0;
    for (size_t c = 0; c < H; c++) nu_max = std::max( nu_max, mixture_rates[c] );

    std::vector<double> Q, Q_scale, p, M;
    std::vector<std::vector<double> > steps, fwd, bwd;
    double lp;

    // M[c*H + e]: from component e at age t down to c at age s <= t, up to a common factor
    auto stepMatrix = [&]( double s, double t, std::vector<double> &m )
    {
        mixtureSensitivity( s, t, all, Q, Q_scale );
        double top = *std::max_element( Q_scale.begin(), Q_scale.end() );
        m.resize( H*H );
        for (size_t c = 0; c < H; c++) for (size_t e = 0; e < H; e++) m[c*H + e] = Q[c*H + e] * exp( Q_scale[c] - top );
    };
    auto normalize = []( std::vector<double> &v )
    {
        double mx = *std::max_element( v.begin(), v.end() );
        if ( mx > 0.0 ) for (double &x : v) x /= mx;
    };

    const size_t panels = 32;
    for (size_t k = 0; k < spans.size(); k++)
    {
        const Stem &st = stems[k];
        std::vector<double> end_vec( H, st.end_component >= 0 ? 0.0 : 1.0 );
        if ( st.end_component >= 0 ) end_vec[ st.end_component ] = 1.0;

        double t0 = spans[k].first, t1 = spans[k].second;
        for (size_t i = findIndex(t0); t0 < t1; i++)
        {
            double t_end = ( i + 1 < num_intervals ) ? std::min( t1, times[i+1] ) : t1;
            double len = t_end - t0;

            // a stiff component's transient at either end of the stem: a mesh graded from 0.02 of its decay time
            double rate = birth[i] + death[i] + fossil[i] * nu_max + anagenetic[i];
            double h_max = len / panels, h = std::min( h_max, 0.02 / rate );
            std::vector<double> left( 1, 0.0 );
            while ( left.back() + h < 0.5 * len && h < h_max ) { left.push_back( left.back() + h ); h *= 2.0; }
            std::vector<double> x = left;
            double rest = len - 2.0 * left.back();
            size_t n_mid = std::max( size_t( 1 ), size_t( std::ceil( rest / h_max ) ) );
            for (size_t j = 1; j <= n_mid; j++) x.push_back( left.back() + rest * j / n_mid );
            for (size_t j = left.size() - 1; j-- > 0; ) x.push_back( len - left[j] );

            // the ages, youngest first, with each cell's midpoint for Simpson
            std::vector<double> y;
            for (size_t j = 0; j + 1 < x.size(); j++) { y.push_back( t0 + x[j] ); y.push_back( t0 + 0.5 * ( x[j] + x[j+1] ) ); }
            y.push_back( t_end );
            size_t N = y.size();

            steps.resize( N - 1 );
            for (size_t j = 0; j + 1 < N; j++) stepMatrix( y[j], y[j+1], steps[j] );

            // a from the stem's top down, then step by step to the youngest age
            fwd.assign( N, std::vector<double>( H, 0.0 ) );
            if ( y[N-1] < st.top )
            {
                stepMatrix( y[N-1], st.top, M );
                for (size_t c = 0; c < H; c++) for (size_t e = 0; e < H; e++) fwd[N-1][c] += st.source[e] * M[c*H + e];
            }
            else fwd[N-1] = st.source;
            normalize( fwd[N-1] );
            for (size_t j = N - 1; j-- > 0; )
            {
                for (size_t c = 0; c < H; c++) for (size_t e = 0; e < H; e++) fwd[j][c] += fwd[j+1][e] * steps[j][c*H + e];
                normalize( fwd[j] );
            }

            // g from the stem's end up, then step by step to the oldest age
            bwd.assign( N, std::vector<double>( H, 0.0 ) );
            if ( y[0] > st.end )
            {
                stepMatrix( st.end, y[0], M );
                for (size_t e = 0; e < H; e++) for (size_t c = 0; c < H; c++) bwd[0][e] += M[c*H + e] * end_vec[c];
            }
            else bwd[0] = end_vec;
            normalize( bwd[0] );
            for (size_t j = 1; j < N; j++)
            {
                for (size_t e = 0; e < H; e++) for (size_t c = 0; c < H; c++) bwd[j][e] += steps[j-1][c*H + e] * bwd[j-1][c];
                normalize( bwd[j] );
            }

            // the integrand by mode; a . g normalizes each point, so the vectors' scales cancel
            std::vector<std::array<double,3> > f( N );
            for (size_t j = 0; j < N; j++)
            {
                mixtureAt( y[j], p, lp );
                double pbar = 0.0, ag = 0.0, a_sum = 0.0, a_p = 0.0, wg = 0.0;
                for (size_t c = 0; c < H; c++)
                {
                    pbar  += mixture_weights[c] * p[c];
                    ag    += fwd[j][c] * bwd[j][c];
                    a_sum += fwd[j][c];
                    a_p   += fwd[j][c] * p[c];
                    wg    += mixture_weights[c] * bwd[j][c];
                }
                double scale = ( ag > 0.0 ? wg / ag : 0.0 );
                f[j] = { scale * birth[i] * ( 1.0 - symmetric[i] ) * a_p,
                         scale * 2.0 * birth[i] * symmetric[i] * pbar * a_sum,
                         scale * anagenetic[i] * a_sum };
            }
            for (size_t j = 0; j + 2 < N; j += 2)
            {
                double wgt = ( y[j+2] - y[j] ) / 6.0;
                budding[k]     += wgt * ( f[j][0] + 4.0 * f[j+1][0] + f[j+2][0] );
                bifurcation[k] += wgt * ( f[j][1] + 4.0 * f[j+1][1] + f[j+2][1] );
                anagenesis[k]  += wgt * ( f[j][2] + 4.0 * f[j+1][2] + f[j+2][2] );
            }
            t0 = t_end;
        }
    }
}


/** Each interval's beta for the psi mixture's solution, empty where there is no symmetric speciation. */
void FossilizedBirthDeathSpeciationProcess::mixtureSymmetric( std::vector<double> &beta ) const
{
    beta.clear();
    if ( hasSymmetricSpeciation() == true ) beta = symmetric;
}


/** Each interval's anagenetic rate, which moves a lineage to a fresh component. */
void FossilizedBirthDeathSpeciationProcess::mixtureAnagenetic( std::vector<double> &rate ) const
{
    rate.clear();
    if ( hasAnagenesis() == true ) rate = anagenetic;
}


/**
 * An unsampled range's lineage from its top to the symmetric node that ends it, starting in source: the stem matrix
 * summed over the component it ends in, since nothing is observed of it.
 */
double FossilizedBirthDeathSpeciationProcess::unsampledMixtureLnProb( const UnsampledRange &r, const std::vector<double> &source ) const
{
    size_t H = num_components;
    std::vector<size_t> columns( H );
    for (size_t c = 0; c < H; c++) columns[c] = c;
    std::vector<double> Q, Q_scale;
    mixtureSensitivity( r.end, r.top, columns, Q, Q_scale );

    // sum over the end component c, each column on its own scale
    double top_scale = -std::numeric_limits<double>::infinity();
    for (size_t c = 0; c < H; c++) top_scale = std::max( top_scale, Q_scale[c] );
    double total = 0.0;
    for (size_t c = 0; c < H; c++)
    {
        double s = 0.0;
        for (size_t e = 0; e < H; e++) s += source[e] * Q[c*H + e];
        total += s * exp( Q_scale[c] - top_scale );
    }
    return log( total ) + top_scale;
}


/** The range parameters, with the anagenetic and symmetric rates the base class does not hold. */
std::vector<const DagNode*> FossilizedBirthDeathSpeciationProcess::getRateParameters( void ) const
{
    std::vector<const DagNode*> out;
    for (const DagNode *n : range_parameters) if ( n != NULL ) out.push_back( n );
    for (const DagNode *n : { (const DagNode*) homogeneous_lambda_a, (const DagNode*) heterogeneous_lambda_a,
                              (const DagNode*) homogeneous_beta, (const DagNode*) heterogeneous_beta }) if ( n != NULL ) out.push_back( n );
    return out;
}


/** The species whose lineage runs through this node: follow the continuing children to a tip. */
long FossilizedBirthDeathSpeciationProcess::speciesBelow( const TopologyNode &node ) const
{
    const TopologyNode *n = &node;
    while ( n->isTip() == false )
    {
        const TopologyNode *next = NULL;
        for (size_t c = 0; c < n->getNumberOfChildren(); c++)
        {
            // a sampled ancestor tip names the range, past which only a symmetric node at its age continues
            if ( continuesRange( n->getChild(c) ) && ( next == NULL || n->getChild(c).isSampledAncestorTip() ) ) next = &n->getChild(c);
        }
        if ( next == NULL ) return -1;
        n = next;
    }
    return long( n->getIndex() );
}


/**
 * By node index, the Poisson means of the hidden speciations on each branch and the observed event that began
 * its species. A hidden event is one that changes the observed lineage's species and leaves the rest unsampled;
 * on a sampled range they are counted only above its tau_1, since a sampled species does not change, and on a
 * range no taxon samples along the whole branch. A branch whose node does not continue its parent's species
 * begins with the observed event that made it: an anagenesis where the continuing sibling is a sampled
 * ancestor, a symmetric speciation where no sibling continues, and otherwise a budding.
 */
void FossilizedBirthDeathSpeciationProcess::branchEvents( BranchEvents &e, bool hidden ) const
{
    // six accessors read one tree state, so compute it once
    if ( cached_events_level >= ( hidden ? 2 : 1 ) )
    {
        e = cached_events;
        return;
    }

    const_cast<FossilizedBirthDeathSpeciationProcess *>( this )->updateRanges();

    const Tree &t = *value;
    size_t nodes = t.getNumberOfNodes();
    for (std::vector<double> *v : { &e.hidden_budding, &e.hidden_bifurcation, &e.hidden_anagenetic,
                                    &e.observed_budding, &e.observed_bifurcation, &e.observed_anagenetic })
    {
        v->assign( nodes, 0.0 );
    }

    std::vector<size_t> at;
    std::vector<std::pair<double,double> > spans;
    std::vector<Stem> stems;
    for (const TopologyNode *node : t.getNodes())
    {
        if ( node->isRoot() ) continue;
        size_t idx = node->getIndex();
        const TopologyNode &parent = node->getParent();

        if ( isModeSuccessor( *node ) == true )
        {
            ( node->continuesParentSpecies() ? e.observed_anagenetic : e.observed_bifurcation )[idx] = 1.0;
        }
        else if ( node->continuesParentSpecies() == false )
        {
            bool sa_sibling = false, any_continues = false;
            for (size_t c = 0; c < parent.getNumberOfChildren(); c++)
            {
                const TopologyNode &sib = parent.getChild(c);
                if ( sib.continuesParentSpecies() == false ) continue;
                any_continues = true;
                if ( sib.isSampledAncestorTip() ) sa_sibling = true;
            }
            if ( sa_sibling )               e.observed_anagenetic[idx]  = 1.0;
            else if ( any_continues == false ) e.observed_bifurcation[idx] = 1.0;
            else                            e.observed_budding[idx]     = 1.0;
        }

        long s = speciesBelow( *node );
        if ( s >= long( ranges.size() ) ) continue;
        double lo = ( s < 0 ? node->getAge() : std::max( node->getAge(), ranges[s].first ) );
        double hi = parent.getAge();

        if ( hi > lo )
        {
            at.push_back( idx );
            spans.push_back( std::make_pair( lo, hi ) );
            if ( hidden == true && num_components > 1 ) stems.push_back( stemOf( *node, s ) );
        }
    }

    if ( hidden == false )
    {
        cached_events = e;
        cached_events_level = 1;
        return;
    }

    std::vector<double> hb, hs, ha;
    if ( num_components > 1 ) mixtureHiddenEvents( spans, stems, hb, hs, ha );
    else                      hiddenEventIntegrals( spans, hb, hs, ha );
    for (size_t k = 0; k < at.size(); k++)
    {
        e.hidden_budding[ at[k] ]     = hb[k];
        e.hidden_bifurcation[ at[k] ] = hs[k];
        e.hidden_anagenetic[ at[k] ]  = ha[k];
    }

    cached_events = e;
    cached_events_level = 2;
}


bool FossilizedBirthDeathSpeciationProcess::isExtended( void ) const { return true; }


/**
 * The sampled tree an extended tree projects to. Each extinct species is cut back to its last appearance: its tip
 * moves there, or, where its lineage leads on to something below it, a sampled ancestor sits there and the lineage
 * carries on to the youngest node it leads to, the node at the old tip suppressed.
 */
Tree FossilizedBirthDeathSpeciationProcess::sampledTree( void ) const
{
    Tree t = *value;

    size_t unsampled_count = 0;
    labelSpecies( t.getRoot(), rangeLabel( t.getRoot(), unsampled_count ), unsampled_count );

    std::vector<TopologyNode*> tips;
    for (size_t i = 0; i < taxa.size(); i++) tips.push_back( &t.getNode(i) );
    TopologyNode *root = &t.getRoot();

    for (size_t i = 0; i < taxa.size(); i++)
    {
        TopologyNode *tip = tips[i];
        if ( ranges[i].survived == true || tip == root ) continue;

        double tk = ranges[i].last;

        // the edge holding the last appearance, on the species' own lineage; above the root, on its stem
        TopologyNode *n = tip;
        while ( n != root && n->getParent().getAge() < tk ) n = &n->getParent();

        if ( n == tip )
        {
            tip->setAge( tk );
            continue;
        }

        // below its sample the lineage leads to a successor: the child it changed into, or a range no taxon samples
        // where it ended at a symmetric speciation. Its path takes the successor's label.
        TopologyNode *p0 = &tip->getParent();
        TopologyNode *c  = ( &p0->getChild(0) == tip ? &p0->getChild(1) : &p0->getChild(0) );
        std::string mine = StringUtilities::to_string( i );
        std::string successor = *c->getNodeParameter( "sp" );
        if ( successor == mine ) successor = "u" + StringUtilities::to_string( unsampled_count++ );
        for (TopologyNode *k = n; k != p0 && k != NULL; )
        {
            k->setNodeParameter( "sp", successor );
            TopologyNode *next = NULL;
            for (size_t j = 0; j < k->getNumberOfChildren(); j++)
            {
                if ( *k->getChild(j).getNodeParameter( "sp" ) == mine ) next = &k->getChild(j);
            }
            k = next;
        }
        c->setNodeParameter( "sp", successor );

        // drop the old tip, and the node it hung from
        p0->removeChild( tip );
        p0->removeChild( c );
        if ( p0 == root )
        {
            p0->setParent( NULL, false );
            c->setParent( NULL, false );
            root = c;
        }
        else
        {
            TopologyNode *gp = &p0->getParent();
            gp->addChild( c );
            c->setParent( gp );
        }
        if ( n == p0 ) n = c;
        // which takes it off its parent; setRoot deletes an old root itself
        if ( p0 != &t.getRoot() ) delete p0;

        // and a sampled ancestor at the last appearance, above what the lineage leads to
        TopologyNode *a = new TopologyNode();
        a->setTree( &t );
        if ( n == root )
        {
            root = a;
        }
        else
        {
            TopologyNode *q = &n->getParent();
            q->removeChild( n );
            q->addChild( a );
            a->setParent( q );
        }
        a->setAge( tk );
        a->setNodeParameter( "sp", StringUtilities::to_string( i ) );
        a->addChild( tip );
        tip->setParent( a );
        a->addChild( n );
        n->setParent( a );
        tip->setAge( tk );
        tip->setSampledAncestor( true );
    }

    // tips by taxon, as the process reads them, then the internal nodes; the tree's own node list is stale
    size_t next = taxa.size();
    std::function<void(TopologyNode&)> index = [&]( TopologyNode &node ) {
        if ( node.isTip() == true ) node.setIndex( std::find( taxa.begin(), taxa.end(), node.getTaxon() ) - taxa.begin() );
        else                        node.setIndex( next++ );
        for (size_t c = 0; c < node.getNumberOfChildren(); c++) index( node.getChild(c) );
    };
    index( *root );
    t.setRoot( root, false );

    return t;
}


void FossilizedBirthDeathSpeciationProcess::executeMethod( const std::string &n, const std::vector<const DagNode *> &args, Tree &rv ) const
{
    if ( n != "getSampledTree" ) throw RbException( "dnFBDSP has no tree method " + n + "." );

    const_cast<FossilizedBirthDeathSpeciationProcess *>( this )->updateRanges();
    rv = sampledTree();
}


/** The per-branch speciation events, then the range methods of the base class. */
void FossilizedBirthDeathSpeciationProcess::executeMethod( const std::string &n, const std::vector<const DagNode *> &args, RbVector<double> &rv ) const
{
    if ( n == "getHiddenSpeciation" || n == "getObservedSpeciation" )
    {
        // the observed events need no hidden integral, so they work under any mixture
        bool hidden = ( n == "getHiddenSpeciation" );
        const std::string &mode = static_cast<const TypedDagNode<std::string> *>( args[0] )->getValue();
        BranchEvents e;
        branchEvents( e, hidden );

        bool bud = ( mode == "all" || mode == "budding" || mode == "cladogenesis" );
        bool bif = ( mode == "all" || mode == "bifurcation" || mode == "cladogenesis" );
        bool ana = ( mode == "all" || mode == "anagenesis" );
        const std::vector<double> &b = ( hidden ? e.hidden_budding : e.observed_budding );
        const std::vector<double> &s = ( hidden ? e.hidden_bifurcation : e.observed_bifurcation );
        const std::vector<double> &a = ( hidden ? e.hidden_anagenetic : e.observed_anagenetic );
        std::vector<double> want( b.size(), 0.0 );
        for (size_t i = 0; i < want.size(); i++) want[i] = ( bud ? b[i] : 0.0 ) + ( bif ? s[i] : 0.0 ) + ( ana ? a[i] : 0.0 );

        // indexed as dnPhyloCTMC's branchRates are, with no entry for the root
        if ( value->getRoot().getIndex() + 1 != want.size() )
        {
            throw RbException("dnFBDSP: the root is not the last node, so a per-branch vector cannot drop it.");
        }
        rv.clear();
        for (size_t i = 0; i + 1 < want.size(); i++) rv.push_back( want[i] );
    }
    else
    {
        AbstractFossilizedBirthDeathRangeProcess::executeMethod( n, args, rv );
    }
}


bool FossilizedBirthDeathSpeciationProcess::hasAnagenesis( void ) const
{
    if ( homogeneous_lambda_a != NULL )
    {
        return homogeneous_lambda_a->getValue() > 0.0;
    }
    if ( heterogeneous_lambda_a != NULL )
    {
        const RbVector<double>& rates = heterogeneous_lambda_a->getValue();
        for (size_t i = 0; i < rates.size(); ++i)
        {
            if ( rates[i] > 0.0 ) return true;
        }
    }

    return false;
}


/**
 * The non-tip child of a sampled ancestor's node other than the symmetric node at its age that
 * the range carries on into. Its range begins there either way, and its flag says whether its lineage bifurcated away
 * from the ancestor's: false for a symmetric speciation, true for an anagenesis.
 */
bool FossilizedBirthDeathSpeciationProcess::isModeSuccessor( const TopologyNode &n ) const
{
    if ( n.isRoot() == true || n.isSampledAncestorTip() == true ) return false;

    const TopologyNode &p = n.getParent();
    bool sa = false;
    for (size_t c = 0; c < p.getNumberOfChildren(); c++)
    {
        if ( p.getChild(c).isSampledAncestorTip() == true ) sa = true;
    }

    return sa == true && ( n.isTip() == true || n.getAge() != p.getAge() );
}


/** Whether n carries its parent's range: its flag, except at a successor, which begins its own range. */
bool FossilizedBirthDeathSpeciationProcess::continuesRange( const TopologyNode &n ) const
{
    return n.continuesParentSpecies() == true && isModeSuccessor( n ) == false;
}


/** Whether both successor modes are possible, so its flag is free. */
bool FossilizedBirthDeathSpeciationProcess::hasFreeSuccessorModes( void ) const
{
    return hasAnagenesis() == true && hasSymmetricSpeciation() == true;
}


bool FossilizedBirthDeathSpeciationProcess::defaultSuccessorAnagenetic( void ) const
{
    return hasAnagenesis() == true || hasSymmetricSpeciation() == false;
}


/** Gibbs-redraw the budding topology with the ranges fixed, or return false and leave it. */
bool FossilizedBirthDeathSpeciationProcess::redrawTopology( void )
{
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // the candidate test cannot reach an anagenetic attachment
    if ( hasAnagenesis() == true ) return false;

    // pure budding, so not a Gibbs step under beta > 0
    if ( hasSymmetricSpeciation() == true ) return false;

    updateRanges();

    size_t n = taxa.size();
    if ( n < 2 ) return false;

    size_t root_lineage = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if ( ranges[i].divergence > ranges[root_lineage].divergence ) root_lineage = i;
    }

    // choose every attachment before touching the tree
    std::vector<size_t> parent(n, n);
    for (size_t k = 0; k < n; ++k)
    {
        if ( k == root_lineage ) continue;

        std::vector<size_t> cand;
        for (size_t j = 0; j < n; ++j)
        {
            if ( j != k && ranges[j].divergence > ranges[k].divergence && ranges[j].death < ranges[k].divergence ) cand.push_back( j );
        }

        if ( cand.empty() == true ) return false;

        size_t pick = size_t( rng->uniform01() * cand.size() );
        if ( pick >= cand.size() ) pick = cand.size() - 1;

        parent[k] = cand[pick];
    }

    std::vector<TopologyNode*> top(n);
    for (size_t i = 0; i < n; ++i)
    {
        TopologyNode* tip = new TopologyNode( taxa[i], i );
        tip->setTipAgeUnconstrained( true );
        tip->setAge( ranges[i].death );
        top[i] = tip;
    }

    std::vector<size_t> byBirth(n);
    for (size_t i = 0; i < n; ++i) byBirth[i] = i;
    std::sort( byBirth.begin(), byBirth.end(), [this](size_t a, size_t b){ return ranges[a].divergence < ranges[b].divergence; } );

    for (size_t idx = 0; idx < n; ++idx)
    {
        size_t k = byBirth[idx];
        if ( parent[k] == n ) continue;

        size_t j = parent[k];

        // budding node at b_k: j continues and k is new, recorded on the children
        TopologyNode* node = new TopologyNode();
        node->setAge( ranges[k].divergence );
        node->addChild( top[j] );
        node->addChild( top[k] );
        top[j]->setParent( node );
        top[k]->setParent( node );
        top[j]->setContinuesParentSpecies( true );
        top[k]->setContinuesParentSpecies( false );
        top[j] = node;
    }

    Tree *psi = new Tree();
    psi->setRooted( true );
    psi->setRoot( top[root_lineage], true );

    delete this->value;
    this->value = psi;

    const std::vector<TopologyNode*> nodes = this->getValue().getNodes();
    for( size_t i = 0; i < this->getValue().getNumberOfTips(); i++)
    {
        size_t j = find(taxa.begin(), taxa.end(), nodes[i]->getTaxon()) - taxa.begin();
        nodes[i]->setIndex(j);
    }
    this->getValue().orderNodesByIndex();

    // exactly one continuation per node
    normalizeContinuationFlags();

    return true;
}

/**
 *
 */
void FossilizedBirthDeathSpeciationProcess::simulateClade(std::vector<TopologyNode *> &n, double age, double present, bool alwaysReturn)
{

    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // get the minimum birth age
    std::vector<double> first_occurrences;

    double current_age = RbConstants::Double::inf;
    double minimum_age = 0.0;
    double max_age = getOriginAge();

    for (size_t i = 0; i < n.size(); ++i)
    {
        // make sure the tip age is equal to the last occurrence
        if( n[i]->isTip() )
        {
            bool extinct = n[i]->getTaxon().isExtinct();

            double present = times.front();

            n[i]->setAge( extinct * rng->uniform01() * (n[i]->getTaxon().getAgeRange()->getMin() - present) + present );

            size_t j = find(taxa.begin(), taxa.end(), n[i]->getTaxon()) - taxa.begin();

            double minmax = std::max(ranges[i].first_min, n[i]->getAge());
            ranges[j].first = GLOBAL_RNG->uniform01()*(std::min(max_age, taxa[i].getAgeRange()->getMax()) - minmax) + minmax;
        }

        double first_occurrence = getMaxTaxonAge( *n[i] );

        if( first_occurrence > minimum_age )
        {
            minimum_age = first_occurrence;
        }

        first_occurrences.push_back( first_occurrence );

        if ( current_age > n[i]->getAge() )
        {
            current_age = n[i]->getAge();
        }
    }

    if( minimum_age > max_age )
    {
        throw RbException() << "Tree age is " << max_age << " but minimum fossil origin is " << minimum_age;
    }


    if ( age <= minimum_age )
    {
        age = rng->uniform01() * ( max_age - minimum_age ) + minimum_age;
    }


    std::vector<double> ages;
    while ( n.size() > 2 && current_age < age )
    {

        // get all the nodes with first occurrences younger than the current age
        std::vector<TopologyNode*> active_nodes;
        std::vector<TopologyNode*> active_right_nodes;
        for (size_t i = 0; i < n.size(); ++i)
        {
            if ( current_age >= n[i]->getAge() )
            {
                active_nodes.push_back( n[i] );

                if( current_age >= first_occurrences[i] )
                {
                    active_right_nodes.push_back( n[i] );
                }
            }

        }

        // we need to get the next node age older than the current age
        double next_node_age = age;
        for (size_t i = 0; i < n.size(); ++i)
        {
            if ( current_age < n[i]->getAge() && n[i]->getAge() < next_node_age )
            {
                next_node_age = n[i]->getAge();
            }
            if ( current_age < first_occurrences[i] && first_occurrences[i] < next_node_age )
            {
                next_node_age = first_occurrences[i];
            }

        }

        // only simulate if there are at least two valid/active nodes and one active right node
        if ( active_nodes.size() <= 2 || active_right_nodes.empty() )
        {
            current_age = next_node_age;
        }
        else
        {
            // now we simulate new ages
            double next_sim_age = simulateNextAge(active_nodes.size()-2, age, present, current_age, alwaysReturn);

            if ( next_sim_age < next_node_age )
            {
                // randomly pick two nodes
                size_t index_left = static_cast<size_t>( floor(rng->uniform01()*active_nodes.size()) );
                TopologyNode* left_child = active_nodes[index_left];

                size_t index_right = static_cast<size_t>( floor(rng->uniform01()*active_right_nodes.size()) );
                TopologyNode* right_child = active_right_nodes[index_right];

                while( left_child == right_child )
                {
                    index_left = static_cast<size_t>( floor(rng->uniform01()*active_nodes.size()) );
                    left_child = active_nodes[index_left];

                    index_right = static_cast<size_t>( floor(rng->uniform01()*active_right_nodes.size()) );
                    right_child = active_right_nodes[index_right];
                }

                // erase the nodes also from the origin nodes vector
                std::vector<TopologyNode *>::iterator child_it_left = std::find( n.begin(), n.end(), left_child );
                std::vector<double>::iterator fa_it_left = first_occurrences.begin() + std::distance( n.begin(), child_it_left );
                double fa_left = *fa_it_left;

                first_occurrences.erase(fa_it_left);
                n.erase(child_it_left);

                std::vector<TopologyNode *>::iterator child_it_right = std::find( n.begin(), n.end(), right_child );
                std::vector<double>::iterator fa_it_right = first_occurrences.begin() + std::distance( n.begin(), child_it_right );
                double fa_right = *fa_it_right;

                first_occurrences.erase(fa_it_right);
                n.erase(child_it_right);


                // create a parent for the two
                TopologyNode *parent = new TopologyNode();
                parent->addChild( left_child );
                parent->addChild( right_child );
                left_child->setParent( parent );
                right_child->setParent( parent );
                left_child->setContinuesParentSpecies( true );
                right_child->setContinuesParentSpecies( false );
                parent->setAge( next_sim_age );

                // insert the parent to our list
                n.push_back( parent );
                first_occurrences.push_back( std::max(fa_left, fa_right) );

                current_age = next_sim_age;
                ages.push_back( next_sim_age );
            }
            else
            {
                current_age = next_node_age;
            }

        }

        if ( n.size() > 2 && current_age >= age  )
            throw RbException() << "Unexpected number of taxa (remaining #taxa was " << n.size() << " and age was " << current_age << " with maximum age of " << age << ") in tree simulation";

    }


    if ( n.size() == 2 )
    {

        // pick two nodes
        TopologyNode* left_child = n[0];
        TopologyNode* right_child = n[1];

        // make sure the speciation event is older than the new species first occurrence
        if( first_occurrences[1] > age )
        {
            if( first_occurrences[0] > age )
            {
                throw RbException() << "Cannot simulate clade of age " << age << ", minimum age is " << minimum_age;
            }
            else
            {
                std::swap( left_child, right_child );
            }
        }
        else if( age > first_occurrences[0] )
        {
            if( rng->uniform01() < 0.5 )
            {
                std::swap( left_child, right_child );
            }
        }

        // erase the nodes also from the origin nodes vector
        n.clear();

        // create a parent for the two
        TopologyNode *parent = new TopologyNode();
        parent->addChild( left_child );
        parent->addChild( right_child );
        left_child->setParent( parent );
        right_child->setParent( parent );
        // name a continuing lineage; mvRotateNode explores the choice
        left_child->setContinuesParentSpecies( true );
        right_child->setContinuesParentSpecies( false );
        parent->setAge( age );

        // insert the parent to our list
        n.push_back( parent );
    }
    else
    {
        throw RbException() << "Unexpected number of taxa (" << n.size() << ") in tree simulation";
    }


}

/**
 * @param alwaysReturn whether the simulation can return times which are not valid draws from the distribution (for initial values)
*/
std::vector<double> FossilizedBirthDeathSpeciationProcess::simulateDivergenceTimes(size_t n, double origin, double present, double min, bool alwaysReturn) const
{
    if(!alwaysReturn) throw RbException("Impossible to simulate under a true FossilizedBirthDeathSpeciation process");

    std::vector<double> t(n, 0.0);

    for (size_t i = 0; i < n; ++i)
    {
        t[i] = simulateDivergenceTime(origin, min);
    }

    // finally sort the times
    std::sort(t.begin(), t.end());

    return t;
}

/**
 * Simulate new speciation times.
 */
double FossilizedBirthDeathSpeciationProcess::simulateDivergenceTime(double origin, double present) const
{
    // incorrect placeholder for constant SSBDP
    // direct forward simulation under FBDRP is not feasible

    // Get the rng
    RandomNumberGenerator* rng = GLOBAL_RNG;

    size_t i = findIndex(present);

    // get the parameters
    double age = origin - present;
    double b = birth[i];
    double d = death[i];
    double p_e = i == 0 ? homogeneous_rho->getValue() : 0.0;


    // get a random draw
    double u = rng->uniform01();

    // compute the time for this draw
    // see Hartmann et al. 2010 and Stadler 2011
    double t = 0.0;
    if ( b > d )
    {
        if( p_e > 0.0 )
        {
            t = ( log( ( (b-d) / (1 - (u)*(1-((b-d)*exp((d-b)*age))/(p_e*b+(b*(1-p_e)-d)*exp((d-b)*age) ) ) ) - (b*(1-p_e)-d) ) / (p_e * b) ) )  /  (b-d);
        }
        else
        {
            t = log( 1 - u * (exp(age*(d-b)) - 1) / exp(age*(d-b)) ) / (b-d);
        }
    }
    else
    {
        if( p_e > 0.0 )
        {
            t = ( log( ( (b-d) / (1 - (u)*(1-(b-d)/(p_e*b*exp((b-d)*age)+(b*(1-p_e)-d) ) ) ) - (b*(1-p_e)-d) ) / (p_e * b) ) )  /  (b-d);
        }
        else
        {
            t = log( 1 - u * (1 - exp(age*(b-d)))  ) / (b-d);
        }
    }

    return present + t;
}


void FossilizedBirthDeathSpeciationProcess::normalizeContinuationFlags( void )
{
    // the repair reads beta per interval
    prepareProbComputation();

    normalizeContinuationFlags( getValue().getRoot() );
}


/** The taxon whose range this node belongs to. */
int FossilizedBirthDeathSpeciationProcess::continuingSpecies( const TopologyNode &node ) const
{
    if ( node.isTip() == true ) return int( node.getIndex() );

    // at a sampled ancestor the tip names the range, whichever slot it is in
    const std::vector<TopologyNode*> &children = node.getChildren();
    for (size_t c = 0; c < children.size(); c++)
    {
        if ( continuesRange( *children[c] ) == false ) continue;

        int s = continuingSpecies( *children[c] );
        if ( s >= 0 ) return s;
    }

    return -1;
}


void FossilizedBirthDeathSpeciationProcess::normalizeContinuationFlags( const TopologyNode &node )
{
    if ( node.isTip() == true ) return;

    std::vector<TopologyNode*> children = node.getChildren();

    for (size_t c = 0; c < children.size(); c++)
    {
        normalizeContinuationFlags( *children[c] );
    }

    // a sampled ancestor tip is a sample of this node's species and always continues
    bool sa_tip = false;
    for (size_t c = 0; c < children.size(); c++)
    {
        if ( children[c]->isSampledAncestorTip() == true )
        {
            sa_tip = true;
            children[c]->setContinuesParentSpecies( true );
        }
    }

    size_t n_cont = 0;
    for (size_t c = 0; c < children.size(); c++)
    {
        if ( children[c]->isSampledAncestorTip() == false && continuesRange( *children[c] ) == true ) ++n_cont;
    }

    // flags from outside default to budding, and to anagenesis at a sampled ancestor; only a label or a
    // move makes a symmetric speciation
    if ( sa_tip == true )
    {
        for (size_t c = 0; c < children.size(); c++)
        {
            if ( children[c]->isSampledAncestorTip() == true ) continue;
            children[c]->setContinuesParentSpecies( isModeSuccessor( *children[c] ) && defaultSuccessorAnagenetic() );
        }
        return;
    }

    if ( n_cont != 1 )
    {
        // the lowest admissible index, so presentations agree
        size_t keep = children.size();

        for (size_t c = 0; c < children.size(); c++)
        {
            if ( children[c]->isSampledAncestorTip() == true ) continue;

            bool feasible = true;
            for (size_t o = 0; o < children.size() && feasible == true; o++)
            {
                if ( o == c || children[o]->isSampledAncestorTip() == true ) continue;

                int s = continuingSpecies( *children[o] );
                if ( s < 0 || ranges[s].first_min >= node.getAge() ) feasible = false;
            }

            if ( feasible == false ) continue;
            if ( keep == children.size() || children[c]->getIndex() < children[keep]->getIndex() ) keep = c;
        }

        // nothing fits: take the lowest index and let the density reject it
        for (size_t c = 0; keep == children.size() && c < children.size(); c++)
        {
            if ( children[c]->isSampledAncestorTip() == true ) continue;
            keep = c;
        }

        for (size_t c = 0; c < children.size(); c++)
        {
            if ( children[c]->isSampledAncestorTip() == true ) continue;
            children[c]->setContinuesParentSpecies( c == keep );
        }
    }
}


FossilizedBirthDeathSpeciationProcess::RangeFlow FossilizedBirthDeathSpeciationProcess::updateRanges( const TopologyNode& node )
{
    if( node.isTip() )
    {
        return RangeFlow{ int(node.getIndex()), 0.0 };
    }

    RangeFlow flow = { -1, 0.0 };

    std::vector<TopologyNode* > children = node.getChildren();

    bool sa = node.isSampledAncestorTipOrParent();

    // one continuing child is budding and none is symmetric speciation
    size_t n_cont = 0;
    bool   sa_tip = false;
    const TopologyNode *past_sa = NULL;

    for (size_t c = 0; c < children.size(); c++)
    {
        if ( children[c]->isSampledAncestorTip() == true )
        {
            sa_tip = true;
            if ( children[c]->continuesParentSpecies() == false ) invalid_continuation = true;
        }
        else if ( continuesRange( *children[c] ) == true )
        {
            ++n_cont;
            past_sa = children[c];
        }
    }

    if ( n_cont > 1 ) invalid_continuation = true;

    if ( sa_tip == true && n_cont == 1 )
    {
        // a range carries on past its sampled ancestor only into the symmetric speciation that ends it, at its age
        if ( hasSymmetricSpeciation() == false || past_sa->isTip() == true || past_sa->getAge() != node.getAge() ) invalid_continuation = true;
        for (size_t c = 0; c < past_sa->getNumberOfChildren(); c++)
        {
            if ( past_sa->getChild(c).continuesParentSpecies() == true ) invalid_continuation = true;
        }
    }

    // no child carries the species, so it ends here by symmetric speciation
    if ( sa_tip == false && n_cont == 0 && hasSymmetricSpeciation() == false ) invalid_continuation = true;

    // stop before indexing with species == -1
    if ( invalid_continuation == true ) return flow;

    if ( sa_tip == false )
    {
        ( n_cont == 1 ? budding_ages : symmetric_ages ).push_back( node.getAge() );
    }

    // the range sampled at this node, which a species beginning here continues
    int sa_species = -1;
    for (size_t c = 0; c < children.size(); c++)
    {
        if ( children[c]->isSampledAncestorTip() == true ) sa_species = int( children[c]->getIndex() );
    }

    // a range that ends below without a taxon reaches this node unnamed
    bool pending = false;
    double pending_end  = 0.0;

    for(int c = 0; c < children.size(); c++)
    {
        const TopologyNode& child = *children[c];

        RangeFlow sub = updateRanges(child);

        // propagate a failed subtree
        if ( invalid_continuation == true ) return flow;

        if ( sub.species < 0 )
        {
            if ( continuesRange( child ) == true )
            {
                pending = true;
                pending_end = sub.end_age;
            }
            else
            {
                // a range no taxon samples begins here
                RangeStart start = ( sa == false ? RangeStart::SPECIATION : RangeStart::SAMPLED_ANCESTOR );
                unsampled.push_back( UnsampledRange{ node.getAge(), sub.end_age, start, sa && child.continuesParentSpecies() } );
            }
            continue;
        }

        int i = sub.species;

        // is child a new species?
        // set start time at this node
        if( continuesRange( child ) == false )
        {
            double age = node.getAge(); // y_{a(i)}

            if ( age != ranges[i].divergence )
            {
                ranges[i].divergence = age;
                dirty_taxa[i] = true;
            }

            next_species[i].parent_is_sa = sa;
            next_species[i].sa_parent    = sa ? sa_species : -1;
            next_species[i].sa_anagenetic = sa && child.continuesParentSpecies();
        }
        // child is the ancestral species
        else
        {
            // one range reaching two tips
            if ( flow.species >= 0 )
            {
                invalid_continuation = true;
                return flow;
            }

            // propagate species index
            flow.species = i;

            // its range ends here and the lineage carries on below
            if ( sa == true ) next_species[i].ends_at_sa = true;

            // if this is the root
            // set the start time to the origin
            if( node.isRoot() )
            {
                double age = getOriginAge();

                if ( age != ranges[i].divergence )
                {
                    ranges[i].divergence = age;
                    origin = age;
                    dirty_taxa[i] = true;
                }
            }
        }
    }

    if ( pending == true )
    {
        if ( sa_species >= 0 )
        {
            // the sampled ancestor here names the range that a symmetric speciation ended below
            if ( ranges[sa_species].death != pending_end )
            {
                ranges[sa_species].death = pending_end;
                dirty_taxa[sa_species] = true;
            }
            next_species[sa_species].ends_symmetric = true;
        }
        else
        {
            flow.end_age = pending_end;
        }
    }

    // no child carries this node's species, so it ends here
    if ( sa_tip == false && n_cont == 0 )
    {
        flow.end_age = node.getAge();
    }

    return flow;
}

/**
 *
 *
 */
void FossilizedBirthDeathSpeciationProcess::prepareProbComputation( void ) const 
{
    // beta first, since the base class solves the psi mixture with it
    if ( homogeneous_lambda_a != NULL )
    {
        anagenetic = std::vector<double>(num_intervals, homogeneous_lambda_a->getValue() );
    }
    else
    {
        anagenetic = heterogeneous_lambda_a->getValue();
    }
    if ( homogeneous_beta != NULL )
    {
        symmetric = std::vector<double>(num_intervals, homogeneous_beta->getValue() );
    }
    else
    {
        symmetric = heterogeneous_beta->getValue();
    }

    AbstractFossilizedBirthDeathRangeProcess::prepareProbComputation();

    for (size_t i = 0; i < num_intervals; i++)
    {
        if ( i < num_intervals-1 )
        {
            double dt = times[i+1] - times[i];

            q_tilde_i[i] = - anagenetic[i] * dt - symmetric[i] * (birth[i] + death[i] + fossil[i]) * dt + (1.0 - symmetric[i]) * q_tilde_i[i];
        }
    }
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
/** Re-take from each tip the ages it holds. */
void FossilizedBirthDeathSpeciationProcess::repairRanges( void )
{
    // a tip is an extinction, which may sit below its bins
    for (size_t i = 0; i < getValue().getNumberOfNodes(); i++)
    {
        TopologyNode &node = getValue().getNode(i);
        // an extant tip is at the present, outside its fossil ages
        if ( node.isTip() )
        {
            size_t ti = node.getIndex();
            if ( ti >= taxa.size() ) continue;

            node.setTipAgeUnconstrained( true );

            double age = node.getAge();

            if ( ranges[ti].death != age )
            {
                ranges[ti].death = age;
                dirty_taxa[ti] = true;
            }
        }
    }
}


void FossilizedBirthDeathSpeciationProcess::updateRanges( void )
{
    // a stale flag would freeze the chain
    invalid_continuation = false;

    repairRanges();

    // a fresh table: the root lineage never reaches the assignment
    next_species.assign( taxa.size(), SpeciesEntry() );
    budding_ages.clear();
    symmetric_ages.clear();
    unsampled.clear();

    const TopologyNode &root = getValue().getRoot();

    RangeFlow top = updateRanges(root);

    // the range that begins at the origin may be one no taxon samples
    origin_unsampled = ( invalid_continuation == false && root.isTip() == false && top.species < 0 );
    if ( origin_unsampled == true ) unsampled.push_back( UnsampledRange{ getOriginAge(), top.end_age, RangeStart::SPECIATION, false } );

    // an entry that moved is dirty
    for (size_t i = 0; i < taxa.size(); ++i)
    {
        if ( next_species[i] == species[i] ) continue;

        species[i] = next_species[i];
        dirty_taxa[i] = true;
    }

    // a lone lineage is both root and tip, so the recursion sets neither of its times
    if ( root.isTip() )
    {
        size_t i = root.getIndex();

        if ( getOriginAge() != ranges[i].divergence )
        {
            ranges[i].divergence = getOriginAge();
            dirty_taxa[i] = true;
        }
    }

    max_birth = 0;
    for (size_t i = 0; i < taxa.size(); i++)
    {
        if ( ranges[i].divergence > ranges[max_birth].divergence ) max_birth = i;
    }

    origin = getOriginAge();
}


/** The tree's keep would otherwise reach the tree base only, and the range base must clear its stored proposal. */
void FossilizedBirthDeathSpeciationProcess::keepSpecialization(const DagNode *toucher)
{
    AbstractFossilizedBirthDeathRangeProcess::keepSpecialization(toucher);
}


void FossilizedBirthDeathSpeciationProcess::restoreSpecialization(const DagNode *toucher)
{
    species = stored_species;
    cached_events_level = 0;

    AbstractFossilizedBirthDeathRangeProcess::restoreSpecialization(toucher);
}


void FossilizedBirthDeathSpeciationProcess::touchSpecialization(const DagNode *toucher, bool touchAll)
{
    bool first_touch = ( touched == false );

    if ( toucher == dag_node )
    {
        if ( first_touch )
        {
            storeRangeTerms();
            stored_ranges     = ranges;
        }

        touched = true;
    }
    else
    {
        AbstractBirthDeathProcess::touchSpecialization(toucher, touchAll);
        AbstractFossilizedBirthDeathRangeProcess::touchSpecialization(toucher, touchAll);
    }

    if ( first_touch == true ) stored_species = species;

    // a parameter or the tree moved, so the per-branch events are stale
    cached_events_level = 0;

    prepareProbComputation();
    updateRanges();
}


/**
 * Swap the parameters held by this distribution.
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void FossilizedBirthDeathSpeciationProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == heterogeneous_lambda_a)
    {
        heterogeneous_lambda_a = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == heterogeneous_beta)
    {
        heterogeneous_beta = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == homogeneous_lambda_a)
    {
        homogeneous_lambda_a = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == homogeneous_beta)
    {
        homogeneous_beta = static_cast<const TypedDagNode<double>* >( newP );
    }
    else
    {
        AbstractBirthDeathProcess::swapParameterInternal(oldP, newP);
        AbstractFossilizedBirthDeathRangeProcess::swapParameterInternal(oldP, newP);
    }
}
