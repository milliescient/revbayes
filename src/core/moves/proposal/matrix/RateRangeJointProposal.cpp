#include "RateRangeJointProposal.h"

#include <cmath>
#include <iostream>
#include <limits>

#include "AbstractFossilizedBirthDeathRangeProcess.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TreeDistributionProperties.h"
#include "DistributionNormal.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RbMathLogic.h"

using namespace RevBayesCore;


RateRangeJointProposal::RateRangeJointProposal( StochasticNode<double> *r, StochasticNode<MatrixReal> *m, StochasticNode<Tree> *tr, double l, bool rj, double sd ) : Proposal(),
    rate( r ),
    shift( NULL ),
    ranges( m ),
    tree( tr ),
    lambda( l ),
    rj_switch( rj ),
    sigma( sd ),
    stored_rate( 0.0 ),
    stored_shift( 0 )
{
    // the ranges are touched first, as for a plain matrix or tree move
    addValueNode();
    addNode( rate );
}


RateRangeJointProposal::RateRangeJointProposal( StochasticNode<std::int64_t> *s, StochasticNode<MatrixReal> *m, StochasticNode<Tree> *tr ) : Proposal(),
    rate( NULL ),
    shift( s ),
    ranges( m ),
    tree( tr ),
    lambda( 0.0 ),
    rj_switch( false ),
    sigma( 0.0 ),
    stored_rate( 0.0 ),
    stored_shift( 0 )
{
    addValueNode();
    addNode( shift );
}


void RateRangeJointProposal::addValueNode( void )
{
    if ( tree != NULL )
    {
        // only an extended tree has extinction times as tip ages
        const TreeDistributionProperties *props = dynamic_cast<const TreeDistributionProperties *>( &tree->getDistribution() );
        if ( props == NULL || props->isExtended() == false )
        {
            throw RbException("mvRateRangeJoint requires an extended tree, whose tips are extinction times.");
        }
        addNode( tree );
    }
    else
    {
        addNode( ranges );
    }
}


AbstractFossilizedBirthDeathRangeProcess* RateRangeJointProposal::process( void ) const
{
    Distribution &d = ( tree != NULL ? unconstrained( tree->getDistribution() ) : unconstrained( ranges->getDistribution() ) );
    AbstractFossilizedBirthDeathRangeProcess* p = dynamic_cast<AbstractFossilizedBirthDeathRangeProcess*>( &d );
    if ( p == NULL )
    {
        throw RbException("mvRateRangeJoint: the ranges must be drawn from a fossilized birth-death range process.");
    }
    return p;
}


double RateRangeJointProposal::extinctionTime( size_t i ) const
{
    return ( tree != NULL ) ? tree->getValue().getNode( i ).getAge() : ranges->getValue()[i][0];
}


void RateRangeJointProposal::setExtinctionTime( size_t i, double d )
{
    if ( tree != NULL )
    {
        // on an extended tree the tip age is the extinction time
        tree->getValue().getNode( i ).setAge( d );
    }
    else
    {
        ranges->getValue()[i][0] = d;
        ranges->addTouchedElementIndex( i * 4 );
    }
}


void RateRangeJointProposal::cleanProposal( void )
{
    stored_ends.clear();
}


RateRangeJointProposal* RateRangeJointProposal::clone( void ) const
{
    return new RateRangeJointProposal( *this );
}


const std::string& RateRangeJointProposal::getProposalName( void ) const
{
    static std::string name = "RateRangeJoint";

    return name;
}


double RateRangeJointProposal::getProposalTuningParameter( void ) const
{
    return ( rate != NULL && rj_switch == false ) ? lambda : RbConstants::Double::nan;
}


void RateRangeJointProposal::readRates( AbstractFossilizedBirthDeathRangeProcess *process, std::vector<double> &t, std::vector<double> &mu, std::vector<double> &psi ) const
{
    process->currentExtinctionAndSampling( t, mu, psi );
}


void RateRangeJointProposal::probeRate( double value, AbstractFossilizedBirthDeathRangeProcess *process, std::vector<double> &t, std::vector<double> &mu, std::vector<double> &psi )
{
    // setValue, not the reference: a dnReversibleJumpMixture keeps which branch it is on beside the value
    rate->setValue( new double( value ) );
    readRates( process, t, mu, psi );
    rate->restore();
}


/**
 * The density mu(x) exp(-H) on (lo, hi], with H the hazard mu + psi accumulated down from hi, the last
 * fossil. Interval j is [t[j], t[j+1]), the last one open.
 */
double RateRangeJointProposal::lnDensity( double x, double lo, double hi, const std::vector<double> &t, const std::vector<double> &mu, const std::vector<double> &psi ) const
{
    if ( !( x > lo && x <= hi ) ) return RbConstants::Double::neginf;

    double hazard = 0.0;
    double ln_z = RbConstants::Double::neginf;
    double ln_f = RbConstants::Double::neginf;
    size_t n = t.size();

    for ( size_t s = 0; s < n; s++ )
    {
        size_t j = n - 1 - s;
        double a = std::max( t[j], lo );
        double c = std::min( ( j + 1 < n ) ? t[j+1] : RbConstants::Double::inf, hi );
        if ( c <= a ) continue;

        double h = mu[j] + psi[j];
        double L = c - a;
        if ( mu[j] > 0.0 )
        {
            double ln_mass = log( mu[j] ) - hazard + ( h > 0.0 ? log( -expm1( -h * L ) ) - log( h ) : log( L ) );
            ln_z = ( ln_z == RbConstants::Double::neginf ) ? ln_mass : std::max( ln_z, ln_mass ) + log1p( exp( -std::fabs( ln_z - ln_mass ) ) );

            if ( x > a && x <= c )
            {
                ln_f = log( mu[j] ) - hazard - h * ( c - x );
            }
        }
        hazard += h * L;
    }

    return ln_f - ln_z;
}


double RateRangeJointProposal::draw( double lo, double hi, const std::vector<double> &t, const std::vector<double> &mu, const std::vector<double> &psi ) const
{
    RandomNumberGenerator* rng = GLOBAL_RNG;

    std::vector<double> ln_mass, top_of, length, h;
    double hazard = 0.0;
    size_t n = t.size();

    for ( size_t s = 0; s < n; s++ )
    {
        size_t j = n - 1 - s;
        double a = std::max( t[j], lo );
        double c = std::min( ( j + 1 < n ) ? t[j+1] : RbConstants::Double::inf, hi );
        if ( c <= a ) continue;

        double hj = mu[j] + psi[j];
        double L = c - a;
        ln_mass.push_back( mu[j] > 0.0 ? log( mu[j] ) - hazard + ( hj > 0.0 ? log( -expm1( -hj * L ) ) - log( hj ) : log( L ) ) : RbConstants::Double::neginf );
        top_of.push_back( c );
        length.push_back( L );
        h.push_back( hj );

        hazard += hj * L;
    }

    double top = RbConstants::Double::neginf;
    for ( double m : ln_mass ) top = std::max( top, m );
    if ( top == RbConstants::Double::neginf ) return RbConstants::Double::nan;

    double total = 0.0;
    for ( double m : ln_mass ) total += exp( m - top );

    double u = rng->uniform01() * total;
    size_t k = 0;
    for ( ; k + 1 < ln_mass.size(); k++ )
    {
        u -= exp( ln_mass[k] - top );
        if ( u < 0.0 ) break;
    }

    // the cell's exponential in the distance down from its top, truncated to the cell
    double v = rng->uniform01();
    double dist = ( h[k] > 0.0 ) ? -log1p( v * expm1( -h[k] * length[k] ) ) / h[k] : v * length[k];
    return top_of[k] - dist;
}


double RateRangeJointProposal::lnLognormal( double x, double centre ) const
{
    double z = ( log( x ) - log( centre ) ) / sigma;
    return -log( x ) - log( sigma ) - 0.5 * log( 2.0 * M_PI ) - 0.5 * z * z;
}


double RateRangeJointProposal::doProposal( void )
{
    RandomNumberGenerator* rng = GLOBAL_RNG;

    AbstractFossilizedBirthDeathRangeProcess* process = this->process();

    std::vector<double> t, mu_old, psi_old;
    readRates( process, t, mu_old, psi_old );
    double present = t.front();

    // every extinction time free to move, which is one not on the present, below its last fossil. On a
    // tree it is also below its parent, and a sampled ancestor tip or a tip at its parent's age is a point
    // mass that a continuous draw could not return to.
    std::vector<size_t> candidates;
    std::vector<double> lo, hi;
    for ( size_t i = 0; i < process->getTaxa().size(); i++ )
    {
        double d = extinctionTime( i );
        if ( process->hasSurvived(i) == true || d <= present ) continue;

        double top = process->getLastAppearance( i );
        if ( tree != NULL )
        {
            const TopologyNode &tip = tree->getValue().getNode( i );
            if ( tip.isRoot() == true || tip.isSampledAncestorTip() == true || d >= tip.getParent().getAge() ) continue;
            top = std::min( top, tip.getParent().getAge() );
        }
        candidates.push_back( i );
        lo.push_back( present );
        hi.push_back( top );
    }

    // set the new rate and touch it only to read the rates through the DAG, then restore, so the
    // move's own touches start from a clean DAG with the matrix first, as for a plain matrix move
    double ln_hastings = 0.0;
    std::vector<double> t_new, mu_new, psi_new;
    if ( rate != NULL && rj_switch == true )
    {
        stored_rate = rate->getValue();
        if ( stored_rate == 0.0 )
        {
            // a probe value no cell holds marks the cells this jump governs, which share one rate now
            double probe = 1.0;
            for ( double m : mu_old ) probe = std::max( probe, 2.0 * m + 1.0 );
            std::vector<double> mu_probe;
            probeRate( probe, process, t_new, mu_probe, psi_new );

            double centre = RbConstants::Double::nan;
            for ( size_t j = 0; j < mu_probe.size(); j++ )
            {
                if ( mu_probe[j] != mu_old[j] ) { centre = mu_old[j]; break; }
            }
            if ( !( centre > 0.0 ) )
            {
                rate->setValue( new double( stored_rate ) );
                return RbConstants::Double::neginf;
            }

            double x = centre * exp( sigma * RbStatistics::Normal::rv( *rng ) );
            probeRate( x, process, t_new, mu_new, psi_new );
            ln_hastings = -lnLognormal( x, centre );
        }
        else
        {
            probeRate( 0.0, process, t_new, mu_new, psi_new );

            // the rate the governed cells fall back to centres the switch-on this switch-off reverses
            double centre = RbConstants::Double::nan;
            for ( size_t j = 0; j < mu_new.size(); j++ )
            {
                if ( mu_new[j] != mu_old[j] ) { centre = mu_new[j]; break; }
            }
            if ( !( centre > 0.0 ) ) return RbConstants::Double::neginf;
            ln_hastings = lnLognormal( stored_rate, centre );
        }
    }
    else if ( rate != NULL )
    {
        stored_rate = rate->getValue();
        double m = lambda * ( rng->uniform01() - 0.5 );
        probeRate( stored_rate * exp( m ), process, t_new, mu_new, psi_new );
        ln_hastings = m;
    }
    else
    {
        stored_shift = shift->getValue();
        shift->setValue( new std::int64_t( 1 - stored_shift ) );
        readRates( process, t_new, mu_new, psi_new );
        shift->restore();
    }

    std::vector<bool> changed( mu_new.size(), false );
    bool any = false;
    for ( size_t j = 0; j < mu_new.size(); j++ )
    {
        changed[j] = ( mu_new[j] != mu_old[j] );
        any = any || changed[j];
    }

    stored_ends.clear();
    if ( any == true )
    {
        for ( size_t k = 0; k < candidates.size(); k++ )
        {
            size_t i = candidates[k];

            // only an end whose window meets an interval where the rate changed
            bool hit = false;
            for ( size_t j = 0; j < changed.size() && hit == false; j++ )
            {
                double top = ( j + 1 < t.size() ) ? t[j+1] : RbConstants::Double::inf;
                hit = changed[j] && t[j] < hi[k] && top > lo[k];
            }
            if ( hit == false ) continue;

            double x_old = extinctionTime( i );
            double x_new = draw( lo[k], hi[k], t_new, mu_new, psi_new );

            stored_ends.push_back( std::make_pair( i, x_old ) );
            if ( RbMath::isFinite( x_new ) == false )
            {
                return RbConstants::Double::neginf;
            }
            setExtinctionTime( i, x_new );

            ln_hastings += lnDensity( x_old, lo[k], hi[k], t, mu_old, psi_old ) - lnDensity( x_new, lo[k], hi[k], t_new, mu_new, psi_new );
        }
    }

    return ln_hastings;
}


void RateRangeJointProposal::prepareProposal( void )
{

}


void RateRangeJointProposal::printParameterSummary( std::ostream &o, bool name_only ) const
{
    if ( rate != NULL && rj_switch == true )
    {
        o << "sigma = ";
        if ( name_only == false ) o << sigma;
    }
    else if ( rate != NULL )
    {
        o << "lambda = ";
        if ( name_only == false ) o << lambda;
    }
}


void RateRangeJointProposal::undoProposal( void )
{
    for ( std::vector< std::pair<size_t,double> >::reverse_iterator it = stored_ends.rbegin(); it != stored_ends.rend(); ++it )
    {
        setExtinctionTime( it->first, it->second );
    }
    stored_ends.clear();
    if ( ranges != NULL ) ranges->clearTouchedElementIndices();

    // through setValue, as ScaleProposalContinuous does, not by writing back through the reference
    if ( rate != NULL )
    {
        rate->setValue( new double( stored_rate ) );
    }
    else
    {
        shift->setValue( new std::int64_t( stored_shift ) );
    }
}


void RateRangeJointProposal::swapNodeInternal( DagNode *oldN, DagNode *newN )
{
    if ( ranges != NULL && oldN == ranges )
    {
        ranges = static_cast< StochasticNode<MatrixReal>* >( newN );
    }
    else if ( tree != NULL && oldN == tree )
    {
        tree = static_cast< StochasticNode<Tree>* >( newN );
    }
    else if ( oldN == rate )
    {
        rate = static_cast< StochasticNode<double>* >( newN );
    }
    else if ( oldN == shift )
    {
        shift = static_cast< StochasticNode<std::int64_t>* >( newN );
    }
}


void RateRangeJointProposal::setProposalTuningParameter( double tp )
{
    if ( rate != NULL && rj_switch == false ) lambda = tp;
}


void RateRangeJointProposal::tune( double r )
{
    if ( rate == NULL || rj_switch == true ) return;

    double p = this->targetAcceptanceRate;
    if ( r > p )
    {
        lambda *= ( 1.0 + ( (r - p) / (1.0 - p) ) );
    }
    else
    {
        lambda /= ( 2.0 - r / p );
    }
    lambda = fmin( 10000, lambda );
}
