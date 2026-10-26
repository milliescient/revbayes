#include "AbstractFossilizedBirthDeathRangeProcess.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ostream>
#include <sstream>
#include <vector>


#include "DagNode.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbBoolean.h"
#include "RbMathCombinatorialFunctions.h"
#include "RbMathLogic.h"
#include "RbException.h"
#include "RbMathFunctions.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "RlUserInterface.h"
#include "Taxon.h"
#include "TimeInterval.h"
#include "TypedDagNode.h"

using namespace RevBayesCore;

/**
 * Constructor.
 *
 * \param[in]    s              Speciation rates.
 * \param[in]    e              Extinction rates.
 * \param[in]    p              Fossil sampling rates.
 * \param[in]    r              Instantaneous sampling probabilities.
 * \param[in]    t              Rate change times.
 * \param[in]    cdt            Condition of the process (time/sampling/survival).
 * \param[in]    tn             Taxa.
 * \param[in]    c              Complete sampling?
 */
AbstractFossilizedBirthDeathRangeProcess::AbstractFossilizedBirthDeathRangeProcess(const DagNode *inspeciation,
                                                                         const DagNode *inextinction,
                                                                         const DagNode *inpsi,
                                                                         const TypedDagNode<double> *inrho,
                                                                         const TypedDagNode< RbVector<double> > *intimes,
                                                                         const std::string &incondition,
                                                                         const std::vector<Taxon> &intaxa,
                                                                         bool comp,
                                                                         const TypedDagNode<double> *inorigin,
                                                                         double inpresent,
                                                                         const TypedDagNode<RbVector<double> > *inrates,
                                                                         const TypedDagNode<Simplex> *inweights) :
    taxa(intaxa),
    condition(incondition),
    homogeneous_rho(inrho),
    timeline( intimes ),
    origin_age( inorigin ),
    psi_mixture_rates( inrates ),
    psi_mixture_weights( inweights ),
    num_components( inrates == NULL ? 1 : inrates->getValue().size() ),
    origin(0.0),
    max_birth(0),
    resampled(false),
    touched(false),
    present_time(inpresent)
{
    // initialize all the pointers to NULL
    homogeneous_lambda             = NULL;
    homogeneous_mu                 = NULL;
    homogeneous_psi                = NULL;
    heterogeneous_lambda           = NULL;
    heterogeneous_mu               = NULL;
    heterogeneous_psi              = NULL;

    // cast the pointers from their input parameters
    heterogeneous_lambda    = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inspeciation);
    homogeneous_lambda      = dynamic_cast<const TypedDagNode<double >*>(inspeciation);
    heterogeneous_mu        = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inextinction);
    homogeneous_mu          = dynamic_cast<const TypedDagNode<double >*>(inextinction);
    heterogeneous_psi       = dynamic_cast<const TypedDagNode<RbVector<double> >*>(inpsi);
    homogeneous_psi         = dynamic_cast<const TypedDagNode<double >*>(inpsi);

    // add the parameters to the model
    range_parameters.push_back( timeline );
    range_parameters.push_back( origin_age );
    range_parameters.push_back( homogeneous_rho );
    range_parameters.push_back( homogeneous_lambda );
    range_parameters.push_back( heterogeneous_lambda );
    range_parameters.push_back( homogeneous_mu );
    range_parameters.push_back( heterogeneous_mu );
    range_parameters.push_back( homogeneous_psi );
    range_parameters.push_back( heterogeneous_psi );
    range_parameters.push_back( psi_mixture_rates );
    range_parameters.push_back( psi_mixture_weights );

    // setup the timeline
    if ( timeline == NULL )
    {
        num_intervals = 1;
    }
    else
    {
        num_intervals = timeline->getValue().size() + (timeline->getValue().front() != present_time);
    }

    if ( num_intervals > 1 )
    {
        std::vector<double> times = timeline->getValue();
        std::vector<double> times_sorted_ascending = times;

        sort(times_sorted_ascending.begin(), times_sorted_ascending.end() );

        if ( times != times_sorted_ascending )
        {
            throw(RbException("Interval times must be provided in ascending order"));
        }

        if ( times.front() < present_time )
        {
            throw(RbException("An interval time is younger than present"));
        }
    }

    RbException no_timeline_err = RbException("No time intervals provided for heterogeneous fossilized birth death process");

    RbException inconsistent_rates = RbException("Inconsistent number of rates in fossilized birth death process.");

    size_t num_rates = 0;

    if( heterogeneous_lambda != NULL || heterogeneous_mu != NULL || heterogeneous_psi != NULL)
    {
        if ( timeline == NULL ) throw(no_timeline_err);

        if( heterogeneous_lambda != NULL )
        {
            num_rates = heterogeneous_lambda->getValue().size();
        }
        if( heterogeneous_mu != NULL )
        {
            if ( num_rates == 0 ) num_rates = heterogeneous_mu->getValue().size();

            if ( heterogeneous_mu->getValue().size() != num_rates ) throw(inconsistent_rates);
        }
        if( heterogeneous_psi != NULL )
        {
            if ( num_rates == 0 ) num_rates = heterogeneous_psi->getValue().size();

            if ( heterogeneous_psi->getValue().size() != num_rates ) throw(inconsistent_rates);
        }
    }
    else
    {
        num_rates = 1;
    }

    if ( num_rates != num_intervals )
    {
        std::stringstream ss;
        ss << "Number of rates does not match number of time intervals in fossilized birth death process.";
        throw(RbException(ss.str()));
    }

    if ( num_components == 0 )
    {
        throw RbException("psiMixtureRates must have at least one component.");
    }
    if ( psi_mixture_weights != NULL && psi_mixture_weights->getValue().size() != num_components )
    {
        throw RbException("psiMixtureWeights must have one weight per psiMixtureRates component.");
    }

    ranges = std::vector<RangeEntry>(taxa.size());

    // start every species in the component nearest the mean rate
    size_t start = 0;
    if ( num_components > 1 )
    {
        const std::vector<double> &rates = psi_mixture_rates->getValue();
        for (size_t h = 1; h < num_components; h++)
        {
            if ( fabs( log( rates[h] ) ) < fabs( log( rates[start] ) ) ) start = h;
        }
    }
    for (size_t i = 0; i < taxa.size(); i++)
    {
        ranges[i].component = start;
    }

    mixture_p_i       = std::vector<std::vector<double> >(num_intervals);
    lambda_pbar_i     = std::vector<double>(num_intervals, 0.0);
    mixture_q_i       = std::vector<std::vector<double> >(num_intervals);
    mixture_q_scale_i = std::vector<std::vector<double> >(num_intervals);

    p_i         = std::vector<double>(num_intervals, 1.0);
    pS_i        = std::vector<double>(num_intervals, 1.0);
    q_i         = std::vector<double>(num_intervals, 0.0);
    q_tilde_i   = std::vector<double>(num_intervals, 0.0);

    birth       = std::vector<double>(num_intervals, 0.0);
    death       = std::vector<double>(num_intervals, 0.0);
    fossil      = std::vector<double>(num_intervals, 0.0);
    times       = std::vector<double>(num_intervals, 0.0);

    partial_likelihood = std::vector<double>(taxa.size(), 0.0);


    dirty_taxa  = std::vector<bool>(taxa.size(), true);

    record_reporting = comp ? COMPLETE : FIRSTLAST;

    updateRecord();

    stored_ranges = ranges;

    prepareProbComputation();

    if ( times.front() > max_present_age )
    {
        throw(RbException("Timeline start time is older than youngest fossil occurrence."));
    }
}


/** Derive from the occurrence record the reporting model and the bounds on the appearances. */
void AbstractFossilizedBirthDeathRangeProcess::updateRecord( void )
{
    // the record fixes the bounds and the reporting model; the sampled ages survive a re-derive
    for ( size_t i = 0; i < taxa.size(); i++ )
    {
        ranges[i].first_min = 0.0;
        ranges[i].first_max = RbConstants::Double::inf;
        ranges[i].last_min  = 0.0;
        ranges[i].last_max  = RbConstants::Double::inf;
        ranges[i].singleton = true;
    }

    max_present_age = RbConstants::Double::inf;

    for ( size_t i = 0; i < taxa.size(); i++ )
    {
        const std::map<TimeInterval, size_t> &record = taxa[i].getOccurrences();

        size_t count = 0;
        for ( std::map<TimeInterval, size_t>::const_iterator Fi = record.begin(); Fi != record.end(); Fi++ )
        {
            const TimeInterval &bin = Fi->first;

            count += Fi->second;

            // find the oldest minimum age
            ranges[i].first_min = std::max(bin.getMin(), ranges[i].first_min);
            // find the youngest maximum age
            ranges[i].last_max = std::min(bin.getMax(), ranges[i].last_max);

            max_present_age = std::min(max_present_age, ranges[i].last_max);
        }

        ranges[i].singleton = ( count < 2 );

        ranges[i].first_max = taxa[i].getMaxAge();
        ranges[i].last_min  = taxa[i].getMinAge();

        // only used when the record has two extremes to order
        ranges[i].last = ranges[i].last_max;
    }

    checkPresenceIntervals();
}


/** The reporting model, which a presence record has to meet. */
void AbstractFossilizedBirthDeathRangeProcess::setReporting( Reporting r )
{
    record_reporting = r;
    checkPresenceIntervals();
}


/** Presence counts the occurrences inside each interval, so a taxon's intervals may not share any. */
void AbstractFossilizedBirthDeathRangeProcess::checkPresenceIntervals( void ) const
{
    if ( record_reporting != PRESENCE ) return;

    for ( size_t i = 0; i < taxa.size(); i++ )
    {
        const auto &ages = taxa[i].getOccurrences();

        for ( auto a = ages.begin(); a != ages.end(); a++ )
        {
            for ( auto b = std::next(a); b != ages.end(); b++ )
            {
                // a shared edge is adjacency, but an exact age inside an interval is in it
                double lo = std::max( a->first.getMin(), b->first.getMin() );
                double hi = std::min( a->first.getMax(), b->first.getMax() );
                bool inside = ( a->first.getMin() < lo && lo < a->first.getMax() ) || ( b->first.getMin() < lo && lo < b->first.getMax() );

                if ( lo < hi || ( lo == hi && inside ) )
                {
                    std::ostringstream ss;
                    ss << "reporting=\"presence\" needs disjoint intervals, and " << taxa[i].getName()
                       << " reports [" << a->first.getMin() << "," << a->first.getMax() << "] and ["
                       << b->first.getMin() << "," << b->first.getMax() << "]";
                    throw RbException( ss.str() );
                }
            }
        }
    }
}


/** Adopt a clamped record's occurrences as the data and re-derive the record model. */
void AbstractFossilizedBirthDeathRangeProcess::setOccurrences( const std::vector<Taxon> &t )
{
    // the ordinary case, which must leave every sampled age exactly as it is
    bool changed = ( t.size() != taxa.size() );

    for (size_t i = 0; i < taxa.size() && changed == false; ++i)
    {
        if ( taxa[i].getOccurrences() != t[i].getOccurrences() ) changed = true;
        if ( taxa[i].isExtinct()     != t[i].isExtinct()     ) changed = true;
    }

    if ( changed == false ) return;

    taxa = t;

    updateRecord();

    prepareProbComputation();

    if ( times.front() > max_present_age )
    {
        throw(RbException("dnFossilRecord: the clamped record has an occurrence younger than the timeline start."));
    }

    dirty_taxa = std::vector<bool>(taxa.size(), true);

}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double AbstractFossilizedBirthDeathRangeProcess::computeLnProbabilityRanges( bool force )
{
    // prepare the probability computation
    prepareProbComputation();

    // only the oldest range begins at the origin
    if ( origin_age != NULL )
    {
        if ( ranges[max_birth].divergence != origin ) return RbConstants::Double::neginf;

        for (size_t i = 0; i < taxa.size(); ++i)
        {
            if ( i != max_birth && ranges[i].divergence >= origin ) return RbConstants::Double::neginf;
        }
    }

    // variable declarations and initialization
    double lnProb = 0.0;

    // add the range terms
    for (size_t i = 0; i < taxa.size(); ++i)
    {
        double b = ranges[i].divergence;
        double d = ranges[i].death;

        double o = ranges[i].first;

        double max_age = taxa[i].getMaxAge();
        double min_age = taxa[i].getMinAge();

        double present = times.front();

        // a move on the value can put an age anywhere, so enforce the support for every taxon
        if ( ranges[i].isOrdered( present ) == false )
        {
            return RbConstants::Double::neginf;
        }
        // a survivor ends at the present, so a continuous move that walks d off it is rejected here
        if ( ranges[i].survived == true && d != present )
        {
            return RbConstants::Double::neginf;
        }
        if ( dirty_taxa[i] == true || force )
        {
            partial_likelihood[i] = rangeLnProb( i );
        }

        lnProb += partial_likelihood[i];
    }

    lnProb += originLnProb();

    lnProb += conditionLnProb();

    if ( RbMath::isFinite(lnProb) == false )
    {
        return RbConstants::Double::neginf;
    }

    return lnProb;
}


/** Birth-death density of range i. */
double AbstractFossilizedBirthDeathRangeProcess::rangeLnProb( size_t i )
{
    if ( num_components > 1 )
    {
        std::vector<double> ln_probs;
        mixtureRangeLnProbs( i, ln_probs, ranges[i].component );

        return ln_probs[ ranges[i].component ];
    }

    double b = ranges[i].divergence;
    double d = ranges[i].death;
    double o = ranges[i].first;

    double present = times.front();

    size_t bi = findIndex(b);
    size_t oi = findIndex(o);
    size_t di = findIndex(d);

    // include speciation density
    double lnProb = log( birth[bi] );

    // multiply by q at the birth time
    lnProb += q(bi, b);

    // include intermediate q terms
    for (size_t j = oi; j < bi; j++)
    {
        lnProb += q_i[j];
    }

    // skip the rest for extant taxa with no fossil samples
    if ( taxa[i].getMaxAge() == present )
    {
        return lnProb;
    }

    // replace q terms at oldest occurrence
    lnProb += q(oi, o, true) - q(oi, o);

    // include intermediate q_tilde terms
    for (size_t j = di; j < oi; j++)
    {
        lnProb += q_tilde_i[j];
    }

    // divide by q_tilde at the death time
    lnProb -= q( di, d, true);

    // the extinction density, or what closes a marginalized range
    if ( ranges[i].survived == false )
    {
        lnProb += rangeEndTerm( i, di, d );
    }

    return lnProb;
}


/** The oldest birth starts the record, so its speciation rate divides out. */
double AbstractFossilizedBirthDeathRangeProcess::originLnProb( void )
{
    return -log( birth[ findIndex(origin) ] );
}


// The count at the present: geometric with an atom at alpha.
void AbstractFossilizedBirthDeathRangeProcess::numTaxaGeometric( double &alpha, double &beta ) const
{
    // (A z + B) / (C z + D), composed from the present up to the origin
    double A = 1.0, B = 0.0, C = 0.0, D = 1.0;

    size_t ori = findIndex(origin);

    for (size_t i = 0; i <= ori; i++)
    {
        double b  = birth[i];
        double d  = death[i];
        double dt = ( i < ori ? times[i+1] : origin ) - times[i];

        if ( dt <= 0.0 )
        {
            continue;
        }

        double diff = b - d;
        double a_i, b_i;

        // b == d leaves the closed form singular
        if ( fabs(diff) * dt < 1E-10 )
        {
            a_i = d * dt / (1.0 + d * dt);
            b_i = b * dt / (1.0 + b * dt);
        }
        else
        {
            double e   = -expm1(-diff * dt);
            double den = diff - d * expm1(-diff * dt);
            a_i = d * e / den;
            b_i = b * e / den;
        }

        // this interval's map, a + (1-a)(1-b) z / (1 - b z)
        double fA = (1.0 - a_i) * (1.0 - b_i) - a_i * b_i;
        double fB = a_i;
        double fC = -b_i;

        double nA = fA * A + fB * C;
        double nB = fA * B + fB * D;
        double nC = fC * A + C;
        double nD = fC * B + D;

        // defined up to scale, so normalise to keep the entries bounded
        double s = std::max( std::max(fabs(nA), fabs(nB)), std::max(fabs(nC), fabs(nD)) );
        if ( s > 0.0 )
        {
            nA /= s; nB /= s; nC /= s; nD /= s;
        }

        A = nA; B = nB; C = nC; D = nD;
    }

    // rho thins each survivor: z -> (1-rho) + rho z
    double r = homogeneous_rho->getValue();

    double tB = A * (1.0 - r) + B;
    double tC = C * r;
    double tD = C * (1.0 - r) + D;

    alpha = tB / tD;
    beta  = -tC / tD;
}


// The taxa the record reports as extant.
size_t AbstractFossilizedBirthDeathRangeProcess::getNumberOfTaxaAtPresent( void ) const
{
    size_t n = 0;
    for (size_t i = 0; i < taxa.size(); i++)
    {
        n += ( taxa[i].isExtinct() == false );
    }

    return n;
}


// log P(n taxa at the present).
double AbstractFossilizedBirthDeathRangeProcess::lnProbNumTaxa( size_t n ) const
{
    double alpha = 0.0, beta = 0.0;
    numTaxaGeometric(alpha, beta);

    if ( n == 0 )
    {
        return ( alpha > 0.0 ? log(alpha) : RbConstants::Double::neginf );
    }

    if ( alpha >= 1.0 || beta >= 1.0 || beta < 0.0 )
    {
        return RbConstants::Double::neginf;
    }

    double lnProb = log(1.0 - alpha) + log(1.0 - beta);

    // beta == 0 forbids more than one survivor, so 0 * log(0) is never evaluated
    if ( n > 1 )
    {
        lnProb += (n - 1) * log(beta);
    }

    return lnProb;
}


/** Conditioning that applies to the process as a whole rather than range by range. */
double AbstractFossilizedBirthDeathRangeProcess::conditionLnProb( void ) const
{
    size_t ori = findIndex(origin);

    if ( condition == "sampling" && num_components > 1 )
    {
        std::vector<double> p;
        double lambda_pbar = 0.0;
        mixtureAt( origin, p, lambda_pbar );

        double pbar = 0.0;
        for (size_t h = 0; h < num_components; h++) pbar += mixture_weights[h] * p[h];

        return -log( 1.0 - pbar );
    }
    if ( condition == "sampling" )
    {
        return -log( 1.0 - p(ori, origin, false) );
    }
    if ( condition == "survival" )
    {
        return -log( 1.0 - p(ori, origin, true) );
    }
    if ( condition == "nTaxa" )
    {
        return -lnProbNumTaxa( getNumberOfTaxaAtPresent() );
    }

    return 0.0;
}



// Total fossil-occurrence log-density, for a dnFossilRecord node.
double AbstractFossilizedBirthDeathRangeProcess::computeLnFossilTotal()
{
    prepareProbComputation();

    size_t num_rho_sampled = 0;
    size_t num_rho_unsampled = 0;

    double lnProb = 0.0;
    for ( size_t i = 0; i < taxa.size(); ++i )
    {
        // status is data, so it is scored here and not by the bare range process
        if ( taxa[i].isExtinct() == false && ranges[i].survived == false )
        {
            return RbConstants::Double::neginf;
        }

        num_rho_sampled   += ( taxa[i].isExtinct() == false );          // l
        num_rho_unsampled += ( taxa[i].isExtinct() && ranges[i].survived && marginalizesExtinction() == false );    // n - m - l

        // an extant taxon with no fossil sample carries no reporting term
        if ( taxa[i].getMaxAge() == times.front() ) continue;

        double r = computeLnFossilRecord(i);
        if ( r == RbConstants::Double::neginf )
        {
            return RbConstants::Double::neginf;
        }
        lnProb += r;
    }

    double rho = homogeneous_rho->getValue();

    if ( num_rho_sampled > 0 )                          // seen at the present, so rho > 0
    {
        if ( rho == 0.0 ) return RbConstants::Double::neginf;
        lnProb += num_rho_sampled * log( rho );
    }
    if ( num_rho_unsampled > 0 )                        // survived unseen, so rho < 1
    {
        if ( rho == 1.0 ) return RbConstants::Double::neginf;
        lnProb += num_rho_unsampled * log( 1.0 - rho );
    }

    return lnProb;
}


/** Fossil-occurrence log-term for taxon i. */
double AbstractFossilizedBirthDeathRangeProcess::computeLnFossilRecord( size_t i ) const
{
    return computeLnFossilRecord( i, num_components > 1 ? mixture_rates[ ranges[i].component ] : 1.0 );
}


/** Fossil-occurrence log-term for taxon i, with psi scaled by a mixture rate. */
double AbstractFossilizedBirthDeathRangeProcess::computeLnFossilRecord( size_t i, double mixture_rate ) const
{
    // a component that never samples cannot carry a fossil, and the terms below would give NaN
    if ( mixture_rate == 0.0 ) return RbConstants::Double::neginf;

    if ( record_reporting == PRESENCE ) return lnPresenceRecord( i, mixture_rate );

    double d = ranges[i].death;
    double o = ranges[i].first;                             // oldest appearance
    size_t oi = findIndex(o);
    double y  = ranges[i].last;                             // youngest appearance
    size_t yi = findIndex(y);

    // the interior runs down to tau_K, or to the range end where the youngest is unreported
    double min_age = ranges[i].singleton == false ? y : d;
    double max_age = o;

    double result = 0.0;

    const std::map<TimeInterval, size_t> &ages = taxa[i].getOccurrences();

    // if there is a range of fossil ages
    if ( min_age != max_age )
    {
        double Lambda = 0.0;                            // interior sampling rate over (last, first)

        std::vector<double> psi(ages.size(), 0.0);

        for (size_t j = 0; j < num_intervals; j++)
        {
            double t_0 = ( j < num_intervals-1 ? times[j+1] : RbConstants::Double::inf );

            if ( t_0 <= min_age )
            {
                continue;
            }
            if ( times[j] >= max_age )
            {
                break;
            }

            // the kappa interior span
            double int_hi = std::min(max_age, t_0);
            double int_lo = std::max(min_age, times[j]);
            double dti = int_hi - int_lo;
            if ( dti > 0.0 ) Lambda += mixture_rate*fossil[j]*dti;

            size_t k = 0;
            // increase running psi total for each observation
            for ( std::map<TimeInterval, size_t>::const_iterator Fi = ages.begin(); Fi != ages.end(); Fi++,k++ )
            {
                if ( Fi->first.getMin() < t_0 && Fi->first.getMax() > times[j] )
                {
                    double dt = 1.0;

                    // only compute dt if this is a non-singleton
                    if ( Fi->first.getMin() != Fi->first.getMax() )
                    {
                        // observed occurrence sampling rate over its interval
                        dt = std::min(Fi->first.getMax(), int_hi) - std::max(Fi->first.getMin(), int_lo);
                    }

                    psi[k] += mixture_rate*fossil[j] * dt;
                }
            }
        }

         // instantaneous rate of the oldest specimen at tau1
        result = log(mixture_rate*fossil[oi]);

        size_t count = 0;
        for ( std::map<TimeInterval, size_t>::const_iterator Fi = ages.begin(); Fi != ages.end(); Fi++ ) count += Fi->second;
        double recip_old = 0.0;                         // sum_{i: tau1 in F_i} count_i / Psi(F_i)
        double recip_young = 0.0;                       // sum_{i: tauK in F_i} count_i / Psi(F_i)
        double diag = 0.0;

        size_t k = 0;
        // compute factors of the sum over each possible oldest/youngest observation
        for ( std::map<TimeInterval, size_t>::const_iterator Fi = ages.begin(); Fi != ages.end(); Fi++,k++ )
        {
            bool eligible_oldest = ( Fi->first.getMax() >= o );
            bool eligible_youngest = ( Fi->first.getMin() <= y );

            // compute sum of reciprocal oldest ranges
            if ( eligible_oldest )
            {
                recip_old += Fi->second / psi[k];
            }

            // compute sum of reciprocal youngest ranges
            if ( eligible_youngest )
            {
                recip_young += Fi->second / psi[k];
            }

            // intervals that could supply both extremes contribute the diagonal term
            if ( eligible_oldest && eligible_youngest )
            {
                diag += Fi->second / (psi[k]*psi[k]);
            }

            // compute product of ranges
            result += log(psi[k]) * Fi->second;
        }

        // sum over which observation supplies the oldest specimen at tau1
        result += log(recip_old) - RbMath::lnFactorial(int(count));

        if ( ranges[i].singleton == false )
        {
            if ( !( o >= y && y >= d && y <= ranges[i].last_max ) )
            {
                return RbConstants::Double::neginf;
            }
            // excluding the diagonal, where a single occurrence supplies both extremes
            result += log(mixture_rate*fossil[yi]) + log(recip_young - diag/recip_old);

            // the first/last rule leaves the interior count unreported, so it is marginalized
            if ( record_reporting == FIRSTLAST )
            {
                double S1 = 0.0, f = 1.0;
                for ( size_t kap = 0; kap < 200; kap++ )
                {
                    S1 += f;
                    f *= Lambda / double(count - 1 + kap);
                    if ( f < 1e-16 * S1 ) break;
                }
                result += log(S1);
            }
        }
    }
    // only one fossil age
    else
    {
        // include instantaneous sampling density
        result = ages.begin()->second * log(mixture_rate*fossil[oi]);
    }

    return result;
}


/** The sampling rate integrated over [lo, hi), which the rate timeline makes piecewise. */
double AbstractFossilizedBirthDeathRangeProcess::samplingExposure( double lo, double hi ) const
{
    double total = 0.0;

    for ( size_t j = 0; j < num_intervals; j++ )
    {
        double top = ( j < num_intervals-1 ? times[j+1] : RbConstants::Double::inf );
        if ( top <= lo ) continue;
        if ( times[j] >= hi ) break;

        total += fossil[j] * ( std::min(hi, top) - std::max(lo, times[j]) );
    }

    return total;
}


namespace {

    /** log of the sum over k >= m of x^k / k!, which is e^x when m <= 0. */
    double lnPoissonTail( double x, long m )
    {
        if ( m <= 0 ) return x;
        if ( x <= 0.0 ) return RbConstants::Double::neginf;

        // above the count the tail is most of e^x, so subtract the head instead
        if ( x > m )
        {
            double head = 0.0, f = exp(-x);
            for ( long k = 0; k < m; k++ )
            {
                head += f;
                f *= x / double(k + 1);
            }
            return x + log1p( -head );
        }

        double sum = 1.0, ratio = 1.0;
        for ( long k = m + 1; ratio >= 1e-16 * sum; k++ )
        {
            ratio *= x / double(k);
            sum += ratio;
        }

        return m * log(x) - RbMath::lnFactorial( int(m) ) + log(sum);
    }
}


/**
 * Presence log-term for taxon i, relative to a process that already charges every stretch of the
 * range for holding no occurrence. A reported interval holds at least its reported count, so it
 * refunds the Poisson tail over its exposure inside [tau_K, tau_1], and an exactly dated occurrence
 * is a density, as under the complete rule. tau_1 and tau_K are occurrences, each taken from an
 * interval holding it, and a singleton is one occurrence that is both.
 */
double AbstractFossilizedBirthDeathRangeProcess::lnPresenceRecord( size_t i, double mixture_rate ) const
{
    bool single = ranges[i].singleton;
    double o = ranges[i].first;
    double y = single ? o : ranges[i].last;

    if ( !( o >= y && y >= ranges[i].death ) ) return RbConstants::Double::neginf;

    std::vector<double> exposure;
    std::vector<long>   need;
    std::vector<bool>   exact, holds_first, holds_last;

    for ( const auto &F : taxa[i].getOccurrences() )
    {
        double lo = std::max( F.first.getMin(), y );
        double hi = std::min( F.first.getMax(), o );
        bool point = ( F.first.getMin() == F.first.getMax() );

        if ( point && lo > hi ) return RbConstants::Double::neginf;

        exposure.push_back( mixture_rate * ( point ? fossil[findIndex(lo)] : ( hi > lo ? samplingExposure(lo, hi) : 0.0 ) ) );
        need.push_back( long(F.second) );
        exact.push_back( point );
        holds_first.push_back( F.first.getMin() <= o && o <= F.first.getMax() );
        holds_last.push_back( F.first.getMin() <= y && y <= F.first.getMax() );
    }

    size_t n = need.size();
    std::vector<double> terms;

    for ( size_t a = 0; a < n; a++ )
    {
        if ( holds_first[a] == false ) continue;

        for ( size_t b = 0; b < n; b++ )
        {
            if ( single == true ? b != a : holds_last[b] == false ) continue;

            double term = 0.0;
            for ( size_t j = 0; j < n; j++ )
            {
                long m = need[j] - ( j == a ) - ( single == false && j == b );
                if ( exact[j] == false ) term += lnPoissonTail( exposure[j], m );
                else if ( m > 0 )        term += m * log( exposure[j] ) - RbMath::lnFactorial( int(m) );
            }
            if ( RbMath::isFinite(term) ) terms.push_back( term );
        }
    }

    if ( terms.empty() ) return RbConstants::Double::neginf;

    double result = log( mixture_rate * fossil[findIndex(o)] ) + RbMath::log_sum_exp( terms );
    if ( single == false ) result += log( mixture_rate * fossil[findIndex(y)] );

    return result;
}


/**
 * return the index i so that t_{i-1} > t >= t_i
 * where t_i is the instantaneous sampling time (i = 0,...,l)
 * t_0 is origin
 * t_l = 0.0
 */
size_t AbstractFossilizedBirthDeathRangeProcess::findIndex(double t) const
{
    return std::prev(std::upper_bound( times.begin(), times.end(), t)) - times.begin();
}


/** log q_tilde from the present up to t. */
double AbstractFossilizedBirthDeathRangeProcess::lnQTildeFromPresent( size_t i, double t ) const
{
    double s = q( i, t, true );

    for ( size_t j = 0; j < i; ++j ) s += q_tilde_i[j];

    return s;
}


/** log q from the present up to t. */
double AbstractFossilizedBirthDeathRangeProcess::lnQFromPresent( size_t i, double t ) const
{
    double s = q( i, t );

    for ( size_t j = 0; j < i; ++j ) s += q_i[j];

    return s;
}


/**
 * log P(species i originated above t), for t between its oldest appearance and its birth.
 *
 * A range's birth time is where it attaches to the nearest range the record samples, not where
 * the species itself began: ancestors that left no occurrence were never named, so a chain of
 * them may stand in between. The stretch holding that chain is the regime in which either
 * daughter of a budding event may carry the species identity, and below the species' own
 * origination the identity is pinned. Those two regimes are q and q_tilde, so the ratio of
 * their segment terms is the survival function
 *
 *     S(t) = q(o) q_tilde(t) / ( q_tilde(o) q(t) ),
 *
 * falling from one at the oldest appearance o. The remainder S(b) is the atom on the
 * attachment itself, where no unnamed ancestor intervenes.
 *
 * Symmetric speciation and anagenesis reach an unnamed ancestor by routes this does not count,
 * so a process admitting either must override it. So must one under complete lineage sampling,
 * where nothing collapses and the origination is the birth time itself.
 */
double AbstractFossilizedBirthDeathRangeProcess::lnOriginationSurvival( size_t i, double t ) const
{
    double o = ranges[i].first;

    size_t oi = findIndex(o);
    size_t ti = findIndex(t);

    return lnQFromPresent(oi, o) - lnQFromPresent(ti, t)
           + lnQTildeFromPresent(ti, t) - lnQTildeFromPresent(oi, o);
}


/**
 * Draw species i's own origination time.
 *
 * Nothing else in the model conditions on it: the occurrences all sit at or below the
 * oldest appearance, and gamma attaches to the birth time. So this is the exact full
 * conditional and the draw always stands, rather than a proposal needing a Hastings
 * ratio. Inverting the survival function needs no quadrature, as it is a ratio of
 * kernels the process already computes.
 */
double AbstractFossilizedBirthDeathRangeProcess::drawOrigination( size_t i ) const
{
    double o = ranges[i].first;
    double b = ranges[i].divergence;

    double ln_u = log( GLOBAL_RNG->uniform01() );

    // the atom at the attachment, where the sampled ancestor is the parent
    if ( ln_u <= lnOriginationSurvival(i, b) ) return b;

    // S decreases over (o, b), so bisect it
    double lo = o;
    double hi = b;

    for (size_t n = 0; n < 200 && hi - lo > 1E-12*(1.0 + b); ++n)
    {
        double mid = 0.5 * (lo + hi);

        if ( lnOriginationSurvival(i, mid) > ln_u ) lo = mid; else hi = mid;
    }

    return 0.5 * (lo + hi);
}



/**
 * The conditional mean of species i's origination time.
 *
 * The origination is confined to the oldest appearance and the divergence, so its mean is
 * the oldest appearance plus the integral of the survival function across them. That
 * carries none of the simulation noise drawOrigination does, and the atom needs no term of
 * its own because the survival function already holds it.
 *
 * The kernels are smooth within a rate interval and kinked only at its boundaries, so the
 * quadrature runs interval by interval. It inherits any override of the survival function.
 */
double AbstractFossilizedBirthDeathRangeProcess::meanOrigination( size_t i ) const
{
    double o = ranges[i].first;
    double b = ranges[i].divergence;

    size_t oi = findIndex(o);
    size_t bi = findIndex(b);

    double integral = 0.0;

    for ( size_t j = oi; j <= bi; ++j )
    {
        double lo = ( j == oi ) ? o : times[j];
        double hi = ( j == bi ) ? b : times[j+1];

        if ( hi <= lo ) continue;

        // composite Simpson over the piece
        size_t steps = 64;
        double h = (hi - lo) / steps;

        double sum = exp( lnOriginationSurvival(i, lo) ) + exp( lnOriginationSurvival(i, hi) );

        for ( size_t k = 1; k < steps; ++k )
        {
            sum += ( k % 2 == 1 ? 4.0 : 2.0 ) * exp( lnOriginationSurvival(i, lo + k*h) );
        }

        integral += sum * h / 3.0;
    }

    return o + integral;
}


/** p_i(t) */
double AbstractFossilizedBirthDeathRangeProcess::p( size_t i, double t, bool survival ) const
{
    // get the parameters
    double b = birth[i];
    double d = death[i];
    double f = survival ? 0.0 : fossil[i];
    double pi = survival ? pS_i[i] : p_i[i];
    double r = (i == 0 ? homogeneous_rho->getValue() : 0.0);
    double ti = times[i];
    
    double diff = b - d - f;
    double dt   = t - ti;

    double A = sqrt( diff*diff + 4.0*b*f);

    // survival takes f = 0, so b == d leaves A = 0, where the closed form is singular
    if ( A < 1E-10 )
    {
        double s = 1.0 - (1.0 - r)*pi;

        return 1.0 - s/(1.0 + b*s*dt);
    }

    double B = ( (1.0 - 2.0*(1.0-r)*pi )*b + d + f ) / A;

    double ln_e = -A*dt;

    double tmp = (1.0 + B) + exp(ln_e)*(1.0 - B);

    return (b + d + f - A * ((1.0+B)-exp(ln_e)*(1.0-B))/tmp)/(2.0*b);
}


/**
 * q_i(t)
 */
double AbstractFossilizedBirthDeathRangeProcess::q( size_t i, double t, bool tilde ) const
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

    if (tilde) q = 0.5 * (q - (b+d+f)*dt);
    
    return q;
}


/**
 *
 *
 */
/** The appearances and the origination and extinction times. */
void AbstractFossilizedBirthDeathRangeProcess::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, RbVector<double> &rv) const
{
    if ( n == "getFirstAppearances" || n == "getLastAppearances" || n == "getDivergenceTimes" || n == "getExtinctionTimes" )
    {
        AbstractFossilizedBirthDeathRangeProcess *self = const_cast<AbstractFossilizedBirthDeathRangeProcess *>( this );
        self->updateRanges();

        if ( n == "getExtinctionTimes" && marginalizesExtinction() == true )
        {
            throw RbException("getExtinctionTimes is unavailable when extended=false: the extinction times are marginalized out rather than sampled. Use getLastAppearances for the youngest occurrence ages, or extended=true to sample extinction times.");
        }

        rv.clear();
        for (size_t i = 0; i < ranges.size(); i++)
        {
            rv.push_back( n == "getFirstAppearances" ? ranges[i].first :
                        ( n == "getLastAppearances"  ? ranges[i].last  :
                        ( n == "getDivergenceTimes"   ? ranges[i].divergence : ranges[i].death ) ) );
        }
    }
    else if ( n == "getExtantStatus" || n == "getSingletonStatus" )
    {
        AbstractFossilizedBirthDeathRangeProcess *self = const_cast<AbstractFossilizedBirthDeathRangeProcess *>( this );
        self->updateRanges();

        // one indicator per sample, so it is the posterior mean that gives a taxon's probability of
        // having survived to the present without being sampled, or of tau_1 and tau_K being one age
        rv.clear();
        for (size_t i = 0; i < ranges.size(); i++)
        {
            bool flag = ( n == "getExtantStatus" ? ranges[i].survived : ranges[i].singleton );
            rv.push_back( flag ? 1.0 : 0.0 );
        }
    }
    else if ( n == "getOriginationTimes" )
    {
        AbstractFossilizedBirthDeathRangeProcess *self = const_cast<AbstractFossilizedBirthDeathRangeProcess *>( this );
        self->updateRanges();
        prepareProbComputation();

        // mean=TRUE integrates the conditional rather than drawing from it, which costs a
        // quadrature and returns the same quantity carrying none of the simulation noise
        bool mean = args.size() > 0 && static_cast<const TypedDagNode<Boolean>* >( args[0] )->getValue();

        rv.clear();
        for (size_t i = 0; i < ranges.size(); i++)
        {
            rv.push_back( mean ? meanOrigination(i) : drawOrigination(i) );
        }
    }
    else if ( n == "getPsiMixtureComponents" )
    {
        rv.clear();
        for (size_t i = 0; i < ranges.size(); i++)
        {
            rv.push_back( ranges[i].component + 1 );
        }
    }
    else
    {
        throw RbException() << "The fossilized birth death range process does not have a member method called '" << n << "'.";
    }
}




/**
 * The origin of the process, which is the oldest birth.
 */
void AbstractFossilizedBirthDeathRangeProcess::executeMethod(const std::string &n, const std::vector<const DagNode *> &args, double &rv) const
{
    if ( n == "getOrigin" )
    {
        AbstractFossilizedBirthDeathRangeProcess *self = const_cast<AbstractFossilizedBirthDeathRangeProcess *>( this );
        self->updateRanges();

        rv = origin;
    }
    else
    {
        throw RbException() << "The fossilized birth death range process does not have a member method called '" << n << "'.";
    }
}


/** Draw taxon i's appearances, nested within its range and its reported bins. */
void AbstractFossilizedBirthDeathRangeProcess::initializeFirstLast(size_t i)
{
    RandomNumberGenerator* rng = GLOBAL_RNG;

    bool augment_youngest = ranges[i].singleton == false;

    double lo, hi;

    // within its reported bin, and above the death unless that is marginalized
    if ( augment_youngest )
    {
        lo = marginalizesExtinction() ? taxa[i].getMinAge() : std::max( ranges[i].death, taxa[i].getMinAge() );
        hi = std::min( ranges[i].last_max, ranges[i].divergence );
        ranges[i].last = ( hi > lo ) ? rng->uniform01()*(hi - lo) + lo : lo;
    }
    else
    {
        ranges[i].last = ranges[i].death;
    }

    // both appearances sit under the birth
    lo = std::max( ranges[i].last, ranges[i].first_min );
    hi = std::min( ranges[i].first_max, ranges[i].divergence );
    ranges[i].first = ( hi > lo ) ? rng->uniform01()*(hi - lo) + lo : lo;

    // a single occurrence is its own youngest
    if ( augment_youngest == false ) ranges[i].last = ranges[i].first;

    // a non-extended extinct range ends at its youngest appearance
    if ( marginalizesExtinction() == true && taxa[i].isExtinct() == true )
    {
        ranges[i].death = ranges[i].last;
        ranges[i].survived = ( ranges[i].last == times.front() );
    }

    dirty_taxa[i] = true;
}

/** Without mvStratigraphicRange the appearances never move; warn once. */
void AbstractFossilizedBirthDeathRangeProcess::warnIfNoResampleMove( void ) const
{
    if ( has_resample_move == false && warned_no_resample == false )
    {
        warned_no_resample = true;
        RBOUT("Warning: no mvStratigraphicRange move. Appearances will not be sampled.");
    }
}


/** Without a dnFossilRecord the sampling is unscored and psi has no data; warn once. */
void AbstractFossilizedBirthDeathRangeProcess::warnIfNoReportingNode( void ) const
{
    if ( has_reporting_node == false && warned_no_reporting == false )
    {
        warned_no_reporting = true;
        RBOUT("Warning: no dnFossilRecord node. Fossil sampling is not scored and psi has no data.");
    }
}


// At rho = 0, condition="nTaxa" is condition="time".
void AbstractFossilizedBirthDeathRangeProcess::warnIfInertNTaxaCondition( void ) const
{
    if ( condition == "nTaxa" && homogeneous_rho->getValue() == 0.0 && warned_inert_ntaxa == false )
    {
        warned_inert_ntaxa = true;
        RBOUT("Warning: condition=\"nTaxa\" is inert at rho=0, where no survivor is observed.");
    }
}


/** Below rho = 1 only mvExtinctionRJSwitch moves a taxon on or off the present; warn once. */
void AbstractFossilizedBirthDeathRangeProcess::warnIfNoSurvivalMove( void ) const
{
    if ( has_survival_move == true || warned_no_survival == true ) return;
    if ( marginalizesExtinction() == true || homogeneous_rho->getValue() == 1.0 ) return;

    for (size_t i = 0; i < taxa.size(); i++)
    {
        if ( taxa[i].isExtinct() == true )
        {
            warned_no_survival = true;
            RBOUT("Warning: rho < 1 and no mvExtinctionRJSwitch move. Survival will not be sampled.");
            return;
        }
    }
}


/** Under presence one reported occurrence may stand for more, which only a move over singleton status samples; warn once. */
void AbstractFossilizedBirthDeathRangeProcess::warnIfNoSingletonMove( void ) const
{
    if ( has_singleton_move == true || warned_no_singleton == true ) return;
    if ( record_reporting != PRESENCE ) return;

    for ( size_t i = 0; i < taxa.size(); i++ )
    {
        const auto &ages = taxa[i].getOccurrences();

        // an exact age leaves nothing to separate
        if ( ranges[i].singleton == true && ages.empty() == false && ages.begin()->first.getMin() != ages.begin()->first.getMax() )
        {
            warned_no_singleton = true;
            RBOUT("Warning: no mvSingletonRJSwitch move. Single occurrences stay single ages.");
            return;
        }
    }
}



/** Flip taxon i between surviving and dying. The caller sets the death time to match. */
void AbstractFossilizedBirthDeathRangeProcess::switchSurvived(size_t i)
{
    switched_range   = i;
    stored_survived  = ranges[i].survived;
    stored_singleton = ranges[i].singleton;
    switched         = true;

    ranges[i].survived = !ranges[i].survived;

    dirty_taxa[i] = true;   // the status is in no value
}


void AbstractFossilizedBirthDeathRangeProcess::resampleFirstLast(size_t i)
{
    stored_range = i;
    stored_first = ranges[i].first;
    stored_last  = ranges[i].last;
    resampled = true;

    dirty_taxa[i] = true;   // a tree does not carry these ages

    // a non-extended extinct range ends at its tip
    if ( marginalizesExtinction() == true && taxa[i].isExtinct() == true && ranges[i].singleton )
    {
        return;
    }

    // data-fixed support
    double hi = ranges[i].first_max;
    ranges[i].first = GLOBAL_RNG->uniform01()*(ranges[i].first_min - hi) + hi;

    // that tip is tau_K
    if ( marginalizesExtinction() == true && taxa[i].isExtinct() == true )
    {
        return;
    }

    // a single reported occurrence is its own youngest
    if ( ranges[i].singleton == false )
    {
        ranges[i].last = GLOBAL_RNG->uniform01()*(ranges[i].last_max - taxa[i].getMinAge()) + taxa[i].getMinAge();
    }
    else
    {
        ranges[i].last = ranges[i].first;
    }

}


/**
 * Expand taxon i's tau_1 uniformly above tau_K, within its interval and below its divergence, or
 * collapse it onto tau_K. tau_K stays, so the move is the same whether or not a value carries it.
 *
 * \return The log hastings ratio.
 */
double AbstractFossilizedBirthDeathRangeProcess::switchSingleton(size_t i)
{
    double y     = ranges[i].last;
    double width = std::min( taxa[i].getMaxAge(), ranges[i].divergence ) - y;

    if ( width <= 0.0 ) return RbConstants::Double::neginf;

    stored_range = i;
    stored_first = ranges[i].first;
    stored_last  = y;
    resampled    = true;

    switched_range   = i;
    stored_survived  = ranges[i].survived;
    stored_singleton = ranges[i].singleton;
    switched         = true;

    bool single = ranges[i].singleton;
    ranges[i].singleton = !single;
    ranges[i].first = single ? y + GLOBAL_RNG->uniform01() * width : y;

    dirty_taxa[i] = true;

    return single ? log( width ) : -log( width );
}


void AbstractFossilizedBirthDeathRangeProcess::drawRanges()
{
    RandomNumberGenerator* rng = GLOBAL_RNG;

    // an unbounded max_age would push every birth to infinity, so fall back to the minima
    double max = 0;
    for (size_t i = 0; i < taxa.size(); i++)
    {
        double o = taxa[i].getMaxAge();
        if ( RbMath::isFinite(o) == false )
        {
            o = 0.0;
            const std::map<TimeInterval, size_t> &record = taxa[i].getOccurrences();
            for ( std::map<TimeInterval, size_t>::const_iterator Fi = record.begin(); Fi != record.end(); Fi++ )
            {
                o = std::max( o, Fi->first.getMin() );
            }
        }
        if ( o > max ) max = o;
    }
    max *= 1.1;
    if ( max == 0.0 ) max = 1.0;

    double present = times.front();

    // The origin comes first and every appearance is drawn under it.
    if ( origin_age != NULL )
    {
        max = origin_age->getValue();
    }

    for (size_t i = 0; i < taxa.size(); i++)
    {
        // Draw d over its range, then the appearances NESTED (d <= last <= first)
        if ( taxa[i].isExtinct() == false )
        {
            ranges[i].survived = true;
            ranges[i].death = present;
        }
        else
        {
            // survives to the present with exp(-mu L) and is missed there with 1-rho
            double cap = std::min( ranges[i].last_max, max );
            double L = cap - present;
            double m = death[ findIndex( cap ) ];
            bool   exponential = ( m > 0.0 && RbMath::isFinite(L) == true );

            double p_survived = 0.0;
            if ( exponential == true )
            {
                double reach = exp( -m * L );
                p_survived = ( 1.0 - homogeneous_rho->getValue() ) * reach;
                p_survived /= p_survived + 1.0 - reach;
            }

            if ( rng->uniform01() < p_survived )
            {
                ranges[i].survived = true;
                ranges[i].death = present;
            }
            else if ( exponential == true )
            {
                // the same exponential, truncated to the deaths the record allows
                double u = rng->uniform01();
                ranges[i].survived = false;
                ranges[i].death = cap + log( 1.0 - u*( 1.0 - exp( -m * L ) ) ) / m;
            }
            else
            {
                ranges[i].survived = false;
                ranges[i].death = rng->uniform01()*L + present;
            }
        }
        ranges[i].divergence = max;

        initializeFirstLast(i);
    }

    // place births oldest-first, each inside a lineage already placed and still alive at it
    std::vector<size_t> order( taxa.size() );
    for (size_t i = 0; i < taxa.size(); i++) order[i] = i;

    std::sort( order.begin(), order.end(), [&](size_t a, size_t b) { return ranges[a].first > ranges[b].first; } );

    for (size_t k = 0; k < order.size(); k++)
    {
        size_t i = order[k];

        // the oldest birth is the origin, already drawn above and bounding every appearance
        if ( k == 0 )
        {
            ranges[i].divergence = ( origin_age != NULL )
                              ? max
                              : ranges[i].first + rng->uniform01()*(max - ranges[i].first);

            continue;
        }

        size_t pick = size_t( rng->uniform01()*k );
        if ( pick >= k ) pick = k - 1;

        size_t a = order[pick];

        // sorted oldest-first, so the window is never empty
        double lo = std::max( ranges[i].first, ranges[a].death );

        ranges[i].divergence = lo + rng->uniform01()*(ranges[a].divergence - lo);
    }
}


void AbstractFossilizedBirthDeathRangeProcess::keepSpecialization(const DagNode *toucher)
{
    dirty_taxa = std::vector<bool>(taxa.size(), false);

    resampled = false;
    switched = false;
    reallocated = false;
    touched = false;
    weld_origin = true;
}


void AbstractFossilizedBirthDeathRangeProcess::restoreSpecialization(const DagNode *toucher)
{
    partial_likelihood = stored_likelihood;

    ranges = stored_ranges;

    if ( resampled )
    {
        ranges[stored_range].first = stored_first;
        ranges[stored_range].last  = stored_last;
    }

    if ( switched )
    {
        ranges[switched_range].survived  = stored_survived;
        ranges[switched_range].singleton = stored_singleton;
    }
    if ( reallocated )
    {
        ranges[reallocated_range].component = stored_component;
    }

    dirty_taxa = std::vector<bool>(taxa.size(), false);

    resampled = false;
    switched = false;
    reallocated = false;
    touched = false;
    weld_origin = true;
}


void AbstractFossilizedBirthDeathRangeProcess::touchSpecialization(const DagNode *toucher, bool touchAll)
{
    if ( touched == false )
    {
        stored_likelihood = partial_likelihood;
        stored_ranges     = ranges;

        dirty_taxa = std::vector<bool>(taxa.size(), true);
    }

    touched = true;
}


/**
 *
 *
 */
void AbstractFossilizedBirthDeathRangeProcess::prepareProbComputation( void ) const
{
    if ( homogeneous_lambda != NULL )
    {
        birth = std::vector<double>(num_intervals, homogeneous_lambda->getValue() );
    }
    else
    {
        birth = heterogeneous_lambda->getValue();
    }
    if ( homogeneous_mu != NULL )
    {
        death = std::vector<double>(num_intervals, homogeneous_mu->getValue() );
    }
    else
    {
        death = heterogeneous_mu->getValue();
    }
    if ( homogeneous_psi != NULL )
    {
        fossil = std::vector<double>(num_intervals, homogeneous_psi->getValue() );
    }
    else
    {
        fossil = heterogeneous_psi->getValue();
    }

    if ( timeline != NULL )
    {
        times = timeline->getValue();
    }
    else
    {
        times.clear();
    }

    if ( times.size() < num_intervals )
    {
        times.insert(times.begin(), present_time);
    }

    for (size_t i = 0; i < num_intervals; i++)
    {
        double ti = times[i];
        double b = birth[i];
        double d = death[i];
        double f = fossil[i];

        if (i < num_intervals-1)
        {
            double r = (i == 0 ? homogeneous_rho->getValue() : 0.0);
            double t = times[i+1];

            double diff = b - d - f;
            double dt   = t - ti;

            double A = sqrt( diff*diff + 4.0*b*f);
            double B = ( (1.0 - 2.0*(1.0-r)*p_i[i] )*b + d + f ) / A;

            double ln_e = -A*dt;

            double tmp = (1.0 + B) + exp(ln_e)*(1.0 - B);

            q_i[i]       = log(4.0) + ln_e - 2.0*log(tmp);
            q_tilde_i[i] = 0.5 * ( q_i[i] - (b+d+f)*dt );
            p_i[i+1]       = (b + d + f - A * ((1.0+B)-exp(ln_e)*(1.0-B))/tmp)/(2.0*b);

            if ( condition == "survival" )
            {
                // survival ignores fossil sampling: the f = 0 case of the recursion above
                diff = b - d;

                A = fabs( diff );

                double s = 1.0 - (1.0 - r)*pS_i[i];

                if ( A < 1E-10 )
                {
                    // b == d leaves A = 0, where the closed form is singular but its limit is not
                    pS_i[i+1] = 1.0 - s/(1.0 + b*s*dt);
                }
                else
                {
                    B = ( (1.0 - 2.0*(1.0-r)*pS_i[i] )*b + d ) / A;

                    ln_e = -A*dt;

                    tmp = (1.0 + B) + exp(ln_e)*(1.0 - B);

                    pS_i[i+1] = (b + d - A * ((1.0+B)-exp(ln_e)*(1.0-B))/tmp)/(2.0*b);
                }
            }
        }
    }

    if ( num_components > 1 )
    {
        size_t H = num_components;

        mixture_rates = psi_mixture_rates->getValue();
        if ( mixture_rates.size() != H )
        {
            throw RbException("psiMixtureRates cannot change its number of components.");
        }
        mixture_weights = ( psi_mixture_weights != NULL ) ? std::vector<double>( psi_mixture_weights->getValue() ) : std::vector<double>(H, 1.0 / H);

        mixture_p_i[0]   = std::vector<double>(H, 1.0 - homogeneous_rho->getValue());
        lambda_pbar_i[0] = 0.0;

        for (size_t i = 0; i + 1 < num_intervals; i++)
        {
            std::vector<double> p = mixture_p_i[i];
            std::vector<double> Q(H*H, 0.0), Q_scale(H, 0.0);
            for (size_t c = 0; c < H; c++) Q[c*H + c] = 1.0;

            double lambda_pbar = lambda_pbar_i[i];
            propagateMixture( i, p, times[i+1] - times[i], lambda_pbar, &Q, &Q_scale );

            mixture_p_i[i+1]     = p;
            lambda_pbar_i[i+1]   = lambda_pbar;
            mixture_q_i[i]       = Q;
            mixture_q_scale_i[i] = Q_scale;
        }
    }
}


/**
 * Carry p of each component, int lambda pbar and the columns of Q a time dt up interval j.
 * With p = y / z, y' = mu z - (lambda + mu + nu_h psi) y and z' = -lambda sum_h w_h y_h are linear,
 * and Q = dp(t)/dp(s). Q holds columns of length H, each scaled by exp(Q_scale).
 */
void AbstractFossilizedBirthDeathRangeProcess::propagateMixture( size_t j, std::vector<double> &p, double dt, double &lambda_pbar, std::vector<double> *Q, std::vector<double> *Q_scale ) const
{
    if ( dt <= 0.0 ) return;

    size_t H = num_components;
    size_t N = H + 1;

    double b = birth[j];
    double d = death[j];
    double f = fossil[j];

    double min_rate = *std::min_element( mixture_rates.begin(), mixture_rates.end() );

    // the two slowest eigenvalues lie within lambda + mu + min_rate psi of each other, which bounds
    // the cancellation in each step's Q; faster components only decay within a step
    size_t n = std::max( size_t(1), size_t( ceil( dt * (b + 2.0*d + min_rate*f) ) ) );
    double dt_step = dt / n;

    // the generator of (y, z) over one step, with z last: diagonal on y, a column mu and a row -lambda w
    std::vector<double> diag(H, 0.0), row(H, 0.0);
    double norm = 0.0, row_norm = 0.0;
    for (size_t h = 0; h < H; h++)
    {
        diag[h] = -(b + d + mixture_rates[h]*f) * dt_step;
        row[h]  = -b * mixture_weights[h] * dt_step;
        norm = std::max( norm, fabs(diag[h]) + d * dt_step );
        row_norm += fabs(row[h]);
    }
    norm = std::max( norm, row_norm );
    double column = d * dt_step;

    // scale the generator to norm at most 1/2, so the series converges fast, and square back after
    int squarings = ( norm > 0.5 ) ? int( ceil( log2( norm / 0.5 ) ) ) : 0;
    double scale = ldexp( 1.0, -squarings );
    for (size_t h = 0; h < H; h++)
    {
        diag[h] *= scale;
        row[h]  *= scale;
    }
    column *= scale;

    // its exponential by the Taylor series, row-major, each term costing O(H^2) by the sparsity
    std::vector<double> E(N*N, 0.0), term(N*N, 0.0), next(N*N, 0.0);
    for (size_t a = 0; a < N; a++)
    {
        E[a*N + a]    = 1.0;
        term[a*N + a] = 1.0;
    }
    for (size_t m = 1; m < 30; m++)
    {
        double largest = 0.0;
        for (size_t a = 0; a < N; a++)
        {
            const double *ta = &term[a*N];
            double to_z = 0.0;
            for (size_t c = 0; c < H; c++)
            {
                double x = ( ta[c] * diag[c] + ta[H] * row[c] ) / m;
                next[a*N + c] = x;
                to_z += ta[c];
                largest = std::max( largest, fabs(x) );
            }
            next[a*N + H] = to_z * column / m;
            largest = std::max( largest, fabs(next[a*N + H]) );
        }
        term.swap(next);
        for (size_t a = 0; a < N*N; a++) E[a] += term[a];
        if ( largest < 1E-18 ) break;
    }
    for (int k = 0; k < squarings; k++)
    {
        for (size_t a = 0; a < N; a++)
        {
            for (size_t c = 0; c < N; c++)
            {
                double x = 0.0;
                for (size_t e = 0; e < N; e++) x += E[a*N + e] * E[e*N + c];
                next[a*N + c] = x;
            }
        }
        E.swap(next);
    }

    size_t num_columns = ( Q == NULL ) ? 0 : Q->size() / H;

    std::vector<double> v(N, 0.0), S(H*H, 0.0), col(H, 0.0);

    for (size_t step = 0; step < n; step++)
    {
        for (size_t a = 0; a < N; a++)
        {
            double x = E[a*N + H];
            for (size_t e = 0; e < H; e++) x += E[a*N + e] * p[e];
            v[a] = x;
        }

        double z = v[H];
        lambda_pbar -= log(z);
        for (size_t a = 0; a < H; a++) p[a] = v[a] / z;

        if ( num_columns == 0 ) continue;

        for (size_t a = 0; a < H; a++)
        {
            for (size_t e = 0; e < H; e++)
            {
                // Q is nonnegative, so a negative entry is rounding in a fast component's column
                S[a*H + e] = std::max( 0.0, ( E[a*N + e] - p[a] * E[H*N + e] ) / z );
            }
        }

        for (size_t c = 0; c < num_columns; c++)
        {
            double *Qc = &(*Q)[c*H];
            double largest = 0.0;
            for (size_t a = 0; a < H; a++)
            {
                double x = 0.0;
                for (size_t e = 0; e < H; e++) x += S[a*H + e] * Qc[e];
                col[a] = x;
                largest = std::max( largest, fabs(x) );
            }
            if ( largest > 0.0 )
            {
                for (size_t a = 0; a < H; a++) Qc[a] = col[a] / largest;
                (*Q_scale)[c] += log(largest);
            }
            else
            {
                std::fill( Qc, Qc + H, 0.0 );
            }
        }
    }
}


/** p of each component and int lambda pbar from the present up to t. */
void AbstractFossilizedBirthDeathRangeProcess::mixtureAt( double t, std::vector<double> &p, double &lambda_pbar ) const
{
    size_t j = findIndex(t);

    p           = mixture_p_i[j];
    lambda_pbar = lambda_pbar_i[j];

    propagateMixture( j, p, t - times[j], lambda_pbar, NULL, NULL );
}


/**
 * Birth-death density of range i given its component h, for every h or for h = only.
 * Between b_i and tau_1 the lineage may pass through unsampled species of any component, so Q
 * enters as sum_e w_e Q[e][h]. Below tau_1, Q_tilde = exp( -(lambda + mu + nu_h psi) t + int lambda pbar ).
 */
void AbstractFossilizedBirthDeathRangeProcess::mixtureRangeLnProbs( size_t i, std::vector<double> &ln_probs, size_t only ) const
{
    size_t H = num_components;

    double b = ranges[i].divergence;
    double d = ranges[i].death;
    double o = ranges[i].first;

    size_t bi = findIndex(b);
    size_t oi = findIndex(o);
    size_t di = findIndex(d);

    std::vector<size_t> columns;
    for (size_t h = 0; h < H; h++)
    {
        if ( only >= H || h == only ) columns.push_back(h);
    }
    size_t num_columns = columns.size();

    std::vector<double> p;
    double lambda_pbar_o = 0.0;
    mixtureAt( o, p, lambda_pbar_o );

    std::vector<double> Q(H*num_columns, 0.0), Q_scale(num_columns, 0.0);
    for (size_t c = 0; c < num_columns; c++) Q[c*H + columns[c]] = 1.0;

    double lambda_pbar = 0.0;
    if ( oi == bi )
    {
        propagateMixture( oi, p, b - o, lambda_pbar, &Q, &Q_scale );
    }
    else
    {
        propagateMixture( oi, p, times[oi+1] - o, lambda_pbar, &Q, &Q_scale );

        // whole intervals, whose columns carry their own scales
        std::vector<double> col(H, 0.0);
        for (size_t j = oi+1; j < bi; j++)
        {
            const std::vector<double> &Qj = mixture_q_i[j];
            const std::vector<double> &sj = mixture_q_scale_i[j];

            for (size_t c = 0; c < num_columns; c++)
            {
                double *Qc = &Q[c*H];

                double top = RbConstants::Double::neginf;
                for (size_t e = 0; e < H; e++)
                {
                    if ( Qc[e] != 0.0 ) top = std::max( top, sj[e] );
                }
                if ( top == RbConstants::Double::neginf ) continue;

                std::fill( col.begin(), col.end(), 0.0 );
                for (size_t e = 0; e < H; e++)
                {
                    if ( Qc[e] == 0.0 ) continue;
                    double x = Qc[e] * exp( sj[e] - top );
                    for (size_t a = 0; a < H; a++) col[a] += Qj[e*H + a] * x;
                }

                double largest = 0.0;
                for (size_t a = 0; a < H; a++) largest = std::max( largest, fabs(col[a]) );
                if ( largest == 0.0 )
                {
                    std::fill( Qc, Qc + H, 0.0 );
                    continue;
                }
                for (size_t a = 0; a < H; a++) Qc[a] = col[a] / largest;
                Q_scale[c] += top + log(largest);
            }
        }

        p = mixture_p_i[bi];
        propagateMixture( bi, p, b - times[bi], lambda_pbar, &Q, &Q_scale );
    }

    ln_probs = std::vector<double>(H, RbConstants::Double::neginf);

    for (size_t c = 0; c < num_columns; c++)
    {
        double s = 0.0;
        for (size_t e = 0; e < H; e++) s += mixture_weights[e] * Q[c*H + e];

        ln_probs[ columns[c] ] = log( birth[bi] ) + log(s) + Q_scale[c];
    }

    // extant taxa with no fossil samples end at their oldest appearance, the present
    if ( taxa[i].getMaxAge() == times.front() ) return;

    double lambda_pbar_d = 0.0;
    mixtureAt( d, p, lambda_pbar_d );

    // the exposure from d_i up to tau_1, without and with psi
    double birth_death = 0.0, sampling = 0.0;
    for (size_t j = di; j <= oi; j++)
    {
        double lo = std::max( d, times[j] );
        double hi = ( j + 1 < num_intervals ) ? std::min( o, times[j+1] ) : o;
        if ( hi <= lo ) continue;

        birth_death += ( birth[j] + death[j] ) * (hi - lo);
        sampling    += fossil[j] * (hi - lo);
    }

    for (size_t c = 0; c < num_columns; c++)
    {
        size_t h = columns[c];

        ln_probs[h] += -( birth_death + mixture_rates[h]*sampling ) + ( lambda_pbar_o - lambda_pbar_d );

        if ( ranges[i].survived == false ) ln_probs[h] += log( death[di] );
    }
}


/** Taxon i's log-density under each component: its range, and its record where one is scored. */
std::vector<double> AbstractFossilizedBirthDeathRangeProcess::psiMixtureComponentLnProbs( size_t i ) const
{
    prepareProbComputation();

    std::vector<double> ln_probs;
    mixtureRangeLnProbs( i, ln_probs, num_components );

    if ( has_reporting_node == true && taxa[i].getMaxAge() != times.front() )
    {
        for (size_t h = 0; h < num_components; h++)
        {
            ln_probs[h] += computeLnFossilRecord( i, mixture_rates[h] );
        }
    }

    return ln_probs;
}


/** Move taxon i to component h. The caller touches the node. */
void AbstractFossilizedBirthDeathRangeProcess::setPsiMixtureComponent( size_t i, size_t h )
{
    reallocated_range = i;
    stored_component  = ranges[i].component;
    reallocated       = true;

    ranges[i].component = h;

    dirty_taxa[i] = true;   // the component is in no value
}


/**
 * Swap the parameters held by this distribution.
 * 
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void AbstractFossilizedBirthDeathRangeProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if (oldP == heterogeneous_lambda)
    {
        heterogeneous_lambda = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == heterogeneous_mu)
    {
        heterogeneous_mu = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == heterogeneous_psi)
    {
        heterogeneous_psi = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == homogeneous_lambda)
    {
        homogeneous_lambda = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == homogeneous_mu)
    {
        homogeneous_mu = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == homogeneous_psi)
    {
        homogeneous_psi = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == homogeneous_rho)
    {
        homogeneous_rho = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == timeline)
    {
        timeline = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == origin_age)
    {
        origin_age = static_cast<const TypedDagNode<double>* >( newP );
    }
    else if (oldP == psi_mixture_rates)
    {
        psi_mixture_rates = static_cast<const TypedDagNode< RbVector<double> >* >( newP );
    }
    else if (oldP == psi_mixture_weights)
    {
        psi_mixture_weights = static_cast<const TypedDagNode<Simplex>* >( newP );
    }
}
