#include "CompleteFossilizedBirthDeathRangeProcess.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "RbConstants.h"
#include "RbException.h"
#include "RbMathLogic.h"
#include "Taxon.h"

using namespace RevBayesCore;


CompleteFossilizedBirthDeathRangeProcess::CompleteFossilizedBirthDeathRangeProcess(const DagNode *inspeciation,
                                                                                         const DagNode *inextinction,
                                                                                         const DagNode *inpsi,
                                                                                         const TypedDagNode<double> *inrho,
                                                                                         const TypedDagNode< RbVector<double> > *intimes,
                                                                                         const std::string &incondition,
                                                                                         const std::vector<Taxon> &intaxa,
                                                                                         bool complete_record,
                                                                                         bool coexistence,
                                                                                         const TypedDagNode<double> *inorigin,
                                                                                         double inpresent) :
    FossilizedBirthDeathRangeProcess(inspeciation, inextinction, inpsi, inrho, intimes, incondition, intaxa, complete_record, coexistence, inorigin, inpresent)
{

}


CompleteFossilizedBirthDeathRangeProcess* CompleteFossilizedBirthDeathRangeProcess::clone( void ) const
{
    return new CompleteFossilizedBirthDeathRangeProcess( *this );
}


/**
 * Nothing collapses, so a range's birth time is the species' own origination and the
 * conditional on it is the point mass there.
 */
double CompleteFossilizedBirthDeathRangeProcess::lnOriginationSurvival( size_t i, double t ) const
{
    return 0.0;
}


/**
 * With p zero the probability of leaving a sample anywhere is one, so there is nothing for a
 * process-level normalization to be about and the ranges carry it instead. Survival and the
 * count at the present still describe the process that generated them, so they are left alone.
 */
double CompleteFossilizedBirthDeathRangeProcess::conditionLnProb( void ) const
{
    if ( condition == "sampling" ) return 0.0;

    return FossilizedBirthDeathRangeProcess::conditionLnProb();
}


/** Only a given start divides a speciation rate back out, as PyRate's oldest interval does. */
double CompleteFossilizedBirthDeathRangeProcess::originLnProb( void )
{
    if ( origin_age == NULL ) return 0.0;

    return FossilizedBirthDeathRangeProcess::originLnProb();
}


/**
 * With no unsampled descendants to allow for there is no p(t), so the kernel is plain
 * exponential waiting: a speciation rate where the range begins, an extinction rate where it
 * ends, and the three rates integrated over its duration.
 */
double CompleteFossilizedBirthDeathRangeProcess::rangeLnProb( size_t i )
{
    double b = ranges[i].divergence;
    double d = ranges[i].death;

    size_t bi = findIndex(b);
    size_t di = findIndex(d);

    // include speciation density
    double lnProb = log( birth[bi] );

    // include extinction density
    if ( ranges[i].survived == false ) lnProb += log( death[di] );

    for ( size_t j = di; j <= bi; j++ )
    {
        double top = ( j < num_intervals-1 ? times[j+1] : RbConstants::Double::inf );

        double dt = std::min(b, top) - std::max(d, times[j]);

        lnProb -= (birth[j] + death[j] + fossil[j])*dt;
    }

    return lnProb + samplingConditionLnProb( i );
}


/**
 * Range i's own conditioning on having been sampled.
 *
 * Nothing is integrated out here, so the conditioning has to be declared. Under
 * condition="time" the ranges are taken as the complete set and no term applies. Under
 * condition="sampling" a range reached the data only by leaving a fossil, and dropping the
 * ranges that left none changes no other range's density, so the conditioning is a division by
 * this range's own sampling probability. A species known only from an extant sample is in no
 * fossil record and takes no term either way.
 */
double CompleteFossilizedBirthDeathRangeProcess::samplingConditionLnProb( size_t i ) const
{
    if ( condition != "sampling" ) return 0.0;

    if ( taxa[i].getMaxAge() == times.front() ) return 0.0;

    double b = ranges[i].divergence;
    double d = ranges[i].death;

    size_t bi = findIndex(b);
    size_t di = findIndex(d);

    double sampling = 0.0;

    for ( size_t j = di; j <= bi; j++ )
    {
        double top = ( j < num_intervals-1 ? times[j+1] : RbConstants::Double::inf );

        sampling += fossil[j] * ( std::min(b, top) - std::max(d, times[j]) );
    }

    return -log( -expm1(-sampling) );
}
