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
                                                                                         bool require_parent,
                                                                                         const TypedDagNode<double> *inorigin,
                                                                                         double inpresent) :
    FossilizedBirthDeathRangeProcess(inspeciation, inextinction, inpsi, inrho, intimes, incondition, intaxa, complete_record, coexistence, require_parent, inorigin, inpresent)
{

}


CompleteFossilizedBirthDeathRangeProcess* CompleteFossilizedBirthDeathRangeProcess::clone( void ) const
{
    return new CompleteFossilizedBirthDeathRangeProcess( *this );
}


/** A range's birth time is the species' own origination. */
double CompleteFossilizedBirthDeathRangeProcess::lnOriginationSurvival( size_t i, double t ) const
{
    return 0.0;
}


/** The ranges carry the sampling conditioning; survival and the count at the present apply as usual. */
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


/** With p = 0 the kernel is plain exponential waiting. */
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


/** Under condition="sampling", divide range i by its own probability of leaving a fossil. */
double CompleteFossilizedBirthDeathRangeProcess::samplingConditionLnProb( size_t i ) const
{
    if ( condition != "sampling" ) return 0.0;

    if ( taxa[i].getAgeRange()->getMax() == times.front() ) return 0.0;

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
