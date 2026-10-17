#include "CompleteFossilizedBirthDeathRangeProcess.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "RbConstants.h"
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
                                                                                   const TypedDagNode<double> *inorigin,
                                                                                   double inpresent) :
    FossilizedBirthDeathRangeProcess(inspeciation, inextinction, inpsi, inrho, intimes, incondition, intaxa, complete_record, inorigin, inpresent)
{

}


CompleteFossilizedBirthDeathRangeProcess* CompleteFossilizedBirthDeathRangeProcess::clone( void ) const
{
    return new CompleteFossilizedBirthDeathRangeProcess( *this );
}


/** Without a given origin every birth is a speciation, as in PyRate. */
double CompleteFossilizedBirthDeathRangeProcess::originLnProb( void )
{
    return origin_age != NULL ? AbstractFossilizedBirthDeathRangeProcess::originLnProb() : 0.0;
}


/** Lineages are independent, so there is no coexistence factor. */
double CompleteFossilizedBirthDeathRangeProcess::computeLnProbability( void )
{
    return computeLnProbabilityRanges();
}


/** Birth-death density of range i, with no p(t). */
double CompleteFossilizedBirthDeathRangeProcess::rangeLnProb( size_t i )
{
    double b = ranges[i].birth;
    double d = ranges[i].death;

    double present = times.front();

    size_t bi = findIndex(b);
    size_t di = findIndex(d);

    // include speciation density
    double lnProb = log( birth[bi] );

    // include extinction density
    if ( ranges[i].survived == false ) lnProb += log( death[di] );

    double psi_b_d = 0.0;

    for ( size_t j = di; j <= bi; j++ )
    {
        double t_0 = ( j < num_intervals-1 ? times[j+1] : RbConstants::Double::inf );

        double dt = std::min(b, t_0) - std::max(d, times[j]);

        lnProb -= (birth[j] + death[j])*dt;

        psi_b_d += fossil[j]*dt;
    }

    // fossil non-sampling normalization over [d,b]
    lnProb -= psi_b_d;

    // an extant taxon with no fossil is observed by the present sample
    if ( condition == "sampling" && taxa[i].getMaxAge() != present )
    {
        lnProb -= log(-expm1(-psi_b_d));
    }

    return lnProb;
}
