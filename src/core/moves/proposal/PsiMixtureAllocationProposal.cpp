#include "PsiMixtureAllocationProposal.h"

#include <cstddef>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

#include "AbstractFossilizedBirthDeathRangeProcess.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RbMathLogic.h"

using namespace RevBayesCore;


PsiMixtureAllocationProposal::PsiMixtureAllocationProposal( StochasticNode<MatrixReal> *n ) : Proposal(),
    matrix( n )
{
    addNode( matrix );

    if ( dynamic_cast<AbstractFossilizedBirthDeathRangeProcess *>( &matrix->getDistribution() ) == NULL )
    {
        throw RbException("mvPsiMixtureAllocation requires a fossilized birth-death range process.");
    }
    if ( rangeProcess().getNumberOfPsiMixtureComponents() < 2 )
    {
        throw RbException("mvPsiMixtureAllocation requires a range process with psiMixtureRates.");
    }
}


/** Resolved on each use, so a swapped node leaves no stale pointer. */
AbstractFossilizedBirthDeathRangeProcess& PsiMixtureAllocationProposal::rangeProcess( void ) const
{
    return dynamic_cast<AbstractFossilizedBirthDeathRangeProcess &>( matrix->getDistribution() );
}


void PsiMixtureAllocationProposal::cleanProposal( void )
{
    ; // do nothing
}


PsiMixtureAllocationProposal* PsiMixtureAllocationProposal::clone( void ) const
{
    return new PsiMixtureAllocationProposal( *this );
}


const std::string& PsiMixtureAllocationProposal::getProposalName( void ) const
{
    static std::string name = "PsiMixtureAllocation";

    return name;
}


double PsiMixtureAllocationProposal::getProposalTuningParameter( void ) const
{
    return RbConstants::Double::nan;
}


/**
 * Draw one taxon's component from its full conditional.
 *
 * \return The hastings ratio, which cancels the posterior ratio.
 */
double PsiMixtureAllocationProposal::doProposal( void )
{
    AbstractFossilizedBirthDeathRangeProcess &process = rangeProcess();

    size_t n = process.getTaxa().size();
    if ( n == 0 ) return 0.0;

    RandomNumberGenerator* rng = GLOBAL_RNG;

    size_t i = size_t( rng->uniform01() * n );
    if ( i >= n ) i = n - 1;

    std::vector<double> ln_probs = process.psiMixtureComponentLnProbs( i );
    size_t H = ln_probs.size();

    for (size_t h = 0; h < H; h++)
    {
        if ( std::isnan( ln_probs[h] ) ) ln_probs[h] = RbConstants::Double::neginf;
    }

    double top = *std::max_element( ln_probs.begin(), ln_probs.end() );
    if ( RbMath::isFinite(top) == false ) return RbConstants::Double::neginf;

    std::vector<double> probs(H, 0.0);
    double total = 0.0;
    for (size_t h = 0; h < H; h++)
    {
        probs[h] = exp( ln_probs[h] - top );
        total += probs[h];
    }

    double u = rng->uniform01() * total;
    size_t h_new = 0;
    while ( h_new + 1 < H && u >= probs[h_new] )
    {
        u -= probs[h_new];
        h_new++;
    }

    size_t h_old = process.getPsiMixtureComponent( i );

    process.setPsiMixtureComponent( i, h_new );

    return ln_probs[h_old] - ln_probs[h_new];
}


void PsiMixtureAllocationProposal::prepareProposal( void )
{

}


void PsiMixtureAllocationProposal::printParameterSummary(std::ostream &o, bool name_only) const
{

}


void PsiMixtureAllocationProposal::undoProposal( void )
{
    // the distribution restores the component
}


void PsiMixtureAllocationProposal::setProposalTuningParameter(double tp)
{

}


void PsiMixtureAllocationProposal::tune( double rate )
{

}


void PsiMixtureAllocationProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    matrix = static_cast<StochasticNode<MatrixReal>* >(newN);
}
