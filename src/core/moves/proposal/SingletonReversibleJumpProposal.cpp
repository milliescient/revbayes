#include "SingletonReversibleJumpProposal.h"

#include <cstddef>
#include <string>
#include <vector>

#include "AbstractFossilizedBirthDeathRangeProcess.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"

using namespace RevBayesCore;


SingletonReversibleJumpProposal::SingletonReversibleJumpProposal( StochasticNode<MatrixReal> *n ) : Proposal(),
    matrix( n )
{
    addNode( matrix );
    collectEligible();
}


SingletonReversibleJumpProposal::SingletonReversibleJumpProposal( StochasticNode<Tree> *n ) : Proposal(),
    tree( n )
{
    addNode( tree );
    collectEligible();
}


/** An exact age leaves nothing to separate. */
void SingletonReversibleJumpProposal::collectEligible( void )
{
    Distribution &d = ( matrix != NULL ? (Distribution&)matrix->getDistribution()
                                       : (Distribution&)tree->getDistribution() );

    if ( dynamic_cast<AbstractFossilizedBirthDeathRangeProcess *>( &d ) == NULL )
    {
        throw RbException("mvSingletonRJSwitch requires a fossilized birth-death range process.");
    }

    const std::vector<Taxon> &taxa = rangeProcess().getTaxa();
    for (size_t i = 0; i < taxa.size(); ++i)
    {
        size_t count = 0;
        for ( const auto &F : taxa[i].getOccurrences() ) count += F.second;

        if ( count == 1 && taxa[i].getMinAge() < taxa[i].getMaxAge() ) eligible.push_back( i );
    }

    rangeProcess().setHasSingletonMove();
}


/** Resolved on each use, so a swapped node leaves no stale pointer. */
AbstractFossilizedBirthDeathRangeProcess& SingletonReversibleJumpProposal::rangeProcess( void ) const
{
    Distribution &d = ( matrix != NULL ? (Distribution&)matrix->getDistribution()
                                       : (Distribution&)tree->getDistribution() );

    return dynamic_cast<AbstractFossilizedBirthDeathRangeProcess &>( d );
}


void SingletonReversibleJumpProposal::cleanProposal( void )
{
    ; // do nothing
}


SingletonReversibleJumpProposal* SingletonReversibleJumpProposal::clone( void ) const
{
    return new SingletonReversibleJumpProposal( *this );
}


const std::string& SingletonReversibleJumpProposal::getProposalName( void ) const
{
    static std::string name = "SingletonReversibleJump";

    return name;
}


double SingletonReversibleJumpProposal::getProposalTuningParameter( void ) const
{
    return RbConstants::Double::nan;
}


/**
 * Expand tau_1 above tau_K or collapse it onto tau_K. The distribution holds both ages, and a
 * matrix also carries tau_1 in its value.
 *
 * \return The hastings ratio.
 */
double SingletonReversibleJumpProposal::doProposal( void )
{
    if ( eligible.empty() == true || rangeProcess().samplesSingletons() == false ) return RbConstants::Double::neginf;

    stored_index = eligible[ (size_t)( GLOBAL_RNG->uniform01() * eligible.size() ) ];
    stored_first = rangeProcess().getFirstAppearance( stored_index );

    double ln_hastings = rangeProcess().switchSingleton( stored_index );

    if ( matrix != NULL )
    {
        // columns are death, youngest, oldest, divergence
        MatrixReal& v = matrix->getValue();
        v[stored_index][2] = rangeProcess().getFirstAppearance( stored_index );

        matrix->addTouchedElementIndex( stored_index * v.getNumberOfColumns() + 2 );
    }

    return ln_hastings;
}


void SingletonReversibleJumpProposal::prepareProposal( void )
{

}


void SingletonReversibleJumpProposal::printParameterSummary(std::ostream &o, bool name_only) const
{

}


void SingletonReversibleJumpProposal::undoProposal( void )
{
    // the distribution restores the ages and the status
    if ( matrix != NULL )
    {
        MatrixReal& v = matrix->getValue();

        v[stored_index][2] = stored_first;

        matrix->addTouchedElementIndex( stored_index * v.getNumberOfColumns() + 2 );
    }
}


void SingletonReversibleJumpProposal::setProposalTuningParameter(double tp)
{

}


void SingletonReversibleJumpProposal::tune( double rate )
{

}


void SingletonReversibleJumpProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    if ( matrix != NULL ) matrix = static_cast<StochasticNode<MatrixReal>* >(newN);
    else                  tree   = static_cast<StochasticNode<Tree>* >(newN);

    // an mcmc swaps in its own copy of the model, which may predate this move
    rangeProcess().setHasSingletonMove();
}
