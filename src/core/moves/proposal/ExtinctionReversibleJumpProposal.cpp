#include "ExtinctionReversibleJumpProposal.h"

#include <cstddef>
#include <string>
#include <vector>

#include "AbstractFossilizedBirthDeathRangeProcess.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"

using namespace RevBayesCore;


ExtinctionReversibleJumpProposal::ExtinctionReversibleJumpProposal( StochasticNode<MatrixReal> *n ) : Proposal(),
    matrix( n )
{
    addNode( matrix );

    if ( dynamic_cast<AbstractFossilizedBirthDeathRangeProcess *>( &matrix->getDistribution() ) == NULL )
    {
        throw RbException("mvExtinctionRJSwitch requires a fossilized birth-death range process.");
    }

    collectEligible();
    rangeProcess().setHasSurvivalMove();
}



/** Only taxa reported extinct can jump. */
void ExtinctionReversibleJumpProposal::collectEligible( void )
{
    const std::vector<Taxon> &taxa = rangeProcess().getTaxa();
    for (size_t i = 0; i < taxa.size(); ++i)
    {
        if ( taxa[i].isExtinct() == true ) eligible.push_back( i );
    }
}


/** Resolved on each use, so a swapped node leaves no stale pointer. */
AbstractFossilizedBirthDeathRangeProcess& ExtinctionReversibleJumpProposal::rangeProcess( void ) const
{
    return dynamic_cast<AbstractFossilizedBirthDeathRangeProcess &>( matrix->getDistribution() );
}


void ExtinctionReversibleJumpProposal::cleanProposal( void )
{
    ; // do nothing
}


ExtinctionReversibleJumpProposal* ExtinctionReversibleJumpProposal::clone( void ) const
{
    return new ExtinctionReversibleJumpProposal( *this );
}


const std::string& ExtinctionReversibleJumpProposal::getProposalName( void ) const
{
    static std::string name = "ExtinctionReversibleJump";

    return name;
}


double ExtinctionReversibleJumpProposal::getProposalTuningParameter( void ) const
{
    return RbConstants::Double::nan;
}


/**
 * Jump one taxon's extinction time between the present and a uniform draw below its youngest
 * appearance.
 *
 * \return The hastings ratio.
 */
double ExtinctionReversibleJumpProposal::doProposal( void )
{
    if ( eligible.empty() == true ) return 0.0;

    RandomNumberGenerator* rng = GLOBAL_RNG;

    stored_index = eligible[ (size_t)( rng->uniform01() * eligible.size() ) ];

    double present = rangeProcess().getPresent();
    // tau_K, the youngest appearance
    double last    = matrix->getValue()[stored_index][1];
    double width   = last - present;

    if ( width <= 0.0 ) return RbConstants::Double::neginf;

    bool survived = rangeProcess().hasSurvived( stored_index );
    double d = survived ? present + rng->uniform01() * width : present;

    rangeProcess().switchSurvived( stored_index );

    MatrixReal& v = matrix->getValue();
    stored_death = v[stored_index][0];
    v[stored_index][0] = d;

    matrix->addTouchedElementIndex( stored_index * v.getNumberOfColumns() );

    return survived ? log( width ) : -log( width );
}


void ExtinctionReversibleJumpProposal::prepareProposal( void )
{

}


void ExtinctionReversibleJumpProposal::printParameterSummary(std::ostream &o, bool name_only) const
{

}


void ExtinctionReversibleJumpProposal::undoProposal( void )
{
    // the distribution restores the status
    MatrixReal& v = matrix->getValue();

    v[stored_index][0] = stored_death;

    matrix->addTouchedElementIndex( stored_index * v.getNumberOfColumns() );
}


void ExtinctionReversibleJumpProposal::setProposalTuningParameter(double tp)
{

}


void ExtinctionReversibleJumpProposal::tune( double rate )
{

}


void ExtinctionReversibleJumpProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    matrix = static_cast<StochasticNode<MatrixReal>* >(newN);
}
