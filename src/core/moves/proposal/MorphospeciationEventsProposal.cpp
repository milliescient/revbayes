#include "MorphospeciationEventsProposal.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "PhyloMorphospeciationProcess.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "RbMathLogic.h"

using namespace RevBayesCore;


MorphospeciationEventsProposal::MorphospeciationEventsProposal( StochasticNode<AbstractHomologousDiscreteCharacterData> *n ) : Proposal(),
    characters( n )
{
    addNode( characters );

    if ( dynamic_cast<PhyloMorphospeciationProcess *>( &characters->getDistribution() ) == NULL )
    {
        throw RbException("mvMorphospeciationEvents requires characters drawn from dnPhyloMorphospeciation.");
    }
    if ( process().allowsMultipleChanges() == false )
    {
        throw RbException("mvMorphospeciationEvents needs dnPhyloMorphospeciation with epsilon. Without it, the events on a branch are integrated out.");
    }
}


/** Resolved on each use, so a swapped node leaves no stale pointer. */
PhyloMorphospeciationProcess& MorphospeciationEventsProposal::process( void ) const
{
    return dynamic_cast<PhyloMorphospeciationProcess &>( characters->getDistribution() );
}


void MorphospeciationEventsProposal::cleanProposal( void )
{
}


MorphospeciationEventsProposal* MorphospeciationEventsProposal::clone( void ) const
{
    return new MorphospeciationEventsProposal( *this );
}


const std::string& MorphospeciationEventsProposal::getProposalName( void ) const
{
    static std::string name = "MorphospeciationEvents";
    return name;
}


double MorphospeciationEventsProposal::getProposalTuningParameter( void ) const
{
    return RbConstants::Double::nan;
}


/**
 * Draw the count of hidden events of one kind on one branch from its full conditional.
 *
 * \return The hastings ratio, which cancels the posterior ratio.
 */
double MorphospeciationEventsProposal::doProposal( void )
{
    PhyloMorphospeciationProcess &p = process();
    if ( p.getNumberOfCharacters() == 0 ) return 0.0;

    // which branches can hold events depends on the tree, so the distribution chooses one after preparing it
    return p.resampleEvents();
}


void MorphospeciationEventsProposal::prepareProposal( void )
{
}


void MorphospeciationEventsProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
}


void MorphospeciationEventsProposal::undoProposal( void )
{
    // the distribution restores the index
}


void MorphospeciationEventsProposal::setProposalTuningParameter(double tp)
{
}


void MorphospeciationEventsProposal::tune( double rate )
{
}


void MorphospeciationEventsProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    characters = static_cast<StochasticNode<AbstractHomologousDiscreteCharacterData>* >(newN);
}
