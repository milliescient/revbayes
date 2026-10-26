#include "MorphospeciationAllocationProposal.h"

#include "PhyloMorphospeciationProcess.h"
#include "RbConstants.h"
#include "RbException.h"

using namespace RevBayesCore;


MorphospeciationAllocationProposal::MorphospeciationAllocationProposal( StochasticNode<AbstractHomologousDiscreteCharacterData> *n ) : Proposal(),
    characters( n )
{
    addNode( characters );

    if ( dynamic_cast<PhyloMorphospeciationProcess *>( &characters->getDistribution() ) == NULL )
    {
        throw RbException("mvMorphospeciationAllocation requires characters drawn from dnPhyloMorphospeciation.");
    }
    if ( process().allowsMultipleChanges() )
    {
        throw RbException("mvMorphospeciationAllocation needs dnPhyloMorphospeciation without epsilon. With it, the characters that changed at each speciation are summed over.");
    }
}


/** Resolved on each use, so a swapped node leaves no stale pointer. */
PhyloMorphospeciationProcess& MorphospeciationAllocationProposal::process( void ) const
{
    return dynamic_cast<PhyloMorphospeciationProcess &>( characters->getDistribution() );
}


void MorphospeciationAllocationProposal::cleanProposal( void )
{
}


MorphospeciationAllocationProposal* MorphospeciationAllocationProposal::clone( void ) const
{
    return new MorphospeciationAllocationProposal( *this );
}


const std::string& MorphospeciationAllocationProposal::getProposalName( void ) const
{
    static std::string name = "MorphospeciationAllocation";
    return name;
}


double MorphospeciationAllocationProposal::getProposalTuningParameter( void ) const
{
    return RbConstants::Double::nan;
}


/**
 * Draw one species' changed character from its full conditional.
 *
 * \return The hastings ratio, which cancels the posterior ratio.
 */
double MorphospeciationAllocationProposal::doProposal( void )
{
    PhyloMorphospeciationProcess &p = process();
    if ( p.getNumberOfCharacters() == 0 ) return 0.0;

    // which speciations there are depends on the tree, so the distribution chooses one after preparing it
    return p.resampleChangedCharacter();
}


void MorphospeciationAllocationProposal::prepareProposal( void )
{
}


void MorphospeciationAllocationProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
}


void MorphospeciationAllocationProposal::undoProposal( void )
{
    // the distribution restores the index
}


void MorphospeciationAllocationProposal::setProposalTuningParameter(double tp)
{
}


void MorphospeciationAllocationProposal::tune( double rate )
{
}


void MorphospeciationAllocationProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    characters = static_cast<StochasticNode<AbstractHomologousDiscreteCharacterData>* >(newN);
}
