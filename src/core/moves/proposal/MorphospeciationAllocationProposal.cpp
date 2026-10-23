#include "MorphospeciationAllocationProposal.h"

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


MorphospeciationAllocationProposal::MorphospeciationAllocationProposal( StochasticNode<AbstractHomologousDiscreteCharacterData> *n ) : Proposal(),
    characters( n )
{
    addNode( characters );

    if ( dynamic_cast<PhyloMorphospeciationProcess *>( &characters->getDistribution() ) == NULL )
    {
        throw RbException("mvMorphospeciationAllocation requires characters drawn from dnPhyloMorphospeciation.");
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

    // with epsilon, the distribution chooses the species after preparing the tree once, and draws its whole set
    if ( p.allowsMultipleChanges() ) return p.resampleChangedSet();

    // which species have a speciation event depends on the topology, so it is read each time
    std::vector<size_t> eligible;
    for (size_t i = 0; i < p.getNumberOfTaxa(); i++)
    {
        if ( p.hasSpeciationEvent(i) ) eligible.push_back( i );
    }
    if ( eligible.empty() ) return 0.0;

    RandomNumberGenerator* rng = GLOBAL_RNG;
    size_t j = eligible[ std::min( eligible.size() - 1, size_t( rng->uniform01() * eligible.size() ) ) ];

    std::vector<double> ln_probs = p.changedCharacterLnProbs( j );
    size_t N = ln_probs.size();
    for (double &x : ln_probs) if ( std::isnan(x) ) x = RbConstants::Double::neginf;

    double top = *std::max_element( ln_probs.begin(), ln_probs.end() );
    if ( RbMath::isFinite(top) == false ) return RbConstants::Double::neginf;

    std::vector<double> probs( N );
    double total = 0.0;
    for (size_t c = 0; c < N; c++)
    {
        probs[c] = std::exp( ln_probs[c] - top );
        total += probs[c];
    }

    double u = rng->uniform01() * total;
    size_t c_new = 0;
    while ( c_new + 1 < N && u >= probs[c_new] )
    {
        u -= probs[c_new];
        c_new++;
    }

    size_t c_old = p.getChangedCharacter( j );
    p.setChangedCharacter( j, c_new );

    return ln_probs[c_old] - ln_probs[c_new];
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
