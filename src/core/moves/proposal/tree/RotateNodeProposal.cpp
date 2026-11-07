#include "RotateNodeProposal.h"

#include <cstddef>
#include <ostream>
#include <vector>

#include "FossilizedBirthDeathSpeciationProcess.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "TopologyNode.h"
#include "Tree.h"

using namespace RevBayesCore;


RotateNodeProposal::RotateNodeProposal( StochasticNode<Tree> *n ) : Proposal(),
    tree( n ),
    node_index( 0 ),
    symmetric( false ),
    failed( false )
{
    addNode( tree );
}


void RotateNodeProposal::cleanProposal( void )
{
}


RotateNodeProposal* RotateNodeProposal::clone( void ) const
{
    return new RotateNodeProposal( *this );
}


const std::string& RotateNodeProposal::getProposalName( void ) const
{
    static std::string name = "SpeciationMode";

    return name;
}


double RotateNodeProposal::getProposalTuningParameter( void ) const
{
    return RbConstants::Double::nan;
}


/**
 * Pick uniformly among the internal nodes that can flip and the successors of sampled ancestors. A node moves its
 * continuation to another child, or with symmetric speciation to none; a successor flips between anagenesis and
 * symmetric speciation. A flip leaves both sets unchanged, so their sizes cancel.
 */
double RotateNodeProposal::doProposal( void )
{
    RandomNumberGenerator* rng = GLOBAL_RNG;

    Tree& tau = tree->getValue();

    failed = false;

    FossilizedBirthDeathSpeciationProcess *fbdsp = dynamic_cast<FossilizedBirthDeathSpeciationProcess *>( &tree->getDistribution() );
    symmetric = ( fbdsp != NULL && fbdsp->hasSymmetricSpeciation() == true );

    // a sampled ancestor node's continuation is forced onto the tip, but the successor's flag holds whether it arose
    // by anagenesis or a symmetric speciation
    bool modes = ( fbdsp != NULL && fbdsp->hasFreeSuccessorModes() == true );
    std::vector<size_t> candidates, mode_nodes;
    for (size_t i = 0; i < tau.getNumberOfNodes(); i++)
    {
        const TopologyNode& n = tau.getNode(i);
        if ( modes == true && fbdsp->isModeSuccessor( n ) == true ) mode_nodes.push_back( i );
        if ( n.isTip() == true || n.getNumberOfChildren() < 2 ) continue;
        if ( n.isSampledAncestorTipOrParent() == true ) continue;

        size_t n_cont = 0;
        for (size_t c = 0; c < n.getNumberOfChildren(); c++)
        {
            if ( n.getChild(c).continuesParentSpecies() == true ) ++n_cont;
        }

        if ( n_cont == 1 || ( symmetric == true && n_cont == 0 ) ) candidates.push_back( i );
    }

    size_t total = candidates.size() + mode_nodes.size();
    if ( total == 0 )
    {
        failed = true;
        // reject, so a move that could not act is not counted as accepted
        return RbConstants::Double::neginf;
    }

    size_t k = size_t( rng->uniform01() * total );
    if ( k >= candidates.size() )
    {
        TopologyNode &successor = tau.getNode( mode_nodes[ k - candidates.size() ] );
        stored_children.assign( 1, &successor );
        stored_flags.assign( 1, successor.continuesParentSpecies() );
        successor.setContinuesParentSpecies( stored_flags[0] == false );
        return 0.0;
    }

    node_index = candidates[k];
    TopologyNode& node = tau.getNode( node_index );

    stored_children = node.getChildren();
    stored_flags.clear();
    size_t holder = stored_children.size();
    for (size_t c = 0; c < stored_children.size(); c++)
    {
        stored_flags.push_back( stored_children[c]->continuesParentSpecies() );
        if ( stored_flags[c] == true ) holder = c;
    }

    // the states are one per child plus, with symmetric speciation, none; draw one of the others uniformly,
    // so the reverse draw has the same chance
    size_t states = stored_children.size() + ( symmetric ? 1 : 0 );
    size_t pick = size_t( rng->uniform01() * ( states - 1 ) );
    if ( pick >= holder ) ++pick;

    for (size_t c = 0; c < stored_children.size(); c++)
    {
        stored_children[c]->setContinuesParentSpecies( c == pick );
    }

    return 0.0;
}


void RotateNodeProposal::prepareProposal( void )
{
}


void RotateNodeProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
}


void RotateNodeProposal::undoProposal( void )
{
    if ( failed == false )
    {
        for (size_t c = 0; c < stored_children.size(); c++)
        {
            stored_children[c]->setContinuesParentSpecies( stored_flags[c] );
        }
    }
}


void RotateNodeProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    if ( oldN == tree )
    {
        tree = static_cast<StochasticNode<Tree>* >(newN);
    }
}


void RotateNodeProposal::setProposalTuningParameter(double tp)
{
}


void RotateNodeProposal::tune( double rate )
{
}
