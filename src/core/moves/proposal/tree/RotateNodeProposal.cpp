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
    failed( false ),
    redrew( false )
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


/** The named successor of each sampled ancestor: follow its sibling's continuation down to a tip. */
static std::vector<size_t> namedSuccessors( const Tree &tau )
{
    std::vector<size_t> s;
    for (size_t i = 0; i < tau.getNumberOfTips(); i++)
    {
        const TopologyNode &tip = tau.getNode(i);
        if ( tip.isSampledAncestorTip() == false || tip.isRoot() == true ) continue;

        const TopologyNode &sa_node = tip.getParent();
        const TopologyNode *n = &sa_node.getChild( &sa_node.getChild(0) == &tip ? 1 : 0 );
        while ( n->isTip() == false )
        {
            const TopologyNode *next = NULL;
            for (size_t c = 0; c < n->getNumberOfChildren(); c++)
            {
                if ( n->getChild(c).continuesParentSpecies() ) next = &n->getChild(c);
            }
            if ( next == NULL ) break;
            n = next;
        }
        if ( n->isTip() == true ) s.push_back( n->getIndex() );
    }
    return s;
}


/**
 * Pick uniformly among the internal nodes that can flip, the successors of sampled ancestors on an extended tree and,
 * once dnFBDSP holds them, the named successors of a sampled tree. A node moves its continuation to another child, or
 * with symmetric speciation to none; an extended successor flips between anagenesis and symmetric speciation; a
 * named successor redraws its first change.
 */
double RotateNodeProposal::doProposal( void )
{
    RandomNumberGenerator* rng = GLOBAL_RNG;

    Tree& tau = tree->getValue();

    failed = false;
    redrew = false;

    FossilizedBirthDeathSpeciationProcess *fbdsp = dynamic_cast<FossilizedBirthDeathSpeciationProcess *>( &tree->getDistribution() );
    symmetric = ( fbdsp != NULL && fbdsp->hasSymmetricSpeciation() == true );
    bool latents = ( fbdsp != NULL && fbdsp->hasFirstChanges() == true );
    std::vector<size_t> successors;
    if ( latents == true ) successors = namedSuccessors( tau );

    // a sampled ancestor node's continuation is forced onto the tip, but on an extended tree the successor's flag
    // holds whether it arose by anagenesis or a symmetric speciation. A flip leaves both sets unchanged, so their
    // sizes cancel.
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

    size_t total = candidates.size() + mode_nodes.size() + successors.size();
    if ( total == 0 )
    {
        failed = true;
        // reject, so a move that could not act is not counted as accepted
        return RbConstants::Double::neginf;
    }

    size_t k = size_t( rng->uniform01() * total );
    if ( k >= candidates.size() + mode_nodes.size() )
    {
        redrew = true;
        return fbdsp->redrawFirstChange( successors[ k - candidates.size() - mode_nodes.size() ] );
    }

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

    // a flip can name a successor or unname one, which changes the reverse pick's chance
    if ( latents == false ) return 0.0;
    return log( double( total ) ) - log( double( candidates.size() + mode_nodes.size() + namedSuccessors( tau ).size() ) );
}


void RotateNodeProposal::prepareProposal( void )
{
}


void RotateNodeProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
}


void RotateNodeProposal::undoProposal( void )
{
    if ( failed == false && redrew == false )
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
