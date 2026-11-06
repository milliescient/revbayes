#include "WilsonBaldingProposal.h"

#include <cmath>
#include <vector>

#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "TopologyNode.h"
#include "TreeDistributionProperties.h"

using namespace RevBayesCore;


WilsonBaldingProposal::WilsonBaldingProposal( StochasticNode<Tree> *n, bool sym ) : Proposal(),
    tree( n ),
    allow_sampled_ancestors( false ),
    budding( false ),
    symmetric( sym ),
    failed( false ),
    node( NULL ),
    parent( NULL ),
    sibling( NULL ),
    grandparent( NULL ),
    target( NULL ),
    target_parent( NULL )
{
    addNode( tree );

    const TreeDistributionProperties *props = dynamic_cast<const TreeDistributionProperties *>( &tree->getDistribution() );
    allow_sampled_ancestors = ( props != NULL && props->allowsSampledAncestors() == true );
    budding                 = ( props != NULL && props->tracksSpeciesContinuation() == true );
}


void WilsonBaldingProposal::cleanProposal( void )
{
}


WilsonBaldingProposal* WilsonBaldingProposal::clone( void ) const
{
    return new WilsonBaldingProposal( *this );
}


const std::string& WilsonBaldingProposal::getProposalName( void ) const
{
    static std::string name = "WilsonBalding";
    return name;
}


double WilsonBaldingProposal::getProposalTuningParameter( void ) const
{
    return RbConstants::Double::nan;
}


/** Nodes that can be pruned: not the root, not a sampled-ancestor tip, not a child of the root. */
size_t WilsonBaldingProposal::numberOfPrunable( const Tree &t ) const
{
    size_t n = 0;
    for (size_t j = 0; j < t.getNumberOfNodes(); ++j)
    {
        const TopologyNode &k = t.getNode(j);
        if ( k.isRoot() == false && k.isSampledAncestorTip() == false && k.getParent().isRoot() == false ) ++n;
    }
    return n;
}


double WilsonBaldingProposal::doProposal( void )
{
    failed = false;

    RandomNumberGenerator* rng = GLOBAL_RNG;
    Tree &t = tree->getValue();

    std::vector<TopologyNode*> prunable;
    for (size_t j = 0; j < t.getNumberOfNodes(); ++j)
    {
        TopologyNode *k = &t.getNode(j);
        if ( k->isRoot() == false && k->isSampledAncestorTip() == false && k->getParent().isRoot() == false ) prunable.push_back( k );
    }
    if ( prunable.empty() )
    {
        failed = true;
        return RbConstants::Double::neginf;
    }

    node        = prunable[ size_t( rng->uniform01() * prunable.size() ) ];
    parent      = &node->getParent();
    sibling     = &parent->getChild(0) == node ? &parent->getChild(1) : &parent->getChild(0);
    grandparent = &parent->getParent();

    bool from_sampled_ancestor = sibling->isSampledAncestorTip();

    // a range carried past its sampled ancestor into a symmetric node: the regraft onto a fossil makes only an
    // anagenesis, so it could never come back
    if ( from_sampled_ancestor == true && budding == true && node->continuesParentSpecies() == true )
    {
        failed = true;
        return RbConstants::Double::neginf;
    }

    // the parent's states: either child continuing, and with symmetric speciation neither
    double flag_states = ( symmetric == true ? 3.0 : 2.0 );

    // the pruned subtree
    std::vector<bool> in_subtree( t.getNumberOfNodes(), false );
    std::vector<TopologyNode*> stack( 1, node );
    while ( stack.empty() == false )
    {
        TopologyNode *k = stack.back();
        stack.pop_back();
        in_subtree[ k->getIndex() ] = true;
        for (size_t c = 0; c < k->getNumberOfChildren(); ++c) stack.push_back( &k->getChild(c) );
    }

    // every place in the pruned tree the subtree could go: a branch (false) or a fossil tip (true)
    double node_age = node->getAge();
    std::vector<TopologyNode*> targets;
    std::vector<bool> onto_tip;
    for (size_t j = 0; j < t.getNumberOfNodes(); ++j)
    {
        TopologyNode *k = &t.getNode(j);
        if ( k == parent || in_subtree[ k->getIndex() ] == true || k->isRoot() == true ) continue;

        TopologyNode *k_parent = ( k == sibling ? grandparent : &k->getParent() );
        bool k_sa = ( k == sibling ? false : k->isSampledAncestorTip() );
        if ( k_sa == true ) continue;

        if ( k_parent->getAge() > std::max( node_age, k->getAge() ) )
        {
            targets.push_back( k );
            onto_tip.push_back( false );
        }
        if ( allow_sampled_ancestors == true && k->isTip() == true && k->getAge() > node_age )
        {
            targets.push_back( k );
            onto_tip.push_back( true );
        }
    }
    if ( targets.empty() )
    {
        failed = true;
        return RbConstants::Double::neginf;
    }

    size_t pick = size_t( rng->uniform01() * targets.size() );
    target = targets[pick];
    bool to_tip = onto_tip[pick];

    double ln_hastings = std::log( double( prunable.size() ) );

    // the reverse move draws the old age along the old branch, or lands on the old sampled ancestor
    if ( from_sampled_ancestor == false )
    {
        ln_hastings -= std::log( grandparent->getAge() - std::max( node_age, sibling->getAge() ) );
        if ( budding == true ) ln_hastings -= std::log( flag_states );
    }

    stored_age               = parent->getAge();
    stored_node_continues    = node->continuesParentSpecies();
    stored_sibling_continues = sibling->continuesParentSpecies();
    stored_parent_continues  = parent->continuesParentSpecies();
    stored_sibling_sa        = sibling->isSampledAncestorTip();
    stored_target_continues  = target->continuesParentSpecies();
    stored_target_sa         = target->isSampledAncestorTip();

    // prune: the sibling takes the parent's place and its species flag
    pos_parent  = grandparent->removeChild( parent );
    pos_sibling = parent->removeChild( sibling );
    grandparent->addChild( sibling, pos_parent );
    sibling->setParent( grandparent );
    if ( budding == true ) sibling->setContinuesParentSpecies( stored_parent_continues );
    if ( from_sampled_ancestor == true ) sibling->setSampledAncestor( false );

    // regraft above the target, which hands its species flag to the new parent
    target_parent = &target->getParent();
    bool target_continues = target->continuesParentSpecies();
    pos_target = target_parent->removeChild( target );
    target_parent->addChild( parent, pos_target );
    parent->setParent( target_parent );
    parent->addChild( target, pos_sibling );
    target->setParent( parent );
    if ( budding == true ) parent->setContinuesParentSpecies( target_continues );

    double new_age;
    if ( to_tip == true )
    {
        // the fossil becomes a sampled ancestor that carries its species on, and the subtree buds off
        new_age = target->getAge();
        target->setSampledAncestor( true );
        if ( budding == true )
        {
            target->setContinuesParentSpecies( true );
            node->setContinuesParentSpecies( false );
        }
    }
    else
    {
        double lo = std::max( node_age, target->getAge() );
        double hi = target_parent->getAge();
        new_age = lo + rng->uniform01() * (hi - lo);
        ln_hastings += std::log( hi - lo );
        if ( budding == true && symmetric == false )
        {
            bool node_continues = rng->uniform01() < 0.5;
            node->setContinuesParentSpecies( node_continues );
            target->setContinuesParentSpecies( node_continues == false );
            ln_hastings += std::log( flag_states );
        }
        else if ( budding == true )
        {
            size_t state = size_t( rng->uniform01() * flag_states );
            node->setContinuesParentSpecies( state == 0 );
            target->setContinuesParentSpecies( state == 1 );
            ln_hastings += std::log( flag_states );
        }
    }
    parent->setAge( new_age );

    ln_hastings -= std::log( double( numberOfPrunable( t ) ) );

    return ln_hastings;
}


void WilsonBaldingProposal::prepareProposal( void )
{
}


void WilsonBaldingProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
}


/** Undo the regraft, then the prune, in reverse order. */
void WilsonBaldingProposal::undoProposal( void )
{
    if ( failed == true ) return;

    parent->removeChild( target );
    target_parent->removeChild( parent );
    target_parent->addChild( target, pos_target );
    target->setParent( target_parent );
    target->setSampledAncestor( stored_target_sa );
    target->setContinuesParentSpecies( stored_target_continues );

    grandparent->removeChild( sibling );
    parent->addChild( sibling, pos_sibling );
    sibling->setParent( parent );
    grandparent->addChild( parent, pos_parent );
    parent->setParent( grandparent );

    parent->setAge( stored_age );
    sibling->setSampledAncestor( stored_sibling_sa );
    sibling->setContinuesParentSpecies( stored_sibling_continues );
    node->setContinuesParentSpecies( stored_node_continues );
    parent->setContinuesParentSpecies( stored_parent_continues );
}


void WilsonBaldingProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    if ( oldN == tree ) tree = static_cast<StochasticNode<Tree>* >( newN );
}
