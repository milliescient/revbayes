#include <cmath>
#include <iostream>
#include <vector>

#include "CollapseExpandBifurcationProposal.h"
#include "AbstractFossilizedBirthDeathRangeProcess.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "TopologyNode.h"
#include "Tree.h"
#include "TreeDistributionProperties.h"

using namespace RevBayesCore;


CollapseExpandBifurcationProposal::CollapseExpandBifurcationProposal( StochasticNode<Tree> *n ) : Proposal(),
    tree( n ),
    present( 0.0 ),
    failed( false ),
    stored_b_age( 0.0 ),
    stored_tip_age( 0.0 ),
    stored_sa( false )
{
    addNode( tree );

    const TreeDistributionProperties *props = dynamic_cast<const TreeDistributionProperties *>( &tree->getDistribution() );
    if ( props == NULL || props->isExtended() == false || props->tracksSpeciesContinuation() == false )
    {
        throw RbException("mvCollapseExpandBifurcation needs an extended tree that tracks which child continues each species, as dnFBDSP does.");
    }

    const AbstractFossilizedBirthDeathRangeProcess *range = dynamic_cast<const AbstractFossilizedBirthDeathRangeProcess *>( &unconstrained( tree->getDistribution() ) );
    if ( range != NULL ) present = range->getPresent();
}


void CollapseExpandBifurcationProposal::cleanProposal( void )
{
}


CollapseExpandBifurcationProposal* CollapseExpandBifurcationProposal::clone( void ) const
{
    return new CollapseExpandBifurcationProposal( *this );
}


const std::string& CollapseExpandBifurcationProposal::getProposalName( void ) const
{
    static std::string name = "CollapseExpandBifurcation";
    return name;
}


double CollapseExpandBifurcationProposal::getProposalTuningParameter( void ) const
{
    return RbConstants::Double::nan;
}


/** The single continuing child of a node with no sampled ancestor, or NULL. */
static TopologyNode* onlyContinuing( TopologyNode &n )
{
    if ( n.isTip() == true || n.getNumberOfChildren() != 2 ) return NULL;

    TopologyNode *found = NULL;
    for (size_t c = 0; c < 2; c++)
    {
        TopologyNode &k = n.getChild(c);
        if ( k.isSampledAncestorTip() == true ) return NULL;
        if ( k.continuesParentLineage() == true )
        {
            if ( found != NULL ) return NULL;
            found = &k;
        }
    }
    return found;
}


static TopologyNode& otherChild( TopologyNode &n, const TopologyNode &k )
{
    return ( &n.getChild(0) == &k ? n.getChild(1) : n.getChild(0) );
}


double CollapseExpandBifurcationProposal::doProposal( void )
{
    RandomNumberGenerator* rng = GLOBAL_RNG;
    Tree &t = tree->getValue();

    failed = true;

    std::vector<TopologyNode*> fossils;
    for (size_t i = 0; i < t.getNumberOfTips(); ++i)
    {
        TopologyNode &n = t.getNode(i);
        if ( n.isFossil() == true && n.isRoot() == false ) fossils.push_back( &n );
    }
    if ( fossils.empty() == true ) return RbConstants::Double::neginf;

    TopologyNode &i = *fossils[ size_t( rng->uniform01() * fossils.size() ) ];

    return ( i.isSampledAncestorTip() == true ? expand( i ) : collapse( i ) );
}


/** Replace the children of a and of b, in the given order. */
void CollapseExpandBifurcationProposal::rewire( TopologyNode &a, const std::vector<TopologyNode*> &a_kids, TopologyNode &b, const std::vector<TopologyNode*> &b_kids )
{
    std::vector<TopologyNode*> old_a = a.getChildren();
    std::vector<TopologyNode*> old_b = b.getChildren();
    for (TopologyNode *k : old_a) a.removeChild( k );
    for (TopologyNode *k : old_b) b.removeChild( k );
    for (TopologyNode *k : a_kids) { a.addChild( k ); k->setParent( &a ); }
    for (TopologyNode *k : b_kids) { b.addChild( k ); k->setParent( &b ); }
}


/** i at the bottom of two buddings, each continuing the range down to i's tip: put both ages at the upper. */
double CollapseExpandBifurcationProposal::collapse( TopologyNode &i )
{
    if ( i.continuesParentLineage() == false ) return RbConstants::Double::neginf;

    TopologyNode &b = i.getParent();
    if ( b.isRoot() == true || onlyContinuing( b ) != &i || b.continuesParentLineage() == false ) return RbConstants::Double::neginf;

    TopologyNode &a = b.getParent();
    if ( onlyContinuing( a ) != &b ) return RbConstants::Double::neginf;

    TopologyNode &d1 = otherChild( a, b );
    TopologyNode &d2 = otherChild( b, i );

    double t_a = a.getAge();
    double t_b = b.getAge();

    stored_tip = &i; stored_a = &a; stored_b = &b;
    stored_a_kids = a.getChildren(); stored_b_kids = b.getChildren();
    stored_b_age = t_b; stored_tip_age = i.getAge(); stored_sa = false;

    // A[ i , B[ D1 , D2 ] ], with i the sampled ancestor and B the symmetric node, both at t_A
    rewire( a, { &i, &b }, b, { &d1, &d2 } );
    b.setAge( t_a );
    i.setAge( t_a );
    i.setSampledAncestor( true );

    failed = false;

    // the reverse expand picks D2 of two and draws t_B and t_s uniformly
    return -( log( 2.0 ) + log( t_a - d2.getAge() ) + log( t_b - present ) );
}


/** i a sampled ancestor whose range carries on into a symmetric node at its age: split it into two buddings and a tip. */
double CollapseExpandBifurcationProposal::expand( TopologyNode &i )
{
    RandomNumberGenerator* rng = GLOBAL_RNG;

    TopologyNode &a = i.getParent();
    TopologyNode &b = otherChild( a, i );
    if ( b.isTip() == true || b.continuesParentLineage() == false || b.hasZeroLengthBranch() == false ) return RbConstants::Double::neginf;
    if ( b.getNumberOfChildren() != 2 ) return RbConstants::Double::neginf;
    for (size_t c = 0; c < 2; c++)
    {
        if ( b.getChild(c).continuesParentLineage() == true || b.getChild(c).isSampledAncestorTip() == true ) return RbConstants::Double::neginf;
    }

    size_t pick = size_t( rng->uniform01() * 2 );
    TopologyNode &d2 = b.getChild( pick );
    TopologyNode &d1 = b.getChild( 1 - pick );

    double t_a = a.getAge();
    double lo_b = d2.getAge();
    if ( t_a <= lo_b || t_a <= present ) return RbConstants::Double::neginf;

    double t_b = lo_b + rng->uniform01() * ( t_a - lo_b );
    if ( t_b <= present ) return RbConstants::Double::neginf;
    double t_s = present + rng->uniform01() * ( t_b - present );

    stored_tip = &i; stored_a = &a; stored_b = &b;
    stored_a_kids = a.getChildren(); stored_b_kids = b.getChildren();
    stored_b_age = b.getAge(); stored_tip_age = i.getAge(); stored_sa = true;

    // A[ D1 , B[ D2 , i ] ]: i buds D1 at t_A, D2 at t_B, and ends at t_s
    rewire( a, { &d1, &b }, b, { &d2, &i } );
    i.setSampledAncestor( false );
    b.setAge( t_b );
    i.setAge( t_s );
    d1.setContinuesParentLineage( false );
    d2.setContinuesParentLineage( false );
    b.setContinuesParentLineage( true );
    i.setContinuesParentLineage( true );

    failed = false;

    return log( 2.0 ) + log( t_a - lo_b ) + log( t_b - present );
}


void CollapseExpandBifurcationProposal::prepareProposal( void )
{
}


void CollapseExpandBifurcationProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
}


void CollapseExpandBifurcationProposal::undoProposal( void )
{
    if ( failed == true ) return;

    TopologyNode &a = *stored_a;
    TopologyNode &b = *stored_b;
    TopologyNode &i = *stored_tip;

    i.setSampledAncestor( false );

    rewire( a, stored_a_kids, b, stored_b_kids );
    b.setAge( stored_b_age );
    i.setAge( stored_tip_age );
    i.setSampledAncestor( stored_sa );
}


void CollapseExpandBifurcationProposal::setProposalTuningParameter(double tp)
{
}


void CollapseExpandBifurcationProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    if ( oldN == tree ) tree = static_cast<StochasticNode<Tree>* >( newN );
}


void CollapseExpandBifurcationProposal::tune( double rate )
{
}
