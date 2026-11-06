#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "Move_RotateNode.h"
#include "RotateNodeProposal.h"
#include "RlBoolean.h"
#include "Probability.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlTimeTree.h"
#include "TypedDagNode.h"
#include "TypeSpec.h"
#include "Move.h"
#include "RevPtr.h"
#include "RevVariable.h"
#include "RlMove.h"
#include "StochasticNode.h"
#include "Tree.h"


using namespace RevLanguage;

Move_RotateNode::Move_RotateNode() : Move()
{

}


Move_RotateNode* Move_RotateNode::clone(void) const
{

    return new Move_RotateNode(*this);
}


void Move_RotateNode::constructInternalObject( void )
{
    // we free the memory first
    delete value;

    // now allocate a new move
    RevBayesCore::TypedDagNode<RevBayesCore::Tree> *tmp = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *t = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );

    double we = static_cast<const RealPos &>( weight->getRevObject() ).getValue();

    bool sym = static_cast<const RlBoolean &>( symmetric->getRevObject() ).getValue();

    double sf = static_cast<const Probability &>( successor_fraction->getRevObject() ).getValue();

    RevBayesCore::RotateNodeProposal *p = new RevBayesCore::RotateNodeProposal( t, sym, sf );

    value = new RevBayesCore::MetropolisHastingsMove(p, we, false);
}


/** Get Rev type of object */
const std::string& Move_RotateNode::getClassType(void)
{

    static std::string rev_type = "Move_RotateNode";

    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Move_RotateNode::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 */
std::string Move_RotateNode::getMoveName( void ) const
{
    std::string c_name = "SpeciationMode";

    return c_name;
}


std::vector<std::string> Move_RotateNode::getMoveAliases( void ) const
{
    return std::vector<std::string>( 1, "RotateNode" );
}


/** Return member rules (no members) */
const MemberRules& Move_RotateNode::getParameterRules(void) const
{

    static MemberRules move_member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {
        move_member_rules.push_back( new ArgumentRule( "tree", TimeTree::getClassTypeSpec(), "The tree on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "symmetric", RlBoolean::getClassTypeSpec(), "Whether a node may also be left with no child continuing its species, a symmetric speciation. Use with beta > 0.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false ) ) );
        move_member_rules.push_back( new ArgumentRule( "successorFraction", Probability::getClassTypeSpec(), "The fraction of proposals that redraw a sampled ancestor's successor's first change of species, on a sampled tree whose event accessors exist.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new Probability( 0.5 ) ) );

        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );

        rules_set = true;
    }

    return move_member_rules;
}


/** Get type spec */
const TypeSpec& Move_RotateNode::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}


/** Get type spec */
void Move_RotateNode::printValue(std::ostream &o) const
{

    o << "Move_RotateNode(";
    if (tree != NULL)
    {
        o << tree->getName();
    }
    else
    {
        o << "?";
    }
    o << ")";
}


/** Set a member variable */
void Move_RotateNode::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "symmetric" )
    {
        symmetric = var;
    }
    else if ( name == "successorFraction" )
    {
        successor_fraction = var;
    }
    else
    {
        Move::setConstParameter(name, var);
    }
}
