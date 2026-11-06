#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "Move_ExtinctionTimeGibbs.h"
#include "ExtinctionTimeGibbsProposal.h"
#include "RlBoolean.h"
#include "RealPos.h"
#include "OptionRule.h"
#include "RlString.h"
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

Move_ExtinctionTimeGibbs::Move_ExtinctionTimeGibbs() : Move()
{

}


Move_ExtinctionTimeGibbs* Move_ExtinctionTimeGibbs::clone(void) const
{

    return new Move_ExtinctionTimeGibbs(*this);
}


void Move_ExtinctionTimeGibbs::constructInternalObject( void )
{
    // we free the memory first
    delete value;

    // now allocate a new move
    RevBayesCore::TypedDagNode<RevBayesCore::Tree> *tmp = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::Tree> *t = static_cast<RevBayesCore::StochasticNode<RevBayesCore::Tree> *>( tmp );

    double we = static_cast<const RealPos &>( weight->getRevObject() ).getValue();

    const std::string &tg = static_cast<const RlString &>( target->getRevObject() ).getValue();
    RevBayesCore::ExtinctionTimeGibbsProposal::Target which = ( tg == "node" ? RevBayesCore::ExtinctionTimeGibbsProposal::NODE :
                                                                ( tg == "range" ? RevBayesCore::ExtinctionTimeGibbsProposal::RANGE :
                                                                ( tg == "ancestor" ? RevBayesCore::ExtinctionTimeGibbsProposal::ANCESTOR : RevBayesCore::ExtinctionTimeGibbsProposal::TIP ) ) );
    RevBayesCore::ExtinctionTimeGibbsProposal *p = new RevBayesCore::ExtinctionTimeGibbsProposal( t, which );

    value = new RevBayesCore::MetropolisHastingsMove(p, we, false);
}


/** Get Rev type of object */
const std::string& Move_ExtinctionTimeGibbs::getClassType(void)
{

    static std::string rev_type = "Move_ExtinctionTimeGibbs";

    return rev_type;
}


/** Get class type spec describing type of object */
const TypeSpec& Move_ExtinctionTimeGibbs::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );

    return rev_type_spec;
}


/**
 * Get the Rev name for the constructor function.
 */
std::string Move_ExtinctionTimeGibbs::getMoveName( void ) const
{
    std::string c_name = "ExtinctionTimeGibbs";

    return c_name;
}


/** Return member rules (no members) */
const MemberRules& Move_ExtinctionTimeGibbs::getParameterRules(void) const
{

    static MemberRules move_member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {
        move_member_rules.push_back( new ArgumentRule( "tree", TimeTree::getClassTypeSpec(), "The tree on which this move operates.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        std::vector<std::string> targets = { "tip", "node", "range", "ancestor" };
        move_member_rules.push_back( new OptionRule( "target", new RlString("tip"), targets, "What moves with the extinctions it redraws: one tip alone, a node slid with the tips below it, a taxon's first and last appearance with its own, or a tip turned into a sampled ancestor at its parent and back." ) );

        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );

        rules_set = true;
    }

    return move_member_rules;
}


/** Get type spec */
const TypeSpec& Move_ExtinctionTimeGibbs::getTypeSpec( void ) const
{

    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}


/** Get type spec */
void Move_ExtinctionTimeGibbs::printValue(std::ostream &o) const
{

    o << "Move_ExtinctionTimeGibbs(";
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
void Move_ExtinctionTimeGibbs::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "target" )
    {
        target = var;
    }

    else
    {
        Move::setConstParameter(name, var);
    }
}
