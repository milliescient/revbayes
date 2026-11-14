#include "Move_RateRangeJoint.h"

#include <cstdint>
#include <string>
#include <vector>

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "MetropolisHastingsMove.h"
#include "Natural.h"
#include "Probability.h"
#include "RateRangeJointProposal.h"
#include "RbException.h"
#include "RealPos.h"
#include "RevObject.h"
#include "RlBoolean.h"
#include "RlMatrixReal.h"
#include "RlString.h"
#include "StochasticNode.h"
#include "TypeSpec.h"
#include "TypedDagNode.h"

using namespace RevLanguage;


Move_RateRangeJoint::Move_RateRangeJoint() : Move()
{

}


Move_RateRangeJoint* Move_RateRangeJoint::clone(void) const
{
    return new Move_RateRangeJoint(*this);
}


void Move_RateRangeJoint::constructInternalObject( void )
{
    delete value;

    double w  = static_cast<const RealPos &>( weight->getRevObject() ).getValue();
    double l  = static_cast<const RealPos &>( lambda->getRevObject() ).getValue();
    bool   t  = static_cast<const RlBoolean &>( tune->getRevObject() ).getValue();
    double tt = static_cast<const Probability &>( tuneTarget->getRevObject() ).getValue();
    bool   rj = static_cast<const RlBoolean &>( rj_switch->getRevObject() ).getValue();
    double sd = static_cast<const RealPos &>( sigma->getRevObject() ).getValue();
    const std::string &sd_side = static_cast<const RlString &>( side->getRevObject() ).getValue();
    if ( sd_side != "extinction" && sd_side != "origination" ) throw RbException("mvRateRangeJoint: side must be \"extinction\" or \"origination\".");
    bool births = ( sd_side == "origination" );

    RevBayesCore::TypedDagNode<RevBayesCore::MatrixReal>* tm = static_cast<const MatrixReal &>( ranges->getRevObject() ).getDagNode();
    RevBayesCore::StochasticNode<RevBayesCore::MatrixReal>* m = dynamic_cast<RevBayesCore::StochasticNode<RevBayesCore::MatrixReal>*>( tm );
    if ( m == NULL )
    {
        throw RbException("mvRateRangeJoint: ranges must be a stochastic range matrix.");
    }

    RevBayesCore::Proposal *p = NULL;
    if ( x->getRevObject().isType( Natural::getClassTypeSpec() ) )
    {
        RevBayesCore::TypedDagNode<std::int64_t>* tx = static_cast<const Natural &>( x->getRevObject() ).getDagNode();
        RevBayesCore::StochasticNode<std::int64_t>* s = dynamic_cast<RevBayesCore::StochasticNode<std::int64_t>*>( tx );
        if ( s == NULL ) throw RbException("mvRateRangeJoint: x must be a stochastic node.");
        if ( rj == true ) throw RbException("mvRateRangeJoint: rjSwitch needs a dnReversibleJumpMixture jump, not a 0/1 indicator.");
        p = new RevBayesCore::RateRangeJointProposal( s, m, births );
        t = false;
    }
    else
    {
        RevBayesCore::TypedDagNode<double>* tx = static_cast<const RealPos &>( x->getRevObject() ).getDagNode();
        RevBayesCore::StochasticNode<double>* r = dynamic_cast<RevBayesCore::StochasticNode<double>*>( tx );
        if ( r == NULL ) throw RbException("mvRateRangeJoint: x must be a stochastic node.");
        p = new RevBayesCore::RateRangeJointProposal( r, m, l, rj, sd, births );
        if ( rj == true ) t = false;
    }
    p->setTargetAcceptanceRate( tt );

    value = new RevBayesCore::MetropolisHastingsMove( p, w, t );
}


const std::string& Move_RateRangeJoint::getClassType(void)
{
    static std::string rev_type = "Move_RateRangeJoint";

    return rev_type;
}


const TypeSpec& Move_RateRangeJoint::getClassTypeSpec(void)
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Move::getClassTypeSpec() ) );

    return rev_type_spec;
}


std::string Move_RateRangeJoint::getMoveName( void ) const
{
    return "RateRangeJoint";
}


const MemberRules& Move_RateRangeJoint::getParameterRules(void) const
{
    static MemberRules move_member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {
        // Natural first: a Natural also converts to RealPos, which would hand the move a converter, not the node
        std::vector<TypeSpec> x_types;
        x_types.push_back( Natural::getClassTypeSpec() );
        x_types.push_back( RealPos::getClassTypeSpec() );
        move_member_rules.push_back( new ArgumentRule( "x", x_types, "The rate, or a 0/1 shift indicator, that the move changes.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "ranges", MatrixReal::getClassTypeSpec(), "The range matrix whose extinction times follow the rate.", ArgumentRule::BY_REFERENCE, ArgumentRule::STOCHASTIC ) );
        move_member_rules.push_back( new ArgumentRule( "lambda", RealPos::getClassTypeSpec(), "The scaling factor of a rate proposal.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(1.0) ) );
        move_member_rules.push_back( new ArgumentRule( "tune", RlBoolean::getClassTypeSpec(), "Should we tune the scaling factor during burnin?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( true ) ) );
        move_member_rules.push_back( new ArgumentRule( "rjSwitch", RlBoolean::getClassTypeSpec(), "Switch a dnReversibleJumpMixture jump between 0 and a rate, instead of scaling it.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( false ) ) );
        move_member_rules.push_back( new ArgumentRule( "side", RlString::getClassTypeSpec(), "Which range ends follow the rate: \"extinction\" times with an extinction rate, \"origination\" times with a speciation rate.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlString("extinction") ) );
        move_member_rules.push_back( new ArgumentRule( "sigma", RealPos::getClassTypeSpec(), "The log-scale sd of a switched-on rate about the rate its cells had.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RealPos(0.5) ) );

        const MemberRules& inheritedRules = Move::getParameterRules();
        move_member_rules.insert( move_member_rules.end(), inheritedRules.begin(), inheritedRules.end() );

        rules_set = true;
    }

    return move_member_rules;
}


const TypeSpec& Move_RateRangeJoint::getTypeSpec( void ) const
{
    static TypeSpec type_spec = getClassTypeSpec();

    return type_spec;
}


void Move_RateRangeJoint::printValue(std::ostream &o) const
{
    o << "Move_RateRangeJoint(" << ( x != NULL ? x->getName() : "?" ) << ")";
}


void Move_RateRangeJoint::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if      ( name == "x" )      x = var;
    else if ( name == "ranges" ) ranges = var;
    else if ( name == "lambda" ) lambda = var;
    else if ( name == "tune" )   tune = var;
    else if ( name == "rjSwitch" ) rj_switch = var;
    else if ( name == "sigma" )  sigma = var;
    else if ( name == "side" )   side = var;
    else                         Move::setConstParameter( name, var );
}
