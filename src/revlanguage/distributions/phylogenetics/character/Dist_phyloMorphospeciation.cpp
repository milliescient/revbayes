#include "Dist_phyloMorphospeciation.h"

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DistributionMemberFunction.h"
#include "MethodTable.h"
#include "ModelVector.h"
#include "Natural.h"
#include "RealPos.h"
#include "RevNullObject.h"
#include "RlDistributionMemberFunction.h"
#include "RlTimeTree.h"
#include "StochasticNode.h"

using namespace RevLanguage;


Dist_phyloMorphospeciation::Dist_phyloMorphospeciation() : TypedDistribution< AbstractHomologousDiscreteCharacterData >()
{
}


Dist_phyloMorphospeciation* Dist_phyloMorphospeciation::clone( void ) const
{
    return new Dist_phyloMorphospeciation(*this);
}


RevBayesCore::PhyloMorphospeciationProcess* Dist_phyloMorphospeciation::createDistribution( void ) const
{
    RevBayesCore::TypedDagNode<RevBayesCore::Tree>* t = static_cast<const TimeTree &>( tree->getRevObject() ).getDagNode();

    std::vector<size_t> k;
    if ( num_states->getRevObject() != RevNullObject::getInstance() )
    {
        for (long x : static_cast<const ModelVector<Natural> &>( num_states->getRevObject() ).getValue()) k.push_back( size_t( x ) );
    }

    return new RevBayesCore::PhyloMorphospeciationProcess( t, k );
}


const std::string& Dist_phyloMorphospeciation::getClassType( void )
{
    static std::string rev_type = "Dist_phyloMorphospeciation";
    return rev_type;
}


const TypeSpec& Dist_phyloMorphospeciation::getClassTypeSpec( void )
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution< AbstractHomologousDiscreteCharacterData >::getClassTypeSpec() ) );
    return rev_type_spec;
}


std::string Dist_phyloMorphospeciation::getDistributionFunctionName( void ) const
{
    return "PhyloMorphospeciation";
}


MethodTable Dist_phyloMorphospeciation::getDistributionMethods( void ) const
{
    MethodTable methods = TypedDistribution< AbstractHomologousDiscreteCharacterData >::getDistributionMethods();

    ArgumentRules* changed_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_phyloMorphospeciation, ModelVector<RealPos> >( "getChangedCharacters", variable, changed_rules, true ) );

    return methods;
}


const MemberRules& Dist_phyloMorphospeciation::getParameterRules(void) const
{
    static MemberRules dist_member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {
        dist_member_rules.push_back( new ArgumentRule( "tree", TimeTree::getClassTypeSpec(), "A tree drawn from dnFBDSP, whose speciation events the characters carry.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::STOCHASTIC ) );
        dist_member_rules.push_back( new ArgumentRule( "nStates", ModelVector<Natural>::getClassTypeSpec(), "The number of states of each character. Left out, each character has as many as its largest observed state, and at least two.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );

        rules_set = true;
    }

    return dist_member_rules;
}


const TypeSpec& Dist_phyloMorphospeciation::getTypeSpec( void ) const
{
    static TypeSpec ts = getClassTypeSpec();
    return ts;
}


void Dist_phyloMorphospeciation::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if ( name == "tree" )
    {
        tree = var;
    }
    else if ( name == "nStates" )
    {
        num_states = var;
    }
    else
    {
        TypedDistribution< AbstractHomologousDiscreteCharacterData >::setConstParameter(name, var);
    }
}
