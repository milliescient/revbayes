#include "Dist_phyloMorphospeciation.h"

#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "DistributionMemberFunction.h"
#include "MethodTable.h"
#include "ModelVector.h"
#include "RlBoolean.h"
#include "OptionRule.h"
#include "RlString.h"
#include "Natural.h"
#include "Probability.h"
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

    // one number is K for every character, a vector K for each
    std::vector<size_t> k;
    size_t shared = 0;
    const RevObject &ns = num_states->getRevObject();
    if ( ns != RevNullObject::getInstance() )
    {
        if ( ns.isType( Natural::getClassTypeSpec() ) )
        {
            shared = size_t( static_cast<const Natural &>( ns ).getValue() );
        }
        else
        {
            for (long x : static_cast<const ModelVector<Natural> &>( ns ).getValue()) k.push_back( size_t( x ) );
        }
    }

    RevBayesCore::TypedDagNode<double>* eps = NULL;
    if ( epsilon->getRevObject() != RevNullObject::getInstance() )
    {
        eps = static_cast<const Probability &>( epsilon->getRevObject() ).getDagNode();
    }

    RevBayesCore::TypedDagNode<double>* eps_a = NULL;
    if ( epsilon_a->getRevObject() != RevNullObject::getInstance() )
    {
        eps_a = static_cast<const Probability &>( epsilon_a->getRevObject() ).getDagNode();
    }

    RevBayesCore::TypedDagNode<double>* ret = NULL;
    if ( retention->getRevObject() != RevNullObject::getInstance() )
    {
        ret = static_cast<const Probability &>( retention->getRevObject() ).getDagNode();
    }

    const std::string &name = static_cast<const RlString &>( coding->getRevObject() ).getValue();
    RevBayesCore::PhyloMorphospeciationProcess::Coding c = RevBayesCore::PhyloMorphospeciationProcess::ALL;
    if ( name == "variable" )    c = RevBayesCore::PhyloMorphospeciationProcess::VARIABLE;
    if ( name == "informative" ) c = RevBayesCore::PhyloMorphospeciationProcess::INFORMATIVE;

    return new RevBayesCore::PhyloMorphospeciationProcess( t, k, eps, eps_a, c, ret, shared );
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

    ArgumentRules* events_rules = new ArgumentRules();
    events_rules->push_back( new ArgumentRule( "mean", RlBoolean::getClassTypeSpec(), "Above a tip coded for no character the counts are summed out: give their mean there rather than drawing them?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
    methods.addFunction( new DistributionMemberFunction<Dist_phyloMorphospeciation, ModelVector<RealPos> >( "getHiddenEventCounts", variable, events_rules, true ) );

    ArgumentRules* modes_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_phyloMorphospeciation, ModelVector<RealPos> >( "getSpeciationModes", variable, modes_rules, true ) );

    ArgumentRules* probs_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_phyloMorphospeciation, ModelVector<RealPos> >( "getChangeProbabilities", variable, probs_rules, true ) );

    ArgumentRules* means_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<Dist_phyloMorphospeciation, ModelVector<RealPos> >( "getHiddenEventMeans", variable, means_rules, true ) );

    return methods;
}


const MemberRules& Dist_phyloMorphospeciation::getParameterRules(void) const
{
    static MemberRules dist_member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {
        dist_member_rules.push_back( new ArgumentRule( "tree", TimeTree::getClassTypeSpec(), "A tree drawn from dnFBDSP, whose speciation events the characters carry.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::STOCHASTIC ) );
        std::vector<TypeSpec> num_states_types = { ModelVector<Natural>::getClassTypeSpec(), Natural::getClassTypeSpec() };
        dist_member_rules.push_back( new ArgumentRule( "nStates", num_states_types, "The number of states: one per character, or one number shared by every character. Left out, each character has as many as its largest observed state, and at least two.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "epsilon", Probability::getClassTypeSpec(), "The probability that a speciation changes each character, independently, so that it may change several or none. Left out, every speciation changes exactly one character.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        dist_member_rules.push_back( new ArgumentRule( "epsilon_a", Probability::getClassTypeSpec(), "The probability that an anagenetic speciation changes each character. Needs epsilon, which then applies to budding alone. Left out, anagenesis changes with probability epsilon too.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
        std::vector<std::string> coding_options = { "all", "variable", "informative" };
        dist_member_rules.push_back( new OptionRule( "coding", new RlString("all"), coding_options, "Which characters were kept: all of them, those that vary among the taxa scored for them, or those that are parsimony informative. All but \"all\" need epsilon." ) );
        dist_member_rules.push_back( new ArgumentRule( "retention", Probability::getClassTypeSpec(), "Under coding=\"informative\", the probability that a character varying only by autapomorphies was kept too. Left out, none was.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );

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
    else if ( name == "epsilon" )
    {
        epsilon = var;
    }
    else if ( name == "epsilon_a" )
    {
        epsilon_a = var;
    }
    else if ( name == "coding" )
    {
        coding = var;
    }
    else if ( name == "retention" )
    {
        retention = var;
    }
    else
    {
        TypedDistribution< AbstractHomologousDiscreteCharacterData >::setConstParameter(name, var);
    }
}
