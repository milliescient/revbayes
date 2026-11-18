#include "Dist_FossilOccurrences.h"

#include <string>
#include <vector>

#include "FossilOccurrenceDistribution.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelVector.h"
#include "OptionRule.h"
#include "RealPos.h"
#include "RevNullObject.h"
#include "RlString.h"
#include "RlMatrixReal.h"
#include "RlTaxon.h"
#include "Taxon.h"
#include "TypeSpec.h"

using namespace RevLanguage;


Dist_FossilOccurrences::Dist_FossilOccurrences() : TypedDistribution< ModelVector<Taxon> >()
{
}


Dist_FossilOccurrences* Dist_FossilOccurrences::clone( void ) const
{
    return new Dist_FossilOccurrences(*this);
}


RevBayesCore::FossilOccurrenceDistribution* Dist_FossilOccurrences::createDistribution( void ) const
{
    // the range process node
    RevBayesCore::DagNode* rn = ranges->getRevObject().getDagNode();

    const std::string &rule = static_cast<const RlString &>( reporting->getRevObject() ).getValue();

    typedef RevBayesCore::AbstractFossilizedBirthDeathRangeProcess Ranges;
    Ranges::Reporting rep = rule == "complete" ? Ranges::COMPLETE : ( rule == "presence" ? Ranges::PRESENCE : Ranges::FIRSTLAST );

    std::vector<double> b;
    if ( bins->getRevObject() != RevNullObject::getInstance() )
    {
        const RevBayesCore::RbVector<double> &v = static_cast<const ModelVector<RealPos> &>( bins->getRevObject() ).getValue();
        for ( size_t i = 0; i < v.size(); ++i ) b.push_back( v[i] );
    }

    RevBayesCore::FossilOccurrenceDistribution* d = new RevBayesCore::FossilOccurrenceDistribution( rn, rep, b );

    return d;
}


const std::string& Dist_FossilOccurrences::getClassType( void )
{
    static std::string rev_type = "Dist_FossilOccurrences";
    return rev_type;
}


const TypeSpec& Dist_FossilOccurrences::getClassTypeSpec( void )
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution< ModelVector<Taxon> >::getClassTypeSpec() ) );
    return rev_type_spec;
}


std::string Dist_FossilOccurrences::getDistributionFunctionName( void ) const
{
    std::string d_name = "FossilOccurrences";
    return d_name;
}


const MemberRules& Dist_FossilOccurrences::getParameterRules(void) const
{
    static MemberRules dist_member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {
        std::vector<TypeSpec> rangesTypes;
        rangesTypes.push_back( MatrixReal::getClassTypeSpec() );  // dnFBDRP
        dist_member_rules.push_back( new ArgumentRule( "ranges",  rangesTypes, "The range process (a dnFBDRP variable) whose species produced the record.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::STOCHASTIC ) );
        std::vector<std::string> optionsReporting;
        optionsReporting.push_back( "complete" );
        optionsReporting.push_back( "firstlast" );
        optionsReporting.push_back( "presence" );
        dist_member_rules.push_back( new OptionRule( "reporting", new RlString("firstlast"), optionsReporting, "Which sampled occurrences are reported: all of them (complete), the oldest and youngest (firstlast), or at least one in each interval (presence)." ) );

        dist_member_rules.push_back( new ArgumentRule( "bins", ModelVector<RealPos>::getClassTypeSpec(), "Time interval boundaries for simulated occurrence age uncertainty. If omitted, simulated occurrences have exact ages.", ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL ) );

        rules_set = true;
    }

    return dist_member_rules;
}


const TypeSpec& Dist_FossilOccurrences::getTypeSpec( void ) const
{
    static TypeSpec ts = getClassTypeSpec();
    return ts;
}


void Dist_FossilOccurrences::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if ( name == "ranges" )
    {
        ranges = var;
    }
    else if ( name == "reporting" )
    {
        reporting = var;
    }
    else if ( name == "bins" )
    {
        bins = var;
    }
    else
    {
        TypedDistribution< ModelVector<Taxon> >::setConstParameter(name, var);
    }
}
