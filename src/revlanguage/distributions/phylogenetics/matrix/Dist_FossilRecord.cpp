#include "Dist_FossilRecord.h"

#include <string>
#include <vector>

#include "FossilRecordProcess.h"
#include "ArgumentRule.h"
#include "ArgumentRules.h"
#include "ModelVector.h"
#include "OptionRule.h"
#include "RlString.h"
#include "RlMatrixReal.h"
#include "RlTimeTree.h"
#include "RlTaxon.h"
#include "Taxon.h"
#include "TypeSpec.h"

using namespace RevLanguage;


Dist_FossilRecord::Dist_FossilRecord() : TypedDistribution< ModelVector<Taxon> >()
{
}


Dist_FossilRecord* Dist_FossilRecord::clone( void ) const
{
    return new Dist_FossilRecord(*this);
}


RevBayesCore::FossilRecordProcess* Dist_FossilRecord::createDistribution( void ) const
{
    // the range process stochastic node (a dnFBDRP or dnFBDSP node)
    RevBayesCore::DagNode* rn = ranges->getRevObject().getDagNode();

    const std::string &rule = static_cast<const RlString &>( reporting->getRevObject() ).getValue();

    typedef RevBayesCore::AbstractFossilizedBirthDeathRangeProcess Ranges;
    Ranges::Reporting rep = rule == "complete" ? Ranges::COMPLETE : ( rule == "presence" ? Ranges::PRESENCE : Ranges::FIRSTLAST );

    RevBayesCore::FossilRecordProcess* d = new RevBayesCore::FossilRecordProcess( rn, rep );

    return d;
}


const std::string& Dist_FossilRecord::getClassType( void )
{
    static std::string rev_type = "Dist_FossilRecord";
    return rev_type;
}


const TypeSpec& Dist_FossilRecord::getClassTypeSpec( void )
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution< ModelVector<Taxon> >::getClassTypeSpec() ) );
    return rev_type_spec;
}


std::string Dist_FossilRecord::getDistributionFunctionName( void ) const
{
    std::string d_name = "FossilRecord";
    return d_name;
}


const MemberRules& Dist_FossilRecord::getParameterRules(void) const
{
    static MemberRules dist_member_rules;
    static bool rules_set = false;

    if ( rules_set == false )
    {
        std::vector<TypeSpec> rangesTypes;
        rangesTypes.push_back( MatrixReal::getClassTypeSpec() );  // dnFBDRP
        rangesTypes.push_back( TimeTree::getClassTypeSpec() );    // dnFBDSP
        dist_member_rules.push_back( new ArgumentRule( "ranges",  rangesTypes, "The FBD range process (a dnFBDRP or dnFBDSP node) supplying b/d, tau, psi and the timeline.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::STOCHASTIC ) );
        std::vector<std::string> optionsReporting;
        optionsReporting.push_back( "complete" );
        optionsReporting.push_back( "firstlast" );
        optionsReporting.push_back( "presence" );
        dist_member_rules.push_back( new OptionRule( "reporting", new RlString("firstlast"), optionsReporting, "Which sampled occurrences are reported: all of them (complete), the oldest and youngest of each species with the interior count unknown (firstlast), or only that each reported interval holds at least its reported count (presence)." ) );

        rules_set = true;
    }

    return dist_member_rules;
}


const TypeSpec& Dist_FossilRecord::getTypeSpec( void ) const
{
    static TypeSpec ts = getClassTypeSpec();
    return ts;
}


void Dist_FossilRecord::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{
    if ( name == "ranges" )
    {
        ranges = var;
    }
    else if ( name == "reporting" )
    {
        reporting = var;
    }
    else
    {
        TypedDistribution< ModelVector<Taxon> >::setConstParameter(name, var);
    }
}
