#include "Dist_BDS.h"

#include <string>
#include <vector>

#include "CompleteFossilizedBirthDeathRangeProcess.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "RlDistribution.h"
#include "TypeSpec.h"

using namespace RevLanguage;


Dist_BDS::Dist_BDS() : Dist_FBDRP()
{

}


/** Clone the object */
Dist_BDS* Dist_BDS::clone( void ) const
{
    return new Dist_BDS(*this);
}


/** Get Rev type of object */
const std::string& Dist_BDS::getClassType( void )
{
    static std::string rev_type = "Dist_BDS";

    return rev_type;
}


/** Get class type spec describing type of an object from this class */
const TypeSpec& Dist_BDS::getClassTypeSpec( void )
{
    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( Dist_FBDRP::getClassTypeSpec() ) );

    return rev_type_spec;
}


/** The name the model is cited under, and the name that says what it assumes. */
std::vector<std::string> Dist_BDS::getDistributionFunctionAliases( void ) const
{
    std::vector<std::string> a_names;
    a_names.push_back( "IndependentRanges" );

    return a_names;
}


std::string Dist_BDS::getDistributionFunctionName( void ) const
{
    std::string d_name = "BDS";

    return d_name;
}


const TypeSpec& Dist_BDS::getTypeSpec( void ) const
{
    static TypeSpec ts = getClassTypeSpec();

    return ts;
}


/** The same arguments, handed to the process under complete lineage sampling. */
RevBayesCore::FossilizedBirthDeathRangeProcess* Dist_BDS::makeProcess(RevBayesCore::DagNode *l,
                                                                      RevBayesCore::DagNode *m,
                                                                      RevBayesCore::DagNode *p,
                                                                      RevBayesCore::TypedDagNode<double> *r,
                                                                      RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> > *rt,
                                                                      const std::string &cond,
                                                                      const std::vector<RevBayesCore::Taxon> &t,
                                                                      RevBayesCore::TypedDagNode<double> *oa,
                                                                      double pres,
                                                                      RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> > *mr,
                                                                      RevBayesCore::TypedDagNode<RevBayesCore::Simplex> *mw,
                                                                      bool coexist,
                                                                      bool need_parent) const
{
    return new RevBayesCore::CompleteFossilizedBirthDeathRangeProcess(l, m, p, r, rt, cond, t, true, coexist, need_parent, oa, pres);
}


/** The same rules, except that the coexistence terms default off here and there is no psi mixture. */
const MemberRules& Dist_BDS::getParameterRules(void) const
{
    static MemberRules dist_member_rules;
    static bool rules_set = false;

    if ( !rules_set )
    {
        const MemberRules &parentRules = Dist_FBDRP::getParameterRules();

        for (size_t i = 0; i < parentRules.size(); ++i)
        {
            const std::string &label = parentRules[i].getArgumentLabel();
            if ( label == "coexistence" || label == "requireParent" || label == "psiMixtureRates" || label == "psiMixtureWeights" ) continue;

            dist_member_rules.push_back( parentRules[i].clone() );
        }

        dist_member_rules.push_back( coexistenceRule( false ) );
        dist_member_rules.push_back( requireParentRule() );

        rules_set = true;
    }

    return dist_member_rules;
}
