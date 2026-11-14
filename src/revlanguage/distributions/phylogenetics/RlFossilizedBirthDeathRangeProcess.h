#ifndef RlFossilizedBirthDeathRangeProcess_H
#define RlFossilizedBirthDeathRangeProcess_H

#include "StochasticNode.h"
#include "TypedDistribution.h"
#include "DistributionMemberFunction.h"
#include "ModelVector.h"
#include "RealPos.h"
#include "RlDistributionMemberFunction.h"
#include "OptionRule.h"
#include "RlDistribution.h"
#include "RevNullObject.h"
#include "RlBoolean.h"
#include "RlString.h"
#include "RlTaxon.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {
    
    /**
     * The RevLanguage wrapper of the FossilizedBirthDeathRangeProcess
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-26, version 1.0
     *
     */
    template <typename rlType>
    class FossilizedBirthDeathRangeProcess : public TypedDistribution<rlType> {
        
    public:
        virtual                                             ~FossilizedBirthDeathRangeProcess(void);                                          //!< Destructor
        
        // Basic utility functions you have to overwrite
        virtual FossilizedBirthDeathRangeProcess<rlType>*   clone(void) const = 0;                                                              //!< Clone the object

        // Basic utility functions you may want to overwrite
        const MemberRules&                                  getParameterRules(void) const;                                                      //!< Get member rules (const)
        MethodTable                                         getDistributionMethods(void) const;                                                 //!< The range times and appearances
        
        // Basic utility functions
        static const std::string&                           getClassType(void);                                                                 //!< Get Rev type
        static const TypeSpec&                              getClassTypeSpec(void);                                                             //!< Get class type spec

    protected:
        FossilizedBirthDeathRangeProcess( void );

        const MemberRules&                                  getCoreParameterRules(void) const;                                              //!< Get member rules (const)
        static void                                         appendParameterRules(MemberRules &rules);                                           //!< Build the member rules
        static ArgumentRule*                                originRule(void);                                                                   //!< The optional origin argument
        static ArgumentRule*                                coexistenceRule(bool on_by_default);                                                //!< The coexistence argument
        static ArgumentRule*                                requireParentRule(void);                                                            //!< The requireParent argument
        RevBayesCore::TypedDagNode<double>*                 createOrigin(void) const;                                                           //!< The given origin node, or NULL if none was given
        std::vector<const RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >*> createRateTimelines(void) const;                    //!< lambdaTimeline, muTimeline and psiTimeline, each NULL if not given

        
        void                                                setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);   //!< Set member variable

        //!< Whether a birth time is a divergence rather than an origination.
        virtual bool                                        hasUnsampledLineages(void) const { return true; }

        // members        
        RevPtr<const RevVariable>                           lambda;                                                                             //!< The speciation rate(s)
        RevPtr<const RevVariable>                           mu;                                                                                 //!< The extinction rate(s)
        RevPtr<const RevVariable>                           psi;                                                                                //!< The fossilization rate(s)
        RevPtr<const RevVariable>                           rho;                                                                                //!< The extant sampling proportion
        RevPtr<const RevVariable>                           timeline;                                                                           //!< The interval times
        RevPtr<const RevVariable>                           lambda_timeline;                                                                    //!< The speciation rate's own interval times
        RevPtr<const RevVariable>                           mu_timeline;                                                                        //!< The extinction rate's own interval times
        RevPtr<const RevVariable>                           psi_timeline;                                                                       //!< The sampling rate's own interval times
        RevPtr<const RevVariable>                           present;                                                                            //!< Where the process stops
        RevPtr<const RevVariable>                           taxa;                                                                               //!< The taxa
        RevPtr<const RevVariable>                           condition;                                                                          //!< The condition of the process
        RevPtr<const RevVariable>                           origin;                                                                             //!< The optional start of the process

    };
    
}


/**
 * Default constructor.
 *
 * The default constructor does nothing except allocating the object.
 */
template <typename rlType>
RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::FossilizedBirthDeathRangeProcess() : TypedDistribution<rlType>()
{

}


/**
 * Default destructor.
 *
 * The default destructor does nothing.
 */
template <typename rlType>
RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::~FossilizedBirthDeathRangeProcess()
{

}


/**
 * Get Rev type of object
 *
 * \return The class' name.
 */
template <typename rlType>
const std::string& RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::getClassType(void)
{

    static std::string rev_type = "FossilizedBirthDeathRange";

    return rev_type;
}


/**
 * Get class type spec describing type of an object from this class (static).
 *
 * \return TypeSpec of this class.
 */
template <typename rlType>
const TypeSpec& RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::getClassTypeSpec(void)
{

    static TypeSpec rev_type_spec = TypeSpec( getClassType(), new TypeSpec( TypedDistribution<rlType>::getClassTypeSpec() ) );

    return rev_type_spec;
}



/**
 * Build the member rules of the range process.
 *
 * \param[out]   rules              The rule list to append to.
 */
template <typename rlType>
void RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::appendParameterRules(MemberRules &rules)
{
    std::vector<TypeSpec> paramTypes;
    paramTypes.push_back( RealPos::getClassTypeSpec() );
    paramTypes.push_back( ModelVector<RealPos>::getClassTypeSpec() );
    rules.push_back( new ArgumentRule( "lambda",  paramTypes, "The speciation rate(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );
    rules.push_back( new ArgumentRule( "mu",      paramTypes, "The extinction rate(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
    rules.push_back( new ArgumentRule( "psi",     paramTypes, "The fossil sampling rate(s).", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );
    rules.push_back( new ArgumentRule( "rho",     Probability::getClassTypeSpec(), "The extant sampling fraction.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(1.0) ) );

    rules.push_back( new ArgumentRule( "timeline",   ModelVector<RealPos>::getClassTypeSpec(), "The rate interval change times of the piecewise constant process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
    rules.push_back( new ArgumentRule( "lambdaTimeline", ModelVector<RealPos>::getClassTypeSpec(), "The rate interval change times of the speciation rate, in place of timeline.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
    rules.push_back( new ArgumentRule( "muTimeline",     ModelVector<RealPos>::getClassTypeSpec(), "The rate interval change times of the extinction rate, in place of timeline.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
    rules.push_back( new ArgumentRule( "psiTimeline",    ModelVector<RealPos>::getClassTypeSpec(), "The rate interval change times of the fossil sampling rate, in place of timeline.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL ) );
    rules.push_back( new ArgumentRule( "present",    RealPos::getClassTypeSpec(), "The time defining the present. Minimum age of the process.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, new RealPos(0.0) ) );

    std::vector<std::string> optionsCondition;
    optionsCondition.push_back( "time" );
    optionsCondition.push_back( "sampling" );
    optionsCondition.push_back( "survival" );
    optionsCondition.push_back( "nTaxa" );
    rules.push_back( new OptionRule( "condition", new RlString("time"), optionsCondition, "The condition of the process." ) );
    rules.push_back( new ArgumentRule( "taxa"  , ModelVector<Taxon>::getClassTypeSpec(), "The taxa with fossil occurrence information.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY ) );

}


template <typename rlType>
ArgumentRule* RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::coexistenceRule(bool on_by_default)
{
    return new ArgumentRule( "coexistence", RlBoolean::getClassTypeSpec(),
                             "Weight each birth by the number of sampled species alive to have produced it?",
                             ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean( on_by_default ) );
}


template <typename rlType>
ArgumentRule* RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::requireParentRule(void)
{
    return new ArgumentRule( "requireParent", RlBoolean::getClassTypeSpec(),
                             "Give zero probability to a birth with no sampled species alive to have produced it? Left out, it follows coexistence. Turn it off when some species were left out of the analysis.",
                             ArgumentRule::BY_VALUE, ArgumentRule::ANY, NULL );
}


template <typename rlType>
ArgumentRule* RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::originRule(void)
{
    return new ArgumentRule( "origin", RealPos::getClassTypeSpec(), "The start of the process. Left out, the process starts at the oldest birth.", ArgumentRule::BY_CONSTANT_REFERENCE, ArgumentRule::ANY, NULL );
}


template <typename rlType>
std::vector<const RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >*> RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::createRateTimelines(void) const
{
    std::vector<const RevBayesCore::TypedDagNode<RevBayesCore::RbVector<double> >*> out;
    for ( const RevPtr<const RevVariable> *v : { &lambda_timeline, &mu_timeline, &psi_timeline } )
    {
        if ( *v == NULL || (*v)->getRevObject() == RevNullObject::getInstance() )
        {
            out.push_back( NULL );
        }
        else
        {
            out.push_back( static_cast<const ModelVector<RealPos> &>( (*v)->getRevObject() ).getDagNode() );
        }
    }
    return out;
}


template <typename rlType>
RevBayesCore::TypedDagNode<double>* RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::createOrigin(void) const
{
    if ( origin == NULL || origin->getRevObject() == RevNullObject::getInstance() )
    {
        return NULL;
    }

    return static_cast<const RealPos &>( origin->getRevObject() ).getDagNode();
}


/**
 * Get the member rules used to create the constructor of this object.
 *
 * \return The member rules.
 */
template <typename rlType>
const MemberRules& RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::getParameterRules(void) const
{

    static MemberRules memberRules;
    static bool rules_set = false;

    if ( !rules_set )
    {
        appendParameterRules( memberRules );

        rules_set = true;
    }

    return memberRules;
}


/**
 * Get the member rules of the range process.
 *
 * \return The member rules.
 */
template <typename rlType>
const MemberRules& RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::getCoreParameterRules(void) const
{

    static MemberRules memberRules;
    static bool rules_set = false;

    if ( !rules_set )
    {
        appendParameterRules( memberRules );

        rules_set = true;
    }

    return memberRules;
}


/**
 * Set a member variable.
 *
 * Sets a member variable with the given name and store the pointer to the variable.
 * The value of the variable might still change but this function needs to be called again if the pointer to
 * the variable changes. The current values will be used to create the distribution object.
 *
 * \param[in]    name     Name of the member variable.
 * \param[in]    var      Pointer to the variable.
 */
template <typename rlType>
void RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var)
{

    if ( name == "lambda" )
    {
        lambda = var;
    }
    else if ( name == "mu" )
    {
        mu = var;
    }
    else if ( name == "psi" )
    {
        psi = var;
    }
    else if ( name == "rho" )
    {
        rho = var;
    }
    else if ( name == "present" )
    {
        present = var;
    }
    else if ( name == "timeline" )
    {
        timeline = var;
    }
    else if ( name == "lambdaTimeline" )
    {
        lambda_timeline = var;
    }
    else if ( name == "muTimeline" )
    {
        mu_timeline = var;
    }
    else if ( name == "psiTimeline" )
    {
        psi_timeline = var;
    }
    else if ( name == "taxa" )
    {
        taxa = var;
    }
    else if ( name == "condition" )
    {
        condition = var;
    }
    else if ( name == "origin" )
    {
        origin = var;
    }
    else
    {
        TypedDistribution<rlType>::setConstParameter(name,var);
    }

}


/** The range times and appearances. */
template <typename rlType>
RevLanguage::MethodTable RevLanguage::FossilizedBirthDeathRangeProcess<rlType>::getDistributionMethods( void ) const
{
    MethodTable methods = TypedDistribution<rlType>::getDistributionMethods();

    ArgumentRules* first_app_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<FossilizedBirthDeathRangeProcess<rlType>, ModelVector<RealPos> >( "getFirstAppearances", this->variable, first_app_arg_rules, true ) );

    ArgumentRules* last_app_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<FossilizedBirthDeathRangeProcess<rlType>, ModelVector<RealPos> >( "getLastAppearances", this->variable, last_app_arg_rules, true ) );

    ArgumentRules* origin_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<FossilizedBirthDeathRangeProcess<rlType>, RealPos >( "getOrigin", this->variable, origin_arg_rules, true ) );

    if ( hasUnsampledLineages() )
    {
        ArgumentRules* divergence_arg_rules = new ArgumentRules();
        methods.addFunction( new DistributionMemberFunction<FossilizedBirthDeathRangeProcess<rlType>, ModelVector<RealPos> >( "getDivergenceTimes", this->variable, divergence_arg_rules, true ) );
    }

    ArgumentRules* extinction_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<FossilizedBirthDeathRangeProcess<rlType>, ModelVector<RealPos> >( "getExtinctionTimes", this->variable, extinction_arg_rules, true ) );

    ArgumentRules* origination_arg_rules = new ArgumentRules();
    origination_arg_rules->push_back( new ArgumentRule( "mean", RlBoolean::getClassTypeSpec(), "Return each origination time's conditional mean instead of a draw?", ArgumentRule::BY_VALUE, ArgumentRule::ANY, new RlBoolean(false) ) );
    methods.addFunction( new DistributionMemberFunction<FossilizedBirthDeathRangeProcess<rlType>, ModelVector<RealPos> >( "getOriginationTimes", this->variable, origination_arg_rules, true ) );

    ArgumentRules* extant_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<FossilizedBirthDeathRangeProcess<rlType>, ModelVector<RealPos> >( "getExtantStatus", this->variable, extant_arg_rules, true ) );

    ArgumentRules* singleton_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<FossilizedBirthDeathRangeProcess<rlType>, ModelVector<RealPos> >( "getSingletonStatus", this->variable, singleton_arg_rules, true ) );

    ArgumentRules* component_arg_rules = new ArgumentRules();
    methods.addFunction( new DistributionMemberFunction<FossilizedBirthDeathRangeProcess<rlType>, ModelVector<RealPos> >( "getPsiMixtureComponents", this->variable, component_arg_rules, true ) );

    return methods;
}

#endif
