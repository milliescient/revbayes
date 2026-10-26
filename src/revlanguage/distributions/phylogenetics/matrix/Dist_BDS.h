#ifndef Dist_BDS_H
#define Dist_BDS_H

#include "Dist_FBDRP.h"

namespace RevLanguage {

    /**
     * The RevLanguage wrapper of the birth-death-with-rateshifts process.
     *
     * The same arguments as dnFBDRP, handed to the process under complete lineage sampling and
     * with the coexistence terms off by default. Under condition="sampling" that is the BDS
     * model of Silvestro et al. (2019), whose likelihood is PyRate's. See
     * CompleteFossilizedBirthDeathRangeProcess for what the assumption costs.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-01-26, version 1.0
     *
     */
    class Dist_BDS : public Dist_FBDRP {

    public:
        Dist_BDS( void );

        // Basic utility functions
        Dist_BDS*                                               clone(void) const;                                                                      //!< Clone the object
        static const std::string&                               getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                                  getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::vector<std::string>                                getDistributionFunctionAliases(void) const;                                             //!< Get the alternative names used for the constructor function in Rev.
        std::string                                             getDistributionFunctionName(void) const;                                                //!< Get the Rev-name for this distribution.
        const TypeSpec&                                         getTypeSpec(void) const;                                                                //!< Get the type spec of the instance
        const MemberRules&                                      getParameterRules(void) const;                                                          //!< Get member rules (const), with the coexistence terms off by default

    protected:
        //!< Nothing is dropped, so nothing collapses and there is no divergence to report.
        bool                                                    hasUnsampledLineages(void) const override { return false; }

        RevBayesCore::FossilizedBirthDeathRangeProcess*         makeProcess(RevBayesCore::DagNode *l,
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
                                                                            bool coexist) const override;
    };

}

#endif
