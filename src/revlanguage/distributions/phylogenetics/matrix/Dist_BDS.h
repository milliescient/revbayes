#ifndef Dist_BDS_H
#define Dist_BDS_H

#include "Dist_FBDRP.h"

namespace RevLanguage {

    /**
     * The RevLanguage wrapper of the range process under complete lineage sampling, dnBDS.
     *
     * Under condition="sampling" this is the BDS model of Silvestro et al. (2019).
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (June Walker)
     * @since 2026-10-24, version 1.0
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
                                                                            bool coexist,
                                                                            bool need_parent) const override;
    };

}

#endif
