#ifndef CompleteFossilizedBirthDeathRangeProcess_H
#define CompleteFossilizedBirthDeathRangeProcess_H

#include "FossilizedBirthDeathRangeProcess.h"

namespace RevBayesCore {

    /**
     * @brief The fossilized birth-death range process under complete lineage sampling.
     *
     * p is zero, and a range's birth time is the species' own origination.
     */
    class CompleteFossilizedBirthDeathRangeProcess : public FossilizedBirthDeathRangeProcess {

    public:
        CompleteFossilizedBirthDeathRangeProcess(const DagNode *speciation,
                                                 const DagNode *extinction,
                                                 const DagNode *psi,
                                                 const TypedDagNode<double>* rho,
                                                 const TypedDagNode<RbVector<double> > *times,
                                                 const std::string &condition,
                                                 const std::vector<Taxon> &taxa,
                                                 bool complete_record,
                                                 bool coexistence,
                                                 bool require_parent,
                                                 const TypedDagNode<double>* origin = NULL,
                                                 double present = 0.0);

        CompleteFossilizedBirthDeathRangeProcess*       clone(void) const override;

    protected:
        double                                          rangeLnProb(size_t i) override;                             //!< Exponential waiting and a Poisson record, with no p(t).
        bool                                            chargesBirth(size_t i) const override;                      //!< The oldest birth pays too, unless the origin is given.
        double                                          conditionLnProb(void) const override;                       //!< The ranges carry the sampling conditioning themselves.
        double                                          lnOriginationSurvival(size_t i, double t) const override;   //!< Nothing collapses, so the origination is the birth time.

        double                                          samplingConditionLnProb(size_t i) const;                    //!< Range i's own conditioning on having left a fossil.
    };
}

#endif
