#ifndef CompleteFossilizedBirthDeathRangeProcess_H
#define CompleteFossilizedBirthDeathRangeProcess_H

#include "FossilizedBirthDeathRangeProcess.h"

namespace RevBayesCore {

    /**
     * @brief The fossilized birth-death range process under complete lineage sampling.
     *
     * No lineage goes unsampled, so no range has unsampled descendants to allow for: p is zero
     * and the kernel collapses from q to plain exponential waiting. Nothing is dropped, so
     * nothing collapses, and a range's birth time is the species' own origination rather than a
     * divergence from the nearest sampled range.
     *
     * Nor is anything conditioned here. The reconstructed process has to integrate over the
     * lineages it drops, and that integration is simultaneously its p(t) factors and the
     * collapse of origination into divergence; here the conditioning is declared instead.
     * condition="time" takes the ranges as the complete set, condition="sampling" divides each
     * by its own probability of having left a fossil.
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
        double                                          originLnProb(void) override;                                //!< Applies only with a given origin.
        double                                          conditionLnProb(void) const override;                       //!< The ranges carry the sampling conditioning themselves.
        double                                          lnOriginationSurvival(size_t i, double t) const override;   //!< Nothing collapses, so the origination is the birth time.

        double                                          samplingConditionLnProb(size_t i) const;                    //!< Range i's own conditioning on having left a fossil.
    };
}

#endif
