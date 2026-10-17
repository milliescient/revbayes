#ifndef CompleteFossilizedBirthDeathRangeProcess_H
#define CompleteFossilizedBirthDeathRangeProcess_H

#include "FossilizedBirthDeathRangeProcess.h"

namespace RevBayesCore {

    /**
     * @brief The complete fossilized birth-death range process: every lineage is sampled.
     *
     * Lineages are independent, with no unsampled lineages. This is the birth-death-with-rateshifts
     * (BDS) model of Silvestro et al. (2019).
     */
    class CompleteFossilizedBirthDeathRangeProcess : public FossilizedBirthDeathRangeProcess {

    public:
        CompleteFossilizedBirthDeathRangeProcess (const DagNode *speciation,
                                                  const DagNode *extinction,
                                                  const DagNode *psi,
                                                  const TypedDagNode<double>* rho,
                                                  const TypedDagNode<RbVector<double> > *times,
                                                  const std::string &condition,
                                                  const std::vector<Taxon> &taxa,
                                                  bool complete_record,
                                                  const TypedDagNode<double>* origin = NULL,
                                                  double present = 0.0);

        CompleteFossilizedBirthDeathRangeProcess*                       clone(void) const override;

    protected:
        double                                          computeLnProbability(void) override;
        double                                          rangeLnProb(size_t i) override;                 //!< Waiting times and a Poisson record, with no p(t).
        double                                          originLnProb(void) override;                    //!< Applies only with a given origin.
        //!< Sampling is conditioned per range, in rangeLnProb.
        double                                          conditionLnProb(void) const override
                                                        {
                                                            if ( condition == "nTaxa" ) return -lnProbNumTaxa( getNumberOfTaxaAtPresent() );
                                                            return 0.0;
                                                        }
    };
}

#endif
