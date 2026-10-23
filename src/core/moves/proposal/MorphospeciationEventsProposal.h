#ifndef MorphospeciationEventsProposal_H
#define MorphospeciationEventsProposal_H

#include "AbstractHomologousDiscreteCharacterData.h"
#include "Proposal.h"
#include "StochasticNode.h"

namespace RevBayesCore {

    class PhyloMorphospeciationProcess;

    /**
     * Gibbs draw of the count of hidden speciation events on one branch, under
     * dnPhyloMorphospeciation with epsilon.
     *
     * The counts are in no value, so this is the only move that reaches them.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (June Walker)
     * @since 2026-10-22, version 1.0
     *
     */
    class MorphospeciationEventsProposal : public Proposal {

    public:
        MorphospeciationEventsProposal( StochasticNode<AbstractHomologousDiscreteCharacterData> *n );

        void                                    cleanProposal(void);
        MorphospeciationEventsProposal*     clone(void) const;
        double                                  doProposal(void);
        const std::string&                      getProposalName(void) const;
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;
        void                                    setProposalTuningParameter(double tp);
        void                                    tune(double r);
        void                                    undoProposal(void);
        bool                                    allowClamped(void) const override { return true; }  //!< Resamples latent counts of the clamped characters, not their states; see #600

    protected:
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);

    private:
        PhyloMorphospeciationProcess&           process(void) const;

        StochasticNode<AbstractHomologousDiscreteCharacterData>* characters = NULL;
    };

}

#endif
