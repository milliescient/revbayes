#ifndef MorphospeciationAllocationProposal_H
#define MorphospeciationAllocationProposal_H

#include "AbstractHomologousDiscreteCharacterData.h"
#include "Proposal.h"
#include "StochasticNode.h"

namespace RevBayesCore {

    class PhyloMorphospeciationProcess;

    /**
     * Gibbs draw of the character that changed at one species' speciation.
     *
     * The index is in no value, so this is the only move that reaches it.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (June Walker)
     * @since 2026-10-19, version 1.0
     *
     */
    class MorphospeciationAllocationProposal : public Proposal {

    public:
        MorphospeciationAllocationProposal( StochasticNode<AbstractHomologousDiscreteCharacterData> *n );

        void                                    cleanProposal(void);
        MorphospeciationAllocationProposal*     clone(void) const;
        double                                  doProposal(void);
        const std::string&                      getProposalName(void) const;
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;
        void                                    setProposalTuningParameter(double tp);
        void                                    tune(double r);
        void                                    undoProposal(void);
        bool                                    allowClamped(void) const override { return true; }  //!< Resamples a latent index of the clamped characters, not their states; see #600

    protected:
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);

    private:
        PhyloMorphospeciationProcess&           process(void) const;

        StochasticNode<AbstractHomologousDiscreteCharacterData>* characters = NULL;
    };

}

#endif
