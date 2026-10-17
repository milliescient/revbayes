#ifndef StratigraphicRangeProposal_H
#define StratigraphicRangeProposal_H

#include "Proposal.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore {

    class FossilizedBirthDeathSpeciationProcess;

    /**
     * Resample one taxon's first and last appearances within the bins its occurrences reported.
     *
     * For dnFBDSP, whose tree does not carry the appearances.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (June Walker)
     * @since 2026-08-01, version 1.0
     *
     */
    class StratigraphicRangeProposal : public Proposal {

    public:
        StratigraphicRangeProposal( StochasticNode<Tree> *n );                                                         //!< Constructor

        // Basic utility functions
        bool                                    allowClamped() const override { return true; }                         //!< Valid on a clamped tree, since it samples the appearances. See #600.
        void                                    cleanProposal(void) override;                                          //!< Clean up proposal
        StratigraphicRangeProposal*             clone(void) const override;                                            //!< Clone object
        double                                  doProposal(void) override;                                             //!< Perform proposal
        const std::string&                      getProposalName(void) const override;                                  //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const override;
        void                                    prepareProposal(void) override;                                        //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const override; //!< Print the parameter summary
        void                                    undoProposal(void) override;                                           //!< Reject the proposal

    protected:

        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN) override;               //!< Swap the DAG nodes the Proposal is working on

    private:

        FossilizedBirthDeathSpeciationProcess&  rangeProcess(void) const;                                              //!< The distribution, which owns the appearances.

        StochasticNode<Tree>*                   variable;                                                              //!< The variable the Proposal is working on
    };

}

#endif
