#ifndef SingletonReversibleJumpProposal_H
#define SingletonReversibleJumpProposal_H

#include "MatrixReal.h"
#include "Proposal.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore {

    class AbstractFossilizedBirthDeathRangeProcess;

    /**
     * Jump one taxon between a single appearance age and two, by expanding tau_1 above tau_K or
     * collapsing it onto tau_K.
     *
     * Under reporting="presence" one reported occurrence may stand for more, so tau_1 and tau_K
     * may differ within its interval. The continuous moves keep their dimension.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (June Walker)
     * @since 2026-10-26, version 1.0
     *
     */
    class SingletonReversibleJumpProposal : public Proposal {

    public:
        SingletonReversibleJumpProposal( StochasticNode<MatrixReal> *n );                           //!< Constructor
        SingletonReversibleJumpProposal( StochasticNode<Tree> *n );                                 //!< Constructor

        // Basic utility functions
        void                                    cleanProposal(void);                                //!< Clean up proposal
        SingletonReversibleJumpProposal*        clone(void) const;                                  //!< Clone object
        double                                  doProposal(void);                                   //!< Perform proposal
        const std::string&                      getProposalName(void) const;                        //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);                              //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const; //!< Print the parameter summary
        void                                    setProposalTuningParameter(double tp);
        void                                    tune(double r);                                     //!< Tune the proposal
        void                                    undoProposal(void);                                 //!< Reject the proposal

    protected:

        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);     //!< Swap the DAG nodes the Proposal is working on

    private:

        AbstractFossilizedBirthDeathRangeProcess& rangeProcess(void) const;                         //!< The distribution, which owns the taxa and the status.

        void                                    collectEligible(void);                              //!< Taxa with one occurrence reported over an interval.

        StochasticNode<MatrixReal>*             matrix = NULL;                                      //!< One of these is the variable, the other stays null.
        StochasticNode<Tree>*                   tree   = NULL;

        std::vector<size_t>                     eligible;                                           //!< Taxa with one occurrence reported over an interval.
        size_t                                  stored_index = 0;                                   //!< The row touched, for the restore.
        double                                  stored_first = 0.0;
    };

}

#endif
