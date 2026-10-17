#ifndef ExtinctionReversibleJumpProposal_H
#define ExtinctionReversibleJumpProposal_H

#include "MatrixReal.h"
#include "Proposal.h"
#include "StochasticNode.h"

namespace RevBayesCore {

    class AbstractFossilizedBirthDeathRangeProcess;

    /**
     * Move one taxon's extinction time between the present and a time above it.
     *
     * Under rho < 1 a taxon reported extinct has a point mass at the present, which the
     * continuous element moves cannot reach.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (June Walker)
     * @since 2026-08-03, version 1.0
     *
     */
    class ExtinctionReversibleJumpProposal : public Proposal {

    public:
        ExtinctionReversibleJumpProposal( StochasticNode<MatrixReal> *n );                          //!< Constructor

        // Basic utility functions
        void                                    cleanProposal(void);                                //!< Clean up proposal
        ExtinctionReversibleJumpProposal*       clone(void) const;                                  //!< Clone object
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

        AbstractFossilizedBirthDeathRangeProcess& rangeProcess(void) const;                         //!< The distribution, which owns the present and the taxa.
        void                                    collectEligible(void);                              //!< Taxa the record reports extinct, fixed for the run.

        StochasticNode<MatrixReal>*             matrix = NULL;                                      //!< The range matrix.

        std::vector<size_t>                     eligible;                                           //!< Taxa reported extinct, the only ones with a point mass.
        size_t                                  stored_index = 0;                                   //!< The row touched, for the restore.
        double                                  stored_death = 0.0;
    };

}

#endif
