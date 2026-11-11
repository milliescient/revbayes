#ifndef PsiMixtureAllocationProposal_H
#define PsiMixtureAllocationProposal_H

#include "MatrixReal.h"
#include "Proposal.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore {

    class AbstractFossilizedBirthDeathRangeProcess;

    /**
     * @brief Gibbs draw of one taxon's psi mixture component, which is in no value and so reached by no other move.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (June Walker)
     * @since 2026-10-17, version 1.0
     *
     */
    class PsiMixtureAllocationProposal : public Proposal {

    public:
        PsiMixtureAllocationProposal( StochasticNode<MatrixReal> *n );                          //!< Constructor
        PsiMixtureAllocationProposal( StochasticNode<Tree> *n );                                //!< Constructor

        // Basic utility functions
        void                                    cleanProposal(void);                                //!< Clean up proposal
        PsiMixtureAllocationProposal*       clone(void) const;                                  //!< Clone object
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

        AbstractFossilizedBirthDeathRangeProcess& rangeProcess(void) const;


        void                                    check(void) const;                                  //!< Refuse a distribution without a mixture.

        StochasticNode<MatrixReal>*             matrix = NULL;                                      //!< One of these is the variable, the other stays null.
        StochasticNode<Tree>*                   tree   = NULL;
    };

}

#endif
