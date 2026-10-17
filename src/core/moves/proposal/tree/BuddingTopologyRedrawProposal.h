#ifndef BuddingTopologyRedrawProposal_H
#define BuddingTopologyRedrawProposal_H

#include <iosfwd>

#include "Proposal.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore {

    /**
     * Redraw the whole budding topology of a dnFBDSP tree with the ranges fixed.
     *
     * Each lineage buds off one drawn uniformly from those alive at its birth. Every compatible
     * tree has the same density, so the draw is from the exact conditional.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (June Walker)
     * @since 2026-07-19, version 1.0
     *
     */
    class BuddingTopologyRedrawProposal : public Proposal {

    public:
        BuddingTopologyRedrawProposal( StochasticNode<Tree> *n );                                      //!< constructor

        // Basic utility functions
        void                                    cleanProposal(void);                                            //!< Clean up proposal
        BuddingTopologyRedrawProposal*           clone(void) const;                                              //!< Clone object
        double                                  doProposal(void);                                               //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                    //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);                                          //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;   //!< Print the parameter summary
        void                                    setProposalTuningParameter(double tp);
        void                                    tune(double r);                                                 //!< Tune the proposal to achieve a better acceptance/rejection ratio
        void                                    undoProposal(void);                                             //!< Reject the proposal

    protected:

        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);                 //!< Swap the DAG nodes on which the Proposal is working on

    private:

        StochasticNode<Tree>*                   variable;                                                       //!< The variable the Proposal is working on

        Tree                                    stored_tree;
        bool                                    failed;
    };

}

#endif
