#ifndef RotateNodeProposal_H
#define RotateNodeProposal_H

#include <cstddef>
#include <iosfwd>
#include <vector>

#include "Proposal.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore {

    class TopologyNode;

    /**
     * Move the continuation of the parent's species to another child of a random internal node,
     * leaving every age and the topology alone.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (June Walker)
     * @since 2026-07-19, version 1.0
     *
     */
    class RotateNodeProposal : public Proposal {

    public:
        RotateNodeProposal( StochasticNode<Tree> *n );                                                //!< constructor

        // Basic utility functions
        void                                    cleanProposal(void);                                            //!< Clean up proposal
        RotateNodeProposal*                     clone(void) const;                                              //!< Clone object
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


        StochasticNode<Tree>*                   tree;                                                           //!< The variable the Proposal is working on

        size_t                                  node_index;
        std::vector<TopologyNode*>              stored_children;
        TopologyNode*                           stored_continuer = NULL;                                 //!< The child that held the continuation before the swap
        TopologyNode*                           moved_continuer = NULL;                                  //!< The child it moved to
        bool                                    failed;
    };

}

#endif
