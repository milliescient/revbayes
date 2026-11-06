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
     * leaving every age and the topology alone. With symmetric speciation, a node's three states are
     * either child continuing or neither, and the move proposes one of the other two uniformly.
     *
     * Once dnFBDSP holds the first changes of species on a sampled tree, a fraction of the proposals instead
     * redraw one taxon's from its density given the tree, which only the characters then judge.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (June Walker)
     * @since 2026-07-19, version 1.0
     *
     */
    class RotateNodeProposal : public Proposal {

    public:
        RotateNodeProposal( StochasticNode<Tree> *n, bool sym = false, double sf = 0.5 );                              //!< constructor

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
        bool                                    symmetric;                                                      //!< Whether a node may also be left with no continuing child
        std::vector<TopologyNode*>              stored_children;
        std::vector<bool>                       stored_flags;                                                   //!< Each child's continuation before the proposal
        bool                                    failed;
        double                                  successor_fraction;                                             //!< The share of proposals that redraw a first change, once there are any
        bool                                    redrew;                                                         //!< The distribution undoes a redrawn first change itself
    };

}

#endif
