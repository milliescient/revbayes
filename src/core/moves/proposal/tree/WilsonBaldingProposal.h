#ifndef WilsonBaldingProposal_H
#define WilsonBaldingProposal_H

#include <string>
#include <vector>

#include "Proposal.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore {

    /**
     * @brief The Wilson-Balding proposal on a time tree that may hold sampled ancestors.
     *
     * A subtree and its parent node are pruned and regrafted either onto a branch, at an age
     * drawn uniformly along it, or onto an older fossil tip, which becomes a sampled ancestor
     * (Wilson and Balding 1998, Gavryushkina et al. 2014). Species continuation flags travel
     * with the branches. The root is never pruned or regrafted above.

     */
    class WilsonBaldingProposal : public Proposal {

    public:
        WilsonBaldingProposal( StochasticNode<Tree> *n );

        void                                    cleanProposal(void);
        WilsonBaldingProposal*                  clone(void) const;
        double                                  doProposal(void);
        const std::string&                      getProposalName(void) const;
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;
        void                                    undoProposal(void);

    protected:

        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);

    private:

        size_t                                  numberOfPrunable(const Tree &t) const;

        StochasticNode<Tree>*                   tree;
        bool                                    allow_sampled_ancestors;
        bool                                    budding;
        bool                                    symmetric;                                  //!< Whether a regrafted node may also be symmetric, read at each proposal

        // stored state to undo the proposal
        bool                                    failed = false;
        TopologyNode*                           node = NULL;               // i
        TopologyNode*                           parent = NULL;             // iP
        TopologyNode*                           sibling = NULL;            // CiP
        TopologyNode*                           grandparent = NULL;        // PiP
        TopologyNode*                           target = NULL;             // k
        TopologyNode*                           target_parent = NULL;      // k's parent in the pruned tree
        size_t                                  pos_parent = 0;
        size_t                                  pos_sibling = 0;
        size_t                                  pos_target = 0;
        double                                  stored_age = 0.0;
        bool                                    stored_node_continues = false;
        bool                                    stored_sibling_continues = false;
        bool                                    stored_parent_continues = false;
        bool                                    stored_target_continues = false;
        bool                                    stored_sibling_sa = false;
        bool                                    stored_target_sa = false;
    };

}

#endif
