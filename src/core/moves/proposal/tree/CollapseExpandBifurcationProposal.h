#ifndef CollapseExpandBifurcationProposal_H
#define CollapseExpandBifurcationProposal_H

#include <cstddef>
#include <iosfwd>
#include <vector>

#include "Proposal.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore {

    class TopologyNode;

    /**
     * @brief On an extended tree, swap a sampled species ending at a symmetric speciation for two buddings and a tip.
     *
     *     collapsed:  A(t_A)[ i(t_A, sampled ancestor) , B(t_A)[ D1 , D2 ] ]
     *     expanded:   A(t_A)[ D1 , B(t_B)[ D2 , i(t_s) ] ]       t_A > t_B > t_s
     *
     * The expand picks D2 of B's two children and draws t_B uniformly above D2 and t_s uniformly above the present.
     * The collapse is determined. Every fossil tip is a candidate, so the choice of i cancels.

     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (June Walker)
     * @since 2026-11-05, version 1.0
     */
    class CollapseExpandBifurcationProposal : public Proposal {

    public:
        CollapseExpandBifurcationProposal( StochasticNode<Tree> *n );                                     //!< constructor

        void                                    cleanProposal(void);
        CollapseExpandBifurcationProposal*        clone(void) const;
        double                                  doProposal(void);
        const std::string&                      getProposalName(void) const;
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);
        void                                    printParameterSummary(std::ostream &o, bool name_only) const;
        void                                    setProposalTuningParameter(double tp);
        void                                    tune(double r);
        void                                    undoProposal(void);

    protected:
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);

    private:
        double                                  collapse(TopologyNode &i);
        double                                  expand(TopologyNode &i);
        void                                    rewire(TopologyNode &a, const std::vector<TopologyNode*> &a_kids, TopologyNode &b, const std::vector<TopologyNode*> &b_kids);

        StochasticNode<Tree>*                   tree;
        double                                  present;
        bool                                    failed;

        TopologyNode*                           stored_tip = NULL;
        TopologyNode*                           stored_a = NULL;
        TopologyNode*                           stored_b = NULL;
        std::vector<TopologyNode*>              stored_a_kids;
        std::vector<TopologyNode*>              stored_b_kids;
        double                                  stored_b_age;
        double                                  stored_tip_age;
        bool                                    stored_sa;
    };

}

#endif
