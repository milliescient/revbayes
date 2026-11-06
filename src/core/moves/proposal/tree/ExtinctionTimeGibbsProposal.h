#ifndef ExtinctionTimeGibbsProposal_H
#define ExtinctionTimeGibbsProposal_H

#include <iosfwd>

#include "Proposal.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore {

    /**
     * On an extended dnFBDSP tree, redraw one fossil tip's extinction time from its conditional given the rest:
     * q~(tau_K -> d) times the range's end, between the present and the younger of its last appearance and its
     * parent. With one psi that is exact and the move always accepts; otherwise it is an independence proposal.
     * The candidates are every extinct fossil tip, a number no move changes, so it cancels; a sampled ancestor
     * among them is rejected.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (June Walker)
     * @since 2026-11-06, version 1.0
     */
    class ExtinctionTimeGibbsProposal : public Proposal {

    public:
        ExtinctionTimeGibbsProposal( StochasticNode<Tree> *n );                                         //!< constructor

        void                                    cleanProposal(void);
        ExtinctionTimeGibbsProposal*            clone(void) const;
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
        StochasticNode<Tree>*                   tree;
        size_t                                  tip_index;
        double                                  stored_age;
        bool                                    failed;
    };

}

#endif
