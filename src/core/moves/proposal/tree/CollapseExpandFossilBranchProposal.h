#ifndef CollapseExpandFossilBranchProposal_H
#define CollapseExpandFossilBranchProposal_H

#include <string>

#include "Proposal.h"
#include "StochasticNode.h"
#include "Tree.h"

namespace RevBayesCore {

    class FossilizedBirthDeathSpeciationProcess;
    
    /**
     * @brief Collapse a fossil tip's branch to make it a sampled ancestor, or expand a sampled ancestor onto a branch of its own.
     *
     * The expanded parent age is drawn uniformly between the tip and the grandparent, or the origin at the root.
     *

     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2012-07-12, version 1.0
     *
     */
    class CollapseExpandFossilBranchProposal : public Proposal {
        
    public:
        CollapseExpandFossilBranchProposal( StochasticNode<Tree> *n, TypedDagNode<double>* o);                //!< Constructor
        
        // Basic utility functions
        void                                    cleanProposal(void);                                          //!< Clean up proposal
        CollapseExpandFossilBranchProposal*     clone(void) const;                                            //!< Clone object
        double                                  doProposal(void);                                             //!< Perform proposal
        const std::string&                      getProposalName(void) const;                                  //!< Get the name of the proposal for summary printing
        double                                  getProposalTuningParameter(void) const;
        void                                    prepareProposal(void);                                        //!< Prepare the proposal
        void                                    printParameterSummary(std::ostream &o, bool name_only) const; //!< Print the parameter summary
        void                                    undoProposal(void);                                           //!< Reject the proposal
        
    protected:
        
        void                                    swapNodeInternal(DagNode *oldN, DagNode *newN);               //!< Swap the DAG nodes the Proposal is working on
        void                                    registerMove(void);                                            //!< Tell a dnFBDSP this move attaches
        
        double                                  collapseBranch(TopologyNode &n);
        double                                  expandBranch(TopologyNode &n);
        
    private:
        
        // parameters
        StochasticNode<Tree>*                   tau;                                                          //!< The variable the Proposal is working on
        TypedDagNode<double>*                   origin;                                                       //!< The variable the Proposal is working on
        
        // stored objects to undo proposal
        TopologyNode*                           storedNode = NULL;
        double                                  storedAge = 0.0;
        bool                                    failed = false;
        bool                                    orients = false;                                                //!< On a dnFBDSP tree, the flags pair a budding's orientation with a successor's mode
        bool                                    storedFlag = false;
        const FossilizedBirthDeathSpeciationProcess* fbdsp = NULL;
    
    };
    
}

#endif

