#ifndef PhyloMorphospeciationProcess_H
#define PhyloMorphospeciationProcess_H

#include <vector>

#include "AbstractHomologousDiscreteCharacterData.h"
#include "MemberObject.h"
#include "RbVector.h"
#include "Tree.h"
#include "TypedDagNode.h"
#include "TypedDistribution.h"

namespace RevBayesCore {

    class FossilizedBirthDeathSpeciationProcess;
    class TopologyNode;

    /**
     * @brief Morphological characters whose changes are the speciation events of a dnFBDSP tree.
     *
     * Every speciation changes the new species at one character, chosen uniformly, to a state drawn
     * from that character's Mk kernel. Conditional on the tree, a character evolves along a species'
     * stem as a Markov chain with generator J - I over the length (1/n) int (lambda p + lambda_a) dt,
     * does not change over the species' sampled range, and starts from the kernel's stationary
     * distribution at the root. Which character changed at each speciation is a latent index per
     * species, sampled by mvMorphospeciationAllocation.
     */
    class PhyloMorphospeciationProcess : public TypedDistribution<AbstractHomologousDiscreteCharacterData>, public MemberObject< RbVector<double> > {

    public:
        PhyloMorphospeciationProcess(const TypedDagNode<Tree> *tree, const std::vector<size_t> &num_states);
        virtual ~PhyloMorphospeciationProcess() {}

        PhyloMorphospeciationProcess*   clone(void) const override;
        double                          computeLnProbability(void) override;
        void                            redrawValue(void) override;
        void                            setValue(AbstractHomologousDiscreteCharacterData *v, bool force=false) override;

        void                            executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const override;
        std::string                     getHiddenStateString(void) const override;              //!< The latent variables, for checkpointing clamped characters.
        void                            setHiddenStateFromString(const std::string &s) override;

        size_t                          getNumberOfCharacters(void) const { return num_states.size(); }
        size_t                          getNumberOfTaxa(void) const { return changed.size(); }
        bool                            hasSpeciationEvent(size_t taxon) const;                 //!< False for the species that begins at the origin.
        size_t                          getChangedCharacter(size_t taxon) const { return changed[taxon]; }
        void                            setChangedCharacter(size_t taxon, size_t c);            //!< Arms the undo. Only mvMorphospeciationAllocation may.
        std::vector<double>             changedCharacterLnProbs(size_t taxon);                  //!< log L for each character changing at this species' speciation.

    protected:
        void                            keepSpecialization(const DagNode *toucher) override;
        void                            restoreSpecialization(const DagNode *toucher) override;
        void                            swapParameterInternal(const DagNode *oldP, const DagNode *newP) override;

    private:
        enum Forcing { NONE, CHANGED, UNCHANGED };

        // a node's partial at the top of its branch, with everything it was computed from
        struct Slot {
            size_t                      version = 0;                                            //!< Unique per computation, 0 when empty.
            size_t                      epoch = 0;
            double                      length = 0.0;
            size_t                      kept = 0;
            long                        budded = -1;
            size_t                      xi = 0;
            size_t                      child_version[2] = {0, 0};
            std::vector<double>         v;
            std::vector<double>         ln_scale;                                               //!< Per character, summed over the subtree.
        };

        void                            resolveProcess(void);
        bool                            prepare(void);
        void                            prepareCache(void);
        const Slot&                     evaluate(const TopologyNode &node);
        void                            combine(const Slot &a, const Slot &b, long j, Forcing forcing, Slot &out) const;
        std::vector<double>             rootLnLikelihoods(const Slot &root) const;
        std::vector<double>             forcedLnLikelihoods(size_t taxon, Forcing forcing);
        void                            carry(std::vector<double> &v, double length) const;
        long                            speciesOf(const TopologyNode &node) const;
        long                            assignSpecies(const TopologyNode &node);
        void                            buildTipPartials(void);
        void                            simulateDown(const TopologyNode &node, std::vector<size_t> &state, std::vector< std::vector<size_t> > &tips);

        const TypedDagNode<Tree>*       tree;
        FossilizedBirthDeathSpeciationProcess* process;
        std::vector<size_t>             num_states;                                             //!< K for each character.
        bool                            states_given;                                           //!< K stated, not read off the data.

        std::vector<size_t>             changed;                                                //!< The character that changed at each species' speciation.
        size_t                          stored_taxon = 0;
        size_t                          stored_character = 0;
        bool                            reallocated = false;

        std::vector<double>             first;                                                  //!< tau_1 of each taxon, refreshed per computation.
        std::vector<double>             stem;                                                   //!< Markov chain length of each node's branch, by node index.
        std::vector<long>               species;                                                //!< The species running through each node.
        std::vector<size_t>             kept_child;                                             //!< The continuing child of each internal node.
        std::vector<long>               budded;                                                 //!< The species each internal node buds.

        // every character's states laid end to end, character c at offset[c]
        std::vector<size_t>             offset;
        size_t                          total_states = 0;
        std::vector< std::vector<double> > tip_partials;                                        //!< The observed states, by tip index.
        bool                            tips_stale = true;

        // two slots per node, so the state before a rejected proposal is usually still there
        std::vector< std::vector<Slot> > cache;
        std::vector<size_t>             current;                                                //!< The slot each node last used.
        size_t                          next_version = 1;
        size_t                          epoch = 1;                                              //!< Bumped when the tip data change.
    };
}

#endif
