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
     * @brief Morphological characters whose changes happen at the speciation events of a dnFBDSP tree.
     *
     * Without epsilon, every speciation changes the new species at one character, chosen uniformly, to
     * a state drawn from that character's Mk kernel. Conditional on the tree, a character evolves along
     * a species' stem as a Markov chain with generator J - I over the length (1/n) int (lambda p + lambda_a) dt,
     * does not change over the species' sampled range, and starts from the kernel's stationary
     * distribution at the root. Which character changed at each speciation is a latent index per
     * species, sampled by mvMorphospeciationAllocation.
     *
     * With epsilon, a speciation changes each character independently, so it can change several or
     * none. Each branch carries latent counts of hidden budding and hidden anagenetic events, Poisson
     * with means int lambda p dt and int lambda_a dt, and given the counts the characters are
     * independent, each carried by ((1 - e) I + e J)^N for each kind of event. The counts are sampled by
     * mvMorphospeciationEvents, and the characters that changed at each observed speciation by
     * mvMorphospeciationAllocation.
     *
     * With epsilon_a, anagenetic speciations change each character with probability epsilon_a and
     * buddings with epsilon, so epsilon / epsilon_a is how much change a budding carries relative to an
     * anagenesis.
     */
    class PhyloMorphospeciationProcess : public TypedDistribution<AbstractHomologousDiscreteCharacterData>, public MemberObject< RbVector<double> > {

    public:
        PhyloMorphospeciationProcess(const TypedDagNode<Tree> *tree, const std::vector<size_t> &num_states, const TypedDagNode<double> *epsilon = NULL, const TypedDagNode<double> *epsilon_a = NULL);
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
        bool                            allowsMultipleChanges(void) const { return epsilon != NULL; }
        bool                            hasSpeciationEvent(size_t taxon) const;                 //!< False for the species that begins at the origin.
        size_t                          getChangedCharacter(size_t taxon) const { return changed[taxon]; }
        void                            setChangedCharacter(size_t taxon, size_t c);            //!< Arms the undo. Only mvMorphospeciationAllocation may.
        std::vector<double>             changedCharacterLnProbs(size_t taxon);                  //!< log L for each character changing at this species' speciation.

        // with epsilon: Gibbs draws of the latent variables, each returning its hastings ratio and arming the undo
        double                          resampleChangedSet(void);                               //!< The characters that changed at one species' speciation, chosen uniformly.
        double                          resampleEvents(void);                                   //!< One kind of hidden event on one branch, both chosen uniformly.

    protected:
        void                            keepSpecialization(const DagNode *toucher) override;
        void                            restoreSpecialization(const DagNode *toucher) override;
        void                            swapParameterInternal(const DagNode *oldP, const DagNode *newP) override;
        void                            touchSpecialization(const DagNode *toucher, bool touchAll) override;

    private:
        enum Forcing { NONE, CHANGED, UNCHANGED };
        enum BranchForcing { AS_IS, IDENTITY, RANDOMIZED };
        enum Mode { BUDDING = 0, ANAGENESIS = 1 };

        // a node's partial at the top of its branch, with everything it was computed from
        struct Slot {
            size_t                      version = 0;                                            //!< Unique per computation, 0 when empty.
            size_t                      epoch = 0;
            double                      length = 0.0;
            long                        events[2] = {0, 0};
            double                      change_prob[2] = {0.0, 0.0};
            size_t                      kept = 0;
            long                        budded = -1;
            size_t                      xi = 0;                                                 //!< The changed character, or with epsilon the changed set's version.
            size_t                      child_version[2] = {0, 0};
            std::vector<double>         v;
            std::vector<double>         ln_scale;                                               //!< Per character, summed over the subtree.
        };

        // a likelihood with one species' change or one branch's matrix overridden
        struct Forced {
            long            taxon = -1;
            Forcing         forcing = NONE;
            long            node = -1;
            BranchForcing   branch = AS_IS;
        };

        void                            resolveProcess(void);
        bool                            prepare(void);
        void                            hiddenIntegrals(const std::vector<long> &sp, const std::vector<double> &tau, std::vector<double> &budding, std::vector<double> &anagenetic) const;  //!< int lambda p and int lambda_a above tau_1 on each node's branch.
        bool                            changeProbabilities(double &e_budding, double &e_anagenesis) const;                     //!< False outside [0, 1].
        size_t                          branchKey(const TopologyNode &node) const;
        long                            assignSpecies(const TopologyNode &node);
        void                            initializeLatents(void);
        bool                            ready(void);                                            //!< Rates, tree, latent variables and cache brought up to date. False if impossible.
        double                          drawChangedSet(size_t taxon);
        double                          drawEvents(size_t node);
        double                          lnSetProbability(const std::vector<char> &set, double e) const;
        void                            prepareCache(void);
        void                            signature(const TopologyNode &node, Slot &sig) const;
        const Slot&                     evaluate(const TopologyNode &node);
        void                            combine(const Slot &a, const Slot &b, long j, Forcing forcing, Slot &out) const;
        std::vector<double>             rootLnLikelihoods(const Slot &root) const;
        std::vector<double>             forcedLnLikelihoods(const Forced &f);
        void                            carry(std::vector<double> &v, size_t node, BranchForcing branch) const;
        void                            mix(std::vector<double> &v, size_t c, double e) const;
        double                          decay(size_t c, double e) const;                        //!< 1 - e K/(K-1), the eigenvalue of one event off the stationary direction.
        bool                            changesAt(long taxon, size_t c) const;
        long                            speciesOf(const TopologyNode &node) const;
        void                            buildTipPartials(void);
        void                            simulateDown(const TopologyNode &node, std::vector<size_t> &state, std::vector< std::vector<size_t> > &tips);

        const TypedDagNode<Tree>*       tree;
        const TypedDagNode<double>*     epsilon;                                                //!< NULL: exactly one change per speciation.
        const TypedDagNode<double>*     epsilon_a;                                              //!< NULL: anagenesis changes with probability epsilon too.
        FossilizedBirthDeathSpeciationProcess* process;
        std::vector<size_t>             num_states;                                             //!< K for each character.
        bool                            states_given;                                           //!< K stated, not read off the data.

        std::vector<size_t>             changed;                                                //!< The character that changed at each species' speciation.
        size_t                          stored_taxon = 0;
        size_t                          stored_character = 0;
        bool                            reallocated = false;

        std::vector< std::vector<char> > changed_set;                                           //!< With epsilon: the characters that changed at each species' speciation.
        std::vector<long>               events[2];                                              //!< With epsilon: hidden events of each kind on each branch, by branch key.
        bool                            latents_ready = false;
        enum Undo { UNDO_NONE, UNDO_SET, UNDO_EVENTS };
        Undo                            undo = UNDO_NONE;
        size_t                          undo_index = 0;
        Mode                            undo_mode = BUDDING;
        std::vector<char>               undo_set;
        size_t                          undo_version = 0;
        std::vector<size_t>             set_version;                                            //!< With epsilon: bumped whenever a species' changed set changes, for the cache.
        size_t                          next_set_version = 1;
        bool                            fresh = false;                                          //!< Prepared by a Gibbs draw, with nothing since but the node's own touch.
        long                            undo_events = 0;

        double                          change_prob[2] = { 0.0, 0.0 };                          //!< With epsilon: each character's chance to change at a budding and at an anagenesis.
        std::vector<double>             first;                                                  //!< tau_1 of each taxon, refreshed per computation.
        std::vector<double>             stem;                                                   //!< Markov chain length of each node's branch, by node index.
        std::vector<long>               species;                                                //!< The species running through each node.
        std::vector<double>             mean_events[2];                                         //!< With epsilon: the Poisson mean of each kind of event on each node's branch.
        std::vector<size_t>             key_of;                                                 //!< The branch key of each node.
        std::vector<size_t>             kept_child;                                             //!< The continuing child of each internal node.
        std::vector<long>               budded;                                                 //!< The species each internal node buds.
        std::vector<char>               anagenetic_species;                                     //!< Each species arose by observed anagenesis.

        // every character's states laid end to end, character c at offset[c]
        std::vector<size_t>             offset;
        size_t                          total_states = 0;
        size_t                          max_states = 2;
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
