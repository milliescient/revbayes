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
     * Without epsilon, each speciation changes one character, a latent index per speciation. With
     * epsilon, each speciation changes every character independently, and each branch carries latent
     * counts of hidden budding and anagenetic events. Latent variables are kept by node index.
     */
    class PhyloMorphospeciationProcess : public TypedDistribution<AbstractHomologousDiscreteCharacterData>, public MemberObject< RbVector<double> > {

    public:
        PhyloMorphospeciationProcess(const TypedDagNode<Tree> *tree, const std::vector<size_t> &num_states, const TypedDagNode<double> *epsilon = NULL, const TypedDagNode<double> *epsilon_a = NULL, bool variable = false, size_t shared_states = 0);
        virtual ~PhyloMorphospeciationProcess() {}

        PhyloMorphospeciationProcess*   clone(void) const override;
        double                          computeLnProbability(void) override;
        void                            redrawValue(void) override;
        void                            setValue(AbstractHomologousDiscreteCharacterData *v, bool force=false) override;

        void                            executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const override;
        std::string                     getHiddenStateString(void) const override;              //!< The latent variables, for checkpointing clamped characters.
        void                            setHiddenStateFromString(const std::string &s) override;

        size_t                          getNumberOfCharacters(void) const { return num_states.size(); }
        bool                            allowsMultipleChanges(void) const { return epsilon != NULL; }

        // Gibbs draws of the latent variables, each returning its hastings ratio and arming the undo
        double                          resampleChangedCharacter(void);                         //!< Without epsilon: the character that changed at one speciation, chosen uniformly.
        double                          resampleEvents(void);                                   //!< With epsilon: one kind of hidden event on one branch, both chosen uniformly.

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
            size_t                      xi = 0;                                                 //!< The changed character, or with epsilon whether the speciation is an anagenesis.
            size_t                      child_version[2] = {0, 0};
            std::vector<double>         v;
            std::vector<double>         ln_scale;                                               //!< Per character, summed over the subtree.
        };

        // a likelihood with the change at one node's speciation, or the matrix of the branch above it, overridden
        struct Forced {
            long            node = -1;
            Forcing         forcing = NONE;
            BranchForcing   branch = AS_IS;
        };

        void                            resolveProcess(void);
        bool                            prepare(void);
        void                            hiddenIntegrals(const std::vector<long> &sp, const std::vector<double> &tau, std::vector<double> &budding, std::vector<double> &anagenetic) const;  //!< int lambda p and int lambda_a above tau_1 on each node's branch.
        bool                            changeProbabilities(double &e_budding, double &e_anagenesis) const;                     //!< False outside [0, 1].
        long                            assignSpecies(const TopologyNode &node);
        void                            initializeLatents(void);
        bool                            ready(void);                                            //!< Rates, tree, latent variables and cache brought up to date. False if impossible.
        double                          drawEvents(size_t node);
        void                            prepareCache(void);
        void                            markCoded(void);                                        //!< Which nodes have a coded tip below them.
        std::vector<char>               codedNodes(const std::vector<char> &tip_coded) const;
        const Slot&                     slotOf(const TopologyNode &node) const;                 //!< A node's current partial, all ones where nothing below it is coded.
        bool                            summedOut(const TopologyNode &node) const;              //!< A tip coded for no character, whose branch's hidden counts never touch the likelihood.
        void                            signature(const TopologyNode &node, Slot &sig) const;
        const Slot&                     evaluate(const TopologyNode &node);
        void                            combine(const Slot &a, const Slot &b, long j, Forcing forcing, Slot &out) const;
        std::vector<double>             rootLnLikelihoods(const Slot &root) const;
        std::vector<double>             forcedLnLikelihoods(const Forced &f);
        void                            carry(std::vector<double> &v, size_t node, BranchForcing branch) const;
        void                            mix(std::vector<double> &v, size_t c, double e) const;
        double                          decay(size_t c, double e) const;                        //!< 1 - e K/(K-1), the eigenvalue of one event off the stationary direction.
        double                          lnVariable(size_t c, double ln_all_zero) const;         //!< log P(character c varies), from the log probability its scored tips are all state 0.
        bool                            changesAt(size_t node, size_t c) const;
        long                            speciesOf(const TopologyNode &node) const;
        void                            buildTipPartials(void);
        void                            simulateDown(const TopologyNode &node, std::vector<size_t> &state, std::vector< std::vector<size_t> > &tips);

        const TypedDagNode<Tree>*       tree;
        const TypedDagNode<double>*     epsilon;                                                //!< NULL: exactly one change per speciation.
        const TypedDagNode<double>*     epsilon_a;                                              //!< NULL: anagenesis changes with probability epsilon too.
        FossilizedBirthDeathSpeciationProcess* process;
        std::vector<size_t>             num_states;                                             //!< K for each character.
        bool                            states_given;                                           //!< K stated, not read off the data.
        size_t                          shared_states;                                          //!< K shared by every character, the count read off the data; 0 if not.
        bool                            variable;                                               //!< Only variable characters were kept, so each is conditioned on varying.

        std::vector<size_t>             changed;                                                //!< The character that changed at each interior node's speciation.
        size_t                          stored_node = 0;
        size_t                          stored_character = 0;
        bool                            reallocated = false;

        std::vector<long>               events[2];                                              //!< With epsilon: hidden events of each kind on the branch above each node.
        bool                            latents_ready = false;
        enum Undo { UNDO_NONE, UNDO_EVENTS };
        Undo                            undo = UNDO_NONE;
        size_t                          undo_index = 0;
        Mode                            undo_mode = BUDDING;
        bool                            fresh = false;                                          //!< Prepared by a Gibbs draw, with nothing since but the node's own touch.
        long                            undo_events = 0;

        double                          change_prob[2] = { 0.0, 0.0 };                          //!< With epsilon: each character's chance to change at a budding and at an anagenesis.
        std::vector<double>             first;                                                  //!< tau_1 of each taxon, refreshed per computation.
        std::vector<double>             stem;                                                   //!< Markov chain length of each node's branch, by node index.
        std::vector<long>               species;                                                //!< The species running through each node.
        std::vector<double>             mean_events[2];                                         //!< With epsilon: the Poisson mean of each kind of event on each node's branch.
        std::vector<size_t>             kept_child;                                             //!< The continuing child of each internal node.
        std::vector<long>               budded;                                                 //!< The species each internal node buds.
        std::vector<char>               anagenetic;                                             //!< Each interior node's speciation is an observed anagenesis.

        // every column's states laid end to end, column c at offset[c]. The columns are the characters and,
        // under variable coding, after them one per K and set of scored tips, with every scored tip in state 0
        std::vector<size_t>             columns;                                                //!< K for each column.
        std::vector<size_t>             all_zero;                                               //!< Under variable coding, the all-0 column each character shares.
        std::vector<size_t>             offset;
        size_t                          total_states = 0;
        size_t                          max_states = 2;
        std::vector< std::vector<double> > tip_partials;                                        //!< The observed states, by tip index.
        std::vector<char>               tip_coded;                                              //!< Each tip has a state at some character, by tip index.
        std::vector<char>               coded;                                                  //!< Each node has a coded tip below it, by node index. Elsewhere every
                                                                                                //!< partial is all ones and the hidden counts above are drawn from their prior.
        Slot                            ones;                                                   //!< The partial of a subtree with nothing coded.
        bool                            tips_stale = true;

        // two slots per node, so the state before a rejected proposal is usually still there
        std::vector< std::vector<Slot> > cache;
        std::vector<size_t>             current;                                                //!< The slot each node last used.
        size_t                          next_version = 1;
        size_t                          epoch = 1;                                              //!< Bumped when the tip data change.
    };
}

#endif
