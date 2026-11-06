#ifndef FossilizedBirthDeathSpeciationProcess_H
#define FossilizedBirthDeathSpeciationProcess_H

#include "AbstractFossilizedBirthDeathRangeProcess.h"
#include "AbstractBirthDeathProcess.h"

namespace RevBayesCore {

    /**
     * @brief What the topology says about one taxon's range: how it begins, and how it ends.
     */
    struct SpeciesEntry {

        bool parent_is_sa   = false;    //!< The range this one budded from was sampled at that node.
        bool ends_at_sa     = false;    //!< The range ends at a sampled ancestor, so its lineage carries on.
        bool ends_symmetric = false;    //!< The range ends at an observed symmetric speciation, which its node pays for.
        int  sa_parent      = -1;       //!< With parent_is_sa, the range sampled at that node.

        bool operator==(const SpeciesEntry &e) const
                {
                    return parent_is_sa == e.parent_is_sa && ends_at_sa == e.ends_at_sa
                           && ends_symmetric == e.ends_symmetric && sa_parent == e.sa_parent;
                }

    };

    /**
     * @brief Piecewise-constant fossilized birth-death species distribution of extended trees.
     *
     * The piecewise-constant fossilized birth-death species process has constant rates for each time interval.
     * At the end of each time interval there may be an abrupt rate-shift (jump) for each
     * of the rates. Additionally, there may be sampling at the end of each interval.
     * Finally, fossils are sampled with rate psi, the others (fossils and extant taxa) are
     * sampled at sampling times (including the present).
     *
     * We assume that the rate vectors have one more element than the rate-change vectors.
     * Thus, one rate-change means always two interval, two rate-changes three interval, and so on.
     *
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (June Walker)
     * @since 2014-03-18, version 1.0
     *
     */
    class FossilizedBirthDeathSpeciationProcess : public AbstractBirthDeathProcess, public AbstractFossilizedBirthDeathRangeProcess {
        
        using AbstractFossilizedBirthDeathRangeProcess::taxa;

    public:
        FossilizedBirthDeathSpeciationProcess (const TypedDagNode<double>* ra,
                                      const DagNode *speciation,
                                      const DagNode *extinction,
                                      const DagNode *psi,
                                      const TypedDagNode<double>* rho,
                                      const DagNode *lambda_a,
                                      const DagNode *beta,
                                      const TypedDagNode<RbVector<double> > *times,
                                      const std::string &condition,
                                      const std::vector<Taxon> &taxa,
                                      bool complete_record,
                                      bool extended = true,
                                      double present = 0.0,
                                      const TypedDagNode<RbVector<double> >* psi_mixture_rates = NULL,
                                      const TypedDagNode<Simplex>* psi_mixture_weights = NULL);  //!< Constructor
        
        // public member functions
        FossilizedBirthDeathSpeciationProcess*          clone(void) const override;                                //!< Create an independent clone

        void                                            setValue(Tree *v, bool force=false) override;                       //!< Adopt a tree's annotations, or draw the state it does not carry.

        std::string                                     getHiddenStateString(void) const override;                          //!< The tree annotated with the range labels, appearances and statuses a Newick drops.
        void                                            setHiddenStateFromString(const std::string &s) override;
        void                                            redrawValue(void) override;
        void                                            redrawValue(SimulationCondition c) override;                        //!< Redraw without reaching the inherited simulator
        bool                                            redrawTopology(void);                                               //!< Gibbs-redraw the budding topology with the ranges fixed. Pure budding only.
        bool                                            hasAnagenesis(void) const;                                          //!< True when any interval's anagenetic rate is positive.
        bool                                            hasSymmetricSpeciation(void) const;                                 //!< True when any interval's symmetric speciation probability is positive.
        void                                            simulateClade(std::vector<TopologyNode *> &n, double age, double present, bool alwaysReturn) override;
        bool                                            allowsSampledAncestors(void) const override { return true; }        //!< The density rejects them where lambda_a = 0.
        bool                                            tracksSpeciesContinuation(void) const override { return true; }
        double                                          hiddenSpeciationIntegral(double t0, double t1) const;               //!< int (lambda p + lambda_a) dt from age t0 up to t1.
        std::vector<double>                             hiddenSpeciationIntegrals(const std::vector<std::pair<double,double> > &spans) const;  //!< The same over many spans, preparing the rates once.
        std::vector<double>                             anageneticIntegrals(const std::vector<std::pair<double,double> > &spans) const;        //!< int lambda_a dt over each span, the anagenetic part of the above.
        std::vector<const DagNode*>                     getRateParameters(void) const;                                      //!< Every node the rates and p(t) depend on.
        void                                            executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const override;  //!< The per-branch speciation events, then the base's range methods.
        //!< By node index: the Poisson means of the hidden speciations on each branch, and the observed event that began the species.
        struct BranchEvents { std::vector<double> hidden_budding, hidden_bifurcation, hidden_anagenetic, observed_budding, observed_bifurcation, observed_anagenetic; };
        void                                            branchEvents(BranchEvents &e, bool hidden = true) const;
        mutable BranchEvents                            cached_events;                                           //!< The last branchEvents, which every accessor of one tree state shares.
        mutable int                                     cached_events_level = 0;                                 //!< 0 stale, 1 observed only, 2 hidden too
        double                                          drawExtinction(size_t i) const;                                  //!< On a sampled tree, a draw of species i's extinction from its conditional.
        void                                            hiddenEventIntegrals(const std::vector<std::pair<double,double> > &spans, std::vector<double> &budding, std::vector<double> &bifurcation, std::vector<double> &anagenesis) const;  //!< The hidden events that change a lineage's species, by mode.

    protected:
        bool                                            marginalizesExtinction(void) const override { return extended == false; }
        double                                          rangeEndTerm(size_t i, size_t di, double d, double p_d = -1.0, double pbar_d = -1.0) const override;
        void                                            mixtureSymmetric(std::vector<double> &beta) const override;
        void                                            mixtureExtraLnProbs(size_t i, std::vector<double> &ln_probs) const override;
        void                                            mixtureStemSource(size_t i, std::vector<double> &w) const override;  //!< An anagenetic species starts in its ancestor's component.
        std::vector<size_t>                             mixtureDependents(size_t i) const override;
        double                                          observedRangeRate(size_t j) const override { return anagenetic[j]; }  //!< No anagenesis inside an observed range.
        void                                            appendIntervalRates(size_t j, std::vector<double> &v) const override;

        void                                            repairRanges(void) override;                             //!< Re-take from the tip what a tree move carries there.
        void                                            updateRanges(void) override;
        void                                            setOccurrences(const std::vector<Taxon> &t) override;    //!< Sync the tree-side taxon copy, then adopt as the base does.
        double                                          symmetricAt(double age) const;                           //!< beta in the interval containing age, read from the parameter.
        void                                            normalizeContinuationFlags(void);                        //!< Repair the whole tree, refreshing the interval cache first.
        void                                            normalizeContinuationFlags(const TopologyNode &node);    //!< Make each node name exactly one continuing child.
        int                                             continuingSpecies(const TopologyNode &node) const;       //!< The taxon whose range this node belongs to, or -1 where none is named.

        //!< The range running through a subtree, or -1 while it is unnamed.
        struct RangeFlow { int species; double end_age; };
        RangeFlow                                       updateRanges(const TopologyNode & );

        //!< A range no taxon samples, from where it began to the symmetric speciation that ended it. At a sampled
        //!< ancestor it begins by anagenesis in an extended tree, and in a sampled tree it carries on that lineage.
        enum class RangeStart { SPECIATION, ANAGENESIS, SAMPLED_ANCESTOR };
        struct UnsampledRange { double top; double end; RangeStart start; int sa_species; };  //!< sa_species: the sampled ancestor it begins at, else -1
        double                                          speciationModeLnProb(void) const;                        //!< The mode of each speciation, and the ranges no taxon samples.
        double                                          unsampledMixtureLnProb(const UnsampledRange &r, const std::vector<double> &source) const;  //!< An unsampled range's path under the mixture, starting in source.

        double                                          pSurvival(double start, double end) const override;             //!< Compute the probability of survival of the process (without incomplete taxon sampling).

        // Parameter management functions
        double                                          computeLnProbabilityTimes(void) const override;                            //!< Compute the log-transformed probability of the current value.
        double                                          computeLnProbabilityDivergenceTimes(void) const override;            //!< Compute the log-transformed probability of the current value.

        bool                                            tipAgeConstrainedToRange(const Taxon &t) const override { return extended == false && t.isExtinct(); } //!< A non-extended extinct tip is its youngest appearance.
        bool                                            validatesTipAgesOnSet(void) const override { return false; }        //!< The density's per-taxon check knows which tips are exempt.

        void                                            setMcmcMode(bool tf) override;

        bool                                            isExtended(void) const override { return extended; }     //!< An extended tree ends each range at an extinction time.

        double                                          lnProbNumTaxa(size_t n, double start, double end, bool MRCA) const override { throw RbException("Cannot compute P(nTaxa)."); }
        double                                          lnProbTreeShape(void) const override;

        double                                          q(size_t i, double t, bool tilde = false) const override;

        double                                          simulateDivergenceTime(double origin, double present) const override;    //!< Simulate a speciation event.
        std::vector<double>                             simulateDivergenceTimes(size_t n, double origin, double present, double min, bool alwaysReturn) const override;                 //!< Simulate n speciation events.

        void                                            keepSpecialization(const DagNode *toucher) override;
        void                                            restoreSpecialization(const DagNode *toucher) override;
        void                                            touchSpecialization(const DagNode *toucher, bool touchAll) override;

        // Parameter management functions
        void                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP) override;                //!< Swap a parameter

        void                                            prepareProbComputation(void) const override;

    private:
        
        // helper functions
        double                                          getMaxTaxonAge( const TopologyNode& ) const;

        void                                            labelSpecies(TopologyNode &node, const std::string &s, size_t &unsampled_count) const;  //!< Stamp the range each node belongs to, top down.
        std::string                                     rangeLabel(const TopologyNode &node, size_t &unsampled_count) const;                  //!< A taxon's index, or a fresh id for a range no taxon samples.
        bool                                            adoptSpeciesLabels(void);                               //!< Take the continuation flags from a tree's range labels.
        bool                                            adoptAppearances(void);                                 //!< Take the appearances from the tips' FAD and LAD.

        bool                                            extended;                                                //!< Tips are extinction times, else those are marginalized out.

        std::vector<SpeciesEntry>                       species;                                                 //!< One per taxon, as the last pull left it.
        std::vector<SpeciesEntry>                       next_species;                                            //!< Where the tree pass builds the replacement.
        std::vector<SpeciesEntry>                       stored_species;                                          //!< The table as the proposal found it.

        std::vector<double>                             budding_ages;                                            //!< Each budding event, which pays 1-beta.
        std::vector<double>                             symmetric_ages;                                          //!< Each symmetric event, which pays 2 beta for both daughters' births.
        std::vector<UnsampledRange>                     unsampled;                                               //!< The ranges no taxon samples.
        mutable bool                                    invalid_continuation = false;                            //!< Some node's children do not name exactly one continuation.

        mutable std::vector<double>                     anagenetic;                                              //!< The sorted anagenetic speciation rates.
        mutable std::vector<double>                     symmetric;                                               //!< The sorted symmetric speciation probabilities.

        const TypedDagNode<double >*                    homogeneous_lambda_a;                                    //!< The homogeneous anagenetic speciation rates.
        const TypedDagNode<RbVector<double> >*          heterogeneous_lambda_a;                                  //!< The heterogeneous anagenetic speciation rates.
        const TypedDagNode<double >*                    homogeneous_beta;                                        //!< The homogeneous symmetric speciation probability.
        const TypedDagNode<RbVector<double> >*          heterogeneous_beta;                                      //!< The heterogeneous symmetric speciation probabilities.

    };
}

#endif
