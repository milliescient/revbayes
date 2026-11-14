#ifndef AbstractFossilizedBirthDeathRangeProcess_H
#define AbstractFossilizedBirthDeathRangeProcess_H

#include <algorithm>
#include <memory>
#include <set>

#include "RbConstants.h"
#include "FossilizedBirthDeathRangeRateClassClosedForm.h"
#include "FossilizedBirthDeathRangeRateClassNumeric.h"
#include "MemberObject.h"
#include "RbVector.h"
#include "Simplex.h"
#include "Taxon.h"
#include "TimeInterval.h"
#include "TypedDistribution.h"
#include "TypedDagNode.h"


namespace RevBayesCore {

    /**
     * @brief One taxon's stratigraphic range: the ages the model samples and the bounds its
     * occurrence record fixes them within.
     */
    struct RangeEntry {

        double  death = 0.0;                                //!< The extinction time.
        double  last  = 0.0;                                //!< The youngest appearance, tau_K.
        double  first = 0.0;                                //!< The oldest appearance, tau_1.
        double  divergence = 0.0;                           //!< Where the range attaches to the nearest sampled range.

        double  first_min = 0.0;                            //!< The oldest reported minimum.
        double  first_max = RbConstants::Double::inf;       //!< The oldest reported maximum.
        double  last_min  = 0.0;                            //!< The youngest reported minimum.
        double  last_max  = RbConstants::Double::inf;       //!< The youngest reported maximum.

        bool    singleton = true;                           //!< Fewer than two occurrences, so the appearances are one age.

        bool    survived = false;                           //!< Alive at the present. Only mvExtinctionRJSwitch changes it.

        size_t  component = 0;                              //!< Its psi mixture component. Only mvPsiMixtureAllocation changes it.

        //! Ages never decrease into the past, and each appearance sits within the bins reporting it.
        bool    isOrdered(double present) const
                {
                    if ( !( divergence > first && first >= last && last >= death && death >= present ) ) return false;

                    return first >= first_min && first <= first_max
                           && last >= last_min && last <= last_max;
                }

    };


    /**
     * @brief Abstract piecewise-constant fossilized birth-death range process.
     *
     * Rates are constant within each time interval and may jump at its boundaries.
     * Fossils are sampled at rate psi, and extant taxa with probability rho at the present only.
     * Each rate vector has one element per interval, and the timeline starts at the present.
     *
     * Each derived process provides updateRanges(), which reads each taxon's origination and
     * extinction times from its matrix or tree value. The likelihood is otherwise shared.
     *
     * @copyright Copyright 2009-
     * @author The RevBayes Development Core Team (Sebastian Hoehna)
     * @since 2014-03-18, version 1.0
     *
     */
    class AbstractFossilizedBirthDeathRangeProcess : public MemberObject< RbVector<double> >, public MemberObject<double> {
        
    public:
        AbstractFossilizedBirthDeathRangeProcess(const DagNode *speciation,
                                            const DagNode *extinction,
                                            const DagNode *psi,
                                            const TypedDagNode<double>* rho,
                                            const TypedDagNode<RbVector<double> > *times,
                                            const std::string &condition,
                                            const std::vector<Taxon> &taxa,
                                            bool complete_record,
                                            const TypedDagNode<double>* origin = NULL,
                                            double present = 0.0,
                                            const TypedDagNode<RbVector<double> >* psi_mixture_rates = NULL,
                                            const TypedDagNode<Simplex>* psi_mixture_weights = NULL,
                                            const TypedDagNode<RbVector<double> >* speciation_timeline = NULL,
                                            const TypedDagNode<RbVector<double> >* extinction_timeline = NULL,
                                            const TypedDagNode<RbVector<double> >* sampling_timeline = NULL);  //!< Constructor

        virtual ~AbstractFossilizedBirthDeathRangeProcess(){};

        void                                            setHasResampleMove(void) { has_resample_move = true; }       //!< Registered by mvStratigraphicRange.
        void                                            setHasReportingNode(void) { has_reporting_node = true; }     //!< Registered by dnFossilRecord.
        void                                            setHasSurvivalMove(void) { has_survival_move = true; }       //!< Registered by mvExtinctionRJSwitch.
        void                                            setHasAllocationMove(void) { has_allocation_move = true; }   //!< Registered by mvPsiMixtureAllocation.
        void                                            setHasSingletonMove(void) { has_singleton_move = true; }     //!< Registered by mvSingletonRJSwitch.

        void                                            executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const;   //!< Expose the appearances and the origination/extinction times.
        void                                            executeMethod(const std::string &n, const std::vector<const DagNode*> &args, double &rv) const;              //!< Expose the origin.
        double                                          getPresent() const { return times.front(); }               //!< The youngest interval boundary, which extinction times cannot precede.
        void                                            currentExtinctionAndSampling(std::vector<double> &t, std::vector<double> &mu, std::vector<double> &psi) const; //!< Interval starts, present first, and mu and psi as the rate nodes hold them now.
        void                                            currentRates(std::vector<double> &t, std::vector<double> &lambda, std::vector<double> &mu, std::vector<double> &psi) const; //!< Interval starts, present first, and all three rates as the rate nodes hold them now.
        const std::vector<Taxon>&                       getTaxa() const { return taxa; }
        void                                            resampleFirstLast(size_t i);                                //!< Redraw taxon i's appearances within its reported bins.
        bool                                            hasSurvived(size_t i) const { return ranges[i].survived; }  //!< Taxon i's latent status: alive at the present.
        double                                          getLastAppearance(size_t i) const { return ranges[i].last; } //!< tau_K.
        void                                            switchSurvived(size_t i);                                   //!< Flip that status, arming the undo. Only mvExtinctionRJSwitch may.
        void                                            warnIfNoResampleMove(void) const;                           //!< Warn once when the appearances have no move.
        double                                          getFirstAppearance(size_t i) const { return ranges[i].first; } //!< tau_1.
        bool                                            samplesSingletons(void) const { return record_reporting == PRESENCE; } //!< Only presence leaves a one-occurrence record free to hold more.
        double                                          switchSingleton(size_t i);                                  //!< Expand tau_1 above tau_K or collapse it onto tau_K, arming the undo. Only mvSingletonRJSwitch may.
        size_t                                          getNumberOfPsiMixtureComponents(void) const { return num_components; }
        size_t                                          getPsiMixtureComponent(size_t i) const { return ranges[i].component; }
        void                                            setPsiMixtureComponent(size_t i, size_t h);                 //!< Set taxon i's component, arming the undo. Only mvPsiMixtureAllocation may.
        std::vector<double>                             psiMixtureComponentLnProbs(size_t i) const;                 //!< Taxon i's log-density under each component, up to a constant.
        void                                            warnIfNoReportingNode(void) const;                          //!< Warn once when psi has no data.
        void                                            warnIfImpossibleInputs(const TypedDagNode<double> *start = NULL) const;  //!< Warn once when fixed inputs make the density zero. start stands in for a missing origin_age.
        void                                            warnIfInertNTaxaCondition(void) const;                      //!< Warn once when rho leaves nTaxa conditioning inert.
        void                                            warnIfNoSurvivalMove(void) const;                           //!< Warn once when rho < 1 and no move samples survival.
        void                                            warnIfNoAllocationMove(void) const;                         //!< Warn once when a psi mixture's components have no move.
        void                                            warnIfNoSingletonMove(void) const;                          //!< Warn once when a presence record has one occurrence and no move samples singleton status.
        void                                            initializeFirstLast(size_t i);                              //!< Draw taxon i's appearances, nested within the range and its reported bins.
        void                                            drawRanges();                                               //!< Draw an initial range and appearances for every taxon.
        double                                          computeLnFossilTotal();                                     //!< Total fossil-record log-density, for a dnFossilRecord node.
        enum Reporting { COMPLETE, FIRSTLAST, PRESENCE };                                                           //!< Which occurrences the record reports.
        void                                            setReporting(Reporting r);                                  //!< The reporting model, pushed here by dnFossilRecord.

        virtual void                                    setOccurrences(const std::vector<Taxon> &t);                //!< Adopt a clamped record's occurrences as the data.
        const std::vector<const DagNode*>&              getRangeParameters(void) const { return range_parameters; } //!< The rate and timeline nodes a dnFossilRecord adopts as parents.

    protected:
        virtual double                                  rangeEndTerm(size_t i, size_t di, double d, double p_d = -1.0, double pbar_d = -1.0) const { return log( death[di] ); } //!< Log term closing range i at d. Under a mixture p_d and pbar_d are its component's and a fresh draw's p at d.
        virtual void                                    mixtureSymmetric(std::vector<double> &beta) const { beta.clear(); }  //!< Each interval's chance of symmetric speciation, empty where there is none.
        virtual void                                    mixtureAnagenetic(std::vector<double> &rate) const { rate.clear(); }  //!< Each interval's anagenetic rate, empty where there is none.
        virtual void                                    mixtureExtraLnProbs(size_t i, std::vector<double> &ln_probs) const {}  //!< Terms beside the ranges that depend on taxon i's component, by component.
        const FossilizedBirthDeathRangeRateClassSolver&                         mixtureSolver(void) const { if ( mixture_numeric_on == true ) return mixture_numeric; return mixture_solution; }  //!< The solver holding the mixture, closed form at beta = 0 and no anagenesis.
        void                                            mixtureAt(double t, std::vector<double> &p, double &lambda_pbar) const;  //!< The mixture solution's p and lambda integral.
        void                                            mixtureSensitivity(double s, double t, const std::vector<size_t> &columns, std::vector<double> &Q, std::vector<double> &Q_scale) const;
        virtual double                                  rangeLnProb(size_t i);                                     //!< Birth-death density of range i.
        virtual bool                                    chargesBirth(size_t i) const;                              //!< Whether range i began with a speciation, whose rate it pays.
        virtual double                                  conditionLnProb(void) const;                               //!< The conditioning normalization, where it applies to the process as a whole.

        virtual void                                    repairRanges(void) {}                                       //!< Re-establish the invariants this process's own moves can break.
        virtual void                                    updateRanges() = 0;                                         //!< Pull the ranges from the value.

        virtual double                                  computeLnProbabilityRanges(bool force = false);
        double                                          computeLnFossilRecord(size_t i) const;                      //!< Fossil-record log-term for taxon i.
        double                                          computeLnFossilRecord(size_t i, double mixture_rate) const; //!< The same, with psi scaled by a mixture rate.
        void                                            mixtureRangeLnProbs(size_t i, std::vector<double> &ln_probs, size_t only, const std::vector<double> *source = NULL) const;  //!< Birth-death density of range i under each component, or under component only, its lineage starting in source.
        virtual void                                    mixtureStemSource(size_t i, std::vector<double> &w) const;  //!< The probability that range i's lineage starts in each component.
        virtual std::vector<size_t>                     mixtureDependents(size_t i) const { return std::vector<size_t>(); }  //!< The ranges whose lineage starts in range i's component.
        virtual double                                  observedRangeRate(size_t j) const { return 0.0; }           //!< A rate in interval j at which an observed range must have no event, beside lambda and mu.
        virtual void                                    appendIntervalRates(size_t j, std::vector<double> &v) const;  //!< Every rate of interval j a range term reads.
        void                                            storeRangeTerms(void);                                      //!< Keep the range terms and what they were computed under, for a rejected proposal.
        void                                            dirtyChangedIntervals(void);                                //!< Dirty the taxa born in or above the youngest interval whose rates or top edge changed.

        // Parameter management functions
        void                                            swapParameterInternal(const DagNode *oldP, const DagNode *newP);                //!< Swap a parameter

        // helper functions
        size_t                                          findIndex(double t) const;                             //!< The index so that times[index] <= t < times[index+1]
        double                                          p(size_t i, double t, bool survival = false) const;
        void                                            numTaxaGeometric(double &alpha, double &beta) const;       //!< The count at the present: geometric with an atom at alpha.
        size_t                                          getNumberOfTaxaAtPresent(void) const;                      //!< The taxa the record reports as extant.
        double                                          lnProbNumTaxa(size_t n) const;                             //!< log P(n taxa at the present).
        virtual double                                  q(size_t i, double t, bool tilde = false) const;
        double                                          lnQTildeFromPresent(size_t i, double t) const;              //!< log q_tilde from the present up to t.
        double                                          lnQFromPresent(size_t i, double t) const;                   //!< log q from the present up to t.

        virtual double                                  lnOriginationSurvival(size_t i, double t) const;            //!< log P(species i originated above t).
        double                                          drawOrigination(size_t i) const;                     //!< Draw species i's own origination from its conditional.
        double                                          meanOrigination(size_t i) const;                      //!< The conditional mean of species i's origination.

        virtual void                                    prepareProbComputation(void) const;
        void                                            prepareTimeline(std::vector<double> &global_timeline, std::vector<double> &lambda, std::vector<double> &mu, std::vector<double> &psi) const; //!< With per-rate timelines, as in BDSTP: their union, and each rate on it.
        void                                            addTimesToGlobalTimeline(std::multiset<double> &event_times, const TypedDagNode<RbVector<double> > *par_times) const; //!< Adds a rate's boundaries to those the global timeline is made from.
        void                                            expandNonGlobalRateParameterVector(std::vector<double> &par, const std::vector<double> &par_times, const std::vector<double> &global_timeline) const; //!< Puts par, given on its own timeline, on the global one.
        std::vector<double>                             intervalStarts(const TypedDagNode<RbVector<double> > *par_times) const; //!< The present, then each boundary above it.
        size_t                                          findIndex(double t, const std::vector<double> &timeline) const; //!< The index so that timeline[index] <= t < timeline[index+1].

        void                                            keepSpecialization(const DagNode *toucher);                  /* NOT VIRTUAL */
        void                                            restoreSpecialization(const DagNode *toucher);               /* NOT VIRTUAL */
        void                                            touchSpecialization(const DagNode *toucher, bool touchAll);  /* NOT VIRTUAL */

        std::vector<Taxon>                              taxa;                                                  //!< Taxa that will be attached to new simulated trees.
        std::string                                     condition;
        Reporting                                       record_reporting;                                       //!< Which occurrences the record reports.
        double                                          lnPresenceRecord(size_t i, double mixture_rate) const;  //!< Each reported interval holds at least its count.
        double                                          samplingExposure(double lo, double hi) const;           //!< The sampling rate integrated over [lo, hi).
        void                                            checkPresenceIntervals(void) const;                     //!< Under PRESENCE, refuse a taxon whose intervals overlap.
        double                                          max_present_age;                                        //!< Youngest occurrence maximum. The timeline may not start above it.

        void                                            updateRecord(void);                                     //!< Refresh the reported bins and the bounds they put on the appearances.

        size_t                                          num_intervals;

        // members
        const TypedDagNode<double >*                    homogeneous_lambda;                                    //!< The homogeneous speciation rate.
        const TypedDagNode<RbVector<double> >*          heterogeneous_lambda;                                  //!< The heterogeneous speciation rates.
        const TypedDagNode<double >*                    homogeneous_mu;                                        //!< The homogeneous extinction rate.
        const TypedDagNode<RbVector<double> >*          heterogeneous_mu;                                      //!< The heterogeneous extinction rates.
        const TypedDagNode<double >*                    homogeneous_psi;                                       //!< The homogeneous fossil sampling rate.
        const TypedDagNode<RbVector<double> >*          heterogeneous_psi;                                     //!< The heterogeneous fossil sampling rates.
        const TypedDagNode<double >*                    homogeneous_rho;                                       //!< The sampling probability at the present.
        const TypedDagNode<RbVector<double> >*          timeline;                                              //!< The interval boundaries.
        const TypedDagNode<RbVector<double> >*          interval_times_speciation;                              //!< lambda's own boundaries, in place of the global timeline.
        const TypedDagNode<RbVector<double> >*          interval_times_extinction;                              //!< mu's own boundaries.
        const TypedDagNode<RbVector<double> >*          interval_times_sampling;                                //!< psi's own boundaries.
        bool                                            using_rate_timelines = false;                           //!< Whether the rates have their own boundaries instead of the global ones.
        const TypedDagNode<double >*                    origin_age;                                            //!< The optional process start (NULL = the oldest birth, under no prior).
        const TypedDagNode<RbVector<double> >*          psi_mixture_rates;                                     //!< The relative psi of each mixture component (NULL = one component).
        const TypedDagNode<Simplex>*                    psi_mixture_weights;                                   //!< The mixture weights (NULL = equal).
        size_t                                          num_components;

        std::vector<const DagNode*>                     range_parameters;

        std::vector<RangeEntry>                         ranges;                                                 //!< One per taxon.
        
        double                                          origin;                                                 //!< The origin time (oldest birth time)
        bool                                            origin_unsampled = false;                               //!< No taxon samples the range that begins at the origin.
        size_t                                          max_birth;                                              //!< Index of the taxon holding the oldest divergence time.

        // the following vectors are used internally for more efficient likelihood calculations and are filled by 'prepareProbComputation'
        mutable std::vector<double>                     birth;                                                  //!< The sorted speciation rates
        mutable std::vector<double>                     death;                                                  //!< The sorted extinction rates
        mutable std::vector<double>                     fossil;                                                 //!< The sorted fossil sampling rates
        mutable std::vector<double>                     times;                                                  //!< The sorted interval times
                        
        mutable std::vector<double>                     q_i;                                                    //!< Probability of no speciation or extinction event in each time interval
        mutable std::vector<double>                     q_tilde_i;                                              //!< Probability of no change in species identity in each time interval
        mutable std::vector<double>                     p_i;                                                    //!< Probability of leaving no sampled descendants from the end of each time interval
        mutable std::vector<double>                     pS_i;                                                   //!< Probability of leaving no descendants from the end of each time interval

        // used only with more than one psi mixture component
        mutable std::vector<double>                     mixture_rates;                                          //!< The relative fossil sampling rate of each component
        mutable std::vector<double>                     mixture_weights;                                        //!< The probability of each component
        mutable FossilizedBirthDeathRangeRateClassNumeric                       mixture_numeric;                                        //!< Used when beta or the anagenetic rate is above 0
        mutable bool                                    mixture_numeric_on = false;                             //!< Whether the mixture is solved numerically.
        mutable std::vector<double>                     mix_beta;                                               //!< scratch: each interval's beta, for the solver
        mutable std::vector<double>                     mix_ana;                                                //!< scratch: each interval's anagenetic rate, for the solver
        mutable FossilizedBirthDeathRangeRateClassClosedForm                      mixture_solution;                                       //!< p, int lambda pbar and Q of every component
        mutable std::vector<double>                     mix_start, mix_p, mix_Q, mix_Q_scale, mix_ln_probs;     //!< scratch for the mixture range density
        mutable std::vector<size_t>                     mix_columns;

                                
        size_t                                          stored_range = 0;                                       //!< The range mvStratigraphicRange drew.
        double                                          stored_first = 0.0;                                     //!< The ages it replaced.
        double                                          stored_last  = 0.0;

        size_t                                          switched_range = 0;                                     //!< The range a switch flipped.
        bool                                            stored_survived = false;                                //!< The statuses before the switch.
        bool                                            stored_singleton = true;
        bool                                            switched = false;                                       //!< That undo is armed.
        size_t                                          reallocated_range = 0;                                  //!< The range mvPsiMixtureAllocation moved, and the component it replaced.
        size_t                                          stored_component = 0;
        bool                                            reallocated = false;

        bool                                            weld_origin = true;                                     //!< Only the origin node drags the oldest birth.

        std::vector<RangeEntry>                         stored_ranges;                                          //!< The table as the proposal found it, restored when it is rejected.

        std::vector<double>                             partial_likelihood;                                     //!< Partial likelihood for each taxon
        std::vector<double>                             stored_likelihood;                                      //!< Stored partial likelihood for each taxon
        std::vector<double>                             range_shared;                                           //!< What every interval shares, as partial_likelihood was computed under it
        std::vector<double>                             range_rates;                                            //!< Each interval's rates, likewise
        std::vector<double>                             stored_range_shared;
        std::vector<double>                             stored_range_rates;
                                
        std::vector<bool>                               dirty_taxa;                                             //!< Indicates whether partial likelihood needs updating

        //! A taxon's record term with the arguments it was computed from
        struct RecordTerm {
            double  death = 0.0, first = 0.0, last = 0.0, rate = 0.0, value = 0.0;
            bool    singleton = false, valid = false;
        };
        std::vector<RecordTerm>                         record_terms;
        std::vector<double>                             record_terms_fossil;                                    //!< The sampling rates and times every record term was computed under
        std::vector<double>                             record_terms_times;
        
        bool                                            touched;                                               //!< Indicates whether any terms need updating
        bool                                            has_resample_move = false;                              //!< Set by mvStratigraphicRange when it attaches.
        bool                                            has_reporting_node = false;                             //!< Set by dnFossilRecord when it attaches.
        bool                                            has_survival_move = false;                              //!< Set by mvExtinctionRJSwitch when it attaches.
        bool                                            has_singleton_move = false;                             //!< Set by mvSingletonRJSwitch when it attaches.
        bool                                            has_allocation_move = false;                            //!< Set by mvPsiMixtureAllocation when it attaches.
        mutable bool                                    warned_no_resample = false;
        mutable bool                                    warned_no_reporting = false;
        mutable bool                                    warned_inert_ntaxa = false;
        mutable bool                                    warned_no_survival = false;
        mutable bool                                    warned_no_singleton = false;
        mutable bool                                    warned_no_allocation = false;
        mutable bool                                    warned_impossible = false;
        bool                                            resampled;                                              //!< The appearance undo is armed.
        double                                          present_time = 0.0;                                     //!< The present, where the timeline starts.
    };
}

#endif
