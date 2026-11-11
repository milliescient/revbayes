#ifndef AbstractFossilizedBirthDeathRangeProcess_H
#define AbstractFossilizedBirthDeathRangeProcess_H

#include <algorithm>
#include <memory>

#include "RbConstants.h"
#include "MemberObject.h"
#include "RbVector.h"
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
                                            double present = 0.0);  //!< Constructor

        virtual ~AbstractFossilizedBirthDeathRangeProcess(){};

        void                                            setHasResampleMove(void) { has_resample_move = true; }       //!< Registered by mvStratigraphicRange.
        void                                            setHasReportingNode(void) { has_reporting_node = true; }     //!< Registered by dnFossilRecord.
        void                                            setHasSurvivalMove(void) { has_survival_move = true; }       //!< Registered by mvExtinctionRJSwitch.
        void                                            setHasSingletonMove(void) { has_singleton_move = true; }     //!< Registered by mvSingletonRJSwitch.

        void                                            executeMethod(const std::string &n, const std::vector<const DagNode*> &args, RbVector<double> &rv) const;   //!< Expose the appearances and the origination/extinction times.
        void                                            executeMethod(const std::string &n, const std::vector<const DagNode*> &args, double &rv) const;              //!< Expose the origin.
        double                                          getPresent() const { return times.front(); }               //!< The youngest interval boundary, which extinction times cannot precede.
        const std::vector<Taxon>&                       getTaxa() const { return taxa; }
        void                                            resampleFirstLast(size_t i);                                //!< Redraw taxon i's appearances within its reported bins.
        bool                                            hasSurvived(size_t i) const { return ranges[i].survived; }  //!< Taxon i's latent status: alive at the present.
        double                                          getLastAppearance(size_t i) const { return ranges[i].last; } //!< tau_K.
        void                                            switchSurvived(size_t i);                                   //!< Flip that status, arming the undo. Only mvExtinctionRJSwitch may.
        void                                            warnIfNoResampleMove(void) const;                           //!< Warn once when the appearances have no move.
        double                                          getFirstAppearance(size_t i) const { return ranges[i].first; } //!< tau_1.
        bool                                            samplesSingletons(void) const { return record_reporting == PRESENCE; } //!< Only presence leaves a one-occurrence record free to hold more.
        double                                          switchSingleton(size_t i);                                  //!< Expand tau_1 above tau_K or collapse it onto tau_K, arming the undo. Only mvSingletonRJSwitch may.
        void                                            warnIfNoReportingNode(void) const;                          //!< Warn once when psi has no data.
        void                                            warnIfImpossibleInputs(const TypedDagNode<double> *start = NULL) const;  //!< Warn once when fixed inputs make the density zero. start stands in for a missing origin_age.
        void                                            warnIfInertNTaxaCondition(void) const;                      //!< Warn once when rho leaves nTaxa conditioning inert.
        void                                            warnIfNoSurvivalMove(void) const;                           //!< Warn once when rho < 1 and no move samples survival.
        void                                            warnIfNoSingletonMove(void) const;                          //!< Warn once when a presence record has one occurrence and no move samples singleton status.
        void                                            initializeFirstLast(size_t i);                              //!< Draw taxon i's appearances, nested within the range and its reported bins.
        void                                            drawRanges();                                               //!< Draw an initial range and appearances for every taxon.
        double                                          computeLnFossilTotal();                                     //!< Total fossil-record log-density, for a dnFossilRecord node.
        enum Reporting { COMPLETE, FIRSTLAST, PRESENCE };                                                           //!< Which occurrences the record reports.
        void                                            setReporting(Reporting r);                                  //!< The reporting model, pushed here by dnFossilRecord.

        virtual void                                    setOccurrences(const std::vector<Taxon> &t);                //!< Adopt a clamped record's occurrences as the data.
        const std::vector<const DagNode*>&              getRangeParameters(void) const { return range_parameters; } //!< The rate and timeline nodes a dnFossilRecord adopts as parents.

    protected:
        virtual double                                  rangeEndTerm(size_t i, size_t di, double d) const { return log( death[di] ); } //!< Log term closing range i at d.
        virtual double                                  rangeLnProb(size_t i);                                     //!< Birth-death density of range i.
        double                                          sightingLnProb(size_t i) const;                            //!< Log density of range i's first and last appearances.
        virtual bool                                    chargesBirth(size_t i) const;                              //!< Whether range i began with a speciation, whose rate it pays.
        virtual double                                  conditionLnProb(void) const;                               //!< The conditioning normalization, where it applies to the process as a whole.

        virtual void                                    repairRanges(void) {}                                       //!< Re-establish the invariants this process's own moves can break.
        virtual void                                    updateRanges() = 0;                                         //!< Pull the ranges from the value.

        virtual double                                  computeLnProbabilityRanges(bool force = false);
        double                                          computeLnFossilRecord(size_t i) const;                      //!< Fossil-record log-term for taxon i.
        virtual double                                  observedRangeRate(size_t j) const { return 0.0; }           //!< A rate in interval j at which an observed range must have no event, beside lambda and mu.
        virtual void                                    appendIntervalRates(size_t j, std::vector<double> &v) const;  //!< Every rate of interval j a range term reads.
        void                                            storeRangeTerms(void);                                      //!< Keep the range terms and what they were computed under, for a rejected proposal.
        void                                            dirtyChangedIntervals(void);                                //!< Dirty the taxa born in or above the youngest interval whose rates changed.

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

        void                                            keepSpecialization(const DagNode *toucher);                  /* NOT VIRTUAL */
        void                                            restoreSpecialization(const DagNode *toucher);               /* NOT VIRTUAL */
        void                                            touchSpecialization(const DagNode *toucher, bool touchAll);  /* NOT VIRTUAL */

        std::vector<Taxon>                              taxa;                                                  //!< Taxa that will be attached to new simulated trees.
        std::string                                     condition;
        Reporting                                       record_reporting;                                       //!< Which occurrences the record reports.
        double                                          lnPresenceRecord(size_t i) const;             //!< Each reported interval holds at least its count.
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
        const TypedDagNode<double >*                    origin_age;                                            //!< The optional process start (NULL = the oldest birth, under no prior).

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


                                
        size_t                                          stored_range = 0;                                       //!< The range mvStratigraphicRange drew.
        double                                          stored_first = 0.0;                                     //!< The ages it replaced.
        double                                          stored_last  = 0.0;

        size_t                                          switched_range = 0;                                     //!< The range a switch flipped.
        bool                                            stored_survived = false;                                //!< The statuses before the switch.
        bool                                            stored_singleton = true;
        bool                                            switched = false;                                       //!< That undo is armed.

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
            double  death = 0.0, first = 0.0, last = 0.0, value = 0.0;
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
        mutable bool                                    warned_no_resample = false;
        mutable bool                                    warned_no_reporting = false;
        mutable bool                                    warned_inert_ntaxa = false;
        mutable bool                                    warned_no_survival = false;
        mutable bool                                    warned_no_singleton = false;
        mutable bool                                    warned_impossible = false;
        bool                                            resampled;                                              //!< The appearance undo is armed.
        double                                          present_time = 0.0;                                     //!< The present, where the timeline starts.
    };
}

#endif
