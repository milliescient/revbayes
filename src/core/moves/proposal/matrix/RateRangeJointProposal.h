#ifndef RateRangeJointProposal_H
#define RateRangeJointProposal_H

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <utility>
#include <vector>

#include "MatrixReal.h"
#include "Proposal.h"
#include "StochasticNode.h"

namespace RevBayesCore {

    class AbstractFossilizedBirthDeathRangeProcess;

    /**
     * A rate move that carries the range ends with it.
     *
     * Scales a rate, or flips a shift indicator, and redraws every range end that can fall in an
     * interval where that rate changed: on the extinction side each extinction time on
     * (present, tauK] from mu(d) exp(-int_d^tauK (mu + psi)), on the origination side each origination
     * time on (tau1, origin] from lambda(b) exp(-int_tau1^b (mu + psi)). Those densities omit the
     * process's other terms, so they are only proposals, corrected by the Hastings ratio. The oldest
     * origination is welded to the origin and never redrawn.
     *
     * With rjSwitch, x is a dnReversibleJumpMixture jump that is 0 for no shift: a jump switched on
     * is proposed lognormally about the rate its cells already have, and one switched off returns
     * them to the rate they inherit.
     */
    class RateRangeJointProposal : public Proposal {

    public:
        RateRangeJointProposal( StochasticNode<double> *r, StochasticNode<MatrixReal> *m, double l, bool rj = false, double sd = 0.5, bool births = false );
        RateRangeJointProposal( StochasticNode<std::int64_t> *s, StochasticNode<MatrixReal> *m, bool births = false );

        void                                    cleanProposal(void);
        RateRangeJointProposal*                 clone(void) const;
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
        // the target rate (lambda or mu) and the hazard rates, read through the DAG
        void                                    readRates( AbstractFossilizedBirthDeathRangeProcess *process, std::vector<double> &t, std::vector<double> &w, std::vector<double> &mu, std::vector<double> &psi ) const;
        // set the rate, read the rates, and restore, leaving the new value in place
        void                                    probeRate( double value, AbstractFossilizedBirthDeathRangeProcess *process, std::vector<double> &t, std::vector<double> &w, std::vector<double> &mu, std::vector<double> &psi );
        // log density of a range end under the proposal, and a draw from it, on (lo, hi]
        double                                  lnDensity( double x, double lo, double hi, const std::vector<double> &t, const std::vector<double> &w, const std::vector<double> &mu, const std::vector<double> &psi ) const;
        double                                  draw( double lo, double hi, const std::vector<double> &t, const std::vector<double> &w, const std::vector<double> &mu, const std::vector<double> &psi ) const;
        double                                  lnLognormal( double x, double centre ) const;

        StochasticNode<double>*                 rate;
        StochasticNode<std::int64_t>*           shift;
        StochasticNode<MatrixReal>*             ranges;
        double                                  lambda;
        bool                                    rj_switch;
        double                                  sigma;
        bool                                    births;                 // origination side: lambda and column 3, not mu and column 0

        double                                  stored_rate;
        std::int64_t                            stored_shift;
        std::vector< std::pair<size_t,double> > stored_ends;
    };

}

#endif
