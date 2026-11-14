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
#include "Tree.h"

namespace RevBayesCore {

    class AbstractFossilizedBirthDeathRangeProcess;

    /**
     * An extinction-rate move that carries the extinction times with it.
     *
     * Scales an extinction rate, or flips a shift indicator, and redraws every extinction time on
     * (present, tauK] that can fall in an interval where the rate changed, from
     * mu(d) exp(-int_d^tauK (mu + psi)). That density omits the process's other terms, so it is only a
     * proposal, corrected by the Hastings ratio. The ranges are a range matrix or an extended tree.
     *
     * With rjSwitch, x is a dnReversibleJumpMixture jump that is 0 for no shift: a jump switched on
     * is proposed lognormally about the rate its cells already have, and one switched off returns
     * them to the rate they inherit.
     */
    class RateRangeJointProposal : public Proposal {

    public:
        RateRangeJointProposal( StochasticNode<double> *r, StochasticNode<MatrixReal> *m, StochasticNode<Tree> *t, double l, bool rj = false, double sd = 0.5 );
        RateRangeJointProposal( StochasticNode<std::int64_t> *s, StochasticNode<MatrixReal> *m, StochasticNode<Tree> *t );

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
        // the range matrix, or an extended tree whose tip ages are the extinction times
        void                                    addValueNode(void);
        AbstractFossilizedBirthDeathRangeProcess* process(void) const;
        double                                  extinctionTime( size_t i ) const;
        void                                    setExtinctionTime( size_t i, double d );
        // the extinction and sampling rates, read through the DAG
        void                                    readRates( AbstractFossilizedBirthDeathRangeProcess *process, std::vector<double> &t, std::vector<double> &mu, std::vector<double> &psi ) const;
        // set the rate, read the rates, and restore, leaving the new value in place
        void                                    probeRate( double value, AbstractFossilizedBirthDeathRangeProcess *process, std::vector<double> &t, std::vector<double> &mu, std::vector<double> &psi );
        // log density of an extinction time under the proposal, and a draw from it, on (lo, hi]
        double                                  lnDensity( double x, double lo, double hi, const std::vector<double> &t, const std::vector<double> &mu, const std::vector<double> &psi ) const;
        double                                  draw( double lo, double hi, const std::vector<double> &t, const std::vector<double> &mu, const std::vector<double> &psi ) const;
        double                                  lnLognormal( double x, double centre ) const;

        StochasticNode<double>*                 rate;
        StochasticNode<std::int64_t>*           shift;
        StochasticNode<MatrixReal>*             ranges;
        StochasticNode<Tree>*                   tree;
        double                                  lambda;
        bool                                    rj_switch;
        double                                  sigma;

        double                                  stored_rate;
        std::int64_t                            stored_shift;
        std::vector< std::pair<size_t,double> > stored_ends;
    };

}

#endif
