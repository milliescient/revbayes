#ifndef FossilizedBirthDeathRangeRateClassClosedForm_H
#define FossilizedBirthDeathRangeRateClassClosedForm_H

#include <cstddef>
#include <vector>

#include "FossilizedBirthDeathRangeRateClassSolver.h"

namespace RevBayesCore {

    /**
     * @brief Extinction probabilities under rate classes that scale psi, in closed form.
     *
     * With p_h = y_h / z the system is linear, so each interval keeps its eigensystem and p,
     * int lambda pbar and Q = dp(t)/dp(s) are closed forms.
     */
    class FossilizedBirthDeathRangeRateClassClosedForm : public FossilizedBirthDeathRangeRateClassSolver {

    public:
        FossilizedBirthDeathRangeRateClassClosedForm(void) = default;

        void                    update(const std::vector<double> &times, const std::vector<double> &lambda, const std::vector<double> &mu,
                                       const std::vector<double> &psi, const std::vector<double> &rates, const std::vector<double> &weights,
                                       double rho);                                                         //!< Solve the intervals whose rates changed and those older. times are the interval starts, the present first.
        size_t                  getNumberOfComponents(void) const override { return num_components; }
        void                    at(double t, std::vector<double> &p, double &lambda_pbar) const override;            //!< p of each component at age t, and int lambda pbar from the present up to t.
        double                  lambdaPbarIntegral(double t0, double t1) const override;                             //!< int lambda pbar between ages t0 <= t1.
        void                    sensitivity(double s, double t, const std::vector<size_t> &columns,
                                            std::vector<double> &Q, std::vector<double> &Q_scale) const override;    //!< Columns of Q = dp(t)/dp(s) for ages s <= t, each scaled by exp(Q_scale).

    private:

        //!< One interval's eigensystem: mode k decays at rate[k], with right vector (vy, vz), left vector (uy, uz)
        //!< already divided by u.v, and coefficient c[k] from p at the interval's start.
        struct Interval {
            double              start = 0.0;
            double              lambda_pbar_start = 0.0;
            std::vector<double> p_start;                        //!< p at the start, from the younger intervals
            std::vector<double> key;                            //!< its start, end and rates, as last solved
            std::vector<double> rate, vy, vz, uy, uz, c;
            size_t              top = 0, second = 0;            //!< the two slowest modes
            std::vector<double> pair;                           //!< pair[(k*K + j)*H + a] = v_ka v_jz - v_ja v_kz
            std::vector<double> full_Q;                         //!< Q over the whole interval, unless it is the oldest
            double              full_Q_scale = 0.0;
        };

        size_t                  findIndex(double t) const;
        void                    solve(Interval &I, double lambda, double mu, double psi, const std::vector<double> &p0) const;
        void                    addMode(Interval &I, double s, const std::vector<double> &vy, double vz, const std::vector<double> &uy, double uz) const;
        void                    pAt(const Interval &I, double dt, std::vector<double> &p, double &ln_z) const;
        void                    wholeInterval(Interval &I, double dt) const;                                //!< Q over a whole interval, on one scale.
        void                    localApply(const Interval &I, double from, double to, const std::vector<size_t> *unit, const double *in,
                                           size_t C, double *out, double *scale) const;     //!< dp(to)/dp(from) applied to C columns, each scaled.

        size_t                  num_components = 0;
        std::vector<double>     nu;                             //!< the relative rate of each component
        std::vector<double>     w;                              //!< the weight of each component
        std::vector<Interval>   intervals;
        std::vector<double>     solved_shared;                  //!< the arguments every interval shares, as last solved

        mutable std::vector<double> cw, e0, e1, col, proj, prev, step_scale, next;  //!< scratch
        mutable std::vector<double> a, vy, uy, pole_a, pole_nu, pole_W;
        mutable std::vector<size_t> order, pole_start;
    };

}

#endif
