#ifndef FossilizedBirthDeathRangeRateClassSolver_H
#define FossilizedBirthDeathRangeRateClassSolver_H

#include <array>
#include <cstddef>
#include <map>
#include <vector>

namespace RevBayesCore {

    /**
     * @brief Base for the psi-mixture solvers: p, the lambda integral, the stem matrix Q, and the hidden species changes on a stem.
     */
    class FossilizedBirthDeathRangeRateClassSolver {

    public:
        virtual                ~FossilizedBirthDeathRangeRateClassSolver(void) {}

        virtual size_t          getNumberOfComponents(void) const = 0;
        virtual void            at(double t, std::vector<double> &p, double &lambda_pbar) const = 0;            //!< p of each component at age t, and the lambda integral from the present up to t.
        virtual void            sensitivity(double s, double t, const std::vector<size_t> &columns,
                                            std::vector<double> &Q, std::vector<double> &Q_scale) const = 0;    //!< Columns of Q = dp(t)/dp(s) for ages s <= t, each scaled by exp(Q_scale).

        virtual void            gridAges(double lo, double hi, std::vector<double> &ages) const { ages.clear(); }  //!< The solver's own grid points strictly between lo and hi, where it caches its steps.

        void                    stemSwitches(double end, double top, const std::vector<double> &source, long end_component,
                                             double lo, double hi, double out[3]) const;                       //!< Expected hidden species changes on [lo, hi] of a stem: budding, symmetric, anagenetic.

    protected:
        void                    setRates(const std::vector<double> &times, const std::vector<double> &lambda, const std::vector<double> &mu,
                                         const std::vector<double> &psi, const std::vector<double> &rates, const std::vector<double> &weights,
                                         double rho, const std::vector<double> &beta, const std::vector<double> &anagenetic);  //!< The arguments stemSwitches reads. An empty beta or anagenetic is zero.

    private:
        size_t                  rateIndex(double t) const;

        std::vector<double>     r_times, r_lambda, r_mu, r_psi, r_nu, r_w, r_beta, r_ana;
        double                  r_rho = -1.0;
        mutable std::map<std::vector<double>, std::array<double,3> > switch_cache;    //!< stemSwitches by its arguments, until the rates change
        mutable std::vector<double>              Q, Q_scale, p, M;
        mutable std::vector<std::vector<double> > steps, fwd, bwd;
    };

}

#endif
