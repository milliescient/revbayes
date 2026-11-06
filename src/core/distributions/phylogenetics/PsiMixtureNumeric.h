#ifndef PsiMixtureNumeric_H
#define PsiMixtureNumeric_H

#include <cstddef>
#include <vector>

#include <boost/numeric/odeint/stepper/runge_kutta_fehlberg78.hpp>

namespace RevBayesCore {

    /**
     * @brief Extinction probabilities under a mixture of psi multiples drawn at birth, integrated with fixed eighth-order steps.
     *
     * The same quantities as PsiMixtureSolution, for any beta: with symmetric speciation both daughters draw a
     * component, which adds lambda beta pbar^2 to dp_h/dt and breaks the linearisation the closed form rests on.
     *
     *     dp_h/dt = mu - (lambda + mu + psi nu_h) p_h + lambda (1-beta) p_h pbar + lambda beta pbar^2
     *
     * The lambda integral is int lambda (1-beta) pbar, a budding whose daughter leaves no sample, which is
     * int lambda pbar at beta = 0. The stem matrix Q = dp(t)/dp(s) solves the variational equation.
     *
     * As the closed form caches each interval, this caches a grid: every interval start is a grid point, steps
     * are at most `step` long, and each point keeps p and the integral, each step its propagator dp/dp over it.
     * A change of arguments drops the grid from the youngest interval it touches. A query integrates only the
     * pieces between its ages and the nearest grid points, and multiplies the cached step propagators between.
     */
    class PsiMixtureNumeric {

    public:
        PsiMixtureNumeric(void) = default;

        void                    update(const std::vector<double> &times, const std::vector<double> &lambda, const std::vector<double> &mu,
                                       const std::vector<double> &psi, const std::vector<double> &rates, const std::vector<double> &weights,
                                       double rho, const std::vector<double> &beta);                         //!< As PsiMixtureSolution::update, with each interval's beta.
        size_t                  getNumberOfComponents(void) const { return num_components; }
        void                    at(double t, std::vector<double> &p, double &lambda_pbar) const;            //!< p of each component at age t, and int lambda (1-beta) pbar from the present up to t.
        double                  lambdaPbarIntegral(double t0, double t1) const;                             //!< int lambda (1-beta) pbar between ages t0 <= t1.
        void                    sensitivity(double s, double t, const std::vector<size_t> &columns,
                                            std::vector<double> &Q, std::vector<double> &Q_scale) const;    //!< Columns of Q = dp(t)/dp(s) for ages s <= t, each scaled by exp(Q_scale).

    private:
        size_t                  findIndex(double t) const;
        size_t                  gridBelow(double t) const;                                                  //!< The last grid point at or below age t, after extending the grid to it.
        void                    extendTo(double t) const;
        void                    advance(std::vector<double> &state, double from, double to, size_t columns) const;  //!< p, the integral and columns of Q, within one interval.
        double                  stepLength(size_t i) const;                                                 //!< The step in interval i, short against its fastest rate.
        void                    derivative(size_t i, const double *x, double *dx, size_t columns) const;
        void                    rk4(size_t i, std::vector<double> &x, double h, size_t columns) const;      //!< one fixed high-order step

        size_t                  num_components = 0;
        std::vector<double>     starts, lam, mu_, psi_, beta_, nu, w;
        double                  rho = 1.0;
        double                  step = 0.5;                    //!< the longest step
        double                  step_rate = 0.5;                //!< the longest step times the interval's fastest rate

        mutable std::vector< boost::numeric::odeint::runge_kutta_fehlberg78< std::vector<double> > > steppers;  //!< one per number of Q columns, so none resizes
        mutable std::vector<double>               grid;         //!< grid ages, the present first
        mutable std::vector<double>               grid_state;   //!< p and the integral at each grid point, H + 1 apiece
        mutable std::vector<double>               grid_step;    //!< dp(g_{k+1})/dp(g_k), row-major H x H, for each step
        mutable std::vector<double>               step_buffer, ext_state, at_state, sens_state, sens_last, sens_next, sens_M;  //!< scratch
    };

}

#endif
