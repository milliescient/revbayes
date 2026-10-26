#ifndef PsiMixtureSolution_H
#define PsiMixtureSolution_H

#include <cstddef>
#include <vector>

namespace RevBayesCore {

    /**
     * @brief The extinction probabilities of a birth-death process whose species each sample fossils
     * at their own multiple of psi, drawn at birth from a mixture.
     *
     * With p_h = y_h / z, the system y_h' = mu z - (lambda + mu + nu_h psi) y_h, z' = -lambda sum_h w_h y_h
     * is linear, and its generator, a diagonal bordered by one row and one column, has real eigenvalues
     * that are the roots of a secular equation. Each interval keeps its eigensystem, so p, int lambda pbar
     * and Q = dp(t)/dp(s) at any age cost a few exponentials, with no stepping. Q is a sum over pairs of
     * modes in which the pair's own term vanishes exactly, so it carries no cancellation however long
     * the interval.
     */
    class PsiMixtureSolution {

    public:
        PsiMixtureSolution(void) = default;

        void                    update(const std::vector<double> &times, const std::vector<double> &lambda, const std::vector<double> &mu,
                                       const std::vector<double> &psi, const std::vector<double> &rates, const std::vector<double> &weights,
                                       double rho);                                                         //!< Solve every interval for these rates, unless they are the last ones. times are the interval starts, the present first.
        size_t                  getNumberOfComponents(void) const { return num_components; }
        void                    at(double t, std::vector<double> &p, double &lambda_pbar) const;            //!< p of each component at age t, and int lambda pbar from the present up to t.
        double                  lambdaPbarIntegral(double t0, double t1) const;                             //!< int lambda pbar between ages t0 <= t1.
        void                    sensitivity(double s, double t, const std::vector<size_t> &columns,
                                            std::vector<double> &Q, std::vector<double> &Q_scale) const;    //!< Columns of Q = dp(t)/dp(s) for ages s <= t, each scaled by exp(Q_scale).

    private:

        //!< One interval's eigensystem: mode k decays at rate[k], with right vector (vy, vz), left vector (uy, uz)
        //!< already divided by u.v, and coefficient c[k] from p at the interval's start.
        struct Interval {
            double              start = 0.0;
            double              lambda_pbar_start = 0.0;
            std::vector<double> rate, vy, vz, uy, uz, c;
            size_t              top = 0, second = 0;            //!< the two slowest modes
            std::vector<double> full_Q;                         //!< Q over the whole interval, unless it is the oldest
            double              full_Q_scale = 0.0;
        };

        size_t                  findIndex(double t) const;
        void                    solve(Interval &I, double lambda, double mu, double psi, const std::vector<double> &p0) const;
        void                    addMode(Interval &I, double s, const std::vector<double> &vy, double vz, const std::vector<double> &uy, double uz) const;
        void                    pAt(const Interval &I, double dt, std::vector<double> &p, double &ln_z) const;
        void                    localSensitivity(const Interval &I, double from, double to, std::vector<double> &M, double &scale) const;

        size_t                  num_components = 0;
        std::vector<double>     nu;                             //!< the relative rate of each component
        std::vector<double>     w;                              //!< the weight of each component
        std::vector<Interval>   intervals;
        std::vector<double>     solved_for;                     //!< the arguments of the last update, which an unchanged call skips

        mutable std::vector<double> cw, e0, e1, col;            //!< scratch
        mutable std::vector<double> a, vy, uy, pole_a, pole_nu, pole_W;
        mutable std::vector<size_t> order, pole_start;
    };

}

#endif
