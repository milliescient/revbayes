#include <algorithm>
#include <array>
#include <cmath>

#include "PsiMixtureSolver.h"

using namespace RevBayesCore;


void PsiMixtureSolver::setRates(const std::vector<double> &times, const std::vector<double> &lambda, const std::vector<double> &mu,
                                const std::vector<double> &psi, const std::vector<double> &rates, const std::vector<double> &weights,
                                const std::vector<double> &beta, const std::vector<double> &anagenetic)
{
    r_times = times; r_lambda = lambda; r_mu = mu; r_psi = psi; r_nu = rates; r_w = weights;
    r_beta = beta.empty() ? std::vector<double>( times.size(), 0.0 ) : beta;
    r_ana  = anagenetic.empty() ? std::vector<double>( times.size(), 0.0 ) : anagenetic;
}


size_t PsiMixtureSolver::rateIndex(double t) const
{
    size_t i = 0;
    while ( i + 1 < r_times.size() && r_times[i+1] <= t ) ++i;
    return i;
}


/**
 * A stem from age top down to end can pass through unsampled ancestors, each in a fresh component. By forward-backward
 * over the components, with a(t) from source at top and g(t) to end_component at end (summed over all where it is -1),
 * E[N] = int a(t) r(t) (w . g(t)) dt / (a . g). A budding hands the lineage on from a parent left in its own component k,
 * at lambda (1-beta) p_k; a symmetric speciation's other daughter is fresh, at 2 lambda beta pbar; anagenesis at its rate.
 */
void PsiMixtureSolver::stemSwitches(double end, double top, const std::vector<double> &source, long end_component,
                                    double lo, double hi, double out[3]) const
{
    size_t H = getNumberOfComponents();
    std::vector<size_t> all( H );
    for (size_t h = 0; h < H; h++) all[h] = h;
    out[0] = out[1] = out[2] = 0.0;

    double nu_max = *std::max_element( r_nu.begin(), r_nu.end() );
    std::vector<double> end_vec( H, end_component >= 0 ? 0.0 : 1.0 );
    if ( end_component >= 0 ) end_vec[ end_component ] = 1.0;

    // M[c*H + e]: from component e at age t down to c at age s <= t, up to a common factor
    auto stepMatrix = [&]( double s, double t, std::vector<double> &m )
    {
        sensitivity( s, t, all, Q, Q_scale );
        double mx = *std::max_element( Q_scale.begin(), Q_scale.end() );
        m.resize( H*H );
        for (size_t c = 0; c < H; c++) for (size_t e = 0; e < H; e++) m[c*H + e] = Q[c*H + e] * exp( Q_scale[c] - mx );
    };
    auto normalize = []( std::vector<double> &v )
    {
        double mx = *std::max_element( v.begin(), v.end() );
        if ( mx > 0.0 ) for (double &x : v) x /= mx;
    };

    const size_t panels = 32;
    double t0 = lo;
    for (size_t i = rateIndex(t0); t0 < hi; i++)
    {
        double t_end = ( i + 1 < r_times.size() ) ? std::min( hi, r_times[i+1] ) : hi;
        double len = t_end - t0;

        // a stiff component's transient at either end of the stem: a mesh graded from 0.02 of its decay time
        double rate = r_lambda[i] + r_mu[i] + r_psi[i] * nu_max + r_ana[i];
        double h_max = len / panels, h = std::min( h_max, 0.02 / rate );
        std::vector<double> left( 1, 0.0 );
        while ( left.back() + h < 0.5 * len && h < h_max ) { left.push_back( left.back() + h ); h *= 2.0; }
        std::vector<double> x = left;
        double rest = len - 2.0 * left.back();
        size_t n_mid = std::max( size_t( 1 ), size_t( std::ceil( rest / h_max ) ) );
        for (size_t j = 1; j <= n_mid; j++) x.push_back( left.back() + rest * j / n_mid );
        for (size_t j = left.size() - 1; j-- > 0; ) x.push_back( len - left[j] );

        // the ages, youngest first, with each cell's midpoint for Simpson
        std::vector<double> y;
        for (size_t j = 0; j + 1 < x.size(); j++) { y.push_back( t0 + x[j] ); y.push_back( t0 + 0.5 * ( x[j] + x[j+1] ) ); }
        y.push_back( t_end );
        size_t N = y.size();

        steps.resize( N - 1 );
        for (size_t j = 0; j + 1 < N; j++) stepMatrix( y[j], y[j+1], steps[j] );

        // a from the top down, then step by step to the youngest age
        fwd.assign( N, std::vector<double>( H, 0.0 ) );
        if ( y[N-1] < top )
        {
            stepMatrix( y[N-1], top, M );
            for (size_t c = 0; c < H; c++) for (size_t e = 0; e < H; e++) fwd[N-1][c] += source[e] * M[c*H + e];
        }
        else fwd[N-1] = source;
        normalize( fwd[N-1] );
        for (size_t j = N - 1; j-- > 0; )
        {
            for (size_t c = 0; c < H; c++) for (size_t e = 0; e < H; e++) fwd[j][c] += fwd[j+1][e] * steps[j][c*H + e];
            normalize( fwd[j] );
        }

        // g from the end up, then step by step to the oldest age
        bwd.assign( N, std::vector<double>( H, 0.0 ) );
        if ( y[0] > end )
        {
            stepMatrix( end, y[0], M );
            for (size_t e = 0; e < H; e++) for (size_t c = 0; c < H; c++) bwd[0][e] += M[c*H + e] * end_vec[c];
        }
        else bwd[0] = end_vec;
        normalize( bwd[0] );
        for (size_t j = 1; j < N; j++)
        {
            for (size_t e = 0; e < H; e++) for (size_t c = 0; c < H; c++) bwd[j][e] += steps[j-1][c*H + e] * bwd[j-1][c];
            normalize( bwd[j] );
        }

        // the integrand by mode; a . g normalizes each point, so the vectors' scales cancel
        double lp;
        std::vector<std::array<double,3> > f( N );
        for (size_t j = 0; j < N; j++)
        {
            at( y[j], p, lp );
            double pbar = 0.0, ag = 0.0, a_sum = 0.0, a_p = 0.0, wg = 0.0;
            for (size_t c = 0; c < H; c++)
            {
                pbar  += r_w[c] * p[c];
                ag    += fwd[j][c] * bwd[j][c];
                a_sum += fwd[j][c];
                a_p   += fwd[j][c] * p[c];
                wg    += r_w[c] * bwd[j][c];
            }
            double scale = ( ag > 0.0 ? wg / ag : 0.0 );
            f[j] = { scale * r_lambda[i] * ( 1.0 - r_beta[i] ) * a_p,
                     scale * 2.0 * r_lambda[i] * r_beta[i] * pbar * a_sum,
                     scale * r_ana[i] * a_sum };
        }
        for (size_t j = 0; j + 2 < N; j += 2)
        {
            double wgt = ( y[j+2] - y[j] ) / 6.0;
            for (size_t m = 0; m < 3; m++) out[m] += wgt * ( f[j][m] + 4.0 * f[j+1][m] + f[j+2][m] );
        }
        t0 = t_end;
    }
}
