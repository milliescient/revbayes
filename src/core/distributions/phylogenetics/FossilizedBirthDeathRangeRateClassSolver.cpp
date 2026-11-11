#include <algorithm>
#include <array>
#include <cmath>

#include "FossilizedBirthDeathRangeRateClassSolver.h"

using namespace RevBayesCore;


void FossilizedBirthDeathRangeRateClassSolver::setRates(const std::vector<double> &times, const std::vector<double> &lambda, const std::vector<double> &mu,
                                const std::vector<double> &psi, const std::vector<double> &rates, const std::vector<double> &weights,
                                double rho, const std::vector<double> &beta, const std::vector<double> &anagenetic)
{
    std::vector<double> b = beta.empty() ? std::vector<double>( times.size(), 0.0 ) : beta;
    std::vector<double> a = anagenetic.empty() ? std::vector<double>( times.size(), 0.0 ) : anagenetic;

    // the cached stem switches hold until the rates change
    if ( times == r_times && lambda == r_lambda && mu == r_mu && psi == r_psi && rates == r_nu && weights == r_w
         && rho == r_rho && b == r_beta && a == r_ana ) return;

    r_times = times; r_lambda = lambda; r_mu = mu; r_psi = psi; r_nu = rates; r_w = weights; r_rho = rho;
    r_beta = b; r_ana = a;
    switch_cache.clear();
}


size_t FossilizedBirthDeathRangeRateClassSolver::rateIndex(double t) const
{
    size_t i = 0;
    while ( i + 1 < r_times.size() && r_times[i+1] <= t ) ++i;
    return i;
}


/**
 * By forward-backward over the components, from source at top down to end_component at end (all components if it is -1).
 */
void FossilizedBirthDeathRangeRateClassSolver::stemSwitches(double end, double top, const std::vector<double> &source, long end_component,
                                    double lo, double hi, double out[3]) const
{
    std::vector<double> key = { end, top, double( end_component ), lo, hi };
    key.insert( key.end(), source.begin(), source.end() );
    std::map<std::vector<double>, std::array<double,3> >::const_iterator hit = switch_cache.find( key );
    if ( hit != switch_cache.end() )
    {
        for (size_t m = 0; m < 3; m++) out[m] = hit->second[m];
        return;
    }

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

        // quadrature nodes y and weights wq: Simpson on the solver's grid, with cells graded next to the stem's ends
        double rate = r_lambda[i] + r_mu[i] + r_psi[i] * nu_max + r_ana[i];
        std::vector<double> y, wq;
        auto add = [&]( double t, double w )
        {
            if ( y.empty() == false && y.back() == t ) wq.back() += w;
            else { y.push_back( t ); wq.push_back( w ); }
        };
        auto cell = [&]( double u, double v )
        {
            double w = ( v - u ) / 6.0;
            add( u, w ); add( 0.5 * ( u + v ), 4.0 * w ); add( v, w );
        };
        // cells from u to v, graded from u, or from v when from_u is false
        auto graded = [&]( double u, double v, bool from_u, double h_max )
        {
            std::vector<double> x( 1, 0.0 );
            double len = v - u, h = std::min( h_max, 0.02 / rate );
            while ( x.back() + h < len && h < h_max ) { x.push_back( x.back() + h ); h *= 2.0; }
            double rest = len - x.back();
            size_t n = std::max( size_t( 1 ), size_t( std::ceil( rest / h_max ) ) );
            double last = x.back();
            for (size_t j = 1; j <= n; j++) x.push_back( last + rest * j / n );
            if ( from_u ) for (size_t j = 0; j + 1 < x.size(); j++) cell( u + x[j], u + x[j+1] );
            else for (size_t j = x.size() - 1; j > 0; j--) cell( v - x[j], v - x[j-1] );
        };

        std::vector<double> G;
        gridAges( t0, t_end, G );
        if ( G.size() < 3 )
        {
            double half = 0.5 * ( t_end - t0 );
            graded( t0, t0 + half, true, half / ( panels / 2 ) );
            graded( t0 + half, t_end, false, half / ( panels / 2 ) );
        }
        else
        {
            double g = ( G.back() - G.front() ) / ( G.size() - 1 );
            if ( t0 == end ) graded( t0, G.front(), true, g );
            else cell( t0, G.front() );

            // composite Simpson on the equal cells, the last three by Simpson's 3/8 when their number is odd
            size_t n = G.size() - 1, even = ( n % 2 == 0 ? n : n - 3 );
            for (size_t j = 0; j < even; j += 2)
            {
                add( G[j], g / 3.0 ); add( G[j+1], 4.0 * g / 3.0 ); add( G[j+2], g / 3.0 );
            }
            if ( even < n )
            {
                double w = 3.0 * g / 8.0;
                add( G[n-3], w ); add( G[n-2], 3.0 * w ); add( G[n-1], 3.0 * w ); add( G[n], w );
            }

            if ( t_end == top ) graded( G.back(), t_end, false, g );
            else cell( G.back(), t_end );
        }
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

        // the integrand by mode, normalized by a . g at each point
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
        for (size_t j = 0; j < N; j++) for (size_t m = 0; m < 3; m++) out[m] += wq[j] * f[j][m];
        t0 = t_end;
    }

    if ( switch_cache.size() > 100000 ) switch_cache.clear();
    switch_cache[key] = { out[0], out[1], out[2] };
}
