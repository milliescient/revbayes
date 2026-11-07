#include "PsiMixtureNumeric.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>

#include <boost/numeric/odeint.hpp> // IWYU pragma: keep

using namespace RevBayesCore;


void PsiMixtureNumeric::update( const std::vector<double> &times, const std::vector<double> &lambda, const std::vector<double> &mu,
                                const std::vector<double> &psi, const std::vector<double> &rates, const std::vector<double> &weights,
                                double r, const std::vector<double> &beta, const std::vector<double> &anagenetic )
{
    setRates( times, lambda, mu, psi, rates, weights, beta, anagenetic );
    std::vector<double> b = beta.empty() ? std::vector<double>( times.size(), 0.0 ) : beta;
    std::vector<double> a = anagenetic.empty() ? std::vector<double>( times.size(), 0.0 ) : anagenetic;

    // a change in what every interval shares drops the whole grid; otherwise only from the youngest changed interval
    bool all = ( num_components != rates.size() || starts != times || nu != rates || w != weights || rho != r || grid.empty() );
    size_t first = all ? 0 : times.size();
    for (size_t i = 0; first == times.size() && i < times.size(); i++)
    {
        if ( lam[i] != lambda[i] || mu_[i] != mu[i] || psi_[i] != psi[i] || beta_[i] != b[i] || ana_[i] != a[i] ) first = i;
    }

    num_components = rates.size();
    starts = times; lam = lambda; mu_ = mu; psi_ = psi; nu = rates; w = weights; rho = r; beta_ = b; ana_ = a;
    etd_cache.clear(); etd_interval = size_t(-1);

    size_t H = num_components;

    if ( first == times.size() ) return;

    if ( first == 0 )
    {
        grid.assign( 1, times.front() );
        grid_state.assign( H + 1, 1.0 - rho );
        grid_state[H] = 0.0;
        grid_step.clear();
        return;
    }

    // keep the points up to the changed interval's start, which only younger intervals reach
    size_t keep = std::upper_bound( grid.begin(), grid.end(), starts[first] ) - grid.begin();
    grid.resize( keep );
    grid_state.resize( keep * (H + 1) );
    grid_step.resize( (keep - 1) * H * H );
}


size_t PsiMixtureNumeric::findIndex( double t ) const
{
    size_t i = std::upper_bound( starts.begin(), starts.end(), t ) - starts.begin();
    return i == 0 ? 0 : i - 1;
}


bool PsiMixtureNumeric::isStiff( size_t i ) const
{
    for (size_t h = 0; h < num_components; h++) if ( ( lam[i] + mu_[i] + psi_[i]*nu[h] + ana_[i] ) * step > stiff_limit ) return true;
    return false;
}


/** The step length in interval i: short against its fastest rate, or in a stiff interval against its fastest slow one. */
double PsiMixtureNumeric::stepLength( size_t i ) const
{
    bool stiff = isStiff( i );
    double fastest = 0.0;
    for (size_t h = 0; h < num_components; h++)
    {
        double k = lam[i] + mu_[i] + psi_[i]*nu[h] + ana_[i];
        if ( stiff == false || k * step <= stiff_limit ) fastest = std::max( fastest, k );
    }
    if ( stiff ) fastest = std::max( fastest, lam[i] + mu_[i] + ana_[i] );
    double h_max = std::min( step, ( stiff ? etd_step_rate : step_rate ) / std::max( fastest, 1E-12 ) );

    // equal steps between the interval's ends; the oldest runs on without end, in steps of h_max
    if ( i + 1 < starts.size() )
    {
        double len = starts[i+1] - starts[i];
        return len / std::ceil( len / h_max );
    }
    return h_max;
}


/** The derivative of p (H), the lambda integral, then C columns of Q (H each), in interval i. */
void PsiMixtureNumeric::derivative( size_t i, const double *x, double *dx, size_t C ) const
{
    size_t H = num_components;
    double l = lam[i], m = mu_[i], f = psi_[i], b = beta_[i], a = ana_[i];

    double pbar = 0.0;
    for (size_t h = 0; h < H; h++) pbar += w[h] * x[h];
    for (size_t h = 0; h < H; h++)
    {
        dx[h] = m - ( l + m + f*nu[h] + a ) * x[h] + l * (1.0 - b) * x[h] * pbar + l * b * pbar * pbar + a * pbar;
    }
    dx[H] = l * (1.0 - b) * pbar;

    // the variational equation: J_ak = delta_ak (-(l+m+f nu_a+a) + l(1-b) pbar) + (l(1-b) p_a + 2 l b pbar + a) w_k
    for (size_t c = 0; c < C; c++)
    {
        const double *q = &x[H + 1 + c*H];
        double *dq = &dx[H + 1 + c*H];
        double wq = 0.0;
        for (size_t k = 0; k < H; k++) wq += w[k] * q[k];
        for (size_t k = 0; k < H; k++)
        {
            dq[k] = ( -( l + m + f*nu[k] + a ) + l * (1.0 - b) * pbar ) * q[k] + ( l * (1.0 - b) * x[k] + 2.0 * l * b * pbar + a ) * wq;
        }
    }
}


/** One fixed eighth-order (Runge-Kutta-Fehlberg 7(8)) step of length h in interval i, in place, by a stepper kept per state size. */
void PsiMixtureNumeric::rk4( size_t i, std::vector<double> &x, double h, size_t C ) const
{
    if ( isStiff( i ) ) { etd4( i, x, h, C ); return; }
    if ( steppers.size() <= C ) steppers.resize( C + 1 );
    auto system = [&]( const std::vector<double> &y, std::vector<double> &dy, double /*t*/ ) { derivative( i, y.data(), dy.data(), C ); };
    steppers[C].do_step( system, x, 0.0, h );
}


/**
 * Cox and Matthews' ETDRK4 with the diagonal linear part L = -(lambda + mu + psi nu_h + a) of each p_h and each Q row,
 * so a fast-sampled component's decay is solved exactly and the step need only resolve the slow rates. The
 * coefficients are Kassam and Trefethen's contour means, accurate as hL goes to zero.
 */
void PsiMixtureNumeric::etd4( size_t i, std::vector<double> &x, double h, size_t C ) const
{
    // a fast transient (each Q column starts from the identity) has local error ~ exp(-k t) (k s)^5, so the substeps
    // start at eps^(1/5)/k and grow as it decays, at most doubling
    double k_max = 0.0, k_slow = std::numeric_limits<double>::infinity();
    for (size_t c = 0; c < num_components; c++)
    {
        double k = lam[i] + mu_[i] + psi_[i]*nu[c] + ana_[i];
        k_max = std::max( k_max, k );
        if ( k * step > stiff_limit ) k_slow = std::min( k_slow, k );
    }
    double sub = 0.0, done = 0.0;
    while ( done < h )
    {
        double allowed = 0.025 / k_max * std::exp( std::min( k_slow * done / 5.0, 600.0 ) );
        double next = std::min( { sub > 0.0 ? 2.0 * sub : allowed, allowed, h - done } );
        if ( h - done - next < 0.5 * next ) next = h - done;
        etdStep( i, x, next, C );
        done += next;
        sub = next;
    }
}


void PsiMixtureNumeric::etdStep( size_t i, std::vector<double> &x, double h, size_t C ) const
{
    size_t H = num_components, n = H + 1 + C*H;

    if ( etd_interval != i || etd_cache.size() > 4096 ) { etd_cache.clear(); etd_interval = i; }
    EtdCoefficients &etd = etd_cache[h];
    if ( etd.E.empty() )
    {
        etd.E.assign( H + 1, 0.0 ); etd.E2 = etd.Q = etd.f1 = etd.f2 = etd.f3 = etd.E;
        const size_t M = 32;
        for (size_t h_ = 0; h_ <= H; h_++)
        {
            double L = ( h_ < H ) ? -( lam[i] + mu_[i] + psi_[i]*nu[h_] + ana_[i] ) : 0.0;
            double z = L * h;
            etd.E[h_] = std::exp( z ); etd.E2[h_] = std::exp( z / 2.0 );
            std::complex<double> q = 0.0, a = 0.0, b = 0.0, c = 0.0;
            for (size_t j = 0; j < M; j++)
            {
                std::complex<double> r = z + std::polar( 1.0, M_PI * ( j + 0.5 ) / M * 2.0 );
                std::complex<double> e = std::exp( r ), r3 = r*r*r;
                q += ( std::exp( r / 2.0 ) - 1.0 ) / r;
                a += ( -4.0 - r + e * ( 4.0 - 3.0*r + r*r ) ) / r3;
                b += ( 2.0 + r + e * ( -2.0 + r ) ) / r3;
                c += ( -4.0 - 3.0*r - r*r + e * ( 4.0 - r ) ) / r3;
            }
            etd.Q[h_] = h * q.real() / M;
            etd.f1[h_] = h * a.real() / M; etd.f2[h_] = h * b.real() / M; etd.f3[h_] = h * c.real() / M;
        }
    }

    // which coefficient each entry takes: p_h and Q row h take component h, the integral the last
    auto slot = [H]( size_t k ) { return k <= H ? k : ( k - H - 1 ) % H; };
    etd_L.resize( n );
    for (size_t k = 0; k < n; k++)
    {
        size_t s = slot( k );
        etd_L[k] = ( s < H ) ? -( lam[i] + mu_[i] + psi_[i]*nu[s] + ana_[i] ) : 0.0;
    }
    auto nonlinear = [&]( const std::vector<double> &y, std::vector<double> &N ) {
        N.resize( n );
        derivative( i, y.data(), N.data(), C );
        for (size_t k = 0; k < n; k++) N[k] -= etd_L[k] * y[k];
    };

    nonlinear( x, etd_N );
    etd_a.resize( n ); etd_b.resize( n ); etd_c.resize( n );
    for (size_t k = 0; k < n; k++) etd_a[k] = etd.E2[slot(k)] * x[k] + etd.Q[slot(k)] * etd_N[k];
    nonlinear( etd_a, etd_Na );
    for (size_t k = 0; k < n; k++) etd_b[k] = etd.E2[slot(k)] * x[k] + etd.Q[slot(k)] * etd_Na[k];
    nonlinear( etd_b, etd_Nb );
    for (size_t k = 0; k < n; k++) etd_c[k] = etd.E2[slot(k)] * etd_a[k] + etd.Q[slot(k)] * ( 2.0 * etd_Nb[k] - etd_N[k] );
    nonlinear( etd_c, etd_Nc );
    for (size_t k = 0; k < n; k++)
    {
        size_t s = slot( k );
        x[k] = etd.E[s] * x[k] + etd.f1[s] * etd_N[k] + 2.0 * etd.f2[s] * ( etd_Na[k] + etd_Nb[k] ) + etd.f3[s] * etd_Nc[k];
    }
}


/**
 * Advance the state from one age to an older one within one interval, in steps no longer than the interval's,
 * all equal, so the result is a smooth function of the ages.
 */
void PsiMixtureNumeric::advance( std::vector<double> &state, double from, double to, size_t C ) const
{
    if ( to <= from ) return;
    size_t i = findIndex( from );
    double h_max = stepLength( i );
    size_t n = size_t( std::ceil( ( to - from ) / h_max - 1E-9 ) );
    if ( n == 0 ) n = 1;
    double h = ( to - from ) / n;
    for (size_t k = 0; k < n; k++) rk4( i, state, h, C );
}


/** Extend the grid out to age t, one interval step at a time, each step's propagator kept. */
void PsiMixtureNumeric::extendTo( double t ) const
{
    size_t H = num_components;
    std::vector<double> &s = ext_state;
    s.resize( H + 1 + H*H );
    while ( grid.back() < t )
    {
        double from = grid.back();
        size_t i = findIndex( from );
        double to = from + stepLength( i );
        if ( i + 1 < starts.size() ) to = std::min( to, starts[i+1] );

        std::fill( s.begin(), s.end(), 0.0 );
        std::copy( grid_state.end() - (H + 1), grid_state.end(), s.begin() );
        for (size_t c = 0; c < H; c++) s[H + 1 + c*H + c] = 1.0;
        rk4( i, s, to - from, H );

        grid.push_back( to );
        grid_state.insert( grid_state.end(), s.begin(), s.begin() + H + 1 );

        // row-major M[a*H + e] = dp_a(to)/dp_e(from), from the columns e
        size_t base = grid_step.size();
        grid_step.resize( base + H*H );
        for (size_t e = 0; e < H; e++) for (size_t a = 0; a < H; a++) grid_step[base + a*H + e] = s[H + 1 + e*H + a];
    }
}


size_t PsiMixtureNumeric::gridBelow( double t ) const
{
    extendTo( t );
    return std::upper_bound( grid.begin(), grid.end(), t ) - grid.begin() - 1;
}


void PsiMixtureNumeric::at( double t, std::vector<double> &p, double &lambda_pbar ) const
{
    size_t k = gridBelow( t );
    size_t H = num_components;
    std::vector<double> &s = at_state;
    s.assign( grid_state.begin() + k*(H + 1), grid_state.begin() + (k + 1)*(H + 1) );
    advance( s, grid[k], t, 0 );

    p.assign( s.begin(), s.begin() + num_components );
    lambda_pbar = s[num_components];
}


double PsiMixtureNumeric::lambdaPbarIntegral( double t0, double t1 ) const
{
    std::vector<double> p;
    double l0, l1;
    at( t0, p, l0 );
    at( t1, p, l1 );
    return l1 - l0;
}


void PsiMixtureNumeric::sensitivity( double s, double t, const std::vector<size_t> &columns, std::vector<double> &Q, std::vector<double> &Q_scale ) const
{
    size_t H = num_components, C = columns.size();

    size_t kt = gridBelow( t );
    size_t ks = gridBelow( s );

    // from s up to the next grid point, or straight to t inside one step
    // p alone up to s, then p with the columns
    std::vector<double> &head = at_state;
    head.assign( grid_state.begin() + ks*(H + 1), grid_state.begin() + (ks + 1)*(H + 1) );
    advance( head, grid[ks], s, 0 );
    std::vector<double> &state = sens_state;
    state.assign( H + 1 + C*H, 0.0 );
    std::copy( head.begin(), head.end(), state.begin() );
    for (size_t c = 0; c < C; c++) state[H + 1 + c*H + columns[c]] = 1.0;
    double first_to = ( kt == ks ) ? t : grid[ks + 1];
    advance( state, s, first_to, C );

    Q.assign( state.begin() + H + 1, state.end() );
    Q_scale.assign( C, 0.0 );

    std::vector<double> &next = sens_next;
    next.resize( H );
    auto rescale = [&]() {
        for (size_t c = 0; c < C; c++)
        {
            double *q = &Q[c*H];
            double largest = 0.0;
            for (size_t a = 0; a < H; a++) largest = std::max( largest, std::fabs( q[a] ) );
            if ( largest > 0.0 )
            {
                for (size_t a = 0; a < H; a++) q[a] /= largest;
                Q_scale[c] += std::log( largest );
            }
        }
    };
    // apply a row-major H x H propagator to every column
    auto apply = [&]( const double *M ) {
        for (size_t c = 0; c < C; c++)
        {
            double *q = &Q[c*H];
            for (size_t a = 0; a < H; a++)
            {
                double x = 0.0;
                for (size_t e = 0; e < H; e++) x += M[a*H + e] * q[e];
                next[a] = x;
            }
            std::copy( next.begin(), next.end(), q );
        }
    };
    if ( kt == ks )
    {
        rescale();
        return;
    }

    // the cached whole steps between, rescaled now and then
    for (size_t k = ks + 1; k < kt; k++)
    {
        apply( &grid_step[k*H*H] );
        if ( ( k - ks ) % 8 == 0 ) rescale();
    }

    // and from the last grid point to t, carrying the columns themselves
    if ( t > grid[kt] )
    {
        std::vector<double> &last = sens_last;
        last.assign( H + 1 + C*H, 0.0 );
        std::copy( grid_state.begin() + kt*(H + 1), grid_state.begin() + (kt + 1)*(H + 1), last.begin() );
        std::copy( Q.begin(), Q.end(), last.begin() + H + 1 );
        advance( last, grid[kt], t, C );
        std::copy( last.begin() + H + 1, last.end(), Q.begin() );
    }
    rescale();
}
