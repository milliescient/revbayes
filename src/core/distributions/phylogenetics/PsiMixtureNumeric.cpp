#include "PsiMixtureNumeric.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include <boost/numeric/odeint.hpp> // IWYU pragma: keep

using namespace RevBayesCore;


void PsiMixtureNumeric::update( const std::vector<double> &times, const std::vector<double> &lambda, const std::vector<double> &mu,
                                const std::vector<double> &psi, const std::vector<double> &rates, const std::vector<double> &weights,
                                double r, const std::vector<double> &beta )
{
    std::vector<double> b = beta.empty() ? std::vector<double>( times.size(), 0.0 ) : beta;

    // a change in what every interval shares drops the whole grid; otherwise only from the youngest changed interval
    bool all = ( num_components != rates.size() || starts != times || nu != rates || w != weights || rho != r || grid.empty() );
    size_t first = all ? 0 : times.size();
    for (size_t i = 0; first == times.size() && i < times.size(); i++)
    {
        if ( lam[i] != lambda[i] || mu_[i] != mu[i] || psi_[i] != psi[i] || beta_[i] != b[i] ) first = i;
    }

    num_components = rates.size();
    starts = times; lam = lambda; mu_ = mu; psi_ = psi; nu = rates; w = weights; rho = r; beta_ = b;

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


/** The step length in interval i: short against its fastest rate, so a fourth-order step is accurate. */
double PsiMixtureNumeric::stepLength( size_t i ) const
{
    double fastest = 0.0;
    for (size_t h = 0; h < num_components; h++) fastest = std::max( fastest, lam[i] + mu_[i] + psi_[i]*nu[h] );
    double h_max = std::min( step, step_rate / std::max( fastest, 1E-12 ) );

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
    double l = lam[i], m = mu_[i], f = psi_[i], b = beta_[i];

    double pbar = 0.0;
    for (size_t h = 0; h < H; h++) pbar += w[h] * x[h];
    for (size_t h = 0; h < H; h++)
    {
        dx[h] = m - ( l + m + f*nu[h] ) * x[h] + l * (1.0 - b) * x[h] * pbar + l * b * pbar * pbar;
    }
    dx[H] = l * (1.0 - b) * pbar;

    // the variational equation: J_ak = delta_ak (-(l+m+f nu_a) + l(1-b) pbar) + (l(1-b) p_a + 2 l b pbar) w_k
    for (size_t c = 0; c < C; c++)
    {
        const double *q = &x[H + 1 + c*H];
        double *dq = &dx[H + 1 + c*H];
        double wq = 0.0;
        for (size_t k = 0; k < H; k++) wq += w[k] * q[k];
        for (size_t a = 0; a < H; a++)
        {
            dq[a] = ( -( l + m + f*nu[a] ) + l * (1.0 - b) * pbar ) * q[a] + ( l * (1.0 - b) * x[a] + 2.0 * l * b * pbar ) * wq;
        }
    }
}


/** One fixed eighth-order (Runge-Kutta-Fehlberg 7(8)) step of length h in interval i, in place, by a stepper kept per state size. */
void PsiMixtureNumeric::rk4( size_t i, std::vector<double> &x, double h, size_t C ) const
{
    if ( steppers.size() <= C ) steppers.resize( C + 1 );
    auto system = [&]( const std::vector<double> &y, std::vector<double> &dy, double /*t*/ ) { derivative( i, y.data(), dy.data(), C ); };
    steppers[C].do_step( system, x, 0.0, h );
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
