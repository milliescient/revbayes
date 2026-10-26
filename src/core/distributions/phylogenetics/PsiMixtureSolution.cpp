#include "PsiMixtureSolution.h"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace RevBayesCore;


namespace {

    const size_t NONE = size_t(-1);

    /**
     * The secular function f(s) = s + lambda mu sum_q W_q / (s + a_q) in the coordinate x = s + a_ref,
     * or x = s when ref is NONE. Measuring from the nearest pole keeps s + a_q exact to its last digit
     * however close two poles lie.
     */
    struct Secular {
        double lm, psi;
        const std::vector<double> *a, *nu, *W;
        size_t ref;

        double off( void ) const { return ref == NONE ? 0.0 : (*a)[ref]; }
        double gap( size_t q ) const { return ref == NONE ? (*a)[q] : psi * ( (*nu)[q] - (*nu)[ref] ); }

        void eval( double x, double &f, double &df ) const
        {
            f  = x - off();
            df = 1.0;
            for (size_t q = 0; q < W->size(); q++)
            {
                double y = 1.0 / ( x + gap(q) );
                double r = (*W)[q] * y;
                f  += lm * r;
                df -= lm * r * y;
            }
        }

        void evalDerivative( double x, double &df, double &d2f ) const
        {
            df  = 1.0;
            d2f = 0.0;
            for (size_t q = 0; q < W->size(); q++)
            {
                double y = 1.0 / ( x + gap(q) );
                double r = lm * (*W)[q] * y * y;
                df  -= r;
                d2f += 2.0 * r * y;
            }
        }
    };


    /** The root of g in (lo, hi), where g changes sign once: Newton's steps, bisection when one leaves the bracket. */
    template <class G>
    double safeNewton( const G &g, double lo, double hi, bool positive_at_lo )
    {
        double x = 0.5 * (lo + hi);
        for (size_t iter = 0; iter < 200; iter++)
        {
            double v, dv;
            g( x, v, dv );
            if ( v == 0.0 ) return x;

            if ( (v > 0.0) == positive_at_lo ) lo = x; else hi = x;

            double next = x - v / dv;
            if ( !( next > lo && next < hi ) ) next = 0.5 * (lo + hi);

            if ( fabs( next - x ) <= 4.0 * std::numeric_limits<double>::epsilon() * fabs( next ) || next == lo || next == hi ) return next;
            x = next;
        }
        return x;
    }

}


size_t PsiMixtureSolution::findIndex( double t ) const
{
    size_t i = 0;
    while ( i + 1 < intervals.size() && intervals[i+1].start <= t ) i++;
    return i;
}


/** Store a mode, dividing its left vector by u.v so that the modes sum to the identity. */
void PsiMixtureSolution::addMode( Interval &I, double s, const std::vector<double> &vy, double vz, const std::vector<double> &uy, double uz ) const
{
    double n = vz * uz;
    for (size_t h = 0; h < num_components; h++) n += vy[h] * uy[h];

    I.rate.push_back( s );
    I.vy.insert( I.vy.end(), vy.begin(), vy.end() );
    I.vz.push_back( vz );
    for (size_t h = 0; h < num_components; h++) I.uy.push_back( uy[h] / n );
    I.uz.push_back( uz / n );
}


/** The eigensystem of one interval's generator, and the coefficients of p0 in it. */
void PsiMixtureSolution::solve( Interval &I, double lambda, double mu, double psi, const std::vector<double> &p0 ) const
{
    size_t H = num_components;

    I.rate.clear(); I.vy.clear(); I.vz.clear(); I.uy.clear(); I.uz.clear();

    a.resize(H); vy.resize(H); uy.resize(H);
    for (size_t h = 0; h < H; h++) a[h] = lambda + mu + nu[h] * psi;

    if ( lambda == 0.0 || mu == 0.0 )
    {
        // the border loses a side, so each component decays on its own beside one mode at rate 0
        for (size_t h = 0; h < H; h++)
        {
            vy[h] = ( lambda == 0.0 && a[h] > 0.0 ) ? mu / a[h] : 0.0;
            uy[h] = ( lambda == 0.0 ) ? 0.0 : -lambda * w[h] / a[h];
        }
        addMode( I, 0.0, vy, 1.0, uy, 1.0 );

        for (size_t h = 0; h < H; h++)
        {
            std::fill( vy.begin(), vy.end(), 0.0 );
            std::fill( uy.begin(), uy.end(), 0.0 );
            vy[h] = 1.0;
            uy[h] = 1.0;
            double vz = ( lambda == 0.0 ) ? 0.0 : lambda * w[h] / a[h];
            double uz = ( lambda == 0.0 && a[h] > 0.0 ) ? -mu / a[h] : 0.0;
            addMode( I, -a[h], vy, vz, uy, uz );
        }
    }
    else
    {
        // the distinct poles, ascending; components with equal a share one, and order[pole_start[q]...] are its members
        order.resize(H);
        for (size_t h = 0; h < H; h++) order[h] = h;
        std::sort( order.begin(), order.end(), [&](size_t x, size_t y){ return a[x] < a[y] || ( a[x] == a[y] && x < y ); } );

        pole_a.clear(); pole_nu.clear(); pole_W.clear(); pole_start.clear();
        for (size_t k = 0; k < H; k++)
        {
            size_t h = order[k];
            if ( pole_a.empty() || a[h] != pole_a.back() )
            {
                pole_start.push_back( k );
                pole_nu.push_back( nu[h] );
                pole_a.push_back( a[h] );
                pole_W.push_back( 0.0 );
            }
            pole_W.back() += w[h];
        }
        size_t num_poles = pole_a.size();
        pole_start.push_back( H );

        Secular f;
        f.lm = lambda * mu; f.psi = psi; f.a = &pole_a; f.nu = &pole_nu; f.W = &pole_W;
        auto in = [&](size_t ref) { f.ref = ref; return [&f](double x, double &v, double &dv){ f.eval(x, v, dv); }; };

        size_t num_bordered = 1;
        for (size_t q = 0; q < num_poles; q++) num_bordered += ( pole_W[q] > 0.0 );

        // a root is kept as x in the coordinate of ref, so that s + a_q = x + gap(q) is exact
        std::vector<std::pair<size_t,double> > roots;
        roots.reserve( num_bordered );

        // one root between each pair of neighbouring poles with weight, found from the nearer one
        size_t r_prev = NONE;
        for (size_t q = 0; q < num_poles; q++)
        {
            if ( pole_W[q] == 0.0 ) continue;
            if ( r_prev != NONE )
            {
                f.ref = r_prev;
                double gap = f.gap(q);
                double v, dv;
                f.eval( -0.5 * gap, v, dv );
                size_t ref = ( v > 0.0 ) ? r_prev : q;
                double lo  = ( v > 0.0 ) ? -0.5 * gap : 0.0;
                double hi  = ( v > 0.0 ) ? 0.0 : 0.5 * gap;
                roots.push_back( std::make_pair( ref, safeNewton( in(ref), lo, hi, true ) ) );
            }
            r_prev = q;
        }

        // two roots above the lowest pole with weight, either side of the minimum of the convex f there
        size_t q0 = 0;
        while ( pole_W[q0] == 0.0 ) q0++;
        double a0 = pole_a[q0];
        f.ref = q0;
        double m = safeNewton( [&f](double x, double &df, double &d2f){ f.evalDerivative(x, df, d2f); }, 0.0, a0, false );
        double fm, dfm;
        f.eval( m, fm, dfm );
        if ( fm < 0.0 )
        {
            roots.push_back( std::make_pair( q0, safeNewton( in(q0), 0.0, m, true ) ) );

            if ( m < 0.5 * a0 )
            {
                roots.push_back( std::make_pair( q0, safeNewton( in(q0), m, a0, false ) ) );
            }
            else
            {
                roots.push_back( std::make_pair( NONE, safeNewton( in(NONE), m - a0, 0.0, false ) ) );
            }
        }
        else
        {
            // a double root to rounding: lambda = mu with nearly all weight at nu psi = 0. Splitting it
            // costs about sqrt(eps) in p, where the stepped solution would lose as much
            double h = 1E-8 * m;
            roots.push_back( std::make_pair( q0, m - h ) );
            roots.push_back( std::make_pair( q0, m + h ) );
        }

        // the bordered modes: v = ( mu / (s + a), 1 ), u = ( -lambda w / (s + a), 1 )
        for (const std::pair<size_t,double> &root : roots)
        {
            f.ref = root.first;
            double x = root.second;
            for (size_t q = 0; q < num_poles; q++)
            {
                double den = x + f.gap(q);
                for (size_t j = pole_start[q]; j < pole_start[q+1]; j++)
                {
                    size_t h = order[j];
                    vy[h] = mu / den;
                    uy[h] = -lambda * w[h] / den;
                }
            }
            addMode( I, x - f.off(), vy, 1.0, uy, 1.0 );
        }

        for (size_t q = 0; q < num_poles; q++)
        {
            size_t first = pole_start[q], last = pole_start[q+1];

            if ( pole_W[q] > 0.0 )
            {
                // components sharing a pole also differ by modes that leave z alone: v = e_i - (w_i / w_anchor) e_anchor,
                // u = e_i - w / W over the group
                size_t anchor = order[first];
                for (size_t j = first; j < last; j++) if ( w[order[j]] > w[anchor] ) anchor = order[j];

                for (size_t j = first; j < last; j++)
                {
                    size_t i = order[j];
                    if ( i == anchor ) continue;
                    std::fill( vy.begin(), vy.end(), 0.0 );
                    std::fill( uy.begin(), uy.end(), 0.0 );
                    vy[i] = 1.0;
                    vy[anchor] = -w[i] / w[anchor];
                    for (size_t g = first; g < last; g++) uy[order[g]] = -w[order[g]] / pole_W[q];
                    uy[i] += 1.0;
                    addMode( I, -pole_a[q], vy, 0.0, uy, 0.0 );
                }
            }
            else
            {
                // a component without weight never feeds z, so it decays on its own: v = e_i, and u is
                // what the bordered modes leave of e_i
                for (size_t j = first; j < last; j++)
                {
                    size_t i = order[j];
                    std::fill( vy.begin(), vy.end(), 0.0 );
                    vy[i] = 1.0;
                    for (size_t h = 0; h < H; h++) uy[h] = ( h == i ) ? 1.0 : 0.0;
                    double uz = 0.0;
                    for (size_t b = 0; b < num_bordered; b++)
                    {
                        double v = I.vy[b*H + i];
                        for (size_t h = 0; h < H; h++) uy[h] -= v * I.uy[b*H + h];
                        uz -= v * I.uz[b];
                    }
                    addMode( I, -pole_a[q], vy, 0.0, uy, uz );
                }
            }
        }
    }

    size_t K = I.rate.size();

    I.top = 0;
    for (size_t k = 1; k < K; k++) if ( I.rate[k] > I.rate[I.top] ) I.top = k;
    I.second = ( I.top == 0 ) ? 1 : 0;
    for (size_t k = 0; k < K; k++) if ( k != I.top && I.rate[k] > I.rate[I.second] ) I.second = k;

    I.c.resize(K);
    for (size_t k = 0; k < K; k++)
    {
        double x = I.uz[k];
        for (size_t h = 0; h < H; h++) x += I.uy[k*H + h] * p0[h];
        I.c[k] = x;
    }
}


/** p a time dt above the interval's start, and ln z, with z = 1 at the start. */
void PsiMixtureSolution::pAt( const Interval &I, double dt, std::vector<double> &p, double &ln_z ) const
{
    size_t H = num_components, K = I.rate.size();

    p.assign( H, 0.0 );
    double Z = 0.0;
    for (size_t k = 0; k < K; k++)
    {
        double x = I.c[k] * exp( ( I.rate[k] - I.rate[I.top] ) * dt );
        Z += x * I.vz[k];
        for (size_t h = 0; h < H; h++) p[h] += x * I.vy[k*H + h];
    }
    for (size_t h = 0; h < H; h++) p[h] /= Z;

    ln_z = I.rate[I.top] * dt + log( Z );
}


/**
 * M[a*H + b] = dp_a(to)/dp_b(from) / exp(scale), for offsets from <= to within one interval.
 * Q_ab = sum_{k != j} e^{(s_k + s_j) dt} c_j (v_ka v_jz - v_ja v_kz) u_kb / z^2, whose largest pair
 * is the two slowest modes, so scale starts at (s_second - s_top) dt and every other factor is at most 1.
 */
void PsiMixtureSolution::localSensitivity( const Interval &I, double from, double to, std::vector<double> &M, double &scale ) const
{
    size_t H = num_components, K = I.rate.size();
    double top = I.rate[I.top], second = I.rate[I.second];

    // the coefficients at from, where z = 1
    cw.assign( K, 0.0 );
    double z_from = 0.0;
    for (size_t j = 0; j < K; j++)
    {
        cw[j] = I.c[j] * exp( ( I.rate[j] - top ) * from );
        z_from += cw[j] * I.vz[j];
    }
    for (size_t j = 0; j < K; j++) cw[j] /= z_from;

    double dt = to - from;
    e0.assign( K, 0.0 );
    e1.assign( K, 0.0 );
    double Z = 0.0;
    for (size_t k = 0; k < K; k++)
    {
        e0[k] = exp( ( I.rate[k] - top ) * dt );
        e1[k] = ( k == I.top ) ? 0.0 : exp( ( I.rate[k] - second ) * dt );
        Z += cw[k] * e0[k] * I.vz[k];
    }

    // G[k*H + a] = sum_{j != k} factor c_j (v_ka v_jz - v_ja v_kz), then M_ab = sum_k G_ka u_kb
    col.assign( K*H, 0.0 );
    for (size_t k = 0; k < K; k++)
    {
        for (size_t j = 0; j < K; j++)
        {
            if ( j == k ) continue;
            double factor = ( k == I.top ) ? e1[j] : ( ( j == I.top ) ? e1[k] : e0[k] * e1[j] );
            double x = factor * cw[j];
            if ( x == 0.0 ) continue;
            for (size_t a = 0; a < H; a++) col[k*H + a] += x * ( I.vy[k*H + a] * I.vz[j] - I.vy[j*H + a] * I.vz[k] );
        }
    }

    M.assign( H*H, 0.0 );
    double largest = 0.0;
    for (size_t a = 0; a < H; a++)
    {
        for (size_t b = 0; b < H; b++)
        {
            double x = 0.0;
            for (size_t k = 0; k < K; k++) x += col[k*H + a] * I.uy[k*H + b];
            // Q is nonnegative, so a negative entry is rounding
            x = std::max( 0.0, x / (Z * Z) );
            M[a*H + b] = x;
            largest = std::max( largest, x );
        }
    }

    scale = ( second - top ) * dt;
    if ( largest > 0.0 )
    {
        for (size_t x = 0; x < H*H; x++) M[x] /= largest;
        scale += log( largest );
    }
}


void PsiMixtureSolution::update( const std::vector<double> &times, const std::vector<double> &lambda, const std::vector<double> &mu,
                                 const std::vector<double> &psi, const std::vector<double> &rates, const std::vector<double> &weights,
                                 double rho )
{
    std::vector<double> &args = cw;
    args.assign( times.begin(), times.end() );
    for (const std::vector<double> *v : { &lambda, &mu, &psi, &rates, &weights }) args.insert( args.end(), v->begin(), v->end() );
    args.push_back( rho );
    args.push_back( double( rates.size() ) );
    if ( args == solved_for ) return;
    solved_for = args;

    num_components = rates.size();
    nu = rates;
    w  = weights;

    size_t N = times.size();
    intervals.resize( N );

    std::vector<double> p( num_components, 1.0 - rho );
    double lambda_pbar = 0.0;

    for (size_t i = 0; i < N; i++)
    {
        Interval &I = intervals[i];
        I.start = times[i];
        I.lambda_pbar_start = lambda_pbar;
        solve( I, lambda[i], mu[i], psi[i], p );

        if ( i + 1 < N )
        {
            double dt = times[i+1] - times[i];
            double ln_z;
            pAt( I, dt, p, ln_z );
            lambda_pbar -= ln_z;
            localSensitivity( I, 0.0, dt, I.full_Q, I.full_Q_scale );
        }
    }
}


void PsiMixtureSolution::at( double t, std::vector<double> &p, double &lambda_pbar ) const
{
    const Interval &I = intervals[ findIndex(t) ];

    double ln_z;
    pAt( I, t - I.start, p, ln_z );
    lambda_pbar = I.lambda_pbar_start - ln_z;
}


double PsiMixtureSolution::lambdaPbarIntegral( double t0, double t1 ) const
{
    std::vector<double> p;
    double l0, l1;
    at( t0, p, l0 );
    at( t1, p, l1 );
    return l1 - l0;
}


/** Columns of Q = dp(t)/dp(s), one per component in columns, Q[c*H + a] = dp_a(t)/dp_columns[c](s) / exp(Q_scale[c]). */
void PsiMixtureSolution::sensitivity( double s, double t, const std::vector<size_t> &columns, std::vector<double> &Q, std::vector<double> &Q_scale ) const
{
    size_t H = num_components, C = columns.size();
    size_t i = findIndex(s), j = findIndex(t);

    std::vector<double> M;
    double scale;
    const Interval &first = intervals[i];
    localSensitivity( first, s - first.start, ( i == j ? t : intervals[i+1].start ) - first.start, M, scale );

    Q.assign( H*C, 0.0 );
    Q_scale.assign( C, scale );
    for (size_t c = 0; c < C; c++)
    {
        for (size_t a = 0; a < H; a++) Q[c*H + a] = M[a*H + columns[c]];
    }

    std::vector<double> next(H);
    auto apply = [&](const std::vector<double> &A, double A_scale) {
        for (size_t c = 0; c < C; c++)
        {
            double *Qc = &Q[c*H];
            double largest = 0.0;
            for (size_t a = 0; a < H; a++)
            {
                double x = 0.0;
                for (size_t e = 0; e < H; e++) x += A[a*H + e] * Qc[e];
                next[a] = x;
                largest = std::max( largest, x );
            }
            if ( largest > 0.0 )
            {
                for (size_t a = 0; a < H; a++) Qc[a] = next[a] / largest;
                Q_scale[c] += A_scale + log( largest );
            }
            else
            {
                std::fill( Qc, Qc + H, 0.0 );
            }
        }
    };

    for (size_t k = i + 1; k < j; k++) apply( intervals[k].full_Q, intervals[k].full_Q_scale );

    if ( j > i )
    {
        const Interval &last = intervals[j];
        localSensitivity( last, 0.0, t - last.start, M, scale );
        apply( M, scale );
    }
}
