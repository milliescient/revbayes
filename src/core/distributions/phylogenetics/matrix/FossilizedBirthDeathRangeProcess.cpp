#include "FossilizedBirthDeathRangeProcess.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <set>
#include <string>
#include <vector>

#include "AbstractFossilizedBirthDeathRangeProcess.h"
#include "DistributionExponential.h"
#include "MatrixReal.h"
#include "RbMathCombinatorialFunctions.h"
#include "RbMathLogic.h"
#include "RbMathFunctions.h"
#include "RbMathLogic.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbException.h"
#include "StochasticNode.h"
#include "TypedDistribution.h"
#include "RbVector.h"
#include "RbVectorImpl.h"
#include "Taxon.h"
#include "TimeInterval.h"

namespace RevBayesCore { class DagNode; }
namespace RevBayesCore { template <class valueType> class TypedDagNode; }

using namespace RevBayesCore;


/**
 * Constructor. 
 * We delegate most parameters to the base class and initialize the members.
 *
 * \param[in]    s              Speciation rates.
 * \param[in]    e              Extinction rates.
 * \param[in]    p              Fossil sampling rates.
 * \param[in]    c              Fossil observation counts.
 * \param[in]    r              Instantaneous sampling probabilities.
 * \param[in]    t              Rate change times.
 * \param[in]    cdt            Condition of the process (time/sampling/survival).
 * \param[in]    tn             Taxa.
 * \param[in]    c              Complete sampling?
 */
FossilizedBirthDeathRangeProcess::FossilizedBirthDeathRangeProcess(const DagNode *inspeciation,
                                                                     const DagNode *inextinction,
                                                                     const DagNode *inpsi,
                                                                     const TypedDagNode<double> *inrho,
                                                                     const TypedDagNode< RbVector<double> > *intimes,
                                                                     const std::string &incondition,
                                                                     const std::vector<Taxon> &intaxa,
                                                                     bool complete_record,
                                                                     bool coexistence,
                                                                     bool require_parent,
                                                                     const TypedDagNode<double> *inorigin,
                                                                     double inpresent,
                                                                     const TypedDagNode<RbVector<double> > *inrates,
                                                                     const TypedDagNode<Simplex> *inweights) :
    TypedDistribution<MatrixReal>(new MatrixReal(intaxa.size(), 4)),
    AbstractFossilizedBirthDeathRangeProcess(inspeciation, inextinction, inpsi, inrho, intimes, incondition, intaxa, complete_record, inorigin, inpresent, inrates, inweights),
    coexistence( coexistence ),
    require_parent( require_parent )
{

    dirty_gamma = std::vector<bool>(taxa.size(), true);
    gamma_i     = std::vector<size_t>(taxa.size(), 0);
    gamma_links = std::vector<std::vector<bool> >(taxa.size(), std::vector<bool>(taxa.size(), false));

    for(std::vector<const DagNode*>::iterator it = range_parameters.begin(); it != range_parameters.end(); it++)
    {
        addParameter(*it);
    }

    redrawValue();
    updateGamma(true);
}


/**
 * The clone function is a convenience function to create proper copies of inherited objected.
 * E.g. a.clone() will create a clone of the correct type even if 'a' is of derived type 'B'.
 *
 * \return A new copy of myself 
 */

void FossilizedBirthDeathRangeProcess::setMcmcMode(bool tf)
{
    TypedDistribution<MatrixReal>::setMcmcMode(tf);
    // the appearances live in the value, where the element moves reach them
    if ( tf == true ) warnIfNoReportingNode();
    if ( tf == true ) warnIfInertNTaxaCondition();
    if ( tf == true ) warnIfNoSurvivalMove();
    if ( tf == true ) warnIfNoSingletonMove();
}

FossilizedBirthDeathRangeProcess* FossilizedBirthDeathRangeProcess::clone( void ) const
{
    return new FossilizedBirthDeathRangeProcess( *this );
}


/** Set the matrix value, clipping any appearance the new ranges leave out of support. */
void FossilizedBirthDeathRangeProcess::setValue(MatrixReal *v, bool force)
{
    TypedDistribution<MatrixReal>::setValue(v, force);

    for (size_t i = 0; i < taxa.size(); i++)
    {
        double d  = (*this->value)[i][0];
        double b  = (*this->value)[i][3];

        ranges[i].last  = (*this->value)[i][1];
        ranges[i].first = (*this->value)[i][2];

        // a value from outside carries no status
        ranges[i].survived = ( d == getPresent() );

        // oldest age: valid range [max(first_min,d), min(max_age,b))
        double lo = std::max( ranges[i].first_min, d );
        double hi = std::min( taxa[i].getMaxAge(), b );
        if ( hi > lo && ( ranges[i].first < lo || ranges[i].first >= b ) )
        {
            ranges[i].first = 0.5 * ( lo + hi );
        }

        // youngest age: valid range [max(d,min_age), min(first,last_max)]
        double lo_y = std::max( d, taxa[i].getMinAge() );
        double hi_y = std::min( ranges[i].first, ranges[i].last_max );
        if ( hi_y > lo_y && ( ranges[i].last < lo_y || ranges[i].last > hi_y ) )
        {
            ranges[i].last = 0.5 * ( lo_y + hi_y );
        }

        (*this->value)[i][1] = ranges[i].last;
        (*this->value)[i][2] = ranges[i].first;
    }

    // a node already touched keeps its cache, which no longer matches
    dirty_taxa  = std::vector<bool>(taxa.size(), true);
    dirty_gamma = std::vector<bool>(taxa.size(), true);

    repairRanges();

    updateRanges();
}


/** The statuses and, with more than one, the psi mixture components are in no value, so they ride with the ranges. */
std::string FossilizedBirthDeathRangeProcess::getHiddenStateString( void ) const
{
    json j;
    j["ranges"]    = this->value->toJSON();
    j["singleton"] = json::array();
    j["survived"]  = json::array();

    for (size_t i = 0; i < taxa.size(); ++i)
    {
        j["singleton"].push_back( ranges[i].singleton );
        j["survived"].push_back( ranges[i].survived );
    }

    if ( getNumberOfPsiMixtureComponents() > 1 )
    {
        j["component"] = json::array();
        for (size_t i = 0; i < taxa.size(); ++i) j["component"].push_back( getPsiMixtureComponent( i ) );
    }

    return j.dump();
}


void FossilizedBirthDeathRangeProcess::setHiddenStateFromString( const std::string &s )
{
    json j = json::parse( s );

    MatrixReal *m = this->value->clone();
    m->initFromString( j["ranges"].dump() );

    // the singleton status ties the appearances when the value is set
    for (size_t i = 0; i < taxa.size(); ++i) ranges[i].singleton = j["singleton"][i].get<bool>();

    setValue( m );

    for (size_t i = 0; i < taxa.size(); ++i) ranges[i].survived = j["survived"][i].get<bool>();

    if ( getNumberOfPsiMixtureComponents() > 1 )
    {
        if ( j.contains("component") == false ) throw RbException("dnFBDRP: the checkpoint holds no psi mixture components.");

        for (size_t i = 0; i < taxa.size(); ++i)
        {
            size_t k = j["component"][i].get<size_t>();
            if ( k >= getNumberOfPsiMixtureComponents() ) throw RbException("dnFBDRP: a psi mixture component in the checkpoint is out of range.");
            ranges[i].component = k;
        }
    }
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
double FossilizedBirthDeathRangeProcess::computeLnProbability( void )
{
    double lnProb = computeLnProbabilityRanges();

    // without them, a birth carries no factor for how many lineages could have produced it,
    // which is exact only where the standing diversity does not vary
    if ( coexistence == false && require_parent == false ) return lnProb;

    updateGamma();

    for( size_t i = 0; i < taxa.size(); i++ )
    {
        if ( gamma_i[i] == 0 )
        {
            // only the origin buds from nothing
            if ( i == max_birth ) continue;

            // Where ranges were left out, the lineage a birth came from may be one of them, so
            // none alive among those kept is an observation rather than an impossibility. The
            // weight still counts the attachments there are, and log(0) is not one of them.
            if ( require_parent == false ) continue;

            return RbConstants::Double::neginf;
        }

        // multiply by the number of possible birth locations
        if ( coexistence == true ) lnProb += log( gamma_i[i] );
    }

    return lnProb;
}




/**
 * Compute the number of ranges that intersect with range i
 *
 * \param[in]    i      index of range for which to compute gamma
 *
 * \return Small gamma
 */
void FossilizedBirthDeathRangeProcess::updateGamma(bool force)
{
    for (size_t i = 0; i < taxa.size(); i++)
    {
        if ( dirty_gamma[i] || force )
        {
            double bi = (*this->value)[i][3];
            double di = (*this->value)[i][0];

            if ( force == true ) gamma_i[i] = 0;

            for (size_t j = 0; j < taxa.size(); j++)
            {
                if (i == j) continue;

                double bj = (*this->value)[j][3];
                double dj = (*this->value)[j][0];

                bool linki = ( bi < bj && bi > dj );
                bool linkj = ( bj < bi && bj > di );

                if ( gamma_links[i][j] != linki && force == false )
                {
                    gamma_i[i] += linki ? 1 : -1;
                }
                if ( gamma_links[j][i] != linkj && force == false )
                {
                    gamma_i[j] += linkj ? 1 : -1;
                }

                if ( force == true ) gamma_i[i] += linki;

                gamma_links[i][j] = linki;
                gamma_links[j][i] = linkj;
            }
        }
    }
}


/**
 * Compute the log-transformed probability of the current value under the current parameter values.
 *
 */
void FossilizedBirthDeathRangeProcess::updateRanges( void )
{
    max_birth = 0;

    for (size_t i = 0; i < taxa.size(); i++)
    {
        ranges[i].divergence = (*this->value)[i][3];
        ranges[i].death = (*this->value)[i][0];
        ranges[i].first = (*this->value)[i][2];
        ranges[i].last  = (*this->value)[i][1];

        if ( ranges[i].divergence > ranges[max_birth].divergence ) max_birth = i;
    }

    if ( origin_age == NULL )
    {
        origin = ranges[max_birth].divergence;
        return;
    }

    origin = origin_age->getValue();

    // the oldest birth follows the origin
    if ( weld_origin == true && ranges[max_birth].divergence != origin )
    {
        if ( touched == true )
        {
            stored_repairs.push_back( std::make_pair( max_birth*4 + 3, (*this->value)[max_birth][3] ) );
        }

        (*this->value)[max_birth][3] = origin;
        ranges[max_birth].divergence      = origin;

        dirty_gamma[max_birth] = true;
        dirty_taxa[max_birth]  = true;
    }
}


/**
 * Simulate new speciation times.
 */
void FossilizedBirthDeathRangeProcess::redrawValue(void)
{
    drawRanges();

    for (size_t i = 0; i < taxa.size(); i++)
    {
        (*this->value)[i][0] = ranges[i].death;
        (*this->value)[i][1] = ranges[i].last;
        (*this->value)[i][2] = ranges[i].first;
        (*this->value)[i][3] = ranges[i].divergence;
    }

    repairRanges();

    updateRanges();
}


/** Below two occurrences tau_1 and tau_K are one age, so the column a move did not write follows the one it did. */
void FossilizedBirthDeathRangeProcess::repairRanges( const std::set<size_t> &touched )
{
    for ( std::set<size_t>::const_iterator it = touched.begin(); it != touched.end(); it++ )
    {
        size_t i = (*it) / 4;
        size_t c = (*it) % 4;

        if ( ranges[i].singleton == false ) continue;

        size_t sib;
        if      ( c == 1 ) sib = 2;
        else if ( c == 2 ) sib = 1;
        else continue;

        stored_repairs.push_back( std::make_pair( i*4 + sib, (*this->value)[i][sib] ) );
        (*this->value)[i][sib] = (*this->value)[i][c];
    }
}


void FossilizedBirthDeathRangeProcess::repairRanges( void )
{
    for (size_t i = 0; i < taxa.size(); i++)
    {
        if ( ranges[i].singleton ) (*this->value)[i][1] = (*this->value)[i][2];
    }
}


void FossilizedBirthDeathRangeProcess::keepSpecialization(const DagNode *toucher)
{
    stored_repairs.clear();

    dirty_gamma = std::vector<bool>(taxa.size(), false);

    AbstractFossilizedBirthDeathRangeProcess::keepSpecialization(toucher);
}

void FossilizedBirthDeathRangeProcess::restoreSpecialization(const DagNode *toucher)
{
    // newest first
    for ( std::vector<std::pair<size_t,double> >::reverse_iterator it = stored_repairs.rbegin(); it != stored_repairs.rend(); ++it )
    {
        (*this->value)[ it->first / 4 ][ it->first % 4 ] = it->second;
    }
    stored_repairs.clear();

    AbstractFossilizedBirthDeathRangeProcess::restoreSpecialization(toucher);
}


void FossilizedBirthDeathRangeProcess::touchSpecialization(const DagNode *toucher, bool touchAll)
{
    weld_origin = ( toucher != dag_node );

    if ( toucher == dag_node )
    {
        if ( touched == false )
        {
            stored_likelihood = partial_likelihood;
            stored_ranges     = ranges;

            std::set<size_t> touched_indices = dag_node->getTouchedElementIndices();

            for ( std::set<size_t>::iterator it = touched_indices.begin(); it != touched_indices.end(); it++)
            {
                size_t i = (*it) / 4; // N x 4 row-major, so the linear index over the columns is the taxon

                dirty_gamma[i] = true;
                dirty_taxa[i]  = true;

            }

            repairRanges( touched_indices );
        }

        touched = true;
    }
    else
    {
        AbstractFossilizedBirthDeathRangeProcess::touchSpecialization(toucher, touchAll);
    }

    updateRanges();
}


/**
 * Swap the parameters held by this distribution.
 *
 * \param[in]    oldP      Pointer to the old parameter.
 * \param[in]    newP      Pointer to the new parameter.
 */
void FossilizedBirthDeathRangeProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    AbstractFossilizedBirthDeathRangeProcess::swapParameterInternal(oldP, newP);
}

/**
 * Row i gives the probability that each range parented range i, and its diagonal entry the
 * probability that its parent was a species the record never sampled, which is what is left
 * when the origination does not sit at the divergence. Rows sum to one.
 *
 * Where the topology is marginalized, the density does not depend on which coexisting range
 * budded the daughter, so given the state every candidate is equally likely and the sampled
 * share divides evenly among them. Where the topology is sampled there is one candidate, so
 * the share lands on it. A range with no candidate at all is the origin, which has no parent,
 * and its whole row sits on the diagonal.
 *
 * A single row carries little, being flat across a set. The ancestry appears only in the
 * posterior average, as the divergences move and the candidate sets change with them.
 */
