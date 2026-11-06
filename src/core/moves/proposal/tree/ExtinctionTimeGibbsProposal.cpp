#include "ExtinctionTimeGibbsProposal.h"

#include <cmath>
#include <vector>

#include "FossilizedBirthDeathSpeciationProcess.h"
#include "RandomNumberFactory.h"
#include "RandomNumberGenerator.h"
#include "RbConstants.h"
#include "RbException.h"
#include "TopologyNode.h"
#include "TreeDistributionProperties.h"

using namespace RevBayesCore;


ExtinctionTimeGibbsProposal::ExtinctionTimeGibbsProposal( StochasticNode<Tree> *n ) : Proposal(),
    tree( n ),
    tip_index( 0 ),
    stored_age( 0.0 ),
    failed( false )
{
    addNode( tree );

    const FossilizedBirthDeathSpeciationProcess *d = dynamic_cast<const FossilizedBirthDeathSpeciationProcess *>( &tree->getDistribution() );
    if ( d == NULL || static_cast<const TreeDistributionProperties *>( d )->isExtended() == false )
    {
        throw RbException("mvExtinctionTimeGibbs needs an extended dnFBDSP tree.");
    }
}


void ExtinctionTimeGibbsProposal::cleanProposal( void )
{
}


ExtinctionTimeGibbsProposal* ExtinctionTimeGibbsProposal::clone( void ) const
{
    return new ExtinctionTimeGibbsProposal( *this );
}


const std::string& ExtinctionTimeGibbsProposal::getProposalName( void ) const
{
    static std::string name = "ExtinctionTimeGibbs";
    return name;
}


double ExtinctionTimeGibbsProposal::getProposalTuningParameter( void ) const
{
    return RbConstants::Double::nan;
}


double ExtinctionTimeGibbsProposal::doProposal( void )
{
    RandomNumberGenerator* rng = GLOBAL_RNG;
    Tree &t = tree->getValue();
    FossilizedBirthDeathSpeciationProcess &d = dynamic_cast<FossilizedBirthDeathSpeciationProcess &>( tree->getDistribution() );

    failed = true;

    std::vector<size_t> fossils;
    for (size_t i = 0; i < t.getNumberOfTips(); ++i)
    {
        const TopologyNode &n = t.getNode(i);
        if ( n.isFossil() == true && n.isRoot() == false ) fossils.push_back( i );
    }
    if ( fossils.empty() == true ) return RbConstants::Double::neginf;

    tip_index = fossils[ size_t( rng->uniform01() * fossils.size() ) ];

    d.prepareForProposal();
    double lo, hi;

    TopologyNode &tip = t.getNode( tip_index );
    stored_age = tip.getAge();

    double ln_old = d.extinctionTimeLnDensity( tip_index, stored_age );
    double age = d.drawExtinctionTime( tip_index );
    double ln_new = d.extinctionTimeLnDensity( tip_index, age );

    tip.setAge( age );
    failed = false;

    return ln_old - ln_new;
}


void ExtinctionTimeGibbsProposal::prepareProposal( void )
{
}


void ExtinctionTimeGibbsProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
}


void ExtinctionTimeGibbsProposal::undoProposal( void )
{
    if ( failed == false ) tree->getValue().getNode( tip_index ).setAge( stored_age );
}


void ExtinctionTimeGibbsProposal::setProposalTuningParameter(double tp)
{
}


void ExtinctionTimeGibbsProposal::swapNodeInternal(DagNode *oldN, DagNode *newN)
{
    if ( oldN == tree ) tree = static_cast<StochasticNode<Tree>* >( newN );
}


void ExtinctionTimeGibbsProposal::tune( double rate )
{
}
