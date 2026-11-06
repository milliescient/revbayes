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


ExtinctionTimeGibbsProposal::ExtinctionTimeGibbsProposal( StochasticNode<Tree> *n, Target t ) : Proposal(),
    tree( n ),
    target( t ),
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
    failed = true;
    moved.clear();
    stored_ages.clear();
    flipped = NULL;

    dynamic_cast<FossilizedBirthDeathSpeciationProcess &>( tree->getDistribution() ).prepareForProposal();

    double hr = ( target == TIP ? proposeTip() : ( target == NODE ? proposeNode() : ( target == RANGE ? proposeRange() : proposeAncestor() ) ) );
    return hr;
}


/** Set tip's extinction to a draw from its conditional as the tree now stands, given its log density before. */
double ExtinctionTimeGibbsProposal::redraw( TopologyNode &tip, double ln_old )
{
    FossilizedBirthDeathSpeciationProcess &d = dynamic_cast<FossilizedBirthDeathSpeciationProcess &>( tree->getDistribution() );
    size_t i = tip.getIndex();

    double age = d.drawExtinctionTime( i );
    if ( std::isnan( age ) == true ) return RbConstants::Double::neginf;

    moved.push_back( &tip );
    stored_ages.push_back( tip.getAge() );
    tip.setAge( age );

    return ln_old - d.extinctionTimeLnDensity( i, age );
}


/** One fossil tip's extinction. The candidates are every extinct fossil tip, a number no move changes. */
double ExtinctionTimeGibbsProposal::proposeTip( void )
{
    RandomNumberGenerator* rng = GLOBAL_RNG;
    Tree &t = tree->getValue();
    FossilizedBirthDeathSpeciationProcess &d = dynamic_cast<FossilizedBirthDeathSpeciationProcess &>( tree->getDistribution() );

    std::vector<size_t> fossils;
    for (size_t i = 0; i < t.getNumberOfTips(); ++i)
    {
        const TopologyNode &n = t.getNode(i);
        if ( n.isFossil() == true && n.isRoot() == false ) fossils.push_back( i );
    }
    if ( fossils.empty() == true ) return RbConstants::Double::neginf;

    size_t i = fossils[ size_t( rng->uniform01() * fossils.size() ) ];
    double lo, hi;
    if ( d.extinctionWindow( i, lo, hi ) == false ) return RbConstants::Double::neginf;

    TopologyNode &tip = t.getNode(i);
    failed = false;
    return redraw( tip, d.extinctionTimeLnDensity( i, tip.getAge() ) );
}


/**
 * Slide an internal node uniformly between its parent and the oldest child whose age it keeps, redrawing each
 * extinct tip below it: those bound the slide only through their own windows, so the slide moves as though
 * their extinctions were integrated out.
 */
double ExtinctionTimeGibbsProposal::proposeNode( void )
{
    RandomNumberGenerator* rng = GLOBAL_RNG;
    Tree &t = tree->getValue();
    FossilizedBirthDeathSpeciationProcess &d = dynamic_cast<FossilizedBirthDeathSpeciationProcess &>( tree->getDistribution() );
    double present = d.getPresent();

    // a sampled ancestor at a symmetric speciation pins its tip, its parent and the symmetric node to one age,
    // which slides as one; no other proposal of this move changes which nodes those are
    std::vector<TopologyNode*> ok;
    for (size_t k = 0; k < t.getNumberOfNodes(); ++k)
    {
        TopologyNode &n = t.getNode(k);
        if ( n.isRoot() == true || n.isTip() == true ) continue;
        if ( n.isSampledAncestorTipOrParent() == false || collapsedSymmetric( n ) != NULL ) ok.push_back( &n );
    }
    if ( ok.empty() == true ) return RbConstants::Double::neginf;

    TopologyNode &node = *ok[ size_t( rng->uniform01() * ok.size() ) ];

    // the nodes that share its age, and the node whose children bound it from below
    std::vector<TopologyNode*> unit( 1, &node );
    TopologyNode *below = &node;
    double hi = node.getParent().getAge();
    TopologyNode *sym = collapsedSymmetric( node );
    if ( sym != NULL )
    {
        TopologyNode *sa = ( &node.getChild(0) == sym ? &node.getChild(1) : &node.getChild(0) );
        unit.push_back( sym );
        unit.push_back( sa );
        below = sym;
        hi = std::min( hi, d.getLastAppearance( sa->getIndex() ) );
    }

    double lo = present, ln_old = 0.0;
    std::vector<TopologyNode*> tips;
    for (size_t c = 0; c < below->getNumberOfChildren(); c++)
    {
        TopologyNode &child = below->getChild(c);
        double wlo, whi;
        if ( child.isTip() == true && d.extinctionWindow( child.getIndex(), wlo, whi ) == true )
        {
            tips.push_back( &child );
            ln_old += d.extinctionTimeLnDensity( child.getIndex(), child.getAge() );
        }
        else
        {
            lo = std::max( lo, child.getAge() );
        }
    }

    if ( hi <= lo ) return RbConstants::Double::neginf;

    double age = lo + rng->uniform01() * ( hi - lo );
    for (TopologyNode *n : unit)
    {
        moved.push_back( n );
        stored_ages.push_back( n->getAge() );
    }
    for (TopologyNode *n : unit) n->setAge( age );
    failed = false;

    double hr = ln_old;
    for (TopologyNode *tip : tips)
    {
        double r = redraw( *tip, 0.0 );
        if ( std::isfinite( r ) == false ) return RbConstants::Double::neginf;
        hr += r;
    }
    return hr;
}


/**
 * Make a fossil tip a sampled ancestor at its parent by raising its extinction to the parent's age, which must not
 * be older than its last appearance, or drop a sampled ancestor's extinction back below the node, drawn from its
 * conditional. The candidates are every extinct fossil tip, a number no move changes.
 */
double ExtinctionTimeGibbsProposal::proposeAncestor( void )
{
    RandomNumberGenerator* rng = GLOBAL_RNG;
    Tree &t = tree->getValue();
    FossilizedBirthDeathSpeciationProcess &d = dynamic_cast<FossilizedBirthDeathSpeciationProcess &>( tree->getDistribution() );

    std::vector<size_t> fossils;
    for (size_t i = 0; i < t.getNumberOfTips(); ++i)
    {
        const TopologyNode &n = t.getNode(i);
        if ( n.isFossil() == true && n.isRoot() == false ) fossils.push_back( i );
    }
    if ( fossils.empty() == true ) return RbConstants::Double::neginf;

    size_t i = fossils[ size_t( rng->uniform01() * fossils.size() ) ];
    TopologyNode &tip = t.getNode(i);
    TopologyNode &parent = tip.getParent();
    TopologyNode &sibling = ( &parent.getChild(0) == &tip ? parent.getChild(1) : parent.getChild(0) );
    // the parent's range ends at the sampled ancestor, so the tip continues it and the successor beside it does not
    if ( parent.getNumberOfChildren() != 2 || sibling.isSampledAncestorTip() == true ) return RbConstants::Double::neginf;
    if ( tip.continuesParentSpecies() == false || sibling.continuesParentSpecies() == true ) return RbConstants::Double::neginf;

    double lo, hi;
    if ( tip.isSampledAncestorTip() == true )
    {
        // expand: the extinction from its conditional below the node, as an ordinary tip there
        tip.setSampledAncestor( false );
        if ( d.extinctionWindow( i, lo, hi ) == false ) { tip.setSampledAncestor( true ); return RbConstants::Double::neginf; }
        moved.push_back( &tip );
        stored_ages.push_back( tip.getAge() );
        double age = d.drawExtinctionTime( i );
        if ( std::isnan( age ) == true ) { tip.setSampledAncestor( true ); return RbConstants::Double::neginf; }
        tip.setAge( age );
        flipped = &tip;
        failed = false;
        return -d.extinctionTimeLnDensity( i, age );
    }

    // collapse: the reverse expand draws the current extinction from the window it has now
    if ( d.extinctionWindow( i, lo, hi ) == false || parent.getAge() > d.getLastAppearance( i ) ) return RbConstants::Double::neginf;
    double ln_old = d.extinctionTimeLnDensity( i, tip.getAge() );
    moved.push_back( &tip );
    stored_ages.push_back( tip.getAge() );
    tip.setAge( parent.getAge() );
    tip.setSampledAncestor( true );
    flipped = &tip;
    failed = false;
    return ln_old;
}


/** For the parent of a sampled ancestor at a symmetric speciation, that symmetric node; otherwise NULL. */
TopologyNode* ExtinctionTimeGibbsProposal::collapsedSymmetric( TopologyNode &a ) const
{
    if ( a.getNumberOfChildren() != 2 ) return NULL;
    TopologyNode *sa = NULL, *b = NULL;
    for (size_t c = 0; c < 2; c++)
    {
        TopologyNode &k = a.getChild(c);
        if ( k.isSampledAncestorTip() == true ) sa = &k; else b = &k;
    }
    if ( sa == NULL || b == NULL || b->isTip() == true || b->getAge() != a.getAge() || b->continuesParentSpecies() == false ) return NULL;
    for (size_t c = 0; c < b->getNumberOfChildren(); c++)
    {
        if ( b->getChild(c).continuesParentSpecies() == true ) return NULL;
    }
    return b;
}


/** Resample a taxon's first and last appearance as mvStratigraphicRange does, and redraw its extinction below them. */
double ExtinctionTimeGibbsProposal::proposeRange( void )
{
    RandomNumberGenerator* rng = GLOBAL_RNG;
    Tree &t = tree->getValue();
    FossilizedBirthDeathSpeciationProcess &d = dynamic_cast<FossilizedBirthDeathSpeciationProcess &>( tree->getDistribution() );

    size_t i = size_t( rng->uniform01() * d.AbstractFossilizedBirthDeathRangeProcess::getTaxa().size() );

    double lo, hi;
    bool carries = d.extinctionWindow( i, lo, hi );
    double ln_old = ( carries ? d.extinctionTimeLnDensity( i, t.getNode(i).getAge() ) : 0.0 );

    // marks taxon i dirty and arms the range's own undo
    d.resampleFirstLast( i );
    failed = false;

    if ( carries == false ) return 0.0;
    if ( d.extinctionWindow( i, lo, hi ) == false ) return RbConstants::Double::neginf;
    return redraw( t.getNode(i), ln_old );
}


void ExtinctionTimeGibbsProposal::prepareProposal( void )
{
}


void ExtinctionTimeGibbsProposal::printParameterSummary(std::ostream &o, bool name_only) const
{
}


void ExtinctionTimeGibbsProposal::undoProposal( void )
{
    if ( flipped != NULL ) flipped->setSampledAncestor( !flipped->isSampledAncestorTip() );
    for (size_t k = moved.size(); k-- > 0; ) moved[k]->setAge( stored_ages[k] );
    moved.clear();
    flipped = NULL;
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
