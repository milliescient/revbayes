#include "FossilRecordProcess.h"

#include "AbstractFossilizedBirthDeathRangeProcess.h"
#include "StochasticNode.h"
#include "MatrixReal.h"
#include "Tree.h"
#include "TimeInterval.h"
#include "RbException.h"

using namespace RevBayesCore;


/** Downcast a ranges stochastic node (MatrixReal-valued dnFBDRP or Tree-valued dnFBDSP) to the shared range base. */
AbstractFossilizedBirthDeathRangeProcess* FossilRecordProcess::rangesOf(const DagNode *n)
{
    const StochasticNode<MatrixReal> *mn = dynamic_cast<const StochasticNode<MatrixReal> *>( n );
    if ( mn != NULL )
    {
        return const_cast<AbstractFossilizedBirthDeathRangeProcess *>(
                   dynamic_cast<const AbstractFossilizedBirthDeathRangeProcess *>( &mn->getDistribution() ) );
    }
    const StochasticNode<Tree> *tn = dynamic_cast<const StochasticNode<Tree> *>( n );
    if ( tn != NULL )
    {
        return const_cast<AbstractFossilizedBirthDeathRangeProcess *>(
                   dynamic_cast<const AbstractFossilizedBirthDeathRangeProcess *>( &tn->getDistribution() ) );
    }
    return NULL;
}


/** The value is the range process's taxa, which carry the occurrences. */
FossilRecordProcess::FossilRecordProcess(const DagNode *rn,
                                         AbstractFossilizedBirthDeathRangeProcess::Reporting rule,
                                         const std::vector<double> &grid) :
    TypedDistribution< RbVector<Taxon> >( new RbVector<Taxon>() ),
    ranges_node( rn ),
    ranges( rangesOf( rn ) ),
    reporting( rule ),
    bins( grid ),
    taxa()
{
    if ( ranges == NULL )
    {
        throw RbException("dnFossilRecord requires a fossilized birth-death range process (e.g. dnFBDRP) as its first argument.");
    }

    addParameter( ranges_node );

    // psi and the timeline are read off the range process
    const std::vector<const DagNode*> &rp = ranges->getRangeParameters();
    for (std::vector<const DagNode*>::const_iterator it = rp.begin(); it != rp.end(); ++it)
    {
        addParameter( *it );
    }

    taxa = ranges->getTaxa();
    *this->value = RbVector<Taxon>( taxa );

    ranges->setReporting( reporting );
    ranges->setReportingBins( bins );
    ranges->setHasReportingNode();
}


FossilRecordProcess* FossilRecordProcess::clone( void ) const
{
    return new FossilRecordProcess( *this );
}


/** Downcast the ranges node's distribution to the range base. */
void FossilRecordProcess::resolveRanges( void )
{
    ranges = rangesOf( ranges_node );
}


/** Clamp a record, which must name the range process's taxa. */
void FossilRecordProcess::setValue(RbVector<Taxon> *v, bool force)
{
    const std::vector<Taxon> &ref = ranges->getTaxa();

    if ( v->size() != ref.size() )
    {
        throw RbException("dnFossilRecord: the clamped record has a different number of taxa than the range process.");
    }

    for (size_t i = 0; i < v->size(); ++i)
    {
        const Taxon &c = (*v)[i];

        const Taxon *r = NULL;
        for (size_t j = 0; j < ref.size(); ++j)
        {
            if ( ref[j].getName() == c.getName() ) { r = &ref[j]; break; }
        }
        if ( r == NULL )
        {
            throw RbException("dnFossilRecord: taxon '" + c.getName() + "' is not in the range process.");
        }
    }

    // the range process adopts the clamped occurrences
    std::vector<Taxon> adopted;
    for (size_t i = 0; i < ref.size(); ++i)
    {
        for (size_t j = 0; j < v->size(); ++j)
        {
            if ( (*v)[j].getName() == ref[i].getName() ) { adopted.push_back( (*v)[j] ); break; }
        }
    }
    ranges->setOccurrences( adopted );

    taxa = adopted;

    TypedDistribution< RbVector<Taxon> >::setValue( v, force );
}


double FossilRecordProcess::computeLnProbability( void )
{
    return ranges->computeLnFossilTotal();
}


void FossilRecordProcess::redrawValue( void )
{
    (*this->value) = taxa;
}


void FossilRecordProcess::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    if ( oldP == ranges_node )
    {
        ranges_node = newP;
        resolveRanges();
    }
}
