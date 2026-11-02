#include "LnProbabilityFunction.h"
#include "Cloneable.h"
#include "StochasticNode.h"

namespace RevBayesCore { class DagNode; }

using namespace RevBayesCore;
/**
 * LnProbabilityFunction of a StochasticNode containing a value of type double
 *
 * @param x the value of type double
 *
 */
LnProbabilityFunction::LnProbabilityFunction(const DagNode *x) : ContinuousFunction( new double(0.0) ),
a( x )
{
    addParameter( a );
    // the probability changes with the parents too, not only with the value
    for (const DagNode *p : a->getParents()) addParameter( p );
}


LnProbabilityFunction* LnProbabilityFunction::clone( void ) const
{
    return new LnProbabilityFunction(*this);
}


void LnProbabilityFunction::swapParameterInternal(const DagNode *oldP, const DagNode *newP)
{
    
    if (oldP == a)
    {
        a = newP;
    }
    
}

void LnProbabilityFunction::update( void )
{
//    *value = a->
    *value = const_cast<DagNode* >(a)->getLnProbability();
}


