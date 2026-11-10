#ifndef TreeDistributionProperties_H
#define TreeDistributionProperties_H

#include "Distribution.h"
#include "Taxon.h"

namespace RevBayesCore {

    /**
     * @brief Properties of a tree distribution that a tree move can ask about.
     *
     * The tip-age policy defaults to isExtended, for processes that answer per taxon.
     */
    class TreeDistributionProperties {

    public:
        virtual ~TreeDistributionProperties(void) {}

        virtual bool    isExtended(void) const { return false; }                                     //!< Tips are extinction times
        virtual bool    allowsSampledAncestors(void) const { return false; }                         //!< Checks if distribution is compatible with sampled ancestors
        virtual bool    tracksSpeciesContinuation(void) const { return false; }                      //!< Node flags say which child continues its parent's species
        virtual bool    tipAgeConstrainedToRange(const Taxon &t) const { return isExtended() == false; }
        virtual const Distribution* getWrappedDistribution(void) const { return NULL; }              //!< The distribution a wrapper such as a topology constraint holds
    };

    //!< The distribution behind any wrappers around d.
    inline Distribution& unconstrained( Distribution &d )
    {
        Distribution *out = &d;
        for (const TreeDistributionProperties *p = dynamic_cast<const TreeDistributionProperties *>( out ); p != NULL && p->getWrappedDistribution() != NULL; p = dynamic_cast<const TreeDistributionProperties *>( out ))
        {
            out = const_cast<Distribution *>( p->getWrappedDistribution() );
        }
        return *out;
    }

    inline const Distribution& unconstrained( const Distribution &d ) { return unconstrained( const_cast<Distribution&>( d ) ); }

}

#endif
