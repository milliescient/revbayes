#ifndef TreeDistributionProperties_H
#define TreeDistributionProperties_H

#include "Taxon.h"

namespace RevBayesCore {

    /**
     * @brief Properties of a tree distribution that a tree move can ask about.
     *
     * The tip-age policies default to isExtended, for processes that answer per taxon.
     */
    class TreeDistributionProperties {

    public:
        virtual ~TreeDistributionProperties(void) {}

        virtual bool    isExtended(void) const { return false; }                                     //!< Tips are extinction times
        virtual bool    allowsSampledAncestors(void) const { return false; }                         //!< Checks if distribution is compatible with sampled ancestors
        virtual bool    tipAgeConstrainedToRange(const Taxon &t) const { return isExtended() == false; }
        virtual bool    validatesTipAgesOnSet(void) const { return isExtended() == false; }
    };

}

#endif
