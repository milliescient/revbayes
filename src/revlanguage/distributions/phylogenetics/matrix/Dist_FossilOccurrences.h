#ifndef Dist_FossilOccurrences_H
#define Dist_FossilOccurrences_H

#include "FossilOccurrenceDistribution.h"
#include "ModelVector.h"
#include "RlTaxon.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {

    /**
     * RevLanguage wrapper of the fossil record conditioned on a range process, dnFossilOccurrences.
     *
     *     rec ~ dnFossilOccurrences( ranges=bd, reporting="firstlast" )
     *     rec.clamp( taxa )
     */
    class Dist_FossilOccurrences : public TypedDistribution< ModelVector<Taxon> > {

    public:
        Dist_FossilOccurrences( void );

        Dist_FossilOccurrences*                         clone(void) const;
        static const std::string&                       getClassType(void);
        static const TypeSpec&                          getClassTypeSpec(void);
        std::string                                     getDistributionFunctionName(void) const;
        const TypeSpec&                                 getTypeSpec(void) const;
        const MemberRules&                              getParameterRules(void) const;

        RevBayesCore::FossilOccurrenceDistribution*     createDistribution(void) const;

    protected:
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);

        RevPtr<const RevVariable>                       ranges;
        RevPtr<const RevVariable>                       reporting;
        RevPtr<const RevVariable>                       bins;
    };

}

#endif
