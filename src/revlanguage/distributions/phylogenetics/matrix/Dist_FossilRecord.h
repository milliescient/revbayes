#ifndef Dist_FossilRecord_H
#define Dist_FossilRecord_H

#include "FossilRecordProcess.h"
#include "ModelVector.h"
#include "RlTaxon.h"
#include "RlTypedDistribution.h"

namespace RevLanguage {

    /**
     * RevLanguage wrapper of the fossil record conditioned on a range process, dnFossilRecord.
     *
     *     rec ~ dnFossilRecord( ranges=bd, reporting="firstlast" )
     *     rec.clamp( taxa )
     */
    class Dist_FossilRecord : public TypedDistribution< ModelVector<Taxon> > {

    public:
        Dist_FossilRecord( void );

        Dist_FossilRecord*                              clone(void) const;
        static const std::string&                       getClassType(void);
        static const TypeSpec&                          getClassTypeSpec(void);
        std::string                                     getDistributionFunctionName(void) const;
        const TypeSpec&                                 getTypeSpec(void) const;
        const MemberRules&                              getParameterRules(void) const;

        RevBayesCore::FossilRecordProcess*              createDistribution(void) const;

    protected:
        void                                            setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);

        RevPtr<const RevVariable>                       ranges;
        RevPtr<const RevVariable>                       bins;   //!< The reporting grid, which is not the rate timeline.
        RevPtr<const RevVariable>                       reporting;
    };

}

#endif
