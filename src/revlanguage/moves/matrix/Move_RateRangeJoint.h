#ifndef Move_RateRangeJoint_H
#define Move_RateRangeJoint_H

#include <ostream>

#include "RlMove.h"
#include "RevObject.h"
#include "RevPtr.h"
#include "RevVariable.h"

namespace RevLanguage {

    class TypeSpec;

    /** mvRateRangeJoint: a rate scale, shift flip or RJ shift switch that redraws the extinction times it affects. */
    class Move_RateRangeJoint : public Move {

    public:
        Move_RateRangeJoint(void);

        virtual Move_RateRangeJoint*                clone(void) const;
        void                                        constructInternalObject(void);
        static const std::string&                   getClassType(void);
        static const TypeSpec&                      getClassTypeSpec(void);
        std::string                                 getMoveName(void) const;
        const MemberRules&                          getParameterRules(void) const;
        virtual const TypeSpec&                     getTypeSpec(void) const;
        virtual void                                printValue(std::ostream& o) const;

    protected:
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);

        RevPtr<const RevVariable>                   x;
        RevPtr<const RevVariable>                   ranges;
        RevPtr<const RevVariable>                   lambda;
        RevPtr<const RevVariable>                   tune;
        RevPtr<const RevVariable>                   rj_switch;
        RevPtr<const RevVariable>                   sigma;
        RevPtr<const RevVariable>                   side;
    };

}

#endif
