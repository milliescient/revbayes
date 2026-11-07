#ifndef Move_CollapseExpandBifurcation_H
#define Move_CollapseExpandBifurcation_H

#include <ostream>
#include <string>
#include <vector>

#include "RlMove.h"
#include "RevPtr.h"
#include "RevVariable.h"

namespace RevLanguage {
    class TypeSpec;

    /**
     * @brief Rev wrapper class for the CollapseExpandBifurcation move.
     *
     * Maps a sampled ancestor that ends at a symmetric speciation to two buddings and a tip, and back.
     */
    class Move_CollapseExpandBifurcation : public Move {

    public:

        Move_CollapseExpandBifurcation(void);                                                                                          //!< Default constructor

        virtual Move_CollapseExpandBifurcation*        clone(void) const;                                                                      //!< Clone object
        void                                        constructInternalObject(void);                                                          //!< We construct the a new internal Move.
        static const std::string&                   getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                 getMoveName(void) const;                                                                //!< Get the name used for the constructor function in Rev.
        const MemberRules&                          getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object
        virtual void                                printValue(std::ostream& o) const;                                                      //!< Print value (for user)

    protected:

        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable

        RevPtr<const RevVariable>                   tree;
    };

}

#endif
