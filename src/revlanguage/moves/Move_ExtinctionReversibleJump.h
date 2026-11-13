#ifndef Move_ExtinctionReversibleJump_H
#define Move_ExtinctionReversibleJump_H

#include "RlMove.h"
#include "TypedDagNode.h"

#include <ostream>
#include <string>

namespace RevLanguage {
    
    
    /**
     * @brief Rev wrapper for the extinction reversible-jump move.
     *
     * Moves one taxon's extinction time between the present and a time above it.
     *
     * @author The RevBayes Development Core Team (June Walker)
     * @copyright GPL version 3
     */
    class Move_ExtinctionReversibleJump : public Move {
        
    public:
        
        Move_ExtinctionReversibleJump(void);                                                                                                                    //!< Default constructor
        
        // Basic utility functions
        virtual Move_ExtinctionReversibleJump*                   clone(void) const;                                                                      //!< Clone object
        void                                        constructInternalObject(void);                                                          //!< We construct the a new internal Move.
        static const std::string&                   getClassType(void);                                                                     //!< Get Rev type
        static const TypeSpec&                      getClassTypeSpec(void);                                                                 //!< Get class type spec
        std::string                                 getMoveName(void) const;                                                                //!< Get the name used for the constructor function in Rev.
        const MemberRules&                          getParameterRules(void) const;                                                          //!< Get member rules (const)
        virtual const TypeSpec&                     getTypeSpec(void) const;                                                                //!< Get language type of the object
        virtual void                                printValue(std::ostream& o) const;                                                      //!< Print value (for user)
        
    protected:
        
        void                                        setConstParameter(const std::string& name, const RevPtr<const RevVariable> &var);       //!< Set member variable
        
        RevPtr<const RevVariable>                   x;
        
    };
    
}

#endif
