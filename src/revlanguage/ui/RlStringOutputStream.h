#ifndef RlStringOutputStream_h
#define RlStringOutputStream_h

#include <iosfwd>

#include "RbConsole.h"

class StringOutputStream : public RevBayesCore::ConsoleStream {
    
public:
    
    StringOutputStream(void) {}
    virtual ~StringOutputStream(void) {}
    
    void                    emptyOutputString(void);
    const std::string&      getOutputString(void) const;
    void                    output(const std::string &o) const;
    void                    outputEndOfLine(void) const;
    
private:
    
    mutable std::string     current_output;
    
};


#endif /* RlStringOutputStream_hpp */
