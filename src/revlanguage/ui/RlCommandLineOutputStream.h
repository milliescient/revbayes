#ifndef RlCommandLineOutputStream_h
#define RlCommandLineOutputStream_h

#include <iosfwd>

#include "RbConsole.h"

class CommandLineOutputStream : public RevBayesCore::ConsoleStream {
    
public:
    
    CommandLineOutputStream(void) {}
    virtual ~CommandLineOutputStream(void) {}
    
    void        output(const std::string &o) const;
    void        outputEndOfLine(void) const;
    
};


#endif /* RlCommandLineOutputStream_hpp */
