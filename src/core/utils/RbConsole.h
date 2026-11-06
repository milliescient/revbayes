#ifndef RbConsole_H
#define RbConsole_H

#include <sstream>
#include <string>

#define RBOUT(m) RevBayesCore::RbConsole::instance().output((m))

namespace RevBayesCore {

    /**
     * @brief Where RbConsole sends its text. A program installs its own with RbConsole::setOutputStream.
     */
    class ConsoleStream {

    public:
        virtual            ~ConsoleStream(void) {}

        virtual void        output(const std::string &o) const = 0;
        virtual void        outputEndOfLine(void) const = 0;
    };


    /**
     * @brief Core's messages to the user, written through RBOUT.
     *
     * Text goes to std::cout until a program installs its own stream.
     * Under MPI only the first process writes.
     */
    class RbConsole {

    public:
        static RbConsole&           instance(void);                                             //!< Get the one RbConsole

        void                        output(std::string msg);                                    //!< Write a message, padded and wrapped, and a newline
        void                        output(std::string msg, const bool hasPadding);             //!< Write a message and a newline, with or without padding
        void                        output(std::ostringstream msg);                             //!< Write a message from a stringstream
        void                        setOutputStream(ConsoleStream *o);                          //!< Send output to o, which the caller keeps

    private:
                                    RbConsole(void);
                                    RbConsole(const RbConsole& x) = delete;
        RbConsole&                  operator=(const RbConsole& x) = delete;

        int                         process_id;
        ConsoleStream*              output_stream;
    };

}

#endif
