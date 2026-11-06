#include <iostream>
#include <string>

#include "RbConsole.h"
#include "RbSettings.h"
#include "RbUtil.h"
#include "StringUtilities.h"

#if defined (RB_MPI)
#include <mpi.h>
#endif

using namespace RevBayesCore;


namespace {

    class StandardStream : public ConsoleStream {

    public:
        void output(const std::string &o) const { std::cout << o; }
        void outputEndOfLine(void) const { std::cout << std::endl; }
    };

}


RbConsole::RbConsole( void ) :
    process_id( 0 ),
    output_stream( NULL )
{
    static StandardStream standard_output;
    output_stream = &standard_output;

#if defined (RB_MPI)
    MPI_Comm_rank(MPI_COMM_WORLD, &process_id);
#endif
}


RbConsole& RbConsole::instance( void )
{
    static RbConsole the_console;
    return the_console;
}


/** Print a message and a newline */
void RbConsole::output(std::string msg)
{

    if ( process_id == 0 )
    {
        std::string tmp = StringUtilities::formatStringForScreen( msg, RbUtils::PAD, RbUtils::PAD, RbSettings::userSettings().getLineWidth() );
        output_stream->output(tmp);
    }
    
}


/** Print a message and a newline without the padding */
void RbConsole::output(std::string msg, const bool hasPadding)
{

    if ( process_id == 0 )
    {
        if (hasPadding == true)
        {
            output(msg);
        }
        else
        {
            output_stream->output(msg);
            output_stream->outputEndOfLine();
        }
        
    }
}


/** Convert to string and then call output to print message string */
void RbConsole::output(std::ostringstream msg)
{

    output( msg.str() );
}


void RbConsole::setOutputStream(ConsoleStream *o)
{
    output_stream = o;
}
