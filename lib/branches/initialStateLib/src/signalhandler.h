//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/signalhandler.h $
//$LastChangedDate: 2014-04-30 20:28:28 +0200 (Mi, 30. Apr 2014) $
//$LastChangedRevision: 1661 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file 
 * @brief Declarations of some default signal handler
 *
 * Here we define some signal handler, which may be used to catch
 * errors like #SIGBUS, #SIGSEV etc. Its main feature is a traceback
 * with demangled C++ names.
 */


#ifndef SIGNALHANDLER_H
#define SIGNALHANDLER_H

#include <stdexcept>

/** 
 * @brief exception class for handling unexpected critical behaviour
 * due to #signal_handler
 */
class eSignal_error : public std::runtime_error
{
  public:
    explicit eSignal_error(const std::string& what) : std::runtime_error(what) {};
    
    virtual ~eSignal_error() throw() {};
};

/**
 * @brief Handler for system signals received by the program
 * 
 * This handles system signals such as KILL, SEGV, etc. It needs to be
 * registered in the main routine via 
 * `signal(<SIGNAME>,signal_handler)`
 * as e.g.
 * ~~~
 * signal(SIGBUS, signal_handler);
 * signal(SIGILL, signal_handler);
 * signal(SIGSEGV, signal_handler);
 * signal(SIGABRT, signal_handler);
 * ~~~
 * You may also use the routine #registerSignalHandler() for this.

 * cf. http://stackoverflow.com/questions/77005/how-to-generate-a-stacktrace-when-my-gcc-c-app-crashes
 * 
 * @param[in] signo ID of the received signal
 */ 
void signal_handler( const int signo );

/**
 * @brief routine to register `signal_handler` as the handler for various signals
 **/
void registerSignalHandler(void);

#endif
