//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/signalhandler.cpp $
//$LastChangedDate: 2014-05-25 19:25:40 +0200 (So, 25. Mai 2014) $
//$LastChangedRevision: 1710 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include "signalhandler.h"

#include <exception>
#include <iostream>

#include <cxxabi.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>     /* malloc, calloc, realloc, free */

void signal_handler( const int signo )
{
  void *stack[20];
  int count;
  
  std::cerr <<__PRETTY_FUNCTION__<< ": Caught signal " << signo << std::endl;
  
  count = backtrace(stack, 20);

  // The short way:
  //  backtrace_symbols_fd(stack, count, 2);

  // The long way:
  
  char ** messages = backtrace_symbols(stack, count);

  // skip first stack frame (points here)
  for (int i = 1; i < count && messages != NULL; ++i)
  {
    char *mangled_name = 0, 
      *offset_begin = 0, *offset_end = 0,
      *adress_begin = 0, *adress_end = 0;

    //    std::cerr << messages[i] << std::endl;
    
    // find parantheses and +address offset surrounding mangled name
    for (char *p = messages[i]; *p; ++p)
    {
      switch (*p)
      {
      case '(': mangled_name = p; break;
      case '+': offset_begin = p; break;
      case ')': offset_end = p; break;
      case '[': adress_begin = p; break;
      case ']': adress_end = p; break;
      }
    }

    // if the line could be processed, attempt to demangle the symbol
    if (mangled_name && offset_begin && offset_end && 
        mangled_name < offset_begin)
    {
      *mangled_name++ = '\0';
      *offset_begin++ = '\0';
      *offset_end++ = '\0';

      if (adress_begin && adress_end)
      {
        *adress_begin++ = '\0';
        *adress_end++ = '\0';
      }
      
      int status;
      char * real_name = abi::__cxa_demangle(mangled_name, 0, 0, &status);
      
      // if demangling is successful, output the demangled function name
      if (status == 0)
      {    
        std::cerr << "(" << i << ") " << messages[i] << " : " 
                  << real_name << "+" << offset_begin << offset_end;
        
        if (adress_begin && adress_end)
        {
          std::cerr << "[" << adress_begin << "]"
                    << " # addr2line -e " << messages[i]
                    << " -a " << adress_begin << " -pfC";
        }
        std::cerr << std::endl;

      }
      // otherwise, output the mangled function name
      else
      {
        std::cerr << "(" << i << ") " << messages[i] << " : " 
                  << mangled_name << "+" << offset_begin << offset_end;

        if (adress_begin && adress_end)
        {
          std::cerr << "[" << adress_begin << "]"
                    << " # addr2line -e " << messages[i] 
                    << " -a " << adress_begin << " -pfC";
        }
        std::cerr << std::endl;
      }
      free(real_name);
    }
    // otherwise, print the whole line
    else
    {
      std::cerr << "(" << i << ") " << messages[i] << std::endl;
    }
  }
  std::cerr << "..." << std::endl;
  std::cerr << std::endl;
  
  free(messages);

  //  std::string errMsg = "Error thrown by signal_handler";
  //  throw eSignal_error( errMsg );

  exit( signo );
}


void registerSignalHandler(void)
{
  signal(SIGBUS, signal_handler);
  signal(SIGILL, signal_handler);
  signal(SIGSEGV, signal_handler);
  signal(SIGABRT, signal_handler);
}


