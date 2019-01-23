//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/logger.h $
//$LastChangedDate: 2015-11-26 18:04:02 +0100 (Do, 26. Nov 2015) $
//$LastChangedRevision: 2234 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file 
 * @brief Declaration for routines implementing a logger
 **/

#ifndef LOGGER_H
#define LOGGER_H

#include "configBAMPS.h"

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifdef LOG4CXX_FOUND

//---------- with Log4cxx: ---------------

#include "log4cxx/logger.h"
#include "log4cxx/basicconfigurator.h"
#include "log4cxx/propertyconfigurator.h"
#include "log4cxx/helpers/exception.h"

#define LOGGER_t log4cxx::LoggerPtr
#define GETLOGGER log4cxx::Logger::getLogger

#define LOG_TRACE  LOG4CXX_TRACE
#define LOG_DEBUG  LOG4CXX_DEBUG
#define LOG_INFO   LOG4CXX_INFO 
#define LOG_WARN   LOG4CXX_WARN 
#define LOG_ERROR  LOG4CXX_ERROR
#define LOG_FATAL  LOG4CXX_FATAL

#define LOG_CONFIGURE_FILE(file) log4cxx::PropertyConfigurator::configure(file)
#define LOG_CONFIGURE() log4cxx::BasicConfigurator::configure()

//----------------------------------------

#else

//---------- no Log4cxx: -----------------

//cf. http://www.drdobbs.com/parallel/logging-in-c/201804215?pgno=1

#include <sstream>
#include <string>
#include <stdio.h>
#include <sys/time.h>

enum class tLogLevel { FATAL, ERROR, WARN, INFO, DEBUG, TRACE };

inline std::string NowTime()
{
  char buffer[11];
  time_t t;
  time(&t);
  tm r = {0};
  strftime(buffer, sizeof(buffer), "[%X]", localtime_r(&t, &r));
  return buffer;
  // struct timeval tv;
  // gettimeofday(&tv, 0);
  // char result[100] = {0};
  // sprintf(result, "%s.%03ld", buffer, (long)tv.tv_usec / 1000); 
  // return result;
}

class tLogMaster
{
public:
  tLogMaster(const std::string &_name): name(_name),level(tLogLevel::TRACE) {};
  tLogMaster(const char * _name): name(_name),level(tLogLevel::TRACE) {};
  std::string name;
  tLogLevel level;
};

class tLog
{
public:
  tLog() {};
  tLog(const tLogMaster & M) {};
  virtual ~tLog();
  std::ostringstream& Get(const char* _file, const int _line, const std::string &_name,tLogLevel level = tLogLevel::INFO);
  static std::string ToString(tLogLevel level);
protected:
  std::ostringstream os;
private:
  tLog(const tLog&) {}; ///< copy constructor, private!
  tLog& operator =(const tLog&); ///< asignment operator, private!
};



inline std::ostringstream& tLog::Get(const char* _file, const int _line, const std::string &_name, tLogLevel level)
{
  os << ToString(level)
     << " " << NowTime()
     << " " << _name
     << " (" << _line << ")"
     << " - ";
  return os;
}

inline tLog::~tLog()
{
  os << std::endl;
  fprintf(stderr, "%s", os.str().c_str());
  fflush(stderr);
}

inline std::string tLog::ToString(tLogLevel level)
{
  switch (level)
  {
    case tLogLevel::FATAL: return "FATAL"; break;
    case tLogLevel::ERROR: return "ERROR"; break;
    case tLogLevel::WARN:  return "WARN "; break;
    case tLogLevel::INFO:  return "INFO "; break;
    case tLogLevel::DEBUG: return "DEBUG"; break;
    case tLogLevel::TRACE: return "TRACE"; break;
    default: return "#####";
  }
}

#ifndef LOG_MAX_LEVEL
//#define LOG_MAX_LEVEL tLogLevel::DEBUG
#define LOG_MAX_LEVEL tLogLevel::TRACE
#endif

#define _LOG(_logger,_level)                                    \
  if (_level > LOG_MAX_LEVEL) ;                                 \
  else if (_level > logger.level) ;                             \
  else tLog().Get(__FILE__,__LINE__,_logger.name,_level)

#define LOGGER_t tLogMaster
#define GETLOGGER(name) name

#define LOG_TRACE(logger, message) { _LOG(logger,tLogLevel::TRACE) << message; }
#define LOG_DEBUG(logger, message) { _LOG(logger,tLogLevel::DEBUG) << message; }
#define LOG_INFO(logger, message)  { _LOG(logger,tLogLevel::INFO) << message; }
#define LOG_WARN(logger, message)  { _LOG(logger,tLogLevel::WARN) << message; }
#define LOG_ERROR(logger, message) { _LOG(logger,tLogLevel::ERROR) << message; }
#define LOG_FATAL(logger, message) { _LOG(logger,tLogLevel::FATAL) << message; }

#define LOG_CONFIGURE_FILE(file) { }
#define LOG_CONFIGURE() { }

//----------------------------------------
#endif

#endif
