//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/TimerClock.h $
//$LastChangedDate: 2016-09-21 12:14:54 +0200 (Mi, 21. Sep 2016) $
//$LastChangedRevision: 2433 $
//$LastChangedBy: senzel $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file
 * @brief This file implements some timing routines.
 **/

#ifndef TIMERCLOCK_H
#define TIMERCLOCK_H

#include <time.h>
#include <sys/time.h>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

// struct timespec {
// 	time_t tv_sec; /* seconds */
// 	long tv_nsec; /* nanoseconds */
// };


/**
 * @brief class to implement a timer
 *
 * This class uses clock_gettime() for 'measuring' the time. You may
 * select, which kind of timer you would like to use:
 *   - CLOCK_REALTIME: the actual wall clock time
 *   - CLOCK_MONOTONIC: for measuring relative real time. It advances
 *     at the same rate as the actual flow of time but it's not
 *     subject to discontinuities from manual or automatic (NTP)
 *     adjustments to the system clock. 
 *   - CLOCK_PROCESS_CPUTIME_ID: amount of CPU time consumed by the
 *     process 
 *   - CLOCK_THREAD_CPUTIME_ID: amount of CPU time consumed by the
 *     thread 
 *
 * see also the man pages: "man clock_gettime"
 **/
class TimerClock
{
public:

#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
  TimerClock() :
    time1(),
    time2()
  {
    Start();
  };
#else  
  TimerClock() :
    time1(),
    time2(),
    clk_id( CLOCK_MONOTONIC )
  {
    Start();
  };

  TimerClock( clockid_t _clk_id ) :
    time1(),
    time2(),
    clk_id( _clk_id )
  { 
    Start();
  };
#endif

  void Start()
  {
    #ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    time1.tv_sec = mts.tv_sec;
    time1.tv_nsec = mts.tv_nsec;
    #else
    clock_gettime(clk_id, &time1);
    #endif
  }

  void Stop()
  {
    #ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
    clock_serv_t cclock;
    mach_timespec_t mts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &mts);
    mach_port_deallocate(mach_task_self(), cclock);
    time2.tv_sec = mts.tv_sec;
    time2.tv_nsec = mts.tv_nsec;
    #else
    clock_gettime(clk_id, &time2);
    #endif
  }

  long long Nanosecs() const
  {
    return (time2.tv_nsec-time1.tv_nsec)
      + (long long)1000000000*(long long)(time2.tv_sec-time1.tv_sec);
  }

#ifndef __MACH__ // OS X does not have clock_gettime, use clock_get_time
  timespec RES_PROCESS_CPUTIME_ID(void) 
  { 
    timespec time1;
    clock_getres(CLOCK_PROCESS_CPUTIME_ID, &time1);
    return time1;
  }
  timespec RES_MONOTONIC(void) 
  { 
    timespec time1;
    clock_getres(CLOCK_MONOTONIC, &time1);
    return time1;
  }
#endif

private:
  timespec time1, time2;

#ifndef __MACH__ // OS X does not have clock_gettime, use clock_get_time 
  clockid_t clk_id;
#endif
};



/* TimerClock timerClock;*/

#endif
