//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/TSC.h $
//$LastChangedDate: 2015-02-07 22:53:32 +0100 (Sa, 07. Feb 2015) $
//$LastChangedRevision: 2078 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file 
 * @brief Usage of Time Stamp Counter for timing purposes
 *
 * please read the relevant literature when using this class, at least
 *  * http://download.intel.com/embedded/software/IA/324264.pdf ("How to Benchmark Code Execution Times on Intel IA-32 and IA-64 Instruction Set Architectures [White Paper]")
 *  * http://www.ccsl.carleton.ca/~jamuir/rdtscpm1.pdf ("Using the RDTSC Instruction for Performance Monitoring")
 *  * http://en.wikipedia.org/wiki/Time_Stamp_Counter
 *       and references therein
 */

#ifndef TSC_H
#define TSC_H

#include "configBAMPS.h"
//#define HAVE_RTDSCP 1


// #ifdef _MSC_VER
// #warning "TSC: MSC"
// #else
// #if HAVE_RDTSCP
// #warning "TSC: HAVE_RDTSCP = yes"
// #else
// #warning "TSC: HAVE_RDTSCP = no"
// #endif
// #endif

/** 
 * @brief class to get access to the time stamp counter
 */
class TimeStampCounter
{
public:
  void Start();
  void Stop();
  unsigned long long Cycles() const;
  
private:
  union Data {
    unsigned long long a;
    unsigned int b[2];
  } m_start, m_end;
};

inline void TimeStampCounter::Start()
{
#ifdef __arm__
  #warning "Running on ARM. No TSC available in userland"
  m_start.a = 0ULL;
#else

#ifdef _MSC_VER
  unsigned int tmp;
  m_start.a = __rdtscp(&tmp);
#else
#if HAVE_RDTSCP
#ifdef __i386
  asm volatile("xorl %%eax,%%eax \n"
               "cpuid"      // serialize
               ::: "%eax", "%ebx", "%ecx", "%edx");
#else
  asm volatile("xorl %%eax,%%eax \n"
               "cpuid"      // serialize
               ::: "%rax", "%rbx", "%rcx", "%rdx");
#endif
  asm volatile("rdtscp" : "=a"(m_start.b[0]), "=d"(m_start.b[1]) :: "ecx" );
#else
#ifdef __i386
  asm volatile("xorl %%eax,%%eax \n"
               "cpuid"      // serialize
               ::: "%eax", "%ebx", "%ecx", "%edx");
#else
  asm volatile("xorl %%eax,%%eax \n"
               "cpuid"      // serialize
               ::: "%rax", "%rbx", "%rcx", "%rdx");
#endif
  asm volatile("rdtsc" : "=a"(m_start.b[0]), "=d"(m_start.b[1]));
#endif
#endif
#endif
}

inline void TimeStampCounter::Stop()
{
#ifdef __arm__
  m_end.a = 0ULL;
#else

#ifdef _MSC_VER
  unsigned int tmp;
  m_end.a = __rdtscp(&tmp);
#else
#if HAVE_RDTSCP
  asm volatile("rdtscp" : "=a"(m_end.b[0]), "=d"(m_end.b[1]) :: "ecx" );
#ifdef __i386
  asm volatile("xorl %%eax,%%eax \n"
               "cpuid"      // serialize
               ::: "%eax", "%ebx", "%ecx", "%edx");
#else
  asm volatile("xorl %%eax,%%eax \n"
               "cpuid"      // serialize
               ::: "%rax", "%rbx", "%rcx", "%rdx");
#endif
#else
#ifdef __i386
  asm volatile("xorl %%eax,%%eax \n"
               "cpuid"      // serialize
               ::: "%eax", "%ebx", "%ecx", "%edx");
#else
  asm volatile("xorl %%eax,%%eax \n"
               "cpuid"      // serialize
               ::: "%rax", "%rbx", "%rcx", "%rdx");
#endif  
  asm volatile("rdtsc" : "=a"(m_end.b[0]), "=d"(m_end.b[1]));
#endif
#endif

#endif
}

inline unsigned long long TimeStampCounter::Cycles() const
{
  return m_end.a - m_start.a;
}

TimeStampCounter timer;


#endif
