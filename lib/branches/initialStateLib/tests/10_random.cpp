//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/tests/10_random.cpp $
//$LastChangedDate: 2016-09-21 12:14:54 +0200 (Mi, 21. Sep 2016) $
//$LastChangedRevision: 2433 $
//$LastChangedBy: senzel $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/*
  This is some simple executable for testing.
  Here we test the possibilities to generate random numbers.
*/

#include <iostream>

#if _OPENMP
#include <omp.h>
#endif

#include "random.h"
#include "TimerClock.h"

using namespace std;

#if _OPENMP
#define ITHREAD omp_get_thread_num()
#else
#define ITHREAD 0
#endif

void test0(void)
{
  std::cout << "ran2: " << ran2.Name() << endl;  

  unsigned int nThreads=1;
  unsigned int nThreadsMax=0;

#if _OPENMP
#pragma omp parallel
  {
    nThreadsMax = omp_get_max_threads();
    nThreads    = omp_get_num_threads();
  }
#endif

  std::cout << "   nThreadsMax = " << nThreadsMax
            << "   nThreads    = " << nThreads
            << endl;

  
  ran2.setSeed( 13 );

#pragma omp parallel for ordered
  for (unsigned int i=0; i < nThreads*2; i++)
  {
#pragma omp ordered
    std::cout << i << "  " << ITHREAD << "  " << ran2() << std::endl;
  }

ran2.setSeed( 113 );

#pragma omp parallel for ordered
  for (unsigned int i=0; i < nThreads*2; i++)
  {
#pragma omp ordered
    std::cout << i << "  " << ITHREAD << "  " << ran2() << std::endl;
  }

  ran2.setSeed( 13 );

#pragma omp parallel for ordered
  for (unsigned int i=0; i < nThreads*2; i++)
  {
#pragma omp ordered
    std::cout << i << "  " << ITHREAD << "  " << ran2() << std::endl;
  }



}

void test0_int(void)
{
  std::cout << "ran2: " << ran2.Name() << endl;  

  unsigned int nThreads=1;
  unsigned int nThreadsMax=0;

#if _OPENMP
#pragma omp parallel
  {
    nThreadsMax = omp_get_max_threads();
    nThreads    = omp_get_num_threads();
  }
#endif

  std::cout << "   nThreadsMax = " << nThreadsMax
            << "   nThreads    = " << nThreads
            << endl;

  
  std::cout << std::endl << "setSeed: 13" << std::endl;
  ran2.setSeed( 13 );

#pragma omp parallel for ordered
  for (unsigned int i=0; i < nThreads*2; i++)
  {
#pragma omp ordered
    std::cout << i << "  " << ITHREAD << "  " << ran2.Int(99) << std::endl;
  }

  std::cout << std::endl << "setSeed: 113" << std::endl;
  ran2.setSeed( 113 );

#pragma omp parallel for ordered
  for (unsigned int i=0; i < nThreads*2; i++)
  {
#pragma omp ordered
    std::cout << i << "  " << ITHREAD << "  " << ran2.Int(99) << std::endl;
  }

  std::cout << std::endl << "setSeed: 13" << std::endl;
  ran2.setSeed( 13 );

#pragma omp parallel for ordered
  for (unsigned int i=0; i < nThreads*2; i++)
  {
#pragma omp ordered
    std::cout << i << "  " << ITHREAD << "  " << ran2.Int(99) << std::endl;
  }



}


void test1(void)
{
  std::cout << "ran2: " << ran2.Name() << endl;  

  TimerClock t1;
#ifdef  __MACH__ // OS X does not have clock_gettime, use clock_get_time
  TimerClock t0;
#else  
  TimerClock t0(CLOCK_PROCESS_CPUTIME_ID);
#endif

  double average = 0;
  unsigned int N = 100000000;

  int nThreads=0;
  int nThreadsMax=0;

#if _OPENMP
#pragma omp parallel
  {
    nThreadsMax = omp_get_max_threads();
    nThreads    = omp_get_num_threads();
  }
#endif

#pragma omp parallel for reduction(+:average)
  for (unsigned int i=0;i<N;i++)
  {
    average += ran2();
  }

  t1.Stop();
  t0.Stop();


  std::cout << "ave: " << average/N 
            << "   t1: " << t1.Nanosecs()
            << "   t0: " << t0.Nanosecs()
            << "   nThreadsMax = " << nThreadsMax
            << "   nThreads    = " << nThreads
            << "   speedUp = " << (1.0*t0.Nanosecs())/t1.Nanosecs()
            << endl;

}

int main(int argc, char **argv)
{

  try
  {
    //    test0();
    test0_int();
    //    test1();
  }
  catch (int e)
  {
    cout << "An exception occurred. Exception Nr. " << e << endl;
    return e;
  }
  return 0;  
}
