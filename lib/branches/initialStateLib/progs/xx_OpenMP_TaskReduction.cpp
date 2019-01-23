//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/progs/xx_OpenMP_TaskReduction.cpp $
//$LastChangedDate: 2015-01-09 14:43:42 +0100 (Fr, 09. Jan 2015) $
//$LastChangedRevision: 2036 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/*
  Here we try to check some OpenMP issues
*/

#include <iostream>     // std::cout, std::fixed, std::scientific
#include <cstdlib>

#if _OPENMP
#include <omp.h>
#endif

#if _OPENMP
#define ITHREAD omp_get_thread_num()
#else
#define ITHREAD 0
#endif


void test1(void)
{
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
            << std::endl;

  int sum1 = 0, sum2 = 0, sum3 = 0;

  std::cout << "&sum1 = " << &sum1 << std::endl;

#pragma omp parallel reduction( + : sum1,sum2,sum3 )
  {
#pragma omp single
    {
      std::cout << "&sum1 = " << &sum1 << std::endl;

      for(int i=0; i < 100; ++i)
        //#pragma omp task firstprivate(i)
#pragma omp task firstprivate(i) shared(sum1,sum2,sum3)
        //#pragma omp task firstprivate(i) threadprivate(sum1,sum2,sum3)
      {

#pragma omp critical
        std::cout << "task " << i << "  thread " << ITHREAD 
                  << "  &sum1 = " << &sum1
                  << std::endl;

        sum1 += 1;
        sum2 += i;
        sum3 += ITHREAD;
      }

#pragma omp taskwait
    }
  }

  std::cout << "  sum1 = " << sum1
            << "  sum2 = " << sum2
            << "  sum3 = " << sum3
            << std::endl;

}

/**
 * @brief The main function of the testing programs
 **/
int main(int argc, char **argv)
{

  try
  {
    test1();
  }
  catch (int e)
  {
    std::cout << "An exception occurred. Exception Nr. " << e << std::endl;
    return e;
  }
  return 0;  
}
