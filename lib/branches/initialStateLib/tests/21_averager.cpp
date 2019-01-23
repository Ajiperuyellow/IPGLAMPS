//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/tests/21_averager.cpp $
//$LastChangedDate: 2014-12-11 22:00:00 +0100 (Do, 11. Dez 2014) $
//$LastChangedRevision: 2015 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/**
 * @brief This is a zero content example for the testing stuff.
 *
 * We do not test anything here, but just provide the basic skeleton
 * of such a file.
 *
 * Running this file as an executable with a return value ==0 is
 * interpreted as *success*. You may indicate the internal number of
 * the failing test by any value != 0.
 **/

#include <iostream>

#include "FPT_compare.h"
#include "averager.h"

/**
 * @brief The function doing the actual work
 *
 * Here we throw an integer in order to indicate, which internal test
 * failed.
 **/

void test1(void)
{
  tAverager Ave;

  if (!FPT_COMP_E(Ave.Mean(),0.0))
  {
    std::cout << "Ave.Mean() = " << Ave.Mean() << std::endl;
    throw 1;
  }

  if (!FPT_COMP_E(Ave.Var(),0.0))
  {
    std::cout << "Ave.Var() = " << Ave.Var() << std::endl;
    throw 2;
  }

  if (!FPT_COMP_E(Ave.VarA(),0.0))
  {
    std::cout << "Ave.VarA() = " << Ave.VarA() << std::endl;
    throw 3;
  }

  

  Ave.add(2.0);
  Ave.add(4.0);
  Ave.add(6.0);

  if (!FPT_COMP_E(Ave.Mean(),4.0))
  {
    std::cout << "Ave.Mean() = " << Ave.Mean() << std::endl;
    throw 11;
  }

  if (!FPT_COMP_E(Ave.Var(),4.0))
  {
    std::cout << "Ave.Var() = " << Ave.Var() << std::endl;
    throw 12;
  }

  if (!FPT_COMP_E(Ave.VarA(),56.0/3-144.0/9))
  {
    std::cout << "Ave.VarA() = " << Ave.VarA() << std::endl;
    throw 13;
  }

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
