//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/tests/99_centralityclass.cpp $
//$LastChangedDate: 2014-12-10 22:53:55 +0100 (Mi, 10. Dez 2014) $
//$LastChangedRevision: 2011 $
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

/**
 * @brief The function doing the actual work
 *
 * Here we throw an integer in order to indicate, which internal test
 * failed.
 **/

void test1(void)
{

  ///
  /// Here you have to add your own code...
  ///

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
