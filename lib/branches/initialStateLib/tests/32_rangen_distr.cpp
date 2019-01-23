//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/tests/32_rangen_distr.cpp $
//$LastChangedDate: 2015-03-15 21:58:26 +0100 (So, 15. Mär 2015) $
//$LastChangedRevision: 2133 $
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

#include <iomanip>
#include <iostream>
//using namespace std;

#include "FPT_compare.h"
#include "rangen_distr.h"
#include "tools.h"

/**
 * @brief The function doing the actual work
 *
 * Here we throw an integer in order to indicate, which internal test
 * failed.
 **/

void test1(void)
{
  ran2.setSeed( 13 );

  std::vector<double> x = create_vector<double>(0.0)(1.0)(2.0)(3.0)(4.0);
  std::vector<double> y = create_vector<double>(0.0)(0.1)(0.3)(0.4)(1.0);

  ranGen_Distr rDist( x, y, interp_cspline );

  // std::cout << std::setprecision(12) << rDist() << std::endl;
  // std::cout << std::setprecision(12) << rDist() << std::endl;
  // std::cout << std::setprecision(12) << rDist() << std::endl;
  // std::cout << std::setprecision(12) << rDist() << std::endl;

  if (!FPT_COMP_E( rDist(), 4.42787837579 ) ) throw 111;
  if (!FPT_COMP_E( rDist(), 4.28010064219 ) ) throw 112;
  if (!FPT_COMP_E( rDist(), 1.63342050911 ) ) throw 113;
  if (!FPT_COMP_E( rDist(), 4.31876801830 ) ) throw 114;

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
