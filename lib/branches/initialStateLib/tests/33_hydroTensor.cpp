//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/tests/33_hydroTensor.cpp $
//$LastChangedDate: 2015-10-29 22:06:02 +0100 (Do, 29. Okt 2015) $
//$LastChangedRevision: 2227 $
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
#include "hydroTensor.h"
#include "tools.h"

/**
 * @brief The function doing the actual work
 *
 * Here we throw an integer in order to indicate, which internal test
 * failed.
 **/

void test1(void)
{
  // numbers are calculated with Mathematica
  
  //  std::cout << std::setprecision(10) << tTmunuNmu::Boltzmann_nDens( 0.200, 0.0, 16) << std::endl;

  if (!FPT_COMP_E( tTmunuNmu::Boltzmann_nDens( 0.200, 0.0, 16), 1.696334698937415 )) throw 111;

  if (!FPT_COMP_E( tTmunuNmu::Boltzmann_nDens( 0.200, 0.138, 16), 1.5260034818492125 )) throw 112;

  

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
