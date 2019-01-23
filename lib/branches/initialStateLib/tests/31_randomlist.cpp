//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/tests/31_randomlist.cpp $
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

#include <iostream>

#include "FPT_compare.h"
#include "random.h"
#include "randomlist.h"

/**
 * @brief The function doing the actual work
 *
 * Here we throw an integer in order to indicate, which internal test
 * failed.
 **/

void test1(void)
{
  ran2.setSeed( 13 );
  randomlist<int> rList(10,3);

  // std::cout << "# slices: " << rList.numberSlice() << std::endl;
  // for (unsigned int iSlice=0; iSlice < rList.numberSlice(); iSlice++)
  // {
  //   std::cout << "#" << iSlice << " : ";
  //   for (auto it=rList.beginSlice(iSlice); it!=rList.endSlice(iSlice); ++it) 
  //   {
  //     std::cout << (*it) << "  ";
  //   }
  //   std::cout << std::endl;
  // }

  if (rList.numberSlice() != 3) throw 100+1;
  if (*rList.beginSlice(0) != 2) throw 110+0;
  if (*rList.beginSlice(1) != 4) throw 110+1;
  if (*rList.beginSlice(2) != 8) throw 110+2;

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
