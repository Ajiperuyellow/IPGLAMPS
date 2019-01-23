//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/tests/00_cmake.cpp $
//$LastChangedDate: 2014-12-10 22:53:55 +0100 (Mi, 10. Dez 2014) $
//$LastChangedRevision: 2011 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/*
  This is some simple executable for testing, whether CMake/CTest works
  correctly.

  It is a little bit more complicated than necessary, but also shows,
  how one can implement mor complicated tests.

*/


#include <iostream>
using namespace std;

void test1(void)
{
  // do nothing
}

int main(int argc, char **argv)
{

  try
  {
    test1();
  }
  catch (int e)
  {
    cout << "An exception occurred. Exception Nr. " << e << endl;
    return e;
  }
  return 0;  
}
