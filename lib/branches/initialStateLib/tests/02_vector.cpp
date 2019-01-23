//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/tests/02_vector.cpp $
//$LastChangedDate: 2015-10-29 22:08:25 +0100 (Do, 29. Okt 2015) $
//$LastChangedRevision: 2228 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/*
  This is some simple executable for testing.
  Here we test the basic features of vectors.
*/

#include <iostream>
using namespace std;

#include "FPT_compare.h"
#include "lorentz.h"

#include "bampsvector.h"


// Some constants:

const double A[4] = { 1,2,3,4 };
const double B[4] = { 4,3,2,1 };


void test1(void)
{
  VectorEPxPyPz v0;
  double Res;

  v0.SetTXYZ( A );

  Res = v0.M2();
  if (!FPT_COMP_E(Res, -28.)) throw 101;

  Res = v0.vec2();
  if (!FPT_COMP_E(Res, 29.)) throw 102;
  
}

void test2(void)
{
  VectorEPxPyPz v0,v1,v2;
  v1.SetTXYZ( A );
  v2.SetTXYZ( B );
  v0 = min( v1,v2-2.0 );
  for (unsigned int i = 0; i<4; ++i)
  {
    if (v0(i) != min( v1(i),v2(i)-2.0 )) throw 201;
  }
}

void test3(void)
{
  VectorEPxPyPz v0,v1,v2;
  v1.SetTXYZ( A );
  v2.SetTXYZ( B );

  v0 = v1+v2;
  for (unsigned int i = 0; i<4; ++i)
  {
    if ( v0(i) != v1(i)+v2(i) ) throw 3100+i;
  }
  v0 = v1;
  v0 += v2;
  for (unsigned int i = 0; i<4; ++i)
  {
    if ( v0(i) != v1(i)+v2(i) ) throw 3110+i;
  }


  
  v0 = v1-v2;
  for (unsigned int i = 0; i<4; ++i)
  {
    if ( v0(i) != v1(i)-v2(i) ) throw 3200+i;
  }
  v0 = v1;
  v0 -= v2;
  for (unsigned int i = 0; i<4; ++i)
  {
    if ( v0(i) != v1(i)-v2(i) ) throw 3210+i;
  }

  

}


int main(int argc, char **argv)
{

  try
  {
    test1();
    test2();
    test3();
  }
  catch (int e)
  {
    cout << "An exception occurred. Exception Nr. " << e << endl;
    return e;
  }
  return 0;  
}
