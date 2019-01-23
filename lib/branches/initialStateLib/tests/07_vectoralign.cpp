//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/tests/07_vectoralign.cpp $
//$LastChangedDate: 2014-12-10 22:53:55 +0100 (Mi, 10. Dez 2014) $
//$LastChangedRevision: 2011 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/*
  This is some simple executable for testing.
  Here we test the basic features of vectors.
*/

#include <iostream>
#include <vector>
using namespace std;

#include "FPT_compare.h"
#include "lorentz.h"

#include "allocator.h"
#include "bampsvector.h"

// Some constants:

const double A[4] = { 1,2,3,4 };
const double B[4] = { 4,3,2,1 };


void test1(void)
{
  cout << __PRETTY_FUNCTION__ << endl;

  VectorEPxPyPz v0,*p1,*p2;
  cout << " ..." << endl;

  p1 = new VectorEPxPyPz();
  cout << "align1: " << is_alignedSSE( p1 ) << endl;
  cout << " ..." << endl;

  double *dd;

  dd = new double( 1.0 );
  cout << "alignd: " << is_alignedSSE( dd ) << endl;
  cout << " ..." << endl;

  p2 = new VectorEPxPyPz();
  cout << "align2: " << is_alignedSSE( p2 ) << endl;

  p1->SetTXYZ( A );
  cout << "set 1" << endl;
  p2->SetTXYZ( B );
  cout << "set 2" << endl;


  v0 = min( (*p1),(*p2)-2.0 );
  cout << "v0" << endl;

  for (unsigned int i = 0; i<4; ++i)
  {
    if (v0(i) != min( (*p1)(i),(*p2)(i)-2.0 )) throw 101;
  }
}

void test2(void)
{
  cout << __PRETTY_FUNCTION__ << endl;
  std::vector<VectorEPxPyPz> vv;
  VectorEPxPyPz v0;

  v0.SetTXYZ( A );
  vv.push_back(v0);
  v0.SetTXYZ( B );
  vv.push_back(v0);

  v0 = min( vv[0],vv[1]-2.0 );
  cout << "v0" << endl;

  for (unsigned int i = 0; i<4; ++i)
  {
    if (v0(i) != min( (vv[0])(i),(vv[1])(i)-2.0 )) throw 201;
  }
}

int main(int argc, char **argv)
{

  try
  {
    test1();
    test2();
  }
  catch (int e)
  {
    cout << "An exception occurred. Exception Nr. " << e << endl;
    return e;
  }
  return 0;  
}
