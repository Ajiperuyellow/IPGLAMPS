//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/tests/06_outputformat.cpp $
//$LastChangedDate: 2014-12-10 22:53:55 +0100 (Mi, 10. Dez 2014) $
//$LastChangedRevision: 2011 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

// At the moment, this is not a 'real' test in the sense of the 'known
// answer tests', where we compare numerical outputs against their
// known answers, but just a compile time test, where we guarantee,
// that the used output modifiers work. The output has still to be
// compared by eyes. (One may imagine to redirect the output to
// strings and then to compare these strings --- to be done).

#include <iostream>
using namespace std;

#include "FPT_compare.h"
#include "lorentz.h"

#include "bampsvector.h"
#include "particleprototype.h"


// Some constants:

const double A[4] = { 1,2,3,4 };

void test1(void)
{
  VectorEPxPyPz v0;
  // double Res;

  v0.SetTXYZ( A );

  std::cout << " v0 = "                  << v0 << std::endl;
  std::cout << " v0 = " << plainvector   << v0 << std::endl;
  std::cout << " v0 = "                  << v0 << std::endl;
  std::cout << " v0 = " << bracketvector << v0 << std::endl;
  //  std::cout << " v0 = " << set_val(4) << v0 << std::endl;
  std::cout << " v0 = "                  << v0 << std::endl;
  

  // Res = v0.M2();
  // if (!FPT_COMP_E(Res, -28.)) throw 101;

  // Res = v0.vec2();
  // if (!FPT_COMP_E(Res, 29.)) throw 102;
  
}

void test2(void)
{
  ParticlePrototype part;
  std::cout << "part = " << part << std::endl;
  std::cout << "part = " << plainvector << part << std::endl;
  std::cout << "part = " << bracketvector << part << std::endl;
  std::cout << "part = " << particlestyle1 << part << std::endl;
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
