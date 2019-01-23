//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/tests/01_lorentz.cpp $
//$LastChangedDate: 2015-02-26 17:06:51 +0100 (Do, 26. Feb 2015) $
//$LastChangedRevision: 2093 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/*
  This is some simple executable for testing.
  Here we test the possibilities to boost some vectors.
*/

#include <iostream>
using namespace std;

#include "FPT_compare.h"
#include "lorentz.h"

#include "bampsvector.h"

// Some constants:

const double A[4] = { 1,2,3,4 };

const double beta1[4] = { 0,0.1,0,0 };
const double beta2[4] = { 7,0.1,0,0 };
const double beta3[4] = { 0,0.8,0.1,0.2 };

const double LA1[4] = {0.804030252207,1.90957184899,3,4};
const double LA2[4] = {0.804030252207,1.90957184899,3,4};
const double LA3[4] = {-3.05329013446,3.05514964723,3.1318937059,4.26378741181};

// void test1_old(void)
// {
//   double LA[4];
  
//   //  cout << setprecision(12) << "LA: " << LA[0] << "," << LA[1] << "," << LA[2] << "," << LA[3] << endl;


//   lorentz(beta1,A,LA);
//   for (int i=0;i<4;i++)
//   {
//     if (!FPT_COMP_E(LA[i],LA1[i])) throw 10+i;
//   }

//   lorentz(beta2,A,LA);
//   for (int i=0;i<4;i++)
//   {
//     if (!FPT_COMP_E(LA[i],LA2[i])) throw 20+i;
//   }

//   lorentz(beta3,A,LA);
//   for (int i=0;i<4;i++)
//   {
//     if (!FPT_COMP_E(LA[i],LA3[i])) throw 30+i;
//   }

// }

void test1_new(void)
{
  VectorEPxPyPz v0,v1,h1a,h1b,h1c,h2a,h2b,h2c;
  VectorTXYZ beta;
  lorentz LL;

  cout << "aligned: " << is_alignedSSE( &beta )
       << " " << is_alignedSSE( &LL ) 
       << "     " << VectorAlignment << endl;
  v0.SetTXYZ( A );
  //  cout << "v0: " << v0 << endl ;
  beta.SetTXYZ( beta1 );
  //  cout << "beta1: " << beta << endl ;
  LL.setBeta( beta );
  v1 = LL * v0;

  //  cout << "LL: " << LL << endl;
  //  cout << "v1: " << v1 << endl ;

  for (int i=0;i<4;i++)
  {
    if (!FPT_COMP_E(v1(i),LA1[i])) 
    {
      cout << "LL: " << LL << endl;
      cout << "v1: " << v1 << endl ;
      throw 110+i;
    }
  }

  LL.boost(v0,v0, h1a,h1b);
  for (int i=0;i<4;i++)
  {
    if (!FPT_COMP_E(h1a(i),LA1[i])) throw 120+i;
    if (!FPT_COMP_E(h1b(i),LA1[i])) throw 130+i;
  }
  
  beta.SetTXYZ( beta3 );
  //  cout << "beta3: " << beta << endl ;
  LL.setBeta( beta );
  v1 = LL * v0;

  //  cout << "v1: " << v1 << endl ;

  for (int i=0;i<4;i++)
  {
    if (!FPT_COMP_E(v1(i),LA3[i])) 
    {
      cout << "LL: " << LL << endl;
      cout << "v1: " << v1 << endl ;
      throw 310+i;
    }
  }

  LL.boost(v0, h1a);
  LL.boostInv(h1a,h2a);
  for (int i=0;i<4;i++)
  {
    if (!FPT_COMP_E(v0(i),h2a(i))) throw 410+i;
  }

  LL.boost(v0,v0, h1a,h1b);
  LL.boostInv(h1a,h1b, h2a,h2b);
  for (int i=0;i<4;i++)
  {
    if (!FPT_COMP_E(v0(i),h2a(i))) throw 510+i;
    if (!FPT_COMP_E(v0(i),h2b(i))) throw 520+i;
  }

  LL.boost(v0,v0,v0, h1a,h1b,h1c);
  LL.boostInv(h1a,h1b,h1c, h2a,h2b,h2c);
  for (int i=0;i<4;i++)
  {
    if (!FPT_COMP_E(v0(i),h2a(i))) throw 610+i;
    if (!FPT_COMP_E(v0(i),h2b(i))) throw 620+i;
    if (!FPT_COMP_E(v0(i),h2c(i))) throw 630+i;
  }

}

void test2_new(void)
{
  VectorEPxPyPz P1,P2;

  //  BetaVector beta;

  //  beta = BoostToCM( P1 );
  

    

}


int main(int argc, char **argv)
{

  try
  {
    //    test1_old();
    test1_new();
  }
  catch (int e)
  {
    cout << "An exception occurred. Exception Nr. " << e << endl;
    return e;
  }
  return 0;  
}
