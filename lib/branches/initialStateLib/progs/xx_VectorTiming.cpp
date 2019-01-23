//---------------------------------------------
//provided by subversion
//---------------------------------------------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/progs/xx_VectorTiming.cpp $
//$LastChangedDate: 2015-03-05 00:01:25 +0100 (Do, 05. Mär 2015) $
//$LastChangedRevision: 2100 $
//$LastChangedBy: gallmei $
//---------------------------------------------
//---------------------------------------------

/*
  Here we try to implement some timing routines for vector4D operations.

*/


#include <iostream>     // std::cout, std::fixed, std::scientific
#include <cstdlib>

#include "configBAMPS.h"
#include "TSC.h"
#include "bampsvector.h"
#include "lorentz.h"

using namespace std;

TimeStampCounter timer0;
TimeStampCounter timer1;


enum { nArr = 10000000 };

vector4D vec1[nArr];
vector4D vec2[nArr];
vector4D vec3[nArr];
long i1[nArr], i2[nArr], i3[nArr], i4[nArr];
long *ii1, *ii2, *ii3, *ii4;

vector4D sumVi0, sumVi1, sumVj0, sumVj1;
vector4D sumWi0, sumWi1, sumWj0, sumWj1;

double TimerDiff(const long nLoop)
{
  if (timer1.Cycles()>timer0.Cycles())
  {
    return (timer1.Cycles()-timer0.Cycles())*(1.0/nLoop);
  }
  else
  {
    return(timer0.Cycles()-timer1.Cycles())*(-1.0/nLoop);
  }
}

void Init(void)
{
  double maxval = 10.0;
  double fak = 1.0/RAND_MAX;

  vector4D sumV;

  for(long iArr=0;iArr < nArr; ++iArr)
  {
    double X=maxval*(2.0*std::rand()*fak-1.0);
    double Y=maxval*(2.0*std::rand()*fak-1.0);
    double Z=maxval*(2.0*std::rand()*fak-1.0);
    double E=maxval*std::rand()*fak;
    E=sqrt(E*E + X*X + Y*Y + Z*Z);
    vec1[iArr].SetTXYZ(E,X,Y,Z);

    sumV += vec1[iArr];

  }

  std::cout << "# Init:   sum = " << sumV*(1.0/nArr) << endl;

}



void test1(void)
{
  //  long nnLoop = 100;
  long nnLoop = 20;
  long nLoop = 1000000;
  double fak = (nArr*1.0)/RAND_MAX;

  lorentz LL0, LL1;

  Init();

  for(long iiLoop=0;iiLoop < nnLoop; ++iiLoop)
  {
    for(long iLoop=0;iLoop < nLoop; ++iLoop)
    {
      i1[iLoop] = (long)(std::rand()*fak);
      i2[iLoop] = (long)(std::rand()*fak);
      i3[iLoop] = (long)(std::rand()*fak);
      i4[iLoop] = (long)(std::rand()*fak);
    }
    
    ii1 = i1;
    ii2 = i2;
    ii3 = i3;
    ii4 = i4;

    timer0.Start();
    for(long iLoop=0;iLoop < nLoop; ++iLoop)
    {
      LL0.setBetaCMScalar(vec1[*ii1],vec1[*ii2]);
      //      LL0.setBetaCM(vec1[*ii1],vec1[*ii2]);

      //      LL0.boostScalar(vec1[*ii3], vec2[*ii1]);
      //      LL0.boostScalar(vec1[*ii4], vec3[*ii2]);

      LL0.boostScalar(vec1[*ii3], vec1[*ii4], vec2[*ii1], vec3[*ii2]);
      //      LL0.boostInvScalar(vec1[*ii3], vec1[*ii4], vec2[*ii1], vec3[*ii2]);


      ii1++;
      ii2++;
      ii3++;
      ii4++;

    }
    timer0.Stop();

    ii1 = i1;
    ii2 = i2;
    ii3 = i3;
    ii4 = i4;
    timer1.Start();
    for(long iLoop=0;iLoop < nLoop; ++iLoop)
    {
      //      LL0.setBetaCMScalar(vec1[*ii1],vec1[*ii2]);
      LL0.setBetaCM(vec1[*ii1],vec1[*ii2]);

      //      LL0.boostScalar(vec1[*ii3], vec2[*ii1]);
      //      LL0.boostScalar(vec1[*ii4], vec3[*ii2]);

      // LL0.boost(vec1[*ii3], vec2[*ii1]);
      // LL0.boost(vec1[*ii4], vec3[*ii2]);

      //      LL0.boostScalar(vec1[*ii3], vec1[*ii4], vec2[*ii1], vec3[*ii2]);

      LL0.boost(vec1[*ii3], vec1[*ii4], vec2[*ii1], vec3[*ii2]);
      //      LL0.boostInv(vec1[*ii3], vec1[*ii4], vec2[*ii1], vec3[*ii2]);

      ii1++;
      ii2++;
      ii3++;
      ii4++;
    }
    timer1.Stop();
    double timePerLoop = TimerDiff( nLoop );
    cout << "timer: " 
	 << timer0.Cycles() << " "
	 << timer1.Cycles() << " -> "
	 << timePerLoop  << endl;

    
    
  }

}

void testAll(void)
{
  //  long nnLoop = 100;
  long nnLoop = 20;
  long nLoop = 1000000;
  double fak = (nArr*1.0)/RAND_MAX;

  double times[10],timesA[10];

  lorentz LL0, LL1;

  Init();

  for (int i=0;i<10;i++) timesA[i] = 0.0;
  
  cout << "# results are given in CPU cycles (TSC)" << endl;
  cout << "#  0 = SIMD-Scalar, Total, 2 vectors" << endl;
  cout << "#  1 = SetBeta, SIMD" << endl;
  cout << "#  2 = SetBeta, Scalar" << endl;
  cout << "#  3 = SetBeta, SIMD (possible cache)" << endl;
  cout << "#  4 = SetBeta, Scalar (possible cache)" << endl;
  cout << "#  5 = Boost, SIMD (2 vectors)" << endl;
  cout << "#  6 = Boost, Scalar (2 vectors)" << endl;
  cout << "#  7 = Boost, SIMD (2 x 1 vector)" << endl;
  cout << "#  8 = Boost, Scalar (2 x 1 vector)" << endl;
  
  cout << setw(3) << "#  "
       << setw(12) << "    0   "
       << setw(12) << "    1   "
       << setw(12) << "    2   "
       << setw(12) << "    3   "
       << setw(12) << "    4   "
       << setw(12) << "    5   "
       << setw(12) << "    6   "
       << setw(12) << "    7   "
       << setw(12) << "    8   " << endl;

  for(long iiLoop=0;iiLoop < nnLoop; ++iiLoop)
  {
    for(long iLoop=0;iLoop < nLoop; ++iLoop)
    {
      i1[iLoop] = (long)(std::rand()*fak);
      i2[iLoop] = (long)(std::rand()*fak);
      i3[iLoop] = (long)(std::rand()*fak);
      i4[iLoop] = (long)(std::rand()*fak);
    }

    ////////////////// TOTAL

    ii1 = i1;    ii2 = i2;    ii3 = i3;    ii4 = i4;

    timer0.Start();
    for(long iLoop=0;iLoop < nLoop; ++iLoop)
    {
      LL0.setBetaCMScalar(vec1[*ii1],vec1[*ii2]);
      LL0.boostScalar(vec1[*ii3], vec1[*ii4], vec2[*ii1], vec3[*ii2]);

      ii1++;      ii2++;      ii3++;      ii4++;
    }
    timer0.Stop();

    ii1 = i1;    ii2 = i2;    ii3 = i3;    ii4 = i4;

    timer1.Start();
    for(long iLoop=0;iLoop < nLoop; ++iLoop)
    {
      LL0.setBetaCM(vec1[*ii1],vec1[*ii2]);
      LL0.boost(vec1[*ii3], vec1[*ii4], vec2[*ii1], vec3[*ii2]);


      ii1++;      ii2++;      ii3++;      ii4++;
    }
    timer1.Stop();

    times[0] = TimerDiff( nLoop );

    
    
    ////////////////// setBetaCM

    ii1 = i1;    ii2 = i2;    ii3 = i3;    ii4 = i4;

    timer0.Start();
    for(long iLoop=0;iLoop < nLoop; ++iLoop)
    {
      LL0.setBetaCMScalar(vec1[*ii1],vec1[*ii2]);

      ii1++;      ii2++;      ii3++;      ii4++;
    }
    timer0.Stop();

    ii1 = i1;    ii2 = i2;    ii3 = i3;    ii4 = i4;

    timer1.Start();
    for(long iLoop=0;iLoop < nLoop; ++iLoop)
    {
      LL0.setBetaCMScalar(vec1[*ii1],vec1[*ii2]);
      LL0.setBetaCM(vec1[*ii1],vec1[*ii2]);

      ii1++;      ii2++;      ii3++;      ii4++;
    }
    timer1.Stop();

    times[1] = TimerDiff( nLoop );

    /// ---------------------------------------


    ii1 = i1;    ii2 = i2;    ii3 = i3;    ii4 = i4;

    timer0.Start();
    for(long iLoop=0;iLoop < nLoop; ++iLoop)
    {
      LL0.setBetaCM(vec1[*ii1],vec1[*ii2]);

      ii1++;      ii2++;      ii3++;      ii4++;
    }
    timer0.Stop();

    ii1 = i1;    ii2 = i2;    ii3 = i3;    ii4 = i4;

    timer1.Start();
    for(long iLoop=0;iLoop < nLoop; ++iLoop)
    {
      LL0.setBetaCMScalar(vec1[*ii1],vec1[*ii2]);
      LL0.setBetaCM(vec1[*ii1],vec1[*ii2]);

      ii1++;      ii2++;      ii3++;      ii4++;
    }
    timer1.Stop();

    times[2] = TimerDiff( nLoop );

    /// ---------------------------------------

    ii1 = i1;    ii2 = i2;    ii3 = i3;    ii4 = i4;

    timer0.Start();
    for(long iLoop=0;iLoop < nLoop; ++iLoop)
    {
      LL0.setBetaCM(vec1[*ii1],vec1[*ii2]);

      ii1++;      ii2++;      ii3++;      ii4++;
    }
    timer0.Stop();

    ii1 = i1;    ii2 = i2;    ii3 = i3;    ii4 = i4;

    timer1.Start();
    for(long iLoop=0;iLoop < nLoop; ++iLoop)
    {
      LL0.setBetaCM(vec1[*ii1],vec1[*ii2]);
      LL0.setBetaCM(vec1[*ii1],vec1[*ii2]);

      ii1++;      ii2++;      ii3++;      ii4++;
    }
    timer1.Stop();

    times[3] = TimerDiff( nLoop );

    /// ---------------------------------------


    ii1 = i1;    ii2 = i2;    ii3 = i3;    ii4 = i4;

    timer0.Start();
    for(long iLoop=0;iLoop < nLoop; ++iLoop)
    {
      LL0.setBetaCMScalar(vec1[*ii1],vec1[*ii2]);

      ii1++;      ii2++;      ii3++;      ii4++;
    }
    timer0.Stop();

    ii1 = i1;    ii2 = i2;    ii3 = i3;    ii4 = i4;

    timer1.Start();
    for(long iLoop=0;iLoop < nLoop; ++iLoop)
    {
      LL0.setBetaCMScalar(vec1[*ii1],vec1[*ii2]);
      LL0.setBetaCMScalar(vec1[*ii1],vec1[*ii2]);

      ii1++;      ii2++;      ii3++;      ii4++;
    }
    timer1.Stop();

    times[4] = TimerDiff( nLoop );


    ////////////////// Boost2

    ii1 = i1;    ii2 = i2;    ii3 = i3;    ii4 = i4;

    timer0.Start();
    for(long iLoop=0;iLoop < nLoop; ++iLoop)
    {
      LL0.setBetaCM(vec1[*ii1],vec1[*ii2]);

      ii1++;      ii2++;      ii3++;      ii4++;
    }
    timer0.Stop();

    ii1 = i1;    ii2 = i2;    ii3 = i3;    ii4 = i4;

    timer1.Start();
    for(long iLoop=0;iLoop < nLoop; ++iLoop)
    {
      LL0.setBetaCM(vec1[*ii1],vec1[*ii2]);
      LL0.boost(vec1[*ii3], vec1[*ii4], vec2[*ii1], vec3[*ii2]);

      ii1++;      ii2++;      ii3++;      ii4++;
    }
    timer1.Stop();

    times[5] = TimerDiff( nLoop );

    /// ---------------------------------------


    ii1 = i1;    ii2 = i2;    ii3 = i3;    ii4 = i4;

    timer0.Start();
    for(long iLoop=0;iLoop < nLoop; ++iLoop)
    {
      LL0.setBetaCM(vec1[*ii1],vec1[*ii2]);

      ii1++;      ii2++;      ii3++;      ii4++;
    }
    timer0.Stop();

    ii1 = i1;    ii2 = i2;    ii3 = i3;    ii4 = i4;

    timer1.Start();
    for(long iLoop=0;iLoop < nLoop; ++iLoop)
    {
      LL0.setBetaCM(vec1[*ii1],vec1[*ii2]);
      LL0.boostScalar(vec1[*ii3], vec1[*ii4], vec2[*ii1], vec3[*ii2]);

      ii1++;      ii2++;      ii3++;      ii4++;
    }
    timer1.Stop();

    times[6] = TimerDiff( nLoop );

    ////////////////// Boost

    ii1 = i1;    ii2 = i2;    ii3 = i3;    ii4 = i4;

    timer0.Start();
    for(long iLoop=0;iLoop < nLoop; ++iLoop)
    {
      LL0.setBetaCM(vec1[*ii1],vec1[*ii2]);

      ii1++;      ii2++;      ii3++;      ii4++;
    }
    timer0.Stop();

    ii1 = i1;    ii2 = i2;    ii3 = i3;    ii4 = i4;

    timer1.Start();
    for(long iLoop=0;iLoop < nLoop; ++iLoop)
    {
      LL0.setBetaCM(vec1[*ii1],vec1[*ii2]);
      LL0.boost(vec1[*ii3], vec2[*ii1]);
      LL0.boost(vec1[*ii4], vec3[*ii2]);

      ii1++;      ii2++;      ii3++;      ii4++;
    }
    timer1.Stop();

    times[7] = TimerDiff( nLoop );

    /// ---------------------------------------


    ii1 = i1;    ii2 = i2;    ii3 = i3;    ii4 = i4;

    timer0.Start();
    for(long iLoop=0;iLoop < nLoop; ++iLoop)
    {
      LL0.setBetaCM(vec1[*ii1],vec1[*ii2]);

      ii1++;      ii2++;      ii3++;      ii4++;
    }
    timer0.Stop();

    ii1 = i1;    ii2 = i2;    ii3 = i3;    ii4 = i4;

    timer1.Start();
    for(long iLoop=0;iLoop < nLoop; ++iLoop)
    {
      LL0.setBetaCM(vec1[*ii1],vec1[*ii2]);
      LL0.boostScalar(vec1[*ii3], vec2[*ii1]);
      LL0.boostScalar(vec1[*ii4], vec3[*ii2]);

      ii1++;      ii2++;      ii3++;      ii4++;
    }
    timer1.Stop();

    times[8] = TimerDiff( nLoop );

    ////////////////// 

    for (int i=0;i<10;i++) timesA[i] += times[i];
    
    cout.precision(3);
    cout.width(12);
    cout << setw(3) << iiLoop << std::fixed
         << setw(12) << times[0] 
         << setw(12) << times[1] 
         << setw(12) << times[2] 
         << setw(12) << times[3] 
         << setw(12) << times[4] 
         << setw(12) << times[5] 
         << setw(12) << times[6] 
         << setw(12) << times[7] 
         << setw(12) << times[8] << endl;  

  }

  cout << endl
       << setw(3) << "Ave" << std::fixed
       << setw(12) << timesA[0]/nnLoop
       << setw(12) << timesA[1]/nnLoop
       << setw(12) << timesA[2]/nnLoop
       << setw(12) << timesA[3]/nnLoop
       << setw(12) << timesA[4]/nnLoop
       << setw(12) << timesA[5]/nnLoop
       << setw(12) << timesA[6]/nnLoop
       << setw(12) << timesA[7]/nnLoop
       << setw(12) << timesA[8]/nnLoop << endl;  

}


int main(int argc, char **argv)
{

  try
  {
    //    test1();
    testAll();
  }
  catch (int e)
  {
    cout << "An exception occurred. Exception Nr. " << e << endl;
    return e;
  }
  return 0;  
}
