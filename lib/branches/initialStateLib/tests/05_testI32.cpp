//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/tests/05_testI32.cpp $
//$LastChangedDate: 2014-12-10 22:53:55 +0100 (Mi, 10. Dez 2014) $
//$LastChangedRevision: 2011 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/*
  This is some simple executable for testing.
  Here we test the possibilities of the 3->2 integration
*/

#include <iostream>
#include <fstream>
#include <string>

#include "configurationbase.h"
#include "random.h"
#include "scattering32.h"
#include "particleprototype.h"
#include "TimerClock.h"

using namespace std;

configBase theConfig;
scattering32 scatt32_object;

double I32_input, I32_output;
TimerClock timerClock;

int SetParameters(std::string line)
{
  double I32;
  double P1[4],P2[4],P3[4];
  int F1,F2,F3;
  double sqrts, md2g_scaled, lambda_scaled,_alpha_s;
  bool _md2_use_simplified_GB, _matrixElement32_22qt;
  int _Ng;
  
  std::istringstream iss(line);
  iss.clear();

  iss >> I32 
      >> P1[0] >> P1[1] >> P1[2] >> P1[3] 
      >> P2[0] >> P2[1] >> P2[2] >> P2[3] 
      >> P3[0] >> P3[1] >> P3[2] >> P3[3]
      >> F1 >> F2 >> F3
      >> sqrts >> md2g_scaled >> lambda_scaled >> _alpha_s
      >> _md2_use_simplified_GB >> _matrixElement32_22qt
      >> _Ng; 

  if (iss.fail()) return 0; // Zero, if failure



  I32_input = I32;
  VectorXYZ v(0,0,0);
  VectorEPxPyPz P1v(P1[0], P1[1],P1[2],P1[3] );
  VectorEPxPyPz P2v(P2[0], P2[1],P2[2],P2[3] );
  VectorEPxPyPz P3v(P3[0], P3[1],P3[2],P3[3] );
  

  scatt32_object.setParameter( v, P1v,P2v,P3v, 
                               static_cast<FLAVOR_TYPE>( F1 ),
                               static_cast<FLAVOR_TYPE>( F2 ),
                               static_cast<FLAVOR_TYPE>( F3 ),
                               sqrts, md2g_scaled, lambda_scaled,_alpha_s,
                               _md2_use_simplified_GB, _matrixElement32_22qt,
                               _Ng);

  return 1; // One, if success
  
}


int main(int argc, char **argv)
{
  try
  {
    // Reading some input options:

    theConfig.readAndProcessProgramOptions( argc, argv );

    // Setting the seed of the random number generator:

    uint32_t seed;
    
    if ( (seed = theConfig.getSeed()) == 0 )
    {
      seed = ran2.findSeed();
    }
    ran2.setSeed( seed );
    cout << "seed: " << seed << endl;

    // Reading the actual calculational parameters:
    string line;
    ifstream myfile ("testNf3_32configurations.dat");
    if (!(myfile.is_open())) throw( 11 );

    int iLine = 0;

    while ( myfile.good() )
    {
      getline(myfile,line);
      if (line[0] != '#')
      {
        if (SetParameters(line))
        {
          iLine++;

          int initialStateIndex;
          I32_COMPUTATION_TYPE I32compType = MONTE_CARLO_INTEGRATION;
          //          I32_COMPUTATION_TYPE I32compType = FAST;

          timerClock.Start();

          I32_output = scatt32_object.getIntegral32_withPrefactors( initialStateIndex, I32compType );

          timerClock.Stop();

          cout << iLine << " "
               << I32_input << " "
               << I32_output << " "
               << timerClock.Nanosecs()
               << endl;

        }
      }
    }
    myfile.close();
  }
  catch (int e)
  {
    cout << "An exception occurred. Exception Nr. " << e << endl;
    return e;
  }
  return 0;  
}
