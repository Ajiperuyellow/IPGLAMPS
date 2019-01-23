//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/cascade.cpp $
//$LastChangedDate: 2018-12-21 20:20:59 +0100 (Fr, 21. Dez 2018) $
//$LastChangedRevision: 2913 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <exception>
#include <stdint.h>
#include <execinfo.h>
#include <signal.h>


#include "cascade.h"
#include "random.h"
#include "heavyIonCollision.h"
#include "mfp_data.h"
#include "global.h"
#include "signalhandler.h"

using std::cout;
using std::endl;
using std::vector;


/**
 * @brief The main routine of BAMPS
 */
int main(int argc, char *argv[])
{
  ns_casc::printSVN( cout );

  registerSignalHandler();

  
  // try-catch block to handle undesired program behaviour
  try
  { 
    //--------------------------------------------------------------
    // create and initialize the main objects needed for
    // configuration,execution and analysis of the simulation
    config theConfig;
    theConfig.readAndProcessProgramOptions( argc, argv );
    mfpForHeavyIonCollision theMFP( &theConfig );
    heavyIonCollision theHIC( &theConfig, &theMFP );
    analysis theAnalysis( &theConfig, &theHIC );
    //--------------------------------------------------------------
   
    //--------------------------------------------------------------
    // initialize the random number generator
    // defined globally in random.h and random.cpp

    uint32_t seed = theConfig.getSeed();
    if (seed == 0) seed = ran2.findSeed();
    ran2.setSeed( seed );
    theAnalysis.setSeed( seed );
    cout << "seed: " << seed << endl;

//     for(int i=1; i< 10;i++)
//     {
//       cout << ran2() << endl;
//     }
//     for(int i=1; i< 10;i++)
//     {
//       cout << ran2() << endl;
//     }
//     ran2.setSeed( seed );
//     for(int i=1; i< 10;i++)
//     {
//       cout << ran2() << endl;
//     }
//     
//     
//     
//     
//     exit(0);
    //--------------------------------------------------------------
    
    //--------cascade-------------------------------------
    cout << "=============start===============" << endl;
    
    theHIC.initialize( theAnalysis );
    theHIC.mainFramework( theAnalysis );
    
    cout << "==============end================" << endl;
    //--------end of cascade------------------------------
  }
  /**
  *
  * handle exceptions
  *
  */
  catch (int e) // provided to handle program termination after displaying usage and help messages
  {
    cout << "Signal caught: " << e << endl;
    if ( e == HELP_MESSAGE_REQUESTED )
    {
      return EXIT_SUCCESS;
    }
    else
    {
      return EXIT_FAILURE;
    }    
  }  
  catch (std::exception& e)
  {
    // output of the error information provided by the exception class
    cout << e.what() << endl; 
    cout << "Program terminated." << endl;
    return EXIT_FAILURE;
  }
  catch (...)
  {
    cout << "Unhandled exception. Program terminated." << endl;
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}
