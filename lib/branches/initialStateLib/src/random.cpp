//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/random.cpp $
//$LastChangedDate: 2014-09-25 21:24:05 +0200 (Do, 25. Sep 2014) $
//$LastChangedRevision: 1884 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <time.h>
#include <boost/random.hpp>

#include "random.h"

using std::cout;
using std::endl;


/** 
 * @brief definition of global RNG object, defined extern in random.h 
 **/
#if (USE_RNG == 1)
#if _OPENMP
randomNumberGeneratorBoost_OpenMP ran2;
#else
randomNumberGeneratorBoost ran2;
#endif
#elif (USE_RNG == 2)
#if _OPENMP
#error "USE_RNG == 2 && _OPENMP is not implemented yet."
#else
randomNumberGeneratordSFMT ran2;
#endif
#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/**
 * This routine computes and sets a (hopefully) unique seed from either
 * /dev/urandom provided by the system or the system time as a
 * fallback. 
 * Using /dev/urandom is prefered and should work on all *nix
 * system. /dev/urandom is a special device that holds (pseudo) random
 * bits generated from user input, network traffic and black
 * voodoo. Use /dev/random for even more unpredictable voodoo. 
 *
 * @return The computed seed.
 */
uint32_t randomNumberGenerator::findSeed()
{
  uint32_t s;
  
  // open /dev/urandom as an ifstream
  std::ifstream urand("/dev/urandom", std::ios::binary );
  
  // get the seed from /dev/urandom, fall back to a seed generated from the current time if /dev/urandom could not have been opened.
  if ( urand.is_open() && urand.good() )
  {
    // read random characters from the input stream assigned to /dev/urandom
    // for this the integer type s is reinterpreted as a character pointer and sizeof(unit32_t) characters are read into s
    urand.read( reinterpret_cast<char *>(&s), sizeof(uint32_t) );
  }
  else
  {
    s = static_cast<uint32_t>( time(NULL) );
  }
  
  // return the seed for the user to have fun with
  return s;
}



//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


