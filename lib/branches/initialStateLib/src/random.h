//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/random.h $
//$LastChangedDate: 2014-11-27 20:38:22 +0100 (Do, 27. Nov 2014) $
//$LastChangedRevision: 1977 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


#ifndef RANDOM_NG_H
#define RANDOM_NG_H

#include <stdint.h>
#include <vector>
#include <boost/random.hpp>
#if _OPENMP
#include <omp.h>
#endif

#include "configBAMPS.h"
#include "dSFMT.hpp"

#ifndef USE_RNG
#define USE_RNG 1
#endif

/** 
 * @brief Abstract class for random number generators
 *
 * This class provides a consistent interface to possibly different
 * random number generators (RNG) and also handles the setting of
 * seeds for these generators. 
 */
class randomNumberGenerator
{
public:
  /** 
   * @brief Standard constructor 
   */
  randomNumberGenerator() : 
    mySeed(),
    myName("Base class")
  {};

  /** 
   * @brief constructor 
   */
  randomNumberGenerator( std::string name ) : 
    mySeed(),
    myName(name)
  {};

  /** 
   * @brief Compute a nice seed
   **/
  uint32_t findSeed();

  /** 
   * @brief Set seed provided by user 
   **/
  uint32_t setSeed( const uint32_t s )  
  { 
    return mySeed = s; 
  }
  
  /** 
   * @brief Get random number
   * @return random floating point number uniformly distributed in
   * (0,1)
   *
   * The endpoints are excludes, i.e. we have an open-open intervall
   */
  //    double getRan() = 0;
  
  /**
   * @brief Get random number via overloaded ()-operator
   * @return random floating point number uniformly distributed in
   * (0,1)
   *
   * The endpoints are excludes, i.e. we have an open-open intervall
   */
  //    double operator()() = 0;

  /**
   * @brief return the name (and other info) as string
   **/
  std::string Name(void) const { return myName; };
  
protected:
  /** @brief Stores the seed for internal use */
  uint32_t mySeed;

  /** @brief Stores the Name and other info */
  std::string myName;
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** 
 * @brief This class encapsulates the BOOST generation of random number
 *
 * This class provides a consistent interface to possibly different
 * random number generators (RNG) and also handles the setting of
 * seeds for these generators. 
 */
class randomNumberGeneratorBoost : public randomNumberGenerator
{
public:
  /** 
   * @brief Standard constructor 
   *
   * This standard constructor initializes the boost::uniform_01
   * object with the RNG and sets a standard seed.
   */
  randomNumberGeneratorBoost() : 
    randomNumberGenerator("Boost"),
    ran( rng ) 
  { 
    setSeed( findSeed() ); 
  }
  
  /**
   * @brief Constructor that already sets the seed
   *
   * This constructor initializes the boost::uniform_01 object with
   * the RNG and sets a user provided seed. 
   */
  randomNumberGeneratorBoost( const uint32_t s ) : 
    randomNumberGenerator("Boost"),
    ran( rng ) 
  { 
    setSeed( s ); 
  }
  
  /** 
   * @brief Set seed provided by user 
   **/
  uint32_t setSeed( const uint32_t s )  
  { 
    rng.seed( s );  
    return randomNumberGenerator::setSeed( s );
  }
  
  /** 
   * @brief Get random number
   * @return random floating point number uniformly distributed in (0,1) (i.e. endpoints excluded)
   */
  double getRan() { return this->operator()(); }
  
  /**
   * @brief Get random number via overloaded ()-operator
   * @return random floating point number uniformly distributed in (0,1) (i.e. endpoints excluded)
   */
  double operator()()
  {
    double r = ran();
    while (r == 0.0) r = ran();
    return r;
  }

  /** 
   * @brief Get integer valued random number
   * @param[in] n the maximal value
   * @return an integer in the range 0...(n-1)
   **/
  unsigned int Int(const unsigned int n)
  {
    boost::random::uniform_int_distribution<> dist(0, n-1);
    return dist( rng );

    // unsigned int r = n * ran();
    // while (r == n) r = n * ran();
    // return r;
  }
  
private:
  
  /**
   * @brief The type of the RNG to use
   **/
  typedef boost::mt19937 tRNG;
  //  typedef  boost::rand48 tRNG;

  /**
   * @brief The type of the distribution
   **/
  typedef boost::uniform_01<tRNG&, double> tDist;

  /**
   * @brief A random number generator (RNG) from the Boost libraries 
   *
   * The Boost RNG is used to obtain random numbers (integer type). It
   * must be seeded prior to use, c.f. setSeed(). 
   * boost::mt19937 uses a Mersenne twister algorithm.
   * Other reasonable choices for the RNG would be
   * boost::lagged_fibonacci607 or boost::kreutzer1986. 
   * See http://www.boost.org/doc/libs/release/libs/random
   */
  tRNG rng;
  
  /** 
   * @brief A distribution generator from the Boost libraries 
   *
   * This maps the random numbers from a Boost RNG to uniformly
   * distributed floating point numbers in the interval [0,1). 
   * 
   * The pointer needs to be initialized in the constructors with the
   * actual RNG. 
   */
  tDist ran;
  
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if _OPENMP

/** 
 * @brief This class encapsulates the BOOST generation of random
 * number while using OpenMP
 *
 * This class provides a consistent interface to possibly different
 * random number generators (RNG) and also handles the setting of
 * seeds for these generators. 
 *
 * If we are using openMP, we define multiple RNG
 */
class randomNumberGeneratorBoost_OpenMP : public randomNumberGenerator
{
private:
  
  /**
   * @brief The type of the RNG to use
   **/
  typedef boost::mt19937 tRNG;
  //  typedef  boost::rand48 tRNG;

  /**
   * @brief The type of the distribution
   **/
  typedef boost::uniform_01<tRNG&, double> tDist;

  /**
   * @brief Array of random number generator (RNG) from the Boost libraries 
   **/
  std::vector<tRNG> rng;
  
  /** 
   * @brief Array of distribution generator from the Boost libraries 
   */
  std::vector<tDist> ran;

  /**
   * @brief maximal number of threads
   **/
  unsigned int nThread;

public:
  /** 
   * @brief Standard constructor 
   *
   * This standard constructor initializes the boost::uniform_01
   * object with the RNG and sets a standard seed.
   */
  randomNumberGeneratorBoost_OpenMP() : 
    randomNumberGenerator("Boost_OpenMP")
  { 
    nThread = omp_get_max_threads();
    rng.resize( nThread );
    for (unsigned int i=0;i<nThread;i++)
    {
      tDist* pRan = new tDist( rng[i] );
      ran.push_back( *pRan );
    }
    setSeed( findSeed() ); 
  }
  
  /**
   * @brief Constructor that already sets the seed
   *
   * This constructor initializes the boost::uniform_01 object with
   * the RNG and sets a user provided seed. 
   **/
  randomNumberGeneratorBoost_OpenMP( const uint32_t s ) : 
    randomNumberGenerator("Boost_OpenMP")
  { 
    nThread = omp_get_max_threads();
    rng.resize( nThread );
    for (unsigned int i=0;i<nThread;i++)
    {
      tDist* pRan = new tDist( rng[i] );
      ran.push_back( *pRan );
    }
    setSeed( s ); 
  }

  /** 
   * @brief Set seed provided by user 
   **/
  uint32_t setSeed( const uint32_t s )  
  { 
    std::cout << "set seed: " << s << "  nThreads=" << nThread << std::endl;
    for (unsigned int i=0; i<nThread; i++)
    {
      rng[i].seed( s + i );  
    }
    return randomNumberGenerator::setSeed( s );
  }
  
  /** 
   * @brief Get random number
   * @return random floating point number uniformly distributed in (0,1) (i.e. endpoints excluded)
   */
  double getRan() { return this->operator()(); }
  
  /**
   * @brief Get random number via overloaded ()-operator
   * @return random floating point number uniformly distributed in (0,1) (i.e. endpoints excluded)
   */
  double operator()()
  {
    unsigned int iThread = omp_get_thread_num();
    double r = ran[iThread]();
    while (r == 0.0) r = ran[iThread]();
    return r;
  }

  /** 
   * @brief Get integer valued random number
   * @param[in] n the maximal value
   * @return an integer in the range 0...(n-1)
   **/
  unsigned int Int(const unsigned int n)
  {
    boost::random::uniform_int_distribution<> dist(0, n-1);
    unsigned int iThread = omp_get_thread_num();
    return dist( rng[iThread] );
  }
  
private:


};
#endif // _OPENMP

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


/** 
 * @brief This class encapsulates the dSFMT generation of random number
 *
 * This class provides a consistent interface to possibly different
 * random number generators (RNG) and also handles the setting of
 * seeds for these generators. 
 */
class randomNumberGeneratordSFMT : public randomNumberGenerator
{
public:
  /** 
   * @brief Standard constructor 
   *
   * This standard constructor initializes the object and sets a
   * standard seed. 
   */
  randomNumberGeneratordSFMT() : 
    randomNumberGenerator("dSFMT"),
    rng( dsfmt_19937_generator() )
  {  
    setSeed( findSeed() ); 
  }
  
  /**
   * @brief Constructor that already sets the seed
   *
   * This constructor initializes the object and sets a user provided
   * seed. 
   */
  randomNumberGeneratordSFMT( const uint32_t s ) : 
    randomNumberGenerator("dSFMT"),
    rng( dsfmt_19937_generator() )
  { 
    setSeed( s ); 
  }
  
  
  /** 
   * @brief Set seed provided by user 
   **/
  uint32_t setSeed( const uint32_t s )  
  { 
    rng.seed( s );  
    return randomNumberGenerator::setSeed( s );
  }
  
  /** 
   * @brief Get random number
   * @return random floating point number uniformly distributed in [0,1)
   */
  double getRan() { return this->operator()(); }
  
  /**
   * @brief Get random number via overloaded ()-operator
   * @return random floating point number uniformly distributed in [0,1)
   */
  double operator()() { return rng.rand_close_open(); };

  /** 
   * @brief Get integer valued random number
   * @param[in] n the maximal value
   * @return an integer in the range 0...(n-1)
   **/
  unsigned int Int(const unsigned int n)
  {
    unsigned int r = n * getRan();
    while (r == n) r = n * getRan();
    return r;
  }
  
private:
  /** 
   * @brief the actual implementation of the dSMFT generator
   **/
  dsfmt_19937_generator rng;
  
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


/** @brief Random number generator (RNG) object, GLOBAL 
*
* Declared as extern here. There MUST be ONE definition in a cpp-file that includes this header in the global 
* scope (it's in random.cpp).
* All files that include this header have then access to the random number generator.
*/
#if (USE_RNG == 1)
#if _OPENMP
extern randomNumberGeneratorBoost_OpenMP ran2;
#else
extern randomNumberGeneratorBoost ran2;
#endif
#elif (USE_RNG == 2)
#if _OPENMP
#error "USE_RNG == 2 && _OPENMP is not implemented yet."
#else
extern randomNumberGeneratordSFMT ran2;
#endif
#else
#error "no RNG selected: set USE_RNG correctly"
#endif

#endif
