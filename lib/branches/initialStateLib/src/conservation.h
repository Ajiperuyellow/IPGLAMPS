
//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/conservation.h $
//$LastChangedDate: 2015-06-22 19:22:38 +0200 (Mo, 22. Jun 2015) $
//$LastChangedRevision: 2181 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file 
 * @brief This file implements conservation checkers, as e.g. the
 * Energy-Momentum conservation checks
 *
 * You may think about additional chechers, as e.g. 'charge' or so...
 *
 **/

#ifndef CONSERVATION_H
#define CONSERVATION_H

#include <stdexcept>
#include <string>
#include <vector>

#include "bampsvector.h"
#include "particleprototype.h"

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** 
 * @brief exception class for handling unexpected critical behaviour
 * within conservation checks
 **/
class eConservation_error : public std::runtime_error
{
public:
  explicit eConservation_error(const std::string& what) : std::runtime_error(what) {};

  virtual ~eConservation_error() throw() {};
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** 
 * @brief class to implement energy and momentum conservation checks
 *
 * Is intended to be run at every time step to indicate possible
 * problems 
 *
 * \tparam T The type of the particles, i.e. the class the particle
 * vector is build of. Since we iterate over a std::vector<T>, so it
 * may be ParticlePrototype or some derived classes. 
 * 
 **/
template<typename T>
class tConservation4Mom
{
public:
  /** 
   * @brief Constructor
   **/
  tConservation4Mom( ):
    Mom()
  {};

  /**
   * @brief initialize or do the check
   *
   * If the stored values are invalid, a call to this routine
   * initializes the class. Any subsequent call of this routine
   * calculates the actual values and compares them with those stored
   * in the first call. Throws an exception, if values differ.
   *
   * By purpose we do not use here FPT_COMP_E():
   * * we want to have the ability to set the accuracy
   * * there positive and negative values are treated differently
   **/
  void check( const std::vector<T> & _particles )
  {
    VectorEPxPyPz _Mom;
    for (auto it = begin(_particles); it != end(_particles); ++it )
    {
      _Mom += it->Mom;
    }
    
    if ( Mom.E() == 0.0 ) // we have to initialize...
    {
      Mom = _Mom;
    }
    else // we have to compare...
    {
      static const double eps = 1e-9; // relative allowed error
      
      if ( (fabs( Mom.E() -_Mom.E() ) > std::max(eps,eps*fabs( Mom.E() )) )||
           (fabs( Mom.Px() -_Mom.Px() ) > std::max(eps,eps*fabs( Mom.Px() )) ) ||
           (fabs( Mom.Py() -_Mom.Py() ) > std::max(eps,eps*fabs( Mom.Py() )) ) ||
           (fabs( Mom.Pz() -_Mom.Pz() ) > std::max(eps,eps*fabs( Mom.Pz() )) ) )
      {
        std::cout << "actual value: " << _Mom << std::endl
                  << "stored value: " << Mom << std::endl
                  << "check 0: " 
                  << (fabs( Mom.E() -_Mom.E() ) > std::max(eps,eps*fabs( Mom.E() ))) 
                  << std::endl
                  << "check x: " 
                  << (fabs( Mom.Px() -_Mom.Px() ) > std::max(eps,eps*fabs( Mom.Px() ))) 
                  << std::endl
                  << "check y: " 
                  << (fabs( Mom.Py() -_Mom.Py() ) > std::max(eps,eps*fabs( Mom.Py() ))) 
                  << std::endl
                  << "check z: " 
                  << (fabs( Mom.Pz() -_Mom.Pz() ) > std::max(eps,eps*fabs( Mom.Pz() ))) 
                  << std::endl;

        throw eConservation_error( "Energy/Momentum conservation violated!" );
      }
    }
  };
    
protected:
  VectorEPxPyPz Mom; ///< the energy/momentum to conserve
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif
