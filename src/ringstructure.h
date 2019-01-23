//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/ringstructure.h $
//$LastChangedDate: 2018-04-28 19:57:48 +0200 (Sa, 28. Apr 2018) $
//$LastChangedRevision: 2742 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

// at revision 902, this file is identical to full/branches/vector4D/src/ringstructure.h
//
// the method 'setLongitudinalGeometry' is new !!!

#ifndef RINGSTRUCTURE_H
#define RINGSTRUCTURE_H

#include <vector>
#include <string>
#include <stdexcept>

#include "bampsvector.h"
#include "particle.h"
#include "ratesmanager.h"
#include "ringcontainer.h"


/**
 * @brief class to store information in multiple rings
 **/
class ringStructure
{
public:
  /** 
   * @brief Constructor
   */
  ringStructure() : 
    centralRadius( 0 ), 
    totalRadius( 0 ), 
    deltaR( 0 ),
    delta_z(0), 
    y_left(0), 
    y_right(0), 
    timenow(0)
  { 
    rings.resize(0); 
  }
  
  /** 
   * @brief Constructor
   */
  ringStructure( const int _nRings, const double _centralRadius, const double _deltaR );

  /** 
   * @brief Destructor
   */
  ~ringStructure() {};
  
  /**
   * @brief Reset the configuration of the rings
   **/
  void resize( const int _nRings, const double _centralRadius, const double _deltaR );
  
  /**
   * @brief return the ring corresponding to given radius
   **/
  ringContainer& getRing( const double _xt ) { return rings[ getIndex( _xt ) ]; }

  /**
   * @brief return the ring given by index
   **/
  const ringContainer& operator[]( const unsigned int _index ) const;

  /**
   * @brief add some information about the longitudinal position
   **/
  void setLongitudinalGeometry( const double _y_left, const double _y_right, const double _t );

  /**
   * @brief calculate the index corresponding to the given radius
   **/
  unsigned int getIndex( const double _xt ) const 
  {
    return std::min( getIndexPure( _xt ), static_cast<unsigned int>(rings.size() - 1) );
  };

  /**
   * @brief calculate the index corresponding to the given radius, no
   * checks
   *
   * This routine does not check, whether the index may be too large
   **/
  unsigned int getIndexPure( const double _xt ) const
  {
    return ( _xt > centralRadius ) ? static_cast<unsigned int>( ( _xt - centralRadius ) / deltaR ) + 1 : 0;    
  };

  /**
   * @brief calculate the index the given particle belongs to
   **/
  unsigned int getIndex( const Particle& _particle ) const { return getIndex( _particle.Pos.Perp() ); };

  /**
   * @brief calculate the index the given particle belongs to, no
   * checks 
   *
   * This routine does not check, whether the index may be too large
   **/
  unsigned int getIndexPure( const Particle& _particle ) const { return getIndexPure( _particle.Pos.Perp() ); };

  /**
   * @brief return the number of rings
   **/
  unsigned int size() const { return rings.size(); }

  /**
   * @brief return the value of the cantral ring radius 
   **/
  double getCentralRadius() const { return centralRadius; }

  /**
   * @brief return the value of the ring width
   **/
  double getDeltaR() const { return deltaR; }
  
  /**
   * @brief reset all rings to its default values
   **/
  void clear();  
  
  /**
   * @brief add information of some particle
   **/
  void addParticle( const double _xt, const Particle& _particle );

  /**
   * @brief add information of some particle
   **/
  void addParticle( const Particle& _particle );

  /**
   * @brief add information of some particle, which is in formation
   * after a geometrical collision
   **/
  void addParticleInFormGeom( const double _xt, const Particle& _particle, const double _time );

  /**
   * @brief add information of some particle, which is in formation
   * after a geometrical collision
   **/
  void addParticleInFormGeom( const Particle& _particle, const double _time );

  /**
   * @brief add the rates given by the particle
   **/
  void addRates( const double _xt, const Particle& _particle );

  /**
   * @brief add the rates given by the particle
   **/
  void addRates( const Particle& _particle );
  
  /** 
   * @brief Do the actual calculation of the averaged values in all
   * rings
   **/
  void prepareAverages( const double _dz, const int _Ntest );
 
  void setValuesAsWorkaround(const int ringIndex, const double Edens,const double Pdens,const double _gamma, const int _numberOfParticles,const double Gdens,const double Qdens);
  
  
  /**
   * @brief give some value for md2g
   *
   * This routine first looks at the stored value at the given
   * transverse position, or looks at the smaller rings for the first
   * non-zero value.
   **/
  double find_md2g( const Particle& _particle ) const;

  /**
   * @brief give some value for md2q
   *
   * This routine first looks at the stored value at the given
   * transverse position, or looks at the smaller rings for the first
   * non-zero value.
   **/
  double find_md2q( const Particle& _particle ) const;

  /**
   * @brief set rate22v of particle from info stored in rings
   **/
  void getRateV( Particle * pPart ) const;

  /**
   * @brief set rate22 and md2g of particle from info stored in rings
   **/
  void getRate( Particle * pPart ) const;

  /**
   * @brief return the value for the energy density
   *
   * if the transverse radius is larger than the stored ring
   * structure, a value of zero is returned.
   **/
  double getEnergyDensity( const VectorEPxPyPz & Pos ) const;

  /**
   * @brief return the value for gamma
   *
   * if the transverse radius is larger than the stored ring
   * structure, a value of zero is returned.
   **/
  double getGamma( const VectorEPxPyPz & Pos ) const
  {
    return rings[ getIndex( Pos.Perp() ) ].getGamma();
  }


private:
  /**
   * @brief Set up all rings
   **/
  void initialize(void);
  
  std::vector<ringContainer> rings; ///< the storage of the rings
  
  double centralRadius; ///< the central radius
  double totalRadius; ///< the maximal radius
  double deltaR; ///< the ring thickness

  double delta_z; ///< longitudinal information
  double y_left; ///< longitudinal information
  double y_right; ///< longitudinal information
  double timenow; ///< the actual time
};




/** @brief exception class */
class eRingStructure_error : public std::runtime_error
{
public:
  explicit eRingStructure_error(const std::string& what) : std::runtime_error(what) {};
  
  virtual ~eRingStructure_error() throw() {};
};


#endif // RINGSTRUCTURE_H
// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on; 
