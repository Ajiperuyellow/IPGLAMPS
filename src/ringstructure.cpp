//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/ringstructure.cpp $
//$LastChangedDate: 2018-04-28 19:57:48 +0200 (Sa, 28. Apr 2018) $
//$LastChangedRevision: 2742 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

// at revision 902, this file is identical to full/branches/vector4D/src/...
//
// the method 'setLongitudinalGeometry' is new!!!
//
// The original version was written for 'ParticleOffline', but since
// this is derived from 'Particle', we can match it to the latter.

#include <vector>
#include <string>
#include <math.h>

#include "configuration.h"
#include "FPT_compare.h"
#include "ringstructure.h"
#include "ringcontainer.h"


ringStructure::ringStructure( const int _nRings, const double _centralRadius, const double _deltaR ) : 
  rings( _nRings ),
  centralRadius( _centralRadius ),
  totalRadius( _centralRadius + ( _nRings - 1 ) * _deltaR ),
  deltaR( _deltaR )
{
  initialize();
}

void ringStructure::initialize(void)
{
  rings[0].relocate( 0, centralRadius );
  for ( unsigned int i = 1; i < rings.size(); i++ )
  {
    rings[i].relocate( rings[i-1].getMaxRadius(), rings[i-1].getMaxRadius() + deltaR );
  }
}

void ringStructure::resize( const int _nRings, const double _centralRadius, const double _deltaR )
{
  centralRadius = _centralRadius;
  totalRadius = _centralRadius + ( _nRings - 1 ) * _deltaR;
  deltaR = _deltaR;
  
  rings.resize( _nRings );
  rings.clear();
  initialize();
}

void ringStructure::setLongitudinalGeometry(const double _y_left, const double _y_right, const double _t)
{
  y_left = _y_left;
  y_right = _y_right;
  timenow = _t;
  delta_z = timenow * ( tanh( y_right ) - tanh( y_left ) );
}


const ringContainer& ringStructure::operator[]( const unsigned int index ) const
{
  if ( index >= rings.size() )
  {
    std::string errMsg = "index out of range in ringStructure";
    throw eRingStructure_error( errMsg );
  }
  return rings[ index ];
}



void ringStructure::addRates( const Particle& _particle )
{
  double xt = _particle.Pos.Perp();
  addRates( xt, _particle );
}


void ringStructure::addRates( const double _xt, const Particle& _particle )
{
  int index = getIndexPure( _xt );
  if ( index < static_cast<int>( rings.size() ) )
  {
    rings[ index ].addRates( _particle );
  }
}


void ringStructure::addParticle( const Particle& _particle )
{
  double xt = _particle.Pos.Perp();
  addParticle( xt, _particle );
}


void ringStructure::addParticle( const double _xt, const Particle& _particle )
{
  if( _particle.FLAVOR > 2*Particle::max_N_light_flavor )
  {
    std::string errMsg = "Heavy quark added to ringStructure which cannot deal with massive particles";
    throw eRingStructure_error( errMsg );
  }
  else
  {
    int index = getIndexPure( _xt );
    if ( index < static_cast<int>( rings.size() ) )
    {
      rings[ index ].addParticle( _particle );
    }
  }
}

void ringStructure::addParticleInFormGeom( const double _xt, const Particle& _particle, const double _time )
{
  if( _particle.FLAVOR > 2*Particle::max_N_light_flavor )
  {
    std::string errMsg = "Heavy quark added to ringStructure which cannot deal with massive particles";
    throw eRingStructure_error( errMsg );
  }
  else
    rings[ getIndex( _xt ) ].addParticleInFormGeom( _particle, _time );
}


void ringStructure::addParticleInFormGeom( const Particle& _particle, const double _time )
{
  if( _particle.FLAVOR > 2*Particle::max_N_light_flavor )
  {
    std::string errMsg = "Heavy quark added to ringStructure which cannot deal with massive particles";
    throw eRingStructure_error( errMsg );
  }
  else
    rings[ getIndex( _particle ) ].addParticleInFormGeom( _particle, _time );
}


void ringStructure::clear()
{
  for ( unsigned int i = 0; i < rings.size(); i++ )
  {
    rings[i].clear();
  }
}


void ringStructure::prepareAverages( const double _dz, const int _Ntest )
{
  for ( unsigned int i = 0; i < rings.size(); i++ )
  {
    rings[i].prepareAverages( _dz, _Ntest );
  }
}

void ringStructure::setValuesAsWorkaround(const int ringIndex, const double Edens, const double Pdens, const double _gamma, const int _numberOfParticles, const double Gdens, const double Qdens)
{
  rings[ringIndex].setValuesAsWorkaround(Edens,Pdens,_gamma,_numberOfParticles,Gdens,Qdens);
}


double ringStructure::find_md2g( const Particle& _particle ) const
{
  int nc = getIndex( _particle );

  if (!rings[nc].getAveragesPrepared())
  {
    // std::cout << "ringStructure::find_md2g: averages not prepared" << std::endl;
    // std::cout << bracketvector << particlestyle2 << _particle << std::endl;
    return ( ns_casc::Ncolor + Particle::N_light_flavor ) * 8 / M_PI * pow( 0.3 , 2.0 );
  }

  double md2g = rings[nc].getAveraged_md2g();
  while ( FPT_COMP_E( md2g , 0.0 ) )
  {
    nc--;
    if ( nc >= 0 )
    {
      // Debye mass from next smallest ring:
      md2g = rings[nc].getAveraged_md2g();
    }
    else
    {
      // was already smallest ring
      // Debye mass in medium with T=300MeV
      md2g = ( ns_casc::Ncolor + Particle::N_light_flavor ) * 8 / M_PI * pow( 0.3 , 2.0 );
    }
  }
  return md2g;
}

double ringStructure::find_md2q( const Particle& _particle ) const
{
  int nc = getIndex( _particle );

  if (!rings[nc].getAveragesPrepared())
  {
    // std::cout << "ringStructure::find_md2q: averages not prepared" << std::endl;
    // std::cout << bracketvector << particlestyle2 << _particle << std::endl;
    return 16 / ( 3 * M_PI ) * pow( 0.3 , 2.0 );
  }

  double md2q = rings[nc].getAveraged_md2q();
  while ( FPT_COMP_E( md2q , 0.0 ) )
  {
    nc--;
    if ( nc >= 0 )
    {
      // Debye mass from next smallest ring:
      md2q = rings[nc].getAveraged_md2q();
    }
    else
    {
      // was already smallest ring
      // Debye mass in medium with T=300MeV
      md2q = 16 / ( 3 * M_PI ) * pow( 0.3 , 2.0 );
    }
  }
  return md2q;
}

void ringStructure::getRateV( Particle * pPart ) const
{
  unsigned int nc = getIndex( *pPart );
          
  pPart->rate22v = rings[nc].rates.getRate( pPart->FLAVOR, c22, GeV );
  pPart->rate23v = rings[nc].rates.getRate( pPart->FLAVOR, c23, GeV );
  pPart->rate32v = rings[nc].rates.getRate( pPart->FLAVOR, c32, GeV );
  
  pPart->rate23 = 0;
  pPart->rate32 = 0;
  pPart->rate22 = 0;
}

void ringStructure::getRate( Particle * pPart ) const
{
  unsigned int nc = getIndex( *pPart );

  pPart->rate22 = rings[nc].rates.getRate( pPart->FLAVOR, c22, GeV );
  pPart->rate23 = rings[nc].rates.getRate( pPart->FLAVOR, c23, GeV );
  pPart->rate32 = rings[nc].rates.getRate( pPart->FLAVOR, c32, GeV );
  pPart->md2g = rings[nc].getAveraged_md2g();
  pPart->md2q = rings[nc].getAveraged_md2q();

}

double ringStructure::getEnergyDensity( const VectorEPxPyPz & Pos ) const
{
  unsigned int nc = getIndexPure( Pos.Perp() );
  return ( nc < rings.size() ) ? rings[nc].getEnergyDensity() : 0.0 ;
}

// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;
