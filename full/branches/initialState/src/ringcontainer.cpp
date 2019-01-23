//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/ringcontainer.cpp $
//$LastChangedDate: 2018-04-28 19:57:48 +0200 (Sa, 28. Apr 2018) $
//$LastChangedRevision: 2742 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

// at revision 912, this is identical to full/branches/vector4D/src/ringcontainer.cpp

#include <math.h>

#include "ringcontainer.h"
#include "particle.h"
#include "configuration.h"
#include "FPT_compare.h"


ringContainer::ringContainer() : 
  minRadius( 0 ), maxRadius( 0 ), deltaR( 0 ), 
  rates(), 
  averagesPrepared( false ),
  numberOfParticles( 0 ), 
  numberOfGluons( 0 ), numberOfQuarks( 0 ),
  numberOfActiveParticles( 0 ), 
  md2g( 0 ), md2q( 0 ),
  v( ),
  v_r( 0 ), 
  E( 0 ), inverseE_gluons( 0 ), inverseE_quarks( 0 ), 
  p_z( 0 ), p_r( 0 ), p_t( 0 ), 
  pr2_over_E( 0 ), pz2_over_E( 0 ), pr_pz_over_E( 0 ),
  gamma( 0 ), energyDensity( 0 ), particleDensity( 0 ), 
  gluonDensity( 0 ), quarkDensity( 0 ), volume( 0 ),
  numberOfCollectedRateObjects( 0 )
{
  rates.normalizeRates();
}



ringContainer::ringContainer( const double _minR, const double _maxR ) : 
  minRadius( _minR ), 
  maxRadius( _maxR ), 
  deltaR( _maxR - _minR ),
  rates(),
  averagesPrepared( false ),
  numberOfParticles( 0 ), 
  numberOfGluons( 0 ), numberOfQuarks( 0 ), 
  numberOfActiveParticles( 0 ), 
  md2g( 0 ), md2q( 0 ),
  v( ),
  v_r( 0 ), 
  E( 0 ), inverseE_gluons( 0 ), inverseE_quarks( 0 ),
  p_z( 0 ), p_r( 0 ), p_t( 0 ), 
  pr2_over_E( 0 ), pz2_over_E( 0 ), pr_pz_over_E( 0 ),
  gamma( 0 ), energyDensity( 0 ), particleDensity( 0 ), 
  gluonDensity( 0 ), quarkDensity( 0 ), volume( 0 ),
  numberOfCollectedRateObjects( 0 )
{
  rates.normalizeRates();
}


void ringContainer::clear()
{
  averagesPrepared = false;
  numberOfParticles = 0;
  numberOfGluons = 0;
  numberOfQuarks = 0;
  particleDensity = 0;
  gluonDensity = 0;
  quarkDensity = 0;
  md2g = 0;
  md2q = 0;
  v = VectorXYZ( );
  v_r = 0;
  inverseE_gluons = 0;
  inverseE_quarks = 0;
  E = 0;
  p_r = 0;
  p_z = 0;
  p_t = 0;
  pr2_over_E = 0;
  pz2_over_E = 0;
  pr_pz_over_E = 0;

  rates.clear();
  rates.normalizeRates();
  numberOfCollectedRateObjects = 0;
}



void ringContainer::addParticle( const Particle& _particle )
{
  const double xt = _particle.Pos.Perp();
  const double oneE = 1 / _particle.Mom.E();

  ++numberOfParticles;
  if ( !_particle.free || FPT_COMP_G( _particle.rate22, 0 ) )
  {
    ++numberOfActiveParticles;
  }

  E += _particle.Mom.E();
  if ( _particle.FLAVOR == gluon )
  {
    ++numberOfGluons;
    inverseE_gluons += oneE;
  }
  else
  {
    ++numberOfQuarks;
    inverseE_quarks += oneE;
  }

  VectorEPxPyPz vv = _particle.Mom;
  vv.NormalizeToE();

  v.Z() += vv.Z();

  double pr;
  if ( xt < 1.0e-5 )
  {
    pr = _particle.Mom.Pt();
    v.X() += vv.X();
    v.Y() += vv.Y();
  }
  else
  {
    double h_x =  vv.X() * _particle.Pos.X() / xt;
    double h_y =  vv.Y() * _particle.Pos.Y() / xt;
    pr = (h_x + h_y) *  _particle.Mom.E();
    v.X() += h_x;
    v.Y() += h_y;
  }
  v_r += pr * oneE;

  p_r += pr;
  p_z += _particle.Mom.Pz();
  p_t += _particle.Mom.Perp();
  pr2_over_E += pow( pr, 2 ) * oneE;
  pz2_over_E += pow( _particle.Mom.Pz(), 2 ) * oneE;
  pr_pz_over_E += pr * _particle.Mom.Pz() * oneE;
}



void ringContainer::addParticleInFormGeom( const Particle& _particle, const double _time )
{
  ++numberOfParticles;
  if ( _particle.FLAVOR == gluon )
  {
    ++numberOfGluons;
  }
  else
  {
    ++numberOfQuarks;
  }

  double Eold = _particle.Old.E();
  double oneE = 1 / Eold;

  double cc = ( _time - _particle.Pos.T() ) * oneE;

  double xx = _particle.Pos.X() + _particle.Old.Px() * cc;
  double yy = _particle.Pos.Y() + _particle.Old.Py() * cc;
  double xt = sqrt( pow( xx, 2 ) + pow( yy, 2 ) );

  E += Eold;
  if ( _particle.FLAVOR == gluon )
  {
    inverseE_gluons += oneE;
  }
  else
  {
    inverseE_quarks += oneE;
  }
  v.Z() += _particle.Old.Pz() * oneE;

  double pr;
  if ( xt < 1.0e-5 )
  {
    pr = _particle.Old.Pt();
    v.X() += _particle.Old.Px() * oneE;
    v.Y() += _particle.Old.Py() * oneE;
  }
  else
  {
    double h_x = _particle.Old.Px() * xx / xt;
    double h_y = _particle.Old.Py() * yy / xt;
    pr = h_x + h_y;
    v.X() += h_x * oneE;
    v.Y() += h_y * oneE;
  }
  v_r += pr * oneE;

  p_r += pr;
  p_z += _particle.Old.Pz();
  p_t += _particle.Mom.Perp();
  pr2_over_E += pow( pr, 2 ) * oneE;
  pz2_over_E += pow( _particle.Old.Pz(), 2 ) * oneE;
  pr_pz_over_E += pr * _particle.Old.Pz() * oneE;
}


void ringContainer::addRates( const Particle& _particle )
{
  rates.addParticleBasedRates( _particle, GeV );
  ++numberOfCollectedRateObjects;
}


void ringContainer::relocate( const double _minR, const double _maxR )
{
  minRadius = _minR;
  maxRadius = _maxR;
  deltaR = _maxR - _minR;
  clear();
}


double ringContainer::getVolume( const double _dz ) const
{
  return ( M_PI * ( pow( maxRadius, 2 ) - pow( minRadius, 2 ) ) * _dz );   //fm^3
}



double ringContainer::getEnergyDensity() const
{
  if ( !averagesPrepared )
  {
    std::string errMsg = "ringContainer::getEnergyDensity: averaged quantity requested without prior call to prepareAverages()";
    throw eRingContainer_error( errMsg );
  }

  return energyDensity;
}




double ringContainer::getParticleDensity() const
{
  if ( !averagesPrepared )
  {
    std::string errMsg = "ringContainer::getParticleDensity: averaged quantity requested without prior call to prepareAverages()";
    throw eRingContainer_error( errMsg );
  }

  return particleDensity;
}



double ringContainer::getGluonDensity() const
{
  if ( !averagesPrepared )
  {
    std::string errMsg = "ringContainer::getGluonDensity: averaged quantity requested without prior call to prepareAverages()";
    throw eRingContainer_error( errMsg );
  }
  
  return gluonDensity;
}



double ringContainer::getQuarkDensity() const
{
  if ( !averagesPrepared )
  {
    std::string errMsg = "ringContainer::getQuarkDensity: averaged quantity requested without prior call to prepareAverages()";
    throw eRingContainer_error( errMsg );
  }
  
  return quarkDensity;
}



double ringContainer::getGamma() const
{
  if ( !averagesPrepared )
  {
    std::string errMsg = "ringContainer::getGamma: averaged quantity requested without prior call to prepareAverages()";
    throw eRingContainer_error( errMsg );
  }

  return gamma;
}



double ringContainer::getAveraged_md2g() const
{
  if ( !averagesPrepared )
  {
    std::string errMsg = "ringContainer::getAveraged_md2g: averaged quantity requested without prior call to prepareAverages()";
    throw eRingContainer_error( errMsg );
  }

  return md2g;
}



double ringContainer::getAveraged_md2q() const
{
  if ( !averagesPrepared )
  {
    std::string errMsg = "ringContainer::getAveraged_md2q: averaged quantity requested without prior call to prepareAverages()";
    throw eRingContainer_error( errMsg );
  }

  return md2q;
}



void ringContainer::prepareAverages( const double _dz, const int _Ntest )
{
  volume = getVolume( _dz );
  averagesPrepared = true;

  if ( numberOfCollectedRateObjects > 0 )
  {
    rates.prepareParticleBasedAverages();
  }

  if ( numberOfParticles > 0 )
  {
    double invEg = inverseE_gluons / ( ns_casc::gG * _Ntest );
    double invEq = ( Particle::N_light_flavor == 0 )?0:inverseE_quarks / ( 2.0 * ns_casc::gQ  * Particle::N_light_flavor * _Ntest );

    md2g = pow( 0.197, 3 ) * 16 * M_PI / volume * ( ns_casc::Ncolor * invEg + Particle::N_light_flavor * invEq );
    md2q = pow( 0.197, 3 ) * 2 * M_PI / volume * 8.0 / 3.0 * ( invEg + invEq );

    gamma = 1 / sqrt( 1.0 - pow( getAveraged_v_z(), 2 ) - pow( getAveraged_v_r(), 2 ) );
    particleDensity = numberOfParticles / ( _Ntest * volume * gamma );
    gluonDensity = numberOfGluons / ( _Ntest * volume * gamma );
    quarkDensity = numberOfQuarks / ( _Ntest * volume * gamma );

    energyDensity = ( E - ( 2 * getAveraged_v_r() * p_r ) - ( 2 * getAveraged_v_z() * p_z )
                      + ( pow( getAveraged_v_r(), 2 ) * pr2_over_E ) + ( pow( getAveraged_v_z(), 2 ) * pz2_over_E )
                      + ( 2 * getAveraged_v_r() * getAveraged_v_z() * pr_pz_over_E ) ) * pow( gamma, 2 ) / ( _Ntest * volume );                //GeV/fm^3
  }
}


double ringContainer::getEffectiveTemperature() const
{
  if ( !averagesPrepared )
  {
    std::string errMsg = "ringContainer::getEffectiveTemperature: averaged quantity requested without prior call to prepareAverages()";
    throw eRingContainer_error( errMsg );
  }
  
  return energyDensity / (3 * particleDensity);
}


void ringContainer::setValuesAsWorkaround(const double Edens, const double  Pdens,const double _gamma,const int _numberOfParticles,const double Gdens,const double Qdens)
{
  averagesPrepared=true;
  energyDensity=Edens;
  particleDensity=Pdens;
  gamma=_gamma;
  numberOfParticles=_numberOfParticles;
  gluonDensity=Gdens;
  quarkDensity=Qdens;
}


double ringContainer::transformEnergyToComovingFrame(VectorEPxPyPz & P) const
{
  if ( !averagesPrepared )
  {
    std::string errMsg = "ringContainer::transformEnergyToComovingFrame: averaged quantity requested without prior call to prepareAverages()";
    throw eRingContainer_error( errMsg );
  }
  
  return ( numberOfParticles > 0 )? (gamma * ( P.E() - Dot3( v, P ) ) / numberOfParticles ):0;
}


// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;
