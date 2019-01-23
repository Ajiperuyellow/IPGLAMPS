//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/cellcontainer.cpp $
//$LastChangedDate: 2018-04-30 12:04:03 +0200 (Mo, 30. Apr 2018) $
//$LastChangedRevision: 2743 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <math.h>
#include <string>
#include <vector>
#include <list>

#include "cellcontainer.h"
#include "FPT_compare.h"
#include "configuration.h"
#include "lorentz.h"

using namespace ns_casc;
using namespace std;

cellContainer::cellContainer() :
  particleList(),
  averagesPrepared( 0 ),
  nCollectedAll2223( 0 ),
  nCollected22( 0 ),
  nCollected23( 0 ),
  alpha_s_22( 0 ),
  alpha_s_23( 0 ),
  md2g_scaled_22( 0 ),
  md2q_scaled_22( 0 ),
  md2g_scaled_23( 0 ),
  md2q_scaled_23( 0 ),
  sigma_22( 0 ),
  sigma_23( 0 ),
  lambdaScaled( 0 )
{
  //  particleList.clear();
}



void cellContainer::clear()
{
  particleList.clear();
  rates.clear();
  corner.clear();

  nCollectedAll2223 = 0;
  nCollected22 = 0;
  nCollected23 = 0;
  alpha_s_22 = 0;
  alpha_s_23 = 0;
  sigma_22 = 0;
  sigma_23 = 0;
  md2g_scaled_22 = 0;
  md2q_scaled_22 = 0;
  md2g_scaled_23 = 0;
  md2q_scaled_23 = 0;
  lambdaScaled = 0;
  averagesPrepared = false;
}


void cellContainer::resetStoredValues()
{
  nCollectedAll2223 = 0;
  nCollected22 = 0;
  nCollected23 = 0;
  alpha_s_22 = 0;
  alpha_s_23 = 0;
  sigma_22 = 0;
  sigma_23 = 0;
  md2g_scaled_22 = 0;
  md2q_scaled_22 = 0;
  md2g_scaled_23 = 0;
  md2q_scaled_23 = 0;
  lambdaScaled = 0;
  rates.clear();
  averagesPrepared = false;
}



void cellContainer::setCoordinates( const int _index, const double _dx, const int _nx, const double _sizeX, const double _dy, const int _ny, const double _sizeY )
{
  const int nxny = _nx * _ny;
  const int indexEta = _index / nxny;
  const int indexY = _index - ( _index / nxny ) * nxny;
  const int indexX = indexY - ( indexY / _nx ) * _nx;

  double leftY = -( _sizeY / 2.0 ) + _dy * ( indexY / _nx );
  double leftX = -( _sizeX / 2.0 ) + _dx * indexX;

  corner.setCorners( leftX, leftX + _dx, leftY, leftY + _dy, indexEta );
  index = _index;
}

void cellContainer::clearThermodynamicsAndRates()
{
  averagesPreparedSave = false;
  ratesSave.clear();
  ratesSave.normalizeRates();
  numberOfCollectedRateSaveObjects = 0;
}

void cellContainer::defineAsEmpty()
{
  averagesPreparedSave = true;
  ratesSave.clear();
  ratesSave.normalizeRates();
  numberOfCollectedRateSaveObjects = 0;
  energyDensity=0;
  gamma=1.0;
  md2g=1000;
  md2q=1000;
  particleDensity=0.;
  gluonDensity=0.;
  quarkDensity=0.;
}

void cellContainer::prepareThermodynamicsWithRadialBoost(double _volume, int _Ntest)
{
  int numberOfGluons=0;
  int numberOfQuarks=0;
  
  averagesPreparedSave=true;
  numberOfParticles=0;
  
  double E=0;
  VectorXYZ v;
  double v_r=0;
  double inverseE_gluons=0;
  double inverseE_quarks=0;
  double p_z=0;
  double p_r=0;
  double p_t=0;
  double pr2_over_E=0;
  double pz2_over_E=0;
  double pr_pz_over_E=0;
  
  list<int>::const_iterator iIt;
  for ( iIt = particleList.begin(); iIt != particleList.end(); iIt++ )
  {
    const double xt   = particles[( *iIt )].Pos.Perp();
    const double oneE = 1. / particles[( *iIt )].Mom.E();
    //std::cout << xt << std::endl;

    ++numberOfParticles;
    
    E += particles[( *iIt )].Mom.E();
    if ( particles[( *iIt )].FLAVOR == gluon )
    {
      ++numberOfGluons;
      inverseE_gluons += oneE;
    }
    else
    {
      ++numberOfQuarks;
      inverseE_quarks += oneE;
    }

    VectorEPxPyPz vv = particles[( *iIt )].Mom;
    vv.NormalizeToE();

    v.Z() += vv.Z();

    double pr;
    if ( xt < 1.0e-5 )
    {
      pr = particles[( *iIt )].Mom.Pt();
      v.X() += vv.X();
      v.Y() += vv.Y();
    }
    else
    {
      double h_x =  vv.X() * particles[( *iIt )].Pos.X() / xt;
      double h_y =  vv.Y() * particles[( *iIt )].Pos.Y() / xt;
      pr = (h_x + h_y) *  particles[( *iIt )].Mom.E();
      v.X() += h_x;
      v.Y() += h_y;
    }
    v_r += pr * oneE;

    p_r += pr;
    p_z += particles[( *iIt )].Mom.Pz();
    p_t += particles[( *iIt )].Mom.Perp();
    pr2_over_E += pow( pr, 2 ) * oneE;
    pz2_over_E += pow( particles[( *iIt )].Mom.Pz(), 2 ) * oneE;
    pr_pz_over_E += pr * particles[( *iIt )].Mom.Pz() * oneE;  

  }  
  

  volume = _volume;   //fm^3;

  if ( numberOfParticles > 5 )
  {
    double invEg = inverseE_gluons / ( ns_casc::gG * _Ntest );
    double invEq = ( Particle::N_light_flavor == 0 )?0:inverseE_quarks / ( 2.0 * ns_casc::gQ  * Particle::N_light_flavor * _Ntest );

    md2g = pow( 0.197, 3 ) * 16 * M_PI / volume * ( ns_casc::Ncolor * invEg + Particle::N_light_flavor * invEq );
    md2q = pow( 0.197, 3 ) * 2 * M_PI / volume * 8.0 / 3.0 * ( invEg + invEq );

    double Averaged_v_z=v.Z() / numberOfParticles;
    double Averaged_v_r=v_r   / numberOfParticles;
    
    gammaCellRadial = 1 / sqrt( 1.0 - pow( Averaged_v_z, 2 ) - pow( Averaged_v_r, 2 ) );
    
    particleDensityCellRadial = numberOfParticles / ( _Ntest * volume * gammaCellRadial );
    gluonDensityCellRadial = numberOfGluons / ( _Ntest * volume * gammaCellRadial );
    quarkDensityCellRadial = numberOfQuarks / ( _Ntest * volume * gammaCellRadial );

    energyDensityCellRadial = ( E - ( 2 * Averaged_v_r * p_r ) - ( 2 * Averaged_v_z * p_z )
                      + ( pow( Averaged_v_r, 2 ) * pr2_over_E ) + ( pow( Averaged_v_z, 2 ) * pz2_over_E )
                      + ( 2 * Averaged_v_r * Averaged_v_z * pr_pz_over_E ) ) * pow( gammaCellRadial, 2 ) / ( _Ntest * volume );                //GeV/fm^3
  } 
  else
  {
    //keep md2g, gamma, densities
    energyDensityCellRadial=0.;
    gammaCellRadial=1.0;
  }
 
}


void cellContainer::prepareThermodynamicsWithRadialBoostWNeighbors(double _volume, int _Ntest, std::vector< int >& _allParticlesListWNeighbors)
{
  int numberOfGluons=0;
  int numberOfQuarks=0;
  
  averagesPreparedSave=true;
  numberOfParticles=0;
  
  double E=0;
  VectorXYZ v;
  double v_r=0;
  double inverseE_gluons=0;
  double inverseE_quarks=0;
  double p_z=0;
  double p_r=0;
  double p_t=0;
  double pr2_over_E=0;
  double pz2_over_E=0;
  double pr_pz_over_E=0;
  
  vector<int>::const_iterator iIt;
  for ( iIt = _allParticlesListWNeighbors.begin(); iIt != _allParticlesListWNeighbors.end(); iIt++ )
  {
    if(particles[( *iIt )].free)
    {
      continue;      
    }
    const double xt   = particles[( *iIt )].Pos.Perp();
    const double oneE = 1. / particles[( *iIt )].Mom.E();
    //std::cout << xt << std::endl;

    ++numberOfParticles;
    
    E += particles[( *iIt )].Mom.E();
    if ( particles[( *iIt )].FLAVOR == gluon )
    {
      ++numberOfGluons;
      inverseE_gluons += oneE;
    }
    else
    {
      ++numberOfQuarks;
      inverseE_quarks += oneE;
    }

    VectorEPxPyPz vv = particles[( *iIt )].Mom;
    vv.NormalizeToE();

    v.Z() += vv.Z();

    double pr;
    if ( xt < 1.0e-5 )
    {
      pr = particles[( *iIt )].Mom.Pt();
      v.X() += vv.X();
      v.Y() += vv.Y();
    }
    else
    {
      double h_x =  vv.X() * particles[( *iIt )].Pos.X() / xt;
      double h_y =  vv.Y() * particles[( *iIt )].Pos.Y() / xt;
      pr = (h_x + h_y) *  particles[( *iIt )].Mom.E();
      v.X() += h_x;
      v.Y() += h_y;
    }
    v_r += pr * oneE;

    p_r += pr;
    p_z += particles[( *iIt )].Mom.Pz();
    p_t += particles[( *iIt )].Mom.Perp();
    pr2_over_E += pow( pr, 2 ) * oneE;
    pz2_over_E += pow( particles[( *iIt )].Mom.Pz(), 2 ) * oneE;
    pr_pz_over_E += pr * particles[( *iIt )].Mom.Pz() * oneE;  

  }  
  

  volume = _volume;   //fm^3;

  if ( numberOfParticles > 5 )
  {
    double invEg = inverseE_gluons / ( ns_casc::gG * _Ntest );
    double invEq = ( Particle::N_light_flavor == 0 )?0:inverseE_quarks / ( 2.0 * ns_casc::gQ  * Particle::N_light_flavor * _Ntest );

    md2g = pow( 0.197, 3 ) * 16 * M_PI / volume * ( ns_casc::Ncolor * invEg + Particle::N_light_flavor * invEq );
    md2q = pow( 0.197, 3 ) * 2 * M_PI / volume * 8.0 / 3.0 * ( invEg + invEq );

    double Averaged_v_z=v.Z() / numberOfParticles;
    double Averaged_v_r=v_r   / numberOfParticles;
    
    boostvector.SetTXYZ(sqrt(pow(Averaged_v_r,2.0)+pow(Averaged_v_z,2.0)),Averaged_v_r,0,Averaged_v_z);//a HACK since never used, only for analysis.
    
    gammaCellRadial = 1 / sqrt( 1.0 - pow( Averaged_v_z, 2 ) - pow( Averaged_v_r, 2 ) );
    
    particleDensityCellRadial = numberOfParticles / ( _Ntest * volume * gammaCellRadial );
    gluonDensityCellRadial = numberOfGluons / ( _Ntest * volume * gammaCellRadial );
    quarkDensityCellRadial = numberOfQuarks / ( _Ntest * volume * gammaCellRadial );

    energyDensityCellRadial = ( E - ( 2 * Averaged_v_r * p_r ) - ( 2 * Averaged_v_z * p_z )
                      + ( pow( Averaged_v_r, 2 ) * pr2_over_E ) + ( pow( Averaged_v_z, 2 ) * pz2_over_E )
                      + ( 2 * Averaged_v_r * Averaged_v_z * pr_pz_over_E ) ) * pow( gammaCellRadial, 2 ) / ( _Ntest * volume );                //GeV/fm^3
  } 
  else
  {
    //keep md2g, gamma, densities
    energyDensityCellRadial=0.;
    gammaCellRadial=1.0;
  }
 
}

void cellContainer::setValuesToRadialValues()
{
    gamma=gammaCellRadial;
    
    energyDensity = energyDensityCellRadial;
    particleDensity =particleDensityCellRadial;
    gluonDensity = gluonDensityCellRadial;
    quarkDensity = quarkDensityCellRadial; 
}

void cellContainer::prepareThermodynamicsCellBoost(double _volume, int _Ntest)
{
  int numberOfGluons=0;
  int numberOfQuarks=0;
  
  averagesPreparedSave=true;
  
  double E=0;
  double inverseE_gluons=0;
  double inverseE_quarks=0;
  
  double px=0;
  double py=0;
  double pz=0;
  
  double pxy=0;
  double pyz=0;
  double pxz=0;
  
  double pxx=0;
  double pyy=0;
  double pzz=0;
  
  VectorEPxPyPz sumOfMom;
  numberOfParticles=0;
  
  list<int>::const_iterator iIt;
  for ( iIt = particleList.begin(); iIt != particleList.end(); iIt++ )
  {    
    if(particles[( *iIt )].free)
    {
      continue;      
    }
    const double oneE = 1./ particles[( *iIt )].Mom.E();

    ++numberOfParticles;
    
    sumOfMom+=particles[( *iIt )].Mom;
    
    E += particles[( *iIt )].Mom.E();
    px+= particles[( *iIt )].Mom.Px();
    py+= particles[( *iIt )].Mom.Py();
    pz+= particles[( *iIt )].Mom.Pz();

    pxx+= particles[( *iIt )].Mom.Px()*particles[( *iIt )].Mom.Px()/particles[( *iIt )].Mom.E();
    pyy+= particles[( *iIt )].Mom.Py()*particles[( *iIt )].Mom.Py()/particles[( *iIt )].Mom.E();
    pzz+= particles[( *iIt )].Mom.Pz()*particles[( *iIt )].Mom.Pz()/particles[( *iIt )].Mom.E();

    pxy+= particles[( *iIt )].Mom.Px()*particles[( *iIt )].Mom.Py()/particles[( *iIt )].Mom.E();
    pxz+= particles[( *iIt )].Mom.Px()*particles[( *iIt )].Mom.Pz()/particles[( *iIt )].Mom.E();
    pyz+= particles[( *iIt )].Mom.Py()*particles[( *iIt )].Mom.Pz()/particles[( *iIt )].Mom.E();   
    
    if ( particles[( *iIt )].FLAVOR == gluon )
    {
      ++numberOfGluons;
      inverseE_gluons += oneE;
    }
    else
    {
      ++numberOfQuarks;
      inverseE_quarks += oneE;
    }
  }  
  

  volume = _volume;   //fm^3;

  if ( numberOfParticles > 2 )
  {
    //LAB TMUNU [GeV/fm^3]
    double T00,T01,T02,T03,T12,T13,T23,T11,T22,T33;
    T00=E   / ( _Ntest * volume );
    T01=px  /( _Ntest * volume );
    T02=py  /( _Ntest * volume );
    T03=pz  /( _Ntest * volume );
    T12=pxy /( _Ntest * volume );
    T13=pxz /( _Ntest * volume );
    T23=pyz /( _Ntest * volume );
    T11=pxx /( _Ntest * volume );
    T22=pyy /( _Ntest * volume );
    T33=pzz /( _Ntest * volume );
    
    VectorEPxPyPz boostbetaVector =  sumOfMom.NormalizeToE();
    lorentz LorentzBoost;
    LorentzBoost.setBeta( boostbetaVector );
    double gammaCell      =  LorentzBoost.gammaVal();
    gamma=gammaCell;
    boostvector = boostbetaVector;
    double L00,L01,L02,L03;
    L00=gamma;
    L01=-gamma*boostbetaVector.X();
    L02=-gamma*boostbetaVector.Y();
    L03=-gamma*boostbetaVector.Z();

    
    energyDensity = L00*L00*T00       +   L01*L01*T11     +   L02*L02*T22   +   L03*L03*T33
                    + 2*L00*L01*T01   +   2*L00*L02*T02   +   2*L00*L03*T03
                    + 2*L01*L02*T12   +   2*L01*L03*T13   +   2*L02*L03*T23;//GeV/fm^3

    double invEg = inverseE_gluons / ( ns_casc::gG * _Ntest );
    double invEq = ( Particle::N_light_flavor == 0 )?0:inverseE_quarks / ( 2.0 * ns_casc::gQ  * Particle::N_light_flavor * _Ntest );

    md2g = pow( 0.197, 3 ) * 16 * M_PI / volume * ( ns_casc::Ncolor * invEg + Particle::N_light_flavor * invEq );
    md2q = pow( 0.197, 3 ) * 2 * M_PI / volume * 8.0 / 3.0 * ( invEg + invEq );

    
    particleDensity = numberOfParticles / ( _Ntest * volume * gamma );
    gluonDensity = numberOfGluons / ( _Ntest * volume * gamma );
    quarkDensity = numberOfQuarks / ( _Ntest * volume * gamma );
    
  }else
  {
    energyDensity=0;
    gamma=1.0;
    //keep all others
  }
}


void cellContainer::prepareThermodynamicsCellBoostWNeighbors(double _volume, int _Ntest, std::vector< int >& _allParticlesListWNeighbors)
{
  int numberOfGluons=0;
  int numberOfQuarks=0;
  
  averagesPreparedSave=true;
  
  double E=0;
  double inverseE_gluons=0;
  double inverseE_quarks=0;
  
  double px=0;
  double py=0;
  double pz=0;
  
  double pxy=0;
  double pyz=0;
  double pxz=0;
  
  double pxx=0;
  double pyy=0;
  double pzz=0;
  
  VectorEPxPyPz sumOfMom;
  
  numberOfParticles=0;
  
  vector<int>::const_iterator iIt;
  for ( iIt = _allParticlesListWNeighbors.begin(); iIt != _allParticlesListWNeighbors.end(); iIt++ )
  {
    if(particles[( *iIt )].free)
    {
      continue;      
    }
    
    const double oneE = 1./ particles[( *iIt )].Mom.E();

    ++numberOfParticles;
    
    sumOfMom+=particles[( *iIt )].Mom;
    
    E += particles[( *iIt )].Mom.E();
    px+= particles[( *iIt )].Mom.Px();
    py+= particles[( *iIt )].Mom.Py();
    pz+= particles[( *iIt )].Mom.Pz();

    pxx+= particles[( *iIt )].Mom.Px()*particles[( *iIt )].Mom.Px()/particles[( *iIt )].Mom.E();
    pyy+= particles[( *iIt )].Mom.Py()*particles[( *iIt )].Mom.Py()/particles[( *iIt )].Mom.E();
    pzz+= particles[( *iIt )].Mom.Pz()*particles[( *iIt )].Mom.Pz()/particles[( *iIt )].Mom.E();

    pxy+= particles[( *iIt )].Mom.Px()*particles[( *iIt )].Mom.Py()/particles[( *iIt )].Mom.E();
    pxz+= particles[( *iIt )].Mom.Px()*particles[( *iIt )].Mom.Pz()/particles[( *iIt )].Mom.E();
    pyz+= particles[( *iIt )].Mom.Py()*particles[( *iIt )].Mom.Pz()/particles[( *iIt )].Mom.E();   
    
    if ( particles[( *iIt )].FLAVOR == gluon )
    {
      ++numberOfGluons;
      inverseE_gluons += oneE;
    }
    else
    {
      ++numberOfQuarks;
      inverseE_quarks += oneE;
    }
  }  
  

  volume = _volume;   //fm^3;
  volumeFromAVG = volume;
  
  if ( numberOfParticles > 2 )
  {
    //LAB TMUNU [GeV/fm^3]
    double T00,T01,T02,T03,T12,T13,T23,T11,T22,T33;
    T00=E   / ( _Ntest * volume );
    T01=px  / ( _Ntest * volume );
    T02=py  / ( _Ntest * volume );
    T03=pz  / ( _Ntest * volume );
    T12=pxy / ( _Ntest * volume );
    T13=pxz / ( _Ntest * volume );
    T23=pyz / ( _Ntest * volume );
    T11=pxx / ( _Ntest * volume );
    T22=pyy / ( _Ntest * volume );
    T33=pzz / ( _Ntest * volume );
    
    VectorEPxPyPz boostbetaVector =  sumOfMom.NormalizeToE();
    boostvector = boostbetaVector;
    lorentz LorentzBoost;
    LorentzBoost.setBeta( boostbetaVector );
    double gammaCell      =  LorentzBoost.gammaVal();
    //if(gammaCell!=1)cout << gammaCell << endl;
    gamma=gammaCell;
    double L00,L01,L02,L03;
    L00=gamma;
    L01=-gamma*boostbetaVector.X();
    L02=-gamma*boostbetaVector.Y();
    L03=-gamma*boostbetaVector.Z();

    
    energyDensity = L00*L00*T00       +   L01*L01*T11     +   L02*L02*T22   +   L03*L03*T33
                    + 2*L00*L01*T01   +   2*L00*L02*T02   +   2*L00*L03*T03
                    + 2*L01*L02*T12   +   2*L01*L03*T13   +   2*L02*L03*T23;//GeV/fm^3

    double invEg = inverseE_gluons / ( ns_casc::gG * _Ntest );
    double invEq = ( Particle::N_light_flavor == 0 )?0:inverseE_quarks / ( 2.0 * ns_casc::gQ  * Particle::N_light_flavor * _Ntest );

    md2g = pow( 0.197, 3 ) * 16 * M_PI / volume * ( ns_casc::Ncolor * invEg + Particle::N_light_flavor * invEq );
    md2q = pow( 0.197, 3 ) * 2 * M_PI / volume * 8.0 / 3.0 * ( invEg + invEq );

    
    particleDensity = numberOfParticles / ( _Ntest * volume * gamma );
    gluonDensity = numberOfGluons / ( _Ntest * volume * gamma );
    quarkDensity = numberOfQuarks / ( _Ntest * volume * gamma );
    
  }else
  {
    energyDensity=0;
    gamma=1.0;
    //keep all others
  }
}





void cellContainer::prepareAverageRates(double volume, int Ntest)
{
  list<int>::const_iterator iIt;
  for ( iIt = particleList.begin(); iIt != particleList.end(); iIt++ )
  {
    
    ratesSave.addParticleBasedRates( particles[( *iIt )], GeV );
    ++numberOfCollectedRateSaveObjects;
  }

  if ( numberOfCollectedRateSaveObjects > 0 )
  {
    ratesSave.prepareParticleBasedAverages();
  }   
}

void cellContainer::prepareAverageRatesWNeighbors(double volume, int Ntest, std::vector< int >& _allParticlesListWNeighbors)
{
  vector<int>::const_iterator iIt;
  for ( iIt = _allParticlesListWNeighbors.begin(); iIt != _allParticlesListWNeighbors.end(); iIt++ )
  {
    
    ratesSave.addParticleBasedRates( particles[( *iIt )], GeV );
    ++numberOfCollectedRateSaveObjects;
  }

  if ( numberOfCollectedRateSaveObjects > 0 )
  {
    ratesSave.prepareParticleBasedAverages();
  }   
}

double cellContainer::getEnergyDensity() const
{
  if ( !averagesPreparedSave )
  {
    std::string errMsg = "Cell::getEnergyDensity: averaged quantity requested without prior call to prepareAverages()";
    throw eCell_error( errMsg );
  }

  return energyDensity;
}

double cellContainer::getEnergyDensityCellRadial() const
{
  if ( !averagesPreparedSave )
  {
    std::string errMsg = "Cell::getEnergyDensity: averaged quantity requested without prior call to prepareAverages()";
    throw eCell_error( errMsg );
  }

  return energyDensityCellRadial;
}


double cellContainer::getParticleDensity() const
{
  if ( !averagesPreparedSave )
  {
    std::string errMsg = "ringContainer::getParticleDensity: averaged quantity requested without prior call to prepareAverages()";
    throw eCell_error( errMsg );
  }

  return particleDensity;
}



double cellContainer::getGluonDensity() const
{
  if ( !averagesPreparedSave )
  {
    std::string errMsg = "ringContainer::getGluonDensity: averaged quantity requested without prior call to prepareAverages()";
    throw eCell_error( errMsg );
  }
  
  return gluonDensity;
}

int cellContainer::getNumberOfParticles() const
{
  if ( !averagesPreparedSave )
  {
    std::string errMsg = "ringContainer::getGluonDensity: averaged quantity requested without prior call to prepareAverages()";
    throw eCell_error( errMsg );
  }
  
  return numberOfParticles;
}


double cellContainer::getQuarkDensity() const
{
  if ( !averagesPreparedSave )
  {
    std::string errMsg = "ringContainer::getQuarkDensity: averaged quantity requested without prior call to prepareAverages()";
    throw eCell_error( errMsg );
  }
  
  return quarkDensity;
}

double cellContainer::getGamma() const
{
  if ( !averagesPreparedSave )
  {
    std::string errMsg = "ringContainer::getGamma: averaged quantity requested without prior call to prepareAverages()";
    throw eCell_error( errMsg );
  }

  return gamma;
}

void cellContainer::getBoostvector(VectorXYZ & boostbeta) const
{
  if ( !averagesPreparedSave )
  {
    std::string errMsg = "ringContainer::getGamma: averaged quantity requested without prior call to prepareAverages()";
    throw eCell_error( errMsg );
  }
  boostbeta = boostvector;
}

void cellContainer::prepareAverages()
{
  if ( !averagesPrepared )
  {
    averagesPrepared = true;
    if ( nCollectedAll2223 > 0 )
    {
      sigma_22 /= static_cast<double>( nCollectedAll2223 );
      sigma_23 /= static_cast<double>( nCollectedAll2223 );
    }
    else
    {
      sigma_22 = 0;
      sigma_23 = 0;
    }

    if ( nCollected22 > 0 )
    {
      md2g_scaled_22 /= static_cast<double>( nCollected22 );
      md2q_scaled_22 /= static_cast<double>( nCollected22 );
      alpha_s_22 /= static_cast<double>( nCollected22 );
    }
    else
    {
      md2g_scaled_22 = 0;
      md2q_scaled_22 = 0;
      alpha_s_22 = 0;
    }

    if ( nCollected23 > 0 )
    {
      md2g_scaled_23 /= static_cast<double>( nCollected23 );
      md2q_scaled_23 /= static_cast<double>( nCollected23 );
      alpha_s_23 /= static_cast<double>( nCollected23 );
      lambdaScaled /= static_cast<double>( nCollected23 );
    }
    else
    {
      md2g_scaled_23 = 0;
      md2q_scaled_23 = 0;
      alpha_s_23 = 0;
      lambdaScaled = 0;
    }
  }
  else
  {
    std::string errMsg = "prepareAverages called for cell that has already been averaged";
    throw eCell_error( errMsg );
  }
}

// The following routine should not be used in full/offlineAnalysis !!!

void cellContainer::writeAveragesToParticle( Particle& _particle ) const
{
  if ( !averagesPrepared )
  {
    std::string errMsg = "writeAveragesToParticle(..) called without prior call to prepareAverages()";
    throw eCell_error( errMsg );
  }
  
  _particle.cs22 = sigma_22;               //1/GeV^2
  _particle.cs23 = sigma_23;               //1/GeV^2
  _particle.md2g_scaled_22 = md2g_scaled_22;
  _particle.md2q_scaled_22 = md2q_scaled_22;
  _particle.md2g_scaled_23 = md2g_scaled_23;
  _particle.md2q_scaled_23 = md2q_scaled_23;
  _particle.as22 = alpha_s_22;
  _particle.as23 = alpha_s_23;
  _particle.lambda_scaled = lambdaScaled;
  
}





cornerCoordinates::cornerCoordinates() :
    x_min( 0 ),
    x_max( 0 ),
    y_min( 0 ),
    y_max( 0 ),
    etaIndex( -1 )
{

}


cornerCoordinates::cornerCoordinates( const double _x_min, const double _x_max, const double _y_min, const double _y_max, const int _etaIndex ) :
    x_min( _x_min ),
    x_max( _x_max ),
    y_min( _y_min ),
    y_max( _y_max ),
    etaIndex( _etaIndex )
{
}



void cornerCoordinates::setCorners( const double _x_min, const double _x_max, const double _y_min, const double _y_max,  const int _etaIndex )
{
  x_min = _x_min;
  x_max = _x_max;
  y_min = _y_min;
  y_max = _y_max;
  etaIndex = _etaIndex;
}


double cornerCoordinates::getVolume( const coordinateEtaBins& _etaBins, const double _time ) const
{
  double deltaZ = _time * ( tanh( _etaBins[etaIndex].right ) - tanh( _etaBins[etaIndex].left ) );
  return (( x_max - x_min ) * ( y_max - y_min ) * deltaZ );
}

// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;
