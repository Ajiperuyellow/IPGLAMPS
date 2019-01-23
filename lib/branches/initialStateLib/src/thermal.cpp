//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/thermal.cpp $
//$LastChangedDate: 2016-09-16 16:31:08 +0200 (Fr, 16. Sep 2016) $
//$LastChangedRevision: 2422 $
//$LastChangedBy: senzel $
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------


#include <math.h>
#include "thermal.h"
#include "FPT_compare.h"

void initProjectile( ParticlePrototype & particle, const double E, const FLAVOR_TYPE _F, const double angle_to_y_axis )
{
  particle.m = ParticlePrototype::getMass( _F );
  particle.FLAVOR = _F;
  
  if( particle.m > E )
  {
    std::string errMsg = "Mass of jet is larger than its energy.";
    throw eThermal_error( errMsg );
  }
  
  const double pp = sqrt(pow(E,2.0)-pow(particle.m,2.0));
  particle.Mom = VectorEPxPyPz(E, pp * sin( angle_to_y_axis ), pp * cos( angle_to_y_axis ), 0.0 );
  particle.Pos = VectorTXYZ( 0.0, 0.0, 0.0, 0.0 );
}

double getGluonDensity( const double T, const double fugacity, const bool quantum_statistics, const double vz )
{
  const double gamma = 1.0 / sqrt( 1.0 - pow( vz , 2.0 ) );
  
  if( quantum_statistics )
  {
    // Bose
    const double zeta_3 = 1.20206; // Riemann Zeta function: zeta(3)
    return ( gamma * fugacity * 16 * zeta_3 * pow( T, 3.0 ) / pow( M_PI, 2.0 ) / pow( 0.197, 3.0 ) ); 
  }
  else
  {
    // Boltzmann
    return ( gamma * fugacity * 16 * pow( T, 3.0 ) / pow( M_PI, 2.0 ) / pow( 0.197, 3.0 ) );
  }
}


double getQuarkDensity( const double T, const double fugacity, const bool quantum_statistics, const double vz )
{
  const double gamma = 1.0 / sqrt( 1.0 - pow( vz , 2.0 ) );
  
  if( quantum_statistics )
  {
    // Fermi Dirac
    const double zeta_3 = 1.20206; // Riemann Zeta function: zeta(3)
    return ( gamma * fugacity * 12 * ParticlePrototype::N_light_flavor * 3.0 / 4.0 * zeta_3 * pow( T, 3.0 ) / pow( M_PI, 2.0 ) / pow( 0.197, 3.0 ) ); 
  }
  else
  {
    // Boltzmann
    return ( gamma * fugacity * 12 * ParticlePrototype::N_light_flavor * pow( T, 3.0 ) / pow( M_PI, 2.0 ) / pow( 0.197, 3.0 ) );
  }
}


void initThermal( ParticlePrototype& particle, const double T, const double gluonFugacity, const double quarkFugacity, const bool quantum_statistics )
{
  double gluon_density = getGluonDensity( T, gluonFugacity, quantum_statistics );
  double quark_density = getQuarkDensity( T, quarkFugacity, quantum_statistics );
  double ratio = gluon_density / ( gluon_density + quark_density );

  // determine flavor
  int flav;
  if ( ran2() < ratio )
  {
    particle.FLAVOR = gluon;
  }
  else
  {
    flav = int( ParticlePrototype::N_light_flavor * ran2() ) + 1;
    if ( flav > ParticlePrototype::N_light_flavor )
    {
      flav = ParticlePrototype::N_light_flavor;
    }
    flav = 2 * flav;     // this gives anti-quarks ( anti_up == 2, anti_down == 4 etc. )
    
    if ( ran2() < 0.5 )
    {
      flav -= 1;   // make a quark of the previously sampled anti-quark
    }
    
    particle.FLAVOR = static_cast<FLAVOR_TYPE>( flav );
  }
  
  // sample mometum distribution according to massless particle
  sampleThermalMomenta( particle, T, particle.FLAVOR, quantum_statistics );
}

void initThermal( ParticlePrototype& particle, const FLAVOR_TYPE _flav, const double T, const double gluonFugacity, const double quarkFugacity, const bool quantum_statistics )
{
  particle.FLAVOR = _flav;
  // sample mometum distribution according to massless particle
  sampleThermalMomenta( particle, T, particle.FLAVOR, quantum_statistics );  
}


/** @brief samples energy/momentum-distribution according to Boltzmann, Bose-Einstein or Fermi-Dirac distribution
 */
void sampleThermalMomenta( ParticlePrototype & particle, const double T, const FLAVOR_TYPE flav, const bool quantum_statistics )
{
  double E;
  double costheta, sintheta;
  double phi;

  particle.m = ParticlePrototype::getMass( particle.FLAVOR );
  
  if( quantum_statistics )
  {
    if( flav == gluon )
      E = sample_energy( T, bose );
    else
      E = sample_energy( T, fermi );
  }
  else
  {
    E = sample_energy( T, boltzmann, particle.m );
//     E = boltzmann_rejection( T );
  }

  phi = 2 * M_PI * ran2();
  costheta = 2.0 * ran2() - 1.0;
  sintheta = sqrt( 1.0 - costheta * costheta );

  if( particle.m > E )
  {
    std::string errMsg = "Mass of parton is larger than its energy.";
    throw eThermal_error( errMsg );
  }
  
  const double pp = sqrt(pow(E,2.0)-pow(particle.m,2.0));
  particle.Mom = VectorEPxPyPz(E,
                               pp * sintheta * cos( phi ),
                               pp * sintheta * sin( phi ),
                               pp * costheta);
}


double getEnergydistribution( const double E, const double T, const QUANTUM_STATISTICS statistics, const double m )
{
  double a = 0;
  switch( statistics )
  {
    case boltzmann : 
      a = 0;
      break;
    case bose :
      a = -1.0;
      break;
    case fermi :
      a = 1.0;
      break;
  }
  
  // energy distribution for massive Boltzmann particles f(E) = sqrt(E^2-m^2) * E * exp( -E/T )
  return ( sqrt( E*E - m*m ) * E / ( exp( E / T ) + a ) );
}

double sample_energy( const double T, const QUANTUM_STATISTICS statistics, const double m )
{
  double E;
  double E_new;
  double g = 0, g_new = 0;
  double ratio;
  double r;
  
  const double E_min = m;
  const double E_max = 50.0 * T;
  
  // select initial values of E
  do
  {
    E = E_min + ran2()*(E_max - E_min);
    g = getEnergydistribution( E, T, statistics, m );
  } 
  while( FPT_COMP_E( g, 0.0 ) );
  
  // number of steps in the Markov chain
  const int n_steps = 100;
  
  // do n_steps steps
  // the current location is (t)
  // one steps consists of
  //  - proposing a new point (t') according to a certain proposal distribution
  //  - calculate the matrix element g(t') at this new point
  //  - if g(t') > g(t) accept the new point (t')
  //  - else accept the new point (t') with a probability g(t') / g(t) (last point is modified since function for sampling t is not symmetric
  //
  for(int i=0; i<n_steps; i++)
  {
    do
    {
      E_new = E_min + ran2()*(E_max - E_min);          // propose new u using a uniform distribution over the entire range
    } 
    while( E_new < E_min || E_new > E_max);
      
    g_new = getEnergydistribution( E_new, T, statistics, m );     // calculate the matrix element at the proposed point
    
    ratio = g_new / g;                                    // ratio of g(u',phi') / g(u,phi)
    
    if ( FPT_COMP_GE(ratio,1.0) || ran2() < ratio )       // accept if g(t') > g(t) or with probability "ratio"
    {
      E = E_new;
      g = g_new;
    }
  }
 
  return E;
}


double getMeanThermalEnergy( const double T, const QUANTUM_STATISTICS statistics, const double m)
{
  const double Nsampl = 10000;
  
  double sumE = 0.0;
  
  for( int i = 0; i < Nsampl; i++ )
  {
    sumE += sample_energy( T, statistics, m );
  }
  
  return( sumE / Nsampl );
}


double boltzmann_rejection( const double T )
{
  double E, f_E, p_E, ratio;
  
  do
  {
    double temp = ran2();
    // To be sure that temp != 0.0
    while ( temp == 0.0)
    {
      temp = ran2();
    }
    
    E = -2 * T * log( temp );  //sample area in [0,8] (see above) and obtain corresponding E -> combined in one step
    f_E = 4 / T * exp( -E / ( 2 * T ) );     //comparison function evaluated at E
    p_E = pow( E, 2.0 ) / ( 2 * pow( T, 3.0 ) ) * exp( -E / T );
    ratio = p_E / f_E;
  }
  while ( ran2() > ratio ); //accept if random value [0,f(E)] is in [0,p(E)]
  
  return E;
}


void sampleThermalMomentaWithVelocity(ParticlePrototype & particle, const double T, const double vz)
//samples energy/momentum-distribution according to p(E)dE=E^2/(2T^3)*exp(-E/T)dE
//using the rejection method and a comparison function f(E) = 4/T*exp(-E/(2*T))
//which has an area of 8
{
  double E,f_E,p_E,ratio, gamma;
  double costheta, sintheta;
  double phi;
  double vz_abs;
  
  vz_abs = fabs(vz); // sampling can be only performed for positive vz, take sign into account below
  
  gamma = 1.0/sqrt(1.0-pow(vz,2.0));

  do
  {
    E = -2.0*T/(1.0-vz_abs)/ gamma * log(ran2());  //sample area in [0,8] and obtain corresponding E -> combined in one step
    ratio = (1.0-vz_abs)*gamma*E/T*sinh(gamma*vz_abs*E/T)*exp(-gamma*(1.0+vz_abs)/(2.0*T)*E);// p_E / f_E;
  } while(ran2() > ratio); //accept if random value [0,f(E)] is in [0,p(E)]

  do
  {
    costheta = cos(M_PI*ran2());
    sintheta = sqrt(1.0-costheta*costheta);
    ratio = sintheta * exp(-gamma*vz_abs/T*E*(1.0+costheta))     ;// p_E / f_E;
  } while(ran2() > ratio); //accept if random value [0,f(E)] is in [0,p(E)]

  phi = 2.0*M_PI*ran2();

  particle.m = 0;
  particle.Mom = VectorEPxPyPz(E,
                               E * sintheta * cos( phi ),
                               E * sintheta * sin( phi ),
                               E * costheta * (-vz / vz_abs)); // correct sign for z direction (smapling was done only for positive vz), I do not know why there should be a minus sign!! Probably in the sampling above costheta is differently defined than I would do this... But checked: here must be a minus sign, than the mean velocity of the medium is equal to vz.
}

// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on; 
