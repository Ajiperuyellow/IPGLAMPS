//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/tests/99_ratesmanager.cpp $
//$LastChangedDate: 2017-04-06 15:21:43 +0200 (Do, 06. Apr 2017) $
//$LastChangedRevision: 2558 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/**
 * @brief This is a zero content example for the testing stuff.
 *
 * We do not test anything here, but just provide the basic skeleton
 * of such a file.
 *
 * Running this file as an executable with a return value ==0 is
 * interpreted as *success*. You may indicate the internal number of
 * the failing test by any value != 0.
 **/

#include <iostream>

#include "FPT_compare.h"
#include "thermal.h"
#include "ratesmanager.h"
#include "globalsettings.h"
#include "interpolation22.h"
#include "interpolation23.h"
#include "binary_cross_sections.h"
#include "scattering22.h"
#include "scattering23.h"
#include "scattering32.h"

// Some constants:
const double T = 0.4;
const double dt = 0.005;
const double particles_per_cell = 100;
const double dv = particles_per_cell / ( getGluonDensity( T ) + getQuarkDensity( T ) );

const int Nf_light = 3;
const int Nf_heavy = 0;
const bool running_coupling = false;
const double maxRunningCoupling = 1.0;
const double lambda_gluon = 0.289622 / 0.197; //GeV^-1
const double lambda_quark = 0.617791 / 0.197; //GeV^-1

const int runs = 10000;

const std::string sep = "\t";

using namespace std;
using namespace ns_casc;

// Returns value of Binomial Coefficient C(n, k)
int binomialCoeff(int n, int k)
{
    int res = 1;
 
    // Since C(n, k) = C(n, n-k)
    if ( k > n - k )
        k = n - k;
 
    // Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
    for (int i = 0; i < k; ++i)
    {
        res *= (n - i);
        res /= (i + 1);
    }
 
    return res;
}

void countFlavors( const ParticlePrototype _particle, int &Ng, int &Nq, int &Nqbar, int &Nc, int &Ncbar, int &Nb, int &Nbbar )
{
  switch( ParticlePrototype::mapToGenericFlavorType( _particle.FLAVOR ) )
  {
    case gluon:
      Ng++;
      break;
    case light_quark:
      Nq++;
      break;
    case anti_light_quark:
      Nqbar++;
      break;
    case charm:
      Nc++;
      break;
    case anti_charm:
      Ncbar++;
      break;
    case bottom:
      Nb++;
      break;
    case anti_bottom:
      Nbbar++;
      break;
    default:
      string errMsg = "Unknown flavor type.";
      throw eRatesManager_error( errMsg );
  }
}

ratesManager calcRatesLib( interpolation23& theI23_massless, interpolation23& theI23_charm_m1, interpolation23& theI23_charm_m2, interpolation23& theI23_bottom_m1, interpolation23& theI23_bottom_m2, interpolation22& theI22 )
{
  int Ng = 0, Nq = 0, Nqbar = 0, Nc = 0, Ncbar = 0, Nb = 0, Nbbar = 0;
  
  double md2g_wo_as = ( Ncolor + ParticlePrototype::N_light_flavor ) * 8 / M_PI * pow( T, 2.0 ); //GeV^2
  double md2q_wo_as = 16 / ( 3 * M_PI ) * pow( T, 2.0 ); //GeV^2
  
  int initialStateIndex;
  int nGluons = -1;

  ParticlePrototype part1, part2, part3;
  FLAVOR_TYPE F1, F2, F3;
  double M1, M2, M3;

  double lambda_scaled, s, Vrel, beta, theta, velocity;
  double cs22, cs23, I32;

  scattering22 scatt22_object( &theI22 );
  scattering23 scatt23_object( &theI23_massless, &theI23_charm_m1, &theI23_charm_m2, &theI23_bottom_m1, &theI23_bottom_m2 );
  scattering32 scatt32_object;

  ratesManager vSigma; //in the case of 3->2 vSigma32 is  I32 / (E1*E2*E3)
  
  vector<ParticlePrototype> particles;
  for( int i = 0; i < particles_per_cell; i++ )
  {
    ParticlePrototype sampledParticle;
    initThermal( sampledParticle, T );
    particles.push_back( sampledParticle );
    countFlavors( sampledParticle, Ng, Nq, Nqbar, Nc, Ncbar, Nb, Nbbar );
  }
  
  cout << setw(10) << left << "precise" << sep << setw(10) << left << "calculated" << endl;
  cout << setw(10) << left << getGluonDensity( T ) << sep << setw(10) << left << Ng / dv << endl;
  cout << setw(10) << left << getQuarkDensity( T ) / 2.0 << sep << setw(10) << left << Nq / dv << endl;
  cout << setw(10) << left << getQuarkDensity( T ) / 2.0 << sep << setw(10) << left << Nqbar / dv << endl << endl << endl;

  cout << setw(10) << left << "# of particles" << sep << setw(10) << left << "Ng" << sep << setw(10) << left << "Nq" << sep << setw(10) << left << "Nqbar" << sep << setw(10) << left << "Nc" << sep << setw(10) << left << "Ncbar" << sep << setw(10) << left << "Nb" << sep << setw(10) << left << "Nbbar" << endl;
  cout << setw(10) << left << "" << sep << setw(10) << left << Ng << sep << setw(10) << left << Nq << sep << Nqbar << sep << setw(10) << left << Nc << sep << setw(10) << left << Ncbar << sep << setw(10) << left << Nb << sep << Nbbar << endl;

  int N22 = 0;
  int N32 = 0;
  
  vector<int> Nproc22, Nproc23, Nproc32;
  Nproc22.resize( 9, 0 );
  Nproc23.resize( 9, 0 );
  Nproc32.resize( 9, 0 );
  
  vSigma.clear();
  for ( vector<ParticlePrototype>::const_iterator iIt = particles.begin(); iIt != particles.end(); iIt++ )
  {
    part1 = *iIt;
    vector<ParticlePrototype>::const_iterator start_j = iIt;
    start_j++;
    
    for ( vector<ParticlePrototype>::const_iterator jIt = start_j; jIt != particles.end(); jIt++ )
    {
      part2 = *jIt;

      F1 = part1.FLAVOR;
      M1 = part1.m;

      F2 = part2.FLAVOR;
      M2 = part2.m;

      int type = interactionType::getInclusiveProcessType( F1, F2, c22 );
      Nproc22[ interactionType::getIndexFromProcessType(type, c22 ) ]++;
      
      type = interactionType::getInclusiveProcessType( F1, F2, c23 );
      Nproc23[ interactionType::getIndexFromProcessType(type, c23 ) ]++;

      // 2<->2  ------------------------
      s = ( part1.Mom + part2.Mom ).M2();
      Vrel = VelRel( part1.Mom, part2.Mom, M1, M2 ); // general relative velocity

      if ( s > 0.1 )
      {
        scatt22_object.setParameter( part1.Mom, part2.Mom, F1, F2, M1, M2, s, md2g_wo_as , md2q_wo_as,
                                      1.0, 1.0, 1.0, false,
                                      10.0, false,
                                      1.0 ); // md2g, md2q are debye masses without the factor alpha_s which is multiplied in scattering22.cpp
        cs22 = scatt22_object.getXSection22( initialStateIndex );
      }
      else
      {
        cs22 = 0.0;
      }
      
      double probab22 = pow( 0.197, 2.0 ) * cs22 * Vrel * ( dt / dv );
      vSigma.add( c22, F1, F2, probab22 );
      N22++;
      //-------------------------------

      // 2->3 ------------------------
      s = ( part1.Mom + part2.Mom ).M2();
      Vrel = VelRel( part1.Mom, part2.Mom, M1, M2 ); // general relative velocity

      double lambdaA;
      if( part1.FLAVOR == gluon )
        lambdaA = lambda_gluon;
      else
        lambdaA = lambda_quark;

      double lambdaB;
      if( part2.FLAVOR == gluon )
        lambdaB = lambda_gluon;
      else
        lambdaB = lambda_quark;
      
      double lambda = ( lambdaA + lambdaB ) / 2.0;
      
      lambda_scaled = ( lambda * sqrt( s ) );
      if ( log( lambda_scaled ) >= 15.0 )
      {
        cout << log( lambda_scaled ) << "   " << lambda << "  " << sqrt( s ) <<  endl;
      }

      initialStateIndex = -1;
      scatt23_object.setParameter( VectorXYZ(), part1.Mom, part2.Mom, F1, F2, M1, M2, sqrt( s ), 
                        md2g_wo_as / s, lambda_scaled, 
                        1.0, 1.0,
                        1.0, 1.0,
                        false,
                        "bampsOrgExtended", "GBimproved", false, 0.3,
                        nGluons, false );
      cs23 = scatt23_object.getXSection23( initialStateIndex ); //1/GeV^2 
      
      double probab23 = pow( 0.197, 2.0 ) * cs23 * Vrel * ( dt / dv );
      vSigma.add( c23, F1, F2, probab23 );
    }
  }
  
  //-------------------------------
  for ( vector<ParticlePrototype>::const_iterator iIt = particles.begin(); iIt != particles.end(); iIt++ )
  {
    part1 = *iIt;
    vector<ParticlePrototype>::const_iterator start_j = iIt;
    start_j++;
    
    for ( vector<ParticlePrototype>::const_iterator jIt = start_j; jIt != particles.end(); jIt++ )
    {
      part2 = *jIt;
      vector<ParticlePrototype>::const_iterator start_k = jIt;
      start_k++;
    
      for ( vector<ParticlePrototype>::const_iterator kIt = start_k; kIt != particles.end(); kIt++ )
      {
        part3 = *kIt;

        N32++;
        
        F1 = part1.FLAVOR;
        M1 = part1.m;

        F2 = part2.FLAVOR;
        M2 = part2.m;

        F3 = part3.FLAVOR;
        M3 = part3.m;

        double lambdaA;
        if( part1.FLAVOR == gluon )
          lambdaA = lambda_gluon;
        else
          lambdaA = lambda_quark;

        double lambdaB;
        if( part2.FLAVOR == gluon )
          lambdaB = lambda_gluon;
        else
          lambdaB = lambda_quark;
      
        double lambdaC;
        if( part3.FLAVOR == gluon )
          lambdaC = lambda_gluon;
        else
          lambdaC = lambda_quark;

        double lambda = ( lambdaA + lambdaB + lambdaC ) / 3.0;

        if ( !( F1 == gluon || F2 == gluon || F3 == gluon ) )
        {
          continue;
        }

        int type = interactionType::getInclusiveProcessType( F1, F2, F3, c32 );
        Nproc32[ interactionType::getIndexFromProcessType(type, c32 ) ]++;

        // 3->2 ------------------------
        s = ( part1.Mom + part2.Mom + part3.Mom ).M2();
        lambda_scaled = lambda * sqrt( s );

        scatt32_object.setParameter( VectorXYZ(), part1.Mom, part2.Mom, part3.Mom, F1, F2, F3, sqrt( s ), md2g_wo_as * coupling::get_constant_coupling() / s, lambda_scaled, coupling::get_constant_coupling(), false, false, 0.3, -1 );
        I32 = scatt32_object.getIntegral32_withPrefactors( initialStateIndex );

        double tempI32 = I32 / pow( 0.197, 5.0 );
        
        double probab32 = I32 * dt / pow( dv, 2.0 );
        vSigma.add( c32, F1, F2, F3, probab32 );
        //-------------------------------
      }
    }
  }

  vector<int> Nfactors22, Nfactors23, Nfactors32;
  Nfactors22.resize( 9, 0 );
  Nfactors23.resize( 9, 0 );
  Nfactors32.resize( 9, 0 );

  Nfactors22[0] = Nfactors23[0] = 1.0 / 2.0 * pow( Ng, 2.0 );
  Nfactors22[1] = Nfactors23[1] = Ng * Nq;
  Nfactors22[2] = Nfactors23[2] = Ng * Nqbar;
  Nfactors22[3] = Nfactors23[3] = 1.0 / 3.0 * Nq * Nqbar;
  Nfactors22[4] = Nfactors23[4] = 1.0 / 6.0 * pow( Nq, 2.0 );
  Nfactors22[5] = Nfactors23[5] = 1.0 / 6.0 * pow( Nqbar, 2.0 );
  Nfactors22[6] = Nfactors23[6] = 1.0 / 3.0 * pow( Nq, 2.0 );
  Nfactors22[7] = Nfactors23[7] = 2.0 / 3.0 * Nq * Nqbar;
  Nfactors22[8] = Nfactors23[8] = 1.0 / 3.0 * pow( Nqbar, 2.0 );
  
  Nfactors32[0] = 1.0 / 6.0 * pow( Ng, 3.0 );
  Nfactors32[1] = 1.0 / 2.0 * pow( Ng, 2.0 ) * Nq;
  Nfactors32[2] = 1.0 / 2.0 * pow( Ng, 2.0 ) * Nqbar;
  Nfactors32[3] = 1.0 / 3.0 * Ng * Nq * Nqbar;
  Nfactors32[4] = 1.0 / 6.0 * pow( Nq, 2.0 ) * Ng;
  Nfactors32[5] = 1.0 / 6.0 * pow( Nqbar, 2.0 ) * Ng;
  Nfactors32[6] = 1.0 / 3.0 * pow( Nq, 2.0 ) * Ng;
  Nfactors32[7] = 2.0 / 3.0 * Ng * Nq * Nqbar;
  Nfactors32[8] = 1.0 / 3.0 * pow( Nqbar, 2.0 ) * Ng;

  for( int i = 0; i < 9; i++ )
  {
    cout << i << sep << Nproc22[i] << sep << Nfactors22[i] << sep << Nproc23[i] << sep << Nfactors23[i] << sep << Nproc32[i] << sep << Nfactors32[i] << endl;
  }
  
  vSigma.normalizeRates( Ng, Nq, Nqbar, Nc, Ncbar, Nb, Nbbar, dt, dv );
  
  return vSigma;
}

ratesManager calcRates_mfp( interpolation23& theI23_massless, interpolation23& theI23_charm_m1, interpolation23& theI23_charm_m2, interpolation23& theI23_bottom_m1, interpolation23& theI23_bottom_m2, interpolation22& theI22 )
{
  int Ng = 0, Nq = 0, Nqbar = 0, Nc = 0, Ncbar = 0, Nb = 0, Nbbar = 0;

  double md2g_wo_as = ( Ncolor + ParticlePrototype::N_light_flavor ) * 8 / M_PI * pow( T, 2.0 ); //GeV^2
  double md2q_wo_as = 16 / ( 3 * M_PI ) * pow( T, 2.0 ); //GeV^2
  
  int initialStateIndex;
  int nGluons = -1;

  ParticlePrototype part1, part2, part3;
  FLAVOR_TYPE F1, F2, F3;
  double M1, M2, M3;

  double lambda_scaled, s, Vrel;
  double cs22, cs23, I32;

  scattering22 scatt22_object( &theI22 );
  scattering23 scatt23_object( &theI23_massless, &theI23_charm_m1, &theI23_charm_m2, &theI23_bottom_m1, &theI23_bottom_m2 );
  scattering32 scatt32_object;

  ratesManager vSigma; //in the case of 3->2 vSigma32 is  I32 / (E1*E2*E3)
  
  vSigma.clear();
  for( int i = 0; i < runs; i++ )
  {
    initThermal( part1, T );
    countFlavors( part1, Ng, Nq, Nqbar, Nc, Ncbar, Nb, Nbbar );
    initThermal( part2, T );
    countFlavors( part1, Ng, Nq, Nqbar, Nc, Ncbar, Nb, Nbbar );
    if( part1.FLAVOR != gluon && part2.FLAVOR != gluon )
    {
      initThermal( part3, gluon, T );
    }
    else
    {
      initThermal( part3, T );
    }
    countFlavors( part3, Ng, Nq, Nqbar, Nc, Ncbar, Nb, Nbbar );

    F1 = part1.FLAVOR;
    M1 = part1.m;

    F2 = part2.FLAVOR;
    M2 = part2.m;

    F3 = part3.FLAVOR;
    M3 = part3.m;
    
    double lambdaA;
    if( part1.FLAVOR == gluon )
      lambdaA = lambda_gluon;
    else
      lambdaA = lambda_quark;

    double lambdaB;
    if( part2.FLAVOR == gluon )
      lambdaB = lambda_gluon;
    else
      lambdaB = lambda_quark;
  
    double lambdaC;
    if( part3.FLAVOR == gluon )
      lambdaC = lambda_gluon;
    else
      lambdaC = lambda_quark;

    double lambda = ( lambdaA + lambdaB + lambdaC ) / 3.0;

    // 2<->2  ------------------------
    s = ( part1.Mom + part2.Mom ).M2();
    Vrel = VelRel( part1.Mom, part2.Mom, M1, M2 ); // general relative velocity

    if ( s > 0.1 )
    {
      scatt22_object.setParameter( part1.Mom, part2.Mom, F1, F2, M1, M2, s, md2g_wo_as , md2q_wo_as,
                                    1.0, 1.0, 1.0, false,
                                    10.0, false,
                                    1.0 ); // md2g, md2q are debye masses without the factor alpha_s which is multiplied in scattering22.cpp
      cs22 = scatt22_object.getXSection22( initialStateIndex );
    }
    else
    {
      cs22 = 0.0;
    }
    
    vSigma.add( c22, F1, F2, Vrel * cs22 );
    //-------------------------------


    // 2->3 ------------------------
    s = ( part1.Mom + part2.Mom ).M2();
    Vrel = VelRel( part1.Mom, part2.Mom, M1, M2 ); // general relative velocity

    lambda_scaled = ( lambda * sqrt( s ) );
    if ( log( lambda_scaled ) >= 15.0 )
    {
      cout << log( lambda_scaled ) << "   " << lambda << "  " << sqrt( s ) <<  endl;
    }

    initialStateIndex = -1;
    scatt23_object.setParameter( VectorXYZ(), part1.Mom, part2.Mom, F1, F2, M1, M2, sqrt( s ), 
                      md2g_wo_as / s, lambda_scaled, 
                      1.0, 1.0,
                      1.0, 1.0,
                      false,
                      "bampsOrgExtended", "GBimproved", false, 0.3,
                      nGluons, false );
    cs23 = scatt23_object.getXSection23( initialStateIndex ); //1/GeV^2 
    
    vSigma.add( c23, F1, F2, Vrel * cs23 );

    //-------------------------------

    // 3->2 ------------------------
    s = ( part1.Mom + part2.Mom + part3.Mom ).M2();
    lambda_scaled = lambda * sqrt( s );

    scatt32_object.setParameter( VectorXYZ(), part1.Mom, part2.Mom, part3.Mom, F1, F2, F3, sqrt( s ), md2g_wo_as * coupling::get_constant_coupling() / s, lambda_scaled, coupling::get_constant_coupling(), false, false, 0.3, -1 );
    I32 = scatt32_object.getIntegral32_withPrefactors( initialStateIndex );

    double tempI32 = I32 / pow( 0.197, 5.0 );
    
    vSigma.add( c32, F1, F2, F3, tempI32 );
    //-------------------------------
  }  
  
  vSigma.normalizeRates( getGluonDensity( T ), getQuarkDensity( T ) / 2, getQuarkDensity( T ) / 2 );
  
  return vSigma;
}

/**
 * @brief The main function of the testing programs
 **/
int main(int argc, char **argv)
{
  try
  {
    coupling::set_isRunning( running_coupling );

    ParticlePrototype::set_N_light_flavor( Nf_light );
    ParticlePrototype::set_N_heavy_flavor( Nf_heavy );
  
    ParticlePrototype::setCharmMass( 1.3 );
    ParticlePrototype::setBottomMass( 4.6 );
   
    interpolation22 theI22;
    theI22.configure( running_coupling, ParticlePrototype::N_light_flavor, ParticlePrototype::N_heavy_flavor, ParticlePrototype::Mcharm, ParticlePrototype::Mbottom, 1.0, coupling::get_constant_coupling() );
  
    /** @brief  interpolation23 object that provides access to tabulated values for the cross section of all 2->3 processes with running coupling */
    // do not load data files right at construction, but after configure() has been called below
    interpolation23 theI23_massless( false );
    interpolation23 theI23_charm_m1( false );
    interpolation23 theI23_charm_m2( false );
    interpolation23 theI23_bottom_m1( false );
    interpolation23 theI23_bottom_m2( false );
    // load 2->3 cross section interpolation data
    theI23_massless.configure( false, 1, 0.0, 1.0, "bamps_org_extended", "GBimproved", false, 0.3, log_interpol );
    if( ParticlePrototype::N_heavy_flavor > 0 )
    {
      theI23_charm_m1.configure( false, 1, ParticlePrototype::Mcharm, 1.0, "bampsOrgExtended", "GBimproved", false, 0.3, log_interpol );
      theI23_charm_m2.configure( false, 2, ParticlePrototype::Mcharm, 1.0, "bampsOrgExtended", "GBimproved", false, 0.3, log_interpol );
    }
    if( ParticlePrototype::N_heavy_flavor > 1 )
    {
      theI23_bottom_m1.configure( false, 1, ParticlePrototype::Mbottom, 1.0, "bampsOrgExtended", "GBimproved", false, 0.3, log_interpol );
      theI23_bottom_m2.configure( false, 2, ParticlePrototype::Mbottom, 1.0, "bampsOrgExtended", "GBimproved", false, 0.3, log_interpol );
    }

    cout << "dt = " << dt << sep << "particles_per_cell = " << particles_per_cell << sep << "=> dv = " << dv << endl;

    const double Navg = 1.0;
    
    ratesManager rates_lib_average;
    rates_lib_average.normalizeRates();
    for( int i = 0; i < 10.0 * Navg; i++ )
    {
      ratesManager rates_lib = calcRatesLib( theI23_massless, theI23_charm_m1, theI23_charm_m2, theI23_bottom_m1, theI23_bottom_m2, theI22 );
      rates_lib_average += rates_lib;
    }
    rates_lib_average /= ( 10.0 * Navg );
    
    ratesManager rates_mfp_average;
    rates_mfp_average.normalizeRates();
    for( int i = 0; i < Navg; i++ )
    {
      ratesManager rates_mfp = calcRates_mfp( theI23_massless, theI23_charm_m1, theI23_charm_m2, theI23_bottom_m1, theI23_bottom_m2, theI22 );
      rates_mfp_average += rates_mfp;
    }
    rates_mfp_average /= Navg;

    cout << endl;
    cout << "#############################" << endl << endl;
    cout << "ratio of inclusive rates (R_mfp / R_lib):" << endl << endl;

    cout << setw( 12 ) << left << "process" << sep;
    for( int i = 0; i < 9; i++ )
    {
      cout <<  setw( 12 ) << left << i << sep;
    }
    cout << endl;
    
    cout << setw( 12 ) << left << "2->2" << sep;
    for( int j = 0; j < 9; j++ )
    {
      int index = interactionType::getInclusiveProcessTypeFromIndex( j, c22 );
      cout << setw( 12 ) << left << rates_mfp_average.getRateInclusive( index, c22, fm ) / rates_lib_average.getRateInclusive( index, c22, fm ) << sep;
    }
    cout << endl;

    cout << setw( 12 ) << left << "2->3" << sep;
    for( int j = 0; j < 9; j++ )
    {
      int index = interactionType::getInclusiveProcessTypeFromIndex( j, c23 );
      cout << setw( 12 ) << left << rates_mfp_average.getRateInclusive( index, c23, fm ) / rates_lib_average.getRateInclusive( index, c23, fm ) << sep;
    }
    cout << endl;

    cout << setw( 12 ) << left << "3->2" << sep;
    for( int j = 0; j < 9; j++ )
    {
      int index = interactionType::getInclusiveProcessTypeFromIndex( j, c32 );
      cout << setw( 12 ) << left << rates_mfp_average.getRateInclusive( index, c32, fm ) / rates_lib_average.getRateInclusive( index, c32, fm ) << sep;
    }
    cout << endl << endl;

    cout << "#############################" << endl << endl;
    cout << "ratio of exclusive rates (R_mfp / R_lib):" << endl << endl;
    
    cout << setw(7) << left << "flav" << setw(10) << left << "2->2" << sep << setw(10) << left << "2->3" << sep << setw(10) << left << "3->2" << endl;
    for( int i = 0; i < 10; i++ )
    {
      cout << setw(7) << left << i << setw(10) << left << rates_mfp_average.getRate( static_cast<FLAVOR_TYPE>( i ), c22, fm ) / rates_lib_average.getRate( static_cast<FLAVOR_TYPE>( i ), c22, fm ) << sep << setw(10) << left << rates_mfp_average.getRate( static_cast<FLAVOR_TYPE>( i ), c23, fm ) / rates_lib_average.getRate( static_cast<FLAVOR_TYPE>( i ), c23, fm ) << sep << setw(10) << left << rates_mfp_average.getRate( static_cast<FLAVOR_TYPE>( i ), c32, fm ) / rates_lib_average.getRate( static_cast<FLAVOR_TYPE>( i ), c32, fm ) << endl;
    }
  }
  catch (int e)
  {
    std::cout << "An exception occurred. Exception Nr. " << e << std::endl;
    return e;
  }
  return 0;  
}
