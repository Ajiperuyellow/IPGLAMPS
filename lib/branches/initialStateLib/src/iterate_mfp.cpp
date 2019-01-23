//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL: file:///home/fochler/svnrepos/cascade_misc/createMFPtables/trunk/src/iteration.cpp $
//$LastChangedDate: 2011-01-08 00:05:01 +0100 (Sat, 08 Jan 2011) $
//$LastChangedRevision: 264 $
//$LastChangedBy: fochler $
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------


#include <iostream>
#include <fstream>
#include <deque>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>

#include "thermal.h"
#include "iterate_mfp.h"
#include "particleprototype.h"
#include "lorentz.h"
#include "scattering22.h"
#include "scattering23.h"
#include "scattering32.h"
#include "ratesmanager.h"
#include "binary_cross_sections.h"
#include "globalsettings.h"

using namespace std;

/**
* Determines mean free path by iterative procedure
*
* Use bisection of intervals
*
* @param[in] theI23_m1
* @param[in] theI23_m2
* @param[in] proj particle vector of projectile
* @param[in] lambdaStart start value for mean free path in fm
* @param[in] T temperature in GeV
* @param[in] sampling_per_iter number of runs for determining one mean free path value in the iterative procedure
* @param[in] theConfig
* @return mean free path in fm
*/
double iterate_mfp_bisection( const double T, const interpolation23& theI23_massless, const interpolation23& theI23_charm_m1, const interpolation23& theI23_charm_m2, const interpolation23& theI23_bottom_m1, const interpolation23& theI23_bottom_m2, const interpolation22& theI22, ParticlePrototype & proj,
                              const configBase& theConfig, std::vector<ratesManager>& ratesArray, ratesManager& iteratedRates, const double gluonFugacity, const double quarkFugacity, const bool quantumStatistics, int sampling_per_iter )
{
  // Debye mass not yet multiplied with coupling alpha_s because this is done in the scattering routines since the scale depends on the process
  double md2g_wo_as, md2q_wo_as;
  if( quantumStatistics )
  {
    // Bose definition
    md2g_wo_as = ( 1.0 + ParticlePrototype::N_light_flavor / 6.0 ) * 4.0 * M_PI * pow( T, 2.0 ); //GeV^2 TODO: how to consider here fugacities out of eq.? 
    md2q_wo_as = 0; //GeV^2, What is the quantum defintion of quark screening mass?
  }
  else
  {
    // Boltzmann definition
    md2g_wo_as = ( gluonFugacity * ns_casc::Ncolor + quarkFugacity * ParticlePrototype::N_light_flavor ) * 8 / M_PI * pow( T, 2.0 ); //GeV^2
    md2q_wo_as = 16 / ( 3 * M_PI ) * ( ( gluonFugacity * ns_casc::Ncolor + quarkFugacity * ParticlePrototype::N_light_flavor ) / 2.0 ) * pow( T, 2.0 ); //GeV^2
  }      
  
  ParticlePrototype part2, part3;

  int initialStateIndex;

  deque<double> R22Array, R23Array, R32Array;

  FLAVOR_TYPE F1, F2, F3;
  double m1, m2, m3;
  double lambda_scaled, s, vRel, velocity;
  double cs22, cs23, I32;

  VectorEPxPyPz P1, P2, P3;
  
  double lambda,lambdaRes; // fm
  bool lambda_converged, converged, lambda_lambdaRes_converged, lambdaRes_converged; // fm
  double lambdaAvr, lambdaResAvr; // fm
  deque<double> lambdaArray, lambdaResArray; //fm
  
  double epsilon;
  int nIterationsMax;
  
  if( !theConfig.I23onlineIntegrationIsSet() ) // I23 readout is really fast -> more runs
  {
    if( theConfig.getInterpolation23Mode() == mix_interpol )
    {
      nIterationsMax = 40;
      epsilon = 0.05;
      if( sampling_per_iter == 0 )
      {
        sampling_per_iter = 600;
      }
    }
    else
    {
//       nIterationsMax = 60;
//       sampling_per_iter = 1000;
//       epsilon = 0.01;
      nIterationsMax = 100;
      epsilon = 0.005;
      if( sampling_per_iter == 0 )
      {
        sampling_per_iter = 5000;
      }
    }
  }
  else
  {
    nIterationsMax = 30;
    epsilon = 0.05;
    if( sampling_per_iter == 0 )
    {
      sampling_per_iter = 400; // default 400
    }
// Old parameter for de/dx testruns
// if( theConfig.getNRuns() <= 5000 ) // no high precision for dE/dx, just a testrun. So do not compute mfp to a high accuracy
// {
//   sampling_per_iter = 50; // default 400
//   epsilon = 0.1;
//   nIterationsMax = 30;
// }
  }
  
  double lambdaMax = 2.0;  //fm
  double lambdaMin = 0.0;  //fm
  lambda = (lambdaMax - lambdaMin) / 2.0 + lambdaMin;  //fm

  int nIterations = 0;

  scattering22 scatt22_object( &theI22 );
  scattering23 scatt23_object( &theI23_massless, &theI23_charm_m1, &theI23_charm_m2, &theI23_bottom_m1, &theI23_bottom_m2 );
  scattering32 scatt32_object;

  ratesManager vSigma;
  ratesArray.clear();
  iteratedRates.clear();
  iteratedRates.isNormalized = true;
  
  do
  {
    vSigma.clear();
    
    double cs_test_sum = 0.0;

    for ( int i = 1;i <= sampling_per_iter;i++ )
    {
      initThermal( part2, T, gluonFugacity, quarkFugacity );
      if ( proj.FLAVOR != gluon && part2.FLAVOR != gluon )
      {
        initThermal( part3, gluon, T );
      }
      else
      {
        initThermal( part3, T);
      }
      
      if ( ran2() < 0.5 )
      {
        F1 = proj.FLAVOR;
        m1 = proj.m;
        P1 = proj.Mom;
        
        F2 = part2.FLAVOR;
        m2 = part2.m;
        P2 = part2.Mom;
      }
      else
      {
        F2 = proj.FLAVOR;
        m2 = proj.m;
        P2 = proj.Mom;
        
        F1 = part2.FLAVOR;
        m1 = part2.m;
        P1 = part2.Mom;
      }
      F3 = part3.FLAVOR;
      m3 = part3.m;
      P3 = part3.Mom;

      s = (P1+P2).M2();
      vRel = VelRel(P1,P2, m1,m2); // general relative velocity

      // 2<->2  ------------------------
      if ( s > 0.1 )
      {
        scatt22_object.setParameter( P1, P2, F1, F2, m1, m2, s, vRel, md2g_wo_as , md2q_wo_as,
                                       theConfig.getKggQQb(), theConfig.getKgQgQ(), theConfig.getKappa_gQgQ(), 
                                       theConfig.isConstantCrossSecGQ(),
                                       theConfig.getConstantCrossSecValueGQ(), theConfig.isIsotropicCrossSecGQ(), theConfig.getKfactor_light() ); // md2g_wo_as, md2q_wo_as are debye masses without the factor alpha_s which is multiplied in scattering22.cpp
        cs22 = scatt22_object.getXSection22();
      }
      else
      {
        cs22 = 0.0;
      }

      vSigma.add( c22, F1, F2, vRel * cs22 );
      //-------------------------------

      // 2->3 ------------------------
      lambda_scaled = ( lambda / 0.197 * sqrt( s ) ); // lambda in fm, sqrt(s) in GeV, lambda_scaled dimensionless
      if ( log( lambda_scaled ) >= 15.0 )
      {
        cout << log( lambda_scaled ) << "   " << lambda << " fm  " << sqrt( s ) <<  endl;
      }

      initialStateIndex = -1;
      scatt23_object.setParameter( VectorXYZ(), P1, P2, F1, F2, m1, m2, sqrt( s ), 
                        md2g_wo_as / s, lambda_scaled, 
                        theConfig.getK23LightPartons(), theConfig.getK23HeavyQuarks(),
                        theConfig.getKappa23LightPartons(), theConfig.getKappa23HeavyQuarks(),
                        theConfig.I23onlineIntegrationIsSet(),
                        theConfig.get23GluonFormationTimeTyp(), theConfig.getMatrixElement23(), theConfig.isMd2CounterTermInI23(), theConfig.get23FudgeFactorLpm() );    
      cs23 = scatt23_object.getXSection23( initialStateIndex ); //1/GeV^2

      vSigma.add( c23, F1, F2, vRel * cs23 );
      
      cs_test_sum += vRel * cs23;
      //-------------------------------

      // 3->2 ------------------------
      s = ( proj.Mom + part2.Mom + part3.Mom ).M2();
      lambda_scaled = ( lambda / 0.197 * sqrt( s ) ); // lambda in fm, sqrt(s) in GeV, lambda_scaled dimensionless

      initialStateIndex = -1;
      scatt32_object.setParameter( VectorXYZ(), proj.Mom, part2.Mom, part3.Mom, F1, F2, F3, sqrt( s ), md2g_wo_as * coupling::get_constant_coupling() / s, lambda_scaled, coupling::get_constant_coupling(), theConfig.isMd2CounterTermInI23(), theConfig.isMatrixElement23_22qt(), theConfig.get23FudgeFactorLpm(), -1 );
      I32 = scatt32_object.getIntegral32_withPrefactors( initialStateIndex );

      double tempI32 = I32 / pow( 0.197, 5.0 );
      vSigma.add( c32, F1, F2, F3, tempI32 );
      //-------------------------------
    }

    vSigma.normalizeRates( getGluonDensity( T, gluonFugacity, quantumStatistics ), getQuarkDensity( T, quarkFugacity, quantumStatistics ) / 2.0, getQuarkDensity( T, quarkFugacity, quantumStatistics ) / 2.0 );
    
    velocity = sqrt( 1.0 - pow( proj.m ,2.0) / pow( proj.Mom.E() ,2.0) );
    // lambda which ones obtain with the previous employed lambda for the cross section
    lambdaRes = vSigma.getLambda( proj.FLAVOR, velocity, fm ); //fm
    
//     cout << "min = " << lambdaMin << "  max = " << lambdaMax << endl;
//     cout << "lambda = " << lambda << "  lambdaRes = " << lambdaRes << endl;
//     cout << "lambda = " << vSigma.getLambda( proj.FLAVOR, velocity, fm ) << " fm    22: " << vSigma.getLambda( proj.FLAVOR, c22, velocity, fm ) << " fm    23: " << vSigma.getLambda( proj.FLAVOR, c23, velocity, fm ) << endl;
    
    // constrain for bisection: we know that the function lambdaRes( lambda ) is monotonly falling
    if( lambdaRes < lambda )
      lambdaMax = lambda;
    else
      lambdaMin = lambda;
    
    lambdaArray.push_back( lambda );
    lambdaResArray.push_back( lambdaRes );

    ratesArray.push_back( vSigma );
    
    if ( lambdaArray.size() >= 4 )
    {
      lambdaAvr = 0;
      for ( int m = ( int( lambdaArray.size() ) - 1 ); m >= ( int( lambdaArray.size() ) - 4 ); m-- )
        lambdaAvr += lambdaArray[m];
      lambdaAvr = lambdaAvr / 4;
      
      lambdaResAvr = 0;
      for ( int m = ( int( lambdaResArray.size() ) - 1 ); m >= ( int( lambdaResArray.size() ) - 4 ); m-- )
        lambdaResAvr += lambdaResArray[m];
      lambdaResAvr = lambdaResAvr / 4;
      
      lambda_lambdaRes_converged = fabs( lambdaAvr - lambdaResAvr ) / lambdaAvr < epsilon;
      
      lambda_converged = ( ( fabs( lambdaArray[lambdaArray.size()-1] - lambdaAvr ) / lambdaAvr < epsilon ) &&
                    ( fabs( lambdaArray[lambdaArray.size()-2] - lambdaAvr ) / lambdaAvr < epsilon ) &&
                    ( fabs( lambdaArray[lambdaArray.size()-3] - lambdaAvr ) / lambdaAvr < epsilon ) &&
                    ( fabs( lambdaArray[lambdaArray.size()-4] - lambdaAvr ) / lambdaAvr < epsilon ) );
                    
      lambdaRes_converged = ( ( fabs( lambdaResArray[lambdaResArray.size()-1] - lambdaResAvr ) / lambdaResAvr < epsilon ) &&
                    ( fabs( lambdaResArray[lambdaResArray.size()-2] - lambdaResAvr ) / lambdaResAvr < epsilon ) &&
                    ( fabs( lambdaResArray[lambdaResArray.size()-3] - lambdaResAvr ) / lambdaResAvr < epsilon ) &&
                    ( fabs( lambdaResArray[lambdaResArray.size()-4] - lambdaResAvr ) / lambdaResAvr < epsilon ) );
                    
      converged = lambda_lambdaRes_converged && lambda_converged;
                    
      if( lambda_converged && !converged )
      {
        if( lambdaResAvr > lambdaMax )
        {
          lambdaMax = lambdaResAvr;
//           cout << "rechange max = " << lambdaMax << endl;
          lambdaArray.clear();
          lambdaResArray.clear();
        }
        else if( lambdaResAvr < lambdaMin )
        {
          lambdaMin = lambdaResAvr;
//           cout << "rechange min = " << lambdaMin << endl;
          lambdaArray.clear();
          lambdaResArray.clear();
        }
      }
    }
    else
      converged = false;
    
    lambda = (lambdaMax - lambdaMin) / 2.0 + lambdaMin; // fm
    
    nIterations++;
  }
  while ( !converged && nIterations < nIterationsMax );

  lambda = ( lambdaAvr + lambdaResAvr ) / 2.0; // fm
  
  if(converged)
    cout << "converged: lambda = " << lambda << " fm  flavor = " << proj.FLAVOR << " E = " << proj.Mom.E() << " GeV" << endl;
  else
    cout << "not converged: converged: lambda = " << lambda << " fm  flavor = " << proj.FLAVOR << " E = " << proj.Mom.E() << " GeV" << endl;
  
  if ( lambdaArray.size() >= 4 )
  {  
//     cout << lambdaAvr << "  " << lambdaArray[lambdaArray.size()-1] << "  " << lambdaArray[lambdaArray.size()-2] << "  " << lambdaArray[lambdaArray.size()-3] << "  " << lambdaArray[lambdaArray.size()-4] << endl;
//     cout << lambdaResAvr << "  " << lambdaResArray[lambdaResArray.size()-1] << "  " << lambdaResArray[lambdaResArray.size()-2] << "  " << lambdaResArray[lambdaResArray.size()-3] << "  " << lambdaResArray[lambdaResArray.size()-4] << endl;
    
    int stop = static_cast<int>( lambdaArray.size() ) - 4;
    for ( int m = lambdaArray.size() - 1; m >= stop; m-- )
    {
      iteratedRates += ratesArray[m];
    }
    iteratedRates /= 4;
  }

  return lambda; // fm
}



/**
* Determines mean free path by iterative procedure
*
* @param[in] theI23_m1 
* @param[in] theI23_m2
* @param[in] proj particle vector of projectile
* @param[in] lambdaStart start value for mean free path in fm
* @param[in] T temperature in GeV
* @param[in] sampling_per_iter number of runs for determining one mean free path value in the iterative procedure
* @param[in] theConfig
* @return mean free path in fm
*/
double iterate_mfp( const double T, const interpolation23& theI23_massless, const interpolation23& theI23_charm_m1, const interpolation23& theI23_charm_m2, const interpolation23& theI23_bottom_m1, const interpolation23& theI23_bottom_m2, const interpolation22& theI22, ParticlePrototype &proj, const double lambdaStart,
                    const configBase& theConfig, std::vector<ratesManager>& ratesArray, ratesManager& iteratedRates, const double gluonFugacity, const double quarkFugacity, const bool quantumStatistics, const int samplings_per_iter )
{
  int nGluons = -1;
  
  // Debye mass not yet multiplied with coupling alpha_s because this is done in the scattering routines since the scale depends on the process
  double md2g_wo_as, md2q_wo_as;
  if( quantumStatistics )
  {
    // Bose definition
    md2g_wo_as = ( 1.0 + ParticlePrototype::N_light_flavor / 6.0 ) * 4.0 * M_PI * pow( T, 2.0 ); //GeV^2 TODO: how to consider here fugacities out of eq.? 
    md2q_wo_as = 0; //GeV^2, What is the quantum defintion of quark screening mass?
  }
  else
  {
    // Boltzmann definition
    md2g_wo_as = ( gluonFugacity * ns_casc::Ncolor + quarkFugacity * ParticlePrototype::N_light_flavor ) * 8 / M_PI * pow( T, 2.0 ); //GeV^2
    md2q_wo_as = 16 / ( 3 * M_PI ) * ( ( gluonFugacity * ns_casc::Ncolor + quarkFugacity * ParticlePrototype::N_light_flavor ) / 2.0 ) * pow( T, 2.0 ); //GeV^2
  }      
  
  ParticlePrototype part2, part3;

  const double epsilon = 0.01;
  double lambda = lambdaStart; // fm
  double lambdaAvr; // fm

  int initialStateIndex;
  const int nIterationsMax = 30;
  bool converged = false;

  deque<double> lambdaArray, R22Array, R23Array, R32Array;

  FLAVOR_TYPE F1, F2, F3;
  double m1, m2, m3;
  double lambda_scaled, s, vRel, velocity;
  double cs22, cs23, I32;

  VectorEPxPyPz P1, P2, P3;

  scattering22 scatt22_object( &theI22 );
  scattering23 scatt23_object( &theI23_massless, &theI23_charm_m1, &theI23_charm_m2, &theI23_bottom_m1, &theI23_bottom_m2 );
  scattering32 scatt32_object;

  ratesManager vSigma;
  ratesArray.clear();
  iteratedRates.clear();
  iteratedRates.isNormalized = true;
  
  do
  {
    vSigma.clear();

    for ( int i = 1;i <= samplings_per_iter;i++ )
    {
      initThermal( part2, T, gluonFugacity, quarkFugacity );
      if ( proj.FLAVOR != gluon && part2.FLAVOR != gluon )
      {
        initThermal( part3, gluon, T );
      }
      else
      {
        initThermal( part3, T, gluonFugacity, quarkFugacity );
      }

      if ( ran2() < 0.5 )
      {
        F1 = proj.FLAVOR;
        m1 = proj.m;
        P1 = proj.Mom;
        
        F2 = part2.FLAVOR;
        m2 = part2.m;
        P2 = part2.Mom;
      }
      else
      {
        F2 = proj.FLAVOR;
        m2 = proj.m;
        P2 = proj.Mom;
        
        F1 = part2.FLAVOR;
        m1 = part2.m;
        P1 = part2.Mom;
      }
      F3 = part3.FLAVOR;
      m3 = part3.m;
      P3 = part3.Mom;

      s = (P1+P2).M2();
      vRel = VelRel(P1,P2, m1,m2); // general relative velocity
      
      // 2<->2  ------------------------
      if( theConfig.doScattering_22() )
      {
	if ( s > 0.1 )
	{
	  scatt22_object.setParameter( P1, P2, F1, F2, m1, m2, s, vRel, md2g_wo_as , md2q_wo_as,
				      theConfig.getKggQQb(), theConfig.getKgQgQ(), theConfig.getKappa_gQgQ(), theConfig.isConstantCrossSecGQ(),
				      theConfig.getConstantCrossSecValueGQ(), theConfig.isIsotropicCrossSecGQ(),
				      theConfig.getKfactor_light() ); // md2g, md2q are debye masses without the factor alpha_s which is multiplied in scattering22.cpp
	  cs22 = scatt22_object.getXSection22();
	}
	else
	{
	  cs22 = 0.0;
	}

	vSigma.add( c22, F1, F2, vRel * cs22 );
      }
      //-------------------------------

      // 2->3 ------------------------
      if( theConfig.doScattering_23() )
      {
	lambda_scaled = ( lambda / 0.197 * sqrt( s ) ); // lambda in fm, sqrt(s) in GeV, lambda_scaled dimensionless
	if ( log( lambda_scaled ) >= 15.0 )
	{
	  cout << log( lambda_scaled ) << "   " << lambda << " fm  " << sqrt( s ) <<  endl;
	}

	initialStateIndex = -1;
	scatt23_object.setParameter( VectorXYZ(), P1, P2, F1, F2, m1, m2, sqrt( s ), 
			  md2g_wo_as / s, lambda_scaled, 
			  theConfig.getK23LightPartons(), theConfig.getK23HeavyQuarks(),
			  theConfig.getKappa23LightPartons(), theConfig.getKappa23HeavyQuarks(),
			  theConfig.I23onlineIntegrationIsSet(),
			  theConfig.get23GluonFormationTimeTyp(), theConfig.getMatrixElement23(), theConfig.isMd2CounterTermInI23(), theConfig.get23FudgeFactorLpm(), 
			  nGluons, theConfig.isMatrixElement23_22qt() );
	cs23 = scatt23_object.getXSection23( initialStateIndex ); //1/GeV^2

	vSigma.add( c23, F1, F2, vRel * cs23 );
      }
      //-------------------------------

      // 3->2 ------------------------
      if( theConfig.doScattering_32() )
      {
	s = ( proj.Mom + part2.Mom + part3.Mom ).M2();
	lambda_scaled = ( lambda / 0.197 * sqrt( s ) ); // lambda in fm, sqrt(s) in GeV, lambda_scaled dimensionless

	initialStateIndex = -1;
	scatt32_object.setParameter( VectorXYZ(), proj.Mom, part2.Mom, part3.Mom, F1, F2, F3, sqrt( s ), md2g_wo_as * coupling::get_constant_coupling() / s, lambda_scaled, coupling::get_constant_coupling(), theConfig.isMd2CounterTermInI23(), theConfig.isMatrixElement23_22qt(), theConfig.get23FudgeFactorLpm(), -1 );
	I32 = scatt32_object.getIntegral32_withPrefactors( initialStateIndex );

	double tempI32 = I32 / pow( 0.197, 5.0 );
	vSigma.add( c32, F1, F2, F3, tempI32 );
      }
      //-------------------------------
    }

    vSigma.normalizeRates( getGluonDensity( T, gluonFugacity, quantumStatistics ), getQuarkDensity( T, quarkFugacity, quantumStatistics ) / 2.0, getQuarkDensity( T, quarkFugacity, quantumStatistics ) / 2.0 );
    
    velocity = sqrt( 1.0 - pow( proj.m ,2.0) / pow( proj.Mom.E() ,2.0) );
    lambda = vSigma.getLambda( proj.FLAVOR, velocity, fm ); //fm

//     cout << lambda << " fm    " <<  lambda / 0.197  << " GeV^-1  " << endl;

    lambdaArray.push_back( lambda );
    ratesArray.push_back( vSigma );

    if ( lambdaArray.size() >= 4 )
    {
      lambdaAvr = 0;
      for ( int m = ( int( lambdaArray.size() ) - 1 ); m >= ( int( lambdaArray.size() ) - 4 ); m-- )
        lambdaAvr += lambdaArray[m];
      lambdaAvr = lambdaAvr / 4;

      if (( fabs( lambdaArray[lambdaArray.size()-1] - lambdaAvr ) / lambdaAvr < epsilon ) &&
          ( fabs( lambdaArray[lambdaArray.size()-2] - lambdaAvr ) / lambdaAvr < epsilon ) &&
          ( fabs( lambdaArray[lambdaArray.size()-3] - lambdaAvr ) / lambdaAvr < epsilon ) &&
          ( fabs( lambdaArray[lambdaArray.size()-4] - lambdaAvr ) / lambdaAvr < epsilon ) )
        converged = true;
      else
        converged = false;
    }
    else
      converged = false;

    if ( lambdaArray.size() >= 6 )
      lambda = lambdaAvr; // fm
  }
  while ( !converged && lambdaArray.size() < nIterationsMax );

  if(converged)
    cout << "converged: lambda = " << lambda << " flavor = " << proj.FLAVOR << " E = " << proj.Mom.E() << " GeV" << endl;
  else
    cout << "not converged: converged: lambda = " << lambda << " flavor = " << proj.FLAVOR << " E = " << proj.Mom.E() << " GeV" << endl;
  
  lambdaAvr = 0;
  int stop = static_cast<int>( lambdaArray.size() ) - 4;
  for ( int m = lambdaArray.size() - 1; m >= stop; m-- )
  {
    lambdaAvr += lambdaArray[m];
    iteratedRates += ratesArray[m];
  }
  lambda = lambdaAvr / 4;
  iteratedRates /= 4;

  return lambda; // fm
}

/**
* Determines thermal mean free path by iterative procedure
*
* @param[in] T temperatur in GeV
* @param[in] flav flavor for which MFP is calculated
* @param[in] theI23_m1 
* @param[in] theI23_m2
* @param[in] proj particle vector of projectile
* @param[in] lambdaStart start value for mean free path in fm
* @param[in] T temperature in GeV
* @param[in] sampling_per_iter number of runs for determining one mean free path value in the iterative procedure
* @param[in] theConfig
* @return mean free path in fm
*/
double iterate_thermal_mfp( const double T, const FLAVOR_TYPE flav, const interpolation23& theI23_massless, const interpolation23& theI23_charm_m1, const interpolation23& theI23_charm_m2, const interpolation23& theI23_bottom_m1, const interpolation23& theI23_bottom_m2, const interpolation22& theI22, const double lambdaStart,
                    const configBase& theConfig, std::vector<ratesManager>& ratesArray, ratesManager& iteratedRates, const double gluonFugacity, const double quarkFugacity, const bool quantumStatistics, const int samplings_per_iter )
{
  int nGluons = -1;
  
  // Debye mass not yet multiplied with coupling alpha_s because this is done in the scattering routines since the scale depends on the process
  double md2g_wo_as, md2q_wo_as;
  if( quantumStatistics )
  {
    // Bose definition
    md2g_wo_as = ( 1.0 + ParticlePrototype::N_light_flavor / 6.0 ) * 4.0 * M_PI * pow( T, 2.0 ); //GeV^2 TODO: how to consider here fugacities out of eq.? 
    md2q_wo_as = 0; //GeV^2, What is the quantum defintion of quark screening mass?
  }
  else
  {
    // Boltzmann definition
    md2g_wo_as = ( gluonFugacity * ns_casc::Ncolor + quarkFugacity * ParticlePrototype::N_light_flavor ) * 8 / M_PI * pow( T, 2.0 ); //GeV^2
    md2q_wo_as = 16 / ( 3 * M_PI ) * ( ( gluonFugacity * ns_casc::Ncolor + quarkFugacity * ParticlePrototype::N_light_flavor ) / 2.0 ) * pow( T, 2.0 ); //GeV^2
  }      
  
  ParticlePrototype part1, part2, part3;
  
  const double epsilon = 0.01;
  double lambda = lambdaStart; // fm
  double lambdaAvr; // fm

  int initialStateIndex;
  const int nIterationsMax = 30;
  bool converged = false;

  deque<double> lambdaArray, R22Array, R23Array, R32Array;

  FLAVOR_TYPE F1, F2, F3;
  double m1, m2, m3;
  double lambda_scaled, s, vRel, velocity;
  double cs22, cs23, I32;

  VectorEPxPyPz P1, P2, P3;

  scattering22 scatt22_object( &theI22 );
  scattering23 scatt23_object( &theI23_massless, &theI23_charm_m1, &theI23_charm_m2, &theI23_bottom_m1, &theI23_bottom_m2 );
  scattering32 scatt32_object;

  ratesManager vSigma; //in the case of 3->2 vSigma32 is  I32 / (E1*E2*E3)
  ratesArray.clear();
  iteratedRates.clear();
  iteratedRates.isNormalized = true;
  
  do
  {
    vSigma.clear();

    for ( int i = 1;i <= samplings_per_iter;i++ )
    {
      initThermal( part1, flav, T, gluonFugacity, quarkFugacity );
      initThermal( part2, T, gluonFugacity, quarkFugacity );
      if ( flav != gluon && part2.FLAVOR != gluon )
      {
        initThermal( part3, gluon, T );
      }
      else
      {
        initThermal( part3, T, gluonFugacity, quarkFugacity );
      }
        
      if ( ran2() < 0.5 )
      {
        F1 = part1.FLAVOR;
        m1 = part1.m;
        P1 = part1.Mom;

        F2 = part2.FLAVOR;
        m2 = part2.m;
        P2 = part2.Mom;
      }
      else
      {
        F2 = part1.FLAVOR;
        m2 = part1.m;
        P2 = part1.Mom;

        F1 = part2.FLAVOR;
        m1 = part2.m;
        P1 = part2.Mom;
      }
      F3 = part3.FLAVOR;
      m3 = part3.m;
      P3 = part3.Mom;

      s = (P1+P2).M2();
      vRel = VelRel(P1,P2, m1,m2); // general relative velocity

      // 2<->2  ------------------------
      if( theConfig.doScattering_22() )
      {
	if ( s > 0.1 )
	{
	  scatt22_object.setParameter( P1, P2, F1, F2, m1, m2, s, vRel, md2g_wo_as , md2q_wo_as,
				      theConfig.getKggQQb(), theConfig.getKgQgQ(), theConfig.getKappa_gQgQ(), theConfig.isConstantCrossSecGQ(),
				      theConfig.getConstantCrossSecValueGQ(), theConfig.isIsotropicCrossSecGQ(),
				      theConfig.getKfactor_light() ); // md2g, md2q are debye masses without the factor alpha_s which is multiplied in scattering22.cpp
	  cs22 = scatt22_object.getXSection22();
	}
	else
	{
	  cs22 = 0.0;
	}
	
	vSigma.add( c22, F1, F2, vRel * cs22 );
      }
      //-------------------------------

      // 2->3 ------------------------
      if( theConfig.doScattering_23() )
      {
	lambda_scaled = ( lambda / 0.197 * sqrt( s ) ); // lambda in fm, sqrt(s) in GeV, lambda_scaled dimensionless
	if ( log( lambda_scaled ) >= 15.0 )
	{
	  cout << log( lambda_scaled ) << "   " << lambda << " fm  " << sqrt( s ) <<  endl;
	}

	initialStateIndex = -1;
	scatt23_object.setParameter( VectorXYZ(), P1, P2, F1, F2, m1, m2, sqrt( s ), 
			  md2g_wo_as / s, lambda_scaled, 
			  theConfig.getK23LightPartons(), theConfig.getK23HeavyQuarks(),
			  theConfig.getKappa23LightPartons(), theConfig.getKappa23HeavyQuarks(),
			  theConfig.I23onlineIntegrationIsSet(),
			  theConfig.get23GluonFormationTimeTyp(), theConfig.getMatrixElement23(), theConfig.isMd2CounterTermInI23(), theConfig.get23FudgeFactorLpm(),
			  nGluons, theConfig.isMatrixElement23_22qt() );
	cs23 = scatt23_object.getXSection23( initialStateIndex ); //1/GeV^2 
	
	vSigma.add( c23, F1, F2, vRel * cs23 );
      }
      //-------------------------------

      // 3->2 ------------------------
      if( theConfig.doScattering_32() )
      {
	s = ( part1.Mom + part2.Mom + part3.Mom ).M2();
	lambda_scaled = ( lambda / 0.197 * sqrt( s ) ); // lambda in fm, sqrt(s) in GeV, lambda_scaled dimensionless

	initialStateIndex = -1;
	scatt32_object.setParameter( VectorXYZ(), part1.Mom, part2.Mom, part3.Mom, F1, F2, F3, sqrt( s ), md2g_wo_as * coupling::get_constant_coupling() / s, lambda_scaled, coupling::get_constant_coupling(), theConfig.isMd2CounterTermInI23(), theConfig.isMatrixElement23_22qt(), theConfig.get23FudgeFactorLpm(), -1 );
	I32 = scatt32_object.getIntegral32_withPrefactors( initialStateIndex );

	double tempI32 = I32 / pow( 0.197, 5.0 );
	vSigma.add( c32, F1, F2, F3, tempI32 );
      }
      //-------------------------------
    }

    vSigma.normalizeRates( getGluonDensity( T, gluonFugacity, quantumStatistics ), getQuarkDensity( T, quarkFugacity, quantumStatistics ) / 2.0, getQuarkDensity( T, quarkFugacity, quantumStatistics ) / 2.0 );
    
    // in former versions the energy of a sampled particle was used. Since this gives a bias now the mean energy of a thermal particle is used. Influences only the HF.
    velocity = sqrt( 1.0 - pow( ParticlePrototype::getMass( flav ), 2.0 ) / pow( getMeanThermalEnergy( T, ( quantumStatistics? fermi: boltzmann ), ParticlePrototype::getMass( flav ) ), 2.0 ) );
    // lambda which ones obtain with the previous employed lambda for the cross section
    lambda = vSigma.getLambda( flav, velocity, fm ); //fm

//     cout << lambda << " fm    " <<  lambda / 0.197  << " GeV^-1  " << endl;

    lambdaArray.push_back( lambda );
    ratesArray.push_back( vSigma );

    if ( lambdaArray.size() >= 4 )
    {
      lambdaAvr = 0;
      for ( int m = ( int( lambdaArray.size() ) - 1 ); m >= ( int( lambdaArray.size() ) - 4 ); m-- )
        lambdaAvr += lambdaArray[m];
      lambdaAvr = lambdaAvr / 4;

      if (( fabs( lambdaArray[lambdaArray.size()-1] - lambdaAvr ) / lambdaAvr < epsilon ) &&
          ( fabs( lambdaArray[lambdaArray.size()-2] - lambdaAvr ) / lambdaAvr < epsilon ) &&
          ( fabs( lambdaArray[lambdaArray.size()-3] - lambdaAvr ) / lambdaAvr < epsilon ) &&
          ( fabs( lambdaArray[lambdaArray.size()-4] - lambdaAvr ) / lambdaAvr < epsilon ) )
        converged = true;
      else
        converged = false;
    }
    else
      converged = false;

    if ( lambdaArray.size() >= 6 )
      lambda = lambdaAvr; // fm
  }
  while ( !converged && lambdaArray.size() < nIterationsMax );

  // in former versions the energy of a sampled particle was used. Since this gives a bias now the mean energy of a thermal particle is used. Influences only the HF.
  if(converged)
    cout << "converged thermal MFP: lambda = " << lambda << " flavor = " << flav << " thermal E = " << getMeanThermalEnergy( T, ( quantumStatistics? fermi: boltzmann ), ParticlePrototype::getMass( flav ) ) << " GeV" << endl;
  else
    cout << "not converged thermal MFP: converged: lambda = " << lambda << " flavor = " << flav << " thermal E = " << getMeanThermalEnergy( T, ( quantumStatistics? fermi: boltzmann ), ParticlePrototype::getMass( flav ) ) << " GeV" << endl;
  
  lambdaAvr = 0;
  int stop = static_cast<int>( lambdaArray.size() ) - 4;
  for ( int m = lambdaArray.size() - 1; m >= stop; m-- )
  {
    lambdaAvr += lambdaArray[m];
    iteratedRates += ratesArray[m];
  }
  lambda = lambdaAvr / 4;
  iteratedRates /= 4;
  
  return lambda; // fm
}


/**
* Determines mean free path of a gluon with only binary interactions
*
* @param[in] theI22
* @param[in] T temperature in GeV
* @param[in] theConfig
* @return mean free path in fm
*/
double get22mfp( const double T, interpolation22& theI22, const configBase& theConfig, const int nRuns22, const double gluonFugacity, const double quarkFugacity, const bool quantumStatistics )
{
  ParticlePrototype partA, partB;
  
  double vSigma22_sum = 0;
  double vSigma22_g_sum = 0;
  double vSigma22_q_sum = 0;
  
  // Debye mass not yet multiplied with coupling alpha_s because this is done in the scattering routines since the scale depends on the process
  double md2g_wo_as, md2q_wo_as;
  if( quantumStatistics )
  {
    // Bose definition
    md2g_wo_as = ( 1.0 + ParticlePrototype::N_light_flavor / 6.0 ) * 4.0 * M_PI * pow( T, 2.0 ); //GeV^2 TODO: how to consider here fugacities out of eq.? 
    md2q_wo_as = 0; //GeV^2, What is the quantum defintion of quark screening mass?
  }
  else
  {
    // Boltzmann definition
    md2g_wo_as = ( gluonFugacity * ns_casc::Ncolor + quarkFugacity * ParticlePrototype::N_light_flavor ) * 8 / M_PI * pow( T, 2.0 ); //GeV^2
    md2q_wo_as = 16 / ( 3 * M_PI ) * ( ( gluonFugacity * ns_casc::Ncolor + quarkFugacity * ParticlePrototype::N_light_flavor ) / 2.0 ) * pow( T, 2.0 ); //GeV^2
  }      
  
  double s;
  VectorEPxPyPz P1, P2;
  double M1, M2;
  double vRel, cs22, vSigma22;
  FLAVOR_TYPE F1, F2;
  double mfp, velocity, rate, rate_g, rate_q;
  int n_coll_g = 0, n_coll_q = 0;

  int initialStateIndex = -1;
  
  scattering22 scatt22_object( &theI22 );

  for ( int i = 0; i < nRuns22; i++ )
  {
    initThermal( partA, gluon, T, quantumStatistics );
    initThermal( partB, T, gluonFugacity, quarkFugacity, quantumStatistics );

    F1 = partA.FLAVOR;
    M1 = partA.m;
    P1 = partA.Mom;

    F2 = partB.FLAVOR;
    M2 = partB.m;
    P2 = partB.Mom;

    s = (P1+P2).M2();
    vRel = VelRel(P1,P2, M1,M2); // general relative velocity

    if ( s > 0.1 )
    {
      scatt22_object.setParameter( P1, P2, F1, F2, M1, M2, s, vRel, md2g_wo_as , md2q_wo_as,
                                       theConfig.getKggQQb(), theConfig.getKgQgQ(), theConfig.getKappa_gQgQ(), 
                                       theConfig.isConstantCrossSecGQ(),
                                       theConfig.getConstantCrossSecValueGQ(), theConfig.isIsotropicCrossSecGQ() ); // md2g_wo_as, md2q_wo_as are debye masses without the factor alpha_s which is multiplied in scattering22.cpp
      cs22 = scatt22_object.getXSection22( initialStateIndex );
    }
    else
      cs22 = 0.0;

    vSigma22 = vRel * cs22;
    
    if( std::isnan( vSigma22 ) )
      cout << "vsigma22 is nan  " << s << endl;
    
    vSigma22_sum += vSigma22;
    
    if( partB.FLAVOR == gluon )
    {
      vSigma22_g_sum += vSigma22;
      n_coll_g++;
    }
    else if( ParticlePrototype::mapToGenericFlavorType( partB.FLAVOR ) == light_quark || ParticlePrototype::mapToGenericFlavorType( partB.FLAVOR ) == anti_light_quark )
    {
      vSigma22_q_sum += vSigma22;
      n_coll_q++;
    }
    else
    {
      string errMsg = "Wrong flavor sampled.";
//       throw eSimulate_error( errMsg );
    }
  }

  if( n_coll_g + n_coll_q != nRuns22 )
  {
    cout << n_coll_g << "  " << n_coll_q << "  " << n_coll_g+n_coll_q << "  " << nRuns22 << endl;
    string errMsg = "Numbers of collisions do not match.";
//     throw eSimulate_error( errMsg );
  }
  
  if( n_coll_g != 0 )
    rate_g = pow(0.197,2.0) * getGluonDensity( T, gluonFugacity, quantumStatistics ) * vSigma22_g_sum / n_coll_g;  // fm^-1
  else
    rate_g = 0.0;
  if( n_coll_q != 0 )
    rate_q = pow(0.197,2.0) * getQuarkDensity( T, quarkFugacity, quantumStatistics ) * vSigma22_q_sum / n_coll_q;  // fm^-1
  else
    rate_q = 0.0;
  rate = rate_g + rate_q;  // fm^-1
  
  velocity = sqrt( 1.0 - pow( partA.m ,2.0) / pow( partA.Mom.E() ,2.0) );
  // lambda which ones obtain with the previous employed lambda for the cross section 
  mfp = velocity / rate; // fm
  
  cout << "thermal mfp = " << mfp << " fm" << endl;

  return mfp; // fm
}
