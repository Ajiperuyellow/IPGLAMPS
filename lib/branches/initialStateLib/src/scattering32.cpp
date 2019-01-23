//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/scattering32.cpp $
//$LastChangedDate: 2014-12-11 22:00:00 +0100 (Do, 11. Dez 2014) $
//$LastChangedRevision: 2015 $
//$LastChangedBy: gallmei $
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------

//-----------------------------------------------------------
// parameters for the integration via the Cuba routines
#define EPSREL 1e-2
#define EPSABS 1e-3
#define VERBOSE 0
//bits 0&1: output level, from 0 (no output) to 3 (lots of information)
//bit 2: 0 = all sets of samples collected on a subregion contribute, 1 = only the last and largest sample contributes
//bit 3: 0 = Sobol quasi-random numbers are used, 1 = Mersenne-Twister pseudo-random numbers are used

#define MINEVAL 0
#define MAXEVAL 600

// Vegas
#define NSTART 50
#define NINCREASE 100
#define NBATCH 1000

// Suave
#define NNEW 200
#define FLATNESS 50

// Divonne
#define KEY1 -1
#define KEY2 -1
#define KEY3 0
#define MAXPASS 2
#define BORDER 0
#define MAXCHISQ 10
#define MINDEVIATION 0.25
#define NEXTRA 0
#define PEAKFINDER 0
//-----------------------------------------------------------


#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include <boost/bind.hpp>
#include <boost/function.hpp>

#include "scattering32.h"
#include "integrate.h"
#include "lorentz.h"
#include "random.h"
#include "FPT_compare.h"
#include "prefactors23.h"
#include "coupling.h"
#include "globalsettings.h"
//#include "interactiontype.h"  //needed for the create_vector<T> routines used in initializing mapOrderForRandomizedInitialState
//#include "distributions.h"  //only needed when a gaussian proposal function is used in scattering32::getMomenta32_metropolis(..)

// include config file that is created at CMake configuration time 
#include "configBAMPS.h"


//-----------------------------------------------------------
// Select the method for numerical integration
#ifndef CUBA_FOUND  // default fallback in case the Cuba library is not available
  #define NR_VEGAS  
#else               // select a method in case the Cuba library is available
  // Only one of the choices below must be defined!
  #define NR_VEGAS
  // #define CUBA_VEGAS
  // #define CUBA_SUAVE
  // #define CUBA_DIVONNE
#endif
//-----------------------------------------------------------


using namespace std;



/**
 * This standard constructor sets all pointers to NULL and sets the standard values
 */
scattering32::scattering32() : 
  beta_abs(),
  F1(), F2(), F3(),
  alpha(), md2g_scaled(), sqrtS(), lambda_scaled(),
  numberOfGluons( 0 ),
  E1_selected(), E3_selected(),
  cos_gamma(), 
  absorpedGluon(),
  md2_use_simplified_GB(),
  matrixElement32_22qt(),
  fudge_factor_lpm()
{
}



/**
 *
 * @param v vector of the collective velocity of the cell
 * @param P1_arg momentum vector of particle 1
 * @param P2_arg momentum vector of particle 2
 * @param P3_arg momentum vector of particle 3
 * @param sqrtS_arg center of mass energy, sqrt(s)
 * @param md2g_scaled_arg squared debye mass, scaled by s
 * @param lambda_scaled_arg mean free path, scaled by sqrt(s)
 */
scattering32::scattering32( const VectorXYZ & v,
                            const VectorEPxPyPz & P1_arg, 
                            const VectorEPxPyPz & P2_arg, 
                            const VectorEPxPyPz & P3_arg,
                            const FLAVOR_TYPE F1_arg, 
                            const FLAVOR_TYPE F2_arg, 
                            const FLAVOR_TYPE F3_arg,
                            const double sqrtS_arg, 
                            const double md2g_scaled_arg, 
                            const double lambda_scaled_arg, 
                            const double _alpha_s, 
                            const bool _md2_use_simplified_GB,
                            const bool _matrixElement32_22qt, 
                            const double fudge_factor_lpm_arg, 
                            const int _Ng ) : 
  P1( P1_arg ), P2( P2_arg ), P3( P3_arg ), 
  F1( F1_arg ), F2( F2_arg ), F3( F3_arg ), 
  alpha( _alpha_s ),
  md2g_scaled( md2g_scaled_arg ), 
  sqrtS( sqrtS_arg ), 
  lambda_scaled( lambda_scaled_arg ), 
  numberOfGluons( _Ng ),
  md2_use_simplified_GB( _md2_use_simplified_GB), 
  matrixElement32_22qt( _matrixElement32_22qt ), 
  fudge_factor_lpm( fudge_factor_lpm_arg )
{
  do
  {
    double selectGluon = ran2() * 3;
    
    if ( selectGluon < 1.0 )
    {
      F3 = F1_arg;
      absorpedGluon = 1;
    }
    else if ( selectGluon < 2.0 )
    {
      F3 = F2_arg;
      absorpedGluon = 2;
    }
    else
    {
      F3 = F3_arg;
      absorpedGluon = 3;
    }
  } while( F3 != gluon );
  
  switch( absorpedGluon )
  {
  case 1:
    F1 = F2_arg;
    F2 = F3_arg;
    P1 = P2_arg;
    P2 = P3_arg;
    P3 = P1_arg;
    break;
  case 2:
    F1 = F1_arg;
    F2 = F3_arg;
    P1 = P1_arg;
    P2 = P3_arg;
    P3 = P2_arg;
    break;
  case 3:
    F1 = F1_arg;
    F2 = F2_arg;
    P1 = P1_arg;
    P2 = P2_arg;
    P3 = P3_arg;
    break;
  default:
    std::string errMsg = "Error in scattering32::scattering32(..). absorpedGluon did not resolve to valid value.";
    throw eScatt32_error( errMsg );
  }
    
  E1_selected = 0.0;
  E3_selected = 0.0;
  cos_gamma = 0.0;

  // boost momenta to rest frame of the cell
  LL_cell.setBeta( v );
  LL_cell.boost(P1,P2,P3, P1cell,P2cell,P3cell);

  // boost momenta to the centre of mass system of the colliding particles
  LL_CM.setBetaCM(P1cell,P2cell,P3cell);
  LL_CM.boost(P1cell,P2cell,P3cell, P1cm,P2cm,P3cm);

  beta_abs = LL_CM.betaVal();

  // transform (rotate) beta into the refernce frame given by p1 and p3
  rotate_beta( P1cm, P3cm, rotated_beta );
}


scattering32::~scattering32()
{
}



/**
 * This method sets all necessary parameters for a given particle triplet. Previous values are deleted or overwritten.
 * Using this method an scattering32 object can be re-used for multiple particle triplets, thus reducing the need to constantly
 * creating new objects.
 *
 * Either this method or the constructor taking the same arguments MUST be called prior to any other methods of the class!
 *
 * @param v vector of the collective velocity of the cell
 * @param P1_arg momentum vector of particle 1
 * @param P2_arg momentum vector of particle 2
 * @param P3_arg momentum vector of particle 3
 * @param sqrtS_arg center of mass energy, sqrt(s)
 * @param md2g_scaled_arg squared debye mass, scaled by s
 * @param lambda_scaled_arg mean free path, scaled by sqrt(s)
 * @return The absolute value of the boost velocity between cell frame and CMS (#beta_vec).
 */
double scattering32::setParameter( const VectorXYZ & v,
                                   const VectorEPxPyPz & P1_arg, 
                                   const VectorEPxPyPz & P2_arg, 
                                   const VectorEPxPyPz & P3_arg,
                                   const FLAVOR_TYPE F1_arg, 
                                   const FLAVOR_TYPE F2_arg, 
                                   const FLAVOR_TYPE F3_arg,
                                   const double sqrtS_arg, 
                                   const double md2g_scaled_arg, 
                                   const double lambda_scaled_arg,
                                   const double _alpha_s, 
                                   const bool _md2_use_simplified_GB,
                                   const bool _matrixElement32_22qt, 
                                   const double fudge_factor_lpm_arg, 
                                   const int _Ng )
{
  alpha = _alpha_s;
  numberOfGluons = _Ng;

  md2_use_simplified_GB = _md2_use_simplified_GB;
  matrixElement32_22qt = _matrixElement32_22qt;
  fudge_factor_lpm = fudge_factor_lpm_arg;
  
  do
  {
    double selectGluon = ran2() * 3;
    
    if ( selectGluon < 1.0 )
    {
      F3 = F1_arg;
      absorpedGluon = 1;
    }
    else if ( selectGluon < 2.0 )
    {
      F3 = F2_arg;
      absorpedGluon = 2;
    }
    else
    {
      F3 = F3_arg;
      absorpedGluon = 3;
    }
  } while( F3 != gluon );
  
  switch( absorpedGluon )
  {
    case 1:
      F1 = F2_arg;
      F2 = F3_arg;
      P1 = P2_arg;
      P2 = P3_arg;
      P3 = P1_arg;
      break;
    case 2:
      F1 = F1_arg;
      F2 = F3_arg;
      P1 = P1_arg;
      P2 = P3_arg;
      P3 = P2_arg;
      break;
    case 3:
      F1 = F1_arg;
      F2 = F2_arg;
      P1 = P1_arg;
      P2 = P2_arg;
      P3 = P3_arg;
      break;
    default:
      std::string errMsg = "Error in scattering32::scattering32(..). absorpedGluon did not resolve to valid value.";
      throw eScatt32_error( errMsg );
  }

  sqrtS = sqrtS_arg;
  md2g_scaled = md2g_scaled_arg;
  lambda_scaled = lambda_scaled_arg;

  E1_selected = 0.0;
  E3_selected = 0.0;
  cos_gamma = 0.0;

  // boost momenta to rest frame of the cell
  LL_cell.setBeta( v );
  LL_cell.boost(P1,P2,P3, P1cell,P2cell,P3cell);


  // boost momenta to the centre of mass system of the colliding particles
  LL_CM.setBetaCM(P1cell,P2cell,P3cell);
  LL_CM.boost(P1cell,P2cell,P3cell, P1cm,P2cm,P3cm);

  beta_abs = LL_CM.betaVal();

  // transform (rotate) beta into the refernce frame given by p1 and p3
  rotate_beta( P1cm, P3cm, rotated_beta );

  return beta_abs;
}


/**
 * scattering32::getIntegral32_bare(..) provides the integrated matrix element for 3->2 processes WITHOUT any prefactors
 * It is:
 * I32bare = \f$\int_{0}^{1} d\cos(\theta) \int_{0}^{\pi} d\phi
 *         \frac{ q_t^2 }{ (q_t^2 + m_D^2)^2  k_t^2 [k_t^2 + q_t^2 - 2 \vec{k}_t\vec{q}_t + m_D^2 ] }\f$        (1)
 * 
 * I32bare is either computed via numerical (Monte Carlo) integration or via an effective method that is based on sampling a 
 * estimation function and comparing to the true matrix element (essentially some sort of rejection sampling).
 * 
 * @param[out] initialStateIndex index that characterizes the initial state, i.e. whether it is gg -> X, qq -> X, etc.
 * @param[in] _compType Which computation method to use
 * @return Integral over the 3->2 matrix element without any prefactors
 */
double scattering32::getIntegral32_bare(int& initialStateIndex, const I32_COMPUTATION_TYPE _compType)
{
  double result = 0;
  
  switch( _compType )
  {
  case FAST:
    result = getIntegral32_bare_fast( initialStateIndex );
    break;
  case MONTE_CARLO_INTEGRATION:
    result = getIntegral32_bare_vegas( initialStateIndex );
    break;
  default:
    std::string errMsg = "Error in scattering32::getIntegral32_bare(). _compType did not resolve to valid value.";
    throw eScatt32_error( errMsg );
  }   
  
  return result;
}



/**
 * scattering32::getIntegral32_withPrefactors(..) provides the integrated matrix element I32 for 3->2 process (see scattering32::getIntegral32_bare) 
 * INCLUDING all the prefactors. This is done in such a way that the probability for a given 3->2 process can be obtained from I32 via
 * p32 = I32 * dT / ( dV^2 * Ntest )
 * 
 * @param[out] initialStateIndex index that characterizes the initial state, i.e. whether it is gg -> X, qq -> X, etc.
 * @param[in] _compType Which computation method to use
 * @return Integral over the 3->2 matrix element including all prefactors
 */
double scattering32::getIntegral32_withPrefactors(int& initialStateIndex, const I32_COMPUTATION_TYPE _compType)
{
  double result = 0;
  double prefactors = 0;
  double I32_bare = getIntegral32_bare( initialStateIndex, _compType );
  
//   prefactors = pow( 0.197, 5.0 ) * 9.0 * M_PI * ns_casc::Ncolor * pow( alpha, 3.0 ) / ( 2.0 * ( pow( ns_casc::Ncolor, 2.0 ) - 1 ) * pow( sqrtS, 2.0 ) * P1[0] * P2[0] * P3[0] ) ;
  prefactors = pow( 0.197, 5.0 ) * 27.0 * M_PI / ( 16.0 * pow( sqrtS, 2.0 ) * P1.E() * P2.E() * P3.E() ) ; // optimized version
  
  result = prefactors * I32_bare;
  
  return result;
}



/**
 * scattering32::getIntegral32(..) provides the integration of the matrix element for 3->2 processes
 * It is:
 * I32 = \f$\int_{0}^{1} d\cos(\theta) \int_{0}^{\pi} d\phi
 *         \frac{ q_t^2 }{ (q_t^2 + m_D^2)^2  k_t^2 [k_t^2 + q_t^2 - 2 \vec{k}_t\vec{q}_t + m_D^2 ] }\f$        (1)
 *
 *
 * The following angles are relevant for the calculation of I32. See notes for more details!
 * - p1,p2,p3 are the initial momentum vectors, p1' and p2' are the final momentum vectors
 * - vector p1 parallel to z-axis
 * - vectors p1,p2,p3 are located in the x-z-plane
 * - theta is the angle between p1' and p1
 * - gamma is the angle between p3 and p1
 * - phi is the angle between e_{x} and p1'*e_{x}
 * - the vector p1' can be decomposed into p1' = sqrt(s)/2 * ( sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta) )
 * - p1'*p3 = sqrt(s)/2*E3*cos(delta)
 * - delta is the angle between p1' and p3, cos(delta) = sin(gamma)*sin(theta)*cos(phi) + cos(gamma)*cos(theta)
 *
 * The LPM cutoff requires that
 *
 * kt^2 > 1/lambda * E3 * ( 1 + beta' * sqrt(s)/2 * cos(delta) * cos(Theta) )    (2)
 *
 * where Theta (capital!) is the angle between the boost-vector beta' (here: beta) and the vector p1'. For comparison with 2->3 note that
 * cosh(y) = E3 / kt,  tanh(y) = p1'*p3/E3 = sqrt(s)/2*cos(delta)  and  cos(Theta) = beta'*p1'/abs(beta'*p1').
 * Furthermore E3 > 1/lambda must be fulfilled, as E3 > kt and the gamma-factor (not the angle gamma) is larger than 1.
 * See notes for more details!
 *
 * theIntegrand is a function object of type integrand32 derived from type integrand. It can be passed to the VEGAS integration routine
 * simply via vegas(2, theIntegrand, &tgral, &sd, &chi2a).
 * Auxiliary members of integrand32 (E1_int, E3_int, cos_gamma_int, beta_vec, beta_abs, md2_int, sqrtS_int) are set
 * via integrand32::set_XX(..) prior to integration via vegas.
 * Results are stored in I32[..].
 *
 * @return Integral over the 3->2 matrix element I32.
 */
double scattering32::getIntegral32_bare_vegas( int& initialStateIndex )
{
  if( F1 > 2 * ParticlePrototype::max_N_light_flavor || F2 > 2 * ParticlePrototype::max_N_light_flavor || F3 > 2 * ParticlePrototype::max_N_light_flavor ) // heavy quark involved
  {
    return 0.0;
  }
  
  theIntegrand.set_md2( md2g_scaled );
  theIntegrand.set_sqrtS( sqrtS );
  theIntegrand.set_matrixElement32_22qt( matrixElement32_22qt );
  theIntegrand.set_md2_use_simplified_GB( md2_use_simplified_GB );

  // set lambda for use in the integrand
  theIntegrand.set_lambda( lambda_scaled );
  theIntegrand.set_fudge_factor_lpm( fudge_factor_lpm );
  
  E1_selected = P1cm(0) / sqrtS;          // CMS energy of particle 1 scaled by sqrt(s)
  E3_selected = P3cm(0) / sqrtS;          // CMS energy of particle 3 scaled by sqrt(s)

  // calculate cos(gamma) = p1*p3 / (abs(p1)*abs(p3))
  cos_gamma = CosTheta(P1cm,P3cm);


  //--------------------------------------------------------------
  // parameters needed for calls to the integration routines
  int neval;  // actual number of integrand evaluations needed
  int fail;   // 0 = desired accuracy was reached, 1 = accuracy not reached, -1 = dimension out of range
  double intResult[NCOMP_32]; // result of the integration over the unit hypercube, NCOMP_32 = #components, 1 in our case
  double error[NCOMP_32];     // presumed absolute error of integral
  double prob[NCOMP_32];      // xi^2 probability that error is NOT a reliable estimate of the true integration error
  //--------------------------------------------------------------

  //--------------------------------------------------------------
  // create the functionoid that handles the integration
  // will be called later with: integrate( theIntegrand, neval, fail, intResult, error, prob );
  // the integration routine is chosen via precompiler defined switches
#ifdef CUBA_VEGAS
  integrate_vegas integrate( NDIM_32, NCOMP_32, EPSREL, EPSABS, VERBOSE, MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH );
#endif

#ifdef CUBA_SUAVE
  integrate_suave integrate( NDIM_32, NCOMP_32, EPSREL, EPSABS, VERBOSE, MINEVAL, MAXEVAL, NNEW, FLATNESS );
#endif

#ifdef CUBA_DIVONNE
  const int ngiven = 0;
  const int ldxgiven = 0;
  double * xgiven;
  int nregions = -1;
  integrate_divonne integrate( NDIM_32, NCOMP_32, EPSREL, EPSABS, VERBOSE, MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
                               BORDER, MAXCHISQ, MINDEVIATION, ngiven, ldxgiven, xgiven, NEXTRA, PEAKFINDER, nregions );
#endif

#ifdef NR_VEGAS
  integrate_nr_vegas integrate( NDIM_32 );
#endif
  //--------------------------------------------------------------

  double I32 = 0;

  theIntegrand.set_E1( E1_selected );
  theIntegrand.set_E3( E3_selected );
  theIntegrand.set_cos_gamma( cos_gamma );
  theIntegrand.set_beta_vec( rotated_beta );

  integrate( theIntegrand, neval, fail, intResult, error, prob );
  I32 = intResult[0];

  // sort flavors of outgoing particles for the processing of prefactors for different processes
  unsigned int _F1 = std::min( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
  unsigned int _F2 = std::max( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );

  double total_prefactor = 0;
  initialStateIndex = -1;

  if (( _F1 + _F2 ) == 0 ) // g+g+g -> g+g, g+g+g -> q+qbar
  {
    prefactor23_gg_ggg pObj1;
    prefactor23_gg_qqbarg pObj2;
    
    initialStateIndex = 0;

    total_prefactor = pObj1.prefactor() * pObj1.symFactor32() + pObj2.prefactor() * pObj2.symFactor32();
  }
  else if ( _F1 == _F2 )  // q+q+g -> q+q, qbar+qbar+g -> qbar+qbar
  {
    prefactor23_qq_qqg pObj1;
    
    initialStateIndex = 1;

    total_prefactor = pObj1.prefactor() * pObj1.symFactor32();
  }
  else if (( _F1 * _F2 ) == 0 ) // g+q+g -> g+q, g+qbar+g -> g+qbar
  {
    prefactor23_qg_qgg pObj1;
    
    initialStateIndex = 2;

    total_prefactor = pObj1.prefactor() * pObj1.symFactor32();
  }
  else if (( _F2 - _F1 ) == 1 && ( _F2 % 2 ) == 0 )  // q+qbar+g -> q+qbar, q+qbar+g -> q'+qbar', q+qbar+g -> g+g
  {
    prefactor23_qqbar_qqbarg pObj1;
    prefactor23_qqbar_qqbarDashg pObj2;
    prefactor23_qqbar_ggg pObj3;
    
    initialStateIndex = 3;

    total_prefactor = pObj1.prefactor() * pObj1.symFactor32() + pObj2.prefactor() * pObj2.symFactor32() + pObj3.prefactor() * pObj3.symFactor32();
  }
  else // q+q'+g -> q+q', q+qbar'+g -> q+qbar'
  {
    prefactor23_qqdash_qqdashg pObj1;
    
    initialStateIndex = 4;

    total_prefactor = pObj1.prefactor() * pObj1.symFactor32();
  }

  return I32 * total_prefactor;
}




/**
 * The integral over the matrix element, I32, is estimated using an upper limit of the real matrix element, which can be
 * integrated analytically. The theta-function is dealt with by finding areas of integration such that the resulting integral
 * is larger than the integral with the theta-function. See notes.
 *
 * Using this method should be faster, though likely less accurate, than employing direct numerical integration
 * via #getIntegral32_vegas.
 *
 * Thus the estimated I32 is larger than the real I32. Using a number of sampled points the method #get_I32estimate_ratio returns
 * an estimate of the ratio of the real to the estimated I32. This is used to correct I32 before it is returned.
 *
 * @return An estimate for I32, corrected by the estimated ratio to the real result
 */
double scattering32::getIntegral32_bare_fast( int& initialStateIndex )
{
  if( F1 > 2 * ParticlePrototype::max_N_light_flavor || F2 > 2 * ParticlePrototype::max_N_light_flavor || F3 > 2 * ParticlePrototype::max_N_light_flavor ) // heavy quark involved
  {
    return 0.0;
  }
  
  const double lower_limit = 1.0e-4;

  E1_selected = P1cm(0) / sqrtS;          // CMS energy of particle 1 scaled by sqrt(s)
  E3_selected = P3cm(0) / sqrtS;          // CMS energy of particle 3 scaled by sqrt(s)

  // calculate cos(gamma) = p1*p3 / (abs(p1)*abs(p3))
  cos_gamma = CosTheta(P1cm,P3cm);

  double I32_estimate = 0;   // will hold the estimated integration result
  double ratio = 0;          // will hold the estimated ratio of the estimated restult to the real result
  double I32 = 0;            // will hold the corrected estimate of the integration result

  I32_estimate = get_I32estimate();
  if ( I32_estimate > lower_limit )
  {
    ratio = get_I32estimate_ratio();
  }
  
  I32 = I32_estimate * ratio * pow( coupling::get_maximum_coupling() , 3.0 );

  double total_prefactor = 0;
  if ( I32 > 0 )
  {
    // sort flavors of outgoing particles for the processing of prefactors for different processes
    unsigned int _F1 = std::min( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
    unsigned int _F2 = std::max( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );

    if (( _F1 + _F2 ) == 0 ) // g+g+g -> g+g, g+g+g -> q+qbar
    {
      prefactor23_gg_ggg pObj1;
      prefactor23_gg_qqbarg pObj2;
      
      initialStateIndex = 0;

      total_prefactor = pObj1.prefactor() * pObj1.symFactor32() + pObj2.prefactor() * pObj2.symFactor32();
    }
    else if ( _F1 == _F2 )  // q+q+g -> q+q, qbar+qbar+g -> qbar+qbar
    {
      prefactor23_qq_qqg pObj1;
      
      initialStateIndex = 1;

      total_prefactor = pObj1.prefactor() * pObj1.symFactor32();
    }
    else if (( _F1 * _F2 ) == 0 ) // g+q+g -> g+q, g+qbar+g -> g+qbar
    {
      prefactor23_qg_qgg pObj1;
      
      initialStateIndex = 2;

      total_prefactor = pObj1.prefactor() * pObj1.symFactor32();
    }
    else if (( _F2 - _F1 ) == 1 && ( _F2 % 2 ) == 0 )  // q+qbar+g -> q+qbar, q+qbar+g -> q'+qbar', q+qbar+g -> g+g
    {
      prefactor23_qqbar_qqbarg pObj1;
      prefactor23_qqbar_qqbarDashg pObj2;
      prefactor23_qqbar_ggg pObj3;
      
      initialStateIndex = 3;
      
      total_prefactor = pObj1.prefactor() * pObj1.symFactor32() + pObj2.prefactor() * pObj2.symFactor32() + pObj3.prefactor() * pObj3.symFactor32();
    }
    else // q+q'+g -> q+q', q+qbar'+g -> q+qbar', qbar+qbar'+g-> qbar+qbar'
    {
      prefactor23_qqdash_qqdashg pObj1;
      
      initialStateIndex = 4;

      total_prefactor = pObj1.prefactor() * pObj1.symFactor32();
    }
  }
  else
  {
    initialStateIndex = -1;
  }
  
  return I32 * total_prefactor;
}



/**
 * For a given choice (p1,p3) of particles 1 and 3 this routine calculates an upper estimate of the integral I32.
 *
 * @return Estimated I32 for given choice of particle 1 and 3
 */
double scattering32::get_I32estimate() const
{
  double cosgamma = fabs( cos_gamma );
  double singamma = sqrt( 1 - pow( cos_gamma, 2.0 ) );
  double C = 1 / ( lambda_scaled / fudge_factor_lpm * E3_selected ) * sqrt( ( 1 - beta_abs ) / ( 1 + beta_abs ) );
//   double C = 1 / ( pow( lambda_scaled, 2 ) * pow( E3_selected, 2 ) * ( 1 + beta_abs ) );

  double G1 = sqrt( 1 - C ) / cosgamma;
  double G3 = singamma;
  double G4 = sqrt( C ) * singamma + sqrt( 1 - C ) * cosgamma;

  double a = sqrt( 1 + md2g_scaled / pow( E1_selected, 2.0 ) );
  double alpha = E3_selected * fudge_factor_lpm / lambda_scaled * sqrt( ( 1 - beta_abs ) / ( 1 + beta_abs ) );  // some abbreviation
  double Q = 0;
  double F = 0;
  
  boost::function< double( const double, const double ) > H;
  
  if ( md2_use_simplified_GB )  // estimate the matrix element where terms have been cancelled out prior to applying the screening with md2
  {
    Q = ( 1 - 2 * sqrt( alpha ) + alpha ) / ( pow( E1_selected, 2 ) * md2g_scaled * alpha );
    //   double Q = 1 / ( pow( E1_selected, 2 ) * md2g_scaled * alpha );
    H = boost::bind( &scattering32::H_simplifiedGB, this, _1, _2 );
  }
  else // estimate the matrix element where the screening with md2 is applied without cancelling out terms
  {
    Q = ( 1 - 2 * sqrt( alpha ) + alpha );
    H = boost::bind( &scattering32::H_actualGB, this, _1, _2 );
  }
  
  if ( 1.0 > sqrt( C ) )
  {
    F += H( a, G3 ) - H( a, 0 );
  }

  if ( G1 > G3 )
  {
    F += H( a, min( 1.0, G1 ) ) - H( a, G3 );
  }

  if ( 1.0 > G1 && G4 > max( G1, G3 ) )
  {
    F += H( a, G4 ) - H( a, max( G1, G3 ) );
  }

  return M_PI * Q * F;
}


/**
 * Utility method used by # get_I32estimate. It calculates the antiderivative of the estimated matrix element at given u.
 * This version is used for estimating the simplified Gunion-Bertsch version where some terms are cancelled prior to applying
 * a screening via the Debye mass. 
 * 
 * @param a
 * @param u
 * @return Value of the antiderivative of the estimated matrix element at given u.
 */
double scattering32::H_simplifiedGB( const double a, const double u ) const
{
  double a2 = pow( a, 2.0 );
  double h = ( a2 + 1 ) / ( 4 * pow( a, 3.0 ) ) * log(( a + u ) / ( a - u ) ) - ( a2 - 1 ) * u / ( 2 * a2 * ( a2 - pow( u, 2.0 ) ) ) ;
  return h;
}


/**
 * Utility method used by # get_I32estimate. It calculates the antiderivative of the estimated matrix element at given u.
 * This version is used for estimating the actual Gunion-Bertsch result where the terms are NOT cancelled prior to applying
 * a screening via the Debye mass. 
 * 
 * @param a
 * @param u
 * @return Value of the antiderivative of the estimated matrix element at given u.
 */
double scattering32::H_actualGB( const double a, const double u ) const
{
  double md4 = pow( md2g_scaled, 2 );
  double E12 = pow( E1_selected, 2 );
  
  double A = 1 / ( md4 * E12 );
  double alpha = E3_selected * fudge_factor_lpm / lambda_scaled * sqrt( ( 1 - beta_abs ) / ( 1 + beta_abs ) );  // some abbreviation
  double B = ( alpha + md2g_scaled ) / ( alpha * md4 * pow( E1_selected, 3 ) );
  double C = ( alpha / 4 + md4 ) / ( alpha * md4 * pow( E12, 2 ) );

  double result = A * H1( a, u ) + B * H2( a, u ) + C * H3( a, u );
  
  return result;
}


/**
 * Utility method used by # get_I32estimate. It calculates the antiderivative of the second term of the estimated matrix element at given u.
 *
 * @param a
 * @param u
 * @return Value of the antiderivative of the second term of the estimated matrix element at given u.
 */
double scattering32::H2( const double a, const double u ) const
{
  double a2 = pow( a, 2 );
  double u2 = pow( u, 2 );
     
  double t1 = u * sqrt( 1 - u2 ) / ( 2 * a2 * (a2 - u2) );
  
  double t2;
  if ( FPT_COMP_E( u, 1.0 ) ) // catch special case where u = 1 and thus A -> inf which implies that ln( (1+iA)/(1-iA) ) = i * pi
  {
    t2 = M_PI / ( 4 * pow( a, 3 ) * sqrt( a2 - 1 ) );
  }
  else
  {
    double A = u * sqrt( ( a2 - 1 ) / ( 1 - u2 ) ) / a;
    double A2 = pow( A, 2 );
    double Re_z = ( 1 - A2 );   // actually it is Re_z = ( 1 - A2 ) / ( 1 + A2 ) but the common denominator of Re(z) and Im(z) does not matter for the computation of atan2
    double Im_z = 2 * A;        // actually it is Im_z = 2 * A / ( 1 + A2 ) but the common denominator Re(z) and Im(z) does not matter for the computation of atan2
    t2 = atan2( Im_z, Re_z ) / ( 4 * pow( a, 3 ) * sqrt( a2 - 1 ) );  // this computes ln(z) as ln(z) = ln(|z|) + i Arg(z) where |z| = 1 in this case and the i is cancelled via sqrt( 1 - a2 ) = i sqrt( a2 - 1 )
  }
  
  return t1 + t2;
}


/**
 * Utility method used by # get_I32estimate. It calculates the antiderivative of the third term of the estimated matrix element at given u.
 *
 * @param a
 * @param u
 * @return Value of the antiderivative of the third term of the estimated matrix element at given u.
 */
double scattering32::H3( const double a, const double u ) const
{
  double a2 = pow( a, 2 );
  double u2 = pow( u, 2 );
  
  double t1 = u / ( 2 * a2 * ( a2 - u2 ) );
  double t2 = log( ( a + u ) / ( a - u ) ) / ( 4 * pow( a, 3 ) );
  
  return t1 + t2;
}


/**
 * The upper estimate of the integrand when simplifying the GB terms prior to applying the screening mass. 
 * This is the function denoted as F-tilde in the handwritten notes. The prefactor Q is not contained here!
 * 
 * @return Upper estimate of the integrand
 */
double scattering32::F_estimate_simplifiedGB( const double a2, const double u ) const
{
  return ( 1 - pow( u, 2 ) ) / pow( a2 - pow( u, 2 ) , 2 );
}


/**
 * The upper estimate of the integrand when NOT simplifying the GB terms prior to applying the screening mass. 
 * This is the function denoted as F-tilde in the handwritten notes. The prefactor Q is not contained here!
 * 
 * @return Upper estimate of the integrand
 */
double scattering32::F_estimate_actualGB( const double a2, const double u ) const
{
  double md4 = pow( md2g_scaled, 2 );
  double E12 = pow( E1_selected, 2 );
  
  double A = 1 / ( md4 * E12 );
  double alpha = E3_selected * fudge_factor_lpm / lambda_scaled * sqrt( ( 1 - beta_abs ) / ( 1 + beta_abs ) );  // some abbreviation
  double B = ( alpha + md2g_scaled ) / ( alpha * md4 * pow( E1_selected, 3 ) );
  double C = ( alpha / 4 + md4 ) / ( alpha * md4 * pow( E12, 2 ) );
  
  double u2 = pow( u, 2 );
  
  double result = ( A * ( 1 - u2 ) + B * sqrt( 1 - u2 ) + C )  / pow( a2 - u2 , 2 );
  
  return result;
}



/**
 * Estimates the ratio of the real result for I32 to the estimated result for I32. This is done by sampling a certain number of points
 * according to the estimated matrix element (which is always larger than the real) and taking the ratio of the function values
 * at these points. An average over the so obtained values for the ratio should give a somewhat reliable result.
 *
 * @return Estimated ratio between real I32 and estimated I32
 */
double scattering32::get_I32estimate_ratio()
{
  const int samples = 60;
  
  double a2 = 1 + md2g_scaled / pow( E1_selected, 2.0 );
  double alpha = E3_selected * fudge_factor_lpm / lambda_scaled * sqrt( ( 1 - beta_abs ) / ( 1 + beta_abs ) );  // some abbreviation
  
  double u, phi , us;
  double g_estimate, g_compare;
  double ratio = 0;

  // the ranges in which the variables u and phi need to be sampled
  const double u_min = 0.0;
  const double u_max = 1.0;
  const double phi_min = 0.0;
  const double phi_max = M_PI;

  const double cosgamma = fabs( cos_gamma );
  const double singamma = sqrt( 1 - pow( cosgamma, 2.0 ) );
  double C = 1 / ( lambda_scaled / fudge_factor_lpm * E3_selected ) * sqrt( ( 1 - beta_abs ) / ( 1 + beta_abs ) );
//   double C = 1 / ( pow( lambda_scaled, 2 ) * pow( E3_selected, 2 ) * ( 1 + beta_abs ) );

  double Q = 0;
  double G = 0; // the upper estimate for the function F-tilde from which u is to be sampled
  
  boost::function< double( const double, const double ) > F;
  if ( md2_use_simplified_GB )  // estimate the matrix element where terms have been cancelled out prior to applying the screening with md2
  {
    Q = ( 1 - 2 * sqrt( alpha ) + alpha ) / ( pow( E1_selected, 2 ) * md2g_scaled * alpha );
    //   double Q = 1 / ( pow( E1_selected, 2 ) * md2g_scaled * alpha );
    
    // F-tilde = Q * (1-u^2) / (a^2 - u^2), see notes
    // the upper estimate G depends on the value of a
    if ( a2 <= 2 )
    {
      G = Q / ( 4 * ( a2 - 1 ) );
    }
    else
    {
      G = Q / pow( a2, 2 );
    }
    
    F = boost::bind( &scattering32::F_estimate_simplifiedGB, this, _1, _2 );
  }
  else // estimate the matrix element where the screening with md2 is applied without cancelling out terms
  {
    Q = ( 1 - 2 * sqrt( alpha ) + alpha );
    
    double md4 = pow( md2g_scaled, 2 );
    double E12 = pow( E1_selected, 2 );
    double A = 1 / ( md4 * E12 );
    double B = ( alpha + md2g_scaled ) / ( alpha * md4 * pow( E1_selected, 3 ) );
    double C = ( alpha / 4 + md4 ) / ( alpha * md4 * pow( E12, 2 ) );
    
    // it is difficult to find the local maximum for this F-tilde, therefore a simpler version F-tilde = Q * ( A + B + C ) / ( u^2 - a^2 )^2 is used that gives an upper limit
    // the maximum of this upper limit is at u = 0 and given by G = Q * ( A + B + C ) / a^4
    G = Q * ( A + B + C ) / pow( a2, 2 ); 
    
    F = boost::bind( &scattering32::F_estimate_actualGB, this, _1, _2 );
  }
  
  
  double sin_delta_2_limit;

  for ( int i = 0; i < samples; i++ )
  {
    do
    {
      u = u_min + ran2() * ( u_max - u_min );
      phi = phi_min + ran2() * ( phi_max - phi_min );
      us = sqrt( 1.0 - pow( u, 2.0 ) );

      // the upper estimate for sin(delta)^2 depends on the value of u = cos(theta) compared to sin(gamma)
      if ( u < singamma )
      {
        sin_delta_2_limit = 1;
      }
      else
      {
        sin_delta_2_limit = 1 - pow(( u * cosgamma - us * singamma ), 2.0 );
      }

      if ( sin_delta_2_limit < C )  // in this case the estimated Theta-function from the LPM cutoff can not be fulfilled
      {
        g_estimate = 0;
      }
      else
      {
        g_estimate = Q * F( a2, u );
      }

      g_compare = G * ran2();
    }
    while ( g_estimate < g_compare );

    double g_real = getMatrixElement( u, phi );
    // alpha_s^2(t) alpha_s(kt) is multiplied in getMatrixElement(), therefore, also multiply the maximum value alpha_s_max^3 tfor the estimate
    g_estimate = g_estimate * pow( coupling::get_maximum_coupling() , 3.0 );

    ratio += g_real / g_estimate;
  }

  ratio = ratio / samples;

  return ratio;
}




/**
 * Samples new momenta of the outgoing particles according to the matrix element.
 *
 * Wrapper, calls a certain implementation of the sampling in u and phi.
 *
 * @param[out] u u=cos(theta), theta = angle between p1' and p1
 * @param[out] phi azimuthal angle between e_x and (p1' e_x)
 * @return The "number" of the absorped gluon (1, 2 or 3) as counted internally corresponding to  F1, F2, F3 passed to setParameter or the constructor
 */
int scattering32::getMomenta32( double& u, double& phi, int& typ, FLAVOR_TYPE& F1arg, FLAVOR_TYPE& F2arg )
{
  getMomenta32_metropolis( u, phi ); // fast

  // slow (as the comparison function in use might be off by some orders of magnitude, leading to high rejection probabilities), use is discouraged
  // error = getMomenta32_rejection(u,phi);
  
  F1arg = F1;
  F2arg = F2;
   
  // sort F1 and F2 such that comparisons below are easier
  // convert to unsigned int to allow for some math that speeds up the decisions
  unsigned int _F1 = std::min( static_cast<unsigned int>( F1arg ), static_cast<unsigned int>( F2arg ) );
  unsigned int _F2 = std::max( static_cast<unsigned int>( F1arg ), static_cast<unsigned int>( F2arg ) );

  if (( _F1 + _F2 ) == 0 ) // g+g -> g+g+g, g+g -> q+qbar+g
  {
    prefactor23_gg_ggg pObj1;
    prefactor23_gg_qqbarg pObj2;
    double pf1 = pObj1.prefactor() * pObj1.symFactor32();
    double pf2 = pObj2.prefactor() * pObj2.symFactor32();

    if ( ran2() < pf1 / ( pf1 + pf2 ) )
    {
      typ = 321; // g+g+g -> g+g
    }
    else
    {
      typ = 322; // g+g+g -> q+qbar
      sampleFlavor23( F1arg, F2arg );
    }
  }
  else if ( _F1 == _F2 )  // q+q+g -> q+q, qbar+qbar+g -> qbar+qbar
  {
    typ = 327;
  }
  else if (( _F1 * _F2 ) == 0 ) // g+q+g -> g+q, g+qbar+g -> g+qbar
  {
    typ = 323;
  }
  else if (( _F2 - _F1 ) == 1 && ( _F2 % 2 ) == 0 )  // q+qbar+g -> q+qbar, q+qbar+g -> q'+qbar', q+qbar+g -> g+g
  {
    prefactor23_qqbar_qqbarg pObj1;
    prefactor23_qqbar_qqbarDashg pObj2;
    prefactor23_qqbar_ggg pObj3;

    double pf1 = pObj1.prefactor() * pObj1.symFactor32();
    double pf2 = pObj2.prefactor() * pObj2.symFactor32();
    double pf3 = pObj3.prefactor() * pObj3.symFactor32();

    double select = ran2() * ( pf1 + pf2 + pf3 );

    if ( select < pf1 )  // q+qbar+g -> q+qbar
    {
      typ = 324;
    }
    else if ( select < ( pf1 + pf2 ) )  // q+qbar+g -> q'+qbar'
    {
      typ = 325;
      sampleFlavor23( F1arg, F2arg, F1 );  // sample flavor excluding F1 (or anti-F1 = F2)
    }
    else  // q+qbar+g -> g+g
    {
      typ = 326;
      F1arg = gluon;
      F2arg = gluon;
    }
  }
  else // q+q'+g -> q+q', q+qbar'+g -> q+qbar'
  {
    typ = 328;
  }

  return absorpedGluon;
}



/**
 * Samples of new momenta according to the matrix element using the rejection method. Slow, since comparison function might be significantly
 * larger than the sampled function at some points. Use #getMomenta32_metropolis instead.
 *
 * @param[out] u_arg u=cos(theta), theta = angle between p1' and p1
 * @param[out] phi_arg azimuthal angle between e_x and (p1' e_x)
 * @return Count how many time the comparison function was below the actual function. Only needed for rejection method. Should be zero.
 */
int scattering32::getMomenta32_rejection( double& u_arg, double& phi_arg ) const
{
  int error = 0;  // see return value

  double u, phi;
  double g = 0, gr = 0;

  // the ranges in which the variables u and phi need to be sampled
  const double u_min = 0.0;
  const double u_max = 1.0;
  const double phi_min = 0.0;
  const double phi_max = M_PI;

  //double beta_abs = 0;   // for testing consistency with old implementation (beta_abs = 0 needs then to be set in integrand32.cpp also)
  double max = lambda_scaled / fudge_factor_lpm / E3_selected * 1 / ( 4 * pow( md2g_scaled, 4.0 ) ) * sqrt( 1 - pow( beta_abs, 2.0 ) ) / ( 1 - beta_abs ) * pow( coupling::get_maximum_coupling() , 3.0 );

  do
  {
    u = u_min + ran2() * ( u_max - u_min );
    phi = phi_min + ran2() * ( phi_max - phi_min );

    g = getMatrixElement( u, phi );

    if ( g > max )
    {
      //cout << "error in get32    g = " << g << "   max = " << max <<  endl;
      error++;
    }

    gr = max * ran2();

  }
  while ( g < gr );

  u_arg = u;
  phi_arg = phi;

  return error; // only necessary when using the rejection method, counts how many times the comparison function is below the actual function
}



/**
 * Routine for sampling u = cos(theta) and phi according to the matrix element (see scattering32::getMatrixElement(..)).
 * An implementation of the metropolis algorithm is used for sampling the distribution. This algorithm does not depend on
 * the knowledge of an absolute normalisation of the distribution
 *
 * Furthermore one needs no comparison function which makes this routine a lot faster than the version implemented in
 * scattering32::getMomenta32(..) as the comparison function in use there can be off the actual distribution by several orders
 * of magnitude, leading to high rejection probabilites for the sampled points.
 *
 * The computational expense is fixed by number of steps in the Markov chain (#n_steps).
 *
 * @param[out] u_arg u=cos(theta), theta = angle between p1' and p1
 * @param[out] phi_arg azimuthal angle between e_x and (p1' e_x)
 * @return By definition 0, provided for compability with #getMomenta32_rejection
 */
int scattering32::getMomenta32_metropolis( double& u_arg, double& phi_arg ) const
{
  double u, phi;
  double u_new, phi_new;
  double g = 0;

  // the ranges in which the variables u and phi need to be sampled
  const double u_min = 0.0;
  const double u_max = 1.0;
  const double phi_min = 0.0;
  const double phi_max = M_PI;

  // only needed when using gaussian distributions for proposing new steps in the markov chain
  // width
  //double sigma_u = (u_max - u_min) / 3.0;
  //double sigma_phi = (phi_max - phi_min) / 3.0;

  // randomly select initial values of u and phi, such that
  do
  {
    u = u_min + ran2() * ( u_max - u_min );
    phi = phi_min + ran2() * ( phi_max - phi_min );

    g = getMatrixElement( u, phi );
  }
  while ( FPT_COMP_E( g, 0.0 ) );


  // number of steps in the Markov chain
  const int n_steps = 50;

  // do n_steps steps
  // the current location is (u,phi)
  // one steps consists of
  //  - proposing a new point (u',phi') according to a certain proposal distribution
  //  - calculate the matrix element g(u',phi') at this new point
  //  - if g(u',phi') > g(u,phi) accept the new point (u',phi')
  //  - else accept the new point (u',phi') with a probability g(u',phi') / g(u, phi)
  //
  for ( int i = 0; i < n_steps; i++ )
  {
    do
    {
      //u_new = distributions::gaussian(u,sigma_u);    // propose new u using a gaussian with width sigma_u around the current value of u
      u_new = u_min + ran2() * ( u_max - u_min );      // propose new u using a uniform distribution over the entire range

    }
    while ( u_new < u_min || u_new > u_max );

    do
    {
      //phi_new = distributions::gaussian(phi,sigma_phi);   // propose new phi using a gaussian with width sigma_phi around the current value of phi
      phi_new = phi_min + ran2() * ( phi_max - phi_min );   // propose new phi using a uniform distribution over the entire range
    }
    while ( phi_new < phi_min || phi_new > phi_max );     // check that the new values are in range

    double g_new = getMatrixElement( u_new, phi_new );           // calculate the matrix element at the proposed point

    double ratio = g_new / g;                                    // ratio of g(u',phi') / g(u,phi)

    if ( FPT_COMP_GE( ratio, 1.0 ) || ran2() < ratio )    // accept if g(u',phi') > g(u,phi) or with probability g(u',phi') / g(u,phi)
    {
      u = u_new;
      phi = phi_new;
      g = g_new;
    }
  }

  u_arg = u;
  phi_arg = phi;

  return 0;  // only necessary when using the rejection method, always 0 when using Metropolis
}



/**
 * Calculates the Matrix elment for 3->2 scatterings at given u and phi. The theta function stemming from the LPM-cutoff is taken into
 * account.
 *
 * @param[in] u u=cos(theta), theta = angle between p1' and p1
 * @param[in] phi azimuthal angle between e_x and (p1' e_x)
 * @return Matrix element for 3->2
 */
double scattering32::getMatrixElement( const double u, const double phi ) const
{
  // E1_selected, E3_selected, cos_gamma and N are private members of scattering32 and have been set in scattering32::getIntegral32_vegas()
  double kt2;
  double us, v, constraint;
  double cos_delta, sin_delta_2;
  double cos_Theta;
  double g;

  //double beta_abs = 0;   // for testing consistency with old implementation (beta_abs = 0 needs then to be set in integrand32.cpp as well)
  double sin_gamma = sqrt( 1.0 - pow( cos_gamma, 2.0 ) );

  us = sqrt( 1.0 - pow( u, 2.0 ) );
  v = cos( phi );

  cos_delta = sin_gamma * us * v + cos_gamma * u;     // cos(delta) = sin(gamma)*sin(theta)*cos(phi) + cos(gamma)*cos(theta)
  sin_delta_2 = 1.0 - pow( cos_delta, 2.0 );          // sin(delta)^2 = 1 - cos(delta)^2
  
  if ( E3_selected * cos_delta + E1_selected * u < 0 )
  {
    return 0;
  }

  kt2 = pow( E3_selected, 2.0 ) * sin_delta_2;               // kt^2 = E3^2 * sin(delta)^2

  // Theta (capital!) is the angle between the boost-vector beta' (here: beta_vec) and the vector p1'
  if ( FPT_COMP_E( beta_abs, 0.0 ) )
  {
    cos_Theta = 0;
  }
  else
  {
    cos_Theta = 1 / beta_abs * ( rotated_beta.X() * us * v + rotated_beta.Y() * us * sqrt( 1 - pow( v, 2.0 ) ) + rotated_beta.Z() * u );
  }

  // the constraint for kt^2 depending on E3, Theta, delta, beta', sqrt(s) and lambda
  constraint = fudge_factor_lpm / lambda_scaled * E3_selected / sqrt( 1 - pow( beta_abs, 2.0 ) ) * ( 1 + ( beta_abs * cos_delta * cos_Theta ) );

  if ( kt2 > constraint )
  {
    double qt2 = pow( E1_selected, 2.0 ) * pow( us, 2.0 );            // qt^2 = E1^2 * sin(theta)^2
    //vector product qt*kt = -E1*E3*sin(theta)*( cos(gamma)*sin(theta) - sin(gamma)*cos(theta)*cos(phi) )
    double qtkt = -E1_selected * E3_selected * us * ( cos_gamma * us - sin_gamma * u * v );
    double qt_kt_2 = qt2 + kt2 - 2 * qtkt;
    double qtkt_kt2 = qtkt - kt2;
    
    // The term ( 1 - xbar )^2 where xbar = E3 / sqrt(s) * exp(|y|).
    // With the help of some analytic rearrangements this can be cast into the form used below,
    // that is easier (faster) to compute than using xbar = E3 / sqrt(s) * exp(arccosh(1/sin(delta))) directly.
    double factor_1minusx_squared = pow( 1 - E3_selected * ( 1 + sqrt( 1 - sin_delta_2 ) ) , 2 );
//     double factor_1minusx_squared = 1;
    
    double matrixElement22_factor, Pg;
    
    // if y < 0 switch from p3 = p1 - q and p4 = p2 + q - k (where p1 + p2 -> p3 + p4 + p5) to
    // p3 = p1 + q - k and p4 = p2 - q. This amounts to the replacements below.
    if ( cos_delta < 0 ) // corresponds to y < 0
    {
      double qt2_temp = qt2;
      qt2 = qt_kt_2;
      qt_kt_2 = qt2_temp;
      
      qtkt_kt2 = -qtkt;
    }
    
    // this computes mandelstam t as t = -q^2 (definition of the vectors as given in the BAMPS notes)
    double mandelstam_t;
    if ( cos_delta > 0 ) // corresponds to y > 0
    {
      mandelstam_t = E1_selected * ( u - 1 );   //  t = (p1' - p1)^2 (as four vectors)
    }
    else
    {
      double E2 = sqrt( pow( E1_selected, 2 ) + 2 * E1_selected * E3_selected * cos_gamma + pow( E3_selected, 2 ) );  
      mandelstam_t = E3_selected * cos_delta + E1_selected * u - E2; //  t = (p2' - p2)^2 (as four vectors) 
    }
    
    if ( matrixElement32_22qt )
    {
      const double md2_22 = md2g_scaled / coupling::get_constant_coupling() * coupling::get_coupling( - qt2 * pow( sqrtS , 2.0 ) ); //!! scale of running alpha_s?
      matrixElement22_factor = 1 / pow(( qt2 + md2_22 ), 2.0 );
    }
    else
    {
      const double md2_22 = md2g_scaled / coupling::get_constant_coupling() * coupling::get_coupling( mandelstam_t * pow( sqrtS , 2.0 ) ); //!! scale of running alpha_s?
      matrixElement22_factor = 1.0 / pow( ( mandelstam_t - md2_22 ) , 2.0 );
    }
    
    const double md2_ktqt = md2g_scaled / coupling::get_constant_coupling() * coupling::get_coupling( -qt_kt_2 * pow( sqrtS , 2.0 ) ); //!!  what scale?
    
    if ( md2_use_simplified_GB )
    {
      Pg = qt2 / ( kt2 * ( qt_kt_2 + md2_ktqt ) );  
    }
    else
    {    
      Pg = 1 / kt2 + qt_kt_2 / pow( qt_kt_2 + md2_ktqt, 2 ) + 2 * ( qtkt_kt2 ) / ( kt2 * ( qt_kt_2 + md2_ktqt ) );
    }
    
    double alpha_s_prefactor = coupling::get_coupling( mandelstam_t * pow( sqrtS , 2.0 ) ) * coupling::get_coupling( mandelstam_t * pow( sqrtS , 2.0 ) ) * coupling::get_coupling( kt2 * pow( sqrtS , 2.0 ) ); //!! at which scale should alpha_s be evaluated?
    
    g = alpha_s_prefactor * factor_1minusx_squared * matrixElement22_factor * Pg;  // the matrix element
  }
  else
  {
    g = 0;
  }

  return g;
}



/**
 * Sets new momenta for the outgoing particle according to u and phi. The results are written to the vectors P1_arg and P2_arg.
 *
 * @param[out] P1_arg[] Momentum vector of outgoing particle 1
 * @param[out] P2_arg[] Momentum vector of outgoing particle 2
 * @param[in] u u=cos(theta), theta = angle between p1' and p1
 * @param[in] phi azimuthal angle between e_x and (p1' e_x)
 */
void scattering32::setNewMomenta32( VectorEPxPyPz & P1_arg, VectorEPxPyPz & P2_arg, const double u, const double phi ) const
{
  VectorEPxPyPz P1_selected, P3_selected;

  switch ( absorpedGluon )
  {
  case 1:
    P1_selected = P2cm;
    P3_selected = P1cm;
    break;
  case 2:
    P1_selected = P1cm;
    P3_selected = P2cm;
    break;
  case 3:
    P1_selected = P1cm;
    P3_selected = P3cm;
    break;
  default:
    std::string errMsg = "Error in scattering32::scattering32(..). absorpedGluon did not resolve to valid value.";
    throw eScatt32_error( errMsg );
    break;
  }

  VectorEPxPyPz ex, ey, ez;

  // unit vector ez[] in z-direction
  // z-direction is given by p1cm[]
  ez = P1_selected * (1./P1_selected(0));

  // unit vector ex[] in x-direction
  // x-direction roughly points in the direction of p3cm[] (see notes)
  // construction is via: ex[] = p3cm[] - ( (p1cm[]*p3cm[]) / (p1cm*p1cm) ) * p1cm[]
  //
  double p1p3 = Dot3( P1_selected, P3_selected ) / pow( P1_selected(0), 2 );

  ex = P3_selected - P1_selected * p1p3;
  ex *= 1./sqrt( ex.vec2() );

  // unit vector ey[] in y-direction
  // ey[] is contructed via vector product ez[] x ex[]

  ey = Cross( ez, ex );


  // u = cos(theta)
  double us;
  if ( u > 1.0 )
    us = 0.0;
  else
    us = sqrt( 1.0 - pow( u, 2.0 ) );  // us = sin(theta)

  // set new momentum vectors according to:
  //         p1'cm[] = sqrt(s)/2 * ( sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)
  //         p3'cm[] = - p1'cm[]
  // (see notes).

  P1_selected = ( ex * ( us * cos( phi ) ) + ey * ( us * sin( phi ) ) + ez * ( u ) ) * ( sqrtS / 2.0 );
  P3_selected = -P1_selected;
  P1_selected(0) = P3_selected(0) = sqrtS / 2.0;


  // boost new momentum vectors back to original frame

  VectorEPxPyPz P1cellTMP,P2cellTMP;

  LL_CM.boostInv(P1_selected,P3_selected, P1cellTMP,P2cellTMP);
  LL_cell.boostInv(P1cellTMP,P2cellTMP, P1_arg,P2_arg);

}



/**
 * Rotates the boost vector #beta_vec to a reference frame given by P1_arg and P3_arg.
 *
 * @param[in] P1_arg
 * @param[in] P3_arg
 * @param[out] beta_new The rotated vector is returned as the argument beta_new[]
 */
void scattering32::rotate_beta( const VectorEPxPyPz & P1_arg, const VectorEPxPyPz & P3_arg, VectorEPxPyPz & beta_new)
{
  //------------------
  // unit vectors for the reference frame used in 3->2 routines are constructed
  //
  VectorEPxPyPz ex, ey, ez;

  // unit vector ez[] in z-direction
  // z-direction is given by p1_arg[]
  ez = P1_arg * (1.0 / P1_arg(0));


  // unit vector ex[] in x-direction
  // x-direction roughly points in the direction of p3_arg[] (see notes)
  // construction is via: ex[] = p3_arg[] - ( (p1_arg[]*p3_arg[]) / (p1_arg*p1_arg) ) * p1_arg[]
  double p1p3 = Dot3( P1_arg, P3_arg ) / pow( P1_arg(0), 2 );

  ex = P3_arg - P1_arg * p1p3;
  ex *= 1./sqrt( ex.vec2() );

  // unit vector ey[] in y-direction
  // ey[] is contructed via vector product ez[] x ex[]

  ey = Cross( ez, ex );
  //
  //------------------


  // the boost vector (given as beta_vec, a private member of scattering32) is transformed (rotated) to the new reference frame
  // given by ex, ey and ez
  // construction is via: a'[i] = V[i,j]*a[j] (sum over j) with V[i,j] = scalar_product(e'_[i],e_[j])
  // e'_[1] for example corresponds to the vector ex[], whereas e_[1] would be the normal unit vector (1,0,0)
  //

 VectorEPxPyPz beta_vec = LL_CM.beta();

  beta_new(0) = 1.0;
  beta_new(1) = beta_vec(1) * ex(1) + beta_vec(2) * ex(2) + beta_vec(3) * ex(3);
  beta_new(2) = beta_vec(1) * ey(1) + beta_vec(2) * ey(2) + beta_vec(3) * ey(3);
  beta_new(3) = beta_vec(1) * ez(1) + beta_vec(2) * ez(2) + beta_vec(3) * ez(3);

}
// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;
