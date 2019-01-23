//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/scattering23.cpp $
//$LastChangedDate: 2014-12-12 20:54:03 +0100 (Fr, 12. Dez 2014) $
//$LastChangedRevision: 2018 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

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
#include <math.h>
#include <algorithm>

#include "scattering23.h"
#include "coupling.h"
#include "lorentz.h"
#include "FPT_compare.h"
#include "random.h"
#include "integrand23.h"
#include "integrate.h"
#include "prefactors23.h"

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



// <<-------------------------------------------------
// some wrapper for old routines to make new routines for heavy quarks with running coupling and more choices in formationTimeTyp, matrix element, etc. compatible with the old ones, try not to use anymore
/**
 * **Old routine. Do not use any more. See constructor below!**
 * 
 * This basic constructor sets all pointers to NULL and only sets up the interpolation23 object.
 * It MUST be used with scattering23::setParameter!
 *
 * @param[in] theI23_massless_arg Pointer to an interpolation23 object for matrix element with only massless particles
 */
scattering23::scattering23(const interpolation23 * const theI23_massless_arg) :
  beta(), cos_theta(), theta(),
  P1(), P2(), 
  P1cell(), P2cell(), 
  P1cm(), P2cm(), 
  F1(), F2(),
  m1(), m2(),
  md2g_wo_as_scaled(),
  sqrtS(), 
  lambda_scaled(),
  numberOfGluons( 0 ), 
  randomizeInput( true ),
  randomizedConfiguration(),
  theI23_massless( theI23_massless_arg ), 
  theI23_charm_m1( NULL ), 
  theI23_charm_m2( NULL ), 
  theI23_bottom_m1( NULL ), 
  theI23_bottom_m2( NULL ),
  I23onlineIntegration(),
  matrixElement23_22qt(),
  md2_counter_term_in_I23(),
  fudge_factor_lpm(),
  K_factor(),
  kappa()
{
}

/**
 * **Old routine. Do not use any more. See constructor below!**
 * 
 * @param[in] theI23_massless_arg Pointer to an interpolation23 object for matrix element with only massless particles
 * @param[in] v vector of the collective velocity of the cell
 * @param[in] P1_arg momentum vector of particle 1
 * @param[in] P2_arg momentum vector of particle 2
 * @param[in] F1_arg flavor of ingoing particle 1
 * @param[in] F2_arg flavor of ingoing particle 2
 * @param[in] sqrtS_arg center of mass energy, sqrt(s)
 * @param[in] md2g_wo_as_scaled_arg squared debye mass, scaled by 1/s and divided by the coupling alpha_s 
 * @param[in] lambda_scaled_arg mean free path, scaled by sqrt(s)
 * @param[in] _alpha_s coupling constant (obsolete and not needed anymore, for constant coupling coupling::get_constant_coupling() is used)
 * @param[in] _Ng Number of gluons in cell (currently not used)
 **/
scattering23::scattering23( const interpolation23 * const theI23_massless_arg, 
                            const VectorXYZ & v, 
                            const VectorEPxPyPz & P1_arg, 
                            const VectorEPxPyPz & P2_arg, 
			    const FLAVOR_TYPE F1_arg, 
			    const FLAVOR_TYPE F2_arg,
			    const double sqrtS_arg, 
			    const double md2g_scaled_arg, 
			    const double lambda_scaled_arg, 
			    const double _alpha_s, 
			    const int _Ng ) :
  P1( P1_arg ), P2( P2_arg ),
  F1( F1_arg ), F2( F2_arg ), 
  m1( 0 ), m2( 0 ), 
  md2g_wo_as_scaled( md2g_scaled_arg / coupling::get_constant_coupling() ), 
  sqrtS( sqrtS_arg ),
  lambda_scaled( lambda_scaled_arg ), 
  numberOfGluons( _Ng ), 
  randomizeInput( true ),
  theI23_massless( theI23_massless_arg ), 
  theI23_charm_m1( NULL ), 
  theI23_charm_m2( NULL ), 
  theI23_bottom_m1( NULL ), 
  theI23_bottom_m2( NULL ),
  I23onlineIntegration( false ), 
  formationTimeTyp( "bamps_org" ), 
  matrixElement23( "GBimproved" ), 
  matrixElement23_22qt(),
  md2_counter_term_in_I23( true ), 
  fudge_factor_lpm( 1 ),
  K_factor ( 1 ),
  kappa ( 1 )
{
  if( F1 > 2 * ParticlePrototype::max_N_light_flavor || F2 > 2 * ParticlePrototype::max_N_light_flavor ) // heavy quark involved
  {
      std::string errMsg = "Error in scattering23. Old constructor which handles only light partons is called with a heavy quark.";
      throw eScatt23_error( errMsg );
  }

  randomizedConfiguration = 0;
  
  if ( randomizeInput )
  {
    if ( ran2() < 0.5 )
    {
      // see above
    }
    else
    {
      randomizedConfiguration = 1;
      
      P1 = P2_arg;
      P2 = P1_arg;
      F1 = F2_arg;
      F2 = F1_arg;
    }
  }
  
  // boost momenta to rest frame of the cell
  LL_cell.setBeta( v );
  LL_cell.boost(P1,P2, P1cell,P2cell);

  LL_CM.setBetaCM(P1cell,P2cell);
  LL_CM.boost(P1cell,P2cell, P1cm,P2cm);

  beta = LL_CM.betaVal();
  cos_theta = CosTheta(LL_CM.beta(), P1cm);

  // absolute value of cos(theta) is taken since the direction of the second boost (from CMS to frame where gluon is emitted transversally)
  // is given by the sign of tanh(y), see notes for details

  theta = acos( cos_theta );
}

/**
 * This method sets all necessary parameters for a given particle
 * pair. Previous values are deleted or overwritten. Using this method
 * an scattering32 object can be re-used for multiple particle pair,
 * thus reducing the need to constantly creating new objects. 
 *
 * Either this method or the constructor taking the same arguments
 * MUST be called prior to any other methods of the class! 
 *
 * @param[in] v vector of the collective velocity of the cell
 * @param[in] P1_arg momentum vector of particle 1
 * @param[in] P2_arg momentum vector of particle 2
 * @param[in] sqrtS_arg center of mass energy, sqrt(s)
 * @param[in] md2g_scaled_arg squared debye mass, scaled by s
 * @param[in] lambda_scaled_arg mean free path, scaled by sqrt(s)
 * @return The absolute value of the boost velocity between cell frame and CMS (#beta_vec).
 **/
double scattering23::setParameter( const VectorXYZ & v, 
                                   const VectorEPxPyPz & P1_arg, 
                                   const VectorEPxPyPz & P2_arg, 
				   const FLAVOR_TYPE F1_arg, 
				   const FLAVOR_TYPE F2_arg, 
				   const double sqrtS_arg,
				   const double md2g_scaled_arg, 
				   const double lambda_scaled_arg, 
				   const double _alpha_s, 
				   const int _Ng ) 
{
  if( F1_arg > 2 * ParticlePrototype::max_N_light_flavor || F2_arg > 2 * ParticlePrototype::max_N_light_flavor ) // heavy quark involved
  {
      std::string errMsg = "Error in scattering23. Old setParameter() function which handles only light partons is called with a heavy quark.";
      throw eScatt23_error( errMsg );
  }
  double mass = 0;
  double md2g_wo_as_scaled_arg = md2g_scaled_arg / coupling::get_constant_coupling();
  return setParameter( v, P1_arg, P2_arg, F1_arg, F2_arg, mass, mass, sqrtS_arg, md2g_wo_as_scaled_arg, lambda_scaled_arg,
		       1.0, 1.0, 1.0, 1.0, false, "bamps_org", "GBimproved", true, 1.0, _Ng );
};
// --------------------------------------------------->>

/**
 * This basic constructor sets all pointers to NULL and only sets up
 * the interpolation23 object. It MUST be used with
 * scattering23::setParameter! 
 *
 * @param[in] theI23_massless_arg Pointer to an interpolation23 object for matrix element with only massless particles
 * @param[in] theI23_charm_m1_arg Pointer to an interpolation23 object for matrix element with heavy charm quark as particle 1.
 * @param[in] theI23_charm_m2_arg Pointer to an interpolation23 object for matrix element with heavy charm quark as particle 2.
 * @param[in] theI23_bottom_m1_arg Pointer to an interpolation23 object for matrix element with heavy bottom quark as particle 1.
 * @param[in] theI23_bottom_m2_arg Pointer to an interpolation23 object for matrix element with heavy bottom quark as particle 2.
 **/
scattering23::scattering23( const interpolation23 * const theI23_massless_arg,
			    const interpolation23 * const theI23_charm_m1_arg, 
			    const interpolation23 * const theI23_charm_m2_arg, 
			    const interpolation23 * const theI23_bottom_m1_arg, 
			    const interpolation23 * const theI23_bottom_m2_arg ) :
  beta(), cos_theta(), theta(),
  P1( ), P2( ), 
  P1cell( ), P2cell( ), 
  P1cm( ), P2cm( ), 
  F1(), F2(),
  m1(), m2(),
  md2g_wo_as_scaled(),
  sqrtS(), 
  lambda_scaled(),
  numberOfGluons( 0 ), 
  randomizeInput( true ),
  randomizedConfiguration(),
  theI23_massless( theI23_massless_arg ), 
  theI23_charm_m1( theI23_charm_m1_arg ), 
  theI23_charm_m2( theI23_charm_m2_arg ), 
  theI23_bottom_m1( theI23_bottom_m1_arg ), 
  theI23_bottom_m2( theI23_bottom_m2_arg ),
  I23onlineIntegration(),
  matrixElement23_22qt(),
  md2_counter_term_in_I23(),
  fudge_factor_lpm(),
  K_factor(),
  kappa()
{
}


/**
 *
 * @param[in] theI23_massless_arg Pointer to an interpolation23 object for matrix element with only massless particles
 * @param[in] theI23_charm_m1_arg Pointer to an interpolation23 object for matrix element with heavy charm quark as particle 1.
 * @param[in] theI23_charm_m2_arg Pointer to an interpolation23 object for matrix element with heavy charm quark as particle 2.
 * @param[in] theI23_bottom_m1_arg Pointer to an interpolation23 object for matrix element with heavy bottom quark as particle 1.
 * @param[in] theI23_bottom_m2_arg Pointer to an interpolation23 object for matrix element with heavy bottom quark as particle 2.
 * @param[in] v vector of the collective velocity of the cell
 * @param[in] P1_arg momentum vector of particle 1
 * @param[in] P2_arg momentum vector of particle 2
 * @param[in] F1_arg flavor of ingoing particle 1
 * @param[in] F2_arg flavor of ingoing particle 1
 * @param[in] m1_arg mass of ingoing particle 1
 * @param[in] m2_arg mass of ingoing particle 1
 * @param[in] sqrtS_arg center of mass energy, sqrt(s)
 * @param[in] md2g_wo_as_scaled_arg squared debye mass, scaled by 1/s and divided by the coupling alpha_s 
 * @param[in] lambda_scaled_arg mean free path, scaled by sqrt(s)
 * @param[in] K_light_parton_arg K factor for process 2->3 with only light quarks
 * @param[in] K_heavy_quark_arg K factor for process 2->3 with one heavy quark
 * @param[in] kappa_light_parton_arg Kappa for Debye screening for process 2->3 with only light quarks
 * @param[in] K_heavy_quark_arg Kappa for Debye screening for process 2->3 with only light quarks
 * @param[in] I23onlineIntegration_arg  true if the integration of the total cross section for 2->3 processes is performed online and not read in from table
 * @param[in] formationTimeTyp_arg type of formation time of radiated gluon in 2->3
 * possibilities are (just give string as input, eg. "bamps_org" ):
 * - **bamps_org**: as implemented in BAMPS, formation time in frame in which gluon is 
 *   purely transverse (Sigma'') is tau'' = 1/kt, boost lambda to this frame and compare
 * - **bamps_org_extended**: bamps_org extended to massive sector, formation time in frame 
 *   in which gluon is purely transverse is tau'' = kt / ( kt^2 + x^2 M^2), with x = k+ / P_hq+
 * - **bamps_org_extended_xwE**: as bamps_org_extended, but with definition x = omega / E_hq
 * - **compare_cm_frame**: compare formation time in center of mass frame, boost lamda to cm frame, 
 *   here tau' = (2) omega / ( kt^2 + x^2 M^2)
 * - **compare_lab_frame**: compare formation time in lab frame, boost kt, omega, x to lab frame, 
 *   here tau = (2) omega / ( kt^2 + x^2 M^2)
 * .
 * all these formulas miss a factor 2 in the formation time since it is omitted in the original 
 * BAMPS version, but according to most literature it should be there... It can be added by 
 * setting fudge_factor to 2
 * @param[in] matrixElement23_arg matrix element which is used for 2->3 scattering
 * possibilities are (just give string as input, eg. "GBimproved" ):
 * - **GBimproved**: Gunion-Bertsch improved version which is correct also at forward and mid-rapidity
 * - **GBorg**: original Gunion-Bertsch version included in BAMPS (is not correct)
 * - **KPRapprox**: KPR matrix element in the limit k5->0, see notes and http://arxiv.org/abs/1109.5539
 * @param[in] md2_counter_term_in_I23_arg Whether a counter term is applied which resembles the old prescription of Xu for masssless particles
 * @param[in] fudge_factor_lpm_arg factor for LPM effect
 * this factor is used to play around with the cutoff
 * instead of k_t > gamma / lambda a modified cutoff k_t > X gamma / lambda
 * is used, where X is the fudge factor
 * @param[in] _Ng Number of gluons in cell (currently not used)
 **/
scattering23::scattering23( const interpolation23 * const theI23_massless_arg, 
			    const interpolation23 * const theI23_charm_m1_arg, 
			    const interpolation23 * const theI23_charm_m2_arg, 
			    const interpolation23 * const theI23_bottom_m1_arg,
			    const interpolation23 * const theI23_bottom_m2_arg, 
                            const VectorXYZ & v, 
                            const VectorEPxPyPz & P1_arg, 
                            const VectorEPxPyPz & P2_arg, 
			    const FLAVOR_TYPE F1_arg, 
			    const FLAVOR_TYPE F2_arg, 
			    const double m1_arg, 
			    const double m2_arg,
			    const double sqrtS_arg, 
			    const double md2g_wo_as_scaled_arg, 
			    const double lambda_scaled_arg,
			    const double K_light_parton_arg, 
			    const double K_heavy_quark_arg, 
			    const double kappa_light_parton_arg, 
			    const double kappa_heavy_quark_arg, 
			    const bool I23onlineIntegration_arg, 
			    const string & formationTimeTyp_arg, 
			    const string & matrixElement23_arg, 
			    const bool md2_counter_term_in_I23_arg,
			    const double fudge_factor_lpm_arg,  
			    const int _Ng, 
                            const bool matrixElement23_22qt_arg ) : 
  P1( P1_arg ), P2( P2_arg ),
  F1( F1_arg ), F2( F2_arg ), 
  m1( m1_arg ), m2( m2_arg ), 
  md2g_wo_as_scaled( md2g_wo_as_scaled_arg ), 
  sqrtS( sqrtS_arg ), 
  lambda_scaled( lambda_scaled_arg ), 
  numberOfGluons( _Ng ), 
  randomizeInput( true ),
  theI23_massless( theI23_massless_arg ), 
  theI23_charm_m1( theI23_charm_m1_arg ), 
  theI23_charm_m2( theI23_charm_m2_arg ), 
  theI23_bottom_m1( theI23_bottom_m1_arg ), 
  theI23_bottom_m2( theI23_bottom_m2_arg ),
  I23onlineIntegration( I23onlineIntegration_arg ), 
  formationTimeTyp( formationTimeTyp_arg ), 
  matrixElement23( matrixElement23_arg ), 
  matrixElement23_22qt( matrixElement23_22qt_arg ), 
  md2_counter_term_in_I23( md2_counter_term_in_I23_arg ), 
  fudge_factor_lpm( fudge_factor_lpm_arg )
{
  if( matrixElement23 != "GBimproved" && matrixElement23 != "GBorg" && matrixElement23 != "KPRapprox" )
  {
    std::string errMsg = "Error in scattering23 constructor. Value for matrixElement23 is not valid.";
    throw eScatt23_error( errMsg );
  }

  if( F1 > 2 * ParticlePrototype::max_N_light_flavor || F2 > 2 * ParticlePrototype::max_N_light_flavor ) // heavy quark involved
  {
    K_factor = K_heavy_quark_arg;
    kappa = kappa_heavy_quark_arg;
  }
  else
  {
    K_factor = K_light_parton_arg;
    kappa = kappa_light_parton_arg;
  }

  randomizedConfiguration = 0;
  
  if ( randomizeInput )
  {
    if ( ran2() < 0.5 )
    {
      // see above
    }
    else
    {
      randomizedConfiguration = 1;
      
      P1 = P2_arg;
      P2 = P1_arg;
      F1 = F2_arg;
      F2 = F1_arg;
      m1 = m2_arg;
      m2 = m1_arg;
    }
  }
  
  // boost momenta to rest frame of the cell
  LL_cell.setBeta( v );
  LL_cell.boost(P1,P2, P1cell,P2cell);

  LL_CM.setBetaCM(P1cell,P2cell);
  LL_CM.boost(P1cell,P2cell, P1cm,P2cm);

  beta = LL_CM.betaVal();
  cos_theta = CosTheta(LL_CM.beta(), P1cm);

  // absolute value of cos(theta) is taken since the direction of the
  // second boost (from CMS to frame where gluon is emitted
  // transversally) is given by the sign of tanh(y), see notes for
  // details 

  theta = acos( cos_theta );
}


scattering23::~scattering23()
{
}


/**
 * This method sets all necessary parameters for a given particle
 * pair. Previous values are deleted or overwritten. 
 * Using this method an scattering32 object can be re-used for multiple
 * particle pair, thus reducing the need to constantly creating new
 * objects. 
 *
 * Either this method or the constructor taking the same arguments MUST
 * be called prior to any other methods of the class! 
 *
 * @param[in] v vector of the collective velocity of the cell
 * @param[in] P1_arg momentum vector of particle 1
 * @param[in] P2_arg momentum vector of particle 2
 * @param[in] F1_arg flavor of ingoing particle 1
 * @param[in] F2_arg flavor of ingoing particle 1
 * @param[in] m1_arg mass of ingoing particle 1
 * @param[in] m2_arg mass of ingoing particle 1
 * @param[in] sqrtS_arg center of mass energy, sqrt(s)
 * @param[in] md2g_wo_as_scaled_arg squared debye mass, scaled by 1/s and divided by the coupling alpha_s 
 * @param[in] lambda_scaled_arg mean free path, scaled by sqrt(s)
 * @param[in] K_light_parton_arg K factor for process 2->3 with only light quarks
 * @param[in] K_heavy_quark_arg K factor for process 2->3 with one heavy quark
 * @param[in] kappa_light_parton_arg Kappa for Debye screening for process 2->3 with only light quarks
 * @param[in] K_heavy_quark_arg Kappa for Debye screening for process 2->3 with only light quarks
 * @param[in] I23onlineIntegration_arg  true if the integration of the total cross section for 2->3 processes is performed online and not read in from table
 * @param[in] formationTimeTyp_arg type of formation time of radiated gluon in 2->3
 * possibilities are (just give string as input, eg. "bamps_org" ):
 * - **bamps_org**: as implemented in BAMPS, formation time in frame in which gluon is 
 *   purely transverse (Sigma'') is tau'' = 1/kt, boost lambda to this frame and compare
 * - **bamps_org_extended**: bamps_org extended to massive sector, formation time in frame 
 *   in which gluon is purely transverse is tau'' = kt / ( kt^2 + x^2 M^2), with x = k+ / P_hq+
 * - **bamps_org_extended_xwE**: as bamps_org_extended, but with definition x = omega / E_hq
 * - **compare_cm_frame**: compare formation time in center of mass frame, boost lamda to cm frame, 
 *   here tau' = (2) omega / ( kt^2 + x^2 M^2)
 * - **compare_lab_frame**: compare formation time in lab frame, boost kt, omega, x to lab frame, 
 *   here tau = (2) omega / ( kt^2 + x^2 M^2)
 * .
 * all these formulas miss a factor 2 in the formation time since it is omitted in the original 
 * BAMPS version, but according to most literature it should be there... It can be added by 
 * setting fudge_factor to 2
 * @param[in] matrixElement23_arg matrix element which is used for 2->3 scattering
 * possibilities are (just give string as input, eg. "GBimproved" ):
 * - **GBimproved**: Gunion-Bertsch improved version which is correct also at forward and mid-rapidity
 * - **GBorg**: original Gunion-Bertsch version included in BAMPS (is not correct)
 * - **KPRapprox**: KPR matrix element in the limit k5->0, see notes and http://arxiv.org/abs/1109.5539
 * @param[in] md2_counter_term_in_I23_arg Whether a counter term is applied which resembles the old prescription of Xu for masssless particles
 * @param[in] fudge_factor_lpm_arg factor for LPM effect
 * this factor is used to play around with the cutoff
 * instead of k_t > gamma / lambda a modified cutoff k_t > X gamma / lambda
 * is used, where X is the fudge factor
 * @param[in] _Ng Number of gluons in cell (currently not used)
 * @return The absolute value of the boost velocity between cell frame and CMS (#beta_vec).
 **/
double scattering23::setParameter( const VectorXYZ & v, 
                                   const VectorEPxPyPz & P1_arg, 
                                   const VectorEPxPyPz & P2_arg,
                                   const FLAVOR_TYPE F1_arg, 
                                   const FLAVOR_TYPE F2_arg, 
                                   const double m1_arg, 
                                   const double m2_arg,
                                   const double sqrtS_arg, 
                                   const double md2g_wo_as_scaled_arg, 
                                   const double lambda_scaled_arg,
                                   const double K_light_parton_arg, 
                                   const double K_heavy_quark_arg, 
                                   const double kappa_light_parton_arg, 
                                   const double kappa_heavy_quark_arg, 
                                   const bool I23onlineIntegration_arg, 
                                   const string & formationTimeTyp_arg, 
                                   const string & matrixElement23_arg, 
                                   const bool md2_counter_term_in_I23_arg, 
                                   const double fudge_factor_lpm_arg,  
                                   const int _Ng, 
                                   const bool matrixElement23_22qt_arg )
{
  sqrtS = sqrtS_arg;
  md2g_wo_as_scaled = md2g_wo_as_scaled_arg;
  lambda_scaled = lambda_scaled_arg;
  I23onlineIntegration = I23onlineIntegration_arg;
  numberOfGluons = _Ng;
  formationTimeTyp = formationTimeTyp_arg;
  matrixElement23 = matrixElement23_arg;
  matrixElement23_22qt = matrixElement23_22qt_arg;
  md2_counter_term_in_I23 = md2_counter_term_in_I23_arg;
  fudge_factor_lpm = fudge_factor_lpm_arg; 
  
  if( matrixElement23 != "GBimproved" && matrixElement23 != "GBorg" && matrixElement23 != "KPRapprox" )
  {
    std::string errMsg = "Error in scattering23::setParameter(). Value for matrixElement23 is not valid.";
    throw eScatt23_error( errMsg );
  }

  randomizedConfiguration = 0;

  P1 = P1_arg;
  P2 = P2_arg;
  F1 = F1_arg;
  F2 = F2_arg;
  m1 = m1_arg;
  m2 = m2_arg;

  if ( randomizeInput )
  {
    if ( ran2() < 0.5 )
    {
      // see above
    }
    else
    {
      randomizedConfiguration = 1;
      
      P1 = P2_arg;
      P2 = P1_arg;
      F1 = F2_arg;
      F2 = F1_arg;
      m1 = m2_arg;
      m2 = m1_arg;
    }
  }

  if( F1 > 2 * ParticlePrototype::max_N_light_flavor || F2 > 2 * ParticlePrototype::max_N_light_flavor ) // heavy quark involved
  {
    K_factor = K_heavy_quark_arg;
    kappa = kappa_heavy_quark_arg;
  }
  else
  {
    K_factor = K_light_parton_arg;
    kappa = kappa_light_parton_arg;
  }

  // boost momenta to rest frame of the cell
  LL_cell.setBeta( v );
  LL_cell.boost(P1,P2, P1cell,P2cell);

  LL_CM.setBetaCM(P1cell,P2cell);
  LL_CM.boost(P1cell,P2cell, P1cm,P2cm);

  beta = LL_CM.betaVal();
  cos_theta = CosTheta(LL_CM.beta(), P1cm);

  // absolute value of cos(theta) is taken since the direction of the second boost (from CMS to frame where gluon is emitted transversally)
  // is given by the sign of tanh(y), see notes for details

  theta = acos( cos_theta );

  return beta;
}


/** @brief Returns total cross section in 1/GeV^2. */
double scattering23::getXSection23( int& initialStateIndex ) const 
{
  const double s = pow( sqrtS , 2.0 );
  const double m1_2 = pow( m1 , 2.0 );
  const double m2_2 = pow( m2 , 2.0 );
  const double prefactor = 1.0 / ( sqrt( pow( s - m1_2 - m2_2 , 2.0 ) - 4.0 * m1_2 * m2_2 ) );
  
  // check if s is too small and pQCD is not applicable
  if( s < 1.1 * ns_casc::lambda2 )
    return 0;
  else
    return ( prefactor * ns_casc::Ncolor * getIntegral23( initialStateIndex ) );
};


/**
 * @return Integral over the 2->3 matrix element, I23.
 * multiply it with 
 *  1 / s * Ncolor * scatt23_object.getIntegral23( initialStateIndex ); //1/GeV^2
 * to get the total cross section
 */
double scattering23::getIntegral23( int& initialStateIndex ) const
{
  double I23_gg_ggg = 0;

  // sort F1 and F2 such that comparisons below are easier
  unsigned int _F1 = std::min( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
  unsigned int _F2 = std::max( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );

  double total_prefactor = 0;
  initialStateIndex = -1;
  
  if ( _F2 > 10 ) // One of the scattering partners is not a parton. These cross sections here are only for parton scatterings.
  {
    initialStateIndex = 6;

    total_prefactor = 0;
  }
  else if (( _F1 + _F2 ) == 0 ) // g+g -> g+g+g, g+g -> q+qbar+g
  {
    prefactor23_gg_ggg pObj1;
    prefactor23_gg_qqbarg pObj2;
    
    initialStateIndex = 0;

    total_prefactor = pObj1.prefactor() * pObj1.symFactor23() + pObj2.prefactor() * pObj2.symFactor23();
  }
  else if ( _F1 >= 7 && _F2 >= 7 )  // heavy quark scattering with another heavy quark is not allowed (and also negligible)
  {
    initialStateIndex = 5;

    total_prefactor = 0;
  }
  else if ( _F1 == _F2 )  // q+q -> q+q+g, qbar+qbar -> qbar+qbar+g
  {
    prefactor23_qq_qqg pObj1;
    
    initialStateIndex = 1;

    total_prefactor = pObj1.prefactor() * pObj1.symFactor23();
  }
  else if (( _F1 * _F2 ) == 0 ) // g+q -> g+q+g, g+qbar -> g+qbar+g
  {
    prefactor23_qg_qgg pObj1;
    
    initialStateIndex = 2;

    total_prefactor = pObj1.prefactor() * pObj1.symFactor23();
  }
  else if (( _F2 - _F1 ) == 1 && ( _F2 % 2 ) == 0 )  // q+qbar -> q+qbar+g, q+qbar -> q'+qbar'+g, q+qbar -> g+g+g
  {
    prefactor23_qqbar_qqbarg pObj1;
    prefactor23_qqbar_qqbarDashg pObj2;
    prefactor23_qqbar_ggg pObj3;
    
    initialStateIndex = 3;

    total_prefactor = pObj1.prefactor() * pObj1.symFactor23() + pObj2.prefactor() * pObj2.symFactor23() + pObj3.prefactor() * pObj3.symFactor23();
  }
  else // q+q' -> q+q'+g, q+qbar' -> q+qbar'+g
  {
    prefactor23_qqdash_qqdashg pObj1;
    
    initialStateIndex = 4;

    total_prefactor = pObj1.prefactor() * pObj1.symFactor23();
  }
  
  total_prefactor = total_prefactor * K_factor;


  if ( total_prefactor > 0 )
  {
    bool I23_result_trustable = true;

    if ( !I23onlineIntegration )
    {
      if( m1 == 0.0 && m2 == 0.0 ) //massless
        I23_gg_ggg = theI23_massless->getI23( log( md2g_wo_as_scaled ), log( lambda_scaled ), beta, fabs( cos_theta ), log( sqrtS*sqrtS ), I23_result_trustable );
      else if( m1 > 0.0 && m2 == 0.0 ) // heavy quark is particle 1
      {
        if( m1 == ParticlePrototype::Mcharm )
          I23_gg_ggg = theI23_charm_m1->getI23( log( md2g_wo_as_scaled ), log( lambda_scaled ), beta, fabs( cos_theta ), log( sqrtS*sqrtS ), I23_result_trustable );
        else if ( m1 == ParticlePrototype::Mbottom )
          I23_gg_ggg = theI23_bottom_m1->getI23( log( md2g_wo_as_scaled ), log( lambda_scaled ), beta, fabs( cos_theta ), log( sqrtS*sqrtS ), I23_result_trustable );
        else
          cout << "Error in getIntegral23(): " << m1 << "  " << m2 << endl;
      }
      else if( m2 > 0.0 && m1 == 0.0 ) // heavy quark is particle 2
      {
        if( m2 == ParticlePrototype::Mcharm )
          I23_gg_ggg = theI23_charm_m2->getI23( log( md2g_wo_as_scaled ), log( lambda_scaled ), beta, fabs( cos_theta ), log( sqrtS*sqrtS ), I23_result_trustable );
        else if ( m2 == ParticlePrototype::Mbottom )
          I23_gg_ggg = theI23_bottom_m2->getI23( log( md2g_wo_as_scaled ), log( lambda_scaled ), beta, fabs( cos_theta ), log( sqrtS*sqrtS ), I23_result_trustable );
        else
          cout << "Error in getIntegral23(): " << m1 << "  " << m2 << endl;
      }
      else
        cout << "Error in getIntegral23() both particles massive: " << m1 << "  " << m2 << endl;
      
      
      if( std::isinf( I23_gg_ggg ) )
      {
        cout << "interpolated I23 is inf. " << I23_gg_ggg << "\t" << md2g_wo_as_scaled << "\t" << lambda_scaled << "\t" << beta << "\t" << cos_theta << "\t" << sqrtS*sqrtS << "\t" << I23_result_trustable << "\t" << F1 << "\t" << F2 <<  endl;
        I23_result_trustable = false;
      }
    }
    
    if( I23onlineIntegration || !I23_result_trustable )
    {
      integrand23 theIntegrand;
      theIntegrand.set_md2_wo_as( md2g_wo_as_scaled );
      theIntegrand.set_lambda( lambda_scaled );
      if( matrixElement23 == "GBorg" )
        theIntegrand.set_cos_theta( fabs(cos_theta) );
      else
        theIntegrand.set_cos_theta( cos_theta );
      theIntegrand.set_beta( beta );
      theIntegrand.set_m1( m1 / sqrtS );
      theIntegrand.set_m2( m2 / sqrtS );
      theIntegrand.set_s( sqrtS*sqrtS );
      theIntegrand.set_kappa( kappa );
      theIntegrand.set_formationTimeTyp( formationTimeTyp );
      theIntegrand.set_matrixElement23( matrixElement23 );
      theIntegrand.set_md2_counter_term_in_I23( md2_counter_term_in_I23 );
      theIntegrand.set_fudge_factor_lpm( fudge_factor_lpm );
      theIntegrand.set_matrixElement23_22qt( matrixElement23_22qt );

      //--------------------------------------------------------------
      // parameters needed for calls to the integration routines
      int neval;  // actual number of integrand evaluations needed
      int fail;   // 0 = desired accuracy was reached, 1 = accuracy not reached, -1 = dimension out of range
      double intResult[NCOMP_23]; // result of the integration over the unit hypercube, NCOMP_23 = #components, 1 in our case
      double error[NCOMP_23];     // presumed absolute error of integral
      double prob[NCOMP_23];      // xi^2 probability that error is NOT a reliable estimate of the true integration error
      //--------------------------------------------------------------

      //--------------------------------------------------------------
      // create the functionoid that handles the integration
      // will be called later with: integrate( theIntegrand, neval, fail, intResult, error, prob );
      // the integration routine is chosen via precompiler defined switches
      #ifdef CUBA_VEGAS
      integrate_vegas integrate( NDIM_23, NCOMP_23, EPSREL, EPSABS, VERBOSE, MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH );
      #endif

      #ifdef CUBA_SUAVE
      integrate_suave integrate( NDIM_23, NCOMP_23, EPSREL, EPSABS, VERBOSE, MINEVAL, MAXEVAL, NNEW, FLATNESS );
      #endif

      #ifdef CUBA_DIVONNE
      const int ngiven = 0;
      const int ldxgiven = 0;
      double * xgiven;
      int nregions = -1;
      integrate_divonne integrate( NDIM_23, NCOMP_23, EPSREL, EPSABS, VERBOSE, MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
                                   BORDER, MAXCHISQ, MINDEVIATION, ngiven, ldxgiven, xgiven, NEXTRA, PEAKFINDER, nregions );
      #endif

      #ifdef NR_VEGAS
      integrate_nr_vegas integrate( NDIM_23 );
      #endif
      //--------------------------------------------------------------

      integrate( theIntegrand, neval, fail, intResult, error, prob );

      I23_gg_ggg = ( 9.0 / ( 2.0 * M_PI ) * intResult[0] );
      
      if( I23_gg_ggg < 0.0 )
        return 0.0;
    }
    
    if( std::isinf( I23_gg_ggg ) )
    {
      cout << "computed I23 is inf. " << md2g_wo_as_scaled << "\t" << lambda_scaled << "\t" << beta << "\t" << cos_theta << "\t" << sqrtS*sqrtS << "\t" << I23_result_trustable << "\t" << F1 << "\t" << F2 <<  endl;
      return 0.0;
    }
  }

  return ( total_prefactor * I23_gg_ggg );
}



/**
 * Calculate the integral I23 of the process directly. This function is used to produce the table I23 for the interpolation.
 * 
 * @param[in] _md2_wo_as_scaled squared debye mass, scaled by 1/s and divided by the coupling alpha_s 
 * @param[in] _lambda_scaled mean free path, scaled by sqrt(s)
 * @param[in] _beta boost beta
 * @param[in] _cos_theta cosine of angle theta
 * @param[in] _m1_scaled mass of particle 1 scaled
 * @param[in] _m2_scaled mass of particle 2 scaled
 * @param[in] _I23onlineIntegration true if the integration of the total cross section for 2->3 processes is performed online and not read in from table
 * @param[in] _s Mandelstam s (needed for the running coupling since qt and kt are dimensionless without a scale)
 * @param[in] _kappa Kappa for Debye screening
 * @param[in] _formationTimeTyp type of formation time of radiated gluon in 2->3
 * possibilities are (just give string as input, eg. "bamps_org" ):
 * - **bamps_org**: as implemented in BAMPS, formation time in frame in which gluon is 
 *   purely transverse (Sigma'') is tau'' = 1/kt, boost lambda to this frame and compare
 * - **bamps_org_extended**: bamps_org extended to massive sector, formation time in frame 
 *   in which gluon is purely transverse is tau'' = kt / ( kt^2 + x^2 M^2), with x = k+ / P_hq+
 * - **bamps_org_extended_xwE**: as bamps_org_extended, but with definition x = omega / E_hq
 * - **compare_cm_frame**: compare formation time in center of mass frame, boost lamda to cm frame, 
 *   here tau' = (2) omega / ( kt^2 + x^2 M^2)
 * - **compare_lab_frame**: compare formation time in lab frame, boost kt, omega, x to lab frame, 
 *   here tau = (2) omega / ( kt^2 + x^2 M^2)
 * .
 * all these formulas miss a factor 2 in the formation time since it is omitted in the original 
 * BAMPS version, but according to most literature it should be there... It can be added by 
 * setting fudge_factor to 2
 * @param[in] _matrixElement23 matrix element which is used for 2->3 scattering
 * possibilities are (just give string as input, eg. "GBimproved" ):
 * - **GBimproved**: Gunion-Bertsch improved version which is correct also at forward and mid-rapidity
 * - **GBorg**: original Gunion-Bertsch version included in BAMPS (is not correct)
 * - **KPRapprox**: KPR matrix element in the limit k5->0, see notes and http://arxiv.org/abs/1109.5539
 * @param[in] _md2_counter_term_in_I23 Whether a counter term is applied which resembles the old prescription of Xu for masssless particles
 * @param[in] _fudge_factor_lpm factor for LPM effect
 * this factor is used to play around with the cutoff
 * instead of k_t > gamma / lambda a modified cutoff k_t > X gamma / lambda
 * is used, where X is the fudge factor
 * @param[out] I23_result_trustable Whether the interpolated value is trustable (currently not applied since interpolation cannot be done with this routine)
 * @return The absolute value of the boost velocity between cell frame and CMS (#beta_vec).
 * @return Integral over the 2->3 matrix element, I23.
 **/ 
double scattering23::getIntegral23_special( const double _md2_wo_as_scaled, 
                                            const double _lambda_scaled, 
                                            const double _beta, 
                                            const double _cos_theta, 
                                            const double _m1_scaled, 
                                            const double _m2_scaled, 
                                            const bool _I23onlineIntegration, 
                                            const double _s, 
                                            const double _kappa, 
                                            const string & _formationTimeTyp, 
                                            const string & _matrixElement23, 
                                            const bool _matrixElement23_22qt, 
                                            const bool _md2_counter_term_in_I23, 
                                            const double _fudge_factor_lpm, 
                                            bool & I23_result_trustable  )
{
  double I23_gg_ggg = 0.0;
  
  if ( !_I23onlineIntegration )
  {
//     if( _m1_scaled == 0.0 && _m2_scaled == 0.0 ) //massless
//       I23_gg_ggg = theI23_massless->getI23( log( _md2_wo_as_scaled ), log( _lambda_scaled ), _beta, _cos_theta, log( _s ), I23_result_trustable );
//     else if( _m1_scaled > 0.0 && _m2_scaled == 0.0 ) // heavy quark is particle 1
//     {
//       if( FPT_COMP_E( _m1_scaled, ParticlePrototype::Mcharm/sqrt(_s) ) )
//         I23_gg_ggg = theI23_charm_m1->getI23( log( _md2_wo_as_scaled ), log( _lambda_scaled ), _beta, _cos_theta, log( _s ), I23_result_trustable );
//       else if ( FPT_COMP_E( _m1_scaled, ParticlePrototype::Mbottom/sqrt(_s) ) )
//         I23_gg_ggg = theI23_bottom_m1->getI23( log( _md_wo_as2_scaled ), log( _lambda_scaled ), _beta, _cos_theta, log( _s ), I23_result_trustable );
//       else
//         cout << "Error in getIntegral23(): " << _m1_scaled << "  " << _m2_scaled << endl;
//     }
//     else if( _m2_scaled > 0.0 && _m1_scaled == 0.0 ) // heavy quark is particle 2
//     {
//       if( FPT_COMP_E( _m2_scaled, ParticlePrototype::Mcharm/sqrt(_s) ) )
//         I23_gg_ggg = theI23_charm_m2->getI23( log( _md2_wo_as_scaled ), log( _lambda_scaled ), _beta, _cos_theta, log( _s ), I23_result_trustable );
//       else if ( FPT_COMP_E( _m2_scaled, ParticlePrototype::Mbottom/sqrt(_s) ) )
//         I23_gg_ggg = theI23_bottom_m2->getI23( log( _md2_wo_as_scaled ), log( _lambda_scaled ), _beta, _cos_theta, log( _s ), I23_result_trustable );
//       else
//         cout << "Error in getIntegral23(): " << _m1_scaled << "  " << _m2_scaled << endl;
//     }
//     else
//     {
//       cout << "Error in getIntegral23() both particles massive: " << _m1_scaled << "  " << _m2_scaled << endl;
//       std::string errMsg = "Error in getIntegral23() both particles massive";
//       throw eScatt23_error( errMsg );
//     }
  }
  else
  {
    integrand23 theIntegrand;
    theIntegrand.set_md2_wo_as( _md2_wo_as_scaled );
    theIntegrand.set_lambda( _lambda_scaled );
    theIntegrand.set_cos_theta( _cos_theta );
    theIntegrand.set_beta( _beta );
    theIntegrand.set_m1( _m1_scaled );
    theIntegrand.set_m2( _m2_scaled );
    theIntegrand.set_s( _s );
    theIntegrand.set_kappa( _kappa );
    theIntegrand.set_formationTimeTyp( _formationTimeTyp );
    theIntegrand.set_matrixElement23( _matrixElement23 );
    theIntegrand.set_matrixElement23_22qt( _matrixElement23_22qt );
    theIntegrand.set_md2_counter_term_in_I23( _md2_counter_term_in_I23 );
    theIntegrand.set_fudge_factor_lpm( _fudge_factor_lpm );

    //--------------------------------------------------------------
    // parameters needed for calls to the integration routines
    int neval;  // actual number of integrand evaluations needed
    int fail;   // 0 = desired accuracy was reached, 1 = accuracy not reached, -1 = dimension out of range
    double intResult[NCOMP_23]; // result of the integration over the unit hypercube, NCOMP_23 = #components, 1 in our case
    double error[NCOMP_23];     // presumed absolute error of integral
    double prob[NCOMP_23];      // xi^2 probability that error is NOT a reliable estimate of the true integration error
    //--------------------------------------------------------------

    //--------------------------------------------------------------
    // create the functionoid that handles the integration
    // will be called later with: integrate( theIntegrand, neval, fail, intResult, error, prob );
    // the integration routine is chosen via precompiler defined switches
    #ifdef CUBA_VEGAS
    integrate_vegas integrate( NDIM_23, NCOMP_23, EPSREL, EPSABS, VERBOSE, MINEVAL, MAXEVAL, NSTART, NINCREASE );
    #endif

    #ifdef CUBA_SUAVE
    integrate_suave integrate( NDIM_23, NCOMP_23, EPSREL, EPSABS, VERBOSE, MINEVAL, MAXEVAL, NNEW, FLATNESS );
    #endif

    #ifdef CUBA_DIVONNE
    const int ngiven = 0;
    const int ldxgiven = 0;
    double * xgiven;
    int nregions = -1;
    integrate_divonne integrate( NDIM_23, NCOMP_23, EPSREL, EPSABS, VERBOSE, MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
                                 BORDER, MAXCHISQ, MINDEVIATION, ngiven, ldxgiven, xgiven, NEXTRA, PEAKFINDER, nregions );
    #endif

    #ifdef NR_VEGAS
    integrate_nr_vegas integrate( NDIM_23 );
    #endif
    //--------------------------------------------------------------

    integrate( theIntegrand, neval, fail, intResult, error, prob );

    I23_gg_ggg = ( 9.0 / ( 2.0 * M_PI ) * intResult[0] );
    
    if( I23_gg_ggg < 0.0 )
      return 0.0;
  }
  
  return I23_gg_ggg;
}


/**
 * Samples new momenta of the outgoing particles according to the matrix element using the rejection method.
 *
 * @param[out] pt1 Sampled qt (in GeV)
 * @param[out] pt3 Sampled kt (in GeV)
 * @param[out] y Sampled rapidity of emmited gluon.
 * @param[out] phi Sampled angle of the outgoing gluon.
 * @param[out] pz1 Sampled longitudinal momentum transfer.
 * @param[out] typ the type of collision, used for analysis purposes
 * @param[out] F1arg Sampled flavor of outgoing particle 1
 * @param[out] F2arg Sampled flavor of outgoing particle 2
 * @return Count how many time the comparison function was below the actual function. Only needed for rejection method. Should be zero.
 **/
int scattering23::getMomenta23( double& pt1, double& pt3, double& y, double& phi, double& pz1, int& typ, FLAVOR_TYPE& F1arg, FLAVOR_TYPE& F2arg ) const
{
  int error = 0;   //0 is returned if no error occurs, 1 otherwise

  double max;
  double qt2, qt, kt2, kt, g, gr;
  double y_min, y_max;
  double V_y, integrand_23;
  double M2,m1_2,m2_2;
  double pz11, pz12, dTf1, dTf2;
  double x_gluon_estimate;
  
  m1_2 = pow( m1/sqrtS ,2.0); // scaled by s
  m2_2 = pow( m2/sqrtS ,2.0); // scaled by s
  
  if( m1_2 >= m2_2 && m2_2 == 0.0 )
    M2 = m1_2;
  else if( m2_2 >= m1_2 && m1_2 == 0.0 )
    M2 = m2_2;
  else
  {
    M2 = 99.9; // Dummy value
    cout << "error in integrand23::operator(). Both masses are larger than 0: " << m1_2 << "  " << m2_2 << endl;
  }
  
  // To estimate comparison function and so forth calculated Debye mass with kappa and a const coupling since the scale of the running coupling cannot be fixed yet
  const double md2g_w_constAs_kappa_scaled = kappa * coupling::get_constant_coupling() * md2g_wo_as_scaled;

  const double lambda_scaled2 = pow( lambda_scaled , 2.0 );
  // check whether the lower boundary of kt integration is larger than upper boundary. Applies only for formationTimeTyp == "bamps_org", since the lower boundary is zero for other formationTimeTypes
  if ( pow( fudge_factor_lpm , 2.0) / lambda_scaled2 >= 0.25 && ( formationTimeTyp == "bamps_org" || M2 == 0.0 ) )
  {
    cout << "lambda_scaled2 = " << lambda_scaled2 << "  sqrtS = " << sqrtS  << " lambda = " << lambda_scaled/sqrtS*0.197 << "fm" << endl;
    cout << "A=" << beta * cos_theta << "   B_lq=" << lambda_scaled * sqrt(1.0-pow(beta,2.0)) / fudge_factor_lpm  << endl;
    std::string errMsg = "1 / lambda_scaled2 >= 0.25";
    throw eScatt23_error( errMsg );
  }
  
  // get maximum y range
  double ymin_kt_independent = 0, ymax_kt_independent = 0;
  maximal_Y_range_23( ymin_kt_independent, ymax_kt_independent );
  double maxV_y = ymax_kt_independent - ymin_kt_independent;
  
  // estimate a maximum constant for the comparison function
  const double deltaTrans_estimate = 0.001;
  double max0 = maxV_y / ( md2g_w_constAs_kappa_scaled * deltaTrans_estimate ) * pow( coupling::get_constant_coupling() , 3.0 );
  
  // estimate x_gluon
  if( formationTimeTyp == "bamps_org_extended" )
    x_gluon_estimate = 0.5;
  else if( formationTimeTyp == "bamps_org_extended_xwE" )
    x_gluon_estimate = 0.0001;
  else
    x_gluon_estimate = 0.0;
  const double x_gluon_2_estimate = pow( x_gluon_estimate , 2.0 );
  
  // calcalute kt and qt ranges needed in their sampling
  double maxV_kt2_x2M2_formTimeBampsOrg = log( ( pow( 1.0 - M2 , 2.0 ) /  4.0  + x_gluon_2_estimate * M2 ) / ( pow( fudge_factor_lpm , 2.0) / lambda_scaled2 + x_gluon_2_estimate * M2 )  );
  double maxV_kt2_x2M2_formTimeHQ = log( pow( 1.0 - M2 , 2.0 ) / ( 4.0 * x_gluon_2_estimate * M2 ) + 1.0 );
  double maxV_qt2 = log( pow( 1.0 - M2 , 2.0 ) / ( 4.0 * md2g_w_constAs_kappa_scaled ) + 1.0 );

  if ( randomizedConfiguration == 0 )
  {
    F1arg = F1;
    F2arg = F2;
  }
  else
  {
    F1arg = F2;
    F2arg = F1;   
  }
  
  // sort F1 and F2 such that comparisons below are easier
  // convert to unsigned int to allow for some math that speeds up the decisions
  unsigned int _F1 = std::min( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );
  unsigned int _F2 = std::max( static_cast<unsigned int>( F1 ), static_cast<unsigned int>( F2 ) );


  //TODO: add extra types for heavy quarks !!
  if (( _F1 + _F2 ) == 0 ) // g+g -> g+g+g, g+g -> q+qbar+g
  {
    prefactor23_gg_ggg pObj1;
    prefactor23_gg_qqbarg pObj2;
    double pf1 = pObj1.prefactor() * pObj1.symFactor23();
    double pf2 = pObj2.prefactor() * pObj2.symFactor23();
    
    if ( ran2() < pf1 / (pf1 + pf2) )
    {
      typ = 231; // g+g -> g+g+g
    }
    else
    {
      typ = 232; // g+g -> q+qbar+g
      sampleFlavor23( F1arg, F2arg );
    }
  }
  else if ( _F1 == _F2 )  // q+q -> q+q+g, qbar+qbar -> qbar+qbar+g
  {
    typ = 237;
  }
  else if ( (_F1 * _F2) == 0 ) // g+q -> g+q+g, g+qbar -> g+qbar+g
  {
    typ = 233;
  }
  else if (( _F2 - _F1 ) == 1 && ( _F2 % 2 ) == 0 )  // q+qbar -> q+qbar+g, q+qbar -> q'+qbar'+g, q+qbar -> g+g+g
  {
    prefactor23_qqbar_qqbarg pObj1;
    prefactor23_qqbar_qqbarDashg pObj2;
    prefactor23_qqbar_ggg pObj3;

    double pf1 = pObj1.prefactor() * pObj1.symFactor23();
    double pf2 = pObj2.prefactor() * pObj2.symFactor23();
    double pf3 = pObj3.prefactor() * pObj3.symFactor23();
    
    double select = ran2() * ( pf1 + pf2 + pf3 );
    
    if ( select < pf1 )  // q+qbar -> q+qbar+g
    {
      typ = 234;
    }
    else if ( select < ( pf1 + pf2 ) )  // q+qbar -> q'+qbar'+g
    {
      typ = 235;
      sampleFlavor23( F1arg, F2arg, F1 );  // sample flavor excluding F1 (or anti-F1 = F2)
    }
    else  // q+qbar -> g+g+g
    {
      typ = 236;
      F1arg = gluon;
      F2arg = gluon;
    }
  }
  else // q+q'-> q+q'+g, q+qbar' -> q+qbar'+g
  {
    typ = 238; // qq' -> qq', qqbar' -> qqbar'
  }
  

  int while_counter = 0;
  int while_count_max = 6000000;
  const int while_count_break = 10000000;
  do
  {
    if( formationTimeTyp == "bamps_org" || M2 == 0.0 )
      kt2 = - x_gluon_2_estimate * M2 + ( pow( fudge_factor_lpm , 2.0) / lambda_scaled2 + x_gluon_2_estimate * M2 ) * exp( ran2() * maxV_kt2_x2M2_formTimeBampsOrg );
    else if( ( formationTimeTyp == "bamps_org_extended" || formationTimeTyp == "bamps_org_extended_xwE" ) && M2 != 0.0 )
      kt2 = x_gluon_2_estimate * M2 * ( exp( ran2() * maxV_kt2_x2M2_formTimeHQ ) - 1.0 );
    else
      cout << "error in scattering23::getMomenta23()" << endl;
    
    qt2 = md2g_w_constAs_kappa_scaled * ( exp( ran2() * maxV_qt2 ) - 1.0 );

    qt = sqrt( qt2 );
    kt = sqrt( kt2 );
    
    bool solvable = y_range_23( y_min, y_max, kt, M2, m1_2, m2_2, lambda_scaled, beta, cos_theta, formationTimeTyp, fudge_factor_lpm );

    if ( solvable )
    {
      V_y = y_max - y_min;
      y = y_min + ( y_max - y_min ) * ran2();
      
//       if( y_min < ymin_kt_independent )
//       {
//         cout << "error in getMomenta23: ymin too small, ymin=" << y_min << "  ymin_ktin=" << ymin_kt_independent << endl;
//         cout << "sqrtS = " << sqrtS  << "   A=" << beta * cos_theta << "   B_lq=" << lambda_scaled * sqrt(1.0-pow(beta,2.0)) / fudge_factor << "   kt=" << kt  << endl;
//       }
//       if( y_max > ymax_kt_independent )
//       {
//         cout << "error in getMomenta23: ymax too large, ymax=" << y_max << "  ymax_ktin=" << ymax_kt_independent << ", ymin=" << y_min << "  ymin_ktin=" << ymin_kt_independent << endl;
//         cout << "sqrtS = " << sqrtS  << "   A=" << beta * cos_theta << "   B_lq=" << lambda_scaled * sqrt(1.0-pow(beta,2.0)) / fudge_factor << "   kt=" << kt  << endl;
//       }
    }
    else
    {
      V_y = 0;
      y = 0;
    }

    phi = M_PI * ran2();
    
    // Now that q_t, k_t, y and phi have been fixed, get
    // the argument of the sigma_23 integral at these values.
    integrand_23 = get_integrand_23( pz11, pz12, dTf1, dTf2, qt, kt, y, phi, M2, m1_2, m2_2, md2g_wo_as_scaled, lambda_scaled, beta, cos_theta, sqrtS*sqrtS, kappa, formationTimeTyp, matrixElement23, matrixElement23_22qt, md2_counter_term_in_I23, fudge_factor_lpm );

    g = V_y * integrand_23;
    
    max = max0 / ( qt2 + md2g_w_constAs_kappa_scaled ) / ( kt2 + x_gluon_2_estimate * M2 );
    
    if ( g > max )
    {
//       double x_gluon = kt * exp(y);
//       cout << "error in getMomenta23: result larger than comparison function: " << g << " > " << max << endl;
//       cout << "x=" << x_gluon << "  x_est=" << x_gluon_2_estimate << "   V_y=" << V_y << "   V_y_max=" << maxV_y << "   dT=" << dTf1+dTf2 << "   dT_est=" << 1.0/deltaTrans_estimate << endl;
      ++error;
    }

    gr = max * ran2();
    
    while_counter++;
    if(while_counter >= while_count_max)
    {
      cout << "=====================================" << endl;
      cout << while_counter << "\t" << g << "\t" << max << "\t" << gr << endl;
      cout << max << " = " << max0 << " * " << 1.0/( qt2 + md2g_w_constAs_kappa_scaled ) << " * " << 1.0/( kt2 + x_gluon_2_estimate * M2 ) << endl;
      cout << qt2 << "\t" << kt2 << endl;
      cout << md2g_wo_as_scaled << "\t" << lambda_scaled << "\t" << cos_theta << "\t" << beta << "\t" <<   beta * cos_theta << endl;
      int dummy;
      cout << getIntegral23( dummy ) << endl;
      while_count_max = while_count_max * 2;
    }
    
    if( while_counter >= while_count_break) // something is wrong -> break
    {
      cout << "======================================" << endl;
      cout << "== break while loop in getMomenta23 ==" << endl;
      cout << "======================================" << endl;
      
      qt = 0.0;
      kt = 0.0001;
      y = 0.0;
      phi = 0.0;
      
      dTf1 = 1.0;
      dTf2 = 0.0;
      pz11 = sqrt( pow( 1.0 - m1_2 - m2_2 , 2.0 ) - 4.0 * m1_2 * m2_2 ) / 2.0 ;
      pz12 = 0.0;
      
      break;
    }

  }
  while ( g < gr );
  

  if ( ran2() < dTf1 / ( dTf1 + dTf2 ) )
    pz1 = pz11 * sqrtS;
  else
    pz1 = pz12 * sqrtS;

  pt1 = qt * sqrtS;
  pt3 = kt * sqrtS;

  return error;
}



/**
 * Integration limits for y are either determined from the kinematic
 * constraint \f$ \cosh(y) < (s-M^2) / (\sqrt{s} 2 k_t) \f$
 * or from the requirement that \f$ \cosh(y)  + A \sinh(y)  < B \f$
 * stemming from the modelling of the LPM effect, where
 * \f$ A = \beta' \cos(\theta) \f$
 * and
 * \f$ B = k_t \lambda \sqrt{1-\beta'^2} \f$.
 *
 * It is: \f$ B > 0 \f$ and \f$0 \leq A \leq 1\f$. 
 * For the physical meaning of beta' and theta see above.
 *
 * @param[out] y_min Returns the minmial y value.
 * @param[out] y_max Returns the maximal y value.
 * @param[in] kt Momentum transfer kt scaled by s
 * @param[in] M2 mass of the heavy quark squared scaled by s
 * @param[in] lambda_scaled_yrange Mean free path scaled by s
 * @param[in] beta_yrange boost beta
 * @param[in] cos_theta_yrange cosine of angle theta
 * @param[in] formationTimeTyp_arg type of formation time of radiated gluon in 2->3
 * @param[in] fudge_factor_lpm_arg factor for LPM effect
 * @return Whether a valid range could be successfully determined.
 **/
bool scattering23::y_range_23( double& y_min, 
                               double& y_max, 
                               const double kt, 
                               const double M2, 
                               const double m1_2, 
                               const double m2_2, 
                               const double lambda_scaled_yrange, 
                               const double beta_yrange, 
                               const double cos_theta_yrange, 
                               const string & formationTimeTyp_yrange, 
                               const double fudge_factor_lpm_yrange )
{
  double A, B;
  bool success = true;
  double y_min_kin, y_max_kin, y_min_lpm, y_max_lpm;
  
  const double y_infinity = 10.0; // large value for rapidity, basically infinity
  
  //<<-----------------------------------------------
  // check for the constraints on y from kinematics 
  
  const double kin_y_constraint = ( 1.0 - M2 ) / ( 2.0 * kt );
  
  if ( FPT_COMP_LE(kin_y_constraint,1.0) )  // no solution in this case
  {
    success = false;
    return success;
  }
  else
  {
    y_max_kin = acosh(kin_y_constraint);
    y_min_kin = -acosh(kin_y_constraint);
  }
  //----------------------------------------------->>
  
  
  //<<-----------------------------------------------
  // check for the constraints on y from LPM effect
  if( ( formationTimeTyp_yrange == "bamps_org"  || M2 == 0.0 ) && fudge_factor_lpm_yrange != 0.0 )
  {
    A = beta_yrange * cos_theta_yrange;
    B = kt * lambda_scaled_yrange * sqrt(1.0-pow(beta_yrange,2.0)) / fudge_factor_lpm_yrange;

    // constrains from LPM effect
    // as B > 0 and -1 <= A <= 1 this should not evaluate to "true"
    if ( FPT_COMP_L(A,-1.0) || FPT_COMP_G(A,1.0) || FPT_COMP_L(B,0.0) )
    {
      cout << "encountered unexpected values when determining limits for y integration" << endl;
      cout << "A = " << A << endl;
      cout << "B = " << B << endl;

      success = false;
      return success;
    }
    else
    {  
      //<<---------------------------------------------
      // some special cases...
      if ( FPT_COMP_E(B,0.0) )                 // B = 0
      {
        success = false;
        return success;
      }
    
      if ( FPT_COMP_E(A,0.0) )                 // A = 0
      {
        if ( FPT_COMP_LE(B,1.0) )            // A = 0  and  B <= 1
        {
          success = false;
          return success;
        }
        else                                 // A = 0  and  B > 1
        {
          y_min_lpm = -acosh(B);
          y_max_lpm = acosh(B);
        }
      }
      //----------------------------------------------->>
      //<<-----------------------------------------------
      // the general case...
      else
      {
        if ( FPT_COMP_LE(B,(sqrt(1+A)*sqrt(1-A))) ) //B <= (sqrt(1+A)*sqrt(1-A)) -> no solution 
        {
          success = false;
          return success;
        }
        else
        {
          // constraints from LPM effect
          y_min_lpm = log( (B - sqrt(pow(A,2.0)+pow(B,2.0)-1))/(1+A) );
          y_max_lpm = log( (B + sqrt(pow(A,2.0)+pow(B,2.0)-1))/(1+A) );
        }
      }
      //----------------------------------------------->>
    }
  }
  else if( formationTimeTyp_yrange == "bamps_org_extended" && fudge_factor_lpm_yrange != 0.0 )
  {
    double a,b,c,d;
    int nmbSol;
    double sol[3];
    double u_min, u_max;
    bool hq_particle1;
    
    A = beta_yrange * cos_theta_yrange;
    B = kt * lambda_scaled_yrange * sqrt(1.0-pow(beta_yrange,2.0)) / fudge_factor_lpm_yrange;
    
    y_max_lpm = y_infinity;
    y_min_lpm = -y_infinity;
    
    if( m1_2 > 0 )
      hq_particle1 = true;
    else
      hq_particle1 = false;
    
    if( hq_particle1 )
    {
      a = - 2.0 * B * M2;
      b = 1.0 + A;
      c = - 2.0 * B;
      d = 1.0 - A;
    }
    else
    {
      a = 1.0 + A;
      b = - 2.0 * B;
      c = 1.0 - A;
      d = - 2.0 * B * M2;
    }
    
//     if ( FPT_COMP_E(a,0.0) ) // a = 0, since M2 != 0, B must be zero -> no solution
    if ( fabs( a ) < 1E-4 ) // a = 0, since M2 != 0, B must be zero -> no solution
    {
      success = false;
      return success;
    }
    
    nmbSol = rootsCubicPolynomial( a, b, c, d, sol );
    
    if ( nmbSol == 0 ) // if no solution there was a problem in rootsCubicPolynomial()     
    {
      cout << "problem in scattering23::y_range_23() - no solution in rootsCubicPolynomial" << endl;
      
      cout << a << "  " << b << "  " << c << "  " << d << "  " 
      << A << "  " << B << "  " << M2 << endl;
      
//      std::string errMsg = "error in y_range_23()";
//      throw eScatt23_error( errMsg );      
      success = false;
      return success;
    }
    else if ( nmbSol == 1 && sol[0] < 1E-20 ) // due to numerical errors the solution can get zero or even negative if alpha in rootsCubicPolynomial() gets very small. Set the solution to a small number to avoid these errors
    {
      double p = b / a;
      double q = c / a;
      double r = d / a;
      
      // helpful abbreviations
      double alpha=(1.0/3.0)*(3.0*q - p*p);
      double beta=(1.0/27.0)*(2.0*p*p*p - 9.0*p*q + 27.0*r);
      double root=(beta*beta/4.0)+(alpha*alpha*alpha/27.0);

      double a1=(-beta/2.0)+sqrt(root);
      double b1=(-beta/2.0)-sqrt(root);
//       if( alpha < 0.1 )
//         sol[0] = 1E-20;
//       else
//         cout << "1 solution in y_range_23() which is very small: " << sol[0] << "  alpha = " << alpha << endl;
      cout << "1 solution in y_range_23() which is very small: " << sol[0] << "  " << a << "  " << b << "  " << c << "  " << d << "  " 
      << p << "  " << q << "  " << r << "  " << alpha << "  " << beta << "  " << root << "  " << a1 << "  " << b1 << endl;
      sol[0] = 1E-20;
    }

    
//     cout << nmbSol << "  " <<  sol[0] << "  " << sol[1]  << "  " << sol[2] << endl;

    // umin is the smallest solution
    u_min = sol[0];
    for(int i=1; i < nmbSol; i++ )
      if( sol[i] < u_min )
        u_min = sol[i];
      
    // umax is the largest solution
    u_max = sol[0];
    for(int i=1; i < nmbSol; i++ )
      if( sol[i] > u_max )
        u_max = sol[i];
      

// cout << "test " << u_min << "  " << nmbSol << "  " << sol[0] << "  " << sol[1] << "  " << sol[2] << endl;
//       cout <<  A << "  " << B << "  " << beta_yrange << "  " << cos_theta_yrange << "  " << lambda_scaled_yrange << "  " << kt << "  " << M2 << "  " << a << "  " << b << "  " << c << "  " << d  << endl;
      
    if( u_min <= 0.0 || u_max <= 0.0 )
    {
      cout << "error in y_range_23(), formationTimeTyp == bamps_org_extended, " << u_min << "  " << nmbSol << "  " << sol[0] << "  " << sol[1] << "  " << sol[2] << endl;
      cout <<  A << "  " << B << "  " << beta_yrange << "  " << cos_theta_yrange << "  " << lambda_scaled_yrange << "  " << kt << "  " << M2 << "  " << a << "  " << b << "  " << c << "  " << d  << endl;
      std::string errMsg = "error in y_range_23()";
      throw eScatt23_error( errMsg );
//       success = false;
//       return success;
    }
    else
    {
      if( hq_particle1 )
        y_min_lpm = log( u_min );
      else
        y_max_lpm = log( u_max );
    }
    
    // special case of 3 solutions. Check whether largest (or smallest for hq_particle1==false) solution is larger (smaller) than y_max_kin (y_min_kin). y_max_lpm (y_min_lpm) must be the intermediate solution
    if( nmbSol == 3 )
    {
      double u_intermediate;

      int iSol = 0;
      do
      {
        if( iSol >= nmbSol )
        {
          std::string errMsg = "error in y_range_23(): determining intermediate solution";
          throw eScatt23_error( errMsg );
        }
        u_intermediate = sol[iSol];
        iSol++;
      }
      while( u_intermediate == u_min || u_intermediate == u_max );
      
      if( hq_particle1 )
      {
        if( log( u_max ) > y_max_kin )
          y_max_lpm = log( u_intermediate );
      }
      else
      {
        if( log( u_min ) < y_min_kin )
          y_min_lpm = log( u_intermediate );
      }
    }
  }  
  else if( formationTimeTyp_yrange == "bamps_org_extended_xwE" )
  {
    // no analytical solution obtained yet for the LPM constrain. So constrains due to LPM effect are checked in getIntegrand23() which causes more numerical effort ( the integration boundaries for y are too large especially at small kt, and the integrand is mostly 0. Though it gives the right results much more calls for Vegas are needed to obtain the same statistical significants as if the y limits could be set analytically ).
    
    y_min_lpm = -y_infinity;
    y_max_lpm = y_infinity;
  }
  else
    cout << "error in scattering23::y_range_23()" << endl;
  //----------------------------------------------->>

  //<<-----------------------------------------------
  // apply the stricter of both constrains
  y_max = min( y_max_kin, y_max_lpm );
  y_min = max( y_min_kin, y_min_lpm );
  
  // for numerics also check that the values are not exorbitantly large
  y_max = min( y_max, y_infinity );
  y_min = max( y_min, -y_infinity );

  
  if ( y_min > y_max )
  {
    success = false;
    return success;
  }
  //----------------------------------------------->>

  return success;
}

/**
* This function calculates the real roots of the equation
*    a x^3 + b x^2 + c x + d = 0
*
* based on formula of Cardano, see also http://www.mathe.tu-freiberg.de/~hebisch/seminar1/kubik.html
*
* @param[in] a coefficient of x^3
* @param[in] b coefficient of x^2
* @param[in] c coefficient of x
* @param[in] d constant coefficient
* @param[out] sol array with real solutions
* @return number of real solutions (maximum 3, minimum 1)
*/
int rootsCubicPolynomial( const double a, const double b, const double c, const double d, double sol[3] )
{
  double p,q,r;
  double root;
  double alpha,beta;
  double A,B;
  double x1,x2,x3;
  double a1,b1;
  double angle;
  
  if( a == 0.0 )
  {
    cout << "error in rootsCubicPolynomial()" << endl;
    return 0;
  }

  // transform to form y^3 + p*y^2 + q*y + r
  p = b / a;
  q = c / a;
  r = d / a;
  
  // helpful abbreviations
  alpha=(1.0/3.0)*(3.0*q - p*p);
  beta=(1.0/27.0)*(2.0*p*p*p - 9.0*p*q + 27.0*r);
  root=(beta*beta/4.0)+(alpha*alpha*alpha/27.0);
  
  const double largeNumber = 1E25;
  if( (beta*beta/4.0) > largeNumber && fabs(alpha*alpha*alpha/27.0) > largeNumber && root <= 0.0 ) // subtract 2 exorbitantly large numbers leads to numberical errors
  {
    if( a > 0.01 ) // this should only happen for small a
      cout << "error in rootsCubicPolynomial()  " << a << "  " << b << "  " << c << "  " << d  << "  " << (beta*beta/4.0) << "  " << fabs(alpha*alpha*alpha/27.0)   << "  " << (beta*beta/4.0) - fabs(alpha*alpha*alpha/27.0)  << endl;

    return 0;
  }

  a1=(-beta/2.0)+sqrt(root);
  b1=(-beta/2.0)-sqrt(root);
  
  // get cubic root, fabs(..) since no negative number in root allowed, if statement to get the sign right
  A = pow( fabs(a1) ,(1.0/3.0));
  if( a1 < 0.0 )
    A = -A;
  B = pow( fabs(b1) ,(1.0/3.0));
  if( b1 < 0.0 )
    B = -B;

  if( root > 0.0 ) // Only one real solution (and two complex)
  {
    x1 = A+B;
    sol[0] = x1 - (p/3.0);
    if( std::isnan(sol[0]) || std::isinf(sol[0]) )
      cout << "error in rootsCubicPolynomial() - 1 solution: " << sol[0] << "  " << p << "  " << q << "  " << r << "  " << a << "  " << b << "  " << c << "  " << d << endl;
    return 1;
  }
  else if( root < 0.0 ) // 3 real solutions
  {
    angle=acos((-beta/2.0)/sqrt(-alpha*alpha*alpha/27.0));
    x1=2.0*sqrt(-alpha/3.0)*cos(angle/3.0);
    x2=2.0*sqrt(-alpha/3.0)*cos((angle/3.0) + 2.0*M_PI/3.0);
    x3=2.0*sqrt(-alpha/3.0)*cos((angle/3.0) + 4.0*M_PI/3.0);
    
    sol[0] = x1 - (p/3.0);
    sol[1] = x2 - (p/3.0);
    sol[2] = x3 - (p/3.0);
    if( std::isnan(sol[0]) || std::isinf(sol[0]) || std::isnan(sol[1]) || std::isinf(sol[1]) || std::isnan(sol[2]) || std::isinf(sol[2]) )
      cout << "error in rootsCubicPolynomial() - 3 solution: " << sol[0] << "  " << sol[1] << "  " << sol[2] << "  " << p << "  " << q << "  " << r << endl;
    return 3;
  }
  else // root == 0, 2 real solutions
  {
    x1=A+B;
    x2=-x1/2.0;
    sol[0] = x1 - (p/3.0);
    sol[1] = x2 - (p/3.0);
    if( std::isnan(sol[0]) || std::isinf(sol[0]) || std::isnan(sol[1]) || std::isinf(sol[1]) )
      cout << "error in rootsCubicPolynomial() - 2 solution: " << sol[0] << "  " << sol[1] << "  " << p << "  " << q << "  " << r << "  " << a << "  " << b << "  " << c << "  " << d << endl;
    return 2;
  }
}



/**
 *
 * The first 4 parameters are only needed in getMomenta23() to sample the longitudinal momentum of outgoing particle 1. 
 * It is not needed to calculate the total cross section.
 * @param[out] pz11 solution 1 of longitudinal momentum of outgoing particle 1
 * @param[out] pz12 solution 2 of longitudinal momentum of outgoing particle 1
 * @param[out] dTf1 deltaTrans_factor of solution 1
 * @param[out] dTf2 deltaTrans_factor of solution 2
 * @param[in] qt Momentum transfer qt scaled
 * @param[in] kt Momentum transfer kt scaled
 * @param[in] y Rapidity of emitted gluon.
 * @param[in] phi azimuthal angle
 * @param[in] M2 mass of the heavy quark squared scaled
 * @param[in] m1_2 mass of particle 1 squared scaled
 * @param[in] m2_2 mass of particle 2 squared scaled
 * @param[in] md2_wo_as_int Debye mass scaled by 1/s and divided by the coupling alpha_s 
 * @param[in] lambda_scaled_int Mean free path scaled by s
 * @param[in] beta_int boost beta
 * @param[in] cos_theta_int cosine of angle theta
 * @param[in] s_int Mandelstam s (needed for the running coupling since qt and kt are dimensionless without a scale)
 * @param[in] kappa_int Kappa for Debye screening
 * @param[in] formationTimeTyp_int type of formation time of radiated gluon in 2->3
 * possibilities are (just give string as input, eg. "bamps_org" ):
 * - **bamps_org**: as implemented in BAMPS, formation time in frame in which gluon is 
 *   purely transverse (Sigma'') is tau'' = 1/kt, boost lambda to this frame and compare
 * - **bamps_org_extended**: bamps_org extended to massive sector, formation time in frame 
 *   in which gluon is purely transverse is tau'' = kt / ( kt^2 + x^2 M^2), with x = k+ / P_hq+
 * - **bamps_org_extended_xwE**: as bamps_org_extended, but with definition x = omega / E_hq
 * - **compare_cm_frame**: compare formation time in center of mass frame, boost lamda to cm frame, 
 *   here tau' = (2) omega / ( kt^2 + x^2 M^2)
 * - **compare_lab_frame**: compare formation time in lab frame, boost kt, omega, x to lab frame, 
 *   here tau = (2) omega / ( kt^2 + x^2 M^2)
 * .
 * all these formulas miss a factor 2 in the formation time since it is omitted in the original 
 * BAMPS version, but according to most literature it should be there... It can be added by 
 * setting fudge_factor to 2
 * @param[in] matrixElement23_int matrix element which is used for 2->3 scattering
 * possibilities are (just give string as input, eg. "GBimproved" ):
 * - **GBimproved**: Gunion-Bertsch improved version which is correct also at forward and mid-rapidity
 * - **GBorg**: original Gunion-Bertsch version included in BAMPS (is not correct)
 * - **KPRapprox**: KPR matrix element in the limit k5->0, see notes and http://arxiv.org/abs/1109.5539
 * @param[in] md2_counter_term_in_I23_int Whether a counter term is applied which resembles the old prescription of Xu for masssless particles
 * @param[in] fudge_factor_lpm_int factor for LPM effect
 * this factor is used to play around with the cutoff
 * instead of k_t > gamma / lambda a modified cutoff k_t > X gamma / lambda
 * is used, where X is the fudge factor
 * @return integrand of the total 2->3 cross section for the given parameters
 **/
double scattering23::get_integrand_23( double& pz11, 
                                       double& pz12, 
                                       double& dTf1, 
                                       double& dTf2, 
                                       const double qt, 
                                       const double kt, 
                                       const double y, 
                                       const double phi, 
                                       const double M2, 
                                       const double m1_2, 
                                       const double m2_2, 
                                       const double md2_wo_as_int, 
                                       const double lambda_scaled_int, 
                                       const double beta_int, 
                                       const double cos_theta_int, 
                                       const double s_int, 
                                       const double kappa_int, 
                                       const string & formationTimeTyp_int, 
                                       const string & matrixElement23_int, 
                                       const bool matrixElement23_22qt_int, 
                                       const bool md2_counter_term_in_I23_int, 
                                       const double fudge_factor_lpm_int )
{
  const double epsilon=1.0e-10;
  double md2_counter_term;
  double delta,delta1,delta2;
  double E1,E3,pz3,qtkt;
  double pz21, pz22;
  double x_gluon, x_gluon_heavy_quark, x2M2, qt_kt_2, matrix_element_gQgQ, Pg;
  double mandelstam_t;
  double A, B;
  double kt2, qt2;
  double alpha_s_prefactor, md2_22;
  
  
  // test restriction due to LPM effect:
  A = beta_int * cos_theta_int;
  
  if( formationTimeTyp_int == "bamps_org"  || M2 == 0.0 )
  {
    B = kt * lambda_scaled_int * sqrt(1.0-pow(beta_int,2.0)) / fudge_factor_lpm_int;
  }
  else if( formationTimeTyp_int == "bamps_org_extended" )
  {
    if( m1_2 > 0 ) // particle 1 is heavy quark, dead cone is in positive y direction
      B = kt * lambda_scaled_int * sqrt(1.0-pow(beta_int,2.0)) / fudge_factor_lpm_int * ( 1.0 + M2 * exp(2.0 * y) );
    else
      B = kt * lambda_scaled_int * sqrt(1.0-pow(beta_int,2.0)) / fudge_factor_lpm_int * ( 1.0 + M2 * exp(-2.0 * y) );
  }
  else if( formationTimeTyp_int == "bamps_org_extended_xwE" )
  {
    B = kt * lambda_scaled_int * sqrt(1.0-pow(beta_int,2.0)) / fudge_factor_lpm_int * ( 1.0 + M2 * pow( 2.0 * cosh(y) / (1.0+M2) , 2.0 ) );
  }
  else
  {
    B = 99.9; // dummy value
    cout << "error in scattering23::get_integrand_23()" << endl;
  }
  
  //  info for scatter plot
//   cout << "A=" << A << "   B_lq=" << lambda_scaled_int * sqrt(1.0-pow(beta_int,2.0)) / fudge_factor_lpm_int  << endl;

  // constrains from LPM effect
  // as B > 0 and -1 <= A <= 1 this should not evaluate to "true"
  if ( FPT_COMP_L(A,-1.0) || FPT_COMP_G(A,1.0) || FPT_COMP_L(B,0.0) )
  {
    cout << "encountered unexpected values when determining limits for y integration" << endl;
    cout << "A = " << A << endl;
    cout << "B = " << B << endl;

    return 0.0;
  }
  
  // LPM constrain cosh(y) + A*sinh(y) < B
  if( ( cosh(y) + A * sinh(y) ) >= B )
    return 0.0;

  //some short hand notations
  kt2 = kt*kt;
  qt2 = qt*qt;
  E3 = kt*cosh(y);
  pz3 = kt*sinh(y);
  qtkt = qt*kt*cos(phi);
  
  if( matrixElement23_int == "GBimproved" )
  {
    x_gluon = kt * exp( fabs( y ) );

    if( m1_2 > 0 ) // particle 1 is heavy quark, dead cone is in positive y direction
      x_gluon_heavy_quark = kt * exp( y );
    else if( m2_2 > 0 ) // particle 2 is heavy quark, dead cone is in negative y direction
      x_gluon_heavy_quark = kt * exp( -y );
  }
  else
  {
    x_gluon = kt * exp(y);
    x_gluon_heavy_quark = x_gluon;
  }
  
  // To evaluate the integral, firstly the term 1/|dF/dy'_1| has to be calculated where F/s = 0:
  // Determine p_{1z}' from the condition F/s = 0:
  // After squaring one obtains the equation
  //
  //                 A*p_{1z}'^2 + B*p_{1z}' + C = 0                                      (1)
  // 
  // with
  //                 A = (1-E_{3}')^2 - p_{3z}' = 1 - 2*E_{3}' + k_t^2
  //
  //                 B = P_{3z}' * ( 1 - 2*E_{3}' + 2*q_t*k_t*cos(phi) + m_1^2 - m_2^2 )
  //
  //                 C = (1-E_{3}')^2 * ( q_t^2 + m_1^2 ) - 1/4 * ( 1 - 2*E_{3}' + 2*q_t*k_t*cos(phi) + m_1^2 - m_2^2 )^2
  //
  // The general solution then reads
  //
  //                 p_{1z}' = ( -B +- sqrt( B^2 - 4*A*C ) ) / 2*A
  //
  // wheras in the case of A = 0 the solution is just
  //
  //                 p_{1z}' = -C/B.
  //
  // Additionally the condition
  //
  //                 1 - 2*E_{3}' + 2*q_t*k_t*cos(theta) - 2*p_{3z}'*p_{1z}' + m_1^2 - m_2^2 >= 0          (2)
  //
  // (equivalent to   B/p_{3z}' - 2*p_{3z}'*p_{1z}' >= 0)
  //
  // has to be fulfilled (see notes).
  //

  double c;
  double deltaTrans_factor = 0;   // this is the variable that will hold 1/|dF/dy'_1| where F=0
                                  //(the factor arising from the transformation of the delta function, hence the name)                               
  
  A = 1.0 - 2.0*E3 + kt2;
  B = pz3 * ( 1.0 - 2.0*E3 + 2.0*qtkt + m1_2 - m2_2 );
  double C = pow((1.0 - E3),2.0) * ( qt2 + m1_2 ) - 0.25*pow( ( 1.0 - 2.0*E3 + 2.0*qtkt + m1_2 - m2_2 ) ,2.0);
  
  if ( FPT_COMP_E(A,0.0) )        // check for special case A = 0
  {
    if ( FPT_COMP_E( B, 0.0 ) ) // if A = 0 and B = 0 there is no solution
    {
      return 0;
    }
    else
    {
      pz11 = -C / B;            // the solution for p_{z1}' (pz of particle 1) in the case of A = 0
      pz21 = - ( pz11 + pz3 );  // pz of particle 2 for this solution, obtained from momentum conservation
      delta1 = B / pz3 - 2.0 * pz3 * pz11;
      dTf1 = 0;
      dTf2 = 0;
      // correct
//       if ( FPT_COMP_GE( delta1, 0.0 ) && pz11 > 0 && pz21 < 0 )  // check whether condition (2) is fulfilled, // why pz11>0 and pz21<0 ? no backward scattering for particles 1 and 2 allowed since small angle scattering approximation is employed which would not make much sense for backward scattering (very large angle, although qt is also small in this case)
      // incorrect
//       if ( FPT_COMP_GE( delta1, 0.0 ) && pz11 > 0 )  // check whether condition (2) is fulfilled, // why pz11 >0 ? no backward emmission allowed since small angle scattering approximation is employed which would not make much sense for backward scattering (very large angle)
      if ( ( matrixElement23_int != "GBorg" && FPT_COMP_GE( delta1, 0.0 ) && pz11 > 0 && pz21 < 0 ) || ( matrixElement23_int == "GBorg" && FPT_COMP_GE( delta1, 0.0 ) && pz11 > 0 ) )
      {
        E1 = sqrt( m1_2 + qt2 + pow( pz11, 2.0 ) );
        c = 2 * fabs( pz11 - E3 * pz11 + E1 * pz3 ); // 1/c is the argument of the sum stemming from the transformation of the delta function
        // i.e., fabs(dF/dy_1') = 2*c where F=0
        if ( c < epsilon )
          deltaTrans_factor = 1.0 / epsilon;
        else
          deltaTrans_factor = 1.0 / c;
        dTf1 = deltaTrans_factor;
      }
      else
      {
        return 0;
      }
    }
  }
  else                                        // the general case, i.e. A != 0 
  {
    delta = pow(B,2.0) - 4.0*A*C;           // the argument of sqrt(..) in the general solution for p_{z1}'
    if(delta < 0.0) 
    {
      return 0;
    }
    else
    {
      // solution 1
      pz11 = ( -B + sqrt( delta ) ) / ( 2.0 * A );
      delta1 = B / pz3 - 2.0 * pz3 * pz11;    // condition (2) requires delta1 >= 0 in order for pz11 to be a valid solution
      pz21 = - ( pz11 + pz3 );  // pz of particle 2 for this solution, obtained from momentum conservation
      
      // solution 2
      pz12 = ( -B - sqrt( delta ) ) / ( 2.0 * A );
      delta2 = B / pz3 - 2.0 * pz3 * pz12;    // condition (2) requires delta2 >= 0 in order for pz12 to be a valid solution
      pz22 = - ( pz12 + pz3 );  // pz of particle 2 for this solution, obtained from momentum conservation
      
      deltaTrans_factor = 0;
      dTf1 = 0;
      dTf2 = 0;
      // correct
//       if ( delta1 >= 0.0 && pz11 > 0 && pz21 < 0 )
      // incorrect
//       if ( delta1 >= 0.0 && pz11 > 0 )
      if( ( matrixElement23_int != "GBorg" && delta1 >= 0.0 && pz11 > 0 && pz21 < 0 ) || ( matrixElement23_int == "GBorg" && delta1 >= 0.0 && pz11 > 0 ) )
      {
        E1 = sqrt( m1_2 + qt2 + pow( pz11, 2.0 ) );
        c = 2 * fabs( pz11 - E3 * pz11 + E1 * pz3 );
        if ( c < epsilon )
          dTf1 = 1.0 / epsilon;
        else
          dTf1 = 1.0 / c;
        deltaTrans_factor += dTf1;
      }
      // correct
//       if ( delta2 >= 0.0 && pz12 > 0 && pz22 < 0 )
      // incorrect
//       if ( delta2 >= 0.0 && pz12 > 0 )
      if( ( matrixElement23_int != "GBorg" && delta2 >= 0.0 && pz12 > 0 && pz22 < 0 ) || ( matrixElement23_int == "GBorg" && delta2 >= 0.0 && pz12 > 0 ) )
      {
        E1 = sqrt( m1_2 + qt2 + pow( pz12, 2.0 ) );
        c = 2 * fabs( pz12 - E3 * pz12 + E1 * pz3 );
        if ( c < epsilon )
          dTf2 = 1.0 / epsilon;
        else
          dTf2 = 1.0 / c;
        deltaTrans_factor += dTf2;
      }
    }
  }
  
  // some abbreviations
  x2M2 = pow( x_gluon_heavy_quark ,2.0) * M2;
  qt_kt_2 = qt2 + kt2 - 2.0*qtkt;
  
  if( matrixElement23_int == "GBimproved" )
  {
    double pz1_scaled;
    
    if ( ran2() < dTf1 / ( dTf1 + dTf2 ) )
      pz1_scaled = pz11;
    else
      pz1_scaled = pz12;//!!!
    
    // calculate Mandelstam t from GB coordinates
    const double pz1_in = ( 1.0 - M2 ) / 2.0; //incoming particle 1 z momentum scaled by sqrt(s)
    const double pz2_in = - ( 1.0 - M2 ) / 2.0; //incoming particle 2 z momentum scaled by sqrt(s)
    const double pz1_out = pz1_scaled; //incoming particle 2 z momentum scaled by sqrt(s)
    const double pz2_out = - ( pz1_scaled + pz3 ); //incoming particle 2 z momentum scaled by sqrt(s)
    double E1_in; //incoming particle 1 energy scaled by sqrt(s)
    double E2_in; //incoming particle 2 energy scaled by sqrt(s)
    double E1_out; //outgoing particle 1 energy scaled by sqrt(s)
    double E2_out; //outgoing particle 2 energy scaled by sqrt(s)
    double mass_term, mass_term_dash; // used for P_i squared + P_i_dash squared ( which is twice of the squared mass of the particle)
    if( m1_2 > 0 ) // particle 1 is heavy quark
    {
      E1_in = ( 1.0 + M2 ) / 2.0;
      E2_in = ( 1.0 - M2 ) / 2.0;
      E1_out = sqrt( qt2     + pow( pz1_out , 2.0 ) + M2 );
      E2_out = sqrt( qt_kt_2 + pow( pz2_out , 2.0 ) );
      mass_term = 2.0 * M2;
      mass_term_dash = 0;
    }
    else
    {
      E1_in = ( 1.0 - M2 ) / 2.0;
      E2_in = ( 1.0 + M2 ) / 2.0;
      E1_out = sqrt( qt2     + pow( pz1_out , 2.0 ) );
      E2_out = sqrt( qt_kt_2 + pow( pz2_out , 2.0 ) + M2 );
      mass_term = 0.0;
      mass_term_dash = 2.0 * M2;
    }
    
    // mandelstam t scaled by s
    mandelstam_t      = mass_term      - 2.0 * ( E1_in * E1_out - pz1_in * pz1_out );
    double mandelstam_t_dash 
                      = mass_term_dash - 2.0 * ( E2_in * E2_out - pz2_in * pz2_out );
    
    // for negative y switch pt of particles 3 and 4 
    if( y < 0 )
    {
      double qt2_temp = qt2;
      qt2 = qt_kt_2;
      qt_kt_2 = qt2_temp;
      
      qtkt = - qt * kt * cos(phi) + kt2;

      mandelstam_t = mandelstam_t_dash;
    }
  }
  
  
  
  if( matrixElement23_int == "GBimproved" )
  {
    // set running alpha_s prefactor
    alpha_s_prefactor = coupling::get_coupling( mandelstam_t * s_int ) * coupling::get_coupling( mandelstam_t * s_int ) * coupling::get_coupling( kt2*s_int ); //!! at which scale should alpha_s be evaluated?
    
    // get matrix_element of 2->2 without prefactor
    md2_22 = kappa_int * md2_wo_as_int * coupling::get_coupling( mandelstam_t * s_int ); //!! scale of running alpha_s?
    
    if( matrixElement23_22qt_int )
      matrix_element_gQgQ = 1.0 / pow(( qt2 + md2_22 ),2.0) * pow( 1.0 - M2 , 2.0 ); // scaled by s without prefactor
    else
      matrix_element_gQgQ = 1.0 / pow( ( mandelstam_t - md2_22 ) , 2.0 ) * pow( 1.0 - M2 , 2.0 ); // scaled by s without prefactor
  }
  else
  {
    // set running alpha_s prefactor
    alpha_s_prefactor = coupling::get_coupling( -qt2*s_int ) * coupling::get_coupling( -qt2*s_int ) * coupling::get_coupling( kt2*s_int ); //!! at which scale should alpha_s be evaluated?
    
    // get matrix_element of 2->2 without prefactor
    md2_22 = kappa_int * md2_wo_as_int * coupling::get_coupling( -qt2*s_int ); //!! scale of running alpha_s?
    
    matrix_element_gQgQ = 1.0 / pow(( qt2 + md2_22 ),2.0) * pow( 1.0 - M2 , 2.0 ); // scaled by s without prefactor
  }
  
  const double md2_ktqt = md2_wo_as_int * coupling::get_coupling( -qt_kt_2*s_int ); //!!  what scale?
  
  if( matrixElement23_int == "GBimproved" )
    Pg = pow( 1.0 - x_gluon , 2.0 ) * ( kt2 / pow( kt2 + x2M2 ,2.0) + qt_kt_2 / pow( qt_kt_2 + md2_ktqt + x2M2 ,2.0) + 2.0 * (qtkt - kt2 ) / (kt2 + x2M2) / (qt_kt_2 + md2_ktqt + x2M2) );
  else if( matrixElement23_int == "GBorg" )
    Pg = kt2 / pow( kt2 + x2M2 ,2.0) + qt_kt_2 / pow( qt_kt_2 + md2_ktqt + x2M2 ,2.0) + 2.0 * (qtkt - kt2 ) / (kt2 + x2M2) / (qt_kt_2 + md2_ktqt + x2M2);
  else if( matrixElement23_int == "KPRapprox" )
    Pg = 1.0 / kt2 * ( ( 1.0 - 2.0 * M2 ) - qt2 * ( 10.0 / 9.0 + 16.0 / 9.0 / ( 1.0 - M2 ) - 1.0 / pow( 1.0 - M2 , 2.0 ) ) );
  else if( matrixElement23_int == "KPRapprox_screen" )
    Pg = kt2 / pow( kt2 + x2M2 ,2.0) * ( ( 1.0 - 2.0 * M2 ) - qt2 * ( 10.0 / 9.0 + 16.0 / 9.0 / ( 1.0 - M2 ) - 1.0 / pow( 1.0 - M2 , 2.0 ) ) );
  else
  {
    Pg = 99.9; // dummy value
    cout << "error get_integrand_23()" << endl;
  }
  
  
  // for M=0 this term does not simplify to Xu's formula for light particles since he screened with a Debye mass after some simplifications which cannot be performed here. However, I think my way of screening is more appropiate (and also Xu agreed). But to compare to older calculations one needs this artificial counter-term which makes my result the same as Xu's for the massless case. The more correct treatment is md2_counter_term = 0.0. However the impact of this counter-term should be small -> check this!!
  if( md2_counter_term_in_I23_int )
  {
    if( matrixElement23_int == "GBimproved" )
      md2_counter_term = pow( 1.0 - x_gluon , 2.0 ) * ( ( md2_ktqt ) / ( kt2 + x2M2 ) / ( qt_kt_2 + md2_ktqt + x2M2 ) - ( md2_ktqt ) / pow( qt_kt_2 + md2_ktqt + x2M2 ,2.0) );
    else if( matrixElement23_int == "GBorg" )
      md2_counter_term = ( md2_ktqt ) / ( kt2 + x2M2 ) / ( qt_kt_2 + md2_ktqt + x2M2 ) - ( md2_ktqt ) / pow( qt_kt_2 + md2_ktqt + x2M2 ,2.0);
    else
    {
      std::string errMsg = "Error in scattering23. md2_counter_term only allowed for GBimproved and GBorg.";
      throw eScatt23_error( errMsg );
    }
  }
  else
    md2_counter_term = 0.0;
  
  return ( alpha_s_prefactor * deltaTrans_factor *  matrix_element_gQgQ * ( Pg - md2_counter_term ) );
}








/**
 * Sets new momenta for the outgoing particle according to the values sampled in getNewMomenta23.
 * The results are written to the vectors P1, P2 and P3.
 *
 * @param[out] P1 Momentum vector of outgoing particle 1
 * @param[out] P2 Momentum vector of outgoing particle 2
 * @param[out] P3 Momentum vector of outgoing particle 3
 * @param[in] R1 Space-Time vector of ingoing particle 1
 * @param[in] R2 Space-Time vector of ingoing particle 2
 * @param[in] PT1 Momentum transfer qt
 * @param[in] PT3 Momentum transfer kt
 * @param[in] y3 Rapidity of emitted gluon.
 * @param[in] phi azimuthal angle
 * @param[in] PZ1 Longitudinal momentum transfer.
 **/
void scattering23::setNewMomenta23( VectorEPxPyPz & P1, 
                                    VectorEPxPyPz & P2,
                                    VectorEPxPyPz & P3,
                                    const VectorTXYZ & R1,
                                    const VectorTXYZ & R2,
                                    const double PT1, const double PT3, 
                                    const double y3, const double phi, const double PZ1 )
{
  double PZ3, sinus, cosinus;

  VectorTXYZ R1cell, R2cell, R1cm, R2cm;
  VectorEPxPyPz P3cell, P3cm;
  VectorEPxPyPz PP, TT, transv;

  LL_cell.boost(R1, R2, R1cell, R2cell);
  LL_CM.boost( R1cell, R2cell, R1cm, R2cm);

  rotation( P1cm, R1cm, R2cm, PP, TT );

  P1cm = TT * PT1 + PP * PZ1;
  P1cm(0) = sqrt( P1cm.vec2() + pow( m1 , 2.0 ) );

  TT = -TT;

  sinus = sin( phi );
  cosinus = cos( phi );
  transv = TT * cosinus + Cross( PP, TT ) * sinus;

  TT = transv;

  PZ3 = PT3 * sinh( y3 );

  P3cm = TT * PT3 + PP * PZ3;
  P3cm(0) = sqrt( P3cm.vec2() );

  P2cm = -( P1cm + P3cm );
  P2cm(0) = sqrt( P2cm.vec2() + pow( m2 , 2.0 ) );

  // boost new momentum vectors back to original frame

  LL_CM.boostInv(P1cm,P2cm,P3cm, P1cell,P2cell,P3cell);

  if ( randomizedConfiguration == 0 )
  {
    LL_cell.boostInv(P1cell,P2cell,P3cell, P1,P2,P3);
  }
  else
  {
    LL_cell.boostInv(P2cell,P1cell,P3cell, P1,P2,P3);
  }
}



/**
 * Returns the maximal and minimal values of the rapidity y possible
 * for ANY k_t (within the appropriate ranges). These values are
 * returned via the ymin and ymax passed to maximal_Y_range_23 as
 * references. 
 *
 * for formation time which do not depend on x and therefore not on y
 * (formationTimeTyp == "bamps_org"), the solution of the 2 curves  
 * \f$ g(y) = kt = ( s - M^2 ) / ( 2 \sqrt{s} \cosh y ) \f$
 * \f$ f(y) = kt = ( \cosh y + A \sinh y ) / B \f$
 * are
 * \f$ y_{\rm{min/max}} = 1/2 \ln( ( C \pm \sqrt( C^2 + A^2 - 1 ) ) / ( 1 + A ) ) \f$
 *
 * for other formation time types it is a bit more complivated, see notes.
 *
 * @param[out] ymin Returns the minmial y value.
 * @param[out] ymax Returns the maximal y value.
 * @return Whether a valid range could be successfully determined.
 **/
bool scattering23::maximal_Y_range_23( double& ymin, double& ymax ) const
{
  double m1_2 = pow( m1/sqrtS ,2.0); // scaled by s
  double m2_2 = pow( m2/sqrtS ,2.0); // scaled by s
  
  double M2;
  if( m1_2 >= m2_2 && m2_2 == 0.0 )
    M2 = m1_2;
  else if( m2_2 >= m1_2 && m1_2 == 0.0 )
    M2 = m2_2;
  else
  {
    M2 = 99.9; // dummy value
    cout << "error in maximal_Y_range_23(). Both masses are larger than 0: " << m1_2 << "  " << m2_2 << endl;
  }
  
  double A = beta * cos_theta;
  double B_lq_wo_kt = lambda_scaled * sqrt(1.0-pow(beta,2.0)) / fudge_factor_lpm; // B for light quarks divided by kt
  
  if( formationTimeTyp == "bamps_org"  || M2 == 0.0 )
  {
    double C = ( 1.0 - M2 ) * B_lq_wo_kt - 1.0;
    double delta = C*C + A*A - 1.0;
    
    if( delta < 0.0 )  // no solution, LPM curve is always larger than kinematic curve
    {
      ymin = ymax = 0.0;
      return false; 
    }

    ymin = 0.5 * log( ( C - sqrt( delta ) ) / ( 1.0 + A ) );
    ymax = 0.5 * log( ( C + sqrt( delta ) ) / ( 1.0 + A ) );
    return true;
  }
  else if( formationTimeTyp == "bamps_org_extended" || formationTimeTyp == "bamps_org_extended_xwE"  )
  {
    double a, b, c;

    // intersection points are solutions to a u^2 + b u + c = 0 with
    if( formationTimeTyp == "bamps_org_extended" )
    {
      if( m1_2 > 0 )
      {
        a = 1.0 + A - 2.0 * B_lq_wo_kt * ( M2 - M2*M2 );
        b = 2.0 * ( 1.0 - B_lq_wo_kt * ( 1.0 - M2 ) );
        c = 1.0 - A;
      }
      else
      {
        a = 1.0 + A;
        b = 2.0 * ( 1.0 - B_lq_wo_kt * ( 1.0 - M2 ) );
        c = 1.0 - A - 2.0 * B_lq_wo_kt * ( M2 - M2*M2 );
      }
    }
    else if( formationTimeTyp == "bamps_org_extended_xwE" )
    {
      a = 1.0 + A - 2.0 * B_lq_wo_kt * ( M2 - M2*M2 ) / pow( 1.0 + M2 , 2.0 );
      b = 2.0 * ( 1.0 - B_lq_wo_kt * ( 1.0 - M2 ) * ( 1.0 + 2.0 * M2 / pow( 1.0 + M2 , 2.0 ) ) );
      c = 1.0 - A - 2.0 * B_lq_wo_kt * ( M2 - M2*M2 ) / pow( 1.0 + M2 , 2.0 );
    }
    else
    {
      // you should never get here ! ;)
      a = b = c = 99.9; // dummy value
    }
    
    if( FPT_COMP_E(a,0.0) )
    {
      cout << "error in maximal_Y_range_23(). a = " << a << endl;
    }
    
    double delta = b*b - 4.0 * a * c;
    
    if( delta < 0.0 )  // no solution, LPM curve is always larger than kinematic curve
    {
      ymin = ymax = 0.0;
      return false;
    }

    // solutions u_{1,2} = \frac{-b \pm \sqrt{b^2-4ac}}{2a}
    double umin = ( - b - sqrt(delta) ) / ( 2.0 * a );
    double umax = ( - b + sqrt(delta) ) / ( 2.0 * a );


    const double y_infinity = 10.0; // large value for rapidity, basically infinity
    
    if( umin > 0.0 && umax > 0.0 ) // u is positive for both solutions: two solutions for y = 0.5 log( u ), the general case
    {
      ymin = 0.5 * log( umin );
      ymax = 0.5 * log( umax );
      
      if( ymin > ymax )
      {
        cout << "problem in scattering23::maximal_Y_range_23()" << endl;
        cout << ymin << "  " << ymax << "  " << A << "  " << B_lq_wo_kt << endl;
        double ytemp = ymin;
        ymin = ymax;
        ymax = ytemp;
      }
    }
    else if( umin > 0.0 && umax <= 0.0 ) // umax is smaller than zero, no solution for ymax
    {
      ymin = 0.5 * log( umin );
      ymax = y_infinity;
      
      if( formationTimeTyp == "bamps_org_extended_xwE" || ( formationTimeTyp == "bamps_org_extended" && m2_2 > 0.0 ) )
        cout << "error in maximal_Y_range_23(). Case 2 should not happen for " << formationTimeTyp << endl;
    }
    else if( umin <= 0.0 && umax > 0.0 ) // umin is smaller than zero, no solution for ymin
    {
      ymin = -y_infinity;
      ymax = 0.5 * log( umax );
      
      if( formationTimeTyp == "bamps_org_extended" && m1_2 > 0.0 )
        cout << "error in maximal_Y_range_23(). Case 3 should not happen for " << formationTimeTyp << endl;
    }
    else // no intersection point at all, LPM curve is always smaller than kinematic curve -> all rapidities are allowed
    {
      ymin = -y_infinity;
      ymax =  y_infinity;
      
      if( formationTimeTyp == "bamps_org_extended" && fudge_factor_lpm != 0.0 )
        cout << "error in maximal_Y_range_23(). Case 4 should not happen for " << formationTimeTyp << ". But happens if i23 interpolation gives non-zero value even if i23 should be zero." << endl;
    }
  }

//   cout << "ymin = " << ymin << "   ymax = " << ymax << "  sqrtS = " << sqrtS << endl;
  
  return true;
}




/**
 * Auxiliary routine. Determines a direction for the transverse momtentum vector in the azimuthal plane.
 * Used by #scattering23::setNewMomenta22
 *
 * @param[in] P Momentum vector of ingoing particle 1
 * @param[in] R1 Space-time vector of ingoing particle 1
 * @param[in] R2 Space-time vector of ingoing particle 2
 * @param[out] PP Direction of momentum (i.e. normalized P)
 * @param[out] TT Direction in the azimuthal plane between the two
 * particles.
 *
 * This routine is identical with scattering22::rotation
 */

void scattering23::rotation( const VectorEPxPyPz & P, const VectorTXYZ & R1, const VectorTXYZ & R2, VectorEPxPyPz & PP, VectorEPxPyPz & TT ) const
{
  double c1 = P.vec2();

  VectorTXYZ dR = R2 - R1;
  double c2 = Dot3(P, dR);

  PP = P * (1.0/sqrt( c1 ));
  TT = P * c2 - dR * c1;

  double c3 = sqrt( TT.vec2() );

  // r2-r1 // p1 => TT[] = 0
  if ( c3 < 1.0e-8 )
  {
    double min = 1.0;
    int mini = -1;

    for ( int i = 1;i <= 3;i++ )
    {
      if ( fabs( PP(i) ) < min )
      {
        min = fabs( PP(i) );
        mini = i;
      }
    }
    switch ( mini )
    {
    case 1:
      TT(1) = 0.0;
      TT(2) =  PP(3);
      TT(3) = -PP(2);
      break;
    case 2:
      TT(1) = -PP(3);
      TT(2) = 0.0;
      TT(3) =  PP(1);
      break;
    case 3:
      TT(1) =  PP(2);
      TT(2) = -PP(1);
      TT(3) = 0.0;
      break;
    default:
      cout << "Error in rotation()" << endl;
    }
  }
  //----------------------

  double phi = 2.0 * M_PI * ran2();
  double sinus = sin( phi );
  double cosinus = cos( phi );

  VectorEPxPyPz transv = TT * cosinus + Cross(PP, TT) * sinus;

  c3 = transv.vec2();
  TT = transv * (1.0/sqrt( c3 ));

}

