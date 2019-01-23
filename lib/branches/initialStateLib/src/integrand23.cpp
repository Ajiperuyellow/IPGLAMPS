//---------------------------------------------------------------------------------------
//provided by subversion
//---------------------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/integrand23.cpp $
//$LastChangedDate: 2013-11-22 16:49:57 +0100 (Fr, 22. Nov 2013) $
//$LastChangedRevision: 1551 $
//$LastChangedBy: uphoff $
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------


#include <math.h>
#include <iostream>

#include "integrand23.h"
#include "FPT_compare.h"
#include "scattering23.h"

using std::cout;
using std::endl;



// overloaded operator() that makes integrand23 a function object - for use with CUBA-vegas
void integrand23::operator()(const int *ndim, const double xx[], const int *ncomp, double ff[]) const
{
  double wgt = 0;
  ff[0] = this->operator()(xx, wgt);  
}


/**
 * This function operator returns the matrix element 2->3 for given qt2, kt2, y and phi (qt2 and kt2 scaled by s to be dimensionless).
 * Arguments are expected to be in the range [0,1], thus indicating the relative position e.g. between qt2_min and qt2_max.
 * This scaled treatment is needed for use with numerical integration routines.
 *
 * @param[in] xx[] unit-offset array with function arguments (xx[1] = qt2, xx[2] = kt2, xx[3] = y, xx[4] = phi)
 * @param[in] wgt there for compatibility reasons - not used at the moment
 */
double integrand23::operator()(const double xx[], double wgt) const
{
  double integrand_23;

  
  double substitution_prefactor = 1.0;
  
  double m1_2 = m1*m1; // scaled by s
  double m2_2 = m2*m2; // scaled by s

  double M2 = ( m1_2 > m2_2 ) ? m1_2 : m2_2;

  if ( m1_2 > 0.0 && m2_2 > 0.0 )
  {
    std::string errMsg("Error in integrand23::operator(). Both masses are larger than 0: ");
    errMsg += m1_2;
    errMsg += " ";
    errMsg += m2_2;
    throw eInt23_error( errMsg );
  }

  //--------------------------------------------
  // q_t^2 is integrated from 0 to 1/4*(1-M^2)^2
  // integrate in ln(qt2) to cure the divergence -> numerically easier 

  double ln_qt2_min = -20.0; // small value
  double ln_qt2_max = log( 0.25 * pow( 1.0 - M2 ,2.0) );
  if (ln_qt2_min >= ln_qt2_max) { return 0; }

  double V = (ln_qt2_max - ln_qt2_min);
  double ln_qt2 = ln_qt2_min + (ln_qt2_max - ln_qt2_min)*xx[1];
  double qt2 = exp(ln_qt2); // substitution 
  double qt = sqrt(qt2);
  
  substitution_prefactor = substitution_prefactor * qt2;

//   qt2_min = 0.0;
//   qt2_max = 0.25 * pow( 1.0 - M2 ,2.0); // is this factor 0.25 really correct?? yes it is, see notes
//   if (qt2_min >= qt2_max)
//     return 0;
//       
//   V = (qt2_max - qt2_min);
//   qt2 = qt2_min + (qt2_max - qt2_min)*xx[1];
//   qt = sqrt(qt2);
  //--------------------------------------------


  //--------------------------------------------
  // k_t^2 is integrated from 1/lambda^2 to 1/4*(1-M^2)^2
  // integrate in ln(qt2) to cure the divergence. This substitution is
  // useful since the 1/kt2 factor in the integrand vanishes. Then the
  // numerical integration is faster (poles are difficult to
  // integrate).  

  double ln_kt2_min = -20.0; // small value
  if( ( formationTimeTyp_int == "bamps_org"  || M2 == 0.0 ) && fudge_factor_lpm_int != 0.0 )
  {
    ln_kt2_min = log( pow( fudge_factor_lpm_int / lambda_scaled_int , 2.0) );
  }
  else if( formationTimeTyp_int == "bamps_org_extended" || formationTimeTyp_int == "bamps_org_extended_xwE" || fudge_factor_lpm_int == 0.0 )
  {
    ln_kt2_min = -20.0; // small value
  }
  else
  {
    std::string errMsg("Error in integrand23::operator(). kt!!! ");
    throw eInt23_error( errMsg );
  }

  double ln_kt2_max = log( 0.25 * pow( 1.0 - M2 ,2.0) );
  if (ln_kt2_min >= ln_kt2_max)
    return 0;

  V = V*(ln_kt2_max - ln_kt2_min);
  double ln_kt2 = ln_kt2_min + (ln_kt2_max - ln_kt2_min)*xx[2];
  double kt2 = exp(ln_kt2); // substitution 
  double kt = sqrt(kt2);
  
  substitution_prefactor = substitution_prefactor * kt2;
  
//   if( formationTimeTyp_int == "bamps_org"  || M2 == 0.0 )
//     kt2_min = pow( fudge_factor_lpm_int / lambda_scaled_int , 2.0);
//   else if( formationTimeTyp_int == "bamps_org_extended" || formationTimeTyp_int == "bamps_org_extended_xwE" )
//     kt2_min = 0.0;
//   else
//     cout << "error in integrand23::operator()" << endl;
//   kt2_max = 0.25 * pow( 1.0 - M2 ,2.0);
//   if (kt2_min >= kt2_max)
//     return 0;
// 
//   V = V*(kt2_max - kt2_min);
//   kt2 = kt2_min + (kt2_max - kt2_min)*xx[2];
//   kt = sqrt(kt2);
  //--------------------------------------------
  
  
  //--------------------------------------------
  // get y range
  
  double y_min, y_max;
  bool solvable = scattering23::y_range_23( y_min, y_max, kt, M2, m1_2, m2_2, lambda_scaled_int, beta_int, cos_theta_int, formationTimeTyp_int, fudge_factor_lpm_int );

  if ( !solvable ) return 0;

  V = V*(y_max - y_min);
  double y = y_min + (y_max - y_min)*xx[3];

  //------------------------------------------>>


  //<<------------------------------------------
  // phi is integrated from 0 to Pi
  double phi_min = 0;
  double phi_max = M_PI;
  if (phi_min >= phi_max)
    return 0;

  V = V*(phi_max - phi_min);
  double phi = phi_min + (phi_max - phi_min)*xx[4];
  //------------------------------------------>>


  // Now that q_t, k_t, y and phi have been fixed, get
  // the argument of the sigma_23 integral at these values.

  double dummy1, dummy2, dummy3, dummy4;

  integrand_23 = scattering23::get_integrand_23( dummy1, dummy2, dummy3, dummy4, qt, kt, y, phi, M2, m1_2, m2_2, md2_wo_as_int, lambda_scaled_int, beta_int, cos_theta_int, s_int, kappa_int, formationTimeTyp_int, matrixElement23_int, matrixElement23_22qt_int, md2_counter_term_in_I23_int, fudge_factor_lpm_int );
  
  return (substitution_prefactor * V * integrand_23);
}

