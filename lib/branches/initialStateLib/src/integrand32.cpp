//---------------------------------------------------------------------------------------
//provided by subversion
//---------------------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/integrand32.cpp $
//$LastChangedDate: 2013-10-29 22:10:47 +0100 (Di, 29. Okt 2013) $
//$LastChangedRevision: 1530 $
//$LastChangedBy: gallmei $
//---------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------


#include <math.h>

#include "integrand32.h"
#include "FPT_compare.h"
#include "coupling.h"


void integrand32::set_beta_vec(const vector4D & x_)
{
  beta_vec[0] = 0;
  for(int i=1;i<=3;i++)
  {
    beta_vec[i] = x_(i);
    //beta_vec[i] = 0;  // for testing consistency with old implementation (also see scattering32.cpp)
  }
  
  beta_abs = sqrt ( pow(beta_vec[1],2.0) + pow(beta_vec[2],2.0) + pow(beta_vec[3],2.0) );
}


// overloaded operator() that makes integrand32 a function object - for use with CUBA-vegas
void integrand32::operator()(const int *ndim, const double xx[], const int *ncomp, double ff[]) const
{
  double wgt = 0;
  ff[0] = this->operator()(xx, wgt);  
}


// overloaded operator() that makes integrand32 a function object - for use with NR-vegas
double integrand32::operator()(const double xx[], double wgt) const
{
  double sin_gamma = sqrt( 1.0 - pow(cos_gamma_int,2.0) );   //sin(gamma)
  
  double u = xx[1];                       // cos(theta)                                   
  double us = sqrt( 1.0 - pow(u,2.0) );  // sin(theta)
  
  double V_phi = M_PI;  // volume for phi integration, as phi is sampled from 0 to M_PI rather than from 0 to 1
  double phi = M_PI*xx[2];
  double v = cos(phi);
  
  double cos_delta = sin_gamma * us * v + cos_gamma_int * u;   // cos(delta) = sin(gamma)*sin(theta)*cos(phi) + cos(gamma)*cos(theta)
  double sin_delta_2 = 1.0 - pow( cos_delta,2.0 );             // sin(delta)^2 = 1 - cos(delta)^2
  
  if ( E3_int * cos_delta + E1_int * u < 0 )
  {
    return 0;
  }
  
  double kt2 = pow(E3_int,2.0) * sin_delta_2;                  // kt^2 = E3^2 * sin(delta)^2 
  
  // Theta (capital!) is the angle between the boost-vector beta' (here: beta_vec) and the vector p1'
  double cos_Theta;
  if ( FPT_COMP_E(beta_abs,0.0) )
  {
    cos_Theta = 0;
  }
  else
  {
    cos_Theta = 1/beta_abs * ( beta_vec[1]*us*v + beta_vec[2]*us*sqrt(1-pow(v,2.0)) + beta_vec[3]*u );
  }
  
  // the constraint for kt^2 depending on E3, Theta, delta, beta', sqrt(s) and lambda
  double constraint = fudge_factor_lpm_int/lambda_scaled_int * E3_int/sqrt(1 - pow(beta_abs,2.0)) * (1 + (beta_abs * cos_delta * cos_Theta)); 
  
  if (kt2 > constraint)
  {
    double qt2 = pow(E1_int,2.0) * pow(us,2.0);                  // qt^2 = E1^2 * sin(theta)^2 
    //vector product qt*kt = -E1*E3*sin(theta)*( cos(gamma)*sin(theta) - sin(gamma)*cos(theta)*cos(phi) )
    double qtkt = -E1_int * E3_int * us * ( cos_gamma_int*us - sin_gamma*u*v );  
    double qt_kt_2 = qt2 + kt2 - 2 * qtkt;
    double qtkt_kt2 = qtkt - kt2;
    
    // The term ( 1 - xbar )^2 where xbar = E3 / sqrt(s) * exp(|y|).
    // With the help of some analytic rearrangements this can be cast into the form used below,
    // that is easier (faster) to compute than using xbar = E3 / sqrt(s) * exp(arccosh(1/sin(delta))) directly.
    double factor_1minusx_squared = pow( 1 - E3_int * ( 1 + sqrt( 1 - sin_delta_2 ) ) , 2 );
    
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
      mandelstam_t = E1_int * ( u - 1 );   //  t = (p1' - p1)^2 (as four vectors)
    }
    else
    {
      double E2 = sqrt( pow( E1_int, 2 ) + 2 * E1_int * E3_int * cos_gamma_int + pow( E3_int, 2 ) );  
      mandelstam_t = E3_int * cos_delta + E1_int * u - E2; //  t = (p2' - p2)^2 (as four vectors)
    }
    
    if ( matrixElement32_22qt )
    {
      const double md2_22 = md2_int / coupling::get_constant_coupling() * coupling::get_coupling( - qt2 * pow( sqrtS_int , 2.0 ) ); //!! scale of running alpha_s?
      matrixElement22_factor = 1 / pow(( qt2 + md2_22 ), 2.0 );
    }
    else
    {
      const double md2_22 = md2_int / coupling::get_constant_coupling() * coupling::get_coupling( mandelstam_t * pow( sqrtS_int , 2.0 ) ); //!! scale of running alpha_s?
      matrixElement22_factor = 1.0 / pow( ( mandelstam_t - md2_22 ) , 2.0 );
    }
    
    const double md2_ktqt = md2_int / coupling::get_constant_coupling() * coupling::get_coupling( -qt_kt_2 * pow( sqrtS_int , 2.0 ) ); //!!  what scale?
    
    if ( md2_use_simplified_GB )
    {
      Pg = qt2 / ( kt2 * ( qt_kt_2 + md2_ktqt ) );  
    }
    else
    {
      Pg = 1 / kt2 + qt_kt_2 / pow( qt_kt_2 + md2_ktqt, 2 ) + 2 * ( qtkt_kt2 ) / ( kt2 * ( qt_kt_2 + md2_ktqt ) );
    }
    
    double alpha_s_prefactor = coupling::get_coupling( mandelstam_t * pow( sqrtS_int , 2.0 ) ) * coupling::get_coupling( mandelstam_t * pow( sqrtS_int , 2.0 ) ) * coupling::get_coupling( kt2 * pow( sqrtS_int , 2.0 ) ); //!! at which scale should alpha_s be evaluated?
    
    return alpha_s_prefactor * factor_1minusx_squared * V_phi * matrixElement22_factor * Pg;  // the matrix element
  }
  else
  {
    return 0;
  }
}
