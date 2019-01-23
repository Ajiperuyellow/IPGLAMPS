//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/integrand32.h $
//$LastChangedDate: 2014-05-30 15:30:33 +0200 (Fr, 30. Mai 2014) $
//$LastChangedRevision: 1733 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


#ifndef INTEGRAND32_H
#define INTEGRAND32_H

#include "vegas.h"
#include "bampsvector.h"

#define NDIM_32 2
#define NCOMP_32 1


// Function object derived from prototype integrand (defined in vegas.h).
// It calculates (via the overloaded ()-operator) the matrix element for 3->2 scattering.

class integrand32 : public integrand
{
  public:
    integrand32() : 
      beta_abs(),
      beta_vec(),
      E1_int(),
      E3_int(),
      cos_gamma_int(),
      md2_int(),
      sqrtS_int(),
      lambda_scaled_int(),
      md2_use_simplified_GB( false ), 
      matrixElement32_22qt( false ),
      fudge_factor_lpm_int()
    { };

  ~integrand32() {};
  
  // overloaded operator() that makes integrand32 a function object - for use with NR-vegas
  double operator()(const double [], double) const;
  
  // verloaded operator() that makes integrand32 a function object - for use with CUBA-vegas
  void operator()(const int*, const double [], const int*, double []) const;
  
  void set_E1(const double x_) {E1_int = x_;}
  void set_E3(const double x_) {E3_int = x_;}
  void set_cos_gamma(const double x_) {cos_gamma_int = x_;}
  void set_md2(const double x_) {md2_int = x_;}
  void set_sqrtS(const double x_) {sqrtS_int = x_;}
  void set_lambda(const double x_) {lambda_scaled_int = x_;}
  
  void set_beta_vec(const vector4D & x_);
  
  void set_md2_use_simplified_GB( const bool x_ ) { md2_use_simplified_GB = x_; }
  void set_matrixElement32_22qt( const bool x_ ) { matrixElement32_22qt = x_; }

  /** @brief Utility routine to set factor for LPM effect to play around
  * @param[in] x_ fudge_factor_lpm_int
  **/
  void set_fudge_factor_lpm(const double x_) {fudge_factor_lpm_int = x_;}
  
private:
  // auxiliary member variables needed for the calculation of the matrix element
  // set via the public routines integrand32::set_???(???)
  double beta_abs;
  double beta_vec[4];
  
  double E1_int;
  double E3_int;
  double cos_gamma_int;
  double md2_int;
  double sqrtS_int;
  double lambda_scaled_int;
  
  /** @brief Whether the version of GB is used that is simplified prior to applying the screening with md2
   * This is equivalent to the option md2_counter_term_in_I23 in scattering23. Setting it to true is equivalent to using the "old" version of GB (but including the (1-x) term).
   */
  bool md2_use_simplified_GB;
    
  /** @brief Whether the gluon propagator in the 2->2 part of the radiative cross section is approximated by 1/qt^4. If false the more exact expression 1/t^2 is employed. */
  bool matrixElement32_22qt;    

  /** @brief factor for LPM effect to play around */
  double fudge_factor_lpm_int;
};


#endif
