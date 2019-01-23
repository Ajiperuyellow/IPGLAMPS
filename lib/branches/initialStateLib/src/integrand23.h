//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/integrand23.h $
//$LastChangedDate: 2014-05-30 15:30:33 +0200 (Fr, 30. Mai 2014) $
//$LastChangedRevision: 1733 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file 
* @brief Declarations for class integrand23.
*/

#ifndef INTEGRAND23_H
#define INTEGRAND23_H

#include <stdexcept>
#include <string>
#include "vegas.h"

#define NDIM_23 4
#define NCOMP_23 1


/**
 * @brief Used for calculating (and integrating over) the matrix element for 2->3 scatterings
 *
 * This class serves as a function object and is derived from the prototype integrand (see vegas.h).
 * Via the overloaded ()-operator an object of this class can be "called" just as it was a function.
 * Utility routines are provided to set certain parameters.
 */
class integrand23 : public integrand
{
  public:
    /** @brief standard constructor */
    integrand23() : 
      md2_wo_as_int(),
      lambda_scaled_int(),
      beta_int(),
      cos_theta_int(),
      m1( 0.0 ),
      m2( 0.0 ),
      s_int(),
      kappa_int(),
      formationTimeTyp_int(),
      matrixElement23_int(),
      matrixElement23_22qt_int(),
      md2_counter_term_in_I23_int(),
      fudge_factor_lpm_int()
    { };
    /** @brief standard destructor */
    ~integrand23() {};
    
    /** @brief Overloaded operator() that makes integrand23 a function object - for use with NR-vegas */
    double operator()(const double [], double) const;
    
    /** @brief Overloaded operator() that makes integrand23 a function object - for use with CUBA-vegas */
    void operator()(const int*, const double [], const int*, double []) const;
    
    /** @brief Utility routine to set the Debye mass (squared) scaled by 1/s
     * @param[in] x_ Debye mass squared in GeV^2 divided by mandelstam s (GeV^2)
     **/
    void set_md2_wo_as(const double x_) {md2_wo_as_int = x_;}
    
    /** @brief Utility routine to set the mean free path
     * @param[in] x_ sqrt(s) in GeV * mean free path in 1/GeV (dimensionless)
     **/
    void set_lambda(const double x_) {lambda_scaled_int = x_;}
    
    /** @brief Utility routine to set the boost velocity from lab frame to CMS
     * @param[in] x_ boost velocity beta from lab frame to CMS
     **/
    void set_beta(const double x_) {beta_int = x_;}
    
    /** @brief Utility routine to set the angle between beta and CMS-axis
     * @param[in] x_ cos(theta) where theta is the angle between the boost velocity vector and the axis of the CMS system
     **/
    void set_cos_theta(const double x_) {cos_theta_int = x_;}
    
    /** @brief Utility routine to set mass of particle one
     * @param[in] x_ mass of particle one divided by sqrt(s) (dimensionless)
     **/
    void set_m1(const double x_) {m1 = x_;}
    
    /** @brief Utility routine to set mass of particle one
     * @param[in] x_ mass of particle one divided by sqrt(s) (dimensionless)
     **/
    void set_m2(const double x_) {m2 = x_;}
    
    /** @brief Utility routine to set Mandelstam s
     * @param[in] x_ Mandelstam s
     **/
    void set_s(const double x_) {s_int = x_;}
    
    /** @brief Utility routine to set kappa for Debye mass
     * @param[in] x_ kappa
     **/
    void set_kappa(const double x_) {kappa_int = x_;}
    
    /** @brief Utility routine to set formation time type for radiated gluon
     * @param[in] x_ formationTimeTyp
     **/
    void set_formationTimeTyp(const std::string & x_) {formationTimeTyp_int = x_;}
    
    /** @brief Utility routine to set matrix element which is used for 2->3
     * @param[in] x_ matrixElement23
     **/
    void set_matrixElement23(const std::string & x_) {matrixElement23_int = x_;}
    
    /** @brief Utility routine to set 2->2 part of matrix element which is used for 2->3
     * @param[in] x_ matrixElement23
     **/
    void set_matrixElement23_22qt(const bool x_) {matrixElement23_22qt_int = x_;}
    
    /** @brief Utility routine to set whether a counter term is applied which resembles the old prescription of Xu for masssless particles
     * @param[in] x_ md2_counter_term_in_I23
     **/
    void set_md2_counter_term_in_I23(const bool x_) {md2_counter_term_in_I23_int = x_;}
    
    /** @brief Utility routine to set factor for LPM effect to play around
     * @param[in] x_ fudge_factor_lpm_int
     **/
    void set_fudge_factor_lpm(const double x_) {fudge_factor_lpm_int = x_;}
    
    
  private:
    /** @brief debye mass squared scaled by 1/s (dimensionless) */
    double md2_wo_as_int;
    /** @brief mean free path scaled by sqrt(s) (dimensionless) */
    double lambda_scaled_int;
    /** @brief boost velocity from lab frame to CMS of 2->3 scattering */
    double beta_int;
    /** @brief cos(theta) where theta is the angle between the boost velocity vector and the axis of the CMS system */
    double cos_theta_int;
    /** @brief mass of particle one, scaled by sqrt(s) (dimensionless) */
    double m1;
    /** @brief mass of particle two, scaled by sqrt(s) (dimensionless) */
    double m2;
    /** @brief Mandelstam s */
    double s_int;
    /** @brief kappa for Debye mass */
    double kappa_int;
    /** @brief type of formation time of radiated gluon in 2->3 */
    std::string formationTimeTyp_int;
    /** @brief matrix element which is used for 2->3 heavy quark scattering */
    std::string matrixElement23_int;
    /** 
     * @brief Whether the gluon propagator in the 2->2 part of the radiative cross section is approximated by 1/qt^4. 
     *
     * If false the more exact expression 1/t^2 is employed. Applies
     * only to matrixElement23 = GBimproved. For GBorg always the
     * approximated propagator is used. 
     **/ 
    bool matrixElement23_22qt_int;
    /** @brief Whether a counter term is applied which resembles the old prescription of Xu for masssless particles */
    bool md2_counter_term_in_I23_int;
    /** @brief factor for LPM effect to play around */
    double fudge_factor_lpm_int;
};

/** 
 * @brief exception class for handling unexpected critical behaviour within integrand23
 **/
class eInt23_error : public std::runtime_error
{
  public:
    explicit eInt23_error(const std::string& what) : std::runtime_error(what) {};
    
    virtual ~eInt23_error() throw() {};
};


#endif
