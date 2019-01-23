//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/scattering23.h $
//$LastChangedDate: 2014-11-16 22:00:01 +0100 (So, 16. Nov 2014) $
//$LastChangedRevision: 1938 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file 
 * @brief Declarations for class integrand23.
 */

#ifndef SCATTERING23_H
#define SCATTERING23_H

#include <stdexcept>
#include <string>
#include "interpolation23.h"
#include "particleprototype.h"
#include "lorentz.h"

int rootsCubicPolynomial( const double a, const double b, const double c, const double d, double sol[3] );


/** @brief Class encapsulating the calculation of 2->3 scatterings */
class scattering23
{
public:
  // <<-------------------------------------------------
  // some wrapper for old routines to make new routines for heavy quarks with running coupling and more choices in formationTimeTyp, matrix element, etc. compatible with the old ones, try not to use anymore. work only with massless particles
  
  /** @brief Basic constructor. scattering23::setParameter MUST be used when using this constructor. DEPRECATED */
  scattering23( const interpolation23 * const theI23_arg );
  
  /** @brief Constructor setting all necessary parameters. 
   * Set default value of _Ng to -1 since it is currently not used in routines. If it is going to be used, do not set a default value and modify all programmes that do not compile anymore. DEPRECATED
   */
  scattering23( const interpolation23 * const theI23_arg, 
                const VectorXYZ & v, 
                const VectorEPxPyPz & P1_arg, 
                const VectorEPxPyPz & P2_arg, 
                const FLAVOR_TYPE F1_arg, 
                const FLAVOR_TYPE F2_arg,
                const double sqrtS_arg, 
                const double md2g_scaled_arg, 
                const double lambda_scaled_arg, 
                const double _alpha_s, 
                const int _Ng = -1 );
  
  /** @brief Sets the internal parameter needed for a specific particle doublet. 
   * Set default value of _Ng to -1 since it is currently not used in routines. If it is going to be used, do not set a default value and modify all programmes that do not compile anymore. DEPRECATED
   */
  double setParameter( const VectorXYZ & v, 
                       const VectorEPxPyPz & P1_arg, 
                       const VectorEPxPyPz & P2_arg, 
                       const FLAVOR_TYPE F1_arg, 
                       const FLAVOR_TYPE F2_arg, 
                       const double sqrtS_arg,
                       const double md2g_scaled_arg, 
                       const double lambda_scaled_arg, 
                       const double _alpha_s,
                       const int _Ng = -1 );
  // --------------------------------------------------->>




  /** @brief Basic constructor. scattering23::setParameter MUST be used when using this constructor. */
  scattering23( const interpolation23 * const theI23_massless_arg,
                const interpolation23 * const theI23_charm_m1_arg, 
                const interpolation23 * const theI23_charm_m2_arg, 
                const interpolation23 * const theI23_bottom_m1_arg, 
                const interpolation23 * const theI23_bottom_m2_arg );

  /** @brief Constructor setting all necessary parameters. 
   * Set default value of _Ng to -1 since it is currently not used in routines. If it is going to be used, do not set a default value and modify all programmes that do not compile anymore.
   */
  scattering23( const interpolation23 * const theI23_massless_arg, 
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
                const double K_light_parton_arg = 1.0, 
                const double K_heavy_quark_arg = 1.0, 
                const double kappa_light_parton_arg = 1.0, 
                const double kappa_heavy_quark_arg = 1.0, 
                const bool I23onlineIntegration_arg = false, 
                const std::string & formationTimeTyp_arg = "bamps_org", 
                const std::string & matrixElement23_arg = "GBimproved", 
                const bool md2_counter_term_in_I23_arg = true, 
                const double fudge_factor_lpm_arg = 1.0,
                const int _Ng = -1, 
                const bool matrixElement23_22qt_arg = false );



  /** @brief Destructor */
  ~scattering23();


  /** @brief Sets the internal parameter needed for a specific particle doublet. 
   * Set default value of _Ng to -1 since it is currently not used in routines. If it is going to be used, do not set a default value and modify all programmes that do not compile anymore.
   */
  double setParameter( const VectorXYZ & v, 
                       const VectorEPxPyPz & P1_arg, 
                       const VectorEPxPyPz & P2_arg, 
                       const FLAVOR_TYPE F1_arg, 
                       const FLAVOR_TYPE F2_arg, 
                       const double m1_arg, 
                       const double m2_arg,
                       const double sqrtS_arg, 
                       const double md2g_wo_as_scaled_arg, 
                       const double lambda_scaled_arg,
                       const double K_light_parton_arg = 1.0, 
                       const double K_heavy_quark_arg = 1.0, 
                       const double kappa_light_parton_arg = 1.0, 
                       const double kappa_heavy_quark_arg = 1.0, 
                       const bool I23onlineIntegration_arg = false, 
                       const std::string & formationTimeTyp_arg = "bamps_org", 
                       const std::string & matrixElement23_arg = "GBimproved", 
                       const bool md2_counter_term_in_I23_arg = true, 
                       const double fudge_factor_lpm_arg = 1.0,
                       const int _Ng = -1, 
                       const bool matrixElement23_22qt_arg = false );
    
  /** @brief Returns total cross section in 1/GeV^2. Wrapper */
  double getXSection23() const { int temp; return getXSection23( temp ); }

  /** @brief Returns total cross section in 1/GeV^2. */
  double getXSection23( int& initialStateIndex ) const;

  /** @brief For calculation of the tables. */
  static double getIntegral23_special( const double _md2_wo_as_scaled, 
                                       const double _lambda_scaled, 
                                       const double _beta, 
                                       const double _cos_theta, 
                                       const double _m1_scaled, 
                                       const double _m2_scaled, 
                                       const bool _I23onlineIntegration, 
                                       const double _s, 
                                       const double _kappa, 
                                       const std::string & _formationTimeTyp, 
                                       const std::string & _matrixElement23, 
                                       const bool _matrixElement23_22qt, 
                                       const bool _md2_counter_term_in_I23, 
                                       const double _fudge_factor_lpm, 
                                       bool & I23_result_trustable_arg );
    
  /** @brief Samples new momenta.*/ 
  int getMomenta23( double& pt1, double& pt3, double& y, double& phi, double& pz1, int& typ, FLAVOR_TYPE& F1arg, FLAVOR_TYPE& F2arg) const;
    
  /** @brief Sets new momenta.*/
  void setNewMomenta23(VectorEPxPyPz & P1, 
                       VectorEPxPyPz & P2,
                       VectorEPxPyPz & P3,
                       const VectorTXYZ & R1,
                       const VectorTXYZ & R2,
                       const double PT1, const double PT3,
                       const double y3, const double phi, const double PZ1);

  /** @brief Utility routine giving access to cos_theta. */
  double getCosTheta() const {return cos_theta;}
    
  /** @brief Get integration range for rapidity integration */
  static bool y_range_23( double& y_min, 
                          double& y_max, 
                          const double kt, 
                          const double M2, 
                          const double m1_2, 
                          const double m2_2, 
                          const double lambda_scaled_yrange, 
                          const double beta_yrange, 
                          const double cos_theta_yrange, 
                          const std::string & formationTimeTyp_yrange, 
                          const double fudge_factor_lpm_yrange );
    
  /** @brief Get the integrand */
  static double get_integrand_23( double& pz11, 
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
                                  const std::string & formationTimeTyp_int, 
                                  const std::string & matrixElement23_int, 
                                  const bool matrixElement23_22qt_int, 
                                  const bool md2_counter_term_in_I23_int, 
                                  const double fudge_factor_lpm_int );
    
  bool isRandomized() const { return ( randomizedConfiguration != 0 ); };
    
private:
  /** @brief Returns the integral over the 2->3 matrix element. Wrapper */
  double getIntegral23() const { int temp; return getIntegral23( temp ); }
  /** @brief Returns the integral over the 2->3 matrix element. */
  double getIntegral23( int& initialStateIndex ) const;  
    
  /** @brief Computes the maximal range in rapidity y */
  bool maximal_Y_range_23(double& ymin, double& ymax) const;
    
  void rotation( const VectorEPxPyPz & P, const VectorTXYZ & R1, const VectorTXYZ & R2, VectorEPxPyPz & PP, VectorEPxPyPz & TT ) const;
    
  /** @brief Velocity with which the computational cell moves (if any). */
  lorentz LL_cell;

  /** @brief Boost from lab frame to CMS. */
  lorentz LL_CM;

  /** @brief absolute value of beta of #LL_CM */
  double beta;
    
  /** @brief cos(theta) where theta is the angle between the boost velocity vector and the axis of the CMS system */
  double cos_theta;

  /** @brief angle between the boost velocity vector and the axis of the CMS system */
  double theta;
    
  /** @brief Original momentum of particle 1 in lab frame. */
  VectorEPxPyPz P1;
  /** @brief Original momentum of particle 2 in lab frame. */
  VectorEPxPyPz P2;
    
  /** @brief Momentum of particle 1 in rest frame of the computational cell. */
  VectorEPxPyPz P1cell;  
  /** @brief Momentum of particle 2 in rest frame of the computational cell. */
  VectorEPxPyPz P2cell;
    
  /** @brief Momentum of particle 1 in CM frame of the scattering. */
  VectorEPxPyPz P1cm;
  /** @brief Momentum of particle 2 in CM frame of the scattering. */
  VectorEPxPyPz P2cm;
    
  /** @brief flavor of particle 1 */
  FLAVOR_TYPE F1;
  /** @brief flavor of particle 2 */
  FLAVOR_TYPE F2;

  /** @brief mass of particle 1 */
  double m1;
  /** @brief mass of particle 2 */
  double m2;
    
  /** @brief Debye mass squared, scaled by 1/s */
  double md2g_wo_as_scaled;
  
  /** @brief Center of mass energy, sqrt(s) */
  double sqrtS;
  
  /** @brief Mean free path, scaled by sqrt(s) */
  double lambda_scaled;
    
  /** @brief number of gluons in given cell, possibly needed to rescale prefactors */
  int numberOfGluons;
    
  /** @brief flag that determines whether the ordering of the input vectors is randomized or not */
  bool randomizeInput;

  /** @brief the actual ordering of the input vectors that might be randomized (0 is default) */ 
  int randomizedConfiguration;
    
  /** @brief Pointer to an #interpolation23 object that is used for interpolating the integrated 2->3 matrix element */
  const interpolation23 * const theI23_massless;
  const interpolation23 * const theI23_charm_m1;
  const interpolation23 * const theI23_charm_m2;
  const interpolation23 * const theI23_bottom_m1;
  const interpolation23 * const theI23_bottom_m2;
  
  /** @brief true if I23 integration is performed online. */
  bool I23onlineIntegration;
  
  /** @brief type of formation time of radiated gluon in 2->3 */
  std::string formationTimeTyp;
  
  /** @brief matrix element which is used for 2->3 heavy quark scattering */
  std::string matrixElement23;
  
  /** @brief Whether the gluon propagator in the 2->2 part of the radiative cross section is approximated by 1/qt^4. If false the more exact expression 1/t^2 is employed. Applies only to matrixElement23 = GBimproved. For GBorg always the approximated propagator is used. */
  bool matrixElement23_22qt;
  
  /** @brief Whether a counter term is applied which resembles the old prescription of Xu for masssless particles */
  bool md2_counter_term_in_I23;
  
  /** @brief factor for LPM effect to play around 
   * this factor is used to play around with the cutoff
   * instead of k_t > gamma / lambda a modified cutoff k_t > X gamma / lambda
   * is used, where X is the fudge factor
   */
  double fudge_factor_lpm;
  
  /** @brief K factor */
  double K_factor;
  
  /** @brief kappa factor for Debye screening */
  double kappa;
};


/** @brief exception class for handling unexpected critical behaviour within 2->3 routines  */
class eScatt23_error : public std::runtime_error
{
public:
  explicit eScatt23_error(const std::string& what) : std::runtime_error(what) {};

  virtual ~eScatt23_error() throw() {};
};


#endif
