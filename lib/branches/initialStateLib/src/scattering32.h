//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/scattering32.h $
//$LastChangedDate: 2014-11-16 22:00:01 +0100 (So, 16. Nov 2014) $
//$LastChangedRevision: 1938 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file 
 * @brief Declarations for the class scattering32.
 */


#ifndef SCATTERING32_H
#define SCATTERING32_H

#include <vector>
#include <stdexcept>
#include <string>

#include "integrand32.h"
#include "particleprototype.h"
#include "lorentz.h"


enum I32_COMPUTATION_TYPE { FAST, MONTE_CARLO_INTEGRATION };


/**
 * @brief Provides routines and data structures for 3->2 scatterings.
 *
 * The class scattering32 encapsulates all routines and data structures needed for 3->2 scattering processes. It provides methods to
 * calculate the integral over the matrix element I32, #getIntegral32_vegas (using numerical integration) or #getIntegral32_fast 
 * (using an estimation procedure). Sampling of new momenta after a 3->2 collision is provided by #getMomenta32 and #setNewMomenta32
 * then sets the particle momenta to the new values.
 *
 * A scattering32 object can be used for more than one particle triplet, just call #setParameter with the appropriate parameters
 * before the evaluation of each triplet.
 */
class scattering32
{
  public:
    /** @brief Standard constructor. When used, setParameter needs to be called prior to other methods! */
    scattering32();
    /** @brief Constructor taking the same arguments as setParamter 
     * Set default value of _Ng to -1 since it is currently not used in routines. If it is going to be used, do not set a default value and modify all programmes that do not compile anymore.
     */
    scattering32( const VectorXYZ & v,
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
                  const bool _md2_use_simplified_GB = false, 
                  const bool _matrixElement23_22qt = false, 
                  const double fudge_factor_lpm_arg = 1.0,
                  const int _Ng = -1 );
    /** @brief standard destructor */
    ~scattering32();
    
    /** @brief Sets the internal parameter needed for a specific particle triplet. 
     * Set default value of _Ng to -1 since it is currently not used in routines. If it is going to be used, do not set a default value and modify all programmes that do not compile anymore.
     */
    double setParameter( const VectorXYZ & v,
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
                         const bool _md2_use_simplified_GB = false, 
                         const bool _matrixElement23_22qt = false, 
                         const double fudge_factor_lpm_arg = 1.0, 
                         const int _Ng = -1 );
    
    /** @brief Returns bare I32 phase space integral with no prefactors */
    double getIntegral32_bare( int& initialStateIndex, const I32_COMPUTATION_TYPE _compType = FAST );
    /** @brief Returns bare I32 phase space integral with no prefactors. Wrapper. */
    double getIntegral32_bare( const I32_COMPUTATION_TYPE _compType = FAST ) { int temp; return getIntegral32_bare( temp, _compType ); }
    
    /** @brief Returns I32 phase space integral including all prefactors */
    double getIntegral32_withPrefactors( int& initialStateIndex, const I32_COMPUTATION_TYPE _compType = FAST );
    /** @brief Returns I32 phase space integral including all prefactors. Wrapper. */
    double getIntegral32_withPrefactors( const I32_COMPUTATION_TYPE _compType = FAST ) { int temp; return getIntegral32_withPrefactors( temp, _compType ); }
    
    /** @brief Samples new momenta. Wrapper method.*/ 
    int getMomenta32( double& u, double& phi, int& typ, FLAVOR_TYPE& F1arg, FLAVOR_TYPE& F2arg);   
    /** @brief Sets new momenta.*/
    void setNewMomenta32( VectorEPxPyPz & P1_arg, 
                          VectorEPxPyPz & P2_arg, 
                          const double u, const double phi) const;
    
    
    /** @brief Returns the number indicating the actual order of P1, P2, P3. */
  //    int getOrder() const;
        
    
    static int choice1_qq;
    static int choice2_qq;
    static int choice3_qq;
    static int events_qq;
    
    static int choice1_gg;
    static int choice2_gg;
    static int choice3_gg;
    static int events_gg;


  private:
    /** @brief Returns bare I32 using numerical integration. */
    double getIntegral32_bare_vegas( int& initialStateIndex );
    /** @brief Returns bare I32 using numerical integration. Wrapper.*/
    double getIntegral32_bare_vegas() { int temp; return getIntegral32_bare_vegas( temp ); }
    /** @brief Returns bare I32 using an estimation procedure. */
    double getIntegral32_bare_fast( int& initialStateIndex );
    /** @brief Returns bare I32 using an estimation procedure. Wrapper. */    
    double getIntegral32_bare_fast() { int temp; return getIntegral32_bare_fast( temp ); }
    
    /** @brief Transforms the boost vector #beta_vec to the reference frame given by P1_arg and P3_arg*/
    void rotate_beta( const vector4D & P1_arg, 
                      const vector4D & P3_arg, 
                      vector4D & beta_new);
    
    /** @brief Get new momenta via the rejection method. Called by wrapper getMomenta32(double& u, double& phi)*/
    int getMomenta32_rejection(double& u, double& phi) const;
    /** @brief Get new momenta via the metropolis algorithm. Called by wrapper getMomenta32(double& u, double& phi)*/
    int getMomenta32_metropolis(double& u_arg, double& phi_arg) const;
    
    /** @brief Returns the matrix element at given u and phi*/
    double getMatrixElement(const double u, const double phi) const;
    
    /** @brief Returns the estimated I32.*/
    double get_I32estimate() const;
    /** @brief Returns the ratio of the "real" I32 to the estimated I32.*/
    double get_I32estimate_ratio();
    
    /** @brief Utility method needed by #get_I32estimate.*/
    double H_simplifiedGB(const double a, const double u) const;
    /** @brief Utility method needed by #get_I32estimate.*/
    double H_actualGB(const double a, const double u) const;    
    
    /** @brief Utility method needed by #H_actualGB.*/
    double H1(const double a, const double u) const { return H_simplifiedGB( a, u ); }
    /** @brief Utility method needed by #H_actualGB.*/
    double H2(const double a, const double u) const;
    /** @brief Utility method needed by #H_actualGB.*/
    double H3(const double a, const double u) const;
    
    /** @brief The estimate of the integrand for the GB version where terms are simplified prior to applying the screening */
    double F_estimate_simplifiedGB( const double a2, const double u ) const;
    /** @brief The estimate of the integrand for the GB version where terms are not simplified prior to applying the screening */
    double F_estimate_actualGB( const double a2, const double u ) const;
        
    /** @brief Function object representing the integrand.
     *
     * integrand32 overloads the () operator, can be called like a normal function.
     */
  integrand32 theIntegrand;
    
    /** @brief Collective velocity of the cell. */
  lorentz LL_cell;

    /** @brief Boost velocity from cell rest frame to CMS of particle triplet.*/
  lorentz LL_CM;

  /** @brief absolute value of beta of #LL_CM */
    double beta_abs;

    /** @brief Boost vector #beta_vec rotated into the frame given by p1 and p3. */
  vector4D rotated_beta;

      /** @brief Original momentum vector P1 in the lab frame. */
  VectorEPxPyPz P1;
  /** @brief Original momentum vector P2 in the lab frame. */  
  VectorEPxPyPz P2;
  /** @brief Original momentum vector P3 in the lab frame. */
  VectorEPxPyPz P3;
    
  /** @brief Momentum vector P1 boosted to the cell frame.*/
  VectorEPxPyPz P1cell;
  /** @brief Momentum vector P2 boosted to the cell frame.*/  
  VectorEPxPyPz P2cell;
  /** @brief Momentum vector P3 boosted to the cell frame.*/
  VectorEPxPyPz P3cell;
    
  /** @brief Momentum vector P1 boosted to the CMS of the colliding particle triplet.*/
  VectorEPxPyPz P1cm;
  /** @brief Momentum vector P2 boosted to the CMS of the colliding particle triplet.*/
  VectorEPxPyPz P2cm;
  /** @brief Momentum vector P3 boosted to the CMS of the colliding particle triplet.*/
  VectorEPxPyPz P3cm;
    
    /** @brief flavor of particle 1 */
    FLAVOR_TYPE F1;
    /** @brief flavor of particle 2 */
    FLAVOR_TYPE F2;
    /** @brief flavor of particle 3 */
    FLAVOR_TYPE F3;
    
    /** @brief alpha_s of the considered process */
    double alpha;
    /** @brief Debye mass squared, scaled by 1/s */
    double md2g_scaled;
    /** @brief Center of mass energy, sqrt(s) */
    double sqrtS;
    /** @brief Mean free path, scaled by sqrt(s) */
    double lambda_scaled;
    
    /** @brief number of gluons in given cell, possibly needed to rescale prefactors */
    int numberOfGluons;
    
    /** @brief Energy of incoming particle 1.
     * The suffix "selected" is due to historic reasons. It is kept for the sake of better comparability with legacy versions of the code
     */
    double E1_selected;
    /** @brief Energy of incoming particle 3.
     * The suffix "selected" is due to historic reasons. It is kept for the sake of better comparability with legacy versions of the code
     */
    double E3_selected;
    /** @brief cos(gamma) for the pair (p1,p3). */
    double cos_gamma;
    /** @brief Indicates which of the gluons is going to be absorped */      
    unsigned short int absorpedGluon;  
    
    
    /** @brief Whether the version of GB is used that is simplified prior to applying the screening with md2
     * This is equivalent to the option md2_counter_term_in_I23 in scattering23. Setting it to true is equivalent to using the "old" version of GB (but including the (1-x) term).
     */
    bool md2_use_simplified_GB;
    
    /** @brief Whether the gluon propagator in the 2->2 part of the radiative cross section is approximated by 1/qt^4. If false the more exact expression 1/t^2 is employed. */
    bool matrixElement32_22qt;    
    
    /** @brief factor for LPM effect to play around 
     * this factor is used to play around with the cutoff
     * instead of k_t > gamma / lambda a modified cutoff k_t > X gamma / lambda
     * is used, where X is the fudge factor
     */
    double fudge_factor_lpm;
};



/** @brief exception class for handling unexpected critical behaviour within 2->3 routines  */
class eScatt32_error : public std::runtime_error
{
  public:
    explicit eScatt32_error(const std::string& what) : std::runtime_error(what) {};

    virtual ~eScatt32_error() throw() {};
};

#endif
