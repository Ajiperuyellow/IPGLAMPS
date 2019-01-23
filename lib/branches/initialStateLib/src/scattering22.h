//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/scattering22.h $
//$LastChangedDate: 2017-05-24 18:07:15 +0200 (Mi, 24. Mai 2017) $
//$LastChangedRevision: 2601 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file 
 * @brief definitions for the scattering22 class
 */

#ifndef SCATTERING22_H
#define SCATTERING22_H

#include <stdexcept>
#include "particleprototype.h"
#include "interpolation22.h"
#include "lorentz.h"
#include "hydroParticleType.h"

/**
 * @brief Provides routines and data structures for 2->2 scatterings.
 *
 * The class scattering22 encapsulates all routines and data structures needed for 2->2 scattering processes. It provides methods
 * to compute the cross section for a given particle pair (#getXSection22), sample the outgoing momenta (#getMomenta22) and to
 * asigne these momenta to the outgoing particles (#setNewMomenta22).
 *
 * A scattering22 object can be used for more than one particle doublet, just call #setParameter with the appropriate parameters
 * before the evaluation of each particle pair.
 */
class scattering22
{
  public:
    /** @brief Standard constructor for only light parton scattering. When used, setParameter needs to be called prior to other methods! */
    scattering22();
    
    /** @brief Standard constructor for heavy quark scattering. When used, setParameter needs to be called prior to other methods! */
    scattering22( const interpolation22 * const theI22 );
    
    /** @brief Constructor taking the same arguments as setParameter for light partons */
    scattering22( const VectorEPxPyPz & P1_arg, 
                  const VectorEPxPyPz & P2_arg, 
                  const FLAVOR_TYPE F1_arg, 
                  const FLAVOR_TYPE F2_arg, 
                  const double s_arg, 
                  const double md2g_arg, 
                  const double md2q_arg, 
                  const double Kfactor_light_arg = 1.0);
    
    /** @brief Constructor taking the same arguments as setParameter for heavy quarks, with optional temperature and Jpsi parameter */
    scattering22( const interpolation22 * const theI22,
                  const VectorEPxPyPz & P1_arg, 
                  const VectorEPxPyPz & P2_arg, 
                  const FLAVOR_TYPE F1_arg, 
                  const FLAVOR_TYPE F2_arg, 
                  const double M1_arg, 
                  const double M2_arg, 
                  const double s_arg, 
                  const double md2g_wo_as_arg, 
                  const double md2q_wo_as_arg, 
                  const double KggQQbar_arg, 
                  const double KgQgQ_arg, 
                  const double kappa_gQgQ_arg,
                  const bool isConstantCrossSecGQ_arg, 
                  const double constantCrossSecValueGQ_arg, 
                  const bool isotropicCrossSecGQ_arg,
                  const double Kfactor_light_arg = 1.0,
                  const double temperature_arg = 0.0, 
                  const double jpsi_dissociation_temperature_arg = 0.0, 
                  const bool isConstantCrossSecJpsi_arg = false, 
                  const double constantCrossSecValueJpsi_arg = 0.0 );
    
    /** @brief standard destructor */
    ~scattering22();
    
    /** @brief Sets the internal parameter needed for a specific particle pair for light partons. */    
    void setParameter( const VectorEPxPyPz & P1_arg, 
                       const VectorEPxPyPz & P2_arg,
                       const FLAVOR_TYPE F1_arg, 
                       const FLAVOR_TYPE F2_arg,
                       const double s_arg, 
                       const double vrel_arg,
                       const double md2g_arg, 
                       const double md2q_arg,
                       const double Kfactor_light_arg = 1.0);
    
    /** @brief Sets the internal parameter needed for a specific particle pair for heavy quarks, with optional temperature and Jpsi parameter  */    
    void setParameter( const VectorEPxPyPz & P1_arg, 
                       const VectorEPxPyPz & P2_arg,
                       const FLAVOR_TYPE F1_arg, 
                       const FLAVOR_TYPE F2_arg, 
                       const double M1_arg, 
                       const double M2_arg, 
                       const double s_arg, 
                       const double vrel_arg,
                       const double md2g_wo_as_arg, 
                       const double md2q_wo_as_arg, 
                       const double KggQQbar_arg, 
                       const double KgQgQ_arg, 
                       const double kappa_gQgQ_arg,
                       const bool isConstantCrossSecGQ_arg, 
                       const double constantCrossSecValueGQ_arg, 
                       const bool isotropicCrossSecGQ_arg,
                       const double Kfactor_light_arg = 1.0,
                       const double temperature_arg = 0.0, 
                       const double jpsi_dissociation_temperature_arg = 0.0, 
                       const bool isConstantCrossSecJpsi_arg = false, 
                       const double constantCrossSecValueJpsi_arg = 0.0 );
    
    
    /** @brief Returns the cross section for the given particle pair (set with #setParameter or in the constructor). Wrapper */
    double getXSection22() const { int temp; return getXSection22( temp ); }
    
    /** @brief Returns the cross section for the given particle pair (set with #setParameter or in the constructor). */
    double getXSection22( int& initialStateIndex ) const;
    
    /** @brief Returns the ELASTIC cross section for the given particle pair (set with #setParameter or in the constructor). Wrapper */
    double getXSectionElastic() const { int temp; return getXSectionElastic( temp ); }
    
    /** @brief Returns the ELASTIC cross section for the given particle pair (set with #setParameter or in the constructor). */
    double getXSectionElastic( int& initialStateIndex ) const; 
    
    /** @brief computes the transport cross section */
    double getTransportXSection22() const;
    
    /** @brief returns the mixture cross section */
  double getMixtureXSection22( hydroParticleType & theHydroParticleType, const std::vector<hydroParticleTypeProperties> & vecType, const int F1, const int F2, const double scalingFactor, std::string & info);
  
    /** @brief samples new momenta only for massless partons */
    void getMomenta22( double& t_hat, int& typ, FLAVOR_TYPE & F1arg, FLAVOR_TYPE & F2arg );
    
    /** @brief samples new momenta */
    void getMomentaAndMasses22(FLAVOR_TYPE& F1arg, FLAVOR_TYPE& F2arg, double& M1_out, double& M2_out, double& t_hat, int& typ);
    
//    /** @brief Jannis 22-isotropic routines  */
//    void getMomentaAndMasses22_isotropic( FLAVOR_TYPE& F1arg, FLAVOR_TYPE& F2arg, double& M1_out, double& M2_out, double& t_hat );
    
    /** @brief get the masses for the particles*/
    void getOnlyMasses22(FLAVOR_TYPE& F1_out, FLAVOR_TYPE& F2_out, double& M1_out, double& M2_out); 
    
    /** @brief samples new momenta for ELASTIC processes */
    void getMomentaElastic( double& PT2, int& typ ); 
    
//    /** @brief sample new momenta isotropically */
//    void getMomenta22_isotropic(double& PT2, int& typ);
    
    // /** @brief sample new momenta isotropically */
    void getMomentaAndMasses22_isotropic(FLAVOR_TYPE& F1arg, FLAVOR_TYPE& F2arg, double& M1_out, double& M2_out, double& t_hat );
    
    /** @brief sets new momenta for outgoing particles */
    void setNewMomenta22( VectorEPxPyPz & P1, 
                          VectorEPxPyPz & P2,
                          const VectorTXYZ & R1,
                          const VectorTXYZ & R2,
                          const double t_hat);
    double getS()
    {
      return s;     
    };
    
    double getVrel()
    {
      return vrel;
    };
    
  protected:
    /** @brief Auxiliary routine. Determines a direction for the transverse momtentum vector in the azimuthal plane. */
    void rotation( const VectorEPxPyPz & P, 
                   const VectorTXYZ & R1, 
                   const VectorTXYZ & R2, 
                   VectorEPxPyPz & PP, 
                   VectorEPxPyPz & TT ) const;
    
    /** @brief Boost from lab frame to CMS of particles.*/
    lorentz LL_CM;
    
    /** @brief Original momentum vector P1 in the lab frame. */
    VectorEPxPyPz P1;
    
    /** @brief Original momentum vector P2 in the lab frame. */ 
    VectorEPxPyPz P2;
    
    /** @brief Momentum vector P1 boosted to the CMS of the colliding particles.*/
    VectorEPxPyPz P1cm;
    
    /** @brief Momentum vector P2 boosted to the CMS of the colliding particles.*/
    VectorEPxPyPz P2cm;
    
    /** @brief flavor of particle 1 */
    FLAVOR_TYPE F1;
    
    /** @brief flavor of particle 2 */
    FLAVOR_TYPE F2;
    
    /** @brief masses of particles 1,2 -> 3,4 */
    double M1,M2,M3,M4;
    
    /** @brief gluonic debye mass (squared) divided by alpha_s */
    double md2_gluon_wo_as;
    
    /** @brief debye mass (squared) for quarks divided by alpha_s */
    double md2_quark_wo_as;
    
    /** @brief mandelstam variable s for the given particle pair */
    double s;

    /** @brief relative velocity for the given particle pair */
    double vrel;
    
    /** @brief smallest temperature of the 2 incoming heavy quarks. Is used for J/psi production, which is not allowed if temperature is above Td */
    double temperature;
    
    /** @brief Dissociation temperature of Jpsi. At higher temperature Jpsi production from QQbar is not possible. */
    double jpsi_dissociation_temperature;
    
    /** @brief Whether a constant cross section is employed for process g+Jpsi -> ccb. */
    bool isConstantCrossSecJpsi;
    
    /** @brief Value of constant cross section for process g+Jpsi -> ccb. */
    double constantCrossSecValueJpsi;
    
    /** @brief K factor for light parton processes */
    double Kfactor_light;
    
    /** @brief K factor for process g + g -> Q + Qbar */
    double KggQQbar;
    
    /** @brief K factor for process g + Q -> g + Q */
    double KgQgQ;
    
    /** @brief Kappa for Debye screening for process g + Q -> g + Q, usually 0.2 (Peshier,Gossiaux) */
    double kappa_gQgQ;
    
    /** @brief Whether an isotropic momentum sampling is employed for process g + Q -> g + Q */
    bool isotropicCrossSecGQ;
    
    /** @brief Whether a constant cross section is employed for process g + Q -> g + Q */
    bool isConstantCrossSecGQ;
    
    /** @brief Value of constant cross section for process g + Q -> g + Q */
    double constantCrossSecValueGQ;
    
    /** @brief Pointer to an #interpolation22 object that is used for interpolating the integrated cross sections for most processes */
    const interpolation22 * const theI22;    
};



/** @brief exception class for handling unexpected critical behaviour within 2->2 routines  */
class eScatt22_error : public std::runtime_error
{
  public:
    explicit eScatt22_error(const std::string& what) : std::runtime_error(what) {};
    
    virtual ~eScatt22_error() throw() {};
};


#endif
