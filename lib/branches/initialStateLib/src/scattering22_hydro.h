//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/scattering22_hydro.h $
//$LastChangedDate: 2017-11-12 20:59:44 +0100 (So, 12. Nov 2017) $
//$LastChangedRevision: 2641 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file 
 * @brief definitions for the scattering22_hydro class
 */

#ifndef SCATTERING22HYDRO_H
#define SCATTERING22HYDRO_H

#include <stdexcept>
#include "hydroParticleType.h"

/**
 * @brief Provides routines for 2->2 scatterings and hydro conditions, but not pQCD.
 *
 * The class scattering22 encapsulates all routines and data
 * structures needed for 2->2 scattering processes in the
 * hydrodynamics framework. 
 */
class scattering22_hydro
{
  public:
    /** @brief Standard constructor */
    scattering22_hydro();
    /** @brief Standard destructor */
    ~scattering22_hydro();
    
    //routines for scatt22 and isotropic
    //----------------------------------------
    /** @brief returns a special variable cross section for one component hydro reproducing multi-component ( only for experimental case, with Andrej El */
    double getVarTimeDepCS( const double time );
    
    /** @brief returns the mixture cross section */
    double getMixtureXSection22( hydroParticleType & theHydroParticleType, const std::vector<hydroParticleTypeProperties> & vecType,
                                 const int F1, const int F2, const double scalingFactor, std::string & info);

    /** @brief returns the cross section for consant mean free path */
    double getConstMFP_XSection( hydroParticleType & theHydroParticleType, const double mfp, const double T00, const double T11, const double T22, const double T33, const double T10, const double T20, const double T30, const double T21, const double T31, const double T32, const double N0, const double N1, const double N2, const double N3, const int N );
    
    
    /** @brief returns the cross section for consant eta over s */
    double getConstEtaOverS_XSection( hydroParticleType & theHydroParticleType, const std::vector<hydroParticleTypeProperties> & vecType, const double etaOverS,
                                      const std::vector<double> & T00, const std::vector<double> & T11, const std::vector<double> & T22,
                                      const std::vector<double> & T33, const std::vector<double> & T10, const std::vector<double> & T20,
                                      const std::vector<double> & T30, const std::vector<double> & T21, const std::vector<double> & T31,
                                      const std::vector<double> & T32, const std::vector<double> & N0, const std::vector<double> & N1,
                                      const std::vector<double> & N2, const std::vector<double> & N3, const std::vector<int> & NN , double & extractedT, double & extractedN, double & extractedFug );
    
    /** 
     * @brief Calculates the scaling factor according to lambda*T^2, which gets into the cross section
     */        
    double getScalingFactorLambdaT2( hydroParticleType & theHydroParticleType, const std::vector<hydroParticleTypeProperties> & vecType,
                                     const std::vector<double> & T00, const std::vector<double> & T11, const std::vector<double> & T22,
                                     const std::vector<double> & T33, const std::vector<double> & T10, const std::vector<double> & T20,
                                     const std::vector<double> & T30, const std::vector<double> & T21, const std::vector<double> & T31,
                                     const std::vector<double> & T32, const std::vector<double> & N0, const std::vector<double> & N1,
                                     const std::vector<double> & N2, const std::vector<double> & N3, const std::vector<int> & NN );          
};


/** @brief exception class for handling unexpected critical behaviour
    within 2->2 hydro routines  */
class eScatt22Hydro_error : public std::runtime_error
{
  public:
    explicit eScatt22Hydro_error(const std::string& what) : std::runtime_error(what) {};
    
    virtual ~eScatt22Hydro_error() throw() {};
};


#endif
