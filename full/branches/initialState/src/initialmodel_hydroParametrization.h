//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/initialmodel_hydroParametrization.h $
//$LastChangedDate: 2014-11-16 21:18:55 +0100 (日, 16 11月 2014) $
//$LastChangedRevision: 1935 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifndef INITIALMODEL_HYDROPARAMETRIZATION_H
#define INITIALMODEL_HYDROPARAMETRIZATION_H

#include <stdexcept>
#include <vector>
#include <string>

#include "initialmodel.h"
#include "configuration.h"
#include "particle.h"
#include "hydroParticleType.h"
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>


/**
 * @brief Class to provide the hydro parametrization initialization
 */
class initialModel_HydroParametrization : public initialModelWS
{
  public:
  initialModel_HydroParametrization( const config& _config, string &infoInitialParameters );
    ~initialModel_HydroParametrization() {};
    
    void populateParticleVector( std::vector<Particle>& _particles );
    
    double theNiemiParametrization();
    double nuclearThicknessFunctionA();
    double nuclearThicknessFunctionB();    
    
    //parameter for integration
    double intPar_x, intPar_y, intPar_z, intPar_pT, intPar_rap, intPar_NTF_z;
    
  protected:    

  private:    
  double integration_HIC_TotalNumber();  
  double integration_HIC_TotalEnergy();
  double integration_HIC_TotalDensityAtZeroPoint( const double x_0, const double x_max);
  void writeInfoInOutput( string &infoInitialParameters );
  void theNTFcalc();
  
  void sampleParticles( std::vector< Particle >& _particles );
  void initialisationNumberForEveryParticleSpecies(hydroParticleType& hp, std::vector<hydroParticleTypeProperties> vecType, std::vector<int> &numberParticles);
  void samplePositions( std::vector< Particle >& _particles, std::vector<hydroParticleTypeProperties> vecType, const int type, const int startValue, const int endValue, double &numAcceptReject);
  void sampleMomenta( std::vector< Particle >& _particles, std::vector<hydroParticleTypeProperties> vecType, const int type, const int startValue, const int endValue, double &numAcceptReject);
  void sampleProperties( std::vector< Particle >& _particles, std::vector<hydroParticleTypeProperties> vecType, const int type, const int startValue, const int endValue, double &numAcceptReject);
  void setLightQuarkFlavor(std::vector<Particle>& _particles, const int j);
  
  double theNiemiParametrization_pT();
  double theNiemiParametrization_rap();
  double theNiemiParametrization_z();
  void theNiemiParametrization_TA_TB(double &TA, double &TB); 
  
  double sampling_TA_TB(double &numAcceptReject);
  double sampling_z(double &numAcceptReject);
  double sampling_pT(double &numAcceptReject, const double maximumPT);
  double sampling_rap(double &numAcceptReject); 
  double findMaximumPT();
  
  double integration_NuclearThicknessFunctionA();
  double integration_NuclearThicknessFunctionB();  
  double nuclearThicknessFunction_woodSaxon(const double nuc, const double s_x, const double s_y);
  
  int nIntegration,nIntegrationNTF;

  //Nnuclear Thickness functions
  std::vector<std::vector<double> > T_A;  
  std::vector<std::vector<double> > T_B;    

  double theNTFStepSize;
  int theNTFColumnsX,theNTFColumnsY;
  
  //variables
  int testparticles;
  double numberOfParticlesToGenerate;
  double boxLengthX,boxLengthY,boxLengthZ;
  double x_0,y_0,z_0,pT_0,rap_0,x_max,y_max,z_max,pT_max,rap_max,NTF_z0,NTF_zMax;
  double p_K,p_Q,p_n,p_m,p_sigRap,p_sigZ,p_A,p_B,p_bx,p_by,p_ksi,p_rho0;    

};



/** @brief exception class for handling unexpected critical behaviour within generation of hydro parametrization initial distributions  */
class eHydroParametrization_error : public std::runtime_error
{
  public:
    explicit eHydroParametrization_error(const std::string& what) : std::runtime_error(what) {};
    
    virtual ~eHydroParametrization_error() throw() {};
};




#endif // INITIALMODEL_PYTHIA_H
