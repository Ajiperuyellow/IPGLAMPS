//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/hydroParticleType.h $
//$LastChangedDate: 2017-11-07 17:10:49 +0100 (Di, 07. Nov 2017) $
//$LastChangedRevision: 2639 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


#ifndef HYDROPARTICLETYPE_H
#define HYDROPARTICLETYPE_H

#include <math.h>
#include <string>
#include <vector>
#include <gsl/gsl_sf_bessel.h>
#include <complex>

#include "FPT_compare.h"

/**
 * @brief Class to define particle classes and store its properties
 **/
class hydroParticleTypeProperties
{
public:
  /**
   * @brief constructor
   **/
  hydroParticleTypeProperties( const double _degeneracy, 
                               const double _mass, 
                               const bool _massless, 
                               const int _flavStart, 
                               const int _flavEnd,
                               const std::string & _name) :
    degeneracyFactor( _degeneracy ),
    mass( _mass ),
    masslessON( _massless ),
    flavorTypeStart( _flavStart ),
    flavorTypeEnd( _flavEnd ),
    nameOfParticle( _name )
  { };
    
  double degeneracyFactor; ///< the degeneracy
  double mass; ///< the mass 
  bool masslessON; ///< treat it as massless
  int flavorTypeStart; ///< the first flavour
  int flavorTypeEnd; ///< the last flavour
  std::string nameOfParticle; ///< the name of the class
    
private:
};  

/**
 * @brief class to handle hydrodynamical calculations
 **/
class hydroParticleType
{
public:
  /** @brief Constructor */
  hydroParticleType() {};

  /** @brief Destructor */
  ~hydroParticleType() {};  
    
  //-----------------------------------------//  
  /**
   * @brief The hydro particle vector with informations for all
   * particle species used in the simulation 
   *
   * This first vector is the general particle info - should be
   * activated anyway neede for some routines, i.e. cross section
   *
   * You may e.g. set it via
   * @code
   * hydroParticleTypeVectorCommon = 
   * create_vector<hydroParticleTypeProperties>
   * ( hydro_gluon )( hydro_lightQuarks );
   * @endcode
   */    
  static std::vector<hydroParticleTypeProperties> hydroParticleTypeVectorCommon;
    
  /**
   * @brief The hydro particle vector with informations for all
   * particle species in the initialisation 
   *
   * This second vector is needed for particle initialization,
   * especially important for the box.
   *
   * @see hydroParticleTypeVectorCommon 
   */    
  static std::vector<hydroParticleTypeProperties> hydroParticleTypeVectorInit;  
    
  /**
   * @brief The hydro particle vector with informations for all
   * particle species used for the final analysis 
   *
   * This third vector is for the analysis in the end
   *
   * @see hydroParticleTypeVectorCommon 
   */    
  static std::vector<hydroParticleTypeProperties> hydroParticleTypeVectorAnalysis;  
  //-----------------------------------------// 

  /** @brief properties of the class 'gluon' */
  static const hydroParticleTypeProperties hydro_gluon;
  /** @brief properties of the class 'up' */
  static const hydroParticleTypeProperties hydro_up;
  /** @brief properties of the class 'antiup' */
  static const hydroParticleTypeProperties hydro_antiup;
  /** @brief properties of the class 'down' */
  static const hydroParticleTypeProperties hydro_down;
  /** @brief properties of the class 'antidown' */
  static const hydroParticleTypeProperties hydro_antidown;
  /** @brief properties of the class 'strange' */
  static const hydroParticleTypeProperties hydro_strange;
  /** @brief properties of the class 'antistrange' */
  static const hydroParticleTypeProperties hydro_antistrange;
  /** @brief properties of the class 'charm' */
  static const hydroParticleTypeProperties hydro_charm;
  /** @brief properties of the class 'anticharm' */
  static const hydroParticleTypeProperties hydro_anticharm;
  /** @brief properties of the class 'lightQuarks' */
  static const hydroParticleTypeProperties hydro_lightQuarks;
  /** @brief properties of the class 'charmQuarks' */
  static const hydroParticleTypeProperties hydro_charmQuarks;
  /** @brief properties of the class 'massA' */
  static const hydroParticleTypeProperties hydro_massA;
  /** @brief properties of the class 'massB' */
  static const hydroParticleTypeProperties hydro_massB;
  /** @brief properties of the class 'masslessC' */
  static const hydroParticleTypeProperties hydro_masslessC;
  /** @brief properties of the class 'masslessD' */
  static const hydroParticleTypeProperties hydro_masslessD;
  
    
  //-----------------------------------------//  
  /**
   * @brief Gives you the according particle type from the hydro particle vector. 
   *
   * As input you give the particle flavor. According to this you
   * get the related "new" particle flavor 
   */      
  int getParticleType(const int particleFlavor, std::vector<hydroParticleTypeProperties> vecType);
  //-----------------------------------------//
    
  //-----------------------------------------//  
  /**
   * @brief Gives the equilibrium temperature for a given energy density and particle density of a particle species.
   */  
  double chooseEqTemperature( const double eDensity, const double nDensity, std::vector<hydroParticleTypeProperties> vecType, const int type, bool & solution);
    
  /**
   * @brief Gives the equilibrium energy density for a given temperature and fugacity of a particle species.
   */  
  double chooseEqEnergyDensity( const double T, const double fug, std::vector<hydroParticleTypeProperties> vecType, const int type);
    
  /**
   * @brief Gives the equilibrium entropy density for a given temperature and fugacity of a particle species.
   */  
  double chooseEqEntropyDensity( const double T, const double fug, std::vector<hydroParticleTypeProperties> vecType, const int type);
    
  /**
   * @brief Gives the equilibrium particle density for a given temperature and fugacity of a particle species.
   */  
  double chooseEqParticleDensity( const double T, const double fug, std::vector<hydroParticleTypeProperties> vecType, const int type);
    
  /**
   * @brief Gives the equilibrium pressure for a given temperature and fugacity of a particle species.
   */    
  double chooseEqPressure( const double T, const double fug, std::vector<hydroParticleTypeProperties> vecType, const int type);
  //-----------------------------------------//      
    
  //-----------------------------------------//
  /**
   * @brief Gives the average temperature or fugacity of a multicomponent system. 
   *
   * Please be patient: 
   * You give vectors with the information of the values of each
   * seperate particle species, but the last component of the vector
   * has to be reserved for the average value. Also, please be sure
   * which particle vector information you give. This is needed for
   * the fugacity.  
   */    
  void lanEck_getAverageTemperatureAndFugacity( const std::vector<hydroParticleTypeProperties> & vecType, std::vector<double>& T, std::vector<double>& fug, const std::vector<double> & n );  
  //-----------------------------------------//  
    
  //-----------------------------------------// 
  /**
   * @brief Gives the fugacity in the Eckart or Landau frame. 
   * 
   * As input you give the particle density and temperature. Because
   * the fugacity depends on the given degeneracy factor, you have
   * to give the correct particle vector information 
   */     
  double lanEck_fugacity(const double nDensity, const double temperature, std::vector<hydroParticleTypeProperties> vecType, const int type);
    
  /**
   * @brief Gives the bulk pressure in the Eckart or Landau frame. 
   *
   * As input you give the isotropic pressure and equilibrium pressure.
   */         
  double lanEck_bulkPressure(const double isoPressure, const double eqPressure);
    
  /**
   * @brief Gives the energy density in the Landau frame. 
   *
   * As input you give components of the energy momentum tensor.
   */     
  double landau_energyDensity(const double T00, const double T11, const double T22, const double T33, const double T10, const double T20, const double T30, const double T21, const double T31, const double T32, bool & solution);
    
  /**
   * @brief Gives the velocities in the Landau frame. 
   * 
   * As input you give components of the energy momentum tensor and
   * the energy density in the Landau frame. 
   */  
  void landau_velocity(const double , const double T11,const double T22, const double T33, const double T10, const double T20, const double T30, const double T21, const double T31, const double T32, const double energyDensity, double& vx, double& vy, double& vz);
    
  /**
   * @brief Gives the velocities in the Eckart frame. 
   *
   * As input you give components of the particle 4-vector.
   */       
  void eckart_velocity(const double N0, const double N1, const double N2, const double N3, double& vx, double& vy, double& vz);
    
  /**
   * @brief Gives the energy density in the Eckart frame. 
   *
   * As input you give components of the energy momentum tensor and
   * the velocities in the Eckart frame. 
   */       
  double eckart_energyDensity(const double T00, const double T11, const double T22, const double T33, const double T10, const double T20,const double T30, const double T21, const double T31, const double T32, const double vx, const double vy, const double vz);
    
  /**
   * @brief Gives the particle density in Eckart or Landau frame. 
   * 
   * As input you give components of the particle 4-vector and the
   * velocities in the according frame. 
   */      
  double lanEck_particleDensity( const double N0, const double N1, const double N2, const double N3, const double vx, const double vy,const double vz);

  /**
   * @brief Gives the isotropic pressure in Eckart or Landau frame. 
   *
   * components of the energy momentum tensor and the velocities in the according frame.
   */     
  double lanEck_isotropicPressure( const double T00, const double T11, const double T22, const double T33, const double T10, const double T20, const double T30, const double T21, const double T31, const double T32, const double v_x, const double v_y, const double v_z);  
  
  /**
   * @brief Gives the components of the shear stress tensor in Eckart or Landau frame. 
   * 
   * As the input you give components of the energy momentum tensor
   * and the velocities in the according frame. 
   */        
  void lanEck_shearStress( const double T00, const double T11, const double T22, const double T33, const double T10, const double T20, const double T30, const double T21, const double T31, const double T32, const double v_x, const double v_y, const double v_z, double & PI00, double & PI11, double & PI22, double & PI33, double & PI10, double & PI20, double & PI30, double & PI21, double & PI31, double & PI32);

  /**
   * @brief Gives the components of the deltamunu in Eckart or Landau frame. 
   *
   * As input you give the velocities in the according frame.
   */     
  void lanEck_deltamunu(const double  v_x, const double  v_y, const double v_z, double & D_0_0, double & D$0_0, double & D_0$0, double & D$0$0, double & D_1_1, double & D$1$1, double & D_2_2, double & D$2$2, double & D_3_3, double & D$3$3, double & D$1_1, double & D$2_2, double & D$3_3, double & D$1_0, double & D$1$0, double & D$0$1, double & D$2_0, double & D$2$0, double & D$0$2, double & D$3_0, double & D$3$0, double & D$0$3, double & D$0_1, double & D_1_0, double & D_0_1, double & D$0_2, double & D_2_0, double & D_0_2, double & D$0_3, double & D_3_0, double & D_0_3, double & D_2_1, double & D$2$1, double & D_1_2, double & D$1$2, double & D_3_1, double & D$3$1, double & D_1_3, double & D$1$3, double & D_3_2, double & D$3$2, double & D_2_3, double & D$2$3, double & D$2_1, double & D$1_2, double & D$3_1, double & D$1_3, double &  D$3_2, double & D$2_3 );
  
  /**
   * @brief Gives the gamma factor in Eckart or Landau frame. 
   *
   * As input you give the velocity in the according frame.
   */    
  double lanEck_gamma(const double vx, const double vy, const double vz);
  
  /**
   * @brief Gives the absolute velocity in Eckart or Landau frame. 
   *
   * As input you give the velocities in the according frame.
   */   
  double lanEck_velocityAbsolute(const double vx, const double vy, const double vz);
  //-----------------------------------------//  

  /**
   * @brief Gives the components of the 4-velocity in Eckart or Landau frame. 
   *
   * As the input you give the velocities in the according frame.
   */       
  void lanEck_fourVelocity(const double vx, const double vy, const double vz, double& u$0, double& u_0, double& u$1, double& u_1, double& u$2, double& u_2, double& u$3, double& u_3 );  

  /**
   * @brief Gives the energy momentum flow in Eckart or Landau frame. 
   *
   * As the input you give components of the energy momentum tensor
   * and the velocities in the according frame. 
   */       
  void lanEck_energyMomentumFlow(const double T00, const double T11, const double T22, const double T33, const double T10, const double T20, const double T30, const double T21, const double T31, const double T32, const double v_x, const double v_y, const double v_z, double& W0, double& W1, double& W2, double& W3);

  /**
   * @brief Gives the particle flow in Eckart or Landau frame. 
   *
   * As input you give components of the particle 4-vector, the
   * velocities in the according frame. 
   */         
  void lanEck_particleFlow(const double N0, const double N1, const double N2, const double N3, const double v_x, const double v_y, const double v_z, double& q0, double& q1, double& q2, double& q3);
 
  /**
   * @brief Gives the heat flow in Eckart or Landau frame. 
   *
   * As input you give the energy momentum flow and particle flow,
   * the energy density, particle density and pressure. 
   */           
  void lanEck_heatFlow(const double W0, const double W1, const double W2, const double W3, const double V0, const double V1, const double V2, const double V3, const double energyDensity, const double pressure, const double particleDensity, double& q0, double& q1, double& q2, double& q3);
    
  /**
   * @brief Gives the PImunuPImunu in Eckart or Landau frame. 
   * 
   * As input you give the shear stress tensor components and the
   * velocities in the according frame. 
   */   
  double lanEck_PImunuPImunu( const double PI00, const double PI11, const double PI22, const double PI33, const double PI10, const double PI20, const double PI30, const double PI21, const double PI31, const double PI32, const double vx, const double vy, const double vz);  
    
  /**
   * @brief Gives all observables in Landau frame. 
   *
   * This is an example for the way to calculate the hydro
   * observables. The arrays contain all seperated species + the
   * whole system as the last one. Also, you have to give explicitly
   * the particle type vector for additional information. 
   */     
  void giveHydroObservablesInLandauFrame( const std::vector<hydroParticleTypeProperties> & vecType,
                                          const std::vector<double> & array_T00, const std::vector<double> & array_T11, const std::vector<double> & array_T22,
                                          const std::vector<double> & array_T33, const std::vector<double> & array_T10, const std::vector<double> & array_T20,
                                          const std::vector<double> & array_T30, const std::vector<double> & array_T21, const std::vector<double> & array_T31,
                                          const std::vector<double> & array_T32, const std::vector<double> & array_N0, const std::vector<double> & array_N1,
                                          const std::vector<double> & array_N2, const std::vector<double> & array_N3,
                                          std::vector<double> & array_energyDensity, std::vector<double> & array_particleDensity, std::vector<double> & array_entropyDensity,
                                          std::vector<double> & array_temperature,   std::vector<double> & array_fugacity,        std::vector<double> & array_eq_particleDensity,
                                          std::vector<double> & array_isoPressure,   std::vector<double> & array_eqPressure,      std::vector<double> & array_bulkPressure,
                                          std::vector<double> & array_vx,            std::vector<double> & array_vy,              std::vector<double> & array_vz,
                                          std::vector<double> & array_vv,            std::vector<double> & array_gamma,
                                          std::vector<double> & array_PI00, std::vector<double> & array_PI11, std::vector<double> & array_PI22,
                                          std::vector<double> & array_PI33, std::vector<double> & array_PI10, std::vector<double> & array_PI20,
                                          std::vector<double> & array_PI30, std::vector<double> & array_PI21, std::vector<double> & array_PI31,
                                          std::vector<double> & array_PI32, std::vector<double> & array_PImunuPImunu,
                                          std::vector<double> & array_q0, std::vector<double> & array_q1,
                                          std::vector<double> & array_q2, std::vector<double> & array_q3,
                                          std::vector<double> & array_W0, std::vector<double> & array_W1,
                                          std::vector<double> & array_W2, std::vector<double> & array_W3,
                                          std::vector<double> & array_V0, std::vector<double> & array_V1,
                                          std::vector<double> & array_V2, std::vector<double> & array_V3);
    
  /**
   * @brief Gives all observables in Eckart frame. 
   *
   * This is an example for the way to calculate the hydro
   * observables. The arrays containt all seperated species + the
   * whole system as the last one. Also, you have to give explicitly
   * the particle type vector for additional information. 
   */     
  void giveHydroObservablesInEckartFrame( const std::vector<hydroParticleTypeProperties> & vecType,
                                          const std::vector<double> & array_T00, const std::vector<double> & array_T11, const std::vector<double> & array_T22,
                                          const std::vector<double> & array_T33, const std::vector<double> & array_T10, const std::vector<double> & array_T20,
                                          const std::vector<double> & array_T30, const std::vector<double> & array_T21, const std::vector<double> & array_T31,
                                          const std::vector<double> & array_T32, const std::vector<double> & array_N0, const std::vector<double> & array_N1,
                                          const std::vector<double> & array_N2, const std::vector<double> & array_N3,
                                          std::vector<double> & array_energyDensity, std::vector<double> & array_particleDensity, std::vector<double> & array_entropyDensity,
                                          std::vector<double> & array_temperature,   std::vector<double> & array_fugacity,        std::vector<double> & array_eq_particleDensity,
                                          std::vector<double> & array_isoPressure,   std::vector<double> & array_eqPressure,      std::vector<double> & array_bulkPressure,
                                          std::vector<double> & array_vx,            std::vector<double> & array_vy,              std::vector<double> & array_vz,
                                          std::vector<double> & array_vv,            std::vector<double> & array_gamma,
                                          std::vector<double> & array_PI00, std::vector<double> & array_PI11, std::vector<double> & array_PI22,
                                          std::vector<double> & array_PI33, std::vector<double> & array_PI10, std::vector<double> & array_PI20,
                                          std::vector<double> & array_PI30, std::vector<double> & array_PI21, std::vector<double> & array_PI31,
                                          std::vector<double> & array_PI32, std::vector<double> & array_PImunuPImunu,
                                          std::vector<double> & array_q0, std::vector<double> & array_q1,
                                          std::vector<double> & array_q2, std::vector<double> & array_q3,
                                          std::vector<double> & array_W0, std::vector<double> & array_W1,
                                          std::vector<double> & array_W2, std::vector<double> & array_W3,
                                          std::vector<double> & array_V0, std::vector<double> & array_V1,
                                          std::vector<double> & array_V2, std::vector<double> & array_V3);  
private:
};


#endif
