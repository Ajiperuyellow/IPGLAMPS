//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/thermal.h $
//$LastChangedDate: 2016-09-16 16:31:08 +0200 (Fr, 16. Sep 2016) $
//$LastChangedRevision: 2422 $
//$LastChangedBy: senzel $
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------


#ifndef THERMAL_H
#define THERMAL_H

#include "random.h"
#include "particleprototype.h"

enum QUANTUM_STATISTICS{ boltzmann, bose, fermi };

void initThermal( ParticlePrototype& particle, const double T, const double gluonFugacity = 1, const double quarkFugacity = 1, const bool quantum_statistics = false );
void initThermal( ParticlePrototype& particle, const FLAVOR_TYPE _flav, const double T, const double gluonFugacity = 1, const double quarkFugacity = 1, const bool quantum_statistics = false );
void sampleThermalMomenta( ParticlePrototype & particle, const double T, const FLAVOR_TYPE flav = gluon, const bool quantum_statistics = false );
void sampleThermalMomentaWithVelocity(ParticlePrototype & particle, const double T, const double vz);
void initProjectile( ParticlePrototype & particle, const double E, const FLAVOR_TYPE _F = gluon, const double angle_to_y_axis = M_PI/2.0 );

double getGluonDensity( const double T, const double fugacity = 1, const bool quantum_statistics = false, const double vz = 0.0 );
/** @brief total quark density = quark density + anti-quark density */
double getQuarkDensity( const double T, const double fugacity = 1, const bool quantum_statistics = false, const double vz = 0.0 ); 

double sample_energy( const double T, const QUANTUM_STATISTICS statistics, const double m = 0.0 );
double getEnergydistribution( const double E, const double T, const QUANTUM_STATISTICS statistics, const double m = 0.0 );
double getMeanThermalEnergy( const double T, const QUANTUM_STATISTICS statistics, const double m = 0.0 );
double boltzmann_rejection( const double T );


/** @brief exception class for handling unexpected critical behaviour within the configuration of the BAMPS run  */
class eThermal_error : public std::runtime_error
{
  public:
    explicit eThermal_error(const std::string& what) : std::runtime_error(what) {};
    
    virtual ~eThermal_error() throw() {};
};


#endif
