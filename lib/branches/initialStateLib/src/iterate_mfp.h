//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL: file:///home/fochler/svnrepos/cascade_misc/createMFPtables/trunk/src/iteration.h $
//$LastChangedDate: 2011-01-08 00:05:01 +0100 (Sat, 08 Jan 2011) $
//$LastChangedRevision: 264 $
//$LastChangedBy
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------


#ifndef ITER_H
#define ITER_H

#include "interpolation23.h"
#include "interpolation22.h"
#include "configurationbase.h"
#include "ratesmanager.h"

double iterate_mfp_bisection( const double T, const interpolation23& theI23_massless, const interpolation23& theI23_charm_m1, const interpolation23& theI23_charm_m2, const interpolation23& theI23_bottom_m1, const interpolation23& theI23_bottom_m2, const interpolation22& theI22, ParticlePrototype & proj,
                              const configBase& theConfig, std::vector<ratesManager>& ratesArray, ratesManager& iteratedRates, const double gluonFugacity = 1.0, const double quarkFugacity = 1.0, const bool quantumStatistics = false, int samplings_per_iter = 0 );

/** @brief wrapper function for iterate mfp bi-section function if iterated rates should not be returned */
inline double iterate_mfp_bisection( const double T, interpolation23& theI23_massless, interpolation23& theI23_charm_m1, interpolation23& theI23_charm_m2, interpolation23& theI23_bottom_m1, interpolation23& theI23_bottom_m2, interpolation22& theI22, ParticlePrototype & proj,
                const configBase& theConfig, const double gluonFugacity = 1.0, const double quarkFugacity = 1.0, const bool quantumStatistics = false, int samplings_per_iter = 0 )
{
  std::vector<ratesManager> ratesArray;
  ratesManager iteratedRates;
  
  return iterate_mfp_bisection( T, theI23_massless, theI23_charm_m1, theI23_charm_m2, theI23_bottom_m1, theI23_bottom_m2, theI22, proj, theConfig, ratesArray, iteratedRates, gluonFugacity, quarkFugacity, quantumStatistics, samplings_per_iter );
}


double iterate_mfp( const double T, const interpolation23& theI23_massless, const interpolation23& theI23_charm_m1, const interpolation23& theI23_charm_m2, const interpolation23& theI23_bottom_m1, const interpolation23& theI23_bottom_m2, const interpolation22& theI22, ParticlePrototype &proj, const double lambdaStart,
                    const configBase& theConfig, std::vector<ratesManager>& ratesArray, ratesManager& iteratedRates, const double gluonFugacity = 1.0, const double quarkFugacity = 1.0, const bool quantumStatistics = false, const int samplings_per_iter = 400 );

/** @brief wrapper function for old iterate mfp function if iterated rates should not be returned */
inline double iterate_mfp( const double T, const interpolation23& theI23_massless, const interpolation23& theI23_charm_m1, const interpolation23& theI23_charm_m2, const interpolation23& theI23_bottom_m1, const interpolation23& theI23_bottom_m2, const interpolation22& theI22, ParticlePrototype &proj, const double lambdaStart,
                    const configBase& theConfig, const double gluonFugacity = 1.0, const double quarkFugacity = 1.0, const bool quantumStatistics = false, const int samplings_per_iter = 400  )
{
  std::vector<ratesManager> ratesArray;
  ratesManager iteratedRates;
  
  return iterate_mfp( T, theI23_massless, theI23_charm_m1, theI23_charm_m2, theI23_bottom_m1, theI23_bottom_m2, theI22, proj, lambdaStart, theConfig, ratesArray, iteratedRates, gluonFugacity, quarkFugacity, quantumStatistics, samplings_per_iter );
}

double iterate_thermal_mfp( const double T, const FLAVOR_TYPE flav, const interpolation23& theI23_massless, const interpolation23& theI23_charm_m1, const interpolation23& theI23_charm_m2, const interpolation23& theI23_bottom_m1, const interpolation23& theI23_bottom_m2, const interpolation22& theI22, const double lambdaStart, const configBase& theConfig, std::vector< ratesManager >& ratesArray, ratesManager& iteratedRates, const double gluonFugacity = 1.0, const double quarkFugacity = 1.0, const bool quantumStatistics = false, const int samplings_per_iter = 400 );

/** @brief wrapper function for iterate thermal mfp function if iterated rates should not be returned */
inline double iterate_thermal_mfp( const double T, const FLAVOR_TYPE flav, const interpolation23& theI23_massless, const interpolation23& theI23_charm_m1, const interpolation23& theI23_charm_m2, const interpolation23& theI23_bottom_m1, const interpolation23& theI23_bottom_m2, const interpolation22& theI22, const double lambdaStart,
                    const configBase& theConfig, const double gluonFugacity = 1.0, const double quarkFugacity = 1.0, const bool quantumStatistics = false, const int samplings_per_iter = 400 )
{
  std::vector<ratesManager> ratesArray;
  ratesManager iteratedRates;
  
  return iterate_thermal_mfp( T, flav, theI23_massless, theI23_charm_m1, theI23_charm_m2, theI23_bottom_m1, theI23_bottom_m2, theI22, lambdaStart, theConfig, ratesArray, iteratedRates, gluonFugacity, quarkFugacity, quantumStatistics, samplings_per_iter );
}

double get22mfp( const double T, interpolation22& theI22, const configBase& theConfig, const int nRuns22, const double gluonFugacity = 1.0, const double quarkFugacity = 1.0, const bool quantumStatistics = false );

#endif
