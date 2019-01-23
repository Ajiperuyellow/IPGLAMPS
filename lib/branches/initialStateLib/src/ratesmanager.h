//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/ratesmanager.h $
//$LastChangedDate: 2016-09-13 17:53:26 +0200 (Di, 13. Sep 2016) $
//$LastChangedRevision: 2407 $
//$LastChangedBy: senzel $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


#ifndef RATESMANAGER_H
#define RATESMANAGER_H


#include <vector>

#include "particleprototype.h"
#include "interactiontype.h"
#include "unit_enum.h"

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @brief exception class for handling unexpected behaviour when managing interaction rates */
class eRatesManager_error : public std::runtime_error
{
public:
  explicit eRatesManager_error( const std::string& what ) : std::runtime_error( what ) {};

  virtual ~eRatesManager_error() throw() {};
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


/**
 * @brief Class encapsulating the management of scattering rates
 *
 */
class ratesManager
{
public:
  /** @brief Constructor */
  ratesManager() :
    isNormalized( false ), 
    rates22( interactionType::indexProcessesInclusive22.size(), 0 ),
    rates23( interactionType::indexProcessesInclusive23.size(), 0 ), 
    rates32( interactionType::indexProcessesInclusive32.size(), 0 ),
    simulatedNumber22( interactionType::indexProcessesInclusive22.size(), 0 ),
    simulatedNumber23( interactionType::indexProcessesInclusive23.size(), 0 ),
    simulatedNumber32( interactionType::indexProcessesInclusive32.size(), 0 ),
    ratesFlav22( 7, 0 ),
    ratesFlav23( 7, 0 ),
    ratesFlav32( 7, 0 ),
    ratesFlavTotal( 7, 0 ),
    nCollectedParticleBasedRates( 7, 0 )
  {
#pragma omp critical
    if ( !faksInitialized ) initializeFaks();
  }

  /** @brief Destructor */
  ~ratesManager() {};
  
  
  ratesManager& operator+=( const ratesManager& rhs );
  ratesManager& operator/=( const double& arg );
  ratesManager& operator*=( const double& arg );
  
  const ratesManager operator+(const ratesManager& rhs) const { return ratesManager(*this) += rhs; }
  const ratesManager operator/(const double& arg) const { return ratesManager(*this) /= arg; }
  const ratesManager operator*(const double& arg) const { return ratesManager(*this) *= arg; }
  
  

  /** @brief Clears all stored rates etc. */
  void clear();
  
  /** @brief Add sampled probability for 2->2 and 2->3 processes */
  int add( const GENERIC_COLL_TYPE genericType, const FLAVOR_TYPE _F1, const FLAVOR_TYPE _F2, const double _P );

  /** @brief Add sampled probability for 3->2 processes */
  int add( const GENERIC_COLL_TYPE genericType, const FLAVOR_TYPE _F1, const FLAVOR_TYPE _F2, const FLAVOR_TYPE _F3, double _P );
  
  /** @brief Add entry for calculation of particle based rates */
  void addParticleBasedRates( const ParticlePrototype& _particle, const UNIT_TYPE unit = GeV );

  /** @brief Finalize calculation of particle based rates */
  void prepareParticleBasedAverages();
  
  /** @brief Normalize the stored rates such that the interaction rate per particle and per time interval can be given (applicable also with heavy quarks) */
  int normalizeRates( const int Ng, const int Nq, const int Nqbar, const int Nc, const int Ncbar, const int Nb, const int Nbbar, const double dt, const double dV );
  /** @brief Normalize the stored rates such that the interaction rate per particle and per time interval can be given (this function can be used if no heavy quarks are present) */
  int normalizeRates( const int Ng, const int Nq, const int Nqbar, const double dt, const double dV ) { return normalizeRates( Ng, Nq, Nqbar, 0, 0, 0, 0, dt, dV ); };
  /** @brief Set all rates to 0 and set normalized flag */
  int normalizeRates();
  /** @brief Normalize the stored rates such that the interaction rate per particle and per time interval can be given (applicable also with heavy quarks), used in routines where vSigma was stored (energyloss, create_mfp_tables...) */
  int normalizeRates( const double ng, const double nq, const double nqbar );
  
  
  /** @brief Get the (total) rate for a given particle type */
  double getRate( const FLAVOR_TYPE particleType, const UNIT_TYPE unit = fm ) const;
  /** @brief Get the rate for a given particle type in given collision types (22, 23, 32) */
  double getRate( const FLAVOR_TYPE particleType, const GENERIC_COLL_TYPE collType, const UNIT_TYPE unit = fm ) const;
  /** @brief Get the rate per inclusive process per volume */
  double getRateInclusive( const int inclCollType, const GENERIC_COLL_TYPE genCollType, const UNIT_TYPE unit = fm ) const;

  /** @brief Get the mean free path for a given particle type */
  double getLambda( const FLAVOR_TYPE particleType, const UNIT_TYPE unit = fm ) const;
  /** @brief Get the mean free path for a given particle type and also specifying the velocity which is important for heavy quarks */
  double getLambda( const FLAVOR_TYPE particleType, const double velocity, const UNIT_TYPE unit = fm ) const;
  /** @brief Get the mean free path for a given particle type in given collision types (22, 23, 32) */
  double getLambda( const FLAVOR_TYPE particleType, const GENERIC_COLL_TYPE collType, const UNIT_TYPE unit = fm ) const;
  /** @brief Get the mean free path for a given particle type in given collision types (22, 23, 32) and also specifying the velocity which is important for heavy quarks */
  double getLambda( const FLAVOR_TYPE particleType, const GENERIC_COLL_TYPE collType, const double velocity, const UNIT_TYPE unit = fm ) const;


  /** @brief Switch indicating whether the collection phase of probabilities has been completed and the rates have been normalized */
  bool isNormalized;

  /** @brief Calculate the index frome the generic flavor */
  unsigned int flavIndex( const FLAVOR_TYPE particleType ) const
  {
    switch ( particleType )
    {
      case gluon:
        return 0;

      case light_quark:
        return 1;

      case anti_light_quark:
        return 2;

      case charm:
      case anti_charm:
      case bottom:
      case anti_bottom:
        return static_cast<int>(particleType) - 4;

      default:
        std::string errMsg = "Particle type not found in ratesManager::getRate. Unrecoverable error.";
        throw eRatesManager_error( errMsg );
        break;
    }
  }

private:

  /** @brief Collect the probabilities for all 2->2 interaction types (each index corresponds to a interaction type as set in class interactionType) */
  std::vector<double> rates22;

  /** @brief Collect the probabilities for all 2->3 interaction types (each index corresponds to a interaction type as set in class interactionType) */
  std::vector<double> rates23;

  /** @brief Collect the probabilities for all 2->3 interaction types (each index corresponds to a interaction type as set in class interactionType) */
  std::vector<double> rates32;

  /** @brief array to store values of e.g. interactionType::getInvolvedInitialInclusive( gluon, c22 ) */
  static std::vector< std::vector< double > > faks22;

  /** @brief array to store values of e.g. interactionType::getInvolvedInitialInclusive( gluon, c23 ) */
  static std::vector< std::vector< double > > faks23;

  /** @brief array to store values of e.g. interactionType::getInvolvedInitialInclusive( gluon, c32 ) */
  static std::vector< std::vector< double > > faks32;
  
  /** @brief flag to indicate, whether faks22 etc are initialized */
  static bool faksInitialized;

  /** @brief Stores 2->2 rates for every flavor after ratesManager::normalizeRates has been called.  */
  std::vector<double> ratesFlav22;

  /** @brief Stores 2->3 rates for every flavor after ratesManager::normalizeRates has been called.  */
  std::vector<double> ratesFlav23;

  /** @brief Stores 3->2 rates for every flavor after ratesManager::normalizeRates has been called.  */
  std::vector<double> ratesFlav32;
  
  /** @brief Stores total (22 + 23 + 32) rate after ratesManager::normalizeRates has been called. */
  std::vector<double> ratesFlavTotal;

  /** @brief The number of collisions stored for each type (each index corresponds to a interaction type as set in class interactionType) */
  std::vector<double> simulatedNumber22;
  /** @brief The number of collisions stored for each type (each index corresponds to a interaction type as set in class interactionType) */
  std::vector<double> simulatedNumber23;
  /** @brief The number of collisions storeed for each type. (each index corresponds to a interaction type as set in class interactionType) */
  std::vector<double> simulatedNumber32;

  /** 
   * @brief number of encountered particles while particle based
   * calculation of rates 
   **/
  std::vector<int> nCollectedParticleBasedRates;

  /** @brief initialize faks */
  void initializeFaks(void);

};



#endif // RATESMANAGER_H
// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;
