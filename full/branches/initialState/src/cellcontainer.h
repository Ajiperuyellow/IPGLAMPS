//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/cellcontainer.h $
//$LastChangedDate: 2018-04-28 19:57:48 +0200 (Sa, 28. Apr 2018) $
//$LastChangedRevision: 2742 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


#ifndef CELLCONTAINER_H
#define CELLCONTAINER_H

#include <list>
#include <stdexcept>

#include "ratesmanager.h"
#include "coordinateBins.h"
#include "particle.h"
#include "configuration.h"


class cornerCoordinates
{
public:
  cornerCoordinates();
  cornerCoordinates( const double _x_min, const double _x_max, const double _y_min, const double _y_max, const int _etaIndex );
  ~cornerCoordinates() {};

  void setCorners( const double _x_min, const double _x_max, const double _y_min, const double _y_max, const int _etaIndex );
  void clear() { setCorners( 0, 0, 0, 0, 0 ); }
  double getVolume( const coordinateEtaBins& _etaBins, const double _time ) const;


  double x_min;
  double x_max;
  double y_min;
  double y_max;
//   double z_min;
//   double z_max;
//   double eta_min;
//   double eta_max;
  int etaIndex;
};








class cellContainer
{
public:
  cellContainer();
  ~cellContainer() {};

  std::list<int> particleList;

  ratesManager rates,ratesSave;

  void clear();
  void resetStoredValues();
  bool empty() const { return particleList.empty(); }
  int size() const { return particleList.size(); }
  void setCoordinates( const int _index, const double _dx, const int _nx, const double _sizeX, const double _dy, const int _ny, const double _sizeY );
  void prepareAverages();
  void writeAveragesToParticle( Particle& _particle ) const;

  void prepareThermodynamicsWithRadialBoost(double volume, int Ntest);
  void prepareThermodynamicsWithRadialBoostWNeighbors(double _volume, int _Ntest, std::vector< int >& _allParticlesListWNeighbors);
  void prepareThermodynamicsCellBoost(double volume, int Ntest);
  void prepareThermodynamicsCellBoostWNeighbors(double _volume, int _Ntest, std::vector< int >& _allParticlesListWNeighbors);
  void prepareAverageRates(double volume, int Ntest);
  void prepareAverageRatesWNeighbors(double volume, int Ntest,std::vector< int >& _allParticlesListWNeighbors);
  void clearThermodynamicsAndRates();
  void defineAsEmpty();
  
  int index;
  cornerCoordinates corner;

  bool averagesPrepared,averagesPreparedSave;

  int nCollectedAll2223;
  int nCollected22;
  int nCollected23;

  double alpha_s_22;
  double alpha_s_23;

  double md2g_scaled_22;
  double md2q_scaled_22;
  double md2g_scaled_23;
  double md2q_scaled_23;

  double sigma_22;
  double sigma_23;

  double lambdaScaled;

  void setValuesToRadialValues();
  double volume;  
  int numberOfParticles,numberOfGluons,numberOfQuarks,numberOfCollectedRateObjects,numberOfCollectedRateSaveObjects; 
  double getAveraged_md2g() const {return md2g;}
  double getAveraged_md2q() const {return md2q;}

  /** @brief return the averaged value **/
  double getParticleDensity() const;

  /** @brief return the averaged value **/
  double getGluonDensity() const;

  int getNumberOfParticles() const;
  
  /** @brief return the averaged value **/
  double getQuarkDensity() const;
  
  /** @brief return the averaged value **/
  double getVolumeFromAVG() const {return volumeFromAVG;};  
  
  /** @brief return the averaged value **/
  double getEnergyDensity() const;
  double getEnergyDensityCellRadial() const;
  
  /** @brief return the averaged value of gamma and boostvector**/
  double getGamma() const;  
  void getBoostvector(VectorXYZ & boostbeta) const;
  
  /** @brief return the averaged value, epsilon/(3*n) **/
  double getEffectiveTemperature() const;
  
  /** @brief return the averaged energy in a comoving frame **/
  double transformEnergyToComovingFrame( VectorEPxPyPz & P ) const;
  
      //Instead of ring-thermodynamic, we use cells:
  double md2g; ///< The averaged value of the gluon Debye mass
  double md2q; ///< The averaged value of the quark Debye mass
  
  double volumeFromAVG;
  double gamma;
  VectorXYZ boostvector;
  double energyDensity;  
  double particleDensity;
  double gluonDensity;
  double quarkDensity;
  
  double gammaCellRadial;
  double energyDensityCellRadial;
  double particleDensityCellRadial;
  double gluonDensityCellRadial;
  double quarkDensityCellRadial;
private:



  
  
};


/** @brief exception class for handling unexpected critical behaviour within cell objects  */
class eCell_error : public std::runtime_error
{
public:
  explicit eCell_error( const std::string& what ) : std::runtime_error( what ) {};

  virtual ~eCell_error() throw() {};
};

#endif // CELLCONTAINER_H
// kate: indent-mode cstyle; space-indent on; indent-width 0;
