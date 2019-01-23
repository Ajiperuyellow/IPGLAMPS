//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/mfp_data.h $
//$LastChangedDate: 2016-09-20 16:39:07 +0200 (Di, 20. Sep 2016) $
//$LastChangedRevision: 2429 $
//$LastChangedBy: senzel $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


#ifndef MFP_DATA_H
#define MFP_DATA_H

#include <vector>
#include <string>
#include <stdexcept>
#include <boost/lexical_cast.hpp>
#include "particleprototype.h"
#include "configurationbase.h"
#include "unit_enum.h"

enum MFP_ERROR_CODE {MFP_SUCCESS, MFP_READOUT_ERROR, MFP_FILE_ERROR, MFP_SIZE_ERROR};  


class mfp_data
{
public:
  mfp_data();
  mfp_data( const double _T, const bool read_file_arg = false );
  ~mfp_data();
  
  double getLambda( const double E, const FLAVOR_TYPE F = gluon ) const;
  double getLambdaTherm( FLAVOR_TYPE F = gluon ) const;
  double getR22(const double E, const FLAVOR_TYPE F = gluon ) const;
  double getR23(const double E, const FLAVOR_TYPE F = gluon ) const;
  double getR32(const double E, const FLAVOR_TYPE F = gluon ) const;
  void displayData() const;   //for debugging purposes
  
private:
  MFP_ERROR_CODE readTables( const std::string & filename, const FLAVOR_TYPE _F );
  int findIndices( const double _E, const std::vector<double> data_E ) const;
  void polint(const double xa[], const double ya[], const int n, const double x,double *y,double *dy) const;
  
  std::vector<double> data_E_gluon;
  std::vector<double> data_E_quark;
  std::vector<double> data_E_charm;
  std::vector<double> data_E_bottom;
  std::vector<double> data_lambda_gluon_fm;
  std::vector<double> data_lambda_gluon_gev;
  std::vector<double> data_lambda_quark_fm;
  std::vector<double> data_lambda_quark_gev;
  std::vector<double> data_lambda_charm_fm;
  std::vector<double> data_lambda_charm_gev;
  std::vector<double> data_lambda_bottom_fm;
  std::vector<double> data_lambda_bottom_gev;
  std::vector<double> data_R22_gluon;
  std::vector<double> data_R23_gluon;
  std::vector<double> data_R32_gluon;
  std::vector<double> data_R22_quark;
  std::vector<double> data_R23_quark;
  std::vector<double> data_R32_quark;
  std::vector<double> data_R22_charm;
  std::vector<double> data_R23_charm;
  std::vector<double> data_R32_charm;
  std::vector<double> data_R22_bottom;
  std::vector<double> data_R23_bottom;
  std::vector<double> data_R32_bottom;
  
  double data_lambda_therm_gluon_fm;
  double data_lambda_therm_gluon_gev;
  double data_lambda_therm_quark_fm;
  double data_lambda_therm_quark_gev;
  double data_lambda_therm_charm_fm;
  double data_lambda_therm_charm_gev;
  double data_lambda_therm_bottom_fm;
  double data_lambda_therm_bottom_gev;
  
  bool includeQuarkMFP;
  bool includeCharmMFP;
  bool includeBottomMFP;
  
  /** @brief indicates whether tables are loaded. Set by default to false since MFP tables are still with the original GB and not suited for realistic runs. */
  bool read_file;
  
  /** @brief Return T as string [MeV] */
  std::string getString_T( const double _T ) const
  {
    return boost::lexical_cast<std::string>( int(_T*1000 + 0.001) ); //*1000 to convert to MeV
  };

};


class mfpForHeavyIonCollision
{
public:
  mfpForHeavyIonCollision( configBase* const _c );
  ~mfpForHeavyIonCollision();
  
  double getMeanFreePath( const double _E, const FLAVOR_TYPE _F, const double _T, const double _ngTest, const double _nqTest, const UNIT_TYPE _unit ) const;
  
  void loadData();
  
  
private:  
  double fugacityDependence( const double _gluonFugacity, const double _quarkFugacity ) const;
  double getGluonDensity( const double T, const double fugacity = 1 ) const ;
  double getQuarkDensity( const double T, const double fugacity = 1 ) const;

  int getStartIndexForTemperatureValues( const double _T, const unsigned int _nValues ) const;
  
  /** @brief interpolation routine */
  void polint( const double xa[], const double ya[], const int n, const double x, double *y, double *dy ) const;
  
  std::vector<mfp_data> mfpData;
  std::vector<double> temperaturesForMfpData;
  
  double fitParameterGluon;
  double fitParameterQuark;
  double fugacityDependenceEvaluatedAtFugacityOne;
  
  bool data_loaded;
  
  configBase * const theConfig;
};


/** @brief exception class for handling unexpected behaviour when reading the files containing the MPF data */
class eMFP_read_error : public std::runtime_error
{
public:
  explicit eMFP_read_error( const std::string& what ) : std::runtime_error( what ) {};

  virtual ~eMFP_read_error() throw() {};
};

#endif
