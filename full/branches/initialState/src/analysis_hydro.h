//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/analysis_hydro.h $
//$LastChangedDate: 2014-02-12 20:15:39 +0100 (水, 12  2月 2014) $
//$LastChangedRevision: 1620 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifndef ANALYSIS_HYDRO_H
#define ANALYSIS_HYDRO_H

#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <list>

#include "particle.h"
#include "configuration.h"
#include "hydroParticleType.h"
#include "FPT_compare.h"

#include "boost/multi_array.hpp"
#include <cassert> 

class config; 
 
using namespace std;
using namespace ns_casc;


class anaHydroData
{
public:
  anaHydroData( const int _nTypes, const int _realTimeSteps, const int _numberColumnsX, const int _numberColumnsY,
                const int _numberColumnsZ, const double _columnWidthX, const double _columnWidthY, const double _columnWidthZ, const string _name);

  typedef boost::multi_array<double, 7> array_typeTmunu;    
  typedef boost::multi_array<double, 6> array_typeNmu;    
  typedef boost::multi_array<double, 2> array_typeMidRapOb;    
 
  array_typeTmunu enDistTmunu;
  array_typeNmu enDistNmu;
  
  array_typeMidRapOb midRapOb_dN_dy;
  array_typeMidRapOb midRapOb_dEt_dy;
  array_typeMidRapOb midRapOb_v2;
  array_typeMidRapOb midRapOb_v4;
  
  const static int mu = 4;
  const static int nu = 4;  
  
  double columnWidthX, columnWidthY, columnWidthZ;
  int numberColumnsX, numberColumnsY, numberColumnsZ;
  string name;
};


//--------------------------------------
//Initialize for hydro analysis
class anaHydro
{
  public:
    anaHydro( config * const theConfig, const int _NNtimeStepsForAnaHydro );
   ~anaHydro();    
    
    hydroParticleType theHydroParticleType;
    
    //----------------------------------//  
    //analysis for hydro    
    void hydroDistribution(const int nn, anaHydroData * ad);
    void hydroDistributionMidRap(const int nn, const double time, const double timeshift, anaHydroData * ad);
    int particleListgetFinalNumber();
    //----------------------------------//
    
    double boxLengthX, boxLengthY, boxLengthZ;
    int nTypes, realTimeSteps, finalNumberInSimulation;
    
    const static int mu = 4;
    const static int nu = 4;
       
    anaHydroData * anaHydroProfileNormal;
    anaHydroData * anaHydroProfileArrow;
    anaHydroData * anaHydroProfileMidRapNormal;
    anaHydroData * anaHydroProfileMidRapArrow;
    
  private:
    //---------------------------------  
    config * const theConfig ;
    //---------------------------------    
    
    void showInfosAtMidRapidity(const int nn, const double volumeBin, const double rapidityRangeInZ, const anaHydroData * ad);
    int getRealNumberOfTimeStepsForHydro();    
  };

  
#endif
