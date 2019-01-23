//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/analysis.cpp $
//$LastChangedDate: 2019-01-05 18:02:58 +0100 (Sa, 05. Jan 2019) $
//$LastChangedRevision: 2917 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <iostream>
#include <list>
#include <math.h>
#include <sstream>
#include <stdlib.h>

#include "analysis.h"
#include "binning.h"
#include "binning2.h"
#include "global.h"
#include "particle.h"
#include "scattering22_hydro.h"
#include "cellcontainer.h"
#include "lorentz.h"
#include "random.h"

using namespace std;
using namespace ns_casc;


void analysis::handle_output_studies( OUTPUT_SCHEME _outputScheme )
{
  

  // By default switches are set off (please do not edit here, but define output scheme and change below, for more information see https://th.physik.uni-frankfurt.de/~bamps/cgi-bin/trac/wiki/RestrictedWiki/AnalysisOutputScheme):
  studyParticleSpectra = theConfig->doOutput_particleSpectraOutput();  // write out particle spectra of all species at the beginning and end
  studyParticleSpectraAllSteps = theConfig->doOutput_particleSpectraOutputAllSteps();  // write out particle spectra of all species at each time step
  studyDetailedParticleOutput = theConfig->doOutput_detailedParticleOutput();  // print all particles at the beginning and end
  studyDetailedParticleOutputAllSteps = theConfig->doOutput_detailedParticleOutputAllSteps();  // print all particles at each time step
  studyHydro = theConfig->getAnalyseForHydro(); // Ioannis hydro analysis routines
  studyDetailedParticlePosMomOutput = theConfig->doOutput_detailedParticlePosMomOutput(); // prints tables of CURRENTLY ACTIVE particle position and momenta, used for Bennys correlation search program (extern/offline)
  
  studyV2 = false;  // write v2(pt) and mean v2(y) and mean pt (y) for various times
  studyBoostDistribution = false; // print the boost distributions for 2->3 and 3->2
  studyHeavyQuarkProduction = false;  // study heavy quark production
  studyJetTracking = false; // track jets
  studyParticleAndCollisionNumbers = false; // print the numbers of all particles, the numbers of the collisions and many more as a function of time
  studyDndy = false; // study dn/dy and dEt/dy
  studyCorrelations = false;
  studyCellwiseVideo = false;
  noAnalysis = false;
  
  
  //---- defining standard rapidity ranges ----
  // only use positiv ranges since the investigated collision systems usually are symmetric in +-y and we therefore only compare the absolute value of y
  analysisRapidityRange yRange;
  yRange.reset( 0, 0.5 );
  rapidityRanges.push_back(yRange);
  yRange.reset( 0, 0.8 );
  rapidityRanges.push_back(yRange);
  yRange.reset( 0, 1.0 );
  rapidityRanges.push_back(yRange);
  yRange.reset( 0, 1.5 );
  rapidityRanges.push_back(yRange);
  yRange.reset( 0, 2.0 );
  rapidityRanges.push_back(yRange);
  yRange.reset( 0, 2.5 );
  rapidityRanges.push_back(yRange);  
  yRange.reset( 0, infinity );
  rapidityRanges.push_back(yRange);  
  //---- defining rapidity ranges ----
  
  // add here a new case for your outpute scheme (which you can create in configuration.h)
  switch ( _outputScheme )
  {
    case allOutputs:
      studyParticleSpectra = true;
      studyParticleSpectraAllSteps = true;
      studyDetailedParticleOutput = true;
      studyDetailedParticleOutputAllSteps = true;
      studyHydro = true;
      studyV2 = true;
      studyBoostDistribution = true;
      studyHeavyQuarkProduction = true;
      studyJetTracking = true;
      studyParticleAndCollisionNumbers = true;
      break;
      
    case testrun_observables:
      studyParticleSpectra = true;
      studyV2 = true;
      studyParticleAndCollisionNumbers = true;
      studyDndy = true;
      break;
    
    case heavy_quark_production:
      studyHeavyQuarkProduction = true;
      break;

    case correlationsearch:
      studyDetailedParticlePosMomOutput = true;
      break;
      
    case phenix_v2:
      studyV2 = true;
      rapidityRanges.clear();
      yRange.reset( 0, 0.35 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 0.5 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 1.0 );
      rapidityRanges.push_back(yRange);
      yRange.reset( 0, 2.0 );
      rapidityRanges.push_back(yRange);
      break;

    case correlations:
      studyCorrelations = 1;
      studyDndy = 0;
      studyParticleAndCollisionNumbers = 0;
      rapidityRanges.clear();
      yRange.reset( 0, 0.5 );
      rapidityRanges.push_back(yRange);
      break;  
    case onlyVideo:
      studyCellwiseVideo=true;
      yRange.reset( 0, 0.5 );
      rapidityRanges.push_back(yRange);
      break;  
    case nothingToAnalyze:
      noAnalysis = true;
      rapidityRanges.clear();
      yRange.reset( 0, 0.5 );
      rapidityRanges.push_back(yRange);
    default:
      break;
  }

}


analysis::analysis( config * const _theConfig, heavyIonCollision * const _theHIC ):
    theConfig( _theConfig ), theHIC( _theHIC ),
    filename_prefix( theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() ),
    outputScheme( theConfig->getOutputScheme() )
{
  //---get time and date---
  time( &start );
  //-----------------------
  
  handle_output_studies( outputScheme );

  //-----------------------  
  switch ( theConfig->getAnaTimeStepsMode() )
  {
  case generalHICAnaTimeSteps:
    anaTimeStep_generalHIC();
    break;
  case detailedHICAnaTimeSteps:
    anaTimeStep_detailed();
    break;
  case heavyQuarkAnaTimeSteps:
    anaTimeStep_heavyQuark();
    break;
  case hydroParametrizationAnaTimeSteps:
    anaTimeStep_hydroParametrization();    
    break;    
  case correlationAnaTimeSteps:
    //anaTimeStep_heavyQuark();
    //anaTimeStep_correlations();
    anaTimeStepVideo();
    break;  
  default:
    cout << "Error in choosing timestep mode for analysis: Default one is chosen" << endl;    
    anaTimeStep_generalHIC();
  }
  //----------------------- 
  
  
  
  //-----------------------
  int NNtimeStepsForAnaHydro = getRealNumberOfTimeStepsForAnaHydro();
  theAnaHydro = new anaHydro( theConfig, NNtimeStepsForAnaHydro );
  //-----------------------  
  
  //--------------------------------------  
  jetTracking_PT = 10.0;

  //---- initialisation of PT-binning ----
  minPT = 1.4;
  maxPT = 34.4;
  binWidthPT = 1.0;
  numberBinsPT = int(( maxPT - minPT + 0.001 ) / binWidthPT );
  ptBinsDY1_gluons = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  ptBinsDY16_gluons = new vector<double>[nTimeSteps+2];
  ptBinsDY2_gluons = new vector<double>[nTimeSteps+2];
  ptBinsDY3_gluons = new vector<double>[nTimeSteps+2];
  ptBinsDY4_gluons = new vector<double>[nTimeSteps+2];
  ptBinsDY5_gluons = new vector<double>[nTimeSteps+2];
  ptBinsAll_gluons = new vector<double>[nTimeSteps+2];
  
  ptBinsDY1_quarks = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  ptBinsDY16_quarks = new vector<double>[nTimeSteps+2];
  ptBinsDY2_quarks = new vector<double>[nTimeSteps+2];
  ptBinsDY3_quarks = new vector<double>[nTimeSteps+2];
  ptBinsDY4_quarks = new vector<double>[nTimeSteps+2];
  ptBinsDY5_quarks = new vector<double>[nTimeSteps+2];
  ptBinsAll_quarks = new vector<double>[nTimeSteps+2];
  
  ptBinsDY1_ups = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  ptBinsDY16_ups = new vector<double>[nTimeSteps+2];
  ptBinsDY2_ups = new vector<double>[nTimeSteps+2];
  ptBinsDY3_ups = new vector<double>[nTimeSteps+2];
  ptBinsDY4_ups = new vector<double>[nTimeSteps+2];
  ptBinsDY5_ups = new vector<double>[nTimeSteps+2];
  ptBinsAll_ups = new vector<double>[nTimeSteps+2];
  
  ptBinsDY1_downs = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  ptBinsDY16_downs = new vector<double>[nTimeSteps+2];
  ptBinsDY2_downs = new vector<double>[nTimeSteps+2];
  ptBinsDY3_downs = new vector<double>[nTimeSteps+2];
  ptBinsDY4_downs = new vector<double>[nTimeSteps+2];
  ptBinsDY5_downs = new vector<double>[nTimeSteps+2];
  ptBinsAll_downs = new vector<double>[nTimeSteps+2];
  
  ptBinsDY1_stranges = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  ptBinsDY16_stranges = new vector<double>[nTimeSteps+2];
  ptBinsDY2_stranges = new vector<double>[nTimeSteps+2];
  ptBinsDY3_stranges = new vector<double>[nTimeSteps+2];
  ptBinsDY4_stranges = new vector<double>[nTimeSteps+2];
  ptBinsDY5_stranges = new vector<double>[nTimeSteps+2];
  ptBinsAll_stranges = new vector<double>[nTimeSteps+2];
  
  ptBinsDY1_anti_ups = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  ptBinsDY16_anti_ups = new vector<double>[nTimeSteps+2];
  ptBinsDY2_anti_ups = new vector<double>[nTimeSteps+2];
  ptBinsDY3_anti_ups = new vector<double>[nTimeSteps+2];
  ptBinsDY4_anti_ups = new vector<double>[nTimeSteps+2];
  ptBinsDY5_anti_ups = new vector<double>[nTimeSteps+2];
  ptBinsAll_anti_ups = new vector<double>[nTimeSteps+2];
  
  ptBinsDY1_anti_downs = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  ptBinsDY16_anti_downs = new vector<double>[nTimeSteps+2];
  ptBinsDY2_anti_downs = new vector<double>[nTimeSteps+2];
  ptBinsDY3_anti_downs = new vector<double>[nTimeSteps+2];
  ptBinsDY4_anti_downs = new vector<double>[nTimeSteps+2];
  ptBinsDY5_anti_downs = new vector<double>[nTimeSteps+2];
  ptBinsAll_anti_downs = new vector<double>[nTimeSteps+2];
  
  ptBinsDY1_anti_stranges = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  ptBinsDY16_anti_stranges = new vector<double>[nTimeSteps+2];
  ptBinsDY2_anti_stranges = new vector<double>[nTimeSteps+2];
  ptBinsDY3_anti_stranges = new vector<double>[nTimeSteps+2];
  ptBinsDY4_anti_stranges = new vector<double>[nTimeSteps+2];
  ptBinsDY5_anti_stranges = new vector<double>[nTimeSteps+2];
  ptBinsAll_anti_stranges = new vector<double>[nTimeSteps+2];
  
  ptBinsDY1_all = new vector<double>[nTimeSteps+2];          //+2 because of initial and final timesteps
  ptBinsDY16_all = new vector<double>[nTimeSteps+2];
  ptBinsDY2_all = new vector<double>[nTimeSteps+2];
  ptBinsDY3_all = new vector<double>[nTimeSteps+2];
  ptBinsDY4_all = new vector<double>[nTimeSteps+2];
  ptBinsDY5_all = new vector<double>[nTimeSteps+2];
  ptBinsAll_all = new vector<double>[nTimeSteps+2];
  
  for ( int j = 0; j < nTimeSteps + 2; j++ )          //+2 because of initial and final timesteps
  {
    for ( int i = 0; i <= numberBinsPT; i++ )
    {
      ptBinsDY1_gluons[j].push_back( 0 );
      ptBinsDY16_gluons[j].push_back( 0 );
      ptBinsDY2_gluons[j].push_back( 0 );
      ptBinsDY3_gluons[j].push_back( 0 );
      ptBinsDY4_gluons[j].push_back( 0 );
      ptBinsDY5_gluons[j].push_back( 0 );
      ptBinsAll_gluons[j].push_back( 0 );
      
      ptBinsDY1_quarks[j].push_back( 0 );
      ptBinsDY16_quarks[j].push_back( 0 );
      ptBinsDY2_quarks[j].push_back( 0 );
      ptBinsDY3_quarks[j].push_back( 0 );
      ptBinsDY4_quarks[j].push_back( 0 );
      ptBinsDY5_quarks[j].push_back( 0 );
      ptBinsAll_quarks[j].push_back( 0 );
      
      ptBinsDY1_ups[j].push_back( 0 );
      ptBinsDY16_ups[j].push_back( 0 );
      ptBinsDY2_ups[j].push_back( 0 );
      ptBinsDY3_ups[j].push_back( 0 );
      ptBinsDY4_ups[j].push_back( 0 );
      ptBinsDY5_ups[j].push_back( 0 );
      ptBinsAll_ups[j].push_back( 0 );
      
      ptBinsDY1_downs[j].push_back( 0 );
      ptBinsDY16_downs[j].push_back( 0 );
      ptBinsDY2_downs[j].push_back( 0 );
      ptBinsDY3_downs[j].push_back( 0 );
      ptBinsDY4_downs[j].push_back( 0 );
      ptBinsDY5_downs[j].push_back( 0 );
      ptBinsAll_downs[j].push_back( 0 );
      
      ptBinsDY1_stranges[j].push_back( 0 );
      ptBinsDY16_stranges[j].push_back( 0 );
      ptBinsDY2_stranges[j].push_back( 0 );
      ptBinsDY3_stranges[j].push_back( 0 );
      ptBinsDY4_stranges[j].push_back( 0 );
      ptBinsDY5_stranges[j].push_back( 0 );
      ptBinsAll_stranges[j].push_back( 0 );
      
      ptBinsDY1_anti_ups[j].push_back( 0 );
      ptBinsDY16_anti_ups[j].push_back( 0 );
      ptBinsDY2_anti_ups[j].push_back( 0 );
      ptBinsDY3_anti_ups[j].push_back( 0 );
      ptBinsDY4_anti_ups[j].push_back( 0 );
      ptBinsDY5_anti_ups[j].push_back( 0 );
      ptBinsAll_anti_ups[j].push_back( 0 );
      
      ptBinsDY1_anti_downs[j].push_back( 0 );
      ptBinsDY16_anti_downs[j].push_back( 0 );
      ptBinsDY2_anti_downs[j].push_back( 0 );
      ptBinsDY3_anti_downs[j].push_back( 0 );
      ptBinsDY4_anti_downs[j].push_back( 0 );
      ptBinsDY5_anti_downs[j].push_back( 0 );
      ptBinsAll_anti_downs[j].push_back( 0 );
      
      ptBinsDY1_anti_stranges[j].push_back( 0 );
      ptBinsDY16_anti_stranges[j].push_back( 0 );
      ptBinsDY2_anti_stranges[j].push_back( 0 );
      ptBinsDY3_anti_stranges[j].push_back( 0 );
      ptBinsDY4_anti_stranges[j].push_back( 0 );
      ptBinsDY5_anti_stranges[j].push_back( 0 );
      ptBinsAll_anti_stranges[j].push_back( 0 );
      
      ptBinsDY1_all[j].push_back( 0 );
      ptBinsDY16_all[j].push_back( 0 );
      ptBinsDY2_all[j].push_back( 0 );
      ptBinsDY3_all[j].push_back( 0 );
      ptBinsDY4_all[j].push_back( 0 );
      ptBinsDY5_all[j].push_back( 0 );
      ptBinsAll_all[j].push_back( 0 );
    }
  }
  //write the bin labels
  for ( int i = 0; i < numberBinsPT; i++ )
  {
    ptBinLabels.push_back( minPT + ( i * binWidthPT ) + ( binWidthPT / 2 ) );
  }
  cout << "number of bins: " << numberBinsPT << "  binWidth: " << binWidthPT << endl;
  
  //---- initialisation of softPT-binning ----
  maxPTSoft = 3.0;
  binWidthSoftPT = 0.1;
  numberBinsSoftPT = int(( maxPTSoft + 0.001 ) / binWidthSoftPT );
  ptBinsSoftAll_gluons = new vector<double>[nTimeSteps+1];
  ptBinsSoftAll_quarks = new vector<double>[nTimeSteps+1];
  ptBinsSoftAll_ups = new vector<double>[nTimeSteps+1];
  ptBinsSoftAll_downs = new vector<double>[nTimeSteps+1];
  ptBinsSoftAll_stranges = new vector<double>[nTimeSteps+1];
  ptBinsSoftAll_anti_ups = new vector<double>[nTimeSteps+1];
  ptBinsSoftAll_anti_downs = new vector<double>[nTimeSteps+1];
  ptBinsSoftAll_anti_stranges = new vector<double>[nTimeSteps+1];
  ptBinsSoftAll_all = new vector<double>[nTimeSteps+1];
  for ( int j = 0; j < nTimeSteps + 1; j++ )          //+1 because of initial timestep
  {
    for ( int i = 0; i < numberBinsSoftPT; i++ )
    {
      ptBinsSoftAll_gluons[j].push_back( 0 );
      ptBinsSoftAll_quarks[j].push_back( 0 );
      ptBinsSoftAll_ups[j].push_back( 0 );
      ptBinsSoftAll_downs[j].push_back( 0 );
      ptBinsSoftAll_stranges[j].push_back( 0 );
      ptBinsSoftAll_anti_ups[j].push_back( 0 );
      ptBinsSoftAll_anti_downs[j].push_back( 0 );
      ptBinsSoftAll_anti_stranges[j].push_back( 0 );
      ptBinsSoftAll_all[j].push_back( 0 );
    }
  }
  //write the bin labels
  for ( int i = 0; i < numberBinsSoftPT; i++ )
  {
    ptSoftBinLabels.push_back(( i * binWidthSoftPT ) + ( binWidthSoftPT / 2 ) );
  }
  //--------------------------------------
  
  cout << "number of test particles = " << theConfig->getTestparticles() << endl;
  
  //-------- differential v2 bins ---------------------
  min_v2_pt = 0.0;
  max_v2_pt = 5.0;
  binWidth_v2_pt = 0.2;
  nBins_v2_pt = int(( max_v2_pt - min_v2_pt + 0.001 ) / binWidth_v2_pt );
  
  rapidityBins.clear();
  rapidityBins.push_back( 30000 );
  rapidityBins.push_back( 2.5 );
  rapidityBins.push_back( 2.0 );
  rapidityBins.push_back( 1.5 );
  rapidityBins.push_back( 1.0 );
  rapidityBins.push_back( 0.8 );
  rapidityBins.push_back( 0.5 );
  nRapidityBins = rapidityBins.size();
  
  differential_v2_values = new vector< vector<double> >[nTimeSteps+1];
  differential_v2_counts = new vector< vector<double> >[nTimeSteps+1];
  mean_PT = new vector<double>[nTimeSteps+1];
  integrated_v2 = new vector<double>[nTimeSteps+1];
  integrated_v2_counts = new vector<double>[nTimeSteps+1];
  
  differential_v2_bins_labels.clear();
  for ( int i = 0; i < nBins_v2_pt; i++ )
  {
    differential_v2_bins_labels.push_back( min_v2_pt + ( ( i + 0.5 ) * binWidth_v2_pt ) );
  }
  
  vector<double> tempVec( nBins_v2_pt, 0 );
  for ( int j = 0; j <= nTimeSteps; j++ )
  {
    differential_v2_values[j].clear();
    differential_v2_values[j].resize( nRapidityBins, tempVec );
    differential_v2_counts[j].clear();
    differential_v2_counts[j].resize( nRapidityBins, tempVec );
    mean_PT[j].clear();
    mean_PT[j].resize( nRapidityBins, 0 );
    integrated_v2[j].clear();
    integrated_v2[j].resize( nRapidityBins, 0 );
    integrated_v2_counts[j].clear();
    integrated_v2_counts[j].resize( nRapidityBins, 0 );
  }
  //-------- differential v2 bins ---------------------


  //---- initialisation of boost binning -----
  numberBinsBoost = 100;
  for ( int i = 0; i < numberBinsBoost; i++ )
  {
    boostDistribution23.push_back( 0 );
    boostDistribution32.push_back( 0 );
  }


  if ( _theConfig->doOutput_progressLog() )
  {
    string filename = _theConfig->getStandardOutputDirectoryName() + "/" + _theConfig->getJobName() + "_progressLog";
    progressLogFile.open( filename.c_str(), ios::out | ios::trunc );
  }
  
  if( studyParticleAndCollisionNumbers )
  {
    string filename = _theConfig->getStandardOutputDirectoryName() + "/" + _theConfig->getJobName() + "_particleAndCollisionNumbers";
    printParticleAndCollisionNumbers.open( filename.c_str(), ios::out | ios::trunc );
  }
  
  //---- initialisation of charm quark counting ----
  charmNmb = new int[nTimeSteps+1];
  charmNmb_eta20 = new int[nTimeSteps+1];
  charmNmb_eta15 = new int[nTimeSteps+1];
  charmNmb_eta10 = new int[nTimeSteps+1];
  charmNmb_eta075 = new int[nTimeSteps+1];
  charmNmb_eta05 = new int[nTimeSteps+1];
  charmNmb_eta035 = new int[nTimeSteps+1];
  charmNmb_strap05 = new int[nTimeSteps+1];

  for ( int i = 0;i <= nTimeSteps;i++ )
  {
    charmNmb[i] = 0;
    charmNmb_eta20[i] = 0;
    charmNmb_eta15[i] = 0;
    charmNmb_eta10[i] = 0;
    charmNmb_eta075[i] = 0;
    charmNmb_eta05[i] = 0;
    charmNmb_eta035[i] = 0;
    charmNmb_strap05[i] = 0;
  }
  //--------------------------------------
  
  //---- initialisation of bottom quark counting ----
  bottomNmb = new int[nTimeSteps+1];
  bottomNmb_eta20 = new int[nTimeSteps+1];
  bottomNmb_eta15 = new int[nTimeSteps+1];
  bottomNmb_eta10 = new int[nTimeSteps+1];
  bottomNmb_eta075 = new int[nTimeSteps+1];
  bottomNmb_eta05 = new int[nTimeSteps+1];
  bottomNmb_eta035 = new int[nTimeSteps+1];

  for ( int i = 0;i <= nTimeSteps;i++ )
  {
    bottomNmb[i] = 0;
    bottomNmb_eta20[i] = 0;
    bottomNmb_eta15[i] = 0;
    bottomNmb_eta10[i] = 0;
    bottomNmb_eta075[i] = 0;
    bottomNmb_eta05[i] = 0;
    bottomNmb_eta035[i] = 0;
  }
  //--------------------------------------


  
  
  
  
  
  
  
  
  
  
  
  
  
  //

}


analysis::~analysis()
{
  finalOutput();
  
  delete theAnaHydro;
  
  delete[] ptBinsDY1_gluons;
  delete[] ptBinsDY16_gluons;
  delete[] ptBinsDY2_gluons;
  delete[] ptBinsDY3_gluons;
  delete[] ptBinsDY4_gluons;
  delete[] ptBinsDY5_gluons;
  delete[] ptBinsAll_gluons;
  
  delete[] ptBinsDY1_quarks;
  delete[] ptBinsDY16_quarks;
  delete[] ptBinsDY2_quarks;
  delete[] ptBinsDY3_quarks;
  delete[] ptBinsDY4_quarks;
  delete[] ptBinsDY5_quarks;
  delete[] ptBinsAll_quarks;
  
  delete[] ptBinsDY1_ups;
  delete[] ptBinsDY16_ups;
  delete[] ptBinsDY2_ups;
  delete[] ptBinsDY3_ups;
  delete[] ptBinsDY4_ups;
  delete[] ptBinsDY5_ups;
  delete[] ptBinsAll_ups;
  
  delete[] ptBinsDY1_downs;
  delete[] ptBinsDY16_downs;
  delete[] ptBinsDY2_downs;
  delete[] ptBinsDY3_downs;
  delete[] ptBinsDY4_downs;
  delete[] ptBinsDY5_downs;
  delete[] ptBinsAll_downs;
  
  delete[] ptBinsDY1_stranges;
  delete[] ptBinsDY16_stranges;
  delete[] ptBinsDY2_stranges;
  delete[] ptBinsDY3_stranges;
  delete[] ptBinsDY4_stranges;
  delete[] ptBinsDY5_stranges;
  delete[] ptBinsAll_stranges;
  
  delete[] ptBinsDY1_anti_ups;
  delete[] ptBinsDY16_anti_ups;
  delete[] ptBinsDY2_anti_ups;
  delete[] ptBinsDY3_anti_ups;
  delete[] ptBinsDY4_anti_ups;
  delete[] ptBinsDY5_anti_ups;
  delete[] ptBinsAll_anti_ups;
  
  delete[] ptBinsDY1_anti_downs;
  delete[] ptBinsDY16_anti_downs;
  delete[] ptBinsDY2_anti_downs;
  delete[] ptBinsDY3_anti_downs;
  delete[] ptBinsDY4_anti_downs;
  delete[] ptBinsDY5_anti_downs;
  delete[] ptBinsAll_anti_downs;
  
  delete[] ptBinsDY1_anti_stranges;
  delete[] ptBinsDY16_anti_stranges;
  delete[] ptBinsDY2_anti_stranges;
  delete[] ptBinsDY3_anti_stranges;
  delete[] ptBinsDY4_anti_stranges;
  delete[] ptBinsDY5_anti_stranges;
  delete[] ptBinsAll_anti_stranges;
  
  delete[] ptBinsDY1_all;
  delete[] ptBinsDY16_all;
  delete[] ptBinsDY2_all;
  delete[] ptBinsDY3_all;
  delete[] ptBinsDY4_all;
  delete[] ptBinsDY5_all;
  delete[] ptBinsAll_all;
  
  delete[] ptBinsSoftAll_gluons;
  delete[] ptBinsSoftAll_quarks;
  delete[] ptBinsSoftAll_ups;
  delete[] ptBinsSoftAll_downs;
  delete[] ptBinsSoftAll_stranges;
  delete[] ptBinsSoftAll_anti_ups;
  delete[] ptBinsSoftAll_anti_downs;
  delete[] ptBinsSoftAll_anti_stranges;
  delete[] ptBinsSoftAll_all;
  
  delete[] mean_PT;
  delete[] integrated_v2;
  delete[] integrated_v2_counts;
  delete[] differential_v2_values;
  delete[] differential_v2_counts;
  
  delete[] charmNmb;
  delete[] charmNmb_eta20;
  delete[] charmNmb_eta15;
  delete[] charmNmb_eta10;
  delete[] charmNmb_eta075;
  delete[] charmNmb_eta05;
  delete[] charmNmb_strap05;
  delete[] charmNmb_eta035;
  delete[] bottomNmb;
  delete[] bottomNmb_eta20;
  delete[] bottomNmb_eta15;
  delete[] bottomNmb_eta10;
  delete[] bottomNmb_eta075;
  delete[] bottomNmb_eta05;
  delete[] bottomNmb_eta035;
}


void analysis::ChargedHadronCorrelations(double time)
{
  
  double pt,px,py,pz,x,y;
  double phi;
  
  double rB0ref = 0.0;
  double iB0ref = 0.0;
  double rB1ref = 0.0;
  double iB1ref = 0.0;
  double rB2ref = 0.0;
  double iB2ref = 0.0;
  double rB3ref = 0.0;
  double iB3ref = 0.0;
  double rB4ref = 0.0;
  double iB4ref = 0.0;
  
  double Testparticles = theConfig->getTestparticles();
  
  double pTReferenceMin = theConfig->getPTRefMin();
  double pTReferenceMax = theConfig->getPTRefMax();
  
  stringstream ss;
  ss << time;
  //creates filename, for example: "./output/run32_corr_time0.5.dat"
  string Fname = "HerwigChHadrons";
  
  string filename = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + "_corr_" + Fname + "_time" + ss.str() + ".dat";
  string filenameCumulant = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + "_cumulant_" + Fname + "_time" + ss.str() + ".dat";
  //cout << "Do the correlation analysis with " << Testparticles << " Testparticles." << endl;
  //cout << particles.size() << endl;
  //cout << particles.size()/Testparticles << endl;
  //cout << "          file: " << filename << endl;  
  
  const double minpT=0.0;
  const double maxpT=8.0;
  const int numberBins = 16;
  
  binningValues rB0(filename, minpT, maxpT, numberBins);
  binningValues iB0(filename, minpT, maxpT, numberBins);
  binningValues rB1(filename, minpT, maxpT, numberBins);
  binningValues iB1(filename, minpT, maxpT, numberBins);
  binningValues rB2(filename, minpT, maxpT, numberBins);
  binningValues iB2(filename, minpT, maxpT, numberBins);
  binningValues rB3(filename, minpT, maxpT, numberBins);
  binningValues iB3(filename, minpT, maxpT, numberBins);
  binningValues rB4(filename, minpT, maxpT, numberBins);
  binningValues iB4(filename, minpT, maxpT, numberBins);
  
  const double minEta = 0.0;
  const double maxEta = 2.0;
  const int numberBinsEta = 20;  
  
  binningValues rB0RapidityDifferential(filenameCumulant,minEta,maxEta,numberBinsEta);
  binningValues rB2RapidityDifferential(filenameCumulant,minEta,maxEta,numberBinsEta);
  binningValues iB2RapidityDifferential(filenameCumulant,minEta,maxEta,numberBinsEta);
  
  double rB0Cumulant=0.;
  double rB2Cumulant=0.;
  double iB2Cumulant=0.;
  
  double r = 0.;
  double offsetX = 0.;
  double offsetY = 0.;
  double eccentricity2 = 0.;
  double angleEccentricity = 0.;
  double eccentricity3 = 0.;
  double eccentricity4 = 0.;
  
  double avrN1 = 0.;
  double avrN2 = 0.;
  double avrN3 = 0.;
  double avrN4 = 0.;
  
  double avcosN1r1 = 0.;
  double avsinN1r1 = 0.;
  double avcosN1 = 0.;
  double avsinN1 = 0.;
  
  double avcosN2 = 0.;
  double avsinN2 = 0.;
  
  double avcosN3 = 0.;
  double avsinN3 = 0.;  

  double avcosN4 = 0.;
  double avsinN4 = 0.;  
  
  //real and imaginary parts of eN= Sum_{i}    r[i]^n Exp[ I n phi[i] ] / Sum_{i}  r[i]^n
  double ReEN1 = 0.;
  double ImEN1 = 0.;
  double ReEN1r1 =  0.;
  double ImEN1r1 =  0.;
  double ReEN2 =  0.;
  double ImEN2 =  0.;  
  double ReEN3 = 0.;
  double ImEN3 = 0.;
  double ReEN4 = 0.;
  double ImEN4 = 0.;
  
  //energy weighted
  double avcosN1_energyweight = 0.;
  double avsinN1_energyweight = 0.;
  double avrN1_energyweight = 0.;
  
  double avcosN1r1_energyweight = 0.;
  double avsinN1r1_energyweight = 0.;
  
  double avcosN2_energyweight = 0.;
  double avsinN2_energyweight = 0.;
  double avrN2_energyweight = 0.;
  
  double avrN3_energyweight = 0.;
  double avsinN3_energyweight = 0.;
  double avcosN3_energyweight = 0.;
  
  double avrN4_energyweight = 0.;
  double avsinN4_energyweight = 0.;
  double avcosN4_energyweight = 0.;
  
  //energy weighted real and imaginary parts of eN= Sum_{i}   e[i] r[i]^n Exp[ I n phi[i] ] / Sum_{i}  e[i] r[i]^n
  double ReEN1_energyweight =  0.;
  double ImEN1_energyweight =  0.;
  double ReEN1r1_energyweight =  0.;
  double ImEN1r1_energyweight =  0.;
  double ReEN2_energyweight =  0.;
  double ImEN2_energyweight =  0.;
  
  double phiSPACE=0.;
  
  //Get the center for this flavor
  int NFromFlavor=0; 
  for(unsigned int i=0;i< hadrons.size(); i++ )
  {    
    //Space anisotropy
    x = hadrons[i].Pos.X();
    y = hadrons[i].Pos.Y();
    offsetX += x;
    offsetY += y;
    NFromFlavor++;
  }
  offsetX /= hadrons.size();
  offsetY /= hadrons.size();
//   cout << Fname + " OFFSET X = " << offsetX << endl;
//   cout << Fname + " OFFSET Y = " << offsetY << endl;
  
  //Get the correlations
  for ( unsigned int i = 0; i < hadrons.size(); i++ )
  {
    double eta = hadrons[i].Mom.Rapidity();
    
    //Cumulants, differential in Delta-eta
    {
      //Momentum anisotropy
      px = hadrons[i].Mom.Px(); 
      py = hadrons[i].Mom.Py();
      pz = hadrons[i].Mom.Pz();
      pt = hadrons[i].Mom.Pt();  
      phi = getAngle2pi(px,py); 
      if( (minEta < eta) && (eta < maxEta) )
      { 
        rB0Cumulant+=1.0;
        rB2Cumulant+=cos(2.0*phi);
        iB2Cumulant+=sin(2.0*phi);
        //B0       
        rB0RapidityDifferential.add(abs(eta),1.0);           
        //B2
        rB2RapidityDifferential.add(abs(eta),cos(2.0*phi));
        iB2RapidityDifferential.add(abs(eta),sin(2.0*phi));    
      }      
    }
    //V2 Analyse
    //for ( int yRangeIndex = 0; yRangeIndex < rapidityRanges.size(); yRangeIndex++ )
    {
      int yRangeIndex = 0;//Midrapidity
      if ( fabs( eta ) >= rapidityRanges[yRangeIndex].yleft && fabs( eta ) <= rapidityRanges[yRangeIndex].yright )
      {
        //Momentum anisotropy
        px = hadrons[i].Mom.Px(); 
        py = hadrons[i].Mom.Py();
        pz = hadrons[i].Mom.Pz();
        pt = hadrons[i].Mom.Pt();  
        phi = getAngle2pi(px,py);
        
        x = hadrons[i].Pos.X()-offsetX,
        y = hadrons[i].Pos.Y()-offsetY;
        
        r = sqrt(pow(x,2.0)+pow(y,2.0));
        phiSPACE =  getAngle2pi(x,y);
        
        if( (pTReferenceMin < pt) && (pt < pTReferenceMax) )
        {
          //B0ref
          rB0ref +=  1.0;
          iB0ref +=  0.0;
          //B1ref
          rB1ref += cos(phi);
          iB1ref += sin(phi);      
          //B2ref
          rB2ref += cos(2.0*phi);
          iB2ref += sin(2.0*phi);  
          //B3ref
          rB3ref += cos(3.0*phi);
          iB3ref += sin(3.0*phi); 
          //B4ref
          rB4ref += cos(4.0*phi);
          iB4ref += sin(4.0*phi); 
        }
        
        //pT dependent
        rB0.add(pt,1.0);
        iB0.add(pt,0.0);
        
        rB1.add(pt,cos(phi));
        iB1.add(pt,sin(phi));
        
        rB2.add(pt,cos(2.0*phi));
        iB2.add(pt,sin(2.0*phi));
        
        rB3.add(pt,cos(3.0*phi));
        iB3.add(pt,sin(3.0*phi));

        rB4.add(pt,cos(4.0*phi));
        iB4.add(pt,sin(4.0*phi));
        
        //Eccentricities
        avcosN1 += pow(r,2.0)*cos(1.0*phiSPACE);
        avsinN1 += pow(r,2.0)*sin(1.0*phiSPACE); 
        
        avrN1 += r;
        avcosN1r1 += pow(r,1.0)*cos(1.0*phiSPACE);
        avsinN1r1 += pow(r,1.0)*sin(1.0*phiSPACE);  
        
        avrN2 += pow(r,2.0);
        avcosN2 += pow(r,2.0)*cos(2.0*phiSPACE);
        avsinN2 += pow(r,2.0)*sin(2.0*phiSPACE);
        
        avrN3 += pow(r,3.0);
        avcosN3 += pow(r,3.0)*cos(3.0*phiSPACE); 
        avsinN3 += pow(r,3.0)*sin(3.0*phiSPACE);
        
        avrN4 += pow(r,4.0);
        avcosN4 += pow(r,4.0)*cos(4.0*phiSPACE); 
        avsinN4 += pow(r,4.0)*sin(4.0*phiSPACE);  
        
        avcosN1_energyweight += pt * pow(r,2.0)*cos(1.0*phiSPACE);
        avsinN1_energyweight += pt * pow(r,2.0)*sin(1.0*phiSPACE); 
        
        avrN1_energyweight += pt * r;
        avcosN1r1_energyweight += pt * pow(r,1.0)*cos(1.0*phiSPACE);
        avsinN1r1_energyweight += pt * pow(r,1.0)*sin(1.0*phiSPACE);  
        
        avrN2_energyweight += pt * pow(r,2.0);
        avcosN2_energyweight += pt * pow(r,2.0)*cos(2.0*phiSPACE);
        avsinN2_energyweight += pt * pow(r,2.0)*sin(2.0*phiSPACE);          
        
        avcosN3_energyweight += pt * pow(r,3.0)*cos(3.0*phiSPACE); 
        avsinN3_energyweight += pt * pow(r,3.0)*sin(3.0*phiSPACE);          
        
        avcosN4_energyweight += pt * pow(r,4.0)*cos(4.0*phiSPACE); 
        avsinN4_energyweight += pt * pow(r,4.0)*sin(4.0*phiSPACE);        
      }
    }
  }  
  
  
  //Normalize
  avrN1 /= rB0ref;
  avrN2 /= rB0ref;
  avrN3 /= rB0ref;
  avcosN1 /= rB0ref; 
  avsinN1 /= rB0ref; 
  avcosN2 /= rB0ref; 
  avsinN2 /= rB0ref; 
  avcosN3 /= rB0ref; 
  avsinN3 /= rB0ref;
  avcosN4 /= rB0ref; 
  avsinN4 /= rB0ref; 
  
  avrN1_energyweight /= rB0ref;
  avrN2_energyweight /= rB0ref;
  avrN3_energyweight /= rB0ref;
  avcosN1_energyweight /= rB0ref; 
  avsinN1_energyweight /= rB0ref; 
  avcosN2_energyweight /= rB0ref; 
  avsinN2_energyweight /= rB0ref; 
  avcosN3_energyweight /= rB0ref; 
  avsinN3_energyweight /= rB0ref;    
  avcosN4_energyweight /= rB0ref; 
  avsinN4_energyweight /= rB0ref;    
  
  eccentricity2 = sqrt(pow(avcosN2,2.0)+pow(avsinN2,2.0))/avrN2;
  eccentricity3 = sqrt(pow(avcosN3,2.0)+pow(avsinN3,2.0))/avrN3; 
  eccentricity4 = sqrt(pow(avcosN4,2.0)+pow(avsinN4,2.0))/avrN4; 

  angleEccentricity = 1./2. * (atan(iB2ref/rB2ref) + M_PI);
  
  //Print out the correlations for vnAnalysis
  {
    fstream file(filename.c_str(),ios::out);
    file  << "#REFERENCE " 
          << rB0ref/Testparticles
          << ' ' 
          << iB0ref/Testparticles 
          << ' ' 
          << rB1ref/Testparticles 
          << ' ' 
          << iB1ref/Testparticles
          << ' ' 
          << rB2ref/Testparticles 
          << ' ' 
          << iB2ref/Testparticles
          << ' ' 
          << rB3ref/Testparticles 
          << ' ' 
          << iB3ref/Testparticles
          << ' ' 
          << rB4ref/Testparticles 
          << ' ' 
          << iB4ref/Testparticles
          << ' '
          << endl;
          
          
    for (int i=0; i<rB0.getNBins(); i++)
    {
      file  << rB0.getBinLabel(i) 
            << "\t" << rB0.getBinValue(i)/Testparticles 
            << "\t" << iB0.getBinValue(i)/Testparticles  
            << "\t" << rB1.getBinValue(i)/Testparticles 
            << "\t" << iB1.getBinValue(i)/Testparticles 
            << "\t" << rB2.getBinValue(i)/Testparticles 
            << "\t" << iB2.getBinValue(i)/Testparticles
            << "\t" << rB3.getBinValue(i)/Testparticles 
            << "\t" << iB3.getBinValue(i)/Testparticles
            << "\t" << rB4.getBinValue(i)/Testparticles 
            << "\t" << iB4.getBinValue(i)/Testparticles
            << endl;
    }    
    file << "#ECCENTRICITY"     //1
          << ' '
          << eccentricity2       //2
          << ' '
          << eccentricity3       //3
          << ' '
          << eccentricity4       //4 new
          << ' '
          << angleEccentricity   //5
          << ' '
          << avcosN1r1           //6
          << ' '
          << avsinN1r1           //7
          << ' '
          << avcosN1             //8
          << ' '
          << avsinN1             //9
          << ' '
          << avrN1               //10
          << ' '
          << avcosN2             //11
          << ' '
          << avsinN2             //12
          << ' '       
          << avrN2               //13
          << ' '
          << avcosN3             //14
          << ' '
          << avsinN3             //15
          << ' '
          << avrN3               //16
          << ' '
          << avcosN4             //17
          << ' '
          << avsinN4             //18
          << ' '
          << avrN4               //19
          << ' '
          << avcosN1r1_energyweight  //20
          << ' '
          << avsinN1r1_energyweight  //21
          << ' '
          << avcosN1_energyweight    //22
          << ' '
          << avsinN1_energyweight    //23
          << ' '
          << avrN1_energyweight      //24
          << ' '
          << avcosN2_energyweight    //25
          << ' '
          << avsinN2_energyweight    //26
          << ' '       
          << avrN2_energyweight      //27
          << ' '
          << avcosN3_energyweight    //28
          << ' '
          << avsinN3_energyweight    //29
          << ' '
          << avrN3_energyweight      //30
          << ' '
          << avcosN4_energyweight    //31
          << ' '
          << avsinN4_energyweight    //32
          << ' '
          << avrN4_energyweight      //33
          << ' '
          << rB0ref                  //34
          << endl;
      
    file.close();   
  }
}




void analysis::ParticleCorrelations(double time, FLAVOR_TYPE thisFlavor)
{
  
  double pt,px,py,pz,x,y;
  double phi;
  
  double rB0ref = 0.0;
  double iB0ref = 0.0;
  double rB1ref = 0.0;
  double iB1ref = 0.0;
  double rB2ref = 0.0;
  double iB2ref = 0.0;
  double rB3ref = 0.0;
  double iB3ref = 0.0;
  double rB4ref = 0.0;
  double iB4ref = 0.0;
  
  double Testparticles = theConfig->getTestparticles();
  
  double pTReferenceMin = theConfig->getPTRefMin();
  double pTReferenceMax = theConfig->getPTRefMax();
  
  stringstream ss;
  ss << time;
  //creates filename, for example: "./output/run32_corr_time0.5.dat"
  string Fname = ParticlePrototype::getName(thisFlavor);
  
  string filename = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + "_corr_" + Fname + "_time" + ss.str() + ".dat";
  string filenameCumulant = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + "_cumulant_" + Fname + "_time" + ss.str() + ".dat";
  //cout << "Do the correlation analysis with " << Testparticles << " Testparticles." << endl;
  //cout << particles.size() << endl;
  //cout << particles.size()/Testparticles << endl;
  //cout << "          file: " << filename << endl;  
  
  const double minpT=0.0;
  const double maxpT=8.0;
  const int numberBins = 16;
  
  binningValues rB0(filename, minpT, maxpT, numberBins);
  binningValues iB0(filename, minpT, maxpT, numberBins);
  binningValues rB1(filename, minpT, maxpT, numberBins);
  binningValues iB1(filename, minpT, maxpT, numberBins);
  binningValues rB2(filename, minpT, maxpT, numberBins);
  binningValues iB2(filename, minpT, maxpT, numberBins);
  binningValues rB3(filename, minpT, maxpT, numberBins);
  binningValues iB3(filename, minpT, maxpT, numberBins);
  binningValues rB4(filename, minpT, maxpT, numberBins);
  binningValues iB4(filename, minpT, maxpT, numberBins);
  
  const double minEta = 0.0;
  const double maxEta = 2.0;
  const int numberBinsEta = 20;  
  
  binningValues rB0RapidityDifferential(filenameCumulant,minEta,maxEta,numberBinsEta);
  binningValues rB2RapidityDifferential(filenameCumulant,minEta,maxEta,numberBinsEta);
  binningValues iB2RapidityDifferential(filenameCumulant,minEta,maxEta,numberBinsEta);
  
  double rB0Cumulant=0.;
  double rB2Cumulant=0.;
  double iB2Cumulant=0.;
  
  double r = 0.;
  double offsetX = 0.;
  double offsetY = 0.;
  double eccentricity2 = 0.;
  double angleEccentricity = 0.;
  double eccentricity3 = 0.;
  double eccentricity4 = 0.;
  
  double avrN1 = 0.;
  double avrN2 = 0.;
  double avrN3 = 0.;
  double avrN4 = 0.;
  
  double avcosN1r1 = 0.;
  double avsinN1r1 = 0.;
  double avcosN1 = 0.;
  double avsinN1 = 0.;
  
  double avcosN2 = 0.;
  double avsinN2 = 0.;
  
  double avcosN3 = 0.;
  double avsinN3 = 0.;  

  double avcosN4 = 0.;
  double avsinN4 = 0.;  
  
  //real and imaginary parts of eN= Sum_{i}    r[i]^n Exp[ I n phi[i] ] / Sum_{i}  r[i]^n
  double ReEN1 = 0.;
  double ImEN1 = 0.;
  double ReEN1r1 =  0.;
  double ImEN1r1 =  0.;
  double ReEN2 =  0.;
  double ImEN2 =  0.;  
  double ReEN3 = 0.;
  double ImEN3 = 0.;
  double ReEN4 = 0.;
  double ImEN4 = 0.;
  
  //energy weighted
  double avcosN1_energyweight = 0.;
  double avsinN1_energyweight = 0.;
  double avrN1_energyweight = 0.;
  
  double avcosN1r1_energyweight = 0.;
  double avsinN1r1_energyweight = 0.;
  
  double avcosN2_energyweight = 0.;
  double avsinN2_energyweight = 0.;
  double avrN2_energyweight = 0.;
  
  double avrN3_energyweight = 0.;
  double avsinN3_energyweight = 0.;
  double avcosN3_energyweight = 0.;
  
  double avrN4_energyweight = 0.;
  double avsinN4_energyweight = 0.;
  double avcosN4_energyweight = 0.;
  
  //energy weighted real and imaginary parts of eN= Sum_{i}   e[i] r[i]^n Exp[ I n phi[i] ] / Sum_{i}  e[i] r[i]^n
  double ReEN1_energyweight =  0.;
  double ImEN1_energyweight =  0.;
  double ReEN1r1_energyweight =  0.;
  double ImEN1r1_energyweight =  0.;
  double ReEN2_energyweight =  0.;
  double ImEN2_energyweight =  0.;
  
  double phiSPACE=0.;
  
  //Get the center for this flavor
  int NFromFlavor=0; 
  for(unsigned int i=0;i< particles.size(); i++ )
  {
    if(ParticlePrototype::mapQuarkGluon( particles[i].FLAVOR ) != thisFlavor)
    {
      continue;
    }      
    //Space anisotropy
    x = particles[i].Pos.X();
    y = particles[i].Pos.Y();
    offsetX += x;
    offsetY += y;
    NFromFlavor++;
  }
  offsetX /= particles.size();
  offsetY /= particles.size();
//   cout << Fname + " OFFSET X = " << offsetX << endl;
//   cout << Fname + " OFFSET Y = " << offsetY << endl;
  
  //Get the correlations
  for ( unsigned int i = 0; i < particles.size(); i++ )
  {
    if(ParticlePrototype::mapQuarkGluon(particles[i].FLAVOR) != thisFlavor)
    {
      continue;
    }
    double eta = particles[i].Mom.Rapidity();
    
    //Cumulants, differential in Delta-eta
    {
      //Momentum anisotropy
      px = particles[i].Mom.Px(); 
      py = particles[i].Mom.Py();
      pz = particles[i].Mom.Pz();
      pt = particles[i].Mom.Pt();  
      phi = getAngle2pi(px,py); 
      if( (minEta < eta) && (eta < maxEta) )
      { 
        rB0Cumulant+=1.0;
        rB2Cumulant+=cos(2.0*phi);
        iB2Cumulant+=sin(2.0*phi);
        //B0       
        rB0RapidityDifferential.add(abs(eta),1.0);           
        //B2
        rB2RapidityDifferential.add(abs(eta),cos(2.0*phi));
        iB2RapidityDifferential.add(abs(eta),sin(2.0*phi));    
      }      
    }
    //V2 Analyse
    //for ( int yRangeIndex = 0; yRangeIndex < rapidityRanges.size(); yRangeIndex++ )
    {
      int yRangeIndex = 0;//Midrapidity
      if ( fabs( eta ) >= rapidityRanges[yRangeIndex].yleft && fabs( eta ) <= rapidityRanges[yRangeIndex].yright )
      {
        //Momentum anisotropy
        px = particles[i].Mom.Px(); 
        py = particles[i].Mom.Py();
        pz = particles[i].Mom.Pz();
        pt = particles[i].Mom.Pt();  
        phi = getAngle2pi(px,py);
        
        x = particles[i].Pos.X()-offsetX,
        y = particles[i].Pos.Y()-offsetY;
        
        r = sqrt(pow(x,2.0)+pow(y,2.0));
        phiSPACE =  getAngle2pi(x,y);
        
        if( (pTReferenceMin < pt) && (pt < pTReferenceMax) )
        {
          //B0ref
          rB0ref +=  1.0;
          iB0ref +=  0.0;
          //B1ref
          rB1ref += cos(phi);
          iB1ref += sin(phi);      
          //B2ref
          rB2ref += cos(2.0*phi);
          iB2ref += sin(2.0*phi);  
          //B3ref
          rB3ref += cos(3.0*phi);
          iB3ref += sin(3.0*phi); 
          //B4ref
          rB4ref += cos(4.0*phi);
          iB4ref += sin(4.0*phi); 
        }
        
        //pT dependent
        rB0.add(pt,1.0);
        iB0.add(pt,0.0);
        
        rB1.add(pt,cos(phi));
        iB1.add(pt,sin(phi));
        
        rB2.add(pt,cos(2.0*phi));
        iB2.add(pt,sin(2.0*phi));
        
        rB3.add(pt,cos(3.0*phi));
        iB3.add(pt,sin(3.0*phi));

        rB4.add(pt,cos(4.0*phi));
        iB4.add(pt,sin(4.0*phi));
        
        //Eccentricities
        avcosN1 += pow(r,2.0)*cos(1.0*phiSPACE);
        avsinN1 += pow(r,2.0)*sin(1.0*phiSPACE); 
        
        avrN1 += r;
        avcosN1r1 += pow(r,1.0)*cos(1.0*phiSPACE);
        avsinN1r1 += pow(r,1.0)*sin(1.0*phiSPACE);  
        
        avrN2 += pow(r,2.0);
        avcosN2 += pow(r,2.0)*cos(2.0*phiSPACE);
        avsinN2 += pow(r,2.0)*sin(2.0*phiSPACE);
        
        avrN3 += pow(r,3.0);
        avcosN3 += pow(r,3.0)*cos(3.0*phiSPACE); 
        avsinN3 += pow(r,3.0)*sin(3.0*phiSPACE);
        
        avrN4 += pow(r,4.0);
        avcosN4 += pow(r,4.0)*cos(4.0*phiSPACE); 
        avsinN4 += pow(r,4.0)*sin(4.0*phiSPACE);  
        
        avcosN1_energyweight += pt * pow(r,2.0)*cos(1.0*phiSPACE);
        avsinN1_energyweight += pt * pow(r,2.0)*sin(1.0*phiSPACE); 
        
        avrN1_energyweight += pt * r;
        avcosN1r1_energyweight += pt * pow(r,1.0)*cos(1.0*phiSPACE);
        avsinN1r1_energyweight += pt * pow(r,1.0)*sin(1.0*phiSPACE);  
        
        avrN2_energyweight += pt * pow(r,2.0);
        avcosN2_energyweight += pt * pow(r,2.0)*cos(2.0*phiSPACE);
        avsinN2_energyweight += pt * pow(r,2.0)*sin(2.0*phiSPACE);          
        
        avcosN3_energyweight += pt * pow(r,3.0)*cos(3.0*phiSPACE); 
        avsinN3_energyweight += pt * pow(r,3.0)*sin(3.0*phiSPACE);          
        
        avcosN4_energyweight += pt * pow(r,4.0)*cos(4.0*phiSPACE); 
        avsinN4_energyweight += pt * pow(r,4.0)*sin(4.0*phiSPACE);        
      }
    }
  }  
  
  
  //Normalize
  avrN1 /= rB0ref;
  avrN2 /= rB0ref;
  avrN3 /= rB0ref;
  avcosN1 /= rB0ref; 
  avsinN1 /= rB0ref; 
  avcosN2 /= rB0ref; 
  avsinN2 /= rB0ref; 
  avcosN3 /= rB0ref; 
  avsinN3 /= rB0ref;
  avcosN4 /= rB0ref; 
  avsinN4 /= rB0ref; 
  
  avrN1_energyweight /= rB0ref;
  avrN2_energyweight /= rB0ref;
  avrN3_energyweight /= rB0ref;
  avcosN1_energyweight /= rB0ref; 
  avsinN1_energyweight /= rB0ref; 
  avcosN2_energyweight /= rB0ref; 
  avsinN2_energyweight /= rB0ref; 
  avcosN3_energyweight /= rB0ref; 
  avsinN3_energyweight /= rB0ref;    
  avcosN4_energyweight /= rB0ref; 
  avsinN4_energyweight /= rB0ref;    
  
  eccentricity2 = sqrt(pow(avcosN2,2.0)+pow(avsinN2,2.0))/avrN2;
  eccentricity3 = sqrt(pow(avcosN3,2.0)+pow(avsinN3,2.0))/avrN3; 
  eccentricity4 = sqrt(pow(avcosN4,2.0)+pow(avsinN4,2.0))/avrN4; 

  angleEccentricity = 1./2. * (atan(iB2ref/rB2ref) + M_PI);
  
  //Print out the correlations for vnAnalysis
  {
    fstream file(filename.c_str(),ios::out);
    file  << "#REFERENCE " 
          << rB0ref/Testparticles
          << ' ' 
          << iB0ref/Testparticles 
          << ' ' 
          << rB1ref/Testparticles 
          << ' ' 
          << iB1ref/Testparticles
          << ' ' 
          << rB2ref/Testparticles 
          << ' ' 
          << iB2ref/Testparticles
          << ' ' 
          << rB3ref/Testparticles 
          << ' ' 
          << iB3ref/Testparticles
          << ' ' 
          << rB4ref/Testparticles 
          << ' ' 
          << iB4ref/Testparticles
          << ' '
          << endl;
          
          
    for (int i=0; i<rB0.getNBins(); i++)
    {
      file  << rB0.getBinLabel(i) 
            << "\t" << rB0.getBinValue(i)/Testparticles 
            << "\t" << iB0.getBinValue(i)/Testparticles  
            << "\t" << rB1.getBinValue(i)/Testparticles 
            << "\t" << iB1.getBinValue(i)/Testparticles 
            << "\t" << rB2.getBinValue(i)/Testparticles 
            << "\t" << iB2.getBinValue(i)/Testparticles
            << "\t" << rB3.getBinValue(i)/Testparticles 
            << "\t" << iB3.getBinValue(i)/Testparticles
            << "\t" << rB4.getBinValue(i)/Testparticles 
            << "\t" << iB4.getBinValue(i)/Testparticles
            << endl;
    }    
    file << "#ECCENTRICITY"     //1
          << ' '
          << eccentricity2       //2
          << ' '
          << eccentricity3       //3
          << ' '
          << eccentricity4       //4 new
          << ' '
          << angleEccentricity   //5
          << ' '
          << avcosN1r1           //6
          << ' '
          << avsinN1r1           //7
          << ' '
          << avcosN1             //8
          << ' '
          << avsinN1             //9
          << ' '
          << avrN1               //10
          << ' '
          << avcosN2             //11
          << ' '
          << avsinN2             //12
          << ' '       
          << avrN2               //13
          << ' '
          << avcosN3             //14
          << ' '
          << avsinN3             //15
          << ' '
          << avrN3               //16
          << ' '
          << avcosN4             //17
          << ' '
          << avsinN4             //18
          << ' '
          << avrN4               //19
          << ' '
          << avcosN1r1_energyweight  //20
          << ' '
          << avsinN1r1_energyweight  //21
          << ' '
          << avcosN1_energyweight    //22
          << ' '
          << avsinN1_energyweight    //23
          << ' '
          << avrN1_energyweight      //24
          << ' '
          << avcosN2_energyweight    //25
          << ' '
          << avsinN2_energyweight    //26
          << ' '       
          << avrN2_energyweight      //27
          << ' '
          << avcosN3_energyweight    //28
          << ' '
          << avsinN3_energyweight    //29
          << ' '
          << avrN3_energyweight      //30
          << ' '
          << avcosN4_energyweight    //31
          << ' '
          << avsinN4_energyweight    //32
          << ' '
          << avrN4_energyweight      //33
          << ' '
          << rB0ref                  //34
          << endl;
      
    file.close();   
  }
  
  //TEST
  //Print out the correlations for DeltaEta-differential cumulant analysis
//   {
//     fstream file2(filenameCumulant.c_str(),ios::out);
//     file2  << "#REFERENCE " 
//           << rB0Cumulant/Testparticles
//           << ' ' 
//           << rB2Cumulant/Testparticles 
//           << ' ' 
//           << iB2Cumulant/Testparticles   
//           << endl;
//     for (int i=0; i<rB0RapidityDifferential.getNBins(); i++)
//     {
//       file2   << rB0RapidityDifferential.getBinLabel(i) 
//               << "\t" << rB0RapidityDifferential.getBinValue(i)/Testparticles 
//               << "\t" << rB2RapidityDifferential.getBinValue(i)/Testparticles  
//               << "\t" << iB2RapidityDifferential.getBinValue(i)/Testparticles 
//               << endl;
//     }  
//     file2.close(); 
//   }
  

//   if(SquaredOn)
//   {
//     for ( unsigned int i = 0; i < particles.size(); i++ )
//     {
//       double eta = particles[i].Mom.Rapidity();
//       //for ( int yRangeIndex = 0; yRangeIndex < rapidityRanges.size(); yRangeIndex++ )
//       {
//         int yRangeIndex = 0;//Midrapidity
//         if ( fabs( eta ) >= rapidityRanges[yRangeIndex].yleft && fabs( eta ) <= rapidityRanges[yRangeIndex].yright )
//         {
//           px = particles[i].Mom.Px(); 
//           py = particles[i].Mom.Py();
//           pz = particles[i].Mom.Pz();
//           pt = particles[i].Mom.Pt();  
//           phi = getAngle2pi(px,py);
//           //B0ref
//           rB0ref +=  1.0;
//           iB0ref +=  0.0;
//           rB0refsquare +=  1.0;
//           iB0refsquare +=  0.0;
//           //B1ref
//           rB1ref += cos(phi);
//           iB1ref += sin(phi);
//           rB1refsquare += cos(phi)*cos(phi);
//           iB1refsquare += sin(phi)*sin(phi);        
//           //B2ref
//           rB2ref += cos(2.0*phi);
//           iB2ref += sin(2.0*phi);  
//           rB2refsquare += cos(2.0*phi)*cos(2.0*phi);
//           iB2refsquare += sin(2.0*phi)*sin(2.0*phi);
//           
//           //pT dependent
//           rB0.add(pt,1.0);
//           iB0.add(pt,0.0);
//           
//           rB1.add(pt,cos(phi));
//           iB1.add(pt,sin(phi));
//           
//           rB2.add(pt,cos(2.0*phi));
//           iB2.add(pt,sin(2.0*phi));
//           
//           rB0square.add(pt,1.0);
//           iB0square.add(pt,0.0);
//           
//           rB1square.add(pt,cos(phi)*cos(phi));
//           iB1square.add(pt,sin(phi)*sin(phi));
//           
//           rB2square.add(pt,cos(2.0*phi)*cos(2.0*phi));
//           iB2square.add(pt,sin(2.0*phi)*sin(2.0*phi));        
//           
//         }
//       }
//     }  
//     
// 
//     fstream file(filename.c_str(),ios::out);
//     file  << "#REFERENCE " 
//           << rB0ref/Testparticles
//           << ' ' 
//           << rB0refsquare/Testparticles
//           << ' ' 
//           << iB0ref/Testparticles 
//           << ' ' 
//           << iB0refsquare/Testparticles
//           << ' ' 
//           << rB1ref/Testparticles 
//           << ' ' 
//           << rB1refsquare/Testparticles
//           << ' ' 
//           << iB1ref/Testparticles
//           << ' ' 
//           << iB1refsquare/Testparticles
//           << ' ' 
//           << rB2ref/Testparticles 
//           << ' '
//           << rB2refsquare/Testparticles
//           << ' ' 
//           << iB2ref/Testparticles
//           << ' ' 
//           << iB2refsquare/Testparticles
//           << endl;
//           
//           
//     for (int i=0; i<rB0.getNBins(); i++)
//     {
//       file  << rB0.getBinLabel(i) 
//             << "\t" << rB0.getBinValue(i)/Testparticles << "\t" << rB0square.getBinValue(i)/Testparticles 
//             << "\t" << iB0.getBinValue(i)/Testparticles << "\t" << iB0square.getBinValue(i)/Testparticles 
//             << "\t" << rB1.getBinValue(i)/Testparticles << "\t" << rB1square.getBinValue(i)/Testparticles 
//             << "\t" << iB1.getBinValue(i)/Testparticles << "\t" << iB1square.getBinValue(i)/Testparticles 
//             << "\t" << rB2.getBinValue(i)/Testparticles << "\t" << rB2square.getBinValue(i)/Testparticles 
//             << "\t" << iB2.getBinValue(i)/Testparticles  << "\t"<< iB2square.getBinValue(i)/Testparticles << endl;
//     }
//   }
//   else
//   {
//     //TODO: QUARK and GLUONS SEPARATELY
//     //Get the center
//     for(unsigned int i=0;i< particles.size(); i++ )
//     {
//       //Space anisotropy
//       x = particles[i].Pos.X();
//       y = particles[i].Pos.Y();
//       offsetX += x;
//       offsetY += y;
//     }
//     offsetX /= particles.size();
//     offsetY /= particles.size();
//     cout << "OFFSET X = " << offsetX << endl;
//     cout << "OFFSET Y = " << offsetY << endl;
//     for ( unsigned int i = 0; i < particles.size(); i++ )
//     {
//       double eta = particles[i].Mom.Rapidity();
//       //for ( int yRangeIndex = 0; yRangeIndex < rapidityRanges.size(); yRangeIndex++ )
//       {
//         int yRangeIndex = 0;//Midrapidity
//         if ( fabs( eta ) >= rapidityRanges[yRangeIndex].yleft && fabs( eta ) <= rapidityRanges[yRangeIndex].yright )
//         {
//           //Momentum anisotropy
//           px = particles[i].Mom.Px(); 
//           py = particles[i].Mom.Py();
//           pz = particles[i].Mom.Pz();
//           pt = particles[i].Mom.Pt();  
//           phi = getAngle2pi(px,py);
//           x = particles[i].Pos.X()-offsetX,
//           y = particles[i].Pos.Y()-offsetY;
//           r = sqrt(pow(x,2.0)+pow(y,2.0));
//           phiSPACE =  getAngle2pi(x,y);
//           
//           if( (pTReferenceMin < pt) && (pt < pTReferenceMax) )
//           {
//             //B0ref
//             rB0ref +=  1.0;
//             iB0ref +=  0.0;
//             //B1ref
//             rB1ref += cos(phi);
//             iB1ref += sin(phi);      
//             //B2ref
//             rB2ref += cos(2.0*phi);
//             iB2ref += sin(2.0*phi);  
//             //B3ref
//             rB3ref += cos(3.0*phi);
//             iB3ref += sin(3.0*phi); 
//             //B4ref
//             rB4ref += cos(4.0*phi);
//             iB4ref += sin(4.0*phi); 
//           }
//           
//           //pT dependent
//           rB0.add(pt,1.0);
//           iB0.add(pt,0.0);
//           
//           rB1.add(pt,cos(phi));
//           iB1.add(pt,sin(phi));
//           
//           rB2.add(pt,cos(2.0*phi));
//           iB2.add(pt,sin(2.0*phi));
//           
//           rB3.add(pt,cos(3.0*phi));
//           iB3.add(pt,sin(3.0*phi));
// 
//           rB4.add(pt,cos(4.0*phi));
//           iB4.add(pt,sin(4.0*phi));
//           
//           //Eccentricities
//           avcosN1 += pow(r,2.0)*cos(1.0*phiSPACE);
//           avsinN1 += pow(r,2.0)*sin(1.0*phiSPACE); 
//           
//           avrN1 += r;
//           avcosN1r1 += pow(r,1.0)*cos(1.0*phiSPACE);
//           avsinN1r1 += pow(r,1.0)*sin(1.0*phiSPACE);  
//           
//           avrN2 += pow(r,2.0);
//           avcosN2 += pow(r,2.0)*cos(2.0*phiSPACE);
//           avsinN2 += pow(r,2.0)*sin(2.0*phiSPACE);
//           
//           avrN3 += pow(r,3.0);
//           avcosN3 += pow(r,3.0)*cos(3.0*phiSPACE); 
//           avsinN3 += pow(r,3.0)*sin(3.0*phiSPACE);
//           
//           avrN4 += pow(r,4.0);
//           avcosN4 += pow(r,4.0)*cos(4.0*phiSPACE); 
//           avsinN4 += pow(r,4.0)*sin(4.0*phiSPACE);  
//           
//           avcosN1_energyweight += pt * pow(r,2.0)*cos(1.0*phiSPACE);
//           avsinN1_energyweight += pt * pow(r,2.0)*sin(1.0*phiSPACE); 
//           
//           avrN1_energyweight += pt * r;
//           avcosN1r1_energyweight += pt * pow(r,1.0)*cos(1.0*phiSPACE);
//           avsinN1r1_energyweight += pt * pow(r,1.0)*sin(1.0*phiSPACE);  
//           
//           avrN2_energyweight += pt * pow(r,2.0);
//           avcosN2_energyweight += pt * pow(r,2.0)*cos(2.0*phiSPACE);
//           avsinN2_energyweight += pt * pow(r,2.0)*sin(2.0*phiSPACE);          
//           
//           avcosN3_energyweight += pt * pow(r,3.0)*cos(3.0*phiSPACE); 
//           avsinN3_energyweight += pt * pow(r,3.0)*sin(3.0*phiSPACE);          
//           
//           avcosN4_energyweight += pt * pow(r,4.0)*cos(4.0*phiSPACE); 
//           avsinN4_energyweight += pt * pow(r,4.0)*sin(4.0*phiSPACE);        
//         }
//       }
//     }  
//     avrN1 /= rB0ref;
//     avrN2 /= rB0ref;
//     avrN3 /= rB0ref;
//     avcosN1 /= rB0ref; 
//     avsinN1 /= rB0ref; 
//     avcosN2 /= rB0ref; 
//     avsinN2 /= rB0ref; 
//     avcosN3 /= rB0ref; 
//     avsinN3 /= rB0ref;
//     avcosN4 /= rB0ref; 
//     avsinN4 /= rB0ref; 
//     
//     avrN1_energyweight /= rB0ref;
//     avrN2_energyweight /= rB0ref;
//     avrN3_energyweight /= rB0ref;
//     avcosN1_energyweight /= rB0ref; 
//     avsinN1_energyweight /= rB0ref; 
//     avcosN2_energyweight /= rB0ref; 
//     avsinN2_energyweight /= rB0ref; 
//     avcosN3_energyweight /= rB0ref; 
//     avsinN3_energyweight /= rB0ref;    
//     avcosN4_energyweight /= rB0ref; 
//     avsinN4_energyweight /= rB0ref;    
//     
//     eccentricity2 = sqrt(pow(avcosN2,2.0)+pow(avsinN2,2.0))/avrN2;
//     eccentricity3 = sqrt(pow(avcosN3,2.0)+pow(avsinN3,2.0))/avrN3; 
//     eccentricity4 = sqrt(pow(avcosN4,2.0)+pow(avsinN4,2.0))/avrN4; 
//     
//     //1/n 
//     angleEccentricity = 1./2. * (atan(iB2ref/rB2ref) + M_PI);
//     cout << "Eccentricity = " << eccentricity2 << endl;
//     fstream file(filename.c_str(),ios::out);
//     file  << "#REFERENCE " 
//           << rB0ref/Testparticles
//           << ' ' 
//           << iB0ref/Testparticles 
//           << ' ' 
//           << rB1ref/Testparticles 
//           << ' ' 
//           << iB1ref/Testparticles
//           << ' ' 
//           << rB2ref/Testparticles 
//           << ' ' 
//           << iB2ref/Testparticles
//           << ' ' 
//           << rB3ref/Testparticles 
//           << ' ' 
//           << iB3ref/Testparticles
//           << ' ' 
//           << rB4ref/Testparticles 
//           << ' ' 
//           << iB4ref/Testparticles
//           << ' '
//           << endl;
//           
//           
//     for (int i=0; i<rB0.getNBins(); i++)
//     {
//       file  << rB0.getBinLabel(i) 
//             << "\t" << rB0.getBinValue(i)/Testparticles 
//             << "\t" << iB0.getBinValue(i)/Testparticles  
//             << "\t" << rB1.getBinValue(i)/Testparticles 
//             << "\t" << iB1.getBinValue(i)/Testparticles 
//             << "\t" << rB2.getBinValue(i)/Testparticles 
//             << "\t" << iB2.getBinValue(i)/Testparticles
//             << "\t" << rB3.getBinValue(i)/Testparticles 
//             << "\t" << iB3.getBinValue(i)/Testparticles
//             << "\t" << rB4.getBinValue(i)/Testparticles 
//             << "\t" << iB4.getBinValue(i)/Testparticles
//             << endl;
//     }    
//     file << "#ECCENTRICITY"     //1
//          << ' '
//          << eccentricity2       //2
//          << ' '
//          << eccentricity3       //3
//          << ' '
//          << eccentricity4       //4 new
//          << ' '
//          << angleEccentricity   //5
//          << ' '
//          << avcosN1r1           //6
//          << ' '
//          << avsinN1r1           //7
//          << ' '
//          << avcosN1             //8
//          << ' '
//          << avsinN1             //9
//          << ' '
//          << avrN1               //10
//          << ' '
//          << avcosN2             //11
//          << ' '
//          << avsinN2             //12
//          << ' '       
//          << avrN2               //13
//          << ' '
//          << avcosN3             //14
//          << ' '
//          << avsinN3             //15
//          << ' '
//          << avrN3               //16
//          << ' '
//          << avcosN4             //17
//          << ' '
//          << avsinN4             //18
//          << ' '
//          << avrN4               //19
//          << ' '
//          << avcosN1r1_energyweight  //20
//          << ' '
//          << avsinN1r1_energyweight  //21
//          << ' '
//          << avcosN1_energyweight    //22
//          << ' '
//          << avsinN1_energyweight    //23
//          << ' '
//          << avrN1_energyweight      //24
//          << ' '
//          << avcosN2_energyweight    //25
//          << ' '
//          << avsinN2_energyweight    //26
//          << ' '       
//          << avrN2_energyweight      //27
//          << ' '
//          << avcosN3_energyweight    //28
//          << ' '
//          << avsinN3_energyweight    //29
//          << ' '
//          << avrN3_energyweight      //30
//          << ' '
//          << avcosN4_energyweight    //31
//          << ' '
//          << avsinN4_energyweight    //32
//          << ' '
//          << avrN4_energyweight      //33
//          << ' '
//          << rB0ref                  //34
//          << endl;
//       
//     file.close(); 
//   }
//   //cout << "done" << endl;
}




void analysis::FragmentedCorrelations(double time)
{

//   string filenameX = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + "_Ds.dat";
//   double z_loop=0.05;
//   fstream fileX(filenameX.c_str(),ios::out);
//   
//   while(z_loop<1)
//   {   
//     z_loop += 0.0001;
//    fileX << z_loop << "\t" << KKPFragmentation(7,1,z_loop,M_SQRT2,gluon) << "\t" << KKPFragmentation(7,1,z_loop,2.0,gluon) << "\t" << KKPFragmentation(7,1,z_loop,4.,gluon)  << endl;
// 
//   }
//   
//   fileX.close();
//   exit(0);

  double pt,px,py,pz,x,y;
  double phi;
  
  double TotalMultInAcceptance=0; //CMS: p_T^h>0.4 GeV
  double TotalMult =0;
  double TotalMultInAcceptanceGluon=0; //CMS: p_T^h>0.4 GeV
  double TotalMultGluon =0;  
  
  
  double rB0ref = 0.0;
  double iB0ref = 0.0;
  double rB1ref = 0.0;
  double iB1ref = 0.0;
  double rB2ref = 0.0;
  double iB2ref = 0.0;
  double rB3ref = 0.0;
  double iB3ref = 0.0;
  double rB4ref = 0.0;
  double iB4ref = 0.0;
  
  double Testparticles = theConfig->getTestparticles();
  
  double pTReferenceMin = theConfig->getPTRefMin();
  double pTReferenceMax = theConfig->getPTRefMax();
  
  double pTHadronReferenceMin = theConfig->getPTHRefMin(); // 0.3 GeV (CMS)
  double pTHadronReferenceMax = theConfig->getPTHRefMax(); // 3.0 GeV (CMS)
  
  
  stringstream ss;
  ss << time;
  //creates filename, for example: "./output/run32_corr_Hadron_time0.5.dat"
  string Fname = "Hadron" ; //ParticlePrototype::getName(thisFlavor);
  
  string filename = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + "_corr_" + Fname + "_time" + ss.str() + ".dat";
  
  //Hadronic pt bins
  const double minpT=0.0;
  const double maxpT=6.0;
  const int numberBins = 20;
  
  binningValues rB0(filename, minpT, maxpT, numberBins);
  binningValues iB0(filename, minpT, maxpT, numberBins);
  binningValues rB1(filename, minpT, maxpT, numberBins);
  binningValues iB1(filename, minpT, maxpT, numberBins);
  binningValues rB2(filename, minpT, maxpT, numberBins);
  binningValues iB2(filename, minpT, maxpT, numberBins);
  binningValues rB3(filename, minpT, maxpT, numberBins);
  binningValues iB3(filename, minpT, maxpT, numberBins);
  binningValues rB4(filename, minpT, maxpT, numberBins);
  binningValues iB4(filename, minpT, maxpT, numberBins);
  
  const double minEta = 0.0;
  const double maxEta = 2.0;
  const int numberBinsEta = 20;  
  
  double rB0Cumulant=0.;
  double rB2Cumulant=0.;
  double iB2Cumulant=0.;
  
  double r = 0.;
  double offsetX = 0.;
  double offsetY = 0.;
  double eccentricity2 = 0.;
  double angleEccentricity = 0.;
  double eccentricity3 = 0.;
  double eccentricity4 = 0.;
  
  double avrN1 = 0.;
  double avrN2 = 0.;
  double avrN3 = 0.;
  double avrN4 = 0.;
  
  double avcosN1r1 = 0.;
  double avsinN1r1 = 0.;
  double avcosN1 = 0.;
  double avsinN1 = 0.;
  
  double avcosN2 = 0.;
  double avsinN2 = 0.;
  
  double avcosN3 = 0.;
  double avsinN3 = 0.;  

  double avcosN4 = 0.;
  double avsinN4 = 0.;  
  
  //real and imaginary parts of eN= Sum_{i}    r[i]^n Exp[ I n phi[i] ] / Sum_{i}  r[i]^n
  double ReEN1 = 0.;
  double ImEN1 = 0.;
  double ReEN1r1 =  0.;
  double ImEN1r1 =  0.;
  double ReEN2 =  0.;
  double ImEN2 =  0.;  
  double ReEN3 = 0.;
  double ImEN3 = 0.;
  double ReEN4 = 0.;
  double ImEN4 = 0.;
  
  //energy weighted
  double avcosN1_energyweight = 0.;
  double avsinN1_energyweight = 0.;
  double avrN1_energyweight = 0.;
  
  double avcosN1r1_energyweight = 0.;
  double avsinN1r1_energyweight = 0.;
  
  double avcosN2_energyweight = 0.;
  double avsinN2_energyweight = 0.;
  double avrN2_energyweight = 0.;
  
  double avrN3_energyweight = 0.;
  double avsinN3_energyweight = 0.;
  double avcosN3_energyweight = 0.;
  
  double avrN4_energyweight = 0.;
  double avsinN4_energyweight = 0.;
  double avcosN4_energyweight = 0.;
  
  //energy weighted real and imaginary parts of eN= Sum_{i}   e[i] r[i]^n Exp[ I n phi[i] ] / Sum_{i}  e[i] r[i]^n
  double ReEN1_energyweight =  0.;
  double ImEN1_energyweight =  0.;
  double ReEN1r1_energyweight =  0.;
  double ImEN1r1_energyweight =  0.;
  double ReEN2_energyweight =  0.;
  double ImEN2_energyweight =  0.;
  
  double phiSPACE=0.;
  
  
  double zmin=0.05;
  double zmax=1.0;
  double dz = zmin/20.;
  
  double Q;
  double ptH,ptP;
  
  FLAVOR_TYPE SPECIES[2] = {gluon, quark}; 
  FLAVOR_TYPE thisFlavor;
  
  //Get the offset, ignore the species for that
  for(unsigned int i=0;i< particles.size(); i++ )
  {     
      /*if(ParticlePrototype::mapQuarkGluon( particles[i].FLAVOR ) != thisFlavor)
      {
        continue;
      }  */    
      //Space anisotropy
      x = particles[i].Pos.X();
      y = particles[i].Pos.Y();
      offsetX += x;
      offsetY += y;
  }
  offsetX /= particles.size();
  offsetY /= particles.size();
   
//   for(int flavNum = 0; flavNum<2;flavNum++)
//   {
//     thisFlavor = SPECIES[flavNum];
//     cout << "FLAVOR = " << thisFlavor << endl;
//     //Get the center for every flavor
//     int NFromFlavor=0; 

  //Loop all particles
  for(unsigned int i=0;i< particles.size(); i++ )
  {
    //TEST progress:
    //cout << double(i)/particles.size()*100.00 << "%" << endl;
    //Parton
    px = particles[i].Mom.Px(); 
    py = particles[i].Mom.Py();
    pz = particles[i].Mom.Pz();
    phi = getAngle2pi(px,py);
    
    x = particles[i].Pos.X()-offsetX,
    y = particles[i].Pos.Y()-offsetY;
      
    r = sqrt(pow(x,2.0)+pow(y,2.0));
    phiSPACE =  getAngle2pi(x,y);      
    
    double eta = particles[i].Mom.Rapidity();
    
    //PartonPT
    ptP=particles[i].Mom.Pt();

    
    int yRangeIndex = 0;//Midrapidity
    if ( fabs( eta ) >= rapidityRanges[yRangeIndex].yleft && fabs( eta ) <= rapidityRanges[yRangeIndex].yright )
    {
      if(ptP>0.4)
      {
        TotalMultInAcceptanceGluon += 1.0;
      }
      TotalMultGluon += 1.0;
      //Do the Integral over z
      double z_loop=zmin;
      while(z_loop<zmax)
      {    
        //cout << z_loop << "\t" << KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR) << endl;
        ptH=z_loop*ptP;
        Q=ptP;
        
        //Total Multiplicity
        if(ptH>0.4)
        {
          TotalMultInAcceptance += 1.0 * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR);
        }
        TotalMult += 1.0 * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR); 
        
        //REFERENCE  (all results times dz because simple integration here)
        if( (pTHadronReferenceMin < ptH) && (ptH < pTHadronReferenceMax) )
        {
          //B0ref
          rB0ref +=  1.0 * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR);
          iB0ref +=  0.0 * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR);
          //B1ref
          rB1ref += cos(phi) * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR);
          iB1ref += sin(phi) * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR);      
          //B2ref
          rB2ref += cos(2.0*phi) * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR);
          iB2ref += sin(2.0*phi) * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR);  
          //B3ref
          rB3ref += cos(3.0*phi) * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR);
          iB3ref += sin(3.0*phi) * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR); 
          //B4ref
          rB4ref += cos(4.0*phi) * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR);
          iB4ref += sin(4.0*phi) * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR); 
        }
        
        //pT dependent, now for HADRONS!
        pt = ptH;
        
        rB0.add(pt,1.0 * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR));
        iB0.add(pt,0.0 * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR));
        
        rB1.add(pt,cos(phi) * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR));
        iB1.add(pt,sin(phi) * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR));
        
        rB2.add(pt,cos(2.0*phi) * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR));
        iB2.add(pt,sin(2.0*phi) * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR));
        
        rB3.add(pt,cos(3.0*phi) * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR));
        iB3.add(pt,sin(3.0*phi) * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR));

        rB4.add(pt,cos(4.0*phi) * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR));
        iB4.add(pt,sin(4.0*phi) * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR));
        
        //Eccentricities --> doesn't change with Fragmentation
        avcosN1 += pow(r,2.0)*cos(1.0*phiSPACE);
        avsinN1 += pow(r,2.0)*sin(1.0*phiSPACE); 
        
        avrN1 += r;
        avcosN1r1 += pow(r,1.0)*cos(1.0*phiSPACE);
        avsinN1r1 += pow(r,1.0)*sin(1.0*phiSPACE);  
        
        avrN2 += pow(r,2.0);
        avcosN2 += pow(r,2.0)*cos(2.0*phiSPACE);
        avsinN2 += pow(r,2.0)*sin(2.0*phiSPACE);
        
        avrN3 += pow(r,3.0);
        avcosN3 += pow(r,3.0)*cos(3.0*phiSPACE); 
        avsinN3 += pow(r,3.0)*sin(3.0*phiSPACE);
        
        avrN4 += pow(r,4.0);
        avcosN4 += pow(r,4.0)*cos(4.0*phiSPACE); 
        avsinN4 += pow(r,4.0)*sin(4.0*phiSPACE);  
        
        //Energyweighted --> changes with Fragmentation
        avcosN1_energyweight += pt * pow(r,2.0)*cos(1.0*phiSPACE) * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR);
        avsinN1_energyweight += pt * pow(r,2.0)*sin(1.0*phiSPACE) * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR); 
        
        avrN1_energyweight += pt * r * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR);
        avcosN1r1_energyweight += pt * pow(r,1.0)*cos(1.0*phiSPACE) * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR);
        avsinN1r1_energyweight += pt * pow(r,1.0)*sin(1.0*phiSPACE) * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR);  
        
        avrN2_energyweight += pt * pow(r,2.0);
        avcosN2_energyweight += pt * pow(r,2.0)*cos(2.0*phiSPACE) * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR);
        avsinN2_energyweight += pt * pow(r,2.0)*sin(2.0*phiSPACE) * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR);          
        
        avcosN3_energyweight += pt * pow(r,3.0)*cos(3.0*phiSPACE) * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR); 
        avsinN3_energyweight += pt * pow(r,3.0)*sin(3.0*phiSPACE) * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR);          
        
        avcosN4_energyweight += pt * pow(r,4.0)*cos(4.0*phiSPACE) * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR); 
        avsinN4_energyweight += pt * pow(r,4.0)*sin(4.0*phiSPACE) * dz * KKPFragmentation(7,1,z_loop,Q,particles[i].FLAVOR); 
        
        //Increase dz
        z_loop+=dz;        
      }  
    }        
  }
  
  //Normalize
  avrN1 /= rB0ref;
  avrN2 /= rB0ref;
  avrN3 /= rB0ref;
  avcosN1 /= rB0ref; 
  avsinN1 /= rB0ref; 
  avcosN2 /= rB0ref; 
  avsinN2 /= rB0ref; 
  avcosN3 /= rB0ref; 
  avsinN3 /= rB0ref;
  avcosN4 /= rB0ref; 
  avsinN4 /= rB0ref; 
  
  avrN1_energyweight /= rB0ref;
  avrN2_energyweight /= rB0ref;
  avrN3_energyweight /= rB0ref;
  avcosN1_energyweight /= rB0ref; 
  avsinN1_energyweight /= rB0ref; 
  avcosN2_energyweight /= rB0ref; 
  avsinN2_energyweight /= rB0ref; 
  avcosN3_energyweight /= rB0ref; 
  avsinN3_energyweight /= rB0ref;    
  avcosN4_energyweight /= rB0ref; 
  avsinN4_energyweight /= rB0ref;    
  
  eccentricity2 = sqrt(pow(avcosN2,2.0)+pow(avsinN2,2.0))/avrN2;
  eccentricity3 = sqrt(pow(avcosN3,2.0)+pow(avsinN3,2.0))/avrN3; 
  eccentricity4 = sqrt(pow(avcosN4,2.0)+pow(avsinN4,2.0))/avrN4; 

  angleEccentricity = 1./2. * (atan(iB2ref/rB2ref) + M_PI);
  
  //Print out the correlations for vnAnalysis
  {
    fstream file(filename.c_str(),ios::out);
    file  << "#REFERENCE " 
          << rB0ref/Testparticles
          << ' ' 
          << iB0ref/Testparticles 
          << ' ' 
          << rB1ref/Testparticles 
          << ' ' 
          << iB1ref/Testparticles
          << ' ' 
          << rB2ref/Testparticles 
          << ' ' 
          << iB2ref/Testparticles
          << ' ' 
          << rB3ref/Testparticles 
          << ' ' 
          << iB3ref/Testparticles
          << ' ' 
          << rB4ref/Testparticles 
          << ' ' 
          << iB4ref/Testparticles
          << ' '
          << endl;
         
          
    for (int i=0; i<rB0.getNBins(); i++)
    {
      file  << rB0.getBinLabel(i) 
            << "\t" << rB0.getBinValue(i)/Testparticles 
            << "\t" << iB0.getBinValue(i)/Testparticles  
            << "\t" << rB1.getBinValue(i)/Testparticles 
            << "\t" << iB1.getBinValue(i)/Testparticles 
            << "\t" << rB2.getBinValue(i)/Testparticles 
            << "\t" << iB2.getBinValue(i)/Testparticles
            << "\t" << rB3.getBinValue(i)/Testparticles 
            << "\t" << iB3.getBinValue(i)/Testparticles
            << "\t" << rB4.getBinValue(i)/Testparticles 
            << "\t" << iB4.getBinValue(i)/Testparticles
            << endl;

    }    
    file << "#ECCENTRICITY"     //1
          << ' '
          << eccentricity2       //2
          << ' '
          << eccentricity3       //3
          << ' '
          << eccentricity4       //4 new
          << ' '
          << angleEccentricity   //5
          << ' '
          << avcosN1r1           //6
          << ' '
          << avsinN1r1           //7
          << ' '
          << avcosN1             //8
          << ' '
          << avsinN1             //9
          << ' '
          << avrN1               //10
          << ' '
          << avcosN2             //11
          << ' '
          << avsinN2             //12
          << ' '       
          << avrN2               //13
          << ' '
          << avcosN3             //14
          << ' '
          << avsinN3             //15
          << ' '
          << avrN3               //16
          << ' '
          << avcosN4             //17
          << ' '
          << avsinN4             //18
          << ' '
          << avrN4               //19
          << ' '
          << avcosN1r1_energyweight  //20
          << ' '
          << avsinN1r1_energyweight  //21
          << ' '
          << avcosN1_energyweight    //22
          << ' '
          << avsinN1_energyweight    //23
          << ' '
          << avrN1_energyweight      //24
          << ' '
          << avcosN2_energyweight    //25
          << ' '
          << avsinN2_energyweight    //26
          << ' '       
          << avrN2_energyweight      //27
          << ' '
          << avcosN3_energyweight    //28
          << ' '
          << avsinN3_energyweight    //29
          << ' '
          << avrN3_energyweight      //30
          << ' '
          << avcosN4_energyweight    //31
          << ' '
          << avsinN4_energyweight    //32
          << ' '
          << avrN4_energyweight      //33
          << ' '
          << rB0ref                  //34
          << ' '
          << TotalMultInAcceptance/Testparticles    //35
          << ' ' 
          << TotalMult/Testparticles                //36
          << endl;
      
    file.close();   
  }
  cout << "dN^Hadron/dy(pTHadron>0.4 GeV) = " << TotalMultInAcceptance/Testparticles << endl;
  cout << "dN^Parton/dy(pTParton>0.4 GeV) = " << TotalMultInAcceptanceGluon/Testparticles << endl;
}

void analysis::PlotCellwise(double time,int number, int IX, int IY, int centralEtaIndex,const vector<cellContainer> &cells )
{
  if(time < 0.0001)
  {
    time=0.0;
  }
  double Testparticles = theConfig->getTestparticles();  
  stringstream ss;
  ss << number;
  //creates filename, for example: "./output/run32_corr_time0.5.dat"
  string filename           = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + "_cell_parameters_time" + ss.str() + ".dat";
  string filenameParticleNo = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + "_particleNo_time" + ss.str() + ".dat"; 
/*  
  double transLen = theConfig->getTransLength();  
  double transBinSize = 0.08;//fm
  int nTransBins = transLen/transBinSize;
  binningValues2d theBins;
  binningValues2d theBinsParticleNumbers;
  theBins.setMinMaxN(-transLen/2., transLen/2., nTransBins, -transLen/2.,transLen/2., nTransBins);
  theBinsParticleNumbers.setMinMaxN(-transLen/2., transLen/2., nTransBins, -transLen/2.,transLen/2., nTransBins);
  double minRap = theConfig->getMinimumEtaInitial();
  double maxRap = theConfig->getMaximumEtaInitial();
  double rapScaleFactor = std::abs(maxRap-minRap);
  
  
  int id;
  int IXY=IX*IY;
  int n_z = cellindex / IXY; 
  int n_y = (cellindex-n_z*IXY) / IX;
  int n_x = (cellindex-n_z*IXY) - (n_y*IX);  
  double T_AMY=0.;
  lorentz LorentzBoost;
  int centralEtaIndex = static_cast<int>( etaBins.size() ) / 2;
  string full_filename = theConfig->getStandardOutputDirectoryName() + "/" +theConfig->getJobName() + "_T_AMY_Central" +  ".dat";
  T_AMY=0;

  fstream outfile( full_filename.c_str(), ios::out | ios::app);
  outfile.precision( 8 );
  
  //Central: IX/2-1*/
 
  int IXY = IX * IY;
  
  fstream file(filename.c_str(),ios::out);
  
  for ( int cellindex = IXY * centralEtaIndex; cellindex < IXY * ( centralEtaIndex + 1 ); cellindex++ )
  { 
    int id;
    int IXY=IX*IY;
    int n_z = cellindex / IXY; 
    int n_y = (cellindex-n_z*IXY) / IX;
    int n_x = (cellindex-n_z*IXY) - (n_y*IX);  
    
    VectorXYZ  velocity;
    
    cells[cellindex].getBoostvector(velocity);
    
    double dx = velocity.X();
    double dy = velocity.Y();
    double lengthTransverse = sqrt(pow(dx,2.0)+pow(dy,2.0));
    dx/=lengthTransverse*10;
    dy/=lengthTransverse*10;
    double gammaVal = cells[cellindex].getGamma();
    double vol = cells[cellindex].getVolumeFromAVG();
    
    file << (cells[cellindex].corner.x_min+cells[cellindex].corner.x_max)/2+cells[cellindex].corner.x_min << "\t" << (cells[cellindex].corner.y_min+cells[cellindex].corner.y_max)/2+cells[cellindex].corner.y_min << "\t" << cells[cellindex].getParticleDensity() << "\t" << cells[cellindex].getEnergyDensity()<< "\t" << dx << "\t" << dy << "\t" << gammaVal << "\t" << vol << endl;
    
  } 

  file.close();   


}

void analysis::Energydensityplot(double time)
{
  if(time < 0.0001)
  {
    time=0.0;
  }
  double Testparticles = theConfig->getTestparticles();  
  stringstream ss;
  ss << time;
  //creates filename, for example: "./output/run32_corr_time0.5.dat"
  string filename           = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + "_energydensity_time" + ss.str() + ".dat";
  string filenameParticleNo = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + "_particleNo_time" + ss.str() + ".dat"; 
  double transLen = theConfig->getTransLength();  
  double transBinSize = 0.08;//fm
  int nTransBins = transLen/transBinSize;
  binningValues2d theBins;
  binningValues2d theBinsParticleNumbers;
  theBins.setMinMaxN(-transLen/2., transLen/2., nTransBins, -transLen/2.,transLen/2., nTransBins);
  theBinsParticleNumbers.setMinMaxN(-transLen/2., transLen/2., nTransBins, -transLen/2.,transLen/2., nTransBins);
  double minRap = theConfig->getMinimumEtaInitial();
  double maxRap = theConfig->getMaximumEtaInitial();
  double rapScaleFactor = std::abs(maxRap-minRap);
  //cout << "rapidity scale factor " << rapScaleFactor << endl;
  for ( unsigned int i = 0; i < particles.size(); i++ )
  {
    if(particles[i].Pos.T()>time)
    { 
      continue;
    }
    double zpos = particles[i].Pos.Z();
    if((pow(time,2.0)>pow(zpos,2.0)))
    {
      double tau = particles[i].Pos.Eigentime();
      double xpos = particles[i].Pos.X();
      double ypos = particles[i].Pos.Y();
//       if(xpos >=  1.28 || xpos <=  -1.28)
//       {
//         cout << "ana:          Particle: " << i << "  Particle momentum: " << particles[i].Mom << '\t' << "Particle position: " << particles[i].Pos << endl; 
//       }
      double edens = particles[i].Mom.Pt()/rapScaleFactor/Testparticles/tau;///tau/rapScaleFactor/Testparticles;
    //cout << particles[i].Mom.Pt()/tau/rapScaleFactor/Testparticles << endl;//particles[i].Mom.Pt()/tau/rapScaleFactor/Testparticles << endl;
      double particleNo = 1.0/Testparticles;
      theBins.add(xpos,ypos,edens);
      theBinsParticleNumbers.add(xpos,ypos,particleNo);
    }
    //.add(particles[i].Pos.X(),particles[i].Pos.Y(),particles[i].Mom.Pt()/tau/rapScaleFactor/Testparticles);
    //cout << "add" << endl;
  }  
  
//Either automatic print:  
//   theBins.setFilename(filename); 
//   theBins.print();
//   theBinsParticleNumbers.setFilename(filenameParticleNo); 
//   theBinsParticleNumbers.print();  
   
//...or print energydensity and particle-number, both
  fstream file(filename.c_str(),ios::out);
  file << "#Testparticles: " << Testparticles << "   bin size: " << transBinSize <<  "   NBins: " << nTransBins << endl;
  for (int i=0; i<nTransBins; i++)
  {
    for (int j=0; j<nTransBins; j++)
    {      
      file << theBins.getBinLabelX(i);
      file << "\t";
      file << theBins.getBinLabelX(j);
      file << "\t";
      file << theBins.getBin (  i, j );
      file << "\t";
      file << theBinsParticleNumbers.getBin(i,j);
      file << endl;
    }
  }   
  file.close();  

}


void analysis::addJetEvent_initial( const int jetID )
{
  jetTrackerSingleEvent tempEvent;
  tempEvent.jet_ID_in = -1;
  tempEvent.jet_ID_out = jetID;
  tempEvent.coll_type = initial_jet;

  tempEvent.R_proj = particles[jetID].Pos;
  tempEvent.P_proj_out = particles[jetID].Mom;

  vector<jetTrackerSingleEvent> tempVec;
  tempVec.push_back( tempEvent );
  jetTracker.push_back( tempVec );
}


void analysis::addJetEvents_final()
{
  double pt;
  for ( unsigned int i = 0; i < particles.size(); i++ )
  {
    pt = particles[i].Mom.Pt();
    if ( pt > jetTracking_PT )
    {
      jetTrackerSingleEvent tempEvent;
      tempEvent.jet_ID_in = i;
      tempEvent.jet_ID_out = -1;

      tempEvent.R_proj = particles[i].Pos;
      tempEvent.P_proj_in = particles[i].Mom;

      tempEvent.coll_type = final_jet;


      unsigned int entity_index = 0;
      while ( entity_index < jetTracker.size() && jetTracker[entity_index].back().jet_ID_out != i )
      {
        entity_index++;
      }

      if ( entity_index < jetTracker.size() )
      {
        jetTracker[entity_index].push_back( tempEvent );
      }
    }
  }

}


int analysis::addJetEvent_in( const int ID_1, const int ID_2, const int ID_3, const jetTrackerCollType coll_type,
                              const double cross_section, const int cell_ID, const double lambda )
{
  double pt1 = particles[ID_1].Mom.Pt();
  double pt2 = particles[ID_2].Mom.Pt();

  double pt3 = -1;
  if ( ID_3 > 0 )
  {
    pt3 = particles[ID_3].Mom.Pt();
  }

  int jetID;
  int partner1 = -1, partner2 = -1;

  if ( pt1 > pt2 )
  {
    jetID = ID_1;
    partner1 = ID_2;
    partner2 = ID_3;
  }
  else
  {
    jetID = ID_2;
    partner1 = ID_1;
    partner2 = ID_3;
  }
  if ( coll_type == c3to2 && pt3 > pt1 && pt3 > pt2 )
  {
    jetID = ID_3;
    partner1 = ID_1;
    partner2 = ID_2;
  }

  jetTrackerSingleEvent tempEvent;
  tempEvent.jet_ID_in = jetID;

  tempEvent.R_proj = particles[jetID].Pos;
  tempEvent.P_proj_in = particles[jetID].Mom;
  tempEvent.P1_in = particles[partner1].Mom;

  if ( coll_type == c3to2 )
  {
    tempEvent.P2_in = particles[partner2].Mom;
  }

  tempEvent.lambda = lambda;
  tempEvent.xSection = cross_section;
  tempEvent.cell_ID = cell_ID;
  tempEvent.coll_type = coll_type;


  unsigned int entity_index = 0;
  while ( entity_index < jetTracker.size() && jetTracker[entity_index].back().jet_ID_out != jetID )
  {
    entity_index++;
  }

  if ( entity_index < jetTracker.size() )
  {
    jetTracker[entity_index].push_back( tempEvent );
  }
  else
  {
    vector< jetTrackerSingleEvent > tempEntity;
    tempEntity.push_back( tempEvent );
    jetTracker.push_back( tempEntity );
  }

  return entity_index;
}


void analysis::addJetEvent_out( const int entity_ID, const int ID_1, const int ID_2, const int ID_3 )
{
  double pt1 = particles[ID_1].Mom.Pt();
  double pt2 = particles[ID_2].Mom.Pt();

  double pt3 = -1;
  if ( ID_3 > 0 )
  {
    pt3 = particles[ID_3].Mom.Pt();
  }

  int jetID;
  int partner1 = -1, partner2 = -1;

  if ( pt1 > pt2 )
  {
    jetID = ID_1;
    partner1 = ID_2;
    partner2 = ID_3;
  }
  else
  {
    jetID = ID_2;
    partner1 = ID_1;
    partner2 = ID_3;
  }
  if ( ID_3 > 0 && pt3 > pt1 && pt3 > pt2 )
  {
    jetID = ID_3;
    partner1 = ID_1;
    partner2 = ID_2;
  }

  if ( entity_ID != -1 )
  {
    jetTracker[entity_ID].back().jet_ID_out = jetID;

    jetTracker[entity_ID].back().P_proj_out = particles[jetID].Mom;
    jetTracker[entity_ID].back().P1_out = particles[partner1].Mom;

    if ( ID_3 > 0 )
    {
      jetTracker[entity_ID].back().P2_out = particles[partner2].Mom;
    }
  }
  else if ( particles[jetID].Mom.Pt() > jetTracking_PT )
  {
    jetTrackerSingleEvent tempEvent;
    tempEvent.jet_ID_in = -1;
    tempEvent.jet_ID_out = jetID;
    tempEvent.coll_type = production;

    tempEvent.R_proj = particles[jetID].Pos;

    tempEvent.P_proj_out = particles[jetID].Mom;
    tempEvent.P1_out = particles[partner1].Mom;

    if ( ID_3 > 0 )
    {
      tempEvent.P2_out = particles[partner2].Mom;
    }

    vector< jetTrackerSingleEvent > tempEntity;
    tempEntity.push_back( tempEvent );
    jetTracker.push_back( tempEntity );
  }

  double pt_partner1 = particles[partner1].Mom.Pt();
  double pt_partner2 = -1;
  if ( ID_3 > 0 )
  {
    pt_partner2 = particles[partner2].Mom.Pt();
  }

  if ( pt_partner1 > jetTracking_PT )
  {
    jetTrackerSingleEvent tempEvent;
    tempEvent.jet_ID_in = jetID;
    tempEvent.jet_ID_out = partner1;
    tempEvent.coll_type = production;

    tempEvent.R_proj = particles[partner1].Pos;
    tempEvent.P_proj_out = particles[partner1].Mom;
    tempEvent.P1_out = particles[jetID].Mom;

    if ( ID_3 > 0 )
    {
      tempEvent.P2_out = particles[partner2].Mom;
    }

    vector< jetTrackerSingleEvent > tempEntity;
    tempEntity.push_back( tempEvent );
    jetTracker.push_back( tempEntity );
  }
  if ( pt_partner2 > jetTracking_PT )
  {
    jetTrackerSingleEvent tempEvent;
    tempEvent.jet_ID_in = jetID;
    tempEvent.jet_ID_out = partner2;
    tempEvent.coll_type = production;
    tempEvent.R_proj = particles[partner2].Pos;
    tempEvent.P_proj_out = particles[partner2].Mom;
    tempEvent.P1_out = particles[jetID].Mom;
    tempEvent.P2_out = particles[partner1].Mom;

    vector< jetTrackerSingleEvent > tempEntity;
    tempEntity.push_back( tempEvent );
    jetTracker.push_back( tempEntity );
  }


}


void analysis::removeJetEvent_in( const int entity_ID )
{
  jetTracker[entity_ID].pop_back();
}


void analysis::makeJetTrackerCopy()
{
  jetTracker_copy = jetTracker;
}


void analysis::restoreJetTracker()
{
  jetTracker = jetTracker_copy;
}


void analysis::exchangeJetID( const int oldID, const int newID )
{
  unsigned int entity_index = 0;
  while ( entity_index < jetTracker.size() && jetTracker[entity_index].back().jet_ID_out != oldID )
  {
    entity_index++;
  }

  if ( entity_index < jetTracker.size() )
  {
    jetTracker[entity_index].back().jet_ID_out = newID;
  }
  else
  {
    cout << "error in exchangeJetID()" << endl;
  }

}


void analysis::collectPtDataInitial()
{
  ptDistribution( gluon, particles, particles.size(), 0 );
  ptDistribution( up, particles, particles.size(), 0 );
  ptDistribution( down, particles, particles.size(), 0 );
  ptDistribution( strange, particles, particles.size(), 0 );
  ptDistribution( anti_up, particles, particles.size(), 0 );
  ptDistribution( anti_down, particles, particles.size(), 0 );
  ptDistribution( anti_strange, particles, particles.size(), 0 );
  ptDistribution( light_quark, particles, particles.size(), 0 );
  ptDistribution( allFlavors, particles, particles.size(), 0 );

  ptSoftDistribution( gluon, particles, particles.size(), 0 );
  ptSoftDistribution( up, particles, particles.size(), 0 );
  ptSoftDistribution( down, particles, particles.size(), 0 );
  ptSoftDistribution( strange, particles, particles.size(), 0 );
  ptSoftDistribution( anti_up, particles, particles.size(), 0 );
  ptSoftDistribution( anti_down, particles, particles.size(), 0 );
  ptSoftDistribution( anti_strange, particles, particles.size(), 0 );
  ptSoftDistribution( light_quark, particles, particles.size(), 0 );
  ptSoftDistribution( allFlavors, particles, particles.size(), 0 );
}



void analysis::collectPtData( const int step )
{
  //are issued starting with 0
  //0 is reserved for initial output, thus step is incremented
  ptDistribution( gluon, particles, particles.size(), step + 1 );
  ptDistribution( up, particles, particles.size(), step + 1 );
  ptDistribution( down, particles, particles.size(), step + 1 );
  ptDistribution( strange, particles, particles.size(), step + 1 );
  ptDistribution( anti_up, particles, particles.size(), step + 1 );
  ptDistribution( anti_down, particles, particles.size(), step + 1 );
  ptDistribution( anti_strange, particles, particles.size(), step + 1 );
  ptDistribution( light_quark, particles, particles.size(), step + 1 );
  ptDistribution( allFlavors, particles, particles.size(), step + 1 );
  
  ptSoftDistribution( gluon, particles, particles.size(), step + 1 );
  ptSoftDistribution( up, particles, particles.size(), step + 1 );
  ptSoftDistribution( down, particles, particles.size(), step + 1 );
  ptSoftDistribution( strange, particles, particles.size(), step + 1 );
  ptSoftDistribution( anti_up, particles, particles.size(), step + 1 );
  ptSoftDistribution( anti_down, particles, particles.size(), step + 1 );
  ptSoftDistribution( anti_strange, particles, particles.size(), step + 1 );
  ptSoftDistribution( light_quark, particles, particles.size(), step + 1 );
  ptSoftDistribution( allFlavors, particles, particles.size(), step + 1 );
}


void analysis::generalOutput( const int step, const double timeshift )
{

  
  if ( ( studyParticleSpectra && ( step == -1  ||  step == (nTimeSteps - 1) ) ) 
    || studyParticleSpectraAllSteps )
  {
    collectPtData( step );
  }
  
  //calls from rhic::rhic_cell(analysis &) are issued starting with 0
  //0 is reserved for initial output, thus step is incremented
  if( studyHeavyQuarkProduction )
  {
    charmNumber( step + 1 );
    bottomNumber( step + 1 );
  }

  if ( ( studyDetailedParticleOutput && ( step == -1  ||  step == (nTimeSteps - 1) ) ) 
    || studyDetailedParticleOutputAllSteps )
  {
    particleOutput( step + 1 );
  }

  if ( studyDetailedParticlePosMomOutput )
  {
    ParticlePosMomOutput( step + 1, timeshift);
  }
  
  if ( studyDndy && ( step == -1  ||  step == (nTimeSteps - 1) ) )
    print_dndy( step + 1 );
  if ( studyDndy &&  (step != -1)  &&  (step != (nTimeSteps - 1) ) )
  {
    print_dndy( step );
  }
  
  if( studyV2 && ( step == -1 ) )
    computeV2RAA( "initial", 0  );
  
  if( studyV2 && ( step == (nTimeSteps - 1) ) )
    computeV2RAA( "final", tstep[nTimeSteps - 1] );


  //for hydroAnalysis
  //----------------------------------------------//
  if( studyHydro )
  {
    int nn = step + 1;
    if(nn < getRealNumberOfTimeStepsForAnaHydro())
      {
        theAnaHydro->hydroDistribution( nn, theAnaHydro->anaHydroProfileNormal );
        theAnaHydro->hydroDistribution( nn, theAnaHydro->anaHydroProfileArrow );
        
        if(nn > 0)
        {
          double time = tstep[step];
          theAnaHydro->hydroDistributionMidRap( nn, time, timeshift, theAnaHydro->anaHydroProfileMidRapNormal );
          theAnaHydro->hydroDistributionMidRap( nn, time, timeshift, theAnaHydro->anaHydroProfileMidRapArrow );  
        }
      }
  }
  //----------------------------------------------//  
  
  string name;
  stringstream ss;

  //HACK WARNING
  if( studyCorrelations )
  {
    if( step == -1 )
    {
      //ParticleCorrelations(0,gluon);
      //ParticleCorrelations(0,quark);
      //cout << "Fragment correlations..." << endl;
      //FragmentedCorrelations(0);
      //Energydensityplot(0);  
      //computeV2RAA( "initial", 0  );
    }  
    else 
    {
      double time = tstep[step];
      ss << time;
      name = "time_" + ss.str();
      cout << "ChargedHadronCorrelations(time)" << endl;
      //ChargedHadronCorrelations(time);
      //ParticleCorrelations(time,gluon);
      //ParticleCorrelations(time,quark);
      //cout << "Fragment correlations..." << endl;
      //FragmentedCorrelations(time);  
      //Energydensityplot(time);    
      //computeV2RAA( name, time  );
     }
  }
  if(studyCellwiseVideo)
  {
    cout << "Study cellwise video.." << endl;
  }
  
  if(noAnalysis)
  {
    cout << "NOANALYSIS" << endl;
    //do nothing;
  }
}



void analysis::generalOutput( const anaType at, const double timeshift )
{
  switch ( at )
  {
  case initial:
    generalOutput( -1, timeshift );            //see +1 in generalOutput(const int)
    break;
  case final:
    generalOutput( nTimeSteps - 1, timeshift );    //see +1 in generalOutput(const int)
    break;
  default:
    cout << "error in call to generalOutput(const anaType)" << endl;
  }
}

void analysis::writeParticleList(const int step)
{
  string name;
  stringstream ss;
  double px,py,pz;
  double time = tstep[step];
  ss << time;
  name = "time_" + ss.str();
  string filename, filename2;
  filename = theConfig->getStandardOutputDirectoryName() + "/" +  theConfig->getJobName() + "_momenta_time_" +name;
  fstream file( filename.c_str(), ios::out | ios::app );
  for ( unsigned int i = 0; i < particles.size(); i++ )
  {
    double eta = particles[i].Mom.Rapidity();
    //for ( int yRangeIndex = 0; yRangeIndex < rapidityRanges.size(); yRangeIndex++ )
    int yRangeIndex = 0;//Midrapidity
    if ( fabs( eta ) >= rapidityRanges[yRangeIndex].yleft && fabs( eta ) <= rapidityRanges[yRangeIndex].yright )
    {
      px = particles[i].Mom.Px(); 
      py = particles[i].Mom.Py();
      pz = particles[i].Mom.Pz();
      file << px << "\t" << py << "\t" << pz << endl;
    }    
  }
  file.close();
}



void cluster::calcCOM()
{
  std::vector<int>::iterator a;
  VectorTXYZ zero(0,0,0,0);
  centerofmass=zero;centerofmassMom=zero;
  centerofmassVel=zero;
  VectorEPxPyPz addVel;
 
  for (  a = particleList.begin(); a != particleList.end(); a++ )
  {
    int myIndex = *a;
    addVel = particles[myIndex].Mom;
    addVel.NormalizeToE();
    centerofmass=centerofmass+particles[myIndex].Pos;
    centerofmassMom = centerofmassMom + particles[myIndex].Mom;
    centerofmassVel = centerofmassVel + addVel;
  }
  centerofmassVel   *= 1./particleList.size();
  centerofmass      *= 1./particleList.size();
  centerofmassMom   *= 1./particleList.size();
}


void cluster::calcVolume()
{
  double D;
  lorentz L;
  S=0;
  volume=0;
  double meanR2=0;
  double minDistance=infinity;
  std::vector<int>::iterator it2,it1;
  for ( it1 = particleList.begin(); it1 != particleList.end(); ++it1 )
  {
    for ( it2 = particleList.begin(); it2 != particleList.end(); ++it2 )
    {
      if (*it1==*it2) continue;        
      
      //D=L.getSpatialDistance(particles[*it1].Pos,particles[*it2].Pos);
      D=L.getDCAsquared(particles[*it1].Pos,particles[*it2].Pos,particles[*it1].Mom,particles[*it2].Mom);

      S+=D;
      meanR2+=D;
      if (D>volume)
      {  
        volume=D;       
      }
      if(minDistance<D)
      {
        minDistance=D;
        MotherA = *it1;
        MotherB = *it2;
      }
      
    }   
  }
  if(particleList.size()>0)
  {
    meanR2 /= (  particleList.size()*(particleList.size()-1)  );
    S/=   (  particleList.size()*(particleList.size()-1)  );
  }
  S*=50.; //works but shrinks cluster
  
  
//   cout << "max dist " << volume << "\t\t mean dist S " << S << endl;
  
  
  //S=S + pow((pow(meanR2,3./2.)-pow(2.01,3.0)),4./3.);
  
//   meanR2=volume;
  
//   cout << "Vol mean " << pow(meanR2,3./2.) << "\t\t Vol wunsch " <<  pow(0.2,3.0)  << " S = " << pow((pow(meanR2,3./2.)-pow(0.2,3.0)),2.0) << endl;
  //works with defined radius
  S=  200*pow(meanR2-0.001,2.0);
//   cout << "max dist " << volume << "\t\t mean dist S " << meanR2 << "\t\t wunsch dist " << 0.05 << "\t\t\tS=" << S << endl;
  
  
  
//   S=10*pow(particleList.size()-200.,2.0)/particleList.size();
//   cout << "S = " << S << endl;
  //cout << "RATIO " << 100*S/(0.001*pow(particleList.size()-100.,2.0)) << endl;
  
  
  //S=100*S; //+0.001*pow(particleList.size()-100.,2.0);
  
  
}

void cluster::calcInvMass()
{
//ToDo
}

void cluster::getDistToCenter()
{
  lorentz L;
  vector<int>::iterator it1;
  distToCenterOfMass.clear();
  for ( it1 = particleList.begin(); it1 != particleList.end(); ++it1 )
  {
    double AbstandSchwerpunkt=L.getSpatialDistance(particles[(*it1)].Pos,centerofmass);
    distToCenterOfMass.push_back(AbstandSchwerpunkt);
  }
}

double analysis::getTotalAction(vector< cluster > allClusters)
{
  double STot=0;
  vector<cluster>::iterator it;
  for(it=allClusters.begin(); it != allClusters.end(); it++)
  {
    STot+=(*it).S;
  }
  return STot;
}

int analysis::GetNextAvailableColor()
{
  for(int i=0;i<Colors.size()-1;i++)
  {
    if(Colors[i+1]-Colors[i]!=1)
    {
    return Colors[i] + 1; 
    }
    
  }
  // ELSE //
  return Colors[Colors.size()-1]+1; 
}

void analysis::SplitRandom(cluster Cluster, cluster &NewClusterOne, cluster &NewClusterTwo, int NCluster)
{
  lorentz L;
  int iRandOne,iRandTwo;
  
  double randomA = ran2();
  double randomB = ran2();
  if(randomA!=randomB)
  {
    iRandOne = std::floor(randomA*Cluster.particleList.size());
    iRandTwo = std::floor(randomB*Cluster.particleList.size());    
  }else
  {
    while(randomA==randomB)
    {
      double randomA = ran2();
      double randomB = ran2();    
      iRandOne = std::floor(randomA*Cluster.particleList.size());
      iRandTwo = std::floor(randomB*Cluster.particleList.size());    
    }
  }  

//   cout << "split random " << iRandOne << "\t" << iRandTwo << endl;
  
  //Select widest particle of iOne: iRandTwo
  double distIOneMax=0;
  for ( int i = 0; i < Cluster.particleList.size(); i++ )
  {
   if(i==iRandOne) continue;
   double distOne = L.getDCAsquared( particles[Cluster.particleList[iRandOne]].Pos, particles[Cluster.particleList[i]].Pos, particles[Cluster.particleList[iRandOne]].Mom, particles[Cluster.particleList[i]].Mom   ) ;
   if(distIOneMax<distOne)
   {
     iRandTwo=i;
     distIOneMax=distOne;
   }
  }
  
  NewClusterOne.particleList.clear();
  NewClusterTwo.particleList.clear();
  NewClusterOne.particleList.push_back(Cluster.particleList[iRandOne]);  
  NewClusterTwo.particleList.push_back(Cluster.particleList[iRandTwo]);
  
  for ( int i = 0; i < Cluster.particleList.size(); i++ )
  {
//     double distOne = L.getSpatialDistance( particles[Cluster.particleList[iRandOne]].Pos, particles[Cluster.particleList[i]].Pos  ) ;
//     double distTwo = L.getSpatialDistance( particles[Cluster.particleList[iRandTwo]].Pos, particles[Cluster.particleList[i]].Pos  ) ;
    
    if(i==iRandOne || i==iRandTwo) continue;
      
    
    double distOne = L.getDCAsquared( particles[Cluster.particleList[iRandOne]].Pos, particles[Cluster.particleList[i]].Pos, particles[Cluster.particleList[iRandOne]].Mom, particles[Cluster.particleList[i]].Mom   ) ;
    double distTwo = L.getDCAsquared( particles[Cluster.particleList[iRandTwo]].Pos, particles[Cluster.particleList[i]].Pos, particles[Cluster.particleList[iRandTwo]].Mom, particles[Cluster.particleList[i]].Mom   ) ;
    
    
    
    if( FPT_COMP_L(distOne,distTwo) )
    {
      NewClusterOne.particleList.push_back(Cluster.particleList[i]);        
    }
    else
    {
      NewClusterTwo.particleList.push_back(Cluster.particleList[i]);       
    }  
  
  }
  // CALCULATE CLUSTER PROPERTIES //

  NewClusterOne.calcCOM();
  NewClusterOne.calcVolume();
  NewClusterOne.calcInvMass();
  NewClusterOne.getDistToCenter();

  NewClusterTwo.calcCOM();
  NewClusterTwo.calcVolume();
  NewClusterTwo.calcInvMass();
  NewClusterTwo.getDistToCenter();
  
//   cout <<  NewClusterOne.particleList.size() << "\t" << NewClusterTwo.particleList.size() << endl;
  
//   cout << ".................. split: Cluster before: " << Cluster.S << "\t\t\t\t Cluster new: " << NewClusterOne.S << "\t" << NewClusterTwo.S << "\t\t" << NewClusterOne.S+NewClusterTwo.S << "\t" << exp(-100*(NewClusterOne.S+NewClusterTwo.S-Cluster.S)) <<   endl;
  
}

void analysis::Split(cluster Cluster, cluster &NewClusterOne, cluster &NewClusterTwo, int NCluster)
{
  std::vector<double>::iterator result;
  result = std::max_element(Cluster.distToCenterOfMass.begin(), Cluster.distToCenterOfMass.end());
  int indexMaxDistance=std::distance(Cluster.distToCenterOfMass.begin(), result);
  
  
//   cout << "\n\n\n\n" << endl;
//   cout << indexMaxDistance << "\t" << Cluster.distToCenterOfMass[indexMaxDistance] <<endl;
  
  lorentz L;
  
  
  for ( int i = 0; i < Cluster.particleList.size(); i++ )
  {
    double distCOM = Cluster.distToCenterOfMass[i];
    double distNewCenter = L.getSpatialDistance( particles[Cluster.particleList[i]].Pos, particles[Cluster.particleList[indexMaxDistance]].Pos  ) ;
    if( distCOM <= distNewCenter )
    {
      NewClusterOne.particleList.push_back(Cluster.particleList[i]);        
    }
    else
    {
      NewClusterTwo.particleList.push_back(Cluster.particleList[i]);       
    }  
  
  }
  // CALCULATE CLUSTER PROPERTIES //

  NewClusterOne.calcCOM();
  NewClusterOne.calcVolume();
  NewClusterOne.calcInvMass();
  NewClusterOne.getDistToCenter();

  NewClusterTwo.calcCOM();
  NewClusterTwo.calcVolume();
  NewClusterTwo.calcInvMass();
  NewClusterTwo.getDistToCenter();
  
  
  // BIGGER CLUSTER ALWAYS CLUSTER ONE OTHERWISE SWAP NAMES //
//   cluster tempC;
//   if(NewClusterOne.volume < NewClusterTwo.volume)
//   {
//     tempC = NewClusterTwo;
//     NewClusterTwo = NewClusterOne;
//     NewClusterOne = tempC;
//   }

  
//   cout << ".................. split: Cluster before: " << Cluster.S << "\t\t\t\t Cluster new: " << NewClusterOne.S << "\t" << NewClusterTwo.S << "\t\t" << NewClusterOne.S+NewClusterTwo.S <<  endl;
  
}

void analysis::ClusterSplitStep(double beta, vector<cluster> &theClusters, bool isWarmUp)
{
  
  
  // CLUSTERS TO MBE CONSIDERED FOR SPLITTING //
  int iOne;

  // RANDOMLY PICK A CLUSTER //
  double A=ran2();
//   cout << A << endl;
  iOne=std::floor(A*theClusters.size());
  
  
  //OR TAKE THE BIGGEST CLUSTER //
  double volMax=0;
  int MaxVolumeIndex=0;
  for(vector<cluster>::iterator it=theClusters.begin();it!=theClusters.end();it++)
  {
    if(volMax<(*it).volume)
    {
      volMax=(*it).volume;
      MaxVolumeIndex=std::distance( theClusters.begin(), it );
    }
  }
  //cout << "Vol Max = " << volMax << "   ind " << MaxVolumeIndex << endl;
  //iOne=MaxVolumeIndex;
  
  
  
  
  if (theClusters[iOne].particleList.size()>2)
  {
    cluster NewClusterOne;
    cluster NewClusterTwo;
    
    double SOld=theClusters[iOne].S;
    
    SplitRandom(theClusters[iOne],NewClusterOne,NewClusterTwo,theClusters.size());
    
    //cout << "after split: " << NewClusterOne.S << "\t" << NewClusterTwo.S << endl;
    
    
      // ACCEPT / REJECT //
    double SNew= (NewClusterOne.S+NewClusterTwo.S ) /2.;
    
    if(!isWarmUp) NumberOfSplitAttempts++;
    //sumOfDeltaS+=SNew-SOld;
    
//     cout << "DeltaS = " <<SNew-SOld << "\t\t\t\t\t" << std::exp(-beta*(SNew-SOld)) << endl;
    
    
    
    // ACCEPT CASE //
    double b=ran2();
//     cout << "B " << b << endl;
    if(FPT_COMP_L(b,std::exp(-beta*(SNew-SOld)) ) || FPT_COMP_L(SNew-SOld,0.0) )
    {
//       cout << "Split No. " << iOne << endl;
      // COLORS //
      NewClusterOne.myColor=theClusters[iOne].myColor;
      NewClusterTwo.myColor=GetNextAvailableColor();
//       cout << "GetNextAvailableColor "<< NewClusterTwo.myColor << endl;
      
      Colors.push_back(NewClusterTwo.myColor);
      std::sort(Colors.begin(), Colors.end());
      
      // FIX NEW CLUSTERS //
      theClusters.erase(theClusters.begin() + iOne);
      theClusters.push_back(NewClusterOne);
      theClusters.push_back(NewClusterTwo);

      
      // ACCEPTANCE COUNTER //
      if(!isWarmUp) NumberOfAcceptedSplits++;

    }else
    {
//       cout << "no split" << endl;
    }
    
      
  } 

  
}

void analysis::ClusterMergeStep(double beta, vector<cluster> &theClusters, bool isWarmUp)
{
  lorentz L;
  
  // CLUSTERS TO MBE CONSIDERED FOR MERGING //
  int iOne; int iTwo;

  // RANDOMLY PICK A CLUSTER //
  double C=ran2();
//   cout << C << endl;
  iOne=std::floor(C*theClusters.size());

  // CHECK FOR CLOSEST CLUSTER //
  double MinDistance=std::numeric_limits<double>::max();

  if(theClusters.size()>1)
  {
    for(int j=0;j<theClusters.size();j++)
    {
      if(j!=iOne)
      { 
//         double CurrentDistance=L.getSpatialDistance(theClusters[iOne].centerofmass, theClusters[j].centerofmass); 
        double CurrentDistance=L.getDCAsquared(theClusters[iOne].centerofmass, theClusters[j].centerofmass,theClusters[iOne].centerofmassVel, theClusters[j].centerofmassVel); 
        
        if(CurrentDistance<MinDistance)
        {
          iTwo=j; 
          MinDistance=CurrentDistance;
        }
      }
    } 

    // CREATAE NEW CLUSTER //
    cluster NewCluster=cluster_merge(theClusters[iOne],theClusters[iTwo]);

    
    
    // ACCEPT / REJECT //
    double SNew=NewCluster.S;
    double SOld= ( theClusters[iOne].S+theClusters[iTwo].S ) /2.;

//     cout << ".................. merge: Cluster before: " << SOld << "\t\t\t\t Cluster new: " << SNew <<  endl;
    
    if(!isWarmUp) NumberOfMergeAttempts++;
//     if(FPT_COMP_G(SNew-SOld,0.0) )
//     {  
//       cout << "Merge Probab: " << std::exp(-beta*(SNew-SOld)) << endl;
//     }
    
//     cout << "DeltaS = " <<SNew-SOld << "\t\t\t\t\t" << std::exp(-beta*(SNew-SOld)) << endl;
    
    // ACCEPT CASE //
    double D=ran2();
    
    if(D<std::exp(-beta*(SNew-SOld)) || FPT_COMP_L(SNew-SOld,0.0) )
    {

      //  REPLACE OLD CLUSTERS BY NEW ONES //
      int iS=std::min(iOne,iTwo);
      int iL=std::max(iOne,iTwo);
//       cout << "Merge!\t" << iOne << " with " << iTwo << endl;
      
      int LowestColor=std::min(theClusters[iOne].myColor,theClusters[iTwo].myColor);
      int HighestColor=std::max(theClusters[iOne].myColor,theClusters[iTwo].myColor);  
      
      theClusters.erase(theClusters.begin() + iL);
      theClusters.erase(theClusters.begin() + iS);
//       cout << "erase...\t" << theClusters.size() << endl;
      NewCluster.myColor=LowestColor;
      theClusters.push_back(NewCluster);
//       cout << "push back...\t" << theClusters.size() << endl;
      // COLORING //
      
      
      // COLOR ASSIGNMENT //
      std::vector<int>::iterator it;
      it = find (Colors.begin(), Colors.end(), HighestColor);
      if (it != Colors.end())
      {
        int indexToErase = std::distance(Colors.begin(),it);
//         std::cout << "Color index to erase = " << indexToErase  << '\n';
//         cout << "Color to erase " << Colors[indexToErase] << endl;
        Colors.erase(Colors.begin()+indexToErase);
      }
      else
      {
//         std::cout << "Element not found in myvector\n";
      }
      
//       for(int j=0;j<Colors.size();j++)
//       {
//         if(Colors[j]==theClusters[iTwo].myColor)
//         {
//           Colors.erase(Colors.begin() + j);
//         }
//       }
      std::sort(Colors.begin(), Colors.end());
      // ACCEPTANCE COUNTER //
      if(!isWarmUp) NumberOfAcceptedMerges++;

    }
    
    //  REJECT CASE //
    else
    {
//       cout << "No Merge!" << endl;
    }
  }
}

cluster analysis::cluster_merge(cluster ClusterOne, cluster ClusterTwo)
{

  // CREATE NEW CLUSTER //
  cluster NewCluster;

  // SET PARTICLE LIST //
  (NewCluster.particleList).insert(NewCluster.particleList.end(), ClusterOne.particleList.begin(), ClusterOne.particleList.end());
  (NewCluster.particleList).insert(NewCluster.particleList.end(), ClusterTwo.particleList.begin(), ClusterTwo.particleList.end());

  // CALCULATE CLUSTER PROPERTIES //
  NewCluster.calcCOM();
  NewCluster.calcVolume();
  NewCluster.calcInvMass();
  NewCluster.getDistToCenter();

  return NewCluster;
}

void analysis::printClusters(vector<cluster> allClusters, int no)
{
   string filename,name;
   stringstream ss;
   ss << no;
//   //+ ss.str() 
   filename = filename_prefix + "_" + ss.str() + ".dat";
   fstream file_central( filename.c_str(), ios::out | ios::trunc  );
    
  cout << "Number Clusters: " << allClusters.size() << endl; 
  vector<cluster>::iterator it;
  for(it=allClusters.begin(); it != allClusters.end(); it++)
  {  
//     if ( (*it).myColor > 12) continue;
//     cout << "CLUSTER " << std::distance(allClusters.begin(),it) << "\t\t" << (*it).myColor << "\t\tS= " << (*it).S << "\t\t" << (*it).particleList.size() << endl;
    for(vector<int>::iterator it2=(*it).particleList.begin(); it2 != (*it).particleList.end(); it2++)
    {
      //cout << particles[(*it2)].Pos.X() << "\t" << particles[(*it2)].Pos.Y() << "\t" << (*it).myColor << endl;
      file_central << particles[(*it2)].Pos.X() << "\t" << particles[(*it2)].Pos.Y() << "\t" << particles[(*it2)].Pos.Z() << "\t" << (*it).myColor << endl;
//      cout << particles[(*it2)].Pos << endl;
//       cout << particles[(*it2)].Mom.E() << endl;
    }
    
    
  }
   
  file_central.close(); 
    
}

void analysis::writeAllClustersInLesHouches(vector< cluster > allClusters, int no)
{
  vector<cluster>::iterator it;
  for(it=allClusters.begin(); it != allClusters.end(); it++)
  {  
    cout << "Write LesHouches for Cluster " << std::distance(allClusters.begin(),it) << endl;
    //WARNING
    //generateLesHouchesEvent((*it), std::distance(allClusters.begin(),it)); 
    GenerateDCATimeOrderedColorsLesHouchesEvent((*it), std::distance(allClusters.begin(),it));
    //HACK
//     break;
  }
}


void analysis::Cluster()
{
  
  ran2.setSeed( 17 );
    
  
  lorentz L;
  
  double beta = 1;  // unit ???
  
  theClusters.clear();
  int NCluster=1;
  theClusters.resize( NCluster );

  Colors.clear();
  
  //INSERT ALL PARTICLES IN CLUSTER ZERO//
  for ( int i = 0; i < particles.size(); i++ )
  {
    theClusters[0].particleList.push_back(i);
  }

  //Insert randomly some particles //
//   int numberToInsert = 7000;
//   vector<int> alreadyIn;
//   for(int k=0; k< numberToInsert; k++)
//   {
//     int item = std::floor(ran2()*particles.size());
//     std::vector<int>::iterator it = std::find(theClusters[0].particleList.begin(),theClusters[0].particleList.end(),item);
//     if(it!=theClusters[0].particleList.end())
//     {  
//       //Element found, go back for choosing a new one
//       k--;
//     }
//     else
//     {
//       //Element not yet found
//       theClusters[0].particleList.push_back(item);
//     }
//   }

  // CALCULATE CLUSTER PROPERTIES //
  theClusters[0].calcCOM();
  theClusters[0].calcVolume();
  theClusters[0].calcInvMass();
  theClusters[0].getDistToCenter();  
  theClusters[0].myColor=0;
  Colors.push_back(theClusters[0].myColor);
  //Print
  {
    cout << "CLUSTER initial "  << "\t\t" << theClusters[0].myColor << "\t\tS= " << theClusters[0].S << "\t\t" << theClusters[0].particleList.size() << endl;
  }
  

  //INIT COUNTERS
  NumberOfMergeAttempts=0;
  NumberOfAcceptedMerges=0;
  NumberOfSplitAttempts=0;
  NumberOfAcceptedSplits=0;
  sumOfDeltaS=0.;
 
  bool isWarmUp=true;
  
  for(int k=0; k< 100;k++)
  {
//     cout << "##############  " << k << " ############" << endl;
//     cout << "Colors: ";
//     for(vector<int>::iterator it= Colors.begin();it!=Colors.end();it++)
//     {
//       cout << *it;
//     }
//     cout << endl;
    
    if (k>50) isWarmUp=false;
//     printClusters(theClusters,k); 
//        cout << "Start Split Step" << endl;
       ClusterSplitStep(beta, theClusters, isWarmUp); 
//     printClusters(theClusters,k); 
//        cout << "Start Merge Step" << endl;
       //TODO not quite time invariant
       ClusterMergeStep(beta, theClusters, isWarmUp);
    
    
//     vector<cluster>::iterator it;
//     for(it=theClusters.begin(); it != theClusters.end(); it++)
//     {  
//       cout << "CLUSTER " << std::distance(theClusters.begin(),it) << "\t\t" << (*it).myColor << "\t\tS= " << (*it).S << "\t\t" << (*it).particleList.size() << endl;
//     } 
    
//     cout << "Colors: ";
//     for(vector<int>::iterator it= Colors.begin();it!=Colors.end();it++)
//     {
//       cout << *it;
//     }
//     cout << endl;  
// //     printClusters(theClusters,k); 

    
  }    
  cout << "Number Clusters: " << theClusters.size() << endl; 
  NumberCluster = theClusters.size();
  cout << "SplitRate " << NumberOfAcceptedSplits/double(NumberOfSplitAttempts) << endl;
  cout << "MergeRate " << NumberOfAcceptedMerges/double(NumberOfMergeAttempts) << endl;
  
//   for( int i = 0; i< particles.size();i++)
//   {
// //     particles[i].Propagate( timenext );
//      cout << particles[i].Mom << endl;
// //     cout <<   particles[i].Mom << endl;
//   }
/*  
  for(int i=0; i< theClusters[0].particleList.size(); i++)
  {
    cout << particles[theClusters[0].particleList[i]].Mom << endl;
  }*/
  
  
  printClusters(theClusters,0); 
  writeAllClustersInLesHouches(theClusters,0);
  
  //cout << "AVG DeltaS=" <<  sumOfDeltaS/NumberOfSplitAttempts << endl;
  
  
  
 
  
//   string name;
//   stringstream ss;
//   double px,py,pz;
//   double time = tstep[step];
//   ss << time;
//   name = "time_" + ss.str();
//   string filename, filename2;
//   filename = theConfig->getStandardOutputDirectoryName() + "/" +  theConfig->getJobName() + "_momenta_time_" +name;
//   filename2 = theConfig->getStandardOutputDirectoryName() + "/" +  theConfig->getJobName() + "_momenta_time_" +name;
//   fstream file( filename.c_str(), ios::out | ios::app );
  

  
  double S=0.,S_old=infinity;
  
  double AbstandSchwerpunkt=0;
  double largestDistanceToSchwerpunkt=0;
  int indexLargestDistParticle=0;
  int indexBiggestCluster=0; 
  double largestDistance=0;
  double D=0;
  VectorXYZ schwerpunkt;
  
}  
 
void analysis::writeMomenta(VectorEPxPyPz aVector, string & thestring)
{
  stringstream ss1,ss2,ss3,ss4;
  double E,px,py,pz;
  E= aVector.E();
  px=aVector.Px();
  py=aVector.Py();
  pz=aVector.Pz();
  ss1 << E;
  ss2 << px;
  ss3 << py;
  ss4 << pz;
  thestring = thestring +   ss2.str() + "\t" + ss3.str() + "\t" + ss4.str() + "\t" + ss1.str() + "\t";
}
 
 
void analysis::getBetaToCoM(cluster ACluster, VectorEPxPyPz & beta)
{
  VectorEPxPyPz sumAll(0.,0.,0.,0.);
  for(int i=0;i<ACluster.particleList.size();i++)
  {    
    if(i+1==ACluster.particleList.size())
    { 
      break;
    }
    sumAll+=particles[ACluster.particleList[i]].Mom;
    i++;
    sumAll+=particles[ACluster.particleList[i]].Mom;
  }
  beta = sumAll.NormalizeToE();
}
 
bool analysis::goodQCDCombination(int i, int j)
{
  if (particles[i].FLAVOR==gluon && particles[j].FLAVOR==gluon )
    return true;
  else if( (particles[i].FLAVOR==up && particles[j].FLAVOR==anti_up) || (particles[i].FLAVOR==down && particles[j].FLAVOR==anti_down) || (particles[i].FLAVOR==strange && particles[j].FLAVOR==anti_strange) )
    return true;
  else 
    return false;
}
 
partonCombinationType analysis::whichQCDCombination(int i, int j)
{
  if (particles[i].FLAVOR==gluon && particles[j].FLAVOR==gluon )
    return glue_glue;
  else if( (particles[i].FLAVOR==up && particles[j].FLAVOR==anti_up) || (particles[i].FLAVOR==down && particles[j].FLAVOR==anti_down) || (particles[i].FLAVOR==strange && particles[j].FLAVOR==anti_strange) )
    return quark_antiquark;
  else 
    return glue_quark;
}
 
 
 
 
void analysis::GenerateDCATimeOrderedColorsLesHouchesEvent(cluster oneCluster, int no)
{
  double EnergyBeam=0.;
  stringstream ss1,ss2,ss3,ss4,ss5,ss6,ss7;
  int N = oneCluster.particleList.size();
  int color1,color2;
  cout << "Original Cluster Size " << N << endl;
  
  int Nparts=0;
  
  double TotE=0.;
  VectorEPxPyPz TotMom_Q(0.,0.,0.,0.);
  VectorEPxPyPz TotMom_AQ(0.,0.,0.,0.);  
  VectorEPxPyPz TotMom(0.,0.,0.,0.);
  
  string momentaString;
  string endOfLine="0.0\t9.0\n" ;
  
  
  
  vector<int> ColorStream1,ColorStream2,FinalParticleList,singleGluonList;

  
  
  double D;
  std::vector<int>::iterator it2,it1;
  bool badchoice=true;
  bool gluonQuark=false;
  partonCombinationType thisCombination;
  
  int countColors=0;
  lorentz L;
  do
  {   
    double minDist=infinity;
    //find nearest combination
    for ( int i = 0; i < oneCluster.particleList.size(); i++ )
    {
      for ( int j = i; j < oneCluster.particleList.size(); j++ )
      {
        if (i==j) continue;        
        
        
        //D=L.getSpatialDistance(particles[*it1].Pos,particles[*it2].Pos);
        D=L.getDCAsquared(particles[oneCluster.particleList[i]].Pos,particles[oneCluster.particleList[j]].Pos,particles[oneCluster.particleList[i]].Mom,particles[oneCluster.particleList[j]].Mom);
        if(minDist>D)
        {
          minDist=D;
          oneCluster.MotherA=i;
          oneCluster.MotherB=j;
         
        }
      }
    } 
    
    thisCombination = whichQCDCombination(oneCluster.particleList[oneCluster.MotherA],oneCluster.particleList[oneCluster.MotherB]);
    
    //cout << oneCluster.particleList[oneCluster.MotherA] << "\t\t\t" << oneCluster.particleList[oneCluster.MotherB] <<  "    switch " << thisCombination << endl;
    
    switch(thisCombination)
    {
      case  glue_glue: 
      {
        //cout << oneCluster.MotherA << "\t\t" << oneCluster.MotherB << "\t\t" << oneCluster.MotherB << endl;
        
        FinalParticleList.push_back(oneCluster.particleList[oneCluster.MotherA]);
        ColorStream1.push_back(501+countColors);
        ColorStream2.push_back(502+countColors);
        
        FinalParticleList.push_back(oneCluster.particleList[oneCluster.MotherB]);
        ColorStream1.push_back(502+countColors);
        ColorStream2.push_back(501+countColors);
        
        
        if(oneCluster.MotherA>oneCluster.MotherB)
        { 
          oneCluster.particleList.erase(oneCluster.particleList.begin()+oneCluster.MotherA);
          oneCluster.particleList.erase(oneCluster.particleList.begin()+oneCluster.MotherB);
        }
        else
        {
          oneCluster.particleList.erase(oneCluster.particleList.begin()+oneCluster.MotherB);
          oneCluster.particleList.erase(oneCluster.particleList.begin()+oneCluster.MotherA);          
        }
        

        
        countColors+=2;
        //cout << "GG " << FinalParticleList.size() << "\t\t\t" <<  oneCluster.particleList.size() << endl;
      }; 
      break;
      case quark_antiquark:
      {
        
        if(  ParticlePrototype::mapToGenericFlavorType(particles[oneCluster.particleList[oneCluster.MotherA]].FLAVOR)==light_quark)
        {
          FinalParticleList.push_back(oneCluster.particleList[oneCluster.MotherA]);
          ColorStream1.push_back(501+countColors);
          ColorStream2.push_back(0);
        
        
          FinalParticleList.push_back(oneCluster.particleList[oneCluster.MotherB]);
          ColorStream1.push_back(0);
          ColorStream2.push_back(501+countColors);    
          
        }else
        {
          FinalParticleList.push_back(oneCluster.particleList[oneCluster.MotherA]);
          ColorStream1.push_back(0);
          ColorStream2.push_back(501+countColors);
        
        
          FinalParticleList.push_back(oneCluster.particleList[oneCluster.MotherB]);
          ColorStream1.push_back(501+countColors);
          ColorStream2.push_back(0);  
        }
        
        
        
        if(oneCluster.MotherA>oneCluster.MotherB)
        { 
          oneCluster.particleList.erase(oneCluster.particleList.begin()+oneCluster.MotherA);
          oneCluster.particleList.erase(oneCluster.particleList.begin()+oneCluster.MotherB);
        }
        else
        {
          oneCluster.particleList.erase(oneCluster.particleList.begin()+oneCluster.MotherB);
          oneCluster.particleList.erase(oneCluster.particleList.begin()+oneCluster.MotherA);          
        }
        

        countColors+=2;
        cout << "QQBAR " << oneCluster.particleList.size() << endl;
      }; 
      break;
      case glue_quark:
      {
        if(particles[oneCluster.particleList[oneCluster.MotherA]].FLAVOR==gluon)
        {
          singleGluonList.push_back(oneCluster.particleList[oneCluster.MotherA]);
          oneCluster.particleList.erase(oneCluster.particleList.begin()+oneCluster.MotherA);
        }
        else
        {
          singleGluonList.push_back(oneCluster.particleList[oneCluster.MotherB]);
          oneCluster.particleList.erase(oneCluster.particleList.begin()+oneCluster.MotherB);       
        }
        //cout << "GQ " << oneCluster.particleList.size() << endl;
      }; 
      break;
    }
  }
  while(oneCluster.particleList.size()>1);    
  
  
  //Sort single gluons to nearest Antenna
  
  for(int i=0;i<singleGluonList.size();i++)
  { 
    double minDist=infinity;
    int closestDCA=0;
    for(int j=0;i<FinalParticleList.size();i++)
    {
      D=L.getDCAsquared(particles[singleGluonList[i]].Pos,particles[FinalParticleList[j]].Pos,particles[singleGluonList[i]].Mom,particles[FinalParticleList[j]].Mom);
      if(minDist>D)
      {
        minDist=D;
        closestDCA=j;       
      }
    }
    //Attach i'th gluon to j'th parton in the FinalParticleList but write it simply at the end of the list
    FinalParticleList.push_back(singleGluonList[i]);
    //check colors of parton j
    color1 = ColorStream1[j];
    color2 = ColorStream2[j];
    if(color1>0 && color2>0)
    {
      //TODO
      //Parton j is a GLUON
//       ColorStream1.push_back(ColorStream1[j]);  //write colors for the gluon
//       ColorStream2.push_back(ColorStream1[j]+1); 
//       
//       ColorStream1[j] = 
//       ColorStream2[j] = 
    }
    else if (color1>0)
    {
      //Parton j is a Quark
      ColorStream1.push_back(ColorStream1[j]);  //write colors for the gluon
      ColorStream2.push_back(ColorStream1[j]+1); 
      ColorStream1[j]=ColorStream1[j]+1;        // change quark color
      
      
    }
    else if (color2>0)
    {
      //Parton j is a AntiQuark
            
      ColorStream1.push_back(ColorStream2[j]);  //write colors for the gluon
      ColorStream2.push_back(ColorStream2[j]+1); 
      ColorStream2[j]=ColorStream2[j]+1;        // change antiquark color
    }
    
  }
  
  
 

  VectorEPxPyPz beta;
  cluster newCluster;
  newCluster.particleList.reserve(FinalParticleList.size());
  newCluster.particleList = FinalParticleList;
  getBetaToCoM(newCluster,  beta);
  lorentz boostToCoM;
  boostToCoM.setBeta(beta);

  
  for(int i=0;i<FinalParticleList.size();i++)
  {        
    particles[FinalParticleList[i]].Mom = boostToCoM.boost(particles[FinalParticleList[i]].Mom); 

    int species = ParticlePrototype::mapToPDG(particles[FinalParticleList[i]].FLAVOR);
    ss7.str("");
    ss7 << species;

    ss5.str("");
    ss6.str("");
    ss5 << ColorStream1[i];
    ss6 << ColorStream2[i];
    
    momentaString = momentaString + ss7.str() + "\t+1\t3\t3\t"+ss5.str()+"\t"+ss6.str()+"\t";
    writeMomenta(particles[FinalParticleList[i]].Mom,momentaString);
    momentaString = momentaString +  "0" +"\t"+endOfLine;
    
    Nparts++;
    
    TotMom  += particles[FinalParticleList[i]].Mom;
    TotE    += particles[FinalParticleList[i]].Mom.E();    
    
  }
  
  cout << "Tot MOM = " << TotMom << endl;

  Nparts += 3; // the incoming and intermediate particle count as well
  
  
  double invariantMass=sqrt(pow(TotMom.E(),2.0)-TotMom.vec2());
  ss3 << invariantMass;

  cout << "Final Cluster Size " << Nparts-3 << endl;
   
  ss2 << Nparts;
  string NpartsStr=ss2.str();
  
   
  string lesHouchesString="<LesHouchesEvents version=\"1.0\">\n<headers>\n</headers>\n<init>\n11 -11 " + ss3.str() + " " + ss3.str() + " -1 -1 -1 -1 1 1\n0 0  3 1001\n</init>\n<event>\n";
  lesHouchesString+=NpartsStr+ "\t1001\t1\t-1\t0.00729927007\t0.3\n";


  VectorEPxPyPz TotEPlusMom(TotE/2.,TotE/2.,0,0);
  VectorEPxPyPz TotEMinusMom(TotE/2.,-TotE/2.,0,0);

  lesHouchesString=lesHouchesString+"11\t-1\t0\t0\t0\t0\t";
  writeMomenta(TotEPlusMom,lesHouchesString);   
  lesHouchesString =  lesHouchesString + "0" +"\t"+endOfLine;
  
  
  lesHouchesString=lesHouchesString+"-11\t-1\t0\t0\t0\t0\t";
  writeMomenta(TotEMinusMom,lesHouchesString);
  lesHouchesString = lesHouchesString +"0" +"\t"+endOfLine;
   
     
  VectorEPxPyPz Photon = TotMom;
  
  lesHouchesString=lesHouchesString+"22\t+2\t1\t2\t0\t0\t";
  writeMomenta(Photon,lesHouchesString);
  //invariant mass = energy beam
  lesHouchesString=lesHouchesString+ ss3.str()  +"\t"+endOfLine;
  
  lesHouchesString=lesHouchesString + momentaString;
    

   lesHouchesString=lesHouchesString+"</event>\n</LesHouchesEvents>";

   
   string filename;

   ss4 << no;

   filename = filename_prefix + "_generatedEvent_" + ss4.str() + ".lhe";
   fstream file_central( filename.c_str(), ios::out | ios::trunc  );
   file_central << lesHouchesString;
   file_central.close();
   
}
 
 
 
void analysis::generateLesHouchesEvent(cluster oneCluster, int no)
{
  double EnergyBeam=0.;
  stringstream ss1,ss2,ss3,ss4,ss5,ss6,ss7;
  int N = oneCluster.particleList.size();
  
  cout << "Original Cluster Size " << N << endl;
  
  int Nparts=0;
  
  double TotE=0.;
  VectorEPxPyPz TotMom_Q(0.,0.,0.,0.);
  VectorEPxPyPz TotMom_AQ(0.,0.,0.,0.);  
  VectorEPxPyPz TotMom(0.,0.,0.,0.);
  
  string momentaString;
  string endOfLine="0.0\t9.0\n" ;
  
  VectorEPxPyPz beta;
  getBetaToCoM(oneCluster,  beta);
  lorentz boostToCoM;
  boostToCoM.setBeta(beta);
  //cout << beta << endl;
  
  for(int i=0;i<oneCluster.particleList.size();i++)
  {    
    
    if(i+1==oneCluster.particleList.size())
    { 
      break;
    }
    particles[oneCluster.particleList[i]].FLAVOR=up;
    particles[oneCluster.particleList[i+1]].FLAVOR=anti_up;
    
    particles[oneCluster.particleList[i]].Mom = boostToCoM.boost(particles[oneCluster.particleList[i]].Mom); 
//     
    
//     cout << particles[oneCluster.particleList[i]].Mom.E() << endl;
    ss5.str("");
    ss6.str("");
    int species = 1;
//     ss5 << species;
    int colorstream = 501+i;
    ss6 << colorstream;
    momentaString = momentaString + "1" + "\t+1\t3\t3\t"+ss6.str()+"\t0"+"\t";
    writeMomenta(particles[oneCluster.particleList[i]].Mom,momentaString);
    momentaString = momentaString +  "0" +"\t"+endOfLine;
    Nparts++;
    TotMom  += particles[oneCluster.particleList[i]].Mom;
    TotE    += particles[oneCluster.particleList[i]].Mom.E();    
    if(ParticlePrototype::mapToGenericFlavorType(particles[oneCluster.particleList[i]].FLAVOR)==light_quark)    TotMom_Q+=particles[oneCluster.particleList[i]].Mom;
    else TotMom_AQ+=particles[oneCluster.particleList[i]].Mom;
    ss5.str("");
    ss6.str("");

    
    i++;
    
    particles[oneCluster.particleList[i]].Mom = boostToCoM.boost(particles[oneCluster.particleList[i]].Mom);
    
    species = -1;
//     ss5 << species;
//     colorstream = 501;
    ss6 << colorstream;
    momentaString = momentaString + "-1" + "\t+1\t3\t3\t"+ "0\t" + ss6.str() +"\t";
    writeMomenta(particles[oneCluster.particleList[i]].Mom,momentaString);
    momentaString = momentaString + "0" + "\t" + endOfLine;
    Nparts++;
    TotMom  += particles[oneCluster.particleList[i]].Mom;
    TotE    += particles[oneCluster.particleList[i]].Mom.E();
    if(ParticlePrototype::mapToGenericFlavorType(particles[oneCluster.particleList[i]].FLAVOR)==light_quark)    TotMom_Q+=particles[oneCluster.particleList[i]].Mom;
    else TotMom_AQ+=particles[oneCluster.particleList[i]].Mom;    
  }

  Nparts += 3; // the incoming and intermediate particle count as well
  
  ss1 << TotMom_Q.E();
  string EnergyBeamStr = ss1.str();
  ss7 << TotMom_AQ.E();
  string EnergyBeamAntiStr = ss7.str();
  
  double invariantMass=sqrt(pow(TotMom.E(),2.0)-TotMom.vec2());
  ss3 << invariantMass;
  cout << "Tot Mom = " << TotMom << endl;
  
  cout << "Final Cluster Size " << Nparts-3 << endl;
  
  //cout << "TotE = " <<  TotE << "\t\t" << TotMom_AQ.E() +TotMom_Q.E() << endl;
  
  ss2 << Nparts;
  string NpartsStr=ss2.str();
  
//   cout << NpartsStr << endl;
  
//   double EnergyBeam=TotEQuark+TotEAQuark;
//   
   
   string lesHouchesString="<LesHouchesEvents version=\"1.0\">\n<headers>\n</headers>\n<init>\n11 -11 " + ss3.str() + " " + ss3.str() + " -1 -1 -1 -1 1 1\n0 0  3 1001\n</init>\n<event>\n";
   lesHouchesString+=NpartsStr+ "\t1001\t1\t-1\t0.00729927007\t0.3\n";
//

   VectorEPxPyPz TotEPlusMom;
   VectorEPxPyPz TotEMinusMom;
   
   TotEPlusMom=TotMom_Q;
   TotEMinusMom=TotMom_AQ;
   
//   EP=N.zeros(4)
//   EM=N.zeros(4)
//   
//   EP[2]=TotEQuark
//   EM[2]=-TotEAQuark
//   
//   EP[3]=TotEQuark
//   EM[3]=TotEAQuark
//   
//   //
//   def writeMomenta(FourVector):  
//   string=""
//   for momComp in FourVector:
//     string+=str(momComp)+"\t"
//   return string
//   //
//   
//   
//   invM=0
     lesHouchesString=lesHouchesString+"11\t-1\t0\t0\t0\t0\t";
     writeMomenta(TotEPlusMom,lesHouchesString);   
     lesHouchesString =  lesHouchesString + "0" +"\t"+endOfLine;
     
     
     lesHouchesString=lesHouchesString+"-11\t-1\t0\t0\t0\t0\t";
     writeMomenta(TotEMinusMom,lesHouchesString);
     lesHouchesString = lesHouchesString +"0" +"\t"+endOfLine;
//   
     
    VectorEPxPyPz Photon = TotMom;
    
    lesHouchesString=lesHouchesString+"22\t+2\t1\t2\t0\t0\t";
    writeMomenta(Photon,lesHouchesString);
    //invariant mass = energy beam
    lesHouchesString=lesHouchesString+ ss3.str()  +"\t"+endOfLine;
    
    lesHouchesString=lesHouchesString + momentaString;
    
    
//   invM=0
//   
//   #random color streams
//   #colorstream1=N.zeros(int(Nquarks/2))
//   #colorstream2=N.zeros(int(Nquarks/2))
//   
//   #for j in N.arange(int(Nquarks/2)):
//   #  colorstream1[j]=501+j
//   #  colorstream2[j]=501+j
//   
//   #random.shuffle(colorstream1)
//   #random.shuffle(colorstream2)
   
//   #Generate pair
//   for j in N.arange(Nquarks):
//   
//     Quark=N.zeros(4)
//     Quark[0]=pxQ[j]        #px
//     Quark[1]=pyQ[j]        #py
//     Quark[2]=pzQ[j]        #pz
//     Quark[3]=N.sqrt(pow(Quark[0],2.0)+pow(Quark[1],2.)+pow(Quark[2],2.))
//     
//     AQuark=N.zeros(4)
//     AQuark[0]=pxAQ[j]         #px
//     AQuark[1]=pyAQ[j]         #py
//     AQuark[2]=pzAQ[j]         #pz
//     AQuark[3]=N.sqrt(pow(AQuark[0],2.0)+pow(AQuark[1],2.)+pow(AQuark[2],2.))
//     
//     Antiquark=AQuark
//     
//     if(sampleIsotropic==True):
//       Quark=rotateMomentum(Quark)
//       Antiquark=BackToBack(Quark)
//     
//     species=1
//     colorstream=501+j         # back-to-back colors
//     #colorstream=colorstream1[j]
//     
//     lesHouchesString+=str(species)+"\t+1\t3\t3\t"+str(colorstream)+"\t"+str(0)+"\t"+writeMomenta(Quark)+str(invM) +"\t"+endOfLine
//     
//     
//     species=-1 
//     
//     colorstream=501+j
//     #colorstream=colorstream2[j]
//     lesHouchesString+=str(species)+"\t+1\t3\t3\t"+str(0)+"\t"+str(colorstream)+"\t"+writeMomenta(Antiquark)+str(invM) +"\t"+endOfLine
//     
//     
//   
//   
   lesHouchesString=lesHouchesString+"</event>\n</LesHouchesEvents>";
//   
//   
//   lesHouchesFormatFile   = "/home/greif/HERWIGNEW/generatedEvent.lhe" 
//   f = codecs.open(lesHouchesFormatFile,'w')
//   f.write(lesHouchesString)
//   f.close()  
   
   string filename;

   ss4 << no;
//   //+ ss.str() 
   filename = filename_prefix + "_generatedEvent_" + ss4.str() + ".lhe";
   fstream file_central( filename.c_str(), ios::out | ios::trunc  );
   file_central << lesHouchesString;
   file_central.close();
   
}
 
vector<string> analysis::split(string const &input) 
{ 
    std::istringstream buffer(input);
    vector<string> ret((std::istream_iterator<string>(buffer)), std::istream_iterator<std::string>());
    return ret;
} 
 
bool analysis::isUncharged(string const particleID)
{
  bool Uncharged=false;
  if(particleID=="22" || particleID=="111" || particleID=="113" || particleID=="-311" || particleID=="223" || particleID=="221")
  {
    Uncharged=true;
  }  
  return Uncharged;
}

bool analysis::isCharged(string particleID)
{
  if(isUncharged(particleID)) 
    return false;
  else
    return true; 
}
 
void analysis::readHadrons()
{
  
  for(int i= 0; i<ClusterNumbers.size();i++)
  {
    stringstream ss1;
    ss1 << ClusterNumbers[i];
    string filename = "myEventSaverun_" + ss1.str() + ".log";

//     //Count entries
//     std::ifstream countParticles( filename.c_str() );
//     if ( countParticles.good() )
//     {
//       int numberOfParticlesToGenerate=0;
//       numberOfParticlesToGenerate = std::count( std::istreambuf_iterator<char>( countParticles ), std::istreambuf_iterator<char>(), '\n' ); // counts number of lines in  particle data file  
//     }
   
    
    std::ifstream readParticles( filename.c_str() );
    cout << "Read hadronic particle file of single hadronized Cluster" << filename.c_str() << endl;
    int StartReading=0;
    int noEvents=0;
    std::string line;
    while( std::getline( readParticles, line ) )   
    {
      std::istringstream iss( line );

      std::string result;
      if( std::getline( iss, result) )
      { 
        if( (StartReading==3 || StartReading==2) && result == "------------------------------------------------------------------------------" )
        {
          StartReading=0;
          noEvents++;
          break; //break while loop of this Cluster
        } 
        if(StartReading==3)
        {
          vector<string> HadronLine;
          HadronLine = split(result); 
//           cout << "Momenta: ";
//           for(int i=0; i< HadronLine.size();i++)
//           {
//             cout << HadronLine[i] << "\t";
//           }
//           cout << endl;
          //-> Save momenta here

          
          Particle aNewHadron;
          
          
          aNewHadron.Mom.E()  = atof(HadronLine[3].c_str());
          aNewHadron.Mom.Px() = atof(HadronLine[0].c_str());
          aNewHadron.Mom.Py() = atof(HadronLine[1].c_str());
          aNewHadron.Mom.Pz() = atof(HadronLine[2].c_str());
          aNewHadron.m = atof(HadronLine[4].c_str());
          hadrons.push_back(aNewHadron);
          
          
          
          
          StartReading=2;
          continue;
        }
        
        if(StartReading==2)
        {
          vector<string> HadronLine;
          HadronLine = split(result); 
          if(HadronLine.size()==4)//line indicating the species
          {
            if(isCharged(HadronLine[2]))
            {
              StartReading=3;
//               cout << "found charged species " << HadronLine[1] << endl;
              continue;
  //             for(int i=0; i< HadronLine.size();i++)
  //             {
  //               cout << HadronLine[i] << "\t";
  //             }
  //             cout << endl;
            }else
            {
//               cout << "found uncharged species " << HadronLine[1] << endl;
              continue;
            }
          }else
          continue; //uncharged species momenta: skip line
        }
        

        

        if(StartReading==1 && result == "--- final:" )
        {
          cout << "--- final:" << endl;
          StartReading=2;
          continue;
        }
       
        if(StartReading==0 && result == "Step 4 performed by DecayHandler" )
        {
          cout << "Step 4 performed by DecayHandler" << endl;
          StartReading=1;
          continue;
        }
      }        
    }
    cout << "collected hadrons =  " << hadrons.size() << endl;
  }
}

//     for ( int i = 0; i < numberOfParticlesToGenerate; i++ )
//     {
//       int flavTemp;
//       double pX, pY, pZ, E;
// 
//       // structure of file
//       // number of pythia event  is hard?   flavour  energy  px  py  pz  m
//       readParticles >> flavTemp >> E >> pX >> pY >> pZ;
//       //cout <<  E << pX << pY << pZ << endl;
//     }
 
 
 
//   std::ifstream file( "sample.txt" );
// 
//     std::string line;
//     while( std::getline( file, line ) )   
//     {
//         std::istringstream iss( line );
// 
//         std::string result;
//         if( std::getline( iss, result , '=') )
//         {
//             if( result == "foo" )
//             {
//                 std::string token;
//                 while( std::getline( iss, token, ',' ) )
//                 {
//                     std::cout << token << std::endl;
//                 }
//             }
//             if( result == "bar" )
//             {
//                //...
//     }
//  
 
/*
def readHERWIGOutput(Seed):
  global eventNo,chargedHadrons
  
  filename = "myEPLUSEMINUS-S"+str(Seed)+".log"
  zeilen = []
  try:
    z = codecs.open(filename, mode='r')
    File_loaded = True
  except:
    print('File does not exist. Skip ' + filename)
    File_loaded = False
  if File_loaded == True:
    for line in z.readlines():
        zeilen.append(line)
    z.close()

    #print('Start the loop...')

    StartReading=0
    for line in zeilen:   
      
      if StartReading==3: 
        StartReading=4
        col = line.split()
        #print col
        if (HadronSpecies=="ChHadrons"):
          if not speciesID in unchargedHadrons:
            BinHerwigInSpectrum(float(col[3]),HadronSpecies)
        elif (HadronSpecies=="Pions"):  
          if speciesID=='211' or speciesID=='-211':
            #print("FoundPion E= " + col[3] )
            BinHerwigInSpectrum(float(col[3]),HadronSpecies)
            
      if StartReading==2:  
        col = line.split()
        StartReading=3
        try:
          speciesID=col[2]
        except:
          StartReading=0
          eventNo+=1

      if not line.find('Step 4 performed by DecayHandler')==-1:
        StartReading=1   
      if not line.find('--- final:')==-1 and StartReading==1:
        StartReading=2 
      if StartReading==4:
        StartReading=2
  print eventNo
  
  */
  
  
  
//   if ( file )
//   {
//     string line;
//     stringstream inStream;
//     double content;
// 
//     do
//     {
//       getline( file, line );
//       if ( line.find( "#", 0 ) == string::npos && line.length() != 0 )    //ignore empty lines and lines with leading "#"
//       {
//         if ( line.find( "-inf", 0 ) == string::npos )  //the string "-inf" can't be pushed into a vector<double>
//         {
//           inStream.str( "" );
//           inStream.clear();
// 
//           inStream.str( line );
//           inStream >> content;
// 
//           if ( !inStream.fail() )
//           {
//             I23data.push_back( content );
//           }
//           else
//           {
//             cout << "severe error in interpolation23::readTables() - ill-formated data file" << endl;
//             return I23_READOUT_ERROR;
//           }
//         }
//         else     //i.e. when string "-inf" is found
//         {
//           I23data.push_back( minus_infinity );  //push_back(minus_infinity) instead of "-inf"
//         }
// 
//       }
//     }
//     while ( !file.eof() );
//   }
//   else
//   {
//     cout << "No I23 data file existent for md2_counter_term_in_I23 = " << md2_counter_term_in_I23 << "  formationTimeTyp = " << formationTimeTyp << "  fudge_factor = " << fudge_factor_lpm << "  matrixElement23 = " << matrixElement23 << "  maxRunningCoupling = " << maxRunningCoupling << endl;
//     // throw an exception if the file could not have been read
//     throw eI23_read_error( "error in interpolation23::readTables() - could not open I23 file. File does not exist." );
//     return I23_FILE_ERROR;
//   }


 
void analysis::runHERWIG()
{
//   string url,wget="wget ";
//   cout << "bitte URL eingeben: ";
//   cin >> url;
//   wget += url;
//   system(wget.c_str());
  stringstream ss4,ss5;
  string filename,filenameFolder;
  ss5 << int(ran2()*100.);
  for(int clusterNumber=0; clusterNumber<NumberCluster;clusterNumber++)
  {
    if(theClusters[clusterNumber].particleList.size()>100 || theClusters[clusterNumber].particleList.size()<5) continue;
    cout << "RUN HERWIG FOR CLUSTER " << clusterNumber << endl;
    ClusterNumbers.push_back(clusterNumber);
    cout << ClusterNumbers.size() << endl;
    ss4.str("");
    ss4 << clusterNumber;
    filename = theConfig->getJobName() + "_generatedEvent_" + ss4.str() + ".lhe";
    filenameFolder = filename_prefix + "_generatedEvent_" + ss4.str() + ".lhe";

    string befehl="./GO.sh " + filename + " " + "myEventSaverun_" + ss4.str() + " " + filenameFolder;
    system(befehl.c_str());
    
    //WARNING: HACK
//     break;
  }
  
  cout << "D O N E - read in Herwig results now" << endl;
  
  
//             p = subprocess.Popen(["./GO.sh",str(Seed)])
//           p.communicate() 
//           readHERWIGOutput(Seed)
  
  
  
  
  
  
  
  
}
  
  
//   while(1)
//   {  
// 
// 
//       
// //       for ( it1 = theClusters[i].particleList.begin(); it1 != theClusters[i].particleList.end(); ++it1 )
// //       {
// //         AbstandSchwerpunkt=L.getSpatialDistance(particles[*it1].Pos,schwerpunkt);
// //         
// //         if (AbstandSchwerpunkt>largestDistanceToSchwerpunkt) 
// //         {
// //           largestDistanceToSchwerpunkt=AbstandSchwerpunkt;
// //           indexLargestDistParticle=*it1;       
// //         }
//         
// //         for ( it2 = theClusters[i].particleList.begin(); it2 != theClusters[i].particleList.end(); ++it2 )
// //         {
// //           if (*it1==*it2) continue;        
// //           D=L.getSpatialDistance(particles[*it1].Pos,particles[*it2].Pos);
// //           S+=D;
// //           if (D>largestDistance)
// //           {  
// //             largestDistance=D;
// //             indexBiggestCluster=i;
// //           }
// //           //largestDistance=L.getSpatialDistance(particles[*it1].Pos,particles[*it2].Pos);
// //   //         cout << L.getSpatialDistance(particles[*it1].Pos,particles[*it2].Pos) << endl;
// //         }
//       }
//       //S/=theClusters[i].particleList.size();
//       
//       //cout << particles[indexLargestDistParticle].Pos << endl;
//       //S += theClusters[i].particleList.size(); 
//     }
//     
//     
//     
//     
//     
//     cout << "S new = " << S << "     S_old = " << S_old << endl;
//     if (S_old<S)
//     {
//       break;
//       cout << "BREAK LOOP" << endl;
//     }
//     else
//     {
//       S_old=S;
//       S=0;
//       cout << "largestDistance = " << largestDistanceToSchwerpunkt << endl;
//       cout << "indexBiggestCluster = " << indexBiggestCluster << endl;
//       cout << "indexLargestDistParticle = " << indexLargestDistParticle << endl;
//       //Split Cluster
//       cluster newCluster;
//       std::list<int>::iterator it3;
//       for ( it3 = theClusters[indexBiggestCluster].particleList.begin(); it3 != theClusters[indexBiggestCluster].particleList.end(); ++it3 )
//       {
//         double dist_Schwerpkt = L.getSpatialDistance(schwerpunkt,particles[*it3].Pos);
//         double dist_NewCenter = L.getSpatialDistance(particles[indexLargestDistParticle].Pos,particles[*it3].Pos);
//         if (dist_NewCenter < dist_Schwerpkt) //Sort in new cluster
//         {
//           newCluster.particleList.push_back(static_cast<int>(*it3));
//           //theClusters[indexBiggestCluster].particleList.remove(static_cast<int>(*it3));
//         }
//       
//       } 
//       cout << "New Cluster with " << newCluster.particleList.size() << " particles" << endl;
//       theClusters.push_back(newCluster);  
//     }
//   } 
//   
//   for ( unsigned int i = 0; i < particles.size(); i++ )
//   {
// 
//     double eta = particles[i].Mom.Rapidity();
//     //for ( int yRangeIndex = 0; yRangeIndex < rapidityRanges.size(); yRangeIndex++ )
//     int yRangeIndex = 0;//Midrapidity
//     if ( fabs( eta ) >= rapidityRanges[yRangeIndex].yleft && fabs( eta ) <= rapidityRanges[yRangeIndex].yright )
//     {
//       px = particles[i].Mom.Px(); 
//       py = particles[i].Mom.Py();
//       pz = particles[i].Mom.Pz();
//       file << px << "\t" << py << "\t" << pz << endl;
//     }    
//   }
//   file.close();





// // writes temperature and velocity of all specified cells in a file, used by Moritz
// void analysis::writeTempCustom( const int step  )
// {
//   
//   int nx, ny, nz, cell_id;
//   double time, pr, XT;
//  
//   const int minNmbTemp = 50; // minimum number of particles to calculate temperature from
// //   const int minNmbTemp = 2; // minimum number of particles to calculate temperature from
//   const int minNmbTempCell = 10; // minimum number of particles in one cell to calculate a Temperature. If the number is below this value, the cell is taken as empty (no temperature)
//   
//   binning binNumber("output/number.dat", -0.5, 49.5, 50);
//   binning binTempCells("output/tempCells.dat", 0.0, 3., 70);
//   binning binTempRings("output/tempRings.dat", 0.0, 2.2, 70);
//   
//   int count_tempCell = 0;
//   int count_tempNeighborCells = 0;
//   int count_noTemp = 0;
//   
//   int IXY = IX * IY;
// 
//   if ( step != 0 && step != nTimeSteps )
//     time = tstep[step-1];
//   else
//     return;
//   
//  
//   // total length of grid system
//   const double xlength = 25;//24.6;
//   const double ylength = 25;//24.6;
//   const double zlength = 13;//12.3;
// 
//   // number of cells in given direction
//   // in cascade IX=40   IY=40   IZ=47
//   const int nCellsx = 10;
//   const int nCellsy = 10;
//   const int nCellsz = 50;
//   nCells = nCellsx * nCellsy * nCellsz;
//   
//   dv = xlength * ylength * zlength / nCells; // volume of each cell
// 
// //   int numberInCell[nCells]; // number of all particles in cell
// //   double vx_cell[nCells]; // total x-velocity of all particles in cell
// //   double vy_cell[nCells];  
// //   double vz_cell[nCells];
// //   double vr_cell[nCells];
// //   double em_cell[nCells];// total energy of all particles in cell
// //   double prm_cell[nCells];
// //   double pzm_cell[nCells];
// //   double pr2em_cell[nCells];
// //   double pz2em_cell[nCells];
// //   double przem_cell[nCells];
// //   double densn_cell[nCells];
// //   double gama_cell[nCells];
// //   double temp_cell[nCells];
// //   double tempWithQuarks_cell[nCells];
// 
// 
//   numberInCell = new int[nCells]; // number of all particles in cell
//   numberInCellQuarks = new int[nCells]; // number of all particles in cell
//   numberInCellGluons = new int[nCells]; // number of all particles in cell
//   fugacityGluons = new double[nCells]; 
//   fugacityQuarks = new double[nCells]; 
//   
//   
//   temp_numberInCell = new int[nCells]; // number of all particles in cell
//   vx_cell = new double[nCells]; // total x-velocity of all particles in cell
//   vy_cell = new double[nCells];  
//   vz_cell = new double[nCells];
//   vr_cell = new double[nCells];
//   em_cell = new double[nCells];// total energy of all particles in cell
//   prm_cell = new double[nCells];
//   pzm_cell = new double[nCells];
//   pr2em_cell = new double[nCells];
//   pz2em_cell = new double[nCells];
//   przem_cell = new double[nCells];
//   densn_cell = new double[nCells];
//   gama_cell = new double[nCells];
//   temp_cell = new double[nCells];
//   tempWithQuarks_cell = new double[nCells];
//   
//   // for copy. If too few particles in one cell, the particles from the surrounding cells are added from _org. Otherwise one would double or triple add particles from cells
//   numberInCell_org = new int[nCells]; // number of all particles in cell
//   vx_cell_org = new double[nCells]; // total x-velocity of all particles in cell
//   vy_cell_org = new double[nCells];  
//   vz_cell_org = new double[nCells];
//   vr_cell_org = new double[nCells];
//   em_cell_org = new double[nCells];// total energy of all particles in cell
//   prm_cell_org = new double[nCells];
//   pzm_cell_org = new double[nCells];
//   pr2em_cell_org = new double[nCells];
//   pz2em_cell_org = new double[nCells];
//   przem_cell_org = new double[nCells];
// // 
// 
// //   for ( int i = 0; i < nCells; i++ )
// //   {
// //     const int nxny = nCellsx * nCellsy;
// //     const int indexZ = i / nxny;
// //     const int indexY = (i -  indexZ  * nxny) / nCellsx;
// //     const int indexX = i -  indexZ  * nxny - indexY * nCellsx;
// //     
// //     energy[i] = 0.0;
// //     numberInCell[i] = indexX + nCellsx * indexY + nCellsx * nCellsy * indexZ;
// //     vx[i] = double(indexX) * xlength / nCellsx;
// //     vy[i] = double(indexY) * ylength / nCellsy;
// //     vz[i] = double(indexZ) * zlength / nCellsz;
// //   }
// 
//   // set all properties to 0
//   for ( int i = 0; i < nCells; i++ )
//   {
//     numberInCell[i] = 0;
//     numberInCellQuarks[i] = 0;
//     numberInCellGluons[i] = 0;
//     fugacityQuarks[i]=0;
//     fugacityGluons[i]=0;
//     temp_numberInCell[i] = 0;
//     vx_cell[i] = 0.0; 
//     vy_cell[i] = 0.0;  
//     vz_cell[i] = 0.0;
//     vr_cell[i] = 0.0;
//     em_cell[i] = 0.0;
//     prm_cell[i] = 0.0;
//     pzm_cell[i] = 0.0;
//     pr2em_cell[i] = 0.0;
//     pz2em_cell[i] = 0.0;
//     przem_cell[i] = 0.0;
//     densn_cell[i] = 0.0;
//     gama_cell[i] = 0.0;
//     temp_cell[i] = 0.0;
//     tempWithQuarks_cell[i] = 0.0;
//     
//     numberInCell_org[i] = 0;
//     vx_cell_org[i] = 0.0; 
//     vy_cell_org[i] = 0.0;  
//     vz_cell_org[i] = 0.0;
//     vr_cell_org[i] = 0.0;
//     em_cell_org[i] = 0.0;
//     prm_cell_org[i] = 0.0;
//     pzm_cell_org[i] = 0.0;
//     pr2em_cell_org[i] = 0.0;
//     pz2em_cell_org[i] = 0.0;
//     przem_cell_org[i] = 0.0;
//   }
//   
//   // sum over all particles
//   for ( int i = 0; i < particles_atTimeNow.size(); i++ )
//   {
//     if ( FPT_COMP_E( particles_atTimeNow[i].Pos.T(), time ) && particles_atTimeNow[i].FLAVOR < 7 ) // only gluons and light quarks
//     {
//       // determine cell id
//       if ( fabs( particles_atTimeNow[i].Pos.X() - xlength / 2.0 ) < 1.0e-6 )
//         nx = nCellsx - 1;
//       else
//         nx = int(( particles_atTimeNow[i].Pos.X() / xlength + 0.5 ) * nCellsx );
// 
//       if ( fabs( particles_atTimeNow[i].Pos.Y() - ylength / 2.0 ) < 1.0e-6 )
//         ny = nCellsy - 1;
//       else
//         ny = int(( particles_atTimeNow[i].Pos.Y() / ylength + 0.5 ) * nCellsy );
// 
//       if ( fabs( particles_atTimeNow[i].Pos.Z() - zlength / 2.0 ) < 1.0e-6 )
//         nz = nCellsz - 1;
//       else
//         nz = int(( particles_atTimeNow[i].Pos.Z() / zlength + 0.5 ) * nCellsz );
// 
// 
//       if (( nx >= nCellsx ) || ( nx < 0 ) || ( ny >= nCellsy ) || ( ny < 0 ) || ( nz >= nCellsz ) || ( nz < 0 ) )
//       {
//         /*cout << "err cell_ID in temp output" << endl;
//         cout << particles_atTimeNow[i].Pos.T() << "\t" << particles_atTimeNow[i].Pos.X() << "\t" << particles_atTimeNow[i].Pos.Y();
//         cout << "\t" << particles_atTimeNow[i].Pos.Z() << endl;
//         cout << nx << "\t" << ny << "\t" << nz << endl;*/
//       }
//       else
//       {
//         cell_id = nx + nCellsx * ny + nCellsx * nCellsy * nz;
//         
//         if(particles_atTimeNow[i].FLAVOR != gluon )
//         {
//           ++numberInCellQuarks[cell_id];
//         }
//         if(particles_atTimeNow[i].FLAVOR == gluon)
//         {
//           ++numberInCellGluons[cell_id];
//         }        
//         ++numberInCell[cell_id];
//         
//         XT = particles_atTimeNow[i].Pos.Perp();
//         if ( XT < 1.0e-5 )
//         {
//           pr = particles_atTimeNow[i].Mom.Pt();
//         }
//         else
//         {
//           pr = ( particles_atTimeNow[i].Mom.Px() * particles_atTimeNow[i].Pos.X()
//                  + particles_atTimeNow[i].Mom.Py() * particles_atTimeNow[i].Pos.Y() ) / XT;
//         }
//         vr_cell[cell_id] += pr / particles_atTimeNow[i].Mom.E();
//         vx_cell[cell_id] += particles_atTimeNow[i].Mom.Px() / particles_atTimeNow[i].Mom.E();
//         vy_cell[cell_id] += particles_atTimeNow[i].Mom.Py() / particles_atTimeNow[i].Mom.E();
//         vz_cell[cell_id] += particles_atTimeNow[i].Mom.Pz() / particles_atTimeNow[i].Mom.E();
//         
//         em_cell[cell_id] += particles_atTimeNow[i].Mom.E();
//         prm_cell[cell_id] += pr;
//         pzm_cell[cell_id] += particles_atTimeNow[i].Mom.Pz();
//         pr2em_cell[cell_id] += pr * pr / particles_atTimeNow[i].Mom.E();
//         pz2em_cell[cell_id] += particles_atTimeNow[i].Mom.Pz() * particles_atTimeNow[i].Mom.Pz() / particles_atTimeNow[i].Mom.E();
//         przem_cell[cell_id] += pr * particles_atTimeNow[i].Mom.Pz() / particles_atTimeNow[i].Mom.E();
//         
//         // temp of particles summed
//         tempWithQuarks_cell[cell_id] += particles_atTimeNow[i].temperature;
//         if( particles_atTimeNow[i].temperature >= 0.1)
//           ++temp_numberInCell[cell_id];
//       }
//     }
//   }
//   
//   // duplicate properties of cell
//   for ( int i = 0; i < nCells; i++ )
//   {
//     numberInCell_org[i] = numberInCell[i];
//     vr_cell_org[i] = vr_cell[i];
//     vx_cell_org[i] = vx_cell[i];
//     vy_cell_org[i] = vy_cell[i];
//     vz_cell_org[i] = vz_cell[i];
//     
//     em_cell_org[i] = em_cell[i];
//     prm_cell_org[i] = prm_cell[i];
//     pzm_cell_org[i] = pzm_cell[i];
//     pr2em_cell_org[i] = pr2em_cell[i];
//     pz2em_cell_org[i] = pz2em_cell[i];
//     przem_cell_org[i] = przem_cell[i];
//   }
// 
//   
//   // determine cells which do not have enough particles
//   for ( int i = 0; i < nCells; i++ )
//   {
//     // If less than minNmbTempCell particles are in the cell it is taken as empty
//     if(numberInCell[i] >= minNmbTempCell)
//     {
//       // If less than minNmbTemp particles are in the cell the temperature the 6 neighbor cells are also taken into account
//       if(numberInCell[i] < minNmbTemp)
//       {
//         int cell_id_neighbor;
//         
//         //neighbors in x direction
//         cell_id_neighbor = i-1;
//         addNeighborCells( i, cell_id_neighbor );
//         cell_id_neighbor = i+1;
//         addNeighborCells( i, cell_id_neighbor );
//         //neighbors in y direction
//         cell_id_neighbor = i+IX;
//         addNeighborCells( i, cell_id_neighbor );
//         cell_id_neighbor = i-IX;
//         addNeighborCells( i, cell_id_neighbor );
//         //neighbors in z direction
//         cell_id_neighbor = i-IXY;
//         addNeighborCells( i, cell_id_neighbor );
//         cell_id_neighbor = i+IXY;
//         addNeighborCells( i, cell_id_neighbor );
//       }
//     }
//   }
//   
//   // calculate temperature for each cell
//   for ( int i = 0; i < nCells; i++ )
//   {
//     // If less than minNmbTempCell particles are in the cell it is taken as empty
//     if(numberInCell[i] >= minNmbTempCell)
//     {
//       // If more than minNmbTemp particles are in the cell the temperature is calculated by them
//       if(numberInCell[i] >= minNmbTemp)
//       {
//         calcTempCell( i );
//         
//         tempWithQuarks_cell[i] = tempWithQuarks_cell[i]/temp_numberInCell[i];
// 
//         binTempCells.add(temp_cell[i]);
//         count_tempCell++;
//       }
//       else // take also 6 neighbor cells into account
//       {
//         count_noTemp++;
//         
//         temp_cell[i] = 0.0; // not enough particles in surrounding to calculate temperature
//         
//         vx_cell[i] = 0.0;
//         vy_cell[i] = 0.0;
//         vz_cell[i] = 0.0;
//       }
//     }
//     else
//     {
//       if(numberInCell[i] == 0)
//       {
//         temp_cell[i] = -2.0; // completely empty
//         vx_cell[i] = 0.0;
//         vy_cell[i] = 0.0;
//         vz_cell[i] = 0.0;
//       }
//       else
//       {
//         temp_cell[i] = -1.0; // very few particles
//         vx_cell[i] = 0.0;
//         vy_cell[i] = 0.0;
//         vz_cell[i] = 0.0;
//       }
//     }
//   }
// 
//   double avFQ = 0.;
//   double avFG = 0.;
//   double avT = 0.;
//   for(int i=1;i<=nCells;i++)
//   {
//     avFQ += fugacityQuarks[i];  
//     avFG += fugacityGluons[i];
//     avT += temp_cell[i];
//   }
//   avFQ /=count_tempCell; 
//   avFG /=count_tempCell;
//   avT /=count_tempCell;
//   
//   
//   string filename,name;
//   stringstream ss;
//   ss << time*10;
//   //+ ss.str() 
//   filename = filename_prefix + "_central" + ".dat";
//   fstream file_central( filename.c_str(), ios::out | ios::app  );
//   
//   const int middle_cell_id = nCellsx/2 + nCellsx * nCellsy/2 +nCellsx*nCellsy*nCellsz/2;
//   cout  << time <<'\t'<< "T   " << '\t' << "FUgQ   " << '\t'<< "FugG  " << '\t' << "nQUARKS  " << '\t'<< "nGluons  "  << '\t' << "nTot   " << endl;
//   cout  << time <<'\t'<< temp_cell[middle_cell_id] << '\t' << fugacityQuarks[middle_cell_id] << '\t'<<fugacityGluons[middle_cell_id] << '\t' << numberInCellQuarks[middle_cell_id] << '\t'<< numberInCellGluons[middle_cell_id]<< '\t' << numberInCell[middle_cell_id] << '\t'<< numberInCellQuarks[middle_cell_id]+numberInCellGluons[middle_cell_id]  << endl;
// 
//   file_central << time <<'\t'<<avT<<'\t'<<avFQ<<'\t'<<avFG<<'\t'<< temp_cell[middle_cell_id] << '\t' << fugacityQuarks[middle_cell_id] << '\t'<<fugacityGluons[middle_cell_id] << '\t' << numberInCellQuarks[middle_cell_id] << '\t'<< numberInCellGluons[middle_cell_id]<< '\t' << numberInCell[middle_cell_id] << endl;
// 
//   
//   delete[] numberInCell; 
//   delete[] vx_cell; 
//   delete[] vy_cell;  
//   delete[] vz_cell;
//   delete[] vr_cell;
//   delete[] em_cell;
//   delete[] prm_cell;
//   delete[] pzm_cell;
//   delete[] pr2em_cell;
//   delete[] pz2em_cell;
//   delete[] przem_cell;
//   delete[] densn_cell;
//   delete[] gama_cell;
//   delete[] temp_cell;
//   delete[] tempWithQuarks_cell;
// }


void analysis::v2_output()
{
  //---- get time and calculate real time needed by simulation ----
  time_t end;
  time( &end );
  int secs = ( int )difftime( end, start );
  int hours = secs / 3600, minutes = ( secs % 3600 ) / 60, seconds = secs - 3600 * hours - 60 * minutes;
  cout << hours << "h" << minutes << "m" << seconds << "s" <<  endl;
  //---------------------------------------------------------------
  
  string filename, filename2;
  filename = theConfig->getStandardOutputDirectoryName() + "/" +  theConfig->getJobName() + "_v2_pt";
  filename2 = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + "_mean_v2";
  
  //-------------------------------------------------------------------
  fstream file( filename.c_str(), ios::out | ios::trunc );

  //---- print header if file is empty ----
  file.seekp( 0, ios::end );
  long size = file.tellp();
  file.seekp( 0, ios::beg );
  if ( size == 0 )
    printHeader( file, v2output, end );
  //---------------------------------------
  
  double temp_v2;
  
  for ( int y_index = nRapidityBins - 1; y_index >= 0; --y_index )
  {
    file << "#y in [ -" << rapidityBins[y_index] << ", " << rapidityBins[y_index] << " ]" << endl;
    
    for ( int pt_index = 0; pt_index < nBins_v2_pt; pt_index++ )
    {
      file << differential_v2_bins_labels[pt_index] << sep;
      for ( int time_index = 0; time_index <= nTimeSteps; time_index++ )
      {
        if ( differential_v2_counts[time_index][y_index][pt_index] > 0 )
        {
          temp_v2 = differential_v2_values[time_index][y_index][pt_index] / differential_v2_counts[time_index][y_index][pt_index];
        }
        else
        {
          temp_v2 = 0;
        }
        file << temp_v2 << sep;
      }
      file << endl;
    }
    file << endl << endl;
  }
  
  file.close();
  //-------------------------------------------------------------------
  
  
  //-------------------------------------------------------------------
  fstream file2( filename2.c_str(), ios::out | ios::trunc );

  //---- print header if file is empty ----
  file2.seekp( 0, ios::end );
  size = file2.tellp();
  file2.seekp( 0, ios::beg );
  if ( size == 0 )
    printHeader( file2, meanv2output, end );
  //---------------------------------------
  
  double temp_int_v2, temp_mean_pt;
  
  file2 << "# [-y,y]    <pt>(t_i)   <v2>(t_i)" << endl;
  file2 << "# [-y,y]" << sep;
  for ( int time_index = 0; time_index <= nTimeSteps; time_index++ )
  {
    file2 << "t_" << time_index << " = " << tstep[time_index] << sep;
  }
  file2 << endl;
  
  for ( int y_index = nRapidityBins - 1; y_index >= 0; --y_index )
  {
    file2 << rapidityBins[y_index] << sep;
    for ( int time_index = 0; time_index <= nTimeSteps; time_index++ )
    {
      if ( integrated_v2_counts[time_index][y_index] > 0 )
      {
        temp_int_v2 = integrated_v2[time_index][y_index] / integrated_v2_counts[time_index][y_index];
        temp_mean_pt = mean_PT[time_index][y_index] / integrated_v2_counts[time_index][y_index];
      }
      else
      {
        temp_int_v2 = 0;
        temp_mean_pt = 0;
      }
    
      file2 << temp_mean_pt << sep << temp_int_v2 << sep; 
    }
    file2 << endl;
  }
  
  file2.close();
  //-------------------------------------------------------------------
  
}






// compute v2 
void analysis::computeV2RAA( string name, const double _outputTime )
{
  v2RAA theV2RAA( theConfig, name, filename_prefix, rapidityRanges );
  
  theV2RAA.setPtBinProperties( 0, 8, 20 ); 
  theV2RAA.setPtBinsLinear(0, 8, 16 );
  theV2RAA.computeFor( gluon, particles, particles.size(), "bulk", _outputTime );
  theV2RAA.computeFor( light_quark, particles, particles.size(), "bulk", _outputTime );
  
  /*
  theV2RAA.computeFor( up, particles, particles.size(), "bulk", _outputTime );
  theV2RAA.computeFor( down, particles, particles.size(), "bulk", _outputTime );
  theV2RAA.computeFor( strange, particles, particles.size(), "bulk", _outputTime );
  theV2RAA.computeFor( anti_up, particles, particles.size(), "bulk", _outputTime );
  theV2RAA.computeFor( anti_down, particles, particles.size(), "bulk", _outputTime );
  theV2RAA.computeFor( anti_strange, particles, particles.size(), "bulk", _outputTime );
  theV2RAA.computeFor( light_quark, particles, particles.size(), "bulk", _outputTime );
  */
  
  
}



v2RAA::v2RAA( config * const c, string name_arg, string filename_prefix_arg, std::vector<analysisRapidityRange> rapidityRanges_arg, const double pt_min_arg, const double pt_max_arg, const int n_g_arg ):
    theConfig( c ), 
    pt_min( pt_min_arg ), pt_max( pt_max_arg ), 
    n_g( n_g_arg ),
    name( name_arg ), 
    filename_prefix( filename_prefix_arg ), 
    rapidityRanges( rapidityRanges_arg )
{
  eta_bins = rapidityRanges.size();
}




void v2RAA::computeFor( const FLAVOR_TYPE _flavTypeToComputeFor, vector<Particle>& _particles, const int n_particles, string additionalNameTag, const double _outputTime )
{
  double eta, pt, v2, xt;
  double sinAlpha, alpha;
  int dummy,dummy_lin, flavor, n_bins;
  int alphaIndex;
  
  double _pt_min, _pt_max;

  string filename_v2, filename_v2_summed, filename_v2_tot, filename_yield, filename_yield_lin, filename_pt_angleDependence, type,filename_averages;
  
  type = Particle::getName( _flavTypeToComputeFor );
  
  n_bins = n_g;
  _pt_max = pt_max;
  _pt_min = pt_min;
  
  // avoid problem with binning pt logaritmitically: cannot deal with pt = 0
  if( _pt_min < 0.05 )
    _pt_min = 0.05;

  
  const double d_ln_pt = ( log( _pt_max ) - log( _pt_min ) ) / n_bins;
  const double d_pt = ( pt_max_lin -  pt_min_lin )/n_g_lin;

  
  
  double v2sum[eta_bins];
  int NmbInRange[eta_bins];
  double MeanPT[eta_bins];
  int NmbInnerRegion = 0;
  for ( int j = 0;j < eta_bins;j++ )
  {
    v2sum[j] = 0.0;
    NmbInRange[j] = 0;
    MeanPT[j]=0;
  }

  
  double ptBinsV2[eta_bins][n_bins+1];
  
  int ptBinsNmb[eta_bins][n_bins+1];
  int ptBinsNmbLin[eta_bins][n_g_lin+1];
  
  int ptBinsInnerRegion[n_bins+1];
  for ( int j = 0;j < n_bins + 1;j++ )
  {
    ptBinsInnerRegion[j] = 0;
    for ( int i = 0;i < eta_bins;i++ )
    {
      ptBinsV2[i][j] = 0.0;
      ptBinsNmb[i][j] = 0.0;
      ptBinsNmbLin[i][j] = 0.0;
    }
  }
  
  const double deltaAlpha = 15; // degrees
  const int nAlphaBins = 6;  // 90° / 15°
  double ptBinsAngleDep[eta_bins][nAlphaBins][n_bins+1];
  for ( int i = 0; i < eta_bins; i++ )
  {
    for ( int j = 0; j < nAlphaBins; j++ )
    {
      for ( int k = 0; k < n_bins + 1; k++ )
      {
        ptBinsAngleDep[i][j][k] = 0;
      }
    }
  }
  

  // compute v2 and bin it into pt bins
  for ( int i = 0; i < n_particles; i++ )
  {
    pt = _particles[i].Mom.Pt();
    xt = _particles[i].Pos.Perp();

    sinAlpha = _particles[i].Mom.Py() / pt;
    alpha = asin( fabs( sinAlpha ) );
    alpha = alpha * 180 / M_PI;
    
    alphaIndex = static_cast<int>( alpha / deltaAlpha );
    if ( alphaIndex >= nAlphaBins )
    {
      alphaIndex = nAlphaBins - 1;
    }
        
    // for most particle species we are interested in the pseudorapidity (for massless particle it does not matter anyhow)
    eta = _particles[i].Mom.Pseudorapidity( _particles[i].m );
    
    // for some scenarios however explicitly the rapidity is measured. So substitute eta by the rapidity:
    if( Particle::mapToGenericFlavorType( _flavTypeToComputeFor ) == dmeson_gen ||
        Particle::mapToGenericFlavorType( _flavTypeToComputeFor ) == bmeson_gen ||
        Particle::mapToGenericFlavorType( _flavTypeToComputeFor ) == jpsi
    )
      eta = _particles[i].Mom.Rapidity();

    v2 = ( pow( _particles[i].Mom.Px(), 2.0 ) - pow( _particles[i].Mom.Py(), 2.0 ) ) / pow( pt, 2.0 );

    flavor = _particles[i].FLAVOR;
    FLAVOR_TYPE genFlavor = Particle::mapToGenericFlavorType( static_cast<FLAVOR_TYPE>( flavor ) );
    
    if( ( _flavTypeToComputeFor == flavor ) || 
          ( _flavTypeToComputeFor == light_quark && ( genFlavor == light_quark || genFlavor == anti_light_quark ) ) ||
          ( _flavTypeToComputeFor == genFlavor )||
          ( _flavTypeToComputeFor == charm && ( flavor == anti_charm ) ) ||
          ( _flavTypeToComputeFor == bottom && ( flavor == anti_bottom ) ) ||
          ( _flavTypeToComputeFor == heavy_quark && ( flavor == charm || flavor == bottom || flavor == anti_charm || flavor == anti_bottom ) )
        ) 
    {
      // individually check for each rapidity range whether this particle needs to be binned
      for ( int yRangeIndex = 0; yRangeIndex < eta_bins; yRangeIndex++ )
      {
        if ( fabs( eta ) >= rapidityRanges[yRangeIndex].yleft && fabs( eta ) <= rapidityRanges[yRangeIndex].yright )
        {
          v2sum[yRangeIndex] += v2;
          NmbInRange[yRangeIndex]++;
          MeanPT[yRangeIndex]+=pt;
          if (pt <= pt_max_lin && pt > pt_min_lin)
          {
            dummy_lin = int((pt-pt_min_lin)/d_pt); //for linear binning
            ptBinsNmbLin[yRangeIndex][dummy_lin]++;            
          }
          if ( pt <= _pt_max && pt > _pt_min )
          {
            dummy = int(( log( pt ) - log( _pt_min ) ) / d_ln_pt );
            ptBinsV2[yRangeIndex][dummy] += v2;                       
            ptBinsNmb[yRangeIndex][dummy]++;           
            ptBinsAngleDep[yRangeIndex][alphaIndex][dummy]++;
          }
        }
      }
    }
  }

  int binMax = 0;
  int binMin = n_particles;
  for ( int k = 0; k < n_bins + 1; k++ )
  {
    if ( ptBinsNmb[0][k] > binMax )
    {
      binMax = ptBinsNmb[0][k];
    }
    if ( ptBinsNmb[0][k] < binMin && ptBinsNmb[0][k] != 0 )
    {
      binMin = ptBinsNmb[0][k];
    }
  }

  // file output
  double pt_out;

  filename_v2_summed = filename_prefix + "_" + type + "_" + additionalNameTag + "_v2_pt_summed_" + name;
  filename_v2_tot = filename_prefix + "_" + type + "_" + additionalNameTag + "_v2_tot_" + name;
  filename_yield = filename_prefix + "_" + type + "_" + additionalNameTag + "_yield_pt_" + name;
  filename_yield_lin = filename_prefix + "_" + type + "_" + additionalNameTag + "_yield_pt_lin_" + name;  
  filename_pt_angleDependence = filename_prefix + "_" + type + "_" + additionalNameTag + "_pt_angular_dependence_" + name;
  filename_averages = filename_prefix + "_" + type + "_" + additionalNameTag + "_averages";// + name;
  
  
  fstream print_v2_summed( filename_v2_summed.c_str(), ios::out | ios::trunc );
  fstream print_v2_tot( filename_v2_tot.c_str(), ios::out | ios::trunc );
  fstream print_yield( filename_yield.c_str(), ios::out | ios::trunc );
  fstream print_yield_lin( filename_yield_lin.c_str(), ios::out | ios::trunc );
  fstream print_pt_angleDependence( filename_pt_angleDependence.c_str(), ios::out | ios::trunc );
  fstream print_averages( filename_averages.c_str(), ios::out | ios::app );
  
  //print _averages  
  //print_averages << "Mean-PT for different rapidity bins and dN/dy." << endl;
  print_averages << _outputTime;
  for ( int yRangeIndex = 0;yRangeIndex < eta_bins;yRangeIndex++ )
  { 
    double delta_eta = 2.0 * ( rapidityRanges[yRangeIndex].yright - rapidityRanges[yRangeIndex].yleft );
    if ( NmbInRange[yRangeIndex] > 0 )
    {
      MeanPT[yRangeIndex] /= NmbInRange[yRangeIndex];
    }
    else
    {
      MeanPT[yRangeIndex] = 0.;
    }
    print_averages << "\t" << MeanPT[yRangeIndex] << "\t" << NmbInRange[yRangeIndex]/(theConfig->getTestparticles())/delta_eta;  
  }  
  print_averages << endl;
  
  
  
  // print total v2
  print_v2_tot << "# total v2 of " << type << endl;
  print_v2_tot << "# t = " << _outputTime << endl;
  print_v2_tot << "# bin statistics for 0.35 mid-rapidity:  Avg per bin=" << double( NmbInRange[0] ) / n_c << "   Min=" << binMin << "   Max=" << binMax << endl;
  print_v2_tot << "# total v2, v2_sum and number in range for different rapidity bins" << endl;

  print_v2_tot << _pt_min;
  print_v2_tot.width( 15 );
  for ( int i = 0;i < eta_bins;i++ )
  {
    if ( NmbInRange[i] > 0 )
    {
      print_v2_tot <<  v2sum[i] / NmbInRange[i];
    }
    else
    {
      print_v2_tot <<  0;
    }
    print_v2_tot.width( 15 );
  }
  for ( int i = 0;i < eta_bins;i++ )
  {
    print_v2_tot <<  v2sum[i];
    print_v2_tot.width( 15 );
    print_v2_tot <<  NmbInRange[i];
    print_v2_tot.width( 15 );
  }
  print_v2_tot << endl;

  print_v2_tot << _pt_max;
  print_v2_tot.width( 15 );
  for ( int i = 0;i < eta_bins;i++ )
  {
    if ( NmbInRange[i] > 0 )
    {
      print_v2_tot <<  v2sum[i] / NmbInRange[i];
    }
    else
    {
      print_v2_tot <<  0;
    }
    print_v2_tot.width( 15 );
  }
  for ( int i = 0;i < eta_bins;i++ )
  {
    print_v2_tot <<  v2sum[i];
    print_v2_tot.width( 15 );
    print_v2_tot <<  NmbInRange[i];
    print_v2_tot.width( 15 );
  }
  print_v2_tot << endl;


  // print summed output, v2 is not computed, but summed v2 and the number in one bin
  print_v2_summed << "# summed v2 of " << type << endl;
  print_v2_summed << "# t = " << _outputTime << endl;
  print_v2_summed << "#";
  print_v2_summed.width( 14 );
  print_v2_summed << "pt       summed v_2 and number in bin for different rapidity bins" << endl;

  for ( int k = 0;k < n_bins + 1;k++ )
  {
    print_v2_summed.width( 15 );
    pt_out = exp( double( k ) * d_ln_pt + log( _pt_min ) + d_ln_pt / 2.0 ); // +dpt/2.0 to shift value in the middle of the intervall, important for small number of bins
    print_v2_summed << pt_out;
    for ( int i = 0;i < eta_bins;i++ )
    {
      print_v2_summed.width( 15 );
      print_v2_summed << ptBinsV2[i][k];
      print_v2_summed.width( 10 );
      print_v2_summed << ptBinsNmb[i][k];
    }
    print_v2_summed << endl;
  }

  // print yield log bins
  print_yield << "# " << type << " yield distribution" << endl;
  print_yield << "# t = " << _outputTime << endl;
  print_yield << "#";
  print_yield.width( 14 );
  print_yield << "pt       yield for different rapidity bins" << endl;

  for ( int k = 0;k < n_bins + 1;k++ )
  {
    print_yield.width( 15 );
    pt_out = exp( double( k ) * d_ln_pt + log( _pt_min ) + d_ln_pt / 2.0 ); // +dpt/2.0 to shift value in the middle of the intervall, important for small number of bins
    print_yield << pt_out;
    const double dpt = exp( double( k ) * d_ln_pt + log( _pt_min ) + d_ln_pt ) - exp( double( k ) * d_ln_pt + log( _pt_min ) );
    for ( int i = 0;i < eta_bins;i++ )
    {
      const double delta_eta = 2.0 * ( rapidityRanges[i].yright - rapidityRanges[i].yleft );
      
      double nInBin = double( ptBinsNmb[i][k] ) / theConfig->getTestparticles() / dpt / delta_eta;
      
      print_yield.width( 15 );
      print_yield << nInBin;
    }
    print_yield << endl;
  }

  // print yield linear bins
  print_yield_lin << "# " << type << " yield distribution" << endl;
  print_yield_lin << "# t = " << _outputTime << endl;
  print_yield_lin << "#pt\t#yield for different rapidity bins" << endl;

  for ( int k = 0;k < n_g_lin;k++ )
  {
    pt_out = double( k ) * d_pt + pt_min_lin + d_pt/2.0;
    print_yield_lin << pt_out;
    for ( int i = 0;i < eta_bins;i++ )
    {
      const double delta_eta = 2.0 * ( rapidityRanges[i].yright - rapidityRanges[i].yleft );
      
      double nInBin = double( ptBinsNmbLin[i][k] ) / theConfig->getTestparticles() / d_pt / delta_eta;
      
      print_yield_lin << "\t";
      print_yield_lin << nInBin;
    }
    print_yield_lin << endl;
  }  
  
  // print yield for RAA for different angles with respect to the reaction plane
  print_pt_angleDependence << "# " << type << " yield distribution for different angles (alpha) with respect to the reaction plane" << endl;
  print_pt_angleDependence << "# t = " << _outputTime << endl;
  print_pt_angleDependence << "#";
  print_pt_angleDependence.width( 14 );
  print_pt_angleDependence << "pt       yield for different rapidity bins" << endl;
  for ( int j = 0; j < nAlphaBins; j++ )
  {
    print_pt_angleDependence << "#alpha in [ " << j * deltaAlpha << "°, " << (j+1)*deltaAlpha << "° ] "<< endl;
    for ( int k = 0;k < n_bins + 1;k++ )
    {
      print_pt_angleDependence.width( 15 );
      pt_out = exp( double( k ) * d_ln_pt + log( _pt_min ) + d_ln_pt / 2.0 ); // +dpt/2.0 to shift value in the middle of the intervall, important for small number of bins
      print_pt_angleDependence << pt_out;
      const double dpt = exp( double( k ) * d_ln_pt + log( _pt_min ) + d_ln_pt ) - exp( double( k ) * d_ln_pt + log( _pt_min ) );
      for ( int i = 0;i < eta_bins;i++ )
      {
        const double delta_eta = 2.0 * ( rapidityRanges[i].yright - rapidityRanges[i].yleft );
        
        print_pt_angleDependence.width( 15 );
        print_pt_angleDependence << double( ptBinsAngleDep[i][j][k] ) / theConfig->getTestparticles() / dpt / delta_eta;
      }
      print_pt_angleDependence << endl;
    }    
    print_pt_angleDependence << endl;
    print_pt_angleDependence << endl;
  }
}






void analysis::finalOutput()
{
  //--------------------//
  if( studyV2 )
    v2_output();
  //--------------------//
  
  //--------------------//
  //write hydro analysis output
  if( studyHydro )
  {
    hydroDistributionOutput();
    hydroDistributionArrowOutput(); 
    hydroDistributionMidRapOutput(); 
    hydroDistributionMidRapArrowOutput();   
    particleListOutput();
    hydroMidRapidityObservablesOutput();
  }
  //--------------------//  
  
  if( studyJetTracking )
    jetTrackerOutput();

  if ( studyParticleSpectra || studyDetailedParticleOutputAllSteps )
  {
    bool printIntemediate = studyParticleSpectraAllSteps;
    printPtSpectra( gluon, printIntemediate );
    printPtSpectra( up, printIntemediate );
    printPtSpectra( down, printIntemediate );
    printPtSpectra( strange, printIntemediate );
    printPtSpectra( anti_up, printIntemediate );
    printPtSpectra( anti_down, printIntemediate );
    printPtSpectra( anti_strange, printIntemediate );
    printPtSpectra( light_quark, printIntemediate );
    printPtSpectra( allFlavors, printIntemediate );
    
    printSoftPtSpectra( gluon, printIntemediate );
    printSoftPtSpectra( up, printIntemediate );
    printSoftPtSpectra( down, printIntemediate );
    printSoftPtSpectra( strange, printIntemediate );
    printSoftPtSpectra( anti_up, printIntemediate );
    printSoftPtSpectra( anti_down, printIntemediate );
    printSoftPtSpectra( anti_strange, printIntemediate );
    printSoftPtSpectra( light_quark, printIntemediate );
    printSoftPtSpectra( allFlavors, printIntemediate );
  }
  
  //---- get time and calculate real time needed by simulation ----
  time_t end;
  time( &end );
  int secs = ( int )difftime( end, start );
  int hours = secs / 3600, minutes = ( secs % 3600 ) / 60, seconds = secs - 3600 * hours - 60 * minutes;
  cout << hours << "h" << minutes << "m" << seconds << "s" <<  endl;
  //---------------------------------------------------------------


  //---- write output file for pt-spectra -------------------------
  string filename3, filename_charmStep, filename_bottomStep;
  if( studyBoostDistribution )
    filename3 = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + "_boostDist";
  else
    filename3 = "/dev/null";
  
  filename_charmStep = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + "_charmNmbStep";
  filename_bottomStep = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + "_bottomNmbStep";
  
  if( Particle::N_heavy_flavor < 2 || !studyHeavyQuarkProduction  )
  {
    filename_bottomStep = "/dev/null";
    if( Particle::N_heavy_flavor < 1 || !studyHeavyQuarkProduction  )
      filename_charmStep = "/dev/null";
  }

  //---- write output file for distributions of boost velocities ------
  fstream file3( filename3.c_str(), ios::out | ios::trunc );
  double num32 = 0, num23 = 0;
  for ( int i = 0; i < numberBinsBoost; i++ )
  {
    num32 += boostDistribution32[i];
    num23 += boostDistribution23[i];
  }

  file3 << "#beta      3->2 (abs.)     3->2 (rel.)    2->3 (abs.)    2->3 (rel.)" << endl;
  file3 << "#Number of 3->2 events: " << num32 << endl;
  file3 << "#Number of 2->3 events: " << num23 << endl;
  file3 << "#" << endl;

  for ( int i = 0; i < numberBinsBoost; i++ )
  {
    file3 << ( i + 0.5 )*( 1.0 / numberBinsBoost ) << sep << boostDistribution32[i] << sep << boostDistribution32[i] / num32
    << sep << boostDistribution23[i] << sep << boostDistribution23[i] / num23 << endl;
  }
  file3.close();
  
  
  //---- write output file for charm number -------------------------
  fstream file_charm( filename_charmStep.c_str(), ios::out | ios::trunc );

  //---- print header if file is empty ----
  file_charm.seekp( 0, ios::end );
  long size = file_charm.tellp();
  file_charm.seekp( 0, ios::beg );
  if ( size == 0 )
    printHeader( file_charm, CharmNumber, end );
  //---------------------------------------
  
  int test_particles = theConfig->getTestparticles();
  for ( int j = 0; j <= nTimeSteps; j++ )
  {
    double time;
    
    if ( j == 0 )
      time = 0.0;
    else
      time = tstep[j-1];
    
    if( time <= theConfig->getRuntime() )
      file_charm << time << sep << 
        double( charmNmb[j] ) / 2.0 << sep << 
        double( charmNmb[j] ) / 2.0 / test_particles << sep <<
        // dN / dY
        double( charmNmb_eta20[j] ) / 2.0 / test_particles / ( 2.0*2.0 ) << sep << 
        double( charmNmb_eta15[j] ) / 2.0 / test_particles / ( 2.0*1.5 ) << sep <<
        double( charmNmb_eta10[j] ) / 2.0 / test_particles / ( 2.0*1.0 ) << sep << 
        double( charmNmb_eta075[j] ) / 2.0 / test_particles / ( 2.0*0.75 ) << sep <<
        double( charmNmb_eta05[j] ) / 2.0 / test_particles / ( 2.0*0.5 ) << sep << 
        double( charmNmb_eta035[j] ) / 2.0 / test_particles / ( 2.0*0.35 ) << sep << 
        double( charmNmb_strap05[j] ) / 2.0 / test_particles / ( 2.0*0.5 )  << endl;
    
  }
  file_charm.close();
  //------------------------------------------------------------
  
   //---- write output file for bottom number -------------------------
  fstream file_bottom( filename_bottomStep.c_str(), ios::out | ios::trunc );

  //---- print header if file is empty ----
  file_bottom.seekp( 0, ios::end );
  size = file_bottom.tellp();
  file_bottom.seekp( 0, ios::beg );
  if ( size == 0 )
    printHeader( file_bottom, BottomNumber, end );
  //---------------------------------------
  
  for ( int j = 0; j <= nTimeSteps; j++ )
  {
    if ( j == 0 )
      file_bottom << "0" << sep << double( bottomNmb[j] ) / 2.0 << sep << double( bottomNmb[j] ) / 2.0 / test_particles << sep <<
      // dN / dY
      double( bottomNmb_eta20[j] ) / 2.0 / test_particles / ( 2.0*2.0 ) << sep << double( bottomNmb_eta15[j] ) / 2.0 / test_particles / ( 2.0*1.5 ) << sep <<
      double( bottomNmb_eta10[j] ) / 2.0 / test_particles / ( 2.0*1.0 ) << sep << double( bottomNmb_eta075[j] ) / 2.0 / test_particles / ( 2.0*0.75 ) << sep <<
      double( bottomNmb_eta05[j] ) / 2.0 / test_particles / ( 2.0*0.5 ) << sep << double( bottomNmb_eta035[j] ) / 2.0 / test_particles / ( 2.0*0.35 ) << endl;
    else if ( tstep[j-1] <= theConfig->getRuntime() )
      file_bottom << tstep[j-1] << sep << double( bottomNmb[j] ) / 2.0 << sep << double( bottomNmb[j] ) / 2.0 / test_particles << sep <<
      // dN / dY
      double( bottomNmb_eta20[j] ) / 2.0 / test_particles / ( 2.0*2.0 ) << sep << double( bottomNmb_eta15[j] ) / 2.0 / test_particles / ( 2.0*1.5 ) << sep <<
      double( bottomNmb_eta10[j] ) / 2.0 / test_particles / ( 2.0*1.0 ) << sep << double( bottomNmb_eta075[j] ) / 2.0 / test_particles / ( 2.0*0.75 ) << sep <<
      double( bottomNmb_eta05[j] ) / 2.0 / test_particles / ( 2.0*0.5 ) << sep << double( bottomNmb_eta035[j] ) / 2.0 / test_particles / ( 2.0*0.35 ) << endl;
  }
  file_bottom.close();
  //------------------------------------------------------------
}



void analysis::printPtSpectra( const FLAVOR_TYPE _flavTypeToComputeFor, const bool printIntemediate )
{
  time_t end;
  time( &end );
  
  vector<double> * _ptBinsAll;
  vector<double> * _ptBinsDY5;
  vector<double> * _ptBinsDY4;
  vector<double> * _ptBinsDY3;
  vector<double> * _ptBinsDY2;
  vector<double> * _ptBinsDY16;
  vector<double> * _ptBinsDY1;
  
  if ( _flavTypeToComputeFor == gluon )
  {
    _ptBinsAll = ptBinsAll_gluons;
    _ptBinsDY5 = ptBinsDY5_gluons;
    _ptBinsDY4 = ptBinsDY4_gluons;
    _ptBinsDY3 = ptBinsDY3_gluons;
    _ptBinsDY2 = ptBinsDY2_gluons;
    _ptBinsDY16 = ptBinsDY16_gluons;
    _ptBinsDY1 = ptBinsDY1_gluons;
  }
  else if ( _flavTypeToComputeFor == light_quark )
  {
    _ptBinsAll = ptBinsAll_quarks;
    _ptBinsDY5 = ptBinsDY5_quarks;
    _ptBinsDY4 = ptBinsDY4_quarks;
    _ptBinsDY3 = ptBinsDY3_quarks;
    _ptBinsDY2 = ptBinsDY2_quarks;
    _ptBinsDY16 = ptBinsDY16_quarks;
    _ptBinsDY1 = ptBinsDY1_quarks;
  }
  else if ( _flavTypeToComputeFor == up )
  {
    _ptBinsAll = ptBinsAll_ups;
    _ptBinsDY5 = ptBinsDY5_ups;
    _ptBinsDY4 = ptBinsDY4_ups;
    _ptBinsDY3 = ptBinsDY3_ups;
    _ptBinsDY2 = ptBinsDY2_ups;
    _ptBinsDY16 = ptBinsDY16_ups;
    _ptBinsDY1 = ptBinsDY1_ups;
  }
  else if ( _flavTypeToComputeFor == down )
  {
    _ptBinsAll = ptBinsAll_downs;
    _ptBinsDY5 = ptBinsDY5_downs;
    _ptBinsDY4 = ptBinsDY4_downs;
    _ptBinsDY3 = ptBinsDY3_downs;
    _ptBinsDY2 = ptBinsDY2_downs;
    _ptBinsDY16 = ptBinsDY16_downs;
    _ptBinsDY1 = ptBinsDY1_downs;
  }
  else if ( _flavTypeToComputeFor == strange )
  {
    _ptBinsAll = ptBinsAll_stranges;
    _ptBinsDY5 = ptBinsDY5_stranges;
    _ptBinsDY4 = ptBinsDY4_stranges;
    _ptBinsDY3 = ptBinsDY3_stranges;
    _ptBinsDY2 = ptBinsDY2_stranges;
    _ptBinsDY16 = ptBinsDY16_stranges;
    _ptBinsDY1 = ptBinsDY1_stranges;
  }
  else if ( _flavTypeToComputeFor == anti_up )
  {
    _ptBinsAll = ptBinsAll_anti_ups;
    _ptBinsDY5 = ptBinsDY5_anti_ups;
    _ptBinsDY4 = ptBinsDY4_anti_ups;
    _ptBinsDY3 = ptBinsDY3_anti_ups;
    _ptBinsDY2 = ptBinsDY2_anti_ups;
    _ptBinsDY16 = ptBinsDY16_anti_ups;
    _ptBinsDY1 = ptBinsDY1_anti_ups;
  }
  else if ( _flavTypeToComputeFor == anti_down )
  {
    _ptBinsAll = ptBinsAll_anti_downs;
    _ptBinsDY5 = ptBinsDY5_anti_downs;
    _ptBinsDY4 = ptBinsDY4_anti_downs;
    _ptBinsDY3 = ptBinsDY3_anti_downs;
    _ptBinsDY2 = ptBinsDY2_anti_downs;
    _ptBinsDY16 = ptBinsDY16_anti_downs;
    _ptBinsDY1 = ptBinsDY1_anti_downs;
  }
  else if ( _flavTypeToComputeFor == anti_strange )
  {
    _ptBinsAll = ptBinsAll_anti_stranges;
    _ptBinsDY5 = ptBinsDY5_anti_stranges;
    _ptBinsDY4 = ptBinsDY4_anti_stranges;
    _ptBinsDY3 = ptBinsDY3_anti_stranges;
    _ptBinsDY2 = ptBinsDY2_anti_stranges;
    _ptBinsDY16 = ptBinsDY16_anti_stranges;
    _ptBinsDY1 = ptBinsDY1_anti_stranges;
  }
  else if ( _flavTypeToComputeFor == allFlavors )
  {
    _ptBinsAll = ptBinsAll_all;
    _ptBinsDY5 = ptBinsDY5_all;
    _ptBinsDY4 = ptBinsDY4_all;
    _ptBinsDY3 = ptBinsDY3_all;
    _ptBinsDY2 = ptBinsDY2_all;
    _ptBinsDY16 = ptBinsDY16_all;
    _ptBinsDY1 = ptBinsDY1_all;
  }
  else
  {
    string errMsg = "error in ptDistribution, flavor not specified";
    throw eAnalysis_error( errMsg );
  }
  
  string type;
  if ( _flavTypeToComputeFor == gluon )
  {
    type = "gluon";
  }
  else if ( _flavTypeToComputeFor == light_quark || _flavTypeToComputeFor == anti_light_quark )
  {
    type = "quark";
  }
  if ( _flavTypeToComputeFor == up )
  {
    type = "up";
  }
  if ( _flavTypeToComputeFor == down )
  {
    type = "down";
  }
  if ( _flavTypeToComputeFor == strange )
  {
    type = "strange";
  }
  if ( _flavTypeToComputeFor == anti_up )
  {
    type = "anti_up";
  }
  if ( _flavTypeToComputeFor == anti_down )
  {
    type = "anti_down";
  }
  if ( _flavTypeToComputeFor == anti_strange )
  {
    type = "anti_strange";
  }
  else if ( _flavTypeToComputeFor == allFlavors )
  {
    type = "allFlavors";
  }
  
  
  string filename = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + "_" + type + "_spectra.f2";
  fstream file( filename.c_str(), ios::out | ios::trunc );
  
  //---- print header if file is empty ----
  file.seekp( 0, ios::end );
  long size = file.tellp();
  file.seekp( 0, ios::beg );
  if ( size == 0 )
  {
    printHeader( file, ptSpectrum, end );
  }
  //---------------------------------------
  
  //---------------------- y in [-0.5,0.5] ---------------------
  file << "#y in [-0.5,0.5]" << endl;
  for ( int i = 0; i < numberBinsPT; i++ )
  {
    file << ptBinLabels[i] << sep;
    file << _ptBinsDY1[0][i] << sep;
    if ( printIntemediate )
    {
      for ( int j = 1; j < nTimeSteps && tstep[j-1] <= theConfig->getRuntime(); j++ )
      {
        file << _ptBinsDY1[j][i] << sep;
      }
    }
    file << _ptBinsDY1[nTimeSteps][i] << sep;
    file << endl;
  }
  //------------------------------------------------------------
  file << endl << endl;
  
  //---------------------- y in [-0.8,0.8] ---------------------
  file << "#y in [-0.8,0.8]" << endl;
  for ( int i = 0; i < numberBinsPT; i++ )
  {
    file << ptBinLabels[i] << sep;
    file << _ptBinsDY16[0][i] << sep;
    if ( printIntemediate )
    {
      for ( int j = 1; j < nTimeSteps && tstep[j-1] <= theConfig->getRuntime(); j++ )
      {
        file << _ptBinsDY16[j][i] << sep;
      }
    }
    file << _ptBinsDY16[nTimeSteps][i] << sep;
    file << endl;
  }
  //------------------------------------------------------------
  file << endl << endl;
  
  //---------------------- y in [-1.0,1.0] ---------------------
  file << "#y in [-1.0,1.0]" << endl;
  for ( int i = 0;i < numberBinsPT;i++ )
  {
    file << ptBinLabels[i] << sep;
    file << _ptBinsDY2[0][i] << sep;
    if ( printIntemediate )
    {
      for ( int j = 1; j < nTimeSteps && tstep[j-1] <= theConfig->getRuntime(); j++ )
      {
        file << _ptBinsDY2[j][i] << sep;
      }
    }
    file << _ptBinsDY2[nTimeSteps][i] << sep;
    file << endl;
  }
  //------------------------------------------------------------
  file << endl << endl;
  
  
  //---------------------- y in [-1.5,1.5] ---------------------
  file << "#y in [-1.5,1.5]" << endl;
  for ( int i = 0; i < numberBinsPT; i++ )
  {
    file << ptBinLabels[i] << sep;
    file << _ptBinsDY3[0][i] << sep;
    if ( printIntemediate )
    {
      for ( int j = 1; j < nTimeSteps && tstep[j-1] <= theConfig->getRuntime(); j++ )
      {
        file << _ptBinsDY3[j][i] << sep;
      }
    }
    file << _ptBinsDY3[nTimeSteps][i] << sep;
    file << endl;
  }
  //------------------------------------------------------------
  file << endl << endl;
  
  //---------------------- y in [-2.0,2.0] ---------------------
  file << "#y in [-2.0,2.0]" << endl;
  for ( int i = 0; i < numberBinsPT; i++ )
  {
    file << ptBinLabels[i] << sep;
    file << _ptBinsDY4[0][i] << sep;
    if ( printIntemediate )
    {
      for ( int j = 1; j < nTimeSteps && tstep[j-1] <= theConfig->getRuntime(); j++ )
      {
        file << _ptBinsDY4[j][i] << sep;
      }
    }
    file << _ptBinsDY4[nTimeSteps][i] << sep;
    file << endl;
  }
  //------------------------------------------------------------
  file << endl << endl;
  
  //---------------------- y in [-2.5,2.5] ---------------------
  file << "#y in [-2.5,2.5]" << endl;
  for ( int i = 0;i < numberBinsPT;i++ )
  {
    file << ptBinLabels[i] << sep;
    file << _ptBinsDY5[0][i] << sep;
    if ( printIntemediate )
    {
      for ( int j = 1; j < nTimeSteps && tstep[j-1] <= theConfig->getRuntime(); j++ )
      {
        file << _ptBinsDY5[j][i] << sep;
      }
    }
    file << _ptBinsDY5[nTimeSteps][i] << sep;
    file << endl;
  }
  //------------------------------------------------------------
  file << endl << endl;
  
  //---------------------- y in [-inf,inf] ---------------------
  file << "#y in [-inf,inf]" << endl;
  for ( int i = 0;i < numberBinsPT;i++ )
  {
    file << ptBinLabels[i] << sep;
    file << _ptBinsAll[0][i] << sep;
    if ( printIntemediate )
    {
      for ( int j = 1; j < nTimeSteps && tstep[j-1] <= theConfig->getRuntime(); j++ )
      {
        file << _ptBinsAll[j][i] << sep;
      }
    }
    file << _ptBinsAll[nTimeSteps][i] << sep;
    file << endl;
  }
  //------------------------------------------------------------
  
  file.close();
  //---------------------------------------------------------------
  
}



void analysis::printSoftPtSpectra( const FLAVOR_TYPE _flavTypeToComputeFor, const bool printIntemediate )
{
  time_t end;
  time( &end );
  
  vector<double> * _ptBinsSoftAll;
  
  if ( _flavTypeToComputeFor == gluon )
  {
    _ptBinsSoftAll = ptBinsSoftAll_gluons;
  }
  else if ( _flavTypeToComputeFor == light_quark )
  {
    _ptBinsSoftAll = ptBinsSoftAll_quarks;
  }
  else if ( _flavTypeToComputeFor == up )
  {
    _ptBinsSoftAll = ptBinsSoftAll_ups;
  }
  else if ( _flavTypeToComputeFor == down )
  {
    _ptBinsSoftAll = ptBinsSoftAll_downs;
  }
  else if ( _flavTypeToComputeFor == strange )
  {
    _ptBinsSoftAll = ptBinsSoftAll_stranges;
  }
  else if ( _flavTypeToComputeFor == anti_up )
  {
    _ptBinsSoftAll = ptBinsSoftAll_anti_ups;
  }
  else if ( _flavTypeToComputeFor == anti_down )
  {
    _ptBinsSoftAll = ptBinsSoftAll_anti_downs;
  }
  else if ( _flavTypeToComputeFor == anti_strange )
  {
    _ptBinsSoftAll = ptBinsSoftAll_anti_stranges;
  }
  else if ( _flavTypeToComputeFor == allFlavors )
  {
    _ptBinsSoftAll = ptBinsSoftAll_all;
  }
  else
  {
    string errMsg = "error in printSoftPtSpectra, flavor not specified";
    throw eAnalysis_error( errMsg );
  }
  
  string type;
  if ( _flavTypeToComputeFor == gluon )
  {
    type = "gluon";
  }
  else if ( _flavTypeToComputeFor == light_quark || _flavTypeToComputeFor == anti_light_quark )
  {
    type = "quark";
  }
  if ( _flavTypeToComputeFor == up )
  {
    type = "up";
  }
  if ( _flavTypeToComputeFor == down )
  {
    type = "down";
  }
  if ( _flavTypeToComputeFor == strange )
  {
    type = "strange";
  }
  if ( _flavTypeToComputeFor == anti_up )
  {
    type = "anti_up";
  }
  if ( _flavTypeToComputeFor == anti_down )
  {
    type = "anti_down";
  }
  if ( _flavTypeToComputeFor == anti_strange )
  {
    type = "anti_strange";
  }
  else if ( _flavTypeToComputeFor == allFlavors )
  {
    type = "allFlavors";
  }
  
  string filename = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + "_" + type + "_spectra.f3";
  fstream file( filename.c_str(), ios::out | ios::trunc );
  
  //---- print header if file is empty ----
  file.seekp( 0, ios::end );
  long size = file.tellp();
  file.seekp( 0, ios::beg );
  if ( size == 0 )
  {
    printHeader( file, ptSpectrum, end );
  }
  //---------------------------------------
  
  
  //---------------------- y in [-inf,inf] ---------------------
  file << "#y in [-inf,inf]" << endl;
  for ( int i = 0;i < numberBinsSoftPT;i++ )
  {
    file << ptSoftBinLabels[i] << sep;
    file << _ptBinsSoftAll[0][i] << sep;
    if ( printIntemediate )
    {
      for ( int j = 1; j < nTimeSteps && tstep[j-1] <= theConfig->getRuntime(); j++ )
      {
        file << _ptBinsSoftAll[j][i] << sep;
      }
    }
    file << _ptBinsSoftAll[nTimeSteps][i] << sep;
    file << endl;
  }
  //------------------------------------------------------------
  
  file.close();
  //---------------------------------------------------------------   
}


void analysis::jetTrackerOutput()
{
  //---- get time and calculate real time needed by simulation ----
  time_t end;
  time( &end );
  int secs = ( int )difftime( end, start );
  int hours = secs / 3600, minutes = ( secs % 3600 ) / 60, seconds = secs - 3600 * hours - 60 * minutes;
  cout << hours << "h" << minutes << "m" << seconds << "s" <<  endl;
  //---------------------------------------------------------------


  //---- write output file for pt-spectra -------------------------
  string filename = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + ".f4";
  fstream file( filename.c_str(), ios::out | ios::trunc );

  //---- print header if file is empty ----
  file.seekp( 0, ios::end );
  long size = file.tellp();
  file.seekp( 0, ios::beg );
  if ( size == 0 )
    printHeader( file, jets, end );
  //---------------------------------------

  for ( unsigned int jet = 0; jet < jetTracker.size(); jet++ )
  {
    for ( unsigned int event = 0; event < jetTracker[jet].size(); event++ )
    {
      file << jetTracker[jet][event].jet_ID_in << sep << jetTracker[jet][event].jet_ID_out << sep
      << jetTracker[jet][event].coll_type << sep;
      // file << jetTracker[jet][event].P_proj_in << sep;
      // file << jetTracker[jet][event].P_proj_out << sep;
      // file << jetTracker[jet][event].R_proj << sep;
      for ( int i = 0; i < 4; i++ )
      {
        file << jetTracker[jet][event].P_proj_in(i) << sep;
      }
      for ( int i = 0; i < 4; i++ )
      {
        file << jetTracker[jet][event].P_proj_out(i) << sep;
      }
      for ( int i = 0; i < 4; i++ )
      {
        file << jetTracker[jet][event].R_proj(i) << sep;
      }
      file << jetTracker[jet][event].xSection << sep << jetTracker[jet][event].lambda << endl;
    }
    file << endl << endl;
  }


}



void analysis::ptDistribution( const FLAVOR_TYPE _flavTypeToComputeFor, vector<Particle>& _particles, const int n_particles, const int step )
{
  double pt, y, v2;
  vector<double> * _ptBinsAll;
  vector<double> * _ptBinsDY5;
  vector<double> * _ptBinsDY4;
  vector<double> * _ptBinsDY3;
  vector<double> * _ptBinsDY2;
  vector<double> * _ptBinsDY16;
  vector<double> * _ptBinsDY1;
  
  if ( _flavTypeToComputeFor == gluon )
  {
    _ptBinsAll = ptBinsAll_gluons;
    _ptBinsDY5 = ptBinsDY5_gluons;
    _ptBinsDY4 = ptBinsDY4_gluons;
    _ptBinsDY3 = ptBinsDY3_gluons;
    _ptBinsDY2 = ptBinsDY2_gluons;
    _ptBinsDY16 = ptBinsDY16_gluons;
    _ptBinsDY1 = ptBinsDY1_gluons;
  }
  else if ( _flavTypeToComputeFor == light_quark )
  {
    _ptBinsAll = ptBinsAll_quarks;
    _ptBinsDY5 = ptBinsDY5_quarks;
    _ptBinsDY4 = ptBinsDY4_quarks;
    _ptBinsDY3 = ptBinsDY3_quarks;
    _ptBinsDY2 = ptBinsDY2_quarks;
    _ptBinsDY16 = ptBinsDY16_quarks;
    _ptBinsDY1 = ptBinsDY1_quarks;
  }
  else if ( _flavTypeToComputeFor == up )
  {
    _ptBinsAll = ptBinsAll_ups;
    _ptBinsDY5 = ptBinsDY5_ups;
    _ptBinsDY4 = ptBinsDY4_ups;
    _ptBinsDY3 = ptBinsDY3_ups;
    _ptBinsDY2 = ptBinsDY2_ups;
    _ptBinsDY16 = ptBinsDY16_ups;
    _ptBinsDY1 = ptBinsDY1_ups;
  }
  else if ( _flavTypeToComputeFor == down )
  {
    _ptBinsAll = ptBinsAll_downs;
    _ptBinsDY5 = ptBinsDY5_downs;
    _ptBinsDY4 = ptBinsDY4_downs;
    _ptBinsDY3 = ptBinsDY3_downs;
    _ptBinsDY2 = ptBinsDY2_downs;
    _ptBinsDY16 = ptBinsDY16_downs;
    _ptBinsDY1 = ptBinsDY1_downs;
  }
  else if ( _flavTypeToComputeFor == strange )
  {
    _ptBinsAll = ptBinsAll_stranges;
    _ptBinsDY5 = ptBinsDY5_stranges;
    _ptBinsDY4 = ptBinsDY4_stranges;
    _ptBinsDY3 = ptBinsDY3_stranges;
    _ptBinsDY2 = ptBinsDY2_stranges;
    _ptBinsDY16 = ptBinsDY16_stranges;
    _ptBinsDY1 = ptBinsDY1_stranges;
  }
  else if ( _flavTypeToComputeFor == anti_up )
  {
    _ptBinsAll = ptBinsAll_anti_ups;
    _ptBinsDY5 = ptBinsDY5_anti_ups;
    _ptBinsDY4 = ptBinsDY4_anti_ups;
    _ptBinsDY3 = ptBinsDY3_anti_ups;
    _ptBinsDY2 = ptBinsDY2_anti_ups;
    _ptBinsDY16 = ptBinsDY16_anti_ups;
    _ptBinsDY1 = ptBinsDY1_anti_ups;
  }
  else if ( _flavTypeToComputeFor == anti_down )
  {
    _ptBinsAll = ptBinsAll_anti_downs;
    _ptBinsDY5 = ptBinsDY5_anti_downs;
    _ptBinsDY4 = ptBinsDY4_anti_downs;
    _ptBinsDY3 = ptBinsDY3_anti_downs;
    _ptBinsDY2 = ptBinsDY2_anti_downs;
    _ptBinsDY16 = ptBinsDY16_anti_downs;
    _ptBinsDY1 = ptBinsDY1_anti_downs;
  }
  else if ( _flavTypeToComputeFor == anti_strange )
  {
    _ptBinsAll = ptBinsAll_anti_stranges;
    _ptBinsDY5 = ptBinsDY5_anti_stranges;
    _ptBinsDY4 = ptBinsDY4_anti_stranges;
    _ptBinsDY3 = ptBinsDY3_anti_stranges;
    _ptBinsDY2 = ptBinsDY2_anti_stranges;
    _ptBinsDY16 = ptBinsDY16_anti_stranges;
    _ptBinsDY1 = ptBinsDY1_anti_stranges;
  }
  else if ( _flavTypeToComputeFor == allFlavors )
  {
    _ptBinsAll = ptBinsAll_all;
    _ptBinsDY5 = ptBinsDY5_all;
    _ptBinsDY4 = ptBinsDY4_all;
    _ptBinsDY3 = ptBinsDY3_all;
    _ptBinsDY2 = ptBinsDY2_all;
    _ptBinsDY16 = ptBinsDY16_all;
    _ptBinsDY1 = ptBinsDY1_all;
  }
  else
  {
    string errMsg = "error in ptDistribution, flavor not specified";
    throw eAnalysis_error( errMsg );
  }
  
  
  FLAVOR_TYPE genFlavor;
  for ( int j = 0; j < n_particles; j++ )
  {
    pt = _particles[j].Mom.Pt();
    y  = _particles[j].Mom.Rapidity();
    v2 = ( pow( particles[j].Mom.Px(), 2.0 ) - pow( particles[j].Mom.Py(), 2.0 ) ) / 
      pow( pt, 2.0 );

        
    genFlavor = Particle::mapToGenericFlavorType( _particles[j].FLAVOR );
    if ( _flavTypeToComputeFor == allFlavors || _particles[j].FLAVOR == _flavTypeToComputeFor || 
      ( _flavTypeToComputeFor == light_quark && ( genFlavor == light_quark || genFlavor == anti_light_quark ) ) )
    {
      //------------------------ y in [-inf,inf] -----------------------
      if ( pt <= maxPT && pt >= minPT )
      {
        if ( pt == maxPT )
        {
          ++_ptBinsAll[step][numberBinsPT - 1];
        }
        else
        {
          if ( pt == minPT )
            ++_ptBinsAll[step][0];
          else
            ++_ptBinsAll[step][int(( pt - minPT )/binWidthPT )];
        }
        
        //----------------------------------------------------------------
        
        //------------------------ y in [-2.5,2.5] -----------------------
        if ( fabs( y ) <= 2.5 )
        {
          if ( pt == maxPT )
          {
            ++_ptBinsDY5[step][numberBinsPT - 1];
          }
          else
          {
            if ( pt == minPT )
              ++_ptBinsDY5[step][0];
            else
              ++_ptBinsDY5[step][int(( pt - minPT )/binWidthPT )];
          }
        }
        //----------------------------------------------------------------
        
        //------------------------ y in [-2.0,2.0] -----------------------
        if ( fabs( y ) <= 2.0 )
        {
          if ( pt == maxPT )
          {
            ++_ptBinsDY4[step][numberBinsPT - 1];
          }
          else
          {
            if ( pt == minPT )
              ++_ptBinsDY4[step][0];
            else
              ++_ptBinsDY4[step][int(( pt - minPT )/binWidthPT )];
          }
        }
        //----------------------------------------------------------------
        
        //------------------------ y in [-1.5,1.5] -----------------------
        if ( fabs( y ) <= 1.5 )
        {
          if ( pt == maxPT )
          {
            ++_ptBinsDY3[step][numberBinsPT - 1];
          }
          else
          {
            if ( pt == minPT )
              ++_ptBinsDY3[step][0];
            else
              ++_ptBinsDY3[step][int(( pt - minPT )/binWidthPT )];
          }
        }
        //----------------------------------------------------------------
        
        //------------------------ y in [-1.0,1.0] -----------------------
        if ( fabs( y ) <= 1.0 )
        {
          if ( pt == maxPT )
          {
            ++_ptBinsDY2[step][numberBinsPT - 1];
          }
          else
          {
            if ( pt == minPT )
              ++_ptBinsDY2[step][0];
            else
              ++_ptBinsDY2[step][int(( pt - minPT )/binWidthPT )];
          }
        }
        //----------------------------------------------------------------
        
        
        //------------------------ y in [-0.8,0.8] -----------------------
        if ( fabs( y ) <= 0.8 )
        {
          if ( pt == maxPT )
          {
            ++_ptBinsDY16[step][numberBinsPT - 1];
          }
          else
          {
            if ( pt == minPT )
              ++_ptBinsDY16[step][0];
            else
              ++_ptBinsDY16[step][int(( pt - minPT )/binWidthPT )];
          }
        }
        //----------------------------------------------------------------
        
        
        //------------------------ y in [-0.5,0.5] -----------------------
        if ( fabs( y ) <= 0.5 )
        {
          if ( pt == maxPT )
          {
            ++_ptBinsDY1[step][numberBinsPT - 1];
          }
          else
          {
            if ( pt == minPT )
              ++_ptBinsDY1[step][0];
            else
              ++_ptBinsDY1[step][int(( pt - minPT )/binWidthPT )];
          }
        }
      }
      //----------------------------------------------------------------
    }
    
    
    //---------------- some v2-binning -----------------------------
    int pt_index = static_cast<int>(( pt - min_v2_pt ) / binWidth_v2_pt );
    int y_bin_max_index = 0;
    
    // rapidity bins are sorted from large bins (index 0 = infinity) to smaller bins
    while ( y_bin_max_index + 1 < nRapidityBins && rapidityBins[y_bin_max_index+1] >= fabs(y) )
    {
      ++y_bin_max_index;
    }
    
    for( int y_index = 0; y_index <= y_bin_max_index; y_index++ )
    {
      if ( pt < max_v2_pt && pt > min_v2_pt )
      {
        differential_v2_values[step][y_index][pt_index] += v2;
        differential_v2_counts[step][y_index][pt_index] += 1;
      }

      if ( _flavTypeToComputeFor == allFlavors )
      {
        mean_PT[step][y_index] += pt;
        integrated_v2[step][y_index] += v2;
        integrated_v2_counts[step][y_index] += 1;
      }
    }
    //---------------- some v2-binning -----------------------------
  }
}









void analysis::ptSoftDistribution( const FLAVOR_TYPE _flavTypeToComputeFor, vector<Particle>& _particles, const int n_particles, const int step )
{
  double pt;
  
  vector<double>* _ptBinsSoftAll;
  
  if ( _flavTypeToComputeFor == gluon )
  {
    _ptBinsSoftAll = ptBinsSoftAll_gluons;
  }
  else if ( _flavTypeToComputeFor == light_quark )
  {
    _ptBinsSoftAll = ptBinsSoftAll_quarks;
  }
  else if ( _flavTypeToComputeFor == up )
  {
    _ptBinsSoftAll = ptBinsSoftAll_ups;
  }
  else if ( _flavTypeToComputeFor == down )
  {
    _ptBinsSoftAll = ptBinsSoftAll_downs;
  }
  else if ( _flavTypeToComputeFor == strange )
  {
    _ptBinsSoftAll = ptBinsSoftAll_stranges;
  }
  else if ( _flavTypeToComputeFor == anti_up )
  {
    _ptBinsSoftAll = ptBinsSoftAll_anti_ups;
  }
  else if ( _flavTypeToComputeFor == anti_down )
  {
    _ptBinsSoftAll = ptBinsSoftAll_anti_downs;
  }
  else if ( _flavTypeToComputeFor == anti_strange )
  {
    _ptBinsSoftAll = ptBinsSoftAll_anti_stranges;
  }
  else if ( _flavTypeToComputeFor == allFlavors )
  {
    _ptBinsSoftAll = ptBinsSoftAll_all;
  }
  else
  {
    string errMsg = "error in ptSoftDistribution, flavor not specified";
    throw eAnalysis_error( errMsg );
  }
  
  FLAVOR_TYPE genFlavor;
  for ( int j = 0; j < n_particles; j++ )
  {
    pt = _particles[j].Mom.Pt();
    
    genFlavor = Particle::mapToGenericFlavorType( _particles[j].FLAVOR );
    if ( _flavTypeToComputeFor == allFlavors || _particles[j].FLAVOR == _flavTypeToComputeFor || 
      ( _flavTypeToComputeFor == light_quark && ( genFlavor == light_quark || genFlavor == anti_light_quark ) ) )
    {
      //------------------------ y in [-inf,inf] -----------------------
      if ( pt <= maxPTSoft && pt >= 0 )
      {
        if ( pt == maxPTSoft )
        {
          ++_ptBinsSoftAll[step][numberBinsSoftPT - 1];
        }
        else
        {
          if ( pt == 0 )
            ++_ptBinsSoftAll[step][0];
          else
            ++_ptBinsSoftAll[step][int( pt/binWidthSoftPT )];
        }
      }
      //----------------------------------------------------------------
    }
  }
}





void analysis::charmNumber( const int step )
{
  double eta;

  for ( unsigned int i = 0; i < particles.size(); i++ )
  {
    if ( particles[i].FLAVOR == 7 || particles[i].FLAVOR == 8 )
    {
      charmNmb[step] += 1; // count charm in whole collision

      eta = particles[i].Mom.Pseudorapidity( particles[i].m );

      if ( fabs( eta ) <= 2.0 )
      {
        charmNmb_eta20[step] += 1;

        if ( fabs( eta ) <= 1.5 )
        {
          charmNmb_eta15[step] += 1;

          if ( fabs( eta ) <= 1.0 )
          {
            charmNmb_eta10[step] += 1;

            if ( fabs( eta ) <= 0.75 )
            {
              charmNmb_eta075[step] += 1;

              if ( fabs( eta ) <= 0.5 )
              {
                charmNmb_eta05[step] += 1;

                if ( fabs( eta ) <= 0.35 )
                {
                  charmNmb_eta035[step] += 1;
                }
              }
            }
          }
        }
      }
      
      
      // for Kai: Just in spacetime rapregion
      // spacetime rapidity
      if ( fabs( particles[i].Pos.Rapidity() ) <= 0.5 )
        charmNmb_strap05[step] += 1;
      
    }
  }
}


void analysis::bottomNumber( const int step )
{
  double eta;

  for ( unsigned int i = 0; i < particles.size(); i++ )
  {
    if ( particles[i].FLAVOR == 9 || particles[i].FLAVOR == 10 )
    {
      bottomNmb[step] += 1; // count bottom in whole collision

      eta = particles[i].Mom.Pseudorapidity( particles[i].m );

      if ( fabs( eta ) <= 2.0 )
      {
        bottomNmb_eta20[step] += 1;

        if ( fabs( eta ) <= 1.5 )
        {
          bottomNmb_eta15[step] += 1;

          if ( fabs( eta ) <= 1.0 )
          {
            bottomNmb_eta10[step] += 1;

            if ( fabs( eta ) <= 0.75 )
            {
              bottomNmb_eta075[step] += 1;

              if ( fabs( eta ) <= 0.5 )
              {
                bottomNmb_eta05[step] += 1;

                if ( fabs( eta ) <= 0.35 )
                {
                  bottomNmb_eta035[step] += 1;
                }
              }
            }
          }
        }
      }
    }
  }
}






void analysis::particleOutput( const int step )
{
  string name;
  stringstream ss;

  if ( step == 0 )
    name = "initial";
  else if ( step == nTimeSteps )
    name = "final";
  else
  {
    ss << step;
    name = "step" + ss.str();
  }

  //creates filename, for example: "./output/run32_step1.f1", "./output/run32_initial.f1" or the likes
  string filename = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + "_" + name + ".f1";
  fstream file( filename.c_str(), ios::out | ios::trunc );

  //---- print header if file is empty ----
  time_t end;
  time( &end );

  file.seekp( 0, ios::end );
  long size = file.tellp();
  file.seekp( 0, ios::beg );
  if ( size == 0 )
    printHeader( file, all, end );
  //---------------------------------------

  for ( unsigned int i = 0; i < particles.size(); i++ )
  {
    file << i << sep << particles[i].unique_id << sep << particles[i].cell_id << sep << particles[i].FLAVOR << sep 
         << particles[i].Pos.T() << sep 
         << particles[i].Pos.X() << sep 
         << particles[i].Pos.Y() << sep 
         << particles[i].Pos.Z() << sep 
         << particles[i].Mom.E() << sep 
         << particles[i].Mom.Px() << sep 
         << particles[i].Mom.Py() << sep
         << particles[i].Mom.Pz() << sep 
         << particles[i].md2g << sep << particles[i].md2q << endl;
  }
  file.close();
}


void analysis::ParticlePosMomOutput( const int step, const double timeshift)
{
  string name;
  stringstream ss;
  
  ss.precision(1);
  ss.setf(ios::fixed,ios::floatfield);
  
//   double pt, p, eta, phi;

  if ( step == 0 )
    name = "initial";
  else if ( step == nTimeSteps )
    name = "final";
  else
  {
    ss.precision(int(tstep[step-1]*10000.0)%1000?4:1);
    ss << tstep[step-1];
    name = "time_" + ss.str() + "fm";
  }
  
  //creates filename, for example: "./output/run32_step1.f1", "./output/run32_initial.f1" or the likes
  string filename = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + "_" + name + ".f0";
  fstream file( filename.c_str(), ios::out | ios::trunc );
  
  //---- print header if file is empty ----
  time_t end;
  time( &end );
  
  file.seekp( 0, ios::end );
  long size = file.tellp();
  file.seekp( 0, ios::beg );
  if ( size == 0 )
  {
    file << "#position and momentum of particles existing at current timestep" << endl;
    file << "#time in this output: " << tstep[step-1] << " fm/c" << endl;
    printHeader( file, pos_mom, end );
  }
  //---------------------------------------
  
  for ( unsigned int i = 0; i < particles.size(); i++ )
  {
    if ( FPT_COMP_LE( particles[i].Pos.T()-timeshift, tstep[step-1] ) || !step ) // initial files contains ALL particles, even not yet sampled, otherwise it would be empty
    {
      //         pt = sqrt(particles[i].PX*particles[i].PX + particles[i].PY*particles[i].PY);
      //         p = sqrt(particles[i].PX*particles[i].PX + particles[i].PY*particles[i].PY + particles[i].PZ*particles[i].PZ);
      //         if (particles[i].PZ == p)
      //             eta = ( particles[i].PZ * 1000 ) / fabs( particles[i].PZ ); // eta +-infinity
      //           else
      //             eta = 0.5 * log( ( p + particles[i].PZ ) / ( p - particles[i].PZ ) ); // eta
      //             if ( pt == 0.0 ) // |p_t| = 0
      //         phi = 1000.0;
      //         else
      //         {
      //           phi = acos( particles[i].PX / pt ) / 3.1415926535897 * 180.0;
      //           if (particles[i].PY < 0)
      //             phi *= (-1.0);
      //         }
      file << particles[i].Pos.X() << sep 
           << particles[i].Pos.Y() << sep 
           << particles[i].Pos.Z() << sep
           << particles[i].Mom.Px() << sep 
           << particles[i].Mom.Py() << sep
           << particles[i].Mom.Pz() << endl; 
    }
  }
  file.close();
}


void analysis::printHeader( fstream & f, const anaType mode, const time_t end )
{
  switch ( mode )
  {
  case ptSpectrum:
    f << "#pt-spectra " << endl;
    break;
  case ptSpectrumSoft:
    f << "#soft pt-spectra " << endl;
    break;
  case jets:
    f << "#jet tracking informatino " << endl;
    break;
  case all:
    f << "#information on particle locations, momenta etc. at step indicated by filename" << endl;
    f << "#initial = 0 fm/c" << endl;
    for ( int i = 0;i < nTimeSteps - 1;i++ )
      f << "#step " << ( i + 1 ) << " = " << tstep[i] << " fm/c" << endl;
    f << "#final = " << theConfig->getRuntime() << " fm/c" << endl;
    break;
  case v2output:
    f << "#differential v2" << endl;
    break;
  case meanv2output:
    f << "#<v2> and <pt>" << endl;
    break;
  case CharmNumber:
    f << "#numbers of charm quarks";
    break;
  case BottomNumber:
    f << "#numbers of bottom quarks";
    break;
  case pos_mom:
    f << "#initial = 0 fm/c" << endl;
    f << "#final = " << theConfig->getRuntime() << " fm/c" << endl;
    break;
  default: f << "#undefined ";
  }

  f << "#start: " << ctime( &start );
  f << "#end: " << ctime( &end );
  f << "#" << endl;
  f << "#simulation parameter:" << endl;
  f << "#testparticles= " << theConfig->getTestparticles() << endl;
  f << "#runtime= " << theConfig->getRuntime() << endl;
  f << "#sqrtS= " << theConfig->getSqrtS() << " GeV" << endl;
  f << "#P0= " << theConfig->getPtCutoff() << " GeV" << endl;
  f << "#b= " << theConfig->getImpactParameter() << " fm" << endl;
  f << "#(" << theConfig->getA() << "," << theConfig->getAatomic() << ") on ("
  << theConfig->getB() << "," << theConfig->getBatomic() << ")" << endl;
  f << "#seed for random generator ran2(): " << seed << endl;
  f << "#" << endl;

  stringstream ss;

  switch ( mode )
  {
  case ptSpectrum:
    f << "#numbers NOT yet corrected for width of bins and number of testparticles!" << endl;
    f << "#binWidth= " << binWidthPT << endl;
    f << "#testparticles= " << theConfig->getTestparticles() << endl;
    f << "#" << endl;
    f << "#pt [GeV]" << sep;
    for ( int j = 0; j <= nTimeSteps; j++ )
    {
      if ( j == 0 )
      {
        f << "initial" << sep;
      }
      else if ( j == nTimeSteps )
      {
        f << "final = " << theConfig->getRuntime() << " fm/c";
      }
      else if ( studyParticleSpectraAllSteps && tstep[j-1] <= theConfig->getRuntime() )
      {
        ss << tstep[j-1];
        f << ss.str() << " fm/c" << sep;
        ss.str( "" );
      }
    }
    f << endl;
    break;
  case ptSpectrumSoft:
    f << "#numbers NOT yet corrected for width of bins and number of testparticles!" << endl;
    f << "#binWidth= " << binWidthSoftPT << endl;
    f << "#testparticles= " << theConfig->getTestparticles() << endl;
    f << "#" << endl;
    f << "#pt [GeV]" << sep;
    for ( int j = 0; j <= nTimeSteps; j++ )
    {
      if ( j == 0 )
        f << "initial" << sep;
      else if ( j == nTimeSteps )
        f << "final = " << theConfig->getRuntime() << " fm/c";
      else if ( tstep[j-1] <= theConfig->getRuntime() )
      {
        ss << tstep[j-1];
        f << ss.str() << " fm/c" << sep;
        ss.str( "" );
      }
    }
    f << endl;
    break;
  case jets:
    f << "#jetID_in" << sep << "jetID_out" << sep << "coll_type" << sep
    << "P_in[0]" << sep << "P_in[1]" << sep << "P_in[2]" << sep  << "P_in[3]" << sep
    << "P_out[0]" << sep << "P_out[1]" << sep << "P_out[2]" << sep  << "P_out[3]" << sep
    << "R[0]" << sep << "R[1]" << sep << "R[2]" << sep << "R[3]" << sep
    << "xSection" << sep << "lambda" << endl;
    break;
  case all:
    f << "#ID" << sep << "uniqueID" << sep << "Cell ID" << sep << "Flavor" << sep << "t" << sep << "x" << sep << "y" << sep << "z" << sep
    << "E"  << sep << "Px"  << sep << "Py"  << sep << "Pz" << sep << "md2g/alpha_s [GeV^2]" << sep
    << "md2q/alpha_s [GeV^2]" << endl;
    break;
  case v2output:
    f << "#testparticles= " << theConfig->getTestparticles() << endl;
    f << "#" << endl;
    f << "#pt [GeV]" << sep;
    for ( int j = 0; j <= nTimeSteps; j++ )
    {
      if ( j == 0 )
        f << "initial" << sep;
      else if ( j == nTimeSteps )
        f << "final = " << theConfig->getRuntime() << " fm/c";
      else
      {
        ss << tstep[j-1];
        f << ss.str() << " fm/c" << sep;
        ss.str( "" );
      }
    }
    f << endl;
    break;
  case CharmNumber:
    f << "#time  # all charm pairs  # all without testparticles  dN_cpair/deta at mid.rap. betw. +-eta for  eta=2  eta=1.5  eta=1  eta=0.75  eta=0.5  eta=0.35  eta_spacetime = 0.5" << endl;
    break;
  case BottomNumber:
    f << "#time  # all bottom pairs  # all without testparticles  dN_bpair/deta at mid.rap. betw. +-eta for  eta=2  eta=1.5  eta=1  eta=0.75  eta=0.5  eta=0.35" << endl;
    break;
  case pos_mom:
    f << "#position in fm and momentum in GeV/c" << endl;
    f << "#X" << sep << "Y" << sep << "Z" << sep << "PX" << sep << "PY" << sep << "PZ" << endl;
    break;
  default: f << endl;
  }
  f << "#" << endl;
}


void analysis::addBoost32( const double beta )
{
  if ( beta <= 1.0 && beta >= 0.0 )
  {
    if ( FPT_COMP_E( beta, 1.0 ) )
    {
      ++boostDistribution32[numberBinsBoost - 1];
    }
    else
    {
      if ( FPT_COMP_E( beta, 1.0 ) )
        ++boostDistribution32[0];
      else
        ++boostDistribution32[int( beta*numberBinsBoost )];  // should be: beta / ( (1.0-0.0)/numberBinsBoost )
    }
  }
  else
    cout << "beta out of range" << endl;
}


void analysis::addBoost23( const double beta )
{
  if ( beta <= 1.0 && beta >= 0.0 )
  {
    if ( FPT_COMP_E( beta, 1.0 ) )
    {
      ++boostDistribution23[numberBinsBoost - 1];
    }
    else
    {
      if ( FPT_COMP_E( beta, 1.0 ) )
        ++boostDistribution23[0];
      else
        ++boostDistribution23[int( beta*numberBinsBoost )];  // should be: beta / ( (1.0-0.0)/numberBinsBoost )
    }
  }
  else
    cout << "beta out of range" << endl;
}


void analysis::registerProgressInfoForOutput( const double _time, const double _dt, const int _nTotalParticles, const int _nParticlesInFormationIntialProduction, const int _nParticlesInFormationGeometricCollisions, const int _nParticlesActive, const int _nParticlesInActiveCells, const vector< int >& _edgeCellSizes, const int _nFreeParticles, const int _nColl, const int _nColl22, const int _nColl23, const int _nColl32, const int _nColle )
{
  if ( progressLogFile.good() )
  {
    string sep = "\t";
    progressLogFile << _time << sep;
    progressLogFile << _dt << sep;
    progressLogFile << _nTotalParticles << sep;
    progressLogFile << _nParticlesInFormationIntialProduction << sep;
    progressLogFile << _nParticlesInFormationGeometricCollisions << sep;
    progressLogFile << _nParticlesActive << sep;
    progressLogFile << _nParticlesInActiveCells << sep;
    for ( unsigned int i = 0; i < _edgeCellSizes.size(); i++ )
    {
      progressLogFile << _edgeCellSizes[i] << sep;
    }
    progressLogFile << _nFreeParticles << sep;
    progressLogFile << _nColl << sep;
    progressLogFile << _nColl22 << sep;
    progressLogFile << _nColl23 << sep;
    progressLogFile << _nColl32 << sep;
    progressLogFile << _nColle << endl;
  }
  else
  {
    string errMsg = "error when attempting to write progress information log";
    throw eAnalysis_error( errMsg );
  }
}


void analysis::writeParticleAndCollisionNumbers(const double _time, const int _nTotalParticles, const int ncoll22, const int ncoll23, const int ncoll32, const int nGet32Errors, const int nGet23Errors, const int _nParticlesInFormationIntialProduction, const int _nParticlesInFormationGeometricCollisions, const int _nParticlesActive, const int _nParticlesInActiveCells, const std::vector< int >& _edgeCellSizes, const int _nFreeParticles, const double sin2Theta, const double vrelSum, const double sin2ThetaVrel, const double ncoll22MIDRAP)
{
  string sep = "\t";
  printParticleAndCollisionNumbers << _time << sep; //1
  printParticleAndCollisionNumbers << _nTotalParticles << sep; //2
  printParticleAndCollisionNumbers << ncoll22 << sep; //3
  printParticleAndCollisionNumbers << ncoll23 << sep;
  printParticleAndCollisionNumbers << ncoll32 << sep;
  printParticleAndCollisionNumbers << nGet32Errors << sep;
  printParticleAndCollisionNumbers << nGet23Errors << sep;
  printParticleAndCollisionNumbers << _nParticlesInFormationIntialProduction << sep;
  printParticleAndCollisionNumbers << _nParticlesInFormationGeometricCollisions << sep;
  printParticleAndCollisionNumbers << _nParticlesActive << sep;  //10
  printParticleAndCollisionNumbers << _nParticlesInActiveCells << sep;
   printParticleAndCollisionNumbers << _nFreeParticles << sep;  //12
  for ( unsigned int i = 0; i < _edgeCellSizes.size(); i++ )
  {
    printParticleAndCollisionNumbers << _edgeCellSizes[i] << sep;
  }
  printParticleAndCollisionNumbers << sin2ThetaVrel << sep; //21
  printParticleAndCollisionNumbers << vrelSum << sep;   //22
  printParticleAndCollisionNumbers << sin2Theta << sep; //23
  printParticleAndCollisionNumbers << ncoll22MIDRAP << sep; //24
  printParticleAndCollisionNumbers << endl;
  //cout << "Average sin^2 Theta = " << sin2Theta/ncoll22MIDRAP << endl;
  //cout << "Average vrel = " << vrelSum/ncoll22MIDRAP << endl;
  
}


// print dndy, detdy of light partons
void analysis::print_dndy( const int step )
{
  string subfix;
  stringstream ss;
  ss << tstep[step];
  
  double rapidityMax = theConfig->getMaximumEtaInitial() + 1;
  double rapidityMin = theConfig->getMinimumEtaInitial() - 1;
  
  if ( step == 0 )
    subfix = "initial";
  else if ( step == nTimeSteps )
    subfix = "final";
  else
    subfix = ss.str();
  
  string filename_prefix = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName();
  string filename;

  filename = filename_prefix + "_dndy_gluon_" + subfix;
  binning dndy_gluon(filename, rapidityMin, rapidityMax, 30);
  filename = filename_prefix + "_dedy_gluon_" + subfix;
  binningValues dedy_gluon(filename, rapidityMin, rapidityMax, 100);
  filename = filename_prefix + "_dndpt_gluon_" + subfix;
  binning dndpt_gluon(filename, 0.0, 20.0, 200);
  
  filename = filename_prefix + "_dndy_quark_" + subfix;
  binning dndy_quark(filename, rapidityMin, rapidityMax, 100);
  filename = filename_prefix + "_dedy_quark_" + subfix;
  binningValues dedy_quark(filename, rapidityMin, rapidityMax, 100);
  filename = filename_prefix + "_dndpt_quark_" + subfix;
  binning dndpt_quark(filename, 0.0, 20.0, 200);
  
  for ( unsigned int i = 0; i < particles.size(); i++ )
  {
    if(particles[i].Pos.T() > tstep[step])
    {
      continue;
    }
    const double eta = particles[i].Pos.Rapidity();
    const double pt = particles[i].Mom.Pt();;

    if( particles[i].FLAVOR == gluon )
    {
      dndy_gluon.add( eta );
      dedy_gluon.add( eta, pt );
      if( fabs( pt ) < 0.5 )
        dndpt_gluon.add( pt );
    }
    else if( particles[i].FLAVOR < 2 * Particle::max_N_light_flavor )
    {
      dndy_quark.add( eta );
      dedy_quark.add( eta, pt );
      if( fabs( pt ) < 0.5 )
        dndpt_quark.add( pt );
    }
  }
  
  dndy_gluon.print();
  //dedy_gluon.print();
  //dndpt_gluon.print();
  
  dndy_quark.print();
  //dedy_quark.print();
  //dndpt_quark.print();
}







void analysis::anaTimeStep_generalHIC()
{
  cout << "RHIC or LHC collision. Use general HIC output." << endl;     
  //---- times for output of data --------
  tstep[0] = 0.5;        //fm/c
  tstep[1] = 1.0;        //fm/c
  tstep[2] = 2.0;        //fm/c
  tstep[3] = 3.0;        //fm/c
  tstep[4] = 4.0;        //fm/c
  tstep[5] = 5.0;        //fm/c
  tstep[6] = 6.0;        //fm/c
  tstep[7] = 7.0;        //fm/c
  tstep[8] = 8.0;        //fm/c
  tstep[9] = 9.0;        //fm/c
  tstep[10] = 10.0;        //fm/c
  tstep[11] = ns_casc::infinity;   //fm/c
  nTimeSteps = 12;      //6 timesteps during the cascade, plus 1 before
  //--------------------------------------
}

void analysis::anaTimeStep_detailed()
{
  // by using this mode "generalOutput()" is called more often, output for more timesteps will be produced
  cout << "RHIC or LHC collision. Use detailed HIC output." << endl;
  //---- times for output of data --------
  
  for ( int ii = 0; ii < 10; ii++ ) // create time steps between [0.2:1.0], 0.2 step
  {
    tstep[ii] = ii*0.2 + 0.2; // fm/c
  }

  // max energy in midrapidity should be below freeze-out energy at 6.5 fm/c at RHIC energys
  // max energy in midrapidity should be below freeze-out energy at 8 fm/c at LHC energys
  for ( int ii = 0; ii < 34; ii++ ) // create time steps between [1.2:7.8], 0.2 step
  {
    tstep[ii+10] = ii*0.2 + 1.2; // fm/c
  }

  for ( int ii = 0; ii < 3; ii++ ) // create time steps between [8.0:10.0], 1.0 step
  {
    tstep[ii+10+34] = ii*1.0 + 8.0; // fm/c
  }

  // max freeze-out time for rapiditys < 1.5 for RHIC and LHC energys -> 16 fm/c
  // max freeze-out time for rapiditys < 3.0 for RHIC and LHC energys -> 32 fm/c
  for ( int ii = 0; ii < 15; ii++ ) // create time steps between [12.0:40.0], 2.0 step
  {
    tstep[ii+10+34+3] = ii*2.0 + 12.0; // fm/c
  }
  
  tstep[10+34+3+15] = ns_casc::infinity;   //fm/c
  nTimeSteps = 10+34+3+15+1;      //62 timesteps during the cascade, plus 1 before
  //--------------------------------------
}

void analysis::anaTimeStep_correlations()
{
//     tstep[0] = 0.1;      //fm/c
//     tstep[1] = 0.5;      //fm/c
//     tstep[2] = 1.0;      //fm/c
//     tstep[3] = 2.0;      //fm/c
//     tstep[4] = infinity;      //fm/c
//     nTimeSteps = 5;      //last tstep + 1
//     //tstep[4] = infinity;      //fm/c
//   
//   
  
//     cout << "Correlation analysis. Use 26 timesteps maximal." << endl;
//     // investigate charm yield at RHIC
//     //---- times for output of data --------
    tstep[0] = 0.5;      //fm/c
    tstep[1] = 1;      //fm/c
    tstep[2] = 1.5;      //fm/c
    tstep[3] = 2;      //fm/c
    tstep[4] = 2.5;      //fm/c
    tstep[5] = 3.;      //fm/c
    tstep[6] = 3.5;      //fm/c
    tstep[7] = 4;      //fm/c
    tstep[8] = 4.5;      //fm/c
    tstep[9] = 5.0;      //fm/c
    tstep[10] = 5.5;      //fm/c
    tstep[11] = 6.0;      //fm/c
    tstep[12] = 6.5;      //fm/c
    tstep[13] = 7.5;      //fm/c
    tstep[14] = infinity; //fm/c
    nTimeSteps = 15;      //last tstep + 1
    //--------------------------------------
    
    
  }

void analysis::anaTimeStepVideo()
{
  for ( int ii = 0; ii < 10000; ii++ ) // create time steps between [0.01:100.0], 0.01 step
  {
    tstep[ii] = ii*0.01 + 0.3; // fm/c
  }
  nTimeSteps = 10000+1; 
}
  
  
void analysis::anaTimeStep_heavyQuark()
{
    cout << "Central RHIC collision. Use 47 timesteps for heavy quark production output." << endl;
    // investigate charm yield at RHIC
    //---- times for output of data --------
    tstep[0] = 0.505;      //fm/c
    tstep[1] = 0.51;      //fm/c
    tstep[2] = 0.515;      //fm/c
    tstep[3] = 0.52;      //fm/c
    tstep[4] = 0.525;      //fm/c
    tstep[5] = 0.5275;      //fm/c
    tstep[6] = 0.53;      //fm/c
    tstep[7] = 0.5325;      //fm/c
    tstep[8] = 0.535;      //fm/c
    tstep[9] = 0.5375;      //fm/c
    tstep[10] = 0.54;      //fm/c
    tstep[11] = 0.5425;      //fm/c
    tstep[12] = 0.545;      //fm/c
    tstep[13] = 0.5475;      //fm/c
    tstep[14] = 0.55;      //fm/c
    tstep[15] = 0.55;      //fm/c
    tstep[16] = 0.6;      //fm/c
    tstep[17] = 0.65;      //fm/c
    tstep[18] = 0.7;      //fm/c
    tstep[19] = 0.75;      //fm/c
    tstep[20] = 0.8;      //fm/c
    tstep[21] = 0.85;      //fm/c
    tstep[22] = 0.9;      //fm/c
    tstep[23] = 0.95;      //fm/c
    tstep[24] = 1.0;      //fm/c
    tstep[25] = 1.1;      //fm/c
    tstep[26] = 1.2;      //fm/c
    tstep[27] = 1.3;      //fm/c
    tstep[28] = 1.4;      //fm/c
    tstep[29] = 1.5;      //fm/c
    tstep[30] = 1.6;      //fm/c
    tstep[31] = 1.7;      //fm/c
    tstep[32] = 1.8;      //fm/c
    tstep[33] = 1.9;      //fm/c
    tstep[34] = 2.0;      //fm/c
    tstep[35] = 2.5;      //fm/c
    tstep[36] = 3.0;      //fm/c
    tstep[37] = 3.5;      //fm/c
    tstep[38] = 4.0;      //fm/c
    tstep[39] = 4.3;      //fm/c
    tstep[40] = 4.6;      //fm/c
    tstep[41] = 5.0;      //fm/c
    tstep[42] = 5.3;      //fm/c
    tstep[43] = 5.6;      //fm/c
    tstep[44] = 6.0;      //fm/c
    tstep[45] = 6.3;      //fm/c
    tstep[46] = 6.6;      //fm/c
    tstep[47] = 7.0;      //fm/c
    tstep[48] = 7.3;      //fm/c
    tstep[49] = 7.6;
    tstep[50] = 8.0;     //fm/c
    tstep[51] = 8.3;      //fm/c
    tstep[52] = 8.6;      //fm/c
    tstep[53] = 9.0;      //fm/c
    tstep[54] = 9.3;      //fm/c
    tstep[55] = 9.6;      //fm/c
    tstep[56] = 10.0;      //fm/c
    tstep[57] = infinity; //fm/c
    nTimeSteps = 58;      //last tstep + 1
    //--------------------------------------
  }
 
void analysis::anaTimeStep_hydroParametrization()
{
  cout << "RHIC or LHC collision for hydro parametrization." << endl;  
  //---- times for output of data --------
  tstep[0] = 0.1;        //fm/c
  tstep[1] = 0.2;        //fm/c
  tstep[2] = 0.3;        //fm/c
  tstep[3] = 0.4;        //fm/c
  tstep[4] = 0.5;        //fm/c
  tstep[5] = 0.6;        //fm/c
  tstep[6] = 0.7;        //fm/c
  tstep[7] = 0.8;        //fm/c
  tstep[8] = 0.9;        //fm/c
  tstep[9] = 1.0;        //fm/c
  tstep[10] = 1.1;        //fm/c
  tstep[11] = 1.2;        //fm/c
  tstep[12] = 1.3;        //fm/c
  tstep[13] = 1.4;        //fm/c
  tstep[14] = 1.5;        //fm/c
  tstep[15] = 1.6;        //fm/c
  tstep[16] = 1.7;        //fm/c
  tstep[17] = 1.8;        //fm/c
  tstep[18] = 1.9;        //fm/c
  tstep[19] = 2.0;        //fm/c  
  tstep[20] = 2.1;        //fm/c
  tstep[21] = 2.2;        //fm/c
  tstep[22] = 2.3;        //fm/c
  tstep[23] = 2.4;        //fm/c
  tstep[24] = 2.5;        //fm/c
  tstep[25] = 2.6;        //fm/c
  tstep[26] = 2.7;        //fm/c
  tstep[27] = 2.8;        //fm/c
  tstep[28] = 2.9;        //fm/c
  tstep[29] = 3.0;        //fm/c   
  tstep[30] = 3.1;        //fm/c
  tstep[31] = 3.2;        //fm/c
  tstep[32] = 3.3;        //fm/c
  tstep[33] = 3.4;        //fm/c
  tstep[34] = 3.5;        //fm/c
  tstep[35] = 3.6;        //fm/c
  tstep[36] = 3.7;        //fm/c
  tstep[37] = 3.8;        //fm/c
  tstep[38] = 3.9;        //fm/c
  tstep[39] = 4.0;        //fm/c     
  tstep[40] = ns_casc::infinity;   //fm/c
  nTimeSteps = 41;      //6 timesteps during the cascade, plus 1 before
  //--------------------------------------    
}

//-----------------------------
//-----------------------------



void analysis::printHeaderForHydro(fstream & f, const anaType mode, const time_t end )
{
  vector<hydroParticleTypeProperties> vecTypeAna(theHydroParticleType.hydroParticleTypeVectorAnalysis);
  vector<hydroParticleTypeProperties> vecTypeCommon(theHydroParticleType.hydroParticleTypeVectorCommon);        
 
        //----------------------------//
  switch (mode)
  {
    case hydroNormal: f << "#fullHydroOutput " << endl;
      break;
    case hydroArrow: f << "#fullHydroArrowOutput " << endl;
      break;   
    case hydroParticleOutput: f << "#particleOutput " << endl;
      break;
    case hydroMidRap: f << "#fullHydroMidRapOutput " << endl;
      break;      
    case hydroMidRapArrow: f << "#fullHydroMidRapArrowOutput " << endl;
      break;         
    case midRapObservables: f << "#midRapObservables " << endl;
      break;      
    default: f << "#undefined "<< endl;
  }
    
  //----------------------------//
  
  f << "#\n";  
  f << "#runTime " << theConfig->getRuntime() << endl;
  f << "#nTimeSteps " << getRealNumberOfTimeStepsForAnaHydro() << endl;
  f << "#testparticles " << theConfig->getTestparticles() << endl;
  f << "#sideLengthX " << theConfig->getBoxLengthForHydroX() << endl;
  f << "#sideLengthY " << theConfig->getBoxLengthForHydroY() << endl;
  f << "#sideLengthZ " << theConfig->getBoxLengthForHydroZ() << endl;
  
  switch (mode)
  {
    case hydroNormal:
      f << "#columnWidthX " << theAnaHydro->anaHydroProfileNormal->columnWidthX << endl; 
      f << "#columnWidthY " << theAnaHydro->anaHydroProfileNormal->columnWidthY << endl; 
      f << "#columnWidthZ " << theAnaHydro->anaHydroProfileNormal->columnWidthZ << endl; 
      f << "#numberColumnX " << theConfig->getNumberColumnForHydroX() << endl; 
      f << "#numberColumnY " << theConfig->getNumberColumnForHydroY() << endl; 
      f << "#numberColumnZ " << theConfig->getNumberColumnForHydroZ() << endl;  
      break;
    case hydroArrow:
      f << "#columnWidthX " << theAnaHydro->anaHydroProfileArrow->columnWidthX << endl; 
      f << "#columnWidthY " << theAnaHydro->anaHydroProfileArrow->columnWidthY << endl; 
      f << "#columnWidthZ " << theAnaHydro->anaHydroProfileArrow->columnWidthZ << endl;    
      f << "#numberColumnX " << theConfig->getNumberColumnArrowForHydroX() << endl; 
      f << "#numberColumnY " << theConfig->getNumberColumnArrowForHydroY() << endl; 
      f << "#numberColumnZ " << theConfig->getNumberColumnArrowForHydroZ() << endl;       
      break;
    case hydroMidRap:
      f << "#columnWidthX " << theAnaHydro->anaHydroProfileMidRapNormal->columnWidthX << endl; 
      f << "#columnWidthY " << theAnaHydro->anaHydroProfileMidRapNormal->columnWidthY << endl; 
      f << "#columnWidthZ " << theAnaHydro->anaHydroProfileMidRapNormal->columnWidthZ << endl; 
      f << "#numberColumnX " << theConfig->getNumberColumnMidRapForHydroX() << endl; 
      f << "#numberColumnY " << theConfig->getNumberColumnMidRapForHydroY() << endl; 
      f << "#numberColumnZ " << 1 << endl;
      break;      
    case hydroMidRapArrow:
      f << "#columnWidthX " << theAnaHydro->anaHydroProfileMidRapArrow->columnWidthX << endl; 
      f << "#columnWidthY " << theAnaHydro->anaHydroProfileMidRapArrow->columnWidthY << endl; 
      f << "#columnWidthZ " << theAnaHydro->anaHydroProfileMidRapArrow->columnWidthZ << endl; 
      f << "#numberColumnX " << theConfig->getNumberColumnMidRapArrowForHydroX() << endl; 
      f << "#numberColumnY " << theConfig->getNumberColumnMidRapArrowForHydroY() << endl; 
      f << "#numberColumnZ " << 1 << endl;
      break;      
    default: break;
  }
  
  f << "#anaDim " << 3 << endl;//has to be 3
  f << "#sideDirection " << 3 << endl;//has to be 3
  f << "#\n";  
  f << "#FinalNInSimulation " << theAnaHydro->particleListgetFinalNumber() << endl;  
  f << "#\n"; 
  f << "#nParticleType " << vecTypeAna.size() << endl;
  f << "#\n";
  f << "#gc_correlationLength " << 1 << endl;  //has to be 1
  f << "#\n"; 
  f << "#ncoll22 " << ncoll22 << endl; 
  f << "#ncoll23 " << ncoll23 << endl;
  f << "#ncoll32 " << ncoll32 << endl;  
  f << "#\n";   
  f << "#startParticleTypeInfo" << endl;   
  for( unsigned int i = 0; i < vecTypeAna.size(); i++ )
  {
    f << "#particleTypeInfo" << i;
    f << " " << vecTypeAna[i].nameOfParticle;
    f << " " << vecTypeAna[i].degeneracyFactor;
    f << " " << vecTypeAna[i].mass;
    f << " " << vecTypeAna[i].masslessON;
    f << " " << vecTypeAna[i].flavorTypeStart;
    f << " " << vecTypeAna[i].flavorTypeEnd;
    f << endl;    
  }
  f << "#endParticleTypeInfo" << endl;  
  f << "#\n";    
  //-----------------------------//
  //some infos 
  f << "#\n";  
  f << "#startInfo" << endl;
  f << "#\n";     
  f << "#versionBAMPS: ";
  ns_casc::printSVNshort(f); //svn number 
  f << "#\n";    
  f << "#start: " << ctime(&start);
  f << "#finish: " << ctime(&end);  //will be replaced by end-time when destructor is called  
  f << "#\n";
  switch ( theConfig->getCrossSectionMethod() )
  {
    case csMethod_pQCD:
      f << "#pQCD is used: " << endl;
      //f << "#alpha_s: " << alpha_s << endl;
      break;
    case csMethod_constCS:
      f << "#a constant cross section is used: " << endl;
      f << "#crossSection [mb] " << theConfig->getInputCrossSectionValue() << endl;
      break;
    case csMethod_constEtaOverS:
      f << "#a constant eta/s is used: " << endl;
      f << "#etaOverS " << theConfig->getInputCrossSectionValue() << endl;
      break;    
    case csMethod_constMFP:
      f << "#a constant mean free path is used: " << endl;
      f << "#meanFreePath [fm] " << theConfig->getInputCrossSectionValue() << endl;
      break;    
    case csMethod_constMixtureCS:
      f << "#mixture cross sections are used: " << endl;                
      for(int i = 0; i <= 16; i++)
      {
        if(theHydroParticleType.getParticleType(i,vecTypeCommon) >= 0)
        {
          for(int j = i; j <= 16; j++)
          {
            if(theHydroParticleType.getParticleType(j,vecTypeCommon) >= 0)
            {
              scattering22_hydro collisions;
              string str;
              double cs22 = collisions.getMixtureXSection22(theHydroParticleType,vecTypeCommon,i,j,1.0,str);
              f << "#" << str << " " << cs22 * 10.0 * pow(0.197,2) << " mb  --  " << cs22 << " GeV^(-2)" << endl;
            }
            if(j == 1){j = 2;}//to skip the many quarks and show only once
            if(j == 3){j = 12;}//to skip the many quarks and show only once
          }
        }
        if(i == 1){i = 2;}//to skip the many quarks and show only once
        if(i == 3){i = 12;}//to skip the many quarks and show only once
      }
      break;
    case csMethod_constMixtureCS_scaledWithLambdaT2:    
      f << "#constant mixture cross section, but with a scaling factor fug*T^2 " << endl;
      for(int i = 0; i <= 16; i++)
      {
        if(theHydroParticleType.getParticleType(i,vecTypeCommon) >= 0)
        {
          for(int j = i; j <= 16; j++)
          {
            if(theHydroParticleType.getParticleType(j,vecTypeCommon) >= 0)
            {
              scattering22_hydro collisions;
              string str;
              double cs22 = collisions.getMixtureXSection22(theHydroParticleType,vecTypeCommon,i,j,1.0,str);
              f << "#" << str << " " << cs22 * 10.0 * pow(0.197,2) << " mb  --  " << cs22 << " GeV^(-2)" << endl;
            }
            if(j == 1){j = 2;}//to skip the many quarks and show only once
            if(j == 3){j = 12;}//to skip the many quarks and show only once
          }
        }
        if(i == 1){i = 2;}//to skip the many quarks and show only once
        if(i == 3){i = 12;}//to skip the many quarks and show only once
      }    
      break; 
    case csMethod_variableCS_scaledWithLambdaT2:    
      f << "  Special variable cross section with a scaling factor fug*T^2 " << endl;
      break;      
    default:
      cout << "#Error in method of cross section! ";
  }
  f << "#\n";   
  f << infoInitialParameters;
  f << "#\n"; 
  f << "#endInfo" << endl;
  f << "#\n";   
  f << "#end"<< endl;  
  f << "#initial seed for random generator ran2(): " << seed << endl;
  switch (mode)
  {
    case hydroNormal: f << "# t-1 \t x-2 \t y-3 \t z-4 \t T00-5 \t T11-6 \t T22-7 \t T33-8 \t T10-9 \t T20-10 \t T30-11 \t T21-12 \t T31-13 \t T32-14 \t N0-15 \t N1-16 \t N2-17 \t N3-18" << endl;
      break;
    case hydroArrow: f << "# t-1 \t x-2 \t y-3 \t z-4 \t T00-5 \t T11-6 \t T22-7 \t T33-8 \t T10-9 \t T20-10 \t T30-11 \t T21-12 \t T31-13 \t T32-14 \t N0-15 \t N1-16 \t N2-17 \t N3-18" << endl;
      break;      
    case hydroParticleOutput: f << "# X-1 \t Y-2 \t Z-3 \t PX-4 \t PY-5 \t PZ-6 \t E-7 \t FLAVOR-8" << endl;
      break;
    case hydroMidRap: f << "# t-1 \t x-2 \t y-3 \t rapZ-4 \t T00-5 \t T11-6 \t T22-7 \t T33-8 \t T10-9 \t T20-10 \t T30-11 \t T21-12 \t T31-13 \t T32-14 \t N0-15 \t N1-16 \t N2-17 \t N3-18" << endl;
      break;      
    case hydroMidRapArrow: f << "# t-1 \t x-2 \t y-3 \t rapZ-4 \t T00-5 \t T11-6 \t T22-7 \t T33-8 \t T10-9 \t T20-10 \t T30-11 \t T21-12 \t T31-13 \t T32-14 \t N0-15 \t N1-16 \t N2-17 \t N3-18" << endl;
      break;    
    case midRapObservables: f << "# t-1 \t dN/dy-2 \t dEt/dy-3 \t v2-4 \t v4-5" << endl;
      break;        
    default: f << endl;      
  }
  f << "#startReadIn\n";
}


void analysis::hydroDistributionOutput()
{
  string filename;
  
  if(theConfig->getAnalyseForHydro()){filename = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + ".f11";}
  else{filename = "/dev/null";}
  
  fstream file( filename.c_str(), ios::out | ios::trunc );
  
  time_t end;
  time( &end );  
  
  printHeaderForHydro(file,hydroNormal,end);

  double positionX,positionY,positionZ;
  
  anaHydroData * ad;  
  ad = new anaHydroData(*theAnaHydro->anaHydroProfileNormal);
  
  int nTypes = theAnaHydro->nTypes;
  int realTimeSteps = theAnaHydro->realTimeSteps;
  int numberColumnsX = ad->numberColumnsX;
  int numberColumnsY = ad->numberColumnsY;
  int numberColumnsZ = ad->numberColumnsZ; 
  double columnWidthX = ad->columnWidthX;
  double columnWidthY = ad->columnWidthY;
  double columnWidthZ = ad->columnWidthZ;  
  double boxLengthX = theAnaHydro->boxLengthX;
  double boxLengthY = theAnaHydro->boxLengthY;
  double boxLengthZ = theAnaHydro->boxLengthZ;   

  
  for ( int type = 0; type < nTypes; type++ )
  {    
    for (int nn = 0 ; nn < realTimeSteps; nn++)
      {
        for (int i = 0 ; i < numberColumnsX ; i++)
        {
          positionX = (double(i) / numberColumnsX) * boxLengthX - boxLengthX/2.0 + 0.5 * columnWidthX;  
          
          for (int j = 0 ; j < numberColumnsY ; j++)
          {
            positionY = (double(j) / numberColumnsY) * boxLengthY - boxLengthY/2.0 + 0.5 * columnWidthY;    
            
            for (int k = 0 ; k < numberColumnsZ ; k++)
            { 
              positionZ = (double(k) / numberColumnsZ) * boxLengthZ - boxLengthZ/2.0 + 0.5 * columnWidthZ; 
            
              int shift_nn = nn - 1;
              
              if(shift_nn<0){file << 0.0 << "\t";}
              else{file << tstep[shift_nn] << "\t";}
              file << positionX << "\t";
              file << positionY << "\t";
              file << positionZ << "\t";              
              file << ad->enDistTmunu[0][0][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[1][1][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[2][2][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[3][3][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[1][0][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[2][0][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[3][0][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[2][1][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[3][1][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[3][2][type][nn][i][j][k] << "\t";
              file << ad->enDistNmu[0][type][nn][i][j][k] << "\t";
              file << ad->enDistNmu[1][type][nn][i][j][k] << "\t";
              file << ad->enDistNmu[2][type][nn][i][j][k] << "\t";
              file << ad->enDistNmu[3][type][nn][i][j][k] << "\t";
              file << endl;
            }
          }
        }   
      file << endl;
      }
    file << endl;   
  }

  file.flush();
  file.close();
  
  //do I have to do that???
  delete ad;
  ad = NULL;
}

void analysis::hydroDistributionArrowOutput()
{
  string filename;
  
  if(theConfig->getAnalyseForHydro()){filename = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + ".f12";}
  else{filename = "/dev/null";}
  
  fstream file( filename.c_str(), ios::out | ios::trunc );
  
  time_t end;
  time( &end );  
  
  printHeaderForHydro(file,hydroArrow,end);

  double positionX,positionY,positionZ;
  anaHydroData * ad;  
  ad = new anaHydroData(*theAnaHydro->anaHydroProfileArrow);
  
  int nTypes = theAnaHydro->nTypes;
  int realTimeSteps = theAnaHydro->realTimeSteps;
  int numberColumnsX = ad->numberColumnsX;
  int numberColumnsY = ad->numberColumnsY;
  int numberColumnsZ = ad->numberColumnsZ; 
  double columnWidthX = ad->columnWidthX;
  double columnWidthY = ad->columnWidthY;
  double columnWidthZ = ad->columnWidthZ;  
  double boxLengthX = theAnaHydro->boxLengthX;
  double boxLengthY = theAnaHydro->boxLengthY;
  double boxLengthZ = theAnaHydro->boxLengthZ;     

  
  for ( int type = 0; type < nTypes; type++ )
  {    
    for (int nn = 0 ; nn < realTimeSteps; nn++)
      {
        for (int i = 0 ; i < numberColumnsX ; i++)
        {
          positionX = (double(i) / numberColumnsX) * boxLengthX - boxLengthX/2.0 + 0.5 * columnWidthX;  
          
          for (int j = 0 ; j < numberColumnsY ; j++)
          {
            positionY = (double(j) / numberColumnsY) * boxLengthY - boxLengthY/2.0 + 0.5 * columnWidthY;    
            
            for (int k = 0 ; k < numberColumnsZ ; k++)
            { 
              positionZ = (double(k) / numberColumnsZ) * boxLengthZ - boxLengthZ/2.0 + 0.5 * columnWidthZ; 
            
              int shift_nn = nn - 1;
              
              if(shift_nn<0){file << 0.0 << "\t";}
              else{file << tstep[shift_nn] << "\t";}
              file << positionX << "\t";
              file << positionY << "\t";
              file << positionZ << "\t";              
              file << ad->enDistTmunu[0][0][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[1][1][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[2][2][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[3][3][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[1][0][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[2][0][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[3][0][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[2][1][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[3][1][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[3][2][type][nn][i][j][k] << "\t";
              file << ad->enDistNmu[0][type][nn][i][j][k] << "\t";
              file << ad->enDistNmu[1][type][nn][i][j][k] << "\t";
              file << ad->enDistNmu[2][type][nn][i][j][k] << "\t";
              file << ad->enDistNmu[3][type][nn][i][j][k] << "\t";
              file << endl;
            }
          }
        }   
      file << endl;
      }
    file << endl;   
  }

  file.flush();
  file.close();
  
  //do I have to do that???
  delete ad;
  ad = NULL;
}



void analysis::hydroDistributionMidRapOutput()
{
  string filename;
  
  if(theConfig->getAnalyseForHydro()){filename = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + ".f16";}
  else{filename = "/dev/null";}
  
  fstream file( filename.c_str(), ios::out | ios::trunc );
  
  time_t end;
  time( &end );  
  
  printHeaderForHydro(file,hydroMidRap,end);

  double positionX,positionY,positionZ;
  
  anaHydroData * ad;  
  ad = new anaHydroData(*theAnaHydro->anaHydroProfileMidRapNormal);
  
  int nTypes = theAnaHydro->nTypes;
  int realTimeSteps = theAnaHydro->realTimeSteps;
  int numberColumnsX = ad->numberColumnsX;
  int numberColumnsY = ad->numberColumnsY;
  int numberColumnsZ = ad->numberColumnsZ; 
  double columnWidthX = ad->columnWidthX;
  double columnWidthY = ad->columnWidthY;
  double columnWidthZ = ad->columnWidthZ;  
  double boxLengthX = theAnaHydro->boxLengthX;
  double boxLengthY = theAnaHydro->boxLengthY;
  double boxLengthZ = theAnaHydro->boxLengthZ;   

  
  for ( int type = 0; type < nTypes; type++ )
  {    
    for (int nn = 0 ; nn < realTimeSteps; nn++)
      {
        for (int i = 0 ; i < numberColumnsX ; i++)
        {
          positionX = (double(i) / numberColumnsX) * boxLengthX - boxLengthX/2.0 + 0.5 * columnWidthX;  
          
          for (int j = 0 ; j < numberColumnsY ; j++)
          {
            positionY = (double(j) / numberColumnsY) * boxLengthY - boxLengthY/2.0 + 0.5 * columnWidthY;    
            
            for (int k = 0 ; k < numberColumnsZ ; k++)
            { 
              positionZ = 0.0; 
            
              int shift_nn = nn - 1;
              
              if(shift_nn<0){file << 0.0 << "\t";}
              else{file << tstep[shift_nn] << "\t";}
              file << positionX << "\t";
              file << positionY << "\t";
              file << positionZ << "\t";              
              file << ad->enDistTmunu[0][0][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[1][1][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[2][2][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[3][3][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[1][0][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[2][0][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[3][0][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[2][1][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[3][1][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[3][2][type][nn][i][j][k] << "\t";
              file << ad->enDistNmu[0][type][nn][i][j][k] << "\t";
              file << ad->enDistNmu[1][type][nn][i][j][k] << "\t";
              file << ad->enDistNmu[2][type][nn][i][j][k] << "\t";
              file << ad->enDistNmu[3][type][nn][i][j][k] << "\t";
              file << endl;
            }
          }
        }   
      file << endl;
      }
    file << endl;   
  }

  file.flush();
  file.close();
  
  //do I have to do that???
  delete ad;
  ad = NULL;
}


void analysis::hydroDistributionMidRapArrowOutput()
{
  string filename;
  
  if(theConfig->getAnalyseForHydro()){filename = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + ".f17";}
  else{filename = "/dev/null";}
  
  fstream file( filename.c_str(), ios::out | ios::trunc );
  
  time_t end;
  time( &end );  
  
  printHeaderForHydro(file,hydroMidRapArrow,end);

  
  double positionX,positionY,positionZ;
  anaHydroData * ad;  
  ad = new anaHydroData(*theAnaHydro->anaHydroProfileMidRapArrow);
  
  int nTypes = theAnaHydro->nTypes;
  int realTimeSteps = theAnaHydro->realTimeSteps;
  int numberColumnsX = ad->numberColumnsX;
  int numberColumnsY = ad->numberColumnsY;
  int numberColumnsZ = ad->numberColumnsZ; 
  double columnWidthX = ad->columnWidthX;
  double columnWidthY = ad->columnWidthY;
  double columnWidthZ = ad->columnWidthZ;  
  double boxLengthX = theAnaHydro->boxLengthX;
  double boxLengthY = theAnaHydro->boxLengthY;
  double boxLengthZ = theAnaHydro->boxLengthZ;   

  
  for ( int type = 0; type < nTypes; type++ )
  {    
    for (int nn = 0 ; nn < realTimeSteps; nn++)
      {
        for (int i = 0 ; i < numberColumnsX ; i++)
        {
          positionX = (double(i) / numberColumnsX) * boxLengthX - boxLengthX/2.0 + 0.5 * columnWidthX;  
          
          for (int j = 0 ; j < numberColumnsY ; j++)
          {
            positionY = (double(j) / numberColumnsY) * boxLengthY - boxLengthY/2.0 + 0.5 * columnWidthY;    
            
            for (int k = 0 ; k < numberColumnsZ ; k++)
            { 
              positionZ = 0.0; 
            
              int shift_nn = nn - 1;
              
              if(shift_nn<0){file << 0.0 << "\t";}
              else{file << tstep[shift_nn] << "\t";}
              file << positionX << "\t";
              file << positionY << "\t";
              file << positionZ << "\t";              
              file << ad->enDistTmunu[0][0][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[1][1][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[2][2][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[3][3][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[1][0][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[2][0][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[3][0][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[2][1][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[3][1][type][nn][i][j][k] << "\t";
              file << ad->enDistTmunu[3][2][type][nn][i][j][k] << "\t";
              file << ad->enDistNmu[0][type][nn][i][j][k] << "\t";
              file << ad->enDistNmu[1][type][nn][i][j][k] << "\t";
              file << ad->enDistNmu[2][type][nn][i][j][k] << "\t";
              file << ad->enDistNmu[3][type][nn][i][j][k] << "\t";
              file << endl;
            }
          }
        }   
      file << endl;
      }
    file << endl;   
  }

  file.flush();
  file.close();
  
  //do I have to do that???
  delete ad;
  ad = NULL;
}


void analysis::hydroMidRapidityObservablesOutput()
{
  string filename;
  
  if(theConfig->getAnalyseForHydro()){filename = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + ".f51";}
  else{filename = "/dev/null";}
  
  fstream file( filename.c_str(), ios::out | ios::trunc );
  
  time_t end;
  time( &end );  
  
  printHeaderForHydro(file,midRapObservables,end);

  anaHydroData * ad;  
  ad = new anaHydroData(*theAnaHydro->anaHydroProfileMidRapNormal);
  
  int nTypes = theAnaHydro->nTypes;
  int realTimeSteps = theAnaHydro->realTimeSteps;
  
  for ( int type = 0; type <= nTypes; type++ )
  {    
    for (int nn = 0 ; nn < realTimeSteps; nn++)
      {
              int shift_nn = nn - 1;
              
              if(shift_nn<0){file << 0.0 << "\t";}
              else{file << tstep[shift_nn] << "\t";}      
              file << ad->midRapOb_dN_dy[type][nn] << "\t";
              file << ad->midRapOb_dEt_dy[type][nn] << "\t";
              file << ad->midRapOb_v2[type][nn] << "\t";
              file << ad->midRapOb_v4[type][nn] << "\t";
              file << endl;
      }
    file << endl;   
  }

  file.flush();
  file.close();
  
  //do I have to do that???
  delete ad;
  ad = NULL;
}



void analysis::particleListOutput()
{
  string filename;
  
  if(theConfig->getAnalyseForHydro()){filename = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + ".f13";}
  else{filename = "/dev/null";}
  
  fstream file( filename.c_str(), ios::out | ios::trunc );
  
  time_t end;
  time( &end );  
  
  printHeaderForHydro(file,hydroParticleOutput,end);
  
  vector<hydroParticleTypeProperties> vecTypeAna(theHydroParticleType.hydroParticleTypeVectorAnalysis);    
  int finalNumberInSimulation = 0;
  
  for( unsigned int i=0; i< particles.size(); i++ )
    {
      //check whether particles are in the unvisible box :-)
      if( fabs(particles[i].Pos.X()) < 0.5 * theAnaHydro->boxLengthX &&
          fabs(particles[i].Pos.Y()) < 0.5 * theAnaHydro->boxLengthY &&
          fabs(particles[i].Pos.Z()) < 0.5 * theAnaHydro->boxLengthZ  )
      { 
        //in case you have more particle species  
        int type = theHydroParticleType.getParticleType(particles[i].FLAVOR,vecTypeAna);
        
        if(type >= 0)
        {
          finalNumberInSimulation++;
          file << particles[i].Pos.X() <<"\t";
          file << particles[i].Pos.Y() <<"\t";  
          file << particles[i].Pos.Z() <<"\t";
          file << particles[i].Mom.Px() <<"\t";
          file << particles[i].Mom.Py() <<"\t";
          file << particles[i].Mom.Pz() <<"\t";
          file << particles[i].Mom.E() <<"\t";
          file << particles[i].FLAVOR <<"\t";
          file << endl;
        }
      }
    }
    
    if(finalNumberInSimulation != theAnaHydro->finalNumberInSimulation)
    {
     cout << "Error in final number for particle list output in hydro!" << endl; 
    }

  file.flush();
  file.close();
}



int analysis::getRealNumberOfTimeStepsForAnaHydro()
  {
    //calculate the number of timesteps for hydro analysis
    int nn = 0;
    
    nn += 1;//because fo initial step
    
    for(int i = 0; i < nTimeSteps; i++)
    {
      double tStepInBAMPS = tstep[i];
      
      if( theConfig->getRuntime() > tStepInBAMPS )
      {
        nn++;
      }
      else
      {
        break;
      }
  
    }
//     cout << "------------------------------------------- " << endl;    
//     cout << "  Number of time steps for Hydro analysis: " << nn << endl;
//     cout << "------------------------------------------- " << endl;      
//     
    return nn;
  }
  

 //**************************************************************************
// Fragmentation functions
// using KKP from http://www.desy.de/~poetter/kkp.html
// rewritten into C++ by Bjoern Schenke 2013
//**************************************************************************

//=====================================================================
//
//     ------------------------------------------------------------
//     Fragmentation functions for: Pions, Kaons, Protons, Neutrons
//            (includes mass-threshholds for c and b quarks)
//
//     Reference: B.A.Kniehl, G.Kramer, B.Potter, NPB582 (2000) 514
//     ------------------------------------------------------------
//
//     ih, iset, x, qs are input; dh is output.
//     ih   = 1 : (pi^+ + pi^-)  /2
//     ih   = 2 : (K^+ + K^-)    /2
//     ih   = 3 : (K^0 + K^0_bar)/2
//     ih   = 4 : (p + p_bar)    /2
//     ih   = 5 : (pi^0)
//     ih   = 6 : (n + n_bar)    /2
//     ih   = 7 : (h^+ + h^-)         [as sum of pions, kaons and protons]
//
//     iset = 0 : LO
//     iset = 1 : NLO
//
//     x    = longitudinal-momentum fraction
//     qs   = fragmentation scale (in GeV)
//
//     Parton label:
//     0    1    2    3    4    5    6    7    8     9    10
//     g    u   ubar  d   dbar  s   sbar  c   cbar   b   bbar
//
//     Lambda_QCD (in GeV):
//     0.088 in LO
//     0.213 in NLO
//
//=====================================================================

double analysis::KKPFragmentation(int ih, int iset, double x, double qs, FLAVOR_TYPE species)
{
    //(a-h,o-z)
    
    double rlam, s, sc, sb;
    // --- Mass-thresholds:
    double rmcc = 2.9788;
    double rmbb = 9.46037;
    // --- Q_0 (in GeV):
    double q0= sqrt(2.);
    // --- BAK
    if(qs < q0)
        qs = q0;
    
    double b[4];
    double a[12];
    double dh[11];
    
    double dpg, dpu, dps, dpc, dpb;
    double dkg, dku, dks, dkc, dkb;
    double dprg, dpru, dprs, dprc, dprb;
    double dkd;
    
    // --- LO FFs
    if (iset == 0)
    {
        rlam= 0.088;
        s=  log(log(qs*qs/(rlam*rlam))/log(q0*q0/(rlam*rlam)));
        sc= log(log(qs*qs/(rlam*rlam))/log(rmcc*rmcc/(rlam*rlam)));
        sb= log(log(qs*qs/(rlam*rlam))/log(rmbb*rmbb/(rlam*rlam)));
        
        //      cout << "s=" << s << endl;
        //       cout << "sc=" << sc << endl;
        //       cout << "sb=" << sb << endl;
        
        
        // ---------------------- LO PION ------------------------------
        b[1]=                            6.04510;
        b[2]=                           -0.71378;
        b[3]=                            2.92133;
        a[1]=                           -6.61523;
        a[2]=                           -1.64978;
        a[3]=                            2.68223;
        a[4]=                            0.14705;
        a[5]=                           -1.08423;
        a[6]=                           -0.43182;
        a[7]=                            1.48429;
        a[8]=                            1.32887;
        a[9]=                           -1.78696;
        a[10]=                            0.23086;
        a[11]=                           -0.29182;
        dpg = (b[1] +a[1]*s + a[2]*s*s +a[3]*s*s*s)*pow(x,(b[2] +a[4]*s + a[5]*s*s +a[6]*s*s*s))
        *pow((1.-x),(b[3] +a[7]*s + a[8]*s*s +a[9]*s*s*s))*(1.+(a[10]*s +a[11]*s*s)/x);
        b[1]=                           0.54610;
        b[2]=                          -1.46616;
        b[3]=                           1.01864;
        a[1]=                           -0.22946;
        a[2]=                           -0.22594;
        a[3]=                            0.21119;
        a[4]=                           -0.45404;
        a[5]=                           -0.12684;
        a[6]=                            0.27646;
        a[7]=                            0.95367;
        a[8]=                           -1.09835;
        a[9]=                            0.74657;
        a[10]=                           -0.01877;
        a[11]=                            0.02949;
        dpu = (b[1] +a[1]*s + a[2]*s*s +a[3]*s*s*s)*pow(x,(b[2] +a[4]*s + a[5]*s*s +a[6]*s*s*s))
        *pow((1.-x),(b[3] +a[7]*s + a[8]*s*s +a[9]*s*s*s))*(1.+(a[10]*s +a[11]*s*s)/x);
        //      dpu= (b[1] +a[1]*s + a[2]*s**2 +a[3]*s**3)*x**
        //  .        (b[2] +a[4]*s + a[5]*s**2 +a[6]*s**3)*(1d0-x)**
        //  .        (b[3] +a[7]*s + a[8]*s**2 +a[9]*s**3)*
        //  .        (1d0+(a[10]*s +a[11]*s**2)/x)
        b[1]=                            22.2815;
        b[2]=                            0.12732;
        b[3]=                            6.13697;
        a[1]=                           -20.8125;
        a[2]=                           -11.5725;
        a[3]=                            15.5372;
        a[4]=                            0.23075;
        a[5]=                           -2.71424;
        a[6]=                            1.72456;
        a[7]=                            2.18849;
        a[8]=                           -5.04475;
        a[9]=                            3.29117;
        a[10]=                            0.09044;
        a[11]=                           -0.07589;
        dps = (b[1] +a[1]*s + a[2]*s*s +a[3]*s*s*s)*pow(x,(b[2] +a[4]*s + a[5]*s*s +a[6]*s*s*s))
        *pow((1.-x),(b[3] +a[7]*s + a[8]*s*s +a[9]*s*s*s))*(1.+(a[10]*s +a[11]*s*s)/x);
        //          dps= (b[1] +a[1]*s + a[2]*s**2 +a[3]*s**3)*x**
        //      .        (b[2] +a[4]*s + a[5]*s**2 +a[6]*s**3)*(1d0-x)**
        //      .        (b[3] +a[7]*s + a[8]*s**2 +a[9]*s**3)*
        //      .        (1d0+(a[10]*s +a[11]*s**2)/x)
        b[1]=                            8.75500;
        b[2]=                           -0.38611;
        b[3]=                            5.61846;
        a[1]=                           -9.32277;
        a[2]=                            1.80600;
        a[3]=                            2.02179;
        a[4]=                           -0.41190;
        a[5]=                           -0.48496;
        a[6]=                            0.42525;
        a[7]=                            0.74035;
        a[8]=                           -0.64929;
        a[9]=                            0.66788;
        a[10]=                            0.06652;
        a[11]=                           -0.05531;
        dpc = (b[1] +a[1]*sc + a[2]*sc*sc +a[3]*sc*sc*sc)*pow(x,(b[2] +a[4]*sc + a[5]*sc*sc +a[6]*sc*sc*sc))
        *pow((1.-x),(b[3] +a[7]*sc + a[8]*sc*sc +a[9]*sc*sc*sc))*(1.+(a[10]*sc +a[11]*sc*sc)/x);
        //       dpc= (b[1] +a[1]*sc + a[2]*sc**2 +a[3]*sc**3)*x**
        //      .        (b[2] +a[4]*sc + a[5]*sc**2 +a[6]*sc**3)*(1d0-x)**
        //      .        (b[3] +a[7]*sc + a[8]*sc**2 +a[9]*sc**3)*
        //      .        (1d0+(a[10]*sc +a[11]*sc**2)/x)
        b[1]=                            0.31147;
        b[2]=                           -1.92993;
        b[3]=                            3.47086;
        a[1]=                           -0.19319;
        a[2]=                           -0.10487;
        a[3]=                            0.18824;
        a[4]=                           -0.44692;
        a[5]=                           -0.08271;
        a[6]=                            0.30441;
        a[7]=                            0.79775;
        a[8]=                           -0.28091;
        a[9]=                            0.39504;
        a[10]=                           -0.04887;
        a[11]=                            0.03212;
        dpb = (b[1] +a[1]*sb + a[2]*sb*sb +a[3]*sb*sb*sb)*pow(x,(b[2] +a[4]*sb + a[5]*sb*sb +a[6]*sb*sb*sb))
        *pow((1.-x),(b[3] +a[7]*sb + a[8]*sb*sb +a[9]*sb*sb*sb))*(1.+(a[10]*sb +a[11]*sb*sb)/x);
        //          dpb= (b[1] +a[1]*sb + a[2]*sb**2 +a[3]*sb**3)*x**
        //      .        (b[2] +a[4]*sb + a[5]*sb**2 +a[6]*sb**3)*(1d0-x)**
        //      .        (b[3] +a[7]*sb + a[8]*sb**2 +a[9]*sb**3)
        //      .        *(1d0+(a[10]*sb +a[11]*sb**2)/x)
        
        // ---------------------- LO KAON ------------------------------
        b[1]=                            0.02862;
        b[2]=                           -2.94091;
        b[3]=                            2.73474;
        a[1]=                           -0.02113;
        a[2]=                            0.00389;
        a[3]=                            0.00901;
        a[4]=                            0.66881;
        a[5]=                           -0.29670;
        a[6]=                            0.20574;
        a[7]=                           -0.58222;
        a[8]=                            0.04329;
        a[9]=                            0.78033;
        a[10]=                            0.03586;
        a[11]=                           -0.01220;
        dkg = (b[1] +a[1]*s + a[2]*s*s +a[3]*s*s*s)*pow(x,(b[2] +a[4]*s + a[5]*s*s +a[6]*s*s*s))
        *pow((1.-x),(b[3] +a[7]*s + a[8]*s*s +a[9]*s*s*s))*(1.+(a[10]*s +a[11]*s*s)/x);
        //       dkg= (b[1] +a[1]*s + a[2]*s**2 +a[3]*s**3)*x**
        //      .        (b[2] +a[4]*s + a[5]*s**2 +a[6]*s**3)*(1d0-x)**
        //      .        (b[3] +a[7]*s + a[8]*s**2 +a[9]*s**3)
        //      .        *(1d0+(a[10]*s +a[11]*s**2)/x)
        b[1]=                           0.25937;
        b[2]=                          -0.61925;
        b[3]=                           0.85946;
        a[1]=                           -0.10502;
        a[2]=                            0.00572;
        a[3]=                           -0.00269;
        a[4]=                            0.09956;
        a[5]=                            0.07389;
        a[6]=                           -0.00070;
        a[7]=                            0.57965;
        a[8]=                            0.26397;
        a[9]=                           -0.12764;
        a[10]=                            0.15303;
        a[11]=                            0.14807;
        dku = (b[1] +a[1]*s + a[2]*s*s +a[3]*s*s*s)*pow(x,(b[2] +a[4]*s + a[5]*s*s +a[6]*s*s*s))
        *pow((1.-x),(b[3] +a[7]*s + a[8]*s*s +a[9]*s*s*s))*(1.+(a[10]*s +a[11]*s*s)/x);
        //      dku= (b[1] +a[1]*s + a[2]*s**2 +a[3]*s**3)*x**
        //  .        (b[2] +a[4]*s + a[5]*s**2 +a[6]*s**3)*(1d0-x)**
        //  .        (b[3] +a[7]*s + a[8]*s**2 +a[9]*s**3)
        //  .        *(1d0+(a[10]*s +a[11]*s**2)/x)
        b[1]=                            5.38115;
        b[2]=                           -0.00321;
        b[3]=                            3.07632;
        a[1]=                           -3.05084;
        a[2]=                           -1.10056;
        a[3]=                            1.31207;
        a[4]=                           -0.25889;
        a[5]=                           -0.18494;
        a[6]=                            0.13994;
        a[7]=                            1.13745;
        a[8]=                           -0.90413;
        a[9]=                            0.56581;
        a[10]=                            0.05141;
        a[11]=                           -0.00697;
        dkd = (b[1] +a[1]*s + a[2]*s*s +a[3]*s*s*s)*pow(x,(b[2] +a[4]*s + a[5]*s*s +a[6]*s*s*s))
        *pow((1.-x),(b[3] +a[7]*s + a[8]*s*s +a[9]*s*s*s))*(1.+(a[10]*s +a[11]*s*s)/x);
        //       dkd= (b[1] +a[1]*s + a[2]*s**2 +a[3]*s**3)*x**
        //  .        (b[2] +a[4]*s + a[5]*s**2 +a[6]*s**3)*(1d0-x)**
        //  .        (b[3] +a[7]*s + a[8]*s**2 +a[9]*s**3)
        //  .        *(1d0+(a[10]*s +a[11]*s**2)/x)
        b[1]=                            5.18266;
        b[2]=                           -0.17751;
        b[3]=                            4.30306;
        a[1]=                           -3.48519;
        a[2]=                           -1.00982;
        a[3]=                            1.17996;
        a[4]=                            0.02309;
        a[5]=                           -0.61327;
        a[6]=                           -0.03532;
        a[7]=                            1.00547;
        a[8]=                           -0.51779;
        a[9]=                            0.20683;
        a[10]=                            0.13514;
        a[11]=                           -0.17778;
        dkc = (b[1] +a[1]*sc + a[2]*sc*sc +a[3]*sc*sc*sc)*pow(x,(b[2] +a[4]*sc + a[5]*sc*sc +a[6]*sc*sc*sc))
        *pow((1.-x),(b[3] +a[7]*sc + a[8]*sc*sc +a[9]*sc*sc*sc))*(1.+(a[10]*sc +a[11]*sc*sc)/x);
        //       dkc= (b[1] +a[1]*sc + a[2]*sc**2 +a[3]*sc**3)*x**
        //      .        (b[2] +a[4]*sc + a[5]*sc**2 +a[6]*sc**3)*(1d0-x)**
        //      .        (b[3] +a[7]*sc + a[8]*sc**2 +a[9]*sc**3)
        //      .        *(1d0+(a[10]*sc +a[11]*sc**2)/x)
        b[1]=                            1.57044;
        b[2]=                           -0.84143;
        b[3]=                            6.01488;
        a[1]=                           -1.78340;
        a[2]=                            0.57100;
        a[3]=                            0.15469;
        a[4]=                           -0.43448;
        a[5]=                           -0.05314;
        a[6]=                           -0.36621;
        a[7]=                            0.72953;
        a[8]=                           -0.64433;
        a[9]=                            0.92351;
        a[10]=                            0.01024;
        a[11]=                           -0.06160;
        dkb = (b[1] +a[1]*sb + a[2]*sb*sb +a[3]*sb*sb*sb)*pow(x,(b[2] +a[4]*sb + a[5]*sb*sb +a[6]*sb*sb*sb))
        *pow((1.-x),(b[3] +a[7]*sb + a[8]*sb*sb +a[9]*sb*sb*sb))*(1.+(a[10]*sb +a[11]*sb*sb)/x);
        //       dkb= (b[1] +a[1]*sb + a[2]*sb**2 +a[3]*sb*sb*sb)*x**
        //      .        (b[2] +a[4]*sb + a[5]*sb**2 +a[6]*sb**3)*(1d0-x)**
        //      .        (b[3] +a[7]*sb + a[8]*sb**2 +a[9]*sb**3)
        //      .        *(1d0+(a[10]*sb +a[11]*sb**2)/x)
        
        // ---------------------- LO PROTON -----------------------------
        b[1]=                            0.73953;
        b[2]=                           -0.76986;
        b[3]=                            7.69079;
        a[1]=                           -1.64519;
        a[2]=                            1.01189;
        a[3]=                           -0.10175;
        a[4]=                           -3.58787;
        a[5]=                            13.8025;
        a[6]=                           -13.8902;
        a[7]=                           -2.84470;
        a[8]=                           -0.36719;
        a[9]=                           -2.21825;
        a[10]=                            1.26515;
        a[11]=                           -1.96117;
        dprg = (b[1] +a[1]*s + a[2]*s*s +a[3]*s*s*s)*pow(x,(b[2] +a[4]*s + a[5]*s*s +a[6]*s*s*s))
        *pow((1.-x),(b[3] +a[7]*s + a[8]*s*s +a[9]*s*s*s))*(1.+(a[10]*s +a[11]*s*s + 0.54769*s*s*s)/x);
        //          dprg= (b[1] +a[1]*s + a[2]*s**2 +a[3]*s**3)*x**
        //      .        (b[2] +a[4]*s + a[5]*s**2 +a[6]*s**3)*(1d0-x)**
        //      .        (b[3] +a[7]*s + a[8]*s**2 +a[9]*s**3)
        //      .        *(1d0+(a[10]*s +a[11]*s**2 + 0.54769*s**3)/x)
        b[1]=                           0.40211;
        b[2]=                          -0.85973;
        b[3]=                           2.80160;
        a[1]=                           -0.21633;
        a[2]=                           -0.07045;
        a[3]=                            0.07831;
        a[4]=                            0.13987;
        a[5]=                           -0.82412;
        a[6]=                            0.43114;
        a[7]=                            0.78923;
        a[8]=                           -0.05344;
        a[9]=                            0.01460;
        a[10]=                            0.05198;
        a[11]=                           -0.04623;
        dpru = (b[1] +a[1]*s + a[2]*s*s +a[3]*s*s*s)*pow(x,(b[2] +a[4]*s + a[5]*s*s +a[6]*s*s*s))
        *pow((1.-x),(b[3] +a[7]*s + a[8]*s*s +a[9]*s*s*s))*(1.+(a[10]*s +a[11]*s*s)/x);
        //          dpru= (b[1] +a[1]*s + a[2]*s**2 +a[3]*s**3)*x**
        //      .        (b[2] +a[4]*s + a[5]*s**2 +a[6]*s**3)*(1d0-x)**
        //      .        (b[3] +a[7]*s + a[8]*s**2 +a[9]*s**3)
        //      .        *(1d0+(a[10]*s +a[11]*s**2]/x)
        b[1]=                            4.07885;
        b[2]=                           -0.09735;
        b[3]=                            4.99191;
        a[1]=                           -2.97392;
        a[2]=                           -0.92973;
        a[3]=                            1.23517;
        a[4]=                            0.25834;
        a[5]=                           -1.52246;
        a[6]=                            0.77060;
        a[7]=                            1.14379;
        a[8]=                           -0.85320;
        a[9]=                            0.45607;
        a[10]=                            0.07174;
        a[11]=                           -0.08321;
        dprs = (b[1] +a[1]*s + a[2]*s*s +a[3]*s*s*s)*pow(x,(b[2] +a[4]*s + a[5]*s*s +a[6]*s*s*s))
        *pow((1.-x),(b[3] +a[7]*s + a[8]*s*s +a[9]*s*s*s))*(1.+(a[10]*s +a[11]*s*s)/x);
        //          dprs= (b[1] +a[1]*s + a[2]*s**2 +a[3]*s**3)*x**
        //      .        (b[2] +a[4]*s + a[5]*s**2 +a[6]*s**3)*(1d0-x)**
        //      .        (b[3] +a[7]*s + a[8]*s**2 +a[9]*s**3)
        //      .        *(1d0+(a[10]*s +a[11]*s**2]/x)
        b[1]=                            0.11061;
        b[2]=                           -1.54340;
        b[3]=                            2.20681;
        a[1]=                           -0.07726;
        a[2]=                            0.05422;
        a[3]=                           -0.03364;
        a[4]=                           -0.20804;
        a[5]=                            0.29038;
        a[6]=                           -0.23662;
        a[7]=                            0.62274;
        a[8]=                            0.29713;
        a[9]=                           -0.21861;
        a[10]=                            0.00831;
        a[11]=                            0.00065;
        dprc = (b[1] +a[1]*sc + a[2]*sc*sc +a[3]*sc*sc*sc)*pow(x,(b[2] +a[4]*sc + a[5]*sc*sc +a[6]*sc*sc*sc))
        *pow((1.-x),(b[3] +a[7]*sc + a[8]*sc*sc +a[9]*sc*sc*sc))*(1.+(a[10]*sc +a[11]*sc*sc)/x);
        //       dprc= (b[1] +a[1]*sc + a[2]*sc**2 +a[3]*sc**3)*x**
        //      .        (b[2] +a[4]*sc + a[5]*sc**2 +a[6]*sc**3)*(1d0-x)**
        //      .        (b[3] +a[7]*sc + a[8]*sc**2 +a[9]*sc**3)
        //      .        *(1d0+(a[10]*sc +a[11]*sc**2]/x)
        b[1]=                            40.0971;
        b[2]=                            0.74249;
        b[3]=                            12.3729;
        a[1]=                           -123.531;
        a[2]=                            128.666;
        a[3]=                           -29.1808;
        a[4]=                           -1.29639;
        a[5]=                           -3.65003;
        a[6]=                            3.05340;
        a[7]=                           -1.04932;
        a[8]=                            0.34662;
        a[9]=                           -1.34412;
        a[10]=                           -0.04290;
        a[11]=                           -0.30359;
        dprb = (b[1] +a[1]*sb + a[2]*sb*sb +a[3]*sb*sb*sb)*pow(x,(b[2] +a[4]*sb + a[5]*sb*sb +a[6]*sb*sb*sb))
        *pow((1.-x),(b[3] +a[7]*sb + a[8]*sb*sb +a[9]*sb*sb*sb))*(1.+(a[10]*sb +a[11]*sb*sb)/x);
        //          dprb= (b[1] +a[1]*sb + a[2]*sb**2 +a[3]*sb**3)*x**
        //      .        (b[2] +a[4]*sb + a[5]*sb**2 +a[6]*sb**3)*(1d0-x)**
        //      .        (b[3] +a[7]*sb + a[8]*sb**2 +a[9]*sb**3)
        //      .        *(1d0+(a[10]*sb +a[11]*sb**2)/x)
    }
    else
    {
        // --- NLO FFs
        if(iset != 1) std::cerr << "ERROR [Fragmentation::kkp]: iset must be 0 (LO) or 1 (NLO)" << std::endl;
        rlam= 0.213;
        s=  log(log(qs*qs/(rlam*rlam))/log(q0*q0/(rlam*rlam)));
        sc= log(log(qs*qs/(rlam*rlam))/log(rmcc*rmcc/(rlam*rlam)));
        sb= log(log(qs*qs/(rlam*rlam))/log(rmbb*rmbb/(rlam*rlam)));
        
        // ---------------------- NLO PION ------------------------------
        b[1]=                            3.73331;
        b[2]=                           -0.74159;
        b[3]=                            2.33092;
        a[1]=                           -3.16946;
        a[2]=                           -0.47683;
        a[3]=                            0.70270;
        a[4]=                           -0.51377;
        a[5]=                           -0.19705;
        a[6]=                           -0.17917;
        a[7]=                            2.03394;
        a[8]=                           -0.50764;
        a[9]=                           -0.08565;
        a[10]=                            0.09466;
        a[11]=                           -0.10222;
        dpg= (b[1] + a[1]*s + a[2]*s*s +a[3]*s*s*s)*pow(x,(b[2] +a[4]*s + a[5]*s*s +a[6]*s*s*s))
        *pow((1.-x),(b[3] +a[7]*s + a[8]*s*s +a[9]*s*s*s))*(1.+(a[10]*s +a[11]*s*s)/x);
        //       dpg= (b[1] +a[1]*s + a[2]*s**2 +a[3]*s**3)*x**
        //      .        (b[2] +a[4]*s + a[5]*s**2 +a[6]*s**3)*(1d0-x)**
        //      .        (b[3] +a[7]*s + a[8]*s**2 +a[9]*s**3)*
        //      .        (1d0+(a[10]*s +a[11]*s**2)/x)
        b[1]=                           0.44809;
        b[2]=                          -1.47598;
        b[3]=                           0.91338;
        a[1]=                           -0.13828;
        a[2]=                           -0.06951;
        a[3]=                            0.01354;
        a[4]=                           -0.30498;
        a[5]=                           -0.01863;
        a[6]=                           -0.12529;
        a[7]=                            0.64145;
        a[8]=                            0.07270;
        a[9]=                           -0.16989;
        a[10]=                            0.07396;
        a[11]=                           -0.07757;
        dpu= (b[1] + a[1]*s + a[2]*s*s +a[3]*s*s*s)*pow(x,(b[2] +a[4]*s + a[5]*s*s +a[6]*s*s*s))
        *pow((1.-x),(b[3] +a[7]*s + a[8]*s*s +a[9]*s*s*s))*(1.+(a[10]*s +a[11]*s*s)/x);
        //       dpu= (b[1] +a[1]*s + a[2]*s**2 +a[3]*s**3)*x**
        //      .        (b[2] +a[4]*s + a[5]*s**2 +a[6]*s**3)*(1d0-x)**
        //      .        (b[3] +a[7]*s + a[8]*s**2 +a[9]*s**3)*
        //      .        (1d0+(a[10]*s +a[11]*s**2)/x)
        b[1]=                            16.5987;
        b[2]=                            0.13345;
        b[3]=                            5.89903;
        a[1]=                           -18.3856;
        a[2]=                            2.44225;
        a[3]=                            2.13225;
        a[4]=                            0.22712;
        a[5]=                           -0.83625;
        a[6]=                            0.38526;
        a[7]=                           -0.16911;
        a[8]=                            0.59886;
        a[9]=                           -0.25630;
        a[10]=                          -0.18619;
        a[11]=                           0.87362;
        dps= (b[1] + a[1]*s + a[2]*s*s +a[3]*s*s*s)*pow(x,(b[2] +a[4]*s + a[5]*s*s +a[6]*s*s*s))
        *pow((1.-x),(b[3] +a[7]*s + a[8]*s*s +a[9]*s*s*s))*(1.+(a[10]*s +a[11]*s*s)/x);
        //          dps= (b[1] +a[1]*s + a[2]*s**2 +a[3]*s**3)*x**
        //      .        (b[2] +a[4]*s + a[5]*s**2 +a[6]*s**3)*(1d0-x)**
        //      .        (b[3] +a[7]*s + a[8]*s**2 +a[9]*s**3)*
        //      .        (1d0+(a[10]*s +a[11]*s**2)/x)
        b[1]=                            6.17173;
        b[2]=                           -0.53618;
        b[3]=                            5.60108;
        a[1]=                           -4.82450;
        a[2]=                           -1.30844;
        a[3]=                            1.95527;
        a[4]=                           -0.27879;
        a[5]=                           -0.51337;
        a[6]=                            0.10900;
        a[7]=                            0.83571;
        a[8]=                           -1.15141;
        a[9]=                            0.77027;
        a[10]=                            0.09268;
        a[11]=                           -0.11267;
        dpc= (b[1] + a[1]*sc + a[2]*sc*sc +a[3]*sc*sc*sc)*pow(x,(b[2] +a[4]*sc + a[5]*sc*sc +a[6]*sc*sc*sc))
        *pow((1.-x),(b[3] +a[7]*sc + a[8]*sc*sc +a[9]*sc*sc*sc))*(1.+(a[10]*sc +a[11]*sc*sc)/x);
        //          dpc= (b[1] +a[1]*sc + a[2]*sc**2 +a[3]*sc**3)*x**
        //      .        (b[2] +a[4]*sc + a[5]*sc**2 +a[6]*sc**3)*(1d0-x)**
        //      .        (b[3] +a[7]*sc + a[8]*sc**2 +a[9]*sc**3)*
        //      .        (1d0+(a[10]*sc +a[11]*sc**2)/x)
        b[1]=                            0.25944;
        b[2]=                           -1.98713;
        b[3]=                            3.52857;
        a[1]=                           -0.11449;
        a[2]=                            0.03733;
        a[3]=                           -0.18028;
        a[4]=                           -0.35858;
        a[5]=                            0.22277;
        a[6]=                           -0.66413;
        a[7]=                            0.72303;
        a[8]=                            0.46260;
        a[9]=                           -0.99235;
        a[10]=                           -0.02701;
        a[11]=                           -0.02089;
        dpb= (b[1] + a[1]*sb + a[2]*sb*sb +a[3]*sb*sb*sb)*pow(x,(b[2] +a[4]*sb + a[5]*sb*sb +a[6]*sb*sb*sb))
        *pow((1.-x),(b[3] +a[7]*sb + a[8]*sb*sb +a[9]*sb*sb*sb))*(1.+(a[10]*sb +a[11]*sb*sb)/x);
        //          dpb= (b[1] +a[1]*sb + a[2]*sb**2 +a[3]*sb**3)*x**
        //      .        (b[2] +a[4]*sb + a[5]*sb**2 +a[6]*sb**3)*(1d0-x)**
        //      .        (b[3] +a[7]*sb + a[8]*sb**2 +a[9]*sb**3)*
        //      .        (1d0+(a[10]*sb +a[11]*sb**2)/x)
        
        // ---------------------- NLO KAON ------------------------------
        b[1]=                            0.23140;
        b[2]=                           -1.36400;
        b[3]=                            1.79761;
        a[1]=                           -0.33644;
        a[2]=                            0.16204;
        a[3]=                           -0.02598;
        a[4]=                            0.97182;
        a[5]=                           -0.02908;
        a[6]=                           -0.43195;
        a[7]=                            1.57116;
        a[8]=                            0.71847;
        a[9]=                           -0.68331;
        a[10]=                            0.36906;
        a[11]=                            2.39060;
        dkg= (b[1] + a[1]*s + a[2]*s*s +a[3]*s*s*s)*pow(x,(b[2] +a[4]*s + a[5]*s*s +a[6]*s*s*s))
        *pow((1.-x),(b[3] +a[7]*s + a[8]*s*s +a[9]*s*s*s))*(1.+(a[10]*s +a[11]*s*s)/x);
        //          dkg= (b[1] +a[1]*s + a[2]*s**2 +a[3]*s**3)*x**
        //      .        (b[2] +a[4]*s + a[5]*s**2 +a[6]*s**3)*(1d0-x)**
        //      .        (b[3] +a[7]*s + a[8]*s**2 +a[9]*s**3)
        //      .        *(1d0+(a[10]*s +a[11]*s**2)/x)
        b[1]=                            0.17806;
        b[2]=                           -0.53733;
        b[3]=                            0.75940;
        a[1]=                           -0.10988;
        a[2]=                           -0.02524;
        a[3]=                            0.03142;
        a[4]=                           -0.60058;
        a[5]=                            0.07863;
        a[6]=                            0.13276;
        a[7]=                            0.61356;
        a[8]=                           -0.43886;
        a[9]=                            0.23942;
        a[10]=                            0.10742;
        a[11]=                            0.12800;
        dku= (b[1] + a[1]*s + a[2]*s*s +a[3]*s*s*s)*pow(x,(b[2] +a[4]*s + a[5]*s*s +a[6]*s*s*s))
        *pow((1.-x),(b[3] +a[7]*s + a[8]*s*s +a[9]*s*s*s))*(1.+(a[10]*s +a[11]*s*s)/x);
        //          dku= (b[1] +a[1]*s + a[2]*s**2 +a[3]*s**3)*x**
        //      .        (b[2] +a[4]*s + a[5]*s**2 +a[6]*s**3)*(1d0-x)**
        //      .        (b[3] +a[7]*s + a[8]*s**2 +a[9]*s**3)
        //      .        *(1d0+(a[10]*s +a[11]*s**2)/x)
        b[1]=                            4.96269;
        b[2]=                            0.05562;
        b[3]=                            2.79926;
        a[1]=                            1.54098;
        a[2]=                           -9.06376;
        a[3]=                            4.94791;
        a[4]=                            1.88660;
        a[5]=                           -2.94350;
        a[6]=                            1.04227;
        a[7]=                            3.02991;
        a[8]=                           -4.14807;
        a[9]=                            1.91494;
        a[10]=                            0.85450;
        a[11]=                           -0.61016;
        dkd= (b[1] + a[1]*s + a[2]*s*s +a[3]*s*s*s)*pow(x,(b[2] +a[4]*s + a[5]*s*s +a[6]*s*s*s))
        *pow((1.-x),(b[3] +a[7]*s + a[8]*s*s +a[9]*s*s*s))*(1.+(a[10]*s +a[11]*s*s)/x);
        //          dkd= (b[1] +a[1]*s + a[2]*s**2 +a[3]*s**3)*x**
        //      .        (b[2] +a[4]*s + a[5]*s**2 +a[6]*s**3)*(1d0-x)**
        //      .        (b[3] +a[7]*s + a[8]*s**2 +a[9]*s**3)
        //      .        *(1d0+(a[10]*s +a[11]*s**2)/x)
        b[1]=                            4.25954;
        b[2]=                           -0.24144;
        b[3]=                            4.21265;
        a[1]=                           -5.44309;
        a[2]=                            6.11031;
        a[3]=                           -3.13973;
        a[4]=                           -1.07757;
        a[5]=                            1.52364;
        a[6]=                           -0.74308;
        a[7]=                            0.25590;
        a[8]=                            0.98423;
        a[9]=                           -0.52839;
        a[10]=                           -0.04000;
        a[11]=                            0.08695;
        dkc= (b[1] + a[1]*sc + a[2]*sc*sc +a[3]*sc*sc*sc)*pow(x,(b[2] +a[4]*sc + a[5]*sc*sc +a[6]*sc*sc*sc))
        *pow((1.-x),(b[3] +a[7]*sc + a[8]*sc*sc +a[9]*sc*sc*sc))*(1.+(a[10]*sc +a[11]*sc*sc)/x);
        //       dkc= (b[1] +a[1]*sc + a[2]*sc**2 +a[3]*sc**3)*x**
        //      .        (b[2] +a[4]*sc + a[5]*sc**2 +a[6]*sc**3)*(1d0-x)**
        //      .        (b[3] +a[7]*sc + a[8]*sc**2 +a[9]*sc**3)
        //      .        *(1d0+(a[10]*sc +a[11]*sc**2)/x)
        b[1]=                            1.32443;
        b[2]=                           -0.88351;
        b[3]=                            6.15221;
        a[1]=                           -1.41156;
        a[2]=                           -0.04809;
        a[3]=                            0.79066;
        a[4]=                           -0.44818;
        a[5]=                           -0.60073;
        a[6]=                            0.45526;
        a[7]=                            0.46679;
        a[8]=                           -0.50792;
        a[9]=                            0.67006;
        a[10]=                           -0.00477;
        a[11]=                           -0.05503;
        dkb= (b[1] + a[1]*sb + a[2]*sb*sb +a[3]*sb*sb*sb)*pow(x,(b[2] +a[4]*sb + a[5]*sb*sb +a[6]*sb*sb*sb))
        *pow((1.-x),(b[3] +a[7]*sb + a[8]*sb*sb +a[9]*sb*sb*sb))*(1.+(a[10]*sb +a[11]*sb*sb)/x);
        //          dkb= (b[1] +a[1]*sb + a[2]*sb**2 +a[3]*sb**3)*x**
        //      .        (b[2] +a[4]*sb + a[5]*sb**2 +a[6]*sb**3)*(1d0-x)**
        //      .        (b[3] +a[7]*sb + a[8]*sb**2 +a[9]*sb**3)
        //      .        *(1d0+(a[10]*sb +a[11]*sb**2)/x)
        
        // ---------------------- NLO PROTON -----------------------------
        b[1]=                            1.56255;
        b[2]=                            0.01567;
        b[3]=                            3.57583;
        a[1]=                           -1.48158;
        a[2]=                           -0.39439;
        a[3]=                            0.51249;
        a[4]=                           -2.16232;
        a[5]=                            2.47127;
        a[6]=                           -0.93259;
        a[7]=                            3.33958;
        a[8]=                           -3.05265;
        a[9]=                            1.21042;
        a[10]=                          -0.84816;
        a[11]=                           1.23583;
        dprg= (b[1] + a[1]*s + a[2]*s*s +a[3]*s*s*s)*pow(x,(b[2] +a[4]*s + a[5]*s*s +a[6]*s*s*s))
        *pow((1.-x),(b[3] +a[7]*s + a[8]*s*s +a[9]*s*s*s))*(1.+(a[10]*s +a[11]*s*s)/x);
        //          dprg= (b[1] +a[1]*s + a[2]*s**2 +a[3]*s**3)*x**
        //      .        (b[2] +a[4]*s + a[5]*s**2 +a[6]*s**3)*(1d0-x)**
        //      .        (b[3] +a[7]*s + a[8]*s**2 +a[9]*s**3)
        //      .        *(1d0+(a[10]*s +a[11]*s**2)/x)
        b[1]=                           1.25946;
        b[2]=                           0.07124;
        b[3]=                           4.12795;
        a[1]=                           -1.17505;
        a[2]=                            0.37550;
        a[3]=                           -0.01416;
        a[4]=                           -0.29533;
        a[5]=                           -0.24540;
        a[6]=                            0.16543;
        a[7]=                            0.98867;
        a[8]=                           -0.46846;
        a[9]=                            0.20750;
        a[10]=                            0.18957;
        a[11]=                           -0.01116;
        dpru= (b[1] + a[1]*s + a[2]*s*s +a[3]*s*s*s)*pow(x,(b[2] +a[4]*s + a[5]*s*s +a[6]*s*s*s))
        *pow((1.-x),(b[3] +a[7]*s + a[8]*s*s +a[9]*s*s*s))*(1.+(a[10]*s +a[11]*s*s)/x);
        //          dpru= (b[1] +a[1]*s + a[2]*s**2 +a[3]*s**3)*x**
        //      .        (b[2] +a[4]*s + a[5]*s**2 +a[6]*s**3)*(1d0-x)**
        //      .        (b[3] +a[7]*s + a[8]*s**2 +a[9]*s**3)
        //      .        *(1d0+(a[10]*s +a[11]*s**2)/x)
        b[1]=                            4.01135;
        b[2]=                            0.17258;
        b[3]=                            5.20766;
        a[1]=                            8.67124;
        a[2]=                           -22.7888;
        a[3]=                            11.4720;
        a[4]=                            4.57608;
        a[5]=                           -9.64835;
        a[6]=                            4.61792;
        a[7]=                            7.25144;
        a[8]=                           -12.6313;
        a[9]=                            6.07314;
        a[10]=                            0.16931;
        a[11]=                           -0.09541;
        dprs= (b[1] + a[1]*s + a[2]*s*s +a[3]*s*s*s)*pow(x,(b[2] +a[4]*s + a[5]*s*s +a[6]*s*s*s))
        *pow((1.-x),(b[3] +a[7]*s + a[8]*s*s +a[9]*s*s*s))*(1.+(a[10]*s +a[11]*s*s)/x);
        //          dprs= (b[1] +a[1]*s + a[2]*s**2 +a[3]*s**3)*x**
        //      .        (b[2] +a[4]*s + a[5]*s**2 +a[6]*s**3)*(1d0-x)**
        //      .        (b[3] +a[7]*s + a[8]*s**2 +a[9]*s**3)
        //      .        *(1d0+(a[10]*s +a[11]*s**2)/x)
        b[1]=                            0.08250;
        b[2]=                           -1.61290;
        b[3]=                            2.01255;
        a[1]=                           -0.04512;
        a[2]=                           -0.00565;
        a[3]=                            0.00900;
        a[4]=                           -0.38012;
        a[5]=                           -0.06840;
        a[6]=                            0.08888;
        a[7]=                            0.63782;
        a[8]=                           -0.14146;
        a[9]=                            0.06083;
        a[10]=                           -0.02958;
        a[11]=                            0.01130;
        dprc= (b[1] + a[1]*sc + a[2]*sc*sc +a[3]*sc*sc*sc)*pow(x,(b[2] +a[4]*sc + a[5]*sc*sc +a[6]*sc*sc*sc))
        *pow((1.-x),(b[3] +a[7]*sc + a[8]*sc*sc +a[9]*sc*sc*sc))*(1.+(a[10]*sc +a[11]*sc*sc)/x);
        //          dprc= (b[1] +a[1]*sc + a[2]*sc**2 +a[3]*sc**3)*x**
        //      .        (b[2] +a[4]*sc + a[5]*sc**2 +a[6]*sc**3)*(1d0-x)**
        //      .        (b[3] +a[7]*sc + a[8]*sc**2 +a[9]*sc**3)
        //      .        *(1d0+(a[10]*sc +a[11]*sc**2)/x)
        b[1]=                            24.2916;
        b[2]=                            0.57939;
        b[3]=                            12.1207;
        a[1]=                          -88.3524;
        a[2]=                            93.1056;
        a[3]=                           -17.4089;
        a[4]=                           -0.80783;
        a[5]=                           -5.07200;
        a[6]=                           -2.45377;
        a[7]=                           -3.27370;
        a[8]=                            1.21188;
        a[9]=                           -5.50374;
        a[10]=                            0.14628;
        a[11]=                           -0.78634;
        dprb= (b[1] + a[1]*sb + a[2]*sb*sb +a[3]*sb*sb*sb)*pow(x,(b[2] +a[4]*sb + a[5]*sb*sb +a[6]*sb*sb*sb))
        *pow((1.-x),(b[3] +a[7]*sb + a[8]*sb*sb +a[9]*sb*sb*sb))*(1.+(a[10]*sb +a[11]*sb*sb)/x);
        //          dprb= (b[1] +a[1]*sb + a[2]*sb**2 +a[3]*sb**3)*x**
        //      .        (b[2] +a[4]*sb + a[5]*sb**2 +a[6]*sb**3)*(1d0-x)**
        //      .        (b[3] +a[7]*sb + a[8]*sb**2 +a[9]*sb**3)
        //      .        *(1d0+(a[10]*sb +a[11]*sb**2)/x)
    }
    
    // --- Evaluate different contributions
    
    double dpd =      dpu;
    dks =      dku;
    double dprd= 0.5*dpru;
    if (qs < rmbb)
    {
        dpb  = 0.;
        dkb  = 0.;
        dprb = 0.;
    }
    if (qs < rmcc)
    {
        dpc  = 0.;
        dkc  = 0.;
        dprc = 0.;
    }
    if (ih == 1)
    {
        dh[0]= dpg/2.;
        dh[1]= dpu/2.;
        dh[2]= dpu/2.;
        dh[3]= dpd/2.;
        dh[4]= dpd/2.;
        dh[5]= dps/2.;
        dh[6]= dps/2.;
        dh[7]= dpc/2.;
        dh[8]= dpc/2.;
        dh[9]= dpb/2.;
        dh[10]=dpb/2.;
    }
    else if (ih == 2)
    {
        dh[0]= dkg/2.;
        dh[1]= dku/2.;
        dh[2]= dku/2.;
        dh[3]= dkd/2.;
        dh[4]= dkd/2.;
        dh[5]= dks/2.;
        dh[6]= dks/2.;
        dh[7]= dkc/2.;
        dh[8]= dkc/2.;
        dh[9]= dkb/2.;
        dh[10]=dkb/2.;
    }
    else if (ih == 3)
    {
        dh[0]= dkg/2.;
        dh[1]= dkd/2.;
        dh[2]= dkd/2.;
        dh[3]= dku/2.;
        dh[4]= dku/2.;
        dh[5]= dks/2.;
        dh[6]= dks/2.;
        dh[7]= dkc/2.;
        dh[8]= dkc/2.;
        dh[9]= dkb/2.;
        dh[10]=dkb/2.;
    }
    else if (ih == 4)
    {
        dh[0]= dprg/2.;
        dh[1]= dpru/2.;
        dh[2]= dpru/2.;
        dh[3]= dprd/2.;
        dh[4]= dprd/2.;
        dh[5]= dprs/2.;
        dh[6]= dprs/2.;
        dh[7]= dprc/2.;
        dh[8]= dprc/2.;
        dh[9]= dprb/2.;
        dh[10]=dprb/2.;
    }
    else if (ih == 5)
    {
        dh[0]= dpg/2.;
        dh[1]= dpu/2.;
        dh[2]= dpu/2.;
        dh[3]= dpd/2.;
        dh[4]= dpd/2.;
        dh[5]= dps/2.;
        dh[6]= dps/2.;
        dh[7]= dpc/2.;
        dh[8]= dpc/2.;
        dh[9]= dpb/2.;
        dh[10]=dpb/2.;
    }
    else if (ih == 6)
    {
        dh[0]= dprg/2.;
        dh[1]= dpru/4.;
        dh[2]= dpru/4.;
        dh[3]= dprd;
        dh[4]= dprd;
        dh[5]= dprs/2.;
        dh[6]= dprs/2.;
        dh[7]= dprc/2.;
        dh[8]= dprc/2.;
        dh[9]= dprb/2.;
        dh[10]=dprb/2.;
    }      
    else
    {  
        dh[0]= dpg+dkg+dprg;
        dh[1]= dpu+dku+dpru;
        dh[2]= dpu+dku+dpru;
        dh[3]= dpd+dkd+dprd;
        dh[4]= dpd+dkd+dprd;
        dh[5]= dps+dks+dprs;
        dh[6]= dps+dks+dprs;
        dh[7]= dpc+dkc+dprc;
        dh[8]= dpc+dkc+dprc;
        dh[9]= dpb+dkb+dprb;
        dh[10]=dpb+dkb+dprb;
    }  
    
  //     Parton label:
  //     0    1    2    3    4    5    6    7    8     9    10
  //     g    u   ubar  d   dbar  s   sbar  c   cbar   b   bbar
  //    
  switch ( species )
  {
    case gluon: return dh[0]; // return gluon part
    break;
    case up: return dh[1]; // return up part
    break;
    case anti_up: return dh[2]; // return anti-up part
    break;
    case down: return dh[3]; // return down part
    break;
    case anti_down: return dh[4]; // return anti-down part
    break;
    case strange: return dh[5]; // return strange part
    break;    
    case anti_strange: return dh[6]; // return anti-strange part
    break;
    default:
     return 0.; 
     break;
  } 
}


///////////////////
// FRAGMENTATION //
///////////////////

// void analysis::PerformFragmentation(int Nz,double zMin,double zMax,int HadronSpecies,InterpolatedSpectrum *GluonSpectrum,InterpolatedSpectrum *HadronSpectrum)
// {
//     
//     // PERFORM CONVOLUTION WITH FRAGMENTATION FUNCTION //
//     for(int ikg=0;ikg<GluonSpectrum->Nk;ikg++){
//         for(int iphig=0;iphig<GluonSpectrum->NPhi;iphig++){
//             
//             for(int iz=0;iz<Nz;iz++){
//                 
//                 // GET dz AND z //
//                 double dz=(zMax-zMin)/DOUBLE(Nz);   double z=(iz+0.5)*dz+zMin;
//                 
//                 if(z<zMin || z>zMax){
//                     std::cerr << "# ERROR" << std::endl;
//                 }
//                 
//                 // GET MOMENTUM kT AND ANGLE phi OF THE GLUON //
//                 double kg=GluonSpectrum->GetK(ikg);    double phig=GluonSpectrum->GetPhi(iphig);
//                 
//                 // COMPUTE g->h FRAGMENTATION FUNCTIONS //
//                 double Dhz=KKPFragmentation(HadronSpecies,1,z,kg);
//                 
//                 // COMPUTE N_h -- NOTE THAT 1/z^2 FACTOR IS CANCELED BY INTEGRATION OVER pT BIN //
//                 double Nh=Dhz*GluonSpectrum->Values[ikg][iphig]*dz;
//                 
//                 if(Dhz<0){
//                     std::cerr << "## WARNING -- NEGATIVE FRAGMENTATION FUNCTION Dhz=" << Dhz << " at z=" << z << std::endl;
//                 }
//                 
//                 // GET MOMENTUM pT AND ANGLE phi OF THE HADRON //
//                 double kh=z*kg; double phih=phig;
//                 
//                 // GET ASSOCIATED MOMENTUM BIN //
//                 int ikh=HadronSpectrum->GetKBin(kh);    int iphih=HadronSpectrum->GetPhiBin(phih);
//                 
//                 if(iphih!=iphig){
//                     
//                     std::cerr << "ERROR IN FRAGMENTATION" << std::endl;
//                     
//                     exit(0);
//                 }
//                 
//                 // INCREASE HADRON SPECTRUM //
//                 if(ikh>=0 && ikh<=HadronSpectrum->Nk-1 && iphih>=0 && iphih<=HadronSpectrum->NPhi-1){
//                     HadronSpectrum->Values[ikh][iphih]+=Nh; HadronSpectrum->dNOverdY+=Nh;
//                 }
//                 
//             }
//         }
//         
//     }
//     
//     // COMMANDLINE OUTPUT //
//     std::cerr << "#CHARGED HADRON dN/dY= " << HadronSpectrum->dNOverdY << std::endl;
//     
// }

