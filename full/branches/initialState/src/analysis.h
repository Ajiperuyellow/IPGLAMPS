//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/analysis.h $
//$LastChangedDate: 2019-03-21 17:38:54 +0100 (Do, 21. Mär 2019) $
//$LastChangedRevision: 2966 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <fstream>
#include <stdint.h>
#include <string>
#include <time.h>
#include <vector>
#include <math.h>

#include "analysis_hydro.h"
#include "bampsvector.h"
#include "configuration.h"
#include "FPT_compare.h"
#include "heavyIonCollision.h"
#include "cellcontainer.h"

using std::vector;
using std::fstream;

class heavyIonCollision;
class config;

enum jetTrackerCollType {initial_jet, final_jet, c2to2, c2to3, c3to2, production};
enum anaType {ptSpectrum, ptSpectrumSoft, all, initial, final, jets, v2output, meanv2output,CharmNumber,BottomNumber,hydroNormal,hydroArrow,hydroMidRap,hydroMidRapArrow,hydroParticleOutput,midRapObservables, pos_mom};
enum partonCombinationType { glue_glue, quark_antiquark, glue_quark, quark_quark };
enum doubleAntennaCombinationType { gluonAntenna_gluonAntenna, gluonAntenna_quarkAntenna, quarkAntenna_quarkAntenna};


const string sep = "\t";
class cluster
{
public:
  cluster():
  centerofmass(0,0,0,0),
  centerofmassMom(0,0,0,0),
  S(0),
  invMass(0),
  volume(0),
  indexBiggestCluster(0),
  myColor(0)
  {};
  ~cluster() {};

  std::vector<int> particleList;
  
  void calcInvMass();
  void calcVolume();
  void getDistToCenter();
  void calcCOM();
  
  VectorTXYZ centerofmass;
  VectorEPxPyPz centerofmassMom;
  VectorEPxPyPz centerofmassVel;
  double S;
  double invMass;
  double volume; // maybe one-dimensional
  int indexBiggestCluster;
  std::vector<double> distToCenterOfMass;
  int myColor;
  
  int MotherA,MotherB;
    
  int streamID1,streamID2;
  
private:
 


};


class analysisRapidityRange
{
  public:
    analysisRapidityRange() : yleft(0), yright(0) {};
    analysisRapidityRange( const double _left, const double _right ) : yleft(_left), yright(_right) {};
  
    void reset( const double _left, const double _right ) { yleft = _left; yright = _right; }
    double ydelta() const { return (yright - yleft); }
    
    double yleft;
    double yright;
};

class jetTrackerSingleEvent
{
public:
  jetTrackerSingleEvent() : 
    R_proj(), 
    P_proj_in(),
    P_proj_out(), 
    P1_in(), P2_in(), P1_out(), P2_out(),
    lambda(-1), 
    xSection(-1), 
    cell_ID(-1)
  { };

  ~jetTrackerSingleEvent() {};
  
  VectorTXYZ R_proj;

  VectorEPxPyPz P_proj_in;
  VectorEPxPyPz P_proj_out;

  VectorEPxPyPz P1_in;
  VectorEPxPyPz P2_in;
  VectorEPxPyPz P1_out;
  VectorEPxPyPz P2_out;

  jetTrackerCollType coll_type;
  double lambda;
  double xSection;
  int cell_ID;
  int jet_ID_in, jet_ID_out;
};

class analysis
{
 public:
  analysis( config * const theConfig, heavyIonCollision * const theHIC );
  ~analysis();

  void generalOutput( const int step, const double timeshift );
  void generalOutput( const anaType at, const double timeshift );
  void output();
  
  void collectPtData( const int step );
  void collectPtDataInitial();

  void addBoost32( const double beta );
  void addBoost23( const double beta );

  void setSeed( uint32_t _s ) { seed = _s; }
  uint32_t getSeed() const { return seed; }

  int addJetEvent_in( const int ID_1, const int ID_2, const int ID_3, 
                      const jetTrackerCollType coll_type,
                      const double cross_section, const int cell_ID, 
                      const double lambda );
  void addJetEvent_initial( const int jetID );
  void addJetEvents_final();
  void addJetEvent_out( const int entity_ID, 
                        const int ID_1, const int ID_2, const int ID_3 );
  void removeJetEvent_in( const int entity_ID );
  void exchangeJetID( const int oldID, const int newID );
  void makeJetTrackerCopy();
  void restoreJetTracker();

  void ParticleCorrelations(double time, FLAVOR_TYPE thisFlavor);
  void ChargedHadronCorrelations(double time);
  void FragmentedCorrelations(double time);
  
  void Energydensityplot(double time);
  void PlotCellwise(double time, int number, int IX, int IY, int centralEtaIndex,const vector<cellContainer> &cells );
  
  double getJetTracking_PT() const {return jetTracking_PT;}

  void registerProgressInfoForOutput( 
           const double _time, const double _dt, const int _nTotalParticles, const int _nParticlesInFormationIntialProduction, const int _nParticlesInFormationGeometricCollisions, const int _nParticlesActive, const int _nParticlesInActiveCells, const vector< int >& _edgeCellSizes, const int _nFreeParticles, const int _nColl, const int _nColl22, const int _nColl23, const int _nColl32, const int _nColle );
  
  void writeParticleAndCollisionNumbers(
           const double _time, 
           const int _nTotalParticles, 
           const int ncoll22, const int ncoll23, const int ncoll32, 
           const int nGet32Errors, const int nGet23Errors, 
           const int _nParticlesInFormationIntialProduction, 
           const int _nParticlesInFormationGeometricCollisions, 
           const int _nParticlesActive, 
           const int _nParticlesInActiveCells, 
           const std::vector< int >& _edgeCellSizes, 
           const int _nFreeParticles,
           const double sin2Theta,
           const double vrelSum, 
           const double sin2ThetaVrel,
           const double ncoll22MIDRAP
                                       );

  void writeParticleList(const int step);
  
  
  double tstep[100000];

  string infoInitialParameters;
  
  //number of collisions
  int ncoll22, ncoll23, ncoll32, ncoll22MIDRAP;  
  double sin2theta,vrelSum,sin2ThetaVrel;

  vector<cluster> theClusters;
  void Cluster();
  void printClusters(vector<cluster> allClusters, int no);
  void writeAllClustersInLesHouches(vector<cluster> allClusters, int no);
  void runHERWIG();
  void generateLesHouchesEvent(cluster oneCluster, int no);
  bool noMoreSingletGluons(vector<vector<int>> ParticleListWColors);
  void GiveNextColors(int& color, int& anticolor, int& recentColor, int & recentAnticolor, int whichParticle, int & countColors);
  void GenerateDCATimeOrderedColorsLesHouchesEvent(cluster oneCluster, int no);
  void distributeEqualQuarks(cluster & oneCluster);
  void GenerateShuffledColorsLesHouchesEvent(cluster oneCluster, int no);
  void writeMomenta(VectorEPxPyPz aVector, string & thestring);
  int NumberCluster=0;
  std::vector<int> ClusterNumbers{};
  void getBetaToCoM(cluster ACluster, VectorEPxPyPz & beta);
  void readHadrons();
  vector<string> split(string const &input); 
  bool isUncharged(string particleID);
  bool isCharged(string particleID);
  partonCombinationType whichQCDCombination(int i,int j);
  void quark_quark_handler(int i,int j);
  void octet_handler(int A1, int A2, int B1, int B2, int & colorA1, int & anticolorA1, int & colorA2, int & anticolorA2, int & colorB1, int & anticolorB1, int & colorB2, int & anticolorB2);
  
private:
  //------------------------
  hydroParticleType theHydroParticleType;
  //------------------------

  void ptDistribution( const FLAVOR_TYPE _flavTypeToComputeFor, 
                       vector< Particle >& _particles, 
                       const int n_particles, 
                       const int step );
  void ptSoftDistribution( const FLAVOR_TYPE _flavTypeToComputeFor, 
                           vector< Particle >& _particles, 
                           const int n_particles, 
                           const int step );
  void particleOutput( const int step );
  void ParticlePosMomOutput( const int step, const double timeshift);
  void finalOutput();
  void v2_output();
  void jetTrackerOutput();
  void charmNumber(const int step);
  void bottomNumber(int arg1);
  void print_dndy( const int step );

  void printPtSpectra( const FLAVOR_TYPE _flavTypeToComputeFor, 
                       const bool printIntemediate );
  void printSoftPtSpectra( const FLAVOR_TYPE _flavTypeToComputeFor, 
                           const bool printIntemediate );
  void printHeader( fstream & f, const anaType mode, const time_t end );
  
  //------------------------------//
  void anaTimeStep_generalHIC();
  void anaTimeStep_heavyQuark();
  void anaTimeStep_hydroParametrization();
  void anaTimeStep_detailed();
  void anaTimeStep_correlations();
  void anaTimeStepVideo();
  //------------------------------//  

  //---------------------------------  
  config * const theConfig ;
  heavyIonCollision * const theHIC;
  anaHydro * theAnaHydro;
  //---------------------------------
  
  void computeV2RAA( string, const double _outputTime );
  
  string filename_prefix;
  
  std::vector<analysisRapidityRange> rapidityRanges;

  time_t start;
  
  OUTPUT_SCHEME outputScheme;
  
  void handle_output_studies( OUTPUT_SCHEME _outputScheme );
  
  bool studyParticleSpectra;
  bool studyParticleSpectraAllSteps;
  bool studyHeavyQuarkProduction;
  bool studyDetailedParticleOutput;
  bool studyDetailedParticleOutputAllSteps;
  bool studyHydro;
  bool studyDetailedParticlePosMomOutput;
  bool studyCorrelations;
  bool studyCellwiseVideo;
  bool noAnalysis;
  
  bool studyV2;
  bool studyBoostDistribution;
  bool studyJetTracking;
  bool studyParticleAndCollisionNumbers;
  bool studyDndy;

  vector<int> boostDistribution32, boostDistribution23;
  int numberBinsBoost;
  
  std::vector<double> *ptBinsDY1_gluons, *ptBinsDY16_gluons, *ptBinsDY2_gluons, *ptBinsDY3_gluons, *ptBinsDY4_gluons, *ptBinsDY5_gluons, *ptBinsAll_gluons;
  std::vector<double> *ptBinsDY1_quarks, *ptBinsDY16_quarks, *ptBinsDY2_quarks, *ptBinsDY3_quarks, *ptBinsDY4_quarks, *ptBinsDY5_quarks, *ptBinsAll_quarks;
  std::vector<double> *ptBinsDY1_ups, *ptBinsDY16_ups, *ptBinsDY2_ups, *ptBinsDY3_ups, *ptBinsDY4_ups, *ptBinsDY5_ups, *ptBinsAll_ups;
  std::vector<double> *ptBinsDY1_downs, *ptBinsDY16_downs, *ptBinsDY2_downs, *ptBinsDY3_downs, *ptBinsDY4_downs, *ptBinsDY5_downs, *ptBinsAll_downs;
  std::vector<double> *ptBinsDY1_stranges, *ptBinsDY16_stranges, *ptBinsDY2_stranges, *ptBinsDY3_stranges, *ptBinsDY4_stranges, *ptBinsDY5_stranges, *ptBinsAll_stranges;
  std::vector<double> *ptBinsDY1_anti_ups, *ptBinsDY16_anti_ups, *ptBinsDY2_anti_ups, *ptBinsDY3_anti_ups, *ptBinsDY4_anti_ups, *ptBinsDY5_anti_ups, *ptBinsAll_anti_ups;
  std::vector<double> *ptBinsDY1_anti_downs, *ptBinsDY16_anti_downs, *ptBinsDY2_anti_downs, *ptBinsDY3_anti_downs, *ptBinsDY4_anti_downs, *ptBinsDY5_anti_downs, *ptBinsAll_anti_downs;
  std::vector<double> *ptBinsDY1_anti_stranges, *ptBinsDY16_anti_stranges, *ptBinsDY2_anti_stranges, *ptBinsDY3_anti_stranges, *ptBinsDY4_anti_stranges, *ptBinsDY5_anti_stranges, *ptBinsAll_anti_stranges;
  std::vector<double> *ptBinsDY1_all, *ptBinsDY16_all, *ptBinsDY2_all, *ptBinsDY3_all, *ptBinsDY4_all, *ptBinsDY5_all, *ptBinsAll_all;
  std::vector<double> ptBinLabels;
  double minPT, maxPT, binWidthPT;
  double maxPTSoft;
  int numberBinsPT;
  std::vector<double> *ptBinsSoftAll_gluons, *ptBinsSoftAll_quarks, *ptBinsSoftAll_all;
  std::vector<double> *ptBinsSoftAll_ups, *ptBinsSoftAll_downs, *ptBinsSoftAll_stranges, *ptBinsSoftAll_anti_ups, *ptBinsSoftAll_anti_downs, *ptBinsSoftAll_anti_stranges;
  std::vector<double> ptSoftBinLabels;
  double binWidthSoftPT;
  int numberBinsSoftPT;
  
  int nRapidityBins;
  vector<double> rapidityBins;
  double min_v2_pt, max_v2_pt, binWidth_v2_pt;
  int nBins_v2_pt;
  vector< vector<double> > *differential_v2_values;
  vector< vector<double> > *differential_v2_counts;
  vector<double> differential_v2_bins_labels;
  vector<double> *mean_PT;
  vector<double> *integrated_v2;
  vector<double> *integrated_v2_counts;
 
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------    

  double KKPFragmentation(int ih, int iset, double x, double qs, FLAVOR_TYPE species);
//   void   PerformFragmentation(int Nz,double zMin,double zMax,int HadronSpecies,InterpolatedSpectrum *GluonSpectrum,InterpolatedSpectrum *HadronSpectrum);
  
  
  
  
  //CLUSTERING
  double getTotalAction(vector<cluster> allClusters);
  
  void ClusterSplitStep(double beta, vector<cluster> &theClusters, bool isWarmUp);
  
  
  cluster cluster_merge(cluster ClusterOne, cluster ClusterTwo); 
  
  void ClusterMergeStep(double beta, vector<cluster> &theClusters, bool isWarmUp);
  int NumberOfMergeAttempts;
  int NumberOfAcceptedMerges;
  int NumberOfSplitAttempts,NumberOfAcceptedSplits;
  double sumOfDeltaS;
  vector<int> Colors;
  int GetNextAvailableColor();
  void Split(cluster Cluster, cluster &NewClusterOne, cluster &NewClusterTwo,int NCluster);
  void SplitRandom(cluster Cluster, cluster &NewClusterOne, cluster &NewClusterTwo,int NCluster);
  //--------------------//
  //for hydro analysis
  void printHeaderForHydro( fstream & f, const anaType mode, const time_t end ); 
  void hydroDistributionOutput();
  void hydroDistributionArrowOutput(); 
  void particleListOutput();
  void hydroDistributionMidRapOutput();      
  void hydroDistributionMidRapArrowOutput();
  void hydroMidRapidityObservablesOutput();   
  int getRealNumberOfTimeStepsForAnaHydro();    
  //--------------------//      

  uint32_t seed;

  vector< vector<jetTrackerSingleEvent> > jetTracker;
  vector< vector<jetTrackerSingleEvent> > jetTracker_copy;

  double jetTracking_PT;

  int nTimeSteps;

  fstream progressLogFile;
  fstream printParticleAndCollisionNumbers;
  
  int *charmNmb;
  int *charmNmb_eta20;
  int *charmNmb_eta15;
  int *charmNmb_eta10;
  int *charmNmb_eta075;
  int *charmNmb_eta05;
  int *charmNmb_eta035;
  int *charmNmb_strap05;
  
  int *bottomNmb;
  int *bottomNmb_eta20;
  int *bottomNmb_eta15;
  int *bottomNmb_eta10;
  int *bottomNmb_eta075;
  int *bottomNmb_eta05;
  int *bottomNmb_eta035;
  
  //---------------------//
  double getAngle2pi(double x, double y) const
  {
    double a = atan2(y,x);
    if ((a*180./M_PI)<0)
    {  
      return ((360.-abs(a*180./M_PI)))*M_PI/180.;
    }
    else
    {
      return a;
    }
  }
  
};

class v2RAA
{
public:
  v2RAA( config * const c, string name_arg, string filename_prefix_arg, std::vector<analysisRapidityRange> rapidityRanges_arg, const double pt_min_arg = 0.0, const double pt_max_arg = 10.0, const int n_g_arg = 25);
  
  void setPtBinProperties( const double pt_min_arg, const double pt_max_arg, const int n_g_arg )
  {
    pt_min = pt_min_arg;
    pt_max = pt_max_arg;
    n_g = n_g_arg;
  }
  void setPtBinsLinear( const double pt_min_arg, const double pt_max_arg, const int n_g_arg )
  {
    pt_min_lin = pt_min_arg;
    pt_max_lin = pt_max_arg;
    n_g_lin = n_g_arg;
  }
  void computeFor( const FLAVOR_TYPE _flavTypeToComputeFor, std::vector< Particle >& _particles, const int n_particles, string additionalNameTag, const double _outputTime );

private:
  
  config * const theConfig ;
  
  double pt_min,pt_max;  
  double pt_min_lin,pt_max_lin; 
  int n_c,n_b,n_g,eta_bins,n_g_lin;
 
  string name,filename_prefix;
  
  std::vector<analysisRapidityRange> rapidityRanges;
  

};


/** @brief exception class for handling unexpected critical behaviour within simulations of heavy ion collisions  */
class eAnalysis_error : public std::runtime_error
{
 public:
  explicit eAnalysis_error(const std::string& what) : std::runtime_error(what) {};
    
  virtual ~eAnalysis_error() throw() {};
};



#endif
