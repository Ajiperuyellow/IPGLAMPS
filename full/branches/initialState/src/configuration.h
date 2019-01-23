//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/configuration.h $
//$LastChangedDate: 2018-12-25 21:53:51 +0100 (Di, 25. Dez 2018) $
//$LastChangedRevision: 2914 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


/** @file
 * @brief This file provides global objects and a configuration interface
 */


#ifndef CONFIGURATION_H
#define CONFIGURATION_H

#include <string>
#include <vector>
#include <stdexcept>
#include <iostream>

// #define BOOST_FILESYSTEM_VERSION 3
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "configurationbase.h"
#include "particle.h"


/** @brief Enumeration type for possible initial state models */
enum INITIAL_STATE_TYPE { miniJetsInitialState, pythiaOfflineInitialState, cgcInitialState, hydroParametrizationInitialState, pythiaOnlineInitialState, simpleInitialState, CYMInitialState, AttractorInitialState };

/** @brief Enumeration type for PDF sources */
enum PDF_SOURCE_TYPE { builtInGRV, LHAPDF };

/** @brief Enumeration type for different variants of computing the mean free path of high-pt particles */
enum JET_MFP_COMPUTATION_TYPE { computeMfpDefault, computeMfpIteration, computeMfpInterpolation };

/** @brief Enumeration type for different time steps modes for the analysis */
enum ANATIMESTEP_MODE { generalHICAnaTimeSteps, heavyQuarkAnaTimeSteps, hydroParametrizationAnaTimeSteps, detailedHICAnaTimeSteps, correlationAnaTimeSteps };



/** @brief Enumeration type for different output schemes to decide which kind of output is printed 
 * set output studies according to the output scheme here in analysis::handle_output_studies( OUTPUT_SCHEME _outputScheme )
*/
enum OUTPUT_SCHEME { no_output = 0, 
allOutputs = 1,
testrun_observables = 10,
correlations = 2,
onlyVideo = 4,
nothingToAnalyze = 3,
// RHIC
phenix_v2 = 21,

// heavy quarks:
heavy_quark_production = 100,

doMovie = 200, // does not work yet

correlationsearch = 500 // detailed timesteps and position momenta output off all active particles
};

/** @brief This extends the namespace ns_casc (see globalsettings.h) to hold a global particle vector */
namespace ns_casc
{  
  /** @brief Vector to hold all particles in the simulation, GLOBAL 
  *
  * Declared as extern here. There MUST be ONE definition in a cpp-file that includes this header in the global 
  * scope (it's in configuration.cpp).
  * All files that include this header can then access the vector.
  */
  extern std::vector<Particle> particles;
  
  extern std::vector<Particle> hadrons;
}
//--------------------------------------------------------//

using std::string;


/**
 * @brief interface for runtime configuration
 *
 * This class provides a basic interface for settings made by the user at runtime. It is based on boost::program_options 
 * and can handle command line arguments as well as INI-type configuration files.
 * 
 * It derives from and extends the base class configBase that is provided in the BAMPS lib. configBase provides the basic
 * functionality, some helper routines (e.g. printing of used parameters) and also some basic boost::program_options::options_description groups
 * and options. The structure of configBase allows this derived class to extend the functionality by 
 *  a) adding to existing boost::program_options::options_description objects
 *  b) defining new boost::program_options::options_description objects and adding those to either the command line or the 
 *    configuration file option groups
 * 
 * This is done in config::initializeProgramOptions() and in config::groupProgramOptions(). In principle these two steps 
 * would be sufficient as they eventually only add to the groups configBase::command_line_options and configBase::config_file_options
 * which are members of the base class. Therefore reading and processing - including the extended options - could then be
 * handled by the base class.
 * However, in order to provide the flexibility to add more sophisticated checking, processing, etc. to the extended options
 * this class explicitly re-implements the processing routines of the base class configBase. Mostly by simply calling 
 * the appropriate member routines of the base class. 
 * 
 * The program flow is therefore identical to the one descriped in the documentation of configBase
 */
class config : public configBase
{
 public:
  /** @brief Constructor that internally reads the provided input file */
  config();
  /** @brief Standard constructor */
  ~config() {};
  
  /** @brief processes command line arguments and settings provided via a configuration file */
  void readAndProcessProgramOptions(const int argc, char* argv[]);

  /** ---- collision parameters ---- */ 
  /** @brief Interface for config::A */
  double getA() const { return A; }
  
  /** @brief Interface for config::Aatomic */
  double getAatomic() const { return Aatomic; }
  
  /** @brief Interface for config::B */
  double getB() const { return B; }
  
  /** @brief Interface for config::Batomic */
  double getBatomic() const { return Batomic; }
  
  /** @brief Interface for config::sqrtS */
  double getSqrtS() const { return sqrtS; }
  
  /** @brief Interface for config::impact */
  double getImpactParameter() const { return impact; }  
  /** ------------------------------- */
  
  /** ---- simulation parameters ---- */ 
  /** @brief Interface for config::freezeOutEnergyDensity */
  double getFreezeOutEnergyDensity() const { return freezeOutEnergyDensity; }

  
  /** @brief Interface for initial state options. */
  bool isSet_transCellSize() const { return isSetTransCellSize; }
  double getTransCellSize() const { return transCellSize; }  
  double getTransLength() const { return transLength; }  
  double getRequiredNinCell() const { return requiredNinCell; }
  double getMinimumEtaInitial() const {return minimumEtaInitial;}
  double getMaximumEtaInitial() const {return maximumEtaInitial;}
  int getNEtaBinsInitial() const {return nEtaBinsInitial;}
  bool doRandomisationAzimuth()  const {return randomizeAzimuth;}
  bool doHomogenSpace() const {return homogenSpace;}
  double getTau0() const {return tau0;}
  int getInitialAttractorCondition() const {return InitialAttractorCondition;}
  double getPTRefMax() const {return pTRefMax;}
  double getPTRefMin() const {return pTRefMin;}
  double getPTHRefMax() const {return pTHRefMax;}
  double getPTHRefMin() const {return pTHRefMin;} 
  /** ------------------------------- */

  /** ---- initial state options ---- */ 
  /** @brief Interface for config::initialStateType */
  INITIAL_STATE_TYPE getInitialStateType() const { return initialStateType; }
  
  /** @brief Interface for config::PDFsource */
  PDF_SOURCE_TYPE getPDFsource() const { return PDFsource; }
  
  /** @brief Interface for config::LHAPDFdatasetName */
  string getLHAPDFdatasetName() const { return LHAPDFdatasetName; }
  
  /** @brief Interface for config::LHAPDFmember */
  unsigned short int getLHAPDFmember() const { return LHAPDFmember; }
  
  /** @brief Interface for config::LHAPDFuseGrid */
  bool getLHAPDFuseGrid() const { return LHAPDFuseGrid; }
  
  /** @brief Interface for config::nuclearPDFs */
  bool useNuclearPDFs() const { return nuclearPDFs; }
  
  /** @brief Interface for config::nuclearPDFdatasetName */
  string getNuclearPDFdatasetName() const { return nuclearPDFdatasetName; }
  
  /** @brief Interface for config::pythiaParticleFile */
  string getPythiaParticleFile() const { return pythiaParticleFile; }

  /** @brief Interface for config::cgcParticleFile */
  string getCgcParticleFile() const { return cgcParticleFile; }
  
  /** @brief Interface for config::P0 */
  double getPtCutoff() const { return P0; }
  
  /** @brief Interfacte to weather to use rings instead of cells for Thermodynamics **/
  bool useRings() const { return useRingsForAVG;  }
  
  /** ------------------------------- */

  /** -------- output options ------- */   
  /** @brief Interface for config::outputSwitch_offlineReconstructionData */
  bool doOutput_offlineReconstructionData() const { return outputSwitch_offlineReconstructionData; }
  
  /** @brief Interface for config::pathdirOfflineData */
  string getPathdirOfflineData() const { return pathdirOfflineData; }

  /** @brief Interface for config::outputSwitch_progressLog */
  bool doOutput_progressLog() const { return outputSwitch_progressLog; }  
  
  /** @brief Interface for config::outputScheme */
  OUTPUT_SCHEME getOutputScheme() const { return outputScheme; }

  /** @brief Interface for config::outputSwitch_detailedParticleOutput */
  bool doOutput_detailedParticleOutput() const { return outputSwitch_detailedParticleOutput; }
  
  /** @brief Interface for config::outputSwitch_detailedParticleOutputAllSteps */
  bool doOutput_detailedParticleOutputAllSteps() const { return outputSwitch_detailedParticleOutputAllSteps; }

  /** @brief Interface for config::outputSwitch_detailedParticleOutputAllSteps */
  bool doOutput_detailedParticlePosMomOutput() const { return outputSwitch_detailedParticlePosMomOutput; }

  /** @brief Interface for config::outputSwitch_particleSpectra */
  bool doOutput_particleSpectraOutput() const { return outputSwitch_particleSpectra; }

  /** @brief Interface for config::outputSwitch_particleSpectraOutputAllSteps */
  bool doOutput_particleSpectraOutputAllSteps() const { return outputSwitch_particleSpectraAllSteps; }
  /** ------------------------------- */
  
  /** -------- heavy quark options ------- */  
  
  /** @brief Interface for config::initialHeavyQuarkFile */
  string getHeavyQuarkParticleFile() const {return initialHeavyQuarkFile; }
  
  /** @brief Is the filename of a seperate heavy quark distribution file given as input */
  bool isSeperateHeavyQuarkParticleFile() const 
  {
    if( initialHeavyQuarkFile == "-" )
      return false;
    else
      return true;
  }
  
  /** @brief Interface for config::initialHeavyQuarkFileTestparticles */
  double getHeavyQuarkParticleFileTestparticles() const {return initialHeavyQuarkFileTestparticles; }
  /** ------------------------------------ */
  
  /** -------- miscellaneous options ------- */ 
  bool repeatTimesteps() const { return switch_repeatTimesteps; }

  /** @brief interface to geometric collisions */  
  bool geometricCollisions() const {return switch_geometricCollisions;}
  
  /** @brief Interface for config::interpolationBorder */
  double getMFPInterpolationBorder() const {return interpolationBorder;}
  
  /** @brief Interface for config::jetMfpComputationSwitch */
  JET_MFP_COMPUTATION_TYPE getJetMfpComputationType() const {return jetMfpComputationSwitch;}
  
  /**@brief Interface for config::scaleTimesteps */
  double getScaleTimesteps() const {return scaleTimesteps;}
  
  /** @brief Interface for config::ANATIMESTEP_MODE */
  ANATIMESTEP_MODE getAnaTimeStepsMode() const { return anaTimeStepsMode; }    
  /** ------------------------------------ */
  
  /** -------- hydro analysis options ------- */  
  /** @brief interface to the total length of the box for the hydro analysis in x-direction */  
  double getBoxLengthForHydroX() const {return boxLengthForHydroX;}
 
  /** @brief interface to the total length of the box for the hydro analysis in y-direction */  
  double getBoxLengthForHydroY() const {return boxLengthForHydroY;}
  
  /** @brief interface to the total length of the box for the hydro analysis in z-direction */  
  double getBoxLengthForHydroZ() const {return boxLengthForHydroZ;}
  
  /** @brief interface to the number of column for midrapidity for the hydro analysis in X-direction */    
  double getNumberColumnMidRapForHydroX() const {return numberColumnMidRapForHydroX;}
  
  /** @brief interface to the number of column for midrapidity for the hydro analysis in Y-direction */    
  double getNumberColumnMidRapForHydroY() const {return numberColumnMidRapForHydroY;}  
  
  /** @brief interface to the number of column for the arrow for midrapidity for the hydro analysis in X-direction */    
  double getNumberColumnMidRapArrowForHydroX() const {return numberColumnMidRapArrowForHydroX;}
  
  /** @brief interface to the number of columns for the arrow for midrapidity for the hydro analysis in Y-direction */    
  double getNumberColumnMidRapArrowForHydroY() const {return numberColumnMidRapArrowForHydroY;}    
  
  /** @brief interface to the length of the space rapidity for midrapidity calculations */   
  double getSpaceRapZRange() const {return spaceRapZRange;}

 protected:
   /** ----- auxiliary routines ----- */
   /** @brief Sort the options into groups */   
   void groupProgramOptions();
   
   /** @brief Actually define all the options */
   void initializeProgramOptions();
   
   /** @brief Read command line arguments and settings provided via a configuration file and via the commmand line */
   void readProgramOptions(const int argc, char* argv[]);
   
   /** @brief Do some processing on program options that have previously been read via configPrototype::readProgramOptions */
   void processProgramOptions();
   
   /** @brief Print a complete configuration file using all current parameter values */
   void printUsedConfigurationParameters();   
   
   /** @brief Do some checks on user-provided options and parameters */
   void checkOptionsForSanity();
   
   /** @brief Some processing of heavy quark options */
   void processHeavyQuarkOptions();
   /** ------------------------------ */
   
   /** ------ boost::program_options objects ------- */ 
   // base class provides:
   //    po::options_description command_line_options;
   //    po::options_description config_file_options;
   //    po::options_description visible_options;
   //    po::positional_options_description pos_options;
   //    po::variables_map vm;
   //    po::options_description usage_information;
   //    po::options_description hidden_options;
   // 
   //    po::options_description general_options;
   //    po::options_description simulation_parameters;
   //    po::options_description output_options;
   //    po::options_description misc_options;
   //    po::options_description heavy_quark_options;
   
   po::options_description collision_parameters;
   po::options_description initial_state_options;   
   /** ------ boost::program_options objects ------- */ 
   
   
  /** ------ general options ------- */  
  // base class provides:
  // string jobname;
  // long seed;
  /** ------------------------------ */
  
  /** ---- collision parameters ---- */  
  /** @brief Mass number of nucleus A  */
  double A;  
  
  /** @brief Atomic number, i.e. number of protons, of nucleus A */
  double Aatomic;
  
  /** @brief Mass number of nucleus B  */
  double B;
  
  /** @brief Atomic number, i.e. number of protons, of nucleus B */
  double Batomic;
  
  /** @brief Center of mass energy per NN pair in GeV */
  double sqrtS;   
  
  /** @brief Impact parameter in fm */
  double impact;  
  /** ------------------------------ */
  
  /** ---- simulation parameters ---- */ 
  // base class provides:
  // double runtime  
  // int testparticles
  // int N_light_flavors_input
  // int N_heavy_flavors_input
  
  /** @brief Energy density for freeze out (in GeV/fm^3)
   * Particles in regions with energy densities below this threshold will be freely streaming
   */
  double freezeOutEnergyDensity;

  /** @brief Switch to set the cell size in transverse direction*/ 
  bool isSetTransCellSize;
  
  /** @brief Set the cell size in transverse direction. If no input, then default value is taken. */
  double transCellSize;     
  
  /** @brief Set the full size in transverse direction */
  double transLength;   
  
  /** @brief Use only some particles from particle input list. */
  int  requiredNinCell; 
  
  /** @brief Boost invariant initial sampling */
  double  minimumEtaInitial;
  double  maximumEtaInitial;
  int  nEtaBinsInitial;
  bool randomizeAzimuth;
  bool homogenSpace;
  double tau0;
  /** ------------------------------- */
  /** @brief For Attractor-Tests */
  int InitialAttractorCondition;
  
    
  /** @brief if 1, then use rings instead of cells **/
  bool useRingsForAVG;
  
  /** @brief Parameters for Correlation output **/
  double pTRefMax;
  double pTRefMin;
  double pTHRefMax;
  double pTHRefMin;
  
  /** ---- initial state options ---- */ 
  /** @brief Which type of initial state to use */
  INITIAL_STATE_TYPE initialStateType;
  
  /** @brief Which source to use for the parton distribution functions
   * 0 = built-in GRV parametrization
   * 1 = PDF parametrizations provided by the LHAPDF library
   */
  PDF_SOURCE_TYPE PDFsource;
  
  /** @brief Name of the LHAPDF data set that should be used */
  string LHAPDFdatasetName;
  
  /** @brief Which member of the given LHAPDF data set should be used */
  unsigned short int LHAPDFmember;
  
  /** @brief Whether to use the grid version of the given LHAPDF set */
  bool LHAPDFuseGrid;
  
  /** @brief Whether to use nPDFs (only available in combination with LHAPDF and minijets) */
  bool nuclearPDFs;
  
  /** @brief Name of the nPDF set that is to be used */
  string nuclearPDFdatasetName;
                     
  /** @brief Relative or full path (including filename) of PYTHIA output file with particle data
   * Declare as "-" if particle momenta should be sampled via glauber method
   */
  string pythiaParticleFile;
  
  /** @brief Relative or full path (including filename) of color glass condensate output file with particle data
   * Declare as "-" if particle momenta should be sampled via glauber method
   */
  string cgcParticleFile; 
  
  /** @brief Lower PT-cutoff [GeV] used for minijet initial conditions */
  double P0;  
  
  
  /** -------- output options ------- */ 
  // provided by base class:
  // string standardOutputDirectoryName
  
  /** @brief Specify whether data needed for later offline reconstruction should be written */
  bool outputSwitch_offlineReconstructionData;

  /** @brief Specify whether progress output should be written to a file */
  bool outputSwitch_progressLog;
  
  /** @brief Output schemes to decide which kind of output is printed */
  OUTPUT_SCHEME outputScheme;
  
  /** @brief Specify whether detailed particle output should be written (initial and final) */
  bool outputSwitch_detailedParticleOutput;
  
  /** @brief Specify whether detailed particle output should be written (all steps) */
  bool outputSwitch_detailedParticleOutputAllSteps;

  /** @brief Specify whether detailed active particle position and momenta should be written (all steps) */
  bool outputSwitch_detailedParticlePosMomOutput;

  /** @brief Specify whether particle spectra should be collected and written (initial and final) */
  bool outputSwitch_particleSpectra;
  
  /** @brief Specify whether particle spectra should be collected and written (all steps) */
  bool outputSwitch_particleSpectraAllSteps;
  
  /** @brief Path to which the output needed for offline analysis is written */
  string pathdirOfflineData;
  /** ------------------------------- */

  /** -------- heavy quark options ------- */ 
  // provided by base class:
  // double KggQQbar;
  // double KgQgQ;
  // double kappa_gQgQ;
  // bool couplingRunning;
  // bool isotropicCrossSecGQ;
  // bool constantCrossSecGQ;
  // double constantCrossSecValueGQ;
  // double Mcharm_input;
  // double Mbottom_input;
  
  /** @brief Filename of file which contains initial heavy quark momentum distribution */
  string initialHeavyQuarkFile;
  
  /** @brief Number of testparticles in initial heavy quark momentum distribution file */
  double initialHeavyQuarkFileTestparticles; 
  /** ------------------------------------ */

  /** -------- hydro analysis options ------- */
  /** @brief total length of box for hydro analysis in x-direction */
  double boxLengthForHydroX;
  
  /** @brief total length of box for hydro analysis in y-direction */ 
  double boxLengthForHydroY;
  
  /** @brief total length of box for hydro analysis in z-direction */
  double boxLengthForHydroZ;  
  
  /** @brief number of the columns in X for Mid Rap hydro analysis */  
  int numberColumnMidRapForHydroX;
  
  /** @brief number of the columns in Y for Mid Rap hydro analysis */  
  int numberColumnMidRapForHydroY; 
  
  /** @brief number of the columns in X for Mid Rap hydro analysis */  
  int numberColumnMidRapArrowForHydroX;
  
  /** @brief number of the columns in Y for Mid Rap hydro analysis */  
  int numberColumnMidRapArrowForHydroY;   
  
  /** @brief total length of box for hydro analysis in z-direction */
  double spaceRapZRange;    
  /** ------------------------------------ */  
  
  /** -------- miscellaneous options ------- */ 
  // provided by base class:
  // bool localCluster
  
  /** @brief X where interpolation of MFP is done for E > X*T */
  double interpolationBorder; 
   
  /** @brief How to compute the mean free path high energy particles?
   *
   * 0 = computeMfpDefault = default, i.e. no special treatment
   * 1 = computeMfpIteration = iterative computation
   * 2 = computeMfpInterpolation = use tabulated mfp data and interpolation functions
   */
  JET_MFP_COMPUTATION_TYPE jetMfpComputationSwitch;
  
  /** @brief Whether timesteps are repeated in cases where the probability has been > 1 */
  bool switch_repeatTimesteps;
  
  /** @brief Whether geometric collisions in the edge cells are performed */  
  bool switch_geometricCollisions;
  
  /**  @brief Factor with which timesteps as computed when building the eta-cells are scaled */
  double scaleTimesteps;
  
  /** @brief Which mode for the analyse time steps */
  ANATIMESTEP_MODE anaTimeStepsMode;   
  /** ------------------------------------ */
};

#endif
