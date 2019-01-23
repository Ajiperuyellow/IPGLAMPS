//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/configurationbase.h $
//$LastChangedDate: 2016-09-21 12:14:54 +0200 (Mi, 21. Sep 2016) $
//$LastChangedRevision: 2433 $
//$LastChangedBy: senzel $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


/** @file
 * @brief This file provides global objects and a configuration interface
 */


#ifndef CONFIGURATIONPROTOTYPE_H
#define CONFIGURATIONPROTOTYPE_H

#include <string>
#include <vector>
#include <stdexcept>
#include <iostream>

// #define BOOST_FILESYSTEM_VERSION 3
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "globalsettings.h"
#include "interpolation23.h"

#define HELP_MESSAGE_REQUESTED 1

/** @brief Enumeration type for different methods of cross sections */
enum CROSS_SECTION_METHOD { csMethod_pQCD, csMethod_constCS, csMethod_constEtaOverS, csMethod_constMFP, csMethod_constMixtureCS, csMethod_constMixtureCS_scaledWithLambdaT2, csMethod_variableCS_scaledWithLambdaT2 };


namespace po = boost::program_options;

/**
 * @brief Basic interface for runtime configuration
 *
 * This class provides a basic interface for settings made by the user at runtime. It is based on boost::program_options 
 * and can handle command line arguments as well as INI-type configuration files.
 * 
 * The class configBase provides the basic functionality as well as the handling of some basic options that are 
 * common to (almost) all variants of BAMPS. It can therefore be used as a standalone interface but is rather meant to
 * be extended in the actual BAMPS code by deriving from it.
 * 
 * In order to provide the means of extending its functionality in derived classes, the handling of program options is split
 * into several separate steps such that each of these routines can in principle be extended or modified in a derived
 * version.
 * The usage is as follows:
 *  1. Call constructor configBase::configBase()
 *   The constructor sets default values to options and calls configBase::initializeProgramOptions() that actually adds
 *   various options to the boost::program_options::options_description objects that are used to categorize the user interface.
 *  2. Use the created object to call the member routine configBase::readAndProcessProgramOptions(const int argc, char* argv[])
 *   This in turn calls several subroutines
 *   1. configBase::groupProgramOptions()  This routine adds the boost::program_options::options_description objects to several 
 *    groups that define which options are to be read from the command line, from an input file, etc.
 *   2. configBase::readProgramOptions(argc, argv)  This routine actually parses the command line and / or configuration file
 *   3. configBase::processProgramOptions()  This is where post-processing of options happens that have been read from user provided parameters
 *   4. configBase::checkOptionsForSanity()  Sanity checks of options (out-of-range checks etc.) could go here. Not used at the moment.
 *   5. configBase::processHeavyQuarkOptions()  Post-processing of options that are special to heavy flavor
 *   6. configBase::printUsedConfigurationParameters()  This writes out the current values of all options into an INI-format file
 *    that could be directly used as an input file to run the program with the very same settings again.
 */
class configBase
{
 public:
  /** @brief Constructor that internally reads the provided input file */
  configBase();
  /** @brief Standard constructor */
  ~configBase() {};
  
  /** @brief processes command line arguments and settings provided via a configuration file */
  virtual void readAndProcessProgramOptions(const int argc, char* argv[]);

  /** ------ general options ------- */  
  /** @brief Interface for configPrototype::jobname */
  std::string getJobName() const { return jobName; }
  
  /** @brief Interface for configPrototype::seed */
  long getSeed() const { return seed; }
  /** ------------------------------- */
  
  /** ---- simulation parameters ---- */ 
  /** @brief Interface for configPrototype::runtime */
  double getRuntime() const { return runtime; }
  
  /** @brief Interface for configPrototype::testparticles */
  double getTestparticles() const { return testparticles; }
  
  /** @brief Interface for config::couplingRunning */
  bool isCouplingRunning() const {return couplingRunning;}

  /** @brief Interface for config::maxRunningCoupling */
  double getMaxRunningCoupling() const {return maxRunningCoupling;}
  
  /** @brief Interface for config::fixedCouplingValue */
  double getfixedCouplingValue() const {return fixedCouplingValue;}
  /** ------------------------------- */

  /** -------- output options ------- */   
  /** @brief Interface for configPrototype::standardOutputDirectoryName */
  std::string getStandardOutputDirectoryName() const { return standardOutputDirectoryName; }
  /** ------------------------------- */

  /** -------- 2->3 parameters ------ */
  /** @brief Interface for configPrototype::I23onlineIntegration */
  bool I23onlineIntegrationIsSet() const {return I23onlineIntegration;}
  
  /** @brief Interface for configPrototype::gluonFormationTimeTyp23 */
  std::string get23GluonFormationTimeTyp() const {return  gluonFormationTimeTyp23; }
  
  /** @brief Interface for configPrototype::matrixElement23 */
  std::string getMatrixElement23() const {return  matrixElement23; }
  
  /** @brief Interface for configPrototype::matrixElement23_22qt */
  bool isMatrixElement23_22qt() const {return  matrixElement23_22qt; }
  
  /** @brief Interface for configPrototype::md2_counter_term_in_I23 */
  bool isMd2CounterTermInI23() const {return  md2_counter_term_in_I23; }
  
  /** @brief Interface for configPrototype::fudge_factor_lpm_23 */
  double get23FudgeFactorLpm() const {return  fudge_factor_lpm_23; }
  
  /** @brief Interface for configPrototype::interpolation23_mode */
  INTERPOLATION_MODE getInterpolation23Mode() const {return  interpolation23_mode; }
  
  /** @brief Interface for configPrototype::kappa23LightPartons */
  double getKappa23LightPartons() const {return  kappa23LightPartons; }
  
  /** @brief Interface for configPrototype::kappa23HeavyQuarks */
  double getKappa23HeavyQuarks() const {return  kappa23HeavyQuarks; }
  
  /** @brief Interface for configPrototype::K23LightPartons */
  double getK23LightPartons() const {return  K23LightPartons; }
  
  /** @brief Interface for configPrototype::K23HeavyQuarks */
  double getK23HeavyQuarks() const {return  K23HeavyQuarks; }
  
  /** ------------------------------- */

  /** -------- miscellaneous options ------- */  
  /** @brief Interface for configPrototype::switch_22 */
  bool doScattering_22() const { return switch_22; }
  /** @brief Interface for configPrototype::switch_23 */
  bool doScattering_23() const { return switch_23; }
  /** @brief Interface for configPrototype::switch_32 */
  bool doScattering_32() const { return switch_32; }
  /** @brief Interface for configPrototype::Kfactor_light */
  double getKfactor_light() const {return Kfactor_light;}
  /** @brief Interface for N_light_flavors_input */
  int getN_light_flavors_input() const {return N_light_flavors_input;}

  /** ------------------------------------ */
  
  /** -------- hydro analysis options ------- */  
  /** @brief interface to the hydro analysis switch */  
  bool getAnalyseForHydro() const {return analyseForHydro;}
  
  /** @brief interface to the number of columns in X for the hydro analysis */
  double getNumberColumnForHydroX() const {return numberColumnForHydroX;}
  
  /** @brief interface to the number of columns in Y for the hydro analysis */
  double getNumberColumnForHydroY() const {return numberColumnForHydroY;}
 
  /** @brief interface to the number of columns in Z for the hydro analysis */
  double getNumberColumnForHydroZ() const {return numberColumnForHydroZ;} 

  /** @brief interface to the number of columns in X for arrow for the hydro analysis */
  double getNumberColumnArrowForHydroX() const {return numberColumnArrowForHydroX;}  

  /** @brief interface to the number of columns in Y for arrow for the hydro analysis */
  double getNumberColumnArrowForHydroY() const {return numberColumnArrowForHydroY;}  
  
  /** @brief interface to the number of columns in Z for arrow for the hydro analysis */
  double getNumberColumnArrowForHydroZ() const {return numberColumnArrowForHydroZ;}    
  
  double getNumberColumnMidRapForHydroX() const {return NumberColumnMidRapForHydroX;}  

  double getNumberColumnMidRapForHydroY() const {return NumberColumnMidRapForHydroY;}  

  double getNumberColumnMidRapForHydroZ() const {return NumberColumnMidRapForHydroZ;}      

  double getNumberColumnMidRapArrowForHydroX() const {return NumberColumnMidRapArrowForHydroX;}  

  double getNumberColumnMidRapArrowForHydroY() const {return NumberColumnMidRapArrowForHydroY;}  

  double getNumberColumnMidRapArrowForHydroZ() const {return NumberColumnMidRapArrowForHydroZ;}       
 
  

  /** ------------------------------------ */  
  
  /** -------- cross section options ------- */ 
  /** @brief Interface for config::CROSS_SECTION_METHOD */
  CROSS_SECTION_METHOD getCrossSectionMethod() const { return crossSectionMethod; }
  
  /** @brief Interface for config::isotropicCrossSection */  
  bool isIsotropicCrossSection() const {return isotropicCrossSection;}

  /** @brief Interface for config::inputCrossSectionValue */  
  double getInputCrossSectionValue() const {return inputCrossSectionValue;}
  
  /** ------------------------------------ */   
  
  /** -------- heavy quark options ------- */ 
  /** @brief Interface for config::KggQQbar */
  double getKggQQb() const {return KggQQbar;}
  
  /** @brief Interface for config::KgQgQ */
  double getKgQgQ() const {return KgQgQ;}
  
  /** @brief Interface for config::kappa_gQgQ */
  double getKappa_gQgQ() const {return kappa_gQgQ;}
  
  /** @brief Interface for config::isotropicCrossSecGQ */
  bool isIsotropicCrossSecGQ() const {return isotropicCrossSecGQ;}
  
  /** @brief Interface for config::constantCrossSecGQ */
  bool isConstantCrossSecGQ() const {return constantCrossSecGQ;}
  
  /** @brief Interface for config::constantCrossSecValueGQ */
  double getConstantCrossSecValueGQ() const {return constantCrossSecValueGQ; }
  /** ------------------------------------ */

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
   
   /** @brief Convert the option descriptions into the INI format of configuration files */
   void printOptionsDescriptionToIniFormat( const boost::program_options::options_description& _desc, boost::filesystem::ofstream& _file );
   
   /** @brief Explicitly convert boost::any data as used in boost::program_options to strings */
   std::string convertBoostAnyToString( const boost::program_options::variable_value& _arg );
  
  /** @brief Do some checks on user-provided options and parameters */
  void checkOptionsForSanity();
  
  /** @brief Create output directory if necessary */
  void checkAndCreateOutputDirectory( boost::filesystem::path& _dir );
  
  /** @brief Some processing of heavy quark options */
  void processHeavyQuarkOptions();
  /** ------------------------------ */
  
  /** ------ boost::program_options objects ------- */ 
  /** @brief Category that groups all options that can be provided via the command line */
  po::options_description command_line_options;
  /** @brief Category that groups all options that can be provided via the configuration file */
  po::options_description config_file_options;
  /** @brief Category that groups all options that are visible in a detailed help message */
  po::options_description visible_options;
  
  /** @brief Positional option used for the configuration file name
   * 
   * Positional options are options that can be given (on the command line) without specifying a parameter name, e.g.
   * the name of the configuration file can be given as ./my_executable inputfile instead of ./my_executable --config-file=inputfile
   */
  po::positional_options_description pos_options;

  /** @brief Variable map that stores all parameter values after command line and configuration file have been parsed */
  po::variables_map vm;

  /** @brief Group options that provide usage information */
  po::options_description usage_information;
  /** @brief Group options that are "hidden" from detailed help messages */
  po::options_description hidden_options;
  
  /** @brief Group general options */
  po::options_description general_options;
  /** @brief Group simulation parameters */
  po::options_description simulation_parameters;
  /** @brief Group output options */
  po::options_description output_options;
  /** @brief Group options for 2->3 processes */
  po::options_description parameters23;
  /** @brief Group miscellaneous options */
  po::options_description misc_options;
  /** @brief Group hydro analysis options */
  po::options_description hydroAnalysis_options;  
  /** @brief Group hydro simulation options */
  po::options_description crossSection_options;    
  /** @brief Group heavy quark options */
  po::options_description heavy_quark_options;
  /** ------ boost::program_options objects ------- */ 
  
  
  /** ------ general options ------- */  
  /** @brief The name of the current job - assigned to output files */
  std::string jobName;
  
  /** @brief Initial seed for the random number generator. seed = 0 forces the random generation of a seed */
  long seed;
  /** ------------------------------ */

  
  /** ---- simulation parameters ---- */ 
  /** @brief Total simulated runtime in fm/c */
  double runtime;  

  /** @brief Number of testparticles per real particle */
  double testparticles;
  
  /** @brief number of active light flavors 
   * ( 0: only gluons, 1: including up, 2: including down, 3: including strange)
   */
  int N_light_flavors_input;
  /** @brief number of active heavy flavors 
   * ( 0: no charm and bottom, 1: only charm, 2: charm and bottom), see particleprototype.h, global parameter 
   */
  int N_heavy_flavors_input;
  
  /** @brief Whether a running coupling is employed for all process */
  bool couplingRunning;
  
  /** @brief value of fixed coupling constant within whole simulation */
  double fixedCouplingValue;
  
  /** @brief maximum value of runnning coupling */
  double maxRunningCoupling;
  
  /** ------------------------------- */

  
  /** -------- output options ------- */ 
  /** @brief Directory to which general output should be written */
  std::string standardOutputDirectoryName;
  
  /** @brief Directory to which general output should be written, which is always local */
  std::string localOutputDirectoryName;  
  /** ------------------------------- */
  
  /** -------- 2->3 parameters ------ */
  /** @brief true if the integration of the total cross section for 2->3 processes is performed online and not read in from table. */
  bool I23onlineIntegration;
  
  /** @brief type of formation time of radiated gluon in 2->3 
   * possibilities are (just give string as input, eg. "bamps_org" ):
   * - bamps_org: as implemented in BAMPS, formation time in frame in which gluon is 
   *   purely transverse (Sigma'') is tau'' = 1/kt, boost lambda to this frame and compare
   * - bamps_org_extended: bamps_org extended to massive sector, formation time in frame 
   *   in which gluon is purely transverse is tau'' = kt / ( kt^2 + x^2 M^2), with x = k+ / P_hq+
   * - bamps_org_extended_xwE: as bamps_org_extended, but with definition x = omega / E_hq
   * - compare_cm_frame: compare formation time in center of mass frame, boost lamda to cm frame, 
   *   here tau' = (2) omega / ( kt^2 + x^2 M^2)
   * - compare_lab_frame: compare formation time in lab frame, boost kt, omega, x to lab frame, 
   *   here tau = (2) omega / ( kt^2 + x^2 M^2)
   * 
   * all these formulas miss a factor 2 in the formation time since it is omitted in the original 
   * BAMPS version, but according to most literature it should be there... It can be added by 
   * setting fudge_factor to 2
   */
  std::string gluonFormationTimeTyp23;
  
  /** @brief matrix element which is used for 2->3 scattering
   * possibilities are (just give string as input, eg. "GBimproved" ):
   * - GBimproved: Gunion-Bertsch improved version which is correct also at forward and mid-rapidity
   * - GBorg: original Gunion-Bertsch version included in BAMPS (is not correct)
   * - KPRapprox: KPR matrix element in the limit k5->0, see notes and http://arxiv.org/abs/1109.5539
   */
  std::string matrixElement23;
  
  /** @brief Whether the gluon propagator in the 2->2 part of the radiative cross section is approximated by 1/qt^4. If false the more exact expression 1/t^2 is employed. Applies only to matrixElement23 = GBimproved. For GBorg always the approximated propagator is used. */
  bool matrixElement23_22qt;
  
  /** @brief Whether a counter term is applied which resembles the old prescription of Xu for masssless particles */
  bool md2_counter_term_in_I23;
  
  /** @brief factor for LPM effect to play around 
    * this factor is used to play around with the cutoff
    * instead of k_t > gamma / lambda a modified cutoff k_t > X gamma / lambda
    * is used, where X is the fudge factor
    */
  double fudge_factor_lpm_23;
  
  /** @brief Stores the mode of interpolation
   * possibilities are (just give number as input )
   * 0 - "log_interpol": (linear) interpolation in log(I23) (original mode),
   * 1 - "lin_interpol": linear interpolation in I23,
   * 2 - "mix_interpol": superposition of both methods
   */
  INTERPOLATION_MODE interpolation23_mode;
  
  /** @brief Kappa for Debye screening for 2->3 process with only light partons (e.g. g + q -> g + q + g) */
  double kappa23LightPartons;
  
  /** @brief Kappa for Debye screening for 2->3 process involving a heavy quark (e.g. g + Q -> g + Q + g) */
  double kappa23HeavyQuarks;
  
  /** @brief K factor to scale the total cross section for 2->3 process with only light partons (e.g. g + q -> g + q + g) */
  double K23LightPartons;
  
  /** @brief K factor to scale the total cross section for 2->3 process involving a heavy quark (e.g. g + Q -> g + Q + g) */
  double K23HeavyQuarks;
  /** ------------------------------- */

  
  /** -------- miscellaneous options ------- */   
  /** @brief Whether 2->2 processes should be used */
  bool switch_22;
  
  /** @brief Whether 2->3 processes should be used */
  bool switch_23;
  
  /** @brief Whether 3->2 processes should be used */
  bool switch_32;

  /** @brief K factor for light parton processes */
  double Kfactor_light;

  /** ------------------------------------ */
  
  /** -------- hydro analysis options ------- */ 
  /** @brief switches hydro analysis ON or OFF*/
  bool analyseForHydro;    

  /** @brief number of columns in X for hydro analysis */
  int numberColumnForHydroX;  
  
  /** @brief number of columns in Y for hydro analysis */
  int numberColumnForHydroY;   
  
  /** @brief number of columns in Z for hydro analysis */
  int numberColumnForHydroZ;     
  
  /** @brief number of columns in X for arrow for hydro analysis */
  int numberColumnArrowForHydroX; 
  
  /** @brief number of columns in Y for arrow for hydro analysis */
  int numberColumnArrowForHydroY; 
  
  /** @brief number of columns in Z for arrow for hydro analysis */
  int numberColumnArrowForHydroZ;
  
  double NumberColumnMidRapForHydroX; 

  double NumberColumnMidRapForHydroY; 

  double NumberColumnMidRapForHydroZ;  

  double NumberColumnMidRapArrowForHydroX; 

  double NumberColumnMidRapArrowForHydroY; 
  
  double NumberColumnMidRapArrowForHydroZ;  
  

  /** ------------------------------------ */  
  
  /** -------- cross section options ------- */ 
  /** @brief Which method of cross sections you want to use */
  CROSS_SECTION_METHOD crossSectionMethod;  
        
  /** @brief whether all cross sections are isotropic */
  bool isotropicCrossSection;
        
  /** @brief an arbitrary value for a kind of cross section. It can be a cross section in [mb], a constant eta over s value or a constant mean free path in [fm].
   * The value it is chosen depends on the method of cross section.
   */
  double inputCrossSectionValue;
  
  /** ------------------------------------ */   
  
  /** -------- heavy quark options ------- */ 
  /** @brief K factor for process g + g -> Q + Qbar */
  double KggQQbar;
  
  /** @brief K factor for process g + Q -> g + Q */
  double KgQgQ;
  
  /** @brief Kappa for Debye screening for process g + Q -> g + Q, usually 0.2 (Peshier,Gossiaux) */
  double kappa_gQgQ;
  
  /** @brief Whether an isotropic momentum sampling is employed for process g + Q -> g + Q */
  bool isotropicCrossSecGQ;
  
  /** @brief Whether a constant cross section is employed for process g + Q -> g + Q */
  bool constantCrossSecGQ;
  
  /** @brief Value of constant cross section for process g + Q -> g + Q */
  double constantCrossSecValueGQ;
  
  /** @brief Mass of charm quarks */
  double Mcharm_input;
  
  /** @brief Mass of bottom quarks */
  double Mbottom_input;
  /** ------------------------------------ */
};


/**
 * This routine checks whether the requested output directory exists and creates it in case it does not. 
 */
inline void configBase::checkAndCreateOutputDirectory(boost::filesystem::path& _dir)
{
  if ( boost::filesystem::exists( _dir ) )
  {
    if ( boost::filesystem::is_directory( _dir ) )
    {
      return;
    }
    else // this handles the rare case where a file (!) with the requested name exists
    {
      boost::filesystem::path renamePath( _dir.string() + ".backup" );
      std::cout << "File with name " << _dir.string() << " blocks the creation of an output folder for offline reconstruction." << std::endl;
      std::cout << "It is renamed to " << renamePath.string() << std::endl;
      boost::filesystem::rename( _dir, renamePath );
      boost::filesystem::create_directory( _dir );       
    }
  }
  else
  {
//     std::cout << "Creating output folder: " << _dir.string() << std::endl;
    boost::filesystem::create_directory( _dir );    
  }
}



/** @brief exception class for handling unexpected critical behaviour within the configuration of the BAMPS run  */
class eConfig_error : public std::runtime_error
{
  public:
    explicit eConfig_error(const std::string& what) : std::runtime_error(what) {};
    
    virtual ~eConfig_error() throw() {};
};


#endif
