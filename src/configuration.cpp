//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/configuration.cpp $
//$LastChangedDate: 2018-12-25 21:53:51 +0100 (Di, 25. Dez 2018) $
//$LastChangedRevision: 2914 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <string>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "configBAMPS.h"
#include "configuration.h"
#include "particle.h"
#include "coupling.h"


using std::cout;
using std::endl;
using std::fstream;
using std::ios;
using std::string;
using std::vector;
namespace po = boost::program_options;



/** @brief definition of particles vector, defined extern in configuration.h */
std::vector<Particle> ns_casc::particles;

/** @brief definition of particles vector, defined extern in configuration.h */
std::vector<Particle> ns_casc::hadrons;


/**
 * The only constructor provided for this class. Defaults for all parameters are set in case there is
 * no input file or it cannot be read.
 */
config::config() :
  configBase(),
  //---- program_options groups ----
  collision_parameters("Parameters of the heavy ion collision"),
  initial_state_options("Options and parameters for the initial state used by the BAMPS simulation"),
  // ---- collision parameters ---- 
  A(197),
  Aatomic(79),
  B(197),
  Batomic(79),
  sqrtS(200),
  impact(0),
  // ---- simulation parameters ---- 
  freezeOutEnergyDensity(0.6),
  isSetTransCellSize(false),
  transCellSize(0.5),
  transLength(20),
  minimumEtaInitial(-1.),
  maximumEtaInitial(1.),
  nEtaBinsInitial(1),
  requiredNinCell(20),  //750 at 2x2 cells
  randomizeAzimuth(false),
  homogenSpace(false),
  tau0(0.2005),
  pTRefMax(8.0),
  pTRefMin(0.0),
  pTHRefMax(3.0),
  pTHRefMin(0.4),  
  // ---- initial state options ----
  initialStateType(miniJetsInitialState),
#ifdef LHAPDF_FOUND
  PDFsource( LHAPDF ),
#else
  PDFsource( builtInGRV ),
#endif
  LHAPDFdatasetName("cteq6l"),
  LHAPDFmember(0),
  LHAPDFuseGrid(false),
  nuclearPDFs(false),
  nuclearPDFdatasetName("EPS09"),
  pythiaParticleFile("-"),
  cgcParticleFile("-"),
  P0(1.4),
  InitialAttractorCondition(0),
  // ---- output options ----
  outputSwitch_offlineReconstructionData(false),
  outputSwitch_progressLog( true ),
  outputScheme(no_output),
  outputSwitch_detailedParticleOutput(false),
  outputSwitch_detailedParticleOutputAllSteps(false),
  outputSwitch_detailedParticlePosMomOutput(false),
  outputSwitch_particleSpectra(false),
  outputSwitch_particleSpectraAllSteps(false),
  pathdirOfflineData("offline_output"),
  // ---- heavy quark options ----
  initialHeavyQuarkFile("-"),
  initialHeavyQuarkFileTestparticles(0.0),
  // ---- hydro analysis options ----
  boxLengthForHydroX(1.0),
  boxLengthForHydroY(1.0),
  boxLengthForHydroZ(1.0),
  numberColumnMidRapForHydroX(1),
  numberColumnMidRapForHydroY(1),
  numberColumnMidRapArrowForHydroX(1),
  numberColumnMidRapArrowForHydroY(1),    
  spaceRapZRange(0.2),
  // ---- miscellaneous options ----
  interpolationBorder(50),
  jetMfpComputationSwitch(computeMfpDefault),
  switch_repeatTimesteps(true),
  switch_geometricCollisions(true),    
  scaleTimesteps(0.2),
  useRingsForAVG(0),
  anaTimeStepsMode(generalHICAnaTimeSteps)
{
  // populate the program_options::options_description groups with names and types of the possible parameters
  // this can add to the options_description groups provided by the base class
  initializeProgramOptions();
}


/**
 * @param[in] argc number of command line arguments, passed from the calling process
 * @param[in] argv[] command line arguments, passed from the calling process
 */ 
void config::readProgramOptions ( const int argc, char* argv[] )
{
  // since parsing the command line and the configuration file only relies on the groups "command_line_options",
  // "pos_options" and "config_file_options" to which the derived class adds but which are members of the base class
  // it is sufficient to call the routine of the base class here
  configBase::readProgramOptions( argc, argv );
}


/**
 * This provides the interface for the actual reading and processing of program options. It needs to be explicitly called
 * after an instance of the configPrototype class has been created.
 * 
 * @param[in] argc number of command line arguments, passed from the calling process
 * @param[in] argv[] command line arguments, passed from the calling process
 */
void config::readAndProcessProgramOptions ( const int argc, char* argv[] )
{
  groupProgramOptions();
  readProgramOptions(argc, argv);
  processProgramOptions();
    
  checkOptionsForSanity();
  
  if( ParticlePrototype::N_heavy_flavor > 0 )
  {
    processHeavyQuarkOptions();
  }
  
  printUsedConfigurationParameters();
}



/**
 * This is the place where post-processing of values that have been stored in the boost::program_options::variables_map
 * should go. For example conversion from integer to enum types etc.
 */
void config::processProgramOptions()
{
  configBase::processProgramOptions();
  
  // some special conversions from integer type values to enum values
  if ( vm.count("initial_state.type") )
  {
    if ( vm["initial_state.type"].as<int>() < 8 && vm["initial_state.type"].as<int>() >= 0 )
    {
      initialStateType = static_cast<INITIAL_STATE_TYPE>( vm["initial_state.type"].as<int>() );
    }
    else
    {
      string errMsg = "parameter \"initial_state.type\" out of range";
      throw eConfig_error( errMsg );
    }
  }
  
  if ( vm.count("initial_state.PDFsource") )
  {
    if ( vm["initial_state.PDFsource"].as<unsigned short int>() < 2 && vm["initial_state.PDFsource"].as<unsigned short int>() >= 0 )
    {
      PDFsource = static_cast<PDF_SOURCE_TYPE>( vm["initial_state.PDFsource"].as<unsigned short int>() );
    }
    else
    {
      string errMsg = "parameter \"initial_state.PDFsource\" out of range";
      throw eConfig_error( errMsg );
    }
  }
  
  if ( vm.count("misc.jet_mfp_computation") )
  {
    if ( vm["misc.jet_mfp_computation"].as<int>() < 3 && vm["misc.jet_mfp_computation"].as<int>() >= 0 )
    {
      jetMfpComputationSwitch = static_cast<JET_MFP_COMPUTATION_TYPE>( vm["misc.jet_mfp_computation"].as<int>() );
    }
    else
    {
      string errMsg = "parameter \"misc.jet_mfp_computation\" out of range";
      throw eConfig_error( errMsg );      
    }
  }
  
  if ( vm.count("misc.anaTimeStepsMode") )
  {
    if ( vm["misc.anaTimeStepsMode"].as<int>() < 5 && vm["misc.anaTimeStepsMode"].as<int>() >= 0 )
    {
      anaTimeStepsMode = static_cast<ANATIMESTEP_MODE>( vm["misc.anaTimeStepsMode"].as<int>() );
    }
    else
    {
      string errMsg = "parameter \"misc.anaTimeStepsMode\" out of range";
      throw eConfig_error( errMsg );
    }
  }  
  
  // see if transCellSize set: this is the case if not a default value or no value given (latter happens if no input file is given at all)
  if( !( vm["simulation.transCellSize"].defaulted() || vm["simulation.transCellSize"].empty() ) )
  {
    isSetTransCellSize = true;
  }


  if ( vm.count("output.outputScheme") )
  {
    outputScheme = static_cast<OUTPUT_SCHEME>( vm["output.outputScheme"].as<int>() );
  }
  
  // for heavy quark production automatically change time step mode for analysis
  if( outputScheme == heavy_quark_production )
  {
    anaTimeStepsMode = heavyQuarkAnaTimeSteps;
  }

  if ( outputScheme == correlationsearch )
  {
    anaTimeStepsMode = detailedHICAnaTimeSteps;
  }

}


/**
 * This routine initializes the boost::program_options::options_description objects that give structure to the handling of
 * user provided input.
 * By using a name pattern <group_name>.<option_name> the options can be given in a configuration file in INI format using 
 * group names such as
 * [group_name]
 * option_name1 = value1
 * option_name2 = value2
 */
void config::initializeProgramOptions()
{
  // Group some collision parameters
  collision_parameters.add_options()
  ("collision.A_mass", po::value<double>( &A )->default_value( A ), "mass number of nucleus A" )
  ("collision.A_atomic", po::value<double>( &Aatomic )->default_value( Aatomic ), "atomic number (number of protons) of nucleus A" )
  ("collision.B_mass", po::value<double>( &B )->default_value( B ), "mass number of nucleus B" )
  ("collision.B_atomic", po::value<double>( &Batomic )->default_value( Batomic ), "atomic number (number of protons) of nucleus B" )
  ("collision.impact", po::value<double>( &impact )->default_value( impact ), "impact parameter (fm)")
  ("collision.sqrt_s", po::value<double>( &sqrtS )->default_value( sqrtS ), "center of mass energy per nucleon-nucleon pair (GeV)")
  ;
  
  // Add some simulation parameters
  simulation_parameters.add_options()
  ("simulation.e_freeze", po::value<double>( &freezeOutEnergyDensity )->default_value( freezeOutEnergyDensity ), "freeze out energy density (GeV/fm^3)")
  ("simulation.transCellSize", po::value<double>( &transCellSize )->default_value( transCellSize ), "transverse size of the cells (fm)")  
  ("simulation.transSize", po::value<double>( &transLength )->default_value( transLength ), "full transverse size (fm)")  
  ;
  
  // Group some options related to the initial state
  initial_state_options.add_options()
  ("initial_state.type", po::value<int>()->default_value( static_cast<int>(initialStateType) ), "initial state type (0 = mini-jets, 1 = pythia offline, 2 = cgc, 3 = hydroParametrization, 4 = pythia online, 5 = simple, 6 = CYM, 7 = Attractor)")
  ("initial_state.PDFsource", po::value<unsigned short int>()->default_value( static_cast<unsigned short int>(PDFsource) ), "which source to use for the PDFs ( 0 = built-in GRV, 1 = PDFs from LHAPDF )")
  ("initial_state.LHAPDFset", po::value<string>( &LHAPDFdatasetName )->default_value( LHAPDFdatasetName ), "name of the LHAPDF data set that should be used")
  ("initial_state.LHAPDFmember", po::value<unsigned short int>( &LHAPDFmember )->default_value( LHAPDFmember ), "which member of the LHAPDF set should be used (does not work with Pythia online initial distribution)")
  ("initial_state.LHAPDFgrid", po::value<bool>( &LHAPDFuseGrid )->default_value( LHAPDFuseGrid ), "whether a grid version of the LHAPDF set should be used")
  ("initial_state.nuclearPDF", po::value<bool>( &nuclearPDFs )->default_value( nuclearPDFs ), "whether to use nuclear PDFs (only available together with LHAPDF and mini-jets)")
  ("initial_state.nuclearPDFname", po::value<string>( &nuclearPDFdatasetName )->default_value( nuclearPDFdatasetName ), "name of the nPDF dataset to use (EPS09, EPS09LO, EPS09NLO, EPS08, EKS98)")
  ("initial_state.minijet_P0", po::value<double>( &P0 )->default_value( P0 ), "lower pT cutoff for minijet initial conditions")
  ("initial_state.pythia_file", po::value<string>( &pythiaParticleFile )->default_value( pythiaParticleFile ), "input file providing pythia particle information, needed when initial_state.type = 1")
  ("initial_state.cgc_file", po::value<string>( &cgcParticleFile )->default_value( cgcParticleFile ), "input file providing cgc particle information, needed when initial_state.type = 2")
  ("initial_state.minimumEta", po::value<double>( &minimumEtaInitial )->default_value( minimumEtaInitial ), "minimum eta value for initial sampling.")
  ("initial_state.maximumEta", po::value<double>( &maximumEtaInitial )->default_value( maximumEtaInitial ), "maximum eta value for initial sampling.")
  ("initial_state.nEtaBins", po::value<int>( &nEtaBinsInitial )->default_value( nEtaBinsInitial ), "number of eta bins for initial sampling.")
  ("initial_state.randomizeAzimuth", po::value<bool>( &randomizeAzimuth )->default_value( randomizeAzimuth ), "whether to randomize the momenta initially. Will naturally give a v2 of 0.")  
  ("initial_state.homogenSpace", po::value<bool>( &homogenSpace )->default_value( homogenSpace ), "whether to randomize transverse position. Will naturally give an epsilon2 of 0.")    
  ("initial_state.tau0", po::value<double>( &tau0 )->default_value( tau0 ), "Tau_0 for boost invariant initial sampling. Default is 0.2005 fm.")
  ("initial_state.AttractorInitialState", po::value<int>( &InitialAttractorCondition )->default_value( InitialAttractorCondition ), "If Attractor tests are done, specifiy the Initial state.")

  
  ;
  
  // Add some options related to the program output  
  output_options.add_options()
  ("output.offline", po::value<bool>( &outputSwitch_offlineReconstructionData )->default_value( outputSwitch_offlineReconstructionData ), "write data for later reconstruction via \"offline\" routines")
  ("output.offline_output_dir", po::value<string>( &pathdirOfflineData )->default_value( pathdirOfflineData ), "directory to which the output for \"offline\" data is written (if output.offline true)")
  ("output.progress_log", po::value<bool>( &outputSwitch_progressLog )->default_value( outputSwitch_progressLog ), "write progress information" )
  ("output.outputScheme", po::value<int>()->default_value( static_cast<int>(outputScheme) ), "output scheme id which configures the analysis routines and decides which output is written. The integer for the desired output scheme is given in the OUTPUT_SCHEME enum in configuration.h (for instance 0 = print nothing, 1 = print everything, ...).")
  ("output.particles", po::value<bool>( &outputSwitch_detailedParticleOutput )->default_value( outputSwitch_detailedParticleOutput ), "write detailed particle output (initial and final)")
  ("output.particles_all_steps", po::value<bool>( &outputSwitch_detailedParticleOutputAllSteps )->default_value( outputSwitch_detailedParticleOutputAllSteps ), "write detailed particle output (all steps)")
  ("output.particles_position_momenta", po::value<bool>( &outputSwitch_detailedParticlePosMomOutput )->default_value( outputSwitch_detailedParticlePosMomOutput ), "write detailed active particle position and momenta output (all steps)")
  ("output.spectra", po::value<bool>( &outputSwitch_particleSpectra )->default_value( outputSwitch_particleSpectra ), "collect and write particle spectra (initial and final)")
  ("output.spectra_all_steps", po::value<bool>( &outputSwitch_particleSpectraAllSteps )->default_value( outputSwitch_particleSpectraAllSteps ), "collect and write particle spectra (all steps)")
  ("output.pTRefMax", po::value<double>( &pTRefMax )->default_value( pTRefMax ), "2 particle correlation output reference bin maximum momentum.")
  ("output.pTRefMin", po::value<double>( &pTRefMin )->default_value( pTRefMin ), "2 particle correlation output reference bin minimum momentum.")
  ;
  
  // Add some miscellaneous options
  misc_options.add_options()
  ("misc.repeat_timesteps", po::value<bool>( &switch_repeatTimesteps )->default_value( switch_repeatTimesteps ), "repeat timesteps in cases where the probability has been > 1" ) 
  ("misc.geometricCollisions", po::value<bool>( &switch_geometricCollisions )->default_value( switch_geometricCollisions ), "do geometric collisions in the edge cells" )
  ("misc.interpolation_border", po::value<double>( &interpolationBorder )->default_value( interpolationBorder ), "X where interpolation of MFP is done for E > X*T")
  ("misc.jet_mfp_computation", po::value<int>()->default_value( jetMfpComputationSwitch ), "special treatment for the mean free path of high energy particles")
  ("misc.scale_dt", po::value<double>( &scaleTimesteps )->default_value( scaleTimesteps ), "factor with which timesteps dt are scaled")
  ("misc.useRings", po::value<bool>( &useRingsForAVG )->default_value( useRingsForAVG ), "If true, use rings instead of cells for boosts and m_debye." )  
  ("misc.anaTimeStepsMode", po::value<int>()->default_value( static_cast<int>(anaTimeStepsMode) ), "choose the mode for the timestep analysis ( 0 = general, 1 = for heavy quarks, 2 = for initial hydro phase, 3 = detailed, 4=Correlations )")
  ;

  // Add options for the hydro analysis
  hydroAnalysis_options.add_options()  
  ("hydro.boxLengthForHydroX", po::value<double>( &boxLengthForHydroX )->default_value( boxLengthForHydroX ), "gives the size for the box for the hydro analysis in X direction " ) 
  ("hydro.boxLengthForHydroY", po::value<double>( &boxLengthForHydroY )->default_value( boxLengthForHydroY ), "gives the size for the box for the hydro analysis in Y direction " )
  ("hydro.boxLengthForHydroZ", po::value<double>( &boxLengthForHydroZ )->default_value( boxLengthForHydroZ ), "gives the size for the box for the hydro analysis in Z direction " ) 
  ("hydro.numberColumnMidRapForHydroX", po::value<int>( &numberColumnMidRapForHydroX )->default_value( numberColumnMidRapForHydroX ), "gives the size number of column for the midRapidity hydro analysis in X direction " ) 
  ("hydro.numberColumnMidRapForHydroY", po::value<int>( &numberColumnMidRapForHydroY )->default_value( numberColumnMidRapForHydroY ), "gives the size number of column for the midRapidity hydro analysis in Y direction " ) 
  ("hydro.numberColumnMidRapArrowForHydroX", po::value<int>( &numberColumnMidRapArrowForHydroX )->default_value( numberColumnMidRapArrowForHydroX ), "gives the number of arrow column for the midRapidity hydro analysis in X directions " )
  ("hydro.numberColumnMidRapArrowForHydroY", po::value<int>( &numberColumnMidRapArrowForHydroY )->default_value( numberColumnMidRapArrowForHydroY ), "gives the number of arrow column for the midRapidity hydro analysis in Y directions " )   
  ("hydro.spaceRapZRange", po::value<double>( &spaceRapZRange )->default_value( spaceRapZRange ), "gives the size for the space rapidity range in Z direction for midrapidity analysis " ) 
  ;
  
  // Add heavy quark options
  heavy_quark_options.add_options()
  ("heavy_quark.initialHeavyQuarkFile", po::value<string>( &initialHeavyQuarkFile )->default_value( initialHeavyQuarkFile ), "filename of seperate initial heavy quark momentum distribution")
  ("heavy_quark.initialHeavyQuarkFileTestparticles", po::value<double>( &initialHeavyQuarkFileTestparticles )->default_value( initialHeavyQuarkFileTestparticles ), "number of testparticles in initial heavy quark momentum distribution file")
  ;
}


/**
 * This routine groups the options_description groups into categories that can control which parameters are accessible via
 * the command line or the configuration file, which parameters are visible when requesting detailed help messages etc.
 */
void config::groupProgramOptions()
{
  configBase::groupProgramOptions(); // first add the options already contained in the base class
  
  // Group options that are meant to be provided via the command line
  command_line_options.add(initial_state_options).add(simulation_parameters);  
  
  // Add some groups that are meant to be provided via a configuration file
  config_file_options.add(collision_parameters).add(initial_state_options);
  
  // Add option groups that are to be printed in a detailed help message
  visible_options.add(collision_parameters).add(initial_state_options);;
}



/**
 * This routine is provided in case future implementations should contain sanity checks for the given values of parameters and
 * options.
 */
void config::checkOptionsForSanity()
{
  configBase::checkOptionsForSanity();  // perform sanity check for base class options
  // sanity checks of parameters and options provided by the derived class can go here 
  
  if( outputScheme == heavy_quark_production && N_heavy_flavors_input == 0 )
  {
    string errMsg = "Cannot study heavy quarks if they are not switched on.";
    throw eConfig_error( errMsg );
  }
}



void config::processHeavyQuarkOptions()
{
  configBase::processHeavyQuarkOptions();  
}



/**
 * Print a file in INI format with all current option values. This file could be used a configuration filename
 * in a subsequent run to reproduce the very same option settings.
 * As it is somewhat easier this way, the routine explicitly implements the output of all option groups without calling
 * the corresponding base class routine.
 */
void config::printUsedConfigurationParameters()
{
  string filename = standardOutputDirectoryName + "/" + jobName + "_used_configuration";
  boost::filesystem::path outputFile( filename );
  boost::filesystem::ofstream output( outputFile, ios::trunc );
 
  printOptionsDescriptionToIniFormat( general_options, output );
  printOptionsDescriptionToIniFormat( simulation_parameters, output );
  printOptionsDescriptionToIniFormat( collision_parameters, output );
  printOptionsDescriptionToIniFormat( initial_state_options, output );
  printOptionsDescriptionToIniFormat( output_options, output );
  printOptionsDescriptionToIniFormat( parameters23, output );
  printOptionsDescriptionToIniFormat( misc_options, output );
  printOptionsDescriptionToIniFormat( hydroAnalysis_options, output );
  printOptionsDescriptionToIniFormat( crossSection_options, output );   
  printOptionsDescriptionToIniFormat( heavy_quark_options, output );
}
