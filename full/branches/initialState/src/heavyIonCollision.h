//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/heavyIonCollision.h $
//$LastChangedDate: 2019-01-05 18:02:58 +0100 (Sa, 05. Jan 2019) $
//$LastChangedRevision: 2917 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** 
 * @file 
 * @brief Declarations for classes and routines that constitute the main program flow of BAMPS 
 */


#ifndef RHIC_H
#define RHIC_H

#include <algorithm>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <math.h>
#include <stdexcept>
#include <vector>


#include "analysis.h"
#include "cellcontainer.h"
#include "configuration.h"
#include "coordinateBins.h"
#include "interpolation22.h"
#include "interpolation23.h"
#include "mfp_data.h"
#include "offlineoutput.h"
#include "ratesmanager.h"
#include "ringstructure.h"
#include "scattering22.h"
#include "scattering23.h"
#include "scattering32.h"


class config;
class analysis;

using std::fstream;


enum simulationType {initialState,simulationState,finalState};

/** 
 * @brief This class provides the routines that constitute the main program flow of BAMPS  
 *  
 * This class provides an interface for the simulation. Usage:
 *    -# Create an object of type heavyIonCollision
 *    -# Call heavyIonCollision::initialize()
 *    -# Call heavyIonCollision::mainFramework
 */
class heavyIonCollision
{
public:
  /** 
   * @brief Constructor
   */
  heavyIonCollision( config * const _theConfig, mfpForHeavyIonCollision* const _theMFP );

  /** 
   * @brief Destructor
   */
  ~heavyIonCollision();

  /** 
   * @brief Initialize particle momenta, positions, etc. 
   */
  void initialize( analysis& aa );

  /** 
   * @brief The main framework of the simulation 
   */
  void mainFramework( analysis& aa );
  /** 
   * @brief Fix the initial screening mass
   */
  void FixInitialScreeningMass(double & md2q,double & md2g);

private:
  /** 
   * @brief A pointer to the global config object that encapsulates
   * user defined settings etc. 
   */ 
  config* const theConfig;

  /** 
   * @brief A pointer to a \p mfpForHeavyIonCollision object that
   * provides interpolated mean free path data for jet particles
   *
   * This is only needed when \sa config::jetMfpComputationSwitch is set
   * to \sa JET_MFP_COMPUTATION_TYPE::computeMfpInterpolation. 
   */
  mfpForHeavyIonCollision* const theMFP;

  /** @brief  interpolation23 object that provides access to tabulated values for the cross section of all 2->3 processes with running coupling */
  interpolation23 theI23_massless;
  interpolation23 theI23_charm_m1;
  interpolation23 theI23_charm_m2;
  interpolation23 theI23_bottom_m1;
  interpolation23 theI23_bottom_m2;
  
  /** @brief  interpolation22 object that provides access to tabulated values for the cross section of all 2->2 processes with running coupling */
  interpolation22 theI22;

  //---- parameters that are copied from the config object for convenience ----
  /** @brief Runtime to simulate [in fm/c] */
  double stop;
  /** @brief Mass number of nucleus A */
  double A;
  /** @brief Atomic number, i.e. number of protons, of nucleus A */
  double Aatomic;
  /** @brief Mass number of nucleus A */
  double B;
  /** @brief Atomic number, i.e. number of protons, of nucleus B */
  double Batomic;
  /** @brief Center of momentum energy per NN pair [GeV] */
  double sqrtS;
  /** @brief Lower pT cutoff [GeV] for mini-jet initial conditions, see config::initialStateType */
  double P0;
  /** @brief Impact parameter [in fm] */
  double Bimp;
  /** @brief Number of test particles per real particle */
  int testpartcl;
  //--------------------------------------------------------------------------

  //--------------------------------------------------------------------------
  /** 
   * @brief Some typical radius
   *
   * was formerly given by WoodSaxonParameter.RA
   */
  double typicalRadius;
  //--------------------------------------------------------------------------


  //--------------------------------------------------------------------------  
  /**
  * @brief Time shifting of particles after initialisation 
  */
  void timeShiftingOfParticles(const double eta_max);

  /**
  * @brief Formation Time for particles 
  */
  void formationTimeForParticles();  
  
  /**
  * @brief Delete Particles from list after initialisation
  */
  void deleteParticlesFromList();

  /**
  * @brief Lists the number of particles of every flavor initialised
  */
  void listParticleNumbersForAllFlavors(const simulationType st, const double time);
  //--------------------------------------------------------------------------
  
  /** 
   * @brief Accumulated number of sampling errors in
   * scattering23::getMomenta23
   *
   * Can give hints on the validity of the envelope function. Should
   * be "small".
   */
  int nGet23Errors;

  /** 
   * @brief Accumulated number of sampling errors in
   * scattering32::getMomenta32
   *
   * Only needed for rejection sampling, 0 for Metropolis sampling per
   * construction. For rejection sampling: can give hints on the
   * validity of the envelope function. Should be "small". 
   */
  int nGet32Errors;
  
  /** 
   * @brief Time shift needed to start simulation at t = 0 fm/c
   *
   * In heavyIonCollision::initialize() particles are first
   * initialized with t = 0 fm/c being the moment when both nuclei are
   * at z = 0, thus particles are also created at t < 0 fm/c (when the
   * nuclei start to overlap). This is later on rectified by virtue of
   * timeshift 
   */
  double timeshift;

  /** 
   * @brief Fine grained Delta-eta used for construction of actual eta
   * bins in coordinateEtaBins::populateEtaBins 
   */  
  double deltaEta_fine;

  /** 
   * @brief The rings that are used for averaging several variables
   * (Debye masses, rates, etc.) 
   */ 
  ringStructure rings;
  
  /** 
   * @brief The sizes and numbers of the cells in transvers direction
   * (Debye masses, rates, etc.) 
   */   
  void buildTransverseCellConfiguration(double &_dx, double &_dy, int &IX, int &IY, double &_translen); 
  
  /** 
  * @brief Build the Eta-Bins
  */  
  void buildLongitudinalCellConfiguration(int NinCell);
  
  /** 
  * @brief Organise the cells and normalize thier rates.
  */ 
  void normalizeCells();
  
   /** 
  * @brief Set some properties of particle vector to zero.
  */  
  void initializeAndResetParticles(analysis& aa);
  
   /** 
  * @brief Check for freezing out criterium and set particles free or not.
  */  
  void findFreeParticles(bool & free, int cellindex);
  void findFreeParticlesCellwise(bool & free, int cellindex);
  
  void handleAllParticlesAreFree(int &Nfree2,int cellindex);
  void handleParticleNumbersRatesDebyeMasses(int & cellindex, int & NinActCell,int & nGluons,bool & free, std::vector< int >&   allParticlesList, std::vector< int >&   gluonList, int & nAllLightQuarks, int & nAllAntiLightQuarks, int & nCharmQuarks, int & nBottomQuarks, int & nAntiCharmQuarks, int & nAntiBottomQuarks);
   /** 
  * @brief Update the particle rates - previous and current - 
  */    
  void updateParticleRatesRings(int cellindex);
  void updateParticleRatesCellwise(int cellindex);
  
  /** 
  * @brief Handle the cells which become edge cells
  */ 
  void handleCellcutCellsRingwise(int cellindex);
  void handleCellcutCellsCellwise(int cellindex);
  
  /** @brief Interface for the output of data needed for later offline reconstruction */
  offlineOutputInterface offlineInterface;
  
  /** @brief Stores the rates of gluons for output needed for offline reconstruction */
  vector< vector<double> > rateForOfflineOutput_gluons;

  /** 
   * @brief Stores the rates of light quarks for output needed for
   * offline reconstruction 
   */
  vector< vector<double> > rateForOfflineOutput_quarks;

  /** 
   * @brief Stores the rates of light anti-quarks for output needed
   * for offline reconstruction 
   */
  vector< vector<double> > rateForOfflineOutput_antiQuarks;
  
  /** 
   * @brief The number of edge cells, possible choices are e.g. 1 for one global edge cell or 8 for an edge cell divided into 8 sectors
   *
   * In order to save computation time, the list of particles in the
   * edge cell (i.e. the list of particles that scatter geometrically)
   * can be split into multiple lists. 
   * The framework is tested for numberOfEdgeCells = 8, i.e. one edge
   * cell for each sector (x > 0, y > 0, z > 0 and x > 0, y > 0, z < 0
   * etc.), it should be rather flexible though. Probably only the
   * indexing in heavyIonCollision::getEdgeCellIndex would need to be
   * adapted.  
   */
  int numberOfEdgeCells; 

  /** 
   * @brief Routine that determines to which cell every particle
   * currently belongs
   */
  void cell_ID( int& NinAct, int& Nfree1, int& NinFormInit, int& NinFormGeom, std::map< int, cellContainer >& _cellsBackup, int& NfrozenOut  );

  /** 
   * @brief 
   */  
  void transversePeriodicHack(vector<Particle>::iterator it);
  void transverseReflectingHack(vector<Particle>::iterator it);
    
  /** 
   * @brief Routine that controls the number of particles in the cells.
   */  
  void controlParticlesInCell();

  /** 
   * @brief Routine for analysis of central cell(s)
   */    
  void analyzeCentralCell( int cellindex ,std::vector< int >& _allParticlesList, int etaSliceIndex, int IX, int IY, const double volume, int nCellsAVG, double time , ringStructure & rings);
  
  /** 
   * @brief Routine that handles the scattering of particles in all
   * cells for a given time step
   */
  void scattering( int& NinActCell, int& Nfree2, bool& again, analysis& aa, std::map< int, cellContainer >& _cellsBackup );

  /** 
   * @brief Routine that initializes form geom.
   */  
  void cleanFormGeom();
  
  /** 
   * @brief Routine that initializes the averages.
   */   
  void populateRingAndCellAverages(double dz, int i);
  /** 
   * @brief Collect particles around a cell in case of less than 30 particles, which would make a shit boost.
   */     
  void collectParticlesInCellWithNeighbors(int cellindex, std::vector< int >& _allParticlesListWNeighbors, int etaSliceIndex, int IX, int IY, int & nCellsAVG);
  
  /** 
   * @brief Handles 3->2 scatterings of particles in a given cell
   */
  void scatt32( cellContainer& _cell, std::vector< int >& allParticlesList, std::vector< int >& gluonList, int& n32, bool& again, analysis& aa );

  /** 
   * @brief Handles 2->2 and 2->3 scatterings of particles in a given
   * cell 
   */
  void scatt2322( cellContainer& _cell, std::vector< int >& allParticlesList, std::vector< int >& gluonList, const double scaleFactor, bool& again, analysis& aa );

  /** 
   * @brief Handles the setting of new momenta etc. for a given 3 -> 2
   * interaction
   */
  void scatt32_utility( scattering32& scatt32_obj, std::list<int>& _cellMembers, std::vector<int>& _allParticlesList, std::vector<int>& _gluonList, const int iscat, const int jscat, const int kscat, int& n32 );
  
  /** 
   * @brief Handles the setting of new momenta etc. for a given 2 -> 3
   * interaction
   */
  int scatt23_utility( scattering23& scatt23_obj, cellContainer& _cell, int iscat, const int jscat );
  
  /** 
   * @brief Handles the setting of new momenta etc. for a given 2 -> 2
   * interaction
   */
  void scatt22_utility( scattering22& scatt22_obj, const int, const int, int& );

  /** 
   * @brief This routine restores the previous state of the cell
   * configurations in case a time step needs to be undone
   */
  void restoreCellBackup( std::map< int, cellContainer >& _cellsBackup );
  
  /** 
   * @brief Iteratively compute the mean free path for a cell
   * containing a jet particle
   */
  double iterateMFP( std::vector< int >& _allParticlesList, std::vector< int >& _gluonList, const int jetID, const double dt, const double dv, const double time );
  
  /** 
   * @brief Cleans up the particle list after possible particle
   * destruction in 3->2 processes
   */
  void removeDeadParticles( analysis& _aa );
  
    /** 
   * @brief Transform some gluons into quarks for tests
   */
  void makeSomeFakeQuarkPairs();
  
  /** 
   * @brief Compute the Delta-t for the first time step 
   */
  double getFirstTimestep( double& _dt, double& _min );
  double getFirstTimestepMinimum( double& _dt, double& _min);
  
  // ----- Auxiliary routines for handling of edge cells -----
  /** 
   * @brief Add given particle to the appropriate edge cell
   */
  void addParticleToEdgeCell( const int _particleID, vector <std::list <int > >& _edgeCell );

  /** 
   * @brief Add given particle to a given edge cell 
   */
  void addParticleToEdgeCell( const int _particleID, const int _edgeCellID, vector <std::list <int > >& _edgeCell );

  /** 
   * @brief Removes given particle from the appropriate edge cell 
   */
  void removeParticleFromEdgeCell( const int _particleID, vector <std::list <int > >& _edgeCell );

  /** 
   * @brief Removes given particle from a given edge cell 
   */
  void removeParticleFromEdgeCell( const int _particleID, const int _edgeCellID, vector <std::list <int > >& _edgeCell );
  
  /** 
   * @brief Get the index of the appropriate edge cell for a given
   * particle ID
   */
  int getEdgeCellIndex( const int _particleID ) const;
  
  /** 
   * @brief Get the index of the appropriate edge cell for a given
   * particle 
   */
  int getEdgeCellIndex( const Particle& _particle ) const;
  
  /** 
   * @brief Move particle from one edge cell to another 
   */
  void swapParticleToNewEdgeCell( const int _particleID, const int _oldEdgeCellID, const int _newEdgeCellID, vector <std::list <int > >& _edgeCell );

  // ---------------------------------------------------------
  
  
  // ---- Routines for hydro ----//
  /** 
   * @brief Prepares the different cross section methods. In case you have pQCD, nothing happens.
   */  
  double getSpecialCS(cellContainer& _cell, const double cellVolume, hydroParticleType & theHydroParticleType, const vector<hydroParticleTypeProperties> & vecType, double & scalingFactor, double & extractedTeff );
 
  /** 
   * @brief Shows the method of the cross section is used.
   */    
  void infoCrossSectionMethod();
 
  /** 
   * @brief Calculates the Tmunu for every particle species within one cell 
   */  
  void getTmunuInCell(cellContainer& _cell, const double cellVolume, hydroParticleType & theHydroParticleType, const vector<hydroParticleTypeProperties> & vecType,
         vector<double>& T00,  vector<double>& T11,  vector<double>& T22,
         vector<double>& T33,  vector<double>& T10,  vector<double>& T20,
         vector<double>& T30,  vector<double>& T21,  vector<double>& T31,
         vector<double>& T32,  vector<double>& N0,  vector<double>& N1,
         vector<double>& N2,  vector<double>& N3,  vector<int>& NN );

    /** 
   * @brief Calculates the Tmunu for every particle species within one cell 
   */  
  void getTmunuAtMidrapidity(double time_, const double cellVolume, hydroParticleType & theHydroParticleType, const vector<hydroParticleTypeProperties> & vecType,
         vector<double>& T00,  vector<double>& T11,  vector<double>& T22,
         vector<double>& T33,  vector<double>& T10,  vector<double>& T20,
         vector<double>& T30,  vector<double>& T21,  vector<double>& T31,
         vector<double>& T32,  vector<double>& N0,  vector<double>& N1,
         vector<double>& N2,  vector<double>& N3,  vector<int>& NN );
  
  /** 
   * @brief Calculates the Tmunu of a single particle within one cell 
   */  
  void getTmunuSingleParticle(int particleIndex, const double cellVolume, hydroParticleType & theHydroParticleType, const vector<hydroParticleTypeProperties> & vecType,
         vector<double>& T00,  vector<double>& T11,  vector<double>& T22,
         vector<double>& T33,  vector<double>& T10,  vector<double>& T20,
         vector<double>& T30,  vector<double>& T21,  vector<double>& T31,
         vector<double>& T32,  vector<double>& N0,  vector<double>& N1,
         vector<double>& N2,  vector<double>& N3,  vector<int>& NN, double eta );
  
  
  /** 
   * @brief Calculates the temperature with different method
   */    
  void getTemperatureAlternative(cellContainer& _cell, const double cellVolume, hydroParticleType & theHydroParticleType, const vector<hydroParticleTypeProperties> & vecType, double & Tsq, double & sumUp, double & EdensGeV4);
  
  /** 
   * @brief Calculates the Tmunu in the central cell
   */  
  void getCentralCellValues(double time);
  
  /** 
   * @brief Calculate hydro observabels in mid-rapidity
   */   
  void getCentralValues(double time);
  
  // ----- Routines for geometric scatterings -----
  /** 
   * @brief Routine that handles the (geometric) scatterings of
   * particles in the edge cell(s) 
   */
  void scatterEdgeParticles( vector< std::list< int > >& _particleList );
  void scatterEdgeParticlesCellwise( vector< std::list< int > >& _particleList );
  
  /** 
   * @brief Perform geometric collisions of particles in a given
   * particle list
   */
  void doGeometricCollisions( std::list<int>& _particleList, const int firstScatteredParticle );

  /** 
   * @brief Computes some values needed to carry out the geometric
   * collisions of two given particles 
   */
  void prepareGeometricCollision( const int iscat, const int jscat, double& M1, double& M2, FLAVOR_TYPE& F1, FLAVOR_TYPE& F2, double& s, double& ct_i, double& ct_j );
  
  /** 
   * @brief Updates the properties of two given particles after a
   * geometric collision 
   */
  void updateAfterGeometricCollision( const int iscat, const int jscat, const VectorEPxPyPz & P1, const VectorEPxPyPz & P2, const double ct_i, const double ct_j );
  
  /** 
   * @brief Computes which will be the first geometric collision in a
   * given list of particles (discards collision time information)
   */
  int updateGeometricCollisionTimesA( std::list<int>& _particleList ) { double dummy = 0; return updateGeometricCollisionTimesA( _particleList, dummy ); }
  
  /** 
   * @brief Computes which will be the first geometric collision in a
   * given list of particles
   */
  int updateGeometricCollisionTimesA( std::list<int>& _particleList, double& _nextCollTime );

  /** 
   * @brief Computes which will be the next geometric collision in a
   * given list of particles (discards collision time information) 
   */
  int updateGeometricCollisionTimesB( std::list<int>& _particleList, int lastUpdatedParticle ) { double dummy = 0; return updateGeometricCollisionTimesB( _particleList, lastUpdatedParticle, dummy ); }
  
  /** 
   * @brief Computes which will be the next geometric collision in a
   * given list of particles
   */
  int updateGeometricCollisionTimesB( std::list< int >& _particleList, int lastUpdatedParticle, double& _nextCollTime );  
  
  /** 
   * @brief Compute the geometric collision time of two given
   * particles 
   */
  double getGeometricCollisionTime( const int, const int );
  // --------------------------------------------
  
  /** 
   * @brief Binomial coefficient (N, k) 
   */
  int binomial( const int N, const int k ) const;
  
  /** 
   * @brief To use the functions of hydroParticleType in the whole heavyIonCollision class
   */  
  hydroParticleType theHydroParticleType;    
};




/** 
 * @brief exception class for handling unexpected critical behaviour
 * within simulations of heavy ion collisions
 */
class eHIC_error : public std::runtime_error
{
  public:
    explicit eHIC_error(const std::string& what) : std::runtime_error(what) {};
    
    virtual ~eHIC_error() throw() {};
};



/** 
 * @brief Utility function that finds an element in a vector and
 * subsequently removes it 
 *
 * Usage: removeElementFromVector<int> ( vector, element)
 *
 * @param[in,out] _vec Vector from which _elementToRemove needs to be erased
 * @param[in] _elementToRemove The element that should be removed
 */
template <class T>
void removeElementFromVector( std::vector<T>& _vec, const T _elementToRemove )
{
  typename std::vector<T>::iterator findIter;
  findIter = std::find( _vec.begin(), _vec.end(), _elementToRemove );
  
  if ( findIter != _vec.end() )
  {
    _vec.erase( findIter );
  }
  else
  {
    std::string errMsg = "Removal of element from vector failed. Unrecoverable error.";
    std::cout << "Removal of element from vector failed." << _elementToRemove << std::endl;
    throw eHIC_error( errMsg );
  }
}


/** 
 * @brief Custom implementation of sign(x)
 * @param[in] _x The value to check
 * @return 0 for _x <= 0, 1 for _x > 0
 */
inline int signFunction( const double _x )
{
  return ( _x > 0 ? 1 : 0 );
}

/**
 * The correct edge cell to which the particle should be added is
 * automatically determined.  
 * @param[in] _particleID ID of the particle to add
 * @param[in] _edgeCell Vector of all possible edge cells
 */
inline void heavyIonCollision::addParticleToEdgeCell( const int _particleID, vector <std::list <int > >& _edgeCell )
{
  int _edgeCellID = ( signFunction( ns_casc::particles[_particleID].Pos.X() ) * 1 ) 
    + ( signFunction( ns_casc::particles[_particleID].Pos.Y() ) * 2 )
    + ( signFunction( ns_casc::particles[_particleID].Pos.Z() ) * 4 );
  _edgeCell[_edgeCellID].push_back( _particleID );
  ns_casc::particles[_particleID].edge = _edgeCellID;
}

/**
 * The ID of the edge cell to which the particle should be added is
 * given via the interface. 
 * @param[in] _particleID ID of the particle to add
 * @param[in] _edgeCellID ID of the edge cell to which the particle
 * should be added (index in _edgeCell vector) 
 * @param[in] _edgeCell Vector of all possible edge cells
 */
inline void heavyIonCollision::addParticleToEdgeCell( const int _particleID, const int _edgeCellID, vector <std::list <int > >& _edgeCell )
{
  _edgeCell[_edgeCellID].push_back( _particleID );
  ns_casc::particles[_particleID].edge = _edgeCellID;
}

/**
 * The correct edge cell from which the particle should be removed is
 * automatically determined. 
 * @param[in] _particleID ID of the particle to remove
 * @param[in] _edgeCell Vector of all possible edge cells
 */
inline void heavyIonCollision::removeParticleFromEdgeCell( const int _particleID, vector <std::list <int > >& _edgeCell )
{
  int _edgeCellID = ns_casc::particles[_particleID].edge;
  _edgeCell[_edgeCellID].remove( _particleID );
}

/**
 * The ID of the edge cell from which the particle should be removed
 * is given via the interface. 
 * @param[in] _particleID ID of the particle to remove
 * @param[in] _edgeCellID ID of the edge cell from which the particle
 * should be removed (index in _edgeCell vector) 
 * @param[in] _edgeCell Vector of all possible edge cells
 */
inline void heavyIonCollision::removeParticleFromEdgeCell( const int _particleID, const int _edgeCellID, vector <std::list <int > >& _edgeCell )
{
  _edgeCell[_edgeCellID].remove( _particleID );
}

/**
 * The ID of the edge cell corresponding to the given particle is
 * computed from its geometrical position 
 * @param[in] _particleID ID of the particle for which the
 * corresponding edge cell ID is to be computed 
 */
inline int heavyIonCollision::getEdgeCellIndex( const int _particleID ) const
{
  return ( ( signFunction( ns_casc::particles[_particleID].Pos.X() ) * 1 ) 
           + ( signFunction( ns_casc::particles[_particleID].Pos.Y() ) * 2 )
           + ( signFunction( ns_casc::particles[_particleID].Pos.Z() ) * 4 ) );
}

/**
 * The ID of the edge cell corresponding to the given particle is
 * computed from its geometrical position 
 * @param[in] _particle Particle for which the corresponding edge cell
 * ID is to be computed 
 */
inline int heavyIonCollision::getEdgeCellIndex( const Particle& _particle ) const
{
  return ( ( signFunction( _particle.Pos.X() ) * 1 ) 
           + ( signFunction( _particle.Pos.Y() ) * 2 ) 
           + ( signFunction( _particle.Pos.Z() ) * 4 ) );
}

/**
 * The IDs of the old and the new edgeCells need to be given explicitly
 * @param[in] _particleID ID of the particle to move from one edge
 * cell to another 
 * @param[in] _oldEdgeCellID ID of the edge cell from which the
 * particle should be removed (index in _edgeCell vector) 
 * @param[in] _newEdgeCellID ID of the edge cell to which the particle
 * should be added (index in _edgeCell vector) 
 * @param[in] _edgeCell Vector of all possible edge cells
 */
inline void heavyIonCollision::swapParticleToNewEdgeCell( const int _particleID, const int _oldEdgeCellID, const int _newEdgeCellID, vector <std::list <int > >& _edgeCell )
{
  _edgeCell[_oldEdgeCellID].remove( _particleID );
  _edgeCell[_newEdgeCellID].push_back( _particleID );
  ns_casc::particles[_particleID].edge = _newEdgeCellID;
}



#endif
// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on; 
