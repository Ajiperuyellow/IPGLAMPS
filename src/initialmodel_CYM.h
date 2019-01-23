//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/initialmodel_CYM.h $
//$LastChangedDate: 2018-06-07 15:37:19 +0200 (Do, 07. Jun 2018) $
//$LastChangedRevision: 2813 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


//--- FILE INCLUDES ATTRACTOR AND CYM INITIAL STATE ---



#ifndef INITIALMODEL_CYM_H
#define INITIALMODEL_CYM_H

#include <stdexcept>
#include <vector>
#include <string>

#include "initialmodel.h"
#include "configuration.h"
#include "particle.h"
#include "vegas.h"
//#include "initialmodel_pythia.h"


/**
 * @brief Class to provide the  offline initialization, which means that files from outside are read in
 */
class initialModel_CYM : public initialModelWS
{
  public:
    initialModel_CYM( const config& _config );
    ~initialModel_CYM() {};
    
  void populateParticleVector( std::vector<Particle>& _particles );
    
  protected:
    /**
     * @brief Set the momenta of the particles
     */
    void sampleMomenta( std::vector<Particle>& _particles );
    /**
     * @brief Set the positions of the particles
     */
    void samplePositions( std::vector<Particle>& _particles );
    /** 
     * @brief Sampling of time and positions of one parton 
     **/
    void sample_TXYZ_one_partcl( Particle& _particle, bool& soft ); 
    /** 
     * @brief Nuclear overlap function
     **/
    double Tab;
    
    /** @brief Typical Radius when not given by WS parameter */
    double Radius();
    
    /** @brief Sample homogeneous in lattice */
    void sampleHomogeneousPositions( double posXmin, double posXmax, double posYmin, double posYmax, double &posX, double &posY );
  
    /** @brief Set the momenta and positions according to the input. */
    int PositionsAndMomenta( std::vector< Particle >& _particles );

    /** @brief Set the momenta and positions according to the input. */
    int PositionsAndMomentaThermal( std::vector< Particle >& _particles );
    
    /** @brief Set the momenta and positions in special way. */
    int PositionsAndMomentaHomogeneousGlasmaFixEDens( std::vector< Particle >& _particles );
    
    /** @brief Parameters of the CYM Model */
    int numberOfTestparticles;
    int numberOfParticlesToGenerate;
    double transCellSize;
    int nEtaBins;
    double minEta;
    double maxEta;
    bool randomiseAzimuth;
    bool homogeniseSpace;
    double InitialTau0;
    
    
 
    std::string filename_ParticleFile;
};

/** @brief exception class for handling unexpected critical behaviour within generation of mini-jet initial distributions  */
class eCYM_error : public std::runtime_error
{
  public:
    explicit eCYM_error(const std::string& what) : std::runtime_error(what) {};
    
    virtual ~eCYM_error() throw() {};
};



/**
 * @brief Class to provide the  offline initialization, which means that files from outside are read in
 */
class initialModel_Attractor : public initialModelWS
{
  public:
    initialModel_Attractor( const config& _config );
    ~initialModel_Attractor() {};
    
  void populateParticleVector( std::vector<Particle>& _particles );
    
  protected:
    
    /** @brief Sample homogeneous in lattice */
    void sampleHomogeneousPositions( double posXmin, double posXmax, double posYmin, double posYmax, double &posX, double &posY );
  
//     /** @brief Set the momenta and positions according to the input. */
//     int PositionsAndMomenta( std::vector< Particle >& _particles );

    /** @brief Set the momenta and positions according to the input. */
    int PositionsAndMomentaThermal( std::vector< Particle >& _particles );
    int PositionsAndMomentaBox( std::vector< Particle >& _particles, double  maxMomentum, double T );
    int PositionsAndMomentaThermalAngularDisturbed( std::vector< Particle >& _particles, double  width );
    
    /** @brief Parameters of the CYM Model */
    int numberOfTestparticles;
    int numberOfParticlesToGenerate;
    double transCellSize;
    int nEtaBins;
    double minEta;
    double maxEta;
    bool randomiseAzimuth;
    bool homogeniseSpace;
    int initialAttractorCondition;
    
 
    std::string filename_ParticleFile;
};

/** @brief exception class for handling unexpected critical behaviour within generation of mini-jet initial distributions  */
class eAttractor_error : public std::runtime_error
{
  public:
    explicit eAttractor_error(const std::string& what) : std::runtime_error(what) {};
    
    virtual ~eAttractor_error() throw() {};
};


#endif // INITIALMODEL_CYM_H
