//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/initialmodel_simple.h $
//$LastChangedDate: 2016-11-30 22:57:48 +0100 (水, 30 11月 2016) $
//$LastChangedRevision: 2476 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifndef INITIALMODEL_SIMPLE_H
#define INITIALMODEL_SIMPLE_H

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
class initialModel_simple : public initialModelWS
{
  public:
    initialModel_simple( const config& _config );
    ~initialModel_simple() {};
    
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
    
    int numberOfTestparticles;
    int numberOfParticlesToGenerate;
 
    std::string filename_ParticleFile;
};

/** @brief exception class for handling unexpected critical behaviour within generation of mini-jet initial distributions  */
class eSimple_error : public std::runtime_error
{
  public:
    explicit eSimple_error(const std::string& what) : std::runtime_error(what) {};
    
    virtual ~eSimple_error() throw() {};
};


#endif // INITIALMODEL_SIMPLE_H
