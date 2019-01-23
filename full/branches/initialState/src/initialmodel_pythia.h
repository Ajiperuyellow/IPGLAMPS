//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/initialmodel_pythia.h $
//$LastChangedDate: 2013-10-19 00:21:06 +0200 (土, 19 10月 2013) $
//$LastChangedRevision: 1500 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifndef INITIALMODEL_PYTHIA_H
#define INITIALMODEL_PYTHIA_H

#include <stdexcept>
#include <vector>
#include <string>

#include "initialmodel.h"
#include "configuration.h"
#include "particle.h"
#include "vegas.h"


/**
 * @brief Class to provide the Pythia initialization
 */
class initialModel_Pythia : public initialModelWS
{
  public:
    initialModel_Pythia( const config& _config );
    ~initialModel_Pythia() {};
    
    void populateParticleVector( std::vector<Particle>& _particles );
    
        
  protected:
    /**
     * @brief Set the positions of the particles
     */
    void samplePositions( std::vector<Particle>& _particles );

    /**
     * @brief Set the momenta of the particles
     *
     * Pure virtual function. Any derived class must at least specify
     * this routine.
     */
    virtual void sampleMomenta( std::vector<Particle>& _particles ) = 0;
    
    /** 
     * @brief Sampling of time and positions of one parton 
     **/
    void sample_TXYZ_one_partcl( Particle& _particle, bool& soft );
    
    /** 
     * @brief Nuclear overlap function
     **/
    double Tab;
    
    int numberOfTestparticles;
    int numberOfParticlesToGenerate;
 
    std::string filename_pythiaParticleFile;
};



/** @brief exception class for handling unexpected critical behaviour within generation of mini-jet initial distributions  */
class ePythia_error : public std::runtime_error
{
  public:
    explicit ePythia_error(const std::string& what) : std::runtime_error(what) {};
    
    virtual ~ePythia_error() throw() {};
};




#endif // INITIALMODEL_PYTHIA_H
