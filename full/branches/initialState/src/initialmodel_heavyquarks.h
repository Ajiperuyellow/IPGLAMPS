//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/initialmodel_heavyquarks.h $
//$LastChangedDate: 2013-10-19 00:21:06 +0200 (土, 19 10月 2013) $
//$LastChangedRevision: 1500 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifndef HEAVY_QUARK_INITIAL_DISTRIBUTION_H
#define HEAVY_QUARK_INITIAL_DISTRIBUTION_H

#include <vector>
#include <string>

#include "initialmodel.h"
#include "configuration.h"
#include "particle.h"

enum SEPERATE_HEAVY_QUARK_INITIAL_MODEL { none, pythia, mcAtNLO };


class initialModel_heavyQuarks : public initialModelWS
{
  public:
    initialModel_heavyQuarks ( const config& _config, SEPERATE_HEAVY_QUARK_INITIAL_MODEL _seperateHeavyQuarkInitialModel = none );
    ~initialModel_heavyQuarks() {};
 
    void populateParticleVector( std::vector<Particle>& _particles );
    
  protected:
    /**
     * @brief Set the positions of the particles
     */
    void samplePositions( std::vector<Particle>& _particles );
    
    /**
     * @brief Set the momenta of the particles
     */
    void sampleMomenta( std::vector<Particle>& _particles );
    
  private:
    SEPERATE_HEAVY_QUARK_INITIAL_MODEL seperateHeavyQuarkInitialModel;
    std::string filename_heavyQuarkParticleFile;
    
    double numberOfTestparticles_heavyQuarkParticleFile;
    int numberOfTestparticles_BAMPS;
    
    int numberOfParticlesToGenerate;
};

#endif 
