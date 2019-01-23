//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/initialmodel_cgc.h $
//$LastChangedDate: 2015-06-09 12:29:40 +0200 (火, 09  6月 2015) $
//$LastChangedRevision: 2168 $
//$LastChangedBy: senzel $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifndef INITIALMODEL_CGC_H
#define INITIALMODEL_CGC_H

#include <string>

#include "initialmodel.h"
#include "configuration.h"
#include "particle.h"


class initialModel_CGC : public initialModelWS
{
  public:
    initialModel_CGC( const config& _config);
    ~initialModel_CGC() {};
    
    void populateParticleVector( std::vector<Particle>& _particles );
  
        
  private:
    std::string filename_cgcParticleFile;
    int nParticlesToGenerate;
    int nTestparticles;
};


#endif // INITIALMODEL_CGC_H
