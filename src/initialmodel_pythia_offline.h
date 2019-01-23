//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/initialmodel_pythia_offline.h $
//$LastChangedDate: 2013-10-19 00:21:06 +0200 (土, 19 10月 2013) $
//$LastChangedRevision: 1500 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifndef INITIALMODEL_PYTHIA_OFFLINE_H
#define INITIALMODEL_PYTHIA_OFFLINE_H

#include <stdexcept>
#include <vector>
#include <string>

#include "initialmodel.h"
#include "configuration.h"
#include "particle.h"
#include "vegas.h"
#include "initialmodel_pythia.h"


/**
 * @brief Class to provide the Pythia offline initialization, which means that files from Pythia are read in
 */
class initialModel_Pythia_offline : public initialModel_Pythia
{
  public:
    initialModel_Pythia_offline( const config& _config );
    ~initialModel_Pythia_offline() {};
       
  protected:
    /**
     * @brief Set the momenta of the particles
     */
    void sampleMomenta( std::vector<Particle>& _particles );
};




#endif // INITIALMODEL_PYTHIA_OFFLINE_H
