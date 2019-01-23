//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/initialmodel_pythia_online.h $
//$LastChangedDate: 2013-10-19 00:21:06 +0200 (土, 19 10月 2013) $
//$LastChangedRevision: 1500 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifndef INITIALMODEL_PYTHIA_ONLINE_H
#define INITIALMODEL_PYTHIA_ONLINE_H

#include <stdexcept>
#include <vector>
#include <string>

#include "initialmodel.h"
#include "configuration.h"
#include "particle.h"
#include "vegas.h"
#include "initialmodel_pythia.h"


/**
 * @brief Class to provide the Pythia online initialization, which means that the event generator is run at the start of the BAMPS event
 */
class initialModel_Pythia_online : public initialModel_Pythia
{
  public:
    initialModel_Pythia_online( const config& _config );
    ~initialModel_Pythia_online() {};
    
  protected:
    /**
     * @brief Set the momenta of the particles
     */
    void sampleMomenta( std::vector<Particle>& _particles );
    
  private:
    /**
     * @brief Auxiliary function to convert pythia id for parton flavor to BAMPS flavor scheme
     */
    FLAVOR_TYPE convertPythiaFlavor( const int _pythia_id );
    
    /**
     * @brief Name of PDF set
     */
    string name_pdfset;
    
    /**
     * @brief Whether to use LHAPDF for the parton distribution function
     */
    bool use_lhapdf;
    
    /**
     * @brief Which shadowing parametrization to use, eg. "EPS09", leave empty "" for none
     */
    string shadowing_parametrization_name;
};




#endif // INITIALMODEL_PYTHIA_ONLINE_H
