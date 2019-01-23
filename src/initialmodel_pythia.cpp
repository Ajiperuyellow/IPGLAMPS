//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/initialmodel_pythia.cpp $
//$LastChangedDate: 2014-02-12 19:12:33 +0100 (水, 12  2月 2014) $
//$LastChangedRevision: 1619 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <list>
#include <algorithm> // for std::count
#include <sstream>


#include "configuration.h"
#include "initialmodel_pythia.h"
#include "random.h"
#include "particle.h"
#include "FPT_compare.h"

using std::cout;
using std::endl;
using namespace ns_casc;


initialModel_Pythia::initialModel_Pythia( const config& _config ) :
  initialModelWS(_config),
  numberOfTestparticles( _config.getTestparticles() ),
  numberOfParticlesToGenerate( 0 )
{
  impactParameter = _config.getImpactParameter();
  sqrtS_perNN = _config.getSqrtS();
  filename_pythiaParticleFile = _config.getPythiaParticleFile();
  if (!WoodSaxonParameter.Calculate( A, impactParameter, sqrtS_perNN))
  {
    std::string errMsg = "Impact parameter b too large. b > 2 R_A0";
    throw ePythia_error( errMsg );
  }

  cout << "======= Generating data sets for sampling of initial state =======" << endl;
  generateTimeDistributionWS(Tab);
  cout << "++++  Tab = " << Tab << "1/mb" << endl;
  cout << "==================================================================" << endl;
}



void initialModel_Pythia::populateParticleVector( std::vector< Particle >& _particles )
{
  /**
  * Reserve memory for the Particle vector. Due to possible particle creation this size should be chosen somewhat larger
  * than the initial number of particles. Doing this, push_back operations to add new particles won't lead to internal
  * memory re-allocation in the vector, which could possibly be time consuming.
  * However push_back operations are always possible, even when the reserved size is used up (as long as there is physical
  * memory available of course). And modern day compilers probably optimize to the extent where it doesn't really matter -
  * but still, it's good practice.
  */
  _particles.reserve( static_cast<int>( numberOfParticlesToGenerate * 1.2 ) );
  /**
  * Now the particle vector is re-sized to hold the initial number of particles (this is NOT done when reserving memory!).
  * The particles are initialised with the standard constructor, actual attributes such as momentum etc. MUST be set later!
  */
  _particles.resize( numberOfParticlesToGenerate );

  sampleMomenta( _particles );
  samplePositions( _particles );
}




void initialModel_Pythia::samplePositions( std::vector< Particle >& _particles )
{
  cout << "Start sampling of particle positions for particles from PYTHIA." << endl;
  // particles are already read from file in momentum()
  int nmb_of_events = _particles.back().N_EVENT_pp;
  int event_tmp = 0;
  VectorTXYZ Pos_tmp;

  bool *soft_event = 0;
  soft_event = new bool[nmb_of_events + 1]; // array to store which events are soft and which not (true if event is soft)
  for ( int i = 1; i <= nmb_of_events; i++ )
  {
    soft_event[i] = false;
  }
  int *partcl_soft_event = 0;
  partcl_soft_event = new int[nmb_of_events + 1]; // array to store one particle number of each event to obtain the position of the event afterwards for soft particles

  // sample positions for binary collisions/events where hard partons are produced. Sample also if there are also soft particles at this event (reason: scaling behavior of soft and hard partons differ)
  for ( unsigned int j = 0; j < _particles.size(); j++ )
  {
    if ( _particles[j].HARD )
    {
      if ( _particles[j].N_EVENT_pp != event_tmp ) // first hard particle of an event
      {
        sample_TXYZ_one_partcl( _particles[j], soft_event[_particles[j].N_EVENT_pp] );
        Pos_tmp = _particles[j].Pos;
        partcl_soft_event[_particles[j].N_EVENT_pp] = j;
        event_tmp = _particles[j].N_EVENT_pp;
      }
      else // all other hard particles of same event get same positions
      {
        _particles[j].Pos = Pos_tmp;
      }
    }
  }

  // start insert for non vanishing b=Bimp
  if ( impactParameter != 0.0 ) // soft particles are deleted for b!=0, because treatment of these is written for b=0 and would cause an error
  {
    // delete soft particles
    for(unsigned int j = 0; j < _particles.size(); j++ )
    {
      if( !_particles[j].HARD )
      {
        // delete last particle if soft
        while( ( !_particles.back().HARD ) && ( j != _particles.size() - 1  ) ) // if particle j is the last particle in the particle list it is deleted here and the then last in the list below as well, which would be wrong.
        {
          _particles.pop_back();
        }
        _particles[j] = _particles.back();
        _particles.pop_back();
      }
    }
    cout << "Soft particles have been deleted since b!=0 and this would cause an error!" << endl;
  }
  // end of insert
  else
  {
    int sum_soft_particle = 0;
    int event = 1;
    // get positions for soft partons from events which have besides hard also soft partons
    int soft_break = 0;
    for ( unsigned int j = 0; j < _particles.size(); j++ )
    {
      if ( !_particles[j].HARD )
      {
        sum_soft_particle++;
      }
    }
    cout <<  "Soft particles:" << sum_soft_particle << endl;

    int sum_soft_events = 0;
    for ( int i = 1; i <= nmb_of_events; i++ )
    {
      if ( soft_event[i] )
      {
        sum_soft_events++;
      }
    }
    if ( sum_soft_events == 0 && sum_soft_particle != 0 )
    {
      string errMsg = "ERROR: No soft events! This will result in an infinite loop!";
      throw ePythia_error( errMsg );
    }
    cout << "Soft events: " << sum_soft_events << endl;

    for ( unsigned int s = 0; s < _particles.size(); s++ )
    {
      if ( !_particles[s].HARD )
      {
        while ( !soft_event[event] )
        {
          if ( event < nmb_of_events )
          {
            event++;
          }
          else
          {
            soft_break = s;   // there are more soft particles than events for it
            event = 1;
          }
        }

        // event is soft
        _particles[s].Pos = _particles[partcl_soft_event[event]].Pos;
        if ( event < nmb_of_events )
        {
          event++;
        }
        else
        {
          soft_break = s;   // there are more soft particles than events for it
          event = 1;
        }
      }
    }
  }
  cout << "Finished sampling of particle positions for particles from PYTHIA." << endl;
}


// sample position for only one particle with id=number
void initialModel_Pythia::sample_TXYZ_one_partcl( Particle& _particle, bool& soft )
{
  double densityA_max;
  double L_z;
  double p_soft;

  sample_TXYZ_singleParticle( _particle );

  // sample if there are also soft partons at this position:
  const double sigma = 40.0 * 0.1; // p+p cross section in mb, converted to fm^2  
  densityA_max = densityA( impactParameter / 2.0, WoodSaxonParameter.velocity * _particle.Pos.T() );

  L_z = 2.0 / WoodSaxonParameter.gamma * sqrt( pow( WoodSaxonParameter.RA, 2.0 ) - _particle.Pos.Perp2() ); // only valid for central collision b(=impactParameter)=0
  p_soft = 1.0 / (sigma * densityA_max * L_z);

//   if(p_soft > 1.0)
//   {
//     cout << "error, p_soft>1 in init_pos() (but doen't matter; so it's clearly a soft event), p_soft=" << p_soft << endl;
//   }

//   if(impactParameter != 0.0) // b != 0 would cause an error since consideration above is only for b=0
//     soft = false;
//   else

  soft = (ran2() < p_soft);

}







// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;
