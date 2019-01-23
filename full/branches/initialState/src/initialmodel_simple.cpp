//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/initialmodel_simple.cpp $
//$LastChangedDate: 2016-12-02 12:09:23 +0100 (金, 02 12月 2016) $
//$LastChangedRevision: 2477 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <list>
#include <algorithm> // for std::count

#include <sstream>
#include <list>

#include "configuration.h"
#include "initialmodel_simple.h"
#include "random.h"
#include "particle.h"
#include "FPT_compare.h"

using std::cout;
using std::endl;
using namespace ns_casc;


initialModel_simple::initialModel_simple( const config& _config ) :
  initialModelWS(_config),
  numberOfTestparticles( _config.getTestparticles() ),
  numberOfParticlesToGenerate( 0 )
{
  impactParameter = _config.getImpactParameter();
  sqrtS_perNN = _config.getSqrtS();
  filename_ParticleFile = _config.getPythiaParticleFile();
  cout << "simple Initial model: " << filename_ParticleFile << endl;
  
  cout << "WOOD SAXON:" << endl;
  if (!WoodSaxonParameter.Calculate( A, impactParameter, sqrtS_perNN))
  {
    std::string errMsg = "Impact parameter b too large. b > 2 R_A0";
    throw eSimple_error( errMsg );
  }
  
  cout << "======= Generating data sets for sampling of initial state =======" << endl;
  generateTimeDistributionWS(Tab);
  cout << "++++  Tab = " << Tab << "1/mb" << endl;
  cout << "==================================================================" << endl;

  std::ifstream countParticles( filename_ParticleFile.c_str() );
  if ( countParticles.good() )
  {
    numberOfParticlesToGenerate = std::count( std::istreambuf_iterator<char>( countParticles ), std::istreambuf_iterator<char>(), '\n' ); // counts number of lines in  particle data file
  }
  else
  {
    string errMsg = "Error at opening particle data file.";
    throw eSimple_error( errMsg );
  }
  cout << "Number of Particles in Datafile: " << numberOfParticlesToGenerate << endl;
  countParticles.close();
}


void initialModel_simple::sampleMomenta( std::vector< Particle >& _particles )
{

  std::ifstream readParticles( filename_ParticleFile.c_str() );
  cout << "Read particle momentum from PYTHIA data file." << endl;
  for ( int i = 0; i < numberOfParticlesToGenerate; i++ )
  {
    int flavTemp;
    double pX, pY, pZ, E;

    // structure of file
    // number of pythia event  is hard?   flavour  energy  px  py  pz  m
    readParticles >> flavTemp >> E >> pX >> pY >> pZ >> _particles[i].m;
    //cout << flavTemp << endl;
    _particles[i].Mom = VectorEPxPyPz(E, pX, pY, pZ);
    _particles[i].FLAVOR = static_cast<FLAVOR_TYPE>( flavTemp );
    
    if ( flavTemp <= 2 * Particle::max_N_light_flavor )
    {
      _particles[i].m = 0;
    }
    
    // to avoid rounding errors compute energy from momenta and mass E^2=p^2+m^2, for light partons this can be different from Pythia's value if in Pythia this parton had a mass
     _particles[i].Mom.E() = sqrt( _particles[i].Mom.vec2() + pow( _particles[i].m, 2.0 ) );
  }
  
  // charm quarks in Pythia have a mass of 1.5 GeV
  // make charm quarks from pythia lighter, if Mcharm is not 1.5 GeV
  if( !FPT_COMP_E( Particle::Mcharm, 1.5 ) )
  {
    const double M_old = 1.5; // charm mass in PYTHIA
    const double M_new = Particle::Mcharm;
    if(M_new > M_old)
    {
      cout << "problem in rhic::init(), charm mass to high. Leads to nan GeV energy." << endl;
    }
    cout << "Make charm quarks lighter. Old mass=" << M_old << "  new mass=" << M_new << endl;
    for( unsigned int j = 0; j < _particles.size(); j++ )
    {
      if( _particles[j].FLAVOR == charm || _particles[j].FLAVOR ==  anti_charm )
      {
        const double pp_old = sqrt( _particles[j].Mom.vec2() );
        const double pp_new = sqrt ( pow ( pp_old, 2.0 ) + pow ( M_old, 2.0 ) - pow ( M_new, 2.0 ) );

        // scaling in order to conserve the energy when making quarks massles
        _particles[j].Mom *= pp_new / pp_old;
        _particles[j].Mom.E() = sqrt( _particles[j].Mom.vec2() + pow( M_new, 2.0 ) );
        _particles[j].m = M_new;
      }
    }
  }
  
  // make bottom quarks from pythia lighter, if Mbottom is not 4.8 GeV
  if( !FPT_COMP_E( Particle::Mbottom, 4.8 ) )
  {
    const double M_old = 4.8;  // bottom mass from pythia
    const double M_new = Particle::Mbottom; // our bottom mass
    if(M_new > M_old)
    {
      cout << "problem in rhic::init(), bottom mass to high. Leads to nan GeV energy." << endl;
    }
    cout << "Make bottom quarks lighter. Old mass=" << M_old << "  new mass=" << M_new << endl;
    for( unsigned int j = 0; j < _particles.size(); j++ )
    {
      if( _particles[j].FLAVOR == bottom || _particles[j].FLAVOR == anti_bottom )
      {
        const double pp_old = sqrt( _particles[j].Mom.vec2() );
        const double pp_new = sqrt ( pow ( pp_old, 2.0 ) + pow ( M_old, 2.0 ) - pow ( M_new, 2.0 ) );

        // scaling in order to conserve the energy when making quarks massles
        _particles[j].Mom *= pp_new / pp_old;
        _particles[j].Mom.E() = sqrt( _particles[j].Mom.vec2() + pow( M_new, 2.0 ) );
        _particles[j].m = M_new;
      }
    }
  }
}

void initialModel_simple::populateParticleVector( std::vector< Particle >& _particles )
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



void initialModel_simple::samplePositions( std::vector< Particle >& _particles )
{
  cout << "Start sampling of particle positions for particles." << endl;
  // particles are already read from file in momentum()

  VectorTXYZ Pos_tmp;

  // sample positions
  for ( unsigned int j = 0; j < _particles.size(); j++ )
  {
    sample_TXYZ_singleParticle( _particles[j] );
    Pos_tmp = _particles[j].Pos;
    _particles[j].Pos = Pos_tmp;
  }

  cout << "Finished simple sampling of particle positions." << endl;
}


// sample position for only one particle with id=number
void initialModel_simple::sample_TXYZ_one_partcl( Particle& _particle, bool& soft )
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

double initialModel_simple::Radius() 
{
  return WoodSaxonParameter.RA;
}

