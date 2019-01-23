//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/initialmodel_pythia_offline.cpp $
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
#include <list>

#include "configuration.h"
#include "initialmodel_pythia_offline.h"
#include "random.h"
#include "particle.h"
#include "FPT_compare.h"

using std::cout;
using std::endl;
using namespace ns_casc;


initialModel_Pythia_offline::initialModel_Pythia_offline( const config& _config ) :
  initialModel_Pythia(_config)
{
  cout << "PYTHIA particle data file: " << filename_pythiaParticleFile << endl;

  std::ifstream countPythiaParticles( filename_pythiaParticleFile.c_str() );
  if ( countPythiaParticles.good() )
  {
    numberOfParticlesToGenerate = std::count( std::istreambuf_iterator<char>( countPythiaParticles ), std::istreambuf_iterator<char>(), '\n' ); // counts number of lines in PYTHIA particle data file
  }
  else
  {
    string errMsg = "Error at opening PYTHIA particle data file.";
    throw ePythia_error( errMsg );
  }

  countPythiaParticles.close();
}

void initialModel_Pythia_offline::sampleMomenta( std::vector< Particle >& _particles )
{

  std::ifstream readPythiaParticles( filename_pythiaParticleFile.c_str() );
  cout << "Read particle momentum from PYTHIA data file." << endl;
  for ( int i = 0; i < numberOfParticlesToGenerate; i++ )
  {
    int flavTemp;
    double pX, pY, pZ, E;

    // structure of file
    // number of pythia event  is hard?   flavour  energy  px  py  pz  m
    readPythiaParticles >> _particles[i].N_EVENT_pp >> _particles[i].HARD >> flavTemp >> E >> pX >> pY >> pZ >> _particles[i].m;
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
