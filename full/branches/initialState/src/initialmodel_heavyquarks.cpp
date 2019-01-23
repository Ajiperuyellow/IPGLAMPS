//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/initialmodel_heavyquarks.cpp $
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

#include "configuration.h"
#include "initialmodel_heavyquarks.h"
#include "particle.h"

using std::cout;
using std::endl;
using namespace ns_casc;


initialModel_heavyQuarks::initialModel_heavyQuarks ( const config& _config, SEPERATE_HEAVY_QUARK_INITIAL_MODEL _seperateHeavyQuarkInitialModel ) : 
    initialModelWS( _config ),
    seperateHeavyQuarkInitialModel( _seperateHeavyQuarkInitialModel ), 
    numberOfParticlesToGenerate( 0 )
{
  impactParameter = _config.getImpactParameter();
  sqrtS_perNN = _config.getSqrtS();
  if ( !WoodSaxonParameter.Calculate( A, impactParameter, sqrtS_perNN) )
  {
    std::string errMsg = "Impact parameter b too large. b > 2 R_A0";
    throw eInitialModel_error( errMsg );
  }

  double Tab;
  generateTimeDistributionWS(Tab);  // The time distribution should have been already computed. However, call the routine anyways to be safe (it will handle multiple calls automatically)

  filename_heavyQuarkParticleFile = _config.getHeavyQuarkParticleFile();
  numberOfTestparticles_heavyQuarkParticleFile = _config.getHeavyQuarkParticleFileTestparticles();
  
  std::ifstream countHeavyQuarks( filename_heavyQuarkParticleFile.c_str() );
  if ( countHeavyQuarks.good() )
  {
    cout << "Using heavy quark particle file: " << filename_heavyQuarkParticleFile << endl;
    numberOfParticlesToGenerate = std::count( std::istreambuf_iterator<char>( countHeavyQuarks ), std::istreambuf_iterator<char>(), '\n' ); // counts number of lines in particle data file
  }
  else
  {
    string errMsg = "Error when trying to open heavy quark particle file: " + filename_heavyQuarkParticleFile;
    throw eInitialModel_error( errMsg );
  }
  
  countHeavyQuarks.close();
}


void initialModel_heavyQuarks::populateParticleVector( std::vector< Particle >& _particles )
{
  // delete all previously sampled heavy quarks
  for(unsigned int j = 0; j < _particles.size(); j++ )
  {
    // if this flavor is not active in the current run delete it or convert it to gluon for PYTHIA case to obtain the same initial energy density
    if( _particles[j].FLAVOR > 2 * Particle::max_N_light_flavor && _particles[j].FLAVOR <= 2 * ( Particle::max_N_light_flavor + Particle::N_heavy_flavor ) )
    {
      // delete last particle if also not active otherwise switch position with particle to be deleted
      while( _particles.back().FLAVOR > 2 * Particle::max_N_light_flavor && _particles.back().FLAVOR <= 2 * ( Particle::max_N_light_flavor + Particle::N_heavy_flavor ) )
      {
        _particles.pop_back();
      }
      _particles[j] = _particles.back();
      _particles.pop_back();
    }
  }

  // sample momenta and position of new heavy quarks
  sampleMomenta( _particles );
  samplePositions( _particles );
}



void initialModel_heavyQuarks::sampleMomenta( std::vector< Particle >& _particles )
{
  std::ifstream readHeavyQuarks( filename_heavyQuarkParticleFile.c_str() );
  cout << "Read particle momentum of heavy quark data file." << endl;
  
  int N_EVENT_tmp, HARD_tmp, FLAVOR_tmp;
  double E_tmp, PX_tmp, PY_tmp, PZ_tmp, MASS_tmp;
  
  int numberOfHeavyQuarks = 0, newNumberOfHeavyQuarks;
  const int nmbOfPartonsWithoutCharmQuarks = _particles.size();
  
  bool charm_mass_change_cout = false, bottom_mass_change_cout = false;
  
  for ( int i = 0; i < numberOfParticlesToGenerate; i++ )
  {
    if( seperateHeavyQuarkInitialModel == pythia )
    {
      // structure of file
      // number of pythia event  is hard?   flavour  energy  px  py  pz  m
      readHeavyQuarks >> N_EVENT_tmp >> HARD_tmp >> FLAVOR_tmp >> E_tmp >> PX_tmp >> PY_tmp >> PZ_tmp >> MASS_tmp;
    }
    else
    {
      string errMsg = "Error. The requested model for the generation of initial heavy quark distributions has not yet been implemented.";
      throw eInitialModel_error( errMsg );
    }

    if( FLAVOR_tmp > 2 * Particle::max_N_light_flavor && FLAVOR_tmp <= 2 * ( Particle::max_N_light_flavor + Particle::N_heavy_flavor ) ) // heavy quarks
    {
      numberOfHeavyQuarks++;
      
      Particle tempParticle;
      tempParticle.FLAVOR = static_cast<FLAVOR_TYPE>( FLAVOR_tmp );;
      tempParticle.Mom = VectorEPxPyPz(0,PX_tmp,PY_tmp,PZ_tmp);
      tempParticle.N_EVENT_pp = N_EVENT_tmp;
      tempParticle.HARD = HARD_tmp;
      
      if( ( ( tempParticle.FLAVOR == charm || tempParticle.FLAVOR == anti_charm ) && Particle::Mcharm != MASS_tmp ) ||
          ( ( tempParticle.FLAVOR == bottom || tempParticle.FLAVOR == anti_bottom ) && Particle::Mbottom != MASS_tmp ) )
      {
        const double M_old = MASS_tmp;
        double M_new;
        if( tempParticle.FLAVOR == charm || tempParticle.FLAVOR == anti_charm )
        {
          M_new = Particle::Mcharm;
          if( !charm_mass_change_cout )
          {
            cout << "Make charm quarks lighter. Old mass=" << M_old << "  new mass=" << M_new << endl;
            charm_mass_change_cout = true;
          }
        }
        else if( tempParticle.FLAVOR == bottom || tempParticle.FLAVOR == anti_bottom )
        {
          M_new = Particle::Mbottom;
          if( !bottom_mass_change_cout )
          {
            cout << "Make bottom quarks lighter. Old mass=" << M_old << "  new mass=" << M_new << endl;
            bottom_mass_change_cout = true;
          }
        }
        if(M_new > M_old)
          cout << "problem in heavyQuarkInitialDistribution::sampleMomenta(), heavy quark mass to high. Leads to nan GeV energy." << endl;

        const double pp_old=sqrt( tempParticle.Mom.vec2() );
        const double pp_new=sqrt(pow(pp_old,2.0)+pow(M_old,2.0)-pow(M_new,2.0));

        // scaling in order to conserve the energy when changing quark mass
        tempParticle.Mom *= pp_new / pp_old;
        tempParticle.m = M_new;
      }
      else
      {
        tempParticle.m = MASS_tmp;
      }
      
      // to avoid rounding errors compute energy from momenta and mass E^2=p^2+m^2, for light partons this can be different from Pythia's value if in Pythia this parton had a mass
      tempParticle.Mom.E() = sqrt( tempParticle.Mom.vec2() + pow(tempParticle.m,2.0) );
  
      _particles.push_back( tempParticle );
    }
  }
  
  // if number of testparticles in heavy quark data file is larger than overall testparticle number, discard some of them
  if ( int( numberOfTestparticles_heavyQuarkParticleFile ) > numberOfTestparticles_BAMPS )
  {
    newNumberOfHeavyQuarks = int( double( numberOfHeavyQuarks ) / numberOfTestparticles_heavyQuarkParticleFile * numberOfTestparticles_BAMPS );
    if ( newNumberOfHeavyQuarks % 2 > 0 ) // newNumberOfHeavyQuarks is odd, but charm quarks always in pairs
    {
      newNumberOfHeavyQuarks++;
    } 

    cout << "NumberOfHeavyQuarks=" << numberOfHeavyQuarks << "  newNumberOfHeavyQuarks=" << newNumberOfHeavyQuarks << "  without rounding=" << double( numberOfHeavyQuarks ) / numberOfTestparticles_heavyQuarkParticleFile * numberOfTestparticles_BAMPS << endl;
    cout << "number before: " << _particles.size();
    cout << "    number after deletion should be: " << _particles.size() - ( numberOfHeavyQuarks - newNumberOfHeavyQuarks );
    for ( int i = 0; i < ( numberOfHeavyQuarks - newNumberOfHeavyQuarks ); i++ )
    {
      _particles.pop_back();
    }
    cout << "    number after deletion: " << _particles.size() << endl;
  }
  // if number of testparticles in pythia data file is less than overall testparticle number, duplicate some of them
  else if ( int( numberOfTestparticles_heavyQuarkParticleFile ) < numberOfTestparticles_BAMPS )
  {
    newNumberOfHeavyQuarks = int( double( numberOfHeavyQuarks ) / numberOfTestparticles_heavyQuarkParticleFile * numberOfTestparticles_BAMPS );
    if ( newNumberOfHeavyQuarks % 2 > 0 ) // newNumberOfHeavyQuarks is odd, but charm quarks always in pairs
    {
      newNumberOfHeavyQuarks++;
    }

    cout << "NumberOfCharmQuarks=" << numberOfHeavyQuarks << "  newNumberOfHeavyQuarks=" << newNumberOfHeavyQuarks << "  without rounding=" << double( numberOfHeavyQuarks ) / numberOfTestparticles_heavyQuarkParticleFile * numberOfTestparticles_BAMPS << endl;
    cout << "number before: " << _particles.size();

    int j = 1;
    for ( int i = 1; i <= ( newNumberOfHeavyQuarks - numberOfHeavyQuarks ) ;i++ ) 
    {
      Particle tempParticle = _particles[nmbOfPartonsWithoutCharmQuarks+j-1];
      _particles.push_back( tempParticle );

      j++;
      if ( j > numberOfHeavyQuarks )
      {
        j = 1;
      }
    }

    cout << "    number after duplication: " << _particles.size() << endl;
  }
}


void initialModel_heavyQuarks::samplePositions( std::vector< Particle >& _particles )
{
  cout << "Start sampling of particle positions for heavy quarks." << endl;
  int event_tmp = 0;
  VectorTXYZ Pos_tmp;


  // sample positions for binary collisions/events where hard partons are produced. Sample also if there are also soft particles at this event (reason: scaling behavior of soft and hard partons differ)
  for ( unsigned int j = 0; j < _particles.size(); j++ )
  {
    if ( _particles[j].FLAVOR > 2 * Particle::max_N_light_flavor && _particles[j].FLAVOR <= 2 * ( Particle::max_N_light_flavor + Particle::N_heavy_flavor ) )
    {
      if ( _particles[j].N_EVENT_pp != event_tmp ) // first hard particle of an event
      {
        sample_TXYZ_singleParticle( _particles[j] );
        Pos_tmp =  _particles[j].Pos;
        event_tmp = _particles[j].N_EVENT_pp;
      }
      else // all other hard particles of same event get same positions
      {
        _particles[j].Pos = Pos_tmp;
      }
    }
  }
}
