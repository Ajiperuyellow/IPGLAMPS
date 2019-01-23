//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/initialmodel_cgc.cpp $
//$LastChangedDate: 2015-06-10 11:23:50 +0200 (水, 10  6月 2015) $
//$LastChangedRevision: 2173 $
//$LastChangedBy: senzel $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <math.h>
#include <string>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <algorithm> // for std::count

#include "configuration.h"
#include "initialmodel_cgc.h"
#include "random.h"
#include "particle.h"

using std::cout;
using std::endl;


initialModel_CGC::initialModel_CGC( const config& _config ) : 
initialModelWS(_config),
nTestparticles( _config.getTestparticles() ),
filename_cgcParticleFile( _config.getCgcParticleFile() )
{
  impactParameter = _config.getImpactParameter();
  sqrtS_perNN = _config.getSqrtS();

  if (!WoodSaxonParameter.Calculate( A, impactParameter, sqrtS_perNN))
  {
    std::string errMsg = "Impact parameter b too large. b > 2 R_A0";
    throw eInitialModel_error( errMsg );
  }

  cout << "CGC particle data file: " << filename_cgcParticleFile << endl;
  std::ifstream countCGCParticles( filename_cgcParticleFile.c_str() );
  if ( countCGCParticles.good() )
  {
    nParticlesToGenerate = 0;
    do
    {
      string line;
      getline ( countCGCParticles, line );
      if( line.length() != 0 && line.find ( "#", 0 ) == string::npos )   //ignore lines with leading "#"
      {
        nParticlesToGenerate++;
      }
    }
    while( !countCGCParticles.eof() );
  }
  else
  {
    string errMsg = "Error at opening CGC particle data file.";
    throw eInitialModel_error( errMsg );
  }

  countCGCParticles.close();
}



void initialModel_CGC::populateParticleVector( std::vector< Particle >& _particles )
{
  /**
  * Reserve memory for the Particle vector. Due to possible particle creation this size should be chosen somewhat larger
  * than the initial number of particles. Doing this, push_back operations to add new particles won't lead to internal
  * memory re-allocation in the vector, which could possibly be time consuming.
  * However push_back operations are always possible, even when the reserved size is used up (as long as there is physical
  * memory available of course). And modern day compilers probably optimize to the extent where it doesn't really matter -
  * but still, it's good practice.
  */
  _particles.reserve( static_cast<int>( nParticlesToGenerate * 1.2 ) );
  /**
  * Now the particle vector is re-sized to hold the initial number of particles (this is NOT done when reserving memory!).
  * The particles are initialised with the standard constructor, actual attributes such as momentum etc. MUST be set later!
  */
  _particles.resize( nParticlesToGenerate );

  //----
  const double tau0 = 0.1; // intitial time 1/1.4 GeV^-1 = 0.14 fm (~1.4 GeV ist saturation momentum)
  //----

  std::ifstream cgcParticles( filename_cgcParticleFile.c_str() );
  int n = 0;
  int nEvents = 0;
  
  do
  {
    string line;
    getline ( cgcParticles, line );
    
    if( line.length() == 0 )
    {
      nEvents++;
    }
    else if ( line.find ( "#", 0 ) == string::npos )   //ignore lines with leading "#"
    {
      std::vector<double> lineContent;
      std::stringstream inStream;
      inStream.str ( line );
      while ( !inStream.eof() )
      {
        double dummy;
        inStream >> dummy;
        lineContent.push_back ( dummy );
      }
      if ( lineContent.size() == 4 )
      {
        double rx = lineContent[0];
        double ry = lineContent[1];
        double y = lineContent[2];
        double pt = lineContent[3];
        
        if ( sqrt( rx*rx + ry*ry ) <= 100.0 ) // just to check if rx, ry are from actual data
        {
          n++;

          _particles[n].Pos = VectorTXYZ(tau0 * cosh( y ),rx,ry,tau0 * sinh( y ));

          double phi = 2.0 * M_PI * ran2();
          _particles[n].Mom.Pz() = pt * sinh( y );
          _particles[n].Mom.Px() = pt * cos( phi );
          _particles[n].Mom.Py() = pt * sin( phi );
          _particles[n].Mom.E() = sqrt( _particles[n].Mom.vec2() );

          _particles[n].m = 0.0;
          _particles[n].FLAVOR = gluon;
        }
      }
      else
      {
          string errMsg = "Error in initialModel_CGC::populateParticleVector - line too short or too long";
          throw eInitialModel_error ( errMsg );
      }
    }
  }
  while ( !cgcParticles.eof() );

  if ( n != nParticlesToGenerate )
  {
    string errMsg = "Error in initialModel_CGC::populateParticleVector, wrong number of particles";
    throw eInitialModel_error( errMsg );
  }
  
  if( nTestparticles != nEvents )
  {
    string errMsg = "Error in initialModel_CGC::populateParticleVector, wrong number of events in comparison with testparticles";
    throw eInitialModel_error( errMsg );
  }
}


