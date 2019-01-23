//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/initialmodel_CYM.cpp $
//$LastChangedDate: 2018-12-21 20:20:59 +0100 (Fr, 21. Dez 2018) $
//$LastChangedRevision: 2913 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <list>
#include <algorithm> // for std::count
#include "lorentz.h"
#include <stdlib.h>



#include <sstream>
#include <list>

#include "configuration.h"
#include "initialmodel_CYM.h"
#include "random.h"
#include "particle.h"
#include "FPT_compare.h"

using std::cout;
using std::endl;
using namespace ns_casc;


initialModel_CYM::initialModel_CYM( const config& _config ) :
  initialModelWS(_config),
  numberOfTestparticles( _config.getTestparticles() ),
  numberOfParticlesToGenerate( 0 )
{
  impactParameter = _config.getImpactParameter();
  sqrtS_perNN = _config.getSqrtS();
  filename_ParticleFile = _config.getPythiaParticleFile();
  transCellSize = _config.getTransCellSize();
  nEtaBins = _config.getNEtaBinsInitial();
  minEta = _config.getMinimumEtaInitial();
  maxEta = _config.getMaximumEtaInitial();
  randomiseAzimuth = _config.doRandomisationAzimuth();
  homogeniseSpace  = _config.doHomogenSpace();
  InitialTau0             = _config.getTau0();
  cout << "CYM Initial model: " << filename_ParticleFile << endl;
  cout << "           tau_0 = " << InitialTau0 << endl;
  
//   cout << "WOOD SAXON:" << endl;
//   if (!WoodSaxonParameter.Calculate( A, impactParameter, sqrtS_perNN))
//   {
//     std::string errMsg = "Impact parameter b too large. b > 2 R_A0";
//     throw eSimple_error( errMsg );
//   }
//   
//   cout << "======= Generating data sets for sampling of initial state =======" << endl;
//   generateTimeDistributionWS(Tab);
//   cout << "++++  Tab = " << Tab << "1/mb" << endl;
//   cout << "==================================================================" << endl;

//   std::ifstream countParticles( filename_ParticleFile.c_str() );
//   if ( countParticles.good() )
//   {
//     numberOfParticlesToGenerate = std::count( std::istreambuf_iterator<char>( countParticles ), std::istreambuf_iterator<char>(), '\n' ); // counts number of lines in  particle data file
//   }
//   else
//   {
//     string errMsg = "Error at opening particle data file.";
//     throw eSimple_error( errMsg );
//   }
//   cout << "Number of Particles in Datafile: " << numberOfParticlesToGenerate << endl;
//   countParticles.close();
}

initialModel_Attractor::initialModel_Attractor( const config& _config ) :
  initialModelWS(_config),
  numberOfTestparticles( _config.getTestparticles() ),
  numberOfParticlesToGenerate( 0 )
{
  impactParameter = _config.getImpactParameter();
  sqrtS_perNN = _config.getSqrtS();
  filename_ParticleFile = _config.getPythiaParticleFile();
  transCellSize = _config.getTransCellSize();
  nEtaBins = _config.getNEtaBinsInitial();
  minEta = _config.getMinimumEtaInitial();
  maxEta = _config.getMaximumEtaInitial();
  randomiseAzimuth = _config.doRandomisationAzimuth();
  homogeniseSpace  = _config.doHomogenSpace();
  initialAttractorCondition = _config.getInitialAttractorCondition();
  cout << "CYM Initial model: " << filename_ParticleFile << endl;
  
//   cout << "WOOD SAXON:" << endl;
//   if (!WoodSaxonParameter.Calculate( A, impactParameter, sqrtS_perNN))
//   {
//     std::string errMsg = "Impact parameter b too large. b > 2 R_A0";
//     throw eSimple_error( errMsg );
//   }
//   
//   cout << "======= Generating data sets for sampling of initial state =======" << endl;
//   generateTimeDistributionWS(Tab);
//   cout << "++++  Tab = " << Tab << "1/mb" << endl;
//   cout << "==================================================================" << endl;

//   std::ifstream countParticles( filename_ParticleFile.c_str() );
//   if ( countParticles.good() )
//   {
//     numberOfParticlesToGenerate = std::count( std::istreambuf_iterator<char>( countParticles ), std::istreambuf_iterator<char>(), '\n' ); // counts number of lines in  particle data file
//   }
//   else
//   {
//     string errMsg = "Error at opening particle data file.";
//     throw eSimple_error( errMsg );
//   }
//   cout << "Number of Particles in Datafile: " << numberOfParticlesToGenerate << endl;
//   countParticles.close();
}

void initialModel_CYM::sampleMomenta( std::vector< Particle >& _particles )
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

void initialModel_CYM::populateParticleVector( std::vector< Particle >& _particles )
{
  numberOfParticlesToGenerate = 10000000;
  
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

  
  
  int TotalNumberParticles;

  TotalNumberParticles = PositionsAndMomenta( _particles );
  
  //Instead, used for some TEST
  //TotalNumberParticles = PositionsAndMomentaThermal( _particles );
  
  cout << "Total Number of Particles from input data: " << TotalNumberParticles <<endl;
  //WARNING
  _particles.resize(TotalNumberParticles);
  
  
  for ( unsigned int j = 0; j < particles.size(); j++ )
  {
    if(particles[j].Mom.E() <=  0. || particles[j].Pos.T() <= 0.)
    {
      particles[j].dead = true;      
    }
  }
  //   sampleMomenta( _particles );
  //   samplePositions( _particles );
}

void initialModel_Attractor::populateParticleVector( std::vector< Particle >& _particles )
{
  numberOfParticlesToGenerate = 30000000;
  
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

  
  
  int TotalNumberParticles;

  
  
  switch(initialAttractorCondition)
  {
    case 0:
      TotalNumberParticles =   PositionsAndMomentaThermalAngularDisturbed( _particles ,  0.1 );
    break;
    case 1:
      TotalNumberParticles =   PositionsAndMomentaThermalAngularDisturbed( _particles ,  0.2 );
    break;
    case 2:
      TotalNumberParticles =   PositionsAndMomentaThermalAngularDisturbed( _particles ,  0.3 );
      break;
    case 3:
      TotalNumberParticles =   PositionsAndMomentaThermalAngularDisturbed( _particles ,  0.4 );
    break;
    case 4:
      TotalNumberParticles =   PositionsAndMomentaThermalAngularDisturbed( _particles ,  0.5 ); 
    break;
    case 5:
      TotalNumberParticles =   PositionsAndMomentaThermalAngularDisturbed( _particles ,  0.6 );
    break;
    case 6:
      TotalNumberParticles =   PositionsAndMomentaThermalAngularDisturbed( _particles ,  0.7 );
    break;
    case 7:
      TotalNumberParticles =   PositionsAndMomentaThermalAngularDisturbed( _particles ,  0.8 );
    break;
    case 8:
      TotalNumberParticles =   PositionsAndMomentaThermalAngularDisturbed( _particles ,  0.9 ); 
    break;    
    case 9:
      TotalNumberParticles =   PositionsAndMomentaThermalAngularDisturbed( _particles ,  1.0 ); 
      break; 
    case 10:
      TotalNumberParticles = PositionsAndMomentaThermal( _particles ); 
    break;      
    default:
      TotalNumberParticles = PositionsAndMomentaThermal( _particles );
    break;
  }
  
  
 
  
  
  
  
  
  
  
  
  
  cout << "Total Number of Particles from input data: " << TotalNumberParticles <<endl;
  
  _particles.resize(TotalNumberParticles);
  
  for ( unsigned int j = 0; j < particles.size(); j++ )
  {
    if(particles[j].Mom.E() <=  0. || particles[j].Pos.T() <= 0.)
    {
      particles[j].dead = true;      
    }
  }

}


int initialModel_CYM::PositionsAndMomenta( std::vector< Particle >& _particles )
{
  //WARNING: Needs to be set correctly!
  double lattice_spacing_a = 0.02*4.; //fm
  bool randomAzimuth = randomiseAzimuth;
  double DividerParticleNumber = 1;
  int NumberSpacingsX = 32;
  int NumberSpacingsY = 32;
  double transverse_X =0;
  double transverse_Y =0;
  double tau0 = InitialTau0; // 0.2005;
  double transverse_px =0.;
  double transverse_py =0.;
  double p_T =0;
  double minimumEta = minEta;
  double maximumEta = maxEta;
  double Eta_Sampling_bins = nEtaBins;
  double start_x=-(NumberSpacingsX*lattice_spacing_a/2);
  double start_y=-(NumberSpacingsY*lattice_spacing_a/2);
  
  cout << minimumEta << "\t" << maximumEta << "\t" << Eta_Sampling_bins << "\t" << DividerParticleNumber << endl;
  double volume  = transCellSize*transCellSize* tau0*abs(sinh( maximumEta )-sinh( minimumEta ) );
  double density = 10./volume/numberOfTestparticles;
  double cs= 10*0.1; //fm^2
  cout << "Mean free path if cross section was 10 mb: " << 1./(cs*density) << " fm" << endl;
  cout << "Cell size over this mean free path: " << transCellSize*cs*density << " < 2 " << endl;
  
  double dEta = (maximumEta-minimumEta)/Eta_Sampling_bins;
  double y=0; //Rapidity
  int particleIndex=0;
  for( unsigned int etaBin = 0; etaBin < Eta_Sampling_bins; etaBin++ )
  {
    double minimumEta_Bin = minimumEta+etaBin*dEta;
    double maximumEta_Bin = minimumEta+etaBin*dEta+dEta;
    
    
    //cout << "Eta: " << etaBin << endl;
    for ( unsigned int xBin = 0; xBin < NumberSpacingsX; xBin++ )
    {
      for ( unsigned int yBin = 0; yBin < NumberSpacingsY; yBin++ )
      {
        int positionID = xBin*(NumberSpacingsY)+yBin;
        //loadMomenta(positionID)
        std::stringstream ss;
        ss << positionID;
        string datafile = filename_ParticleFile +  "_" + ss.str() + ".dat";
        //Count entries
        std::ifstream countParticles( datafile.c_str() );
        if ( countParticles.good() )
        {
          numberOfParticlesToGenerate = std::count( std::istreambuf_iterator<char>( countParticles ), std::istreambuf_iterator<char>(), '\n' ); // counts number of lines in  particle data file  
          
          numberOfParticlesToGenerate = int(numberOfParticlesToGenerate/DividerParticleNumber);
          cout << "x = " << start_x+lattice_spacing_a*xBin << "   N= " << numberOfParticlesToGenerate << endl;
          
        }
        else
        {
          string errMsg = "Error at opening particle data file.";
          throw eCYM_error( errMsg );
        }
        countParticles.close();
        //cout << "      positionID=" << positionID << " , Datafile = " + datafile + "   Number of Particles in Datafile: " << numberOfParticlesToGenerate << endl; 
        //read in
        std::ifstream readParticles( datafile.c_str() );
        //cout << "Read particle momentum from data file." << endl;
        for ( int n = 0; n < numberOfParticlesToGenerate; n++ )
        {

          readParticles >>  transverse_px >> transverse_py ;
          //cout << flavTemp << endl;
          
          p_T=sqrt(pow(transverse_px,2.0)+pow(transverse_py,2.0));
          
          if(randomiseAzimuth)
          {
            double angle = 2.0*M_PI*ran2();
            double pxBef = transverse_px;
            double pyBef = transverse_py;
            
            transverse_px = pxBef*cos(angle) + pyBef*sin(angle);
            transverse_py = -pxBef*sin(angle) + pyBef*cos(angle); 
          }

          //homogenoeus sampling in Eta within one eta-bin
          y=minimumEta_Bin+ran2()*dEta;   
          _particles[particleIndex+n].Mom.Pz() = p_T * sinh( y );
          _particles[particleIndex+n].Mom.Px() = transverse_px;
          _particles[particleIndex+n].Mom.Py() = transverse_py;
          _particles[particleIndex+n].Mom.E() = sqrt( _particles[particleIndex+n].Mom.vec2() );
          //cout << _particles[particleIndex+n].Mom.E()/(p_T*cosh(y)) << endl;
          _particles[particleIndex+n].m = 0.0;
          _particles[particleIndex+n].FLAVOR = gluon;
          
          
          if(homogeniseSpace)
          {
            sampleHomogeneousPositions(start_x,start_x+lattice_spacing_a*(NumberSpacingsX),start_y,start_y+lattice_spacing_a*(NumberSpacingsY),transverse_X,transverse_Y );
          }
          else
          {
            sampleHomogeneousPositions(start_x+lattice_spacing_a*xBin,start_x+lattice_spacing_a*(xBin+1),start_y+lattice_spacing_a*yBin,start_y+lattice_spacing_a*(yBin+1),transverse_X,transverse_Y );
          }

          //cout << start_x+lattice_spacing_a*xBin << " - " << start_x+lattice_spacing_a*(xBin+1) << "  tX = " << transverse_X <<  endl;
          _particles[particleIndex+n].Pos=VectorTXYZ(tau0 * cosh(y) ,transverse_X,transverse_Y,tau0 * sinh( y ));
//           if(transverse_X >=  1.0 )
//           {
//             cout << "          Particle: " << particleIndex+n << "  Particle momentum: " << _particles[particleIndex+n].Mom << '\t' << "Particle position: " << _particles[particleIndex+n].Pos << endl; 
//           }
          //cout << "          Particle: " << particleIndex+n << "  Particle momentum: " << _particles[n].Mom << '\t' << "Particle position: " << _particles[n].Pos << endl;
        }
        particleIndex +=numberOfParticlesToGenerate;
      }
    }
  } 
  return particleIndex;
}

int initialModel_Attractor::PositionsAndMomentaThermal( std::vector< Particle >& _particles )
{
  //WARNING: Needs to be set correctly!
  double lattice_spacing_a = 0.02*4.; //fm
  
  
  
  bool randomAzimuth = randomiseAzimuth;
  double DividerParticleNumber = 1;
  int NumberSpacingsX = 32;
  int NumberSpacingsY = 32;
  
  //TEST
  NumberSpacingsX = 2; //50;
  NumberSpacingsY = 2; //50;
  lattice_spacing_a = 0.5; //0.02;
  
  
  double transverse_X =0;
  double transverse_Y =0;
  
  //WARNING
  double tau0 = 0.4; //0.2005;
  
  
  double transverse_px =0.;
  double transverse_py =0.;
  double p_T =0;
  double minimumEta = minEta;
  double maximumEta = maxEta;
  double Eta_Sampling_bins = nEtaBins;
  double start_x=-(NumberSpacingsX*lattice_spacing_a/2);
  double start_y=-(NumberSpacingsY*lattice_spacing_a/2);
  double dEta = (maximumEta-minimumEta)/Eta_Sampling_bins;
  cout << minimumEta << "\t" << maximumEta << "\t" << Eta_Sampling_bins << "\t" << DividerParticleNumber  << "\t" << dEta << endl;
  
//   double volume1  = transCellSize*transCellSize* tau0*abs(sinh( maximumEta )-sinh( minimumEta ) );
//    double volume =  transCellSize*transCellSize* tau0* ( tanh( maximumEta  ) - tanh( minimumEta ) );
//  double volumeGallmeister = transCellSize*transCellSize* tau0* (maximumEta-minimumEta)*NumberSpacingsX*NumberSpacingsY;
//   cout << volume1 << "\t" << volume << endl;
  
//   double density = 10./volume/numberOfTestparticles;
//   double cs= 10*0.1; //fm^2
//   cout << "Mean free path if cross section was 10 mb: " << 1./(cs*density) << " fm" << endl;
//   cout << "Cell size over this mean free path: " << transCellSize*cs*density << " < 2 " << endl;
  
  
  double y=0; //Rapidity
  double eta,etaDistributed=0;
  int particleIndex=0;
  
  double pp,f_E,costheta,ratio,phi,sintheta,E,p_E,T,gamma;
  T = 0.5; //GeV
//   numberOfParticlesToGenerate = 100;
   
   //int gluonNumber = round(  16 * pow ( T, 3.0 ) / pow ( M_PI, 2.0 ) * volumeGallmeister / pow ( 0.197, 3.0 ) * numberOfTestparticles );
//   cout << "number per cell = " << gluonNumber << "   TotNumber = " << gluonNumber*NumberSpacingsX*NumberSpacingsY << endl;
//   
   //numberOfParticlesToGenerate=gluonNumber;
  
  //ENERGY DENSITY = 19.6337 GeV^4
  
  //numberOfParticlesToGenerate=int(100000/NumberSpacingsX/NumberSpacingsY);
  
  
 
  int NumberRapBins = 1001;
  double dRap = (maximumEta-minimumEta)/NumberRapBins;
  cout << "MUST be 9 = " << transCellSize*transCellSize*NumberSpacingsX*NumberSpacingsY << endl;
  
  //Overall particle index
  particleIndex=0;
  for(unsigned int rapBin = 0; rapBin < NumberRapBins; rapBin++)
  {
      
  //     cout << dRap <<  "\t" << transCellSize*transCellSize*NumberSpacingsX*NumberSpacingsY*tau0 << "\t" << sinh( eta+dRap/2. ) - sinh( eta-dRap/2. ) << "\t" <<  endl;
      eta=minimumEta + rapBin * dRap + dRap/2.;
      
      
      
      {
        double t=tau0 * cosh( eta );
        double z=tau0 * sinh( eta );
        double etaNew=0.5*log( (t+z)/(t-z) );
        
        cout << "eta: " << eta << " etaNew= " << etaNew << "\t" << " tau0 * sinh( eta ) = z : " <<  tau0 * sinh( eta ) << endl;
      }
    
    
    
    
    
    gamma=cosh(eta);
    double volumeHere = transCellSize*transCellSize*NumberSpacingsX*NumberSpacingsY*tau0*  (  sinh( eta+dRap/2. )-  sinh( eta-dRap/2. )  );
    int gluonNumber = round(  16. * pow ( T, 3.0 ) / pow ( M_PI, 2.0 ) * volumeHere / pow ( 0.197, 3.0 ) * numberOfTestparticles );
    numberOfParticlesToGenerate=gluonNumber;
    
    if(fabs(eta) < 0.001)
      cout << "CENTRAL DENSITY = " << gluonNumber/volumeHere/numberOfTestparticles << " 1/fm^3   = " << gluonNumber/volumeHere/numberOfTestparticles*pow ( 0.197, 3.0 ) << "  GeV^3 " << endl;
    
//     cout << rapBin << "\t" << eta-dRap/2. << "\t" << eta+dRap/2.  << "   " << volumeHere << " fm^3   " << gluonNumber  << " DENSITY " << gluonNumber/volumeHere/numberOfTestparticles << " 1/fm^3   = " << gluonNumber/volumeHere/numberOfTestparticles*pow ( 0.197, 3.0 ) << "  GeV^3 " << endl;
    
    //Generate numberOfParticlesToGenerate particles at Eta

    double pzSum=0.;
    
    
    for ( int n = 0; n < numberOfParticlesToGenerate; n++ )
    {
      double X,Y;
      sampleHomogeneousPositions( start_x,start_x+ transCellSize*NumberSpacingsX,start_y,start_y+ transCellSize*NumberSpacingsY,X,Y);
      etaDistributed=eta-dRap/2. + ran2() * dRap;
      _particles[particleIndex].Pos=VectorTXYZ(tau0 * cosh( etaDistributed ),X,Y,tau0 * sinh( etaDistributed ));
      //CHECK
      double t=_particles[particleIndex].Pos.T();
      double z=_particles[particleIndex].Pos.Z();
      double etaNew=0.5*log( (t+z)/(t-z) );
//       cout << t << "\t" << z << "\t" << eta << endl;
//       cout << eta << " etaNew: " <<  etaNew << endl;
      
      //Sample LRF
      do{
        E = - 2.0 * T * log ( ran2() ); //sample area and obtain corresponding E -> combined in one step
        f_E = 4.0 / T * exp ( -E / ( 2.0 * T ) ); //comparison function evaluated at E
        p_E = E * sqrt ( pow ( E , 2.0 ) ) / ( 2.0 * pow ( T , 3.0 ) ) * exp ( -E / T );
        ratio = p_E / f_E;
      }
      while ( ran2() > ratio ); //accept if random value [0,f(E)] is in [0,p(E)]

      pp = sqrt ( pow ( E , 2 )  );
      phi = 2 * M_PI * ran2();
      costheta = 2.0 * ran2() - 1.0;
      sintheta = sqrt ( 1.0 - pow ( costheta, 2.0 ) );

      _particles[particleIndex].Mom.E()  = E;
      _particles[particleIndex].Mom.Pz() = pp * costheta;                  
      _particles[particleIndex].Mom.Py() = pp * sintheta * sin ( phi );
      _particles[particleIndex].Mom.Px() = pp * sintheta * cos ( phi );
      
      pzSum+=pow(_particles[particleIndex].Mom.Pz(),2.0)/E;
      
      p_T=sqrt(pow(particles[particleIndex].Mom.Px(),2.0)+pow(particles[particleIndex].Mom.Py(),2.0));
      _particles[particleIndex].m = 0.0;
      _particles[particleIndex].FLAVOR = gluon;
      //           //Boost
      
      VectorTXYZ boostvector;
      double beta;
      if(eta>0.00000000000000000000)
      {
        beta=-sqrt(pow(gamma,2.0)-1.)/gamma;
      }else
      {
        beta=sqrt(pow(gamma,2.0)-1.)/gamma;
      }
      
      boostvector.SetTXYZ(1.,0,0,beta);      
      lorentz L( boostvector );
      
      
      
//       lorentz L( VectorEPxPyPz( 1.0, 0.0, 0.0, - tanh(eta) ) );
      _particles[particleIndex].Mom = L.boost( _particles[particleIndex].Mom );
      particleIndex++;
    }
    cout << "ETA=" <<eta << "\tParticleNo=" << numberOfParticlesToGenerate << "\tPZSUM=" << (pzSum/volumeHere)*(pow(0.197,3.0))/numberOfTestparticles << endl;
  } 
    
//     for ( unsigned int xBin = 0; xBin < NumberSpacingsX; xBin++ )
//     {
//       for ( unsigned int yBin = 0; yBin < NumberSpacingsY; yBin++ )
//       {
//         for ( int n = 0; n < numberOfParticlesToGenerate; n++ )
//         {
//           sampleHomogeneousPositions(start_x+lattice_spacing_a*xBin,start_x+lattice_spacing_a*(xBin+1),start_y+lattice_spacing_a*yBin,start_y+lattice_spacing_a*(yBin+1),transverse_X,transverse_Y );
//           //eta=eta-dRap/2. + ran2() * dRap;
//           _particles[particleIndex].Pos=VectorTXYZ(tau0 * cosh( eta ),transverse_X,transverse_Y,tau0 * sinh( eta ));
//           //gamma=cosh(eta);
//         
//           //Sample LRF
//           do{
//             E = - 2.0 * T * log ( ran2() ); //sample area and obtain corresponding E -> combined in one step
//             f_E = 4.0 / T * exp ( -E / ( 2.0 * T ) ); //comparison function evaluated at E
//             p_E = E * sqrt ( pow ( E , 2.0 ) ) / ( 2.0 * pow ( T , 3.0 ) ) * exp ( -E / T );
//             ratio = p_E / f_E;
//           }
//           while ( ran2() > ratio ); //accept if random value [0,f(E)] is in [0,p(E)]
// 
//           pp = sqrt ( pow ( E , 2 )  );
//           phi = 2 * M_PI * ran2();
//           costheta = 2.0 * ran2() - 1.0;
//           sintheta = sqrt ( 1.0 - pow ( costheta, 2.0 ) );
// 
//           _particles[particleIndex].Mom.E()  = E;
//           _particles[particleIndex].Mom.Pz() = pp * costheta;                  
//           _particles[particleIndex].Mom.Py() = pp * sintheta * sin ( phi );
//           _particles[particleIndex].Mom.Px() = pp * sintheta * cos ( phi );
//         
//           //Boost
// //           lorentz L( VectorEPxPyPz( 1.0, 0.0, 0.0, - tanh(eta) ) );
// //           _particles[particleIndex].Mom = L.boost( _particles[particleIndex].Mom );
//           
//           p_T=sqrt(pow(particles[particleIndex].Mom.Px(),2.0)+pow(particles[particleIndex].Mom.Py(),2.0));
// 
//           _particles[particleIndex].m = 0.0;
//           _particles[particleIndex].FLAVOR = gluon;
//         
//           particleIndex++;
//         }
//       }
//     }
 
    
  
   
   
//   for ( unsigned int xBin = 0; xBin < NumberSpacingsX; xBin++ )
//   {
//     for ( unsigned int yBin = 0; yBin < NumberSpacingsY; yBin++ )
//     {
//       for ( int n = 0; n < numberOfParticlesToGenerate; n++ )
//       {
//         
//         sampleHomogeneousPositions(start_x+lattice_spacing_a*xBin,start_x+lattice_spacing_a*(xBin+1),start_y+lattice_spacing_a*yBin,start_y+lattice_spacing_a*(yBin+1),transverse_X,transverse_Y );
//         eta=minimumEta + ran2() * dEta;
//         _particles[particleIndex].Pos=VectorTXYZ(tau0 * cosh( eta ),transverse_X,transverse_Y,tau0 * sinh( eta ));
//         gamma=cosh(eta);
//         
//         
//         do{
//           E = - 2.0 * T * log ( ran2() ); //sample area and obtain corresponding E -> combined in one step
//           f_E = 4.0 / T * exp ( -E / ( 2.0 * T ) ); //comparison function evaluated at E
//           p_E = E * sqrt ( pow ( E , 2.0 ) ) / ( 2.0 * pow ( T , 3.0 ) ) * exp ( -E / T );
//           ratio = p_E / f_E;
//         }
//         while ( ran2() > ratio ); //accept if random value [0,f(E)] is in [0,p(E)]
// 
//         pp = sqrt ( pow ( E , 2 )  );
//         phi = 2 * M_PI * ran2();
//         costheta = 2.0 * ran2() - 1.0;
//         sintheta = sqrt ( 1.0 - pow ( costheta, 2.0 ) );
// 
//         _particles[particleIndex].Mom.E()  = E;
//         _particles[particleIndex].Mom.Pz() = pp * costheta;
//         
//         
//         _particles[particleIndex].Mom.Py() = pp * sintheta * sin ( phi );
//         _particles[particleIndex].Mom.Px() = pp * sintheta * cos ( phi );
//         
//          //TEST
//          lorentz L( VectorEPxPyPz( 1.0, 0.0, 0.0, - tanh(eta) ) );
//          _particles[particleIndex].Mom = L.boost( _particles[particleIndex].Mom );
//         
//          p_T=sqrt(pow(particles[particleIndex].Mom.Px(),2.0)+pow(particles[particleIndex].Mom.Py(),2.0));
//          
// //          _particles[particleIndex].Mom.Pz() = p_T * sinh( eta );
// //          _particles[particleIndex].Mom.E()=sqrt( _particles[particleIndex].Mom.vec2() );
// /*        _particles[particleIndex].Mom.Pz()=pp * costheta;
//         _particles[particleIndex].Mom.E()=E;
//         
//         double boostbeta = ((eta > 0) ? -1 : ((eta < 0) ? 1 : 0))* sqrt(1.0-1.0/(gamma*gamma));
//         double E_transformed,pz_transformed;
//         //BOOST:
//         E_transformed = _particles[particleIndex].Mom.E()*gamma- (boostbeta*gamma*_particles[particleIndex].Mom.Pz());
//         pz_transformed= (-boostbeta*gamma*_particles[particleIndex].Mom.E()) + (gamma*_particles[particleIndex].Mom.Pz());
//             */   
//         
// //         _particles[particleIndex].Mom.Pz() = pz_transformed;
// //         _particles[particleIndex].Mom.E() = E_transformed;
// //         
// //         _particles[particleIndex].Mom.Px() = sqrt(pow(E_transformed,2.0)-pow(pz_transformed,2.0)) *  cos ( phi );
// //         _particles[particleIndex].Mom.Py() = sqrt(pow(E_transformed,2.0)-pow(pz_transformed,2.0)) *  sin ( phi );
// //    
// //         
//         
// //         p_T=sqrt(pow(particles[particleIndex].Mom.Px(),2.0)+pow(particles[particleIndex].Mom.Py(),2.0));
//         
//         //homogenoeus sampling in Eta within one eta-bin
// //         y=minimumEta + ran2() * dEta;   
// //         _particles[particleIndex].Mom.Pz() = p_T * sinh( y );
// 
// //         _particles[particleIndex].Mom.E() = sqrt( _particles[particleIndex].Mom.vec2() );
//         //cout << _particles[particleIndex+n].Mom.E()/(p_T*cosh(y)) << endl;
//         _particles[particleIndex].m = 0.0;
//         _particles[particleIndex].FLAVOR = gluon;
//       
// 
// //         cout << _particles[particleIndex].Mom << endl;
//         particleIndex++;
//         
//        
//       }
//     }
//   }
  cout << particleIndex << endl;
  particleIndex--;
  return particleIndex;
}  
  
 
double sampleHomogeneous(double minimum, double maximum)
{
  return minimum+ran2()*(maximum-minimum);
} 
 
double getNormalDistro(double x, double width)
{
  return 1./sqrt(2*M_PI*pow(width,2.0))*exp(-pow(x,2.0)/(2*pow(width,2.0)));
}

double getNormalDistroMaximum(double width)
{
  return getNormalDistro(0., width);
}
 
 
int initialModel_Attractor::PositionsAndMomentaThermalAngularDisturbed( std::vector< Particle >& _particles, double width )
{

  bool randomAzimuth = randomiseAzimuth;
  double DividerParticleNumber = 1;
  
  //TEST
  int NumberSpacingsX = 2; //50;
  int NumberSpacingsY = 2; //50;
  double lattice_spacing_a = 0.5; //0.02;
  
  
  double transverse_X =0;
  double transverse_Y =0;
  
  //WARNING
  double tau0 = 0.4; //0.2005;
  
  double transverse_px =0.;
  double transverse_py =0.;
  double p_T =0;
  double minimumEta = minEta;
  double maximumEta = maxEta;
  double Eta_Sampling_bins = nEtaBins;
  double start_x=-(NumberSpacingsX*lattice_spacing_a/2);
  double start_y=-(NumberSpacingsY*lattice_spacing_a/2);
  double dEta = (maximumEta-minimumEta)/Eta_Sampling_bins;
  cout << minimumEta << "\t" << maximumEta << "\t" << Eta_Sampling_bins << "\t" << DividerParticleNumber  << "\t" << dEta << endl;
  
  double y=0; //Rapidity
  double eta,etaDistributed=0;
  int particleIndex=0;
  
  double pp,f_E,costheta,ratio,phi,sintheta,E,p_E,T,gamma;
  T = 0.5; //GeV

  int NumberRapBins = 1001;
  double dRap = (maximumEta-minimumEta)/NumberRapBins;
  cout << "MUST be 9 = " << transCellSize*transCellSize*NumberSpacingsX*NumberSpacingsY << endl;
  
  //Overall particle index
  particleIndex=0;
  for(unsigned int rapBin = 0; rapBin < NumberRapBins; rapBin++)
  {  
  //     cout << dRap <<  "\t" << transCellSize*transCellSize*NumberSpacingsX*NumberSpacingsY*tau0 << "\t" << sinh( eta+dRap/2. ) - sinh( eta-dRap/2. ) << "\t" <<  endl;
      eta=minimumEta + rapBin * dRap + dRap/2.;
      {
        double t=tau0 * cosh( eta );
        double z=tau0 * sinh( eta );
        double etaNew=0.5*log( (t+z)/(t-z) );        
        cout << "eta: " << eta << " etaNew= " << etaNew << "\t" << " tau0 * sinh( eta ) = z : " <<  tau0 * sinh( eta ) << endl;
      }

    gamma=cosh(eta);
    double volumeHere = transCellSize*transCellSize*NumberSpacingsX*NumberSpacingsY*tau0*  (  sinh( eta+dRap/2. )-  sinh( eta-dRap/2. )  );
    int gluonNumber = round(  16. * pow ( T, 3.0 ) / pow ( M_PI, 2.0 ) * volumeHere / pow ( 0.197, 3.0 ) * numberOfTestparticles );
    numberOfParticlesToGenerate=gluonNumber;
    
    if(fabs(eta) < 0.001)
      cout << "CENTRAL DENSITY = " << gluonNumber/volumeHere/numberOfTestparticles << " 1/fm^3   = " << gluonNumber/volumeHere/numberOfTestparticles*pow ( 0.197, 3.0 ) << "  GeV^3 " << endl;
    
    //cout << rapBin << "\t" << eta-dRap/2. << "\t" << eta+dRap/2.  << "   " << volumeHere << " fm^3   " << gluonNumber  << " DENSITY " << gluonNumber/volumeHere/numberOfTestparticles << " 1/fm^3   = " << gluonNumber/volumeHere/numberOfTestparticles*pow ( 0.197, 3.0 ) << "  GeV^3 " << endl;
    
    //Generate numberOfParticlesToGenerate particles at Eta
    double pzSum=0.;
    
    for ( int n = 0; n < numberOfParticlesToGenerate; n++ )
    {
      double X,Y;
      sampleHomogeneousPositions( start_x,start_x+ transCellSize*NumberSpacingsX,start_y,start_y+ transCellSize*NumberSpacingsY,X,Y);
      etaDistributed=eta-dRap/2. + ran2() * dRap;
      _particles[particleIndex].Pos=VectorTXYZ(tau0 * cosh( etaDistributed ),X,Y,tau0 * sinh( etaDistributed ));
      //CHECK
      double t=_particles[particleIndex].Pos.T();
      double z=_particles[particleIndex].Pos.Z();
      double etaNew=0.5*log( (t+z)/(t-z) );
//       cout << t << "\t" << z << "\t" << eta << endl;
//       cout << eta << " etaNew: " <<  etaNew << endl;
      //Sample Gaussian weighted cosTheta
      double cosThetaGaussian;
      while(1)
      {
        cosThetaGaussian = sampleHomogeneous(-1,1);
        if(  ran2() <=  getNormalDistro(cosThetaGaussian,width)/getNormalDistroMaximum(width) )  
        {
          break;
        } 
      }
      
      //Sample LRF
      do{
        E = - 2.0 * T * log ( ran2() ); //sample area and obtain corresponding E -> combined in one step
        f_E = 4.0 / T * exp ( -E / ( 2.0 * T ) ); //comparison function evaluated at E
        p_E = E * sqrt ( pow ( E , 2.0 ) ) / ( 2.0 * pow ( T , 3.0 ) ) * exp ( -E / T );
        ratio = p_E / f_E;
      }
      while ( ran2() > ratio ); //accept if random value [0,f(E)] is in [0,p(E)]

      pp = sqrt ( pow ( E , 2 )  );
      phi = 2 * M_PI * ran2();
      costheta = cosThetaGaussian;
      sintheta = sqrt ( 1.0 - pow ( costheta, 2.0 ) );

      _particles[particleIndex].Mom.E()  = E;
      _particles[particleIndex].Mom.Pz() = pp * costheta;                  
      _particles[particleIndex].Mom.Py() = pp * sintheta * sin ( phi );
      _particles[particleIndex].Mom.Px() = pp * sintheta * cos ( phi );
      
      pzSum+=pow(_particles[particleIndex].Mom.Pz(),2.0)/E;
      
      p_T=sqrt(pow(particles[particleIndex].Mom.Px(),2.0)+pow(particles[particleIndex].Mom.Py(),2.0));
      _particles[particleIndex].m = 0.0;
      _particles[particleIndex].FLAVOR = gluon;
      //           //Boost
      
      VectorTXYZ boostvector;
      double beta;
      if(eta>0.00000000000000000000)
      {
        beta=-sqrt(pow(gamma,2.0)-1.)/gamma;
      }else
      {
        beta=sqrt(pow(gamma,2.0)-1.)/gamma;
      }
      
      boostvector.SetTXYZ(1.,0,0,beta);      
      lorentz L( boostvector );
      
      
      
//       lorentz L( VectorEPxPyPz( 1.0, 0.0, 0.0, - tanh(eta) ) );
      _particles[particleIndex].Mom = L.boost( _particles[particleIndex].Mom );
      particleIndex++;
    }
    cout << "ETA=" <<eta << "\tParticleNo=" << numberOfParticlesToGenerate << "\tPZSUM=" << (pzSum/volumeHere)*(pow(0.197,3.0))/numberOfTestparticles << endl;
  } 
    
//     for ( unsigned int xBin = 0; xBin < NumberSpacingsX; xBin++ )
//     {
//       for ( unsigned int yBin = 0; yBin < NumberSpacingsY; yBin++ )
//       {
//         for ( int n = 0; n < numberOfParticlesToGenerate; n++ )
//         {
//           sampleHomogeneousPositions(start_x+lattice_spacing_a*xBin,start_x+lattice_spacing_a*(xBin+1),start_y+lattice_spacing_a*yBin,start_y+lattice_spacing_a*(yBin+1),transverse_X,transverse_Y );
//           //eta=eta-dRap/2. + ran2() * dRap;
//           _particles[particleIndex].Pos=VectorTXYZ(tau0 * cosh( eta ),transverse_X,transverse_Y,tau0 * sinh( eta ));
//           //gamma=cosh(eta);
//         
//           //Sample LRF
//           do{
//             E = - 2.0 * T * log ( ran2() ); //sample area and obtain corresponding E -> combined in one step
//             f_E = 4.0 / T * exp ( -E / ( 2.0 * T ) ); //comparison function evaluated at E
//             p_E = E * sqrt ( pow ( E , 2.0 ) ) / ( 2.0 * pow ( T , 3.0 ) ) * exp ( -E / T );
//             ratio = p_E / f_E;
//           }
//           while ( ran2() > ratio ); //accept if random value [0,f(E)] is in [0,p(E)]
// 
//           pp = sqrt ( pow ( E , 2 )  );
//           phi = 2 * M_PI * ran2();
//           costheta = 2.0 * ran2() - 1.0;
//           sintheta = sqrt ( 1.0 - pow ( costheta, 2.0 ) );
// 
//           _particles[particleIndex].Mom.E()  = E;
//           _particles[particleIndex].Mom.Pz() = pp * costheta;                  
//           _particles[particleIndex].Mom.Py() = pp * sintheta * sin ( phi );
//           _particles[particleIndex].Mom.Px() = pp * sintheta * cos ( phi );
//         
//           //Boost
// //           lorentz L( VectorEPxPyPz( 1.0, 0.0, 0.0, - tanh(eta) ) );
// //           _particles[particleIndex].Mom = L.boost( _particles[particleIndex].Mom );
//           
//           p_T=sqrt(pow(particles[particleIndex].Mom.Px(),2.0)+pow(particles[particleIndex].Mom.Py(),2.0));
// 
//           _particles[particleIndex].m = 0.0;
//           _particles[particleIndex].FLAVOR = gluon;
//         
//           particleIndex++;
//         }
//       }
//     }
 
    
  
   
   
//   for ( unsigned int xBin = 0; xBin < NumberSpacingsX; xBin++ )
//   {
//     for ( unsigned int yBin = 0; yBin < NumberSpacingsY; yBin++ )
//     {
//       for ( int n = 0; n < numberOfParticlesToGenerate; n++ )
//       {
//         
//         sampleHomogeneousPositions(start_x+lattice_spacing_a*xBin,start_x+lattice_spacing_a*(xBin+1),start_y+lattice_spacing_a*yBin,start_y+lattice_spacing_a*(yBin+1),transverse_X,transverse_Y );
//         eta=minimumEta + ran2() * dEta;
//         _particles[particleIndex].Pos=VectorTXYZ(tau0 * cosh( eta ),transverse_X,transverse_Y,tau0 * sinh( eta ));
//         gamma=cosh(eta);
//         
//         
//         do{
//           E = - 2.0 * T * log ( ran2() ); //sample area and obtain corresponding E -> combined in one step
//           f_E = 4.0 / T * exp ( -E / ( 2.0 * T ) ); //comparison function evaluated at E
//           p_E = E * sqrt ( pow ( E , 2.0 ) ) / ( 2.0 * pow ( T , 3.0 ) ) * exp ( -E / T );
//           ratio = p_E / f_E;
//         }
//         while ( ran2() > ratio ); //accept if random value [0,f(E)] is in [0,p(E)]
// 
//         pp = sqrt ( pow ( E , 2 )  );
//         phi = 2 * M_PI * ran2();
//         costheta = 2.0 * ran2() - 1.0;
//         sintheta = sqrt ( 1.0 - pow ( costheta, 2.0 ) );
// 
//         _particles[particleIndex].Mom.E()  = E;
//         _particles[particleIndex].Mom.Pz() = pp * costheta;
//         
//         
//         _particles[particleIndex].Mom.Py() = pp * sintheta * sin ( phi );
//         _particles[particleIndex].Mom.Px() = pp * sintheta * cos ( phi );
//         
//          //TEST
//          lorentz L( VectorEPxPyPz( 1.0, 0.0, 0.0, - tanh(eta) ) );
//          _particles[particleIndex].Mom = L.boost( _particles[particleIndex].Mom );
//         
//          p_T=sqrt(pow(particles[particleIndex].Mom.Px(),2.0)+pow(particles[particleIndex].Mom.Py(),2.0));
//          
// //          _particles[particleIndex].Mom.Pz() = p_T * sinh( eta );
// //          _particles[particleIndex].Mom.E()=sqrt( _particles[particleIndex].Mom.vec2() );
// /*        _particles[particleIndex].Mom.Pz()=pp * costheta;
//         _particles[particleIndex].Mom.E()=E;
//         
//         double boostbeta = ((eta > 0) ? -1 : ((eta < 0) ? 1 : 0))* sqrt(1.0-1.0/(gamma*gamma));
//         double E_transformed,pz_transformed;
//         //BOOST:
//         E_transformed = _particles[particleIndex].Mom.E()*gamma- (boostbeta*gamma*_particles[particleIndex].Mom.Pz());
//         pz_transformed= (-boostbeta*gamma*_particles[particleIndex].Mom.E()) + (gamma*_particles[particleIndex].Mom.Pz());
//             */   
//         
// //         _particles[particleIndex].Mom.Pz() = pz_transformed;
// //         _particles[particleIndex].Mom.E() = E_transformed;
// //         
// //         _particles[particleIndex].Mom.Px() = sqrt(pow(E_transformed,2.0)-pow(pz_transformed,2.0)) *  cos ( phi );
// //         _particles[particleIndex].Mom.Py() = sqrt(pow(E_transformed,2.0)-pow(pz_transformed,2.0)) *  sin ( phi );
// //    
// //         
//         
// //         p_T=sqrt(pow(particles[particleIndex].Mom.Px(),2.0)+pow(particles[particleIndex].Mom.Py(),2.0));
//         
//         //homogenoeus sampling in Eta within one eta-bin
// //         y=minimumEta + ran2() * dEta;   
// //         _particles[particleIndex].Mom.Pz() = p_T * sinh( y );
// 
// //         _particles[particleIndex].Mom.E() = sqrt( _particles[particleIndex].Mom.vec2() );
//         //cout << _particles[particleIndex+n].Mom.E()/(p_T*cosh(y)) << endl;
//         _particles[particleIndex].m = 0.0;
//         _particles[particleIndex].FLAVOR = gluon;
//       
// 
// //         cout << _particles[particleIndex].Mom << endl;
//         particleIndex++;
//         
//        
//       }
//     }
//   }
  cout << particleIndex << endl;
  particleIndex--;
  return particleIndex;
}  
 
 
 
 
int initialModel_Attractor::PositionsAndMomentaBox( std::vector< Particle >& _particles, double  maxMomentum, double T )
{
  //WARNING: Needs to be set correctly!
  double lattice_spacing_a = 0.02*4.; //fm
  
  
  
  bool randomAzimuth = randomiseAzimuth;
  double DividerParticleNumber = 1;
  int NumberSpacingsX = 32;
  int NumberSpacingsY = 32;
  
  //TEST
  NumberSpacingsX = 6;
  NumberSpacingsY = 6;
  lattice_spacing_a = 0.5;
  
  
  double transverse_X =0;
  double transverse_Y =0;
  double tau0 = 0.2005;
  double transverse_px =0.;
  double transverse_py =0.;
  double p_T =0;
  double minimumEta = minEta;
  double maximumEta = maxEta;
  double Eta_Sampling_bins = nEtaBins;
  double start_x=-(NumberSpacingsX*lattice_spacing_a/2);
  double start_y=-(NumberSpacingsY*lattice_spacing_a/2);
  double dEta = (maximumEta-minimumEta)/Eta_Sampling_bins;
  cout << minimumEta << "\t" << maximumEta << "\t" << Eta_Sampling_bins << "\t" << DividerParticleNumber  << "\t" << dEta << endl;
  
//   double volume1  = transCellSize*transCellSize* tau0*abs(sinh( maximumEta )-sinh( minimumEta ) );
//    double volume =  transCellSize*transCellSize* tau0* ( tanh( maximumEta  ) - tanh( minimumEta ) );
  double volumeGallmeister = transCellSize*transCellSize* tau0* (maximumEta-minimumEta)*NumberSpacingsX*NumberSpacingsY;
//   cout << volume1 << "\t" << volume << endl;
  
//   double density = 10./volume/numberOfTestparticles;
//   double cs= 10*0.1; //fm^2
//   cout << "Mean free path if cross section was 10 mb: " << 1./(cs*density) << " fm" << endl;
//   cout << "Cell size over this mean free path: " << transCellSize*cs*density << " < 2 " << endl;
  
  
  double y=0; //Rapidity
  double eta=0;
  int particleIndex=0;
  
   double pp,f_E,costheta,ratio,phi,sintheta,E,p_E,gamma;
//   T = 1.0;
//   numberOfParticlesToGenerate = 100;
   
   //int gluonNumber = round(  16 * pow ( T, 3.0 ) / pow ( M_PI, 2.0 ) * volumeGallmeister / pow ( 0.197, 3.0 ) * numberOfTestparticles );
//   cout << "number per cell = " << gluonNumber << "   TotNumber = " << gluonNumber*NumberSpacingsX*NumberSpacingsY << endl;
//   
   //numberOfParticlesToGenerate=gluonNumber;
  
  //ENERGY DENSITY = 19.6337 GeV^4
  
  //numberOfParticlesToGenerate=int(100000/NumberSpacingsX/NumberSpacingsY);
  
  
 
  int NumberRapBins = 251;
  double dRap = (maximumEta-minimumEta)/NumberRapBins;
  cout << "MUST be 9 = " << transCellSize*transCellSize*NumberSpacingsX*NumberSpacingsY << endl;
  
  //Overall particle index
  particleIndex=0;
  for(unsigned int rapBin = 0; rapBin < NumberRapBins; rapBin++)
  {
      
  //     cout << dRap <<  "\t" << transCellSize*transCellSize*NumberSpacingsX*NumberSpacingsY*tau0 << "\t" << sinh( eta+dRap/2. ) - sinh( eta-dRap/2. ) << "\t" <<  endl;
      eta=minimumEta + rapBin * dRap + dRap/2.;
      
      
      
      {
      double t=tau0 * cosh( eta );
      double z=tau0 * sinh( eta );
      double etaNew=0.5*log( (t+z)/(t-z) );
      
      cout << "eta: " << eta << " etaNew= " << etaNew << "\t" << " tau0 * sinh( eta ) = z : " <<  tau0 * sinh( eta ) << endl;
    }
    
    
    
    
    
    gamma=cosh(eta);
    double volumeHere = transCellSize*transCellSize*NumberSpacingsX*NumberSpacingsY*tau0*  (  sinh( eta+dRap/2. )-  sinh( eta-dRap/2. )  );
    int gluonNumber = round(gamma*  16. * pow ( T, 3.0 ) / pow ( M_PI, 2.0 ) * volumeHere / pow ( 0.197, 3.0 ) * numberOfTestparticles );
    numberOfParticlesToGenerate=gluonNumber;
    
    if(fabs(eta) < 0.001)
      cout << "CENTRAL DENSITY = " << gluonNumber/volumeHere/numberOfTestparticles << " 1/fm^3   = " << gluonNumber/volumeHere/numberOfTestparticles*pow ( 0.197, 3.0 ) << "  GeV^3 " << endl;
    
//     cout << rapBin << "\t" << eta-dRap/2. << "\t" << eta+dRap/2.  << "   " << volumeHere << " fm^3   " << gluonNumber  << " DENSITY " << gluonNumber/volumeHere/numberOfTestparticles << " 1/fm^3   = " << gluonNumber/volumeHere/numberOfTestparticles*pow ( 0.197, 3.0 ) << "  GeV^3 " << endl;
    
    //Generate numberOfParticlesToGenerate particles at Eta

    
    for ( int n = 0; n < numberOfParticlesToGenerate; n++ )
    {
      double X,Y;
      sampleHomogeneousPositions( start_x,start_x+ transCellSize*NumberSpacingsX,start_y,start_y+ transCellSize*NumberSpacingsY,X,Y);
      eta=eta-dRap/2. + ran2() * dRap;
      _particles[particleIndex].Pos=VectorTXYZ(tau0 * cosh( eta ),X,Y,tau0 * sinh( eta ));
     
      E=maxMomentum*ran2(); 
      pp = sqrt ( pow ( E , 2 )  );
      
      
      phi = 2 * M_PI * ran2();
      costheta = 2.0 * ran2() - 1.0;
      sintheta = sqrt ( 1.0 - pow ( costheta, 2.0 ) );

      _particles[particleIndex].Mom.E()  = E;
      _particles[particleIndex].Mom.Pz() = pp * costheta;                  
      _particles[particleIndex].Mom.Py() = pp * sintheta * sin ( phi );
      _particles[particleIndex].Mom.Px() = pp * sintheta * cos ( phi );
      p_T=sqrt(pow(particles[particleIndex].Mom.Px(),2.0)+pow(particles[particleIndex].Mom.Py(),2.0));
      _particles[particleIndex].m = 0.0;
      _particles[particleIndex].FLAVOR = gluon;
      //           //Boost
      lorentz L( VectorEPxPyPz( 1.0, 0.0, 0.0, - tanh(eta) ) );
      _particles[particleIndex].Mom = L.boost( _particles[particleIndex].Mom );
      particleIndex++;
    }
    cout << eta << "\t" << numberOfParticlesToGenerate << endl;
  } 
    
//     for ( unsigned int xBin = 0; xBin < NumberSpacingsX; xBin++ )
//     {
//       for ( unsigned int yBin = 0; yBin < NumberSpacingsY; yBin++ )
//       {
//         for ( int n = 0; n < numberOfParticlesToGenerate; n++ )
//         {
//           sampleHomogeneousPositions(start_x+lattice_spacing_a*xBin,start_x+lattice_spacing_a*(xBin+1),start_y+lattice_spacing_a*yBin,start_y+lattice_spacing_a*(yBin+1),transverse_X,transverse_Y );
//           //eta=eta-dRap/2. + ran2() * dRap;
//           _particles[particleIndex].Pos=VectorTXYZ(tau0 * cosh( eta ),transverse_X,transverse_Y,tau0 * sinh( eta ));
//           //gamma=cosh(eta);
//         
//           //Sample LRF
//           do{
//             E = - 2.0 * T * log ( ran2() ); //sample area and obtain corresponding E -> combined in one step
//             f_E = 4.0 / T * exp ( -E / ( 2.0 * T ) ); //comparison function evaluated at E
//             p_E = E * sqrt ( pow ( E , 2.0 ) ) / ( 2.0 * pow ( T , 3.0 ) ) * exp ( -E / T );
//             ratio = p_E / f_E;
//           }
//           while ( ran2() > ratio ); //accept if random value [0,f(E)] is in [0,p(E)]
// 
//           pp = sqrt ( pow ( E , 2 )  );
//           phi = 2 * M_PI * ran2();
//           costheta = 2.0 * ran2() - 1.0;
//           sintheta = sqrt ( 1.0 - pow ( costheta, 2.0 ) );
// 
//           _particles[particleIndex].Mom.E()  = E;
//           _particles[particleIndex].Mom.Pz() = pp * costheta;                  
//           _particles[particleIndex].Mom.Py() = pp * sintheta * sin ( phi );
//           _particles[particleIndex].Mom.Px() = pp * sintheta * cos ( phi );
//         
//           //Boost
// //           lorentz L( VectorEPxPyPz( 1.0, 0.0, 0.0, - tanh(eta) ) );
// //           _particles[particleIndex].Mom = L.boost( _particles[particleIndex].Mom );
//           
//           p_T=sqrt(pow(particles[particleIndex].Mom.Px(),2.0)+pow(particles[particleIndex].Mom.Py(),2.0));
// 
//           _particles[particleIndex].m = 0.0;
//           _particles[particleIndex].FLAVOR = gluon;
//         
//           particleIndex++;
//         }
//       }
//     }
 
    
  
   
   
//   for ( unsigned int xBin = 0; xBin < NumberSpacingsX; xBin++ )
//   {
//     for ( unsigned int yBin = 0; yBin < NumberSpacingsY; yBin++ )
//     {
//       for ( int n = 0; n < numberOfParticlesToGenerate; n++ )
//       {
//         
//         sampleHomogeneousPositions(start_x+lattice_spacing_a*xBin,start_x+lattice_spacing_a*(xBin+1),start_y+lattice_spacing_a*yBin,start_y+lattice_spacing_a*(yBin+1),transverse_X,transverse_Y );
//         eta=minimumEta + ran2() * dEta;
//         _particles[particleIndex].Pos=VectorTXYZ(tau0 * cosh( eta ),transverse_X,transverse_Y,tau0 * sinh( eta ));
//         gamma=cosh(eta);
//         
//         
//         do{
//           E = - 2.0 * T * log ( ran2() ); //sample area and obtain corresponding E -> combined in one step
//           f_E = 4.0 / T * exp ( -E / ( 2.0 * T ) ); //comparison function evaluated at E
//           p_E = E * sqrt ( pow ( E , 2.0 ) ) / ( 2.0 * pow ( T , 3.0 ) ) * exp ( -E / T );
//           ratio = p_E / f_E;
//         }
//         while ( ran2() > ratio ); //accept if random value [0,f(E)] is in [0,p(E)]
// 
//         pp = sqrt ( pow ( E , 2 )  );
//         phi = 2 * M_PI * ran2();
//         costheta = 2.0 * ran2() - 1.0;
//         sintheta = sqrt ( 1.0 - pow ( costheta, 2.0 ) );
// 
//         _particles[particleIndex].Mom.E()  = E;
//         _particles[particleIndex].Mom.Pz() = pp * costheta;
//         
//         
//         _particles[particleIndex].Mom.Py() = pp * sintheta * sin ( phi );
//         _particles[particleIndex].Mom.Px() = pp * sintheta * cos ( phi );
//         
//          //TEST
//          lorentz L( VectorEPxPyPz( 1.0, 0.0, 0.0, - tanh(eta) ) );
//          _particles[particleIndex].Mom = L.boost( _particles[particleIndex].Mom );
//         
//          p_T=sqrt(pow(particles[particleIndex].Mom.Px(),2.0)+pow(particles[particleIndex].Mom.Py(),2.0));
//          
// //          _particles[particleIndex].Mom.Pz() = p_T * sinh( eta );
// //          _particles[particleIndex].Mom.E()=sqrt( _particles[particleIndex].Mom.vec2() );
// /*        _particles[particleIndex].Mom.Pz()=pp * costheta;
//         _particles[particleIndex].Mom.E()=E;
//         
//         double boostbeta = ((eta > 0) ? -1 : ((eta < 0) ? 1 : 0))* sqrt(1.0-1.0/(gamma*gamma));
//         double E_transformed,pz_transformed;
//         //BOOST:
//         E_transformed = _particles[particleIndex].Mom.E()*gamma- (boostbeta*gamma*_particles[particleIndex].Mom.Pz());
//         pz_transformed= (-boostbeta*gamma*_particles[particleIndex].Mom.E()) + (gamma*_particles[particleIndex].Mom.Pz());
//             */   
//         
// //         _particles[particleIndex].Mom.Pz() = pz_transformed;
// //         _particles[particleIndex].Mom.E() = E_transformed;
// //         
// //         _particles[particleIndex].Mom.Px() = sqrt(pow(E_transformed,2.0)-pow(pz_transformed,2.0)) *  cos ( phi );
// //         _particles[particleIndex].Mom.Py() = sqrt(pow(E_transformed,2.0)-pow(pz_transformed,2.0)) *  sin ( phi );
// //    
// //         
//         
// //         p_T=sqrt(pow(particles[particleIndex].Mom.Px(),2.0)+pow(particles[particleIndex].Mom.Py(),2.0));
//         
//         //homogenoeus sampling in Eta within one eta-bin
// //         y=minimumEta + ran2() * dEta;   
// //         _particles[particleIndex].Mom.Pz() = p_T * sinh( y );
// 
// //         _particles[particleIndex].Mom.E() = sqrt( _particles[particleIndex].Mom.vec2() );
//         //cout << _particles[particleIndex+n].Mom.E()/(p_T*cosh(y)) << endl;
//         _particles[particleIndex].m = 0.0;
//         _particles[particleIndex].FLAVOR = gluon;
//       
// 
// //         cout << _particles[particleIndex].Mom << endl;
//         particleIndex++;
//         
//        
//       }
//     }
//   }
  cout << particleIndex << endl;
  particleIndex--;
  return particleIndex;
}  
   
 
int initialModel_CYM::PositionsAndMomentaHomogeneousGlasmaFixEDens( std::vector< Particle >& _particles )
{
  //WARNING: Needs to be set correctly!
  double lattice_spacing_a = 0.02*4.; //fm
  
  
  
  bool randomAzimuth = randomiseAzimuth;
  double DividerParticleNumber = 1;
  int NumberSpacingsX = 32;
  int NumberSpacingsY = 32;
  
  //TEST
  NumberSpacingsX = 6;
  NumberSpacingsY = 6;
  lattice_spacing_a = 0.5;
  
  
  double transverse_X =0;
  double transverse_Y =0;
  double tau0 = 0.2005;
  double transverse_px =0.;
  double transverse_py =0.;
  double p_T =0;
  double minimumEta = minEta;
  double maximumEta = maxEta;
  double Eta_Sampling_bins = nEtaBins;
  double start_x=-(NumberSpacingsX*lattice_spacing_a/2);
  double start_y=-(NumberSpacingsY*lattice_spacing_a/2);
  double dEta = (maximumEta-minimumEta)/Eta_Sampling_bins;
  cout << minimumEta << "\t" << maximumEta << "\t" << Eta_Sampling_bins << "\t" << DividerParticleNumber  << "\t" << dEta << endl;
  
//   double volume1  = transCellSize*transCellSize* tau0*abs(sinh( maximumEta )-sinh( minimumEta ) );
  double volume =  transCellSize*transCellSize* tau0* ( tanh( maximumEta  ) - tanh( minimumEta ) );
  
//   cout << volume1 << "\t" << volume << endl;
  
  double density = 10./volume/numberOfTestparticles;
  double cs= 10*0.1; //fm^2
  cout << "Mean free path if cross section was 10 mb: " << 1./(cs*density) << " fm" << endl;
  cout << "Cell size over this mean free path: " << transCellSize*cs*density << " < 2 " << endl;
  
  
  double y=0; //Rapidity
  double eta=0;
  int particleIndex=0;
  
  double pp,f_E,costheta,ratio,phi,sintheta,E,p_E,T,gamma;
  
  numberOfParticlesToGenerate = 100;
  T = 0.5;
  int gluonNumber = round (  16 * pow ( T, 3.0 ) / pow ( M_PI, 2.0 ) * volume / pow ( 0.197, 3.0 ) * numberOfTestparticles );
  cout << "number per cell = " << gluonNumber << "   TotNumber = " << gluonNumber*NumberSpacingsX*NumberSpacingsY << endl;
  
  numberOfParticlesToGenerate=gluonNumber*NumberSpacingsX*NumberSpacingsY;
  
  
  for ( unsigned int xBin = 0; xBin < NumberSpacingsX; xBin++ )
  {
    for ( unsigned int yBin = 0; yBin < NumberSpacingsY; yBin++ )
    {
      for ( int n = 0; n < gluonNumber; n++ )
      {
        
        sampleHomogeneousPositions(start_x+lattice_spacing_a*xBin,start_x+lattice_spacing_a*(xBin+1),start_y+lattice_spacing_a*yBin,start_y+lattice_spacing_a*(yBin+1),transverse_X,transverse_Y );
        eta=minimumEta + ran2() * dEta;
        _particles[particleIndex].Pos=VectorTXYZ(tau0 * cosh( eta ),transverse_X,transverse_Y,tau0 * sinh( eta ));
        gamma=cosh(eta);
        
        
        do{
          E = - 2.0 * T * log ( ran2() ); //sample area and obtain corresponding E -> combined in one step
          f_E = 4.0 / T * exp ( -E / ( 2.0 * T ) ); //comparison function evaluated at E
          p_E = E * sqrt ( pow ( E , 2.0 ) ) / ( 2.0 * pow ( T , 3.0 ) ) * exp ( -E / T );
          ratio = p_E / f_E;
        }
        while ( ran2() > ratio ); //accept if random value [0,f(E)] is in [0,p(E)]

        pp = sqrt ( pow ( E , 2 )  );
        phi = 2 * M_PI * ran2();
        costheta = 2.0 * ran2() - 1.0;
        sintheta = sqrt ( 1.0 - pow ( costheta, 2.0 ) );

        _particles[particleIndex].Mom.Pz()=pp * costheta;
        _particles[particleIndex].Mom.E()=E;
        
        double boostbeta = ((eta > 0) ? -1 : ((eta < 0) ? 1 : 0))* sqrt(1.0-1.0/(gamma*gamma));
        double E_transformed,pz_transformed;
        //BOOST:
        E_transformed = _particles[particleIndex].Mom.E()*gamma- (boostbeta*gamma*_particles[particleIndex].Mom.Pz());
        pz_transformed= (-boostbeta*gamma*_particles[particleIndex].Mom.E()) + (gamma*_particles[particleIndex].Mom.Pz());
        
        _particles[particleIndex].Mom.Pz() = pz_transformed;
        _particles[particleIndex].Mom.E() = E_transformed;
        _particles[particleIndex].Mom.Px() = sqrt(pow(E_transformed,2.0)-pow(pz_transformed,2.0)) *  cos ( phi );
        _particles[particleIndex].Mom.Py() = sqrt(pow(E_transformed,2.0)-pow(pz_transformed,2.0)) *  sin ( phi );
   
        
        
        p_T=sqrt(pow(particles[particleIndex].Mom.Px(),2.0)+pow(particles[particleIndex].Mom.Py(),2.0));
        
        //homogenoeus sampling in Eta within one eta-bin
        y=minimumEta + ran2() * dEta;   
//         _particles[particleIndex].Mom.Pz() = p_T * sinh( y );

//         _particles[particleIndex].Mom.E() = sqrt( _particles[particleIndex].Mom.vec2() );
        //cout << _particles[particleIndex+n].Mom.E()/(p_T*cosh(y)) << endl;
        _particles[particleIndex].m = 0.0;
        _particles[particleIndex].FLAVOR = gluon;
      

//         cout << _particles[particleIndex].Mom << endl;
        particleIndex++;
        
       
      }
    }
  }
  cout << particleIndex << endl;
  particleIndex--;
  return particleIndex;
}  
   
 
void initialModel_CYM::sampleHomogeneousPositions( double posXmin, double posXmax, double posYmin, double posYmax, double &posX, double &posY )
{
  double random1 = ran2();
  double random2 = ran2();
  posX = random1*(posXmax-posXmin)+posXmin;
  posY = random2*(posYmax-posYmin)+posYmin;
}

void initialModel_Attractor::sampleHomogeneousPositions( double posXmin, double posXmax, double posYmin, double posYmax, double &posX, double &posY )
{
  double random1 = ran2();
  double random2 = ran2();
  posX = random1*(posXmax-posXmin)+posXmin;
  posY = random2*(posYmax-posYmin)+posYmin;
}

void initialModel_CYM::samplePositions( std::vector< Particle >& _particles )
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
void initialModel_CYM::sample_TXYZ_one_partcl( Particle& _particle, bool& soft )
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

double initialModel_CYM::Radius() 
{
  return WoodSaxonParameter.RA;
}

