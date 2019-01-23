//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/initialmodel_pythia_online.cpp $
//$LastChangedDate: 2014-07-15 11:00:55 +0200 (火, 15  7月 2014) $
//$LastChangedRevision: 1798 $
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
#include "configBAMPS.h"
#include "hadroncrosssection.h"
#include "initialmodel_pythia_online.h"
#include "random.h"
#include "particle.h"
#include "FPT_compare.h"

#ifdef Pythia_FOUND
#include "Pythia8/Pythia.h"
#include "Pythia8/LHAPDFInterface.h"
#endif

using std::cout;
using std::endl;
using namespace ns_casc;


initialModel_Pythia_online::initialModel_Pythia_online( const config& _config ) :
  initialModel_Pythia(_config)
{
  if( _config.getPDFsource() == LHAPDF )
    use_lhapdf = true;
  else
    use_lhapdf = false;
  
  if( use_lhapdf )
  {
    name_pdfset = _config.getLHAPDFdatasetName();
    if( _config.getLHAPDFuseGrid() )
      name_pdfset += ".LHgrid";
    else
      name_pdfset += ".LHpdf";

    if( _config.useNuclearPDFs() )
      shadowing_parametrization_name = _config.getNuclearPDFdatasetName();
    else
      shadowing_parametrization_name = "";
  }
}

void initialModel_Pythia_online::sampleMomenta( std::vector< Particle >& _particles )
{
#ifdef Pythia_FOUND
  // give correct xml directory to pythia object, only possible via constructor
  string xmlpath = PYTHIA_XML_DIR;
  
  // Generator; shorthand for event and particleData.                           
  Pythia8::Pythia pythia( xmlpath );
  Pythia8::Event& event      = pythia.event;
  Pythia8::ParticleData& pdt = pythia.particleData;

  // Key requirement: switch off HadronLevel to prevent hadronization and only get partons
  pythia.readString("HadronLevel:all = off");
  
  // switch on soft QCD processes, which also includes min. bias processes, which also describe hard processes (bad labeling in Pythia). Form Pythia's manual concerning soft QCD processes: "As a rule, the processes in this class should not be mixed with the simulation of other processes. All by themselves, they are intended to represent the total cross section of hadron collisions, with the exception of the "rare processes" that one wishes to study separately. In particular, jet physics at all scales occurs as part of the minimum-bias description. "
  pythia.readString("SoftQCD:all = on");
  
  // suppress unnecessary output
  pythia.readString("Main:timesToShow = 0"); //                ! show how far along run is this many times
  pythia.readString("Main:showChangedSettings = off"); //      ! print changed flags/modes/parameters
  pythia.readString("Main:showChangedParticleData = off"); //  ! print changed particle and decay data
  pythia.readString("Next:numberShowEvent = 0"); //            ! suppress full listing of first events
  
  // use LHAPDF
  if( use_lhapdf )
  {
    pythia.readString("PDF:useLHAPDF = on");
    pythia.readString("PDF:LHAPDFset = " + name_pdfset);
    
    // Pythia8_bamps_link::shadowing_parametrization = shadowing_parametrization_name;
    
    // if( shadowing_parametrization_name != "" )
    // {
    //   if( A == B )
    //   {
    //     Pythia8_bamps_link::A = A;
    //   }
    //   else
    //   {
    //     std::string errMsg = "Mass number of both nucleii is not the same!";
    //     throw ePythia_error( errMsg );
    //   }
    // }
  }

  // Set the random seed of the Pythia Random Number Generator.
  // Please note: Here we use a random number genarted by the
  // BAMPS-RNG, since every call of pythia.init should produce
  // different events.
  pythia.readString("Random:setSeed = on");
  std::stringstream buffer;
  buffer << "Random:seed = " << static_cast<long>((ran2()*899999998.0)+1.0) ;
  pythia.readString(buffer.str());
  
  int proton_code = 2212;
  // Initialize.
  pythia.init( proton_code, proton_code, sqrtS_perNN );
  
  XSgenerator_pp_inelast_Pythia sigma_inelast;
  double sigma_inel_pp = sigma_inelast(sqrtS_perNN) ; // mb
  
  // Number of binary collisions times testparticles
  int N_bin = int( sigma_inel_pp * Tab * numberOfTestparticles );
  cout << "number of binary collisions: " << N_bin /numberOfTestparticles << endl;
  
  // run pythia events
  for (int iEvent = 1; iEvent <= N_bin; ++iEvent) 
  {
    if( !pythia.next() )
      continue;
    
//     pythia.event.list();
//     pythia.process.list();
//     pythia.info.list();
      
    // go through particle list of pythia event
    for( int i = 0; i < pythia.event.size(); ++i) 
    {
      if( pythia.event[i].isFinal() )
      {
        if(  pythia.event[i].isGluon() || pythia.event[i].isQuark() )
        {
          Particle tempParticle;
          tempParticle.N_EVENT_pp = iEvent;
          tempParticle.HARD = true; // regard all particles as HARD to keep all of them (soft particles are deleted for impact parameter b!=0 later since position sampling is not written for soft particles)

          tempParticle.Mom = VectorEPxPyPz( pythia.event[i].e(),
                                            pythia.event[i].px(),
                                            pythia.event[i].py(),
                                            pythia.event[i].pz() );

          tempParticle.m = pythia.event[i].m();
          tempParticle.FLAVOR = convertPythiaFlavor( pythia.event[i].id() );
    
          // light partons are massless
          if ( tempParticle.FLAVOR <= 2 * Particle::max_N_light_flavor )
          {
            tempParticle.m = 0;
          }
    
          // to avoid rounding errors compute energy from momenta and mass E^2=p^2+m^2, for light partons this can be different from Pythia's value if in Pythia this parton had a mass
          tempParticle.Mom.E() = sqrt( tempParticle.Mom.vec2() + pow( tempParticle.m, 2.0 ) );
          
          _particles.push_back( tempParticle );
        }
      }
    }
  }
#else
  std::string errMsg = "Pythia 8 not available to sample initial distributions";
  throw ePythia_error( errMsg );
#endif
  
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


FLAVOR_TYPE initialModel_Pythia_online::convertPythiaFlavor( const int _pythia_flavor_id )
{
  // quark notation    g     u     d     s     c     b     t
  // Pythia (PDG)      21    2     1     3     4     5     6      (antipartcle *(-1))
  // cascade           0     1     3     5     7     9     11     (antipartcle +1)
  
  FLAVOR_TYPE flav;
  switch ( _pythia_flavor_id )
  {
    case 21:
      flav = gluon;
      break;
    case 2:
      flav = up;
      break;
    case -2:
      flav = anti_up;
      break;
    case 1:
      flav = down;
      break;
    case -1:
      flav = anti_down;
      break;
    case 3:
      flav = strange;
      break;
    case -3:
      flav = anti_strange;
      break;
    case 4:
      flav = charm;
      break;
    case -4:
      flav = anti_charm;
      break;
    case 5:
      flav = bottom;
      break;
    case -5:
      flav = anti_bottom;
      break;    
    default:
      std::string errMsg = "Error in convertPythiaFlavor!";
      throw ePythia_error( errMsg );
      break;
  }
  
  return flav;
}
