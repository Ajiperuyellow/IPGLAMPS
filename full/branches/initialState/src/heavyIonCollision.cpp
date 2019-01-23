//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/heavyIonCollision.cpp $
//$LastChangedDate: 2018-04-19 08:34:36 +0200 (Do, 19. Apr 2018) $
//$LastChangedRevision: 2723 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <iostream>
#include <math.h>
#include <string>
#include <deque>
#include <stdlib.h>
#include <iterator>

#include "binary_cross_sections.h"
#include "configuration.h"
#include "particle.h"
#include "random.h"
#include "heavyIonCollision.h"
#include "lorentz.h"
#include "FPT_compare.h"
#include "configuration.h"
#include "coordinateBins.h"
#include "coupling.h"
#include "initialmodel.h"
#include "initialmodel_minijets.h"
#include "initialmodel_pythia_offline.h"
#include "initialmodel_pythia_online.h"
#include "initialmodel_cgc.h"
#include "initialmodel_heavyquarks.h"
#include "initialmodel_hydroParametrization.h"
#include "initialmodel_simple.h"
#include "initialmodel_CYM.h"

#include "debug.h"


using namespace std;
using namespace ns_casc;


/**
 * The constructor for the heavyIonCollision class. It sets various
 * parameters, the name of output files, etc. 
 *
 * @param[in] _theConfig Pointer to a config object. This is used to
 * obtain user provided settings from the input file 
 * @param[in] _theMFP Pointer to a mfpForHeavyIonCollision
 * object. This is needed in case the computation of the mfp for high
 * energy particles via interpolation from tabulated values is
 * requested. See ::JET_MFP_COMPUTATION_TYPE and
 * config::getJetMfpComputationType() 
 */
heavyIonCollision::heavyIonCollision( config*const _theConfig, mfpForHeavyIonCollision*const _theMFP ) : 
  theConfig( _theConfig ), theMFP( _theMFP ),
  theI23_massless( false ), theI23_charm_m1( false ), theI23_charm_m2( false ), theI23_bottom_m1( false ), theI23_bottom_m2( false ), // do not load data files right at construction, but after configure() has been called below
  rings( 10, 1.5, 1.0 ),  // 10 rings, inner ring with radius 1.5 fm, outer rings with radius 1.0 fm
  offlineInterface( _theConfig->getPathdirOfflineData().c_str(), _theConfig->doOutput_offlineReconstructionData() )
{
  offlineInterface.setAdditionalFilenameTag( theConfig->getJobName() );
  
  stop = theConfig->getRuntime();
  A = theConfig->getA();
  Aatomic = theConfig->getAatomic();
  B = theConfig->getB();
  Batomic = theConfig->getBatomic();
  sqrtS = theConfig->getSqrtS();
  P0 = theConfig->getPtCutoff();
  Bimp = theConfig->getImpactParameter();
  testpartcl = theConfig->getTestparticles();
  
  // load 2->2 cross section interpolation data
  if( theConfig->doScattering_22() )
  {
    theI22.configure( theConfig->isCouplingRunning(), Particle::N_light_flavor, Particle::N_heavy_flavor, Particle::Mcharm, Particle::Mbottom, theConfig->getMaxRunningCoupling(), theConfig->getfixedCouplingValue() );
  }
  
  // load 2->3 cross section interpolation data
  if( theConfig->doScattering_23() )
  {
    theI23_massless.configure( theConfig->I23onlineIntegrationIsSet(), 1, 0.0, theConfig->getKappa23LightPartons(), theConfig->get23GluonFormationTimeTyp(), theConfig->getMatrixElement23(), theConfig->isMd2CounterTermInI23(), theConfig->get23FudgeFactorLpm(), theConfig->getInterpolation23Mode(), theConfig->isMatrixElement23_22qt() );
    if( Particle::N_heavy_flavor > 0 )
    {
      theI23_charm_m1.configure( theConfig->I23onlineIntegrationIsSet(), 1, Particle::Mcharm, theConfig->getKappa23HeavyQuarks(), theConfig->get23GluonFormationTimeTyp(), theConfig->getMatrixElement23(), theConfig->isMd2CounterTermInI23(), theConfig->get23FudgeFactorLpm(), theConfig->getInterpolation23Mode(), theConfig->isMatrixElement23_22qt() );
      theI23_charm_m2.configure( theConfig->I23onlineIntegrationIsSet(), 2, Particle::Mcharm, theConfig->getKappa23HeavyQuarks(), theConfig->get23GluonFormationTimeTyp(), theConfig->getMatrixElement23(), theConfig->isMd2CounterTermInI23(), theConfig->get23FudgeFactorLpm(), theConfig->getInterpolation23Mode(), theConfig->isMatrixElement23_22qt() );
    }
    if( Particle::N_heavy_flavor > 1 )
    {
      theI23_bottom_m1.configure( theConfig->I23onlineIntegrationIsSet(), 1, Particle::Mbottom, theConfig->getKappa23HeavyQuarks(), theConfig->get23GluonFormationTimeTyp(), theConfig->getMatrixElement23(), theConfig->isMd2CounterTermInI23(), theConfig->get23FudgeFactorLpm(), theConfig->getInterpolation23Mode(), theConfig->isMatrixElement23_22qt() );
      theI23_bottom_m2.configure( theConfig->I23onlineIntegrationIsSet(), 2, Particle::Mbottom, theConfig->getKappa23HeavyQuarks(), theConfig->get23GluonFormationTimeTyp(), theConfig->getMatrixElement23(), theConfig->isMd2CounterTermInI23(), theConfig->get23FudgeFactorLpm(), theConfig->getInterpolation23Mode(), theConfig->isMatrixElement23_22qt() );
      
    }
  }

  timeshift = 0;
  deltaEta_fine = 0.001;
  
  nGet23Errors = 0;
  nGet32Errors = 0;

  numberOfEdgeCells = 8;
  
  
}

heavyIonCollision::~heavyIonCollision()
{
}



/**
 * This routine initializes: particle momenta, positions, masses,
 * flavors and formation times. 
 * The type of initial state model that is used is fetched from the
 * user settings via config::getInitialStateType() 
 */
void heavyIonCollision::initialize( analysis& aa )
{
  //---------------  
  bool bool_deleteParticlesFromList, bool_formationTimeForParticles;
  //---------------    
  
  //---------------    
  //needed for the time shifitng of the particles
  const double eta_max = 5.0;  
  //---------------    
  
  //---------------  
  string infoInitialParameters;//string to write the initial conditions in the output file
  //---------------  
  
  //---------------  
  initialModel *initialmodel;

  switch ( theConfig->getInitialStateType() )
  {
    case miniJetsInitialState:
      initialmodel = new initialModel_minijets( *theConfig );
      bool_deleteParticlesFromList = true;      
      bool_formationTimeForParticles = true;
      break;
    case pythiaOfflineInitialState:
      initialmodel = new initialModel_Pythia_offline( *theConfig );
      bool_deleteParticlesFromList = true;      
      bool_formationTimeForParticles = true;    
      break;
    case pythiaOnlineInitialState:
      initialmodel = new initialModel_Pythia_online( *theConfig );
      bool_deleteParticlesFromList = true;      
      bool_formationTimeForParticles = true;    
      break;
    case cgcInitialState:
      initialmodel = new initialModel_CGC( *theConfig );
      bool_deleteParticlesFromList = true;      
      bool_formationTimeForParticles = true;     
      break;
    case hydroParametrizationInitialState:
      initialmodel = new initialModel_HydroParametrization( *theConfig, infoInitialParameters );
      bool_deleteParticlesFromList = false;      
      bool_formationTimeForParticles = false;   
      break;     
    case simpleInitialState:
      initialmodel = new initialModel_simple( *theConfig );
      bool_deleteParticlesFromList = true;      
      bool_formationTimeForParticles = true;   
      break;     
    case CYMInitialState:
      initialmodel = new initialModel_CYM( *theConfig );
      bool_deleteParticlesFromList = false;      
      bool_formationTimeForParticles = false;   
      break; 
    case AttractorInitialState:
      initialmodel = new initialModel_Attractor( *theConfig );
      bool_deleteParticlesFromList = false;      
      bool_formationTimeForParticles = false;   
      break; 
    default:
      std::string errMsg = "Model for sampling the initial state not implemented yet!";
      throw eInitialModel_error( errMsg );
      break;
  }
  //---------------  
  
  //---------------
  initialmodel->populateParticleVector( particles );
  initialmodel->setUniqueID( particles );

  typicalRadius = initialmodel->Radius();
  cout << "Typical Radius: " << typicalRadius <<  endl;
  //---------------  
  
  //---------------
  aa.infoInitialParameters = infoInitialParameters;
  //---------------    
  
  //---------------
  //Special for heavy quarks
  if( theConfig->isSeperateHeavyQuarkParticleFile() )
  {
    initialModel_heavyQuarks _heavyQuarkInitialDistribution( *theConfig, pythia );
    _heavyQuarkInitialDistribution.populateParticleVector( particles );

    typicalRadius = _heavyQuarkInitialDistribution.Radius();
  }
  //---------------

  //---------------
  if(bool_deleteParticlesFromList){deleteParticlesFromList();}
  //---------------  
  
  //---------------   
  timeShiftingOfParticles(eta_max);  
  //---------------   
  
  //---------------    
  if(bool_formationTimeForParticles){formationTimeForParticles();}
  //---------------    

  //---------------   
  listParticleNumbersForAllFlavors(initialState,0.0);
  //--------------- 
  
  delete initialmodel;
}


/**
* @brief Delete Particles from list after initialisation
*/
void heavyIonCollision::deleteParticlesFromList()
{
  int Nbefore = particles.size();
  
  for(unsigned int j = 0; j < particles.size(); j++ )
  {
    // if this flavor is not active in the current run delete it or convert it to gluon for PYTHIA case to obtain the same initial energy density
    if( ( particles[j].FLAVOR > 2 * Particle::N_light_flavor ) && 
        ( ( particles[j].FLAVOR <= 2 * Particle::max_N_light_flavor ) || 
          ( particles[j].FLAVOR > 2 * ( Particle::max_N_light_flavor + Particle::N_heavy_flavor ) ) ) )
    {
      if( theConfig->getInitialStateType() == pythiaOnlineInitialState || theConfig->getInitialStateType() == pythiaOfflineInitialState )
      {
        // transform non-active quarks to gluons
        double pp = particles[j].Mom.vec2();
        double ee = pp + pow( particles[j].m , 2.0 );

        // scaling in order to conserve the energy when making quarks massles
        double fak = sqrt( ee / pp );
        particles[j].Mom *= fak;
        particles[j].Mom(0) = sqrt( particles[j].Mom.vec2() );
        particles[j].m = 0.0;
        particles[j].FLAVOR = gluon;   
      }
      else
      {
        // delete last particle if also not active otherwise switch position with particle to be deleted
        while( ( particles.back().FLAVOR > 2 * Particle::N_light_flavor ) && 
               ( ( particles.back().FLAVOR <= 2 * Particle::max_N_light_flavor ) || 
                 ( particles.back().FLAVOR > 2 * ( Particle::max_N_light_flavor + Particle::N_heavy_flavor ) ) ) &&
               ( j != particles.size() - 1  ) ) // if particle j is the last particle in the particle list it is deleted here and the then last in the list below as well, which would be wrong.
        {
          particles.pop_back();
        }
        particles[j] = particles.back();
        particles.pop_back();
      }
    }
  }
  cout << "#### " << particles.size() << " out of " << Nbefore << " particles kept for simulation ( N_f = N_f_light_quarks + N_f_heavy_quarks = " << Particle::N_light_flavor << " + " <<   Particle::N_heavy_flavor << " )." << endl;
}


/**
* @brief Time shifting of particles after initialisation 
*/
void heavyIonCollision::timeShiftingOfParticles(const double eta_max)
{
  double dt, max = 0;
  
  for ( unsigned int i = 0; i < particles.size(); i++ )
  {
    dt = fabs( particles[i].Pos.Z() ) / tanh( eta_max ) - particles[i].Pos.T();  //tanh(eta)=z/t
    if ( dt > max ){max = dt;}
  }
  
  timeshift = max;

  for ( unsigned int i = 0; i < particles.size(); i++ )
  {
    particles[i].Pos.T() += timeshift;
  }
  
  cout << endl;
  cout << "timeshift for the particles is set to " << timeshift << " fm/c" << endl;
  cout << endl;
}

/**
* @brief Formation time for particles
*
* formation time: 1/sqrt( p_T^2 + m^2) = 1/m_T
*/
void heavyIonCollision::formationTimeForParticles()
{
  double dt, y, MT;
  
  for ( unsigned int i = 0; i < particles.size(); i++ )
  {
    y = particles[i].Mom.Rapidity();
    MT = particles[i].Mom.Mt( particles[i].m );
    dt = 1 / MT * cosh( y ) * 0.197327;  //fm/c   //cosh(y) = gamma  (of that particle wrt motion in z-direction)
    if ( MT < 0.01 ) // ugly hack
    {
      dt = infinity;
    }

    particles[i].Propagate( particles[i].Pos.T() + dt );
  }
}


/**
* @brief Lists the number of particles of every flavor at different states
*/
void heavyIonCollision::listParticleNumbersForAllFlavors(const simulationType st, const double time)
{  
  // List particle numbers for all flavors
  cout << "==========================" << endl;
  cout << "List of Particles: " << endl;  
  //------------------------
  switch ( st )
    {
    case initialState:
      cout << "Initial State " << endl;
      break;
    case simulationState:
      cout << "time " << time << endl;
      break;
    case finalState:
      cout << "Final State " << endl;
      break;
    default:
      cout << "Error in simulation state!" << endl;
    }
  //------------------------  
  cout << "# of all = " << particles.size() << endl;
  
  int maxFlavor = 16;
  
  if( !theConfig->getAnalyseForHydro() )
  {
    for( int i = 0; i <= 2 * Particle::N_light_flavor; i++ )
    {
      int fl_sum = 0;
      for( unsigned int j = 0; j < particles.size(); j++ )
        if(particles[j].FLAVOR == i)
          fl_sum++;
      cout << "Flavor = " << i << "   # = " << fl_sum << endl;
    }
    for( int i = 2 * Particle::max_N_light_flavor + 1; i <= 2 * ( Particle::max_N_light_flavor + Particle::N_heavy_flavor ); i++ )
    {
      int fl_sum = 0;
      for( unsigned int j = 0; j < particles.size(); j++ )
        if(particles[j].FLAVOR == i)
          fl_sum++;
      cout << "Flavor = " << i << "   # = " << fl_sum << endl;
    }
  }
  else
  {
    for( int i = 0; i <= maxFlavor; i++ )
    {
      int fl_sum = 0;
      for( unsigned int j = 0; j < particles.size(); j++ )
        if(particles[j].FLAVOR == i)
          fl_sum++;
      cout << "Flavor = " << i << "   # = " << fl_sum << endl;
    }
  }
  
  double sum = 0.0;
  for( unsigned int j = 0; j < particles.size(); j++ )
    sum += particles[j].Mom.E();
  //------------------------
  switch ( st )
    {
    case initialState:
      cout << "E(init) = " << sum << " GeV" << endl;
      break;
    case simulationState:
      cout << "E(" << time << ") = " << sum << " GeV" << endl;
      break;
    case finalState:
      cout << "E(final) = " << sum << " GeV" << endl;
      cout << "E(final)/Testparticles = " << sum/testpartcl << " GeV" << endl;
      break;
    default:
      cout << "Error in simulation state!" << endl;
    }
  //------------------------  
  for( unsigned int j = 0; j < particles.size(); j++ )
  {
    if(particles[j].FLAVOR > maxFlavor)
    {
      cout << "There are particles with larger flavor then listed here; Please check! " << endl;
      break;
    }
  }  
  cout << "==========================" << endl;
  

}



void heavyIonCollision::prepareGeometricCollision( const int iscat, const int jscat, 
                                                   double& M1, double& M2, 
                                                   FLAVOR_TYPE& F1, FLAVOR_TYPE& F2, double& s,
                                                   double& ct_i, double& ct_j )
{
  double a, b, c, d, ee;

  M1 = particles[iscat].m;
  F1 = particles[iscat].FLAVOR;

  M2 = particles[jscat].m;
  F2 = particles[jscat].FLAVOR;

  s = (particles[iscat].Mom + particles[jscat].Mom).M2();

  a = Dot( particles[iscat].Mom, particles[jscat].Pos - particles[iscat].Pos );
  b = Dot( particles[jscat].Mom, particles[jscat].Pos - particles[iscat].Pos );

  c = M1 * M1;
  d = M2 * M2;

  ee = Dot( particles[iscat].Mom, particles[jscat].Mom );

  ct_i = particles[iscat].Pos.T() - ( a * d - b * ee ) / ( ee * ee - c * d ) * particles[iscat].Mom.E();
  ct_j = particles[jscat].Pos.T() + ( b * c - a * ee ) / ( ee * ee - c * d ) * particles[jscat].Mom.E();
}




/**
* @brief Update partons after geometrical 2->2 collision
*/
void heavyIonCollision::updateAfterGeometricCollision( const int iscat, const int jscat, const VectorEPxPyPz & P1, const VectorEPxPyPz & P2, const double ct_i, const double ct_j )
{
  double m1,m2;

  particles[iscat].Old = particles[iscat].Mom;
  m1 = particles[iscat].m;

  particles[jscat].Old = particles[jscat].Mom;
  m2 = particles[jscat].m ;

  // set new outgoing particle 1:
  particles[iscat].Propagate( ct_i );
  particles[iscat].Mom = P1;
  particles[iscat].Mom(0) = sqrt(m1*m1 + P1.vec2());

  // set new outgoing particle 2:
  particles[jscat].Propagate( ct_j );
  particles[jscat].Mom = P2;
  particles[jscat].Mom(0) = sqrt(m2*m2 + P2.vec2());
}


/**
 * @param[in] ii ID of particle 1
 * @param[in] jj ID of particle 2
 */
double heavyIonCollision::getGeometricCollisionTime( const int ii, const int jj )
{
  FLAVOR_TYPE F1, F2;
  double M1, M2;
  double s, a, b, c, d, ee, f;
  double h, AB, sd, cs, ct1, ct2, min, max, ot, md2g_wo_as, md2q_wo_as;

  if (( particles[ii].coll_id == particles[jj].coll_id ) && ( particles[ii].coll_id != 0 ) )
    return infinity;

  M1 = particles[ii].m;
  F1 = particles[ii].FLAVOR;

  M2 = particles[jj].m;
  F2 = particles[jj].FLAVOR;

  s = (particles[ii].Mom + particles[jj].Mom).M2();


  //   if(s <= (-2.0*tcut)) return infinity;//s>=-(t+u) if t<=tcut, u<=tcut

  a = Dot( particles[ii].Mom, particles[jj].Pos - particles[ii].Pos );
  b = Dot( particles[jj].Mom, particles[jj].Pos - particles[ii].Pos );

  c = M1 * M1;
  d = M2 * M2;

  ee = Dot( particles[ii].Mom, particles[jj].Mom );

  f = ( particles[jj].Pos - particles[ii].Pos ).M2();

  h = a + b;
  if ( h > 0.0 )
    AB = -c * b + ee * a;
  else
    AB = d * a - ee * b;
  if ( AB >= 0.0 )
    return infinity;

  if ((( ee*ee - c*d )/(particles[ii].Mom.E() * particles[jj].Mom.E()) ) < 1.0e-20 )
  {
    cout << "parallel ee*ee-c*d= " << ee*ee - c*d << endl;
    return infinity;
  }

  sd = -f - ( a * a * d + b * b * c - 2.0 * a * b * ee ) / ( ee * ee - c * d );//fm^2
  sd = sd * 10.0;//fm^2 -> mb
  if ( sd < 0.0 )
  {
    std::string errMsg = "sd < 0 in getGeometricCollisionTime(..)";
    throw eHIC_error( errMsg );
  }
  
  if ( particles[ii].md2g < particles[jj].md2g )
  {
    md2g_wo_as = particles[ii].md2g;
  }
  else
  {
    md2g_wo_as = particles[jj].md2g;
  }
  if ( particles[ii].md2q < particles[jj].md2q )
  {
    md2q_wo_as = particles[ii].md2q;
  }
  else
  {
    md2q_wo_as = particles[jj].md2q;
  }

  // Ugly hack. For md2g all ELASTIC total cross sections become infinite, thus the two particles will collide.
  // Needs a better solution..
  if ( md2g_wo_as != 0 )
  {
    //------ isotropic scattering -----
    //   cs=30.0;//mb
    //   cs=10.0;//mb
    //   cs=5.0;//mb
    //   cs=2.0;//mb
    //------ isotropic scattering -----
    scattering22 scatt22_object( &theI22 );
    scatt22_object.setParameter( particles[ii].Mom, particles[jj].Mom, 
                                 F1, F2, M1, M2, s, md2g_wo_as , md2q_wo_as,
                                 theConfig->getKggQQb(), theConfig->getKgQgQ(), 
                                 theConfig->getKappa_gQgQ(), 
                                 theConfig->isConstantCrossSecGQ(),
                                 theConfig->getConstantCrossSecValueGQ(), 
                                 theConfig->isIsotropicCrossSecGQ(),
                                 theConfig->getKfactor_light() ); // md2g, md2q are debye masses without the factor alpha_s which is multiplied in scattering22.cpp
    cs = scatt22_object.getXSectionElastic(); //1/GeV^2
    cs = cs * 0.3894; //1/GeV^2->mb
    cs = cs / testpartcl;
    
    if (( M_PI*sd ) > cs )
    {
      return infinity;
    }
  }

  ct1 = particles[ii].Pos.T() - ( a * d - b * ee ) / ( ee * ee - c * d ) * particles[ii].Mom.E();
  ct2 = particles[jj].Pos.T() + ( b * c - a * ee ) / ( ee * ee - c * d ) * particles[jj].Mom.E();

  if ( ct1 < ct2 )
  {
    min = ct1;
    max = ct2;
  }
  else
  {
    min = ct2;
    max = ct1;
  }

  if (( min < particles[ii].Pos.T() ) || ( min < particles[jj].Pos.T() ) || ( max > stop ) )
    return infinity;
  //   if(max > stop) return infinity;
  else
  {
    ot = min;
    //     ot=max;
    //     ot=0.5*(ct1+ct2);

    return ot;
  }
}




int heavyIonCollision::binomial( const int N, const int k ) const
{
  if ( N < k )
  {
    return 0;
  }
  else if ( k == 3 )
  {
    return ( N * ( N - 1 ) * ( N - 2 ) / 6 );
  }
  else if ( k == 2 )
  {
    return ( N * ( N - 1 ) / 2 );
  }
  else if ( k <= 0 )
  {
    return 0;
  }
  else if ( N == k )
  {
    return 1;
  }
  else
  {
    int nominator = 1;
    for ( int i = 0; i < k; i++ )
    {
      nominator *= ( N - i );
    }

    int denominator = 1;
    for ( int i = 1; i <= k; i++ )
    {
      denominator *= i;
    }

    return ( nominator / denominator );
  }

  return 0;
}



/**
 * In case config::jetMfpComputationSwitch requests that the mean free path of jet particles should be computed iteratively,
 * this routine does the iterative procedure based on the particles that are in the same cell as the jet particle.
 *
 * Not written for heavy quarks!
 *
 * @param[in] _allParticlesList List of all particles in the given cell
 * @param[in] _gluonList List of all gluons in the given cell
 * @param[in] jetID ID of the "jet" particle for which the mean free path should be iteratively computed
 * @param[in] dt Size of the current time step
 * @param[in] dv Volume of the current cell
 * @param[in] time Current time of the simulation
 */
double heavyIonCollision::iterateMFP( std::vector< int >& _allParticlesList, std::vector< int >& _gluonList, const int jetID, const double dt, const double dv, const double time )
{
  int iscat, jscat;
  int n32 = 0, n22 = 0, n23 = 0;
  double lambda_scaled, s;
  double probab22 = 0, probab23 = 0, probab32 = 0;
  double cs22, cs23, I32;
  double R22, R23, R32;
  double Vrel, md2g_wo_as, md2q_wo_as;
  double lambda, lambdaAvr;
  int iter = 0;
  double avS2x = 0, avS32 = 0;
  deque<double> lambdaArray;
  double betaDistEntry;

  int initialStateIndex = -1;
  FLAVOR_TYPE F1, F2, F3;

  const double epsilon = 0.01;
  const int nIterationsMax = 8;
  bool converged = false;


  //----------------------------- to which ring does the jet belong? ---------------------------
  int nc = rings.getIndex( particles[jetID] );
  //--------------------------------------------------------------------------------------------

  // === Particle 1: jetID
  F1 = particles[jetID].FLAVOR;
  
  const int nTotal = _allParticlesList.size();
  const int nGluons = _gluonList.size();
  int nJetGluons = 0;
  if ( particles[jetID].FLAVOR == gluon )
  {
    nJetGluons = 1; 
  }
  const int nAllQuarks = nTotal - nGluons;
  const int nJetQuarks = 1 - nJetGluons;
  
  const int allTriplets = ( binomial( nTotal, 2 ) )  - ( binomial( nAllQuarks, 2 ) * nJetQuarks );

  //------------------------------ determine initial value of lambda ---------------------------
  double csgg = 0;
  const double small = 1.0e-4;
  double jetRate = ( particles[jetID].rate22 + particles[jetID].rate23 + particles[jetID].rate32 +
                     particles[jetID].rate22v + particles[jetID].rate23v + particles[jetID].rate32v ) / 2.0;
  if ( jetRate > small )
  {
    lambda = 1 / jetRate;
  }
  else
  {
    unsigned int selected = 0;
    do  // pick random particle (that is neihter dead nor the jet) from the _particleList
    {
      selected = static_cast<int>( ran2() * _allParticlesList.size() );
      if ( selected >= _allParticlesList.size() )
      {
        continue;
      }      
      iscat = _allParticlesList[selected];
    }
    while ( particles[iscat].dead || iscat == jetID || particles[iscat].isHeavyQuark() );
    
    // === Particle 2 : iscat
    
    s = (particles[jetID].Mom + particles[iscat].Mom).M2();
    md2g_wo_as = ( particles[jetID].md2g + particles[iscat].md2g ) / 2.0;
    md2q_wo_as = ( particles[jetID].md2q + particles[iscat].md2q ) / 2.0;
    
    xsection_gg_gg csObj( s, md2g_wo_as, md2q_wo_as, &theI22, theConfig->getKfactor_light() );
    csgg = csObj.totalCrossSection();
    lambda = ( dv * rings[nc].getGamma() * testpartcl   ) / ( pow( 0.197, 3.0 ) * _allParticlesList.size() * csgg );  
  }
  //--------------------------------------------------------------------------------------------
  
  scattering22 scatt22_object( &theI22 );
  scattering32 scatt32_object;
  scattering23 scatt23_object( &theI23_massless, &theI23_charm_m1, &theI23_charm_m2, &theI23_bottom_m1, &theI23_bottom_m2 );
  
  do
  {
    //------------------------ 2<->2 & 2->3-----------------------
    for ( unsigned int m1 = 0; m1 < _allParticlesList.size(); m1++ )
    {
      jscat = _allParticlesList[m1];
      if (( jscat != jetID ) && ( !particles[jscat].dead ) && ( !particles[jscat].isHeavyQuark() ) )
      {
        // === Particle 2: jscat
        F2 = particles[jscat].FLAVOR;
        
        s = (particles[jetID].Mom + particles[jscat].Mom).M2();
        avS2x += s;
        if ( s > 1.1*lambda2 )
        {
          n22++;
          Vrel = VelRel(particles[jetID].Mom,particles[jscat].Mom, 0.0,0.0);
          
          md2g_wo_as = ( particles[jetID].md2g + particles[jscat].md2g ) / 2.0;
          md2q_wo_as = ( particles[jetID].md2q + particles[jscat].md2q ) / 2.0;
          
          if ( theConfig->doScattering_22() )
          {
            double mass = 0;
            scatt22_object.setParameter( particles[jetID].Mom, particles[jscat].Mom,
                                         F1, F2, mass, mass, s, md2g_wo_as , md2q_wo_as,
                                         theConfig->getKggQQb(), theConfig->getKgQgQ(), 
                                         theConfig->getKappa_gQgQ(), 
                                         theConfig->isConstantCrossSecGQ(),
                                         theConfig->getConstantCrossSecValueGQ(), 
                                         theConfig->isIsotropicCrossSecGQ(), 
                                         theConfig->getKfactor_light() ); // md2g_wo_as, md2q_wo_as are debye masses without the factor alpha_s which is multiplied in scattering22.cpp
            cs22 = scatt22_object.getXSection22();
            probab22 += pow( 0.197, 2.0 ) * cs22 * Vrel * dt / ( dv * testpartcl );
          }
          else
          {
            probab22 += 0;
          }
          
          n23++;
          lambda_scaled = lambda * sqrt( s );
          
          if ( theConfig->doScattering_23() )
          {
            double mass = 0;
            betaDistEntry = scatt23_object.setParameter( rings[nc].getAveraged_v(),
                                                         particles[jetID].Mom, particles[jscat].Mom,
                                                         F1, F2, mass, mass, sqrt( s ), 
                                                         md2g_wo_as / s, lambda_scaled, 
                                                         theConfig->getK23LightPartons(), 
                                                         theConfig->getK23HeavyQuarks(),
                                                         theConfig->getKappa23LightPartons(), 
                                                         theConfig->getKappa23HeavyQuarks(),
                                                         theConfig->I23onlineIntegrationIsSet(),
                                                         theConfig->get23GluonFormationTimeTyp(), 
                                                         theConfig->getMatrixElement23(), 
                                                         theConfig->isMd2CounterTermInI23(), 
                                                         theConfig->get23FudgeFactorLpm(), 
                                                         _gluonList.size() );    
            cs23 = scatt23_object.getXSection23( initialStateIndex ); //1/GeV^2
            probab23 += pow( 0.197, 2.0 ) * cs23 * Vrel * dt / ( dv * testpartcl );
          }
          else
          {
            probab23 += 0;
          }
        }
      }
    }
    R22 = probab22 / dt * rings[nc].getGamma();
    R23 = probab23 / dt * rings[nc].getGamma();
    //-------------------------------------------------------------

    
    //---------------------------- 3->2 ---------------------------
    const int consideredTriplets = 5;
    double scaleForSelectedTriplets = 1;
    if ( theConfig->doScattering_32() )
    {
      if ( allTriplets > 20 )
      {      
        unsigned int m2, m3;
        scaleForSelectedTriplets = static_cast<double>( allTriplets ) / static_cast<double>( consideredTriplets );
        
        for ( int i = 0; i < consideredTriplets && _gluonList.size() > 0; i++ )
        {
          do
          {
            do
            {
              m2 = int ( _allParticlesList.size() * ran2() );
              if ( m2 == _allParticlesList.size() )
              {
                m2 = _allParticlesList.size() - 1;
              }
              iscat = _allParticlesList[m2];
            }
            while ( particles[iscat].dead || iscat == jetID || particles[iscat].isHeavyQuark() );
            F2 = particles[iscat].FLAVOR;
            
            do
            {
              m3 = int ( _allParticlesList.size() * ran2() );
              if ( m3 == _allParticlesList.size() )
              {
                m3 = _allParticlesList.size() - 1;
              }
              jscat = _allParticlesList[m3];
            }
            while ( iscat == jscat || particles[jscat].dead || jscat == jetID || particles[jscat].isHeavyQuark() );
            F3 = particles[jscat].FLAVOR;          
          }
          while ( !( F1 == gluon || F2 == gluon || F3 == gluon ) );
          
          // === Particle 2: iscat
          // === Particle 3: jscat
          
          s = ( particles[jetID].Mom + particles[iscat].Mom + particles[jscat].Mom ).M2();
          avS32 += s;
          
          if ( s > 1.1*lambda2 )
          {
            n32++;
            lambda_scaled = lambda * sqrt( s );
            
            md2g_wo_as = ( particles[iscat].md2g + particles[jscat].md2g + particles[jetID].md2g ) / 3.0;
            md2g_wo_as = ( particles[iscat].md2q + particles[jscat].md2q + particles[jetID].md2q ) / 3.0;
            
            // create scattering32 object for the given 3 particles
            betaDistEntry = scatt32_object.setParameter( rings[nc].getAveraged_v(),
                                                         particles[jetID].Mom, particles[iscat].Mom, particles[jscat].Mom,
                                                         F1, F2, F3, 
                                                         sqrt( s ), md2g_wo_as * coupling::get_constant_coupling() / s, lambda_scaled, coupling::get_constant_coupling(), theConfig->isMd2CounterTermInI23(), theConfig->isMatrixElement23_22qt(), theConfig->get23FudgeFactorLpm(), nGluons );  // create scattering32 object for the given 3 particles
            I32 = scatt32_object.getIntegral32_withPrefactors();                        // get the integral I32 for the given 3 particles
            
            probab32 += I32 * dt / ( pow( dv, 2.0 ) * pow( testpartcl, 2.0 ) );
          }
        }
      }
      else
      {
        //---------------------------- 3->2 ---------------------------
        for ( unsigned int m1 = 0; m1 < _allParticlesList.size() - 1; m1++ )
        {
          iscat = _allParticlesList[m1];
          if ( !particles[iscat].dead && iscat != jetID && !particles[iscat].isHeavyQuark() )
          {   
            for ( unsigned int m2 = m1 + 1; m2 < _allParticlesList.size(); m2++ )
            {
              jscat = _allParticlesList[m2];
              if ( !particles[jscat].dead && jscat != jetID && !particles[jscat].isHeavyQuark() )
              {
                F2 = particles[iscat].FLAVOR;
                F3 = particles[jscat].FLAVOR;

                // at least one of the particles must be a gluon
                // otherwise go to next step in the loop
                if ( !( F1 == gluon || F2 == gluon || F3 == gluon ) )
                {
                  continue;
                }
              
                // === Particle 2: iscat
                // === Particle 3: jscat
              
                s = ( particles[jetID].Mom + particles[iscat].Mom + particles[jscat].Mom ).M2();
                avS32 += s;
              
                if ( s > 1.1*lambda2 )
                {
                  n32++;
                  lambda_scaled = lambda * sqrt( s );

                  md2g_wo_as = ( particles[iscat].md2g + particles[jscat].md2g + particles[jetID].md2g ) / 3.0;
                  md2g_wo_as = ( particles[iscat].md2q + particles[jscat].md2q + particles[jetID].md2q ) / 3.0;
                  
                  // create scattering32 object for the given 3 particles
                  betaDistEntry = scatt32_object.setParameter( rings[nc].getAveraged_v(), 
                                                               particles[jetID].Mom, particles[iscat].Mom, particles[jscat].Mom,
                                                               F1, F2, F3, sqrt( s ), md2g_wo_as * coupling::get_constant_coupling() / s, lambda_scaled, coupling::get_constant_coupling(), theConfig->isMd2CounterTermInI23(), theConfig->isMatrixElement23_22qt(), theConfig->get23FudgeFactorLpm(), nGluons );  // create scattering32 object for the given 3 particles
                  I32 = scatt32_object.getIntegral32_withPrefactors();                        // get the integral I32 for the given 3 particles

                  
                  probab32 += I32 * dt / ( pow( dv, 2.0 ) * pow( testpartcl, 2.0 ) );
                }
              }
            }
          }
        }
      }
      R32 = probab32 / dt * rings[nc].getGamma();  // fm^-1
    }
    else
    {
      R32 = 0;
    }
    //-------------------------------------------------------------
    
    lambda = 1 / ( R22 + R23 + R32 ); //fm
    lambda = lambda / 0.197;  //GeV^-1

    if ( lambdaArray.size() < 4 )
      lambdaArray.push_back( lambda );
    else
    {
      lambdaArray.push_back( lambda );
      lambdaArray.pop_front();
    }

    if ( lambdaArray.size() == 4 )
    {
      lambdaAvr = 0;
      for ( unsigned int m = 0; m < lambdaArray.size(); m++ )
      {
        lambdaAvr += lambdaArray[m];
      }
      lambdaAvr = lambdaAvr / lambdaArray.size();

      if (( fabs( lambdaArray[3] - lambdaAvr ) / lambdaAvr < epsilon ) &&
          ( fabs( lambdaArray[2] - lambdaAvr ) / lambdaAvr < epsilon ) &&
          ( fabs( lambdaArray[1] - lambdaAvr ) / lambdaAvr < epsilon ) &&
          ( fabs( lambdaArray[0] - lambdaAvr ) / lambdaAvr < epsilon ) )
      {
        converged = true;
      }
      else
      {
        converged = false;
      }
    }
    else
    {
      converged = false;
    }

    ++iter;
    probab22 = probab23 = probab32 = 0;
  }
  while ( !converged && iter < nIterationsMax );

  lambdaAvr = 0;
  for ( unsigned int m = 0; m < lambdaArray.size(); m++ )
    lambdaAvr += lambdaArray[m];
  lambdaAvr = lambdaAvr / lambdaArray.size();

  return lambdaAvr;
}


// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;  replace-tabs on;
