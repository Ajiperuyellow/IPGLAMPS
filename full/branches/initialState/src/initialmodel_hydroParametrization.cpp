//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/initialmodel_hydroParametrization.cpp $
//$LastChangedDate: 2014-11-16 21:18:55 +0100 (日, 16 11月 2014) $
//$LastChangedRevision: 1935 $
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
#include "initialmodel_hydroParametrization.h"
#include "random.h"
#include "particle.h"
#include "FPT_compare.h"

using std::cout;
using std::endl;
using std::vector;
//using namespace ns_casc;


initialModel_HydroParametrization::initialModel_HydroParametrization( const config& _config, string &infoInitialParameters ) :   initialModelWS(_config) 
{

//-------------------------
testparticles = _config.getTestparticles();
//-------------------------
  
//parameters fixed, do not change them
boxLengthX = 20.0;//fm
boxLengthY = 20.0;//fm
boxLengthZ = 8.0;//fm
//-------------------------

//-------------------------
x_0 = -boxLengthX/2.0;
y_0 = -boxLengthY/2.0;
z_0 = -boxLengthZ/2.0;
x_max = boxLengthX/2.0;
y_max = boxLengthY/2.0;
z_max = boxLengthZ/2.0;
//-------------------------

//-------------------------
//values which are fine for our conditions
pT_0 = 0.0;//has to be zero
rap_0 = -50.0;//has to be exact -50!!!
NTF_z0 = -100.0;//has to be at least -100!!!

pT_max = 10.0;//has to be 10
rap_max = - rap_0;//has to be exact 50!!!
NTF_zMax = - NTF_z0;//has to be at least 100!!!
//-------------------------

//-------------------------
theNTFStepSize = 0.5;//fm //1 is a good value
theNTFColumnsX = boxLengthX/theNTFStepSize + 0.5;
theNTFColumnsY = boxLengthY/theNTFStepSize + 0.5;

T_A.resize( theNTFColumnsX, vector<double>( theNTFColumnsY,0.0));    
T_B.resize( theNTFColumnsX, vector<double>( theNTFColumnsY,0.0));     
//-------------------------

//-------------------------
//number of integration
nIntegration = 500000;
nIntegrationNTF = 1000;//1000 is enough
//-------------------------

//-------------------------
//paramters for the heavy-ion collisions
p_A = _config.getA();//Nucleon number of A
p_B = _config.getB();//Nucleon number of B
p_bx = _config.getImpactParameter();
p_by = 0.0;//this is chosen in such a way, that the y-component is always zero
//-------------------------

//-------------------------
//values for andrej
p_K = 0.0135;
p_Q = 1.3;
p_n = 4.0;
p_m = 1.5;
p_sigRap = 1.0;
p_sigZ = 0.13;

p_ksi = 0.54;//fm
p_rho0 = 0.17;



//values for Gabriel and harri
// p_K = 0.6;
// p_Q = 0.6;
// p_n = 1.4;
// p_m = 5.5;
// p_sigRap = 2.25;
// p_sigZ = 0.065;
// 
// p_ksi = 0.54;//fm
// p_rho0 = 0.17;

//-------------------------


//-------------------------
//double sqrtS_perNN = _config.getSqrtS();
//_WoodSaxonParameter.Calculate( p_A, p_bx, sqrtS_perNN);
//-------------------------  

//-------------------------
writeInfoInOutput(infoInitialParameters);
//-------------------------

//-------------------------
//Calculation of the nuclear thickness functions TA and TB
theNTFcalc();
//-------------------------
  
//-------------------------
//Calculation of the total particle number
double totalNumberOfParticlesIntegrated = integration_HIC_TotalNumber() / pow(0.197,3);
  
//Calculation of the total energy 
double totalEnergyIntegrated = integration_HIC_TotalEnergy() / pow(0.197,3); 

cout << endl;
cout << "---------------------------------------------" << endl;
cout << "Number of real particles in HIC: " << totalNumberOfParticlesIntegrated << endl;
cout << "Energy in the system:            " << totalEnergyIntegrated << endl;
cout << "---------------------------------------------" << endl;
cout << endl;
//-------------------------

//--------------------------------------
double iTDAZP_x0 = - 1.0;//should be negative
double iTDAZP_xMax = - iTDAZP_x0;
double iTDAZP_volume = pow( fabs(iTDAZP_x0 - iTDAZP_xMax) , 3 ); 


//Calculation of the total particle density at zero position 
double totalParticleDensityAtZeroPosition = integration_HIC_TotalDensityAtZeroPoint(iTDAZP_x0,iTDAZP_xMax) / iTDAZP_volume / pow(0.197,3);

cout << endl;
cout << "---------------------------------------------" << endl;
cout << "Particle density at volume x0 = y0 = z0 = " << iTDAZP_x0 << " and xM = yM = zM = " << iTDAZP_xMax << ": " << totalParticleDensityAtZeroPosition << " 1/fm^3" << endl;
cout << "---------------------------------------------" << endl;
cout << endl;

//-------------------------
numberOfParticlesToGenerate = totalNumberOfParticlesIntegrated * _config.getTestparticles();
//-------------------------

cout << endl;
cout << "---------------------------------------------" << endl;
cout << "Number of particles with testparticles in HIC: " << numberOfParticlesToGenerate << endl;
cout << "---------------------------------------------" << endl;
cout << endl;
//-------------------------
  
}


void initialModel_HydroParametrization::populateParticleVector( std::vector< Particle >& _particles )
{
  //-----------------------//  
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
  //-----------------------//  

  //-----------------------//   
  sampleParticles( _particles );
  //-----------------------//   
  
}

void initialModel_HydroParametrization::sampleParticles( std::vector< Particle >& _particles )
{
  //-----------------------//  
  hydroParticleType theHydroParticleType;  
  vector<hydroParticleTypeProperties> vecTypeInit(theHydroParticleType.hydroParticleTypeVectorInit);
    
  int nTypes = vecTypeInit.size();
  int nTypesPlusComplete = nTypes + 1; 
  
  vector<int> numberParticles(nTypesPlusComplete,0);  
  vector<int> sizeParticles(nTypesPlusComplete,0);  
  //-----------------------//
  
  //-----------------------//  
  initialisationNumberForEveryParticleSpecies(theHydroParticleType, vecTypeInit, numberParticles); 
  //-----------------------//
  
  //-----------------------//  
  for(int i = 0; i < nTypesPlusComplete; i++)
  {
    if(i == 0){sizeParticles[i] = 0;}
    else{sizeParticles[i] = sizeParticles[i-1] + numberParticles[i-1];}
  } 
  //-----------------------//  
  
  //-----------------------//
  //sampling the momenta, positions and porperties of the particles
  double numAcceptReject = 0.0;
  
  for(int type = 0; type < nTypes; type++)
    { 
      sampleMomenta( _particles, vecTypeInit, type, sizeParticles[type], sizeParticles[type + 1], numAcceptReject );
      samplePositions( _particles, vecTypeInit, type, sizeParticles[type], sizeParticles[type + 1], numAcceptReject );
      sampleProperties( _particles, vecTypeInit, type, sizeParticles[type], sizeParticles[type + 1], numAcceptReject );
    }
  //-----------------------//  
  
  cout << endl;
  cout << "--------------------------------" << endl;
  cout << "Number of tries for reject/accept: " << numAcceptReject << endl;
  cout << "--------------------------------" << endl;
  cout << endl;  
  
}


void initialModel_HydroParametrization::initialisationNumberForEveryParticleSpecies(hydroParticleType& hp , vector<hydroParticleTypeProperties> vecType, vector<int> &numberParticles)
{  
  int nTypes = vecType.size();
  int nTypesPlusComplete = nTypes + 1;   
  
  vector<double> numberParticlesAsDouble(nTypesPlusComplete,0.0);
  
  double volume = boxLengthX * boxLengthY * boxLengthZ;
  double volNorTestp = volume / pow(0.197,3.0) * testparticles;
  double T = 1.0;//the temperature is the same for every particle species
  double fug = 1.0;
  
  //here you can set whether you set the number of particles automatically or according to a manual setting.
  //be careful, because you can have errors here if you do it manually
  bool manualSetting = false;
  
  if(manualSetting == false)
  {
    cout << "------------------------------------------" << endl;
    cout << "Automatic Setting for number of particles!" << endl;
    for (int i = 0; i < nTypes; i++)
    { 
      numberParticlesAsDouble[i] = hp.chooseEqParticleDensity( T,fug,vecType,i ) * volNorTestp;
      numberParticlesAsDouble[nTypes] += numberParticlesAsDouble[i];       
    }
  }
  else
  {
    cout << "------------------------------------------" << endl;
    cout << "Manual Setting for number of particles!" << endl;    
    int upQuark = hp.getParticleType( up, vecType);
    int downQuark = hp.getParticleType( down, vecType);
    
    numberParticlesAsDouble[upQuark] = hp.chooseEqParticleDensity( T,fug,vecType,upQuark ) * volNorTestp;
    numberParticlesAsDouble[downQuark] = hp.chooseEqParticleDensity( T,fug,vecType,downQuark ) * volNorTestp / 1.0;    
    
    for (int i = 0; i < nTypes; i++)
    {     
    numberParticlesAsDouble[nTypes] += numberParticlesAsDouble[i];   
    }
  }

  //RENORMALIZE to the number of particles really generated
  for (int i = 0; i <= nTypes; i++)
  { 
    numberParticlesAsDouble[i] = numberParticlesAsDouble[i] * double(numberOfParticlesToGenerate) / numberParticlesAsDouble[nTypes];      
  }

  //transfer the number of particles from double to int
  for(int i = 0; i < nTypes; i++)
  {
    numberParticles[i] = round( numberParticlesAsDouble[i] );
    numberParticles[nTypes] += numberParticles[i];

    cout << "Number of " << vecType[i].nameOfParticle << ": " << numberParticles[i] << endl;    
  }
}



void initialModel_HydroParametrization::samplePositions( std::vector< Particle >& _particles, vector<hydroParticleTypeProperties> vecType, const int type, const int startValue, const int endValue, double &numAcceptReject)
{
  for(int i = startValue; i < endValue; i++)
    {   
      //-----------------------//
      //sampling the coordinates in space
      sampling_TA_TB(numAcceptReject);
      sampling_z(numAcceptReject);
      //-----------------------//

      _particles[i].Pos = VectorTXYZ(0.0,intPar_x,intPar_y,intPar_z);
    }
}


void initialModel_HydroParametrization::sampleMomenta( std::vector< Particle >& _particles, vector<hydroParticleTypeProperties> vecType, const int type, const int startValue, const int endValue, double &numAcceptReject)
{
  double maximumPT = findMaximumPT();  
  
  for(int i = startValue; i < endValue; i++)
    {     
      //-----------------------//
      //sampling the coordinates in momentum and space
      sampling_rap(numAcceptReject);
      sampling_pT(numAcceptReject,maximumPT);
      //-----------------------//
      
      double phi = 2.0 * M_PI * ran2();    
      _particles[i].Mom = VectorEPxPyPz(intPar_pT*cosh(intPar_rap),
                                        intPar_pT*cos(phi),
                                        intPar_pT*sin(phi),
                                        intPar_pT*sinh(intPar_rap));
    }
}



void initialModel_HydroParametrization::sampleProperties( std::vector< Particle >& _particles, vector<hydroParticleTypeProperties> vecType, const int type, const int startValue, const int endValue, double &numAcceptReject)
{   
    //-----------------------//
    //sampling the properties
    double mass = vecType[type].mass; 
    int flavorStart = vecType[type].flavorTypeStart;
    int flavorEnd = vecType[type].flavorTypeEnd;      

    //set mass
    for(int i = startValue; i < endValue; i++)
      {     
        _particles[i].m = mass;   
      }
    
    //set flavors
    if(flavorStart == flavorEnd)
      {
        for(int i = startValue; i < endValue; i++)
        {         
          _particles[i].FLAVOR = static_cast<FLAVOR_TYPE> (flavorStart);
        }
      }  
    else
      {
      //set flavor for light quarks!!
      //if these are no light quarks, you should loook what the fuck you are doing
      //not written for light quarks
      for(int i = startValue; i < endValue; i+=2)
        {
          setLightQuarkFlavor(_particles, i);
        }    
      }
}


void initialModel_HydroParametrization::setLightQuarkFlavor(std::vector< Particle >& _particles, const int j)
{
  int flav;
  int Nflavor = ParticlePrototype::N_light_flavor;
  
  flav = int( Nflavor * ran2() ) + 1;
  if(flav > Nflavor) {flav = Nflavor;}
  flav = 2*flav;
  _particles[j].FLAVOR = static_cast<FLAVOR_TYPE> (flav);
  _particles[j+1].FLAVOR = static_cast<FLAVOR_TYPE> (flav-1);
}

//-----------------------------//
//-----------------------------//

double initialModel_HydroParametrization::theNiemiParametrization()
{
  double TA,TB;
  
  double part0 = p_K;
  double part1 = theNiemiParametrization_pT();
  double part2 = theNiemiParametrization_rap();
  double part3 = theNiemiParametrization_z();
  theNiemiParametrization_TA_TB(TA,TB);//part4

  double result = part0 * part1 * part2 * part3 * TA * TB; 
  
  return result;
}


double initialModel_HydroParametrization::theNiemiParametrization_pT()
{
  double value = pow( pow(p_Q,p_n) / ( pow(p_Q,p_n) + pow(intPar_pT,p_n) ) ,p_m);
  
  return value;
}

double initialModel_HydroParametrization::theNiemiParametrization_rap()
{
  double value = exp(- pow(intPar_rap,2.0) / (2.0*pow(p_sigRap,2.0)) );
  
  return value;
}

double initialModel_HydroParametrization::theNiemiParametrization_z()
{
  double value = exp(- pow(intPar_z,2.0) / (2.0*pow(p_sigZ,2.0)) );
  
  return value;
}

void initialModel_HydroParametrization::theNiemiParametrization_TA_TB(double &TA, double &TB)
{
  int i = theNTFColumnsX*intPar_x/boxLengthX + 0.5*theNTFColumnsX;
  int j = theNTFColumnsY*intPar_y/boxLengthY + 0.5*theNTFColumnsY;    
  
  TA = T_A[i][j];
  TB = T_B[i][j];
}



void initialModel_HydroParametrization::theNTFcalc()
{
  for(int i = 0; i < theNTFColumnsX; i++){
    for(int j = 0; j < theNTFColumnsY; j++){
        
      intPar_x = double(i)/theNTFColumnsX*boxLengthX - boxLengthX/2.0 + theNTFStepSize/2.0;
      intPar_y = double(j)/theNTFColumnsY*boxLengthY - boxLengthY/2.0 + theNTFStepSize/2.0;       

      T_A[i][j] = integration_NuclearThicknessFunctionA();
      T_B[i][j] = integration_NuclearThicknessFunctionB();
      } 
  }
}

//------------------------------------------//
double initialModel_HydroParametrization::nuclearThicknessFunctionA()
{
  double T_A = nuclearThicknessFunction_woodSaxon(p_A, intPar_x + 0.5*p_bx, intPar_y + 0.5*p_by );

  return T_A;
}

double initialModel_HydroParametrization::nuclearThicknessFunctionB()
{
  double T_B = nuclearThicknessFunction_woodSaxon(p_B, intPar_x - 0.5*p_bx, intPar_y - 0.5*p_by );
  
  return T_B;
}

double initialModel_HydroParametrization::nuclearThicknessFunction_woodSaxon(const double nuc, const double s_x, const double s_y)
{
double p_R0 =1.12*pow(nuc,1.0/3.0) - 0.86*pow(nuc,-1.0/3.0);

double r = sqrt( pow(s_x,2.0) + pow(s_y,2.0) + pow(intPar_NTF_z,2.0) );

double rho_nuc = p_rho0 / ( 1.0 + exp( (r-p_R0) / p_ksi ) );

return rho_nuc;
}
//------------------------------------------//

//----------------------------------
//----------------------------------

double initialModel_HydroParametrization::integration_HIC_TotalDensityAtZeroPoint( const double x_0, const double x_max)
{
  double g_integration_HIC_TotalDensityAtZeroPoint(double *k, size_t dim, void *params);
  //--------------------------------
  //Vegas
  //--------------------------------  
  double res, err;

  double xl[5] = { x_0,x_0,x_0,pT_0,rap_0};
  double xu[5] = { x_max,x_max,x_max,pT_max,rap_max};

  const gsl_rng_type *TT;
  gsl_rng *r;
     
  gsl_monte_function G = { &g_integration_HIC_TotalDensityAtZeroPoint, 5, this};//this

  size_t calls = nIntegration;
     
  gsl_rng_env_setup ();
     
  TT = gsl_rng_default;
  r = gsl_rng_alloc (TT);

  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (5);
    
  gsl_monte_vegas_integrate (&G, xl, xu, 5, calls, r, s, &res, &err);

  gsl_monte_vegas_free (s);

  gsl_rng_free (r);
  //--------------------------------
  //Vegas
  //--------------------------------  
  
  return res;
}

double g_integration_HIC_TotalDensityAtZeroPoint(double *k, size_t dim, void *params)
{  
  initialModel_HydroParametrization* parameter = (initialModel_HydroParametrization*) params;
  
  parameter->intPar_x = k[0];
  parameter->intPar_y = k[1]; 
  parameter->intPar_z = k[2];
  parameter->intPar_pT = k[3];
  parameter->intPar_rap = k[4];  

  double result = parameter->theNiemiParametrization() * 2.0 * M_PI * pow(parameter->intPar_pT,1);
       
  return result;
}



double initialModel_HydroParametrization::integration_HIC_TotalNumber()
{
  double g1(double *k, size_t dim, void *params);
  //--------------------------------
  //Vegas
  //--------------------------------  
  double res, err;

  double xl[5] = { x_0,y_0,z_0,pT_0,rap_0};
  double xu[5] = { x_max,y_max,z_max, pT_max,rap_max};

  const gsl_rng_type *TT;
  gsl_rng *r;
     
  gsl_monte_function G = { &g1, 5, this};//this

  size_t calls = nIntegration;
     
  gsl_rng_env_setup ();
     
  TT = gsl_rng_default;
  r = gsl_rng_alloc (TT);

  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (5);
    
  gsl_monte_vegas_integrate (&G, xl, xu, 5, calls, r, s, &res, &err);

  gsl_monte_vegas_free (s);

  gsl_rng_free (r);
  //--------------------------------
  //Vegas
  //--------------------------------  
  
  return res;
}

double g1(double *k, size_t dim, void *params)
{  
  initialModel_HydroParametrization* parameter = (initialModel_HydroParametrization*) params;
  
  parameter->intPar_x = k[0];
  parameter->intPar_y = k[1]; 
  parameter->intPar_z = k[2];
  parameter->intPar_pT = k[3];
  parameter->intPar_rap = k[4];  

  double result = parameter->theNiemiParametrization() * 2.0 * M_PI * pow(parameter->intPar_pT,1);
       
  return result;
}



double initialModel_HydroParametrization::integration_HIC_TotalEnergy()
{
  double gE(double *k, size_t dim, void *params);
  //--------------------------------
  //Vegas
  //--------------------------------  
  double res, err;

  double xl[5] = { x_0,y_0,z_0,pT_0,rap_0};
  double xu[5] = { x_max,y_max,z_max, pT_max,rap_max};


  const gsl_rng_type *TT;
  gsl_rng *r;
     
  gsl_monte_function G = { &gE, 5, this};//this

  size_t calls = nIntegration;
     
  gsl_rng_env_setup ();
     
  TT = gsl_rng_default;
  r = gsl_rng_alloc (TT);

  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (5);
    
  gsl_monte_vegas_integrate (&G, xl, xu, 5, calls, r, s, &res, &err);

  gsl_monte_vegas_free (s);

  gsl_rng_free (r);
  //--------------------------------
  //Vegas
  //--------------------------------  
  
  return res;
}

double gE(double *k, size_t dim, void *params)
{  
  initialModel_HydroParametrization* parameter = (initialModel_HydroParametrization*) params;
  
  parameter->intPar_x = k[0];
  parameter->intPar_y = k[1]; 
  parameter->intPar_z = k[2];
  parameter->intPar_pT = k[3];
  parameter->intPar_rap = k[4];  
 

  double result = parameter->theNiemiParametrization() * 2.0 * M_PI * pow(parameter->intPar_pT,2) * cosh(parameter->intPar_rap);//the pT^2 has to be there!!!
       
  return result;
}


double initialModel_HydroParametrization::integration_NuclearThicknessFunctionA()
{
  double g2(double *k, size_t dim, void *params);
  //--------------------------------
  //Vegas
  //--------------------------------  
  double res, err;

  double xl[1] = { NTF_z0 };
  double xu[1] = { NTF_zMax };

  const gsl_rng_type *TT;
  gsl_rng *r;
     
  gsl_monte_function G = { &g2, 1, this};//this

  size_t calls = nIntegrationNTF;
     
  gsl_rng_env_setup ();
     
  TT = gsl_rng_default;
  r = gsl_rng_alloc (TT);

  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (1);
    
  gsl_monte_vegas_integrate (&G, xl, xu, 1, calls, r, s, &res, &err);

  gsl_monte_vegas_free (s);

  gsl_rng_free (r);
  //--------------------------------
  //Vegas
  //--------------------------------  
  
  return res;
}

double g2(double *k, size_t dim, void *params)
{ 
  initialModel_HydroParametrization* parameter = (initialModel_HydroParametrization*) params;
  
  parameter->intPar_NTF_z = k[0];  

  double result = parameter->nuclearThicknessFunctionA();

  return result;
} 


double initialModel_HydroParametrization::integration_NuclearThicknessFunctionB()
{
  double g3(double *k, size_t dim, void *params);
  //--------------------------------
  //Vegas
  //--------------------------------  
  double res, err;

  double xl[1] = { NTF_z0 };
  double xu[1] = { NTF_zMax };

  const gsl_rng_type *TT;
  gsl_rng *r;
     
  gsl_monte_function G = { &g3, 1, this};//this

  size_t calls = nIntegrationNTF;
     
  gsl_rng_env_setup ();
     
  TT = gsl_rng_default;
  r = gsl_rng_alloc (TT);

  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (1);
     
  gsl_monte_vegas_integrate (&G, xl, xu, 1, calls, r, s, &res, &err);
  
  gsl_monte_vegas_free (s);

  gsl_rng_free (r);
  //--------------------------------
  //Vegas
  //--------------------------------  
  
  return res;
}

double g3(double *k, size_t dim, void *params)
{  
  initialModel_HydroParametrization* parameter = (initialModel_HydroParametrization*) params;
  
  parameter->intPar_NTF_z = k[0]; 

  double result = parameter->nuclearThicknessFunctionB();
  
  return result;
} 

double initialModel_HydroParametrization::sampling_rap(double &numAcceptReject)
  {
    double ratio,reference,probability;  
      
    intPar_rap = 0.0; 
    reference = theNiemiParametrization_rap(); //das cosh(y) muss hier weg.... /cosh(intPar_rap)

      do{
        intPar_rap = (rap_max - rap_0) * ran2() + rap_0;
        
        probability = theNiemiParametrization_rap();
        
        ratio = probability/reference;
        if(ratio > 1.0){cout << "Ratio for rapidity is " << ratio << endl;}
        numAcceptReject++;
      }while( (ran2() > ratio) ); //accept if random value [0,f(E)] is in [0,p(E)] 

      return 0;
  }


double initialModel_HydroParametrization::sampling_pT(double &numAcceptReject, const double maximum)
  {
    double ratio,reference,probability;  

    reference = maximum;

      do{
        intPar_pT = ran2()*pT_max;
        
        probability = theNiemiParametrization_pT() * intPar_pT;
        
        ratio = probability/reference;
        if(ratio > 1.0){cout << "Ratio for pT is " << ratio << endl;}
        numAcceptReject++;      
      }while( (ran2() > ratio) ); //accept if random value [0,f(E)] is in [0,p(E)] 

  return 0;
  }


double initialModel_HydroParametrization::sampling_z(double &numAcceptReject)
  {
    double ratio,reference,probability;  
      
    intPar_z = 0.0; 
    reference = theNiemiParametrization_z(); 

      do{
        intPar_z = boxLengthZ*ran2() - boxLengthZ/2.0;
        
        probability = theNiemiParametrization_z();
        
        ratio = probability/reference;
        if(ratio > 1.0){cout << "Ratio for z is " << ratio << endl;}
        numAcceptReject++;      
      }while( (ran2() > ratio) ); //accept if random value [0,f(E)] is in [0,p(E)] 

  return 0;
  }
    
    
double initialModel_HydroParametrization::sampling_TA_TB(double &numAcceptReject)
  {
    double ratio,TA,TB,reference,probability;  
      
    intPar_x = 0.0;
    intPar_y = 0.0; 
    theNiemiParametrization_TA_TB(TA,TB); 
    
    reference = TA*TB;

      do{
        intPar_x = boxLengthX*ran2() - boxLengthX/2.0;
        intPar_y = boxLengthY*ran2() - boxLengthY/2.0;
        
        theNiemiParametrization_TA_TB(TA,TB);
        
        probability = TA*TB;
        
        ratio = probability/reference;
        
        if(ratio > 1.01){cout << "Ratio for TA and TB is " << ratio << endl;}
        numAcceptReject++;      
      }while( (ran2() > ratio) ); //accept if random value [0,f(E)] is in [0,p(E)] 
 
      return 0;
 }
  

double initialModel_HydroParametrization::findMaximumPT()
{
  cout << endl;
  cout << "---------------------------" << endl;
  cout << "Try to find the maximum value for the pt function:" << endl;
  
  //find the maximum value
  double stepSize = 0.000001 * pT_max;//size of the step
  double value_start = 0.0;//start value has to be zero in this case
  double value_go;//temporary value
  bool found = false;//has to be false in the beginning
  double result;//the final result
  intPar_pT = 0.0;
  
  do{
      value_go = theNiemiParametrization_pT() * intPar_pT;
      
      if(value_go<value_start)
        {
          result = value_start + 0.01 * value_start;
          found = true;
        }
      else
        {
          value_start = value_go;
          intPar_pT += stepSize;
        }
      
    }while(!found);

    cout << "maximum found at: " << result << endl;   
    cout << "---------------------------" << endl; 
    cout << endl;

    return result;
}      

 
void initialModel_HydroParametrization::writeInfoInOutput( string &infoInitialParameters )
{
 
//-------------------------//
//write in the output
infoInitialParameters += "#  The initial parameters for the hydro parametrization: \n";

string str_pT_0 = "#  pT_0 = ";
string str_pT_max = "#  pT_max = ";    
string str_rap_0 = "#  rap_0 = ";
string str_rap_max = "#  rap_max = ";    
string str_NTF_z0  = "#  NTF_z0 = ";
string str_NTF_zMax  = "#  NTF_zMax = "; 
string str_theNTFStepSize = "#  theNTFStepSize = ";
string str_nIntegration = "#  nIntegration = ";
string str_nIntegrationNTF = "#  nIntegrationNTF = ";    
string str_p_K = "#  p_K = ";
string str_p_Q = "#  p_Q = ";    
string str_p_n  = "#  p_n = ";
string str_p_m  = "#  p_m = "; 
string str_p_sigRap = "#  p_sigRap = ";
string str_p_sigZ = "#  p_sigZ = ";
string str_p_A = "#  p_A = ";    
string str_p_B  = "#  p_B = ";
string str_p_bx  = "#  p_bx = "; 
string str_p_by = "#  p_by = ";
string str_p_ksi  = "#  p_ksi = "; 
string str_p_rho0 = "#  p_rho0 = ";
    
std::stringstream valueStr_pT_0;
std::stringstream valueStr_pT_max;    
std::stringstream valueStr_rap_0;
std::stringstream valueStr_rap_max;    
std::stringstream valueStr_NTF_z0;
std::stringstream valueStr_NTF_zMax; 
std::stringstream valueStr_theNTFStepSize;  
std::stringstream valueStr_nIntegration;
std::stringstream valueStr_nIntegrationNTF;    
std::stringstream valueStr_p_K;
std::stringstream valueStr_p_Q;    
std::stringstream valueStr_p_n;
std::stringstream valueStr_p_m; 
std::stringstream valueStr_p_sigRap;    
std::stringstream valueStr_p_sigZ;
std::stringstream valueStr_p_A;    
std::stringstream valueStr_p_B;
std::stringstream valueStr_p_bx;    
std::stringstream valueStr_p_by;
std::stringstream valueStr_p_ksi; 
std::stringstream valueStr_p_rho0;
    
valueStr_pT_0 << pT_0;
valueStr_pT_max << pT_max;    
valueStr_rap_0 << rap_0;
valueStr_rap_max << rap_max;     
valueStr_NTF_z0 << NTF_z0;
valueStr_NTF_zMax << NTF_zMax;
valueStr_theNTFStepSize << theNTFStepSize;
valueStr_nIntegration << nIntegration;
valueStr_nIntegrationNTF << nIntegrationNTF;    
valueStr_p_K << p_K;
valueStr_p_Q << p_Q;     
valueStr_p_n << p_n;
valueStr_p_m << p_m;
valueStr_p_sigRap << p_sigRap;
valueStr_p_sigZ << p_sigZ;
valueStr_p_A << p_A;    
valueStr_p_B << p_B;
valueStr_p_bx << p_bx;     
valueStr_p_by << p_by;
valueStr_p_ksi << p_ksi;
valueStr_p_rho0 << p_rho0;     


//-----------------------------------//
//save the particle distribution in the string
infoInitialParameters += str_pT_0 + valueStr_pT_0.str() + "\n"
+ str_pT_max + valueStr_pT_max.str() + "\n"  
+ str_rap_0 + valueStr_rap_0.str() + "\n"
+ str_rap_max + valueStr_rap_max.str() + "\n"    
+ str_NTF_z0 + valueStr_NTF_z0.str() + "\n"
+ str_NTF_zMax + valueStr_NTF_zMax.str() + "\n"
+ str_theNTFStepSize + valueStr_theNTFStepSize.str() + "\n"
+ str_nIntegrationNTF + valueStr_nIntegration.str() + "\n"
+ str_nIntegrationNTF + valueStr_nIntegrationNTF.str() + "\n"  
+ str_p_K + valueStr_p_K.str() + "\n"
+ str_p_Q + valueStr_p_Q.str() + "\n"    
+ str_p_n + valueStr_p_n.str() + "\n"
+ str_p_m + valueStr_p_m.str() + "\n"
+ str_p_sigRap + valueStr_p_sigRap.str() + "\n"
+ str_p_sigZ + valueStr_p_sigZ.str() + "\n"
+ str_p_A + valueStr_p_A.str() + "\n"  
+ str_p_B + valueStr_p_B.str() + "\n"
+ str_p_bx + valueStr_p_bx.str() + "\n"    
+ str_p_by + valueStr_p_by.str() + "\n"
+ str_p_ksi + valueStr_p_ksi.str() + "\n"
+ str_p_rho0 + valueStr_p_rho0.str() + "\n";
//-----------------------------------//
 
}
