//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/hydroParticleType.cpp $
//$LastChangedDate: 2017-12-14 20:59:36 +0100 (Do, 14. Dez 2017) $
//$LastChangedRevision: 2645 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <string>

#include "hydroParticleType.h"
#include "particleprototype.h"
#include "tools.h"

using namespace std;


/** 
 * These Analysis Routines are for obtaining the T_munu and
 * N_mu. They are not well documentated. Ask
 * gallmeister@th.physik.uni-frankfurt.de.  
 **/ 

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

const hydroParticleTypeProperties hydroParticleType::hydro_gluon( 16.0, 0.0, true, gluon, gluon, "gluon" );
const hydroParticleTypeProperties hydroParticleType::hydro_up( 6.0, 0.0, true, up, up, "up" );
const hydroParticleTypeProperties hydroParticleType::hydro_antiup( 6.0, 0.0, true, anti_up, anti_up, "antiup" );
const hydroParticleTypeProperties hydroParticleType::hydro_down( 6.0, 0.0, true, down, down, "down" );
const hydroParticleTypeProperties hydroParticleType::hydro_antidown( 6.0, 0.0, true, anti_down, anti_down, "antidown" );
const hydroParticleTypeProperties hydroParticleType::hydro_strange( 6.0, 0.101, false, strange, strange, "strange" );
const hydroParticleTypeProperties hydroParticleType::hydro_antistrange( 6.0, 0.101, false, anti_strange, anti_strange, "antistrange" );
const hydroParticleTypeProperties hydroParticleType::hydro_charm( 6.0, 1.27, false, charm, charm, "charm" );
const hydroParticleTypeProperties hydroParticleType::hydro_anticharm( 6.0, 1.27, false, anti_charm, anti_charm, "anticharm" );
const hydroParticleTypeProperties hydroParticleType::hydro_lightQuarks( ParticlePrototype::N_light_flavor * 2.0 * 2.0 * 3.0, 0.0, true, up, anti_strange, "lightQuarks" );
const hydroParticleTypeProperties hydroParticleType::hydro_charmQuarks( 1.0 * 2.0 * 2.0 * 3.0, ParticlePrototype::Mcharm, false, charm, anti_charm, "charmQuark" );
const hydroParticleTypeProperties hydroParticleType::hydro_massA( 16.0, ParticlePrototype::MmassA, false, massA, massA, "massA" );
const hydroParticleTypeProperties hydroParticleType::hydro_massB( 16.0, ParticlePrototype::MmassB, false, massB, massB, "massB" );
const hydroParticleTypeProperties hydroParticleType::hydro_masslessC( 16.0, ParticlePrototype::MmasslessC, true, masslessC, masslessC, "masslessC" );
const hydroParticleTypeProperties hydroParticleType::hydro_masslessD( 16.0, ParticlePrototype::MmasslessD, true, masslessD, masslessD, "masslessD" );
  
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


std::vector<hydroParticleTypeProperties> hydroParticleType::hydroParticleTypeVectorCommon = create_vector<hydroParticleTypeProperties>
  ( hydro_gluon );
  
//   ( hydro_up )( hydro_antiup )( hydro_down )( hydro_antidown)( hydro_strange )( hydro_antistrange );

std::vector<hydroParticleTypeProperties> hydroParticleType::hydroParticleTypeVectorInit = create_vector<hydroParticleTypeProperties>
  ( hydro_gluon );
//   ( hydro_up )( hydro_antiup )( hydro_down )( hydro_antidown)( hydro_strange )( hydro_antistrange );

std::vector<hydroParticleTypeProperties> hydroParticleType::hydroParticleTypeVectorAnalysis = create_vector<hydroParticleTypeProperties>
  ( hydro_gluon );
//   ( hydro_up )( hydro_antiup )( hydro_down )( hydro_antidown)( hydro_strange )( hydro_antistrange );


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------



int hydroParticleType::getParticleType(const int particleFlavor, vector<hydroParticleTypeProperties> vecType)
{
  //---------//
  int type = -1; // in case it is not in list, give -1 back
  int nn = -1;
  //---------//

  vector<hydroParticleTypeProperties>::iterator it;
  
  for ( it = vecType.begin(); it != vecType.end(); ++it )
  {
    nn++;
    
    if ( particleFlavor >= ( *it ).flavorTypeStart && particleFlavor <= ( *it ).flavorTypeEnd )
    {
      type = nn;
      break;
    }
  }

  return type;  
}


double hydroParticleType::chooseEqParticleDensity( const double T, const double fug, vector<hydroParticleTypeProperties> vecType, const int type)
{

  double nDensity;
  double mass = vecType[type].mass; 
  double deg = vecType[type].degeneracyFactor;
  bool masslessON = vecType[type].masslessON;
  
  if(masslessON)
  {
    nDensity = fug * deg / pow( M_PI,2 ) * pow(T,3.0);
  }
  else
  {
    nDensity = fug * deg * 0.5 / pow( M_PI,2 ) * pow(mass,2.0) * T * gsl_sf_bessel_Kn(2,mass/T);  
  }
  
  return nDensity;
}
  
double hydroParticleType::chooseEqPressure( const double T, const double fug, vector<hydroParticleTypeProperties> vecType, const int type)
{
  double Pressure = chooseEqParticleDensity(T,fug,vecType,type) * T;

  return Pressure;
}

double hydroParticleType::chooseEqEnergyDensity( const double T, const double fug, vector<hydroParticleTypeProperties> vecType, const int type)
{

  double eDensity;
  double mass = vecType[type].mass; 
  double deg = vecType[type].degeneracyFactor;
  bool masslessON = vecType[type].masslessON;
  
  if(masslessON)
  {
    eDensity = 3.0 * chooseEqPressure(T,fug,vecType,type);
  }
  else
  {
    eDensity = 3.0 * chooseEqPressure(T,fug,vecType,type) + fug * deg * 0.5 / pow( M_PI,2 ) * pow(mass,3.0) * T * gsl_sf_bessel_Kn(1,mass/T);   
  }
  
  return eDensity;
}
  
double hydroParticleType::chooseEqEntropyDensity( const double T, const double fug, vector<hydroParticleTypeProperties> vecType, const int type)
{
  double entropyDensity = (chooseEqEnergyDensity(T,fug,vecType,type) + chooseEqPressure(T,fug,vecType,type) - T * log(fug) * chooseEqParticleDensity(T,fug,vecType,type)) / T;

  return entropyDensity;
} 
  


double hydroParticleType::chooseEqTemperature(const double eDensity, const double nDensity, vector<hydroParticleTypeProperties> vecType, const int type, bool & solution)
{
  double temperature;
  double mass = vecType[type].mass; 
  bool masslessON = vecType[type].masslessON;
  const double massOverTBorder = pow(10.0,10);

  solution = true;

  //eDensity unit GeV/fm³
  //nDensity unit 1/fm³
  //T unit GeV
  
  
  if(masslessON)
  {
    temperature = eDensity / (3.0 * nDensity);
  }
  else
  {
    //-----------------------------//
    //Sekantenverfahren
    //-----------------------------//
    static const double precision = 0.001;
    double T,T_min,T_max;
    double F_min,F = 0.0;
      
    T_min = 1.0/3.0*(eDensity / nDensity - mass);//minimum value for the temperature
    if(FPT_COMP_LE(T_min,0.0))
    {
      T_min = 0.001;
      cout << "Error in T_min for getting temperature" << endl;
    }  //in case T_min is negative
    T_max = eDensity / (3.0 * nDensity);//the maximum value is the energy density for an ultrarelativistic gas
    T = T_max;//its start value
    
    int counter = 0;
    do
    {
      if(mass/T > massOverTBorder || mass/T < 1.0/massOverTBorder)
      {
        temperature = 0.0;
        solution = false;      
        break;
      }
      else
      {
        F_min = 3.0 * T_min + mass * gsl_sf_bessel_Kn(1,mass/T_min) / gsl_sf_bessel_Kn(2,mass/T_min) - eDensity / nDensity;
        F = 3.0 * T + mass * gsl_sf_bessel_Kn(1,mass/T) / gsl_sf_bessel_Kn(2,mass/T) - eDensity / nDensity;
      }
          
      temperature = T - (T - T_min) / (F - F_min) * F;

      T = temperature;//set T to the new value
          
      //-------------------
      counter++;
      if(counter > 1000){break;}
      //-------------------  
          
    } while(fabs(F) > precision);
        
    //      if(fabs(F)> precision)
    //        {cout << "Error in iterating the temperature: Value not reached!" << endl;}
    //      
    if(temperature > T_max || temperature < T_min)
    {
      //cout << "Error in iterating the temperature: Not in range of T_min or T_max!" << endl;
      temperature = 0.0;//set it to the value of a massless gas
      solution = false;
    }
  }

  return temperature;      
}



double hydroParticleType::landau_energyDensity(const double T00, const double T11, const double T22, const double T33, const double T10, const double T20, const double T30, const double T21, const double T31, const double T32, bool & solution)
{
  
  //------------------------------------------------------------------------------------------------
  //energy density in landau frame
  //this is valid for a common equation of 4. degree
  //still in beta-phase!!
  //-----------------------------------------
    
  double energyDensity = 0.0;
  double c_A,c_B,c_C,c_D,c_E;
  double d_B,d_C,d_D,d_E;
  double e_C,e_D,e_E;
  double e_R,e_S,e_T;
  bool e_DposSign;
    
  complex<double> y1(0.0, 0.0 );
  complex<double> y2(0.0, 0.0 );
  complex<double> y3(0.0, 0.0 );
  complex<double> y4(0.0, 0.0 );
  complex<double> e1(0.0, 0.0 );
  complex<double> e2(0.0, 0.0 );
  complex<double> e3(0.0, 0.0 );
  complex<double> e4(0.0, 0.0 ); 
  complex<double> z1(0.0, 0.0 );
  complex<double> z2(0.0, 0.0 );
  complex<double> z3(0.0, 0.0 );
  complex<double> sqrt_z1(0.0, 0.0 );
  complex<double> sqrt_z2(0.0, 0.0 );
  complex<double> sqrt_z3(0.0, 0.0 );
    
  //-----------------------------------------
  //Calculating the energydensity with analytical method
  //-----------------------------------------
  //it has to be of 4th degree, otherwise the solution is not correct. Be careful!!

  c_A = -1.0;
  c_B = T00 - T11 - T22 - T33;
    
  c_C = (T00*(T11+T22+T33) - T11*T22 - T11*T33 - T22*T33 + pow(T32,2) + pow(T21,2) + pow(T31,2) - pow(T10,2) - pow(T20,2) - pow(T30,2));

  c_D = (T00*(T11*T22 + T11*T33 + T22*T33 - pow(T32,2) - pow(T21,2) - pow(T31,2))
	 - T11*T22*T33 - 2.0*T21*T32*T31 + T11*pow(T32,2) + T33*pow(T21,2) + T22*pow(T31,2)
	 - pow(T10,2)*(T22 + T33) + 2.0*T21*T10*T20 + 2.0*T31*T10*T30
	 - pow(T20,2)*(T11 + T33) + 2.0*T20*T32*T30
	 - pow(T30,2)*(T11 + T22));

  c_E = (T00*T11*T22*T33 + 2.0*T00*T21*T32*T31 + pow(T10,2)*pow(T32,2) + pow(T31,2)*pow(T20,2) + pow(T21,2)*pow(T30,2)
	 + 2.0*T21*T10*T20*T33 + 2.0*T20*T11*T32*T30 + 2*T31*T10*T30*T22
	 - 2.0*T21*T10*T32*T30 - 2.0*T31*T10*T20*T32 - 2.0*T31*T21*T30*T20
	 - T00*T11*pow(T32,2) - T00*T33*pow(T21,2) - T00*T22*pow(T31,2) - T22*T33*pow(T10,2) - T11*T33*pow(T20,2) - T11*T22*pow(T30,2));

  //switch to general form
  //x^4 + d_B*y^3 + d_C*y^2 + d_D*y + d_E = 0
  d_B = c_B/c_A;
  d_C = c_C/c_A;
  d_D = c_D/c_A;
  d_E = c_E/c_A;
    
  //substitution: x = y -d_A; still 4 degree, but without y^3
  //y^4 + e_C*y^2 + e_D*y + e_E = 0
  e_C = -3.0/8.0*pow( d_B,2 ) + d_C;
  e_D = 1.0/8.0*pow( d_B,3 ) - d_B*d_C/2.0 + d_D;
  e_E = -3.0*pow( d_B,4 )/256.0 + pow( d_B,2 )*d_C/16.0 - d_B*d_D/4.0 + d_E;
    
  //transfrom to a cubic equation;
  //z^3 + e_R*z^2 + e_S*z + e_T = 0
  e_R = -2.0 * e_C;
  e_S = pow( e_C,2 ) - 4.0*e_E;
  e_T = pow( e_D,2 );

  //-----------------------------------------
  //the solutions of the cubic equation
  //-----------------------------------------
  do
  {
    double p,q,R;

    //------------------------------------------------------------------
    p = e_S - pow( e_R,2 ) / 3.0;
    q = 2.0/27.0*pow( e_R,3 ) - e_R*e_S/3.0 + e_T;
    R = pow(q/2.0,2) + pow(p/3.0,3);
      
    if(R < 0.0)
    {
      double a,b,r;
      double phi;
        
      a = q/2.0;
      b = - sqrt(-R);
      r = fabs(sqrt( pow(a,2) - R));

      phi = atan(b/a);

      if(phi > 0.0)
      {
	z1 = + 2.0*pow(r,1.0/3.0)*cos(phi/3.0);
	z2 = - 2.0*pow(r,1.0/3.0)*cos(phi/3.0 - M_PI/3.0);
	z3 = - 2.0*pow(r,1.0/3.0)*cos(phi/3.0 + M_PI/3.0);   
      }
      else
      {
	phi = phi + M_PI;

	z1 = + 2.0*pow(r,1.0/3.0)*cos(phi/3.0);
	z2 = - 2.0*pow(r,1.0/3.0)*cos(phi/3.0 - M_PI/3.0);
	z3 = - 2.0*pow(r,1.0/3.0)*cos(phi/3.0 + M_PI/3.0);              
      }
    }
    else
    {
      double v,u,u$3,v$3;
                        
      u$3 = -q/2.0 + sqrt(R);
      v$3 = -q/2.0 - sqrt(R);
                        
      if(v$3 < 0.0)
      {
	v$3 = -v$3;
	v = -pow(v$3,1.0/3.0);
      }
      else
      {
	v = pow(v$3,1.0/3.0);
      }

      if(u$3 < 0.0)
      {
	u$3 = -u$3;
	u = -pow(u$3,1.0/3.0);
      }
      else
      {
	u = pow(u$3,1.0/3.0);
      }
            
      complex<double> z2_temp(-(u+v)/2.0,-(u-v)/2.0*sqrt(3.0));
      complex<double> z3_temp(-(u+v)/2.0, (u-v)/2.0*sqrt(3.0));   
          
      z1 = u + v;
      z2 = z2_temp;
      z3 = z3_temp;
    }
  } while(false);

  //Substituin back with z = z - e_R/3
  z1 = z1 - e_R/3.0;
  z2 = z2 - e_R/3.0;
  z3 = z3 - e_R/3.0;

  //define whether -e_D is positiv or negativ; need for next step...
  if( -e_D >= 0.0){e_DposSign = true;}
  else{e_DposSign = false;}
    
  //define the square roots, they have to be negativ except
  //for sqrt(z2), whiich sign depends on e_D
  sqrt_z1 = - sqrt( -z1 );
  if(e_DposSign){sqrt_z2 = sqrt( -z2);}
  else{sqrt_z2 = - sqrt( -z2);}
  sqrt_z3 = - sqrt( -z3 );  

  //find solution for y by adding all z
  y1 = ( sqrt_z1 + sqrt_z2 + sqrt_z3 )/2.0; 
  y2 = ( sqrt_z1 - sqrt_z2 - sqrt_z3 )/2.0;
  y3 = ( -sqrt_z1 + sqrt_z2 - sqrt_z3 )/2.0; 
  y4 = ( -sqrt_z1 - sqrt_z2 + sqrt_z3 )/2.0;    
    
  //substituin back x = y-d_B/4
  e1 = y1 - d_B/4.0;
  e2 = y2 - d_B/4.0;
  e3 = y3 - d_B/4.0;
  e4 = y4 - d_B/4.0;
    
  //in general there is one positive solution, which represents the energy density
  //and three negative one, which are the pressure of the energy momentum tensor;
  //having only 2 particles, there is a possibility for two positive solution, so
  //the largest one has to be chosen;
  //in the case there is only one particle, the solution is definitve not exact!!
  //then one should set it to T00, as approximation!!
    
  int nn_e = 0;
  if(e1.real()> 0.0){energyDensity = e1.real(); nn_e++;}
  if(e2.real()> 0.0){energyDensity = e2.real(); nn_e++;}
  if(e3.real()> 0.0){energyDensity = e3.real(); nn_e++;}
  if(e4.real()> 0.0){energyDensity = e4.real(); nn_e++;}
    
  if(nn_e == 0)
  {
    //cout << "No positive solution in eDensity Landau frame! " << endl;
    energyDensity = T00;
    solution = false;
  }
    
  if(nn_e == 1){solution = true;}    
    
  if(nn_e == 2)
  {
    //cout << "2 positive solution in eDensity Landau frame! " << endl;      
    if( e1.real() > energyDensity  ) { energyDensity = e1.real(); }
    if( e2.real() > energyDensity  ) { energyDensity = e2.real(); }
    if( e3.real() > energyDensity  ) { energyDensity = e3.real(); }
    if( e4.real() > energyDensity  ) { energyDensity = e4.real(); }
    solution = true;      
  }
    
  if(nn_e >= 3)
  {
    //cout << "3 or more positive solution in eDensity Landau frame! " << endl;
    energyDensity = T00;
    solution = false;       
  }
    
  /*
    if(!(nn_e == 1))
    {
    //------------------------------------------------------------------
    //Energiedichte im 1-dimensionalen Fall, dient als Abschätzung
    double root_F;
    double p_OneD,q_OneD,energyDensityOneD;

    p_OneD = T33 - T00;
    q_OneD = - T33*T00 + pow(T30,2);
    energyDensityOneD = -(p_OneD/2.0) + sqrt( pow((p_OneD/2.0),2) - q_OneD); // GeV / fm^3
    //------------------------------------------------------------------

    //------------------------------------------------------------------
    //Newton-Verfahren //übergangsweise
    double initial = energyDensityOneD + 1.0;
    //------------------------
              
    cout << initial << endl;
              
    for(int j = 1; j<100000; j++)
    {
    double x_nv,F,dF;
    x_nv = initial;

    F = pow(x_nv,4) + e_C*pow(x_nv,2) + e_D*x_nv + e_E;
    dF = 4.0*pow(x_nv,3) + 2.0*e_C*x_nv + e_D;

    root_F = initial - F/dF;
    initial = root_F;
    }

    energyDensity = root_F;

    cout << "Error in energy density calculations" << endl;
    cout << "Starting Newton-Verfahren for energy density:" << endl;
    cout << "Solutions for e: " << e1.real() << " " << e2.real() << " " << e3.real() << " " << e4.real() << endl;
    cout << "root_F Newton:   " <<  root_F << endl;
    cout << "eDensity old:    " <<  energyDensity << endl;
    cout << "difference:      " <<  energyDensity - root_F << endl;
    }
  */

  return energyDensity;
}


void hydroParticleType::eckart_velocity(const double N0, const double N1, const double N2, const double N3, double& vx, double& vy, double& vz)
{
  vx = N1/N0;
  vy = N2/N0;
  vz = N3/N0;
}

double hydroParticleType::lanEck_velocityAbsolute(const double vx, const double vy, const double vz)
{
  double vv;

  vv = sqrt( pow( vx,2 ) + pow( vy,2 ) + pow( vz,2 ) );

  return vv;
}

void hydroParticleType::landau_velocity(const double , const double T11,const double T22, const double T33, const double T10, const double T20, const double T30, const double T21, const double T31, const double T32, const double energyDensity, double& vx, double& vy, double& vz)
{
  double denominator;
  denominator = (T11+energyDensity)*( (T22+energyDensity)*(T33+energyDensity) - pow(T32,2)) - T21*(T21*(T33+energyDensity) - T32*T31) + T31*(T21*T32 - T31*(T22+energyDensity));
  vx = (T10*((T22+energyDensity)*(T33+energyDensity) - pow(T32,2))- T21*(T20*(T33+energyDensity) - T32*T30) + T31*(T20*T32 - T30*(T22+energyDensity))) / denominator;
  vy = ((T11+energyDensity)*(T20*(T33+energyDensity) - T32*T30) - T10*(T21*(T33+energyDensity) - T32*T31) + T31*(T21*T30 - T20*T31)) / denominator;
  vz = ((T11+energyDensity)*(T30*(T22+energyDensity) - T20*T32) - T21*(T21*T30 - T20*T31) + T10*(T21*T32 - T31*(T22+energyDensity))) / denominator;
}

double hydroParticleType::lanEck_gamma(const double vx, const double vy, const double vz)
{
  double vv, gamma;

  vv = sqrt( pow( vx,2 ) + pow( vy,2 ) + pow( vz,2 ) );
  gamma = 1.0 / sqrt( 1.0 - pow( vv, 2 ) );

  return gamma;
}


double hydroParticleType::lanEck_fugacity(const double nDensity, const double temperature,
					  vector<hydroParticleTypeProperties> vecType, const int type)
{
  double fugacity, eqParticleDensity;

  //nDensity is in/fm^3 
  eqParticleDensity = chooseEqParticleDensity(temperature,1.0,vecType,type) / pow(0.197,3);//1/fm^3 
  fugacity = nDensity / eqParticleDensity;

  return fugacity;
}



double hydroParticleType::lanEck_bulkPressure(const double isoPressure, const double eqPressure)
{
  double bulkPressure = isoPressure - eqPressure;
   
  return bulkPressure;
}

double hydroParticleType::eckart_energyDensity(const double T00, const double T11, const double T22, const double T33, const double T10, const double T20,const double T30, const double T21, const double T31, const double T32, const double vx, const double vy, const double vz)
{
  double eDensity, gamma;
  double u_0,u_1,u_2,u_3;  
  double T01,T02,T03,T12,T13,T23;
 
  T01 = T10;
  T02 = T20;
  T03 = T30;
  T12 = T21;
  T13 = T31;
  T23 = T32;

  gamma = lanEck_gamma( vx, vy, vz);
  u_0 = gamma;
  u_1 = - gamma * vx;
  u_2 = - gamma * vy;
  u_3 = - gamma * vz;
  // double u$0 = gamma;
  // double u$1 = gamma * vx;
  // double u$2 = gamma * vy;
  // double u$3 = gamma * vz;


  eDensity = u_0 * T00 * u_0 + u_0 * T01 * u_1 + u_0 * T02 * u_2 + u_0 * T03 * u_3
    + u_1 * T10 * u_0 + u_1 * T11 * u_1 + u_1 * T12 * u_2 + u_1 * T13 * u_3
    + u_2 * T20 * u_0 + u_2 * T21 * u_1 + u_2 * T22 * u_2 + u_2 * T23 * u_3
    + u_3 * T30 * u_0 + u_3 * T31 * u_1 + u_3 * T32 * u_2 + u_3 * T33 * u_3;
           
  return eDensity;
}



double hydroParticleType::lanEck_particleDensity( const double N0, const double N1, const double N2, const double N3, const double vx, const double vy,const double vz)
{
  double particleDensity, gamma;
  gamma = lanEck_gamma( vx, vy, vz);
  particleDensity = gamma * (N0 - vx*N1 - vy*N2 - vz*N3);
  
  return particleDensity;
}


void hydroParticleType::lanEck_shearStress( const double T00, const double T11, const double T22, const double T33, const double T10, const double T20, const double T30, const double T21, const double T31, const double T32, const double v_x, const double v_y, const double v_z, double & PI00, double & PI11, double & PI22, double & PI33, double & PI10, double & PI20, double & PI30, double & PI21, double & PI31, double & PI32)
{
  //------------------------------------------------------------------------------------------------
  // shear stress tensor und co
  //-----------------------------------------
  double D_0_0,D$0_0,D_0$0,D$0$0;
  double D_1_1,D$1$1,D_2_2,D$2$2,D_3_3,D$3$3;
  double D$1_1,D$2_2,D$3_3;
  double D$1_0,D$1$0,D$0$1,D$2_0,D$2$0,D$0$2,D$3_0,D$3$0,D$0$3;
  double D$0_1,D_1_0,D_0_1,D$0_2,D_2_0,D_0_2,D$0_3,D_3_0,D_0_3;
  double D_2_1,D$2$1,D_1_2,D$1$2,D_3_1,D$3$1,D_1_3,D$1$3,D_3_2,D$3$2,D_2_3,D$2$3;
  double D$2_1,D$1_2,D$3_1,D$1_3,D$3_2,D$2_3;

  //-----------
  lanEck_deltamunu( v_x, v_y, v_z, D_0_0, D$0_0, D_0$0, D$0$0, D_1_1, D$1$1, D_2_2, D$2$2, D_3_3, D$3$3, D$1_1, D$2_2, D$3_3, D$1_0, D$1$0, D$0$1, D$2_0, D$2$0, D$0$2, D$3_0, D$3$0, D$0$3, D$0_1, D_1_0, D_0_1, D$0_2, D_2_0, D_0_2, D$0_3, D_3_0, D_0_3, D_2_1, D$2$1, D_1_2, D$1$2, D_3_1, D$3$1, D_1_3, D$1$3, D_3_2, D$3$2, D_2_3, D$2$3, D$2_1, D$1_2, D$3_1, D$1_3,  D$3_2, D$2_3 );      
  //------------     

  //-----------------------------------------
  // p_mn - Tensor
  //-----------------------------------------
  PI00 = (D$0_0 *D$0_0 - 1.0/3.0*D$0$0*D_0_0)*T00
    + (D$0_1 *D$0_1 - 1.0/3.0*D$0$0*D_1_1)*T11 + (D$0_2 *D$0_2 - 1.0/3.0*D$0$0*D_2_2)*T22 + (D$0_3 *D$0_3 - 1.0/3.0*D$0$0*D_3_3)*T33
    + 2.0*((D$0_1 *D$0_0 - 1.0/3.0*D$0$0*D_1_0)*T10 + (D$0_2 *D$0_0 - 1.0/3.0*D$0$0*D_2_0)*T20 + (D$0_3 *D$0_0 - 1.0/3.0*D$0$0*D_3_0)*T30)
    + 2.0*((D$0_2 *D$0_1 - 1.0/3.0*D$0$0*D_2_1)*T21 + (D$0_3 *D$0_1 - 1.0/3.0*D$0$0*D_3_1)*T31 + (D$0_3 *D$0_2 - 1.0/3.0*D$0$0*D_3_2)*T32);

  PI11 = (D$1_0 *D$1_0 - 1.0/3.0*D$1$1*D_0_0)*T00
    + (D$1_1 *D$1_1 - 1.0/3.0*D$1$1*D_1_1)*T11 + (D$1_2 *D$1_2 - 1.0/3.0*D$1$1*D_2_2)*T22 + (D$1_3 *D$1_3 - 1.0/3.0*D$1$1*D_3_3)*T33
    + 2.0*((D$1_1 *D$1_0 - 1.0/3.0*D$1$1*D_1_0)*T10 + (D$1_2 *D$1_0 - 1.0/3.0*D$1$1*D_2_0)*T20 + (D$1_3 *D$1_0 - 1.0/3.0*D$1$1*D_3_0)*T30)
    + 2.0*((D$1_2 *D$1_1 - 1.0/3.0*D$1$1*D_2_1)*T21 + (D$1_3 *D$1_1 - 1.0/3.0*D$1$1*D_3_1)*T31 + (D$1_3 *D$1_2 - 1.0/3.0*D$1$1*D_3_2)*T32);

  PI22 = (D$2_0 *D$2_0 - 1.0/3.0*D$2$2*D_0_0)*T00
    + (D$2_1 *D$2_1 - 1.0/3.0*D$2$2*D_1_1)*T11 + (D$2_2 *D$2_2 - 1.0/3.0*D$2$2*D_2_2)*T22 + (D$2_3 *D$2_3 - 1.0/3.0*D$2$2*D_3_3)*T33
    + 2.0*((D$2_1 *D$2_0 - 1.0/3.0*D$2$2*D_1_0)*T10 + (D$2_2 *D$2_0 - 1.0/3.0*D$2$2*D_2_0)*T20 + (D$2_3 *D$2_0 - 1.0/3.0*D$2$2*D_3_0)*T30)
    + 2.0*((D$2_2 *D$2_1 - 1.0/3.0*D$2$2*D_2_1)*T21 + (D$2_3 *D$2_1 - 1.0/3.0*D$2$2*D_3_1)*T31 + (D$2_3 *D$2_2 - 1.0/3.0*D$2$2*D_3_2)*T32);

  PI33 = (D$3_0 *D$3_0 - 1.0/3.0*D$3$3*D_0_0)*T00
    + (D$3_1 *D$3_1 - 1.0/3.0*D$3$3*D_1_1)*T11 + (D$3_2 *D$3_2 - 1.0/3.0*D$3$3*D_2_2)*T22 + (D$3_3 *D$3_3 - 1.0/3.0*D$3$3*D_3_3)*T33
    + 2.0*((D$3_1 *D$3_0 - 1.0/3.0*D$3$3*D_1_0)*T10 + (D$3_2 *D$3_0 - 1.0/3.0*D$3$3*D_2_0)*T20 + (D$3_3 *D$3_0 - 1.0/3.0*D$3$3*D_3_0)*T30)
    + 2.0*((D$3_2 *D$3_1 - 1.0/3.0*D$3$3*D_2_1)*T21 + (D$3_3 *D$3_1 - 1.0/3.0*D$3$3*D_3_1)*T31 + (D$3_3 *D$3_2 - 1.0/3.0*D$3$3*D_3_2)*T32);

  PI10 = (D$1_0 *D$0_0 - 1.0/3.0*D$1$0*D_0_0)*T00
    + (D$1_1 *D$0_1 - 1.0/3.0*D$1$0*D_1_1)*T11 + (D$1_2 *D$0_2 - 1.0/3.0*D$1$0*D_2_2)*T22 + (D$1_3 *D$0_3 - 1.0/3.0*D$1$0*D_3_3)*T33
    + (D$1_0 *D$0_1 - 1.0/3.0*D$1$0*D_0_1)*T10 + (D$1_0 *D$0_2 - 1.0/3.0*D$1$0*D_0_2)*T20 + (D$1_0 *D$0_3 - 1.0/3.0*D$1$0*D_0_3)*T30
    + (D$1_1 *D$0_0 - 1.0/3.0*D$1$0*D_1_0)*T10 + (D$1_2 *D$0_0 - 1.0/3.0*D$1$0*D_2_0)*T20 + (D$1_3 *D$0_0 - 1.0/3.0*D$1$0*D_3_0)*T30
    + (D$1_2 *D$0_1 - 1.0/3.0*D$1$0*D_2_1)*T21 + (D$1_3 *D$0_1 - 1.0/3.0*D$1$0*D_3_1)*T31 + (D$1_3 *D$0_2 - 1.0/3.0*D$1$0*D_3_2)*T32
    + (D$1_1 *D$0_2 - 1.0/3.0*D$1$0*D_1_2)*T21 + (D$1_1 *D$0_3 - 1.0/3.0*D$1$0*D_1_3)*T31 + (D$1_2 *D$0_3 - 1.0/3.0*D$1$0*D_2_3)*T32;

  PI20 = (D$2_0 *D$0_0 - 1.0/3.0*D$2$0*D_0_0)*T00
    + (D$2_1 *D$0_1 - 1.0/3.0*D$2$0*D_1_1)*T11 + (D$2_2 *D$0_2 - 1.0/3.0*D$2$0*D_2_2)*T22 + (D$2_3 *D$0_3 - 1.0/3.0*D$2$0*D_3_3)*T33
    + (D$2_0 *D$0_1 - 1.0/3.0*D$2$0*D_0_1)*T10 + (D$2_0 *D$0_2 - 1.0/3.0*D$2$0*D_0_2)*T20 + (D$2_0 *D$0_3 - 1.0/3.0*D$2$0*D_0_3)*T30
    + (D$2_1 *D$0_0 - 1.0/3.0*D$2$0*D_1_0)*T10 + (D$2_2 *D$0_0 - 1.0/3.0*D$2$0*D_2_0)*T20 + (D$2_3 *D$0_0 - 1.0/3.0*D$2$0*D_3_0)*T30
    + (D$2_2 *D$0_1 - 1.0/3.0*D$2$0*D_2_1)*T21 + (D$2_3 *D$0_1 - 1.0/3.0*D$2$0*D_3_1)*T31 + (D$2_3 *D$0_2 - 1.0/3.0*D$2$0*D_3_2)*T32
    + (D$2_1 *D$0_2 - 1.0/3.0*D$2$0*D_1_2)*T21 + (D$2_1 *D$0_3 - 1.0/3.0*D$2$0*D_1_3)*T31 + (D$2_2 *D$0_3 - 1.0/3.0*D$2$0*D_2_3)*T32;

  PI30 = (D$3_0 *D$0_0 - 1.0/3.0*D$3$0*D_0_0)*T00
    + (D$3_1 *D$0_1 - 1.0/3.0*D$3$0*D_1_1)*T11 + (D$3_2 *D$0_2 - 1.0/3.0*D$3$0*D_2_2)*T22 + (D$3_3 *D$0_3 - 1.0/3.0*D$3$0*D_3_3)*T33
    + (D$3_0 *D$0_1 - 1.0/3.0*D$3$0*D_0_1)*T10 + (D$3_0 *D$0_2 - 1.0/3.0*D$3$0*D_0_2)*T20 + (D$3_0 *D$0_3 - 1.0/3.0*D$3$0*D_0_3)*T30
    + (D$3_1 *D$0_0 - 1.0/3.0*D$3$0*D_1_0)*T10 + (D$3_2 *D$0_0 - 1.0/3.0*D$3$0*D_2_0)*T20 + (D$3_3 *D$0_0 - 1.0/3.0*D$3$0*D_3_0)*T30
    + (D$3_2 *D$0_1 - 1.0/3.0*D$3$0*D_2_1)*T21 + (D$3_3 *D$0_1 - 1.0/3.0*D$3$0*D_3_1)*T31 + (D$3_3 *D$0_2 - 1.0/3.0*D$3$0*D_3_2)*T32
    + (D$3_1 *D$0_2 - 1.0/3.0*D$3$0*D_1_2)*T21 + (D$3_1 *D$0_3 - 1.0/3.0*D$3$0*D_1_3)*T31 + (D$3_2 *D$0_3 - 1.0/3.0*D$3$0*D_2_3)*T32;

  PI21 = (D$2_0 *D$1_0 - 1.0/3.0*D$2$1*D_0_0)*T00
    + (D$2_1 *D$1_1 - 1.0/3.0*D$2$1*D_1_1)*T11 + (D$2_2 *D$1_2 - 1.0/3.0*D$2$1*D_2_2)*T22 + (D$2_3 *D$1_3 - 1.0/3.0*D$2$1*D_3_3)*T33
    + (D$2_0 *D$1_1 - 1.0/3.0*D$2$1*D_0_1)*T10 + (D$2_0 *D$1_2 - 1.0/3.0*D$2$1*D_0_2)*T20 + (D$2_0 *D$1_3 - 1.0/3.0*D$2$1*D_0_3)*T30
    + (D$2_1 *D$1_0 - 1.0/3.0*D$2$1*D_1_0)*T10 + (D$2_2 *D$1_0 - 1.0/3.0*D$2$1*D_2_0)*T20 + (D$2_3 *D$1_0 - 1.0/3.0*D$2$1*D_3_0)*T30
    + (D$2_2 *D$1_1 - 1.0/3.0*D$2$1*D_2_1)*T21 + (D$2_3 *D$1_1 - 1.0/3.0*D$2$1*D_3_1)*T31 + (D$2_3 *D$1_2 - 1.0/3.0*D$2$1*D_3_2)*T32
    + (D$2_1 *D$1_2 - 1.0/3.0*D$2$1*D_1_2)*T21 + (D$2_1 *D$1_3 - 1.0/3.0*D$2$1*D_1_3)*T31 + (D$2_2 *D$1_3 - 1.0/3.0*D$2$1*D_2_3)*T32;

  PI31 = (D$3_0 *D$1_0 - 1.0/3.0*D$3$1*D_0_0)*T00
    + (D$3_1 *D$1_1 - 1.0/3.0*D$3$1*D_1_1)*T11 + (D$3_2 *D$1_2 - 1.0/3.0*D$3$1*D_2_2)*T22 + (D$3_3 *D$1_3 - 1.0/3.0*D$3$1*D_3_3)*T33
    + (D$3_0 *D$1_1 - 1.0/3.0*D$3$1*D_0_1)*T10 + (D$3_0 *D$1_2 - 1.0/3.0*D$3$1*D_0_2)*T20 + (D$3_0 *D$1_3 - 1.0/3.0*D$3$1*D_0_3)*T30
    + (D$3_1 *D$1_0 - 1.0/3.0*D$3$1*D_1_0)*T10 + (D$3_2 *D$1_0 - 1.0/3.0*D$3$1*D_2_0)*T20 + (D$3_3 *D$1_0 - 1.0/3.0*D$3$1*D_3_0)*T30
    + (D$3_2 *D$1_1 - 1.0/3.0*D$3$1*D_2_1)*T21 + (D$3_3 *D$1_1 - 1.0/3.0*D$3$1*D_3_1)*T31 + (D$3_3 *D$1_2 - 1.0/3.0*D$3$1*D_3_2)*T32
    + (D$3_1 *D$1_2 - 1.0/3.0*D$3$1*D_1_2)*T21 + (D$3_1 *D$1_3 - 1.0/3.0*D$3$1*D_1_3)*T31 + (D$3_2 *D$1_3 - 1.0/3.0*D$3$1*D_2_3)*T32;

  PI32 = (D$3_0 *D$2_0 - 1.0/3.0*D$3$2*D_0_0)*T00
    + (D$3_1 *D$2_1 - 1.0/3.0*D$3$2*D_1_1)*T11 + (D$3_2 *D$2_2 - 1.0/3.0*D$3$2*D_2_2)*T22 + (D$3_3 *D$2_3 - 1.0/3.0*D$3$2*D_3_3)*T33
    + (D$3_0 *D$2_1 - 1.0/3.0*D$3$2*D_0_1)*T10 + (D$3_0 *D$2_2 - 1.0/3.0*D$3$2*D_0_2)*T20 + (D$3_0 *D$2_3 - 1.0/3.0*D$3$2*D_0_3)*T30
    + (D$3_1 *D$2_0 - 1.0/3.0*D$3$2*D_1_0)*T10 + (D$3_2 *D$2_0 - 1.0/3.0*D$3$2*D_2_0)*T20 + (D$3_3 *D$2_0 - 1.0/3.0*D$3$2*D_3_0)*T30
    + (D$3_2 *D$2_1 - 1.0/3.0*D$3$2*D_2_1)*T21 + (D$3_3 *D$2_1 - 1.0/3.0*D$3$2*D_3_1)*T31 + (D$3_3 *D$2_2 - 1.0/3.0*D$3$2*D_3_2)*T32
    + (D$3_1 *D$2_2 - 1.0/3.0*D$3$2*D_1_2)*T21 + (D$3_1 *D$2_3 - 1.0/3.0*D$3$2*D_1_3)*T31 + (D$3_2 *D$2_3 - 1.0/3.0*D$3$2*D_2_3)*T32;

}




void hydroParticleType::lanEck_deltamunu(const double  v_x, const double  v_y, const double v_z, double & D_0_0, double & D$0_0, double & D_0$0, double & D$0$0, double & D_1_1, double & D$1$1, double & D_2_2, double & D$2$2, double & D_3_3, double & D$3$3, double & D$1_1, double & D$2_2, double & D$3_3, double & D$1_0, double & D$1$0, double & D$0$1, double & D$2_0, double & D$2$0, double & D$0$2, double & D$3_0, double & D$3$0, double & D$0$3, double & D$0_1, double & D_1_0, double & D_0_1, double & D$0_2, double & D_2_0, double & D_0_2, double & D$0_3, double & D_3_0, double & D_0_3, double & D_2_1, double & D$2$1, double & D_1_2, double & D$1$2, double & D_3_1, double & D$3$1, double & D_1_3, double & D$1$3, double & D_3_2, double & D$3$2, double & D_2_3, double & D$2$3, double & D$2_1, double & D$1_2, double & D$3_1, double & D$1_3, double &  D$3_2, double & D$2_3 )
{
  double gamma = lanEck_gamma(v_x, v_y, v_z);
                
  D_0_0 = D$0_0 = D_0$0 = D$0$0 = 1.0 - pow(gamma,2);

  D_1_1 = D$1$1 = -1.0 - pow(gamma,2) * pow(v_x,2);
  D_2_2 = D$2$2 = -1.0 - pow(gamma,2) * pow(v_y,2);
  D_3_3 = D$3$3 = -1.0 - pow(gamma,2) * pow(v_z,2);

  D$1_1 = 1.0 + pow(gamma,2) * pow(v_x,2);
  D$2_2 = 1.0 + pow(gamma,2) * pow(v_y,2);
  D$3_3 = 1.0 + pow(gamma,2) * pow(v_z,2);

  D$1_0 = D$1$0 = D$0$1 = - pow(gamma,2) * v_x;
  D$2_0 = D$2$0 = D$0$2 = - pow(gamma,2) * v_y;
  D$3_0 = D$3$0 = D$0$3 = - pow(gamma,2) * v_z;

  D$0_1 = D_1_0 = D_0_1 = pow(gamma,2) * v_x;
  D$0_2 = D_2_0 = D_0_2 = pow(gamma,2) * v_y;
  D$0_3 = D_3_0 = D_0_3 = pow(gamma,2) * v_z;

  D_2_1 = D$2$1 = D_1_2 = D$1$2 = - pow(gamma,2) * v_x * v_y;
  D_3_1 = D$3$1 = D_1_3 = D$1$3 = - pow(gamma,2) * v_x * v_z;
  D_3_2 = D$3$2 = D_2_3 = D$2$3 = - pow(gamma,2) * v_y * v_z;

  D$2_1 = D$1_2 = pow(gamma,2) * v_x * v_y;
  D$3_1 = D$1_3 = pow(gamma,2) * v_x * v_z;
  D$3_2 = D$2_3 = pow(gamma,2) * v_y * v_z;
}
              
double hydroParticleType::lanEck_isotropicPressure( const double T00, const double T11, const double T22, const double T33, const double T10, const double T20, const double T30, const double T21, const double T31, const double T32, const double v_x, const double v_y, const double v_z)
{
  //calculate isotropic pressure in Landau or Eckart frame          
  double isoPressure;
    
  //------------------------------------------------------------------------------------------------
  //-----------------------------------------
  double D_0_0,D$0_0,D_0$0,D$0$0;
  double D_1_1,D$1$1,D_2_2,D$2$2,D_3_3,D$3$3;
  double D$1_1,D$2_2,D$3_3;
  double D$1_0,D$1$0,D$0$1,D$2_0,D$2$0,D$0$2,D$3_0,D$3$0,D$0$3;
  double D$0_1,D_1_0,D_0_1,D$0_2,D_2_0,D_0_2,D$0_3,D_3_0,D_0_3;
  double D_2_1,D$2$1,D_1_2,D$1$2,D_3_1,D$3$1,D_1_3,D$1$3,D_3_2,D$3$2,D_2_3,D$2$3;
  double D$2_1,D$1_2,D$3_1,D$1_3,D$3_2,D$2_3;

  //-----------
  lanEck_deltamunu( v_x, v_y, v_z, D_0_0, D$0_0, D_0$0, D$0$0, D_1_1, D$1$1, D_2_2, D$2$2, D_3_3, D$3$3, D$1_1, D$2_2, D$3_3, D$1_0, D$1$0, D$0$1, D$2_0, D$2$0, D$0$2, D$3_0, D$3$0, D$0$3, D$0_1, D_1_0, D_0_1, D$0_2, D_2_0, D_0_2, D$0_3, D_3_0, D_0_3, D_2_1, D$2$1, D_1_2, D$1$2, D_3_1, D$3$1, D_1_3, D$1$3, D_3_2, D$3$2, D_2_3, D$2$3, D$2_1, D$1_2, D$3_1, D$1_3,  D$3_2, D$2_3 );      
  //------------     
  
  //-----------------------------------------
  //pressure
  //-----------------------------------------

  isoPressure = - 1.0/3.0*(D_0_0*T00 + D_1_1*T11 + D_2_2*T22 + D_3_3*T33) 
    - 2.0/3.0*( D_1_0*T10 + D_2_0*T20 + D_3_0*T30 + D_2_1*T21 + D_3_1*T31 + D_3_2*T32);
  //----------------------------------------- 
              
  return isoPressure;
}

void hydroParticleType::lanEck_fourVelocity(const double vx, const double vy, const double vz, double& u$0, double& u_0, double& u$1, double& u_1, double& u$2, double& u_2, double& u$3, double& u_3 )
{
  double gamma;
    
  gamma = lanEck_gamma( vx, vy, vz);
  u$0 = u_0 = gamma;
  u$1 = gamma * vx;
  u$2 = gamma * vy;
  u$3 = gamma * vz;
  u_1 = - gamma * vx;
  u_2 = - gamma * vy;
  u_3 = - gamma * vz;
}
              

void hydroParticleType::lanEck_energyMomentumFlow(const double T00, const double T11, const double T22, const double T33, const double T10, const double T20, const double T30, const double T21, const double T31, const double T32, const double v_x, const double v_y, const double v_z, double& W0, double& W1, double& W2, double& W3)
{
  double T01,T02,T03,T12,T13,T23;
            
  T01 = T10;
  T02 = T20;
  T03 = T30;
  T12 = T21;
  T13 = T31;
  T23 = T32;  
  
  //------------------------------------------------------------------------------------------------
  //-----------------------------------------
  double D_0_0,D$0_0,D_0$0,D$0$0;
  double D_1_1,D$1$1,D_2_2,D$2$2,D_3_3,D$3$3;
  double D$1_1,D$2_2,D$3_3;
  double D$1_0,D$1$0,D$0$1,D$2_0,D$2$0,D$0$2,D$3_0,D$3$0,D$0$3;
  double D$0_1,D_1_0,D_0_1,D$0_2,D_2_0,D_0_2,D$0_3,D_3_0,D_0_3;
  double D_2_1,D$2$1,D_1_2,D$1$2,D_3_1,D$3$1,D_1_3,D$1$3,D_3_2,D$3$2,D_2_3,D$2$3;
  double D$2_1,D$1_2,D$3_1,D$1_3,D$3_2,D$2_3;

  //-----------
  lanEck_deltamunu( v_x, v_y, v_z, D_0_0, D$0_0, D_0$0, D$0$0, D_1_1, D$1$1, D_2_2, D$2$2, D_3_3, D$3$3, D$1_1, D$2_2, D$3_3, D$1_0, D$1$0, D$0$1, D$2_0, D$2$0, D$0$2, D$3_0, D$3$0, D$0$3, D$0_1, D_1_0, D_0_1, D$0_2, D_2_0, D_0_2, D$0_3, D_3_0, D_0_3, D_2_1, D$2$1, D_1_2, D$1$2, D_3_1, D$3$1, D_1_3, D$1$3, D_3_2, D$3$2, D_2_3, D$2$3, D$2_1, D$1_2, D$3_1, D$1_3,  D$3_2, D$2_3 );      
  //------------
              
  double u$0,u$1,u$2,u$3,u_0,u_1,u_2,u_3;
  //------------      
  lanEck_fourVelocity(v_x, v_y, v_z, u$0, u_0, u$1, u_1, u$2, u_2, u$3, u_3 );
  //------------      

  W0 = D$0_0 * T00 * u_0 + D$0_0 * T01 * u_1 + D$0_0 * T02 * u_2 + D$0_0 * T03 * u_3
    +D$0_1 * T10 * u_0 + D$0_1 * T11 * u_1 + D$0_1 * T12 * u_2 + D$0_1 * T13 * u_3
    +D$0_2 * T20 * u_0 + D$0_2 * T21 * u_1 + D$0_2 * T22 * u_2 + D$0_2 * T23 * u_3
    +D$0_3 * T30 * u_0 + D$0_3 * T31 * u_1 + D$0_3 * T32 * u_2 + D$0_3 * T33 * u_3;
              
  W1 = D$1_0 * T00 * u_0 + D$1_0 * T01 * u_1 + D$1_0 * T02 * u_2 + D$1_0 * T03 * u_3
    +D$1_1 * T10 * u_0 + D$1_1 * T11 * u_1 + D$1_1 * T12 * u_2 + D$1_1 * T13 * u_3
    +D$1_2 * T20 * u_0 + D$1_2 * T21 * u_1 + D$1_2 * T22 * u_2 + D$1_2 * T23 * u_3
    +D$1_3 * T30 * u_0 + D$1_3 * T31 * u_1 + D$1_3 * T32 * u_2 + D$1_3 * T33 * u_3;

  W2 = D$2_0 * T00 * u_0 + D$2_0 * T01 * u_1 + D$2_0 * T02 * u_2 + D$2_0 * T03 * u_3
    +D$2_1 * T10 * u_0 + D$2_1 * T11 * u_1 + D$2_1 * T12 * u_2 + D$2_1 * T13 * u_3
    +D$2_2 * T20 * u_0 + D$2_2 * T21 * u_1 + D$2_2 * T22 * u_2 + D$2_2 * T23 * u_3
    +D$2_3 * T30 * u_0 + D$2_3 * T31 * u_1 + D$2_3 * T32 * u_2 + D$2_3 * T33 * u_3;
                                      
  W3 = D$3_0 * T00 * u_0 + D$3_0 * T01 * u_1 + D$3_0 * T02 * u_2 + D$3_0 * T03 * u_3
    +D$3_1 * T10 * u_0 + D$3_1 * T11 * u_1 + D$3_1 * T12 * u_2 + D$3_1 * T13 * u_3
    +D$3_2 * T20 * u_0 + D$3_2 * T21 * u_1 + D$3_2 * T22 * u_2 + D$3_2 * T23 * u_3
    +D$3_3 * T30 * u_0 + D$3_3 * T31 * u_1 + D$3_3 * T32 * u_2 + D$3_3 * T33 * u_3;
                                       
}

void hydroParticleType::lanEck_particleFlow(const double N0, const double N1, const double N2, const double N3, const double v_x, const double v_y, const double v_z, double& V0, double& V1, double& V2, double& V3)
{
  //------------------------------------------------------------------------------------------------
  //-----------------------------------------
  double D_0_0,D$0_0,D_0$0,D$0$0;
  double D_1_1,D$1$1,D_2_2,D$2$2,D_3_3,D$3$3;
  double D$1_1,D$2_2,D$3_3;
  double D$1_0,D$1$0,D$0$1,D$2_0,D$2$0,D$0$2,D$3_0,D$3$0,D$0$3;
  double D$0_1,D_1_0,D_0_1,D$0_2,D_2_0,D_0_2,D$0_3,D_3_0,D_0_3;
  double D_2_1,D$2$1,D_1_2,D$1$2,D_3_1,D$3$1,D_1_3,D$1$3,D_3_2,D$3$2,D_2_3,D$2$3;
  double D$2_1,D$1_2,D$3_1,D$1_3,D$3_2,D$2_3;

  //-----------
  lanEck_deltamunu( v_x, v_y, v_z, D_0_0, D$0_0, D_0$0, D$0$0, D_1_1, D$1$1, D_2_2, D$2$2, D_3_3, D$3$3, D$1_1, D$2_2, D$3_3, D$1_0, D$1$0, D$0$1, D$2_0, D$2$0, D$0$2, D$3_0, D$3$0, D$0$3, D$0_1, D_1_0, D_0_1, D$0_2, D_2_0, D_0_2, D$0_3, D_3_0, D_0_3, D_2_1, D$2$1, D_1_2, D$1$2, D_3_1, D$3$1, D_1_3, D$1$3, D_3_2, D$3$2, D_2_3, D$2$3, D$2_1, D$1_2, D$3_1, D$1_3,  D$3_2, D$2_3 );      
  //------------

  V0 = (N0*D$0$0 - N1*D$1$0 - N2*D$2$0 - N3*D$3$0);
  V1 = (N0*D$0$1 - N1*D$1$1 - N2*D$2$1 - N3*D$3$1);
  V2 = (N0*D$0$2 - N1*D$1$2 - N2*D$2$2 - N3*D$3$2);
  V3 = (N0*D$0$3 - N1*D$1$3 - N2*D$2$3 - N3*D$3$3);
}

void hydroParticleType::lanEck_heatFlow(const double W0, const double W1, const double W2, const double W3, const double V0, const double V1, const double V2, const double V3, const double energyDensity, const double pressure, const double particleDensity, double& q0, double& q1, double& q2, double& q3)
{
  double h = (energyDensity + pressure) / particleDensity;//enthalpy per particle

  q0 = W0 - h * V0;
  q1 = W1 - h * V1;
  q2 = W2 - h * V2;
  q3 = W3 - h * V3;  
}


double hydroParticleType::lanEck_PImunuPImunu( const double PI00, const double PI11, const double PI22, const double PI33, const double PI10, const double PI20, const double PI30, const double PI21, const double PI31, const double PI32, const double vx, const double vy, const double vz)
{
  double PI01,PI02,PI03,PI12,PI13,PI23;
  PI01 = PI10;
  PI02 = PI20;
  PI03 = PI30;
  PI12 = PI21;
  PI13 = PI31;
  PI23 = PI32;  
  
  double PImunuPImunu = pow( PI00,2 ) - pow( PI10,2 ) - pow( PI20,2 ) - pow( PI30,2 )
    - pow( PI01,2 ) + pow( PI11,2 ) + pow( PI21,2 ) + pow( PI31,2 )
    - pow( PI02,2 ) + pow( PI12,2 ) + pow( PI22,2 ) + pow( PI32,2 )
    - pow( PI03,2 ) + pow( PI13,2 ) + pow( PI23,2 ) + pow( PI33,2 ); 
           
  return PImunuPImunu;
}

void hydroParticleType::lanEck_getAverageTemperatureAndFugacity(  const vector<hydroParticleTypeProperties> & vecType, vector<double>& T, vector<double>& fug, const vector<double> & n )
{  
  int nTypes = vecType.size();

  //temperature
  do
  {
    double sum_nT = 0.0;

    for(int i = 0; i < nTypes; i++)
    {
      sum_nT += n[i] * T[i];
    }

    T[nTypes] = sum_nT / n[nTypes];   
    
  } while(false);
  
  if(FPT_COMP_G(T[nTypes],0.0))
  {
    //fugacity
    do
    {
      double sum_nEq = 0.0;

      for(int i = 0; i < nTypes; i++)
      {
        sum_nEq += chooseEqParticleDensity(T[nTypes],1.0,vecType,i) / pow(0.197,3);
      }
      
      fug[nTypes] = n[nTypes] / sum_nEq;
      
    } while(false);  
  }

}


void hydroParticleType::giveHydroObservablesInLandauFrame( const vector<hydroParticleTypeProperties> & vecType,
							   const vector<double> & array_T00, const vector<double> & array_T11, const vector<double> & array_T22,
							   const vector<double> & array_T33, const vector<double> & array_T10, const vector<double> & array_T20,
							   const vector<double> & array_T30, const vector<double> & array_T21, const vector<double> & array_T31,
							   const vector<double> & array_T32, const vector<double> & array_N0, const vector<double> & array_N1,
							   const vector<double> & array_N2, const vector<double> & array_N3,
							   vector<double> & array_energyDensity, vector<double> & array_particleDensity, vector<double> & array_entropyDensity,
							   vector<double> & array_temperature,   vector<double> & array_fugacity,        vector<double> & array_eq_particleDensity,
							   vector<double> & array_isoPressure,   vector<double> & array_eqPressure,      vector<double> & array_bulkPressure,
							   vector<double> & array_vx,            vector<double> & array_vy,              vector<double> & array_vz,
							   vector<double> & array_vv,            vector<double> & array_gamma,
							   vector<double> & array_PI00, vector<double> & array_PI11, vector<double> & array_PI22,
							   vector<double> & array_PI33, vector<double> & array_PI10, vector<double> & array_PI20,
							   vector<double> & array_PI30, vector<double> & array_PI21, vector<double> & array_PI31,
							   vector<double> & array_PI32, vector<double> & array_PImunuPImunu,
							   vector<double> & array_q0, vector<double> & array_q1,
							   vector<double> & array_q2, vector<double> & array_q3,
							   vector<double> & array_W0, vector<double> & array_W1,
							   vector<double> & array_W2, vector<double> & array_W3,
							   vector<double> & array_V0, vector<double> & array_V1,
							   vector<double> & array_V2, vector<double> & array_V3)
{
  int nTypes = vecType.size(); 
  
  //--------------------------------// 
  //this is the summation of all components
  //--------------------------------// 
  int completeType = nTypes;

  for(int type = completeType; type >= 0; type--)
  { 

    double T00, T11, T22, T33, T10, T20, T30, T21, T31, T32;
    double N0, N1, N2, N3;
    //------------------------------------------------------------------------------------------------//
    T00 = array_T00[type];//GeV/fm^3 from getTmunuinCell
    T11 = array_T11[type];
    T22 = array_T22[type];
    T33 = array_T33[type];

    T10 = array_T10[type];
    T20 = array_T20[type];
    T30 = array_T30[type];

    T21 = array_T21[type];
    T31 = array_T31[type];
    T32 = array_T32[type];

    N0 = array_N0[type]; //1/fm^3 from getTmunuinCell
    N1 = array_N1[type];
    N2 = array_N2[type];
    N3 = array_N3[type];
    //------------------------------------------------------------------------------------------------  
        
    do
    {
      //------------------------------------------------------------------------------------------------
      // Landau Frame
      //-----------------------------------------
      double energyDensity, particleDensity, isoPressure, temperature, fugacity;
      double vx, vy, vz, vv, gamma;
      double PI00, PI11, PI22, PI33, PI10, PI20, PI30, PI21, PI31, PI32;
      double W0, W1, W2, W3, V0, V1, V2, V3, q0, q1, q2, q3;
      double PImunuPImunu;
    
      //-------------------------------//
      if(type == completeType)
      {
        bool solution;
        energyDensity = landau_energyDensity(T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, solution);//GeV/fm^3
                
        //---------
        if( std::isinf( energyDensity ) || std::isnan( energyDensity ) || FPT_COMP_LE(energyDensity,0.0) || solution == false )
        {
          break;          
        }
        //---------

        landau_velocity(T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, energyDensity, vx, vy, vz);
        vv = lanEck_velocityAbsolute(vx, vy, vz);
        gamma = lanEck_gamma(vx, vy, vz);
                
        //---------
        if ( std::isinf( vx ) || std::isnan( vx ) || FPT_COMP_GE(vx,1.0)
            || std::isinf( vy ) || std::isnan( vy ) || FPT_COMP_GE(vy,1.0)
            || std::isinf( vz ) || std::isnan( vz ) || FPT_COMP_GE(vz,1.0)
            || std::isinf( vv ) || std::isnan( vv ) || FPT_COMP_GE(vv,1.0) )
        {
          break;         
        }  
        //---------
      }
      else
      {
        vx = array_vx[completeType];
        vy = array_vy[completeType];
        vz = array_vz[completeType];
                
        //this one is the eckart definition, but for multi component system we do not calculate
        //the energy density using the landau method. Using the velocities from the whole system,
        //we get the energy density via the normal "eckart" definition
        energyDensity = eckart_energyDensity( T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, vx, vy, vz );
              
        //---------  
        if( std::isinf( energyDensity ) || std::isnan( energyDensity ) || FPT_COMP_LE(energyDensity,0.0) ){break;}
        //---------
      }
      //-------------------------------//
        
        

      //-------------------------------//      
      particleDensity = lanEck_particleDensity( N0, N1, N2, N3, vx, vy, vz);//1/fm^3
      isoPressure = lanEck_isotropicPressure( T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, vx, vy, vz);
      lanEck_shearStress( T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, vx, vy, vz, PI00, PI11, PI22, PI33, PI10, PI20, PI30, PI21, PI31, PI32); // Gev/fm^3
      PImunuPImunu = lanEck_PImunuPImunu( PI00, PI11, PI22, PI33, PI10, PI20, PI30, PI21, PI31, PI32, vx, vy, vz);    
      lanEck_energyMomentumFlow( T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, vx, vy, vz, W0, W1, W2, W3);
      lanEck_particleFlow( N0, N1, N2, N3, vx, vy, vz, V0, V1, V2, V3);
      lanEck_heatFlow( W0, W1, W2, W3, V0, V1, V2, V3, energyDensity, isoPressure, particleDensity, q0, q1, q2, q3);  
      //-------------------------------//
        
      //-------------------------------//
      //get the temperatures only for each species
      if(type == completeType){}//do nothing
      else
      {
        bool solutionTemperature;       
        temperature = chooseEqTemperature(energyDensity,particleDensity,vecType,type, solutionTemperature);//GeV
                
        if(solutionTemperature == true)
        {
          fugacity = lanEck_fugacity(particleDensity, temperature, vecType, type);
        }
        else
        {
          temperature = 0.0;
          fugacity = 0.0; 
        }
      }
      //-------------------------------//
        

      //----------------------------------------------------------------------//
      //transfer the velocity components
      if(type == completeType)
      {
        array_vx[type] = vx;
        array_vy[type] = vy;
        array_vz[type] = vz;
        array_vv[type] = vv;
        array_gamma[type] = gamma;   
      }
      else
      {
        double vx_oneSpecies, vy_oneSpecies, vz_oneSpecies, vv_oneSpecies, gamma_oneSpecies;
        landau_velocity(T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, energyDensity, vx_oneSpecies, vy_oneSpecies, vz_oneSpecies);
        vv_oneSpecies = lanEck_velocityAbsolute(vx_oneSpecies, vy_oneSpecies, vz_oneSpecies);
        gamma_oneSpecies = lanEck_gamma(vx_oneSpecies, vy_oneSpecies, vz_oneSpecies);
                
        if ( std::isinf( vx_oneSpecies ) || std::isnan( vx_oneSpecies ) || FPT_COMP_GE(vx_oneSpecies,1.0)
            || std::isinf( vy_oneSpecies ) || std::isnan( vy_oneSpecies ) || FPT_COMP_GE(vy_oneSpecies,1.0)
            || std::isinf( vz_oneSpecies ) || std::isnan( vz_oneSpecies ) || FPT_COMP_GE(vz_oneSpecies,1.0)
            || std::isinf( vv_oneSpecies ) || std::isnan( vv_oneSpecies ) || FPT_COMP_GE(vv_oneSpecies,1.0) )
        {
          //do nothing
        }
        else
        {
          array_vx[type] = vx_oneSpecies;
          array_vy[type] = vy_oneSpecies;
          array_vz[type] = vz_oneSpecies;
          array_vv[type] = vv_oneSpecies;
          array_gamma[type] = gamma_oneSpecies;
        }
      }
      //----------------------------------------------------------------------//

      //----------------------------------------------------------------------//
      //transfer the temperatures
      if(type == completeType){}//do nothing
      else
      {
        array_temperature[type] = temperature;
        array_fugacity[type] = fugacity;  
      }
      //----------------------------------------------------------------------//  
        
      //----------------------------------------------------------------------//
      //transfer all others
      array_energyDensity[type] = energyDensity;
      array_particleDensity[type] = particleDensity;
      array_isoPressure[type] = isoPressure;
        
      array_PI00[type] = PI00;
      array_PI11[type] = PI11;
      array_PI22[type] = PI22;
      array_PI33[type] = PI33;
      array_PI10[type] = PI10;
      array_PI20[type] = PI20;
      array_PI30[type] = PI30;
      array_PI21[type] = PI21;
      array_PI31[type] = PI31;
      array_PI32[type] = PI32;
        
      array_PImunuPImunu[type] = PImunuPImunu;
        
      array_V0[type] = V0;
      array_V1[type] = V1;
      array_V2[type] = V2;
      array_V3[type] = V3;
        
      array_W0[type] = W0;
      array_W1[type] = W1;
      array_W2[type] = W2;
      array_W3[type] = W3;
        
      array_q0[type] = q0;
      array_q1[type] = q1;
      array_q2[type] = q2;
      array_q3[type] = q3;    
      //----------------------------------------------------------------------//
    } while(false);
  }


  // get the temperature and fugacity for the whole system
  int nTypesPlusComplete = nTypes + 1;
          
  do
  {
    //------------------------------------------------------------------------------------------------
    // Landau Frame
    //-----------------------------------------
    vector<double> temperature(nTypesPlusComplete,0.0); 
    vector<double> particleDensity(nTypesPlusComplete,0.0);  
    vector<double> fugacity(nTypesPlusComplete,0.0);   
            
    for(int type = 0; type <= nTypes; type++)
    {
      temperature[type] = array_temperature[type];
      particleDensity[type] = array_particleDensity[type];
      fugacity[type] = array_fugacity[type];
    }
    //-----------------------------------------------//    
            
    //-----------------------------------------------//
    // calculate energy density and particle density
    if( FPT_COMP_G(particleDensity[nTypes],0.0) )
    {
      lanEck_getAverageTemperatureAndFugacity( vecType, temperature, fugacity, particleDensity );
    }
    //-----------------------------------------------//

    //-----------------------------------------------//  
    // transfer
    array_temperature[completeType] = temperature[completeType];
    array_fugacity[completeType] = fugacity[completeType];
    //-----------------------------------------------//             
  }while(false);
  
  
  
  
  //for each particle species
  for ( int type = 0; type <= nTypes; type++ )
  { 
    //------------------------------------------------------------------------------------------------
    // Landau Frame
    //-----------------------------------------
    double entropyDensity = 0.0, eqPressure = 0.0, bulkPressure = 0.0, eqDensity=0.0;
        
    //--------------------------
    //-- these are the individual ones
    //      double particleDensity = array_particleDensity[type];
    double isoPressure = array_isoPressure[type];  
        
    //-- here we use only the complete type ones
    double temperature = array_temperature[completeType];
    double fugacity = array_fugacity[completeType];
        
    if ( std::isinf( temperature ) || std::isnan( temperature ) || FPT_COMP_LE(temperature,0.0)
          || std::isinf( fugacity ) || std::isnan( fugacity )    || FPT_COMP_GE(fugacity,1000.0) )
    {
      break;
    }

    if(type == completeType)
    {
      for(int jj = 0; jj < nTypes; jj++)
      {
        double tempIsoPressure = array_isoPressure[jj];
        double tempEqPressure = chooseEqPressure(temperature,fugacity,vecType,jj) / pow(0.197,3);//GeV/fm^3

        entropyDensity += chooseEqEntropyDensity(temperature,fugacity,vecType,jj) / pow(0.197,3);//1/fm^3
        eqPressure += chooseEqPressure(temperature,fugacity,vecType,jj) / pow(0.197,3);//GeV/fm^3
        bulkPressure += lanEck_bulkPressure( tempIsoPressure, tempEqPressure);
        eqDensity += chooseEqParticleDensity(temperature,1.0,vecType,jj) / pow(0.197,3);//1/fm^3
      }   
    }
    else
    {
      entropyDensity = chooseEqEntropyDensity(temperature,fugacity,vecType,type) / pow(0.197,3);//1/fm^3
      eqPressure = chooseEqPressure(temperature,fugacity,vecType,type) / pow(0.197,3);//GeV/fm^3
      bulkPressure = lanEck_bulkPressure(isoPressure, eqPressure);
      eqDensity = chooseEqParticleDensity(temperature,1.0,vecType,type) / pow(0.197,3);//1/fm^3
    }
        
    //transfer
    array_entropyDensity[type] = entropyDensity;
    array_eqPressure[type] = eqPressure;  
    array_bulkPressure[type] = bulkPressure;
    array_eq_particleDensity[type] = eqDensity;//1/fm^3
  }
}



void hydroParticleType::giveHydroObservablesInEckartFrame( const vector<hydroParticleTypeProperties> & vecType,
							   const vector<double> & array_T00, const vector<double> & array_T11, const vector<double> & array_T22,
							   const vector<double> & array_T33, const vector<double> & array_T10, const vector<double> & array_T20,
							   const vector<double> & array_T30, const vector<double> & array_T21, const vector<double> & array_T31,
							   const vector<double> & array_T32, const vector<double> & array_N0, const vector<double> & array_N1,
							   const vector<double> & array_N2, const vector<double> & array_N3,
							   vector<double> & array_energyDensity, vector<double> & array_particleDensity, vector<double> & array_entropyDensity,
							   vector<double> & array_temperature,   vector<double> & array_fugacity,        vector<double> & array_eq_particleDensity,
							   vector<double> & array_isoPressure,   vector<double> & array_eqPressure,      vector<double> & array_bulkPressure,
							   vector<double> & array_vx,            vector<double> & array_vy,              vector<double> & array_vz,
							   vector<double> & array_vv,            vector<double> & array_gamma,
							   vector<double> & array_PI00, vector<double> & array_PI11, vector<double> & array_PI22,
							   vector<double> & array_PI33, vector<double> & array_PI10, vector<double> & array_PI20,
							   vector<double> & array_PI30, vector<double> & array_PI21, vector<double> & array_PI31,
							   vector<double> & array_PI32, vector<double> & array_PImunuPImunu,
							   vector<double> & array_q0, vector<double> & array_q1,
							   vector<double> & array_q2, vector<double> & array_q3,
							   vector<double> & array_W0, vector<double> & array_W1,
							   vector<double> & array_W2, vector<double> & array_W3,
							   vector<double> & array_V0, vector<double> & array_V1,
							   vector<double> & array_V2, vector<double> & array_V3)
{
  int nTypes = vecType.size(); 
   
  //--------------------------------// 
  //this is the summation of all components
  //--------------------------------// 
  int completeType = nTypes;

  for(int type = nTypes; type >= 0; type--){ 

    double T00, T11, T22, T33, T10, T20, T30, T21, T31, T32;
    double N0, N1, N2, N3;
    //------------------------------------------------------------------------------------------------//
    T00 = array_T00[type];//1/fm^3
    T11 = array_T11[type];
    T22 = array_T22[type];
    T33 = array_T33[type];

    T10 = array_T10[type];
    T20 = array_T20[type];
    T30 = array_T30[type];

    T21 = array_T21[type];
    T31 = array_T31[type];
    T32 = array_T32[type];

    N0 = array_N0[type]; //1/fm^3
    N1 = array_N1[type];
    N2 = array_N2[type];
    N3 = array_N3[type];
    //------------------------------------------------------------------------------------------------  

        
    do
    {
      //------------------------------------------------------------------------------------------------
      // Eckart Frame
      //-----------------------------------------
      double energyDensity, particleDensity, isoPressure, temperature, fugacity;
      double vx, vy, vz, vv, gamma;
      double PI00, PI11, PI22, PI33, PI10, PI20, PI30, PI21, PI31, PI32;
      double W0, W1, W2, W3, V0, V1, V2, V3, q0, q1, q2, q3;
      double PImunuPImunu;
        
      //-------------------------------//      
      if(type == completeType)
      {
        eckart_velocity( N0, N1, N2, N3, vx, vy, vz);
        vv = lanEck_velocityAbsolute(vx, vy, vz);
        gamma = lanEck_gamma(vx, vy, vz);
        
        //---------
	if ( std::isinf( vx ) || std::isnan( vx ) || FPT_COMP_GE(vx,1.0)
	     || std::isinf( vy ) || std::isnan( vy ) || FPT_COMP_GE(vy,1.0)
	     || std::isinf( vz ) || std::isnan( vz ) || FPT_COMP_GE(vz,1.0)
	     || std::isinf( vv ) || std::isnan( vv ) || FPT_COMP_GE(vv,1.0) )
	{break;}   
        //---------  
          
      }
      else
      {
	vx = array_vx[completeType];
	vy = array_vy[completeType];
	vz = array_vz[completeType];
      }       
      //-------------------------------//     

      //-------------------------------//
      energyDensity = eckart_energyDensity( T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, vx, vy, vz ); 
      if( std::isinf( energyDensity ) || std::isnan( energyDensity ) || FPT_COMP_LE(energyDensity,0.0) ){break;}
      //-------------------------------//
        

      //-------------------------------//      
      particleDensity = lanEck_particleDensity( N0, N1, N2, N3, vx, vy, vz);
      isoPressure = lanEck_isotropicPressure( T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, vx, vy, vz);
      lanEck_shearStress( T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, vx, vy, vz, PI00, PI11, PI22, PI33, PI10, PI20, PI30, PI21, PI31, PI32); 
      PImunuPImunu = lanEck_PImunuPImunu( PI00, PI11, PI22, PI33, PI10, PI20, PI30, PI21, PI31, PI32, vx, vy, vz);    
      lanEck_energyMomentumFlow( T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, vx, vy, vz, W0, W1, W2, W3);
      lanEck_particleFlow( N0, N1, N2, N3, vx, vy, vz, V0, V1, V2, V3);
      lanEck_heatFlow( W0, W1, W2, W3, V0, V1, V2, V3, energyDensity, isoPressure, particleDensity, q0, q1, q2, q3);  
      //-------------------------------//
        
      //-------------------------------//
      //get the temperatures only for each species
      if(type == completeType){}//do nothing
      else
      {
	bool solutionTemperature;       
	temperature = chooseEqTemperature(energyDensity,particleDensity,vecType,type, solutionTemperature);//GeV
          
	if(solutionTemperature == true)
	{fugacity = lanEck_fugacity(particleDensity, temperature, vecType, type);}
	else
	{
	  temperature = 0.0;
	  fugacity = 0.0; 
	}
      }
      //-------------------------------//
        

      //----------------------------------------------------------------------//
      //transfer the velocity components
      if(type == completeType)
      {
	array_vx[type] = vx;
	array_vy[type] = vy;
	array_vz[type] = vz;
	array_vv[type] = vv;
	array_gamma[type] = gamma;   
      }
      else
      {
	double vx_oneSpecies, vy_oneSpecies, vz_oneSpecies, vv_oneSpecies, gamma_oneSpecies;
	eckart_velocity( N0, N1, N2, N3, vx_oneSpecies, vy_oneSpecies, vz_oneSpecies);
	vv_oneSpecies = lanEck_velocityAbsolute(vx_oneSpecies, vy_oneSpecies, vz_oneSpecies);
	gamma_oneSpecies = lanEck_gamma(vx_oneSpecies, vy_oneSpecies, vz_oneSpecies);
          
	if ( std::isinf( vx_oneSpecies ) || std::isnan( vx_oneSpecies ) || FPT_COMP_GE(vx_oneSpecies,1.0)
	     || std::isinf( vy_oneSpecies ) || std::isnan( vy_oneSpecies ) || FPT_COMP_GE(vy_oneSpecies,1.0)
	     || std::isinf( vz_oneSpecies ) || std::isnan( vz_oneSpecies ) || FPT_COMP_GE(vz_oneSpecies,1.0)
	     || std::isinf( vv_oneSpecies ) || std::isnan( vv_oneSpecies ) || FPT_COMP_GE(vv_oneSpecies,1.0) )
	{
	  //do nothing
	}
	else
	{
	  array_vx[type] = vx_oneSpecies;
	  array_vy[type] = vy_oneSpecies;
	  array_vz[type] = vz_oneSpecies;
	  array_vv[type] = vv_oneSpecies;
	  array_gamma[type] = gamma_oneSpecies;
	}
      }
      //----------------------------------------------------------------------//

      //----------------------------------------------------------------------//
      //transfer the temperatures
      if(type == completeType){}//do nothing
      else
      {
	array_temperature[type] = temperature;
	array_fugacity[type] = fugacity;  
      }
      //----------------------------------------------------------------------//  
        
      //----------------------------------------------------------------------//
      //transfer all others
      array_energyDensity[type] = energyDensity;
      array_particleDensity[type] = particleDensity;
      array_isoPressure[type] = isoPressure;
        
      array_PI00[type] = PI00;
      array_PI11[type] = PI11;
      array_PI22[type] = PI22;
      array_PI33[type] = PI33;
      array_PI10[type] = PI10;
      array_PI20[type] = PI20;
      array_PI30[type] = PI30;
      array_PI21[type] = PI21;
      array_PI31[type] = PI31;
      array_PI32[type] = PI32;
        
      array_PImunuPImunu[type] = PImunuPImunu;
        
      array_V0[type] = V0;
      array_V1[type] = V1;
      array_V2[type] = V2;
      array_V3[type] = V3;
        
      array_W0[type] = W0;
      array_W1[type] = W1;
      array_W2[type] = W2;
      array_W3[type] = W3;
        
      array_q0[type] = q0;
      array_q1[type] = q1;
      array_q2[type] = q2;
      array_q3[type] = q3;    
      //----------------------------------------------------------------------//
    } while(false);
  }


  //get the temperature and fugacity for the whole system
  int nTypesPlusComplete = nTypes + 1;
          
  do
  {
    //------------------------------------------------------------------------------------------------
    //Eckart Frame
    //-----------------------------------------
    vector<double> temperature(nTypesPlusComplete,0.0); 
    vector<double> particleDensity(nTypesPlusComplete,0.0);  
    vector<double> fugacity(nTypesPlusComplete,0.0);   
            
    for(int type = 0; type <= nTypes; type++)
    {
      temperature[type] = array_temperature[type];
      particleDensity[type] = array_particleDensity[type];
      fugacity[type] = array_fugacity[type];
    }
    //-----------------------------------------------//    
            
    //-----------------------------------------------//
    //calculate energy density and particle density
    if( FPT_COMP_G(particleDensity[nTypes],0.0) )
    {
      lanEck_getAverageTemperatureAndFugacity( vecType, temperature, fugacity, particleDensity );
    }
    //-----------------------------------------------//

    //-----------------------------------------------//  
    //transfer
    array_temperature[completeType] = temperature[completeType];
    array_fugacity[completeType] = fugacity[completeType];
    //-----------------------------------------------//             
  } while(false);
  
  
  
  
  //for each particle species
  for(int type = 0; type <= nTypes; type++)
  { 
    do
    {
      //------------------------------------------------------------------------------------------------
      // Landau Frame
      //-----------------------------------------
      double entropyDensity = 0.0,eqPressure = 0.0,bulkPressure = 0.0;
          
      //--------------------------
      //-- these are the individual ones
      //      double particleDensity = array_particleDensity[type];
      double isoPressure = array_isoPressure[type];  
          
      //-- here we use only the complete type ones
      double temperature = array_temperature[completeType];
      double fugacity = array_fugacity[completeType];
          
      if ( std::isinf( temperature ) || std::isnan( temperature ) || FPT_COMP_LE(temperature,0.0)
           || std::isinf( fugacity ) || std::isnan( fugacity )    || FPT_COMP_GE(fugacity,1000.0) )
      {
        break;
      }

      if(type == completeType)
      {
        for(int jj = 0; jj < nTypes; jj++)
        {
          double tempIsoPressure = array_isoPressure[jj];
          double tempEqPressure = chooseEqPressure(temperature,fugacity,vecType,jj) / pow(0.197,3);//GeV/fm^3
          
          entropyDensity += chooseEqEntropyDensity(temperature,fugacity,vecType,jj) / pow(0.197,3);//1/fm^3
          eqPressure += chooseEqPressure(temperature,fugacity,vecType,jj) / pow(0.197,3);//GeV/fm^3
          bulkPressure += lanEck_bulkPressure( tempIsoPressure, tempEqPressure);
        }   
      }
      else
      {
        entropyDensity = chooseEqEntropyDensity(temperature,fugacity,vecType,type) / pow(0.197,3);//1/fm^3
        eqPressure = chooseEqPressure(temperature,fugacity,vecType,type) / pow(0.197,3);//GeV/fm^3
        bulkPressure = lanEck_bulkPressure(isoPressure, eqPressure);
      }
          
      //transfer
      array_entropyDensity[type] = entropyDensity;
      array_eqPressure[type] = eqPressure;  
      array_bulkPressure[type] = bulkPressure;
      
      
    } while(false);
  }
        
}
