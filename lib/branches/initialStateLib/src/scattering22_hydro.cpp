//-----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/scattering22_hydro.cpp $
//$LastChangedDate: 2018-04-27 19:16:12 +0200 (Fr, 27. Apr 2018) $
//$LastChangedRevision: 2741 $
//$LastChangedBy: greif $
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


#include "scattering22_hydro.h"
#include "particleprototype.h"
#include <iostream>
using namespace std;

using std::vector;


scattering22_hydro::scattering22_hydro(){}

scattering22_hydro::~scattering22_hydro(){}

//time - timeShift
double scattering22_hydro::getVarTimeDepCS( const double time )
{
/*  
  //needed for function
  double nDensity = 1000.0;
  
  //---------------------------------------
  //---------------------------------------
  double cs11,cs12,cs22;
  double n1,n2;
  double a,b,c,d;
  double A,B,C,D;
  double lambda1,lambda2;
  double r,r0,fRatio;
  double sigma;
  
  cs11 = 2.0;
  cs12 = 1.0;
  cs22 = 0.5;

  fRatio = 1.0;
  n1 = nDensity * fRatio / (1.0 + fRatio);
  n2 = nDensity / (1.0 + fRatio);  
  
  r0 = n1/n2;
  
  a = - 5.0/9.0 * n1 * cs11 - 7.0/9.0 * n2 * cs12;
  b = + 2.0/9.0 * n1 * cs12;
  c = - 5.0/9.0 * n2 * cs22 - 7.0/9.0 * n1 * cs12;
  d = + 2.0/9.0 * n2 * cs12;
  
  lambda1 = ( a + c + sqrt( pow(a-c,2.0) + 4.0 * b * d ) ) / 2.0;
  lambda2 = ( a + c - sqrt( pow(a-c,2.0) + 4.0 * b * d ) ) / 2.0;
  
  A = ( ( lambda2 - a ) * r0 - b ) / ( lambda2 - lambda1 );
  B = ( ( a - lambda1 ) * r0 + b ) / ( lambda2 - lambda1 );
  C = A * ( ( lambda1 - a ) / b );
  D = B * ( ( lambda2 - a ) / b );
  
  r = ( A * exp( lambda1 * time ) + B * exp( lambda2 * time ) ) / ( C * exp( lambda1 * time ) + D * exp( lambda2 * time ) );
  
  sigma = 1.0 / nDensity * ( r / ( 1.0 + r) * ( n1 * cs11 + n2 * cs12 ) + 1.0 / ( 1.0 + r) * ( n2 * cs22 + n1 * cs12 ) ); 
  
  cout << "sigma " << sigma << endl;
  
  return sigma;
  */

  double sigma;

  sigma = 0.225 * tanh( -0.85 * time ) + 1.125;

  return sigma;
  
}

/**
 * Get the mixutre cross sections for ELASTIC 2->2 processes
 * 
 * @return cross section in units 1/GeV^2
 */
double scattering22_hydro::getMixtureXSection22( hydroParticleType & theHydroParticleType, const vector<hydroParticleTypeProperties> & vecType, const int F1, const int F2, const double scalingFactor, std::string & info)
{
  //here are the cross section for the mixture case and different cross section
  //all cross section are supposed to be in mb. A converter gives it back to 1/GeV²
  //if not, you have to deactivate the conversion
  double XSectionMixture;
  double cs = -1.0;
  
  double cs_GluonGluon, cs_GluonQuarkUp, cs_GluonQuarkDown, cs_GluonMassA, cs_GluonMassB, cs_GluonMasslessC, cs_GluonMasslessD;
  double cs_QuarkUpGluon, cs_QuarkUpQuarkUp, cs_QuarkUpQuarkDown, cs_QuarkUpMassA, cs_QuarkUpMassB, cs_QuarkUpMasslessC, cs_QuarkUpMasslessD;
  double cs_QuarkDownGluon, cs_QuarkDownQuarkUp, cs_QuarkDownQuarkDown, cs_QuarkDownMassA, cs_QuarkDownMassB, cs_QuarkDownMasslessC, cs_QuarkDownMasslessD;     
  double cs_MassAGluon, cs_MassAQuarkUp, cs_MassAQuarkDown, cs_MassAMassA, cs_MassAMassB, cs_MassAMasslessC, cs_MassAMasslessD;  
  double cs_MassBGluon, cs_MassBQuarkUp, cs_MassBQuarkDown, cs_MassBMassA, cs_MassBMassB, cs_MassBMasslessC, cs_MassBMasslessD;
  double cs_MasslessCGluon, cs_MasslessCQuarkUp, cs_MasslessCQuarkDown , cs_MasslessCMassA, cs_MasslessCMassB, cs_MasslessCMasslessC, cs_MasslessCMasslessD;  
  double cs_MasslessDGluon, cs_MasslessDQuarkUp, cs_MasslessDQuarkDown , cs_MasslessDMassA, cs_MasslessDMassB, cs_MasslessDMasslessC, cs_MasslessDMasslessD;  
  
  //temporäre Werte zum Arbeiten
  //Vorsicht, manchmal 
  double UpzuUp = 2.0;
  double UpzuDown = UpzuUp / 2.0;
  double DownzuDown = UpzuUp / 4.0;
  
  
  
//   double UpzuUp = 1.531;
//   double UpzuDown = 1.531;
//   double DownzuDown = 1.531;
//         
        //------------------------
  double CzuC = 2.0;
  double CzuD = CzuC / 2.0;
  double DzuD = CzuC / 4.0;     
  
  //----------------------------------------//  
  //gluons with X
  //----------------------------------------//  
  cs_GluonGluon = 1.125;
  cs_GluonQuarkUp = 2.0;
  cs_GluonQuarkDown = 2.0;      
  cs_GluonMassA = 3.0;
  cs_GluonMassB = 4.0;
  cs_GluonMasslessC = 5.0;
  cs_GluonMasslessD = 6.0;  
  //----------------------------------------//
  
  //----------------------------------------//  
  //quarks Up with X
  //----------------------------------------//
  cs_QuarkUpGluon = cs_GluonQuarkUp;//set above
  //----------------------------------------//  
  cs_QuarkUpQuarkUp = UpzuUp;
  cs_QuarkUpQuarkDown = UpzuDown;
  cs_QuarkUpMassA = 3.0;
  cs_QuarkUpMassB = 4.0;
  cs_QuarkUpMasslessC = 5.0;
  cs_QuarkUpMasslessD = 6.0;
  //----------------------------------------//  
        
  //----------------------------------------//  
  //quarks Down with X
  //----------------------------------------//
  cs_QuarkDownGluon = cs_GluonQuarkDown;//set above
  cs_QuarkDownQuarkUp = cs_QuarkUpQuarkDown;//set above  
  //----------------------------------------//  
  cs_QuarkDownQuarkDown = DownzuDown;
  cs_QuarkDownMassA = 3.0;
  cs_QuarkDownMassB = 4.0;
  cs_QuarkDownMasslessC = 5.0;
  cs_QuarkDownMasslessD = 6.0;
  //----------------------------------------//          
  
  //massA with X
  //----------------------------------------//  
  cs_MassAGluon = cs_GluonMassA;//set above
  cs_MassAQuarkUp = cs_QuarkUpMassA;//set above
  cs_MassAQuarkDown = cs_QuarkDownMassA;//set above  
  //----------------------------------------//
  cs_MassAMassA = 3.0;
  cs_MassAMassB = 4.0;  
  cs_MassAMasslessC = 5.0;
  cs_MassAMasslessD = 6.0;  
  //----------------------------------------//  

  //----------------------------------------//
  //mass B with X
  //----------------------------------------//  
  cs_MassBGluon = cs_GluonMassB;//set above
  cs_MassBQuarkUp = cs_QuarkUpMassB;//set above
  cs_MassBQuarkDown = cs_QuarkDownMassB;//set above  
  cs_MassBMassA = cs_MassAMassB;//set above
  //----------------------------------------//  
  cs_MassBMassB = 4.0; 
  cs_MassBMasslessC = 5.0;
  cs_MassBMasslessD = 6.0;
  //----------------------------------------//  

  //----------------------------------------//
  //masslessC with X
  //----------------------------------------//  
  cs_MasslessCGluon = cs_GluonMasslessC;//set above
  cs_MasslessCQuarkUp = cs_QuarkUpMasslessC;//set above
  cs_MasslessCQuarkDown = cs_QuarkDownMasslessC;//set above  
  cs_MasslessCMassA = cs_MassAMasslessC;//set above
  cs_MasslessCMassB = cs_MassBMasslessC;//set above
  //----------------------------------------//  
  cs_MasslessCMasslessC = CzuC;
  cs_MasslessCMasslessD = CzuD;  
  //----------------------------------------//  

  //----------------------------------------//
  //massless D with X
  //----------------------------------------//  
  cs_MasslessDGluon = cs_GluonMasslessD;//set above
  cs_MasslessDQuarkUp = cs_QuarkUpMasslessD;//set above
  cs_MasslessDQuarkDown = cs_QuarkDownMasslessD;//set above  
  cs_MasslessDMassA = cs_MassAMasslessD;//set above
  cs_MasslessDMassB = cs_MassBMasslessD;//set above
  cs_MasslessDMasslessC = cs_MasslessCMasslessD;//set above
  //----------------------------------------// 
  cs_MasslessDMasslessD = DzuD;    
  //----------------------------------------//  

  switch ( F1 )
  {
  case gluon:
    switch ( F2 )
    {
    case gluon:     cs = cs_GluonGluon;  break;    // gluon <-> gluon
    case up:        cs = cs_GluonQuarkUp; break;   // gluon <-> quarkUp
    case down:      cs = cs_GluonQuarkDown; break; // gluon <-> quarkDown
    case massA:     cs = cs_GluonMassA; break;     // gluon <-> massA   
    case massB:     cs = cs_GluonMassB; break;     // gluon <-> massB
    case masslessC: cs = cs_GluonMasslessC; break; // gluon <-> masslessC  
    case masslessD: cs = cs_GluonMasslessD; break; // gluon <-> masslessD  
    } break;

  case up:
    switch ( F2 )
    {
    case gluon:     cs = cs_QuarkUpGluon;  break;    // quarkUp <-> gluon
    case up:        cs = cs_QuarkUpQuarkUp; break;   // quarkUp <-> quarkUp
    case down:      cs = cs_QuarkUpQuarkDown; break; // quarkUp <-> quarkDown
    case massA:     cs = cs_QuarkUpMassA; break;     // quarkUp <-> massA   
    case massB:     cs = cs_QuarkUpMassB; break;     // quarkUp <-> massB
    case masslessC: cs = cs_QuarkUpMasslessC; break; // quarkUp <-> masslessC  
    case masslessD: cs = cs_QuarkUpMasslessD; break; // quarkUp <-> masslessD  
    } break;

  case down:
    switch ( F2 )
    {
    case gluon:     cs = cs_QuarkDownGluon;  break;    // quarkDown <-> gluon
    case up:        cs = cs_QuarkDownQuarkUp; break;   // quarkDown <-> quarkUp
    case down:      cs = cs_QuarkDownQuarkDown; break; // quarkDown <-> quarkDown
    case massA:     cs = cs_QuarkDownMassA; break;     // quarkDown <-> massA   
    case massB:     cs = cs_QuarkDownMassB; break;     // quarkDown <-> massB
    case masslessC: cs = cs_QuarkDownMasslessC; break; // quarkDown <-> masslessC  
    case masslessD: cs = cs_QuarkDownMasslessD; break; // quarkDown <-> masslessD  
    } break;

  case massA:
    switch ( F2 )
    {
    case gluon:     cs = cs_MassAGluon;  break;    // massA <-> gluon
    case up:        cs = cs_MassAQuarkUp; break;   // massA <-> quarkUp
    case down:      cs = cs_MassAQuarkDown; break; // massA <-> quarkDown
    case massA:     cs = cs_MassAMassA; break;     // massA <-> massA   
    case massB:     cs = cs_MassAMassB; break;     // massA <-> massB
    case masslessC: cs = cs_MassAMasslessC; break; // massA <-> masslessC  
    case masslessD: cs = cs_MassAMasslessD; break; // massA <-> masslessD  
    } break;

  case massB:
    switch ( F2 )
    {
    case gluon:     cs = cs_MassBGluon;  break;    // massB <-> gluon
    case up:        cs = cs_MassBQuarkUp; break;   // massB <-> quarkUp
    case down:      cs = cs_MassBQuarkDown; break; // massB <-> quarkDown
    case massA:     cs = cs_MassBMassA; break;     // massB <-> massA   
    case massB:     cs = cs_MassBMassB; break;     // massB <-> massB
    case masslessC: cs = cs_MassBMasslessC; break; // massB <-> masslessC  
    case masslessD: cs = cs_MassBMasslessD; break; // massB <-> masslessD  
    } break;

  case masslessC:
    switch ( F2 )
    {
    case gluon:     cs = cs_MasslessCGluon;  break;    // masslessC <-> gluon
    case up:        cs = cs_MasslessCQuarkUp; break;   // masslessC <-> quarkUp
    case down:      cs = cs_MasslessCQuarkDown; break; // masslessC <-> quarkDown
    case massA:     cs = cs_MasslessCMassA; break;     // masslessC <-> massA   
    case massB:     cs = cs_MasslessCMassB; break;     // masslessC <-> massB
    case masslessC: cs = cs_MasslessCMasslessC; break; // masslessC <-> masslessC  
    case masslessD: cs = cs_MasslessCMasslessD; break; // masslessC <-> masslessD 
    } break;

  case masslessD:
    switch ( F2 )
    {
    case gluon:     cs = cs_MasslessDGluon;  break;    // masslessD <-> gluon
    case up:        cs = cs_MasslessDQuarkUp; break;   // masslessD <-> quarkUp
    case down:      cs = cs_MasslessDQuarkDown; break; // masslessD <-> quarkDown
    case massA:     cs = cs_MasslessDMassA; break;     // masslessD <-> massA   
    case massB:     cs = cs_MasslessDMassB; break;     // masslessD <-> massB
    case masslessC: cs = cs_MasslessDMasslessC; break; // masslessD <-> masslessC  
    case masslessD: cs = cs_MasslessDMasslessD; break; // masslessD <-> masslessD 
    } break;
  }

  
  //--------------------------------------------------------------------//
  //depends if you want to convert mb to 1/GeV^2
  //XSectionMixture = cs / scalingFactor / (10.0 * pow(0.197,2) );//convert from mb to 1/GeV²
  XSectionMixture = cs / scalingFactor;
  //--------------------------------------------------------------------//
  

  //--------------------------------------------------------------------//      
  int typeF1 = theHydroParticleType.getParticleType(F1,vecType);
  int typeF2 = theHydroParticleType.getParticleType(F2,vecType);

  if( typeF1 >= 0 && typeF2 >= 0)
  {
    info += vecType[typeF1].nameOfParticle + " <-> " + vecType[typeF2].nameOfParticle;
  }
  else
  {
    std::string errMsg = "Error in mixture cross section: No correct flavor!";
    throw eScatt22Hydro_error( errMsg );
  }
  //--------------------------------------------------------------------//

  //--------------------------------------------------------------------//
  if(cs < 0.0)
  {
    std::stringstream aa;
    aa << cs;
    std::string errMsg = "Error in mixture cross section: Value of cross section is " + aa.str();
    throw eScatt22Hydro_error( errMsg );
  }
  //--------------------------------------------------------------------//
  
  return XSectionMixture;    
}


double scattering22_hydro::getConstMFP_XSection( hydroParticleType & theHydroParticleType, const double mfp, const double T00, const double T11, const double T22, const double T33, const double T10, const double T20, const double T30, const double T21, const double T31, const double T32, const double N0, const double N1, const double N2, const double N3, const int NN )
{  
  double particleDensity, vx, vy, vz;
  
  if(NN > 2)
  {
    bool solution;
    double energyDensity = theHydroParticleType.landau_energyDensity(T00,T11,T22,T33,T10,T20,T30,T21,T31,T32,solution);
    theHydroParticleType.landau_velocity(T00,T11,T22,T33,T10,T20,T30,T21,T31,T32,energyDensity,vx,vy,vz);    
    particleDensity = theHydroParticleType.lanEck_particleDensity(N0,N1,N2,N3,vx,vy,vz); // 1 / / fm^3
  }
  else
  {
    theHydroParticleType.eckart_velocity(N0,N1,N2,N3,vx,vy,vz);
    //    double energyDensity = theHydroParticleType.eckart_energyDensity(T00,T11,T22,T33,T10,T20,T30,T21,T31,T32,vx,vy,vz);
    particleDensity = theHydroParticleType.lanEck_particleDensity(N0,N1,N2,N3,vx,vy,vz); // 1 / / fm^3               
  }
  //---------
  double cs = 1.0 / ( particleDensity * mfp) / pow( 0.197,2 ); //1 / GeV^2
  //---------
  
  return cs;
}


//---------------------------------//
//only valid for massless particles, otherwise many corrections have to be done
double scattering22_hydro::getConstEtaOverS_XSection( hydroParticleType & theHydroParticleType, const vector<hydroParticleTypeProperties> & vecType, const double etaOverS,
                                const vector<double> & T00, const vector<double> & T11, const vector<double> & T22,
                                const vector<double> & T33, const vector<double> & T10, const vector<double> & T20,
                                const vector<double> & T30, const vector<double> & T21, const vector<double> & T31,
                                const vector<double> & T32, const vector<double> & N0, const vector<double> & N1,
                                const vector<double> & N2, const vector<double> & N3, const vector<int> & NN , double & extractedT, double & extractedN, double & extractedFug)
{  
  int nTypes = vecType.size();  
  int nTypesPlusComplete = nTypes + 1; 
  
  vector<double> energyDensity(nTypesPlusComplete,0); 
  vector<double> particleDensity(nTypesPlusComplete,0); 
  vector<double> entropyDensity(nTypesPlusComplete,0); 
  vector<double> temperature(nTypesPlusComplete,0);  
  vector<double> fugacity(nTypesPlusComplete,0); 
  vector<double> eq_particleDensity(nTypesPlusComplete,0);  
  vector<double> isoPressure(nTypesPlusComplete,0); 
  vector<double> eqPressure(nTypesPlusComplete,0); 
  vector<double> bulkPressure(nTypesPlusComplete,0);
  vector<double> vx(nTypesPlusComplete,0);    
  vector<double> vy(nTypesPlusComplete,0); 
  vector<double> vz(nTypesPlusComplete,0); 
  vector<double> vv(nTypesPlusComplete,0);
  vector<double> gamma(nTypesPlusComplete,0);
  vector<double> PI00(nTypesPlusComplete,0);    
  vector<double> PI11(nTypesPlusComplete,0); 
  vector<double> PI22(nTypesPlusComplete,0); 
  vector<double> PI33(nTypesPlusComplete,0);   
  vector<double> PI10(nTypesPlusComplete,0);    
  vector<double> PI20(nTypesPlusComplete,0); 
  vector<double> PI30(nTypesPlusComplete,0); 
  vector<double> PI21(nTypesPlusComplete,0); 
  vector<double> PI31(nTypesPlusComplete,0);    
  vector<double> PI32(nTypesPlusComplete,0); 
  vector<double> PImunuPImunu(nTypesPlusComplete,0); 
  vector<double> q0(nTypesPlusComplete,0); 
  vector<double> q1(nTypesPlusComplete,0);
  vector<double> q2(nTypesPlusComplete,0);
  vector<double> q3(nTypesPlusComplete,0);
  vector<double> W0(nTypesPlusComplete,0); 
  vector<double> W1(nTypesPlusComplete,0);
  vector<double> W2(nTypesPlusComplete,0);
  vector<double> W3(nTypesPlusComplete,0);
  vector<double> V0(nTypesPlusComplete,0); 
  vector<double> V1(nTypesPlusComplete,0);
  vector<double> V2(nTypesPlusComplete,0);
  vector<double> V3(nTypesPlusComplete,0);  
 
//   cout << T00[0] << endl;
  
  theHydroParticleType.giveHydroObservablesInLandauFrame( vecType, T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, N0, N1, N2, N3,
                                energyDensity, particleDensity, entropyDensity, temperature, fugacity, eq_particleDensity,
                                isoPressure, eqPressure, bulkPressure, vx, vy, vz, vv, gamma,
                                PI00, PI11, PI22, PI33, PI10, PI20, PI30, PI21, PI31, PI32, PImunuPImunu, q0, q1, q2, q3,
                                W0, W1, W2, W3, V0, V1, V2, V3);
  
  //particleDensity[nTypes] in 1/fm^3 !!!!!!!!!!!!!!!!!!!
  
  double cs,vRel;
  vRel = 1.0;
//cs = 1.267 * temperature[nTypes] / entropyDensity[nTypes] / etaOverS / vRel / pow(0.197,3); // 1/GeV^2
//   cs = 1.267 * temperature[nTypes] / ( 4.0 * particleDensity[nTypes] ) / etaOverS / vRel / pow(0.197,3); // 1/GeV^2  
//   cout << " w/o fugacity " <<  cs / pow(0.197,2.0) * 10. << " mb " << endl;
  double eDensity = energyDensity[nTypes];
  extractedT    = temperature[nTypes];
  extractedN    = particleDensity[nTypes];
  extractedFug  = fugacity[nTypes];
  
//   std::cout << extractedT << "\t\t" << eDensity/(3.*extractedN) << std::endl;
  
  
  //CAUTION
//   cs = 1.267 * temperature[nTypes] / ( particleDensity[nTypes]*pow(0.197,3) * ( 4.0 - log(fugacity[nTypes]) ) * etaOverS * vRel ) ; // 1/GeV^2  
  
  cs = 1.2 * temperature[nTypes] / ( particleDensity[nTypes]*pow(0.197,3) * ( 4.0 - log(fugacity[nTypes]) ) * etaOverS * vRel ) ; // 1/GeV^2  


  
//   TESTcase:
//   cs = 1.267 / ( 4. * 16./pow(M_PI,2.0)*0.4*0.4 * 0.2  ) ; // 1/GeV^2  
  
  
//   cout << cs * pow(0.197,2.0) * 10. << " mb " << endl;
  
//   cout << cs / pow(0.197,2.0) * 10. << " mb " << endl;
 
  //if(1.0/(particleDensity[nTypes] * cs) / pow(0.197,2.0) <  0.05)
//   {
//    cout << "particleDensity[nTypes] = " << particleDensity[nTypes]* pow(0.197,3.0)/(pow(temperature[nTypes],3.0))*pow(M_PI,2.0) << endl;
//    cout << "                               " << 16./pow(M_PI,2.0)*pow(temperature[nTypes],3.0)*fugacity[nTypes] << endl;
//    cout << "                               " << particleDensity[nTypes]* pow(0.197,3.0) << endl;
//    cout << "                               " << 16./pow(M_PI,2.0)*pow(temperature[nTypes],3.0)*fugacity[nTypes] / (particleDensity[nTypes]* pow(0.197,3.0)) << endl;
//    cout << "                               " << fugacity[nTypes] << endl;
//    cout << endl;
//    cout << "log(lambda) = " << log(fugacity[nTypes]) << endl; 
//    cout << "n*(4-ln lambda)/T^3 = " <<  particleDensity[nTypes]* pow(0.197,3.0)*(4.-log(fugacity[nTypes])) / (pow(temperature[nTypes],3.0)) << endl;
//    cout << "was kommt raus?" << endl;
//    cout << "NN " << NN[nTypes] << endl;
//    cout << "cs " << cs << endl; 
//    cout << "T " << temperature[nTypes] << endl;
//    cout << "fug " << fugacity[nTypes] << endl;   
//    cout << "s " << entropyDensity[nTypes] << endl;
//    cout << "etaS " << etaOverS << endl;
//    cout << "T00 " << T00[nTypes] << endl;  
//    cout << "T33 " << T33[nTypes] << endl;    
//    cout << "N0 " << N0[nTypes] << endl;
//    cout << "v " << vv[nTypes] << endl;  
//    cout << "lambda = " << 1.0/(particleDensity[nTypes] * cs) / pow(0.197,2.0) << endl;
//   }

  return cs;
}

double scattering22_hydro::getScalingFactorLambdaT2( hydroParticleType & theHydroParticleType, const vector<hydroParticleTypeProperties> & vecType,
                                const vector<double> & T00, const vector<double> & T11, const vector<double> & T22,
                                const vector<double> & T33, const vector<double> & T10, const vector<double> & T20,
                                const vector<double> & T30, const vector<double> & T21, const vector<double> & T31,
                                const vector<double> & T32, const vector<double> & N0, const vector<double> & N1,
                                const vector<double> & N2, const vector<double> & N3, const vector<int> & NN )
{
  int nTypes = vecType.size();  
  int nTypesPlusComplete = nTypes + 1;  
  double scalingFactorLambdaT2;
  
  vector<double> energyDensity(nTypesPlusComplete,0); 
  vector<double> particleDensity(nTypesPlusComplete,0); 
  vector<double> entropyDensity(nTypesPlusComplete,0); 
  vector<double> temperature(nTypesPlusComplete,0);  
  vector<double> fugacity(nTypesPlusComplete,0); 
  vector<double> eq_particleDensity(nTypesPlusComplete,0);  
  vector<double> isoPressure(nTypesPlusComplete,0); 
  vector<double> eqPressure(nTypesPlusComplete,0); 
  vector<double> bulkPressure(nTypesPlusComplete,0);
  vector<double> vx(nTypesPlusComplete,0);    
  vector<double> vy(nTypesPlusComplete,0); 
  vector<double> vz(nTypesPlusComplete,0); 
  vector<double> vv(nTypesPlusComplete,0);
  vector<double> gamma(nTypesPlusComplete,0);
  vector<double> PI00(nTypesPlusComplete,0);    
  vector<double> PI11(nTypesPlusComplete,0); 
  vector<double> PI22(nTypesPlusComplete,0); 
  vector<double> PI33(nTypesPlusComplete,0);   
  vector<double> PI10(nTypesPlusComplete,0);    
  vector<double> PI20(nTypesPlusComplete,0); 
  vector<double> PI30(nTypesPlusComplete,0); 
  vector<double> PI21(nTypesPlusComplete,0); 
  vector<double> PI31(nTypesPlusComplete,0);    
  vector<double> PI32(nTypesPlusComplete,0); 
  vector<double> PImunuPImunu(nTypesPlusComplete,0); 
  vector<double> q0(nTypesPlusComplete,0); 
  vector<double> q1(nTypesPlusComplete,0);
  vector<double> q2(nTypesPlusComplete,0);
  vector<double> q3(nTypesPlusComplete,0);
  vector<double> W0(nTypesPlusComplete,0); 
  vector<double> W1(nTypesPlusComplete,0);
  vector<double> W2(nTypesPlusComplete,0);
  vector<double> W3(nTypesPlusComplete,0);
  vector<double> V0(nTypesPlusComplete,0); 
  vector<double> V1(nTypesPlusComplete,0);
  vector<double> V2(nTypesPlusComplete,0);
  vector<double> V3(nTypesPlusComplete,0);  


  theHydroParticleType.giveHydroObservablesInLandauFrame( vecType, T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, N0, N1, N2, N3,
                                energyDensity, particleDensity, entropyDensity, temperature, fugacity, eq_particleDensity,
                                isoPressure, eqPressure, bulkPressure, vx, vy, vz, vv, gamma,
                                PI00, PI11, PI22, PI33, PI10, PI10, PI30, PI21, PI31, PI32, PImunuPImunu, q0, q1, q2, q3,
                                W0, W1, W2, W3, V0, V1, V2, V3);

/*  
  theHydroParticleType.giveHydroObservablesInEckartFrame( vecType, T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, N0, N1, N2, N3,
                                energyDensity, particleDensity, entropyDensity, temperature, fugacity, eq_particleDensity,
                                isoPressure, eqPressure, bulkPressure, vx, vy, vz, vv, gamma,
                                PI00, PI11, PI22, PI33, PI10, PI10, PI30, PI21, PI31, PI32, PImunuPImunu, q0, q1, q2, q3,
                                W0, W1, W2, W3, V0, V1, V2, V3);  
 */
  //-------------------//
  //for the massless case this is valid
  double T = energyDensity[nTypes] / (3.0 * particleDensity[nTypes]);
  scalingFactorLambdaT2 = fugacity[nTypes] * pow( T,2 );   
  
  //this is for the general case
//   scalingFactorLambdaT2 = fugacity[nTypes] * pow( temperature[nTypes],2 ); 


  if( FPT_COMP_LE(scalingFactorLambdaT2, 0.0) || NN[nTypes] <= 2)
    {scalingFactorLambdaT2 = 100000000000000.0;}//some arbitrary large value
  //-------------------//    

  return scalingFactorLambdaT2; 
}
