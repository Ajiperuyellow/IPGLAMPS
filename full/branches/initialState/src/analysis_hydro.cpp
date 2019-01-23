//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/analysis_hydro.cpp $
//$LastChangedDate: 2014-02-12 19:12:33 +0100 (水, 12  2月 2014) $
//$LastChangedRevision: 1619 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include "analysis_hydro.h"
#include "analysis.h"
#include "particle.h"

using namespace std;
using namespace ns_casc;



anaHydro::anaHydro( config * const _theConfig, const int _NNtimeStepsForAnaHydro):
  theConfig( _theConfig )
{
  //-------------------------------------------
  //---- Initialisation the hydro analysis ----
  //-------------------------------------------
  
  //----------------------------------------- 
  double columnWidthX, columnWidthY, columnWidthZ, numberColumnsX, numberColumnsY, numberColumnsZ;
  string name;
  
  realTimeSteps = _NNtimeStepsForAnaHydro;
  
  nTypes = theHydroParticleType.hydroParticleTypeVectorAnalysis.size();
  boxLengthX = theConfig->getBoxLengthForHydroX();
  boxLengthY = theConfig->getBoxLengthForHydroY();
  boxLengthZ = theConfig->getBoxLengthForHydroZ();
  //----------------------------------------
  
  //-----------------------------------------
  //normal profile  
  numberColumnsX = theConfig->getNumberColumnForHydroX();
  numberColumnsY = theConfig->getNumberColumnForHydroY();
  numberColumnsZ = theConfig->getNumberColumnForHydroZ();
  
  columnWidthX = boxLengthX / numberColumnsX;
  columnWidthY = boxLengthY / numberColumnsY;
  columnWidthZ = boxLengthZ / numberColumnsZ;
  
  name = "normal profile";
  
  anaHydroProfileNormal = new anaHydroData( nTypes,realTimeSteps,numberColumnsX,numberColumnsY,numberColumnsZ,columnWidthX,columnWidthY,columnWidthZ,name );
  //-----------------------------------------
  
  //-----------------------------------------
  //arrow profile
  numberColumnsX = theConfig->getNumberColumnArrowForHydroX();
  numberColumnsY = theConfig->getNumberColumnArrowForHydroY();
  numberColumnsZ = theConfig->getNumberColumnArrowForHydroZ();
  
  columnWidthX = boxLengthX / numberColumnsX;
  columnWidthY = boxLengthY / numberColumnsY;
  columnWidthZ = boxLengthZ / numberColumnsZ;

  name = "arrow profile";  
  
  anaHydroProfileArrow = new anaHydroData( nTypes,realTimeSteps,numberColumnsX,numberColumnsY,numberColumnsZ,columnWidthX,columnWidthY,columnWidthZ,name );  
  //-----------------------------------------  
  
  //-----------------------------------------
  //nromal profile at midrapidity
  numberColumnsX = theConfig->getNumberColumnMidRapForHydroX();
  numberColumnsY = theConfig->getNumberColumnMidRapForHydroY();
  numberColumnsZ = 1;
  
  columnWidthX = boxLengthX / numberColumnsX;
  columnWidthY = boxLengthY / numberColumnsY;
  columnWidthZ = -1.0;
  
  name = "normal midrapidity profile";  
  
  anaHydroProfileMidRapNormal = new anaHydroData( nTypes,realTimeSteps,numberColumnsX,numberColumnsY,numberColumnsZ,columnWidthX,columnWidthY,columnWidthZ,name );    
  //-----------------------------------------  
  
  //-----------------------------------------
  //arrow profile at midrapidity 
  numberColumnsX = theConfig->getNumberColumnMidRapArrowForHydroX();
  numberColumnsY = theConfig->getNumberColumnMidRapArrowForHydroY();
  numberColumnsZ = 1;
  
  columnWidthX = boxLengthX / numberColumnsX;
  columnWidthY = boxLengthY / numberColumnsY;
  columnWidthZ = -1.0;
  
  name = "arrow midrapidity profile";   
  
  anaHydroProfileMidRapArrow = new anaHydroData( nTypes,realTimeSteps,numberColumnsX,numberColumnsY,numberColumnsZ,columnWidthX,columnWidthY,columnWidthZ,name );    
  //-----------------------------------------
}     


anaHydro::~anaHydro()
{
  delete anaHydroProfileNormal;
  delete anaHydroProfileArrow;
  delete anaHydroProfileMidRapNormal;
  delete anaHydroProfileMidRapArrow;
}


anaHydroData::anaHydroData( const int _nTypes, const int _realTimeSteps, const int _numberColumnsX, const int _numberColumnsY,
                            const int _numberColumnsZ, const double _columnWidthX, const double _columnWidthY, const double _columnWidthZ,
                            const string _name)
{
  
  numberColumnsX = _numberColumnsX;
  numberColumnsY = _numberColumnsY;
  numberColumnsZ = _numberColumnsZ;
  columnWidthX = _columnWidthX;
  columnWidthY = _columnWidthY;
  columnWidthZ = _columnWidthZ;
  name = _name;
  
  //------------------------------------
  array_typeTmunu::extent_gen extentsTmunu; 
  array_typeNmu::extent_gen extentsNmu;
  array_typeMidRapOb::extent_gen extentsMidRapOb;
    
  enDistNmu.resize(extentsNmu[mu][_nTypes][_realTimeSteps][numberColumnsX][numberColumnsY][numberColumnsZ]);
  enDistTmunu.resize(extentsTmunu[mu][nu][_nTypes][_realTimeSteps][numberColumnsX][numberColumnsY][numberColumnsZ]);
  
  int nTypesComplete = _nTypes + 1;
  midRapOb_dN_dy.resize(extentsMidRapOb[nTypesComplete][_realTimeSteps]);
  midRapOb_dEt_dy.resize(extentsMidRapOb[nTypesComplete][_realTimeSteps]);
  midRapOb_v2.resize(extentsMidRapOb[nTypesComplete][_realTimeSteps]);
  midRapOb_v4.resize(extentsMidRapOb[nTypesComplete][_realTimeSteps]);
  
  //--------------------------------------------------------  
  for(int type = 0; type < _nTypes; type++){  
    for(int nn = 0; nn < _realTimeSteps; nn++){    
      for(int i = 0; i < numberColumnsX; i++){
        for(int j = 0; j < numberColumnsY; j++){
          for(int k = 0; k < numberColumnsZ; k++){
            for(int indexMu = 0; indexMu < mu; indexMu++){
              for(int indexNu = 0; indexNu < nu; indexNu++){      
                
                enDistTmunu[indexMu][indexNu][type][nn][i][j][k] = 0.0;
                
              }
            }
          }
        }
      }
    }
  }
  
  for(int type = 0; type < _nTypes; type++){  
    for(int nn = 0; nn < _realTimeSteps; nn++){    
      for(int i = 0; i < numberColumnsX; i++){
        for(int j = 0; j < numberColumnsY; j++){
          for(int k = 0; k < numberColumnsZ; k++){
            for(int indexMu = 0; indexMu < mu; indexMu++){   
              
              enDistNmu[indexMu][type][nn][i][j][k] = 0.0;
                
            }
          }
        }
      }
    }
  }
  
  for(int type = 0; type <= _nTypes; type++){  
    for(int nn = 0; nn < _realTimeSteps; nn++){
      midRapOb_dN_dy[type][nn] = 0.0;
      midRapOb_dEt_dy[type][nn] = 0.0;
      midRapOb_v2[type][nn] = 0.0;
      midRapOb_v4[type][nn] = 0.0;
    }
  }
  //--------------------------------------------
}



void anaHydro::hydroDistribution(const int nn, anaHydroData * ad)
{
  vector<double> P(mu,0.0);//momentum vector 

  int locationX,locationY,locationZ;
  int type;
  int testpartcl = theConfig->getTestparticles(); 
  
  double volumeBin = ad->columnWidthX * ad->columnWidthY * ad->columnWidthZ;
  if( volumeBin < 0.0){cout << "Error in hydroDistribution - negative volume!" << endl;}
  
  double vt = volumeBin * testpartcl;  
  
  vector<hydroParticleTypeProperties> vecTypeAna(theHydroParticleType.hydroParticleTypeVectorAnalysis);
  
  //   int count_dN_dSpaceRap = 0;//counting the number of particles in that region
  //   double count_dEt_dSpaceRap = 0.0;//counting the the transverse energy in that region    
 
  //---------------    
  for(unsigned int i = 0;i < particles.size(); i++)
  {      
    //check whether particles are in the unvisible box :-)
    if( fabs(particles[i].Pos.X()) < 0.5 * boxLengthX &&
        fabs(particles[i].Pos.Y()) < 0.5 * boxLengthY &&
        fabs(particles[i].Pos.Z()) < 0.5 * boxLengthZ  )
    {
      
      //searching for the location of the paricle
      locationX = ( particles[i].Pos.X() + 0.5 * boxLengthX ) / ad->columnWidthX;      
      locationY = ( particles[i].Pos.Y() + 0.5 * boxLengthY ) / ad->columnWidthY;      
      locationZ = ( particles[i].Pos.Z() + 0.5 * boxLengthZ ) / ad->columnWidthZ;

      //in case you have more particle species      
      type = theHydroParticleType.getParticleType(particles[i].FLAVOR,vecTypeAna);
        
      if(type >= 0)
      {    
        //cartesian coordinates
        P[0] = particles[i].Mom.E(); // P_E
        P[1] = particles[i].Mom.Px(); // P_X
        P[2] = particles[i].Mom.Py(); // P_Y
        P[3] = particles[i].Mom.Pz(); // P_Z
        
        //      //falls 2D
        //      count_dN_dSpaceRap++;
        //      count_dEt_dSpaceRap += sqrt(pow(P[1],2) + pow(P[2],2) + pow(mass,2));   
        
        //--------------- 
        for(int indexMu = 0; indexMu < mu; indexMu++){
            
          ad->enDistNmu[indexMu][type][nn][locationX][locationY][locationZ] += ( P[indexMu] / P[0] ) / vt;
             
          for(int indexNu = 0; indexNu < nu; indexNu++){   
            ad->enDistTmunu[indexMu][indexNu][type][nn][locationX][locationY][locationZ] += ( P[indexMu] * P[indexNu] / P[0] ) / vt;
          }
        }
        //---------------
      }
    }
  }
  //---------------
  
  //   //----------------------
  //   //this is only when having 2D
  //   showInfosAtMidRapidity(nn,volumeBin, ad->columnWidthZ,count_dN_dSpaceRap,count_dEt_dSpaceRap,ad);
  //   //----------------------   
}



void anaHydro::hydroDistributionMidRap(const int nn, const double time, const double timeshift, anaHydroData * ad)
{
  vector<double> P(mu,0.0);//momentum vector
  vector<int> count_N(nTypes+ 1,0.0);
  
  //for the case of midrapidity
  double spaceRapZRange = theConfig->getSpaceRapZRange();
  
  double tau = time;//in midrapidity tau = t_lab
  int testpartcl = theConfig->getTestparticles();
  double spaceRapidityOnlyArea = ad->columnWidthX * ad->columnWidthY * ( tanh(spaceRapZRange/2.0) - tanh(-spaceRapZRange/2.0) );
  double volumeBin = spaceRapidityOnlyArea * tau;
  double vt = volumeBin * testpartcl;
  
  vector<hydroParticleTypeProperties> vecTypeAna(theHydroParticleType.hydroParticleTypeVectorAnalysis); 

  //-----------------------------
  for(unsigned int i = 0;i < particles.size(); i++)
  {
    int locationX,locationY,locationZ;
    int type;  
    double rapZ, realTime;

    //check whether particles are in the unvisible box :-)
    if( fabs(particles[i].Pos.X()) < 0.5 * boxLengthX &&
        fabs(particles[i].Pos.Y()) < 0.5 * boxLengthY  )
    {            

      //searching for the location of the paricle       
      locationX = ( particles[i].Pos.X() + 0.5 * boxLengthX ) / ad->columnWidthX;      
      locationY = ( particles[i].Pos.Y() + 0.5 * boxLengthY ) / ad->columnWidthY;      

      realTime = particles[i].Pos.T() - timeshift;//could also be "time", but this is for error check!!!

      //to calculate if particle is in the space rapidity region
      rapZ = 0.5 * log( ( realTime + particles[i].Pos.Z() ) / ( realTime - particles[i].Pos.Z() ) );
      locationZ = 0;//this has to be zero

      //in case you have more particle species  
      type = theHydroParticleType.getParticleType(particles[i].FLAVOR,vecTypeAna);
      
      if(type >= 0)
      {
        if(fabs(rapZ) < spaceRapZRange / 2.0 )
        {     
          double mass, p_T;

          //cartesian coordinates
          P[0] = particles[i].Mom.E(); // P_E
          P[1] = particles[i].Mom.Px(); // P_X
          P[2] = particles[i].Mom.Py(); // P_Y
          P[3] = particles[i].Mom.Pz(); // P_Z
          mass = particles[i].m; //mass of particle

            p_T  = sqrt( pow(P[1],2) + pow(P[2],2) );

            //------------------
            count_N[type] ++;
            //------------------
            
            //observables averaged over the whole midrapidity region        
            ad->midRapOb_dN_dy[type][nn]++;
            ad->midRapOb_dEt_dy[type][nn]+= sqrt( pow(P[1],2) + pow(P[2],2) + pow(mass,2) );
            ad->midRapOb_v2[type][nn] += ( pow(P[1],2) - pow(P[2],2) ) / pow( p_T,2 );
            ad->midRapOb_v4[type][nn] += ( pow(P[1],4) - 6.0 * pow(P[1],2) * pow(P[2],2) + pow(P[2],4) ) / pow( p_T,4 );    
            //-----------------------------------------------    

          //--------------- 
          for(int indexMu = 0; indexMu < mu; indexMu++){
              
            ad->enDistNmu[indexMu][type][nn][locationX][locationY][locationZ] += ( P[indexMu] / P[0] ) / vt;
              
            for(int indexNu = 0; indexNu < nu; indexNu++){   
              ad->enDistTmunu[indexMu][indexNu][type][nn][locationX][locationY][locationZ] += ( P[indexMu] * P[indexNu] / P[0] ) / vt;
            }
          }
          //---------------
        }
      }
    }
  }

    //-----------------------------    
    //observables averaged over the whole midrapidity region
    
    //get complete type
    for(int type = 0; type < nTypes; type++){  
      count_N[nTypes] += count_N[type];
      ad->midRapOb_dN_dy[nTypes][nn] += ad->midRapOb_dN_dy[type][nn];
      ad->midRapOb_dEt_dy[nTypes][nn] += ad->midRapOb_dEt_dy[type][nn];
      ad->midRapOb_v2[nTypes][nn] += ad->midRapOb_v2[type][nn];
      ad->midRapOb_v4[nTypes][nn] += ad->midRapOb_v4[type][nn];
    }    
    
    //normalize correct
    for(int type = 0; type <= nTypes; type++){    
      ad->midRapOb_dN_dy[type][nn] /= double(testpartcl) * spaceRapZRange ;
      ad->midRapOb_dEt_dy[type][nn] /= double(testpartcl) * spaceRapZRange ;
      ad->midRapOb_v2[type][nn] /= double(count_N[type]);
      ad->midRapOb_v4[type][nn] /= double(count_N[type]);
      //-----------------------------------------------------
    }
    
  //----------------------
  showInfosAtMidRapidity(nn,volumeBin,spaceRapZRange,ad);
  //----------------------  
}


void anaHydro::showInfosAtMidRapidity(const int nn, const double volumeBin, const double rapidityRangeInZ, const anaHydroData * ad )
{
  vector<hydroParticleTypeProperties> vecTypeAna(theHydroParticleType.hydroParticleTypeVectorAnalysis); 
  int nTypesPlusComplete = nTypes + 1;   
  
  //values at zero position
  double T_central, fug_central, n_central, e_central, T00_central;
  
  //-----------------------------
  //-----------------------------
  //this calculates the entropy and other observables per rapidity
  typedef boost::multi_array<double, 4> array_typeTemp;
  array_typeTemp::extent_gen extentsTemp;
  
  array_typeTemp array_entropyDensity;
  array_typeTemp array_particleDensity;
  array_typeTemp array_gamma; 
  
  array_entropyDensity.resize(extentsTemp[nTypesPlusComplete][ad->numberColumnsX][ad->numberColumnsY][ad->numberColumnsZ]);
  array_particleDensity.resize(extentsTemp[nTypesPlusComplete][ad->numberColumnsX][ad->numberColumnsY][ad->numberColumnsZ]);
  array_gamma.resize(extentsTemp[nTypesPlusComplete][ad->numberColumnsX][ad->numberColumnsY][ad->numberColumnsZ]);  
  
  for(int type = 0; type <= nTypes; type++){    
    for(int i = 0; i < ad->numberColumnsX; i++){
      for(int j = 0; j < ad->numberColumnsY; j++){
        for(int k = 0; k < ad->numberColumnsZ; k++){
          array_entropyDensity[type][i][j][k] = 0.0;  
          array_particleDensity[type][i][j][k] = 0.0;
          array_gamma[type][i][j][k] = 0.0;
        }
      }
    }
  }
  //-----------------------------
 
  //-----------------------------
  for (int i = 0 ; i < ad->numberColumnsX ; i++)
  {
    for (int j = 0 ; j < ad->numberColumnsY ; j++)
    {
      for (int k = 0 ; k < ad->numberColumnsZ ; k++)
      {    
        vector<double> T00(nTypesPlusComplete,0);    
        vector<double> T11(nTypesPlusComplete,0); 
        vector<double> T22(nTypesPlusComplete,0); 
        vector<double> T33(nTypesPlusComplete,0);   
        vector<double> T10(nTypesPlusComplete,0);    
        vector<double> T20(nTypesPlusComplete,0); 
        vector<double> T30(nTypesPlusComplete,0); 
        vector<double> T21(nTypesPlusComplete,0); 
        vector<double> T31(nTypesPlusComplete,0);    
        vector<double> T32(nTypesPlusComplete,0);
        vector<double> N0(nTypesPlusComplete,0); 
        vector<double> N1(nTypesPlusComplete,0);
        vector<double> N2(nTypesPlusComplete,0);
        vector<double> N3(nTypesPlusComplete,0);      
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
  
        for(int type = 0; type < nTypes; type++){
          int completeType = nTypes;  
              
          T00[type] = ad->enDistTmunu[0][0][type][nn][i][j][k];
          T11[type] = ad->enDistTmunu[1][1][type][nn][i][j][k];
          T22[type] = ad->enDistTmunu[2][2][type][nn][i][j][k];
          T33[type] = ad->enDistTmunu[3][3][type][nn][i][j][k];

          T10[type] = ad->enDistTmunu[1][0][type][nn][i][j][k];
          T20[type] = ad->enDistTmunu[2][0][type][nn][i][j][k];
          T30[type] = ad->enDistTmunu[3][0][type][nn][i][j][k];
                  
          T21[type] = ad->enDistTmunu[2][1][type][nn][i][j][k];
          T31[type] = ad->enDistTmunu[3][1][type][nn][i][j][k];
          T32[type] = ad->enDistTmunu[3][2][type][nn][i][j][k];
                  
          N0[type] = ad->enDistNmu[0][type][nn][i][j][k];
          N1[type] = ad->enDistNmu[1][type][nn][i][j][k];
          N2[type] = ad->enDistNmu[2][type][nn][i][j][k];
          N3[type] = ad->enDistNmu[3][type][nn][i][j][k];
            
          T00[completeType] += T00[type];
          T11[completeType] += T11[type];
          T22[completeType] += T22[type];
          T33[completeType] += T33[type];

          T10[completeType] += T10[type];
          T20[completeType] += T20[type];
          T30[completeType] += T30[type];
                  
          T21[completeType] += T21[type];
          T31[completeType] += T31[type];
          T32[completeType] += T32[type];
                  
          N0[completeType] += N0[type];
          N1[completeType] += N1[type];
          N2[completeType] += N2[type];
          N3[completeType] += N3[type];
        }
            
        theHydroParticleType.giveHydroObservablesInLandauFrame( vecTypeAna, T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, N0, N1, N2, N3,
                                                                energyDensity, particleDensity, entropyDensity, temperature, fugacity, eq_particleDensity,
                                                                isoPressure, eqPressure, bulkPressure, vx, vy, vz, vv, gamma,
                                                                PI00, PI11, PI22, PI33, PI10, PI10, PI30, PI21, PI31, PI32, PImunuPImunu, q0, q1, q2, q3,
                                                                W0, W1, W2, W3, V0, V1, V2, V3);
            
            
            
        for(int type = 0; type <= nTypes; type++){    
          array_entropyDensity[type][i][j][k] = entropyDensity[type];  
          array_particleDensity[type][i][j][k] = particleDensity[type];
          array_gamma[type][i][j][k] = gamma[type];
        }
            
        if( i == ( (ad->numberColumnsX - 1) / 2 ) && j == ( (ad->numberColumnsY - 1) / 2) )
        {
          T_central = temperature[nTypes];
          fug_central = fugacity[nTypes];
          n_central = particleDensity[nTypes];
          e_central = energyDensity[nTypes];
          T00_central = T00[nTypes];
        } 
      }
    }
  }
  //-----------------------------      

  //-----------------------------
  double wholeEntropyAtMidRapidity = 0.0;
  double wholeParticlesAtMidRapidity = 0.0;
  double wholeT00AtMidRapidity = 0.0;
  double wholeN0AtMidRapidity = 0.0;  
    
  for(int type = 0; type < nTypes; type++){    
    for(int i = 0; i < ad->numberColumnsX; i++){
      for(int j = 0; j < ad->numberColumnsY; j++){
        for(int k = 0; k < ad->numberColumnsZ; k++){
          wholeEntropyAtMidRapidity += array_entropyDensity[type][i][j][k] * array_gamma[type][i][j][k] * volumeBin;
          wholeParticlesAtMidRapidity += array_particleDensity[type][i][j][k] * array_gamma[type][i][j][k] * volumeBin;
          wholeT00AtMidRapidity += ad->enDistTmunu[0][0][type][nn][i][j][k] * volumeBin;
          wholeN0AtMidRapidity += ad->enDistNmu[0][type][nn][i][j][k] * volumeBin;
        }
      }
    }
  }
  //-----------------------------    

  cout << "-----------------------------------------------------" << endl;
  cout << "  Values at midrapidity for the " << ad->name << endl;
  cout << "  dNN/dy        = " << ad->midRapOb_dN_dy[nTypes][nn] << endl;
  cout << "  dE_t/dy       = " << ad->midRapOb_dEt_dy[nTypes][nn]  << endl;   
  cout << "  dS/dy         = " << wholeEntropyAtMidRapidity / rapidityRangeInZ << endl;
  cout << "  d(n*gamma)/dy = " << wholeParticlesAtMidRapidity / rapidityRangeInZ << endl;  
  cout << "  dN0/dy        = " << wholeN0AtMidRapidity / rapidityRangeInZ << endl;      
  cout << "  dT00/dy       = " << wholeT00AtMidRapidity / rapidityRangeInZ << endl;      
  cout << "-----------------------------------------------------" << endl;
  
  cout << "-----------------------------------------------------" << endl;
  cout << "  Values at zero bin for the " << ad->name << endl;
  cout << "  Temperature : " << T_central << " GeV" << endl;
  cout << "  fugacity    : " << fug_central << endl;
  cout << "  nDensity    : " << n_central << " 1/fm^3 " << endl;
  cout << "  eDensity    : " << e_central << " GeV/fm^3" << endl;  
  cout << "  T00         : " << T00_central << " GeV/fm^3" << endl;   
  cout << "-----------------------------------------------------" << endl;  
  
}




int anaHydro::particleListgetFinalNumber()
{
  vector<hydroParticleTypeProperties> vecTypeAna(theHydroParticleType.hydroParticleTypeVectorAnalysis);   
  finalNumberInSimulation = 0;  
  
  for(unsigned int i=0; i< particles.size(); i++)
  {
    //check whether particles are in the unvisible box :-)
    if( fabs(particles[i].Pos.X()) < 0.5 * boxLengthX &&
        fabs(particles[i].Pos.Y()) < 0.5 * boxLengthY &&
        fabs(particles[i].Pos.Z()) < 0.5 * boxLengthZ  )
    { 
      //in case you have more particle species  
      int type = theHydroParticleType.getParticleType(particles[i].FLAVOR,vecTypeAna);
        
      if(type >= 0)
      {
        finalNumberInSimulation++;
      }
    }
  }  
    
  return finalNumberInSimulation;
}
