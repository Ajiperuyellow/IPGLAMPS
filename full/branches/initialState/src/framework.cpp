//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/framework.cpp $
//$LastChangedDate: 2019-01-05 18:02:58 +0100 (Sa, 05. Jan 2019) $
//$LastChangedRevision: 2917 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <list>
#include <map>
#include <utility>
#include <numeric>
#include <sstream>
#include <algorithm>

#include "coupling.h"
#include "particle.h"
#include "random.h"
#include "heavyIonCollision.h"
#include "FPT_compare.h"
#include "lorentz.h"
#include "configuration.h"
#include "ringstructure.h"
#include "cellcontainer.h"
#include "ratesmanager.h"
#include "binary_cross_sections.h"
#include "cellcontainer.h"
#include "unit_enum.h"
#include "scattering22_hydro.h"

#include "debug.h"

using std::cout;
using std::endl;
using std::vector;
using std::list;
using std::map;
using std::pair;
using namespace ns_casc;


namespace
{
  const int cellcut = 0;  //WARNING
  int IX, IY, IZ;
  double timenow, timenext;
  double dx, dy, dv, transLen, dt;

  vector<cellContainer> cells;
  vector< list<int> > edgeCell;
  list<int> formGeom;

  list<int> deadParticleList;
  coordinateBins dNdEta;
  coordinateEtaBins etaBins;

  long ncoll, ncoll32, ncoll23, ncoll22, ncolle, ncoll22MIDRAP;
  double sin2Theta, vrelSum, sin2ThetaVrel;
  double randomShiftEta, randomShiftX, randomShiftY;
  bool dodo1;
  int nn_ana1;
  
  int countMFP;
  double MFPRatio,MFPfm;
  
}


/**
 * This routine provides the framework for the actual simulation. Its central element is a do-while-loop that
 * loops over all time steps until the desired run time is reached and calls the scattering routines.
 *
 * @param[in] aa A pointer to the analysis object that is used for output, analysis and data collection. This pointer is passed on to various subroutines.
 */
void heavyIonCollision::mainFramework( analysis& aa )
{
  vector<Particle> particlesCopy;
  map< int, cellContainer > cellsBackup;
  vector< list<int> > edgeCellCopy;

  int NinCell, nn_ana2, NfrozenOut=0;
  int NinAct, NinActCell, NinFormInit, NinFormGeom, Nfree1, Nfree2;
  double md2g, md2q, sum, sum1, min, max, throughtime;
  bool again;


  const string sep = "\t";

  dodo1 = false;
  nn_ana1 = nn_ana2 = 0;

  
  //-------------------- build cell configuration --------------------  
  buildTransverseCellConfiguration(dx,dy,IX,IY,transLen);
  
  //Set average number of particles per cell, homogeneously averaged over transverse area, looking at total particles/transverse area
  NinCell = theConfig->getRequiredNinCell(); //usually 10-20
  cout  << "NinCell required: " << NinCell << endl;
  buildLongitudinalCellConfiguration(NinCell);
  
  normalizeCells();
  //------------------------------------------------------------------

  initializeAndResetParticles(aa);

  //---------- get the first timestep ----------
  timenow = 0.0;
  throughtime = 0.3;    // fm/c for p0=1.4, b=0,2,4,6 fm
  
  switch ( theConfig->getInitialStateType() )
  {
    case miniJetsInitialState:
      timenext = getFirstTimestep( dt, min );
      dt = dt * 0.001;
      break;
    case pythiaOnlineInitialState:
      timenext = getFirstTimestep( dt, min );
      dt = dt * 0.001;
      break;
    case pythiaOfflineInitialState:
      timenext = getFirstTimestep( dt, min );
      dt = dt * 0.001;
      break;
    case cgcInitialState:
      timenext = getFirstTimestep( dt, min );
      dt = dt * 0.001;
      break;
    case hydroParametrizationInitialState:
      timenext = timeshift + 0.01;
      dt = 0.001;
      break;
    case simpleInitialState:
      timenext = getFirstTimestep( dt, min );
      dt = dt * 0.001;
      break;
    case CYMInitialState:
      timenext = getFirstTimestepMinimum( dt, min );
      dt = dt * 0.001; //WARNING -> sensitive to "again"
      break;
    case AttractorInitialState:
      timenext = 0.4 ;//0.2005; //getFirstTimestepMinimum( dt, min );
      dt=0.01;
//       dt = dt * 0.001; //WARNING -> sensitive to "again"
      break;
    default:
      std::string errMsg = "Error in getting first timestep. No correct modus is chosen!";
      throw eHIC_error( errMsg );
      break;
  }
  cout << "Timestep = " << dt << endl;
  cout << "timenext = " << timenext << endl;
  

  //---------- get the first timestep ----------
  
  // Avoid analysis before first timestep because it leads to an earlier beginning of cascade 
  // by going back in time to the first analysis time step later on when the analysis routine is called.
  // Only for movie output (=6) analysis output is needed even at very early time steps.
  if( !theConfig->getOutputScheme() != doMovie )
    while(timenext >= aa.tstep[nn_ana1] + timeshift)
      nn_ana1++;
  //--------------------  

  DDD;

  infoCrossSectionMethod();    

  FixInitialScreeningMass(md2q,md2g);

  DDD;

  //------------------------------------------------------------------------------------------
  //---------------------------------- beginning of cascade ----------------------------------
  //------------------------------------------------------------------------------------------
  aa.generalOutput( initial, timeshift );

  DDD;
  
  ncoll = ncoll32 = ncoll23 = ncoll22 = ncolle = 0;
  sin2Theta = 0.;
  sin2ThetaVrel = 0.;
  vrelSum = 0.;
  ncoll22MIDRAP = 0;
  
  //---------- register data for offline reconstruction ----------
  if ( theConfig->doOutput_offlineReconstructionData() )
  { 
    DDD;
    offlineDataSimulationParameters collectSimulationParameters( aa.getSeed(), sqrtS, Bimp, A, Aatomic, B, Batomic,
                                                                 testpartcl, particles.size(), timenext, timeshift, 
                                                                 theConfig->getFreezeOutEnergyDensity(), rings.size(),
                                                                 rings.getCentralRadius(), rings.getDeltaR(), dx, dy,
                                                                 transLen, IX, IY, IZ, Particle::N_light_flavor, Particle::N_heavy_flavor );
    offlineInterface.submitOfflineDataForOutput( &collectSimulationParameters );

    offlineDataInitialParticles collectInitialParticleConfiguration( &particles );
    offlineInterface.submitOfflineDataForOutput( &collectInitialParticleConfiguration );
  }
  //---------- register data for offline reconstruction ----------
  
  DDD;

  if(theConfig->useRings())
  {
    scatterEdgeParticles( edgeCell );
  }else
  {
    //TODO at the moment only propagate
    scatterEdgeParticlesCellwise( edgeCell );
  }
  
  
  timenow = timenext; //fm
  int analysisNumber=0;
  DDD;

  cout << "=====================" << endl;
  cout << "=====================" << endl;
  cout << "Start parton cascade:" << endl;
  
 
  do
  {
    //HACK WORKS
//     for( int i = 0; i< particles.size();i++)
//     {
//       particles[i].Propagate( timenext );
//     }
  
    
  //HACK  
//     for ( int j = 0; j < particles.size(); j++ )
//     {
//       particles[j].Pos.T() =   timenext*0.999;   
//     }
    
    
    
    if ( timenow < throughtime )
    {
      randomShiftEta = 0;
    }
    else
    {
      randomShiftEta = ( etaBins[ etaBins.getCentralIndex()].right - etaBins[ etaBins.getCentralIndex()].left ) / 2  * ran2();
    }
    randomShiftX = dx * ran2();
    randomShiftY = dy * ran2();
    
   
    

    //---------- register data for offline reconstruction ----------   
    if ( theConfig->doOutput_offlineReconstructionData() )
    {
      offlineInterface.submitOfflineDataForOutput( event_newTimestep );
    }
    //---------- register data for offline reconstruction ----------
    
    //Generates the left-and right coordinates of the eta-bins.
    etaBins.populateEtaBins( particles, dNdEta, randomShiftEta, timenow, dt, dx, deltaEta_fine );
  
    // make copy ---------------
    particlesCopy = particles;
    edgeCellCopy = edgeCell;
    cellsBackup.clear();
    //--------------------------

    //--------------------------    
    timenext = timenow + dt;
    //--------------------------    
    //cout << timenow << "  time: " << timenext << " fm/c, timestep = " << dt << "  timeshift " << timeshift << "   aa.tstep[nn_ana1] " << aa.tstep[nn_ana1] << endl;
    
    if(theConfig->getAnalyseForHydro() == true)
    {cout << "time: " << timenext << " fm/c" << endl;}    
    
    if ( timenext >= (aa.tstep[nn_ana1] + timeshift) )
    {
      timenext = aa.tstep[nn_ana1] + timeshift;
      if(!FPT_COMP_E(timenow+dt,timenext))
      {
        dt = timenext - timenow; //AVOID ROUNDING ERROR
      }
      cout << dt << endl;
      dodo1 = true;
//       cout << "-------------------------------" << endl;
//       cout << "real time:          " << aa.tstep[nn_ana1] << endl;
//       cout << "time in simulation: " << timenext << endl;
      cout << "time: " << timenext << " ncoll22: " << ncoll22 << endl ;
    }

    
    
    cell_ID( NinAct, Nfree1, NinFormInit, NinFormGeom, cellsBackup, NfrozenOut );
    
//     cout << "CONTROL " << endl;
//     cout << NinAct << "\t" << NinActCell+Nfree2 << "\t" << Nfree1 << "\t" << NinFormInit << "\t" << NinFormGeom << "\t" << endl;
//     
    
//WARNING
//     aa.makeJetTrackerCopy();
    
    double ncoll22_backup = ncoll22;
    double ncoll23_backup = ncoll23;
    double ncoll32_backup = ncoll32;
    double ncolle_backup = ncolle;
    double sin2Theta_backup = sin2Theta;
    double sin2ThetaVrel_backup = sin2ThetaVrel;
    double vRel_backup = vrelSum;
    double ncoll22MIDRAP_backup = ncoll22MIDRAP;
    for( int ii = 0; ii < IZ; ii++ )
    {
      rateForOfflineOutput_gluons[ii].assign( rings.size(), 0 );
      rateForOfflineOutput_quarks[ii].assign( rings.size(), 0 );
      rateForOfflineOutput_antiQuarks[ii].assign( rings.size(), 0 );
    }
  
    deadParticleList.clear();
  
    

    
//     cout << "Original: " << particles[17].Pos << endl;
//     int i = 17;
//     double x=particles[i].Pos.X();
//     double y=particles[i].Pos.Y();
//     double z=particles[i].Pos.Z();  
//     double px=particles[i].Mom.Px();
//     double py=particles[i].Mom.Py();
//     double pz=particles[i].Mom.Pz(); 
//     double E=particles[i].Mom.E();
    
       scattering( NinActCell, Nfree2, again, aa, cellsBackup );
    
//    //TEST PROPAGATE ALL PARTICLES MYSELF
//       for(int i=0; i<particles.size();i++)
//       {
//         double x=particles[i].Pos.X();
//         double y=particles[i].Pos.Y();
//         double z=particles[i].Pos.Z();
//         
//         double px=particles[i].Mom.Px();
//         double py=particles[i].Mom.Py();
//         double pz=particles[i].Mom.Pz();
//         
//         double E=particles[i].Mom.E();
//         
//         
//         
//         particles[i].Pos=VectorTXYZ(timenext,x+px/E*dt,y+py/E*dt,z+pz/E*dt);
//         
//       }

//    cout << "After BAMPS  Propagation: " << particles[17].Pos << endl;
// 
//    double mydt = ( timenext - timenow );
//    VectorTXYZ particle17=VectorTXYZ(timenext,x+px/E*mydt,y+py/E*mydt,z+pz/E*mydt);
//    cout << "After my own Propagation: " << particle17 << endl;
//    //particles[i].Pos=VectorTXYZ(timenext,x+px/E*dt,y+py/E*dt,z+pz/E*dt);
    
    
    
    
    
    
    
    
    
    if ( theConfig->repeatTimesteps() )
    {
      while ( again )
      {
        dodo1 = false;
        
        aa.restoreJetTracker();
        
        particles = particlesCopy;
        edgeCell = edgeCellCopy;
        restoreCellBackup( cellsBackup );
        
        deadParticleList.clear();
        
        cout << "again" << "\t" << timenow << endl;
        
        timenext = timenow + dt;
     
        cell_ID( NinAct, Nfree1, NinFormInit, NinFormGeom, cellsBackup, NfrozenOut );
        
        vrelSum = vRel_backup;
        sin2Theta = sin2Theta_backup;
        sin2ThetaVrel = sin2ThetaVrel_backup;
        ncoll22MIDRAP = ncoll22MIDRAP_backup;
        ncoll22 = ncoll22_backup;
        ncoll23 = ncoll23_backup;
        ncoll32 = ncoll32_backup;
        ncolle = ncolle_backup;
        ncoll = ncoll22 + ncoll23 + ncoll32 + ncolle;
        
        for( int ii = 0; ii < IZ; ii++ )
        {
          rateForOfflineOutput_gluons[ii].assign( rings.size(), 0 );
          rateForOfflineOutput_quarks[ii].assign( rings.size(), 0 );
          rateForOfflineOutput_antiQuarks[ii].assign( rings.size(), 0 );
        }
        
        scattering( NinActCell, Nfree2, again, aa, cellsBackup );
      }
    }

    //---------- register data for offline reconstruction ----------
    if ( theConfig->doOutput_offlineReconstructionData() )
    {
      offlineDataCellConfiguration collectCellConfiguration( timenow, timenext, randomShiftX, randomShiftY, randomShiftEta, etaBins );
      offlineDataInteractionRates collectInteractionRates( rateForOfflineOutput_gluons, rateForOfflineOutput_quarks, rateForOfflineOutput_antiQuarks );
                                                          
      offlineInterface.submitOfflineDataForOutput( &collectCellConfiguration );
      offlineInterface.submitOfflineDataForOutput( &collectInteractionRates );
    }
    //---------- register data for offline reconstruction ----------


    // SCATTER EDGES //
    if(theConfig->useRings())
    {
      scatterEdgeParticles( edgeCell );
    }else
    {
      //TODO at the moment only propagate
       scatterEdgeParticlesCellwise( edgeCell );
    }
    
    removeDeadParticles( aa );

    if ( theConfig->doOutput_progressLog() )
    {
      std::vector<int> edgeCellSizes( edgeCell.size() );
      for ( unsigned int i = 0; i < edgeCell.size(); i++ )
      {
        edgeCellSizes[i] = edgeCell[i].size();
      } 
      aa.registerProgressInfoForOutput( timenow, dt, particles.size(), NinFormInit, NinFormGeom,
                                        NinAct, NinActCell, edgeCellSizes, Nfree1+Nfree2,
                                        ncoll, ncoll22, ncoll23, ncoll32, ncolle );
    }
          

    // // just error checking if masses and flavors are correct
    // for ( int j = 0; j < particles.size(); j++ )
    // {
    //   double E_check = sqrt( particles[j].Mom.vec2() + pow( particles[j].m, 2 ) );
    //   if( !FPT_COMP_E( particles[j].Mom.E(), E_check ) )
    //   {
    // 	cout << "Error! Particle " << j << " does not fulfill E^2 = p^2 + m^2. " << timenow << endl;
    // 	cout << " " << particles[j].FLAVOR << " " << particles[j].m << particles[j].Mom << particles[j].Pos << endl;
    // 	throw 123;
    //   }

    //   if( !FPT_COMP_E( particles[j].m, Particle::getMass( particles[j].FLAVOR ) ) )
    //   {
    // 	cout << "Error! Particle " << j << " with flavor " << particles[j].FLAVOR << " has wrong mass: " << particles[j].m << "  " << Particle::getMass( particles[j].FLAVOR )  << endl;
    // 	throw 123;
    //   }
      
    //   if( ( particles[j].FLAVOR > 2 * Particle::N_light_flavor ) && ( ( particles[j].FLAVOR <= 2 * 3 ) || ( particles[j].FLAVOR > 2 * ( 3 + Particle::N_heavy_flavor ) ) ) )
    //   {
    // 	cout << "Error! Particle " << j << " with flavor " << particles[j].FLAVOR << " should not exist: Nf = " << Particle::N_light_flavor << " + " <<   Particle::N_heavy_flavor <<  endl;
    // 	throw 123;
    //   }
    // }

      //TEST PROPAGATE ALL PARTICLES MYSELF
//       for(int i=0; i<particles.size();i++)
//       {
//         double x=particles[i].Pos.X();
//         double y=particles[i].Pos.Y();
//         double z=particles[i].Pos.Z();
//         
//         double px=particles[i].Mom.Px();
//         double py=particles[i].Mom.Py();
//         double pz=particles[i].Mom.Pz();
//         
//         double E=particles[i].Mom.E();
//         
//         dt = .01;
//         
//         particles[i].Pos=VectorTXYZ(timenext,x+px/E*dt,y+py/E*dt,z+pz/E*dt);
//         
//       }


    timenow = timenext;
   
    if ( dodo1 )
    {
      //TEST for AttractorInitialState
      //getCentralValues(timenow);
      if(theConfig->getOutputScheme()==onlyVideo)
      {
        int centralEtaIndex = static_cast<int>( etaBins.size() ) / 2;
        aa.PlotCellwise(timenow, analysisNumber,  IX,  IY, centralEtaIndex,cells);
        analysisNumber++;
      }
      aa.generalOutput( nn_ana1, timeshift );
      
      std::vector<int> edgeCellSizes( edgeCell.size() );
      for ( unsigned int i = 0; i < edgeCell.size(); i++ )
      {
        edgeCellSizes[i] = edgeCell[i].size();
      } 
      
//       //NEW
//       //Find Particles at Midrapidity
//       int particleSizeMidrap = 0;
//       for ( int j = 0; j < particles.size(); j++ )
//       {
//         if( fabs( particles[j].Mom.Rapidity()  ) >= -0.5 && fabs( particles[j].Mom.Rapidity()  ) <= 0.5 )
//         { 
//           particleSizeMidrap++;
//         }       
//       }
//       
      //TEST
//       aa.Cluster(nn_ana1);
//       

      
//       aa.Cluster(nn_ana1);
      //aa.writeParticleAndCollisionNumbers( timenow, particleSizeMidrap, ncoll22, ncoll23, ncoll32, nGet32Errors, nGet23Errors, NinFormInit, NinFormGeom, NinAct, NinActCell, edgeCellSizes, Nfree1+Nfree2, sin2Theta, vrelSum, sin2ThetaVrel, ncoll22MIDRAP );
      
      //TEST
      
      cout << NfrozenOut << "\t\t" << double(NfrozenOut)/particles.size() << endl;
      
      nn_ana1++;
      dodo1 = false;
    }
    
    
    
  }
  while ( timenow < stop + timeshift && (  double(NfrozenOut)/particles.size()  < 0.9  )  );

  // FINAL DCA CLUSTERING //
  cout << "Timestep - " << nn_ana1 << endl;
  makeSomeFakeQuarkPairs();
  aa.Cluster();
  aa.runHERWIG();
  aa.readHadrons();
  aa.ChargedHadronCorrelations(timenow);

  
  
  
  //---------- register data for offline reconstruction ----------
  if ( theConfig->doOutput_offlineReconstructionData() )
  {
    offlineInterface.submitOfflineDataForOutput( event_endOfCascade );
  }
  //---------- register data for offline reconstruction ----------
  aa.addJetEvents_final();
  aa.generalOutput( final, timeshift );
  cout << "number of errors in get32(...) = " << nGet32Errors << endl;
  cout << "number of errors in get23(...) = " << nGet23Errors << endl;

  //show and save the number of collisions
  cout << endl;  
  cout << "--------------------------------" << endl;
  cout << "Number of collisions in simulation: " << endl;
  cout << "ncoll22 " << ncoll22 << endl;
  cout << "ncoll23 " << ncoll23 << endl;
  cout << "ncoll32 " << ncoll32 << endl;
  cout << "Sum sin^2(theta) " << sin2Theta << endl;
  cout << "Sum sin^2(theta)*v_rel " << sin2ThetaVrel << endl;
  cout << "vrelSum " << vrelSum << endl;
  cout << "ncoll22MIDRAP " << ncoll22MIDRAP << endl;
  cout << "ncoll22/Nparticles " << ncoll22/particles.size() << endl;
  
  //get the number of collisions for the analysis
  aa.ncoll22 = ncoll22;
  aa.ncoll23 = ncoll23;
  aa.ncoll32 = ncoll32;  
  aa.sin2theta = sin2Theta; 
  aa.sin2ThetaVrel = sin2ThetaVrel;
  aa.vrelSum = vrelSum;
  aa.ncoll22MIDRAP = ncoll22MIDRAP;
  cout << "--------------------------------" << endl;
  cout << endl;
  

  
  //------------------------------------------------------------------------------------------
  //------------------------------------- end of cascade -------------------------------------
  //------------------------------------------------------------------------------------------
  
  // List particle numbers for all flavors
  listParticleNumbersForAllFlavors(finalState,0.0);


  for ( unsigned int i = 0; i < cells.size(); i++ )
  {
    cells[i].clear();
  }
  cellsBackup.clear();
  edgeCell.clear();
  edgeCellCopy.clear();
}

void heavyIonCollision::makeSomeFakeQuarkPairs()
{ 
  int indexOne,indexTwo;
  int fractionToTransform = floor((particles.size()/2.)*0.1);
  cout << "Transform " << 2*fractionToTransform << " of " << particles.size() << " particles to q and anti-qs" << endl;
  for(int i=0; i<fractionToTransform; i++)
  {
    for(int j=0;j<particles.size()-1;j++)
    {
      if(particles[j].FLAVOR!=gluon)continue;
      indexOne=j;
      for(int k=j+1;j<particles.size();k++)
      {
        if(particles[k].FLAVOR!=gluon)continue;
        indexTwo=k;
        break;
      }
      break;
    }
    
    
//     do
//     {
//       indexOne = floor(particles.size()*ran2());
//       indexTwo = floor(particles.size()*ran2());
//     }while(particles[indexOne].FLAVOR==gluon || particles[indexTwo].FLAVOR == gluon);
    
    int aflavor=floor(3*ran2());
    
    switch(aflavor)
    {
      case 0: {particles[indexOne].FLAVOR=up; particles[indexTwo].FLAVOR=anti_up;}break;
      case 1: {particles[indexOne].FLAVOR=down; particles[indexTwo].FLAVOR=anti_down;}break;
      case 2: {particles[indexOne].FLAVOR=strange; particles[indexTwo].FLAVOR=anti_strange;}break;
    }
  }
}

void heavyIonCollision::FixInitialScreeningMass(double & md2q,double & md2g)
{
  double sum,sum1,max;
  //---------- initial screening mass ----------
  sum = sum1 = max = 0.0;
  for ( unsigned int j = 0; j < particles.size(); j++ )
  {
    if ( particles[j].Pos.T() < timenext && particles[j].FLAVOR <= 2 * Particle::N_light_flavor ) // Heavy quarks do not contribute to screening mass, since routines are written for massless particles and the number of heavy quarks is negligible anyhow.
    {
      if ( particles[j].FLAVOR == 0 )
      {
        sum += 1.0 / particles[j].Mom.E();
      }
      else
      {
        sum1 += 1.0 / particles[j].Mom.E();
      }

      if ( fabs( particles[j].Pos.Z() ) > max )
      {
        max = fabs( particles[j].Pos.Z() );
      }
    }
  }
  sum = sum / gG / testpartcl;
  
  if ( Particle::N_light_flavor > 0 )
  {
    sum1 = sum1 / gQ / ( 2.0 * Particle::N_light_flavor ) / testpartcl;
  }
  else
  {
    sum1 = 0;
  }
  
  dv = M_PI * 36.0 * 2.0 * max;
  md2q = pow( 0.197, 3.0 ) * 2.0 * M_PI / dv;//GeV^3
  md2g = 8.0 * md2q * ( Ncolor * sum + Particle::N_light_flavor * sum1 );//GeV^2
  md2q = md2q * gG / gQ * ( sum + sum1 );//GeV^2

  DDD;

  for ( unsigned int k = 0; k < particles.size(); k++ )
  {
    if ( particles[k].Pos.T() < timenext )
    {
      particles[k].md2g = md2g;
      particles[k].md2q = md2q;

      addParticleToEdgeCell( k, edgeCell );
      particles[k].init = false;
    }
  }
  //---------- initial screening mass ----------
}

void heavyIonCollision::buildLongitudinalCellConfiguration(int NinCell)
{
  //Choose how to construct the eta bins
  if(theConfig->getInitialStateType()==CYMInitialState)
  {
    IZ = etaBins.constructEtaBinsPPb( NinCell, transLen, dx, dy, typicalRadius, theConfig->getTestparticles(), particles.size() );
  }else if (theConfig->getInitialStateType()==AttractorInitialState)
  {
    IZ = etaBins.constructEtaBinsPPb( NinCell, transLen, dx, dy, typicalRadius, theConfig->getTestparticles(), particles.size() );
  }
  else
  {
    IZ = etaBins.constructEtaBins( NinCell, Bimp, dx, dy, typicalRadius, theConfig->getTestparticles(), particles.size() );
  }

  etaBins.setTimestepScaling( theConfig->getScaleTimesteps() );
}

void heavyIonCollision::normalizeCells()
{
  vector<double> tempVec( rings.size(), 0 );
  rateForOfflineOutput_gluons.assign( IZ, tempVec );
  rateForOfflineOutput_quarks.assign( IZ, tempVec );
  rateForOfflineOutput_antiQuarks.assign( IZ, tempVec );

  edgeCell.resize( numberOfEdgeCells );  // handle 8 different edge cells
  for ( int i = 0; i < numberOfEdgeCells; i++ )
  {
    edgeCell[i].clear();
  }
  
  int ncell = IX * IY * IZ;
  cells.resize( ncell );
  for ( unsigned int i = 0; i < cells.size(); i++ )
  {
    cells[i].setCoordinates( i, dx, IX, transLen, dy, IY, transLen );
    cells[i].rates.normalizeRates();
  }

  cout << "IX=" << IX << "\t" << "IY=" << IY << "\t" << "IZ=" << IZ << endl;
  cout << "number of cells=" << ncell << endl;
  cout << "critical energy density: e_crit = " << theConfig->getFreezeOutEnergyDensity() << endl;  
}

void heavyIonCollision::initializeAndResetParticles(analysis & aa)
{
  //---------- initialize some variables ----------
  double pt_initial;
  for ( unsigned int j = 0; j < particles.size(); j++ )
  {
    particles[j].cell_id = -1;//-1:unformed
    particles[j].edge = -1;
    particles[j].dead = false;
    particles[j].init = true;
    particles[j].free = false;
    particles[j].coll_id = 0;
    particles[j].collisionTime = infinity;//fm
    particles[j].collisionPartner = -1;
    particles[j].rate23 = particles[j].rate32 = particles[j].rate22 = 0.0;//GeV
    particles[j].rate23v = particles[j].rate32v = particles[j].rate22v = 0.0;//GeV
    particles[j].as22 = particles[j].as23 = 0.0;
    particles[j].cs22 = particles[j].cs23;//1/GeV^2
    particles[j].md2g_scaled_22 = 0;
    particles[j].md2q_scaled_22 = 0;
    particles[j].md2g_scaled_23 = 0;
    particles[j].md2q_scaled_23 = 0;
    particles[j].lambda_scaled = 0;

    pt_initial = particles[j].Mom.Pt();
    if ( pt_initial > aa.getJetTracking_PT() )
    {
      aa.addJetEvent_initial( j );
    }
  }  
}

void heavyIonCollision::findFreeParticles(bool & free, int cellindex)
{
  free = true;
  list<int>::const_iterator iIt;
  for ( iIt = cells[cellindex].particleList.begin(); iIt != cells[cellindex].particleList.end(); iIt++ )
  {
    int id = *iIt;
    double xt = particles[id].Pos.Perp();

    int nc = rings.getIndexPure( xt );
  
    if ( nc < rings.size() )
    {
      if ( rings[nc].getEnergyDensity() < theConfig->getFreezeOutEnergyDensity() )
      {
        particles[id].free = true;
      }
      else
      {
        particles[id].free = false;
      }
    }
    else
    {
      particles[id].free = true;
    }

    if ( !particles[id].free )
    {
      free = false;
    }
  }
}

void heavyIonCollision::findFreeParticlesCellwise(bool & free, int cellindex)
{
  free = true;
  list<int>::const_iterator iIt;
  for ( iIt = cells[cellindex].particleList.begin(); iIt != cells[cellindex].particleList.end(); iIt++ )
  {
    int id = *iIt;
  
    if ( cells[cellindex].getEnergyDensity() < theConfig->getFreezeOutEnergyDensity() )
    {
      particles[id].free = true;
    }
    else
    {
      particles[id].free = false;
    }

    if ( !particles[id].free )
    {
      free = false;
    }
  }
}

void heavyIonCollision::updateParticleRatesRings(int cellindex)
{
  list<int>::const_iterator iIt;
  for ( iIt = cells[cellindex].particleList.begin(); iIt != cells[cellindex].particleList.end(); iIt++ )
  {
    int id = *iIt;

    double xt = particles[id].Pos.Perp();
    int nc = rings.getIndex( xt );

    particles[id].rate23v = particles[id].rate23;   //GeV
    particles[id].rate32v = particles[id].rate32;   //GeV
    particles[id].rate22v = particles[id].rate22;   //GeV

    particles[id].rate22 = cells[cellindex].rates.getRate( particles[id].FLAVOR, c22, GeV ) * rings[nc].getGamma();
    particles[id].rate23 = cells[cellindex].rates.getRate( particles[id].FLAVOR, c23, GeV ) * rings[nc].getGamma();
    particles[id].rate32 = cells[cellindex].rates.getRate( particles[id].FLAVOR, c32, GeV ) * rings[nc].getGamma();
    
    cells[cellindex].writeAveragesToParticle( particles[id] );
  }
}

void heavyIonCollision::updateParticleRatesCellwise(int cellindex)
{
  list<int>::const_iterator iIt;
  for ( iIt = cells[cellindex].particleList.begin(); iIt != cells[cellindex].particleList.end(); iIt++ )
  {
    int id = *iIt;

    double xt = particles[id].Pos.Perp();
    int nc = rings.getIndex( xt );

    particles[id].rate23v = particles[id].rate23;   //GeV
    particles[id].rate32v = particles[id].rate32;   //GeV
    particles[id].rate22v = particles[id].rate22;   //GeV

    particles[id].rate22 = cells[cellindex].rates.getRate( particles[id].FLAVOR, c22, GeV ) * cells[cellindex].getGamma();
    particles[id].rate23 = cells[cellindex].rates.getRate( particles[id].FLAVOR, c23, GeV ) * cells[cellindex].getGamma();
    particles[id].rate32 = cells[cellindex].rates.getRate( particles[id].FLAVOR, c32, GeV ) * cells[cellindex].getGamma();
    
    cells[cellindex].writeAveragesToParticle( particles[id] );
  }
}

void heavyIonCollision::handleCellcutCellsRingwise(int cellindex)
{
  list<int>::const_iterator iIt;
  for ( iIt = cells[cellindex].particleList.begin(); iIt != cells[cellindex].particleList.end(); iIt++ )
  {
    int id = *iIt;

    double xt = particles[id].Pos.Perp();
    int nc = rings.getIndex( xt );
    
    particles[id].rate22v = rings[nc].rates.getRate( particles[id].FLAVOR, c22, GeV );
    particles[id].rate23v = rings[nc].rates.getRate( particles[id].FLAVOR, c23, GeV );
    particles[id].rate32v = rings[nc].rates.getRate( particles[id].FLAVOR, c32, GeV );

    particles[id].rate23 = 0;        //GeV
    particles[id].rate32 = 0;        //GeV
    particles[id].rate22 = 0;        //GeV
    particles[id].as22 = 0;
    particles[id].as23 = 0;
    particles[id].cs22 = 0;          //1/GeV^2
    particles[id].cs23 = 0;          //1/GeV^2
    particles[id].md2g_scaled_22 = 0;
    particles[id].md2q_scaled_22 = 0;
    particles[id].md2g_scaled_23 = 0;
    particles[id].md2q_scaled_23 = 0;
    particles[id].lambda_scaled = 0;

    if ( particles[id].edge < 0 )      // new member in edge cell
    {
      addParticleToEdgeCell( id, edgeCell );
      
      particles[id].md2g = rings[nc].getAveraged_md2g();
      particles[id].md2q = rings[nc].getAveraged_md2q();
    }
    else
    {
      int newEdgeCellIndex = getEdgeCellIndex( id );
      if ( particles[id].edge != newEdgeCellIndex )  // check if the edge particle needs to be put into another edge cell
      {
        swapParticleToNewEdgeCell( id, particles[id].edge, newEdgeCellIndex, edgeCell );
      }
    }
  }
  cells[cellindex].rates.normalizeRates();  
}

void heavyIonCollision::handleCellcutCellsCellwise(int cellindex)
{
  list<int>::const_iterator iIt;
  for ( iIt = cells[cellindex].particleList.begin(); iIt != cells[cellindex].particleList.end(); iIt++ )
  {
    int id = *iIt;
    particles[id].rate22v = cells[cellindex].ratesSave.getRate( particles[id].FLAVOR, c22, GeV );
    particles[id].rate23v = cells[cellindex].ratesSave.getRate( particles[id].FLAVOR, c23, GeV );
    particles[id].rate32v = cells[cellindex].ratesSave.getRate( particles[id].FLAVOR, c32, GeV );

    particles[id].rate23 = 0;        //GeV
    particles[id].rate32 = 0;        //GeV
    particles[id].rate22 = 0;        //GeV
    particles[id].as22 = 0;
    particles[id].as23 = 0;
    particles[id].cs22 = 0;          //1/GeV^2
    particles[id].cs23 = 0;          //1/GeV^2
    particles[id].md2g_scaled_22 = 0;
    particles[id].md2q_scaled_22 = 0;
    particles[id].md2g_scaled_23 = 0;
    particles[id].md2q_scaled_23 = 0;
    particles[id].lambda_scaled = 0;

    particles[id].Propagate( timenext ); //Propagate particle
    
    if ( particles[id].edge < 0 )      // new member in edge cell
    {
      addParticleToEdgeCell( id, edgeCell );
      
      particles[id].md2g = cells[cellindex].getAveraged_md2g();
      particles[id].md2q = cells[cellindex].getAveraged_md2q();
    }
    else
    {
      int newEdgeCellIndex = getEdgeCellIndex( id );
      if ( particles[id].edge != newEdgeCellIndex )  // check if the edge particle needs to be put into another edge cell
      {
        swapParticleToNewEdgeCell( id, particles[id].edge, newEdgeCellIndex, edgeCell );
      }
    }
  }
  cells[cellindex].rates.normalizeRates();  
}

/**
 * This routine goes through all edge cells and performs the geometric scatterings of particles in each of these edge cells.
 * @param[in] _particleList A vector of edge cells. Each edge cell is a list of particle IDs contained in this edge cell.
 */
void heavyIonCollision::scatterEdgeParticles( vector< list<int> >& _particleList )
{
  int index;
  double xt;
  int nc;
  
  DDD;

  for ( unsigned int _edgeCellID = 0 ; _edgeCellID < _particleList.size(); _edgeCellID++ )
  {
    list<int>::iterator iIt;

    DDD;
    
    for ( iIt = _particleList[_edgeCellID].begin(); iIt != _particleList[_edgeCellID].end(); )
    {
      index = *iIt;
      
      if ( particles[index].md2g < 0 )  //new member
      {
        xt = particles[index].Pos.Perp();
        nc = rings.getIndex( xt );

        particles[index].md2g = rings[nc].getAveraged_md2g();
        // if new member is added to ring that has no particles in it, the "averaged" md2g is just zero, which leads to problems. Thus give it md2g from next ring that has particles in it.
        if( FPT_COMP_E( particles[index].md2g , 0.0 ) )
        {
          do
          {
            nc--;
            if( nc < 0 ) // was already smallest ring
              particles[index].md2g = ( Ncolor + Particle::N_light_flavor ) * 8 / M_PI * pow( 0.3 , 2.0 ); //GeV^2, Debye mass in medium with T=300MeV
            else
              particles[index].md2g = rings[nc].getAveraged_md2g(); // Debye mass from next smallest ring
          } while( FPT_COMP_E( particles[index].md2g , 0.0 ) );
        }
        
        nc = rings.getIndex( xt );

        particles[index].md2q = rings[nc].getAveraged_md2q();
        // if new member is added to ring that has no particles in it, the "averaged" md2q is just zero, which leads to problems. Thus give it md2q from next ring that has particles in it.
        if( FPT_COMP_E( particles[index].md2q , 0.0 ) )
        {
          do
          {
            nc--;
            if( nc < 0 ) // was already smallest ring
              particles[index].md2q = 16 / ( 3 * M_PI ) * pow( 0.3 , 2.0 ); //GeV^2, quark Debye mass in medium with T=300MeV
            else
              particles[index].md2q = rings[nc].getAveraged_md2q(); // Debye mass from next smallest ring
          } while( FPT_COMP_E( particles[index].md2q , 0.0 ) );
        }
      }
      
      if ( particles[index].edge != _edgeCellID )
      {
        iIt = _particleList[_edgeCellID].erase( iIt );
      }
      else
      {
        ++iIt;
      }
    }

    DDD;
    
    if(theConfig->geometricCollisions())
    {
      DDD;
      if ( _particleList[_edgeCellID].size() > 1 )
      {
        int particleIndexOfFirstCollision = updateGeometricCollisionTimesA( _particleList[_edgeCellID] );
        DDD;
        doGeometricCollisions( _particleList[_edgeCellID], particleIndexOfFirstCollision );
      }
    }

    DDD;
    
    for ( iIt = _particleList[_edgeCellID].begin(); iIt != _particleList[_edgeCellID].end(); iIt++ )
    {
      index = *iIt;
      if ( particles[index].Pos.T() < timenext )
      {
        // t = timenext, pos = pos + c*mom
        particles[index].Propagate( timenext );
      }
    }

    DDD;
    
  }
}

/**
 * This routine goes through all edge cells and performs the geometric scatterings of particles in each of these edge cells.
 * @param[in] _particleList A vector of edge cells. Each edge cell is a list of particle IDs contained in this edge cell.
 */
void heavyIonCollision::scatterEdgeParticlesCellwise( vector< list<int> >& _particleList )
{
  int index;
  double xt;
  int nc;
  
  for ( unsigned int _edgeCellID = 0 ; _edgeCellID < _particleList.size(); _edgeCellID++ )
  {
    list<int>::iterator iIt;
    for ( iIt = _particleList[_edgeCellID].begin(); iIt != _particleList[_edgeCellID].end(); iIt++ )
    {
      index = *iIt;
      if ( particles[index].Pos.T() < timenext )
      {
        // t = timenext, pos = pos + c*mom
        particles[index].Propagate( timenext );
      }
    }
  }
}


/**
 * This routine goes through all etaBins and prints its particle list. 
 */
void heavyIonCollision::controlParticlesInCell()
{
  int particleNumber = 0;
  int IXY = IX * IY;
  double cellCenter_x=0;
  double cellCenter_y=0;
  for ( int i = etaBins.min_index(); i <= etaBins.max_index(); i++ )
  {
    if (i==etaBins.getCentralIndex())
    {
      for ( int j = IXY * i; j < IXY * ( i + 1 ); j++ )
      {
        cellCenter_x=(cells[j].corner.x_max-cells[j].corner.x_min)/2.+cells[j].corner.x_min;
        cellCenter_y=(cells[j].corner.y_max-cells[j].corner.y_min)/2.+cells[j].corner.y_min;
        
        if(abs(cellCenter_x) < 1 && abs(cellCenter_y) < 1)
        cout << cellCenter_x << "\t" << cellCenter_y << "\t" <<    cells[j].particleList.size()<<endl;
        /*list<int>::const_iterator iIt;
        for ( iIt = cells[j].particleList.begin(); iIt != cells[j].particleList.end(); iIt++ )
        {
        }  
        */
      }
    }
  }
}



/**
 * @param[out] NinAct The number of active particles, i.e. particles for which T < timenext
 * @param[out] Nfree1 The number of particles that are free because they are outside the transversal grid
 * @param[out] NinFormInit The number of particles that are not active (T > timenext) because they are still within the initial formation time
 * @param[out] NinFormGeom The number of particles that are not active (T > timenext) because their time has been set within geometrical collision routines
 * @param[out] NfrozenOut The number of particles that are not active (T > timenext) because their time has been set within geometrical collision routines
 * @param[in,out] _cellsBackup The backup structure for cells whose content has changed within the current time step
 */
/*
Cell-ID:
-1: unformed
-2: edge cell, if particle position is outside the box
-100: free particle, does not scatter anymore
Otherwise:
cell_id = nx + IX * ny + IX * IY * nz;
Only particles which are neiter FREE nor EDGE are in the list:
cells[cell_id].particleList


*/

void heavyIonCollision::cell_ID( int& NinAct, int& Nfree1, int& NinFormInit, int& NinFormGeom, map< int, cellContainer >& _cellsBackup, int& NfrozenOut )
{
  int nx, ny, nz, cell_id;
  double halfsize;

  halfsize = transLen / 2.0;

  NfrozenOut=0;
  
  //OV NinAct=NinForm1=0;
  NinAct = Nfree1 = NinFormInit = NinFormGeom = 0;

  formGeom.clear();

  vector<Particle>::iterator it;
  for ( it = particles.begin(); it != particles.end(); ++it )
  {
    
    if ((*it).free) NfrozenOut++;
    
    //
    //transversePeriodicHack(it);
    //transverseReflectingHack(it);
    if (( *it ).Pos.T() < timenext )
    {
      NinAct++;
      ( *it ).init = false;

      if ( fabs(( *it ).Pos.X() - randomShiftX - halfsize ) < 1.0e-6 )
      {
        nx = IX - 1;
      }
      else
      {
        nx = static_cast<int>(((( *it ).Pos.X() - randomShiftX ) / transLen + 0.5 ) * IX );
      }

      if ( fabs(( *it ).Pos.Y() - randomShiftY - halfsize ) < 1.0e-6 )
      {
        ny = IY - 1;
      }
      else
      {
        ny = int(((( *it ).Pos.Y() - randomShiftY ) / transLen + 0.5 ) * IY );
      }

      nz = etaBins.getIndex(( *it ).eta );

      if (( nz >= IZ ) || ( nz < 0 ) )
      {
        cout << "eee  " << nz << "  " << IZ << "  " << (*it).eta << "  " << (*it).unique_id << endl;
        cout << "eee  " << etaBins.getCentralIndex() << "  " << etaBins.min_index_limit() << "  " << etaBins.max_index_limit() << "  " << etaBins.size() << endl;
        cout << "eee  " << (*it).Mom << endl;
        cout << "eee  " << (*it).Pos << endl;
        cout << "eee  " << timenow << "  " << (*it).free << "  " << (*it).edge  << "  " << (*it).dead << "  " << (*it).cell_id << endl;
        for ( int ii = 0; ii < etaBins.size(); ii++ )
        {
          cout << "eee# " << ii << "  " << etaBins[ii].left << "  " << etaBins[ii].right << "  " << etaBins[ii].content << endl;
        }
        std::string errMsg = "error in cell_id, nz out of range";
        throw eHIC_error( errMsg );
      }

      if (( nz < etaBins.min_index() ) || ( nz > etaBins.max_index() ) || ( nx < 0 ) || ( nx >= IX ) || ( ny < 0 ) || ( ny >= IY ) )
      {
        cell_id = -2;// -2:edge

        ( *it ).rate23 = 0;//GeV
        ( *it ).rate32 = 0;//GeV
        ( *it ).rate22 = 0;//GeV
        ( *it ).as22 = ( *it ).as23 = 0;
        ( *it ).cs22 = 0;//1/GeV^2
        ( *it ).cs23 = 0;//1/GeV^2
        ( *it ).md2g_scaled_22 = 0;
        ( *it ).md2q_scaled_22 = 0;
        ( *it ).md2g_scaled_23 = 0;
        ( *it ).md2q_scaled_23 = 0;
        ( *it ).lambda_scaled = 0;
      }
      else if (( nx < 0 ) || ( nx >= IX ) || ( ny < 0 ) || ( ny >= IY ) )
      {
        cell_id = -100;// -100:free
        Nfree1++;

        ( *it ).free = true;
        ( *it ).edge = -1;

        ( *it ).rate23 = 0;//GeV
        ( *it ).rate32 = 0;//GeV
        ( *it ).rate22 = 0;//GeV
        ( *it ).as22 = ( *it ).as23 = 0;
        ( *it ).cs22 = 0;//1/GeV^2
        ( *it ).cs23 = 0;//1/GeV^2
        ( *it ).md2g_scaled_22 = 0;
        ( *it ).md2q_scaled_22 = 0;
        ( *it ).md2g_scaled_23 = 0;
        ( *it ).md2q_scaled_23 = 0;
        ( *it ).lambda_scaled = 0;

        ( *it ).Propagate( timenext );
      }
      else
      {
        cell_id = nx + IX * ny + IX * IY * nz;
      }
      
      //Check if particle has changed the cell
      if ( cell_id != (*it).cell_id )
      {
        int old_cell_id = (*it).cell_id;
        //If old cell was a normal cell and no edge or free particle
        if ( old_cell_id >= 0 )
        {
          if ( theConfig->repeatTimesteps() )
          {
            auto lb=_cellsBackup.find(old_cell_id);
            if ( lb == _cellsBackup.end() )
            {
              _cellsBackup.insert( std::pair< int, cellContainer >( old_cell_id, cells[old_cell_id] ) );
            }
          }
          cells[old_cell_id].particleList.remove( static_cast<int>( it - particles.begin() ) );
        }
        //Set the new cell-ID
        ( *it ).cell_id = cell_id;
        //If this new cell-ID is a normal ID, insert the particle in the particle list:
        if ( cell_id >= 0 )
        {
          if ( theConfig->repeatTimesteps() )
          {
            auto lb=_cellsBackup.find(cell_id);
            if ( lb == _cellsBackup.end() )
            {
              _cellsBackup.insert( std::pair< int, cellContainer >( cell_id, cells[cell_id] ) );
            }
          }
          cells[cell_id].particleList.push_back( static_cast<int>( it - particles.begin() ) );
        }
        //BUT if this new cell is a Edge Cell, add or update edge cell:
        else if ( cell_id == -2 )
        {
          //IF the particle was not an edge cell member before, add it to the edge cells:
          if ( ( *it ).edge < 0 ) // new member
          {
            addParticleToEdgeCell( static_cast<int>( it - particles.begin() ), edgeCell );
            ( *it ).collisionTime = infinity;
            ( *it ).collisionPartner = -1;
            ( *it ).md2g = -1.0;
            ( *it ).md2q = -1.0;
          }
          //Otherwise, update the edge cell it belongs to
          else
          {
            int newEdgeCellIndex = getEdgeCellIndex( (*it) );
            if ( ( *it ).edge != newEdgeCellIndex )  // check if the edge particle needs to be put into another edge cell
            {
              int id = static_cast<int>( it - particles.begin() );
              swapParticleToNewEdgeCell( id, particles[id].edge, newEdgeCellIndex, edgeCell );
            } 
          }
        }
        //If the particle is free (cell-id=-100), remove it from any edge cell if it had belonged to one before:
        else if ( cell_id == -100 )
        {
          if ( (*it).edge >= 0 )
          {
            int id = static_cast<int>( it - particles.begin() );
            removeParticleFromEdgeCell( id, particles[id].edge, edgeCell );
            (*it).edge = -1;
          }
        }
        else
        {
          std::string errMsg = "error in cell_id, cell_id not defined";
          throw eHIC_error( errMsg );
        }
      }
    }
    else
    {
      if (( *it ).init )
      {
        NinFormInit++;
      }
      else
      {
        NinFormGeom++;
        formGeom.push_back( static_cast<int>( it - particles.begin() ) );
      }

      if ( (*it).cell_id != -1 )
      {
        int old_cell_id = (*it).cell_id;
        if ( old_cell_id >= 0 )
        {
          if ( theConfig->repeatTimesteps() )
          {
            auto lb=_cellsBackup.find(old_cell_id);
            if ( lb == _cellsBackup.end() )
            {
              _cellsBackup.insert( std::pair< int, cellContainer >( old_cell_id, cells[old_cell_id] ) );
            }
          }
          cells[old_cell_id].particleList.remove( static_cast<int>( it - particles.begin() ) );
        }
      }
      
      if ( (*it).edge >= 0 )
      {
        int id = static_cast<int>( it - particles.begin() );
        removeParticleFromEdgeCell( id, particles[id].edge, edgeCell );
      }

      ( *it ).cell_id = -1;
      ( *it ).edge = -1;

      ( *it ).rate23 = 0;//GeV
      ( *it ).rate32 = 0;//GeV
      ( *it ).rate22 = 0;//GeV
      ( *it ).as22 = ( *it ).as23 = 0;
      ( *it ).cs22 = 0;//1/GeV^2
      ( *it ).cs23 = 0;//1/GeV^2
      ( *it ).md2g_scaled_22 = 0;
      ( *it ).md2q_scaled_22 = 0;
      ( *it ).md2g_scaled_23 = 0;
      ( *it ).md2q_scaled_23 = 0;
      ( *it ).lambda_scaled = 0;
    }
  }
}

void heavyIonCollision::transversePeriodicHack(vector<Particle>::iterator it)
{
  double halfsize;

  halfsize = transLen / 2.0;
  
  if ( ( *it ).Pos.X()  >  halfsize )
  {
    ( *it ).Pos.X() = -halfsize;
  }
  
  if ( ( *it ).Pos.X()  < -halfsize )
  {
    ( *it ).Pos.X() = halfsize;
  }
  if ( ( *it ).Pos.Y()  >  halfsize )
  {
    ( *it ).Pos.Y() = -halfsize;
  }
  
  if ( ( *it ).Pos.Y()  < -halfsize )
  {
    ( *it ).Pos.Y() = halfsize;
  }  
//   cout << ( *it ).Pos.X() << endl;
}

void heavyIonCollision::transverseReflectingHack(vector<Particle>::iterator it)
{
  double halfsize;

  halfsize = transLen / 2.0;
  
  if ( ( *it ).Pos.X()  >  halfsize )
  {
    ( *it ).Mom.Px() *= -1.;
  }
  
  if ( ( *it ).Pos.X()  < -halfsize )
  {
    ( *it ).Mom.Px() *= -1.;
  }
  if ( ( *it ).Pos.Y()  >  halfsize )
  {
    ( *it ).Mom.Py() *= -1.;
  }
  
  if ( ( *it ).Pos.Y()  < -halfsize )
  {
    ( *it ).Mom.Py() *= -1.;
  }  
//   cout << ( *it ).Pos.X() << endl;
}

void heavyIonCollision::cleanFormGeom()
{ 
  double cc, eta;
  bool free;
  double zz;
  
  if ( !formGeom.empty() )
  {
    list<int>::iterator iIt;
    int id = -1;
    for ( iIt = formGeom.begin(); iIt != formGeom.end(); )
    {
      id = *iIt;
      cc = ( timenow - particles[id].Pos.T() ) / particles[id].Old.E();

      zz = particles[id].Pos.Z() + particles[id].Old.Pz() * cc;
      eta = 0.5 * log(( timenow + zz ) / ( timenow - zz ) );

      if (( eta < etaBins[ etaBins.min_index()].left ) || ( eta > etaBins[ etaBins.max_index()].right ) )
      {
        iIt = formGeom.erase( iIt );
      }
      else
      {
        ++iIt;
      }
    }
  }
}


void heavyIonCollision::collectParticlesInCellWithNeighbors(int cellindex, std::vector< int >& _allParticlesListWNeighbors, int etaSliceIndex, int IX, int IY,int & nCellsAVG)
{
  list<int>::const_iterator iIt;
  int totalNumber=0;
  int sizeTot=0;
  int id;
  int IXY=IX*IY;
  int n_z = cellindex / IXY; 
  int n_y = (cellindex-n_z*IXY) / IX;
  int n_x = (cellindex-n_z*IXY) - (n_y*IX);
  bool testWithSpecificNumbers =false;//default: false

  int rowCells=3;
  
  if( ( n_x > (rowCells-1)/2 ) && ( n_x < (IX-(rowCells-1)/2) ) && (n_y >(rowCells-1)/2) && (n_y < (IY-(rowCells-1)/2)) )
  {
    nCellsAVG=rowCells*rowCells;

    //cells in the middle
    for(int loop_y=0;loop_y<=(rowCells-1);loop_y++)
    {
      for(int loop_x=0;loop_x<=(rowCells-1);loop_x++)
      {
        int n_loop_x=(n_x-(rowCells-1)/2)+loop_x;
        int n_loop_y=(n_y-(rowCells-1)/2)+loop_y;
        int CELLINDEX_LOOP= n_loop_x + IX*n_loop_y + IXY*etaSliceIndex;
        
//         cout << loop_x << "\t" << loop_y<<"\t" << CELLINDEX_LOOP << endl;
//         totalNumber+=cells[CELLINDEX_LOOP].particleList.size();
//         cout << cells[CELLINDEX_LOOP].particleList.size() << endl;
        totalNumber+=cells[CELLINDEX_LOOP].size();
        for ( iIt = cells[CELLINDEX_LOOP].particleList.begin(); iIt != cells[CELLINDEX_LOOP].particleList.end(); iIt++ )
        {        
          id = *iIt;                   
          _allParticlesListWNeighbors.push_back( id );                  
        }
      }
    }     
  }
  else
  {
    nCellsAVG=0;
  }
  if(totalNumber<30)
  {
    nCellsAVG=0;   
  }  
}


void heavyIonCollision::populateRingAndCellAverages(double dz, int i)
{
  double cc,zz,eta;
  rings.clear();
  int IXY = IX * IY;
  
  theConfig->useRings();
  
  for ( int j = IXY * i; j < IXY * ( i + 1 ); j++ )
  { 
    if(theConfig->useRings())
    {
      //use ring-average for DEBYE mass and RATES
      list<int>::const_iterator iIt;
      for ( iIt = cells[j].particleList.begin(); iIt != cells[j].particleList.end(); iIt++ )
      {
        // calculation of energy density & co. is only written for massless particles, only add those -> for the medium properties heavy quarks do not contribute anyhow
        if( particles[( *iIt )].FLAVOR <= 2*Particle::max_N_light_flavor ) 
        {
          rings.addParticle( particles[( *iIt )] );
        }
        
        rings.addRates( particles[(*iIt)] );
      } 
    }else
    {
      //use cell-average for DEBYE mass and RATES    
      double volume = dx * dy * dz;//cells[j].corner.getVolume(etaBins,timenow);
      
      //HACK
      //if(cells[j].size() < 30)//WARNING 30
      {
        
        std::vector< int > listWithNeighbors;
        int nCellsAVG=0;
        collectParticlesInCellWithNeighbors(j, listWithNeighbors, i,IX,  IY, nCellsAVG);
        if(nCellsAVG>1)
        {
          double volumeCluster=nCellsAVG*volume;
          cells[j].clearThermodynamicsAndRates();      
          cells[j].prepareThermodynamicsCellBoostWNeighbors(volumeCluster,theConfig->getTestparticles(),listWithNeighbors);
          
          //cells[j].prepareThermodynamicsWithRadialBoostWNeighbors(volumeCluster,theConfig->getTestparticles(),listWithNeighbors);
          //cells[j].setValuesToRadialValues();
          cells[j].prepareAverageRatesWNeighbors(volumeCluster,theConfig->getTestparticles(),listWithNeighbors);        
//           cout << "Cluster" << "\t\t\t" << cells[j].corner.x_min << "\t" << cells[j].corner.y_min << endl;
        }else
        {
          cells[j].defineAsEmpty();
//           cout << "EMPTY" << "\t\t\t" << cells[j].corner.x_min << "\t" << cells[j].corner.y_min << endl;
        }
      }
//       else
//       {
// //         cout << "Single" << "\t\t\t" << cells[j].corner.x_min << "\t" << cells[j].corner.y_min << endl;
//         cells[j].clearThermodynamicsAndRates();      
//         cells[j].prepareThermodynamicsCellBoost(volume,theConfig->getTestparticles());
//         cells[j].prepareAverageRates(volume,theConfig->getTestparticles());
//         //Use cells, but boost radially only
//         //cells[j].prepareThermodynamicsWithRadialBoost(volume,theConfig->getTestparticles());
//       }
      
      
      

    }
  }
  
  if(theConfig->useRings())
  {  
    list<int>::iterator iIt;
    for ( iIt = formGeom.begin(); iIt != formGeom.end(); iIt++ )
    {
      int id = *iIt;

      cc = ( timenow - particles[id].Pos.T() ) / particles[id].Old.E();
      zz = particles[id].Pos.Z() + particles[id].Old.Pz() * cc;
      eta = 0.5 * log(( timenow + zz ) / ( timenow - zz ) );

      if (( eta >= etaBins[i].left ) && ( eta <= etaBins[i].right ) )
      {
        // calculation of energy density & co. is only written for massless particles, only add those -> for the medium properties heavy quarks do not contribute anyhow
        if( particles[( *iIt )].FLAVOR <= 2*Particle::max_N_light_flavor ) 
        {
          rings.addParticleInFormGeom( particles[id], timenow );
        }
        
        iIt = formGeom.erase( iIt );    // erase this particle such that it needs not be looped over for the next rapidity slab
      }
    }
    
    rings.prepareAverages( dz, testpartcl );  
    
    //TODO: offlineOutput uses only rings at the moment. Could be changed in future.
    for ( int jj = 0; jj < rings.size(); ++jj )
    {
      rateForOfflineOutput_gluons[i][jj] = rings[jj].rates.getRate( gluon, GeV );
      rateForOfflineOutput_quarks[i][jj] = rings[jj].rates.getRate( light_quark, GeV );
      rateForOfflineOutput_antiQuarks[i][jj] = rings[jj].rates.getRate( anti_light_quark, GeV );
    }
    //---------- populate the ring structure for averages ----------
  }
}


void heavyIonCollision::analyzeCentralCell( int cellindex ,std::vector< int >& _allParticlesList, int etaSliceIndex, int IX, int IY, const double volume, int nCellsAVG, double time , ringStructure&_rings)
{
  int id;
  int IXY=IX*IY;
  int n_z = cellindex / IXY; 
  int n_y = (cellindex-n_z*IXY) / IX;
  int n_x = (cellindex-n_z*IXY) - (n_y*IX);  
  double T_AMY=0.;
  lorentz LorentzBoost;
  int centralEtaIndex = static_cast<int>( etaBins.size() ) / 2;
  string full_filename = theConfig->getStandardOutputDirectoryName() + "/" +theConfig->getJobName() + "_T_AMY_Central" +  ".dat";
  T_AMY=0;

  fstream outfile( full_filename.c_str(), ios::out | ios::app);
  outfile.precision( 8 );
  
  //Central: IX/2-1
 
    if( (n_x==(IX/2)) &&  (n_y==IY/2-1) && (etaSliceIndex==centralEtaIndex) )
    {
      if(volume>0)
      {         
        if(!cells[cellindex].empty())
        {
          cout << cells[cellindex].gamma << "\t" << cells[cellindex].gammaCellRadial << endl;
          cout << "Edens c= " << cells[cellindex].getEnergyDensity() << endl;
          cout << "Edens r= " << cells[cellindex].getEnergyDensityCellRadial() << endl;
        }        
        
         cout << volume << "\t" << _allParticlesList.size() << endl;
         cout << cells[cellindex].getAveraged_md2g() << "\t\t" << _rings[0].getAveraged_md2g() <<endl;
         cout << cells[cellindex].getAveraged_md2q() << "\t\t" << _rings[0].getAveraged_md2q() <<endl;
         cout << "DENS= " << cells[cellindex].numberOfParticles/(cells[cellindex].volume*testpartcl) << "\t\t" << _rings[0].getGluonDensity() << endl;
         cout << "DENS= " << cells[cellindex].gluonDensity << "\t\t" << _rings[0].getGluonDensity() << endl;
          //cout << "Particles: " << _allParticlesList.size() << "\t\tT central= " << T_AMY << "\t\t\t nCellsAVG = " << nCellsAVG << endl;
          //outfile << time << "\t" << T_AMY << "\t" << LorentzBoost.gammaVal();
     
      }
        //outfile << time << "\t#\t#";
    }
  
  outfile.close();
}



/**
 * This routine first sets and computes several average values (ring structure). Then it loops over all cells and 
 * calls the scattering routines for 2->2, 2->3, 3->2 for each of these cells.
 *
 * @param[out] NinActCell Number of particles in "active" cells, i.e. in cells whose particle content is larger than ::cellcut
 * @param[out] Nfree2 Number of particles that are in cells in which all particles are "free"
 * @param aa A reference to the analysis object that handles output, data collection etc.
 * @param _cellsBackup The backup structure for cells whose content has changed within the current time step
 */
void heavyIonCollision::scattering( int& NinActCell, int& Nfree2, bool& again, analysis& aa, map< int, cellContainer >& _cellsBackup )
{
  int n32, nc;
  double dz, xt, cc, eta;
  bool free;
  double zz;
  vector<int> gluonList, allParticlesList;
  const string sep = "\t";
  countMFP = 0;
  MFPRatio = 0.;
  MFPfm = 0.;
  int nGluons = 0;
  int nAllLightQuarks = 0;
  int nAllAntiLightQuarks = 0;
  again = false;
  int IXY = IX * IY;
  Nfree2 = NinActCell = 0;

  cleanFormGeom();
  
  
  //TEST
//   int nTot=0;
//   for ( int i = etaBins.min_index(); i <= etaBins.max_index(); i++ )
//   {
//     for ( int j = IXY * i; j < IXY * ( i + 1 ); j++ )
//     { 
//       nTot+=cells[j].particleList.size();
//       for(std::list<int>::iterator iIt = cells[j].particleList.begin(); iIt != cells[j].particleList.end(); iIt++ )
//       {
//          particles[*iIt].Propagate( timenext );
//       }
//     }
//   }
//   cout << "nTOT from Cells: " << nTot << endl;
//   for( int i = 0; i< particles.size();i++)
//   {
// //     particles[i].Propagate( timenext );
//     cout << particles[i].Mom << endl;
//   }
  
  
  for ( int i = etaBins.min_index(); i <= etaBins.max_index(); i++ )
  {
    //---------- populate the ring structure for averages ----------
    dz = timenow * ( tanh( etaBins[i].right ) - tanh( etaBins[i].left ) );
    dv = dx * dy * dz;  
    populateRingAndCellAverages(dz,i);

    for ( int j = IXY * i; j < IXY * ( i + 1 ); j++ )
    {   
      if ( !cells[j].empty() )
      {
        cells[j].resetStoredValues();
        gluonList.clear();
        allParticlesList.clear();
        gluonList.reserve( cells[j].size() );
        allParticlesList.reserve( cells[j].size() );

        //WARNING CELLWISE OR RINGWISE FREEZEOUT ?!
        if(theConfig->useRings())
        {
          findFreeParticles(free,j);
        }else
        {
          findFreeParticlesCellwise(free,j);
        }

        
        if ( free )    // ALL particles in this cell are free
        {
          handleAllParticlesAreFree(Nfree2,j);
        }
        else
        {
          //normal loop through cells
          if ( cells[j].size() >= cellcut )
          {       
            int nCharmQuarks = 0;
            int nAntiCharmQuarks = 0;
            int nBottomQuarks = 0;
            int nAntiBottomQuarks = 0;
            handleParticleNumbersRatesDebyeMasses(j,NinActCell,nGluons,free,allParticlesList,gluonList,nAllLightQuarks,nAllAntiLightQuarks,nCharmQuarks,nAntiCharmQuarks,nBottomQuarks,nAntiBottomQuarks);

            n32 = 0;
            if ( theConfig->doScattering_32() )
            {
              scatt32( cells[j], allParticlesList, gluonList, n32, again, aa );
            }
            if ( theConfig->repeatTimesteps() && again )
            {
              if ( theConfig->doOutput_offlineReconstructionData() )
              { 
                offlineInterface.resetTemporaryStorage();     
              }
              return;
            }
            
            double scaleFactor;
            if( ( nGluons - n32 ) > 0 )
            {
              scaleFactor = static_cast<double>( nGluons ) / static_cast<double>( nGluons - n32 );
            }
            else if( nGluons - n32 == 0) // no gluons in the cell anymore (all of them have been absorbed in 3->2 collisions). In this case set scaleFactor to 1. But it does not really matter since the scaleFactor is only used if a gluon is involved in a 2->x processes which cannot happen since no gluons are present anymore.
            {
              scaleFactor = 1.0;
            }
            else 
            {
              string errMsg = "n32 > nGluons!";
              throw eHIC_error( errMsg );
            }
           
            //analyzeCentralCell(j,allParticlesList,i,IX,IY,dv,1,timenext,rings);
            scatt2322( cells[j], allParticlesList, gluonList, scaleFactor, again, aa );
            
            if ( theConfig->repeatTimesteps() && again )
            {
              if ( theConfig->doOutput_offlineReconstructionData() )
              {
                offlineInterface.resetTemporaryStorage();
              }
              return;
            }

            cells[j].rates.normalizeRates( nGluons, nAllLightQuarks, nAllAntiLightQuarks, nCharmQuarks, nAntiCharmQuarks, nBottomQuarks, nAntiBottomQuarks, dt, dv );
            cells[j].prepareAverages();
              
            if(theConfig->useRings())
            {
              updateParticleRatesRings(j);
            }else
            {
              updateParticleRatesCellwise(j);
            }
          }
          else   //belongs to if(nn >= cellcut) -> these particles go into the edge cells
          {
            cout << "Cellcut!" << endl;
            //WARNING CELLWISE OR RINGWISE
            if(theConfig->useRings())
            {
              handleCellcutCellsRingwise(j);
            }else
            {
              handleCellcutCellsCellwise(j);
            }
            
          }
        }
      }
      else
      {
//         cells[j].rates.normalizeRates();
      }
    }
  }  
  
  if(theConfig->getCrossSectionMethod()==csMethod_constEtaOverS)
  {
    cout << "MFP/fm = " << MFPfm/countMFP  << "\tmfp/dz " << MFPRatio/countMFP << "\t\t\t dz[fm] " << dz <<endl;
  }
  
  if ( theConfig->doOutput_offlineReconstructionData() )
  {
    offlineInterface.outputAndResetTemporaryStorage();
  }  
}

void heavyIonCollision::handleParticleNumbersRatesDebyeMasses(int & cellindex, int & NinActCell,int & nGluons,bool & free, std::vector< int >&   allParticlesList, std::vector< int >&   gluonList, int & nAllLightQuarks, int & nAllAntiLightQuarks, int & nCharmQuarks, int & nBottomQuarks, int & nAntiCharmQuarks, int & nAntiBottomQuarks)
{
  int nc;
  double xt;
  NinActCell += cells[cellindex].size();
  nGluons = 0;
  free = false;
  vector<int> nLightQuarks( Particle::N_light_flavor , 0 );
  vector<int> nAntiLightQuarks( Particle::N_light_flavor, 0 );

  list<int>::const_iterator iIt;
  for ( iIt = cells[cellindex].particleList.begin(); iIt != cells[cellindex].particleList.end(); iIt++ )
  {
    int id = *iIt;

    particles[id].edge = -1;
    allParticlesList.push_back( id );

    if ( particles[id].FLAVOR == 0 )
    {
      nGluons++;
      gluonList.push_back( *iIt );
    }
    else
    {
      switch ( particles[id].FLAVOR )
      {
      case up:
        ++nLightQuarks[0];
        break;
      case down:
        ++nLightQuarks[1];
        break;
      case strange:
        ++nLightQuarks[2];
        break;
      case charm:
        ++nCharmQuarks;
        break;
      case bottom:
        ++nBottomQuarks;
        break;
      case anti_up:
        ++nAntiLightQuarks[0];
        break;
      case anti_down:
        ++nAntiLightQuarks[1];
        break;
      case anti_strange:
        ++nAntiLightQuarks[2];
        break;
      case anti_charm:
        ++nAntiCharmQuarks;
        break;
      case anti_bottom:
        ++nAntiBottomQuarks;
        break;
      default:
        break;
      }
    }

    xt = particles[id].Pos.Perp();
    nc = rings.getIndex( xt );

    if(theConfig->useRings())
    { 
      //WARNING RINGS ARE USED:
      particles[id].rate22 = rings[nc].rates.getRate( particles[id].FLAVOR, c22, GeV );
      particles[id].rate23 = rings[nc].rates.getRate( particles[id].FLAVOR, c23, GeV );
      particles[id].rate32 = rings[nc].rates.getRate( particles[id].FLAVOR, c32, GeV );             
      particles[id].md2g = rings[nc].getAveraged_md2g();
      particles[id].md2q = rings[nc].getAveraged_md2q();                
    }else
    {
      //WARNING this uses cell-based-averages for RATES and DEBYE-MASSES instead of rings.
      particles[id].md2g =    cells[cellindex].getAveraged_md2g();
      particles[id].md2q =    cells[cellindex].getAveraged_md2q();             
      particles[id].rate22 =  cells[cellindex].ratesSave.getRate( particles[id].FLAVOR, c22, GeV );
      particles[id].rate23 =  cells[cellindex].ratesSave.getRate( particles[id].FLAVOR, c23, GeV );
      particles[id].rate32 =  cells[cellindex].ratesSave.getRate( particles[id].FLAVOR, c32, GeV );                   
    }
  }
  
  nAllLightQuarks = std::accumulate( nLightQuarks.begin(), nLightQuarks.end(), 0 );
  nAllAntiLightQuarks = std::accumulate( nAntiLightQuarks.begin(), nAntiLightQuarks.end(), 0 );  
}

void heavyIonCollision::handleAllParticlesAreFree(int &Nfree2, int cellindex)
{
  double xt;
  int nc;
  Nfree2 += cells[cellindex].size();

  list<int>::const_iterator iIt;
  for ( iIt = cells[cellindex].particleList.begin(); iIt != cells[cellindex].particleList.end(); iIt++ )
  {
    int id = *iIt;
    particles[id].edge = -1;

    if(theConfig->useRings())
    {
     //IF RINGS ARE USED:
      xt = particles[id].Pos.Perp();
      nc = rings.getIndex( xt );
      particles[id].rate22v = rings[nc].rates.getRate( particles[id].FLAVOR, c22, GeV );
      particles[id].rate23v = rings[nc].rates.getRate( particles[id].FLAVOR, c23, GeV );
      particles[id].rate32v = rings[nc].rates.getRate( particles[id].FLAVOR, c32, GeV );
    }
    else
    {
      //WARNING HACK this uses cell-based-averages for RATES instead of rings.           
      particles[id].rate22v =  cells[cellindex].ratesSave.getRate( particles[id].FLAVOR, c22, GeV );
      particles[id].rate23v =  cells[cellindex].ratesSave.getRate( particles[id].FLAVOR, c23, GeV );
      particles[id].rate32v =  cells[cellindex].ratesSave.getRate( particles[id].FLAVOR, c32, GeV );         
    }
    particles[id].rate23 = 0;   //GeV
    particles[id].rate32 = 0;   //GeV
    particles[id].rate22 = 0;   //GeV
    particles[id].Propagate( timenext );
  }            
  cells[cellindex].rates.normalizeRates();
  cells[cellindex].defineAsEmpty();
}

/**
 * Given a certain cell and time step, this routine computes does the Monte Carlo sampling of 3->2 interactions of 
 * particle triplets within this cell.
 * 
 * @param _cell The cell for which the 3->2 interactions should be computed
 * @param[in] allParticlesList The list of all particles in the cell (TODO Redundant information)
 * @param[in] gluonList The list of all gluons in the cell
 * @param[out] n32 The number of 3->2 processes that actually took place
 * @param[out] again Whether the time step needs to be undone because a probability in this routine was found to be >1
 * @param aa A reference to the analysis object that handles output, data collection etc.
 */
void heavyIonCollision::scatt32( cellContainer& _cell, vector<int>& allParticlesList, vector<int>& gluonList, int& n32, bool& again, analysis& aa )
{
  const double epsilon = 1.0e-4;
  int iscat, jscat, kscat;
  unsigned int m1, m2, m3;
  int ringIndex;
  int consideredTriplets;  //number of triplets to consider

  FLAVOR_TYPE F1, F2, F3;

  double averagedRate;
  double s, csgg, I32, probab32, xt;
  double md2g_wo_as, md2q_wo_as;
  double lambda_scaled;
  double ran2out;
  double betaDistEntry;

  //n32=0;

  const int nTotal = _cell.particleList.size();

  list<int>::const_iterator iIt;
  const int nGluons = gluonList.size();
  const int nAllQuarks = allParticlesList.size() - nGluons;
  const int allTriplets = binomial( nTotal, 3 )  - binomial( nAllQuarks, 3 );

  scattering32 scatt32_object;

  if ( allTriplets > 20 )
  {
    if ( nTotal >= 20 )
    {
      consideredTriplets = nTotal;  // number of triplets to consider for 3 -> 2 processes
    }
    else
    {
      consideredTriplets = 20;  //20 is an empirical value
    }

    double scaleForSelectedTriplets = static_cast<double>( allTriplets ) / static_cast<double>( consideredTriplets );

    for ( int i = 0; i < consideredTriplets && gluonList.size() > 0; i++ )
    {
      // the first particle is forced to be a gluon in order to reduce unnecessary samplings
      m1 = int ( gluonList.size() * ran2() );
      
      if ( m1 == gluonList.size() )
        m1 = gluonList.size() - 1;
      iscat = gluonList[m1];
      F1 = particles[iscat].FLAVOR;
      
      do
      {
        m2 = int ( allParticlesList.size() * ran2() );
        
        if ( m2 == allParticlesList.size() )
          m2 = allParticlesList.size() - 1;
        
        jscat = allParticlesList[m2];
      }
      while ( jscat == iscat );
      F2 = particles[jscat].FLAVOR;
      
      do
      {
        m3 = int ( allParticlesList.size() * ran2() );
        
        if ( m3 == allParticlesList.size() )
          m3 = allParticlesList.size() - 1;
        
        kscat = allParticlesList[m3];
      }
      while (( kscat == iscat ) || ( kscat == jscat ) );
      F3 = particles[kscat].FLAVOR;

      s = (particles[iscat].Mom + particles[jscat].Mom + particles[kscat].Mom).M2();

      averagedRate = ( particles[iscat].rate22 + particles[iscat].rate23 + particles[iscat].rate32 +
                       particles[jscat].rate22 + particles[jscat].rate23 + particles[jscat].rate32 +
                       particles[kscat].rate22 + particles[kscat].rate23 + particles[kscat].rate32 +
                       particles[iscat].rate22v + particles[iscat].rate23v + particles[iscat].rate32v +
                       particles[jscat].rate22v + particles[jscat].rate23v + particles[jscat].rate32v +
                       particles[kscat].rate22v + particles[kscat].rate23v + particles[kscat].rate32v ) / ( 3.0 * 2.0 );

      if ( s < 1.1*lambda2 )
        probab32 = -1.0;
      else
      {
        //factor as (alpha_s) is not included in definitions of partcl[jscat].md2g
        md2g_wo_as = ( particles[iscat].md2g + particles[jscat].md2g + particles[kscat].md2g ) / 3.0;
        md2q_wo_as = ( particles[iscat].md2q + particles[jscat].md2q + particles[kscat].md2q ) / 3.0;

        xt = ( particles[iscat].Pos.Perp() +
               particles[jscat].Pos.Perp() +
               particles[kscat].Pos.Perp() ) / 3.0;

        ringIndex = rings.getIndex( xt );
        
        //WORKAROUND: set rings to value of cell
        if(!theConfig->useRings())
        {
          rings.setValuesAsWorkaround(ringIndex, _cell.getEnergyDensity(),_cell.getParticleDensity(),_cell.getGamma(),_cell.getNumberOfParticles(),_cell.getGluonDensity(),_cell.getQuarkDensity());
        }

        // should the jet mean path computed more accurately for light partons?
        double interpolLimit = theConfig->getMFPInterpolationBorder() * rings[ringIndex].getEffectiveTemperature();
        if ( ( F1 <= 2*Particle::max_N_light_flavor && F2 <= 2*Particle::max_N_light_flavor  && F3 <= 2*Particle::max_N_light_flavor ) && theConfig->getJetMfpComputationType() != computeMfpDefault )
        {
          double E1_dash = rings[ringIndex].transformEnergyToComovingFrame( particles[iscat].Mom );
          double E2_dash = rings[ringIndex].transformEnergyToComovingFrame( particles[jscat].Mom );
          double E3_dash = rings[ringIndex].transformEnergyToComovingFrame( particles[kscat].Mom );
          
          if ( !( E1_dash < interpolLimit && E2_dash < interpolLimit && E3_dash < interpolLimit ) )
          {
            int tempJetID;
            double tempJetEnergy;
            FLAVOR_TYPE tempJetFlavor;
            if ( E1_dash >= E2_dash && E1_dash >= E3_dash )
            {
              tempJetID = iscat;
              tempJetEnergy = E1_dash;
              tempJetFlavor = F1;
            }
            else if ( E2_dash > E3_dash ) 
            {
              tempJetID = jscat;
              tempJetEnergy = E2_dash;
              tempJetFlavor = F2;
            }
            else
            {
              tempJetID = kscat;
              tempJetEnergy = E3_dash;
              tempJetFlavor = F3;
            }
            
            if ( theConfig->getJetMfpComputationType() == computeMfpInterpolation )
            {
              averagedRate = 1 / ( theMFP->getMeanFreePath( tempJetEnergy, tempJetFlavor, rings[ringIndex].getEffectiveTemperature(),
                                                            rings[ringIndex].getGluonDensity(), rings[ringIndex].getQuarkDensity(), GeV ) );
            }
            else if ( theConfig->getJetMfpComputationType() == computeMfpIteration )
            {
              averagedRate = 1 / iterateMFP( allParticlesList, gluonList, tempJetID, dt, dv, timenow );
            }
            else
            {
              string errMsg = "jetMfpComputationType not valid";
              throw eHIC_error( errMsg );
            }
          }
        }

        if ( averagedRate > epsilon )
        {
          lambda_scaled = sqrt( s ) / averagedRate;
        }
        else
        {
          xsection_gg_gg csObj( s, md2g_wo_as, md2q_wo_as, &theI22, theConfig->getKfactor_light() );
          csgg = csObj.totalCrossSection();
          lambda_scaled = ( dv * rings[ringIndex].getGamma() * testpartcl * sqrt( s ) ) / ( nTotal * csgg * pow( 0.197, 3.0 ) );
        }

        // create scattering32 object for the given 3 particles
        betaDistEntry = scatt32_object.setParameter( rings[ringIndex].getAveraged_v(), 
                                                     particles[iscat].Mom, particles[jscat].Mom, particles[kscat].Mom,
                                                     F1, F2, F3, sqrt( s ), md2g_wo_as * coupling::get_constant_coupling() / s, lambda_scaled, coupling::get_constant_coupling(), theConfig->isMd2CounterTermInI23(), theConfig->isMatrixElement23_22qt(), theConfig->get23FudgeFactorLpm(), nGluons );  // create scattering32 object for the given 3 particles
        aa.addBoost32( betaDistEntry );
        I32 = scatt32_object.getIntegral32_withPrefactors();                        // get the integral I32 for the given 3 particles

        probab32 = scaleForSelectedTriplets * I32 * dt / ( pow( dv, 2 ) * pow( testpartcl, 2 ) );
        
        _cell.rates.add( c32, F1, F2, F3, probab32 );
      }

      if ( probab32 > 1.0 )
      {
        if ( theConfig->repeatTimesteps() )
        {
          cout << "P32 = " << probab32 << " > 1" << endl;
          cout << "dt (old) = " << dt << endl;
          again = true;
          dt = 0.5 / ( probab32 / dt );
          cout << "dt (new) = " << dt << endl;
          return;
        }
        else
        {
          probab32 = 1;
        }
      }
      

      ran2out = ran2();
      if ( ran2out < probab32 )
      {        
        double pt_iscat = particles[iscat].Mom.Pt();
        double pt_jscat = particles[jscat].Mom.Pt();
        double pt_kscat = particles[kscat].Mom.Pt();

        int jetEventIndex = -1;
        if ( pt_iscat > aa.getJetTracking_PT() || pt_jscat > aa.getJetTracking_PT() || pt_kscat > aa.getJetTracking_PT() )
        {
          jetEventIndex = aa.addJetEvent_in( iscat, jscat, kscat, c3to2, I32, _cell.index, lambda_scaled / sqrt( s ) );
        }
        
        scatt32_utility( scatt32_object, _cell.particleList, allParticlesList, gluonList, iscat, jscat, kscat, n32 );

        // TO BE CHECKED: Do we really have to recalculate pT in the
        // following ????
        // Also in a copy of this code somwhere below

        if ( particles[iscat].dead )
        {
          pt_jscat = particles[jscat].Mom.Pt();
          pt_kscat = particles[kscat].Mom.Pt();
          if ( jetEventIndex != -1 || pt_jscat > aa.getJetTracking_PT() || pt_kscat > aa.getJetTracking_PT() )
          {
            aa.addJetEvent_out( jetEventIndex, jscat, kscat, -1 );
          }
        }
        else if ( particles[jscat].dead )
        {
          pt_iscat = particles[iscat].Mom.Pt();
          pt_kscat = particles[kscat].Mom.Pt();
          if ( jetEventIndex != -1 || pt_iscat > aa.getJetTracking_PT() || pt_kscat > aa.getJetTracking_PT() )
          {
            aa.addJetEvent_out( jetEventIndex, iscat, kscat, -1 );
          }
        }
        else
        {
          pt_iscat = particles[iscat].Mom.Pt();
          pt_jscat = particles[jscat].Mom.Pt();
          if ( jetEventIndex != -1 || pt_iscat > aa.getJetTracking_PT() || pt_jscat > aa.getJetTracking_PT() )
          {
            aa.addJetEvent_out( jetEventIndex, iscat, jscat, -1 );
          }
        }
      }
    }
  }
  else
  {
    for ( unsigned int m1 = 0; m1 < allParticlesList.size() - 2; m1++ )
    {
      iscat = allParticlesList[m1];

      for ( unsigned int m2 = m1 + 1; m2 < allParticlesList.size() - 1; m2++ )
      {
        jscat = allParticlesList[m2];

        for ( unsigned int m3 = m2 + 1; m3 < allParticlesList.size(); m3++ )
        {
          kscat = allParticlesList[m3];

          F1 = particles[iscat].FLAVOR;
          F2 = particles[jscat].FLAVOR;
          F3 = particles[kscat].FLAVOR;

          // at least one of the particles must be a gluon
          // otherwise go to next step in the loop
          if ( !( F1 == gluon || F2 == gluon || F3 == gluon ) )
          {
            continue;
          }

          s = (particles[iscat].Mom + particles[jscat].Mom + particles[kscat].Mom).M2();

          averagedRate = ( particles[iscat].rate22 + particles[iscat].rate23 + particles[iscat].rate32 +
                           particles[jscat].rate22 + particles[jscat].rate23 + particles[jscat].rate32 +
                           particles[kscat].rate22 + particles[kscat].rate23 + particles[kscat].rate32 +
                           particles[iscat].rate22v + particles[iscat].rate23v + particles[iscat].rate32v +
                           particles[jscat].rate22v + particles[jscat].rate23v + particles[jscat].rate32v +
                           particles[kscat].rate22v + particles[kscat].rate23v + particles[kscat].rate32v ) / ( 3.0 * 2.0 );

          if ( s < 1.1*lambda2 )
          {
            probab32 = -1.0;
          }
          else
          {
            //factor as (alpha_s) is not included in definitions of partcl[jscat].md2g
            md2g_wo_as = ( particles[iscat].md2g + particles[jscat].md2g + particles[kscat].md2g ) / 3.0;
            md2q_wo_as = ( particles[iscat].md2q + particles[jscat].md2q + particles[kscat].md2q ) / 3.0;

            xt = ( particles[iscat].Pos.Perp() +
                   particles[jscat].Pos.Perp() +
                   particles[kscat].Pos.Perp() ) / 3.0;

            ringIndex = rings.getIndex( xt );
            //WORKAROUND: set rings to value of cell
            if(!theConfig->useRings())
            {
              rings.setValuesAsWorkaround(ringIndex, _cell.getEnergyDensity(),_cell.getParticleDensity(),_cell.getGamma(),_cell.getNumberOfParticles(),_cell.getGluonDensity(),_cell.getQuarkDensity());
            }
            
            // should the jet mean path computed more accurately for light partons?
            double interpolLimit = theConfig->getMFPInterpolationBorder() * rings[ringIndex].getEffectiveTemperature();
            if ( ( F1 <= 2*Particle::max_N_light_flavor && F2 <= 2*Particle::max_N_light_flavor  && F3 <= 2*Particle::max_N_light_flavor ) && theConfig->getJetMfpComputationType() != computeMfpDefault )
            { 
              double E1_dash = rings[ringIndex].transformEnergyToComovingFrame( particles[iscat].Mom );
              double E2_dash = rings[ringIndex].transformEnergyToComovingFrame( particles[jscat].Mom );
              double E3_dash = rings[ringIndex].transformEnergyToComovingFrame( particles[kscat].Mom );
              
              if ( !( E1_dash < interpolLimit && E2_dash < interpolLimit && E3_dash < interpolLimit ) )
              {
                int tempJetID;
                double tempJetEnergy;
                FLAVOR_TYPE tempJetFlavor;
                if ( E1_dash >= E2_dash && E1_dash >= E3_dash )
                {
                  tempJetID = iscat;
                  tempJetEnergy = E1_dash;
                  tempJetFlavor = F1;
                }
                else if ( E2_dash > E3_dash ) 
                {
                  tempJetID = jscat;
                  tempJetEnergy = E2_dash;
                  tempJetFlavor = F2;
                }
                else
                {
                  tempJetID = kscat;
                  tempJetEnergy = E3_dash;
                  tempJetFlavor = F3;
                }
                
                if ( theConfig->getJetMfpComputationType() == computeMfpInterpolation )
                {
                  averagedRate = 1 / ( theMFP->getMeanFreePath( tempJetEnergy, tempJetFlavor, rings[ringIndex].getEffectiveTemperature(),
                                                                rings[ringIndex].getGluonDensity(), rings[ringIndex].getQuarkDensity(), GeV ) );
                }
                else if ( theConfig->getJetMfpComputationType() == computeMfpIteration )
                {
                  averagedRate = 1 / iterateMFP( allParticlesList, gluonList, tempJetID, dt, dv, timenow );
                }
                else
                {
                  string errMsg = "jetMfpComputationType not valid";
                  throw eHIC_error( errMsg );
                }
                
              }
            }

            if ( averagedRate > epsilon )
            {
              lambda_scaled = sqrt( s ) / averagedRate;
            }
            else
            {
              xsection_gg_gg csObj( s, md2g_wo_as, md2q_wo_as, &theI22, theConfig->getKfactor_light() );
              csgg = csObj.totalCrossSection();
              lambda_scaled = ( dv * rings[ringIndex].getGamma() * testpartcl * sqrt( s ) ) / ( nTotal * csgg * pow( 0.197, 3.0 ) );
            }
            
            // create scattering32 object for the given 3 particles
            betaDistEntry = scatt32_object.setParameter( rings[ringIndex].getAveraged_v(), 
                                                         particles[iscat].Mom, particles[jscat].Mom, particles[kscat].Mom,
                                                         F1, F2, F3, sqrt( s ), md2g_wo_as * coupling::get_constant_coupling() / s, lambda_scaled, coupling::get_constant_coupling(), theConfig->isMd2CounterTermInI23(), theConfig->isMatrixElement23_22qt(), theConfig->get23FudgeFactorLpm(), nGluons );  // create scattering32 object for the given 3 particles
            aa.addBoost32( betaDistEntry );
            I32 = scatt32_object.getIntegral32_withPrefactors();                        // get the integral I32 for the given 3 particles

            probab32 = I32 * dt / ( pow( dv, 2 ) * pow( testpartcl, 2 ) );
            
            _cell.rates.add( c32, F1, F2, F3, probab32 );
          }

          if ( probab32 > 1.0 )
          {
            if ( theConfig->repeatTimesteps() )
            {
              cout << "P32 = " << probab32 << " > 1" << endl;
              cout << "dt (old) = " << dt << endl;
              again = true;
              dt = 0.5 / ( probab32 / dt );
              cout << "dt (new) = " << dt << endl;
              return;
            }
            else
            {
              probab32 = 1;
            }
          }

          ran2out = ran2();
          if ( ran2out < probab32 )
          {            
            double pt_iscat = particles[iscat].Mom.Pt();
            double pt_jscat = particles[jscat].Mom.Pt();
            double pt_kscat = particles[kscat].Mom.Pt();

            int jetEventIndex = -1;
            if ( pt_iscat > aa.getJetTracking_PT() || pt_jscat > aa.getJetTracking_PT() || pt_kscat > aa.getJetTracking_PT() )
            {
              jetEventIndex = aa.addJetEvent_in( iscat, jscat, kscat, c3to2, I32, _cell.index, lambda_scaled / sqrt( s ) );
            }
            
            scatt32_utility( scatt32_object, _cell.particleList, allParticlesList, gluonList, iscat, jscat, kscat, n32 );

            if ( particles[iscat].dead )
            {
              pt_jscat = particles[jscat].Mom.Pt();
              pt_kscat = particles[kscat].Mom.Pt();
              if ( jetEventIndex != -1 || pt_jscat > aa.getJetTracking_PT() || pt_kscat > aa.getJetTracking_PT() )
              {
                aa.addJetEvent_out( jetEventIndex, jscat, kscat, -1 );
              }
                
              m2 = m1 + 1;
              m3 = m2; // +1 is added in next step of the loop
              iscat = allParticlesList[m1]; // old particle is removed, update jscat with the new particle id at the same place in the list
              jscat = allParticlesList[m2]; // update jscat
            }
            else if ( particles[jscat].dead )
            {
              pt_iscat = particles[iscat].Mom.Pt();
              pt_kscat = particles[kscat].Mom.Pt();
              if ( jetEventIndex != -1 || pt_iscat > aa.getJetTracking_PT() || pt_kscat > aa.getJetTracking_PT() )
              {
                aa.addJetEvent_out( jetEventIndex, iscat, kscat, -1 );
              }
                
              m3 = m2; // +1 is added in next step of the loop
              jscat = allParticlesList[m2];  // old particle is removed, update jscat with the new particle id at the same place in the list
            }
            else if ( particles[kscat].dead )
            {                
              pt_iscat = particles[iscat].Mom.Pt();
              pt_jscat = particles[jscat].Mom.Pt();
              if ( jetEventIndex != -1 || pt_iscat > aa.getJetTracking_PT() || pt_jscat > aa.getJetTracking_PT() )
              {
                aa.addJetEvent_out( jetEventIndex, iscat, jscat, -1 );
              }
                
              m3--;  // reduce m3 by one since in next step of loop it is increased by +1
            }
          }
        }
      }
    }
  }
}



/**
 * Given a certain cell and time step, this routine computes does the Monte Carlo sampling of 2->2 and 2->3 interactions of 
 * particle pairs within this cell.
 * 
 * @param _cell The cell for which the 3->2 interactions should be computed
 * @param[in] allParticlesList The list of all particles in the cell (TODO Redundant information)
 * @param[in] gluonList The list of all gluons in the cell
 * @param[in] scaleFactor As 3->2 interactions are computed first in heavyIonCollision::scattering, the probability for a given particle pair might need to be adjusted
 * @param[out] again Whether the time step needs to be undone because a probability in this routine was found to be >1
 * @param aa A reference to the analysis object that handles output, data collection etc.
 */
void heavyIonCollision::scatt2322( cellContainer& _cell, std::vector< int >& allParticlesList, std::vector< int >& gluonList, const double scaleFactor, bool& again, analysis& aa )
{
  const double epsilon = 1.0e-4;
  int iscat, jscat, typ;
  double s, csgg, cs23, cs22, Vrel, lambda_scaled;
  double M1, M2;
  double probab22, probab23, probab2322;
  double averagedRate;
  double xt;
  double betaDistEntry;
  double md2g_wo_as, md2q_wo_as;
  int ringIndex;
  int initialStateIndex = -1;
  FLAVOR_TYPE F1, F2;
  double specialCS, scalingFactor;  

  const int nTotal = _cell.particleList.size();

  int nGluons = 0;
//   vector<int> nQuarks( Nflavor, 0 );
//   vector<int> nAntiQuarks( Nflavor, 0 );

  list<int>::const_iterator iIt;
  list<int>::const_iterator jIt;

  //---------------------------------------------
  scattering23 scatt23_object( &theI23_massless, &theI23_charm_m1, &theI23_charm_m2, &theI23_bottom_m1, &theI23_bottom_m2 );
  scattering22 scatt22_object( &theI22 );
  scattering22_hydro scatt22_hydro;
  //---------------------------------------------  

//  const int nAllQuarks = std::accumulate( nQuarks.begin(), nQuarks.end(), 0 ) + std::accumulate( nAntiQuarks.begin(), nAntiQuarks.end(), 0 );
  const int allPairs = binomial( nTotal, 2 );
  const int consideredPairs = 5000;
  unsigned int m1, m2;
  double scaleForSelectedPairs;
  scaleForSelectedPairs =  static_cast<double>( allPairs ) / static_cast<double>( consideredPairs );
  //---------------------------------//
  
  //static definition! hydroParticleTypeVectorCommon
  
  vector<hydroParticleTypeProperties> vecTypeCommon(theHydroParticleType.hydroParticleTypeVectorCommon);
  
  //cout << vecTypeCommon[0].degeneracyFactor << endl;
  double TeffLocalCell =0.;
  specialCS = getSpecialCS( _cell, dv, theHydroParticleType, vecTypeCommon, scalingFactor, TeffLocalCell );
  //---------------------------------//  
  
  if ( allPairs > consideredPairs && (theConfig->getInitialStateType()==CYMInitialState) )  //only a certain number consideredPairs of all pairs will be considered. Only necessary for huge gradients, as in pPb with CYM input.
  {    
    for ( int i = 0; i < consideredPairs; i++ )
    {
      //------- randomly pick first particle -------
      iIt =  _cell.particleList.begin();
      m1 = int ( _cell.particleList.size() * ran2() );
      if(m1==_cell.particleList.size())
      {
        m1=_cell.particleList.size()-1;
      }
      std::advance( iIt, m1 );
      iscat = *iIt;
      //------- randomly pick second particle -------
      do
      {
        jIt =  _cell.particleList.begin();
        m2 = int ( _cell.particleList.size() * ran2() );
        if(m2==_cell.particleList.size())
        {
          m2=_cell.particleList.size()-1;
        }
        std::advance( jIt, m2 );
        jscat = *jIt;      
      }
      while ( jscat == iscat );
      
      //Do the scattering of the selected particles:
      {
        F1 = particles[iscat].FLAVOR;
        M1 = particles[iscat].m;
        F2 = particles[jscat].FLAVOR;
        M2 = particles[jscat].m;

        s = (particles[iscat].Mom + particles[jscat].Mom).M2();

        xt = ( particles[iscat].Pos.Perp() +
               particles[jscat].Pos.Perp() ) / 2.0;

        ringIndex = rings.getIndex( xt );
        
        //WORKAROUND: set rings to value of cell
        if(!theConfig->useRings())
        {
          rings.setValuesAsWorkaround(ringIndex, _cell.getEnergyDensity(),_cell.getParticleDensity(),_cell.getGamma(),_cell.getNumberOfParticles(),_cell.getGluonDensity(),_cell.getQuarkDensity());
        }

        
        if ( ( F1 <= 2*Particle::max_N_light_flavor && F2 <= 2*Particle::max_N_light_flavor ) && theConfig->getJetMfpComputationType() != computeMfpDefault )
        { 
          double interpolLimit = theConfig->getMFPInterpolationBorder() * rings[ringIndex].getEffectiveTemperature();
          double E1_dash = rings[ringIndex].transformEnergyToComovingFrame( particles[iscat].Mom );
          double E2_dash = rings[ringIndex].transformEnergyToComovingFrame( particles[jscat].Mom );
          
          if ( E1_dash < interpolLimit && E2_dash < interpolLimit )
          {
            averagedRate = ( particles[iscat].rate22 + particles[iscat].rate23 + particles[iscat].rate32 +
            particles[jscat].rate22 + particles[jscat].rate23 + particles[jscat].rate32 +
            particles[iscat].rate22v + particles[iscat].rate23v + particles[iscat].rate32v +
            particles[jscat].rate22v + particles[jscat].rate23v + particles[jscat].rate32v ) / ( 2.0 * 2.0 );            
          }
          else
          {
            int tempJetID;
            double tempJetEnergy;
            FLAVOR_TYPE tempJetFlavor;
            
            if ( E1_dash >= interpolLimit && E1_dash >= E2_dash )
            {
              tempJetID = iscat;
              tempJetEnergy = E1_dash;
              tempJetFlavor = F1;
            }
            else
            {
              tempJetID = jscat;
              tempJetEnergy = E2_dash;
              tempJetFlavor = F2;
            }
            
            if ( theConfig->getJetMfpComputationType() == computeMfpInterpolation )
            {
              averagedRate = 1 / ( theMFP->getMeanFreePath( tempJetEnergy, tempJetFlavor, rings[ringIndex].getEffectiveTemperature(),
                                                            rings[ringIndex].getGluonDensity(), rings[ringIndex].getQuarkDensity(), GeV ) );
            }
            else if ( theConfig->getJetMfpComputationType() == computeMfpIteration )
            {
              averagedRate = 1 / iterateMFP( allParticlesList, gluonList, tempJetID, dt, dv, timenow );
            }
            else
            {
              string errMsg = "jetMfpComputationType not valid";
              throw eHIC_error( errMsg );
            }
            
          }
        }
        else
        {
          averagedRate = ( particles[iscat].rate22 + particles[iscat].rate23 + particles[iscat].rate32 +
                         particles[jscat].rate22 + particles[jscat].rate23 + particles[jscat].rate32 +
                         particles[iscat].rate22v + particles[iscat].rate23v + particles[iscat].rate32v +
                         particles[jscat].rate22v + particles[jscat].rate23v + particles[jscat].rate32v ) / ( 2.0 * 2.0 );
        }

        if ( s < 1.1*lambda2 )
        {
          probab2322 = -1.0;
          cs22 = cs23 = 0.0; //1/GeV^2
        }
        else
        {
          Vrel = VelRel( particles[iscat].Mom,particles[jscat].Mom, M1,M2 ); // general relative velocity

          //factor as (alpha_s) is not included in definitions of partcl[jscat].md2g
          md2g_wo_as = ( particles[iscat].md2g + particles[jscat].md2g ) / 2.0;
          md2q_wo_as = ( particles[iscat].md2q + particles[jscat].md2q ) / 2.0;

          //Reconstruction of the temperature for Nf=3 only.
          //double effectiveLRFTempGeVNf3 = sqrt( md2g/8.0*M_PI / (ns_casc::Ncolor + 3.0) ); 
          
          //Reconstruction of the temperature for Nf=3 only.
          double effectiveLRFTempGeVNf0 = sqrt( md2g_wo_as/8.0*M_PI / (ns_casc::Ncolor + 0.0) );
//           cout << effectiveLRFTempGeVNf0/TeffLocalCell<< endl;
          
          
          if ( theConfig->doScattering_22() )
          {
            scatt22_object.setParameter( particles[iscat].Mom, particles[jscat].Mom,
                                         F1, F2, M1, M2, s, Vrel, md2g_wo_as , md2q_wo_as,
                                        theConfig->getKggQQb(), theConfig->getKgQgQ(), theConfig->getKappa_gQgQ(), 
                                        theConfig->isConstantCrossSecGQ(),
                                        theConfig->getConstantCrossSecValueGQ(), theConfig->isIsotropicCrossSecGQ(),
                                        theConfig->getKfactor_light(),Vrel ); // md2g, md2q are debye masses without the factor alpha_s which is multiplied in scattering22.cpp
            
            //---------------------------------//
            string dummyStr;
            switch ( theConfig->getCrossSectionMethod() )
            {
            case csMethod_pQCD:
              cs22 = scatt22_object.getXSection22( initialStateIndex );
              break;
            case csMethod_constCS:
              cs22 = specialCS;
              break;
            case csMethod_constEtaOverS:
              cs22 = specialCS;  
              break;    
            case csMethod_constMFP:
              cs22 = specialCS;    
              break;    
            case csMethod_constMixtureCS:
              cs22 = scatt22_hydro.getMixtureXSection22(theHydroParticleType,vecTypeCommon,F1,F2,1.0,dummyStr);
              break; 
            case csMethod_constMixtureCS_scaledWithLambdaT2:  
              cs22 = scatt22_hydro.getMixtureXSection22(theHydroParticleType,vecTypeCommon,F1,F2,scalingFactor,dummyStr);
              break;  
            case csMethod_variableCS_scaledWithLambdaT2:
              cs22 = specialCS / scalingFactor;
              break;
            default:
              string errMsg = "Error in method of cross section";
              throw eHIC_error( errMsg );
            }
            //---------------------------------//
            
            probab22 = pow( 0.197, 2.0 ) * cs22 * Vrel * dt * scaleForSelectedPairs / ( dv * testpartcl );
//             cout << probab22 << endl;
          }
          else
          {
            cs22 = 0;
            probab22 = 0;
          }


          if ( theConfig->doScattering_23() )
          {
            if ( averagedRate > epsilon )
            {
              lambda_scaled = sqrt( s ) / averagedRate;  //dimensionless
            }
            else
            {
              xsection_gg_gg csObj( s, md2g_wo_as, md2q_wo_as, &theI22, theConfig->getKfactor_light() );
              csgg = csObj.totalCrossSection();
              lambda_scaled = ( dv * rings[ringIndex].getGamma() * testpartcl * sqrt( s ) ) / ( nTotal * csgg * pow( 0.197, 3.0 ) );
            }
            
            betaDistEntry = scatt23_object.setParameter( rings[ringIndex].getAveraged_v(),
                                                         particles[iscat].Mom, particles[jscat].Mom,
                                                         F1, F2, M1, M2, sqrt( s ), 
                        md2g_wo_as / s, lambda_scaled, 
                        theConfig->getK23LightPartons(), theConfig->getK23HeavyQuarks(),
                        theConfig->getKappa23LightPartons(), theConfig->getKappa23HeavyQuarks(),
                        theConfig->I23onlineIntegrationIsSet(),
                        theConfig->get23GluonFormationTimeTyp(), theConfig->getMatrixElement23(), theConfig->isMd2CounterTermInI23(),
                        theConfig->get23FudgeFactorLpm(), nGluons, theConfig->isMatrixElement23_22qt() );
            cs23 = scatt23_object.getXSection23();
            aa.addBoost23( betaDistEntry );

            probab23 = pow( 0.197, 2.0 ) * cs23 * Vrel * dt * scaleForSelectedPairs / ( dv * testpartcl );
          }
          else
          {
            cs23 = 0;
            probab23 = 0;
          }
          
          if ( F1 == gluon )
          {
            probab22 *= scaleFactor;
            probab23 *= scaleFactor;
          }
          if ( F2 == gluon )
          {
            probab22 *= scaleFactor;
            probab23 *= scaleFactor;
          }

          _cell.rates.add( c22, F1, F2, probab22 );
          _cell.rates.add( c23, F1, F2, probab23 );

          if ( cs22 > 0.0 )
          {
            ++_cell.nCollected22;
            _cell.md2g_scaled_22 += md2g_wo_as * coupling::get_constant_coupling() / s;
            _cell.md2q_scaled_22 += md2q_wo_as * coupling::get_constant_coupling() / s;
            _cell.alpha_s_22 +=  coupling::get_constant_coupling();
          }
          if ( cs23 > 0.0 )
          {
            ++_cell.nCollected23;
            _cell.md2g_scaled_23 += md2g_wo_as * coupling::get_constant_coupling() / s;
            _cell.md2q_scaled_23 += md2q_wo_as * coupling::get_constant_coupling() / s;
            _cell.alpha_s_23 +=  coupling::get_constant_coupling();
            _cell.lambdaScaled += lambda_scaled;
          }

          probab2322 = probab22 + probab23;
        }
        //--------------------------------------//

        ++_cell.nCollectedAll2223;
        _cell.sigma_22 += cs22;               //1/GeV^2
        _cell.sigma_23 += cs23;               //1/GeV^2

        if ( probab2322 > 1.0 )
        {
          if ( theConfig->repeatTimesteps() )
          {
            cout << "P2322=" << probab2322 << ">1" << endl;
            again = true;
            cout << "dt (old) = " << dt << endl;
            dt = 0.5 / ( probab2322 / dt );
            cout << "dt (new) = " << dt << endl;
            return;
          }
          else
          {
            cout << "P2322=" << probab2322 << ">1" << endl;             
            probab2322 = 1;
          }
        }

        if ( ran2() < probab2322 )
        {
          double pt_iscat = particles[iscat].Mom.Pt();
          double pt_jscat = particles[jscat].Mom.Pt();
          double pt_nmb;

          if ( ran2() * probab2322 < probab23 )
          {      
            //------------------------------------            
            int jetEventIndex = -1;
            if(theConfig->getAnalyseForHydro() == false)
              {                  
                if ( pt_iscat > aa.getJetTracking_PT() || pt_jscat > aa.getJetTracking_PT() )
                {
                  jetEventIndex = aa.addJetEvent_in( iscat, jscat, -1, c2to3, cs23, _cell.index, lambda_scaled / sqrt( s ) );
                }
              }
              
            //------------------------------------
            int newIndex = scatt23_utility( scatt23_object, _cell, iscat, jscat );

            if(theConfig->getAnalyseForHydro() == false)
            {  
              pt_iscat = particles[iscat].Mom.Pt();
              pt_jscat = particles[jscat].Mom.Pt();
              pt_nmb = particles[newIndex].Mom.Pt();
              if ( jetEventIndex != -1 || pt_iscat > aa.getJetTracking_PT() || pt_jscat > aa.getJetTracking_PT() || pt_nmb > aa.getJetTracking_PT() )
              {
                aa.addJetEvent_out( jetEventIndex, iscat, jscat, newIndex );
              }

            }
            //------------------------------------              
          }
          else
          {
            //------------------------------------            
            int jetEventIndex = -1;
            if(theConfig->getAnalyseForHydro() == false)
              {                  
                if ( pt_iscat > aa.getJetTracking_PT() || pt_jscat > aa.getJetTracking_PT() )
                {
                  jetEventIndex = aa.addJetEvent_in( iscat, jscat, -1, c2to2, cs22, _cell.index, lambda_scaled / sqrt( s ) );
                }
              }

            //------------------------------------
            scatt22_utility( scatt22_object, iscat, jscat, typ );
            //------------------------------------            

            if(theConfig->getAnalyseForHydro() == false)
            {
              pt_iscat = particles[iscat].Mom.Pt();
              pt_jscat = particles[jscat].Mom.Pt();
              if ( jetEventIndex != -1 || pt_iscat > aa.getJetTracking_PT() || pt_jscat > aa.getJetTracking_PT() )
              {
                aa.addJetEvent_out( jetEventIndex, iscat, jscat, -1 );
              }
            }
          }
        }
      }
    }
  }
  else
  {
    for ( iIt =  _cell.particleList.begin(); iIt != _cell.particleList.end(); iIt++ )
    {
      iscat = *iIt;

      list<int>::const_iterator start_j = iIt;
      start_j++;

      for ( jIt = start_j; jIt != _cell.particleList.end(); jIt++ )
      {
        jscat = *jIt;

        F1 = particles[iscat].FLAVOR;
        M1 = particles[iscat].m;
        F2 = particles[jscat].FLAVOR;
        M2 = particles[jscat].m;

        s = (particles[iscat].Mom + particles[jscat].Mom).M2();

        xt = ( particles[iscat].Pos.Perp() +
               particles[jscat].Pos.Perp() ) / 2.0;

        ringIndex = rings.getIndex( xt );
        
        //WORKAROUND: set rings to value of cell
        if(!theConfig->useRings())
        {
          rings.setValuesAsWorkaround(ringIndex, _cell.getEnergyDensity(),_cell.getParticleDensity(),_cell.getGamma(),_cell.getNumberOfParticles(),_cell.getGluonDensity(),_cell.getQuarkDensity());
        }
        
        
        if ( ( F1 <= 2*Particle::max_N_light_flavor && F2 <= 2*Particle::max_N_light_flavor ) && theConfig->getJetMfpComputationType() != computeMfpDefault )
        { 
          double interpolLimit = theConfig->getMFPInterpolationBorder() * rings[ringIndex].getEffectiveTemperature();
          double E1_dash = rings[ringIndex].transformEnergyToComovingFrame( particles[iscat].Mom );
          double E2_dash = rings[ringIndex].transformEnergyToComovingFrame( particles[jscat].Mom );
          
          if ( E1_dash < interpolLimit && E2_dash < interpolLimit )
          {
            averagedRate = ( particles[iscat].rate22 + particles[iscat].rate23 + particles[iscat].rate32 +
            particles[jscat].rate22 + particles[jscat].rate23 + particles[jscat].rate32 +
            particles[iscat].rate22v + particles[iscat].rate23v + particles[iscat].rate32v +
            particles[jscat].rate22v + particles[jscat].rate23v + particles[jscat].rate32v ) / ( 2.0 * 2.0 );            
          }
          else
          {
            int tempJetID;
            double tempJetEnergy;
            FLAVOR_TYPE tempJetFlavor;
            
            if ( E1_dash >= interpolLimit && E1_dash >= E2_dash )
            {
              tempJetID = iscat;
              tempJetEnergy = E1_dash;
              tempJetFlavor = F1;
            }
            else
            {
              tempJetID = jscat;
              tempJetEnergy = E2_dash;
              tempJetFlavor = F2;
            }
            
            if ( theConfig->getJetMfpComputationType() == computeMfpInterpolation )
            {
              averagedRate = 1 / ( theMFP->getMeanFreePath( tempJetEnergy, tempJetFlavor, rings[ringIndex].getEffectiveTemperature(),
                                                            rings[ringIndex].getGluonDensity(), rings[ringIndex].getQuarkDensity(), GeV ) );
            }
            else if ( theConfig->getJetMfpComputationType() == computeMfpIteration )
            {
              averagedRate = 1 / iterateMFP( allParticlesList, gluonList, tempJetID, dt, dv, timenow );
            }
            else
            {
              string errMsg = "jetMfpComputationType not valid";
              throw eHIC_error( errMsg );
            }
            
          }
        }
        else
        {
          averagedRate = ( particles[iscat].rate22 + particles[iscat].rate23 + particles[iscat].rate32 +
                         particles[jscat].rate22 + particles[jscat].rate23 + particles[jscat].rate32 +
                         particles[iscat].rate22v + particles[iscat].rate23v + particles[iscat].rate32v +
                         particles[jscat].rate22v + particles[jscat].rate23v + particles[jscat].rate32v ) / ( 2.0 * 2.0 );
        }

        if ( s < 1.1*lambda2 )
        {
          probab2322 = -1.0;
          cs22 = cs23 = 0.0; //1/GeV^2
        }
        else
        {
          Vrel = VelRel( particles[iscat].Mom,particles[jscat].Mom, M1,M2 ); // general relative velocity

          //factor as (alpha_s) is not included in definitions of partcl[jscat].md2g
          md2g_wo_as = ( particles[iscat].md2g + particles[jscat].md2g ) / 2.0;
          md2q_wo_as = ( particles[iscat].md2q + particles[jscat].md2q ) / 2.0;
          
          //Reconstruction of the temperature for Nf=3 only.
          //double effectiveLRFTempGeVNf3 = sqrt( md2g/8.0*M_PI / (ns_casc::Ncolor + 3.0) ); 
          
          //Reconstruction of the temperature for Nf=3 only.
          double effectiveLRFTempGeVNf0 = sqrt( md2g_wo_as/8.0*M_PI / (ns_casc::Ncolor + 0.0) );
//           cout << effectiveLRFTempGeVNf0/TeffLocalCell<< endl;
          
          if ( theConfig->doScattering_22() )
          {
            scatt22_object.setParameter( particles[iscat].Mom, particles[jscat].Mom,
                                         F1, F2, M1, M2, s, Vrel, md2g_wo_as , md2q_wo_as,
                                        theConfig->getKggQQb(), theConfig->getKgQgQ(), theConfig->getKappa_gQgQ(), 
                                        theConfig->isConstantCrossSecGQ(),
                                        theConfig->getConstantCrossSecValueGQ(), theConfig->isIsotropicCrossSecGQ(),
                                        theConfig->getKfactor_light() ); // md2g, md2q are debye masses without the factor alpha_s which is multiplied in scattering22.cpp
            
            //---------------------------------//
            string dummyStr;
            switch ( theConfig->getCrossSectionMethod() )
            {
            case csMethod_pQCD:
              cs22 = scatt22_object.getXSection22( initialStateIndex );
              break;
            case csMethod_constCS:
              cs22 = specialCS;
              break;
            case csMethod_constEtaOverS:
              cs22 = specialCS;  
              break;    
            case csMethod_constMFP:
              cs22 = specialCS;    
              break;    
            case csMethod_constMixtureCS:
              cs22 = scatt22_hydro.getMixtureXSection22(theHydroParticleType,vecTypeCommon,F1,F2,1.0,dummyStr);
              break; 
            case csMethod_constMixtureCS_scaledWithLambdaT2:  
              cs22 = scatt22_hydro.getMixtureXSection22(theHydroParticleType,vecTypeCommon,F1,F2,scalingFactor,dummyStr);
              break;  
            case csMethod_variableCS_scaledWithLambdaT2:
              cs22 = specialCS / scalingFactor;
              break;
            default:
              string errMsg = "Error in method of cross section";
              throw eHIC_error( errMsg );
            }
            //---------------------------------//
            
            probab22 = pow( 0.197, 2.0 ) * cs22 * Vrel * dt / ( dv * testpartcl );
          }
          else
          {
            cs22 = 0;
            probab22 = 0;
          }


          if ( theConfig->doScattering_23() )
          {
            if ( averagedRate > epsilon )
            {
              lambda_scaled = sqrt( s ) / averagedRate;  //dimensionless
            }
            else
            {
              xsection_gg_gg csObj( s, md2g_wo_as, md2q_wo_as, &theI22, theConfig->getKfactor_light() );
              csgg = csObj.totalCrossSection();
              lambda_scaled = ( dv * rings[ringIndex].getGamma() * testpartcl * sqrt( s ) ) / ( nTotal * csgg * pow( 0.197, 3.0 ) );
            }
            
            betaDistEntry = scatt23_object.setParameter( rings[ringIndex].getAveraged_v(),
                                                         particles[iscat].Mom, particles[jscat].Mom,
                                                         F1, F2, M1, M2, sqrt( s ), 
                        md2g_wo_as / s, lambda_scaled, 
                        theConfig->getK23LightPartons(), theConfig->getK23HeavyQuarks(),
                        theConfig->getKappa23LightPartons(), theConfig->getKappa23HeavyQuarks(),
                        theConfig->I23onlineIntegrationIsSet(),
                        theConfig->get23GluonFormationTimeTyp(), theConfig->getMatrixElement23(), theConfig->isMd2CounterTermInI23(),
                        theConfig->get23FudgeFactorLpm(), nGluons, theConfig->isMatrixElement23_22qt() );
            cs23 = scatt23_object.getXSection23();
            aa.addBoost23( betaDistEntry );

            probab23 = pow( 0.197, 2.0 ) * cs23 * Vrel * dt / ( dv * testpartcl );
          }
          else
          {
            cs23 = 0;
            probab23 = 0;
          }
          
          if ( F1 == gluon )
          {
            probab22 *= scaleFactor;
            probab23 *= scaleFactor;
          }
          if ( F2 == gluon )
          {
            probab22 *= scaleFactor;
            probab23 *= scaleFactor;
          }

          _cell.rates.add( c22, F1, F2, probab22 );
          _cell.rates.add( c23, F1, F2, probab23 );

          if ( cs22 > 0.0 )
          {
            ++_cell.nCollected22;
            _cell.md2g_scaled_22 += md2g_wo_as * coupling::get_constant_coupling() / s;
            _cell.md2q_scaled_22 += md2q_wo_as * coupling::get_constant_coupling() / s;
            _cell.alpha_s_22 +=  coupling::get_constant_coupling();
          }
          if ( cs23 > 0.0 )
          {
            ++_cell.nCollected23;
            _cell.md2g_scaled_23 += md2g_wo_as * coupling::get_constant_coupling() / s;
            _cell.md2q_scaled_23 += md2q_wo_as * coupling::get_constant_coupling() / s;
            _cell.alpha_s_23 +=  coupling::get_constant_coupling();
            _cell.lambdaScaled += lambda_scaled;
          }

          probab2322 = probab22 + probab23;
        }
        //--------------------------------------//

        ++_cell.nCollectedAll2223;
        _cell.sigma_22 += cs22;               //1/GeV^2
        _cell.sigma_23 += cs23;               //1/GeV^2

        if ( probab2322 > 1.0 )
        {
          if ( theConfig->repeatTimesteps() )
          {
            cout << "P2322=" << probab2322 << ">1" << endl;
            again = true;
            cout << "dt (old) = " << dt << endl;
            dt = 0.5 / ( probab2322 / dt );
            cout << "dt (new) = " << dt << endl;
            return;
          }
          else
          {
            cout << "P2322=" << probab2322 << ">1" << endl;             
            probab2322 = 1;
          }
        }

        if ( ran2() < probab2322 )
        {
          double pt_iscat = particles[iscat].Mom.Pt();
          double pt_jscat = particles[jscat].Mom.Pt();
          double pt_nmb;

          if ( ran2() * probab2322 < probab23 )
          {      
            //------------------------------------            
            int jetEventIndex = -1;
            if(theConfig->getAnalyseForHydro() == false)
              {                  
                if ( pt_iscat > aa.getJetTracking_PT() || pt_jscat > aa.getJetTracking_PT() )
                {
                  jetEventIndex = aa.addJetEvent_in( iscat, jscat, -1, c2to3, cs23, _cell.index, lambda_scaled / sqrt( s ) );
                }
              }
              
            //------------------------------------
            int newIndex = scatt23_utility( scatt23_object, _cell, iscat, jscat );

            if(theConfig->getAnalyseForHydro() == false)
            {  
              pt_iscat = particles[iscat].Mom.Pt();
              pt_jscat = particles[jscat].Mom.Pt();
              pt_nmb = particles[newIndex].Mom.Pt();
              if ( jetEventIndex != -1 || pt_iscat > aa.getJetTracking_PT() || pt_jscat > aa.getJetTracking_PT() || pt_nmb > aa.getJetTracking_PT() )
              {
                aa.addJetEvent_out( jetEventIndex, iscat, jscat, newIndex );
              }

            }
            //------------------------------------              
          }
          else
          {
            //------------------------------------            
            int jetEventIndex = -1;
            if(theConfig->getAnalyseForHydro() == false)
              {                  
                if ( pt_iscat > aa.getJetTracking_PT() || pt_jscat > aa.getJetTracking_PT() )
                {
                  jetEventIndex = aa.addJetEvent_in( iscat, jscat, -1, c2to2, cs22, _cell.index, lambda_scaled / sqrt( s ) );
                }
              }

            //------------------------------------
            scatt22_utility( scatt22_object, iscat, jscat, typ );
            //------------------------------------            

            if(theConfig->getAnalyseForHydro() == false)
            {
              pt_iscat = particles[iscat].Mom.Pt();
              pt_jscat = particles[jscat].Mom.Pt();
              if ( jetEventIndex != -1 || pt_iscat > aa.getJetTracking_PT() || pt_jscat > aa.getJetTracking_PT() )
              {
                aa.addJetEvent_out( jetEventIndex, iscat, jscat, -1 );
              }
            }
          }
        }
      }
    }
  }


  for( iIt = _cell.particleList.begin(); iIt != _cell.particleList.end(); iIt++ )
  {
    particles[*iIt].Propagate( timenext );
  }
  
  
  
  
}


/**
 * When a 3->2 interaction takes place, this routine handles the setting of new momenta, flavors, the removal of annihilated particles etc.
 * @param[in] scatt32_obj A reference to the scattering32 object currently associated with the considered particle triplet
 * @param _cellMembers A reference to the particle list of the current cell
 * @param[in] _allParticlesList List of all particles in the cell (TODO redundant information)
 * @param[in] _gluonList List of all gluons in the cell
 * @param[in] iscat ID of particle 1
 * @param[in] jscat ID of particle 2
 * @param[in] kscat ID of particle 3
 * @param[out] n32 Number of 3->2 processes, n32 is increased in this routine
 */
void heavyIonCollision::scatt32_utility( scattering32& scatt32_obj, std::list< int >& _cellMembers, std::vector< int >& _allParticlesList, std::vector< int >& _gluonList, const int iscat, const int jscat, const int kscat, int& n32 )
{
  double Tmax, TT, u, phi;
  FLAVOR_TYPE F1, F2, F3;
  int typ;

  F1 = particles[iscat].FLAVOR;
  F2 = particles[jscat].FLAVOR;
  F3 = particles[kscat].FLAVOR;
  
  // these routines are not written for heavy flavors
  if( F1 > 2*Particle::max_N_light_flavor || F2 > 2*Particle::max_N_light_flavor || F3 > 2*Particle::max_N_light_flavor )
  {
    std::string errMsg = "Heavy flavor in scatt32_utility.";
    throw eHIC_error( errMsg );
  }

  ncoll++;
  ncoll32++;
  n32++;
  particles[iscat].coll_id = particles[jscat].coll_id = particles[kscat].coll_id = ncoll;

  if ( particles[iscat].Pos.T() > particles[jscat].Pos.T() )
    Tmax = particles[iscat].Pos.T();
  else
    Tmax = particles[jscat].Pos.T();
  if ( particles[kscat].Pos.T() > Tmax )
    Tmax = particles[kscat].Pos.T();
  TT = ( timenext - Tmax ) * ran2() + Tmax;

  particles[iscat].Propagate( TT );
  particles[jscat].Propagate( TT );
  particles[kscat].Propagate( TT );

  VectorEPxPyPz P1_new, P2_new;
  int absorbedGluon = scatt32_obj.getMomenta32( u, phi, typ, F1, F2 );
  scatt32_obj.setNewMomenta32( P1_new, P2_new, u, phi );
  P1_new(0) = sqrt( P1_new.vec2() );
  P2_new(0) = sqrt( P2_new.vec2() );

  int i1, i2, i3;
  switch ( absorbedGluon )
  {
    case 1:
    {
      i1 = jscat;
      i2 = kscat;
      i3 = iscat;
    } break;
    case 2:
    {
      i1 = iscat;
      i2 = kscat;
      i3 = jscat;
    } break;
    case 3:
    {
      i1 = iscat;
      i2 = jscat;
      i3 = kscat;
    } break;
    
  }

  particles[i1].FLAVOR = F1;
  particles[i1].Mom = P1_new;
  particles[i2].FLAVOR = F2;
  particles[i2].Mom = P2_new;
  
  if ( particles[i3].FLAVOR != gluon )
  {
    cout << "xx " << iscat << "  " << jscat << "  " << kscat;
    cout << " : " << i1 << "  " << i2 << "  " << i3 << endl;
    
    string errMsg = "error in scatt32_utility: dead particle i3 is not a gluon";
    throw eHIC_error( errMsg );
  }
  
  //mark absorbed gluon for removal from global particle list
  deadParticleList.push_back( i3 );
  particles[i3].dead = true;
  
  //erase absorbed gluon from the local list of gluons in this cell
  removeElementFromVector<int>( _gluonList, i3 );
  //erase absorbed gluon from the local list of all particles in this cell
  removeElementFromVector<int>( _allParticlesList, i3 );
  
  //erase absorbed gluon from the list of all particles in this cell
  _cellMembers.remove( i3 );
  
  if ( theConfig->doOutput_offlineReconstructionData() )
  {
    tPointerToOfflineData _tmpPtr( new offlineDataInteraction32( i1, i2, i3, TT, 
                                                                 particles[i1].Mom,
                                                                 particles[i2].Mom,
                                                                 F1, F2 ) );
    
    offlineInterface.registerOfflineDataForTemporaryStorage( event_interaction32, _tmpPtr );
  }
}




/**
 * When a 2->3 interaction takes place, this routine handles the setting of new momenta, flavors, the adding of the radiated particle etc.
 * @param[in] scatt23_obj A reference to the scattering23 object currently associated with the considered particle pair
 * @param _cell A reference to the current cell
 * @param[in] iscat ID of particle 1
 * @param[in] jscat ID of particle 2
 */
int heavyIonCollision::scatt23_utility( scattering23& scatt23_obj, cellContainer& _cell, int iscat, const int jscat )
{
  FLAVOR_TYPE F1, F2;
  //  int nc, IXY;
  //  double cx, cy, cz, dz;
  double Tmax, TT, pt1, pt3, y, phi, pz1;
  int typ;

  int etaIndex = etaBins.getIndex( particles[iscat].eta );
  double leftZ = timenow * tanh( etaBins[etaIndex].left );
  double deltaZ = timenow * tanh( etaBins[etaIndex].right ) - leftZ;
  double leftY = _cell.corner.y_min + randomShiftY;
  double leftX = _cell.corner.x_min + randomShiftX;
  
  F1 = particles[iscat].FLAVOR;
  F2 = particles[jscat].FLAVOR;

  ncoll++;
  ncoll23++;
  particles[iscat].coll_id = particles[jscat].coll_id = ncoll;

  Tmax = std::max(particles[iscat].Pos.T(),particles[jscat].Pos.T() );
  TT = ( timenext - Tmax ) * ran2() + Tmax;

  particles[iscat].Propagate( TT );
  particles[jscat].Propagate( TT );

  nGet23Errors += scatt23_obj.getMomenta23( pt1, pt3, y, phi, pz1, typ, F1, F2 );

  VectorEPxPyPz P1new, P2new, P3new;
  scatt23_obj.setNewMomenta23( P1new, P2new, P3new, 
                               particles[iscat].Pos, particles[jscat].Pos, 
                               pt1, pt3, y, phi, pz1 );
  
  particles[iscat].FLAVOR = F1;
  particles[iscat].Mom = P1new;

  particles[jscat].FLAVOR = F2;
  particles[jscat].Mom = P2new;

  Particle tempParticle;
  tempParticle.FLAVOR = gluon;
  tempParticle.Mom = P3new;
  tempParticle.Pos = VectorTXYZ(TT,leftX + dx * ran2(),leftY + dy * ran2(),leftZ + deltaZ * ran2());
  tempParticle.m = 0.0;
  tempParticle.md2g = ( particles[iscat].md2g + particles[jscat].md2g ) / 2.0;
  tempParticle.md2q = ( particles[iscat].md2q + particles[jscat].md2q ) / 2.0;
  tempParticle.edge = -1;
  tempParticle.dead = false;
  tempParticle.free = false;
  tempParticle.cell_id = -1;
  tempParticle.coll_id = ncoll;
  tempParticle.collisionTime = infinity;
  tempParticle.collisionPartner = -1;
  tempParticle.rate23 = particles[iscat].rate23;    //GeV
  tempParticle.rate32 = particles[iscat].rate32;    //GeV
  tempParticle.rate22 = particles[iscat].rate22;    //GeV
  tempParticle.rate23v = particles[iscat].rate23v;  //GeV
  tempParticle.rate32v = particles[iscat].rate32v;  //GeV
  tempParticle.rate22v = particles[iscat].rate22v;  //GeV
  tempParticle.as22 = tempParticle.as23 = 0;
  tempParticle.cs22 = tempParticle.cs23 = 0;        //1/GeV^2
  tempParticle.md2g_scaled_22 = 0;
  tempParticle.md2q_scaled_22 = 0;
  tempParticle.md2g_scaled_23 = 0;
  tempParticle.md2q_scaled_23 = 0;
  tempParticle.lambda_scaled = 0;
  tempParticle.unique_id = Particle::unique_id_counter;
  ++Particle::unique_id_counter;
  
  particles.push_back( tempParticle );
  int newIndex = particles.size() - 1;

  if ( theConfig->doOutput_offlineReconstructionData() )
  {
    tPointerToOfflineData _tmpPtr( new offlineDataInteraction23( iscat, jscat, newIndex, TT, 
                                                                 particles[iscat].Mom,
                                                                 particles[jscat].Mom,
                                                                 particles[newIndex].Pos,
                                                                 particles[newIndex].Mom,
                                                                 F1, F2, gluon ) );
    offlineInterface.registerOfflineDataForTemporaryStorage( event_interaction23, _tmpPtr );
  }

  particles[newIndex].Propagate( timenext );
  
  return newIndex;
}



/**
 * When a 2->2 interaction takes place, this routine handles the setting of new momenta, flavors, etc.
 * @param[in] scatt22_obj A reference to the scattering22 object currently associated with the considered particle pair
 * @param[in] iscat ID of particle 1
 * @param[in] jscat ID of particle 2
 * @param[out] typ Type of the 2->2 interaction (221 = gg -> gg, 222 = gg -> qqbar, etc.)
 */
void heavyIonCollision::scatt22_utility( scattering22& scatt22_obj, const int iscat, const int jscat, int& typ )
{
  FLAVOR_TYPE F1, F2;
  double Tmax, TT;
  double M1, M2;
  double t_hat,s, vrel;

  s=scatt22_obj.getS();
  vrel=scatt22_obj.getVrel();
  
  F1 = particles[iscat].FLAVOR;
  M1 = particles[iscat].m;

  F2 = particles[jscat].FLAVOR;
  M2 = particles[jscat].m;
  
  if(particles[iscat].Pos.T() > timenext)
    cout << "Should not happen, scatt22: " << particles[iscat].Pos.T() << endl;
  if(particles[jscat].Pos.T() > timenext)
    cout << "Should not happen, scatt22:" <<  particles[iscat].Pos.T() << endl;  
  
  //  double s_initial = (particles[iscat].Mom + particles[jscat].Mom).M2();

  
  ncoll++;
  ncoll22++;
  particles[iscat].coll_id = particles[jscat].coll_id = ncoll;

  Tmax = std::max(particles[iscat].Pos.T(),particles[jscat].Pos.T() );
  TT = ( timenext - Tmax ) * ran2() + Tmax;

  particles[iscat].Propagate( TT );
  particles[jscat].Propagate( TT );

  if( theConfig->isIsotropicCrossSection() == false)
  {  
    // determine type of scattering, momentum transfer t_hat, new flavor and new masses
    scatt22_obj.getMomentaAndMasses22(F1, F2, M1, M2, t_hat, typ);
  }
  else
  {
    // determine type of scattering, momentum transfer t_hat, new flavor and new masses for isotropic processes
    scatt22_obj.getMomentaAndMasses22_isotropic( F1, F2, M1, M2, t_hat );
  }
  
  //Compute sin^2 theta at midrapidity
  double eta_i = particles[iscat].Mom.Rapidity();
  double eta_j = particles[jscat].Mom.Rapidity();
  double eta = std::min(std::abs(eta_i),std::abs(eta_j));
  
  if ( fabs( eta ) >= -0.5 && fabs( eta ) <= 0.5 )
  {
    sin2ThetaVrel += vrel*(- 4.*t_hat/s - 4.*pow(t_hat/s,2.0) );
    sin2Theta += (- 4.*t_hat/s - 4.*pow(t_hat/s,2.0) );
    vrelSum += vrel;
    ncoll22MIDRAP +=1;
  }
  // translate momemtum transfer t_hat into actual momenta of outgoing particles
  scatt22_obj.setNewMomenta22( particles[iscat].Mom, particles[jscat].Mom,
                               particles[iscat].Pos, particles[jscat].Pos,
                               t_hat );
  
  // double s_final = (particles[iscat].Mom + particles[jscat].Mom).M2();
  // if (!FPT_COMP_E(s_initial,s_final))
  // {
  //   std::string errMsg = "++++++++++++   2->2 violation of energy conservation   +++++++++++++"; 
  //   throw eHIC_error( errMsg ); 
  // }
  
  particles[iscat].FLAVOR = F1;
  particles[iscat].m = M1;

  particles[jscat].FLAVOR = F2;
  particles[jscat].m = M2;
  
  if ( theConfig->doOutput_offlineReconstructionData() )
  {
    tConstPointerToOfflineData _tmpPtr( new offlineDataInteraction22( iscat, jscat, TT, 
                                                                      particles[iscat].Mom,
                                                                      particles[jscat].Mom,
                                                                      F1, F2 ) );
    offlineInterface.registerOfflineDataForTemporaryStorage( event_interaction22, _tmpPtr );
  }
}


/**
 * @param[in] _particleList The list of particles for which the geometric collisions should be performed (i.e. one of the edge cells)
 * @param[in] firstScatteredParticle The ID of the particle which is the first to collide geometrically within the given _particleList
 */
void heavyIonCollision::doGeometricCollisions( std::list< int >& _particleList, const int firstScatteredParticle )
{
  FLAVOR_TYPE F1, F2;
  int jscat, typ;
  double M1, M2, md2g_wo_as, md2q_wo_as;
  double s, t_hat, time, ct_i, ct_j;
  //  double t;
  
  time = timenow;
  if ( firstScatteredParticle != -1 )
  {
    int iscat = firstScatteredParticle;
    while ( iscat != -1 && particles[iscat].collisionTime < timenext )
    {
      if ( particles[iscat].collisionTime < time )
      {
        int collP = particles[iscat].collisionPartner;
        cout << "tttt " << iscat << "  " << particles[iscat].unique_id << "  " << particles[iscat].collisionTime << "  " << time << "  " << particles[iscat].coll_id << endl;
        cout << "tttt " << particles[iscat].free << "   " << particles[iscat].dead << "  " << particles[iscat].edge << "  " << particles[iscat].FLAVOR << endl;
        cout << "tttt " << particles[iscat].Pos << endl;
        cout << "tttt " << particles[iscat].Mom << endl;
        cout << "tttt p  " << collP << "  " << particles[collP].unique_id << "  " << particles[collP].collisionTime << "  " << time << "  " << particles[collP].coll_id << endl;
        cout << "tttt p  " << particles[collP].free << "   " << particles[collP].dead << "  " << particles[collP].edge << "  " << particles[collP].FLAVOR << endl;
        cout << "tttt p  " << particles[collP].Pos << endl;
        cout << "tttt p  " << particles[collP].Mom << endl;
        cout << "tttt particle list: " << _particleList.size() << endl;
        string errMsg = "backward propagation in geometricCollisionForSpecifiedParticle";
        throw eHIC_error( errMsg );
      }
      
      jscat = particles[iscat].collisionPartner;
      
      ncoll++;
      ncolle++;
      particles[iscat].coll_id = particles[jscat].coll_id = ncoll;
      
      prepareGeometricCollision( iscat, jscat, M1, M2, F1, F2, s, ct_i, ct_j );
      
      if ( ct_i < ct_j )
      {
        time = ct_i;
      }
      else
      {
        time = ct_j;
      }
      
      if ( particles[iscat].md2g < particles[jscat].md2g )
      {
        md2g_wo_as = particles[iscat].md2g;
      }
      else
      {
        md2g_wo_as = particles[jscat].md2g;
      }
      
      if ( particles[iscat].md2q < particles[jscat].md2q )
      {
        md2q_wo_as = particles[iscat].md2q;
      }
      else
      {
        md2q_wo_as = particles[jscat].md2q;
      }
            
      scattering22 scatt22_object( &theI22 );
      scatt22_object.setParameter( particles[iscat].Mom, particles[jscat].Mom,
                                   F1, F2, M1, M2, s, 1.0, md2g_wo_as , md2q_wo_as,
                                   theConfig->getKggQQb(), theConfig->getKgQgQ(), 
                                   theConfig->getKappa_gQgQ(), 
                                   theConfig->isConstantCrossSecGQ(),
                                   theConfig->getConstantCrossSecValueGQ(), 
                                   theConfig->isIsotropicCrossSecGQ(),
                                   theConfig->getKfactor_light() ); // md2g, md2q are debye masses without the factor alpha_s which is multiplied in scattering22.cpp
      // determine type of scattering and momentum transfer t_hat for elastic collision
      scatt22_object.getMomentaElastic( t_hat, typ );

      VectorEPxPyPz P1new = particles[iscat].Mom;
      VectorEPxPyPz P2new = particles[jscat].Mom;

      scatt22_object.setNewMomenta22( P1new, P2new,
                                      particles[iscat].Pos, particles[jscat].Pos,
                                      t_hat );

      updateAfterGeometricCollision( iscat, jscat, P1new, P2new, ct_i, ct_j );

      //---------- register data for offline reconstruction ----------
      if ( theConfig->doOutput_offlineReconstructionData() )
      {
        offlineInterface.submitOfflineDataForOutput( event_interactionElastic );
        offlineDataInteractionElastic collectInteractionElastic( iscat, jscat, ct_i, ct_j,
                                                                 particles[iscat].Mom,
                                                                 particles[jscat].Mom );
        offlineInterface.submitOfflineDataForOutput( &collectInteractionElastic );
      }
      //---------- register data for offline reconstruction ----------
      int nextParticle = updateGeometricCollisionTimesB( _particleList, iscat );
      iscat = nextParticle;

    }
  }
}



/**
 * @param[in] _particleList The list of particles for which the geometric collisions should be computed (i.e. one of the edge cells)
 * @param[out] _nextCollTime The time of the next collision
 * @return ID of the particle which is next (first) to collide geometrically
 */
int heavyIonCollision::updateGeometricCollisionTimesA( list<int>& _particleList, double& _nextCollTime )
{
  int mi, mj;
  double ot_ij, min;

  int indexFirstCollisionParticle = -1;
  particles[ _particleList.front()].collisionPartner = -1;
  min = particles[ _particleList.front()].collisionTime = infinity;

  list<int>::const_iterator iIt;
  list<int>::const_iterator jIt;

  for ( iIt =  _particleList.begin(); iIt != _particleList.end(); iIt++ )
  {
    mi = *iIt;
    particles[mi].collisionPartner = -1;
    particles[mi].collisionTime = infinity;

    list<int>::const_iterator start_j = iIt;
    start_j++;

    for ( jIt = start_j; jIt != _particleList.end(); jIt++ )
    {
      mj = *jIt;
      ot_ij = getGeometricCollisionTime( mi, mj );
      if ( ot_ij < particles[mi].collisionTime )
      {
        particles[mi].collisionTime = ot_ij;
        particles[mi].collisionPartner = mj;
      }
    }

    if ( particles[mi].collisionTime < min )
    {
      min = particles[mi].collisionTime;
      indexFirstCollisionParticle = mi;
    }
  }
  _nextCollTime = min;
  return indexFirstCollisionParticle;
}



/**
 * @param[in] _particleList The list of particles for which the geometric collisions should be computed (i.e. one of the edge cells)
 * @param[in] lastUpdatedParticle The ID of the particle for which the previous geometric collision had been carried out
 * @param[out] _nextCollTime The time of the next collision
 * @return ID of the particle which is next to collide geometrically
 */
int heavyIonCollision::updateGeometricCollisionTimesB( std::list<int>& _particleList, int lastUpdatedParticle, double& _nextCollTime )
{
  int mi, mj;
  double ot_ij, ot_ij1, ot_ij2, min;
  
  min = infinity;
  int indexNextCollisionParticle = -1;
  
  const int updated1 = lastUpdatedParticle;
  const int updated2 = particles[ updated1 ].collisionPartner;

  _particleList.remove( updated1 );
  _particleList.remove( updated2 );

  list<int>::const_iterator iIt;
  list<int>::const_iterator jIt;

  for ( iIt = _particleList.begin(); iIt != _particleList.end(); iIt++ )
  {
    mi = *iIt;
    if (( particles[mi].collisionPartner == updated1 ) || ( particles[mi].collisionPartner == updated2 ) )
    {
      particles[mi].collisionPartner = -1;
      particles[mi].collisionTime = infinity;

      list<int>::const_iterator start_j = iIt;
      start_j++;

      for ( jIt = start_j; jIt != _particleList.end(); jIt++ )
      {
        mj = *jIt;

        ot_ij = getGeometricCollisionTime( mi, mj );

        if ( ot_ij < particles[mi].collisionTime )
        {
          particles[mi].collisionTime = ot_ij;
          particles[mi].collisionPartner = mj;
        }
      }
      
      ot_ij = getGeometricCollisionTime( mi, updated1 );
      if ( ot_ij < particles[mi].collisionTime )
      {
        particles[mi].collisionTime = ot_ij;
        particles[mi].collisionPartner = updated1;
      }
      
      ot_ij = getGeometricCollisionTime( mi, updated2 );
      if ( ot_ij < particles[mi].collisionTime )
      {
        particles[mi].collisionTime = ot_ij;
        particles[mi].collisionPartner = updated2;
      }     
    }
    else
    {
      ot_ij1 = getGeometricCollisionTime( mi, updated1 );
      ot_ij2 = getGeometricCollisionTime( mi, updated2 );

      if ( ot_ij1 < particles[mi].collisionTime )
      {
        particles[mi].collisionTime = ot_ij1;
        particles[mi].collisionPartner = updated1;
      }
      if ( ot_ij2 < particles[mi].collisionTime )
      {
        particles[mi].collisionTime = ot_ij2;
        particles[mi].collisionPartner = updated2;
      }
    }

    if ( particles[mi].collisionTime < min )
    {
      min = particles[mi].collisionTime;
      indexNextCollisionParticle = mi;
    }
  }

  particles[updated1].collisionTime = particles[updated2].collisionTime = infinity;
  particles[updated1].collisionPartner = particles[updated2].collisionPartner = -1;
  _particleList.push_back( updated1 );
  _particleList.push_back( updated2 );

  _nextCollTime = min;
  return indexNextCollisionParticle;
}




/**
 * This routines takes care that particles that have been absorbed in 3->2 processes are removed properly, handling cell list,
 * edge lists etc.
 * @param aa A reference to the analysis object that handles output, data collection etc.
 */
void heavyIonCollision::removeDeadParticles( analysis& _aa )
{
  double pt_new;
  int lastIndex = -1;
  deadParticleList.sort();
  while ( particles.back().dead )
  {
    if ( particles.back().edge >= 0 )
    {
      std::string errMsg = "error in removeDeadParticles, removed edge particle";
      throw eHIC_error( errMsg );
    }
    
    int _CellID_deadParticle = particles.back().cell_id;
    if ( _CellID_deadParticle >= 0 )
    {
      cells[_CellID_deadParticle].particleList.remove( particles.size() - 1 );
    }
    
    particles.pop_back();
    deadParticleList.pop_back();
  }

  list<int>::iterator it;
  for ( it = deadParticleList.begin(); it != deadParticleList.end(); it++ )
  {
    lastIndex = particles.size() - 1;
    pt_new = particles.back().Mom.Pt();
    if ( pt_new > _aa.getJetTracking_PT() )
    {
      _aa.exchangeJetID( lastIndex, ( *it ) );
    }

    //---------- register data for offline reconstruction ----------
    if ( theConfig->doOutput_offlineReconstructionData() )
    {
      offlineDataParticleIdSwap collectParticleIdSwap( *it, lastIndex );
      offlineInterface.submitOfflineDataForOutput( &collectParticleIdSwap );
      offlineInterface.submitOfflineDataForOutput( event_particleIdSwap );
    }
    //---------- register data for offline reconstruction ----------

    int _cellID_deadParticle = particles[(*it)].cell_id;
    int _cellID_swappedParticle = particles.back().cell_id;
    if ( _cellID_deadParticle >= 0 )
    {
      cells[_cellID_deadParticle].particleList.remove( *it );
    }
    if ( _cellID_swappedParticle >= 0 )
    {
      cells[_cellID_swappedParticle].particleList.remove( lastIndex );
      cells[_cellID_swappedParticle].particleList.push_back( (*it) );    
    }

    particles[( *it )] = particles.back();
    if ( particles.back().edge >= 0 )
    {
      int swappedParticleEdgeID = particles.back().edge;
      edgeCell[swappedParticleEdgeID].remove( lastIndex );
//       edgeCell[swappedParticleEdgeID].remove( *it ); // just to be sure it's not added twice
      edgeCell[swappedParticleEdgeID].push_back( *it );

      list<int>::iterator jIt;
      for ( jIt = edgeCell[swappedParticleEdgeID].begin(); jIt != edgeCell[swappedParticleEdgeID].end() ; jIt++ )
      {
        if ( particles[*jIt].collisionPartner == lastIndex )
        {
          particles[*jIt].collisionPartner = *it;
        }
      }
    }
    particles.pop_back();

    while ( particles.back().dead )
    {
      if ( particles.back().edge >= 0 )
      {
        std::string errMsg = "error in removeDeadParticles, removed edge particle";
        throw eHIC_error( errMsg );
      }
      
      int _CellID_deadParticle = particles.back().cell_id;
      if ( _CellID_deadParticle >= 0 )
      {
        cells[_CellID_deadParticle].particleList.remove( particles.size() - 1 );
      }
      
      particles.pop_back();
      deadParticleList.pop_back();
    }
  }

}

/**
 * @param[out] _dt The size \Delta t of the first time step
 * @param[out] _min The minimum time of all initial particles (according to their formation time)
 * @return The time t of the first time step 
 */
double heavyIonCollision::getFirstTimestepMinimum( double& _dt, double& _min)
{
  double PartOfParticlesFormed = 0.01;
  
  double min = infinity;       // fm/c
  double max = 0;            // fm/c
  int tempN = 0;
  double _timenext = infinity;

  for ( unsigned int i = 0; i < particles.size(); i++ )
  {
    if(particles[i].Pos.T() > 0)
    {
      if ( particles[i].Pos.T() < min )
      {
        min = particles[i].Pos.T();
      }
      if ( particles[i].Pos.T() > max )
      {
        max = particles[i].Pos.T();
      }
    }
  }
  cout << "minimal time of particles= " << min << "\t" << "maximum time of particles= " << max << endl;

  int upperLimit = etaBins.getNinEtaBin() + etaBins.getNinEtaBin() / 10;       //  10% surplus
  _timenext = min;
  _dt = max - min;
  _dt = _dt / 100.0;
  do
  {
    _timenext += _dt;
    tempN = 0;
    for ( unsigned int j = 0; j < particles.size(); j++ )
    {
      if ( particles[j].Pos.T() < _timenext )
      {
        ++tempN;
      }
    }
  }while(tempN<particles.size()*PartOfParticlesFormed);
  
  cout << "Time when " << PartOfParticlesFormed*100 <<"% of particles formed: " << _timenext << " fm/c \t and dt=" << _dt << " fm/c" << endl;
 
  _min = min;
  return _timenext;
}



/**
 * @param[out] _dt The size \Delta t of the first time step
 * @param[out] _min The minimum time of all initial particles (according to their formation time)
 * @return The time t of the first time step 
 */
double heavyIonCollision::getFirstTimestep( double& _dt, double& _min)
{
  double min = infinity;       // fm/c
  double max = 0;            // fm/c
  int tempN = 0;
  double _timenext = infinity;

  for ( unsigned int i = 0; i < particles.size(); i++ )
  {
    if(particles[i].Pos.T() > 0)
    {
      if ( particles[i].Pos.T() < min )
      {
        min = particles[i].Pos.T();
      }
      if ( particles[i].Pos.T() > max )
      {
        max = particles[i].Pos.T();
      }
    }
  }
  cout << "minimal time of particle= " << min << "\t" << "maximum time of particle== " << max << endl;

  int upperLimit = etaBins.getNinEtaBin() + etaBins.getNinEtaBin() / 10;       //  10% surplus
  _timenext = max;
  _dt = max - min;

  cout << _timenext << '\t' << _dt << endl;
  
  if ( particles.size() > upperLimit )
  {
    tempN = particles.size();
    do
    {
      _dt = _dt / 2.0;
      if ( tempN > upperLimit )
      {
        _timenext -= _dt;
      }
      else
      {
        _timenext += _dt;
      }

      tempN = 0;
      for ( unsigned int j = 0; j < particles.size(); j++ )
      {
        if ( particles[j].Pos.T() < _timenext )
        {
          ++tempN;
        }
      }
    }
    while (( tempN < etaBins.getNinEtaBin() ) || ( tempN > upperLimit ) );
  }

  _min = min;
  return _timenext;
}

/**
 * Here we build the transverse cell structure for the whole simulation. The transverse cell structure is static. Due to the
 * small expansion in transverse direction this is fine.
 */
void heavyIonCollision::buildTransverseCellConfiguration(double &_dx, double &_dy, int &_IX, int &_IY, double &_transLen)
{
  if ( theConfig->isSet_transCellSize() )
  {
    _dx = _dy = theConfig->getTransCellSize();//set by input file
    _transLen = theConfig->getTransLength();
    //_transLen = 20.0;//is default value    
  }
  else
    {
      //-------------------- build cell configuration --------------------
      _dx = _dy = 0.5;  // default values
      _transLen = 20.0; 
      
      // assuming ~170000 total particles and IZ ~ 40 these values are good choices for dx and dy
      // see estimateCellSizes.nb
      if ( theConfig->getImpactParameter() <= 5.5 )
      {
        _dx = _dy = 0.5;
        _transLen = 20.0;
      }
      else if ( theConfig->getImpactParameter() <= 7.8 )
      {
        _dx = _dy = 0.4;
        _transLen = 20.0;
      }
      else if ( theConfig->getImpactParameter() <= 9.2 )
      {
        _dx = _dy = 0.3;
        _transLen = 20.1;
      }
      else if ( theConfig->getImpactParameter() <= 11.5 )
      {
        _dx = _dy = 0.2;
        _transLen = 20.0;
      }
      else if ( theConfig->getImpactParameter() <= 12.4 )
      {
        _dx = _dy = 0.1;
        _transLen = 20.0;
      }
      else if ( theConfig->getImpactParameter() <= 14.0 )
      {
        _dx = _dy = 0.1;
    //     dx = dy = 0.05;
        _transLen = 20.0;
      }
    }
  
  //--------------------------------------
  _IX = _IY = int( _transLen / dx );
  //--------------------------------------    
  cout << "Geometry:" << endl;
  cout << "dx=" << dx  << " fm" << endl;
  cout << "dy=" << dy  << " fm" << endl;
  cout << "Transverse length = " << transLen << " fm" << endl;
}

/**
 * In case a time step needs to be undone, this routine restores cells to their state at the start of the time step that
 * have actually been modified (as opposed to making a backup of ALL cells, which is very resource demanding).
 * 
 * @param[in,out] _cellsBackup The backup structure for cells whose content has changed within the current time step
 */
void heavyIonCollision::restoreCellBackup( std::map< int, cellContainer >& _cellsBackup )
{
  map<int,cellContainer>::const_iterator myIt;
  for ( myIt = _cellsBackup.begin(); myIt != _cellsBackup.end(); myIt++ )
  {
    cells[ (*myIt).first ] = (*myIt).second;
  }
  
  _cellsBackup.clear();
}

/** Scheis auf Zellen*/
void heavyIonCollision::getCentralValues(double time)
{
  vector<hydroParticleTypeProperties> vecTypeCommon(theHydroParticleType.hydroParticleTypeVectorCommon);
  
  int nTypes = vecTypeCommon.size();
  int nTypesPlusComplete = nTypes + 1;  
  
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
  
  vector<int> NN(nTypesPlusComplete,0);
  
 
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

  int countMidRapParticles=0;
  double pzSum=0.;
  double DeltaEta=0.02;
  
  cout << "ANALYSIS:  particle size " << particles.size() << "\t\t\t TIME = " << time << endl;
  
/*  for(int i=0; i<etaBins.size();i++)
  {
    cout << etaBins[i].left << "\t" <<  etaBins[i].right << endl;
  } */ 
  
//WARNING sinh or cosh?
  double  dz=time*tanh(DeltaEta*2.);
          dz=time*sinh(DeltaEta*2.);
  int IXY = IX*IY;
  double volume  = dz*dx*dy*IXY; //total slice volume
  
  for(int i =0; i<particles.size();i++)
  {
//     double eta=particles[i].Pos.Rapidity();
    
    
    double t=particles[i].Pos.T();
    double zz=particles[i].Pos.Z();
    double eta=0.5*log( fabs((t+zz)/(t-zz)) );    

    if (fabs(eta)<DeltaEta && FPT_COMP_LE(t, time))
    {
//       cout << t << "\t" << zz << "\t" << eta << endl;
      countMidRapParticles++;
      pzSum += particles[i].Mom.Pz()*particles[i].Mom.Pz()/particles[i].Mom.E();
      getTmunuSingleParticle( i, volume, theHydroParticleType, vecTypeCommon, T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, N0, N1, N2, N3, NN, eta );
      
    }
  }
  if(countMidRapParticles!=NN[nTypes]) 
    cout << "SOMETHING WENT WRONG" << endl;
  
  
  theHydroParticleType.giveHydroObservablesInLandauFrame( vecTypeCommon, T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, N0, N1, N2, N3,
                                                          energyDensity, particleDensity, entropyDensity, temperature, fugacity, eq_particleDensity,
                                                          isoPressure, eqPressure, bulkPressure, vx, vy, vz, vv, gamma,
                                                          PI00, PI11, PI22, PI33, PI10, PI20, PI30, PI21, PI31, PI32, PImunuPImunu, q0, q1, q2, q3,
                                                          W0, W1, W2, W3, V0, V1, V2, V3); 
  
  cout << "countMidRapParticles:          " << countMidRapParticles     << endl;
  cout << "pzSum/volume/testpartcl:       " << pzSum/volume/testpartcl  << endl;
  
  
  cout << endl;
  cout << "TESTPARTICLES                = " << testpartcl << endl;
  cout << "Volume central slice [fm^3]  = " << volume << endl;
  cout << "Volume central cell [fm^3]   = " << volume/IXY << endl;
  cout << "Dens. from part/vol [GeV^3]  = " << (NN[nTypes]/volume/testpartcl)*pow(0.197,3.0) << endl;
  cout << "Manualcount/vol [GeV^3]      = " << (countMidRapParticles/volume/testpartcl)*pow(0.197,3.0) << endl;
  cout << "Pressure [GeV^4]             = " << isoPressure[nTypes]*pow(0.197,3.0)     << endl;
  cout << "Density [GeV^3]              = " << particleDensity[nTypes]*pow(0.197,3.0)  << endl;           
  cout << "T=P/n [GeV]                  = " << (isoPressure[nTypes]*pow(0.197,3.0) )/(particleDensity[nTypes]*pow(0.197,3.0) )     << endl;
  cout << "T Janni [GeV]                = " << temperature[nTypes]      << endl;
  cout << "16/PI^2T^3 [GeV^3]           = " << pow(temperature[nTypes],3.0)*16./pow(M_PI,2.0) << endl;
  cout << "eq_particleDensity [GeV^3]   = " << eq_particleDensity[nTypes]*pow(0.197,3.0) << endl;
  cout << "Fuacity = n/neq              = " << particleDensity[nTypes]/eq_particleDensity[nTypes]   << endl;
  cout << "Fugacity Janni               = " << fugacity[nTypes] << endl;
  cout << "ENTROPY JANNI GeV^3          = " << entropyDensity[nTypes]*pow(0.197,3.0)  << endl;
  cout << "( 4- mu/T )*n                = " << (4.-log(fugacity[nTypes]))*particleDensity[nTypes]*pow(0.197,3.0) << endl;
  cout << endl;
  cout << "ENERGY DENSITY Landau        = " << energyDensity[nTypes]*pow(0.197,3.0) << " GeV^4" << endl;   
//   cout << "ENERGY DENSITY simple        = " << EdensTot << " GeV^4" << endl;   
  cout << "ENERGY DENSITY T00           = " << T00[nTypes]*pow(0.197,3.0) << " GeV^4" << endl;   
// 
//   string filenameHydroCentral3 = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + "_HydroLandau_" + ".dat";
//   fstream file3(filenameHydroCentral3.c_str(),ios::app);
//   file3 << time << "\t" << time/0.197 * pow(temperature[nTypes],3.0) << endl;
//   file3.close();
  
  double pi = (PI11[nTypes]*2.0+PI22[nTypes]*2.0+ (-PI33[nTypes])) /3.;  //GeV/fm^3
//          pi = pi*pow(0.197,3.0)/(isoPressure[nTypes]*pow(0.197,3.0)); // scaled by temperature
  double p  = isoPressure[nTypes];
  double anisotropy = (p-pi)/(p+pi/2);
  double anisotropy2 = T33[nTypes]/((T22[nTypes]+T11[nTypes])/2.);
  double pt = ((T22[nTypes]+T11[nTypes])/2.); //GeV/fm³
//          pt = pt*pow(0.197,3.0)/pow(temperature[nTypes],4.0);// scaled by temperature
  double pl = T33[nTypes];                    //GeV/fm³
//          pl = pl*pow(0.197,3.0)/pow(temperature[nTypes],4.0);// scaled by temperature
  double enthalpie = energyDensity[nTypes] + isoPressure[nTypes]; //GeV/fm³
  cout << "SHOULD ALL BE THE SAME: "<< endl;
  cout << PI11[nTypes]*2.0 << "\t" << PI22[nTypes]*2.0 << "\t" << -PI33[nTypes] << endl;
  cout << "PL/T^4 = " << pl << endl;
  cout << "PT/T^4 = " << pt << endl;
  cout << "pi/T^4 = " << pi << endl;  
  
  
  double AndrejPLPT = (p-pi)/(p+pi/2.);
  
//   
// /*  cout << anisotropy << "\t\t" << anisotropy2 << "\t\t\t" << anisotropy/anisotropy2 << endl; 
//   cout << (-PI33[nTypes]) << endl;
//   cout << "PL= " << T33[nTypes]/pow(temperature[nTypes],4.0) << "  PT= " << (T22[nTypes]+T11[nTypes])/(2.*pow(temperature[nTypes],4.0)) << endl;
//   cout << gamma[nTypes]   << "\t" <<   vx[nTypes] << "\t" << vy[nTypes] << "\t" << vz[nTypes] << endl;
//  */ 
//   cout << sqrt(sumAll)<<"\t\t" << sqrt(Tsq/count) << "\t\t" << temperature[nTypes] << endl;
// //   cout << (pt-pl)/enthalpie << "\t\t" << NN[nTypes] << endl;
// //   cout << time/0.197 * pow(temperature[nTypes],3.0) << " GeV²" << endl;
// //   cout << 4*particleDensity[nTypes]*(-pow(temperature[nTypes],2.0)*pi)/(  )
//   
// //   cell_id = nx + IX * ny + IX * IY * nz;
//   
  string filenameHydroCentral = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + "_Hydro2_" + ".dat";
  fstream file(filenameHydroCentral.c_str(),ios::app);
  
//   1: time [fm]
//   2: T [GeV]
//   3: s [Gev^3]
//   4: lambda
//   5: n [Gev^3]
//   6: pi [Gev^4]
//   7: p  [Gev^4]
//   8: e  [Gev^4]
//   9: pl [Gev^4]
//   10:pt [Gev^4]
  file << time << "\t" << temperature[nTypes] << "\t" << entropyDensity[nTypes]*pow(0.197,3.0) << "\t" << fugacity[nTypes] << "\t" << particleDensity[nTypes]*pow(0.197,3.0) << "\t" << pi*pow(0.197,3.0) << "\t" << isoPressure[nTypes]*pow(0.197,3.0)<< "\t" << energyDensity[nTypes]*pow(0.197,3.0) << "\t" << (p-pi)*pow(0.197,3.0) << "\t" << (p+pi/2.)*pow(0.197,3.0) << endl;
 
//   file << time << "\t" << pi*pow(0.197,3.0)/(isoPressure[nTypes]*pow(0.197,3.0)) << "\t" << (pi*pow(0.197,3.0))/(temperature[nTypes]*(entropyDensity[nTypes]*pow(0.197,3.0)) ) << "\t" <<  (time/0.197*temperature[nTypes])*(pi*pow(0.197,3.0))/(temperature[nTypes]*entropyDensity[nTypes]*pow(0.197,3.0)) << "\t" << temperature[nTypes] << "\t" << entropyDensity[nTypes]*pow(0.197,3.0) << "\t" << fugacity[nTypes] << "\t" << particleDensity[nTypes]*pow(0.197,3.0) << endl;
  
//   << (time/0.197*temperature[nTypes])*(pi*pow(0.197,3.0))/(temperature[nTypes]*entropyDensity[nTypes]*pow(0.197,3.0)) << "\t" << entropyDensity[nTypes]*pow(0.197,3.0) << "\t" << p << "\t" << pi/(p) << "\t" << NN[nTypes] << "\t" << (pt-pl)/enthalpie <<  "\t" << time/0.197 * pow(temperature[nTypes],3.0) << "\t" << pl/pt << "\t" <<  pow(energyDensity[nTypes],3./4.)*time/0.197 << "\t" << particleDensity[nTypes]*time/0.197 << "\t" << pl << "\t" << fugacity[nTypes] << endl;
  
  
//   file << time << "\t" << temperature[nTypes] << "\t" << p << "\t" << pi/(p) << "\t" << NN[nTypes] << "\t" << (pt-pl)/enthalpie <<  "\t" << time/0.197 * pow(temperature[nTypes],3.0) << "\t" << pl/pt << "\t" <<  pow(energyDensity[nTypes],3./4.)*time/0.197 << "\t" << particleDensity[nTypes]*time/0.197 << "\t" << pl << "\t" << fugacity[nTypes] << endl;
  file.close(); 
  
  string filenameHydroCentral2 = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + "_Hydro3_" + ".dat";
  fstream file2(filenameHydroCentral2.c_str(),ios::app);
  file2 << time << "\t" << AndrejPLPT << "\t" << temperature[nTypes] << "\t" << pi/energyDensity[nTypes] << "\t" << endl;
  file2.close();
 }

void heavyIonCollision::getCentralCellValues(double time)
{
  vector<hydroParticleTypeProperties> vecTypeCommon(theHydroParticleType.hydroParticleTypeVectorCommon);
  
//   int x_ind_half,y_ind_half,z_ind_half;
//   int x_ind_half2,y_ind_half2,z_ind_half2;
//   bool twocells_x=false;
//   bool twocells_y=false;
//   bool twocells_z=false;
//   int avgSizeL = 5;
//   int avgSizeT = 7;
//   std::vector<int>  cellsForLoop_x(avgSizeT,0);
//   std::vector<int>  cellsForLoop_y(avgSizeT,0);
//   std::vector<int>  cellsForLoop_z(avgSizeL,0);
//   
//   if(IX%2==0)
//   {
//     x_ind_half=IX/2;
//     x_ind_half2=IX/2+1;
//     twocells_x=true;
//   }else
//   {
//     x_ind_half=IX/2;
//     x_ind_half2=-1;
//     std::vector<int> temp={IX/2-3,IX/2-2,IX/2-1,IX/2,IX/2+1,IX/2+2,IX/2+3};
//     cellsForLoop_x = temp;
//   }
//   
//   if(IY%2==0)
//   {
//     y_ind_half=IY/2;
//     y_ind_half2=IY/2+1;
//     twocells_y=true;
//   }else
//   {
//     y_ind_half=IX/2;
//     y_ind_half2=-1;
//     std::vector<int> temp={IY/2-3,IY/2-2,IY/2-1,IY/2,IY/2+1,IY/2+2,IY/2+3};
//     cellsForLoop_y = temp;
//   }
//   
//   if(IZ%2==0)
//   {
//     z_ind_half=IZ/2;
//     z_ind_half2=IZ/2+1;
//     twocells_z=true;
//   }else
//   {
//     z_ind_half=IZ/2;
//     z_ind_half2=-1;
//     std::vector<int> temp={IZ/2-2,IZ/2-1,IZ/2,IZ/2+1,IZ/2+2};
//     cellsForLoop_z = temp;
//   }  
//   
//   

  
  
/*  int cellsForLoop_x[2]={x_ind_half,x_ind_half2};
  int cellsForLoop_y[2]={y_ind_half,y_ind_half2};
  int cellsForLoop_z[2]={z_ind_half,z_ind_half2};
 */ 
  int cell_id;
  
  int nTypes = vecTypeCommon.size();
  int nTypesPlusComplete = nTypes + 1;  
  
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
  
  vector<int> NN(nTypesPlusComplete,0);
  
 
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
  
//   cout << cellsForLoop_x[0] << "\t" << cellsForLoop_x[1] << endl;
//   cout << cellsForLoop_y[0] << "\t" << cellsForLoop_y[1]  << endl;
//   cout << cellsForLoop_z[0] << "\t" << cellsForLoop_z[1] << endl;
  int IXY = IX*IY;
  
  int etaIndex = etaBins.getCentralIndex();
  double dz  = time * ( tanh( etaBins[etaIndex].right ) - tanh( etaBins[etaIndex].left ) );
  double dz2 = 0.2  * ( sinh( etaBins[etaIndex].right ) - sinh( etaBins[etaIndex].left ) );
  
//   for(int i=0; i<etaBins.size();i++)
//   {
//     cout << etaBins[i].left << "\t" <<  etaBins[i].right << endl;
//   }
//   cout << "COMPARE dz:  tanh: " << dz << "    sinh: " << dz2 << "    and transverse area: " << dx*dy*IXY << " fm^2" << endl;
//   
  //dz = time  *  (  sinh( eta+dRap/2. )-  sinh( eta-dRap/2. )  )
  
//   double dz = time * (etaBins[etaIndex].right - etaBins[etaIndex].left);
  double volume  = dz*dx*dy*IXY; //total slice volume
  
  //HACK
  //volume = IXY*0.07147;
  
  double Tsq = 0.;
  double EdensTot=0.;
  int count=0;
  double sumAll=0.;
  int countManually=0;
  for ( int cell_id = IXY * etaIndex; cell_id < IXY * ( etaIndex + 1 ); cell_id++ )
  {
    //TEST
//    list<int>::const_iterator iIt; 
//    for ( iIt = cells[cell_id].particleList.begin(); iIt != cells[cell_id].particleList.end(); iIt++ )
//    {
//      int mm = *iIt;
//      if(particles[mm].Pos.T()<=time)
//      {
//        countManually++;
//      }
//    }
    
//     double Tsqtemp=0.;
//     getTemperatureAlternative( cells[cell_id], dz*dx*dy, theHydroParticleType, vecTypeCommon, Tsqtemp, sumAll);
//     Tsq += Tsqtemp;
//     count++;
    
//     if(FPT_COMP_GE(cells[cell_id].particleList.size(),15)) //WARNING hard particle cut
//     {
//       getTmunuInCell( cells[cell_id], dz*dx*dy, theHydroParticleType, vecTypeCommon, T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, N0, N1, N2, N3, NN );
//     }
    
//     if(cells[cell_id].corner.x_min > -1.1 && cells[cell_id].corner.x_max < 1.1 && cells[cell_id].corner.y_min > -1.1 && cells[cell_id].corner.y_max < 1.1)
//     {
//       if(FPT_COMP_LE(cells[cell_id].particleList.size(),20))
//       {
//         continue;
//       } 
//       cout << cells[cell_id].particleList.size() << endl;  
      count++;
//       cout << cells[cell_id].particleList.size() << endl;
      getTmunuInCell( cells[cell_id], volume, theHydroParticleType, vecTypeCommon, T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, N0, N1, N2, N3, NN );
      
      //Tmunu GeV/fm³
      //Numu 1/fm³
      
      
      double Tsqtemp=0.;
      double EdensHere=0.;
      getTemperatureAlternative( cells[cell_id], volume, theHydroParticleType, vecTypeCommon, Tsqtemp, sumAll, EdensHere);
      Tsq += Tsqtemp;
      EdensTot+=EdensHere; 
//     }
    
  }
  cout << "CELL No: = " << count << "   NN = " << NN[nTypes] << "     formed=" << countManually << endl;
  sumAll /= volume*testpartcl*gG/(2.*pow(M_PI,2.0));
  
  theHydroParticleType.giveHydroObservablesInLandauFrame( vecTypeCommon, T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, N0, N1, N2, N3,
                                                          energyDensity, particleDensity, entropyDensity, temperature, fugacity, eq_particleDensity,
                                                          isoPressure, eqPressure, bulkPressure, vx, vy, vz, vv, gamma,
                                                          PI00, PI11, PI22, PI33, PI10, PI20, PI30, PI21, PI31, PI32, PImunuPImunu, q0, q1, q2, q3,
                                                          W0, W1, W2, W3, V0, V1, V2, V3); 
  cout << endl;
  cout << "TESTPARTICLES                = " << testpartcl << endl;
  cout << "Volume central slice [fm^3]  = " << volume << endl;
  cout << "Volume central cell [fm^3]   = " << volume/IXY << endl;
  cout << "Dens. from part/vol [GeV^3]  = " << (NN[nTypes]/volume/testpartcl)*pow(0.197,3.0) << endl;
  cout << "Manualcount/vol [GeV^3]      = " << (countManually/volume/testpartcl)*pow(0.197,3.0) << endl;
  cout << "Pressure [GeV^4]             = " << isoPressure[nTypes]*pow(0.197,3.0)     << endl;
  cout << "Density [GeV^3]              = " << particleDensity[nTypes]*pow(0.197,3.0)  << endl;           
  cout << "T=P/n [GeV]                  = " << (isoPressure[nTypes]*pow(0.197,3.0) )/(particleDensity[nTypes]*pow(0.197,3.0) )     << endl;
  cout << "T Janni [GeV]                = " << temperature[nTypes]      << endl;
  cout << "16/PI^2T^3 [GeV^3]           = " << pow(temperature[nTypes],3.0)*16./pow(M_PI,2.0) << endl;
  cout << "eq_particleDensity [GeV^3]   = " << eq_particleDensity[nTypes]*pow(0.197,3.0) << endl;
  cout << "Fuacity = n/neq              = " << particleDensity[nTypes]/eq_particleDensity[nTypes]   << endl;
  cout << "Fugacity Janni               = " << fugacity[nTypes] << endl;
  
  cout << endl;
  cout << "ENERGY DENSITY Landau        = " << energyDensity[nTypes]*pow(0.197,3.0) << " GeV^4" << endl;   
  cout << "ENERGY DENSITY simple        = " << EdensTot << " GeV^4" << endl;   
  cout << "ENERGY DENSITY T00           = " << T00[nTypes]*pow(0.197,3.0) << " GeV^4" << endl;   
// 
//   string filenameHydroCentral3 = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + "_HydroLandau_" + ".dat";
//   fstream file3(filenameHydroCentral3.c_str(),ios::app);
//   file3 << time << "\t" << time/0.197 * pow(temperature[nTypes],3.0) << endl;
//   file3.close();
  
  double pi = (PI11[nTypes]*2.0+PI22[nTypes]*2.0+ (-PI33[nTypes])) /3.;
  
  cout << "SHOULD ALL BE THE SAME: "<< endl;
  cout << PI11[nTypes]*2.0 << "\t" << PI22[nTypes]*2.0 << "\t" << -PI33[nTypes] << endl;
  
  
  double p  = isoPressure[nTypes];
  double anisotropy = (p-pi)/(p+pi/2);
  double anisotropy2 = T33[nTypes]/((T22[nTypes]+T11[nTypes])/2.);
  
  double pt = ((T22[nTypes]+T11[nTypes])/2.); //GeV/fm³
  double pl = T33[nTypes];                    //GeV/fm³
  double enthalpie = energyDensity[nTypes] + isoPressure[nTypes]; //GeV/fm³
  
  cout << "PL = " << pl << endl;
  cout << "PT = " << pt << endl;
  
  
  
  
  
  
//   
// /*  cout << anisotropy << "\t\t" << anisotropy2 << "\t\t\t" << anisotropy/anisotropy2 << endl; 
//   cout << (-PI33[nTypes]) << endl;
//   cout << "PL= " << T33[nTypes]/pow(temperature[nTypes],4.0) << "  PT= " << (T22[nTypes]+T11[nTypes])/(2.*pow(temperature[nTypes],4.0)) << endl;
//   cout << gamma[nTypes]   << "\t" <<   vx[nTypes] << "\t" << vy[nTypes] << "\t" << vz[nTypes] << endl;
//  */ 
//   cout << sqrt(sumAll)<<"\t\t" << sqrt(Tsq/count) << "\t\t" << temperature[nTypes] << endl;
// //   cout << (pt-pl)/enthalpie << "\t\t" << NN[nTypes] << endl;
// //   cout << time/0.197 * pow(temperature[nTypes],3.0) << " GeV²" << endl;
// //   cout << 4*particleDensity[nTypes]*(-pow(temperature[nTypes],2.0)*pi)/(  )
//   
// //   cell_id = nx + IX * ny + IX * IY * nz;
//   
  string filenameHydroCentral = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + "_Hydro_" + ".dat";
  fstream file(filenameHydroCentral.c_str(),ios::app);
  
//   1: time
//   2: Number
//   3: pt-pl/enthalpie
//   4: T^3tau
//   5: pl/pt
//   6: e^(3/4)*tau
//   7: n*tau
//   8: pl
  
  
  file << time << "\t" << temperature[nTypes] << "\t" << p << "\t" << pi/(p) << "\t" << NN[nTypes] << "\t" << (pt-pl)/enthalpie <<  "\t" << time/0.197 * pow(temperature[nTypes],3.0) << "\t" << pl/pt << "\t" <<  pow(energyDensity[nTypes],3./4.)*time/0.197 << "\t" << particleDensity[nTypes]*time/0.197 << "\t" << pl << endl;
  file.close(); 
  
//   //CentralCell
//   int halfX_index = IX/2;
//   int HalfCellNumberCentralBox = 5; 
//   cout << "Central Cells:" << endl;
//   for ( int cell_id = IXY * etaIndex; cell_id < IXY * ( etaIndex + 1 ); cell_id++ )
//   {
//     if(cells[cell_id].corner.x_min > -1. && cells[cell_id].corner.x_max < 1 && cells[cell_id].corner.y_min > -1. && cells[cell_id].corner.y_max < 1)
//     {
//       cout << cells[cell_id].corner.x_min << "\t" <<  cells[cell_id].corner.x_max << "\t" << cells[cell_id].corner.y_min << "\t" <<  cells[cell_id].corner.y_max << "\t" << cells[cell_id].particleList.size() << endl;
//       
//     }
//   }  
//    
//   etaIndex = etaBins.getCentralIndex()-1;
//   dz = time * ( tanh( etaBins[etaIndex].right ) - tanh( etaBins[etaIndex].left ) );
//   volume  = dz*dx*dy*IXY; //total slice volume
//   
//   for ( int cell_id = IXY * etaIndex; cell_id < IXY * ( etaIndex + 1 ); cell_id++ )
//   {
//     getTmunuInCell( cells[cell_id], dv, theHydroParticleType, vecTypeCommon, T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, N0, N1, N2, N3, NN );
//   }
//   
//   etaIndex = etaBins.getCentralIndex()+1;
//   dz = time * ( tanh( etaBins[etaIndex].right ) - tanh( etaBins[etaIndex].left ) );
//   volume  = dz*dx*dy*IXY; //total slice volume
//   
//   for ( int cell_id = IXY * etaIndex; cell_id < IXY * ( etaIndex + 1 ); cell_id++ )
//   {
//     getTmunuInCell( cells[cell_id], dv, theHydroParticleType, vecTypeCommon, T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, N0, N1, N2, N3, NN );
//   }  
  
  
/*  for(int x=0;x<avgSizeT;x++)
  {
    if(cellsForLoop_x[x]<0)
      continue;
      for(int y=0;y<avgSizeT;y++)
      { 
        if(cellsForLoop_y[y]<0)
          continue;
        for(int z=0;z<avgSizeL;z++)
        {
          if(cellsForLoop_z[z]<0)
            continue;
          cell_id = cellsForLoop_x[x] + IX * cellsForLoop_y[y] + IX * IY * cellsForLoop_z[z];
//           cells[cell_id];
          getTmunuInCell( cells[cell_id], dv, theHydroParticleType, vecTypeCommon, T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, N0, N1, N2, N3, NN );

        }    
      }
  } */ 
  
  
  
  
//   for(int x=0;x<avgSizeT;x++)
//   {
//     if(cellsForLoop_x[x]<0)
//       continue;
//       for(int y=0;y<avgSizeT;y++)
//       { 
//         if(cellsForLoop_y[y]<0)
//           continue;
//         for(int z=0;z<avgSizeL;z++)
//         {
//           if(cellsForLoop_z[z]<0)
//             continue;
//           cell_id = cellsForLoop_x[x] + IX * cellsForLoop_y[y] + IX * IY * cellsForLoop_z[z];
// //           cells[cell_id];
//           getTmunuInCell( cells[cell_id], dv, theHydroParticleType, vecTypeCommon, T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, N0, N1, N2, N3, NN );
// 
//         }    
//       }
//   }

//   theHydroParticleType.giveHydroObservablesInEckartFrame(vecTypeCommon, T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, N0, N1, N2, N3,
//                                 energyDensity, particleDensity, entropyDensity, temperature, fugacity, eq_particleDensity,
//                                 isoPressure, eqPressure, bulkPressure, vx, vy, vz, vv, gamma,
//                                 PI00, PI11, PI22, PI33, PI10, PI20, PI30, PI21, PI31, PI32, PImunuPImunu, q0, q1, q2, q3,
//                                 W0, W1, W2, W3, V0, V1, V2, V3);

//   string filenameHydroCentral2 = theConfig->getStandardOutputDirectoryName() + "/" + theConfig->getJobName() + "_HydroEckart_" + ".dat";
//   fstream file2(filenameHydroCentral2.c_str(),ios::app);
//   file2 << time << "\t" << time/0.197 * pow(temperature[nTypes],3.0) << "\t" << temperature[nTypes] <<  endl;
//   file2.close();
 

 }
/**
 * This routine prepares the the calculation of the cross sections for non pQCD calculations. In case for pQCD nothing happens.
 * In case another method is chosen, several routines are started via a switch method.
 */
double heavyIonCollision::getSpecialCS(cellContainer& _cell, const double cellVolume, hydroParticleType & theHydroParticleType, const vector<hydroParticleTypeProperties> & vecType, double & scalingFactor, double & extractedTeff)
{ 
  int nTypes = vecType.size();
  int nTypesPlusComplete = nTypes + 1;  
  
  
  double specialCS;
    
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
  
  vector<int> NN(nTypesPlusComplete,0);
  
  //-----------------------------//
  scattering22_hydro sc22;
  double inputValueCS = theConfig->getInputCrossSectionValue();
  //-----------------------------//  
  
  //-----------------------------//
  //get the real time
  double realTime = timenext - timeshift;
  //-----------------------------//
  
  //---------------------------------//
  scalingFactor = 1.0;//standard values, if not changed
  specialCS = 0.0;//standard values, if not changed
  
  double extractedN, extractedFug;
  
  switch ( theConfig->getCrossSectionMethod() )
  {
    case csMethod_pQCD:
      //do nothing
      break;
    case csMethod_constCS:
      specialCS = inputValueCS / pow(0.197,2) / 10.0;//1/GeV^2
      break;
    case csMethod_constEtaOverS:
      if(FPT_COMP_GE(_cell.particleList.size(),10)) //WARNING hard particle cut
      {
        getTmunuInCell( _cell, dv, theHydroParticleType, vecType, T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, N0, N1, N2, N3, NN );
        //getTmunuAtMidrapidity( realTime, dv, theHydroParticleType, vecType, T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, N0, N1, N2, N3, NN );
        
        specialCS = sc22.getConstEtaOverS_XSection( theHydroParticleType, vecType, inputValueCS, T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, N0, N1, N2, N3, NN, extractedTeff, extractedN, extractedFug );      
        //cout << _cell.corner.etaIndex << endl;
        int etaIndex = _cell.corner.etaIndex;
        double  dz  = timenext * ( tanh( etaBins[etaIndex].right ) - tanh( etaBins[etaIndex].left ) );
        double  mfp = 1./(extractedN*(specialCS*pow(0.197,2.0)));
        
        if(etaIndex == etaBins.getCentralIndex())
        {
          countMFP++;
          MFPRatio+=mfp/dz;
          MFPfm+=mfp;
          cout << "T/GeV=" << extractedTeff << "   MFP = " << mfp << " fm,   dz = " << dz << "   " << "Cell = " << etaBins[etaIndex].left << " --- " << etaBins[etaIndex].right << endl;
        }
      }
      else
      {
        specialCS = 0.;
      }
      break;    
    case csMethod_constMFP:
      getTmunuInCell( _cell, dv, theHydroParticleType, vecType, T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, N0, N1, N2, N3, NN );
      specialCS = sc22.getConstMFP_XSection( theHydroParticleType, inputValueCS, T00[nTypes], T11[nTypes], T22[nTypes], T33[nTypes], T10[nTypes], T20[nTypes], T30[nTypes], T21[nTypes], T31[nTypes], T32[nTypes], N0[nTypes], N1[nTypes], N2[nTypes], N3[nTypes], NN[nTypes]);
      break;    
    case csMethod_constMixtureCS:
      //do nothing
      break;
    case csMethod_constMixtureCS_scaledWithLambdaT2:
      getTmunuInCell( _cell, dv, theHydroParticleType, vecType, T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, N0, N1, N2, N3, NN );
      scalingFactor = sc22.getScalingFactorLambdaT2( theHydroParticleType, vecType, T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, N0, N1, N2, N3, NN );
      break;  
    case csMethod_variableCS_scaledWithLambdaT2:
      getTmunuInCell( _cell, dv, theHydroParticleType, vecType, T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, N0, N1, N2, N3, NN );
      scalingFactor = sc22.getScalingFactorLambdaT2( theHydroParticleType, vecType, T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, N0, N1, N2, N3, NN );
      specialCS = sc22.getVarTimeDepCS( realTime );
      break;      
    default:
      string errMsg = "Error in method of cross section";
      throw eHIC_error( errMsg );
  }

return specialCS;
}


void heavyIonCollision::getTemperatureAlternative(cellContainer& _cell, const double cellVolume, hydroParticleType & theHydroParticleType, const vector<hydroParticleTypeProperties> & vecType, double & Tsq , double & sumUp, double & EdensGeV4)
{
  int nTypes = vecType.size();
  
  list<int>::const_iterator iIt;  

  //Total number of particles:
  int countParticlesInCell = *(_cell.particleList.end());
  
  double sum=0;
  double sumE = 0.;
  double md2q=0;
  double md2g = 0.;
//       int type = theHydroParticleType.getParticleType(particles[mm].FLAVOR,vecType);
    
  double vt = cellVolume * testpartcl;
  for ( iIt = _cell.particleList.begin(); iIt != _cell.particleList.end(); iIt++ )
  {    
    int mm = *iIt;
    sum += 1.0 / particles[mm].Mom.E();
    sumE+=particles[mm].Mom.E();
  }
  
  sumUp += sum;
  sum = sum / gG / vt;
  
  Tsq = 2.*pow(M_PI,2.0)*sum*pow(0.197,3.0); //GeV^2
  EdensGeV4 = sumE/vt*pow(0.197,3.0);
}

/**
 * This routine calculates the Tmunu and Nmu for a specific cell.
 */
void heavyIonCollision::getTmunuInCell(cellContainer& _cell, const double cellVolume, hydroParticleType & theHydroParticleType, const vector<hydroParticleTypeProperties> & vecType,
         vector<double>& T00,  vector<double>& T11,  vector<double>& T22,
         vector<double>& T33,  vector<double>& T10,  vector<double>& T20,
         vector<double>& T30,  vector<double>& T21,  vector<double>& T31,
         vector<double>& T32,  vector<double>& N0,  vector<double>& N1,
         vector<double>& N2,  vector<double>& N3,  vector<int>& NN )
{  
  int nTypes = vecType.size();
  
  list<int>::const_iterator iIt;  

  //Total number of particles:
  int countParticlesInCell = *(_cell.particleList.end());
  
  
  
  for ( iIt = _cell.particleList.begin(); iIt != _cell.particleList.end(); iIt++ )
  {
    
    int mm = *iIt;

    int type = theHydroParticleType.getParticleType(particles[mm].FLAVOR,vecType);
    
    double vt = cellVolume * testpartcl;

    
    if(type < 0)
    {
      cout << "Particle not in particle list for Tmunu" << endl;
    }
    else
    {
      T00[type] += particles[mm].Mom.E()  / vt; //GeV/fm^3
      T11[type] += pow( particles[mm].Mom.Px(),2 ) / particles[mm].Mom.E()  / vt; //GeV/fm^3
      T22[type] += pow( particles[mm].Mom.Py(),2 ) / particles[mm].Mom.E()  / vt; //GeV/fm^3
      T33[type] += pow( particles[mm].Mom.Pz(),2 ) / particles[mm].Mom.E()  / vt; //GeV/fm^3

      T10[type] += particles[mm].Mom.Px()  / vt; //GeV/fm^3
      T20[type] += particles[mm].Mom.Py()  / vt; //GeV/fm^3
      T30[type] += particles[mm].Mom.Pz()  / vt; //GeV/fm^3

      T21[type] += particles[mm].Mom.Py() * particles[mm].Mom.Px() / particles[mm].Mom.E()  / vt; //GeV/fm^3
      T31[type] += particles[mm].Mom.Pz() * particles[mm].Mom.Px() / particles[mm].Mom.E()  / vt; //GeV/fm^3
      T32[type] += particles[mm].Mom.Pz() * particles[mm].Mom.Py() / particles[mm].Mom.E()  / vt; //GeV/fm^3

      N0[type] += 1.0 / vt; // 1/fm^3
      N1[type] += particles[mm].Mom.Px() / particles[mm].Mom.E() / vt; // 1/fm^3
      N2[type] += particles[mm].Mom.Py() / particles[mm].Mom.E() / vt; // 1/fm^3
      N3[type] += particles[mm].Mom.Pz() / particles[mm].Mom.E() / vt; // 1/fm^3
      NN[type]++;
      
      //for all particles
      T00[nTypes] += particles[mm].Mom.E()  / vt; //GeV/fm^3
      T11[nTypes] += pow( particles[mm].Mom.Px(),2 ) / particles[mm].Mom.E()  / vt; //GeV/fm^3
      T22[nTypes] += pow( particles[mm].Mom.Py(),2 ) / particles[mm].Mom.E()  / vt; //GeV/fm^3
      T33[nTypes] += pow( particles[mm].Mom.Pz(),2 ) / particles[mm].Mom.E()  / vt; //GeV/fm^3

      T10[nTypes] += particles[mm].Mom.Px()  / vt; //GeV/fm^3
      T20[nTypes] += particles[mm].Mom.Py()  / vt; //GeV/fm^3
      T30[nTypes] += particles[mm].Mom.Pz()  / vt; //GeV/fm^3

      T21[nTypes] += particles[mm].Mom.Py() * particles[mm].Mom.Px() / particles[mm].Mom.E()  / vt; //GeV/fm^3
      T31[nTypes] += particles[mm].Mom.Pz() * particles[mm].Mom.Px() / particles[mm].Mom.E()  / vt; //GeV/fm^3
      T32[nTypes] += particles[mm].Mom.Pz() * particles[mm].Mom.Py() / particles[mm].Mom.E()  / vt; //GeV/fm^3

      N0[nTypes] += 1.0 / vt; // 1/fm^3
      N1[nTypes] += particles[mm].Mom.Px() / particles[mm].Mom.E() / vt; // 1/fm^3
      N2[nTypes] += particles[mm].Mom.Py() / particles[mm].Mom.E() / vt; // 1/fm^3
      N3[nTypes] += particles[mm].Mom.Pz() / particles[mm].Mom.E() / vt; // 1/fm^3
      NN[nTypes]++;
      //-------------------------         
    }
  }
  
}

/**
 * This routine calculates the Tmunu and Nmu for a specific cell.
 */
void heavyIonCollision::getTmunuAtMidrapidity(double time_, const double cellVolume, hydroParticleType & theHydroParticleType, const vector<hydroParticleTypeProperties> & vecType,
         vector<double>& T00,  vector<double>& T11,  vector<double>& T22,
         vector<double>& T33,  vector<double>& T10,  vector<double>& T20,
         vector<double>& T30,  vector<double>& T21,  vector<double>& T31,
         vector<double>& T32,  vector<double>& N0,  vector<double>& N1,
         vector<double>& N2,  vector<double>& N3,  vector<int>& NN )
{  
  //TODO  
  vector<hydroParticleTypeProperties> vecTypeCommon(theHydroParticleType.hydroParticleTypeVectorCommon);
  
  int nTypes = vecTypeCommon.size();
  int nTypesPlusComplete = nTypes + 1;  
  
  int countMidRapParticles=0;
  double pzSum=0.;
  double DeltaEta=0.01;

  
  //WARNING sinh or cosh?
  double dz=time_*sinh(DeltaEta*2.);
  int IXY = IX*IY;
  double volume  = dz*dx*dy*IXY; //total slice volume
  
  for(int i =0; i<particles.size();i++)
  {
//     double eta=particles[i].Pos.Rapidity();
    double t=particles[i].Pos.T();
    double zz=particles[i].Pos.Z();
    double eta=0.5*log( fabs((t+zz)/(t-zz)) );    

    if (fabs(eta)<DeltaEta && FPT_COMP_LE(t, time_))
    {
      getTmunuSingleParticle( i, volume, theHydroParticleType, vecTypeCommon, T00, T11, T22, T33, T10, T20, T30, T21, T31, T32, N0, N1, N2, N3, NN, eta );
      
    }
  }
  
}


/**
 * This routine calculates the Tmunu for a particle with index i
 */
void heavyIonCollision::getTmunuSingleParticle(int particleIndex, const double cellVolume, hydroParticleType & theHydroParticleType, const vector<hydroParticleTypeProperties> & vecType,
         vector<double>& T00,  vector<double>& T11,  vector<double>& T22,
         vector<double>& T33,  vector<double>& T10,  vector<double>& T20,
         vector<double>& T30,  vector<double>& T21,  vector<double>& T31,
         vector<double>& T32,  vector<double>& N0,  vector<double>& N1,
         vector<double>& N2,  vector<double>& N3,  vector<int>& NN, double eta )
{  
  int nTypes = vecType.size();
  
  int mm = particleIndex;

  int type = theHydroParticleType.getParticleType(particles[mm].FLAVOR,vecType);
    
  double vt = cellVolume * testpartcl;

//   //BOOST mit eta_i
//   // gamma = cosh(eta)
//   VectorTXYZ boostvector;
//   double beta;
//   double gamma = cosh(eta);
//   if(eta>0.00000000000000000000)
//   {
//     beta=sqrt(pow(gamma,2.0)-1.)/gamma;
//   }else
//   {
//     beta=-sqrt(pow(gamma,2.0)-1.)/gamma;
//   }
//   
//   boostvector.SetTXYZ(1.,0,0,beta);
//   
//   //lorentz L( VectorEPxPyPz( 1.0, 0.0, 0.0, - tanh(eta) ) );
//   
//   
//   lorentz L( boostvector);
//   VectorEPxPyPz BoostedMomentum = L.boost(particles[mm].Mom);
    
  VectorEPxPyPz BoostedMomentum = particles[mm].Mom;
  
  if(type < 0)
  {
    cout << "Particle not in particle list for Tmunu" << endl;
  }
  else
  {
    T00[type] += BoostedMomentum.E()  / vt; //GeV/fm^3
    T11[type] += pow( BoostedMomentum.Px(),2 ) / BoostedMomentum.E()  / vt; //GeV/fm^3
    T22[type] += pow( BoostedMomentum.Py(),2 ) / BoostedMomentum.E()  / vt; //GeV/fm^3
    T33[type] += pow( BoostedMomentum.Pz(),2 ) / BoostedMomentum.E()  / vt; //GeV/fm^3

    T10[type] += BoostedMomentum.Px()  / vt; //GeV/fm^3
    T20[type] += BoostedMomentum.Py()  / vt; //GeV/fm^3
    T30[type] += BoostedMomentum.Pz()  / vt; //GeV/fm^3

    T21[type] += BoostedMomentum.Py() * BoostedMomentum.Px() / BoostedMomentum.E()  / vt; //GeV/fm^3
    T31[type] += BoostedMomentum.Pz() * BoostedMomentum.Px() / BoostedMomentum.E()  / vt; //GeV/fm^3
    T32[type] += BoostedMomentum.Pz() * BoostedMomentum.Py() / BoostedMomentum.E()  / vt; //GeV/fm^3

    N0[type] += 1.0 / vt; // 1/fm^3
    N1[type] += BoostedMomentum.Px() / BoostedMomentum.E() / vt; // 1/fm^3
    N2[type] += BoostedMomentum.Py() / BoostedMomentum.E() / vt; // 1/fm^3
    N3[type] += BoostedMomentum.Pz() / BoostedMomentum.E() / vt; // 1/fm^3
    NN[type]++;
    
    //for all particles
    T00[nTypes] += BoostedMomentum.E()  / vt; //GeV/fm^3
    T11[nTypes] += pow( BoostedMomentum.Px(),2 ) / BoostedMomentum.E()  / vt; //GeV/fm^3
    T22[nTypes] += pow( BoostedMomentum.Py(),2 ) / BoostedMomentum.E()  / vt; //GeV/fm^3
    T33[nTypes] += pow( BoostedMomentum.Pz(),2 ) / BoostedMomentum.E()  / vt; //GeV/fm^3

    T10[nTypes] += BoostedMomentum.Px()  / vt; //GeV/fm^3
    T20[nTypes] += BoostedMomentum.Py()  / vt; //GeV/fm^3
    T30[nTypes] += BoostedMomentum.Pz()  / vt; //GeV/fm^3

    T21[nTypes] += BoostedMomentum.Py() * BoostedMomentum.Px() / BoostedMomentum.E()  / vt; //GeV/fm^3
    T31[nTypes] += BoostedMomentum.Pz() * BoostedMomentum.Px() / BoostedMomentum.E()  / vt; //GeV/fm^3
    T32[nTypes] += BoostedMomentum.Pz() * BoostedMomentum.Py() / BoostedMomentum.E()  / vt; //GeV/fm^3

    N0[nTypes] += 1.0 / vt; // 1/fm^3
    N1[nTypes] += BoostedMomentum.Px() / BoostedMomentum.E() / vt; // 1/fm^3
    N2[nTypes] += BoostedMomentum.Py() / BoostedMomentum.E() / vt; // 1/fm^3
    N3[nTypes] += BoostedMomentum.Pz() / BoostedMomentum.E() / vt; // 1/fm^3
    NN[nTypes]++;
    //-------------------------         
  }
  
  
}


/**
 * This routine gives an info which method of cross section is used and show some additional informations of the used cross sections.
 */
void heavyIonCollision::infoCrossSectionMethod()
{
  hydroParticleType theHydroParticleType;  
  vector<hydroParticleTypeProperties> vecTypeCommon(theHydroParticleType.hydroParticleTypeVectorCommon);    
  
  scattering22_hydro collisions;   
  
  cout << "---------------------------------------------" << endl;
  cout << "  Type of cross section: ";  
  switch ( theConfig->getCrossSectionMethod() )
  {
  case csMethod_pQCD:
    cout << "  pQCD " << endl;    
    break;
  case csMethod_constCS:
    cout << "  constant cross section " << endl;
    cout << "  cs[mb] = " << theConfig->getInputCrossSectionValue() << endl;
    cout << "  cs[1/GeV^2] = " << theConfig->getInputCrossSectionValue() / pow(0.197,2) / 10.0  << endl;    
    break;
  case csMethod_constEtaOverS:
    cout << "  constant eta over s " << endl;
    cout << "  eta/s = " << theConfig->getInputCrossSectionValue() << endl;   
    break;    
  case csMethod_constMFP:
    cout << "  constant mean free path " << endl;
    cout << "  mfp[fm] = " << theConfig->getInputCrossSectionValue() << endl;    
    break;    
  case csMethod_constMixtureCS:
    cout << "  constant mixture cross section " << endl;   
    for(int i = 0; i <= 16; i++)
      {
        if(theHydroParticleType.getParticleType(i,vecTypeCommon) >= 0)
        {
          for(int j = i; j <= 16; j++)
            {
              if(theHydroParticleType.getParticleType(j,vecTypeCommon) >= 0)
              {
                string str;
                double cs22 = collisions.getMixtureXSection22(theHydroParticleType,vecTypeCommon,i,j,1.0,str);
                cout << "  " << str << " " << cs22 * 10.0 * pow(0.197,2) << " mb  --  " << cs22 << " GeV^(-2)" << endl;
              }
            if(j == 1){j = 2;}//to skip the many quarks and show only once
            if(j == 3){j = 12;}//to skip the many quarks and show only once
            }
        }
        if(i == 1){i = 2;}//to skip the many quarks and show only once
        if(i == 3){i = 12;}//to skip the many quarks and show only once
      }    
    break;
  case csMethod_constMixtureCS_scaledWithLambdaT2:
    cout << "  constant mixture cross section, but with a scaling factor fug*T^2 " << endl;
    for(int i = 0; i <= 16; i++)
      {
        if(theHydroParticleType.getParticleType(i,vecTypeCommon) >= 0)
        {
          for(int j = i; j <= 16; j++)
            {
              if(theHydroParticleType.getParticleType(j,vecTypeCommon) >= 0)
              {
                string str;
                double cs22 = collisions.getMixtureXSection22(theHydroParticleType,vecTypeCommon,i,j,1.0,str);
                cout << "  " << str << " " << cs22 * 10.0 * pow(0.197,2) << " mb  --  " << cs22 << " GeV^(-2)" << endl;
              }
            if(j == 1){j = 2;}//to skip the many quarks and show only once
            if(j == 3){j = 12;}//to skip the many quarks and show only once
            }
        }
        if(i == 1){i = 2;}//to skip the many quarks and show only once
        if(i == 3){i = 12;}//to skip the many quarks and show only once
      }        
    break; 
  case csMethod_variableCS_scaledWithLambdaT2:    
    cout << "  Special variable cross section with a scaling factor fug*T^2 " << endl;
    break;
  default:
    string errMsg = "Error in method of cross section";
    throw eHIC_error( errMsg );
  }
  cout << "---------------------------------------------" << endl;
}
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;  replace-tabs on;
