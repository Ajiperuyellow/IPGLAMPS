//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/mfp_data.cpp $
//$LastChangedDate: 2016-09-20 16:39:07 +0200 (Di, 20. Sep 2016) $
//$LastChangedRevision: 2429 $
//$LastChangedBy: senzel $
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>

#include "mfp_data.h"
#include "particleprototype.h"

using namespace std;

mfp_data::mfp_data() : 
  data_lambda_therm_gluon_fm( -1 ), data_lambda_therm_gluon_gev( -1 ),
  data_lambda_therm_quark_fm( -1 ), data_lambda_therm_quark_gev( -1 ),
  data_lambda_therm_charm_fm( -1 ), data_lambda_therm_charm_gev( -1 ),
  data_lambda_therm_bottom_fm( -1 ), data_lambda_therm_bottom_gev( -1 ),
  includeQuarkMFP( false ),
  includeCharmMFP( false ),
  includeBottomMFP( false ),
  read_file( false )
{
  if( read_file )
  {
    if ( readTables( "mfp_table.dat", gluon ) == MFP_SUCCESS )
    {
      cout << "MFP readout successful" << endl;
    }
    else
    {
      string errMsg = "Readout of MFP table (mfp_table.dat) failed. Unrecoverable error.";
      throw eMFP_read_error( errMsg );
    }
  }
}



mfp_data::mfp_data( const double _T, const bool read_file_arg ):
  data_lambda_therm_gluon_fm( -1 ), data_lambda_therm_gluon_gev( -1 ),
  data_lambda_therm_quark_fm( -1 ), data_lambda_therm_quark_gev( -1 ),
  data_lambda_therm_charm_fm( -1 ), data_lambda_therm_charm_gev( -1 ),
  data_lambda_therm_bottom_fm( -1 ), data_lambda_therm_bottom_gev( -1 ),
  includeQuarkMFP( false ),
  includeCharmMFP( false ),
  includeBottomMFP( false ),
  read_file( read_file_arg )
{
  if( read_file )
  {
    string prefix = "data/MFP/mfpTable_";
    string basename;
    if( coupling::isRunning() )
    {
      basename = "asRun_asMax1_0_womdct_ff0_3_GBimproved_nf3_svn2422_T" + getString_T( _T ) + "MeV";
    }
    else
    {
      basename = "as0_3_womdct_ff0_3_GBimproved_nf3_svn2422_T" + getString_T( _T ) + "MeV";
    }

    string filename = prefix + "m0_" + basename + "_gluon.dat";
    if ( readTables( filename, gluon ) == MFP_SUCCESS )
    {
      cout << "MFP readout (gluons) successful: " << filename << endl;
    }
    else
    {
      string errMsg = "Readout of MFP table (" + filename + ") failed. Unrecoverable error.";
      throw eMFP_read_error( errMsg );
    }
    
    filename = prefix + "m0_" + basename + "_quark.dat";
    if ( readTables( filename, light_quark ) == MFP_SUCCESS )
    {
      includeQuarkMFP = true;
      cout << "MFP readout (quarks) successful: "  << filename << endl;
    }
    
    if( ParticlePrototype::N_heavy_flavor > 0 )
    {
      filename = prefix + "m13_" + basename + "_biSectionIter_charm.dat";
      if ( readTables( filename, charm ) == MFP_SUCCESS )
      {
        includeCharmMFP = true;
        cout << "MFP readout (charm) successful: "  << filename << endl;
      }
    }

    if( ParticlePrototype::N_heavy_flavor > 1 )
    {
      filename = prefix + "m46_" + basename + "_biSectionIter_bottom.dat";
      if ( readTables( filename, bottom ) == MFP_SUCCESS )
      {
        includeBottomMFP = true;
        cout << "MFP readout (bottom) successful: "  << filename << endl;
      }
    }
    
//     displayData();
  }
}



mfp_data::~mfp_data()
{
}



MFP_ERROR_CODE mfp_data::readTables( const string & filename, const FLAVOR_TYPE _F )
{
  fstream file( filename.c_str(), ios::in );

  if ( file )
  {
    string line;
    stringstream inStream;
    double E, lambda_fm, lambda_gev, R22, R23, R32;
    double lambda_therm_fm = -1, lambda_therm_gev = -1;
    double dummy;
    vector<double> lineContent;

    do
    {
      getline( file, line );
      if ( line.find( "#", 0 ) == string::npos && line.length() != 0 )    //ignore empty lines and lines with leading "#"
      {
        inStream.str( "" );
        inStream.clear();
        lineContent.clear();

        inStream.str( line );
        while ( !inStream.eof() )
        {
          inStream >> dummy;
          lineContent.push_back( dummy );
        }

        if ( lineContent.size() >= 3 )
        {
          E = lineContent[0];
          
          lambda_fm = lineContent[1];
          lambda_gev = lineContent[2];
          if( lambda_therm_fm < 0 )
          {
            lambda_therm_fm = lineContent[3];
            lambda_therm_gev = lambda_therm_fm / 0.197;
          }
          R22 = lineContent[6];
          R23 = lineContent[7];
          R32 = lineContent[8];
        }
        else
        {
          string errMsg = "severe error in interpolation23::readTables() - line too short";
          throw eMFP_read_error( errMsg );
        }


        if ( !inStream.fail() )
        {
          if ( _F == gluon )
          {
            data_E_gluon.push_back( E );
            data_lambda_gluon_fm.push_back( lambda_fm );
            data_lambda_gluon_gev.push_back( lambda_gev );
            data_R22_gluon.push_back( R22 );
            data_R23_gluon.push_back( R23 );
            data_R32_gluon.push_back( R32 );
            data_lambda_therm_gluon_fm = lambda_therm_fm;
            data_lambda_therm_gluon_gev = lambda_therm_gev;
          }
          else if ( _F == light_quark )
          {
            data_E_quark.push_back( E );
            data_lambda_quark_fm.push_back( lambda_fm );
            data_lambda_quark_gev.push_back( lambda_gev );
            data_R22_quark.push_back( R22 );
            data_R23_quark.push_back( R23 );
            data_R32_quark.push_back( R32 );
            data_lambda_therm_quark_fm = lambda_therm_fm;
            data_lambda_therm_quark_gev = lambda_therm_gev;
          }
          else if ( _F == charm )
          {
            data_E_charm.push_back( E );
            data_lambda_charm_fm.push_back( lambda_fm );
            data_lambda_charm_gev.push_back( lambda_gev );
            data_R22_charm.push_back( R22 );
            data_R23_charm.push_back( R23 );
            data_R32_charm.push_back( R32 );
            data_lambda_therm_charm_fm = lambda_therm_fm;
            data_lambda_therm_charm_gev = lambda_therm_gev;
          }
          else if ( _F == bottom )
          {
            data_E_bottom.push_back( E );
            data_lambda_bottom_fm.push_back( lambda_fm );
            data_lambda_bottom_gev.push_back( lambda_gev );
            data_R22_bottom.push_back( R22 );
            data_R23_bottom.push_back( R23 );
            data_R32_bottom.push_back( R32 );
            data_lambda_therm_bottom_fm = lambda_therm_fm;
            data_lambda_therm_bottom_gev = lambda_therm_gev;
          }
        }
        else
        {
          string errMsg = "severe error in interpolation23::readTables() - ill-formated data file";
          throw eMFP_read_error( errMsg );
        }
      }
    }
    while ( !file.eof() );
  }
  else
  {
    string errMsg = "error in mfp_data::readTables() - could not open " + filename;
    throw eMFP_read_error( errMsg );
  }

  return MFP_SUCCESS;
}



double mfp_data::getLambdaTherm(FLAVOR_TYPE F) const
{
  if( !read_file )
  {
    string errMsg = "error in mfp_data - Tables not loaded since they are calculated with the original Gunion-Bertsch matrix element.";
    throw eMFP_read_error( errMsg );
  }
  
  if ( F == gluon )
  {
    if ( data_lambda_therm_gluon_gev < 0 )
    {
      string errMsg = "error in mfp_data::readTables() - thermal MFP not read";
      throw eMFP_read_error( errMsg );
    }
    else
    {
      return data_lambda_therm_gluon_gev;
    }
  }
  else if ( F == up || F == down || F == strange || F == anti_up || F == anti_down || F == anti_strange || F == light_quark || F == anti_light_quark )
  {
    if ( data_lambda_therm_quark_gev < 0 )
    {
      string errMsg = "error in mfp_data::readTables() - thermal MFP not read";
      throw eMFP_read_error( errMsg );
    }
    else
    {
      return data_lambda_therm_quark_gev;
    }
  }
  else if ( F == charm || F == anti_charm )
  {
    if ( data_lambda_therm_charm_gev < 0 )
    {
      string errMsg = "error in mfp_data::readTables() - thermal MFP not read";
      throw eMFP_read_error( errMsg );
    }
    else
    {
      return data_lambda_therm_charm_gev;
    }
  }
  else if ( F == bottom || F == anti_bottom )
  {
    if ( data_lambda_therm_bottom_gev < 0 )
    {
      string errMsg = "error in mfp_data::readTables() - thermal MFP not read";
      throw eMFP_read_error( errMsg );
    }
    else
    {
      return data_lambda_therm_bottom_gev;
    }
  }
  else
  {
    string errMsg = "error in mfp_data::readTables() - flavor type could not be resolved";
    throw eMFP_read_error( errMsg );
  }
}


int mfp_data::findIndices( const double _E, const std::vector<double> data_E ) const
{
  const int data_size = data_E.size();
  
  //find j such that E is between data_E[j] and data_E[j+1]
  //bisection-algorithm according to "numerical recipes in c++", chapter 3.4
  int ju, jl, j;
  bool ascend;

  jl = -1;
  ju = data_size;
  ascend = ( data_E[data_size-1] >= data_E[0] );

  while ( ju - jl > 1 )
  {
    int jm = ( ju + jl ) >> 1;  //>> 1 bitwise shift corresponding to division by 2
    if ( ( _E >= data_E[jm] ) == ascend )
      jl = jm;
    else
      ju = jm;
  }
  if ( _E == data_E[0] )
    j = 0;
  else if ( _E == data_E[data_size-1] )
    j = data_size - 2;
  else
    j = jl;
  //-----------------------------------------------------
  
  return j;
}


// in 1/GeV
double mfp_data::getLambda( const double E, const FLAVOR_TYPE F ) const
{
  if( !read_file )
  {
    string errMsg = "error in mfp_data - Tables not loaded since they are calculated with the original Gunion-Bertsch matrix element.";
    throw eMFP_read_error( errMsg );
  }
  
  const std::vector<double> *dataArray = NULL;
  const std::vector<double> *data_E = NULL;
  
  if ( F == gluon )
  {
    dataArray = &data_lambda_gluon_gev;
    data_E = &data_E_gluon;
  }
  else if ( F == up || F == down || F == strange || F == anti_up || F == anti_down || F == anti_strange || F == light_quark || F == anti_light_quark )
  {
    if ( includeQuarkMFP )
    {
      dataArray = &data_lambda_quark_gev;
      data_E = &data_E_quark;
    }
    else
    {
      string errMsg = "error in mfp_data::readTables() - quark MFP requested but not provided by data file";
      throw eMFP_read_error( errMsg );
    }
  }
  else if ( F == charm || F == anti_charm )
  {
    if ( includeCharmMFP )
    {
      dataArray = &data_lambda_charm_gev;
      data_E = &data_E_charm;      
    }
    else
    {
      string errMsg = "error in mfp_data::readTables() - charm MFP requested but not provided by data file";
      throw eMFP_read_error( errMsg );
    }
  }
  else if ( F == bottom || F == anti_bottom )
  {
    if ( includeBottomMFP )
    {
      dataArray = &data_lambda_bottom_gev;
      data_E = &data_E_bottom;      
    }
    else
    {
      string errMsg = "error in mfp_data::readTables() - bottom MFP requested but not provided by data file";
      throw eMFP_read_error( errMsg );
    }
  }
  else
  {
    string errMsg = "error in mfp_data::readTables() - flavor type could not be resolved";
    throw eMFP_read_error( errMsg );
  }
  
  const int j = findIndices( E, *data_E );
  
  //create sub-arrays and interpolate    ----------------
  double result, err;
  int ord = 4;
  double x[ord+1];
  double y[ord+1];
  int k = min( max( j - ( ord - 1 ) / 2, 0 ), static_cast<int>((*data_E).size()) - ord );  //see "numerical recipes", chapter 3.4

  for ( int i = 1;i <= ord;++i )
  {
    x[i] = (*data_E)[k+( i-1 )];
    y[i] = (*dataArray)[k+( i-1 )];
  }
  polint( x, y, ord, E, &result, &err );
  //-----------------------------------------------------

  return result;
}


double mfp_data::getR22(const double E, const FLAVOR_TYPE F) const
{
  if( !read_file )
  {
    string errMsg = "error in mfp_data - Tables not loaded since they are calculated with the original Gunion-Bertsch matrix element.";
    throw eMFP_read_error( errMsg );
  }
  
  const std::vector<double> *dataArray = NULL;
  const std::vector<double> *data_E = NULL;
  
  if ( F == gluon )
  {
    dataArray = &data_R22_gluon;
    data_E = &data_E_gluon;
  }
  else if ( F == up || F == down || F == strange || F == anti_up || F == anti_down || F == anti_strange || F == light_quark || F == anti_light_quark )
  {
    if ( includeQuarkMFP )
    {
      dataArray = &data_R22_quark;
      data_E = &data_E_quark;
    }
    else
    {
      string errMsg = "error in mfp_data::readTables() - quark MFP requested but not provided by data file";
      throw eMFP_read_error( errMsg );
    }
  }
  else if ( F == charm || F == anti_charm )
  {
    if ( includeCharmMFP )
    {
      dataArray = &data_R22_charm;
      data_E = &data_E_charm;      
    }
    else
    {
      string errMsg = "error in mfp_data::readTables() - charm MFP requested but not provided by data file";
      throw eMFP_read_error( errMsg );
    }
  }
  else if ( F == bottom || F == anti_bottom )
  {
    if ( includeBottomMFP )
    {
      dataArray = &data_R22_bottom;
      data_E = &data_E_bottom;      
    }
    else
    {
      string errMsg = "error in mfp_data::readTables() - bottom MFP requested but not provided by data file";
      throw eMFP_read_error( errMsg );
    }
  }
  else
  {
    string errMsg = "error in mfp_data::readTables() - flavor type could not be resolved";
    throw eMFP_read_error( errMsg );
  }
  
  const int j = findIndices( E, *data_E );
  
  //create sub-arrays and interpolate    ----------------
  double result, err;
  int ord = 4;
  double x[ord+1];
  double y[ord+1];
  int k = min( max( j - ( ord - 1 ) / 2, 0 ), static_cast<int>((*data_E).size()) - ord );  //see "numerical recipes", chapter 3.4

  for ( int i = 1;i <= ord;++i )
  {
    x[i] = (*data_E)[k+( i-1 )];
    y[i] = (*dataArray)[k+( i-1 )];
  }
  polint( x, y, ord, E, &result, &err );
  //-----------------------------------------------------
  
  return result;
}



double mfp_data::getR23(const double E, const FLAVOR_TYPE F ) const
{
  if( !read_file )
  {
    string errMsg = "error in mfp_data - Tables not loaded since they are calculated with the original Gunion-Bertsch matrix element.";
    throw eMFP_read_error( errMsg );
  }
  
  const std::vector<double> *dataArray = NULL;
  const std::vector<double> *data_E = NULL;
  
  if ( F == gluon )
  {
    dataArray = &data_R23_gluon;
    data_E = &data_E_gluon;
  }
  else if ( F == up || F == down || F == strange || F == anti_up || F == anti_down || F == anti_strange || F == light_quark || F == anti_light_quark )
  {
    if ( includeQuarkMFP )
    {
      dataArray = &data_R23_quark;
      data_E = &data_E_quark;
    }
    else
    {
      string errMsg = "error in mfp_data::readTables() - quark MFP requested but not provided by data file";
      throw eMFP_read_error( errMsg );
    }
  }
  else if ( F == charm || F == anti_charm )
  {
    if ( includeCharmMFP )
    {
      dataArray = &data_R23_charm;
      data_E = &data_E_charm;      
    }
    else
    {
      string errMsg = "error in mfp_data::readTables() - charm MFP requested but not provided by data file";
      throw eMFP_read_error( errMsg );
    }
  }
  else if ( F == bottom || F == anti_bottom )
  {
    if ( includeBottomMFP )
    {
      dataArray = &data_R23_bottom;
      data_E = &data_E_bottom;      
    }
    else
    {
      string errMsg = "error in mfp_data::readTables() - bottom MFP requested but not provided by data file";
      throw eMFP_read_error( errMsg );
    }
  }
  else
  {
    string errMsg = "error in mfp_data::readTables() - flavor type could not be resolved";
    throw eMFP_read_error( errMsg );
  }
  
  const int j = findIndices( E, *data_E );
  
  //create sub-arrays and interpolate    ----------------
  double result, err;
  int ord = 4;
  double x[ord+1];
  double y[ord+1];
  int k = min( max( j - ( ord - 1 ) / 2, 0 ), static_cast<int>((*data_E).size()) - ord );  //see "numerical recipes", chapter 3.4

  for ( int i = 1;i <= ord;++i )
  {
    x[i] = (*data_E)[k+( i-1 )];
    y[i] = (*dataArray)[k+( i-1 )];
  }
  polint( x, y, ord, E, &result, &err );
  //-----------------------------------------------------
  
  return result;
}

double mfp_data::getR32( const double E, const FLAVOR_TYPE F) const
{
  if( !read_file )
  {
    string errMsg = "error in mfp_data - Tables not loaded since they are calculated with the original Gunion-Bertsch matrix element.";
    throw eMFP_read_error( errMsg );
  }
  
  const std::vector<double> *dataArray = NULL;
  const std::vector<double> *data_E = NULL;
  
  if ( F == gluon )
  {
    dataArray = &data_R32_gluon;
    data_E = &data_E_gluon;
  }
  else if ( F == up || F == down || F == strange || F == anti_up || F == anti_down || F == anti_strange || F == light_quark || F == anti_light_quark )
  {
    if ( includeQuarkMFP )
    {
      dataArray = &data_R32_quark;
      data_E = &data_E_quark;
    }
    else
    {
      string errMsg = "error in mfp_data::readTables() - quark MFP requested but not provided by data file";
      throw eMFP_read_error( errMsg );
    }
  }
  else if ( F == charm || F == anti_charm )
  {
//     if ( includeCharmMFP )
//     {
//       dataArray = &data_R32_charm;
//       data_E = &data_E_charm;      
//     }
//     else
//     {
//       string errMsg = "error in mfp_data::readTables() - charm MFP requested but not provided by data file";
//       throw eMFP_read_error( errMsg );
//     }
    string errMsg = "error in mfp_data::readTables() - charm R32 requested but not provided by data file";
    throw eMFP_read_error( errMsg );
  }
  else if ( F == bottom || F == anti_bottom )
  {
//     if ( includeBottomMFP )
//     {
//       dataArray = &data_R32_bottom;
//       data_E = &data_E_bottom;      
//     }
//     else
//     {
//       string errMsg = "error in mfp_data::readTables() - bottom MFP requested but not provided by data file";
//       throw eMFP_read_error( errMsg );
//     }
    string errMsg = "error in mfp_data::readTables() - bottom R32 requested but not provided by data file";
    throw eMFP_read_error( errMsg );
  }
  else
  {
    string errMsg = "error in mfp_data::readTables() - flavor type could not be resolved";
    throw eMFP_read_error( errMsg );
  }
  
  const int j = findIndices( E, *data_E );
  
  //create sub-arrays and interpolate    ----------------
  double result, err;
  int ord = 4;
  double x[ord+1];
  double y[ord+1];
  int k = min( max( j - ( ord - 1 ) / 2, 0 ), static_cast<int>((*data_E).size()) - ord );  //see "numerical recipes", chapter 3.4

  for ( int i = 1;i <= ord;++i )
  {
    x[i] = (*data_E)[k+( i-1 )];
    y[i] = (*dataArray)[k+( i-1 )];
  }
  polint( x, y, ord, E, &result, &err );
  //-----------------------------------------------------
  
  return result;
}


void mfp_data::polint(const double xa[], const double ya[], const int n, const double x,double *y,double *dy) const
{
  int ns=1;
  double den,dif,ho,hp,w;
  double *c,*d;

  dif=fabs(x-xa[1]);

  c=new double[n+1];
  d=new double[n+1];

  for(int i=1;i<=n;i++)
  {
    double dift;
    if((dift=fabs(x-xa[i])) < dif)
    {
      ns=i;
      dif=dift;
    }
    c[i]=ya[i];
    d[i]=ya[i];
  }

  *y=ya[ns--];

  for(int m=1;m<n;m++)
  {
    for(int i=1;i<=n-m;i++)
    {
      ho=xa[i]-x;
      hp=xa[i+m]-x;
      w=c[i+1]-d[i];
      if ((den=ho-hp) == 0.0)
        cout << "error in routine polint" << endl;
      den=w/den;
      d[i]=hp*den;
      c[i]=ho*den;
    }

    *y += (*dy=(2*ns <(n-m) ? c[ns+1] : d[ns--]));
  }

  delete[] c;
  delete[] d;
}


void mfp_data::displayData() const
{
  if( !read_file )
  {
    string errMsg = "error in mfp_data - Tables not loaded since they are calculated with the original Gunion-Bertsch matrix element.";
    throw eMFP_read_error( errMsg );
  }
 
  cout << "number of entries (gluon): " << data_E_gluon.size() << endl;
  for (int i=0; i < data_E_gluon.size(); i++)
  {
    cout << data_E_gluon[i] << "\t" << data_lambda_gluon_fm[i] << "\t" << data_lambda_gluon_gev[i] << endl; 
  }
  cout << endl << endl;
  
  if( includeQuarkMFP )
  {
    cout << "number of entries (quark): " << data_E_quark.size() << endl;
    for (int i=0; i < data_E_quark.size(); i++)
    {
      cout << data_E_quark[i] << "\t" << data_lambda_quark_fm[i] << "\t" << data_lambda_quark_gev[i] << endl; 
    }
    cout << endl << endl;
  }

  if( includeCharmMFP )
  {
    cout << "number of entries (charm): " << data_E_charm.size() << endl;
    for (int i=0; i < data_E_charm.size(); i++)
    {
      cout << data_E_charm[i] << "\t" << data_lambda_charm_fm[i] << "\t" << data_lambda_charm_gev[i] << endl; 
    }
    cout << endl << endl;
  }

  if( includeBottomMFP )
  {
    cout << "number of entries (bottom): " << data_E_bottom.size() << endl;
    for (int i=0; i < data_E_bottom.size(); i++)
    {
      cout << data_E_bottom[i] << "\t" << data_lambda_bottom_fm[i] << "\t" << data_lambda_bottom_gev[i] << endl; 
    }
    cout << endl << endl;
  }
}


mfpForHeavyIonCollision::mfpForHeavyIonCollision( configBase *const _c ) : 
  fitParameterGluon(-0.569485),
  fitParameterQuark(-1.07081),
  data_loaded( false ),
  theConfig( _c )
{
  //just to be sure
  temperaturesForMfpData.clear();
  mfpData.clear();
}




mfpForHeavyIonCollision::~mfpForHeavyIonCollision()
{
  temperaturesForMfpData.clear();
  mfpData.clear();
}


void mfpForHeavyIonCollision::loadData()
{
  //write temperatures for which the mfp data should be read from tables into a vector
  // T in GeV
  temperaturesForMfpData.push_back( 0.1 );
  temperaturesForMfpData.push_back( 0.2 );
  temperaturesForMfpData.push_back( 0.3 );
  temperaturesForMfpData.push_back( 0.4 );
  temperaturesForMfpData.push_back( 0.5 );
  temperaturesForMfpData.push_back( 0.6 );
  temperaturesForMfpData.push_back( 0.8 );
  temperaturesForMfpData.push_back( 1.0 );
  temperaturesForMfpData.push_back( 1.5 );
  temperaturesForMfpData.push_back( 2.0 );
  
  //read temperature from the above created vector, initialize corresponding mfp_data object and push it into vector
  for ( unsigned int i = 0; i < temperaturesForMfpData.size(); i++ )
  {
    mfp_data tempMFP( temperaturesForMfpData[i], true );
    
    mfpData.push_back( tempMFP );
  }
  
  fugacityDependenceEvaluatedAtFugacityOne = fugacityDependence( 1, 1 );
  
  data_loaded = true;
}


double mfpForHeavyIonCollision::getMeanFreePath(const double _E, const FLAVOR_TYPE _F, const double _T, const double _ngTest, const double _nqTest, const UNIT_TYPE _unit ) const
{
  if( !data_loaded )
  {
    string errMsg = "MFP data not loaded yet.";
    throw eMFP_read_error( errMsg );
  }
  
  const double gluonFugacity = _ngTest / ( getGluonDensity( _T ) );
  double quarkFugacity;
  if( ParticlePrototype::N_light_flavor != 0 )
    quarkFugacity = _nqTest / ( getQuarkDensity( _T ) );
  else
    quarkFugacity = 0.0;
  
  const double scaleForFugacity = fugacityDependence( gluonFugacity, quarkFugacity ) / fugacityDependenceEvaluatedAtFugacityOne;
  
  const int nInterpolValues = 4;
  int startIndex = getStartIndexForTemperatureValues( _T, nInterpolValues );
  
  double xa[nInterpolValues+1], ya[nInterpolValues+1];
  
  for ( int i = 0; i < nInterpolValues; i++ )
  {
    xa[i+1] = temperaturesForMfpData[ startIndex + i ];
    ya[i+1] = mfpData[ startIndex + i ].getLambda( _E, _F );  //1/GeV
  }
  
  double interpolResult = 0, interpolError = 0;
  polint( xa, ya, nInterpolValues, _T, &interpolResult, &interpolError );
  
  double conversionFactor = 1;
  switch( _unit )
  {
    case GeV:
      conversionFactor = 1;
      break;
    case fm:
      conversionFactor = 0.197;
      break;
  }
  interpolResult *= conversionFactor;
  
  interpolResult *= scaleForFugacity;
  
  return interpolResult;
}



double mfpForHeavyIonCollision::fugacityDependence( const double _gluonFugacity, const double _quarkFugacity ) const
{
  double termGluon = sqrt( 2 * _gluonFugacity - fitParameterGluon * pow( _gluonFugacity, 2 ) );
  double termQuark = sqrt( 2 * _quarkFugacity - fitParameterQuark * pow( _quarkFugacity, 2 ) );
  return ( 1 / ( termGluon + termQuark ) );
}




double mfpForHeavyIonCollision::getGluonDensity( const double T, const double fugacity ) const
{
  return fugacity * 16 * pow( T, 3.0 ) / pow( M_PI, 2.0 ) / pow( 0.197, 3.0 );
}



double mfpForHeavyIonCollision::getQuarkDensity( const double T, const double fugacity ) const
{
  return fugacity * 12 * ParticlePrototype::N_light_flavor * pow( T, 3.0 ) / pow( M_PI, 2.0 ) / pow( 0.197, 3.0 );
}



int mfpForHeavyIonCollision::getStartIndexForTemperatureValues(const double _T, const unsigned int _nValues) const
{
  if ( _nValues > temperaturesForMfpData.size() )
  {
    string errMsg = "requested number of interpolation points larger than temperature array size";
    throw eMFP_read_error( errMsg );
  }
  
  unsigned int nn = 0;
  while ( temperaturesForMfpData[nn] < _T && nn < temperaturesForMfpData.size() )
  {
    ++nn;
  }
  
  nn = nn - _nValues / 2;
  if ( nn < 0 )
  {
    nn = 0;
  }
  else if ( nn + _nValues > temperaturesForMfpData.size() ) 
  {
    while ( nn + _nValues > temperaturesForMfpData.size() )
    {
      --nn;
    }
  }

  return nn;
}




/**
* Routine that linearly interpolates. Attention: Offset = 1 !
*
* @param[in] xa[] array of x-values
* @param[in] ya[] array of y-values
* @param[in] n number of x- and y-values used for the interpolation
* @param[in] x value at which the interpolated y should be evaluated
* @param[out] y result
* @param[out] dy some sort of error
*/
void mfpForHeavyIonCollision::polint(const double xa[], const double ya[], const int n, const double x, double* y, double* dy) const
{
  
  int ns = 1;
  double den, dif, dift, ho, hp, w;
  double *c, *d;
  
  dif = fabs( x - xa[1] );
  
  c = new double[n+1];
  d = new double[n+1];
  
  for ( int i = 1;i <= n;i++ )
  {
    if (( dift = fabs( x - xa[i] ) ) < dif )
    {
      ns = i;
      dif = dift;
    }
    c[i] = ya[i];
    d[i] = ya[i];
  }
  
  *y = ya[ns--];
  
  for ( int m = 1;m < n;m++ )
  {
    for ( int i = 1;i <= n - m;i++ )
    {
      ho = xa[i] - x;
      hp = xa[i+m] - x;
      w = c[i+1] - d[i];
      if (( den = ho - hp ) == 0.0 )
      {
        cout << "Error in routine polint" << endl;
        cout << x << endl;
        cout << xa[1] << "\t" << ya[1] << endl;
        cout << xa[2] << "\t" << ya[2] << endl;
        cout << xa[3] << "\t" << ya[3] << endl;
        cout << xa[4] << "\t" << ya[4] << endl;
      }
      den = w / den;
      d[i] = hp * den;
      c[i] = ho * den;
    }
    
    *y += ( *dy = ( 2 * ns < ( n - m ) ? c[ns+1] : d[ns--] ) );
  }
  
  delete[] c;
  delete[] d;
}


