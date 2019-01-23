//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/binning3.cpp $
//$LastChangedDate: 2016-04-19 20:13:24 +0200 (Di, 19. Apr 2016) $
//$LastChangedRevision: 2328 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


#include <math.h>

#include "binning3.h"

using namespace std;

//----------------------------------------------------------------------------

binning3d::binning3d() :
  x_min(1e20),
  x_max(-1e20),
  y_min(1e20),
  y_max(-1e20),
  z_min(1e20),
  z_max(-1e20),
  binWidth_x(),
  binWidth_y(),
  binWidth_z(),
  nBins_x(10),           // default number of bins
  nBins_y(10),           // default number of bins
  nBins_z(10),           // default number of bins
  nEntries(0)       // number of collected events=0
{
}

// binning2d::binning2d ( const string name, const double min_x_arg, const double max_x_arg, const int n_x, const double min_y_arg, const double max_y_arg, const int n_y )
// {
//   x_min = min_x_arg;
//   x_max = max_x_arg;
//   nBins_x = n_x;
//   y_min = min_y_arg;
//   y_max = max_y_arg;
//   nBins_y = n_y;
// 
//   // initialize 2 dimensional array with 0
//   for ( int i = 0; i < nBins_x; i++ )
//   {
//     theBins.push_back ( vector<int>() );
// 
//     for ( int j = 0; j < nBins_y; j++ )
//       theBins[i].push_back ( 0 );
//   }
// 
//   binWidth_x = ( x_max - x_min ) / nBins_x;
//   binWidth_y = ( y_max - y_min ) / nBins_y;
//   btype = fixed_min_max_number;
//   nEntries = 0;
//   finished = true;
//   setFilename ( name );
// }

void binning3d::setMinMaxN( const double min_x_arg, const double max_x_arg, const int n_x, const double min_y_arg, const double max_y_arg, const int n_y, const double min_z_arg, const double max_z_arg, const int n_z )
{
  x_min = min_x_arg;
  x_max = max_x_arg;
  nBins_x = n_x;
  y_min = min_y_arg;
  y_max = max_y_arg;
  nBins_y = n_y;
  z_min = min_z_arg;
  z_max = max_z_arg;
  nBins_z = n_z;

  // initialize 2 dimensional array with 0
  for ( int i = 0; i < nBins_x; i++ )
  {
    theBins.push_back ( vector< vector<int> >() );
    for ( int j = 0; j < nBins_y; j++ )
    {
      theBins[i].push_back ( vector<int>() );
      for ( int k = 0; k < nBins_z; k++ )
      {
        theBins[i][j].push_back( 0 );
      }
    }
  }
  binWidth_x = ( x_max - x_min ) / nBins_x;
  binWidth_y = ( y_max - y_min ) / nBins_y;
  binWidth_z = ( z_max - z_min ) / nBins_z;
  nEntries = 0;
}

void binning3d::setMinMaxWidth( const double min_x_arg, const double max_x_arg, const double width_x, const double min_y_arg, const double max_y_arg, const double width_y, const double min_z_arg, const double max_z_arg, const double width_z )
{
  x_min = min_x_arg;
  x_max = max_x_arg;
  binWidth_x = width_x;
  
  y_min = min_y_arg;
  y_max = max_y_arg;
  binWidth_y = width_y;
  
  z_min = min_z_arg;
  z_max = max_z_arg;
  binWidth_z = width_z;

  nBins_x = int ( ( x_max - x_min ) / binWidth_x );
  nBins_y = int ( ( y_max - y_min ) / binWidth_y );
  nBins_z = int ( ( z_max - z_min ) / binWidth_z );
  
  // initialize 2 dimensional array with 0
  for ( int i = 0; i < nBins_x; i++ )
  {
    theBins.push_back ( vector< vector<int> >() );
    for ( int j = 0; j < nBins_y; j++ )
    {
      theBins[i].push_back ( vector<int>() );
      for ( int k = 0; k < nBins_z; k++ )
      {
        theBins[i][j].push_back( 0 );
      }
    }
  }
  
  nEntries = 0;
}



//add element x to the binned data
void binning3d::add( const double x, const double y, const double z )
{

  int index_x = int ( ( x - x_min ) / binWidth_x );
  if ( index_x >= nBins_x )
    index_x = nBins_x - 1;
  else if ( index_x < 0 )
    index_x = 0;
  int index_y = int ( ( y - y_min ) / binWidth_y );
  if ( index_y >= nBins_y )
    index_y = nBins_y - 1;
  else if ( index_y < 0 )
    index_y = 0;
  int index_z = int ( ( z - z_min ) / binWidth_z );
  if ( index_z >= nBins_z )
    index_z = nBins_z - 1;
  else if ( index_z < 0 )
    index_z = 0;
  ++theBins[index_x][index_y][index_z];

  ++nEntries;
}

vector<int> binning3d::getIndex( const double x, const double y, const double z )
{
  vector<int> indices( 3, 0 );
  
  int index_x = int ( ( x - x_min ) / binWidth_x );
  if ( index_x >= nBins_x )
    index_x = nBins_x - 1;
  else if ( index_x < 0 )
    index_x = 0;
  int index_y = int ( ( y - y_min ) / binWidth_y );
  if ( index_y >= nBins_y )
    index_y = nBins_y - 1;
  else if ( index_y < 0 )
    index_y = 0;
  int index_z = int ( ( z - z_min ) / binWidth_z );
  if ( index_z >= nBins_z )
    index_z = nBins_z - 1;
  else if ( index_z < 0 )
    index_z = 0;
  
  indices[0] = index_x;
  indices[1] = index_y;
  indices[2] = index_z;
  
  return indices;
}

//sort into bins if necessary and print results to file
// void binning3d::print()
// {
//   vector<int> sumX ( nBins_x, 0 );
//   for ( int i = 0; i < nBins_x; i++ )
//   {
//     for ( int j = 0; j < nBins_y; j++ )
//     {
//       sumX[i] += theBins[i][j];
//     }
//   }
// 
//   vector<int> sumY ( nBins_y, 0 );
//   for ( int j = 0; j < nBins_y; j++ )
//   {
//     for ( int i = 0; i < nBins_x; i++ )
//     {
//       sumY[j] += theBins[i][j];
//     }
//   }
// 
//   vector<int> sumZ ( nBins_z, 0 );
//   for ( int k = 0; k < nBins_y; j++ )
//   {
//     for ( int i = 0; i < nBins_x; i++ )
//     {
//       sumY[j] += theBins[i][j];
//     }
//   }
// 
//   fstream file ( filename.c_str(), ios::out );
//   file << "#bin width_x: " << binWidth_x << "   bin width_y: " << binWidth_y <<  "   number of binned events: " << nEntries << endl;
//   for ( int i = 0; i < nBins_x; i++ )
//   {
//     for ( int j = 0; j < nBins_y; j++ )
//     {
//       file.width ( 15 );
//       file << x_min + ( i + 0.5 ) *binWidth_x;
//       file.width ( 15 );
//       file << y_min + ( j + 0.5 ) *binWidth_y;
//       file.width ( 15 );
//       file << theBins[i][j];
//       file.width ( 15 );
// //      file << double(theBins[i][j])/sumX[i]/sumY[j]/nEntries/binWidth_x/binWidth_y << endl;
//       file << double ( theBins[i][j] ) / nEntries / binWidth_x / binWidth_y << endl;
//     }
//     file << endl;
//   }
//   file.close();
// }

//----------------------------------------------------------------------------
