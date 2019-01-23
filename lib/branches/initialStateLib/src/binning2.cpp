//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/binning2.cpp $
//$LastChangedDate: 2017-04-06 14:55:30 +0200 (Do, 06. Apr 2017) $
//$LastChangedRevision: 2556 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


#include <math.h>
#include <iostream>
#include "binning2.h"

using namespace std;

//----------------------------------------------------------------------------

binning2d::binning2d( ) :
  x_min(1e20),
  x_max(-1e20),
  y_min(1e20),
  y_max(-1e20),
  binWidth_x(),
  binWidth_y(),
  nBins_x(100),           // default number of bins
  nBins_y(100),           // default number of bins
  nEntries(0)       // number of collected events=0
{
}

binning2d::binning2d( const string & name, 
                      const double min_x_arg, const double max_x_arg, const int n_x, 
                      const double min_y_arg, const double max_y_arg, const int n_y ) :
  x_min(min_x_arg),
  x_max(max_x_arg),
  y_min(min_y_arg),
  y_max(max_y_arg),
  binWidth_x( ( x_max - x_min ) / n_x ),
  binWidth_y( ( y_max - y_min ) / n_y ),
  nBins_x(n_x),
  nBins_y(n_y),
  nEntries(0)       // number of collected events=0
{
  // initialize 2 dimensional array with 0
  for ( int i = 0; i < nBins_x; i++ )
  {
    theBins.push_back ( vector<int>() );
    for ( int j = 0; j < nBins_y; j++ )
    {
      theBins[i].push_back ( 0 );
    }
  }

  setFilename ( name );
}

binning2d::binning2d( const string & name, 
                      const double min_x_arg, const double max_x_arg, const double width_x, 
                      const double min_y_arg, const double max_y_arg, const double width_y ) :
  x_min(min_x_arg),
  x_max(max_x_arg),
  y_min(min_y_arg),
  y_max(max_y_arg),
  binWidth_x(width_x),
  binWidth_y(width_y),
  nBins_x( ( x_max - x_min ) / binWidth_x ),
  nBins_y( ( y_max - y_min ) / binWidth_y ),
  nEntries(0)       // number of collected events=0
{
  // initialize 2 dimensional array with 0
  for ( int i = 0; i < nBins_x; i++ )
  {
    theBins.push_back ( vector<int>() );
    for ( int j = 0; j < nBins_y; j++ )
    {
      theBins[i].push_back ( 0 );
    }
  }

  setFilename ( name );
}


void binning2d::setMinMaxN( const double min_x_arg, const double max_x_arg, const int n_x, 
                            const double min_y_arg, const double max_y_arg, const int n_y )
{
  x_min = min_x_arg;
  x_max = max_x_arg;
  nBins_x = n_x;
  y_min = min_y_arg;
  y_max = max_y_arg;
  nBins_y = n_y;
  
  // initialize 2 dimensional array with 0
  for ( int i = 0; i < nBins_x; i++ ) 
  {
    theBins.push_back ( vector<int>() );

    for ( int j = 0; j < nBins_y; j++ )
    {
      theBins[i].push_back ( 0 );
    }
  }
  
  binWidth_x = (x_max - x_min) / nBins_x;
  binWidth_y = (y_max - y_min) / nBins_y;
  nEntries = 0; 
}

void binning2d::setMinMaxWidth( const double min_x_arg, const double max_x_arg, const double width_x, 
                                const double min_y_arg, const double max_y_arg, const double width_y )
{
  x_min = min_x_arg;
  x_max = max_x_arg;
  binWidth_x = width_x;

  y_min = min_y_arg;
  y_max = max_y_arg;
  binWidth_y = width_y;

  nBins_x = int ( ( x_max - x_min ) / binWidth_x );
  nBins_y = int ( ( y_max - y_min ) / binWidth_y );

  // initialize 2 dimensional array with 0
  for ( int i = 0; i < nBins_x; i++ )
  {
    theBins.push_back ( vector<int>() );

    for ( int j = 0; j < nBins_y; j++ )
    {
      theBins[i].push_back ( 0 );
    }
  }

  nEntries = 0;
}


void binning2d::add(const double x, const double y)
{
  int index_x = int((x - x_min)/binWidth_x);
  if (index_x >= nBins_x)
    index_x = nBins_x-1;
  else if (index_x < 0)
    index_x = 0;

  int index_y = int((y - y_min)/binWidth_y);
  if (index_y >= nBins_y)
    index_y = nBins_y-1;
  else if (index_y < 0)
    index_y = 0;

  ++theBins[index_x][index_y];
  ++nEntries;
}

vector<int> binning2d::getIndex( const double x, const double y )
{
  vector<int> indices( 2, 0 );
  
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
  
  indices[0] = index_x;
  indices[1] = index_y;
  
  return indices;
}

void binning2d::print()
{
  vector<int> sumX ( nBins_x, 0 );
  for ( int i = 0; i < nBins_x; i++ )
  {
    for ( int j = 0; j < nBins_y; j++ )
    {
      sumX[i] += theBins[i][j];
    }
  }

  vector<int> sumY ( nBins_y, 0 );
  for ( int j = 0; j < nBins_y; j++ )
  {
    for ( int i = 0; i < nBins_x; i++ )
    {
      sumY[j] += theBins[i][j];
    }
  }

  fstream file ( filename.c_str(), ios::out );
  file << "#bin width_x: " << binWidth_x << "   bin width_y: " << binWidth_y <<  "   number of binned events: " << nEntries << endl;
  for (int i=0; i<nBins_x; i++)
  {
    for (int j=0; j<nBins_y; j++)
    {
      file.width(15);
      file << x_min + (i+0.5)*binWidth_x;
      file.width(15);
      file << y_min + (j+0.5)*binWidth_y;
      file.width(15);
      file << theBins[i][j];
      file.width ( 15 );
      //      file << double(theBins[i][j])/sumX[i]/sumY[j]/nEntries/binWidth_x/binWidth_y << endl;
      file << double ( theBins[i][j] ) / nEntries / binWidth_x / binWidth_y << endl;
    }
    file << endl;
  }   

  file<<endl;

  for ( int i = 0; i < nBins_x; i++ )
  {
    file.width(15);
    file << x_min + (i+0.5)*binWidth_x;
    file.width(15);
    file << sumX[i] / ( nEntries * binWidth_x );
    file << endl;
  }
  file << endl;
  file << endl;

  for ( int i = 0; i < nBins_y; i++ )
  {
    file.width(15);
    file << y_min + (i+0.5)*binWidth_y;
    file.width(15);
    file << sumY[i] / ( nEntries * binWidth_y );
    file << endl;
  }
  file << endl;
  file << endl;

  file.close();
}

/**
 * This is mainly needed for OpenMP reduction purposes. 
 */
binning2d & binning2d::operator += (const binning2d & q)
{
  for (int i=0; i<nBins_x; i++)
  {
    for (int j=0; j<nBins_y; j++)
    {
      theBins[i][j] += q.theBins[i][j];
    }
  }
  nEntries += q.nEntries;
  return *this;
}


//----------------------------------------------------------------------------

binningValues2d::binningValues2d()
{
}

binningValues2d::binningValues2d( const string & name, const double min_x_arg, const double max_x_arg, const int n_x, const double min_y_arg, const double max_y_arg, const int n_y )
{
  x_min = min_x_arg;
  x_max = max_x_arg;
  nBins_x = n_x;
  y_min = min_y_arg;
  y_max = max_y_arg;
  nBins_y = n_y;
  
  // initialize 2 dimensional array with 0
  for ( int i = 0; i < nBins_x; i++ ) 
  {
    theBins.push_back ( vector<double>() );
    theBinsN.push_back ( vector<int>() );

    for ( int j = 0; j < nBins_y; j++ )
    {
      theBins[i].push_back ( 0.0 );
      theBinsN[i].push_back ( 0 );
    }
  }
  
  binWidth_x = (x_max - x_min) / nBins_x;
  binWidth_y = (y_max - y_min) / nBins_y;
  nEntries = 0; 
  setFilename(name);
}

void binningValues2d::setMinMaxN( const double min_x_arg, const double max_x_arg, const int n )
{
  x_min = y_min = min_x_arg;
  x_max = y_max = max_x_arg;
  nBins_x = nBins_y = n;

  // initialize 2 dimensional array with 0
  for ( int i = 0; i < nBins_x; i++ )
  {
    theBins.push_back ( vector<double>() );
    theBinsN.push_back ( vector<int>() );

    for ( int j = 0; j < nBins_y; j++ )
    {
      theBins[i].push_back ( 0.0 );
      theBinsN[i].push_back ( 0 );
    }
  }

  binWidth_x = ( x_max - x_min ) / nBins_x;
  binWidth_y = ( y_max - y_min ) / nBins_y;
  nEntries = 0;
}


void binningValues2d::setMinMaxN( const double min_x_arg, const double max_x_arg, const int n_x, const double min_y_arg, const double max_y_arg, const int n_y )
{
  x_min = min_x_arg;
  x_max = max_x_arg;
  nBins_x = n_x;

  y_min = min_y_arg;
  y_max = max_y_arg;
  nBins_y = n_y;

  // initialize 2 dimensional array with 0
  for ( int i = 0; i < nBins_x; i++ )
  {
    theBins.push_back ( vector<double>() );
    theBinsN.push_back ( vector<int>() );

    for ( int j = 0; j < nBins_y; j++ )
    {
      theBins[i].push_back ( 0.0 );
      theBinsN[i].push_back ( 0 );
    }
  }

  binWidth_x = ( x_max - x_min ) / nBins_x;
  binWidth_y = ( y_max - y_min ) / nBins_y;
  nEntries = 0;
}

void binningValues2d::setMinMaxWidth( const double min_x_arg, const double max_x_arg, const double width_x, const double min_y_arg, const double max_y_arg, const double width_y)
{
  x_min = min_x_arg;
  x_max = max_x_arg;
  binWidth_x = width_x;

  y_min = min_y_arg;
  y_max = max_y_arg;
  binWidth_y = width_y;

  nBins_x = int ( ( x_max - x_min ) / binWidth_x );
  nBins_y = int ( ( y_max - y_min ) / binWidth_y );

  // initialize 2 dimensional array with 0
  for ( int i = 0; i < nBins_x; i++ )
  {
    theBins.push_back ( vector<double>() );
    theBinsN.push_back ( vector<int>() );

    for ( int j = 0; j < nBins_y; j++ )
    {
      theBins[i].push_back ( 0.0 );
      theBinsN[i].push_back ( 0 );
    }
  }

  nEntries = 0;
}


void binningValues2d::add( const double x, const double y, const double v)
{
  
  int index_x = int((x - x_min)/binWidth_x);
  if (index_x >= nBins_x)
    index_x = nBins_x-1;
  else if (index_x < 0)
    index_x = 0;
  int index_y = int((y - y_min)/binWidth_y);
  if (index_y >= nBins_y)
    index_y = nBins_y-1;
  else if (index_y < 0)
    index_y = 0;
    
  theBins[index_x][index_y] += v;
  ++theBinsN[index_x][index_y];
  
  ++nEntries;
}


//sort into bins if necessary and print results to file
void binningValues2d::print()
{
   
  fstream file(filename.c_str(),ios::out);
  file << "#bin width_x: " << binWidth_x << "   bin width_y: " << binWidth_y <<  "   number of binned events: " << nEntries << endl;
  for (int i=0; i<nBins_x; i++)
  {
    for (int j=0; j<nBins_y; j++)
    {
      file.width(15);
      file << x_min + (i+0.5)*binWidth_x;
      file.width(15);
      file << y_min + (j+0.5)*binWidth_y;
      file.width(15);
      if(theBinsN[i][j] != 0)
        file << theBins[i][j];
      else
        file << 0;
      file.width(15);
      file << theBinsN[i][j] << endl;
    }
    file << endl;
  }   
  file.close();
}

//sort into bins if necessary and print results to file
// Print in slightly other format than print(). Useful for 2D image plotting  
void binningValues2d::print2DPlot()
{
   
  fstream file(filename.c_str(),ios::out);
  file << "#bin width_x: " << binWidth_x << "   bin width_y: " << binWidth_y <<  "   number of binned events: " << nEntries << endl;
  for (int i=0; i<nBins_x; i++)
  {
    for (int j=0; j<nBins_y; j++)
    {      
      file << x_min + (i+0.5)*binWidth_x;
      file << "\t";
      file << y_min + (j+0.5)*binWidth_y;
      file << "\t";
      if(theBinsN[i][j] != 0)
        file << theBins[i][j];
      else
        file << 0;
      file << endl;
    }
  }   
  file.close();
}

double binningValues2d::getBinLabelX ( const int i )
{
  if ( i >= 0 && i < nBins_x )
    return x_min + ( i + 0.5 ) * binWidth_x;
  else
    return 0;
}


double binningValues2d::getBinLabelY ( const int i )
{
  if ( i >= 0 && i < nBins_y )
    return y_min + ( i + 0.5 ) * binWidth_y;
  else
    return 0;
}


double binningValues2d::getBin ( const int i, const int j )
{
  if ( i >= 0 && i < nBins_x && j >= 0 && j < nBins_y )
  {
    if(theBinsN[i][j] != 0)
      return ( theBins[i][j] );
    else
      return 0.;     
  }
  else
    return 0.;
}

binningValues2d::~binningValues2d()
{
}

//----------------------------------------------------------------------------

