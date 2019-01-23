//------------------------------------------------------------------------------------------
//provided by subversion
//------------------------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/interpolation_n_dimensions.cpp $
//$LastChangedDate: 2013-10-16 21:44:56 +0200 (Mi, 16. Okt 2013) $
//$LastChangedRevision: 1488 $
//$LastChangedBy: gallmei $
//------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <math.h>

#include "interpolation_n_dimensions.h"

using namespace std;


interpolation_n_D::~interpolation_n_D()
{
  data.clear();
}


DATA_ERROR_CODE interpolation_n_D::readTables()
{
  fstream file( filename.c_str(), ios::in );

  if ( file )
  {
    string line;
    stringstream inStream;
    double content;

    do
    {
      getline( file, line );
      if ( line.find( "#", 0 ) == string::npos && line.length() != 0 )    //ignore empty lines and lines with leading "#"
      {
        if ( line.find( "-inf", 0 ) == string::npos )  //the string "-inf" can't be pushed into a vector<double>
        {
          inStream.str( "" );
          inStream.clear();

          inStream.str( line );
          inStream >> content;

          if ( !inStream.fail() )
          {
            data.push_back( content );
          }
          else
          {
            cout << "severe error in interpolation1d::readTables() - ill-formated data file" << endl;
            return DATA_READOUT_ERROR;
          }
        }
        else     //i.e. when string "-inf" is found
        {
          cout << "error  in interpolation1d::readTables() - inf entries" << endl;
          return DATA_READOUT_ERROR;
        }

      }
    }
    while ( !file.eof() );
  }
  else
  {
    // throw an exception if the file could not have been read
    throw eData_read_error( "error in interpolation1d::readTables() - could not open " + filename );
    return DATA_FILE_ERROR;
  }

  // The number of tabulated values read from file must be equal to the expected number, otherwise
  // the interpolation routine will fail to work correctly!
  if ( data.size() != expected_entries )
  {
    cout << data.size() << " instead of " << expected_entries << "entries" << endl;
    cout << "Filename: " << filename << endl;
    // throw an exception if the number of entries in the file is not right
    throw eData_read_error( "error in interpolation_n_D::readTables() - unexpected number of entries in data table" );
    return DATA_SIZE_ERROR;
  }
  else
    return DATA_SUCCESS;
}



/**
 * Configure 1-dimensional interpolation
 *
 * @param[in] n_i_arg total numbers of tabulated points on the a "axes"
 * @param[in] delta_a_arg spacings of tabulated values in a-direction
 * @param[in] a_start_arg the smallest tabulated value of a
 * @return data, the value of data at (a) found via interpolation from the tabulated values of f
 */
void interpolation1d::configure( const int n_i_arg, const double delta_a_arg, const double a_start_arg, string filename_arg )
{
  // total numbers of tabulated points on the a "axes"
  n_i = n_i_arg;
  
  // spacings of tabulated values in a-direction
  delta_a = delta_a_arg;

  // the smallest tabulated value of a,
  a_start = a_start_arg;
  
  filename = filename_arg;
  
  loadData();
}


/**
 * Reserve memory and load the data.
 */
void interpolation1d::loadData()
{
  expected_entries = n_i;
  //allocate memory based on the known size (# of entries) of the table
  //a matter of efficieny - not necessary
  data.reserve( expected_entries + 1 );

  if ( readTables() == DATA_SUCCESS )
  {
    cout << "data readout successful: " << filename << endl;
  }
}


/**
 * 1-dimensional interpolation / extrapolation of the tabulated values at a given point (a)
 * This routine serves as a wrapper to different methods of interpolation (linear, polynomial, cubic splines). 
 *
 * @param[in] a point first dimension
 * @return data, the value of data at (a) found via interpolation from the tabulated values of f
 */
double interpolation1d::getInterpolatedData( const double a) const
{
  return linearInterpolation( a );
}


/**
 * 1-dimensional interpolation / extrapolation of the tabulated values at a given point (a)
 * This routine is closely tailored to the known structure of the tabulated data (stored in interpolation1d::data).
 * See below for a short documentation of the required data structure.
 *
 * @param[in] a point first dimension
 * @return data, the value of data at (a) found via interpolation from the tabulated values of f
 */
double interpolation1d::linearInterpolation( const double a ) const
{
  // find indices i such that
  // a[i] <= a < a[i+1]
  const int i = int(( a - a_start ) / delta_a ); 

  if ( i < 0)
  {
    cout << "unexpected - i not in range" << endl;
    cout << n_i << "  " << endl;
    cout << i << "  " << endl;
    cout << a << "  " <<endl;
    cout << a_start << "  " << endl;

    return 0;
  }

  // temporary vector objects used to initialize the multi-dimensional arrays based on STL vectors
  vector<double> f( 2, 0 );

  int delta_ex_a;
  int index_i;
  
  if ( i < (n_i - 1) ) // interpolated linearly in a direction
  {
    index_i = i+1;
    delta_ex_a = 1;
  }
  else
  {
    index_i = (n_i - 1);
    delta_ex_a = 2;
  }
  
  for ( int p = 0; p <= 1; p++ )
    f[p] = data[ get_index( index_i + (( p-1 )*delta_ex_a ) ) ]; // without interpolation this becomes simply data[ get_index( i+p, j+q )]

  const double a_i = a_start + ( index_i - delta_ex_a ) * delta_a; // without interpolation this is simply a_i = a_start + i * delta_a
  const double delta_a_local = delta_a * delta_ex_a; // without interpolation this is simply delta_a

  // f_a[q] holds the values of f[p][q] linearly interpolated in a-direction.
  // Similar to the notation above, f_a[0] corresponds to the value at (a,b[j])
  // and f_a[1] corresponds to the value at (a,b[j+1]).
  // Note that these points are at a not at a[i+p]!
  double f_a = f[0] + ( a - a_i ) / delta_a_local * ( f[1] - f[0] );

  return f_a;
}


/**
 * Utility routine to get the index number (line number - 1) of the point P=(a[i]) 
 * such that data[index] is the tabulated value corresponding to the point P.
 *
 * @param[in] i point first dimension
 * @return index number (line number - 1)
 */
int interpolation1d::get_index( const int i) const
{
  return i;
  //add 1 to get the line number of the desired value in the data file (data vector has offset 0, line numbers have offset 1)
}



/**
 * Configure 2-dimensional interpolation
 *
 * @param[in] n_i_arg total numbers of tabulated points on the a "axes"
 * @param[in] n_j_arg total numbers of tabulated points on the b "axes"
 * @param[in] delta_a_arg spacings of tabulated values in a-direction
 * @param[in] delta_b_arg spacings of tabulated values in b-direction
 * @param[in] a_start_arg the smallest tabulated value of a
 * @param[in] b_start_arg the smallest tabulated value of b
 * @return data, the value of data at (a) found via interpolation from the tabulated values of f
 */
void interpolation2d::configure( const int n_i_arg, const int n_j_arg, const double delta_a_arg, const double delta_b_arg, const double a_start_arg, const double b_start_arg, string filename_arg )
{
  // total numbers of tabulated points on the a, b "axes"
  n_i = n_i_arg;
  n_j = n_j_arg;
  
  // spacings of tabulated values in a, b-direction
  delta_a = delta_a_arg;
  delta_b = delta_b_arg;

  // the smallest tabulated values of a, b, c and d
  a_start = a_start_arg;
  b_start = b_start_arg;
  
  filename = filename_arg;
  
  loadData();
}


/**
 * Reserve memory and load the data.
 */
void interpolation2d::loadData()
{
  expected_entries = n_i * n_j;
  //allocate memory based on the known size (# of entries) of the table
  //a matter of efficieny - not necessary
  data.reserve( expected_entries + 1 );

  if ( readTables() == DATA_SUCCESS )
  {
    cout << "data readout successful: " << filename << endl;
  }
}


/**
 * 2-dimensional interpolation / extrapolation of the tabulated values at a given point (a,b)
 * This routine serves as a wrapper to different methods of interpolation (linear, polynomial, cubic splines). 
 *
 * @param[in] a point first dimension
 * @param[in] b point second dimension
 * @return data = exp(f), the value of data at (a,b,c,d) found via interpolation from the tabulated values of f
 */
double interpolation2d::getInterpolatedData( const double a, const double b) const
{
  return linearInterpolation( a, b);
}



/**
 * 2-dimensional interpolation / extrapolation of the tabulated values at a given point (a,b)
 * This routine is closely tailored to the known structure of the tabulated data (stored in interpolation2d::data).
 * See below for a short documentation of the required data structure.
 *
 * @param[in] a point first dimension
 * @param[in] b point second dimension
 * @return data, the value of data at (a,b) found via interpolation from the tabulated values
 */
double interpolation2d::linearInterpolation( const double a, const double b) const
{
  // find indices i,j such that
  // a[i] <= a < a[i+1], b[j] <= b < b[j+1]
  const int i = int(( a - a_start ) / delta_a ); 
  const int j = int(( b - b_start ) / delta_b );

  if ( i < 0 || j < 0)
  {
    cout << "unexpected - (i,j) not in range" << endl;
    cout << n_i << "  " << n_j << "  "<< endl;
    cout << i << "  " << j   << endl;
    cout << a << "  " << b  << endl;
    cout << a_start << "  " << b_start  << endl;

    return 0;
  }

  // temporary vector objects used to initialize the multi-dimensional arrays based on STL vectors
  vector<double> t1( 2, 0 );
  vector< vector<double> > f( 2, t1 );

  int delta_ex_a,delta_ex_b;
  int index_i, index_j;
  
  if ( i < (n_i - 1) ) // interpolated linearly in a direction
  {
    index_i = i+1;
    delta_ex_a = 1;
  }
  else
  {
    index_i = (n_i - 1);
    delta_ex_a = 2;
  }
  
  if ( j < (n_j - 1) ) // interpolated linearly in b direction
  {
    index_j = j+1;
    delta_ex_b = 1;
  }
  else
  {
    index_j = (n_j - 1);
    delta_ex_b = 2;
  }

  for ( int p = 0; p <= 1; p++ )
    for ( int q = 0; q <= 1; q++ )
      f[p][q] = data[ get_index( index_i + (( p-1 )*delta_ex_a ), index_j + (( q-1 )*delta_ex_b ) )]; // without interpolation this becomes simply data[ get_index( i+p, j+q )]

  const double a_i = a_start + ( index_i - delta_ex_a ) * delta_a; // without interpolation this is simply a_i = a_start + i * delta_a
  const double delta_a_local = delta_a * delta_ex_a; // without interpolation this is simply delta_a
  
  const double b_j = b_start + ( index_j - delta_ex_b ) * delta_b; // without interpolation this is simply b_i = b_start + i * delta_b
  const double delta_b_local = delta_b * delta_ex_b; // without interpolation this is simply delta_b

  // f_a[q] holds the values of f[p][q] linearly interpolated in a-direction.
  // Similar to the notation above, f_a[0] corresponds to the value at (a,b[j])
  // and f_a[1] corresponds to the value at (a,b[j+1]).
  // Note that these points are at a not at a[i+p]!
  vector<double> f_a( 2, 0 );
  for ( int q = 0; q <= 1; q++ )
    f_a[q] = f[0][q] + ( a - a_i ) / delta_a_local * ( f[1][q] - f[0][q] );

  // f_ab holds the values of f_a[q] linearly interpolated in b-direction.
  // thus the final result!
  double f_ab = f_a[0] + ( b - b_j ) / delta_b_local * ( f_a[1] - f_a[0] );

  return f_ab;
}


/**
 * Utility routine to get the index number (line number - 1) of the point P=(a[i],b[j]) 
 * such that data[index] is the tabulated value corresponding to the point P.
 *
 * @param[in] i point first dimension
 * @param[in] j point second dimension
 * @return index number (line number - 1)
 */
int interpolation2d::get_index( const int i, const int j) const
{
  return (( i * n_j ) + ( j ) );
  //add 1 to get the line number of the desired value in the data file (data vector has offset 0, line numbers have offset 1)
}





/**
 * Configure 3-dimensional interpolation
 *
 * @param[in] n_i_arg total numbers of tabulated points on the a "axes"
 * @param[in] n_j_arg total numbers of tabulated points on the b "axes"
 * @param[in] n_k_arg total numbers of tabulated points on the c "axes"
 * @param[in] delta_a_arg spacings of tabulated values in a-direction
 * @param[in] delta_b_arg spacings of tabulated values in b-direction
 * @param[in] delta_c_arg spacings of tabulated values in c-direction
 * @param[in] a_start_arg the smallest tabulated value of a
 * @param[in] b_start_arg the smallest tabulated value of b
 * @param[in] c_start_arg the smallest tabulated value of c
 * @return data, the value of data at (a) found via interpolation from the tabulated values of f
 */
void interpolation3d::configure( const int n_i_arg, const int n_j_arg, const int n_k_arg, const double delta_a_arg, const double delta_b_arg, const double delta_c_arg, const double a_start_arg, const double b_start_arg, const double c_start_arg, string filename_arg )
{
  // total numbers of tabulated points on the a, b, c "axes"
  n_i = n_i_arg;
  n_j = n_j_arg;
  n_k = n_k_arg;
  
  // spacings of tabulated values in a, b, c-direction
  delta_a = delta_a_arg;
  delta_b = delta_b_arg;
  delta_c = delta_c_arg;

  // the smallest tabulated values of a, b, c 
  a_start = a_start_arg;
  b_start = b_start_arg;
  c_start = c_start_arg;
  
  filename = filename_arg;
  
  loadData();
}


/**
 * Reserve memory and load the data.
 */
void interpolation3d::loadData()
{
  expected_entries = n_i * n_j * n_k;
  //allocate memory based on the known size (# of entries) of the table
  //a matter of efficieny - not necessary
  data.reserve( expected_entries + 1 );

  if ( readTables() == DATA_SUCCESS )
  {
    cout << "data readout successful: " << filename << endl;
  }
}


/**
 * 3-dimensional interpolation / extrapolation of the tabulated values at a given point (a,b,c)
 * This routine serves as a wrapper to different methods of interpolation (linear, polynomial, cubic splines). 
 *
 * @param[in] a point first dimension
 * @param[in] b point second dimension
 * @param[in] c point second dimension
 * @return data = f, the value of data at (a,b,c) found via interpolation from the tabulated values of f
 */
double interpolation3d::getInterpolatedData( const double a, const double b, const double c) const
{
  return linearInterpolation( a, b, c);
}



/**
 * 3-dimensional interpolation / extrapolation of the tabulated values at a given point (a,b,c)
 * This routine is closely tailored to the known structure of the tabulated data (stored in interpolation3d::data).
 * See below for a short documentation of the required data structure.
 *
 * @param[in] a point first dimension
 * @param[in] b point second dimension
 * @param[in] c point second dimension
 * @return data, the value of data at (a,b,c) found via interpolation from the tabulated values
 */
double interpolation3d::linearInterpolation( const double a, const double b, const double c) const
{
  // find indices i,j,k such that
  // a[i] <= a < a[i+1], b[j] <= b < b[j+1], c[k] <= c < c[k+1]
  const int i = int(( a - a_start ) / delta_a ); 
  const int j = int(( b - b_start ) / delta_b );
  const int k = int(( c - c_start ) / delta_c );

  if ( i < 0 || j < 0 || k < 0)
  {
    cout << "unexpected - (i,j,k) not in range" << endl;
    cout << n_i << "  " << n_j << "  " << n_k << "  " << endl;
    cout << i << "  " << j << "  " << k << endl;
    cout << a << "  " << b << "  " << c << endl;
    cout << a_start << "  " << b_start << "  " << c_start << endl;

    return 0;
  }

  // temporary vector objects used to initialize the multi-dimensional arrays based on STL vectors
  vector<double> t1( 2, 0 );
  vector< vector<double> > t2( 2, t1 );
  vector< vector< vector<double> > > f( 2, t2 );

  int delta_ex_a,delta_ex_b,delta_ex_c;
  int index_i, index_j, index_k;
  
  if ( i < (n_i - 1) ) // interpolated linearly in a direction
  {
    index_i = i+1;
    delta_ex_a = 1;
  }
  else
  {
    index_i = (n_i - 1);
    delta_ex_a = 2;
  }
  
  if ( j < (n_j - 1) ) // interpolated linearly in b direction
  {
    index_j = j+1;
    delta_ex_b = 1;
  }
  else
  {
    index_j = (n_j - 1);
    delta_ex_b = 2;
  }
    
  if ( k < (n_k - 1) ) // interpolated linearly in c direction
  {
    index_k = k+1;
    delta_ex_c = 1;
  }
  else
  {
    index_k = (n_k - 1);
    delta_ex_c = 2;
  }

  for ( int p = 0; p <= 1; p++ )
    for ( int q = 0; q <= 1; q++ )
      for ( int r = 0; r <= 1; r++ )
        f[p][q][r] = data[ get_index( index_i + (( p-1 )*delta_ex_a ), index_j + (( q-1 )*delta_ex_b ), index_k + (( r-1 )*delta_ex_c ) )]; // without interpolation this becomes simply data[ get_index( i+p, j+q, k+r )]

  const double a_i = a_start + ( index_i - delta_ex_a ) * delta_a; // without interpolation this is simply a_i = a_start + i * delta_a
  const double delta_a_local = delta_a * delta_ex_a; // without interpolation this is simply delta_a
  
  const double b_j = b_start + ( index_j - delta_ex_b ) * delta_b; // without interpolation this is simply b_i = b_start + i * delta_b
  const double delta_b_local = delta_b * delta_ex_b; // without interpolation this is simply delta_b
  
  const double c_k = c_start + ( index_k - delta_ex_c ) * delta_c; // without interpolation this is simply c_k = c_start + k * delta_c
  const double delta_c_local = delta_c * delta_ex_c; // without interpolation this is simply delta_c
      
  // f_a[q][r] holds the values of f[p][q][r] linearly interpolated in a-direction.
  // Similar to the notation above, f_a[0][0] corresponds to the value at (a,b[j],c[k])
  // and f_a[1],[1] corresponds to the value at (a,b[j+1],c[k+1]).
  // Note that these points are at a not at a[i+p]!
  vector< vector<double> > f_a( 2, t1 );
  for ( int q = 0; q <= 1; q++ )
    for ( int r = 0; r <= 1; r++ )
      f_a[q][r] = f[0][q][r] + ( a - a_i ) / delta_a_local * ( f[1][q][r] - f[0][q][r] );

  // f_ab[q] holds the values of f_a[q][r] linearly interpolated in b-direction.
  vector<double> f_ab( 2, 0 );
  for ( int r = 0; r <= 1; r++ )
    f_ab[r] = f_a[0][r] + ( b - b_j ) / delta_b_local * ( f_a[1][r] - f_a[0][r] );

  // f_abc holds the values of f_ab[r] linearly interpolated in c-direction.
  // thus the final result!
  double f_abc = f_ab[0] + ( c - c_k ) / delta_c_local * ( f_ab[1] - f_ab[0] );

  return f_abc;
}


/**
 * Utility routine to get the index number (line number - 1) of the point P=(a[i],b[j],c[k]) 
 * such that data[index] is the tabulated value corresponding to the point P.
 *
 * @param[in] i point first dimension
 * @param[in] j point second dimension
 * @param[in] k point second dimension
 * @return index number (line number - 1)
 */
int interpolation3d::get_index( const int i, const int j, const int k) const
{
  return (( i * n_j * n_k ) + ( j * n_k ) + k );
  //add 1 to get the line number of the desired value in the data file (data vector has offset 0, line numbers have offset 1)
}

