//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/coordinateBins.cpp $
//$LastChangedDate: 2018-04-19 08:34:36 +0200 (Do, 19. Apr 2018) $
//$LastChangedRevision: 2723 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <stdexcept>
#include <vector>
#include <string>
#include <math.h>
#include <iostream>

#include "coordinateBins.h"
#include "configuration.h"

using std::vector;
using std::cout;
using std::endl;

using namespace ns_casc;


coordinateBins::coordinateBins( const int _size, const double _min, const double _max ) :
  _min_real( _min ), _max_real( _max ), 
  _min_index_limit( 0 ), _max_index_limit( _size - 1 ),
  _min_index_active( 0 ), _max_index_active( _size - 1 ), 
  _negative_indices( false )
{
  bins.resize( _size );

  if ( _size > 0 )
  {
    _delta_x = ( _max_real - _min_real ) / _size;
  }
  else
  {
    _delta_x = 0;
  }

  double tmp = _min_real;
  for ( unsigned int i = _min_index_limit; i <= _max_index_limit; i++ )
  {
    bins[i].left = tmp;
    tmp += _delta_x;
    bins[i].right = tmp;
  }
}



void coordinateBins::reshape( const int _size, const double _min, const double _max )
{
  _negative_indices = false;
  _min_real = _min;
  _max_real = _max;
  _min_index_limit = 0;
  _max_index_limit = _size - 1;
  _min_index_active = 0;
  _max_index_active = _size - 1;

  bins.clear();
  bins.resize( _size );

  if ( _size > 0 )
  {
    _delta_x = ( _max_real - _min_real ) / _size;
  }
  else
  {
    _delta_x = 0;
  }

  double tmp = _min_real;
  for ( unsigned int i = _min_index_limit; i <= _max_index_limit; i++ )
  {
    bins[i].left = tmp;
    tmp += _delta_x;
    bins[i].right = tmp;
    bins[i].content = 0;
  }
}


int coordinateBins::getIndex( const double _x ) const
{
  int index = static_cast<int>(( _x - _min_real ) / _delta_x );

  index -= _min_index_limit;

  if ( index < _min_index_limit )
  {
    index = _min_index_limit;
  }
  if ( index > _max_index_limit )
  {
    index = _max_index_limit;
  }

  return index;
}



const coordinateSubBin& coordinateBins::operator[]( const int _index ) const
{
  int tempIndex = _index - _min_index_limit;
  return bins[ tempIndex ];
}


//todo
void coordinateEtaBins::populateEtaBinsManually( std::vector<Particle> & parts, coordinateBins& _dNdEta, const double _etaShift, double _timenow, double& _dt, const double _dx, const double Netabins, const double eta_min_fixed, const double eta_max_fixed )
{
//   const double n_eta = 2 * static_cast<int>( eta_max_fixed / _dEta_fine );
//   _dNdEta.reshape( n_eta, eta_min_fixed, eta_max_fixed );
//   for ( unsigned int i = 0; i < parts.size(); i++ )
//   {    
//     parts[i].eta = parts[i].Pos.Rapidity();
// 
//     if ( parts[i].Pos.T() < ( _timenow + _dt ) )
//     {
//       _dNdEta.increase( parts[i].eta );
//     }
//   }  
//   int centralIndex = this->getCentralIndex();
  
  double deltaEta = (eta_max_fixed-eta_min_fixed)/Netabins;
  
  for ( unsigned int i = 0; i < (bins.size()); ++i )
  {
    double leftCoordinate = eta_min_fixed + i*deltaEta;
    double rightCoordinate = eta_min_fixed + (i+1)*deltaEta;
    
    
    bins[i].setLeftRight( leftCoordinate, rightCoordinate);
  }  
  _min_index_active = 0;
  _max_index_active = bins.size()-1;
  _min_index_limit = 0;
  _max_index_limit =  bins.size()-1;
  
}



void coordinateEtaBins::populateEtaBins( std::vector<Particle> & parts, coordinateBins& _dNdEta, const double _etaShift, double _timenow, double& _dt, const double _dx, const double _dEta_fine )
{
  const double delta_eta_fine = _dEta_fine;
  //WARNING was hardcoded to 8. i put 2.
  const double eta_max = 2.0;
  const double n_eta = 2 * static_cast<int>( eta_max / delta_eta_fine );

  _dNdEta.reshape( n_eta, -eta_max, + eta_max );

  for ( unsigned int i = 0; i < parts.size(); i++ )
  {    
    parts[i].eta = parts[i].Pos.Rapidity();

    if ( parts[i].Pos.T() < ( _timenow + _dt ) )
    {
      _dNdEta.increase( parts[i].eta );
    }
  }
  
  int centralIndex = this->getCentralIndex();
  
  for ( unsigned int i = 0; i < bins.size(); ++i )
  {
    bins[i].setLeftRight( -infinity, infinity );
  }


  double dz;

  _dt = _dx;

  int np, npr, npl;
  unsigned int nsum1 = 0, nsum2 = 0;
  unsigned int nMax = NinEtaBin;
  const int dNdEta_centralIndex = static_cast<int>( _dNdEta.size() ) / 2;
  const double dEta = _dNdEta.get_dx();

  //---------- populate central eta bin ----------
  np = dNdEta_centralIndex + static_cast<int>( _etaShift / dEta );
  npr = np;
  npl = np - 1;

  nsum1 = _dNdEta[npr].content;
  nsum2 = _dNdEta[npl].content;
  while ( nsum1 + nsum2 < nMax && ( npr < _dNdEta.max_index() || npl > _dNdEta.min_index() ) )
  {
    if ( npr < _dNdEta.max_index() )
    {
      ++npr;
      nsum1 += _dNdEta[npr].content;
    }
    if ( npl > _dNdEta.min_index() )
    {
      --npl;
      nsum2 += _dNdEta[npl].content;
    }
  }
  
  bins[centralIndex].right = _dNdEta[npr].right;
  bins[centralIndex].left = _dNdEta[npl].left;
  bins[centralIndex].content = nsum1 + nsum2;

  if ( bins.size() >= 3 )
  {
    bins[centralIndex + 1].left = bins[centralIndex].right;
    bins[centralIndex - 1].right = bins[centralIndex].left;
  }
  
  dz = _timenow * ( tanh( bins[centralIndex].right ) - tanh( bins[centralIndex].left ) );
  if ( dz < _dt )
  {
    _dt = dz;
  }
  
   //---------- populate central eta bin ----------
  

  unsigned int nRest;
  unsigned int index = centralIndex;
  //---------- populate bins with eta > 0 ----------
  nRest = 0;
  for ( int i = npr + 1; i < _dNdEta.size(); i++ )
  {
    nRest += _dNdEta[i].content;
  }

  if ( nRest < nMax )
  {
    nMax = NinEtaBin / 2;
  }
  while ( nRest >= nMax && index < _max_index_limit - 1 )
  {
    ++index;
    nsum1 = 0;
    while ( nsum1 < nMax && npr < _dNdEta.max_index() )
    {
      ++npr;
      nsum1 += _dNdEta[npr].content;
    }
    
    bins[index].right = _dNdEta[npr].right;
    bins[index].content = nsum1;
    if ( index + 1 < bins.size() )
    {
      bins[index + 1].left = bins[index].right;
    }
    
    dz = _timenow * ( tanh( bins[index].right ) - tanh( bins[index].left ) );
    if ( dz < _dt )
    {
      _dt = dz;
    }
    
    nRest -= nsum1;
    if ( nRest < nMax )
    {
      nMax = NinEtaBin / 2;
    }
  }
  _max_index_active = index;
  //---------- populate bins with eta > 0 ----------
  
  

  //---------- populate bins with eta < 0 ----------
  nRest = 0;
  index = centralIndex;
  
  for ( int i = 0; i < npl; i++ )
  {
    nRest += _dNdEta[i].content;
  }

  nMax = NinEtaBin;
  if ( nRest < nMax )
  {
    nMax = NinEtaBin / 2;
  }

  while ( nRest >= nMax && npl >= _dNdEta.min_index() && index > _min_index_limit + 1 )
  {
    index--;
    nsum2 = 0;
    while ( nsum2 < nMax )
    {
      --npl;
      nsum2 += _dNdEta[npl].content;
    }
    bins[index].left = _dNdEta[npl].left;
    bins[index].content = nsum2;
    if ( index - 1 >= 0 )
    {
      bins[index - 1].right = bins[index].left;
    }
    
    dz = _timenow * ( tanh( bins[index].right ) - tanh( bins[index].left ) );
    if ( dz < _dt )
    {
      _dt = dz;
    }

    nRest -= nsum2;
    if ( nRest < nMax )
    {
      nMax = NinEtaBin / 2;
    }
  }
  _min_index_active = index;
  //---------- populate bins with eta < 0 ----------
  if ( _min_index_active < 0 || _max_index_active >= bins.size() )
  {
    std::string errMsg = "Index out of range in populateEtaBins. EstimatedMaxNumber in constructEtaBins might need adjustment.";
    throw eCoordBins_error( errMsg );
  }
  
  _dt = timestepScaling * _dt;
  
}
/**
 * @param[in] particlesSize Total number of particles
 * @param[in] _nTest testparticles, currentlx not used here
 * @param[in] _RA not used here
 * @param[in] _dy y-length of cells 
 * @param[in] _dx x-length of cells
 * @param[in] _transLen transverse length of the collision. Will be fixed during the simulation.
 * @param[in] _NperCell Given particles per cell value. Rapidity bins will be adjusted.
 */
int coordinateEtaBins::constructEtaBinsManually( const int _NperCell, const double _transLen, const double _dx, const double _dy, const double _RA, const int _nTest, const int particlesSize )
{
    int _IZ = 101;

    bins.resize( _IZ );
    _min_index_limit = 0;
    _max_index_limit = bins.size() - 1;

    cout << "Manual adjustment of eta-bins:  IZ = " << _IZ << endl;
    
//     if ( particlesSize <= ( 8 * NinEtaBin ) )
//     {
//         cout << "N = " << particlesSize << "  NinEtaBin = " << NinEtaBin << "  IZ = " << _IZ << endl;
//         cout << "recommended number of test particles > " << ( ( 8 * NinEtaBin * _nTest ) / particlesSize ) << endl;
//         std::string errMsg = "insufficient number of test particles";
//         throw eCoordBins_error( errMsg );
//     }

    return _IZ;
}


/**
 * @param[in] particlesSize Total number of particles
 * @param[in] _nTest testparticles, currentlx not used here
 * @param[in] _RA not used here
 * @param[in] _dy y-length of cells 
 * @param[in] _dx x-length of cells
 * @param[in] _transLen transverse length of the collision. Will be fixed during the simulation.
 * @param[in] _NperCell Given particles per cell value. Rapidity bins will be adjusted.
 */
int coordinateEtaBins::constructEtaBinsPPb( const int _NperCell, const double _transLen, const double _dx, const double _dy, const double _RA, const int _nTest, const int particlesSize )
{
    // The following line relies on Woods-Saxon density distribution:
    double initialArea = _transLen*_transLen;
    NinEtaBin = int( initialArea / ( _dx * _dy ) * _NperCell ); //Total Number per eta bin

    int estimatedMaxNumber = particlesSize * 1.01; //WARNING //HACK
    
    int _IZ = ( estimatedMaxNumber / 2 / NinEtaBin  + 2 ) * 2 + 1;

    bins.resize( _IZ );
    _min_index_limit = 0;
    _max_index_limit = bins.size() - 1;

    cout << "N = " << particlesSize << "  NinEtaBin = " << NinEtaBin << "  IZ = " << _IZ << endl;
    
//     if ( particlesSize <= ( 8 * NinEtaBin ) )
//     {
//         cout << "N = " << particlesSize << "  NinEtaBin = " << NinEtaBin << "  IZ = " << _IZ << endl;
//         cout << "recommended number of test particles > " << ( ( 8 * NinEtaBin * _nTest ) / particlesSize ) << endl;
//         std::string errMsg = "insufficient number of test particles";
//         throw eCoordBins_error( errMsg );
//     }

    return _IZ;
}

/**
 * The new parameter RA corresponds to the former variable WoodSaxon&
 * _param.RA, i.e. to the radius parameter
 *
 * The new parameter particlesSize corresponds to the former particles.size()
 */
int coordinateEtaBins::constructEtaBins( const int _NperCell, const double _b, const double _dx, const double _dy, const double _RA, const int _nTest, const int particlesSize )
{
  // The following line relies on Woods-Saxon density distribution:
  double initialArea = 4 * ( _RA - _b / 2 ) * sqrt(( _RA - _b / 2 ) * ( _RA + _b / 2 ) ); 
  NinEtaBin = int( initialArea / ( _dx * _dy ) * _NperCell ); //number of transverse cells times 10 - required value

  int estimatedMaxNumber = particlesSize * 1.3;
  if ( Particle::N_light_flavor == 0 )
  {
    estimatedMaxNumber *= 1.4;  // ugly fix, empirical value
  }
  int _IZ = ( estimatedMaxNumber / 2 / NinEtaBin  + 2 ) * 2 + 1;

  bins.resize( _IZ );
  _min_index_limit = 0;
  _max_index_limit = bins.size() - 1;

  cout << "number of cells in one etabin=" << int( pow(( 2 * _RA / _dx ), 2 ) )
       << "\t" << "particle number in one Etabin=" << NinEtaBin << endl;
  if ( particlesSize <= ( 8 * NinEtaBin ) )
  {
    cout << "N = " << particlesSize << "  should be = " << 8*NinEtaBin << "  IZ = " << _IZ << endl;
    cout << "recommended number of test particles > " << ( ( 8 * NinEtaBin * _nTest ) / particlesSize ) << endl;
    std::string errMsg = "insufficient number of test particles";
    throw eCoordBins_error( errMsg );
  }

  return _IZ;
}



int coordinateEtaBins::getCentralIndex() const
{
  if ( bins.size() % 2 != 1 )
  {
    string errMsg = "number of eta bins must be odd";
    throw eCoordBins_error( errMsg );
  }
  else
  {
    return ( static_cast<int>( bins.size() ) / 2 );
  }
}



int coordinateEtaBins::getIndex( const double _eta ) const
{
  unsigned int index = this->getCentralIndex();
  if ( _eta >= 0 )
  {
    while ( _eta >= bins[index].right && index < _max_index_limit )
    {
      ++index;
    }
  }
  else
  {
    while ( _eta < bins[index].left && index >= _min_index_limit )
    {
      --index;
    }
  }
  
  return index;
}



// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on;  replace-tabs on;
