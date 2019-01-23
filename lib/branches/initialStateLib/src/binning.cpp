//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/binning.cpp $
//$LastChangedDate: 2018-04-27 19:16:12 +0200 (Fr, 27. Apr 2018) $
//$LastChangedRevision: 2741 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


#include "binning.h"

#include <iostream>
#include <math.h>

using namespace std;

/**
 * set all variables to its defaults
 */
binning::binning() :
  x_min(1e20),
  x_max(-1e20),
  binWidth(),
  nBins(100),
  nEntries(),
  theList(),
  theBins(nBins,0),
  finished(false),
  bintype( BINTYPE::dummy ),
  filename("defaultBinning.dat")
{
}

/**
 * set all variables to its defaults
 *
 * @param[in] name file name for output
 */
binning::binning(const string & name) :
  x_min(1e20),
  x_max(-1e20),
  binWidth(),
  nBins(100),
  nEntries(),
  theList(),
  theBins(nBins,0),
  finished(false),
  bintype( BINTYPE::dummy ),
  filename( name )
{
}

/**
 * @param[in] n number of bins
 *
 * first do the storage into a list
 */
binning::binning(const int n) :
  x_min(1e20),
  x_max(-1e20),
  binWidth(),
  nBins(n),
  nEntries(),
  theList(),
  theBins(nBins,0),
  finished(false),
  bintype( BINTYPE::fixed_number ),
  filename("defaultBinning.dat")
{
}

/**
 * @param[in] name file name for output
 * @param[in] n number of bins
 *
 * first do the storage into a list
 */
binning::binning(const string & name, const int n) :
  x_min(1e20),
  x_max(-1e20),
  binWidth(),
  nBins(n),
  nEntries(),
  theList(),
  theBins(nBins,0),
  finished(false),
  bintype( BINTYPE::fixed_number ),
  filename(name)
{
}

/**
 * @param[in] width The bin width
 *
 * first do the storage into a list
 */
binning::binning(const double width) :
  x_min(1e20),
  x_max(-1e20),
  binWidth(width),
  nBins(),
  nEntries(),
  theList(),
  theBins(nBins,0),
  finished(false),
  bintype( BINTYPE::fixed_width ),
  filename("defaultBinning.dat")
{
}

/**
 * @param[in] name file name for output
 * @param[in] width The bin width
 *
 * first do the storage into a list
 */
binning::binning(const string & name, const double width) :
  x_min(1e20),
  x_max(-1e20),
  binWidth(width),
  nBins(),
  nEntries(),
  theList(),
  theBins(nBins,0),
  finished(false),
  bintype( BINTYPE::fixed_width ),
  filename(name)
{
}

/**
 * @param[in] min_arg lower bound
 * @param[in] max_arg upper bound
 * @param[in] n number of bins
 *
 * All parameters for the binning are given and entries are added
 * directly to the correspondin bin (no intermediate storage in a list)
 */
binning::binning(const double min_arg, const double max_arg, const int n) :
  x_min( min_arg ),
  x_max( max_arg ),
  binWidth( (x_max - x_min) / n ),
  nBins( n ),
  nEntries(),
  theList(),
  theBins(nBins,0),
  finished(true),
  bintype( BINTYPE::fixed_ALL ),
  filename("defaultBinning.dat")
{
}

/**
 * @param[in] name file name for output
 * @param[in] min_arg lower bound
 * @param[in] max_arg upper bound
 * @param[in] n number of bins
 *
 * All parameters for the binning are given and entries are added
 * directly to the correspondin bin (no intermediate storage in a list)
 */
binning::binning(const string & name, const double min_arg, const double max_arg, const int n) :
  x_min( min_arg ),
  x_max( max_arg ),
  binWidth( (x_max - x_min) / n ),
  nBins( n ),
  nEntries(),
  theList(),
  theBins(nBins,0),
  finished(true),
  bintype( BINTYPE::fixed_ALL ),
  filename( name )
{
}


binning::~binning()
{
}

/**
 * @param[in] min_arg lower bound
 * @param[in] max_arg upper bound
 * @param[in] width_arg bin width
 *
 * The number of bins is calculated.
 */
void binning::setMinMaxWidth(const double min_arg, const double max_arg, const double width_arg)
{
  x_min = min_arg;
  x_max = max_arg;
  binWidth = width_arg;
  nBins = int((x_max-x_min)/binWidth);
  if( (x_min + nBins*binWidth) < x_max )
    nBins++;
  theBins.resize(nBins,0);
  bintype = BINTYPE::fixed_ALL;
  nEntries = 0;
  finished = true;

}

/**
 * @param[in] min_arg lower bound
 * @param[in] max_arg upper bound
 * @param[in] n number of bins
 */
void binning::setMinMaxN(const double min_arg, const double max_arg, const int n)
{
  x_min = min_arg;
  x_max = max_arg;
  nBins = n;
  theBins.resize(nBins,0);
  binWidth = (x_max - x_min) / nBins;
  bintype = BINTYPE::fixed_ALL;
  nEntries = 0;
  finished = true;
}


/**
 * @param[in] n number of bins
 *
 * **NOTE:** This routine does not adjust the bins or allocates
 * memory. *Use with care!*
 */
void binning::setNumber(const int n)
{
  nBins = n;
  bintype = BINTYPE::fixed_number;
}

/**
 * @param[in] w bin width
 *
 * **NOTE:** This routine does not adjust the bins or allocates
 * memory. *Use with care!*
 */
void binning::setWidth(const double w)
{
  binWidth = w;
  bintype = BINTYPE::fixed_width;
}


/**
 * @param[in] x value of element to add
 *
 * If bins are known, element can be inserted directly into the
 * desired bin. Otherwise, the value is added to a list and the
 * binning is done at the very end by #finish.
 */
void binning::add(const double x)
{
  //if min and max are known AND data collection hasn't started yet, x can be sorted into the appropriate bin immediately
  if (bintype == BINTYPE::fixed_ALL && theList.size() == 0)
  {
    int index = int((x - x_min)/binWidth);
    if (index >= nBins)
      index = nBins-1;
    else if (index < 0)
    {
      index = 0;
      //       cout << "smaller than bin: " << x << endl;
    }
    ++theBins[index];
  }
  //otherwise all data has to be collected before sorting into bins starts
  else
    theList.push_back(x);   //push x into the temporary list

  ++nEntries;
}

/**
 * if necessary, sort into bins by calling #finish
 */
void binning::print()
{
  fstream file(filename.c_str(),ios::out);
  print(file);
  file.close();
}

/**
 * if necessary, sort into bins by calling #finish
 */
void binning::print(std::fstream & file)
{
  doFinish();

  file << "#bin width: " << binWidth
       << "   number of binned events: " << nEntries << endl;

  for (int i=0; i<nBins; i++)
  {
    file << x_min + (i+0.5)*binWidth << "\t"
         << theBins[i] << "\t"
         << theBins[i]/binWidth << "\t"
         << theBins[i]/(binWidth*nEntries) << endl;
  }
}


void binning::finish()
{
  if (bintype != BINTYPE::fixed_ALL)
  {
    if (theList.size() == 0)
    {
      std::string errMsg = "Error in 'finish' of "+filename+"! size=0";
      throw eBINNING_error( errMsg );
    }

    x_min = x_max = theList[0];
    for(unsigned int i=1; i<theList.size(); i++)
    {
      if(theList[i] < x_min)
        x_min = theList[i];
      if(theList[i] > x_max)
        x_max = theList[i];
    }
  }

  if (bintype == BINTYPE::fixed_width)
  {
    nBins = int((x_max - x_min) / binWidth) + 1;
    theBins.resize(nBins,0);
  }
  else
  {
    binWidth = (x_max - x_min) / nBins;
  }

  for(unsigned int i=0; i<theList.size(); i++)
  {
    int index = int((theList[i] - x_min)/binWidth);
    if (index >= nBins)
      index = nBins-1;
    ++theBins[index];
  }

  finished = true;

}


void binning::getVectors(vector<double>& x, vector<double>& y, double& width)
{
  doFinish();

  x.clear(); x.reserve(nBins);
  y.clear(); y.reserve(nBins);
  width = binWidth;

  for ( int i = 0; i < nBins; i++ )
  {
    x.push_back( x_min + (i+0.5) * binWidth );
    y.push_back( theBins[i] );
  }
}



double binning::getBinLabel(const int i)
{
  doFinish();

  if (i >= 0 && i < nBins)
    return x_min + (i+0.5)*binWidth;
  else
    return 0;
}


double binning::getBin(const int i)
{
  doFinish();

  if (i >= 0 && i < nBins)
    return theBins[i]/binWidth;
  else
    return 0;
}


int binning::getBinRaw(const int i)
{
  doFinish();

  if (i >= 0 && i < nBins)
    return theBins[i];
  else
    return 0;
}

double binning::getBinRelative(const int i)
{
  doFinish();

  if (i >= 0 && i < nBins)
    return double(theBins[i])/(double(nEntries)*binWidth);
  else
    return 0;
}

//----------------------------------------------------------------------------

/*
  default constructor
*/
binningLog::binningLog() :
  x_min(1e20),
  x_max(-1e20),
  binWidth(),
  nBins(100),            // default number of bins
  nEntries(),
  theList(),
  theBins(nBins,0),
  finished(false),
  bintype( BINTYPE::dummy ),
  filename("defaultBinningLog.dat")
{
}


binningLog::binningLog(const double min_arg, const double max_arg, const int n) :
  x_min( min_arg ),
  x_max( max_arg ),
  binWidth( (log( x_max ) - log( x_min )) / n ),
  nBins( n ),
  nEntries(),
  theList(),
  theBins(nBins,0),
  finished(true),
  bintype( BINTYPE::fixed_ALL ),
  filename("defaultBinningLog.dat")
{
}

binningLog::binningLog(const std::string name, const double min_arg, const double max_arg, const int n) :
  x_min( min_arg ),
  x_max( max_arg ),
  binWidth( (log( x_max ) - log( x_min )) / n ),
  nBins( n ),
  nEntries(),
  theList(),
  theBins(nBins,0),
  finished(true),
  bintype( BINTYPE::fixed_ALL ),
  filename( name )
{
}


binningLog::~binningLog()
{
}


void binningLog::setMinMaxN(const double min_arg, const double max_arg, const int n)
{
  x_min = min_arg;
  x_max = max_arg;
  nBins = n;
  theBins.resize(nBins,0);
  binWidth = (log( x_max ) - log( x_min )) / nBins;
  bintype = BINTYPE::fixed_ALL;
  nEntries = 0;
  finished = true;
}


void binningLog::setNumber(const int n)
{
  nBins = n;
  bintype = BINTYPE::fixed_number;
}


void binningLog::add(const double x)
{
  //if min and max are known AND data collection hasn't started yet, x can be sorted into the appropriate bin immediately
  if (bintype == BINTYPE::fixed_ALL && theList.size() == 0)
  {
    int index = static_cast<int>(( log( x ) - log( x_min ) ) / binWidth );
    if (index >= nBins)
      index = nBins-1;
    else if (index < 0)
    {
      index = 0;
      //       cout << "smaller than bin: " << x << endl;
    }
    ++theBins[index];
  }
  //otherwise all data has to be collected before sorting into bins starts
  else
    theList.push_back(x);   //push x into the temporary list

  ++nEntries;
}


void binningLog::print()
{
  doFinish();

  fstream file(filename.c_str(),ios::out);
  file << "#bin width: " << binWidth
       << "   number of binned events: " << nEntries << endl;

  for (int i=0; i<nBins; i++)
  {
    double label = exp( double( i ) * binWidth + log( x_min ) + binWidth / 2.0 );
    double dx = exp( double( i ) * binWidth + log( x_min ) + binWidth ) - exp( double( i ) * binWidth + log( x_min ) );
    file << label << "\t"
         << theBins[i] << "\t"
         << theBins[i]/dx << "\t"
         << theBins[i]/(dx*nEntries) << endl;
  }
  file.close();
}


void binningLog::finish()
{
  if (bintype != BINTYPE::fixed_ALL)
  {
    x_min = x_max = theList[0];
    for(unsigned int i=1; i<theList.size(); i++)
    {
      if(theList[i] < x_min)
        x_min = theList[i];
      if(theList[i] > x_max)
        x_max = theList[i];
    }
  }

  if (bintype == BINTYPE::fixed_width)
  {
    nBins = int((x_max - x_min) / binWidth) + 1;
  }
  else
  {
    binWidth = (log( x_max ) - log( x_min )) / nBins;
  }

  int index;
  for(unsigned int i=0; i<theList.size(); i++)
  {
    index = static_cast<int>(( log( theList[i] ) - log( x_min ) ) / binWidth );
    if (index >= nBins)
      index = nBins-1;
    ++theBins[index];
  }

  finished = true;
}


double binningLog::getBinLabel(const int i)
{
  doFinish();

  if (i >= 0 && i < nBins)
    return exp( double( i ) * binWidth + log( x_min ) + binWidth / 2.0 );
  else
    return 0;
}


double binningLog::getBin(const int i)
{
  doFinish();

  double dx = exp( double( i ) * binWidth + log( x_min ) + binWidth ) - exp( double( i ) * binWidth + log( x_min ) );

  if (i >= 0 && i < nBins)
    return theBins[i]/dx;
  else
    return 0;
}



int binningLog::getBinRaw(const int i)
{
  doFinish();

  if (i >= 0 && i < nBins)
    return theBins[i];
  else
    return 0;
}




double binningLog::getBinRelative(const int i)
{
  doFinish();

  double dx = exp( double( i ) * binWidth + log( x_min ) + binWidth ) - exp( double( i ) * binWidth + log( x_min ) );

  if (i >= 0 && i < nBins)
    return double(theBins[i])/double(nEntries)/dx;
  else
    return 0;
}

//----------------------------------------------------------------------------

binningValues::binningValues() :
  x_min(1e20),
  x_max(-1e20),
  binWidth(),
  nBins(100),
  nEntries(),
  theBins(nBins,0),
  theBinsSquared(nBins,0),
  theBinsN(nBins,0),
  filename("defaultBinningValues.dat")
{
}


/**
 * @param[in] name filename
 * @param[in] min_arg lower bound
 * @param[in] max_arg upper bound
 * @param[in] width width of bins
 */
binningValues::binningValues(const string & name, const double min_arg, const double max_arg, const double width) :
  x_min(min_arg),
  x_max(max_arg),
  binWidth(width),
  nBins(),
  nEntries(),
  filename(name)
{
  nBins = int ( ( x_max - x_min ) / binWidth );
  if ( ( x_min + nBins * binWidth ) < x_max )
    nBins++;

  theBins.resize(nBins,0.0);
  theBinsSquared.resize(nBins,0.0);
  theBinsN.resize(nBins,0);
}

/**
 * @param[in] name filename
 * @param[in] min_arg lower bound
 * @param[in] max_arg upper bound
 * @param[in] n number of bins
 */
binningValues::binningValues(const string & name, const double min_arg, const double max_arg, const int n) :
  x_min(min_arg),
  x_max(max_arg),
  binWidth(),
  nBins(n),
  nEntries(),
  filename(name)
{
  theBins.resize(nBins,0.0);
  theBinsSquared.resize(nBins,0.0);
  theBinsN.resize(nBins,0);
  binWidth = (x_max - x_min) / nBins;
}

/**
 * @param[in] x x-value indicating the bin
 * @param[in] v value to add
 */
void binningValues::add(const double x, const double v)
{
  int index = int((x - x_min)/binWidth);
  if (index < nBins)
  {
    //index = nBins-1;
    if (index < 0)
    {  
      //do nothing!
      //index = 0;
    }
    if (index >= 0)
    {
      theBins[index] += v;
      theBinsSquared[index] += v * v;
      ++theBinsN[index];

      ++nEntries;
    }
  }
}

void binningValues::print(std::fstream &file)
{
  file << "#bin width: " << binWidth << "   number of binned events: " << nEntries << endl;
  for (int i=0; i<nBins; i++)
  {
    file.width(15);
    file << x_min + (i+0.5)*binWidth;
    file.width(15);
    if(theBinsN[i] != 0)
    {
      file << theBins[i];
      file.width(15);
      file << theBinsN[i];
      file.width(15);
      file << theBins[i] / theBinsN[i] << endl;
    }
    else
    {
      file << 0;
      file.width(15);
      file << 0;
      file.width(15);
      file << 0 << endl;
    }
  }
}

void binningValues::print()
{
  fstream file(filename.c_str(),ios::out);
  print(file);
  file.close();
}

void binningValues::setMinMaxWidth( const double min_arg, const double max_arg, const double width_arg )
{
  x_min = min_arg;
  x_max = max_arg;
  binWidth = width_arg;
  nBins = int ( ( x_max - x_min ) / binWidth );
  if ( ( x_min + nBins * binWidth ) < x_max )
    nBins++;
  theBins.resize ( nBins, 0 );
  theBinsSquared.resize( nBins, 0 );
  nEntries = 0;
  theBinsN.resize ( nBins, 0 );
}

void binningValues::setMinMaxN( const double min_arg, const double max_arg, const int n )
{
  x_min = min_arg;
  x_max = max_arg;
  nBins = n;
  theBins.resize ( nBins, 0 );
  theBinsSquared.resize( nBins, 0 );
  binWidth = ( x_max - x_min ) / nBins;
  nEntries = 0;
  theBinsN.resize ( nBins, 0 );
}

double binningValues::getBinLabel( const int i )
{
  if ( i >= 0 && i < nBins )
    return x_min + ( i + 0.5 ) * binWidth;
  else
    return 0;
}


int binningValues::getBinCount( const int i )
{
  if ( i >= 0 && i < nBins )
    return theBinsN[i];
  else
    return 0;
}


long double binningValues::getBinValue( const int i )
{
  if ( i >= 0 && i < nBins )
    return theBins[i];
  else
    return 0;
}

long double binningValues::getBinValueSquared( const int i )
{
  if ( i >= 0 && i < nBins )
    return theBinsSquared[i];
  else
    return 0;
}

double binningValues::getBinMean( const int i )
{
  if ( i >= 0 && i < nBins )
  {
    if ( theBins[i] > 0 )
    {
      return theBins[i] / theBinsN[i];
    }
  }
  return 0.0;
}


double binningValues::getBinError( const int i )
{
  if ( theBinsN[i] > 0 )
  {
    double sigma = sqrt( theBinsSquared[i] / theBinsN[i] - pow( theBins[i] / theBinsN[i], 2.0 ) );
    return ( sigma / sqrt( theBinsN[i] ) );
  }

  return 0.0;
}

double binningValues::getSigma(const int i )
{
  if ( theBinsN[i] > 0 )
  {
    double sigma = sqrt( theBinsSquared[i] / theBinsN[i] - pow( theBins[i] / theBinsN[i], 2.0 ) );
    return ( sigma );
  }

  return 0.0;
}

/**
 * see source code for details
 */
double binningValues::getAreaNormalization( const int i )
{
  if ( i >= 0 && i < nBins )
  {
    double rmin = i * binWidth + x_min;
    double rmax = ( i + 1 ) * binWidth + x_min;
    return ( M_PI * ( pow ( rmax, 2.0 ) - pow ( rmin, 2.0 ) ) );
  }
  return 0.0;
}

//----------------------------------------------------------------------------

binningLogValues::binningLogValues() :
  x_min(1e20),
  x_max(-1e20),
  binWidth(),
  nBins(100),
  nEntries(),
  theBins(nBins,0),
  theBinsSquared(nBins,0),
  theBinsN(nBins,0),
  filename("defaultBinningValues.dat")
{
}


/**
 * @param[in] name filename
 * @param[in] min_arg lower bound
 * @param[in] max_arg upper bound
 * @param[in] n number of bins
 */
binningLogValues::binningLogValues(const string & name, const double min_arg, const double max_arg, const int n) :
  x_min(min_arg),
  x_max(max_arg),
  binWidth(),
  nBins(n),
  nEntries(),
  filename(name)
{
  theBins.resize(nBins,0.0);
  theBinsSquared.resize(nBins,0.0);
  theBinsN.resize(nBins,0);
  binWidth = (log( x_max ) - log( x_min )) / nBins;
  setFilename(name);
}

/**
 * @param[in] x x-value indicating the bin
 * @param[in] v value to add
 */
void binningLogValues::add(const double x, const double v)
{

  int index = static_cast<int>(( log( x ) - log( x_min ) ) / binWidth );
  if (index >= nBins)
    index = nBins-1;
  else if (index < 0)
    index = 0;

  theBins[index] += v;
  theBinsSquared[index] += v * v;
  ++theBinsN[index];

  ++nEntries;
}

void binningLogValues::print(std::fstream &file)
{
  file << "#bin width: " << binWidth << "   number of binned events: " << nEntries << endl;
  for (int i=0; i<nBins; i++)
  {
    double label = exp( double( i ) * binWidth + log( x_min ) + binWidth / 2.0 ); // +dpt/2.0 to shift value in the middle of the intervall, important for small number of bins
    double dx = exp( double( i ) * binWidth + log( x_min ) + binWidth ) - exp( double( i ) * binWidth + log( x_min ) );
    file.width(15);
    file << label;
    file.width(15);
    if(theBinsN[i] != 0)
    {
      file << theBins[i];
      file.width(15);
      file << theBins[i] / dx;
      file.width(15);
      file << theBinsN[i] << endl;
    }
    else
    {
      file << 0;
      file.width(15);
      file << 0;
      file.width(15);
      file << 0 << endl;
    }
  }
}

void binningLogValues::print()
{
  fstream file(filename.c_str(),ios::out);
  print(file);
  file.close();
}

void binningLogValues::setMinMaxN( const double min_arg, const double max_arg, const int n )
{
  x_min = min_arg;
  x_max = max_arg;
  nBins = n;
  theBins.resize ( nBins, 0 );
  theBinsSquared.resize( nBins, 0 );
  theBinsN.resize ( nBins, 0 );
  binWidth = ( log( x_max ) - log( x_min ) ) / nBins;
  nEntries = 0;
}

double binningLogValues::getBinLabel(const int i)
{
  if (i >= 0 && i < nBins)
    return exp( double( i ) * binWidth + log( x_min ) + binWidth / 2.0 );
  else
    return 0;
}


int binningLogValues::getBinCount( const int i )
{
  if ( i >= 0 && i < nBins )
    return theBinsN[i];
  else
    return 0;
}


long double binningLogValues::getBinValue( const int i )
{
  if ( i >= 0 && i < nBins )
    return theBins[i];
  else
    return 0;
}

long double binningLogValues::getBinValueSquared( const int i )
{
  if ( i >= 0 && i < nBins )
    return theBinsSquared[i];
  else
    return 0;
}

double binningLogValues::getBinMean( const int i )
{
  if ( i >= 0 && i < nBins )
  {
    if ( theBins[i] > 0 )
    {
      return theBins[i] / theBinsN[i];
    }
  }
  return 0.0;
}


double binningLogValues::getBinError( const int i )
{
  if ( theBinsN[i] > 0 )
  {
    double sigma = sqrt( theBinsSquared[i] / theBinsN[i] - pow( theBins[i] / theBinsN[i], 2.0 ) );
    return ( sigma / sqrt( theBinsN[i] ) );
  }

  return 0.0;
}

double binningLogValues::getSigma(const int i )
{
  if ( theBinsN[i] > 0 )
  {
    double sigma = sqrt( theBinsSquared[i] / theBinsN[i] - pow( theBins[i] / theBinsN[i], 2.0 ) );
    return ( sigma );
  }

  return 0.0;
}

//----------------------------------------------------------------------------

binningInRanges::binningInRanges()
{
}

binningInRanges::binningInRanges( const int _dim )
{
  theBinnings.resize( _dim );
  lowValues.resize( _dim );
  highValues.resize( _dim );
}

binningInRanges::~binningInRanges()
{

}

void binningInRanges::setNumberRanges(const int _dim)
{
  theBinnings.resize( _dim );
  lowValues.resize( _dim );
  highValues.resize( _dim );
}

void binningInRanges::setRangevalues(const int _dim, const double minValue, const double maxValue)
{
  if ( _dim < 0 )
  {
    throw eBINNING_error ("Unrecoverable error in binningInRanges...");
  }
  else if ( _dim > static_cast< int >(theBinnings.size()) - 1 )
  {
    throw eBINNING_error ("Unrecoverable error in binningInRanges...");
  }
  else
  {
    lowValues[_dim] = minValue;
    highValues[_dim] = maxValue;
  }
}

void binningInRanges::setMinMaxN(const double min_arg, const double max_arg, const int n)
{
  for ( unsigned int i = 0; i < theBinnings.size(); i++ )
  {
    theBinnings[i].setMinMaxN( min_arg, max_arg, n);
  }
}

void binningInRanges::setMinMaxWidth(const double min_arg, const double max_arg, const double width_arg)
{
  for ( unsigned int i = 0; i < theBinnings.size(); i++ )
  {
    theBinnings[i].setMinMaxWidth( min_arg, max_arg, width_arg );
  }
}

void binningInRanges::add(const double rangeParameter, const double x)
{
  for ( unsigned int i = 0; i < theBinnings.size(); i++ )
  {
    if ( rangeParameter >= lowValues[i] && rangeParameter <= highValues[i] )
    {
      theBinnings[i].add( x );
    }
  }
}

void binningInRanges::getBorders(const int _dim, double& lowBorder, double& highBorder)
{
  lowBorder = lowValues[_dim];
  highBorder = highValues[_dim];
}


binningLogInRanges::binningLogInRanges()
{
}

binningLogInRanges::binningLogInRanges( const int _dim )
{
  theBinnings.resize( _dim );
  lowValues.resize( _dim );
  highValues.resize( _dim );
}

binningLogInRanges::~binningLogInRanges()
{

}

void binningLogInRanges::setNumberRanges(const int _dim)
{
  theBinnings.resize( _dim );
  lowValues.resize( _dim );
  highValues.resize( _dim );
}

void binningLogInRanges::setRangevalues(const int _dim, const double minValue, const double maxValue)
{
  if ( _dim < 0 )
  {
    throw eBINNING_error ("Unrecoverable error in binningInRanges...");
  }
  else if ( _dim > static_cast< int >(theBinnings.size()) - 1 )
  {
    throw eBINNING_error ("Unrecoverable error in binningInRanges...");
  }
  else
  {
    lowValues[_dim] = minValue;
    highValues[_dim] = maxValue;
  }
}

void binningLogInRanges::setMinMaxN(const double min_arg, const double max_arg, const int n)
{
  for ( unsigned int i = 0; i < theBinnings.size(); i++ )
  {
    theBinnings[i].setMinMaxN( min_arg, max_arg, n);
  }
}

void binningLogInRanges::add(const double rangeParameter, const double x)
{
  for ( unsigned int i = 0; i < theBinnings.size(); i++ )
  {
    if ( rangeParameter >= lowValues[i] && rangeParameter <= highValues[i] )
    {
      theBinnings[i].add( x );
    }
  }
}

void binningLogInRanges::getBorders(const int _dim, double& lowBorder, double& highBorder)
{
  lowBorder = lowValues[_dim];
  highBorder = highValues[_dim];
}

//----------------------------------------------------------------------------
