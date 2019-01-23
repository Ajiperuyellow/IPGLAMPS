//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/binning3.h $
//$LastChangedDate: 2016-04-19 20:13:24 +0200 (Di, 19. Apr 2016) $
//$LastChangedRevision: 2328 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file 
 * @brief Declaration of 3dim histogramms
 *
 */

#ifndef BINNING3_H
#define BINNING3_H

#include "binning.h"

#include <fstream>
#include <vector>
#include <string>

//----------------------------------------------------------------------------

/** 
 * @brief class to implement 3dim histograms
 *
 * Contrary to the 1dim case (see class #binning), we have not the
 * possibility to first store all values and do the binning at the
 * end. Here we have to know the binning area and number of bins (or
 * width) from the very beginning.
 */

class binning3d
{
public:
  //     binning();
  //     binning(const std::string);
  //
  //     binning(const double, const double, const int);
  //     binning3d ( const std::string name, const double min_x_arg, const double max_x_arg, const int n_x, const double min_y_arg, const double max_y_arg, const int n_y );
  binning3d();

  //     binning(const int);
  //     binning(const std::string, const int);
  //
  //     binning(const double);
  //     binning(const std::string, const double);
  //
  //     ~binning();

  //     void setWidth(const double);
  //     void setNumber(const int);
  void setFilename ( const std::string & name ) { filename = name; };

  void setMinMaxN( const double min_x_arg, const double max_x_arg, const int n_x, const double min_y_arg, const double max_y_arg, const int n_y, const double min_z_arg, const double max_z_arg, const int n_z );

  void setMinMaxWidth( const double min_x_arg, const double max_x_arg, const double width_x, const double min_y_arg, const double max_y_arg, const double width_y, const double min_z_arg, const double max_z_arg, const double width_z );

  //     double getWidth() const {return binWidth;}
  //     int getNBins() const {return nBins;}

  //     double getBinLabel(const int i);
  //     double getBin(const int i);
    
  std::vector< int > getIndex( const double x, const double y, const double z );
  int getBinRaw( const int i, const int j, const int k ) { return theBins[i][j][k]; };

  void add( const double x, const double y, const double z );
  void print();

private:
  void setFilename() { setFilename( "defaultBinning3.dat" ); };

  std::vector< std::vector< std::vector<int> > > theBins;

  double x_min, x_max, y_min, y_max, z_min, z_max, binWidth_x, binWidth_y, binWidth_z;
  int nBins_x, nBins_y, nBins_z;
  int nEntries;
  
  std::string filename;
};

//----------------------------------------------------------------------------

#endif
