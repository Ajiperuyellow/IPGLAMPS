//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/binning2.h $
//$LastChangedDate: 2017-04-06 14:55:30 +0200 (Do, 06. Apr 2017) $
//$LastChangedRevision: 2556 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file 
 * @brief Declaration of 2dim histogramms
 *
 */

#ifndef BINNING2_H
#define BINNING2_H

#include "binning.h"

#include <fstream>
#include <vector>
#include <string>

//----------------------------------------------------------------------------

/** 
 * @brief class to implement 2dim histograms
 *
 * Contrary to the 1dim case (see class #binning), we have not the
 * possibility to first store all values and do the binning at the
 * end. Here we have to know the binning area and number of bins (or
 * width) from the very beginning.
 */
class binning2d
{
public:
  
  /** @brief Default constructor */
  binning2d();
  
  /** @brief Constructor for specifying the binning area and the number of bins + file name */
  binning2d(const std::string & name, 
            const double min_x_arg, const double max_x_arg, const int n_x, 
            const double min_y_arg, const double max_y_arg, const int n_y);

  /** @brief Constructor for specifying the binning area and the width of bins + file name */
  binning2d(const std::string & name, 
            const double min_x_arg, const double max_x_arg, const double width_x,
            const double min_y_arg, const double max_y_arg, const double width_y);
    
    
  /** 
   * @brief Set the filename for output
   * 
   * @param[in] name filename
   */
  void setFilename(const std::string & name) {filename = name;}


  /** @brief Set the binning area and the number of bins and the type to #fixed_min_max_number */
  void setMinMaxN( const double min_x_arg, const double max_x_arg, const int n_x, 
                   const double min_y_arg, const double max_y_arg, const int n_y );

  /** @brief Set the binning area and the width of the bins and the type to #fixed_min_max_number */
  void setMinMaxWidth( const double min_x_arg, const double max_x_arg, const double width_x, 
                       const double min_y_arg, const double max_y_arg, const double width_y );
    
  /** @brief Return the value of some bin without scaling with the bin width */
  int getBinRaw( const int i, const int j ) { return theBins[i][j]; };
  
  /** @brief Add an element to the binned data */
  void add(const double x, const double y);

  /** @brief Print the histogram 
   *
   * at the end, also the projections on the x- and the y-axis are
   * printed 
   **/
  void print();
    
  /** @brief Return the index of some bin */
  std::vector< int > getIndex( const double x, const double y );

  /** 
   * @brief The Increment Operator 
   *
   * This routine just addds the contents of the bins. No checking
   * of dimensions etc.
   **/
  binning2d & operator += (const binning2d &);

private:
  /** 
   * @brief Set the filename for output
   *
   * Sets the value to a default value.
   */
  void setFilename() {setFilename("defaultBinning2.dat");}
  
  std::vector< std::vector<int> > theBins;
  
  double x_min, x_max, y_min, y_max, binWidth_x, binWidth_y;
  int nBins_x, nBins_y;
  int nEntries;
  
  std::string filename;    
};

//----------------------------------------------------------------------------

/** 
 * @brief class to implement 2dim histograms, where for every bin the
 * average of added values is reported.
 **/
class binningValues2d
{
public:
  binningValues2d();
  //     binning(const std::string);
  //     
  //     binning(const double, const double, const int);
  binningValues2d(const std::string & name, const double min_x_arg, const double max_x_arg, const int n_x, const double min_y_arg, const double max_y_arg, const int n_y);
    
  /** @brief Destructor */
  ~binningValues2d();
  
  

  //     binning(const int);
  //     binning(const std::string, const int);
  //     
  //     binning(const double);
  //     binning(const std::string, const double);
  //     
  //     ~binning();
    
  //     void setWidth(const double);
  //     void setNumber(const int);
  void setFilename(const std::string & name) {filename = name;}
  void setMinMaxN( const double x_min, const double x_max, const int n );
  void setMinMaxN( const double min_x_arg, const double max_x_arg, const int n_x, const double min_y_arg, const double max_y_arg, const int n_y );
  void setMinMaxWidth( const double min_x_arg, const double max_x_arg, const double width_x, const double min_y_arg, const double max_y_arg, const double width_y );
  
  //     double getWidth() const {return binWidth;}
  //     int getNBins() const {return nBins;}

  double getBinLabelX( const int i );
  double getBinLabelY( const int i );
  double getBin( const int i, const int j );
  int getBinRaw( const int i, const int j );
  
  void add( const double x, const double y, const double v );
  void print();
  void print2DPlot();
  
private:
  void setFilename() {setFilename("defaultBinningValues2.dat");}
    
  //     std::vector<double> theList;
  //     std::vector<int> theBins;
  std::vector< std::vector<double> > theBins;
  std::vector< std::vector<int> > theBinsN;
    
  double x_min, x_max, y_min, y_max, binWidth_x, binWidth_y;
  int nBins_x, nBins_y;
  int nEntries;
  
  std::string filename;
};




#endif
