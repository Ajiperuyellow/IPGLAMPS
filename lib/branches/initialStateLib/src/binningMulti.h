//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/binningMulti.h $
//$LastChangedDate: 2015-01-28 18:21:52 +0100 (Mi, 28. Jan 2015) $
//$LastChangedRevision: 2060 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file 
 * @brief Declaration of binning complex classes
 *
 */

#ifndef BINNINGMULTI_H
#define BINNINGMULTI_H

#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>

/**
 * @brief class to implement a templated 'binning'
 *
 * The class T has to provide:
 * - T += T
 * - T *= a
 * - operator<<   which prints the object in one single line
 *
 * Instead of T/a we use T * (1.0/a) -- due to speed
 *
 **/
template<typename T>
class tBinningMulti
{  
protected:
  double x_min, x_max, binWidth;
  unsigned int nBins;
  unsigned long nEntries;
  std::string filename;
  std::vector<T> theBins;
public:

  /** @brief Default constructor */
  tBinningMulti():
    x_min(1e20),
    x_max(-1e20),
    binWidth(),
    nBins(100),            // default number of bins
    nEntries(0),
    filename("tBinningMulti.dat"),
    theBins(nBins, T() )  // fill with 0
  { };

  /** 
   * @brief Constructor for specifying the binning area and the number
   * of bins
   *
   * @param[in] min_arg lower bound
   * @param[in] max_arg upper bound
   * @param[in] n number of bins
   **/
  tBinningMulti(const double min_arg, const double max_arg, const int n) :
    x_min(min_arg),
    x_max(max_arg),
    binWidth( (max_arg-min_arg)/n ),
    nBins(n),            // default number of bins
    nEntries(0),
    filename("tBinningMulti.dat"),
    theBins(nBins, T() )  // fill with 0
  { };

  /** 
   * @brief Constructor for specifying the binning area and the number
   * of bins + file name 
   *
   * @param[in] name file name for output
   * @param[in] min_arg lower bound
   * @param[in] max_arg upper bound
   * @param[in] n number of bins
   **/
  tBinningMulti(const std::string & name, const double min_arg, const double max_arg, const int n):
    x_min(min_arg),
    x_max(max_arg),
    binWidth( (max_arg-min_arg)/n ),
    nBins(n),            // default number of bins
    nEntries(0),
    filename(name),
    theBins(nBins, T() )  // fill with 0
  { };

  /** 
   * @brief Constructor for specifying the binning area and the bin width
   *
   * @param[in] min_arg lower bound
   * @param[in] max_arg upper bound
   * @param[in] w width of the bins
   **/
  tBinningMulti(const double min_arg, const double max_arg, const double w) :
    x_min(min_arg),
    x_max(max_arg),
    binWidth( w ),
    nBins( ceil( (max_arg-min_arg)/w ) ),            // default number of bins
    nEntries(0),
    filename("tBinningMulti.dat"),
    theBins(nBins, T() )  // fill with 0
  { };

  /** 
   * @brief Constructor for specifying the binning area and the bin width
   *
   * @param[in] name file name for output
   * @param[in] min_arg lower bound
   * @param[in] max_arg upper bound
   * @param[in] w width of the bins
   **/
  tBinningMulti(const std::string & name, const double min_arg, const double max_arg, const double w) :
    x_min(min_arg),
    x_max(max_arg),
    binWidth( w ),
    nBins( ceil( (max_arg-min_arg)/w ) ),            // default number of bins
    nEntries(0),
    filename(name),
    theBins(nBins, T() )  // fill with 0
  { };

  /** @brief Add an element to the binned data */
  void add(const double x, const T val)
  {
    int index = int((x - x_min)/binWidth);
    if (index >= nBins)
    {
      index = nBins-1;
    }
    else if (index < 0)
    {
      index = 0;
    }
    theBins[index] += val;
    ++nEntries;
  };

  /** @brief Print the histogram */
  void print(std::fstream &file)
  {
    file << "# bin width: " << binWidth 
         << "   number entries: " << nEntries << std::endl;

    const double fak = 1.0/binWidth;
    for (int i=0; i<nBins; i++)
    {
      file << x_min + (i+0.5)*binWidth << "\t " 
           << theBins[i]*fak << std::endl;
    } 
  };

  /** @brief Print the histogram */
  void print()
  {
    std::fstream file(filename.c_str(),std::ios::out);
    print(file);
    file.close();
  }

  /** 
   * @brief Print the histogram 
   *
   * here all bins are additionally normalize to number of entries
   **/
  void printNormed(std::fstream &file)
  {
    file << "# bin width: " << binWidth 
         << "   number entries: " << nEntries << std::endl;

    const double fak = 1.0/(binWidth*nEntries);
    for (int i=0; i<nBins; i++)
    {
      file << x_min + (i+0.5)*binWidth << "\t " 
           << theBins[i]*fak << std::endl;
    } 
  };

  /** 
   * @brief Print the histogram 
   *
   * here all bins are additionally normalize to number of entries
   **/
  void printNormed()
  {
    std::fstream file(filename.c_str(),std::ios::out);
    printNormed(file);
    file.close();
  }

};

#endif
