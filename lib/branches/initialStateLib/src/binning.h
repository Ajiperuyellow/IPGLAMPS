//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/binning.h $
//$LastChangedDate: 2017-10-27 18:21:49 +0200 (Fr, 27. Okt 2017) $
//$LastChangedRevision: 2637 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file
 * @brief Declaration of 1dim histogramms
 *
 */

#ifndef BINNING_H
#define BINNING_H

#include <algorithm>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>

//----------------------------------------------------------------------------

/**
 * @brief Type of binning
 */
enum class BINTYPE
{
  dummy, dynamic, fixed_number, fixed_width, fixed_ALL
};

//----------------------------------------------------------------------------

/**
 * @brief class to implement 1dim histograms
 *
 * You may select range and bin width already at construction, or
 * leave this free til the end. In this case, all entries are first
 * stored into a list and min/max are calculated at the end and
 * allocation and filling of the bins is also done at the end
 */
class binning
{
public:

  /** @brief Default constructor */
  binning();

  /** @brief Default constructor + file name */
  binning(const std::string & name);

  /** @brief Constructor for specifying the number of bins */
  binning(const int n);

  /** @brief Constructor for specifying the number of bins + file name */
  binning(const std::string & name, const int n);

  /** @brief Constructor for specifying the bin width */
  binning(const double width);

  /** @brief Constructor for specifying the bin width + file name */
  binning(const std::string & name, const double width);

  /** @brief Constructor for specifying the binning area and the number of bins */
  binning(const double min_arg, const double max_arg, const int n);

  /** @brief Constructor for specifying the binning area and the number of bins + file name */
  binning(const std::string & name, const double min_arg, const double max_arg, const int n);


  /** @brief Destructor */
  ~binning();

  /**
   * @brief Set the filename for output
   *
   * @param[in] name filename
   */
  void setFilename(const std::string & name) { filename = name; }

  /** @brief Set the bin width and the type to #fixed_width */
  void setWidth(const double);

  /** @brief Set the number of bins and the type to #fixed_number */
  void setNumber(const int);

  /** @brief Set the binning area and the number of bins and the type to #fixed_min_max_number */
  void setMinMaxN(const double x_min, const double x_max, const int n);

  /** @brief Set the binning area and the width of the bins and the type to #fixed_min_max_number */
  void setMinMaxWidth(const double x_min, const double x_max, const double width_arg);

  /** @brief Return the bin width */
  double getWidth() const {return binWidth;}

  /** @brief Return the number of bins */
  int getNBins() const {return nBins;}

  /** @brief Return the label of some bin */
  double getBinLabel(const int i);

  /** @brief Return the value of some bin */
  double getBin(const int i);

  /** @brief Return the value of some bin without scaling with the bin width */
  int getBinRaw(const int i);
  
  /** @brief Return the value of some bin, divided by number of all binned events */
  double getBinRelative( const int i);

  /** @brief Add an element to the binned data */
  void add(const double);
  
  /** @brief Print the histogram */
  void print();

  /** @brief Print the histogram */
  void print(std::fstream &file);

  /** @brief Return the binned data as vectors */
  void getVectors(std::vector<double>& x, std::vector<double>& y, double& width);

  /**
   * @brief The Increment Operator: adding two binning objects
   *
   * every component is added, min/max values are corrected.
   *
   * It is *not* checked, whether the parameters (e.g. binwidth,
   * nBins,...) are the same
   **/
  inline binning & operator += (const binning & a)
  {
    x_min = std::min(x_min, a.x_min);
    x_max = std::max(x_max, a.x_max);
    // binwidth
    // nBins
    nEntries += a.nEntries;

    // http://stackoverflow.com/questions/3177241/what-is-the-best-way-to-concatenate-two-vectors
    if (a.theList.size() > 0)
    {
      theList.reserve( theList.size() + a.theList.size() );
      theList.insert( theList.end(), a.theList.begin(), a.theList.end() );
    }

    // http://stackoverflow.com/questions/3376124/how-to-add-element-by-element-of-two-stl-vectors
    std::transform(theBins.begin( ), theBins.end( ), a.theBins.begin( ), theBins.begin( ),std::plus<double>( ));

    // finished
    // btype
    // filename

  }

private:

  /**
   * @brief Do the actual binning
   */
  void finish();

  /**
   * @brief Do the actual binning, if necessary
   */
  inline void doFinish(void) { if (!finished) finish(); };

private:
  double x_min;
  double x_max;
  double binWidth;
  int nBins;
  int nEntries;

  std::vector<double> theList;
  std::vector<int> theBins;
  
  bool finished;

  BINTYPE bintype;

  std::string filename;
};

//----------------------------------------------------------------------------
/**
 * @brief class to implement 1dim histograms with logarithmic x-scale
 *
 * You may select range and bin width already at construction, or
 * leave this free til the end. In this case, all entries are first
 * stored into a list and min/max are calculated at the end and
 * allocation and filling of the bins is also done at the end
 */
class binningLog
{
public:
  /** @brief Default constructor */
  binningLog();

  /** @brief Default constructor + file name */
  binningLog(const std::string);

  /** @brief Constructor for specifying the binning area and the number of bins */
  binningLog(const double min_arg, const double max_arg, const int n);

  /** @brief Constructor for specifying the binning area and the number of bins + file name*/
  binningLog(const std::string name, const double min_arg, const double max_arg, const int n);

  /** @brief Destructor */
  ~binningLog();

  void setNumber(const int);
  void setFilename(const std::string name) {filename = name;}

  void setMinMaxN(const double x_min, const double x_max, const int n);

  int getNBins() const {return nBins;}

  double getBinLabel(const int i);
  double getBin(const int i);
  int getBinRaw(const int i);
  double getBinRelative( const int i);

  void add(const double);
  void print();
  void getVectors(std::vector<double>& x, std::vector<double>& y, double& width);

private:
  void setFilename() {setFilename("defaultBinning.dat");}
  void finish();

  /**
   * @brief Do the actual binning, if necessary
   */
  inline void doFinish(void) { if (!finished) finish(); };


private:
  double x_min;
  double x_max;
  double binWidth;
  int nBins;
  int nEntries;

  std::vector<double> theList;
  std::vector<int> theBins;

  bool finished;

  BINTYPE bintype;

  std::string filename;
};

//----------------------------------------------------------------------------

/**
 * @brief class to implement a 1dim histogram, where for every bin the
 * average of added values is reported.
 *
 * This class reports the mean value and the statistical error for
 * each bin.
 *
 * Contrary to the (usual) 1dim case (see class #binning), we have not the
 * possibility to first store all values and do the binning at the
 * end. Here we have to know the binning area and number of bins (or
 * width) from the very beginning.
 */
class binningValues
{
public:

  /** @brief Default constructor */
  binningValues();

  /** @brief Constructor for specifying the binning area and the number of bins + file name */
  binningValues( const std::string & name, const double min_arg, const double max_arg, const int n);

  /** @brief Constructor for specifying the binning area and the width of bins + file name */
  binningValues( const std::string & name, const double min_arg, const double max_arg, const double width);

  /**
   * @brief Set the filename for output
   *
   * @param[in] name filename
   */
  void setFilename(const std::string & name) { filename = name; }

  /** @brief Set the binning area and the number of bins and the type to #fixed_min_max_number */
  void setMinMaxN( const double x_min, const double x_max, const int n );

  /** @brief Set the binning area and the width of the bins and the type to #fixed_min_max_number */
  void setMinMaxWidth( const double x_min, const double x_max, const double width_arg );

  /** @brief Return the bin width */
  double getWidth() const {return binWidth;}

  /** @brief Return the number of bins */
  int getNBins() const {return nBins;}

  /** @brief Return the label of some bin */
  double getBinLabel(const int i);

  /** @brief Return the number of values inserted in some bin */
  int getBinCount( const int i );

  /** @brief Return the sum of all values inserted in some bin */
  long double getBinValue( const int i );

  /** @brief Return the sum of all squared values inserted in some bin */
  long double getBinValueSquared( const int i );

  /** @brief Return the value of the statistical error of some bin */
  double getBinError( const int i );

  /** @brief Return the sqrt of the variance of some bin */
  double getSigma( const int i );

  /** @brief Return the sqrt of the variance divided by number of entries of some bin */
  double getBinMean( const int i );

  /** @return some correction factor */
  double getAreaNormalization( const int i );

  /** @brief Add an element to the binned data */
  void add(const double x, const double v);
  
  /** @brief Print the histogram */
  void print();

  /** @brief Print the histogram */
  void print(std::fstream &file);

  /**
   * @brief The Increment Operator: adding two binning objects
   *
   * every component is added, min/max values are corrected.
   *
   * It is *not* checked, whether the parameters (e.g. binwidth,
   * nBins,...) are the same
   **/
  inline binningValues & operator += (const binningValues & a)
  {
    x_min = std::min(x_min, a.x_min);
    x_max = std::max(x_max, a.x_max);
    // binwidth
    // nBins
    nEntries += a.nEntries;

    // http://stackoverflow.com/questions/3376124/how-to-add-element-by-element-of-two-stl-vectors
    std::transform(theBins.begin( ), theBins.end( ), a.theBins.begin( ), theBins.begin( ),std::plus<double>( ));
    std::transform(theBinsSquared.begin( ), theBinsSquared.end( ), a.theBinsSquared.begin( ), theBinsSquared.begin( ),std::plus<double>( ));
    std::transform(theBinsN.begin( ), theBinsN.end( ), a.theBinsN.begin( ), theBinsN.begin( ),std::plus<int>( ));

    // finished
    // btype
    // filename

  }


private:
  /**
   * @brief Set the filename for output
   *
   * Sets the value to a default value.
   */
  void setFilename() {setFilename("defaultBinningValues.dat");}

private:
  double x_min;
  double x_max;
  double binWidth;
  int nBins;
  int nEntries;

  std::vector<long double > theBins;
  std::vector<long double > theBinsSquared;
  std::vector<int > theBinsN;

  std::string filename;
};


/**
 * @brief class to implement a 1dim histogram, where for every bin the
 * average of added values is reported.
 *
 * This class reports the mean vaule and the statistical error for
 * each bin.
 *
 * Contrary to the (usual) 1dim case (see class #binning), we have not the
 * possibility to first store all values and do the binning at the
 * end. Here we have to know the binning area and number of bins (or
 * width) from the very beginning.
 */
class binningLogValues
{
public:

  /** @brief Default constructor */
  binningLogValues();

  /** @brief Constructor for specifying the binning area and the number of bins + file name */
  binningLogValues( const std::string & name, const double min_arg, const double max_arg, const int n);

  /** @brief Constructor for specifying the binning area and the width of bins + file name */
  binningLogValues( const std::string & name, const double min_arg, const double max_arg, const double width);

  /**
   * @brief Set the filename for output
   *
   * @param[in] name filename
   */
  void setFilename(const std::string & name) { filename = name; }

  /** @brief Set the binning area and the number of bins and the type to #fixed_min_max_number */
  void setMinMaxN( const double x_min, const double x_max, const int n );

  /** @brief Return the number of bins */
  int getNBins() const {return nBins;}

  /** @brief Return the label of some bin */
  double getBinLabel(const int i);

  /** @brief Return the number of values inserted in some bin */
  int getBinCount ( const int i );

  /** @brief Return the sum of all values inserted in some bin */
  long double getBinValue( const int i );

  /** @brief Return the sum of all squared values inserted in some bin */
  long double getBinValueSquared( const int i );

  /** @brief Return the value of the statistical error of some bin */
  double getBinError( const int i );

  /** @brief Return the sqrt of the variance of some bin */
  double getSigma( const int i );

  /** @brief Return the sqrt of the variance divided by number of entries of some bin */
  double getBinMean( const int i );

  /** @return some correction factor */
  double getAreaNormalization( const int i );

  /** @brief Add an element to the binned data */
  void add(const double x, const double v);

  /** @brief Print the histogram */
  void print();

  /** @brief Print the histogram */
  void print(std::fstream &file);

private:
  /**
   * @brief Set the filename for output
   *
   * Sets the value to a default value.
   */
  void setFilename() {setFilename("defaultBinningValues.dat");}

  double x_min;
  double x_max;
  double binWidth;
  int nBins;
  int nEntries;

  std::vector<long double > theBins;
  std::vector<long double > theBinsSquared;
  std::vector<int > theBinsN;

  std::string filename;  
};

//----------------------------------------------------------------------------

class binningInRanges
{
public:

  binningInRanges();
  binningInRanges( const int _dim );
  ~binningInRanges();

  void setNumberRanges( const int _dim );
  void setMinMaxN( const double min_arg, const double max_arg, const int n);
  void setMinMaxWidth( const double min_arg, const double max_arg, const double width_arg );
  void setFilename( const std::string & name ) {filename = name;}

  void setRangevalues( const int _dim, const double minValue, const double maxValue );

  double getBinLabel(const int i){ return theBinnings[0].getBinLabel( i ); };
  double getBin(const int _dim, const int i){ return theBinnings[_dim].getBin( i ); };
  int getBinRaw(const int _dim, const int i){ return theBinnings[_dim].getBinRaw( i ); };
  double getBinRelative(const int _dim, const int i){ return theBinnings[_dim].getBinRelative( i ); };

  double getWidth() const { return theBinnings[0].getWidth(); }
  int getNBins() const { return theBinnings[0].getNBins(); }

  int getDimension() { return theBinnings.size(); };
  void getBorders( const int _dim, double &lowBorder, double &highBorder );
  void add( const double rangeParameter, const double x );

  void print();

private:

  std::vector< binning > theBinnings;
  std::vector< double > lowValues;
  std::vector< double > highValues;

  std::string filename;
};


class binningLogInRanges
{
public:

  binningLogInRanges();
  binningLogInRanges( const int _dim );
  ~binningLogInRanges();

  void setNumberRanges( const int _dim );
  void setMinMaxN( const double min_arg, const double max_arg, const int n);
  void setFilename( const std::string & name ) {filename = name;}

  void setRangevalues( const int _dim, const double minValue, const double maxValue );

  double getBinLabel(const int i){ return theBinnings[0].getBinLabel( i ); };
  double getBin(const int _dim, const int i){ return theBinnings[_dim].getBin( i ); };
  int getBinRaw(const int _dim, const int i){ return theBinnings[_dim].getBinRaw( i ); };
  double getBinRelative(const int _dim, const int i){ return theBinnings[_dim].getBinRelative( i ); };

  int getNBins() const { return theBinnings[0].getNBins(); }

  int getDimension() { return theBinnings.size(); };
  void getBorders( const int _dim, double &lowBorder, double &highBorder );
  void add( const double rangeParameter, const double x );

  void print();

private:

  std::vector< binningLog > theBinnings;
  std::vector< double > lowValues;
  std::vector< double > highValues;

  std::string filename;
};




//----------------------------------------------------------------------------

/** @brief exception class for handling unexpected critical behaviour within the binning class  */
class eBINNING_error : public std::runtime_error
{
public:
  explicit eBINNING_error ( const std::string& what ) : std::runtime_error ( what ) {};

  virtual ~eBINNING_error() throw() {};
};

//----------------------------------------------------------------------------

// WARNING: if 'initializer( omp_priv = omp_orig )' is used, make sure
// that the histogram is initialized, but empty at the beginning of
// the reduction. Otherwise, the already existing entries are
// multiplied by the number of threads!

#pragma omp declare                          \
  reduction(+ : binning :                    \
            omp_out += omp_in )              \
  initializer( omp_priv = omp_orig )


//initializer( omp_priv = binning() )

#pragma omp declare                          \
  reduction(+ : binningValues :              \
            omp_out += omp_in )              \
  initializer( omp_priv = omp_orig )


//initializer( omp_priv = binningValues() )



#endif
