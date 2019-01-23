//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/interpolationGSL.h $
//$LastChangedDate: 2015-02-25 21:03:59 +0100 (Mi, 25. Feb 2015) $
//$LastChangedRevision: 2089 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifndef INTERPOLATIONGSL_H
#define INTERPOLATIONGSL_H

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>


/** @brief Enumeration type for possible interpolation types */
enum INTERP_TYPE { 
  interp_linear,
  interp_polynomial,
  interp_cspline
};

/** 
 * @brief This class encapsulates the interpolation routines of the GSL.
 *
 * We are using the interpolation object directly, since we store the
 * x- and y-values directly in our class. Therefore, the GSL spline
 * type would double the storage.
 */
class interpolationGSL
{
public:
  /** 
   * @brief Default Constructor 
   */
  interpolationGSL() : acc(), interp(), X(), Y(), nPoints(0), name("noname") {};
  
  /** 
   * @brief Constructor 
   */
  interpolationGSL(const std::vector<double>&, const std::vector<double>&, const INTERP_TYPE interp_type = interp_linear);
  
  /** 
   * @brief Constructor 
   */
  interpolationGSL(const double x[], const double y[], const int n, const INTERP_TYPE interp_type = interp_linear);
  
  /** 
   * @brief Constructor 
   */
  interpolationGSL(const std::string & Name, const std::vector<double>&, const std::vector<double>&, const INTERP_TYPE interp_type = interp_linear);
  
  /** 
   * @brief Constructor 
   */
  interpolationGSL(const std::string & Name, const double x[], const double y[], const int n, const INTERP_TYPE interp_type = interp_linear);
  

  /** 
   * @brief Destructor
   */
  ~interpolationGSL();
  
  /**
   * @brief access routine
   **/
  double eval(const double x) const;
  
  /**
   * @brief overloading the () operator
   **/
  double operator()(const double x) const { return eval(x); }
  
  /**
   * @brief returning the minimal x-value
   **/
  double minX() const { return X[0]; }
  
  /**
   * @brief returning the maximal x-value
   **/
  double maxX() const { return X[nPoints-1]; }
  
  /**
   * @brief return the derivative at point x
   **/
  double Deriv(const double x) const;

  /**
   * @brief return the integral over [x[0], x]
   **/
  double Integ(const double x) const;
  
  /**
   * @brief print the data points to some file
   **/
  void PrintDataPoints(std::ostream &f, const int xPrecision = 6);

  /**
   * @brief print the data points to some file
   *
   * Opens/closes the the file with the given name and print to it.
   *
   * The file name argument is not defined to be a reference, since we
   * want to have auto-conversion from 'const char [...]'. This allows
   * calling this method as '...PrintDataPoints("aaa.dat")'.
   **/
  void PrintDataPoints(std::string fName, const int xPrecision = 6)
  {
    std::fstream file;
    file.open(fName, std::fstream::out );
    PrintDataPoints(file,xPrecision);
    file.close();
  }

    
  /**
   * @brief print the interpolation curve to some file
   **/
  void PrintInterpolation(std::ostream &f, const int N, const int xPrecision = 6);
  
  
protected:
  /**
   * @brief set some internal GSL parameters during construction
   **/
  void setParamGSL(const INTERP_TYPE interp_type);
  
  
private:
  gsl_interp_accel *acc;
  gsl_interp *interp;
  
  double *X;
  double *Y;
  unsigned int nPoints;
  std::string name;
};

/** 
 * @brief exception class for handling unexpected critical behaviour within interpolation
 **/
class interpolationGSL_error : public std::runtime_error
{
public:
  explicit interpolationGSL_error(const std::string& what) : std::runtime_error(what) {};
  
  virtual ~interpolationGSL_error() throw() {};
};


#endif
