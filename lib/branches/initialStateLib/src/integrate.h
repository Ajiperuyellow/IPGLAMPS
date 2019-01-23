//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/integrate.h $
//$LastChangedDate: 2014-05-30 15:30:33 +0200 (Fr, 30. Mai 2014) $
//$LastChangedRevision: 1733 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifndef INTEGRATE_H
#define INTEGRATE_H

// include config file that is created at CMake configuration time 
#include "configBAMPS.h"
#ifdef CUBA_FOUND
  #include <cuba.h>
#endif

#include "vegas.h"
#include "integrand32.h"

/** @brief typedef for pointer to function that finds peaks of integrands, needed for Divonne algorithm */
typedef void (*peakfinder_t)(const int *, const double [], int *, double []);


/**
 * @brief Wrapper class for use with CUBA integration routines
 *
 * Static declarations so that the variables and the member function can be used without actually creating a wrapper object
 *- > little trick to get rid of global scope variables, implement object oriented structure and allow for splitting into different
 * translation units
 * 
 * Usage: 1. Set integrand_cuba_wrapper::theIntegrand_in_wrapper to the actual integrand object.
 * 2. Use integrand_cuba_wrapper::_integrand where a function pointer of type integrand_t (see cuba.h) is expected 
 */
class integrand_cuba_wrapper
{
public:
  /** @brief Pointer to the actual integrand object */
  static integrand* theIntegrand_in_wrapper;
    
  /** @brief Wrapper function around the integrand functionoid, use where integrand_t (see cuba.h) is expected */
  static int _integrand(const int *ndim, const double xx[], const int *ncomp, double ff[], void* userdata)  
  { 
    (*theIntegrand_in_wrapper)(ndim,xx,ncomp,ff); 
    return 0;
  }
    
};


#ifdef CUBA_FOUND
/**
* @brief Generic prototype for function objects to be used with Cuba library routines
*
* The basic idea is the following:
* In order to use this additional wrapper layer, create an object of a class derived from integrate_cuba_generic with the
* appropriate constructor arguments. Then just call the overloaded operator() (as if it was a function) to get the integration
* result (as integral[0]). The interface operator() should be the same for each method of numerical integration, only the constructor
* should change.
*/
class integrate_cuba_generic
{
public:
  integrate_cuba_generic( const int _ndim, const int _ncomp, const double _epsrel, const double _epsabs, const int _flags, 
                          const int _mineval, const int _maxeval ) :
    ndim(_ndim), ncomp(_ncomp), epsrel(_epsrel), epsabs(_epsabs), flags(_flags), mineval(_mineval), maxeval(_maxeval) {};
                        
  void operator()( integrand& _integrand, int& neval, int& fail, double integral[], double error[], double prob[] ) {};   
  //     virtual void operator()( integrand_t integrand, int& neval, int& fail, double integral[], double error[], double prob[] ); 
    
    
protected:
  int ndim;
  int ncomp;
  double epsrel;
  double epsabs;
  int flags;
  int mineval;
  int maxeval;
};


/** @brief Integration with the Vegas routine from the Cuba library */
class integrate_vegas : public integrate_cuba_generic
{
public:
  integrate_vegas( const int _ndim, const int _ncomp, const double _epsrel, const double _epsabs, const int _flags,
                   const int _mineval, const int _maxeval, const int _nstart, const int _nincrease, const int _nbatch ) :
    integrate_cuba_generic( _ndim, _ncomp, _epsrel, _epsabs, _flags, _mineval, _maxeval),
    nbatch(_nbatch), nstart(_nstart), nincrease(_nincrease) {};
  
  void operator()( integrand& _integrand, int& neval, int& fail, double integral[], double error[], double prob[] )
  {
    const int seed = 0;
    const int gridno = 0;
    integrand_cuba_wrapper::theIntegrand_in_wrapper = &_integrand;
    Vegas( ndim, ncomp, integrand_cuba_wrapper::_integrand, 0, epsrel, epsabs, flags, seed,
           mineval, maxeval, nstart, nincrease, 
           nbatch, gridno, "",
           &neval, &fail, integral, error, prob );
  }
  
    
private:
  int nbatch;
  int nstart;
  int nincrease;   
};



/** @brief Integration with the Suave routine from the Cuba library */
class integrate_suave : public integrate_cuba_generic
{
public:
  integrate_suave( const int _ndim, const int _ncomp, const double _epsrel, const double _epsabs, const int _flags,
                   const int _mineval, const int _maxeval, const int _nnew, const double _flatness ) :
    integrate_cuba_generic( _ndim, _ncomp, _epsrel, _epsabs, _flags, _mineval, _maxeval),
    nregions(),
    nnew(_nnew), 
    flatness(_flatness) {};
  
  void operator()( integrand& _integrand, int& neval, int& fail, double integral[], double error[], double prob[] )
  {
    const int seed = 0;
    integrand_cuba_wrapper::theIntegrand_in_wrapper = &_integrand;
    Suave( ndim, ncomp, integrand_cuba_wrapper::_integrand, 0, epsrel, epsabs, flags, seed,
           mineval, maxeval, nnew, flatness, &nregions, &neval, &fail, integral, error, prob );
  }
  
  int nregions;
  
private:
  int nnew;
  double flatness;
};



/** @brief Integration with the Divonne routine from the Cuba library */
class integrate_divonne : public integrate_cuba_generic
{
public:
  integrate_divonne( const int _ndim, const int _ncomp, const double _epsrel, const double _epsabs, const int _flags,
                     const int _mineval, const int _maxeval, const int _key1, const int _key2, const int _key3,
                     const int _maxpass, const int _border, const double _maxchisq, const double _mindeviation,
                     const int _ngiven, const int _ldxgiven, double _xgiven[], const int _nextra,
                     peakfinder_t _peakfinder, const int _nregions ) :
    integrate_cuba_generic( _ndim, _ncomp, _epsrel, _epsabs, _flags, _mineval, _maxeval),
    nregions(_nregions),
    key1(_key1), key2(_key2), key3(_key3),
    maxpass(_maxpass), border(_border), maxchisq(_maxchisq), mindeviation(_mindeviation),
    ngiven(_ngiven), ldxgiven(_ldxgiven), xgiven(_xgiven), nextra(_nextra),
    peakfinder(_peakfinder) {};
  
  void operator()( integrand& _integrand, int& neval, int& fail, double integral[], double error[], double prob[] )
  {
    const int seed = 0;
    integrand_cuba_wrapper::theIntegrand_in_wrapper = &_integrand;
    Divonne(ndim, ncomp, integrand_cuba_wrapper::_integrand, 0, epsrel, epsabs, flags, seed, 
            mineval, maxeval, key1, key2, key3, maxpass, border, maxchisq, mindeviation, 
            ngiven, ldxgiven, xgiven, nextra, peakfinder, &nregions, &neval, &fail, integral, error, prob);
  }
  
  int nregions;
  
private:
  int key1, key2, key3;
  int maxpass;
  double border;
  double maxchisq;
  double mindeviation;
  int ngiven; 
  int ldxgiven;
  double * xgiven;
  int nextra;
  peakfinder_t peakfinder;
};
#endif  //#ifdef CUBA_FOUND


/** @brief Integration with the Vegas routine from Numerical Recipes */
class integrate_nr_vegas
{
public:
  integrate_nr_vegas(const int _ndim) : ndim(_ndim) {};
  
  void operator()( integrand& _integrand , int& neval, int& fail, double integral[], double error[], double prob[] )
  {
    double tgral, sd, chi2a;
    
    vegas( ndim, _integrand, &tgral, &sd, &chi2a );
    
    neval = -1;
    fail = 0;
    integral[0] = tgral;
    error[0] = sd;
    prob[0] = chi2a;
  }
  
                       
private:
  int ndim;
};


#endif
