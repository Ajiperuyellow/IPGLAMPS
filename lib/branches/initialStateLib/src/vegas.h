//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/vegas.h $
//$LastChangedDate: 2014-02-10 09:50:50 +0100 (Mo, 10. Feb 2014) $
//$LastChangedRevision: 1612 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifndef VEGAS_H
#define VEGAS_H

/**
 * @brief abstract class for integrands to pass to integration routines.
 *
 * The operator () is overloaded such that we have a functor object. 
 * Here we implement:
 * * void operator(): for usage with CUBA-Vegas 
 * * double operator(): for usage with NR-Vegas
 **/
class integrand
{
  public:
    /** 
     * @brief Overloaded operator() that makes integrand a functor object
     * - for use with CUBA-vegas
     *
     * just calls double operator() const
     **/ 
    virtual void operator()(const int *ndim, const double xx[], const int *ncomp, double ff[]) const
    {
      double wgt = 0;
      ff[0] = this->operator()( xx, wgt );
    };
    
    /** 
     * @brief Overloaded operator() that makes integrand a functor object 
     * - for use with NR-vegas
     **/ 
    virtual double operator()(const double [], double) const = 0;

    /**
     * @brief Destructor
     **/
    virtual ~integrand() {};
};


/**
 * @brief takes a function object derived from integrand as argument
 **/
void vegas(int ndim, integrand& fxn, double *tgral, double *sd, double *chi2a);

/**
 * @brief ???
 **/
void rebin(double rc, int nd, double r[], double xin[], double xi[]);

#endif
