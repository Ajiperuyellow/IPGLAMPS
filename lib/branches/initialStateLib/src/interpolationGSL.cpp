//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/interpolationGSL.cpp $
//$LastChangedDate: 2015-02-25 21:03:59 +0100 (Mi, 25. Feb 2015) $
//$LastChangedRevision: 2089 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include "interpolationGSL.h"

#include <iomanip>      // std::setprecision
#include <gsl/gsl_errno.h>

using std::cout;
using std::endl;



interpolationGSL::interpolationGSL(const std::vector<double>& x, const std::vector<double>& y, const INTERP_TYPE interp_type) : acc(), interp(), X(), Y(), nPoints(0), name("noname")
{
  if ( ( nPoints  = x.size() ) != y.size() )
  {
    std::string errMsg = "Size of arrays not equal.";
    throw interpolationGSL_error( errMsg );
  }

  setParamGSL(interp_type);

  X = new double[nPoints];
  Y = new double[nPoints];

  std::copy(x.begin(),x.end(), X);
  std::copy(y.begin(),y.end(), Y);

  gsl_interp_init(interp, X, Y, nPoints);
}

interpolationGSL::interpolationGSL(const double x[], const double y[], const int N, const INTERP_TYPE interp_type) : acc(), interp(), X(), Y(), nPoints(N), name("noname")
{
  setParamGSL(interp_type);

  X = new double[nPoints];
  Y = new double[nPoints];

  std::copy(x,x+nPoints, X);
  std::copy(y,y+nPoints, Y);

  gsl_interp_init(interp, X, Y, nPoints);
}

interpolationGSL::interpolationGSL(const std::string & Name, const std::vector<double>& x, const std::vector<double>& y, const INTERP_TYPE interp_type) : acc(), interp(), X(), Y(), nPoints(0), name(Name)
{
  if ( ( nPoints  = x.size() ) != y.size() )
  {
    std::string errMsg = "Size of arrays not equal.";
    throw interpolationGSL_error( errMsg );
  }

  setParamGSL(interp_type);

  X = new double[nPoints];
  Y = new double[nPoints];

  std::copy(x.begin(),x.end(), X);
  std::copy(y.begin(),y.end(), Y);

  gsl_interp_init(interp, X, Y, nPoints);
}

interpolationGSL::interpolationGSL(const std::string & Name, const double x[], const double y[], const int N, const INTERP_TYPE interp_type) : acc(), interp(), X(), Y(), nPoints(N), name(Name)
{
  setParamGSL(interp_type);

  X = new double[nPoints];
  Y = new double[nPoints];

  std::copy(x,x+nPoints, X);
  std::copy(y,y+nPoints, Y);

  gsl_interp_init(interp, X, Y, nPoints);
}


void interpolationGSL::setParamGSL(const INTERP_TYPE interp_type)
{
  acc = gsl_interp_accel_alloc();
  switch (interp_type)
  {
    case interp_linear: 
    {
      interp = gsl_interp_alloc(gsl_interp_linear, nPoints);
    }  break;
    case interp_polynomial: 
    {
      interp = gsl_interp_alloc(gsl_interp_polynomial, nPoints); 
    } break;
    case interp_cspline: 
    {
      interp = gsl_interp_alloc(gsl_interp_cspline, nPoints);  
    } break;
    default:
    {
      std::string errMsg = "interp_type unknown.";
      throw interpolationGSL_error( errMsg );
    }
      
  }
}

interpolationGSL::~interpolationGSL()
{
  gsl_interp_free(interp);
  gsl_interp_accel_free(acc);
  delete[] X;
  delete[] Y;
}

double interpolationGSL::eval(const double x) const
{
  double y;

  if ( (gsl_interp_eval_e(interp, X, Y, x, acc, &y)) )
  {
    y = 0.0; // Error occured, most probably x out of range.
  }

  return y;
}

double interpolationGSL::Deriv(const double x) const
{
  double y;

  if (gsl_interp_eval_deriv_e(interp, X, Y, x, acc, &y))
  {
    y = 0.0; // Error occured, most probably x out of range.
  }

  return y;
}

double interpolationGSL::Integ(const double x) const
{
  double y;

  if (gsl_interp_eval_integ_e(interp, X, Y, X[0], x, acc, &y))
  {
    y = 0.0; // Error occured, most probably x out of range.
  }

  return y;
}

void interpolationGSL::PrintDataPoints(std::ostream &f, const int xPrecision)
{
  f << "# " << name << endl;
  f << "# Data Points: " << nPoints << endl;
  //  f << " precision:: "<< f.precision() << endl;
  for (unsigned int i=0; i<nPoints;i++)
  {
    f << std::scientific << std::setprecision(xPrecision) 
      << X[i] << "  " << Y[i] << endl;
  }
}

void interpolationGSL::PrintInterpolation(std::ostream &f, const int N, const int xPrecision)
{
  double dx = (X[nPoints-1]-X[0])/N;
  double x;
  int j;

  f << "# Interpolation: " << N << " Subpoints" << endl;

  for (unsigned int i=0; i<nPoints-1;i++)
  {
    dx = (X[i+1]-X[i])/N;
    for (j=0, x=X[i]; j<N; j++, x+=dx)
    {
      f << std::scientific << std::setprecision(xPrecision) 
	<< x << "  " << eval(x) << endl;
    }
  }
  f << std::scientific << std::setprecision(xPrecision) 
    << X[nPoints-1] << "  " << eval(X[nPoints-1]) << endl;
}
