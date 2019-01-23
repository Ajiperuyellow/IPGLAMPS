//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/initialmodel_minijets.h $
//$LastChangedDate: 2013-10-19 00:21:06 +0200 (土, 19 10月 2013) $
//$LastChangedRevision: 1500 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


#ifndef INITIALMODEL_MINIJETS_H
#define INITIALMODEL_MINIJETS_H

#include <stdexcept>
#include <vector>
#include <string>

#include "initialmodel.h"
#include "configuration.h"
#include "particle.h"
#include "vegas.h"
#include "rangen_distr.h"
#include "pdfinterface.h"



/**
 * @brief Class to provide the minijet initialization
 */
class initialModel_minijets : public initialModelWS
{
public:
  initialModel_minijets( const config& _config );
  ~initialModel_minijets();
  
  void populateParticleVector( std::vector<Particle>& _particles );
  
private:
  
  /** 
   * @brief sampling of pT, y1,y2, and the flavours of the outgoing
   * particles.
   *
   * pT is set according dsigma/d(pT). Then, also pX and pY are set.
   *
   * In the second step, y1 and y2 ar eset according
   * dsigma/d(pT)d(y1)d(y2) for fixed pT via a rejection
   * method. Knowing y1 and y2, also pZ and E are known.
   *
   * In a third step, also the flavours of the outgoing partons are
   * set by a flavour decomposition of dsigma/d(pT)d(y1)d(y2).
   **/
  void sample_PXYZE_FLAV( std::vector<Particle>& _particles ) const;
  
  /** 
   * @brief sampling of collision times and positions of the parton
   *        pair
   **/
  void sample_TXYZ( std::vector<Particle>& _particles );
  
  /** 
   * @brief Plot the initial pT distribution
   * 
   * This routine is for testing purposes. It just prints out the
   * values of dsigma/d(pT)
   **/
  void Plot(void);
  
  /**
   * @brief Initilization of Time and pT distributions
   *
   * This routine is called once by the constructor.
   **/
  void generateSamplingDataSets();
  
  /** 
   * @brief Tabulate dsigma/d(pT)
   * 
   * This routines calculates dsigma/d(pT) and sets up the system for
   * random generation according this distribution.
   **/
  void generatePtDistribution( double& sigma_jet );
  

  double P0; /**< lower PT-cutoff (in GeV) */

  int nTestparticles;

  /** @brief Number of (high-pt) particles that is added on top of the reconstructed medium, using it as a background. */
  int nParticlesToGenerate;
  
  ranGen_Distr *distrPT; /**< the random generator for the PT distribution */
  interpolationGSL *maxIntegrandPT; /**< the maximum value of the integrand per pT */

  interfacePDF_generic *PDF; /**< the parton distribution generator */

};


/**
 * @brief Class to provide the function d(sigma)/d(pT)d(y1)d(y2) for
 * the integration over y1 and y2
 *
 * This class is used by the integration routines in order to
 * calculate d(sigma)/d(pT) by integrating over y1 and y2. While
 * sqrt(s) and pT are class variables, the values of y1 and y2 have to
 * be given as function argments to the operator().
 */
class integrand_distPT : public integrand
{
public:
  integrand_distPT() :
    pT(0.0), sqrtS(0.0), xT(0.0), PDF(0) {};

  integrand_distPT(const double _pT, const double _sqrtS, interfacePDF_generic *_PDF) : 
    pT(_pT), sqrtS(_sqrtS), PDF(_PDF) 
  {
    xT=(sqrtS==0.0)?0.0:2.0 * pT / sqrtS;
  };

  ~integrand_distPT() {};
  
  /** 
   * @brief Overloaded operator() that makes integrand_distPT a function
   * object - for use with CUBA-vegas 
   **/
  void operator()(const int*, const double [], const int*, double []) const;

  /** 
   * @brief Overloaded operator() that makes integrand_distPT a function
   * object - for use with NR-vegas 
   **/
  double operator()(const double [], double) const;
  
  /** 
   * @brief Set the pT value
   */
  void setPT(const double _pT) {pT = _pT; xT=(sqrtS==0.0)?0.0:2.0 * pT / sqrtS;}

  /** 
   * @brief Set the sqrtS value
   */
  void setSqrtS(const double _sqrtS) {sqrtS = _sqrtS; xT=(sqrtS==0.0)?0.0:2.0 * pT / sqrtS;}

  /** 
   * @brief Set the pointer to the PDF calss
   */
  void setPDF( interfacePDF_generic *_PDF) {PDF = _PDF;}

  /**
   * @brief the cross section d(sigma)/d(PT2)d(y1)d(y2)
   **/
  double calculate(const double Y1, const double Y2, const int n1, const int n2, double arr[13][13]) const;

  /**
   * @brief the maximum value of d(sigma)/d(PT2)d(y1)d(y2) for fixed pT
   **/
  double maxIntegrand( double& y ) const;

  /**
   * @brief the maximum value of d(sigma)/d(PT2)d(y1)d(y2) for fixed pT
   **/
  double getY1max(void) const;
  double getXT(void) const { return xT;};

private:
  double pT; /**< The stored value of pT */
  double sqrtS; /**< The stored value of sqrt(s) */
  double xT; /**< = 2.0 * pT / sqrtS */
  interfacePDF_generic *PDF; /**< pointer to the PDF generator */

  
  /**
   * @brief massless pQCD cross section: g g -> g g
   **/
  double gggg(const double t, const double u) const
  { 
    return (9./2.*(3.-t*u-u/(t*t)-t/(u*u))); 
  }

  /**
   * @brief massless pQCD cross section: g g -> q qbar
   */
  double ggqq(const double t, const double u) const
  { 
    return (1./6.*(t/u+u/t)-3./8.*(t*t+u*u)); 
  }

  /**
   * @brief massless pQCD cross section: q qbar -> g g
   */
  double qqgg(const double t, const double u) const
  { 
    return (64./9.*ggqq(t,u)); 
  }

  /**
   * @brief massless pQCD cross section: g q -> g q, g qbar -> g qbar
   **/
  double gqgq(const double t, const double u) const
  { 
    return (-4./9.*(1./u+u)+(1.+u*u)/(t*t)); 
  }

  /**
   * @brief massless pQCD cross section: q_a q_b -> q_a q_b
   * (t-channel gluon exchange)
   **/
  double qq12(const double t, const double u) const
  { 
    return (4./9.*(1.+u*u)/(t*t)); 
  }

  /**
   * @brief massless pQCD cross section: q_a q_a -> q_a q_a
   **/
  double qqqq(const double t, const double u) const
  { 
    return (4./9.*((1.+u*u)/(t*t)+(1.+t*t)/(u*u))-8./27./t/u); 
  }

  /**
   * @brief massless pQCD cross section: q_a qbar_a -> q_a qbar_a
   * (same flavour and antiflavour)
   **/
  double qaqs(const double t, const double u) const
  { 
    return (4./9.*((1.+u*u)/(t*t)+u*u+t*t)-8./27.*u*u/t); 
  }

  /**
   * @brief massless pQCD cross section: q_a qbar_a -> q_b qbar_b
   * (annihilation process)
   **/
  double qaqd(const double t, const double u) const
  { 
    return (4./9.*(t*t+u*u)); 
  }

};





/** 
 * @brief exception class for handling unexpected critical behaviour
 * within generation of mini-jet initial distributions 
 */
class eMiniJet_error : public std::runtime_error
{
  public:
    explicit eMiniJet_error(const std::string& what) : std::runtime_error(what) {};
    
    virtual ~eMiniJet_error() throw() {};
};


#endif // INITIALMODEL_MINIJETS_H
