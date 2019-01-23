//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/initialmodel.h $
//$LastChangedDate: 2013-10-19 00:21:06 +0200 (土, 19 10月 2013) $
//$LastChangedRevision: 1500 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------



#ifndef INITIALMODEL_H
#define INITIALMODEL_H

#include <vector>
#include <iostream>
#include <boost/smart_ptr.hpp>

#include "particle.h"
#include "configuration.h"
#include "woodsaxon.h"
#include "rangen_distr.h"
#include "vegas.h"


typedef boost::shared_ptr<ranGen_Distr> tPointerToRanGen;


/** 
 * @brief abstract base class for all models that generate initial
 * particle distributions 
 */
class initialModel
{
public:
  /** 
   * @brief constructor 
   */
  initialModel( ) {};

  /** 
   * @brief Reserve memory for the Particle vector and sample momenta
   * and positions. 
   *
   * Pure virtual function. Any derived class must at least specify
   * this routine.
   */
  virtual void populateParticleVector( std::vector<Particle>& _particles ) = 0;

  /**
   * @brief routine to set all Particle::unique_id
   *
   * Loops over all particles and sets the unique ID of every particle.
   */
  void setUniqueID( std::vector<Particle>& _particles );

  /**
   * @brief Return a 'typical' radius of the model
   */
  virtual double Radius();

};



/** 
 * @brief abstract base class for all models that generate initial
 * particle distributions, extended by Wood-Saxon information
 */
class initialModelWS : public initialModel
{
public:
  WoodSaxon WoodSaxonParameter;

  /** 
   * @brief constructor 
   */
  initialModelWS( const config& _config);

  /** 
   * @brief destructor 
   */
  ~initialModelWS();

  /**
   * @brief Return a 'typical' radius of the model
   */
  virtual double Radius();

  /** 
   * @brief Reserve memory for the Particle vector and sample momenta
   * and positions. 
   */
  virtual void populateParticleVector( std::vector<Particle>& _particles ) = 0;

  /** 
   * @brief nuclear density (Woods-Saxon distribution) 
   *
   * This is just a shortcut to WoodSaxonParameter.densityA
   */
  double densityA(double b, double z) const
  { return WoodSaxonParameter.densityA(b,z); };

  /** 
   * @brief Sample XYZ according Wood-Saxon (for one particle) 
   *
   * This routine samples the X,Y,Z coordinates according the
   * Wood-Saxon parameters
   *
   * @param[in] T The time coordinate
   * @param[out] X,Y,Z The spatial coordinates
   **/
  void sample_XYZ_WoodSaxon(double T, double &X, double &Y, double &Z) const;

  /** 
   * @brief Sample time T according Wood-Saxon (for one particle) 
   * 
   * @param[out] T The time coordinate
   **/
  void sample_T_WoodSaxon(double &T) const;

  /** 
   * @brief Sample time and coordinate position for a given single particle
   * 
   * @param[in,out] _particle The particle for which t,x,y,z should be sampled
   **/
  void sample_TXYZ_singleParticle( Particle& _particle );

  /**
   * @brief Generate the distribution function for the time variable
   *
   * @param[put] T_AB the overlap
   **/
  void generateTimeDistributionWS(double & T_AB);

protected:
  double A;              /**< mass number of nucleus A */   
  double Aatomic;        /**< atomic number, i.e. number of protons, of nucleus A */
  double B;              /**< mass number of nucleus B */   
  double Batomic;        /**< atomic number, i.e. number of protons, of nucleus B */
  double impactParameter;/**< impact parameter (in fm)*/
  double sqrtS_perNN;    /**< Energy (in GeV)*/

  tPointerToRanGen distrTime; /**< the random generator for the time distribution */
};


/** 
 * @brief exception class for handling unexpected critical behaviour
 * within generation of initial distributions  
 */
class eInitialModel_error : public std::runtime_error
{
public:
  explicit eInitialModel_error(const std::string& what) : std::runtime_error(what) {};
  
  virtual ~eInitialModel_error() throw() {};
};


class integrand_time : public integrand
{
public:
  integrand_time() {};
  integrand_time(const double b, const double t) : bImp(b),time(t) {};
  ~integrand_time() {};
  
  /** 
   * @brief Overloaded operator() that makes integrand_time a functor
   * object - for use with CUBA-vegas 
   **/ 
  void operator()(const int*, const double [], const int*, double []) const;

  /** 
   * @brief Overloaded operator() that makes integrand_time a
   * functor object - for use with NR-vegas 
   **/  
  double operator()(const double [], double) const;
  
  void setB(const double b) { bImp = b; }
  void setTime(const double t) { time = t; }
  void setWoodSaxonParameter( const WoodSaxon& _w ) { woodSaxonParameter = _w; }
  
private:
  double bImp;
  double time;
  WoodSaxon woodSaxonParameter;
};




#endif // INITIALMODEL_H


// kate: indent-mode cstyle; space-indent on; indent-width 2; replace-tabs on; 
