#include <math.h>

#include "coupling.h"
#include "particleprototype.h"
#include "globalsettings.h"

int coupling::Nflavor = ParticlePrototype::N_light_flavor;
bool coupling::isRunning_global = false;
double coupling::fixedCoupling_global = 0.3;
double coupling::alpha_max_global = 1.0; // maximum value of alpha_s

/**
* @brief Get the running coupling
* @param[in] Q2 scale for running coupling, e.g. momentum transfer of process
* @return coupling strength at this scale
*/
double coupling::get_running_coupling( const double Q2 )
{
  const double beta0 = 11.0 - 2.0 / 3.0 * Nflavor;
  const double Lambda2 = ns_casc::lambda2; // = (200 MeV)^2

  double alpha;
  
  if(Q2 < -Lambda2)
  {
    alpha = 4.0*M_PI/beta0 * 1.0/log(-Q2/Lambda2);
    if(alpha > alpha_max_global)
      alpha = alpha_max_global;
  }
  else if(Q2 > 0.0)
  {
    alpha = 4.0*M_PI/beta0 *  (0.5 - atan( log(Q2/Lambda2) / M_PI) / M_PI) ;
    if(alpha > alpha_max_global)
      alpha = alpha_max_global;
  }
  else
    alpha = alpha_max_global;
  
  return alpha;
}


/**
* @brief returns coupling strength for given scale
* @param[in] Q2 scale for running coupling, e.g. momentum transfer of process
* @param[in] isRunning whether the coupling should run or be constant
* @return coupling strength
*/
double coupling::get_coupling( const double Q2, const bool isRunning )
{
  double alpha;
  if( isRunning )
    alpha = get_running_coupling( Q2 );
  else
    alpha = get_constant_coupling();
  
  return alpha;
}


/**
* @brief returns coupling strength for given scale
* @param[in] Q2 scale for running coupling, e.g. momentum transfer of process
* @return coupling strength
*/
double coupling::get_coupling( const double Q2 )
{
  double alpha;
  if( isRunning_global )
    alpha = get_running_coupling( Q2 );
  else
    alpha = get_constant_coupling();
  
  return alpha;
}
