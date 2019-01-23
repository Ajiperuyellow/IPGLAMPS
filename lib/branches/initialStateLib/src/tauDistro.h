//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/tauDistro.h $
//$LastChangedDate: 2016-03-23 21:04:01 +0100 (Mi, 23. Mär 2016) $
//$LastChangedRevision: 2319 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file
 * @brief Declaration of calculation distributions for constant
 * eigentime tau
 *
 * At the moment, we implement the generation of
 *   dN/deta_s, dE_T/deta_s, dN/dy, dE_T/dy,
 * which are all binned simultanously in one loop over the particle
 * vector
 */

#ifndef TAUDISTRO_H
#define TAUDISTRO_H

#include <cmath>
#include <iostream>
#include <vector>
#include <string>

#include "binning.h"
#include "binning2.h"

/**
 * @brief class to implement array of histograms for binning as
 * function of eta_s at constant tau
 */
class tTauDistro
{
public:

  /**
   * @brief Constructor
   *
   * @param[in] tstep the values of tau, when stuff should be binned
   * @param[in] fName0 path and base name of output filenames
   * @param[in] vMax maximum value of etas/y to consider in binning
   * @param[in] vDelta bin width of etas/y to consider in binning
   */
  tTauDistro(const std::vector<double> & tstep, const std::string & fName0,
             const double vMax, const double vDelta);

  /**
   * @brief Destructor
   *
   * The destructor automatically calls the print routine.
   */
  ~tTauDistro()
  {
    if(toPrint) print();
  };

  /**
   * @brief do the binning at a specific time
   *
   * @param[in] t value of t
   * @param[in] dt value of timestep size
   * @param[in] parts particle vector
   *
   * we use templates here, since the particle vector may be a vector
   * of different classes in the different branches/codes
   **/
  template<typename T>
  void add(const double t, const double dt, std::vector<T> particles);

  /**
   * @brief write out all the stuff
   **/
  void print(void);

  /**
   * @brief generate a filename for a given timestep
   *
   * for step==0 and step==nTimeSteps, it uses the strings 'initial'
   * and 'final', otherwise it uses the numerical value
   **/
  std::string filename_step( const int step ) const;


protected:

  std::vector<double> taustep; ///< the values of the timestep

  std::string fNamePre; ///< contains path and base name of output filenames

  std::vector<binningValues> dNdeta; ///< the binnings of dN/deta_s
  std::vector<binningValues> dEdeta; ///< the binnings of dE_T/deta_s
  std::vector<binningValues> dNdy; ///< the binnings of dN/dy
  std::vector<binningValues> dEdy; ///< the binnings of dE_T/dy


private:
  bool toPrint; ///< internal flag: whether print in destructor


};

template<typename T>
void tTauDistro::add(const double t, const double dt, std::vector<T> particles)
{

#pragma omp parallel for
  for( auto it = begin(particles); it < end(particles); ++it )
  {
    const double tau = it->Pos.Eigentime();
    const double etas = it->Pos.Rapidity();
    const double dTau = dt*sqrt(1-pow(tanh(etas),2));

    for (int iTau = 0; iTau < taustep.size(); ++iTau)
    {
      const double & Tau = taustep[iTau];

      if (tau <= Tau-dTau) continue;
      if (tau > Tau) continue;

      const double pt = it->Mom.Pt();
      const double y = it->Mom.Rapidity();

#pragma omp critical(add_dNdeta)
      dNdeta[iTau].add(etas, 1.0);
#pragma omp critical(add_dEdeta)
      dEdeta[iTau].add(etas, pt);

#pragma omp critical(add_dNdy)
      dNdy[iTau].add(y, 1.0);
#pragma omp critical(add_dEdy)
      dEdy[iTau].add(y, pt);

      break; // loop finished!
    }

  }
  toPrint = true;
}

#endif
