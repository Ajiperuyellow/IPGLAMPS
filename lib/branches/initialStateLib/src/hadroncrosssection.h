//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/hadroncrosssection.h $
//$LastChangedDate: 2015-01-28 18:19:25 +0100 (Mi, 28. Jan 2015) $
//$LastChangedRevision: 2058 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file 
 * @brief Parametrization of hadron-hadron cross sections
 *
 * Some of these parametrizations are explicietely given in
 * literature, as e.g. the PDG versions HPR1R2 in the 2012 version:
 *   - J. Beringer et al. (Particle Data Group), Phys. Rev. D86, 010001(2012) 
 *
 * The second way is to parametrize more complicated parametrizations
 * in a simplier way, as e.g. with
 * \f$ y = A + B\ln(x) + C\ln^2(x) \f$ with \f$ x = \sqrt{s} \f$.
 * This is done here for the Pythia cross section parametrizations,
 * i.e. Schuler:1993wr:
 *   - G. Schuler and T. Sjostrand, Phys. Rev. D49, 2257 (1994).
 * 
 * Here all contributions were listed for the range 200GeV to 7000GeV
 * and fitted with gnuplot.
 *
 *
 * Another way would be to use the old, p_lab dependent
 * parametrisations of PDG(1994):
 *    - L. Montanet et al., Phys. Rev., D 50 (1994) 1173-1826 
 */

#ifndef HADRONCROSSSECTION_H
#define HADRONCROSSSECTION_H

//#include <iostream>
#include <math.h>


/**
 * @brief class to implement a pomeron-reggeon parametrization of a
 * cross section 
 *
 * This is according PDG, 
 * J. Beringer et al. (Particle Data Group), Phys. Rev. D86, 010001(2012) 
 *
 * The operator () is overloaded such that we have a functor object. 
 *
 * The template parameters have to be integers:
 *   - mA,mB,M: masses in MeV
 *   - P,R1,R2: in mb/100
 *   - eta1,eta2: 1/10000 (i.e. 4 digits)
 **/

template<unsigned int mA, unsigned int mB, unsigned int M, unsigned int P, unsigned int R1, signed int R2, unsigned int eta1, unsigned int eta2>
class XSgenerator_HPR1R2
{
  public:
    double operator()(const double sqrts)
    {
      const double H = M_PI*0.197*0.197*1000*1000*10/(M*M); // in GeV^2
      const double sM = pow( 0.001*(mA+mB+M), 2); // in GeV^2 !!!
      const double h = sqrts*sqrts/sM;

      // std::cout << "H = " << H << std::endl;
      // std::cout << "sM= " << sM << std::endl;
      // std::cout << "h = " << h << std::endl;

      // std::cout << "1 : " <<H*pow(log(h),2) << std::endl;
      // std::cout << "2 : " <<P*0.01 << std::endl;
      // std::cout << "3 : " <<R1 * 0.01 * pow( h, -(eta1*0.0001) ) << std::endl;
      // std::cout << "4 : " <<R2 * 0.01 * pow( h, -(eta2*0.0001) ) << std::endl;
      // std::cout << "eta1 = " << -(eta1*0.0001) << std::endl;

      return H*pow(log(h),2) + P*0.01 + R1 * 0.01 * pow( h, -(eta1*0.0001) ) + R2 * 0.01 * pow( h, -(eta2*0.0001) );

    }
};


/**
 * @brief class to implement some quadratic (fit-)function of log(sqrt(s))
 * 
 * This template class implements a parametric fit according
 * \f$ y = A + B\ln(x) + C\ln^2(x) \f$ with \f$ x = \sqrt{s} \f$.
 * 
 * The template parameters have to be integers:
 *    - A : 
 *
 * The fit parameters for the inelastic contributions are:
 * - 'non-diffractive':
 *   * a               = 29.1873          +/- 0.1814       (0.6214%)
 *   * b               = -3.86686         +/- 0.04913      (1.271%)
 *   * c               = 0.71339          +/- 0.003292     (0.4615%)
 *
 * - single diffractive (A B -> X B, A B -> A X), has to be taken twice
 *   * a               = 1.48719          +/- 0.006384     (0.4292%)
 *   * b               = 0.693689         +/- 0.001729     (0.2493%)
 *   * c               = -0.018363        +/- 0.0001159    (0.631%)
 *
 * - double diffractive (A B -> X X)
 *   * a               = -1.27515         +/- 0.009422     (0.7389%)
 *   * b               = 0.992487         +/- 0.002552     (0.2572%)
 *   * c               = 0.00758181       +/- 0.000171     (2.256%)
 *
 **/
template<long A, long B, long C>
class XSgenerator_fit
{
  public:
    double operator()(const double sqrts)
    {
      const double h = log(sqrts);
      return A*0.0001 + B*0.0001 * h + C*0.0001 * h * h;
    }
};

/**
 * @brief cross section parametrisation:  PDG, HPR1R2 for sqrt(s)>7GeV, pp_total
 **/
typedef XSgenerator_HPR1R2<938,938,2076,3373,1367,777,4120,5626> XSgenerator_pp_tot_PDG;

/**
 * @brief cross section parametrisation: Schuler/Sjostrand, pp_total
 *
 *   * a               = 38.6041          +/- 0.2427       (0.6286%)
 *   * b               = -2.63411         +/- 0.06574      (2.496%)
 *   * c               = 0.962335         +/- 0.004405     (0.4577%)
 **/
typedef XSgenerator_fit<386041,-26341,9623> XSgenerator_pp_tot_Pythia;

/**
 * @brief cross section parametrisation: Schuler/Sjostrand, pp_elast
 *
 *   * a               = 7.71752          +/- 0.06227      (0.8068%)
 *   * b               = -1.14711         +/- 0.01687      (1.47%)
 *   * c               = 0.278089         +/- 0.00113      (0.4064%)
 **/
typedef XSgenerator_fit<77175,-11471,2781> XSgenerator_pp_elast_Pythia;

/**
 * @brief cross section parametrisation: Schuler/Sjostrand, pp_inelast
 *
 *   * a               = 30.8865          +/- 0.184        (0.5959%)
 *   * b               = -1.48699         +/- 0.04986      (3.353%)
 *   * c               = 0.684246         +/- 0.003341     (0.4882%)
 **/
typedef XSgenerator_fit<308865,-14870,6842> XSgenerator_pp_inelast_Pythia;


#endif
