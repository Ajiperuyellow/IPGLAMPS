//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/interpolation23.h $
//$LastChangedDate: 2014-12-11 22:00:00 +0100 (Do, 11. Dez 2014) $
//$LastChangedRevision: 2015 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file 
 * @brief Declarations for the class interpolation23.
 */


#ifndef INTERPOLATION23_H
#define INTERPOLATION23_H

#include <vector>
#include <stdexcept>
#include <string.h>


enum I23_ERROR_CODE {I23_SUCCESS, I23_READOUT_ERROR, I23_FILE_ERROR, I23_SIZE_ERROR}; 
enum INTERPOLATION_MODE {log_interpol, lin_interpol, mix_interpol};


class interpolation23
{
public:
  // <<-------------------------------------------------
  // a wrapper for old routines to make new routines for heavy quarks with running coupling and more choices in formationTimeTyp, matrix element, etc. compatible with the old ones, try not to use anymore
  interpolation23( const bool read_I23_data_file = true ) { configure( !read_I23_data_file, 1, 0, 1, "bamps_org", "GBimproved", true, 1.0, log_interpol ) ; };
  // --------------------------------------------------->>

  interpolation23( const bool I23onlineIntegration_arg, 
                   const int id_hq_arg, 
                   const double M_hq_arg, 
                   const double kappa_arg, 
                   const std::string & formationTimeTyp_arg, 
                   const std::string & matrixElement23_arg, 
                   const bool md2_counter_term_in_I23_arg, 
                   const double fudge_factor_lpm_arg, 
                   const INTERPOLATION_MODE interpolation_mode_arg, 
                   const bool matrixElement23_22qt_arg = false );// wrapper

  ~interpolation23();
    
  /** @brief initializises framework, stores information about size and starting points of data files */
  void configure( const bool I23onlineIntegration_arg, 
                  const int id_hq_arg, 
                  const double M_hq_arg, 
                  const double kappa_arg, 
                  const std::string & formationTimeTyp_arg, 
                  const std::string & matrixElement23_arg, 
                  const bool md2_counter_term_in_I23_arg, 
                  const double fudge_factor_lpm_arg, 
                  const INTERPOLATION_MODE interpolation_mode_arg, 
                  const bool matrixElement23_22qt_arg = false );
  
  /** @brief find I23 via interpolation from tabulated values */
  double getI23( const double a, const double b, const double ca, const double d, const double e, bool & trustable ) const;
  /** @brief returns true when interpolation can be performed with the given parameters, false otherwise */
  //     bool canInterpolate(const double, const double, const double, const double) const;
  
  void setInterpolationMode( INTERPOLATION_MODE _x ) { interpolation_mode = _x; };
  
private:
  /** @brief This routine reads the I23 tables from a file */
  I23_ERROR_CODE readTables();
  int get_index( const int i, const int j, const int k, const int l, const int m ) const;
  int get_index_old_massless_i23_table( const int i, const int j, const int k, const int l ) const;
  bool inRange_old_massless_i23_table(const double, const double, const double, const double) const;
  //     int getShift(const std::vector<double>&) const;
  
  /** @brief get I23 using linear interpolation */
  double getI23_linearInterpolation( const double a, const double b, const double c, const double d, const double e, bool & trustable ) const;
  /** @brief get I23 using linear interpolation (works only for old massless i23 table for constant coupling) */
  double getI23_linearInterpolation_old_massless_i23_table( const double a, const double b, const double c, const double d ) const;
  
  /** @brief get I23 using polynomial interpolation */
  //     double getI23_polyInterpolation(const double, const double, const double, const double) const;
    
  /** @brief get I23 using spline interpolation */
  //     double getI23_splineInterpolation(const double, const double, const double, const double) const;
    
  /** @brief polynomial interpolation and extrapolation as found in numerical recipes */
  //     void polynomialInterpolation(const std::vector<double>& xa, const std::vector<double>& ya, const int n, const double x, 
  //                                  double& y, double& dy) const;
    
  /** @brief the data structure holding the tabulated values */
  std::vector<double> I23data;        
    
    
  /** @brief number of entries in a = ln(md2) direction */
  int n_i;
  /** @brief number of entries in b = ln(lambda) direction */
  int n_j;
  /** @brief number of entries in c = ln(gamma) direction for lower values of ln(gamma) = 0 .. 3.4 */
  int n_k_low;
  /** @brief number of entries in c = ln(gamma) direction for higher values of ln(gamma) = 4.4 .. 7.4 */
  int n_k_high;
  /** @brief total number of entries in c = ln(gamma) direction */
  int n_k;   
  /** @brief number of entries in d = cos(theta) direction */
  int n_l;
  /** @brief number of entries in e = M direction */
  int n_m;
    
  /** @brief total number of entries to be expected in the input file */
  unsigned int expected_entries;
    
  /** @brief spacing of tabulated a values */
  double delta_a;
  /** @brief spacing of tabulated b values */
  double delta_b;
  /** @brief spacing of tabulated c values for small c */
  double delta_c_low;
  /** @brief spacing of tabulated c values for large c*/
  double delta_c_high;
  /** @brief spacing of tabulated d values */
  double delta_d;
  /** @brief spacing of tabulated e values */
  double delta_e;
    
  /** @brief lowest tabulated value in a-direction */
  double a_start;
  /** @brief lowest tabulated value in b-direction */
  double b_start;
  /** @brief lowest tabulated value in c-direction */
  double c_start;
  /** @brief lowest tabulated value in d-direction */
  double d_start;
  /** @brief lowest tabulated value in e-direction */
  double e_start;
    
  /** @brief boundary between small and high c-values with a different spacing */
  double c_separator;
    
  /** @brief a 'very low' negative number */
  int minus_infinity;
    
  /** @brief true if I23 integration is performed online. In this case there is no need to read in table. */
  bool I23onlineIntegration;
    
  /** @brief whether heavy quark is particle 1 or 2 */
  int id_hq;
    
  /** @brief Mass of incoming particle (can be 0 or finite), the other incoming particle is finite anyhow */
  double M_hq;
    
  /** @brief Stores the mode of interpolation: "lin_interpol": linear interpolation in I23, "log_interpol": (linear) interpolation in log(I23) (original mode), "mix_interpol": superposition of both methods */
  INTERPOLATION_MODE interpolation_mode;
    
  /** @brief type of formation time of radiated gluon in 2->3 */
  std::string formationTimeTyp;
    
  /** @brief matrix element which is used for 2->3 heavy quark scattering */
  std::string matrixElement23;
    
  /** @brief Whether the gluon propagator in the 2->2 part of the radiative cross section is approximated by 1/qt^4. If false the more exact expression 1/t^2 is employed. Applies only to matrixElement23 = GBimproved. For GBorg always the approximated propagator is used. */
  bool matrixElement23_22qt;
    
  /** @brief Whether a counter term is applied which resembles the old prescription of Xu for masssless particles */
  bool md2_counter_term_in_I23;
    
  /** @brief factor for LPM effect to play around 
   * this factor is used to play around with the cutoff
   * instead of k_t > gamma / lambda a modified cutoff k_t > X gamma / lambda
   * is used, where X is the fudge factor
   */
  double fudge_factor_lpm;
    
  /** @brief kappa factor for Debye screening */
  double kappa;
    
  /** @brief whether the old version (before SVN) was used for creating the table (this had a different prefactor) */
  bool version_before_svn;
    
  /** @brief SVN revision which was used to create table */
  std::string svn_revision;
    
  /** @brief whether the old version for i23 table for massless particles with a constant coupling of 0.3 should be used */
  bool use_old_massless_i23_table;
    
};



/** @brief exception class for handling unexpected behaviour when reading the file containing the I23 tables */
class eI23_read_error : public std::runtime_error
{
public:
  explicit eI23_read_error(const std::string& what) : std::runtime_error(what) {};

  virtual ~eI23_read_error() throw() {};
};

#endif
