//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/interpolation_n_dimensions.h $
//$LastChangedDate: 2014-11-16 22:00:01 +0100 (So, 16. Nov 2014) $
//$LastChangedRevision: 1938 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file 
 * @brief Declarations for the class interpolation1d.
 */


#ifndef INTERPOLATION_n_D_H
#define INTERPOLATION_n_D_H

#include <vector>
#include <stdexcept>
#include <string>


enum DATA_ERROR_CODE {DATA_SUCCESS, DATA_READOUT_ERROR, DATA_FILE_ERROR, DATA_SIZE_ERROR}; 


class interpolation_n_D
{
public:
  ~interpolation_n_D();
    
  /** @brief Set filename of data file */
  void setFilename( std::string filename_arg ) { filename = filename_arg; };

  /** @brief Configure the interpolation routine by giving the main parameters */
  void configure();
    
  /** @brief find value via interpolation from tabulated values */
  double getInterpolatedData() const;
    
    
protected:
  /** @brief This routine reads the table from a file */
  DATA_ERROR_CODE readTables();

  /** @brief the data structure holding the tabulated values */
  std::vector<double> data;
    
  /** @brief total number of entries to be expected in the input file */
  unsigned int expected_entries;

  /** @brief data file holding the tabulated values, see below for required structure */
  std::string filename;
};


class interpolation1d : public interpolation_n_D
{
public:
  /** @brief Configure the interpolation routine by giving the main parameters */
  void configure( const int n_i_arg, const double delta_a_arg,const double a_start_arg, std::string filename_arg = "default.dat" );
    
  /** @brief find value via interpolation from tabulated values */
  double getInterpolatedData( const double a ) const;
    
protected:
  /** @brief Utility routine to get the index number */
  int get_index(const int) const;
    
  /** @brief Reserve memory and load the data. */
  void loadData();
    
  /** @brief get interpolated data using linear interpolation */
  double linearInterpolation(const double a) const;   

  /** @brief number of entries in a-direction */
  int n_i;
    
  /** @brief spacing of tabulated a values */
  double delta_a;

  /** @brief lowest tabulated value in a-direction */
  double a_start;
};


class interpolation2d : public interpolation_n_D
{
public:
  /** @brief Configure the interpolation routine by giving the main parameters */
  void configure( const int n_i_arg, const int n_j_arg, const double delta_a_arg, const double delta_b_arg, const double a_start_arg, const double b_start_arg, std::string filename_arg = "default.dat" );
    
  /** @brief find value via interpolation from tabulated values */
  double getInterpolatedData( const double a, const double b ) const;
    
protected:
  /** @brief Utility routine to get the index number */
  int get_index(const int, const int) const;
    
  /** @brief Reserve memory and load the data. */
  void loadData();
    
  /** @brief get value using linear interpolation */
  double linearInterpolation(const double, const double) const;  

  /** @brief number of entries in a direction */
  int n_i;
  /** @brief number of entries in b direction */
  int n_j;
    
  /** @brief spacing of tabulated a values */
  double delta_a;
  /** @brief spacing of tabulated b values */
  double delta_b;

  /** @brief lowest tabulated value in a-direction */
  double a_start;
  /** @brief lowest tabulated value in b-direction */
  double b_start;
};


class interpolation3d : public interpolation_n_D
{
public:
  /** @brief Configure the interpolation routine by giving the main parameters */
  void configure( const int n_i_arg, const int n_j_arg, const int n_k_arg, const double delta_a_arg, const double delta_b_arg, const double delta_c_arg, const double a_start_arg, const double b_start_arg, const double c_start_arg, std::string filename_arg = "default.dat" );
  
  /** @brief find value via interpolation from tabulated values */
  double getInterpolatedData( const double a, const double b, const double c ) const;
  
protected:
  /** @brief Utility routine to get the index number */
  int get_index(const int, const int, const int) const;
  
  /** @brief Reserve memory and load the data. */
  void loadData();
  
  /** @brief get value using linear interpolation */
  double linearInterpolation(const double, const double, const double) const;  
  
  /** @brief number of entries in a direction */
  int n_i;
  /** @brief number of entries in b direction */
  int n_j;
  /** @brief number of entries in c direction */
  int n_k;
  
  /** @brief spacing of tabulated a values */
  double delta_a;
  /** @brief spacing of tabulated b values */
  double delta_b;
  /** @brief spacing of tabulated c values */
  double delta_c;
  
  /** @brief lowest tabulated value in a-direction */
  double a_start;
  /** @brief lowest tabulated value in b-direction */
  double b_start;
  /** @brief lowest tabulated value in c-direction */
  double c_start;
};



/** @brief exception class for handling unexpected behaviour when reading the file containing the IgQ tables */
class eData_read_error : public std::runtime_error
{
public:
  explicit eData_read_error(const std::string& what) : std::runtime_error(what) {};

  virtual ~eData_read_error() throw() {};
};



#endif
