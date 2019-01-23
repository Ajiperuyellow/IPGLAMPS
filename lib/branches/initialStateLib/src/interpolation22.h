//--------------------------------------------------------- -*- c++ -*- ------
//$HeadURL: file:///home/bamps/svn/lib/branches/initialStateLib/src/interpolation22.h $
//$LastChangedDate: 2015-06-09 12:35:17 +0200 (Di, 09. Jun 2015) $
//$LastChangedRevision: 2169 $
//$LastChangedBy: senzel $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/** @file 
 * @brief Declarations for the class interpolation22.
 */


#ifndef INTERPOLATION22_H
#define INTERPOLATION22_H

#include <math.h>
#include <boost/lexical_cast.hpp>

#include "globalsettings.h"
#include "interpolation_n_dimensions.h"
#include "particleprototype.h"
#include "FPT_compare.h"


/** @brief exception class for handling unexpected behaviour when reading the file containing the IgQ tables */
class eI22_read_error : public std::runtime_error
{
public:
  explicit eI22_read_error(const std::string& what) : std::runtime_error(what) {};

  virtual ~eI22_read_error() throw() {};
};


class interpolation22_process
{
public:
  /** @brief Standard constructor. When used, configure() needs to be called prior to other methods! */
  interpolation22_process( const std::string & _process_name = "default" ) :
    process_name( _process_name ),
    // standard values for massless particles
    n_ln_s ( 108 ),
    n_ln_md2g_wo_as ( 104 ),
    n_ln_md2q_wo_as ( 34 ),
    delta_ln_s ( 0.1 ),
    delta_ln_md2g_wo_as ( 0.1 ),
    delta_ln_md2q_wo_as ( 0.3 ),
    ln_s_start ( log( 1.1 * ns_casc::lambda2 ) ),
    ln_md2g_wo_as_start( log( 0.001 / 0.3 ) ),
    ln_md2q_wo_as_start( log( 0.0002 / 0.3 ) )
  { };
  ~interpolation22_process() { };
    
  /** @brief get the interpolated cross section of the process interactionType, which are tabled in interactiontype.h */
  double get_cs( const double s, const double md2g_wo_as_arg = 0.0, const double md2q_wo_as_arg = 0.0 ) const { return 0.0; };
    
protected:
  /** @brief get filename of data file of the given process. */
  std::string getFilename( const bool asRun_arg, const int N_light_flav, const bool is_heavy_quark_involved = false, const double maximum_coupling = 1.0 ) const
  {
    std::string filename;
    if( is_heavy_quark_involved )
      filename = "data/cs22/external_data/cs_" + process_name;
    else
      filename = "data/cs22/cs_" + process_name;
    
    if( asRun_arg )
    {
      std::string nameAsMax = "";
      if( FPT_COMP_E( maximum_coupling, 1.0 ) )
        nameAsMax = "";
      else if( FPT_COMP_E( maximum_coupling, 0.3 ) )
        nameAsMax = "_asMax_03";
      else if( FPT_COMP_E( maximum_coupling, 0.6 ) )
        nameAsMax = "_asMax_06";
      else if( FPT_COMP_E( maximum_coupling, 0.8 ) )
        nameAsMax = "_asMax_08";
      else if( FPT_COMP_E( maximum_coupling, 1.5 ) )
        nameAsMax = "_asMax_15";
      else if( FPT_COMP_E( maximum_coupling, 2.0 ) )
        nameAsMax = "_asMax_2";
      else
        throw eI22_read_error( "I22 tables for given maximum coupling not calculated.");
        
      filename = filename + "_asRun_Nf_" + boost::lexical_cast<std::string>( N_light_flav ) + nameAsMax;
    }
    else
      filename = filename + "_asConst";
      
    filename = filename + "_table.dat";
      
    return filename;
  }

  /** @brief Name of of the given process, eg. "gg_gg". */
  std::string process_name;
    
  /** @brief number of entries in ln(s) direction */
  int n_ln_s;
  /** @brief number of entries in ln(md2g_wo_as) direction */
  int n_ln_md2g_wo_as;
  /** @brief number of entries in ln(md2q_wo_as) direction */
  int n_ln_md2q_wo_as;
    
  /** @brief spacing of tabulated ln(s) values */
  double delta_ln_s;
  /** @brief spacing of tabulated ln(md2g_wo_as) values */
  double delta_ln_md2g_wo_as;
  /** @brief spacing of tabulated ln(md2q_wo_as) values */
  double delta_ln_md2q_wo_as;

  /** @brief lowest tabulated value in ln(s) direction */
  double ln_s_start;
  /** @brief lowest tabulated value in ln(md2g_wo_as) direction */
  double ln_md2g_wo_as_start;
  /** @brief lowest tabulated value in ln(md2q_wo_as) direction */
  double ln_md2q_wo_as_start;
};


class interpolation22_process_1d : public interpolation22_process
{
public:
  /** @brief get the interpolated cross section of the process interactionType, which are tabled in interactiontype.h */
  double get_cs( const double s, const double md2g_wo_as_arg = 0.0, const double md2q_wo_as_arg = 0.0 ) const;
    
  /** @brief Load table */
  void loadTable( const bool asRun, const int N_light_flav, const bool is_heavy_quark_involved = false )
  {
    I22.configure( n_ln_s, delta_ln_s, ln_s_start, getFilename( asRun, N_light_flav, is_heavy_quark_involved ) );
  };
    
protected:
  /** @brief Actual interpolation routine with stored data. */
  interpolation1d I22;
};


class interpolation22_process_2d : public interpolation22_process
{
public:
  interpolation22_process_2d( const std::string & _process_name = "default" ) : interpolation22_process( _process_name ), I22() {};
  /** @brief get the interpolated cross section of the process interactionType, which are tabled in interactiontype.h */
  double get_cs( const double s, const double md2g_wo_as_arg, const double md2q_wo_as_arg = 0.0 ) const;
    
  /** @brief Load table */
  void loadTable( const bool asRun, const int N_light_flav, const bool is_heavy_quark_involved = false, const double maximum_coupling = 1.0 )
  {
    I22.configure( n_ln_md2g_wo_as, n_ln_s, delta_ln_md2g_wo_as, delta_ln_s, ln_md2g_wo_as_start, ln_s_start, getFilename( asRun, N_light_flav, is_heavy_quark_involved, maximum_coupling ) );
  };
    
protected:
  /** @brief Actual interpolation routine with stored data. */
  interpolation2d I22;
};


class interpolation22_process_3d : public interpolation22_process
{
public:
  interpolation22_process_3d( const std::string & _process_name = "default" ) : interpolation22_process( _process_name ), I22() {};
  /** @brief get the interpolated cross section of the process interactionType, which are tabled in interactiontype.h */
  double get_cs( const double s, const double md2g_wo_as_arg, const double md2q_wo_as_arg ) const;
  
  /** @brief Load table */
  void loadTable( const bool asRun, const int N_light_flav, const bool is_heavy_quark_involved = false, const double maximum_coupling = 1.0 )
  {
    I22.configure( n_ln_md2q_wo_as, n_ln_md2g_wo_as, n_ln_s, delta_ln_md2q_wo_as, delta_ln_md2g_wo_as, delta_ln_s, ln_md2q_wo_as_start, ln_md2g_wo_as_start, ln_s_start, getFilename( asRun, N_light_flav, is_heavy_quark_involved, maximum_coupling ) );
  };
  
protected:
  /** @brief Actual interpolation routine with stored data. */
  interpolation3d I22;
};

// light quarks
class interpolation22_gg_gg : public interpolation22_process_2d
{
public:
  interpolation22_gg_gg() : interpolation22_process_2d( "gg_gg" ) {};
};

class interpolation22_gg_qqbar : public interpolation22_process_3d
{
public:
  interpolation22_gg_qqbar() : interpolation22_process_3d( "gg_qqbar" ) {};
};

class interpolation22_qg_qg : public interpolation22_process_3d
{
public:
  interpolation22_qg_qg() : interpolation22_process_3d( "qg_qg" ) {};
};

class interpolation22_qqbar_qqbar : public interpolation22_process_2d
{
public:
  interpolation22_qqbar_qqbar() : interpolation22_process_2d( "qqbar_qqbar" ) {};
};

class interpolation22_qqbar_qqbarDash : public interpolation22_process_2d
{
public:
  interpolation22_qqbar_qqbarDash() : interpolation22_process_2d( "qqbar_qqbarDash" ) {};
};

class interpolation22_qq_qq : public interpolation22_process_2d
{
public:
  interpolation22_qq_qq() : interpolation22_process_2d( "qq_qq" ) {};
};

class interpolation22_qqdash_qqdash : public interpolation22_process_2d
{
public:
  interpolation22_qqdash_qqdash() : interpolation22_process_2d( "qqDash_qqDash" ) {};
};


class interpolation22_gg_ccbar : public interpolation22_process_3d
{
public:
  /** @brief Set parameters like n_ln_s, delta_ln_s and ln_s_start. */
  void setParametersAndLoadTable( const bool asRun_arg, const int N_light_flav, const double M_charm );
};


class interpolation22_gg_bbbar : public interpolation22_process_3d
{
public:
  /** @brief Set parameters like n_ln_s, delta_ln_s and ln_s_start. */
  void setParametersAndLoadTable( const bool asRun_arg, const int N_light_flav, const double M_bottom );
};


class interpolation22_qqbar_ccbar : public interpolation22_process_2d
{
public:
  /** @brief Set parameters like n_ln_s, delta_ln_s and ln_s_start. */
  void setParametersAndLoadTable( const bool asRun_arg, const int N_light_flav, const double M_charm );
};


class interpolation22_qqbar_bbbar : public interpolation22_process_2d
{
public:
  /** @brief Set parameters like n_ln_s, delta_ln_s and ln_s_start. */
  void setParametersAndLoadTable( const bool asRun_arg, const int N_light_flav, const double M_bottom );
};



class interpolation22_gc_gc : public interpolation22_process_3d
{
public:
  /** @brief Set parameters like n_ln_s, delta_ln_s and ln_s_start. */
  void setParametersAndLoadTable( const bool asRun_arg, const int N_light_flav, const double M_charm );
};


class interpolation22_gb_gb : public interpolation22_process_3d
{
public:
  /** @brief Set parameters like n_ln_s, delta_ln_s and ln_s_start. */
  void setParametersAndLoadTable( const bool asRun_arg, const int N_light_flav, const double M_bottom );
};


class interpolation22_qc_qc : public interpolation22_process_2d
{
public:
  /** @brief Set parameters like n_ln_s, delta_ln_s and ln_s_start. */
  void setParametersAndLoadTable( const bool asRun_arg, const int N_light_flav, const double M_charm );
};


class interpolation22_qb_qb : public interpolation22_process_2d
{
public:
  /** @brief Set parameters like n_ln_s, delta_ln_s and ln_s_start. */
  void setParametersAndLoadTable( const bool asRun_arg, const int N_light_flav, const double M_bottom );
};




class interpolation22
{
public:
  /** @brief Standard constructor. When used, configure() needs to be called prior to other methods! */
  interpolation22() { };
  /** @brief Constructor that also sets the parameters */
  interpolation22( const bool asRun, const int N_light_flav, const int N_heavy_flav = 0 ) { configure( asRun, N_light_flav, N_heavy_flav ); };
  ~interpolation22() { };
    
  /** @brief Sets parameters and and initializes interpolation routines of given processes */
  void configure( const bool asRun, const int N_light_flav, const int N_heavy_flav = 0, const double M_c = ParticlePrototype::Mcharm, const double M_b = ParticlePrototype::Mbottom, const double maximum_coupling = 1.0, const double fixed_coupling_value = 0.3 );
    
  /** @brief get the interpolated cross section of the process interactionType, which are tabled in interactiontype.h */
  double get_cs( const int interactionType, const double s, const double md2g_wo_as_arg = 0.0, const double md2q_wo_as_arg = 0.0 ) const;
    
private:
    
  interpolation22_gg_gg I_gg_gg;
  interpolation22_gg_qqbar I_gg_qqbar;
  interpolation22_qg_qg I_qg_qg;
  interpolation22_qqbar_qqbar I_qqbar_qqbar;
  interpolation22_qqbar_qqbarDash I_qqbar_qqbarDash;
  interpolation22_qq_qq I_qq_qq;
  interpolation22_qqdash_qqdash I_qqdash_qqdash;

  interpolation22_gg_ccbar I_gg_ccbar;
  interpolation22_gg_bbbar I_gg_bbbar;
  interpolation22_qqbar_ccbar I_qqbar_ccbar;
  interpolation22_qqbar_bbbar I_qqbar_bbbar;
    
  interpolation22_gc_gc I_gc_gc;
  interpolation22_gb_gb I_gb_gb;
  interpolation22_qc_qc I_qc_qc;
  interpolation22_qb_qb I_qb_qb;
};




#endif
