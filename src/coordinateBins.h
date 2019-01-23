//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/coordinateBins.h $
//$LastChangedDate: 2018-04-19 08:34:36 +0200 (Do, 19. Apr 2018) $
//$LastChangedRevision: 2723 $
//$LastChangedBy: greif $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


#ifndef COORDINATE_BINS_H
#define COORDINATE_BINS_H

#include <stdexcept>
#include <vector>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/nvp.hpp>

#include "particle.h"

using std::vector;


/**
 * @brief class to implement one single coordinate bin
 */
class coordinateSubBin
{
  public:
    /** @brief default constructor */
    coordinateSubBin() : 
      left( 0 ), 
      right( 0 ), 
      content( 0 ) 
    {};

    /** @brief destructor */
    ~coordinateSubBin() {};

    /** @brief set the values of #left and #right */
    void setLeftRight( const double _l, const double _r ) { left = _l; right = _r; }
  
    /** @brief increase #content, prefix */
    coordinateSubBin& operator++() { ++content; return *this; }

    /** @brief increase #content, postfix 
     * 
     * __very expensive, do not use!__
     */
    coordinateSubBin operator++( int unused ) 
    { 
      coordinateSubBin temp = *this; 
      ++content; 
      return temp; 
    }
  
    /** @brief reset all values to its defaults */
    void clear() { left = 0; right = 0; content = 0; }
    
    double left; ///< ...
    double right; ///< ...
    int content; ///< the number of entries

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & BOOST_SERIALIZATION_NVP( left );
      ar & BOOST_SERIALIZATION_NVP( right );
      ar & BOOST_SERIALIZATION_NVP( content );
    }
};

/**
 * @brief abstract class to implement coordinate bins
 *
 * The bins are implemented as a vector of #coordinateSubBin
 */
class coordinateBins
{
  public:

    /** @brief default constructor */
    coordinateBins() : 
      _min_real( 0 ), 
      _max_real( 0 ), 
      _delta_x( 0 ),
      _min_index_limit( 0 ), 
      _max_index_limit( 0 ),  
      _min_index_active( 0 ), 
      _max_index_active( 0 ),
      _negative_indices( false ) 
    { 
      bins.clear(); 
    }

    /** @brief construtor */
    coordinateBins( const int _size, const double _min, const double _max );

    /** @brief destructor */
    ~coordinateBins() {};

    /** @brief return the value of #_min_real */
    double min_real()  { return _min_real; }

    /** @brief return the value of #_max_real */
    double max_real()  { return _max_real; }

    /** @brief return the value of #_min_index_active */
    int min_index() { return _min_index_active; }

    /** @brief return the value of #_max_index_active */
    int max_index() { return _max_index_active; }
  
    /** @brief return the value of #_min_index_limit */
    int min_index_limit() { return _min_index_limit; }

    /** @brief return the value of #_max_index_limit */
    int max_index_limit() { return _max_index_limit; }

    /** @brief return the length of the bins vector, i.e. its size() */
    int size() { return bins.size(); }

    /** @brief return the value of #_delta_x */
    double get_dx() { return _delta_x; }

    /** @brief reset the object to its default values */
    void clear() { *this = coordinateBins( 0, 0, 0 ); }

    /** @brief ... */
    void reshape( const int _size, const double _min, const double _max );

    /** @brief ... */
    const coordinateSubBin& operator[]( const int _index ) const;

    /** @brief ... */
    int getIndex( const double _x ) const;
  
    /** @brief ... */
    void increase( const double _x ) { ++bins[ getIndex( _x ) ]; }
  
  

  protected:
    vector< coordinateSubBin > bins; ///< the vector of the bins

    double _min_real; ///< ...
    double _max_real; ///< ...
    double _delta_x; ///< ...

    int _min_index_limit; ///< ...
    int _max_index_limit; ///< ...
    
    int _min_index_active; ///< ...
    int _max_index_active; ///< ...

    bool _negative_indices; ///< ...
 
  private:
  
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
      ar & BOOST_SERIALIZATION_NVP( bins );
      ar & BOOST_SERIALIZATION_NVP( _min_real );
      ar & BOOST_SERIALIZATION_NVP( _max_real );
      ar & BOOST_SERIALIZATION_NVP( _delta_x );
      ar & BOOST_SERIALIZATION_NVP( _min_index_limit );
      ar & BOOST_SERIALIZATION_NVP( _max_index_limit );
      ar & BOOST_SERIALIZATION_NVP( _min_index_active );
      ar & BOOST_SERIALIZATION_NVP( _max_index_active );
      ar & BOOST_SERIALIZATION_NVP( _negative_indices );  
    } 
};


/**
 * @brief class to implement coordinate bins in the variable eta
 */
class coordinateEtaBins : public coordinateBins
{
  public:

    /** @brief default constructor */
    coordinateEtaBins() : 
      coordinateBins(), 
      NinEtaBin( 0 ), 
      timestepScaling(0.2) 
    {};

    /** @brief constructor */
    coordinateEtaBins( const int _size, 
                       const double _min, 
                       const double _max, 
                       const double _scaleTimestep = 0.2 ) : 
      coordinateBins( _size, _min, _max ), 
      NinEtaBin(0), 
      timestepScaling(_scaleTimestep) 
    {};

    /** @brief destructor */
    ~coordinateEtaBins() {};
  
    /** @brief ... */
    void populateEtaBins( std::vector<Particle> & parts, 
                          coordinateBins& _dNdEta, 
                          const double _etaShift, 
                          double _timenow, 
                          double& _dt, 
                          const double _dx, 
                          const double _dEta_fine
                          );
    
    void populateEtaBinsManually( std::vector<Particle> & parts, coordinateBins& _dNdEta, const double _etaShift, double _timenow, double& _dt, const double _dx, const double _dEta_fine, const double eta_min_fixed, const double eta_max_fixed );
    /** @brief ... */
    int constructEtaBins( const int _NperCell, 
                          const double _b, 
                          const double _dx, 
                          const double _dy, 
                          const double _RA, 
                          const int _nTest, 
                          const int particlesSize );

    int constructEtaBinsPPb( const int _NperCell, const double _transLen, const double _dx, const double _dy, const double _RA, const int _nTest, const int particlesSize );    
    int constructEtaBinsManually( const int _NperCell, const double _transLen, const double _dx, const double _dy, const double _RA, const int _nTest, const int particlesSize );    
        
    /** @brief ... */
    int getCentralIndex() const;

    /** @brief ... */
    int getIndex( const double _eta ) const;

    /** @brief set the value of #timestepScaling */
    void setTimestepScaling( const double _scaleTimestep ) { timestepScaling = _scaleTimestep; }

    /** @brief return the value of #NinEtaBin */
    int getNinEtaBin() const { return NinEtaBin; }

private:
  unsigned int NinEtaBin; ///< ...
  
  double timestepScaling; ///< ...
  
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( coordinateBins );
    ar & BOOST_SERIALIZATION_NVP( NinEtaBin );
  } 
};



/** 
 * @brief exception class for handling unexpected critical behaviour within coordinate bin errors 
 * 
 */
class eCoordBins_error : public std::runtime_error
{
  public:
    explicit eCoordBins_error(const std::string& what) : std::runtime_error(what) {};
    
    virtual ~eCoordBins_error() throw() {};
};


#endif
