//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/offlineoutput.h $
//$LastChangedDate: 2014-11-16 21:18:55 +0100 (日, 16 11月 2014) $
//$LastChangedRevision: 1935 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

// at revision 902, this file is identical to full/branches/vector4D/src/offlineoutput.h

#ifndef OFFLINE_OUTPUT_H
#define OFFLINE_OUTPUT_H

// #define BINARY_OFFLINE_OUTPUT 0
#define BINARY_OFFLINE_OUTPUT 1

// if BINARY_OFFLINE_OUTPUT==0, i.e. with text output,
// you may switch to XML format (for output only)
#define XML_OFFLINE_OUTPUT 0
// #define XML_OFFLINE_OUTPUT 1

#define DEFAULT_EVENT_NUMBER_ESTIMATE 1000000
#define DEFAULT_TEMPORARY_EVENT_NUMBER_ESTIMATE 100

#define BAMPS_RECONSTRUCT_ACTIVATE_READOUT

#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <typeindex>
#include <vector>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/level.hpp>
#include <boost/smart_ptr.hpp>

#ifdef BAMPS_RECONSTRUCT_ACTIVATE_READOUT
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#endif 

// #define BOOST_FILESYSTEM_VERSION 3
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "coordinateBins.h"
#include "particle.h"
#include "bampsvector.h"

class eOfflineOutput_error;
enum offlineEventType { event_interaction22, event_interaction23, event_interaction32, event_interactionElastic, 
                        event_particleIdSwap, event_newTimestep, event_endOfCascade, event_dummy = 99 };


/**
 * @brief Generic Class for Offline Data
 */

class offlineDataGeneric
{
public:
  offlineDataGeneric() {};
  ~offlineDataGeneric() {};
    
  virtual std::string getFilenameIdentifier() const { return this->filenameIdentifier; }
  virtual size_t getSize() const = 0;
    
  static std::string filenameIdentifier;
    
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
  }
};


typedef boost::shared_ptr<offlineDataGeneric> tPointerToOfflineData;
typedef boost::shared_ptr<const offlineDataGeneric> tConstPointerToOfflineData;


/**
 * @brief Class... 
 */

class offlineOutputInterface
{
public:
  offlineOutputInterface( const char* _outputDir, const bool _active = true ) : 
    outputDirectory( _outputDir ), 
    filenamePrefix( "offline_" ), 
    filenameSuffix( ".dat" ),
    additionalFilenameTag( "" )
  { 
    if ( _active )
    {
      checkAndCreateOutputDirectory( outputDirectory );
      eventList.reserve( DEFAULT_EVENT_NUMBER_ESTIMATE ); 
      temporaryStorage.reserve( DEFAULT_TEMPORARY_EVENT_NUMBER_ESTIMATE );
    }
  }
  offlineOutputInterface( const char* _outputDir, const int eventNumberEstimate, const bool _active = true ) : 
    outputDirectory( _outputDir ), 
    filenamePrefix( "offline_" ), 
    filenameSuffix( ".dat" ),
    additionalFilenameTag( "" )
  { 
    if ( _active )
    {
      checkAndCreateOutputDirectory( outputDirectory );
      eventList.reserve( eventNumberEstimate );
      temporaryStorage.reserve( DEFAULT_TEMPORARY_EVENT_NUMBER_ESTIMATE );
    }
  }
  ~offlineOutputInterface();
    
  void setAdditionalFilenameTag( std::string _tag ) { additionalFilenameTag = _tag + "_"; }
    
  void submitOfflineDataForOutput( const offlineDataGeneric* const _data );
  void submitOfflineDataForOutput( const offlineEventType _eventType )
  { 
    eventList.push_back( _eventType );
  } 
    
  void registerOfflineDataForTemporaryStorage( tConstPointerToOfflineData _data )
  {
    temporaryStorage.push_back( _data );
  }
  void registerOfflineDataForTemporaryStorage( const offlineEventType _eventType, tConstPointerToOfflineData _data );
  void resetTemporaryStorage() { temporaryStorage.clear(); }
  void outputAndResetTemporaryStorage();
    
#ifdef BAMPS_RECONSTRUCT_ACTIVATE_READOUT
  /** read the stuff */
  template<class T>
  boost::shared_ptr<T> readOfflineDataFromArchive();
    
  /** provide ability to undo read operation */
  template<class T>
  void undoLastReadOperation();
    
  template<class T>
  void temporaryStoreData( tPointerToOfflineData _data );
#endif
    
private:
  void checkAndCreateOutputDirectory( boost::filesystem::path& _dir );
    
  typedef boost::shared_ptr<boost::filesystem::ofstream> tPointerToFilestream;
  typedef std::map< std::type_index, tPointerToFilestream > tFileStreamMap;
  tFileStreamMap outputMap;

#if BINARY_OFFLINE_OUTPUT > 0 
  typedef boost::archive::binary_oarchive tArchive;
#else
#if XML_OFFLINE_OUTPUT > 0
  typedef boost::archive::xml_oarchive tArchive;
#else
  typedef boost::archive::text_oarchive tArchive;
#endif
#endif
  typedef boost::shared_ptr<tArchive> tPointerToArchive;
  typedef std::map< std::type_index, tPointerToArchive > tArchiveMap;
  tArchiveMap archiveMap;
    
    
#ifdef BAMPS_RECONSTRUCT_ACTIVATE_READOUT
  /** read the stuff */
  typedef boost::shared_ptr<boost::filesystem::ifstream> tPointerToInputFilestream;
  typedef std::map< std::type_index, tPointerToInputFilestream > tInputFileStreamMap;
  tInputFileStreamMap inputMap;

#if BINARY_OFFLINE_OUTPUT > 0 
  typedef boost::archive::binary_iarchive tInputArchive;
#else
  typedef boost::archive::text_iarchive tInputArchive;
#endif
  typedef boost::shared_ptr<tInputArchive> tPointerToInputArchive;
  typedef std::map< std::type_index, tPointerToInputArchive > tInputArchiveMap;
  tInputArchiveMap inputArchiveMap;
  /** read the stuff */
    
  /** save the last read position */
  typedef std::map< std::type_index, std::streampos > tStreamPositionMap;
  tStreamPositionMap lastStreamPositionMap;
    
  typedef std::map< std::type_index, tPointerToOfflineData > tOfflineDataMap;
  tOfflineDataMap backupLastReadData;     
  /** save the last read position */
#endif
    
  typedef std::vector<tConstPointerToOfflineData> tTemporaryDataStorage;
  tTemporaryDataStorage temporaryStorage;
    
  std::vector<offlineEventType> eventList;
  void collectAndOutputEventList();    
    
  boost::filesystem::path outputDirectory;
  std::string filenamePrefix;
  std::string filenameSuffix;
  std::string additionalFilenameTag;
};


inline void offlineOutputInterface::checkAndCreateOutputDirectory(boost::filesystem::path& _dir)
{
  if ( boost::filesystem::exists( _dir ) )
  {
    if ( boost::filesystem::is_directory( _dir ) )
    {
      return;
    }
    else
    {
      boost::filesystem::path renamePath( _dir.string() + ".backup" );
      std::cout << "File with name " << _dir.string() << " blocks the creation of an output folder for offline reconstruction." << std::endl;
      std::cout << "It is renamed to " << renamePath.string() << std::endl;
      boost::filesystem::rename( _dir, renamePath );
      boost::filesystem::create_directory( _dir );       
    }
  }
  else
  {
    std::cout << "Creating output folder for offline reconstruction data: " << _dir.string() << std::endl;
    boost::filesystem::create_directory( _dir );    
  }
}

#ifdef BAMPS_RECONSTRUCT_ACTIVATE_READOUT
template<class T>
boost::shared_ptr<T> offlineOutputInterface::readOfflineDataFromArchive()
{
  tPointerToInputFilestream stream;
  tPointerToInputArchive archive;
  
  //   check whether an output archive for the given type of data has already been created, re-use if yes, create if no
  tInputArchiveMap::iterator it = inputArchiveMap.find( std::type_index(typeid(T)) );
  if( it != inputArchiveMap.end() )
  {
    tInputFileStreamMap::iterator itFile = inputMap.find( std::type_index(typeid(T)) );
    if( itFile != inputMap.end() )
    {
      stream = itFile->second;    
    }
    archive = it->second;
  }
  else
  {
    stream.reset( new boost::filesystem::ifstream() );
    std::string filename = outputDirectory.string() + "/" + filenamePrefix + additionalFilenameTag + T::filenameIdentifier + filenameSuffix;
    
    if( !boost::filesystem::exists( filename ) )
    {
      std::string errMsg = filename + " does not exist.";
      throw eOfflineOutput_error( errMsg );
    }

#if BINARY_OFFLINE_OUTPUT > 0
    stream->open( filename.c_str(), std::ios::binary );
#else
    stream->open( filename.c_str() );
#endif
    
    inputMap.insert(std::make_pair( std::type_index(typeid(T)), stream));
    archive.reset( new tInputArchive( *stream ) );
    inputArchiveMap.insert(std::make_pair( std::type_index(typeid(T)), archive));
    
    lastStreamPositionMap.insert( std::make_pair( std::type_index(typeid(T)), 0 ));
  }
  
  tOfflineDataMap::iterator itBackupData = backupLastReadData.find( std::type_index(typeid(T)) );
  if( itBackupData != backupLastReadData.end() )
  {
    boost::shared_ptr< T > shrptr = boost::static_pointer_cast< T >( itBackupData->second );
    backupLastReadData.erase( itBackupData );
    return shrptr;
  }
  else
  {
    // the actual reading
    tStreamPositionMap::iterator itLastPos = lastStreamPositionMap.find( std::type_index(typeid(T)) );
    if( itLastPos != lastStreamPositionMap.end() )
    {
      itLastPos->second = stream->tellg();    
    }
    offlineDataGeneric* _ptr = 0;
    (*archive) & _ptr;
    
    T* ptrToDerived = dynamic_cast<T*>(_ptr);
    if ( ptrToDerived == 0 )
    {
      std::string errMsg = "Bad cast. Attempted dynamic_cast from ";
      errMsg += std::type_index(typeid(offlineDataGeneric)).name();
      errMsg += " to ";
      errMsg += std::type_index(typeid(T)).name();
      throw eOfflineOutput_error( errMsg );
    }
    
    boost::shared_ptr<T> shrptr( ptrToDerived );
    return shrptr;
  }
}


template<class T>
void offlineOutputInterface::undoLastReadOperation()
{
  tPointerToInputFilestream stream;
  
  std::streampos lastStreamPosition;
  tStreamPositionMap::iterator itLastPos = lastStreamPositionMap.find( std::type_index(typeid(T)) );
  if( itLastPos != lastStreamPositionMap.end() )
  {
    lastStreamPosition = itLastPos->second;    
  }
  else
  {
    std::string errMsg = "Attempting undo operation on stream that has not been initialized yet";
    throw eOfflineOutput_error( errMsg );
  }
  
  //   check whether an output archive for the given type of data has already been created, re-use if yes, create if no
  tInputFileStreamMap::iterator it = inputMap.find( std::type_index(typeid(T)) );
  if( it != inputMap.end() )
  {
    stream = it->second;
    stream->seekg( lastStreamPosition );
  }
}


template<class T>
void offlineOutputInterface::temporaryStoreData(tPointerToOfflineData _data)
{
  tOfflineDataMap::iterator itBackupData = backupLastReadData.find( std::type_index(typeid(T)) );
  if( itBackupData != backupLastReadData.end() )
  {
    std::string errMsg = "Attempting to store temporary data into a storage that is already in use";
    throw eOfflineOutput_error( errMsg );
  }
  else
  {    
    backupLastReadData.insert(std::make_pair( std::type_index(typeid(T)), _data));
  }
}
#endif


/**
 * @brief Class... 
 */

class offlineDataEventType : public offlineDataGeneric
{
public:
  offlineDataEventType() : offlineDataGeneric(), event( event_interaction22 ) {};
  offlineDataEventType( const offlineEventType _event ) : offlineDataGeneric(), event(_event) {};
  ~offlineDataEventType() {};
    
  size_t getSize() const { return sizeof( offlineDataEventType ); }
  std::string getFilenameIdentifier() const { return this->filenameIdentifier; }
  static std::string filenameIdentifier;
    
  /** @brief event type */
  offlineEventType event;
    
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( offlineDataGeneric );
    ar & BOOST_SERIALIZATION_NVP( event );
  }
};

/**
 * @brief Class for offline output of 22 scattering
 */

class offlineDataInteraction22 : public offlineDataGeneric
{
public:
  offlineDataInteraction22() : 
    offlineDataGeneric(),
    iscat(-1), jscat(-1), time(-1), 
    Mom1(-999,-999,-999,-999),
    Mom2(-999,-999,-999,-999),
    F1(gluon), F2(gluon) {};

  offlineDataInteraction22( const int _iscat, const int _jscat, const double _time, 
                            const VectorEPxPyPz & _Mom1,
                            const VectorEPxPyPz & _Mom2,
                            const FLAVOR_TYPE _F1, const FLAVOR_TYPE _F2 ) : 
    offlineDataGeneric(),
    iscat(_iscat), jscat(_jscat), time(_time),
    Mom1( _Mom1 ),
    Mom2( _Mom2 ),
    F1(_F1), F2(_F2) {};

  ~offlineDataInteraction22() {};
    
  size_t getSize() const { return sizeof( offlineDataInteraction22 ); }
  std::string getFilenameIdentifier() const { return this->filenameIdentifier; }
  static std::string filenameIdentifier;
    
  /** @brief Current ID of particle 1 */
  int iscat;
  /** @brief Current ID of particle 2 */
  int jscat;
  /** @brief Current time of the collision */
  double time;

  /** @brief Momentum of particle 1 AFTER the collison */
  VectorEPxPyPz Mom1;
  /** @brief Momentum of particle 2 AFTER the collison */
  VectorEPxPyPz Mom2;

  /** @brief Flavor of particle 1 AFTER the collison */
  FLAVOR_TYPE F1; 
  /** @brief Flavor of particle 1 AFTER the collison */
  FLAVOR_TYPE F2;
    
  friend class boost::serialization::access;

  template<class Archive>
  void load(Archive & ar, const unsigned int version)
  {
    double x, y, z;

    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( offlineDataGeneric );
    ar & BOOST_SERIALIZATION_NVP( iscat );
    ar & BOOST_SERIALIZATION_NVP( jscat );
    ar & BOOST_SERIALIZATION_NVP( time );

    ar & boost::serialization::make_nvp( "pix", x );
    ar & boost::serialization::make_nvp( "piy", y );
    ar & boost::serialization::make_nvp( "piz", z );
    Mom1.SetTXYZ(0.0,x,y,z);

    ar & boost::serialization::make_nvp( "pjx", x );
    ar & boost::serialization::make_nvp( "pjy", y );
    ar & boost::serialization::make_nvp( "pjz", z );
    Mom2.SetTXYZ(0.0,x,y,z);

    ar & BOOST_SERIALIZATION_NVP( F1 );
    ar & BOOST_SERIALIZATION_NVP( F2 );
  }

  template<class Archive>
  void save(Archive & ar, const unsigned int version) const
  {
    double t, x, y, z;

    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( offlineDataGeneric );
    ar & BOOST_SERIALIZATION_NVP( iscat );
    ar & BOOST_SERIALIZATION_NVP( jscat );
    ar & BOOST_SERIALIZATION_NVP( time );

    Mom1.GetTXYZ(t,x,y,z);
    ar & boost::serialization::make_nvp( "pix", x );
    ar & boost::serialization::make_nvp( "piy", y );
    ar & boost::serialization::make_nvp( "piz", z );

    Mom2.GetTXYZ(t,x,y,z);
    ar & boost::serialization::make_nvp( "pjx", x );
    ar & boost::serialization::make_nvp( "pjy", y );
    ar & boost::serialization::make_nvp( "pjz", z );

    ar & BOOST_SERIALIZATION_NVP( F1 );
    ar & BOOST_SERIALIZATION_NVP( F2 );
  }
  BOOST_SERIALIZATION_SPLIT_MEMBER()
};

/**
 * @brief Class for offline output of 23 scattering
 */

class offlineDataInteraction23 : public offlineDataGeneric
{
public:
  offlineDataInteraction23() : 
    offlineDataGeneric(),
    iscat(-1), jscat(-1), newp(-1), time(-1), 
    Mom1(-999,-999,-999,-999),
    Mom2(-999,-999,-999,-999),
    Pos3(-999,-999,-999,-999),
    Mom3(-999,-999,-999,-999),
    F1(gluon), F2(gluon), F3(gluon) {};

  offlineDataInteraction23( const int _iscat, const int _jscat, const int _newp, 
                            const double _time, 
                            const VectorEPxPyPz & _Mom1,
                            const VectorEPxPyPz & _Mom2,
                            const VectorTXYZ & _Pos3,
                            const VectorEPxPyPz & _Mom3,
                            const FLAVOR_TYPE _F1, const FLAVOR_TYPE _F2, const FLAVOR_TYPE _F3 ) : 
    offlineDataGeneric(),
    iscat(_iscat), jscat(_jscat), newp(_newp), time(_time), 
    Mom1( _Mom1 ),
    Mom2( _Mom2 ),
    Pos3( _Pos3 ),
    Mom3( _Mom3 ),
    F1(_F1), F2(_F2), F3(_F3) {};

  ~offlineDataInteraction23() {};
    
  size_t getSize() const { return sizeof( offlineDataInteraction23 ); }
  std::string getFilenameIdentifier() const { return this->filenameIdentifier; }
  static std::string filenameIdentifier;
    
  /** @brief Current ID of particle 1 */
  int iscat;
  /** @brief Current ID of particle 2 */
  int jscat;
  /** @brief ID of newly created particle */
  int newp;
  /** @brief Current time of the collision */
  double time;

  /** @brief Momentum of particle 1 AFTER the collison */
  VectorEPxPyPz Mom1;
  /** @brief Momentum of particle 2 AFTER the collison */
  VectorEPxPyPz Mom2;

  /** @brief Position of new particle (particle 3) */
  VectorTXYZ Pos3;
  /** @brief Momentum of new particle (particle 3) */
  VectorEPxPyPz Mom3;

  /** @brief Flavor of particle 1 AFTER the collison */
  FLAVOR_TYPE F1; 
  /** @brief Flavor of particle 1 AFTER the collison */
  FLAVOR_TYPE F2;
  /** @brief Flavor of particle 3 AFTER the collison */
  FLAVOR_TYPE F3;
    
  friend class boost::serialization::access;

  template<class Archive>
  void load(Archive & ar, const unsigned int version)
  {
    double x, y, z;

    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( offlineDataGeneric );
    ar & BOOST_SERIALIZATION_NVP( iscat );
    ar & BOOST_SERIALIZATION_NVP( jscat );
    ar & BOOST_SERIALIZATION_NVP( newp );
    ar & BOOST_SERIALIZATION_NVP( time );

    ar & boost::serialization::make_nvp( "pix", x );
    ar & boost::serialization::make_nvp( "piy", y );
    ar & boost::serialization::make_nvp( "piz", z );
    Mom1.SetTXYZ(0.0,x,y,z);

    ar & boost::serialization::make_nvp( "pjx", x );
    ar & boost::serialization::make_nvp( "pjy", y );
    ar & boost::serialization::make_nvp( "pjz", z );
    Mom2.SetTXYZ(0.0,x,y,z);

    ar & boost::serialization::make_nvp( "newx", x );
    ar & boost::serialization::make_nvp( "newy", y );
    ar & boost::serialization::make_nvp( "newz", z );
    Pos3.SetTXYZ(0.0,x,y,z);

    ar & boost::serialization::make_nvp( "newpx", x );
    ar & boost::serialization::make_nvp( "newpy", y );
    ar & boost::serialization::make_nvp( "newpz", z );
    Mom3.SetTXYZ(0.0,x,y,z);

    ar & BOOST_SERIALIZATION_NVP( F1 );
    ar & BOOST_SERIALIZATION_NVP( F2 );
    ar & BOOST_SERIALIZATION_NVP( F3 );
  }

  template<class Archive>
  void save(Archive & ar, const unsigned int version) const
  {
    double t, x, y, z;

    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( offlineDataGeneric );
    ar & BOOST_SERIALIZATION_NVP( iscat );
    ar & BOOST_SERIALIZATION_NVP( jscat );
    ar & BOOST_SERIALIZATION_NVP( newp );
    ar & BOOST_SERIALIZATION_NVP( time );

    Mom1.GetTXYZ(t,x,y,z);
    ar & boost::serialization::make_nvp( "pix", x );
    ar & boost::serialization::make_nvp( "piy", y );
    ar & boost::serialization::make_nvp( "piz", z );

    Mom2.GetTXYZ(t,x,y,z);
    ar & boost::serialization::make_nvp( "pjx", x );
    ar & boost::serialization::make_nvp( "pjy", y );
    ar & boost::serialization::make_nvp( "pjz", z );

    Pos3.GetTXYZ(t,x,y,z);
    ar & boost::serialization::make_nvp( "newx", x );
    ar & boost::serialization::make_nvp( "newy", y );
    ar & boost::serialization::make_nvp( "newz", z );

    Mom3.GetTXYZ(t,x,y,z);
    ar & boost::serialization::make_nvp( "newpx", x );
    ar & boost::serialization::make_nvp( "newpy", y );
    ar & boost::serialization::make_nvp( "newpz", z );

    ar & BOOST_SERIALIZATION_NVP( F1 );
    ar & BOOST_SERIALIZATION_NVP( F2 );
    ar & BOOST_SERIALIZATION_NVP( F3 );
  }
  BOOST_SERIALIZATION_SPLIT_MEMBER()
};

/**
 * @brief Class for offline output of 32 scattering
 */

class offlineDataInteraction32 : public offlineDataGeneric
{
public:
  offlineDataInteraction32() : 
    offlineDataGeneric(),
    iscat(-1), jscat(-1), dead(-1), time(-1), 
    Mom1(-999,-999,-999,-999),
    Mom2(-999,-999,-999,-999),
    F1(gluon), F2(gluon) {};

  offlineDataInteraction32( const int _iscat, const int _jscat, const int _dead, 
                            const double _time, 
                            const VectorEPxPyPz & _Mom1,
                            const VectorEPxPyPz & _Mom2,
                            const FLAVOR_TYPE _F1, const FLAVOR_TYPE _F2 ) : 
    offlineDataGeneric(),
    iscat(_iscat), jscat(_jscat), dead(_dead), 
    time(_time), 
    Mom1( _Mom1 ),
    Mom2( _Mom2 ),
    F1(_F1), F2(_F2) {};

  ~offlineDataInteraction32() {};
    
  size_t getSize() const { return sizeof( offlineDataInteraction32 ); }
  std::string getFilenameIdentifier() const { return this->filenameIdentifier; }
  static std::string filenameIdentifier;
    
  /** @brief Current ID of particle 1 */
  int iscat;
  /** @brief Current ID of particle 2 */
  int jscat;
  /** @brief ID of absorbed ("dead") particle */
  int dead;
  /** @brief Current time of the collision */
  double time;

  /** @brief Momentum of particle 1 AFTER the collison */
  VectorEPxPyPz Mom1;
  /** @brief Momentum of particle 2 AFTER the collison */
  VectorEPxPyPz Mom2;

  /** @brief Flavor of particle 1 AFTER the collison */
  FLAVOR_TYPE F1; 
  /** @brief Flavor of particle 1 AFTER the collison */
  FLAVOR_TYPE F2;
    
  friend class boost::serialization::access;

  template<class Archive>
  void load(Archive & ar, const unsigned int version)
  {
    double x, y, z;

    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( offlineDataGeneric );
    ar & BOOST_SERIALIZATION_NVP( iscat );
    ar & BOOST_SERIALIZATION_NVP( jscat );
    ar & BOOST_SERIALIZATION_NVP( dead );
    ar & BOOST_SERIALIZATION_NVP( time );

    ar & boost::serialization::make_nvp( "pix", x );
    ar & boost::serialization::make_nvp( "piy", y );
    ar & boost::serialization::make_nvp( "piz", z );
    Mom1.SetTXYZ(0.0,x,y,z);

    ar & boost::serialization::make_nvp( "pjx", x );
    ar & boost::serialization::make_nvp( "pjy", y );
    ar & boost::serialization::make_nvp( "pjz", z );
    Mom2.SetTXYZ(0.0,x,y,z);

    ar & BOOST_SERIALIZATION_NVP( F1 );
    ar & BOOST_SERIALIZATION_NVP( F2 );
  }

  template<class Archive>
  void save(Archive & ar, const unsigned int version) const
  {
    double t, x, y, z;

    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( offlineDataGeneric );
    ar & BOOST_SERIALIZATION_NVP( iscat );
    ar & BOOST_SERIALIZATION_NVP( jscat );
    ar & BOOST_SERIALIZATION_NVP( dead );
    ar & BOOST_SERIALIZATION_NVP( time );

    Mom1.GetTXYZ(t,x,y,z);
    ar & boost::serialization::make_nvp( "pix", x );
    ar & boost::serialization::make_nvp( "piy", y );
    ar & boost::serialization::make_nvp( "piz", z );

    Mom2.GetTXYZ(t,x,y,z);
    ar & boost::serialization::make_nvp( "pjx", x );
    ar & boost::serialization::make_nvp( "pjy", y );
    ar & boost::serialization::make_nvp( "pjz", z );

    ar & BOOST_SERIALIZATION_NVP( F1 );
    ar & BOOST_SERIALIZATION_NVP( F2 );
  }
  BOOST_SERIALIZATION_SPLIT_MEMBER()
};

/**
 * @brief Class for offline output of elastic scattering
 */

class offlineDataInteractionElastic : public offlineDataGeneric
{
public:
  offlineDataInteractionElastic( ) : 
    offlineDataGeneric(),
    iscat(-1), jscat(-1), 
    ct_i(-1), ct_j(-1), 
    Mom1(-999,-999,-999,-999),
    Mom2(-999,-999,-999,-999) {};

  offlineDataInteractionElastic( const int _iscat, const int _jscat, 
                                 const double _ct_i, const double _ct_j, 
                                 const VectorEPxPyPz & _Mom1,
                                 const VectorEPxPyPz & _Mom2) :
    offlineDataGeneric(),
    iscat(_iscat), jscat(_jscat), 
    ct_i(_ct_i), ct_j(_ct_j), 
    Mom1( _Mom1 ),
    Mom2( _Mom2 ) {};

  ~offlineDataInteractionElastic() {};
    
  size_t getSize() const { return sizeof( offlineDataInteractionElastic ); }
  std::string getFilenameIdentifier() const { return this->filenameIdentifier; }
  static std::string filenameIdentifier;
    
  /** @brief Current ID of particle 1 */
  int iscat;
  /** @brief Current ID of particle 2 */
  int jscat;
  /** @brief Geometric collision time of particle 1 */
  double ct_i;
  /** @brief Geometric collision time of particle 2 */
  double ct_j;

  /** @brief Momentum of particle 1 AFTER the collison */
  VectorEPxPyPz Mom1;
  /** @brief Momentum of particle 2 AFTER the collison */
  VectorEPxPyPz Mom2;

  friend class boost::serialization::access;

  template<class Archive>
  void load(Archive & ar, const unsigned int version)
  {
    double x, y, z;

    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( offlineDataGeneric );
    ar & BOOST_SERIALIZATION_NVP( iscat );
    ar & BOOST_SERIALIZATION_NVP( jscat );
    ar & BOOST_SERIALIZATION_NVP( ct_i );
    ar & BOOST_SERIALIZATION_NVP( ct_j );
    
    ar & boost::serialization::make_nvp( "pix", x );
    ar & boost::serialization::make_nvp( "piy", y );
    ar & boost::serialization::make_nvp( "piz", z );
    Mom1.SetTXYZ(0.0,x,y,z);

    ar & boost::serialization::make_nvp( "pjx", x );
    ar & boost::serialization::make_nvp( "pjy", y );
    ar & boost::serialization::make_nvp( "pjz", z );
    Mom2.SetTXYZ(0.0,x,y,z);
  }

  template<class Archive>
  void save(Archive & ar, const unsigned int version) const
  {
    double t, x, y, z;

    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( offlineDataGeneric );
    ar & BOOST_SERIALIZATION_NVP( iscat );
    ar & BOOST_SERIALIZATION_NVP( jscat );
    ar & BOOST_SERIALIZATION_NVP( ct_i );
    ar & BOOST_SERIALIZATION_NVP( ct_j );
    
    Mom1.GetTXYZ(t,x,y,z);
    ar & boost::serialization::make_nvp( "pix", x );
    ar & boost::serialization::make_nvp( "piy", y );
    ar & boost::serialization::make_nvp( "piz", z );

    Mom2.GetTXYZ(t,x,y,z);
    ar & boost::serialization::make_nvp( "pjx", x );
    ar & boost::serialization::make_nvp( "pjy", y );
    ar & boost::serialization::make_nvp( "pjz", z );
  }
  BOOST_SERIALIZATION_SPLIT_MEMBER()
};

/**
 * @brief Class... 
 */

class offlineDataCellConfiguration : public offlineDataGeneric
{
public:
  offlineDataCellConfiguration( ) : 
    offlineDataGeneric(), 
    timenow(-1), timenext(-1),
    randomShiftX(0), randomShiftY(0), randomShiftEta(0) {};

  offlineDataCellConfiguration( const double _timenow, const double _timenext, 
                                const double _randomshiftX, 
                                const double _randomShiftY, 
                                const double _randomShiftEta,
                                const coordinateEtaBins& _etaBins ) :
    offlineDataGeneric(),
    timenow(_timenow), timenext(_timenext),
    randomShiftX(_randomshiftX), randomShiftY(_randomShiftY), 
    randomShiftEta(_randomShiftEta),
    etaBins(_etaBins) {};
  ~offlineDataCellConfiguration() {};
    
  size_t getSize() const { return sizeof( offlineDataCellConfiguration ); }
  std::string getFilenameIdentifier() const { return this->filenameIdentifier; }
  static std::string filenameIdentifier;
    
  /** @brief Current simulation time */
  double timenow;
  /** @brief The next time step */
  double timenext;
  /** @brief Random shift of the cell grid in x-direction */
  double randomShiftX;
  /** @brief Random shift of the cell grid in y-direction */
  double randomShiftY;
  /** @brief Random shift of the cell grid in eta-direction (= z-direction) */
  double randomShiftEta;
  /** @brief The current cell configuration in eta-direction */
  coordinateEtaBins etaBins;    
    
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( offlineDataGeneric );
    ar & BOOST_SERIALIZATION_NVP( timenow );
    ar & BOOST_SERIALIZATION_NVP( timenext );
    ar & BOOST_SERIALIZATION_NVP( randomShiftX );
    ar & BOOST_SERIALIZATION_NVP( randomShiftY );
    ar & BOOST_SERIALIZATION_NVP( randomShiftEta );
    ar & BOOST_SERIALIZATION_NVP( etaBins );
  }
};

/**
 * @brief Class... 
 */

class offlineDataParticleIdSwap : public offlineDataGeneric
{
public:
  offlineDataParticleIdSwap() : 
    offlineDataGeneric(),
    removedParticleID(-1), replacingParticleID(-1) {};
  offlineDataParticleIdSwap( const int _removedParticleID, const int _replacingParticleID ) : 
    offlineDataGeneric(),
    removedParticleID(_removedParticleID),
    replacingParticleID(_replacingParticleID) {};
    
  ~offlineDataParticleIdSwap() {};
    
  size_t getSize() const { return sizeof( offlineDataParticleIdSwap ); }
  std::string getFilenameIdentifier() const { return this->filenameIdentifier; }
  static std::string filenameIdentifier;
    
  /** @brief ID of the removed ("dead") particle */
  int removedParticleID;
  /** @brief ID of the particle that takes the place of the removed particle */
  int replacingParticleID;
    
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( offlineDataGeneric );
    ar & BOOST_SERIALIZATION_NVP( removedParticleID );
    ar & BOOST_SERIALIZATION_NVP( replacingParticleID );
  }
};

/**
 * @brief Class... 
 */

class offlineDataParticleNumbers : public offlineDataGeneric
{
public:
  offlineDataParticleNumbers( ) : 
    offlineDataGeneric(),
    time(-1), dt(-1), nTotalParticles(0),
    nParticlesInFormationIntialProduction(0),
    nParticlesInFormationGeometricCollisions(0),
    nParticlesActive(0), nParticlesInActiveCells(0),
    edgeCellSizes(0), nFreeParticles( 0 ) {};
  offlineDataParticleNumbers( const double _time, const double _dt, const int _nTotalParticles,
                              const int _nParticlesInFormationIntialProduction, 
                              const int _nParticlesInFormationGeometricCollisions, const int _nParticlesActive,
                              const int _nParticlesInActiveCells, const std::vector<int>& _edgeCellSizes,
                              const int _nFreeParticles ) : 
    offlineDataGeneric(),
    time(_time), dt(_dt), nTotalParticles(_nTotalParticles),
    nParticlesInFormationIntialProduction(_nParticlesInFormationIntialProduction),
    nParticlesInFormationGeometricCollisions(_nParticlesInFormationGeometricCollisions),
    nParticlesActive(_nParticlesActive), nParticlesInActiveCells(_nParticlesInActiveCells),
    edgeCellSizes(_edgeCellSizes), nFreeParticles( _nFreeParticles ) {};
  ~offlineDataParticleNumbers() {};
    
  size_t getSize() const { return sizeof( offlineDataParticleNumbers ); }
  std::string getFilenameIdentifier() const { return this->filenameIdentifier; }
  static std::string filenameIdentifier;
    
  /** @brief Current simulation time */
  double time;
  /** @brief Current time step \Delta t */
  double dt;
  /** @brief Current size of the particle vector */
  int nTotalParticles;
  //    /** @brief Number of particles that are in formation */
  //     int nParticlesInFormation;
  /** @brief Number of particles that are in formation from initial production */
  int nParticlesInFormationIntialProduction;
  /** @brief Number of particles that are in formation due to geometric collisions */
  int nParticlesInFormationGeometricCollisions;
  /** @brief The number of active particles, i.e. particles for which T < timenext */
  int nParticlesActive;
  /** @brief Number of particles in "active" cells, i.e. in cells whose particle content is larger than ::cellcut */
  int nParticlesInActiveCells;
  /** @brief The sizes of the edge cells (= number of particles in these cells) */
  std::vector<int> edgeCellSizes;
  /** @brief Total number of particles that are free because they are a) outside the transversal grid or b) in free cells */
  int nFreeParticles;
    
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( offlineDataGeneric );
    ar & BOOST_SERIALIZATION_NVP( time );
    ar & BOOST_SERIALIZATION_NVP( dt );
    ar & BOOST_SERIALIZATION_NVP( nTotalParticles );
    //       ar & BOOST_SERIALIZATION_NVP( nParticlesInFormation );
    ar & BOOST_SERIALIZATION_NVP( nParticlesInFormationIntialProduction );
    ar & BOOST_SERIALIZATION_NVP( nParticlesInFormationGeometricCollisions );
    ar & BOOST_SERIALIZATION_NVP( nParticlesActive );
    ar & BOOST_SERIALIZATION_NVP( nParticlesInActiveCells );
    ar & BOOST_SERIALIZATION_NVP( edgeCellSizes );
    ar & BOOST_SERIALIZATION_NVP( nFreeParticles );
  }
};

/**
 * @brief Class... 
 */

class offlineDataCollisionNumbers : public offlineDataGeneric
{
public:
  offlineDataCollisionNumbers( ) : 
    offlineDataGeneric(), time( -1 ),
    nCollisions22( 0 ), nCollisions23( 0 ), nCollisions32( 0 ),
    nCollisionsElastic( 0 ) {};
  offlineDataCollisionNumbers( const double _time, const int _nCollisions22, const int _nCollisions23, 
                               const int _nCollisions32, const int _nCollisionsElastic ) : 
    offlineDataGeneric(), time( _time ),
    nCollisions22( _nCollisions22 ), nCollisions23( _nCollisions23 ), nCollisions32( _nCollisions32 ),
    nCollisionsElastic( _nCollisionsElastic ) {};
  ~offlineDataCollisionNumbers() {};
    
  size_t getSize() const { return sizeof( offlineDataCollisionNumbers ); }
  std::string getFilenameIdentifier() const { return this->filenameIdentifier; }
  static std::string filenameIdentifier;
    
  /** @brief Current simulation time */
  double time;
  /** @brief Number of 2->2 interactions (divided by the test particle number) */
  int nCollisions22;
  /** @brief Number of 2->3 interactions (divided by the test particle number) */
  int nCollisions23;
  /** @brief Number of 3->2 interactions (divided by the test particle number) */
  int nCollisions32;
  /** @brief Number of geometric collisions (divided by the test particle number) */
  int nCollisionsElastic;
    
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( offlineDataGeneric );
    ar & BOOST_SERIALIZATION_NVP( time );
    ar & BOOST_SERIALIZATION_NVP( nCollisions22 );
    ar & BOOST_SERIALIZATION_NVP( nCollisions23 );
    ar & BOOST_SERIALIZATION_NVP( nCollisions32 );
    ar & BOOST_SERIALIZATION_NVP( nCollisionsElastic );
  }
};

/**
 * @brief Class... 
 */

class offlineDataInteractionRates : public offlineDataGeneric
{
public:
  offlineDataInteractionRates( ) : offlineDataGeneric() {};
  offlineDataInteractionRates( std::vector< vector<double> >& _gluonRates,
                               std::vector< vector<double> >& _quarkRates,
                               std::vector< vector<double> >& _antiQuarkRates ) : 
    offlineDataGeneric(),
    gluonRates(_gluonRates),
    quarkRates(_quarkRates), antiQuarkRates(_antiQuarkRates) {};
  ~offlineDataInteractionRates() {};

  size_t getSize() const { return sizeof( offlineDataInteractionRates ); }
  std::string getFilenameIdentifier() const { return this->filenameIdentifier; }
  static std::string filenameIdentifier;
                                 
  /** @brief Stores the rates of gluons for output needed for offline reconstruction */
  std::vector< std::vector<double> > gluonRates;
  /** @brief Stores the rates of light quarks for output needed for offline reconstruction */
  std::vector< std::vector<double> > quarkRates;
  /** @brief Stores the rates of light anti-quarks for output needed for offline reconstruction */
  std::vector< std::vector<double> > antiQuarkRates;
    
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( offlineDataGeneric );
    ar & BOOST_SERIALIZATION_NVP( gluonRates );
    ar & BOOST_SERIALIZATION_NVP( quarkRates );
    ar & BOOST_SERIALIZATION_NVP( antiQuarkRates );
  }
};

/**
 * @brief Class... 
 */

class offlineDataSimulationParameters : public offlineDataGeneric
{
public:
  offlineDataSimulationParameters() : 
    offlineDataGeneric(),
    seed(0), sqrtS(0), impactParameter(0),
    massNumberNucleusA(0), atomicNumberNucleusA(0),massNumberNucleusB(0), atomicNumberNucleusB(0),
    numberOfTestparticles(0), initialNumberOfParticles(0),firstTimeStep(0), timeShift(0),
    freezeOutEnergyDensity(0), ringStructureSize(0), ringStructureCentralRadius(0),
    ringStructureDeltaR(0), cellSizeDeltaX(0), cellSizeDeltaY(0),
    transversalSize(0), gridSizeX(0), gridSizeY(0), gridSizeZ(0) {};
  offlineDataSimulationParameters( const uint32_t _seed, const double _sqrtS, const double _impactParameter,
                                   const double _massNumberNucleusA, const double _atomicNumberNucleusA, const double _massNumberNucleusB,
                                   const double _atomicNumberNucleusB, const int _numberOfTestparticles, const int _initialNumberOfParticles,
                                   const double _firstTimeStep, const double _timeShift, const double _freezeOutEnergyDensity,
                                   const int _ringStructureSize, const double _ringStructureCentralRadius, const double _ringStructureDeltaR,
                                   const double _cellSizeDeltaX, const double _cellSizeDeltaY, const double _transversalSize,
                                   const double _gridSizeX, const double _gridSizeY, const double _gridSizeZ, const int _N_light_flav, const int _N_heavy_flav ) : 
    offlineDataGeneric(),
    seed( _seed ), sqrtS(_sqrtS),
    impactParameter(_impactParameter), massNumberNucleusA(_massNumberNucleusA), atomicNumberNucleusA(_atomicNumberNucleusA),
    massNumberNucleusB(_massNumberNucleusB), atomicNumberNucleusB(_atomicNumberNucleusB),
    numberOfTestparticles(_numberOfTestparticles), initialNumberOfParticles(_initialNumberOfParticles),
    firstTimeStep(_firstTimeStep), timeShift(_timeShift), freezeOutEnergyDensity(_freezeOutEnergyDensity),
    ringStructureSize(_ringStructureSize), ringStructureCentralRadius(_ringStructureCentralRadius),
    ringStructureDeltaR(_ringStructureDeltaR), cellSizeDeltaX(_cellSizeDeltaX), cellSizeDeltaY(_cellSizeDeltaY),
    transversalSize(_transversalSize), gridSizeX(_gridSizeX), gridSizeY(_gridSizeY), gridSizeZ(_gridSizeZ), N_light_flav(_N_light_flav), N_heavy_flav(_N_heavy_flav) {};
  ~offlineDataSimulationParameters() {};

  size_t getSize() const { return sizeof( offlineDataSimulationParameters ); }
  std::string getFilenameIdentifier() const { return this->filenameIdentifier; }
  static std::string filenameIdentifier;
    
  /** @brief initial seed of the original BAMPS run */
  uint32_t seed;
  /** @brief center of mass energy per NN pair [GeV] */
  double sqrtS;
  /** @brief impact parameter [fm] */
  double impactParameter;
  /** @brief Mass number of nucleus A  */
  double massNumberNucleusA;
  /** @brief Atomic number, i.e. number of protons, of nucleus A */
  double atomicNumberNucleusA;
  /** @brief Mass number of nucleus B  */
  double massNumberNucleusB;
  /** @brief Atomic number, i.e. number of protons, of nucleus B */
  double atomicNumberNucleusB;
  /** @brief number of test particles per real particle */
  int numberOfTestparticles;
  /** @brief total number of initial particles of the original BAMPS run */
  int initialNumberOfParticles;
  /** @brief the first time step*/
  double firstTimeStep;
  /** @brief intial time shift after the production of the particles*/
  double timeShift;
  /** @brief energy density for kinetic freeze out [GeV/fm^3] */
  double freezeOutEnergyDensity;
  /** @brief size of the ring structure (= number of rings) */
  int ringStructureSize;
  /** @brief central radius of the ring structure */
  double ringStructureCentralRadius;
  /** @brief delta R for the outer rings */
  double ringStructureDeltaR;
  /** @brief cell size in x-direction (dx) */
  double cellSizeDeltaX;
  /** @brief cell size in y-direction (dy) */
  double cellSizeDeltaY;
  /** @brief */
  double transversalSize;
  /** @brief number of cells in x-direction */
  double gridSizeX;
  /** @brief number of cells in y-direction */
  double gridSizeY;
  /** @brief number of cells in z-direction (= eta-direction) */
  double gridSizeZ;
  /** @brief number of active light flavors */
  int N_light_flav;
  /** @brief number of active heavy flavors */
  int N_heavy_flav;
    
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( offlineDataGeneric );
    ar & BOOST_SERIALIZATION_NVP( seed );
    ar & BOOST_SERIALIZATION_NVP( sqrtS );
    ar & BOOST_SERIALIZATION_NVP( impactParameter );
    ar & BOOST_SERIALIZATION_NVP( massNumberNucleusA );
    ar & BOOST_SERIALIZATION_NVP( atomicNumberNucleusA );
    ar & BOOST_SERIALIZATION_NVP( massNumberNucleusB );
    ar & BOOST_SERIALIZATION_NVP( atomicNumberNucleusB );
    ar & BOOST_SERIALIZATION_NVP( numberOfTestparticles );
    ar & BOOST_SERIALIZATION_NVP( initialNumberOfParticles );
    ar & BOOST_SERIALIZATION_NVP( firstTimeStep );
    ar & BOOST_SERIALIZATION_NVP( timeShift );
    ar & BOOST_SERIALIZATION_NVP( freezeOutEnergyDensity );
    ar & BOOST_SERIALIZATION_NVP( ringStructureSize );
    ar & BOOST_SERIALIZATION_NVP( ringStructureCentralRadius );
    ar & BOOST_SERIALIZATION_NVP( ringStructureDeltaR );
    ar & BOOST_SERIALIZATION_NVP( cellSizeDeltaX );
    ar & BOOST_SERIALIZATION_NVP( cellSizeDeltaY );
    ar & BOOST_SERIALIZATION_NVP( transversalSize );
    ar & BOOST_SERIALIZATION_NVP( gridSizeX );
    ar & BOOST_SERIALIZATION_NVP( gridSizeY );
    ar & BOOST_SERIALIZATION_NVP( gridSizeZ );
    ar & BOOST_SERIALIZATION_NVP( N_light_flav );
    ar & BOOST_SERIALIZATION_NVP( N_heavy_flav );
  }
};


namespace boost {
  namespace serialization {

    template<class Archive>
    void load(Archive & ar, vector4D & g, const unsigned int version)
    {
      // please note: the data member 'Mem' is protected and not
      // public, therefore we have to use the public set-routine.

      vector4D::Scalar Mem[4];
      ar & make_nvp( "Mem", Mem );
      g.SetTXYZ(Mem);
    }

    template<class Archive>
    void save(Archive & ar, const vector4D & g, const unsigned int version)
    {
      vector4D::Scalar Mem[4];
      g.GetTXYZ(Mem);
      ar & make_nvp( "Mem", Mem );
    }

    template<class Archive>
    void load(Archive & ar, ParticlePrototype & g, const unsigned int version)
    {
      ar & make_nvp( "unique_id", g.unique_id );
      ar & make_nvp( "cell_id", g.cell_id );
      ar & make_nvp( "FLAVOR", g.FLAVOR );
      ar & make_nvp( "m", g.m );

      switch (version)
      {
        case 0:
        {
          double x, y, z, t;

          ar & make_nvp( "T", t );
          ar & make_nvp( "X", x );
          ar & make_nvp( "Y", y );
          ar & make_nvp( "Z", z );
          g.Pos.SetTXYZ(t,x,y,z);

          ar & make_nvp( "E", t );
          ar & make_nvp( "PX", x );
          ar & make_nvp( "PY", y );
          ar & make_nvp( "PZ", z );
          g.Mom.SetTXYZ(t,x,y,z);
        } break;

        default:
        {
          ar & make_nvp( "Pos", g.Pos );
          ar & make_nvp( "Mom", g.Mom );
        } break;
      }
    }

    template<class Archive>
    void save(Archive & ar, const ParticlePrototype & g, const unsigned int version)
    {
      ar & make_nvp( "unique_id", g.unique_id );
      ar & make_nvp( "cell_id", g.cell_id );
      ar & make_nvp( "FLAVOR", g.FLAVOR );
      ar & make_nvp( "m", g.m );

      // Normally, when calling 'save' the value of version is the
      // highest possible one. For testing purposes, you also may
      // set it to a former value. see below, BOOST_CLASS_VERSION(...) 

      switch (version)
      {
        case 0:
        {
          double x, y, z, t;

          g.Pos.GetTXYZ(t,x,y,z);
          ar & make_nvp( "T", t );
          ar & make_nvp( "X", x );
          ar & make_nvp( "Y", y );
          ar & make_nvp( "Z", z );

          g.Mom.GetTXYZ(t,x,y,z);
          ar & make_nvp( "E", t );
          ar & make_nvp( "PX", x );
          ar & make_nvp( "PY", y );
          ar & make_nvp( "PZ", z );
        } break;

        default:
        {
          ar & make_nvp( "Pos", g.Pos );
          ar & make_nvp( "Mom", g.Mom );
        } break;
      }
    }

    
    /* template<class Archive> */
    /* void serialize(Archive & ar, ParticlePrototype & g, const unsigned int version) */
    /* { */
    /*   ar & g.unique_id; */
    /*   ar & g.cell_id; */
    /*   ar & g.FLAVOR; */
    /*   ar & g.m; */
    /*   ar & g.Pos; */
    /*   ar & g.Mom; */
    /* } */
    
    
    template<class Archive> 
    void load(Archive & ar, Particle & g, const unsigned int version) 
    { 
      ar >> make_nvp( "ParticlePrototype", boost::serialization::base_object<ParticlePrototype>(g) );
      ar & make_nvp( "md2g", g.md2g );
      ar & make_nvp( "md2q", g.md2q );
    }
    
    template<class Archive> 
    void save(Archive & ar, const Particle & g, const unsigned int version)  
    { 
      ar << make_nvp( "ParticlePrototype", boost::serialization::base_object<ParticlePrototype>(g) );
      ar & make_nvp( "md2g", g.md2g );
      ar & make_nvp( "md2q", g.md2q );
    } 
    
    // template<class Archive> 
    // void serialize(Archive & ar, Particle & g, const unsigned int version) 
    // { 
    //   ar & boost::serialization::base_object<ParticlePrototype>(g); 
    //   ar & g.md2g; 
    //   ar & g.md2q; 
    // }
    
  }
}
BOOST_SERIALIZATION_SPLIT_FREE(vector4D)
BOOST_SERIALIZATION_SPLIT_FREE(ParticlePrototype)
BOOST_SERIALIZATION_SPLIT_FREE(Particle)

/* BOOST_CLASS_VERSION(ParticlePrototype, 1) */
/* BOOST_CLASS_VERSION(Particle, 1) */

/**
 * @brief Class... 
 */

class offlineDataInitialParticles : public offlineDataGeneric
{
public:
  offlineDataInitialParticles() : 
    offlineDataGeneric(),
    pointerToParticleVector(0) {};
  offlineDataInitialParticles( const std::vector< Particle >* const _particles ) : 
    offlineDataGeneric(),
    pointerToParticleVector(_particles) {};
  ~offlineDataInitialParticles() {};

  size_t getSize() const { return sizeof( offlineDataInteractionRates ); }
  std::string getFilenameIdentifier() const { return this->filenameIdentifier; }
  static std::string filenameIdentifier;
                                 
  typedef boost::shared_ptr< const std::vector<Particle> > tPointerToParticleVector;
  /** @brief A shared pointer to the particle vector that is read from the archive */
  tPointerToParticleVector particleVector;
    
  friend class boost::serialization::access;
  template<class Archive>
  void save(Archive & ar, const unsigned int version) const  // split save / load operations to prevent memory leaks when loading
  {
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( offlineDataGeneric );
    ar & BOOST_SERIALIZATION_NVP( pointerToParticleVector );
  }
  template<class Archive>
  void load(Archive & ar, const unsigned int version) // split save / load operations to prevent memory leaks when loading
  {
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( offlineDataGeneric );
    ar & BOOST_SERIALIZATION_NVP( pointerToParticleVector );
    particleVector.reset( pointerToParticleVector );  
    // The outside world can only use this shared pointer to access the restored data, thus when this object goes 
    // out of scope, the memory is automatically released. This would not be the case when only providing the bare
    // pointer.
  }
  BOOST_SERIALIZATION_SPLIT_MEMBER()
    
  private:
  /** @brief A pointer to the particle vector that needs to be archived (for input / save operations) */
  const std::vector<Particle>* pointerToParticleVector;
};

// "Export" the classes derived from offlineDataGeneric in case they need to be archived via boost:serialization
// using pointers to the base class
BOOST_CLASS_EXPORT_KEY( offlineDataInteraction22 )
BOOST_CLASS_EXPORT_KEY( offlineDataInteraction23 )
BOOST_CLASS_EXPORT_KEY( offlineDataInteraction32 )
BOOST_CLASS_EXPORT_KEY( offlineDataInteractionElastic )
BOOST_CLASS_EXPORT_KEY( offlineDataCellConfiguration )
BOOST_CLASS_EXPORT_KEY( offlineDataParticleIdSwap )
BOOST_CLASS_EXPORT_KEY( offlineDataParticleNumbers )
BOOST_CLASS_EXPORT_KEY( offlineDataCollisionNumbers )
BOOST_CLASS_EXPORT_KEY( offlineDataInteractionRates )
BOOST_CLASS_EXPORT_KEY( offlineDataEventType )
BOOST_CLASS_EXPORT_KEY( offlineDataSimulationParameters )
BOOST_CLASS_EXPORT_KEY( offlineDataInitialParticles )

BOOST_CLASS_IMPLEMENTATION( offlineDataGeneric, boost::serialization::object_serializable)
BOOST_CLASS_IMPLEMENTATION( offlineDataInteraction22, boost::serialization::object_serializable)
BOOST_CLASS_IMPLEMENTATION( offlineDataInteraction23, boost::serialization::object_serializable)
BOOST_CLASS_IMPLEMENTATION( offlineDataInteraction32, boost::serialization::object_serializable)
BOOST_CLASS_IMPLEMENTATION( offlineDataInteractionElastic, boost::serialization::object_serializable)
BOOST_CLASS_IMPLEMENTATION( offlineDataCellConfiguration, boost::serialization::object_serializable)
BOOST_CLASS_IMPLEMENTATION( offlineDataParticleIdSwap, boost::serialization::object_serializable)
BOOST_CLASS_IMPLEMENTATION( offlineDataParticleNumbers, boost::serialization::object_serializable)
BOOST_CLASS_IMPLEMENTATION( offlineDataCollisionNumbers, boost::serialization::object_serializable)
BOOST_CLASS_IMPLEMENTATION( offlineDataInteractionRates, boost::serialization::object_serializable)
BOOST_CLASS_IMPLEMENTATION( offlineDataEventType, boost::serialization::object_serializable)
BOOST_CLASS_IMPLEMENTATION( offlineDataSimulationParameters, boost::serialization::object_serializable)
BOOST_CLASS_IMPLEMENTATION( offlineDataInitialParticles, boost::serialization::object_serializable)

BOOST_CLASS_TRACKING(offlineDataGeneric, boost::serialization::track_never)
BOOST_CLASS_TRACKING(offlineDataInteraction22, boost::serialization::track_never)
BOOST_CLASS_TRACKING(offlineDataInteraction23, boost::serialization::track_never)
BOOST_CLASS_TRACKING(offlineDataInteraction32, boost::serialization::track_never)
BOOST_CLASS_TRACKING(offlineDataInteractionElastic, boost::serialization::track_never)
BOOST_CLASS_TRACKING(offlineDataCellConfiguration, boost::serialization::track_never)
BOOST_CLASS_TRACKING(offlineDataParticleIdSwap, boost::serialization::track_never)
BOOST_CLASS_TRACKING(offlineDataParticleNumbers, boost::serialization::track_never)
BOOST_CLASS_TRACKING(offlineDataCollisionNumbers, boost::serialization::track_never)
BOOST_CLASS_TRACKING(offlineDataInteractionRates, boost::serialization::track_never)
BOOST_CLASS_TRACKING(offlineDataEventType, boost::serialization::track_never)
BOOST_CLASS_TRACKING(offlineDataSimulationParameters, boost::serialization::track_never)
BOOST_CLASS_TRACKING(offlineDataInitialParticles, boost::serialization::track_never)

BOOST_CLASS_IS_WRAPPER( offlineDataGeneric* )


/** @brief exception class for handling unexpected critical behaviour within simulations of heavy ion collisions  */
class eOfflineOutput_error : public std::runtime_error
{
public:
  explicit eOfflineOutput_error(const std::string& what) : std::runtime_error(what) {};
    
  virtual ~eOfflineOutput_error() throw() {};
};


#endif
