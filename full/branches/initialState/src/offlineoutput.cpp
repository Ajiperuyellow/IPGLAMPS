//--------------------------------------------------------- -*- c++ -*- ------
//provided by subversion
//----------------------------------------------------------------------------
//$HeadURL: file:///home/bamps/svn/full/branches/initialState/src/offlineoutput.cpp $
//$LastChangedDate: 2014-11-16 21:18:55 +0100 (日, 16 11月 2014) $
//$LastChangedRevision: 1935 $
//$LastChangedBy: gallmei $
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#include <typeinfo>
#include <iostream>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "offlineoutput.h"

using std::type_index;

std::string  offlineDataGeneric::filenameIdentifier = "oflineDataGeneric";
std::string  offlineDataEventType::filenameIdentifier = "eventList";
std::string  offlineDataInteraction22::filenameIdentifier = "interaction22";
std::string  offlineDataInteraction23::filenameIdentifier = "interaction23";
std::string  offlineDataInteraction32::filenameIdentifier = "interaction32";
std::string  offlineDataInteractionElastic::filenameIdentifier = "interactionElastic";
std::string  offlineDataCellConfiguration::filenameIdentifier = "cellConfiguration";
std::string  offlineDataParticleIdSwap::filenameIdentifier = "particleIdSwap";
std::string  offlineDataParticleNumbers::filenameIdentifier = "particleNumbers";
std::string  offlineDataCollisionNumbers::filenameIdentifier = "collisionNumbers";
std::string  offlineDataInteractionRates::filenameIdentifier = "interactionRates";
std::string  offlineDataSimulationParameters::filenameIdentifier = "simulationParameters";
std::string  offlineDataInitialParticles::filenameIdentifier = "initialParticles";


offlineOutputInterface::~offlineOutputInterface()
{
  collectAndOutputEventList();
  
  archiveMap.clear(); // due to the usage of shared pointers this is sufficient to take care of all cleanup
}



void offlineOutputInterface::collectAndOutputEventList()
{
  for ( unsigned int i = 0; i < eventList.size(); i++ )
  {
    offlineDataEventType collectEventType( eventList[i] );
    this->submitOfflineDataForOutput( &collectEventType );
  }
}


void offlineOutputInterface::submitOfflineDataForOutput( const offlineDataGeneric*const _data )
{
  tPointerToFilestream stream;
  tPointerToArchive archive;
 
  //check whether an output archive for the given type of data has already been created, re-use if yes, create if no
  tArchiveMap::iterator it = archiveMap.find( type_index(typeid(*_data)) );
  if( it != archiveMap.end() )
  {
//     stream = it->second;
    archive = it->second;
  }
  else
  {
    stream.reset( new boost::filesystem::ofstream() );
    std::string filename = outputDirectory.string() + "/" + filenamePrefix + additionalFilenameTag + _data->getFilenameIdentifier() + filenameSuffix;
    
    #if BINARY_OFFLINE_OUTPUT > 0
      stream->open( filename.c_str(), std::ios::trunc | std::ios::binary );
    #else
      stream->open( filename.c_str(), std::ios::trunc );
    #endif
    
    outputMap.insert(std::make_pair( type_index(typeid(*_data)), stream));
    archive.reset( new tArchive( *stream ) );
    archiveMap.insert(std::make_pair( type_index(typeid(*_data)), archive));
  }
  
  // the actual output
  (*archive) & BOOST_SERIALIZATION_NVP( _data );
}


void offlineOutputInterface::registerOfflineDataForTemporaryStorage(const offlineEventType _eventType, tConstPointerToOfflineData _data)
{
  this->registerOfflineDataForTemporaryStorage( tConstPointerToOfflineData( new offlineDataEventType(_eventType) ) );
  this->registerOfflineDataForTemporaryStorage( _data );
}


void offlineOutputInterface::outputAndResetTemporaryStorage()
{
  for ( tTemporaryDataStorage::const_iterator it = temporaryStorage.begin(); it != temporaryStorage.end(); it++ )
  {
    //check whether the given data is actually just an offlineEventType and give it some special treatment 
    if ( type_index(typeid(*((*it).get()))) == type_index( typeid(offlineDataEventType) ) )
    {
      const offlineDataEventType* _tmpPtr = dynamic_cast<const offlineDataEventType*>( (*it).get() );
      submitOfflineDataForOutput( _tmpPtr->event );
    }
    else
    {
      this->submitOfflineDataForOutput( (*it).get() ); //submit a pointer to the object stored at the iterator it
    }
  }
  
  temporaryStorage.clear();
}


BOOST_CLASS_EXPORT_IMPLEMENT( offlineDataInteraction22 )
BOOST_CLASS_EXPORT_IMPLEMENT( offlineDataInteraction23 )
BOOST_CLASS_EXPORT_IMPLEMENT( offlineDataInteraction32 )
BOOST_CLASS_EXPORT_IMPLEMENT( offlineDataInteractionElastic )
BOOST_CLASS_EXPORT_IMPLEMENT( offlineDataCellConfiguration )
BOOST_CLASS_EXPORT_IMPLEMENT( offlineDataParticleIdSwap )
BOOST_CLASS_EXPORT_IMPLEMENT( offlineDataParticleNumbers )
BOOST_CLASS_EXPORT_IMPLEMENT( offlineDataCollisionNumbers )
BOOST_CLASS_EXPORT_IMPLEMENT( offlineDataInteractionRates )
BOOST_CLASS_EXPORT_IMPLEMENT( offlineDataEventType )
BOOST_CLASS_EXPORT_IMPLEMENT( offlineDataSimulationParameters )
BOOST_CLASS_EXPORT_IMPLEMENT( offlineDataInitialParticles )
