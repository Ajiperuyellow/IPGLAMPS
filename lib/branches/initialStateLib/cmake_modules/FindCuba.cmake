# Try to find the Cuba library
#
# Oliver Fochler
#
# This package defines
#  Cuba_FOUND - The Cuba library has been found
#  Cuba_INCLUDE_DIRS - The directory in which the Cuba headers reside
#  Cuba_LIBRARIES - The Cuba library

# find the path of the include files, based on two headers
FIND_PATH( Cuba_INCLUDE_DIR 
    NAMES cuba.h
    PATHS $ENV{HOME}/usr/include      # suggest a user based include tree
    PATH_SUFFIXES cuba CUBA Cuba   # suggest some path suffixes in which the headers could be located
)

# find the library
FIND_LIBRARY( Cuba_LIBRARY
    NAMES cuba
    PATHS $ENV{HOME}/usr/lib          # suggest a user based include tree
    PATH_SUFFIXES cuba CUBA Cuba   # suggest some path suffixes in which the headers could be located
)

# adhere to the standard nomenclature for find_package routines and set some variables
SET( Cuba_LIBRARIES ${Cuba_LIBRARY} )
SET( Cuba_INCLUDE_DIRS ${Cuba_INCLUDE_DIR} )

# handle the QUIETLY and REQUIRED arguments and set Cuba_FOUND to TRUE
# if all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS( Cuba DEFAULT_MSG Cuba_LIBRARY Cuba_INCLUDE_DIR)

# display some status information
SET ( Cuba_FOUND ${CUBA_FOUND} CACHE INTERNAL "Provide Cuba_FOUND in addition to CUBA_FOUND" FORCE )
IF( Cuba_FOUND )
  IF (NOT Cuba_FIND_QUIETLY)
    MESSAGE( STATUS "** Cuba library: ${Cuba_LIBRARIES}" )
    MESSAGE( STATUS "** Cuba include: ${Cuba_INCLUDE_DIRS}" )
  ENDIF (NOT Cuba_FIND_QUIETLY)
ELSE (Cuba_FOUND)
  IF (Cuba_FIND_REQUIRED)
    MESSAGE(FATAL_ERROR "Could not find cuba")
  ENDIF (Cuba_FIND_REQUIRED)
ENDIF( Cuba_FOUND )

# the variables will only show up in the GUI in the "advanced" view
MARK_AS_ADVANCED(Cuba_INCLUDE_DIR Cuba_LIBRARY )
