# Try to find the BAMPS library
#
# Oliver Fochler
#
# This package defines
#  BAMPS_FOUND - The BAMPS library has been found
#  BAMPS_INCLUDE_DIRS - The directory in which the BAMPS headers reside
#  BAMPS_LIBRARIES - The BAMPS library

# find the path of the include files, based on two headers
FIND_PATH( BAMPS_INCLUDE_DIR 
    NAMES scattering22.h binary_cross_sections.h
    PATHS $ENV{HOME}/usr/include      # suggest a user based include tree
    PATH_SUFFIXES BAMPS bamps Bamps   # suggest some path suffixes in which the headers could be located
)

# find the library
FIND_LIBRARY( BAMPS_LIBRARY
    NAMES bamps
    PATHS $ENV{HOME}/usr/lib          # suggest a user based include tree
    PATH_SUFFIXES BAMPS bamps Bamps   # suggest some path suffixes in which the headers could be located
)

# adhere to the standard nomenclature for find_package routines and set some variables
SET( BAMPS_LIBRARIES ${BAMPS_LIBRARY} )
SET( BAMPS_INCLUDE_DIRS ${BAMPS_INCLUDE_DIR} )

# handle the QUIETLY and REQUIRED arguments and set BAMPS_FOUND to TRUE
# if all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS( BAMPS DEFAULT_MSG BAMPS_LIBRARY BAMPS_INCLUDE_DIR)

# display some status information
IF( BAMPS_FOUND )
 MESSAGE( STATUS "** BAMPS library: ${BAMPS_LIBRARIES}" )
 MESSAGE( STATUS "** BAMPS include: ${BAMPS_INCLUDE_DIRS}" )
ENDIF( BAMPS_FOUND )

# the variables will only show up in the GUI in the "advanced" view
MARK_AS_ADVANCED(BAMPS_INCLUDE_DIR BAMPS_LIBRARY )
