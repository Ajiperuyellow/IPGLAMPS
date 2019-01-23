#####################################
# cmake module for finding LHAPDF
#
# returns:
#   LHAPDF_FOUND
#   LHAPDF_LIBRARIES
#   LHAPDF_INCLUDE_DIR
# 
#------------------------------------
# adopted from:
# Author: Kamil Sobol
# http://code.google.com/p/winhacplusplus/source/browse/trunk/cmake/modules/FindLHAPDF.cmake
#------------------------------------
#####################################

#MESSAGE(STATUS "Looking for LHAPDF...")

# try to find LHAPDF in user defined path
FIND_LIBRARY(LHAPDF_LIB
	     NAMES LHAPDF
	     PATHS ${LHAPDF_PREFIX}
	     )

# if not try to find LHAPDF in standard instalation paths
IF(${LHAPDF_LIB} MATCHES "LHAPDF_LIB-NOTFOUND")
  FIND_LIBRARY(LHAPDF_LIB
	       NAMES LHAPDF
	       PATHS $ENV{HOME}/usr/lib /usr/lib /usr/local/lib
	       )
ENDIF(${LHAPDF_LIB} MATCHES "LHAPDF_LIB-NOTFOUND")

# if not found try to use lhapdf-config tool
IF(${LHAPDF_LIB} MATCHES "LHAPDF_LIB-NOTFOUND")
  FIND_PROGRAM(LHAPDF_CONFIG_EXECUTABLE NAMES lhapdf-config)

  IF(${LHAPDF_CONFIG_EXECUTABLE} MATCHES "LHAPDF_CONFIG_EXECUTABLE-NOTFOUND")
    MESSAGE(STATUS "Looking for LHAPDF... - lhapdf-config executable not found")
  ELSE(${LHAPDF_CONFIG_EXECUTABLE} MATCHES "LHAPDF_CONFIG_EXECUTABLE-NOTFOUND")
    MESSAGE(STATUS "Looking for LHAPDF... - using lhapdf-config executable")
    EXEC_PROGRAM(${LHAPDF_CONFIG_EXECUTABLE} ARGS "--prefix" OUTPUT_VARIABLE LHAPDF_PREFIX)
    FIND_LIBRARY(LHAPDF_LIB
		 NAMES LHAPDF
		 PATHS ${LHAPDF_PREFIX}
		 )
  ENDIF(${LHAPDF_CONFIG_EXECUTABLE} MATCHES "LHAPDF_CONFIG_EXECUTABLE-NOTFOUND")
ENDIF(${LHAPDF_LIB} MATCHES "LHAPDF_LIB-NOTFOUND")

FIND_PATH( LHAPDF_INCLUDE_DIR 
    NAMES LHAPDF.h
    PATHS $ENV{HOME}/usr/include/LHAPDF ${LHAPDF_PREFIX} /usr/include/LHAPDF /usr/local/include/LHAPDF
)

# remove a trailing 'LHAPDF', since the header files are always adressed as 'LHAPDF/xxx.h'
STRING(REGEX REPLACE "/LHAPDF" "" LHAPDF_INCLUDE_DIR ${LHAPDF_INCLUDE_DIR})


SET( LHAPDF_LIBRARIES ${LHAPDF_LIB} )
SET( LHAPDF_INCLUDE_DIRS ${LHAPDF_INCLUDE_DIR} )

# handle the QUIETLY and REQUIRED arguments and set ..._FOUND to TRUE
# if all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS( LHAPDF DEFAULT_MSG LHAPDF_LIB LHAPDF_INCLUDE_DIR)

## final printout.
IF(${LHAPDF_LIB} MATCHES "LHAPDF_LIB-NOTFOUND")
  MESSAGE(STATUS "Looking for LHAPDF... - not found !")
ELSE(${LHAPDF_LIB} MATCHES "LHAPDF_LIB-NOTFOUND")
  MESSAGE(STATUS "** LHAPDF library: ${LHAPDF_LIBRARIES}" )
  MESSAGE(STATUS "** LHAPDF include: ${LHAPDF_INCLUDE_DIRS}" )
#  MESSAGE(STATUS "** LHAPDF found: ${LHAPDF_FOUND}" )
ENDIF(${LHAPDF_LIB} MATCHES "LHAPDF_LIB-NOTFOUND")

# MESSAGE(STATUS "** LHAPDF found: ${LHAPDF_FOUND}" )

MARK_AS_ADVANCED( LHAPDF_LIB LHAPDF_INCLUDE_DIR )