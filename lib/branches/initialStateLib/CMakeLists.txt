CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules/")

######### Set custom build flags and options for various compilers #########
INCLUDE(SetCompilersAndBuildFlags)

######### Set project name #########
PROJECT(BAMPS_Library)

######### Activate packaging, use as "make package_source" #########
INCLUDE(CMakeCustomPackageOptions.cmake)
INCLUDE(CPack)

######### Setup host properties (Vector, RDTSCP, ...) #########
INCLUDE(HostProperties)
SetVectorImplType()
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_FLAGS_MARCH}")

######### check that C++11 code compiles #########
INCLUDE(CheckCXXSourceCompiles)
CHECK_CXX_SOURCE_COMPILES("
#include <initializer_list>
struct A {};
struct B : public A { using A::A; };
template<typename T = void> int f() { return 1; }
int main() {
  auto l = []() { return f(); };
  return l() - 1;
}
" CXX11_COMPILES)

IF( NOT CXX11_COMPILES )
#   message(FATAL_ERROR "Your compiler did not manage to compile a simple C++11 program. Please get a newer C++ compiler.")
   message(WARNING "Your compiler did not manage to compile a simple C++11 program. Please get a newer C++ compiler."
   "It may be okay for this code. Keep fingers crossed that it manages the C++11 features used.")
ENDIF()

######### check for -fdiagnostics-color compiler option #########
INCLUDE(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-fdiagnostics-color=always" HAS_DIAGNOSTICSCOLOR)
IF( HAS_DIAGNOSTICSCOLOR )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fdiagnostics-color=always")
ENDIF()

######### Find Boost libraries #########

SET(Boost_ADDITIONAL_VERSIONS "1.40" "1.40.0" "1.41" "1.41.0" "1.45" "1.45.0" "1.45.1" "1.46" "1.46.0" "1.46.1")
IF(NOT BOOST_ROOT AND NOT DEFINED ENV{BOOST_ROOT})
  SET(BOOST_ROOT "$ENV{HOME}/usr")
ENDIF()
MESSAGE(STATUS "** Boost Hint: ${BOOST_ROOT}")

FIND_PACKAGE(Boost 1.40 COMPONENTS program_options system filesystem regex REQUIRED)
MESSAGE(STATUS "** Boost Include: ${Boost_INCLUDE_DIR}")
MESSAGE(STATUS "** Boost Libraries: ${Boost_LIBRARY_DIRS}")

######### Find Cuba library #########
FIND_PACKAGE(Cuba)

########## Find GSL library #########
FIND_PACKAGE(GSL REQUIRED)

########## Find ROOT library #########
#FIND_PACKAGE(ROOT REQUIRED)
#MESSAGE(STATUS "** ROOT Include: ${ROOT_INCLUDE_DIR}")

########## This otion controls whther Log4cxx should be searched #########
## Default: ON
## Change on the command line via -DWITH_LOG4CXX=OFF
OPTION( WITH_LOG4CXX "Whether try to use Log4cxx or not" ON)
########## Find Log4cxx library #########
IF( WITH_LOG4CXX )
  IF( NOT LOG4CXX_ROOT AND NOT DEFINED ENV{LOG4CXX_ROOT} )
    SET(LOG4CXX_ROOT "$ENV{HOME}/usr")
  ENDIF()
  FIND_PACKAGE(Log4cxx)
ENDIF()

######### This option controls whether optional goodies (subversion revison number and doxygen generation) should be disabled #########
## Default: FALSE
## Change on the command line via -DDISABLE_BUILD_EXTRAS=TRUE
OPTION( DISABLE_BUILD_EXTRAS "Whether to disable some goodies [svn rev number and doxygen targets]" OFF )

######### This option controls whether OpenMP should be used #########
## Default: FALSE
## Change on the command line via -DWITH_OPENMP=ON
OPTION( WITH_OPENMP "Whether to use OpenMP or not" OFF )

######### Find OpenMP and enable/disable it according flag #########

IF( WITH_OPENMP )
  FIND_PACKAGE( OpenMP )
ENDIF()
IF( OPENMP_FOUND )
  #    MESSAGE(STATUS "** OpenMP found and will be used.")
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -D_GLIBCXX_PARALLEL")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
ELSE()
  #    MESSAGE(STATUS "** OpenMP found, but will *not* be used.")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
ENDIF()

######### Find svn and retrieve revision number #########
IF( NOT DISABLE_BUILD_EXTRAS )
  FIND_PACKAGE(Subversion)
  IF( Subversion_FOUND )
    Subversion_WC_INFO(${CMAKE_CURRENT_SOURCE_DIR} ER)
    SET(SUBVERSION_REVISION ${ER_WC_REVISION})
    SET(CPACK_PACKAGE_VERSION_PATCH "${SUBVERSION_REVISION}")
  ENDIF()
ENDIF()
#SET(Subversion_FOUND FALSE)

######### Handle doxygen files via CMake, adds custom target, use as "make doxygen" #########
IF( NOT DISABLE_BUILD_EXTRAS )
  INCLUDE(ConfigureDoxygenTarget.cmake)
ENDIF( NOT DISABLE_BUILD_EXTRAS )

######### Set include and link directories #########
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
IF( Cuba_FOUND )
  INCLUDE_DIRECTORIES(${Cuba_INCLUDE_DIRS})
ENDIF()
IF( LOG4CXX_FOUND )
  INCLUDE_DIRECTORIES(${LOG4CXX_INCLUDE_DIRS})
ENDIF( )
INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS})
IF( ROOT_FOUND )
  INCLUDE_DIRECTORIES(${ROOT_INCLUDE_DIR})
ENDIF()

######### Set an option to enable the generation of install targets #########
OPTION( GENERATE_BAMPS_LIB_INSTALL_TARGETS "Whether to generate install targets for the BAMPS library or not" ON )

######### Enable Testing  #########
ENABLE_TESTING()

# enable dashboard scripting
INCLUDE(CTest)

######### Add source directories  #########
ADD_SUBDIRECTORY(tests)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(progs)


