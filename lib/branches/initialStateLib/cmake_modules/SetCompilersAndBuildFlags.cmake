######### This sets custom build flags and includes options for various compilers #########

### This file has to be included and called before any calls to 'PROJECT()' or
### 'ENABLE_LANGUAGE()' in CMakeLists.txt, since we set the compiler via the 
### environment variable. This can be done only at the *initial* configuration.
###
### We also set some CMAKE_CXX_<TYPE> variables. They have to be initialized
### with the 'CACHE' keyword, since they would be overwritten by
### the first call of Modules/CMakeCXXInformation.cmake. 
### (The given descriptions are copied from this file.)
### Values given at the prompt are respected.
###
### Changing the value of CMAKE_BUILD_TYPE at the prompt is allowed.

SET(BAMPS_BUILD_TYPE "GCC" CACHE STRING "Possible build types are: GCC, Intel, IntelCSC, IntelLoeweCSC, AMD, AMDCSC, AMDLoeweCSC" )
MARK_AS_ADVANCED( BAMPS_BUILD_TYPE )

IF( NOT CMAKE_BUILD_TYPE )
  SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel" FORCE)
ENDIF( NOT CMAKE_BUILD_TYPE )


#########
## We may define own CMAKE_BUILD_TYPE:

## -DCMAKE_BUILD_TYPE=DebugProf :

SET(CMAKE_CXX_FLAGS_DEBUGPROF "-O3 -g -pg -Wall -Wunused-variable" CACHE STRING 
  "Flags used by the compiler during profiling builds.")
SET(CMAKE_C_FLAGS_DEBUGPROF "-O3 -g -pg -Wall -Wunused-variable" CACHE STRING 
  "Flags used by the compiler during profiling builds.")

## -DCMAKE_BUILD_TYPE=Backtrace :

SET(CMAKE_CXX_FLAGS_BACKTRACE "-O0 -g -rdynamic" CACHE STRING 
  "Flags used by the compiler during builds allowing backtrace.")
SET(CMAKE_C_FLAGS_BACKTRACE "-O0 -g -rdynamic" CACHE STRING 
  "Flags used by the compiler during builds allowing backtrace.")
MARK_AS_ADVANCED( CMAKE_CXX_FLAGS_BACKTRACE CMAKE_C_FLAGS_BACKTRACE)

#########

IF (BAMPS_BUILD_TYPE STREQUAL "Intel")
  MESSAGE(STATUS "** Compiling with Intel settings **")
  SET(ENV{CXX} "icpc")
  SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -no-prec-div -w -fp-model fast=2 -xHost" CACHE STRING 
    "Flags used by the compiler during release builds (/MD /Ob1 /Oi /Ot /Oy /Gs will produce slightly less optimized but smaller files).")
  SET(CMAKE_CXX_FLAGS_DEBUG "-g" CACHE STRING 
    "Flags used by the compiler during debug builds.") 
  ##
ELSEIF (BAMPS_BUILD_TYPE STREQUAL "IntelParallel")
  MESSAGE(STATUS "** Compiling with Intel parallel settings **")
  SET(ENV{CXX} "icpc")
  SET(CMAKE_CXX_FLAGS_RELEASE "-parallel -O3 -no-prec-div -w -fp-model fast=2 -xHost" CACHE STRING 
    "Flags used by the compiler during release builds (/MD /Ob1 /Oi /Ot /Oy /Gs will produce slightly less optimized but smaller files).")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "IntelLoeweCSC")
  MESSAGE(STATUS "** Compiling with Intel settings for LOEWE CSC**")
  SET(ENV{CXX} "icpc")
  SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -simd -msse3 -no-prec-div -w -fp-relaxed -fp-model fast=2 -xHost" CACHE STRING 
    "Flags used by the compiler during release builds (/MD /Ob1 /Oi /Ot /Oy /Gs will produce slightly less optimized but smaller files).")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "IntelCSC")
  MESSAGE(STATUS "** Compiling with IntelFast settings for CSC**")
  SET(ENV{CXX} "icpc")
  SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -no-prec-div -fp-relaxed -w -fp-model fast=2 -xHost" CACHE STRING 
    "Flags used by the compiler during release builds (/MD /Ob1 /Oi /Ot /Oy /Gs will produce slightly less optimized but smaller files).")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "IntelProfGen")
  MESSAGE(STATUS "** Compiling with Intel generating profile information**")
  SET(ENV{CXX} "icpc")
  SET(CMAKE_CXX_FLAGS_RELEASE "-O2 -prof-gen -w" CACHE STRING 
    "Flags used by the compiler during release builds (/MD /Ob1 /Oi /Ot /Oy /Gs will produce slightly less optimized but smaller files).")
  SET(CMAKE_CXX_FLAGS_DEBUG "-g -prof-gen" CACHE STRING 
    "Flags used by the compiler during debug builds.")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "IntelProfUse")
  MESSAGE(STATUS "** Compiling with Intel using profile information**")
  SET(ENV{CXX} "icpc")
  SET(CMAKE_CXX_FLAGS_RELEASE "-prof-use -prof-dir${CMAKE_CURRENT_BINARY_DIR} -O3 -no-prec-div -w -fp-model fast=2 -xHost" CACHE STRING 
    ".Flags used by the compiler during release builds (/MD /Ob1 /Oi /Ot /Oy /Gs will produce slightly less optimized but smaller files)...")
  SET(CMAKE_CXX_FLAGS_DEBUG "-prof-use -prof-dir${CMAKE_CURRENT_BINARY_DIR} -g" CACHE STRING 
    "Flags used by the compiler during debug builds.")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "IntelProfUseParallel")
  MESSAGE(STATUS "** Compiling with Intel using profile information**")
  SET(ENV{CXX} "icpc")
  SET(CMAKE_CXX_FLAGS_RELEASE "-prof-use -prof-dir${CMAKE_CURRENT_BINARY_DIR} -parallel -O3 -no-prec-div -w -fp-model fast=2 -xHost" CACHE STRING 
    "Flags used by the compiler during release builds (/MD /Ob1 /Oi /Ot /Oy /Gs will produce slightly less optimized but smaller files).")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "AMDLoeweCSC")
  MESSAGE(STATUS "** Compiling with AMD settings for LOEWE CSC**")
  SET(ENV{CXX} "openCC")
  SET(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG" CACHE STRING 
    "Flags used by the compiler during release builds (/MD /Ob1 /Oi /Ot /Oy /Gs will produce slightly less optimized but smaller files).")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "AMDCSC")
  MESSAGE(STATUS "** Compiling with AMD settings for LOEWE CSC**")
  SET(ENV{CXX} "openCC")
  SET(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG" CACHE STRING 
    "Flags used by the compiler during release builds (/MD /Ob1 /Oi /Ot /Oy /Gs will produce slightly less optimized but smaller files).")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "AMD")
  MESSAGE(STATUS "** Compiling with AMD settings**")
  SET(ENV{CXX} "openCC")
  SET(CMAKE_CXX_FLAGS_RELEASE "-Ofast -DNDEBUG" CACHE STRING 
    "Flags used by the compiler during release builds (/MD /Ob1 /Oi /Ot /Oy /Gs will produce slightly less optimized but smaller files).")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "GCC49")
  MESSAGE(STATUS "** Compiling with GCC 4.9 (g++) settings**")
  SET(ENV{CXX} "g++-4.9")
  ##    
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "GCC49vtune")
  MESSAGE(STATUS "** Compiling with GCC 4.9 (g++) settings + vtune**")
  SET(ENV{CXX} "g++-4.9")
  ##    
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "GCC49ompp")
  MESSAGE(STATUS "** Compiling with GCC 4.9 (g++) settings + ompp profiling **")
  SET(ENV{CXX} "kinst-ompp g++-4.9")
  ##    
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "GCC")
  MESSAGE(STATUS "** Compiling with GCC (g++) settings**")
  SET(ENV{CXX} "g++")
  ##
ELSEIF(BAMPS_BUILD_TYPE STREQUAL "CLANG")
  MESSAGE(STATUS "** Compiling with CLANG settings**")
  SET(ENV{CC} "clang")
  SET(ENV{CXX} "clang++")
  SET(CMAKE_C_FLAGS                "-Wall -std=c99")
  SET(CMAKE_C_FLAGS_DEBUG          "-g")
  SET(CMAKE_C_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
  SET(CMAKE_C_FLAGS_RELEASE        "-O4 -DNDEBUG")
  SET(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g")
  SET(CMAKE_CXX_FLAGS                "-Wall")
  SET(CMAKE_CXX_FLAGS_DEBUG          "-g")
  SET(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
  SET(CMAKE_CXX_FLAGS_RELEASE        "-v -O4 -DNDEBUG" CACHE STRING 
    "Flags used by the compiler during release builds (/MD /Ob1 /Oi /Ot /Oy /Gs will produce slightly less optimized but smaller files).")
  SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
  SET(CMAKE_AR      "/usr/bin/llvm-ar-3.4" CACHE STRING "")
  SET(CMAKE_LINKER  "/usr/bin/llvm-link-3.4" CACHE STRING "")
  SET(CMAKE_NM      "/usr/bin/llvm-nm-3.4" CACHE STRING "")
  SET(CMAKE_OBJDUMP "/usr/bin/llvm-objdump-3.4" CACHE STRING "")
  SET(CMAKE_RANLIB  "/usr/bin/llvm-ranlib-3.4" CACHE STRING "")
ELSE()
  MESSAGE(STATUS "** Compiling with default settings**")
  MESSAGE(STATUS "** Use -DBAMPS_BUILD_TYPE to specify predefined sets of build options. Possible Values are: GCC, Intel, IntelCSC, IntelLoeweCSC, AMD, AMDCSC, AMDLoeweCSC **")
ENDIF()

#########


macro(SetVectorImplType)

  # set default vector implementation:

  if (HAVE_VECTOR_AVX)
#    MESSAGE(WARNING "The AVX implementation lacks correct alignment. Setting to SSE!")
    SET(VECTOR_IMPL_TYPE "AVX" CACHE STRING "Possible vector implementation types are: Scalar, SSE, AVX" )
#    SET(VECTOR_IMPL_TYPE "SSE" CACHE STRING "Possible vector implementation types are: Scalar, SSE, AVX" )
  elseif(HAVE_VECTOR_SSE3)
    SET(VECTOR_IMPL_TYPE "SSE" CACHE STRING "Possible vector implementation types are: Scalar, SSE, AVX" )
  else()
    SET(VECTOR_IMPL_TYPE "Scalar" CACHE STRING "Possible vector implementation types are: Scalar, SSE, AVX" )
  endif()

  MARK_AS_ADVANCED( VECTOR_IMPL_TYPE )

  MESSAGE(STATUS "** Vector implementation: ${VECTOR_IMPL_TYPE}")

  # set flags for vector implementation:
  
  IF (VECTOR_IMPL_TYPE STREQUAL "AVX")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2 -msse3 -mssse3 -msse4.1 -msse4.2 -mavx -msse2avx -mno-sse4a")
  ELSEIF (VECTOR_IMPL_TYPE STREQUAL "SSE")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2 -msse3")
  ELSEIF (VECTOR_IMPL_TYPE STREQUAL "Scalar")
  ENDIF()

  #########

  # enable C++11 features:
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++0x")
  
  #########


endmacro(SetVectorImplType)

