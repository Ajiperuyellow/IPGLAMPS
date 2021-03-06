#=============================================================================
# Find properties of the host system:
# * SSE/AVX support
# * RDTSCP opcode support (needed for TSC timing)
#
#
# Parts are taken from "Vc" by Matthias Kretz:
# (copyright see below)
# 1) The macros "_my_find", "AutodetectHostArchitecture" 
#    are copied one-to-one (and enhanced by some
#    additional hosts)
# 2) "SetHostVectorSupport" is a boiled down version of
#    "OptimizeForArchitecture"
#=============================================================================

include(CheckCXXCompilerFlag)

# Taken from "Vc":
#=============================================================================
# Copyright 2010-2013 Matthias Kretz <kretz@kde.org>
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#  * Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
#  * The names of Kitware, Inc., the Insight Consortium, or the names of
#    any consortium members, or of any contributors, may not be used to
#    endorse or promote products derived from this software without
#    specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS ``AS IS''
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#=============================================================================

macro(_my_find _list _value _ret)
   list(FIND ${_list} "${_value}" _found)
   if(_found EQUAL -1)
      set(${_ret} FALSE)
   else(_found EQUAL -1)
      set(${_ret} TRUE)
   endif(_found EQUAL -1)
endmacro(_my_find)

macro(AutodetectHostArchitecture)
   set(TARGET_ARCHITECTURE "generic")
   set(_vendor_id)
   set(_cpu_family)
   set(_cpu_model)
   if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
      file(READ "/proc/cpuinfo" _cpuinfo)
      string(REGEX REPLACE ".*vendor_id[ \t]*:[ \t]+([a-zA-Z0-9_-]+).*" "\\1" _vendor_id "${_cpuinfo}")
      string(REGEX REPLACE ".*cpu family[ \t]*:[ \t]+([a-zA-Z0-9_-]+).*" "\\1" _cpu_family "${_cpuinfo}")
      string(REGEX REPLACE ".*model[ \t]*:[ \t]+([a-zA-Z0-9_-]+).*" "\\1" _cpu_model "${_cpuinfo}")
      string(REGEX REPLACE ".*flags[ \t]*:[ \t]+([^\n]+).*" "\\1" _cpu_flags "${_cpuinfo}")
   elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
      exec_program("/usr/sbin/sysctl -n machdep.cpu.vendor" OUTPUT_VARIABLE _vendor_id)
      exec_program("/usr/sbin/sysctl -n machdep.cpu.model"  OUTPUT_VARIABLE _cpu_model)
      exec_program("/usr/sbin/sysctl -n machdep.cpu.family" OUTPUT_VARIABLE _cpu_family)
      exec_program("/usr/sbin/sysctl -n machdep.cpu.features" OUTPUT_VARIABLE _cpu_flags)
      string(TOLOWER "${_cpu_flags}" _cpu_flags)
      string(REPLACE "." "_" _cpu_flags "${_cpu_flags}")
   elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
      get_filename_component(_vendor_id "[HKEY_LOCAL_MACHINE\\Hardware\\Description\\System\\CentralProcessor\\0;VendorIdentifier]" NAME CACHE)
      get_filename_component(_cpu_id "[HKEY_LOCAL_MACHINE\\Hardware\\Description\\System\\CentralProcessor\\0;Identifier]" NAME CACHE)
      mark_as_advanced(_vendor_id _cpu_id)
      string(REGEX REPLACE ".* Family ([0-9]+) .*" "\\1" _cpu_family "${_cpu_id}")
      string(REGEX REPLACE ".* Model ([0-9]+) .*" "\\1" _cpu_model "${_cpu_id}")
   endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
   if(_vendor_id STREQUAL "GenuineIntel")
     if(_cpu_family EQUAL 6)
         if(_cpu_model EQUAL 87)
            set(TARGET_ARCHITECTURE "knl")  # Knights Landing
         elseif(_cpu_model EQUAL 92)
            set(TARGET_ARCHITECTURE "goldmont")
         elseif(_cpu_model EQUAL 90 OR _cpu_model EQUAL 76)
            set(TARGET_ARCHITECTURE "silvermont")
         elseif(_cpu_model EQUAL 102)
            set(TARGET_ARCHITECTURE "cannonlake")
         elseif(_cpu_model EQUAL 85) # 55
            set(TARGET_ARCHITECTURE "skylake-avx512")
         elseif(_cpu_model EQUAL 78 OR _cpu_model EQUAL 94) # 4E, 5E
            set(TARGET_ARCHITECTURE "skylake")
         elseif(_cpu_model EQUAL 61 OR _cpu_model EQUAL 71 OR _cpu_model EQUAL 86)
            set(TARGET_ARCHITECTURE "broadwell")
         elseif(_cpu_model EQUAL 60 OR _cpu_model EQUAL 69 OR _cpu_model EQUAL 70 OR _cpu_model EQUAL 63)
            set(TARGET_ARCHITECTURE "haswell")
	 elseif(_cpu_model EQUAL 79) # E5-2630 
	    set(TARGET_ARCHITECTURE "sandy-bridge-ep")
         elseif(_cpu_model EQUAL 58 OR _cpu_model EQUAL 62)
            set(TARGET_ARCHITECTURE "ivy-bridge")
         elseif(_cpu_model EQUAL 42 OR _cpu_model EQUAL 45)
            set(TARGET_ARCHITECTURE "sandy-bridge")
         elseif(_cpu_model EQUAL 37 OR _cpu_model EQUAL 44 OR _cpu_model EQUAL 47)
            set(TARGET_ARCHITECTURE "westmere")
         elseif(_cpu_model EQUAL 26 OR _cpu_model EQUAL 30 OR _cpu_model EQUAL 31 OR _cpu_model EQUAL 46)
            set(TARGET_ARCHITECTURE "nehalem")
         elseif(_cpu_model EQUAL 23 OR _cpu_model EQUAL 29)
            set(TARGET_ARCHITECTURE "penryn")
         elseif(_cpu_model EQUAL 15)
            set(TARGET_ARCHITECTURE "merom")
         elseif(_cpu_model EQUAL 28)
            set(TARGET_ARCHITECTURE "atom")
         elseif(_cpu_model EQUAL 14)
            set(TARGET_ARCHITECTURE "core")
         elseif(_cpu_model LESS 14)
            message(WARNING "Your CPU (family ${_cpu_family}, model ${_cpu_model}) is not known. Auto-detection of optimization flags failed and will use the generic CPU settings with SSE2.")
            set(TARGET_ARCHITECTURE "generic")
         else()
            message(WARNING "Your CPU (family ${_cpu_family}, model ${_cpu_model}) is not known. Auto-detection of optimization flags failed and will use the 65nm Core 2 CPU settings.")
            set(TARGET_ARCHITECTURE "merom")
         endif()
      elseif(_cpu_family EQUAL 7) # Itanium (not supported)
         message(WARNING "Your CPU (Itanium: family ${_cpu_family}, model ${_cpu_model}) is not supported by OptimizeForArchitecture.cmake.")
      elseif(_cpu_family EQUAL 15) # NetBurst
         list(APPEND _available_vector_units_list "sse" "sse2")
         if(_cpu_model GREATER 2) # Not sure whether this must be 3 or even 4 instead
            list(APPEND _available_vector_units_list "sse" "sse2" "sse3")
         endif(_cpu_model GREATER 2)
      endif(_cpu_family EQUAL 6)
   elseif(_vendor_id STREQUAL "AuthenticAMD")
      if(_cpu_family EQUAL 22) # 16h
         set(TARGET_ARCHITECTURE "AMD 16h")
      elseif(_cpu_family EQUAL 21) # 15h
         if(_cpu_model LESS 2)
            set(TARGET_ARCHITECTURE "bulldozer")
         else()
            set(TARGET_ARCHITECTURE "piledriver")
         endif()
      elseif(_cpu_family EQUAL 20) # 14h
         set(TARGET_ARCHITECTURE "AMD 14h")
      elseif(_cpu_family EQUAL 18) # 12h
         if(_cpu_model EQUAL 1) # Llano
	    set(TARGET_ARCHITECTURE "llano")
         endif()
      elseif(_cpu_family EQUAL 16) # 10h
         set(TARGET_ARCHITECTURE "barcelona")
      elseif(_cpu_family EQUAL 15)
         set(TARGET_ARCHITECTURE "k8")
         if(_cpu_model GREATER 64) # I don't know the right number to put here. This is just a guess from the hardware I have access to
            set(TARGET_ARCHITECTURE "k8-sse3")
         endif(_cpu_model GREATER 64)
      endif()
   endif(_vendor_id STREQUAL "GenuineIntel")
endmacro()


macro(SetHostVectorSupport)
   set(TARGET_ARCHITECTURE "auto" CACHE STRING "CPU architecture to optimize for. Using an incorrect setting here can result in crashes of the resulting binary because of invalid instructions used.\nSetting the value to \"auto\" will try to optimize for the architecture where cmake is called.\nOther supported values are: \"none\", \"generic\", \"core\", \"merom\" (65nm Core2), \"penryn\" (45nm Core2), \"nehalem\", \"westmere\", \"sandy-bridge\", \"sandy-bridge-ep\", \"ivy-bridge\", \"haswell\", \"broadwell\", \"skylake\", \"skylake-avx512\", \"cannonlake\", \"silvermont\", \"goldmont\", \"knl\" (Knights Landing), \"atom\", \"k8\", \"k8-sse3\", \"barcelona\", \"istanbul\", \"magny-cours\", \"bulldozer\", \"interlagos\", \"piledriver\", \"AMD 14h\", \"AMD 16h\".")
   set(_force)
   if(NOT _last_target_arch STREQUAL "${TARGET_ARCHITECTURE}")
#      message(STATUS "target changed from \"${_last_target_arch}\" to \"${TARGET_ARCHITECTURE}\"")
      set(_force FORCE)
   endif()
   set(_last_target_arch "${TARGET_ARCHITECTURE}" CACHE STRING "" FORCE)
   mark_as_advanced(_last_target_arch)
   string(TOLOWER "${TARGET_ARCHITECTURE}" TARGET_ARCHITECTURE)

   set(_march_flag_list)
   set(_available_vector_units_list)
   set(_set_march 1)

   if(TARGET_ARCHITECTURE STREQUAL "auto")
      AutodetectHostArchitecture()
      message(STATUS "Detected CPU: ${TARGET_ARCHITECTURE}")
      set(_set_march 0)
   else()
      message(STATUS "Defined CPU: ${TARGET_ARCHITECTURE}")
   endif(TARGET_ARCHITECTURE STREQUAL "auto")

   macro(_nehalem)
      list(APPEND _march_flag_list "nehalem")
      list(APPEND _march_flag_list "corei7")
      list(APPEND _march_flag_list "core2")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "ssse3" "sse4.1" "sse4.2")
   endmacro()
   macro(_westmere)
      list(APPEND _march_flag_list "westmere")
      _nehalem()
   endmacro()
   macro(_sandybridge)
      list(APPEND _march_flag_list "sandybridge")
      list(APPEND _march_flag_list "corei7-avx")
      _westmere()
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "ssse3" "sse4.1" "sse4.2" "avx")
   endmacro()
   macro(_ivybridge)
      list(APPEND _march_flag_list "ivybridge")
      list(APPEND _march_flag_list "core-avx-i")
      _sandybridge()
      list(APPEND _available_vector_units_list "rdrnd" "f16c")
   endmacro()
   macro(_haswell)
      list(APPEND _march_flag_list "haswell")
      list(APPEND _march_flag_list "core-avx2")
      _ivybridge()
      list(APPEND _available_vector_units_list "avx2" "fma" "bmi" "bmi2")
   endmacro()
   macro(_broadwell)
      list(APPEND _march_flag_list "broadwell")
      _haswell()
   endmacro()
   macro(_skylake)
      list(APPEND _march_flag_list "skylake")
      _broadwell()
   endmacro()
   macro(_skylake_avx512)
      list(APPEND _march_flag_list "skylake-avx512")
      _skylake()
      list(APPEND _available_vector_units_list "avx512f" "avx512cd" "avx512dq" "avx512bw" "avx512vl")
   endmacro()
   macro(_cannonlake)
      list(APPEND _march_flag_list "cannonlake")
      _skylake_avx512()
      list(APPEND _available_vector_units_list "avx512ifma" "avx512vbmi")
   endmacro()
   macro(_knightslanding)
      list(APPEND _march_flag_list "knl")
      _broadwell()
      list(APPEND _available_vector_units_list "avx512f" "avx512pf" "avx512er" "avx512cd")
   endmacro()
   macro(_silvermont)
      list(APPEND _march_flag_list "silvermont")
      _westmere()
      list(APPEND _available_vector_units_list "rdrnd")
   endmacro()
   macro(_goldmont)
      list(APPEND _march_flag_list "goldmont")
      _silvermont()
    endmacro()
    
   if(TARGET_ARCHITECTURE STREQUAL "core")
      list(APPEND _march_flag_list "core2")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3")
   elseif(TARGET_ARCHITECTURE STREQUAL "merom")
      list(APPEND _march_flag_list "merom")
      list(APPEND _march_flag_list "core2")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "ssse3")
   elseif(TARGET_ARCHITECTURE STREQUAL "penryn")
      list(APPEND _march_flag_list "penryn")
      list(APPEND _march_flag_list "core2")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "ssse3")
      message(STATUS "Sadly the Penryn architecture exists in variants with SSE4.1 and without SSE4.1.")
      if(_cpu_flags MATCHES "sse4_1")
         message(STATUS "SSE4.1: enabled (auto-detected from this computer's CPU flags)")
         list(APPEND _available_vector_units_list "sse4.1")
      else()
         message(STATUS "SSE4.1: disabled (auto-detected from this computer's CPU flags)")
       endif()
   elseif(TARGET_ARCHITECTURE STREQUAL "knl")
      _knightslanding()
   elseif(TARGET_ARCHITECTURE STREQUAL "cannonlake")
      _cannonlake()
   elseif(TARGET_ARCHITECTURE STREQUAL "skylake-xeon" OR TARGET_ARCHITECTURE STREQUAL "skylake-avx512")
      _skylake_avx512()
   elseif(TARGET_ARCHITECTURE STREQUAL "skylake")
      _skylake()
   elseif(TARGET_ARCHITECTURE STREQUAL "broadwell")
      _broadwell()
   elseif(TARGET_ARCHITECTURE STREQUAL "haswell")
      _haswell()
   elseif(TARGET_ARCHITECTURE STREQUAL "ivy-bridge")
      _ivybridge()
   elseif(TARGET_ARCHITECTURE STREQUAL "sandy-bridge")
      _sandybridge()
   elseif(TARGET_ARCHITECTURE STREQUAL "westmere")
      _westmere()
   elseif(TARGET_ARCHITECTURE STREQUAL "nehalem")
      _nehalem()
   elseif(TARGET_ARCHITECTURE STREQUAL "goldmont")
      _goldmont()
   elseif(TARGET_ARCHITECTURE STREQUAL "silvermont")
      _silvermont()
    elseif(TARGET_ARCHITECTURE STREQUAL "sandy-bridge-ep")
      list(APPEND _march_flag_list "sandybridge")
      list(APPEND _march_flag_list "core-avx2")
      list(APPEND _march_flag_list "core-avx-i")
      list(APPEND _march_flag_list "corei7-avx")
      list(APPEND _march_flag_list "core2")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "ssse3" "sse4.1" "sse4.2" "avx" "avx2" "rdrnd" "f16c" "fma")
    elseif(TARGET_ARCHITECTURE STREQUAL "atom")
      list(APPEND _march_flag_list "atom")
      list(APPEND _march_flag_list "bonnell")
      list(APPEND _march_flag_list "core2")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "ssse3")
   elseif(TARGET_ARCHITECTURE STREQUAL "k8")
      list(APPEND _march_flag_list "k8")
      list(APPEND _available_vector_units_list "sse" "sse2")
   elseif(TARGET_ARCHITECTURE STREQUAL "k8-sse3")
      list(APPEND _march_flag_list "k8-sse3")
      list(APPEND _march_flag_list "k8")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3")
   elseif(TARGET_ARCHITECTURE STREQUAL "AMD 16h")
      list(APPEND _march_flag_list "btver2")
      list(APPEND _march_flag_list "btver1")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "ssse3" "sse4a" "sse4.1" "sse4.2" "avx" "f16c")
   elseif(TARGET_ARCHITECTURE STREQUAL "AMD 14h")
      list(APPEND _march_flag_list "btver1")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "ssse3" "sse4a")
   elseif(TARGET_ARCHITECTURE STREQUAL "piledriver")
      list(APPEND _march_flag_list "bdver2")
      list(APPEND _march_flag_list "bdver1")
      list(APPEND _march_flag_list "bulldozer")
      list(APPEND _march_flag_list "barcelona")
      list(APPEND _march_flag_list "core2")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "ssse3" "sse4a" "sse4.1" "sse4.2" "avx" "xop" "fma4" "fma" "f16c")
   elseif(TARGET_ARCHITECTURE STREQUAL "interlagos")
      list(APPEND _march_flag_list "bdver1")
      list(APPEND _march_flag_list "bulldozer")
      list(APPEND _march_flag_list "barcelona")
      list(APPEND _march_flag_list "core2")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "ssse3" "sse4a" "sse4.1" "sse4.2" "avx" "xop" "fma4")
   elseif(TARGET_ARCHITECTURE STREQUAL "bulldozer")
      list(APPEND _march_flag_list "bdver1")
      list(APPEND _march_flag_list "bulldozer")
      list(APPEND _march_flag_list "barcelona")
      list(APPEND _march_flag_list "core2")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "ssse3" "sse4a" "sse4.1" "sse4.2" "avx" "xop" "fma4")
   elseif(TARGET_ARCHITECTURE STREQUAL "barcelona")
      list(APPEND _march_flag_list "barcelona")
      list(APPEND _march_flag_list "core2")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "sse4a")
   elseif(TARGET_ARCHITECTURE STREQUAL "istanbul")
      list(APPEND _march_flag_list "barcelona")
      list(APPEND _march_flag_list "core2")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "sse4a")
   elseif(TARGET_ARCHITECTURE STREQUAL "llano")
      list(APPEND _march_flag_list "amdfam10")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "sse4a")
   elseif(TARGET_ARCHITECTURE STREQUAL "magny-cours")
      list(APPEND _march_flag_list "barcelona")
      list(APPEND _march_flag_list "core2")
      list(APPEND _available_vector_units_list "sse" "sse2" "sse3" "sse4a")
   elseif(TARGET_ARCHITECTURE STREQUAL "generic")
      list(APPEND _march_flag_list "generic")
   elseif(TARGET_ARCHITECTURE STREQUAL "none")
      # add this clause to remove it from the else clause
   else(TARGET_ARCHITECTURE STREQUAL "core")
      message(FATAL_ERROR "Unknown target architecture: \"${TARGET_ARCHITECTURE}\". Please set TARGET_ARCHITECTURE to a supported value.")
   endif(TARGET_ARCHITECTURE STREQUAL "core")



   if(NOT TARGET_ARCHITECTURE STREQUAL "none")
      _my_find(_available_vector_units_list "sse2" SSE2_FOUND)
      _my_find(_available_vector_units_list "sse3" SSE3_FOUND)
      _my_find(_available_vector_units_list "ssse3" SSSE3_FOUND)
      _my_find(_available_vector_units_list "sse4.1" SSE4_1_FOUND)
      _my_find(_available_vector_units_list "sse4.2" SSE4_2_FOUND)
      _my_find(_available_vector_units_list "sse4a" SSE4a_FOUND)
      if(DEFINED Vc_AVX_INTRINSICS_BROKEN AND Vc_AVX_INTRINSICS_BROKEN)
         UserWarning("AVX disabled per default because of old/broken compiler")
         set(AVX_FOUND false)
         set(XOP_FOUND false)
         set(FMA4_FOUND false)
      else()
         _my_find(_available_vector_units_list "avx" AVX_FOUND)
         if(DEFINED Vc_FMA4_INTRINSICS_BROKEN AND Vc_FMA4_INTRINSICS_BROKEN)
            UserWarning("FMA4 disabled per default because of old/broken compiler")
            set(FMA4_FOUND false)
         else()
            _my_find(_available_vector_units_list "fma4" FMA4_FOUND)
         endif()
         if(DEFINED Vc_XOP_INTRINSICS_BROKEN AND Vc_XOP_INTRINSICS_BROKEN)
            UserWarning("XOP disabled per default because of old/broken compiler")
            set(XOP_FOUND false)
         else()
            _my_find(_available_vector_units_list "xop" XOP_FOUND)
         endif()
      endif()

      set(HAVE_VECTOR_SSE2   ${SSE2_FOUND}   CACHE BOOL "We have support of SSE2. If SSE2 instructions are not enabled the SSE implementation will be disabled." ${_force})
      set(HAVE_VECTOR_SSE3   ${SSE3_FOUND}   CACHE BOOL "We have support of SSE3. If SSE3 instructions are not enabled they will be emulated." ${_force})
      set(HAVE_VECTOR_SSSE3  ${SSSE3_FOUND}  CACHE BOOL "We have support of SSSE3. If SSSE3 instructions are not enabled they will be emulated." ${_force})
      set(HAVE_VECTOR_SSE4_1 ${SSE4_1_FOUND} CACHE BOOL "We have support of SSE4.1. If SSE4.1 instructions are not enabled they will be emulated." ${_force})
      set(HAVE_VECTOR_SSE4_2 ${SSE4_2_FOUND} CACHE BOOL "We have support of SSE4.2. If SSE4.2 instructions are not enabled they will be emulated." ${_force})
      set(HAVE_VECTOR_SSE4a  ${SSE4a_FOUND}  CACHE BOOL "We have support of SSE4a. If SSE4a instructions are not enabled they will be emulated." ${_force})
      set(HAVE_VECTOR_AVX    ${AVX_FOUND}    CACHE BOOL "We have support of AVX. This will double some of the vector sizes relative to SSE." ${_force})
      set(HAVE_VECTOR_XOP    ${XOP_FOUND}    CACHE BOOL "We have support of XOP." ${_force})
      set(HAVE_VECTOR_FMA4   ${FMA4_FOUND}   CACHE BOOL "We have support of FMA4." ${_force})


      set(_cxx_flags)

      if(_set_march)

        if(CMAKE_CXX_COMPILER MATCHES "/(icpc|icc)$") # ICC (on Linux)
          
          # to be inserted...
          
        else() # not MSVC and not ICC => GCC, Clang, Open64
          foreach(_flag ${_march_flag_list})
            string(REGEX REPLACE "[-.+/:= ]" "_" _flag_esc "${_flag}")
            CHECK_CXX_COMPILER_FLAG("-march=${_flag}" HAS_FLAG_${_flag_esc})
            if(HAS_FLAG_${_flag_esc})
#              message(STATUS "Compiler supports ${_flag}")
              set(_cxx_flags "${_cxx_flags} -march=${_flag}")
            else()
#              message(STATUS "Compiler does not supports ${_flag}")
            endif()
            
          endforeach(_flag)
        endif()

      else()

#        message(STATUS "setting flag -march=native")
#        set(_cxx_flags "${_cxx_flags} -march=native")

      endif(_set_march)


      set(CXX_FLAGS_MARCH ${_cxx_flags} CACHE STRING "additional CXX flags -march=..." FORCE)


  endif ()

  # the variables will only show up in the GUI in the "advanced" view
  MARK_AS_ADVANCED(CXX_FLAGS_MARCH)

endmacro()

macro(SetHostRDTSCP)
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    execute_process(COMMAND cat /proc/cpuinfo
      COMMAND grep flags
      COMMAND grep -c rdtscp
      OUTPUT_VARIABLE HAVE_RDTSCP)
  elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    execute_process(COMMAND /usr/sbin/sysctl -n machdep.cpu.features
      COMMAND grep -c rdtscp
      OUTPUT_VARIABLE HAVE_RDTSCP)
  endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set(HAVE_RDTSCP ${HAVE_RDTSCP} CACHE INTEGER "CPU has RDTSCP opcode (indicates also number of cores)" FORCE)
endmacro()

macro(SetHostRDRAND)
  execute_process(COMMAND cat /proc/cpuinfo
    COMMAND grep flags
    COMMAND grep -c rdrand
    OUTPUT_VARIABLE HAVE_RDRAND)

  set(HAVE_RDRAND ${HAVE_RDRAND} CACHE INTEGER "CPU has RDRAND opcode (indicates also number of cores)" FORCE)
endmacro()

### 
# now DO the macros:
###

SetHostVectorSupport()
SetHostRDTSCP()
SetHostRDRAND()

###
