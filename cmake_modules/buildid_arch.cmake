# FEAST BUILD_ID mechanism:
#       detection of underlying hardware microarchitecture
#
# This module sets the following variable:
#   FEAST_CPU_TYPE
#
# TODO implement more "guessing" macros for CUDA, CL, whatever
#
# Note to developers:
#   - The two macro functions must be kept synchroneous.
#   - See bottom of this file for the actual logic.
#   - Using cmake functions instead of macros doesn't work
#
# \author Dominik Goeddeke
# \author Dirk Ribbrock



# small macro to parse the arch token in a given BUILD_ID
macro (parse_arch_token)

  # Intel CPUs
  if (BUILD_ID MATCHES "^i486-.+|.+-i486-.+|.+-i486$")
    set(FEAST_CPU_TYPE "i486")
  elseif (BUILD_ID MATCHES "^pentium-.+|.+-pentium-.+|.+-pentium$")
    set(FEAST_CPU_TYPE "pentium")
  elseif (BUILD_ID MATCHES "^pentiumpro-.+|.+-pentiumpro-.+|.+-pentiumpro$")
    set(FEAST_CPU_TYPE "pentiumpro")
  elseif (BUILD_ID MATCHES "^pentium2-.+|.+-pentium2-.+|.+-pentium2$")
    set(FEAST_CPU_TYPE "pentium2")
  elseif (BUILD_ID MATCHES "^pentium3-.+|.+-pentium3-.+|.+-pentium3$")
    set(FEAST_CPU_TYPE "pentium3")
  elseif (BUILD_ID MATCHES "^pentiumm-.+|.+-pentiumm-.+|.+-pentiumm$")
    set(FEAST_CPU_TYPE "pentiumm")
  elseif (BUILD_ID MATCHES "^pentium4m-.+|.+-pentium4m-.+|.+-pentium4m$")
    set(FEAST_CPU_TYPE "pentium4m")
  elseif (BUILD_ID MATCHES "^coresolo-.+|.+-coresolo-.+|.+-coresolo$")
    set(FEAST_CPU_TYPE "coresolo")
  elseif (BUILD_ID MATCHES "^coreduo-.+|.+-coreduo-.+|.+-coreduo$")
    set(FEAST_CPU_TYPE "coreduo")
  elseif (BUILD_ID MATCHES "^penryn-.+|.+-penryn-.+|.+-penryn$")
    set(FEAST_CPU_TYPE "penryn")
  elseif (BUILD_ID MATCHES "^nehalem-.+|.+-nehalem-.+|.+-nehalem$")
    set(FEAST_CPU_TYPE "nehalem")
  elseif (BUILD_ID MATCHES "^westmere-.+|.+-westmere-.+|.+-westmere$")
    set(FEAST_CPU_TYPE "westmere")
  elseif (BUILD_ID MATCHES "^itanium-.+|.+-itanium-.+|.+-itanium$")
    set(FEAST_CPU_TYPE "itanium")
  elseif (BUILD_ID MATCHES "^pentium4-.+|.+-pentium4-.+|.+-pentium4$")
    set(FEAST_CPU_TYPE "pentium4")
  elseif (BUILD_ID MATCHES "^nocona-.+|.+-nocona-.+|.+-nocona$")
    set(FEAST_CPU_TYPE "nocona")
  elseif (BUILD_ID MATCHES "^itanium2-.+|.+-itanium2-.+|.+-itanium2$")
    set(FEAST_CPU_TYPE "itanium2")

  # AMD CPUs
  elseif (BUILD_ID MATCHES "^amd486-.+|.+-amd486-.+|.+-amd486$")
    set(FEAST_CPU_TYPE "amd486")
  elseif (BUILD_ID MATCHES "^k5-.+|.+-k5-.+|.+-k5$")
    set(FEAST_CPU_TYPE "k5")
  elseif (BUILD_ID MATCHES "^k6-.+|.+-k6-.+|.+-k6$")
    set(FEAST_CPU_TYPE "k6")
  elseif (BUILD_ID MATCHES "^athlon-.+|.+-athlon-.+|.+-athlon$")
    set(FEAST_CPU_TYPE "athlon")
  elseif (BUILD_ID MATCHES "^athlonxp-.+|.+-athlonxp-.+|.+-athlonxp$")
    set(FEAST_CPU_TYPE "athlonxp")
  elseif (BUILD_ID MATCHES "^opteron-.+|.+-opteron-.+|.+-opteron$")
    set(FEAST_CPU_TYPE "opteron")
  elseif (BUILD_ID MATCHES "^athlon64-.+|.+-athlon64-.+|.+-athlon64$")
    set(FEAST_CPU_TYPE "athlon64")
  elseif (BUILD_ID MATCHES "^opteronx2-.+|.+-opteronx2-.+|.+-opteronx2$")
    set(FEAST_CPU_TYPE "opteronx2")
  elseif (BUILD_ID MATCHES "^turionx2-.+|.+-turionx2-.+|.+-turionx2$")
    set(FEAST_CPU_TYPE "turionx2")
  elseif (BUILD_ID MATCHES "^barcelona-.+|.+-barcelona-.+|.+-barcelona$")
    set(FEAST_CPU_TYPE "barcelona")

  # sparc CPUs
  elseif (BUILD_ID MATCHES "^sparc-.+|.+-sparc-.+|.+-sparc$")
    set(FEAST_CPU_TYPE "sparc")
  else()
    set(FEAST_CPU_TYPE "FEAST_CPU_TYPE-NOTFOUND")
  endif ()
endmacro (parse_arch_token)



# macro to guess the CPU type
# Based on code found on http://code.compeng.uni-frankfurt.de/
# and on guess_id() in the old FEAST
macro (guess_cpu)
  # initialise necessary variables
  set (_vendor_id)
  set (_cpu_family)
  set (_cpu_model)

  # extract CPU information via system calls
  if (CMAKE_SYSTEM_NAME STREQUAL "Linux")

    # on Linux, CPU information is available via /proc/cpuinfo
    file(READ "/proc/cpuinfo" _cpuinfo)
    string(REGEX REPLACE ".*vendor_id[ \t]*:[ \t]+([a-zA-Z0-9_-]+).*" "\\1" _vendor_id "${_cpuinfo}")
    string(REGEX REPLACE ".*cpu family[ \t]*:[ \t]+([a-zA-Z0-9_-]+).*" "\\1" _cpu_family "${_cpuinfo}")
    string(REGEX REPLACE ".*model[ \t]*:[ \t]+([a-zA-Z0-9_-]+).*" "\\1" _cpu_model "${_cpuinfo}")
    string(REGEX REPLACE ".*flags[ \t]*:[ \t]+([a-zA-Z0-9_-]+).*" "\\1" _cpu_flags "${_cpuinfo}")

  elseif (CMAKE_SYSTEM_NAME STREQUAL "SunOS")
    # TODO daraus ein Ticket machen
    # for SunOS, Dominik has no idea how to query CPU information:
    # according to a one-hour google search (and discussions with Christian Becker
    # and Andreas Langer), the tools prtconf, prtdiag and (on some systems)
    # prtinfo, all with -v for verbose output, should give what we need.
    # Unfortunately, these tools do not provide any useful information on the SunOS
    # machines Dominik tested (Apr. 13 2011, peppone, daontooine, naboo), so
    # the current implementation defaults to generic builds for SunOS. One could argue
    # that these machines are not really our target platform, but: For massively
    # shared memory multithreading, we have no other choice.
    set (FEAST_CPU_TYPE "sparc")
    message (STATUS "WARNING: SunOs support is only half-baked")

  elseif (CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    # on Mac/Darvin, CPU information can be queried with the small tool
    # /usr/sbin/sysctl
    exec_program("/usr/sbin/sysctl -n machdep.cpu.vendor" OUTPUT_VARIABLE _vendor_id)
    exec_program("/usr/sbin/sysctl -n machdep.cpu.model"  OUTPUT_VARIABLE _cpu_model)
    exec_program("/usr/sbin/sysctl -n machdep.cpu.family" OUTPUT_VARIABLE _cpu_family)
    exec_program("/usr/sbin/sysctl -n machdep.cpu.features" OUTPUT_VARIABLE _cpu_flags)
    string(TOLOWER "${_cpu_flags}" _cpu_flags)
    string(REPLACE "." "_" _cpu_flags "${_cpu_flags}")
    message (STATUS "WARNING: UNCHARTED BUILD SYSTEM TERRITORY FROM NOW ON")
    # some info by Matthias Moeller, copied from the FEAT2 guess_id() script:
    # #from /usr/include/mach/machine.h
    #         cputype=`sysctl hw.cputype | awk '{print $2 ; exit;}'`
    #         cpusubt=`sysctl hw.cpusubtype | awk '{print $2; exit;}'`
    #         cpufamily=`sysctl hw.cpufamily | awk '{print $2; exit;}'`
    #         case "${cputype}" in
    #             7) cpu=x86
    #                 case "${cpufamily}" in
    #                     1943433984) cpu="coresolo";; # Intel Core Solo / Core Duo (32-bit Pentim-M with SSE3) (Yonah)
    #                     1114597871) cpu="coreduo";;  # Intel Core Duo (Merom)
    #                      2028621756) cpu="penryn" ;;  # Intel Penryn
    #                     1801080018) cpu="nehalem" ;; # Intel Nehalem
    #                     *)
    #                         case "${cpusubt}" in
    #                             3) cpu="i386" ;;
    #                             4) cpu="i486" ;;
    #                             5) cpu="pentium" ;;
    #                             22) cpu="pentiumpro" ;;
    #                             54|86) cpu="pentium2" ;;
    #                             103|119) cpu="celeron" ;;
    #                             8|24|40) cpu="pentium3" ;;
    #                             9) cpu="pentiumm" ;;
    #                             10|26) cpu="pentium4" ;;
    #                             11) cpu="itanium" ;;
    #                             27) cpu="itanium2" ;;
    #                             12) cpu="xeon" ;;
    #                             28) cpu="xeon2" ;;
    #                         esac
    #                         ;;
    #                 esac
    #                 ;;

  elseif (CMAKE_SYSTEM_NAME STREQUAL "Windows")

    # on Windows, CPU information is buried in the system registry database
    get_filename_component(_vendor_id "[HKEY_LOCAL_MACHINE\\Hardware\\Description\\System\\CentralProcessor\\0;VendorIdentifier]" NAME CACHE)
    get_filename_component(_cpu_id "[HKEY_LOCAL_MACHINE\\Hardware\\Description\\System\\CentralProcessor\\0;Identifier]" NAME CACHE)
    mark_as_advanced(_vendor_id _cpu_id)
    string(REGEX REPLACE ".* Family ([0-9]+) .*" "\\1" _cpu_family "${_cpu_id}")
    string(REGEX REPLACE ".* Model ([0-9]+) .*" "\\1" _cpu_model "${_cpu_id}")
    message(AUTHOR_WARNING "If you know how to query the processor capabilities wrt. SSE on Windows, let Dirk know!")
    message (STATUS "WARNING: UNCHARTED BUILD SYSTEM TERRITORY FROM NOW ON")

  else ()
    # unknown arch, so bail out
    # TODO format properly
    message (STATUS "##############################################################")
    message (STATUS "ERROR: The operating system ${CMAKE_SYSTEM_NAME}              ")
    message (STATUS "is not supported.                                             ")
    message (FATAL_ERROR "##############################################################")

  endif (CMAKE_SYSTEM_NAME STREQUAL "Linux")

  # now, map the information we found to a human-readable CPU microarchitecture token
  # first, for Intel processors
  if (_vendor_id STREQUAL "GenuineIntel")
    if (_cpu_family EQUAL 4)
      set (FEAST_CPU_TYPE "i486")
    elseif (_cpu_family EQUAL 5)
      set (FEAST_CPU_TYPE "pentium")
    elseif (_cpu_family EQUAL 6)
      if (_cpu_model LESS 2)
        set(FEAST_CPU_TYPE "pentiumpro")
      elseif (_cpu_model LESS 7)
        set (FEAST_CPU_TYPE "pentium2")
      elseif (_cpu_model LESS 9)
        set (FEAST_CPU_TYPE "pentium3")
      elseif (_cpu_model EQUAL 9)
        set (FEAST_CPU_TYPE "pentiumm")
      elseif (_cpu_model LESS 12)
        set (FEAST_CPU_TYPE "pentium3")
      elseif (_cpu_model EQUAL 13)
        set (FEAST_CPU_TYPE "pentium4m")
      elseif (_cpu_model EQUAL 14)
        set (FEAST_CPU_TYPE "coresolo")
      elseif (_cpu_model EQUAL 15)
        set (FEAST_CPU_TYPE "coreduo")
      elseif (_cpu_model EQUAL 23)
        set (FEAST_CPU_TYPE "penryn")
      elseif (_cpu_model EQUAL 26)
        set (FEAST_CPU_TYPE "nehalem")
      elseif (_cpu_model EQUAL 37)
        set(FEAST_CPU_TYPE "westmere")
      elseif (_cpu_model EQUAL 44)
        set(FEAST_CPU_TYPE "westmere")
      endif (_cpu_model LESS 2)
    elseif (_cpu_family EQUAL 7)
      set (FEAST_CPU_TYPE "itanium")
    elseif (_cpu_family EQUAL 15)
      if (_cpu_model LESS 3)
        set (FEAST_CPU_TYPE "pentium4")
      elseif (_cpu_model EQUAL 3 OR _cpu_model EQUAL 4)
        set (FEAST_CPU_TYPE "nocona")
      endif (_cpu_model LESS 3)
    elseif (_cpu_family EQUAL 32)
      set (FEAST_CPU_TYPE "itanium2")
    endif (_cpu_family EQUAL 4)

  # then, for AMD processors
  elseif (_vendor_id STREQUAL "AuthenticAMD")
    if (_cpu_family EQUAL 4)
      set (FEAST_CPU_TYPE "amd486")
    elseif (_cpu_family EQUAL 5)
      if (_cpu_model LESS 6)
        set (FEAST_CPU_TYPE "k5")
      else ()
        set (FEAST_CPU_TYPE "k6")
      endif (_cpu_model LESS 6)
    elseif (_cpu_family EQUAL 6)
      if (_cpu_model LESS 6)
        set (FEAST_CPU_TYPE "athlon")
      else ()
        set (FEAST_CPU_TYPE "athlonxp")
      endif (_cpu_model LESS 6)
    elseif (_cpu_family EQUAL 15)
      if (_cpu_model LESS 5)
        set (FEAST_CPU_TYPE "athlon64")
      elseif (_cpu_model EQUAL 5)
        set (FEAST_CPU_TYPE "opteron")
      elseif (_cpu_model LESS 37)
        set (FEAST_CPU_TYPE "athlon64")
      elseif (_cpu_model EQUAL 37)
        set (FEAST_CPU_TYPE "opteron")
      elseif (_cpu_model EQUAL 65)
        set (FEAST_CPU_TYPE "opteronx2")
      elseif (_cpu_model EQUAL 72)
        set (FEAST_CPU_TYPE "turionx2")
      endif (_cpu_model LESS 5)
    elseif (_cpu_family EQUAL 16)
      set (FEAST_CPU_TYPE "turionx2")
    endif (_cpu_family EQUAL 4)

  # TODO insert sparc support here once it has been implemented properly

  endif (_vendor_id STREQUAL "GenuineIntel")
endmacro (guess_cpu)





# with these macros defined, we can finally do what we want to do,
# i.e., set the variable FEAST_CPU_TYPE

# case 1: guessing the CPU is requested (because BUILD_ID has been empty
#         and using defaults has been requested
if (GUESS_ARCH)
  set (FEAST_CPU_TYPE "FEAST_CPU_TYPE-NOTFOUND")
  guess_cpu ()
  message (STATUS "##############################################################")
  if (NOT FEAST_CPU_TYPE STREQUAL "FEAST_CPU_TYPE-NOTFOUND")
    message (STATUS "Identified CPU type: ${FEAST_CPU_TYPE}")
  else ()
    message (STATUS "No supported CPU found, using default compiler settings with no arch-aware optimisation.")
    set (FEAST_CPU_TYPE "generic")
  endif ()

else ()

  # shortcut to improve speed
  if (NOT DISPLAY_HELP_ONLY)
    # otherwise, three possible sub-cases: Either the arch token is set, or it is not set,
    # or it is misspelled. There is no way to distinguish between these cases properly,
    # so assume it's there and try to parse it first. If that fails, guess an ID.
    # If parsing succeeds, do not validate the type (by calling guess_cpu again) because
    # this prevents cross-compilation (ok, technically, cross-compilation would still be
    # possible by introducing a "force" token to the build-ID which disables the verification
    # again, but that's an additional level of complexity not considered appropriate by
    # Hilmar and Dominik: if the user sets an architecture token (recall that it's not
    # mandatory in the first place), it is used by force.
    # A possible use-case not involving cross-compilation is that one gets a new CPU,
    # and first wants to test the code with the compiler settings from the generation
    # before (e.g., run nehalem-optimised builds on westmere).
    parse_arch_token ()

    # either it's not given or it's misspelled, so try to guess one
    if (FEAST_CPU_TYPE STREQUAL "FEAST_CPU_TYPE-NOTFOUND")
      message (STATUS "##############################################################")

      guess_cpu ()
      if (NOT FEAST_CPU_TYPE STREQUAL "FEAST_CPU_TYPE-NOTFOUND")
        message (STATUS "Identified CPU type: ${FEAST_CPU_TYPE}")

      else ()
        message (STATUS "No supported CPU found, using default compiler settings with no arch-aware optimisation.")
        set (FEAST_CPU_TYPE "generic")
      endif ()

    endif (FEAST_CPU_TYPE STREQUAL "FEAST_CPU_TYPE-NOTFOUND")
  endif (NOT DISPLAY_HELP_ONLY)
endif (GUESS_ARCH)
