import os
import sys
import shutil
import json
import copy
import platform

__author__ = "Maximilian Esser"
__date__   = "February 2023"

#This script initializes and configures the enviroment for visual studio code.
#The idea is to uses this with the automated configure process of CMakeTools and therefore
#requires the C++Tools and CMakeTools Extensions
#The base idea is, that the required toolchains are detected by the kits component of CMakeTools
#and that this script configures the cmake-variants.yaml file to introduce all necessary CMake options

from build_system.detect_cpu import detect_cpu
from build_system.configure_gcc import configure_gcc
from build_system.configure_icc import configure_icc
from build_system.configure_clang import configure_clang
from build_system.configure_pgi import configure_pgi
from build_system.feat_util import remove_string
from build_system.feat_util import remove_substring
from build_system.feat_util import is_found
from build_system.feat_util import find_exe
from build_system.feat_util import get_output
from build_system.thirdparty_package import *
from build_system.is_subdirectory import *


def add_to_dict(dict, string):
  str_list = string.split(" ")
  for pair in str_list:
    if pair == '':
      continue
    ls = pair.split("=")
    key = ls[0][2:].split(':')[0]
    dict[key] = ls[1]

# TODO help, ? und guess nicht in sys.argv[0] suchen

# output help screen
if len(sys.argv) > 1 and ("help" in " ".join(sys.argv) or "?" in " ".join(sys.argv)):
  print ("Usage: configure_feat [build_id]/[guess]/[help]")
  print ("")
  print ("Configure FEAT for VScode and set up all necessary environment variables.")
  print ("You can execute the script (located in the source directory)  from any folder,")
  print ("that shall contain your object files and executables, and thus become your build directory.")
  print ("")
  print ("Please note that you need to rerun this configure script after you have altered the")
  print ("software or hardware stack of your computer.")
  print ("")
  print ("============= REQUIREMENTS ==================")
  print ("Needs cmake (at least version 2.8) to run properly.")
  print ("")
  print ("============= EXAMPLES ==================")
  print ("")
  print ("Recommended: Let configure_feat choose suitable configuration")
  print ("by inspecting available tools in your $PATH:")
  print ("    %s"%sys.argv[0])
  print ("or")
  print ("    %s guess"%sys.argv[0])
  print
  print ("Use a specific build id:")
  print ("    %s build-id"%sys.argv[0])
  print ("Example:")
  print ("    %s mpi-cuda-umfpack"%sys.argv[0])
  print ("")
  print ("============== BUILD ID ====================")
  print ("The FEAT build-id is a string of various keywords, delimited by hyphes or whitespaces")
  print ("It is used to choose different tools and configurations related to the build process.")
  print ("The keywords, sorted by categories are:")
  print ("  Cluster Parallelisation: mpi")
  #print ("  Compiler frontends: ccache")
  print ("  Backends: cuda, mkl")
  #print ("  Generators: unixmake, ninja, xcode, mingw, msys - defaults to unixmake if none is choosen")
  print ("  Instrumentation: scorep")
  print ("  Debugging: valgrind, cudamemcheck, deathhandler")
  print ("  Precision extensions: quadmath, half, floatx")
  print ("  Compilation process optimisations: cotire")
  print ("  Compiler optimisation steps: lto")
  print ("  Third party source code: alglib, fparser, half, hypre, metis, superlu, umfpack, zfp, zlib, zoltan")
  print ("Note that you must provide exactly one build mode, one parallelisation mode and one compiler.")
  print ("Note that you cannot combine any build id token with the guess token.")
  print ("Note that cudamemcheck is only active if cuda is activated, too.")
  print ("Note that the actual compiler version and mpi implementation is choosen by first come, fist serve")
  print ("corresponding to your $PATH settings.")
#   print ("Note that the source code for third party libraries is automatically downloaded if it is not")
#   print ("provided in the first place.")
  print ("Note that the source code for third party libraries has to be manually downloaded")
  print ("and be available in build or in your $PATH.")
  print ("")
  print ("============== ADDITIONAL OPTIONS ==============")
  print ("--system_host_compiler=/usr/bin/gcc")
  print ("  Selects the system compiler (and header files) used by icc and clang.")
  print ("  Be aware that icc needs the path to a gcc binary and clang needs the path to folder containing the whole toolchain.")
  print ("--cuda_host_compiler=/usr/bin/gcc")
  print ("  Selects the system compiler (and header files) used by nvcc.")
  print ("--cuda_arch=sm_30")
  print ("  Selects the cuda architecture target (defaults to sm_30 if not set).")
  print ("--export_compile_commands")
  print ("  Let CMake export the used compile commands to compile_commands.json.")
  print ("--cpu=cputype")
  print ("  Override cpu auto detection with provided string (see build_system/detect_cpu.py for a list).")
  print ("--cxx=path/to/binary")
  print ("  Override default c++ compiler binary by absolute path to binary. This does not replace the proper compiler statement in the build id.")
  print ("--cc=path/to/binary")
  print ("  Override default c compiler binary by absolute path to binary. This muss 'match' the c++ compiler and the compiler statement in the build id.")
  print ("--cuda_verbose")
  print ("  Be way more verbosive on cuda compilation steps")
  print ("--unroll_banded")
  print ("  Use fancy template meta programming to unroll SparseMatrixBanded kernels.")
  print ("--eickt")
  print ("  Explicitly Instanciate all Common Kernel Templates.")
  print ("--sf_esoteric")
  print ("  Enable uncommon data/index types in solver factory.")
  print ("--ovr_mpi_ops")
  print ("  Override MPI operations with custom implementations.")
  print ("--mpi_thread_multiple.")
  print ("  Use multiple threads for asynchronous mpi calls.")
  print ("--mpi_c=xyz")
  print ("  Set mpi c compiler wrapper to xyz.")
  print ("--mpi_cxx=xyz")
  print ("  Set mpi c++ compiler wrapper to xyz.")
  print ("--mkl_sparse_executor")
  print ("  Use the intel mkl sparse executor interface.")
  print ("--restrict_errors")
  print ("  Abort compilation after an error is encountered.")
  print ("")
  print ("Note that the first gcc(g++) in $PATH is choosen as the cuda host compiler, if no other one is provided.")
  print ("Note that these additional options are only valid in build-id mode and not in guess mode.")
  sys.exit()

#First of all, we set our path variables:
trunk_dirname = os.path.abspath(os.path.dirname(sys.argv[0]))
#print(trunk_dirname)

extra_options = []
# separate extra options from build id
for i in range(1, len(sys.argv)):
  if len(sys.argv[i]) > 2 and sys.argv[i][:2] == "--":
    extra_options.append(sys.argv[i])

# check for cmake binary
if not is_found ("cmake"):
  print("Error: no cmake binary found!")
  sys.exit()

cpu = detect_cpu()

if len(sys.argv) == 1 or "guess" in sys.argv:
  buildid_string = "opt"
  if is_found("mpic++") or is_found("mpicxx"):
    buildid_string += "-mpi"
  if is_found("g++"):
    buildid_string += "-gcc"
  elif is_found("icpc"):
    buildid_string += "-icc"
  elif is_found("clang++"):
    buildid_string += "-clang"
  else:
    print("Error: no comptabile compiler found nor provided!")
    sys.exit()
else:
  # add every parameter to the build id that is no extra_option
  buildid_string = ""
  for i in range(1, len(sys.argv)):
    if not (len(sys.argv[i]) > 2 and sys.argv[i][:2] == "--"):
      buildid_string += sys.argv[i] + "-"

buildid_string = buildid_string.strip("-")
# generate handy build id array of tokens
buildid = buildid_string.split("-")
#erase all optimazation flags
buildid = [id for id in buildid if not id in ["opt", "debug", "fast", "noop"]]
# list of buildid tokens, remove token when used, results in list of unused tokens
unused_tokens = copy.deepcopy(buildid)

# initialize variables
cxxflags = ""
ldflags = ""
cmake_flags = {'BUILD_ID':buildid_string}

debug_cmake = {'FEAT_DEBUG_MODE':'ON'}
opt_cmake = {'FEAT_DEBUG_MODE':'OFF'}

debug_flags= ""
opt_flags = ""
opt_c_flags = ""
debug_c_flags = ""
debug_ldf_flags = ""
opt_ldf_flags = ""

dic_json = dict()
# parallelisation selection
if "mpi" in buildid:
  dic_json['mpi'] = {
    'default':'no-mpi',
    'choices':{
        'no-mpi':{
            'short':'no-mpi',
            'long':'Does not compile with mpi support'
        },
        'mpi':{
            'short':'mpi',
            'long':'Does compile with mpi support',
            'settings':{
                'FEAT_HAVE_MPI':'TRUE'
            }
        }
    }
  }
  remove_string(unused_tokens, "mpi")


# evaluate extra options
cputype = ""
system_host_compiler = ""
cuda_host_compiler = ""
compiler_cxx = ""
compiler_cc = ""
cuda_arch = ""
cuda_verbose = False
restrict_errors = False

unused_extra_options = copy.copy(extra_options)
for option in extra_options:
  if option.startswith("--cpu="):
    unused_extra_options.remove(option)
    if cputype:
      print ("Error: multiple cpu parameters in\n" + "\n".join(extra_options))
      exit(1)
    cputype = option.replace("--cpu=", "", 1)

  if option.startswith("--system_host_compiler="):
    unused_extra_options.remove(option)
    if system_host_compiler:
      print ("Error: multiple system_host_compiler parameters in\n" + "\n".join(extra_options))
      exit(1)
    system_host_compiler = option.replace("--system_host_compiler=", "", 1)
    if not os.path.exists(system_host_compiler):
      print ("Error: system_host_compiler " + system_host_compiler + " does not exist!")
    cmake_flags["FEAT_SYSTEM_HOST_COMPILER"] = system_host_compiler

  if option.startswith("--cuda_host_compiler="):
    unused_extra_options.remove(option)
    if cuda_host_compiler:
      print ("Error: multiple cuda_host_compiler parameters in\n" + "\n".join(extra_options))
      exit(1)
    cuda_host_compiler = option.replace("--cuda_host_compiler=", "", 1)
    if not os.path.exists(cuda_host_compiler):
      print ("Error: cuda_host_compiler " + cuda_host_compiler + " does not exist!")

  if option.startswith("--cxx="):
    unused_extra_options.remove(option)
    if compiler_cxx:
      print ("Error: multiple cxx parameters in\n" + "\n".join(extra_options))
      exit(1)
    compiler_cxx = option.replace("--cxx=", "", 1)
    if not os.path.exists(compiler_cxx):
      print ("Error: compiler " + compiler_cxx + " does not exist!")

  if option.startswith("--cc="):
    unused_extra_options.remove(option)
    if compiler_cc:
      print ("Error: multiple cc parameters in\n" + "\n".join(extra_options))
      exit(1)
    compiler_cc = option.replace("--cc=", "", 1)
    if not os.path.exists(compiler_cc):
      print ("Error: compiler " + compiler_cc + " does not exist!")

  if option.startswith("--cuda_arch="):
    unused_extra_options.remove(option)
    if cuda_arch:
      print ("Error: multiple cuda_arch parameters in\n" + "\n".join(extra_options))
      exit(1)
    cuda_arch = option.replace("--cuda_arch=", "", 1)

  if option.startswith("--cuda_verbose"):
    unused_extra_options.remove(option)
    cuda_verbose = True

  if option.startswith("--export_compile_commands"):
    unused_extra_options.remove(option)
    cmake_flags["CMAKE_EXPORT_COMPILE_COMMANDS"] = 'ON'

  if option.startswith("--unroll_banded"):
    unused_extra_options.remove(option)
    cmake_flags["FEAT_UNROLL_BANDED"]= "TRUE"

  if option.startswith("--eickt"):
    unused_extra_options.remove(option)
    cmake_flags["FEAT_EICKT"] = "TRUE"

  if option.startswith("--sf_esoteric"):
    unused_extra_options.remove(option)
    cmake_flags["FEAT_SF_ESOTERIC"] = "TRUE"

  if option.startswith("--ovr_mpi_ops"):
    unused_extra_options.remove(option)
    cmake_flags["FEAT_OVERRIDE_MPI_OPS"] = "TRUE"

  if option.startswith("--mpi_thread_multiple"):
    unused_extra_options.remove(option)
    cmake_flags["FEAT_MPI_THREAD_MULTIPLE"] = "TRUE"

  if option.startswith("--mpi_c="):
    unused_extra_options.remove(option)
    mpi_c = option.replace("--mpi_c=", "", 1)
    cmake_flags["FEAT_MPI_C"] = mpi_c

  if option.startswith("--mpi_cxx="):
    unused_extra_options.remove(option)
    mpi_cxx = option.replace("--mpi_cxx=", "", 1)
    cmake_flags["FEAT_MPI_CXX"] = mpi_cxx

  if option.startswith("--mkl_sparse_executor"):
    unused_extra_options.remove(option)
    cmake_flags["FEAT_USE_MKL_SPARSE_EXECUTOR"] = "TRUE"

  if option.startswith("--restrict_errors"):
    unused_extra_options.remove(option)
    restrict_errors = True

if "mkl" in buildid:
  remove_string(unused_tokens, "mkl")
  cmake_flags["FEAT_HAVE_MKL"] = "ON"
  cxxflags += " -DMKL_ILP64"

# give cmake on macOS a hint to the openmp lib
if platform.system() == "Darwin":
  ldflags += " -L/usr/local/opt/libomp/lib" #-lomp
  cxxflags += " -I/usr/local/opt/libomp/include"

#before this, the cxx flag has to be finalized, since its used as base for the seperate
#build flags
if not cputype:
  cputype = detect_cpu()
cmake_flags["FEAT_CPU_TYPE"] = cputype

if "gcc" in buildid or "gnu" in buildid or "g++" in buildid:
  if not compiler_cxx:
    compiler_cxx="g++"
  if not is_found(compiler_cxx):
    print ("Error: Chosen cxx binary '" + compiler_cxx +"' not found!")
    sys.exit(1)
  if not compiler_cc:
    compiler_cc="gcc"
  if not is_found(compiler_cc):
    print ("Error: Chosen cc binary '" + compiler_cc +"' not found!")
    sys.exit(1)
  remove_string(unused_tokens, "gcc")
  remove_string(unused_tokens, "gnu")
  os.environ["CC"] = compiler_cc
  os.environ["CXX"] = compiler_cxx
  os.environ["LD"] = compiler_cxx
  cxxflags_temp, cmake_temp = configure_gcc (cputype, ['debug',buildid] , compiler_cxx, restrict_errors)
  debug_flags = cxxflags + cxxflags_temp
  add_to_dict(debug_cmake, cmake_temp)
  cxxflags_temp, cmake_temp = configure_gcc (cputype, ['opt',buildid] , compiler_cxx, restrict_errors)
  opt_flags = cxxflags + cxxflags_temp
  add_to_dict(opt_cmake, cmake_temp)
  if "coverage" in buildid:
    remove_string(unused_tokens, "coverage")
    ldflags += " -fprofile-arcs -ftest-coverage"
  cmake_flags["FEAT_COMPILER_ID"] = "gcc"
elif "clang" in buildid or "llvm" in buildid:
  if not compiler_cxx:
    compiler_cxx = "clang++"
  if not is_found(compiler_cxx):
    print ("Error: Chosen cxx binary '" + compiler_cxx +"' not found!")
    sys.exit(1)
  if not compiler_cc:
    compiler_cc = "clang"
  if not is_found(compiler_cc):
    print ("Error: Chosen cc binary '" + compiler_cc +"' not found!")
    sys.exit(1)
  remove_string(unused_tokens, "clang")
  remove_string(unused_tokens, "llvm")
  os.environ["CC"] = compiler_cc
  os.environ["CXX"] = compiler_cxx
  os.environ["LD"] = compiler_cxx
  #if(platform.system() != "Windows"):
  #  cxxflags += " -pthreads"
  debug_flags = cxxflags + " " + configure_clang (cputype, ['debug',buildid], compiler_cxx, system_host_compiler, restrict_errors)
  opt_flags = cxxflags + " " + configure_clang (cputype, ['opt',buildid], compiler_cxx, system_host_compiler, restrict_errors)
  cmake_flags["FEAT_COMPILER_ID"] = "clang"
elif "pgi" in buildid:
  if not compiler_cxx:
    compiler_cxx = "pgc++"
  if not is_found(compiler_cxx):
    print ("Error: Chosen cxx binary '" + compiler_cxx +"' not found!")
    sys.exit(1)
  if not compiler_cc:
    compiler_cc = "pgcc"
  if not is_found(compiler_cc):
    print ("Error: Chosen cc binary '" + compiler_cc +"' not found!")
    sys.exit(1)
  remove_string(unused_tokens, "pgi")
  os.environ["CC"] = compiler_cc
  os.environ["CXX"] = compiler_cxx
  os.environ["LD"] = compiler_cxx
  debug_flags = cxxflags + " " + configure_pgi (cputype, ['debug',buildid], compiler_cxx, restrict_errors)
  opt_flags = cxxflags + " " + configure_pgi (cputype, ['opt',buildid], compiler_cxx, restrict_errors)
  cmake_flags["FEAT_COMPILER_ID"] = "pgi"
else:
  print ("Error: No supported compiler found in build id:")
  print (buildid_string)
  sys.exit(1)

if "cuda" in buildid:
  if not is_found("nvcc"):
    print ("Error: Choosen backend compiler binary nvcc not found!")
    sys.exit(1)
  remove_string(unused_tokens, "cuda")
  cuda_dict = dict()

  cuda_dict["FEAT_HAVE_CUDA"] = "ON"
  if not cuda_arch:
    cuda_arch = "sm_30"
  cuda_dict["FEAT_CUDA_ARCH"] = cuda_arch

  if cuda_verbose:
    cuda_arch["FEAT_CUDA_VERBOSE"] = "ON"

  if "cudamemcheck" in buildid:
    remove_string(unused_tokens, "cudamemcheck")
    cmake_flags["FEAT_CUDAMEMCHECK"] = "ON"

  if not cuda_host_compiler:
    cuda_host_compiler = find_exe("g++")
  cmake_flags["FEAT_CUDA_HOST_COMPILER"] = cuda_host_compiler
  dic_json['cuda'] = {
    'default':'cuda',
    'choices':{
        'no-cuda':{
            'short':'no-cuda',
            'long':'Does not compile with cuda support',
        },
        'cuda':{
            'short':'cuda',
            'long':'Does compile with cuda support',
            'settings':cuda_dict
        }
    }
  }



if "valgrind" in buildid:
  if not is_found("valgrind"):
    print ("Error: Choosen debugger valgrind not found!")
    sys.exit(1)
  print("valgrind version: " + get_output("valgrind --version")[0])
  remove_string(unused_tokens, "valgrind")
  cmake_flags["FEAT_VALGRIND"] = "ON"

if "quadmath" in buildid:
  remove_string(unused_tokens, "quadmath")
  cmake_flags["FEAT_HAVE_QUADMATH"] = "ON"

if "cotire" in buildid:
  remove_string(unused_tokens, "cotire")
  cmake_flags["FEAT_COTIRE"] = "ON"

# prevent parmetis without mpi enabled (needed by parmetis)
if "mpi" not in buildid and ("parmetis" in buildid or "metis" in buildid):
  print("Error: Using ParMETIS without mpi is not allowed!")
  sys.exit()

# SuperLU needs MPI, because we only use SuperLU_DIST
if ("mpi" not in buildid) and ("superlu" in buildid):
  print("Error: SuperLU can only be used in combination with MPI!")
  sys.exit()

# SuperLU needs MPI, because we only use SuperLU_DIST
if ("mpi" not in buildid) and ("zoltan" in buildid):
  print("Error: Zoltan can only be used in combination with MPI!")
  sys.exit()



# optional third party libraries
for package in available_packages(trunk_dirname+os.sep+"build_system",trunk_dirname+os.sep+"thirdparty"):
  for name in package.names:
    if name in buildid:
      package.add() ##for now the package should be added beforehand
      cmake_flags_tmp = package.cmake_flags
      tmp_dict = dict()
      add_to_dict(tmp_dict, cmake_flags_tmp)
      dic_json[name] = {
        'default':name,
        'choices':{
            "no_"+name:{
                'short':"no_"+name,
                'long':'Does not compile with '+ name +' support',
            },
            name:{
                'short':name,
                'long':'Does compile with ' + name + ' support',
                'settings':tmp_dict
            }
        }
      }
      remove_string(unused_tokens, name)
      continue

#clean up previous config file
clean = [ "CMakeCache.txt", trunk_dirname+os.sep+"CMakeCache.txt" ]
for i in clean:
  if os.path.exists(i):
    os.unlink(i)
clean = [ "feat_config.hpp", trunk_dirname+os.sep+"feat_config.hpp" ]
for i in clean:
  if os.path.exists(i):
    os.unlink(i)
clean = [ "rules.ninja", trunk_dirname+os.sep+"rules.ninja" ]
for i in clean:
  if os.path.exists(i):
    os.unlink(i)
clean = [ "build.ninja", trunk_dirname+os.sep+"build.ninja" ]
for i in clean:
  if os.path.exists(i):
    os.unlink(i)
clean = [ "CMakeFiles", trunk_dirname+os.sep+"CMakeFiles" ]
for i in clean:
  if os.path.exists(i):
    shutil.rmtree(i)

# export compiler/linker flags
os.environ["CXXFLAGS"] = cxxflags
if "lto" in buildid:
  remove_string(unused_tokens, "lto")
  opt_ldf_flags = ldflags + opt_flags
else:
  opt_ldf_flags = ldflags
debug_ldf_flags = ldflags
opt_c_flags = "-O3"
debug_c_flags = "-O3"
# set system host compiler in cflags to pass cmake's c compiler check
if "icc" in buildid or "intel" in buildid or "icpc" in buildid:
  if system_host_compiler:
    opt_c_flags = "-O3  -gcc-name=" + system_host_compiler
    debug_c_flags = "-O3  -gcc-name=" + system_host_compiler

#for now we do not allow to add or retrieve packages through cmake...
#to do add variants for packages
dic_json["cmakeFlags"] = {'default':'generic',
                          'choices': {
                            'generic':{
                              'short':'generic',
                              'long':'Use cmakeFlag as given',
                              'settings': cmake_flags
                            }
                          }
                        }


dic_json['buildType'] = dict()
dic_json['buildType']['default'] = 'opt'
dic_json['buildType'][ 'description'] = 'Build Option'
dic_json['buildType']['choices'] = dict()
dic_json['buildType']['choices']['debug'] = {'short':'Debug',
                                             'long':'Build with debugging info',
                                             'buildType':'Debug'}
dic_json['buildType']['choices']['debug']['settings'] = debug_cmake
dic_json['buildType']['choices']['debug']['env'] = {'CXXFLAGS':debug_flags,
                                                    'LDFLAGS':debug_ldf_flags,
                                                    'CFLAGS':debug_c_flags}
dic_json['buildType']['choices']['opt'] = {'short':'Opt',
                                             'long':'Build without debugging info',
                                             'buildType':'Release'}
dic_json['buildType']['choices']['opt']['settings'] = opt_cmake
dic_json['buildType']['choices']['opt']['env'] = {'CXXFLAGS':opt_flags,
                                                  'LDFLAGS':opt_ldf_flags,
                                                  'CFLAGS':opt_c_flags}

# print out choosen configuration
print ("============== configure_feat ===========")
print ("Build-ID: " + buildid_string)
if len(unused_tokens) == 0:
  #print ("Unused Build-ID tokens: none")
  pass
else:
  print ("")
  print ("Warning: Unused Build-ID tokens: " + "-".join(unused_tokens))
if len(unused_extra_options) == 0:
  #print ("Unused extra option tokens: none")
  pass
else:
  print ("")
  print ("Warning: Unused extra option tokens: " + " ".join(unused_extra_options))
print ("")
print ("Extra Options: " + " ".join(extra_options))
print ("")
print ("CPU Type: " + cputype)
print ("")
print ("cxx: " + os.environ["CXX"])
print ("")
print ("cxxflags_opt: " + opt_flags)
print ("")
print ("cxxflags_debug: " + debug_flags)
print ("")
print ("ldflags: " + opt_ldf_flags)
print ("")
if (cuda_arch != ""):
  print ("cuda_arch: " + cuda_arch)
print ("")
print ("Warning: Reload VScode window before selecting variant!")
print("You can now use CMakeTools to configure your project")

#see if we find mpi in our path, if so, add it as variant
# if(is_found(''))
# print(json.dumps(dic_json, sort_keys=False, indent=4))

print("Writing variants file to " + trunk_dirname + "/.vscode/cmake-variants.json")

with open(trunk_dirname + "/.vscode/cmake-variants.json", 'w') as f:
  f.write(json.dumps(dic_json, sort_keys=False, indent=4))
