
#!/usr/bin/env python
########################################################################################################################
# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.
########################################################################################################################
# FEAT3 Visual Studio 2022 project file generator
# -----------------------------------------------
# This script can be used to generate FEAT application and test project files for Microsoft (R) Visual Studio 2022.
#
# Important:
# This tool MUST be executed from the FEAT root directory, because it will not be able to resolve paths and references
# otherwise.
#
# USAGE: vc17-gen <project-path> <project-name> <files...>
#
#   or
#
# USAGE: vc17-gen -app <app-source-files...>
#
#   or
#
# USAGE: vc17-gen -test <test-source-files...>
#
# If the '-app' option is given, then all following arguments are interpreted as application source
# or header files. This script will generate a corresponding application project file in the same
# directory and with the same name as the first passed application source file.
#
# If the '-test' option is given, then all following arguments are interpreted as test source files.
# This script will generate a corresponding test project file in the 'testing' subdirectory of the
# FEAT root directory with the same path and title as the first passed test source file.
#
# \author Peter Zajac
########################################################################################################################

import sys
import os
import re
import uuid

sys.dont_write_bytecode = True

########################################################################################################################

kernel_project_path = os.path.join(".","build_system","vc17","kernel.vc17.vcxproj")
build_modes_path = os.path.join(".","build_system","vc17","build-modes.xml")

build_ids = []
kernel_guid = None
project_guid = "{" + str(uuid.uuid4()).upper() + "}"
project_name = None
project_path = None
root_path = None

test_mode = False
app_mode = False

cpp_list = []
hpp_list = []
cu_list = []

########################################################################################################################
# MAIN
########################################################################################################################

# no arguments?
args = sys.argv[1:]
if len(args) < 1:
  #      123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
  print("This script can be used to generate FEAT application and test project files for use with")
  print("Microsoft (R) Visual Studio 2022.")
  print("")
  print("Important: This tool MUST be executed from the FEAT root directory, because")
  print("           it will not be able to resolve paths and references otherwise.")
  print("")
  print("USAGE: vc17-gen <project-path> <project-name> <files...>")
  print("\n  or\n")
  print("USAGE: vc17-gen -app <app-source-files...>")
  print("\n  or\n")
  print("USAGE: vc17-gen -test <test-source-files...>")
  print("")
  print("If the '-app' option is given, then all following arguments are interpreted as application source")
  print("or header files. This script will generate a corresponding application project file in the same")
  print("directory and with the same name as the first passed application source file.")
  print("")
  print("If the '-test' option is given, then all following arguments are interpreted as test source files.")
  print("This script will generate a corresponding test project file in the 'testing' subdirectory of the")
  print("FEAT root directory with the same path and title as the first passed test source file.")
  sys.exit(0)

# check kernel project path
if not os.path.exists(kernel_project_path):
  print("ERROR: could not find kernel project file")
  print("Note: You must execute this script from the FEAT root directory!")
  sys.exit(1)

# read kernel guid
ifs = open(kernel_project_path, "rt")
for l in ifs.readlines():
  l = l.strip()
  r = re.search("<ProjectGuid>(.+)</ProjectGuid>", l)
  if r:
    kernel_guid = r.group(1)
    break
ifs.close()

# sanity check
if kernel_guid == None:
  print("ERROR: Failed to query kernel project GUID")
  sys.exit(1)

# read build modes
ifs = open(build_modes_path, "rt")
for l in ifs.readlines():
  l = l.strip()
  r = re.search(r"'\$\(Configuration\)'=='(.+)'", l)
  if r:
    build_ids = build_ids + [r.group(1)]
ifs.close()

# sanity check
if len(build_ids) < 2: # at least 'dbg' and 'opt'
  print("ERROR: Failed to query build modes")
  sys.exit(1)

print("Kernel  GUID:   " + kernel_guid);
print("Project GUID:   " + project_guid);

# app-mode ?
if (args[0] == "-app") or (args[0] == "--app"):
  app_mode = True
  args = args[1:]
# test-mode ?
elif (args[0] == "-test") or (args[0] == "--test"):
  test_mode = True
  args = args[1:]

# invalid argument?
if args[0].startswith('-'):
  print("ERROR: invalid option: " + args[0])
  sys.exit(1)

# app mode, test mode or normal mode?
if app_mode:
  # invalid argument count?
  if len(args) < 1:
    print("ERROR: invalid number of arguments")
    sys.exit(1)
  # split path of first source
  project_path, project_name = os.path.split(os.path.splitext(args[0])[0])
  # interpret all arguments as source or header files
  for a in args:
    if a.endswith(".cpp"):
      cpp_list = cpp_list + [os.path.relpath(a, project_path)]
    elif a.endswith(".hpp"):
      hpp_list = hpp_list + [os.path.relpath(a, project_path)]
    elif a.endswith(".cu"):
      cu_list = cu_list + [os.path.relpath(a, project_path)]
    else:
      print("ERROR: unknown file type for file '" + a + "'")
      sys.exit(1)
elif test_mode:
  # invalid argument count?
  if len(args) < 1:
    print("ERROR: invalid number of arguments")
    sys.exit(1)
  # split path of first source
  ppth, project_name = os.path.split(os.path.splitext(args[0])[0])
  # combine project path
  project_path = os.path.join(".","testing",os.path.relpath(ppth, os.path.join(".","kernel")))
  # interpret all arguments as test source files
  for a in args:
    cpp_list = cpp_list + [os.path.relpath(a, project_path)]
else: # normal mode
  # invalid argument count?
  if len(args) < 2:
    print("ERROR: invalid number of arguments")
    sys.exit(1)
  # get project path and name
  project_path = args[0]
  project_name = args[1]
  # interpret all remaining arguments as source or header files
  for a in args[2:]:
    if a.endswith(".cpp"):
      cpp_list = cpp_list + [os.path.relpath(a, project_path)]
    elif a.endswith(".hpp"):
      hpp_list = hpp_list + [os.path.relpath(a, project_path)]
    elif a.endswith(".cu"):
      cu_list = cu_list + [os.path.relpath(a, project_path)]
    else:
      print("ERROR: unknown file type for file '" + a + "'")
      sys.exit(1)

print("Project Path:   " + project_path)
print("Project Name:   " + project_name)
print("")

# define file paths
file_path_sln = os.path.join(project_path, project_name + ".vc17.sln");
file_path_vcx = os.path.join(project_path, project_name + ".vc17.vcxproj");

# generate relative root path
root_path = os.path.relpath(".", project_path)
#print("Root Path: " + root_path)

# create project path
os.makedirs(project_path, exist_ok=True)

# check if the path exists
if not os.path.isdir(project_path):
  print("ERROR: project path does not exist")
  sys.exit(1)
if os.path.exists(file_path_sln):
  print("ERROR: solution file already exists: " + file_path_sln)
  sys.exit(1)
if os.path.exists(file_path_vcx):
  print("ERROR: project file already exists: " + file_path_vcx)
  sys.exit(1)

########################################################################################################################

# write solution file
print("Writing solution file: '" + file_path_sln + "'...")
ofs = open(file_path_sln, "wt")

# write UTF-8 bom
ofs.write("\xEF\xBB\xBF\n")

# write header
ofs.write("Microsoft Visual Studio Solution File, Format Version 12.00\n")
ofs.write("# Visual Studio Version 17\n")
ofs.write("VisualStudioVersion = 17.3.32929.385\n")
ofs.write("MinimumVisualStudioVersion = 10.0.40219.1\n")

# include app project
ofs.write("Project(\"{8BC9CEB8-8B4A-11D0-8D11-00A0C91BC942}\") = \"" + project_name + ".vc17\",")
ofs.write("\"" + project_name + ".vc17.vcxproj\", \"" + project_guid + "\"\nEndProject\n")

# include kernel project
ofs.write("Project(\"{8BC9CEB8-8B4A-11D0-8D11-00A0C91BC942}\") = \"kernel.vc17\",")
ofs.write("\"" + os.path.relpath(kernel_project_path, project_path) + "\",\"" + kernel_guid + "\"\nEndProject\n")

# write solution configurations
ofs.write("Global\n")
ofs.write("\tGlobalSection(SolutionConfigurationPlatforms) = preSolution\n")
for bid in build_ids:
  ofs.write("\t\t" + bid + "|x64 = " + bid + "|x64\n")
ofs.write("\tEndGlobalSection\n")

# write project configurations
ofs.write("\tGlobalSection(ProjectConfigurationPlatforms) = postSolution\n")
for bid in build_ids:
  ofs.write("\t\t" + project_guid + "." + bid + "|x64.ActiveCfg = " + bid + "|x64\n")
  ofs.write("\t\t" + project_guid + "." + bid + "|x64.Build.0 = "   + bid + "|x64\n")
for bid in build_ids:
  ofs.write("\t\t" + kernel_guid + "." + bid + "|x64.ActiveCfg = " + bid + "|x64\n")
  ofs.write("\t\t" + kernel_guid + "." + bid + "|x64.Build.0 = "   + bid + "|x64\n")
ofs.write("\tEndGlobalSection\n")

# write solution properties
ofs.write("\tGlobalSection(SolutionProperties) = preSolution\n")
ofs.write("\t\tHideSolutionNode = FALSE\n")
ofs.write("\tEndGlobalSection\n")
ofs.write("EndGlobal\n")
ofs.close()

########################################################################################################################

# write project file
print("Writing project file: '" + file_path_vcx + "'...")
ofs = open(file_path_vcx, "wt")

# write UTF-8 bom
ofs.write("\xEF\xBB\xBF")

# write header
ofs.write("<?xml version=\"1.0\" encoding=\"utf-8\"?>\n")
ofs.write("<Project DefaultTargets=\"Build\" xmlns=\"http://schemas.microsoft.com/developer/msbuild/2003\">\n")

# write global project properties
ofs.write("  <!-- global project properties -->\n")
ofs.write("  <PropertyGroup Label=\"Globals\">\n")
ofs.write("    <ProjectGuid>" + project_guid + "</ProjectGuid>\n")
ofs.write("    <FeatAppName>" + project_name + "</FeatAppName>\n")
#ofs.write("    <FeatRootPath>$(MSBuildProjectDirectory)\\" + root_path + "</FeatRootPath>\n")
ofs.write("    <FeatRootPath>" + os.path.abspath(".") + "</FeatRootPath>\n")
ofs.write("  </PropertyGroup>\n")

# write common config import
ofs.write("  <!-- import common config -->\n")
ofs.write("  <Import Project=\"$(FeatRootPath)\\build_system\\vc17\\common-config.xml\" />\n")

# write header inclusion list
ofs.write("  <!-- ********************************************************************* -->\n")
ofs.write("  <!-- Header File List -->\n")
ofs.write("  <!-- ********************************************************************* -->\n")
ofs.write("  <ItemGroup Label=\"Header-Files\">\n")
if test_mode:
  ofs.write("    <ClInclude Include=\"" + root_path + "\\test_system\\test_system.hpp\" />\n")
for fn in hpp_list:
  print("Adding HPP header file: '%s'" % fn)
  ofs.write("    <ClInclude Include=\"" + fn + "\" />\n")
ofs.write("  </ItemGroup>\n")

# write source inclusion list
ofs.write("  <!-- ********************************************************************* -->\n")
ofs.write("  <!-- Source File List -->\n")
ofs.write("  <!-- ********************************************************************* -->\n")
ofs.write("  <ItemGroup Label=\"Source-Files\">\n")
if test_mode:
  ofs.write("    <ClCompile Include=\"" + root_path + "\\test_system\\test_system.cpp\" />\n")
for fn in cpp_list:
  print("Adding CPP source file: '%s'" % fn)
  ofs.write("    <ClCompile Include=\"" + fn + "\" />\n")
ofs.write("  </ItemGroup>\n")

# write CUDA inclusion list
ofs.write("  <!-- ********************************************************************* -->\n")
ofs.write("  <!-- CUDA File List -->\n")
ofs.write("  <!-- ********************************************************************* -->\n")
ofs.write("  <ItemGroup Label=\"CUDA-Files\">\n")
for fn in cu_list:
  print("Adding CUDA source file: '%s'" % fn)
  ofs.write("    <CudaCompile Include=\"" + fn + "\" />\n")
ofs.write("  </ItemGroup>\n")

# write kernel project reference
ofs.write("  <!-- ********************************************************************* -->\n")
ofs.write("  <!-- Final Imports -->\n")
ofs.write("  <!-- ********************************************************************* -->\n")
ofs.write("  <ItemGroup>\n")
ofs.write("    <ProjectReference Include=\"$(FeatRootPath)\\build_system\\vc17\\kernel.vc17.vcxproj\">\n")
ofs.write("      <Project>" + kernel_guid + "</Project>\n")
ofs.write("    </ProjectReference>\n")
ofs.write("  </ItemGroup>\n")

# write app-target import
if test_mode:
  ofs.write("  <Import Project=\"$(FeatRootPath)\\build_system\\vc17\\target-test.xml\" />\n")
else:
  ofs.write("  <Import Project=\"$(FeatRootPath)\\build_system\\vc17\\target-app.xml\" />\n")

# end-of-file
ofs.write("</Project>\n")
ofs.close();

########################################################################################################################

print("\nFinished!")
