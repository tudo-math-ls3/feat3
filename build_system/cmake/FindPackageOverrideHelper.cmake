# Helper function for properly overriding find_package calls by FetchContent_MakeAvailable
# This function
# - copies any <lowercaseName>-extra.cmake file for the given package to the CMAKE_FIND_PACKAGE_REDIRECTS_DIR
# - writes a proper config-version file into the CMAKE_FIND_PACKAGE_REDIRECTS_DIR
#
# Arguments:
# - package: Name of the corresponding FetchContent_Declare entry. Must not be lowercase
# - version: The package version as major.minor.patch
# - compat: The version compatability mode
# See the docs of write_basic_package_version_file for more info on version and compat args.
function(find_package_override_helper package version compat)
  include(CMakePackageConfigHelpers)

  string(TOLOWER ${package} package)
  set(extra_file ${PROJECT_SOURCE_DIR}/build_system/cmake/${package}-extra.cmake)
  set(version_file ${CMAKE_FIND_PACKAGE_REDIRECTS_DIR}/${package}-config-version.cmake)

  # Copy <lowercaseName>-extra file to CMAKE_FIND_PACKAGE_REDIRECTS_DIR if it exists
  if(EXISTS ${extra_file})
    file(COPY ${extra_file} DESTINATION ${CMAKE_FIND_PACKAGE_REDIRECTS_DIR})
  endif()

  # Write out version file for package
  write_basic_package_version_file(${version_file} VERSION ${version} COMPATIBILITY ${compat})
endfunction()
