# Helper function for adding a thirdparty library (TPL) to FEAT3.
#
# Signature:
# feat_get_tpl(
#   PACKAGE_NAME <name>
#   VERSION <version>
#   URL <url>
#   URL_HASH <hash>
#   [SOURCE_SUBDIR <dir>]
#   [PATCH_COMMAND_WINDOWS command...]
#   [PATCH_COMMAND_LINUX command...]
#   [CONFIG]
#   [EXCLUDE_FROM_ALL]
# )
#
# Parameters
# - PACKAGE_NAME, name of the TPL. Must be same name as used by find_package.
# - VERSION, minimum version of TPL if found by find_package, downloaded version otherwise
# - URL, URL of TPL archive
# - URL_HASH, hash of TPL archive. See https://cmake.org/cmake/help/latest/command/string.html#hash
# - SOURCE_SUBDIR, directory of root cmake file of dependency
# - PATCH_COMMAND_WINDOWS, patch command to run for windows builds
# - PATCH_COMMAND_LINUX, patch command to run for linux builds
# - CONFIG, if set dependencies will be searched for in config mode
# - EXCLUDE_FROM_ALL, if set targets defined by the TPL will not be added to the all target
#
# Behavior
#
# feat_define_tpl can either use libraries installed on the system or download dependencies itself, if required.
#
# If FEAT_PREFER_EXTERNAL_TPL is true or <name>_DIR is defined
# feat_define_tpl first tries to find the requested TPL using find_package.
# Parameters <name> and <version> are forwarded to the find_package call
#
# If no find_package call is to be done or the find_package call failed,
# feat_define_tpl will declare the dependency using FetchContent_declare.
#
# If a system library was found feat_define_tpl will set <name>_FOUND in the parent scope.
# If a dependency was declared using FetchContent_declare feat_define_tpl will add it to
# a list called MAKE_AVAIL_LIST and update that list in the parent scope.
function(feat_get_tpl)

  cmake_parse_arguments(
    PARSE_ARGV 0
    TPL
    "CONFIG;EXLUDE_FROM_ALL" # Options
    "PACKAGE_NAME;VERSION;URL;URL_HASH;SOURCE_SUBDIR" # Single value keywords
    "PATCH_COMMAND_WINDOWS;PATCH_COMMAND_LINUX" # Multi value keywords
  )
  message(STATUS "------------------------------------------------------")
  message(STATUS "- Getting TPL: ${TPL_PACKAGE_NAME} (Version: ${TPL_VERSION})")
  message(STATUS "------------------------------------------------------")

  if(FEAT_PREFER_EXTERNAL_TPL OR ${TPL_PACKAGE_NAME}_DIR)
    list(APPEND FIND_PACKAGE_ARGS ${TPL_PACKAGE_NAME} ${TPL_VERSION})

    if(TPL_CONFIG)
        list(APPEND FIND_PACKAGE_ARGS CONFIG)
    endif()

    if(FEAT_NO_EXTERNAL_DOWNLOAD)
        list(APPEND FIND_PACKAGE_ARGS REQUIRED)
    endif()

    message(STATUS "Trying to find ${TPL_PACKAGE_NAME}")
    find_package(${FIND_PACKAGE_ARGS})
  endif()

  if(${TPL_PACKAGE_NAME}_FOUND)
    print_package_info(${TPL_PACKAGE_NAME})
  else()

    # NOTE(mmuegge): Implementing caching logic ourselves, because its simpler.
    # In CMake versions <= 3.29 ExternalProject_Add (which is internally used by the FetchContent module)
    # supports a DOWNLOAD_DIR settings[1] and a DOWNLOAD_NO_EXTRACT option[2].
    # These should keep downloaded archives around in a directory chosen by us.
    # I am unsure if those archives would be reused automatically by CMake or if we
    # would still need to explicitly pass those paths to FetchContent_declare.
    #
    # In CMake versions >= 3.30 the FetchContent module ignores the DOWNLOAD_DIR settings
    # as per policy CMP0168[3].
    #
    # [1]: https://cmake.org/cmake/help/latest/module/ExternalProject.html#directory-options
    # [2]: https://cmake.org/cmake/help/latest/module/ExternalProject.html#url
    # [3]: https://cmake.org/cmake/help/latest/policy/CMP0168.html
    if(FEAT_TPL_CACHE_DIRECTORY)
      message(STATUS "Checking cache at ${FEAT_TPL_CACHE_DIRECTORY}")
      cmake_path(GET TPL_URL FILENAME ARCHIVE_NAME)
      set(ARCHIVE_NAME ${TPL_PACKAGE_NAME}-${ARCHIVE_NAME})
      set(ARCHIVE_PATH ${FEAT_TPL_CACHE_DIRECTORY}/${ARCHIVE_NAME})

      if(EXISTS ${ARCHIVE_PATH})
        file(MD5 ${ARCHIVE_PATH} ARCHIVE_HASH)
      endif()

      if(NOT EXISTS ${ARCHIVE_PATH} OR NOT "MD5=${ARCHIVE_HASH}" STREQUAL ${TPL_URL_HASH})
        message(STATUS "Archive ${ARCHIVE_NAME} not in cache or hash does not match. Downloading...")
        file(DOWNLOAD ${TPL_URL} ${ARCHIVE_PATH} SHOW_PROGRESS EXPECTED_HASH ${TPL_URL_HASH} STATUS DOWNLOAD_STATUS)

        list(GET DOWNLOAD_STATUS 0 STATUS_CODE)
        list(GET DOWNLOAD_STATUS 1 ERROR_MESSAGE)

        if(NOT ${STATUS_CODE} EQUAL 0)
          message(FATAL_ERROR "Error occured when downloading TPL ${TPL_PACKAGE_NAME}: ${ERROR_MESSAGE} (Status Code: ${STATUS_CODE})")
        endif()
      else()
        message(STATUS "Found archive ${ARCHIVE_NAME} in cache. Reusing archive")
      endif()

      set(TPL_URL ${ARCHIVE_PATH})
    else()
      message(STATUS "Fetching TPL from: ${TPL_URL}")
    endif()


    message(STATUS "Adding ${TPL_PACKAGE_NAME} to build")
    list(APPEND DECLARE_ARGS
      ${TPL_PACKAGE_NAME}
      URL ${TPL_URL}
      URL_HASH ${TPL_URL_HASH}
      UPDATE_DISCONNECTED ON
      OVERRIDE_FIND_PACKAGE
    )

    if(WIN32 AND TPL_PATCH_COMMAND_WINDOWS)
      list(JOIN TPL_PATCH_COMMAND_WINDOWS " " TPL_PATCH_COMMAND_WINDOWS_STRING)
      message(STATUS "Patching TPL with command: " ${TPL_PATCH_COMMAND_WINDOWS_STRING})
      list(APPEND DECLARE_ARGS PATCH_COMMAND ${TPL_PATCH_COMMAND_WINDOWS})
    endif()

    if(NOT WIN32 AND TPL_PATCH_COMMAND_LINUX)
      list(JOIN TPL_PATCH_COMMAND_LINUX " " TPL_PATCH_COMMAND_LINUX_STRING)
      message(STATUS "Patching TPL with command: " ${TPL_PATCH_COMMAND_LINUX_STRING})
      list(APPEND DECLARE_ARGS PATCH_COMMAND ${TPL_PATCH_COMMAND_LINUX})
    endif()

    if(TPL_SOURCE_SUBDIR)
      list(APPEND DECLARE_ARGS SOURCE_SUBDIR ${TPL_SOURCE_SUBDIR})
    endif()

    if(TPL_EXCLUDE_FROM_ALL)
      list(APPEND DECLARE_ARGS EXCLUDE_FROM_ALL)
    endif()

    FetchContent_Declare(${DECLARE_ARGS})

    list(APPEND MAKE_AVAIL_LIST ${TPL_PACKAGE_NAME})
    find_package_override_helper(${TPL_PACKAGE_NAME} ${TPL_VERSION} AnyNewerVersion)
  endif()

  return(PROPAGATE ${TPL_PACKAGE_NAME}_FOUND MAKE_AVAIL_LIST)
endfunction()
