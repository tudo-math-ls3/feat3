#pragma once
#ifndef KERNEL_VC10_FEAST_CONFIG_HPP
#define KERNEL_VC10_FEAST_CONFIG_HPP 1

// include user-config header
#include <vc10_user_config.hpp>

/**
 * \file
 * \brief FEAST configuration file for Visual Studio projects.
 *
 * \author Peter Zajac
 */

// Hide the internals of this file from doxygen - it might conflict with other definitions.
/// \cond nodoxy

/* ********************************************************************************************* */
// GENERAL VISUAL STUDIO CONFIGURATIONS
/* ********************************************************************************************* */

// use 'unsigned long long' for 'Index' typedef for 64-bit compilations
//#ifdef _WIN64
//#  define FEAST_INDEX_ULL 1
//#endif

/* ********************************************************************************************* */
// INTEL MKL SUPPORT
/* ********************************************************************************************* */

// check for 32/64 bit MKL path
#ifndef FEAST_VC10_MKL_ROOT_DIR
#  if defined(_WIN64) && defined(FEAST_VC10_MKL_ROOT_DIR_64)
#    define FEAST_VC10_MKL_ROOT_DIR FEAST_VC10_MKL_ROOT_DIR_64
#  elif !defined(_WIN64) && defined(FEAST_VC10_MKL_ROOT_DIR_32)
#    define FEAST_VC10_MKL_ROOT_DIR FEAST_VC10_MKL_ROOT_DIR_32
#  endif
#endif // !defined(FEAST_VC10_MKL_ROOT_DIR)

// check for MKL backend support
#ifdef FEAST_BACKENDS_MKL
#  ifndef FEAST_VC10_MKL_ROOT_DIR
#    error The FEAST_VC10_MKL_ROOT_DIR macro must be defined for MKL support!
#  endif
// define MKL header alias paths
#  ifndef FEAST_VC10_MKL_H_PATH
#    define FEAST_VC10_MKL_H_PATH FEAST_VC10_MKL_ROOT_DIR##/include/mkl.h
#  endif
#  ifndef FEAST_VC10_MKL_SPBLAS_H_PATH
#    define FEAST_VC10_MKL_SPBLAS_H_PATH FEAST_VC10_MKL_ROOT_DIR##/include/mkl_spblas.h
#  endif
#  define FEAST_VC10_MKL_H_PATH_ <FEAST_VC10_MKL_H_PATH>
#  define FEAST_VC10_MKL_SPBLAS_H_PATH_ <FEAST_VC10_MKL_SPBLAS_H_PATH>
// add include-aliases
#  pragma include_alias(<mkl.h>, FEAST_VC10_MKL_H_PATH_)
#  pragma include_alias(<mkl_spblas.h>, FEAST_VC10_MKL_SPBLAS_H_PATH_)
#endif // defined(FEAST_BACKENDS_MKL)

/// \endcond

#endif // KERNEL_VC10_FEAST_CONFIG_HPP
