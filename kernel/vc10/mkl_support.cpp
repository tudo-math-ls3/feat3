#include <kernel/base_header.hpp>
#ifdef FEAST_BACKENDS_MKL

/**
 * \file MKL support file for Visual Studio
 *
 * \author Peter Zajac
 */

/// \cond internal

// ensure that we compile the hand-made Visual Studio project files
#ifndef VISUAL_STUDIO
#  error This source file must not be compiled by the CMake build system.
#endif

// define STRINGISE preproc abomination
#define STRINGISE_(x) #x
#define STRINGISE(x) STRINGISE_(x)

// define MKL LIB-path
#ifndef FEAST_VC10_MKL_LIB_DIR
#  ifdef _WIN64
#    define FEAST_VC10_MKL_LIB_DIR FEAST_VC10_MKL_ROOT_DIR##/lib/intel64/
#  else
#    define FEAST_VC10_MKL_LIB_DIR FEAST_VC10_MKL_ROOT_DIR##/lib/ia32/
#  endif
#endif
#define FEAST_VC10_MKL_LIB_DIR_ STRINGISE(FEAST_VC10_MKL_LIB_DIR)

// link against interface library
#ifdef _WIN64
#  pragma comment(lib, FEAST_VC10_MKL_LIB_DIR_ "mkl_intel_lp64.lib")
#else
#  pragma comment(lib, FEAST_VC10_MKL_LIB_DIR_ "mkl_intel_c.lib")
#endif

// link against threading library
#pragma comment(lib, FEAST_VC10_MKL_LIB_DIR_ "mkl_sequential.lib")

// link against core library
#pragma comment(lib, FEAST_VC10_MKL_LIB_DIR_ "mkl_core.lib")

// include LAFEM-MKL source files
#include <kernel/lafem/axpy_mkl.cpp>
#include <kernel/lafem/component_product_mkl.cpp>
#include <kernel/lafem/defect_mkl.cpp>
#include <kernel/lafem/difference_mkl.cpp>
#include <kernel/lafem/dot_product_mkl.cpp>
#include <kernel/lafem/norm_mkl.cpp>
#include <kernel/lafem/product_matvec_mkl.cpp>
#include <kernel/lafem/scale_mkl.cpp>
#include <kernel/lafem/sum_mkl.cpp>

/// \endcond

#endif // defined(FEAST_BACKENDS_MKL)
