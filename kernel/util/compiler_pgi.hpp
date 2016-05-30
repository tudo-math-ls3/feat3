#pragma once
#ifndef KERNEL_UTIL_COMPILER_PGI_HPP
#define KERNEL_UTIL_COMPILER_PGI_HPP 1

/**
 * \file compiler_pgi.hpp
 *
 * \brief Compiler detection header for PGI C++ compiler.
 *
 * \author Dominik Goeddeke
 */

#if !defined(FEAT_COMPILER) && defined(__PGI)

// PGI compiler does not make its version number available as a
// preprocessor macro, unfortunately

// define FEAT_COMPILER_PGI macro
# define FEAT_COMPILER_PGI 1
# define FEAT_COMPILER "PGI C/C++ compiler"

// PGI compiler does not support "nullptr", do not define HAVE_CPP0X_NULLPTR 1

#endif // !defined(FEAT_COMPILER) && defined(__PGI)

#endif // KERNEL_UTIL_COMPILER_PGI_HPP
