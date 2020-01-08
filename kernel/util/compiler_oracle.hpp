// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_UTIL_COMPILER_ORACLE_HPP
#define KERNEL_UTIL_COMPILER_ORACLE_HPP 1

/**
 * \file compiler_oracle.hpp
 *
 * \brief Compiler detection header for SunStudio/OracleStudio compilers.
 *
 * \author Dominik Goeddeke
 */

#if !defined(FEAT_COMPILER) && defined(__SUNPRO_CC)

// define FEAT_COMPILER_ORACLE macro
#  define FEAT_COMPILER_ORACLE __SUNPRO_CC

// map version to human-readable string and add a few
// C++0x-specific details.
#  if (__SUNPRO_CC == 0x5110)
//   SunStudio does not support nullptr, do not define HAVE_CPP0X_NULLPTR 1
#    define FEAT_COMPILER "SunStudio/OracleStudio C/C++ compiler 12.2"
#  else
#    define FEAT_COMPILER "SunStudio/OracleStudio C/C++ compiler (unknown version)"
#  endif

#endif // !defined(FEAT_COMPILER) && defined(__SUNPRO_CC)

#endif // KERNEL_UTIL_COMPILER_ORACLE_HPP
