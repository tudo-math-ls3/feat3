// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2024 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_UTIL_COMPILER_CRAY_HPP
#define KERNEL_UTIL_COMPILER_CRAY_HPP 1

/**
 * \file compiler_cray.hpp
 *
 * \brief Compiler detection header for Cray C++ compiler.
 *
 * \author Peter Zajac
 */

#if !defined(FEAT_COMPILER) && defined(_CRAYC)

// calc linear sortable Cray version
#define FEAT_COMPILER_CRAY (_RELEASE_MAJOR * 10000 + _RELEASE_MINOR)

#define FEAT_COMPILER "Cray C++ Compiler " // _RELEASE_STRING

#endif // !defined(FEAT_COMPILER) && defined(_CRAYC)

#endif // KERNEL_UTIL_COMPILER_CRAY_HPP
