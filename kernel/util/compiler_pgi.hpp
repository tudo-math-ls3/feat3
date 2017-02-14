#pragma once
#ifndef KERNEL_UTIL_COMPILER_PGI_HPP
#define KERNEL_UTIL_COMPILER_PGI_HPP 1

/**
 * \file compiler_pgi.hpp
 *
 * \brief Compiler detection header for PGI C++ compiler.
 */

#if !defined(FEAT_COMPILER) && defined(__PGI)

// calc linear sortable pgi version
# define _PGI_VER (__PGIC__ * 10000 + __PGIC_MINOR__ * 100 + __PGIC_PATCHLEVEL__)

// define FEAT_COMPILER_PGI macro
# define FEAT_COMPILER_PGI _PGI_VER
# define FEAT_COMPILER "PGI C/C++ compiler"

#endif // !defined(FEAT_COMPILER) && defined(__PGI)

#endif // KERNEL_UTIL_COMPILER_PGI_HPP
