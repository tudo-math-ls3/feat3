#pragma once
#ifndef UTIL_COMPILER_INTEL_HPP
/// Header guard
#define UTIL_COMPILER_INTEL_HPP

// Compiler detection header for Intel C++ compiler.

/// \todo implement compiler detection for Intel Compiler Suite

#if !defined(FEAST_COMPILER) && defined(__INTEL_COMPILER)

// define FEAST_COMPILER_INTEL macro
#  define FEAST_COMPILER_INTEL 1

#  define FEAST_COMPILER "Intel C/C++ compiler"

#endif // !defined(FEAST_COMPILER) && defined(__INTEL_COMPILER)

#endif // UTIL_COMPILER_INTEL_HPP
