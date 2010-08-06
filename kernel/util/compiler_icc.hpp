#pragma once
#ifndef KERNEL_COMPILER_ICC_HPP
#define KERNEL_COMPILER_ICC_HPP 1

// Compiler detection header for Intel C++ compiler.

/// \todo implement compiler detection for ICC

#if !defined(FEAST_COMPILER) && defined(__INTEL_COMPILER)

// define FEAST_COMPILER_ICC macro
#  define FEAST_COMPILER_ICC 1

#  define FEAST_COMPILER "Intel C/C++ compiler"

#endif // !defined(FEAST_COMPILER) && defined(__INTEL_COMPILER)

#endif // KERNEL_COMPILER_ICC_HPP
