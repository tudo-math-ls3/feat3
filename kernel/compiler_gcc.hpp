#pragma once
#ifndef KERNEL_COMPILER_GCC_HPP
#define KERNEL_COMPILER_GCC_HPP 1

// Compiler detection header for GNU C++ compiler.

/// \todo detect implement compiler detection for GCC

#if !defined(FEAST_COMPILER) && defined(__GNUC__)

// define FEAST_COMPILER_GCC
#  define FEAST_COMPILER_GCC 1

#  define FEAST_COMPILER "GNU C++ compiler"

#endif // !defined(FEAST_COMPILER) && defined(__GNUC__)

#endif // KERNEL_COMPILER_GCC_HPP
