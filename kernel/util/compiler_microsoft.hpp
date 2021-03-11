// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_UTIL_COMPILER_MICROSOFT_HPP
#define KERNEL_UTIL_COMPILER_MICROSOFT_HPP 1

/**
 * \file compiler_microsoft.hpp
 *
 * \brief Compiler detection header for Microsoft Visual C++ compiler.
 *
 * \author Peter Zajac
 */
#if !defined(FEAT_COMPILER) && defined(_MSC_VER)

// define FEAT_COMPILER_MICROSOFT macro
#  define FEAT_COMPILER_MICROSOFT _MSC_VER

// detect the compiler verson and define the FEAT_COMPILER macro
#  if (_MSC_VER >= 1920)
#    define FEAT_COMPILER "Microsoft Visual C++ 2019"
#  elif (_MSC_VER >= 1910)
#    define FEAT_COMPILER "Microsoft Visual C++ 2017"
#  elif (_MSC_VER >= 1900)
#    define FEAT_COMPILER "Microsoft Visual C++ 2015"
#  elif (_MSC_VER >= 1800)
#    define FEAT_COMPILER "Microsoft Visual C++ 2013"
#  elif (_MSC_VER >= 1700)
#    define FEAT_COMPILER "Microsoft Visual C++ 2012"
#  elif (_MSC_VER >= 1600)
#    define FEAT_COMPILER "Microsoft Visual C++ 2010"
#  else
  // this compiler version won't be able to compile FEAT anyway...
#    define FEAT_COMPILER "Microsoft C/C++ compiler"
#  endif

#  define FEAT_PRAGMA_IVDEP __pragma(loop(ivdep))

#  define FEAT_DISABLE_WARNINGS __pragma(warning(push, 0))
#  define FEAT_RESTORE_WARNINGS __pragma(warning(pop))

#  define FORCE_INLINE __forceinline

// define the noinline specifier
#  define NOINLINE __declspec(noinline)

// C4061: enumerator 'identifier' in switch of enum 'enumeration' is not explicitly handled by a case label
// This warning is emitted when a 'switch' handles one or more cases using a 'default' block.
// Note: If there are unhandled cases and there is no default block, the compiler emits a C4062 warning.
#  pragma warning(disable: 4061)

// C4127: conditional expression is constant
// This warning arises for instance in an expression like 'if(true)'.
#  pragma warning(disable: 4127)

// C4180: qualifier applied to function type has no meaning; ignored
// This warning arises when a non-pointer return type of a function is declared as 'const'.
#  pragma warning(disable: 4180)

// C4503: 'identifier': decorated name length exceeded, name was truncated
// This warning arises from heavy template nesting, blowing the compiler's limit on maximal name lengths.
// Running into this warning does not affect the correctness of the program, however, it might confuse
// the debugger.
#  pragma warning(disable: 4503)

// C4512: 'class': assignment operator could not be generated
#  pragma warning(disable: 4512)

// C4514: 'function': unreferenced inline function has been removed
// This is an annoying optimization information.
#  pragma warning(disable: 4514)

// C4555: expression has no effect; expected expression with side-effect
#  pragma warning(disable: 4555)

// C4571: Informational: catch(...) semantics changed since Visual C++ 7.1;
//        structured exceptions (SEH) are no longer caught
#  pragma warning(disable: 4571)

// C4625: 'derived class': copy constructor could not be generated because
//        base class copy constructor is inaccessible
// This warning arises from our non-copyable instantiation policy.
#  pragma warning(disable: 4625)

// C4626: 'derived class': assignment operator could not be generated because
//        base class assignment operator is inaccessible
// This warning arises from our non-copyable instantiation policy.
#  pragma warning(disable: 4626)

// C4702: unreachable code
// This warning is emitted in any case where the compiler detects a statement which will never be executed.
// This happens e.g. in for-loops which always perform zero iterations due to the chosen template parameters,
// buzzword 'dead code elimination'.
#  pragma warning(disable: 4702)

// C4710: 'function': function not inlined
// This is an annoying optimization information.
#  pragma warning(disable: 4710)

// C4711: function 'function' selected for inline expansion
// This is an annoying optimization information.
#  pragma warning(disable: 4711)

// C4738: storing 32-bit float result in memory, possible loss of performance
// This is an optimization warning, which arises from the strict fp-model.
#  pragma warning(disable: 4738)

// C4820: 'bytes' bytes padding added after construct 'member_name'
// This warning is disabled as it is mass-produced when compiling standard libraries.
#  pragma warning(disable: 4820)

// C4883: function size suppresses optimizations
// This is an annoying optimization information.
#  pragma warning(disable: 4883)

// C4938: 'var' : Floating point reduction variable may cause inconsistent results
//                under /fp:strict or #pragma fenv_access
// This warning is issued when reducing floating point variables via OpenMP and
// it has been enabled because it only states the obvious effects of multithreading.
#  pragma warning(disable: 4938)

// C5024: 'class' : move constructor was implicitly defined as deleted
// C5025: 'class' : move assignment operator was implicitly defined as deleted
// C5026: 'class' : move constructor was implicitly defined as deleted because
//                  a base class move constructor is inaccessible or deleted
// C5027: 'class' : move assignment operator was implicitly defined as deleted because
//                  a base class move assignment operator is inaccessible or deleted
#  pragma warning(disable: 5024)
#  pragma warning(disable: 5025)
#  pragma warning(disable: 5026)
#  pragma warning(disable: 5027)

// disable CRT security warnings for standard C/C++ library functions
#ifndef _CRT_SECURE_NO_WARNINGS
#  define _CRT_SECURE_NO_WARNINGS 1
#endif

#endif // !defined(FEAT_COMPILER) && defined(_MSC_VER)

#endif // KERNEL_UTIL_COMPILER_MICROSOFT_HPP
