#pragma once
#ifndef UTIL_COMPILER_MICROSOFT_HPP
/// Header guard
#define UTIL_COMPILER_MICROSOFT_HPP 1



/**
 * \file compiler_microsoft.hpp
 *
 * \brief Compiler detection header for Microsoft Visual C++ compiler.
 *
 * \author Peter Zajac
 */

#if !defined(FEAST_COMPILER) && defined(_MSC_VER)

// define FEAST_COMPILER_MICROSOFT macro
#  define FEAST_COMPILER_MICROSOFT _MSC_VER

// detect the compiler verson and define the FEAST_COMPILER macro
#  if (_MSC_VER >= 1600)
#    define FEAST_COMPILER "Microsoft Visual C++ .NET 2010"
#  elif (_MSC_VER >= 1500)
#    define FEAST_COMPILER "Microsoft Visual C++ .NET 2008"
#  elif (_MSC_VER >= 1400)
#    define FEAST_COMPILER "Microsoft Visual C++ .NET 2005"
#  elif (_MSC_VER >= 1310)
#    define FEAST_COMPILER "Microsoft Visual C++ .NET 2003"
#  elif (_MSC_VER >= 1300)
#    define FEAST_COMPILER "Microsoft Visual C++ .NET 2002"
#  elif (_MSC_VER >= 1200)
#    define FEAST_COMPILER "Microsoft Visual C++ 6"
#  elif (_MSC_VER >= 1100)
#    define FEAST_COMPILER "Microsoft Visual C++ 5"
#  elif (_MSC_VER >= 1000)
#    define FEAST_COMPILER "Microsoft Visual C++ 4"
#  elif (_MSC_VER >= 900)
#    define FEAST_COMPILER "Microsoft Visual C++ 2"
#  elif (_MSC_VER >= 800)
#    define FEAST_COMPILER "Microsoft Visual C++"
#  else
  // old MS C/C++ compiler, the time before Visual Studio,
  // this compiler version won't be able to compile Feast anyway...
#    define FEAST_COMPILER "Microsoft C/C++ compiler"
#  endif

// Now let's see what C++0x features the compiler offers
#  if (_MSC_VER >= 1600)
#    define HAVE_CPP0X_NULLPTR 1
#    define HAVE_CPP0X_STATIC_ASSERT 1
#  endif

// disable CRT security warnings for standard C/C++ library functions
#define _CRT_SECURE_NO_WARNINGS 1

#endif // !defined(FEAST_COMPILER) && defined(_MSC_VER)

#endif // UTIL_COMPILER_MICROSOFT_HPP
