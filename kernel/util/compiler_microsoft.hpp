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
#if !defined(FEAST_COMPILER) && defined(_MSC_VER)

// define FEAST_COMPILER_MICROSOFT macro
#  define FEAST_COMPILER_MICROSOFT _MSC_VER

// detect the compiler verson and define the FEAST_COMPILER macro
#  if (_MSC_VER >= 1600)
#    define FEAST_COMPILER "Microsoft Visual C++ .NET 2010 (or newer)"
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

// Now let's see what C++11 features the compiler offers
#  if (_MSC_VER >= 1600)
#    define HAVE_CPP11_NULLPTR 1
#    define HAVE_CPP11_STATIC_ASSERT 1
#  endif

// C4100: 'identifier': unreferenced formal parameter
// This warning arises from unused function parameters. As we are using doxygen, parameters always need
// names to produce useful documentation independent of whether the actual parameter is used or not.
#  pragma warning(disable: 4100)

// C4127: conditional expression is constant
#  pragma warning(disable: 4127)

// C4512: 'class': assignment operator could not be generated
#  pragma warning(disable: 4512)

// C4514: 'function': unreferenced inline function has been removed
// This is an annoying optimisation information.
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

// C4710: 'function': function not inlined
// This is an annoying optimisation information.
#  pragma warning(disable: 4710)

// C4820: 'bytes' bytes padding added after construct 'member_name'
// This warning is disabled as it is mass-produced when compiling standard libraries.
#  pragma warning(disable: 4820)

// C4986: 'function' exception specification does not match previous declaration
// This warning is disabled as it is mass-produced when compiling standard libraries.
#  pragma warning(disable: 4986)

// disable CRT security warnings for standard C/C++ library functions
#  define _CRT_SECURE_NO_WARNINGS 1

#endif // !defined(FEAST_COMPILER) && defined(_MSC_VER)

#endif // KERNEL_UTIL_COMPILER_MICROSOFT_HPP
