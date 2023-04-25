// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_RUNTIME_HPP
#define KERNEL_RUNTIME_HPP 1

#include <kernel/base_header.hpp>

#include <iostream>

#ifdef FEAT_HAVE_CUDA
#define CUDA_SKELETON_VOID(function_cuda, ...)\
  function_cuda(__VA_ARGS__);\
  return;
#define CUDA_SKELETON_VOID_T1(t1, function_cuda, ...)\
  function_cuda<t1>(__VA_ARGS__);\
  return;
#define CUDA_SKELETON_VOID_T2(t1, t2, function_cuda, ...)\
  function_cuda<t1, t2>(__VA_ARGS__);\
  return;
#else
#define CUDA_SKELETON_VOID(function_cuda, ...)\
  [[fallthrough]];
#define CUDA_SKELETON_VOID_T1(t1, function_cuda, ...)\
  [[fallthrough]];
#define CUDA_SKELETON_VOID_T2(t1, t2, function_cuda, ...)\
  [[fallthrough]];
#endif

#ifdef FEAT_HAVE_MKL
#define MKL_SKELETON_VOID(function_mkl, ...)\
  function_mkl(__VA_ARGS__);\
  return;
#define MKL_SKELETON_VOID_T1(t1, function_mkl, ...)\
  function_mkl<t1>(__VA_ARGS__);\
  return;
#define MKL_SKELETON_VOID_T2(t1, t2, function_mkl, ...)\
  function_mkl<t1, t2 >(__VA_ARGS__);\
  return;
#else
#define MKL_SKELETON_VOID(function_mkl, ...)\
  [[fallthrough]];
#define MKL_SKELETON_VOID_T1(t1, function_mkl, ...)\
  [[fallthrough]];
#define MKL_SKELETON_VOID_T2(t1, t2, function_mkl, ...)\
  [[fallthrough]];
#endif

#define BACKEND_SKELETON_VOID(function_cuda, function_mkl, function_generic, ...)\
  switch(Runtime::get_preferred_backend())\
  {\
    case PreferredBackend::cuda:\
      CUDA_SKELETON_VOID(function_cuda, __VA_ARGS__)\
    case PreferredBackend::mkl:\
      MKL_SKELETON_VOID(function_mkl, __VA_ARGS__)\
    case PreferredBackend::generic:\
    default:\
      function_generic(__VA_ARGS__);\
      return;\
  }

#define BACKEND_SKELETON_VOID_T1(t1, function_cuda, function_mkl, function_generic, ...)\
  switch(Runtime::get_preferred_backend())\
  {\
    case PreferredBackend::cuda:\
      CUDA_SKELETON_VOID_T1(t1, function_cuda, __VA_ARGS__)\
    case PreferredBackend::mkl:\
      MKL_SKELETON_VOID_T1(t1, function_mkl, __VA_ARGS__)\
    case PreferredBackend::generic:\
    default:\
      function_generic<t1>(__VA_ARGS__);\
      return;\
  }

#define BACKEND_SKELETON_VOID_T2(t1, t2, function_cuda, function_mkl, function_generic, ...)\
  switch(Runtime::get_preferred_backend())\
  {\
    case PreferredBackend::cuda:\
      CUDA_SKELETON_VOID_T2(t1, t2, function_cuda, __VA_ARGS__)\
    case PreferredBackend::mkl:\
      MKL_SKELETON_VOID_T2(t1, t2, function_mkl, __VA_ARGS__)\
    case PreferredBackend::generic:\
    default:\
      function_generic<t1, t2>(__VA_ARGS__);\
      return;\
  }

#ifdef FEAT_HAVE_CUDA
#define CUDA_SKELETON_RETURN(function_cuda, ...)\
  return function_cuda(__VA_ARGS__);
#else
#define CUDA_SKELETON_RETURN(function_cuda, ...)\
  [[fallthrough]];
#endif

#ifdef FEAT_HAVE_MKL_
#define MKL_SKELETON_RETURN(function_mkl, ...)\
  return function_mkl(__VA_ARGS__);
#else
#define MKL_SKELETON_RETURN(function_mkl, ...)\
  [[fallthrough]];
#endif

#define BACKEND_SKELETON_RETURN(function_cuda, function_mkl, function_generic, ...)\
  switch(Runtime::get_preferred_backend())\
  {\
    case PreferredBackend::cuda:\
      CUDA_SKELETON_RETURN(function_cuda, __VA_ARGS__)\
    case PreferredBackend::mkl:\
      MKL_SKELETON_RETURN(function_mkl, __VA_ARGS__)\
    case PreferredBackend::generic:\
    default:\
      return function_generic(__VA_ARGS__);\
  }

namespace FEAT
{
  /// The backend that shall be used in all compute heavy calculations
  enum class PreferredBackend
  {
    generic = 0, /**< generic c++ code */
    mkl, /**< intel mkl blas library */
    cuda /**< nvidia cuda gpgpu support */
  };

  std::ostream & operator<< (std::ostream & left, FEAT::PreferredBackend backend);

  /// The class Runtime encapsulates various settings and functionality needed to run FEAT properly
  class Runtime
  {
  private:
    /// signals, if initialize was called
    static bool _initialized;

    /// signals, if finalize was called
    static bool _finished;

    /// the currently preferred backend
    static PreferredBackend _preferred_backend;

  public:
    /**
     * \brief FEAT initialization
     *
     * This function performs the basic initialization of the FEAT library.
     *
     * \attention
     * This function should be the first functional called in an application's
     * \p main function.
     *
     * \param[in] argc, argv
     * The argument parameters of the calling \p main function.
     */
    static void initialize(int& argc, char**& argv, PreferredBackend backend = PreferredBackend::generic);

    /**
     * \brief FEAT abortion
     *
     * This function terminates this process and, in a MPI-based run, also
     * all other processes belonging to this group.
     *
     * \param[in] dump_call_stack
     * Specifies whether to dump the call-stack to stderr prior to process termination.\n
     * Note that a call-stack dump may not be available on all platforms.
     */
    [[noreturn]] static void abort(bool dump_call_stack = true);

    /**
     * \brief FEAT finalization
     *
     * This function finalizes the FEAT library.
     *
     * \attention
     * This function should be the last function called in an application's
     * \p main function.
     *
     * \note
     * Internally this functions calls the MemoryPool::finalize function.
     * To get proper warnings one memory that is still in use and not freed correctly
     * one has to make sure, that any FEAT Container has been destructed when calling the finalize method.
     * This is usually achieved by keeping the C++ main function slim and kicking off all the fancy application stuff in a separate function/method.
     * Thus (at most) every FEAT related stuff is destructed when this separate function/method ends.
     *
     * \returns
     * An exit code (<c>EXIT_SUCCESS</c>) that can be returned by the \p main function.
     */
    static int finalize();

    /// set new preferred backend
    static void set_preferred_backend(PreferredBackend preferred_backend);

    /// get current preferred backend
    static PreferredBackend get_preferred_backend();
  };

  /**
   * Output operator for PreferredBackend
   */
  std::ostream & operator<< (std::ostream & left, PreferredBackend value);

} // namespace FEAT

#endif // KERNEL_RUNTIME_HPP
