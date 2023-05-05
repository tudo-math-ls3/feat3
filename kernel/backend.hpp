// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_BACKEND_HPP
#define KERNEL_BACKEND_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>

// includes, system
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
  switch(Backend::get_preferred_backend())\
  {\
    case FEAT::PreferredBackend::cuda:\
      CUDA_SKELETON_VOID(function_cuda, __VA_ARGS__)\
    case FEAT::PreferredBackend::mkl:\
      MKL_SKELETON_VOID(function_mkl, __VA_ARGS__)\
    case FEAT::PreferredBackend::generic:\
    default:\
      function_generic(__VA_ARGS__);\
      return;\
  }

#define BACKEND_SKELETON_VOID_T1(t1, function_cuda, function_mkl, function_generic, ...)\
  switch(Backend::get_preferred_backend())\
  {\
    case FEAT::PreferredBackend::cuda:\
      CUDA_SKELETON_VOID_T1(t1, function_cuda, __VA_ARGS__)\
    case FEAT::PreferredBackend::mkl:\
      MKL_SKELETON_VOID_T1(t1, function_mkl, __VA_ARGS__)\
    case FEAT::PreferredBackend::generic:\
    default:\
      function_generic<t1>(__VA_ARGS__);\
      return;\
  }

#define BACKEND_SKELETON_VOID_T2(t1, t2, function_cuda, function_mkl, function_generic, ...)\
  switch(Backend::get_preferred_backend())\
  {\
    case FEAT::PreferredBackend::cuda:\
      CUDA_SKELETON_VOID_T2(t1, t2, function_cuda, __VA_ARGS__)\
    case FEAT::PreferredBackend::mkl:\
      MKL_SKELETON_VOID_T2(t1, t2, function_mkl, __VA_ARGS__)\
    case FEAT::PreferredBackend::generic:\
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
  switch(Backend::get_preferred_backend())\
  {\
    case FEAT::PreferredBackend::cuda:\
      CUDA_SKELETON_RETURN(function_cuda, __VA_ARGS__)\
    case FEAT::PreferredBackend::mkl:\
      MKL_SKELETON_RETURN(function_mkl, __VA_ARGS__)\
    case FEAT::PreferredBackend::generic:\
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

  /**
   * \brief Backend support class
   *
   * \author Dirk Ribbrock, Peter Zajac
   */
  class Backend
  {
  private:
    /// the currently preferred backend
    static PreferredBackend _preferred_backend;

  public:
    /// set new preferred backend
    static void set_preferred_backend(PreferredBackend preferred_backend);

    /// get current preferred backend
    static PreferredBackend get_preferred_backend();
  }; // class Backend

  std::ostream & operator<< (std::ostream & left, FEAT::PreferredBackend backend);
} // namespace FEAT

#endif // KERNEL_BACKEND_HPP
