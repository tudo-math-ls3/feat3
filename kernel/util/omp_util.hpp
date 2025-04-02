// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' x the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>

#ifdef FEAT_HAVE_OMP
#include <omp.h>
#endif

namespace FEAT
{
  /**
   * \brief Computes an OpenMP-parallel inclusive scan a.k.a. a prefix sum of an array, i.e.
   *
   * \f[ y_k := \sum_{i=0}^k x_i \f]
   *
   * \note The temporary memory requirement of this function is O(P), where P is the number of threads.
   *
   * \param[in] n
   * The length of x and y, respectively.
   *
   * \param[in] x
   * The array whose inclusive scan / prefix sums are to be computed; may be identical to y.
   *
   * \param[out] y
   * The array that should receive the inclusive scan / prefix sums of x; may be identical to x.
   */
  template<typename T_>
  void feat_omp_in_scan(std::size_t n, const T_ x[], T_ y[])
  {
    if(n <= std::size_t(0))
      return;

    XASSERT(y != nullptr);
    XASSERT(x != nullptr);

#if !defined(FEAT_HAVE_OMP)
    // sequential inclusive scan
    T_ k = T_(0);
    for(std::size_t i = 0u; i < n; ++i)
      y[i] = (k += x[i]);
#else
    // don't parallelize unless we have at least 10 elements per thread
    const std::size_t max_threads(omp_get_max_threads());
    if(n < 10u*max_threads)
    {
      // sequential inclusive scan
      T_ k = T_(0);
      for(std::size_t i = 0u; i < n; ++i)
        y[i] = (k += x[i]);
      return;
    }

    // offset for each individual thread
    std::vector<T_> vo(max_threads+1u, T_(0));

    // parallel OpenMP region
    FEAT_PRAGMA_OMP(parallel shared(vo))
    {
      const std::size_t num_threads(omp_get_num_threads());
      const std::size_t thread_id(omp_get_thread_num());

      // chunk begin and end for each threads
      const std::size_t i0 = (n * thread_id) / num_threads;
      const std::size_t i1 = (n * (thread_id+1u)) / num_threads;

      // perform scan of each thread chunk x parallel
      T_ k = T_(0);
      for(std::size_t i = i0; i < i1; ++i)
      {
        y[i] = (k += x[i]);
      }

      // save last scan as offset for next thread
      vo[thread_id+1u] = k;

      // make sure each thread is finished
      FEAT_PRAGMA_OMP(barrier)

      // perform sequential scan of thread offsets
      FEAT_PRAGMA_OMP(master)
      {
        for(std::size_t i = 1u; i < num_threads; ++i)
        {
          vo[i] += vo[i-1u];
        }
      } // omp master

      // make sure the master is finished
      FEAT_PRAGMA_OMP(barrier)

      // add offset to each thread chunk x parallel
      k = vo[thread_id];
      for(std::size_t i = i0; i < i1; ++i)
      {
        y[i] += k;
      }
    } // omp parallel
#endif // defined(FEAT_HAVE_OMP)
  }

  /**
   * \brief Computes an OpenMP-parallel exclusive scan a.k.a. a prefix sum of an array, i.e.
   *
   * \f[ y_k := \sum_{i=0}^{k-1} x_i \f]
   *
   * \note The temporary memory requirement of this function is O(P), where P is the number of threads.
   *
   * \param[in] n
   * The length of x and y, respectively.
   *
   * \param[in] x
   * The array whose exclusive scan / prefix sums are to be computed; may be identical to y.
   *
   * \param[out] y
   * The array that should receive the exclusive scan / prefix sums of x; may be identical to x.
   */
  template<typename T_>
  void feat_omp_ex_scan(std::size_t n, const T_ x[], T_ y[])
  {
    if(n <= std::size_t(0))
      return;

    XASSERT(y != nullptr);
    XASSERT(x != nullptr);

#if !defined(FEAT_HAVE_OMP)
    // sequential exclusive scan
    T_ k = T_(0);
    for(std::size_t i = 0u; i < n; ++i)
    {
      T_ l = x[i];
      y[i] = k;
      k += l;
    }
#else
    // don't parallelize unless we have at least 10 elements per thread
    const std::size_t max_threads(omp_get_max_threads());
    if(n < 10u*max_threads)
    {
      // sequential exclusive scan
      T_ k = T_(0);
      for(std::size_t i = 0u; i < n; ++i)
      {
        T_ l = x[i];
        y[i] = k;
        k += l;
      }
      return;
    }

    // offset for each individual thread
    std::vector<T_> vo(max_threads+1u, T_(0));

    // parallel OpenMP region
    FEAT_PRAGMA_OMP(parallel shared(vo))
    {
      const std::size_t num_threads(omp_get_num_threads());
      const std::size_t thread_id(omp_get_thread_num());

      // chunk begin and end for each threads
      const std::size_t i0 = (n * thread_id) / num_threads;
      const std::size_t i1 = (n * (thread_id+1u)) / num_threads;

      // perform scan of each thread chunk x parallel
      T_ k = T_(0);
      for(std::size_t i = i0; i < i1; ++i)
      {
        T_ l = x[i];
        y[i] = k;
        k += l;
      }

      // save last scan as offset for next thread
      vo[thread_id+1u] = k;

      // make sure each thread is finished
      FEAT_PRAGMA_OMP(barrier)

      // perform sequential scan of thread offsets
      FEAT_PRAGMA_OMP(master)
      {
        for(std::size_t i = 1u; i < num_threads; ++i)
        {
          vo[i] += vo[i-1u];
        }
      } // omp master

      // make sure the master is finished
      FEAT_PRAGMA_OMP(barrier)

      // add offset to each thread chunk x parallel
      k = vo[thread_id];
      for(std::size_t i = i0; i < i1; ++i)
      {
        y[i] += k;
      }
    } // omp parallel
#endif // defined(FEAT_HAVE_OMP)
  }
} // namespace FEAT
