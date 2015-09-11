#pragma once
#ifndef KERNEL_LAFEM_TUNING_HPP
#define KERNEL_LAFEM_TUNING_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/util/cuda_util.hpp>

#include <iostream>
#include <functional>

namespace FEAST
{
  namespace LAFEM
  {
    class Tuning
    {
      private:
      static double _run_bench(std::function<void (void)> func)
      {
        Index iters(1);
        //warmup
        func();
        MemoryPool<Mem::CUDA>::synchronize();

        TimeStamp at, bt;
        at.stamp();
        func();
        MemoryPool<Mem::CUDA>::synchronize();
        bt.stamp();
        double test_run_time(bt.elapsed(at));
        if (test_run_time < 0.1)
          iters = Index(0.1 / test_run_time) + 1;

        std::vector<double> times;
        for (Index i(0) ; i < 3 ; ++i)
        {
          at.stamp();
          for (Index j(0) ; j < iters ; ++j)
          {
            func();
          }
          MemoryPool<Mem::CUDA>::synchronize();
          bt.stamp();
          times.push_back(bt.elapsed(at));
        }

        double mean(0);
        for (auto & time : times)
          mean += time;
        mean /= double(times.size());
        mean /= double(iters); // mean per single func call

        return mean;
      }

      public:
      /**
       * \brief Find optimal cuda runtime parameters
       *
       * \param[in] matrix The matrix to use for tuning.
       * \param[in] vector The vector to use for tuning.
       *
       * Finds optimal cuda runtime parameters (cuda blocksize) for the given matrix and vector.
       * The blocksizes are set globally and are used by all following cuda operations.
       * Thus the given matrix/vector should represent the majority of the following computational effort.
       */
      template <typename DT_, typename IT_>
      static void tune_cuda_blocksize(SparseMatrixELL<Mem::CUDA, DT_, IT_> & matrix, DenseVector<Mem::CUDA, DT_, IT_> & vector)
      {
        DenseVector<Mem::CUDA, DT_, IT_> temp(vector.size(), DT_(0));

        /////// SPMV
        double best_toe = std::numeric_limits<double>::max();
        Index best_blocksize = MemoryPool<Mem::CUDA>::blocksize_spmv;

        for (Index blocksize(64) ; blocksize <= 1280 ; blocksize+=64)
        {
          MemoryPool<Mem::CUDA>::blocksize_spmv = blocksize;
          double toe = std::numeric_limits<double>::max();
          try
          {
            auto func = [&] () { Arch::ProductMatVec<Mem::CUDA>::ell(temp.elements(), matrix.val(), matrix.col_ind(), matrix.cs(), matrix.cl(), vector.elements(), matrix.C(), matrix.rows()); };
            toe = LAFEM::Tuning::_run_bench(func);
            Util::cuda_check_last_error();
          }
          catch (const FEAST::InternalError & e)
          {
            continue;
          }
          if (toe < best_toe)
          {
            best_toe = toe;
            best_blocksize = blocksize;
          }
        }

        MemoryPool<Mem::CUDA>::blocksize_spmv = best_blocksize;

        ///// AXPY
        best_toe = std::numeric_limits<double>::max();
        best_blocksize = MemoryPool<Mem::CUDA>::blocksize_axpy;

        for (Index blocksize(64) ; blocksize <= 1280 ; blocksize+=64)
        {
          MemoryPool<Mem::CUDA>::blocksize_axpy = blocksize;
          double toe = std::numeric_limits<double>::max();
          try
          {
            auto func = [&] () { Arch::Axpy<Mem::CUDA>::dv(temp.elements(), DT_(1.234), vector.elements(), temp.elements(), temp.size()); };
            toe = LAFEM::Tuning::_run_bench(func);
            Util::cuda_check_last_error();
          }
          catch (const FEAST::InternalError & e)
          {
            continue;
          }
          if (toe < best_toe)
          {
            best_toe = toe;
            best_blocksize = blocksize;
          }
        }

        MemoryPool<Mem::CUDA>::blocksize_axpy = best_blocksize;
      }
    };
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_TUNING_HPP
