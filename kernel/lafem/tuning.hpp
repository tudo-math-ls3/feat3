#pragma once
#ifndef KERNEL_LAFEM_TUNING_HPP
#define KERNEL_LAFEM_TUNING_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/util/cuda_util.hpp>

#include <iostream>

#ifdef FEAT_COMPILER_MICROSOFT
// Microsoft's STL implementation is not safe for /Wall "by design", see:
// https://connect.microsoft.com/VisualStudio/feedback/details/1217660/warning-c4265-when-using-functional-header
// We do not disable the warning C4265 globally in the compiler detection header, as
// the warning is kind of useful, so we need to use a few pragmas around the include:
#  pragma warning(push, 0)
#  pragma warning(disable:4265)
#  include <functional>
#  pragma warning(pop)
#else
#  include <functional>
#endif // FEAT_COMPILER_MICROSOFT

namespace FEAT
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
       *
       * Finds optimal cuda runtime parameters (cuda blocksize) for the given matrix and vector.
       * The blocksizes are set globally and are used by all following cuda operations.
       * Thus the given matrix/vector should represent the majority of the following computational effort.
       */
      template <typename DT_, typename IT_>
      static void tune_cuda_blocksize(SparseMatrixELL<Mem::CUDA, DT_, IT_> & matrix)
      {
        DenseVector<Mem::CUDA, DT_, IT_> vector(matrix.columns(), DT_(1));
        DenseVector<Mem::CUDA, DT_, IT_> temp(matrix.rows(), DT_(0));

        /////// SPMV
        double best_toe = std::numeric_limits<double>::max();
        Index best_blocksize = MemoryPool<Mem::CUDA>::blocksize_spmv;

        for (Index blocksize(64) ; blocksize <= 1280 ; blocksize+=64)
        {
          MemoryPool<Mem::CUDA>::blocksize_spmv = blocksize;
          double toe = std::numeric_limits<double>::max();
          try
          {
            auto func = [&] () { Arch::Apply<Mem::CUDA>::ell(temp.elements(), DT_(1), vector.elements(), DT_(0), temp.elements(),
                matrix.val(), matrix.col_ind(), matrix.cs(), matrix.cl(), matrix.C(), matrix.rows()); };
            toe = LAFEM::Tuning::_run_bench(func);
            Util::cuda_check_last_error();
          }
          catch (const FEAT::InternalError & e)
          {
            //did we get at least one valid configuration?
            if (best_toe < std::numeric_limits<double>::max())
              //then break, as we have reached the maximum blocksize already (blocksize loop is counting upwards)
              break;
            else
              //no: we have to increase the blocksize for the first valid configuration
              continue;
          }
          if (toe < best_toe)
          {
            best_toe = toe;
            best_blocksize = blocksize;
          }
        }

        if (best_toe < std::numeric_limits<double>::max())
          MemoryPool<Mem::CUDA>::blocksize_spmv = best_blocksize;
        else
          throw InternalError(__func__, __FILE__, __LINE__, "no valid spmv configuration found!");

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
          catch (const FEAT::InternalError & e)
          {
            //did we get at least one valid configuration?
            if (best_toe < std::numeric_limits<double>::max())
              //then break, as we have reached the maximum blocksize already (blocksize loop is counting upwards)
              break;
            else
              //no: we have to increase the blocksize for the first valid configuration
              continue;
          }
          if (toe < best_toe)
          {
            best_toe = toe;
            best_blocksize = blocksize;
          }
        }

        if (best_toe < std::numeric_limits<double>::max())
          MemoryPool<Mem::CUDA>::blocksize_axpy = best_blocksize;
        else
          throw InternalError(__func__, __FILE__, __LINE__, "no valid axpy configuration found!");
      }
    };
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_TUNING_HPP
