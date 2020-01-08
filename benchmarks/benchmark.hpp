// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef BENCHMARKS_BENCHMARK_HPP
#define BECHNMARKS_BENCHMARK_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/memory_pool.hpp>

#include <functional>

namespace FEAT
{
  /// common benchmarking methods
  namespace Benchmark
  {
    /**
     * Runs a given benchmark and evaluates the timing results
     *
     * \param[in] func The function to evaluate (best given as lambda).
     * \param[in] flops The flop count of a single function call.
     * \param[in] bytes The amount of bytes moved by a single function call.
     *
     **/
    template <typename Mem_>
    void run_bench(std::function<void (void)> func, double flops, double bytes)
    {
#ifdef FEAT_DEBUG_MODE
      throw InternalError("You are running a benchmark in DEBUG mode!");
#endif

      Index iters(1);
      //warmup
      func();
      MemoryPool<Mem_>::synchronise();

      TimeStamp at;
      func();
      MemoryPool<Mem_>::synchronise();
      double test_run_time(at.elapsed_now());
      std::cout<<"test time: "<<test_run_time<<std::endl;
      if (test_run_time < 0.1)
        iters = Index(0.1 / test_run_time) + 1;
      std::cout<<"iters: "<<iters<<std::endl;

      std::vector<double> times;
      for (Index i(0) ; i < 10 ; ++i)
      {
        at.stamp();
        for (Index j(0) ; j < iters ; ++j)
        {
          func();
        }
        MemoryPool<Mem_>::synchronise();
        times.push_back(at.elapsed_now());
      }

      double mean(0);
      /// \todo use KahanSummation from statistics.hpp, see https://stackoverflow.com/questions/10330002/sum-of-small-double-numbers-c/10330857#10330857 for std::accumulate usage
      for (auto & time : times)
        mean += time;
      mean /= double(times.size());
      std::cout<<"TOE: "<<std::fixed<<mean<<"; duration of "<< iters << " function calls, average over " << 10 << " runs."<<std::endl;
      std::cout<<"TOE per function call: "<<std::fixed<<mean/double(iters)<<std::endl;
      flops *= iters;
      flops /= mean;
      flops /= 1000.; // kilo
      flops /= 1000.; // mega
      flops /= 1000.; // giga
      std::cout<<"GFlop/s: "<<flops<<std::endl;
      bytes *= iters;
      bytes /= mean;
      bytes /= 1024.; // kilo
      bytes /= 1024.; // mega
      bytes /= 1024.; // giga
      std::cout<<"GByte/s: "<<bytes<<std::endl;
      std::cout<<"=============================================="<<std::endl;
    }
  }
}

#endif //BENCHMARKS_BENCHMARK_HPP
