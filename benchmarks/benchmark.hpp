// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef BENCHMARKS_BENCHMARK_HPP
#define BECHNMARKS_BENCHMARK_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/kahan_accumulator.hpp>

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
    void run_bench(std::function<void (void)> func, double flops, double bytes)
    {
#ifdef FEAT_DEBUG_MODE
      std::cout << "\n" << String(80, '*') << "\n" << String(80, '*') << "\n\n";
      std::cout << "WARNING: You are running a benchmark in DEBUG mode!\n\n";
      std::cout << "This is generally a really bad idea because all performance metrics will\n";
      std::cout << "be utterly pointless due to missing compile-time optimizations, so only\n";
      std::cout << "do this if you are currently in bug hunting mode.\n\n";
      std::cout << String(80, '*') << "\n" << String(80, '*') << "\n\n";
#endif

      Index iters(1);
      //warmup
      func();
      Runtime::synchronize_devices();

      TimeStamp at;
      func();
      Runtime::synchronize_devices();
      double test_run_time(at.elapsed_now());
      std::cout<<"test time: "<<test_run_time<<"\n";
      if (test_run_time < 0.1)
        iters = Index(0.1 / test_run_time) + 1;
      std::cout<<"iters: "<<iters<<"\n";

      std::vector<double> times;
      for (Index i(0) ; i < 10 ; ++i)
      {
        at.stamp();
        for (Index j(0) ; j < iters ; ++j)
        {
          func();
        }
        Runtime::synchronize_devices();
        times.push_back(at.elapsed_now());
      }

      KahanAccumulator<double> mean_ka;
      for (auto & time : times)
        mean_ka += time;
      mean_ka.value /= double(times.size());
      std::cout<<"TOE: "<<std::fixed<<mean_ka.value<<"; duration of "<< iters << " function calls, average over " << 10 << " runs.\n";
      std::cout<<"TOE per function call: "<<std::fixed<<mean_ka.value/double(iters)<<"\n";
      flops *= double(iters);
      flops /= mean_ka.value;
      flops /= 1000.; // kilo
      flops /= 1000.; // mega
      flops /= 1000.; // giga
      std::cout<<"GFlop/s: "<<flops<<"\n";
      bytes *= double(iters);
      bytes /= mean_ka.value;
      bytes /= 1024.; // kilo
      bytes /= 1024.; // mega
      bytes /= 1024.; // giga
      std::cout<<"GByte/s: "<<bytes<<"\n";
      std::cout << String(80, '=') << "\n";
    }
  }
}

#endif //BENCHMARKS_BENCHMARK_HPP
