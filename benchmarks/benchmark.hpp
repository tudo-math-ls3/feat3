#pragma once
#ifndef BENCHMARKS_BENCHMARK_HPP
#define BECHNMARKS_BENCHMARK_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/lafem/memory_pool.hpp>

#include <functional>

namespace FEAST
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
      Index iters(1);
      TimeStamp at, bt;
      at.stamp();
      func();
      MemoryPool<Mem_>::synchronize();
      bt.stamp();
      double test_run_time(bt.elapsed(at));
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
        MemoryPool<Mem_>::synchronize();
        bt.stamp();
        times.push_back(bt.elapsed(at));
      }

      double mean(0);
      for (auto & time : times)
        mean += time;
      mean /= double(times.size());
      std::cout<<"TOE: "<<std::fixed<<mean<<std::endl;
      flops *= iters;
      flops /= mean;
      flops /= 1000; // kilo
      flops /= 1000; // mega
      flops /= 1000; // giga
      std::cout<<"GFlop/s: "<<flops<<std::endl;
      bytes *= iters;
      bytes /= mean;
      bytes /= 1024; // kilo
      bytes /= 1024; // mega
      bytes /= 1024; // giga
      std::cout<<"GByte/s: "<<bytes<<std::endl;
      std::cout<<"=============================================="<<std::endl;
    }
  }
}

#endif //BENCHMARKS_BENCHMARK_HPP
