#pragma once
#ifndef SCARC_GUARD_SOLVER_PATTERN_HH
#define SCARC_GUARD_SOLVER_PATTERN_HH 1

#include<kernel/scarc/solver_functor.hpp>

using namespace FEAST::Foundation;
using namespace FEAST;

namespace FEAST
{
  namespace ScaRC
  {
    class Richardson
    {
    };

    class RichardsonLocal
    {
    };

    template<typename Pattern_>
    class SolverPatternGeneration
    {
    };

    template<>
    class SolverPatternGeneration<Richardson>
    {
      public:

        static std::shared_ptr<FunctorBase> execute(std::shared_ptr<VectorData>& x,
            std::shared_ptr<FunctorBase>& p)
        {
          p = std::shared_ptr<FunctorBase>(new ProxyPreconApply);
          return std::shared_ptr<FunctorBase>(new ProxyVectorSum<VectorData, ProxyPreconApply>(*(x.get()), *((ProxyPreconApply*)(p.get()))));
        }
    };

    template<>
      class SolverPatternGeneration<RichardsonLocal>
      {
        public:

          static std::shared_ptr<FunctorBase> execute(std::shared_ptr<MatrixData>& A,
              std::shared_ptr<VectorData>& x,
              std::shared_ptr<VectorData>& b,
              std::shared_ptr<FunctorBase>& p)
          {
            std::shared_ptr<FunctorBase> defect(new ProxyDefect<VectorData, MatrixData, VectorData>(*(b.get()), *(A.get()), *(x.get())));
            p = std::shared_ptr<FunctorBase>(new ProxyPreconApply(defect));

            return std::shared_ptr<FunctorBase>(new ProxyVectorSum<VectorData, ProxyPreconApply>(*(x.get()), *((ProxyPreconApply*)(p.get()))));
          }
      };
  }
}

#endif
