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

    template<typename Pattern_>
    class SolverPatternGeneration
    {
    };

    template<>
    class SolverPatternGeneration<Richardson>
    {
      public:

        /// x_{k+1} <- [u] + x_k, [u] represents a preconditioner applied to a vector
        static std::shared_ptr<FunctorBase> execute(std::shared_ptr<VectorData>& x,
            std::shared_ptr<FunctorBase>& p)
        {
          p = std::shared_ptr<FunctorBase>(new ProxyPreconApply);
          return std::shared_ptr<FunctorBase>(new ProxyVectorSum<ProxyPreconApply, VectorData>(*((ProxyPreconApply*)(p.get())), *(x.get())));
        }

        /// x_{k+1} <- P(b-Ax_k) + x_k
        static std::shared_ptr<FunctorBase> execute(std::shared_ptr<MatrixData>& A,
            std::shared_ptr<VectorData>& x,
            std::shared_ptr<VectorData>& b,
            std::shared_ptr<FunctorBase>& p)
        {
          std::shared_ptr<FunctorBase> defect(new ProxyDefect<VectorData, MatrixData, VectorData>(*(b.get()), *(A.get()), *(x.get())));
          p = std::shared_ptr<FunctorBase>(new ProxyPreconApply(defect));

          return std::shared_ptr<FunctorBase>(new ProxyVectorSum<ProxyPreconApply, VectorData>(*((ProxyPreconApply*)(p.get())), *(x.get())));
        }
    };
  }
}

#endif
