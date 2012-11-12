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
          return std::shared_ptr<FunctorBase>(new ProxyVectorSum<ProxyPreconApply, VectorData>(*(reinterpret_cast<std::shared_ptr<ProxyPreconApply>* >(&p)), x));
        }

        /// x_{k+1} <- P(b-Ax_k) + x_k
        static std::shared_ptr<FunctorBase> execute(std::shared_ptr<MatrixData>& A,
            std::shared_ptr<VectorData>& x,
            std::shared_ptr<VectorData>& b,
            std::shared_ptr<MatrixData>& P,
            std::shared_ptr<FunctorBase>& p)
        {
          ///for convenience: add a pseudo expression to p although it is covered by the Richardson functor
          std::shared_ptr<FunctorBase> defect(new ProxyDefect<VectorData, MatrixData, VectorData>(b, A, x));
          std::shared_ptr<FunctorBase> rich(new ProxyRichardson<MatrixData, VectorData, VectorData, MatrixData>(A, x, b, P));

          p = std::shared_ptr<FunctorBase>(new ProxyMatrixVectorProduct<MatrixData, ProxyDefect<VectorData, MatrixData, VectorData> >(P, (*reinterpret_cast<std::shared_ptr<ProxyDefect<VectorData, MatrixData, VectorData> >* >(&defect))));

          return std::shared_ptr<FunctorBase>(new ProxyPreconApply(rich));

          //return std::shared_ptr<FunctorBase>(new ProxyVectorSum<ProxyPreconApply, VectorData>(*(reinterpret_cast<std::shared_ptr<ProxyPreconApply>* >(&p)), x));
        }
    };
  }
}

#endif
