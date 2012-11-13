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

    template<typename Pattern_>
    class SolverPatternGeneration
    {
    };

    template<typename Pattern_>
    class SolverOperatorGeneration
    {
    };

    class ApproximateInverseMultiply
    {
    };

    template<>
    class SolverOperatorGeneration<ApproximateInverseMultiply>
    {
      public:

        template<typename PExpr_, typename XExpr_>
        static std::shared_ptr<ProxyMatrixVectorProduct<PExpr_, XExpr_> > execute(std::shared_ptr<PExpr_>& P, std::shared_ptr<XExpr_>& x)
        {
          return std::shared_ptr<ProxyMatrixVectorProduct<PExpr_, XExpr_> >(new ProxyMatrixVectorProduct<PExpr_, XExpr_>(P, x));
        }
    };

    class Richardson
    {
    };

    template<>
    class SolverPatternGeneration<Richardson>
    {
      public:

        /// x_{k+1} <- [u] + x_k, [u] represents a preconditioner applied to a vector
        template<typename XExpr_>
        static std::shared_ptr<FunctorBase> execute(std::shared_ptr<XExpr_>& x,
            std::shared_ptr<FunctorBase>& p)
        {
          p = std::shared_ptr<FunctorBase>(new ProxyPreconApply);
          return std::shared_ptr<FunctorBase>(new ProxyVectorSum<ProxyPreconApply, XExpr_>(*(reinterpret_cast<std::shared_ptr<ProxyPreconApply>* >(&p)), x));
        }

        /// x_{k+1} <- P(b-Ax_k) + x_k
        template<typename AExpr_,
                 typename XExpr_,
                 typename BExpr_,
                 typename PExpr_>
        static std::shared_ptr<FunctorBase> execute(std::shared_ptr<AExpr_>& A,
                                                    std::shared_ptr<XExpr_>& x,
                                                    std::shared_ptr<BExpr_>& b,
                                                    std::shared_ptr<PExpr_>& P,
                                                    std::shared_ptr<FunctorBase>& p)
        {
          p = *(reinterpret_cast<std::shared_ptr<FunctorBase>* >(&P));
          std::shared_ptr<FunctorBase> rich(new ProxyRichardson<AExpr_, XExpr_, BExpr_, PExpr_>(A, x, b, P));
          return std::shared_ptr<FunctorBase>(new ProxyPreconApply(rich));
        }
    };

  }
}

#endif
