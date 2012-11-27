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

    class RichardsonProxy
    {
    };

    template<typename Pattern_, typename Algo_>
    struct SolverPatternGeneration
    {
    };

    template<typename Algo_>
    struct SolverPatternGeneration<Richardson, Algo_>
    {
      template<typename Tag_,
               typename DataType_,
               template<typename, typename> class VT_,
               template<typename, typename> class MT_>
      static std::shared_ptr<SolverFunctorBase<VT_<Tag_, DataType_> > > execute(VT_<Tag_, DataType_>& y,
                                                                                MT_<Tag_, DataType_>& A,
                                                                                VT_<Tag_, DataType_>& x,
                                                                                VT_<Tag_, DataType_>& b,
                                                                                DataType_& norm_init_storage,
                                                                                DataType_& norm_storage,
                                                                                Index& num_iter_storage,
                                                                                Index max_iter)
      {
        return std::shared_ptr<SolverFunctorBase<VT_<Tag_, DataType_> > >(new DefectFunctor<Algo_, VT_<Tag_, DataType_>, MT_<Tag_, DataType_> >(y, b, A, x));
      }

    };

  }
}

#endif
