#pragma once
#ifndef SCARC_GUARD_SOLVER_PATTERN_HH
#define SCARC_GUARD_SOLVER_PATTERN_HH 1

#include<kernel/scarc/solver_functor.hpp>
#include<kernel/scarc/solver_data.hpp>

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
               template<typename, typename> class MT_,
               template<typename, typename> class ST_>
      static std::shared_ptr<SolverFunctorBase<VT_<Tag_, DataType_> > > execute(VT_<Tag_, DataType_>& y,
                                                                                SolverDataBase<DataType_, Tag_, VT_, MT_, ST_>& data,
                                                                                Index max_iter)
      {
        std::shared_ptr<SolverFunctorBase<VT_<Tag_, DataType_> > > result(new CompoundSolverFunctor<Algo_, VT_<Tag_, DataType_> >());

        ((CompoundSolverFunctor<Algo_, VT_<Tag_, DataType_> >*)(result.get()))->add_functor(new DefectFunctor<Algo_, VT_<Tag_, DataType_>, MT_<Tag_, DataType_> >(y, data.stored_rhs(), data.stored_sys(), data.stored_sol()));

        return result;
      }

    };

  }
}

#endif
