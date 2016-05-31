#pragma once
#ifndef KERNEL_GLOBAL_CONVERT_PRECOND_HPP
#define KERNEL_GLOBAL_CONVERT_PRECOND_HPP 1

// includes, FEAT
#include <kernel/solver/base.hpp>
#include <kernel/global/gate.hpp>

namespace FEAT
{
  namespace Global
  {
    /**
     * \brief Global Solver conversion module.
     *
     * This class implements a simple solver module, which converts data before and after
     * the inner solver call
     *
     * There are two typical applications:
     * - mixed precission iterative refinement
     *   The outer loop calls the ConvertPrecond object, which calls the inner solver,
     *   that is running in lower precission
     *
     * - mixed arch intra solver
     *   For example a v-cycle multigrid can solve the finer levels on the gpu and switch
     *   completly to the cpu for the coarser levels.
     *   This can be easily obtained by passing the ConvertPrecond object as a coarse solver to the
     *   'finer' v cycle and let the let the ConvertPrecond object call the 'coarser' v cycle as its own solver.
     *
     * \author Dirk Ribbrock
     */
    template<typename VectorOuter_, typename VectorInner_>
    class ConvertPrecond :
      public SolverBase<VectorOuter_>
    {
    public:
      using VectorTypeOuter = VectorOuter_;
      using BaseClass = SolverBase<VectorTypeOuter>;
      using GateTypeOuter = Global::Gate<typename VectorTypeOuter::LocalVectorType>;

      using VectorTypeInner = VectorInner_;
      using SolverTypeInner = SolverBase<VectorTypeInner>;
      using GateTypeInner = Global::Gate<typename VectorTypeInner::LocalVectorType>;

    protected:
      std::shared_ptr<SolverTypeInner> _inner_solver;
      GateTypeInner* _gate_inner;
      GateTypeOuter* _gate_outer;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] inner_solver
       * The acutal solver, which shall be executed with converted rhs/sol vectors
       *
       * \param[in] gate_inner
       * The gate corresponding to the inner rhs/sol vectors.
       *
       * \param[in] gate_outer
       * The gate corresponding to the outer rhs/sol vectors.
       */
      explicit ConvertPrecond(std::shared_ptr<SolverTypeInner> inner_solver, GateTypeInner * gate_inner, GateTypeOuter * gate_outer) :
        _inner_solver(inner_solver),
        _gate_inner(gate_inner),
        _gate_outer(gate_outer)
      {
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "Convert";
      }

      virtual void init_branch(String parent = "") override
      {
        BaseClass::init_branch(parent);
        if(_inner_solver)
          _inner_solver->init_branch(parent + "::" + this->name());
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();
        if(_inner_solver)
          _inner_solver->init_symbolic();
      }

      virtual void init_numeric() override
      {
        BaseClass::init_numeric();
        if(_inner_solver)
          _inner_solver->init_numeric();
      }

      virtual void done_numeric() override
      {
        if(_inner_solver)
          _inner_solver->done_numeric();
        BaseClass::done_numeric();
      }

      virtual void done_symbolic() override
      {
        if(_inner_solver)
          _inner_solver->done_symbolic();
        BaseClass::done_symbolic();
      }

      virtual String get_formated_solver_tree() override
      {
        String result;
        result += this->name();
        if(_inner_solver)
        {
          result += " ( ";
          result += _inner_solver->get_formated_solver_tree();
          result += " ) ";
        }
        return result;
      }

      virtual Status apply(VectorTypeOuter& vec_cor, const VectorTypeOuter& vec_def) override
      {
        VectorTypeInner vec_def_inner;
        vec_def_inner.convert(_gate_inner, vec_def);
        VectorTypeInner vec_cor_inner(vec_def_inner.clone(LAFEM::CloneMode::Layout));

        _inner_solver->apply(vec_cor_inner, vec_def_inner);

        vec_cor.convert(_gate_outer, vec_cor_inner);
        return Status::success;
      }
    }; // class ConvertPrecond<...>

    /**
     * \brief Creates a new ConvertPrecond solver object
     *
       * \param[in] inner_solver
       * The acutal solver, which shall be executed with converted rhs/sol vectors
       *
       * \param[in] gate_inner
       * The gate corresponding to the inner rhs/sol vectors.
       *
       * \param[in] gate_outer
       * The gate corresponding to the outer rhs/sol vectors.
       *
       * \returns
       * A shared pointer to a new ConvertPrecond object.
       */
    template<typename VectorOuter_, typename VectorInner_>
    inline std::shared_ptr<ConvertPrecond<VectorOuter_, VectorInner_>> new_convert_precond(std::shared_ptr<SolverBase<VectorInner_> > inner_solver,
      Global::Gate<typename VectorInner_::LocalVectorType> * gate_inner, Global::Gate<typename VectorOuter_::LocalVectorType> * gate_outer)
    {
      return std::make_shared<ConvertPrecond<VectorOuter_, VectorInner_>>(inner_solver, gate_inner, gate_outer);
    }
  } // namespace Global
} // namespace FEAT

#endif // KERNEL_GLOBAL_CONVERT_PRECOND_HPP
