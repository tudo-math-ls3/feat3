// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/solver/base.hpp>
#include <kernel/global/vector.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Solver conversion module.
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
     * \note
     * This class template is specialized for Global::Matrix and Global::Filter instances below.
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

      using VectorTypeInner = VectorInner_;
      using SolverTypeInner = SolverBase<VectorTypeInner>;

    protected:
      std::shared_ptr<SolverTypeInner> _inner_solver;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] inner_solver
       * The actual solver, which shall be executed with converted rhs/sol vectors
       */
      explicit ConvertPrecond(std::shared_ptr<SolverTypeInner> inner_solver) :
        _inner_solver(inner_solver)
      {
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "Convert";
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

      virtual Status apply(VectorTypeOuter& vec_cor, const VectorTypeOuter& vec_def) override
      {
        Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));

        VectorTypeInner vec_def_inner;
        vec_def_inner.convert(vec_def);
        VectorTypeInner vec_cor_inner(vec_def_inner.clone(LAFEM::CloneMode::Layout));

        Statistics::add_solver_expression(std::make_shared<ExpressionCallPrecond>(this->name(), this->_inner_solver->name()));
        Status status = _inner_solver->apply(vec_cor_inner, vec_def_inner);
        if(!status_success(status))
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, 0));
          return status;
        }

        vec_cor.convert(vec_cor_inner);
        Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::success, 0));
        return Status::success;
      }
    }; // class ConvertPrecond<...>

    /**
     * \brief Creates a new ConvertPrecond solver object
     *
     * \param[in] inner_solver
     * The actual solver, which shall be executed with converted rhs/sol vectors
     *
     * \returns
     * A shared pointer to a new ConvertPrecond object.
     */
    template<typename VectorOuter_, typename VectorInner_>
    inline std::shared_ptr<ConvertPrecond<VectorOuter_, VectorInner_>> new_convert_precond(std::shared_ptr<SolverBase<VectorInner_> > inner_solver)
    {
      return std::make_shared<ConvertPrecond<VectorOuter_, VectorInner_>>(inner_solver);
    }


    /**
     * \brief Convert preconditioner specialization for Global::Vector
     */
    template<typename LocalVectorOuter_, typename LocalVectorInner_, typename MirrorOuter_, typename MirrorInner_>
    class ConvertPrecond<Global::Vector<LocalVectorOuter_, MirrorOuter_>, Global::Vector<LocalVectorInner_, MirrorInner_>> :
      public SolverBase<Global::Vector<LocalVectorOuter_, MirrorOuter_>>
    {
    public:
      using VectorTypeOuter = Global::Vector<LocalVectorOuter_, MirrorOuter_>;
      using BaseClass = SolverBase<VectorTypeOuter>;

      using VectorTypeInner = Global::Vector<LocalVectorInner_, MirrorInner_>;
      using SolverTypeInner = SolverBase<VectorTypeInner>;

      //using GateTypeOuter = Gate<VectorOuter_, MirrorOuter_>;
      //using GateTypeInner = Gate<VectorInner_, MirrorInner_>;

    protected:
      std::shared_ptr<SolverTypeInner> _inner_solver;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] inner_solver
       * The actual solver, which shall be executed with converted rhs/sol vectors
       */
      explicit ConvertPrecond(std::shared_ptr<SolverTypeInner> inner_solver) :
        _inner_solver(inner_solver)
      {
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "Convert";
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

      virtual Status apply(VectorTypeOuter& vec_cor, const VectorTypeOuter& vec_def) override
      {
        Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));

        VectorTypeInner vec_def_inner;
        vec_def_inner.convert(vec_def_inner.get_gate(), vec_def);
        VectorTypeInner vec_cor_inner(vec_def_inner.clone(LAFEM::CloneMode::Layout));

        Statistics::add_solver_expression(std::make_shared<ExpressionCallPrecond>(this->name(), this->_inner_solver->name()));
        Status status = _inner_solver->apply(vec_cor_inner, vec_def_inner);
        if(!status_success(status))
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, 0));
          return status;
        }

        vec_cor.convert(vec_cor.get_gate(), vec_cor_inner);
        Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::success, 0));
        return Status::success;
      }
    }; // class ConvertPrecond<...>

    /**
     * \brief Creates a new ConvertPrecond solver object
     *
     * \param[in] inner_solver
     * The actual solver, which shall be executed with converted rhs/sol vectors
     *
     * \returns
     * A shared pointer to a new ConvertPrecond object.
     */
    template<typename LocalVectorOuter_, typename LocalVectorInner_, typename MirrorOuter_, typename MirrorInner_>
    inline std::shared_ptr<ConvertPrecond<Global::Vector<LocalVectorOuter_, MirrorOuter_>, Global::Vector<LocalVectorInner_, MirrorInner_>>>
      new_convert_precond(std::shared_ptr<SolverBase<Global::Vector<LocalVectorInner_, MirrorInner_>>> inner_solver)
    {
      return std::make_shared<ConvertPrecond<Global::Vector<LocalVectorOuter_, MirrorOuter_>, Global::Vector<LocalVectorInner_, MirrorInner_>>>(inner_solver);
    }
  } // namespace Solver
} // namespace FEAT
