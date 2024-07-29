// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2024 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SOLVER_DESCENT_HPP
#define KERNEL_SOLVER_DESCENT_HPP 1

// includes, FEAT
#include <kernel/solver/iterative.hpp>

namespace FEAT
{
  namespace Solver
  {
    enum class DescentVariant
    {
      /// simple fixed Richardson iteration
      fixed,
      /// symmetrically preconditioned steepest descent
      steepest,
      /// left-preconditioned defect minimization
      defect,
      /// minimize residual
      //min_res,
    };

    std::ostream& operator<<(std::ostream& os, DescentVariant dv)
    {
      switch(dv)
      {
      case DescentVariant::fixed:
        return os << "fixed";

      case DescentVariant::steepest:
        return os << "steepest";

      case DescentVariant::defect:
        return os << "defect";

      default:
        return os << "???";
      }
    }

    /**
     * \brief (Preconditioned) Descent solver implementation
     *
     * \tparam Matrix_
     * The matrix class to be used by the solver.
     *
     * \tparam Filter_
     * The filter class to be used by the solver.
     *
     * \see Chapter 5.3.2, Algorithm 5.3 in \cite Saad03
     *
     * \author Peter Zajac
     */
    template<
      typename Matrix_,
      typename Filter_>
    class Descent :
      public PreconditionedIterativeSolver<typename Matrix_::VectorTypeR>
    {
    public:
      typedef Matrix_ MatrixType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::VectorTypeR VectorType;
      typedef typename MatrixType::DataType DataType;
      typedef PreconditionedIterativeSolver<VectorType> BaseClass;

      typedef SolverBase<VectorType> PrecondType;

    protected:
      /// the matrix for the solver
      const MatrixType& _system_matrix;
      /// the filter for the solver
      const FilterType& _system_filter;
      /// temporary vectors
      VectorType _vec_r, _vec_p, _vec_q;//, _vec_z;
      /// the chosen descent type
      DescentVariant _variant;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] matrix
       * A reference to the system matrix.
       *
       * \param[in] filter
       * A reference to the system filter.
       *
       * \param[in] variant
       * The descent variant to use
       *
       * \param[in] precond
       * A pointer to the preconditioner. May be \c nullptr.
       */
      explicit Descent(const MatrixType& matrix, const FilterType& filter,
        DescentVariant variant, std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass("Descent", precond),
        _system_matrix(matrix),
        _system_filter(filter),
        _variant(variant)
      {
        // set communicator by system matrix
        this->_set_comm_by_matrix(matrix);
        this->set_plot_name("Descent[" + stringify(this->_variant) + "]");
      }

      /**
       * \brief Constructor using a PropertyMap
       *
       * \param[in] section_name
       * The name of the config section, which it does not know by itself
       *
       * \param[in] section
       * A pointer to the PropertyMap section configuring this solver
       *
       * \param[in] matrix
       * The system matrix.
       *
       * \param[in] filter
       * The system filter.
       *
       * \param[in] precond
       * The preconditioner. May be \c nullptr.
       *
       * \returns
       * A shared pointer to a new Descent object.
       */
      explicit Descent(const String& section_name, PropertyMap* section,
        const MatrixType& matrix, const FilterType& filter,
        std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass("Descent", section_name, section, precond),
        _system_matrix(matrix),
        _system_filter(filter)
      {
        // set communicator by system matrix
        this->_set_comm_by_matrix(matrix);

        // get our variant
        auto variant_p = section->query("variant");
        if(!variant_p.second)
          throw ParseError("Descent: config section is missing the mandatory variant!");
        if(variant_p.first.compare_no_case("steepest") == 0)
          this->_variant = DescentVariant::steepest;
        else if(variant_p.first.compare_no_case("defect") == 0)
          this->_variant = DescentVariant::defect;
        else
          throw ParseError("Descent: invalid variant: " + variant_p.first);

        this->set_plot_name("Descent[" + stringify(this->_variant) + "]");
      }

      virtual String name() const override
      {
        return "Descent";
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();

        _vec_r = this->_system_matrix.create_vector_r();
        _vec_p = this->_system_matrix.create_vector_r();
        _vec_q = this->_system_matrix.create_vector_r();
        //_vec_z = this->_system_matrix.create_vector_r();
      }

      virtual void done_symbolic() override
      {
        //this->_vec_z.clear();
        this->_vec_q.clear();
        this->_vec_p.clear();
        this->_vec_r.clear();
        BaseClass::done_symbolic();
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        // save defect
        this->_vec_r.copy(vec_def);

        // clear solution vector
        vec_cor.format();

        // apply solver
        this->_status = _apply_intern(vec_cor);

        // plot summary
        this->plot_summary();

        // return status
        return this->_status;
      }

      virtual Status correct(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // compute defect
        this->_system_matrix.apply(this->_vec_r, vec_sol, vec_rhs, -DataType(1));
        this->_system_filter.filter_def(this->_vec_r);

        // apply solver
        this->_status = _apply_intern(vec_sol);

        // plot summary
        this->plot_summary();

        // return status
        return this->_status;
      }

    protected:
      virtual Status _apply_intern(VectorType& vec_sol)
      {
        IterationStats pre_iter(*this);
        Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));

        const MatrixType& matrix(this->_system_matrix);
        const FilterType& filter(this->_system_filter);
        VectorType& vec_r(this->_vec_r);
        VectorType& vec_p(this->_vec_p);
        VectorType& vec_q(this->_vec_q);
        //VectorType& vec_z(this->_vec_z);

        // set initial defect:
        // r[0] := b - A*x[0]
        Status status = this->_set_initial_defect(vec_r, vec_sol);
        if(status != Status::progress)
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
          pre_iter.destroy();
          return status;
        }

        // apply preconditioner to defect vector
        // p[0] := M^{-1} * r[0]
        if(!this->_apply_precond(vec_p, vec_r, filter))
        {
          pre_iter.destroy();
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
          return Status::aborted;
        }

        pre_iter.destroy();

        // start iterating
        while(status == Status::progress)
        {
          IterationStats stat(*this);

          // q[k] := A*p[k]
          matrix.apply(vec_q, vec_p);
          filter.filter_def(vec_q);

          // compute alpha
          DataType gamma = DataType(1), omega = DataType(1);
          switch(this->_variant)
          {
          case DescentVariant::fixed:
            gamma = omega = DataType(1);
            break;

          case DescentVariant::steepest:
            // alpha[k] := < r[k], p[k] > / < q[k], p[k] >
            gamma = vec_r.dot(vec_p);
            omega = vec_q.dot(vec_p);
            break;

          case DescentVariant::defect:
            // alpha[k] := < r[k], q[k] > / < q[k], q[k] >
            gamma = vec_r.dot(vec_q);
            omega = vec_q.dot(vec_q);
            break;
          }

          // compute alpha
          DataType alpha = gamma / omega;

          // update solution vector:
          // x[k+1] := x[k] + alpha[k] * p[k]
          vec_sol.axpy(vec_p, vec_sol, alpha);

          // update defect vector:
          // r[k+1] := r[k] - alpha[k] * q[k]
          vec_r.axpy(vec_q, vec_r, -alpha);

          // compute defect norm
          status = this->_set_new_defect(vec_r, vec_sol);
          if(status != Status::progress)
          {
            stat.destroy();
            Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
            return status;
          }

          // apply preconditioner
          // p[k+1] := M^{-1} * r[k+1]
          if(!this->_apply_precond(vec_p, vec_r, filter))
          {
            stat.destroy();
            Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
            return Status::aborted;
          }
        }

        // we should never reach this point...
        Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::undefined, this->get_num_iter()));
        return Status::undefined;
      }
    }; // class Descent<...>

    /**
     * \brief Creates a new Descent solver object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] variant
     * The desired descent variant.
     *
     * \param[in] precond
     * The preconditioner. May be \c nullptr.
     *
     * \returns
     * A shared pointer to a new Descent object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<Descent<Matrix_, Filter_>> new_descent(
      const Matrix_& matrix, const Filter_& filter, DescentVariant variant,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr)
    {
      return std::make_shared<Descent<Matrix_, Filter_>>(matrix, filter, variant, precond);
    }

    /**
     * \brief Creates a new Descent solver object using a PropertyMap
     *
     * \param[in] section_name
     * The name of the config section, which it does not know by itself
     *
     * \param[in] section
     * A pointer to the PropertyMap section configuring this solver
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] precond
     * The preconditioner. May be \c nullptr.
     *
     * \returns
     * A shared pointer to a new Descent object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<Descent<Matrix_, Filter_>> new_pmr(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr)
    {
      return std::make_shared<Descent<Matrix_, Filter_>>(section_name, section, matrix, filter, precond);
    }
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_DESCENT_HPP
