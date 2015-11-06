#pragma once
#ifndef KERNEL_SOLVER_SCHUR_PRECOND_HPP
#define KERNEL_SOLVER_SCHUR_PRECOND_HPP 1

// includes, FEAST
#include <kernel/solver/base.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>
#include <kernel/lafem/tuple_filter.hpp>
#include <kernel/lafem/tuple_vector.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/global/filter.hpp>

namespace FEAST
{
  namespace Solver
  {
    /**
     * \brief Schur-Complement preconditioner type.
     *
     * This enumeration specifies the various preconditioner types implemented in the
     * SchurPrecond class template.
     */
    enum class SchurType
    {
      /// diagonal Schur-Complement preconditioner
      diagonal,
      /// lower-diagonal Schur-Complement preconditioner
      lower,
      /// upper-diagonal Schur-Complement preconditioner
      upper,
      /// full Schur-Complement preconditioner
      full
    };

    /**
     * \brief Schur-Complement preconditioner
     *
     * This class implements the Schur-Complement preconditioner, which is a special preconditioner for
     * saddle-point systems of the form
     * \f[\begin{bmatrix} A & B\\D & 0\end{bmatrix} \cdot \begin{bmatrix} x_u\\x_p\end{bmatrix} = \begin{bmatrix} f_u\\f_p\end{bmatrix}\f]
     *
     * Let \f$ S \approx -DA^{-1}B\f$ be an approximation of the Schur-complement matrix, then this
     * class supports a total of four different preconditioners for the saddle-point system above:
     * - SchurType::diagonal: Solves the system
     * \f[\begin{bmatrix} A & 0\\0 & S\end{bmatrix} \cdot \begin{bmatrix} x_u\\x_p\end{bmatrix} = \begin{bmatrix} f_u\\f_p\end{bmatrix}\f]
     *
     * - SchurType::lower: Solves the system
     * \f[\begin{bmatrix} A & 0\\D & S\end{bmatrix} \cdot \begin{bmatrix} x_u\\x_p\end{bmatrix} = \begin{bmatrix} f_u\\f_p\end{bmatrix}\f]
     *
     * - SchurType::upper: Solves the system
     * \f[\begin{bmatrix} A & B\\0 & S\end{bmatrix} \cdot \begin{bmatrix} x_u\\x_p\end{bmatrix} = \begin{bmatrix} f_u\\f_p\end{bmatrix}\f]
     *
     * - SchurType::full: Solves the system
     * \f[\begin{bmatrix} I & 0\\-DA^{-1} & I\end{bmatrix} \cdot \begin{bmatrix} A & B\\0 & S\end{bmatrix} \cdot \begin{bmatrix} x_u\\x_p\end{bmatrix} = \begin{bmatrix} f_u\\f_p\end{bmatrix}\f]
     *
     * The required solution steps of \f$ A^{-1} \f$ and \f$ S^{-1} \f$ are performed by two sub-solvers,
     * which have to be created by the user and supplied to the constructor of this object.
     *
     * \tparam MatrixA_, MatrixB_, MatrixD_
     * The types of the sub-matrices A, B and D.
     *
     * \tparam FilterV_, FilterP_
     * The types of the filters for the velocity and pressure components, respectively.
     *
     * \note
     * This class template is specialised for Global::Matrix and Global::Filter instances.
     *
     * \author Peter Zajac
     */
    template<typename MatrixA_, typename MatrixB_, typename MatrixD_, typename FilterV_, typename FilterP_>
    class SchurPrecond :
      public SolverBase<LAFEM::TupleVector<typename MatrixB_::VectorTypeL, typename MatrixD_::VectorTypeL> >
    {
    public:
      /// our data type
      typedef typename MatrixA_::DataType DataType;
      /// velocity vector type
      typedef typename MatrixB_::VectorTypeL VectorTypeV;
      /// pressure vector type
      typedef typename MatrixD_::VectorTypeL VectorTypeP;
      /// system vector type
      typedef LAFEM::TupleVector<VectorTypeV, VectorTypeP> VectorType;
      /// base-class typedef
      typedef SolverBase<VectorType> BaseClass;

      /// A-block solver type
      typedef SolverBase<VectorTypeV> SolverA;
      /// S-block solver type
      typedef SolverBase<VectorTypeP> SolverS;

    private:
      const MatrixA_& _matrix_a;
      const MatrixB_& _matrix_b;
      const MatrixD_& _matrix_d;
      /// our system filter
      const FilterV_& _filter_v;
      const FilterP_& _filter_p;
      /// our A-block solver
      std::shared_ptr<SolverA> _solver_a;
      /// our S-block solver
      std::shared_ptr<SolverS> _solver_s;
      /// our Schur-Complement type
      SchurType _schur_type;
      /// auto-initialise of S-solver
      const bool _auto_init_s;
      /// a temporary defec vector
      VectorTypeV _vec_tmp_v;
      VectorTypeP _vec_tmp_p;

    public:
      /**
       * \brief Constructs a Schur-Complement preconditioner
       *
       * \param[in] matrix
       * The saddle-point system matrix.
       *
       * \param[in] filter
       * The system filter.
       *
       * \param[in] solver_a
       * The solver representing \f$A^{-1}\f$.
       *
       * \param[in] solver_s
       * The solver representing \f$S^{-1}\f$.
       *
       * \param[in] type
       * Specifies the type of the preconditioner. See this class' documentation for details.
       *
       * \param[in] auto_init_s
       * Specifies whether this solver object should call the init/done functions of the \p solver_s
       * object. If set to \c false, then the caller is responsible for the initialisation and
       * finalisation of the S-block solver object.
       */
      explicit SchurPrecond(
        const MatrixA_& matrix_a,
        const MatrixB_& matrix_b,
        const MatrixD_& matrix_d,
        const FilterV_& filter_v,
        const FilterP_& filter_p,
        std::shared_ptr<SolverA> solver_a,
        std::shared_ptr<SolverS> solver_s,
        SchurType type = SchurType::diagonal,
        bool auto_init_s = true) :
        _matrix_a(matrix_a),
        _matrix_b(matrix_b),
        _matrix_d(matrix_d),
        _filter_v(filter_v),
        _filter_p(filter_p),
        _solver_a(solver_a),
        _solver_s(solver_s),
        _schur_type(type),
        _auto_init_s(auto_init_s)
      {
        ASSERT_(solver_a != nullptr);
        ASSERT_(solver_s != nullptr);
      }

      virtual ~SchurPrecond()
      {
      }

      virtual String name() const override
      {
        return "Schur";
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();
        _solver_a->init_symbolic();
        if(_auto_init_s)
        {
          _solver_s->init_symbolic();
        }

        // create a temporary vector
        if(_schur_type != SchurType::diagonal)
        {
          _vec_tmp_v = _matrix_b.create_vector_l();
          _vec_tmp_p = _matrix_d.create_vector_l();
        }
      }

      virtual void init_numeric() override
      {
        BaseClass::init_numeric();
        _solver_a->init_numeric();
        if(_auto_init_s)
        {
          _solver_s->init_numeric();
        }
      }

      virtual void init_branch(String parent = "") override
      {
        BaseClass::init_branch(parent);
        _solver_a->init_branch(parent + "::" + this->name());
        _solver_s->init_branch(parent + "::" + this->name());
      }

      virtual void done_numeric() override
      {
        if(_auto_init_s)
        {
          _solver_s->done_numeric();
        }
        _solver_a->done_numeric();
        BaseClass::done_numeric();
      }

      virtual void done_symbolic() override
      {
        if(_schur_type != SchurType::diagonal)
        {
          _vec_tmp_p.clear();
          _vec_tmp_v.clear();
        }
        if(_auto_init_s)
        {
          _solver_s->done_symbolic();
        }
        _solver_a->done_symbolic();
        BaseClass::done_symbolic();
      }

      virtual String get_formated_solver_tree() override
      {
        String result;
        result += this->name();
        result += " ( ";
        result += "S: " + _solver_s->get_formated_solver_tree() + " /  A: " + _solver_a->get_formated_solver_tree();
        result += " ) ";
        return result;
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        // fetch the references
        VectorTypeV& tmp_v = this->_vec_tmp_v;
        VectorTypeP& tmp_p = this->_vec_tmp_p;
        VectorTypeV& sol_v = vec_cor.template at<Index(0)>();
        VectorTypeP& sol_p = vec_cor.template at<Index(1)>();
        const VectorTypeV& rhs_v = vec_def.template at<Index(0)>();
        const VectorTypeP& rhs_p = vec_def.template at<Index(1)>();

        // now let's check the preconditioner type
        switch(_schur_type)
        {
        case SchurType::diagonal:
          // solve A*u_v = f_v
          if(!status_success(_solver_a->apply(sol_v, rhs_v)))
            return Status::aborted;

          // solve S*u_p = f_p
          if(!status_success(_solver_s->apply(sol_p, rhs_p)))
            return Status::aborted;

          // okay
          return Status::success;

        case SchurType::lower:
          // solve A*u_v = f_v
          if(!status_success(_solver_a->apply(sol_v, rhs_v)))
            return Status::aborted;

          // compute g_p := f_p - D*u_v
          this->_matrix_d.apply(tmp_p, sol_v, rhs_p, -DataType(1));

          // apply pressure defect filter
          this->_filter_p.filter_def(tmp_p);

          // solve S*u_p = g_p
          if(!status_success(_solver_s->apply(sol_p, tmp_p)))
            return Status::aborted;

          // okay
          return Status::success;

        case SchurType::upper:
          // solve S*u_p = f_p
          if(!status_success(_solver_s->apply(sol_p, rhs_p)))
            return Status::aborted;

          // compute g_v := f_v - B*u_p
          this->_matrix_b.apply(tmp_v, sol_p, rhs_v, -DataType(1));

          // apply velocity defect filter
          this->_filter_v.filter_def(tmp_v);

          // solve A*u_v = g_v
          if(!status_success(_solver_a->apply(sol_v, tmp_v)))
            return Status::aborted;

          // okay
          return Status::success;

        case SchurType::full:
          // Note: We will use the first component of the solution vector here.
          //       It will be overwritten by the third solution step below.
          // solve A*u_v = f_v
          if(!status_success(_solver_a->apply(sol_v, rhs_v)))
            return Status::aborted;

          // compute g_p := f_p - D*u_v
          this->_matrix_d.apply(tmp_p, sol_v, rhs_p, -DataType(1));

          // apply pressure defect filter
          this->_filter_p.filter_def(tmp_p);

          // solve S*u_p = g_p
          if(!status_success(_solver_s->apply(sol_p, tmp_p)))
            return Status::aborted;

          // compute g_v := f_v - B*u_p
          this->_matrix_b.apply(tmp_v, sol_p, rhs_v, -DataType(1));

          // apply velocity defect filter
          this->_filter_v.filter_def(tmp_v);

          // solve A*u_v = g_v
          if(!status_success(_solver_a->apply(sol_v, tmp_v)))
            return Status::aborted;

          // okay
          return Status::success;
        }

        // we should never come out here...
        return Status::aborted;
      }
    }; // class SchurPrecond<...>

    /**
     * \brief SchurPrecond specialisation for Global matrices
     *
     * \author Peter Zajac
     */
    template<typename MatrixA_, typename MatrixB_, typename MatrixD_, typename FilterV_, typename FilterP_>
    class SchurPrecond
      <
        Global::Matrix<MatrixA_>, Global::Matrix<MatrixB_>, Global::Matrix<MatrixD_>,
        Global::Filter<FilterV_>, Global::Filter<FilterP_>
      > :
      public Solver::SolverBase<Global::Vector<LAFEM::TupleVector<typename MatrixB_::VectorTypeL, typename MatrixD_::VectorTypeL>>>
    {
    public:
      typedef LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_> LocalMatrixType;
      typedef LAFEM::TupleFilter<FilterV_, FilterP_> LocalFilterType;
      typedef LAFEM::TupleVector<typename MatrixB_::VectorTypeL, typename MatrixD_::VectorTypeL> LocalVectorType;

      typedef Global::Matrix<LocalMatrixType> GlobalMatrixType;
      typedef Global::Filter<LocalFilterType> GlobalFilterType;
      typedef Global::Vector<LocalVectorType> GlobalVectorType;

      typedef typename GlobalMatrixType::DataType DataType;

      typedef Solver::SolverBase<GlobalVectorType> BaseClass;

      typedef Global::Matrix<MatrixA_> GlobalMatrixTypeA;
      typedef Global::Matrix<MatrixB_> GlobalMatrixTypeB;
      typedef Global::Matrix<MatrixD_> GlobalMatrixTypeD;

      typedef typename MatrixB_::VectorTypeL LocalVectorTypeV;
      typedef typename MatrixD_::VectorTypeL LocalVectorTypeP;

      typedef Global::Vector<LocalVectorTypeV> GlobalVectorTypeV;
      typedef Global::Vector<LocalVectorTypeP> GlobalVectorTypeP;

      typedef Global::Filter<FilterV_> GlobalFilterTypeV;
      typedef Global::Filter<FilterP_> GlobalFilterTypeP;

      typedef Global::Gate<LocalVectorTypeV> GateTypeV;
      typedef Global::Gate<LocalVectorTypeP> GateTypeP;

      typedef Solver::SolverBase<GlobalVectorTypeV> SolverA;
      typedef Solver::SolverBase<GlobalVectorTypeP> SolverS;

    protected:
      const GlobalMatrixTypeA& _matrix_a;
      const GlobalMatrixTypeB& _matrix_b;
      const GlobalMatrixTypeD& _matrix_d;
      const GlobalFilterTypeV& _filter_v;
      const GlobalFilterTypeP& _filter_p;
      GlobalVectorTypeV _vec_rhs_v, _vec_sol_v, _vec_def_v;
      GlobalVectorTypeP _vec_rhs_p, _vec_sol_p, _vec_def_p;
      std::shared_ptr<SolverA> _solver_a;
      std::shared_ptr<SolverS> _solver_s;
      SchurType _schur_type;
      const bool _auto_init_s;

    public:
      explicit SchurPrecond(
        const GlobalMatrixTypeA& matrix_a,
        const GlobalMatrixTypeB& matrix_b,
        const GlobalMatrixTypeD& matrix_d,
        const GlobalFilterTypeV& filter_v,
        const GlobalFilterTypeP& filter_p,
        std::shared_ptr<SolverA> solver_a,
        std::shared_ptr<SolverS> solver_s,
        SchurType type = SchurType::diagonal,
        bool auto_init_s = true
        ) :
        _matrix_a(matrix_a),
        _matrix_b(matrix_b),
        _matrix_d(matrix_d),
        _filter_v(filter_v),
        _filter_p(filter_p),
        _solver_a(solver_a),
        _solver_s(solver_s),
        _schur_type(type),
        _auto_init_s(auto_init_s)
      {
      }

      virtual String name() const override
      {
        return "Schur";
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();
        _solver_a->init_symbolic();
        if(_auto_init_s)
        {
          _solver_s->init_symbolic();
        }

        _vec_rhs_v = _matrix_b.create_vector_l();
        _vec_sol_v = _matrix_b.create_vector_l();
        _vec_rhs_p = _matrix_d.create_vector_l();
        _vec_sol_p = _matrix_d.create_vector_l();

        if(_schur_type != SchurType::diagonal)
        {
          _vec_def_v = _matrix_b.create_vector_l();
          _vec_def_p = _matrix_d.create_vector_l();
        }
      }

      virtual void init_numeric() override
      {
        BaseClass::init_numeric();
        _solver_a->init_numeric();
        if(_auto_init_s)
        {
          _solver_s->init_numeric();
        }
      }

      virtual void init_branch(String parent = "") override
      {
        BaseClass::init_branch(parent);
        _solver_a->init_branch(parent + "::" + this->name());
        _solver_s->init_branch(parent + "::" + this->name());
      }

      virtual void done_numeric() override
      {
        if(_auto_init_s)
        {
          _solver_s->done_numeric();
        }
        _solver_a->done_numeric();
        BaseClass::done_numeric();
      }

      virtual void done_symbolic() override
      {
        if(_schur_type != SchurType::diagonal)
        {
          _vec_def_p.clear();
          _vec_def_v.clear();
        }
        _vec_sol_p.clear();
        _vec_rhs_p.clear();
        _vec_sol_v.clear();
        _vec_rhs_v.clear();
        if(_auto_init_s)
        {
          _solver_s->done_symbolic();
        }
        _solver_a->done_symbolic();
        BaseClass::done_symbolic();
      }

      virtual String get_formated_solver_tree() override
      {
        String result;
        result += this->name();
        result += " ( ";
        result += "S: " + _solver_s->get_formated_solver_tree() + " /  A: " + _solver_a->get_formated_solver_tree();
        result += " ) ";
        return result;
      }

      virtual Status apply(GlobalVectorType& vec_cor, const GlobalVectorType& vec_def) override
      {
        // first of all, copy RHS
        (*_vec_rhs_v).copy((*vec_def).template at<0>());
        (*_vec_rhs_p).copy((*vec_def).template at<1>());

        // now let's check the preconditioner type
        switch(_schur_type)
        {
        case SchurType::diagonal:
          // solve A*u_v = f_v
          if(!status_success(_solver_a->apply(_vec_sol_v, _vec_rhs_v)))
            return Status::aborted;

          // solve S*u_p = f_p
          if(!status_success(_solver_s->apply(_vec_sol_p, _vec_rhs_p)))
            return Status::aborted;

          // okay
          break;

        case SchurType::lower:
          // solve A*u_v = f_v
          if(!status_success(_solver_a->apply(_vec_sol_v, _vec_rhs_v)))
            return Status::aborted;

          // compute g_p := f_p - D*u_v
          _matrix_d.apply(_vec_def_p, _vec_sol_v, _vec_rhs_p, -DataType(1));

          // apply pressure defect filter
          _filter_p.filter_def(_vec_def_p);

          // solve S*u_p = g_p
          if(!status_success(_solver_s->apply(_vec_sol_p, _vec_def_p)))
            return Status::aborted;

          // okay
          break;

        case SchurType::upper:
          // solve S*u_p = f_p
          if(!status_success(_solver_s->apply(_vec_sol_p, _vec_rhs_p)))
            return Status::aborted;

          // compute g_v := f_v - B*u_p
          _matrix_b.apply(_vec_def_v, _vec_sol_p, _vec_rhs_v, -DataType(1));

          // apply velocity defect filter
          _filter_v.filter_def(_vec_def_v);

          // solve A*u_v = g_v
          if(!status_success(_solver_a->apply(_vec_sol_v, _vec_def_v)))
            return Status::aborted;

          // okay
          break;

        case SchurType::full:
          // Note: We will use the first component of the solution vector here.
          //       It will be overwritten by the third solution step below.
          // solve A*u_v = f_v
          if(!status_success(_solver_a->apply(_vec_sol_v, _vec_rhs_v)))
            return Status::aborted;

          // compute g_p := f_p - D*u_v
          _matrix_d.apply(_vec_def_p, _vec_sol_v, _vec_rhs_p, -DataType(1));

          // apply pressure defect filter
          _filter_p.filter_def(_vec_def_p);

          // solve S*u_p = g_p
          if(!status_success(_solver_s->apply(_vec_sol_p, _vec_def_p)))
            return Status::aborted;

          // compute g_v := f_v - B*u_p
          _matrix_b.apply(_vec_def_v, _vec_sol_p, _vec_rhs_v, -DataType(1));

          // apply velocity defect filter
          _filter_v.filter_def(_vec_def_v);

          // solve A*u_v = g_v
          if(!status_success(_solver_a->apply(_vec_sol_v, _vec_def_v)))
            return Status::aborted;

          // okay
          break;
        }

        // finally, copy sol
        (*vec_cor).template at<0>().copy(*_vec_sol_v);
        (*vec_cor).template at<1>().copy(*_vec_sol_p);

        // okay
        return Status::success;
      }
    }; // class SchurPrecond<...>

    /**
     * \brief Creates a new SchurPrecond solver object
     *
     * \param[in] matrix_a, matrix_b, matrix_d
     * The three sub-matrices of the saddle-point system matrix.
     *
     * \param[in] filter_v, filter_p
     * The two sub-filters of the system.
     *
     * \param[in] solver_a
     * The solver object for the A-matrix block.
     *
     * \param[in] solver_s
     * The solver object for the S-matrix block.
     *
     * \param[in] type
     * Specifies the type of the preconditioner.
     *
     * \param[in] auto_init_s
     * Specifies whether to automatically initialise the S-matrix solver.
     *
     * \returns
     * A shared pointer to a new SchurPrecond object.
     */
     /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAST_COMPILER_GNU) && (FEAST_COMPILER_GNU < 40900)
    template<typename MatrixA_, typename MatrixB_, typename MatrixD_, typename FilterV_, typename FilterP_, typename SolverA_, typename SolverS_>
    inline std::shared_ptr<SchurPrecond<MatrixA_, MatrixB_, MatrixD_, FilterV_, FilterP_>> new_schur_precond(
      const MatrixA_& matrix_a, const MatrixB_& matrix_b, const MatrixD_& matrix_d,
      const FilterV_& filter_v, const FilterP_& filter_p,
      std::shared_ptr<SolverA_> solver_a,
      std::shared_ptr<SolverS_> solver_s,
      SchurType type = SchurType::diagonal,
      bool auto_init_s = true)
    {
      return std::make_shared<SchurPrecond<MatrixA_, MatrixB_, MatrixD_, FilterV_, FilterP_>>
        (matrix_a, matrix_b, matrix_d, filter_v, filter_p, solver_a, solver_s, type, auto_init_s);
    }
#else
    template<typename MatrixA_, typename MatrixB_, typename MatrixD_, typename FilterV_, typename FilterP_>
    inline std::shared_ptr<SchurPrecond<MatrixA_, MatrixB_, MatrixD_, FilterV_, FilterP_>> new_schur_precond(
      const MatrixA_& matrix_a, const MatrixB_& matrix_b, const MatrixD_& matrix_d,
      const FilterV_& filter_v, const FilterP_& filter_p,
      std::shared_ptr<SolverBase<typename MatrixB_::VectorTypeL>> solver_a,
      std::shared_ptr<SolverBase<typename MatrixD_::VectorTypeL>> solver_s,
      SchurType type = SchurType::diagonal,
      bool auto_init_s = true)
    {
      return std::make_shared<SchurPrecond<MatrixA_, MatrixB_, MatrixD_, FilterV_, FilterP_>>
        (matrix_a, matrix_b, matrix_d, filter_v, filter_p, solver_a, solver_s, type, auto_init_s);
    }
#endif
  } // namespace Solver
} // namespace FEAST

#endif // KERNEL_SOLVER_SCHUR_PRECOND_HPP
