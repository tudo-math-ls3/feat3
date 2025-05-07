// FEAT3: Finite Element Analysis Toolbox, Version 3
// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/global/filter.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>
#include <kernel/lafem/tuple_filter.hpp>
#include <kernel/lafem/tuple_mirror.hpp>
#include <kernel/lafem/tuple_vector.hpp>
#include <kernel/solver/base.hpp>

namespace FEAT
{
  namespace Solver
  {

    /**
     * \brief BFBT Schur-complement preconditioner
     *
     * This class implements the BFBT preconditioner, which provides an approximation
     * to the Schur complement \f$ S \approx -DA^{-1}B \f$ of a saddle-point system
     * of the form:
     * \f[
     * \begin{bmatrix} A & B \\ D & 0 \end{bmatrix}
     * \cdot
     * \begin{bmatrix} x_u \\ x_p \end{bmatrix}
     * =
     * \begin{bmatrix} f_u \\ f_p \end{bmatrix}
     * \f]
     *
     * The BFBT preconditioner uses a mass matrix-based approximation of the inverse operator
     * \f$ A^{-1} \f$, i.e., the Schur complement is approximated as:
     * \f[
     * S_{\mathrm{BFBT}} = L_p^{-1} D M^{-1} A M^{-1} B L_p^{-1}
     * \f]
     * where \f$ M \f$ is the lumped velocity mass matrix and \f$ L_p \f$ is the pressure laplace matrix.
     * For more information see \cite{elman2006block}.
     *
     * Separate left and right solvers for \f$ L_p^{-1} \f$ can be provided
     * by the user for greater flexibility.
     *
     * \tparam MatrixA_, MatrixB_, MatrixD_
     * The types of the sub-matrices A, B, and D.
     *
     * \tparam FilterV_, FilterP_
     * The types of the filters for the velocity and pressure components.
     *
     * \note
     * This class template is specialized for Global::Matrix and Global::Filter instances below.
     *
     * \author Pia Ritter
     */
    template <typename MatrixA_, typename MatrixB_, typename MatrixD_, typename FilterV_, typename FilterP_>
    class BFBT :
      public SolverBase<typename MatrixD_::VectorTypeL>
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
      typedef SolverBase<VectorTypeP> BaseClass;

      /// left solver type
      typedef SolverBase<VectorTypeP> SolverLeft;
      /// right solver type
      typedef SolverBase<VectorTypeP> SolverRight;

    private:
      const MatrixA_& _matrix_a;
      const MatrixB_& _matrix_b;
      const MatrixD_& _matrix_d;
      /// our system filter
      const FilterV_& _filter_v;
      const FilterP_& _filter_p;

      /// our left solver
      std::shared_ptr<SolverLeft> _solver_left;
      /// our right solver
      std::shared_ptr<SolverRight> _solver_right;

      /// a temporary defect vector
      VectorTypeV _vec_tmp_v1;
      VectorTypeV _vec_tmp_v2;

      // inverse lumped velocity mass vector
      const VectorTypeV& _lumped_velo_mass_vec;

    public:
      /**
       * \brief Constructs a scaled BFBT preconditioner
       *
       * \param[in] matrix_a, matrix_b, matrix_d
       * The three sub-matrices of the saddle-point matrix.
       *
       * \param[in] filter_v, filter_p
       * The velocity and pressure filter, resp.
       *
       * \param[in] solver_left
       * The solver representing left \f$L_p^{-1}\f$.
       *
       * \param[in] solver_right
       * The solver representing right \f$L_p^{-1}\f$.
       *
       */
      explicit BFBT(
        const MatrixA_& matrix_a,
        const MatrixB_& matrix_b,
        const MatrixD_& matrix_d,
        const FilterV_& filter_v,
        const FilterP_& filter_p,
        std::shared_ptr<SolverLeft> solver_left,
        std::shared_ptr<SolverRight> solver_right,
        const VectorTypeV& lumped_velo_mass_vec) :
        _matrix_a(matrix_a),
        _matrix_b(matrix_b),
        _matrix_d(matrix_d),
        _filter_v(filter_v),
        _filter_p(filter_p),
        _solver_left(solver_left),
        _solver_right(solver_right),
        _lumped_velo_mass_vec(lumped_velo_mass_vec)
      {
        XASSERTM(solver_left != nullptr, "left-solver must be given");
        XASSERTM(solver_right != nullptr, "right-solver must be given");
      }

      virtual ~BFBT()
      {
      }

      virtual String name() const override
      {
        return "BFBT";
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();
        if(_solver_left == _solver_right)
        {
          _solver_left->init_symbolic();
        }
        else
        {
          _solver_left->init_symbolic();
          _solver_right->init_symbolic();
        }
        // create a temporary vector
        _vec_tmp_v1 = _matrix_b.create_vector_l();
        _vec_tmp_v2 = _matrix_b.create_vector_l();
      }

      virtual void init_numeric() override
      {
        BaseClass::init_numeric();
        if(_solver_left == _solver_right)
        {
          _solver_left->init_numeric();
        }
        else
        {
          _solver_left->init_numeric();
          _solver_right->init_numeric();
        }
      }

      virtual void done_numeric() override
      {
        _vec_tmp_v1.clear();
        _vec_tmp_v2.clear();
        if(_solver_left == _solver_right)
        {
          _solver_left->done_numeric();
        }
        else
        {
          _solver_left->done_numeric();
          _solver_right->done_numeric();
        }
        BaseClass::done_numeric();
      }

      virtual void done_symbolic() override
      {
        if(_solver_left == _solver_right)
        {
          _solver_left->done_symbolic();
        }
        else
        {
          _solver_left->done_symbolic();
          _solver_right->done_symbolic();
        }
        BaseClass::done_symbolic();
      }


      virtual Status apply(VectorTypeP& vec_cor, const VectorTypeP& vec_def) override
      {
        Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));

        // fetch references
        VectorTypeV& tmp_v1 = this->_vec_tmp_v1;
        VectorTypeV& tmp_v2 = this->_vec_tmp_v2;
        VectorTypeP vec_cor_right(vec_cor.clone(LAFEM::CloneMode::Layout));
        VectorTypeP vec_def_left(vec_cor.clone(LAFEM::CloneMode::Layout));

        // vec_cor_right = L_right^{-1} * vec_def
        Status status_right = _solver_right->apply(vec_cor_right, vec_def);
        if(!status_success(status_right))
          return Status::aborted;

        // tmp1 = B * vec_cor_right
        _matrix_b.apply(tmp_v1, vec_cor_right);
        this->_filter_v.filter_def(tmp_v1);

        // tmp1 = M^{-1} * tmp1
        tmp_v1.component_product(tmp_v1, _lumped_velo_mass_vec);
        this->_filter_v.filter_def(tmp_v1);

        // tmp2 = A * tmp1
        _matrix_a.apply(tmp_v2, tmp_v1);
        this->_filter_v.filter_def(tmp_v2);

        // tmp2 = M^{-1} * tmp2
        tmp_v2.component_product(tmp_v2, _lumped_velo_mass_vec);
        this->_filter_v.filter_def(tmp_v2);

        // vec_tmp = D * tmp2
        _matrix_d.apply(vec_def_left, tmp_v2);
        this->_filter_p.filter_def(vec_def_left);

        // vec_cor = L_left^{-1} * vec_tmp
        Status status_left = _solver_left->apply(vec_cor, vec_def_left);
        if(!status_success(status_left))
          return Status::aborted;

        vec_cor.scale(vec_cor, -1);
        return Status::success;
      }
    }; // class BFBT<...>

    /**
     * \brief BFBT specialization for Global matrices
     *
     * \author Pia Ritter
     */
    template <
      typename MatrixA_, typename MatrixB_, typename MatrixD_,
      typename FilterV_, typename FilterP_, typename MirrorV_,
      typename MirrorP_>
    class BFBT
      <
        Global::Matrix<MatrixA_, MirrorV_, MirrorV_>,
        Global::Matrix<MatrixB_, MirrorV_, MirrorP_>,
        Global::Matrix<MatrixD_, MirrorP_, MirrorV_>,
        Global::Filter<FilterV_, MirrorV_>,
        Global::Filter<FilterP_, MirrorP_>
      > :
      public Solver::SolverBase<
        Global::Vector<
          LAFEM::TupleVector<
            typename MatrixB_::VectorTypeL,
            typename MatrixD_::VectorTypeL>,
          LAFEM::TupleMirror<
            MirrorV_, MirrorP_>
        >
      >
    {
    public:
      typedef LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_> LocalMatrixType;
      typedef LAFEM::TupleVector<typename MatrixB_::VectorTypeL, typename MatrixD_::VectorTypeL> LocalVectorType;
      typedef LAFEM::TupleMirror<MirrorV_, MirrorP_> MirrorType;

      typedef Global::Matrix<LocalMatrixType, MirrorType, MirrorType> GlobalMatrixType;
      typedef Global::Vector<LocalVectorType, MirrorType> GlobalVectorType;

      typedef typename GlobalVectorType::DataType DataType;

      typedef Solver::SolverBase<GlobalVectorType> BaseClass;

      typedef Global::Matrix<MatrixA_, MirrorV_, MirrorV_> GlobalMatrixTypeA;
      typedef Global::Matrix<MatrixB_, MirrorV_, MirrorP_> GlobalMatrixTypeB;
      typedef Global::Matrix<MatrixD_, MirrorP_, MirrorV_> GlobalMatrixTypeD;

      typedef Global::Filter<FilterV_, MirrorV_> GlobalFilterTypeV;
      typedef Global::Filter<FilterP_, MirrorP_> GlobalFilterTypeP;

      typedef typename MatrixB_::VectorTypeL LocalVectorTypeV;
      typedef typename MatrixD_::VectorTypeL LocalVectorTypeP;

      typedef Global::Vector<LocalVectorTypeV, MirrorV_> GlobalVectorTypeV;
      typedef Global::Vector<LocalVectorTypeP, MirrorP_> GlobalVectorTypeP;

      typedef Solver::SolverBase<GlobalVectorTypeV> SolverA;
      typedef Solver::SolverBase<GlobalVectorTypeP> SolverS;

    protected:
      const GlobalMatrixTypeA& _matrix_a;
      const GlobalMatrixTypeB& _matrix_b;
      const GlobalMatrixTypeD& _matrix_d;
      const GlobalFilterTypeV& _filter_v;
      const GlobalFilterTypeP& _filter_p;
      std::shared_ptr<SolverA> _solver_left;
      std::shared_ptr<SolverS> _solver_right;
      GlobalVectorTypeV _vec_tmp_v1, _vec_tmp_v2;
      const GlobalVectorTypeV& _lumped_velo_mass_vec;

    public:
      explicit BFBT(
        const GlobalMatrixTypeA& matrix_a,
        const GlobalMatrixTypeB& matrix_b,
        const GlobalMatrixTypeD& matrix_d,
        const GlobalFilterTypeV& filter_v,
        const GlobalFilterTypeP& filter_p,
        std::shared_ptr<SolverA> solver_left,
        std::shared_ptr<SolverS> solver_right,
        const GlobalVectorTypeV& lumped_velo_mass_vec
        ) :
        _matrix_a(matrix_a),
        _matrix_b(matrix_b),
        _matrix_d(matrix_d),
        _filter_v(filter_v),
        _filter_p(filter_p),
        _solver_left(solver_left),
        _solver_right(solver_right),
        _lumped_velo_mass_vec(lumped_velo_mass_vec)
      {
      }

      virtual String name() const override
      {
        return "BFBT";
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();
        _solver_left->init_symbolic();
        if(_solver_right != _solver_left)
          _solver_right->init_symbolic();
        _vec_tmp_v1 = _matrix_b.create_vector_l();
        _vec_tmp_v2 = _matrix_b.create_vector_l();
      }

      virtual void init_numeric() override
      {
        BaseClass::init_numeric();
        _solver_left->init_numeric();
        if(_solver_right != _solver_left)
          _solver_right->init_numeric();
      }

      virtual void done_numeric() override
      {
        _vec_tmp_v1.clear();
        _vec_tmp_v2.clear();
        _solver_left->done_numeric();
        if(_solver_right != _solver_left)
          _solver_right->done_numeric();
        BaseClass::done_numeric();
      }

      virtual void done_symbolic() override
      {
        _solver_left->done_symbolic();
        if(_solver_right != _solver_left)
          _solver_right->done_symbolic();
        BaseClass::done_symbolic();
      }

      virtual Status apply(GlobalVectorTypeP& vec_cor, const GlobalVectorTypeP& vec_def) override
      {
        Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));

        GlobalVectorTypeP vec_cor_right(vec_cor.clone(LAFEM::CloneMode::Layout));
        GlobalVectorTypeP vec_def_left(vec_cor.clone(LAFEM::CloneMode::Layout));

        // vec_cor_right = S^{-1} * vec_def
        if(!status_success(_solver_right->apply(vec_cor_right, vec_def)))
          return Status::aborted;

        // tmp_v1 = B * vec_cor_right
        _matrix_b.apply(_vec_tmp_v1, vec_cor_right);
        this->_filter_v.filter_def(_vec_tmp_v1);

        // tmp_v1 = M^{-1} * tmp_v1
        _vec_tmp_v1.component_product(_vec_tmp_v1, _lumped_velo_mass_vec);
        this->_filter_v.filter_def(_vec_tmp_v1);

        // tmp_v2 = A * tmp_v1
        _matrix_a.apply(_vec_tmp_v2, _vec_tmp_v1);
        this->_filter_v.filter_def(_vec_tmp_v2);

        // tmp_v2 = M^{-1} * tmp_v2
        _vec_tmp_v2.component_product(_vec_tmp_v2, _lumped_velo_mass_vec);
        this->_filter_v.filter_def(_vec_tmp_v2);

        // vec_def_left = D * tmp_v2
        _matrix_d.apply(vec_def_left, _vec_tmp_v2);
        this->_filter_p.filter_def(vec_def_left);

        // vec_cor = -A^{-1} * vec_def_left
        if(!status_success(_solver_left->apply(vec_cor, vec_def_left)))
          return Status::aborted;

        vec_cor.scale(vec_cor, -1.0);
        Statistics::add_solver_expression(
          std::make_shared<ExpressionEndSolve>(this->name(), Status::success, 0));
        return Status::success;
      }
    };

    /**
     * \brief Creates a new BFBT solver object
     *
     * \param[in] matrix_a, matrix_b, matrix_d
     * The three sub-matrices of the saddle-point system matrix.
     *
     * \param[in] filter_v, filter_p
     * The two sub-filters of the system.
     *
     * \param[in] solver_left
     * The solver object for the left L_p-matrix.
     *
     * \param[in] solver_right
     * The solver object for the right L_p-matrix.
     *
     * \param[in] lumped_velo_mass_vec
     * The inverse lumped velocity mass matrix vector
     *
     * \returns
     * A shared pointer to a new BFBT object.
     */
    template <typename MatrixA_, typename MatrixB_, typename MatrixD_, typename FilterV_, typename FilterP_>
    inline std::shared_ptr<BFBT<MatrixA_, MatrixB_, MatrixD_, FilterV_, FilterP_>> new_bfbt(
      const MatrixA_& matrix_a, const MatrixB_& matrix_b,
      const MatrixD_& matrix_d, const FilterV_& filter_v,
      const FilterP_& filter_p,
      std::shared_ptr<SolverBase<typename MatrixD_::VectorTypeL>> solver_left,
      std::shared_ptr<SolverBase<typename MatrixD_::VectorTypeL>> solver_right,
      const typename MatrixB_::VectorTypeL& lumped_velo_mass_vec)
    {
      return std::make_shared<BFBT<MatrixA_, MatrixB_, MatrixD_, FilterV_, FilterP_>>
        (matrix_a, matrix_b, matrix_d, filter_v, filter_p, solver_left, solver_right, lumped_velo_mass_vec);
    }
  } // namespace Solver
} // namespace FEAT
