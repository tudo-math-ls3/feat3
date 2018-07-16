#pragma once
#ifndef KERNEL_SOLVER_ADP_SOLVER_BASE_HPP
#define KERNEL_SOLVER_ADP_SOLVER_BASE_HPP 1

// includes, FEAT
#include <kernel/solver/base.hpp>
#include <kernel/global/alg_dof_parti.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/filter.hpp>
#include <kernel/global/vector.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Base-Class for Solvers based on Algebraic-DOF-Partitioning
     *
     * See the documentation of the specialisation for more details.
     *
     * \author Peter Zajac
     */
    template<typename Matrix_, typename Filter_>
    class ADPSolverBase;

    /**
     * \brief Base-Class for Solvers based on Algebraic-DOF-Partitioning
     *
     * See the documentation of the specialisation for more details.
     *
     * \author Peter Zajac
     */
    template<typename LocalMatrix_, typename Mirror_, typename GlobalFilter_>
    class ADPSolverBase<Global::Matrix<LocalMatrix_, Mirror_, Mirror_>, GlobalFilter_> :
      public SolverBase<Global::Vector<typename LocalMatrix_::VectorTypeL, Mirror_>>
    {
    public:
      typedef LocalMatrix_ LocalMatrixType;
      typedef typename LocalMatrixType::VectorTypeL LocalVectorType;
      typedef Mirror_ MirrorType;

      typedef Global::Matrix<LocalMatrixType, MirrorType, MirrorType> GlobalMatrixType;
      typedef Global::Vector<LocalVectorType, MirrorType> GlobalVectorType;
      typedef GlobalFilter_ GlobalFilterType;

      typedef SolverBase<GlobalVectorType> BaseClass;

      typedef GlobalVectorType VectorType;

      typedef Global::AlgDofParti<LocalVectorType, MirrorType> AlgDofPartiType;
      typedef Global::AlgDofPartiMatrix<LocalMatrixType, MirrorType> ADPMatrixType;
      typedef Global::AlgDofPartiVector<LocalVectorType, MirrorType> ADPVectorType;

    protected:
      /// the global system matrix
      const GlobalMatrixType& _matrix;
      /// the global system filter
      const GlobalFilterType& _filter;
      /// the algebraic DOF partitioning
      AlgDofPartiType _alg_dof_parti;
      /// the ADP matrix
      ADPMatrixType _adp_matrix;
      /// two ADP vectors for correction and defect vectors
      ADPVectorType _adp_vec_def, _adp_vec_cor;

      explicit ADPSolverBase(const GlobalMatrixType& matrix, const GlobalFilterType& filter) :
        BaseClass(),
        _matrix(matrix),
        _filter(filter),
        _alg_dof_parti(),
        _adp_matrix(),
        _adp_vec_def(),
        _adp_vec_cor()
      {
      }

      /// \returns A pointer to the underlying communicator
      const Dist::Comm* _get_comm() const
      {
        return this->_matrix.get_comm();
      }

    public:
      /**
       * \brief Symbolic Initialisation
       *
       * This function assembles the algebraic DOF partitioning, creates the ADP matrix,
       * uploads the matrix structure and creates the two ADP vectors.
       */
      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();

        // assemble the algebraic dof partitioning
        _alg_dof_parti.assemble_by_gate(*this->_matrix.get_row_gate());

        // upload the matrix structure
        _adp_matrix = ADPMatrixType(&this->_alg_dof_parti);
        _adp_matrix.upload_symbolic(this->_matrix);

        // create two vectors
        _adp_vec_def = ADPVectorType(&this->_alg_dof_parti);
        _adp_vec_cor = ADPVectorType(&this->_alg_dof_parti);

        // format matrix and vectors
        _adp_matrix.owned().format();
        _adp_vec_def.owned().format();
        _adp_vec_cor.owned().format();
      }

      /**
       * \brief Numeric Initialisation
       *
       * This function uploads the numeric matrix values and
       * applies the filter onto the ADP matrix.
       */
      virtual void init_numeric() override
      {
        BaseClass::init_numeric();

        // upload matrix values
        _adp_matrix.upload_numeric(this->_matrix);

        // filter matrix
        _adp_matrix.filter_matrix(this->_filter);
      }

      virtual void done_symbolic() override
      {
        _adp_vec_cor.clear();
        _adp_vec_def.clear();
        _adp_matrix.clear();
        _alg_dof_parti.clear();

        BaseClass::done_symbolic();
      }
    }; // class ADPSolverBase<Global::Matrix<...>, ...>

  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_ADP_SOLVER_BASE_HPP
