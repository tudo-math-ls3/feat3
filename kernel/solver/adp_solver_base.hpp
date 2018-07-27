#pragma once
#ifndef KERNEL_SOLVER_ADP_SOLVER_BASE_HPP
#define KERNEL_SOLVER_ADP_SOLVER_BASE_HPP 1

// includes, FEAT
#include <kernel/solver/base.hpp>
#include <kernel/global/alg_dof_parti.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/filter.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/global/pmdcdsc_matrix.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Base-Class for Solvers based on Algebraic-DOF-Partitioning
     *
     * This class template specifies a common interface for solvers based on algebraic DOF partitioning.
     * There are two specialisations of this class, which depend on the type of the system matrix.
     *
     * \attention
     * This class is only defined for this doxygen documentation! There exists no generic implementation
     * of this class template, but only the two specialisations for Global::Matrix and Global::PMDCDSCMatrix,
     * however, both these specialisations implement the same protected interface, which is defined in
     * this generic template.
     *
     * \tparam Matrix_
     * The type of the system matrix that is to be used by the solver.
     *
     * \tparam Filter_
     * The type of the system filter.
     *
     * \tparam SolverBase_
     * Specifies the base-class that this class should derive from. Defaults to Solver::SolverBase, but
     * it may also be set to Solver::IterativeSolver if the derived class represents an iterative solver.
     *
     * \author Peter Zajac
     */
    template<typename Matrix_, typename Filter_, typename SolverBase_ = Solver::SolverBase<typename Matrix_::VectorTypeL>>
#ifndef DOXYGEN
    class ADPSolverBase;
#else // DOXYGEN
    class ADPSolverBase :
      public SolverBase_
    {
    protected:
      /// constructor
      explicit ADPSolverBase(const GlobalMatrixType& matrix, const GlobalFilterType& filter);

      /// \returns A pointer to the underlying communicator
      const Dist::Comm* _get_comm() const;

      /// \returns The number of DOFs owned by this process.
      Index _get_num_owned_dofs() const;

      /// \returns The total number of DOFs over all processes.
      Index _get_num_global_dofs() const;

      /// \returns The offset of the first global DOF index of this process.
      Index _get_global_dof_offset() const;

      /// \returns The number of non-zero entries in this process's matrix partition
      Index _get_mat_num_nze() const;

      /// \returns The value array of the matrix partition.
      const DataType* _get_mat_vals() const;

      /// \returns The row-pointer array of the matrix partition.
      const IndexType* _get_mat_row_ptr() const;

      /// \returns The column-index array of the matrix partition.
      const IndexType* _get_mat_col_idx() const;

      /// \returns The data array of the owned defect vector.
      const DataType* _get_vec_def_vals(const VectorType& vec_def);

      /// \returns The data array of the owned defect vector.
      DataType* _get_vec_def_vals(VectorType& vec_def);

      /// \returns The data array of the owned correction vector.
      DataType* _get_vec_cor_vals(VectorType& vec_cor);

      /// \returns Uploads a defect vector.
      void _upload_vec_def(const VectorType& vec_def);

      /// \returns Uploads a correction vector.
      void _upload_vec_cor(const VectorType& vec_cor);

      /// \returns Downloads a defect vector.
      void _download_vec_def(VectorType& vec_def);

      /// \returns Downloads a correction vector.
      void _download_vec_cor(VectorType& vec_cor);
    }; // class ADPSolverBase
#endif // DOXYGEN

    /**
     * \brief Implementation of ADPSolverBase interface for Global::Matrix instances
     *
     * This class implements the ADPSolverBase interface for Global::Matrix instances, which derive
     * from standard (scalar) finite element discretisations. Most of the dirty work is outsourced
     * to the Global::AlgDofParti, Global::AlgDofPartiVector and Global::AlgDoPartiMatrix classes.
     *
     * \note
     * Although this class does not explicitly expect the local matrix to be a LAFEM::SparseMatrixCSR,
     * this is (currently) the only matrix type that is supported by the underlying Global::AlgDofParti
     * class template.
     *
     * \author Peter Zajac
     */
    template<typename LocalMatrix_, typename Mirror_, typename GlobalFilter_, typename SolverBase_>
    class ADPSolverBase<Global::Matrix<LocalMatrix_, Mirror_, Mirror_>, GlobalFilter_, SolverBase_> :
      //public SolverBase<Global::Vector<typename LocalMatrix_::VectorTypeL, Mirror_>>
      public SolverBase_
    {
    public:
      typedef LocalMatrix_ LocalMatrixType;
      typedef typename LocalMatrixType::VectorTypeL LocalVectorType;
      typedef Mirror_ MirrorType;

      typedef Global::Matrix<LocalMatrixType, MirrorType, MirrorType> GlobalMatrixType;
      typedef Global::Vector<LocalVectorType, MirrorType> GlobalVectorType;
      typedef GlobalFilter_ GlobalFilterType;

      //typedef SolverBase<GlobalVectorType> BaseClass;
      typedef SolverBase_ BaseClass;

      typedef GlobalVectorType VectorType;

      typedef Global::AlgDofParti<LocalVectorType, MirrorType> AlgDofPartiType;
      typedef Global::AlgDofPartiMatrix<LocalMatrixType, MirrorType> ADPMatrixType;
      typedef Global::AlgDofPartiVector<LocalVectorType, MirrorType> ADPVectorType;

      typedef typename AlgDofPartiType::DataType DataType;
      typedef typename AlgDofPartiType::IndexType IndexType;

    private:
      /// the algebraic DOF partitioning
      AlgDofPartiType _alg_dof_parti;
      /// the ADP matrix
      ADPMatrixType _adp_matrix;
      /// two ADP vectors for correction and defect vectors
      ADPVectorType _adp_vec_def, _adp_vec_cor;

    protected:
      /// the global system matrix
      const GlobalMatrixType& _system_matrix;
      /// the global system filter
      const GlobalFilterType& _system_filter;

      explicit ADPSolverBase(const GlobalMatrixType& matrix, const GlobalFilterType& filter) :
        BaseClass(),
        _alg_dof_parti(),
        _adp_matrix(),
        _adp_vec_def(),
        _adp_vec_cor(),
        _system_matrix(matrix),
        _system_filter(filter)
      {
      }

    public:
      /// \returns A pointer to the underlying communicator
      const Dist::Comm* _get_comm() const
      {
        return this->_system_matrix.get_comm();
      }

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
        _alg_dof_parti.assemble_by_gate(*this->_system_matrix.get_row_gate());

        // upload the matrix structure
        _adp_matrix = ADPMatrixType(&this->_alg_dof_parti);
        _adp_matrix.upload_symbolic(this->_system_matrix);

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
        _adp_matrix.upload_numeric(this->_system_matrix);

        // filter matrix
        _adp_matrix.filter_matrix(this->_system_filter);
      }

      virtual void done_symbolic() override
      {
        _adp_vec_cor.clear();
        _adp_vec_def.clear();
        _adp_matrix.clear();
        _alg_dof_parti.clear();

        BaseClass::done_symbolic();
      }

    protected:
      Index _get_num_owned_dofs() const
      {
        return this->_alg_dof_parti.get_num_owned_dofs();
      }

      Index _get_num_global_dofs() const
      {
        return this->_alg_dof_parti.get_num_global_dofs();
      }

      Index _get_global_dof_offset() const
      {
        return this->_alg_dof_parti.get_global_dof_offset();
      }

      Index _get_mat_num_nze() const
      {
        return this->_adp_matrix.owned().used_elements();
      }

      const DataType* _get_mat_vals() const
      {
        return this->_adp_matrix.owned().val();
      }

      const IndexType* _get_mat_row_ptr() const
      {
        return this->_adp_matrix.owned().row_ptr();
      }

      const IndexType* _get_mat_col_idx() const
      {
        return this->_adp_matrix.owned().col_ind();
      }

      const DataType* _get_vec_def_vals(const VectorType&)
      {
        return this->_adp_vec_def.owned().elements();
      }

      DataType* _get_vec_def_vals(VectorType&)
      {
        return this->_adp_vec_def.owned().elements();
      }

      DataType* _get_vec_cor_vals(VectorType&)
      {
        return this->_adp_vec_cor.owned().elements();
      }

      void _upload_vec_def(const VectorType& vec_def)
      {
        this->_adp_vec_def.upload(vec_def);
      }

      void _upload_vec_cor(const VectorType& vec_cor)
      {
        this->_adp_vec_cor.upload(vec_cor);
      }

      void _download_vec_def(VectorType& vec_def)
      {
        this->_adp_vec_def.download(vec_def);
      }

      void _download_vec_cor(VectorType& vec_cor)
      {
        this->_adp_vec_cor.download(vec_cor);
      }
    }; // class ADPSolverBase<Global::Matrix<...>, ...>


    /**
     * \brief Implementation of ADPSolverBase interface for Global::PMDCDSCMatrix instances
     *
     * This class implements the ADPSolverBase interface for Global::PMDCDSCMatrix instances,
     * which are used to implement the global Schur-complement matrices appearing in Navier-Stokes
     * equations with discontinuous pressure spaces. The most prominent use of this specialisation
     * is the Pressure-Poisson-Problem solver of the infamous "PP" application.
     *
     * Most of the dirty work is outsourced to the Global::PMDCDSCMatrix class, most notably to
     * its two member functions Global::PMDCDSCMatrix::asm_adp_symbolic() and
     * Global::PMDCDSCMatrix::asm_adp_numeric().
     *
     * \author Peter Zajac
     */
    template<typename MatrixB_, typename MatrixD_, typename GlobalFilter_, typename SolverBase_>
    class ADPSolverBase<Global::PMDCDSCMatrix<MatrixB_, MatrixD_>, GlobalFilter_, SolverBase_> :
      //public SolverBase<typename Global::PMDCDSCMatrix<MatrixB_, MatrixD_>::VectorTypeL>
      public SolverBase_
    {
    public:
      typedef Global::PMDCDSCMatrix<MatrixB_, MatrixD_> GlobalMatrixType;
      typedef GlobalFilter_ GlobalFilterType;

      typedef typename GlobalMatrixType::LocalMatrixTypeS LocalMatrixTypeS;

      typedef typename GlobalMatrixType::VectorTypeL GlobalVectorType;

      //typedef SolverBase<GlobalVectorType> BaseClass;
      typedef SolverBase_ BaseClass;

      typedef GlobalVectorType VectorType;

      typedef typename LocalMatrixTypeS::DataType DataType;
      typedef typename LocalMatrixTypeS::IndexType IndexType;

    private:
      /// our global DOF offset
      Index _glob_dof_offset;
      // the total global DOF count
      Index _glob_dof_count;
      /// our Schur-complement matrix slice
      LocalMatrixTypeS _matrix_slice;

    protected:
      /// the global system matrix
      const GlobalMatrixType& _system_matrix;
      /// the global system filter
      const GlobalFilterType& _system_filter;

      explicit ADPSolverBase(const GlobalMatrixType& matrix, const GlobalFilterType& filter) :
        BaseClass(),
        _glob_dof_offset(0),
        _glob_dof_count(0),
        _system_matrix(matrix),
        _system_filter(filter)
      {
      }

    public:
      /// \returns A pointer to the underlying communicator
      const Dist::Comm* _get_comm() const
      {
        return this->_system_matrix.get_comm();
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();

        // assemble matrix slice structure
        _matrix_slice = _system_matrix.asm_adp_symbolic(_glob_dof_offset, _glob_dof_count);
      }

      virtual void init_numeric() override
      {
        BaseClass::init_numeric();

        // upload matrix values
        _system_matrix.asm_adp_numeric(_matrix_slice);

        // filter matrix
        this->_system_filter.local().filter_mat(_matrix_slice);
      }

      virtual void done_symbolic() override
       {
        _matrix_slice.clear();
        BaseClass::done_symbolic();
      }

    protected:

      Index _get_num_global_dofs() const
      {
        return _glob_dof_count();
      }

      Index _get_global_dof_offset() const
      {
        return _glob_dof_offset;
      }

      Index _get_num_owned_dofs() const
      {
        return _matrix_slice.rows();
      }

      Index _get_mat_num_nze() const
      {
        return _matrix_slice.used_elements();
      }

      const DataType* _get_mat_vals() const
      {
        return _matrix_slice.val();
      }

      const IndexType* _get_mat_row_ptr() const
      {
        return _matrix_slice.row_ptr();
      }

      const IndexType* _get_mat_col_idx() const
      {
        return _matrix_slice.col_ind();
      }

      const DataType* _get_vec_def_vals(const VectorType& vec_def)
      {
        return vec_def.local().elements();
      }

      DataType* _get_vec_def_vals(VectorType& vec_def)
      {
        return vec_def.local().elements();
      }

      DataType* _get_vec_cor_vals(VectorType& vec_cor)
      {
        return vec_cor.local().elements();
      }

      void _upload_vec_def(const VectorType&)
      {
        // nothing to do here
      }

      void _upload_vec_cor(const VectorType&)
      {
        // nothing to do here
      }

      void _download_vec_def(VectorType&)
      {
        // nothing to do here
      }

      void _download_vec_cor(VectorType&)
      {
        // nothing to do here
      }
    }; // class ADPSolverBase<Global::PMDCDSCMatrix<...>, ...>
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_ADP_SOLVER_BASE_HPP
