#pragma once
#ifndef KERNEL_SOLVER_HYPRE_HPP
#define KERNEL_SOLVER_HYPRE_HPP 1

#if defined(FEAT_HAVE_MPI) && defined(FEAT_HAVE_HYPRE)

// includes, FEAT
#include <kernel/solver/adp_solver_base.hpp>

// includes, system
#include <vector>

// includes, HYPRE third-party lib
#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Base-Class for solvers/preconditioners borrowed from HYPRE library
     *
     * This class acts as a base-class for all solvers that we borrow from the HYPRE library.
     * As all our HYPRE solvers use the "IJ"-based ParCSR format, which corresponds to our
     * "algebraic DOF partitioning", this class derives from the ADPSolverBase class, which
     * takes care of the translation between the system matrix and the ADP data structures.
     *
     * This base-class takes care of allocating, initialising and updating the required
     * HYPRE matrix and vector objects.
     *
     * \author Peter Zajac
     */
    template<typename Matrix_, typename Filter_, typename SolverBase_ = Solver::SolverBase<typename Matrix_::VectorTypeL>>
    class HypreSolverBase :
      public ADPSolverBase<Matrix_, Filter_, SolverBase_>
    {
    public:
      /// our base-class
      typedef ADPSolverBase<Matrix_, Filter_, SolverBase_> BaseClass;
      // the vector type
      typedef typename BaseClass::VectorType VectorType;

    protected:
      /// HYPRE global DOF indices array
      std::vector<HYPRE_Int> _hypre_dof_idx;
      /// HYPRE number of non-zero entries per row array
      std::vector<HYPRE_Int> _hypre_num_nze;
      /// HYPRE column indices array
      std::vector<HYPRE_Int> _hypre_col_idx;

      // note: the following members are actually pointers
      /// handle to the HYPRE IJ matrix object
      HYPRE_IJMatrix _hypre_matrix;
      /// handles to the HYPRE IJ vector objects
      HYPRE_IJVector _hypre_vec_def, _hypre_vec_cor;
      /// HYPRE ParCSRMatrix handle to the system matrix
      HYPRE_ParCSRMatrix _hypre_parcsr;
      /// HYPRE ParVector handles to the defect/correction vectors
      HYPRE_ParVector _hypre_par_def, _hypre_par_cor;

      explicit HypreSolverBase(const Matrix_& matrix, const Filter_& filter) :
        BaseClass(matrix, filter),
        _hypre_dof_idx(),
        _hypre_num_nze(),
        _hypre_col_idx(),
        _hypre_matrix(nullptr),
        _hypre_vec_def(nullptr),
        _hypre_vec_cor(nullptr),
        _hypre_parcsr(nullptr),
        _hypre_par_def(nullptr),
        _hypre_par_cor(nullptr)
      {
      }

      /**
       * \brief Uploads the HYPRE defect vector
       *
       * This function first uploads the given defect vector into its ADP defect vector
       * counterpart and afterwards uploads that into the HYPRE vector counterpart.
       *
       * \param[in] vec_def
       * The defect vector to be uploaded from.
       */
      void _upload_def(const VectorType& vec_def)
      {
        // upload defect to ADP vector
        this->_upload_vec_def(vec_def);

        // set HYPRE defect vector values
        HYPRE_IJVectorSetValues(_hypre_vec_def, HYPRE_Int(_hypre_dof_idx.size()), _hypre_dof_idx.data(),
          this->_get_vec_def_vals(vec_def));
      }

      /**
       * \brief Format the HYPRE correction vector
       */
      void _format_cor()
      {
        // format HYPRE correction vector
        HYPRE_ParVectorSetConstantValues(_hypre_par_cor, HYPRE_Real(0));
      }

      /**
       * \brief Downloads the HYPRE correction vector
       *
       * This function first downloads the HYPRE vector into its ADP correction vector
       * counterpart and afterwards downloads that into the given correction vector.
       *
       * \param[out] vec_cor
       * The correction vector to download to.
       */
      void _download_cor(VectorType& vec_cor)
      {
        // get HYPRE correction vector values
        HYPRE_IJVectorGetValues(_hypre_vec_cor, HYPRE_Int(_hypre_dof_idx.size()), _hypre_dof_idx.data(),
          this->_get_vec_cor_vals(vec_cor));

        // download correction from APD vector
        this->_download_vec_cor(vec_cor);

        // apply correction filter
        this->_system_filter.filter_cor(vec_cor);
      }

    public:
      /**
       * \brief Symbolic Initialisation
       *
       * This function creates the HYPRE matrix and vector objects and initialises their
       * structure/layout by using the algebraic DOF partitioning (ADP) that is managed by the
       * base-class. This function also performs an initial upload of the matrix and vector
       * values from the ADP structures (because HYPRE requires this), although these values
       * may be undefined (but existent) at this point.
       */
      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();

        // get our gate and our comm
        const Dist::Comm& comm = *this->_get_comm();

        // get dof offset and count
        const HYPRE_Int dof_offset = HYPRE_Int(this->_get_global_dof_offset());
        const HYPRE_Int num_owned  = HYPRE_Int(this->_get_num_owned_dofs());

        // set up HYPRE dof indices vector
        _hypre_dof_idx.resize((std::size_t)num_owned);
        for(HYPRE_Int i(0); i < num_owned; ++i)
          _hypre_dof_idx[std::size_t(i)] = dof_offset + i;

        // get lower and upper DOF bounds
        const HYPRE_Int ilower = _hypre_dof_idx.front();
        const HYPRE_Int iupper = _hypre_dof_idx.back();

        // create HYPRE matrix
        HYPRE_IJMatrixCreate(comm.mpi_comm(), ilower, iupper, ilower, iupper, &_hypre_matrix);
        HYPRE_IJMatrixSetObjectType(_hypre_matrix, HYPRE_PARCSR);
        HYPRE_IJMatrixInitialize(_hypre_matrix);

        // create HYPRE vectors
        HYPRE_IJVectorCreate(comm.mpi_comm(), ilower, iupper, &_hypre_vec_def);
        HYPRE_IJVectorCreate(comm.mpi_comm(), ilower, iupper, &_hypre_vec_cor);
        HYPRE_IJVectorSetObjectType(_hypre_vec_def, HYPRE_PARCSR);
        HYPRE_IJVectorSetObjectType(_hypre_vec_cor, HYPRE_PARCSR);
        HYPRE_IJVectorInitialize(_hypre_vec_def);
        HYPRE_IJVectorInitialize(_hypre_vec_cor);

        // get matrix dimensions
        const HYPRE_Int num_nze = HYPRE_Int(this->_get_mat_num_nze());

        // get matrix structure arrays
        const auto* row_ptr_a = this->_get_mat_row_ptr();
        const auto* col_idx_a = this->_get_mat_col_idx();

        // create HYPRE index vectors
        _hypre_num_nze.resize((std::size_t)num_owned);
        _hypre_col_idx.resize((std::size_t)num_nze);

        // loop over all owned matrix rows
        for(HYPRE_Int i(0); i < num_owned; ++i)
        {
          _hypre_num_nze[std::size_t(i)] = HYPRE_Int(row_ptr_a[i+1] - row_ptr_a[i]);
        }
        for(HYPRE_Int i(0); i < num_nze; ++i)
        {
          _hypre_col_idx[std::size_t(i)] = HYPRE_Int(col_idx_a[i]);
        }

        // vector for initial vector values
        std::vector<HYPRE_Real> vv((std::size_t)num_owned, 0.0);

        // set matrix structure + dummy values
        HYPRE_IJMatrixSetValues(_hypre_matrix, num_owned, _hypre_num_nze.data(), _hypre_dof_idx.data(),
          _hypre_col_idx.data(), this->_get_mat_vals());

        // set vector structures + dummy values
        HYPRE_IJVectorSetValues(_hypre_vec_def, num_owned, _hypre_dof_idx.data(), vv.data());
        HYPRE_IJVectorSetValues(_hypre_vec_cor, num_owned, _hypre_dof_idx.data(), vv.data());

        // assemble matrix and get ParCSR object pointer
        HYPRE_IJMatrixAssemble(_hypre_matrix);
        HYPRE_IJMatrixGetObject(_hypre_matrix, (void**) &_hypre_parcsr);

        // assemble vectors and get Par object pointers
        HYPRE_IJVectorAssemble(_hypre_vec_def);
        HYPRE_IJVectorAssemble(_hypre_vec_cor);
        HYPRE_IJVectorGetObject(_hypre_vec_def, (void **) &_hypre_par_def);
        HYPRE_IJVectorGetObject(_hypre_vec_cor, (void **) &_hypre_par_cor);
      }

      /**
       * \brief Numeric Initialisation
       *
       * This function uploads the numerical values of the ADP matrix, which
       * is managed by the base-class, to the HYPRE matrix.
       */
      virtual void init_numeric() override
      {
        BaseClass::init_numeric();

        // update matrix values of HYPRE matrix
        HYPRE_IJMatrixSetValues(_hypre_matrix, HYPRE_Int(_hypre_dof_idx.size()), _hypre_num_nze.data(),
          _hypre_dof_idx.data(), _hypre_col_idx.data(), this->_get_mat_vals());
      }

      /**
       * \brief Symbolic Finalisation
       *
       * This function destroys all HYPRE objects managed by this class
       * and resets all auxiliary vectors and pointers.
       */
      virtual void done_symbolic() override
      {
        // these are just pointers to HYPRE's internal objects
        _hypre_par_cor = nullptr;
        _hypre_par_def = nullptr;
        _hypre_parcsr = nullptr;

        // destroy HYPRE objects
        HYPRE_IJVectorDestroy(_hypre_vec_cor);
        HYPRE_IJVectorDestroy(_hypre_vec_def);
        HYPRE_IJMatrixDestroy(_hypre_matrix);

        // reset object handles
        _hypre_vec_cor = nullptr;
        _hypre_vec_def = nullptr;
        _hypre_matrix = nullptr;

        // clear index arrays
        _hypre_col_idx.clear();
        _hypre_num_nze.clear();
        _hypre_dof_idx.clear();

        BaseClass::done_symbolic();
      }
    }; // class HypreSolverBase<...>

    /**
     * \brief HYPRE ParaSails Preconditioner Wrapper class template
     *
     * This class acts as a wrapper around the ParaSails preconditioner from the HYPRE library.
     * ParaSails is some sort of parallel sparse approximate inverse preconditioner, see the
     * documentation of HYPRE for details.
     *
     * \todo support setting of solver parameters
     *
     * \author Peter Zajac
     */
    template<typename Matrix_, typename Filter_>
    class ParaSailsPrecond :
      public HypreSolverBase<Matrix_, Filter_>
    {
    public:
      /// our base-class
      typedef HypreSolverBase<Matrix_, Filter_> BaseClass;
      /// the vector type
      typedef typename BaseClass::VectorType VectorType;

    protected:
      /// the HYPRE solver object
      HYPRE_Solver _precond;

    public:
      explicit ParaSailsPrecond(const Matrix_& matrix, const Filter_& filter) :
        BaseClass(matrix, filter),
        _precond(nullptr)
      {
      }

      explicit ParaSailsPrecond(const String& DOXY(section_name), PropertyMap* DOXY(section),
        const Matrix_& matrix, const Filter_& filter)
         :
        ParaSailsPrecond(matrix, filter)
      {
        /// \todo parse parameters from property map
      }

      virtual String name() const override
      {
        return "ParaSailsPrecond";
      }

      virtual void init_numeric() override
      {
        BaseClass::init_numeric();

        // create ParaSails preconditioner
        HYPRE_ParaSailsCreate(this->_get_comm()->mpi_comm(), &_precond);

        /// \todo set ParaSails parameters
        HYPRE_ParaSailsSetParams(_precond, 0.1, 1);
        HYPRE_ParaSailsSetFilter(_precond, 0.05);
        HYPRE_ParaSailsSetSym(_precond, 2);
        //HYPRE_ParaSailsSetLogging(_precond, 3);

        // setup ParaSails
        // according to the documentation, the two vectors are ignored by this function
        HYPRE_ParaSailsSetup(this->_precond, this->_hypre_parcsr, this->_hypre_par_def, this->_hypre_par_cor);
      }

      virtual void done_numeric() override
      {
        if(_precond != nullptr)
        {
          HYPRE_ParaSailsDestroy(_precond);
          _precond = nullptr;
        }

        BaseClass::done_numeric();
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        // upload defect vector and format correction
        this->_upload_def(vec_def);
        this->_format_cor();

        // apply ParaSails preconditioner
        HYPRE_ParaSailsSolve(this->_precond, this->_hypre_parcsr, this->_hypre_par_def, this->_hypre_par_cor);

        // download correction
        this->_download_cor(vec_cor);

        // okay
        return Status::success;
      }
    }; // class ParaSailsPrecond

    /**
     * \brief Creates a new ParaSailsPrecond solver object
     *
     * \param[in] matrix
     * The global system matrix.
     *
     * \param[in] filter
     * The global system filter.
     *
     * \returns
     * A shared pointer to a new ParaSailsPrecond object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<ParaSailsPrecond<Matrix_, Filter_>> new_parasails_precond(
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<ParaSailsPrecond<Matrix_, Filter_>>(matrix, filter);
    }

    /**
     * \brief Creates a new ParaSailsPrecond solver object based on a PropertyMap
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
     * \returns
     * A shared pointer to a new ParaSailsPrecond object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<ParaSailsPrecond<Matrix_, Filter_>> new_parasails_precond(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<ParaSailsPrecond<Matrix_, Filter_>>(section_name, section, matrix, filter);
    }

    /**
     * \brief HYPRE Euclid Preconditioner Wrapper class template
     *
     * This class acts as a wrapper around the Euclid preconditioner from the HYPRE library.
     * Euclid is some sort of parallel ILU preconditioner, see the documentation of HYPRE
     * for details.
     *
     * \todo support setting of solver parameters
     *
     * \author Peter Zajac
     */
    template<typename Matrix_, typename Filter_>
    class EuclidPrecond :
      public HypreSolverBase<Matrix_, Filter_>
    {
    public:
      /// our base-class
      typedef HypreSolverBase<Matrix_, Filter_> BaseClass;
      /// the vector type
      typedef typename BaseClass::VectorType VectorType;

    protected:
      /// the HYPRE solver object
      HYPRE_Solver _precond;

    public:
      explicit EuclidPrecond(const Matrix_& matrix, const Filter_& filter) :
        BaseClass(matrix, filter),
        _precond(nullptr)
      {
      }

      explicit EuclidPrecond(const String& DOXY(section_name), PropertyMap* DOXY(section),
        const Matrix_& matrix, const Filter_& filter)
         :
        EuclidPrecond(matrix, filter)
      {
        /// \todo parse parameters from property map
      }

      virtual String name() const override
      {
        return "EuclidPrecond";
      }

      virtual void init_numeric() override
      {
        BaseClass::init_numeric();

        // create Euclid preconditioner
        HYPRE_EuclidCreate(this->_get_comm()->mpi_comm(), &_precond);

        /// \todo set Euclid parameters

        // setup Euclid
        // according to the documentation, the two vectors are ignored by this function
        HYPRE_EuclidSetup(this->_precond, this->_hypre_parcsr, this->_hypre_par_def, this->_hypre_par_cor);
      }

      virtual void done_numeric() override
      {
        if(_precond != nullptr)
        {
          HYPRE_EuclidDestroy(_precond);
          _precond = nullptr;
        }
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        // upload defect vector and format correction
        this->_upload_def(vec_def);
        this->_format_cor();

        // apply Euclid preconditioner
        HYPRE_EuclidSolve(this->_precond, this->_hypre_parcsr, this->_hypre_par_def, this->_hypre_par_cor);

        // download correction
        this->_download_cor(vec_cor);

        // okay
        return Status::success;
      }
    }; // class EuclidPrecond

    /**
     * \brief Creates a new EuclidPrecond solver object
     *
     * \param[in] matrix
     * The global system matrix.
     *
     * \param[in] filter
     * The global system filter.
     *
     * \returns
     * A shared pointer to a new EuclidPrecond object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<EuclidPrecond<Matrix_, Filter_>> new_euclid_precond(
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<EuclidPrecond<Matrix_, Filter_>>(matrix, filter);
    }

    /**
     * \brief Creates a new EuclidPrecond solver object based on a PropertyMap
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
     * \returns
     * A shared pointer to a new EuclidPrecond object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<EuclidPrecond<Matrix_, Filter_>> new_euclid_precond(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<EuclidPrecond<Matrix_, Filter_>>(section_name, section, matrix, filter);
    }

    /**
     * \brief HYPRE BoomerAMGWrapper class template
     *
     * This class acts as a wrapper around the BoomerAMG solver/preconditioner from the HYPRE library.
     *
     * \todo support setting of solver parameters
     *
     * \author Peter Zajac
     */
    template<typename Matrix_, typename Filter_>
    class BoomerAMG :
      public HypreSolverBase<Matrix_, Filter_>
    {
    public:
      /// our base-class
      typedef HypreSolverBase<Matrix_, Filter_> BaseClass;
      /// the vector type
      typedef typename BaseClass::VectorType VectorType;

    protected:
      /// the HYPRE solver object
      HYPRE_Solver _precond;

    public:
      explicit BoomerAMG(const Matrix_& matrix, const Filter_& filter) :
        BaseClass(matrix, filter),
        _precond(nullptr)
      {
      }

      explicit BoomerAMG(const String& DOXY(section_name), PropertyMap* DOXY(section),
        const Matrix_& matrix, const Filter_& filter)
         :
        BoomerAMG(matrix, filter)
      {
        /// \todo parse parameters
      }

      virtual String name() const override
      {
        return "BoomerAMG";
      }

      virtual void init_numeric() override
      {
        BaseClass::init_numeric();

        // create BoomerAMG preconditioner
        HYPRE_BoomerAMGCreate(&_precond);

        /// \todo set BoomerAMG parameters

        // set tolerance to 0 and maximum iterations to 1, so that BoomerAMG
        // acts as a preconditioner rather than a real iterative solver
        HYPRE_BoomerAMGSetTol(_precond, 0.0);
        HYPRE_BoomerAMGSetMaxIter(_precond, 1);

        // setup BoomerAMG
        // according to the documentation, the two vectors are ignored by this function
        HYPRE_BoomerAMGSetup(this->_precond, this->_hypre_parcsr, this->_hypre_par_def, this->_hypre_par_cor);
      }

      virtual void done_numeric() override
      {
        if(_precond != nullptr)
        {
          HYPRE_BoomerAMGDestroy(_precond);
          _precond = nullptr;
        }
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        // upload defect vector and format correction
        this->_upload_def(vec_def);
        this->_format_cor();

        // apply BoomerAMG preconditioner
        HYPRE_BoomerAMGSolve(this->_precond, this->_hypre_parcsr, this->_hypre_par_def, this->_hypre_par_cor);

        // download correction
        this->_download_cor(vec_cor);

        // okay
        return Status::success;
      }
    }; // class BoomerAMG

    /**
     * \brief Creates a new BoomerAMG solver object
     *
     * \param[in] matrix
     * The global system matrix.
     *
     * \param[in] filter
     * The global system filter.
     *
     * \returns
     * A shared pointer to a new BoomerAMG object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<BoomerAMG<Matrix_, Filter_>> new_boomeramg(
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<BoomerAMG<Matrix_, Filter_>>(matrix, filter);
    }

    /**
     * \brief Creates a new BoomerAMG solver object based on a PropertyMap
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
     * \returns
     * A shared pointer to a new BoomerAMG object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<BoomerAMG<Matrix_, Filter_>> new_boomeramg(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<BoomerAMG<Matrix_, Filter_>>(section_name, section, matrix, filter);
    }
  } // namespace Solver
} // namespace FEAT

#endif // defined(FEAT_HAVE_MPI) && defined(FEAT_HAVE_HYPRE)
#endif // KERNEL_SOLVER_HYPRE_HPP
