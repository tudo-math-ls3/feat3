// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SOLVER_HYPRE_HPP
#define KERNEL_SOLVER_HYPRE_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/solver/adp_solver_base.hpp>

#if defined(FEAT_HAVE_HYPRE) || defined(DOXYGEN)

// includes, system
#include <vector>

namespace FEAT
{
  namespace Solver
  {
    /// \cond internal
    namespace Hypre
    {
      /**
       * \brief Creates a core wrapper object for HYPRE.
       *
       * The "core wrapper object", which is returned by this function, is basically a collection
       * of all HYPRE data structures, which are required to represent a linear system A*x=b, i.e.
       * a partitioned matrix as well as two vectors.
       *
       * \param[in] comm
       * A pointer to an MPI_Comm object that represents the communicator.
       * This argument is ignored in non-MPI builds.
       *
       * \param[in] dof_offset
       * The global DOF offset for this process.
       *
       * \param[in] num_owned_dofs
       * The number of global DOFs owned by this process.
       *
       * \param[in] row_ptr
       * The row-pointer array of the partitioned CSR matrix.
       *
       * \param[in] col_idx
       * The column-index array of the partitioned CSR matrix.
       *
       * \returns
       * A pointer to a newly allocated core wrapper object.
       */
      void* create_core(const void* comm, Index dof_offset, Index num_owned_dofs,
        const unsigned int* row_ptr, const unsigned int* col_idx);
      void* create_core(const void* comm, Index dof_offset, Index num_owned_dofs,
        const unsigned long* row_ptr, const unsigned long* col_idx);
      void* create_core(const void* comm, Index dof_offset, Index num_owned_dofs,
        const unsigned long long* row_ptr, const unsigned long long* col_idx);

      void destroy_core(void* core);

      void set_matrix_values(void* core, const double* vals);
      void set_vec_cor_values(void* core, const double* vals);
      void set_vec_def_values(void* core, const double* vals);
      void get_vec_cor_values(const void* core, double* vals);
      void get_vec_def_values(const void* core, double* vals);

      // Parasails Wrappers
      void* create_parasails(void* core, int* iparam, double* dparam);
      void destroy_parasails(void* solver);
      void solve_parasails(void* core, void* solver);

      // Euclid Wrappers
      void* create_euclid(void* core, int* iparam, double* dparam);
      void destroy_euclid(void* solver);
      void solve_euclid(void* core, void* solver);

      // BoomerAMG Wrappers
      void* create_boomeramg(void* core, int* iparam, double* dparam);
      void destroy_boomeramg(void* solver);
      void solve_boomeramg(void* core, void* solver);
    } // namespace Hypre
    /// \endcond

    /**
     * \brief Base-Class for solvers/preconditioners borrowed from HYPRE library
     *
     * This class acts as a base-class for all solvers that we borrow from the HYPRE library.
     * As all our HYPRE solvers use the "IJ"-based ParCSR format, which corresponds to our
     * "algebraic DOF partitioning", this class derives from the ADPSolverBase class, which
     * takes care of the translation between the system matrix and the ADP data structures.
     *
     * This base-class takes care of allocating, initialising and updating the required
     * HYPRE matrix and vector objects, which are outsourced into an opaque core wrapper object.
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
      /// a pointer to our opaque core wrapper object
      void* _core;

      explicit HypreSolverBase(const Matrix_& matrix, const Filter_& filter) :
        BaseClass(matrix, filter),
        _core(nullptr)
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
        Hypre::set_vec_def_values(this->_core, this->_get_vec_def_vals(vec_def));
      }

      /**
       * \brief Format the HYPRE correction vector
       */
      void _format_cor()
      {
        // format HYPRE correction vector
        Hypre::set_vec_cor_values(this->_core, nullptr);
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
        Hypre::get_vec_cor_values(this->_core, this->_get_vec_cor_vals(vec_cor));

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

        XASSERT(this->_core == nullptr);

        // create our HYPRE core wrapper object
        this->_core = Hypre::create_core(
#ifdef FEAT_HAVE_MPI
          &this->_get_comm()->mpi_comm(),
#else
          nullptr, // no communicator for non-MPI builds
#endif
          this->_get_global_dof_offset(), this->_get_num_owned_dofs(), this->_get_mat_row_ptr(), this->_get_mat_col_idx());
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

        XASSERT(this->_core != nullptr);

        // update matrix values of HYPRE matrix
        Hypre::set_matrix_values(this->_core, this->_get_mat_vals());
      }

      /**
       * \brief Symbolic Finalisation
       *
       * This function destroys all HYPRE objects managed by this class
       * and resets all auxiliary vectors and pointers.
       */
      virtual void done_symbolic() override
      {
        XASSERT(this->_core != nullptr);

        Hypre::destroy_core(this->_core);
        this->_core = nullptr;

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
      void* _solver;

      /// integer parameters:
      // 0: num levels, third argument for HYPRE_ParaSailsSetParams, default: 1
      // 1: symmetry, second argument for HYPRE_ParaSailsSetSym, default: 2
      int _iparam[2];

      /// double parameters:
      // 0: threshold, second argument for HYPRE_ParaSailsSetParams, default: 0.1
      // 1: filter, second argument for HYPRE_ParaSailsSetFilter, default: 0.05
      double _dparam[2];

    public:
      explicit ParaSailsPrecond(const Matrix_& matrix, const Filter_& filter) :
        BaseClass(matrix, filter),
        _solver(nullptr)
      {
        _iparam[0] = 1;
        _iparam[1] = 2;

        _dparam[0] = 0.1;
        _dparam[1] = 0.05;
      }

      explicit ParaSailsPrecond(const String& section_name, PropertyMap* section,
        const Matrix_& matrix, const Filter_& filter)
         :
        ParaSailsPrecond(matrix, filter)
      {
        auto level_p = section->query("level");
        if(level_p.second && !level_p.first.parse(_iparam[0]))
          throw ParseError(section_name + ".level", level_p.first, "a non-negative integer");

        auto sym_p = section->query("sym");
        if(sym_p.second && !sym_p.first.parse(_iparam[1]))
          throw ParseError(section_name + ".sym", sym_p.first, "a non-negative integer");

        auto thresh_p = section->query("thresh");
        if(thresh_p.second && !thresh_p.first.parse(_dparam[0]))
          throw ParseError(section_name + ".thresh", thresh_p.first, "a positive float");

        auto filter_p = section->query("filter");
        if(filter_p.second && !filter_p.first.parse(_dparam[1]))
          throw ParseError(section_name + ".filter", filter_p.first, "a positive float");
      }

      virtual String name() const override
      {
        return "ParaSailsPrecond";
      }

      virtual void init_numeric() override
      {
        BaseClass::init_numeric();

        // create ParaSails preconditioner
        this->_solver = Hypre::create_parasails(this->_core, this->_iparam, this->_dparam);
      }

      virtual void done_numeric() override
      {
        if(_solver != nullptr)
        {
          Hypre::destroy_parasails(_solver);
          _solver = nullptr;
        }

        BaseClass::done_numeric();
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        // upload defect vector and format correction
        this->_upload_def(vec_def);
        this->_format_cor();

        // apply ParaSails preconditioner
        Hypre::solve_parasails(this->_core, this->_solver);

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
      void* _solver;

      /// integer parameters:
      // 0: factorisation level for ILU(k), default: 1
      int _iparam[1];

      /// double parameters:
      // 0: drop tolerance for ILU(k), default: 0.0
      double _dparam[1];

    public:
      explicit EuclidPrecond(const Matrix_& matrix, const Filter_& filter) :
        BaseClass(matrix, filter),
        _solver(nullptr)
      {
        _iparam[0] = 1;
        _dparam[0] = 0.0;
      }

      explicit EuclidPrecond(const String& section_name, PropertyMap* section,
        const Matrix_& matrix, const Filter_& filter)
         :
        EuclidPrecond(matrix, filter)
      {
        auto level_p = section->query("level");
        if(level_p.second && !level_p.first.parse(_iparam[0]))
          throw ParseError(section_name + ".level", level_p.first, "a non-negative integer");

        auto drop_p = section->query("drop");
        if(drop_p.second && !drop_p.first.parse(_dparam[0]))
          throw ParseError(section_name + ".drop", drop_p.first, "a non-negative float");
      }

      virtual String name() const override
      {
        return "EuclidPrecond";
      }

      virtual void init_numeric() override
      {
        BaseClass::init_numeric();

        // create Euclid preconditioner
        this->_solver = Hypre::create_euclid(this->_core, this->_iparam, this->_dparam);
      }

      virtual void done_numeric() override
      {
        if(_solver != nullptr)
        {
          Hypre::destroy_euclid(_solver);
          _solver = nullptr;
        }
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        // upload defect vector and format correction
        this->_upload_def(vec_def);
        this->_format_cor();

        // apply Euclid preconditioner
        Hypre::solve_euclid(this->_core, this->_solver);

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
      void* _solver;

      /// integer parameters:
      /// 0: maximum iterations, default: 1
      int _iparam[1];

      /// double parameters:
      /// 0: tolerance, default: 0.0
      double _dparam[1];

    public:
      explicit BoomerAMG(const Matrix_& matrix, const Filter_& filter) :
        BaseClass(matrix, filter),
        _solver(nullptr)
      {
        // set tolerance to 0 and maximum iterations to 1, so that BoomerAMG
        // acts as a preconditioner rather than a real iterative solver
        _iparam[0] = 1;
        _dparam[0] = 0.0;
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
        this->_solver = Hypre::create_boomeramg(this->_core, this->_iparam, this->_dparam);
      }

      virtual void done_numeric() override
      {
        if(_solver != nullptr)
        {
          Hypre::destroy_boomeramg(_solver);
          _solver = nullptr;
        }
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        // upload defect vector and format correction
        this->_upload_def(vec_def);
        this->_format_cor();

        // apply BoomerAMG preconditioner
        Hypre::solve_boomeramg(this->_core, this->_solver);

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

#endif // defined(FEAT_HAVE_HYPRE) || defined(DOXYGEN)
#endif // KERNEL_SOLVER_HYPRE_HPP
