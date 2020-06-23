// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>

#ifdef FEAT_HAVE_UMFPACK
#include <kernel/solver/umfpack.hpp>
FEAT_DISABLE_WARNINGS
#include <umfpack.h>
FEAT_RESTORE_WARNINGS

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Umfpack function wrapper helper class template
     *
     * The UMFPACK library contains multiple versions of its functions which
     * offer the same functionality for different data/index types.
     * Most notably, UMFPACK provides functions for two types on indices, which
     * are used for the row-pointer and column-index arrays:
     * - The "*_di_*" functions work with indices of type 'int', which is always
     *   a 32-bit integer type
     * - The "*_dl_*" functions work with indices of type 'SuiteSparse_long', which
     *   is a typedef that coincides to '__int64' for 64-bit Windows systems and
     *   'long' for any other platform.
     *
     * Unfortunately, this leads to a problem on 32-bit Windows systems: both versions
     * work with 32-bit ints and there is no variant for 64-bit ints on 32-bit Windows
     * systems.
     *
     * The task of this wrapper class is to choose the corresponding UMFPACK functions
     * depending on the size of the FEAT::Index type and the current platform.
     *
     * \tparam Idx_
     * The type that is used as an index type for the row-pointer and column-index arrays.
     * This type usually coincides with FEAT::Index.
     *
     * \tparam sidx_
     * The size of Idx_ in bytes.
     *
     * \tparam eq_
     * Specifies whether 'int' and 'SuiteSparse_long' have the same size.
     * If so, the specialisation for 64-bit integers is disabled, as UMFPACK does not offer it.
     *
     * \author Peter Zajac
     */
    template<typename Idx_, int sidx_ = sizeof(Idx_), bool eq_ = (sizeof(int) == sizeof(SuiteSparse_long))>
    struct UmfpackWrapper;

    // specialisation for Idx_ = 32-bit integer
    template<typename Idx_, bool eq_>
    struct UmfpackWrapper<Idx_, sizeof(int), eq_>
    {
      static void init_defaults(double* control)
      {
        ::umfpack_di_defaults(control);
      }

      static int init_symbolic(Idx_ nrows, Idx_ ncols, const Idx_* row_ptr, const Idx_* col_idx,
        const double* data, void** symb, const double* ctrl, double* info)
      {
        static_assert(sizeof(Idx_) == sizeof(int), "invalid index size");
        return int(::umfpack_di_symbolic(static_cast<int>(nrows), static_cast<int>(ncols),
          reinterpret_cast<const int*>(row_ptr), reinterpret_cast<const int*>(col_idx),
          data, symb, ctrl, info));
      }

      static int init_numeric(const Idx_* row_ptr, const Idx_* col_idx,
        const double* data, void* symb, void** nume, const double* ctrl, double* info)
      {
        static_assert(sizeof(Idx_) == sizeof(int), "invalid index size");
        return int(::umfpack_di_numeric(
          reinterpret_cast<const int*>(row_ptr), reinterpret_cast<const int*>(col_idx),
          data, symb, nume, ctrl, info));
      }

      static void free_symbolic(void** symb)
      {
        ::umfpack_di_free_symbolic(symb);
      }

      static void free_numeric(void** nume)
      {
        ::umfpack_di_free_numeric(nume);
      }

      static int solve(int sys, const Idx_* row_ptr, const Idx_* col_idx, const double* data,
        double* x, const double* b, void* nume, const double* control, double* info)
      {
        static_assert(sizeof(Idx_) == sizeof(int), "invalid index size");
        return int(::umfpack_di_solve(static_cast<int>(sys),
          reinterpret_cast<const int*>(row_ptr), reinterpret_cast<const int*>(col_idx),
          data, x, b, nume, control, info));
      }
    };

    /// specialisation for Idx_ = 64-bit integer (only non-32-bit-windows systems)
    template<typename Idx_>
    struct UmfpackWrapper<Idx_, sizeof(SuiteSparse_long), false>
    {
      static void init_defaults(double* control)
      {
        ::umfpack_dl_defaults(control);
      }

      static int init_symbolic(Idx_ nrows, Idx_ ncols, const Idx_* row_ptr, const Idx_* col_idx,
        const double* data, void** symb, const double* ctrl, double* info)
      {
        static_assert(sizeof(Idx_) == sizeof(SuiteSparse_long), "invalid index size");
        return int(::umfpack_dl_symbolic(static_cast<SuiteSparse_long>(nrows), static_cast<SuiteSparse_long>(ncols),
          reinterpret_cast<const SuiteSparse_long*>(row_ptr), reinterpret_cast<const SuiteSparse_long*>(col_idx),
          data, symb, ctrl, info));
      }

      static int init_numeric(const Idx_* row_ptr, const Idx_* col_idx,
        const double* data, void* symb, void** nume, const double* ctrl, double* info)
      {
        static_assert(sizeof(Idx_) == sizeof(SuiteSparse_long), "invalid index size");
        return int(::umfpack_dl_numeric(
          reinterpret_cast<const SuiteSparse_long*>(row_ptr), reinterpret_cast<const SuiteSparse_long*>(col_idx),
          data, symb, nume, ctrl, info));
      }

      static void free_symbolic(void** symb)
      {
        ::umfpack_dl_free_symbolic(symb);
      }

      static void free_numeric(void** nume)
      {
        ::umfpack_dl_free_numeric(nume);
      }

      static int solve(int sys, const Idx_* row_ptr, const Idx_* col_idx, const double* data,
        double* x, const double* b, void* nume, const double* control, double* info)
      {
        static_assert(sizeof(Idx_) == sizeof(SuiteSparse_long), "invalid index size");
        return int(::umfpack_dl_solve(static_cast<SuiteSparse_long>(sys),
          reinterpret_cast<const SuiteSparse_long*>(row_ptr), reinterpret_cast<const SuiteSparse_long*>(col_idx),
          data, x, b, nume, control, info));
      }
    };

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    Umfpack::Umfpack(const Umfpack::MatrixType& system_matrix) :
      _system_matrix(system_matrix),
      _umf_control(new double[UMFPACK_CONTROL]),
      _umf_symbolic(nullptr),
      _umf_numeric(nullptr),
      _sym_peak_size(0),
      _sym_mem_size(0),
      _num_mem_size(0),
      _umf_peak_size(0)
    {
      // initialise default umfpack control values
      UmfpackWrapper<Index>::init_defaults(_umf_control);
    }

    /// virtual destructor
    Umfpack::~Umfpack()
    {
      done();
      if(_umf_control != nullptr)
        delete [] _umf_control;
    }

    void Umfpack::init_symbolic()
    {
      BaseClass::init_symbolic();

      // ensure that we don't have a system matrix assigned
      if(_umf_symbolic != nullptr)
        throw SolverException("UMFPACK: already have symbolic factorisation");

      // umfpack info array
      double info[UMFPACK_INFO];

      // try to perform symbolic factorisation
      int status = UmfpackWrapper<Index>::init_symbolic(
        _system_matrix.rows(),
        _system_matrix.columns(),
        _system_matrix.row_ptr(),
        _system_matrix.col_ind(),
        nullptr,
        &_umf_symbolic,
        _umf_control,
        info
        );

      // check status code
      switch(status)
      {
      case UMFPACK_OK:
        break;
      case UMFPACK_ERROR_out_of_memory:
        throw SolverException("UMFPACK: out of memory");
      case UMFPACK_ERROR_invalid_matrix:
        throw InvalidMatrixStructureException("UMFPACK: invalid matrix structure");
      case UMFPACK_ERROR_internal_error:
        throw SolverException("UMFPACK: internal error");
      default:
        throw SolverException("UMFPACK: unknown error");
      }

      // gather statistics
      _sym_peak_size = std::size_t(info[UMFPACK_SIZE_OF_UNIT] * info[UMFPACK_SYMBOLIC_PEAK_MEMORY]);
      _sym_mem_size = std::size_t(info[UMFPACK_SIZE_OF_UNIT] * info[UMFPACK_SYMBOLIC_SIZE]);
    }

    void Umfpack::done_symbolic()
    {
      if(_umf_symbolic == nullptr)
        return;

      UmfpackWrapper<Index>::free_symbolic(&_umf_symbolic);

      _umf_symbolic = nullptr;

      BaseClass::done_symbolic();
    }

    void Umfpack::init_numeric()
    {
      BaseClass::init_numeric();

      if(_umf_symbolic == nullptr)
        throw SolverException("UFMPACK: symbolic factorisation missing");

      // umfpack info array
      double info[UMFPACK_INFO];

      // try to perform symbolic factorisation
      int status = UmfpackWrapper<Index>::init_numeric(
        _system_matrix.row_ptr(),
        _system_matrix.col_ind(),
        _system_matrix.val(),
        _umf_symbolic,
        &_umf_numeric,
        _umf_control,
        info
        );

      // check status code
      switch(status)
      {
      case UMFPACK_OK:
        break;
      case UMFPACK_ERROR_out_of_memory:
        throw SolverException("UMFPACK: out of memory");
      case UMFPACK_ERROR_invalid_matrix:
        throw InvalidMatrixStructureException("UMFPACK: invalid matrix structure");
      case UMFPACK_WARNING_singular_matrix:
        throw SingularMatrixException("UMFPACK: singular matrix");
      case UMFPACK_ERROR_different_pattern:
        throw SolverException("UMFPACK: different pattern");
      case UMFPACK_ERROR_internal_error:
        throw SolverException("UMFPACK: internal error");
      default:
        throw SolverException("UMFPACK: unknown error");
      }

      // gather statistics
      _umf_peak_size = std::size_t(info[UMFPACK_SIZE_OF_UNIT] * info[UMFPACK_PEAK_MEMORY]);
      _num_mem_size = std::size_t(info[UMFPACK_SIZE_OF_UNIT] * info[UMFPACK_NUMERIC_SIZE]);
    }

    void Umfpack::done_numeric()
    {
      if(_umf_numeric == nullptr)
        return;

      UmfpackWrapper<Index>::free_numeric(&_umf_numeric);

      _umf_numeric = nullptr;
      BaseClass::done_numeric();
    }

    Status Umfpack::apply(VectorType& x, const VectorType& b)
    {
      // umfpack info array
      double info[UMFPACK_INFO];

      // solve
      int status = UmfpackWrapper<Index>::solve(
        UMFPACK_At,
        _system_matrix.row_ptr(),
        _system_matrix.col_ind(),
        _system_matrix.val(),
        x.elements(),
        b.elements(),
        _umf_numeric,
        _umf_control,
        info);

      // check status code
      return (status == UMFPACK_OK) ? Status::success :  Status::aborted;
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    UmfpackMean::UmfpackMean(const MatrixType& system_matrix, const VectorType& weight_vector) :
      BaseClass(),
      _system_matrix(system_matrix),
      _weight_vector(weight_vector),
      _umfpack(_solver_matrix)
    {
    }

    void UmfpackMean::init_symbolic()
    {
      BaseClass::init_symbolic();

      // get the number of rows/columns
      const Index n = _weight_vector.size();
      if((n != _system_matrix.rows()) || (n != _system_matrix.columns()))
        throw InvalidMatrixStructureException("system matrix/weight vector dimension mismatch");

      // get the number of non-zeroes
      const Index nnze = _system_matrix.used_elements();

      // allocate our solver matrix
      _solver_matrix = MatrixType(n+1, n+1, nnze + 2*n);

      // get our input matrix arrays
      const Index* irow_ptr = _system_matrix.row_ptr();
      const Index* icol_idx = _system_matrix.col_ind();

      // get our output matrix arrays
      Index* orow_ptr = _solver_matrix.row_ptr();
      Index* ocol_idx = _solver_matrix.col_ind();

      // assemble the solver matrix structure
      orow_ptr[0] = Index(0);
      for(Index i(0); i < n; ++i)
      {
        Index op = orow_ptr[i];
        // copy input matrix row pattern
        for(Index ip(irow_ptr[i]); ip < irow_ptr[i+1]; ++ip, ++op)
        {
          ocol_idx[op] = icol_idx[ip];
        }
        // append a single entry in the last column
        ocol_idx[op] = n;
        // set next row pointer
        orow_ptr[i+1] = ++op;
      }

      // append an (almost) dense row (length n of n+1)
      Index op(orow_ptr[n]);
      for(Index j(0); j < n; ++j, ++op)
      {
        ocol_idx[op] = j;
      }
      // set last row pointer
      orow_ptr[n+1] = op;

      // Okay, solver matrix structure assembly complete

      // create temporary vectors
      _vec_x = _solver_matrix.create_vector_r();
      _vec_b = _solver_matrix.create_vector_r();

      // initialise umfpack
      _umfpack.init_symbolic();
    }

    void UmfpackMean::done_symbolic()
    {
      _umfpack.done_symbolic();
      _vec_b.clear();
      _vec_x.clear();
      _solver_matrix.clear();

      BaseClass::done_symbolic();
    }

    void UmfpackMean::init_numeric()
    {
      BaseClass::init_numeric();

      const Index n = _system_matrix.rows();

      // get input matrix arrays
      const Index* irow_ptr = _system_matrix.row_ptr();
      const double* idata = _system_matrix.val();

      // get input vector array
      const double* weight = _weight_vector.elements();

      // get output matrix arrays
      const Index* orow_ptr = _solver_matrix.row_ptr();
      double* odata = _solver_matrix.val();

      // loop over all input matrix rows
      for(Index i(0); i < n; ++i)
      {
        Index op(orow_ptr[i]);
        // copy input matrix row data
        for(Index ip(irow_ptr[i]); ip < irow_ptr[i+1]; ++ip, ++op)
        {
          odata[op] = idata[ip];
        }
        // copy weight vector entry as last column entry
        odata[op] = weight[i];
      }

      // copy weight vector into last row
      Index op(orow_ptr[n]);
      for(Index j(0); j < n; ++j, ++op)
      {
        odata[op] = weight[j];
      }

      // okay, solver matrix data assembly completed

      // initialise umfpack
      _umfpack.init_numeric();
    }

    void UmfpackMean::done_numeric()
    {
      _umfpack.done_numeric();

      BaseClass::done_numeric();
    }

    Status UmfpackMean::apply(VectorType& vec_sol, const VectorType& vec_rhs)
    {
      const Index n = vec_sol.size();

      // copy RHS vector
      const double* vr = vec_rhs.elements();
      double* vb = _vec_b.elements();
      for(Index i(0); i < n; ++i)
      {
        vb[i] = vr[i];
      }
      vb[n] = 0.0;

      // solve system
      Status status = _umfpack.apply(_vec_x, _vec_b);

      // copy sol vector
      const double* vx = _vec_x.elements();
      double* vs = vec_sol.elements();
      for(Index i(0); i < n; ++i)
      {
        vs[i] = vx[i];
      }

      // okay
      return status;
    }

  } // namespace Solver
} // namespace FEAT
#else
// insert dummy function to suppress linker warnings
void dummy_function() {}
#endif // FEAT_HAVE_UMFPACK
