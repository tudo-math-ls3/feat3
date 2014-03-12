#include <kernel/lafem/umfpack.hpp>
#ifdef FEAST_HAVE_UMFPACK
#include <umfpack.h>

namespace FEAST
{
  namespace LAFEM
  {
    template<typename Idx_, int sidx = sizeof(Idx_)>
    struct UmfpackWrapper;

    template<typename Idx_>
    struct UmfpackWrapper<Idx_, sizeof(int)>
    {
      static void init_defaults(double* control)
      {
        ::umfpack_di_defaults(control);
      }

      static int init_symbolic(Idx_ nrows, Idx_ ncols, const Idx_* row_ptr, const Idx_* col_idx,
        const double* data, void** symb, const double* ctrl, double* info)
      {
        static_assert(sizeof(Idx_) == sizeof(int), "invalid index size");
        return ::umfpack_di_symbolic(static_cast<int>(nrows), static_cast<int>(ncols),
          reinterpret_cast<const int*>(row_ptr), reinterpret_cast<const int*>(col_idx),
          data, symb, ctrl, info);
      }

      static int init_numeric(const Idx_* row_ptr, const Idx_* col_idx,
        const double* data, void* symb, void** nume, const double* ctrl, double* info)
      {
        static_assert(sizeof(Idx_) == sizeof(int), "invalid index size");
        return ::umfpack_di_numeric(
          reinterpret_cast<const int*>(row_ptr), reinterpret_cast<const int*>(col_idx),
          data, symb, nume, ctrl, info);
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
        return ::umfpack_di_solve(static_cast<int>(sys),
          reinterpret_cast<const int*>(row_ptr), reinterpret_cast<const int*>(col_idx),
          data, x, b, nume, control, info);
      }
    };

    template<typename Idx_>
    struct UmfpackWrapper<Idx_, sizeof(SuiteSparse_long)>
    {
      static void init_defaults(double* control)
      {
        ::umfpack_dl_defaults(control);
      }

      static int init_symbolic(Idx_ nrows, Idx_ ncols, const Idx_* row_ptr, const Idx_* col_idx,
        const double* data, void** symb, const double* ctrl, double* info)
      {
        static_assert(sizeof(Idx_) == sizeof(int), "invalid index size");
        return ::umfpack_dl_symbolic(static_cast<SuiteSparse_long>(nrows), static_cast<SuiteSparse_long>(ncols),
          reinterpret_cast<const SuiteSparse_long*>(row_ptr), reinterpret_cast<const SuiteSparse_long*>(col_idx),
          data, symb, ctrl, info);
      }

      static int init_numeric(const Idx_* row_ptr, const Idx_* col_idx,
        const double* data, void* symb, void** nume, const double* ctrl, double* info)
      {
        static_assert(sizeof(Idx_) == sizeof(int), "invalid index size");
        return ::umfpack_dl_numeric(
          reinterpret_cast<const SuiteSparse_long*>(row_ptr), reinterpret_cast<const SuiteSparse_long*>(col_idx),
          data, symb, nume, ctrl, info);
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
        return ::umfpack_dl_solve(static_cast<SuiteSparse_long>(sys),
          reinterpret_cast<const SuiteSparse_long*>(row_ptr), reinterpret_cast<const SuiteSparse_long*>(col_idx),
          data, x, b, nume, control, info);
      }
    };

    Umfpack::Umfpack() :
      _system_matrix(nullptr),
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

    Umfpack::Umfpack(const MatrixType* system_matrix) :
      Umfpack()
    {
      init(system_matrix);
    }

    /// virtual destructor
    Umfpack::~Umfpack()
    {
      free();
      if(_umf_control != nullptr)
        delete [] _umf_control;
    }

    void Umfpack::init_symbolic(const MatrixType* system_matrix, bool ignore_data)
    {
      // ensure that we don't have a system matrix assigned
      if((_system_matrix != nullptr) || (_umf_symbolic != nullptr))
        throw InternalError("already have symbolic factorisation");

      // store system matrix
      ASSERT_(system_matrix != nullptr);
      _system_matrix = system_matrix;

      // umfpack info array
      double info[UMFPACK_INFO];

      // try to perform symbolic factorisation
      int status = UmfpackWrapper<Index>::init_symbolic(
        _system_matrix->rows(),
        _system_matrix->columns(),
        _system_matrix->row_ptr(),
        _system_matrix->col_ind(),
        ignore_data ? nullptr : _system_matrix->val(),
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
        throw InternalError("UMFPACK: out of memory");
      case UMFPACK_ERROR_invalid_matrix:
        throw InternalError("UMFPACK: invalid matrix");
      case UMFPACK_ERROR_internal_error:
        throw InternalError("UMFPACK: internal error");
      default:
        throw InternalError("UMFPACK: unknown error");
      }

      // gather statistics
      _sym_peak_size = std::size_t(info[UMFPACK_SIZE_OF_UNIT] * info[UMFPACK_SYMBOLIC_PEAK_MEMORY]);
      _sym_mem_size = std::size_t(info[UMFPACK_SIZE_OF_UNIT] * info[UMFPACK_SYMBOLIC_SIZE]);
    }

    void Umfpack::init_numeric()
    {
      if((_system_matrix == nullptr) || (_umf_symbolic == nullptr))
        throw InternalError("symbolic factorisation missing");

      // umfpack info array
      double info[UMFPACK_INFO];

      // try to perform symbolic factorisation
      int status = UmfpackWrapper<Index>::init_numeric(
        _system_matrix->row_ptr(),
        _system_matrix->col_ind(),
        _system_matrix->val(),
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
        throw InternalError("UMFPACK: out of memory");
      case UMFPACK_ERROR_invalid_matrix:
        throw InternalError("UMFPACK: invalid matrix");
      case UMFPACK_ERROR_different_pattern:
        throw InternalError("UMFPACK: different pattern");
      case UMFPACK_WARNING_singular_matrix:
        throw InternalError("UMFPACK: singular matrix");
      case UMFPACK_ERROR_internal_error:
        throw InternalError("UMFPACK: internal error");
      default:
        throw InternalError("UMFPACK: unknown error");
      }

      // gather statistics
      _umf_peak_size = std::size_t(info[UMFPACK_SIZE_OF_UNIT] * info[UMFPACK_PEAK_MEMORY]);
      _num_mem_size = std::size_t(info[UMFPACK_SIZE_OF_UNIT] * info[UMFPACK_NUMERIC_SIZE]);
    }

    void Umfpack::init(const MatrixType* system_matrix)
    {
      init_symbolic(system_matrix);
      init_numeric();
    }

    void Umfpack::free_numeric()
    {
      if(_umf_numeric == nullptr)
        return;

      UmfpackWrapper<Index>::free_numeric(&_umf_numeric);

      _umf_numeric = nullptr;
    }

    void Umfpack::free_symbolic()
    {
      if(_umf_symbolic == nullptr)
        return;

      UmfpackWrapper<Index>::free_symbolic(&_umf_symbolic);

      _umf_symbolic = nullptr;
      _system_matrix = nullptr;
    }

    void Umfpack::free()
    {
      free_numeric();
      free_symbolic();
    }

    void Umfpack::solve(VectorType& x, const VectorType& b)
    {
      // umfpack info array
      double info[UMFPACK_INFO];

      // solve
      int status = UmfpackWrapper<Index>::solve(
        UMFPACK_At,
        _system_matrix->row_ptr(),
        _system_matrix->col_ind(),
        _system_matrix->val(),
        x.elements(),
        b.elements(),
        _umf_numeric,
        _umf_control,
        info);

      // check status code
      switch(status)
      {
      case UMFPACK_OK:
        break;
      default:
        throw InternalError("UMFPACK: error");
      }
    }
  } // namespace LAFEM
} // namespace FEAST
#endif // FEAST_HAVE_UMFPACK
