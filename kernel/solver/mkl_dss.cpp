// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>

#ifdef FEAT_HAVE_MKL
#include <kernel/solver/mkl_dss.hpp>
FEAT_DISABLE_WARNINGS
#include <mkl_dss.h>
FEAT_RESTORE_WARNINGS

namespace FEAT
{
  namespace Solver
  {
    MKLDSS::MKLDSS(const MatrixType& system_matrix) :
      BaseClass(),
      _system_matrix(system_matrix),
      _mkl_dss_handle(nullptr)
    {
    }

    MKLDSS::~MKLDSS()
    {
    }

    void MKLDSS::init_symbolic()
    {
      BaseClass::init_symbolic();
      XASSERT(_mkl_dss_handle == nullptr);

      MKL_INT ret = 0;
      MKL_INT opt = 0;
      MKL_INT neq = _system_matrix.rows();
      MKL_INT nze = _system_matrix.used_elements();

      // create DSS handle
      opt = MKL_DSS_TERM_LVL_ERROR + MKL_DSS_ZERO_BASED_INDEXING;
#ifdef DEBUG
      opt += MKL_DSS_MSG_LVL_WARNING;
#endif
      ret = ::dss_create(_mkl_dss_handle, opt);

      switch(ret)
      {
      case MKL_DSS_SUCCESS:
        break;
      case MKL_DSS_OUT_OF_MEMORY:
        throw SolverException("MKL DSS: out of memory");
      default:
        throw SolverException("MKL DSS: unknown error");
      }

      // get matrix index arrays and convert to MKL_INT if necessary
      std::vector<MKL_INT> row_ptr_v, col_idx_v;
      const Index* a_row_ptr = _system_matrix.row_ptr();
      const Index* a_col_idx = _system_matrix.col_ind();
      const MKL_INT* row_ptr = nullptr;
      const MKL_INT* col_idx = nullptr;
      if constexpr(sizeof(Index) == sizeof(MKL_INT))
      {
        row_ptr = reinterpret_cast<const MKL_INT*>(a_row_ptr);
        col_idx = reinterpret_cast<const MKL_INT*>(a_col_idx);
      }
      else
      {
        row_ptr_v.resize(neq+1u);
        col_idx_v.resize(nze);
        for(Index i(0); i <= Index(neq); ++i)
          row_ptr_v[i] = MKL_INT(a_row_ptr[i]);
        for(Index i(0); i < Index(nze); ++i)
          col_idx_v[i] = MKL_INT(a_col_idx[i]);
        row_ptr = row_ptr_v.data();
        col_idx = col_idx_v.data();
      }

      // set matrix structure
      opt = MKL_DSS_NON_SYMMETRIC;
      ret = ::dss_define_structure(_mkl_dss_handle, opt, row_ptr, neq, neq, col_idx, nze);

      switch(ret)
      {
      case MKL_DSS_SUCCESS:
        break;
      case MKL_DSS_OUT_OF_MEMORY:
        throw SolverException("MKL DSS: out of memory");
      case MKL_DSS_STRUCTURE_ERR:
        throw InvalidMatrixStructureException("MKL DSS: invalid matrix structure");
      default:
        throw SolverException("MKL DSS: unknown error");
      }

      // reorder matrix
      opt = MKL_DSS_AUTO_ORDER;
      ret = ::dss_reorder(_mkl_dss_handle, opt, nullptr);

      switch(ret)
      {
      case MKL_DSS_SUCCESS:
        break;
      case MKL_DSS_OUT_OF_MEMORY:
        throw SolverException("MKL DSS: out of memory");
      default:
        throw SolverException("MKL DSS: unknown error");
      }
    }

    void MKLDSS::done_symbolic()
    {
      if(_mkl_dss_handle)
      {
        MKL_INT opt = MKL_DSS_TERM_LVL_ERROR;
#ifdef DEBUG
        opt += MKL_DSS_MSG_LVL_WARNING;
#endif
        ::dss_delete(_mkl_dss_handle, opt);
        _mkl_dss_handle = nullptr;
      }
      BaseClass::done_symbolic();
    }

    void MKLDSS::init_numeric()
    {
      BaseClass::init_numeric();
      XASSERT(_mkl_dss_handle != nullptr);

      MKL_INT ret = 0;
      MKL_INT opt = MKL_DSS_INDEFINITE;

      ret = ::dss_factor_real(_mkl_dss_handle, opt, _system_matrix.val());

      switch(ret)
      {
      case MKL_DSS_SUCCESS:
        break;
      case MKL_DSS_OUT_OF_MEMORY:
        throw SolverException("MKL DSS: out of memory");
      case MKL_DSS_ZERO_PIVOT:
        throw SingularMatrixException("MKL DSS: zero pivot encountered");
      default:
        throw SolverException("MKL DSS: unknown error");
      }
    }

    void MKLDSS::done_numeric()
    {
      // nothing to do here
      BaseClass::done_numeric();
    }

    Status MKLDSS::apply(VectorType& vec_sol, const VectorType& vec_rhs)
    {
      MKL_INT opt = 0;
      MKL_INT nrhs = 1;
      MKL_INT ret = ::dss_solve_real(_mkl_dss_handle, opt, vec_rhs.elements(), nrhs, vec_sol.elements());
      return (ret == MKL_DSS_SUCCESS ? Status::success : Status::aborted);
    }
  } // namespace Solver
} // namespace FEAT
#else
// insert dummy function to suppress linker warnings
void dummy_function_mkldss() {}
#endif // FEAT_HAVE_MKL
