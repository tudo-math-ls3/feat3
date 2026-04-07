// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/matvecmult_bcsr_dense.hpp>
#include <kernel/util/exception.hpp>

FEAT_DISABLE_WARNINGS
#include <mkl.h>
#include <mkl_spblas.h>
FEAT_RESTORE_WARNINGS

namespace FEAT::LAFEM::Arch
{
  // the code below casts Index* to MKL_INT*
  static_assert(sizeof(Index) == sizeof(MKL_INT), "size mismatch between MKL_INT and FEAT::Index");

  void MatVecMultBCSRDense::exec_mkl_impl(float * r, const float alpha, const float * const x, const float * const y, const float * const val, const Index * const col_idx, const Index * const row_ptr, const Index rows, const Index cols, const Index, const int block_size, const bool transposed)
  {
    MKL_INT mrows = (MKL_INT)rows;
    MKL_INT mcolumns = (MKL_INT)cols;
    MKL_INT mblocksize = (MKL_INT)block_size;
    MKL_INT mcopysize = (transposed ? mcolumns : mrows) * mblocksize;

    if(y == nullptr)
      Memory::memset_main(r, 0, sizeof(float) * std::size_t(mcopysize));
    else if(r != y)
    {
      MKL_INT one = 1;
      scopy(&mcopysize, y, &one, r, &one);
    }

    float beta = 1.0f;

#ifndef FEAT_USE_MKL_LEGACY_SPMV

    sparse_operation_t opt = (transposed ? SPARSE_OPERATION_TRANSPOSE : SPARSE_OPERATION_NON_TRANSPOSE);

    sparse_matrix_t A;
    FEAT_DISABLE_WARNINGS
    mkl_sparse_s_create_bsr(&A, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, mrows, mcolumns, mblocksize, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr + 1), (MKL_INT*)col_idx, (float*) val);
    FEAT_RESTORE_WARNINGS
    matrix_descr md;
    md.type = SPARSE_MATRIX_TYPE_GENERAL;
    mkl_sparse_s_mv(opt, alpha, A, md, x, beta, r);
    mkl_sparse_destroy(A);

#else

    char trans = (transposed ? 'T' : 'N');
    char matdescra[6];
    matdescra[0] = 'G';
    matdescra[1] = 0; //ingored by mkl, but valgrind complains otherwise
    matdescra[2] = 0; //ingored by mkl, but valgrind complains otherwise
    matdescra[3] = 'C';

    FEAT_DISABLE_WARNINGS
    mkl_sbsrmv(&trans, &mrows, &mcolumns, &mblocksize, (const float*)&alpha, matdescra, (const float*) val, (const MKL_INT*)col_idx, (const MKL_INT*)row_ptr, (const MKL_INT*)(row_ptr) + 1, (const float*)x, (const float*)&beta, r);
    FEAT_RESTORE_WARNINGS

#endif
  }

  void MatVecMultBCSRDense::exec_mkl_impl(double * r, const double alpha, const double * const x, const double * const y, const double * const val, const Index * const col_idx, const Index * const row_ptr, const Index rows, const Index cols, const Index, const int block_size, const bool transposed)
  {
    MKL_INT mrows = (MKL_INT)rows;
    MKL_INT mcolumns = (MKL_INT)cols;
    MKL_INT mblocksize = (MKL_INT)block_size;
    MKL_INT mcopysize = (transposed ? mcolumns : mrows) * mblocksize;

    if(y == nullptr)
      Memory::memset_main(r, 0, sizeof(double) * std::size_t(mcopysize));
    else if(r != y)
    {
      MKL_INT one = 1;
      dcopy(&mcopysize, y, &one, r, &one);
    }

    double beta = 1.0;

#ifndef FEAT_USE_MKL_LEGACY_SPMV

    sparse_operation_t opt = (transposed ? SPARSE_OPERATION_TRANSPOSE : SPARSE_OPERATION_NON_TRANSPOSE);

    sparse_matrix_t A;
    FEAT_DISABLE_WARNINGS
    mkl_sparse_d_create_bsr(&A, SPARSE_INDEX_BASE_ZERO, SPARSE_LAYOUT_ROW_MAJOR, mrows, mcolumns, mblocksize, (MKL_INT*)row_ptr, (MKL_INT*)(row_ptr + 1), (MKL_INT*)col_idx, (double*) val);
    FEAT_RESTORE_WARNINGS
    matrix_descr md;
    md.type = SPARSE_MATRIX_TYPE_GENERAL;
    mkl_sparse_d_mv(opt, alpha, A, md, x, beta, r);
    mkl_sparse_destroy(A);

#else

    char trans = (transposed ? 'T' : 'N');
    char matdescra[6];
    matdescra[0] = 'G';
    matdescra[1] = 0; //ingored by mkl, but valgrind complains otherwise
    matdescra[2] = 0; //ingored by mkl, but valgrind complains otherwise
    matdescra[3] = 'C';

    FEAT_DISABLE_WARNINGS
    mkl_dbsrmv(&trans, &mrows, &mcolumns, &mblocksize, (const double*)&alpha, matdescra, (const double*) val, (const MKL_INT*)col_ind, (const MKL_INT*)row_ptr, (const MKL_INT*)(row_ptr) + 1, (const double*)x, (const double*)&beta, r);
    FEAT_RESTORE_WARNINGS

#endif
  }
} // namespace FEAT::LAFEM::Arch
