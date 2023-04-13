// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/lafem/matrix_mirror_buffer.hpp>
#include <kernel/util/binary_stream.hpp>

#include <list>
#include <sstream>
#include <cstdio>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for the matrix mirror buffer class.
 *
 * \test test description missing
 *
 * \tparam Mem_
 * description missing
 *
 * \tparam DT_
 * description missing
 *
 * \author Dirk Ribbrock
 */
template<
  typename DT_,
  typename IT_>
class MatrixMirrorBufferTest
  : public UnitTest
{
public:
  MatrixMirrorBufferTest(PreferredBackend backend)
    : UnitTest("MatrixMirrorBufferTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~MatrixMirrorBufferTest()
  {
  }

  virtual void run() const override
  {
    MatrixMirrorBuffer<DT_, IT_> zero1;

    MatrixMirrorBuffer<DT_, IT_> a(10, 11, 30, 5);
    TEST_CHECK_EQUAL(a.rows(), 10ul);
    TEST_CHECK_EQUAL(a.columns(), 11ul);
    TEST_CHECK_EQUAL(a.used_elements(), 30ul);
    TEST_CHECK_EQUAL(a.entries_per_nonzero(), 5ul);

    MatrixMirrorBuffer<DT_, IT_> b;
    b.convert(a);
    TEST_CHECK_EQUAL(b.val(), a.val());
    TEST_CHECK_EQUAL(b.col_ind(), a.col_ind());
    TEST_CHECK_EQUAL(b.row_ptr(), a.row_ptr());

    auto c = a.clone();
    TEST_CHECK_NOT_EQUAL(c.val(), a.val());
    TEST_CHECK_EQUAL(c.col_ind(), a.col_ind());
    TEST_CHECK_EQUAL(c.row_ptr(), a.row_ptr());
  }
};
MatrixMirrorBufferTest<float, unsigned int> cpu_matrix_mirror_buffer_test_float_uint(PreferredBackend::generic);
MatrixMirrorBufferTest<double, unsigned int> cpu_matrix_mirror_buffer_test_double_uint(PreferredBackend::generic);
MatrixMirrorBufferTest<float, unsigned long> cpu_matrix_mirror_buffer_test_float_ulong(PreferredBackend::generic);
MatrixMirrorBufferTest<double, unsigned long> cpu_matrix_mirror_buffer_test_double_ulong(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
MatrixMirrorBufferTest<float, unsigned long> mkl_matrix_mirror_buffer_test_float_ulong(PreferredBackend::mkl);
MatrixMirrorBufferTest<double, unsigned long> mkl_matrix_mirror_buffer_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
MatrixMirrorBufferTest<__float128, unsigned int> matrix_mirror_buffer_test_float128_uint(PreferredBackend::generic);
MatrixMirrorBufferTest<__float128, unsigned long> matrix_mirror_buffer_test_float128_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
MatrixMirrorBufferTest<Half, unsigned int> matrix_mirror_buffer_test_half_uint(PreferredBackend::generic);
MatrixMirrorBufferTest<Half, unsigned long> matrix_mirror_buffer_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
MatrixMirrorBufferTest<float, unsigned int> cuda_matrix_mirror_buffer_test_float_uint(PreferredBackend::cuda);
MatrixMirrorBufferTest<double, unsigned int> cuda_matrix_mirror_buffer_test_double_uint(PreferredBackend::cuda);
MatrixMirrorBufferTest<float, unsigned long> cuda_matrix_mirror_buffer_test_float_ulong(PreferredBackend::cuda);
MatrixMirrorBufferTest<double, unsigned long> cuda_matrix_mirror_buffer_test_double_ulong(PreferredBackend::cuda);
#endif
