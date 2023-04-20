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
MatrixMirrorBufferTest <float, std::uint32_t> cpu_matrix_mirror_buffer_test_float_uint32(PreferredBackend::generic);
MatrixMirrorBufferTest <double, std::uint32_t> cpu_matrix_mirror_buffer_test_double_uint32(PreferredBackend::generic);
MatrixMirrorBufferTest <float, std::uint64_t> cpu_matrix_mirror_buffer_test_float_uint64(PreferredBackend::generic);
MatrixMirrorBufferTest <double, std::uint64_t> cpu_matrix_mirror_buffer_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
MatrixMirrorBufferTest <float, std::uint64_t> mkl_matrix_mirror_buffer_test_float_uint64(PreferredBackend::mkl);
MatrixMirrorBufferTest <double, std::uint64_t> mkl_matrix_mirror_buffer_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
MatrixMirrorBufferTest <__float128, std::uint32_t> matrix_mirror_buffer_test_float128_uint32(PreferredBackend::generic);
MatrixMirrorBufferTest <__float128, std::uint64_t> matrix_mirror_buffer_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
MatrixMirrorBufferTest <Half, std::uint32_t> matrix_mirror_buffer_test_half_uint32(PreferredBackend::generic);
MatrixMirrorBufferTest <Half, std::uint64_t> matrix_mirror_buffer_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
MatrixMirrorBufferTest <float, std::uint32_t> cuda_matrix_mirror_buffer_test_float_uint32(PreferredBackend::cuda);
MatrixMirrorBufferTest <double, std::uint32_t> cuda_matrix_mirror_buffer_test_double_uint32(PreferredBackend::cuda);
MatrixMirrorBufferTest <float, std::uint64_t> cuda_matrix_mirror_buffer_test_float_uint64(PreferredBackend::cuda);
MatrixMirrorBufferTest <double, std::uint64_t> cuda_matrix_mirror_buffer_test_double_uint64(PreferredBackend::cuda);
#endif
