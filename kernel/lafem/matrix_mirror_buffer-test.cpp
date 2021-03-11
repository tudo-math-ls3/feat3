// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
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
  typename Mem_,
  typename DT_,
  typename IT_>
class MatrixMirrorBufferTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  MatrixMirrorBufferTest()
    : FullTaggedTest<Mem_, DT_, IT_>("MatrixMirrorBufferTest")
  {
  }

  virtual ~MatrixMirrorBufferTest()
  {
  }

  virtual void run() const override
  {
    MatrixMirrorBuffer<Mem_, DT_, IT_> zero1;

    MatrixMirrorBuffer<Mem_, DT_, IT_> a(10, 11, 30, 5);
    TEST_CHECK_EQUAL(a.rows(), 10ul);
    TEST_CHECK_EQUAL(a.columns(), 11ul);
    TEST_CHECK_EQUAL(a.used_elements(), 30ul);
    TEST_CHECK_EQUAL(a.entries_per_nonzero(), 5ul);

    MatrixMirrorBuffer<Mem_, DT_, IT_> b;
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
MatrixMirrorBufferTest<Mem::Main, float, unsigned int> cpu_matrix_mirror_buffer_test_float_uint;
MatrixMirrorBufferTest<Mem::Main, double, unsigned int> cpu_matrix_mirror_buffer_test_double_uint;
MatrixMirrorBufferTest<Mem::Main, float, unsigned long> cpu_matrix_mirror_buffer_test_float_ulong;
MatrixMirrorBufferTest<Mem::Main, double, unsigned long> cpu_matrix_mirror_buffer_test_double_ulong;
#ifdef FEAT_HAVE_QUADMATH
MatrixMirrorBufferTest<Mem::Main, __float128, unsigned int> cpu_matrix_mirror_buffer_test_float128_uint;
MatrixMirrorBufferTest<Mem::Main, __float128, unsigned long> cpu_matrix_mirror_buffer_test_float128_ulong;
#endif
#ifdef FEAT_HAVE_CUDA
MatrixMirrorBufferTest<Mem::CUDA, float, unsigned int> cuda_matrix_mirror_buffer_test_float_uint;
MatrixMirrorBufferTest<Mem::CUDA, double, unsigned int> cuda_matrix_mirror_buffer_test_double_uint;
MatrixMirrorBufferTest<Mem::CUDA, float, unsigned long> cuda_matrix_mirror_buffer_test_float_ulong;
MatrixMirrorBufferTest<Mem::CUDA, double, unsigned long> cuda_matrix_mirror_buffer_test_double_ulong;
#endif
