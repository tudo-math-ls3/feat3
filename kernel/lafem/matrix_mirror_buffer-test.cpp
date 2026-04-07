// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
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
  explicit MatrixMirrorBufferTest(PreferredBackend backend)
    : UnitTest("MatrixMirrorBufferTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    MatrixMirrorBuffer<DT_, IT_> zero1;

    MatrixMirrorBuffer<DT_, IT_> a(10, 11, 30, 5);
    TEST_CHECK_EQUAL(a.num_rows(), 10ul);
    TEST_CHECK_EQUAL(a.num_cols(), 11ul);
    TEST_CHECK_EQUAL(a.num_nzes(), 30ul);
    TEST_CHECK_EQUAL(a.entries_per_nonzero(), 5ul);

    MatrixMirrorBuffer<DT_, IT_> b;
    b.convert(a);
    TEST_CHECK_EQUAL(b.val_arbiter(), a.val_arbiter());
    TEST_CHECK_EQUAL(b.col_idx_arbiter(), a.col_idx_arbiter());
    TEST_CHECK_EQUAL(b.row_ptr_arbiter(), a.row_ptr_arbiter());

    auto c = a.clone();
    TEST_CHECK_NOT_EQUAL(c.val_arbiter(), a.val_arbiter());
    TEST_CHECK_EQUAL(c.col_idx_arbiter(), a.col_idx_arbiter());
    TEST_CHECK_EQUAL(c.row_ptr_arbiter(), a.row_ptr_arbiter());
  }
};
SPAWN_UNIT_TEST_2T_P(MatrixMirrorBufferTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MatrixMirrorBufferTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MatrixMirrorBufferTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MatrixMirrorBufferTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(MatrixMirrorBufferTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(MatrixMirrorBufferTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(MatrixMirrorBufferTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MatrixMirrorBufferTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(MatrixMirrorBufferTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MatrixMirrorBufferTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(MatrixMirrorBufferTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MatrixMirrorBufferTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MatrixMirrorBufferTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MatrixMirrorBufferTest, double, std::uint64_t, PreferredBackend::cuda);
#endif
