// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2024 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/test_matrix_factory.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

template<typename DT_, typename IT_>
class TestMatrixFactoryTest :
  public UnitTest
{
public:
  TestMatrixFactoryTest() :
    UnitTest("TestMatrixFactoryTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name())
  {
  }

  virtual void run() const override
  {
    test_generic_csr();
  }

  void test_generic_csr() const
  {
    SparseMatrixCSR<DT_, IT_> matrix;
    TestMatrixFactory factory;
    factory.seed_rng_by_timer();
    factory.create(matrix,  4u, 6u, 15u, TestMatrixFlags::generic);

    TEST_CHECK_EQUAL(matrix.rows(), 4u);
    TEST_CHECK_EQUAL(matrix.columns(), 6u);
    TEST_CHECK_EQUAL(matrix.used_elements(), 15u);
    check_structure_generic(matrix);
    check_values_generic(matrix);
  }

  void check_structure_generic(const SparseMatrixCSR<DT_, IT_>& matrix) const
  {
    const Index nrows = matrix.rows();
    const Index ncols = matrix.columns();
    const Index nnze = matrix.used_elements();
    const IT_* row_ptr = matrix.row_ptr();
    const IT_* col_idx = matrix.col_ind();

    // check row pointer
    TEST_CHECK(row_ptr[0] >= 0u);
    for(Index i(0); i < nrows; ++i)
    {
      TEST_CHECK(row_ptr[i] < row_ptr[i+1u]);
    }
    TEST_CHECK(row_ptr[nrows] <= nnze);

    // check column indices
    for(Index i(0); i < nrows; ++i)
    {
      for(IT_ j = row_ptr[i]; j < row_ptr[i+1]; ++j)
      {
        TEST_CHECK(col_idx[j] >= IT_(0));
        TEST_CHECK(col_idx[j] < ncols);
        if((j+1) < row_ptr[i+1])
        {
          // column indices must be ascending
          TEST_CHECK(col_idx[j] < col_idx[j+1]);
        }
      }
    }
  }

  void check_values_generic(const SparseMatrixCSR<DT_, IT_>& matrix) const
  {
    const Index nnze = matrix.used_elements();
    const DT_* val = matrix.val();
    for(Index i(0); i < nnze; ++i)
    {
      TEST_CHECK_IN_RANGE(val[i], DT_(-1), DT_(1));
    }
  }
}; // class TestMatrixFactoryTest

TestMatrixFactoryTest<double, std::uint64_t> test_matrix_factory_test_double_uint64;
