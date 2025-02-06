// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
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
    test_symmetric_struct_csr();
    test_symmetric_csr();
    test_non_negative_csr();
    test_non_empty_diag_csr();
    test_non_zero_diag_csr();
    test_diagonal_dominant_csr();
    test_non_negative_diagonal_dominant_csr();
    test_non_negative_symmetric_csr();

    test_generic_bcsr();
    test_symmetric_struct_bcsr();
    test_symmetric_bcsr();
    test_non_negative_bcsr();
    test_non_empty_diag_bcsr();
    test_non_zero_diag_bcsr();
    test_diagonal_dominant_bcsr();
    test_non_negative_diagonal_dominant_bcsr();
    test_non_negative_symmetric_bcsr();
  }

  void test_generic_csr() const
  {
    SparseMatrixCSR<DT_, IT_> matrix;
    TestMatrixFactory factory;
    std::cout << "RNG Seed: " << factory.get_rng_seed() << std::endl;
    factory.create(matrix,  4u, 6u, 15u, TestMatrixFlags::generic);

    TEST_CHECK_EQUAL(matrix.rows(), 4u);
    TEST_CHECK_EQUAL(matrix.columns(), 6u);
    TEST_CHECK_EQUAL(matrix.used_elements(), 15u);
    check_structure_generic(matrix);
    check_values_generic(matrix);
  }

  void test_generic_bcsr() const
  {
    SparseMatrixBCSR<DT_, IT_, 2u, 3u> matrix;
    TestMatrixFactory factory;
    std::cout << "RNG Seed: " << factory.get_rng_seed() << std::endl;
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
      TEST_CHECK(row_ptr[i] <= row_ptr[i+1u]);
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

  template<int BlockHeight_, int BlockWidth_>
  void check_structure_generic(const SparseMatrixBCSR<DT_, IT_, BlockHeight_, BlockWidth_>& matrix) const
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
      TEST_CHECK(row_ptr[i] <= row_ptr[i+1u]);
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

  template<int BlockHeight_, int BlockWidth_>
  void check_values_generic(const SparseMatrixBCSR<DT_, IT_, BlockHeight_, BlockWidth_>& matrix) const
  {
    const Index nnze = matrix.used_elements();
    const auto * val = matrix.val();
    for(Index i(0); i < nnze; ++i)
      for(int j(0); j < BlockHeight_; ++j)
        for(int k(0); k< BlockWidth_; ++k)
          TEST_CHECK_IN_RANGE(val[i][j][k], DT_(-1), DT_(1));
  }

  void test_symmetric_struct_csr() const
  {
    SparseMatrixCSR<DT_, IT_> matrix;
    TestMatrixFactory factory;
    std::cout << "RNG Seed: " << factory.get_rng_seed() << std::endl;
    factory.create(matrix,  6u, 6u, 15u, TestMatrixFlags::symmetric_struct);

    TEST_CHECK_EQUAL(matrix.rows(), 6u);
    TEST_CHECK_EQUAL(matrix.columns(), 6u);
    TEST_CHECK_IN_RANGE(matrix.used_elements(), 14u, 16u);
    check_structure_symmetric(matrix);
    check_values_generic(matrix);
  }

  void test_symmetric_struct_bcsr() const
  {
    SparseMatrixBCSR<DT_, IT_, 2u, 2u> matrix;
    TestMatrixFactory factory;
    std::cout << "RNG Seed: " << factory.get_rng_seed() << std::endl;
    factory.create(matrix,  6u, 6u, 15u, TestMatrixFlags::symmetric_struct);

    TEST_CHECK_EQUAL(matrix.rows(), 6u);
    TEST_CHECK_EQUAL(matrix.columns(), 6u);
    TEST_CHECK_IN_RANGE(matrix.used_elements(), 14u, 16u);
    check_structure_symmetric(matrix);
    check_values_generic(matrix);
  }

  void check_structure_symmetric(const SparseMatrixCSR<DT_, IT_>& matrix) const
  {
    const Index nrows = matrix.rows();
    const Index ncols = matrix.columns();
    const IT_* row_ptr = matrix.row_ptr();
    const IT_* col_idx = matrix.col_ind();
    bool symmetric_entry = false;

    TEST_CHECK(nrows == ncols);

    // Check that the matrix has a symmetric structure
    for(Index i = 0; i < nrows; ++i)
    {
      for(IT_ j = row_ptr[i]; j < row_ptr[i+1]; ++j)
      {
        IT_ col_j = col_idx[j];  // Column index of the entry

        // Ensure that for every entry at (i, col_j) there is an entry at (col_j, i)
        symmetric_entry = false;

        for(IT_ k = row_ptr[col_j]; k < row_ptr[col_j + 1]; ++k)
        {
          if(col_idx[k] == i)
          {
            symmetric_entry = true;
            break;
          }
        }
        TEST_CHECK(symmetric_entry);
      }
    }
  }

  template<int BlockHeight_, int BlockWidth_>
  void check_structure_symmetric(const SparseMatrixBCSR<DT_, IT_, BlockHeight_, BlockWidth_>& matrix) const
  {
    const Index nrows = matrix.rows();
    const Index ncols = matrix.columns();
    const IT_* row_ptr = matrix.row_ptr();
    const IT_* col_idx = matrix.col_ind();
    bool symmetric_entry = false;

    TEST_CHECK(nrows == ncols);
    TEST_CHECK(BlockWidth_ == BlockHeight_);


    // Check that the matrix has a symmetric structure
    for(Index i = 0; i < nrows; ++i)
    {
      for(IT_ j = row_ptr[i]; j < row_ptr[i+1]; ++j)
      {
        IT_ col_j = col_idx[j];  // Column index of the entry

        // Ensure that for every entry at (i, col_j) there is an entry at (col_j, i)
        symmetric_entry = false;

        for(IT_ k = row_ptr[col_j]; k < row_ptr[col_j + 1]; ++k)
        {
          if(col_idx[k] == i)
          {
            symmetric_entry = true;
            break;
          }
        }
        TEST_CHECK(symmetric_entry);
      }
    }
  }

  void test_symmetric_csr() const
  {
    SparseMatrixCSR<DT_, IT_> matrix;
    TestMatrixFactory factory;
    std::cout << "RNG Seed: " << factory.get_rng_seed() << std::endl;
    factory.create(matrix, 6u, 6u, 15u, TestMatrixFlags::symmetric);

    TEST_CHECK_EQUAL(matrix.rows(), 6u);
    TEST_CHECK_EQUAL(matrix.columns(), 6u);
    TEST_CHECK_IN_RANGE(matrix.used_elements(), 14u, 16u);
    check_symmetric(matrix);
    check_values_generic(matrix);
  }

  void test_symmetric_bcsr() const
  {
    SparseMatrixBCSR<DT_, IT_, 2u, 2u> matrix;
    TestMatrixFactory factory;
    std::cout << "RNG Seed: " << factory.get_rng_seed() << std::endl;
    factory.create(matrix, 6u, 6u, 15u, TestMatrixFlags::symmetric);

    TEST_CHECK_EQUAL(matrix.rows(), 6u);
    TEST_CHECK_EQUAL(matrix.columns(), 6u);
    TEST_CHECK_IN_RANGE(matrix.used_elements(), 14u, 16u);
    check_symmetric(matrix);
    check_values_generic(matrix);
  }

  void check_symmetric(const SparseMatrixCSR<DT_, IT_>& matrix) const
  {
    const Index nrows = matrix.rows();
    const Index ncols = matrix.columns();
    const IT_* row_ptr = matrix.row_ptr();
    const IT_* col_idx = matrix.col_ind();
    const DT_* values = matrix.val();
    bool symmetric_entry = false;

    TEST_CHECK(nrows == ncols);

    // Check that the matrix is symmetric
    for(Index i = 0; i < nrows; ++i)
    {
      for(IT_ j = row_ptr[i]; j < row_ptr[i+1]; ++j)
      {
        IT_ col_j = col_idx[j];  // Column index of the entry

        // Ensure that the entry at (i, col_j) = (col_j, i)
        symmetric_entry = false;

        for(IT_ k = row_ptr[col_j]; k < row_ptr[col_j + 1]; ++k)
        {
          if(col_idx[k] == i)
          {
            // Ensure that the values are equal for symmetry
            TEST_CHECK(values[j] == values[k]);
            symmetric_entry = true;
            break;
          }
        }
        TEST_CHECK(symmetric_entry);
      }
    }
  }

  template<int BlockHeight_, int BlockWidth_>
  void check_symmetric(const SparseMatrixBCSR<DT_, IT_, BlockHeight_, BlockWidth_>& matrix) const
  {
    const DT_ eps(Math::pow(Math::eps<DT_>(), DT_(0.8)));
    const Index nrows = matrix.rows();
    const Index ncols = matrix.columns();
    const IT_* row_ptr = matrix.row_ptr();
    const IT_* col_idx = matrix.col_ind();
    const auto* values = matrix.val();
    bool symmetric_entry = false;

    TEST_CHECK(nrows == ncols);
    TEST_CHECK(BlockHeight_ == BlockWidth_);

    // Check that the matrix is symmetric
    for(Index i = 0; i < nrows; ++i)
    {
      for(IT_ j = row_ptr[i]; j < row_ptr[i+1]; ++j)
      {
        IT_ col_j = col_idx[j];  // Column index of the entry

        // Ensure that the entry at (i, col_j) = (col_j, i)

        for(IT_ k = row_ptr[col_j]; k < row_ptr[col_j + 1]; ++k)
        {
          if(col_idx[k] == i)
          {
            for(int l(0); l < BlockHeight_; ++l) {
              for(int m(0); m < BlockWidth_; ++m) {
                // Ensure that the values are equal for symmetry
                TEST_CHECK_EQUAL_WITHIN_EPS(values[k][m][l], values[j][l][m], eps);
              }
            }
            symmetric_entry = true;
            break;
          }
        }
        TEST_CHECK(symmetric_entry);
      }
    }
  }

  void test_non_negative_csr() const
  {
    SparseMatrixCSR<DT_, IT_> matrix;
    TestMatrixFactory factory;
    std::cout << "RNG Seed: " << factory.get_rng_seed() << std::endl;
    factory.create(matrix,  4u, 6u, 15u, TestMatrixFlags::non_negative);

    TEST_CHECK_EQUAL(matrix.rows(), 4u);
    TEST_CHECK_EQUAL(matrix.columns(), 6u);
    TEST_CHECK_EQUAL(matrix.used_elements(), 15u);
    check_structure_generic(matrix);
    check_values_non_negative(matrix);
  }

  void test_non_negative_bcsr() const
  {
    SparseMatrixBCSR<DT_, IT_, 2u, 3u> matrix;
    TestMatrixFactory factory;
    std::cout << "RNG Seed: " << factory.get_rng_seed() << std::endl;
    factory.create(matrix,  4u, 6u, 15u, TestMatrixFlags::non_negative);

    TEST_CHECK_EQUAL(matrix.rows(), 4u);
    TEST_CHECK_EQUAL(matrix.columns(), 6u);
    TEST_CHECK_EQUAL(matrix.used_elements(), 15u);
    check_structure_generic(matrix);
    check_values_non_negative(matrix);
  }

  void check_values_non_negative(const SparseMatrixCSR<DT_, IT_>& matrix) const
  {
    const Index nnze = matrix.used_elements();
    const DT_* val = matrix.val();
    for(Index i(0); i < nnze; ++i)
    {
      TEST_CHECK(val[i] > DT_(0));
    }
  }

  template<int BlockHeight_, int BlockWidth_>
  void check_values_non_negative(const SparseMatrixBCSR<DT_, IT_, BlockHeight_, BlockWidth_>& matrix) const
  {
    const Index nnze = matrix.used_elements();
    const auto* val = matrix.val();
    for(Index i(0); i < nnze; ++i)
      for(int j(0); j < BlockHeight_; ++j)
        for(int k(0); k < BlockWidth_; ++k)
          TEST_CHECK(val[i][j][k] > DT_(0));
  }

  void test_non_empty_diag_csr() const
  {
    SparseMatrixCSR<DT_, IT_> matrix;
    TestMatrixFactory factory;
    std::cout << "RNG Seed: " << factory.get_rng_seed() << std::endl;
    factory.create(matrix,  6u, 6u, 15u, TestMatrixFlags::non_empty_diag);

    TEST_CHECK_EQUAL(matrix.rows(), 6u);
    TEST_CHECK_EQUAL(matrix.columns(), 6u);
    TEST_CHECK_EQUAL(matrix.used_elements(), 15u);
    check_structure_non_empty_diag(matrix);
    check_values_generic(matrix);
  }

  void test_non_empty_diag_bcsr() const
  {
    SparseMatrixBCSR<DT_, IT_, 2u, 2u> matrix;
    TestMatrixFactory factory;
    std::cout << "RNG Seed: " << factory.get_rng_seed() << std::endl;
    factory.create(matrix,  6u, 6u, 15u, TestMatrixFlags::non_empty_diag);

    TEST_CHECK_EQUAL(matrix.rows(), 6u);
    TEST_CHECK_EQUAL(matrix.columns(), 6u);
    TEST_CHECK_EQUAL(matrix.used_elements(), 15u);
    check_structure_non_empty_diag(matrix);
    check_values_generic(matrix);
  }

  void check_structure_non_empty_diag(const SparseMatrixCSR<DT_, IT_>& matrix) const
  {
    const Index nrows = matrix.rows();
    const Index ncols = matrix.columns();
    const IT_* row_ptr = matrix.row_ptr();
    const IT_* col_idx = matrix.col_ind();
    bool diagonal_entry = false;

    TEST_CHECK(nrows == ncols);

    for (Index i = 0; i < nrows; ++i)
    {
      diagonal_entry = false;

      for (IT_ j = row_ptr[i]; j < row_ptr[i+1]; ++j)
      {
        if (col_idx[j] == i) // Check if the column index matches the row index
        {
          diagonal_entry = true;
          break;
        }
      }
      TEST_CHECK(diagonal_entry);
    }
  }

  template<int BlockHeight_, int BlockWidth_>
  void check_structure_non_empty_diag(const SparseMatrixBCSR<DT_, IT_, BlockHeight_, BlockWidth_>& matrix) const
  {
    const Index nrows = matrix.rows();
    const Index ncols = matrix.columns();
    const IT_* row_ptr = matrix.row_ptr();
    const IT_* col_idx = matrix.col_ind();
    bool diagonal_entry = false;

    TEST_CHECK(nrows == ncols);
    TEST_CHECK(BlockHeight_ == BlockWidth_);

    for (Index i = 0; i < nrows; ++i)
    {
      diagonal_entry = false;

      for (IT_ j = row_ptr[i]; j < row_ptr[i+1]; ++j)
      {
        if (col_idx[j] == i) // Check if the column index matches the row index
        {
          diagonal_entry = true;
          break;
        }
      }
      TEST_CHECK(diagonal_entry);
    }
  }

  void test_non_zero_diag_csr() const
  {
    SparseMatrixCSR<DT_, IT_> matrix;
    TestMatrixFactory factory;
    std::cout << "RNG Seed: " << factory.get_rng_seed() << std::endl;
    factory.create(matrix,  6u, 6u, 15u, TestMatrixFlags::non_zero_diag);

    TEST_CHECK_EQUAL(matrix.rows(), 6u);
    TEST_CHECK_EQUAL(matrix.columns(), 6u);
    TEST_CHECK_EQUAL(matrix.used_elements(), 15u);
    check_structure_non_empty_diag(matrix);
    check_values_non_zero_diag(matrix);
  }

  void test_non_zero_diag_bcsr() const
  {
    SparseMatrixBCSR<DT_, IT_, 2u, 2u> matrix;
    TestMatrixFactory factory;
    std::cout << "RNG Seed: " << factory.get_rng_seed() << std::endl;
    factory.create(matrix,  6u, 6u, 15u, TestMatrixFlags::non_zero_diag);

    TEST_CHECK_EQUAL(matrix.rows(), 6u);
    TEST_CHECK_EQUAL(matrix.columns(), 6u);
    TEST_CHECK_EQUAL(matrix.used_elements(), 15u);
    check_structure_non_empty_diag(matrix);
    check_values_non_zero_diag(matrix);
  }

  void check_values_non_zero_diag(const SparseMatrixCSR<DT_, IT_>& matrix) const
  {
    const Index nrows = matrix.rows();
    const Index ncols = matrix.columns();
    const IT_* row_ptr = matrix.row_ptr();
    const IT_* col_idx = matrix.col_ind();
    const DT_* values = matrix.val();

    TEST_CHECK(nrows == ncols);

    for (Index i = 0; i < nrows; ++i)
    {
      for (IT_ j = row_ptr[i]; j < row_ptr[i+1]; ++j)
      {
        if (col_idx[j] == i) // Check if the column index matches the row index
        {
          TEST_CHECK_NOT_EQUAL(values[j], DT_(0));
          break;
        }
      }
    }
  }

  template<int BlockHeight_, int BlockWidth_>
  void check_values_non_zero_diag(const SparseMatrixBCSR<DT_, IT_, BlockWidth_, BlockHeight_>& matrix) const
  {
    const Index nrows = matrix.rows();
    const Index ncols = matrix.columns();
    const IT_* row_ptr = matrix.row_ptr();
    const IT_* col_idx = matrix.col_ind();
    const auto* values = matrix.val();

    TEST_CHECK(nrows == ncols);
    TEST_CHECK(BlockWidth_ == BlockHeight_);

    for (Index i = 0; i < nrows; ++i)
    {
      for (IT_ j = row_ptr[i]; j < row_ptr[i+1]; ++j)
      {
        if (col_idx[j] == i) // Check if the column index matches the row index
        {
          for(int k(0); k < BlockWidth_; ++k)
            TEST_CHECK_NOT_EQUAL(values[j][k][k], DT_(0));
          break;
        }
      }
    }
  }

  void test_diagonal_dominant_csr() const
  {
    SparseMatrixCSR<DT_, IT_> matrix;
    TestMatrixFactory factory;
    std::cout << "RNG Seed: " << factory.get_rng_seed() << std::endl;
    factory.create(matrix,  6u, 6u, 15u, TestMatrixFlags::diagonal_dominant);

    TEST_CHECK_EQUAL(matrix.rows(), 6u);
    TEST_CHECK_EQUAL(matrix.columns(), 6u);
    TEST_CHECK_EQUAL(matrix.used_elements(), 15u);
    check_structure_non_empty_diag(matrix);
    check_values_diagonal_dominant(matrix);
  }

  void test_diagonal_dominant_bcsr() const
  {
    SparseMatrixBCSR<DT_, IT_, 2u, 2u> matrix;
    TestMatrixFactory factory;
    std::cout << "RNG Seed: " << factory.get_rng_seed() << std::endl;
    factory.create(matrix,  6u, 6u, 15u, TestMatrixFlags::diagonal_dominant);

    TEST_CHECK_EQUAL(matrix.rows(), 6u);
    TEST_CHECK_EQUAL(matrix.columns(), 6u);
    TEST_CHECK_EQUAL(matrix.used_elements(), 15u);
    check_structure_non_empty_diag(matrix);
    check_values_diagonal_dominant(matrix);
  }

  void check_values_diagonal_dominant(const SparseMatrixCSR<DT_, IT_>& matrix) const
  {
    const Index nrows = matrix.rows();
    const Index ncols = matrix.columns();
    const IT_* row_ptr = matrix.row_ptr();
    const IT_* col_idx = matrix.col_ind();
    const DT_* values = matrix.val();

    TEST_CHECK(nrows == ncols);

    for(Index i = 0; i < nrows; ++i)
    {
      DT_ off_diagonal = 0;
      for(IT_ j = row_ptr[i]; j < row_ptr[i + 1]; ++j)
      {
        if(col_idx[j] != i) // Skip the diagonal element
        {
          off_diagonal += std::abs(values[j]);
        }
      }

      for(IT_ j = row_ptr[i]; j < row_ptr[i + 1]; ++j)
      {
        if(col_idx[j] == i) // Check if this is the diagonal element
        {
          TEST_CHECK(values[j] > off_diagonal);
          break;
        }
      }
    }
  }

  template<int BlockHeight_, int BlockWidth_>
  void check_values_diagonal_dominant(const SparseMatrixBCSR<DT_, IT_,BlockWidth_, BlockHeight_>& matrix) const
  {
    const Index nrows = matrix.rows();
    const Index ncols = matrix.columns();
    const IT_* row_ptr = matrix.row_ptr();
    const IT_* col_idx = matrix.col_ind();
    const auto* values = matrix.val();

    TEST_CHECK(nrows == ncols);
    TEST_CHECK(BlockWidth_ == BlockHeight_);

    for(Index i = 0; i < nrows; ++i)
    {
      DT_ off_diagonal = 0;
      for(IT_ j = row_ptr[i]; j < row_ptr[i + 1]; ++j)
      {
        if(col_idx[j] != i) // Skip the diagonal block element
          for(int k(0); k < BlockWidth_; ++k)
            for(int l(0); l < BlockHeight_; ++l)
              off_diagonal += std::abs(values[j][k][l]);
        else
          for(int k(0); k < BlockWidth_; ++k)
            for(int l(0); l < BlockHeight_; ++l)
              if(k != l) // Skip the diagonal element
                off_diagonal += std::abs(values[j][k][l]);
      }

      for(IT_ j = row_ptr[i]; j < row_ptr[i + 1]; ++j)
      {
        if(col_idx[j] == i) // Check if this is the diagonal block element
        {
          for(int k(0); k < BlockWidth_; ++k)
            TEST_CHECK(values[j][k][k] > off_diagonal);
          break;
        }
      }
    }
  }

  void test_non_negative_diagonal_dominant_csr() const
  {
    SparseMatrixCSR<DT_, IT_> matrix;
    TestMatrixFactory factory;
    std::cout << "RNG Seed: " << factory.get_rng_seed() << std::endl;
    factory.create(matrix,  6u, 6u, 15u, TestMatrixFlags::non_negative | TestMatrixFlags::diagonal_dominant);

    TEST_CHECK_EQUAL(matrix.rows(), 6u);
    TEST_CHECK_EQUAL(matrix.columns(), 6u);
    TEST_CHECK_EQUAL(matrix.used_elements(), 15u);
    check_values_non_negative(matrix);
    check_structure_non_empty_diag(matrix);
    check_values_diagonal_dominant(matrix);
  }

  void test_non_negative_diagonal_dominant_bcsr() const
  {
    SparseMatrixBCSR<DT_, IT_, 2u, 2u> matrix;
    TestMatrixFactory factory;
    std::cout << "RNG Seed: " << factory.get_rng_seed() << std::endl;
    factory.create(matrix,  6u, 6u, 15u, TestMatrixFlags::non_negative | TestMatrixFlags::diagonal_dominant);

    TEST_CHECK_EQUAL(matrix.rows(), 6u);
    TEST_CHECK_EQUAL(matrix.columns(), 6u);
    TEST_CHECK_EQUAL(matrix.used_elements(), 15u);
    check_values_non_negative(matrix);
    check_structure_non_empty_diag(matrix);
    check_values_diagonal_dominant(matrix);
  }

  void test_non_negative_symmetric_csr() const
  {
    SparseMatrixCSR<DT_, IT_> matrix;
    TestMatrixFactory factory;
    std::cout << "RNG Seed: " << factory.get_rng_seed() << std::endl;
    factory.create(matrix,  6u, 6u, 15u, TestMatrixFlags::non_negative | TestMatrixFlags::symmetric);

    TEST_CHECK_EQUAL(matrix.rows(), 6u);
    TEST_CHECK_EQUAL(matrix.columns(), 6u);
    TEST_CHECK_IN_RANGE(matrix.used_elements(), 14u, 16u);
    check_values_non_negative(matrix);
    check_symmetric(matrix);
  }

  void test_non_negative_symmetric_bcsr() const
  {
    SparseMatrixBCSR<DT_, IT_, 2u, 2u> matrix;
    TestMatrixFactory factory;
    std::cout << "RNG Seed: " << factory.get_rng_seed() << std::endl;
    factory.create(matrix,  6u, 6u, 15u, TestMatrixFlags::non_negative | TestMatrixFlags::symmetric);

    TEST_CHECK_EQUAL(matrix.rows(), 6u);
    TEST_CHECK_EQUAL(matrix.columns(), 6u);
    TEST_CHECK_IN_RANGE(matrix.used_elements(), 14u, 16u);
    check_values_non_negative(matrix);
    check_symmetric(matrix);
  }
}; // class TestMatrixFactoryTest

TestMatrixFactoryTest<double, std::uint64_t> test_matrix_factory_test_double_uint64;
