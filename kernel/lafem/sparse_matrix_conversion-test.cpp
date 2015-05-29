#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_banded.hpp>
#include <kernel/util/random.hpp>
#include <kernel/adjacency/permutation.hpp>

#include <cstdio>
#include <sstream>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::TestSystem;

/**
 * \brief Test class for sparse matrix conversions.
 *
 * \test test description missing
 *
 * \tparam Mem_
 * description missing
 *
 * \tparam DT_
 * description missing
 *
 * \tparam IT_
 * description missing
 *
 * \author Dirk Ribbrock
 */
template<
  typename Mem_,
  typename DT_,
  typename IT_>
class SparseMatrixConversionTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{

public:

  SparseMatrixConversionTest()
    : FullTaggedTest<Mem_, DT_, IT_>("sparse_matrix_conversion_test")
  {
  }

  virtual void run() const
  {

    SparseMatrixCOO<Mem_, DT_, IT_> a(121, 121);
    for (Index row(0) ; row < a.rows() ; ++row)
    {
      for (Index col(0) ; col < a.columns() ; ++col)
      {
        if(row == col)
          a(row, col, DT_(2));
        else if((row == col+1) || (row+1 == col))
          a(row, col, DT_(-1));
      }
    }

    {
      SparseMatrixCOO<Mem_, DT_, IT_> coo_m1;
      coo_m1.convert(a);
      SparseMatrixELL<Mem_, DT_, IT_> ell_m1(coo_m1);
      SparseMatrixCOO<Mem_, DT_, IT_> coo_m2(ell_m1);
      TEST_CHECK_EQUAL(coo_m2, coo_m1);
      SparseMatrixCSR<Mem_, DT_, IT_> csr_m1(coo_m1);
      SparseMatrixCOO<Mem_, DT_, IT_> coo_m3(csr_m1);
      TEST_CHECK_EQUAL(coo_m3, coo_m1);
    }

    {
      SparseMatrixELL<Mem_, DT_, IT_> ell_m1(a);
      SparseMatrixCOO<Mem_, DT_, IT_> coo_m1(ell_m1);
      SparseMatrixELL<Mem_, DT_, IT_> ell_m2(coo_m1);
      TEST_CHECK_EQUAL(ell_m2, ell_m1);
      SparseMatrixCSR<Mem_, DT_, IT_> csr_m1(ell_m1);
      SparseMatrixELL<Mem_, DT_, IT_> ell_m3(csr_m1);
      TEST_CHECK_EQUAL(ell_m3, ell_m1);
    }
    {
      SparseMatrixCSR<Mem_, DT_, IT_> csr_m1(a);
      SparseMatrixCOO<Mem_, DT_, IT_> coo_m1(csr_m1);
      SparseMatrixCSR<Mem_, DT_, IT_> csr_m2(coo_m1);
      TEST_CHECK_EQUAL(csr_m2, csr_m1);
      SparseMatrixELL<Mem_, DT_, IT_> ell_m1(csr_m1);
      SparseMatrixCSR<Mem_, DT_, IT_> csr_m3(ell_m1);
      TEST_CHECK_EQUAL(csr_m3, csr_m1);
    }

  }
};
SparseMatrixConversionTest<Mem::Main, float, unsigned int> sparse_matrix_conversion_test_float_uint;
SparseMatrixConversionTest<Mem::Main, double, unsigned int> sparse_matrix_conversion_test_double_uint;
SparseMatrixConversionTest<Mem::Main, float, unsigned long> sparse_matrix_conversion_test_float_ulong;
SparseMatrixConversionTest<Mem::Main, double, unsigned long> sparse_matrix_conversion_test_double_ulong;


/**
 * \brief Test class for sparse matrix conversions.
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
  typename DT_>
class SparseMatrixCudaConversionTest
  : public TaggedTest<Mem_, DT_>
{

public:

  SparseMatrixCudaConversionTest()
    : TaggedTest<Mem_, DT_>("sparse_matrix_cuda_conversion_test")
  {
  }

  virtual void run() const
  {

    SparseMatrixCOO<Mem::Main, DT_> a(121, 121);
    for (Index row(0) ; row < a.rows() ; ++row)
    {
      for (Index col(0) ; col < a.columns() ; ++col)
      {
        if(row == col)
          a(row, col, DT_(2));
        else if((row == col+1) || (row+1 == col))
          a(row, col, DT_(-1));
      }
    }

    {
      SparseMatrixCOO<Mem::CUDA, DT_> coo_m1;
      coo_m1.convert(a);
      SparseMatrixELL<Mem::CUDA, DT_> ell_m1(coo_m1);
      SparseMatrixCOO<Mem::CUDA, DT_> coo_m2(ell_m1);
      TEST_CHECK_EQUAL(coo_m2, coo_m1);
      SparseMatrixCSR<Mem::CUDA, DT_> csr_m1(coo_m1);
      SparseMatrixCOO<Mem::CUDA, DT_> coo_m3(csr_m1);
      TEST_CHECK_EQUAL(coo_m2, coo_m1);
    }

    {
      SparseMatrixELL<Mem::CUDA, DT_> ell_m1(a);
      SparseMatrixCOO<Mem::CUDA, DT_> coo_m1(ell_m1);
      SparseMatrixELL<Mem::CUDA, DT_> ell_m2(coo_m1);
      TEST_CHECK_EQUAL(ell_m2, ell_m1);
      SparseMatrixCSR<Mem::CUDA, DT_> csr_m1(ell_m1);
      SparseMatrixELL<Mem::CUDA, DT_> ell_m3(csr_m1);
      TEST_CHECK_EQUAL(ell_m3, ell_m1);
    }

    {
      SparseMatrixCSR<Mem::CUDA, DT_> csr_m1(a);
      SparseMatrixCOO<Mem::CUDA, DT_> coo_m1(csr_m1);
      SparseMatrixCSR<Mem::CUDA, DT_> csr_m2(coo_m1);
      TEST_CHECK_EQUAL(csr_m2, csr_m1);
      SparseMatrixELL<Mem::CUDA, DT_> ell_m1(csr_m1);
      SparseMatrixCSR<Mem::CUDA, DT_> csr_m3(ell_m1);
      TEST_CHECK_EQUAL(csr_m3, csr_m1);
    }

    {
      SparseMatrixCOO<Mem::CUDA, DT_> coo_m1;
      coo_m1.convert(a);
      SparseMatrixELL<Mem::Main, DT_> ell_m1(coo_m1);
      SparseMatrixCOO<Mem::CUDA, DT_> coo_m2(ell_m1);
      TEST_CHECK_EQUAL(coo_m2, coo_m1);
      SparseMatrixCSR<Mem::Main, DT_> csr_m1(coo_m1);
      SparseMatrixCOO<Mem::CUDA, DT_> coo_m3(csr_m1);
      TEST_CHECK_EQUAL(coo_m2, coo_m1);
    }

    {
      SparseMatrixELL<Mem::CUDA, DT_> ell_m1(a);
      SparseMatrixCOO<Mem::Main, DT_> coo_m1(ell_m1);
      SparseMatrixELL<Mem::CUDA, DT_> ell_m2(coo_m1);
      TEST_CHECK_EQUAL(ell_m2, ell_m1);
      SparseMatrixCSR<Mem::Main, DT_> csr_m1(ell_m1);
      SparseMatrixELL<Mem::CUDA, DT_> ell_m3(csr_m1);
      TEST_CHECK_EQUAL(ell_m3, ell_m1);
    }

    {
      SparseMatrixCSR<Mem::CUDA, DT_> csr_m1(a);
      SparseMatrixCOO<Mem::Main, DT_> coo_m1(csr_m1);
      SparseMatrixCSR<Mem::CUDA, DT_> csr_m2(coo_m1);
      TEST_CHECK_EQUAL(csr_m2, csr_m1);
      SparseMatrixELL<Mem::Main, DT_> ell_m1(csr_m1);
      SparseMatrixCSR<Mem::CUDA, DT_> csr_m3(ell_m1);
      TEST_CHECK_EQUAL(csr_m3, csr_m1);
    }

  }
};
#ifdef FEAST_BACKENDS_CUDA
SparseMatrixCudaConversionTest<Mem::CUDA, float> sparse_matrix_cuda_conversion_test_float;
SparseMatrixCudaConversionTest<Mem::CUDA, double> sparse_matrix_cuda_conversion_test_double;
#endif

/**
 * \brief Test class for the convertion from SparseMatrixBanded to another matrix-format
 *
 * \test test description missing
 *
 * \tparam MT_
 * description missing
 *
 * \author Christoph Lohmann
 */
template<typename MT_>
class SparseMatrixBandedConversionTest
  : public FullTaggedTest<typename MT_::MemType, typename MT_::DataType, typename MT_::IndexType>
{
public:
  SparseMatrixBandedConversionTest()
    : FullTaggedTest<typename MT_::MemType, typename MT_::DataType, typename MT_::IndexType>("sparse_matrix_banded_conversion_test: " + MT_::name())
  {
  }

  typedef typename MT_::MemType Mem_;
  typedef typename MT_::DataType DT_;
  typedef typename MT_::IndexType IT_;

  virtual void run() const
  {
    Random::SeedType seed(Random::SeedType(time(nullptr)));
    Random random(seed);
    std::cout << "seed: " << seed << std::endl;

    // create random matrix
    const Index tsize(20);
    const Index rows(tsize + random(Index(0), Index(20)));
    const Index columns(tsize + random(Index(0), Index(20)));

    const Index num_of_offsets(5 + random(Index(0), Index(10)));

    DenseVector<Mem::Main, IT_, IT_> tvec_offsets(num_of_offsets);
    DenseVector<Mem_, DT_, IT_> vec_val(num_of_offsets * rows);

    // create random vector of offsets
    FEAST::Adjacency::Permutation permutation(rows + columns - 1, random);
    for (Index i(0); i < num_of_offsets; ++i)
    {
      tvec_offsets(i, IT_(permutation.get_perm_pos()[i]));
    }
    std::sort(tvec_offsets.elements(), tvec_offsets.elements() + num_of_offsets);

    DenseVector<Mem_, IT_, IT_> vec_offsets;
    vec_offsets.convert(tvec_offsets);

    // fill data-array
    for (Index i(0); i < vec_val.size(); ++i)
    {
      vec_val(i, random(DT_(0), DT_(10)));
    }

    // create test-matrix
    SparseMatrixBanded<Mem_, DT_, IT_> sys_banded(rows, columns, vec_val, vec_offsets);

    MT_ sys_other(sys_banded);

    for (Index i(0); i < sys_banded.rows(); ++i)
    {
      for (Index j(0); j < sys_banded.columns(); ++j)
      {
        TEST_CHECK_EQUAL(sys_banded(i, j), sys_other(i, j));
      }
    }
  }
};

SparseMatrixBandedConversionTest<SparseMatrixCSR<Mem::Main, float, unsigned int> > sparse_matrix_banded_csr_conversion_test_float_uint;
SparseMatrixBandedConversionTest<SparseMatrixCSR<Mem::Main, double, unsigned int> > sparse_matrix_banded_csr_conversion_test_double_uint;
SparseMatrixBandedConversionTest<SparseMatrixELL<Mem::Main, float, unsigned int> > sparse_matrix_banded_ell_conversion_test_float_uint;
SparseMatrixBandedConversionTest<SparseMatrixELL<Mem::Main, double, unsigned int> > sparse_matrix_banded_ell_conversion_test_double_uint;
SparseMatrixBandedConversionTest<SparseMatrixCOO<Mem::Main, float, unsigned int> > sparse_matrix_banded_coo_conversion_test_float_uint;
SparseMatrixBandedConversionTest<SparseMatrixCOO<Mem::Main, double, unsigned int> > sparse_matrix_banded_coo_conversion_test_double_uint;

SparseMatrixBandedConversionTest<SparseMatrixCSR<Mem::Main, float, unsigned long> > sparse_matrix_banded_csr_conversion_test_float_ulong;
SparseMatrixBandedConversionTest<SparseMatrixCSR<Mem::Main, double, unsigned long> > sparse_matrix_banded_csr_conversion_test_double_ulong;
SparseMatrixBandedConversionTest<SparseMatrixELL<Mem::Main, float, unsigned long> > sparse_matrix_banded_ell_conversion_test_float_ulong;
SparseMatrixBandedConversionTest<SparseMatrixELL<Mem::Main, double, unsigned long> > sparse_matrix_banded_ell_conversion_test_double_ulong;
SparseMatrixBandedConversionTest<SparseMatrixCOO<Mem::Main, float, unsigned long> > sparse_matrix_banded_coo_conversion_test_float_ulong;
SparseMatrixBandedConversionTest<SparseMatrixCOO<Mem::Main, double, unsigned long> > sparse_matrix_banded_coo_conversion_test_double_ulong;

#ifdef FEAST_BACKENDS_CUDA
SparseMatrixBandedConversionTest<SparseMatrixCSR<Mem::CUDA, float, unsigned int> > cuda_sparse_matrix_banded_csr_conversion_test_float_uint;
SparseMatrixBandedConversionTest<SparseMatrixCSR<Mem::CUDA, double, unsigned int> > cuda_sparse_matrix_banded_csr_conversion_test_double_uint;
SparseMatrixBandedConversionTest<SparseMatrixELL<Mem::CUDA, float, unsigned int> > cuda_sparse_matrix_banded_ell_conversion_test_float_uint;
SparseMatrixBandedConversionTest<SparseMatrixELL<Mem::CUDA, double, unsigned int> > cuda_sparse_matrix_banded_ell_conversion_test_double_uint;
SparseMatrixBandedConversionTest<SparseMatrixCOO<Mem::CUDA, float, unsigned int> > cuda_sparse_matrix_banded_coo_conversion_test_float_uint;
SparseMatrixBandedConversionTest<SparseMatrixCOO<Mem::CUDA, double, unsigned int> > cuda_sparse_matrix_banded_coo_conversion_test_double_uint;

SparseMatrixBandedConversionTest<SparseMatrixCSR<Mem::CUDA, float, unsigned long> > cuda_sparse_matrix_banded_csr_conversion_test_float_ulong;
SparseMatrixBandedConversionTest<SparseMatrixCSR<Mem::CUDA, double, unsigned long> > cuda_sparse_matrix_banded_csr_conversion_test_double_ulong;
SparseMatrixBandedConversionTest<SparseMatrixELL<Mem::CUDA, float, unsigned long> > cuda_sparse_matrix_banded_ell_conversion_test_float_ulong;
SparseMatrixBandedConversionTest<SparseMatrixELL<Mem::CUDA, double, unsigned long> > cuda_sparse_matrix_banded_ell_conversion_test_double_ulong;
SparseMatrixBandedConversionTest<SparseMatrixCOO<Mem::CUDA, float, unsigned long> > cuda_sparse_matrix_banded_coo_conversion_test_float_ulong;
SparseMatrixBandedConversionTest<SparseMatrixCOO<Mem::CUDA, double, unsigned long> > cuda_sparse_matrix_banded_coo_conversion_test_double_ulong;
#endif
