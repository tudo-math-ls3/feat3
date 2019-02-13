// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
//#include <kernel/archs.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_banded.hpp>
#include <kernel/util/random.hpp>
#include <kernel/adjacency/permutation.hpp>
#include <kernel/lafem/sparse_matrix_factory.hpp>

#include <cstdio>
#include <sstream>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for sparse matrix conversions.
 *
 * \test test description missing
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
  typename DT_,
  typename IT_>
class SparseMatrixConversionTest
  : public UnitTest
{

public:
  SparseMatrixConversionTest(PreferredBackend backend)
    : UnitTest("sparse_matrix_conversion_test", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixConversionTest()
  {
  }

  virtual void run() const override
  {
    SparseMatrixFactory<DT_, IT_> a(4, 4);
    for (Index row(0) ; row < a.rows() ; ++row)
    {
      for (Index col(0) ; col < a.columns() ; ++col)
      {
        if(row == col)
          a.add(row, col, DT_(2));
        else if((row == col+1) || (row+1 == col))
          a.add(row, col, DT_(-1));
      }
    }
    SparseMatrixCSR<DT_, IT_> csr_m1(a.make_csr());
    SparseMatrixBanded<DT_, IT_> banded_m1;
    banded_m1.convert(csr_m1);
    SparseMatrixCSR<DT_, IT_> csr_m2;
    csr_m2.convert(banded_m1);
    TEST_CHECK_EQUAL(csr_m1.used_elements(), banded_m1.used_elements());
    TEST_CHECK_EQUAL(csr_m2.used_elements(), banded_m1.used_elements());
    TEST_CHECK_EQUAL(csr_m2.used_elements(), csr_m1.used_elements());
    TEST_CHECK_EQUAL(csr_m2, csr_m1);
    SparseMatrixCSR<DT_, IT_> csr_m3(banded_m1);
    TEST_CHECK_EQUAL(csr_m3, csr_m1);
    csr_m1.format();
    csr_m2.format();
    csr_m3.convert_reverse(csr_m2);
    csr_m2.convert_reverse(csr_m1);
    TEST_CHECK_EQUAL(csr_m1, csr_m3);
  }
};
SparseMatrixConversionTest<float, unsigned int> sparse_matrix_conversion_test_float_uint(PreferredBackend::generic);
SparseMatrixConversionTest<double, unsigned int> sparse_matrix_conversion_test_double_uint(PreferredBackend::generic);
SparseMatrixConversionTest<float, unsigned long> sparse_matrix_conversion_test_float_ulong(PreferredBackend::generic);
SparseMatrixConversionTest<double, unsigned long> sparse_matrix_conversion_test_double_ulong(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixConversionTest<float, unsigned long> mkl_sparse_matrix_conversion_test_float_ulong(PreferredBackend::mkl);
SparseMatrixConversionTest<double, unsigned long> mkl_sparse_matrix_conversion_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixConversionTest<__float128, unsigned long> sparse_matrix_conversion_test_float128_ulong(PreferredBackend::generic);
SparseMatrixConversionTest<__float128, unsigned int> sparse_matrix_conversion_test_float128_uint(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixConversionTest<Half, unsigned int> sparse_matrix_conversion_test_half_uint(PreferredBackend::generic);
SparseMatrixConversionTest<Half, unsigned long> sparse_matrix_conversion_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixConversionTest<float, unsigned long> cuda_sparse_matrix_conversion_test_float_ulong(PreferredBackend::cuda);
SparseMatrixConversionTest<double, unsigned long> cuda_sparse_matrix_conversion_test_double_ulong(PreferredBackend::cuda);
SparseMatrixConversionTest<float, unsigned int> cuda_sparse_matrix_conversion_test_float_uint(PreferredBackend::cuda);
SparseMatrixConversionTest<double, unsigned int> cuda_sparse_matrix_conversion_test_double_uint(PreferredBackend::cuda);
#endif

/**
 * \brief Test class for the conversion from SparseMatrixBanded to another matrix-format
 *
 * \test test description missing
 *
 * \tparam MT_
 * description missing
 *
 * \author Christoph Lohmann
 */
template<
  typename MT_>
class SparseMatrixBandedConversionTest
  : public UnitTest
{
public:
  SparseMatrixBandedConversionTest(PreferredBackend backend)
    : UnitTest("sparse_matrix_banded_conversion_test: " + MT_::name(), Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual ~SparseMatrixBandedConversionTest()
  {
  }

  typedef typename MT_::DataType DT_;
  typedef typename MT_::IndexType IT_;

  virtual void run() const override
  {
    Random::SeedType seed(Random::SeedType(time(nullptr)));
    Random random(seed);
    std::cout << "seed: " << seed << std::endl;

    // create random matrix
    const Index tsize(20);
    const Index rows(tsize + random(Index(0), Index(20)));
    const Index columns(tsize + random(Index(0), Index(20)));

    const Index num_of_offsets(5 + random(Index(0), Index(10)));

    DenseVector<IT_, IT_> vec_offsets(num_of_offsets);
    DenseVector<DT_, IT_> vec_val(num_of_offsets * rows);

    // create random vector of offsets
    FEAT::Adjacency::Permutation permutation(rows + columns - 1, random);
    for (Index i(0); i < num_of_offsets; ++i)
    {
      vec_offsets(i, IT_(permutation.get_perm_pos()[i]));
    }
    std::sort(vec_offsets.elements(), vec_offsets.elements() + num_of_offsets);


    // fill data-array
    for (Index i(0); i < vec_val.size(); ++i)
    {
      vec_val(i, random(DT_(0), DT_(10)));
    }

    // create test-matrix
    SparseMatrixBanded<DT_, IT_> sys_banded(rows, columns, vec_val, vec_offsets);

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

SparseMatrixBandedConversionTest<SparseMatrixCSR<float, unsigned int> > sparse_matrix_banded_csr_conversion_test_float_uint(PreferredBackend::generic);
SparseMatrixBandedConversionTest<SparseMatrixCSR<double, unsigned int> > sparse_matrix_banded_csr_conversion_test_double_uint(PreferredBackend::generic);
SparseMatrixBandedConversionTest<SparseMatrixCSR<float, unsigned long> > sparse_matrix_banded_csr_conversion_test_float_ulong(PreferredBackend::generic);
SparseMatrixBandedConversionTest<SparseMatrixCSR<double, unsigned long> > sparse_matrix_banded_csr_conversion_test_double_ulong(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixBandedConversionTest<SparseMatrixCSR<float, unsigned long> > mkl_sparse_matrix_banded_csr_conversion_test_float_ulong(PreferredBackend::mkl);
SparseMatrixBandedConversionTest<SparseMatrixCSR<double, unsigned long> > mkl_sparse_matrix_banded_csr_conversion_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixBandedConversionTest<SparseMatrixCSR<__float128, unsigned long> > sparse_matrix_banded_csr_conversion_test_float128_ulong(PreferredBackend::generic);
SparseMatrixBandedConversionTest<SparseMatrixCSR<__float128, unsigned int> > sparse_matrix_banded_csr_conversion_test_float128_uint(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixBandedConversionTest<SparseMatrixCSR<Half, unsigned int> > sparse_matrix_banded_csr_conversion_test_half_uint(PreferredBackend::generic);
SparseMatrixBandedConversionTest<SparseMatrixCSR<Half, unsigned long> > sparse_matrix_banded_csr_conversion_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixBandedConversionTest<SparseMatrixCSR<float, unsigned int> > cuda_sparse_matrix_banded_csr_conversion_test_float_uint(PreferredBackend::cuda);
SparseMatrixBandedConversionTest<SparseMatrixCSR<double, unsigned int> > cuda_sparse_matrix_banded_csr_conversion_test_double_uint(PreferredBackend::cuda);
SparseMatrixBandedConversionTest<SparseMatrixCSR<float, unsigned long> > cuda_sparse_matrix_banded_csr_conversion_test_float_ulong(PreferredBackend::cuda);
SparseMatrixBandedConversionTest<SparseMatrixCSR<double, unsigned long> > cuda_sparse_matrix_banded_csr_conversion_test_double_ulong(PreferredBackend::cuda);
#endif
