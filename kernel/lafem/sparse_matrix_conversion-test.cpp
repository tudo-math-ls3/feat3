// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
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
SparseMatrixConversionTest <float, std::uint32_t> sparse_matrix_conversion_test_float_uint32(PreferredBackend::generic);
SparseMatrixConversionTest <double, std::uint32_t> sparse_matrix_conversion_test_double_uint32(PreferredBackend::generic);
SparseMatrixConversionTest <float, std::uint64_t> sparse_matrix_conversion_test_float_uint64(PreferredBackend::generic);
SparseMatrixConversionTest <double, std::uint64_t> sparse_matrix_conversion_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixConversionTest <float, std::uint64_t> mkl_sparse_matrix_conversion_test_float_uint64(PreferredBackend::mkl);
SparseMatrixConversionTest <double, std::uint64_t> mkl_sparse_matrix_conversion_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixConversionTest <__float128, std::uint64_t> sparse_matrix_conversion_test_float128_uint64(PreferredBackend::generic);
SparseMatrixConversionTest <__float128, std::uint32_t> sparse_matrix_conversion_test_float128_uint32(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixConversionTest <Half, std::uint32_t> sparse_matrix_conversion_test_half_uint32(PreferredBackend::generic);
SparseMatrixConversionTest <Half, std::uint64_t> sparse_matrix_conversion_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixConversionTest <float, std::uint64_t> cuda_sparse_matrix_conversion_test_float_uint64(PreferredBackend::cuda);
SparseMatrixConversionTest <double, std::uint64_t> cuda_sparse_matrix_conversion_test_double_uint64(PreferredBackend::cuda);
SparseMatrixConversionTest <float, std::uint32_t> cuda_sparse_matrix_conversion_test_float_uint32(PreferredBackend::cuda);
SparseMatrixConversionTest <double, std::uint32_t> cuda_sparse_matrix_conversion_test_double_uint32(PreferredBackend::cuda);
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
    std::cout << "seed: " << seed << "\n";

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

SparseMatrixBandedConversionTest <SparseMatrixCSR<float, std::uint32_t> > sparse_matrix_banded_csr_conversion_test_float_uint32(PreferredBackend::generic);
SparseMatrixBandedConversionTest <SparseMatrixCSR<double, std::uint32_t> > sparse_matrix_banded_csr_conversion_test_double_uint32(PreferredBackend::generic);
SparseMatrixBandedConversionTest <SparseMatrixCSR<float, std::uint64_t> > sparse_matrix_banded_csr_conversion_test_float_uint64(PreferredBackend::generic);
SparseMatrixBandedConversionTest <SparseMatrixCSR<double, std::uint64_t> > sparse_matrix_banded_csr_conversion_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SparseMatrixBandedConversionTest <SparseMatrixCSR<float, std::uint64_t> > mkl_sparse_matrix_banded_csr_conversion_test_float_uint64(PreferredBackend::mkl);
SparseMatrixBandedConversionTest <SparseMatrixCSR<double, std::uint64_t> > mkl_sparse_matrix_banded_csr_conversion_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SparseMatrixBandedConversionTest <SparseMatrixCSR<__float128, std::uint64_t> > sparse_matrix_banded_csr_conversion_test_float128_uint64(PreferredBackend::generic);
SparseMatrixBandedConversionTest <SparseMatrixCSR<__float128, std::uint32_t> > sparse_matrix_banded_csr_conversion_test_float128_uint32(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SparseMatrixBandedConversionTest <SparseMatrixCSR<Half, std::uint32_t> > sparse_matrix_banded_csr_conversion_test_half_uint32(PreferredBackend::generic);
SparseMatrixBandedConversionTest <SparseMatrixCSR<Half, std::uint64_t> > sparse_matrix_banded_csr_conversion_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SparseMatrixBandedConversionTest <SparseMatrixCSR<float, std::uint32_t> > cuda_sparse_matrix_banded_csr_conversion_test_float_uint32(PreferredBackend::cuda);
SparseMatrixBandedConversionTest <SparseMatrixCSR<double, std::uint32_t> > cuda_sparse_matrix_banded_csr_conversion_test_double_uint32(PreferredBackend::cuda);
SparseMatrixBandedConversionTest <SparseMatrixCSR<float, std::uint64_t> > cuda_sparse_matrix_banded_csr_conversion_test_float_uint64(PreferredBackend::cuda);
SparseMatrixBandedConversionTest <SparseMatrixCSR<double, std::uint64_t> > cuda_sparse_matrix_banded_csr_conversion_test_double_uint64(PreferredBackend::cuda);
#endif
