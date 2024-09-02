// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/power_vector.hpp>
#include <kernel/lafem/tuple_vector.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/mean_filter.hpp>
#include <kernel/lafem/power_filter.hpp>
#include <kernel/lafem/tuple_filter.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

template<
  typename DataType_,
  typename IndexType_>
class MetaFilterTest
  : public UnitTest
{
public:
  typedef DataType_ DataType;
  typedef IndexType_ IndexType;

  typedef DenseVector<DataType, IndexType> ScalarVector;
  typedef PowerVector<ScalarVector, 2> PowerVector2;
  typedef TupleVector<PowerVector2, ScalarVector> MetaVector;

  typedef UnitFilter<DataType, IndexType> ScalarFilter1;
  typedef MeanFilter<DataType, IndexType> ScalarFilter2;
  typedef PowerFilter<ScalarFilter1, 2> PowerFilter2;
  typedef TupleFilter<PowerFilter2, ScalarFilter2> MetaFilter;

   MetaFilterTest(PreferredBackend backend)
    : UnitTest("MetaFilterTest", Type::Traits<DataType>::name(), Type::Traits<IndexType>::name(), backend)
  {
  }

  virtual ~MetaFilterTest()
  {
  }

  static MetaFilter gen_filter(Index m)
  {
    // create a unit-filter
    ScalarVector fv(2);
    fv(0, DataType(1));
    fv(1, DataType(5));
    DenseVector<IndexType, IndexType> idx(2);
    idx(0, 0);
    XASSERT(m > 0);
    idx(1, IndexType(m-1));
    ScalarFilter1 unit_filter(m, fv, idx);

    // create vectors for mean-filter
    ScalarVector mfv(m, DataType(1)), mfw(m, DataType(0));
    DataType* fw(mfw.elements());
    for(Index i(0); i < m; ++i)
      fw[i] = DataType(i+1);

    // create a mean-filter
    ScalarFilter2 mean_filter(std::move(mfv), std::move(mfw), DataType(0), DataType(((m+1)*(m+2))/2));

    // create a power-filer
    PowerFilter2 power_filter;
    power_filter.template at<0>().clone(unit_filter);
    power_filter.template at<1>().clone(unit_filter);

    // return the tuple-filter
    return MetaFilter(std::move(power_filter), std::move(mean_filter));
  }

  static MetaVector gen_vector(Index m)
  {
    PowerVector2 vec;
    vec.template at<0>() = ScalarVector(m, DataType(2));
    vec.template at<1>() = ScalarVector(m, DataType(3));

    return MetaVector(std::move(vec), ScalarVector(m, DataType(1)));
  }

  static MetaVector gen_vector_sol(Index m)
  {
    ScalarVector vx(m, DataType(2));
    ScalarVector vy(m, DataType(3));
    ScalarVector vz(m, DataType(2) / DataType(7));

    DataType* fx(vx.elements());
    DataType* fy(vy.elements());
    fx[0] = fy[0] = DataType(1);
    XASSERT(m > 0);
    fx[m-1] = fy[m-1] = DataType(5);

    // create a power-vector
    PowerVector2 vec;
    vec.template at<0>().convert(vx);
    vec.template at<1>().convert(vy);

    ScalarVector tvz;
    tvz.convert(vz);
    return MetaVector(std::move(vec), std::move(tvz));
  }

  static MetaVector gen_vector_def(Index m)
  {
    ScalarVector vx(m, DataType(2));
    ScalarVector vy(m, DataType(3));
    ScalarVector vz(m, DataType(0));

    DataType* fx(vx.elements());
    DataType* fy(vy.elements());
    DataType* fz(vz.elements());
    XASSERT(m > 0);
    fx[0] = fy[0] = fx[m-1] = fy[m-1] = DataType(0);
    for(Index i(0); i < m; ++i)
    {
      fz[i] = DataType_(32 - 10*int(i)) / DataType_(42);
    }

    // create a power-vector
    PowerVector2 vec;
    vec.template at<0>().convert(vx);
    vec.template at<1>().convert(vy);

    ScalarVector tvz;
    tvz.convert(vz);
    return MetaVector(std::move(vec), std::move(tvz));
  }

  virtual void run() const override
  {
    const DataType tol(Math::pow(Math::eps<DataType>(), DataType(0.8)));

    const Index m(5);

    // create a power-filter
    MetaFilter filter(gen_filter(m));

    // generate two input vector
    MetaVector vec_sol(gen_vector(m));
    MetaVector vec_def(vec_sol.clone());

    // appy sol filter
    filter.filter_sol(vec_sol);
    filter.filter_def(vec_def);

    // generate ref vectors
    const MetaVector ref_sol(gen_vector_sol(m));
    const MetaVector ref_def(gen_vector_def(m));

    // subtract reference
    vec_sol.axpy(ref_sol, -DataType(1));
    vec_def.axpy(ref_def, -DataType(1));

    // check norm
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_sol.norm2(), DataType(0), tol);
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_def.norm2(), DataType(0), tol);
  }
};

MetaFilterTest <float, std::uint64_t> meta_filter_test_float_uint64(PreferredBackend::generic);
MetaFilterTest <double, std::uint64_t> meta_filter_test_double_uint64(PreferredBackend::generic);
MetaFilterTest <float, std::uint32_t> meta_filter_test_float_uint32(PreferredBackend::generic);
MetaFilterTest <double, std::uint32_t> meta_filter_test_double_uint32(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
MetaFilterTest <float, std::uint64_t> mkl_meta_filter_test_float_uint64(PreferredBackend::mkl);
MetaFilterTest <double, std::uint64_t> mkl_meta_filter_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
MetaFilterTest <__float128, std::uint64_t> meta_filter_test_float128_uint64(PreferredBackend::generic);
MetaFilterTest <__float128, std::uint32_t> meta_filter_test_float128_uint32(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
MetaFilterTest <Half, std::uint32_t> meta_filter_test_half_uint32(PreferredBackend::generic);
MetaFilterTest <Half, std::uint64_t> meta_filter_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
MetaFilterTest <float, std::uint64_t> cuda_meta_filter_test_float_uint64(PreferredBackend::cuda);
MetaFilterTest <double, std::uint64_t> cuda_meta_filter_test_double_uint64(PreferredBackend::cuda);
MetaFilterTest <float, std::uint32_t> cuda_meta_filter_test_float_uint32(PreferredBackend::cuda);
MetaFilterTest <double, std::uint32_t> cuda_meta_filter_test_double_uint32(PreferredBackend::cuda);
#endif
