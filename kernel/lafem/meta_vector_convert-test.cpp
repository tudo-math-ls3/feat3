// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>
#include <kernel/util/random.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/power_vector.hpp>
#include <kernel/lafem/tuple_vector.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief MetaVectorConvertTest
 *
 * \author Peter Zajac
 */
template<
  typename DT_,
  typename IT_>
class MetaVectorConvertTest
  : public UnitTest
{
public:
  typedef DenseVector<DT_, IT_> DV;
  typedef DenseVectorBlocked<DT_, IT_, 1> DVB1;
  typedef DenseVectorBlocked<DT_, IT_, 2> DVB2;
  typedef DenseVectorBlocked<DT_, IT_, 3> DVB3;
  typedef PowerVector<DVB3, 1> PV1;
  typedef PowerVector<DVB2, 2> PV2;
  typedef PowerVector<DVB1, 3> PV3;
  typedef TupleVector<PV1, PV2, PV3, DV> TV;

  explicit MetaVectorConvertTest(PreferredBackend backend)
    : UnitTest("MetaVectorConvertTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = TestSystem::tol<DT_>();

    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";

    static const Index block_sizes[] = {3,2,2,1,1,1,1};

    // we have 7 total sub-vectors, so choose 8 random vector sizes between 1 and 11
    Index sizes[7], raw_sizes[7], raw_offsets[8] = {0};
    for(int i = 0; i < 7; ++i)
    {
      sizes[i] = rng(Index(1), Index(11));
      raw_sizes[i] = block_sizes[i] * sizes[i];
      raw_offsets[i+1] = raw_offsets[i] + raw_sizes[i];
    }
    const Index total_raw_size = raw_offsets[7];

    // create a tuple vector and get references to all its scalar/blocked sub-vectors
    TV tv;
    DVB3& tv_dvb3   = tv.template at<0>().template at<0>();
    DVB2& tv_dvb2_0 = tv.template at<1>().template at<0>();
    DVB2& tv_dvb2_1 = tv.template at<1>().template at<1>();
    DVB1& tv_dvb1_0 = tv.template at<2>().template at<0>();
    DVB1& tv_dvb1_1 = tv.template at<2>().template at<1>();
    DVB1& tv_dvb1_2 = tv.template at<2>().template at<2>();
    DV&   tv_dv     = tv.template at<3>();

    // allocate all sub-vectors to their corresponding size
    tv_dvb3   = DVB3(sizes[0]);
    tv_dvb2_0 = DVB2(sizes[1]);
    tv_dvb2_1 = DVB2(sizes[2]);
    tv_dvb1_0 = DVB1(sizes[3]);
    tv_dvb1_1 = DVB1(sizes[4]);
    tv_dvb1_2 = DVB1(sizes[5]);
    tv_dv     = DV  (sizes[6]);

    TEST_CHECK_EQUAL(tv.size_raw(), total_raw_size);

    // allocate a dense vector of total raw size and fill it with random values
    DV ref_vector(rng, total_raw_size, DT_(1), DT_(2));

    // derive sub-vectors of our reference vector
    DVB3 ref_dvb3  (sizes[0], ref_vector.elements_arbiter().attach(raw_offsets[0]*sizeof(DT_), raw_sizes[0]*sizeof(DT_)));
    DVB2 ref_dvb2_0(sizes[1], ref_vector.elements_arbiter().attach(raw_offsets[1]*sizeof(DT_), raw_sizes[1]*sizeof(DT_)));
    DVB2 ref_dvb2_1(sizes[2], ref_vector.elements_arbiter().attach(raw_offsets[2]*sizeof(DT_), raw_sizes[2]*sizeof(DT_)));
    DVB1 ref_dvb1_0(sizes[3], ref_vector.elements_arbiter().attach(raw_offsets[3]*sizeof(DT_), raw_sizes[3]*sizeof(DT_)));
    DVB1 ref_dvb1_1(sizes[4], ref_vector.elements_arbiter().attach(raw_offsets[4]*sizeof(DT_), raw_sizes[4]*sizeof(DT_)));
    DVB1 ref_dvb1_2(sizes[5], ref_vector.elements_arbiter().attach(raw_offsets[5]*sizeof(DT_), raw_sizes[5]*sizeof(DT_)));
    DV   ref_dv    (sizes[6], ref_vector.elements_arbiter().attach(raw_offsets[6]*sizeof(DT_), raw_sizes[6]*sizeof(DT_)));

    // copy values into our tuple vector
    ref_vector.copy_to(tv);

    // now compare all individual sub-vectors
    TEST_CHECK_LESS_THAN(tv_dvb3.max_rel_diff(ref_dvb3), tol);
    TEST_CHECK_LESS_THAN(tv_dvb2_0.max_rel_diff(ref_dvb2_0), tol);
    TEST_CHECK_LESS_THAN(tv_dvb2_1.max_rel_diff(ref_dvb2_1), tol);
    TEST_CHECK_LESS_THAN(tv_dvb1_0.max_rel_diff(ref_dvb1_0), tol);
    TEST_CHECK_LESS_THAN(tv_dvb1_1.max_rel_diff(ref_dvb1_1), tol);
    TEST_CHECK_LESS_THAN(tv_dvb1_2.max_rel_diff(ref_dvb1_2), tol);
    TEST_CHECK_LESS_THAN(tv_dv.max_rel_diff(ref_dv), tol);

    // create a new dense vector, copy from our tuple vector and compare it to the reference
    DV dense_vec;
    dense_vec.convert(tv);

    TEST_CHECK(dense_vec.same_layout(ref_vector));
    TEST_CHECK_LESS_THAN(dense_vec.max_rel_diff(ref_vector), tol);
  }
};

SPAWN_UNIT_TEST_2T_P(MetaVectorConvertTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaVectorConvertTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaVectorConvertTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaVectorConvertTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(MetaVectorConvertTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(MetaVectorConvertTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(MetaVectorConvertTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaVectorConvertTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(MetaVectorConvertTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaVectorConvertTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(MetaVectorConvertTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MetaVectorConvertTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MetaVectorConvertTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MetaVectorConvertTest, double, std::uint64_t, PreferredBackend::cuda);
#endif
