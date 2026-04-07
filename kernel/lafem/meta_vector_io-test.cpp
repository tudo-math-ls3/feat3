// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/lafem/meta_vector_test_base.hpp>
#include <kernel/util/binary_stream.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Meta-Vector scale test class
 *
 * \test The 'write/read' operation of the DenseVector, PowerVector and TupleVector class templates.
 *
 * \author Dirk Ribbrock
 */
template<
  typename DataType_,
  typename IndexType_>
class MetaVectorIOTest
  : public MetaVectorTestBase<DataType_, IndexType_>
{
public:
  typedef DataType_ DataType;
  typedef MetaVectorTestBase<DataType_, IndexType_> BaseClass;
  typedef typename BaseClass::MetaVector MetaVector;

  explicit MetaVectorIOTest(PreferredBackend backend) :
    BaseClass("MetaVectorIOTest", Type::Traits<DataType>::name(), Type::Traits<IndexType_>::name(), backend)
  {
  }

  using BaseClass::fx00;
  using BaseClass::fx01;
  using BaseClass::fx1;

  virtual void run() const override
  {
    const DataType tol = TestSystem::tol<DataType>();

    const Index n00 = 5;
    const Index n01 = 10;
    const Index n1 = 7;

    MetaVector x(this->gen_vector_x(n00, n01, n1));
    MetaVector z(this->gen_vector_null(n00, n01, n1));

    BinaryStream bs;
    x.write_out(FileMode::fm_binary, bs);
    bs.seekg(0);
    z.read_from(FileMode::fm_binary, bs);

    TEST_CHECK_LESS_THAN(z.max_rel_diff(x), tol);
  }
};

SPAWN_UNIT_TEST_2T_P(MetaVectorIOTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaVectorIOTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaVectorIOTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaVectorIOTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(MetaVectorIOTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(MetaVectorIOTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
//#ifdef FEAT_HAVE_QUADMATH
//SPAWN_UNIT_TEST_2T_P(MetaVectorIOTest, __float128, std::uint32_t, PreferredBackend::generic);
//SPAWN_UNIT_TEST_2T_P(MetaVectorIOTest, __float128, std::uint64_t, PreferredBackend::generic);
//#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(MetaVectorIOTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaVectorIOTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(MetaVectorIOTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MetaVectorIOTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MetaVectorIOTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MetaVectorIOTest, double, std::uint64_t, PreferredBackend::cuda);
#endif
