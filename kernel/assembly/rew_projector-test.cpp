// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/assembly/rew_projector.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/space/cro_rav_ran_tur/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/math.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

/**
 * \brief RewProjector test class template
 *
 * \test Tests the Assembly::RewProjector class
 *
 * \tparam DataType_
 * The data type for the test. Shall be either double or float.
 *
 * \tparam IndexType_
 * The index type for the test. Shall be either unsigned int or unsigned long.
 *
 * \author Peter Zajac
 */
template<typename DataType_, typename IndexType_>
class RewProjectorTest :
  public UnitTest
{
  typedef LAFEM::DenseVector<DataType_, IndexType_> VectorType;

  typedef Geometry::ConformalMesh<Shape::Quadrilateral> QuadMesh;

  typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

  typedef Space::Discontinuous::Element<QuadTrafo> QuadSpaceQ0;
  typedef Space::Lagrange1::Element<QuadTrafo> QuadSpaceQ1;
  typedef Space::Lagrange2::Element<QuadTrafo> QuadSpaceQ2;
  typedef Space::CroRavRanTur::Element<QuadTrafo> QuadSpaceRT;

public:
  explicit RewProjectorTest(PreferredBackend backend) :
    UnitTest("RewProjectorTest", Type::Traits<DataType_>::name(), Type::Traits<IndexType_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    test_unit_2d();
  }

  void test_unit_2d() const
  {
    // compute eps
    const DataType_ tol = TestSystem::relaxed_tol<DataType_>();

    // create coarse mesh
    Geometry::RefinedUnitCubeFactory<QuadMesh> unit_factory(3);
    QuadMesh mesh(unit_factory);

    // create trafo
    QuadTrafo trafo(mesh);

    // project Q0
    DataType_ q0_err = project<QuadSpaceQ0>(trafo, "gauss-legendre:2");
    q0_err = q0_err / DataType_(7.973149267259641e-02) - DataType_(1);
    TEST_CHECK_EQUAL_WITHIN_EPS(q0_err, DataType_(0), tol);

    // project Q1
    DataType_ q1_err = project<QuadSpaceQ1>(trafo, "gauss-legendre:3");
    q1_err = q1_err / DataType_(4.143030804263830e-03) - DataType_(1);
    TEST_CHECK_EQUAL_WITHIN_EPS(q1_err, DataType_(0), tol);

    // project Q2
    DataType_ q2_err = project<QuadSpaceQ2>(trafo, "gauss-legendre:4");
    q2_err = q2_err / DataType_(2.157823711367303e-04) - DataType_(1);
    TEST_CHECK_EQUAL_WITHIN_EPS(q2_err, DataType_(0), tol);

    // project RT
    DataType_ rt_err = project<QuadSpaceRT>(trafo, "gauss-legendre:3");
    rt_err = rt_err / DataType_(7.5903151744156896e-03) - DataType_(1);
    TEST_CHECK_EQUAL_WITHIN_EPS(rt_err, DataType_(0), tol);
  }

  template<typename Space_>
  DataType_ project(QuadTrafo& trafo, String cubature_name) const
  {
    // create space
    Space_ space(trafo);

    // define function
    Analytic::Common::SineBubbleFunction<2> function;

    // define a cubature factory
    Cubature::DynamicFactory cubature_factory(cubature_name);

    // project function into FE space
    VectorType vector;
    Assembly::RewProjector::project(vector, function, space, cubature_factory);

    // compute L2-Error
    return Assembly::ScalarErrorComputer<0>::compute(vector, function, space, cubature_factory).norm_h0;
  }

};

SPAWN_UNIT_TEST_2T_P(RewProjectorTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(RewProjectorTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(RewProjectorTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(RewProjectorTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(RewProjectorTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(RewProjectorTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(RewProjectorTest, __float128, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(RewProjectorTest, __float128, std::uint64_t, PreferredBackend::generic);
#endif
//#ifdef FEAT_HAVE_HALFMATH
// Disabled: tol too sharp, error is about 3 times larger than even relaxed tolerance
//SPAWN_UNIT_TEST_2T_P(RewProjectorTest, Half, std::uint32_t, PreferredBackend::generic);
// Disabled: tol too sharp, error is about 3 times larger than even relaxed tolerance
//SPAWN_UNIT_TEST_2T_P(RewProjectorTest, Half, std::uint64_t, PreferredBackend::generic);
//#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(RewProjectorTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(RewProjectorTest, double, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(RewProjectorTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(RewProjectorTest, double, std::uint64_t, PreferredBackend::cuda);
#endif
