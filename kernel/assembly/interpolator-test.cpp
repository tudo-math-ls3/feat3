// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/math.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

/**
 * \brief Interpolator test class template
 *
 * \test Tests the Assembly::Interpolator class
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
class InterpolatorTest :
  public UnitTest
{
  typedef LAFEM::DenseVector<DataType_> VectorType;

  typedef Geometry::ConformalMesh<Shape::Quadrilateral> QuadMesh;

  typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

  typedef Space::Discontinuous::Element<QuadTrafo> QuadSpaceQ0;
  typedef Space::Lagrange1::Element<QuadTrafo> QuadSpaceQ1;

public:
  InterpolatorTest(PreferredBackend backend) :
    UnitTest("InterpolatorTest", Type::Traits<DataType_>::name(), Type::Traits<IndexType_>::name(), backend)
  {
  }

  virtual ~InterpolatorTest()
  {
  }

  virtual void run() const override
  {
    test_unit_2d();
  }

  void test_unit_2d() const
  {
    // create coarse mesh
    Geometry::RefinedUnitCubeFactory<QuadMesh> unit_factory(3);
    QuadMesh mesh(unit_factory);

    // run tests
    test_unit_2d_q1(mesh);
    test_unit_2d_q0(mesh);
  }

  void test_unit_2d_q1(QuadMesh& mesh) const
  {
    // compute eps
    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.8));

    // create trafo
    QuadTrafo trafo(mesh);

    // create space
    QuadSpaceQ1 space(trafo);

    // define function
    Analytic::Common::SineBubbleFunction<2> sine_bubble;

    // interpolate function into FE space
    VectorType vector;
    Assembly::Interpolator::project(vector, sine_bubble, space);

    // get mesh vertex count
    const Index num_verts = mesh.get_num_entities(0);

    // validate vector size
    TEST_CHECK_EQUAL(vector.size(), num_verts);

    // get the vertex set of the mesh
    const typename QuadMesh::VertexSetType& vertex_set(mesh.get_vertex_set());

    // loop over all vertices
    for (Index i(0); i < num_verts; ++i)
    {
      // compute sine-bubble value in vertex position
      DataType_ s = Analytic::Common::SineBubbleStatic<DataType_>
        ::eval(DataType_(vertex_set[i][0]), DataType_(vertex_set[i][1]));

      // validate vector data
      TEST_CHECK_EQUAL_WITHIN_EPS(vector(i), s, eps);
    }
  }

  void test_unit_2d_q0(QuadMesh& mesh) const
  {
    // compute eps
    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.8));

    // create trafo
    QuadTrafo trafo(mesh);

    // create space
    QuadSpaceQ0 space(trafo);

    // define function
    Analytic::Common::SineBubbleFunction<2> sine_bubble;

    // interpolate function into FE space
    VectorType vector;
    Assembly::Interpolator::project(vector, sine_bubble, space);

    // get mesh element count
    const Index num_quads = mesh.get_num_entities(2);

    // validate vector size
    TEST_CHECK_EQUAL(vector.size(), num_quads);

    // get the vertex set of the mesh
    const typename QuadMesh::VertexSetType& vertex_set(mesh.get_vertex_set());

    // get the index set of the mesh
    const Geometry::IndexSet<4>& index_set(mesh.get_index_set<2, 0>());

    // loop over all quads
    for (Index i(0); i < num_quads; ++i)
    {
      // compute quad center
      DataType_ x(0), y(0);
      for (int j(0); j < 4; ++j)
      {
        x += DataType_(vertex_set[index_set[i][j]][0]);
        y += DataType_(vertex_set[index_set[i][j]][1]);
      }

      // compute sine-bubble value in quad center
      DataType_ s = Analytic::Common::SineBubbleStatic<DataType_>::eval(DataType_(0.25) * x, DataType_(0.25) * y);

      // validate vector data
      TEST_CHECK_EQUAL_WITHIN_EPS(vector(i), s, eps);
    }
  }
};

InterpolatorTest <float, std::uint32_t> interpolator_test_float_uint32(PreferredBackend::generic);
InterpolatorTest <double, std::uint32_t> interpolator_test_double_uint32(PreferredBackend::generic);
InterpolatorTest <float, std::uint64_t> interpolator_test_float_uint64(PreferredBackend::generic);
InterpolatorTest <double, std::uint64_t> interpolator_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
InterpolatorTest <float, std::uint64_t> mkl_interpolator_test_float_uint64(PreferredBackend::mkl);
InterpolatorTest <double, std::uint64_t> mkl_interpolator_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
InterpolatorTest <__float128, std::uint32_t> interpolator_test_float128_uint32(PreferredBackend::generic);
InterpolatorTest <__float128, std::uint64_t> interpolator_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
InterpolatorTest <Half, std::uint32_t> interpolator_test_half_uint32(PreferredBackend::generic);
InterpolatorTest <Half, std::uint64_t> interpolator_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
InterpolatorTest <float, std::uint32_t> cuda_interpolator_test_float_uint32(PreferredBackend::cuda);
InterpolatorTest <double, std::uint32_t> cuda_interpolator_test_double_uint32(PreferredBackend::cuda);
InterpolatorTest <float, std::uint64_t> cuda_interpolator_test_float_uint64(PreferredBackend::cuda);
InterpolatorTest <double, std::uint64_t> cuda_interpolator_test_double_uint64(PreferredBackend::cuda);
#endif
