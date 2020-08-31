// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
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
 * \author Peter Zajac
 */
template<typename DataType_>
class InterpolatorTest :
  public TestSystem::TaggedTest<Archs::None, DataType_>
{
  typedef LAFEM::DenseVector<Mem::Main, DataType_> VectorType;

  typedef Geometry::ConformalMesh<Shape::Quadrilateral> QuadMesh;

  typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

  typedef Space::Discontinuous::Element<QuadTrafo> QuadSpaceQ0;
  typedef Space::Lagrange1::Element<QuadTrafo> QuadSpaceQ1;

public:
  InterpolatorTest() :
    TestSystem::TaggedTest<Archs::None, DataType_>("InterpolatorTest")
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
    for(Index i(0); i < num_verts; ++i)
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
    const Geometry::IndexSet<4>& index_set(mesh.get_index_set<2,0>());

    // loop over all quads
    for(Index i(0); i < num_quads; ++i)
    {
      // compute quad center
      DataType_ x(0), y(0);
      for(int j(0); j < 4; ++j)
      {
        x += DataType_(vertex_set[index_set[i][j]][0]);
        y += DataType_(vertex_set[index_set[i][j]][1]);
      }

      // compute sine-bubble value in quad center
      DataType_ s = Analytic::Common::SineBubbleStatic<DataType_>::eval(DataType_(0.25)*x, DataType_(0.25)*y);

      // validate vector data
      TEST_CHECK_EQUAL_WITHIN_EPS(vector(i), s, eps);
    }
  }
};

InterpolatorTest<float> interpolator_test_float;
InterpolatorTest<double> interpolator_test_double;
