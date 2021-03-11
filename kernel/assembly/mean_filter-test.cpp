// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/assembly/mean_filter_assembler.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/mean_filter.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/random.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

template<typename DataType_>
class MeanFilterTest :
  public TestSystem::TaggedTest<Archs::None, DataType_>
{
  typedef LAFEM::MeanFilter<Mem::Main, DataType_, Index> FilterType;
  typedef LAFEM::DenseVector<Mem::Main, DataType_, Index> VectorType;

  typedef Geometry::ConformalMesh<Shape::Quadrilateral> QuadMesh;

  typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

  typedef Space::Lagrange2::Element<QuadTrafo> QuadSpaceQ2;

public:
  MeanFilterTest() :
    TestSystem::TaggedTest<Archs::None, DataType_>("MeanFilterTest")
  {
  }

  virtual ~MeanFilterTest()
  {
  }

  virtual void run() const override
  {
    test_unit_2d();
  }

  void test_unit_2d() const
  {
    const DataType_ tol = Math::pow(Math::eps<DataType_>(), DataType_(0.75));
    // create coarse mesh
    Geometry::RefinedUnitCubeFactory<QuadMesh> unit_factory(3);
    QuadMesh mesh(unit_factory);
    QuadTrafo trafo(mesh);
    QuadSpaceQ2 space(trafo);

    // create a cubature factory
    Cubature::DynamicFactory cubature_factory("gauss-legendre:3");

    FilterType filter;

    VectorType vec_prim, vec_dual;
    Assembly::MeanFilterAssembler::assemble(vec_prim, vec_dual, space, cubature_factory);


    Assembly::MeanFilterAssembler::assemble(filter, space, cubature_factory, DataType_(0.5));

    Random rng;

    VectorType vec_sol(rng, space.get_num_dofs(), DataType_(-2.0), DataType_(1.0));
    VectorType vec_rhs(rng, space.get_num_dofs(), DataType_(-2.0), DataType_(1.0));
    VectorType vec_def(rng, space.get_num_dofs(), DataType_(-2.0), DataType_(1.0));
    VectorType vec_cor(rng, space.get_num_dofs(), DataType_(-2.0), DataType_(1.0));

    // apply filter to all vectors
    filter.filter_sol(vec_sol);
    filter.filter_cor(vec_cor);
    filter.filter_rhs(vec_rhs);
    filter.filter_def(vec_def);

    // check primal integrals
    DataType_ int_sol = vec_sol.dot(vec_dual);
    TEST_CHECK_EQUAL_WITHIN_EPS(int_sol, DataType_(0.5), tol);
    DataType_ int_cor = vec_cor.dot(vec_dual);
    TEST_CHECK_EQUAL_WITHIN_EPS(int_cor, DataType_(0.0), tol);

    DataType_ int_rhs = vec_rhs.dot(vec_prim);
    TEST_CHECK_EQUAL_WITHIN_EPS(int_rhs, DataType_(0.0), tol);
    DataType_ int_def = vec_def.dot(vec_prim);
    TEST_CHECK_EQUAL_WITHIN_EPS(int_def, DataType_(0.0), tol);
  }
};

MeanFilterTest<double> mean_filter_test_double;
