#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/rannacher_turek/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/common_functions.hpp>

using namespace FEAST;

typedef Geometry::ConformalMesh<Shape::Quadrilateral> QuadMesh;
typedef Geometry::RefinedUnitCubeFactory<QuadMesh> QuadMeshFactory;
typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

typedef Space::Lagrange1::Element<QuadTrafo> QuadSpaceQ1;
typedef Space::RannacherTurek::Element<QuadTrafo> QuadSpaceQ1T;
typedef Space::Discontinuous::Element<QuadTrafo> QuadSpaceQ0;

template<typename Space_>
void test_interpolation(Index level)
{
  // define mesh factory
  QuadMeshFactory mesh_factory(level);

  // create mesh
  QuadMesh mesh(mesh_factory);

  // create trafo
  QuadTrafo trafo(mesh);

  // create space
  Space_ space(trafo);

  // define functor
  Assembly::Common::SineBubbleFunction sine_bubble;

  // interpolate functor into FE space
  LAFEM::DenseVector<Mem::Main, double> vector;
  Assembly::Interpolator::project(vector, sine_bubble, space);

  // compute L2-error
  Cubature::DynamicFactory cubature_factory("auto-degree:10");
  double l2err = Assembly::ScalarErrorComputerL2::compute(vector, sine_bubble, space, cubature_factory);
  double h1err = Assembly::ScalarErrorComputerH1::compute(vector, sine_bubble, space, cubature_factory);

  // print error
  std::cout << "Level: " << level <<
    " , L2-Error: " << std::scientific << l2err <<
    " , H1-Error: " << std::scientific << h1err << std::endl;
}

int main(int /*argc*/, char** /*argv*/)
{
  for(Index i(0); i < 5; ++i)
  {
    test_interpolation<QuadSpaceQ1>(i);
  }

/*
  // define mesh factory
  QuadMeshFactory mesh_factory(5);

  // create mesh
  QuadMesh mesh(mesh_factory);

  // create trafo
  QuadTrafo trafo(mesh);

  // create space
  //typedef QuadSpaceQ0 QuadSpace;
  typedef QuadSpaceQ1 QuadSpace;
  //typedef QuadSpaceQ1T QuadSpace;
  QuadSpace space(trafo);

  // define functor
  //typedef Space::DeriveFunctorStaticWrapper<SineBubble> FunctorType;
  typedef Analytic::StaticWrapperFunctor<SineBubble> FunctorType;
  FunctorType functor;

  // interpolate functor into FE space
  VectorType vector;
  Assembly::Interpolator::project(vector, functor, space);

  // project to vertices
  VectorType vertex_vector;
  Assembly::DiscreteVertexProjector::project(vertex_vector, vector, space);

  // project to cells
  VectorType cell_vector;
  Assembly::DiscreteCellProjector::project(cell_vector, vector, space, "auto-degree:5");

  //
  Geometry::ExportVTK<QuadMesh> exporter(mesh);
  exporter.open("test.vtk");
  exporter.write_vertex_scalar("test-vtx", vertex_vector.elements());
  exporter.write_cell_scalar("test-cell", cell_vector.elements());
  exporter.close();
  */

  return 0;
}
