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

using namespace FEAST;

typedef double DataType;

typedef Geometry::ConformalMesh<Shape::Quadrilateral> QuadMesh;
typedef Geometry::RefinedUnitCubeFactory<QuadMesh> QuadMeshFactory;
typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

typedef Space::Lagrange1::Element<QuadTrafo> QuadSpaceQ1;
typedef Space::RannacherTurek::Element<QuadTrafo> QuadSpaceQ1T;
typedef Space::Discontinuous::Element<QuadTrafo> QuadSpaceQ0;

typedef LAFEM::DenseVector<Mem::Main, DataType> VectorType;

template<typename DataType_>
class SineBubble :
  public Analytic::StaticFunction<DataType_>
{
public:
  /// returns the constant pi
  static DataType_ pi()
  {
    static const DataType_ value(DataType_(2) * std::acos(DataType_(0)));
    return value;
  }

  /// 1D: function value
  static DataType_ eval(DataType_ x)
  {
    return std::sin(pi() * x);
  }

  /// 2D: function value
  static DataType_ eval(DataType_ x, DataType_ y)
  {
    return std::sin(pi() * x) * std::sin(pi() * y);
  }

  /// 3D: function value
  static DataType_ eval(DataType_ x, DataType_ y, DataType_ z)
  {
    return std::sin(pi() * x) * std::sin(pi() * y) * std::sin(pi() * z);
  }

  /// 1D: X-derivative
  static DataType_ der_x(DataType_ x)
  {
    return pi() * std::cos(pi() * x);
  }

  /// 2D: X-derivative
  static DataType_ der_x(DataType_ x, DataType_ y)
  {
    return pi() * std::cos(pi() * x) * std::sin(pi() * y);
  }

  /// 2D: Y-derivative
  static DataType_ der_y(DataType_ x, DataType_ y)
  {
    return pi() * std::sin(pi() * x) * std::cos(pi() * y);
  }

  /// 3D: X-derivative
  static DataType_ der_x(DataType_ x, DataType_ y, DataType_ z)
  {
    return pi() * std::cos(pi() * x) * std::sin(pi() * y) * std::sin(pi() * z);
  }

  /// 3D: Y-derivative
  static DataType_ der_y(DataType_ x, DataType_ y, DataType_ z)
  {
    return pi() * std::sin(pi() * x) * std::cos(pi() * y) * std::sin(pi() * z);
  }

  /// 3D: Z-derivative
  static DataType_ der_z(DataType_ x, DataType_ y, DataType_ z)
  {
    return pi() * std::sin(pi() * x) * std::sin(pi() * y) * std::cos(pi() * z);
  }
}; // class SineBubble<...>

int main(int argc, char* argv[])
{
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

  //*
  Geometry::ExportVTK<QuadMesh> exporter(mesh);
  exporter.open("test.vtk");
  exporter.write_vertex_scalar("test-vtx", vertex_vector.elements());
  exporter.write_cell_scalar("test-cell", cell_vector.elements());
  exporter.close();
  //*/

  return 0;
}