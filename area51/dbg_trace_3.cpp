#include <kernel/base_header.hpp>
#include <kernel/util/runtime.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/space/ext_vtk_writer.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/assembly/trace_assembler.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/geometry/mesh_file_writer.hpp>

// We are using FEAT
using namespace FEAT;

namespace DbgTrace3
{
  template<int i_>
  class NiFunctional : public Assembly::LinearFunctional
  {
  public:
    static constexpr TrafoTags trafo_config = TrafoTags::normal;
    static constexpr SpaceTags test_config = SpaceTags::value;
    template<typename AsmTraits_>
    class Evaluator :
      public Assembly::LinearFunctional::template Evaluator<AsmTraits_>
    {
    public:
      /// data type
      typedef typename AsmTraits_::DataType DataType;
      /// trafo data type
      typedef typename AsmTraits_::TrafoData TrafoData;
      /// test function data type
      typedef typename AsmTraits_::TestBasisData TestBasisData;

      typedef DataType ValueType;

      DataType ni;

    public:
      explicit Evaluator(const NiFunctional&) {}
      void set_point(const TrafoData& tau)
      {
        ni = tau.normal[i_];
      }

      ValueType eval(const TestBasisData& psi) const
      {
        return ni * psi.value;
      }
    };
  };

  template<int dim_>
  class NFunctional : public Assembly::LinearFunctional
  {
  public:
    static constexpr TrafoTags trafo_config = TrafoTags::normal;
    static constexpr SpaceTags test_config = SpaceTags::value;
    template<typename AsmTraits_>
    class Evaluator :
      public Assembly::LinearFunctional::template Evaluator<AsmTraits_>
    {
    public:
      /// data type
      typedef typename AsmTraits_::DataType DataType;
      /// trafo data type
      typedef typename AsmTraits_::TrafoData TrafoData;
      /// test function data type
      typedef typename AsmTraits_::TestBasisData TestBasisData;

      typedef Tiny::Vector<DataType, dim_> ValueType;
      ValueType n;

    public:
      explicit Evaluator(const NFunctional&) {}
      void set_point(const TrafoData& tau)
      {
        n = tau.normal;
      }

      ValueType eval(const TestBasisData& psi) const
      {
        return psi.value * n;
      }
    };
  };

  template<typename ShapeType>
  void run(const String shape_name, Index level)
  {
    //typedef Shape::Quadrilateral ShapeType;
    typedef double DataType;
    typedef Index IndexType;
    typedef Mem::Main MemType;

    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Geometry::MeshPart<MeshType> BoundaryType;
    typedef Geometry::RefinedUnitCubeFactory<MeshType> MeshFactoryType;
    typedef Geometry::BoundaryFactory<MeshType> BoundaryFactoryType;

    static constexpr int dim = ShapeType::dimension;

    //Index level = 0;

    std::cout << "Creating Mesh on Level " << level << "..." << std::endl;
    MeshFactoryType mesh_factory(level);
    MeshType mesh(mesh_factory);

    std::cout << "Creating Boundary..." << std::endl;
    BoundaryFactoryType boundary_factory(mesh);
    BoundaryType boundary(boundary_factory);

    std::cout << "Creating Trafo..." << std::endl;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    TrafoType trafo(mesh);

    std::cout << "Creating Space..." << std::endl;
    typedef Space::Discontinuous::ElementP0<TrafoType> SpaceType;
    SpaceType space(trafo);

    //typedef LAFEM::DenseVector<MemType, DataType, IndexType> VectorType;
    typedef LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, dim> BVectorType;


    std::cout << "Allocating and initializing vectors and matrix..." << std::endl;
    BVectorType vec_n(space.get_num_dofs());
    /*VectorType vec_x(space.get_num_dofs());
    VectorType vec_y(space.get_num_dofs());
    VectorType vec_z(space.get_num_dofs());*/
    vec_n.format();
    /*vec_x.format();
    vec_y.format();
    vec_z.format();*/

    Cubature::DynamicFactory cubature_factory("auto-degree:1");


    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    {
      NFunctional<dim> nf;
      /*NiFunctional<0> nix;
      NiFunctional<1> niy;
      NiFunctional<2> niz;*/
      Assembly::TraceAssembler<TrafoType> trace_asm(trafo);
      trace_asm.add_mesh_part(boundary);
      trace_asm.compile();
      trace_asm.assemble_functional_vector(vec_n, nf, space, cubature_factory);
      /*trace_asm.assemble_functional_vector(vec_x, nix, space, cubature_factory);
      trace_asm.assemble_functional_vector(vec_y, niy, space, cubature_factory);
      if(ShapeType::dimension > 2)
        trace_asm.assemble_functional_vector(vec_z, niz, space, cubature_factory);*/
    }

    /*{
      std::ofstream ofs("dbg-trace-3-" + shape_name + "-lvl" + stringify(level) + ".xml");
      Geometry::MeshFileWriter writer(ofs);
      writer.write_mesh(mesh);
    }*/

    // First of all, build the filename string
    String vtk_name(String("./dbg-trace-3-") + shape_name + "-lvl" + stringify(level));

    std::cout << "Writing VTK file '" << vtk_name << ".vtu'..." << std::endl;

    Space::ExtVtkWriter<TrafoType> ext_vtk(trafo, 0);

    ext_vtk.open("dbg-trace-" + shape_name + "-" + stringify(level) + ".vtk");
    ext_vtk.write_values_blocked("n", space, vec_n);
    /*ext_vtk.write_values("nx", space, vec_x.elements());
    ext_vtk.write_values("ny", space, vec_y.elements());
    if(ShapeType::dimension > 2)
      ext_vtk.write_values("nz", space, vec_z.elements());*/
    ext_vtk.close();

    /*
    // Create a VTK exporter for our mesh
    Geometry::ExportVTK<MeshType> exporter(mesh);

    // add the vertex-projection of our solution and rhs vectors
    exporter.add_vertex_scalar("nx", vec_x.elements());
    exporter.add_vertex_scalar("ny", vec_y.elements());
    if(ShapeType::dimension > 2)
      exporter.add_vertex_scalar("nz", vec_z.elements());

    // finally, write the VTK file
    exporter.write(vtk_name);
    */

    // That's all, folks.
    std::cout << "Finished!" << std::endl;
  }
  void main()
  {
    //run<Shape::Quadrilateral>("quad");
    //run<Shape::Triangle>("tria", 0);
    //run<Shape::Triangle>("tria", 1);
    //run<Shape::Triangle>("tria", 2);
    run<Shape::Tetrahedron>("tetra", 0);
    run<Shape::Tetrahedron>("tetra", 1);
    run<Shape::Tetrahedron>("tetra", 2);
    //run<Shape::Hexahedron>("hexa");
  }
}

int main(int argc, char** argv)
{
  FEAT::Runtime::initialize(argc, argv);
  DbgTrace3::main();
  return FEAT::Runtime::finalize();
}