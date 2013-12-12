#include <kernel/geometry/reference_cell_factory.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/ext_vtk_writer.hpp>

// FE spaces
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/rannacher_turek/element.hpp>
#include <kernel/space/bogner_fox_schmit/element.hpp>
#include <kernel/space/hermite3/element.hpp>
#include <kernel/space/argyris/element.hpp>

#include <vector>

using namespace FEAST;

template<typename SpaceType_>
void dump_basis(String vtk_name, Index num_refines = 5);


int main(int, char**)
{
  // 1D Line Elements
  {
    typedef Geometry::ConformalMesh< Shape::Hypercube<1> > MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;

    // Discontinuous-0
    dump_basis< Space::Discontinuous::Element<TrafoType> >("1d_line_discontinuous-0.vtk");
    // Lagrange-1
    dump_basis< Space::Lagrange1::Element<TrafoType> >("1d_line_lagrange-1.vtk");
    // Bogner-Fox-Schmit
    dump_basis< Space::BognerFoxSchmit::Element<TrafoType> >("1d_line_bogner_fox_schmit.vtk");
  }

  // 2D Tria Elements
  {
    typedef Geometry::ConformalMesh< Shape::Simplex<2> > MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;

    // Discontinuous-0
    dump_basis< Space::Discontinuous::Element<TrafoType> >("2d_tria_discontinuous-0.vtk");
    // Lagrange-1
    dump_basis< Space::Lagrange1::Element<TrafoType> >("2d_tria_lagrange-1.vtk");
    // Hermite-3
    dump_basis< Space::Hermite3::Element<TrafoType> >("2d_tria_hermite-3.vtk");
    // Argyris
    dump_basis< Space::Argyris::Element<TrafoType> >("2d_tria_argyris.vtk");
  }

  // 3D Tetra Elements
  {
    //typedef Geometry::ConformalMesh< Shape::Simplex<3> > MeshType;
    //typedef Trafo::Standard::Mapping<MeshType> TrafoType;

    // Discontinuous-0
    //dump_basis< Space::Discontinuous::Element<TrafoType> >("3d_tetra_discontinuous-0.vtk");
    // Lagrange-1
    //dump_basis< Space::Lagrange1::Element<TrafoType> >("3d_tetra_lagrange-1.vtk");
  }

  // 2D Quad Elements
  {
    typedef Geometry::ConformalMesh< Shape::Hypercube<2> > MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;

    // Discontinuous-0
    dump_basis< Space::Discontinuous::Element<TrafoType> >("2d_quad_discontinuous-0.vtk");
    // Lagrange-1
    dump_basis< Space::Lagrange1::Element<TrafoType> >("2d_quad_lagrange-1.vtk");
    // Rannacher-Turek
    dump_basis< Space::RannacherTurek::Element<TrafoType> >("2d_quad_rannacher_turek.vtk");
    // Bogner-Fox-Schmit
    dump_basis< Space::BognerFoxSchmit::Element<TrafoType> >("2d_quad_bogner_fox_schmit.vtk");
    // Hermite-3
    dump_basis< Space::Hermite3::Element<TrafoType> >("2d_quad_hermite-3.vtk");
  }

  // 3D Hexa Elements
  {
    //typedef Geometry::ConformalMesh< Shape::Hypercube<3> > MeshType;
    //typedef Trafo::Standard::Mapping<MeshType> TrafoType;

    // Discontinuous-0
    //dump_basis< Space::Discontinuous::Element<TrafoType> >("3d_hexa_discontinuous-0.vtk");
    // Lagrange-1
    //dump_basis< Space::Lagrange1::Element<TrafoType> >("3d_hexa_lagrange-1.vtk");
  }
}

template<bool _enable>
struct DumpWrapper
{
  template<typename Writer_, typename Space_>
  static void write_values(Writer_&, Space_&) {}
  template<typename Writer_, typename Space_>
  static void write_gradients(Writer_&, Space_&) {}
  template<typename Writer_, typename Space_>
  static void write_hessians(Writer_&, Space_&) {}
};

template<>
struct DumpWrapper<true>
{
  template<typename Writer_, typename Space_>
  static void write_values(Writer_& writer, Space_& space)
  {
    // create a dof-vector
    Index num_dofs = space.get_num_dofs();
    std::vector<double> v(num_dofs, 0.0);

    for(Index i(0); i < num_dofs; ++i)
    {
      if(i > Index(0))
        v[i-1] = 0.0;
      v[i] = 1.0;
      writer.write_values(String("phi_") + (i < 10 ? "0" : "") + stringify(i), space, v.data());
    }
  }

  template<typename Writer_, typename Space_>
  static void write_gradients(Writer_& writer, Space_& space)
  {
    // create a dof-vector
    Index num_dofs = space.get_num_dofs();
    std::vector<double> v(num_dofs, 0.0);

    for(Index i(0); i < num_dofs; ++i)
    {
      if(i > Index(0))
        v[i-1] = 0.0;
      v[i] = 1.0;
      writer.write_gradients(String("phi_grad_") + (i < 10 ? "0" : "") + stringify(i), space, v.data());
    }
  }

  template<typename Writer_, typename Space_>
  static void write_hessians(Writer_& writer, Space_& space)
  {
    // create a dof-vector
    Index num_dofs = space.get_num_dofs();
    std::vector<double> v(num_dofs, 0.0);

    for(Index i(0); i < num_dofs; ++i)
    {
      if(i > Index(0))
        v[i-1] = 0.0;
      v[i] = 1.0;
      writer.write_hessians(String("phi_hess_") + (i < 10 ? "0" : "") + stringify(i), space, v.data());
    }
  }
};

template<typename SpaceType_>
void dump_basis(String vtk_name, Index num_refines)
{
  typedef SpaceType_ SpaceType;
  typedef typename SpaceType_::TrafoType TrafoType;
  typedef typename TrafoType::MeshType MeshType;
  typedef typename MeshType::ShapeType ShapeType;

  // create a mesh, trafo and space
  Geometry::ReferenceCellFactory<ShapeType> factory;
  MeshType mesh(factory);
  TrafoType trafo(mesh);
  SpaceType space(trafo);

  // print a message
  std::cout << "Writing '" << vtk_name << "'..." << std::endl;

  // create an extended VTK writer
  Space::ExtVtkWriter<TrafoType> vtk_writer(trafo, num_refines);
  vtk_writer.open(vtk_name);

  // write basis functions
  typedef typename TrafoType::template Evaluator<>::Type TrafoEval;
  typedef typename SpaceType::template Evaluator<TrafoEval>::Type SpaceEval;
  DumpWrapper<SpaceEval::can_value != 0>::write_values(vtk_writer, space);
  DumpWrapper<SpaceEval::can_grad  != 0>::write_gradients(vtk_writer, space);
  DumpWrapper<SpaceEval::can_hess  != 0>::write_hessians(vtk_writer, space);

  // okay
  vtk_writer.close();
}
