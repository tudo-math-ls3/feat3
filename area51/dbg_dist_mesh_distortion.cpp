// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/mesh_distortion.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/util/dist.hpp>
#include <control/domain/parti_domain_control.hpp>

#include <vector>

namespace DbgDistMeshDisto
{
  using namespace FEAT;

  template<typename MeshType_>
  std::vector<double> analyze_mesh(MeshType_& mesh)
  {
    typedef Trafo::Standard::Mapping<MeshType_> TrafoType;
    TrafoType trafo(mesh);
    typedef typename TrafoType::ShapeType ShapeType;
    static constexpr int shape_dim = ShapeType::dimension;
    const int num_verts = Shape::FaceTraits<ShapeType, 0>::count;

    // create trafo evaluator
    typename TrafoType::template Evaluator<>::Type trafo_eval(trafo);
    const Index num_elems = trafo_eval.get_num_cells();

    Tiny::Matrix<double, shape_dim, shape_dim> jacmat;
    Tiny::Vector<double, shape_dim> vertex;
    std::vector<double> res(num_elems, 0.0);

    // loop over all elements
    for(Index ielem(0); ielem < num_elems; ++ielem)
    {
      trafo_eval.prepare(ielem);

      double min_jac = 1E+99;
      double max_jac = 1E-99;

      // loop over all corner vertices
      for(int iv(0); iv < num_verts; ++iv)
      {
        // set vertex coords
        for(int k(0); k < shape_dim; ++k)
          vertex[k] = Shape::ReferenceCell<ShapeType>::template vertex<double>(iv, k);

        // evaluate jacobian matrix
        trafo_eval.calc_jac_mat(jacmat, vertex);

        // compute determinant = volume
        double jac_def = jacmat.det();

        // update min/max
        Math::minimax(jac_def, min_jac, max_jac);
      }

      if(min_jac > +1E-10)
        res[ielem] = min_jac / max_jac;
      else if(max_jac < -1E-10)
        res[ielem] = max_jac / min_jac;
      else if((min_jac < -1E-10) && (max_jac > +1E-10))
        res[ielem] = -Math::min(-min_jac, max_jac) / Math::max(-min_jac, max_jac);
      else
        res[ielem] = 0.0;

      trafo_eval.finish();
    }

    return res;
  }

  template<typename Shape_>
  void run_reader(Dist::Comm& comm, SimpleArgParser& args, Geometry::MeshFileReader& mesh_reader)
  {
    // define our types
    typedef Geometry::ConformalMesh<Shape_> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange1::Element<TrafoType> SpaceType;
    typedef Geometry::RootMeshNode<MeshType> RootMeshNodeType;

    // create our domain control
    typedef Control::Domain::SimpleDomainLevel<MeshType, TrafoType, SpaceType> DomainLevelType;
    Control::Domain::PartiDomainControl<DomainLevelType> domain(comm, false);

    domain.set_desired_levels(args.query("level")->second);

    domain.create(mesh_reader);

    comm.print("\nDesired Levels: " + domain.format_desired_levels());
    comm.print(  "Chosen  Levels: " + domain.format_chosen_levels());

    const int level_idx = domain.front()->get_level_index();

    double h = 1.0 /  double(1 << level_idx);

    args.parse("h", h);

    comm.print("\nh = " + stringify(h));

    // clone the finest mesh node a few times
    std::unique_ptr<RootMeshNodeType> mesh_node_1 = domain.front()->get_mesh_node()->clone_unique();
    std::unique_ptr<RootMeshNodeType> mesh_node_2 = domain.front()->get_mesh_node()->clone_unique();
    std::unique_ptr<RootMeshNodeType> mesh_node_3 = domain.front()->get_mesh_node()->clone_unique();

    // distort meshes
    Geometry::DistributedMeshDistortion<MeshType> disto_1(comm, *mesh_node_1);
    Geometry::DistributedMeshDistortion<MeshType> disto_2(comm, *mesh_node_2);
    Geometry::DistributedMeshDistortion<MeshType> disto_3(comm, *mesh_node_3);

    disto_1.distort_uniform(0.1 * h);
    disto_2.distort_shortest_edge_uniform(0.1);
    disto_3.distort_shortest_edge_local(0.1);

    Geometry::ExportVTK<MeshType> exporter_1(*mesh_node_1->get_mesh());
    Geometry::ExportVTK<MeshType> exporter_2(*mesh_node_2->get_mesh());
    Geometry::ExportVTK<MeshType> exporter_3(*mesh_node_3->get_mesh());

    // compute jac matrices
    std::vector<double> jacs_1 = analyze_mesh(*mesh_node_1->get_mesh());
    std::vector<double> jacs_2 = analyze_mesh(*mesh_node_2->get_mesh());
    std::vector<double> jacs_3 = analyze_mesh(*mesh_node_3->get_mesh());

    exporter_1.add_cell_scalar("jacs", jacs_1.data());
    exporter_2.add_cell_scalar("jacs", jacs_2.data());
    exporter_3.add_cell_scalar("jacs", jacs_3.data());

    exporter_1.write("dbg-dist-mesh-disto.1", comm);
    exporter_2.write("dbg-dist-mesh-disto.2", comm);
    exporter_3.write("dbg-dist-mesh-disto.3", comm);
  }

  void main(int argc, char* argv [])
  {
    // create world communicator
    Dist::Comm comm(Dist::Comm::world());
    comm.print(String(100u, '*'));

    // dump system call
    {
      String s("Arguments: ");
      s.append(argv[0]);
      for(int i(1); i < argc; ++i)
        s.append(" ").append(argv[i]);
      comm.print(s);
    }

    // print number of processes
    comm.print("Number of Processes: " + stringify(comm.size()));

    // create arg parser
    SimpleArgParser args(argc, argv);

    // check command line arguments
    args.support("mesh");
    args.support("level");
    args.support("h");

    // check for unsupported options
    auto unsupported = args.query_unsupported();
    if (!unsupported.empty())
    {
      // print all unsupported options to cerr
      for (auto it = unsupported.begin(); it != unsupported.end(); ++it)
        comm.print(std::cerr, "ERROR: unknown option '--" + (*it).second + "'");

      comm.print(std::cerr, "Supported Options are:");
      comm.print(std::cerr, args.get_supported_help());

      // abort
      FEAT::Runtime::abort();
    }
    // Our mesh file reader
    Geometry::MeshFileReader mesh_reader;

    // read in the mesh files
    mesh_reader.add_mesh_files(comm, args.query("mesh")->second);

    // read the mesh file root markups
    mesh_reader.read_root_markup();
    String mesh_type = mesh_reader.get_meshtype_string();

    // run the corresponding version
    if(mesh_type == "conformal:hypercube:2:2") run_reader<Shape::Hypercube<2>>(comm, args, mesh_reader); else
    if(mesh_type == "conformal:hypercube:3:3") run_reader<Shape::Hypercube<3>>(comm, args, mesh_reader); else
    if(mesh_type == "conformal:simplex:2:2")   run_reader<Shape::Simplex<2>>  (comm, args, mesh_reader); else
    if(mesh_type == "conformal:simplex:3:3")   run_reader<Shape::Simplex<3>>  (comm, args, mesh_reader); else
    {
      comm.print(std::cerr, "ERROR: unsupported mesh type '" + mesh_type + "'");
      FEAT::Runtime::abort();
    }
  }
} // namespace DbgDistMeshDisto

int main(int argc, char* argv [])
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);
  try
  {
    DbgDistMeshDisto::main(argc, argv);
  }
  catch (const std::exception& exc)
  {
    std::cerr << "ERROR: unhandled exception: " << exc.what() << std::endl;
    FEAT::Runtime::abort();
  }
  catch (...)
  {
    std::cerr << "ERROR: unknown exception" << std::endl;
    FEAT::Runtime::abort();
  }
  return 0;
}
