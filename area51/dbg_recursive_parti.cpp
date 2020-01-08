// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>

#include <control/domain/parti_domain_control.hpp>

namespace DbgRecursiveParti
{
  using namespace FEAT;

  template<typename MeshType_, typename DomainLevel_>
  void write_vtk(const DomainLevel_& dom_lvl, String name)
  {
    const Dist::Comm& comm = dom_lvl.layer().comm();

    String vname = name + "." + stringify(dom_lvl->get_level_index());
    Geometry::ExportVTK<MeshType_> vtk(dom_lvl->get_mesh());

    // add all halos
    std::vector<double> vtx(dom_lvl->get_mesh().get_num_vertices(), 0.0);

    for(int r(0); r < comm.size(); ++r)
    {
      const auto* halo = dom_lvl->get_mesh_node_ptr()->get_halo(r);
      if(halo == nullptr)
      {
        vtk.add_vertex_scalar("halo:" + stringify(r).pad_front(2,'0'), vtx.data());
        continue;
      }
      Index n = halo->get_num_entities(0);
      const auto& v = halo->template get_target_set<0>();
      for(Index i(0); i < n; ++i) vtx[v[i]] = 1.0;
      vtk.add_vertex_scalar("halo:" + stringify(r).pad_front(2,'0'), vtx.data());
      for(Index i(0); i < n; ++i) vtx[v[i]] = 0.0;
    }

    vtk.write(vname, comm);
  }

  template<typename Shape_>
  void run(SimpleArgParser& args, Dist::Comm& comm, Geometry::MeshFileReader& mesh_reader)
  {
    typedef Shape_ ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange1::Element<TrafoType> SpaceType;

    // create our domain control
    typedef Control::Domain::SimpleDomainLevel<MeshType, TrafoType, SpaceType> DomainLevelType;
    Control::Domain::PartiDomainControl<DomainLevelType> domain(comm, true);

    domain.parse_args(args);
    domain.set_desired_levels(args.query("level")->second);

    domain.create(mesh_reader);

    // plot our levels
    comm.print("\nDesired Levels: " + domain.format_desired_levels());

    comm.print("\nChosen Partitioning Info:\n" + domain.get_chosen_parti_info());

    comm.print("\nChosen Layers:");
    comm.allprint(domain.dump_layers());

    comm.print("\nAncestry Info:");
    comm.allprint(domain.dump_ancestry());

    comm.print("\nLayer Levels:");
    comm.allprint(domain.dump_layer_levels());

    comm.print("\nVirtual Levels:");
    comm.allprint(domain.dump_virt_levels());

    comm.print("");
    comm.print("Desired Levels: " + domain.format_desired_levels());
    comm.print("Chosen  Levels: " + domain.format_chosen_levels());

    // write VTKs
    if(args.check("vtk") >= 0)
    {
      String name("dbg-recurse-parti");
      args.parse("vtk", name);
      comm.print("Writing VTK files '" + name + ".*.vtu'...");
      for(Index lvl(0); lvl < domain.size_physical(); ++lvl)
      {
        write_vtk<MeshType>(domain.at(lvl), name);
      }
      /*{
        String vname = name + "." + stringify(lvl);
        Geometry::ExportVTK<MeshType> vtk(domain.at(lvl)->get_mesh());
        vtk.write(vname, domain.at(lvl).layer().comm());
      }*/
    }

  }

  void main(int argc, char** argv)
  {
    // create world communicator
    Dist::Comm comm(Dist::Comm::world());

    // print number of processes
    comm.print("Number of Processes: " + stringify(comm.size()));

    // create arg parser
    SimpleArgParser args(argc, argv);

    // check command line arguments
    Control::Domain::add_supported_pdc_args(args);
    args.support("debug");
    args.support("mesh");
    args.support("level");
    args.support("vtk");

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

    if(args.check("mesh") < 1)
    {
      comm.print(std::cerr, "ERROR: Mandatory option '--mesh <mesh-file>' is missing!");
      FEAT::Runtime::abort();
    }
    if(args.check("level") < 1)
    {
      comm.print(std::cerr, "ERROR: Mandatory option '--level <levels>' is missing!");
      FEAT::Runtime::abort();
    }

#ifdef FEAT_COMPILER_MICROSOFT
    if(args.check("debug") > 0)
    {
      const auto& va = args.query("debug")->second;
      for(const auto& s : va)
      {
        int debug_rank = -1;
        if(s.parse(debug_rank) && (comm.rank() == debug_rank))
        {
          __debugbreak();
          break;
        }
      }
    }
#endif

    // Our mesh file reader
    Geometry::MeshFileReader mesh_reader;

    // read in the mesh files
    mesh_reader.add_mesh_files(comm, args.query("mesh")->second);

    // read the mesh file root markups
    mesh_reader.read_root_markup();
    String mesh_type = mesh_reader.get_meshtype_string();

    // run 2D or 3D ?
    if(mesh_type == "conformal:hypercube:2:2") run<Shape::Hypercube<2>>(args, comm, mesh_reader); else
    //if(mesh_type == "conformal:hypercube:3:3") run<Shape::Hypercube<3>>(args, comm, mesh_reader); else
    //if(mesh_type == "conformal:simplex:2:2")   run<Shape::Simplex<2>>  (args, comm, mesh_reader); else
    //if(mesh_type == "conformal:simplex:3:3")   run<Shape::Simplex<3>>  (args, comm, mesh_reader); else
    {
      comm.print(std::cerr, "ERROR: unsupported mesh type '" + mesh_type + "'");
      FEAT::Runtime::abort();
    }

    comm.print("\nFinished");
  }
} // namespace DbgRecursiveParti

int main(int argc, char** argv)
{
  FEAT::Runtime::initialise(argc, argv);
#ifdef DEBUG
  DbgRecursiveParti::main(argc, argv);
#else
  try
  {
    DbgRecursiveParti::main(argc, argv);
  }
  catch(std::exception& exc)
  {
    std::cout << exc.what() << std::endl;
    FEAT::Runtime::abort();
  }
#endif
  return FEAT::Runtime::finalise();
}
