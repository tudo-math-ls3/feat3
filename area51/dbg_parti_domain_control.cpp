// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/util/dist.hpp>
#include <kernel/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/stop_watch.hpp>
#include <kernel/geometry/export_vtk.hpp>

#include <control/domain/parti_domain_control.hpp>

namespace DbgPartiDomainControl
{
  using namespace FEAT;

  template<typename Shape_>
  void run(Dist::Comm& comm, SimpleArgParser& args, Geometry::MeshFileReader& mesh_reader);

  template<typename DomLvl_>
  void compute_ranks(std::vector<int>& ranks, const typename Control::Domain::PartiDomainControl<DomLvl_>::Ancestor& ancestor);

  void main(int argc, char** argv)
  {
    // create world communicator
    Dist::Comm comm(Dist::Comm::world());

    // print number of processes
    comm.print("Number of Processes: " + stringify(comm.size()));

    // create arg parser
    SimpleArgParser args(argc, argv);
    Control::Domain::add_supported_pdc_args(args);
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


    // Our mesh file reader
    Geometry::MeshFileReader mesh_reader;

    // read in the mesh files
    mesh_reader.add_mesh_files(comm, args.query("mesh")->second);

    // read the mesh file root markups
    mesh_reader.read_root_markup();
    String mesh_type = mesh_reader.get_meshtype_string();

    // run the corresponding version
    if(mesh_type == "conformal:hypercube:2:2") run<Shape::Hypercube<2>>(comm, args, mesh_reader); else
    if(mesh_type == "conformal:hypercube:3:3") run<Shape::Hypercube<3>>(comm, args, mesh_reader); else
    //if(mesh_type == "conformal:simplex:2:2")   run_reader<Shape::Simplex<2>>  (comm, args, mesh_reader); else
    //if(mesh_type == "conformal:simplex:3:3")   run_reader<hape::Simplex<3>>  (comm, args, mesh_reader); else
    {
      comm.print(std::cerr, "ERROR: unsupported mesh type '" + mesh_type + "'");
      FEAT::Runtime::abort();
    }
  }

  template<typename Shape_>
  void run(Dist::Comm& comm, SimpleArgParser& args, Geometry::MeshFileReader& mesh_reader)
  {
    typedef Geometry::ConformalMesh<Shape_> MeshType;
    typedef Control::Domain::DomainLevel<MeshType> DomainLevelType;
    Control::Domain::PartiDomainControl<DomainLevelType> domain(comm, true);

    StopWatch watch_parti;
    domain.parse_args(args);
    domain.set_desired_levels(args.query("level")->second);

    watch_parti.start();
    domain.create(mesh_reader);
    watch_parti.stop();

    comm.print("\nPartitioner runtime: " + watch_parti.elapsed_string());

    comm.print("\nDesired Levels: " + domain.format_desired_levels());
    comm.print(  "Chosen  Levels: " + domain.format_chosen_levels());
    comm.print("\nChosen Partitioning Info:\n" + domain.get_chosen_parti_info());

    // dump domain info if desired
    comm.print("\nDomain Ancestry:");
    comm.allprint(domain.dump_ancestry());
    comm.print("\nDomain Layers:");
    comm.allprint(domain.dump_layers());
    comm.print("\nDomain Layer Levels:");
    comm.allprint(domain.dump_layer_levels());
    comm.print("\nDomain Virtual Levels:");
    comm.allprint(domain.dump_virt_levels());

    const auto& ancestry = domain.get_ancestry();

    if(args.check("vtk") >= 0)
    {
      const MeshType& mesh = domain.front()->get_mesh();
      Geometry::ExportVTK<MeshType> vtk(mesh);

      //for(const auto& ac : ancestry)
      for(std::size_t i(0); i < ancestry.size(); ++i)
      {
        const auto& ac = ancestry.at(i);
        std::vector<int> ranks(mesh.get_num_elements(), ac.progeny_group + ac.progeny_child);
        //compute_ranks(ranks, ac);
        vtk.add_cell_scalar("ac:" + stringify(i), ranks.data());
      }
      vtk.write("dbg-parti-domain-control.n" + stringify(comm.size()), comm);
    }
  }
} // namespace DbgPartiDomainControl

int main(int argc, char** argv)
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);
  DbgPartiDomainControl::main(argc, argv);
  return 0;
}
