// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/analytic/lambda_function.hpp>

#include <control/domain/voxel_domain_control.hpp>
#include <control/scalar_basic.hpp>

namespace DbgAdaptDrop
{
  using namespace FEAT;

  typedef Shape::Hexahedron ShapeType;

  typedef Geometry::ConformalMesh<ShapeType> MeshType;

  typedef Geometry::RootMeshNode<MeshType> MeshNodeType;

  typedef Trafo::Standard::Mapping<MeshType> TrafoType;

  typedef Space::Lagrange1::Element<TrafoType> SpaceType;

  typedef Control::Domain::SimpleDomainLevel<MeshType, TrafoType, SpaceType> DomainLevelSubType;
  typedef Control::Domain::VoxelDomainLevelWrapper<DomainLevelSubType> DomainLevelType;

  typedef Control::Domain::VoxelDomainControl<DomainLevelType> DomainControl;

  typedef Control::ScalarBasicSystemLevel<> SystemLevelType;
}


int main(int argc, char** argv)
{
  using namespace FEAT;
  using namespace DbgAdaptDrop;

  Runtime::ScopeGuard runtime_guard(argc, argv);

  Dist::Comm comm = Dist::Comm::world();

  Control::Domain::VoxelDomainControl<DomainLevelType> domain(comm, true);

  /*if(argc < 2)
  {
    std::cout << "USAGE: " << argv[0] << " <OFF-file>\n";
    return 0;
  }*/

  SimpleArgParser args(argc, argv);

  // create base mesh node
  //domain.create_base_mesh_2d(4, 4, 0.0, 1.0, 0.0, 1.0);
  //domain.create_base_mesh_3d(4, 4, 4, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
  //domain.create_base_mesh_3d(6, 7, 8, -96.0, 96.0, -164.0, 60.0, 0.0, 256.0);
  domain.create_base_mesh_3d(6, 7, 9, -96.0, 96.0, -164.0, 60.0, 0.0, 288.0);
  domain.parse_args(args);

  // todo: add base mesh meshparts

  //domain.set_desired_levels("5 0");
  //domain.set_desired_levels("3 2:1 0");
  //domain.set_desired_levels("6 4:2 2:1 0");
  //domain.set_desired_levels("4 3:4 2:2 1:1 0");
  //domain.set_desired_levels(5, 2, 0);
  domain.set_desired_levels(args.query("level")->second);

  // create slag mask from OFF file
  if(args.check("off") > 0)
  {
    String off_name = args.query("off")->second.front();
    comm.print("Creating voxel map from OFF file '" + off_name + "'...");
    domain.create_voxel_map_from_off(off_name, false, 0.0);
  }
  if(args.check("vxl") > 0)
  {
    String vxl_name = args.query("vxl")->second.front();
    comm.print("Loading voxel map from file '" + vxl_name + "'...");
    domain.read_voxel_map(vxl_name);
  }

  domain.keep_voxel_map();
  domain.create_hierarchy();

  comm.print("Desired Levels: " + domain.format_desired_levels());
  comm.print("Chosen  Levels: " + domain.format_chosen_levels());

  comm.print("Base-Mesh Creation Time: " + domain.get_watch_base_mesh().elapsed_string().pad_front(7));
  comm.print("Voxel-Map Creation Time: " + domain.get_watch_voxel_map().elapsed_string().pad_front(7));
  comm.print("Hierarchy Creation Time: " + domain.get_watch_hierarchy().elapsed_string().pad_front(7));


  if(1)
  {
    comm.print("\nAncestry Dump");
    comm.allprint(domain.dump_ancestry());
    comm.print("\nLayer Dump:");
    comm.allprint(domain.dump_layers());
    //comm.print("\nLayer Level Dump:");
    //comm.allprint(domain.dump_layer_levels());
    comm.print("\nVirtual Level Dump:");
    comm.allprint(domain.dump_virt_levels());
    comm.print("\nSlag Layer Level Dump:");
    comm.allprint(domain.dump_slag_layer_levels());

    String str;
    for(std::size_t i(0); i < domain.size_virtual(); ++i)
    {
      Index nel = 0;
      if(i < domain.size_physical())
        nel = domain.at(i)->get_mesh().get_num_elements();
      comm.allreduce(&nel, &nel, 1, Dist::op_sum);
      if(comm.rank() == 0)
        str = stringify(domain.at(i)->get_level_index()) + ": " + stringify(nel).pad_front(10) + "\n" + str;
    }
    comm.print(String("Element Counts:\n") + str);
  }

  /*if(comm.rank() == 0)
  {
    std::ofstream ofs("slagmask.bin", std::ios_base::binary);
    const std::vector<char>& slagvec = domain.get_voxel_map().get_map(); //domain.get_slag_mask_vector();
    ofs.write(slagvec.data(), std::streamsize(slagvec.size() * 8ull));
  }*/

  const Index num_levels = domain.size_physical();

  std::deque<std::shared_ptr<SystemLevelType>> system(num_levels);
  std::deque<typename SystemLevelType::GlobalSystemVector> vectors(num_levels);

  const String cubature("auto-degree:5");

  // create system levels
  /*for (Index i(0); i < num_levels; ++i)
  {
    domain.at(i)->domain_asm.compile_all_elements();
    system.at(i) = std::make_shared<SystemLevelType>();
    system.at(i)->assemble_gate(domain.at(i));
    if((i+1) < domain.size_virtual())
    {
      system.at(i)->assemble_coarse_muxer(domain.at(i+1));
      //system.at(i)->assemble_transfer(domain.at(i), domain.at(i+1), cubature);
      //assemble_sad_transfer(*system.at(i), domain.at(i), domain.at(i+1), cubature);
      Control::SADTransferAssembler::assemble_scalar_basic_transfer(*system.at(i), domain.at(i), domain.at(i+1), cubature);
    }
    vectors.at(i) = typename SystemLevelType::GlobalSystemVector(&system.at(i)->gate_sys,
      system.at(i)->gate_sys._freqs.clone(LAFEM::CloneMode::Layout));
  }

  auto func_xy = Analytic::create_lambda_function_scalar_2d([](double x, double y) {return x+y;});

  for(Index i(0); i < num_levels; ++i)
  {
    Index ii = num_levels - i - 1u;

    if(ii+1u < domain.size_physical())
    {
      // prolongate
      system.at(ii)->transfer_sys.prol(vectors.at(ii), vectors.at(ii+1u));
    }
    else if(ii+1u < domain.size_virtual())
    {
      system.at(ii)->transfer_sys.prol_recv(vectors.at(ii));
    }
    else if(ii+1u == domain.size_virtual())
    {
      Assembly::Interpolator::project(vectors.at(ii).local(), func_xy, domain.at(ii)->space);
    }
  }*/

  for(Index k(0); k < domain.size_physical(); ++k)
  {
    Geometry::ExportVTK<MeshType> vtk(domain.at(k)->get_mesh());
    std::vector<int> mask = domain.gather_vertex_voxel_map(k);
    std::vector<Real> weights = domain.gather_element_voxel_weights(k);
    vtk.add_vertex_scalar("mask", mask.data());
    //vtk.add_vertex_scalar("test", vectors.at(k).local().elements());
    vtk.add_cell_scalar("color", domain.at(k)->element_coloring.get_coloring());
    vtk.add_cell_scalar("layer", domain.at(k)->element_layering.get_coloring());
    vtk.add_cell_scalar("weights", weights.data());
    vtk.write("dbg-voxel-domain-gendie.lvl" + stringify(domain.at(k)->get_level_index()), domain.at(k).layer().comm());
  }

  return 0;
}
