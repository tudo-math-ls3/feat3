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

namespace DbgVoxelDomain
{
  using namespace FEAT;

  template<int dim_>
  class RingDropper :
    public Geometry::VoxelMasker<double, dim_>
  {
  public:
    typedef Geometry::VoxelMasker<double, dim_> BaseClass;

    virtual void mask_line(std::vector<int>& mask, const double x_min, const double x_max, const Tiny::Vector<double, dim_>& point) override
    {
      const std::size_t nv = mask.size();
      Tiny::Vector<double, dim_> midpoint(0.0);
      Tiny::Vector<double, dim_> pt(point);

      for(std::size_t i(0); i < nv; ++i)
      {
        pt[0] = BaseClass::x_coord(x_min, x_max, i, nv);
        double r = (pt - midpoint).template norm_euclid_n<2>();
        mask[i] = ((r > 1.1) && (r < 2.1) ? 1 : 0);
      }
    }
  };

  //typedef Shape::Quadrilateral ShapeType;
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
  using namespace DbgVoxelDomain;

  Runtime::ScopeGuard runtime_guard(argc, argv);

  Dist::Comm comm = Dist::Comm::world();

  SimpleArgParser args(argc, argv);
  args.support("level");
  args.support("debug");

  XASSERT(args.query_unsupported().empty());

  Control::Domain::VoxelDomainControl<DomainLevelType> domain(comm, true);

  // create base mesh node
  //domain.create_base_mesh_2d(4, 4, 0.0, 1.0, 0.0, 1.0);
  //domain.create_base_mesh_2d(6, 6, -3.0*15./16., 3.0*15./16., -3.0*15./16., 3.0*15./16.);
  domain.create_base_mesh_3d(6, 6, 3, -3.0*13./16., 3.0*13./16., -3.0*13./16., 3.0*13./16., 0., 3.);
  //domain.create_base_mesh_2d(6, 6, -3.0*13./16., 3.0*13./16., -3.0*13./16., 3.0*13./16.);

  // todo: add base mesh meshparts

  domain.parse_args(args);
  //domain.set_desired_levels(args.query("level")->second);
  domain.set_desired_levels("4 0");
  //domain.set_desired_levels("3 1:2 0");
  //domain.set_desired_levels("4 2:4 0");
  //domain.set_desired_levels("6 4:2 2:1 0");
  //domain.set_desired_levels("4 3:4 2:2 1:1 0");
  //domain.set_desired_levels(5, 2, 0);

  domain.keep_voxel_map();
  domain.create_voxel_map(*std::make_unique<RingDropper<ShapeType::dimension>>(), 0.1);
  //domain.create_voxel_map_from_lambda([](auto p) { return (p[1] < 2.375-2.25*p[0] ) && (p[1] < 2.25*p[0]+0.125); }, Real(0));
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
  }

  const Index num_levels = domain.size_physical();

  std::deque<std::shared_ptr<SystemLevelType>> system(num_levels);
  std::deque<typename SystemLevelType::GlobalSystemVector> vectors(num_levels);

  const String cubature("auto-degree:5");

  // create system levels
  for (Index i(0); i < num_levels; ++i)
  {
    domain.at(i)->domain_asm.compile_all_elements();
    system.at(i) = std::make_shared<SystemLevelType>();
    system.at(i)->assemble_gate(domain.at(i));

    {
      String m("Neighbors Level ");
      m += stringify(i) + ": ";
      String n;
      String v(std::size_t(comm.size()), ' ');
      v[std::size_t(domain.at(i).layer().comm().rank())] = ':';
      for(int r : system.at(i)->gate_sys._ranks)
      {
        (n += " ") += stringify(r);
        v[std::size_t(r)] = 'X';
      }
      domain.at(i).layer().comm().allprint(m + v + n);
      //comm.print(m + v);
    }

    /*if((i+1) < domain.size_virtual())
    {
      system.at(i)->assemble_coarse_muxer(domain.at(i+1));
      //system.at(i)->assemble_transfer(domain.at(i), domain.at(i+1), cubature);
      //assemble_sad_transfer(*system.at(i), domain.at(i), domain.at(i+1), cubature);
      Control::VoxelTransferAssembler::assemble_scalar_basic_transfer(*system.at(i), domain.at(i), domain.at(i+1), cubature);
    }
    vectors.at(i) = typename SystemLevelType::GlobalSystemVector(&system.at(i)->gate_sys,
      system.at(i)->gate_sys._freqs.clone(LAFEM::CloneMode::Layout));
    */
  }

  /*
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
  }
  //*/
  for(Index k(0); k < domain.size_physical(); ++k)
  {
    Geometry::ExportVTK<MeshType> vtk(domain.at(k)->get_mesh());
    std::vector<int> mask = domain.gather_vertex_voxel_map(k); //gather_vertex_slag_mask(k);
    std::vector<Real> weights = domain.gather_element_voxel_weights(k); //gather_element_slag_weights(k);
    vtk.add_vertex_scalar("mask", mask.data());
    //vtk.add_vertex_scalar("test", vectors.at(k).local().elements());
    vtk.add_cell_scalar("weights", weights.data());
    vtk.add_cell_scalar("color", domain.at(k)->element_coloring.get_coloring());
    if(!domain.at(k)->element_layering.empty())
      vtk.add_cell_scalar("layer", domain.at(k)->element_layering.get_coloring());
    /*Adjacency::Graph& layering = domain.at(k)->element_layering;
    if(layering.get_num_nodes_domain() > 0u)
    {
      std::vector<Index> vlyr(domain.at(k)->element_coloring.get_num_nodes(), 0u);
      const Index* dom_ptr = layering.get_domain_ptr();
      const Index* img_idx = layering.get_image_idx();
      for(Index i(0); i < layering.get_num_nodes_domain(); ++i)
        for(Index j(dom_ptr[i]); j < dom_ptr[i+1]; ++j)
          vlyr[img_idx[j]] = i;
      vtk.add_cell_scalar("layer", vlyr.data());
    }*/
    vtk.write("dbg-voxel-domain.lvl" + stringify(domain.at(k)->get_level_index()), domain.at(k).layer().comm());
  }
  //for(Index k(0); (k < domain.size_physical() + (domain.has_ghost() ? 1 : 0)); ++k)
  /*{
    int k = 2;
    Geometry::ExportVTK<MeshType> vtk_c(domain.at(k).level_c().get_mesh());
    vtk_c.write("dbg-voxel-domain-c.lvl" + stringify(domain.at(k)->get_level_index()), domain.at(k).layer_c().comm());
  }*/

  return 0;
}
