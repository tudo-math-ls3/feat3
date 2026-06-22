// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include "kernel/geometry/export_vtk.hpp"
#include "kernel/geometry/mesh_file_writer.hpp"
#include "kernel/util/simple_arg_parser.hpp"
#include <filesystem>
#include <kernel/base_header.hpp>
#include <kernel/geometry/cgal.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/intern/coarse_fine_cell_mapping.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/meshopt/mesh_deformer.hpp>
#include <kernel/runtime.hpp>
#include <kernel/shape.hpp>
#include <kernel/solver/gmres.hpp>
#include <kernel/trafo/inverse_mapping.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/arg_parser.hpp>
#include <kernel/util/dist.hpp>

namespace MeshDeformer
{
  using namespace FEAT;

  using MeshType = Geometry::ConformalMesh<Shape::Hexahedron>;
  using DataType = MeshType::CoordType;
  using MeshNodeType = Geometry::RootMeshNode<MeshType>;
  using TrafoType = Trafo::Standard::Mapping<MeshType>;
  using InverseMappingType = Trafo::InverseMapping<TrafoType, DataType>;
  using HInverseMappingType =
    Trafo::HierarchicalInverseMapping<InverseMappingType, Geometry::Intern::CoarseFineCellMapping<MeshType>>;

  static constexpr int dim = 3;

  struct MeshDef
  {
    Tiny::Vector<Index, dim> num_cells;
    Tiny::Vector<Real, dim> min;
    Tiny::Vector<Real, dim> max;

    friend std::ostream& operator<<(std::ostream& s, const MeshDef& d)
    {
      s << "MeshDef{ num_cells: " << d.num_cells << ", min: " << d.min << ", max: " << d.max << " }";
      return s;
    }

    friend std::istream& operator>>(std::istream& s, MeshDef& d)
    {
      for(int i(0); i < dim; i++)
      {
        s >> d.min(i);
      }
      for(int i(0); i < dim; i++)
      {
        s >> d.max(i);
      }
      for(int i(0); i < dim; i++)
      {
        s >> d.num_cells(i);
      }

      return s;
    }
  };

  struct MeshDeformerArgs : public ArgParser
  {
    Parameter<bool> help = parameter(false).short_flag("-h").long_flag("--help").help_string("Show this help text");
    Parameter<String> output_filename =
      parameter(String("mesh-lvl-"))
        .long_flag("--output-filename")
        .property("output-filename")
        .help_string("Prefix for output files. Output files will be written as prefix<level>.xml");
    Parameter<String> surface_mesh =
      parameter(String(""))
        .long_flag("--surface-mesh")
        .property("surface-mesh")
        .required()
        .help_string(
          "The surface mesh that is used to determine target cell sizes.\nCells close to the surface mesh will become "
          "smaller.\nMust be a triangle mesh in the OFF format.");
    Parameter<MeshDef> start_mesh =
      parameter(MeshDef{})
        .long_flag("--start-mesh-definition")
        .property("start-mesh-definition")
        .required()
        .help_string(
          "Definition for the starting mesh. Formatted as <x_min y_min z_min x_max y_max z_max cells_x cells_y "
          "cells_z>.\nThe starting mesh is a cuboid with the given size, divided into the given number of cells along "
          "each axis.\nWill be refined further until it has at least --min-cells cells before deformation "
          "starts.\nNote "
          "that the definition has to be given as a single parameter, i.e. it has to be quoted if passed via a command "
          "line argument.");
    Parameter<Index> min_cells = parameter(Index(10000))
                                   .long_flag("--min-cells")
                                   .property("min-cells")
                                   .help_string(
                                     "Minimum mesh size before deformation begins. Input meshes will be globally "
                                     "refined until this target number is hit.");
    Parameter<DataType> rel_min =
      parameter(DataType(0.2))
        .long_flag("--rel-min")
        .property("rel-min")
        .help_string(
          "Minimum relative taret size for cells. This is the relative size for cells near the surface mesh.");
    Parameter<DataType> rel_max_inner =
      parameter(DataType(0.5))
        .long_flag("--rel-max-inner")
        .property("rel-max-inner")
        .help_string("Maximum relative target size for cells within the surface mesh.");
    Parameter<DataType> rel_max_outer =
      parameter(DataType(1.0))
        .long_flag("--rel-max-outer")
        .property("rel-max-outer")
        .help_string("Maximum relative target size for cells outside the surface mesh.");
    Parameter<DataType> radius = parameter(DataType(1.0))
                                   .long_flag("--radius")
                                   .property("radius")
                                   .help_string(
                                     "Interpolation radius for target cell sizes. Cell sizes are interpolated from "
                                     "rel-min (on the surface) to rel-max-inner or rel-max-outer within this radius.\n"
                                     "All cells beyond this radius get a uniform target cell size.");

    Parameter<String> config_file =
      parameter(String(""))
        .long_flag("--config-file")
        .help_string("Config file in FEAT's PropertyMap format. See help text for property names of all options.");

    Parameter<Index> max_iters = parameter(Index(3))
                                   .long_flag("--max-iters")
                                   .property("max-iters")
                                   .help_string(
                                     "Maximum number of iterations in each blend step. Algorithm will exit early if "
                                     "quality of deformation is good enough.");
    Parameter<Index> blend_steps =
      parameter(Index(1))
        .long_flag("--blend-steps")
        .property("blend-steps")
        .help_string(
          "Maximum number of blend steps. Will cause target cell sizes to slowly increase from current sizes to actual "
          "target sizes. More blend steps cause the deformation to be more stable.");
    Parameter<Index> post_refine =
      parameter(Index(0))
        .long_flag("--post-refine")
        .property("post-refine")
        .help_string("Number of times to refine the deformed mesh further after deformation.");

    Parameter<bool> vtk = parameter(false).long_flag("--vtk").property("vtk").help_string("If true, write vtk of deformed mesh.");

    Parameter<bool> export_undeformed = parameter(false).long_flag("--export-undeformed").property("export-undeformed").help_string("Export mesh before deforming it");
  };

  static int main(int argc, char** argv)
  {
    SimpleArgParser bootstrap(argc, argv);

    MeshDeformerArgs args;
    args.set_description(
      "Mesh deformation program. Deforms a mesh to match target cell sizes. Cell sizes are derived from a surface "
      "mesh. Cells closer to the surface will become smaller.");

    if(bootstrap.check("help") >= 0)
    {
      std::cout << args.help_text() << "\n";
      return 0;
    }

    bool parse_success(false);
    if(bootstrap.check("config-file") > 0)
    {
      String config_file;
      bootstrap.parse("config-file", config_file);
      parse_success = args.parse(argc, argv, config_file);
    }
    else
    {
      parse_success = args.parse(argc, argv);
    }

    if(!parse_success)
    {
      for(const String& s : args.errors())
      {
        std::cerr << s << "\n";
      }

      return 1;
    }

    std::cout << "Running with settings:\n";
    std::cout << args.display() << "\n";

    // Create starting mesh
    std::unique_ptr<MeshNodeType> mesh_node = MeshNodeType::make_unique(
      Geometry::StructUnitCubeFactory<MeshType>::make_unique_from(
        args.start_mesh->num_cells[0],
        args.start_mesh->num_cells[1],
        args.start_mesh->num_cells[2]));

    // Scale and translate to correct bounding box
    Tiny::Vector<Real, dim> translation = args.start_mesh->min;
    Tiny::Vector<Real, dim> scale = args.start_mesh->max - args.start_mesh->min;
    {
      auto& vtx_set = mesh_node->get_mesh()->get_vertex_set();
      const Index n = vtx_set.get_num_vertices();
      for(Index i(0); i < n; ++i)
      {
        for(int k = 0; k < dim; ++k)
        {
          vtx_set[i][k] = translation[k] + vtx_set[i][k] * scale[k];
        }
      }
    }


    while(mesh_node->get_mesh()->get_num_elements() < *args.min_cells)
    {
      mesh_node = mesh_node->refine_unique();
    }

    if(*args.export_undeformed)
    {
      Geometry::ExportVTK<MeshType> exporter(*mesh_node->get_mesh());
      exporter.write("undeformed");
    }

    // Load surface mesh
    Geometry::CGALWrapper<DataType> cgal(*args.surface_mesh, Geometry::CGALFileMode::fm_off);

    // Define monitor function
    const Real radius = *args.radius;
    const auto monitor_fn = [&](auto vertex)
    {
      std::cout << "monitor_fn: vertex = " << vertex << "\n";
      Real dist = Math::sqrt(cgal.squared_distance(vertex[0], vertex[1], vertex[2]));
      bool inside = cgal.point_inside(vertex[0], vertex[1], vertex[2]);
      const Real rel_min = *args.rel_min;
      const Real rel_max = inside ? *args.rel_max_inner : *args.rel_max_outer;

      // Interpolate between rel_min and rel_max within radius around the surface
      return Math::clamp(rel_min + ((dist / radius) * (rel_max - rel_min)), rel_min, rel_max);
    };

    // Create inverse mapping
    TrafoType trafo(*mesh_node->get_mesh());
    InverseMappingType inv_mapping(trafo);

    // Create trivial gate
    Dist::Comm comm(Dist::Comm::world());
    Global::Gate<LAFEM::DenseVector<DataType, Index>, LAFEM::VectorMirror<DataType, Index>> gate(comm);

    // Precondition for deformer
    mesh_node->get_mesh()->reorient_boundary_facets();

    MeshOpt::MeshDeformer<MeshType> deformer(&gate);
    deformer.set_blend_steps(*args.blend_steps);
    deformer.set_max_iters_per_blend_step(*args.max_iters);
    auto quality = deformer.deform(*mesh_node, monitor_fn);

    std::cout << "Quality (0-Norm)  : " << quality.l2norm << "\n";
    std::cout << "Quality (inf-Norm): " << quality.maxnorm << "\n";

    for(Index i(0); i <= *args.post_refine; i++)
    {
      std::ofstream out(*args.output_filename + stringify(i) + ".xml");
      Geometry::MeshFileWriter writer(out);
      writer.write(mesh_node.get());

      if(*args.vtk)
      {
        Geometry::ExportVTK<MeshType> exporter(*mesh_node->get_mesh());
        exporter.write(*args.output_filename + stringify(i));
      }

      if(i < *args.post_refine)
      {
        mesh_node = mesh_node->refine_unique();
      }
    }

    return 0;
  }
} // namespace MeshDeformer

int main(int argc, char** argv)
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  return MeshDeformer::main(argc, argv);
}
