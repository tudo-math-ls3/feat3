// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include "kernel/geometry/common_factories.hpp"
#include "kernel/geometry/templates/sun_zhao_ma_expansion_data.hpp"
#include "kernel/shape.hpp"
#include <kernel/geometry/adaptive_mesh.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/subdivision_levels.hpp>
#include <kernel/util/simple_arg_parser.hpp>

using namespace FEAT;

namespace DbgAdaptiveMesh
{
  /**
   * Produces a vtk file containing all 3D templates of the chosen template set
   */
  template<typename TemplateSet_>
  void export_templates_3d()
  {
    using FoundationMeshType = Geometry::ConformalMesh<Shape::Hexahedron>;

    std::uint64_t num_templates = 256;

    const std::array<Index, 4> mesh_size = {{8 * num_templates, 12 * num_templates, 6 * num_templates, num_templates}};
    FoundationMeshType foundation_mesh(mesh_size.data());

    Geometry::SubdivisionLevels levels(foundation_mesh.get_num_vertices());

    auto& vertex_set = foundation_mesh.get_vertex_set();
    auto& v_at_c = foundation_mesh.get_index_set<3, 0>();
    auto& v_at_f = foundation_mesh.get_index_set<2, 0>();
    auto& v_at_e = foundation_mesh.get_index_set<1, 0>();
    auto& e_at_c = foundation_mesh.get_index_set<3, 1>();
    auto& e_at_f = foundation_mesh.get_index_set<2, 1>();
    auto& f_at_c = foundation_mesh.get_index_set<3, 2>();

    using CellFaceMap = Geometry::Intern::FaceIndexMapping<Shape::Hexahedron, 2, 0>;
    using CellEdgeMap = Geometry::Intern::FaceIndexMapping<Shape::Hexahedron, 1, 0>;
    using FaceEdgeMap = Geometry::Intern::FaceIndexMapping<Shape::Hexahedron, 2, 1>;

    for(std::uint64_t type = 0; type < num_templates; type++)
    {
      Index indices[8];
      for(Index i = 0; i < 8; i++)
      {
        indices[i] = 8 * type + i;
      }

      vertex_set[indices[0]] = {2.0 * (Real)type, 0.0, 0.0};
      vertex_set[indices[1]] = {2.0 * (Real)type + 1.0, 0.0, 0.0};
      vertex_set[indices[2]] = {2.0 * (Real)type, 1.0, 0.0};
      vertex_set[indices[3]] = {2.0 * (Real)type + 1.0, 1.0, 0.0};

      vertex_set[indices[4]] = {2.0 * (Real)type, 0.0, 1.0};
      vertex_set[indices[5]] = {2.0 * (Real)type + 1.0, 0.0, 1.0};
      vertex_set[indices[6]] = {2.0 * (Real)type, 1.0, 1.0};
      vertex_set[indices[7]] = {2.0 * (Real)type + 1.0, 1.0, 1.0};

      for(int i = 0; i < 12; i++)
      {
        e_at_c[type][i] = 12 * type + i;

        v_at_e[12 * type + i][0] = indices[CellEdgeMap::map(i, 0)];
        v_at_e[12 * type + i][1] = indices[CellEdgeMap::map(i, 1)];
      }

      for(int i = 0; i < 6; i++)
      {
        f_at_c[type][i] = 6 * type + i;

        v_at_f[6 * type + i][0] = indices[CellFaceMap::map(i, 0)];
        v_at_f[6 * type + i][1] = indices[CellFaceMap::map(i, 1)];
        v_at_f[6 * type + i][2] = indices[CellFaceMap::map(i, 2)];
        v_at_f[6 * type + i][3] = indices[CellFaceMap::map(i, 3)];

        e_at_f[6 * type + i][0] = 12 * type + FaceEdgeMap::map(i, 0);
        e_at_f[6 * type + i][1] = 12 * type + FaceEdgeMap::map(i, 1);
        e_at_f[6 * type + i][2] = 12 * type + FaceEdgeMap::map(i, 2);
        e_at_f[6 * type + i][3] = 12 * type + FaceEdgeMap::map(i, 3);
      }

      for(int i = 0; i < 8; i++)
      {
        v_at_c[type][i] = indices[i];
        levels[indices[i]] = (type & (1UL << (std::uint64_t) i)) > 0 ? 1 : 0;
      }
    }

    Geometry::AdaptiveMesh<TemplateSet_, Shape::Hexahedron> adaptive_mesh(foundation_mesh);

    std::cout << "Building Adaptive Mesh\n";
    adaptive_mesh.adapt(levels);

    std::cout << "Building ConformalMesh\n";
    FoundationMeshType refined_mesh = adaptive_mesh.to_conformal_mesh(Geometry::Layer{1});

    std::cout << "Exporting Mesh\n";
    Geometry::ExportVTK<FoundationMeshType> exporter(refined_mesh);

    exporter.write("3d_templates");
    std::cout << "Done\n";
  }

  /**
   * Produces a vtk file containing all 2D templates of the chosen template set
   */
  template<typename TemplateSet_>
  void export_templates_2d()
  {
    using FoundationMeshType = Geometry::ConformalMesh<Shape::Quadrilateral>;

    const std::array<Index, 4> mesh_size = {{4 * 16, 4 * 16, 16, 0}};
    FoundationMeshType foundation_mesh(mesh_size.data());

    Geometry::SubdivisionLevels levels(foundation_mesh.get_num_vertices());

    auto& vertex_set = foundation_mesh.get_vertex_set();
    auto& v_at_f = foundation_mesh.get_index_set<2, 0>();
    auto& v_at_e = foundation_mesh.get_index_set<1, 0>();
    auto& e_at_f = foundation_mesh.get_index_set<2, 1>();

    using EdgeMap = Geometry::Intern::FaceIndexMapping<Shape::Quadrilateral, 1, 0>;

    for(std::uint64_t type = 0; type < 16; type++)
    {
      Index indices[4] = {4 * type, 4 * type + 1, 4 * type + 2, 4 * type + 3};
      vertex_set[indices[0]] = {2.0 * (Real)type, 0.0};
      vertex_set[indices[1]] = {2.0 * (Real)type + 1.0, 0.0};
      vertex_set[indices[2]] = {2.0 * (Real)type, 1.0};
      vertex_set[indices[3]] = {2.0 * (Real)type + 1.0, 1.0};

      for(int i = 0; i < 4; i++)
      {
        v_at_f[type][i] = indices[i];
        v_at_e[4 * type + i][0] = indices[EdgeMap::map(i, 0)];
        v_at_e[4 * type + i][1] = indices[EdgeMap::map(i, 1)];
        e_at_f[type][i] = 4 * type + i;
        levels[indices[i]] = (type & (1UL << (std:: uint64_t)i)) > 0 ? 1 : 0;
      }
    }

    Geometry::AdaptiveMesh<TemplateSet_, Shape::Quadrilateral> adaptive_mesh(foundation_mesh);

    std::cout << "Building Adaptive Mesh\n";
    adaptive_mesh.adapt(levels);

    std::cout << "Building ConformalMesh\n";
    FoundationMeshType refined_mesh = adaptive_mesh.to_conformal_mesh(Geometry::Layer{1});

    std::cout << "Exporting Mesh\n";
    Geometry::ExportVTK<FoundationMeshType> exporter(refined_mesh);

    exporter.write("2d_templates");
    std::cout << "Done\n";
  }

  template<typename TemplateSet_>
  void concave_refinement()
  {
    using FoundationMeshType = Geometry::ConformalMesh<Shape::Hexahedron>;

    Geometry::RefinedUnitCubeFactory<FoundationMeshType> foundation_factory(4);
    FoundationMeshType foundation(foundation_factory);

    Geometry::SubdivisionLevels sdls(foundation.get_num_vertices());

    auto& vertices = foundation.get_vertex_set();
    for(Index i(0); i < vertices.get_num_vertices(); ++i)
    {
      auto& vertex = vertices[i];
      if(vertex[0] > 0.5 || vertex[1] > 0.5 || vertex[2] > 0.5)
      {
        sdls[i] = 1;
      }
      else
      {
        sdls[i] = 0;
      }
    }

    Geometry::AdaptiveMesh<TemplateSet_, Shape::Hexahedron> a_mesh(foundation);
    a_mesh.adapt(sdls, Geometry::AdaptiveMesh<TemplateSet_, Shape::Hexahedron>::ImportBehaviour::All);

    std::cout << "Building Adaptive Mesh\n";
    a_mesh.adapt(sdls);

    std::cout << "Building ConformalMesh\n";
    FoundationMeshType refined_mesh = a_mesh.to_conformal_mesh(Geometry::Layer{1});

    std::cout << "Exporting Mesh\n";
    Geometry::ExportVTK<FoundationMeshType> exporter(refined_mesh);

    exporter.write("concave_refinement");
  }

  void main(int argc, char** argv)
  {
    SimpleArgParser args(argc, argv);

    args.support("program");
    args.support("template_set");

    std::deque<std::pair<int, String>> unsupported = args.query_unsupported();
    if(!unsupported.empty())
    {
      std::cerr << "\n";
      for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
        std::cerr << "ERROR: unsupported option #" << (*it).first << " '--" << (*it).second << "'\n";
      return;
    }

    if(args.check("program") != 1)
    {
      std::cout << "ERROR: mandatory program parameter is missing!\n";
      std::cout << "Valid programs are: 2d_templates, 3d_templates\n";
    }

    if(args.check("template_set") != 1)
    {
      std::cout << "ERROR: mandatory template_set parameter is missing!\n";
      std::cout << "Valid template sets are: schneiders, schneiders_new\n";
    }

    String program;
    args.parse("program", program);

    String template_set;
    args.parse("template_set", template_set);

    std::cout << "Using template set: " << template_set << "\n";

    if(template_set.compare_no_case("schneiders") == 0 && program.compare_no_case("2d_templates") == 0)
    {
      export_templates_2d<Geometry::SchneidersTemplates>();
    }
    if(template_set.compare_no_case("schneiders") == 0 && program.compare_no_case("3d_templates") == 0)
    {
      export_templates_3d<Geometry::SchneidersTemplates>();
    }
    if(template_set.compare_no_case("sunzhaoma") == 0 && program.compare_no_case("2d_templates") == 0)
    {
      export_templates_2d<Geometry::SunZhaoMaTemplates>();
    }
    if(template_set.compare_no_case("sunzhaoma") == 0 && program.compare_no_case("3d_templates") == 0)
    {
      export_templates_3d<Geometry::SunZhaoMaTemplates>();
    }
    if(template_set.compare_no_case("sunzhaoma") == 0 && program.compare_no_case("concave_refinement") == 0)
    {
      concave_refinement<Geometry::SunZhaoMaTemplates>();
    }
    if(template_set.compare_no_case("sunzhaoma_expanded") == 0 && program.compare_no_case("concave_refinement") == 0)
    {
      concave_refinement<Geometry::SunZhaoMaExpansionTemplates>();
    }
  }
} // namespace DbgAdaptiveMesh

int main(int argc, char** argv)
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  DbgAdaptiveMesh::main(argc, argv);

  return 0;
}
