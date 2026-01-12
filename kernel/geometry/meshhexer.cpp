// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>

#ifdef FEAT_HAVE_MESHHEXER

#include <kernel/geometry/adaptive_mesh.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/meshhexer.hpp>
#include <kernel/geometry/subdivision_levels.hpp>
#include <kernel/shape.hpp>
#include <kernel/util/assertion.hpp>

#include <meshhexer/meshhexer.hpp>
#include <meshhexer/types.hpp>

namespace FEAT::Geometry
{
  MeshHexerSurfaceMesh::~MeshHexerSurfaceMesh()
  {
    delete static_cast<MeshHexer::SurfaceMesh*>(_surface_mesh);
  }

  bool MeshHexerSurfaceMesh::is_closed() const
  {
    return static_cast<MeshHexer::SurfaceMesh*>(_surface_mesh)->is_closed();
  }

  bool MeshHexerSurfaceMesh::is_wound_consistently() const
  {
    return static_cast<MeshHexer::SurfaceMesh*>(_surface_mesh)->is_wound_consistently();
  }

  bool MeshHexerSurfaceMesh::is_outward_oriented() const
  {
    return static_cast<MeshHexer::SurfaceMesh*>(_surface_mesh)->is_outward_oriented();
  }

  double MeshHexerSurfaceMesh::minimal_aspect_ratio() const
  {
    return static_cast<MeshHexer::SurfaceMesh*>(_surface_mesh)->minimal_aspect_ratio();
  }

  double MeshHexerSurfaceMesh::maximal_aspect_ratio() const
  {
    return static_cast<MeshHexer::SurfaceMesh*>(_surface_mesh)->maximal_aspect_ratio();
  }

  std::vector<Gap> MeshHexerSurfaceMesh::gaps()
  {
    std::vector<MeshHexer::Gap> gaps = static_cast<MeshHexer::SurfaceMesh*>(_surface_mesh)->gaps();

    // Convert to FEAT type
    std::vector<Gap> result;
    result.reserve(gaps.size());
    for(const auto gap : gaps)
    {
      result.emplace_back(gap.diameter, gap.confidence);
    }

    return result;
  }

  Gap MeshHexerSurfaceMesh::min_gap()
  {
    MeshHexer::Gap gap = static_cast<MeshHexer::SurfaceMesh*>(_surface_mesh)->min_gap();
    return {gap.diameter, gap.confidence};
  }

  template<typename DataType_>
  ConformalMesh<Shape::Hexahedron, 3, DataType_> MeshHexerSurfaceMesh::fbm_mesh(
    DataType_ x_min,
    DataType_ y_min,
    DataType_ z_min,
    DataType_ x_max,
    DataType_ y_max,
    DataType_ z_max,
    std::uint64_t levels)
  {
    using AdaptiveMeshType = AdaptiveMesh<SchneidersTemplates, Shape::Hexahedron, 3, DataType_>;

    auto data = prepare_fbm_mesh<DataType_>(x_min, y_min, z_min, x_max, y_max, z_max, levels);

    AdaptiveMeshType amesh(data.first);
    amesh.adapt(data.second, AdaptiveMeshType::ImportBehaviour::All);

    return amesh.to_conformal_mesh(Layer{amesh.num_layers() - 1});
  }

  template<typename DataType_>
  ConformalMesh<Shape::Hexahedron, 3, DataType_> MeshHexerSurfaceMesh::fbm_mesh(std::uint64_t levels)
  {
    const MeshHexer::BoundingBox bb = static_cast<MeshHexer::SurfaceMesh*>(_surface_mesh)->bounding_box();
    return fbm_mesh<DataType_>(
      DataType_(bb.min.x),
      DataType_(bb.min.y),
      DataType_(bb.min.z),
      DataType_(bb.max.x),
      DataType_(bb.max.y),
      DataType_(bb.max.z),
      levels);
  }

  template<typename DataType_>
  std::pair<ConformalMesh<Shape::Hexahedron, 3, DataType_>, SubdivisionLevels> MeshHexerSurfaceMesh::prepare_fbm_mesh(
    DataType_ x_min,
    DataType_ y_min,
    DataType_ z_min,
    DataType_ x_max,
    DataType_ y_max,
    DataType_ z_max,
    std::uint64_t levels)
  {
    using MeshType = ConformalMesh<Shape::Hexahedron, 3, DataType_>;
    using VertexType = typename MeshType::VertexType;

    MeshHexer::BoundingBox bb{
      {double(x_min), double(y_min), double(z_min)},
      {double(x_max), double(y_max), double(z_max)}};
    MeshHexer::FBMMeshSettings settings{bb, levels};

    MeshHexer::VolumeMesh vmesh = static_cast<MeshHexer::SurfaceMesh*>(_surface_mesh)->fbm_mesh(settings);

    SubdivisionLevels sdls(vmesh.num_vertices());
    for(Index i(0); i < vmesh.num_vertices(); i++)
    {
      sdls[i] = vmesh.subdivision_level(i);
    }

    std::array<Index, 4> size{vmesh.num_vertices(), vmesh.num_edges(), vmesh.num_faces(), vmesh.num_cells()};
    ConformalMesh<Shape::Hexahedron, 3, DataType_> mesh(size.data());

    // Copy vertices
    {
      auto& vertices = mesh.get_vertex_set();
      for(Index i(0); i < vmesh.num_vertices(); i++)
      {
        MeshHexer::Point p = vmesh.vertex(i);
        vertices[i] = VertexType{p.x, p.y, p.z};
      }
    }

    // Copy edges
    {
      auto& v_at_e = mesh.template get_index_set<1, 0>();
      for(Index i(0); i < vmesh.num_edges(); i++)
      {
        v_at_e(i, 0) = vmesh.edge(i, 0);
        v_at_e(i, 1) = vmesh.edge(i, 1);
      }
    }

    // Copy faces
    {
      auto& v_at_f = mesh.template get_index_set<2, 0>();
      for(Index i(0); i < vmesh.num_faces(); i++)
      {
        for(int j(0); j < 4; j++)
        {
          v_at_f(i, j) = vmesh.face(i, std::size_t(j));
        }
      }
    }

    // Copy cells
    {
      auto& v_at_c = mesh.template get_index_set<3, 0>();
      for(Index i(0); i < vmesh.num_cells(); i++)
      {
        for(int j(0); j < 8; j++)
        {
          v_at_c(i, j) = vmesh.cell(i, std::size_t(j));
        }
      }
    }

    mesh.deduct_topology_from_top();

    return std::make_pair(std::move(mesh), std::move(sdls));
  }

  template<typename DataType_>
  std::pair<ConformalMesh<Shape::Hexahedron, 3, DataType_>, SubdivisionLevels>
  MeshHexerSurfaceMesh::prepare_fbm_mesh(std::uint64_t levels)
  {
    const MeshHexer::BoundingBox bb = static_cast<MeshHexer::SurfaceMesh*>(_surface_mesh)->bounding_box();
    return prepare_fbm_mesh<DataType_>(
      DataType_(bb.min.x),
      DataType_(bb.min.y),
      DataType_(bb.min.z),
      DataType_(bb.max.x),
      DataType_(bb.max.y),
      DataType_(bb.max.z),
      levels);
  }

  void MeshHexerSurfaceMesh::write_to_file(const std::string& filename)
  {
    static_cast<MeshHexer::SurfaceMesh*>(_surface_mesh)->write_to_file(filename);
  }

  MeshHexerSurfaceMesh load_from_file(const String& filename, bool triangulate)
  {
    auto result = MeshHexer::load_from_file(filename, triangulate);

    if(result.is_err())
    {
      std::cerr << result.err_ref() << "\n";
      XABORTM("Reading MeshHexerSurfaceMesh failed!");
    }

    auto* mesh = new MeshHexer::SurfaceMesh();
    *mesh = std::move(result).take_ok();

    return MeshHexerSurfaceMesh(mesh);
  }

  template ConformalMesh<Shape::Hexahedron, 3, double>
  MeshHexerSurfaceMesh::fbm_mesh(double, double, double, double, double, double, std::uint64_t);
  template ConformalMesh<Shape::Hexahedron, 3, float>
  MeshHexerSurfaceMesh::fbm_mesh(float, float, float, float, float, float, std::uint64_t);

  template std::pair<ConformalMesh<Shape::Hexahedron, 3, double>, SubdivisionLevels>
  MeshHexerSurfaceMesh::prepare_fbm_mesh(double, double, double, double, double, double, std::uint64_t);
  template std::pair<ConformalMesh<Shape::Hexahedron, 3, float>, SubdivisionLevels>
  MeshHexerSurfaceMesh::prepare_fbm_mesh(float, float, float, float, float, float, std::uint64_t);

  template ConformalMesh<Shape::Hexahedron, 3, double> MeshHexerSurfaceMesh::fbm_mesh(std::uint64_t);
  template ConformalMesh<Shape::Hexahedron, 3, float> MeshHexerSurfaceMesh::fbm_mesh(std::uint64_t);

  template std::pair<ConformalMesh<Shape::Hexahedron, 3, double>, SubdivisionLevels>
  MeshHexerSurfaceMesh::prepare_fbm_mesh(std::uint64_t);
  template std::pair<ConformalMesh<Shape::Hexahedron, 3, float>, SubdivisionLevels>
  MeshHexerSurfaceMesh::prepare_fbm_mesh(std::uint64_t);
} // namespace FEAT::Geometry

#endif
