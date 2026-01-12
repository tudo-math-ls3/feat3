// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/geometry/mesh_file_writer.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/meshhexer.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/string.hpp>

#include <fstream>

namespace FBMMeshGenerator
{
  using namespace FEAT;

  namespace {
    int main(int argc, char** argv)
    {
      FEAT::SimpleArgParser parser(argc, argv);

      parser.support("surface-mesh", "Surface mesh to create an FBM mesh for");
      parser.support("level", "Intended maximum level of mesh-hierarchy");
      parser.support("bounding-box", "Manually set bounding box for fbm mesh. Uses bounding box of surface mesh otherwise");
      parser.support("output", "Output file path");

      std::deque<std::pair<int, String>> unsupported_args = parser.query_unsupported();
      if (!unsupported_args.empty())
      {
        for (const std::pair<int, String> &unsupported : parser.query_unsupported())
        {
          std::cerr << "Unsupported argument " << unsupported.second << " in position " << unsupported.first << ".\n";
        }

        return 1;
      }

      String surface_mesh_path;
      if(const auto *arg = parser.query("surface-mesh"))
      {
        if(arg->second.size() != 1)
        {
          std::cerr << "Option --surface-mesh expects exactly one parameter.\n";
          return 1;
        }

        surface_mesh_path = arg->second.front();
      }
      else
      {
        std::cerr << "Option --surface-mesh is required.\n";
        return 1;
      }

      std::uint64_t level(1);
      if(const auto* arg = parser.query("level"))
      {
        if(arg->second.size() != 1)
        {
          std::cerr << "Option --level expects exactly one parameter.\n";
          return 1;
        }

        parser.parse("level", level);
      }
      else
      {
        std::cerr << "Option --level is required.\n";
        return 1;
      }

      bool custom_bb = false;
      Real min_x(0.0);
      Real min_y(0.0);
      Real min_z(0.0);
      Real max_x(0.0);
      Real max_y(0.0);
      Real max_z(0.0);

      if(const auto* arg = parser.query("bounding-box"))
      {
        if(arg->second.size() != 6)
        {
          std::cerr << "Option --bounding-box expects exactly six parameters ([min_x, min_y, min_z, max_x, max_y, max_z]).\n";
          return 1;
        }

        custom_bb = true;
        parser.parse("bounding-box", min_x, min_y, min_z, max_x, max_y, max_z);
      }

      String output("fbm.xml");
      if(const auto* arg = parser.query("output"))
      {
        if(arg->second.size() != 1)
        {
          std::cerr << "Option --output expects exactly one parameter.\n";
          return 1;
        }

        parser.parse("output", output);
      }

      std::ofstream output_stream(output);
        Geometry::MeshFileWriter writer(output_stream);

      Geometry::MeshHexerSurfaceMesh surface_mesh = Geometry::load_from_file(surface_mesh_path);
      if(custom_bb)
      {
        auto mesh = surface_mesh.fbm_mesh<Real>(min_x, min_y, min_z, max_x, max_y, max_z, level);
        Geometry::RootMeshNode<decltype(mesh)> node(std::make_unique<decltype(mesh)>(std::move(mesh)));
        writer.write(&node);
      }
      else
      {
        auto mesh = surface_mesh.fbm_mesh<Real>(level);
        Geometry::RootMeshNode<decltype(mesh)> node(std::make_unique<decltype(mesh)>(std::move(mesh)));
        writer.write(&node);
      }

      return 0;
    }
  }
}

int main(int argc, char** argv)
{
  FBMMeshGenerator::main(argc, argv);
}
