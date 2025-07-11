// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/runtime.hpp>
#include <kernel/geometry/adaptive_mesh.hpp>
#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/mesh_file_writer.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/util/simple_arg_parser.hpp>

#include <cstring>

namespace MeshRefiner
{
  using namespace FEAT;

  template <typename Shape, typename TemplateSet>
  int run(std::size_t full_pre, String sdl_file, std::size_t full_post, String output_file, Geometry::MeshFileReader &reader)
  {
    using MeshType = Geometry::ConformalMesh<Shape>;
    using RootMeshNodeType = Geometry::RootMeshNode<MeshType>;

    Geometry::MeshAtlas<MeshType> mesh_atlas;

    std::unique_ptr<RootMeshNodeType> mesh_node;

    try
    {
      mesh_node = reader.parse(mesh_atlas);
    }
    catch (const std::exception &e)
    {
      std::cerr << "Error parsing mesh files: " << e.what() << "\n";
      return 1;
    }

    for (std::size_t lvl(0); lvl < full_pre; lvl++)
    {
      mesh_node = mesh_node->refine_unique();
    }

    if (!sdl_file.empty())
    {
      LAFEM::DenseVector<std::uint64_t> sdls;
      sdls.read_from(LAFEM::FileMode::fm_mtx, sdl_file);

      mesh_node = mesh_node->template refine_partial_unique<TemplateSet>(sdls);
    }

    for (std::size_t lvl(0); lvl < full_post; lvl++)
    {
      mesh_node = mesh_node->refine_unique();
    }

    std::ofstream output(output_file);

    if (output.fail())
    {
      std::cerr << "Opening output file for writing failed.\n";
      return 1;
    }

    Geometry::MeshFileWriter writer(output);
    writer.write(mesh_node.get(), &mesh_atlas);

    return 0;
  }

  const char *help_string =
      "mesh-refiner: Apply refinements to mesh files\n"
      "\n"
      "Usage:\n"
      "mesh-refiner [-h|--help] <mesh-files> [--template-set <Schneiders|SunZhaoMa>] [--full-pre <n>] [--adaptive <sdl-file>] [--full-post <n>] [--output-file <path>]\n"
      "\n"
      "Options:\n"
      "\t--adaptive <sdl-file>\n"
      "\t\tApply a conformal partial refinement to the mesh.\n"
      "\t\t<sdl-file> must be a vector in the MatrixMarket format that assigns a subdivision level to each vertex.\n"
      "\t--full-pre <n>\n"
      "\t\tApply <n> global 2-refinements to the mesh, before doing any partial refinement.\n"
      "\t--full-post <n>\n"
      "\t\tApply <n> global 2-refinements to the mesh, after doing any partial refinement.\n"
      "\t--output-file <path>\n"
      "\t\tOutput file. Defaults to refined.xml\n"
      "\t--template-set <Schneiders|SuZhaoMa>\n"
      "\t\tChoose the template set for partial refinements\n";

  int main(int argc, char **argv)
  {
    if (argc == 1)
    {
      std::cerr << help_string;
      return 1;
    }

    if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)
    {
      std::cout << help_string;
      return 0;
    }

    SimpleArgParser parser(argc, argv);

    parser.support("adaptive");
    parser.support("full-pre");
    parser.support("full-post");
    parser.support("output-file");
    parser.support("template-set");

    std::deque<std::pair<int, String>> unsupported_args = parser.query_unsupported();
    if (!unsupported_args.empty())
    {
      for (const std::pair<int, String> &unsupported : parser.query_unsupported())
      {
        std::cerr << "Unsupported argument " << unsupported.second << " in position " << unsupported.first << ".\n";
      }

      std::cerr << "\n";
      std::cerr << help_string;
      return 1;
    }

    std::deque<String> mesh_files;
    std::size_t full_pre = 0;
    std::size_t full_post = 0;
    String sdl_file;
    String template_set("Schneiders");
    String output_file("refined.xml");

    for (int i(1); i < parser.num_skipped_args(); i++)
    {
      mesh_files.push_back(String(argv[i]));
    }

    int full_pre_count = parser.check("full-pre");
    if (full_pre_count != -1)
    {
      if (full_pre_count != 1)
      {
        std::cerr << "Option --full-pre expects exactly one parameter.\n";
        return 1;
      }

      parser.parse("full-pre", full_pre);
    }

    int full_post_count = parser.check("full-post");
    if (full_post_count != -1)
    {
      if (full_post_count != 1)
      {
        std::cerr << "Option --full-post expects exactly one parameter.\n";
        return 1;
      }

      parser.parse("full-post", full_post);
    }

    int sdl_file_count = parser.check("adaptive");
    if (sdl_file_count != -1)
    {
      if (sdl_file_count != 1)
      {
        std::cerr << "Option --adaptive expects exactly one parameter.\n";
        return 1;
      }

      parser.parse("adaptive", sdl_file);
    }

    int template_set_count = parser.check("template-set");
    if (template_set_count != -1)
    {
      if (template_set_count != 1)
      {
        std::cerr << "Option --template-set expects exactly one parameter.\n";
        return 1;
      }

      parser.parse("template-set", template_set);

      if (template_set != "Schneiders" && template_set != "SunZhaoMa")
      {
        std::cerr << "Invalid parameters for option --template-set\n";
        std::cerr << "Valid options are Schneiders and SunZhaoMa\n";
        return 1;
      }
    }

    int output_file_count = parser.check("output-file");
    if (output_file_count != -1)
    {
      if (output_file_count != 1)
      {
        std::cerr << "Option --output-file expects exactly one parameter.\n";
        return 1;
      }

      parser.parse("output-file", output_file);
    }

    Geometry::MeshFileReader mesh_reader;

    try
    {
      mesh_reader.add_mesh_files(mesh_files);
    }
    catch (const std::exception &e)
    {
      std::cerr << "Error reading mesh files: " << e.what() << "\n";
      return 1;
    }

    try
    {
      mesh_reader.read_root_markup();
    }
    catch (const std::exception &e)
    {
      std::cerr << "Error parsing mesh file root markup: " << e.what() << "\n";
      return 1;
    }

    const String mesh_type = mesh_reader.get_meshtype_string();

    if (mesh_type.empty())
    {
      std::cerr << "Error: Found no mesh_type in mesh files. Did you supply all required mesh files?\n";
      return 1;
    }

    if (mesh_type != "conformal:hypercube:2:2" && mesh_type != "conformal:hypercube:3:3")
    {
      std::cerr << "Error: Unsupported mesh type " << mesh_type << "\n";
      return 1;
    }

    if (mesh_type == "conformal:hypercube:2:2")
      if (template_set == "Schneiders")
        return run<Shape::Hypercube<2>, Geometry::SchneidersTemplates>(full_pre, sdl_file, full_post, output_file, mesh_reader);
      else
        return run<Shape::Hypercube<2>, Geometry::SunZhaoMaExpansionTemplates>(full_pre, sdl_file, full_post, output_file, mesh_reader);
    else if (template_set == "Schneiders")
      return run<Shape::Hypercube<3>, Geometry::SchneidersTemplates>(full_pre, sdl_file, full_post, output_file, mesh_reader);
    else
      return run<Shape::Hypercube<3>, Geometry::SunZhaoMaExpansionTemplates>(full_pre, sdl_file, full_post, output_file, mesh_reader);

    return 0;
  }
} // namespace MeshRefiner

int main(int argc, char **argv)
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);
  return MeshRefiner::main(argc, argv);
}
