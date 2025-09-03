// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/mesh_file_writer.hpp>
#include <kernel/geometry/parsed_hit_test_factory.hpp>
#include <kernel/geometry/boundary_factory.hpp>

using namespace FEAT;
using namespace FEAT::Geometry;


static void display_help()
{
  std::cout <<
    "\nmesh-part-gen: Creates mesh-parts for FEAT mesh files based on hit test functions\n\n"
#ifndef FEAT_HAVE_FPARSER
    "WARNING: You have configured FEAT without the 'fparser' third-party library, which\n"
    "is required for the --hit-test functionality, so this feature is not available for\n"
    "this build.\n\n"
#endif
    "Mandatory arguments:\n"
    "--------------------\n"
    " --mesh <path to mesh file(s)>\n"
    "Specifies the sequence of input mesh files paths.\n\n"
    " --out <path to mesh-part file>\n"
    "Specifies the output mesh-part file path.\n\n"
    "Optional arguments:\n"
    "--------------------\n"
    " --out-complete\n"
    "Specifies that the output mesh file should also contain the contents of the input\n"
    "mesh file(s) along with the newly created mesh-parts. If not given, the output mesh\n"
    "file will only contain the newly created mesh-parts but nothing that was read in\n"
    "from the input mesh file(s).\n\n"
    " --bnd-whole [<name>]\n"
    "Specifies that a single mesh part for the whole boundary is to be generated.\n"
    "If specified without <name>, the generated mesh-part will be named 'bnd'.\n\n"
    " --hit-test <name1> <formula1> [<name2> <formula2> ...]\n"
    "Specifies a set of name-formula argument pairs which are used to generate mesh-parts\n"
    "by using the Geometry::ParsedHitTestFactory class. The first component of each pair\n"
    "specifies the name for the mesh-part, whereas the second component specifies the formula\n"
    "in x,y,z coordinates, which is to be used for the hit test of the mesh part.\n"
    "A vertex or edge/face/cell will be contained in the mesh-part if the formula evaluates\n"
    "to a positive value in its coordinates or midpoint coordinates, respectively.\n"
    "You can also use the > and < operators to compare coordinates explicitly.\n"
    "It is strongly advised to put all formulae in quotations marks, because otherwise your\n"
    "shell/script interpreter might digest the formulae thus leading to incorrect results.\n\n"
    " --bnd-complement <name1> <parts1> [<name2> <parts2> ...]\n"
    "Specifies that a boundary mesh-part is to be created that complements a set of other\n"
    "given mesh-parts. This option can complement mesh-parts that have been read in from\n"
    "the input mesh file or created via the --hit-test option. If the complement to a set\n"
    "of mesh-parts is to be computed, the names of the mesh-parts have to be given as a\n"
    "single parameter by enclosing them in quotation marks, separated by spaces.\n\n"
    " --help\n"
    "Displays this message\n\n"
    "Example:\n"
    "--------\n"
    "Assume that the input mesh is a 3D mesh that has an inflow at z=0 and an outflow at z=1\n"
    "then we can use the following call to create 3 mesh-parts for the inflow, the outflow\n"
    "as well as the remaining part of the boundary, which we call the 'wall':\n\n"
    "  mesh-part-gen --mesh input_mesh.xml --out output_mesh.xml \\\n"
    "    --hit-test \"bnd:inflow\" \"z<0.0001\" \"bnd:outflow\" \"z>0.9999\" \\\n"
    "    --bnd-complement \"bnd:wall\" \"bnd:inflow bnd:outflow\"\n";
  std::cout.flush();
}

template<typename MeshType_>
bool build_hit_test_meshparts(Geometry::RootMeshNode<MeshType_>& mesh_node_out,
  const Geometry::RootMeshNode<MeshType_>& mesh_node_in, const std::map<String, String>& hit_formulae)
{
#ifdef FEAT_HAVE_FPARSER
  // loop over all hit-test formulae
  for(auto it = hit_formulae.begin(); it != hit_formulae.end(); ++it)
  {
    // get part name and formula
    String name = it->first;
    String formula = it->second;

    std::cout << "Creating mesh-part '" << name << "' by hit-test formula '" << formula << "'...\n";

    try
    {
      // try to parse the formula
      Geometry::ParsedHitTestFactory<MeshType_> factory(*mesh_node_in.get_mesh());

      // parse formula
      factory.parse(formula);

      // create mesh part
      auto mesh_part = factory.make_unique();

      // print target counts
      std::cout << "Target entity counts:";
      for(int i(0); i <= MeshType_::shape_dim; ++i)
        std::cout << ' ' << mesh_part->get_num_entities(i);
      std::cout << std::endl;

      // add to mesh node
      mesh_node_out.add_mesh_part(name, std::move(mesh_part));
    }
    catch(std::exception& exc)
    {
      std::cerr << "ERROR: in boundary function formula '" << formula << "'\n";
      std::cerr << exc.what() << std::endl;
      return false;
    }
  }

  return true;
#else // no FEAT_HAVE_FPARSER
  std::cout << "ERROR: Creation of mesh parts based on hit-test is only available if FEAT is\n";
  std::cout << "configured and build with the 'fparser' third-party library!\n";
  return false;
#endif // FEAT_HAVE_FPARSER
}

template<typename MeshType_>
bool build_bnd_complement_meshparts(Geometry::RootMeshNode<MeshType_>& mesh_node_out,
  const Geometry::RootMeshNode<MeshType_>& mesh_node_in, const std::map<String, String>& bnd_complements)
{
  // loop over boundary complements
  for(auto it = bnd_complements.begin(); it != bnd_complements.end(); ++it)
  {
    // get part name and formula
    String name = it->first;
    std::deque<String> others = it->second.split_by_whitespaces();

    std::cout << "Creating boundary mesh-part '" << name << "' complementary to mesh-parts '"
      << stringify_join(others, "', '") << "'...\n";

    // create masked boundary factory
    Geometry::MaskedBoundaryFactory<MeshType_> factory(*mesh_node_in.get_mesh());

    // add masked parts
    for(const auto& part_name : others)
    {
      const Geometry::MeshPart<MeshType_>* part = nullptr;
      if(part == nullptr)
        part = mesh_node_in.find_mesh_part(part_name);
      if(part == nullptr)
        part = mesh_node_out.find_mesh_part(part_name);
      if(part == nullptr)
      {
        std::cout << "ERROR: mesh-part '" << part_name << "' for boundary complement part '" << name << "' not found!\n";
        return false;
      }
      factory.add_mask_meshpart(*part);
    }

    // compile the factory
    factory.compile();

    // create mesh part
    auto mesh_part = factory.make_unique();

    // print target counts
    std::cout << "Target entity counts:";
    for(int i(0); i <= MeshType_::shape_dim; ++i)
      std::cout << ' ' << mesh_part->get_num_entities(i);
    std::cout << std::endl;

    // add to mesh node
    mesh_node_out.add_mesh_part(name, std::move(mesh_part));
  }

  return true;
}

template<typename Mesh_>
int run_xml(SimpleArgParser& args, Geometry::MeshFileReader& mesh_reader)
{
  // write complete output file?
  const bool out_complete = (args.check("out-complete") < 0);

  // build hit-formulae map
  std::map<String, String> hit_formulae;
  if(args.check("hit-test") > 0)
  {
    std::deque<String> parts = args.query("hit-test")->second;
    if(parts.size() % 2u != 0u)
    {
      std::cerr << "ERROR: invalid number of parameters for option --hit-test, expected an even count\n";
      return 1;
    }

    for(std::size_t i(0); i < parts.size(); i += 2u)
      hit_formulae.emplace(parts[i], parts[i+1u]);
  }

  // build boundary complements map
  std::map<String, String> bnd_complements;
  if(args.check("bnd-complement") > 0)
  {
    std::deque<String> parts = args.query("bnd-complement")->second;
    if(parts.size() % 2u != 0u)
    {
      std::cerr << "ERROR: invalid number of parameters for option --bnd-complement, expected an even count\n";
      return 1;
    }

    for(std::size_t i(0); i < parts.size(); i += 2u)
      bnd_complements.emplace(parts[i], parts[i+1u]);
  }

  // create an empty atlas and a root mesh node
  auto atlas = Geometry::MeshAtlas<Mesh_>::make_unique();
  std::shared_ptr<Geometry::RootMeshNode<Mesh_>> node;
  Geometry::PartitionSet part_set;

  // try to parse the mesh file
#ifndef DEBUG
  try
#endif
  {
    std::cout << "Parsing mesh files...\n";
    node = mesh_reader.parse(*atlas, &part_set);
  }
#ifndef DEBUG
  catch(std::exception& exc)
  {
    std::cerr << "ERROR: " << exc.what() << std::endl;
    return 1;
  }
  catch(...)
  {
    std::cerr << "ERROR: unknown exception\n";
    return 1;
  }
#endif

  // adapt coarse mesh
  node->adapt();

  // create output node unless we want to write the complete mesh
  std::shared_ptr<Geometry::RootMeshNode<Mesh_>> node_out = node;
  if(args.check("out-complete") < 0)
    node_out = Geometry::RootMeshNode<Mesh_>::make_unique(nullptr, nullptr);

  // build hit-test mesh-parts
  if(!build_hit_test_meshparts(*node_out, *node, hit_formulae))
    return false;

  // build complementary mesh-parts
  if(!build_bnd_complement_meshparts(*node_out, *node, bnd_complements))
    return false;

  // export mesh part node to file
  String out_name;
  if(args.parse("out", out_name) > 0)
  {
    std::cout << "Writing " << (out_complete ? "complete mesh" : "mesh-parts") << " to file '" << out_name << "'...\n";
    std::ofstream ofs(out_name, std::ios_base::out);
    Geometry::MeshFileWriter writer(ofs);
    if(out_complete)
      writer.write(node_out.get(), atlas.get(), &part_set);
    else
      writer.write(node_out.get());
    std::cout << "Done!\n";
  }

  return 0;
}

int main(int argc, char* argv[])
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  // This is the list of all supported meshes that could appear in the mesh file
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 2, Real> S2M2D;
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 3, Real> S2M3D;
  typedef Geometry::ConformalMesh<Shape::Simplex<3>, 3, Real> S3M3D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 1, Real> H1M1D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 2, Real> H1M2D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 3, Real> H1M3D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, Real> H2M2D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 3, Real> H2M3D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<3>, 3, Real> H3M3D;

  SimpleArgParser args(argc, argv);

  // need help?
  if((argc < 2) || (args.check("help") > -1))
  {
    display_help();
    return 0;
  }

  args.support("mesh");
  args.support("out");
  args.support("hit-test");
  args.support("bnd-complement");
  args.support("out-complete");

  // check for unsupported options
  auto unsupported = args.query_unsupported();
  if( !unsupported.empty() )
  {
    // print all unsupported options to cerr
    for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
      std::cerr << "ERROR: unsupported option '--" << (*it).second << "'\n";

    display_help();
    return 1;
  }

  int num_mesh_files = args.check("mesh");
  if(num_mesh_files < 1)
  {
    std::cerr << "ERROR: You have to specify at least one mesh file with --mesh <files...>\n";
    display_help();
    return 1;
  }

  // get our filename deque
  auto mpars = args.query("mesh");
  XASSERT(mpars != nullptr);
  const std::deque<String>& filenames = mpars->second;

  // create an empty mesh file reader
  Geometry::MeshFileReader mesh_reader;
  mesh_reader.add_mesh_files(filenames);

  // read root markup
  mesh_reader.read_root_markup();

  // get mesh type
  const String mtype = mesh_reader.get_meshtype_string();

  std::cout << "Mesh Type: " << mtype << std::endl;

  if(mtype == "conformal:hypercube:1:1")
    return run_xml<H1M1D>(args, mesh_reader);
  if(mtype == "conformal:hypercube:1:2")
    return run_xml<H1M2D>(args, mesh_reader);
  if(mtype == "conformal:hypercube:1:3")
    return run_xml<H1M3D>(args, mesh_reader);
  if(mtype == "conformal:hypercube:2:2")
    return run_xml<H2M2D>(args, mesh_reader);
  if(mtype == "conformal:hypercube:2:3")
    return run_xml<H2M3D>(args, mesh_reader);
  if(mtype == "conformal:hypercube:3:3")
    return run_xml<H3M3D>(args, mesh_reader);
  if(mtype == "conformal:simplex:2:2")
    return run_xml<S2M2D>(args, mesh_reader);
  if(mtype == "conformal:simplex:2:3")
    return run_xml<S2M3D>(args, mesh_reader);
  if(mtype == "conformal:simplex:3:3")
    return run_xml<S3M3D>(args, mesh_reader);

  std::cout << "ERROR: unsupported mesh type!\n";

  return 1;
}
