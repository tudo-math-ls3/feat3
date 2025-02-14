// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/stop_watch.hpp>
#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/mesh_file_writer.hpp>
#include <kernel/geometry/parti_iterative.hpp>
#include <kernel/geometry/parti_parmetis.hpp>
#include <kernel/geometry/parti_zoltan.hpp>

using namespace FEAT;
using namespace FEAT::Geometry;

static constexpr int pad_len = 32;

void display_help(const Dist::Comm& comm, SimpleArgParser& args)
{
  //          123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
  comm.print("mesh-partitioner: Creates a partitioning for a given mesh file and saves it to another mesh file");
  comm.print("\nThis tool can be used to create a partitioning file for a given mesh file, which is then saved");
  comm.print("as a separate mesh file and which can be passed to MPI parallel applications in addition to the");
  comm.print("original mesh file to ensure that these applications picks one of the precomputed partitionings");
  comm.print("instead of applying a partitioner at simulation time.\n");

  comm.print("Please note that this tool can only be really useful if FEAT was configured and linked against the");
  comm.print("ParMETIS and/or Zoltan libraries, because FEAT does not provide any 'real' partitioners by its own");
  comm.print("except for a very experimental genetic partitioner that should not be used for real world use cases.");
  comm.print("In consequence, this tool must be compiled with MPI as both of these libraries require MPI.\n");

  comm.print("This tool allows you to specify the number of desired partitions/patches independently of the number");
  comm.print("of MPI processes that this tool is run with, however, not all MPI processes may be used to run the");
  comm.print("partitioner if the workload is too small to distribute over all MPI processes.\n");

  comm.print("This tool allows you to specify the refinement level of the mesh that is to be partitioned, if the");
  comm.print("unrefined mesh stored in the mesh file that you want to partition is too coarse for the number of");
  comm.print("patches that you want to create. The chosen partitioning level is then also written to the output");
  comm.print("file and the PartiDomainControl will also refine the mesh to the required level if it chooses to");
  comm.print("select the corresponding partitioning.\n");

  comm.print("By default, the partitioning is named 'auto', which tells the PartiDomainControl, which is used to");
  comm.print("read in the partition mesh file, that the corresponding partitioning can be chosen automatically");
  comm.print("if the number of MPI processes matches the number of patches in the partitioning. If you choose a");
  comm.print("different name for the partitioning, you will have to tell the PartiDomainControl explicitly that");
  comm.print("you want to use that particular partitioning by supplying its name to the '--part-extern-name'");
  comm.print("command line parameter of the PartDomainControl object. Unless you indent to try out different");
  comm.print("partitionings for the same number of patches, it is recommended to simply leave the partition name");
  comm.print("as 'auto'.\n");

  comm.print("By default, this tool creates the dual graph for the partitioner based on the facet-adjacency of");
  comm.print("the elements, however, you can also specify that the dual graph is to be defined based on the");
  comm.print("vertex-adjacency of the elements. The former one usually results in nicer partitionings, but you");
  comm.print("might want to try out the vertex-based adjacency if the facet-based element adjacency does not yield");
  comm.print("satisfactory results for the mesh that you are trying to partition.\n");

  comm.print("If you want to visualize the partitioning created by this tool, you can simply use the 'mesh2vtk'");
  comm.print("tool and supply both the input mesh file as well as the partitioning mesh file written by this tool");
  comm.print("as input mesh files and the mesh2vtk tool will write out the partitioning as a cell variable in the");
  comm.print("resulting VTK file. If you have created the partitioning on a refinement level greater than 0, you");
  comm.print("also have to tell the mesh2vtk tool to write out the refined mesh on at least the same level to");
  comm.print("visualize the partitioning because it is not available on lower refinement levels, of course.\n");

  comm.print("This tool supports the following command line parameters:\n");
  comm.print(args.get_supported_help());
}

template<typename MeshType_>
int run(const Dist::Comm& comm, SimpleArgParser& args, Geometry::MeshFileReader& mesh_reader)
{
  static constexpr int shape_dim = MeshType_::shape_dim;
  StopWatch watch_total;
  watch_total.start();

  // parse parameters
  Index parti_level = 0;
  Index num_parts = 0;
  //Index min_elems = 1; // per MPI process
  String parti_name = "auto";
  int parti_prio = 1;
  String out_name;
  bool dual_verts = args.check("dual-by-verts") >= 0;
  bool allow_parmetis = args.check("no-parmetis") < 0;
  bool allow_zoltan = args.check("no-zoltan") < 0;
  bool allow_genetic = args.check("genetic") > 0;
  bool success_genetic = false;
  bool success_parmetis = false;
  bool success_zoltan = false;
  bool success = false;
  int subcomm_size_parmetis = 0;
  int subcomm_size_zoltan = 0;
  double genetic_init = 0.0, genetic_mutate = 0.0;

  if(args.parse("parts", num_parts) < 1)
  {
    comm.print(std::cerr, "ERROR: Failed to parse desired number of partitions!");
    return 1;
  }

  if(args.parse("level", parti_level) < 0)
  {
    comm.print(std::cerr, "ERROR: Failed to parse desired partition level!");
    return 1;
  }

  if(args.parse("out", out_name) < 0)
  {
    comm.print(std::cerr, "ERROR: Failed to parse output filename!");
    return 1;
  }
  if(allow_genetic && (args.parse("genetic", genetic_init, genetic_mutate) < 2))
  {
    comm.print(std::cerr, "ERROR: Failed to parse genetic init or mutate times!");
    return 1;
  }

  // print desired configuration
  comm.print(String("Dual Graph Adjacency").pad_back(pad_len, '.') + ": " + stringify(dual_verts ? "by vertices" : "by facets"));
  comm.print(String("Desired Partitions").pad_back(pad_len, '.') + ": " + stringify(num_parts));
  comm.print(String("Partitioning Level").pad_back(pad_len, '.') + ": " + stringify(parti_level));
  comm.print(String("Partitioning Name").pad_back(pad_len, '.') + ": '" + parti_name + "'");
  comm.print(String("Partitioning Priority").pad_back(pad_len, '.') + ": " + stringify(parti_prio));
#ifdef FEAT_HAVE_PARMETIS
  comm.print(String("ParMETIS Partitioner").pad_back(pad_len, '.') + ": " + stringify(allow_parmetis ? "enabled" : "disabled"));
#else
  comm.print(String("ParMETIS Partitioner").pad_back(pad_len, '.') + ": not available");
#endif
#ifdef FEAT_HAVE_ZOLTAN
  comm.print(String("Zoltan Partitioner").pad_back(pad_len, '.') + ": " + stringify(allow_zoltan ? "enabled" : "disabled"));
#else
  comm.print(String("Zoltan Partitioner").pad_back(pad_len, '.') + ": not available");
#endif
  comm.print(String("Genetic Partitioner").pad_back(pad_len, '.') + ": " + stringify(allow_genetic ? "enabled" : "disabled"));
  comm.print(String("Genetic Init Time").pad_back(pad_len, '.') + ": " + stringify(genetic_init));
  comm.print(String("Genetic Mutate Time").pad_back(pad_len, '.') + ": " + stringify(genetic_mutate));

  // print non-fatal warnings
  if(num_parts < Index(comm.size()))
  {
    comm.print("\nWARNING: Requested less partitions than there are MPI ranks here");
  }

  comm.print_flush();

  // read in the mesh
  MeshAtlas<MeshType_> mesh_atlas;
  typedef RootMeshNode<MeshType_> MeshNodeType;
  std::unique_ptr<MeshNodeType> mesh_node = mesh_reader.parse(mesh_atlas);

  // discard all mesh parts, because we don't need them for this application
  mesh_node->remove_all_mesh_parts();

  // the empty partition object
  PartitionSet parti_set;

  // refine mesh up to desired level and print statistics
  comm.print("\nLevel    #GlobalElems     Mesh Size  |  #OptPartElems     Mesh Size");
  comm.print("--------------------------------------------------------------------");
  Index num_elems = 0u;
  for(Index ilevel(0u); ilevel <= parti_level; ++ilevel)
  {
    // refine mesh node
    if(ilevel > 0u)
      mesh_node = mesh_node->refine_unique();

    num_elems = mesh_node->get_mesh()->get_num_elements();
    std::size_t bytes = mesh_node->bytes();

    // print some statistics
    String s;
    s += stringify(ilevel).pad_front(5) + ":";
    s += stringify(num_elems).pad_front(15);
    s += stringify_bytes(bytes, 3, 9).pad_back(15) + " |";
    s += stringify(num_elems/num_parts).pad_front(15);
    s += stringify_bytes(bytes/num_parts, 3, 9).pad_back(15);
    comm.print(s);
  }

  // check number of elements vs number of partitions
  if(num_elems < num_parts)
  {
    comm.print("\nERROR: There are less elements on the desired level than desired partitions!");
    return 1;
  }

  // get vertices-at-element and faces-at-element index sets
  const auto& vertex_set = mesh_node->get_mesh()->get_vertex_set();
  const auto& verts_at_elem = mesh_node->get_mesh()->template get_index_set<shape_dim, 0>();
  const auto& faces_at_elem = mesh_node->get_mesh()->template get_index_set<shape_dim, shape_dim-1>();

  // prevent unused variable warnings
  (void)vertex_set;
  (void)verts_at_elem;
  (void)faces_at_elem;

  // try ParMETIS first if it is available and allowed by the user
  StopWatch watch_parmetis;
#ifdef FEAT_HAVE_PARMETIS
  if(!success && allow_parmetis)
  {
    watch_parmetis.start();
    comm.print("\nApplying ParMETIS partitioner on level " + stringify(parti_level) + "...");
    comm.print("This might take a while, so please be patient...");
    comm.print_flush();

    PartiParMETIS parti_parmetis(comm, 1u, 0);

    std::vector<Real> weights;
    if(dual_verts)
      success_parmetis = parti_parmetis.execute(verts_at_elem, verts_at_elem, vertex_set, num_parts, weights);
    else
      success_parmetis = parti_parmetis.execute(faces_at_elem, verts_at_elem, vertex_set, num_parts, weights);

    subcomm_size_parmetis = parti_parmetis.get_sub_comm_size();

    if(success_parmetis)
    {
      comm.print("ParMETIS claims it found a valid partitioning!");
      parti_set.add_partition(Partition(parti_parmetis.build_elems_at_rank(), parti_name, parti_prio, int(parti_level)));
      success = true;
    }
    else
    {
      comm.print("ParMETIS admits that it failed to find a valid partitioning!");
    }
    comm.print_flush();

    watch_parmetis.stop();
  }
#else // no ParMETIS
  // this is just to avoid unused variable warnings...
  (void)success_parmetis;
  (void)allow_parmetis;
  (void)subcomm_size_parmetis;
#endif // FEAT_HAVE_PARMETIS

  // try Zoltan next if it is available and allowed by the user
  StopWatch watch_zoltan;
#ifdef FEAT_HAVE_ZOLTAN
  if(!success && allow_zoltan)
  {
    watch_zoltan.start();
    comm.print("\nApplying Zoltan Hypergraph partitioner on level " + stringify(parti_level) + "...");
    comm.print("This might take a while, so please be patient...");
    comm.print_flush();

    PartiZoltan parti_zoltan(comm, 1u, 0);

    std::vector<Real> weights;
    if(dual_verts)
      success_zoltan = parti_zoltan.execute(verts_at_elem, num_parts, weights);
    else
      success_zoltan = parti_zoltan.execute(faces_at_elem, num_parts, weights);

    subcomm_size_zoltan = parti_zoltan.get_sub_comm_size();

    // Zoltan writes to cout in debug mode, so flush on all ranks and wait in barrier...
    std::cout.flush();
    comm.barrier();

    if(success_zoltan)
    {
      comm.print("Zoltan claims it found a valid partitioning!");
      parti_set.add_partition(Partition(parti_zoltan.build_elems_at_rank(), parti_name, parti_prio, int(parti_level)));
      success = true;
    }
    else
    {
      comm.print("Zoltan admits that it failed to find a valid partitioning!");
    }
    comm.print_flush();

    watch_zoltan.stop();
  }
#else // no Zoltan
  // this is just to avoid unused variable warnings...
  (void)success_zoltan;
  (void)allow_zoltan;
  (void)subcomm_size_zoltan;
#endif // FEAT_HAVE_ZOLTAN

  // give the genetic partitioner a try if the user really thinks that this is a good idea
  StopWatch watch_genetic;
  if(!success && allow_genetic)
  {
    watch_genetic.start();
    comm.print("\nApplying genetic partitioner on level " + stringify(parti_level) + "...");
    comm.print("This might take a while, so please be patient...");
    comm.print_flush();

    try
    {
      PartiIterative<MeshType_> parti_genetic(*mesh_node->get_mesh(), comm, num_parts, genetic_init, genetic_mutate);
      comm.print("Genetic partitioner claims it found a valid partitioning!");
      parti_set.add_partition(Partition(parti_genetic.build_elems_at_rank(), parti_name, parti_prio, int(parti_level)));
      success = success_genetic = true;
    }
    catch(...)
    {
      comm.print("Genetic partitioner admits that it failed to find a valid partitioning!");
    }
    comm.print_flush();

    watch_genetic.stop();
  }

  // did anyone find a partitioning?
  if(parti_set.get_partitions().size() == 0u)
  {
    comm.print(std::cerr, "\nERROR: no valid partition found!");
    return 1;
  }

  // get the partitioning
  Partition& parti = parti_set.get_partitions().front();

  // get number of elements per patch
  comm.print("\nValidating partitioning...");
  std::vector<Index> nelp(parti.get_num_patches());
  {
    Index n = parti.get_num_patches();
    const auto& patches = parti.get_patches();
    bool failed = false;
    for(Index i(0); i < n; ++i)
    {
      if((nelp[i] = patches.degree(i)) == 0u)
      {
        comm.print(std::cerr, "\nERROR: partitioner returned an empty patch " + stringify(i));
        failed = true;
      }
    }
    if(failed)
      return 1;
  }
  comm.print("Partitioning seems OK!");

  // print information
  comm.print("\nPartitioning Information:");
  comm.print(String("Total Number of Elements").pad_back(pad_len, '.') + ": " + stringify(parti.get_num_elements()));
  comm.print(String("Total Number of Patches").pad_back(pad_len, '.') + ": " + stringify(parti.get_num_patches()));
  comm.print(String("Average #Elements per Patch").pad_back(pad_len, '.') + ": " + stringify(double(parti.get_num_elements())/double(parti.get_num_patches())));
  if(!nelp.empty())
  {
    // compute and show some balancing statistics
    Index min_epp(nelp.front()), max_epp(nelp.front());
    for(Index k : nelp)
      Math::minimax(k, min_epp, max_epp);
    Index n = parti.get_num_patches();
    std::nth_element(nelp.begin(), nelp.begin() + (n/2u), nelp.end());
    comm.print(String("Median  #Elements per Patch").pad_back(pad_len, '.') + ": " + stringify(nelp[n/2u]));
    comm.print(String("Minimum #Elements per Patch").pad_back(pad_len, '.') + ": " + stringify(min_epp));
    comm.print(String("Maximum #Elements per Patch").pad_back(pad_len, '.') + ": " + stringify(max_epp));
    comm.print(String("Partitioning Balance").pad_back(pad_len, '.') + ": " + stringify_fp_fix(double(min_epp)/double(max_epp), 3));
  }
  comm.print(String("ParMETIS Sub-Comm Size").pad_back(pad_len, '.') + ": " + stringify(subcomm_size_parmetis) + " of " + stringify(comm.size()));
  comm.print(String("Zoltan Sub-Comm Size").pad_back(pad_len, '.') + ": " + stringify(subcomm_size_zoltan) + " of " + stringify(comm.size()));
  comm.print(String("ParMETIS Runtime").pad_back(pad_len, '.') + ": " + watch_parmetis.elapsed_string() + " seconds");
  comm.print(String("Zoltan Runtime").pad_back(pad_len, '.') + ": " + watch_zoltan.elapsed_string() + " seconds");
  comm.print(String("Genetic Runtime").pad_back(pad_len, '.') + ": " + watch_genetic.elapsed_string() + " seconds");
  comm.print_flush();

  // write output file
  if(!out_name.empty())
  {
    comm.print("\nWriting Partitioning to file '" + out_name + "'...");
    if(comm.rank() == 0)
    {
      std::ofstream ofs(out_name, std::ios_base::out);
      {
        MeshFileWriter writer(ofs);
        writer.write_only(parti_set);
      }
      ofs.close();
    }
  }

  watch_total.stop();
  comm.print(String("\nTotal Runtime").pad_back(pad_len, '.') + ": " + watch_total.elapsed_string() + " seconds");

  return 0;
}

int main(int argc, char** argv)
{
  Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  Dist::Comm comm(Dist::Comm::world());

  // This is the list of all supported meshes that could appear in the mesh file
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 2, Real> S2M2D;
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 3, Real> S2M3D;
  typedef Geometry::ConformalMesh<Shape::Simplex<3>, 3, Real> S3M3D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 1, Real> H1M1D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 2, Real> H1M2D;
  //typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 3, Real> H1M3D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, Real> H2M2D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 3, Real> H2M3D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<3>, 3, Real> H3M3D;

  SimpleArgParser args(argc, argv);

  args.support("mesh", "<meshfile>\nMandatory: Specifies the name of the input mesh file to be partitioned.");
  args.support("out", "<partition-file>\nMandatory: Specifies the name of the output partition mesh file.");
  args.support("parts", "<n>\nMandatory: Specifies the number of partitions/patches to create.");
  args.support("dual-by-verts", "\nOptional: Specifies that the dual graph is to be computed by vertices rather than by facets.");
  args.support("level", "<n>\nOptional: Specifies the mesh partitioning level; defaults to 0.");
  args.support("name", "<name>\nOptional: Specifies the name of the partitioning; defaults to 'auto'.");
  args.support("prio", "<priority>\nOptional: Specifies the priority of the partitioning; defaults to 1.");
  args.support("no-2lvl", "\nOptional: Specifies that the 2-level partitioner should not be used.");
  args.support("genetic", "<time-init> <time-mutate>\nOptional: Specifies that the genetic partitioner should be used.\n"
    "This partitioner is highly experimental and should not be used by mere mortals like you.");
  args.support("no-parmetis", "\nOptional: Specifies that the ParMETIS partitioner should not be used.");
  args.support("no-zoltan", "\nOptional: Specifies that the Zoltan partitioner should not be used.");
  args.support("help", "\nDisplays this help message.");

  // need help?
  if((argc < 2) || (args.check("help") >= 0))
  {
    display_help(comm, args);
    return 0;
  }

  // let's print some basic information
  comm.print(String("MPI processes").pad_back(pad_len, '.') + ": " + stringify(comm.size()));

  // check for unsupported options
  auto unsupported = args.query_unsupported();
  if( !unsupported.empty() )
  {
    // print all unsupported options to cerr
    for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
      comm.print(std::cerr, "\nERROR: unsupported option '--" + (*it).second);

    display_help(comm, args);
    return 1;
  }

  int num_mesh_files = args.check("mesh");
  if(num_mesh_files < 1)
  {
    comm.print(std::cerr, "\nERROR: You have to specify at least one meshfile with --mesh <filenames...>");
    display_help(comm, args);
    return 1;
  }
  if(args.check("parts") < 1)
  {
    comm.print(std::cerr, "\nERROR: You have to specify the number of desired partitions via --num-parts <n>");
    display_help(comm, args);
    return 1;
  }
  if(args.check("out") < 1)
  {
    comm.print(std::cerr, "\nERROR: You have to specify the output partition file via --out <filename>");
    display_help(comm, args);
    return 1;
  }

  // get our filename deque
  auto mpars = args.query("mesh");
  XASSERT(mpars != nullptr);
  const std::deque<String>& filenames = mpars->second;

  comm.print(String("Input Mesh Files").pad_back(pad_len, '.') + ": '" + stringify_join(filenames, "', '") + "'");

  // create an empty mesh file reader
  Geometry::MeshFileReader mesh_reader;
  mesh_reader.add_mesh_files(comm, filenames);

  // read root markup
  mesh_reader.read_root_markup();

  // get mesh type
  const String mtype = mesh_reader.get_meshtype_string();

  comm.print(String("Mesh Type").pad_back(pad_len, '.') + ": " + mtype);

  if(mtype == "conformal:hypercube:2:2")
    return run<H2M2D>(comm, args, mesh_reader);
  if(mtype == "conformal:hypercube:2:3")
    return run<H2M3D>(comm, args, mesh_reader);
  if(mtype == "conformal:hypercube:3:3")
    return run<H3M3D>(comm, args, mesh_reader);
  if(mtype == "conformal:simplex:2:2")
    return run<S2M2D>(comm, args, mesh_reader);
  if(mtype == "conformal:simplex:2:3")
    return run<S2M3D>(comm, args, mesh_reader);
  if(mtype == "conformal:simplex:3:3")
    return run<S3M3D>(comm, args, mesh_reader);

  comm.print(std::cerr, "\nERROR: unsupported mesh type!\n");
  return 1;
}
