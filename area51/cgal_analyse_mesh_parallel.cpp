// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

//
// Supported Command Line Parameters
// =================================
// This application supports the following command line parameters:
//
// --mesh <filenames...>
// Mandatory: Specifies the filename(s) of at least one input mesh file.
//
// --level <n>
// Optional: Specifies the mesh refinement level.
// If not given, level 0 (unrefined) is used.
//
// --vtk <name>
// Optional: Specifies the filename of the output VTK file.
// If not given, no VTK file output is written.
//
// --off <name>
// Mandatory: Specify the off-file that is to be used for the analysis.
//
// \author Malte Schuh
//

// We start our application with a batch of includes...

#include <kernel/base_header.hpp>
// Misc. FEAT includes
#include <kernel/util/string.hpp>                          // for String
#include <kernel/util/runtime.hpp>                         // for Runtime
#include <kernel/util/simple_arg_parser.hpp>               // for SimpleArgParser

// FEAT-Geometry includes
#include <kernel/geometry/conformal_mesh.hpp>              // for ConformalMesh
#include <kernel/geometry/export_vtk.hpp>                  // for ExportVTK
#include <kernel/geometry/mesh_node.hpp>                   // for RootMeshNode
#include <kernel/geometry/mesh_file_reader.hpp>            // for MeshFileReader

#include <kernel/util/statistics.hpp>
#include <kernel/util/stop_watch.hpp>
#include <control/domain/parti_domain_control.hpp>
#include <control/scalar_basic.hpp>
#include <control/statistics.hpp>

#include <kernel/util/dist_file_io.hpp>                   // to read the mesh in parallel

// not exactly neccesary because we will not do any calculations,
// but the domain control requires a mapping and a space
#include <kernel/trafo/standard/mapping.hpp>
// We use a hack somewhere: For lagrange1, the number of DOF
// and the number of vertices in the mesh is identical!
#include <kernel/space/lagrange1/element.hpp>

// We need to have cgal to be able to include it
#ifdef FEAT_HAVE_CGAL
#include <kernel/geometry/cgal.hpp>


namespace AnalyseMeshCGALParallel
{
  using namespace FEAT;

  // define our arch types
  typedef Mem::Main MemType;
  typedef Real DataType;
  typedef Index IndexType;

  template<typename Shape_>
  void run(SimpleArgParser& args, Dist::Comm& comm)
  {
    // First, define our shape-type based on the template parameter. At this point,
    // 'ShapeType' may refer to any of those types that were specified in the if-else
    // cascade in the main function above.
    typedef Shape_ ShapeType;

    // Let's print out the shape-type name just for convenience:
    comm.print("Shape Type: " + ShapeType::name());

    // Now that we know the shape-type that we want to use in this specialization of the
    // 'run' function template, we can continue with the remaining typedefs as usual:

    // Use the unstructured conformal mesh class
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    // since we do not solve any pde we do not need the
    // trafo-type or the space-type, but the domain-control
    // needs it.
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    // We use a hack later: The number of vertices in the
    // mesh is identical to the dof for the lagrange.
    typedef Space::Lagrange1::Element<TrafoType> SpaceType;

    // create our domain control
    typedef Control::Domain::SimpleDomainLevel<MeshType, TrafoType, SpaceType> DomainLevelType;
    typedef Control::Domain::PartiDomainControl<DomainLevelType> DomainType;
    DomainType domain(comm, true);

    domain.parse_args(args);
    domain.set_desired_levels(args.query("level")->second);
    domain.create(args.query("mesh")->second);
    int lvl_max = domain.get_desired_level_max();
    comm.print("Max Level: " + stringify(lvl_max));

    // print partitioning info
    comm.print(domain.get_chosen_parti_info());

    // define our system level
    typedef Control::ScalarUnitFilterSystemLevel<MemType, DataType, IndexType> SystemLevelType;

    std::deque<std::shared_ptr<SystemLevelType>> system_levels;

    const Index num_levels = domain.size_physical();

    // create system levels
    for (Index i(0); i < num_levels; ++i)
    {
      system_levels.push_back(std::make_shared<SystemLevelType>());
    }

    // We are mostly interested in the gates!
    comm.print("Assembling gates, muxers and transfers...");

    for (Index i(0); i < num_levels; ++i)
    {
      system_levels.at(i)->assemble_gate(domain.at(i));
      if((i+1) < domain.size_virtual())
      {
        system_levels.at(i)->assemble_coarse_muxer(domain.at(i+1));
        //system_levels.at(i)->assemble_transfer(domain.at(i), domain.at(i+1), cubature);
      }
    }

    // Fetch the root-mesh-node from the finest level and the corresponding mesh
    typedef typename DomainType::MeshNodeType RootMeshNodeType;
    // therefore, first fetch the finest level
    DomainLevelType& finest_lvl = *domain.front();
    RootMeshNodeType* rmn = finest_lvl.get_mesh_node();
    MeshType& mesh = *rmn->get_mesh();

    // Now it is somehow identical to the non-parallel:
    // We'll print out the number of mesh elements for the rank 0 just for kicks:
    Index number_of_elements = mesh.get_num_elements();
    Index number_of_vertices = mesh.get_num_vertices();
    if (comm.rank() == 0)
    {
      comm.print("Number of Elements for rank 0: " + stringify(number_of_elements));
      comm.print("Number of Vertices for rank 0: " + stringify(number_of_vertices));
    }

    // Now create the global and local vector
    // Get the local and global vector types
    typedef typename SystemLevelType::LocalSystemVector LocalVectorType;
    typedef typename SystemLevelType::GlobalSystemVector GlobalVectorType;

    // Now the system level gates are assembled, so we can use them to
    // create the global vectors
    // We create them by feet
    GlobalVectorType inside_outside_global(&(system_levels.front()->gate_sys),number_of_vertices);
    GlobalVectorType squared_distance_global(&(system_levels.front()->gate_sys),number_of_vertices);
    GlobalVectorType signed_distance_global(&(system_levels.front()->gate_sys),number_of_vertices);

    // Get the local vectors
    LocalVectorType& inside_outside_local = inside_outside_global.local();
    LocalVectorType& squared_distance_local = squared_distance_global.local();
    LocalVectorType& signed_distance_local = signed_distance_global.local();

    DataType* SignedDist = signed_distance_local.elements();
    DataType* SqDist = squared_distance_local.elements();
    DataType* InOut = inside_outside_local.elements();

    // Load the off-filename
    String off_file_name;
    args.parse("off",off_file_name);
    StopWatch loading_off_file;
    loading_off_file.reset();
    loading_off_file.start();
    // We could now do
    // Geometry::CGALWrapper cw(off_file_name);
    // But in the parallel case this is not good.
    // Better way is to use the dist_file_io to read it into a stream
    // and hand the stream to the constructor:
    std::stringstream ss_off_file;
    DistFileIO::read_common(ss_off_file,off_file_name,comm);
    Geometry::CGALWrapper cw(ss_off_file);
    loading_off_file.stop();
    comm.print("Time loading the Off-file: " + loading_off_file.elapsed_string(TimeFormat::s_m));

    // Get the Vertex-Set of the mesh
    typedef typename MeshType::VertexSetType VertexSetType;
    VertexSetType& vertices = (&mesh)->get_vertex_set();

    // Temp variables as result for the inside-outside-test
    bool is_inside;
    DataType sign_distance;
    DataType sq_distance;

    // Loop through all vertices and analyse if they are
    // inside or outside of the geometrical object defined by the
    // off-file and calculate the distance to it.
    StopWatch distance_calc;
    distance_calc.reset();
    distance_calc.start();
    for(Index ivtx(0); ivtx<number_of_vertices;++ivtx)
    {
      is_inside = cw.point_inside(vertices[ivtx][0],vertices[ivtx][1],vertices[ivtx][2]);
      if (is_inside == true)
      {
        sign_distance = -1.0;
      }
      else
      {
        sign_distance = 1.0;
      }
      InOut[ivtx] = sign_distance;

      sq_distance = cw.squared_distance(vertices[ivtx][0],vertices[ivtx][1],vertices[ivtx][2]);
      SqDist[ivtx] = sq_distance;

      SignedDist[ivtx] = sign_distance * Math::sqrt(sq_distance);
    }
    distance_calc.stop();
    comm.print("Time spent with distance calculation: " + distance_calc.elapsed_string(TimeFormat::s_m));

    // Syncronise
    // sync_0: Add up, sync_1: average value at the boarders
    signed_distance_global.sync_1();
    squared_distance_global.sync_1();
    inside_outside_global.sync_1();

    String vtk_name;
    String lvl;
    args.parse("level",lvl);
    if(args.parse("vtk", vtk_name) > 0)
    {
      vtk_name = vtk_name + "-lvl-"+ stringify(lvl_max) + String("-n-") + stringify(comm.size());
      comm.print("Writing vtk-file: " + vtk_name);
      Geometry::ExportVTK<MeshType> exporter(mesh);
      exporter.add_vertex_scalar("InsideOutside",inside_outside_local.elements());
      exporter.add_vertex_scalar("SquaredDistance",squared_distance_local.elements());
      exporter.add_vertex_scalar("DistanceToOff",signed_distance_local.elements());
      exporter.write(vtk_name,comm);
    }
  }

  void main(int argc, char* argv [])
  {
    // create world communicator
    Dist::Comm comm(Dist::Comm::world());

    // print number of processes
    comm.print("Number of Processes: " + stringify(comm.size()));

    // create arg parser
    SimpleArgParser args(argc, argv);
    // Now let's add the supported options
    // For the sake of simplicity, the actual arguments are restricted to:
    // 1) Mandatory: The input mesh option:
    args.support("mesh", "<filenames...>\n"
      "Mandatory: Specifies the filenames of the mesh files that are to be parsed.\n");
    // 2) Optional: The refinement level option:
    args.support("level", "<n>\n"
      "Optional: Specifies the refinement level. If not given, defaults to 0.\n");

    args.support("off", "<filename>\n"
     "Mandatory: Specifies the off-file that is to be analysed\n");

    // 4) Optional: The VTK output filename option:
    args.support("vtk", "<filename>\n"
      "Optional: Specifies the filename of the output VTU file.\n"
      "If not given, no output VTU file is created.\n");

    // Before checking for unsupported arguments, let's check if the user supplied any arguments
    // at all, as otherwise this application would call Runtime::abort() further below:
    if(args.num_args() <= 2)
    {
      comm.print("");
      comm.print("Info: For this application, you need to specify at least one");
      comm.print("input mesh file via the '--mesh <filenames...>' option.");
      comm.print("and the off-file via the '--off <filename>' option");

      comm.print("Supported Options:");
      comm.print(args.get_supported_help());
      return;
    }

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

    // Now let us check whether all mandatory options have been given.
    // In our case, the only mandatory option is the '--mesh' option, which specifies the path(s)
    // to the mesh file(s) that describe our computational domain.
    if(args.check("mesh") <= 0)
    {
      // This application requires at least one mesh file...
      comm.print(std::cerr,"ERROR: You have to specify at least one mesh file via '--mesh <files...>'");

      // ...so we cannot continue without one.
      Runtime::abort();
    }
    if(args.check("off") <= 0)
    {
      // This application requires at least one off file...
      comm.print(std::cerr,"ERROR: You have to specify one off file via '--off <filename>'");

      // ...so we cannot continue without one.
      Runtime::abort();
    }

    // Now check what kind of mesh we have
    // Due to the above check we know that the user supplied the '--mesh' option, so
    // we can query the deque of strings that represent the parameters of the option:
    const std::deque<String>& filenames = args.query("mesh")->second;

    // For convenience, we'll print the filenames to the console by utilising the
    // 'stringify_join' function that will stringify each item in the deque and
    // concatenate them using a separator string:
    comm.print("Mesh Files: " + stringify_join(filenames, " "));

    // Now we can create the actual reader object. Note that we need only one mesh reader
    // object even if there are multiple files to be read.
    Geometry::MeshFileReader mesh_reader;

    // Next we have to pass the filename deque to the reader's 'add_mesh_files' function.
    // However, if one of the files does not exist, the following function may throw
    // an instance of the 'FEAT::FileNotFound' exception, so we enclose this function
    // call in a try-catch block:
    try
    {
      // try to add the filenames to the mesh reader:
      mesh_reader.add_mesh_files(filenames);
    }
    catch(const std::exception& exc)
    {
      // Something went wrong; probably one of the files could not be opened...
      comm.print(std::cerr, "ERROR: " + stringify(exc.what()));
      Runtime::abort();
    }

    // Okay, if we arrive there, then the reader has successfully opened all
    // input files. However, the reader did not actually start parsing any of
    // those files yet, so we will initiate the actual parsing process by
    // calling the 'read_root_markup' function. Again, this function may
    // throw various types of exceptions -- most notably Xml::Scanner related
    // exceptions that derive from Xml::Error -- so we also put this call
    // inside a try-catch block:
    try
    {
      // Try to read the root markup: at this point the reader analyses the first
      // line(s) of the input files to determine whether the files are actually
      // FEAT mesh files. If not, then the parser will throw an instance of the
      // 'Xml::SyntaxError' (in case a file is not an XML file) or 'Xml::GrammarError'
      // (in case the file is an XML file but not a FEAT mesh file).
      mesh_reader.read_root_markup();
    }
    catch(const std::exception& exc)
    {
      // Oops...
      comm.print(std::cerr , "ERROR: " + stringify(exc.what()));
      Runtime::abort();
    }

    // If we arrive here, then the reader has successfully started parsing
    // the root nodes of all input files and is now able to tell us what
    // type of mesh is stored in the input file(s) (assuming that the mesh
    // files provide this information in the root markup, which is not
    // mandatory). So we can now query the mesh type from the reader:
    const String mesh_type = mesh_reader.get_meshtype_string();

    // Let us first check whether the mesh-file(s) provide us with a mesh-type.
    // If this information is missing, then the user probably did not specify all
    // required files -- or the mesh file simply did not specify the mesh-type.
    if(mesh_type.empty())
    {
      comm.print("ERROR: Mesh file(s) did not provide a mesh-type!");
      comm.print("");
      comm.print("Did you supply all required mesh-files?");
      Runtime::abort();
    }

    // And we'll print the mesh type string to the console:
    comm.print("Mesh Type : " + mesh_type);

    // Run the according run-function that has the mesh-type as a template
    // parameter
    // At the Moment we only allow 3d because cgal allows only 3d, but
    // we leave this here so if 2d or 1d is supported from cgal then we
    // can just comment it in.
    //    if     (mesh_type == "conformal:hypercube:1:1") // 1D mesh
//      run<Shape::Hypercube<1>>(args,comm);
//    else if(mesh_type == "conformal:hypercube:2:2") // 2D quadrilateral mesh
//      run<Shape::Hypercube<2>>(args,comm);
    if(mesh_type == "conformal:hypercube:3:3") // 3D hexahedron mesh
      run<Shape::Hypercube<3>>(args,comm);
//    else if(mesh_type == "conformal:simplex:2:2")   // 2D triangle mesh
//      run<Shape::Simplex<2>>(args, mesh_reader);
    else if(mesh_type == "conformal:simplex:3:3")   // 3D tetrahedron mesh
      run<Shape::Simplex<3>>(args,comm);
    else
    {
      // The mesh-type is either invalid or not supported. In fact, there are some
      // other valid mesh-type specifiers (e.g. sub-dimensional meshes), which this
      // (and most other) application does not support by construction.
      comm.print(std::cerr,"ERROR: unsupported mesh type!");
      Runtime::abort();
    }
  }
} //AnalyseMeshCGALParallel

int main(int argc, char* argv [])
{
  FEAT::Runtime::initialize(argc, argv);
  try
  {
    AnalyseMeshCGALParallel::main(argc, argv);
  }
  catch (const std::exception& exc)
  {
    std::cerr << "ERROR: unhandled exception: " << exc.what() << std::endl;
    FEAT::Runtime::abort();
  }
  catch (...)
  {
    std::cerr << "ERROR: unknown exception" << std::endl;
    FEAT::Runtime::abort();
  }
  return FEAT::Runtime::finalize();
}

#else // not defined(FEAT_HAVE_CGAL)

int main(/*int argc, char* argv []*/)
{
  std::cerr << "ERROR: Not compiled with CGAL!" << std::endl;
  return 1;
}

#endif // FEAT_HAVE_CGAL
