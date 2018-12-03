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


// FEAT-LAFEM includes
#include <kernel/lafem/dense_vector.hpp>                   // for DenseVector
// we need to have cgal to be able to include it
#ifdef FEAT_HAVE_CGAL
#include <kernel/geometry/cgal.hpp>

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// We are using FEAT
using namespace FEAT;

namespace AnalyseMeshCGAL
{
  // Our LAFEM containers work in main memory.
  typedef Mem::Main MemType;
  // Our data arrays should be double precision.
  typedef double DataType;
  // Use the default index type for indexing.
  typedef Index IndexType;

  // Use the standard dense vector
  typedef LAFEM::DenseVector<MemType, DataType, IndexType> VectorType;

  // In most tutorials, we have defined the shape-type of the mesh as a (compile-time
  // constant) typedef at this point prior to defining the actual application code.
  // However, in this application the shape of the mesh is specified by the input mesh file(s) at
  // runtime, so we have a conflict between the code design used in FEAT, which requires the
  // shape-type at compile-time, and the fact that input mesh files are read at runtime.
  //
  // There are two possible solutions to this problem:
  //
  // 1) We could (just as in many tutorials) typedef the shape-type as, say,
  //    Shape::Quadrilateral (or some other type, of course) and then check whether the shape-type
  //    of the input mesh files matches this chosen definition. In case of a mismatch, we would
  //    then simply emit an error message stating the unsupported shape of the input mesh file
  //    and abort the application.
  //
  // 2) We could out-source the actual application code into a separate function template (named e.g.
  //    'run'), which is templatised in the mesh shape-type, and then use an if-else cascade to
  //    call the corresponding shape-type specialisation of that function template based on the
  //    shape-type of the input mesh file(s).
  //
  // In this application, we choose the second solution, which is (slightly) more complex to implement
  // than the first one. For this, we have to split the application code over two functions:
  //
  // 1) The 'main' function, which performs most of the basic initialisation prior to the actual
  //    mesh object construction. In this application, this boils down to setting up the SimpleArgParser
  //    and creating and initialising the MeshFileReader to obtain the shape-type of the input
  //    mesh files.
  //
  // 2) The 'run' function template, which is templatised in the shape-type that is to be used.
  //    This function then contains the 'actual' application code that starts with our well-known
  //    geometry typedefs (MeshType, TrafoType, etc.) as well as the creation of the mesh object.
  //    The remainder of this function is then more or less equivalent to all tutorials.


  // Here comes the forward declaration of the "run" function template; it is templatised in the
  // shape-type, which is chosen from the input mesh file(s), and expects the SimpleArgParser and
  // MeshFileReader objects as arguments, which are created in the "main" function.
  // This 'run' function template is implemented right below the following "main" function.
  template<typename Shape_>
  void run(SimpleArgParser& args, Geometry::MeshFileReader& mesh_reader);


  // Here's our main function
  void main(int argc, char* argv[])
  {
    // First of all, let's create a SimpleArgParser
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
      std::cout << std::endl;
      std::cout << "Info: For this application, you need to specify at least one" << std::endl;
      std::cout << "input mesh file via the '--mesh <filenames...>' option." << std::endl;
      std::cout << "and the off-file via the '--off <filename>' option" <<std::endl;
      std::cout << std::endl;

      std::cout << "Supported Options:" << std::endl;
      std::cout << args.get_supported_help();
      return;
    }

    // As usual, we check for unsupported arguments and abort with an error message if necessary:
    std::deque<std::pair<int,String>> unsupported = args.query_unsupported();
    if(!unsupported.empty())
    {
      // Print unsupported options
      for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
      {
        std::cerr << "ERROR: unsupported option #" << (*it).first << " '--" << (*it).second << "'" << std::endl;
      }

      // And abort
      Runtime::abort();
    }

    // Now let us check whether all mandatory options have been given.
    // In our case, the only mandatory option is the '--mesh' option, which specifies the path(s)
    // to the mesh file(s) that describe our computational domain.
    if(args.check("mesh") <= 0)
    {
      // This tool requires at least one mesh file...
      std::cerr << std::endl << "ERROR: You have to specify at least one mesh file via '--mesh <files...>'" << std::endl;

      // ...so we cannot continue without one.
      Runtime::abort();
    }
    if(args.check("off") <= 0)
    {
      // This application requires at least one off file...
      std::cerr << std::endl << "ERROR: You have to specify one off file via '--off <filename>'" << std::endl;

      // ...so we cannot continue without one.
      Runtime::abort();
    }

    // At this point, you may ask why one should supply more than one mesh file here.
    // The reason for this is that the FEAT mesh file(s) do not only contain a mesh but
    // also possibly other additional important data such as
    // 1) charts for describing the analytic domain boundaries
    // 2) mesh-parts for describing the discrete domain boundary parts
    // 3) partitionings for parallel simulations

    // Now all these parts are usually combined in a single mesh file for the sake of
    // convenience, but it is also possible to split this data among several disjoint
    // files. This is especially interesting if you have several different meshes that
    // discretise the same analytic domain and you want to "outsource" the common
    // analytic domain description into a separate common file.

    // Therefore, we always have to expect that the user does not supply just one filename,
    // but a set of filenames that we have to pass to the MeshFileReader, so thus we are
    // always dealing with a deque of filename strings instead of a single string.

    // Due to the above check we know that the user supplied the '--mesh' option, so
    // we can query the deque of strings that represent the parameters of the option:
    const std::deque<String>& filenames = args.query("mesh")->second;

    // For convenience, we'll print the filenames to the console by utilising the
    // 'stringify_join' function that will stringify each item in the deque and
    // concatenate them using a separator string:
    std::cout << "Mesh Files: " << stringify_join(filenames, " ") << std::endl;

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
      std::cerr << std::endl << "ERROR: " << exc.what() << std::endl;
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
      std::cerr << std::endl << "ERROR: " << exc.what() << std::endl;
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
      std::cout << std::endl;
      std::cout << "ERROR: Mesh file(s) did not provide a mesh-type!" << std::endl;
      std::cout << std::endl;
      std::cout << "Did you supply all required mesh-files?" << std::endl;
      Runtime::abort();
    }

    // And we'll print the mesh type string to the console:
    std::cout << "Mesh Type : " << mesh_type << std::endl;

    // Now comes the interesting part:
    // The "mesh_type" string specifies the mesh-type and therefore the shape-type
    // of the mesh that is stored in the mesh-file(s). At this point, we check all
    // five mesh-types that are supported by this application using the
    // following if-else cascade and call the corresponding "run" function template
    // specialisation for the required shape-type and pass our SimpleArgParser and
    // MeshFileReader objects as parameters to it:
    // Note that everything else than the 3d-meshes is commented out
    // because the cgal-wrapper supports only 3d-situations at the time of coding.
//    if     (mesh_type == "conformal:hypercube:1:1") // 1D mesh
//      run<Shape::Hypercube<1>>(args, mesh_reader);
//    else if(mesh_type == "conformal:hypercube:2:2") // 2D quadrilateral mesh
//      run<Shape::Hypercube<2>>(args, mesh_reader);
    if(mesh_type == "conformal:hypercube:3:3") // 3D hexahedron mesh
      run<Shape::Hypercube<3>>(args, mesh_reader);
//    else if(mesh_type == "conformal:simplex:2:2")   // 2D triangle mesh
//      run<Shape::Simplex<2>>(args, mesh_reader);
    else if(mesh_type == "conformal:simplex:3:3")   // 3D tetrahedron mesh
      run<Shape::Simplex<3>>(args, mesh_reader);
    else
    {
      // The mesh-type is either invalid or not supported. In fact, there are some
      // other valid mesh-type specifiers (e.g. sub-dimensional meshes), which this
      // (and most other) application does not support by construction.
      std::cout << std::endl << "ERROR: unsupported mesh type!" << std::endl;
      Runtime::abort();
    }

    // And that's it for the main function.
  } // void main(...)


  // And this is where the magic happens: our 'run' function template.
  template<typename Shape_>
  void run(SimpleArgParser& args, Geometry::MeshFileReader& mesh_reader)
  {
    // First, define our shape-type based on the template parameter. At this point,
    // 'ShapeType' may refer to any of those types that were specified in the if-else
    // cascade in the main function above.
    typedef Shape_ ShapeType;

    // Let's print out the shape-type name just for convenience:
    std::cout << "Shape Type: " << ShapeType::name() << std::endl;

    // Now that we know the shape-type that we want to use in this specialisation of the
    // 'run' function template, we can continue with the remaining typedefs as usual:

    // Use the unstructured conformal mesh class
    typedef Geometry::ConformalMesh<ShapeType> MeshType;

    // At this point, we have to introduce a new class:
    // the (mesh) atlas, which is merely a class that manages
    // a named set of charts. A chart is a geometric object that is used for the analytic
    // description of the computational domain, such as e.g. a circle, a Bezier-Spline,
    // a sphere or a surface triangulation. At this point, we could provide a typedef
    // for the mesh atlas type, but we can also skip this, as we require the actual
    // type only for the one following variable declaration.
    // The mesh atlas class is templatised in the mesh type and we do not have to pass
    // anything to its constructor:
    Geometry::MeshAtlas<MeshType> mesh_atlas;

    // In addition, we have to define the (root) mesh-node type (which has also been used in
    // tutorial 06). This class is responsible for managing a mesh as well as a named set of
    // mesh parts defined on that particular mesh. As we are going to need this type several
    // times, we'll use a typedef for convenience here:
    typedef Geometry::RootMeshNode<MeshType> RootMeshNodeType;

    // Next, let us create a shared-pointer for the mesh-node, which will be assigned to
    // a new object by the mesh file reader in the next step.
    std::shared_ptr<RootMeshNodeType> mesh_node;

    // Okay, let the parsing begin:
    std::cout << "Parsing mesh files..." << std::endl;

    // At this point, a lot of errors could occur if either one of the files is corrupt/invalid
    // or if the set of different input files is inconsistent/mismatching. Therefore, we enclose
    // the actual call of the 'parse' function in a try-catch block:
    try
    {
      // In its most simple form (precisely: overload), the only input parameter is a reference
      // to the mesh atlas and the function returns a shared-pointer to a newly allocated
      // mesh-node object, which we assign to the previously declared variable:
      mesh_node = mesh_reader.parse(mesh_atlas);
    }
    catch(const std::exception& exc)
    {
      // That's not good...
      std::cerr << std::endl << "ERROR: " << exc.what() << std::endl;
      Runtime::abort();
    }

    // If we come out here, then the user has actually succeeded in choosing a consistent, valid
    // and supported set of input mesh files and the mesh file reader has successfully finished
    // its job. As a first action with our new mesh-node, we'll query the names of all mesh-parts
    // from the mesh-node, so that we can print them to the console.
    // The 'true' parameter for the following function call specifies that we want to query only
    // the names of the mesh-parts that were parsed from the mesh-file and that we are NOT
    // interested in any additional mesh-parts that the management code had to create for internal
    // use (mesh-parts for "internal use" are characterised by a leading underscore in their name).
    std::deque<String> meshpart_names = mesh_node->get_mesh_part_names(true);
    std::cout << "Mesh Parts: " << stringify_join(meshpart_names, " ") << std::endl;

    // In many cases, the mesh-file only contains a relatively "coarse" mesh that first has to be
    // refined a few times to obtain a mesh that is fine enough for a finite element discretisation.
    // The main reason for this is that FEAT is a software package that uses geometric multigrid as
    // one of its core components, and geometric multigrid needs a mesh hierarchy.
    // However, we do not want to deal with multigrid in this application, but we probably still need
    // to refine the parsed mesh, so let us check whether the user has specified a desired level
    // and, if so, try to parse it:
    Index level(0);
    if(args.parse("level", level) < 0)
    {
      // That's not meant to happen...
      std::cerr << std::endl << "ERROR: Failed to parse '--level' parameter" << std::endl;

      // and abort our program
      Runtime::abort();
    }

    // Did the user specify a level > 0?
    if(level > Index(0))
      std::cout << "Refining Mesh to Level " << level << "..." << std::endl;

    // Okay, so let's refine the mesh-node up to the desired level:
    for(Index lvl(0); lvl < level; ++lvl)
    {
      // Refining the mesh-node is really easy: we just have to call the 'refine' function
      // of the mesh-node. Since we are using shared-pointers here, we have to encapsulate
      // the 'naked' pointer returned by the function in a shared_ptr:
      mesh_node = std::shared_ptr<RootMeshNodeType>(mesh_node->refine());
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    // Okay, now that we have refined the mesh-node up to the desired level, we can now obtain a
    // reference to the actual mesh from it by calling the 'get_mesh' function and dereferencing
    // the returned pointer:
    MeshType& mesh = *mesh_node->get_mesh();

    // We'll print out the number of mesh elements just for kicks:
    Index number_of_elements = mesh.get_num_elements();
    Index number_of_vertices = mesh.get_num_vertices();
    std::cout << "Number of Elements: " << number_of_elements << std::endl;
    std::cout << "Number of Vertices: " << number_of_vertices << std::endl;

    // Vectors to save the result
    VectorType distance_to_off(number_of_vertices);
    VectorType squared_distance(number_of_vertices);
    VectorType inside_outside(number_of_vertices);

    DataType* dist_to_off = distance_to_off.elements();
    DataType* sq_dist = squared_distance.elements();
    DataType* in_out = inside_outside.elements();

    // load the off-file with the cgal-wrapper
    std::cout<<"Loading the off-file"<<std::endl;
    String off_file_name;
    args.parse("off",off_file_name);
    Geometry::CGALWrapper cw(off_file_name);

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
      in_out[ivtx] = sign_distance;

      sq_distance = cw.squared_distance(vertices[ivtx][0],vertices[ivtx][1],vertices[ivtx][2]);
      sq_dist[ivtx] = sq_distance;

      dist_to_off[ivtx] = sign_distance * Math::sqrt(sq_distance);
    }



    // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    // Post-Processing: Export to VTK file

    // The filename of our VTK file
    String vtk_name;
    if(args.parse("vtk", vtk_name) > 0)
    {
      std::cout << std::endl;
      std::cout << "Writing VTK file '" << vtk_name << ".vtu'..." << std::endl;

      // Create a VTK exporter for our mesh
      Geometry::ExportVTK<MeshType> exporter(mesh);

      // add the solution vector
      exporter.add_vertex_scalar("InsideOutside", inside_outside.elements());
      exporter.add_vertex_scalar("SquaredDistance", squared_distance.elements());
      exporter.add_vertex_scalar("DistanceToOff", distance_to_off.elements());

      //write the VTK file
      exporter.write(vtk_name);
    }

    // That's it for today.
    std::cout << std::endl << "Finished!" << std::endl;
  } // void run<Shape_>(...)

} // namespace

// Here's our main function
int main(int argc, char* argv[])
{
  // Initialise the runtime
  Runtime::initialise(argc, argv);

  // call the tutorial's main function
  AnalyseMeshCGAL::main(argc, argv);

  // Finalise the runtime
  return Runtime::finalise();
}
#else
int main(/*int argc, char* argv[]*/)
{
  // Print a welcome message
  std::cerr << "Not Compiled with CGAL. Cannot run without CGAL" << std::endl;
  return 1;
}
#endif  // FEAT_HAVE_CGAL