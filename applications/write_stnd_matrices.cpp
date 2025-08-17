#include <kernel/base_header.hpp>
#include <kernel/runtime.hpp>

#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/conformal_mesh.hpp>

#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/domain_assembler.hpp>
#include <kernel/assembly/domain_assembler_basic_jobs.hpp>
#include <kernel/assembly/gpdv_assembler.hpp>

#include <kernel/assembly/symbolic_assembler.hpp>

#include <kernel/space/discontinuous/element.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>

#include <kernel/trafo/standard/mapping.hpp>

#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>

#include <kernel/util/simple_arg_parser.hpp>

using namespace FEAT;
typedef double DataType;
typedef Index IndexType;
// Choose a cubature rule
String cubature_name = "auto-degree:5";

namespace MatrixCreator
{

  template<typename Space_>
  void _write_laplace_matrix(const Space_& space, const String& name_base)
  {
    typedef LAFEM::SparseMatrixCSR<DataType, IndexType> MatrixType;
    constexpr int dim = Space_::shape_dim;

    Assembly::DomainAssembler dom_asm(space.get_trafo());
    dom_asm.compile_all_elements();
    MatrixType matrix;
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix, space);
    Assembly::Common::LaplaceOperator laplace;
    matrix.format();
    Assembly::assemble_bilinear_operator_matrix_1(dom_asm, matrix, laplace, space, cubature_name);

    std::cout << "Matrix size " << name_base << " " << (matrix.bytes()/(1024*1024)) << "MB \n";
    matrix.write_out(LAFEM::FileMode::fm_mtx, name_base + "_dim_" + stringify(dim) + ".mtx");
  }

  template<typename Space1_, typename Space2_>
  void _write_saddlepoint_matrices(const Space1_& space1, const Space2_& space2, const String& name_base)
  {
    constexpr int dim = Space1_::shape_dim;
    typedef LAFEM::SparseMatrixBCSR<DataType, IndexType, dim, dim> MatrixTypeA;
    typedef LAFEM::SparseMatrixBCSR<DataType, IndexType, dim, 1> MatrixTypeB;
    typedef LAFEM::SparseMatrixBCSR<DataType, IndexType, 1, dim> MatrixTypeD;

    Assembly::DomainAssembler dom_asm(space1.get_trafo());
    dom_asm.compile_all_elements();
    {
      MatrixTypeA matrix;
      Assembly::SymbolicAssembler::assemble_matrix_std1(matrix, space1);
      Assembly::Common::DuDvOperatorBlocked<dim> defo;
      matrix.format();
      Assembly::assemble_bilinear_operator_matrix_1(dom_asm, matrix, defo, space1, cubature_name);
      std::cout << "Matrix size " << name_base << " " << (matrix.bytes()/(1024*1024)) << "MB \n";

      matrix.write_out(LAFEM::FileMode::fm_mtx, name_base + "_dim_" + stringify(dim) + "_" +
                                                    "block_" + stringify(dim) + "x" + stringify(dim)+ ".mtx");
    }
    {
      MatrixTypeB mat_b;
      MatrixTypeD mat_d;
      Assembly::GradPresDivVeloAssembler::assemble(mat_b, mat_d, space1, space2, cubature_name);


      mat_b.write_out(LAFEM::FileMode::fm_mtx, name_base + "_dim_" + stringify(dim) + "_" +
                                                    "block_" + stringify(dim) + "x" + stringify(1)+ ".mtx");
      mat_d.write_out(LAFEM::FileMode::fm_mtx, name_base + "_dim_" + stringify(dim) + "_" +
                                                    "block_" + stringify(1) + "x" + stringify(dim)+ ".mtx");

    }
  }

  template<typename Shape_>
  void run(SimpleArgParser& args, Geometry::MeshFileReader& mesh_reader)
  {
    typedef Shape_ ShapeType;
    std::cout << "Shape Type: " << ShapeType::name() << "\n";
    // Use the unstructured conformal mesh class
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    // Use the standard transformation mapping
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    // Use the Lagrange-1 element
    typedef Space::Lagrange1::Element<TrafoType> SpaceTypeQ1;
    typedef Space::Lagrange2::Element<TrafoType> SpaceTypeQ2;
    typedef Space::Discontinuous::ElementP1<TrafoType> SpaceTypeDiscP1;

    Geometry::MeshAtlas<MeshType> mesh_atlas;
    typedef Geometry::RootMeshNode<MeshType> RootMeshNodeType;
    std::unique_ptr<RootMeshNodeType> mesh_node;
    std::cout << "Parsing mesh files...\n";
    try
    {
      mesh_node = mesh_reader.parse(mesh_atlas);
    }
    catch(const std::exception& exc)
    {
      std::cerr << "\nERROR: " << exc.what() << "\n";
      Runtime::abort();
    }
    Index level(0);
    if(args.parse("level", level) < 0)
    {
      // That's not meant to happen...
      std::cerr << "\nERROR: Failed to parse '--level' parameter\n";

      // and abort our program
      Runtime::abort();
    }
    if(level > Index(0))
      std::cout << "Refining Mesh to Level " << level << "...\n";

    // Okay, so let's refine the mesh-node up to the desired level:
    for(Index lvl(0); lvl < level; ++lvl)
    {
      // Refining the mesh-node is really easy: we just have to call the 'refine_unique' function
      // of the mesh-node, which gives us a unique_ptr of the refined mesh node.
      mesh_node = mesh_node->refine_unique();
    }
    MeshType& mesh = *mesh_node->get_mesh();
    std::cout << "Number of Elements: " << mesh.get_num_elements() << "\n";
    TrafoType trafo(mesh);

    // mtx file base
    String mtx_base;
    if(args.check("mtx-name-schema") > 0)
    {
      mtx_base = args.query("mtx-name-schema")->second.front();
    }
    else
    {
      mtx_base = "feat_matrix";
    }

    mtx_base += "_lvl" + stringify(level);



    SpaceTypeQ1 space_q1(trafo);
    SpaceTypeQ2 space_q2(trafo);
    SpaceTypeDiscP1 space_disc(trafo);

    std::cout << "Space Q1 DoFs " << space_q1.get_num_dofs() << "\n";
    std::cout << "Space Q2 DoFs " << space_q2.get_num_dofs() << "\n";

    _write_laplace_matrix(space_q1, mtx_base + "_q1_crs");
    _write_laplace_matrix(space_q2, mtx_base + "_q2_crs");
    _write_saddlepoint_matrices(space_q2, space_disc, mtx_base + "_q2_bcrs");


  }


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
    // 3) Optional: The VTK output filename option:
    args.support("mtx-name-schema", "<filename>\n"
      "Optional: Specifies the filename of the output mtx files.\n"
      "If not given, a default name is used.\n");

    // Before checking for unsupported arguments, let's check if the user supplied any arguments
    // at all, as otherwise this application would call Runtime::abort() further below:
    if(args.num_args() <= 1)
    {
      std::cout << "\n";
      std::cout << "Info: For this tutorial, you need to specify at least one\n";
      std::cout << "input mesh file via the '--mesh <filenames...>' option.\n";
      std::cout << "\n";

      std::cout << "Supported Options:\n";
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
        std::cerr << "ERROR: unsupported option #" << (*it).first << " '--" << (*it).second << "'\n";
      }

      // And abort
      Runtime::abort();
    }

    // Now let us check whether all mandatory options have been given.
    // In our case, the only mandatory option is the '--mesh' option, which specifies the path(s)
    // to the mesh file(s) that describe our computational domain.
    if(args.check("mesh") <= 0)
    {
      // This tutorial requires at least one mesh file...
      std::cerr << "\nERROR: You have to specify at least one mesh file via '--mesh <files...>'\n";

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
    // discretize the same analytic domain and you want to "outsource" the common
    // analytic domain description into a separate common file.

    // Therefore, we always have to expect that the user does not supply just one filename,
    // but a set of filenames that we have to pass to the MeshFileReader, so thus we are
    // always dealing with a deque of filename strings instead of a single string.

    // Due to the above check we know that the user supplied the '--mesh' option, so
    // we can query the deque of strings that represent the parameters of the option:
    const std::deque<String>& filenames = args.query("mesh")->second;

    // For convenience, we'll print the filenames to the console by utilizing the
    // 'stringify_join' function that will stringify each item in the deque and
    // concatenate them using a separator string:
    std::cout << "Mesh Files: " << stringify_join(filenames, " ") << "\n";

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
      std::cerr << "\nERROR: " << exc.what() << "\n";
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
      std::cerr << "\nERROR: " << exc.what() << "\n";
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
      std::cout << "\n";
      std::cout << "ERROR: Mesh file(s) did not provide a mesh-type!\n";
      std::cout << "\n";
      std::cout << "Did you supply all required mesh-files?\n";
      Runtime::abort();
    }

    // And we'll print the mesh type string to the console:
    std::cout << "Mesh Type : " << mesh_type << "\n";

    // Now comes the interesting part:
    // The "mesh_type" string specifies the mesh-type and therefore the shape-type
    // of the mesh that is stored in the mesh-file(s). At this point, we check all
    // five mesh-types that are supported by this tutorial application using the
    // following if-else cascade and call the corresponding "run" function template
    // specialization for the required shape-type and pass our SimpleArgParser and
    // MeshFileReader objects as parameters to it:
    if     (mesh_type == "conformal:hypercube:1:1") // 1D mesh
      run<Shape::Hypercube<1>>(args, mesh_reader);
    else if(mesh_type == "conformal:hypercube:2:2") // 2D quadrilateral mesh
      run<Shape::Hypercube<2>>(args, mesh_reader);
    else if(mesh_type == "conformal:hypercube:3:3") // 3D hexahedron mesh
      run<Shape::Hypercube<3>>(args, mesh_reader);
    else if(mesh_type == "conformal:simplex:2:2")   // 2D triangle mesh
      run<Shape::Simplex<2>>(args, mesh_reader);
    else if(mesh_type == "conformal:simplex:3:3")   // 3D tetrahedron mesh
      run<Shape::Simplex<3>>(args, mesh_reader);
    else
    {
      // The mesh-type is either invalid or not supported. In fact, there are some
      // other valid mesh-type specifiers (e.g. sub-dimensional meshes), which this
      // (and most other) application does not support by construction.
      std::cout << "\nERROR: unsupported mesh type!\n";
      Runtime::abort();
    }

    // And that's it for the main function.
  } // void main(...)
}











// Here's our main function
int main(int argc, char* argv[])
{
  // Initialize the runtime
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  // Print a welcome message
  std::cout << "Welcome to FEAT's tutorial #08: MeshFileReader\n";

  // call the tutorial's main function
  MatrixCreator::main(argc, argv);

  // Finalize the runtime
  return 0;
}