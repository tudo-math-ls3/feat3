// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/mesh_permutation.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/assembly/grid_transfer.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/runtime.hpp>

using namespace FEAT;
using namespace FEAT::Geometry;

typedef LAFEM::SparseMatrixCSR<Mem::Main, double, Index> MatrixType;
typedef LAFEM::DenseVector<Mem::Main, double, Index> VectorType;

String get_file_title(const String& filename)
{
  // find last slash
  std::size_t p = filename.find_last_of("\\/");
  if(p == filename.npos)
    p = 0;
  else
    ++p;

  // fine last dot
  std::size_t q = filename.find_last_of(".");
  if(q == filename.npos)
    return filename.substr(p);
  else
    return filename.substr(p, q-p);
}

template<typename Mesh_>
void process_coloring(Geometry::ExportVTK<Mesh_>& exp, const Mesh_& mesh)
{
  const std::vector<Index>& col = mesh.get_mesh_permutation().get_element_coloring();
  if(col.empty())
    return;

  std::vector<double> v(mesh.get_num_elements(), 0.0);
  for(std::size_t i(0); i+1u < col.size(); ++i)
  {
    Index k0 = col.at(i);
    Index k1 = col.at(i+1u);
    for(Index k(k0); k < k1; ++k)
      v[k] = double(i+1u);
  }

  exp.add_cell_scalar("coloring", v.data());
}

template<typename Mesh_>
void process_layering(Geometry::ExportVTK<Mesh_>& exp, const Mesh_& mesh)
{
  const std::vector<Index>& lay = mesh.get_mesh_permutation().get_element_layering();
  if(lay.empty())
    return;

  std::vector<double> v(mesh.get_num_elements(), 0.0), w(mesh.get_num_elements(), 0.0);
  for(std::size_t i(0); i+1u < lay.size(); ++i)
  {
    Index k0 = lay.at(i);
    Index k1 = lay.at(i+1u);
    for(Index k(k0); k < k1; ++k)
    {
      v[k] = double(i+1u);
      w[k] = double(i%2);
    }
  }

  exp.add_cell_scalar("layering", v.data());
  exp.add_cell_scalar("layer_mod2", w.data());
}

template<typename Mesh_>
void prolongate(VectorType& vec_f1, VectorType& vec_f2, Mesh_& mesh_f, Mesh_& mesh_c)
{
  typedef Trafo::Standard::Mapping<Mesh_> TrafoType;
  typedef Space::Lagrange1::Element<TrafoType> SpaceType;

  TrafoType trafo_f(mesh_f);
  TrafoType trafo_c(mesh_c);
  SpaceType space_f(trafo_f);
  SpaceType space_c(trafo_c);

  vec_f1 = VectorType(space_f.get_num_dofs());
  vec_f2 = VectorType(space_f.get_num_dofs());
  VectorType vec_c = VectorType(space_c.get_num_dofs());

  // interpolate on fine mesh
  Analytic::Common::SineBubbleFunction<Mesh_::world_dim> sine_bubble;
  Assembly::Interpolator::project(vec_f1, sine_bubble, space_f);

  // assemble prolongation
  MatrixType prol;
  Assembly::SymbolicAssembler::assemble_matrix_2lvl(prol, space_f, space_c);

  prol.format();
  Cubature::DynamicFactory cubature("auto-degree:3");
  Assembly::GridTransfer::assemble_prolongation_direct(prol, space_f, space_c, cubature);

  // assemble truncation
  MatrixType trunc = prol.transpose();
  trunc.format();
  Assembly::GridTransfer::assemble_truncation_direct(trunc, space_f, space_c, cubature);

  // truncate f1 to c
  trunc.apply(vec_c, vec_f1);

  // prolongate to f2
  prol.apply(vec_f2, vec_c);
}

template<typename Mesh_>
void run_xml(SimpleArgParser& args, Geometry::MeshFileReader& mesh_reader, const String& filename)
{
  // parse levels
  Index lvl_min(0);
  Index lvl_max(1);
  args.parse("level", lvl_max, lvl_min);

  // parse permutation strategy
  Geometry::PermutationStrategy strat = Geometry::PermutationStrategy::none;
  String sperm;
  args.parse("perm", sperm);
  if(sperm.empty() || (sperm == "none"))
    strat = Geometry::PermutationStrategy::none;
  else if((sperm == "rand") || (sperm == "random"))
    strat = Geometry::PermutationStrategy::random;
  else if((sperm == "lexi"))
    strat = Geometry::PermutationStrategy::lexicographic;
  else if((sperm == "color") || (sperm == "colored"))
    strat = Geometry::PermutationStrategy::colored;
  else if((sperm == "cmk") || (sperm == "acmk"))
    strat = Geometry::PermutationStrategy::cuthill_mckee;
  else if((sperm == "rcmk") || (sperm == "racmk"))
    strat = Geometry::PermutationStrategy::cuthill_mckee_reversed;
  else if((sperm == "gcmk"))
    strat = Geometry::PermutationStrategy::geometric_cuthill_mckee;
  else if((sperm == "rgcmk"))
    strat = Geometry::PermutationStrategy::geometric_cuthill_mckee_reversed;
  else
  {
    std::cout << "ERROR: unknown permutation strategy '" << sperm << "'" << std::endl;
    return;
  }

  std::cout << "Permutation strategy: ";
  switch(strat)
  {
  case Geometry::PermutationStrategy::none:
    std::cout << "none" << std::endl;
    break;
  case Geometry::PermutationStrategy::other:
    std::cout << "other" << std::endl;
    break;
  case Geometry::PermutationStrategy::random:
    std::cout << "random" << std::endl;
    break;
  case Geometry::PermutationStrategy::colored:
    std::cout << "colored" << std::endl;
    break;
  case Geometry::PermutationStrategy::lexicographic:
    std::cout << "lexicographic" << std::endl;
    break;
  case Geometry::PermutationStrategy::cuthill_mckee:
    std::cout << "algebraic Cuthill-McKee" << std::endl;
    break;
  case Geometry::PermutationStrategy::cuthill_mckee_reversed:
    std::cout << "algebraic Cuthill-McKee reversed" << std::endl;
    break;
  case Geometry::PermutationStrategy::geometric_cuthill_mckee:
    std::cout << "geometric Cuthill-McKee" << std::endl;
    break;
  case Geometry::PermutationStrategy::geometric_cuthill_mckee_reversed:
    std::cout << "geometric Cuthill-McKee reversed" << std::endl;
    break;
  default:
    // make picky compilers STFU
    break;
  }


  // create an empty atlas and a root mesh node
  Geometry::MeshAtlas<Mesh_> atlas;
  //Geometry::RootMeshNode<Mesh_>* node = new Geometry::RootMeshNode<Mesh_>(nullptr, atlas);
  auto node = std::make_shared<Geometry::RootMeshNode<Mesh_>>(nullptr, &atlas);

  // try to parse the mesh file
#ifndef DEBUG
  try
#endif
  {
    std::cout << "Parsing mesh files..." << std::endl;
    // Now parse the mesh file
    mesh_reader.parse(*node, atlas, nullptr);
  }
#ifndef DEBUG
  catch(std::exception& exc)
  {
    std::cerr << "ERROR: " << exc.what() << std::endl;
    return;
  }
  catch(...)
  {
    std::cerr << "ERROR: unknown exception" << std::endl;
    return;
  }
#endif

  // adapt coarse mesh
  node->adapt();

  Random rng;

  // refine
  for(Index lvl(1); lvl <= lvl_max; ++lvl)
  {
    std::cout << "Refining up to level " << lvl << "..." << std::endl;
    auto coarse = node;
    node = coarse->refine_shared();

    if(lvl < lvl_min)
    {
      continue;
    }

    //auto* fine = coarse->refine();
    auto fine = node->clone_shared();

    // permute coarse and fine mesh nodes
    coarse->create_permutation(strat);
    fine->create_permutation(strat);

    // get our mesh
    Mesh_& mesh = *fine->get_mesh();

    // validate coloring/layering
    if(!mesh.validate_element_coloring())
      std::cout << "WARNING: invalid coloring!" << std::endl;
    if(!mesh.validate_element_layering())
      std::cout << "WARNING: invalid layering!" << std::endl;

    VectorType vec_f1, vec_f2;
    prolongate(vec_f1, vec_f2, *fine->get_mesh(), *coarse->get_mesh());

    // Create a VTK exporter for our mesh
    FEAT::String vtkname = filename + "." + stringify(lvl);

    std::cout << "Writing file '" << vtkname << ".vtu'..." << std::endl;
    Geometry::ExportVTK<Mesh_> exporter(mesh);

    process_coloring(exporter, mesh);
    process_layering(exporter, mesh);
    exporter.add_vertex_scalar("sine", vec_f1.elements());
    exporter.add_vertex_scalar("prol", vec_f2.elements());

    exporter.write(vtkname);
  }
}

void run(int argc, char* argv[])
{
  // This is the list of all supported meshes that could appear in the mesh file
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 2, Real> S2M2D;
  typedef Geometry::ConformalMesh<Shape::Simplex<3>, 3, Real> S3M3D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, Real> H2M2D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<3>, 3, Real> H3M3D;

  SimpleArgParser args(argc, argv);

  args.support("mesh");
  args.support("vtk");
  args.support("level");
  args.support("perm");

  // check for unsupported options
  auto unsupported = args.query_unsupported();
  if( !unsupported.empty() )
  {
    // print all unsupported options to cerr
    for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
      std::cerr << "ERROR: unsupported option '--" << (*it).second << "'" << std::endl;
    return;
  }

  int num_mesh_files = args.check("mesh");
  if(num_mesh_files < 1)
  {
    std::cerr << "ERROR: You have to specify at least one meshfile with --mesh <files...>" << std::endl;
    return;
  }

  // get our filename deque
  auto mpars = args.query("mesh");
  XASSERT(mpars != nullptr);
  const std::deque<String>& filenames = mpars->second;

  // our VTK output filename
  String vtk_name = get_file_title(filenames.back());
  args.parse("vtk", vtk_name);

  // create an empty mesh file reader
  Geometry::MeshFileReader mesh_reader;
  mesh_reader.add_mesh_files(filenames);

  // read root markup
  mesh_reader.read_root_markup();

  // get mesh type
  const String mtype = mesh_reader.get_meshtype_string();

  std::cout << "Mesh Type: " << mtype << std::endl;

  if(mtype == "conformal:hypercube:2:2") run_xml<H2M2D>(args, mesh_reader, vtk_name); else
  if(mtype == "conformal:hypercube:3:3") run_xml<H3M3D>(args, mesh_reader, vtk_name); else
  if(mtype == "conformal:simplex:2:2") run_xml<S2M2D>(args, mesh_reader, vtk_name); else
  if(mtype == "conformal:simplex:3:3") run_xml<S3M3D>(args, mesh_reader, vtk_name); else
  std::cout << "ERROR: unsupported mesh type!" << std::endl;
}

int main(int argc, char* argv[])
{
  Runtime::initialize(argc, argv);
  run(argc, argv);
  return Runtime::finalize();
}
