// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.


// Misc. FEAT includes
#include <kernel/util/string.hpp>                          // for String
#include <kernel/runtime.hpp>                         // for Runtime
#include <kernel/util/dist.hpp>                            // NEW: for Dist::Comm

// FEAT-Geometry includes
#include <kernel/geometry/conformal_mesh.hpp>              // for ConformalMesh
#include <kernel/geometry/export_vtk.hpp>                  // for ExportVTK
#include <kernel/geometry/mesh_part.hpp>                   // for MeshPart
#include <kernel/geometry/mesh_node.hpp>                   // NEW: for RootMeshNode, MeshNode
#include <kernel/geometry/unit_cube_patch_generator.hpp>   // NEW: for UnitCubePatchGenerator
#include <kernel/geometry/boundary_factory.hpp>

// FEAT-Trafo includes
#include <kernel/trafo/standard/mapping.hpp>               // the standard Trafo mapping

// FEAT-Space includes
#include <kernel/space/lagrange2/element.hpp>              // the Lagrange-1 Element (aka "Q1")
#include <kernel/space/lagrange1/element.hpp>              // the Lagrange-1 Element (aka "Q1")
#include <kernel/space/discontinuous/element.hpp>

// FEAT-Cubature includes
#include <kernel/cubature/dynamic_factory.hpp>             // for DynamicFactory

// FEAT-Analytic includes
#include <kernel/analytic/common.hpp>                      // for SineBubbleFunction
#include <kernel/adjacency/export_tga.hpp>


// FEAT-Assembly includes
#include <kernel/assembly/symbolic_assembler.hpp>          // for SymbolicAssembler
#include <kernel/assembly/unit_filter_assembler.hpp>       // for UnitFilterAssembler
#include <kernel/assembly/error_computer.hpp>              // for L2/H1-error computation
#include <kernel/assembly/bilinear_operator_assembler.hpp> // for BilinearOperatorAssembler
#include <kernel/assembly/linear_functional_assembler.hpp> // for LinearFunctionalAssembler
#include <kernel/assembly/discrete_projector.hpp>          // for DiscreteVertexProjector
#include <kernel/assembly/common_operators.hpp>            // for LaplaceOperator
#include <kernel/assembly/common_functionals.hpp>          // for LaplaceFunctional
#include <kernel/assembly/mirror_assembler.hpp>            // NEW: for MirrorAssembler

// FEAT-LAFEM includes
#include <kernel/lafem/dense_vector.hpp>                   // for DenseVector
#include <kernel/lafem/sparse_matrix_csr.hpp>              // for SparseMatrixCSR
#include <kernel/lafem/mean_filter.hpp>                    // for UnitFilter
#include <kernel/lafem/unit_filter.hpp>                    // for UnitFilter
#include <kernel/lafem/none_filter.hpp>                    // for UnitFilter
#include <kernel/lafem/vector_mirror.hpp>                  // NEW: for VectorMirror
#include <kernel/lafem/tuple_vector.hpp>
#include <kernel/lafem/tuple_filter.hpp>
#include <kernel/lafem/tuple_matrix.hpp>
#include <kernel/lafem/tuple_mirror.hpp>

// FEAT-Global includes
#include <kernel/global/gate.hpp>                          // NEW: for Global::Gate
#include <kernel/global/filter.hpp>                        // NEW: for Global::Filter
#include <kernel/global/matrix.hpp>                        // NEW: for Global::Matrix
#include <kernel/global/vector.hpp>                        // NEW: for Global::Vector
#include <kernel/global/alg_dof_parti_system.hpp>

// FEAT-Solver includes
#include <kernel/solver/pcg.hpp>                           // for PCG
#include <kernel/solver/schwarz_precond.hpp>               // NEW: for SchwarzPrecond
#include <kernel/solver/ilu_precond.hpp>                   // NEW: for ILUPrecond

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// We are using FEAT
using namespace FEAT;


template<typename DT_, typename IT_>
String stringify_buf(const LAFEM::MatrixMirrorBuffer<DT_, IT_>& buf)
{
  String s1;
  s1 += "[";
  const auto* row_ptr = buf.row_ptr();
  const auto* col_idx = buf.col_ind();
  Index n = buf.rows();
  for(Index i(0); i < n; ++i)
  {
    if(i > 0u)
      s1 += " |";
    for(Index j = row_ptr[i]; j < row_ptr[i+1]; ++j)
    {
      s1 += " ";
      s1 += stringify(col_idx[j]);
    }
  }
  s1 += "]";
  return s1;
}

template<typename IT_>
String stringify_struct(Index n, const IT_* row_ptr, const IT_* col_idx)
{
  String s1;
  s1 += "[";
  for(Index i(0); i < n; ++i)
  {
    if(i > 0u)
      s1 += " |";
    for(Index j = row_ptr[i]; j < row_ptr[i+1]; ++j)
    {
      s1 += " ";
      s1 += stringify(col_idx[j]);
    }
  }
  s1 += "]";
  return s1;
}

template<typename DT_, typename IT_>
String stringify_val(Index n, const IT_* row_ptr, const DT_* val)
{
  String s1;
  s1 += "[";
  for(Index i(0); i < n; ++i)
  {
    if(i > 0u)
      s1 += " |";
    for(Index j = row_ptr[i]; j < row_ptr[i+1]; ++j)
    {
      s1 += " ";
      s1 += stringify_fp_fix(val[j], 3);
    }
  }
  s1 += "]";
  return s1;
}

// Here's our main function
int main(int argc, char* argv[])
{
  // Before we can do anything else, we first need to initialize the FEAT runtime environment:
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  // Specify the desired mesh refinement level, defaulted to 5.
  Index level(2);
  int dbg_rank = -1;

  // Now let's see if we have command line parameters: This tutorial supports passing
  // the refinement level as a command line parameter, to investigate the behavior of the L2/H1
  // errors of the discrete solution.
  if(argc > 1)
  {
    // Try to parse the last argument to obtain the desired mesh refinement level.
    int ilevel(0);
    if(!String(argv[1]).parse(ilevel) || (ilevel < 1))
    {
      // Failed to parse
      std::cerr << "ERROR: Failed to parse '" << argv[1] << "' as refinement level.\n";
      std::cerr << "Note: The last argument must be a positive integer.\n";
      // Abort our runtime environment
      Runtime::abort();
    }
    // If parsing was successful, use the given information and notify the user
    level = Index(ilevel);
  }
  if(argc > 2)
  {
    String(argv[2]).parse(dbg_rank);
  }

  typedef Shape::Quadrilateral ShapeType;
  typedef Geometry::ConformalMesh<ShapeType> MeshType;
  typedef Geometry::MeshPart<MeshType> MeshPartType;
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
  typedef Space::Lagrange2::Element<TrafoType> SpaceType1;
  typedef Space::Lagrange1::Element<TrafoType> SpaceType2;
  //typedef Space::Discontinuous::ElementP1<TrafoType> SpaceType2;
  typedef Geometry::RootMeshNode<MeshType> RootMeshNodeType;

  typedef double DataType;
  typedef Index IndexType;

  typedef LAFEM::DenseVectorBlocked<DataType, IndexType, 2> LocalVectorTypeV;
  typedef LAFEM::DenseVector<DataType, IndexType> LocalVectorTypeP;
  typedef LAFEM::TupleVector<LocalVectorTypeV, LocalVectorTypeP> LocalVectorType;

  //typedef LAFEM::DenseVector<DataType, IndexType> BufferVectorType;

  typedef LAFEM::VectorMirror<DataType, IndexType> VectorMirrorTypeV;
  typedef LAFEM::VectorMirror<DataType, IndexType> VectorMirrorTypeP;
  typedef LAFEM::TupleMirror<VectorMirrorTypeV, VectorMirrorTypeP> VectorMirrorType;

  typedef Global::Gate<LocalVectorType, VectorMirrorType> GateType;


  typedef LAFEM::SparseMatrixBCSR<DataType, IndexType, 2, 2> LocalMatrixTypeA;
  typedef LAFEM::SparseMatrixBCSR<DataType, IndexType, 2, 1> LocalMatrixTypeB;
  typedef LAFEM::SparseMatrixBCSR<DataType, IndexType, 1, 2> LocalMatrixTypeD;
  //typedef LAFEM::SparseMatrixCSR<DataType, IndexType> LocalMatrixType11;
  //typedef LAFEM::SparseMatrixCSR<DataType, IndexType> LocalMatrixType12;
  //typedef LAFEM::SparseMatrixCSR<DataType, IndexType> LocalMatrixType21;
  //typedef LAFEM::SparseMatrixCSR<DataType, IndexType> LocalMatrixType22;

  //typedef LAFEM::TupleMatrix<
    //LAFEM::TupleMatrixRow<LocalMatrixType11, LocalMatrixType12>,
    //LAFEM::TupleMatrixRow<LocalMatrixType21, LocalMatrixType22>> LocalMatrixType;
  typedef LAFEM::SaddlePointMatrix<LocalMatrixTypeA, LocalMatrixTypeB, LocalMatrixTypeD> LocalMatrixType;

  typedef LAFEM::UnitFilterBlocked<DataType, IndexType, 2> LocalFilterTypeV;
  //typedef LAFEM::UnitFilter<DataType, IndexType> LocalFilterType1;
  //typedef LAFEM::MeanFilter<DataType, IndexType> LocalFilterType2;
  typedef LAFEM::UnitFilter<DataType, IndexType> LocalFilterTypeP;

  typedef LAFEM::TupleFilter<LocalFilterTypeV, LocalFilterTypeP> LocalFilterType;

  typedef Global::Matrix<LocalMatrixType, VectorMirrorType, VectorMirrorType> GlobalMatrixType;
  typedef Global::Vector<LocalVectorType, VectorMirrorType> GlobalVectorType;
  typedef Global::Filter<LocalFilterType, VectorMirrorType> GlobalFilterType;

  //typedef Global::AlgDofPartiVector<LocalVectorType, VectorMirrorType> AlgDofPartiVectorType;
  //typedef Global::AlgDofPartiMatrix<LocalMatrixType, VectorMirrorType> AlgDofPartiMatrixType;

  Dist::Comm comm = Dist::Comm::world();
  comm.print("\nNumber of processes: " + stringify(comm.size()));
  if((comm.size() != 1) && (comm.size() != 4) && (comm.size() != 16) && (comm.size() != 64))
  {
    comm.print(std::cerr, "ERROR: You must run this tutorial with 1, 4, 16 or 64 processes!");
    Runtime::abort();
  }

  std::unique_ptr<RootMeshNodeType> root_mesh_node;

  std::vector<int> neighbor_ranks;

  Index lvl = Geometry::UnitCubePatchGenerator<MeshType>::create_unique(
    comm.rank(), comm.size(), root_mesh_node, neighbor_ranks);

  {
    Geometry::MaskedBoundaryFactory<MeshType> bnd_factory(*root_mesh_node->get_mesh());
    for(const auto& x : root_mesh_node->get_halo_map())
      bnd_factory.add_mask_meshpart(*x.second);
    bnd_factory.compile();
    root_mesh_node->add_mesh_part("bnd", bnd_factory.make_unique());
  }

  String msg = "Neighbors of process " + stringify(comm.rank()) + ":";
  for(int i : neighbor_ranks)
    msg += " " + stringify(i);
  comm.allprint(msg);

  comm.print("\nBase Mesh Level: " + stringify(lvl));

  if(lvl < level)
  {
    comm.print("Refining Mesh to Level " + stringify(level) + "...");

    for(; lvl < level; ++lvl)
    {
      root_mesh_node = root_mesh_node->refine_unique();
    }
  }

  MeshType& mesh = *root_mesh_node->get_mesh();
  TrafoType trafo(mesh);
  SpaceType1 space_1(trafo);
  SpaceType2 space_2(trafo);


  GateType gate(comm);
  for(auto it = neighbor_ranks.begin(); it != neighbor_ranks.end(); ++it)
  {
    const int neighbor_rank = *it;
    const MeshPartType* neighbor_halo = root_mesh_node->get_halo(neighbor_rank);

    XASSERTM(neighbor_halo != nullptr, "Failed to retrieve neighbor halo!");

    VectorMirrorType neighbor_mirror;
    Assembly::MirrorAssembler::assemble_mirror(neighbor_mirror.at<0>(), space_1, *neighbor_halo);
    Assembly::MirrorAssembler::assemble_mirror(neighbor_mirror.at<1>(), space_2, *neighbor_halo);
    gate.push(neighbor_rank, std::move(neighbor_mirror));
  }

  {
    LocalVectorType tmp;
    tmp.at<0>() = LocalVectorTypeV(space_1.get_num_dofs());
    tmp.at<1>() = LocalVectorTypeP(space_2.get_num_dofs());
    gate.compile(std::move(tmp));
  }

  GlobalMatrixType matrix(&gate, &gate);

  GlobalFilterType filter;

  {
    Assembly::UnitFilterAssembler<MeshType> unit_asm;
    auto* part = root_mesh_node->find_mesh_part("bnd");
    if(part)
      unit_asm.add_mesh_part(*part);
    unit_asm.assemble(filter.local().at<0>(), space_1);
    unit_asm.assemble(filter.local().at<1>(), space_2);
  }

  Assembly::SymbolicAssembler::assemble_matrix_std1(matrix.local().block_a(), space_1);
  Assembly::SymbolicAssembler::assemble_matrix_std2(matrix.local().block_b(), space_1, space_2);
  Assembly::SymbolicAssembler::assemble_matrix_std2(matrix.local().block_d(), space_2, space_1);

  Random rng(std::size_t(17 + 7*comm.rank()));
  matrix.local().block_a().format(rng, 1.0, 9.0);
  matrix.local().block_b().format(rng, 1.0, 9.0);
  matrix.local().block_d().format(rng, 1.0, 9.0);

  //matrix.local().at<0,0>().format(0.11);
  //matrix.local().at<0,1>().format(0.12);
  //matrix.local().at<1,0>().format(0.21);
  //matrix.local().at<1,1>().format(0.22);

  GlobalVectorType vec_x = matrix.create_vector_r();
  GlobalVectorType vec_x2 = matrix.create_vector_r();
  GlobalVectorType vec_b = matrix.create_vector_l();
  GlobalVectorType vec_b2 = matrix.create_vector_l();
  GlobalVectorType vec_d = matrix.create_vector_l();

  vec_x.format(rng, 1.0, 9.0);
  vec_x2.format(0.0);
  vec_b.format(0.0);
  vec_b2.format(0.0);

  filter.filter_cor(vec_x);

  matrix.apply(vec_b, vec_x);

  filter.filter_def(vec_b);

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  comm.print(String(120, '#'));

  comm.print("Assembling AlgDofParti...");

  auto adp = std::make_shared<Global::AlgDofParti<LocalVectorType, VectorMirrorType>>();

  // breakpoint for debugging
#ifdef FEAT_COMPILER_MICROSOFT
  if(comm.rank() == dbg_rank) __debugbreak();
#else
  (void)dbg_rank;
#endif

  adp->assemble_by_gate(gate);
  adp->assemble_allgather(true);

  comm.print("Block Information:");
  comm.allprint(adp->get_block_information());

  String s2;

  comm.print("\nGlobal DOF Offet / Count:");
  comm.allprint(stringify(adp->_global_dof_count) + " / " + stringify(adp->_global_dof_offset)
    + " : " + stringify(adp->get_num_owned_dofs()));

  comm.print("\n> _global_dof_idx");
  comm.allprint(stringify(adp->_global_dof_idx));

  comm.print("\n> _owned_mirror");
  comm.allprint(stringify(adp->_owned_mirror));

  comm.print("\n> _neighbors_owner");
  s2.clear();
  for(Index i = 0; i < adp->get_num_owner_neighbors(); ++i)
    s2 += stringify(adp->get_owner_rank(i)) + ":" + stringify(adp->get_owner_mirror(i)) + " / ";
  comm.allprint(s2);

  comm.print("\n> _neighbors_donee");
  s2.clear();
  for(Index i = 0; i < adp->get_num_donee_neighbors(); ++i)
    s2 += stringify(adp->get_donee_rank(i)) + ":" + stringify(adp->get_donee_mirror(i)) + " / ";
  comm.allprint(s2);

  comm.print_flush();

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  comm.print(String(120, '-'));
  comm.print("Assembling AlgDofPartiSystem...");
  comm.print_flush();

  auto adp_system = std::make_shared<Global::AlgDofPartiSystem<LocalMatrixType, LocalFilterType, VectorMirrorType>>
    (matrix, filter, adp);

  adp_system->init_symbolic();

  typedef std::uint32_t IndexTypeADP;
  typedef double DataTypeADP;

  std::vector<IndexTypeADP> adp_row_ptr(adp_system->get_adp_matrix_rows()+1u);
  std::vector<IndexTypeADP> adp_col_idx(adp_system->get_adp_matrix_nzes());
  std::vector<DataTypeADP> adp_val(adp_system->get_adp_matrix_nzes());

  adp_system->upload_matrix_symbolic(adp_row_ptr.data(), adp_col_idx.data());
  adp_system->upload_matrix_numeric(adp_val.data(), adp_row_ptr.data(), adp_col_idx.data());

  adp_system->upload_filter();

  comm.barrier();

  String s1;

  comm.print("\n> adp_system _owner_bufs:");
  s1.clear();
  for(const auto& x : adp_system->_owner_bufs)
  {
    s1 += stringify_buf(x) + " / ";
  }
  comm.allprint(s1);

  comm.print("\n> adp_system _donee_bufs:");
  s1.clear();
  for(const auto& x : adp_system->_donee_bufs)
  {
    s1 += stringify_buf(x) + " / ";
  }
  comm.allprint(s1);

  comm.print("\n> adp_system _owned_matrix:");
  comm.allprint(stringify_struct(adp_system->get_adp_matrix_rows(), adp_row_ptr.data(), adp_col_idx.data()));

  comm.print("\n> adp_system _donee_data_mirs:");
  s1.clear();
  for(const auto& x : adp_system->_donee_data_mirs)
  {
    s1 += stringify(x) + " / ";
  }
  comm.allprint(s1);

  comm.print("\n> adp_system _owner_data_mirs:");
  s1.clear();
  for(const auto& x : adp_system->_owner_data_mirs)
  {
    s1 += x.dump() + " / ";
  }
  comm.allprint(s1);

  comm.print("\n> adp_system _owned_data_mir:");
  s1.clear();
  comm.allprint(adp_system->_owned_data_mir.dump());

  comm.print("\n> adp_system owned _matrix val:");
  comm.allprint(stringify_val(adp_system->get_adp_matrix_rows(), adp_row_ptr.data(), adp_val.data()));

  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////

  std::vector<DataTypeADP> adp_vec_x2(adp->get_num_owned_dofs());
  std::vector<DataTypeADP> adp_vec_b2(adp->get_num_owned_dofs());

  adp->upload_vector(adp_vec_x2.data(), vec_x.local(), nullptr, nullptr);

  adp_system->filter_vec_cor(adp_vec_x2.data());

  adp->download_vector(adp_vec_x2.data(), vec_x2.local(), nullptr, nullptr);

  vec_d.copy(vec_x);
  vec_d.axpy(vec_x2, -1.0);
  DataType x2_err = vec_d.max_abs_element() / vec_x.max_abs_element();

  comm.print("\n> vector X values:");
  comm.allprint(stringify(vec_x.local()) + " vs\n" + stringify(vec_x2.local()));
  comm.print("> vector X error: " + stringify_fp_sci(x2_err));
  comm.print_flush();

  adp_system->filter_matrix(adp_val.data(), adp_row_ptr.data(), adp_col_idx.data());

  //adp_system->apply(adp_vec_b2, adp_vec_x2);
  adp_system->apply(adp_vec_b2.data(), adp_vec_x2.data(), adp_val.data(), adp_row_ptr.data(), adp_col_idx.data());

  //adp_system->filter_vec_def(adp_vec_b2.data());

  adp->download_vector(adp_vec_b2.data(), vec_b2.local(), nullptr, nullptr);

  vec_d.copy(vec_b);
  vec_d.axpy(vec_b2, -1.0);
  DataType b2_err = vec_d.max_abs_element() / vec_b.max_abs_element();

  comm.print("\n> vector B values:");
  comm.allprint(stringify(vec_b.local()) + " vs\n" + stringify(vec_b2.local()));

  comm.print("> vector B error: " + stringify_fp_sci(b2_err));
  comm.print_flush();

  adp_system.reset();
  adp.reset();

  return 0;
}
