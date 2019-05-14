// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/util/runtime.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/assembly/burgers_assembler.hpp>
#include <kernel/assembly/burgers_assembly_job.hpp>
#include <kernel/assembly/domain_assembler.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/analytic/wrappers.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>

namespace DbgBurgersAsm
{
  using namespace FEAT;

  static constexpr int dim = 2;

  typedef Shape::Hypercube<dim> ShapeType;
  typedef Geometry::ConformalMesh<ShapeType> MeshType;
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
  typedef Space::Lagrange2::Element<TrafoType> SpaceType;

  typedef Mem::Main MemType;
  typedef double DataType;
  typedef Index IndexType;

  typedef LAFEM::SparseMatrixCSR<MemType, DataType, IndexType> ScalarMatrixType;
  typedef LAFEM::DenseVector<MemType, DataType, IndexType> ScalarVectorType;

  typedef LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType, dim, dim> BlockedMatrixType;
  typedef LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, dim> BlockedVectorType;

  DataType diff_max(const BlockedMatrixType& a, const BlockedMatrixType& b)
  {
    const auto* val_a = a.val();
    const auto* val_b = b.val();

    DataType r = DataType(0);
    for(Index k(0); k < a.used_elements(); ++k)
    {
      auto x = (val_a[k] - val_b[k]);
      for(int i(0); i < a.BlockHeight; ++i)
        for(int j(0); j < a.BlockWidth; ++j)
          r = Math::max(r, Math::abs(x(i,j)));
    }
    return r;
  }

  DataType diff_max(const BlockedVectorType& a, const BlockedVectorType& b)
  {
    const auto* val_a = a.elements();
    const auto* val_b = b.elements();

    DataType r = DataType(0);
    for(Index k(0); k < a.used_elements(); ++k)
    {
      auto x = (val_a[k] - val_b[k]);
      for(int i(0); i < a.BlockSize; ++i)
        r = Math::max(r, Math::abs(x[i]));
    }
    return r;
  }

  DataType diff_max(const ScalarMatrixType& a, const ScalarMatrixType& b)
  {
    const auto* val_a = a.val();
    const auto* val_b = b.val();

    DataType r = DataType(0);
    for(Index k(0); k < a.used_elements(); ++k)
    {
      r = Math::max(r, Math::abs(val_a[k] - val_b[k]));
    }
    return r;
  }

  DataType diff_max(const ScalarVectorType& a, const ScalarVectorType& b)
  {
    const auto* val_a = a.elements();
    const auto* val_b = b.elements();

    DataType r = DataType(0);
    for(Index k(0); k < a.used_elements(); ++k)
    {
      r = Math::max(r, Math::abs(val_a[k] - val_b[k]));
    }
    return r;
  }

  void main()
  {
    Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(3);
    MeshType mesh(mesh_factory);
    TrafoType trafo(mesh);
    SpaceType space(trafo);

    // create domain assembler
    Assembly::DomainAssembler<TrafoType> domain_asm(trafo);
    domain_asm.compile_all_elements(0u, Assembly::ThreadingStrategy::single);

    BlockedMatrixType bmatrix_1, bmatrix_2;
    ScalarMatrixType smatrix_1, smatrix_2;

    BlockedVectorType bvec_rhs_1, bvec_rhs_2, bvec_sol, bvec_conv;
    ScalarVectorType svec_rhs_1, svec_rhs_2, svec_sol;

    // assemble matrix structures
    Assembly::SymbolicAssembler::assemble_matrix_std1(bmatrix_1, space);
    Assembly::SymbolicAssembler::assemble_matrix_std1(smatrix_1, space);
    bmatrix_2 = bmatrix_1.clone(LAFEM::CloneMode::Layout);
    smatrix_2 = smatrix_1.clone(LAFEM::CloneMode::Layout);

    // interpolate some convection vector
    {
      Analytic::Common::CosineWaveFunction<dim> cos_func;
      //Analytic::ScalarCurl<decltype(cos_func)> conv_func(cos_func);
      Analytic::Gradient<decltype(cos_func)> conv_func(cos_func);
      Assembly::Interpolator::project(bvec_conv, conv_func, space);
      Assembly::Interpolator::project(bvec_sol, conv_func, space);
      Assembly::Interpolator::project(svec_sol, cos_func, space);
    }

    // create vectors
    bvec_rhs_1 = bmatrix_1.create_vector_l();
    bvec_rhs_2 = bmatrix_1.create_vector_l();
    svec_rhs_1 = smatrix_1.create_vector_l();
    svec_rhs_2 = smatrix_1.create_vector_l();

    // format everything to be assembled
    bmatrix_1.format();
    bmatrix_2.format();
    smatrix_1.format();
    smatrix_2.format();
    bvec_rhs_1.format();
    bvec_rhs_2.format();
    svec_rhs_1.format();
    svec_rhs_2.format();

    const String cubature_name = "gauss-legendre:3";
    Cubature::DynamicFactory cubature_factory(cubature_name);

    // create old burgers assembler
    Assembly::BurgersAssembler<DataType, IndexType, dim> burgers_asm;
    burgers_asm.deformation = true;
    burgers_asm.nu = burgers_asm.theta = burgers_asm.beta = burgers_asm.frechet_beta = 1.0;
    burgers_asm.sd_delta = burgers_asm.sd_nu = 1.0;
    burgers_asm.set_sd_v_norm(bvec_conv);

    // assemble matrix and RHS vector
    burgers_asm.assemble_matrix(bmatrix_1, bvec_conv, space, cubature_factory, 1.0);
    burgers_asm.assemble_vector(bvec_rhs_1, bvec_conv, bvec_sol, space, cubature_factory, 1.0);

    // create burgers jobs and execute them
    Assembly::BurgersBlockedMatrixAssemblyJob<BlockedMatrixType, SpaceType> bbm_job(bmatrix_2, bvec_conv, space, cubature_name);
    Assembly::BurgersBlockedVectorAssemblyJob<BlockedVectorType, SpaceType> bbv_job(bvec_rhs_2, bvec_sol, bvec_conv, space, cubature_name);

    bbm_job.deformation = burgers_asm.deformation;
    bbm_job.nu = bbm_job.theta = bbm_job.beta = burgers_asm.nu;
    bbm_job.frechet_beta = burgers_asm.frechet_beta;
    bbm_job.sd_delta = bbm_job.sd_nu = burgers_asm.sd_delta;
    bbm_job.set_sd_v_norm(bvec_conv);

    bbv_job.deformation = burgers_asm.deformation;
    bbv_job.nu = bbv_job.theta = bbv_job.beta = burgers_asm.nu;
    bbv_job.frechet_beta = 0.0 /*burgers_asm.frechet_beta*/;
    bbv_job.sd_delta = bbv_job.sd_nu = 0.0 /*burgers_asm.sd_delta*/;
    bbv_job.set_sd_v_norm(bvec_conv);

    domain_asm.assemble(bbm_job);
    domain_asm.assemble(bbv_job);

    std::cout << "Diff of Blocked Matrix: " << stringify_fp_sci(diff_max(bmatrix_1, bmatrix_2)) << std::endl;
    std::cout << "Diff of Blocked Vector: " << stringify_fp_sci(diff_max(bvec_rhs_1, bvec_rhs_2)) << std::endl;

    burgers_asm.deformation = false;
    burgers_asm.nu = burgers_asm.theta = burgers_asm.beta = 1.0;
    burgers_asm.frechet_beta = 0.0;
    burgers_asm.sd_delta = burgers_asm.sd_nu = 1.0;
    burgers_asm.set_sd_v_norm(bvec_conv);

    // assemble matrix and RHS vector
    burgers_asm.assemble_scalar_matrix(smatrix_1, bvec_conv, space, cubature_factory, 1.0);
    //burgers_asm.assemble_scalar_vector(svec_rhs_1, bvec_conv, svec_sol, space, cubature_factory, 1.0);
    smatrix_1.apply(svec_rhs_1, svec_sol);

    // create burgers jobs and execute them
    Assembly::BurgersScalarMatrixAssemblyJob<ScalarMatrixType, SpaceType, BlockedVectorType> bsm_job(smatrix_2, bvec_conv, space, cubature_name);
    Assembly::BurgersScalarVectorAssemblyJob<ScalarVectorType, SpaceType, BlockedVectorType> bsv_job(svec_rhs_2, svec_sol, bvec_conv, space, cubature_name);

    bsm_job.deformation = burgers_asm.deformation;
    bsm_job.nu = bsm_job.theta = bsm_job.beta = burgers_asm.nu;
    bsm_job.frechet_beta = burgers_asm.frechet_beta;
    bsm_job.sd_delta = bsm_job.sd_nu = burgers_asm.sd_delta;
    bsm_job.set_sd_v_norm(bvec_conv);

    bsv_job.deformation = burgers_asm.deformation;
    bsv_job.nu = bsv_job.theta = bsv_job.beta = burgers_asm.nu;
    bsv_job.frechet_beta = burgers_asm.frechet_beta;
    bsv_job.sd_delta = bsv_job.sd_nu = burgers_asm.sd_delta;
    bsv_job.set_sd_v_norm(bvec_conv);

    domain_asm.assemble(bsm_job);
    domain_asm.assemble(bsv_job);

    std::cout << "Diff of Scalar  Matrix: " << stringify_fp_sci(diff_max(smatrix_1, smatrix_2)) << std::endl;
    std::cout << "Diff of Scalar  Vector: " << stringify_fp_sci(diff_max(svec_rhs_1, svec_rhs_2)) << std::endl;
  }
}

int main(int argc, char** argv)
{
  FEAT::Runtime::initialize(argc, argv);
  DbgBurgersAsm::main();
  return FEAT::Runtime::finalize();
}
