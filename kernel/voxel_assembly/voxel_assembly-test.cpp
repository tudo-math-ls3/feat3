// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/runtime.hpp>
#include <kernel/backend.hpp>

#include <test_system/test_system.hpp>

#include <kernel/voxel_assembly/poisson_assembler.hpp>
#include <kernel/voxel_assembly/defo_assembler.hpp>
#include <kernel/voxel_assembly/burgers_assembler.hpp>
#include <kernel/voxel_assembly/burgers_velo_material_assembler.hpp>
#include <kernel/voxel_assembly/helper/voxel_coloring.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/burgers_assembler.hpp>
#include <kernel/assembly/burgers_velo_material_assembly_job.hpp>

#include <kernel/geometry/boundary_factory.hpp>            // for BoundaryFactory
#include <kernel/geometry/conformal_mesh.hpp>              // for ConformalMesh
#include <kernel/geometry/common_factories.hpp>            // for RefinedUnitCubeFactory
#include <kernel/geometry/mesh_part.hpp>                   // for MeshPart

#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange2/element.hpp>

#include <kernel/assembly/symbolic_assembler.hpp>          // for SymbolicAssembler
#include <kernel/assembly/domain_assembler.hpp>            // for DomainAssembler
#include <kernel/assembly/domain_assembler_basic_jobs.hpp>    // for Assembly::assemble_***

#include <kernel/cubature/dynamic_factory.hpp>             // for DynamicFactory

#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>

#include <kernel/adjacency/adjactor.hpp>



using namespace FEAT;
using namespace FEAT::TestSystem;

/**
 * \brief Base Test class for the voxel assembly classes.
 *
 * \test test description missing
 *
 * \tparam DT_
 * description missing
 *
 * \tparam IT_
 * description missing
 *
 * \author Maximilian Esser
 */
template<
  typename DT_,
  typename IT_,
  typename Shape_>
class VoxelAssemblyTest
  : public UnitTest
{
public:
  typedef DT_ DataType;
  typedef IT_ IndexType;
  typedef Shape_ ShapeType;
  typedef Geometry::ConformalMesh<ShapeType> MeshType;
  typedef Geometry::MeshPart<MeshType> MeshPartType;
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
  typedef Space::Lagrange2::Element<TrafoType> Lagrange2SpaceType;
  const int level;

  VoxelAssemblyTest(String test_name, int level_, PreferredBackend backend)
    : UnitTest(test_name, Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend), level(level_)
  {
  }

  virtual ~VoxelAssemblyTest()
  {
  }


  const DataType tol = Math::Limits<DataType>::epsilon() * ( (this->_preferred_backend == PreferredBackend::generic) ? DataType(1000) : DataType(100));

  virtual void run_matrix_compare() const = 0;

  virtual void run_defect_compare() const = 0;

  virtual void run() const override
  {
    this->run_matrix_compare();
    this->run_defect_compare();

  }

  template<typename MatrixType_>
  void compare_matrices(const MatrixType_& mat_a, const MatrixType_& mat_b) const
  {
    auto matrix_comp = mat_a.clone(LAFEM::CloneMode::Deep);
    matrix_comp.copy(mat_b);
    matrix_comp.axpy(mat_a, DataType(-1));
    // std::cout << matrix_comp << "\n";
    TEST_CHECK_EQUAL_WITHIN_EPS(matrix_comp.norm_frobenius(), DataType(0), tol);
  }

  template<typename VectorType_>
  void compare_defects(const VectorType_& vec_a, const VectorType_& vec_b) const
  {
    auto vec_comp = vec_a.clone(LAFEM::CloneMode::Deep);
    vec_comp.copy(vec_b);
    vec_comp.axpy(vec_a, DataType(-1));
    // std::cout << matrix_comp << "\n";
    TEST_CHECK_EQUAL_WITHIN_EPS(vec_comp.norm2(), DataType(0), tol*DataType(1.5));
  }

}; // class VoxelAssemblyTest


/**
 * \brief Poisson Voxel Assembly Test class.
 *
 * \test Test poisson assembly for given shapetype
 *
 * \tparam DT_
 * description missing
 *
 * \tparam IT_
 * description missing
 *
 * \tparam Shape_
 * The used shape
 *
 * \author Maximilian Esser
 */
template<
  typename DT_,
  typename IT_,
  typename Shape_>
class VoxelPoissonAssemblyTest
  : public VoxelAssemblyTest<DT_, IT_, Shape_>
{
public:
  typedef VoxelAssemblyTest<DT_, IT_, Shape_> BaseClass;
  // typedef typename BaseClass::DataType DataType;
  // typedef typename BaseClass::IndexType IndexType;
  typedef typename BaseClass::ShapeType ShapeType;
  typedef typename BaseClass::MeshType MeshType;
  typedef typename BaseClass::MeshPartType MeshPartType;
  typedef typename BaseClass::TrafoType TrafoType;
  typedef typename BaseClass::Lagrange2SpaceType Lagrange2SpaceType;

  VoxelPoissonAssemblyTest(int level_, PreferredBackend backend)
    : BaseClass("VoxelPoissonAssemblyTest" + String("dim") + stringify(Shape_::dimension), level_, backend)
  {
  }

  virtual ~VoxelPoissonAssemblyTest()
  {
  }

  typedef DT_ DataType;
  typedef IT_ IndexType;


  virtual void run_matrix_compare() const override
  {
    typedef LAFEM::SparseMatrixCSR<DataType, IndexType> MatrixType;

    Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(Index(this->level));
    MeshType mesh(mesh_factory);

    Geometry::BoundaryFactory<MeshType> boundary_factory(mesh);
    MeshPartType boundary(boundary_factory);
    TrafoType trafo(mesh);
    Lagrange2SpaceType space(trafo);

    MatrixType matrix_ref;
    MatrixType matrix_new;

    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix_ref, space);
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix_new, space);

    matrix_ref.format();
    matrix_new.format();

    Assembly::DomainAssembler<TrafoType> domain_assembler(trafo);
    domain_assembler.compile_all_elements();

    String cubature_name = "auto-degree:5";

    Assembly::Common::LaplaceOperator laplace_operator;
    Assembly::assemble_bilinear_operator_matrix_1(
      domain_assembler, matrix_ref, laplace_operator, space, cubature_name);

    std::vector<int> coloring = VoxelAssembly::UnitCubeColoring<ShapeType>::create_coloring(this->level);
    VoxelAssembly::test_coloring(mesh, coloring);
    VoxelAssembly::VoxelPoissonAssembler<Lagrange2SpaceType, DataType, IndexType> poisson_assembler(space, coloring, Lagrange2SpaceType::world_dim == 3 ? 8 : 4);

    poisson_assembler.assemble_matrix1(matrix_new, space, Cubature::DynamicFactory(cubature_name));

    this->compare_matrices(matrix_new, matrix_ref);

  }

  virtual void run_defect_compare() const override
  {
    return;
  }


}; // class PoissonVoxelAssemblyTest

/**
 * \brief Voxel Defo Assembly Test class.
 *
 * \test Test poisson assembly for given shapetype
 *
 * \tparam DT_
 * description missing
 *
 * \tparam IT_
 * description missing
 *
 * \tparam Shape_
 * The used shape
 *
 * \author Maximilian Esser
 */
template<
  typename DT_,
  typename IT_,
  typename Shape_>
class VoxelDefoAssemblyTest
  : public VoxelAssemblyTest<DT_, IT_, Shape_>
{
public:
  typedef VoxelAssemblyTest<DT_, IT_, Shape_> BaseClass;
  // typedef typename BaseClass::DataType DataType;
  // typedef typename BaseClass::IndexType IndexType;
  typedef typename BaseClass::ShapeType ShapeType;
  typedef typename BaseClass::MeshType MeshType;
  typedef typename BaseClass::MeshPartType MeshPartType;
  typedef typename BaseClass::TrafoType TrafoType;
  typedef typename BaseClass::Lagrange2SpaceType Lagrange2SpaceType;

  VoxelDefoAssemblyTest(int level_, PreferredBackend backend)
    : BaseClass("VoxelDefoAssemblyTest" + String("dim") + stringify(Shape_::dimension), level_, backend)
  {
  }

  virtual ~VoxelDefoAssemblyTest()
  {
  }

  typedef DT_ DataType;
  typedef IT_ IndexType;


  virtual void run_matrix_compare() const override
  {
    typedef LAFEM::SparseMatrixBCSR<DataType, IndexType, ShapeType::dimension, ShapeType::dimension> MatrixType;
    typedef LAFEM::DenseVectorBlocked<DataType, IndexType, ShapeType::dimension> VectorType;

    Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(Index(this->level));
    MeshType mesh(mesh_factory);

    Geometry::BoundaryFactory<MeshType> boundary_factory(mesh);
    MeshPartType boundary(boundary_factory);
    TrafoType trafo(mesh);
    Lagrange2SpaceType space(trafo);

    MatrixType matrix_ref;
    MatrixType matrix_new;

    VectorType dummy_conv;

    DataType nu = DataType(0.78);

    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix_ref, space);
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix_new, space);

    matrix_ref.format();
    matrix_new.format();

    dummy_conv = matrix_ref.create_vector_r();

    Assembly::BurgersAssembler<DataType, IndexType, ShapeType::dimension> burgers_asm;
    burgers_asm.deformation = true;
    burgers_asm.nu = nu;

    String cubature_name = "auto-degree:5";

    burgers_asm.assemble_matrix(matrix_ref, dummy_conv, space, Cubature::DynamicFactory(cubature_name));

    std::vector<int> coloring = VoxelAssembly::UnitCubeColoring<ShapeType>::create_coloring(this->level);
    VoxelAssembly::test_coloring(mesh, coloring);
    VoxelAssembly::VoxelDefoAssembler<Lagrange2SpaceType, DataType, IndexType> defo_assembler(space, coloring, Lagrange2SpaceType::world_dim == 3 ? 8 : 4);
    defo_assembler.nu = nu;

    defo_assembler.assemble_matrix1(matrix_new, space, Cubature::DynamicFactory(cubature_name));

    this->compare_matrices(matrix_new, matrix_ref);

  }

  virtual void run_defect_compare() const override
  {
    return;
  }


}; // class VoxelDefoAssemblyTest

/**
 * \brief Voxel Burgers Assembly Test class.
 *
 * \test Test burgers assembly for given shapetype
 *
 * \tparam DT_
 * description missing
 *
 * \tparam IT_
 * description missing
 *
 * \tparam Shape_
 * The used shape
 *
 * \author Maximilian Esser
 */
template<
  typename DT_,
  typename IT_,
  typename Shape_>
class VoxelBurgersAssemblyTest
  : public VoxelAssemblyTest<DT_, IT_, Shape_>
{
public:
  typedef VoxelAssemblyTest<DT_, IT_, Shape_> BaseClass;
  // typedef typename BaseClass::DataType DataType;
  // typedef typename BaseClass::IndexType IndexType;
  typedef typename BaseClass::ShapeType ShapeType;
  typedef typename BaseClass::MeshType MeshType;
  typedef typename BaseClass::MeshPartType MeshPartType;
  typedef typename BaseClass::TrafoType TrafoType;
  typedef typename BaseClass::Lagrange2SpaceType Lagrange2SpaceType;

  VoxelBurgersAssemblyTest(int level_, PreferredBackend backend)
    : BaseClass("VoxelBurgersAssemblyTest" + String("dim") + stringify(Shape_::dimension), level_, backend)
  {
  }

  virtual ~VoxelBurgersAssemblyTest()
  {
  }

  typedef DT_ DataType;
  typedef IT_ IndexType;


  virtual void run_matrix_compare() const override
  {
    typedef LAFEM::SparseMatrixBCSR<DataType, IndexType, ShapeType::dimension, ShapeType::dimension> MatrixType;
    typedef LAFEM::DenseVectorBlocked<DataType, IndexType, ShapeType::dimension> VectorType;

    Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(Index(this->level));
    MeshType mesh(mesh_factory);

    Geometry::BoundaryFactory<MeshType> boundary_factory(mesh);
    MeshPartType boundary(boundary_factory);
    TrafoType trafo(mesh);
    Lagrange2SpaceType space(trafo);

    MatrixType matrix_ref;
    MatrixType matrix_new;

    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix_ref, space);
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix_new, space);

    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";

    VectorType dummy_conv = VectorType(rng, matrix_ref.rows(), DataType(-1), DataType(1));

    DataType nu = DataType(0.78);
    DataType beta = DataType(0.3);
    DataType frechet_beta = DataType(1.);
    DataType theta = DataType(1.3);
    DataType sd_delta = DataType(0.57);
    bool deformation = true;
    DataType alpha = DataType(0.66);


    matrix_ref.format();
    matrix_new.format();

    Assembly::BurgersAssembler<DataType, IndexType, ShapeType::dimension> burgers_asm;
    burgers_asm.deformation = deformation;
    burgers_asm.nu = nu;
    burgers_asm.sd_nu = nu;
    burgers_asm.beta = beta;
    burgers_asm.frechet_beta = frechet_beta;
    burgers_asm.theta = theta;
    burgers_asm.sd_delta = sd_delta;
    burgers_asm.set_sd_v_norm(dummy_conv);

    String cubature_name = "auto-degree:5";

    burgers_asm.assemble_matrix(matrix_ref, dummy_conv, space, Cubature::DynamicFactory(cubature_name), alpha);

    std::vector<int> coloring = VoxelAssembly::UnitCubeColoring<ShapeType>::create_coloring(this->level);
    VoxelAssembly::test_coloring(mesh, coloring);
    VoxelAssembly::VoxelBurgersAssembler<Lagrange2SpaceType, DataType, IndexType> burgers_assembler(space, coloring);
    burgers_assembler.deformation = deformation;
    burgers_assembler.nu = nu;
    burgers_assembler.sd_nu = nu;
    burgers_assembler.beta = beta;
    burgers_assembler.frechet_beta = frechet_beta;
    burgers_assembler.theta = theta;
    burgers_assembler.sd_delta = sd_delta;
    burgers_assembler.set_sd_v_norm(dummy_conv);

    TEST_CHECK_EQUAL_WITHIN_EPS(burgers_asm.sd_v_norm, burgers_assembler.sd_v_norm, Math::Limits<DataType>::epsilon()*DataType(100));

    burgers_assembler.assemble_matrix1(matrix_new, dummy_conv, space, Cubature::DynamicFactory(cubature_name), alpha);

    // std::cout << "Matrix ref " << matrix_ref << "\nMatrix_new " << matrix_new << "\n";

    this->compare_matrices(matrix_new, matrix_ref);

  }

  virtual void run_defect_compare() const override
  {
    typedef LAFEM::DenseVectorBlocked<DataType, IndexType, ShapeType::dimension> VectorType;

    Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(Index(this->level));
    MeshType mesh(mesh_factory);

    Geometry::BoundaryFactory<MeshType> boundary_factory(mesh);
    MeshPartType boundary(boundary_factory);
    TrafoType trafo(mesh);
    Lagrange2SpaceType space(trafo);

    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";

    VectorType dummy_conv = VectorType(rng, space.get_num_dofs(), DataType(-1), DataType(1));
    VectorType vector_ref = dummy_conv.clone(LAFEM::CloneMode::Deep);
    VectorType vector_new = dummy_conv.clone(LAFEM::CloneMode::Deep);
    VectorType primal_vec = VectorType(rng, space.get_num_dofs(), DataType(-1), DataType(1));

    DataType nu = DataType(-0.78);
    DataType beta = DataType(0.93);
    DataType theta = DataType(1.3);
    bool deformation = true;
    DataType alpha = DataType(0.66);


    vector_ref.format();
    vector_new.format();

    Assembly::BurgersAssembler<DataType, IndexType, ShapeType::dimension> burgers_asm;
    burgers_asm.deformation = deformation;
    burgers_asm.nu = nu;
    burgers_asm.beta = beta;
    burgers_asm.theta = theta;

    String cubature_name = "auto-degree:5";

    burgers_asm.assemble_vector(vector_ref, dummy_conv, primal_vec, space, Cubature::DynamicFactory(cubature_name), alpha);

    std::vector<int> coloring = VoxelAssembly::UnitCubeColoring<ShapeType>::create_coloring(this->level);
    VoxelAssembly::test_coloring(mesh, coloring);
    VoxelAssembly::VoxelBurgersAssembler<Lagrange2SpaceType, DataType, IndexType> burgers_assembler(space, coloring);
    burgers_assembler.deformation = deformation;
    burgers_assembler.nu = nu;
    burgers_assembler.beta = beta;
    burgers_assembler.theta = theta;

    burgers_assembler.assemble_vector(vector_new, dummy_conv, primal_vec, space, Cubature::DynamicFactory(cubature_name), alpha);

    this->compare_defects(vector_new, vector_ref);

  }


}; // class VoxelDefoAssemblyTest

/**
 * \brief Voxel Burgers Assembly Test class.
 *
 * \test Test burgers assembly for given shapetype
 *
 * \tparam DT_
 * description missing
 *
 * \tparam IT_
 * description missing
 *
 * \tparam Shape_
 * The used shape
 *
 * \author Maximilian Esser
 */
template<
  typename DT_,
  typename IT_,
  typename Shape_,
  VoxelAssembly::MaterialType material_type_>
class VoxelBurgersVeloMaterialAssemblyTest
  : public VoxelAssemblyTest<DT_, IT_, Shape_>
{
public:
  typedef VoxelAssemblyTest<DT_, IT_, Shape_> BaseClass;
  // typedef typename BaseClass::DataType DataType;
  // typedef typename BaseClass::IndexType IndexType;
  typedef typename BaseClass::ShapeType ShapeType;
  typedef typename BaseClass::MeshType MeshType;
  typedef typename BaseClass::MeshPartType MeshPartType;
  typedef typename BaseClass::TrafoType TrafoType;
  typedef typename BaseClass::Lagrange2SpaceType Lagrange2SpaceType;

  VoxelBurgersVeloMaterialAssemblyTest(int level_, PreferredBackend backend)
    : BaseClass("VoxelBurgersVeloMaterialAssemblyTest" + String("dim") + stringify(Shape_::dimension), level_, backend)
  {
  }

  virtual ~VoxelBurgersVeloMaterialAssemblyTest()
  {
  }

  typedef DT_ DataType;
  typedef IT_ IndexType;


  virtual void run_matrix_compare() const override
  {
    typedef LAFEM::SparseMatrixBCSR<DataType, IndexType, ShapeType::dimension, ShapeType::dimension> MatrixType;
    typedef LAFEM::DenseVectorBlocked<DataType, IndexType, ShapeType::dimension> VectorType;

    Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(Index(this->level));
    MeshType mesh(mesh_factory);

    Geometry::BoundaryFactory<MeshType> boundary_factory(mesh);
    MeshPartType boundary(boundary_factory);
    TrafoType trafo(mesh);
    Lagrange2SpaceType space(trafo);

    MatrixType matrix_ref;
    MatrixType matrix_new;

    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix_ref, space);
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix_new, space);

    Random::SeedType seed(Random::SeedType(time(nullptr)));
    Random rng(seed);

    VectorType dummy_conv = VectorType(rng, matrix_ref.rows(), DataType(-1), DataType(1));

    DataType nu = DataType(0.78);
    DataType beta = DataType(0.3);
    DataType frechet_beta = DataType(1.);
    DataType theta = DataType(1.3);
    DataType sd_delta = DataType(0.57);
    bool deformation = true;
    DataType reg_eps = DataType(1E-100);
    DataType frechet_material = DataType(0.7);
    DataType alpha = DataType(0.66);
    DataType mu_0 = nu;
    DataType exp = DataType(1.7);
    DataType lambda = DataType(0.8);
    DataType a_T = DataType(1.2);
    DataType yasuda_a = DataType(1.8);
    DataType mu_inf = DataType(3.5);

    String cubature_name = "auto-degree:5";

    matrix_ref.format();
    matrix_new.format();
    Assembly::DomainAssembler<TrafoType> dom_asm(trafo);
    Assembly::BurgersVeloMaterialBlockedMatrixAssemblyJob burgers_asm(matrix_ref, dummy_conv, space, cubature_name,
      VoxelAssembly::Intern::ViscoFunctorNew<DataType, material_type_>{mu_0, exp, lambda},
      VoxelAssembly::Intern::ViscoDFunctorNew<DataType, material_type_>{mu_0, exp, lambda});

    burgers_asm.deformation = deformation;
    burgers_asm.nu = nu;
    burgers_asm.sd_nu = nu;
    burgers_asm.beta = beta;
    burgers_asm.frechet_beta = frechet_beta;
    burgers_asm.theta = theta;
    burgers_asm.sd_delta = sd_delta;
    burgers_asm.set_sd_v_norm(dummy_conv);
    burgers_asm.reg_eps = reg_eps;
    burgers_asm.frechet_material = frechet_material;
    burgers_asm.aT = a_T;

    dom_asm.compile_all_elements();
    dom_asm.assemble(burgers_asm);
    matrix_ref.scale(matrix_ref, alpha);


    std::vector<int> coloring = VoxelAssembly::UnitCubeColoring<ShapeType>::create_coloring(this->level);
    VoxelAssembly::test_coloring(mesh, coloring);
    VoxelAssembly::VoxelBurgersVeloMaterialAssembler<Lagrange2SpaceType, DataType, IndexType> burgers_assembler(space, coloring, Lagrange2SpaceType::world_dim == 3 ? 8 : 4);
    burgers_assembler.deformation = deformation;
    burgers_assembler.nu = nu;
    burgers_assembler.sd_nu = nu;
    burgers_assembler.beta = beta;
    burgers_assembler.frechet_beta = frechet_beta;
    burgers_assembler.theta = theta;
    burgers_assembler.sd_delta = sd_delta;
    burgers_assembler.set_sd_v_norm(dummy_conv);
    burgers_assembler.reg_eps = reg_eps;
    burgers_assembler.frechet_material = frechet_material;
    burgers_assembler.mu_0 = mu_0;
    burgers_assembler.exp = exp;
    burgers_assembler.lambda = lambda;
    burgers_assembler.a_T = a_T;
    burgers_assembler.yasuda_a = yasuda_a;
    burgers_assembler.mu_inf = mu_inf;
    burgers_assembler.material_type = material_type_;

    TEST_CHECK_EQUAL_WITHIN_EPS(burgers_asm.sd_v_norm, burgers_assembler.sd_v_norm, Math::Limits<DataType>::epsilon()*DataType(100));

    burgers_assembler.assemble_matrix1(matrix_new, dummy_conv, space, Cubature::DynamicFactory(cubature_name), alpha);

    this->compare_matrices(matrix_new, matrix_ref);

  }

  virtual void run_defect_compare() const override
  {
    typedef LAFEM::DenseVectorBlocked<DataType, IndexType, ShapeType::dimension> VectorType;

    Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(Index(this->level));
    MeshType mesh(mesh_factory);

    Geometry::BoundaryFactory<MeshType> boundary_factory(mesh);
    MeshPartType boundary(boundary_factory);
    TrafoType trafo(mesh);
    Lagrange2SpaceType space(trafo);

    Random::SeedType seed(Random::SeedType(time(nullptr)));
    Random rng(seed);

    VectorType dummy_conv = VectorType(rng, space.get_num_dofs(), DataType(-1), DataType(1));
    VectorType vector_ref = dummy_conv.clone(LAFEM::CloneMode::Deep);
    VectorType vector_new = dummy_conv.clone(LAFEM::CloneMode::Deep);
    VectorType primal_vec = VectorType(rng, space.get_num_dofs(), DataType(-1), DataType(1));

    DataType nu = DataType(-0.78);
    DataType beta = DataType(0.93);
    DataType theta = DataType(1.3);
    bool deformation = true;
    DataType alpha = DataType(0.66);
    DataType reg_eps = DataType(1E-100);
    DataType frechet_material = DataType(0.7);
    DataType mu_0 = nu;
    DataType exp = DataType(1.7);
    DataType lambda = DataType(0.8);
    DataType a_T = DataType(1.2);
    DataType yasuda_a = DataType(1.8);
    DataType mu_inf = DataType(3.5);


    vector_ref.format();
    vector_new.format();

    String cubature_name = "auto-degree:5";

    Assembly::DomainAssembler<TrafoType> dom_asm(trafo);
    Assembly::BurgersVeloMaterialBlockedVectorAssemblyJob burgers_asm(vector_ref, primal_vec, dummy_conv, space, cubature_name,
      VoxelAssembly::Intern::ViscoFunctorNew<DataType, material_type_>{mu_0, exp, lambda},
      VoxelAssembly::Intern::ViscoDFunctorNew<DataType, material_type_>{mu_0, exp, lambda});

    burgers_asm.deformation = deformation;
    burgers_asm.nu = nu;
    burgers_asm.sd_nu = nu;
    burgers_asm.beta = beta;
    burgers_asm.theta = theta;
    burgers_asm.reg_eps = reg_eps;
    burgers_asm.aT = a_T;

    dom_asm.compile_all_elements();
    dom_asm.assemble(burgers_asm);
    vector_ref.scale(vector_ref, alpha);


    std::vector<int> coloring = VoxelAssembly::UnitCubeColoring<ShapeType>::create_coloring(this->level);
    VoxelAssembly::test_coloring(mesh, coloring);
    VoxelAssembly::VoxelBurgersVeloMaterialAssembler<Lagrange2SpaceType, DataType, IndexType> burgers_assembler(space, coloring, Lagrange2SpaceType::world_dim == 3 ? 8 : 4);
    burgers_assembler.deformation = deformation;
    burgers_assembler.nu = nu;
    burgers_assembler.beta = beta;
    burgers_assembler.theta = theta;
    burgers_assembler.reg_eps = reg_eps;
    burgers_assembler.frechet_material = frechet_material;
    burgers_assembler.mu_0 = mu_0;
    burgers_assembler.exp = exp;
    burgers_assembler.lambda = lambda;
    burgers_assembler.a_T = a_T;
    burgers_assembler.yasuda_a = yasuda_a;
    burgers_assembler.mu_inf = mu_inf;
    burgers_assembler.material_type = material_type_;

    burgers_assembler.assemble_vector(vector_new, dummy_conv, primal_vec, space, Cubature::DynamicFactory(cubature_name), alpha);

    this->compare_defects(vector_new, vector_ref);

  }
}; // class VoxelDefoAssemblyTest

constexpr int _lvl_2d = 3;
constexpr int _lvl_3d = 2;

VoxelPoissonAssemblyTest<double, std::uint32_t, Shape::Hypercube<2>> poisson_vassembly_double_uint32_quadliteral_test_generic(_lvl_2d, PreferredBackend::generic);
VoxelPoissonAssemblyTest<float, std::uint32_t, Shape::Hypercube<2>> poisson_vassembly_float_uint32_quadliteral_test_generic(_lvl_2d, PreferredBackend::generic);
VoxelPoissonAssemblyTest<double, std::uint64_t, Shape::Hypercube<2>> poisson_vassembly_double_uint64_quadliteral_test_generic(_lvl_2d, PreferredBackend::generic);
VoxelPoissonAssemblyTest<float, std::uint64_t, Shape::Hypercube<2>> poisson_vassembly_float_uint64_quadliteral_test_generic(_lvl_2d, PreferredBackend::generic);
VoxelPoissonAssemblyTest<double, std::uint32_t, Shape::Hypercube<3>> poisson_vassembly_double_uint32_hexaedral_test_generic(_lvl_3d, PreferredBackend::generic);
VoxelPoissonAssemblyTest<float, std::uint32_t, Shape::Hypercube<3>> poisson_vassembly_float_uint32_hexaedral_test_generic(_lvl_3d, PreferredBackend::generic);
VoxelPoissonAssemblyTest<double, std::uint64_t, Shape::Hypercube<3>> poisson_vassembly_double_uint64_hexaedral_test_generic(_lvl_3d, PreferredBackend::generic);
VoxelPoissonAssemblyTest<float, std::uint64_t, Shape::Hypercube<3>> poisson_vassembly_float_uint64_hexaedral_test_generic(_lvl_3d, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
VoxelPoissonAssemblyTest<double, std::uint32_t, Shape::Hypercube<2>> poisson_vassembly_double_uint32_quadliteral_test_cuda(_lvl_2d, PreferredBackend::cuda);
VoxelPoissonAssemblyTest<float, std::uint32_t, Shape::Hypercube<2>> poisson_vassembly_float_uint32_quadliteral_test_cuda(_lvl_2d, PreferredBackend::cuda);
VoxelPoissonAssemblyTest<double, std::uint64_t, Shape::Hypercube<2>> poisson_vassembly_double_uint64_quadliteral_test_cuda(_lvl_2d, PreferredBackend::cuda);
VoxelPoissonAssemblyTest<float, std::uint64_t, Shape::Hypercube<2>> poisson_vassembly_float_uint64_quadliteral_test_cuda(_lvl_2d, PreferredBackend::cuda);
VoxelPoissonAssemblyTest<double, std::uint32_t, Shape::Hypercube<3>> poisson_vassembly_double_uint32_hexaedral_test_cuda(_lvl_3d, PreferredBackend::cuda);
VoxelPoissonAssemblyTest<float, std::uint32_t, Shape::Hypercube<3>> poisson_vassembly_float_uint32_hexaedral_test_cuda(_lvl_3d, PreferredBackend::cuda);
VoxelPoissonAssemblyTest<double, std::uint64_t, Shape::Hypercube<3>> poisson_vassembly_double_uint64_hexaedral_test_cuda(_lvl_3d, PreferredBackend::cuda);
VoxelPoissonAssemblyTest<float, std::uint64_t, Shape::Hypercube<3>> poisson_vassembly_float_uint64_hexaedral_test_cuda(_lvl_3d, PreferredBackend::cuda);
#endif

VoxelDefoAssemblyTest<double, std::uint32_t, Shape::Hypercube<2>> defo_vassembly_double_uint32_quadliteral_test_generic(_lvl_2d, PreferredBackend::generic);
#ifndef DEBUG
VoxelDefoAssemblyTest<float, std::uint32_t, Shape::Hypercube<2>> defo_vassembly_float_uint32_quadliteral_test_generic(_lvl_2d, PreferredBackend::generic);
VoxelDefoAssemblyTest<double, std::uint64_t, Shape::Hypercube<2>> defo_vassembly_double_uint64_quadliteral_test_generic(_lvl_2d, PreferredBackend::generic);
VoxelDefoAssemblyTest<float, std::uint64_t, Shape::Hypercube<2>> defo_vassembly_float_uint64_quadliteral_test_generic(_lvl_2d, PreferredBackend::generic);
VoxelDefoAssemblyTest<double, std::uint32_t, Shape::Hypercube<3>> defo_vassembly_double_uint32_hexaedral_test_generic(_lvl_3d, PreferredBackend::generic);
VoxelDefoAssemblyTest<float, std::uint32_t, Shape::Hypercube<3>> defo_vassembly_float_uint32_hexaedral_test_generic(_lvl_3d, PreferredBackend::generic);
VoxelDefoAssemblyTest<double, std::uint64_t, Shape::Hypercube<3>> defo_vassembly_double_uint64_hexaedral_test_generic(_lvl_3d, PreferredBackend::generic);
VoxelDefoAssemblyTest<float, std::uint64_t, Shape::Hypercube<3>> defo_vassembly_float_uint64_hexaedral_test_generic(_lvl_3d, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
VoxelDefoAssemblyTest<double, std::uint32_t, Shape::Hypercube<2>> defo_vassembly_double_uint32_quadliteral_test_cuda(_lvl_2d, PreferredBackend::cuda);
#ifndef DEBUG
VoxelDefoAssemblyTest<float, std::uint32_t, Shape::Hypercube<2>> defo_vassembly_float_uint32_quadliteral_test_cuda(_lvl_2d, PreferredBackend::cuda);
VoxelDefoAssemblyTest<double, std::uint64_t, Shape::Hypercube<2>> defo_vassembly_double_uint64_quadliteral_test_cuda(_lvl_2d, PreferredBackend::cuda);
VoxelDefoAssemblyTest<float, std::uint64_t, Shape::Hypercube<2>> defo_vassembly_float_uint64_quadliteral_test_cuda(_lvl_2d, PreferredBackend::cuda);
#endif
VoxelDefoAssemblyTest<float, std::uint32_t, Shape::Hypercube<3>> defo_vassembly_float_uint32_hexaedral_test_cuda(_lvl_3d, PreferredBackend::cuda);
//only run one 3D gpu test since these have very long runtimes in debug mode...
#ifndef DEBUG
VoxelDefoAssemblyTest<double, std::uint32_t, Shape::Hypercube<3>> defo_vassembly_double_uint32_hexaedral_test_cuda(_lvl_3d, PreferredBackend::cuda);
VoxelDefoAssemblyTest<double, std::uint64_t, Shape::Hypercube<3>> defo_vassembly_double_uint64_hexaedral_test_cuda(_lvl_3d, PreferredBackend::cuda);
VoxelDefoAssemblyTest<float, std::uint64_t, Shape::Hypercube<3>> defo_vassembly_float_uint64_hexaedral_test_cuda(_lvl_3d, PreferredBackend::cuda);
#endif
#endif

VoxelBurgersAssemblyTest<double, std::uint32_t, Shape::Hypercube<2>> burgers_vassembly_double_uint32_quadliteral_test_generic(_lvl_2d, PreferredBackend::generic);
VoxelBurgersAssemblyTest<float, std::uint32_t, Shape::Hypercube<2>> burgers_vassembly_float_uint32_quadliteral_test_generic(_lvl_2d, PreferredBackend::generic);
VoxelBurgersAssemblyTest<double, std::uint64_t, Shape::Hypercube<2>> burgers_vassembly_double_uint64_quadliteral_test_generic(_lvl_2d, PreferredBackend::generic);
VoxelBurgersAssemblyTest<float, std::uint64_t, Shape::Hypercube<2>> burgers_vassembly_float_uint64_quadliteral_test_generic(_lvl_2d, PreferredBackend::generic);
VoxelBurgersAssemblyTest<double, std::uint32_t, Shape::Hypercube<3>> burgers_vassembly_double_uint32_hexaedral_test_generic(_lvl_3d, PreferredBackend::generic);
#ifndef DEBUG
VoxelBurgersAssemblyTest<float, std::uint32_t, Shape::Hypercube<3>> burgers_vassembly_float_uint32_hexaedral_test_generic(_lvl_3d, PreferredBackend::generic);
VoxelBurgersAssemblyTest<double, std::uint64_t, Shape::Hypercube<3>> burgers_vassembly_double_uint64_hexaedral_test_generic(_lvl_3d, PreferredBackend::generic);
VoxelBurgersAssemblyTest<float, std::uint64_t, Shape::Hypercube<3>> burgers_vassembly_float_uint64_hexaedral_test_generic(_lvl_3d, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
VoxelBurgersAssemblyTest<double, std::uint32_t, Shape::Hypercube<2>> burgers_vassembly_double_uint32_quadliteral_test_cuda(_lvl_2d, PreferredBackend::cuda);
VoxelBurgersAssemblyTest<float, std::uint32_t, Shape::Hypercube<2>> burgers_vassembly_float_uint32_quadliteral_test_cuda(_lvl_2d, PreferredBackend::cuda);
VoxelBurgersAssemblyTest<double, std::uint64_t, Shape::Hypercube<2>> burgers_vassembly_double_uint64_quadliteral_test_cuda(_lvl_2d, PreferredBackend::cuda);
VoxelBurgersAssemblyTest<float, std::uint64_t, Shape::Hypercube<2>> burgers_vassembly_float_uint64_quadliteral_test_cuda(_lvl_2d, PreferredBackend::cuda);
VoxelBurgersAssemblyTest<float, std::uint32_t, Shape::Hypercube<3>> burgers_vassembly_float_uint32_hexaedral_test_cuda(_lvl_3d, PreferredBackend::cuda);
#ifndef DEBUG
VoxelBurgersAssemblyTest<double, std::uint32_t, Shape::Hypercube<3>> burgers_vassembly_double_uint32_hexaedral_test_cuda(_lvl_3d, PreferredBackend::cuda);
VoxelBurgersAssemblyTest<double, std::uint64_t, Shape::Hypercube<3>> burgers_vassembly_double_uint64_hexaedral_test_cuda(_lvl_3d, PreferredBackend::cuda);
VoxelBurgersAssemblyTest<float, std::uint64_t, Shape::Hypercube<3>> burgers_vassembly_float_uint64_hexaedral_test_cuda(_lvl_3d, PreferredBackend::cuda);
#endif
#endif

VoxelBurgersVeloMaterialAssemblyTest<double, std::uint32_t, Shape::Hypercube<2>, VoxelAssembly::MaterialType::carreau> burgers_vm_vassembly_double_uint32_quadliteral_test_generic(_lvl_2d, PreferredBackend::generic);
VoxelBurgersVeloMaterialAssemblyTest<float, std::uint32_t, Shape::Hypercube<2>, VoxelAssembly::MaterialType::carreau> burgers_vm_vassembly_float_uint32_quadliteral_test_generic(_lvl_2d, PreferredBackend::generic);
VoxelBurgersVeloMaterialAssemblyTest<double, std::uint64_t, Shape::Hypercube<2>, VoxelAssembly::MaterialType::carreau> burgers_vm_vassembly_double_uint64_quadliteral_test_generic(_lvl_2d, PreferredBackend::generic);
VoxelBurgersVeloMaterialAssemblyTest<float, std::uint64_t, Shape::Hypercube<2>, VoxelAssembly::MaterialType::carreau> burgers_vm_vassembly_float_uint64_quadliteral_test_generic(_lvl_2d, PreferredBackend::generic);
VoxelBurgersVeloMaterialAssemblyTest<double, std::uint32_t, Shape::Hypercube<3>, VoxelAssembly::MaterialType::carreau> burgers_vm_vassembly_double_uint32_hexaedral_test_generic(_lvl_3d, PreferredBackend::generic);
VoxelBurgersVeloMaterialAssemblyTest<float, std::uint32_t, Shape::Hypercube<3>, VoxelAssembly::MaterialType::carreau> burgers_vm_vassembly_float_uint32_hexaedral_test_generic(_lvl_3d, PreferredBackend::generic);
VoxelBurgersVeloMaterialAssemblyTest<double, std::uint64_t, Shape::Hypercube<3>, VoxelAssembly::MaterialType::carreau> burgers_vm_vassembly_double_uint64_hexaedral_test_generic(_lvl_3d, PreferredBackend::generic);
VoxelBurgersVeloMaterialAssemblyTest<float, std::uint64_t, Shape::Hypercube<3>, VoxelAssembly::MaterialType::carreau> burgers_vm_vassembly_float_uint64_hexaedral_test_generic(_lvl_3d, PreferredBackend::generic);
#ifdef FEAT_HAVE_CUDA
VoxelBurgersVeloMaterialAssemblyTest<double, std::uint32_t, Shape::Hypercube<2>, VoxelAssembly::MaterialType::carreau> burgers_vm_vassembly_double_uint32_quadliteral_test_cuda(_lvl_2d, PreferredBackend::cuda);
VoxelBurgersVeloMaterialAssemblyTest<float, std::uint32_t, Shape::Hypercube<2>, VoxelAssembly::MaterialType::carreau> burgers_vm_vassembly_float_uint32_quadliteral_test_cuda(_lvl_2d, PreferredBackend::cuda);
VoxelBurgersVeloMaterialAssemblyTest<double, std::uint64_t, Shape::Hypercube<2>, VoxelAssembly::MaterialType::carreau> burgers_vm_vassembly_double_uint64_quadliteral_test_cuda(_lvl_2d, PreferredBackend::cuda);
VoxelBurgersVeloMaterialAssemblyTest<float, std::uint64_t, Shape::Hypercube<2>, VoxelAssembly::MaterialType::carreau> burgers_vm_vassembly_float_uint64_quadliteral_test_cuda(_lvl_2d, PreferredBackend::cuda);
VoxelBurgersVeloMaterialAssemblyTest<float, std::uint32_t, Shape::Hypercube<3>, VoxelAssembly::MaterialType::carreau> burgers_vm_vassembly_float_uint32_hexaedral_test_cuda(_lvl_3d, PreferredBackend::cuda);
#ifndef DEBUG
VoxelBurgersVeloMaterialAssemblyTest<double, std::uint32_t, Shape::Hypercube<3>, VoxelAssembly::MaterialType::carreau> burgers_vm_vassembly_double_uint32_hexaedral_test_cuda(_lvl_3d, PreferredBackend::cuda);
VoxelBurgersVeloMaterialAssemblyTest<double, std::uint64_t, Shape::Hypercube<3>, VoxelAssembly::MaterialType::carreau> burgers_vm_vassembly_double_uint64_hexaedral_test_cuda(_lvl_3d, PreferredBackend::cuda);
VoxelBurgersVeloMaterialAssemblyTest<float, std::uint64_t, Shape::Hypercube<3>, VoxelAssembly::MaterialType::carreau> burgers_vm_vassembly_float_uint64_hexaedral_test_cuda(_lvl_3d, PreferredBackend::cuda);
#endif
#endif
