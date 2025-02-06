// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/assembly/mean_filter_assembler.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/mean_filter.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/random.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

template<typename DataType_, typename IndexType_>
class MeanFilterAssemblerTest :
  public UnitTest
{
  typedef LAFEM::MeanFilter<DataType_, IndexType_> FilterType;
  typedef LAFEM::DenseVector<DataType_, IndexType_> VectorType;

  typedef Geometry::ConformalMesh<Shape::Quadrilateral> QuadMesh;

  typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

  typedef Space::Lagrange2::Element<QuadTrafo> QuadSpaceQ2;

public:
  MeanFilterAssemblerTest(PreferredBackend backend) :
    UnitTest("MeanFilterAssemblerTest", Type::Traits<DataType_>::name(), Type::Traits<IndexType_>::name(), backend)
  {
  }

  virtual ~MeanFilterAssemblerTest()
  {
  }

  virtual void run() const override
  {
    test_unit_2d();
  }

  void test_unit_2d() const
  {
    //const DataType_ tol = Math::pow(Math::eps<DataType_>(), DataType_(0.75));
    DataType_ tol;
    if (typeid(DataType_) == typeid(double))
    {
      tol = Math::pow(Math::eps<DataType_>(), DataType_(0.75));
    }
    #ifdef FEAT_HAVE_QUADMATH
    else if (typeid(DataType_) == typeid(__float128))
    {
      tol = Math::pow(Math::eps<DataType_>(), DataType_(0.75));
    }
    #endif
    else
    {
      tol = Math::pow(Math::eps<DataType_>(), DataType_(0.4));
    }

    // create coarse mesh
    Geometry::RefinedUnitCubeFactory<QuadMesh> unit_factory(3);
    QuadMesh mesh(unit_factory);
    QuadTrafo trafo(mesh);
    QuadSpaceQ2 space(trafo);

    // create a cubature factory
    Cubature::DynamicFactory cubature_factory("gauss-legendre:3");

    FilterType filter;

    VectorType vec_prim, vec_dual;
    Assembly::MeanFilterAssembler::assemble(vec_prim, vec_dual, space, cubature_factory);


    Assembly::MeanFilterAssembler::assemble(filter, space, cubature_factory, DataType_(0.5));

    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";

    VectorType vec_sol(rng, space.get_num_dofs(), DataType_(-2.0), DataType_(1.0));
    VectorType vec_rhs(rng, space.get_num_dofs(), DataType_(-2.0), DataType_(1.0));
    VectorType vec_def(rng, space.get_num_dofs(), DataType_(-2.0), DataType_(1.0));
    VectorType vec_cor(rng, space.get_num_dofs(), DataType_(-2.0), DataType_(1.0));

    // apply filter to all vectors
    filter.filter_sol(vec_sol);
    filter.filter_cor(vec_cor);
    filter.filter_rhs(vec_rhs);
    filter.filter_def(vec_def);

    // check primal integrals
    DataType_ int_sol = vec_sol.dot(vec_dual);
    TEST_CHECK_EQUAL_WITHIN_EPS(int_sol, DataType_(0.5), tol);
    DataType_ int_cor = vec_cor.dot(vec_dual);
    TEST_CHECK_EQUAL_WITHIN_EPS(int_cor, DataType_(0.0), tol);

    DataType_ int_rhs = vec_rhs.dot(vec_prim);
    TEST_CHECK_EQUAL_WITHIN_EPS(int_rhs, DataType_(0.0), tol);
    DataType_ int_def = vec_def.dot(vec_prim);
    TEST_CHECK_EQUAL_WITHIN_EPS(int_def, DataType_(0.0), tol);
  }
};

MeanFilterAssemblerTest <float, std::uint32_t> mean_filter_asm_test_float_uint32(PreferredBackend::generic);
MeanFilterAssemblerTest <double, std::uint32_t> mean_filter_asm_test_double_uint32(PreferredBackend::generic);
MeanFilterAssemblerTest <float, std::uint64_t> mean_filter_asm_test_float_uint64(PreferredBackend::generic);
MeanFilterAssemblerTest <double, std::uint64_t> mean_filter_asm_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
MeanFilterAssemblerTest <float, std::uint64_t> mkl_mean_filter_asm_test_float_uint64(PreferredBackend::mkl);
MeanFilterAssemblerTest <double, std::uint64_t> mkl_mean_filter_asm_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
MeanFilterAssemblerTest <__float128, std::uint32_t> mean_filter_asm_test_float128_uint32(PreferredBackend::generic);
MeanFilterAssemblerTest <__float128, std::uint64_t> mean_filter_asm_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
MeanFilterAssemblerTest <Half, std::uint32_t> mean_filter_asm_test_half_uint32(PreferredBackend::generic);
MeanFilterAssemblerTest <Half, std::uint64_t> mean_filter_asm_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
MeanFilterAssemblerTest <float, std::uint32_t> cuda_mean_filter_asm_test_float_uint32(PreferredBackend::cuda);
MeanFilterAssemblerTest <double, std::uint32_t> cuda_mean_filter_asm_test_double_uint32(PreferredBackend::cuda);
MeanFilterAssemblerTest <float, std::uint64_t> cuda_mean_filter_asm_test_float_uint64(PreferredBackend::cuda);
MeanFilterAssemblerTest <double, std::uint64_t> cuda_mean_filter_asm_test_double_uint64(PreferredBackend::cuda);
#endif
