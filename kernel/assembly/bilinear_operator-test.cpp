// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/structured_mesh.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/pointstar_structure.hpp>
#include <kernel/space/argyris/element.hpp>
#include <kernel/space/bogner_fox_schmit/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/cro_rav_ran_tur/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/math.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

template<typename Trafo_>
using MyP0 = Space::Discontinuous::Element<Trafo_, Space::Discontinuous::Variant::StdPolyP<0>>;


template<typename DT_, typename IT_>
class BilinearOperatorTest :
  public UnitTest
{
  typedef DT_ DataType_;
  typedef IT_ IndexType_;
  typedef LAFEM::SparseMatrixCSR<DataType_, IndexType_> MatrixType;
  typedef LAFEM::DenseVector<DataType_, IndexType_> VectorType;

  typedef Geometry::ConformalMesh<Shape::Quadrilateral, 2, DataType_> QuadMesh;

  typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

  typedef Space::Discontinuous::Element<QuadTrafo> QuadSpaceQ0;
  typedef Space::Lagrange1::Element<QuadTrafo> QuadSpaceQ1;

public:
  BilinearOperatorTest(PreferredBackend backend) :
    UnitTest("BilinearOperatorTest", Type::Traits<DataType_>::name(), Type::Traits<IndexType_>::name(), backend)
  {
  }

  virtual ~BilinearOperatorTest()
  {
  }

  virtual void run() const override
  {
    test_unit_2d();

    // Argyris has polynomial degree 6, so the integration for the indentity operator has only to be exact to
    // degree 6 + 6 = 12
    test_apply1<Space::Argyris::Element, Shape::Simplex<2>, Assembly::Common::IdentityOperator, 2>(2,12);
    // degree 1 + 1 - 2 = 0
    test_apply1<Space::Lagrange1::Element, Shape::Hypercube<3>, Assembly::Common::LaplaceOperator, 3>(2, 0);

    // RT has degree 1, BFC has degree 3, so the integration for the Laplace operator has only to be exact to
    // degree 1 + 3 - 2 = 2
    test_apply2<Space::CroRavRanTur::Element, Space::BognerFoxSchmit::Element, Shape::Hypercube<2>, Assembly::Common::LaplaceOperator, 2>(3, 2);
    // degree 0 + 1 = 1
    test_apply2<MyP0, Space::Lagrange1::Element, Shape::Simplex<3>, Assembly::Common::IdentityOperator, 3>(1, 1);
  }

  template<template<typename> class SpaceType_, typename ShapeType_, typename OperatorType_, int dim_>
  void test_apply1(const Index level, const int cubature_degree) const
  {
    typedef Geometry::ConformalMesh<ShapeType_, dim_, DataType_> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef SpaceType_<TrafoType> SpaceType;

    Geometry::RefinedUnitCubeFactory<MeshType> unit_factory(level);
    MeshType my_mesh(unit_factory);
    TrafoType my_trafo(my_mesh);
    SpaceType my_space(my_trafo);

    // create a matrix
    MatrixType matrix;
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix, my_space);
    matrix.format();

    // create a cubature factory
    Cubature::DynamicFactory cubature_factory("auto-degree:"+stringify(cubature_degree));

    // Our bilinear operator
    OperatorType_ my_operator;

    // Generate bogus coefficient vector for the FE function
    VectorType x(my_space.get_num_dofs(), DataType_(0));
    // This factor is there so that the l2-norm of x does not get too big
    DataType_ fac(DataType_(1)/DataType_(my_space.get_num_dofs()));
    for(Index i(0); i < x.size(); ++i)
      x(i, fac*(DataType_(2)*DataType_(i) - Math::sqrt(DataType_(i))));

    // Compute reference solution by assembling the matrix and multiplying x to it
    Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix, my_operator, my_space, cubature_factory);
    VectorType ref(my_space.get_num_dofs(), DataType_(0));
    matrix.apply(ref, x);

    // Compute the result given by the apply1() routine
    VectorType res(my_space.get_num_dofs(), DataType_(0));
    Assembly::BilinearOperatorAssembler::apply1(res, x, my_operator, my_space, cubature_factory);

    // Compute || ref - res ||_l2
    res.axpy(ref, -DataType_(1));
    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.8));
    TEST_CHECK_EQUAL_WITHIN_EPS(res.norm2(), DataType_(0), eps);

  }

  template
  <
    template<typename> class TestSpaceType_,
    template<typename> class TrialSpaceType_,
    typename ShapeType_,
    typename OperatorType_,
    int dim_
  >
  void test_apply2(const Index level, const int cubature_degree) const
  {
    typedef Geometry::ConformalMesh<ShapeType_, dim_, DataType_> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef TestSpaceType_<TrafoType> TestSpaceType;
    typedef TrialSpaceType_<TrafoType> TrialSpaceType;

    Geometry::RefinedUnitCubeFactory<MeshType> unit_factory(level);
    MeshType my_mesh(unit_factory);
    TrafoType my_trafo(my_mesh);
    TestSpaceType my_test_space(my_trafo);
    TrialSpaceType my_trial_space(my_trafo);

    // create a matrix
    MatrixType matrix;
    Assembly::SymbolicAssembler::assemble_matrix_std2(matrix, my_test_space, my_trial_space);
    matrix.format();

    // create a cubature factory
    Cubature::DynamicFactory cubature_factory("auto-degree:"+stringify(cubature_degree));

    // Create bogus coefficient vector for the FE function
    VectorType x(my_trial_space.get_num_dofs(), DataType_(0));
    // This factor is there so that the l2-norm of x does not get too big
    DataType_ fac(DataType_(1)/DataType_(my_trial_space.get_num_dofs()));
    for(Index i(0); i < x.size(); ++i)
      x(i, fac*(DataType_(2)*DataType_(i) - Math::sqrt(DataType_(i))));

    // Our bilinear operator
    OperatorType_ my_operator;

    // Compute reference solution by assembling the matrix and multiplying x to it
    Assembly::BilinearOperatorAssembler::assemble_matrix2
      (matrix, my_operator, my_test_space, my_trial_space, cubature_factory);

    VectorType ref(my_test_space.get_num_dofs(), DataType_(0));
    matrix.apply(ref, x);

    // Compute the result given by the apply2() routine
    VectorType res(my_test_space.get_num_dofs(), DataType_(0));
    Assembly::BilinearOperatorAssembler::apply2
      (res, x, my_operator, my_test_space, my_trial_space, cubature_factory);

    // Compute || ref - res ||_l2
    res.axpy(ref, -DataType_(1));

    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.8));
    TEST_CHECK_EQUAL_WITHIN_EPS(res.norm2(), DataType_(0), eps);

  }


  void test_unit_2d() const
  {
    // create coarse mesh
    Geometry::RefinedUnitCubeFactory<QuadMesh> unit_factory(3);
    QuadMesh mesh(unit_factory);

    // run tests
    test_unit_2d_q0(mesh);
    test_unit_2d_q1(mesh);
  }

  void test_unit_2d_q0(QuadMesh& mesh) const
  {
    // compute eps
    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.9));

    // create trafo
    QuadTrafo trafo(mesh);

    // create space
    QuadSpaceQ0 space(trafo);

    // create a matrix
    MatrixType matrix;
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix, space);
    matrix.format();

    // create a cubature factory
    Cubature::DynamicFactory cubature_factory("barycentre");

    // assemble the identity operator with barycentre cubature rule
    Assembly::Common::IdentityOperator operat;
    Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix, operat, space, cubature_factory);

    // create CSR-matrix of matrix
    LAFEM::SparseMatrixCSR<DataType_, IndexType_> tmp_matrix;
    tmp_matrix.convert(matrix);
    // fetch the matrix arrays
    DataType_* data = tmp_matrix.val();

    const DataType_ v = DataType_(1) / DataType_(mesh.get_num_entities(2));

    // loop over all matrix rows
    for(Index i(0); i < tmp_matrix.rows(); ++i)
    {
      // validate entry
      TEST_CHECK_EQUAL_WITHIN_EPS(data[i], v, eps);
    }
  }

  void test_unit_2d_q1(QuadMesh& mesh) const
  {
    // compute eps
    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.9));

    // create trafo
    QuadTrafo trafo(mesh);

    // create space
    QuadSpaceQ1 space(trafo);

    // create two matrices
    MatrixType matrix_1, matrix_2;
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix_1, space);
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix_2, space);
    matrix_1.format();
    matrix_2.format();

    // create a cubature factory
    Cubature::DynamicFactory cubature_factory_trz("trapezoidal");
    Cubature::DynamicFactory cubature_factory_gl2("gauss-legendre:2");

    // assemble the identity operator with barycentre cubature rule
    Assembly::Common::LaplaceOperator operat;
    Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix_1, operat, space, cubature_factory_trz);
    Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix_2, operat, space, cubature_factory_gl2);

    // get mesh element count
    Index num_verts = mesh.get_num_entities(0);
    Index num_quads = mesh.get_num_entities(2);

    const Geometry::IndexSet<4>& vatq(mesh.template get_index_set<2,0>());

    // create a temporary array to count the number of quads adjacent to a vertex
    LAFEM::DenseVector<Index> qatv(num_verts, Index(0));
    for(Index i(0); i < num_quads; ++i)
    {
      for(int j(0); j < 4; ++j)
        qatv(vatq(i,j), qatv(vatq(i,j)) + 1);
    }

    // create a dof-mapping
    typename QuadSpaceQ1::DofMappingType dof_mapping(space);

    // create local matrix data
    Tiny::Matrix<DataType_,4,4> lmd1, lmd2;
    typename MatrixType::GatherAxpy gather1(matrix_1), gather2(matrix_2);

    // some constants
    static const DataType_ zero = DataType_(0);
    static const DataType_ one = DataType_(1);
    static const DataType_ _12 = DataType_(1) / DataType_(2);
    static const DataType_ _13 = DataType_(1) / DataType_(3);
    static const DataType_ _16 = DataType_(1) / DataType_(6);

    // loop over all quads
    for(Index cell(0); cell < num_quads; ++cell)
    {
      // fetch the local matrix
      dof_mapping.prepare(cell);
      lmd1.format();
      lmd2.format();
      gather1(lmd1, dof_mapping, dof_mapping);
      gather2(lmd2, dof_mapping, dof_mapping);

      // loop over all 4x4 entries
      for(int i(0); i < 4; ++i)
      {
        Index nvi = qatv(vatq(cell,i));
        for(int j(0); j < 4; ++j)
        {
          Index nvj = qatv(vatq(cell,j));
          if(i == j)
          {
            // main diagonal entry
            TEST_CHECK_EQUAL_WITHIN_EPS(lmd1(i,j), DataType_(nvi), eps);
            TEST_CHECK_EQUAL_WITHIN_EPS(lmd2(i,j), DataType_(2*nvi)*_13, eps);
          }
          else if((i^j) == 3)
          {
            // off-diagonal entry, two vertices sharing no common edge
            TEST_CHECK_EQUAL_WITHIN_EPS(lmd1(i,j), zero, eps);
            TEST_CHECK_EQUAL_WITHIN_EPS(lmd2(i,j), -_13, eps);
          }
          else if((nvi+nvj) > 4)
          {
            // off-diagonal entry, two vertices sharing an inner edge
            TEST_CHECK_EQUAL_WITHIN_EPS(lmd1(i,j), -one, eps);
            TEST_CHECK_EQUAL_WITHIN_EPS(lmd2(i,j), -_13, eps);
          }
          else
          {
            // off-diagonal entry, two vertices sharing a boundary edge
            TEST_CHECK_EQUAL_WITHIN_EPS(lmd1(i,j), -_12, eps);
            TEST_CHECK_EQUAL_WITHIN_EPS(lmd2(i,j), -_16, eps);
          }
        }
      }

      dof_mapping.finish();
    }
  }

};

BilinearOperatorTest <float, std::uint32_t> bilinear_operator_test_float_uint32(PreferredBackend::generic);
BilinearOperatorTest <float, std::uint64_t> bilinear_operator_test_float_uint64(PreferredBackend::generic);
BilinearOperatorTest <double, std::uint32_t> bilinear_operator_test_double_uint32(PreferredBackend::generic);
BilinearOperatorTest <double, std::uint64_t> bilinear_operator_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
BilinearOperatorTest <float, std::uint64_t> mkl_bilinear_operator_test_float_uint64(PreferredBackend::mkl);
BilinearOperatorTest <double, std::uint64_t> mkl_bilinear_operator_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
BilinearOperatorTest <__float128, std::uint32_t> bilinear_operator_test_float128_uint32(PreferredBackend::generic);
BilinearOperatorTest <__float128, std::uint64_t> bilinear_operator_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
BilinearOperatorTest <Half, std::uint32_t> bilinear_operator_test_half_uint32(PreferredBackend::generic);
BilinearOperatorTest <Half, std::uint64_t> bilinear_operator_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
BilinearOperatorTest <float, std::uint32_t> cuda_bilinear_operator_test_float_uint32(PreferredBackend::cuda);
BilinearOperatorTest <double, std::uint32_t> cuda_bilinear_operator_test_double_uint32(PreferredBackend::cuda);
BilinearOperatorTest <float, std::uint64_t> cuda_bilinear_operator_test_float_uint64(PreferredBackend::cuda);
BilinearOperatorTest <double, std::uint64_t> cuda_bilinear_operator_test_double_uint64(PreferredBackend::cuda);
#endif

template<typename DataType_, typename IndexType_>
class BandedBilinearOperatorTest :
  public UnitTest
{
  typedef LAFEM::SparseMatrixBanded<DataType_, IndexType_> MatrixType;
  typedef LAFEM::DenseVector<DataType_, IndexType_> VectorType;

  typedef Geometry::StructuredMesh<2, 2, DataType_> QuadMesh;

  typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

  typedef Space::Discontinuous::Element<QuadTrafo> QuadSpaceQ0;
  typedef Space::Lagrange1::Element<QuadTrafo> QuadSpaceQ1;

public:
  BandedBilinearOperatorTest(PreferredBackend backend) :
    UnitTest("BandedBilinearOperatorTest", Type::Traits<DataType_>::name(), Type::Traits<IndexType_>::name(), backend)
  {
  }

  virtual ~BandedBilinearOperatorTest()
  {
  }

  virtual void run() const override
  {
    test_unit_2d();
  }

  void test_unit_2d() const
  {
    // create a 4x4 structured mesh
    const Index num_slices[2] = {Index(4), Index(4)};
    QuadMesh mesh(num_slices);

    // initialize our 5x5 vertices
    auto& v = mesh.get_vertex_set();
    for(Index i(0), k(0); i <= num_slices[0]; ++i)
    {
      for(Index j(0); j <= num_slices[1]; ++j, ++k)
      {
        v[k][0] = DataType_(i) / DataType_(num_slices[0]);
        v[k][1] = DataType_(j) / DataType_(num_slices[1]);
      }
    }

    // run tests
    test_unit_2d_q0(mesh);
    test_unit_2d_q1(mesh);
  }

  void test_unit_2d_q0(QuadMesh& mesh) const
  {
    // compute eps
    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.9));

    // create trafo
    QuadTrafo trafo(mesh);

    // create space
    QuadSpaceQ0 space(trafo);

    // create a matrix
    std::vector<IndexType_> nsi;
    nsi.push_back(IndexType_(mesh.get_num_slices(0)));
    nsi.push_back(IndexType_(mesh.get_num_slices(1)));
    MatrixType matrix(LAFEM::PointstarStructureFE::template value<DataType_,IndexType_>(Index(0), nsi));
    matrix.format();

    // create a cubature factory
    Cubature::DynamicFactory cubature_factory("barycentre");

    // assemble the identity operator with barycentre cubature rule
    Assembly::Common::IdentityOperator operat;
    Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix, operat, space, cubature_factory);

    // fetch the matrix diagonal
    DataType_* data = matrix.val();
    const DataType_ v = DataType_(1) / DataType_(mesh.get_num_entities(2));

    // loop over all matrix rows
    for(Index i(0); i < matrix.rows(); ++i)
    {
      // validate entry
      TEST_CHECK_EQUAL_WITHIN_EPS(data[i], v, eps);
    }
  }

  void test_unit_2d_q1(QuadMesh& mesh) const
  {
    // compute eps
    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.9));

    // create trafo
    QuadTrafo trafo(mesh);

    // create space
    QuadSpaceQ1 space(trafo);

    // create two matrices
    std::vector<IndexType_> nsi;
    nsi.push_back(IndexType_(mesh.get_num_slices(0)));
    nsi.push_back(IndexType_(mesh.get_num_slices(1)));
    MatrixType matrix_1(LAFEM::PointstarStructureFE::template value<DataType_,IndexType_>(Index(1), nsi));

    matrix_1.format();
    MatrixType matrix_2(matrix_1.clone());

    // create a cubature factory
    Cubature::DynamicFactory cubature_factory_trz("trapezoidal");
    Cubature::DynamicFactory cubature_factory_gl2("gauss-legendre:2");

    // assemble the identity operator with barycentre cubature rule
    Assembly::Common::LaplaceOperator operat;
    Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix_1, operat, space, cubature_factory_trz);
    Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix_2, operat, space, cubature_factory_gl2);

    // get mesh element count
    Index num_verts = mesh.get_num_entities(0);
    Index num_quads = mesh.get_num_entities(2);

    const auto& vatq(mesh.template get_index_set<2,0>());

    // create a temporary array to count the number of quads adjacent to a vertex
    LAFEM::DenseVector<Index> qatv(num_verts, Index(0));
    for(Index i(0); i < num_quads; ++i)
    {
      for(int j(0); j < 4; ++j)
        qatv(vatq(i,j), qatv(vatq(i,j)) + 1);
    }

    // create a dof-mapping
    typename QuadSpaceQ1::DofMappingType dof_mapping(space);

    // create local matrix data
    Tiny::Matrix<DataType_,4,4> lmd1, lmd2;
    typename MatrixType::GatherAxpy gather1(matrix_1), gather2(matrix_2);

    // some constants
    static const DataType_ zero = DataType_(0);
    static const DataType_ one = DataType_(1);
    static const DataType_ _12 = DataType_(1) / DataType_(2);
    static const DataType_ _13 = DataType_(1) / DataType_(3);
    static const DataType_ _16 = DataType_(1) / DataType_(6);

    // loop over all quads
    for(Index cell(0); cell < num_quads; ++cell)
    {
      // fetch the local matrix
      dof_mapping.prepare(cell);
      lmd1.format();
      lmd2.format();
      gather1(lmd1, dof_mapping, dof_mapping);
      gather2(lmd2, dof_mapping, dof_mapping);

      // loop over all 4x4 entries
      for(int i(0); i < 4; ++i)
      {
        Index nvi = qatv(vatq(cell,i));
        for(int j(0); j < 4; ++j)
        {
          Index nvj = qatv(vatq(cell,j));
          if(i == j)
          {
            // main diagonal entry
            TEST_CHECK_EQUAL_WITHIN_EPS(lmd1(i,j), DataType_(nvi), eps);
            TEST_CHECK_EQUAL_WITHIN_EPS(lmd2(i,j), DataType_(2*nvi)*_13, eps);
          }
          else if((i^j) == 3)
          {
            // off-diagonal entry, two vertices sharing no common edge
            TEST_CHECK_EQUAL_WITHIN_EPS(lmd1(i,j), zero, eps);
            TEST_CHECK_EQUAL_WITHIN_EPS(lmd2(i,j), -_13, eps);
          }
          else if((nvi+nvj) > 4)
          {
            // off-diagonal entry, two vertices sharing an inner edge
            TEST_CHECK_EQUAL_WITHIN_EPS(lmd1(i,j), -one, eps);
            TEST_CHECK_EQUAL_WITHIN_EPS(lmd2(i,j), -_13, eps);
          }
          else
          {
            // off-diagonal entry, two vertices sharing a boundary edge
            TEST_CHECK_EQUAL_WITHIN_EPS(lmd1(i,j), -_12, eps);
            TEST_CHECK_EQUAL_WITHIN_EPS(lmd2(i,j), -_16, eps);
          }
        }
      }

      dof_mapping.finish();
    }
  }
};

BandedBilinearOperatorTest <float, std::uint32_t> banded_bilinear_operator_test_float_uint32(PreferredBackend::generic);
BandedBilinearOperatorTest <float, std::uint64_t> banded_bilinear_operator_test_float_uint64(PreferredBackend::generic);
BandedBilinearOperatorTest <double, std::uint32_t> banded_bilinear_operator_test_double_uint32(PreferredBackend::generic);
BandedBilinearOperatorTest <double, std::uint64_t> banded_bilinear_operator_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
BandedBilinearOperatorTest <float, std::uint64_t> mkl_banded_bilinear_operator_test_float_uint64(PreferredBackend::mkl);
BandedBilinearOperatorTest <double, std::uint64_t> mkl_banded_bilinear_operator_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
BandedBilinearOperatorTest <__float128, std::uint32_t> banded_bilinear_operator_test_float128_uint32(PreferredBackend::generic);
BandedBilinearOperatorTest <__float128, std::uint64_t> banded_bilinear_operator_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
BandedBilinearOperatorTest <Half, std::uint32_t> banded_bilinear_operator_test_half_uint32(PreferredBackend::generic);
BandedBilinearOperatorTest <Half, std::uint64_t> banded_bilinear_operator_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
BandedBilinearOperatorTest <float, std::uint32_t> cuda_banded_bilinear_operator_test_float_uint32(PreferredBackend::cuda);
BandedBilinearOperatorTest <double, std::uint32_t> cuda_banded_bilinear_operator_test_double_uint32(PreferredBackend::cuda);
BandedBilinearOperatorTest <float, std::uint64_t> cuda_banded_bilinear_operator_test_float_uint64(PreferredBackend::cuda);
BandedBilinearOperatorTest <double, std::uint64_t> cuda_banded_bilinear_operator_test_double_uint64(PreferredBackend::cuda);
#endif
