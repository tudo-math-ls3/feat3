#include <test_system/test_system.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/structured_mesh.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/pointstar_structure.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/math.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

template<typename MatrixType_>
class BilinearOperatorTest :
  public TestSystem::FullTaggedTest<typename MatrixType_::MemType, Archs::None, typename MatrixType_::DataType, typename MatrixType_::IndexType>
{
  typedef typename MatrixType_::MemType MemType_;
  typedef typename MatrixType_::DataType DataType_;
  typedef typename MatrixType_::IndexType IndexType_;
  typedef LAFEM::DenseVector<MemType_, DataType_, IndexType_> VectorType;

  typedef Geometry::ConformalMesh<Shape::Quadrilateral, 2, 2, DataType_> QuadMesh;

  typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

  typedef Space::Discontinuous::Element<QuadTrafo> QuadSpaceQ0;
  typedef Space::Lagrange1::Element<QuadTrafo> QuadSpaceQ1;

public:
  BilinearOperatorTest() :
    TestSystem::FullTaggedTest<typename MatrixType_::MemType, Archs::None, DataType_, IndexType_>("BilinearOperatorTest<" + MatrixType_::name() + ">")
  {
  }

  virtual void run() const
  {
    test_unit_2d();
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
    MatrixType_ matrix;
    Assembly::SymbolicMatrixAssembler<>::assemble1(matrix, space);
    matrix.format();

    // create a cubature factory
    Cubature::DynamicFactory cubature_factory("barycentre");

    // assemble the identity operator with barycentre cubature rule
    Assembly::Common::IdentityOperator operat;
    Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix, operat, space, cubature_factory);

    // create CSR-matrix of matrix
    LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_> tmp_matrix;
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
    MatrixType_ matrix_1, matrix_2;
    Assembly::SymbolicMatrixAssembler<>::assemble1(matrix_1, space);
    Assembly::SymbolicMatrixAssembler<>::assemble1(matrix_2, space);
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
    LAFEM::DenseVector<Mem::Main, Index> qatv(num_verts, Index(0));
    for(Index i(0); i < num_quads; ++i)
    {
      for(Index j(0); j < 4; ++j)
        qatv(vatq(i,j), qatv(vatq(i,j)) + 1);
    }

    // create a dof-mapping
    typedef typename QuadSpaceQ1::DofMappingType DofMapping;
    typename QuadSpaceQ1::DofMappingType dof_mapping(space);

    // create local matrix data
    Assembly::LocalMatrixData<Tiny::Matrix<DataType_,4,4>, DofMapping, DofMapping>
      lmd1(dof_mapping,dof_mapping), lmd2(dof_mapping,dof_mapping);
    LAFEM::GatherAxpy<MatrixType_> gather1(matrix_1), gather2(matrix_2);

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
      gather1(lmd1);
      gather2(lmd2);

      // loop over all 4x4 entries
      for(Index i(0); i < 4; ++i)
      {
        Index nvi = qatv(vatq(cell,i));
        for(Index j(0); j < 4; ++j)
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

BilinearOperatorTest<LAFEM::SparseMatrixCSR<Mem::Main, float, unsigned int> > bilinear_operator_test_csr_float_uint;
BilinearOperatorTest<LAFEM::SparseMatrixCSR<Mem::Main, float, unsigned long> > bilinear_operator_test_csr_float_ulong;
BilinearOperatorTest<LAFEM::SparseMatrixCSR<Mem::Main, double, unsigned int> > bilinear_operator_test_csr_double_uint;
BilinearOperatorTest<LAFEM::SparseMatrixCSR<Mem::Main, double, unsigned long> > bilinear_operator_test_csr_double_ulong;

BilinearOperatorTest<LAFEM::SparseMatrixCOO<Mem::Main, float, unsigned int> > bilinear_operator_test_coo_float_uint;
BilinearOperatorTest<LAFEM::SparseMatrixCOO<Mem::Main, float, unsigned long> > bilinear_operator_test_coo_float_ulong;
BilinearOperatorTest<LAFEM::SparseMatrixCOO<Mem::Main, double, unsigned int> > bilinear_operator_test_coo_double_uint;
BilinearOperatorTest<LAFEM::SparseMatrixCOO<Mem::Main, double, unsigned long> > bilinear_operator_test_coo_double_ulong;

BilinearOperatorTest<LAFEM::SparseMatrixELL<Mem::Main, float, unsigned int> > bilinear_operator_test_ell_float_uint;
BilinearOperatorTest<LAFEM::SparseMatrixELL<Mem::Main, float, unsigned long> > bilinear_operator_test_ell_float_ulong;
BilinearOperatorTest<LAFEM::SparseMatrixELL<Mem::Main, double, unsigned int> > bilinear_operator_test_ell_double_uint;
BilinearOperatorTest<LAFEM::SparseMatrixELL<Mem::Main, double, unsigned long> > bilinear_operator_test_ell_double_ulong;

#ifdef FEAST_HAVE_QUADMATH
BilinearOperatorTest<LAFEM::SparseMatrixCSR<Mem::Main, __float128, unsigned int> > bilinear_operator_test_csr_float128_uint;
BilinearOperatorTest<LAFEM::SparseMatrixCSR<Mem::Main, __float128, unsigned long> > bilinear_operator_test_csr_float128_ulong;

BilinearOperatorTest<LAFEM::SparseMatrixCOO<Mem::Main, __float128, unsigned int> > bilinear_operator_test_coo_float128_uint;
BilinearOperatorTest<LAFEM::SparseMatrixCOO<Mem::Main, __float128, unsigned long> > bilinear_operator_test_coo_float128_ulong;

BilinearOperatorTest<LAFEM::SparseMatrixELL<Mem::Main, __float128, unsigned int> > bilinear_operator_test_ell_float128_uint;
BilinearOperatorTest<LAFEM::SparseMatrixELL<Mem::Main, __float128, unsigned long> > bilinear_operator_test_ell_float128_ulong;
#endif


template<typename MemType_, typename DataType_, typename IndexType_>
class BandedBilinearOperatorTest :
  public TestSystem::FullTaggedTest<MemType_, Archs::None, DataType_, IndexType_>
{
  typedef LAFEM::SparseMatrixBanded<Mem::Main, DataType_, IndexType_> MatrixType;
  typedef LAFEM::DenseVector<MemType_, DataType_, IndexType_> VectorType;

  typedef Geometry::StructuredMesh<2, 2, 2, DataType_> QuadMesh;

  typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

  typedef Space::Discontinuous::Element<QuadTrafo> QuadSpaceQ0;
  typedef Space::Lagrange1::Element<QuadTrafo> QuadSpaceQ1;

public:
  BandedBilinearOperatorTest() :
    TestSystem::FullTaggedTest<MemType_, Archs::None, DataType_, IndexType_>("BandedBilinearOperatorTest")
  {
  }

  virtual void run() const
  {
    test_unit_2d();
  }

  void test_unit_2d() const
  {
    // create an 4x4 structured mesh
    const Index num_slices[2] = {Index(4), Index(4)};
    QuadMesh mesh(num_slices);

    // initialise our 5x5 vertices
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
    //std::vector<IndexType_> nsi;
    //nsi.push_back(mesh.get_num_slices(0));
    //nsi.push_back(mesh.get_num_slices(1));
    //MatrixType matrix(LAFEM::PointstarStructureFE<Algo::Generic>::template value<DataType_,IndexType_>(Index(0), nsi));
    Index neq = space.get_num_dofs();
    LAFEM::DenseVector<MemType_,DataType_, IndexType_> vv(neq, DataType_(0));
    LAFEM::DenseVector<MemType_, IndexType_, IndexType_> vi(Index(1), IndexType_(neq) - IndexType_(1));
    MatrixType matrix(neq, neq, vv, vi);

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
    //std::vector<IndexType_> nsi;
    //nsi.push_back(mesh.get_num_slices(0));
    //nsi.push_back(mesh.get_num_slices(1));
    //MatrixType matrix_1(LAFEM::PointstarStructureFE<Algo::Generic>::template value<DataType_,IndexType_>(Index(1), nsi));
    Index nsl = mesh.get_num_slices(0) + Index(1);
    Index neq = space.get_num_dofs();
    Index nze = Index(9)*neq;
    LAFEM::DenseVector<MemType_,DataType_, IndexType_> vv(nze, DataType_(0));
    LAFEM::DenseVector<MemType_, IndexType_, IndexType_> vi(Index(9));
    vi(Index(0), IndexType_(neq - nsl) - IndexType_(2));
    vi(Index(1), IndexType_(neq - nsl) - IndexType_(1));
    vi(Index(2), IndexType_(neq - nsl) - IndexType_(0));
    vi(Index(3), IndexType_(neq      ) - IndexType_(2));
    vi(Index(4), IndexType_(neq      ) - IndexType_(1));
    vi(Index(5), IndexType_(neq      ) - IndexType_(0));
    vi(Index(6), IndexType_(neq + nsl) - IndexType_(2));
    vi(Index(7), IndexType_(neq + nsl) - IndexType_(1));
    vi(Index(8), IndexType_(neq + nsl) - IndexType_(0));
    MatrixType matrix_1(neq, neq, vv, vi);

    MatrixType matrix_2(matrix_1.clone());
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

    const auto& vatq(mesh.template get_index_set<2,0>());

    // create a temporary array to count the number of quads adjacent to a vertex
    LAFEM::DenseVector<Mem::Main, Index> qatv(num_verts, Index(0));
    for(Index i(0); i < num_quads; ++i)
    {
      for(Index j(0); j < 4; ++j)
        qatv(vatq(i,j), qatv(vatq(i,j)) + 1);
    }

    // create a dof-mapping
    typedef typename QuadSpaceQ1::DofMappingType DofMapping;
    typename QuadSpaceQ1::DofMappingType dof_mapping(space);

    // create local matrix data
    Assembly::LocalMatrixData<Tiny::Matrix<DataType_,4,4>, DofMapping, DofMapping>
      lmd1(dof_mapping,dof_mapping), lmd2(dof_mapping,dof_mapping);
    LAFEM::GatherAxpy<MatrixType> gather1(matrix_1), gather2(matrix_2);

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
      gather1(lmd1);
      gather2(lmd2);

      // loop over all 4x4 entries
      for(Index i(0); i < 4; ++i)
      {
        Index nvi = qatv(vatq(cell,i));
        for(Index j(0); j < 4; ++j)
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

BandedBilinearOperatorTest<Mem::Main, float, unsigned int> banded_bilinear_operator_test_float_uint;
BandedBilinearOperatorTest<Mem::Main, float, unsigned long> banded_bilinear_operator_test_float_ulong;
BandedBilinearOperatorTest<Mem::Main, double, unsigned int> banded_bilinear_operator_test_double_uint;
BandedBilinearOperatorTest<Mem::Main, double, unsigned long> banded_bilinear_operator_test_double_ulong;

#ifdef FEAST_HAVE_QUADMATH
BandedBilinearOperatorTest<Mem::Main, __float128, unsigned int> banded_bilinear_operator_test_float128_uint;
BandedBilinearOperatorTest<Mem::Main, __float128, unsigned long> banded_bilinear_operator_test_float128_ulong;
#endif
