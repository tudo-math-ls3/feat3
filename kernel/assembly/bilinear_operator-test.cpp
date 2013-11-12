#include <test_system/test_system.hpp>
#include <kernel/assembly/standard_operators.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/space/dof_adjacency.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/math.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

template<typename DataType_>
class BilinearOperatorTest :
  public TestSystem::TaggedTest<Archs::None, DataType_>
{
  typedef LAFEM::DenseVector<Mem::Main, DataType_> VectorType;
  typedef LAFEM::SparseMatrixCSR<Mem::Main, DataType_> MatrixType;

  typedef Geometry::ConformalMesh<Shape::Quadrilateral> QuadMesh;

  typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

  typedef Space::Discontinuous::Element<QuadTrafo> QuadSpaceQ0;
  typedef Space::Lagrange1::Element<QuadTrafo> QuadSpaceQ1;

public:
  BilinearOperatorTest() :
    TestSystem::TaggedTest<Archs::None, DataType_>("BilinearOperatorTest")
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

  void test_unit_2d_q0(const QuadMesh& mesh) const
  {
    // compute eps
    const DataType_ eps = Math::pow(Math::Limits<DataType_>::epsilon(), DataType_(0.9));

    // create trafo
    QuadTrafo trafo(mesh);

    // create space
    QuadSpaceQ0 space(trafo);

    // create a CSR matrices
    Adjacency::Graph dof_adjacency(Space::DofAdjacency<>::assemble(space));
    MatrixType matrix(dof_adjacency);
    matrix.clear();

    // create a cubature factory
    Cubature::DynamicFactory cubature_factory("barycentre");

    // assemble the identity operator with barycentre cubature rule
    Assembly::BilinearScalarIdentityOperator operat;
    Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix, operat, space, cubature_factory);

    // fetch the matrix arrays
    DataType_* data = matrix.val();

    const DataType_ v = DataType_(1) / DataType_(mesh.get_num_entities(2));

    // loop over all matrix rows
    for(Index i(0); i < matrix.rows(); ++i)
    {
      // validate entry
      TEST_CHECK_EQUAL_WITHIN_EPS(data[i], v, eps);
    }
  }

  void test_unit_2d_q1(const QuadMesh& mesh) const
  {
    // compute eps
    const DataType_ eps = Math::pow(Math::Limits<DataType_>::epsilon(), DataType_(0.9));

    // create trafo
    QuadTrafo trafo(mesh);

    // create space
    QuadSpaceQ1 space(trafo);

    // create two CSR matrices
    Adjacency::Graph dof_adjacency(Space::DofAdjacency<>::assemble(space));
    MatrixType matrix_1(dof_adjacency);
    MatrixType matrix_2(dof_adjacency);
    matrix_1.clear();
    matrix_2.clear();

    // create a cubature factory
    Cubature::DynamicFactory cubature_factory_trz("trapezoidal");
    Cubature::DynamicFactory cubature_factory_gl2("gauss-legendre:2");

    // assemble the identity operator with barycentre cubature rule
    Assembly::BilinearScalarLaplaceOperator operat;
    Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix_1, operat, space, cubature_factory_trz);
    Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix_2, operat, space, cubature_factory_gl2);

    // get mesh element count
    Index num_verts = mesh.get_num_entities(0);
    Index num_quads = mesh.get_num_entities(2);

    const Geometry::IndexSet<4>& vatq(mesh.get_index_set<2,0>());

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
      lmd1.clear();
      lmd2.clear();
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

BilinearOperatorTest<float> bilinear_operator_test_float;
BilinearOperatorTest<double> bilinear_operator_test_double;
