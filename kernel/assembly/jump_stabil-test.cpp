// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/cro_rav_ran_tur/element.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/assembly/trace_assembler.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/runtime.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

// reference jumps computed by FEAT2 (inner edges only)
// parameters:
//   gamma  = 0.01
//   gamma* = 0.0
//   exp    = 2.0
//   mult   = 1.0
//   nu     = 1.0
//   cub.   = 5-point Gauss-Legendre

static const double ref_jumps_q1_grad[] =
{
  0.000000000000E+00,
  2.666666666667E-02,
  6.192881254230E-03,
  9.035685301801E-04,
  1.173890312926E-04,
  1.481576972996E-05,
  1.856439277898E-06,
  2.321947419078E-07,
  2.902871429374E-08,
  3.628726240484E-09,
  4.535957446099E-10
};

static const double ref_jumps_q2_grad[] =
{
  0.000000000000E+00,
  1.164021703892E-03,
  1.073544645324E-05,
  8.723580979868E-08,
  6.881645595531E-10,
  5.389184464480E-12,
  4.144326881337E-14,
  -5.140613578791E-16,
  -2.286102417061E-15,
  -4.246260405313E-15
};

static const double ref_jumps_q1t_grad[] =
{
  0.000000000000E+00,
  3.232381016959E-02, // EL_EM30_2D
  9.774203971881E-03,
  1.601650637608E-03,
  2.208123298499E-04,
  2.872861604293E-05,
  3.655598765458E-06,
  4.607845710407E-07,
  5.783136740318E-08,
  7.243301815328E-09,
  9.063101117509E-10
};


template<typename DT_, typename IT_>
class JumpStabilTest :
  public UnitTest
{
public:
  typedef DT_ DataType;
  typedef IT_ IndexType;

  typedef LAFEM::SparseMatrixCSR<DataType, IndexType> MatrixType;
  typedef LAFEM::DenseVector<DataType, IndexType> VectorType;

  typedef Geometry::ConformalMesh<Shape::Quadrilateral> MeshType;
  typedef Geometry::RefinedUnitCubeFactory<MeshType> MeshFactory;
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;

  const DataType tol;

  JumpStabilTest(PreferredBackend backend) :
    UnitTest("JumpStabilTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend),
    tol(Math::pow(Math::eps<DataType>(), DataType(0.8)))
  {
  }

  virtual void run() const override
  {
    const IT_ level_max = IT_(5);

    // create coarse mesh
    std::shared_ptr<MeshType> mesh;
    Geometry::UnitCubeFactory<MeshType> mesh_factory;

    std::cout << "Q1: Gradient" << "\n";
    mesh = std::make_shared<MeshType>(mesh_factory);
    for(IT_ level = IT_(1); level < level_max; ++level)
    {
      {
        auto coarse_mesh = mesh;
        Geometry::StandardRefinery<MeshType> refinery(*coarse_mesh);
        mesh = std::make_shared<MeshType>(refinery);
      }

      TrafoType trafo(*mesh);
      Space::Lagrange1::Element<TrafoType> space(trafo);

      // test Q1
      test_space(space, level, DataType(level < 11 ? ref_jumps_q1_grad[level] : 0.0));
    }
    std::cout << "\n";

    std::cout << "Q2: Gradient" << "\n";
    mesh = std::make_shared<MeshType>(mesh_factory);
    for(IT_ level = IT_(1); level < level_max; ++level)
    {
      {
        auto coarse_mesh = mesh;
        Geometry::StandardRefinery<MeshType> refinery(*coarse_mesh);
        mesh = std::make_shared<MeshType>(refinery);
      }

      TrafoType trafo(*mesh);
      Space::Lagrange2::Element<TrafoType> space(trafo);

      // test Q2
      test_space(space, level, DataType(level < 10 ? ref_jumps_q2_grad[level] : 0.0));
    }
    std::cout << "\n";

    std::cout << "Q1~: Gradient" << "\n";
    mesh = std::make_shared<MeshType>(mesh_factory);
    for(IT_ level = IT_(1); level < level_max; ++level)
    {
      {
        auto coarse_mesh = mesh;
        Geometry::StandardRefinery<MeshType> refinery(*coarse_mesh);
        mesh = std::make_shared<MeshType>(refinery);
      }

      TrafoType trafo(*mesh);
      Space::CroRavRanTur::Element<TrafoType> space(trafo);

      // test Q1~
      test_space(space, level, DataType(level < 11 ? ref_jumps_q1t_grad[level] : 0.0));
    }
    std::cout << "\n";
  }

  template<typename Space_>
  void test_space(Space_& space, IT_ level, DataType ref_value) const
  {
    // assemble extended matrix structure
    MatrixType matrix;
    Assembly::SymbolicAssembler::assemble_matrix_ext_facet1(matrix, space);

    // allocate vectors
    VectorType vec_x(matrix.create_vector_l());
    VectorType vec_y(matrix.create_vector_l());

    // format
    matrix.format();
    vec_x.format();
    vec_y.format();

    // Assemble EOJ matrix
    Cubature::DynamicFactory cubature_factory("gauss-legendre:5");
    Assembly::TraceAssembler<typename Space_::TrafoType> trace_asm(space.get_trafo());
    trace_asm.compile_all_facets(true, false);
    trace_asm.assemble_jump_stabil_operator_matrix(matrix, space, cubature_factory, 0.01, 2.0, 2.0);

    // interpolate sine bubble
    Analytic::Common::SineBubbleFunction<2> sine_bubble;
    Assembly::Interpolator::project(vec_x, sine_bubble, space);

    // compute y := A*x
    matrix.apply(vec_y, vec_x);

    // compute <x,y>
    DataType jump = vec_x.dot(vec_y);

    // print
    std::cout << "Level " << stringify(level).pad_front(2) << ": "
      << stringify_fp_sci(jump, 12) << " | "
      << stringify_fp_sci(ref_value, 12) << " | "
      << stringify_fp_sci(Math::abs(ref_value-jump), 5) << "\n";

    // check
    TEST_CHECK_EQUAL_WITHIN_EPS(jump, ref_value, tol);
  }
}; // class JumpStabilTest

JumpStabilTest <double, std::uint32_t> jump_stabil_test_double_uint32(PreferredBackend::generic);
JumpStabilTest <double, std::uint64_t> jump_stabil_test_double_uint64(PreferredBackend::generic);
