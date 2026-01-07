// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.


#include <kernel/base_header.hpp>
#include <kernel/runtime.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/analytic/lambda_function.hpp>
#include <test_system/test_system.hpp>

#include <kernel/assembly/surface_integrator.hpp>
#include <kernel/assembly/surface_integrator_basic_jobs.hpp>



using namespace FEAT;
using namespace FEAT::TestSystem;



template<typename DataType_, typename IndexType_>
class SurfaceIntegratorTest :
  public UnitTest
{
public:
  typedef DataType_ DataType;
  typedef IndexType_ IndexType;

public:
  SurfaceIntegratorTest(PreferredBackend backend) :
    UnitTest("SurfaceIntegratorTest", Type::Traits<DataType_>::name(), Type::Traits<IndexType_>::name(), backend)
  {
  }

  virtual ~SurfaceIntegratorTest()
  {
  }

  template<typename DT_>
  static DT_ calc_simple_analytic_simplex_3d(const Tiny::Vector<DT_, 3>& v0, const Tiny::Vector<DT_, 3>& v1, const Tiny::Vector<DT_, 3>& v2)
  {
    Tiny::Matrix<DT_, 3, 2> trafo;
    for(int d = 0; d < 3; ++d)
    {
      trafo[d][0] = v1[d] - v0[d];
      trafo[d][1] = v2[d] - v0[d];
    }

    auto normal = Tiny::orthogonal(trafo).normalize().negate();
    Tiny::Matrix<DT_, 2, 3> trans;
    trans.set_transpose(trafo);
    auto gram = trans * trafo;
    auto vol = Math::sqrt(gram.det());


    DT_ integral = vol * (DT_(1)/DT_(3)*normal[0]*(trafo[0][0] + trafo[0][1] + DT_(3) * v0[0])
                  + DT_(1)/DT_(6)*normal[1]*(trafo[1][0] + trafo[1][1] + DT_(3) * v0[1])
                  - DT_(1)/DT_(2)*normal[2] *(trafo[2][0] + trafo[2][1] + DT_(3) * v0[2]));
    return integral;
  }

  virtual void run() const override
  {
    test_unit_3d();
  }

  void test_unit_3d() const
  {
    const DataType_ eps = TestSystem::tol<DataType_>();
    constexpr int dim = 3;
    typedef Shape::Hypercube<dim> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceType;
    constexpr int lvl = 3;
    Geometry::RefinedUnitCubeFactory<MeshType> unitcube_factory(lvl);
    MeshType mesh(unitcube_factory);

    TrafoType trafo(mesh);
    SpaceType space(trafo);

    LAFEM::DenseVectorBlocked<DataType, IndexType, dim> primal_vec(space.get_num_dofs());
    Analytic::SimplifiedLambdaVectorFunction3D simple_lambda([](const auto& point)
        {return Tiny::Vector<DataType, dim>{
          DataType(2)*point[0], DataType(1)*point[1], -DataType(3)*point[2]
        };});

    Assembly::Interpolator::project(primal_vec, simple_lambda, space);

    std::vector<std::array<DataType, dim>> vtx{
      {DataType(0.1), DataType(0.1), DataType(0.1)},
      {DataType(0.2), DataType(0.18), DataType(0.1)},
      {DataType(0.15), DataType(0.34), DataType(0.2)},
      {DataType(0.3), DataType(0.56), DataType(0.15)}
    };
    std::vector<std::array<IndexType, 3>> indx{
      {IndexType(0), IndexType(1), IndexType(2)}, {IndexType(1), IndexType(3), IndexType(2)}
    };

    const String cubature_string = "auto-degree:3";
    Cubature::DynamicFactory cubature(cubature_string);

    const Index surface_num = indx.size();

    LAFEM::DenseVector<DataType, IndexType> vec_face(surface_num);

    Assembly::SurfaceIntegrator<TrafoType, Shape::Simplex<2>> surf_int(trafo, cubature);
    surf_int.set_vertices(vtx);
    surf_int.set_normal_view(true);
    surf_int.add_face(indx[0]);
    surf_int.add_face(indx[1]);
    surf_int.compile();

    Assembly::NormalValueSurfaceIntegratorJob surface_int_job(vec_face, primal_vec, space, DataType(1));

    surf_int.assemble(surface_int_job);

    for(Index k = 0; k < surface_num; ++k)
    {
      DataType ref_res = calc_simple_analytic_simplex_3d(Tiny::Vector<DataType, 3>{vtx[indx[k][0]][0], vtx[indx[k][0]][1], vtx[indx[k][0]][2]},
        Tiny::Vector<DataType, 3>{vtx[indx[k][1]][0], vtx[indx[k][1]][1], vtx[indx[k][1]][2]},
        Tiny::Vector<DataType, 3>{vtx[indx[k][2]][0], vtx[indx[k][2]][1], vtx[indx[k][2]][2]});
      TEST_CHECK_EQUAL_WITHIN_EPS(vec_face(k), ref_res, eps);
    }

  }

};

SurfaceIntegratorTest<float, std::uint32_t> surface_integrator_test_float_uint32(PreferredBackend::generic);
SurfaceIntegratorTest<double, std::uint32_t> surface_integrator_test_double_uint32(PreferredBackend::generic);
SurfaceIntegratorTest<float, std::uint64_t> surface_integrator_test_float_uint64(PreferredBackend::generic);
SurfaceIntegratorTest<double, std::uint64_t> surface_integrator_test_double_uint64(PreferredBackend::generic);
