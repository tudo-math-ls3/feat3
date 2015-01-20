#include <test_system/test_system.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/math.hpp>
#ifdef FEAST_HAVE_ALGLIB
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/reference_cell_factory.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_smoother.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_smoother_q1hack.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_functional_2d_q1.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_functional_2d_q1_d2.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_functional_2d_q1hack.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_functional_2d_p1.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_functional_2d_p1_d2.hpp>

using namespace FEAST;
using namespace FEAST::TestSystem;

/// \brief Helper class for resizing tests
template<typename ShapeType_>
struct helperclass;

/// \brief Specialisation for hypercubes
template<int shape_dim_>
struct helperclass< FEAST::Shape::Hypercube<shape_dim_> >
{
  /// \brief Sets coordinates so we deal the the reference element
  template<typename VectorType_, typename DataType_>
  static void set_coords(VectorType_& coords_, const DataType_& scaling)
  {
    for(Index i(0); i < Index(1 << shape_dim_); ++i)
    {
      for(Index d(0); d < Index(shape_dim_); ++d)
        coords_[d](i, (DataType_(((i >> d) & 1) << 1) - DataType_(1)) * scaling );
    }
  }

};

/// \brief Specialisation for 2d simplices
template<>
struct helperclass< FEAST::Shape::Simplex<2> >
{
  /// \brief Sets coordinates so we deal the the Rumpf reference element
  template<typename VectorType_, typename DataType_>
  static void set_coords(VectorType_& coords_, const DataType_& scaling)
  {
    coords_[0](0, DataType_(0));
    coords_[1](0, DataType_(0));

    coords_[0](1, DataType_(1) * scaling);
    coords_[1](1, DataType_(0));

    coords_[0](2, DataType_(0.5) * scaling);
    coords_[1](2, DataType_(0.5)*Math::sqrt(DataType_(3))*scaling);
  }
};

/**
 * \brief Test for Rumpf smoothers and functionals
 *
 * When the input mesh consists of a single Rumpf reference element, the mesh optimiser is supposed not to change it,
 * as:
 * 1) The functional value is minimal for the reference element and
 * 2) the functional gradient evaluates to zero meaning that
 * 3) the optimiser stops after one evaluation of the gradient.
 *
 * This means the functional value pre and post optimisation should be the same.
 *
 * \author Jordi Paul
 **/
template
<
  typename DataType_,
  typename MemType_,
  typename ShapeType_,
  template<typename, typename> class FunctionalType_,
  template<typename ... > class RumpfSmootherType_
>
class RumpfSmootherTest_2d
: public TestSystem::TaggedTest<Archs::None, DataType_>
{
  public:
    typedef MemType_ MemType;
    typedef DataType_ DataType;

    typedef ShapeType_ ShapeType;
    typedef Geometry::ConformalMesh<ShapeType, ShapeType::dimension,ShapeType::dimension, DataType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;

    typedef FunctionalType_<DataType, ShapeType> FunctionalType;
    typedef RumpfSmootherType_<DataType, MemType, TrafoType, FunctionalType> RumpfSmootherType;

    RumpfSmootherTest_2d() :
      TestSystem::TaggedTest<Archs::None, DataType_>("rumpf_smoother_test")
      {
      }

    virtual void run() const
    {

      // Mesh and trafo
      Geometry::ReferenceCellFactory<ShapeType, DataType> mesh_factory;
      MeshType mesh(mesh_factory);
      TrafoType trafo(mesh);

      // In 2d, the cofactor matrix is not used
      DataType fac_norm = DataType(1e0),fac_det = DataType(1e0),fac_cof = DataType(0), fac_reg(DataType(1e-8));
      FunctionalType my_functional(fac_norm, fac_det, fac_cof, fac_reg);

      // The smoother in all its template glory
      RumpfSmootherType rumpflpumpfl(trafo, my_functional);

      // Set bdry_id to 0 again so the mesh smoother can resize the single element
      for(Index i(0); i < mesh.get_num_entities(0); ++i)
        rumpflpumpfl._bdry_id[i] = 0;

      // It is possible to scale the reference element, but we do not want that here
      DataType scaling(1);
      // This transforms the unit element to the Rumpf reference element
      helperclass<ShapeType>::set_coords(rumpflpumpfl._coords, scaling);

      // Since we changed the internal _coords, they have to be copied back to the mesh
      rumpflpumpfl.set_coords();
      rumpflpumpfl.init();


      // Compute initial functional value
      DataType fval_pre = rumpflpumpfl.compute_functional();
      // Optimise the mesh
      rumpflpumpfl.optimise();
      // Compute new functional value
      DataType fval_post = rumpflpumpfl.compute_functional();

      const DataType eps = Math::eps<DataType_>();
      TEST_CHECK_EQUAL_WITHIN_EPS(fval_pre, fval_post, eps);

    }
};

template<typename A, typename B, typename C, typename D>
using MySmoother = Geometry::RumpfSmoother<A, B, C, D>;

template<typename A, typename B, typename C, typename D>
using MySmootherQ1Hack = Geometry::RumpfSmootherQ1Hack<A, B, C, D>;

typedef Mem::Main MemType;

RumpfSmootherTest_2d<float, MemType, Shape::Hypercube<2>, Geometry::RumpfFunctional, MySmoother> test_hc_f_1;
RumpfSmootherTest_2d<float, MemType, Shape::Hypercube<2>, Geometry::RumpfFunctional_D2, MySmoother> test_hc_f_2;
RumpfSmootherTest_2d<float, MemType, Shape::Simplex<2>, Geometry::RumpfFunctional, MySmoother> test_s_f_1;
RumpfSmootherTest_2d<float, MemType, Shape::Simplex<2>, Geometry::RumpfFunctional_D2, MySmoother> test_s_f_2;

RumpfSmootherTest_2d<double, MemType, Shape::Hypercube<2>, Geometry::RumpfFunctional, MySmoother> test_hc_d_1;
RumpfSmootherTest_2d<double, MemType, Shape::Hypercube<2>, Geometry::RumpfFunctional_D2, MySmoother> test_hc_d_2;
RumpfSmootherTest_2d<double, MemType, Shape::Simplex<2>, Geometry::RumpfFunctional, MySmoother> test_s_d_1;
RumpfSmootherTest_2d<double, MemType, Shape::Simplex<2>, Geometry::RumpfFunctional_D2, MySmoother> test_s_d_2;

template<typename A, typename B>
using MyFunctionalQ1Hack = Geometry::RumpfFunctionalQ1Hack<A, B, Geometry::RumpfFunctional>;

template<typename A, typename B>
using MyFunctionalQ1Hack_D2 = Geometry::RumpfFunctionalQ1Hack<A, B, Geometry::RumpfFunctional_D2>;

RumpfSmootherTest_2d<float, MemType, Shape::Hypercube<2>, MyFunctionalQ1Hack, MySmootherQ1Hack> test_q1hack_f_1;
RumpfSmootherTest_2d<float, MemType, Shape::Hypercube<2>, MyFunctionalQ1Hack_D2, MySmootherQ1Hack> test_q1hack_f_2;
RumpfSmootherTest_2d<double, MemType, Shape::Hypercube<2>, MyFunctionalQ1Hack, MySmootherQ1Hack> test_q1hack_d_1;
RumpfSmootherTest_2d<double, MemType, Shape::Hypercube<2>, MyFunctionalQ1Hack_D2, MySmootherQ1Hack> test_q1hack_d_2;

#endif
