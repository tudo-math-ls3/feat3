#include <iostream>
#include <kernel/archs.hpp>
#include <kernel/util/math.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_smoother_levelset.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_functional_levelset_2d_p1.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_functional_levelset_2d_q1.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/assembly/common_functions.hpp>
#include <kernel/assembly/discrete_projector.hpp>

using namespace FEAST;

namespace FEAST
{
  namespace Assembly
  {
    namespace Common
    {
      /**
       * @brief Class Template for the gradient of DistanceFunctionSD
       *
       * \tparam ImgPointType_
       * Type for the point in the DistanceFunction
       *
       * \tparam component
       * Component of the gradient
       *
       **/
      template<typename ImgPointType_, int component>
      class DistanceFunctionSD_grad :
        public AnalyticFunction
      {
        public:
          /// Datatype from the point
          typedef typename ImgPointType_::DataType DataType;
          /** \copydoc AnalyticFunction::FunctionCapabilities */
          enum FunctionCapabilities
          {
            can_value = 1,
            can_grad = 0,
            can_hess = 0
          };

          /** \copydoc AnalyticFunction::ConfigTraits */
          template<typename Config_>
          struct ConfigTraits
          {
            /**
             * \brief Trafo configuration tag class
             *
             * \see Trafo::ConfigBase
             */
            struct TrafoConfig :
              public Trafo::ConfigBase
            {
              enum
              {
                need_img_point = 1
              };
            };
          };

          /** \copydoc AnalyticFunction::Evaluator */
          template<typename EvalTraits_>
          class Evaluator :
            public AnalyticFunction::Evaluator<EvalTraits_>
        {
          public:
            /// trafo evaluator data
            typedef typename EvalTraits_::TrafoEvaluator TrafoEvaluator;
            /// trafo data type
            typedef typename EvalTraits_::TrafoData TrafoData;
            /// coefficient data type
            typedef typename EvalTraits_::DataType DataType;
            /// value type
            typedef typename EvalTraits_::ValueType ValueType;
            /// gradient type
            typedef typename EvalTraits_::GradientType GradientType;
            /// hessian type
            typedef typename EvalTraits_::HessianType HessianType;
            /// type for the points the analytic function is evaluated at
            typedef typename EvalTraits_::ImagePointType ImgPointType;

          private:
            /// Function to evaluate
            const DistanceFunctionSD_grad& _function;

          public:
            /// Constructor
            explicit Evaluator(const DistanceFunctionSD_grad& function) :
              _function(function)
              {
              }

            ValueType value(const TrafoData& tau) const
            {
              ImgPointType tmp (tau.img_point - _function._point);
              DataType norm = tmp.norm_euclid();
              if(norm <= Math::eps<DataType>())
              {
                return DataType(0);
              }
              else
              {
                return _function._b*( tau.img_point[component] - _function._point(component))/norm;
              }
            }
        }; // class DistanceFunctionSD::Evaluator<...>

        private:
          /// The point to which the distance to is calculated
          ImgPointType_& _point;
          /// Displacement of the function
          DataType _a;
          /// Scaling factor
          DataType _b;

        public:
          /// Constructor
          explicit DistanceFunctionSD_grad(ImgPointType_& x0_, const DataType a_, const DataType b_) :
            _point(x0_),
            _a(a_),
            _b(b_)
            {
            }

      }; // class DistanceFunctionSD_grad

    } // namespace Common
  } // namespace Assembly
} // namespace FEAST

/**
 * @brief Runs mesh smoother stuff
 *
 **/
  template<typename DataType_, typename ShapeType_ >
  void run()
  {
    typedef DataType_ DataType;
    typedef Mem::Main MemType;
    typedef ShapeType_ ShapeType;

    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange1::Element<TrafoType> SpaceType;
    typedef Tiny::Vector<DataType,2,2> ImgPointType;
    typedef Assembly::Common::DistanceFunctionSD<ImgPointType> AnalyticFunctionType;
    typedef Assembly::Common::DistanceFunctionSD_grad<ImgPointType, 0> AnalyticFunctionGrad0Type;
    typedef Assembly::Common::DistanceFunctionSD_grad<ImgPointType, 1> AnalyticFunctionGrad1Type;
    //    typedef typename Geometry::RumpfFunctionalLevelset<MemType, DataType, TrafoType> FunctionalType;
    typedef typename Geometry::RumpfFunctionalLevelset<MemType, DataType, TrafoType> FunctionalType;

    typedef LAFEM::DenseVector<MemType, DataType> VectorType;


    // Mesh and trafo
    Index level(4);
    Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(level);
    MeshType mesh(mesh_factory);
    TrafoType trafo(mesh);

    // Levelset function stuff
    ImgPointType x0(DataType(0));
    x0.v[0] = DataType(0.5);
    x0.v[1] = DataType(0.85);
    // Levelset function is phi(x) = displacement + scaling * || x - x0 ||_2
    DataType displacement(0.5);
    DataType scaling(-4);
    AnalyticFunctionType analytic_lvlset(x0, displacement, scaling);
    AnalyticFunctionGrad0Type analytic_lvlset_grad0(x0, displacement, scaling);
    AnalyticFunctionGrad1Type analytic_lvlset_grad1(x0, displacement, scaling);

    SpaceType lvlsetspace(trafo);
    VectorType lvlset_vec, lvlset_vec_out;
    Assembly::Interpolator::project(lvlset_vec, analytic_lvlset, lvlsetspace);
    FEAST::Assembly::DiscreteVertexProjector::project(lvlset_vec_out, lvlset_vec, lvlsetspace);

    DataType pi(Math::pi<DataType>());

    DataType deltat(DataType(0.0125));

    DataType fac_norm = DataType(5e-1),fac_det = DataType(1e0),fac_cof = DataType(0), fac_reg(DataType(1e-8));
    bool from_original(false);
    bool align_to_lvlset(false);
    bool r_adaptivity(true);
    FunctionalType my_functional(fac_norm, fac_det, fac_cof, fac_reg);

    // The smoother in all its template glory
    Geometry::RumpfSmootherLevelsetAnalytic
    <
      AnalyticFunctionType,
      AnalyticFunctionGrad0Type,
      AnalyticFunctionGrad1Type,
      SpaceType,
      FunctionalType,
      TrafoType,
      DataType,
      MemType
    > rumpflpumpfl(trafo, my_functional, align_to_lvlset, from_original, r_adaptivity, analytic_lvlset, analytic_lvlset_grad0, analytic_lvlset_grad1);
    rumpflpumpfl.init();

    DataType* func_norm(new DataType[mesh.get_num_entities(2)]);
    DataType* func_det(new DataType[mesh.get_num_entities(2)]);
    DataType* func_det2(new DataType[mesh.get_num_entities(2)]);
    DataType* func_lvlset(new DataType[mesh.get_num_entities(2)]);

    // Evaluates the levelset function and its gradient
    rumpflpumpfl.prepare();
    // Compute initial functional value
    DataType fval(0);
    fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_det2, func_lvlset);
    std::cout << "fval pre optimisation = " << scientify(fval) << std::endl;
    // Compute initial functional gradient
    rumpflpumpfl.compute_gradient();

    // Filename for vtk files
    std::string filename;

    Geometry::ExportVTK<MeshType> writer_pre(mesh);
    writer_pre.add_scalar_cell("lambda", rumpflpumpfl._lambda.elements() );
    writer_pre.add_scalar_cell("h", rumpflpumpfl._h[0].elements() );
    writer_pre.add_scalar_vertex("grad_0", &rumpflpumpfl._grad[0]);
    writer_pre.add_scalar_vertex("grad_1", &rumpflpumpfl._grad[mesh.get_num_entities(0)]);
    writer_pre.add_scalar_vertex("levelset", rumpflpumpfl._lvlset_vtx_vec.elements());
    writer_pre.add_scalar_vertex("lvlset_grad_0", rumpflpumpfl._lvlset_grad_vtx_vec[0].elements());
    writer_pre.add_scalar_vertex("lvlset_grad_1", rumpflpumpfl._lvlset_grad_vtx_vec[1].elements());
    writer_pre.add_scalar_cell("levelset_constraint", func_lvlset );
    writer_pre.write("pre_initial.vtk");

    rumpflpumpfl.optimise();

    fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_det2, func_lvlset);
    std::cout << "fval post optimisation = " << scientify(fval) << std::endl;

    Geometry::ExportVTK<MeshType> writer_post(mesh);
    writer_post.add_scalar_cell("lambda", rumpflpumpfl._lambda.elements() );
    writer_post.add_scalar_cell("h", rumpflpumpfl._h[0].elements() );
    writer_post.add_scalar_vertex("grad_0", &rumpflpumpfl._grad[0]);
    writer_post.add_scalar_vertex("grad_1", &rumpflpumpfl._grad[mesh.get_num_entities(0)]);
    writer_post.add_scalar_vertex("levelset", rumpflpumpfl._lvlset_vtx_vec.elements());
    writer_post.add_scalar_vertex("lvlset_grad_0", rumpflpumpfl._lvlset_grad_vtx_vec[0].elements());
    writer_post.add_scalar_vertex("lvlset_grad_1", rumpflpumpfl._lvlset_grad_vtx_vec[1].elements());
    writer_post.add_scalar_cell("levelset_constraint", func_lvlset );
    writer_post.write("post_initial.vtk");

    DataType time(0);
    Index n(0);

    DataType* mesh_velocity(new DataType[mesh.get_num_entities(0)]);

    // Old mesh coordinates for computing the mesh velocity
    LAFEM::DenseVector<MemType, DataType_> coords_old[MeshType::world_dim];
    for(Index d = 0; d < MeshType::world_dim; ++d)
      coords_old[d]= std::move(LAFEM::DenseVector<MemType, DataType>(mesh.get_num_entities(0)));

    Index outputstep(1);
    deltat /= DataType(outputstep);
    std::cout << "deltat = " << scientify(deltat) << ", outputstep = " << outputstep << std::endl;
    std::cout << "fac_norm = " << scientify(fac_norm) << ", fac_det= " << scientify(fac_det)
    << ", fac_reg = " << scientify(fac_reg) << std::endl;

    while(time < DataType(6))
    {

      std::cout << "timestep " << n << std::endl;
      time+= deltat;

      // Save old vertex coordinates
      for(Index d(0); d < MeshType::world_dim; ++d)
      {
        for(Index i(0); i < mesh.get_num_entities(0); ++i)
          coords_old[d](i, rumpflpumpfl._coords[d](i));
      }

      // update leveset function
      x0.v[0] = DataType(0.5) + 0.35*Math::sin(DataType(2)*pi*time);
      x0.v[1] = DataType(0.5) + 0.35*Math::cos(DataType(2)*pi*time);
      analytic_lvlset.set_point(x0);

      if( n%outputstep || outputstep==1)
      {
        fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_det2, func_lvlset);
        rumpflpumpfl.compute_gradient();
        std::cout << "fval post optimisation = " << scientify(fval) << std::endl;

        filename = "pre_" + stringify(n) + ".vtk";
        Geometry::ExportVTK<MeshType> writer_pre(mesh);

        writer_pre.add_scalar_cell("lambda", rumpflpumpfl._lambda.elements() );
        writer_pre.add_scalar_cell("h", rumpflpumpfl._h[0].elements() );
        writer_pre.add_scalar_vertex("grad_0", &rumpflpumpfl._grad[0]);
        writer_pre.add_scalar_vertex("grad_1", &rumpflpumpfl._grad[mesh.get_num_entities(0)]);
        writer_pre.add_scalar_cell("norm", func_norm);
        writer_pre.add_scalar_cell("det", func_det);
        writer_pre.add_scalar_cell("det2", func_det2);
        writer_pre.add_scalar_vertex("grad_0", &rumpflpumpfl._grad[0]);
        writer_pre.add_scalar_vertex("grad_1", &rumpflpumpfl._grad[mesh.get_num_entities(0)]);
        writer_pre.add_scalar_vertex("levelset", rumpflpumpfl._lvlset_vtx_vec.elements());
        writer_pre.add_scalar_cell("levelset_constraint", func_lvlset );
        //std::cout << "Writing " << filename << std::endl;
        writer_pre.write(filename);
      }

      rumpflpumpfl.optimise();

      DataType max_mesh_velocity(-1e10);
      DataType ideltat = DataType(1)/deltat;
      // Compute grid velocity
      for(Index i(0); i < mesh.get_num_entities(0); ++i)
      {
        mesh_velocity[i] = DataType(0);
        for(Index d(0); d < MeshType::world_dim; ++d)
          mesh_velocity[i] += Math::sqr(ideltat*(coords_old[d](i) - rumpflpumpfl._coords[d](i)));

        mesh_velocity[i] = Math::sqrt(mesh_velocity[i]);
        if(mesh_velocity[i] > max_mesh_velocity)
          max_mesh_velocity = mesh_velocity[i];
      }
      std::cout << "max mesh velocity = " << scientify(max_mesh_velocity) << std::endl;

      if( n%outputstep || outputstep==1)
      {
        fval = rumpflpumpfl.compute_functional(func_norm,func_det,func_det2,func_lvlset);
        rumpflpumpfl.compute_gradient();
        std::cout << "fval post optimisation = " << scientify(fval) << std::endl;

        filename = "post_" + stringify(n) + ".vtk";
        Geometry::ExportVTK<MeshType> writer_post(mesh);

        writer_post.add_scalar_cell("lambda", rumpflpumpfl._lambda.elements() );
        writer_post.add_scalar_cell("h", rumpflpumpfl._h[0].elements() );
        writer_post.add_scalar_vertex("grad_0", &rumpflpumpfl._grad[0]);
        writer_post.add_scalar_vertex("grad_1", &rumpflpumpfl._grad[mesh.get_num_entities(0)]);
        writer_post.add_scalar_cell("norm", func_norm);
        writer_post.add_scalar_cell("det", func_det);
        writer_post.add_scalar_cell("det2", func_det2);
        writer_post.add_scalar_vertex("lvlset_grad_0", rumpflpumpfl._lvlset_grad_vtx_vec[0].elements());
        writer_post.add_scalar_vertex("lvlset_grad_1", rumpflpumpfl._lvlset_grad_vtx_vec[1].elements());
        writer_post.add_scalar_vertex("levelset", rumpflpumpfl._lvlset_vtx_vec.elements());
        writer_post.add_scalar_vertex("mesh_velocity", mesh_velocity);
        writer_post.add_scalar_cell("levelset_constraint", func_lvlset );
        //std::cout << "Writing " << filename << std::endl;
        writer_post.write(filename);

      }

      n++;
    }

    delete func_norm;
    delete func_det;
    delete func_det2;
    delete func_lvlset;

    delete mesh_velocity;
  }

int main()
{
  run<double,Shape::Simplex<2>>();
  return 0;
}
