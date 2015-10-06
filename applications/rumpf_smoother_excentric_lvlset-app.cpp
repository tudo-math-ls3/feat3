#include <kernel/base_header.hpp>
#ifdef FEAST_HAVE_ALGLIB
#include <kernel/archs.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/meshopt/rumpf_smoother_conc.hpp>
#include <kernel/meshopt/rumpf_smoother_lvlset.hpp>
#include <kernel/meshopt/rumpf_smoother_lvlset_q1hack.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_p1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_p1_d2.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1_d2.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1hack.hpp>
#include <kernel/meshopt/rumpf_functionals/lvlset.hpp>
#include <kernel/meshopt/rumpf_functionals/conc_2d_p1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/conc_2d_p1_d2.hpp>
#include <kernel/meshopt/rumpf_functionals/conc_2d_q1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/conc_2d_q1_d2.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/simple_arg_parser.hpp>

using namespace FEAST;

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
      template<int component, int dim_, typename DataType_>
      class DistanceFunctionSD_grad :
        public Analytic::Function
      {
        public:
        static constexpr int domain_dim = dim_;
        typedef Analytic::Image::Scalar ImageType;

        static constexpr bool can_value = true;

        typedef Tiny::Vector<DataType_, dim_> PointType;

          /** \copydoc AnalyticFunction::Evaluator */
          template<typename EvalTraits_>
          class Evaluator :
            public Analytic::Function::Evaluator<EvalTraits_>
        {
          public:
            /// coefficient data type
            typedef typename EvalTraits_::DataType DataType;
            /// evaluation point type
            typedef typename EvalTraits_::PointType PointType;
            /// value type
            typedef typename EvalTraits_::ValueType ValueType;
            /// gradient type
            typedef typename EvalTraits_::GradientType GradientType;
            /// hessian type
            typedef typename EvalTraits_::HessianType HessianType;

          private:
            /// Function to evaluate
            const DistanceFunctionSD_grad& _function;

          public:
            /// Constructor
            explicit Evaluator(const DistanceFunctionSD_grad& function) :
              _function(function)
              {
              }

            void value(ValueType& val, const PointType& point) const
            {
              PointType tmp (point - _function._point);
              DataType norm = tmp.norm_euclid();
              if(norm <= Math::eps<DataType>())
              {
                val = DataType(0);
              }
              else
              {
                val =  _function._b*(point[component] - _function._point(component))/norm;
              }
            }
        }; // class DistanceFunctionSD::Evaluator<...>

        private:
          /// The point to which the distance to is calculated
          PointType& _point;
          /// Displacement of the function
          DataType_ _a;
          /// Scaling factor
          DataType_ _b;

        public:
          /// Constructor
          explicit DistanceFunctionSD_grad(PointType& x0_, const DataType_ a_, const DataType_ b_) :
            _point(x0_),
            _a(a_),
            _b(b_)
            {
            }
          /// Sets _point to x0_
          void set_point(const PointType x0_)
          {
            _point = x0_;
          }

      }; // class DistanceFunctionSD_grad

      template<int component, int derivative_component, int dim_, typename DataType_>
      class PlaneDistanceFunctionSD_grad :
        public Analytic::Function
      {
        public:
        static constexpr int domain_dim = dim_;
        typedef Analytic::Image::Scalar ImageType;

        static constexpr bool can_value = true;

        typedef Tiny::Vector<DataType_, dim_> PointType;

          /** \copydoc AnalyticFunction::Evaluator */
          template<typename EvalTraits_>
          class Evaluator :
            public Analytic::Function::Evaluator<EvalTraits_>
        {
          public:
            /// coefficient data type
            typedef typename EvalTraits_::DataType DataType;
            /// evaluation point type
            typedef typename EvalTraits_::PointType PointType;
            /// value type
            typedef typename EvalTraits_::ValueType ValueType;
            /// gradient type
            typedef typename EvalTraits_::GradientType GradientType;
            /// hessian type
            typedef typename EvalTraits_::HessianType HessianType;

          private:
            /// Function to evaluate
            const PlaneDistanceFunctionSD_grad& _function;

          public:
            /// Constructor
            explicit Evaluator(const PlaneDistanceFunctionSD_grad& function) :
              _function(function)
              {
              }

            void value(ValueType& val, const PointType&) const
            {
              if(component == derivative_component)
                val = _function._b;
              else
                val = DataType(0);
            }

        }; // class PlaneDistanceFunctionSD::Evaluator<...>

        private:
          /// The point to which the distance to is calculated
          PointType& _point;
          /// Scaling factor
          DataType_ _b;

        public:
          /// Constructor
          explicit PlaneDistanceFunctionSD_grad(PointType& x0_, const DataType_ b_) :
            _point(x0_),
            _b(b_)
            {
            }

      }; // class DistanceFunctionSD_grad

      /**
       * \brief Function representing the gradient of the minimum of two analytic functions
       *
       * This is needed i.e. if there are several objects implicitly defined by their zero level sets. The class
       * in general supports function values, gradients and hessians for all dimensions, depending on the two analytic
       * functions supporting these.
       *
       * \Warning As min is non differentiable in general, Bad Things(TM) may happen when computing the gradient
       * and/or hessian where the function values are nearly identical.
       *
       * \tparam AnalyticFunctionType1
       * Type for the first AnalyticFunction
       *
       * \tparam AnalyticFunctionType2
       * Type for the second AnalyticFunction
       *
       * \author Jordi Paul
       */
      template<typename AnalyticFunctionType1, typename AnalyticFunctionType2, typename AnalyticFunctionType1Der, typename AnalyticFunctionType2Der>
      class MinOfTwoFunctionsDer :
        public Analytic::Function
      {
      public:
        /// our domain dimension
        static constexpr int domain_dim = AnalyticFunctionType1::domain_dim;
        /// our image type
        typedef Analytic::Image::Scalar ImageType;

        private:

          /// The first AnalyticFunction
          const AnalyticFunctionType1& _f1;
          /// The second AnalyticFunction
          const AnalyticFunctionType2& _f2;
          /// The first AnalyticFunctions derivative
          const AnalyticFunctionType1Der& _f1_der;
          /// The second AnalyticFunctions derivative
          const AnalyticFunctionType2Der& _f2_der;

        public:
          /// Can compute function values if both AnalyticFunctions can do that
          static constexpr bool can_value = (AnalyticFunctionType1Der::can_value && AnalyticFunctionType2Der::can_value);
          /// Can compute the function gradient if both AnalyticFunctions can do that
          static constexpr bool can_grad = (AnalyticFunctionType1Der::can_grad && AnalyticFunctionType2Der::can_grad);
          /// Can compute the function hessian if both AnalyticFunctions can do that
          static constexpr bool can_hess = (AnalyticFunctionType1Der::can_hess && AnalyticFunctionType2Der::can_hess);

          /** \copydoc AnalyticFunction::Evaluator */
          template<typename EvalTraits_>
          class Evaluator :
            public Analytic::Function::Evaluator<EvalTraits_>
        {
          public:
            /// coefficient data type
            typedef typename EvalTraits_::DataType DataType;
            /// evaluation point type
            typedef typename EvalTraits_::PointType PointType;
            /// value type
            typedef typename EvalTraits_::ValueType ValueType;
            /// gradient type
            typedef typename EvalTraits_::GradientType GradientType;
            /// hessian type
            typedef typename EvalTraits_::HessianType HessianType;

          private:
            /// Function to evaluate
            const MinOfTwoFunctionsDer& _function;
            typename AnalyticFunctionType1::template Evaluator<EvalTraits_> _f1_eval;
            typename AnalyticFunctionType2::template Evaluator<EvalTraits_> _f2_eval;
            typename AnalyticFunctionType1Der::template Evaluator<EvalTraits_> _f1_der_eval;
            typename AnalyticFunctionType2Der::template Evaluator<EvalTraits_> _f2_der_eval;

          public:
            /// Constructor
            explicit Evaluator(const MinOfTwoFunctionsDer& function) :
              _function(function),
              _f1_eval(function._f1),
              _f2_eval(function._f2),
              _f1_der_eval(function._f1_der),
              _f2_der_eval(function._f2_der)
              {
              }

            void value(ValueType& val, const PointType& point) const
            {
              ValueType fval1, fval2;
              _f1_eval.value(fval1, point);
              _f2_eval.value(fval2, point);
              if(Math::abs(fval1-fval2) < Math::eps<DataType>())
                val = DataType(0);
              else if(fval1 < fval2)
                _f1_der_eval.value(val, point);
              else
                _f2_der_eval.value(val, point);
            }

            void gradient(GradientType& grad, const PointType& point) const
            {
              ValueType fval1, fval2;
              _f1_eval.value(fval1, point);
              _f2_eval.value(fval2, point);
              if(Math::abs(fval1-fval2) < Math::eps<DataType>())
                grad.format();
              else if(fval1 < fval2)
                _f1_der_eval.gradient(grad, point);
              else
                _f2_der_eval.gradient(grad, point);
            }

            void hessian(HessianType& hess, const PointType& point) const
            {
              ValueType fval1, fval2;
              _f1_eval.value(fval1, point);
              _f2_eval.value(fval2, point);
              if(Math::abs(fval1 - fval2) < Math::eps<DataType>())
                hess.format();
              else if(fval1 < fval2)
                _f1_der_eval.hessian(hess, point);
              else
                _f2_der_eval.hessian(hess, point);
            }

        }; // class MinOfTwoFunctions::Evaluator<...>

        public:
          /// Constructor
          explicit MinOfTwoFunctionsDer(const AnalyticFunctionType1& f1_, const AnalyticFunctionType2& f2_,
          const AnalyticFunctionType1Der& f1_der_, const AnalyticFunctionType2Der& f2_der_) :
            _f1(f1_),
            _f2(f2_),
            _f1_der(f1_der_),
            _f2_der(f2_der_)
            {
            }

      }; // class MinOfTwoFunctionsDer


template<typename PointType, typename DataType>
void centre_point_outer(PointType& my_point, DataType time)
{
  my_point.v[0] = DataType(0.5) - DataType(0.125)*Math::cos(DataType(2)*Math::pi<DataType>()*DataType(time));
  my_point.v[1] = DataType(0.5) - DataType(0.125)*Math::sin(DataType(2)*Math::pi<DataType>()*DataType(time));
}

template<typename PointType, typename DataType>
void centre_point_inner(PointType& my_point, DataType time)
{
  my_point.v[0] = DataType(0.5) - DataType(0.1875)*Math::cos(DataType(2)*Math::pi<DataType>()*DataType(time));
  my_point.v[1] = DataType(0.5) - DataType(0.1875)*Math::sin(DataType(2)*Math::pi<DataType>()*DataType(time));
}

/**
 * \brief Wrapper struct as functions do not seem to agree with template template parameters
 **/
template
<
  typename DT_,
  typename MeshType_,
  template<typename ... > class RumpfSmootherType_,
  template<typename, typename> class FunctionalType_,
  template<typename, typename> class LevelsetFunctionalType_
  > struct LevelsetApp
{
  /// Precision for meshes etc, everything else uses the same data type
  typedef DT_ DataType;
  /// Rumpf Smoothers are implemented for Mem::Main only
  typedef Mem::Main MemType;
  /// So we use Index
  typedef Index IndexType;
  /// The type of the mesh
  typedef MeshType_ MeshType;
  /// Shape of the mesh cells
  typedef typename MeshType::ShapeType ShapeType;
  /// The corresponding transformation
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
  /// The Rumpf functional variant
  typedef FunctionalType_<DataType, ShapeType> FunctionalType;
  /// The levelset part of the Rumpf functional
  typedef LevelsetFunctionalType_<DataType, ShapeType> LevelsetFunctionalType;
  /// Type for points in the mesh
  typedef Tiny::Vector<DataType,2,2> ImgPointType;

  /**
   * \brief Routine that does the actual work
   *
   * \param[in] my_streamer
   * MeshStreamer that contains the data from the mesh file.
   *
   * \param[in] level
   * Number of refines.
   */
  static int run(MeshStreamer& my_streamer, Index level, DT_ deltat)
  {
    // Read mesh from the MeshStreamer and create the MeshAtlas
    std::cout << "Creating mesh atlas..." << std::endl;
    Geometry::MeshAtlas<MeshType_>* atlas = nullptr;
    try
    {
      atlas = new Geometry::MeshAtlas<MeshType_>(my_streamer);
    }
    catch(std::exception& exc)
    {
      std::cerr << "ERROR: " << exc.what() << std::endl;
      return 1;
    }

    // Create mesh node
    std::cout << "Creating mesh node..." << std::endl;
    Geometry::RootMeshNode<MeshType_>* rmn = nullptr;
    try
    {
      rmn = new Geometry::RootMeshNode<MeshType_>(my_streamer, atlas);
      rmn ->adapt();
    }
    catch(std::exception& exc)
    {
      std::cerr << "ERROR: " << exc.what() << std::endl;
      return 1;
    }

    // refine
    for(Index lvl(1); lvl <= level; ++lvl)
    {
      std::cout << "Refining up to level " << lvl << "..." << std::endl;
      auto* old = rmn;
      rmn = old->refine();
      delete old;
    }

    MeshType* mesh = rmn->get_mesh();

    typedef Analytic::Common::DistanceFunctionSD<2, DataType> AnalyticFunctionType1;
    typedef DistanceFunctionSD_grad<0, 2, DataType> AnalyticFunction1Grad0Type;
    typedef DistanceFunctionSD_grad<1, 2, DataType> AnalyticFunction1Grad1Type;

    typedef Analytic::Common::MinOfTwoFunctions<AnalyticFunctionType1, AnalyticFunctionType1> AnalyticFunctionType;
    typedef MinOfTwoFunctionsDer<AnalyticFunctionType1, AnalyticFunctionType1, AnalyticFunction1Grad0Type, AnalyticFunction1Grad0Type> AnalyticFunctionGrad0Type;
    typedef MinOfTwoFunctionsDer<AnalyticFunctionType1, AnalyticFunctionType1, AnalyticFunction1Grad1Type, AnalyticFunction1Grad1Type> AnalyticFunctionGrad1Type;

    // Reference point for the rotation
    ImgPointType x0(DataType(0));

    centre_point_outer(x0,DataType(0));
    DataType scaling_outer(-1);
    DataType radius_outer(0.35);

    // Analytic function for the distance to the outer circle's boundary
    AnalyticFunctionType1 outer(x0, radius_outer, scaling_outer);
    AnalyticFunction1Grad0Type outer_grad0(x0, radius_outer, scaling_outer);
    AnalyticFunction1Grad1Type outer_grad1(x0, radius_outer, scaling_outer);

    DataType scaling_inner(-scaling_outer);
    DataType radius_inner(-scaling_inner*0.275);
    centre_point_inner(x0,DataType(0));

    // Analytic function for the distance to the inner circle's boundary
    AnalyticFunctionType1 inner(x0, radius_inner, scaling_inner);
    AnalyticFunction1Grad0Type inner_grad0(x0, radius_inner, scaling_inner);
    AnalyticFunction1Grad1Type inner_grad1(x0, radius_inner, scaling_inner);

    AnalyticFunctionType analytic_lvlset(inner, outer);
    AnalyticFunctionGrad0Type analytic_lvlset_grad0(inner, outer, inner_grad0, outer_grad0);
    AnalyticFunctionGrad1Type analytic_lvlset_grad1(inner, outer, inner_grad1, outer_grad1);

    // Now we can define the whole smoother
    typedef RumpfSmootherType_
    <
      AnalyticFunctionType,
      AnalyticFunctionGrad0Type,
      AnalyticFunctionGrad1Type,
      TrafoType,
      FunctionalType,
      LevelsetFunctionalType
    > RumpfSmootherType;


    // Parameters for the functional
    DataType fac_norm = DataType(1e0),fac_det = DataType(1e0),fac_cof = DataType(0), fac_reg(DataType(1e-8));
    // Create the functional
    FunctionalType my_functional(fac_norm, fac_det, fac_cof, fac_reg);
    // Parameters for the Rumpf smoother
    bool align_to_lvlset(false);
    bool r_adaptivity(true);
    DataType r_adapt_reg = DataType(1e-2), r_adapt_pow = DataType(0.5);
    // Create levelset part of the functional
    LevelsetFunctionalType my_levelset_functional;

    // Set slip boundary conditions at these parts of the boundary
    std::deque<String> slip_list;
    slip_list.push_back("left");
    slip_list.push_back("right");

    // Set Dirichlet boundary conditions at these parts of the boundary
    std::deque<String> dirichlet_list;
    dirichlet_list.push_back("bottom");
    dirichlet_list.push_back("top");

    // The smoother in all its template glory
    RumpfSmootherType rumpflpumpfl(rmn, dirichlet_list, slip_list, my_functional, my_levelset_functional,
    align_to_lvlset, r_adaptivity, analytic_lvlset, analytic_lvlset_grad0, analytic_lvlset_grad1);

    // Set the r-adaptivity parameters
    rumpflpumpfl.set_r_adapt_params(r_adapt_reg,r_adapt_pow);

    rumpflpumpfl.init();
    rumpflpumpfl.print();

    // Arrays for saving the contributions of the different Rumpf functional parts
    DataType* func_norm(new DataType[mesh->get_num_entities(MeshType::shape_dim)]);
    DataType* func_det(new DataType[mesh->get_num_entities(MeshType::shape_dim)]);
    DataType* func_rec_det(new DataType[mesh->get_num_entities(MeshType::shape_dim)]);
    DataType* func_lvlset(new DataType[mesh->get_num_entities(MeshType::shape_dim)]);

    // Evaluates the levelset function and its gradient
    rumpflpumpfl.prepare();
    // Compute initial functional value
    DataType fval(0);
    fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_rec_det, func_lvlset);
    std::cout << "fval pre optimisation = " << scientify(fval) << " cell size quality indicator: " << scientify(rumpflpumpfl.cell_size_quality()) << std::endl;
    // Compute initial functional gradient
    rumpflpumpfl.compute_gradient();

    // Filename for vtk files
    std::string filename;

    // Write initial state to file
    Geometry::ExportVTK<MeshType> writer_pre_initial(*mesh);
    writer_pre_initial.add_scalar_cell("lambda", rumpflpumpfl._lambda.elements() );
    writer_pre_initial.add_field_cell_blocked_vector("h", rumpflpumpfl._h );
    writer_pre_initial.add_scalar_vertex("levelset", rumpflpumpfl._lvlset_vtx_vec.elements());
    writer_pre_initial.add_field_vertex("lvlset_grad", rumpflpumpfl._lvlset_grad_vtx_vec[0].elements(), rumpflpumpfl._lvlset_grad_vtx_vec[1].elements());
    writer_pre_initial.add_field_cell("fval", func_norm, func_det, func_rec_det);
    writer_pre_initial.add_scalar_cell("levelset_constraint", func_lvlset );
    writer_pre_initial.add_field_vertex_blocked_vector("grad", rumpflpumpfl._grad);
    writer_pre_initial.write("pre_initial");

    // Optimise the mesh
    rumpflpumpfl.optimise();

    fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_rec_det, func_lvlset);
    std::cout << "fval post optimisation = " << scientify(fval) << " cell size quality indicator: " << scientify(rumpflpumpfl.cell_size_quality()) << std::endl;

    // Call prepare() again because the mesh changed due to the optimisation and it was not called again after the
    // last iteration
    rumpflpumpfl.prepare();
    rumpflpumpfl.compute_gradient();

    // Write optimised initial mesh
    Geometry::ExportVTK<MeshType> writer_post_initial(*mesh);
    writer_post_initial.add_scalar_cell("lambda", rumpflpumpfl._lambda.elements() );
    writer_post_initial.add_field_cell_blocked_vector("h", rumpflpumpfl._h );
    writer_post_initial.add_scalar_vertex("levelset", rumpflpumpfl._lvlset_vtx_vec.elements());
    writer_post_initial.add_field_vertex("lvlset_grad", rumpflpumpfl._lvlset_grad_vtx_vec[0].elements(), rumpflpumpfl._lvlset_grad_vtx_vec[1].elements());
    writer_post_initial.add_field_cell("fval", func_norm, func_det, func_rec_det);
    writer_post_initial.add_scalar_cell("levelset_constraint", func_lvlset );
    writer_post_initial.add_field_vertex_blocked_vector("grad", rumpflpumpfl._grad);
    writer_post_initial.write("post_initial");

    DataType time(0);
    Index n(0);

    // Old mesh coordinates for computing the mesh velocity
    LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, MeshType::world_dim>
      coords_old(mesh->get_num_entities(0),DataType(0));
    LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, MeshType::world_dim>
      mesh_velocity(mesh->get_num_entities(0), DataType(0));

    while(time < DataType(1))
    {
      std::cout << "timestep " << n << std::endl;
      time+= deltat;

      // Save old vertex coordinates
      coords_old.clone(rumpflpumpfl._coords);

      // Update levelset function
      centre_point_outer(x0,time);
      outer.set_point(x0);
      outer_grad0.set_point(x0);
      outer_grad1.set_point(x0);

      // Update reference point
      centre_point_inner(x0,time);
      inner.set_point(x0);
      inner_grad0.set_point(x0);
      inner_grad1.set_point(x0);

      // Compute functional value and gradient
      rumpflpumpfl.prepare();
      fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_rec_det, func_lvlset);
      rumpflpumpfl.compute_gradient();
      std::cout << "fval pre optimisation = " << scientify(fval) << " cell size quality indicator: " <<
        scientify(rumpflpumpfl.cell_size_quality()) << std::endl;

      // Write pre-optimisation mesh
      filename = "pre_" + stringify(n);
      Geometry::ExportVTK<MeshType> writer_pre(*mesh);
      writer_pre.add_scalar_cell("lambda", rumpflpumpfl._lambda.elements() );
      writer_pre.add_field_cell_blocked_vector("h", rumpflpumpfl._h);
      writer_pre.add_scalar_vertex("levelset", rumpflpumpfl._lvlset_vtx_vec.elements());
      writer_pre.add_field_vertex("lvlset_grad", rumpflpumpfl._lvlset_grad_vtx_vec[0].elements(),
      rumpflpumpfl._lvlset_grad_vtx_vec[1].elements());
      writer_pre.add_field_cell("fval", func_norm, func_det, func_rec_det);
      writer_pre.add_scalar_cell("levelset_constraint", func_lvlset );
      writer_pre.add_field_vertex_blocked_vector("grad", rumpflpumpfl._grad);
      writer_pre.write(filename);

      // Optimise the mesh
      rumpflpumpfl.optimise();

      // Compute grid velocity
      DataType max_mesh_velocity(-1e10);
      DataType ideltat = DataType(1)/deltat;
      for(Index i(0); i < mesh->get_num_entities(0); ++i)
      {
        mesh_velocity(i, ideltat*(rumpflpumpfl._coords(i) - coords_old(i)));

        DataType my_mesh_velocity(mesh_velocity(i).norm_euclid());

        if(my_mesh_velocity > max_mesh_velocity)
          max_mesh_velocity = my_mesh_velocity;
      }
      std::cout << "max mesh velocity = " << scientify(max_mesh_velocity) << std::endl;

      rumpflpumpfl.prepare();
      fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_rec_det, func_lvlset);
      rumpflpumpfl.compute_gradient();
      std::cout << "fval post optimisation = " << scientify(fval) << " cell size quality indicator: "
      << scientify(rumpflpumpfl.cell_size_quality()) << std::endl;

      // Write post-optimisation mesh
      filename = "post_" + stringify(n);
      Geometry::ExportVTK<MeshType> writer_post(*mesh);
      writer_post.add_scalar_cell("lambda", rumpflpumpfl._lambda.elements() );
      writer_post.add_field_cell_blocked_vector("h", rumpflpumpfl._h);
      writer_post.add_scalar_vertex("levelset", rumpflpumpfl._lvlset_vtx_vec.elements());
      writer_post.add_field_vertex("lvlset_grad", rumpflpumpfl._lvlset_grad_vtx_vec[0].elements(),
      rumpflpumpfl._lvlset_grad_vtx_vec[1].elements());
      writer_post.add_field_cell("fval", func_norm, func_det, func_rec_det);
      writer_post.add_scalar_cell("levelset_constraint", func_lvlset );
      writer_post.add_field_vertex_blocked_vector("grad", rumpflpumpfl._grad);
      writer_post.add_field_vertex_blocked_vector("mesh_velocity", mesh_velocity);
      writer_post.write(filename);

      n++;

    }

    // Clean up
    delete rmn;
    if(atlas != nullptr)
      delete atlas;

    delete[] func_norm;
    delete[] func_det;
    delete[] func_rec_det;
    delete[] func_lvlset;

    return 0;

  }

}; // struct LevelsetApp

template<typename A, typename B, typename C, typename D, typename E, typename F>
using MySmoother = Meshopt::RumpfSmootherLevelsetConcAnalytic<A, B, C, D, E, F>;

template<typename A, typename B, typename C, typename D, typename E, typename F>
using MySmootherQ1Hack = Meshopt::RumpfSmootherLevelsetAnalyticQ1Hack<A, B, C, D, E, F>;

template<typename A, typename B>
using MyFunctional= Meshopt::RumpfFunctionalConc<A, B>;

template<typename A, typename B>
using MyFunctionalQ1Hack = Meshopt::RumpfFunctionalQ1Hack<A, B, Meshopt::RumpfFunctional>;

int main(int argc, char* argv[])
{
  // Creata a parser for command line arguments.
  SimpleArgParser args(argc, argv);

  if( args.check("help") > -1 || args.num_args()==1)
  {
    std::cout << "Rumpf Smoother Application for Levelset usage: " << std::endl;
    std::cout << "Required arguments: --filename [String]: Path to a FEAST mesh file." << std::endl;
    std::cout << "Optional arguments: --level [unsigned int]: Number of refines, defaults to 0." << std::endl;
    exit(1);
  }
  // Specify supported command line switches
  args.support("level");
  args.support("filename");
  args.support("help");
  // Refinement level
  Index lvl_max(0);
  // Input file name, required
  FEAST::String filename;
  // Get unsupported command line arguments
  std::deque<std::pair<int,String> > unsupported = args.query_unsupported();
  if( !unsupported.empty() )
  {
    // print all unsupported options to cerr
    for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
      std::cerr << "ERROR: unsupported option '--" << (*it).second << "'" << std::endl;
  }

  // Check and parse --filename
  if(args.check("filename") != 1 )
    throw InternalError(__func__, __FILE__, __LINE__, "Invalid option for --filename");
  else
  {
    args.parse("filename", filename);
    std::cout << "Reading mesh from file " << filename << std::endl;
  }

  // Check and parse --level
  if(args.check("level") != 1)
    std::cout << "No refinement level specified, defaulting to 0." << std::endl;
  else
  {
    args.parse("level", lvl_max);
    std::cout << "Refinement level " << lvl_max << std::endl;
  }

  // Create a MeshStreamer and read the mesh file
  MeshStreamer my_streamer;
  my_streamer.parse_mesh_file(filename);

  // This is the raw mesh data my_streamer read from filename
  auto& mesh_data = my_streamer.get_root_mesh_node()->mesh_data;
  // Marker int for the MeshType
  int mesh_type = mesh_data.mesh_type;
  // Marker int for the ShapeType
  int shape_type = mesh_data.shape_type;

  ASSERT(mesh_type == mesh_data.mt_conformal, "This application only works for conformal meshes!");

  typedef double DataType;

  DataType deltat(DataType(1e-2));

  // This is the list of all supported meshes that could appear in the mesh file
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 2, 2, Real> Simplex2Mesh_2d;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, 2, Real> Hypercube2Mesh_2d;

  // Call the run() method of the appropriate wrapper class
  if(shape_type == mesh_data.st_tria)
    return LevelsetApp<DataType, Simplex2Mesh_2d, MySmoother, MyFunctional, Meshopt::RumpfFunctionalLevelset>::
      run(my_streamer, lvl_max, deltat);
  if(shape_type == mesh_data.st_quad)
    return LevelsetApp<DataType, Hypercube2Mesh_2d, MySmoother, MyFunctional, Meshopt::RumpfFunctionalLevelset>::
      run(my_streamer, lvl_max, deltat);

  // If no MeshType from the list was in the file, return 1
  return 1;
}
/// \endcond
#endif // FEAST_HAVE_ALGLIB
