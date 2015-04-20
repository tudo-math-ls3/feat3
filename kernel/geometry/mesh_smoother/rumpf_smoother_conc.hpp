#pragma once
#ifndef KERNEL_GEOMETRY_RUMPF_SMOOTHER_CONC_HPP
#define KERNEL_GEOMETRY_RUMPF_SMOOTHER_CONC_HPP 1

#include <kernel/geometry/mesh_smoother/rumpf_smoother_lvlset.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
// DEBUG
#include <kernel/geometry/export_vtk.hpp>


namespace FEAST
{
  namespace Geometry
  {
    /**
     * \brief Variant with the levelset function given explicitly
     *
     * \tparam AnalyticFunctionType_
     * Analytic function providing level set information
     *
     * \tparam AnalyticFunctionGrad0Type_
     * 1st component function of the AnalyticFunction's gradient
     *
     * \tparam AnalyticFunctionGrad0Type_
     * 2nd component function of the AnalyticFunction's gradient
     *
     * \tparam SpaceType_
     * FE space for the levelset function
     *
     * \tparam FunctionalType_
     * Functional used for determining mesh quality.
     *
     * \tparam TrafoType_
     * Our transformation.
     *
     * \tparam DataType_
     * Our datatype.
     *
     * \tparam MemType_
     * Memory architecture.
     *
     * \tparam H_EvalType
     * Local meshsize evaluator
     *
     * \author Jordi Paul
     *
     **/
    template
    <
      typename AnalyticFunctionType_,
      typename AnalyticFunctionGrad0Type_,
      typename AnalyticFunctionGrad1Type_,
      typename DataType_,
      typename MemType_,
      typename TrafoType_,
      typename FunctionalType_,
      typename LevelsetFunctionalType_,
      typename SpaceType_,
      typename H_EvalType_ = H_Evaluator<TrafoType_, DataType_>
    >
    class RumpfSmootherLevelsetConcAnalytic :
      public RumpfSmootherLevelsetAnalytic
      <
        AnalyticFunctionType_,
        AnalyticFunctionGrad0Type_,
        AnalyticFunctionGrad1Type_,
        DataType_,
        MemType_,
        TrafoType_,
        FunctionalType_,
        LevelsetFunctionalType_,
        SpaceType_,
        H_EvalType_
      >
    {
      public:
        typedef DataType_ DataType;
        /// Transformation type
        typedef TrafoType_ TrafoType;
        /// Baseclass
        typedef RumpfSmootherLevelsetAnalytic
        <
          AnalyticFunctionType_,
          AnalyticFunctionGrad0Type_,
          AnalyticFunctionGrad1Type_,
          DataType_,
          MemType_,
          TrafoType_,
          FunctionalType_,
          LevelsetFunctionalType_,
          SpaceType_,
          H_EvalType_
        > BaseClass;

        /// Vector types for element sizes etc.
        typedef LAFEM::DenseVector<MemType_, DataType_> VectorType;
        /// Our type of mesh
        typedef typename BaseClass::MeshType MeshType;
        /// Our shape type
        typedef typename BaseClass::ShapeType ShapeType;

        /// Vector for saving the mesh concentration, one entry per cell
        VectorType _conc;
        /// The sum of all entries in _conc
        DataType _sum_conc;
        /// Gradient of the mesh concentration wrt. the world coordinates
        VectorType _grad_conc[MeshType::world_dim];
        /// Gradient of the local optimal scales h wrt. the vertex coordinates. Each entry in the DenseVectorBlocked
        /// represents one cell. Each cell's block contains world_dim*(number of local vertices) entries.
        LAFEM::DenseVectorBlocked<MemType_, DataType_, Index, MeshType::world_dim*Shape::FaceTraits<ShapeType,0>::count> _grad_h;

        /**
         * \copydoc RumpfSmootherLevelset()
         *
         * \param[in] align_to_lvlset_
         * Align vertices/edges/faces with 0 levelset?
         *
         * \param[in] r_adaptivity_
         * Use levelset-based r_adaptivity?
         *
         * \param[in] analytic_function_
         * The analytic function representing the levelset function
         *
         * \param[in] analytic_function_grad0_
         * The first component of the gradient
         *
         * \param[in] analytic_function_grad1_
         * The second component of the gradient
         *
         **/
        explicit RumpfSmootherLevelsetConcAnalytic(
          TrafoType& trafo_,
          FunctionalType_& functional_,
          LevelsetFunctionalType_& levelset_functional_,
          bool align_to_lvlset_,
          bool r_adaptivity_,
          AnalyticFunctionType_& analytic_function_,
          AnalyticFunctionGrad0Type_& analytic_function_grad0_,
          AnalyticFunctionGrad1Type_& analytic_function_grad1_)
          : BaseClass(trafo_, functional_, levelset_functional_, align_to_lvlset_, r_adaptivity_,
          analytic_function_, analytic_function_grad0_, analytic_function_grad1_),
          _conc(trafo_.get_mesh().get_num_entities(ShapeType::dimension)),
          _sum_conc(DataType(0)),
          _grad_h(this->_mesh.get_num_entities(ShapeType::dimension),DataType(0))
          {
            for(int d(0); d < MeshType::world_dim; ++d)
              _grad_conc[d]= std::move(VectorType(this->_mesh.get_num_entities(ShapeType::dimension)));
            this->_update_h = true;
          }

        /**
         * \copydoc RumpfSmootherLevelset::prepare()
         *
         * In this case, it evaluates the analytic levelset function on the current mesh and then evaluates it in the
         * current mesh's vertices.
         *
         **/
        virtual void prepare() override
        {
          this->set_coords();
          // Evaluate levelset function
          Assembly::Interpolator::project(this->_lvlset_vec, this->_analytic_lvlset, this->_lvlset_space);
          // Evaluate the gradient of the levelset function
          Assembly::Interpolator::project(this->_lvlset_grad_vec[0], this->_analytic_lvlset_grad0, this->_lvlset_space);
          Assembly::Interpolator::project(this->_lvlset_grad_vec[1], this->_analytic_lvlset_grad1, this->_lvlset_space);
          // Project the levelset function to a grid vector on the new grid and...
          FEAST::Assembly::DiscreteVertexProjector::project(this->_lvlset_vtx_vec, this->_lvlset_vec, this->_lvlset_space);

          // ... do the same for its gradient
          FEAST::Assembly::DiscreteVertexProjector::project(this->_lvlset_grad_vtx_vec[0],this->_lvlset_grad_vec[0], this->_lvlset_space);
          FEAST::Assembly::DiscreteVertexProjector::project(this->_lvlset_grad_vtx_vec[1],this->_lvlset_grad_vec[1], this->_lvlset_space);

          // As the levelset might be used for r-adaptivity and it's vertex vector might have changed, re-compute
          if(this->_r_adaptivity)
          {
            this->compute_lambda();

            // DEBUG
            //for(Index cell(0); cell < this->_mesh.get_num_entities(ShapeType::dimension); cell++ )
            //  this->_mu(cell,this->_lambda(cell));

            this->compute_h();
            compute_grad_h();
          }
        }

        /**
         * \brief Computes the mesh concentration function for each cell
         **/
        void compute_conc()
        {
          // Total number of cells in the mesh
          Index ncells(this->_mesh.get_num_entities(ShapeType::dimension));
          // Index set for local/global numbering
          auto& idx = this->_mesh.template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector <DataType_, MeshType::world_dim> h;

          _sum_conc = DataType(0);
          for(Index cell(0); cell < ncells; ++cell)
          {
            // This will be the average of the levelset values at the vertices
            DataType tmp(0);

            for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
              tmp += this->_lvlset_vtx_vec(idx(cell,j));

            tmp = tmp/DataType(Shape::FaceTraits<ShapeType,0>::count);

            this->_conc(cell, conc_val(tmp));
            _sum_conc+=this->_conc(cell);
          }
        }

        /**
         * \brief Computes the value of the concentration function
         *
         * \param[in] dist
         * The distance.
         *
         **/
        DataType conc_val(DataType dist)
        {
          //return DataType(1);
          if(dist < DataType(0))
            return DataType(2)*Math::pow(this->_r_adapt_reg + Math::abs(dist),this->_r_adapt_pow);

          return Math::pow(this->_r_adapt_reg + Math::abs(dist),this->_r_adapt_pow);
        }

        /**
         * \brief Computes the derivative of the concentration function
         *
         * \param[in] dist
         * The distance.
         *
         **/
        DataType conc_der(DataType dist)
        {
          //return DataType(0);
          if(dist < DataType(0))
            return DataType(2)*this->_r_adapt_pow * Math::pow(this->_r_adapt_reg + Math::abs(dist),this->_r_adapt_pow - DataType(1))*Math::signum(dist);

          return this->_r_adapt_pow * Math::pow(this->_r_adapt_reg + Math::abs(dist),this->_r_adapt_pow - DataType(1))*Math::signum(dist);
        }

        /**
         * \brief Computes the local gradient of the concentration function wrt. the vertices
         *
         * \tparam Tgrad_
         * Type of the local gradient, i.e. Tiny::Matrix
         *
         * \tparam Tl_
         * Type for the vector of levelset values, i.e. Tiny::vector
         *
         * \tparam Tgradl_
         * Type for the local gradient of the levelset values wrt. the vertices, i.e. Tiny::Matrix
         *
         * \param[out] grad_loc
         * The gradient of the concentration wrt. the vertices
         *
         * \param[in] lvlset_vals_
         * The levelset values at the vertices
         *
         * \param[in] lvlset_grad_vals_
         * The grandient of the levelset wrt. the vertices, evaluated at the vertices themselves
         *
         **/
        template<typename Tgrad_, typename Tl_, typename Tgradl_>
        void compute_grad_conc_local(Tgrad_& grad_loc_, const Tl_& lvlset_vals_, const Tgradl_& lvlset_grad_vals_)
        {
          grad_loc_ = DataType(0);

            // This will be the average of the levelset values at the vertices
          DataType val(0);
          for(Index i(0); i < Index(Shape::FaceTraits<ShapeType,0>::count); ++i)
            val += lvlset_vals_(i);

          val = val/DataType(Shape::FaceTraits<ShapeType,0>::count);

          for(Index d(0); d < Index(MeshType::world_dim); ++d)
          {
            for(Index i(0); i < Index(Shape::FaceTraits<ShapeType,0>::count); ++i)
              grad_loc_(d,i) = conc_der(val) * lvlset_grad_vals_(d,i)
                /DataType(Shape::FaceTraits<ShapeType,0>::count);
          }
        }

        /**
         * \copydoc Baseclass::compute_lambda()
         *
         *  In this case,
         *  \f$ \lambda(T) = \frac{\mathrm{conc}(T)}{\sum_{K \in \mathcal{T} \mathrm{conc}(t)}} \f$
         *
         **/
        virtual void compute_lambda() override
        {
          compute_conc();
          for(Index cell(0); cell < this->_mesh.get_num_entities(ShapeType::dimension); ++cell)
          {
            this->_lambda(cell, this->_conc(cell)/this->_sum_conc);
          }
        }

        /**
         * \brief Computes the local gradient of the optimal scales
         *
         **/
        void compute_grad_h()
        {
          VectorType grad_sum_det[MeshType::world_dim];
          for(int d(0); d < MeshType::world_dim; ++d)
            grad_sum_det[d] = std::move(VectorType(this->_mesh.get_num_entities(0), DataType(0)));

          // Index set for local/global numbering
          auto& idx = this->_mesh.template get_index_set<ShapeType::dimension,0>();
          // This will hold the levelset values at the mesh vertices for one element
          FEAST::Tiny::Vector <DataType_, Shape::FaceTraits<ShapeType,0>::count> lvlset_vals;
          // This will hold the levelset gradient values for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> lvlset_grad_vals;
          // This will hold the local gradient values for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> grad_loc(0);
          // This will hold the computed local gradient values for one element for copy assigning to the blocked
          // datatype
          FEAST::Tiny::Vector<DataType, MeshType::world_dim*Shape::FaceTraits<ShapeType,0>::count> tmp(0);

          DataType exponent = DataType(1)/DataType(MeshType::world_dim) - DataType(1);

          DataType sum_det = H_EvalType_::compute_sum_det(this->_coords, this->_trafo);
          H_EvalType_::compute_grad_sum_det(grad_sum_det, this->_coords, this->_trafo);

          for(Index cell(0); cell < this->_mesh.get_num_entities(ShapeType::dimension); ++cell)
          {
            grad_loc.format();
            for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              // Get levelset
              lvlset_vals(j) = this->_lvlset_vec(idx(cell,j));
              for(Index d(0); d < Index(MeshType::world_dim); ++d)
                // Get levelset gradient
                lvlset_grad_vals(d,j) = this->_lvlset_grad_vtx_vec[d](idx(cell,j));
            }

            compute_grad_conc_local(grad_loc, lvlset_vals, lvlset_grad_vals);

            for(Index d(0); d < Index(MeshType::world_dim); ++d)
            {
              for(Index i(0); i < Shape::FaceTraits<ShapeType,0>::count; ++i)
                tmp(d*Shape::FaceTraits<ShapeType,0>::count + i) =
                  DataType(1)/DataType(MeshType::world_dim)*Math::pow(_conc(cell)/_sum_conc*sum_det,exponent)
                  *( _conc(cell)/_sum_conc*grad_sum_det[d](idx(cell,i)) +
                      ( grad_loc(d,i)*this->_sum_conc + this->_conc(cell)*this->_grad_conc[d](idx(cell,i)))
                      / Math::sqr(this->_sum_conc)*sum_det);
            }
            _grad_h(cell, tmp);

          }

        } // compute_grad_h

        /**
         * \brief Computes the gradient of the mesh concentration
         *
         **/
        void compute_grad_conc()
        {
          // Index set for local/global numbering
          auto& idx = this->_mesh.template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector <DataType_,MeshType::world_dim> h;
          // This will hold the local gradient for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> grad_loc;
          // This will hold the levelset values at the mesh vertices for one element
          FEAST::Tiny::Vector <DataType_, Shape::FaceTraits<ShapeType,0>::count> lvlset_vals;
          // This will hold the levelset gradient values for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> lvlset_grad_vals;

          // Compute the functional value for each cell
          for(Index cell(0); cell < this->_mesh.get_num_entities(ShapeType::dimension); ++cell)
          {
            for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              // Get levelset
              lvlset_vals(j) = this->_lvlset_vec(idx(cell,j));
              for(Index d(0); d < Index(MeshType::world_dim); ++d)
                // Get levelset gradient
                lvlset_grad_vals(d,j) = this->_lvlset_grad_vtx_vec[d](idx(cell,j));
            }

            compute_grad_conc_local(grad_loc, lvlset_vals, lvlset_grad_vals);

            for(Index d(0); d < Index(MeshType::world_dim); ++d)
            {
              // Add local contributions to global gradient vector
              for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
                this->_grad_conc[d](idx(cell,j)) += grad_loc(d,j);
            }
          }

          this->_filter_grad();

        } // compute_grad_conc()

        /// \copydoc RumpfSmoother::compute_gradient()
        virtual void compute_gradient() override
        {
          // Total number of vertices in the mesh
          Index nvertices(this->_mesh.get_num_entities(0));
          // Total number of cells in the mesh
          Index ncells(this->_mesh.get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->_mesh.template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector <DataType_,MeshType::world_dim> h;
          // This will hold the local gradient for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> grad_loc;
          // This will hold the levelset values at the mesh vertices for one element
          FEAST::Tiny::Vector <DataType_, Shape::FaceTraits<ShapeType,0>::count> lvlset_vals;
          // This will hold the levelset gradient values for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> lvlset_grad_vals;

          // Clear gradient vector
          for(Index i(0); i < Index(MeshType::world_dim)*nvertices; ++i)
            this->_grad[i] = DataType_(0);

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              // Get levelset
              lvlset_vals(j) = this->_lvlset_vec(idx(cell,j));
              for(Index d(0); d < Index(MeshType::world_dim); ++d)
              {
                h(d) = this->_h[d](cell);
                // Get local coordinates
                x(d,j) = this->_coords[d](idx(cell,j));
                // Get levelset gradient
                lvlset_grad_vals(d,j) = this->_lvlset_grad_vtx_vec[d](idx(cell,j));
              }
            }

            // Compute local gradient and save it to grad_loc
            this->_functional.compute_local_grad(x, h, grad_loc);
            //grad_loc.format();
            this->_functional.add_grad_h_part(grad_loc, x, h, this->_grad_h(cell));
            // This does not contain the weighting yet
            grad_loc *= this->_mu(cell);
            // Add levelset penalty term, which is not weighted with lambda
            this->_lvlset_functional.add_lvlset_penalty_grad(lvlset_vals, lvlset_grad_vals, grad_loc, this->lvlset_constraint_last);

            for(Index d(0); d < Index(MeshType::world_dim); ++d)
            {
              // Add local contributions to global gradient vector
              for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
                this->_grad[d*nvertices + idx(cell,j)] += grad_loc(d,j);
            }
          }

          this->_filter_grad();

        } // compute_gradient
    }; // class RumpfSmootherLevelsetConcAnalytic

  } // namespace Geometry
} // namespace FEAST
#endif // KERNEL_GEOMETRY_RUMPF_SMOOTHER_CONC_HPP
