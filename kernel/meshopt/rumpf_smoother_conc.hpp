#pragma once
#ifndef KERNEL_MESHOPT_RUMPF_SMOOTHER_CONC_HPP
#define KERNEL_MESHOPT_RUMPF_SMOOTHER_CONC_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/meshopt/rumpf_smoother_lvlset.hpp>

namespace FEAST
{
  namespace Meshopt
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
     * \tparam CoordType
     * Our datatype.
     *
     * \tparam MemType
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
      typename TrafoType_,
      typename FunctionalType_,
      typename LevelsetFunctionalType_,
      typename H_EvalType_ = H_Evaluator<TrafoType_, typename TrafoType_::MeshType::CoordType>
    >
    class RumpfSmootherLevelsetConcAnalytic :
      public RumpfSmootherLevelsetAnalytic
      <
        AnalyticFunctionType_,
        AnalyticFunctionGrad0Type_,
        AnalyticFunctionGrad1Type_,
        TrafoType_,
        FunctionalType_,
        LevelsetFunctionalType_,
        H_EvalType_
      >
    {
      public:
        /// Type for the transformation
        typedef TrafoType_ TrafoType;
        /// The mesh the transformation is defined on
        typedef typename TrafoType::MeshType MeshType;
        /// The precision of the mesh coordinates
        typedef typename MeshType::CoordType CoordType;

        /// Type for the functional
        typedef FunctionalType_ FunctionalType;
        /// Type for the levelset part of the functional
        typedef LevelsetFunctionalType_ LevelsetFunctionalType;

        /// Only Mem::Main is supported atm
        typedef Mem::Main MemType;
        /// We always use Index for now
        typedef Index IndexType;

        /// ShapeType of said mesh
        typedef typename MeshType::ShapeType ShapeType;
        /// Vector type for element sizes etc.
        typedef LAFEM::DenseVector<MemType, CoordType, IndexType> ScalarVectorType;
        /// Vector type for element scales etc.
        typedef LAFEM::DenseVectorBlocked<MemType, CoordType, IndexType, MeshType::world_dim> VectorType;
        /// Type for the filter enforcing boundary conditions
        typedef LAFEM::UnitFilterBlocked<MemType, CoordType, IndexType, MeshType::world_dim> FilterType;
        /// Finite Element space for the transformation
        typedef typename Intern::TrafoFE<TrafoType>::Space TrafoSpace;

        /// Baseclass
        typedef RumpfSmootherLevelsetAnalytic
        <
          AnalyticFunctionType_,
          AnalyticFunctionGrad0Type_,
          AnalyticFunctionGrad1Type_,
          TrafoType_,
          FunctionalType_,
          LevelsetFunctionalType_,
          H_EvalType_
        > BaseClass;

        /// Vector for saving the mesh concentration, one entry per cell
        ScalarVectorType _conc;
        /// The sum of all entries in _conc
        CoordType _sum_conc;
        /// Gradient of the mesh concentration wrt. the world coordinates
        VectorType _grad_conc;
        /// Gradient of the local optimal scales h wrt. the vertex coordinates.
        // Each entry in the DenseVectorBlocked represents one cell. Each cell's block contains
        // world_dim*(number of local vertices) entries. Has to be serialised like this because there is no
        // DenseVector that saves a Tiny::Matrix
        LAFEM::DenseVectorBlocked<MemType, CoordType, Index, MeshType::world_dim*Shape::FaceTraits<ShapeType,0>::count> _grad_h;

        /**
         * \copydoc RumpfSmootherLevelsetAnalytic()
         *
         * \note Because the dependence of \f$ h \f$ on the vertex coordinates is taken into account here, \c
         * _update_h is set to \c true.
         *
         **/
        explicit RumpfSmootherLevelsetConcAnalytic(Geometry::RootMeshNode<MeshType>* rmn_,
        std::deque<String>& dirichlet_list_, std::deque<String> slip_list_,
          FunctionalType_& functional_, LevelsetFunctionalType_& lvlset_functional_,
          bool align_to_lvlset_, bool r_adaptivity_,
          AnalyticFunctionType_& analytic_function_,
          AnalyticFunctionGrad0Type_& analytic_function_grad0_, AnalyticFunctionGrad1Type_& analytic_function_grad1_)
          : BaseClass(rmn_, dirichlet_list_, slip_list_, functional_, lvlset_functional_,
          align_to_lvlset_, r_adaptivity_,
          analytic_function_, analytic_function_grad0_, analytic_function_grad1_),
          _conc(rmn_->get_mesh()->get_num_entities(ShapeType::dimension)),
          _sum_conc(CoordType(0)),
          _grad_conc(rmn_->get_mesh()->get_num_entities(0), CoordType(0)),
          _grad_h(rmn_->get_mesh()->get_num_entities(ShapeType::dimension),CoordType(0))
          {
            this->_update_h = true;
          }

        /// \copydoc BaseClass::init()
        virtual void init() override
        {
          // Write any potential changes to the mesh
          this->set_coords();

          // Assemble the homogeneous filter
          this->_slip_asm.assemble(this->_filter.template at<0>(), this->_trafo_space);

          // Call prepare to evaluate the levelset function
          prepare();

          // Compute element weights
          this->compute_mu();
          // Compute desired element size distribution
          this->compute_lambda();
          // Compute target scales
          this->compute_h();
          // Compute their gradient wrt. the vertex coordinates
          compute_grad_h();
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
          //// Evaluate the gradient of the levelset function
          //Assembly::Interpolator::project(this->_lvlset_grad_vec[0], this->_analytic_lvlset_grad0, this->_lvlset_space);
          //Assembly::Interpolator::project(this->_lvlset_grad_vec[1], this->_analytic_lvlset_grad1, this->_lvlset_space);
          //BaseClass::prepare();

          //if(this->_r_adaptivity)
          //{
          //  compute_grad_h();
          //}

          this->set_coords();
          //this->_slip_asm.assemble(this->filter.at<1>(), this->_trafo_space);
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
            this->compute_h();
            compute_grad_h();
          }
        }

        /**
         * \brief Computes the mesh concentration function for each cell
         *
         * This means it computes
         * \f[
         *   \forall K \in \mathcal{T}_h: \underbrace{c(K)}_{ =: \mathtt{conc(cell)}} ~ \mathrm{and} ~ \underbrace{\sum_{K \in \mathcal{T}_h} c(K)}_{\mathtt{\_sum\_conc}}
         * \f]
         **/
        void compute_conc()
        {
          // Total number of cells in the mesh
          Index ncells(this->get_mesh()->get_num_entities(ShapeType::dimension));
          // Index set for local/global numbering
          auto& idx = this->get_mesh()->template get_index_set<ShapeType::dimension,0>();

          _sum_conc = CoordType(0);
          for(Index cell(0); cell < ncells; ++cell)
          {
            // This will be the average of the levelset values at the vertices
            CoordType tmp(0);

            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
              tmp += this->_lvlset_vtx_vec(idx(cell,Index(j)));

            tmp = tmp/CoordType(Shape::FaceTraits<ShapeType,0>::count);

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
         * \returns
         * The mesh concentration for the given distance.
         *
         **/
        CoordType conc_val(CoordType dist)
        {
          return Math::pow(this->_r_adapt_reg + Math::abs(dist),this->_r_adapt_pow);
        }

        /**
         * \brief Computes the derivative of the concentration function
         *
         * \param[in] dist
         * The distance.
         *
         * \returns
         * The derivative of mesh concentration wrt. the distance.
         *
         **/
        CoordType conc_der(CoordType dist)
        {
          return this->_r_adapt_pow * Math::pow(this->_r_adapt_reg + Math::abs(dist),this->_r_adapt_pow - CoordType(1))*Math::signum(dist);
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
         * \param[out] grad_loc_
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
          grad_loc_.format(CoordType(0));

          // This will be the average of the levelset values at the vertices
          CoordType val(0);
          for(int i(0); i < Shape::FaceTraits<ShapeType,0>::count; ++i)
            val += lvlset_vals_(i);

          val = val/CoordType(Shape::FaceTraits<ShapeType,0>::count);

          for(int i(0); i < Shape::FaceTraits<ShapeType,0>::count; ++i)
            grad_loc_[i] = conc_der(val)/CoordType(Shape::FaceTraits<ShapeType,0>::count) * lvlset_grad_vals_[i];
        }

        /**
         * \copydoc Baseclass::compute_lambda()
         *
         *  In this case,
         *  \f$ \forall K_k \in \mathcal{T}_h: \lambda(K_k) = \frac{c(K_k)}{\sum_{K \in \mathcal{T}_h c(K)}}. \f$
         *
         **/
        virtual void compute_lambda() override
        {
          compute_conc();
          for(Index cell(0); cell < this->get_mesh()->get_num_entities(ShapeType::dimension); ++cell)
          {
            this->_lambda(cell, this->_conc(cell)/this->_sum_conc);
          }
        }

        /**
         * \brief Computes the local gradient of the optimal scales
         *
         * The optimal scales \f$ h \f$ depend on the concentration \f$ c \f$ by the relation
         * \f[
         *   \forall K_k \in \mathcal{T}_h: h_k = \left( \frac{c (K_k)}{\sum_{l=1}^N c(K_l) \sum_{l=1}^N
         *   \mathrm{det} \nabla R_{T,l}(\phi)} \right)^{\frac{1}{d}}
         * \f]
         *
         * The concentration function in turn indirectly depends on the vertex locations
         * \f$ \{ x_i : i=1, \dots, n \} \f$ via
         * \f[ c(K_k) = (\alpha + \sum_{x_j \in K_k} \varphi(x_j))^\beta, \f]
         * so that we have to take this dependency into account for the full gradient. Define
         * \f[
         *   s_d := \sum_{l=1}^N \mathrm{det} \nabla R_{T,l}(\Phi), s_c :=  \frac{c (K_k)}{\sum_{l=1}^N c(K_l)}.
         * \f]
         * So for each \f$ x_j \f$ we arrive at

         * \f{align*}
         *   \frac{\partial h(K_k)}{\partial x_j} & = \frac{1}{d} \left( \frac{c(K_k)}{s_d s_c} \right)^
         *   {\frac{1}{d}-1} \left[ c(K_k) ( \frac{\partial s_d}{\partial x_j} s_c + s_d
         *   \frac{\partial s_c}{\partial x_j} ) + \frac{\partial c(K_k)}{\partial x_j} s_d s_c \right]
         *   (s_d s_c)^{-2}
         * \f}
         *
         **/
        void compute_grad_h()
        {
          _grad_h.format();
          VectorType grad_sum_det(this->get_mesh()->get_num_entities(0), CoordType(0));
          CoordType sum_det = H_EvalType_::compute_sum_det(this->_coords, this->_trafo);
          H_EvalType_::compute_grad_sum_det(grad_sum_det, this->_coords, this->_trafo);
          compute_grad_conc();

          // Index set for local/global numbering
          auto& idx = this->get_mesh()->template get_index_set<ShapeType::dimension,0>();
          // This will hold the levelset values at the mesh vertices for one element
          FEAST::Tiny::Vector<CoordType, Shape::FaceTraits<ShapeType,0>::count> lvlset_vals;
          // This will hold the levelset gradient values for one element for passing to other routines
          FEAST::Tiny::Matrix<CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> lvlset_grad_vals;
          // This will hold the local gradient values for one element for passing to other routines
          FEAST::Tiny::Matrix<CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> grad_loc(0);
          // This will hold the computed local gradient values for one element for copy assigning to the blocked
          // datatype
          FEAST::Tiny::Vector<CoordType, MeshType::world_dim*Shape::FaceTraits<ShapeType,0>::count> tmp(0);

          CoordType exponent = CoordType(1)/CoordType(MeshType::world_dim) - CoordType(1);

          for(Index cell(0); cell < this->get_mesh()->get_num_entities(ShapeType::dimension); ++cell)
          {
            grad_loc.format();
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              Index i(idx(cell, Index(j)));
              // Get levelset
              lvlset_vals(j) = this->_lvlset_vec(i);
              for(int d(0); d < MeshType::world_dim; ++d)
                // Get levelset gradient
                lvlset_grad_vals(j,d) = this->_lvlset_grad_vtx_vec[d](i);
            }

            compute_grad_conc_local(grad_loc, lvlset_vals, lvlset_grad_vals);

            for(int d(0); d < MeshType::world_dim; ++d)
            {
              for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
              {
                Index i(idx(cell, Index(j)));

                tmp(j*MeshType::world_dim +d) =
                  CoordType(1)/CoordType(MeshType::world_dim)*Math::pow(_conc(cell)/_sum_conc*sum_det,exponent)
                  *( _conc(cell)*(grad_sum_det(i)(d)*_sum_conc + sum_det*_grad_conc(i)(d) )
                      + grad_loc(j,d) * sum_det *_sum_conc)
                  / Math::sqr(_sum_conc*sum_det);
              }
            }
            _grad_h(cell, tmp + _grad_h(cell));

          }

        } // compute_grad_h

        /**
         * \brief Computes the gradient of the sum of all mesh concentrations
         *
         * \f[
         *   \mathrm{grad\_conc}(k,i) = \frac{\partial}{\partial x_j} \sum_{l=1}^N c(K_l),
         * \f]
         * where \f$ i \f$ is the global index of the local vertex \f$ j \f$.
         *
         **/
        void compute_grad_conc()
        {
          // Clear the old gradient
          _grad_conc.format();

          // Index set for local/global numbering
          auto& idx = this->get_mesh()->template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector <CoordType,MeshType::world_dim> h;
          // This will hold the local gradient for one element for passing to other routines
          FEAST::Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> grad_loc;
          // This will hold the levelset values at the mesh vertices for one element
          FEAST::Tiny::Vector <CoordType, Shape::FaceTraits<ShapeType,0>::count> lvlset_vals;
          // This will hold the levelset gradient values for one element for passing to other routines
          FEAST::Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> lvlset_grad_vals;

          // Compute the functional value for each cell
          for(Index cell(0); cell < this->get_mesh()->get_num_entities(ShapeType::dimension); ++cell)
          {
            // Collect levelset and levelset grad values
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              // Global vertex/dof index
              Index i(idx(cell, Index(j)));
              // Get levelset
              lvlset_vals(j) = this->_lvlset_vec(i);
              for(int d(0); d < MeshType::world_dim; ++d)
                // Get levelset gradient
                lvlset_grad_vals(j,d) = this->_lvlset_grad_vtx_vec[d](i);
            }

            // Compute gradient of the concentration on this cell
            compute_grad_conc_local(grad_loc, lvlset_vals, lvlset_grad_vals);

            // Add local contributions to global gradient vector
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              Index i(idx(cell, Index(j)));
              this->_grad_conc(i, _grad_conc(i) + grad_loc[j]);
            }
          }

        } // compute_grad_conc()

        /// \copydoc RumpfSmootherBase::compute_gradient()
        virtual void compute_gradient() override
        {
          // Clear gradient vector
          this->_grad.format();

          // Total number of cells in the mesh
          Index ncells(this->get_mesh()->get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->get_mesh()->template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix<CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector<CoordType, MeshType::world_dim> h;
          // This will hold the local gradient for one element for passing to other routines
          FEAST::Tiny::Matrix<CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> grad_loc;
          // This will hold the levelset values at the mesh vertices for one element
          FEAST::Tiny::Vector<CoordType, Shape::FaceTraits<ShapeType,0>::count> lvlset_vals;
          // This will hold the levelset gradient values for one element for passing to other routines
          FEAST::Tiny::Matrix<CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> lvlset_grad_vals;

          // Compute contribution from each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            grad_loc.format();
            h = this->_h(cell);

            // Collect levelset and levelset grad values
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              Index i(idx(cell, Index(j)));
              // Get levelset
              lvlset_vals(j) = this->_lvlset_vec(i);
              // Get local coordinates
              x[j] = this->_coords(i);
              for(int d(0); d < MeshType::world_dim; ++d)
              {
                // Get levelset gradient
                lvlset_grad_vals[j](d) = this->_lvlset_grad_vtx_vec[d](i);
              }
            }

            // Compute local gradient and save it to grad_loc
            this->_functional.compute_local_grad(x, h, grad_loc);
            // Add the contribution from the dependence of h on the vertex coordinates
            this->_functional.add_grad_h_part(grad_loc, x, h, this->_grad_h(cell));
            // This does not contain the weighting yet
            grad_loc *= this->_mu(cell);
            // Add levelset penalty term, which is not weighted with mu
            this->_lvlset_functional.add_lvlset_penalty_grad(lvlset_vals, lvlset_grad_vals, grad_loc, this->lvlset_constraint_last);

            // Add local contributions to global gradient vector
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              // Global vertex/dof index
              Index i(idx(cell,Index(j)));
              this->_grad(i, this->_grad(i) + grad_loc[j]);
            }
          }

          this->_filter.filter_cor(this->_grad);

        } // compute_gradient
  }; // class RumpfSmootherLevelsetConcAnalytic

} // namespace Meshopt
} // namespace FEAST
#endif // KERNEL_MESHOPT_RUMPF_SMOOTHER_CONC_HPP
