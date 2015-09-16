#pragma once
#ifndef KERNEL_MESHOPT_RUMPF_SMOOTHER_LVLSET_Q1HACK_HPP
#define KERNEL_MESHOPT_RUMPF_SMOOTHER_LVLSET_Q1HACK_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/meshopt/rumpf_smoother_lvlset.hpp>

namespace FEAST
{
  namespace Meshopt
  {
    /**
     * \brief RumpfSmoother variant for using levelset information together with the Q1 hack
     *
     * \copydoc RumpfSmootherLevelsetAnalytic
     *
     * This variant is for Q1 transformations only and uses an idea by Steffen Basting: Split every hypercube
     * into simplices and evaluate the Rumpf functional for the P1 transformation of those. Do this for every
     * possible way of splitting the hypercube into simplices.
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
      typename H_EvalType_ = H_EvaluatorQ1Hack<TrafoType_, typename TrafoType_::MeshType::CoordType>
    >
    class RumpfSmootherLevelsetAnalyticQ1Hack :
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
        /// Type for the Rumpf functional
        typedef FunctionalType_ FunctionalType;
        /// Type for the levelset part of Rumpf functional
        typedef LevelsetFunctionalType_ LevelsetFunctionalType;

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

        /// Type of the underlying mesh
        typedef typename TrafoType_::MeshType MeshType;
        /// Data type for the vertex coordinates
        typedef typename MeshType::CoordType CoordType;
        /// The shape of the mesh's cells
        typedef typename MeshType::ShapeType ShapeType;

        // The functional has to use Simplex<shape_dim> and the trafo Hypercube<shape_dim>
        //static_assert( std::is_same<ShapeType, Shape::Hypercube<ShapeType::dimension> >::value &&
        //    std::is_same<typename FunctionalType::ShapeType, Shape::Simplex<ShapeType::dimension> >::value,
        //    "ShapeTypes of the transformation / functional have to be Hypercube<d>/Simplex<d> for RumpfSmootherQ1Hack" );


        /// \copydoc BaseClass::RumpfSmootherLevelsetAnalytic()
        explicit RumpfSmootherLevelsetAnalyticQ1Hack(Geometry::RootMeshNode<MeshType>* rmn_,
        std::deque<String>& dirichlet_list_, std::deque<String>& slip_list_,
          TrafoType_& trafo_, FunctionalType_& functional_, LevelsetFunctionalType_& lvlset_functional_,
          bool align_to_lvlset_, bool r_adaptivity_,
          AnalyticFunctionType_& analytic_function_,
          AnalyticFunctionGrad0Type_& analytic_function_grad0_,
          AnalyticFunctionGrad1Type_& analytic_function_grad1_)
          : BaseClass(rmn_, dirichlet_list_, slip_list_,
            trafo_, functional_, lvlset_functional_,
            align_to_lvlset_, r_adaptivity_ ,
            analytic_function_, analytic_function_grad0_, analytic_function_grad1_)
            {
            }

        /// \copydoc BaseClass::compute_functional()
        virtual CoordType compute_functional() override
        {
          CoordType fval(0);
          // Total number of cells in the mesh
          Index ncells(this->_mesh.get_num_entities(ShapeType::dimension));

          // In 2d, each hypercube is split into 2 simplices and there are two possible permutations
          const int n_perms(4);
          const Index perm[4][3] =
          {
            {Index(0), Index(1), Index(2)},
            {Index(1), Index(3), Index(2)},
            {Index(0), Index(3), Index(2)},
            {Index(0), Index(1), Index(3)}
          };
          // Index set for local/global numbering
          auto& idx = this->_mesh.template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector <CoordType, MeshType::world_dim> h;
          // This will hold the levelset values at the mesh vertices for one element
          FEAST::Tiny::Vector <CoordType, Shape::FaceTraits<ShapeType,0>::count> lvlset_vals;

          // Reset last levelset constraint
          this->lvlset_constraint_last = CoordType(0);
          CoordType lvlset_penalty(0);

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            h = this->_h(cell);
            // For every permutation, compute the terms wrt. the Simplex functional on the cell defined by the
            // vertices we pick through the permutation
            for(int p(0); p < n_perms; ++p)
            {
              for(int j(0); j < 3; j++)
                x[j] = this->_coords(idx(cell,perm[p][j]));

              fval += this->_mu(cell)*this->_functional.compute_local_functional(x,h);
            } // permutations

            // The levelset penalty has to be computed for the hypercube case
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; j++)
              lvlset_vals(j) = this->_lvlset_vec(idx(cell,Index(j)));

            lvlset_penalty = this->_lvlset_functional.compute_lvlset_penalty(lvlset_vals);
            // Add local penalty term to global levelset constraint
            this->lvlset_constraint_last += lvlset_penalty;
          } // cell

          fval += this->_lvlset_functional.fac_lvlset/CoordType(2)*Math::sqr(this->lvlset_constraint_last);

          return fval;
        } // compute_functional

        /*
         * This makes the baseclass' compute_functional available, as one might want to compute the functional value
         * without the levelset contribution and now the compiler no longer complains about overloaded virtual
         */
        using BaseClass::compute_functional;
        /**
         * \copydoc RumpfSmootherLevelset::compute_functional()
         *
         * \param[in] func_norm
         * The contribution of the Frobenius norm for each cell
         *
         * \param[in] func_det
         * The contribution of the det term for each cell
         *
         * \param[in] func_rec_det
         * The contribution of the 1/det term for each cell
         *
         * \param[in] func_lvlset
         * The contribution of the levelset term for each cell
         *
         * Debug variant that saves the different contributions for each cell.
         **/
        virtual CoordType compute_functional(CoordType* func_norm, CoordType* func_det,
        CoordType* func_rec_det, CoordType* func_lvlset) override
        {
          CoordType fval(0);
          // Total number of cells in the mesh
          Index ncells(this->_mesh.get_num_entities(ShapeType::dimension));

          // In 2d, each hypercube is split into 2 simplices and there are two possible permutations
          const int n_perms(4);
          const Index perm[4][3] =
          {
            {Index(0), Index(1), Index(2)},
            {Index(1), Index(3), Index(2)},
            {Index(0), Index(3), Index(2)},
            {Index(0), Index(1), Index(3)}
          };
          // Index set for local/global numbering
          auto& idx = this->_mesh.template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector <CoordType, MeshType::world_dim> h;
          // This will hold the levelset values at the mesh vertices for one element
          FEAST::Tiny::Vector <CoordType, Shape::FaceTraits<ShapeType,0>::count> lvlset_vals;

          CoordType norm_A(0), det_A(0), rec_det_A(0), lvlset_penalty(0);

          CoordType func_norm_tot(0);
          CoordType func_det_tot(0);
          CoordType func_rec_det_tot(0);

          // Reset last levelset constraint
          this->lvlset_constraint_last = CoordType(0);

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            func_norm[cell] = CoordType(0);
            func_det[cell] = CoordType(0);
            func_rec_det[cell] = CoordType(0);

            h = this->_h(cell);
            // For every permutation, compute the terms wrt. the Simplex functional on the cell defined by the
            // vertices we pick through the permutation
            for(int p(0); p < n_perms; ++p)
            {
              for(int j(0); j < 3; ++j)
                x[j] = this->_coords(idx(cell,perm[p][j]));

              fval += this->_mu(cell)*this->_functional.compute_local_functional(x,h, norm_A, det_A, rec_det_A);

              func_norm[cell] += this->_mu(cell) * norm_A;
              func_det[cell] += this->_mu(cell) * det_A;
              func_rec_det[cell] += this->_mu(cell) * rec_det_A;
            }
            func_norm_tot += func_norm[cell];
            func_det_tot += func_det[cell];
            func_rec_det_tot += func_rec_det[cell];

            // The levelset penalty has to be computed for the hypercube case
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
              lvlset_vals(j) = this->_lvlset_vec(idx(cell,Index(j)));

            lvlset_penalty = this->_lvlset_functional.compute_lvlset_penalty(lvlset_vals);
            // Add local penalty term to global levelset constraint
            this->lvlset_constraint_last += lvlset_penalty;
            func_lvlset[cell] +=  lvlset_penalty;
          }

          // Scale levelset constraint correctly and add it to the functional value
          fval += this->_lvlset_functional.fac_lvlset/CoordType(2)*Math::sqr(this->lvlset_constraint_last);

          std::cout << "func_norm = " << scientify(func_norm_tot) << ", func_det = " << scientify(func_det_tot) <<
            ", func_rec_det = " << scientify(func_rec_det_tot) << ", func_lvlset = " <<
            scientify(this->_lvlset_functional.fac_lvlset/CoordType(2)*Math::sqr(this->lvlset_constraint_last)) << std::endl;


          return fval;
        } // compute_functional

        /// \copydoc BaseClass::compute_gradient()
        virtual void compute_gradient() override
        {
          // Total number of cells in the mesh
          Index ncells(this->_mesh.get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->_mesh.template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector <CoordType, MeshType::world_dim> h;
          // This will hold the local gradient for one element for passing to other routines
          FEAST::Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> grad_loc;
          // This will hold the levelset values at the mesh vertices for one element
          FEAST::Tiny::Vector <CoordType, Shape::FaceTraits<ShapeType,0>::count> lvlset_vals;
          // This will hold the levelset gradient values for one element for passing to other routines
          FEAST::Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> lvlset_grad_vals;

          // In 2d, each hypercube is split into 2 simplices and there are two possible permutations
          const int n_perms(4);
          const Index perm[4][3] =
          {
            {Index(0), Index(1), Index(2)},
            {Index(1), Index(3), Index(2)},
            {Index(0), Index(3), Index(2)},
            {Index(0), Index(1), Index(3)}
          };

          // Clear gradient vector
          this->_grad.format(CoordType(0));

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            h = this->_h(cell);
            // For every permutation, compute the terms wrt. the Simplex functional on the cell defined by the
            // vertices we pick through the permutation
            for(int p(0); p < n_perms; ++p)
            {
              // Get local coordinates
              for(int j(0); j < 3; ++j)
              {
                x[j] = this->_coords(idx(cell,perm[p][j]));
                // Get levelset gradient
              }

              this->_functional.compute_local_grad(x, h, grad_loc);

              // Add local contributions to global gradient vector
              for(int j(0); j < 3; ++j)
              {
                // Global vertex/dof index
                Index i(idx(cell,perm[p][j]));
                Tiny::Vector<CoordType, MeshType::world_dim, MeshType::world_dim> tmp(this->_grad(i));
                tmp += grad_loc[j];

                this->_grad(i,tmp);
              }

            } // permutations

            grad_loc.format();

            // The levelset penalty has to be computed for the hypercube case
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; j++)
            {
              Index i(idx(cell, Index(j)));

              lvlset_vals(j) = this->_lvlset_vec(i);

              for(int d(0); d < MeshType::world_dim; ++d)
                lvlset_grad_vals(j,d) = this->_lvlset_grad_vtx_vec[d](i);

            }
            // Add levelset penalty term, which is not weighted with lambda
            this->_lvlset_functional.add_lvlset_penalty_grad(lvlset_vals, lvlset_grad_vals, grad_loc, this->lvlset_constraint_last);

            // Add local contributions to global gradient vector
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              // Global vertex/dof index
              Index i(idx(cell,Index(j)));
              Tiny::Vector<CoordType, MeshType::world_dim, MeshType::world_dim> tmp(this->_grad(i));
              tmp += grad_loc[j];

              this->_grad(i,tmp);
            }
          }

          this->_filter.filter_cor(this->_grad);

        } // compute_gradient

    }; // class RumpfSmootherLevelsetAnalyticQ1Hack

  } // namespace Meshopt
} // namespace FEAST
#endif // KERNEL_MESHOPT_RUMPF_SMOOTHER_LVLSET_Q1HACK_HPP
