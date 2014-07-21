#pragma once
#ifndef KERNEL_GEOMETRY_RUMPF_SMOOTHER_LVLSET_Q1HACK_HPP
#define KERNEL_GEOMETRY_RUMPF_SMOOTHER_LVLSET_Q1HACK_HPP 1

#include <kernel/geometry/mesh_smoother/rumpf_smoother_lvlset.hpp>

namespace FEAST
{
  namespace Geometry
  {

    /**
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
      typename SpaceType_,
      typename FunctionalType_,
      typename TrafoType_,
      typename DataType_,
      typename MemType_,
      typename H_EvalType_ = H_EvaluatorQ1Hack<TrafoType_, DataType_>
    >
    class RumpfSmootherLevelsetAnalyticQ1Hack :
      public RumpfSmootherLevelsetAnalytic
      <
        AnalyticFunctionType_,
        AnalyticFunctionGrad0Type_,
        AnalyticFunctionGrad1Type_,
        SpaceType_,
        FunctionalType_,
        TrafoType_,
        DataType_,
        MemType_,
        H_EvalType_
      >
    {
      public:
        /// Our functional type
        typedef FunctionalType_ FunctionalType;
        /// Transformation type
        typedef TrafoType_ TrafoType;
        /// Our datatype
        typedef DataType_ DataType;
        /// Memory architecture
        typedef MemType_ MemType;
        /// ShapeType
        typedef typename TrafoType::ShapeType ShapeType;
        /// Mesh of said ShapeType
        typedef Geometry::ConformalMesh<ShapeType, ShapeType::dimension, ShapeType::dimension, DataType> MeshType;
        // The functional has to use Simplex<shape_dim> and the trafo Hypercube<shape_dim>
        //static_assert( std::is_same<ShapeType, Shape::Hypercube<ShapeType::dimension> >::value &&
        //    std::is_same<typename FunctionalType::ShapeType, Shape::Simplex<ShapeType::dimension> >::value,
        //    "ShapeTypes of the transformation / functional have to be Hypercube<d>/Simplex<d> for RumpfSmootherQ1Hack" );

        /// Baseclass
        typedef RumpfSmootherLevelsetAnalytic
        <
          AnalyticFunctionType_,
          AnalyticFunctionGrad0Type_,
          AnalyticFunctionGrad1Type_,
          SpaceType_,
          FunctionalType_,
          TrafoType_,
          DataType_,
          MemType_,
          H_EvalType_
        > BaseClass;

        /// \copydoc BaseClass::RumpfSmootherLevelsetAnalytic()
        explicit RumpfSmootherLevelsetAnalyticQ1Hack(
          TrafoType& trafo_,
          FunctionalType_& functional_,
          bool align_to_lvlset_,
          bool from_original_,
          bool r_adaptivity_,
          AnalyticFunctionType_& analytic_function_,
          AnalyticFunctionGrad0Type_& analytic_function_grad0_,
          AnalyticFunctionGrad1Type_& analytic_function_grad1_)
          : BaseClass(
            trafo_,
            functional_,
            align_to_lvlset_,
            from_original_,
            r_adaptivity_ ,
            analytic_function_,
            analytic_function_grad0_,
            analytic_function_grad1_)
            {
            }

        /// \copydoc BaseClass::compute_functional()
        virtual DataType compute_functional() override
        {
          DataType_ fval(0);
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
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector <DataType_,MeshType::world_dim> h;
          // This will hold the levelset values at the mesh vertices for one element
          FEAST::Tiny::Vector <DataType_, Shape::FaceTraits<ShapeType,0>::count> lvlset_vals;

          // Reset last levelset constraint
          this->lvlset_constraint_last = DataType(0);
          DataType lvlset_penalty(0);

          // Compute the functional value for each cell ...
          for(Index cell(0); cell < ncells; ++cell)
          {
            // ... and for each simplex in all permutations
            for(Index p(0); p < n_perms; ++p)
            {
              for(Index d(0); d < this->_world_dim; d++)
              {
                h(d) = this->_h[d](cell);
                for(Index j(0); j < 3; j++)
                  x(d,j) = this->_coords[d](idx(cell,perm[p][j]));
              }
              fval += this->_lambda(cell)*this->_functional.compute_local_functional(x,h);
            } // permutations

            // The levelset penalty has to be computed for the hypercube case
            for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; j++)
              lvlset_vals(j) = this->_lvlset_vec(idx(cell,j));

            lvlset_penalty = this->_functional.compute_lvlset_penalty(lvlset_vals);
            // Add local penalty term to global levelset constraint
            this->lvlset_constraint_last += lvlset_penalty;
          } // cell

          fval += this->_functional.fac_lvlset/DataType(2)*Math::sqr(this->lvlset_constraint_last);

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
         * \param[in] func_det2
         * The contribution of the 1/det term for each cell
         *
         * \param[in] func_lvlset
         * The contribution of the levelset term for each cell
         *
         * Debug variant that saves the different contributions for each cell.
         **/
        virtual DataType compute_functional(DataType_* func_norm, DataType_* func_det,
        DataType_* func_det2, DataType* func_lvlset) override
        {
          DataType_ fval(0);
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
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector <DataType_,MeshType::world_dim> h;
          // This will hold the levelset values at the mesh vertices for one element
          FEAST::Tiny::Vector <DataType_, Shape::FaceTraits<ShapeType,0>::count> lvlset_vals;

          DataType_ norm_A(0), det_A(0), det2_A(0), lvlset_penalty(0);

          DataType_ func_norm_tot(0);
          DataType_ func_det_tot(0);
          DataType_ func_det2_tot(0);

          // Reset last levelset constraint
          this->lvlset_constraint_last = DataType(0);

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            func_norm[cell] = DataType(0);
            func_det[cell] = DataType(0);
            func_det2[cell] = DataType(0);

            //std::cout << "cell " << cell << std::endl;
            for(Index p(0); p < n_perms; ++p)
            {
              //std::cout << " permutation " << p << std::endl;
              for(Index d(0); d < this->_world_dim; ++d)
              {
                h(d) = this->_h[d](cell);
                for(Index j(0); j < 3; ++j)
                {
                  x(d,j) = this->_coords[d](idx(cell,perm[p][j]));
                  //std::cout << perm[p][j] << ": " << scientify(x(d,j)) << " " ;
                }
                //std::cout << std::endl;
              }

              fval += this->_lambda(cell)*this->_functional.compute_local_functional(x,h, norm_A, det_A, det2_A);

              func_norm[cell] += this->_lambda(cell) * norm_A;
              func_det[cell] += this->_lambda(cell) * det_A;
              func_det2[cell] += this->_lambda(cell) * det2_A;
              func_norm_tot += func_norm[cell];
              func_det_tot += func_det[cell];
              func_det2_tot += func_det2[cell];
            }

            // The levelset penalty has to be computed for the hypercube case
            for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; j++)
              lvlset_vals(j) = this->_lvlset_vec(idx(cell,j));

            lvlset_penalty = this->_functional.compute_lvlset_penalty(lvlset_vals);
            // Add local penalty term to global levelset constraint
            this->lvlset_constraint_last += lvlset_penalty;
            func_lvlset[cell] +=  lvlset_penalty;
          }

          // Scale levelset constraint correctly and add it to the functional value
          fval += this->_functional.fac_lvlset/DataType(2)*Math::sqr(this->lvlset_constraint_last);

          std::cout << "func_norm = " << scientify(func_norm_tot) << ", func_det = " << scientify(func_det_tot) <<
            ", func_det2 = " << scientify(func_det2_tot) << ", func_lvlset = " <<
            scientify(this->_functional.fac_lvlset/DataType(2)*Math::sqr(this->lvlset_constraint_last)) << std::endl;

          return fval;
        } // compute_functional

        /// \copydoc BaseClass::compute_gradient()
        virtual void compute_gradient()
        {
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
          for(Index i(0); i < this->_world_dim*this->_nk; ++i)
            this->_grad[i] = DataType_(0);

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            for(Index p(0); p < n_perms; ++p)
            {
              for(Index d(0); d < this->_world_dim; ++d)
              {
                h(d) = this->_h[d](cell);
                // Get local coordinates
                for(Index j(0); j < 3; ++j)
                {
                  x(d,j) = this->_coords[d](idx(cell,perm[p][j]));
                  // Get levelset gradient
                }
              }

              this->_functional.compute_local_grad(x, h, grad_loc);

              // Add local contributions to global gradient vector
              for(Index d(0); d < this->_world_dim; ++d)
              {
                for(Index j(0); j < 3; ++j)
                  this->_grad[d*this->_nk + idx(cell,perm[p][j])] += this->_lambda(cell)*grad_loc(d,j);
              }
            } // permutations

            grad_loc.format();
            // The levelset penalty has to be computed for the hypercube case
            for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; j++)
            {
              lvlset_vals(j) = this->_lvlset_vec(idx(cell,j));

              for(Index d(0); d < this->_world_dim; ++d)
                lvlset_grad_vals(d,j) = this->_lvlset_grad_vtx_vec[d](idx(cell,j));

            }
            // Add levelset penalty term, which is not weighted with lambda
            this->_functional.add_lvlset_penalty_grad(lvlset_vals, lvlset_grad_vals, grad_loc, this->lvlset_constraint_last);
            for(Index d(0); d < this->_world_dim; ++d)
            {
              // Add local contributions to global gradient vector
              for(Index j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
                this->_grad[d*this->_nk + idx(cell,j)] += grad_loc(d,j);
            }
          }

          // Set gradient to 0 where bdry_id == -1
          // TODO: Convert this to proper filtering etc. and allow for more general BCs.
          for(Index i(0); i < this->_mesh.get_num_entities(0); ++i)
          {
            if(this->_bdry_id[i] == -1)
            {
              for(Index d(0); d < this->_world_dim; ++d)
                this->_grad[d*this->_nk + i] = DataType(0);
            }
          }

        } // compute_gradient

    }; // class RumpfSmootherLevelsetAnalytic

  } // namespace Geometry
} // namespace FEAST
#endif // KERNEL_GEOMETRY_RUMPF_SMOOTHER_HPP
