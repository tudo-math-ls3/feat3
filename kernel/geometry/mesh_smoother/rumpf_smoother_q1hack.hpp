#pragma once
#ifndef KERNEL_GEOMETRY_RUMPF_SMOOTHER_Q1HACK_HPP
#define KERNEL_GEOMETRY_RUMPF_SMOOTHER_Q1HACK_HPP 1

#include <kernel/geometry/mesh_smoother/rumpf_smoother.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /**
     * \copydoc RumpfSmoother
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
      typename DataType_,
      typename MemType_,
      typename TrafoType_,
      typename FunctionalType_
    >
    class RumpfSmootherQ1Hack :
      public RumpfSmootherBase
      <
        DataType_,
        MemType_,
        TrafoType_,
        FunctionalType_,
        H_EvaluatorQ1Hack<TrafoType_, DataType_>
      >
    {
      public:
        /// Our datatype
        typedef DataType_ DataType;
        /// Memory architecture
        typedef MemType_ MemType;
        /// Transformation type
        typedef TrafoType_ TrafoType;
        /// Our functional type
        typedef FunctionalType_ FunctionalType;
        /// Meshsize evaluator
        typedef H_EvaluatorQ1Hack<TrafoType_, DataType_> H_EvalType;
        /// ShapeType
        typedef typename TrafoType::ShapeType ShapeType;
        /// Mesh of said ShapeType
        typedef Geometry::ConformalMesh<ShapeType, ShapeType::dimension, ShapeType::dimension, DataType> MeshType;
        /// Who's my daddy?
        typedef RumpfSmootherBase<DataType_, MemType_, TrafoType, FunctionalType_, H_EvalType > BaseClass;
        // The functional has to use Simplex<shape_dim> and the trafo Hypercube<shape_dim>
        //static_assert( std::is_same<ShapeType, Shape::Hypercube<ShapeType::dimension> >::value &&
        //std::is_same<typename FunctionalType::ShapeType, Shape::Simplex<ShapeType::dimension> >::value, "ShapeTypes of the transformation / functional have to be Hypercube<d>/Simplex<d> for RumpfSmootherQ1Hack" );

        /// \copydoc RumpfSmoother()
        explicit RumpfSmootherQ1Hack( const TrafoType& trafo_, FunctionalType& functional_)
          : BaseClass(trafo_, functional_)
          {
          }

        /// \copydoc BaseClass::compute_functional()
        virtual DataType compute_functional()
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
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, 3> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector<DataType_, MeshType::world_dim> h;

          // Compute the functional value for each cell...
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
            }
          }

          return fval;
        } // compute_functional

        /// \copydoc BaseClass::compute_functional(DataType_*, DataType_*, DataType_*)
        virtual DataType compute_functional( DataType_* func_norm, DataType_* func_det, DataType_* func_rec_det )
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
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, 3> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector<DataType_, MeshType::world_dim> h;

          DataType_ norm_A(0), det_A(0), rec_det_A(0);

          DataType_ func_norm_tot(0);
          DataType_ func_det_tot(0);
          DataType_ func_rec_det_tot(0);

          // Compute the functional value for each cell...
          for(Index cell(0); cell < ncells; ++cell)
          {
            func_norm[cell] = DataType(0);
            func_det[cell] = DataType(0);
            func_rec_det[cell] = DataType(0);

            //std::cout << "cell " << cell << std::endl;
            // ... and for each simplex in all permutations
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
              fval += this->_lambda(cell)*this->_functional.compute_local_functional(x,h, norm_A, det_A, rec_det_A);

              func_norm[cell] += this->_lambda(cell) * norm_A;
              func_det[cell] += this->_lambda(cell) * det_A;
              func_rec_det[cell] += this->_lambda(cell) * rec_det_A;
            }
            func_norm_tot += func_norm[cell];
            func_det_tot += func_det[cell];
            func_rec_det_tot += func_rec_det[cell];
          }

          std::cout << "func_norm = " << scientify(func_norm_tot) << ", func_det = " << scientify(func_det_tot) <<
            ", func_rec_det = " << scientify(func_rec_det_tot) << std::endl;

          return fval;
        } // compute_functional

        /// \brief Computes the gradient of the functional with regard to the nodal coordinates.
        /// \copydoc BaseClass::compute_gradient()
        virtual void compute_gradient()
        {
          // Total number of cells in the mesh
          Index ncells(this->_mesh.get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->_mesh.template get_index_set<ShapeType::dimension,0>();

          // In 2d, each hypercube is split into 2 simplices and there are two possible permutations
          const int n_perms(4);
          const Index perm[4][3] =
          {
            {Index(0), Index(1), Index(2)},
            {Index(1), Index(3), Index(2)},
            {Index(0), Index(3), Index(2)},
            {Index(0), Index(1), Index(3)}
          };
          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <DataType_, MeshType::world_dim, 3> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector<DataType_, MeshType::world_dim> h;
          // This will hold the local gradient for one element for passing to other routines
          FEAST::Tiny::Matrix<DataType_, MeshType::world_dim, 3 > grad_loc;

          // Clear gradient vector
          for(Index i(0); i < this->_world_dim*this->_nk; ++i)
            this->_grad[i] = DataType_(0);

          // Compute the functional value for each cell...
          for(Index cell(0); cell < ncells; ++cell)
          {
            // ... and for each simplex in all permutations
            for(Index p(0); p < n_perms; ++p)
            {
              for(Index d(0); d < this->_world_dim; ++d)
              {
                h(d) = this->_h[d](cell);
                // Get local coordinates
                for(Index j(0); j < 3; ++j)
                  x(d,j) = this->_coords[d](idx(cell,perm[p][j]));
              }

              this->_functional.compute_local_grad(x, h, grad_loc);

              for(Index d(0); d < this->_world_dim; ++d)
              {
                // Add local contributions to global gradient vector
                for(Index j(0); j < 3; ++j)
                  this->_grad[d*this->_nk + idx(cell,perm[p][j])] += this->_lambda(cell)*grad_loc(d,j);
              }
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

        /// \copydoc MeshSmoother::optimise()
        virtual void optimise()
        {
          ALGLIBWrapper<RumpfSmootherQ1Hack<DataType_, MemType_, TrafoType_, FunctionalType_>>::minimise_functional_cg(*this);
          // Important: Copy back the coordinates the mesh optimiser changed to the original mesh.
          this->set_coords();
        }

    }; // class RumpfSmootherQ1Hack

  } // namespace Geometry
} // namespace FEAST
#endif // KERNEL_GEOMETRY_RUMPF_SMOOTHER_HPP
