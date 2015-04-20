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
        /// And who am I?
        typedef RumpfSmootherQ1Hack<DataType_, MemType_, TrafoType_, FunctionalType_> MyType;
        /// Since the functional contains a ShapeType, these have to be the same
        static_assert( std::is_same<ShapeType, typename FunctionalType::ShapeType>::value, "ShapeTypes of the transformation and functional have to agree" );
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
              for(Index d(0); d < MeshType::world_dim; d++)
              {
                h(d) = this->_h[d](cell);
                for(Index j(0); j < 3; j++)
                  x(d,j) = this->_coords[d](idx(cell,perm[p][j]));
              }
              fval += this->_mu(cell)*this->_functional.compute_local_functional(x,h);
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
              for(Index d(0); d < MeshType::world_dim; ++d)
              {
                h(d) = this->_h[d](cell);
                for(Index j(0); j < 3; ++j)
                {
                  x(d,j) = this->_coords[d](idx(cell,perm[p][j]));
                  //std::cout << perm[p][j] << ": " << scientify(x(d,j)) << " " ;
                }
                //std::cout << std::endl;
              }
              fval += this->_mu(cell)*this->_functional.compute_local_functional(x,h, norm_A, det_A, rec_det_A);

              func_norm[cell] += this->_mu(cell) * norm_A;
              func_det[cell] += this->_mu(cell) * det_A;
              func_rec_det[cell] += this->_mu(cell) * rec_det_A;
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
          // Total number of vertices in the mesh
          Index nvertices(this->_mesh.get_num_entities(0));
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
          for(Index i(0); i < MeshType::world_dim*nvertices; ++i)
            this->_grad[i] = DataType_(0);

          // Compute the functional value for each cell...
          for(Index cell(0); cell < ncells; ++cell)
          {
            // ... and for each simplex in all permutations
            for(Index p(0); p < n_perms; ++p)
            {
              for(Index d(0); d < MeshType::world_dim; ++d)
              {
                h(d) = this->_h[d](cell);
                // Get local coordinates
                for(Index j(0); j < 3; ++j)
                  x(d,j) = this->_coords[d](idx(cell,perm[p][j]));
              }

              this->_functional.compute_local_grad(x, h, grad_loc);

              for(Index d(0); d < MeshType::world_dim; ++d)
              {
                // Add local contributions to global gradient vector
                for(Index j(0); j < 3; ++j)
                  this->_grad[d*nvertices + idx(cell,perm[p][j])] += this->_mu(cell)*grad_loc(d,j);
              }
            }
          }

          this->_filter_grad();

        } // compute_gradient

        /// \copydoc MeshSmoother::optimise()
        virtual void optimise()
        {
          int total_grad_evals(0);
          int total_iterations(0);
          int termination_type(0);

          ALGLIBWrapper<MyType>::minimise_functional_cg(total_grad_evals, total_iterations, termination_type,*this);
          std::cout << total_iterations << " mincg iterations, " << total_grad_evals << " grad evals, terminationtype was " << termination_type << std::endl;
          // Important: Copy back the coordinates the mesh optimiser changed to the original mesh.
          this->set_coords();
        }

    }; // class RumpfSmootherQ1Hack

  } // namespace Geometry
} // namespace FEAST
#endif // KERNEL_GEOMETRY_RUMPF_SMOOTHER_HPP
