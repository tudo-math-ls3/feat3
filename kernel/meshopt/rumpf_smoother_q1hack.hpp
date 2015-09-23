#pragma once
#ifndef KERNEL_MESHOPT_RUMPF_SMOOTHER_Q1HACK_HPP
#define KERNEL_MESHOPT_RUMPF_SMOOTHER_Q1HACK_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/meshopt/rumpf_smoother.hpp>

namespace FEAST
{
  namespace Meshopt
  {
    /**
     * \copydoc RumpfSmoother
     *
     * This variant is for Q1 transformations only and uses an idea by Steffen Basting: Split every hypercube
     * into simplices and evaluate the Rumpf functional for the P1 transformation of those. Do this for every
     * possible way of splitting the hypercube into simplices.
     *
     * \verbatim
        Hypercube:    Splitting 1:   Splitting 2
         2-----3       2-----3        2-----3
         |     |       | \   |        |   / |
         |     |       |   \ |        | /   |
         0-----1       0-----1        0-----1
     * \endverbatim
     * Splitting 1 gives permutations 1 and 2, splitting 2 gives permutations 3 and 4.
     *
     * \author Jordi Paul
     *
     **/
    template
    <
      typename TrafoType_,
      typename FunctionalType_
    >
    class RumpfSmootherQ1Hack :
      public RumpfSmootherBase
      <
        TrafoType_,
        FunctionalType_,
        H_EvaluatorQ1Hack<TrafoType_, typename TrafoType_::MeshType::CoordType>
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
        /// Meshsize evaluator
        typedef H_EvaluatorQ1Hack<TrafoType, CoordType> H_EvalType;

        /// Only Mem::Main is supported atm
        typedef Mem::Main MemType;
        /// We always use Index for now
        typedef Index IndexType;

        /// Our base class
        typedef RumpfSmootherBase<TrafoType_, FunctionalType_, H_EvalType> BaseClass;
        /// Own type for passing to ALGLIB
        typedef RumpfSmootherQ1Hack<TrafoType_, FunctionalType_> MyType;

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
        /// Since the functional contains a ShapeType, these have to be the same
        static_assert( std::is_same<ShapeType, typename FunctionalType::ShapeType>::value, "ShapeTypes of the transformation / functional have to agree" );

        /// \copydoc RumpfSmootherBase()
        explicit RumpfSmootherQ1Hack(Geometry::RootMeshNode<MeshType>* rmn_,
        std::deque<String>& dirichlet_list_, std::deque<String>& slip_list_, FunctionalType_& functional_)
          : BaseClass(rmn_, dirichlet_list_, slip_list_, functional_)
          {
          }

        /// \copydoc BaseClass::compute_functional()
        virtual CoordType compute_functional()
        {
          CoordType fval(0);
          // Total number of cells in the mesh
          Index ncells(this->get_mesh()->get_num_entities(ShapeType::dimension));

          // In 2d, each hypercube is split into 2 simplices and there are two possible permutations
          const int n_perms(4);
          const int perm[4][3] =
          {
            {0, 1, 2},
            {1, 3, 2},
            {0, 3, 2},
            {0, 1, 3}
          };

          // Scale everything with this to be consistent
          CoordType scal(CoordType(1)/CoordType(n_perms));

          // Index set for local/global numbering
          auto& idx = this->get_mesh()->template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <CoordType, 3, MeshType::world_dim> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector<CoordType, MeshType::world_dim> h;

          // Compute the functional value for each cell...
          for(Index cell(0); cell < ncells; ++cell)
          {
            h = this->_h(cell);
            // ... and for each simplex in all permutations
            for(int p(0); p < n_perms; ++p)
            {
              for(int d(0); d < MeshType::world_dim; d++)
              {
                for(int j(0); j < 3; j++)
                  x[j] = this->_coords(idx(cell,Index(perm[p][j])));
              }
              fval += this->_mu(cell)*this->_functional.compute_local_functional(x,h);
            }
          }

          return scal*fval;
        } // compute_functional

        /**
         * \copydoc RumpfSmoother::compute_functional(CoordType*, CoordType*, CoordType*)
         *
         */
        virtual CoordType compute_functional( CoordType* func_norm, CoordType* func_det, CoordType* func_rec_det )
        {
          CoordType fval(0);
          // Total number of cells in the mesh
          Index ncells(this->get_mesh()->get_num_entities(ShapeType::dimension));

          // In 2d, each hypercube is split into 2 simplices and there are two possible permutations
          const int n_perms(4);
          const int perm[4][3] =
          {
            {0, 1, 2},
            {1, 3, 2},
            {0, 3, 2},
            {0, 1, 3}
          };

          // Scale everything with this to be consistent
          CoordType scal(CoordType(1)/CoordType(n_perms));

          // Index set for local/global numbering
          auto& idx = this->get_mesh()->template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <CoordType, 3, MeshType::world_dim> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector<CoordType, MeshType::world_dim> h;

          CoordType norm_A(0), det_A(0), rec_det_A(0);

          CoordType func_norm_tot(0);
          CoordType func_det_tot(0);
          CoordType func_rec_det_tot(0);

          // Compute the functional value for each cell...
          for(Index cell(0); cell < ncells; ++cell)
          {
            func_norm[cell] = CoordType(0);
            func_det[cell] = CoordType(0);
            func_rec_det[cell] = CoordType(0);

            h = this->_h(cell);
            // ... and for each simplex in all permutations
            for(int p(0); p < n_perms; ++p)
            {
              //std::cout << " permutation " << p << std::endl;
              for(int d(0); d < MeshType::world_dim; ++d)
              {
                for(int j(0); j < 3; ++j)
                  x[j] = this->_coords(idx(cell,Index(perm[p][j])));
              }
              fval += this->_mu(cell)*this->_functional.compute_local_functional(x,h, norm_A, det_A, rec_det_A);

              func_norm[cell] += scal*this->_mu(cell) * norm_A;
              func_det[cell] += scal*this->_mu(cell) * det_A;
              func_rec_det[cell] += scal*this->_mu(cell) * rec_det_A;
            }
            func_norm_tot += func_norm[cell];
            func_det_tot += func_det[cell];
            func_rec_det_tot += func_rec_det[cell];
          }

          fval *= scal;

          std::cout << "func_norm = " << scientify(func_norm_tot) << ", func_det = " << scientify(func_det_tot) <<
            ", func_rec_det = " << scientify(func_rec_det_tot) << std::endl;

          return fval;
        } // compute_functional

        /// \brief Computes the gradient of the functional with regard to the nodal coordinates.
        /// \copydoc BaseClass::compute_gradient()
        virtual void compute_gradient()
        {
          // Total number of cells in the mesh
          Index ncells(this->get_mesh()->get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->get_mesh()->template get_index_set<ShapeType::dimension,0>();

          // In 2d, each hypercube is split into 2 simplices and there are two possible permutations
          const int n_perms(4);
          const int perm[4][3] =
          {
            {0, 1, 2},
            {1, 3, 2},
            {0, 3, 2},
            {0, 1, 3}
          };

          // Scale everything with this to be consistent
          CoordType scal(CoordType(1)/CoordType(n_perms));

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix<CoordType, 3, MeshType::world_dim> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector<CoordType, MeshType::world_dim> h;
          // This will hold the local gradient for one element for passing to other routines
          FEAST::Tiny::Matrix<CoordType, 3, MeshType::world_dim> grad_loc;

          // Clear gradient vector
          this->_grad.format();

          // Compute the functional value for each cell...
          for(Index cell(0); cell < ncells; ++cell)
          {
            h = this->_h(cell);
            // ... and for each simplex in all permutations
            for(int p(0); p < n_perms; ++p)
            {
              // Get local coordinates
              for(int j(0); j < 3; ++j)
                x[j] = this->_coords(idx(cell,Index(perm[p][j])));

              this->_functional.compute_local_grad(x, h, grad_loc);

              // Add local contributions to global gradient vector
              for(int j(0); j < 3; ++j)
              {
                Index i(idx(cell,Index(perm[p][j])));
                Tiny::Vector<CoordType, MeshType::world_dim, MeshType::world_dim> tmp(this->_grad(i));
                tmp += (scal*this->_mu(cell))*grad_loc[j];

                this->_grad(i,tmp);
              }
            }
          }

          this->_filter.filter_cor(this->_grad);

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

  } // namespace Meshopt
} // namespace FEAST
#endif // KERNEL_MESHOPT_RUMPF_SMOOTHER_HPP
