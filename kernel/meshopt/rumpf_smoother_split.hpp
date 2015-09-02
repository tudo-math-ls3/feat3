#pragma once
#ifndef KERNEL_MESHOPT_RUMPF_SMOOTHER_SPLIT_HPP
#define KERNEL_MESHOPT_RUMPF_SMOOTHER_SPLIT_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/lafem/unit_filter_blocked.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/mesh_smoother/h_evaluator.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_smoother.hpp>
// ALGLIB includes
FEAST_DISABLE_WARNINGS
#include <thirdparty/alglib/cpp/src/optimization.h>
FEAST_RESTORE_WARNINGS
namespace FEAST
{
  namespace Meshopt
  {
    /**
     * \brief RumpfSmoother variant where only the Frobenius norm term is treated implicitly
     *
     * Both terms containing det are treated explicitly and there is a Richardson iteration around all of this.
     *
     * \note The FunctionalType_ must offer a compute_grad_norm() routine, which only computes the gradient of
     * Frobenius norm term.
     *
     */
    template
    <
      typename TrafoType_,
      typename FunctionalType_,
      typename H_EvalType_ = H_Evaluator<TrafoType_, typename TrafoType_::MeshType::CoordType_>
    >
    class RumpfSmootherSplit:
      public RumpfSmoother<TrafoType_, FunctionalType_, H_EvalType_>
    {
        /// Our base class
        typedef RumpfSmoother<TrafoType_, FunctionalType_, H_EvalType_> BaseClass;

        /// Type for the transformation
        typedef TrafoType_ TrafoType;

        /// The mesh the transformation is defined on
        typedef typename TrafoType::MeshType MeshType;
        /// The precision of the mesh coordinates
        typedef typename MeshType::CoordType CoordType;

        /// Only Mem::Main is supported atm
        typedef Mem::Main MemType;
        /// We always use the precision of the mesh
        typedef CoordType DataType;
        /// We always use Index for now
        typedef Index IndexType;

        /// Type for the functional
        typedef FunctionalType_ FunctionalType;
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
        // Since the functional contains a ShapeType, these have to be the same
        static_assert( std::is_same<ShapeType, typename FunctionalType::ShapeType>::value,
        "ShapeTypes of the transformation / functional have to agree" );
        /// Coordinates of the vertices, as they get changed in the optimisation process
        VectorType _coords_old[MeshType::world_dim];
        /// Relaxation parameter for the Richardson iteration.
        CoordType relaxation_parameter;


      public:
        /**
         * \brief Constructor
         *
         * \param[in] trafo_
         * Reference to the underlying transformation
         *
         * \param[in] functional_
         * Reference to the functional used
         *
         */
        explicit RumpfSmootherSplit(TrafoType& trafo_, FunctionalType& functional_)
          : BaseClass(trafo_, functional_)
          {
          }

        /**
         * \brief Constructor with filter
         *
         * This constructor takes a filter, so that boundary conditions can be specified.
         *
         * \param[in] trafo_
         * Reference to the underlying transformation
         *
         * \param[in] functional_
         * Reference to the functional used
         *
         */
        explicit RumpfSmootherSplit(TrafoType& trafo_, FunctionalType& functional_, FilterType& filter_)
          : BaseClass(trafo_, functional_, filter_)
          {
          }

        /// \brief Destructor
        virtual ~RumpfSmootherSplit()
        {
        };

        /**
         * \brief Computes the functional value on the current mesh.
         *
         * \returns
         * The functional value
         *
         **/
        virtual CoordType compute_functional()
        {
          CoordType_ fval(0);
          // Total number of cells in the mesh
          Index ncells(this->_mesh.get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->_mesh.template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <CoordType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;
          FEAST::Tiny::Matrix <CoordType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x_old;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector<CoordType_, MeshType::world_dim> h;

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            h = this->_h(cell);
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; j++)
            {
              Index i(idx(cell, Index(j)));
              x[j] = this->_coords(i);
              x_old[j] = this->_coords_old(i);
            }

            fval += this->_mu(cell) * (this->_functional._fac_norm*this->_functional.compute_norm_A(x,h)
                + this->_functional._fac_det*this->_functional.compute_det_A(x_old,h)
                +this->_functional._fac_rec_det*this->_functional.compute_rec_det_A(x_old,h));
          }

          return fval;
        } // compute_functional

        /**
         * \brief Computes the functional value on the current mesh.
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
         * \returns
         * The functional value
         *
         * Debug variant that saves the different contributions for each cell.
         **/
        virtual CoordType compute_functional( CoordType_* func_norm, CoordType_* func_det, CoordType_* func_rec_det )
        {
          CoordType_ fval(0);
          // Total number of cells in the mesh
          Index ncells(this->_mesh.get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->_mesh.template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <CoordType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector<CoordType_, MeshType::world_dim> h;

          CoordType_ norm_A(0), det_A(0), rec_det_A(0);

          CoordType_ func_norm_tot(0);
          CoordType_ func_det_tot(0);
          CoordType_ func_rec_det_tot(0);
          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            h = this->_h(cell);
            // Get local coordinates
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; j++)
              x[j] = this->_coords(idx(cell,Index(j)));

            // Scale local functional value with lambda
            fval += this->_mu(cell) * this->_functional.compute_local_functional(x,h, norm_A, det_A, rec_det_A);

            func_norm[cell] = this->_mu(cell) * norm_A;
            func_det[cell] = this->_mu(cell) * det_A;
            func_rec_det[cell] = this->_mu(cell) * rec_det_A;
            func_norm_tot += func_norm[cell];
            func_det_tot += func_det[cell];
            func_rec_det_tot += func_rec_det[cell];
          }

          std::cout << "fval = " << scientify(fval) << " func_norm = " << scientify(func_norm_tot)
          << ", func_det = " << scientify(func_det_tot) << ", func_rec_det = " << scientify(func_rec_det_tot) << std::endl;

          return fval;
        } // compute_functional

        /// \brief Computes the gradient of the functional with regard to the nodal coordinates.
        virtual void compute_gradient()
        {
          // Total number of vertices in the mesh
          Index nvertices(this->_mesh.get_num_entities(0));
          // Total number of cells in the mesh
          Index ncells(this->_mesh.get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->_mesh.template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <CoordType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector<CoordType_, MeshType::world_dim> h;
          // This will hold the local gradient for one element for passing to other routines
          FEAST::Tiny::Matrix<CoordType_, MeshType::world_dim, Shape::FaceTraits<ShapeType,0>::count> grad_loc;

          // Clear gradient vector
          this->_grad.format(CoordType_(0));

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            h = this->_h(cell);
            // Get local coordinates
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
              x[j] = this->_coords(idx(cell,Index(j)));

            this->_functional.compute_grad_norm(x, h, grad_loc);

            // Add local contributions to global gradient vector
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              Index i(idx(cell,Index(j)));
              Tiny::Vector<CoordType, MeshType::world_dim, MeshType::world_dim> tmp(this->_grad(i));
              tmp += this->_mu(cell)*grad_loc[j];

              this->_grad(i,tmp);
            }
          }

          this->_filter.filter_cor(this->_grad);

        } // compute_gradient

        /// \copydoc MeshSmoother::optimise()
        virtual void optimise()
        {
          CoordType diff(1);
          CoordType tol(Math::pow(Math::eps<CoordType_>(),CoordType(0.75) ));

          int total_grad_evals(0);
          int total_iterations(0);
          int termination_type(0);

          // coords_new = relaxation_parameter*coords + (1 - relaxation_parameter)*coords_old
          relaxation_parameter = CoordType(0.1);
          for(Index iter(0); iter < 1000; ++iter)
          {
            diff = CoordType(0);
            int iterations(0);
            int grad_evals(0);
            // Save old coordinates for explicit terms
            this->_coords_old.clone(this->_coords);

            ALGLIBWrapper<RumpfSmoother<CoordType_, MemType_, TrafoType_, FunctionalType_, H_EvalType_>>::minimise_functional_cg(grad_evals, iterations, termination_type,*this);

            total_grad_evals += grad_evals;
            total_iterations += iterations;

            for(Index i(0); i < this->_mesh.get_num_entities(0); ++i)
            {
              for(int d(0); d < MeshType::world_dim; ++d)
                diff+= Math::sqr(this->_coords_old(i)(d) - this->_coords(i)(d));
            }
            diff = Math::sqrt(diff/CoordType(this->_mesh.get_num_entities(0)));

            // Important: Copy back the coordinates the mesh optimiser changed to the original mesh.
            this->_coords_old.template scale(this->_coords_old,CoordType(1)-relaxation_parameter);
            this->_coords.template axpy(this->_coords, this->_coords_old, relaxation_parameter);

            this->set_coords();
            if(diff < CoordType(tol))
            {
              std::cout << iter << " fixed point iterations, " << total_iterations << " mincg iterations, " << total_grad_evals << " grad evals, last terminationtype was " << termination_type << std::endl;
              return;
            }
          }
        }

        /**
         * \brief Prepares the functional for evaluation.
         *
         * Needs to be called whenever any data like the mesh, the levelset function etc. changed.
         *
         **/
        virtual void prepare()
        {
        }
    }; // class RumpfSmootherSplit

  } // namespace Meshopt
} // namespace FEAST

#endif // KERNEL_MESHOPT_RUMPF_SMOOTHER_SPLIT_HPP
