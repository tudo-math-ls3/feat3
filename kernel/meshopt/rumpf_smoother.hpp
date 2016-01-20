#pragma once
#ifndef KERNEL_MESHOPT_RUMPF_SMOOTHER_HPP
#define KERNEL_MESHOPT_RUMPF_SMOOTHER_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/assembly/slip_filter_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/lafem/filter_chain.hpp>
#include <kernel/lafem/slip_filter.hpp>
#include <kernel/lafem/unit_filter_blocked.hpp>
#include <kernel/meshopt/h_evaluator.hpp>
#include <kernel/meshopt/mesh_smoother.hpp>
// ALGLIB includes
FEAST_DISABLE_WARNINGS
#include <optimization.h>
FEAST_RESTORE_WARNINGS

namespace FEAST
{
  namespace Meshopt
  {

    template<typename RumpfSmootherType_>
    struct ALGLIBWrapper;

    /**
     * \brief Baseclass for a family of variational mesh optimisation algorithms.
     *
     * \tparam TrafoType_
     * Type of the underlying transformation.
     *
     * \tparam FunctionalType_
     * Functional used for defining mesh quality. \see RumpfFunctional
     *
     * \tparam H_EvalType
     * Local meshsize evaluator. \see H_Evaluator
     *
     * Mesh optimisation algorithms derived from Martin Rumpf's paper \cite Rum96.
     *
     * \note The evaluation of the nonlinear functional requires operations that are essentially similar to an
     * assembly of an operator into a matrix. This is implemented for Mem::Main only. This in turn means that this
     * family of mesh optimisation algorithms is implemented for Mem::Main only.
     *
     * Assume we have a regular, conforming mesh \f$ \mathcal{T} \f$ and each cell \f$ K \in \mathcal{T} \f$ can
     * be expressed as the image of an (optimal) reference cell \f$ \hat{K} \f$ such that
     * \f[
     *   \forall K \in \mathcal{T}: \exists R_K : \hat{K} \to K: \forall \hat{x} \in \hat{K}: \det \nabla
     *   R_K(\hat{x}) \neq 0,
     * \f]
     * which ensures that the mapping is nonsingular and that the orientation of the reference cell is preserved. We
     * are now looking for a deformation \f$ \Phi: \mathcal{T} \to \Phi(\mathcal{T})\f$ such that
     * \f$ \Phi(\mathcal{T} \f$ is optimal in the sense that it minimises a functional of the form
     *
     * \f[
     *   \mathcal{F}(\Phi) = \int_\Omega \mathcal{L}(\Phi,x) dx = \sum_{K \in \mathcal{T}} \mu_K \int_K
     *   \mathcal{L}_K(\Phi,x) dx.
     * \f]
     *
     * Under the assumptions of frame indifference and translation invariance it can be shown that
     * \f[
     *   \forall K \in \mathcal{T}: \exists L_K \in SO_d \times K: \mathcal{L}_K: \mathcal{L}(\Phi,\cdot) =
     *   L_K(\nabla \Phi, \cdot) = L_K(\nabla R_K, \cdot)
     * \f]
     * and that the local functional is of the form
     * \f[
     *    F(\nabla R_K(\Phi))  := \int_K L(\nabla R_K (\Phi)(x)) dx = \mu_K L( \| \nabla R_T(\Phi) \|_F^2,
     *    \| \mathrm{cof} \nabla R_T(\Phi) \|_F^2, \det \nabla R_T(\Phi) )
     * \f]
     *
     * In the code, \f$ F \f$ is called the RumpfFunctional.
     *
     *
     * \author Jordi Paul
     *
     */
    template
    <
      typename TrafoType_,
      typename FunctionalType_,
      typename H_EvalType_ = H_Evaluator<TrafoType_, typename TrafoType_::MeshType::CoordType>
    >
    class RumpfSmootherBase:
      public MeshSmoother<typename TrafoType_::MeshType>
    {
      public :
        /// Type for the transformation
        typedef TrafoType_ TrafoType;
        /// The mesh the transformation is defined on
        typedef typename TrafoType::MeshType MeshType;
        /// The precision of the mesh coordinates
        typedef typename MeshType::CoordType CoordType;
        /// Type for the functional
        typedef FunctionalType_ FunctionalType;

        /// Only Mem::Main is supported atm
        typedef Mem::Main MemType;
        /// We always use the precision of the mesh
        typedef CoordType DataType;
        /// We always use Index for now
        typedef Index IndexType;

        /// Our base class
        typedef MeshSmoother<MeshType> BaseClass;

        /// ShapeType of said mesh
        typedef typename MeshType::ShapeType ShapeType;

        /// Vector type for element sizes etc.
        typedef LAFEM::DenseVector<MemType, CoordType, IndexType> ScalarVectorType;
        /// Vector type for coordinate vectors etc.
        typedef LAFEM::DenseVectorBlocked<MemType, CoordType, IndexType, MeshType::world_dim> VectorType;
        /// Filter for Dirichlet boundary conditions
        typedef LAFEM::UnitFilterBlocked<MemType, CoordType, IndexType, MeshType::world_dim> DirichletFilterType;
        /// Filter for slip boundary conditions
        typedef LAFEM::SlipFilter<MemType, CoordType, IndexType, MeshType::world_dim> SlipFilterType;
        /// Combined filter
        typedef LAFEM::FilterChain<SlipFilterType, DirichletFilterType> FilterType;

        /// Finite Element space for the transformation
        typedef typename Intern::TrafoFE<TrafoType>::Space TrafoSpace;

        /// Since the functional contains a ShapeType, these have to be the same
        static_assert(std::is_same<ShapeType, typename FunctionalType::ShapeType>::value,
        "ShapeTypes of the transformation / functional have to agree" );

        /// The transformation defining the physical mesh
        TrafoType _trafo;
        /// The FE space for the transformation, needed for filtering
        TrafoSpace _trafo_space;
        /// The filter enforcing boundary conditions
        FilterType _filter;

        /// Global gradient of the functional
        VectorType _grad;
        /// The functional for determining mesh quality
        FunctionalType& _functional;

        /// Assembler for Dirichlet boundary conditions
        Assembly::UnitFilterAssembler<MeshType> _dirichlet_asm;
        /// Assembler for slip boundary conditions
        Assembly::SlipFilterAssembler<MeshType> _slip_asm;

        // This is public for debugging purposes
      public:
        /// Weights for the local contributions to the global functional value.
        ScalarVectorType _mu;
        /// Weights for local mesh size
        // In the optimal case, every cell in the mesh has the size lambda(cell)
        ScalarVectorType _lambda;
        /// Size parameters for the local reference element.
        VectorType _h;
        /// List of boundary identifiers to enforce Dirichlet boundary conditions at
        std::deque<String> _dirichlet_list;
        /// List of boundary identifiers to enforce slip boundary conditions at
        std::deque<String> _slip_list;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] rmn_
         * The RootMeshNode representing the tree of root mesh, all of its MeshParts and Charts
         *
         * \param[in] dirichlet_list_
         * List of boundary identifiers for enforcing Dirichlet boundary conditions, can be empty
         *
         * \param[in] slip_list_
         * List of boundary identifiers for enforcing slip boundary conditions, can be empty
         *
         * \param[in] functional_
         * Reference to the functional used
         *
         */
        explicit RumpfSmootherBase(Geometry::RootMeshNode<MeshType>* rmn_,
        std::deque<String>& dirichlet_list_, std::deque<String>& slip_list_,
         FunctionalType_& functional_)
          : BaseClass(rmn_),
          _trafo(*(rmn_->get_mesh())),
          _trafo_space(_trafo),
          _filter(),
          _grad(rmn_->get_mesh()->get_num_entities(0)),
          _functional(functional_),
          _dirichlet_asm(),
          _slip_asm(*(rmn_->get_mesh())),
          _mu(rmn_->get_mesh()->get_num_entities(ShapeType::dimension)),
          _lambda(rmn_->get_mesh()->get_num_entities(ShapeType::dimension)),
          _h(rmn_->get_mesh()->get_num_entities(ShapeType::dimension)),
          _dirichlet_list(dirichlet_list_),
          _slip_list(slip_list_)
          {
            // Add all specified mesh parts to the Dirichlet filter assember
            for(auto& it : this->_dirichlet_list)
            {
              auto* mpp = this->_mesh_node->find_mesh_part(it);
              if(mpp != nullptr)
                _dirichlet_asm.add_mesh_part(*mpp);
            }

            // Assemble the homogeneous filter
            _dirichlet_asm.assemble(_filter.template at<1>(), _trafo_space);

            // Add all specified mesh parts to the slip filter assember
            for(auto& it : this->_slip_list)
            {
              auto* mpp = this->_mesh_node->find_mesh_part(it);
              if(mpp != nullptr)
                _slip_asm.add_mesh_part(*mpp);
            }

          }

        /// \brief Destructor
        virtual ~RumpfSmootherBase()
        {
        };

        /**
         * \brief The class name
         *
         * \returns String with the class name
         */
        static String name()
        {
          return "RumpfSmootherBase<"+MeshType::name()+">";
        }

        /**
         * \brief Prints some characteristics of the RumpfSmoother object
         **/
        virtual void print()
        {
          if(this->_dirichlet_list.size() > 0)
          {
            std::cout << "Dirichlet boundaries:";
            for(auto& it : this->_dirichlet_list)
              std::cout << " " << it;
            std::cout << std::endl;
          }

          if(this->_slip_list.size() > 0)
          {
            std::cout << "Slip boundaries:";
            for(auto& it : this->_slip_list)
              std::cout << " " << it;
            std::cout << std::endl;
          }

          _functional.print();
        }

        /**
         * \brief Performs one-time initialisations
         *
         * These are not done in the constructor as compute_lambda() etc. could be overwritten in a derived class,
         * potentially using uninitialised members.
         *
         */
        virtual void init() override
        {
          // Write any potential changes to the mesh
          this->set_coords();
          // Assemble the homogeneous filter
          this->_slip_asm.assemble(this->_filter.template at<0>(), this->_trafo_space);

          // Compute desired element size distribution
          this->compute_lambda();
          // Compute target scales
          this->compute_h();
          // Compute element weights
          this->compute_mu();
        }

        /**
         * \brief Prepares the functional for evaluation.
         *
         * Needs to be called whenever any data like the mesh, the levelset function etc. changed.
         *
         **/
        virtual void prepare() override
        {
          this->set_coords();

          // Adapt all slip boundaries
          for(auto& it : _slip_list)
            this->_mesh_node->adapt_by_name(it);

          // The slip filter contains the outer unit normal, so reassemble it
          _slip_asm.assemble(_filter.template at<0>(), _trafo_space);

          this->get_coords();
        }

        /**
         * \brief Computes a quality indicator concerning the cell sizes
         *
         * In a truly optimal mesh (consisting ONLY of Rumpf reference cells of the right size), every cell's volume is
         * exactly lambda(cell). This is especially the goal for r-adaptivity.
         * So in an optimal mesh,
         * \f[
         *   \forall K \in \mathcal{T}_h: \frac{|K|}{\lambda(K)} = 1,
         * \f]
         * so we compute the Euclidean norm of the vector \f$(v)_i = \frac{1}{N}(1 -  \frac{|K_i|}{\lambda(K_i)} \f$.
         * This is scaled by the number of cells so it is independant of the refinement level. Not sure if the
         * scaling part is sensible, though.
         *
         * \returns The relative cell size quality indicator.
         *
         **/
        CoordType cell_size_quality()
        {
          typename LAFEM::DenseVector<Mem::Main, CoordType, Index> tmp(
            this->get_mesh()->get_num_entities(ShapeType::dimension));

          CoordType my_vol(0);

          for(Index cell(0); cell < this->get_mesh()->get_num_entities(ShapeType::dimension); ++cell)
          {
            my_vol = this->_trafo.template compute_vol<ShapeType, CoordType>(cell);
            tmp(cell, Math::abs(CoordType(1) - my_vol/this->_lambda(cell)));
          }

          return tmp.norm2()/Math::sqrt(CoordType(this->get_mesh()->get_num_entities(ShapeType::dimension)));
        }

        /**
         * \brief Computes the functional value on the current mesh.
         *
         * \returns
         * The functional value
         *
         **/
        virtual CoordType compute_functional() = 0;

        /**
         * \brief Computes the gradient of the functional with regard to the nodal coordinates.
         *
         * Usually, prepare() should be called before calling this.
         *
         */
        virtual void compute_gradient() = 0;

        /**
         * \brief Computes the volume of the optimal reference for each cell and saves it to _h.
         *
         **/
        virtual void compute_h()
        {
          H_EvalType_::compute_h(_h, this->_coords, _lambda, this->_trafo);
        }

        /// \brief Computes the weights _lambda.
        virtual void compute_lambda()
        {
          compute_lambda_current();
        }

        /// \brief Computes the uniformly distributed weights _lambda.
        virtual void compute_lambda_current()
        {
          Index ncells(this->get_mesh()->get_num_entities(ShapeType::dimension));

          // As this uses the transformation, it is always carried out in Mem::Main
          typename LAFEM::DenseVector<Mem::Main, CoordType, Index> tmp(
            this->get_mesh()->get_num_entities(ShapeType::dimension));

          CoordType sum_lambda(0);
          for(Index cell(0); cell < ncells; ++cell)
          {
            tmp(cell, this->_trafo.template compute_vol<ShapeType, CoordType>(cell));
            sum_lambda+=tmp(cell);
          }

          // Scale so that sum(lambda) = 1
          tmp.scale(tmp, CoordType(1)/sum_lambda);
          _lambda.convert(tmp);

        }

        /// \brief Computes the uniformly distributed weights _lambda.
        virtual void compute_lambda_uniform()
        {
          Index ncells(this->get_mesh()->get_num_entities(ShapeType::dimension));

          _lambda.format(CoordType(1)/CoordType(ncells));
        }

        /// \brief Computes the weights mu
        virtual void compute_mu()
        {
          Index ncells(this->get_mesh()->get_num_entities(ShapeType::dimension));
          _mu.format(CoordType(1)/CoordType(ncells));
        }

    }; // class RumpfSmootherBase

    /// \copydoc Meshopt::RumpfSmootherBase
    template
    <
      typename TrafoType_,
      typename FunctionalType_,
      typename H_EvalType_ = H_Evaluator<TrafoType_, typename TrafoType_::MeshType::CoordType>
    >
    class RumpfSmoother:
      public RumpfSmootherBase<TrafoType_, FunctionalType_, H_EvalType_>
    {
      public :
        /// Our base class
        typedef RumpfSmootherBase<TrafoType_, FunctionalType_, H_EvalType_> BaseClass;

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
        /// Filter for Dirichlet boundary conditions
        typedef LAFEM::UnitFilterBlocked<MemType, CoordType, IndexType, MeshType::world_dim> DirichletFilterType;
        /// Filter for slip boundary conditions
        typedef LAFEM::SlipFilter<MemType, CoordType, IndexType, MeshType::world_dim> SlipFilterType;
        /// Combined filter
        typedef LAFEM::FilterChain<SlipFilterType, DirichletFilterType> FilterType;
        /// Finite Element space for the transformation
        typedef typename Intern::TrafoFE<TrafoType>::Space TrafoSpace;

        /// Since the functional contains a ShapeType, these have to be the same
        static_assert( std::is_same<ShapeType, typename FunctionalType::ShapeType>::value,
        "ShapeTypes of the transformation / functional have to agree" );

      public:
        /**
         * \copydoc RumpfSmootherBase()
         *
         */
        explicit RumpfSmoother(Geometry::RootMeshNode<MeshType>* rmn_,
        std::deque<String>& dirichlet_list_, std::deque<String>& slip_list_,
        FunctionalType_& functional_)
          : BaseClass(rmn_, dirichlet_list_, slip_list_, functional_)
          {
          }

        /// \brief Destructor
        virtual ~RumpfSmoother()
        {
        };

        /**
         * \brief The class name
         *
         * \returns String with the class name
         */
        static String name()
        {
          return "RumpfSmoother<"+MeshType::name()+">";
        }

        /**
         * \brief Prints some characteristics of the RumpfSmoother object
         */
        virtual void print() override
        {
          std::cout << name() << std::endl;
          BaseClass::print();
        }

        /// \copydoc RumpfSmootherBase::compute_functional()
        virtual CoordType compute_functional() override
        {
          CoordType fval(0);
          // Total number of cells in the mesh
          Index ncells(this->get_mesh()->get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->get_mesh()->template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector<CoordType, MeshType::world_dim> h;

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            h = this->_h(cell);
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; j++)
              x[j] = this->_coords(idx(cell,Index(j)));

            fval += this->_mu(cell) * this->_functional.compute_local_functional(x,h);
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
        virtual CoordType compute_functional(CoordType* func_norm, CoordType* func_det, CoordType* func_rec_det)
        {
          CoordType fval(0);
          // Total number of cells in the mesh
          Index ncells(this->get_mesh()->get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->get_mesh()->template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector<CoordType, MeshType::world_dim> h;

          CoordType norm_A(0), det_A(0), rec_det_A(0);

          CoordType func_norm_tot(0);
          CoordType func_det_tot(0);
          CoordType func_rec_det_tot(0);
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

          std::cout << "fval = " << stringify_fp_sci(fval) << " func_norm = " << stringify_fp_sci(func_norm_tot)
          << ", func_det = " << stringify_fp_sci(func_det_tot) << ", func_rec_det = " << stringify_fp_sci(func_rec_det_tot) << std::endl;

          return fval;
        } // compute_functional

        /// \copydoc BaseClass::compute_gradient()
        virtual void compute_gradient() override
        {
          // Total number of cells in the mesh
          Index ncells(this->get_mesh()->get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->get_mesh()->template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAST::Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;
          // Local cell dimensions for passing to other routines
          FEAST::Tiny::Vector<CoordType, MeshType::world_dim> h;
          // This will hold the local gradient for one element for passing to other routines
          FEAST::Tiny::Matrix<CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> grad_loc;

          // Clear gradient vector
          this->_grad.format();

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            h = this->_h(cell);
            // Get local coordinates
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
              x[j] = this->_coords(idx(cell,Index(j)));

            this->_functional.compute_local_grad(x, h, grad_loc);

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
        virtual void optimise() override
        {
          int total_grad_evals(0);
          int total_iterations(0);
          int termination_type(0);

          ALGLIBWrapper<RumpfSmoother<TrafoType_, FunctionalType_, H_EvalType_>>::minimise_functional_cg(total_grad_evals, total_iterations, termination_type,*this);
          // Important: Copy back the coordinates the mesh optimiser changed to the original mesh.
          this->prepare();
          std::cout << total_iterations << " mincg iterations, " << total_grad_evals << " grad evals, terminationtype was " << termination_type << std::endl;
        }

    }; // class RumpfSmoother

    /**
     * \brief Wrapper struct for calling and passing data to ALGLIB.
     *
     * \tparam RumpfSmootherType_
     * The type of the RumpfSmoother used, needed for safe downcasting.
     *
     * \author Jordi Paul
     *
     */
    template<typename RumpfSmootherType_>
    struct ALGLIBWrapper
    {
      /// Datatype from the RumpfSmoother
      typedef typename RumpfSmootherType_::CoordType CoordType;

      /**
       * \brief Minimises the functional using ALGLIB's nonlinear CG
       *
       * \param[out] grad_eval_count
       * The number of gradient evaluations is returned here
       *
       * \param[out] iteration_count
       * The number of optimiser iterations is returned here
       *
       * \param[out] termination_type
       * The termination type of the optimiser is returned here
       *
       * \param[in] my_smoother
       * Mesh smoother providing data and routines for evaluating the functional and its gradient.
       *
       */
      static void minimise_functional_cg(int& grad_eval_count, int& iteration_count, int& termination_type,
      RumpfSmootherType_& my_smoother)
      {

        const int world_dim(RumpfSmootherType_::MeshType::world_dim);
        const Index num_vert(my_smoother.get_mesh()->get_num_entities(0));

        // Array for the coordinates for passing to ALGLIB's optimiser
        alglib::real_1d_array x;
        // Warning: ae_int_t/Index conversion
        x.setlength(alglib::ae_int_t(world_dim*num_vert));
        for(Index j(0); j < num_vert; ++j)
        {
          for(int d(0); d < world_dim; ++d)
          {
            // Warning: Index/ae_int_t conversion
            x[alglib::ae_int_t(world_dim)*alglib::ae_int_t(j) + alglib::ae_int_t(d)]
              = double(my_smoother._coords(j)(d));
          }
        }

        double epsg = 1.e-10;
        double epsf = 0.;
        double epsx = 0.;
        alglib::ae_int_t maxits = 10000; // 1000
        alglib::mincgstate state;
        alglib::mincgreport rep;

        alglib::mincgcreate(x, state);
        alglib::mincgsetcgtype(state, -1);
        alglib::mincgsetcond(state, epsg, epsf, epsx, maxits);
        // mincgsuggeststep(state, 1.e-4);
        alglib::mincgoptimize(state, ALGLIBWrapper::functional_grad, nullptr, &my_smoother);
        alglib::mincgresults(state, x, rep);

        // This is the variant that uses automatic differentiation, but it is too slow in general
        //alglib::mincgcreatef(x, 1e-6, state);
        //alglib::mincgsetcgtype(state, -1);
        //alglib::mincgsetcond(state, epsg, epsf, epsx, maxits);
        //// mincgsuggeststep(state, 1.e-4);
        //alglib::mincgoptimize(state, ALGLIBWrapper::functional, nullptr, &my_smoother);
        //alglib::mincgresults(state, x, rep);

        std::cout << "mincg: terminationtype " << rep.terminationtype << ", " << rep.iterationscount << " its, " << rep.nfev << " grad evals" << std::endl;

        // Warning: ae_int_t to int conversions
        iteration_count = int(rep.iterationscount);
        grad_eval_count = int(rep.nfev);
        termination_type = int(rep.terminationtype);
        if(rep.terminationtype == -8)
          throw InternalError(__func__, __FILE__, __LINE__, "Optimizer stopped with status -8 (internal ALGLIB error)");

      }

      ///**
      // * \brief Computes the functional value
      // *
      // * The actual work is done by the mesh smoother's routines. The mesh smoother is passed as ptr and then casted
      // *
      // * \param[in] x
      // * Vertex coordinates for the mesh in its current state in the optimisation process.
      // *
      // * \param[out] func
      // * Functional value.
      // *
      // * \param[out] grad
      // * Gradient of the functional.
      // *
      // * \param[in] ptr
      // * The mesh smoother.
      // *
      // */
      //static void functional(const alglib::real_1d_array& x, double& func, void* ptr)
      //{
      //  // Evil downcast, but we know what we are doing, right?
      //  RumpfSmootherType_* my_smoother = reinterpret_cast<RumpfSmootherType_*> (ptr);

      //  const int world_dim(my_smoother->get_world_dim());
      //  const Index num_vert(my_smoother->get_num_vert());

      //  // Copy back the vertex coordinates, needed for computing the gradient on the modified mesh
      //  for(int d(0); d < world_dim; ++d)
      //  {
      //    for(Index j(0); j < num_vert; ++j)
      //    {
      //      // Skip the bounday vertices to preserve their value
      //      if(my_smoother->_bdry_id[j] >=0)
      //        // Warning: Index/ae_int_t conversion
      //        my_smoother->_coords[d](j,CoordType(x[alglib::ae_int_t(d*num_vert + j)]));
      //    }
      //  }

      //  // Let the functional initialise stuff if it needs to, like evaluating the levelset function on the current
      //  // mesh
      //  my_smoother->prepare();
      //  // Compute functional value
      //  func = double(my_smoother->compute_functional());

      //}

      /**
       * \brief Computes the functional value and its gradient
       *
       * The actual work is done by the mesh smoother's routines. The mesh smoother is passed as ptr and then casted
       *
       * \param[in] x
       * Vertex coordinates for the mesh in its current state in the optimisation process.
       *
       * \param[out] func
       * Functional value.
       *
       * \param[out] grad
       * Gradient of the functional.
       *
       * \param[in] ptr
       * The mesh smoother.
       *
       */
      static void functional_grad(const alglib::real_1d_array& x, double& func, alglib::real_1d_array& grad, void* ptr)
      {
        // Evil downcast, but we know what we are doing, right?
        RumpfSmootherType_* my_smoother = reinterpret_cast<RumpfSmootherType_*> (ptr);

        const int world_dim(RumpfSmootherType_::MeshType::world_dim);
        const Index num_vert(my_smoother->get_mesh()->get_num_entities(0));

        // Copy back the vertex coordinates, needed for computing the gradient on the modified mesh
        for(Index j(0); j < num_vert; ++j)
        {
          Tiny::Vector<CoordType, world_dim, world_dim> tmp;
          // Warning: Index/ae_int_t conversion
          for(int d(0); d < world_dim; ++d)
            tmp(d) = CoordType(x[alglib::ae_int_t(world_dim)*alglib::ae_int_t(j) + alglib::ae_int_t(d)]);

          my_smoother->_coords(j,tmp);
        }

        // Let the functional initialise stuff if it needs to, like evaluating the levelset function on the current
        // mesh
        my_smoother->prepare();
        // Compute functional value
        func = double(my_smoother->compute_functional());
        // Compute functional gradient
        my_smoother->compute_gradient();

        // Copy to array for the gradient for passing to ALGLIB's optimiser
        //CoordType norm_grad(0);
        for(Index j(0); j < world_dim*num_vert; ++j)
        {
          grad[alglib::ae_int_t(j)] = my_smoother->_grad.template elements<LAFEM::Perspective::pod>()[j];
        }
      }
    }; // struct ALGLIBWrapper

  } // namespace Meshopt
} // namespace FEAST
#endif // KERNEL_MESHOPT_RUMPF_SMOOTHER_HPP
