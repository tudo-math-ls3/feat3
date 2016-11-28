#pragma once
#ifndef KERNEL_MESHOPT_HYPERELASTICITY_FUNCTIONAL_HPP
#define KERNEL_MESHOPT_HYPERELASTICITY_FUNCTIONAL_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/assembly/slip_filter_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/lafem/filter_chain.hpp>
#include <kernel/lafem/filter_sequence.hpp>
#include <kernel/lafem/slip_filter.hpp>
#include <kernel/lafem/unit_filter_blocked.hpp>
#include <kernel/meshopt/mesh_concentration_function.hpp>
#include <kernel/meshopt/mesh_quality_functional.hpp>
#include <kernel/meshopt/rumpf_trafo.hpp>
#include <kernel/util/dist.hpp>

#include <map>

// This is for the explicit template instantiation below
#ifdef FEAT_EICKT
#include <kernel/meshopt/rumpf_functional.hpp>
#include <kernel/meshopt/rumpf_functionals/p1.hpp>
#include <kernel/meshopt/rumpf_functionals/q1.hpp>
#endif


namespace FEAT
{
  namespace Meshopt
  {

    /**
     * \brief Enum class for different types of scale computations
     *
     * - once_*** means the scales are only computed once and remain the same over all calls to the functional
     * - current_*** means the scales are recomputed with every call to init() (i.e. if they are to be computed for
     *   every timestep in a transient simulation)
     * - iter_*** means the scales are recomputed with every call to prepare(), i.e. in every nonlinear solver
     *   iteration. If the scales explicitly depend on the solution (i.e. the distance of the mesh vertices to some
     *   surface), this is needed.
     */
    enum class ScaleComputation
    {
      undefined = 0,
      once_uniform,
      once_cellsize,
      once_concentration,
      current_uniform,
      current_cellsize,
      current_concentration,
      iter_concentration,
      };

    /// \cond internal
    /**
     * \brief Streaming operator for ScaleComputations
     */
    inline std::ostream& operator<<(std::ostream& os, ScaleComputation scale_computation)
    {
      switch(scale_computation)
      {
        case ScaleComputation::undefined:
          return os << "undefined";
        case ScaleComputation::once_uniform:
          return os << "once_uniform";
        case ScaleComputation::once_cellsize:
          return os << "once_cellsize";
        case ScaleComputation::once_concentration:
          return os << "once_concentration";
        case ScaleComputation::current_uniform:
          return os << "current_uniform";
        case ScaleComputation::current_cellsize:
          return os << "current_cellsize";
        case ScaleComputation::current_concentration:
          return os << "current_concentration";
        case ScaleComputation::iter_concentration:
          return os << "iter_concentration";
        default:
          return os << "-unknown-";
      }
    }

    inline void operator<<(ScaleComputation& scale_computation, const String& sc_name)
    {
      if(sc_name == "undefined")
        scale_computation = ScaleComputation::undefined;
      else if(sc_name == "once_uniform")
        scale_computation = ScaleComputation::once_uniform;
      else if(sc_name == "once_cellsize")
        scale_computation = ScaleComputation::once_cellsize;
      else if(sc_name == "once_concentration")
        scale_computation = ScaleComputation::once_concentration;
      else if(sc_name == "current_uniform")
        scale_computation = ScaleComputation::current_uniform;
      else if(sc_name == "current_cellsize")
        scale_computation = ScaleComputation::current_cellsize;
      else if(sc_name == "current_concentration")
        scale_computation = ScaleComputation::current_concentration;
      else if(sc_name == "iter_concentration")
        scale_computation = ScaleComputation::iter_concentration;
      else
        throw InternalError(__func__, __FILE__, __LINE__, "Unknown ScaleComputation identifier string "
            +sc_name);
    }
    /// \endcond

    /**
     * \brief Baseclass for a family of variational mesh optimisation algorithms.
     *
     * \tparam Mem_
     * Memory architecture for solving.
     *
     * \tparam DT_
     * Floating point type for solving.
     *
     * \tparam IT_
     * Index type for solver vectors etc.
     *
     * \tparam Trafo_
     * Type of the underlying transformation.
     *
     * \tparam FunctionalType_
     * Functional used for defining mesh quality. \see RumpfFunctional
     *
     * \tparam RefCellTrafo_
     * Basically choses the reference cell according to which to optimise the mesh quality for.
     *
     * Mesh optimisation algorithms derived from Martin Rumpf's paper \cite Rum96.
     *
     * \note The evaluation of the nonlinear functional requires operations that are essentially similar to an
     * assembly of an operator into a matrix. This is implemented for Mem::Main only. This in turn means that this
     * family of mesh optimisation algorithms is implemented for Mem::Main only.
     *
     * See the specialisation for Mem::Main.
     *
     */
    template
    <
      typename Mem_,
      typename DT_,
      typename IT_,
      typename Trafo_,
      typename FunctionalType_,
      typename RefCellTrafo_ = RumpfTrafo<Trafo_, typename Trafo_::MeshType::CoordType>
    >
    class HyperelasticityFunctional;
    /**
     * \brief Baseclass for a family of variational mesh optimisation algorithms.
     *
     * \tparam DT_
     * Floating point type for solving.
     *
     * \tparam IT_
     * Index type for solver vectors etc.
     *
     * \tparam Trafo_
     * Type of the underlying transformation.
     *
     * \tparam FunctionalType_
     * Functional used for defining mesh quality. \see RumpfFunctional
     *
     * \tparam RefCellTrafo_
     * Basically choses the reference cell according to which to optimise the mesh quality for.
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
     * \author Jordi Paul
     *
     */
    template
    <
      typename DT_,
      typename IT_,
      typename Trafo_,
      typename FunctionalType_,
      typename RefCellTrafo_
    >
    class HyperelasticityFunctional<Mem::Main, DT_, IT_, Trafo_, FunctionalType_, RefCellTrafo_>:
    public MeshQualityFunctional<typename Trafo_::MeshType>
    {
      public :
        /// Type for the transformation
        typedef Trafo_ TrafoType;
        /// The mesh the transformation is defined on
        typedef typename TrafoType::MeshType MeshType;
        /// The precision of the mesh coordinates
        typedef typename MeshType::CoordType CoordType;
        /// Type for the functional
        typedef FunctionalType_ FunctionalType;

        /// Type of the reference cell trafo for the mesh quality
        typedef RefCellTrafo_ RefCellTrafo;

        /// Only Mem::Main is supported atm
        typedef Mem::Main MemType;
        /// We always use the precision of the mesh
        typedef CoordType DataType;
        /// We always use Index for now
        typedef Index IndexType;

        /// Our base class
        typedef MeshQualityFunctional<MeshType> BaseClass;

        /// Block height of the operator's gradient
        static constexpr int BlockHeight = MeshType::world_dim;
        /// Block Weight of the operator's gradient
        static constexpr int BlockWidth = MeshType::world_dim;

        /// ShapeType of said mesh
        typedef typename MeshType::ShapeType ShapeType;

        /// Vector type for element sizes etc.
        typedef LAFEM::DenseVector<MemType, CoordType, IndexType> ScalarVectorType;
        /// Vector type for coordinate vectors etc.
        typedef LAFEM::DenseVectorBlocked<MemType, CoordType, IndexType, MeshType::world_dim> VectorType;
        /// Filter for Dirichlet boundary conditions
        typedef LAFEM::UnitFilterBlocked<MemType, DT_, IT_, MeshType::world_dim> DirichletFilterType;
        /// Sequence of Dirichlet filters for several different boundary parts
        typedef LAFEM::FilterSequence<DirichletFilterType> DirichletFilterSequence;
        /// Filter for slip boundary conditions
        typedef LAFEM::SlipFilter<MemType, DT_, IT_, MeshType::world_dim> SlipFilterType;
        /// Sequence of Slip filters for several different boundary parts
        typedef LAFEM::FilterSequence<SlipFilterType> SlipFilterSequence;
        /// Combined filter
        typedef LAFEM::FilterChain<SlipFilterSequence, DirichletFilterSequence> FilterType;

        /// Finite Element space for the transformation
        typedef typename Intern::TrafoFE<TrafoType>::Space TrafoSpace;

        /// Output vector type of the operator
        typedef LAFEM::DenseVectorBlocked<MemType, DT_, IT_, MeshType::world_dim> VectorTypeL;
        /// Input vector type of the operator
        typedef LAFEM::DenseVectorBlocked<MemType, DT_, IT_, MeshType::world_dim> VectorTypeR;
        /// Type of the gradient vector
        typedef VectorTypeR GradientType;
        /// Type for exchanging information between state variable and mesh
        typedef typename BaseClass::CoordsBufferType CoordsBufferType;

        /// Since the functional contains a ShapeType, these have to be the same
        static_assert(std::is_same<ShapeType, typename FunctionalType::ShapeType>::value,
        "ShapeTypes of the transformation / functional have to agree" );

        /// The transformation defining the physical mesh
        TrafoType& _trafo;
        /// The FE space for the transformation, needed for filtering
        TrafoSpace& _trafo_space;
        /// Assembler for Dirichlet boundary conditions
        std::map<String, std::shared_ptr<Assembly::UnitFilterAssembler<MeshType>>>& _dirichlet_asm;
        /// Assembler for slip boundary conditions
        std::map<String, std::shared_ptr<Assembly::SlipFilterAssembler<MeshType>>>& _slip_asm;
        /// The functional for determining mesh quality
        std::shared_ptr<FunctionalType> _functional;
        /// The mesh concentration function (if any)
        std::shared_ptr<MeshConcentrationFunctionBase<Trafo_, RefCellTrafo_>> _mesh_conc;
        /// These are the scalars that need to be synchronised in init() or prepare()
        std::map<String, DataType*> sync_scalars;
        /// These are the vectors that need to be synchronised (type-0 to type-1)
        std::set<VectorTypeR*> sync_vecs;

      protected:
        /// Weights for the local contributions to the global functional value.
        ScalarVectorType _mu;
        /// Weights for local mesh size
        // In the optimal case, every cell in the mesh has the size lambda(cell)
        ScalarVectorType _lambda;
        /// Size parameters for the local reference element.
        ScalarVectorType _h;

      private:
        /// This is the number of DoFs in the trial space (= number of mesh vertices)
        const Index _columns;
        /// This is the number of DoFs in the test space (= number of mesh vertices)
        const Index _rows;
        /// How to compute the optimal scales
        ScaleComputation _scale_computation;

        /// Factor for the alignment penalty term
        DataType _penalty_param;
        /// Last computed contraint (violation)
        DataType _alignment_constraint;

        /// This is the sum of the determinant of the transformation of the reference cell to the mesh cell over all
        /// cells
        DataType _sum_det;
        /// The sum of all lambda (= optimal cell size) over all cells before normalisation
        DataType _sum_lambda;
        /// The sum of all mu (= cell weights) over all cells before normalisation
        DataType _sum_mu;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] rmn_
         * The RootMeshNode representing the tree of root mesh, all of its MeshParts and Charts
         *
         * \param[in] trafo_space_
         * A reference to the Finite Element space used for the transformation
         *
         * \param[in] dirichlet_asm_
         * A map of Strings to UnitFilterAssemblers for all Dirichlet boundaries
         *
         * \param[in] slip_asm_
         * A map of Strings to SlipFilterAssemblers for all slip boundaries
         *
         * \param[in] functional_
         * The (cell-local) functional used
         *
         * This is the simple constructor that always sets the scale computation to once_uniform and the mesh
         * concentration function to nullptr.
         *
         */
        explicit HyperelasticityFunctional(
          Geometry::RootMeshNode<MeshType>* rmn_,
          TrafoSpace& trafo_space_,
          std::map<String, std::shared_ptr<Assembly::UnitFilterAssembler<MeshType>>>& dirichlet_asm_,
          std::map<String, std::shared_ptr<Assembly::SlipFilterAssembler<MeshType>>>& slip_asm_,
          std::shared_ptr<FunctionalType_> functional_)
          : BaseClass(rmn_),
          _trafo(trafo_space_.get_trafo()),
          _trafo_space(trafo_space_),
          _dirichlet_asm(dirichlet_asm_),
          _slip_asm(slip_asm_),
          _functional(functional_),
          _mesh_conc(nullptr),
          sync_scalars(),
          sync_vecs(),
          _mu(rmn_->get_mesh()->get_num_entities(ShapeType::dimension), DataType(1)),
          _lambda(rmn_->get_mesh()->get_num_entities(ShapeType::dimension)),
          _h(rmn_->get_mesh()->get_num_entities(ShapeType::dimension)),
          _columns(_trafo_space.get_num_dofs()),
          _rows(_trafo_space.get_num_dofs()),
          _scale_computation(ScaleComputation::once_uniform),
          _penalty_param(0),
          _alignment_constraint(0),
          _sum_det(0),
          _sum_lambda(0),
          _sum_mu(0)
          {

            XASSERTM(functional_ != nullptr, "Cell functional must not be nullptr");

            for(Index i(0); i < _mu.size(); ++i)
            {
              _sum_mu += _mu(i);
            }

            sync_scalars.emplace("_sum_mu",&_sum_mu);

            // Compute desired element size distribution
            _compute_scales_once();
          }

        /**
         * \brief Constructor
         *
         * \param[in] rmn_
         * The RootMeshNode representing the tree of root mesh, all of its MeshParts and Charts
         *
         * \param[in] trafo_space_
         * A reference to the Finite Element space used for the transformation
         *
         * \param[in] dirichlet_asm_
         * A map of Strings to UnitFilterAssemblers for all Dirichlet boundaries
         *
         * \param[in] slip_asm_
         * A map of Strings to SlipFilterAssemblers for all slip boundaries
         *
         * \param[in] functional_
         * The (cell-local) functional used
         *
         * \param[in] scale_computation_
         * The type of scale computation to use
         *
         * \param[in] mesh_conc_
         * The mesh concentration function to use. If no scale computation that uses this is set, this may be
         * nullptr.
         *
         * \param[in] penalty_param_
         * The penalty parameter, e.g. for surface alignment. Defaults to 0, which turns off this functionality.
         *
         */
        explicit HyperelasticityFunctional(
          Geometry::RootMeshNode<MeshType>* rmn_,
          TrafoSpace& trafo_space_,
          std::map<String, std::shared_ptr<Assembly::UnitFilterAssembler<MeshType>>>& dirichlet_asm_,
          std::map<String, std::shared_ptr<Assembly::SlipFilterAssembler<MeshType>>>& slip_asm_,
          std::shared_ptr<FunctionalType_> functional_,
          ScaleComputation scale_computation_,
          std::shared_ptr<MeshConcentrationFunctionBase<Trafo_, RefCellTrafo_>> mesh_conc_ = nullptr,
          DataType penalty_param_ = DataType(0))
          : BaseClass(rmn_),
          _trafo(trafo_space_.get_trafo()),
          _trafo_space(trafo_space_),
          _dirichlet_asm(dirichlet_asm_),
          _slip_asm(slip_asm_),
          _functional(functional_),
          _mesh_conc(nullptr),
          sync_scalars(),
          sync_vecs(),
          _mu(rmn_->get_mesh()->get_num_entities(ShapeType::dimension), DataType(1)),
          _lambda(rmn_->get_mesh()->get_num_entities(ShapeType::dimension)),
          _h(rmn_->get_mesh()->get_num_entities(ShapeType::dimension)),
          _columns(_trafo_space.get_num_dofs()),
          _rows(_trafo_space.get_num_dofs()),
          _scale_computation(scale_computation_),
          _penalty_param(penalty_param_),
          _alignment_constraint(0),
          _sum_det(0),
          _sum_lambda(0),
          _sum_mu(0)
          {

            XASSERTM(functional_ != nullptr, "Cell functional must not be nullptr!\n");
            XASSERTM(_penalty_param >= DataType(0), "penalty_param must be >= 0!\n");

            //_mu(17, DataType(128));
            //_mu(20, DataType(128));
            //_mu(25, DataType(128));
            //_mu(28, DataType(128));

            for(Index i(0); i < _mu.size(); ++i)
            {
              _sum_mu += _mu(i);
            }

            sync_scalars.emplace("_sum_mu",&_sum_mu);

            if(
              ( _scale_computation == ScaleComputation::once_concentration ||
                _scale_computation == ScaleComputation::current_concentration ||
                _scale_computation == ScaleComputation::iter_concentration ) &&
              mesh_conc_ == nullptr)
              {
                throw InternalError(__func__,__FILE__,__LINE__,
                "Scale computation set to "+stringify(_scale_computation)+", but no concentration function was given");
              }

            if(mesh_conc_ != nullptr)
            {
              _mesh_conc = mesh_conc_->create_empty_clone();
              _mesh_conc->set_mesh_node(rmn_);
              _mesh_conc->add_sync_vecs(sync_vecs);
            }

            // Perform one time scal computation
            _compute_scales_once();
          }

        /// Explicitly delete default constructor
        HyperelasticityFunctional() = delete;
        /// Explicitly delete copy constructor
        HyperelasticityFunctional(const HyperelasticityFunctional&) = delete;

        /// \brief Destructor
        virtual ~HyperelasticityFunctional()
        {
        }

        /**
         * \brief Checks if the functional is empty (= the null functional)
         *
         * \returns True if the number of DoFs is zero.)
         */
        bool empty() const
        {
          return (_trafo_space.get_num_dofs() == Index(0));
        }

        /**
         * \brief Creates an L-vector for the functional's gradient
         */
        VectorTypeL create_vector_l() const
        {
          return VectorTypeL(_trafo_space.get_num_dofs());
        }

        /**
         * \brief Creates an R-vector for the functional and its gradient
         */
        VectorTypeR create_vector_r() const
        {
          return VectorTypeR(_trafo_space.get_num_dofs());
        }

        /**
         * \brief Returns the number of columns
         *
         * \returns The number of columns.
         */
        const Index& columns() const
        {
          return _columns;
        }

        /**
         * \brief Returns the number of rows
         *
         * \returns The number of rows.
         */
        const Index& rows() const
        {
          return _rows;
        }

        /**
         * \brief The class name
         *
         * \returns String with the class name
         */
        static String name()
        {
          return "HyperelasticityFunctional<"+MeshType::name()+">";
        }

        /**
         * \brief Prints some characteristics of the HyperelasticityFunctional object
         */
        virtual void print()
        {
          Index pad_width(30);
          Dist::Comm comm_world(Dist::Comm::world());

          String msg(name()+":");
          comm_world.print(msg);

          msg = String("Scale computation").pad_back(pad_width, '.') + String(": ") + stringify(_scale_computation);
          comm_world.print(msg);

          _functional->print();
          if(_mesh_conc != nullptr)
            _mesh_conc->print();

        }

        /**
         * \brief Adds relevant quantities of this object to a VTK exporter
         *
         * \param[in, out] exporter
         * The exporter to add our data to.
         *
         * \note This does not compute and/or add the functional's gradient to the exporter because this requires
         * filtering and synchronisation if there is more than one patch. Compute the gradient at a point where
         * you know these things (e.g. in the application itself, or in Control::Meshopt::HyperelasticityControl).
         *
         */
        virtual void add_to_vtk_exporter(Geometry::ExportVTK<typename Trafo_::MeshType>& exporter) const
        {
          exporter.add_cell_scalar("h", this->_h.elements());
          exporter.add_cell_scalar("lambda", this->_lambda.elements());
          exporter.add_cell_scalar("mu", this->_mu.elements());

          DataType* fval_norm(new DataType[this->get_mesh()->get_num_entities(MeshType::shape_dim)]);
          DataType* fval_det(new DataType[this->get_mesh()->get_num_entities(MeshType::shape_dim)]);
          DataType* fval_rec_det(new DataType[this->get_mesh()->get_num_entities(MeshType::shape_dim)]);

          DataType fval(0);
          eval_fval_cellwise(fval, fval_norm, fval_det, fval_rec_det);
          exporter.add_cell_vector("fval", fval_norm, fval_det, fval_rec_det);

          if(_mesh_conc != nullptr)
          {
            exporter.add_vertex_scalar("dist", _mesh_conc->get_dist().elements());
            exporter.add_vertex_vector("grad_dist", _mesh_conc->get_grad_dist());

            if(this->_penalty_param > DataType(0))
            {
              DataType* constraint_at_cell = new DataType[this->get_mesh()->get_num_entities(MeshType::shape_dim)];

              this->_mesh_conc->compute_constraint(constraint_at_cell);
              exporter.add_cell_scalar("alignment_constraint", constraint_at_cell);

              delete[] constraint_at_cell;
            }
          }

          delete [] fval_norm;
          delete [] fval_det;
          delete [] fval_rec_det;

        }

        /**
         * \brief Gets the penalty parameter
         *
         * \returns The penalty parameter
         */
        DataType get_penalty_param() const
        {
          return _penalty_param;
        }

        /**
         * \brief Sets the penalty parameter
         *
         * \param[in] penalty_param_
         * The new penalty parameter.
         *
         */
        void set_penalty_param(const DataType penalty_param_)
        {
          XASSERTM(penalty_param_ >= DataType(0),"Penalty parameter must be >= 0!");
          XASSERTM(penalty_param_ <= Math::pow(Math::huge<DataType>(), DataType(0.25)),
          "Excessively large penalty parameter.");

          _penalty_param = penalty_param_;
        }

        /**
         * \brief Gets the last computed (alignment) constraint
         *
         * This is useful if you KNOW the constraint was already computed (i.e. by prepare())
         *
         * \returns The (alignment) constraint.
         */
        DataType get_constraint()
        {
          return _alignment_constraint;
        }

        //void set_mu(std::vector<CoordType>& cells, const CoordType& weight)
        //{
        //  for(Index i(0); i < cells.size(); ++i)
        //  {
        //    _mu(cells(i), weight);
        //  }

        //  _sum_mu = CoordType(0);
        //  for(Index i(0); i < _mu.size(); ++i)
        //  {
        //    _sum_mu += _mu(i);
        //  }
        //}

        /**
         * \brief Performs intialisations that cannot be done in the constructor
         *
         * \note This needs some synchronisation, so if this object lives on just one of many patches, use
         * init_pre_sync() and init_post_sync() like the Global::NonlinearFunctional class does.
         */
        virtual void init() override
        {
          init_pre_sync();
          init_post_sync();
        }

        /**
         * \brief Performs initialisations that cannot be done in the constructor, pre synchronisation phase
         *
         * \see init()
         */
        virtual void init_pre_sync()
        {
          // Write any potential changes to the mesh
          this->buffer_to_mesh();

          // Compute distances if necessary
          if((this->_scale_computation == ScaleComputation::current_concentration) ||
              (this->_scale_computation == ScaleComputation::iter_concentration) ||
              (_penalty_param > DataType(0)) )
              {
                if(_mesh_conc == nullptr)
                  throw InternalError(__func__,__FILE__,__LINE__,
                  "Scale computation set to "+stringify(_scale_computation)+"and alignment penalty factor to "+
                  stringify_fp_sci(_penalty_param)+", but no concentration function was given");

                _mesh_conc->compute_dist();
              }

          // Compute desired element size distribution
          this->_compute_scales_init();
        }

        /**
         * \brief Performs initialisations that cannot be done in the constructor, post synchronisation phase
         *
         * \see init()
         */
        virtual void init_post_sync()
        {
          // Scale mu so that || mu ||_1 = 1
          // For this, _sum_mu needs to have been summed up over all patches in the synchronisation phase.
          if(_sum_mu != DataType(1))
          {
            _mu.scale(_mu, DataType(1)/_sum_mu);
            _sum_mu = DataType(1);
          }

          // Scale lambda so that || lambda ||_1 = 1
          // For this, _sum_lambda needs to have been summed up over all patches in the synchronisation phase.
          if(_sum_lambda != DataType(1))
          {
            _lambda.scale(_lambda, DataType(1)/_sum_lambda);
            _sum_lambda = DataType(1);
          }

          // Compute the new optimal scales. For this, lambda needs to be scaled correctly and _sum_det needs to be
          // summed up over all patches in the synchronisation phase.
          RefCellTrafo_::compute_h(_h, _lambda, _sum_det);
        }

        /**
         * \brief Prepares the functional for evaluation.
         *
         * \param[in] vec_state
         * The state to evaluate the functional at.
         *
         * \param[in] filter
         * The filter. We might need to reassemble it.
         *
         * Needs to be called whenever any data like the mesh, the levelset function etc. changed.
         *
         */
        virtual void prepare(const VectorTypeR& vec_state, FilterType& filter)
        {
          prepare_pre_sync(vec_state, filter);
          prepare_post_sync(vec_state, filter);
        }

        /**
         * \brief Prepares the functional for evaluation, pre synchronisation phase
         *
         * \param[in] vec_state
         * The state to evaluate the functional at.
         *
         * \param[in] filter
         * The filter. We might need to reassemble it.
         *
         * \see prepare()
         *
         */
        virtual void prepare_pre_sync(const VectorTypeR& vec_state, FilterType& filter)
        {
          // Download state if necessary
          CoordsBufferType vec_buf;
          vec_buf.convert(vec_state);

          // Copy to buffer and mesh
          this->_coords_buffer.copy(vec_state);
          this->buffer_to_mesh();

          //auto& dirichlet_filters = filter.template at<1>();

          //for(auto& it : dirichlet_filters)
          //{
          //  const auto& assembler = _dirichlet_asm.find(it.first);

          //  if(assembler == _dirichlet_asm.end())
          //    throw InternalError(__func__,__FILE__,__LINE__,
          //    "Could not find unit filter assembler for filter with key "+it.first);

          //  assembler->second->assemble(it.second, _trafo_space);
          //}

          // The slip filter contains the outer unit normal, so reassemble it
          auto& slip_filters = filter.template at<0>();
          for(const auto& it:slip_filters)
            this->_mesh_node->adapt_by_name(it.first);

          // Copy back to the buffer as adaption might have changed the mesh
          this->mesh_to_buffer();

          for(auto& it : slip_filters)
          {
            const auto& assembler = _slip_asm.find(it.first);

            if(assembler == _slip_asm.end())
            {
              throw InternalError(__func__,__FILE__,__LINE__,
              "Could not find slip filter assembler for filter with key "+it.first);
            }

            assembler->second->assemble(it.second, _trafo_space);
          }

          if( (this->_scale_computation == ScaleComputation::iter_concentration) ||
              (_penalty_param > DataType(0)))
              {
                if(_mesh_conc == nullptr)
                {
                  throw InternalError(__func__,__FILE__,__LINE__,
                  "Scale computation set to "+stringify(_scale_computation)+"and alignment penalty factor to "+
                  stringify_fp_sci(_penalty_param)+", but no concentration function was given");
                }

                _mesh_conc->compute_dist();
              }

          _compute_scales_iter();

          if(_penalty_param > DataType(0))
          {
            _alignment_constraint = this->_mesh_conc->compute_constraint();
            this->sync_scalars.emplace("_alignment_constraint",&_alignment_constraint);
          }

        }

        /**
         * \brief Prepares the functional for evaluation, post synchronisation phase
         *
         * \param[in] vec_state
         * The state to evaluate the functional at.
         *
         * \param[in] filter
         * The filter. We might need to reassemble it.
         *
         * \see prepare()
         *
         */
        virtual void prepare_post_sync(const VectorTypeR& DOXY(vec_state), FilterType& DOXY(filter))
        {
          if(_sum_lambda != DataType(1))
          {
            _lambda.scale(_lambda, DataType(1)/_sum_lambda);
            _sum_lambda = DataType(1);
          }

          RefCellTrafo_::compute_h(_h, _lambda, _sum_det);

          if(this->_scale_computation == ScaleComputation::iter_concentration)
          {
            _mesh_conc->compute_grad_h(this->_sum_det);
          }
        }

        /**
         * \brief Computes a quality indicator concerning the cell sizes
         *
         * \param[out] lambda_min
         * Minimum of the optimal cell size lambda over all cells
         *
         * \param[out] lambda_max
         * Maximum of the optimal cell size lambda over all cells
         *
         * \param[out] vol_min
         * Minimum cell volume
         *
         * \param[out] vol_max
         * Maximum cell volume
         *
         * \param[out] vol
         * Total volume of the domain given by the mesh
         *
         * In a truly optimal mesh (consisting ONLY of reference cells of the right size), every cell's volume is
         * exactly lambda(cell). This is especially the goal for r-adaptivity.
         * So in an optimal mesh,
         * \f[
         *   \forall K \in \mathcal{T}_h: |K|/|\Omega| = \lambda(K)
         * \f]
         * so we compute the 1-norm of the vector
         * \f$(v)_i = \left| \frac{|K_i|}{\sum_j |K_j|} - \lambda(K_i) \right| \f$.
         *
         * \returns The relative cell size quality indicator.
         *
         * \note lambda_min, lambda_max, vol_min, and vol_max are all volume fractions.
         *
         * \note As these quantities are global, this function has a pre_sync and post_sync part.
         *
         * \see Control::HyperelasticityFunctionalControl::compute_cell_size_defect()
         *
         */
        virtual CoordType compute_cell_size_defect(CoordType& lambda_min, CoordType& lambda_max, CoordType& vol_min,
        CoordType& vol_max, CoordType& vol) const override
        {

          compute_cell_size_defect_pre_sync(vol_min, vol_max, vol);
          return compute_cell_size_defect_post_sync(lambda_min, lambda_max, vol_min, vol_max, vol);

        }

        /**
         * \brief Computes a quality indicator concerning the cell sizes, pre synchronisation phase
         *
         * \param[out] vol_min
         * Minimum cell volume
         *
         * \param[out] vol_max
         * Maximum cell volume
         *
         * \param[out] vol
         * Total volume of the domain given by the mesh
         *
         * \returns The relative cell size quality indicator.
         *
         * \note lambda_min, lambda_max, vol_min, and vol_max are all volume fractions.
         *
         * \note As these quantities are global, this function has a pre_sync and post_sync part.
         *
         * \see Control::HyperelasticityFunctionalControl::compute_cell_size_defect()
         *
         */
        void compute_cell_size_defect_pre_sync(CoordType& vol_min, CoordType& vol_max, CoordType& vol) const
        {
          vol = CoordType(0);

          vol_min = Math::huge<CoordType>();
          vol_max = CoordType(0);

          for(Index cell(0); cell < this->get_mesh()->get_num_entities(ShapeType::dimension); ++cell)
          {
            CoordType my_vol = this->_trafo.template compute_vol<ShapeType, CoordType>(cell);
            vol_min = Math::min(vol_min, my_vol);
            vol_max = Math::max(vol_min, my_vol);
            vol += my_vol;
          }
        }

        /**
         * \brief Computes a quality indicator concerning the cell sizes, pre synchronisation phase
         *
         * \param[out] lambda_min
         * Minimum of the optimal cell size lambda over all cells
         *
         * \param[out] lambda_max
         * Maximum of the optimal cell size lambda over all cells
         *
         * \param[out] vol_min
         * Minimum cell volume
         *
         * \param[out] vol_max
         * Maximum cell volume
         *
         * \param[out] vol
         * Total volume of the domain given by the mesh
         *
         * \returns The relative cell size quality indicator.
         *
         * \note lambda_min, lambda_max, vol_min, and vol_max are all volume fractions.
         *
         * \note As these quantities are global, this function has a pre_sync and post_sync part.
         *
         * \see Control::HyperelasticityFunctionalControl::compute_cell_size_defect()
         *
         */
        virtual CoordType compute_cell_size_defect_post_sync(CoordType& lambda_min, CoordType& lambda_max,
        CoordType& vol_min, CoordType& vol_max, const CoordType& vol) const
        {
          CoordType size_defect(0);
          lambda_min = Math::huge<CoordType>();
          lambda_max = CoordType(0);

          for(Index cell(0); cell < this->get_mesh()->get_num_entities(ShapeType::dimension); ++cell)
          {
            size_defect += Math::abs(this->_trafo.template
                compute_vol<ShapeType, CoordType>(cell)/vol - this->_lambda(cell));

            lambda_min = Math::min(lambda_min, this->_lambda(cell));
            lambda_max = Math::max(lambda_max, this->_lambda(cell));
          }

          vol_min /= vol;
          vol_max/= vol;

          return size_defect;
        }

        /**
         * \brief Evaluates the functional's value and gradient at the current state
         *
         * \param[out] fval
         * The functional value computed
         *
         * \param[out] grad
         * The functional's gradient computed
         *
         * \param[in] add_penalty_fval
         * Whether to add the penalty term
         *
         * \note The awkward handling of the penalty term is because this needs to be handled in THIS class and not
         * in another class that uses it because the the quadratic penalty method (\see Solver::QPenalty) needs to
         * call an inner solver that is unaware of the existence (or special treatment) of the penalty term.
         *
         */
        virtual void eval_fval_grad(CoordType& fval, VectorTypeL& grad, const bool& add_penalty_fval = true)
        {
          typedef typename TrafoType::template Evaluator<ShapeType, DataType>::Type TrafoEvaluator;
          typedef typename TrafoSpace::template Evaluator<TrafoEvaluator>::Type SpaceEvaluator;

          // Increase number of functional evaluations
          this->_num_func_evals++;
          this->_num_grad_evals++;

          // Total number of cells in the mesh
          Index ncells(this->get_mesh()->get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->get_mesh()->template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAT::Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;
          // This will hold the local gradient for one element for passing to other routines
          FEAT::Tiny::Matrix<CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> grad_loc;

          // Clear gradient vector
          grad.format();
          fval = DataType(0);

          TrafoEvaluator trafo_eval(this->_trafo);
          SpaceEvaluator space_eval(this->_trafo_space);

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            DataType fval_loc(0);
            trafo_eval.prepare(cell);
            space_eval.prepare(trafo_eval);

            // Get local coordinates
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              x[j] = this->_coords_buffer(idx(cell,Index(j)));
            }

            auto mat_tensor = RefCellTrafo_::compute_mat_tensor(x, this->_h(cell));

            this->_functional->eval_fval_grad(
              fval_loc, grad_loc, mat_tensor, trafo_eval, space_eval, x, this->_h(cell));

            // Add the contribution from the dependence of h on the vertex coordinates
            if(this->_mesh_conc != nullptr && this->_mesh_conc->use_derivative())
            {
              const auto& grad_h = this->_mesh_conc->get_grad_h();
              this->_functional->add_grad_h_part(
                grad_loc, mat_tensor, trafo_eval, space_eval, x, this->_h(cell), grad_h(cell));
            }

            fval += this->_mu(cell)*fval_loc;
            // Add local contributions to global gradient vector
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              Index i(idx(cell,Index(j)));
              Tiny::Vector<CoordType, MeshType::world_dim, MeshType::world_dim> tmp(grad(i));
              tmp += this->_mu(cell)*grad_loc[j];

              grad(i,tmp);
            }
          }

          if(add_penalty_fval && (this->_penalty_param > DataType(0)) )
          {
            XASSERTM(this->_mesh_conc != nullptr,
            "You need a mesh concentration function for imposing a mesh alignment constraint!");

            // Add the quadratic penalty term to the functional value
            fval += this->_penalty_param*DataType(0.5)*Math::sqr(this->_alignment_constraint);

            // Add the gradient of the penalty term
            this->_mesh_conc->add_constraint_grad(grad, this->_alignment_constraint, this->_penalty_param);
          }

        }

        /**
         * \brief Computes the functional value parts on every cell
         *
         * \param[in] fval
         * The functional value summed up over all cells.
         *
         * \param[in] fval_norm
         * Array to receive the Frobenius norm part of the functional value on every cell.
         *
         * \param[in] fval_cof
         * Array to receive the cofactor matrix part of the functional value on every cell.
         *
         * \param[in] fval_det
         * Array to receive the det part of the functional value on every cell.
         *
         */
        virtual void eval_fval_cellwise(CoordType& fval, CoordType* fval_norm, CoordType* fval_cof,
        CoordType* fval_det) const
        {
          typedef typename TrafoType::template Evaluator<ShapeType, DataType>::Type TrafoEvaluator;
          typedef typename TrafoSpace::template Evaluator<TrafoEvaluator>::Type SpaceEvaluator;

          // Total number of cells in the mesh
          Index ncells(this->get_mesh()->get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->get_mesh()->template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAT::Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;
          // Local cell dimensions for passing to other routines
          CoordType h(0);
          // This will hold the local gradient for one element for passing to other routines
          FEAT::Tiny::Matrix<CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> grad_loc;

          fval = DataType(0);

          TrafoEvaluator trafo_eval(this->_trafo);
          SpaceEvaluator space_eval(this->_trafo_space);

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            DataType fval_loc(0);
            trafo_eval.prepare(cell);
            space_eval.prepare(trafo_eval);

            h = this->_h(cell);
            // Get local coordinates
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              x[j] = this->_coords_buffer(idx(cell,Index(j)));
            }

            auto mat_tensor = RefCellTrafo_::compute_mat_tensor(x, this->_h(cell));

            this->_functional->eval_fval_cellwise(fval_loc, mat_tensor, trafo_eval, space_eval, x, h,
            fval_norm[cell], fval_cof[cell], fval_det[cell]);

            fval += this->_mu(cell)*fval_loc;
          }

        } // eval_fval_cellwise

      protected:
        /// \brief Computes the weights _lambda according to the current mesh
        virtual void _compute_lambda_cellsize()
        {
          Index ncells(this->get_mesh()->get_num_entities(ShapeType::dimension));

          _sum_lambda = DataType(0);
          for(Index cell(0); cell < ncells; ++cell)
          {
            _lambda(cell, this->_trafo.template compute_vol<ShapeType, CoordType>(cell));
            _sum_lambda += _lambda(cell);
          }

          this->sync_scalars.emplace("_sum_lambda",&_sum_lambda);

        }

        /// \brief Computes the uniformly distributed weights _lambda.
        virtual void _compute_lambda_uniform()
        {
          _lambda.format(CoordType(1));

          _sum_lambda = CoordType(_lambda.size());
          this->sync_scalars.emplace("_sum_lambda",&_sum_lambda);
        }

        /// \brief Computes _lambda according to the concentration function given by _mesh_conc
        virtual void _compute_lambda_conc()
        {
          _mesh_conc->compute_conc();
          this->sync_scalars.emplace("_sum_conc",&(_mesh_conc->get_sum_conc()));

          _mesh_conc->compute_grad_conc();

          _mesh_conc->compute_grad_sum_det(this->_coords_buffer);

          _lambda.copy(_mesh_conc->get_conc());

          _sum_lambda = DataType(0);
          for(Index cell(0); cell < this->get_mesh()->get_num_entities(ShapeType::dimension); ++cell)
          {
            _sum_lambda += _lambda(cell);
          }
          this->sync_scalars.emplace("_sum_lambda",&_sum_lambda);

        }

        /**
         * \brief Recomputes the optimal scales, every call to the solver
         */
        virtual void _compute_scales_init()
        {

          switch(_scale_computation)
          {
            case ScaleComputation::once_uniform:
            case ScaleComputation::once_cellsize:
            case ScaleComputation::once_concentration:
              return;
            case ScaleComputation::current_uniform:
              _compute_lambda_uniform();
              break;
            case ScaleComputation::current_cellsize:
              _compute_lambda_cellsize();
              break;
            case ScaleComputation::current_concentration:
            case ScaleComputation::iter_concentration:
              _compute_lambda_conc();
              break;
            default:
              return;
          }

          _sum_det = RefCellTrafo_::compute_sum_det(this->_coords_buffer, *(this->get_mesh()));
          this->sync_scalars.emplace("_sum_det",&_sum_det);

        }

        /**
         * \brief Recomputes the optimal scales, every iteration
         *
         * To be called from prepare(). This consists of first computing _lambda and then _h.
         */
        virtual void _compute_scales_iter()
        {
          switch(_scale_computation)
          {
            case ScaleComputation::iter_concentration:
              _compute_lambda_conc();
              break;
            default:
              return;
          }

          _sum_det = RefCellTrafo_::compute_sum_det(this->_coords_buffer, *(this->get_mesh()));
          this->sync_scalars.emplace("_sum_det",&_sum_det);
        }

        /**
         * \brief Computes the scales, just once
         *
         * To be called from a constructor. This consists of first computing _lambda and then _h.
         */
        virtual void _compute_scales_once()
        {
          switch(_scale_computation)
          {
            case ScaleComputation::once_uniform:
              _compute_lambda_uniform();
              break;
            case ScaleComputation::once_cellsize:
              _compute_lambda_cellsize();
              break;
            case ScaleComputation::once_concentration:
              _compute_lambda_conc();
              break;
            default:
              return;
          }

          _sum_det = RefCellTrafo_::compute_sum_det(this->_coords_buffer, *(this->get_mesh()));
          this->sync_scalars.emplace("_sum_det",&_sum_det);

        }

    }; // class HyperelasticityFunctional

#ifdef FEAT_EICKT
    /// \cond internal
    extern template class HyperelasticityFunctional
    <
      Mem::Main, double, Index,
      Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Simplex<2>, 2, 2, double >>,
      Meshopt::RumpfFunctional
      <
        double,
        Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Simplex<2>, 2, 2, double >>
      >
    >;

    extern template class HyperelasticityFunctional
    <
      Mem::Main, double, Index,
      Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Simplex<3>, 3, 3, double >>,
      Meshopt::RumpfFunctional
      <
        double,
        Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Simplex<3>, 3, 3, double >>
      >
    >;

    extern template class HyperelasticityFunctional
    <
      Mem::Main, double, Index,
      Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Hypercube<2>, 2, 2, double >>,
      Meshopt::RumpfFunctional
      <
        double,
        Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Hypercube<2>, 2, 2, double >>
      >
    >;

    extern template class HyperelasticityFunctional
    <
      Mem::Main, double, Index,
      Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Hypercube<3>, 3, 3, double >>,
      Meshopt::RumpfFunctional
      <
        double,
        Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Hypercube<3>, 3, 3, double >>
      >
    >;
    /// \endcond
#endif // FEAT_EICKT

  } // namespace Meshopt
} // namespace FEAT
#endif // KERNEL_MESHOPT_HYPERELASTICITY_FUNCTIONAL_HPP
