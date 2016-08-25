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
#include <kernel/meshopt/rumpf_functionals/2d_p1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_p1_d2.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1_d2.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1split.hpp>
#include <kernel/meshopt/rumpf_trafo.hpp>
#include <kernel/util/comm_base.hpp>

#include <map>

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
    class HyperelasticityFunctionalBase;
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
    class HyperelasticityFunctionalBase<Mem::Main, DT_, IT_, Trafo_, FunctionalType_, RefCellTrafo_>:
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

        // This is public for debugging purposes
      protected:
        /// Weights for the local contributions to the global functional value.
        ScalarVectorType _mu;
        /// Weights for local mesh size
        // In the optimal case, every cell in the mesh has the size lambda(cell)
        ScalarVectorType _lambda;
        /// Size parameters for the local reference element.
        VectorType _h;

      private:
        /// This is the number of DoFs in the trial space (= number of mesh vertices)
        const Index _columns;
        /// This is the number of DoFs in the test space (= number of mesh vertices)
        const Index _rows;
        /// How to compute the optimal scales
        ScaleComputation _scale_computation;

      protected:
        /// Factor for the alignment penalty term
        DataType _penalty_param;
        /// Last computed contraint (violation)
        DataType _alignment_constraint;

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
        explicit HyperelasticityFunctionalBase(
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
          _mu(rmn_->get_mesh()->get_num_entities(ShapeType::dimension)),
          _lambda(rmn_->get_mesh()->get_num_entities(ShapeType::dimension)),
          _h(rmn_->get_mesh()->get_num_entities(ShapeType::dimension)),
          _columns(_trafo_space.get_num_dofs()),
          _rows(_trafo_space.get_num_dofs()),
          _scale_computation(ScaleComputation::once_uniform),
          _penalty_param(0),
          _alignment_constraint(0)
          {

            XASSERTM(functional_ != nullptr, "Cell functional must not be nullptr");

            // Compute desired element size distribution
            _compute_scales_once();
            // Compute element weights
            _compute_mu();
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
         */
        explicit HyperelasticityFunctionalBase(
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
          _mu(rmn_->get_mesh()->get_num_entities(ShapeType::dimension)),
          _lambda(rmn_->get_mesh()->get_num_entities(ShapeType::dimension)),
          _h(rmn_->get_mesh()->get_num_entities(ShapeType::dimension)),
          _columns(_trafo_space.get_num_dofs()),
          _rows(_trafo_space.get_num_dofs()),
          _scale_computation(scale_computation_),
          _penalty_param(penalty_param_),
          _alignment_constraint(0)
          {
            if(( _scale_computation == ScaleComputation::once_concentration ||
                  _scale_computation == ScaleComputation::current_concentration ||
                  _scale_computation == ScaleComputation::iter_concentration ) &&
                mesh_conc_ == nullptr)
              throw InternalError(__func__,__FILE__,__LINE__,
              "Scale computation set to "+stringify(_scale_computation)+", but no concentration funtion was given");

            if(mesh_conc_ != nullptr)
            {
              _mesh_conc = mesh_conc_->create_empty_clone();
              _mesh_conc->set_mesh_node(rmn_);
            }

            // Perform one time scal computation
            _compute_scales_once();
            // Compute the cell weights
            _compute_mu();
          }

        /// Explicitly delete default constructor
        HyperelasticityFunctionalBase() = delete;
        /// Explicitly delete copy constructor
        HyperelasticityFunctionalBase(const HyperelasticityFunctionalBase&) = delete;

        /// \brief Destructor
        virtual ~HyperelasticityFunctionalBase()
        {
        };

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

        /// \copydoc BaseClass::name()
        virtual String name() const override
        {
          return "HyperelasticityFunctionalBase<"+MeshType::name()+">";
        }

        /**
         * \brief Prints some characteristics of the HyperelasticityFunctional object
         */
        virtual void print()
        {
          Util::mpi_cout_pad_line("Scale computation",_scale_computation);
          _functional->print();
          //if(_mesh_conc != nullptr)
          //  _mesh_conc->print();
        }

        /// \copydoc BaseClass::add_to_vtk_exporter()
        virtual void add_to_vtk_exporter(Geometry::ExportVTK<typename Trafo_::MeshType>& exporter) const
        {
          exporter.add_cell_vector("h", this->_h);
          exporter.add_cell_scalar("lambda", this->_lambda.elements());

          DataType* func_norm(new DataType[this->get_mesh()->get_num_entities(MeshType::shape_dim)]);
          DataType* func_det(new DataType[this->get_mesh()->get_num_entities(MeshType::shape_dim)]);
          DataType* func_rec_det(new DataType[this->get_mesh()->get_num_entities(MeshType::shape_dim)]);

          compute_func_cellwise(func_norm, func_det, func_rec_det);
          exporter.add_cell_vector("fval", func_norm, func_det, func_rec_det);

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

          delete [] func_norm;
          delete [] func_det;
          delete [] func_rec_det;

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
          XASSERTM(penalty_param_ <= Math::pow(Math::huge<DataType>(), DataType(0.25)), "Excessively large penalty parameter.");

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

        ///**
        // * \brief Computes (alignment) constraint at the current state of _mesh_conc
        // *
        // * This also saves the result to _alignment_constraint.
        // *
        // * \returns The (alignment) constraint.
        // */
        //DataType compute_constraint()
        //{
        //  XASSERT(this->_mesh_conc != nullptr);

        //  _alignment_constraint = this->_mesh_conc->compute_constraint();

        //  return _alignment_constraint;

        //}
        ///**
        // * \brief Computes the (mesh alignment) constraint on every cell
        // *
        // * \param[in] constraint_at_cells
        // * Array to receive the constraint violation on every cell
        // *
        // * \returns The sum of the constraint violation over all cells.
        // */
        //DataType compute_constraint(CoordType* constraint_at_cells)
        //{
        //  XASSERT(this->_mesh_conc != nullptr);

        //  _alignment_constraint = this->_mesh_conc->compute_constraint(constraint_at_cells);

        //  return _alignment_constraint;
        //}

        /// \copydoc BaseClass::init()
        virtual void init() override
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

          // Download state if neccessary
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
              throw InternalError(__func__,__FILE__,__LINE__,
              "Could not find slip filter assembler for filter with key "+it.first);

            assembler->second->assemble(it.second, _trafo_space);
          }

          if( (this->_scale_computation == ScaleComputation::iter_concentration) ||
              (_penalty_param > DataType(0)))
              {
                if(_mesh_conc == nullptr)
                  throw InternalError(__func__,__FILE__,__LINE__,
                  "Scale computation set to "+stringify(_scale_computation)+"and alignment penalty factor to "+
                  stringify_fp_sci(_penalty_param)+", but no concentration function was given");

                _mesh_conc->compute_dist();
              }

          if(_penalty_param > DataType(0))
            _alignment_constraint = this->_mesh_conc->compute_constraint();

          _compute_scales_iter();
        }

        /**
         * \brief Computes a quality indicator concerning the cell sizes
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
         */
        virtual CoordType compute_cell_size_defect(CoordType& lambda_min, CoordType& lambda_max,
        CoordType& vol_min, CoordType& vol_max) const override
        {
          CoordType size_defect(0);
          CoordType vol(0);

          lambda_min = Math::huge<CoordType>();
          lambda_max = CoordType(0);
          vol_min = Math::huge<CoordType>();
          vol_max = CoordType(0);

          for(Index cell(0); cell < this->get_mesh()->get_num_entities(ShapeType::dimension); ++cell)
          {
            CoordType my_vol = this->_trafo.template compute_vol<ShapeType, CoordType>(cell);
            vol_min = Math::min(vol_min, my_vol);
            vol_max = Math::max(vol_min, my_vol);
            vol += my_vol;
          }

#ifdef FEAT_HAVE_MPI
          CoordType vol_snd(vol);
          Util::Comm::allreduce(&vol_snd, &vol, Index(1), Util::CommOperationSum());

          CoordType vol_min_snd(vol_min);
          Util::Comm::allreduce(&vol_min_snd, &vol_min, Index(1), Util::CommOperationMin());

          CoordType vol_max_snd(vol_max);
          Util::Comm::allreduce(&vol_max_snd, &vol_max, Index(1), Util::CommOperationMax());
#endif

          vol_min /= vol;
          vol_max /= vol;

          for(Index cell(0); cell < this->get_mesh()->get_num_entities(ShapeType::dimension); ++cell)
          {
            size_defect += Math::abs(this->_trafo.template compute_vol<ShapeType, CoordType>(cell)/vol - this->_lambda(cell));
            lambda_min = Math::min(lambda_min, this->_lambda(cell));
            lambda_max = Math::max(lambda_max, this->_lambda(cell));
          }

#ifdef FEAT_HAVE_MPI
          CoordType size_defect_snd(size_defect);
          Util::Comm::allreduce(&size_defect_snd, &size_defect, Index(1), Util::CommOperationSum());

          CoordType lambda_min_snd(lambda_min);
          Util::Comm::allreduce(&lambda_min_snd, &lambda_min, Index(1), Util::CommOperationMin());

          CoordType lambda_max_snd(lambda_max);
          Util::Comm::allreduce(&lambda_max_snd, &lambda_max, Index(1), Util::CommOperationMax());
#endif
          return size_defect;
        }

        /**
         * \brief Computes the functional value on the current mesh.
         *
         * \returns
         * The functional value
         *
         */
        virtual CoordType compute_func() = 0;

        /**
         * \brief Computes the functional value parts on every cell
         *
         * \param[in] func_norm
         * Array to receive the Frobenius norm part of the functional value on every cell.
         *
         * \param[in] func_det
         * Array to receive the det part of the functional value on every cell.
         *
         * \param[in] func_rec_det
         * Array to receive the 1/det part of the functional value on every cell.
         *
         * \returns
         * The functional value, summed up over all parts and cells.
         */
        virtual CoordType compute_func_cellwise(CoordType* DOXY(func_norm), CoordType* DOXY(func_det), CoordType* DOXY(func_rec_det)) const = 0;

        /**
         * \brief Computes the gradient of the functional with regard to the nodal coordinates.
         *
         * Usually, prepare() should be called before calling this.
         *
         */
        virtual void compute_grad(VectorTypeL&) = 0;

        /**
         * \brief Computes the gradient of the functional with regard to the nodal coordinates.
         *
         * Usually, prepare() should be called before calling this. This const version does not increase the gradient
         * evaluation counter.
         *
         */
        virtual void compute_grad(VectorTypeL&) const = 0;

        /**
         * \brief Computes the volume of the optimal reference for each cell and saves it to _h.
         *
         */
        virtual void compute_h()
        {
          RefCellTrafo_::compute_h(_h, this->_coords_buffer, _lambda, *(this->get_mesh()));
        }

      protected:

        /// \brief Computes the weights _lambda according to the current mesh
        virtual void _compute_lambda_cellsize()
        {
          Index ncells(this->get_mesh()->get_num_entities(ShapeType::dimension));

          for(Index cell(0); cell < ncells; ++cell)
            _lambda(cell, this->_trafo.template compute_vol<ShapeType, CoordType>(cell));

        }

        /// \brief Computes the uniformly distributed weights _lambda.
        virtual void _compute_lambda_uniform()
        {
          _lambda.format(CoordType(1));
        }

        /// \brief Computes _lambda according to the concentration function given by _mesh_conc
        virtual void _compute_lambda_conc()
        {
          if(_mesh_conc == nullptr)
            throw InternalError(__func__,__FILE__,__LINE__,
            "Scale computation set to "+stringify(_scale_computation)+", but no concentration function was given");

          _mesh_conc->compute_dist();
          _mesh_conc->compute_conc();
          _mesh_conc->compute_grad_h(this->get_coords());

          _lambda.copy(_mesh_conc->get_conc());

        }
        /// \brief Computes the weights mu
        virtual void _compute_mu()
        {
          Index ncells(this->get_mesh()->get_num_entities(ShapeType::dimension));

#ifdef FEAT_HAVE_MPI
          Index ncells_send(ncells);
          Util::Comm::allreduce(&ncells_send, &ncells, 1, Util::CommOperationSum());
#endif
          _mu.format(CoordType(1)/CoordType(ncells));
        }

        /**
         * \brief Recomputes the optimal scales, every call to the solver
         *
         * To be called from init(). This consists of first computing _lambda and then _h.
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

          // Rescale so that sum lambda == 1
          CoordType sum_lambda(0);

          for(Index cell(0); cell < this->get_mesh()->get_num_entities(ShapeType::dimension); ++cell)
            sum_lambda += _lambda(cell);

#ifdef FEAT_HAVE_MPI
          CoordType sum_lambda_send(sum_lambda);
          Util::Comm::allreduce(&sum_lambda_send, &sum_lambda, 1, Util::CommOperationSum());
#endif

          _lambda.scale(_lambda, CoordType(1)/sum_lambda);

          // Compute the scales
          compute_h();
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

          // Rescale so that sum lambda == 1
          CoordType sum_lambda(0);

          for(Index cell(0); cell < this->get_mesh()->get_num_entities(ShapeType::dimension); ++cell)
            sum_lambda += _lambda(cell);

#ifdef FEAT_HAVE_MPI
          CoordType sum_lambda_send(sum_lambda);
          Util::Comm::allreduce(&sum_lambda_send, &sum_lambda, 1, Util::CommOperationSum());
#endif
          _lambda.scale(_lambda, CoordType(1)/sum_lambda);

          compute_h();
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

          // Rescale so that sum lambda == 1
          CoordType sum_lambda(0);

          for(Index cell(0); cell < this->get_mesh()->get_num_entities(ShapeType::dimension); ++cell)
            sum_lambda += _lambda(cell);

#ifdef FEAT_HAVE_MPI
          CoordType sum_lambda_send(sum_lambda);
          Util::Comm::allreduce(&sum_lambda_send, &sum_lambda, 1, Util::CommOperationSum());
#endif
          _lambda.scale(_lambda, CoordType(1)/sum_lambda);

          // Compute the scales
          compute_h();
        }

    }; // class HyperelasticityFunctionalBase

    /// \copydoc Meshopt::HyperelasticityFunctionalBase
    template
    <
      typename Mem_,
      typename DT_,
      typename IT_,
      typename Trafo_,
      typename FunctionalType_,
      typename RefCellTrafo_ = RumpfTrafo<Trafo_, typename Trafo_::MeshType::CoordType>
    >
    class HyperelasticityFunctional:
      public HyperelasticityFunctionalBase<Mem_, DT_, IT_, Trafo_, FunctionalType_, RefCellTrafo_>
    {
      public :
        /// Our base class
        typedef HyperelasticityFunctionalBase<Mem_, DT_, IT_, Trafo_, FunctionalType_, RefCellTrafo_> BaseClass;

        /// Only Mem::Main is supported atm
        typedef Mem_ MemType;
        /// We always use the precision of the mesh
        typedef DT_ DataType;
        /// We always use Index for now
        typedef IT_ IndexType;

        /// Type for the transformation
        typedef Trafo_ TrafoType;
        /// Type for the functional
        typedef FunctionalType_ FunctionalType;
        /// Type of the reference cell trafo for the mesh quality
        typedef RefCellTrafo_ RefCellTrafo;

        /// The mesh the transformation is defined on
        typedef typename TrafoType::MeshType MeshType;
        /// ShapeType of said mesh
        typedef typename MeshType::ShapeType ShapeType;
        /// The precision of the mesh coordinates
        typedef typename MeshType::CoordType CoordType;

        /// Output vector type of the operator
        typedef typename BaseClass::VectorTypeL VectorTypeL;
        /// Input vector type of the operator
        typedef typename BaseClass::VectorTypeL VectorTypeR;
        /// Type of the gradient vector
        typedef typename BaseClass::GradientType GradientType;
        /// Type for exchanging information between state variable and mesh
        typedef typename BaseClass::CoordsBufferType CoordsBufferType;

        /// Filter for Dirichlet boundary conditions
        typedef typename BaseClass::DirichletFilterType DirichletFilterType;
        /// Sequence of Dirichlet filters for several different boundary parts
        typedef typename BaseClass::DirichletFilterSequence DirichletFilterSequence;
        /// Filter for slip boundary conditions
        typedef typename BaseClass::SlipFilterType SlipFilterType;
        /// Sequence of Slip filters for several different boundary parts
        typedef typename BaseClass::SlipFilterSequence SlipFilterSequence;
        /// Combined filter
        typedef typename BaseClass::FilterType FilterType;
        /// Since the functional contains a ShapeType, these have to be the same
        static_assert( std::is_same<ShapeType, typename FunctionalType::ShapeType>::value,
        "ShapeTypes of the transformation / functional have to agree" );

      public:
        /**
         * \copydoc HyperelasticityFunctionalBase()
         */
        using BaseClass::BaseClass;

        /// Explicitly delete default constructor
        HyperelasticityFunctional() = delete;
        /// Explicitly delete copy constructor
        HyperelasticityFunctional(const HyperelasticityFunctional&) = delete;

        /// \brief Destructor
        virtual ~HyperelasticityFunctional()
        {
        };

        /// \copydoc BaseClass::name()
        virtual String name() const override
        {
          return "HyperelasticityFunctional<"+MeshType::name()+">";
        }

        /**
         * \brief Prints some characteristics of the HyperelasticityFunctional object
         */
        virtual void print() override
        {
          Util::mpi_cout(name()+" settings:\n");
          BaseClass::print();
        }

        /// \copydoc HyperelasticityFunctionalBase::compute_func()
        virtual CoordType compute_func() override
        {
          // Increase number of functional evaluations
          this->_num_func_evals++;

          DataType fval(static_cast<const HyperelasticityFunctional*>(this)->_compute_func_intern());

          if(this->_penalty_param > DataType(0))
          {
            XASSERT(this->_mesh_conc != nullptr);
            fval += this->_penalty_param*DataType(0.5)*Math::sqr(this->_alignment_constraint);
          }
          return fval;
        }

        /// \copydoc compute_func(CoordType*,CoordType*,CoordType*)
        virtual CoordType compute_func_cellwise(CoordType* func_norm, CoordType* func_det, CoordType* func_rec_det) const override
        {
          CoordType fval(0);
          // Total number of cells in the mesh
          Index ncells(this->get_mesh()->get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->get_mesh()->template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAT::Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;
          // Local cell dimensions for passing to other routines
          FEAT::Tiny::Vector<CoordType, MeshType::world_dim> h;

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
              x[j] = this->_coords_buffer(idx(cell,Index(j)));

            // Scale local functional value with lambda
            fval += this->_mu(cell) * this->_functional->compute_local_functional(x,h, norm_A, det_A, rec_det_A);

            func_norm[cell] = this->_mu(cell) * norm_A;
            func_det[cell] = this->_mu(cell) * det_A;
            func_rec_det[cell] = this->_mu(cell) * rec_det_A;
            func_norm_tot += func_norm[cell];
            func_det_tot += func_det[cell];
            func_rec_det_tot += func_rec_det[cell];
          }

          return fval;
        } // compute_func


        /// \copydoc BaseClass::compute_grad(VectorTypeL&)
        virtual void compute_grad(VectorTypeL& grad) override
        {
          // Increase number of functional evaluations
          this->_num_grad_evals++;

          static_cast<const HyperelasticityFunctional*>(this)->compute_grad(grad);
        }

        /// \copydoc BaseClass::compute_grad(VectorTypeL&)
        virtual void compute_grad(VectorTypeL& grad) const override
        {
          if(this->_mesh_conc != nullptr && this->_mesh_conc->use_derivative())
            _compute_grad_with_conc(grad);
          else
            _compute_grad_without_conc(grad);

          if(this->_penalty_param > DataType(0))
          {
            XASSERTM(this->_mesh_conc != nullptr,
            "You need a mesh concentration function for imposing a mesh alignment constraint!");
            this->_mesh_conc->add_constraint_grad(grad, this->_alignment_constraint, this->_penalty_param);
          }
        }

      protected:
        /// \copydoc HyperelasticityFunctionalBase::compute_func()
        virtual CoordType _compute_func_intern() const
        {
          CoordType fval(0);
          // Total number of cells in the mesh
          Index ncells(this->get_mesh()->get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          const auto& idx = this->get_mesh()->template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAT::Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;
          // Local cell dimensions for passing to other routines
          FEAT::Tiny::Vector<CoordType, MeshType::world_dim> h;

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            h = this->_h(cell);
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; j++)
              x[j] = this->_coords_buffer(idx(cell,Index(j)));

            fval += this->_mu(cell) * this->_functional->compute_local_functional(x,h);
          }

          return fval;
        } // compute_func


        /// \copydoc BaseClass::compute_grad()
        virtual void _compute_grad_with_conc(VectorTypeL& grad) const
        {
          XASSERTM(this->_mesh_conc != nullptr, "You need a mesh concentration function for this.");

          // Total number of cells in the mesh
          Index ncells(this->get_mesh()->get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->get_mesh()->template get_index_set<ShapeType::dimension,0>();

          const auto& grad_h = this->_mesh_conc->get_grad_h();

          // This will hold the coordinates for one element for passing to other routines
          FEAT::Tiny::Matrix<CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;
          // Local cell dimensions for passing to other routines
          FEAT::Tiny::Vector<CoordType, MeshType::world_dim> h;
          // This will hold the local gradient for one element for passing to other routines
          FEAT::Tiny::Matrix<CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> grad_loc;

          // Clear gradient vector
          grad.format();

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            h = this->_h(cell);
            // Get local coordinates
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
              x[j] = this->_coords_buffer(idx(cell,Index(j)));

            this->_functional->compute_local_grad(x, h, grad_loc);
            // Add the contribution from the dependence of h on the vertex coordinates
            this->_functional->add_grad_h_part(grad_loc, x, h, grad_h(cell));

            // Add local contributions to global gradient vector
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              Index i(idx(cell,Index(j)));
              Tiny::Vector<CoordType, MeshType::world_dim, MeshType::world_dim> tmp(grad(i));
              tmp += this->_mu(cell)*grad_loc[j];

              grad(i,tmp);
            }
          }

        } // compute_grad_with_conc

        /// \copydoc BaseClass::compute_grad()
        virtual void _compute_grad_without_conc(VectorTypeL& grad) const
        {
          // Total number of cells in the mesh
          Index ncells(this->get_mesh()->get_num_entities(ShapeType::dimension));

          // Index set for local/global numbering
          auto& idx = this->get_mesh()->template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAT::Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;
          // Local cell dimensions for passing to other routines
          FEAT::Tiny::Vector<CoordType, MeshType::world_dim> h;
          // This will hold the local gradient for one element for passing to other routines
          FEAT::Tiny::Matrix<CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> grad_loc;

          // Clear gradient vector
          grad.format();

          // Compute the functional value for each cell
          for(Index cell(0); cell < ncells; ++cell)
          {
            h = this->_h(cell);
            // Get local coordinates
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
              x[j] = this->_coords_buffer(idx(cell,Index(j)));

            this->_functional->compute_local_grad(x, h, grad_loc);

            // Add local contributions to global gradient vector
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              Index i(idx(cell,Index(j)));
              Tiny::Vector<CoordType, MeshType::world_dim, MeshType::world_dim> tmp(grad(i));
              tmp += this->_mu(cell)*grad_loc[j];

              grad(i,tmp);
            }
          }

        } // compute_grad_without_conc

    }; // class HyperelasticityFunctional

    /// \cond internal
    extern template class HyperelasticityFunctional
    <
      Mem::Main, double, Index,
      Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Simplex<2>, 2, 2, double >>,
      Meshopt::RumpfFunctional<double,Shape::Simplex<2>>
    >;

    extern template class FEAT::Meshopt::HyperelasticityFunctional
    <
      Mem::Main, double, Index,
      Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2>, 2, 2, double>>,
      Meshopt::RumpfFunctional<double,Shape::Hypercube<2>>
    >;

    extern template class FEAT::Meshopt::HyperelasticityFunctional
    <
      Mem::Main, double, Index,
      Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Simplex<2>, 2, 2, double>>,
      Meshopt::RumpfFunctional_D2<double,Shape::Simplex<2>>
    >;

    extern template class FEAT::Meshopt::HyperelasticityFunctional
    <
      Mem::Main, double, Index,
      Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2>, 2, 2, double>>,
      Meshopt::RumpfFunctional_D2<double,Shape::Hypercube<2>>
    >;

    extern template class FEAT::Meshopt::HyperelasticityFunctional
    <
      Mem::Main, double, Index,
      Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Simplex<2>, 2, 2, double >>,
      Meshopt::RumpfFunctionalQ1Split<double, Shape::Simplex<2>, FEAT::Meshopt::RumpfFunctional>
    >;

    extern template class FEAT::Meshopt::HyperelasticityFunctional
    <
      Mem::Main, double, Index,
      Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Simplex<2>, 2, 2, double >>,
      Meshopt::RumpfFunctionalQ1Split<double, Shape::Simplex<2>, FEAT::Meshopt::RumpfFunctional_D2>
    >;

    extern template class FEAT::Meshopt::HyperelasticityFunctional
    <
      Mem::Main, double, Index,
      Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Hypercube<2>, 2, 2, double >>,
      Meshopt::RumpfFunctionalQ1Split<double, Shape::Hypercube<2>, FEAT::Meshopt::RumpfFunctional>
    >;

    extern template class FEAT::Meshopt::HyperelasticityFunctional
    <
      Mem::Main, double, Index,
      Trafo::Standard::Mapping<Geometry::ConformalMesh< Shape::Hypercube<2>, 2, 2, double >>,
      Meshopt::RumpfFunctionalQ1Split<double, Shape::Hypercube<2>, FEAT::Meshopt::RumpfFunctional_D2>
    >;
    /// \endcond

  } // namespace Meshopt
} // namespace FEAT
#endif // KERNEL_MESHOPT_HYPERELASTICITY_FUNCTIONAL_HPP
