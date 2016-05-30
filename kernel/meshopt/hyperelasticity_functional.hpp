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
#include <kernel/meshopt/mesh_quality_functional.hpp>
#include <kernel/meshopt/rumpf_trafo.hpp>
#include <kernel/util/comm_base.hpp>

#include <map>

namespace FEAT
{
  namespace Meshopt
  {

    /**
     * \brief Enum class for different types of scale computations
     */
    enum class ScaleComputation
    {
      undefined = 0,
      once_uniform,
      once_initial,
      once_chart_distance,
      current_initial,
      current_chart_distance,
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
        case ScaleComputation::once_initial:
          return os << "once_initial";
        case ScaleComputation::once_chart_distance:
          return os << "once_chart_distance";
        case ScaleComputation::current_initial:
          return os << "current_initial";
        case ScaleComputation::current_chart_distance:
          return os << "current_chart_distance";
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
      else if(sc_name == "once_initial")
        scale_computation = ScaleComputation::once_initial;
      else if(sc_name == "once_chart_distance")
        scale_computation = ScaleComputation::once_chart_distance;
      else if(sc_name == "current_initial")
        scale_computation = ScaleComputation::current_initial;
      else if(sc_name == "current_chart_distance")
        scale_computation = ScaleComputation::current_chart_distance;
      else if(sc_name == "iter_concentration")
        scale_computation = ScaleComputation::iter_concentration;
      else
        throw InternalError(__func__, __FILE__, __LINE__, "Unknown ScaleComputation identifier string "
            +sc_name);
    }
    /// \endcond

    template
    <
      typename Mem_,
      typename DT_,
      typename IT_,
      typename TrafoType_,
      typename FunctionalType_,
      typename RefCellTrafo_ = RumpfTrafo<TrafoType_, typename TrafoType_::MeshType::CoordType>
    >
    class HyperelasticityFunctionalBase;
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
     * Local meshsize evaluator. \see RumpfTrafo
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
      typename DT_,
      typename IT_,
      typename TrafoType_,
      typename FunctionalType_,
      typename RefCellTrafo_
    >
    class HyperelasticityFunctionalBase<Mem::Main, DT_, IT_, TrafoType_, FunctionalType_, RefCellTrafo_>:
      public MeshQualityFunctional<typename TrafoType_::MeshType>
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

        // This is public for debugging purposes
      public:
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

        ScaleComputation _scale_computation;
        std::deque<String> _distance_charts;

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
          _mu(rmn_->get_mesh()->get_num_entities(ShapeType::dimension)),
          _lambda(rmn_->get_mesh()->get_num_entities(ShapeType::dimension)),
          _h(rmn_->get_mesh()->get_num_entities(ShapeType::dimension)),
          _columns(_trafo_space.get_num_dofs()),
          _rows(_trafo_space.get_num_dofs()),
          _scale_computation(ScaleComputation::once_uniform),
          _distance_charts()
          {
            // Compute desired element size distribution
            this->compute_scales_once();
            // Compute element weights
            this->compute_mu();
          }

        explicit HyperelasticityFunctionalBase(
          Geometry::RootMeshNode<MeshType>* rmn_,
          TrafoSpace& trafo_space_,
          std::map<String, std::shared_ptr<Assembly::UnitFilterAssembler<MeshType>>>& dirichlet_asm_,
          std::map<String, std::shared_ptr<Assembly::SlipFilterAssembler<MeshType>>>& slip_asm_,
          std::shared_ptr<FunctionalType_> functional_,
          ScaleComputation scale_computation_,
          const std::deque<String>& distance_charts_)
          : BaseClass(rmn_),
          _trafo(trafo_space_.get_trafo()),
          _trafo_space(trafo_space_),
          _dirichlet_asm(dirichlet_asm_),
          _slip_asm(slip_asm_),
          _functional(functional_),
          _mu(rmn_->get_mesh()->get_num_entities(ShapeType::dimension)),
          _lambda(rmn_->get_mesh()->get_num_entities(ShapeType::dimension)),
          _h(rmn_->get_mesh()->get_num_entities(ShapeType::dimension)),
          _columns(_trafo_space.get_num_dofs()),
          _rows(_trafo_space.get_num_dofs()),
          _scale_computation(scale_computation_),
          _distance_charts()
          {
            for(const auto& it:distance_charts_)
              _distance_charts.push_back(it);

            if(_scale_computation == ScaleComputation::once_chart_distance ||
                _scale_computation == ScaleComputation::current_chart_distance)
                {
                  if(_distance_charts.empty())
                    throw InternalError(__func__,__FILE__,__LINE__,
                    "Scale computation set to "+stringify(_scale_computation)+", but no charts were given");
                }
            // Perform one time scal computation
            compute_scales_once();
            // Compute the cell weights
            compute_mu();
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
         * \brief Checks if the functional is empty (= the null functional
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
          return "HyperelasticityFunctionalBase<"+MeshType::name()+">";
        }

        /**
         * \brief Prints some characteristics of the HyperelasticityFunctional object
         */
        virtual void print()
        {
          std::cout << "scale computation:" << _scale_computation << std::endl;
          if(!_dirichlet_asm.empty())
          {
            std::cout << "Dirichlet boundaries:";
              for(const auto& it:_dirichlet_asm)
                std::cout << " " << it.first;
            std::cout << std::endl;
          }
          _functional->print();
        }

        virtual void add_to_vtk_exporter(Geometry::ExportVTK<typename TrafoType_::MeshType>& exporter) const
        {
          exporter.add_cell_vector("h", this->_h);
          exporter.add_cell_scalar("lambda", this->_lambda.elements());

          DataType* func_norm(new DataType[this->get_mesh()->get_num_entities(MeshType::shape_dim)]);
          DataType* func_det(new DataType[this->get_mesh()->get_num_entities(MeshType::shape_dim)]);
          DataType* func_rec_det(new DataType[this->get_mesh()->get_num_entities(MeshType::shape_dim)]);

          compute_func(func_norm, func_det, func_rec_det);
          exporter.add_cell_vector("fval", func_norm, func_det, func_rec_det);

          delete [] func_norm;
          delete [] func_det;
          delete [] func_rec_det;

        }

        void add_distance_chart(const String& chart_name)
        {
          if(this->_mesh_node->get_atlas()->find_chart(chart_name) == nullptr)
            throw InternalError(__func__,__FILE__,__LINE__,"Chart "+chart_name+" not found in atlas");

          _distance_charts.push_back(chart_name);
        }

        virtual void init() override
        {
          // Write any potential changes to the mesh
          this->buffer_to_mesh();

          // Compute desired element size distribution
          this->compute_scales_init();
        }

        /**
         * \brief Prepares the functional for evaluation.
         *
         * Needs to be called whenever any data like the mesh, the levelset function etc. changed.
         *
         */
        virtual void prepare(const VectorTypeR& vec_state, FilterType& filter)
        {
          this->_coords_buffer.convert(vec_state);
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

          for(auto& it : slip_filters)
          {
            const auto& assembler = _slip_asm.find(it.first);

            if(assembler == _slip_asm.end())
              throw InternalError(__func__,__FILE__,__LINE__,
              "Could not find slip filter assembler for filter with key "+it.first);

            assembler->second->assemble(it.second, _trafo_space);
          }

          compute_scales_iter();

          this->mesh_to_buffer();
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
        virtual CoordType compute_func() const = 0;
        virtual CoordType compute_func(CoordType*, CoordType*, CoordType*) const = 0;

        /**
         * \brief Computes the gradient of the functional with regard to the nodal coordinates.
         *
         * Usually, prepare() should be called before calling this.
         *
         */
        virtual void compute_grad(VectorTypeL&) const = 0;

        /**
         * \brief Computes the volume of the optimal reference for each cell and saves it to _h.
         *
         */
        virtual void compute_h()
        {
          RefCellTrafo_::compute_h(_h, this->_coords_buffer, _lambda, this->_trafo);
        }

        virtual void compute_scales_iter()
        {
          switch(_scale_computation)
          {
            case ScaleComputation::iter_concentration:
              break;
            default:
              return;
          }

          // Rescale so that sum lambda == 1
          CoordType sum_lambda(0);

          for(Index cell(0); cell < this->get_mesh()->get_num_entities(ShapeType::dimension); ++cell)
            sum_lambda += _lambda(cell);

          CoordType sum_lambda_send(sum_lambda);
          Util::Comm::allreduce(&sum_lambda, 1, &sum_lambda_send);

          _lambda.scale(_lambda, CoordType(1)/sum_lambda);

          compute_h();
        }

        virtual void compute_scales_init()
        {
          switch(_scale_computation)
          {
            case ScaleComputation::once_uniform:
            case ScaleComputation::once_initial:
            case ScaleComputation::once_chart_distance:
              return;
            case ScaleComputation::current_initial:
              compute_lambda_initial();
              break;
            case ScaleComputation::current_chart_distance:
              compute_lambda_chart_dist();
              break;
            default:
              throw InternalError(__func__,__FILE__,__LINE__,
              "Unhandled scale computation"+stringify(_scale_computation));
          }

          // Rescale so that sum lambda == 1
          CoordType sum_lambda(0);

          for(Index cell(0); cell < this->get_mesh()->get_num_entities(ShapeType::dimension); ++cell)
            sum_lambda += _lambda(cell);

          CoordType sum_lambda_send(sum_lambda);
          Util::Comm::allreduce(&sum_lambda, 1, &sum_lambda_send);

          _lambda.scale(_lambda, CoordType(1)/sum_lambda);

          // Compute the scales
          compute_h();
        }

        /// \brief Computes the weights _lambda.
        virtual void compute_scales_once()
        {
          switch(_scale_computation)
          {
            case ScaleComputation::once_uniform:
              compute_lambda_uniform();
              break;
            case ScaleComputation::once_initial:
              compute_lambda_initial();
              break;
            case ScaleComputation::once_chart_distance:
              compute_lambda_chart_dist();
              break;
            default:
              return;
          }

          // Rescale so that sum lambda == 1
          CoordType sum_lambda(0);

          for(Index cell(0); cell < this->get_mesh()->get_num_entities(ShapeType::dimension); ++cell)
            sum_lambda += _lambda(cell);

          CoordType sum_lambda_send(sum_lambda);
          Util::Comm::allreduce(&sum_lambda, 1, &sum_lambda_send);

          _lambda.scale(_lambda, CoordType(1)/sum_lambda);

          // Compute the scales
          compute_h();
        }

        /// \brief Computes the weights _lambda according to the initial mesh
        virtual void compute_lambda_initial()
        {
          Index ncells(this->get_mesh()->get_num_entities(ShapeType::dimension));

          for(Index cell(0); cell < ncells; ++cell)
            _lambda(cell, this->_trafo.template compute_vol<ShapeType, CoordType>(cell));

        }

        /// \brief Computes the uniformly distributed weights _lambda.
        virtual void compute_lambda_uniform()
        {
          Index ncells(this->get_mesh()->get_num_entities(ShapeType::dimension));

          _lambda.format(CoordType(1)/CoordType(ncells));
        }

        virtual void compute_lambda_chart_dist()
        {
          Index ncells(this->get_mesh()->get_num_entities(ShapeType::dimension));
          //_lambda.format(Math::huge<CoordType>());
          _lambda.format(CoordType(0));

          const auto& idx = this->get_mesh()->template get_index_set<ShapeType::dimension, 0>();
          const auto& vtx = this->get_mesh()->get_vertex_set();

          Tiny::Vector<CoordType, MeshType::world_dim> midpoint(CoordType(0));

          for(Index cell(0); cell < ncells; ++cell)
          {
            midpoint.format(CoordType(0));
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              Index i(idx(cell, Index(j)));
              midpoint += vtx[i];
            }
            midpoint *= (DataType(1))/DataType(Shape::FaceTraits<ShapeType,0>::count);

            CoordType midpoint_dist(0);
            for(const auto& it:_distance_charts)
            {
              auto* chart = this->_mesh_node->get_atlas()->find_mesh_chart(it);
              if(chart == nullptr)
                throw InternalError(__func__,__FILE__,__LINE__,"Could not find chart "+it);

              midpoint_dist += chart->dist(midpoint);

            }
            _lambda(cell, midpoint_dist);
          }
        }

        /// \brief Computes the weights mu
        virtual void compute_mu()
        {
          Index ncells(this->get_mesh()->get_num_entities(ShapeType::dimension));

          Index ncells_send(ncells);
          Util::Comm::allreduce(&ncells, 1, &ncells_send);

          _mu.format(CoordType(1)/CoordType(ncells));
        }

    }; // class HyperelasticityFunctionalBase

    /// \copydoc Meshopt::HyperelasticityFunctionalBase
    template
    <
      typename Mem_,
      typename DT_,
      typename IT_,
      typename TrafoType_,
      typename FunctionalType_,
      typename RefCellTrafo_ = RumpfTrafo<TrafoType_, typename TrafoType_::MeshType::CoordType>
    >
    class HyperelasticityFunctional:
      public HyperelasticityFunctionalBase<Mem_, DT_, IT_, TrafoType_, FunctionalType_, RefCellTrafo_>
    {
      public :
        /// Our base class
        typedef HyperelasticityFunctionalBase<Mem_, DT_, IT_, TrafoType_, FunctionalType_, RefCellTrafo_> BaseClass;

        /// Type for the transformation
        typedef TrafoType_ TrafoType;

        /// The mesh the transformation is defined on
        typedef typename TrafoType::MeshType MeshType;
        /// The precision of the mesh coordinates
        typedef typename MeshType::CoordType CoordType;

        /// Only Mem::Main is supported atm
        typedef Mem_ MemType;
        /// We always use the precision of the mesh
        typedef DT_ DataType;
        /// We always use Index for now
        typedef IT_ IndexType;

        /// Type for the functional
        typedef FunctionalType_ FunctionalType;
        /// ShapeType of said mesh
        typedef typename MeshType::ShapeType ShapeType;

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
         *
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
        virtual void print() override
        {
          std::cout << name() << std::endl;
          BaseClass::print();
        }

        /// \copydoc HyperelasticityFunctionalBase::compute_func()
        virtual CoordType compute_func()
        {
          // Increase number of functional evaluations
          this->_num_func_evals++;

          return static_cast<const HyperelasticityFunctional*>(this)->compute_func();
        }

        /// \copydoc HyperelasticityFunctionalBase::compute_func()
        virtual CoordType compute_func() const override
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
        virtual CoordType compute_func(CoordType* func_norm, CoordType* func_det, CoordType* func_rec_det)
        {
          // Increase number of functional evaluations
          this->_num_func_evals++;

          return static_cast<const HyperelasticityFunctional*>(this)->compute_func(func_norm, func_det, func_rec_det);
        }

        virtual CoordType compute_func(CoordType* func_norm, CoordType* func_det, CoordType* func_rec_det) const override
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

        /// \copydoc BaseClass::compute_grad()
        virtual void compute_grad(VectorTypeL& grad)
        {
          // Increase number of functional evaluations
          this->_num_grad_evals++;

          static_cast<const HyperelasticityFunctional*>(this)->compute_grad(grad);

        }

        /// \copydoc BaseClass::compute_grad()
        virtual void compute_grad(VectorTypeL& grad) const override
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

        } // compute_grad

    }; // class HyperelasticityFunctional

  } // namespace Meshopt
} // namespace FEAT
#endif // KERNEL_MESHOPT_HYPERELASTICITY_FUNCTIONAL_HPP
