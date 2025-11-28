#include <kernel/base_header.hpp>
#include <kernel/runtime.hpp>

#include <applications/gendie/fbm_control_helper.hpp>
#include <applications/gendie/steady_flow_solver.hpp>
#include <applications/gendie/steady_stokes_solver.hpp>
#include <applications/gendie/defect_assembler.hpp>
#include <applications/gendie/system_assembler.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/util/dist_file_io.hpp>
#include <kernel/util/dist.hpp>
#include <applications/gendie/gendie_common.hpp>
#include <kernel/assembly/trace_assembler.hpp>
#include <kernel/assembly/trace_assembler_basic_jobs.hpp>
#include <kernel/assembly/surface_integrator.hpp>
#include <kernel/assembly/surface_integrator_basic_jobs.hpp>
#include <kernel/geometry/facet_flipper.hpp>
#include <kernel/geometry/hit_test_factory.hpp>
#include <kernel/geometry/intern/meshpart_converter.hpp>
#include <kernel/geometry/atlas/bezier.hpp>
#include <kernel/util/xml_scanner.hpp>

#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/dense_vector.hpp>

// open the namespace and define the DomainLevel and SystemLevel classes
namespace FBMTest
{
  using namespace FEAT;

  namespace Intern
  {
    // template for the steady-state inflow function used in bench1 and bench2
    template<typename DT_, int dim_>
    class InflowHelper :
      public Analytic::Function
    {
    public:
      static constexpr int domain_dim = dim_;
      typedef Analytic::Image::Vector<dim_> ImageType;
      static constexpr bool can_value = true;

    protected:
      DT_ _vmax;

    public:
      explicit InflowHelper(DT_ vmax) :
        _vmax(vmax)
      {
      }

      template<typename Traits_>
      class Evaluator :
        public Analytic::Function::Evaluator<Traits_>
      {
      protected:
        typedef typename Traits_::DataType DataType;
        typedef typename Traits_::PointType PointType;
        typedef typename Traits_::ValueType ValueType;

        const DataType _vmax;

      public:
        explicit Evaluator(const InflowHelper& function) :
          _vmax(DataType(function._vmax))
        {
        }

        ValueType value(const PointType& point)
        {
          ValueType val;
          val[0] = _vmax * DataType(4) * (point[1])*(0.41-point[1])/(0.41*0.41);
          for(int l = 1; l < dim_; ++l)
            val[l] = DataType(0);
          return val;
        }
      };
    }; // class InflowHelper

    template<typename FaceVector_, typename FEVector_, typename Space_>
    class FaceValueTraceAssemblyJob
    {
    public:
      typedef typename FEVector_::DataType DataType;
      typedef typename FEVector_::ValueType ValueType;
      typedef typename FaceVector_::ValueType FaceValueType;

      static constexpr TrafoTags trafo_config = TrafoTags::none;
      static constexpr SpaceTags space_config = SpaceTags::value;
      static constexpr TrafoTags facet_trafo_config = TrafoTags::jac_det | TrafoTags::normal;

      class Task :
        public FEAT::Assembly::TraceAssemblyBasicTaskBase1<DataType, Space_, trafo_config, facet_trafo_config, space_config>
      {
      public:
        /// our base-class typedef
        typedef FEAT::Assembly::TraceAssemblyBasicTaskBase1<DataType, Space_, trafo_config, facet_trafo_config, space_config> BaseClass;

        /// our assembly traits
        typedef typename BaseClass::AsmTraits AsmTraits;

        // we do not support pairwise assembly
        static constexpr bool assemble_pairwise = false;
        /// this task needs to scatter
        static constexpr bool need_scatter = false;
        /// this task has no combine
        static constexpr bool need_combine = false;

        using typename BaseClass::TrafoEvalData;
        using typename BaseClass::TrafoFacetEvalData;

      protected:
        FaceVector_& face_vector;
        /// the vector that is to be integrated
        const FEVector_& vector;
        static constexpr int fe_dim = FEVector_::BlockSize;
        /// the cubature rule used for integration
        typename BaseClass::CubatureRuleType cubature_rule;
        /// the local vector to be assembled
        typename AsmTraits::template TLocalVector<ValueType> local_vector;
        /// the vector gather object
        typename FEVector_::GatherAxpy gather_axpy;
        Tiny::Vector<DataType, fe_dim> loc_v;
        FaceValueType loc_face_val;

      public:
        explicit Task(FaceValueTraceAssemblyJob& job) :
          BaseClass(job.space),
          face_vector(job.face_vector),
          vector(job.vector),
          cubature_rule(Cubature::ctor_factory, job.cubature_factory),
          local_vector(),
          gather_axpy(vector),
          loc_v(),
          loc_face_val()
        {
        }

        /// \brief Performs the local assembly.
        void assemble()
        {
          // format local vector
          local_vector.format();

          // gather local vector data
          gather_axpy(local_vector, this->dof_mapping);

          // fetch number of local dofs
          const int num_loc_dofs = this->space_eval.get_num_local_dofs();
          loc_face_val = FaceValueType(DataType(0));

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // prepare point
            this->prepare_point(cubature_rule.get_point(k));

            loc_v.format();
            // calculate local gradient
            for(int l = 0; l < num_loc_dofs; ++l)
            {
              loc_v.axpy(this->space_data.phi[l].value, local_vector[l]);
            }

            // do the dirty work
            Tiny::axpy(loc_face_val, loc_v, cubature_rule.get_weight(k)*this->trafo_facet_data.jac_det);
          } // continue with next cubature point

          // scatter our face value
          face_vector(this->trafo_facet_eval.get_cell_index(), (DataType(1) / this->trafo_facet_eval.volume()) * loc_face_val);
        }

        /// \brief Scatters the local assembly.
        void scatter()
        {
          // nothing to do here
        }

        /// \brief Finalizes the assembly.
        void combine()
        {
        }
      }; // class Task

    protected:
      FaceVector_& face_vector;
      /// a reference to the vector that is to be integrated
      const FEVector_& vector;
      /// a reference to the finite element space to be used as test-/trial-space
      const Space_& space;
      /// the cubature factory to be used for integration
      Cubature::DynamicFactory cubature_factory;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] vector_
       * A \resident reference to the coefficient vector of the discrete function.
       *
       * \param[in] space_
       * A \resident reference to the finite element space used for the discretization
       *
       * \param[in] cubature_
       * The name of the cubature rule to use for integration.
       */
      explicit FaceValueTraceAssemblyJob(FaceVector_& face_vec, const FEVector_& vector_, const Space_& space_, const String& cubature_):
        face_vector(face_vec),
        vector(vector_),
        space(space_),
        cubature_factory(cubature_)
      {
      }

    }; // class FaceValueTraceAssemblyJob<...>

    template<typename FaceVector_, typename FEVector_, typename Space_>
    class FaceGradientTraceAssemblyJob
    {
    public:
      typedef typename FEVector_::DataType DataType;
      typedef typename FEVector_::ValueType ValueType;
      typedef typename FaceVector_::ValueType FaceValueType;

      static constexpr TrafoTags trafo_config = TrafoTags::none;
      static constexpr SpaceTags space_config = SpaceTags::grad;
      static constexpr TrafoTags facet_trafo_config = TrafoTags::jac_det | TrafoTags::normal;

      class Task :
        public FEAT::Assembly::TraceAssemblyBasicTaskBase1<DataType, Space_, trafo_config, facet_trafo_config, space_config>
      {
      public:
        /// our base-class typedef
        typedef FEAT::Assembly::TraceAssemblyBasicTaskBase1<DataType, Space_, trafo_config, facet_trafo_config, space_config> BaseClass;

        /// our assembly traits
        typedef typename BaseClass::AsmTraits AsmTraits;

        // we do not support pairwise assembly
        static constexpr bool assemble_pairwise = false;
        /// this task needs to scatter
        static constexpr bool need_scatter = false;
        /// this task has no combine
        static constexpr bool need_combine = false;

        using typename BaseClass::TrafoEvalData;
        using typename BaseClass::TrafoFacetEvalData;

      protected:
        FaceVector_& face_vector;
        /// the vector that is to be integrated
        const FEVector_& vector;
        static constexpr int fe_dim = FEVector_::BlockSize;
        /// the cubature rule used for integration
        typename BaseClass::CubatureRuleType cubature_rule;
        /// the local vector to be assembled
        typename AsmTraits::template TLocalVector<ValueType> local_vector;
        /// the vector gather object
        typename FEVector_::GatherAxpy gather_axpy;
        Tiny::Matrix<DataType, fe_dim, fe_dim> loc_grad;
        FaceValueType loc_face_val;

      public:
        explicit Task(FaceGradientTraceAssemblyJob& job) :
          BaseClass(job.space),
          face_vector(job.face_vector),
          vector(job.vector),
          cubature_rule(Cubature::ctor_factory, job.cubature_factory),
          local_vector(),
          gather_axpy(vector),
          loc_grad(),
          loc_face_val()
        {
        }

        /// \brief Performs the local assembly.
        void assemble()
        {
          // format local vector
          local_vector.format();

          // gather local vector data
          gather_axpy(local_vector, this->dof_mapping);

          // fetch number of local dofs
          const int num_loc_dofs = this->space_eval.get_num_local_dofs();
          loc_face_val = FaceValueType(DataType(0));

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // prepare point
            this->prepare_point(cubature_rule.get_point(k));

            loc_grad.format();
            // calculate local gradient
            for(int l = 0; l < num_loc_dofs; ++l)
            {
              loc_grad.add_outer_product(local_vector[l], this->space_data.phi[l].grad);
            }

            // construct normal gradient
            auto normal_grad = loc_grad * this->trafo_facet_data.normal;

            // do the dirty work
            Tiny::axpy(loc_face_val, normal_grad, cubature_rule.get_weight(k)*this->trafo_facet_data.jac_det);
          } // continue with next cubature point
          // scatter our face value
          face_vector(this->trafo_facet_eval.get_cell_index(), (DataType(1) / this->trafo_facet_eval.volume()) * loc_face_val);
        }

        /// \brief Scatters the local assembly.
        void scatter()
        {
          // nothing to do here
        }

        /// \brief Finalizes the assembly.
        void combine()
        {
        }
      }; // class Task

    protected:
      FaceVector_& face_vector;
      /// a reference to the vector that is to be integrated
      const FEVector_& vector;
      /// a reference to the finite element space to be used as test-/trial-space
      const Space_& space;
      /// the cubature factory to be used for integration
      Cubature::DynamicFactory cubature_factory;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] vector_
       * A \resident reference to the coefficient vector of the discrete function.
       *
       * \param[in] space_
       * A \resident reference to the finite element space used for the discretization
       *
       * \param[in] cubature_
       * The name of the cubature rule to use for integration.
       */
      explicit FaceGradientTraceAssemblyJob(FaceVector_& face_vec, const FEVector_& vector_, const Space_& space_, const String& cubature_):
        face_vector(face_vec),
        vector(vector_),
        space(space_),
        cubature_factory(cubature_)
      {
      }

    }; // class FaceGradientTraceAssemblyJob<...>

    template<typename FaceVector_, typename Space_, typename FEVector_>
    class FaceValueIntegratorJob
    {
    public:
      typedef Space_ SpaceType;
      typedef FEVector_ FEVector;
      typedef FaceVector_ FaceVector;
      /// the data-type of the vector
      typedef typename FEVector_::DataType DataType;
      /// the value-type of the vector
      typedef typename FEVector_::ValueType ValueType;
      typedef typename FaceVector_::ValueType FaceValueType;

      static constexpr TrafoTags trafo_config_ = TrafoTags::img_point | TrafoTags::dom_point;
      static constexpr SpaceTags space_config_ = SpaceTags::value;

      /// our assembly traits
      typedef Assembly::AsmTraits1<
        DataType,
        Space_,
        trafo_config_,
        space_config_
      > AsmTraits;




      class Task :
        public FEAT::Assembly::FEVectorIntegratorTaskCRTP<Task, AsmTraits, Space_, FEVector_>
      {
      public:
        typedef FEAT::Assembly::FEVectorIntegratorTaskCRTP<Task, AsmTraits, Space_, FEVector_> BaseClass;
        typedef typename BaseClass::IndexType IndexType;
        typedef typename BaseClass::DomainPointType DomainPointType;
        static constexpr int dim = BaseClass::dim;
        /// the vector that is to be assembled
        FaceVector_& vector;
        FaceValueType face_val;
        /// the scatter scaling factor
        DataType scatter_alpha;

      public:
        /**
        * \brief Constructor
        */
        explicit Task(FaceValueIntegratorJob& job) :
                BaseClass(job.space, job.primal_vec),
                vector(job.face_vec),
                face_val(),
                scatter_alpha(job.scatter_alpha)
        {
        }

        void _assemble()
        {
          face_val.format();

          BaseClass::_assemble();

          this->vector(this->_cur_surface_index, (DataType(1) / this->_face_volume) * face_val);
        }

        void _integrate(DataType weight, [[maybe_unused]] IndexType point_idx)
        {
          Tiny::axpy(face_val, this->loc_value_holder.value, weight * scatter_alpha);
        }

      }; // class Task

      FaceVector_& face_vec;
      const FEVector_& primal_vec;
      const Space_& space;
      DataType scatter_alpha;

      explicit FaceValueIntegratorJob(FaceVector_& face_vec_, const FEVector_& primal_vec_, const Space_& space_, DataType scatter_alpha_ = DataType(1)) :
        face_vec(face_vec_),
        primal_vec(primal_vec_),
        space(space_),
        scatter_alpha(scatter_alpha_)
      {
      }
    };

    template<typename FaceVector_, typename Space_, typename FEVector_>
    class FaceGradientIntegratorJob
    {
    public:
      typedef Space_ SpaceType;
      typedef FEVector_ FEVector;
      typedef FaceVector_ FaceVector;
      /// the data-type of the vector
      typedef typename FEVector_::DataType DataType;
      /// the value-type of the vector
      typedef typename FEVector_::ValueType ValueType;
      typedef typename FaceVector_::ValueType FaceValueType;

      static constexpr TrafoTags trafo_config_ = TrafoTags::img_point | TrafoTags::dom_point;
      static constexpr SpaceTags space_config_ = SpaceTags::grad;

      /// our assembly traits
      typedef Assembly::AsmTraits1<
        DataType,
        Space_,
        trafo_config_,
        space_config_
      > AsmTraits;




      class Task :
        public FEAT::Assembly::FEVectorIntegratorTaskCRTP<Task, AsmTraits, Space_, FEVector_>
      {
      public:
        typedef FEAT::Assembly::FEVectorIntegratorTaskCRTP<Task, AsmTraits, Space_, FEVector_> BaseClass;
        typedef typename BaseClass::IndexType IndexType;
        typedef typename BaseClass::DomainPointType DomainPointType;
        static constexpr int dim = BaseClass::dim;
        /// the vector that is to be assembled
        FaceVector_& vector;
        FaceValueType face_val;
        /// the scatter scaling factor
        DataType scatter_alpha;

      public:
        /**
        * \brief Constructor
        */
        explicit Task(FaceGradientIntegratorJob& job) :
                BaseClass(job.space, job.primal_vec),
                vector(job.face_vec),
                face_val(),
                scatter_alpha(job.scatter_alpha)
        {
        }

        void _assemble()
        {
          face_val.format();

          BaseClass::_assemble();

          this->vector(this->_cur_surface_index, (DataType(1) / this->_face_volume) * face_val);
        }

        void _integrate(DataType weight, [[maybe_unused]] IndexType point_idx)
        {
          auto normal_grad = this->loc_value_holder.grad * this->_normal;
          Tiny::axpy(face_val, normal_grad, weight * scatter_alpha);
        }

      }; // class Task

      FaceVector_& face_vec;
      const FEVector_& primal_vec;
      const Space_& space;
      DataType scatter_alpha;

      explicit FaceGradientIntegratorJob(FaceVector_& face_vec_, const FEVector_& primal_vec_, const Space_& space_, DataType scatter_alpha_ = DataType(1)) :
        face_vec(face_vec_),
        primal_vec(primal_vec_),
        space(space_),
        scatter_alpha(scatter_alpha_)
      {
      }
    };

    template<typename FaceVector_, typename Space_, typename FEVector_, typename MaskSpace_, typename MaskVector_>
    class FaceGradientFBMDeltaPeakIntegratorJob
    {
    public:
      typedef Space_ SpaceType;
      typedef MaskSpace_ MaskSpaceType;
      typedef FEVector_ FEVector;
      typedef FaceVector_ FaceVector;
      typedef MaskVector_ MaskVector;
      /// the data-type of the vector
      typedef typename FEVector_::DataType DataType;
      /// the value-type of the vector
      typedef typename FEVector_::ValueType ValueType;
      typedef typename FaceVector_::ValueType FaceValueType;

      static constexpr TrafoTags trafo_config_ = TrafoTags::img_point | TrafoTags::dom_point | TrafoTags::jac_det | TrafoTags::jac_mat;
      static constexpr SpaceTags space_config_ = SpaceTags::grad;
      static constexpr SpaceTags mask_config_ = SpaceTags::ref_grad | SpaceTags::grad;

      /// our assembly traits
      typedef Assembly::AsmTraits1<
        DataType,
        Space_,
        trafo_config_,
        space_config_
      > AsmTraits;

      /// our assembly traits
      typedef Assembly::AsmTraits1<
        DataType,
        MaskSpace_,
        trafo_config_,
        mask_config_
      > MaskAsmTraits;



      class Task :
        public FEAT::Assembly::SurfaceIntegratorTaskBase<AsmTraits>
      {
      public:
        typedef FEAT::Assembly::SurfaceIntegratorTaskBase<AsmTraits> BaseClass;
        typedef typename BaseClass::IndexType IndexType;
        typedef typename BaseClass::DomainPointType DomainPointType;
        static constexpr int dim = BaseClass::dim;
        typedef typename AsmTraits::DataType DataType;
        typedef typename FEVector_::ValueType ValueType;
        static constexpr int image_dim = FEAT::Assembly::Intern::ValueTypeHelper<ValueType>::dim;

        /// the vector that is to be assembled
        FaceVector_& vector;
        const FEVector& fe_vector;
        const MaskVector& mask_vector;

        const Space_& space;
        const MaskSpace_& mask_space;
        typename AsmTraits::TrafoEvaluator trafo_eval;
        /// the space evaluator
        typename AsmTraits::SpaceEvaluator space_eval;
        typename MaskAsmTraits::SpaceEvaluator mask_space_eval;
        /// the space dof-mapping
        typename AsmTraits::DofMapping dof_mapping;
        typename MaskAsmTraits::DofMapping mask_dof_mapping;
        /// the trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;
        /// the space evaluation data
        typename AsmTraits::SpaceEvalData space_data;
        typename MaskAsmTraits::SpaceEvalData mask_space_data;
        /// the local vector to be gathered
        typename AsmTraits::template TLocalVector<ValueType> local_vector;
        typename MaskAsmTraits::template TLocalVector<typename MaskAsmTraits::DataType> local_mask_vector;
        /// the gather object
        typename FEVector_::GatherAxpy gather_axpy;
        typename MaskVector_::GatherAxpy mask_gather_axpy;
        /// local fe point/grad/hess values
        FEAT::Assembly::Intern::LocalFEValueHolder<AsmTraits, DataType, dim, image_dim> loc_value_holder;
        FEAT::Assembly::Intern::LocalFEValueHolder<AsmTraits, DataType, dim, 1> loc_mask_value_holder;

        Cubature::Rule<typename AsmTraits::ShapeType> cubature_rule;

        FaceValueType face_val;
        FaceValueType loc_cell_val;
        /// the scatter scaling factor
        DataType scatter_alpha;

      public:
        /**
        * \brief Constructor
        */
        explicit Task(FaceGradientFBMDeltaPeakIntegratorJob& job) :
                BaseClass(job.space.get_trafo()),
                vector(job.face_vec),
                fe_vector(job.primal_vec),
                mask_vector(job.mask_vec),
                space(job.space),
                mask_space(job.mask_space),
                trafo_eval(this->_trafo),
                space_eval(space),
                mask_space_eval(mask_space),
                dof_mapping(space),
                mask_dof_mapping(mask_space),
                gather_axpy(fe_vector),
                mask_gather_axpy(mask_vector),
                cubature_rule(Cubature::ctor_factory, job.cubature_factory),
                face_val(),
                loc_cell_val(),
                scatter_alpha(job.scatter_alpha)
        {
        }


        template<typename ValueHolder_, typename LocVecType_, typename SpaceEval_, typename SpaceData_>
        void gather_point_values(ValueHolder_& val_holder_, const LocVecType_& loc_vec_, const SpaceEval_& space_eval_, const SpaceData_& space_data_)
        {
          const int num_loc_dofs = space_eval_.get_num_local_dofs();
          if constexpr(ValueHolder_::has_value)
          {
            val_holder_.value = typename ValueHolder_::ValueType(0);
            for(int i = 0; i < num_loc_dofs; ++i)
            {
              Tiny::axpy(val_holder_.value, loc_vec_[i], space_data_.phi[i].value);
            }
          }
          if constexpr(ValueHolder_::has_grad)
          {
            val_holder_.grad.format();
            for(int i = 0; i < num_loc_dofs; ++i)
            {
              if constexpr(ValueHolder_::image_dim == 1)
                Tiny::axpy(val_holder_.grad, space_data_.phi[i].grad, loc_vec_[i]);
              else
                val_holder_.grad.add_outer_product(loc_vec_[i], space_data_.phi[i].grad);
            }
          }
          if constexpr(ValueHolder_::has_hess)
          {
            XABORTM("Hessian not implemented yet");
          }

        }

        template<typename ValueHolder_, typename LocVecType_, typename SpaceEval_, typename SpaceData_>
        void gather_ref_point_values(ValueHolder_& val_holder_, const LocVecType_& loc_vec_, const SpaceEval_& space_eval_, const SpaceData_& space_data_)
        {
          const int num_loc_dofs = space_eval_.get_num_local_dofs();
          if constexpr(ValueHolder_::has_value)
          {
            val_holder_.value = typename ValueHolder_::ValueType(0);
            for(int i = 0; i < num_loc_dofs; ++i)
            {
              Tiny::axpy(val_holder_.value, loc_vec_[i], space_data_.phi[i].ref_value);
            }
          }
          if constexpr(ValueHolder_::has_grad)
          {
            val_holder_.grad.format();
            for(int i = 0; i < num_loc_dofs; ++i)
            {
              if constexpr(ValueHolder_::image_dim == 1)
                Tiny::axpy(val_holder_.grad, space_data_.phi[i].ref_grad, loc_vec_[i]);
              else
                val_holder_.grad.add_outer_product(loc_vec_[i], space_data_.phi[i].ref_grad);
            }
          }
          if constexpr(ValueHolder_::has_hess)
          {
            XABORTM("Hessian not implemented yet");
          }

        }

        void _eval_point_cell_val(DataType point_weight)
        {
          Tiny::axpy(loc_cell_val, loc_value_holder.grad*this->_normal, point_weight);
        }

        void _integrate_local_cell()
        {
          loc_cell_val = FaceValueType(0);

          const DataType cell_vol = trafo_eval.volume();

          for(int pt = 0; pt < cubature_rule.get_num_points(); ++pt)
          {
            auto dom_point = cubature_rule.get_point(pt);
            DataType weight = cubature_rule.get_weight(pt);
            trafo_eval(trafo_data, dom_point);
            space_eval(space_data, trafo_data);
            mask_space_eval(mask_space_data, trafo_data);
            // construct local fe values
            gather_point_values(loc_value_holder, local_vector, space_eval, space_data);
            gather_ref_point_values(loc_mask_value_holder, local_mask_vector, mask_space_eval, mask_space_data);

            // construct our mask gradient norm, in normal direction, integrating our folding has to
            // evaluate to one, which corresponds to one h term of the gradient
            // but since we evaluate a point representing the whole integral, we have to correct for the volume
            // of the d-1 manifold. Its easier to simply evaluate the reference gradient instead and divide by the volume
            // of our d-dimensional cell
            auto gradient_norm = DataType(2) * loc_mask_value_holder.grad.norm_euclid() / cell_vol;

            this->_eval_point_cell_val(weight*gradient_norm*trafo_data.jac_det);

          }
        }

        void _prepare_cell(IndexType cell)
        {
          dof_mapping.prepare(cell);
          mask_dof_mapping.prepare(cell);
          trafo_eval.prepare(cell);
          space_eval.prepare(trafo_eval);
          mask_space_eval.prepare(trafo_eval);
          local_vector.format();
          gather_axpy(local_vector, dof_mapping);
          local_mask_vector.format();
          mask_gather_axpy(local_mask_vector, mask_dof_mapping);
          this->_integrate_local_cell();
        }

        void _assemble_cell(const std::vector<IndexType>& domain_point_idx)
        {
          for(auto pti : domain_point_idx)
          {
            Tiny::axpy(face_val, loc_cell_val, this->_point_weights[pti] / this->_face_volume);
          }

        }


      public:
        void assemble()
        {
          face_val = DataType(0);
          for(IndexType k=0; k < this->_cell_helper.size(); ++k)
          {
            const auto& domain_point_idx = this->_cell_to_domain_point.at(k);
            if(domain_point_idx.empty())
              continue;

            this->_prepare_cell(this->_cell_helper[k]);
            this->_assemble_cell(domain_point_idx);
          }
          this->vector(this->_cur_surface_index, face_val*scatter_alpha);


        }
      }; // class Task

      FaceVector_& face_vec;
      const FEVector_& primal_vec;
      const MaskVector& mask_vec;
      const Space_& space;
      const MaskSpace_& mask_space;
      Cubature::DynamicFactory cubature_factory;
      DataType scatter_alpha;

      explicit FaceGradientFBMDeltaPeakIntegratorJob(FaceVector_& face_vec_, const FEVector_& primal_vec_, const MaskVector_& mask_vec_, const Space_& space_, const MaskSpace_& mask_space_, DataType scatter_alpha_ = DataType(1)) :
        face_vec(face_vec_),
        primal_vec(primal_vec_),
        mask_vec(mask_vec_),
        space(space_),
        mask_space(mask_space_),
        cubature_factory("gauss-legendre:3"),
        scatter_alpha(scatter_alpha_)
      {
      }
    };

    template<typename FaceVector_, typename Space_, typename FEVector_>
    class FaceGradientCellMeanIntegratorJob
    {
    public:
      typedef Space_ SpaceType;
      typedef FEVector_ FEVector;
      typedef FaceVector_ FaceVector;
      /// the data-type of the vector
      typedef typename FEVector_::DataType DataType;
      /// the value-type of the vector
      typedef typename FEVector_::ValueType ValueType;
      typedef typename FaceVector_::ValueType FaceValueType;

      static constexpr TrafoTags trafo_config_ = TrafoTags::img_point | TrafoTags::dom_point | TrafoTags::jac_det | TrafoTags::jac_mat;
      static constexpr SpaceTags space_config_ = SpaceTags::grad;

      /// our assembly traits
      typedef Assembly::AsmTraits1<
        DataType,
        Space_,
        trafo_config_,
        space_config_
      > AsmTraits;

      class Task :
        public FEAT::Assembly::SurfaceIntegratorTaskBase<AsmTraits>
      {
      public:
        typedef FEAT::Assembly::SurfaceIntegratorTaskBase<AsmTraits> BaseClass;
        typedef typename BaseClass::IndexType IndexType;
        typedef typename BaseClass::DomainPointType DomainPointType;
        static constexpr int dim = BaseClass::dim;
        typedef typename AsmTraits::DataType DataType;
        typedef typename FEVector_::ValueType ValueType;
        static constexpr int image_dim = FEAT::Assembly::Intern::ValueTypeHelper<ValueType>::dim;

        /// the vector that is to be assembled
        FaceVector_& vector;
        const FEVector& fe_vector;

        const Space_& space;
        typename AsmTraits::TrafoEvaluator trafo_eval;
        /// the space evaluator
        typename AsmTraits::SpaceEvaluator space_eval;
        /// the space dof-mapping
        typename AsmTraits::DofMapping dof_mapping;
        /// the trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;
        /// the space evaluation data
        typename AsmTraits::SpaceEvalData space_data;
        /// the local vector to be gathered
        typename AsmTraits::template TLocalVector<ValueType> local_vector;
        /// the gather object
        typename FEVector_::GatherAxpy gather_axpy;
        /// local fe point/grad/hess values
        FEAT::Assembly::Intern::LocalFEValueHolder<AsmTraits, DataType, dim, image_dim> loc_value_holder;

        Cubature::Rule<typename AsmTraits::ShapeType> cubature_rule;

        FaceValueType face_val;
        FaceValueType loc_cell_val;
        /// the scatter scaling factor
        DataType scatter_alpha;

      public:
        /**
        * \brief Constructor
        */
        explicit Task(FaceGradientCellMeanIntegratorJob& job) :
                BaseClass(job.space.get_trafo()),
                vector(job.face_vec),
                fe_vector(job.primal_vec),
                space(job.space),
                trafo_eval(this->_trafo),
                space_eval(space),
                dof_mapping(space),
                gather_axpy(fe_vector),
                cubature_rule(Cubature::ctor_factory, job.cubature_factory),
                face_val(),
                loc_cell_val(),
                scatter_alpha(job.scatter_alpha)
        {
        }


        template<typename ValueHolder_, typename LocVecType_, typename SpaceEval_, typename SpaceData_>
        void gather_point_values(ValueHolder_& val_holder_, const LocVecType_& loc_vec_, const SpaceEval_& space_eval_, const SpaceData_& space_data_)
        {
          const int num_loc_dofs = space_eval_.get_num_local_dofs();
          if constexpr(ValueHolder_::has_value)
          {
            val_holder_.value = typename ValueHolder_::ValueType(0);
            for(int i = 0; i < num_loc_dofs; ++i)
            {
              Tiny::axpy(val_holder_.value, loc_vec_[i], space_data_.phi[i].value);
            }
          }
          if constexpr(ValueHolder_::has_grad)
          {
            val_holder_.grad.format();
            for(int i = 0; i < num_loc_dofs; ++i)
            {
              if constexpr(ValueHolder_::image_dim == 1)
                Tiny::axpy(val_holder_.grad, space_data_.phi[i].grad, loc_vec_[i]);
              else
                val_holder_.grad.add_outer_product(loc_vec_[i], space_data_.phi[i].grad);
            }
          }
          if constexpr(ValueHolder_::has_hess)
          {
            XABORTM("Hessian not implemented yet");
          }

        }

        void _eval_point_cell_val(DataType point_weight)
        {
          Tiny::axpy(loc_cell_val, loc_value_holder.grad*this->_normal, point_weight);
        }

        void _integrate_local_cell()
        {
          loc_cell_val = FaceValueType(0);

          const DataType cell_vol = trafo_eval.volume();

          //integrate value over cell
          for(int pt = 0; pt < cubature_rule.get_num_points(); ++pt)
          {
            auto dom_point = cubature_rule.get_point(pt);
            DataType weight = cubature_rule.get_weight(pt);
            trafo_eval(trafo_data, dom_point);
            space_eval(space_data, trafo_data);
            // construct local fe values
            gather_point_values(loc_value_holder, local_vector, space_eval, space_data);

            this->_eval_point_cell_val(weight*trafo_data.jac_det);
          }
          // calculate cell mean value
          loc_cell_val *= DataType(1) / cell_vol;
        }

        void _prepare_cell(IndexType cell)
        {
          dof_mapping.prepare(cell);
          trafo_eval.prepare(cell);
          space_eval.prepare(trafo_eval);
          local_vector.format();
          gather_axpy(local_vector, dof_mapping);
          this->_integrate_local_cell();
        }

        void _assemble_cell(const std::vector<IndexType>& domain_point_idx)
        {
          for(auto pti : domain_point_idx)
          {
            Tiny::axpy(face_val, loc_cell_val, this->_point_weights[pti] / this->_face_volume);
          }

        }


      public:
        void assemble()
        {
          face_val = DataType(0);
          for(IndexType k=0; k < this->_cell_helper.size(); ++k)
          {
            const auto& domain_point_idx = this->_cell_to_domain_point.at(k);
            if(domain_point_idx.empty())
              continue;

            this->_prepare_cell(this->_cell_helper[k]);
            this->_assemble_cell(domain_point_idx);
          }
          this->vector(this->_cur_surface_index, face_val*scatter_alpha);


        }
      }; // class Task

      FaceVector_& face_vec;
      const FEVector_& primal_vec;
      const Space_& space;
      Cubature::DynamicFactory cubature_factory;
      DataType scatter_alpha;

      explicit FaceGradientCellMeanIntegratorJob(FaceVector_& face_vec_, const FEVector_& primal_vec_, const Space_& space_, DataType scatter_alpha_ = DataType(1)) :
        face_vec(face_vec_),
        primal_vec(primal_vec_),
        space(space_),
        cubature_factory("gauss-legendre:3"),
        scatter_alpha(scatter_alpha_)
      {
      }
    };
  }

  /// define application dimension
  #ifndef FEAT_APP_DIM
  #define FEAT_APP_DIM = 2
  #endif
  static constexpr int dim = FEAT_APP_DIM;
  static_assert((dim == 2) || (dim == 3), "invalid dimension");

  /// output padding length
  // static constexpr std::size_t pad_len = 30u;

  /// output padding character
  // static constexpr char pad_char = '.';

  /// helper function to parse arguments
  template<typename T_>
  T_ parse(SimpleArgParser& args, const String& name, T_ default_value)
  {
    args.parse(name, default_value);
    return default_value;
  }

  /// helper to parse backend string
  PreferredBackend parse_backend(const String& backend_string)
  {
    if(backend_string.compare_no_case("cuda") == 0)
    {
      return PreferredBackend::cuda;
    }
    else if(backend_string.compare_no_case("mkl") == 0)
    {
      return PreferredBackend::mkl;
    }
    else if(backend_string.compare_no_case("generic") == 0)
    {
      return PreferredBackend::generic;
    }
    XABORTM("Unknown backend " + backend_string);
    return PreferredBackend::generic;

  }

  /// our one and only index type
  typedef Index IndexType;

  // depending on whether FEAT_CCND_APP_QUADMATH is defined,
  // we use quad precision or double precision
#ifdef FEAT_CCND_APP_QUADMATH
#  define Q_(x) (x##Q)
  typedef __float128 DataType;
  // static constexpr int fp_num_digs = 35;
  // static const char* fp_typename = "quadruple";
#else
#  define Q_(x) x
  typedef double DataType;
  // static constexpr int fp_num_digs = 15;
  // static const char* fp_typename = "double";
#endif

  // our matrix types
  typedef LAFEM::SparseMatrixBCSR<DataType, IndexType, dim, dim> LocalMatrixBlockA;
  typedef LAFEM::SparseMatrixBCSR<DataType, IndexType, dim, 1> LocalMatrixBlockB;
  typedef LAFEM::SparseMatrixBCSR<DataType, IndexType, 1, dim> LocalMatrixBlockD;
  typedef LAFEM::SparseMatrixCSR<DataType, IndexType> LocalScalarMatrix;

  // our vector types
  typedef LAFEM::DenseVectorBlocked<DataType, IndexType, dim> LocalVeloVector;
  typedef LAFEM::DenseVector<DataType, IndexType> LocalPresVector;

  // define our mesh type and other geometry related classes
  typedef Shape::Hypercube<dim> ShapeType;
  typedef Geometry::ConformalMesh<ShapeType, dim, DataType> MeshType;
  typedef Geometry::MeshPart<MeshType> MeshPartType;
  typedef Geometry::RootMeshNode<MeshType> MeshNodeType;
  typedef Geometry::Atlas::ChartBase<MeshType> ChartBaseType;
  typedef Geometry::MeshAtlas<MeshType> MeshAtlasType;
  typedef typename Shape::FaceTraits<ShapeType, dim-1>::ShapeType FaceType;

  // define our trafo type: standard or isoparametric
#ifdef FEAT_APP_ISOPARAM
  // static constexpr bool isoparam = true;
  typedef Trafo::Isoparam::Mapping<MeshType, 2> TrafoType;
#else
  // static constexpr bool isoparam = false;
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
#endif

  // define FE space types
#if defined(FEAT_APP_Q1T_P0)
  typedef Space::CroRavRanTur::Element<TrafoType> SpaceVeloType;
  typedef Space::Discontinuous::ElementP0<TrafoType> SpacePresType;
#elif defined(FEAT_APP_Q1TBNP_P1DC)
  typedef Space::Q1TBNP::Element<TrafoType> SpaceVeloType;
  typedef Space::Discontinuous::ElementP1<TrafoType> SpacePresType;
#else
  typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;
  typedef Space::Discontinuous::ElementP1<TrafoType> SpacePresType;
#endif


  /// our domain level
  typedef Control::Domain::StokesDomainLevel<MeshType, TrafoType, SpaceVeloType, SpacePresType> DomainLevel;
  typedef Control::StokesFBMDomainLevelBase<SpaceVeloType, SpacePresType> FBMDomainLevel;

  /// our system level
  typedef Control::StokesSystemLevelBase<dim, LocalMatrixBlockA, LocalMatrixBlockB, LocalMatrixBlockD, LocalScalarMatrix> SystemLevel;
  typedef Control::StokesFBMSystemLevelBase<dim, LocalMatrixBlockA, LocalMatrixBlockB, LocalMatrixBlockD, LocalScalarMatrix> FBMSystemLevel;

  /// our partidomaincontrol
  template<typename DomainLevel_>
  using PartitionControl = Control::Domain::PartiDomainControl<DomainLevel_>;

  String cubature_matrix_a = "gauss-legendre:3";
  String cubature_matrix_b = "gauss-legendre:3";
  String cubature_matrix_m = "gauss-legendre:3";
  String cubature_transfer = "gauss-legendre:3";
  String cubature_defect = "gauss-legendre:3";

  DataType nu = 0.001;
  DataType beta = 1.;
  bool deformation = false;
  bool navier = true;
  DataType abs_tol = 1E-10;

  template<typename Domain_>
  void read_mesh_file(Domain_& domain, const Dist::Comm& comm, SimpleArgParser& args, Gendie::Logger& logger)
  {
    if(args.check("mesh") < 1)
    {
      XABORTM("You have to provide a meshfile");
    }

    std::deque<String> mesh_files = args.query("mesh")->second;

    logger.print("Reading from files " + std::accumulate(mesh_files.begin(), mesh_files.end(), String(""), [](const FEAT::String& a, const FEAT::String& b){return a + " " + b;}), Gendie::info);

    Geometry::MeshFileReader mesh_reader;
    mesh_reader.add_mesh_files(comm, mesh_files);
    mesh_reader.read_root_markup();

    domain.create(mesh_reader);

    domain.add_trafo_mesh_part_charts();
  }

  template<typename Domain_, typename SystemLevel_, bool fbm_>
  void create_systems(Domain_& domain, std::deque<std::shared_ptr<SystemLevel_>>& system)
  {
    const Index num_levels(domain.size_physical());
    system.resize(num_levels);
    for(Index i(0); i < num_levels; ++i)
    {
      system.at(i).reset(new SystemLevel_());
      // get domain and system levels
      auto& dom_lvl = *domain.at(i);
      auto& sys_lvl = *system.at(i);
      dom_lvl.domain_asm.set_max_worker_threads(std::size_t(0));

      // compile domain assembler
      dom_lvl.domain_asm.compile_all_elements();

      // assemble gates
      sys_lvl.assemble_gates(domain.at(i));
    }

    // assemble muxers and transfers
    for (Index i(0); (i < domain.size_physical()) && ((i+1) < domain.size_virtual()); ++i)
    {
      system.at(i)->assemble_coarse_muxers(domain.at(i+1));
      if((i+1) < domain.size_physical())
        system.at(i)->assemble_transfers(*system.at(i+1), domain.at(i), domain.at(i+1), cubature_transfer, true);
      else
        system.at(i)->assemble_transfers(domain.at(i), domain.at(i+1), cubature_transfer, true);
    }

    for(Index i(0); i < num_levels; ++i)
    {
      // get domain and system levels
      auto& dom_lvl = *domain.at(i);
      auto& sys_lvl = *system.at(i);


      // assemble matrix structures
      sys_lvl.assemble_velo_struct(dom_lvl.space_velo);
      sys_lvl.assemble_pres_struct(dom_lvl.space_pres);

      // assemble matrices
      {
        auto& loc_a = sys_lvl.matrix_a.local();
        loc_a.format();
        Assembly::Common::LaplaceOperatorBlocked<dim> lapl_op;
        Assembly::Common::DuDvOperatorBlocked<dim> dudv_op;
        String cubature = cubature_matrix_a;
        if(cubature.empty())
          cubature = "auto-degree:" + stringify(2*SpaceVeloType::local_degree+2);
        if(deformation)
          Assembly::assemble_bilinear_operator_matrix_1(dom_lvl.domain_asm, loc_a, dudv_op, dom_lvl.space_velo, cubature, nu);
        else
          Assembly::assemble_bilinear_operator_matrix_1(dom_lvl.domain_asm, loc_a, lapl_op, dom_lvl.space_velo, cubature, nu);
      }
      if constexpr(fbm_)
      {
        sys_lvl.velo_mass_matrix = sys_lvl.matrix_a.clone(LAFEM::CloneMode::Weak);
        auto& loc_m = sys_lvl.velo_mass_matrix.local();
        loc_m.format();
        String cubature = cubature_matrix_m;
        if(cubature.empty())
          cubature = "auto-degree:" + stringify(2*SpaceVeloType::local_degree+2);
        Assembly::Common::IdentityOperatorBlocked<dim> id_op;
        Assembly::assemble_bilinear_operator_matrix_1(dom_lvl.domain_asm, loc_m, id_op, dom_lvl.space_velo, cubature);
      }
      sys_lvl.assemble_grad_div_matrices(dom_lvl.domain_asm, dom_lvl.space_velo, dom_lvl.space_pres, cubature_matrix_b);

      // compile the system matrix
      sys_lvl.compile_system_matrix();
    }

  }

  template<typename Domain_, typename SystemLevel_, bool fbm_ = false>
  void assemble_filters(const std::deque<String>& inflows, const std::deque<String>& outflows, std::deque<std::shared_ptr<SystemLevel_>>& system, const Domain_& domain, DataType vmax = 0.3, bool use_q2_fbm = false)
  {
    // the names of the mesh parts on which to assemble
    std::deque<String> part_names = domain.front()->get_mesh_node()->get_mesh_part_names(true);
    for(std::size_t i(0); i < system.size(); ++i)
    {
      // get our local velocity filters
      auto& filter_v_noflow = system.at(i)->get_local_velo_unit_filter_seq().find_or_add("noflow");
      auto& filter_v_inflow = system.at(i)->get_local_velo_unit_filter_seq().find_or_add("inflow");

      // create unit-filter assembler
      Assembly::UnitFilterAssembler<MeshType> /*unit_asm_inflow,*/ unit_asm_noflow;
      // loop over all boundary parts except for the right one, which is outflow
      for(const auto& name : part_names)
      {
        // skip non-boundary mesh-parts
        if((!name.starts_with("bnd:")))
          continue;

        // try to fetch the corresponding mesh part node
        auto* mesh_part_node = domain.at(i)->get_mesh_node()->find_mesh_part_node(name);
        XASSERT(mesh_part_node != nullptr);

        auto* mesh_part = mesh_part_node->get_mesh();
        String name_lower = name.lower();
        if (mesh_part != nullptr)
        {
          if(std::find(inflows.begin(), inflows.end(), name)!=inflows.end())
          {
            Assembly::UnitFilterAssembler<MeshType> loc_unit_asm_inflow;

            // inflow
            loc_unit_asm_inflow.add_mesh_part(*mesh_part);
            loc_unit_asm_inflow.assemble(filter_v_inflow, domain.at(i)->space_velo, Intern::InflowHelper<DataType, dim>(vmax));
          }
          else if(std::find(outflows.begin(), outflows.end(), name) == outflows.end())
          {
            // outflow
            unit_asm_noflow.add_mesh_part(*mesh_part);
          }
        }
      }
      {
        unit_asm_noflow.assemble(filter_v_noflow, domain.at(i)->space_velo);
      }
      // assemble fbm parts
      if constexpr(fbm_)
      {
        auto& fbm_asm = *domain.at(i)->fbm_assembler;

        system.at(i)->assemble_fbm_filters(fbm_asm, domain.at(i)->space_velo, domain.at(i)->space_pres, (i==0u), use_q2_fbm, false);
      }
      // compile system filter
      system.at(i)->compile_system_filter();
    }
  }

  enum surface_integration_type : std::uint8_t
  {
    point_values = 0,
    delta_peak = 1,
    cell_mean = 2,
    unkown = 3
  };

  static inline surface_integration_type parse_integration_type(const FEAT::String& type_string)
  {
    if((type_string.compare_no_case("point") == 0) || (type_string.compare_no_case("point_value") == 0) || (type_string.compare_no_case("point-value") == 0))
    {
      return surface_integration_type::point_values;
    }
    else if((type_string.compare_no_case("delta") == 0) || (type_string.compare_no_case("peak") == 0) || (type_string.compare_no_case("delta-peak") == 0)|| (type_string.compare_no_case("delta_peak") == 0))
    {
      return surface_integration_type::delta_peak;
    }
    else if((type_string.compare_no_case("cell") == 0) || (type_string.compare_no_case("mean") == 0) || (type_string.compare_no_case("cell-mean") == 0)|| (type_string.compare_no_case("cell_mean") == 0))
    {
      return surface_integration_type::cell_mean;
    }

    return unkown;
  }

  static inline FEAT::String print_integration_type(surface_integration_type type)
  {
    switch (type)
    {
      case point_values:
        return FEAT::String("Point Values");
      case delta_peak:
        return FEAT::String("Delta Peak");
      case cell_mean:
        return FEAT::String("Cell Mean");
      case unkown:
        return FEAT::String("unkown");
    }

    return FEAT::String();
  }

  void run_fbm(Dist::Comm& comm, SimpleArgParser& args, Gendie::Logger& logger)
  {
    logger.print("Running fbm case", Gendie::info);
    args.support("reference-data-file", "Read in ref surface integral");
    args.support("nu");
    args.support("stokes");
    args.support("v-max");
    args.support("vtk");
    args.support("use-q2-fbm");
    args.support("face-cub-deg");
    args.support("summed-cub");
    args.support("use-external-fbm-chart");
    args.support("use-fbm-function", "Use an analyic chart for the circle (not recommended)");
    args.support("calc-face-values", "Instead of normal gradients, integrate the actual values of our FE function (which analytical should be zero)");
    args.support("integrator-type", "Integrate surface by integrating either point-value; delta-peak; cell-mean");
    args.support("fbm-chart", "Which chart to use for our fbm assembly");
    args.support("ref-shape-type", "The shape type of the reference face integrals (not used for now)");

    bool face_values = args.check("calc-face-values") >= 0;
    surface_integration_type integrator_type = surface_integration_type::point_values;
    String cubature_face_int = "gauss-legendre:3";
    if(args.check("integrator-type") > 0)
    {
      integrator_type = parse_integration_type(args.query("integrator-type")->second.front());
    }
    String external_fbm_chart_file;
    if(args.check("use-external-fbm-chart") > 0)
    {
      external_fbm_chart_file = args.query("use-external-fbm-chart")->second.front();
    }
    if(args.check("summed-cub")>0)
      cubature_face_int = "refine*" + args.query("summed-cub")->second.front() + ":";
    if(args.check("face-cub-deg")>0)
      cubature_face_int = cubature_face_int + "gauss-legendre:" + args.query("face-cub-deg")->second.front();
    else
      cubature_face_int += "gauss-legendre:3";
    nu = args.parse_default("nu", 0.001);
    navier = !args.parse_default("stokes", false);
    beta = navier ? DataType(1) : DataType(0);
    DataType vmax =args.parse_default("v-max", 0.3);
    bool use_q2_fbm = args.check("use-q2-fbm") >= 0;
    bool use_fbm_function = args.check("use-fbm-function") >= 0;
    logger.print("Running with v-max: " + stringify(vmax), Gendie::info);
    logger.print("Running with nu: " + stringify(nu), Gendie::info);
    logger.print("Running navier: " + stringify(navier), Gendie::info);
    logger.print("Calc face values: " + String(face_values ? "Yes" : "No"), Gendie::info);
    logger.print("Surface Integrator Type: " + print_integration_type(integrator_type), Gendie::info);
    logger.print("Using q2 FBM: " + String(use_q2_fbm ? "Yes" : "No") , Gendie::info);
    logger.print("Setting up domain", Gendie::info);
    PartitionControl<FBMDomainLevel> domain(comm, true);
    {
      auto levels = args.query("level")->second;
      domain.set_desired_levels(levels);
    }

    const String fbm_meshpart_name = "fbm";

    domain.parse_args(args);


    read_mesh_file(domain, comm, args, logger);

    String fbm_chart_name = "fbm";
    if(args.check("fbm-chart") >= 0)
    {
      fbm_chart_name = args.query("fbm-chart")->second.front();
    }
    Tiny::Vector<DataType, dim> mp(DataType(0.2));
    if constexpr (dim == 3)
      mp[0] = DataType(0.5);

    Geometry::SphereHitTestFunction<DataType, dim> fbm_hit_func(mp, 0.05);

    for(Index i = 0; i < domain.size_physical(); ++i)
    {
      auto& cur_dom = *domain.at(i);

      if(use_fbm_function)
      {
        auto* mesh_node = cur_dom.get_mesh_node();
        Geometry::HitTestFactory<decltype(fbm_hit_func), MeshType> hit_factory(fbm_hit_func, *mesh_node->get_mesh());
        mesh_node->add_mesh_part(fbm_meshpart_name, hit_factory.make_unique());
      }
      else
      {
        std::unique_ptr<ChartBaseType> tmp_chart;
        const ChartBaseType* fbm_chart = nullptr;
        if(external_fbm_chart_file.empty())
          fbm_chart = cur_dom.get_mesh_node()->get_atlas()->find_mesh_chart(fbm_chart_name);
        else
        {
          if constexpr(dim == 2)
          {
            typedef Geometry::Atlas::BezierChartParser<MeshType> BezierChartParserType;
            std::ifstream instream;
            instream.open(external_fbm_chart_file);
            XASSERTM(instream.good(), "Error while readin in file from " + external_fbm_chart_file);
            Xml::Scanner scanner(instream);
            std::shared_ptr<BezierChartParserType> bezier_chart_parser = std::make_shared<BezierChartParserType>(tmp_chart);
            scanner.scan(bezier_chart_parser);
            instream.close();

            fbm_chart = tmp_chart.get();
          }
          else
          {
            XABORTM("You shall not arrive here");
          }
        }
        XASSERTM(fbm_chart, "Could not find fbm chart with the name "+ fbm_chart_name);

        cur_dom.get_mesh_node()->add_mesh_part(fbm_meshpart_name, FEAT::Geometry::make_unique_meshpart_by_chart_hit_test(*cur_dom.get_mesh_node()->get_mesh(), *fbm_chart, true, DataType(1E-5)));
      }

      cur_dom.create_fbm_assembler(domain.at(i).layer().comm(), fbm_meshpart_name);
    }


    std::deque<std::shared_ptr<FBMSystemLevel>> system;
    create_systems<decltype(domain), FBMSystemLevel, true>(domain, system);

    // now we need to assemble our filters, first get the names of the inflow and outflow
    if((args.check("inflow") < 1) || (args.check("outflow") < 1))
    {
      comm.print("You have to provide inflow and outflow names");
      Runtime::abort();
    }
    std::deque<String> inflows = args.query("inflow")->second;
    std::deque<String> outflows = args.query("outflow")->second;

    assemble_filters<decltype(domain), FBMSystemLevel, true>(inflows, outflows, system, domain, vmax, use_q2_fbm);

    logger.print("After filter assemble", Gendie::info);

    // system is assembled, now we need to create the components
    auto stokes_solver = std::make_shared<Gendie::MultigridVankaFlowSolver<FBMSystemLevel, true>>(system, domain, nullptr, &logger);
    auto defect_asm = std::make_shared<Gendie::BurgersDefectFlowAssembler<decltype(domain), FBMSystemLevel>>(domain, system.front(), nullptr);
    defect_asm->set_cubature(cubature_matrix_a);
    defect_asm->nu = nu;
    defect_asm->beta = beta;
    auto system_asm = std::make_shared<Gendie::BurgersSystemFlowAssembler<decltype(domain)>>(domain, nullptr);
    system_asm->set_cubature(cubature_matrix_a);
    system_asm->nu = nu;
    system_asm->beta = beta;

    typedef Gendie::MultigridVankaFlowSolver<FBMSystemLevel, true> StS;
    typedef Gendie::BurgersDefectFlowAssembler<decltype(domain), FBMSystemLevel> CaD;
    typedef Gendie::BurgersSystemFlowAssembler<decltype(domain)> CaS;

    auto solver = Gendie::AlPiNeSteadyFlowSolver<StS, CaD, CaS, decltype(system.front()->filter_sys)>(stokes_solver, defect_asm, system_asm, system.front()->filter_sys, &logger);

    solver.min_picard_steps = 2;
    solver._tol_abs = abs_tol;

    auto sol_vec = system.front()->create_global_vector_sys();
    auto rhs_vec = system.front()->create_global_vector_sys();

    sol_vec.format();
    rhs_vec.format();

    system.front()->filter_sys.filter_sol(sol_vec);
    system.front()->filter_sys.filter_rhs(rhs_vec);

    // system.front()->apply_fbm_filter_to_rhs(rhs_vec);

    solver.init();

    solver.init_sol(sol_vec, rhs_vec);
    solver.apply(sol_vec, rhs_vec);

    solver.done();

    if(args.check("vtk") >= 0)
    {
      const MeshType* vtk_mesh = &domain.front()->get_mesh();
      std::unique_ptr<MeshType> refined_mesh;
      Geometry::StandardRefinery<MeshType> refinery(domain.front()->get_mesh());
      refined_mesh = refinery.make_unique();
      vtk_mesh = refined_mesh.get();
      Geometry::ExportVTK<MeshType> vtk_export(*vtk_mesh);
      vtk_export.add_vertex_vector("Vec V", sol_vec.local().template at<0>());
      vtk_export.add_vertex_scalar("fbm_mask", system.front()->fbm_mask_velo.data());
      {
        const auto& filter_vec_velo = system.front()->filter_interface_fbm.get_filter_vector();
        LAFEM::DenseVectorBlocked<DataType, IndexType, dim> filter_vec(filter_vec_velo.size());
        filter_vec.format(-1);
        for(Index k = 0; k < filter_vec_velo.used_elements(); ++k)
        {
          filter_vec(filter_vec_velo.indices()[k], filter_vec_velo.elements()[k]);
        }

        vtk_export.add_vertex_vector("Debug: interface_filter_v", filter_vec);

      }
      {
        const auto& filter_vec_velo = system.front()->get_local_velo_unit_filter_seq().find_or_add("fbm").get_filter_vector();
        LAFEM::DenseVectorBlocked<DataType, IndexType, dim> filter_vec(filter_vec_velo.size());
        filter_vec.format(-1);
        for(Index k = 0; k < filter_vec_velo.used_elements(); ++k)
        {
          filter_vec(filter_vec_velo.indices()[k], filter_vec_velo.elements()[k]);
        }

        vtk_export.add_vertex_vector("Debug: fbm_filter_v", filter_vec);

      }

      String vtk_name = "fbm_test_lvl_" + stringify(domain.size_virtual());
      if(args.check("vtk") > 0)
      {
        vtk_name = args.query("vtk")->second.front();
      }
      vtk_export.write(vtk_name, comm);
    }

    String ref_data_base;
    if(args.check("reference-data-file") > 0)
    {
      ref_data_base = args.query("reference-data-file")->second.front();
    }
    else
    {
      XABORTM("No data file given");
    }

    BinaryStream stream_in;

    // read in through commen mpi data io
    DistFileIO::read_common(stream_in, ref_data_base + "_verts.bin");

    LAFEM::DenseVectorBlocked<DataType, IndexType, dim> verts(FEAT::LAFEM::FileMode::fm_binary, stream_in);


    LAFEM::DenseVectorBlocked<DataType, IndexType, dim> face_val;

    {
      typedef Shape::Hypercube<dim-1> FaceShapeType;
      constexpr int num_verts = Shape::FaceTraits<FaceShapeType, 0>::count;

      FEAT::Assembly::SurfaceIntegrator<TrafoType, FaceShapeType> surf_integrator(domain.front()->trafo, Cubature::DynamicFactory(cubature_face_int));
      std::array<Tiny::Vector<DataType, dim>, num_verts> face_verts;
      for(Index k = 0; k < verts.size(); k += num_verts)
      {
        for(Index l = 0; l < num_verts; ++l)
        {
          face_verts[l] = verts(k+l);
        }
        surf_integrator.add_face_vertices(face_verts);
      }
      surf_integrator.set_mask_vector(domain.front()->fbm_assembler->get_fbm_mask_vector(dim));

      surf_integrator.compile();

      face_val = LAFEM::DenseVectorBlocked<DataType, IndexType, dim>(verts.size()/num_verts);
      face_val.format();

      if(face_values)
      {
        Intern::FaceValueIntegratorJob int_job(face_val, sol_vec.local().template at<0>(), domain.front()->space_velo);

        surf_integrator.assemble(int_job);

      }
      else
      {
        switch (integrator_type)
        {
          case point_values:
          {
            Intern::FaceGradientIntegratorJob int_job(face_val, sol_vec.local().template at<0>(), domain.front()->space_velo);

            surf_integrator.assemble(int_job);
            break;
          }
          case delta_peak:
          {
            FEAT::Space::Lagrange1::Element<TrafoType> mask_space(domain.front()->trafo);
            // for now, use lagrange2 -> lagrange1 ordering to construct our mask vector
            LAFEM::DenseVector<DataType, IndexType> mask_vector(domain.front()->space_velo.get_num_dofs());
            mask_vector.format();
            {
              const auto& loc_fbm_filter_vec = system.front()->get_local_velo_unit_filter_seq().find_or_add("fbm").get_filter_vector();
              const auto& index_array = loc_fbm_filter_vec.indices();
              auto* mask_elem = mask_vector.elements();
              for(Index k = 0; k < loc_fbm_filter_vec.used_elements(); ++k)
              {
                mask_elem[index_array[k]] = DataType(1);
              }
            }
            Intern::FaceGradientFBMDeltaPeakIntegratorJob int_job( face_val, sol_vec.local().template at<0>(), mask_vector, domain.front()->space_velo, mask_space);
            int_job.scatter_alpha = DataType(1);

            surf_integrator.assemble(int_job);
            break;
          }
          case cell_mean:
          {
            Intern::FaceGradientCellMeanIntegratorJob int_job(face_val, sol_vec.local().template at<0>(), domain.front()->space_velo);

            surf_integrator.assemble(int_job);
            break;
          }
          default:
          {
            XABORTM("Unkown Integrator Case");
            break;
          }
        }
      }
    }

    LAFEM::DenseVectorBlocked<DataType, IndexType, dim> all_face_val(comm.rank()==0 ? face_val.size() : 0);
    all_face_val.format();

    {
      comm.reduce(face_val.template elements<FEAT::LAFEM::Perspective::pod>(), all_face_val.template elements<FEAT::LAFEM::Perspective::pod>(),
                   face_val.template size<FEAT::LAFEM::Perspective::pod>(), Dist::op_sum, 0);
    }

    if(comm.rank() == 0)
    {
      LAFEM::DenseVectorBlocked<DataType, IndexType, dim> ref_vals(FEAT::LAFEM::FileMode::fm_binary, ref_data_base + "_values.bin");
      // std::cout << all_face_val << "\n--------------------------------------\n";

      // std::cout << ref_vals << "\n--------------------------------------\n";
      LAFEM::DenseVectorBlocked<DataType, IndexType, dim> diff(ref_vals.clone(LAFEM::CloneMode::Weak));

      diff.axpy(all_face_val, DataType(-1));

      auto err = DataType(1) / Math::sqrt(DataType(diff.size())) * diff.norm2_blocked();

      std::cout << "Error is " << err << "\n-----------------------\n";

      std::cout << "Ref: " << ref_vals << "\n-----------------------\n";
      std::cout << "Calc: " << all_face_val << "\n-------------------\n";
      std::cout << "Diff: " << diff << "\n--------------------------\n";
    }






  }

  void run_no_fbm(Dist::Comm& comm, SimpleArgParser& args, Gendie::Logger& logger)
  {
    logger.print("Running no fbm case", Gendie::info);
    args.support("create-reference-data", "Create reference data and write it out as blocked vector");
    args.support("nu");
    args.support("stokes");
    args.support("v-max");
    nu = args.parse_default("nu", 0.001);
    navier = !args.parse_default("stokes", false);
    beta = navier ? DataType(1) : DataType(0);
    DataType vmax =args.parse_default("v-max", 0.3);
    args.support("calc-face-values", "Instead of normal gradients, integrate the actual values of our FE function (which are zero due to boundary conditions)");
    bool face_values = args.check("calc-face-values") >= 0;
    logger.print("Running with v-max: " + stringify(vmax), Gendie::info);
    logger.print("Running with nu: " + stringify(nu), Gendie::info);
    logger.print("Running navier: " + stringify(navier), Gendie::info);
    logger.print("Calc face values: " + String(face_values ? "Yes" : "No"), Gendie::info);
    logger.print("Setting up domain", Gendie::info);

    args.support("vtk");

    args.support("export-chart");
    const bool export_bezier_chart = args.check("export-chart") >= 0;
    // first of all, create domain
    PartitionControl<DomainLevel> domain(comm, true);
    {
      auto levels = args.query("level")->second;
      domain.set_desired_levels(levels);
      if(export_bezier_chart)
        domain.keep_base_levels();
    }

    domain.parse_args(args);


    read_mesh_file(domain, comm, args, logger);


    std::deque<std::shared_ptr<SystemLevel>> system;
    create_systems<decltype(domain), SystemLevel, false>(domain, system);

    // now we need to assemble our filters, first get the names of the inflow and outflow
    if((args.check("inflow") < 1) || (args.check("outflow") < 1))
    {
      comm.print("You have to provide inflow and outflow names");
      Runtime::abort();
    }
    std::deque<String> inflows = args.query("inflow")->second;
    std::deque<String> outflows = args.query("outflow")->second;

    assemble_filters<decltype(domain), SystemLevel, false>(inflows, outflows, system, domain, vmax);

    logger.print("After filter assemble", Gendie::info);
    // system is assembled, now we need to create the components
    auto stokes_solver = std::make_shared<Gendie::MultigridVankaFlowSolver<SystemLevel, false>>(system, domain, nullptr, &logger);
    auto defect_asm = std::make_shared<Gendie::BurgersDefectFlowAssembler<decltype(domain), SystemLevel>>(domain, system.front(), nullptr);
    defect_asm->set_cubature(cubature_matrix_a);
    defect_asm->nu = nu;
    defect_asm->beta = beta;
    auto system_asm = std::make_shared<Gendie::BurgersSystemFlowAssembler<decltype(domain)>>(domain, nullptr);
    system_asm->set_cubature(cubature_matrix_a);
    system_asm->nu = nu;
    system_asm->beta = beta;

    typedef Gendie::MultigridVankaFlowSolver<SystemLevel, false> StS;
    typedef Gendie::BurgersDefectFlowAssembler<decltype(domain), SystemLevel> CaD;
    typedef Gendie::BurgersSystemFlowAssembler<decltype(domain)> CaS;

    auto solver = Gendie::AlPiNeSteadyFlowSolver<StS, CaD, CaS, decltype(system.front()->filter_sys)>(stokes_solver, defect_asm, system_asm, system.front()->filter_sys, &logger);

    solver.min_picard_steps = 2;
    solver._tol_abs = abs_tol;

    auto sol_vec = system.front()->create_global_vector_sys();
    auto rhs_vec = system.front()->create_global_vector_sys();

    sol_vec.format();
    rhs_vec.format();

    system.front()->filter_sys.filter_sol(sol_vec);
    system.front()->filter_sys.filter_rhs(rhs_vec);

    solver.init();

    solver.init_sol(sol_vec, rhs_vec);
    solver.apply(sol_vec, rhs_vec);

    solver.done();

    if(args.check("vtk")>=0)
    {
      Geometry::ExportVTK<MeshType> vtk_export(domain.front()->get_mesh());
      vtk_export.add_vertex_vector("Vec V", sol_vec.local().template at<0>());
      String vtk_name = "no_fbm_test_lvl_" + stringify(domain.size_virtual());
      if(args.check("vtk") > 0)
      {
        vtk_name = args.query("vtk")->second.front();
      }
      vtk_export.write(vtk_name, comm);
    }

  //  std::cout << "Vec sol " << sol_vec.local().template at<0>() << "\n ---------------------------------- \n";


    bool create_ref_data = false;
    String ref_file = "ref_data";

    create_ref_data = args.check("create-reference-data") >= 0;
    if(args.check("create-reference-data") >= 1)
    {
      ref_file = args.query("create-reference-data")->second.front();
    }

    if(!create_ref_data)
    {
      comm.print("Do not create ref data!");
      return;
    }

    String mesh_part_name = "bnd:c";
    if(args.check("ref-part-name") >= 0)
    {
      mesh_part_name = args.query("ref-part-name")->second.front();
    }

    // at this point, since we hav to guarantee positive oriented normals, use the normal flipper to fix our target set holder
    FEAT::Geometry::FacetFlipper<ShapeType>::reorient(domain.front()->get_mesh().get_index_set_holder());


    // create ref data, so first, get local mesh parts
    const auto* ref_mesh_part_node = domain.front()->get_mesh_node()->find_mesh_part_node(mesh_part_name);
    XASSERTM(ref_mesh_part_node, "Could not find meshpart " + mesh_part_name);

    FEAT::LAFEM::DenseVectorBlocked<DataType, IndexType, dim> part_vertices;
    const auto& face_to_vert = domain.front()->get_mesh().get_index_set<dim-1, 0>();

    FEAT::LAFEM::DenseVectorBlocked<DataType, IndexType, dim> vec_face;

    if(const auto* mesh_part = ref_mesh_part_node->get_mesh())
    {
      const auto& target_set = mesh_part->get_target_set<dim-1>();
      // construct our vertices for the local face
      // get our mapping face to vertices
      const auto& vertex_set = domain.front()->get_mesh().get_vertex_set();
      part_vertices = FEAT::LAFEM::DenseVectorBlocked<DataType, IndexType, dim>(target_set.get_num_entities()*face_to_vert.get_num_indices());
      for(Index k = 0; k < target_set.get_num_entities(); ++k)
      {
        const Index cur_face = target_set[k];

        const auto& loc_face_set = face_to_vert[cur_face];
        for(Index l = 0; l < Index(face_to_vert.get_num_indices()); ++l)
        {
          part_vertices(k*Index(face_to_vert.get_num_indices()) + l, vertex_set[loc_face_set[int(l)]]);
        }
      }

      // now, assemble our local face vector
      FEAT::Assembly::TraceAssembler trace_asm(domain.front()->trafo);
      trace_asm.add_mesh_part(*mesh_part);

      trace_asm.compile();

      FEAT::LAFEM::DenseVectorBlocked<DataType, IndexType, dim> vec_face_whole_mesh(domain.front()->get_mesh().get_num_entities(dim-1));
      vec_face_whole_mesh.format();

      if(face_values)
      {
        Intern::FaceValueTraceAssemblyJob trace_asm_job(vec_face_whole_mesh, sol_vec.local().template at<0>(), domain.front()->space_velo, cubature_matrix_a);
        trace_asm.assemble(trace_asm_job);

      }
      else
      {
        Intern::FaceGradientTraceAssemblyJob trace_asm_job(vec_face_whole_mesh, sol_vec.local().template at<0>(), domain.front()->space_velo, cubature_matrix_a);
        trace_asm.assemble(trace_asm_job);
      }

      vec_face = LAFEM::DenseVectorBlocked<DataType, IndexType, dim>(target_set.get_num_entities());
      vec_face.format();
      for(Index k = 0; k < target_set.get_num_entities(); ++k)
      {
        vec_face(k, vec_face_whole_mesh(target_set[k]));
      }
    }

    // communicate overall number of faces to export
    std::vector<Index> num_faces_per_rank(comm.size());
    {
      Index num_faces = vec_face.size();
      comm.gather(&num_faces, 1, num_faces_per_rank.data(), 1, 0);
    }

    Index num_faces_global = comm.rank() == 0 ? std::accumulate(num_faces_per_rank.begin(), num_faces_per_rank.end(), Index(0)) : Index(0);
    // will be filled on rank 0
    FEAT::LAFEM::DenseVectorBlocked<DataType, IndexType, dim> global_vec_face(num_faces_global);
    FEAT::LAFEM::DenseVectorBlocked<DataType, IndexType, dim> global_part_vertices(num_faces_global*face_to_vert.get_num_indices());


    std::vector<Index> num_faces_offsets(comm.size());
    std::exclusive_scan(num_faces_per_rank.begin(), num_faces_per_rank.end(), num_faces_offsets.begin(), Index(0));

    // now create gather and sendbuffer and get our values
    {
      std::vector<std::vector<DataType>> recv_buffer(comm.rank() == 0 ? comm.size() : 0);
      if(comm.rank() == 0)
      {
        for(Index k = 0; k < num_faces_per_rank.size(); ++k)
        {
          recv_buffer.at(k).resize(num_faces_per_rank[k]*dim*(1+face_to_vert.get_num_indices()));
        }
      }
      // reserve enough space for face values and face vertices
      std::vector<DataType> vtx_send_buffer(vec_face.size()*dim*(1+face_to_vert.get_num_indices()));

      // fill our buffer array, first face_values, then our vertices
      const DataType* vtx_val = part_vertices.template elements<LAFEM::Perspective::pod>();
      const DataType* face_val = vec_face.template elements<LAFEM::Perspective::pod>();
      DataType* gvtx_val = global_part_vertices.template elements<LAFEM::Perspective::pod>();
      DataType* gface_val = global_vec_face.template elements<LAFEM::Perspective::pod>();
      for(Index k = 0; k < vec_face.size()*dim; ++k)
      {
        vtx_send_buffer[k] = face_val[k];
      }
      for(Index k = 0; k < vec_face.size()*dim*(face_to_vert.get_num_indices()); ++k)
      {
        vtx_send_buffer[k+vec_face.size()*dim] = vtx_val[k];
      }

      Dist::RequestVector requests;
      if(comm.rank() == 0)
      {
        for(Index k = 1; k < num_faces_per_rank.size(); ++k)
        {
          if(num_faces_per_rank[k] > 0)
          {
            requests.push_back(comm.irecv(recv_buffer.at(k).data(), num_faces_per_rank[k]*dim*(1+face_to_vert.get_num_indices()), int(k)));
          }
        }

        // copy this rank
        {
          auto rcomm_rank = 0;
          for(Index k = 0; k < num_faces_per_rank.at(rcomm_rank)*dim; ++k)
          {
            gface_val[k] = face_val[k];
          }
          for(Index k = 0; k < num_faces_per_rank.at(rcomm_rank)*dim*(face_to_vert.get_num_indices()); ++k)
          {
            gvtx_val[k] = vtx_val[k];
          }

        }

        for(std::size_t idx(0u); requests.wait_any(idx); )
        {
          #ifdef FEAT_HAVE_MPI
          Index rcomm_rank = Index(requests.get_status(idx).source());
          #else
          Index rcomm_rank = Index(0);
          #endif
          Index cur_offset = num_faces_offsets[rcomm_rank]*dim;
          for(Index k = 0; k < num_faces_per_rank.at(rcomm_rank)*dim; ++k)
          {
            gface_val[cur_offset + k] = recv_buffer.at(rcomm_rank)[k];
          }
          cur_offset = num_faces_offsets[rcomm_rank]*dim*(face_to_vert.get_num_indices());
          Index buffer_offset = num_faces_per_rank.at(rcomm_rank)*dim;
          for(Index k = 0; k < num_faces_per_rank.at(rcomm_rank)*dim*(face_to_vert.get_num_indices()); ++k)
          {
            gvtx_val[cur_offset+k] = recv_buffer.at(rcomm_rank)[k+buffer_offset];
          }
        }
      }
      else if(vec_face.size() > 0)
      {
        comm.send(vtx_send_buffer.data(), vtx_send_buffer.size(), 0);
      }

    }

    // write out our vectors
    logger.print(String("Writing out reference data to ") + ref_file, Gendie::info);
    if(comm.rank() == 0)
    {
      global_part_vertices.write_out(FEAT::LAFEM::FileMode::fm_binary, ref_file + String("_verts.bin"));
      global_vec_face.write_out(FEAT::LAFEM::FileMode::fm_binary, ref_file + String("_values.bin"));
    }

    if constexpr(dim == 2)
    {
      if(comm.rank() == 0 && export_bezier_chart)
      {
        // construct our mesh on rank 0 on our specified resolution
        const auto& base_mesh_node = comm.size() == 1 ? domain.front()->get_mesh_node() : domain.front().level_b().get_mesh_node();

        if(const auto* meshpart = base_mesh_node->find_mesh_part_node(mesh_part_name)->get_mesh())
        {
          String bezier_name = "bezier_" + mesh_part_name;
          if(args.check("export-chart") > 0)
          {
            bezier_name = args.query("export-chart")->second.front();
          }
          std::cout << "Writing bezier chart to " << bezier_name << ".xml\n";

          std::stringstream stream;
          const MeshType* mesh = comm.size() == 1 ? &domain.front()->get_mesh() : &domain.front().level_b().get_mesh();
          auto bezier_chart = Geometry::Intern::MeshPartToChartExporter<ShapeType>::meshpart_to_bezier(*mesh, *meshpart);
          bezier_chart->write(stream, String(" "));

          {
            std::ofstream file;
            file.open(bezier_name + ".xml", std::ios_base::out | std::ios_base::trunc);
            XASSERTM(file.good(), "Something went wrong while writing the file");
            file << stream.rdbuf();
            file.close();
          }
        }
        else
        {
          std::cerr << "Could not find meshpart on combined mesh\n";
        }
      }
    }


  }


  void main(int argc, char** argv)
  {
    Dist::Comm comm(Dist::Comm::world());
    SimpleArgParser args(argc, argv);
    Gendie::Logger logger(comm);

    // cast typenames to ignore warnings
    (void)Gendie::fp_typename;
    (void)Gendie::ix_typename;

    args.support("fbm-mode");
    args.support("mesh");
    args.support("inflow-name");
    args.support("outflow-name");
    args.support("inflow-velo");
    args.support("level");

    if(args.check("level") <= 0)
    {
      comm.print("You have to provide a level string");
      Runtime::abort();
    }

    if(args.query("fbm-mode"))
    {
      run_fbm(comm, args, logger);
    }
    else
    {
      run_no_fbm(comm, args, logger);
    }
  }

}

int main(int argc, char** argv)
{
  FEAT::Runtime::ScopeGuard guard(argc, argv);
  FBMTest::main(argc, argv);
  return 0;
}
