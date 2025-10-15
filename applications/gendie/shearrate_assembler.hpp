#pragma once

#include <kernel/assembly/domain_assembler.hpp>
#include <kernel/assembly/domain_assembler_basic_jobs.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <kernel/util/math.hpp>


namespace Gendie
{
  template<typename Vector_, typename PrimalVector_, typename Space_>
  class JacobianVectorAssemblyJob
  {
  public:
    typedef typename Vector_::DataType DataType;
    typedef typename Vector_::ValueType ValueType;
    typedef typename PrimalVector_::ValueType PrimalValueType;

    typedef Space_ SpaceType;
    typedef typename SpaceType::TrafoType TrafoType;

    static constexpr int dim_ = SpaceType::shape_dim;


    class Task
    {
    public:
      /// this task doesn't need to scatter
      static constexpr bool need_scatter = true;
      /// this task needs to combine
      static constexpr bool need_combine = false;

    protected:
      static constexpr FEAT::SpaceTags space_tags = FEAT::SpaceTags::value | FEAT::SpaceTags::grad;

      /// our assembly traits
      typedef FEAT::Assembly::AsmTraits1<DataType, Space_, FEAT::TrafoTags::jac_det, space_tags> AsmTraits;
      /// the vector that is to be scattered into
      Vector_& vector;
      /// the vector that is to be integrated
      const PrimalVector_& primal_vector;
      /// the finite element space to be used
      const Space_& space;
      /// the cubature factory used for integration
      const typename AsmTraits::TrafoType& trafo;
      /// the trafo evaluator
      typename AsmTraits::TrafoEvaluator trafo_eval;
      /// the space evaluator
      typename AsmTraits::SpaceEvaluator space_eval;
      /// the space dof-mapping
      typename AsmTraits::DofMapping dof_mapping;
      /// the cubature rule used for integration
      typename AsmTraits::CubatureRuleType cubature_rule;
      /// the trafo evaluation data
      typename AsmTraits::TrafoEvalData trafo_data;
      /// the space evaluation data
      typename AsmTraits::SpaceEvalData space_data;
      /// the local vector to be assembled
      typename AsmTraits::template TLocalVector<ValueType> local_vector;
      /// the local vector to be assembled
      typename AsmTraits::template TLocalVector<PrimalValueType> local_primal_vector;
      /// The local gradient type
      FEAT::Tiny::Matrix<DataType, dim_, dim_> grad_v;
      /// the vector gather object
      typename PrimalVector_::GatherAxpy gather_axpy;
      /// the vector scatter object
      typename Vector_::ScatterAxpy scatter_axpy;
      /// scatter scaling factor
      DataType alpha;
      /// maximum number of local dofs
      static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

    public:
      explicit Task(JacobianVectorAssemblyJob& job) :
        vector(job._vector),
        primal_vector(job._primal_vector),
        space(job._space),
        trafo(space.get_trafo()),
        trafo_eval(trafo),
        space_eval(space),
        dof_mapping(space),
        cubature_rule(FEAT::Cubature::ctor_factory, job._cubature_factory),
        trafo_data(),
        space_data(),
        local_vector(),
        local_primal_vector(),
        grad_v(),
        gather_axpy(primal_vector),
        scatter_axpy(vector),
        alpha(job._alpha)
      {
      }

      void prepare(FEAT::Index cell)
      {
        // prepare dof mapping
        dof_mapping.prepare(cell);

        // prepare trafo evaluator
        trafo_eval.prepare(cell);

        // prepare space evaluator
        space_eval.prepare(trafo_eval);

        // format vector to be gathered
        local_primal_vector.format();

        // gather local vector data
        gather_axpy(local_primal_vector, dof_mapping);
      }

      void assemble()
      {
        // format local vector
        local_vector.format();

        // fetch number of local dofs
        const int num_loc_dofs = space_eval.get_num_local_dofs();

        // loop over all quadrature points and integrate
        for(int k(0); k < cubature_rule.get_num_points(); ++k)
        {
          // compute trafo data
          trafo_eval(trafo_data, cubature_rule.get_point(k));

          // compute basis function data
          space_eval(space_data, trafo_data);

          // construct gradient
          grad_v.format();
          for(int i = 0; i < num_loc_dofs; ++i)
          {
            grad_v.add_outer_product(local_primal_vector[i], space_data.phi[i].grad);
          }

          const DataType weight = cubature_rule.get_weight(k) * trafo_data.jac_det;

          // and now integrate with test function
          for(int i = 0;  i < num_loc_dofs; ++i)
          {
            const DataType value = weight * space_data.phi[i].value;
            for(int m = 0; m < dim_; ++m)
            {
              for(int n = 0; n < dim_; ++n)
              {
                FEAT::Tiny::axpy(local_vector[i][m*dim_ + n], grad_v[m][n], value);
              }
            }
          }

        }
      }

      void scatter()
      {
        // nothing to do here
        scatter_axpy(local_vector, dof_mapping, alpha);
      }

      void finish()
      {
        space_eval.finish();
        trafo_eval.finish();
      }

      void combine()
      {
      }
    }; // class Task

    protected:
      /// the vector that is to be assembled
      Vector_& _vector;
      /// the vector that is to be assembled
      const PrimalVector_& _primal_vector;
      /// the finite element space to be used
      const Space_& _space;
      /// the cubature factory
      FEAT::Cubature::DynamicFactory _cubature_factory;
      /// scatter alpha
      DataType _alpha;

    public:
      explicit JacobianVectorAssemblyJob(Vector_& vector, const PrimalVector_& primal_vector, const Space_& space,
        const FEAT::String& cubature, DataType alpha = DataType(1))
        : _vector(vector),
          _primal_vector(primal_vector),
          _space(space),
          _cubature_factory(cubature),
          _alpha(alpha)
      {
      }

  };

  template<typename ShearVector_, typename VeloVector_, typename MassVector_, typename Space_>
  static inline void assemble_shear_vector(FEAT::Assembly::DomainAssembler<typename Space_::TrafoType>& dom_asm, ShearVector_& shear_vec, const VeloVector_& velo_vector, const MassVector_& lumped_inv_mass,
    const Space_& space, const FEAT::String& cubature, typename ShearVector_::DataType min_shear = typename ShearVector_::DataType(1E-2),
    typename ShearVector_::DataType max_shear = typename ShearVector_::DataType(1E+6), bool clamp_shear = false)
  {
    XASSERTM(space.get_num_dofs() == shear_vec.size(), "Shear vector and space do not fit");
    XASSERTM(space.get_num_dofs() == velo_vector.local().size(), "Velocity vector and space do not fit");
    XASSERTM(space.get_num_dofs() == lumped_inv_mass.size(), "Lumped mass vector and space do not fit");

    // first create temporary gradient vector
    constexpr int dim_ = Space_::shape_dim;
    static_assert(dim_ > 1);
    typedef typename ShearVector_::DataType DataType;
    // we require a gate for our "grad" vector
    typedef FEAT::LAFEM::DenseVectorBlocked<typename ShearVector_::DataType, typename ShearVector_::IndexType, dim_*dim_> GradVectorType;
    FEAT::Global::Gate<GradVectorType, typename VeloVector_::MirrorType> grad_gate;
    {
      GradVectorType freq_vec(space.get_num_dofs());
      freq_vec.format();
      grad_gate.convert(*velo_vector.get_gate(), std::move(freq_vec));
    }
    // create global grad vector
    FEAT::Global::Vector<GradVectorType, typename VeloVector_::MirrorType> grad_vector(&grad_gate, space.get_num_dofs());
    grad_vector.format();

    // assemble dual gradient values
    JacobianVectorAssemblyJob grad_asm_job(grad_vector.local(), velo_vector.local(), space, cubature);
    dom_asm.assemble(grad_asm_job);

    // sync vector -> type 0 summed values
    grad_vector.sync_0();

    auto* loc_val_shear = shear_vec.elements();
    const auto* loc_val_grad = grad_vector.local().elements();
    const auto* loc_mass_inv = lumped_inv_mass.elements();

    // and now apply inverted mass matrix, calculate Du and from there the (clamped) shearrate
    FEAT_PRAGMA_OMP(parallel for)
    for(FEAT::Index i = 0; i < space.get_num_dofs(); ++i)
    {
      FEAT::Tiny::Matrix<DataType, dim_, dim_> loc_grad;
      for(int k = 0; k < dim_; ++k)
      {
        for(int l = 0; l < dim_; ++l)
        {
          loc_grad[k][l] = loc_val_grad[i][k*dim_ + l]*loc_mass_inv[i][k];
        }
      }
      FEAT::Tiny::Matrix<DataType, dim_, dim_> loc_grad_t;
      loc_grad_t.set_transpose(loc_grad);
      loc_grad_t.axpy(DataType(1), loc_grad);
      loc_val_shear[i] = FEAT::Math::sqrt(0.5) * loc_grad_t.norm_frobenius();
      if(clamp_shear)
      {
        loc_val_shear[i] = FEAT::Math::clamp(loc_val_shear[i], max_shear, min_shear);
      }
    }
  }

}