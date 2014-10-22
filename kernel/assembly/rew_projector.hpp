#pragma once
#ifndef KERNEL_ASSEMBLY_REW_PROJECTOR_HPP
#define KERNEL_ASSEMBLY_REW_PROJECTOR_HPP 1

// includes, FEAST
#include <kernel/assembly/asm_traits.hpp>
#include <kernel/util/linear_algebra.hpp>

namespace FEAST
{
  namespace Assembly
  {
    /**
     * \brief Restricted Element-Wise Projection operator
     *
     * This class template implements the "restricted element-wise projection" operator
     * which projects an (analytic) function into a Finite-Element space by performing
     * weighted element-wise L2-projections.
     *
     * \author Peter Zajac
     */
    class RewProjector
    {
    private:
      /// \cond internal
      struct AsmSpaceConfig :
        public Space::ConfigBase
      {
        enum
        {
          need_value = 1
        };
      };
      struct AsmTrafoConfig :
        public Trafo::ConfigBase
      {
        enum
        {
          need_img_point = 1,
          need_jac_det = 1
        };
      };
      /// \endcond

    public:
      /**
       * \brief Weighting type enumeration
       */
      enum WeightType
      {
        /// use arithmetic averaging
        wt_arithmetic = 0,
        /// use volume-based averaging
        wt_volume
      };

      /**
       * \brief Projects an analytic function into a finite element space.
       *
       * \param[out] vector
       * A reference to the coefficient vector that is to be assembled.
       *
       * \param[in] function
       * An object implementing the AnalyticFunction interface capable of computing function values.
       *
       * \param[in] space
       * A reference to the space to which to project into.
       *
       * \param[in] cubature_factory
       * A cubature factory to be used for integration.
       *
       * \param[in] weight_type
       * The weighting type to be used.
       */
      template<
        typename Vector_,
        typename Function_,
        typename Space_,
        typename CubatureFactory_>
      static void project(
        Vector_& vector,
        const Function_& function,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        WeightType weight_type = wt_volume)
      {
        typedef Vector_ VectorType;
        typedef Function_ FunctionType;
        typedef Space_ SpaceType;
        typedef typename SpaceType::TrafoType TrafoType;

        static_assert(Function_::can_value != 0, "analytic function can't compute function values");

        // define assembly traits
        typedef AsmTraits1<typename Vector_::DataType, SpaceType, AsmTrafoConfig, AsmSpaceConfig> AsmTraits;
        typedef typename AsmTraits::DataType DataType;

        // define the cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // fetch the trafo and the mesh
        const TrafoType& trafo(space.get_trafo());

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create a functor evaluator
        typename FunctionType::template Evaluator<typename AsmTraits::AnalyticEvalTraits> func_eval(function);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // allocate fine-mesh mass matrix
        typename AsmTraits::LocalMatrixType mass;

        // pivot array for factorisation
        Index pivot[mass.n];

        // create local vector data
        typename AsmTraits::LocalVectorType lvad;

        // create local weight data
        typename AsmTraits::LocalVectorType lwad;

        // fetch the dof count
        const Index num_dofs(space.get_num_dofs());

        // create a clear output vector
        vector = Vector_(num_dofs, DataType(0));

        // allocate weight vector
        Vector_ weight(num_dofs, DataType(0));

        // create vector scatter-axpy
        LAFEM::ScatterAxpy<VectorType> scatter_axpy(vector);

        // create weight scatter-axpy
        LAFEM::ScatterAxpy<VectorType> weight_scatter_axpy(weight);

        // format local weight vector
        if(weight_type == wt_arithmetic)
          lwad.format(DataType(1));

        // loop over all cells of the mesh
        for(Index cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // prepare functor evaluator
          func_eval.prepare(trafo_eval);

          // fetch number of local dofs
          Index num_loc_dofs(space_eval.get_num_local_dofs());

          // format local system
          lvad.format();
          mass.format();

          // cell volume
          DataType volume(DataType(0));

          // loop over all cubature points and integrate
          for(Index k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(k));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // compute function value
            DataType fv(func_eval.value(trafo_data));

            // pre-compute integration weight
            DataType omega(trafo_data.jac_det * cubature_rule.get_weight(k));

            // update volume
            volume += omega;

            // test function loop
            for(Index i(0); i < num_loc_dofs; ++i)
            {
              // assemble RHS entry
              lvad(i) += omega * fv * space_data.phi[i].value;

              // trial function loop
              for(Index j(0); j < num_loc_dofs; ++j)
              {
                mass(i,j) += omega * space_data.phi[i].value * space_data.phi[j].value;
              }
            }
          }

          // finish functor evaluator
          func_eval.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();

          // try to factorise the local mass matrix
          LinAlg::mat_factorise(num_loc_dofs, num_loc_dofs, Index(mass.n), &mass.v[0][0], pivot);

          // solve M*x=f
          LinAlg::mat_solve_vec<false>(num_loc_dofs, &lvad.v[0], Index(mass.n), &mass.v[0][0], pivot);

          // choose weight
          if(weight_type == wt_volume)
          {
            lvad *= volume;
            lwad.format(volume);
          }

          // prepare dof-mapping
          dof_mapping.prepare(cell);

          // incorporate local vector and weights
          scatter_axpy(lvad, dof_mapping);
          weight_scatter_axpy(lwad, dof_mapping);

          // finish dof-mapping
          dof_mapping.finish();
        }

        // finally, scale the vector by the weights
        DataType* vx(vector.elements());
        DataType* wx(weight.elements());
        for(Index i(0); i < num_dofs; ++i)
        {
          if(wx[i] > DataType(0))
            vx[i] /= wx[i];
        }

        // okay
      }
    }; // class RewProjector
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_REW_PROJECTOR_HPP
