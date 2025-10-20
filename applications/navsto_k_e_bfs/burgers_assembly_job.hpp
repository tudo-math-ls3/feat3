// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/assembly/base.hpp>
#include <kernel/assembly/asm_traits.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/global/vector.hpp>

namespace Turb
{
  template<typename MatrixA_, typename VeloVector_, typename TurbVector_, typename SpaceVelo_>
  class TurbBurgersAssemblyJob
  {
  public:
    typedef MatrixA_ MatrixTypeA; // type of matrix_a
    typedef VeloVector_ VeloVectorType; // type of velocity vector component
    typedef TurbVector_ TurbVectorType; // type of turb vector component
    typedef SpaceVelo_ SpaceVeloType;

    typedef typename MatrixTypeA::DataType DataType;

    static constexpr int dim = SpaceVeloType::shape_dim; // system dimension

    /// inout: pointer to A system matrix to be assembled
    MatrixTypeA* matrix_a = nullptr;

    /// inout: pointer to k def vector to be assembled
    VeloVectorType* vec_def_v = nullptr;

    /// input: pointer to last v solution vector
    const VeloVectorType* vec_sol_v = nullptr;

    /// input: pointer to last k solution vector
    const TurbVectorType* vec_sol_k = nullptr;

    /// input: pointer to last k solution vector
    const TurbVectorType* vec_sol_e = nullptr;

    /// input: velocity space type
    const SpaceVeloType* space_velo = nullptr;

    /// input: use deformation tensor?
    bool deformation = true;

    /// enable k-e viscosity?
    bool enable_k_e = false;

    /// input: viscosity parameter nu
    DataType nu = DataType(1);

    /// input: constant C_mu
    DataType C_mu = DataType(0);

    /// input: convection term
    DataType beta = DataType(1);

    /// input: convection term Frechet derivative for matrix
    DataType beta_frechet = DataType(1);

    /// input: reaction term theta for matrix a
    DataType theta = DataType(0);

  public:
    class Task
    {
    public:
      /// this task needs to scatter
      static constexpr bool need_scatter = true;
      /// this task doesn't need to combine
      static constexpr bool need_combine = false;

      /// our assembly traits
      typedef Assembly::AsmTraits1<DataType, SpaceVeloType, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

      /// a reference to our job object
      const TurbBurgersAssemblyJob& job;

      /// the test-/trial-space to be used
      const SpaceVeloType& space_velo;
      /// the cubature factory used for integration
      const typename AsmTraits::TrafoType& trafo;
      /// the trafo evaluator
      typename AsmTraits::TrafoEvaluator trafo_eval;
      /// the space evaluator
      typename AsmTraits::TestEvaluator space_eval_velo;
      /// the space dof-mapping
      typename AsmTraits::TestDofMapping dof_mapping_velo;
      /// the cubature rule used for integration
      typename AsmTraits::CubatureRuleType cubature_rule;
      /// the trafo evaluation data
      typename AsmTraits::TrafoEvalData trafo_data;
      /// the space evaluation data
      typename AsmTraits::TestEvalData space_data_velo;
      /// convection vector gather-axpy object
      typename VeloVectorType::GatherAxpy gather_vec_v;
      typename TurbVectorType::GatherAxpy gather_vec_k;
      typename TurbVectorType::GatherAxpy gather_vec_e;

      /// matrix scatter objects
      typedef typename MatrixTypeA::ScatterAxpy ScatterMatrixA;
      std::unique_ptr<ScatterMatrixA> scatter_matrix_a;
      /// vector scatter objects
      typedef typename VeloVectorType::ScatterAxpy ScatterVeloVector;
      std::unique_ptr<ScatterVeloVector> scatter_vec_def;

      /// maximum number of local dofs
      static constexpr int max_local_dofs_velo = AsmTraits::max_local_test_dofs;

      /// actual number of local dofs on current element
      int num_local_dofs_velo;

      /// local convection field dofs
      Tiny::Vector<Tiny::Vector<DataType, dim>, max_local_dofs_velo> local_dofs_v;

      /// local turbulence dofs
      Tiny::Vector<DataType, max_local_dofs_velo> local_dofs_k, local_dofs_e;

      /// local convection field value
      Tiny::Vector<DataType, dim> value_v;

      // our local velocity gradient
      Tiny::Matrix<DataType, dim, dim> grad_v;

      /// local k and e values
      DataType value_k, value_e;

      /// the local k and e matrices to be assembled
      typedef Tiny::Matrix<DataType, dim, dim> MatrixValue;
      Tiny::Matrix<MatrixValue, max_local_dofs_velo, max_local_dofs_velo> local_matrix_a;

      /// the local k and e rhs vectors to be assembled
      typedef Tiny::Vector<DataType, dim> VectorValue;
      Tiny::Vector<VectorValue, max_local_dofs_velo> local_vec_def;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] job_
       * A \resident reference to the job containing the assembly parameters.
       */
      explicit Task(const TurbBurgersAssemblyJob& job_) :
        job(job_),
        space_velo(*job_.space_velo),
        trafo(space_velo.get_trafo()),
        trafo_eval(trafo),
        space_eval_velo(space_velo),
        dof_mapping_velo(space_velo),
        gather_vec_v(*job.vec_sol_v),
        gather_vec_k(*job.vec_sol_k),
        gather_vec_e(*job.vec_sol_e)

      {
        Cubature::DynamicFactory cubature_factory("gauss-legendre:3");
        cubature_factory.create(cubature_rule);
        if(job.matrix_a)
          scatter_matrix_a = std::make_unique<ScatterMatrixA>(*job.matrix_a);
        if(job.vec_def_v)
          scatter_vec_def = std::make_unique<ScatterVeloVector>(*job.vec_def_v);
      }

      /**
       * \brief Prepares the assembly task for a cell/element.
       *
       * \param[in] cell
       * The index of the element/cell that is to be assembled on
       */
      void prepare(const Index cell)
      {
        // prepare dof mapping
        dof_mapping_velo.prepare(cell);

        // prepare trafo evaluator
        trafo_eval.prepare(cell);

        // prepare space evaluator
        space_eval_velo.prepare(trafo_eval);

        // fetch number of local dofs
        num_local_dofs_velo = space_eval_velo.get_num_local_dofs();

        // gather our local dofs
        local_dofs_v.format();
        local_dofs_k.format();
        local_dofs_e.format();
        gather_vec_v(local_dofs_v, dof_mapping_velo);
        gather_vec_k(local_dofs_k, dof_mapping_velo);
        gather_vec_e(local_dofs_e, dof_mapping_velo);
      }

      /**
       * \brief Prepares the task for a cubature point
       *
       * This function evaluates the trafo and the space in the current cubature point and
       * also evaluates the local convection field value and gradient as well as the streamline
       * diffusion coefficients, if required.
       *
       * \param[in] cubature_point
       * The index of the current cubature point
       */
      void prepare_point(int cubature_point)
      {
        // compute trafo data
        trafo_eval(trafo_data, cubature_rule.get_point(cubature_point));

        // compute basis function data
        space_eval_velo(space_data_velo, trafo_data);

        // compute local v
        value_v.format();
        grad_v.format();
        for(int i(0); i < num_local_dofs_velo; ++i)
        {
          value_v.axpy(space_data_velo.phi[i].value, local_dofs_v[i]);
          grad_v.add_outer_product(local_dofs_v[i], space_data_velo.phi[i].grad);
        }

        // compute local k and e
        value_k = value_e = DataType(0);

        for(int i = 0; i < num_local_dofs_velo; ++i)
        {
          value_k += local_dofs_k[i] * space_data_velo.phi[i].value;
          value_e += local_dofs_e[i] * space_data_velo.phi[i].value;
        }
      }

      /**
       * \brief Performs the assembly of the local matrices and vectors.
       *
       * This is where the magic happens.
       */
      void assemble()
      {
        // fetch coefficients
        const DataType my_beta = this->job.beta;
        const DataType my_beta_frechet = this->job.beta_frechet;
        const DataType my_theta = this->job.theta;

        // format local matrices and vectors
        local_matrix_a.format();
        local_vec_def.format();

        // loop over all quadrature points and integrate
        for(int point(0); point < this->cubature_rule.get_num_points(); ++point)
        {
          // prepare trafo and space for cubature point
          this->prepare_point(point);

          // precompute cubature weight
          const DataType weight = this->trafo_data.jac_det * this->cubature_rule.get_weight(point);

          // compute local nu
          DataType my_nu = this->job.nu;
          if(this->job.enable_k_e)
            my_nu += this->job.C_mu * value_k * value_k / value_e;

          // test function loop
          for(int i(0); i < this->num_local_dofs_velo; ++i)
          {
            // trial function loop
            for(int j(0); j < this->num_local_dofs_velo; ++j)
            {
              // add diffusion
              local_matrix_a[i][j].add_scalar_main_diag(my_nu * weight * Tiny::dot(this->space_data_velo.phi[i].grad, this->space_data_velo.phi[j].grad));
              if(this->job.deformation)
              {
                local_matrix_a[i][j].add_outer_product(this->space_data_velo.phi[j].grad, this->space_data_velo.phi[i].grad, my_nu * weight);
              }

              // add convection
              local_matrix_a[i][j].add_scalar_main_diag(my_beta * weight * this->space_data_velo.phi[i].value * Tiny::dot(value_v, this->space_data_velo.phi[j].grad));

              // add reaction
              local_matrix_a[i][j].add_scalar_main_diag(my_theta * weight *  this->space_data_velo.phi[i].value * this->space_data_velo.phi[j].value);

              // update defect vector
              //local_vec_def[i] -= local_matrix_a[i][j] * local_dofs_v[j];

              // add frechet derivative of convection
              local_matrix_a[i][j].axpy(my_beta_frechet * weight * this->space_data_velo.phi[i].value * this->space_data_velo.phi[j].value, this->grad_v);

            } // next trial function

          } // next test function

        } // continue with next cubature point

        for(int i(0); i < this->num_local_dofs_velo; ++i)
        {
          // trial function loop
          for(int j(0); j < this->num_local_dofs_velo; ++j)
          {
            local_vec_def[i] -= local_matrix_a[i][j] * local_dofs_v[j];
          }
        }
      }

      /// scatters the local matrices and vectors
      void scatter()
      {
        if(this->scatter_matrix_a)
          this->scatter_matrix_a->operator()(this->local_matrix_a, this->dof_mapping_velo, this->dof_mapping_velo);
        if(this->scatter_vec_def)
          this->scatter_vec_def->operator()(this->local_vec_def, this->dof_mapping_velo);
      }

      /// Finishes the assembly on the current cell.
      void finish()
      {
        // finish evaluators
        space_eval_velo.finish();
        trafo_eval.finish();

        // finish dof mapping
        dof_mapping_velo.finish();
      }

      /// Finalizes the assembly.
      void combine()
      {
        // nothing to do here
      }
    }; // class TurbBurgersAssemblyJob<...>::Task
  }; // class TurbBurgersAssemblyJob<...>
} // namespace Turb
