// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include "base.hpp"
#include <kernel/assembly/base.hpp>
#include <kernel/assembly/asm_traits.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/global/vector.hpp>

namespace Turb
{
  template<typename TurbMatrix_, typename TurbVector_, typename VeloVector_, typename SpaceTurb_, typename SpaceVelo_>
  class TurbSystemAssemblyJob
  {
  public:
    typedef TurbMatrix_ TurbMatrixType; // type of matrix_sys_k and matrix_sys_e
    typedef TurbVector_ TurbVectorType; // type of vector_k and vector_e
    typedef VeloVector_ VeloVectorType; // type of velocity vector component
    typedef SpaceTurb_ SpaceTurbType;
    typedef SpaceVelo_ SpaceVeloType;

    typedef typename TurbMatrixType::DataType DataType;

    static constexpr int dim = SpaceVeloType::shape_dim; // system dimension

    /// inout: pointer to k system matrix to be assembled
    TurbMatrixType* matrix_k = nullptr;

    /// inout: pointer to e system matrix to be assembled
    TurbMatrixType* matrix_e = nullptr;

    /// inout: pointer to k rhs vector to be assembled
    TurbVectorType* vec_rhs_k = nullptr;

    /// inout: pointer to e rhs vector to be assembled
    TurbVectorType* vec_rhs_e = nullptr;

    /// input: pointer to last k solution vector
    const TurbVectorType* vec_sol_k = nullptr;

    /// input: pointer to last k solution vector
    const TurbVectorType* vec_sol_e = nullptr;

    /// input: pointer to last v solution vector
    const VeloVectorType* vec_sol_v = nullptr;

    /// input: turbulence space type
    const SpaceTurbType* space_turb = nullptr;

    /// input: velocity space type
    const SpaceVeloType* space_velo = nullptr;
    TurbVectorType* ref = nullptr;

    TurbVectorType* p_k = nullptr;
    TurbVectorType* p_k_bnd = nullptr;
    TurbVectorType* p_k_neumann = nullptr;

    /// input: parameters
    DataType sigma_k = DataType(0);
    DataType sigma_e = DataType(0);
    DataType l_max = DataType(0);
    DataType C_mu = DataType(0);
    DataType nu_min = DataType(0);
    DataType nu = DataType(0);
    String problem = "dirichlet";

    /// input: reaction term theta for matrix k and e
    DataType theta_matrix_k = DataType(0);
    DataType theta_matrix_e = DataType(0);

    /// input: reaction term theta for rhs vectors
    DataType theta_vec_rhs = DataType(0);

    // input: scaling parameters for convection, diffusion, reaction
    DataType scal_diff = DataType(1);
    DataType scal_conv = DataType(1);
    DataType scal_reac = DataType(1);
    DataType scal_p_k = DataType(1);

  public:
    class Task
    {
    public:
      /// this task needs to scatter
      static constexpr bool need_scatter = true;

      /// this task doesn't need to combine
      static constexpr bool need_combine = false;

      /// our assembly traits
      typedef Assembly::AsmTraits2<DataType, SpaceVeloType, SpaceTurbType, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad, SpaceTags::value|SpaceTags::grad> AsmTraits;

      /// a reference to our job object
      const TurbSystemAssemblyJob& job;

      /// the test-/trial-space to be used
      const SpaceVeloType& space_velo;
      const SpaceTurbType& space_turb;

      /// the cubature factory used for integration
      const typename AsmTraits::TrafoType& trafo;

      /// the trafo evaluator
      typename AsmTraits::TrafoEvaluator trafo_eval;

      /// the space evaluator
      typename AsmTraits::TestEvaluator space_eval_velo;
      typename AsmTraits::TrialEvaluator space_eval_turb;

      /// the space dof-mapping
      typename AsmTraits::TestDofMapping dof_mapping_velo;
      typename AsmTraits::TrialDofMapping dof_mapping_turb;

      /// the cubature rule used for integration
      typename AsmTraits::CubatureRuleType cubature_rule;

      /// the trafo evaluation data
      typename AsmTraits::TrafoEvalData trafo_data;

      /// the space evaluation data
      typename AsmTraits::TestEvalData space_data_velo;
      typename AsmTraits::TrialEvalData space_data_turb;

      /// convection vector gather-axpy object
      typename VeloVectorType::GatherAxpy gather_vec_v;
      typename TurbVectorType::GatherAxpy gather_vec_k;
      typename TurbVectorType::GatherAxpy gather_vec_e;
      typename TurbVectorType::GatherAxpy gather_vec_ref;
      typename TurbVectorType::GatherAxpy gather_vec_p_k;
      typename TurbVectorType::GatherAxpy gather_vec_p_k_bnd;
      typename TurbVectorType::GatherAxpy gather_vec_p_k_neumann;

      /// matrix scatter objects
      typename TurbMatrixType::ScatterAxpy scatter_matrix_k;
      typename TurbMatrixType::ScatterAxpy scatter_matrix_e;

      /// vector scatter objects
      typename TurbVectorType::ScatterAxpy scatter_vec_rhs_k;
      typename TurbVectorType::ScatterAxpy scatter_vec_rhs_e;
      typename TurbVectorType::ScatterAxpy scatter_vec_rhs_ref;
      typename TurbVectorType::ScatterAxpy scatter_vec_rhs_p_k;
      typename TurbVectorType::ScatterAxpy scatter_vec_rhs_p_k_bnd;
      typename TurbVectorType::ScatterAxpy scatter_vec_rhs_p_k_neumann;

      /// maximum number of local dofs
      static constexpr int max_local_dofs_velo = AsmTraits::max_local_test_dofs;
      static constexpr int max_local_dofs_turb = AsmTraits::max_local_trial_dofs;

      /// actual number of local dofs on current element
      int num_local_dofs_velo, num_local_dofs_turb;

      /// local convection field dofs
      Tiny::Vector<Tiny::Vector<DataType, dim>, max_local_dofs_velo> local_dofs_v;

      /// local turbulence dofs
      Tiny::Vector<DataType, max_local_dofs_turb> local_dofs_k, local_dofs_e, local_dofs_ref;
      Tiny::Vector<DataType, max_local_dofs_turb>  local_dofs_p_k, local_dofs_p_k_bnd, local_dofs_p_k_neumann;

      /// local convection field value
      Tiny::Vector<DataType, dim> value_v;
      Tiny::Matrix<DataType, dim, dim> grad_v;

      /// local k and e values
      DataType value_k, value_e, velo_norm;

      /// the local k and e matrices to be assembled
      Tiny::Matrix<DataType, max_local_dofs_turb, max_local_dofs_turb> local_matrix_k, local_matrix_e;

      /// the local k and e rhs vectors to be assembled
      Tiny::Vector<DataType, max_local_dofs_turb> local_vec_rhs_k, local_vec_rhs_e, local_ref;
      Tiny::Vector<DataType, max_local_dofs_turb> local_p_k, local_p_k_bnd, local_p_k_neumann;

    public:
      /**
      * \brief Constructor
      *
      * \param[in] job_
      * A \resident reference to the job containing the assembly parameters.
      */
      explicit Task(const TurbSystemAssemblyJob& job_) :
        job(job_),
        space_velo(*job_.space_velo),
        space_turb(*job_.space_turb),
        trafo(space_turb.get_trafo()),
        trafo_eval(trafo),
        space_eval_velo(space_velo),
        space_eval_turb(space_turb),
        dof_mapping_velo(space_velo),
        dof_mapping_turb(space_turb),
        gather_vec_v(*job.vec_sol_v),
        gather_vec_k(*job.vec_sol_k),
        gather_vec_e(*job.vec_sol_e),
        gather_vec_ref(*job.ref),
        gather_vec_p_k(*job.p_k),
        gather_vec_p_k_bnd(*job.p_k_bnd),
        gather_vec_p_k_neumann(*job.p_k_neumann),
        scatter_matrix_k(*job.matrix_k),
        scatter_matrix_e(*job.matrix_e),
        scatter_vec_rhs_k(*job.vec_rhs_k),
        scatter_vec_rhs_e(*job.vec_rhs_e),
        scatter_vec_rhs_ref(*job.ref),
        scatter_vec_rhs_p_k(*job.p_k),
        scatter_vec_rhs_p_k_bnd(*job.p_k_bnd),
        scatter_vec_rhs_p_k_neumann(*job.p_k_neumann)
      {
        Cubature::DynamicFactory cubature_factory("gauss-legendre:3");
        cubature_factory.create(cubature_rule);
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
        dof_mapping_turb.prepare(cell);

        // prepare trafo evaluator
        trafo_eval.prepare(cell);

        // prepare space evaluator
        space_eval_velo.prepare(trafo_eval);
        space_eval_turb.prepare(trafo_eval);

        // fetch number of local dofs
        num_local_dofs_velo = space_eval_velo.get_num_local_dofs();
        num_local_dofs_turb = space_eval_turb.get_num_local_dofs();

        // gather our local dofs
        local_dofs_v.format();
        local_dofs_k.format();
        local_dofs_e.format();
        local_dofs_ref.format();
        local_dofs_p_k.format();
        local_dofs_p_k_bnd.format();
        local_dofs_p_k_neumann.format();
        gather_vec_v(local_dofs_v, dof_mapping_velo);
        gather_vec_k(local_dofs_k, dof_mapping_turb);
        gather_vec_e(local_dofs_e, dof_mapping_turb);
        gather_vec_ref(local_dofs_ref, dof_mapping_turb);
        gather_vec_p_k(local_dofs_p_k, dof_mapping_turb);
        gather_vec_p_k_bnd(local_dofs_p_k_bnd, dof_mapping_turb);
        gather_vec_p_k_neumann(local_dofs_p_k_neumann, dof_mapping_turb);
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
        space_eval_turb(space_data_turb, trafo_data);

        // compute local v
        value_v.format();
        grad_v.format();
        for(int i(0); i < num_local_dofs_velo; ++i)
        {
          value_v.axpy(space_data_velo.phi[i].value, local_dofs_v[i]);
          grad_v.add_outer_product(local_dofs_v[i], space_data_velo.phi[i].grad);
        }

        // compute local k and e
        velo_norm = value_v.norm_euclid();

        // compute local k and e
        value_k = value_e = DataType(0);
        local_ref.format();
        local_p_k.format();
        local_p_k_bnd.format();
        local_p_k_neumann.format();

        for(int i = 0; i < num_local_dofs_turb; ++i)
        {
          value_k += local_dofs_k[i] * space_data_turb.phi[i].value;
          value_e += local_dofs_e[i] * space_data_turb.phi[i].value;
          if(this->job.problem == "neumann")
          {
            local_p_k[i] = local_dofs_p_k[i] * space_data_turb.phi[i].value;
            local_p_k_bnd[i] = local_dofs_p_k_bnd[i] * space_data_turb.phi[i].value;
            local_p_k_neumann[i] = local_dofs_p_k_neumann[i] * space_data_turb.phi[i].value;
          }
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
        const DataType sig_k = this->job.sigma_k;
        const DataType sig_e = this->job.sigma_e;
        const DataType _C_mu = this->job.C_mu;
        const DataType _l_max = this->job.l_max;
        const DataType _nu_min = this->job.nu_min;
        const DataType _nu = this->job.nu;
        const DataType theta_mat_k = this->job.theta_matrix_k;
        const DataType theta_mat_e = this->job.theta_matrix_e;
        const DataType theta_vec = this->job.theta_vec_rhs;
        const DataType scale_diff = this->job.scal_diff;
        const DataType scale_conv = this->job.scal_conv;
        const DataType scale_reac = this->job.scal_reac;
        const DataType scale_p_k = this->job.scal_p_k;

        // format local matrices and vectors
        local_matrix_k.format();
        local_matrix_e.format();
        local_vec_rhs_k.format();
        local_vec_rhs_e.format();
        local_p_k.format();
        local_p_k_bnd.format();
        local_p_k_neumann.format();

        // Compute nu_T according to eq (24)
        const DataType value_nu_bnd = 0.41 * 11.06 * _nu;

        // loop over all quadrature points and integrate
        for(int point(0); point < this->cubature_rule.get_num_points(); ++point)
        {
          // prepare trafo and space for cubature point
          this->prepare_point(point);

          // precompute cubature weight
          const DataType weight = this->trafo_data.jac_det * this->cubature_rule.get_weight(point);

          // Compute nu_T according to (8)
          DataType l_star = _l_max;
          if (_C_mu * value_k * Math::sqrt(value_k) < value_e * _l_max)
          {
            l_star = _C_mu * value_k * Math::sqrt(value_k) / value_e;
          }

          DataType value_nu = Math::max(_nu_min, l_star * Math::sqrt(value_k));

          // Compute gamma according to eq (9)
          const DataType value_gamma = _C_mu*value_k/value_nu;

          // test function loop
          for(int i(0); i < this->num_local_dofs_turb; ++i)
          {
            if(job.problem == "neumann")
            {
              // Compute nu_T on boundary according to eq (24)
              if(Math::abs(local_dofs_ref[i]) > 0.1) value_nu = value_nu_bnd;
              else value_nu = Math::max(_nu_min, l_star * Math::sqrt(value_k));
            }
            // trial function loop
            for(int j(0); j < this->num_local_dofs_turb; ++j)
            {
              // compute scalar value
              const DataType laplace_op = Tiny::dot(this->space_data_turb.phi[i].grad, this->space_data_turb.phi[j].grad);
              const DataType mass_op = this->space_data_turb.phi[i].value * this->space_data_turb.phi[j].value;
              const DataType conv_op = Tiny::dot(this->value_v, this->space_data_turb.phi[j].grad) * this->space_data_turb.phi[i].value;

              // update local matrix k
              local_matrix_k[i][j] += weight * (scale_diff * value_nu / sig_k * laplace_op
                                              + scale_reac * theta_mat_k * value_gamma * mass_op
                                              + scale_conv * conv_op);

              // update local matrix e
              local_matrix_e[i][j] += weight * (scale_diff * value_nu / sig_e * laplace_op
                                              + scale_reac * theta_mat_e * value_gamma * mass_op
                                              + scale_conv * conv_op);
            } // next trial function

            // compute P_k
            DataType P_k = 0.0;
            DataType P_k_bnd = 0.0;
            Tiny::Matrix<DataType, dim, dim> tmp;
            tmp.set_transpose(this->grad_v);
            tmp.axpy(1.0, this->grad_v);
            P_k = value_nu / 2.0 * tmp.norm_frobenius() * tmp.norm_frobenius();

            if(job.problem == "dirichlet")
            {
              // update rhs vector k
              local_vec_rhs_k[i] += weight * (scale_p_k * P_k * this->space_data_turb.phi[i].value );

              // update rhs vector e
              local_vec_rhs_e[i] += weight * (value_gamma * theta_vec * scale_p_k * P_k * this->space_data_turb.phi[i].value);
            }
            if(job.problem == "neumann")
            {
              // Compute P_k according to eqs (22) and (28)
              P_k_bnd = Math::pow(Math::max(Math::pow(_C_mu, 0.25) * Math::sqrt(value_k), velo_norm / 11.06 ), 4 ) / value_nu;

              // Set values p_k
              local_p_k[i]      += P_k * weight * this->space_data_turb.phi[i].value;
              local_p_k_bnd[i]  += P_k_bnd * weight * this->space_data_turb.phi[i].value;
              if (Math::abs(local_dofs_ref[i]) > 0.1)
                local_p_k_neumann[i] += P_k_bnd * weight * this->space_data_turb.phi[i].value;
              else
                local_p_k_neumann[i] += P_k * weight * this->space_data_turb.phi[i].value;

              // update rhs vector k and e
              if (Math::abs(local_dofs_ref[i]) > 0.1)
              {
                local_vec_rhs_k[i] += weight * (scale_p_k * P_k_bnd * this->space_data_turb.phi[i].value );
                local_vec_rhs_e[i] += weight * (value_gamma * theta_vec * scale_p_k * P_k_bnd * this->space_data_turb.phi[i].value);
              }
              else
              {
                local_vec_rhs_k[i] += weight * (scale_p_k * P_k * this->space_data_turb.phi[i].value );
                local_vec_rhs_e[i] += weight * (value_gamma * theta_vec * scale_p_k * P_k * this->space_data_turb.phi[i].value);
              }
            }
          } // next test function
        } // continue with next cubature point
      }

      /// scatters the local matrices and vectors
      void scatter()
      {
        this->scatter_matrix_k(this->local_matrix_k, this->dof_mapping_turb, this->dof_mapping_turb);
        this->scatter_matrix_e(this->local_matrix_e, this->dof_mapping_turb, this->dof_mapping_turb);
        this->scatter_vec_rhs_k(this->local_vec_rhs_k, this->dof_mapping_turb);
        this->scatter_vec_rhs_e(this->local_vec_rhs_e, this->dof_mapping_turb);
        this->scatter_vec_rhs_ref(this->local_ref, this->dof_mapping_turb);
        this->scatter_vec_rhs_p_k(this->local_p_k, this->dof_mapping_turb);
        this->scatter_vec_rhs_p_k_bnd(this->local_p_k_bnd, this->dof_mapping_turb);
        this->scatter_vec_rhs_p_k_neumann(this->local_p_k_neumann, this->dof_mapping_turb);
      }

      /// Finishes the assembly on the current cell.
      void finish()
      {
        // finish evaluators
        space_eval_turb.finish();
        space_eval_velo.finish();
        trafo_eval.finish();

        // finish dof mapping
        dof_mapping_turb.finish();
        dof_mapping_velo.finish();
      }

      /// Finalizes the assembly.
      void combine()
      {
        // nothing to do here
      }
    }; // class TurbSystemAssemblyJob<...>::Task
  }; // class TurbSystemAssemblyJob<...>
} // namespace Turb
