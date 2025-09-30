
// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ASSEMBLY_BURGERS_VELO_MATERIAL_ASSEMBLY_JOB_HPP
#define KERNEL_ASSEMBLY_BURGERS_VELO_MATERIAL_ASSEMBLY_JOB_HPP 1

#include <kernel/assembly/base.hpp>
#include <kernel/assembly/asm_traits.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/global/vector.hpp>

namespace FEAT
{
  namespace Assembly
  {
    /**
     * \brief Base class for Burgers typed assembly jobs, using a function object for the actual viscosity model
     *
     * This base class is used by various Burgers assemblers to outsource common coefficient
     * member variables as well as code related to precomputing quantities required for the
     * streamline diffusion stabilization.
     *
     * This operator works simliar to the Burgers assembler, with one major change.
     * Instead of a constant viscosity, this assembler allows for a general viscosity dependent on
     * the second invariant of the deformation tensor:
     *    \f[\mathbf{\dot{\gamma}}(u) := \sqrt(\frac{1}{2})\sqrt{\nabla u : \nabla u^T } \f]
     *
     * You will have to provide two functions, mapping from \mathbb{R} \to \mathbb{R} :
     *   -  \b \nu(\dot{\gamma})
     *   -  \b \nu^\prime(\dot{\gamma})
     *
     * This can be provided by any object type that provides the () operator, for example a simple lambda.
     * \note The object should be trivial copyable, so probably do not use any captcha shenaningans.
     *
     *  The derivative is only used for the exact material jacobian. If you dont require this, simply pass a trivial function.
     *
     * \tparam DataType_
     * The datatype to be used for the coefficient member variables
     *
     * \tparam Space_
     * The finite element space to be used for assembly
     *
     * \tparam ViscFunc_
     *
     * \tparam ViscDerFunc_
     *
     * \tparam ConvVector_
     * The type of the convection (velocity) vector field
     *
     * \author Maximilian Esser
     */
    template<typename DataType_, typename Space_, typename ViscFunc_, typename ViscDerFunc_, typename ConvVector_>
    class BurgersVeloMaterialAssemblyJobBase
    {
    public:
      /// the datatype we use here
      typedef DataType_ DataType;

      /// the space type
      typedef Space_ SpaceType;

      /// the convection vector type
      typedef ConvVector_ ConvVectorType;

      /// the function object describing the viscosity function depending on the second invariant of the defo tensor
      typedef ViscFunc_ ViscosityFunctionType;

      /// the function object describing the derivitive of the viscosity function
      typedef ViscDerFunc_ ViscosityDerivFunctionType;

      /// Functionobject
      ViscosityFunctionType visco_func;

      /// Der Functionobject
      ViscosityDerivFunctionType visco_der_func;

      /// specifies whether to use the deformation tensor
      bool deformation;

      /// scaling parameter for diffusive operator \b L (aka viscosity)
      DataType_ nu;

      /// scaling parameter for reactive operator \b M
      DataType_ theta;

      /// scaling parameter for convective operator \b K
      DataType_ beta;

      /// scaling parameter for Frechet derivative of convective operator <b>K'</b>
      DataType_ frechet_beta;

      /// Pseudo temperature factor
      DataType aT;

      /// scaling parameter for streamline diffusion stabilization operator \b S
      DataType_ sd_delta;

      /// viscosity parameter nu_S for streamline diffusion (usually equal to nu)
      DataType_ sd_nu;

      /// velocity norm for streamline diffusion
      DataType_ sd_v_norm;

      /// Regularisation parameter for material derivative
      DataType_ reg_eps;

      /// the convection vector to be used
      const ConvVector_& convection_vector;

      /// the test-/trial-space to be used
      const Space_& space;

      /// the cubature factory to use
      Cubature::DynamicFactory cubature_factory;

      /// Assemble full Jacobian?
      DataType frechet_material;

      /// Scatter constant
      DataType alpha;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] conv_vector
       * A \resident reference to the convection vector for the convective term.
       *
       * \param[in] space_
       * A \resident reference to the space to be assembled on.
       *
       * \param[in] cubature
       * The name of the cubature rule to be used for assembly.
       */
      explicit BurgersVeloMaterialAssemblyJobBase(const ConvVector_& conv_vector, const Space_& space_, const String& cubature,
                                              ViscosityFunctionType _visc_fun, ViscosityDerivFunctionType _visc_der_fun) :
        visco_func(_visc_fun),
        visco_der_func(_visc_der_fun),
        deformation(false),
        nu(DataType_(0)),
        theta(DataType_(0)),
        beta(DataType_(0)),
        frechet_beta(DataType_(0)),
        aT(DataType(1)),
        sd_delta(DataType_(0)),
        sd_nu(DataType_(0)),
        sd_v_norm(DataType_(0)),
        reg_eps(DataType_(1E-100)),
        convection_vector(conv_vector),
        space(space_),
        cubature_factory(cubature),
        frechet_material(DataType(0)),
        alpha(DataType(1))
      {
      }

      /**
       * \brief Calculates the convection field norm \f$\|v\|_\Omega\f$ for the local streamline diffusion parameter delta_T.
       *
       * \param[in] convect
       * The \transient (local) convection field vector.
       *
       * \returns
       * The convection vector field norm.
       */
      template<typename IndexType_, int conv_dim_>
      static DataType calc_sd_v_norm(const LAFEM::DenseVectorBlocked<DataType_, IndexType_, conv_dim_>& convect)
      {
        const auto* vals = convect.elements();
        DataType_ r = DataType_(0);
        for(Index i(0); i < convect.size(); ++i)
          r = Math::max(r, vals[i].norm_euclid());
        return r;
      }

      /**
       * \brief Calculates the convection field norm \f$\|v\|_\Omega\f$ for the streamline diffusion parameter delta_T.
       *
       * \note
       * This function automatically syncs the norm over all processes by using the vector's gate.
       *
       * \param[in] convect
       * The \transient (global) convection field vector.
       *
       * \returns
       * The convection vector field norm.
       */
      template<typename LocalVector_, typename Mirror_>
      static DataType calc_sd_v_norm(const Global::Vector<LocalVector_, Mirror_>& convect)
      {
        const auto* gate = convect.get_gate();
        if(gate != nullptr)
          return gate->max(calc_sd_v_norm(convect.local()));
        else
          return calc_sd_v_norm(convect.local());
      }

      /**
       * \brief Sets the convection field norm \f$\|v\|_\Omega\f$ for the local streamline diffusion parameter delta_T.
       *
       * \param[in] convect
       * The \transient (local) convection field vector.
       */
      template<typename VectorType_>
      void set_sd_v_norm(const VectorType_& convect)
      {
        this->sd_v_norm = calc_sd_v_norm(convect);
      }
    }; // class BurgersVeloMaterialAssemblyJobBase<...>

    /**
     * \brief Base class for Burgers assembly tasks
     *
     * This base class is used by various Burgers assemblers to outsource common code related to
     * the computation of the local velocity field value as well as precomputing quantities for
     * the streamline diffusion stabilization.
     *
     * This base class already implements the prepare() and finish() member functions required
     * by the Task interface, so a derived class only needs to implement the assemble() and
     * scatter() functions.
     *
     * This function implements the prepare_point() function, which performs most of the
     * preprocessing that has to be performed for each cubature point and which can be called
     * by the derived class' assemble() function.
     *
     * \tparam Job_
     * The job class that encapsulates the most derived class from task this base class.
     * This Job_ class should derive (directly or indirectly) from BurgersVeloMaterialAssemblyJobBase,
     * because this task accesses a lot of the job's member variables directly.
     *
     * \tparam DataType_
     * The datatype to be used for the coefficient member variables
     *
     * \tparam Function_
     * The Function type. Should be anything that takes in the second invariant of the strain rate tensor and outputs a viscosity, i.e.
     *  Function_:  R -> R
     *
     * \author Maximilian Esser
     */
    template<typename Job_, typename DataType_, typename Function_>
    class BurgersVeloMaterialAssemblyTaskBase
    {
    public:
      /// this task needs to scatter
      static constexpr bool need_scatter = true;
      /// this task doesn't need to combine
      static constexpr bool need_combine = false;

      /// the data type to be used by the assembly
      typedef DataType_ DataType;

      /// the space type
      typedef typename Job_::SpaceType SpaceType;

      /// the convection vector type
      typedef typename Job_::ConvVectorType ConvVectorType;

      /// the shape dimension
      static constexpr int shape_dim = SpaceType::shape_dim;

      /// the convection vector dimension
      static constexpr int conv_dim = ConvVectorType::BlockSize;

      /// our assembly traits
      typedef AsmTraits1<DataType, SpaceType, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

      /// our function type
      typedef Function_ FunctionType;

      /// a reference to our job object
      const Job_& job;

      /// our function to be used
      const FunctionType visco_func;

      /// tolerance for checks == 0: sqrt(eps)
      const DataType tol_eps;

      /// specifies whether to use the deformation tensor
      const bool deformation;

      /// scaling parameter for diffusive operator \b L (aka viscosity)
      const DataType nu;

      /// scaling parameter for reactive operator \b M
      const DataType theta;

      /// scaling parameter for convective operator \b K
      const DataType beta;

      /// scaling parameter for Frechet derivative of convective operator <b>K'</b>
      const DataType frechet_beta;

      /// scaling parameter for streamline diffusion stabilization operator \b S
      const DataType sd_delta;

      /// viscosity parameter nu_S for streamline diffusion (usually equal to nu)
      const DataType sd_nu;

      /// velocity norm for streamline diffusion
      const DataType sd_v_norm;

      /// Scaling operator for Frechet material parameter
      const DataType material_frechet;

      /// Pseudo temperature factor
      const DataType aT;

      /// keep track what we need to assemble
      const bool need_diff, need_conv, need_conv_frechet, need_reac, need_streamdiff, need_material_frechet;

      /// Some local data needed for our material model
      DataType_ gamma_dot, reg_eps, nu_loc;

      /// enable computation of certain quantities
      bool need_loc_v, need_loc_grad_v, need_mean_v, need_local_h;

      /// the test-/trial-space to be used
      const SpaceType& space;
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
      /// convection vector gather-axpy object
      typename ConvVectorType::GatherAxpy gather_conv;

      /// maximum number of local dofs
      static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

      /// actual number of local dofs on current element
      int num_local_dofs;

      /// local convection field dofs
      Tiny::Vector<Tiny::Vector<DataType, conv_dim>, max_local_dofs> local_conv_dofs;

      /// local convection field value
      Tiny::Vector<DataType, conv_dim> loc_v, mean_v;

      /// local convection field gradient
      Tiny::Matrix<DataType, conv_dim, conv_dim> loc_grad_v;

      /// strain rate tensor
      Tiny::Matrix<DataType, conv_dim, conv_dim> strain_rate_tensor_2;

      /// local streamline diffusion coefficients
      Tiny::Vector<DataType, max_local_dofs> streamdiff_coeffs;

      /// reference element barycenter
      Tiny::Vector<DataType, shape_dim> barycenter;

      /// local mesh width for streamline diffusion
      DataType local_h;

      /// local delta for streamline diffusion
      DataType local_delta;

      /// local scatter value
      DataType alpha = DataType(1);

    public:
      /**
       * \brief Constructor
       *
       * \param[in] job_
       * A \resident reference to the job containing the assembly parameters.
       */
      explicit BurgersVeloMaterialAssemblyTaskBase(const Job_& job_) :
        job(job_),
        visco_func(job_.visco_func),
        tol_eps(Math::sqrt(Math::eps<DataType>())),
        deformation(job_.deformation),
        nu(job_.nu),
        theta(job_.theta),
        beta(job_.beta),
        frechet_beta(job_.frechet_beta),
        sd_delta(job_.sd_delta),
        sd_nu(job_.sd_nu),
        sd_v_norm(job_.sd_v_norm),
        material_frechet(job_.frechet_material),
        aT(job_.aT),
        need_diff(true),
        need_conv(Math::abs(beta) > DataType(0)),
        need_conv_frechet(Math::abs(frechet_beta) > DataType(0)),
        need_reac(Math::abs(theta) > DataType(0)),
        need_streamdiff((Math::abs(sd_delta) > DataType(0)) && (sd_v_norm > tol_eps)),
        need_material_frechet(material_frechet > DataType(0)),
        gamma_dot(),
        reg_eps(job_.reg_eps),
        nu_loc(nu),
        need_loc_v(need_conv),
        need_loc_grad_v(true),
        need_mean_v(need_streamdiff),
        need_local_h(need_streamdiff),
        space(job_.space),
        trafo(space.get_trafo()),
        trafo_eval(trafo),
        space_eval(space),
        dof_mapping(space),
        cubature_rule(Cubature::ctor_factory, job_.cubature_factory),
        trafo_data(),
        space_data(),
        gather_conv(job_.convection_vector),
        num_local_dofs(0),
        local_conv_dofs(),
        loc_v(),
        mean_v(),
        loc_grad_v(),
        strain_rate_tensor_2(),
        streamdiff_coeffs(),
        barycenter(),
        local_h(DataType(0)),
        local_delta(DataType(0)),
        alpha(job.alpha)
      {
        // compute reference element barycenter
        for(int i(0); i < shape_dim; ++i)
          barycenter[i] = Shape::ReferenceCell<typename SpaceType::ShapeType>::template centre<DataType>(i);
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
        dof_mapping.prepare(cell);

        // prepare trafo evaluator
        trafo_eval.prepare(cell);

        // prepare space evaluator
        space_eval.prepare(trafo_eval);

        // fetch number of local dofs
        num_local_dofs = space_eval.get_num_local_dofs();

        // gather our local convection dofs
        local_conv_dofs.format();
        gather_conv(local_conv_dofs, dof_mapping);

        // compute mesh width if necessary
        if(need_mean_v || need_local_h || need_streamdiff)
        {
          // reset local h and delta
          local_h = local_delta = DataType(0);

          // evaluate trafo and space at barycenter
          trafo_eval(trafo_data, barycenter);
          space_eval(space_data, trafo_data);

          // compute velocity at barycenter
          mean_v.format();
          for(int i(0); i < num_local_dofs; ++i)
            mean_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);

          // compute norm of mean velocity
          const DataType local_norm_v = mean_v.norm_euclid();

          // do we have a non-zero velocity?
          if(local_norm_v > tol_eps)
          {
            // compute local mesh width w.r.t. mean velocity
            local_h = trafo_eval.width_directed(mean_v) * local_norm_v;

            // compute local Re_T
            const DataType local_re = (local_norm_v * local_h) / this->sd_nu;

            // compute local delta
            local_delta = this->sd_delta * (local_h / this->sd_v_norm) * (DataType(2)*local_re) / (DataType(1) + local_re);
          }
        }
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
        space_eval(space_data, trafo_data);

        // evaluate convection function and its gradient (if required)
        if(need_loc_v || need_conv || need_streamdiff)
        {
          loc_v.format();
          for(int i(0); i < num_local_dofs; ++i)
          {
            // update velocity value
            loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
          }
        }
        if(need_loc_grad_v || need_conv_frechet)
        {
          loc_grad_v.format();
          for(int i(0); i < num_local_dofs; ++i)
          {
            // update velocity gradient
            loc_grad_v.add_outer_product(local_conv_dofs[i], space_data.phi[i].grad);
          }
        }

        // evaluate streamline diffusion coefficients
        if(need_streamdiff)
        {
          for(int i(0); i < num_local_dofs; ++i)
          {
            streamdiff_coeffs[i] = Tiny::dot(loc_v, space_data.phi[i].grad);
          }
        }

        //default to 1, only assemble if carreau is needed...
        strain_rate_tensor_2.set_transpose(loc_grad_v);
        strain_rate_tensor_2.axpy(DataType(1.0), loc_grad_v);
        gamma_dot = Math::sqrt(DataType(0.5))*strain_rate_tensor_2.norm_frobenius();

        //default to nu
        nu_loc = visco_func(gamma_dot, aT);
      }

      /**
       * \brief Finishes the assembly on the current cell.
       */
      void finish()
      {
        // finish evaluators
        space_eval.finish();
        trafo_eval.finish();

        // finish dof mapping
        dof_mapping.finish();
      }

      /**
       * \brief Finalizes the assembly.
       */
      void combine()
      {
        // nothing to do here
      }
    }; // class BurgersVeloMaterialAssemblyTaskBase<...>

    /**
     * \brief Base class for blocked Burgers material assembly tasks
     *
     * This base class is used by various Burgers assemblers to outsource common code related to
     * local matrix assembly based on the chosen parameters.
     *
     * \tparam Job_
     * The job class that encapsulates the most derived class from task this base class.
     * This Job_ class should derive (directly or indirectly) from BurgersVeloMaterialAssemblyJobBase,
     * because this task accesses a lot of the job's member variables directly.
     *
     * \tparam DataType_
     * The datatype to be used for the coefficient member variables
     *
     * \tparam block_size_
     * The size of the matrix/vector blocks
     *
     * \author Maximilian Esser
     */
    template<typename Job_, typename DataType_, int block_size_, typename ViscFun_, typename ViscDerFun_>
    class BurgersVeloMaterialBlockedAssemblyTaskBase :
      public BurgersVeloMaterialAssemblyTaskBase<Job_, DataType_, ViscFun_>
    {
    public:
      typedef BurgersVeloMaterialAssemblyTaskBase<Job_, DataType_, ViscFun_> BaseClass;

      typedef DataType_ DataType;
      typedef ViscDerFun_ ViscosityDerivFunctionType;

      using BaseClass::max_local_dofs;
      using BaseClass::num_local_dofs;

    protected:
      /// our local matrix data
      typedef Tiny::Matrix<DataType_, block_size_, block_size_> MatrixValue;
      typedef Tiny::Matrix<MatrixValue, max_local_dofs, max_local_dofs> LocalMatrixType;
      LocalMatrixType local_matrix;
      /// our Visco derivative func
      ViscosityDerivFunctionType visco_der_func;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] job_
       * A \resident reference to the job containing the assembly parameters.
       */
      explicit BurgersVeloMaterialBlockedAssemblyTaskBase(const Job_& job_) :
        BaseClass(job_),
        local_matrix(),
        visco_der_func(job_.visco_der_func)
      {
      }

      /**
       * \brief Assembles the Burger operator in a single cubature point
       *
       * \param[in] weight
       * The cubature weight, pre-multiplied by the Jacobian determinant
       */
      void assemble_burgers_point(const DataType weight)
      {
        // assemble diffusion matrix?
        if(this->need_diff)
        {
          // assemble gradient-tensor diffusion

          // test function loop
          for(int i(0); i < this->num_local_dofs; ++i)
          {
            // trial function loop
            for(int j(0); j < this->num_local_dofs; ++j)
            {
              // compute scalar value
              const DataType value = this->nu_loc * weight * Tiny::dot(this->space_data.phi[i].grad, this->space_data.phi[j].grad);

              // update local matrix
              local_matrix[i][j].add_scalar_main_diag(value);
            }
          }
          // assemble deformation tensor?
          if(this->deformation)
          {
            // test function loop
            for(int i(0); i < this->num_local_dofs; ++i)
            {
              // trial function loop
              for(int j(0); j < this->num_local_dofs; ++j)
              {
                // add outer product of grad(phi) and grad(psi)
                local_matrix[i][j].add_outer_product(this->space_data.phi[j].grad, this->space_data.phi[i].grad, this->nu_loc * weight);
              }
            }
            // frechet carreau term "only" makes sense in defo formulation
            if(this->need_material_frechet)
            {
              const DataType fac =  this->material_frechet * visco_der_func(this->gamma_dot, this->aT)/(this->gamma_dot + this->reg_eps);
              //std::cout << fac ;
              // test function loop
              for(int i(0); i < this->num_local_dofs; ++i)
              {
                Tiny::Vector<DataType, BaseClass::ConvVectorType::BlockSize> du_grad_phi;
                du_grad_phi.set_mat_vec_mult(this->strain_rate_tensor_2, this->space_data.phi[i].grad);

                // trial function loop
                for(int j(0); j < this->num_local_dofs; ++j)
                {
                  Tiny::Vector<DataType, BaseClass::ConvVectorType::BlockSize> du_grad_psi;
                  du_grad_psi.set_mat_vec_mult(this->strain_rate_tensor_2, this->space_data.phi[j].grad);
                  // add outer product of grad(phi) and grad(psi)
                  local_matrix[i][j].add_outer_product(du_grad_phi, du_grad_psi, fac * weight);
                }
              }
            }
          }
        }

        // assemble convection?
        if(this->need_conv)
        {
          // test function loop
          for(int i(0); i < this->num_local_dofs; ++i)
          {
            // trial function loop
            for(int j(0); j < this->num_local_dofs; ++j)
            {
              // compute scalar value
              const DataType value = this->beta * weight * this->space_data.phi[i].value * Tiny::dot(this->loc_v, this->space_data.phi[j].grad);

              // update local matrix
              local_matrix[i][j].add_scalar_main_diag(value);
            }
          }
        }

        // assemble convection Frechet?
        if(this->need_conv_frechet)
        {
          // test function loop
          for(int i(0); i < this->num_local_dofs; ++i)
          {
            // trial function loop
            for(int j(0); j < this->num_local_dofs; ++j)
            {
              // compute scalar value
              const DataType value = this->frechet_beta * weight * this->space_data.phi[i].value * this->space_data.phi[j].value;

              // update local matrix
              local_matrix[i][j].axpy(value, this->loc_grad_v);
            }
          }
        }

        // assemble reaction?
        if(this->need_reac)
        {
          // test function loop
          for(int i(0); i < this->num_local_dofs; ++i)
          {
            // trial function loop
            for(int j(0); j < this->num_local_dofs; ++j)
            {
              // compute scalar value
              const DataType value = this->theta * weight *  this->space_data.phi[i].value * this->space_data.phi[j].value;

              // update local matrix
              local_matrix[i][j].add_scalar_main_diag(value);
            }
          }
        }
      }

      /**
       * \brief Assembles the Streamline Diffusion stabilization operator in a single cubature point
       *
       * \param[in] weight
       * The cubature weight, pre-multiplied by the Jacobian determinant
       */
      void assemble_streamline_diffusion(const DataType weight)
      {
        // assemble streamline diffusion?
        if(this->need_streamdiff && (this->local_delta > this->tol_eps))
        {
          // test function loop
          for(int i(0); i < this->num_local_dofs; ++i)
          {
            // trial function loop
            for(int j(0); j < this->num_local_dofs; ++j)
            {
              // compute scalar value
              const DataType value = this->local_delta * weight * this->streamdiff_coeffs[i] * this->streamdiff_coeffs[j];

              // update local matrix
              local_matrix[i][j].add_scalar_main_diag(value);
            }
          }
        }
      }

      /**
       * \brief Performs the assembly of the local matrix.
       */
      void assemble()
      {
        local_matrix.format();

        // loop over all quadrature points and integrate
        for(int point(0); point < this->cubature_rule.get_num_points(); ++point)
        {
          // prepare trafo and space for cubature point
          this->prepare_point(point);

          // precompute cubature weight
          const DataType weight = this->trafo_data.jac_det * this->cubature_rule.get_weight(point);

          // assemble the Burgers operator
          this->assemble_burgers_point(weight);

          // assemble the streamline diffusion
          this->assemble_streamline_diffusion(weight);
        } // continue with next cubature point
      }
    }; // class BurgersVeloMaterialBlockedAssemblyTaskBase<...>

    /**
     * \brief Burgers material assembly job for block matrix
     *
     * This assembly job implements the blocked (i.e. vector-valued) Burgers operator assembly,
     * which is assembled into a matrix.
     *
     * See the documentation of the BurgersVeloMaterialAssemblyJobBase class template (which is a base class
     * of this class) for details about the terms and coefficients offered by this assembly job.
     *
     * \tparam Matrix_
     * The type of the blocked matrix that is to be assembled.
     *
     * \tparam Space_
     * The finite element space to be used for assembly
     *
     * \tparam ViscoFunc_
     *
     * \tparam ViscoDerFunc_
     *
     * \tparam ConvVector_
     * The type of the convection (velocity) vector field
     *
     * \author Maximilian Esser
     */
    template<typename Matrix_, typename Space_, typename ViscoFunc_, typename ViscoDerFunc_, typename ConvVector_ = typename Matrix_::VectorTypeR>
    class BurgersVeloMaterialBlockedMatrixAssemblyJob :
      public BurgersVeloMaterialAssemblyJobBase<typename Matrix_::DataType, Space_, ViscoFunc_, ViscoDerFunc_, ConvVector_>
    {
    public:
      /// the matrix type
      typedef Matrix_ MatrixType;
      /// the data type
      typedef typename MatrixType::DataType DataType;

      /// our base class
      typedef BurgersVeloMaterialAssemblyJobBase<DataType, Space_, ViscoFunc_, ViscoDerFunc_, ConvVector_> BaseClass;

      // no nonsense, please
      static_assert(Matrix_::BlockHeight == Matrix_::BlockWidth, "only square matrix blocks are supported here");

      /// the block size
      static constexpr int block_size = Matrix_::BlockHeight;

      /// the matrix to be assembled
      MatrixType& matrix;

    public:
      /**
       * \brief Constructor
       *
       * \param[in,out] matrix_
       * A \resident reference to the matrix that is to be assembled.
       *
       * \param[in] conv_vector
       * A \resident reference to the convection field vector.
       *
       * \param[in] space_
       * A \resident reference to the space that is be used for assembly.
       *
       * \param[in] cubature_name
       * The name of the cubature rule to be used for assembly.
       */
      explicit BurgersVeloMaterialBlockedMatrixAssemblyJob(Matrix_& matrix_, const ConvVector_& conv_vector,
        const Space_& space_, const String& cubature_name, ViscoFunc_ visc_fun, ViscoDerFunc_ visc_der_func) :
        BaseClass(conv_vector, space_, cubature_name, visc_fun, visc_der_func),
        matrix(matrix_)
      {
      }

    public:
      /// the actual assembly task
      class Task :
        public BurgersVeloMaterialBlockedAssemblyTaskBase<BurgersVeloMaterialBlockedMatrixAssemblyJob, DataType, block_size, ViscoFunc_, ViscoDerFunc_>
      {
      public:
        typedef BurgersVeloMaterialBlockedAssemblyTaskBase<BurgersVeloMaterialBlockedMatrixAssemblyJob, DataType, block_size, ViscoFunc_, ViscoDerFunc_> BaseClass;

      protected:
        /// matrix scatter-axpy object
        typename MatrixType::ScatterAxpy scatter_matrix;

      public:
        /// constructor
        explicit Task(const BurgersVeloMaterialBlockedMatrixAssemblyJob& job_) :
          BaseClass(job_),
          scatter_matrix(job_.matrix)
        {
        }

        // prepare, assemble and finish are already implemented in the base classes

        /// scatters the local matrix
        void scatter()
        {
          this->scatter_matrix(this->local_matrix, this->dof_mapping, this->dof_mapping, this->alpha);
        }
      }; // class BurgersVeloMaterialBlockedMatrixAssemblyJob::Task
    }; // class BurgersVeloMaterialBlockedMatrixAssemblyJob

    /**
     * \brief Burgers Material assembly job for block right-hand-side vector
     *
     * This assembly job implements the blocked (i.e. vector-valued) Burgers operator assembly,
     * which is assembled into a right-hand-side vector.
     *
     * See the documentation of the BurgersVeloMaterialAssemblyJobBase class template (which is a base class
     * of this class) for details about the terms and coefficients offered by this assembly job.
     *
     * \tparam Vector_
     * The type of the blocked vector that is to be assembled.
     *
     * \tparam Space_
     * The finite element space to be used for assembly
     *
     * \tparam ConvVector_
     * The type of the convection (velocity) vector field; usually identical to Vector_
     *
     * \author Maximilian Esser
     */
    template<typename Vector_, typename Space_, typename ViscFunc_, typename ViscDerFunc_, typename ConvVector_ = Vector_>
    class BurgersVeloMaterialBlockedVectorAssemblyJob :
      public BurgersVeloMaterialAssemblyJobBase<typename Vector_::DataType, Space_, ViscFunc_, ViscDerFunc_, ConvVector_>
    {
    public:
      /// the vector type
      typedef Vector_ VectorType;
      /// the data type
      typedef typename VectorType::DataType DataType;

      /// our base class
      typedef BurgersVeloMaterialAssemblyJobBase<DataType, Space_, ViscFunc_, ViscDerFunc_, ConvVector_> BaseClass;

      /// the block size
      static constexpr int block_size = Vector_::BlockSize;

      /// the RHS vector to be assembled
      VectorType& vector_rhs;

      /// the solution vector
      const VectorType& vector_sol;

    public:
      /**
       * \brief Constructor
       *
       * \param[in,out] rhs_vector_
       * A \resident reference to the RHS vector that is to be assembled.
       *
       * \param[in] sol_vector_
       * A \resident reference to the solution vector that the matrix is to be multiplied by.
       * May be the same object as \p conv_vector, but must not be the same as \p rhs_vector_.
       *
       * \param[in] conv_vector
       * A \resident reference to the convection field vector.
       * May be the same object as \p sol_vector, but must not be the same as \p rhs_vector_.
       *
       * \param[in] space_
       * A \resident reference to the space that is be used for assembly.
       *
       * \param[in] cubature_name
       * The name of the cubature rule to be used for assembly.
       */
      explicit BurgersVeloMaterialBlockedVectorAssemblyJob(Vector_& rhs_vector_, const Vector_& sol_vector_,
        const ConvVector_& conv_vector, const Space_& space_, const String& cubature_name,
        ViscFunc_ visc_fun, ViscDerFunc_ visc_der_func) :
        BaseClass(conv_vector, space_, cubature_name, visc_fun, visc_der_func),
        vector_rhs(rhs_vector_),
        vector_sol(sol_vector_)
      {
        XASSERTM(&rhs_vector_ != &sol_vector_, "rhs and solution vectors must not be the same object");
        // compare void addresses to avoid warnings in case the classes are different
        XASSERTM((void*)&rhs_vector_ != (void*)&conv_vector, "rhs and convection vectors must not be the same object");
      }

    public:
      /// the actual assembly task
      class Task :
        public BurgersVeloMaterialBlockedAssemblyTaskBase<BurgersVeloMaterialBlockedVectorAssemblyJob, DataType, block_size, ViscFunc_, ViscDerFunc_>
      {
      public:
        typedef BurgersVeloMaterialBlockedAssemblyTaskBase<BurgersVeloMaterialBlockedVectorAssemblyJob, DataType, block_size, ViscFunc_, ViscDerFunc_> BaseClass;

        using BaseClass::max_local_dofs;

      protected:
        /// sol vector gather-axpy object
        typename VectorType::GatherAxpy gather_vec_sol;
        /// rhs vector scatter-axpy object
        typename VectorType::ScatterAxpy scatter_vec_rhs;
        /// specifies whether the solution vector and the convection vector are the same object
        const bool sol_equals_conv;

        /// local solution vector dofs
        Tiny::Vector<Tiny::Vector<DataType, block_size>, max_local_dofs> local_sol_dofs;

        /// local rhs vector dofs
        Tiny::Vector<Tiny::Vector<DataType, block_size>, max_local_dofs> local_vector;

      public:
        /// constructor
        explicit Task(const BurgersVeloMaterialBlockedVectorAssemblyJob& job_) :
          BaseClass(job_),
          gather_vec_sol(job_.vector_sol),
          scatter_vec_rhs(job_.vector_rhs),
          sol_equals_conv((void*)&job_.vector_sol == (void*)&job_.convection_vector),
          local_sol_dofs(),
          local_vector()
        {
        }

        void prepare(const Index cell)
        {
          BaseClass::prepare(cell);

          // gather local solution dofs if required
          if(!sol_equals_conv)
          {
            this->local_sol_dofs.format();
            this->gather_vec_sol(this->local_sol_dofs, this->dof_mapping);
          }
        }

        /**
         * \brief Computes the local vector by multiplying the local matrix by the solution vector
         */
        void compute_local_vector()
        {
          local_vector.format();

          // compute local vector from local matrix and sol vector dofs
          if(!sol_equals_conv)
          {
            for(int i(0); i < this->num_local_dofs; ++i)
              for(int j(0); j < this->num_local_dofs; ++j)
                this->local_vector[i].add_mat_vec_mult(this->local_matrix(i,j), this->local_sol_dofs[j]);
          }
          else
          {
            for(int i(0); i < this->num_local_dofs; ++i)
              for(int j(0); j < this->num_local_dofs; ++j)
                this->local_vector[i].add_mat_vec_mult(this->local_matrix(i,j), this->local_conv_dofs[j]);
          }
        }

        void assemble()
        {
          // assemble the local matrix
          BaseClass::assemble();

          // compute the local vector from the matrix
          compute_local_vector();
        }

        /// scatters the local vector
        void scatter()
        {
          this->scatter_vec_rhs(this->local_vector, this->dof_mapping, this->alpha);
        }
      }; // class BurgersVeloMaterialBlockedVectorAssemblyJob::Task
    }; // class BurgersVeloMaterialBlockedVectorAssemblyJob<...>

    /**
     * \brief Base class for scalar Burgers Material assembly tasks
     *
     * This base class is used by various Burgers assemblers to outsource common code related to
     * local matrix assembly based on the chosen parameters.
     *
     * \tparam Job_
     * The job class that encapsulates the most derived class from task this base class.
     * This Job_ class should derive (directly or indirectly) from BurgersVeloMaterialAssemblyJobBase,
     * because this task accesses a lot of the job's member variables directly.
     *
     * \tparam DataType_
     * The datatype to be used for the coefficient member variables
     *
     * \TODO: The frechet material term is not correctly implemented (yet)
     *
     * \author Maximilian Esser
     */
    template<typename Job_, typename DataType_, typename ViscFunc_, typename ViscDerFunc_>
    class BurgersVeloMaterialScalarAssemblyTaskBase :
      public BurgersVeloMaterialAssemblyTaskBase<Job_, DataType_, ViscFunc_>
    {
    public:
      typedef BurgersVeloMaterialAssemblyTaskBase<Job_, DataType_, ViscFunc_> BaseClass;

      typedef DataType_ DataType;

      using BaseClass::max_local_dofs;
      using BaseClass::num_local_dofs;

    protected:
      /// our local matrix data
      typedef Tiny::Matrix<DataType, max_local_dofs, max_local_dofs> LocalMatrixType;
      LocalMatrixType local_matrix;
      ViscDerFunc_ visco_der_fun;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] job_
       * A \resident reference to the job containing the assembly parameters.
       */
      explicit BurgersVeloMaterialScalarAssemblyTaskBase(const Job_& job_) :
        BaseClass(job_),
        local_matrix(),
        visco_der_fun(job_.visc_der_fun)
      {
      }

      /**
       * \brief Assembles the Burger operator in a single cubature point \todo: Missing frechet term?
       *
       * \param[in] weight
       * The cubature weight, pre-multiplied by the Jacobian determinant
       */
      void assemble_burgers_point(const DataType weight)
      {
        // assemble diffusion matrix?
        if(this->need_diff)
        {
          // assemble gradient-tensor diffusion

          // test function loop
          for(int i(0); i < this->num_local_dofs; ++i)
          {
            // trial function loop
            for(int j(0); j < this->num_local_dofs; ++j)
            {
              local_matrix[i][j] += this->nu_loc * weight
                * Tiny::dot(this->space_data.phi[i].grad, this->space_data.phi[j].grad);
            }
          }

          // \TODO: additional frechet term?
        }

        // assemble convection?
        if(this->need_conv)
        {
          // test function loop
          for(int i(0); i < this->num_local_dofs; ++i)
          {
            // trial function loop
            for(int j(0); j < this->num_local_dofs; ++j)
            {
              local_matrix[i][j] += this->beta * weight
                * this->space_data.phi[i].value * Tiny::dot(this->loc_v, this->space_data.phi[j].grad);
            }
          }
        }

        // assemble reaction?
        if(this->need_reac)
        {
          // test function loop
          for(int i(0); i < this->num_local_dofs; ++i)
          {
            // trial function loop
            for(int j(0); j < this->num_local_dofs; ++j)
            {
              local_matrix[i][j] += this->theta * weight
                * this->space_data.phi[i].value * this->space_data.phi[j].value;
            }
          }
        }
      }

      /**
       * \brief Assembles the Streamline Diffusion stabilization operator in a single cubature point
       *
       * \param[in] weight
       * The cubature weight, pre-multiplied by the Jacobian determinant
       */
      void assemble_streamline_diffusion(const DataType weight)
      {
        // assemble streamline diffusion?
        if(this->need_streamdiff && (this->local_delta > this->tol_eps))
        {
          // test function loop
          for(int i(0); i < this->num_local_dofs; ++i)
          {
            // trial function loop
            for(int j(0); j < this->num_local_dofs; ++j)
            {
              local_matrix[i][j] += this->local_delta * weight
                * this->streamdiff_coeffs[i] * this->streamdiff_coeffs[j];
            }
          }
        }
      }

      /**
       * \brief Performs the assembly of the local matrix.
       */
      void assemble()
      {
        local_matrix.format();

        // loop over all quadrature points and integrate
        for(int point(0); point < this->cubature_rule.get_num_points(); ++point)
        {
          // prepare trafo and space for cubature point
          this->prepare_point(point);

          // precompute cubature weight
          const DataType weight = this->trafo_data.jac_det * this->cubature_rule.get_weight(point);

          // assemble the Burgers operator
          this->assemble_burgers_point(weight);

          // assembles the streamline diffusion stabilization
          this->assemble_streamline_diffusion(weight);
        } // continue with next cubature point
      }
    }; // class BurgersVeloMaterialScalarAssemblyTaskBase<...>

    /**
     * \brief Burgers Material assembly job for scalar matrix
     *
     * This assembly job implements the scalar Burgers operator assembly,
     * which is assembled into a matrix.
     *
     * See the documentation of the BurgersVeloMaterialAssemblyJobBase class template (which is a base class
     * of this class) for details about the terms and coefficients offered by this assembly job.
     *
     * \tparam Matrix_
     * The type of the scalar matrix that is to be assembled.
     *
     * \tparam Space_
     * The finite element space to be used for assembly
     *
     * \tparam ConvVector_
     * The type of the convection (velocity) vector field
     *
     * \author Maximilian Esser
     */
    template<typename Matrix_, typename Space_, typename ViscFunc_, typename ViscDerFunc_, typename ConvVector_>
    class BurgersVeloMaterialScalarMatrixAssemblyJob :
      public BurgersVeloMaterialAssemblyJobBase<typename Matrix_::DataType, Space_, ViscFunc_, ViscDerFunc_, ConvVector_>
    {
    public:
      /// the matrix type
      typedef Matrix_ MatrixType;
      /// the data type
      typedef typename MatrixType::DataType DataType;

      /// our base class
      typedef BurgersVeloMaterialAssemblyJobBase<DataType, Space_, ViscFunc_, ViscDerFunc_, ConvVector_> BaseClass;

      /// the matrix to be assembled
      MatrixType& matrix;

    public:
      /**
       * \brief Constructor
       *
       * \param[in,out] matrix_
       * A \resident reference to the matrix that is to be assembled.
       *
       * \param[in] conv_vector
       * A \resident reference to the convection field vector.
       *
       * \param[in] space_
       * A reference to the space that is be used for assembly.
       *
       * \param[in] cubature_name
       * The name of the cubature rule to be used for assembly.
       */
      explicit BurgersVeloMaterialScalarMatrixAssemblyJob(Matrix_& matrix_, const ConvVector_& conv_vector,
        const Space_& space_, const String& cubature_name,
        ViscFunc_ visc_fun, ViscDerFunc_ visc_der_fun) :
        BaseClass(conv_vector, space_, cubature_name, visc_fun, visc_der_fun),
        matrix(matrix_)
      {
      }

    public:
      /// the actual assembly task
      class Task :
        public BurgersVeloMaterialScalarAssemblyTaskBase<BurgersVeloMaterialScalarMatrixAssemblyJob, DataType, ViscFunc_, ViscDerFunc_>
      {
      public:
        typedef BurgersVeloMaterialScalarAssemblyTaskBase<BurgersVeloMaterialScalarMatrixAssemblyJob, DataType, ViscFunc_, ViscDerFunc_> BaseClass;

      protected:
        /// matrix scatter-axpy object
        typename MatrixType::ScatterAxpy scatter_matrix;

      public:
        /// constructor
        explicit Task(const BurgersVeloMaterialScalarMatrixAssemblyJob& job_) :
          BaseClass(job_),
          scatter_matrix(job_.matrix)
        {
        }

        // prepare, assemble and finish are already implemented in the base classes

        /// scatters the local matrix
        void scatter()
        {
          this->scatter_matrix(this->local_matrix, this->dof_mapping, this->dof_mapping, this->alpha);
        }
      }; // class BurgersVeloMaterialScalarMatrixAssemblyJob::Task
    }; // class BurgersVeloMaterialScalarMatrixAssemblyJob

    /**
     * \brief Burgers Material assembly job for scalar right-hand-side vector
     *
     * This assembly job implements the scalar Burgers operator assembly,
     * which is assembled into a right-hand-side vector.
     *
     * See the documentation of the BurgersAssemblyJobBase class template (which is a base class
     * of this class) for details about the terms and coefficients offered by this assembly job.
     *
     * \tparam Vector_
     * The type of the scalar vector that is to be assembled.
     *
     * \tparam Space_
     * The finite element space to be used for assembly
     *
     * \tparam ConvVector_
     * The type of the convection (velocity) vector field.
     *
     * \author Maximilian Esser
     */
    template<typename Vector_, typename Space_, typename ViscoFunc_, typename ViscoDerFunc_, typename ConvVector_>
    class BurgersVeloMaterialScalarVectorAssemblyJob :
      public BurgersVeloMaterialAssemblyJobBase<typename Vector_::DataType, Space_, ViscoFunc_, ViscoDerFunc_, ConvVector_>
    {
    public:
      /// the vector type
      typedef Vector_ VectorType;
      /// the data type
      typedef typename VectorType::DataType DataType;

      /// our base class
      typedef BurgersVeloMaterialAssemblyJobBase<DataType, Space_, ViscoFunc_, ViscoDerFunc_, ConvVector_> BaseClass;

      /// the RHS vector to be assembled
      VectorType& vector_rhs;

      /// the solution vector
      const VectorType& vector_sol;

    public:
      /**
       * \brief Constructor
       *
       * \param[in,out] rhs_vector_
       * A \resident reference to the RHS vector that is to be assembled.
       * Must not be the same object as \p sol_vector_.
       *
       * \param[in] sol_vector_
       * A \resident reference to the solution vector that the matrix is to be multiplied by.
       * Must not be the same object as \p rhs_vector_.
       *
       * \param[in] conv_vector
       * A \resident reference to the convection field vector.
       *
       * \param[in] space_
       * A \resident reference to the space that is be used for assembly.
       *
       * \param[in] cubature_name
       * The name of the cubature rule to be used for assembly.
       */
      explicit BurgersVeloMaterialScalarVectorAssemblyJob(Vector_& rhs_vector_, const Vector_& sol_vector_,
        const ConvVector_& conv_vector, const Space_& space_, const String& cubature_name) :
        BaseClass(conv_vector, space_, cubature_name),
        vector_rhs(rhs_vector_),
        vector_sol(sol_vector_)
      {
        XASSERTM(&rhs_vector_ != &sol_vector_, "rhs and solution vectors must not be the same object");
      }

    public:
      /// the actual assembly task
      class Task :
        public BurgersVeloMaterialScalarAssemblyTaskBase<BurgersVeloMaterialScalarVectorAssemblyJob, DataType, ViscoFunc_, ViscoDerFunc_>
      {
      public:
        typedef BurgersVeloMaterialScalarAssemblyTaskBase<BurgersVeloMaterialScalarVectorAssemblyJob, DataType, ViscoFunc_, ViscoDerFunc_> BaseClass;

        using BaseClass::max_local_dofs;

      protected:
        /// sol vector gather-axpy object
        typename VectorType::GatherAxpy gather_vec_sol;
        /// rhs vector scatter-axpy object
        typename VectorType::ScatterAxpy scatter_vec_rhs;

        /// local solution vector dofs
        Tiny::Vector<DataType, max_local_dofs> local_sol_dofs;

        /// local rhs vector dofs
        Tiny::Vector<DataType, max_local_dofs> local_vector;

      public:
        /// constructor
        explicit Task(const BurgersVeloMaterialScalarVectorAssemblyJob& job_) :
          BaseClass(job_),
          gather_vec_sol(job_.vector_sol),
          scatter_vec_rhs(job_.vector_rhs),
          local_sol_dofs(),
          local_vector()
        {
        }

        void prepare(const Index cell)
        {
          BaseClass::prepare(cell);

          // gather local solution vector dofs
          this->local_sol_dofs.format();
          this->gather_vec_sol(this->local_sol_dofs, this->dof_mapping);
        }

        /**
         * \brief Computes the local vector by multiplying the local matrix by the solution vector
         */
        void compute_local_vector()
        {
          // compute local vector from local matrix and sol vector dofs
          this->local_vector.format();
          for(int i(0); i < this->num_local_dofs; ++i)
            for(int j(0); j < this->num_local_dofs; ++j)
              this->local_vector[i] += this->local_matrix(i,j) * this->local_sol_dofs[j];
          //this->local_vector.set_mat_vec_mult(this->local_matrix, this->local_sol_dofs);
        }

        void assemble()
        {
          // assemble the local matrix
          BaseClass::assemble();

          // compute local vector from the matrix
          compute_local_vector();
        }

        /// scatters the local vector
        void scatter()
        {
          this->scatter_vec_rhs(this->local_vector, this->dof_mapping, this->alpha);
        }
      }; // class BurgersVeloMaterialScalarVectorAssemblyJob::Task
    }; // class BurgersVeloMaterialScalarVectorAssemblyJob<...>
  } // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_BURGERS_VELO_MATERIAL_ASSEMBLY_JOB_HPP
