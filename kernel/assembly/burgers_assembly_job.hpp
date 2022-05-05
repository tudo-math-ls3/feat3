// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ASSEMBLY_BURGERS_ASSEMBLY_JOB_HPP
#define KERNEL_ASSEMBLY_BURGERS_ASSEMBLY_JOB_HPP

#include <kernel/assembly/base.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/global/vector.hpp>

namespace FEAT
{
  namespace Assembly
  {
    /**
     * \brief Base class for Burgers assembly jobs
     *
     * This base class is used by various Burgers assemblers to outsource common coefficient
     * member variables as well as code related to precomputing quantities required for the
     * streamline diffusion stabilization.
     *
     * The operators assembled by the derived classes are the scalar and vector-valued Burgers operators:
     *
     * \f[\mathbf{N}(v,u,\psi) := \nu \mathbf{L}(u,\psi) + \theta \mathbf{M}(u,\psi) + \beta \mathbf{K}(v,u,\psi) + \beta' \mathbf{K'}(v,u,\psi) + \delta \mathbf{S}(v,u,\psi)\f]
     *
     * where
     * - \b L is the diffusive operator, which is either
     *   - the <em>gradient tensor</em>:
     *     \f[\mathbf{L}(u,\psi) := \int_\Omega \nabla u \cdot \nabla \psi~dx \f]
     *   - or the <em>deformation tensor</em> (vector-valued operators only):
     *     \f[\mathbf{L}(u,\psi) := \frac{1}{4} \int_\Omega (\nabla+\nabla^\top) u : (\nabla+\nabla^\top) \psi~dx\f]
     * - \b M is the reactive operator:
     *   \f[\mathbf{M}(u,\psi) := \int_\Omega u\psi~dx\f]
     * - \b K is the convective operator:
     *   \f[\mathbf{K}(v,u,\psi) := \int_\Omega (v\cdot \nabla) u \psi~dx\f]
     * - <b>K'</b> is the additional term for the Fr&eacute;chet derivative of the convective operator \b K
         (vector-valued operators only; see note below for details):
     *   \f[\mathbf{K'}(v,u,\psi) := \int_\Omega (\nabla v) u \psi~dx\f]
     * - \b S is the Samarskij-style streamline-diffusion stabilization operator (see note below for details):
     *   \f[\mathbf{S}(v,u,\psi) := \sum_{T\in\mathcal{T}_h}\delta_{T}\int_T (v\cdot\nabla u)\cdot(v\cdot\nabla\psi)~dx\f]
     *   where
     *   \f[\delta_T := \frac{h_T}{\|v\|_\Omega}\cdot\frac{2Re_T}{1+Re_T}\qquad\textnormal{and}\qquad Re_T := \frac{\|v\|_T\cdot h_T}{\nu_\mathbf{S}}\f]
     *
     * <b>Notes on the Fr&eacute;chet derivative of the convective term:</b>\n
     * The trilinearform \f$\mathbf{K}(v,u,\psi)\f$ represents the weak formulation of the convective term \f$v\cdot\nabla u\f$,
     * which is the generalization of the non-linear term \f$u\cdot\nabla u\f$. To obtain the Fr&eacute;chet derivative
     * \f$D\mathbf{K}\f$ of this non-linear operator in the case of \f$v=u\f$, which is required for the Newton method
     * on the Burgers equation, by definition of the Fr&eacute;chet derivative we have to compute
     * \f[D\mathbf{K}(v,u,\psi) = \bigg[\lim_{\varepsilon\rightarrow 0} \frac{\mathbf{K}(v + \varepsilon \hat{v}, u+\varepsilon \hat{u}, \psi)}{\varepsilon}\bigg]\bigg|_{\hat{u}=u,\hat{v}=v}\f]
     * with \f$\hat{u},\hat{v}\in V_h\f$ with \f$\|\hat{u}\| = \|\hat{v}\| = 1\f$. So, here goes:
       \f{align*}
       \lim_{\varepsilon\rightarrow 0} \frac{\mathbf{K}(v + \varepsilon \hat{v}, u+\varepsilon \hat{u}, \psi)}{\varepsilon}
       &=\lim_{\varepsilon\rightarrow 0} \int_\Omega \frac{(v + \varepsilon \hat{v})\cdot\nabla(u+\varepsilon \hat{u}) - v\cdot\nabla u}{\varepsilon}\psi~dx\\
       &= \lim_{\varepsilon\rightarrow 0} \int_\Omega\frac{v\cdot\nabla u +\varepsilon v\cdot\nabla \hat{u} + \varepsilon \hat{v}\cdot\nabla u+\varepsilon^2 \hat{v}\cdot\nabla \hat{u} - v\cdot\nabla u}{\varepsilon}\psi~dx\\
       &= \lim_{\varepsilon\rightarrow 0} \int_\Omega\frac{\varepsilon v\cdot\nabla \hat{u} + \varepsilon \hat{v}\cdot\nabla u+\varepsilon^2 \hat{v}\cdot\nabla \hat{u}}{\varepsilon}\psi~dx\\
       &= \lim_{\varepsilon\rightarrow 0}\frac{\varepsilon}{\varepsilon}\int_\Omega (v\cdot\nabla \hat{u})\psi~dx + \lim_{\varepsilon\rightarrow 0}\frac{\varepsilon}{\varepsilon}\int_\Omega (\hat{v}\cdot\nabla u)\psi~dx + \lim_{\varepsilon\rightarrow 0} \frac{\varepsilon^2}{\varepsilon} \int_\Omega (\hat{v}\cdot\nabla \hat{u})\psi~dx\\
       &= \int_\Omega (v\cdot\nabla \hat{u})\psi~dx + \int_\Omega (\hat{v}\cdot\nabla u)\psi~dx + \lim_{\varepsilon\rightarrow 0} \varepsilon \int_\Omega (\hat{v}\cdot\nabla \hat{u})\psi~dx\\
       &= \int_\Omega (v\cdot\nabla \hat{u})\psi~dx + \int_\Omega (\hat{v}\cdot\nabla u)\psi~dx
       \f}
     * Finally, we obtain our Fr&eacute;chet derivative by setting \f$\hat{u} := u\f$ and \f$\hat{v} = v\f$:
     * \f[D\mathbf{K}(v,u,\psi) = \underbrace{\int_\Omega (v\cdot\nabla u)\psi~dx}_{=~\mathbf{K}(v,u,\psi)} + \underbrace{\int_\Omega (v\cdot\nabla u)\psi~dx}_{=:~\mathbf{K}'(v,u,\psi)}\f]
     * Note that in the case of this assembly class, the Fr&eacute;chet derivative \f$D\mathbf{K}(v,u,\psi)\f$ is split
     * into its two additive parts \f$\mathbf{K}(v,u,\psi)\f$ and \f$\mathbf{K'}(v,u,\psi)\f$, which are individually
     * scaled by \f$\beta\f$ and \f$\beta'\f$, respectively, so one has to set both \f$\beta = \beta' = 1\f$ to use the
     * full Fr&eacute;chet derivative operator \f$D\mathbf{K}\f$.
     *
     * <b>Notes on Streamline Diffusion Stabilization:</b>\n
     * The implementation of the streamline diffusion stabilization is based on Turek's CFD book (see \cite TurekCFD,
     * pages 119 -- 123), however, in the formula above, the local stabilization parameter \f$\delta_T\f$ is not pre-
     * multiplied by the global parameter \f$\delta^*\f$, compare eq. (3.75) in the book, as we use delta as a scaling
     * parameter for the whole operator \b S in the definition of the full operator \b N, see above.
     * To enable streamline diffusion stabilization, you have to perform several tasks:
     * -# Set the global stabilization parameter #sd_delta, which corresponds to the delta* parameter in the CFD book.
     *    According to the book (and previous FEAT versions), the value should be in range [0.1, 2], but according
     *    to my personal experience, even smaller values may be a good choice (e.g. 0.02).
     * -# Set the viscosity parameter #sd_nu, which is usually set equal to nu. The only reason for this extra parameter
     *    to exist is so that you can use this assembler class to assemble the stabilization without necessarily also
     *    assembling the diffusion operator, which is what would happen if you set nu to a non-zero value.
     * -# Set the (maximum) velocity norm #sd_v_norm of the convection field, which corresponds to \f$\|v\|_\Omega\f$
     *    and is required for the computation of the local delta_T parameters. You can do this by either specifying the
     *    value directly or, more conveniently, use one of this class's set_sd_v_norm() functions to compute this norm
     *    from a (local or global) convection field vector.
     *
     * \tparam DataType_
     * The datatype to be used for the coefficient member variables
     *
     * \tparam Space_
     * The finite element space to be used for assembly
     *
     * \tparam ConvVector_
     * The type of the convection (velocity) vector field
     *
     * \author Peter Zajac
     */
    template<typename DataType_, typename Space_, typename ConvVector_>
    class BurgersAssemblyJobBase
    {
    public:
      /// the datatype we use here
      typedef DataType_ DataType;

      /// the space type
      typedef Space_ SpaceType;

      /// the convection vector type
      typedef ConvVector_ ConvVectorType;

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

      /// scaling parameter for streamline diffusion stabilization operator \b S
      DataType_ sd_delta;

      /// viscosity parameter nu_S for streamline diffusion (usually equal to nu)
      DataType_ sd_nu;

      /// velocity norm for streamline diffusion
      DataType_ sd_v_norm;

      /// the convection vector to be used
      const ConvVector_& convection_vector;

      /// the test-/trial-space to be used
      const Space_& space;

      /// the cubature factory to use
      Cubature::DynamicFactory cubature_factory;

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
      explicit BurgersAssemblyJobBase(const ConvVector_& conv_vector, const Space_& space_, const String& cubature) :
        deformation(false),
        nu(DataType_(0)),
        theta(DataType_(0)),
        beta(DataType_(0)),
        frechet_beta(DataType_(0)),
        sd_delta(DataType_(0)),
        sd_nu(DataType_(0)),
        sd_v_norm(DataType_(0)),
        convection_vector(conv_vector),
        space(space_),
        cubature_factory(cubature)
      {
      }

      /**
       * \brief Sets the convection field norm \f$\|v\|_\Omega\f$ for the local streamline diffusion parameter delta_T.
       *
       * \param[in] convect
       * The \transient (local) convection field vector.
       */
      template<typename IndexType_, int conv_dim_>
      void set_sd_v_norm(const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, conv_dim_>& convect)
      {
        const auto* vals = convect.elements();
        DataType_ r = DataType_(0);
        for(Index i(0); i < convect.size(); ++i)
          r = Math::max(r, vals[i].norm_euclid());
        this->sd_v_norm = r;
      }

      /**
       * \brief Sets the convection field norm \f$\|v\|_\Omega\f$ for the streamline diffusion parameter delta_T.
       *
       * \note
       * This function automatically syncs the norm over all processes by using the vector's gate.
       *
       * \param[in] convect
       * The \transient (global) convection field vector.
       */
      template<typename LocalVector_, typename Mirror_>
      void set_sd_v_norm(const Global::Vector<LocalVector_, Mirror_>& convect)
      {
        this->set_sd_v_norm(convect.local());
        const auto* gate = convect.get_gate();
        if(gate != nullptr)
          this->sd_v_norm = gate->max(this->sd_v_norm);
      }
    }; // class BurgersAssemblyJobBase<...>

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
     * This Job_ class should derive (directly or indirectly) from BurgersAssemblyJobBase,
     * because this task accesses a lot of the job's member variables directly.
     *
     * \tparam DataType_
     * The datatype to be used for the coefficient member variables
     *
     * \author Peter Zajac
     */
    template<typename Job_, typename DataType_>
    class BurgersAssemblyTaskBase
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

      /// a reference to our job object
      const Job_& job;

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

      /// keep track what we need to assemble
      const bool need_diff, need_conv, need_conv_frechet, need_reac, need_streamdiff;

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

      /// local streamline diffusion coefficients
      Tiny::Vector<DataType, max_local_dofs> streamdiff_coeffs;

      /// reference element barycenter
      Tiny::Vector<DataType, shape_dim> barycenter;

      /// local delta for streamline diffusion
      DataType local_delta;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] job_
       * A \resident reference to the job containing the assembly parameters.
       */
      explicit BurgersAssemblyTaskBase(const Job_& job_) :
        job(job_),
        tol_eps(Math::sqrt(Math::eps<DataType>())),
        deformation(job_.deformation),
        nu(job_.nu),
        theta(job_.theta),
        beta(job_.beta),
        frechet_beta(job_.frechet_beta),
        sd_delta(job_.sd_delta),
        sd_nu(job_.sd_nu),
        sd_v_norm(job_.sd_v_norm),
        need_diff(Math::abs(nu) > DataType(0)),
        need_conv(Math::abs(beta) > DataType(0)),
        need_conv_frechet(Math::abs(frechet_beta) > DataType(0)),
        need_reac(Math::abs(theta) > DataType(0)),
        need_streamdiff((Math::abs(sd_delta) > DataType(0)) && (sd_v_norm > tol_eps)),
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
        streamdiff_coeffs(),
        barycenter(),
        local_delta(DataType(0))
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
        if(need_streamdiff)
        {
          // reset local delta
          local_delta = DataType(0);

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
            const DataType local_h = trafo_eval.width_directed(mean_v) * local_norm_v;

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
        if(need_conv || need_streamdiff)
        {
          loc_v.format();
          for(int i(0); i < num_local_dofs; ++i)
          {
            // update velocity value
            loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
          }
        }
        if(need_conv_frechet)
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
    }; // class BurgersAssemblyTaskBase<...>

    /**
     * \brief Base class for blocked Burgers assembly tasks
     *
     * This base class is used by various Burgers assemblers to outsource common code related to
     * local matrix assembly based on the chosen parameters.
     *
     * \tparam Job_
     * The job class that encapsulates the most derived class from task this base class.
     * This Job_ class should derive (directly or indirectly) from BurgersAssemblyJobBase,
     * because this task accesses a lot of the job's member variables directly.
     *
     * \tparam DataType_
     * The datatype to be used for the coefficient member variables
     *
     * \tparam block_size_
     * The size of the matrix/vector blocks
     *
     * \author Peter Zajac
     */
    template<typename Job_, typename DataType_, int block_size_>
    class BurgersBlockedAssemblyTaskBase :
      public BurgersAssemblyTaskBase<Job_, DataType_>
    {
    public:
      typedef BurgersAssemblyTaskBase<Job_, DataType_> BaseClass;

      typedef DataType_ DataType;

      using BaseClass::max_local_dofs;
      using BaseClass::num_local_dofs;

    protected:
      /// our local matrix data
      typedef Tiny::Matrix<DataType_, block_size_, block_size_> MatrixValue;
      typedef Tiny::Matrix<MatrixValue, max_local_dofs, max_local_dofs> LocalMatrixType;
      LocalMatrixType local_matrix;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] job_
       * A \resident reference to the job containing the assembly parameters.
       */
      explicit BurgersBlockedAssemblyTaskBase(const Job_& job_) :
        BaseClass(job_),
        local_matrix()
      {
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
                const DataType value = this->nu * weight * Tiny::dot(this->space_data.phi[i].grad, this->space_data.phi[j].grad);

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
                  local_matrix[i][j].add_outer_product(this->space_data.phi[j].grad, this->space_data.phi[i].grad, this->nu * weight);
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
        } // continue with next cubature point
      }
    }; // class BurgersBlockedAssemblyTaskBase<...>

    /**
     * \brief Burgers assembly job for block matrix
     *
     * This assembly job implements the blocked (i.e. vector-valued) Burgers operator assembly,
     * which is assembled into a matrix.
     *
     * See the documentation of the BurgersAssemblyJobBase class template (which is a base class
     * of this class) for details about the terms and coefficients offered by this assembly job.
     *
     * \tparam Matrix_
     * The type of the blocked matrix that is to be assembled.
     *
     * \tparam Space_
     * The finite element space to be used for assembly
     *
     * \tparam ConvVector_
     * The type of the convection (velocity) vector field
     *
     * \author Peter Zajac
     */
    template<typename Matrix_, typename Space_, typename ConvVector_ = typename Matrix_::VectorTypeR>
    class BurgersBlockedMatrixAssemblyJob :
      public BurgersAssemblyJobBase<typename Matrix_::DataType, Space_, ConvVector_>
    {
    public:
      /// the matrix type
      typedef Matrix_ MatrixType;
      /// the data type
      typedef typename MatrixType::DataType DataType;

      /// our base class
      typedef BurgersAssemblyJobBase<DataType, Space_, ConvVector_> BaseClass;

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
      explicit BurgersBlockedMatrixAssemblyJob(Matrix_& matrix_, const ConvVector_& conv_vector,
        const Space_& space_, const String& cubature_name) :
        BaseClass(conv_vector, space_, cubature_name),
        matrix(matrix_)
      {
      }

    public:
      /// the actual assembly task
      class Task :
        public BurgersBlockedAssemblyTaskBase<BurgersBlockedMatrixAssemblyJob, DataType, block_size>
      {
      public:
        typedef BurgersBlockedAssemblyTaskBase<BurgersBlockedMatrixAssemblyJob, DataType, block_size> BaseClass;

      protected:
        /// matrix scatter-axpy object
        typename MatrixType::ScatterAxpy scatter_matrix;

      public:
        /// constructor
        explicit Task(const BurgersBlockedMatrixAssemblyJob& job_) :
          BaseClass(job_),
          scatter_matrix(job_.matrix)
        {
        }

        // prepare, assemble and finish are already implemented in the base classes

        /// scatters the local matrix
        void scatter()
        {
          this->scatter_matrix(this->local_matrix, this->dof_mapping, this->dof_mapping);
        }
      }; // class BurgersBlockedMatrixAssemblyJob::Task
    }; // class BurgersBlockedMatrixAssemblyJob

    /**
     * \brief Burgers assembly job for block right-hand-side vector
     *
     * This assembly job implements the blocked (i.e. vector-valued) Burgers operator assembly,
     * which is assembled into a right-hand-side vector.
     *
     * See the documentation of the BurgersAssemblyJobBase class template (which is a base class
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
     * \author Peter Zajac
     */
    template<typename Vector_, typename Space_, typename ConvVector_ = Vector_>
    class BurgersBlockedVectorAssemblyJob :
      public BurgersAssemblyJobBase<typename Vector_::DataType, Space_, ConvVector_>
    {
    public:
      /// the vector type
      typedef Vector_ VectorType;
      /// the data type
      typedef typename VectorType::DataType DataType;

      /// our base class
      typedef BurgersAssemblyJobBase<DataType, Space_, ConvVector_> BaseClass;

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
      explicit BurgersBlockedVectorAssemblyJob(Vector_& rhs_vector_, const Vector_& sol_vector_,
        const ConvVector_& conv_vector, const Space_& space_, const String& cubature_name) :
        BaseClass(conv_vector, space_, cubature_name),
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
        public BurgersBlockedAssemblyTaskBase<BurgersBlockedVectorAssemblyJob, DataType, block_size>
      {
      public:
        typedef BurgersBlockedAssemblyTaskBase<BurgersBlockedVectorAssemblyJob, DataType, block_size> BaseClass;

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
        explicit Task(const BurgersBlockedVectorAssemblyJob& job_) :
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

        void assemble()
        {
          // assemble the local matrix
          BaseClass::assemble();
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

        /// scatters the local vector
        void scatter()
        {
          this->scatter_vec_rhs(this->local_vector, this->dof_mapping);
        }
      }; // class BurgersBlockedVectorAssemblyJob::Task
    }; // class BurgersBlockedVectorAssemblyJob<...>

    /**
     * \brief Base class for scalar Burgers assembly tasks
     *
     * This base class is used by various Burgers assemblers to outsource common code related to
     * local matrix assembly based on the chosen parameters.
     *
     * \tparam Job_
     * The job class that encapsulates the most derived class from task this base class.
     * This Job_ class should derive (directly or indirectly) from BurgersAssemblyJobBase,
     * because this task accesses a lot of the job's member variables directly.
     *
     * \tparam DataType_
     * The datatype to be used for the coefficient member variables
     *
     * \author Peter Zajac
     */
    template<typename Job_, typename DataType_>
    class BurgersScalarAssemblyTaskBase :
      public BurgersAssemblyTaskBase<Job_, DataType_>
    {
    public:
      typedef BurgersAssemblyTaskBase<Job_, DataType_> BaseClass;

      typedef DataType_ DataType;

      using BaseClass::max_local_dofs;
      using BaseClass::num_local_dofs;

    protected:
      /// our local matrix data
      typedef Tiny::Matrix<DataType, max_local_dofs, max_local_dofs> LocalMatrixType;
      LocalMatrixType local_matrix;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] job_
       * A \resident reference to the job containing the assembly parameters.
       */
      explicit BurgersScalarAssemblyTaskBase(const Job_& job_) :
        BaseClass(job_),
        local_matrix()
      {
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
                local_matrix[i][j] += this->nu * weight
                  * Tiny::dot(this->space_data.phi[i].grad, this->space_data.phi[j].grad);
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
        } // continue with next cubature point
      }
    }; // class BurgersScalarAssemblyTaskBase<...>

    /**
     * \brief Burgers assembly job for scalar matrix
     *
     * This assembly job implements the scalar Burgers operator assembly,
     * which is assembled into a matrix.
     *
     * See the documentation of the BurgersAssemblyJobBase class template (which is a base class
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
     * \author Peter Zajac
     */
    template<typename Matrix_, typename Space_, typename ConvVector_>
    class BurgersScalarMatrixAssemblyJob :
      public BurgersAssemblyJobBase<typename Matrix_::DataType, Space_, ConvVector_>
    {
    public:
      /// the matrix type
      typedef Matrix_ MatrixType;
      /// the data type
      typedef typename MatrixType::DataType DataType;

      /// our base class
      typedef BurgersAssemblyJobBase<DataType, Space_, ConvVector_> BaseClass;

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
      explicit BurgersScalarMatrixAssemblyJob(Matrix_& matrix_, const ConvVector_& conv_vector,
        const Space_& space_, const String& cubature_name) :
        BaseClass(conv_vector, space_, cubature_name),
        matrix(matrix_)
      {
      }

    public:
      /// the actual assembly task
      class Task :
        public BurgersScalarAssemblyTaskBase<BurgersScalarMatrixAssemblyJob, DataType>
      {
      public:
        typedef BurgersScalarAssemblyTaskBase<BurgersScalarMatrixAssemblyJob, DataType> BaseClass;

      protected:
        /// matrix scatter-axpy object
        typename MatrixType::ScatterAxpy scatter_matrix;

      public:
        /// constructor
        explicit Task(const BurgersScalarMatrixAssemblyJob& job_) :
          BaseClass(job_),
          scatter_matrix(job_.matrix)
        {
        }

        // prepare, assemble and finish are already implemented in the base classes

        /// scatters the local matrix
        void scatter()
        {
          this->scatter_matrix(this->local_matrix, this->dof_mapping, this->dof_mapping);
        }
      }; // class BurgersScalarMatrixAssemblyJob::Task
    }; // class BurgersScalarMatrixAssemblyJob

    /**
     * \brief Burgers assembly job for scalar right-hand-side vector
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
     * \author Peter Zajac
     */
    template<typename Vector_, typename Space_, typename ConvVector_>
    class BurgersScalarVectorAssemblyJob :
      public BurgersAssemblyJobBase<typename Vector_::DataType, Space_, ConvVector_>
    {
    public:
      /// the vector type
      typedef Vector_ VectorType;
      /// the data type
      typedef typename VectorType::DataType DataType;

      /// our base class
      typedef BurgersAssemblyJobBase<DataType, Space_, ConvVector_> BaseClass;

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
      explicit BurgersScalarVectorAssemblyJob(Vector_& rhs_vector_, const Vector_& sol_vector_,
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
        public BurgersScalarAssemblyTaskBase<BurgersScalarVectorAssemblyJob, DataType>
      {
      public:
        typedef BurgersScalarAssemblyTaskBase<BurgersScalarVectorAssemblyJob, DataType> BaseClass;

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
        explicit Task(const BurgersScalarVectorAssemblyJob& job_) :
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

        void assemble()
        {
          // assemble the local matrix
          BaseClass::assemble();

          // compute local vector from local matrix and sol vector dofs
          this->local_vector.format();
          for(int i(0); i < this->num_local_dofs; ++i)
            for(int j(0); j < this->num_local_dofs; ++j)
              this->local_vector[i] += this->local_matrix(i,j) * this->local_sol_dofs[j];
          //this->local_vector.set_mat_vec_mult(this->local_matrix, this->local_sol_dofs);
        }

        /// scatters the local vector
        void scatter()
        {
          this->scatter_vec_rhs(this->local_vector, this->dof_mapping);
        }
      }; // class BurgersScalarVectorAssemblyJob::Task
    }; // class BurgersScalarVectorAssemblyJob<...>
  } // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_BURGERS_ASSEMBLY_JOB_HPP
