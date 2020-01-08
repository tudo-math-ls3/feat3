// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ASSEMBLY_BURGERS_ASSEMBLER_HPP
#define KERNEL_ASSEMBLY_BURGERS_ASSEMBLER_HPP

#include <kernel/assembly/asm_traits.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/global/vector.hpp>

namespace FEAT
{
  namespace Assembly
  {
    /**
     * \brief Burgers operator assembly class
     *
     * This class is responsible for assembling the scalar and vector-valued Burgers operators:
     *
     * \f[\mathbf{N}(v,u,\psi) := \nu \mathbf{L}(u,\psi) + \theta \mathbf{M}(u,\psi) + \beta \mathbf{K}(v,u,\psi) + \beta' \mathbf{K'(v,u,\psi)} + \delta \mathbf{S}(v,u,\psi)\f]
     *
     * where
     * - \b L is the diffusive operator, which is either
     *   - the <em>gradient tensor</em>:
     *     \f[\mathbf{L}(u,\psi) := \int_\Omega \nabla u \cdot \nabla \psi \f]
     *   - or the <em>deformation tensor</em>:
     *     \f[\mathbf{L}(u,\psi) := \frac{1}{4} \int_\Omega (\nabla+\nabla^\top) u : (\nabla+\nabla^\top) \psi\f]
     * - \b M is the reactive operator:
     *   \f[\mathbf{M}(u,\psi) := \int_\Omega u\psi\f]
     * - \b K is the convective operator:
     *   \f[\mathbf{K}(v,u,\psi) := \int_\Omega v\cdot \nabla u \psi\f]
     * - <b>K'</b> is the Frechet derivative of the convective operator:
     *   \f[\mathbf{K'}(v,u,\psi) := \int_\Omega \nabla v u \psi\f]
     * - \b S is the Samarskij-style streamline-diffusion stabilisation operator:
     *   \f[\mathbf{S}(v,u,\psi) := \sum_{T\in\mathcal{T}_h}\delta_{T}\int_T (v\cdot\nabla u)\cdot(v\cdot\nabla\psi)\f]
     *   where
     *   \f[\delta_T := \frac{h_T}{\|v\|_\Omega}\cdot\frac{2Re_T}{1+Re_T}\qquad\textnormal{and}\qquad Re_T := \frac{\|v\|_T\cdot h_T}{\nu_\mathbf{S}}\f]
     *
     * <b>Notes on Streamline Diffusion Stabilisation:</b>\n
     * The implementation of the streamline diffusion stabilisation is based on Turek's CFD book (see \cite TurekCFD,
     * pages 119 -- 123), however, in the formula above, the local stabilisation parameter \f$\delta_T\f$ is not pre-
     * multiplied by the global parameter \f$\delta^*\f$, compare eq. (3.75) in the book, as we use delta as a scaling
     * parameter for the whole operator \b S in the definition of the full operator \b N, see above.
     * To enable streamline diffusion stabilisation, you have to perform several tasks:
     * -# Set the global stabilisation parameter #sd_delta, which corresponds to the delta* parameter in the CFD book.
     *    According to the book (and previous FEAT versions), the value should be in range [0.1, 2], but according
     *    to my personal experience, even smaller values may be a good choice (e.g. 0.02).
     * -# Set the viscosity parameter #sd_nu, which is usually set equal to nu. The only reason for this extra parameter
     *    is so that you can use this assembler class to assemble the stabilisation without necessarily also assembling
     *    the diffusion operator, which is what would happen if you set nu to a non-zero value.
     * -# Set the (maximum) velocity norm #sd_v_norm of the convection field, which corresponds to \f$\|v\|_\Omega\f$ and
     *    is required for the computation of the local delta_T parameters. You can do this by either specifying the
     *    value directly or, more conveniently, use one of this class's set_sd_v_norm() functions to compute this norm
     *    from a (local or global) convection field vector.
     *
     * \author Peter Zajac
     */
    template<typename DataType_, typename IndexType_, int dim_>
    class BurgersAssembler
    {
    public:
      /// the datatype we use here
      typedef DataType_ DataType;

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

      /// scaling parameter for streamline diffusion stabilisation operator \b S
      DataType_ sd_delta;

      /// viscosity parameter nu_S for streamline diffusion (usually equal to nu)
      DataType_ sd_nu;

      /// velocity norm for streamline diffusion
      DataType_ sd_v_norm;

      /// default constructor
      BurgersAssembler() :
        deformation(false),
        nu(DataType_(1)),
        theta(DataType_(0)),
        beta(DataType_(0)),
        frechet_beta(DataType_(0)),
        sd_delta(DataType_(0)),
        sd_nu(DataType_(0)),
        sd_v_norm(DataType_(0))
      {
      }

      /**
       * \brief Assembles the Burgers operator into a matrix.
       *
       * \param[in,out] matrix
       * The matrix to be assembled.
       *
       * \param[in] convect
       * The transport vector for the convection.
       *
       * \param[in] space
       * The velocity space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] scale
       * A scaling factor for the matrix to be assembled.
       */
      template<typename Space_>
      void assemble_matrix(
        LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, dim_>& matrix,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect,
        const Space_& space,
        const Cubature::DynamicFactory& cubature_factory,
        const DataType_ scale = DataType_(1)
        ) const
      {
        // validate matrix and vector dimensions
        XASSERTM(matrix.rows() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.columns() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;
        typedef LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, dim_> MatrixType;

        // tolerance: sqrt(eps)
        const DataType tol_eps = Math::sqrt(Math::eps<DataType>());

        // first of all, let's see what we have to assemble
        const bool need_diff = (Math::abs(nu) > DataType(0));
        const bool need_conv = (Math::abs(beta) > DataType(0));
        const bool need_conv_frechet = (Math::abs(frechet_beta) > DataType(0));
        const bool need_reac = (Math::abs(theta) > DataType(0));
        const bool need_streamdiff = (Math::abs(sd_delta) > DataType(0)) && (sd_v_norm > tol_eps);

        // define our assembly traits
        typedef AsmTraits1<DataType_, Space_, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

        // fetch our trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        typename MatrixType::ScatterAxpy scatter_matrix(matrix);

        // create convection gather-axpy
        typename VectorType::GatherAxpy gather_conv(convect);

        // get maximum number of local dofs
        static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

        // create local matrix data
        typedef Tiny::Matrix<DataType, dim_, dim_> MatrixValue;
        typedef Tiny::Matrix<MatrixValue, max_local_dofs, max_local_dofs> LocalMatrixType;
        LocalMatrixType local_matrix;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs> LocalVectorType;

        // local convection field dofs
        LocalVectorType local_conv_dofs;

        // our local velocity value
        Tiny::Vector<DataType, dim_> loc_v, mean_v;

        // our local velocity gradient
        Tiny::Matrix<DataType, dim_, dim_> loc_grad_v;

        loc_v.format();
        mean_v.format();
        loc_grad_v.format();

        // compute reference element barycentre
        Tiny::Vector<DataType, dim_> barycentre;
        for(int i(0); i < dim_; ++i)
          barycentre[i] = Shape::ReferenceCell<typename Space_::ShapeType>::template centre<DataType>(i);

        // our local streamline diffusion coefficients
        Tiny::Vector<DataType, max_local_dofs> streamdiff_coeffs;
        streamdiff_coeffs.format();

        // our local delta for streamline diffusion
        DataType local_delta = DataType(0);

        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // initialise dof-mapping
          dof_mapping.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);

          // compute mesh width if necessary
          if(need_streamdiff)
          {
            // reset local delta
            local_delta = DataType(0);

            // evaluate trafo and space at barycentre
            trafo_eval(trafo_data, barycentre);
            space_eval(space_data, trafo_data);

            // compute velocity at barycentre
            mean_v.format();
            for(int i(0); i < num_loc_dofs; ++i)
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

          // format our local matrix and vector
          local_matrix.format();

          // loop over all quadrature points and integrate
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(point));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // pre-compute cubature weight
            const DataType weight = trafo_data.jac_det * cubature_rule.get_weight(point);

            // evaluate convection function and its gradient (if required)
            if(need_conv || need_streamdiff)
            {
              loc_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity value
                loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
              }
            }
            if(need_conv_frechet)
            {
              loc_grad_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity gradient
                loc_grad_v.add_outer_product(local_conv_dofs[i], space_data.phi[i].grad);
              }
            }

            // evaluate streamline diffusion coefficients
            if(need_streamdiff)
            {
              for(int i(0); i < num_loc_dofs; ++i)
              {
                streamdiff_coeffs[i] = Tiny::dot(loc_v, space_data.phi[i].grad);
              }
            }

            // assemble diffusion matrix?
            if(need_diff && !deformation)
            {
              // assemble gradient-tensor diffusion

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = nu * weight * Tiny::dot(space_data.phi[i].grad, space_data.phi[j].grad);

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }
            else if(need_diff && deformation)
            {
              // assemble deformation-tensor diffusion

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute inner product of grad(phi) and grad(psi)
                  const DataType value = nu * weight * Tiny::dot(space_data.phi[j].grad, space_data.phi[i].grad);

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);

                  // add outer product of grad(phi) and grad(psi)
                  local_matrix[i][j].add_outer_product(space_data.phi[j].grad, space_data.phi[i].grad, nu * weight);
                }
              }
            }

            // assemble convection?
            if(need_conv)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = beta * weight * space_data.phi[i].value * Tiny::dot(loc_v, space_data.phi[j].grad);

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }

            // assemble convection Frechet?
            if(need_conv_frechet)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = frechet_beta * weight * space_data.phi[i].value * space_data.phi[j].value;

                  // update local matrix
                  local_matrix[i][j].axpy(value, loc_grad_v);
                }
              }
            }

            // assemble reaction?
            if(need_reac)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = theta * weight *  space_data.phi[i].value * space_data.phi[j].value;

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }

            // assemble streamline diffusion?
            if(need_streamdiff && (local_delta > tol_eps))
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = local_delta * weight * streamdiff_coeffs[i] * streamdiff_coeffs[j];

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }

            // continue with next cubature point
          }

          // scatter into matrix
          scatter_matrix(local_matrix, dof_mapping, dof_mapping, scale);

          // finish dof mapping
          dof_mapping.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();
        }
      }

      /**
       * \brief Assembles the Burgers operator into a scalar matrix.
       *
       * \param[in,out] matrix
       * The scalar matrix to be assembled.
       *
       * \param[in] convect
       * The transport vector for the convection.
       *
       * \param[in] space
       * The velocity space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] scale
       * A scaling factor for the matrix to be assembled.
       */
      template<typename Matrix_, typename Space_>
      void assemble_scalar_matrix(
        Matrix_& matrix,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect,
        const Space_& space,
        const Cubature::DynamicFactory& cubature_factory,
        const DataType_ scale = DataType_(1)
        ) const
      {
        // validate matrix and vector dimensions
        XASSERTM(matrix.rows() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.columns() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;
        typedef Matrix_ MatrixType;

        // tolerance: sqrt(eps)
        const DataType tol_eps = Math::sqrt(Math::eps<DataType>());

        // first of all, let's see what we have to assemble
        const bool need_diff = (Math::abs(nu) > DataType(0));
        const bool need_conv = (Math::abs(beta) > DataType(0));
        const bool need_reac = (Math::abs(theta) > DataType(0));
        const bool need_streamdiff = (Math::abs(sd_delta) > DataType(0)) && (sd_v_norm > tol_eps);

        // deformation tensor is not available for scalar matrices
        XASSERTM(!deformation, "deformation tensor is not available for scalar matrices");
        XASSERTM(frechet_beta == DataType(0), "convection Frechet derivative is not available for scalar matrices");

        // define our assembly traits
        typedef AsmTraits1<DataType_, Space_, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

        // fetch our trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        typename MatrixType::ScatterAxpy scatter_matrix(matrix);

        // create convection gather-axpy
        typename VectorType::GatherAxpy gather_conv(convect);

        // get maximum number of local dofs
        static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

        // create local matrix data
        typedef typename AsmTraits::LocalMatrixType LocalMatrixType;
        LocalMatrixType local_matrix;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs> LocalVectorType;

        // local convection field dofs
        LocalVectorType local_conv_dofs;

        // our local velocity value
        Tiny::Vector<DataType, dim_> loc_v, mean_v;
        loc_v.format();
        mean_v.format();

        // compute reference element barycentre
        Tiny::Vector<DataType, dim_> barycentre;
        for(int i(0); i < dim_; ++i)
          barycentre[i] = Shape::ReferenceCell<typename Space_::ShapeType>::template centre<DataType>(i);

        // our local streamline diffusion coefficients
        Tiny::Vector<DataType, max_local_dofs> streamdiff_coeffs;
        streamdiff_coeffs.format();

        // our local delta for streamline diffusion
        DataType local_delta = DataType(0);

        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // initialise dof-mapping
          dof_mapping.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);

          // compute mesh width if necessary
          if(need_streamdiff)
          {
            // reset local delta
            local_delta = DataType(0);

            // evaluate trafo and space at barycentre
            trafo_eval(trafo_data, barycentre);
            space_eval(space_data, trafo_data);

            // compute velocity at barycentre
            mean_v.format();
            for(int i(0); i < num_loc_dofs; ++i)
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

          // format our local matrix and vector
          local_matrix.format();

          // loop over all quadrature points and integrate
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(point));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // pre-compute cubature weight
            const DataType weight = trafo_data.jac_det * cubature_rule.get_weight(point);

            // evaluate convection function and its gradient (if required)
            if(need_conv || need_streamdiff)
            {
              loc_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity value
                loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
              }
            }

            // evaluate streamline diffusion coefficients
            if(need_streamdiff)
            {
              for(int i(0); i < num_loc_dofs; ++i)
              {
                streamdiff_coeffs[i] = Tiny::dot(loc_v, space_data.phi[i].grad);
              }
            }

            // assemble diffusion?
            if(need_diff)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  local_matrix[i][j] += nu * weight * Tiny::dot(space_data.phi[i].grad, space_data.phi[j].grad);
                }
              }
            }

            // assemble convection?
            if(need_conv)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  local_matrix[i][j] += beta * weight * space_data.phi[i].value * Tiny::dot(loc_v, space_data.phi[j].grad);
                }
              }
            }

            // assemble reaction?
            if(need_reac)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  local_matrix[i][j] += theta * weight *  space_data.phi[i].value * space_data.phi[j].value;
                }
              }
            }

            // assemble streamline diffusion?
            if(need_streamdiff && (local_delta > tol_eps))
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  local_matrix[i][j] += local_delta * weight * streamdiff_coeffs[i] * streamdiff_coeffs[j];
                }
              }
            }

            // continue with next cubature point
          }

          // scatter into matrix
          scatter_matrix(local_matrix, dof_mapping, dof_mapping, scale);

          // finish dof mapping
          dof_mapping.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();
        }
      }

      /**
       * \brief Assembles the Burgers operator into a vector.
       *
       * \param[in,out] vector
       * The vector to be assembled.
       *
       * \param[in] convect
       * The transport vector for the convection.
       *
       * \param[in] primal
       * The primal vector, usually a solution vector.
       *
       * \param[in] space
       * The velocity space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] scale
       * A scaling factor the the vector to be assembled.
       */
      template<typename Space_>
      void assemble_vector(
        LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& vector,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& primal,
        const Space_& space,
        const Cubature::DynamicFactory& cubature_factory,
        const DataType_ scale = DataType_(1)
        ) const
      {
        // validate matrix and vector dimensions
        XASSERTM(vector.size() == space.get_num_dofs(), "invalid vector size");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;

        // first of all, let's see what we have to assemble
        const bool need_diff = (Math::abs(nu) > DataType(0));
        const bool need_conv = (Math::abs(beta) > DataType(0));
        //const bool need_conv_frechet = (frechet_beta > DataType(0));
        const bool need_reac = (Math::abs(theta) > DataType(0));

        // define our assembly traits
        typedef AsmTraits1<DataType_, Space_, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

        // fetch our trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create vector-scatter-axpy (if needed)
        typename VectorType::ScatterAxpy scatter_vector(vector);

        // create convection gather-axpy
        typename VectorType::GatherAxpy gather_conv(convect);

        // create primal gather-axpy
        typename VectorType::GatherAxpy gather_prim(primal);

        // get maximum number of local dofs
        static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs> LocalVectorType;
        LocalVectorType local_vector;

        // local convection field dofs
        LocalVectorType local_conv_dofs;

        // local primal vector dofs
        LocalVectorType local_prim_dofs;

        // our local velocity value
        Tiny::Vector<DataType, dim_> loc_v;

        // our local velocity gradient
        //Tiny::Matrix<DataType, dim_, dim_> loc_grad_v;

        loc_v.format();
        //loc_grad_v.format();

        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // initialise dof-mapping
          dof_mapping.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);

          // gather our local primal dofs
          local_prim_dofs.format();
          gather_prim(local_prim_dofs, dof_mapping);

          // format our local vector
          local_vector.format();

          // loop over all quadrature points and integrate
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(point));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // pre-compute cubature weight
            const DataType weight = trafo_data.jac_det * cubature_rule.get_weight(point);

            // evaluate convection function and its gradient (if required)
            if(need_conv)
            {
              loc_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity value
                loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
              }
            }
            /*if(need_conv_frechet)
            {
              loc_grad_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity gradient
                loc_grad_v.add_outer_product(local_conv_dofs[i], space_data.phi[i].grad);
              }
            }*/

            // assemble diffusion matrix?
            if(need_diff && !deformation)
            {
              // assemble gradient-tensor diffusion

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = nu * weight * Tiny::dot(space_data.phi[i].grad, space_data.phi[j].grad);

                  // update local vector
                  local_vector[i].axpy(value, local_prim_dofs[j]);
                }
              }
            }
            else if(need_diff && deformation)
            {
              // assemble deformation-tensor diffusion

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute inner product of grad(phi) and grad(psi)
                  const DataType value1 = nu * weight * Tiny::dot(space_data.phi[i].grad, space_data.phi[j].grad);

                  // compute outer product of grad(phi) and grad(psi)
                  const DataType value2 = nu * weight * Tiny::dot(local_prim_dofs[j], space_data.phi[i].grad);

                  // update local vector
                  local_vector[i].axpy(value1, local_prim_dofs[j]);
                  local_vector[i].axpy(value2, space_data.phi[j].grad);
                }
              }
            }

            // assemble convection?
            if(need_conv)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = beta * weight * space_data.phi[i].value * Tiny::dot(loc_v, space_data.phi[j].grad);

                  // update local vector
                  local_vector[i].axpy(value, local_prim_dofs[j]);
                }
              }
            }

            // assemble reaction?
            if(need_reac)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = theta * weight *  space_data.phi[i].value * space_data.phi[j].value;

                  // update local vector
                  local_vector[i].axpy(value, local_prim_dofs[j]);
                }
              }
            }

            // continue with next cubature point
          }

          // scatter into vector
          scatter_vector(local_vector, dof_mapping, scale);

          // finish dof mapping
          dof_mapping.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();
        }
      }

      /**
       * \brief Sets the convection field norm \f$\|v\|_\Omega\f$ for the local streamline diffusion parameter delta_T.
       *
       * \param[in] convect
       * The (local) convection field vector.
       */
      void set_sd_v_norm(const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect)
      {
        const auto* vals = convect.elements();
        DataType_ r = DataType(0);
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
       * The (global) convection field vector.
       */
      template<typename LocalVector_, typename Mirror_>
      void set_sd_v_norm(const Global::Vector<LocalVector_, Mirror_>& convect)
      {
        this->set_sd_v_norm(convect.local());
        const auto* gate = convect.get_gate();
        if(gate != nullptr)
          this->sd_v_norm = gate->max(this->sd_v_norm);
      }
    }; // class BurgersAssembler<...>
  } // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_BURGERS_ASSEMBLER_HPP
