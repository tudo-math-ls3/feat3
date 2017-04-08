#pragma once
#ifndef KERNEL_ASSEMBLY_BURGERS_ASSEMBLER_HPP
#define KERNEL_ASSEMBLY_BURGERS_ASSEMBLER_HPP

#include <kernel/assembly/asm_traits.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>

namespace FEAT
{
  namespace Assembly
  {
    /**
     * \brief Burgers operator assembly class
     *
     * This class is responsible for assembling the vector-valued Burgers operator:
     *
     * \f[\mathbf{N}(v,u,\psi) := -\nu \mathbf{L}(u,\psi) + \beta \mathbf{K}(v,u,\psi) + \theta \mathbf{M}(u,\psi)\f]
     *
     * where
     * - \f$\mathbf{L}(u,\psi)\f$ is either the
     *   - gradient tensor: \f$\int_\Omega \nabla u \cdot \nabla \psi \f$
     *   - deformation tensor: \f$\frac{1}{4} \int_\Omega (\nabla+\nabla^\top) u : (\nabla+\nabla^\top) \psi\f$
     * - \f$\mathbf{K}(v,u,\psi)\f$ is the convection operator: \f$\int_\Omega v\cdot \nabla u \psi\f$
     * - \f$\mathbf{M}(u,\psi)\f$ is the reaction operator: \f$\int_\Omega u\psi\f$
     *
     * \author Peter Zajac
     */
    template<typename DataType_, typename IndexType_, int dim_>
    class BurgersAssembler
    {
    public:
      typedef DataType_ DataType;

      // specifies whether to use the deformation tensor
      bool deformation;

      /// viscosity parameter: nu
      DataType_ nu;
      DataType_ nu_scale_matrix;
      DataType_ nu_scale_vector;

      /// convection scaling parameter:
      DataType_ beta;
      DataType_ beta_scale_matrix;
      DataType_ beta_scale_vector;

      /// convection Frechet scaling parameter:
      DataType_ frechet_beta;
      DataType_ frechet_beta_scale_matrix;
      //DataType_ frechet_beta_scale_vector;

      /// reaction scaling parameter:
      DataType_ theta;
      DataType_ theta_scale_matrix;
      DataType_ theta_scale_vector;

      BurgersAssembler() :
        deformation(false),
        nu(DataType_(1)),
        nu_scale_matrix(DataType(1)),
        nu_scale_vector(DataType(1)),
        beta(DataType_(0)),
        beta_scale_matrix(DataType(1)),
        beta_scale_vector(DataType(1)),
        frechet_beta(DataType_(0)),
        frechet_beta_scale_matrix(DataType(1)),
        //frechet_beta_scale_vector(DataType(1)),
        theta(DataType_(0)),
        theta_scale_matrix(DataType(1)),
        theta_scale_vector(DataType(1))
      {
      }

      /**
       * \brief Assembles the Burgers operator into a matrix and/or vector.
       *
       * \param[in] space
       * The velocity space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] convect
       * The transport vector for the convection.
       *
       * \param[in,out] matrix
       * A pointer to the matrix to be assembled. May be \c nullptr if no matrix is to be assembled.
       *
       * \param[in,out] vector
       * A pointer to the vector to be assembled. May be \c nullptr if no vector is to be assembled.
       *
       * \param[in] scale_matrix
       * A scaling factor for the matrix to be assembled.
       *
       * \param[in] scale_vector
       * A scaling factor the the vector to be assembled.
       */
      template<typename Space_>
      void assemble(
        const Space_& space,
        const Cubature::DynamicFactory& cubature_factory,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect,
        LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, dim_>* matrix = nullptr,
        LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>* vector = nullptr,
        const DataType_ scale_matrix = DataType_(1),
        const DataType_ scale_vector = DataType_(1)
        ) const
      {
        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;
        typedef LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, dim_> MatrixType;

        // first of all, let's see what we have to assemble
        const bool need_diff = (nu > DataType(0));
        const bool need_conv = (beta > DataType(0));
        const bool need_conv_frechet = (frechet_beta > DataType(0));
        const bool need_reac = (theta > DataType(0));

        const bool have_matrix = (matrix != nullptr);
        const bool have_vector = (vector != nullptr);

        // nothing to do?
        XASSERTM(have_matrix || have_vector, "Neither matrix nor vector to be assembled");

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

        // create matrix scatter-axpy (if needed)
        std::shared_ptr<typename MatrixType::ScatterAxpy> scatter_matrix =
          (have_matrix ? std::make_shared<typename MatrixType::ScatterAxpy>(*matrix) : nullptr);

        // create vector-scatter-axpy (if needed)
        std::shared_ptr<typename VectorType::ScatterAxpy> scatter_vector =
          (have_vector ? std::make_shared<typename VectorType::ScatterAxpy>(*vector) : nullptr);

        // create velocity gather-axpy
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
        LocalVectorType local_vector;

        // local convection field dofs
        LocalVectorType local_conv_dofs;

        // our local velocity value
        Tiny::Vector<DataType, dim_> loc_v;

        // our local velocity gradient
        Tiny::Matrix<DataType, dim_, dim_> loc_grad_v;

        loc_v.format();
        loc_grad_v.format();

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

          // gather our local velocity dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);

          // format our local matrix and vector
          local_matrix.format();
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

            // evaluate velocity function and its gradient (if required)
            if(need_conv)
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

            // assemble diffusion matrix?
            if(need_diff && !deformation)
            {
              // assemble gradient-tensor diffusion
              const DataType nu_mat = nu * nu_scale_matrix;
              const DataType nu_vec = nu * nu_scale_vector;

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = weight * Tiny::dot(space_data.phi[i].grad, space_data.phi[j].grad);

                  // update local matrix
                  for(int k(0); k < dim_; ++k)
                  {
                    // update matrix
                    local_matrix[i][j][k][k] += nu_mat * value;

                    // update vector
                    local_vector[i][k] += nu_vec * value * local_conv_dofs[j][k];
                  }
                }
              }
            }
            else if(need_diff && deformation)
            {
              /// \todo figure out the correct scaling factor (1/2 ?, 1/4 ?)
              // assemble deformation-tensor diffusion
              const DataType nu_mat = nu * nu_scale_matrix;
              const DataType nu_vec = nu * nu_scale_vector;

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute inner product of  grad(phi) and grad(psi)
                  const DataType value = weight * Tiny::dot(space_data.phi[i].grad, space_data.phi[j].grad);

                  // update local matrix
                  for(int k(0); k < dim_; ++k)
                  {
                    // update matrix
                    local_matrix[i][j][k][k] += nu_mat * value;

                    // update vector
                    local_vector[i][k] += nu_vec * value * local_conv_dofs[j][k];
                  }

                  // add outer product of grad(phi) and grad(psi)
                  local_matrix[i][j].add_outer_product(space_data.phi[i].grad, space_data.phi[j].grad, nu_mat*weight);

                  // update local vector
                  const DataType vecval = nu_vec * weight * Tiny::dot(local_conv_dofs[j], space_data.phi[j].grad);
                  for(int k(0); k < dim_; ++k)
                  {
                    local_vector[i][k] += vecval * space_data.phi[i].grad[k];
                  }
                }
              }
            }

            // assemble convection?
            if(need_conv)
            {
              const DataType beta_mat = beta * beta_scale_matrix;
              const DataType beta_vec = beta * beta_scale_vector;

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = weight * space_data.phi[i].value * Tiny::dot(loc_v, space_data.phi[j].grad);

                  // update local matrix
                  for(int k(0); k < dim_; ++k)
                  {
                    // update matrix
                    local_matrix[i][j][k][k] += beta_mat * value;

                    // update vector
                    local_vector[i][k] += beta_vec * value * local_conv_dofs[j][k];
                  }
                }
              }
            }

            // assemble convection Frechet?
            if(need_conv_frechet)
            {
              const DataType frechet_beta_mat = frechet_beta * frechet_beta_scale_matrix;
              //const DataType frechet_beta_vec = frechet_beta * frechet_beta_scale_vector;

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = weight * space_data.phi[i].value * space_data.phi[j].value;

                  // update local matrix
                  local_matrix[i][j].axpy(frechet_beta_mat * value, loc_grad_v);
                }
              }
            }

            // assemble reaction?
            if(need_reac)
            {
              const DataType theta_mat = theta * theta_scale_matrix;
              const DataType theta_vec = theta * theta_scale_vector;

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = weight *  space_data.phi[i].value * space_data.phi[j].value;

                  // update local matrix
                  for(int k(0); k < dim_; ++k)
                  {
                    // update matrix
                    local_matrix[i][j][k][k] += theta_mat * value;

                    // update vector
                    local_vector[i][k] += theta_vec * value * local_conv_dofs[j][k];
                  }
                }
              }
            }

            // continue with next cubature point
          }

          // scatter into matrix
          if(have_matrix)
          {
            scatter_matrix->operator()(local_matrix, dof_mapping, dof_mapping, scale_matrix);
          }

          // scatter into vector
          if(have_vector)
          {
            scatter_vector->operator()(local_vector, dof_mapping, scale_vector);
          }

          // finish dof mapping
          dof_mapping.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();
        }
      }
    };
  } // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_BURGERS_ASSEMBLER_HPP
