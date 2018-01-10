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
     * \f[\mathbf{N}(v,u,\psi) := \nu \mathbf{L}(u,\psi) + \beta \mathbf{K}(v,u,\psi) + \theta \mathbf{M}(u,\psi)\f]
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

      /// convection scaling parameter:
      DataType_ beta;

      /// convection Frechet scaling parameter:
      DataType_ frechet_beta;

      /// reaction scaling parameter:
      DataType_ theta;

      BurgersAssembler() :
        deformation(false),
        nu(DataType_(1)),
        beta(DataType_(0)),
        frechet_beta(DataType_(0)),
        theta(DataType_(0))
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

        // first of all, let's see what we have to assemble
        const bool need_diff = (Math::abs(nu) > DataType(0));
        const bool need_conv = (Math::abs(beta) > DataType(0));
        const bool need_conv_frechet = (frechet_beta > DataType(0));
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

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);

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
    }; // class BurgersAssembler<...>
  } // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_BURGERS_ASSEMBLER_HPP
