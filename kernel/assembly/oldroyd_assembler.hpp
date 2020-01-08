// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ASSEMBLY_OLDROYD_ASSEMBLER_HPP
#define KERNEL_ASSEMBLY_OLDROYD_ASSEMBLER_HPP 1

#include <kernel/assembly/asm_traits.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>

namespace FEAT
{
  namespace Assembly
  {
    /**
     * \brief Oldroyd-B operator assembly class
     *
     * This class is responsible for assembling the <em>upper convected time derivative</em> operator
     *
     * \f[ \lambda\sigma + \gamma v\cdot\nabla\sigma - \zeta(\nabla v\cdot\sigma + \sigma\cdot\nabla v^\top)\f]
     *
     * which is used in the Oldroyd-B model for the stress tensor \f$\sigma\f$ in a 3-field Stokes formulation.
     *
     * This class can assemble the above mentioned operator for both symmetric and unsymmetric representations of
     * \f$\sigma\f$, where stress tensor is encoded as a block vector in the following fashion:
     *
     * - 2D unsymmetric (4 components): \f$(\sigma_{11}, \sigma_{12}, \sigma_{21}, \sigma_{22})\f$
     * - 2D symmetric (3 components): \f$(\sigma_{11}, \sigma_{22}, \sigma_{12})\f$
     * - 3D unsymmetric (9 components): \f$(\sigma_{11}, \sigma_{12}, \sigma_{13}, \sigma_{21}, \sigma_{22}, \sigma_{23}, \sigma_{31}, \sigma_{32}, \sigma_{33})\f$
     * - 3D symmetric (6 components): \f$(\sigma_{11}, \sigma_{22}, \sigma_{33}, \sigma_{12}, \sigma_{23}, \sigma_{31})\f$
     *
     * \author Peter Zajac
     */
    class OldroydAssembler
    {
    public:
      /**
       * \brief Assembles the Oldroyd operator onto a stress matrix
       *
       * \param[inout] matrix
       * The matrix to be assembled.
       *
       * \param[in] convect
       * The velocity field vector \e v.
       *
       * \param[in] space_v
       * The velocity space.
       *
       * \param[in] space_s
       * The stress space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] lambda
       * The scaling factor \f$\lambda\f$ for the reactive term \f$\sigma\f$.
       *
       * \param[in] gamma
       * The scaling factor \f$\gamma\f$ for the convective term \f$v\cdot\nabla\sigma\f$
       *
       * \param[in] zeta
       * The scaling factor \f$\zeta\f$ for the complicated term \f$\nabla v\cdot\sigma + \sigma\cdot\nabla v^\top\f$.
       */
      template<typename DT_, typename IT_, typename SpaceV_, typename SpaceS_, int dim_, int nsc_>
      static void assemble_matrix(
        LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, nsc_, nsc_>& matrix,
        const LAFEM::DenseVectorBlocked<Mem::Main, DT_, IT_, dim_>& convect,
        const SpaceV_& space_v,
        const SpaceS_& space_s,
        const Cubature::DynamicFactory& cubature_factory,
        const DT_ lambda,
        const DT_ gamma = DT_(1),
        const DT_ zeta = DT_(1)
      )
      {
        // validate matrix and vector dimensions
        XASSERTM(matrix.rows() == space_s.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.columns() == space_s.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(convect.size() == space_v.get_num_dofs(), "invalid vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DT_, IT_, dim_> VectorType;
        typedef LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, nsc_, nsc_> MatrixType;

        // define our assembly traits
        // we will "abuse" the 'trial space' of this assembly space for the velocity space,
        // the stress space is both the test and the 'real' trial space here
        typedef AsmTraits2<DT_, SpaceS_, SpaceV_, TrafoTags::jac_det,
          SpaceTags::value|SpaceTags::grad, SpaceTags::value|SpaceTags::grad> AsmTraits;

        // out datatype
        typedef typename AsmTraits::DataType DataType;

        // fetch our trafo
        const typename AsmTraits::TrafoType& trafo = space_s.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create the space evaluators
        typename AsmTraits::TestEvaluator space_eval_s(space_s);
        typename AsmTraits::TrialEvaluator space_eval_v(space_v);

        // create the dof-mappings
        typename AsmTraits::TestDofMapping dof_mapping_s(space_s);
        typename AsmTraits::TrialDofMapping dof_mapping_v(space_v);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::TestEvalData space_data_s;
        typename AsmTraits::TrialEvalData space_data_v;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        typename MatrixType::ScatterAxpy scatter_matrix(matrix);

        // create convection gather-axpy
        typename VectorType::GatherAxpy gather_conv(convect);

        // get maximum number of local dofs
        static constexpr int max_local_dofs_s = AsmTraits::max_local_test_dofs;
        static constexpr int max_local_dofs_v = AsmTraits::max_local_trial_dofs;

        // create local matrix data
        typedef Tiny::Matrix<DataType, nsc_, nsc_> MatrixValue;
        typedef Tiny::Matrix<MatrixValue, max_local_dofs_s, max_local_dofs_s> LocalMatrixType;
        LocalMatrixType local_matrix;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs_v> LocalVectorType;

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
          space_eval_s.prepare(trafo_eval);
          space_eval_v.prepare(trafo_eval);

          // initialise dof-mapping
          dof_mapping_s.prepare(cell);
          dof_mapping_v.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs_s = space_eval_s.get_num_local_dofs();
          const int num_loc_dofs_v = space_eval_v.get_num_local_dofs();

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping_v);

          // format our local matrix and vector
          local_matrix.format();

          // loop over all quadrature points and integrate
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(point));

            // compute basis function data
            space_eval_s(space_data_s, trafo_data);
            space_eval_v(space_data_v, trafo_data);

            // pre-compute cubature weight
            const DataType weight = trafo_data.jac_det * cubature_rule.get_weight(point);

            // evaluate convection function and its gradient (if required)
            loc_v.format();
            loc_grad_v.format();
            for(int i(0); i < num_loc_dofs_v; ++i)
            {
              loc_v.axpy(space_data_v.phi[i].value, local_conv_dofs[i]);
              loc_grad_v.add_outer_product(local_conv_dofs[i], space_data_v.phi[i].grad);
            }

            // test function loop
            for(int i(0); i < num_loc_dofs_s; ++i)
            {
              // trial function loop
              for(int j(0); j < num_loc_dofs_s; ++j)
              {
                // add reactive term 'lambda*sigma'
                local_matrix[i][j].add_scalar_main_diag(lambda * weight *  space_data_s.phi[i].value * space_data_s.phi[j].value);

                // add scalar term 'gamma*v*grad(sigma)' for all components of sigma
                local_matrix[i][j].add_scalar_main_diag(gamma * weight * space_data_s.phi[i].value * Tiny::dot(loc_v, space_data_s.phi[j].grad));

                // evaluate actual Oldroyd operator
                core_eval(local_matrix[i][j], loc_grad_v, -zeta * weight * space_data_s.phi[i].value * space_data_s.phi[j].value);
              }
            }

            // continue with next cubature point
          }

          // scatter into matrix
          scatter_matrix(local_matrix, dof_mapping_s, dof_mapping_s, DataType(1));

          // finish dof mapping
          dof_mapping_v.finish();
          dof_mapping_s.finish();

          // finish evaluators
          space_eval_v.finish();
          space_eval_s.finish();
          trafo_eval.finish();
        }
      }

    protected:
      /**
       * \brief Auxiliary evaluation function: unsymmetric 2D version
       *
       * This function evaluates the 2D version of the term \f$\nabla v\cdot\sigma + \sigma\cdot\nabla v^\top\f$
       * without exploiting the symmetry of the stress tensor, thus resulting in a 4x4 matrix.
       *
       * The components of the stress tensor are assumed to be ordered as (sigma_11, sigma_12, sigma_21, sigma_22).
       */
      template<typename DT_>
      static void core_eval(Tiny::Matrix<DT_, 4, 4, 4, 4>& M, const Tiny::Matrix<DT_, 2, 2, 2, 2>& G, const DT_ omega)
      {
        // we have
        //
        // X :=       grad(v)        *           sigma            +            sigma           *        grad(v)^T
        //
        //    = [ dx(v1)   dy(v1) ]  *  [ sigma_11   sigma_12 ]   +   [ sigma_11   sigma_12 ]  *  [ dx(v1)   dx(v2) ]
        //      [ dx(v2)   dy(v2) ]     [ sigma_21   sigma_22 ]       [ sigma_21   sigma_22 ]     [ dy(v1)   dy(v2) ]
        //      \_________________/     \_____________________/       \_____________________/     \_________________/
        //              =: G                      =: S                          = S                      = G^T
        //
        //    = [  G_11     G_12  ]  *  [   S_1         S_2   ]   +   [   S_1         S_2   ]  *  [  G_11     G_21  ]
        //      [  G_21     G_22  ]     [   S_3         S_4   ]       [   S_3         S_4   ]     [  G_12     G_22  ]
        //
        //    = [ G_11*S_1 + G_12*S_3     G_11*S_2 + G_12*S_4 ]   +   [ S_1*G_11 + S_2*G_12     S_1*G_21 + S_2*G_22 ]
        //      [ G_21*S_1 + G_22*S_3     G_21*S_2 + G_22*S_4 ]   +   [ S_3*G_11 + S_4*G_12     S_3*G_21 + S_4*G_22 ]
        //
        //    = [ 2*G_11*S_1 + G_12*S_2 + G_12*S_3          G_21*S_1 + (G_11+G_22)*S_2 + G_12*S_4 ]
        //      [ G_21*S_1 + (G_11+G_22)*S_3 + G_12*S_4     G_21*S_2 + G_21*S_3 + 2*G_22*S_4      ]
        //
        //    = [ X_1   X_2 ]
        //      [ X_3   X_4 ]
        //
        // Now re-write definition of X as 4x4 matrix-vector-product of M and S
        //
        // [ X_1 ] = [ M_11  M_12  M_13  M_14 ] * [ S_1 ]
        // [ X_2 ]   [ M_21  M_22  M_23  M_24 ]   [ S_2 ]
        // [ X_3 ]   [ M_31  M_32  M_33  M_34 ]   [ S_3 ]
        // [ X_4 ]   [ M_41  M_42  M_43  M_44 ]   [ S_4 ]

        // X_1 = M_1. * S
        M(0,0) += omega * DT_(2) * G(0,0);        // 2 * G_11 * S_1
        M(0,1) += omega * G(0,1);                 // G_12 * S_2
        M(0,2) += omega * G(0,1);                 // G_12 * S_3
      //M(0,3) += DT_(0);                         // 0 * S_4

        // X_2 = M_2. * S
        M(1,0) += omega * G(1,0);                 // G_21 * S_1
        M(1,1) += omega * (G(0,0) + G(1,1));      // (G_11+G_22) * S_2
      //M(1,2) += DT_(0);                         // 0 * S_3
        M(1,3) += omega * G(0,1);                 // G_12 * S_4

        // X_3 = M_3. * S
        M(2,0) += omega * G(1,0);                 // G_21 * S_1
      //M(2,1) += DT_(0);                         // 0 * S_2
        M(2,2) += omega * (G(0,0) + G(1,1));      // (G_11+G_22) * S_3
        M(2,3) += omega * (G(0,1));               // G_12 * S_4

        // X_4 = M_4. * S
      //M(3,0) += DT_(0);                         // 0 * S_1
        M(3,1) += omega * G(1,0);                 // G_21 * S_2
        M(3,2) += omega * G(1,0);                 // G_21 * S_3
        M(3,3) += omega * DT_(2) * G(1,1);        // 2 * G_22 * S_4
      }

      /**
       * \brief Auxiliary evaluation function: symmetric 2D version
       *
       * This function evaluates the 2D version of the term \f$\nabla v\cdot\sigma + \sigma\cdot\nabla v^\top\f$
       * while exploiting the symmetry of the stress tensor, thus resulting in a (reduced) 3x3 matrix.
       *
       * The components of the stress tensor are assumed to be ordered as (sigma_11, sigma_22, sigma_12).
       */
      template<typename DT_>
      static void core_eval(Tiny::Matrix<DT_, 3, 3, 3, 3>& M, const Tiny::Matrix<DT_, 2, 2, 2, 2>& G, const DT_ omega)
      {
        // we have
        //
        // X :=       grad(v)        *           sigma            +            sigma           *        grad(v)^T
        //
        //    = [ dx(v1)   dy(v1) ]  *  [ sigma_11   sigma_12 ]   +   [ sigma_11   sigma_12 ]  *  [ dx(v1)   dx(v2) ]
        //      [ dx(v2)   dy(v2) ]     [ sigma_21   sigma_22 ]       [ sigma_21   sigma_22 ]     [ dy(v1)   dy(v2) ]
        //      \_________________/     \_____________________/       \_____________________/     \_________________/
        //              =: G                      =: S                          = S                      = G^T
        //
        //    = [  G_11     G_12  ]  *  [   S_1         S_3   ]   +   [   S_1         S_3   ]  *  [  G_11     G_21  ]
        //      [  G_21     G_22  ]     [   S_3         S_2   ]       [   S_3         S_2   ]     [  G_12     G_22  ]
        //
        // At this point we exploit the symmetry of sigma, i.e. we have S_3 := sigma_12 = sigma_21, then:
        //
        //    = [ G_11*S_1 + G_12*S_3     G_11*S_3 + G_12*S_2 ]   +   [ S_1*G_11 + S_3*G_12     S_1*G_21 + S_3*G_22 ]
        //      [ G_21*S_1 + G_22*S_3     G_21*S_3 + G_22*S_2 ]   +   [ S_3*G_11 + S_2*G_12     S_3*G_21 + S_2*G_22 ]
        //
        //    = [ 2*G_11*S_1 + 2*G_12*S_3                   G_21*S_1 + G_12*S_2 + (G_11+G_22)*S_3 ]
        //      [ G_21*S_1 + G_12*S_2 + (G_11+G_22)*S_3     2*G_22*S_2 + 2*G_21*S_3               ]
        //
        //    = [ X_1   X_3 ]
        //      [ X_3   X_2 ]
        //
        // At this point, we exploit the symmetry of X and it as 3x3 matrix-vector-product of M and S
        //
        // [ X_1 ] = [ M_11  M_12  M_13 ] * [ S_1 ]
        // [ X_2 ]   [ M_21  M_22  M_23 ]   [ S_3 ]
        // [ X_3 ]   [ M_31  M_32  M_33 ]   [ S_3 ]

        // X_1 = M_1. * S
        M(0,0) += omega * DT_(2) * G(0,0);        // 2 * G_11 * S_1
      //M(0,1) += DT_(0);                         // 0 * S_2
        M(0,2) += omega * DT_(2) * G(0,1);        // 2 * G_12 * S_3

        // X_2 = M_2. * S
      //M(1,0) += DT_(0);                         // 0 * S_1
        M(1,1) += omega * DT_(2) * G(1,1);        // 2 * G_22 * S_2
        M(1,2) += omega * DT_(2) * G(1,0);        // 2 * G_21 * S_3

        // X_2 = M_2. * S
        M(2,0) += omega * G(1,0);                 // G_21 * S_1
        M(2,1) += omega * G(0,1);                 // G_12 * S_2
        M(2,2) += omega * (G(0,0) + G(1,1));      // (G_11+G_22) * S_3
      }
    }; // class OldroydAssembler
  } // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_OLDROYD_ASSEMBLER_HPP
