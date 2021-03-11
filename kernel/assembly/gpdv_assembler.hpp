// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ASSEMBLY_GPDV_ASSEMBLER_HPP
#define KERNEL_ASSEMBLY_GPDV_ASSEMBLER_HPP 1

#include <kernel/assembly/asm_traits.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>

namespace FEAT
{
  namespace Assembly
  {
    /**
     * \brief Grad(P)/Div(V) assembler class
     *
     * This class assembles the two bilinearforms \f$\nabla p\f$ and \f$\nabla\cdot v\f$
     * of the Stokes equations, from which the weak formulations are obtained by integration
     * by parts, i.e.
     *
     * \f[ b(\psi,\varphi) := - \int_\Omega \psi \nabla \varphi\f]
     * \f[ d(\varphi,\psi) := - \int_\Omega (\nabla \cdot \varphi) \psi\f]
     *
     * \note
     * Both bilinearforms are by default multiplied by -1 to obtain symmetry of the resulting Stokes matrix.
     *
     * \author Peter Zajac
     */
    class GradPresDivVeloAssembler
    {
    public:
      /**
       * \brief Assembles the B and D matrices
       *
       * \param[in,out] matrix_b, matrix_d
       * The two matrices to be assembled. If the matrices are empty, their structure is
       * assembled automatically.
       *
       * \param[in] space_velo
       * The velocity space.
       *
       * \param[in] space_pres
       * The pressure space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] scale_b, scale_d
       * The scaling factors for the matrices. These default to -1 due to integration by parts.
       */
      template<
        typename DataType_, typename IndexType_, int dim_,
        typename SpaceVelo_, typename SpacePres_,
        typename CubatureFactory_>
      static void assemble(
        LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, 1>& matrix_b,
        LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, 1, dim_>& matrix_d,
        const SpaceVelo_& space_velo,
        const SpacePres_& space_pres,
        const CubatureFactory_& cubature_factory,
        const DataType_ scale_b = -DataType_(1),
        const DataType_ scale_d = -DataType_(1)
        )
      {
        typedef DataType_ DataType;

        // ensure that the matrices have the correct dimensions
        static_assert(SpaceVelo_::shape_dim == dim_, "invalid matrix block sizes");

        // make sure that velocity and pressure have the same trafo
        XASSERTM((&space_velo.get_trafo()) == (&space_pres.get_trafo()),
          "Velocity and Pressure spaces must be defined on the same trafo!");

        typedef LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, 1> MatrixB;
        typedef LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, 1, dim_> MatrixD;

        // assembly traits
        typedef AsmTraits2<
          DataType,
          SpaceVelo_,
          SpacePres_,
          TrafoTags::jac_det,
          SpaceTags::grad,
          SpaceTags::value> AsmTraits;

        // fetch the trafo
        const typename AsmTraits::TrafoType& trafo = space_velo.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create space evaluators
        typename AsmTraits::TestEvaluator velo_eval(space_velo);
        typename AsmTraits::TrialEvaluator pres_eval(space_pres);

        // create dof-mappings
        typename AsmTraits::TestDofMapping velo_dof_mapping(space_velo);
        typename AsmTraits::TrialDofMapping pres_dof_mapping(space_pres);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::TestEvalData velo_data;
        typename AsmTraits::TrialEvalData pres_data;

        // maximum number of dofs
        static constexpr int max_velo_dofs = AsmTraits::max_local_test_dofs;
        static constexpr int max_pres_dofs = AsmTraits::max_local_trial_dofs;

        // create local matrix data
        //typename AsmTraits::LocalMatrixType lmd;
        typedef Tiny::Matrix<DataType, dim_, 1> BEntryType;
        typedef Tiny::Matrix<DataType, 1, dim_> DEntryType;
        Tiny::Matrix<BEntryType, max_velo_dofs, max_pres_dofs> local_b;
        Tiny::Matrix<DEntryType, max_pres_dofs, max_velo_dofs> local_d;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // first of all, check whether we need to assemble the matrix structures
        if(matrix_b.empty())
        {
          Adjacency::Graph graph_b(SymbolicAssembler::assemble_graph_std2(space_velo, space_pres));
          matrix_b = MatrixB(graph_b);
        }
        matrix_b.format();
        if(matrix_d.empty())
        {
          Adjacency::Graph graph_d(SymbolicAssembler::assemble_graph_std2(space_pres, space_velo));
          matrix_d = MatrixD(graph_d);
          matrix_d.format();
        }
        matrix_d.format();

        // create matrix scatter-axpy
        typename MatrixB::ScatterAxpy scatter_b(matrix_b);
        typename MatrixD::ScatterAxpy scatter_d(matrix_d);

        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluators
          velo_eval.prepare(trafo_eval);
          pres_eval.prepare(trafo_eval);

          // fetch number of local dofs
          const int num_loc_velo_dofs = velo_eval.get_num_local_dofs();
          const int num_loc_pres_dofs = pres_eval.get_num_local_dofs();

          // format local matrix
          local_b.format();
          local_d.format();

          // loop over all quadrature points and integrate
          for(int pt(0); pt < cubature_rule.get_num_points(); ++pt)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(pt));

            // compute basis function data
            velo_eval(velo_data, trafo_data);
            pres_eval(pres_data, trafo_data);

            // test function loop
            for(int i(0); i < num_loc_velo_dofs; ++i)
            {
              // trial function loop
              for(int j(0); j < num_loc_pres_dofs; ++j)
              {
                const DataType pv = trafo_data.jac_det * cubature_rule.get_weight(pt) * pres_data.phi[j].value;

                // dimension loop
                for(int k(0); k < dim_; ++k)
                {
                  local_b[i][j][k][0] += velo_data.phi[i].grad[k] * pv;
                  local_d[j][i][0][k] += velo_data.phi[i].grad[k] * pv;
                }

                // continue with next trial function
              }
              // continue with next test function
            }
            // continue with next cubature point
          }

          // finish evaluators
          pres_eval.finish();
          velo_eval.finish();
          trafo_eval.finish();

          // initialize dof-mappings
          velo_dof_mapping.prepare(cell);
          pres_dof_mapping.prepare(cell);

          // incorporate local matrix
          scatter_b(local_b, velo_dof_mapping, pres_dof_mapping, scale_b);
          scatter_d(local_d, pres_dof_mapping, velo_dof_mapping, scale_d);

          // finish dof mapping
          pres_dof_mapping.finish();
          velo_dof_mapping.finish();

          // continue with next cell
        }
      }
    };
  } // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_GPDV_ASSEMBLER_HPP
