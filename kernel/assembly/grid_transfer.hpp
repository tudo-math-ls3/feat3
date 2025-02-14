// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/util/exception.hpp>
#include <kernel/lafem/base.hpp>
#include <kernel/assembly/asm_traits.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/math.hpp>
#include <kernel/geometry/intern/coarse_fine_cell_mapping.hpp>
#include <kernel/geometry/mesh_permutation.hpp>
#include <kernel/trafo/inverse_mapping.hpp>

// includes, system
#include <set>
#include <map>
#include <vector>

namespace FEAT
{
  namespace Assembly
  {
    /**
     * \brief Grid-Transfer assembly class template
     *
     * This class template implements the assembly of grid transfer operators.
     *
     * \author Peter Zajac
     */
    class GridTransfer
    {
    public:
      /**
       * \brief Exception for singular mass matrix errors
       *
       * This exception is thrown if the inversion of a local mass matrix fails, thus resulting
       * in invalid entries in a prolongation/truncation/transfer matrix.
       */
      class LocalMassMatrixSingularException :
        public FEAT::Exception
      {
      public:
        LocalMassMatrixSingularException() :
          Exception("local mass matrix inversion error")
        {
        }
      }; // class LocalMassMatrixSingularException

      /**
       * \brief Assembles a prolongation matrix and its corresponding weight vector.
       *
       * To obtain the final prolongation matrix, one needs to invert the weight vector
       * component-wise and scale the matrix rows by the inverted weights afterwards.
       * This can be accomplished by the \c component_invert and \c scale_rows operations
       * of the vector and matrix containers, resp.
       *
       * \param[in,out] matrix
       * A \transient reference to the prolongation matrix that is to be assembled.
       *
       * \param[in,out] vector
       * A \transient reference to the weight vector for the prolongation matrix.
       *
       * \param[in] fine_space
       * A \transient reference to the fine-mesh test-space to be used.
       *
       * \param[in] coarse_space
       * A \transient reference to the coarse-mesh trial-space to be used.
       *
       * \param[in] cubature_name
       * The name of the cubature rule to be used for integration of the mass matrices.
       */
      template<
        typename Matrix_,
        typename Vector_,
        typename FineSpace_,
        typename CoarseSpace_>
      static void assemble_prolongation(
        Matrix_& matrix,
        Vector_& vector,
        const FineSpace_& fine_space,
        const CoarseSpace_& coarse_space,
        const String& cubature_name)
      {
        Cubature::DynamicFactory cubature_factory(cubature_name);
        assemble_prolongation(matrix, vector, fine_space, coarse_space, cubature_factory);
      }

      /**
       * \brief Assembles a prolongation matrix and its corresponding weight vector.
       *
       * To obtain the final prolongation matrix, one needs to invert the weight vector
       * component-wise and scale the matrix rows by the inverted weights afterwards.
       * This can be accomplished by the \c component_invert and \c scale_rows operations
       * of the vector and matrix containers, resp.
       *
       * \param[in,out] matrix
       * A \transient reference to the prolongation matrix that is to be assembled.
       *
       * \param[in,out] vector
       * A \transient reference to the weight vector for the prolongation matrix.
       *
       * \param[in] fine_space
       * A \transient reference to the fine-mesh test-space to be used.
       *
       * \param[in] coarse_space
       * A \transient reference to the coarse-mesh trial-space to be used.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       */
      template<
        typename Matrix_,
        typename Vector_,
        typename FineSpace_,
        typename CoarseSpace_>
      static void assemble_prolongation(
        Matrix_& matrix,
        Vector_& vector,
        const FineSpace_& fine_space,
        const CoarseSpace_& coarse_space,
        const Cubature::DynamicFactory& cubature_factory)
      {
        // validate matrix and vector dimensions
        XASSERTM(matrix.rows() == fine_space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.columns() == coarse_space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(vector.size() == fine_space.get_num_dofs(), "invalid vector size");

        typedef typename Matrix_::DataType DataType;

        // typedefs for trafos, mesh and shape
        typedef typename FineSpace_::TrafoType FineTrafoType;
        typedef typename CoarseSpace_::TrafoType CoarseTrafoType;
        typedef typename CoarseSpace_::ShapeType ShapeType;

        // typedefs for dof-mappings
        typedef typename FineSpace_::DofMappingType FineDofMapping;
        typedef typename CoarseSpace_::DofMappingType CoarseDofMapping;

        // typedefs for trafo evaluators
        typedef typename FineTrafoType::template Evaluator<ShapeType, DataType>::Type FineTrafoEvaluator;
        typedef typename CoarseTrafoType::template Evaluator<ShapeType, DataType>::Type CoarseTrafoEvaluator;

        // typedefs for space evaluators
        typedef typename FineSpace_::template Evaluator<FineTrafoEvaluator>::Type FineSpaceEvaluator;
        typedef typename CoarseSpace_::template Evaluator<CoarseTrafoEvaluator>::Type CoarseSpaceEvaluator;

        // define fine and coarse mesh trafo configurations
        typedef typename FineSpaceEvaluator::template ConfigTraits<SpaceTags::value> FineSpaceConfigTraits;
        typedef typename CoarseSpaceEvaluator::template ConfigTraits<SpaceTags::value> CoarseSpaceConfigTraits;
        static constexpr TrafoTags fine_trafo_config = TrafoTags::jac_det | FineSpaceConfigTraits::trafo_config;
        static constexpr TrafoTags coarse_trafo_config = TrafoTags::jac_det | CoarseSpaceConfigTraits::trafo_config;

        // typedefs for trafo data
        typedef typename FineTrafoEvaluator::template ConfigTraits<fine_trafo_config>::EvalDataType FineTrafoEvalData;
        typedef typename CoarseTrafoEvaluator::template ConfigTraits<coarse_trafo_config>::EvalDataType CoarseTrafoEvalData;

        // typedef for space data
        typedef typename FineSpaceEvaluator::template ConfigTraits<SpaceTags::value>::EvalDataType FineSpaceEvalData;
        typedef typename CoarseSpaceEvaluator::template ConfigTraits<SpaceTags::value>::EvalDataType CoarseSpaceEvalData;

        // typedefs for cubature rule and factory
        typedef typename Intern::CubatureTraits<FineTrafoEvaluator>::RuleType CubatureRuleType;

        // fetch the trafos
        const FineTrafoType& fine_trafo = fine_space.get_trafo();
        const CoarseTrafoType& coarse_trafo = coarse_space.get_trafo();

        // get fine mesh element inverse permutation
        const Adjacency::Permutation& coarse_perm =
          coarse_trafo.get_mesh().get_mesh_permutation().get_perm();
        const Adjacency::Permutation& fine_perm =
          fine_trafo.get_mesh().get_mesh_permutation().get_inv_perm();

        // create the cubature rules
        CubatureRuleType fine_cubature(Cubature::ctor_factory, cubature_factory);
        CubatureRuleType refine_cubature;
        Cubature::RefineFactoryCore::create(refine_cubature, fine_cubature);

        // helper struct to calculate fine mesh cell index
        const Geometry::Intern::CoarseFineCellMapping<
          typename FineSpace_::MeshType, typename CoarseSpace_::MeshType>
          cfmapping(fine_trafo.get_mesh(), coarse_trafo.get_mesh());

        // OpenMP parallel region
        FEAT_PRAGMA_OMP(parallel)
        {
          FineTrafoEvalData fine_trafo_data;
          CoarseTrafoEvalData coarse_trafo_data;

          FineSpaceEvalData fine_space_data;
          CoarseSpaceEvalData coarse_space_data;

          // create matrix scatter-axpy
          typename Matrix_::ScatterAxpy scatter_maxpy(matrix);
          typename Vector_::ScatterAxpy scatter_vaxpy(vector);

          // create DOF-mappings
          FineDofMapping fine_dof_mapping(fine_space);
          CoarseDofMapping coarse_dof_mapping(coarse_space);

          // create trafo evaluators
          FineTrafoEvaluator fine_trafo_eval(fine_trafo);
          CoarseTrafoEvaluator coarse_trafo_eval(coarse_trafo);

          // create space evaluators
          FineSpaceEvaluator fine_space_eval(fine_space);
          CoarseSpaceEvaluator coarse_space_eval(coarse_space);

          // allocate fine-mesh mass matrix
          Tiny::Matrix<DataType, FineSpaceEvaluator::max_local_dofs, FineSpaceEvaluator::max_local_dofs> mass;

          // allocate local matrix data for interlevel-mesh mass matrix
          Tiny::Matrix<DataType, FineSpaceEvaluator::max_local_dofs, CoarseSpaceEvaluator::max_local_dofs> lmd, lid;

          // allocate local vector data for weight vector
          Tiny::Vector<DataType, FineSpaceEvaluator::max_local_dofs> lvd;

          // pivot array for factorization
          int pivot[FineSpaceEvaluator::max_local_dofs];

          // loop over all coarse mesh cells
          FEAT_PRAGMA_OMP(for)
          for(Index ccell = 0; ccell < coarse_trafo_eval.get_num_cells(); ++ccell)
          {
            // prepare coarse trafo evaluator
            coarse_trafo_eval.prepare(ccell);

            // prepare coarse space evaluator
            coarse_space_eval.prepare(coarse_trafo_eval);

            // fetch number of local coarse DOFs
            const int coarse_num_loc_dofs = coarse_space_eval.get_num_local_dofs();

            // prepare coarse mesh dof-mapping
            coarse_dof_mapping.prepare(ccell);

            // get coarse cell index with respect to 2-level ordering
            const Index ccell_2lvl = (coarse_perm.empty() ? ccell : coarse_perm.map(ccell));

            // loop over all child cells
            for(Index child(0); child < cfmapping.get_num_children(); ++child)
            {
              // calculate fine mesh cell index with respect to 2 level ordering
              const Index fcell_2lvl = cfmapping.calc_fcell(ccell_2lvl, child);

              // get fine cell index with respect to potential permutation
              const Index fcell = (fine_perm.empty() ? fcell_2lvl : fine_perm.map(fcell_2lvl));

              // prepare fine trafo evaluator
              fine_trafo_eval.prepare(fcell);

              // prepare fine space evaluator
              fine_space_eval.prepare(fine_trafo_eval);

              // fetch number of local fine DOFs
              const int fine_num_loc_dofs = fine_space_eval.get_num_local_dofs();

              // format local matrices
              mass.format();
              lmd.format();
              lvd.format();

              // loop over all cubature points and integrate
              for(int k(0); k < fine_cubature.get_num_points(); ++k)
              {
                // compute coarse mesh cubature point index
                const int l(int(child) * fine_cubature.get_num_points() + k);

                // compute trafo data
                fine_trafo_eval(fine_trafo_data, fine_cubature.get_point(k));
                coarse_trafo_eval(coarse_trafo_data, refine_cubature.get_point(l));

                // compute basis function data
                fine_space_eval(fine_space_data, fine_trafo_data);
                coarse_space_eval(coarse_space_data, coarse_trafo_data);

                // fine mesh test function loop
                for(int i(0); i < fine_num_loc_dofs; ++i)
                {
                  // fine mesh trial function loop
                  for(int j(0); j < fine_num_loc_dofs; ++j)
                  {
                    mass(i,j) += fine_trafo_data.jac_det * fine_cubature.get_weight(k) *
                      fine_space_data.phi[i].value * fine_space_data.phi[j].value;
                    // go for next fine mesh trial DOF
                  }

                  // coarse mesh trial function loop
                  for(int j(0); j < coarse_num_loc_dofs; ++j)
                  {
                    lmd(i,j) +=
                      fine_trafo_data.jac_det * fine_cubature.get_weight(k) *
                      fine_space_data.phi[i].value * coarse_space_data.phi[j].value;
                    // go for next fine mesh trial DOF
                  }
                  // go for next fine mesh test DOF
                }
                // go for next cubature point
              }

              // finish coarse mesh evaluators
              fine_space_eval.finish();
              fine_trafo_eval.finish();

              // invert fine mesh mass matrix
              Math::invert_matrix(fine_num_loc_dofs, mass.sn, &mass.v[0][0], pivot);

              // Note:
              // Usually, one would check whether the determinant returned by the invert_matrix
              // function is normal. However, this can lead to false alerts when assembling in
              // single precision, as the mass matrix entries are of magnitude h^2 (in 2D), i.e.
              // the determinant can become subnormal or even (numerically) zero although the
              // condition number of the matrix is still fine and the inversion was successful.
              // Therefore, we first multiply the (hopefully) inverted mass matrix by the
              // inter-level mass matrix and check whether the Frobenius norm of the result
              // is normal. If our matrix inversion failed, the result is virtually guaranteed
              // to be garbage, so this should serve well enough as a sanity check.

              // compute X := M^{-1}*N
              lid.set_mat_mat_mult(mass, lmd);

              // sanity check for matrix inversion
              if(!Math::isnormal(lid.norm_frobenius()))
                throw LocalMassMatrixSingularException();

              // prepare fine mesh dof-mapping
              fine_dof_mapping.prepare(fcell);
              lvd.format(DataType(1));

              FEAT_PRAGMA_OMP(critical)
              {
                // incorporate local matrix
                scatter_maxpy(lid, fine_dof_mapping, coarse_dof_mapping);

                // update weights
                scatter_vaxpy(lvd, fine_dof_mapping);
              }

              // finish fine mesh dof-mapping
              fine_dof_mapping.finish();

              // go for next child cell
            }

            // finish coarse mesh evaluators
            coarse_space_eval.finish();
            coarse_trafo_eval.finish();

            // finish coarse mesh dof-mapping
            coarse_dof_mapping.finish();

            // go for next coarse mesh cell
          } // FEAT_PRAGMA_OMP(for)
        } // FEAT_PRAGMA_OMP(parallel)
      }

      /**
       * \brief Assembles a prolongation matrix.
       *
       * \attention
       * This function <b>must not</b> be used to assemble prolongation matrices for
       * parallel (i.e. global) simulations, as it will be scaled incorrectly due to
       * missing weight synchronization!\n
       * Use this function only in serial simulations!
       *
       * \param[in,out] matrix
       * A \transient reference to the prolongation matrix that is to be assembled.
       *
       * \param[in] fine_space
       * A \transient reference to the fine-mesh test-space to be used.
       *
       * \param[in] coarse_space
       * A \transient reference to the coarse-mesh trial-space to be used.
       *
       * \param[in] cubature_name
       * The name of the cubature rule to be used for integration.
       */
      template<
        typename Matrix_,
        typename FineSpace_,
        typename CoarseSpace_>
        static void assemble_prolongation_direct(
          Matrix_& matrix,
          const FineSpace_& fine_space,
          const CoarseSpace_& coarse_space,
          const String& cubature_name)
      {
        Cubature::DynamicFactory cubature_factory(cubature_name);
        assemble_prolongation_direct(matrix, fine_space, coarse_space, cubature_factory);
      }

      /**
       * \brief Assembles a prolongation matrix.
       *
       * \attention
       * This function <b>must not</b> be used to assemble prolongation matrices for
       * parallel (i.e. global) simulations, as it will be scaled incorrectly due to
       * missing weight synchronization!\n
       * Use this function only in serial simulations!
       *
       * \param[in,out] matrix
       * A \transient reference to the prolongation matrix that is to be assembled.
       *
       * \param[in] fine_space
       * A \transient reference to the fine-mesh test-space to be used.
       *
       * \param[in] coarse_space
       * A \transient reference to the coarse-mesh trial-space to be used.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       */
      template<
        typename Matrix_,
        typename FineSpace_,
        typename CoarseSpace_>
      static void assemble_prolongation_direct(
        Matrix_& matrix,
        const FineSpace_& fine_space,
        const CoarseSpace_& coarse_space,
        const Cubature::DynamicFactory& cubature_factory)
      {
        // create a weight vector
        auto weight = matrix.create_vector_l();
        matrix.format();
        weight.format();

        // assemble matrix and weight
        assemble_prolongation(matrix, weight, fine_space, coarse_space, cubature_factory);

        // scale prolongation matrix rows by inverse weights
        weight.component_invert(weight);
        matrix.scale_rows(matrix, weight);
      }

      /**
       * \brief Assembles a truncation matrix and its corresponding weight vector.
       *
       * To obtain the final truncation matrix, one needs to invert the weight vector
       * component-wise and scale the matrix rows by the inverted weights afterwards.
       * This can be accomplished by the \c component_invert and \c scale_rows operations
       * of the vector and matrix containers, resp.
       *
       * \param[in,out] matrix
       * A \transient reference to the truncation matrix that is to be assembled.
       *
       * \param[in,out] vector
       * A \transient reference to the weight vector for the reduction matrix.
       *
       * \param[in] fine_space
       * A \transient reference to the fine-mesh test-space to be used.
       *
       * \param[in] coarse_space
       * A \transient reference to the coarse-mesh trial-space to be used.
       *
       * \param[in] cubature_name
       * The name of the cubature rule to be used for integration.
       */
      template<
        typename Matrix_,
        typename Vector_,
        typename FineSpace_,
        typename CoarseSpace_>
      static void assemble_truncation(
        Matrix_& matrix,
        Vector_& vector,
        const FineSpace_& fine_space,
        const CoarseSpace_& coarse_space,
        const String& cubature_name)
      {
        Cubature::DynamicFactory cubature_factory(cubature_name);
        assemble_truncation(matrix, vector, fine_space, coarse_space, cubature_factory);
      }

      /**
       * \brief Assembles a truncation matrix and its corresponding weight vector.
       *
       * To obtain the final truncation matrix, one needs to invert the weight vector
       * component-wise and scale the matrix rows by the inverted weights afterwards.
       * This can be accomplished by the \c component_invert and \c scale_rows operations
       * of the vector and matrix containers, resp.
       *
       * \param[in,out] matrix
       * A \transient reference to the truncation matrix that is to be assembled.
       *
       * \param[in,out] vector
       * A \transient reference to the weight vector for the reduction matrix.
       *
       * \param[in] fine_space
       * A \transient reference to the fine-mesh test-space to be used.
       *
       * \param[in] coarse_space
       * A \transient reference to the coarse-mesh trial-space to be used.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       */
      template<
        typename Matrix_,
        typename Vector_,
        typename FineSpace_,
        typename CoarseSpace_>
      static void assemble_truncation(
        Matrix_& matrix,
        Vector_& vector,
        const FineSpace_& fine_space,
        const CoarseSpace_& coarse_space,
        const Cubature::DynamicFactory& cubature_factory)
      {
        // validate matrix and vector dimensions
        XASSERTM(matrix.rows() == coarse_space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.columns() == fine_space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(vector.size() == coarse_space.get_num_dofs(), "invalid vector size");

        typedef typename Matrix_::DataType DataType;

        // typedefs for trafos, mesh and shape
        typedef typename FineSpace_::TrafoType FineTrafoType;
        typedef typename CoarseSpace_::TrafoType CoarseTrafoType;
        typedef typename CoarseSpace_::ShapeType ShapeType;

        // typedefs for dof-mappings
        typedef typename FineSpace_::DofMappingType FineDofMapping;
        typedef typename CoarseSpace_::DofMappingType CoarseDofMapping;

        // typedefs for trafo evaluators
        typedef typename FineTrafoType::template Evaluator<ShapeType, DataType>::Type FineTrafoEvaluator;
        typedef typename CoarseTrafoType::template Evaluator<ShapeType, DataType>::Type CoarseTrafoEvaluator;

        // typedefs for space evaluators
        typedef typename FineSpace_::template Evaluator<FineTrafoEvaluator>::Type FineSpaceEvaluator;
        typedef typename CoarseSpace_::template Evaluator<CoarseTrafoEvaluator>::Type CoarseSpaceEvaluator;

        // define fine and coarse mesh trafo configurations
        typedef typename FineSpaceEvaluator::template ConfigTraits<SpaceTags::value> FineSpaceConfigTraits;
        typedef typename CoarseSpaceEvaluator::template ConfigTraits<SpaceTags::value> CoarseSpaceConfigTraits;
        static constexpr TrafoTags fine_trafo_config = TrafoTags::jac_det | FineSpaceConfigTraits::trafo_config;
        static constexpr TrafoTags coarse_trafo_config = TrafoTags::jac_det | CoarseSpaceConfigTraits::trafo_config;

        // typedefs for trafo data
        typedef typename FineTrafoEvaluator::template ConfigTraits<fine_trafo_config>::EvalDataType FineTrafoEvalData;
        typedef typename CoarseTrafoEvaluator::template ConfigTraits<coarse_trafo_config>::EvalDataType CoarseTrafoEvalData;

        // typedef for space data
        typedef typename FineSpaceEvaluator::template ConfigTraits<SpaceTags::value>::EvalDataType FineSpaceEvalData;
        typedef typename CoarseSpaceEvaluator::template ConfigTraits<SpaceTags::value>::EvalDataType CoarseSpaceEvalData;

        // typedefs for cubature rule and factory
        typedef typename Intern::CubatureTraits<FineTrafoEvaluator>::RuleType CubatureRuleType;

        // fetch the trafos
        const FineTrafoType& fine_trafo = fine_space.get_trafo();
        const CoarseTrafoType& coarse_trafo = coarse_space.get_trafo();

        // helper struct to calculate fine mesh cell index
        const Geometry::Intern::CoarseFineCellMapping<
          typename FineSpace_::MeshType, typename CoarseSpace_::MeshType>
          cfmapping(fine_trafo.get_mesh(), coarse_trafo.get_mesh());

        // get fine mesh element inverse permutation
        const Adjacency::Permutation& coarse_perm =
          coarse_trafo.get_mesh().get_mesh_permutation().get_perm();
        const Adjacency::Permutation& fine_perm =
          fine_trafo.get_mesh().get_mesh_permutation().get_inv_perm();

        // create the cubature rules
        CubatureRuleType coarse_cubature(Cubature::ctor_factory, cubature_factory);
        CubatureRuleType fine_cubature(Cubature::ctor_factory, cubature_factory);
        CubatureRuleType refine_cubature;
        Cubature::RefineFactoryCore::create(refine_cubature, fine_cubature);

        FEAT_PRAGMA_OMP(parallel)
        {
          FineTrafoEvalData fine_trafo_data;
          CoarseTrafoEvalData coarse_trafo_data;
          FineSpaceEvalData fine_space_data;
          CoarseSpaceEvalData coarse_space_data;

          // create matrix scatter-axpy
          typename Matrix_::ScatterAxpy scatter_maxpy(matrix);
          typename Vector_::ScatterAxpy scatter_vaxpy(vector);

          // create DOF-mappings
          FineDofMapping fine_dof_mapping(fine_space);
          CoarseDofMapping coarse_dof_mapping(coarse_space);

          // create trafo evaluators
          FineTrafoEvaluator fine_trafo_eval(fine_trafo);
          CoarseTrafoEvaluator coarse_trafo_eval(coarse_trafo);

          // create space evaluators
          FineSpaceEvaluator fine_space_eval(fine_space);
          CoarseSpaceEvaluator coarse_space_eval(coarse_space);

          // allocate fine-mesh mass matrix
          Tiny::Matrix<DataType, CoarseSpaceEvaluator::max_local_dofs, CoarseSpaceEvaluator::max_local_dofs> mass;

          // allocate local matrix data for interlevel-mesh mass matrix
          Tiny::Matrix<DataType, CoarseSpaceEvaluator::max_local_dofs, FineSpaceEvaluator::max_local_dofs> lmd, lid;

          // allocate local vector data for weight vector
          Tiny::Vector<DataType, CoarseSpaceEvaluator::max_local_dofs> lvd;

          // pivot array for factorization
          int pivot[CoarseSpaceEvaluator::max_local_dofs];

          // loop over all coarse mesh cells
          FEAT_PRAGMA_OMP(for)
          for(Index ccell = 0; ccell < coarse_trafo_eval.get_num_cells(); ++ccell)
          {
            // prepare coarse trafo evaluator
            coarse_trafo_eval.prepare(ccell);

            // prepare coarse space evaluator
            coarse_space_eval.prepare(coarse_trafo_eval);

            // fetch number of local coarse DOFs
            const int coarse_num_loc_dofs = coarse_space_eval.get_num_local_dofs();

            // prepare coarse mesh dof-mapping
            coarse_dof_mapping.prepare(ccell);

            // Let's assemble the coarse-mesh mass matrix first
            mass.format();
            for(int k(0); k < coarse_cubature.get_num_points(); ++k)
            {
              // compute trafo and space data
              coarse_trafo_eval(coarse_trafo_data, coarse_cubature.get_point(k));
              coarse_space_eval(coarse_space_data, coarse_trafo_data);

              // coarse mesh test function loop
              for(int i(0); i < coarse_num_loc_dofs; ++i)
              {
                // coarse mesh trial function loop
                for(int j(0); j < coarse_num_loc_dofs; ++j)
                {
                  mass(i,j) += coarse_trafo_data.jac_det * coarse_cubature.get_weight(k) *
                    coarse_space_data.phi[i].value * coarse_space_data.phi[j].value;
                }
              }
            }

            // invert coarse mesh mass matrix
            Math::invert_matrix(coarse_num_loc_dofs, mass.sn, &mass.v[0][0], pivot);

            // get coarse cell index with respect to 2-level ordering
            const Index ccell_2lvl = (coarse_perm.empty() ? ccell : coarse_perm.map(ccell));

            // loop over all child cells
            for(Index child(0); child < cfmapping.get_num_children(); ++child)
            {
              // calculate fine mesh cell index with respect to 2 level ordering
              const Index fcell_2lvl = cfmapping.calc_fcell(ccell_2lvl, child);

              // get fine cell index with respect to potential permutation
              const Index fcell = (fine_perm.empty() ? fcell_2lvl : fine_perm.map(fcell_2lvl));

              // prepare fine trafo evaluator
              fine_trafo_eval.prepare(fcell);

              // prepare fine space evaluator
              fine_space_eval.prepare(fine_trafo_eval);

              // fetch number of local fine DOFs
              const int fine_num_loc_dofs = fine_space_eval.get_num_local_dofs();

              // format local matrices
              lmd.format();

              // loop over all cubature points and integrate
              for(int k(0); k < fine_cubature.get_num_points(); ++k)
              {
                // compute coarse mesh cubature point index
                const int l(int(child) * fine_cubature.get_num_points() + k);

                // compute trafo data
                fine_trafo_eval(fine_trafo_data, fine_cubature.get_point(k));
                coarse_trafo_eval(coarse_trafo_data, refine_cubature.get_point(l));

                // compute basis function data
                fine_space_eval(fine_space_data, fine_trafo_data);
                coarse_space_eval(coarse_space_data, coarse_trafo_data);

                // coarse mesh test function loop
                for(int i(0); i < coarse_num_loc_dofs; ++i)
                {
                  // fine mesh trial function loop
                  for(int j(0); j < fine_num_loc_dofs; ++j)
                  {
                    lmd(i,j) +=
                      fine_trafo_data.jac_det * fine_cubature.get_weight(k) *
                      coarse_space_data.phi[i].value * fine_space_data.phi[j].value;
                    // go for next fine mesh trial DOF
                  }
                  // go for next coarse mesh test DOF
                }
                // go for next cubature point
              }

              // finish coarse mesh evaluators
              fine_space_eval.finish();
              fine_trafo_eval.finish();

              // compute X := M^{-1}*N
              lid.set_mat_mat_mult(mass, lmd);

              // sanity check for matrix inversion
              if(!Math::isnormal(lid.norm_frobenius()))
                throw LocalMassMatrixSingularException();

              // prepare fine mesh dof-mapping
              fine_dof_mapping.prepare(fcell);

              // incorporate local matrix
              FEAT_PRAGMA_OMP(critical)
              {
                scatter_maxpy(lid, coarse_dof_mapping, fine_dof_mapping);
              }

              // finish fine mesh dof-mapping
              fine_dof_mapping.finish();

              // go for next child cell
            }

            // update weights
            lvd.format(DataType(1));
            FEAT_PRAGMA_OMP(critical)
            {
              scatter_vaxpy(lvd, coarse_dof_mapping);
            }

            // finish coarse mesh dof-mapping
            coarse_dof_mapping.finish();

            // finish coarse mesh evaluators
            coarse_space_eval.finish();
            coarse_trafo_eval.finish();

            // go for next coarse mesh cell
          } // FEAT_PRAGMA_OMP(for)
        } // FEAT_PRAGMA_OMP(parallel)
      }

      /**
       * \brief Assembles a truncation matrix.
       *
       * \attention
       * This function <b>must not</b> be used to assemble truncation matrices for
       * parallel (i.e. global) simulations, as it will be scaled incorrectly due to
       * missing weight synchronization!\n
       * Use this function only in serial simulations!
       *
       * \param[in,out] matrix
       * A \transient reference to the truncation matrix that is to be assembled.
       *
       * \param[in] fine_space
       * A \transient reference to the fine-mesh test-space to be used.
       *
       * \param[in] coarse_space
       * A \transient reference to the coarse-mesh trial-space to be used.
       *
       * \param[in] cubature_name
       * The name of the cubature rule to be used for integration.
       */
      template<
        typename Matrix_,
        typename FineSpace_,
        typename CoarseSpace_>
      static void assemble_truncation_direct(
        Matrix_& matrix,
        const FineSpace_& fine_space,
        const CoarseSpace_& coarse_space,
        const String& cubature_name)
      {
        Cubature::DynamicFactory cubature_factory(cubature_name);
        assemble_truncation_direct(matrix, fine_space, coarse_space, cubature_factory);
      }

      /**
       * \brief Assembles a truncation matrix.
       *
       * \attention
       * This function <b>must not</b> be used to assemble truncation matrices for
       * parallel (i.e. global) simulations, as it will be scaled incorrectly due to
       * missing weight synchronization!\n
       * Use this function only in serial simulations!
       *
       * \param[in,out] matrix
       * A \transient reference to the truncation matrix that is to be assembled.
       *
       * \param[in] fine_space
       * A \transient reference to the fine-mesh test-space to be used.
       *
       * \param[in] coarse_space
       * A \transient reference to the coarse-mesh trial-space to be used.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       */
      template<
        typename Matrix_,
        typename FineSpace_,
        typename CoarseSpace_>
      static void assemble_truncation_direct(
        Matrix_& matrix,
        const FineSpace_& fine_space,
        const CoarseSpace_& coarse_space,
        const Cubature::DynamicFactory& cubature_factory)
      {
        // create a weight vector
        auto weight = matrix.create_vector_l();
        matrix.format();
        weight.format();

        // assemble matrix and weight
        assemble_truncation(matrix, weight, fine_space, coarse_space, cubature_factory);

        // scale truncation matrix rows by inverse weights
        weight.component_invert(weight);
        matrix.scale_rows(matrix, weight);
      }

      /**
       * \brief Assembles a generic inter-mesh transfer matrix and its corresponding weight vector.
       *
       * To obtain the final inter-mesh transfer matrix, one needs to invert the weight vector
       * component-wise and scale the matrix rows by the inverted weights afterwards.
       * This can be accomplished by the \c component_invert and \c scale_rows operations
       * of the vector and matrix containers, resp.
       *
       * \note This function is parallelized using OpenMP.
       *
       * \attention
       * It is strongly recommended to use an open cubature rule, i.e. all points of the cubature
       * rule shall be in the interior of the reference cell. This reduces the risk of point
       * unmapping failures
       *
       * \param[in,out] matrix
       * A \transient reference to the transfer matrix that is to be assembled.
       *
       * \param[in,out] vector
       * A \transient reference to the weight vector for the transfer matrix.
       *
       * \param[in] space_target
       * A \transient reference to the target-mesh (test) space to be used.
       *
       * \param[in] space_source
       * A \transient reference to the source-mesh (trial) space to be used.
       *
       * \param[in] target_to_source
       * A \transient reference to a target-to-source mesh adjactor. If the source and/or target meshes
       * are permuted, then the adjactor is assumed to be also permuted consistently.
       *
       * \param[in] cubature_name
       * The name of the cubature rule to be used for integration of the mass matrices.
       * This should be an open cubature rule, i.e. all point shall be in the interior of the
       * reference cell to avoid unmapping errors.
       *
       * \returns
       * The number of cubature points that failed to be unmapped onto the source mesh.
       * If the return value is 0, then all points were unmapped successfully and therefore the
       * transfer matrix is assumed to be valid, otherwise the transfer matrix \e may be invalid.
       */
      template<
        typename Matrix_,
        typename Vector_,
        typename TargetSpace_,
        typename SourceSpace_,
        typename Target2SourceAdjactor_>
      static int assemble_intermesh_transfer(
        Matrix_& matrix,
        Vector_& vector,
        const TargetSpace_& space_target,
        const SourceSpace_& space_source,
        const Target2SourceAdjactor_& trg2src,
        const String& cubature_name)
      {
        // validate matrix and vector dimensions
        XASSERTM(matrix.rows() == space_target.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.columns() == space_source.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(vector.size() == space_target.get_num_dofs(), "invalid vector size");

        typedef typename Matrix_::DataType DataType;

        // typedefs for trafos, mesh and shape
        typedef typename TargetSpace_::TrafoType TargetTrafoType;
        typedef typename SourceSpace_::TrafoType SourceTrafoType;
        typedef typename SourceSpace_::ShapeType ShapeType;

        // typedefs for dof-mappings
        typedef typename TargetSpace_::DofMappingType TargetDofMapping;
        typedef typename SourceSpace_::DofMappingType SourceDofMapping;

        // typedefs for trafo evaluators
        typedef typename TargetTrafoType::template Evaluator<ShapeType, DataType>::Type TargetTrafoEvaluator;
        typedef typename SourceTrafoType::template Evaluator<ShapeType, DataType>::Type SourceTrafoEvaluator;

        // typedefs for space evaluators
        typedef typename TargetSpace_::template Evaluator<TargetTrafoEvaluator>::Type TargetSpaceEvaluator;
        typedef typename SourceSpace_::template Evaluator<SourceTrafoEvaluator>::Type SourceSpaceEvaluator;

        // define target and source mesh trafo configurations
        typedef typename TargetSpaceEvaluator::template ConfigTraits<SpaceTags::value> TargetSpaceConfigTraits;
        typedef typename SourceSpaceEvaluator::template ConfigTraits<SpaceTags::value> SourceSpaceConfigTraits;
        static constexpr TrafoTags target_trafo_config = TrafoTags::img_point | TrafoTags::jac_det | TargetSpaceConfigTraits::trafo_config;
        static constexpr TrafoTags source_trafo_config = TrafoTags::jac_det | SourceSpaceConfigTraits::trafo_config;

        // typedefs for trafo data
        typedef typename TargetTrafoEvaluator::template ConfigTraits<target_trafo_config>::EvalDataType TargetTrafoEvalData;
        typedef typename SourceTrafoEvaluator::template ConfigTraits<source_trafo_config>::EvalDataType SourceTrafoEvalData;

        // typedef for space data
        typedef typename TargetSpaceEvaluator::template ConfigTraits<SpaceTags::value>::EvalDataType TargetSpaceEvalData;
        typedef typename SourceSpaceEvaluator::template ConfigTraits<SpaceTags::value>::EvalDataType SourceSpaceEvalData;

        // typedef for source inverse mapping data
        typedef Trafo::InverseMapping<SourceTrafoType, DataType> SourceInverseMappingType;
        typedef typename SourceInverseMappingType::InvMapDataType SourceInvMapData;

        // typedefs for cubature rule and factory
        typedef typename Intern::CubatureTraits<TargetTrafoEvaluator>::RuleType CubatureRuleType;

        // local target mass matrix
        typedef Tiny::Matrix<DataType, TargetSpaceEvaluator::max_local_dofs, TargetSpaceEvaluator::max_local_dofs> LocTrgMatrixType;

        // local target-source inter-mesh matrix
        typedef Tiny::Matrix<DataType, TargetSpaceEvaluator::max_local_dofs, SourceSpaceEvaluator::max_local_dofs> LocTSIMatrixType;

        // fetch the trafos
        const TargetTrafoType& target_trafo = space_target.get_trafo();
        const SourceTrafoType& source_trafo = space_source.get_trafo();

        // create the cubature rule
        Cubature::DynamicFactory cubature_factory(cubature_name);
        CubatureRuleType target_cubature(Cubature::ctor_factory, cubature_factory);
        const std::size_t num_cub_pts = std::size_t(target_cubature.get_num_points());

        // number of failed cubature points
        int failed_points = 0;

        // we open a OpenMP section here since this is a quite expensive assembly
        FEAT_PRAGMA_OMP(parallel reduction(+:failed_points))
        {
          // create matrix scatter-axpy
          typename Matrix_::ScatterAxpy scatter_maxpy(matrix);
          typename Vector_::ScatterAxpy scatter_vaxpy(vector);

          // create DOF-mappings
          TargetDofMapping target_dof_mapping(space_target);
          SourceDofMapping source_dof_mapping(space_source);

          // create trafo evaluators
          TargetTrafoEvaluator target_trafo_eval(target_trafo);
          SourceTrafoEvaluator source_trafo_eval(source_trafo);

          // create space evaluators
          TargetSpaceEvaluator space_target_eval(space_target);
          SourceSpaceEvaluator space_source_eval(space_source);

          TargetTrafoEvalData target_trafo_data;
          SourceTrafoEvalData source_trafo_data;
          TargetSpaceEvalData space_target_data;
          SourceSpaceEvalData space_source_data;

          // target space basis values, pre-multiplied by cubature weights and jacobian determinants
          typedef Tiny::Vector<DataType, TargetSpaceEvaluator::max_local_dofs> TargetBasisValueVector;
          std::vector<TargetBasisValueVector> target_basis_values;
          target_basis_values.resize(num_cub_pts);

          // create an inverse mapping object for the source
          SourceInverseMappingType source_inv_map(source_trafo);

          // the inverse mapping data object
          std::vector<SourceInvMapData> source_inv_map_data;
          source_inv_map_data.resize(num_cub_pts);

          // [local source cell index] = (inverse map cell index, cubature point index)
          std::vector<std::vector<std::pair<std::size_t, std::size_t>>> source_cell_map;
          source_cell_map.reserve(30);

          // allocate target-mesh mass matrix
          LocTrgMatrixType trg_mass;

          // allocate local matrix data for inter-mesh mass matrix
          LocTSIMatrixType tsi_mass;

          // local transfer matrix
          std::vector<LocTSIMatrixType> loc_trans;
          loc_trans.reserve(30);

          // allocate local vector data for weight vector and set all entries to 1
          Tiny::Vector<DataType, TargetSpaceEvaluator::max_local_dofs> lvd;
          lvd.format(DataType(1));

          // pivot array for factorization
          int pivot[TargetSpaceEvaluator::max_local_dofs];

          // loop over all target mesh cells
          // note: OpenMP (sometimes) requires a signed iteration variable
          FEAT_PRAGMA_OMP(for schedule(dynamic,16))
          for(std::int64_t itrg_cell = 0ull; itrg_cell < std::int64_t(target_trafo_eval.get_num_cells()); ++itrg_cell)
          {
            // cast target cell index
            const Index trg_cell = Index(itrg_cell);

            // collect all source cells that intersect with the current target cell
            std::vector<Index> source_cells;
            for(auto it = trg2src.image_begin(trg_cell); it != trg2src.image_end(trg_cell); ++it)
              source_cells.push_back(*it);

            // allocate source cell map to correct size
            source_cell_map.resize(source_cells.size());

            // allocate local transfer matrices
            loc_trans.resize(source_cells.size());

            // prepare target trafo evaluator
            target_trafo_eval.prepare(trg_cell);

            // prepare target space evaluator
            space_target_eval.prepare(target_trafo_eval);

            // fetch number of local target DOFs
            const int target_num_loc_dofs = space_target_eval.get_num_local_dofs();

            // format local matrices
            trg_mass.format();

            // first of all, let's loop over all target cubature points and unmap them onto the source mesh
            // we will also integrate the target mesh mass matrix while we're at it
            for(std::size_t k(0); k < num_cub_pts; ++k)
            {
              // compute trafo data
              target_trafo_eval(target_trafo_data, target_cubature.get_point(int(k)));

              // compute basis function data
              space_target_eval(space_target_data, target_trafo_data);

              // integrate target mesh mass matrix entries
              for(int i(0); i < target_num_loc_dofs; ++i)
              {
                // save pre-weighted target basis value for later use
                target_basis_values.at(k)[i] = target_trafo_data.jac_det * target_cubature.get_weight(int(k)) * space_target_data.phi[i].value;
                for(int j(0); j < target_num_loc_dofs; ++j)
                {
                  trg_mass(i,j) += target_basis_values.at(k)[i] * space_target_data.phi[j].value;
                }
              }

              // try to unmap cubature point from real coordinates onto source cell reference coordinates
              SourceInvMapData& simd = source_inv_map_data.at(k);
              simd = source_inv_map.unmap_point(target_trafo_data.img_point, source_cells, false);

              // did we fail to unmap the point?
              if(simd.empty())
              {
                ++failed_points;
                continue;
              }

              // now let's loop over all source cells that were found to contain/intersect the current cubature
              // point and update our source-cell-to-cubature-point map
              for(std::size_t i(0); i < simd.size(); ++i)
              {
                // find the local source cell index that corresponds to this source cell index
                std::size_t loc_idx = 0u;
                for(; loc_idx < source_cells.size(); ++loc_idx)
                  if(source_cells.at(loc_idx) == simd.cells.at(i))
                    break;
                XASSERT(loc_idx < source_cells.size());

                // add this cubature point to the set of points for this source cell
                source_cell_map.at(loc_idx).push_back(std::make_pair(i, k));
              }
            }

            // finish target mesh evaluators
            space_target_eval.finish();
            target_trafo_eval.finish();

            // invert target mesh mass matrix
            Math::invert_matrix(target_num_loc_dofs, trg_mass.sn, &trg_mass.v[0][0], pivot);

            // Note:
            // Usually, one would check whether the determinant returned by the invert_matrix
            // function is normal. However, this can lead to false alerts when assembling in
            // single precision, as the mass matrix entries are of magnitude h^2 (in 2D), i.e.
            // the determinant can become subnormal or even (numerically) zero although the
            // condition number of the matrix is still fine and the inversion was successful.
            // Therefore, we first multiply the (hopefully) inverted mass matrix by the
            // inter-level mass matrix and check whether the Frobenius norm of the result
            // is normal. If our matrix inversion failed, the result is virtually guaranteed
            // to be garbage, so this should serve well enough as a sanity check.

            // now let's loop over all source cells that were found
            for(std::size_t is = 0u; is < source_cell_map.size(); ++is)
            {
              // get source cell index
              const Index src_cell = source_cells.at(is);

              // get the cubature points map
              std::vector<std::pair<std::size_t, std::size_t>> cub_pts_map = source_cell_map.at(is);

              // no cubature points on this source cell?
              if(cub_pts_map.empty())
                continue;

              // prepare source trafo evaluator
              source_trafo_eval.prepare(src_cell);

              // prepare source space evaluator
              space_source_eval.prepare(source_trafo_eval);

              // fetch number of local source DOFs
              const int source_num_loc_dofs = space_source_eval.get_num_local_dofs();

              // format our inter-mesh mass matrix
              tsi_mass.format();

              // loop over all cubature points and integrate
              for(auto cit = cub_pts_map.begin(); cit != cub_pts_map.end(); ++cit)
              {
                // get the index of the source cell in the inverse mapping data object
                const std::size_t src_i = cit->first;

                // get the index of the cubature point
                const std::size_t cub_k = cit->second;

                // get the source inverse mapping data for this cubature point
                SourceInvMapData& simd = source_inv_map_data.at(cub_k);

                // ensure that this is the right cell
                XASSERT(simd.cells.at(src_i) == src_cell);

                // compute our averaging weight: we have to normalize in case that this cubature
                // point is intersecting more than one source cell
                DataType weight = DataType(1) / DataType(simd.size());

                // alright, evaluate our source trafo and space
                source_trafo_eval(source_trafo_data, simd.dom_points.at(src_i));

                // compute basis function data
                space_source_eval(space_source_data, source_trafo_data);

                // integrate inter-mesh mass matrix
                for(int i(0); i < target_num_loc_dofs; ++i)
                  for(int j(0); j < source_num_loc_dofs; ++j)
                    tsi_mass(i, j) += weight * target_basis_values.at(cub_k)[i] * space_source_data.phi[j].value;

              } // go for next cubature point

              // finish source mesh evaluators
              space_source_eval.finish();
              source_trafo_eval.finish();

              // compute X := M^{-1}*N
              loc_trans.at(is).set_mat_mat_mult(trg_mass, tsi_mass);

              // sanity check for matrix inversion
              if(!Math::isnormal(loc_trans.at(is).norm_frobenius()))
                throw LocalMassMatrixSingularException();

            } // continue with next source cell

            // prepare target mesh dof-mapping
            target_dof_mapping.prepare(trg_cell);

            // scattering is a potential data race
            FEAT_PRAGMA_OMP(critical)
            {
              // loop over all intersecting source cells
              for(std::size_t is = 0u; is < source_cells.size(); ++is)
              {
                // no cubature points on this source cell?
                if(source_cell_map.at(is).empty())
                  continue;

                // prepare source mesh dof-mapping
                source_dof_mapping.prepare(source_cells.at(is));

                // incorporate local matrix
                scatter_maxpy(loc_trans.at(is), target_dof_mapping, source_dof_mapping);

                // finish source mesh dof-mapping
                source_dof_mapping.finish();

                // clear the cubature point map for this source cell
                source_cell_map.at(is).clear();
              }

              // update target weights
              scatter_vaxpy(lvd, target_dof_mapping);
            } // pragma omp critical

            // finish target mesh dof-mapping
            target_dof_mapping.finish();
          } // go for next target mesh cell
        } // pragma omp parallel

        // return number of failed points
        return failed_points;
      }

      /**
       * \brief Assembles a generic inter-mesh transfer matrix and its corresponding weight vector.
       *
       * \attention
       * This function <b>must not</b> be used to assemble transfer matrices for parallel (i.e. global) simulations,
       * as it will be scaled incorrectly due to missing weight synchronization!\n
       * Use this function only in serial simulations!
       *
       * \note This function is parallelized using OpenMP.
       *
       * \param[in,out] matrix
       * A \transient reference to the transfer matrix that is to be assembled.
       *
       * \param[in] space_target
       * A \transient reference to the target-mesh (test) space to be used.
       *
       * \param[in] space_source
       * A \transient reference to the source-mesh (trial) space to be used.
       *
       * \param[in] target_to_source
       * A \transient reference to a target-to-source mesh adjactor. If the source and/or target meshes
       * are permuted, then the adjactor is assumed to be also permuted consistently.
       *
       * \param[in] cubature_name
       * The name of the cubature rule to be used for integration of the mass matrices.
       *
       * \returns
       * The number of cubature points that failed to be unmapped onto the source mesh.
       * If the return value is 0, then all points were unmapped successfully and therefore the
       * transfer matrix is assumed to be valid, otherwise the transfer matrix \e may be invalid.
       */
      template<
        typename Matrix_,
        typename TargetSpace_,
        typename SourceSpace_,
        typename Target2SourceAdjactor_>
      static int assemble_intermesh_transfer_direct(
        Matrix_& matrix,
        const TargetSpace_& space_target,
        const SourceSpace_& space_source,
        const Target2SourceAdjactor_& trg2src,
        const String& cubature_name)
      {
        // create a weight vector
        auto weight = matrix.create_vector_l();
        matrix.format();
        weight.format();

        // assemble matrix and weight
        int rtn = assemble_intermesh_transfer(matrix, weight, space_target, space_source, trg2src, cubature_name);

        // scale transfer matrix rows by inverse weights
        weight.component_invert(weight);
        matrix.scale_rows(matrix, weight);

        return rtn;
      }

      /**
       * \brief Prolongates a primal vector and assembles a compatible weight vector
       *
       * To obtain the final prolongated vector, one needs to invert the weight vector component-wise
       * and scale the prolongated vector component-wise by the inverted weights afterwards.
       * This can be accomplished by the \c component_invert and \c component_product
       * operations of the vector container, resp.
       *
       * \attention
       * In the case of global vectors, the weight vector has to be synchronized via sync_0 before
       * the component-wise inversion and scaling to obtain the correctly scaled primal vector.
       *
       * \param[in,out] vector_f
       * A \transient reference to the fine-mesh vector that is to be assembled
       * Is assumed to be allocated and formatted to 0.
       *
       * \param[in,out] vector_w
       * A \transient reference to the fine-mesh weight vector that is to be assembled.
       * Is assumed to be allocated and formatted to 0.
       *
       * \param[in] vector_c
       * A \transient reference to the coarse-mesh vector that is to be prolongated.
       *
       * \param[in] fine_space
       * A \transient reference to the fine-mesh test-space to be used.
       *
       * \param[in] coarse_space
       * A \transient reference to the coarse-mesh trial-space to be used.
       *
       * \param[in] cubature_name
       * The name of the cubature rule to be used for integration.
       */
      template<
        typename Vector_,
        typename FineSpace_,
        typename CoarseSpace_>
      static void prolongate_vector(
        Vector_& vector_f,
        Vector_& vector_w,
        const Vector_& vector_c,
        const FineSpace_& fine_space,
        const CoarseSpace_& coarse_space,
        const String& cubature_name)
      {
        // validate matrix and vector dimensions
        XASSERTM(vector_f.size() == vector_w.size(), "invalid vector sizes");
        XASSERTM(vector_f.size() == fine_space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(vector_c.size() == coarse_space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(vector_f.size() == fine_space.get_num_dofs(), "invalid vector size");

        typedef typename Vector_::DataType DataType;
        typedef typename Vector_::ValueType ValueType;

        // typedefs for trafos, mesh and shape
        typedef typename FineSpace_::TrafoType FineTrafoType;
        typedef typename CoarseSpace_::TrafoType CoarseTrafoType;
        typedef typename CoarseSpace_::ShapeType ShapeType;

        // typedefs for dof-mappings
        typedef typename FineSpace_::DofMappingType FineDofMapping;
        typedef typename CoarseSpace_::DofMappingType CoarseDofMapping;

        // typedefs for trafo evaluators
        typedef typename FineTrafoType::template Evaluator<ShapeType, DataType>::Type FineTrafoEvaluator;
        typedef typename CoarseTrafoType::template Evaluator<ShapeType, DataType>::Type CoarseTrafoEvaluator;

        // typedefs for space evaluators
        typedef typename FineSpace_::template Evaluator<FineTrafoEvaluator>::Type FineSpaceEvaluator;
        typedef typename CoarseSpace_::template Evaluator<CoarseTrafoEvaluator>::Type CoarseSpaceEvaluator;

        // define fine and coarse mesh trafo configurations
        typedef typename FineSpaceEvaluator::template ConfigTraits<SpaceTags::value> FineSpaceConfigTraits;
        typedef typename CoarseSpaceEvaluator::template ConfigTraits<SpaceTags::value> CoarseSpaceConfigTraits;
        static constexpr TrafoTags fine_trafo_config = TrafoTags::jac_det | FineSpaceConfigTraits::trafo_config;
        static constexpr TrafoTags coarse_trafo_config = TrafoTags::jac_det | CoarseSpaceConfigTraits::trafo_config;

        // typedefs for trafo data
        typedef typename FineTrafoEvaluator::template ConfigTraits<fine_trafo_config>::EvalDataType FineTrafoEvalData;
        typedef typename CoarseTrafoEvaluator::template ConfigTraits<coarse_trafo_config>::EvalDataType CoarseTrafoEvalData;
        FineTrafoEvalData fine_trafo_data;
        CoarseTrafoEvalData coarse_trafo_data;

        // typedef for space data
        typedef typename FineSpaceEvaluator::template ConfigTraits<SpaceTags::value>::EvalDataType FineSpaceEvalData;
        typedef typename CoarseSpaceEvaluator::template ConfigTraits<SpaceTags::value>::EvalDataType CoarseSpaceEvalData;
        FineSpaceEvalData fine_space_data;
        CoarseSpaceEvalData coarse_space_data;

        // typedefs for cubature rule and factory
        typedef typename Assembly::Intern::CubatureTraits<FineTrafoEvaluator>::RuleType CubatureRuleType;

        // create gather/scatter-axpys
        typename Vector_::GatherAxpy gather_c(vector_c);
        typename Vector_::ScatterAxpy scatter_f(vector_f);
        typename Vector_::ScatterAxpy scatter_w(vector_w);

        // create DOF-mappings
        FineDofMapping fine_dof_mapping(fine_space);
        CoarseDofMapping coarse_dof_mapping(coarse_space);

        // fetch the trafos
        const FineTrafoType& fine_trafo = fine_space.get_trafo();
        const CoarseTrafoType& coarse_trafo = coarse_space.get_trafo();

        // create trafo evaluators
        FineTrafoEvaluator fine_trafo_eval(fine_trafo);
        CoarseTrafoEvaluator coarse_trafo_eval(coarse_trafo);

        // create space evaluators
        FineSpaceEvaluator fine_space_eval(fine_space);
        CoarseSpaceEvaluator coarse_space_eval(coarse_space);

        // create the cubature rules
        Cubature::DynamicFactory cubature_factory(cubature_name);
        CubatureRuleType fine_cubature(Cubature::ctor_factory, cubature_factory);
        CubatureRuleType refine_cubature;
        Cubature::RefineFactoryCore::create(refine_cubature, fine_cubature);

        // allocate fine-mesh mass matrix
        Tiny::Matrix<DataType, FineSpaceEvaluator::max_local_dofs, FineSpaceEvaluator::max_local_dofs> mass;

        // allocate local matrix data for interlevel-mesh mass matrix
        Tiny::Matrix<DataType, FineSpaceEvaluator::max_local_dofs, CoarseSpaceEvaluator::max_local_dofs> lmd, lid;

        // allocate local vector data for weight/coarse/fine vector
        Tiny::Vector<ValueType, FineSpaceEvaluator::max_local_dofs> lv_w, lv_c, lv_f;

        // pivot array for factorization
        int pivot[FineSpaceEvaluator::max_local_dofs];

        // helper struct to calculate fine mesh cell index
        const Geometry::Intern::CoarseFineCellMapping<
          typename FineSpace_::MeshType, typename CoarseSpace_::MeshType>
          cfmapping(fine_trafo.get_mesh(), coarse_trafo.get_mesh());

        // get fine mesh element inverse permutation
        const Adjacency::Permutation& coarse_perm =
          coarse_trafo.get_mesh().get_mesh_permutation().get_perm();
        const Adjacency::Permutation& fine_perm =
          fine_trafo.get_mesh().get_mesh_permutation().get_inv_perm();

        // loop over all coarse mesh cells
        for(Index ccell(0); ccell < coarse_trafo_eval.get_num_cells(); ++ccell)
        {
          // prepare coarse trafo evaluator
          coarse_trafo_eval.prepare(ccell);

          // prepare coarse space evaluator
          coarse_space_eval.prepare(coarse_trafo_eval);

          // fetch number of local coarse DOFs
          const int coarse_num_loc_dofs = coarse_space_eval.get_num_local_dofs();

          // prepare coarse mesh dof-mapping
          coarse_dof_mapping.prepare(ccell);

          // gather local coarse vector
          lv_c.format();
          gather_c(lv_c, coarse_dof_mapping);

          // get coarse cell index with respect to 2-level ordering
          const Index ccell_2lvl = (coarse_perm.empty() ? ccell : coarse_perm.map(ccell));

          // loop over all child cells
          for(Index child(0); child < cfmapping.get_num_children(); ++child)
          {
            // calculate fine mesh cell index with respect to 2 level ordering
            const Index fcell_2lvl = cfmapping.calc_fcell(ccell_2lvl, child);

            // get fine cell index with respect to potential permutation
            const Index fcell = (fine_perm.empty() ? fcell_2lvl : fine_perm.map(fcell_2lvl));

            // prepare fine trafo evaluator
            fine_trafo_eval.prepare(fcell);

            // prepare fine space evaluator
            fine_space_eval.prepare(fine_trafo_eval);

            // fetch number of local fine DOFs
            const int fine_num_loc_dofs = fine_space_eval.get_num_local_dofs();

            // format local matrices
            mass.format();
            lmd.format();
            lv_w.format();
            lv_f.format();

            // loop over all cubature points and integrate
            for(int k(0); k < fine_cubature.get_num_points(); ++k)
            {
              // compute coarse mesh cubature point index
              const int l(int(child) * fine_cubature.get_num_points() + k);

              // compute trafo data
              fine_trafo_eval(fine_trafo_data, fine_cubature.get_point(k));
              coarse_trafo_eval(coarse_trafo_data, refine_cubature.get_point(l));

              // compute basis function data
              fine_space_eval(fine_space_data, fine_trafo_data);
              coarse_space_eval(coarse_space_data, coarse_trafo_data);

              // fine mesh test function loop
              for(int i(0); i < fine_num_loc_dofs; ++i)
              {
                // fine mesh trial function loop
                for(int j(0); j < fine_num_loc_dofs; ++j)
                {
                  mass(i,j) += fine_trafo_data.jac_det * fine_cubature.get_weight(k) *
                    fine_space_data.phi[i].value * fine_space_data.phi[j].value;
                  // go for next fine mesh trial DOF
                }

                // coarse mesh trial function loop
                for(int j(0); j < coarse_num_loc_dofs; ++j)
                {
                  lmd(i,j) +=
                    fine_trafo_data.jac_det * fine_cubature.get_weight(k) *
                    fine_space_data.phi[i].value * coarse_space_data.phi[j].value;
                  // go for next fine mesh trial DOF
                }
                // go for next fine mesh test DOF
              }
              // go for next cubature point
            }

            // finish coarse mesh evaluators
            fine_space_eval.finish();
            fine_trafo_eval.finish();

            // invert fine mesh mass matrix
            Math::invert_matrix(fine_num_loc_dofs, mass.sn, &mass.v[0][0], pivot);

            // Note:
            // Usually, one would check whether the determinant returned by the invert_matrix
            // function is normal. However, this can lead to false alerts when assembling in
            // single precision, as the mass matrix entries are of magnitude h^2 (in 2D), i.e.
            // the determinant can become subnormal or even (numerically) zero although the
            // condition number of the matrix is still fine and the inversion was successful.
            // Therefore, we first multiply the (hopefully) inverted mass matrix by the
            // inter-level mass matrix and check whether the Frobenius norm of the result
            // is normal. If our matrix inversion failed, the result is virtually guaranteed
            // to be garbage, so this should serve well enough as a sanity check.

            // compute X := M^{-1}*N
            lid.set_mat_mat_mult(mass, lmd);

            // sanity check for matrix inversion
            if(!Math::isnormal(lid.norm_frobenius()))
              throw LocalMassMatrixSingularException();

            // compute local fine vector
            lv_f.format();
            for(int i(0); i < fine_num_loc_dofs; ++i)
              for(int j(0); j < coarse_num_loc_dofs; ++j)
                lv_f[i] += lid[i][j] * lv_c[j];

            // prepare fine mesh dof-mapping
            fine_dof_mapping.prepare(fcell);

            // scatter local fine vector
            scatter_f(lv_f, fine_dof_mapping);

            // update weights
            lv_w.format(DataType(1));
            scatter_w(lv_w, fine_dof_mapping);

            // finish fine mesh dof-mapping
            fine_dof_mapping.finish();

            // go for next child cell
          }

          // finish coarse mesh evaluators
          coarse_space_eval.finish();
          coarse_trafo_eval.finish();

          // finish coarse mesh dof-mapping
          coarse_dof_mapping.finish();

          // go for next coarse mesh cell
        }
      }

      /**
       * \brief Prolongates a primal vector directly
       *
       * \attention
       * This function <b>must not</b> be used to prolongate vectors for parallel (i.e. global)
       * simulations, as it will be scaled incorrectly due to missing weight vector synchronization!
       *
       * \param[in,out] vector_f
       * A \transient reference to the fine-mesh vector that is to be assembled
       * Is assumed to be allocated and formatted to 0.
       *
       * \param[in] vector_c
       * A \transient reference to the coarse-mesh vector that is to be prolongated.
       *
       * \param[in] fine_space
       * A \transient reference to the fine-mesh test-space to be used.
       *
       * \param[in] coarse_space
       * A \transient reference to the coarse-mesh trial-space to be used.
       *
       * \param[in] cubature_name
       * The name of the cubature rule to be used for integration.
       */
      template<
        typename Vector_,
        typename FineSpace_,
        typename CoarseSpace_>
      static void prolongate_vector_direct(
        Vector_& vector_f,
        const Vector_& vector_c,
        const FineSpace_& fine_space,
        const CoarseSpace_& coarse_space,
        const String& cubature_name)
      {
        Vector_ vector_w = vector_f.clone(LAFEM::CloneMode::Layout);
        vector_w.format();

        prolongate_vector(vector_f, vector_w, vector_c, fine_space, coarse_space, cubature_name);

        // finally, scale fine mesh vector by inverse weights
        vector_w.component_invert(vector_w);
        vector_f.component_product(vector_f, vector_w);
      }

      /**
       * \brief Performs a generic inter-mesh transfer of a primal vector and assembles a compatible weight vector
       *
       * To obtain the final transferred vector, one needs to invert the weight vector component-wise
       * and scale the prolongated vector component-wise by the inverted weights afterwards.
       * This can be accomplished by the \c component_invert and \c component_product
       * operations of the vector container, resp.
       *
       * \note This function is parallelized using OpenMP.
       *
       * \attention
       * It is strongly recommended to use an open cubature rule, i.e. all points of the cubature
       * rule shall be in the interior of the reference cell. This reduces the risk of point
       * unmapping failures
       *
       * \param[in,out] vector_target
       * A \transient reference to the target vector that is to be assembled.
       * Is assumed to be allocated and formatted to 0.
       *
       * \param[in,out] vector_weight
       * A \transient reference to the weight vector for the transfer matrix.
       * Is assumed to be allocated and formatted to 0.
       *
       * \param[in] vector_source
       * A \transient reference to the source vector that is to be transferred.
       *
       * \param[in] space_target
       * A \transient reference to the target-mesh (test) space to be used.
       *
       * \param[in] space_source
       * A \transient reference to the source-mesh (trial) space to be used.
       *
       * \param[in] target_to_source
       * A \transient reference to a target-to-source mesh adjactor. If the source and/or target meshes
       * are permuted, then the adjactor is assumed to be also permuted consistently.
       *
       * \param[in] cubature_name
       * The name of the cubature rule to be used for integration of the mass matrices.
       * This should be an open cubature rule, i.e. all point shall be in the interior of the
       * reference cell to avoid unmapping errors.
       *
       * \returns
       * The number of cubature points that failed to be unmapped onto the source mesh.
       * If the return value is 0, then all points were unmapped successfully and therefore the
       * transfer matrix is assumed to be valid, otherwise the transfer matrix \e may be invalid.
       */
      template<
        typename Vector_,
        typename TargetSpace_,
        typename SourceSpace_,
        typename Target2SourceAdjactor_>
      static int transfer_intermesh_vector(
        Vector_& vector_target,
        Vector_& vector_weight,
        const Vector_& vector_source,
        const TargetSpace_& space_target,
        const SourceSpace_& space_source,
        const Target2SourceAdjactor_& trg2src,
        const String& cubature_name)
      {
        // validate matrix and vector dimensions
        XASSERTM(vector_target.size() == space_target.get_num_dofs(), "invalid vector size");
        XASSERTM(vector_weight.size() == space_target.get_num_dofs(), "invalid vector size");
        XASSERTM(vector_source.size() == space_source.get_num_dofs(), "invalid vector size");

        typedef typename Vector_::DataType DataType;
        typedef typename Vector_::ValueType ValueType;

        // typedefs for trafos, mesh and shape
        typedef typename TargetSpace_::TrafoType TargetTrafoType;
        typedef typename SourceSpace_::TrafoType SourceTrafoType;
        typedef typename SourceSpace_::ShapeType ShapeType;

        // typedefs for dof-mappings
        typedef typename TargetSpace_::DofMappingType TargetDofMapping;
        typedef typename SourceSpace_::DofMappingType SourceDofMapping;

        // typedefs for trafo evaluators
        typedef typename TargetTrafoType::template Evaluator<ShapeType, DataType>::Type TargetTrafoEvaluator;
        typedef typename SourceTrafoType::template Evaluator<ShapeType, DataType>::Type SourceTrafoEvaluator;

        // typedefs for space evaluators
        typedef typename TargetSpace_::template Evaluator<TargetTrafoEvaluator>::Type TargetSpaceEvaluator;
        typedef typename SourceSpace_::template Evaluator<SourceTrafoEvaluator>::Type SourceSpaceEvaluator;

        // define target and source mesh trafo configurations
        typedef typename TargetSpaceEvaluator::template ConfigTraits<SpaceTags::value> TargetSpaceConfigTraits;
        typedef typename SourceSpaceEvaluator::template ConfigTraits<SpaceTags::value> SourceSpaceConfigTraits;
        static constexpr TrafoTags target_trafo_config = TrafoTags::img_point | TrafoTags::jac_det | TargetSpaceConfigTraits::trafo_config;
        static constexpr TrafoTags source_trafo_config = TrafoTags::jac_det | SourceSpaceConfigTraits::trafo_config;

        // typedefs for trafo data
        typedef typename TargetTrafoEvaluator::template ConfigTraits<target_trafo_config>::EvalDataType TargetTrafoEvalData;
        typedef typename SourceTrafoEvaluator::template ConfigTraits<source_trafo_config>::EvalDataType SourceTrafoEvalData;

        // typedef for space data
        typedef typename TargetSpaceEvaluator::template ConfigTraits<SpaceTags::value>::EvalDataType TargetSpaceEvalData;
        typedef typename SourceSpaceEvaluator::template ConfigTraits<SpaceTags::value>::EvalDataType SourceSpaceEvalData;

        // typedef for source inverse mapping data
        typedef Trafo::InverseMapping<SourceTrafoType, DataType> SourceInverseMappingType;
        typedef typename SourceInverseMappingType::InvMapDataType SourceInvMapData;

        // typedefs for cubature rule and factory
        typedef typename Intern::CubatureTraits<TargetTrafoEvaluator>::RuleType CubatureRuleType;

        // local target mass matrix
        typedef Tiny::Matrix<DataType, TargetSpaceEvaluator::max_local_dofs, TargetSpaceEvaluator::max_local_dofs> LocTrgMatrixType;

        // local target-source inter-mesh matrix
        typedef Tiny::Matrix<DataType, TargetSpaceEvaluator::max_local_dofs, SourceSpaceEvaluator::max_local_dofs> LocTSIMatrixType;

        // local target and source vectors
        typedef Tiny::Vector<ValueType, TargetSpaceEvaluator::max_local_dofs> LocTrgVectorType;
        typedef Tiny::Vector<ValueType, SourceSpaceEvaluator::max_local_dofs> LocSrcVectorType;

        // fetch the trafos
        const TargetTrafoType& target_trafo = space_target.get_trafo();
        const SourceTrafoType& source_trafo = space_source.get_trafo();

        // create the cubature rule
        Cubature::DynamicFactory cubature_factory(cubature_name);
        CubatureRuleType target_cubature(Cubature::ctor_factory, cubature_factory);
        const std::size_t num_cub_pts = std::size_t(target_cubature.get_num_points());

        // number of failed cubature points
        int failed_points = 0;

        // we open a OpenMP section here since this is a quite expensive assembly
        FEAT_PRAGMA_OMP(parallel reduction(+:failed_points))
        {
          // create gather/scatter-axpys
          typename Vector_::GatherAxpy gather_s(vector_source);
          typename Vector_::ScatterAxpy scatter_t(vector_target);
          typename Vector_::ScatterAxpy scatter_w(vector_weight);

          // create DOF-mappings
          TargetDofMapping target_dof_mapping(space_target);
          SourceDofMapping source_dof_mapping(space_source);

          // create trafo evaluators
          TargetTrafoEvaluator target_trafo_eval(target_trafo);
          SourceTrafoEvaluator source_trafo_eval(source_trafo);

          // create space evaluators
          TargetSpaceEvaluator space_target_eval(space_target);
          SourceSpaceEvaluator space_source_eval(space_source);

          TargetTrafoEvalData target_trafo_data;
          SourceTrafoEvalData source_trafo_data;
          TargetSpaceEvalData space_target_data;
          SourceSpaceEvalData space_source_data;

          // target space basis values, pre-multiplied by cubature weights and jacobian determinants
          typedef Tiny::Vector<DataType, TargetSpaceEvaluator::max_local_dofs> TargetBasisValueVector;
          std::vector<TargetBasisValueVector> target_basis_values;
          target_basis_values.resize(num_cub_pts);

          // create an inverse mapping object for the source
          SourceInverseMappingType source_inv_map(source_trafo);

          // the inverse mapping data object
          std::vector<SourceInvMapData> source_inv_map_data;
          source_inv_map_data.resize(num_cub_pts);

          // [local source cell index] = (inverse map cell index, cubature point index)
          std::vector<std::vector<std::pair<std::size_t, std::size_t>>> source_cell_map;
          source_cell_map.reserve(30);

          // allocate local stuff
          LocTrgMatrixType trg_mass;
          LocTSIMatrixType tsi_mass, loc_trans;
          LocTrgVectorType loc_trg, loc_weight;
          LocSrcVectorType loc_src;

          // format local weights to 1
          loc_weight.format(DataType(1));

          // pivot array for factorization
          int pivot[TargetSpaceEvaluator::max_local_dofs];

          // loop over all target mesh cells
          // note: OpenMP (sometimes) requires a signed iteration variable
          FEAT_PRAGMA_OMP(for schedule(dynamic,16))
          for(std::int64_t itrg_cell = 0ull; itrg_cell < std::int64_t(target_trafo_eval.get_num_cells()); ++itrg_cell)
          {
            // cast target cell index
            const Index trg_cell = Index(itrg_cell);

            // collect all source cells that intersect with the current target cell
            std::vector<Index> source_cells;
            for(auto it = trg2src.image_begin(trg_cell); it != trg2src.image_end(trg_cell); ++it)
              source_cells.push_back(*it);

            // allocate source cell map to correct size
            source_cell_map.resize(source_cells.size());

            // allocate local transfer matrices
            //loc_trg_vec.resize(source_cells.size());
            loc_trg.format();

            // prepare target trafo evaluator
            target_trafo_eval.prepare(trg_cell);

            // prepare target space evaluator
            space_target_eval.prepare(target_trafo_eval);

            // fetch number of local target DOFs
            const int target_num_loc_dofs = space_target_eval.get_num_local_dofs();

            // format local matrices
            trg_mass.format();

            // first of all, let's loop over all target cubature points and unmap them onto the source mesh
            // we will also integrate the target mesh mass matrix while we're at it
            for(std::size_t k(0); k < num_cub_pts; ++k)
            {
              // compute trafo data
              target_trafo_eval(target_trafo_data, target_cubature.get_point(int(k)));

              // compute basis function data
              space_target_eval(space_target_data, target_trafo_data);

              // integrate target mesh mass matrix entries
              for(int i(0); i < target_num_loc_dofs; ++i)
              {
                // save pre-weighted target basis value for later use
                target_basis_values.at(k)[i] = target_trafo_data.jac_det * target_cubature.get_weight(int(k)) * space_target_data.phi[i].value;
                for(int j(0); j < target_num_loc_dofs; ++j)
                {
                  trg_mass(i,j) += target_basis_values.at(k)[i] * space_target_data.phi[j].value;
                }
              }

              // try to unmap cubature point from real coordinates onto source cell reference coordinates
              SourceInvMapData& simd = source_inv_map_data.at(k);
              simd = source_inv_map.unmap_point(target_trafo_data.img_point, source_cells, false);

              // did we fail to unmap the point?
              if(simd.empty())
              {
                ++failed_points;
                continue;
              }

              // now let's loop over all source cells that were found to contain/intersect the current cubature
              // point and update our source-cell-to-cubature-point map
              for(std::size_t i(0); i < simd.size(); ++i)
              {
                // find the local source cell index that corresponds to this source cell index
                std::size_t loc_idx = 0u;
                for(; loc_idx < source_cells.size(); ++loc_idx)
                  if(source_cells.at(loc_idx) == simd.cells.at(i))
                    break;
                XASSERT(loc_idx < source_cells.size());

                // add this cubature point to the set of points for this source cell
                source_cell_map.at(loc_idx).push_back(std::make_pair(i, k));
              }
            }

            // finish target mesh evaluators
            space_target_eval.finish();
            target_trafo_eval.finish();

            // invert target mesh mass matrix
            Math::invert_matrix(target_num_loc_dofs, trg_mass.sn, &trg_mass.v[0][0], pivot);

            // Note:
            // Usually, one would check whether the determinant returned by the invert_matrix
            // function is normal. However, this can lead to false alerts when assembling in
            // single precision, as the mass matrix entries are of magnitude h^2 (in 2D), i.e.
            // the determinant can become subnormal or even (numerically) zero although the
            // condition number of the matrix is still fine and the inversion was successful.
            // Therefore, we first multiply the (hopefully) inverted mass matrix by the
            // inter-level mass matrix and check whether the Frobenius norm of the result
            // is normal. If our matrix inversion failed, the result is virtually guaranteed
            // to be garbage, so this should serve well enough as a sanity check.

            // now let's loop over all source cells that were found
            for(std::size_t is = 0u; is < source_cell_map.size(); ++is)
            {
              // get source cell index
              const Index src_cell = source_cells.at(is);

              // get the cubature points map
              std::vector<std::pair<std::size_t, std::size_t>> cub_pts_map = source_cell_map.at(is);

              // no cubature points on this source cell?
              if(cub_pts_map.empty())
                continue;

              // prepare source trafo evaluator
              source_trafo_eval.prepare(src_cell);

              // prepare source space evaluator
              space_source_eval.prepare(source_trafo_eval);

              // fetch number of local source DOFs
              const int source_num_loc_dofs = space_source_eval.get_num_local_dofs();

              // format our inter-mesh mass matrix
              tsi_mass.format();

              // loop over all cubature points and integrate
              for(auto cit = cub_pts_map.begin(); cit != cub_pts_map.end(); ++cit)
              {
                // get the index of the source cell in the inverse mapping data object
                const std::size_t src_i = cit->first;

                // get the index of the cubature point
                const std::size_t cub_k = cit->second;

                // get the source inverse mapping data for this cubature point
                SourceInvMapData& simd = source_inv_map_data.at(cub_k);

                // ensure that this is the right cell
                XASSERT(simd.cells.at(src_i) == src_cell);

                // compute our averaging weight: we have to normalize in case that this cubature
                // point is intersecting more than one source cell
                DataType weight = DataType(1) / DataType(simd.size());

                // alright, evaluate our source trafo and space
                source_trafo_eval(source_trafo_data, simd.dom_points.at(src_i));

                // compute basis function data
                space_source_eval(space_source_data, source_trafo_data);

                // integrate inter-mesh mass matrix
                for(int i(0); i < target_num_loc_dofs; ++i)
                  for(int j(0); j < source_num_loc_dofs; ++j)
                    tsi_mass(i, j) += weight * target_basis_values.at(cub_k)[i] * space_source_data.phi[j].value;

              } // go for next cubature point

                // finish source mesh evaluators
              space_source_eval.finish();
              source_trafo_eval.finish();

              // compute X := M^{-1}*N
              loc_trans.set_mat_mat_mult(trg_mass, tsi_mass);

              // sanity check for matrix inversion
              if(!Math::isnormal(loc_trans.norm_frobenius()))
                throw LocalMassMatrixSingularException();

              // prepare source mesh dof-mapping
              source_dof_mapping.prepare(source_cells.at(is));

              // gather the local vector
              loc_src.format();
              gather_s(loc_src, source_dof_mapping);

              // finish source mesh dof-mapping
              source_dof_mapping.finish();

              // update local target vector
              for(int i(0); i < target_num_loc_dofs; ++i)
                for(int j(0); j < source_num_loc_dofs; ++j)
                  loc_trg[i] += loc_trans[i][j] * loc_src[j];

            } // continue with next source cell

              // prepare target mesh dof-mapping
            target_dof_mapping.prepare(trg_cell);

            // scattering is a potential data race
            FEAT_PRAGMA_OMP(critical)
            {
              // update target weights
              scatter_t(loc_trg, target_dof_mapping);
              scatter_w(loc_weight, target_dof_mapping);
            } // pragma omp critical

              // finish target mesh dof-mapping
            target_dof_mapping.finish();
          } // go for next target mesh cell
        } // pragma omp parallel

          // return number of failed points
        return failed_points;
      }

      /**
       * \brief Performs a generic inter-mesh transfer of a primal vector directly
       *
       * \attention
       * This function <b>must not</b> be used to transfer vectors for parallel (i.e. global) simulations,
       * as it will be scaled incorrectly due to missing weight synchronization!\n
       * Use this function only in serial simulations!
       *
       * \note This function is parallelized using OpenMP.
       *
       * \param[in,out] vector_target
       * A \transient reference to the target vector that is to be assembled.
       * Is assumed to be allocated and formatted to 0.
       *
       * \param[in] vector_source
       * A \transient reference to the source vector that is to be transferred.
       *
       * \param[in] space_target
       * A \transient reference to the target-mesh (test) space to be used.
       *
       * \param[in] space_source
       * A \transient reference to the source-mesh (trial) space to be used.
       *
       * \param[in] target_to_source
       * A \transient reference to a target-to-source mesh adjactor. If the source and/or target meshes
       * are permuted, then the adjactor is assumed to be also permuted consistently.
       *
       * \param[in] cubature_name
       * The name of the cubature rule to be used for integration of the mass matrices.
       *
       * \returns
       * The number of cubature points that failed to be unmapped onto the source mesh.
       * If the return value is 0, then all points were unmapped successfully and therefore the
       * transferred vector is assumed to be valid, otherwise the transferred vector \e may be invalid.
       */
      template<
        typename Vector_,
        typename TargetSpace_,
        typename SourceSpace_,
        typename Target2SourceAdjactor_>
      static int transfer_intermesh_vector_direct(
        Vector_& vector_target,
        const Vector_& vector_source,
        const TargetSpace_& space_target,
        const SourceSpace_& space_source,
        const Target2SourceAdjactor_& trg2src,
        const String& cubature_name)
      {
        // create a weight vector
        Vector_ weight = vector_target.clone(LAFEM::CloneMode::Layout);
        weight.format();

        // assemble matrix and weight
        int rtn = transfer_intermesh_vector_direct(vector_target, weight, vector_source, space_target, space_source, trg2src, cubature_name);

        // scale transfer matrix rows by inverse weights
        weight.component_invert(weight);
        vector_target.component_product(vector_target, weight);

        return rtn;
      }
    }; // class GridTransfer<...>
  } // namespace Assembly
} // namespace FEAT
