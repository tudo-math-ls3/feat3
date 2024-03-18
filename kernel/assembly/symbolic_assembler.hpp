// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ASSEMBLY_SYMBOLIC_ASSEMBLER_HPP
#define KERNEL_ASSEMBLY_SYMBOLIC_ASSEMBLER_HPP 1

// includes, FEAT
#include <kernel/util/assertion.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/space/dof_mapping_renderer.hpp>
#include <kernel/lafem/null_matrix.hpp>
#include <kernel/geometry/intern/coarse_fine_cell_mapping.hpp>

namespace FEAT
{
  namespace Assembly
  {
    /**
     * \brief Symbolic Matrix and Graph assembler class
     *
     * This class is used for the assembly of matrix sparsity patterns, which are either
     * assembled into an Adjacency::Graph structure or directly into a sparse matrix
     * container as e.g. LAFEM::SparseMatrixCSR.
     *
     * This class supports the assembly of the following sparsity patterns:
     * - Standard pattern for identical or different test-/trial-spaces defined on the same mesh.
     *   This pattern is used for standard operator matrices. The corresponding functions are
     *   assemble_matrix_std1 and assemble_matrix_std2.
     *
     * - Extended-facet pattern for identical or different test-/trial-spaces on the same mesh.
     *   This pattern is used for matrices which include operator coupling over neighbor elements,
     *   e.g. for jump-stabilization or discontinuous Galerkin methods. The corresponding functions
     *   are assemble_matrix_ext_facet1 and assemble_matrix_ext_facet2.
     *
     * - Extended-node pattern for identical or different test-/trial-spaces on the same mesh.
     *   This pattern is used for matrices which include operator coupling over neighbor nodes.
     *   The corresponding functions are assemble_matrix_ext_node1 and assemble_matrix_ext-node2.
     *
     * - Standard pattern for test-/trial-spaces defined on a 2-level-refined mesh pair.
     *   This pattern is used for prolongation and restriction matrices in multigrid methods.
     *   The corresponding function is assemble_matrix_2lvl.
     *
     * \author Peter Zajac
     */
    class SymbolicAssembler
    {
    public:
      /**
       * \brief Assembles the standard Dof-Adjacency graph for different test- and trial-spaces.
       *
       * \param[in] test_space, trial_space
       * The \transient test- and trial-spaces to be used for the assembly. Must be defined on the same mesh.
       *
       * \returns
       * The standard Dof-Adjacency graph of the test- and trial-space combination.
       */
      template<typename TestSpace_, typename TrialSpace_>
      static Adjacency::Graph assemble_graph_std2(const TestSpace_& test_space, const TrialSpace_& trial_space)
      {
        // render dof-graphs
        Adjacency::Graph test_dof_graph(Space::DofMappingRenderer::render(test_space));
        Adjacency::Graph trial_dof_graph(Space::DofMappingRenderer::render(trial_space));

        // check dimensions
        XASSERTM(test_dof_graph.get_num_nodes_domain() == trial_dof_graph.get_num_nodes_domain(), "invalid test-/trial-space pair");

        // render transposed test-dof-mapping
        Adjacency::Graph test_dof_support(Adjacency::RenderType::transpose, test_dof_graph);

        // render composite test-dof-mapping/trial-dof-support graph
        Adjacency::Graph dof_adjactor(Adjacency::RenderType::injectify_sorted, test_dof_support, trial_dof_graph);

        // return the graph
        return dof_adjactor;
      }

      /**
       * \brief Assembles the standard Dof-Adjacency graph for identical test- and trial-spaces.
       *
       * \param[in] space
       * The \transient space representing the test- and trial spaces to be used for the assembly.
       *
       * \returns
       * The standard Dof-Adjacency graph of the space.
       */
      template<typename Space_>
      static Adjacency::Graph assemble_graph_std1(const Space_& space)
      {
        // create dof-mapping
        Adjacency::Graph dof_graph(Space::DofMappingRenderer::render(space));

        // render transposed dof-mapping
        Adjacency::Graph dof_support(Adjacency::RenderType::transpose, dof_graph);

        // render composite dof-mapping/dof-support graph
        Adjacency::Graph dof_adjactor(Adjacency::RenderType::injectify_sorted, dof_support, dof_graph);

        // return the graph
        return dof_adjactor;
      }

      /**
       * \brief Assembles the extended-facet Dof-Adjacency graph for different test- and trial-spaces.
       *
       * \param[in] test_space, trial_space
       * The \transient test- and trial-spaces to be used for the assembly. Must be defined on the same mesh.
       *
       * \returns
       * The extended-facet Dof-Adjacency graph of the test- and trial-space combination.
       */
      template<typename TestSpace_, typename TrialSpace_>
      static Adjacency::Graph assemble_graph_ext_facet2(const TestSpace_& test_space, const TrialSpace_& trial_space)
      {
        // render dof-graphs
        Adjacency::Graph test_dof_graph(Space::DofMappingRenderer::render(test_space));
        Adjacency::Graph trial_dof_graph(Space::DofMappingRenderer::render(trial_space));

        // check dimensions
        XASSERTM(test_dof_graph.get_num_nodes_domain() == trial_dof_graph.get_num_nodes_domain(), "invalid test-/trial-space pair");

        // get the shape dimension
        typedef typename TestSpace_::ShapeType ShapeType;
        static constexpr Index shape_dim = Index(ShapeType::dimension);

        // get the facet index set
        const auto& facet_at_shape = test_space.get_trafo().get_mesh().template get_index_set<shape_dim, shape_dim-1>();

        // transpose to get the facet support
        Adjacency::Graph shape_at_facet(Adjacency::RenderType::transpose, facet_at_shape);

        // render transposed test-dof-mapping
        Adjacency::Graph test_dof_support(Adjacency::RenderType::transpose, test_dof_graph);

        // render extended test-dof-support
        Adjacency::Graph test_dof_ext_sup(Adjacency::RenderType::injectify, test_dof_support, facet_at_shape);

        // render extended trial-dof-mapping
        Adjacency::Graph trial_dof_ext_graph(Adjacency::RenderType::injectify, shape_at_facet, trial_dof_graph);

        // render composite test-dof-mapping/trial-dof-support graph
        Adjacency::Graph dof_adjactor(Adjacency::RenderType::injectify_sorted, test_dof_ext_sup, trial_dof_ext_graph);

        // return the graph
        return dof_adjactor;
      }

      /**
       * \brief Assembles the extended-facet Dof-Adjacency graph for identical test- and trial-spaces.
       *
       * \param[in] space
       * The \transient space representing the test- and trial spaces to be used for the assembly.
       *
       * \returns
       * The extended-facet Dof-Adjacency graph of the space.
       */
      template<typename Space_>
      static Adjacency::Graph assemble_graph_ext_facet1(const Space_& space)
      {
        // create dof-mapping
        Adjacency::Graph dof_graph(Space::DofMappingRenderer::render(space));

        // get the shape dimension
        typedef typename Space_::ShapeType ShapeType;
        static constexpr Index shape_dim = Index(ShapeType::dimension);

        // get the facet index set
        const auto& facet_at_shape = space.get_trafo().get_mesh().template get_index_set<shape_dim, shape_dim-1>();

        // transpose to get the facet support
        Adjacency::Graph shape_at_facet(Adjacency::RenderType::transpose, facet_at_shape);

        // render extended dof-mapping
        Adjacency::Graph dof_ext_graph(Adjacency::RenderType::injectify, shape_at_facet, dof_graph);

        // render transposed extended dof-mapping
        Adjacency::Graph dof_ext_sup(Adjacency::RenderType::transpose, dof_ext_graph);

        // render composite dof-mapping/dof-support graph
        Adjacency::Graph dof_adjactor(Adjacency::RenderType::injectify_sorted, dof_ext_sup, dof_ext_graph);

        // return the graph
        return dof_adjactor;
      }

      /**
       * \brief Assembles the extended-node Dof-Adjacency graph for different test- and trial-spaces.
       *
       * \param[in] test_space, trial_space
       * The \transient test- and trial-spaces to be used for the assembly. Must be defined on the same mesh.
       *
       * \returns
       * The extended-node Dof-Adjacency graph of the test- and trial-space combination.
       */
      template<typename TestSpace_, typename TrialSpace_>
      static Adjacency::Graph assemble_graph_ext_node2(const TestSpace_& test_space, const TrialSpace_& trial_space)
      {
        // render dof-graphs
        Adjacency::Graph test_dof_graph(Space::DofMappingRenderer::render(test_space));
        Adjacency::Graph trial_dof_graph(Space::DofMappingRenderer::render(trial_space));

        // check dimensions
        XASSERTM(test_dof_graph.get_num_nodes_domain() == trial_dof_graph.get_num_nodes_domain(), "invalid test-/trial-space pair");

        // get the shape dimension
        typedef typename TestSpace_::ShapeType ShapeType;
        static constexpr Index shape_dim = Index(ShapeType::dimension);

        // get the facet index set
        const auto& facet_at_shape = test_space.get_trafo().get_mesh().template get_index_set<shape_dim, 0>();

        // transpose to get the facet support
        Adjacency::Graph shape_at_facet(Adjacency::RenderType::transpose, facet_at_shape);

        // render transposed test-dof-mapping
        Adjacency::Graph test_dof_support(Adjacency::RenderType::transpose, test_dof_graph);

        // render extended test-dof-support
        Adjacency::Graph test_dof_ext_sup(Adjacency::RenderType::injectify, test_dof_support, facet_at_shape);

        // render extended trial-dof-mapping
        Adjacency::Graph trial_dof_ext_graph(Adjacency::RenderType::injectify, shape_at_facet, trial_dof_graph);

        // render composite test-dof-mapping/trial-dof-support graph
        Adjacency::Graph dof_adjactor(Adjacency::RenderType::injectify_sorted, test_dof_ext_sup, trial_dof_ext_graph);

        // return the graph
        return dof_adjactor;
      }

      /**
       * \brief Assembles the extended-node Dof-Adjacency graph for identical test- and trial-spaces.
       *
       * \param[in] space
       * The \transient space representing the test- and trial spaces to be used for the assembly.
       *
       * \returns
       * The extended-node Dof-Adjacency graph of the space.
       */
      template<typename Space_>
      static Adjacency::Graph assemble_graph_ext_node1(const Space_& space)
      {
        // create dof-mapping
        Adjacency::Graph dof_graph(Space::DofMappingRenderer::render(space));

        // get the shape dimension
        typedef typename Space_::ShapeType ShapeType;
        static constexpr Index shape_dim = Index(ShapeType::dimension);

        // get the facet index set
        const auto& facet_at_shape = space.get_trafo().get_mesh().template get_index_set<shape_dim, 0>();

        // transpose to get the facet support
        Adjacency::Graph shape_at_facet(Adjacency::RenderType::transpose, facet_at_shape);

        // render extended dof-mapping
        Adjacency::Graph dof_ext_graph(Adjacency::RenderType::injectify, shape_at_facet, dof_graph);

        // render transposed extended dof-mapping
        Adjacency::Graph dof_ext_sup(Adjacency::RenderType::transpose, dof_ext_graph);

        // render composite dof-mapping/dof-support graph
        Adjacency::Graph dof_adjactor(Adjacency::RenderType::injectify_sorted, dof_ext_sup, dof_ext_graph);

        // return the graph
        return dof_adjactor;
      }

      /**
       * \brief Assembles the diagonal Dof-Adjacency graph for identical test- and trial-spaces.
       *
       * \param[in] space
       * The \transient space representing the test- and trial spaces to be used for the assembly.
       *
       * \returns
       * The diagonal Dof-Adjacency graph of the space.
       */
      template<typename Space_>
      static Adjacency::Graph assemble_graph_diag(const Space_& space)
      {
        // get number of DOFs
        const Index n = space.get_num_dofs();

        // create a graph
        Adjacency::Graph graph(n, n, n);

        // get and fill the arrays
        Index* dom_ptr = graph.get_domain_ptr();
        Index* img_idx = graph.get_image_idx();
        for(Index i(0); i < n; ++i)
          dom_ptr[i] = img_idx[i] = i;
        dom_ptr[n] = n;

        // return the graph
        return graph;
      }

      /**
      * \brief Assembles the Dof-Adjacency graph where test and trial spaces are defined on different meshes.
      *
      * \param[in] test_space
      * The \transient test-space defined on the test space mesh.
      *
      * \param[in] trial_space
      * The \transient trial-space defined on the trial space mesh.
      *
      * \param[in] trial2test_adjactor
      * A \transient object implementing the adjactor interface, which maps from the set of trial space mesh cells
      * to the test space mesh cells, which iterates over the set of all test space mesh cells intersecting with
      * a given trial space mesh cell.
      *
      * \returns
      * The inter-mesh Dof-Adjacency graph for the test- and trial-space combination defined on different meshes.
      */
      template<typename TestSpace_, typename TrialSpace_, typename TrialToTestAdjator_>
      static Adjacency::Graph assemble_graph_intermesh(
        const TestSpace_& test_space,
        const TrialSpace_& trial_space,
        const TrialToTestAdjator_& trial2test_adjactor)
      {
        // create test- and trial-dof-mappers
        Adjacency::Graph test_dof_mapping(Space::DofMappingRenderer::render(test_space));
        Adjacency::Graph trial_dof_mapping(Space::DofMappingRenderer::render(trial_space));

        // get trial and test mesh permutations
        const Adjacency::Permutation& trial_perm = trial_space.get_trafo().get_mesh().get_mesh_permutation().get_perm();
        const Adjacency::Permutation& test_inv_perm = test_space.get_trafo().get_mesh().get_mesh_permutation().get_inv_perm();
        if(!trial_perm.empty() || !test_inv_perm.empty())
        {
          // render trial2test adjactor
          Adjacency::Graph trial2test_graph(Adjacency::RenderType::as_is, trial2test_adjactor);

          // permute graph indices
          Adjacency::Graph permuted_graph;
          if(trial_perm.empty())
          {
            // no trial permutation, only test mesh permuted
            permuted_graph = std::move(trial2test_graph);
            permuted_graph.permute_indices(test_inv_perm);
          }
          else if(test_inv_perm.empty())
          {
            // no test permutation, only trial mesh permuted
            // create a dummy identity permutation for the test mesh
            Adjacency::Permutation id_perm(trial2test_graph.get_num_nodes_image(), Adjacency::Permutation::type_identity);
            permuted_graph = Adjacency::Graph(trial2test_graph, trial_perm, id_perm);
          }
          else
          {
            // both trial and test permutations exist
            permuted_graph = Adjacency::Graph(trial2test_graph, trial_perm, test_inv_perm);
          }

          // render transposed test-dof-mapping
          Adjacency::Graph test_dof_support(Adjacency::RenderType::injectify_transpose, permuted_graph, test_dof_mapping);

          // render composite test-dof-mapping/trial-dof-support graph
          Adjacency::Graph dof_adjactor(Adjacency::RenderType::injectify_sorted, test_dof_support, trial_dof_mapping);

          // return the graph
          return dof_adjactor;
        }
        else // no permutation
        {
          // render transposed test-dof-mapping
          Adjacency::Graph test_dof_support(Adjacency::RenderType::injectify_transpose, trial2test_adjactor, test_dof_mapping);

          // render composite test-dof-mapping/trial-dof-support graph
          Adjacency::Graph dof_adjactor(Adjacency::RenderType::injectify_sorted, test_dof_support, trial_dof_mapping);

          // return the graph
          return dof_adjactor;
        }
      }

      /**
       * \brief Assembles the 2-level refinement Dof-Adjacency graph.
       *
       * \param[in] fine_space
       * The \transient test-space defined on the refined mesh.
       *
       * \param[in] coarse_space
       * The \transient trial-space defined on the coarse mesh.
       *
       * \attention
       * This function silently assumes that \p fine_space is defined on the
       * mesh that was obtained by applying the standard 2-level refinement
       * algorithm on the underlying mesh of \p coarse_space.
       *
       * \returns
       * The standard Dof-Adjacency graph for the fine-test- and coarse-trial-space combination.
       */
      template<typename FineTestSpace_, typename CoarseTrialSpace_>
      static Adjacency::Graph assemble_graph_2lvl(
        const FineTestSpace_& fine_space,
        const CoarseTrialSpace_& coarse_space)
      {
        // create an refinement adjactor
        Geometry::Intern::CoarseFineCellMapping<typename FineTestSpace_::MeshType, typename CoarseTrialSpace_::MeshType>
            refine_adjactor(fine_space.get_trafo().get_mesh(), coarse_space.get_trafo().get_mesh());

        // call inter-mesh assembler
        return assemble_graph_intermesh(fine_space, coarse_space, refine_adjactor);
      }

      /**
       * \brief Assembles a standard matrix structure from a test-trial-space pair.
       *
       * \param[out] matrix
       * A \transient reference to the matrix to be assembled.
       *
       * \param[in] test_space, trial_space
       * The \transient test- and trial-spaces to be used for the assembly.
       */
      template<typename MatrixType_, typename TestSpace_, typename TrialSpace_>
      static void assemble_matrix_std2(MatrixType_ & matrix,
        const TestSpace_& test_space, const TrialSpace_& trial_space)
      {
        matrix = MatrixType_(assemble_graph_std2(test_space, trial_space));
      }

      /// specialization for NullMatrix
      template<typename DT_, typename IT_, int BH_, int BW_, typename TestSpace_, typename TrialSpace_>
      static void assemble_matrix_std2(
        LAFEM::NullMatrix<DT_, IT_, BH_, BW_>& matrix,
        const TestSpace_& test_space, const TrialSpace_& trial_space)
      {
        // only the dimensions are required here
        matrix.resize(test_space.get_num_dofs(), trial_space.get_num_dofs());
      }

      /**
       * \brief Assembles a standard matrix structure from a single space.
       *
       * \param[out] matrix
       * A \transient reference to the matrix to be assembled.
       *
       * \param[in] space
       * The \transient space to be used for the assembly.
       */
      template<typename MatrixType_, typename Space_>
      static void assemble_matrix_std1(MatrixType_ & matrix, const Space_& space)
      {
        matrix = MatrixType_(assemble_graph_std1(space));
      }

      /// specialization for NullMatrix
      template<typename DT_, typename IT_, int BH_, int BW_, typename Space_>
      static void assemble_matrix_std1(
        LAFEM::NullMatrix<DT_, IT_, BH_, BW_>& matrix, const Space_& space)
      {
        // only the dimensions are required here
        matrix.resize(space.get_num_dofs(), space.get_num_dofs());
      }

      /**
       * \brief Assembles an extended-facet matrix structure from a test-trial-space pair.
       *
       * \param[out] matrix
       * A \transient reference to the matrix to be assembled.
       *
       * \param[in] test_space, trial_space
       * The \transient test- and trial-spaces to be used for the assembly.
       */
      template<typename MatrixType_, typename TestSpace_, typename TrialSpace_>
      static void assemble_matrix_ext_facet2(MatrixType_ & matrix,
        const TestSpace_& test_space, const TrialSpace_& trial_space)
      {
        matrix = MatrixType_(assemble_graph_ext_facet2(test_space, trial_space));
      }

      /**
       * \brief Assembles an extended-facet matrix structure from a single space.
       *
       * \param[out] matrix
       * A \transient reference to the matrix to be assembled.
       *
       * \param[in] space
       * The \transient space to be used for the assembly.
       */
      template<typename MatrixType_, typename Space_>
      static void assemble_matrix_ext_facet1(MatrixType_ & matrix, const Space_& space)
      {
        matrix = MatrixType_(assemble_graph_ext_facet1(space));
      }

      /**
       * \brief Assembles an extended-node matrix structure from a test-trial-space pair.
       *
       * \param[out] matrix
       * A \transient reference to the matrix to be assembled.
       *
       * \param[in] test_space, trial_space
       * The \transient test- and trial-spaces to be used for the assembly.
       */
      template<typename MatrixType_, typename TestSpace_, typename TrialSpace_>
      static void assemble_matrix_ext_node2(MatrixType_ & matrix,
        const TestSpace_& test_space, const TrialSpace_& trial_space)
      {
        matrix = MatrixType_(assemble_graph_ext_node2(test_space, trial_space));
      }

      /**
       * \brief Assembles an extended-node matrix structure from a single space.
       *
       * \param[out] matrix
       * A \transient reference to the matrix to be assembled.
       *
       * \param[in] space
       * The \transient space to be used for the assembly.
       */
      template<typename MatrixType_, typename Space_>
      static void assemble_matrix_ext_node1(MatrixType_ & matrix, const Space_& space)
      {
        matrix = MatrixType_(assemble_graph_ext_node1(space));
      }

      /**
       * \brief Assembles a diagonal matrix structure from a single space.
       *
       * \param[out] matrix
       * A \transient reference to the matrix to be assembled.
       *
       * \param[in] space
       * The \transient space to be used for the assembly.
       */
      template<typename MatrixType_, typename Space_>
      static void assemble_matrix_diag(MatrixType_ & matrix, const Space_& space)
      {
        matrix = MatrixType_(assemble_graph_diag(space));
      }

      /**
       * \brief Assembles a matrix structure from a test-trial-mesh pair defined on different meshes.
       *
       * \param[out] matrix
       * A \transient reference to the matrix to be assembled.
       *
       * \param[in] test_space, trials_space
       * The \transient fine and coarse spaces to be used for the assembly.
       *
       * \param[in] trial2test_adjactor
       * A \transient object implementing the adjactor interface, which maps from the set of trial space mesh cells
       * to the test space mesh cells, which iterates over the set of all test space mesh cells intersecting with
       * a given trial space mesh cell.
       */
      template<typename MatrixType_, typename TestSpace_, typename TrialSpace_, typename Trial2TestAdjactor_>
      static void assemble_matrix_intermesh(MatrixType_ & matrix,
        const TestSpace_& test_space, const TrialSpace_& trial_space,
        const Trial2TestAdjactor_& trial2test_adjactor)
      {
        matrix = MatrixType_(assemble_graph_intermesh(test_space, trial_space, trial2test_adjactor));
      }

      /**
       * \brief Assembles a 2-level matrix structure from a fine-coarse-space pair.
       *
       * \param[out] matrix
       * A \transient reference to the matrix to be assembled.
       *
       * \param[in] fine_space, coarse_space
       * The \transient fine and coarse spaces to be used for the assembly.
       *
       * \attention
       * This function silently assumes that \p fine_space is defined on the
       * mesh that was obtained by applying the standard 2-level refinement
       * algorithm on the underlying mesh of \p coarse_space.
       */
      template<typename MatrixType_, typename FineSpace_, typename CoarseSpace_>
      static void assemble_matrix_2lvl(MatrixType_ & matrix,
        const FineSpace_& fine_space, const CoarseSpace_& coarse_space)
      {
        matrix = MatrixType_(assemble_graph_2lvl(fine_space, coarse_space));
      }
    }; // class SymbolicAssembler
  } // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_SYMBOLIC_ASSEMBLER_HPP
