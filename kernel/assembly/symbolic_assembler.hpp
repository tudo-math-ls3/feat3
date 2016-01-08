#pragma once
#ifndef KERNEL_ASSEMBLY_SYMBOLIC_ASSEMBLER_HPP
#define KERNEL_ASSEMBLY_SYMBOLIC_ASSEMBLER_HPP 1

// includes, FEAST
#include <kernel/space/dof_mapping_renderer.hpp>

// includes, FEAST-LAFEM
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>

namespace FEAST
{
  namespace Assembly
  {
    /**
     * \brief Assembly Stencil namespace
     *
     * This namespace contains several stencil tag classes which are used by the DofAdjacency
     * class template to assemble dof adjacency graphs.
     */
    namespace Stencil
    {
      /**
       * \brief Standard adjacency stencil tag class
       *
       * This is the stencil tag that is used to assemble standard operators.
       */
      class Standard
      {
      };

      /**
       * \brief Extended-Facet adjacency stencil tag class
       *
       * This is the stencil tag that is used to assemble operators which have additional dof couplings
       * over the facets.\n
       * This stencil is necessary for dicontinous Galerkin and edge-oriented jump-stabilisation operators.
       */
      class ExtendedFacet
      {
      };

      /**
       * \brief Standard Refinement adjacency stencil tag class
       *
       * This is the stencil tag that is used to assemble operators with test- and trial-spaces on two
       * different meshes, where the trial-space is defined on the mesh refined from the test-space mesh
       * using the standard refinement algorithm.\n
       * This stencil is necessary for prolongation and restriction operators.
       */
      class StandardRefinement
      {
      };
    } // namespace Stencil

    template<typename Stencil_ = Stencil::Standard>
    class SymbolicGraphAssembler DOXY({});

    template<>
    class SymbolicGraphAssembler<Stencil::Standard>
    {
    public:
      /**
       * \brief Assembles the standard Dof-Adjacency graph for different test- and trial-spaces.
       *
       * \param[in] test_space, trial_space
       * The test- and trial-spaces to be used for the assembly. Must be defined on the same mesh.
       *
       * \returns
       * The standard Dof-Adjacency graph of the test- and trial-space combination.
       */
      template<
        typename TestSpace_,
        typename TrialSpace_>
      static Adjacency::Graph assemble_graph(const TestSpace_& test_space, const TrialSpace_& trial_space)
      {
        // render dof-graphs
        Adjacency::Graph test_dof_graph(Space::DofMappingRenderer::render(test_space));
        Adjacency::Graph trial_dof_graph(Space::DofMappingRenderer::render(trial_space));

        // render transposed test-dof-mapping
        Adjacency::Graph test_dof_support(Adjacency::rt_transpose, test_dof_graph);

        // render composite test-dof-mapping/trial-dof-support graph
        Adjacency::Graph dof_adjactor(Adjacency::rt_injectify, test_dof_support, trial_dof_graph);

        // sort the dof-adjactor graph
        dof_adjactor.sort_indices();

        // return the graph
        return dof_adjactor;
      }

      /**
       * \brief Assembles the standard Dof-Adjacency graph for identical test- and trial-spaces.
       *
       * \param[in] space
       * The space representing the test- and trial spaces to be used for the assembly.
       *
       * \returns
       * The standard Dof-Adjacency graph of the space.
       */
      template<typename Space_>
      static Adjacency::Graph assemble_graph(const Space_& space)
      {
        // create dof-mapping
        Adjacency::Graph dof_graph(Space::DofMappingRenderer::render(space));

        // render transposed dof-mapping
        Adjacency::Graph dof_support(Adjacency::rt_transpose, dof_graph);

        // render composite dof-mapping/dof-support graph
        Adjacency::Graph dof_adjactor(Adjacency::rt_injectify, dof_support, dof_graph);

        // sort the sof-adjactor graph
        dof_adjactor.sort_indices();

        // return the graph
        return dof_adjactor;
      }
    }; // class SymbolicGraphAssembler<Stencil::Standard>

    template<>
    class SymbolicGraphAssembler<Stencil::ExtendedFacet>
    {
    public:
      /**
       * \brief Assembles the extended-facet Dof-Adjacency graph for different test- and trial-spaces.
       *
       * \param[in] test_space, trial_space
       * The test- and trial-spaces to be used for the assembly. Must be defined on the same mesh.
       *
       * \returns
       * The extended-facet Dof-Adjacency graph of the test- and trial-space combination.
       */
      template<
        typename TestSpace_,
        typename TrialSpace_>
      static Adjacency::Graph assemble_graph(const TestSpace_& test_space, const TrialSpace_& trial_space)
      {
        // render dof-graphs
        Adjacency::Graph test_dof_graph(Space::DofMappingRenderer::render(test_space));
        Adjacency::Graph trial_dof_graph(Space::DofMappingRenderer::render(trial_space));

        // get the shape dimension
        typedef typename TestSpace_::ShapeType ShapeType;
        static constexpr Index shape_dim = Index(ShapeType::dimension);

        // get the facet index set
        const auto& facet_at_shape = test_space.get_trafo().get_mesh().template get_index_set<shape_dim, shape_dim-1>();

        // transpose to get the facet support
        Adjacency::Graph shape_at_facet(Adjacency::rt_transpose, facet_at_shape);

        // render transposed test-dof-mapping
        Adjacency::Graph test_dof_support(Adjacency::rt_transpose, test_dof_graph);

        // render extended test-dof-support
        Adjacency::Graph test_dof_ext_sup(Adjacency::rt_injectify, test_dof_support, facet_at_shape);

        // render extended trial-dof-mapping
        Adjacency::Graph trial_dof_ext_graph(Adjacency::rt_injectify, shape_at_facet, trial_dof_graph);

        // render composite test-dof-mapping/trial-dof-support graph
        Adjacency::Graph dof_adjactor(Adjacency::rt_injectify, test_dof_ext_sup, trial_dof_ext_graph);

        // sort the dof-adjactor graph
        dof_adjactor.sort_indices();

        // return the graph
        return dof_adjactor;
      }

      /**
       * \brief Assembles the standard Dof-Adjacency graph for identical test- and trial-spaces.
       *
       * \param[in] space
       * The space representing the test- and trial spaces to be used for the assembly.
       *
       * \returns
       * The standard Dof-Adjacency graph of the space.
       */
      template<typename Space_>
      static Adjacency::Graph assemble_graph(const Space_& space)
      {
        // create dof-mapping
        Adjacency::Graph dof_graph(Space::DofMappingRenderer::render(space));

        // get the shape dimension
        typedef typename Space_::ShapeType ShapeType;
        static constexpr Index shape_dim = Index(ShapeType::dimension);

        // get the facet index set
        const auto& facet_at_shape = space.get_trafo().get_mesh().template get_index_set<shape_dim, shape_dim-1>();

        // transpose to get the facet support
        Adjacency::Graph shape_at_facet(Adjacency::rt_transpose, facet_at_shape);

        // render extended dof-mapping
        Adjacency::Graph dof_ext_graph(Adjacency::rt_injectify, shape_at_facet, dof_graph);

        // render transposed extended dof-mapping
        Adjacency::Graph dof_ext_sup(Adjacency::rt_transpose, dof_ext_graph);

        // render composite dof-mapping/dof-support graph
        Adjacency::Graph dof_adjactor(Adjacency::rt_injectify, dof_ext_sup, dof_ext_graph);

        // sort the sof-adjactor graph
        dof_adjactor.sort_indices();

        // return the graph
        return dof_adjactor;
      }
    }; // class SymbolicGraphAssembler<Stencil::ExtendedFacet>

    /**
     * \brief Standard-Refinement Dof-Adjacency assembler class
     *
     * This class assembles a standard Dof-Adjacency graph for a combination of test- and trial-spaces
     * where the test-space mesh is refined from the trial-space mesh.
     *
     * \author Peter Zajac
     */
    template<>
    class SymbolicGraphAssembler<Stencil::StandardRefinement>
    {
    private:
      /// \cond internal
      class RefinementAdjactor
      {
      public:
        typedef Adjacency::Adjactor::IndexImageIterator ImageIterator;

      private:
        Index _num_elements;
        Index _num_children;

      public:
        explicit RefinementAdjactor(Index num_elements, Index num_children) :
          _num_elements(num_elements),
          _num_children(num_children)
        {
        }

        Index get_num_nodes_domain() const
        {
          return _num_elements;
        }

        Index get_num_nodes_image() const
        {
          return _num_elements * _num_children;
        }

        ImageIterator image_begin(Index domain_node) const
        {
          return ImageIterator(_num_children * domain_node);
        }

        ImageIterator image_end(Index domain_node) const
        {
          return ImageIterator(_num_children * (domain_node + 1));
        }
      };
      /// \endcond

    public:
      /**
       * \brief Assembles the standard Dof-Adjacency graph.
       *
       * \param[in] fine_test_space
       * The test-space defined on the refined mesh.
       *
       * \param[in] coarse_trial_space
       * The trial-space defined on the coarse mesh.
       *
       * \returns
       * The standard Dof-Adjacency graph for the fine-test- and coarse-trial-space combination.
       */
      template<
        typename FineTestSpace_,
        typename CoarseTrialSpace_>
      static Adjacency::Graph assemble_graph(
        const FineTestSpace_& fine_test_space,
        const CoarseTrialSpace_& coarse_trial_space)
      {
        // fetch the shape of the space
        typedef typename CoarseTrialSpace_::ShapeType ShapeType;

        // fetch number of coarse mesh elements
        Index num_elements_coarse = coarse_trial_space.get_trafo().get_mesh().get_num_entities(ShapeType::dimension);

        // fetch number of fine mesh elements
        Index num_elements_fine = fine_test_space.get_trafo().get_mesh().get_num_entities(ShapeType::dimension);

        // calculate child count
        Index num_children = num_elements_fine / num_elements_coarse;

        // create an refinement adjactor
        RefinementAdjactor refine_adjactor(num_elements_coarse, num_children);

        // create test- and trial-dof-mappers
        Adjacency::Graph test_dof_mapping(Space::DofMappingRenderer::render(fine_test_space));
        Adjacency::Graph trial_dof_mapping(Space::DofMappingRenderer::render(coarse_trial_space));

        // render transposed test-dof-mapping
        Adjacency::Graph test_dof_support(Adjacency::rt_injectify_transpose, refine_adjactor, test_dof_mapping);

        // render composite test-dof-mapping/trial-dof-support graph
        Adjacency::Graph dof_adjactor(Adjacency::rt_injectify, test_dof_support, trial_dof_mapping);

        // sort the dof-adjactor graph
        dof_adjactor.sort_indices();

        // return the graph
        return dof_adjactor;
      }
    }; // class SymbolicGraphAssembler<Stencil::StandardRefinement>

    /**
     * \brief Symbolic matrix assembler base-class
     *
     * \author Peter Zajac
     */
    class SymbolicMatrixAssemblerBase
    {
    public:
      /**
       * \brief Assembles a matrix from a Graph.
       *
       * \tparam MatrixType_
       * The type of the lafem matrix to be assembled.
       *
       * \param[out] matrix
       * A reference to the matrix to be assembled.
       *
       * \param[in] graph
       * The graph representing the sparsity pattern.
       */
      template<typename MatrixType_>
      static void assemble(MatrixType_& matrix, const Adjacency::Graph& graph)
      {
        // build the matrix
        matrix = MatrixType_(graph);
      }
    };

    /**
     * \brief Symbolic Matrix assembler class template
     *
     * \tparam Stencil_
     * One of the tag classes defined in the Stencil namespace identifying the type of the matrix
     * stencil to be assembled.
     *
     * \author Peter Zajac
     */
    template<typename Stencil_ = Stencil::Standard>
    class SymbolicMatrixAssembler :
      public SymbolicMatrixAssemblerBase
    {
    public:
      /**
       * \brief Assembles a matrix from a test-trial-space pair.
       *
       * \param[out] matrix
       * A reference to the matrix to be assembled.
       *
       * \param[in] test_space, trial_space
       * The test- and trial-spaces to be used for the assembly.
       */
      template<
        typename MatrixType_,
        typename TestSpace_,
        typename TrialSpace_>
      static void assemble2(MatrixType_ & matrix,
        const TestSpace_& test_space, const TrialSpace_& trial_space)
      {
        // assemble Graph
        SymbolicMatrixAssemblerBase::assemble(matrix,
          SymbolicGraphAssembler<Stencil_>::assemble_graph(test_space, trial_space));
      }

      /**
       * \brief Assembles a matrix from a single space.
       *
       * \param[out] matrix
       * A reference to the matrix to be assembled.
       *
       * \param[in] space
       * The space to be used for the assembly.
       */
      template<
        typename MatrixType_,
        typename Space_>
      static void assemble1(MatrixType_ & matrix, const Space_& space)
      {
        // assemble Graph
        SymbolicMatrixAssemblerBase::assemble(matrix,
          SymbolicGraphAssembler<Stencil_>::assemble_graph(space));
      }
    }; // class SymbolicMatrixAssembler<...>

    /**
     * \brief Symbolic Matrix assembler class template for StandardRefinement stencil
     *
     * \author Peter Zajac
     */
    template<>
    class SymbolicMatrixAssembler<Stencil::StandardRefinement> :
      public SymbolicMatrixAssemblerBase
    {
    public:
      /**
       * \brief Assembles a matrix from a fine-coarse-space pair.
       *
       * \param[out] matrix
       * A reference to the matrix to be assembled.
       *
       * \param[in] fine_space, coarse_space
       * The fine and coarse spaces to be used for the assembly.
       */
      template<
        typename MatrixType_,
        typename FineSpace_,
        typename CoarseSpace_>
      static void assemble(MatrixType_ & matrix,
        const FineSpace_& fine_space, const CoarseSpace_& coarse_space)
      {
        // assemble Graph
        SymbolicMatrixAssemblerBase::assemble(matrix,
          SymbolicGraphAssembler<Stencil::StandardRefinement>::assemble_graph(fine_space, coarse_space));
      }
    }; // class SymbolicMatrixAssembler<...>
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_SYMBOLIC_ASSEMBLER_HPP
