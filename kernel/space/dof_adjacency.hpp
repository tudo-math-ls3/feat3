#pragma once
#ifndef KERNEL_SPACE_DOF_ADJACENCY_HPP
#define KERNEL_SPACE_DOF_ADJACENCY_HPP 1

// includes, FEAST
#include <kernel/util/graph.hpp>

namespace FEAST
{
  namespace Space
  {
    /**
     * \brief Finite-Element Stencil namespace
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

    /**
     * \brief Dof-Adjacency assembler class template
     *
     * This class is used for the assembly of dof-adjacency graphs of operators.
     *
     * \tparam Stencil_
     * The stencil that is to be assembled. See Stencil namespace for details.
     *
     * \author Peter Zajac
     */
    template<typename Stencil_ = Stencil::Standard>
    class DofAdjacency DOXY({});

    /**
     * \brief Standard-Stencil Dof-Adjacency assembler class
     *
     * This class assembles a standard Dof-Adjacency graph for a combination of test- and trial-spaces
     * on the same mesh.
     *
     * \author Peter Zajac
     */
    template<>
    class DofAdjacency<Stencil::Standard>
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
      static Graph assemble(const TestSpace_& test_space, const TrialSpace_& trial_space)
      {
        // create test- and trial-dof-mappers
        typename TestSpace_::DofMappingType test_dof_mapping(test_space);
        typename TrialSpace_::DofMappingType trial_dof_mapping(trial_space);

        // render transposed test-dof-mapping
        Graph test_dof_support(Graph::rt_transpose, test_dof_mapping);

        // render composite test-dof-mapping/trial-dof-support graph
        Graph dof_adjactor(Graph::rt_injectify, test_dof_support, trial_dof_mapping);

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
      static Graph assemble(const Space_& space)
      {
        // create dof-mapping
        typename Space_::DofMappingType dof_mapping(space);

        // render transposed dof-mapping
        Graph dof_support(Graph::rt_transpose, dof_mapping);

        // render composite dof-mapping/dof-support graph
        Graph dof_adjactor(Graph::rt_injectify, dof_support, dof_mapping);

        // sort the sof-adjactor graph
        dof_adjactor.sort_indices();

        // return the graph
        return dof_adjactor;
      }
    }; // class DofAdjacency<Stencil::Standard>

    /**
     * \brief Standard-Refinement Dof-Adjacency assembler class
     *
     * This class assembles a standard Dof-Adjacency graph for a combination of test- and trial-spaces
     * where the test-space mesh is refined from the trial-space mesh.
     *
     * \author Peter Zajac
     */
    template<>
    class DofAdjacency<Stencil::StandardRefinement>
    {
    private:
      /// \cond internal
      class RefinementAdjactor
      {
      public:
        typedef Adjactor::IndexImageIterator ImageIterator;

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
      static Graph assemble(
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
        typename FineTestSpace_::DofMappingType test_dof_mapping(fine_test_space);
        typename CoarseTrialSpace_::DofMappingType trial_dof_mapping(coarse_trial_space);

        // render transposed test-dof-mapping
        Graph test_dof_support(Graph::rt_injectify_transpose, refine_adjactor, test_dof_mapping);

        // render composite test-dof-mapping/trial-dof-support graph
        Graph dof_adjactor(Graph::rt_injectify, test_dof_support, trial_dof_mapping);

        // sort the dof-adjactor graph
        dof_adjactor.sort_indices();

        // return the graph
        return dof_adjactor;
      }
    }; // class DofAdjacency<Stencil::StandardRefinement>
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_DOF_ADJACENCY_HPP
