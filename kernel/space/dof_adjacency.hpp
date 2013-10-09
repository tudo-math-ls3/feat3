#pragma once
#ifndef KERNEL_SPACE_DOF_ADJACENCY_HPP
#define KERNEL_SPACE_DOF_ADJACENCY_HPP 1

// includes, FEAST
#include <kernel/adjacency/graph.hpp>

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
     * \brief Dof-Mapping renderer class
     *
     * \author Peter Zajac
     */
    class DofMappingRenderer
    {
    public:
      /**
       * \brief Renders the dof-mapping of a space into an adjacency graph.
       *
       * \param[in] space
       * A reference to a finite element space whose dof-mapping is to be rendered.
       *
       * \returns
       * The dof-mapping of the space given as an adjacency graph.
       */
      template<typename Space_>
      static Adjacency::Graph render(const Space_& space)
      {
        // create a dof-mapping
        typename Space_::DofMappingType dof_map(space);

        // fetch the number of cells and dofs
        const Index num_cells(space.get_mesh().get_num_entities(Space_::shape_dim));
        const Index num_dofs(space.get_num_dofs());

        // allocate the domain-pointer
        Index* dom_ptr = new Index[num_cells+1];

        // loop over all cells and build the domain pointer array
        dom_ptr[0] = Index(0);
        for(Index i(0); i < num_cells; ++i)
        {
          Index num_adj(0);
          dof_map.prepare(i);
          for(Index j(0); j < dof_map.get_num_local_dofs(); ++j)
            num_adj += dof_map.get_num_contribs(j);
          dof_map.finish();
          dom_ptr[i+1] = dom_ptr[i] + num_adj;
        }

        // allocate the index array
        Index* img_idx = new Index[dom_ptr[num_cells]];

        // loop over all cells and build the image index array
        for(Index i(0); i < num_cells; ++i)
        {
          Index l(dom_ptr[i]);
          dof_map.prepare(i);
          for(Index j(0); j < dof_map.get_num_local_dofs(); ++j)
          {
            for(Index k(0); k < dof_map.get_num_contribs(j); ++k, ++l)
              img_idx[l] = dof_map.get_index(j, k);
          }
          dof_map.finish();
        }

        // create an adjacency graph
        return Adjacency::Graph(num_cells, num_dofs, dom_ptr[num_cells], dom_ptr, nullptr, img_idx, false);
      }
    }; // class DofMappingRenderer

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
      static Adjacency::Graph assemble(const TestSpace_& test_space, const TrialSpace_& trial_space)
      {
        // render dof-graphs
        Adjacency::Graph test_dof_graph(DofMappingRenderer::render(test_space));
        Adjacency::Graph trial_dof_graph(DofMappingRenderer::render(trial_space));

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
      static Adjacency::Graph assemble(const Space_& space)
      {
        // create dof-mapping
        Adjacency::Graph dof_graph(DofMappingRenderer::render(space));

        // render transposed dof-mapping
        Adjacency::Graph dof_support(Adjacency::rt_transpose, dof_graph);

        // render composite dof-mapping/dof-support graph
        Adjacency::Graph dof_adjactor(Adjacency::rt_injectify, dof_support, dof_graph);

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
      static Adjacency::Graph assemble(
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
        Adjacency::Graph test_dof_mapping(DofMappingRenderer::render(fine_test_space));
        Adjacency::Graph trial_dof_mapping(DofMappingRenderer::render(coarse_trial_space));

        // render transposed test-dof-mapping
        Adjacency::Graph test_dof_support(Adjacency::rt_injectify_transpose, refine_adjactor, test_dof_mapping);

        // render composite test-dof-mapping/trial-dof-support graph
        Adjacency::Graph dof_adjactor(Adjacency::rt_injectify, test_dof_support, trial_dof_mapping);

        // sort the dof-adjactor graph
        dof_adjactor.sort_indices();

        // return the graph
        return dof_adjactor;
      }
    }; // class DofAdjacency<Stencil::StandardRefinement>
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_DOF_ADJACENCY_HPP
