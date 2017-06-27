#pragma once
#ifndef KERNEL_SPACE_DOF_MAPPING_RENDERER_HPP
#define KERNEL_SPACE_DOF_MAPPING_RENDERER_HPP 1

// includes, FEAT
#include <kernel/adjacency/graph.hpp>

namespace FEAT
{
  namespace Space
  {
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
          dof_map.prepare(i);
          dom_ptr[i+1] = dom_ptr[i] + Index(dof_map.get_num_local_dofs());
          dof_map.finish();
        }

        // allocate the index array
        Index* img_idx = new Index[dom_ptr[num_cells]];

        // loop over all cells and build the image index array
        for(Index i(0); i < num_cells; ++i)
        {
          Index l(dom_ptr[i]);
          dof_map.prepare(i);
          for(int j(0); j < dof_map.get_num_local_dofs(); ++j, ++l)
          {
            img_idx[l] = dof_map.get_index(j);
          }
          dof_map.finish();
        }

        // create an adjacency graph
        return Adjacency::Graph(num_cells, num_dofs, dom_ptr[num_cells], dom_ptr, img_idx, false);
      }
    }; // class DofMappingRenderer
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_DOF_MAPPING_RENDERER_HPP
