// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

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
        const Index num_cells(space.get_mesh().get_num_entities(Space_::shape_dim)); //num_cells
        const Index num_dofs(space.get_num_dofs()); //num_dofs

        // loop over all cells and build the domain pointer array

        Index num_indices_image(0);
        for(Index i(0); i < num_cells; ++i)
        {
          dof_map.prepare(i);
          num_indices_image += Index(dof_map.get_num_local_dofs());
          dof_map.finish();
        }
        // create an adjacency graph
        Adjacency::Graph graph (num_cells, num_dofs, num_indices_image) ;
        Index* dom_ptr = graph.get_domain_ptr();
        Index* img_idx = graph.get_image_idx();

        //build the domain pointer array
        //loop over all cells and build the image index array
        dom_ptr[0] = Index(0);

        for(Index i(0); i < num_cells; ++i)
        {
          Index l(dom_ptr[i]);
          dof_map.prepare(i);
          dom_ptr[i + 1] = l + Index(dof_map.get_num_local_dofs());
          for(int j(0); j < dof_map.get_num_local_dofs(); ++j, ++l)
          {
            img_idx[l] = dof_map.get_index(j);
          }
          dof_map.finish();
        }

        // create an adjacency graph
        return graph;
      }
    }; // class DofMappingRenderer
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_DOF_MAPPING_RENDERER_HPP
