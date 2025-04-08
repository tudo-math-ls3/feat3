// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

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
        // fetch the number of cells and dofs
        const Index num_cells(space.get_mesh().get_num_entities(Space_::shape_dim));
        const Index num_dofs(space.get_num_dofs());

        // allocate domain pointer vector
        std::vector<Index> dom_ptr(num_cells+1u);
        dom_ptr[0u] = Index(0);

#ifdef FEAT_HAVE_OMP
        FEAT_PRAGMA_OMP(parallel)
        {
          typename Space_::DofMappingType dof_map(space);

          FEAT_PRAGMA_OMP(for)
          for(Index i = 0; i < num_cells; ++i)
          {
            dof_map.prepare(i);
            dom_ptr[i+1u] = Index(dof_map.get_num_local_dofs());
            dof_map.finish();
          }
        }
        feat_omp_in_scan(dom_ptr.size(), dom_ptr.data(), dom_ptr.data());
#else // no FEAT_HAVE_OMP
        {
          typename Space_::DofMappingType dof_map(space);
          for(Index i = 0; i < num_cells; ++i)
          {
            dof_map.prepare(i);
            dom_ptr[i+1u] = dom_ptr[i] + Index(dof_map.get_num_local_dofs());
            dof_map.finish();
          }
        }
#endif // FEAT_HAVE_OMP

        // allocate image index vector
        std::vector<Index> img_idx(dom_ptr[num_cells]);

        FEAT_PRAGMA_OMP(parallel)
        {
          typename Space_::DofMappingType dof_map(space);

          FEAT_PRAGMA_OMP(for)
          for(Index i = 0; i < num_cells; ++i)
          {
            Index l(dom_ptr[i]);
            dof_map.prepare(i);
            for(int j(0); j < dof_map.get_num_local_dofs(); ++j, ++l)
            {
              img_idx[l] = dof_map.get_index(j);
            }
            dof_map.finish();
          }
        }

        // create an adjacency graph
        return Adjacency::Graph(num_dofs, std::move(dom_ptr), std::move(img_idx));
      }
    }; // class DofMappingRenderer
  } // namespace Space
} // namespace FEAT
