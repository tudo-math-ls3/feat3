// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/string.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      /**
       * \brief Wrapper class for computing facet neighbors
       */
      struct FacetNeighbors
      {
        // \brief Descriptive String
        static String name()
        {
          return "FacetNeighbors";
        }

        /**
         * \brief Computes neighbor information for entities sharing a facet
         *
         * \param[out] neighbors
         * The neighbors IndexSet. This will be of the same type as the IndexSet<shape_dim, shape_dim -1> in
         * ConformalMesh.
         *
         * \param[in] facet_idx
         * Facet at cell IndexSet.
         *
         */
        template<typename NeighborIndexSetType_, typename FacetIndexSetType_>
        static void compute(NeighborIndexSetType_& neighbors, const FacetIndexSetType_& facet_idx)
        {
          Index num_cells(facet_idx.get_num_entities());
          Index num_facets(facet_idx.get_index_bound());

          XASSERT(neighbors.get_num_entities() == num_cells);

          // A facet is shared by exactly 2 cells if it is interiour, and is present in exactly one cell if it is at
          // the boundary
          typedef Index SharedBy[2];
          std::vector<Index> shared_by_vec(size_t(2)*num_facets);
          auto shared_by = reinterpret_cast<SharedBy*>(shared_by_vec.data());

          // ~Index(0) is the marker for "no neighbor"
          FEAT_PRAGMA_OMP(parallel for)
          for(Index l(0); l < num_facets; ++l)
          {
            shared_by[l][0] = ~Index(0);
            shared_by[l][1] = ~Index(0);
          }

          // For each facet, find the cells sharing it
          FEAT_PRAGMA_OMP(parallel for)
          for(Index k(0); k < num_cells; ++k)
          {
            for(int j(0); j < facet_idx.num_indices; ++j)
            {
              // Index of the facet
              Index l(facet_idx[k][j]);

              if(shared_by[l][0] == ~Index(0))
                shared_by[l][0] = k;
              else if(shared_by[l][1] == ~Index(0))
                shared_by[l][1] = k;
              else
                XABORTM("Facet "+stringify(l)+" is shared by cells "+stringify(shared_by[l][0])+", "+stringify(shared_by[l][1])+" and again by "+stringify(k));
            }
          }

          // For every cell and for every facet of that cell, the neighbor at a face is the OTHER cell sharing it
          // (if any)
          FEAT_PRAGMA_OMP(parallel for)
          for(Index k(0); k < num_cells; ++k)
          {
            for(int j(0); j < facet_idx.num_indices; ++j)
            {
              // Index of the facet
              Index l(facet_idx[k][j]);

              if(shared_by[l][0] == k)
                neighbors[k][j] = shared_by[l][1];
              else if(shared_by[l][1] == k)
                neighbors[k][j] = shared_by[l][0];
              else
                XABORTM("Facet "+stringify(l)+" found at cell "+stringify(k)+" but is shared by cells "+stringify(shared_by[l][0])+", "+stringify(shared_by[l][1]));
            }
          }

        } // compute()

      }; // struct FacetNeighbors

    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT
