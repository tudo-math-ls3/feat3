#pragma once
#ifndef FEAT_KERNEL_GEOMETRY_INTERN_FACET_NEIGHBOURS_HPP
#define FEAT_KERNEL_GEOMETRY_INTERN_FACET_NEIGHBOURS_HPP 1

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
       * \brief Wrapper class for computing facet neighbours
       */
      struct FacetNeighbours
      {
        // \brief Descriptive String
        static String name()
        {
          return "FacetNeighbours";
        }

        /**
         * \brief Computes neighbour information for entities sharing a facet
         *
         * \param[out] neighbours
         * The neighbours IndexSet. This will be of the same type as the IndexSet<shape_dim, shape_dim -1> in
         * ConformalMesh.
         *
         * \param[in] facet_idx
         * Facet at cell IndexSet.
         *
         */
        template<typename NeighbourIndexSetType_, typename FacetIndexSetType_>
        static void compute(NeighbourIndexSetType_& neighbours, const FacetIndexSetType_& facet_idx)
        {
          Index num_cells(facet_idx.get_num_entities());
          Index num_facets(facet_idx.get_index_bound());

          ASSERT(neighbours.get_num_entities() == num_cells, "neighbours.get_num_entitities() = "+stringify(neighbours.get_num_entities())+ " != " +stringify(num_cells) + " = num_cells");

          // A facet is shared by exactly 2 cells if it is interiour, and is present in exactly one cell if it is at
          // the boundary
          typedef Index SharedBy[2];
          SharedBy* shared_by(new SharedBy[num_facets]);

          // ~Index(0) is the marker for "no neighbour"
          for(Index l(0); l < num_facets; ++l)
          {
            shared_by[l][0] = ~Index(0);
            shared_by[l][1] = ~Index(0);
          }

          // For each facet, find the cells sharing it
          for(Index k(0); k < num_cells; ++k)
          {
            for(int j(0); j < facet_idx.num_indices; ++j)
            {
              // Index of the facet
              Index l(facet_idx[k][Index(j)]);

              if(shared_by[l][0] == ~Index(0))
                shared_by[l][0] = k;
              else if(shared_by[l][1] == ~Index(0))
                shared_by[l][1] = k;
              else
                throw InternalError("Facet "+stringify(l)+" is shared by cells "+stringify(shared_by[l][0])+", "+stringify(shared_by[l][1])+" and again by "+stringify(k));
            }
          }

          // For every cell and for every facet of that cell, the neighbour at a face is the OTHER cell sharing it
          // (if any)
          for(Index k(0); k < num_cells; ++k)
          {
            for(int j(0); j < facet_idx.num_indices; ++j)
            {
              // Index of the facet
              Index l(facet_idx[k][Index(j)]);

              if(shared_by[l][0] == k)
                neighbours[k][j] = shared_by[l][1];
              else if(shared_by[l][1] == k)
                neighbours[k][j] = shared_by[l][0];
              else
                throw InternalError("Facet "+stringify(l)+" found at cell "+stringify(k)+" but is shared by cells "+stringify(shared_by[l][0])+", "+stringify(shared_by[l][1]));
            }
          }

          delete[] shared_by;

        } // compute()

      }; // struct FacetNeighbours

    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // FEAT_KERNEL_GEOMETRY_INTERN_FACET_NEIGHBOURS_HPP
