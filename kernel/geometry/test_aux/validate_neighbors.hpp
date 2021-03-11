// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_TEST_AUX_VALIDATE_NEIGHBOURS_HPP
#define KERNEL_GEOMETRY_TEST_AUX_VALIDATE_NEIGHBOURS_HPP 1
#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/mesh_part.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace TestAux
    {
      /**
       * \brief Checks if the neighbor information is consistent.
       *
       * \tparam MeshType
       * Type of the mesh to validate
       *
       * \param[in] mesh
       * The mesh to check
       *
       * This checks:
       *  - if k1 has a neighbor k2, then k2 has k1 as neighbor as well
       *  - if k1's facet l says there is no neighbor, then l is at the boundary
       *
       */
      template<typename MeshType>
      void validate_neighbors(const MeshType& mesh)
      {
        static constexpr int facet_dim = MeshType::shape_dim-1;
        const auto& neigh = mesh.get_neighbors();
        const auto& facet_idx = mesh.template get_index_set<MeshType::shape_dim, facet_dim>();

        Geometry::BoundaryFactory<MeshType> boundary_factory(mesh);
        Geometry::MeshPart<MeshType> boundary(boundary_factory);

        Index num_facets(mesh.get_num_entities(MeshType::shape_dim-1));

        bool* at_boundary(new bool[num_facets]);

        for(Index l(0); l < num_facets; ++l)
          at_boundary[l] = false;

        // Mark all boundary facets
        const auto& facet_trg = boundary.template get_target_set<MeshType::shape_dim-1>();
        for(Index l(0); l < facet_trg.get_num_entities(); ++l)
        {
          at_boundary[facet_trg[l]] = true;
        }

        // Check the mesh
        for(Index k(0); k < mesh.get_num_entities(MeshType::shape_dim); ++k)
        {
          for(int j(0); j < facet_idx.num_indices; ++j)
          {
            Index other_cell(neigh(k,j));
            // If we have a neighbor...
            if(other_cell != ~Index(0))
            {
              bool ok(false);
              // ... then this should NOT be the boundary
              if(at_boundary[facet_idx(k,j)])
                throw String("Facet "+stringify(facet_idx(k,j))+
                    " is at the boundary, but the _neighborinformation says it is not!");

              // Check for vice versa neighbor relationship
              for(int i(0); i < facet_idx.num_indices; ++i)
              {
                if(neigh(other_cell,i) == k)
                {
                  ok = true;
                  break;
                }
              }
              if(!ok)
                throw String("Cell "+stringify(k)+" has neighbor "+stringify(other_cell)+" but not vice versa!");
            }
            else
            {
              /// So we do not have a neighbor. Then we SHOULD be at the boundary
              if(!at_boundary[facet_idx(k,j)])
               throw String("Facet "+stringify(facet_idx(k,j))+
                   " is not at the boundary, but the _neighbor information says it is!");
            }
          }

        }

        // Clean up
        delete[] at_boundary;

      }
    } // namespace TestAux
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_TEST_AUX_VALIDATE_NEIGHBOURS_HPP
