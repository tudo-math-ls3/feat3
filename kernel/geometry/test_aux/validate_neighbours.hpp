#pragma once
#ifndef KERNEL_GEOMETRY_TEST_AUX_VALIDATE_NEIGHBOURS_HPP
#define KERNEL_GEOMETRY_TEST_AUX_VALIDATE_NEIGHBOURS_HPP 1
#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/mesh_part.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /// \cond internal
    namespace TestAux
    {
      /**
       * \brief Checks if the neighbour information is consistent.
       *
       * \tparam MeshType
       * Type of the mesh to validate
       *
       * \param[in] mesh
       * The mesh to check
       *
       * This checks:
       *  - if k1 has a neighbour k2, then k2 has k1 as neighbour as well
       *  - if k1's facet l says there is no neighbour, then l is at the boundary
       *
       */
      template<typename MeshType>
      void validate_neighbours(const MeshType& mesh)
      {
        static constexpr int facet_dim = MeshType::shape_dim-1;
        const auto& neigh = mesh.get_neighbours();
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
            Index other_cell(neigh[k][Index(j)]);
            bool ok(false);
            // If we have a neighbour...
            if(other_cell != ~Index(0))
            {
              // ... then this should NOT be the boundary
              if(at_boundary[facet_idx[k][j]])
                throw String("Facet "+stringify(facet_idx[k][j])+
                    " is at the boundary, but the _neighbourinformation says it is not!");

              // Check for vice versa neighbour relationship
              for(int i(0); i < facet_idx.num_indices; ++i)
              {
                if(neigh[other_cell][Index(i)] == k)
                {
                  ok = true;
                  break;
                }
              }
              if(!ok)
                throw String("Cell "+stringify(k)+" has neighbour "+stringify(other_cell)+" but not vice versa!");
            }
            else
            {
              /// So we do not have a neighbour. Then we SHOULD be at the boundary
              if(!at_boundary[facet_idx[k][j]])
               throw String("Facet "+stringify(facet_idx[k][j])+
                   " is not at the boundary, but the _neighbour information says it is!");
            }
          }

        }

        // Clean up
        delete[] at_boundary;

      }
    } // namespace TestAux
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_TEST_AUX_VALIDATE_NEIGHBOURS_HPP
