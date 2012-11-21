#pragma once
#ifndef KERNEL_GEOMETRY_INTERN_STRUCT_NUM_ENTITIES_HPP
#define KERNEL_GEOMETRY_INTERN_STRUCT_NUM_ENTITIES_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      /**
       * \brief Num-Entities computer for structured meshes
       *
       * This class computes the number of entities of dimension \p face_dim_ in
       * a structured mesh of dimension \p shape_dim_.
       *
       * \author Peter Zajac
       */
      template<
        int shape_dim_,
        int face_dim_>
      struct StructNumEntities DOXY({});

      template<int shape_dim_>
      struct StructNumEntities<shape_dim_, shape_dim_>
      {
        static Index compute(const Index num_slices[])
        {
          // compute number of elements
          Index count = num_slices[0];
          for(int i(1); i < shape_dim_; ++i)
            count *= (num_slices[i]);
          return count;
        }
      };

      template<int shape_dim_>
      struct StructNumEntities<shape_dim_, 0>
      {
        static Index compute(const Index num_slices[])
        {
          // compute number of vertices
          Index count = num_slices[0] + 1;
          for(int i(1); i < shape_dim_; ++i)
            count *= (num_slices[i] + 1);
          return count;
        }
      };

      template<>
      struct StructNumEntities<2, 1>
      {
        static Index compute(const Index num_slices[])
        {
          // compute number of edges in quad mesh
          return
            num_slices[0] * (num_slices[1] + 1) +
            num_slices[1] * (num_slices[0] + 1);
        }
      };

      template<>
      struct StructNumEntities<3, 1>
      {
        static Index compute(const Index num_slices[])
        {
          // compute number of edges in hexa mesh
          return
            num_slices[0] * (num_slices[1] + 1) * (num_slices[2] + 1) +
            num_slices[1] * (num_slices[0] + 1) * (num_slices[2] + 1) +
            num_slices[2] * (num_slices[0] + 1) * (num_slices[1] + 1);
        }
      };

      template<>
      struct StructNumEntities<3, 2>
      {
        static Index compute(const Index num_slices[])
        {
          // compute number of quads in hexa mesh
          return
            (num_slices[0] + 1) * num_slices[1] * num_slices[2] +
            (num_slices[1] + 1) * num_slices[0] * num_slices[2] +
            (num_slices[2] + 1) * num_slices[0] * num_slices[1];
        }
      };

      /**
       * \brief Computes the num_entities array
       *
       * This class template computes the num_entities array from a structured mesh's num_slices
       * array.
       *
       * \author Peter Zajac
       */
      template<
        int shape_dim_,
        int cell_dim_ = shape_dim_>
      struct StructNumEntitiesWrapper
      {
        static void compute(Index num_entities[], const Index num_slices[])
        {
          StructNumEntitiesWrapper<shape_dim_, cell_dim_ - 1>::compute(num_entities, num_slices);
          num_entities[cell_dim_] = StructNumEntities<shape_dim_, cell_dim_>::compute(num_slices);
        }
      };

      template<int shape_dim_>
      struct StructNumEntitiesWrapper<shape_dim_, 0>
      {
        static void compute(Index num_entities[], const Index num_slices[])
        {
          num_entities[0] = StructNumEntities<shape_dim_, 0>::compute(num_slices);
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_INTERN_STRUCT_NUM_ENTITIES_HPP
