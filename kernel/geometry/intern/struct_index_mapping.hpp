#pragma once
#ifndef KERNEL_GEOMETRY_INTERN_STRUCT_INDEX_MAPPING_HPP
#define KERNEL_GEOMETRY_INTERN_STRUCT_INDEX_MAPPING_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      template<
        int shape_dim_,
        int cell_dim_,
        int face_dim_>
      struct StructIndexMapping;

      /// vert-at-edge for edge mesh
      template<>
      struct StructIndexMapping<1, 1, 0>
      {
        static Index compute(Index i, Index j, const Index num_slices[])
        {
          return i + j;
        }
      };

      /// vert-at-quad for quad mesh
      template<>
      struct StructIndexMapping<2, 2, 0>
      {
        static Index compute(Index i, Index j, const Index num_slices[])
        {
          // compute index
          return ((i / num_slices[0]) + ((j >> 1) & 0x1)) * (num_slices[0] + 1) + (i % num_slices[0]) + (j & 0x1);
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_INTERN_STRUCT_INDEX_MAPPING_HPP
