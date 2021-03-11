// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_INTERN_DUAL_ADAPTOR_HPP
#define KERNEL_GEOMETRY_INTERN_DUAL_ADAPTOR_HPP 1

#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/vertex_set.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      template<typename Mesh_, typename Shape_ = typename Mesh_::ShapeType>
      struct DualAdaptor
      {
        static void adapt(Mesh_&, const Mesh_&)
        {
          // do nothing
        }
      };

      template<typename Mesh_>
      struct DualAdaptor<Mesh_, Shape::Hypercube<1> >
      {
        static void adapt(Mesh_&, const Mesh_&)
        {
          // do nothing
        }
      };

      template<typename Mesh_, int shape_dim_>
      struct DualAdaptor<Mesh_, Shape::Hypercube<shape_dim_> >
      {
        static_assert(shape_dim_ > 1, "invalid shape dimension");

        static void adapt(Mesh_& mesh_f, const Mesh_& mesh_c)
        {
          typedef typename Mesh_::CoordType CoordType;

          // get number of facets
          static constexpr int nfe = Shape::FaceTraits<Shape::Hypercube<shape_dim_>, shape_dim_-1>::count;

          // compute scaling factor
          const CoordType scale = CoordType(1) / CoordType(nfe);

          // get cell count
          Index num_cells = mesh_c.get_num_entities(shape_dim_);

          // compute facet-vertex and cell-vertex offsets
          Index fvo = Index(0);
          for(int i(0); (i+1) < shape_dim_; ++i)
          {
            fvo += mesh_c.get_num_entities(i);
          }
          Index cvo = fvo + mesh_c.get_num_entities(shape_dim_ - 1);

          // get fine mesh vertex set
          auto& vtx = mesh_f.get_vertex_set();

          // get coarse mesh facet index set
          const auto& facet = mesh_c.template get_index_set<shape_dim_, shape_dim_-1>();

          // loop over all coarse mesh quads
          for(Index i(0); i < num_cells; ++i)
          {
            // get the fine mesh vertex
            auto& v = vtx[cvo + i];
            v.format();

            // loop over all faces
            for(int j(0); j < nfe; ++j)
              v += scale * vtx[fvo + facet(i,j)];
          }
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_INTERN_DUAL_ADAPTOR_HPP
