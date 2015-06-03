#pragma once
#ifndef KERNEL_FOUNDATION_HALO_INTERFACE_HPP
#define KERNEL_FOUNDATION_HALO_INTERFACE_HPP

#include <kernel/foundation/mesh.hpp>
#include <kernel/foundation/halo.hpp>
#include <kernel/foundation/halo_control.hpp>

namespace FEAST
{
  namespace Foundation
  {
    template<unsigned delta_, typename Dim_>
    struct HaloInterface
    {
    };

    template<unsigned delta_>
    struct HaloInterface<delta_, Dim1D>
    {
      ///for concrete halos
      template<typename Level_, typename MeshType_, typename WT_, template<typename, typename> class ST_>
      static Geometry::MeshPart<Geometry::ConformalMesh<Shape::Hypercube<1> > > convert(Halo<delta_, Level_, MeshType_, WT_, ST_> source)
      {
        typename MeshType_::index_type_ target_sizes[2];
        HaloControl<Dim1D::tag_value, WT_>::fill_sizes(source, &target_sizes[0]);
        Geometry::MeshPart<Geometry::ConformalMesh<Shape::Hypercube<1> > > result(&target_sizes[0]);
        HaloControl<Dim1D::tag_value, WT_>::fill_target_set(source, result);

        return result;
      }

      ///delegate-overload for anonymous halos
      template<typename MeshType_, typename WT_, template<typename, typename> class ST_>
      static Geometry::MeshPart<Geometry::ConformalMesh<Shape::Hypercube<1> > > convert(HaloBase<MeshType_, WT_, ST_>* source)
      {
        typename MeshType_::index_type_ target_sizes[2];

        switch(source->get_level())
        {
          case pl_vertex:
            HaloControl<Dim1D::tag_value, WT_>::fill_sizes(*static_cast<Halo<delta_, PLVertex, MeshType_, WT_, ST_>*>(source), &target_sizes[0]);
            break;
          case pl_edge:
            HaloControl<Dim1D::tag_value, WT_>::fill_sizes(*static_cast<Halo<delta_, PLEdge, MeshType_, WT_, ST_>*>(source), &target_sizes[0]);
            break;
          default:
            throw InternalError("Undefined polytope level!");
        }

        Geometry::MeshPart<Geometry::ConformalMesh<Shape::Hypercube<1> > > result(&target_sizes[0]);

        switch(source->get_level())
        {
          case pl_vertex:
            HaloControl<Dim1D::tag_value, WT_>::fill_target_set(*static_cast<Halo<delta_, PLVertex, MeshType_, WT_, ST_>*>(source), result);
            break;
          case pl_edge:
            HaloControl<Dim1D::tag_value, WT_>::fill_target_set(*static_cast<Halo<delta_, PLEdge, MeshType_, WT_, ST_>*>(source), result);
            break;
          default:
            throw InternalError("Undefined polytope level!");
        }

        return result;
      }
    };

    template<unsigned delta_>
    struct HaloInterface<delta_, Dim2D>
    {
      ///for concrete halos
      template<typename Level_, typename MeshType_, typename WT_, template<typename, typename> class ST_>
      static Geometry::MeshPart<Geometry::ConformalMesh<Shape::Hypercube<2> > > convert(Halo<delta_, Level_, MeshType_, WT_, ST_> source)
      {
        typename MeshType_::index_type_ target_sizes[3];
        HaloControl<Dim2D::tag_value, WT_>::fill_sizes(source, &target_sizes[0]);
        Geometry::MeshPart<Geometry::ConformalMesh<Shape::Hypercube<2> > > result(&target_sizes[0]);
        HaloControl<Dim2D::tag_value, WT_>::fill_target_set(source, result);

        return result;
      }

      ///delegate-overload for anonymous halos
      template<typename MeshType_, typename WT_, template<typename, typename> class ST_>
      static Geometry::MeshPart<Geometry::ConformalMesh<Shape::Hypercube<2> > > convert(HaloBase<MeshType_, WT_, ST_>* source)
      {
        typename MeshType_::index_type_ target_sizes[3];

        switch(source->get_level())
        {
          case pl_vertex:
            HaloControl<Dim2D::tag_value, WT_>::fill_sizes(*static_cast<Halo<delta_, PLVertex, MeshType_, WT_, ST_>*>(source), &target_sizes[0]);
            break;
          case pl_edge:
            HaloControl<Dim2D::tag_value, WT_>::fill_sizes(*static_cast<Halo<delta_, PLEdge, MeshType_, WT_, ST_>*>(source), &target_sizes[0]);
            break;
          case pl_face:
            HaloControl<Dim2D::tag_value, WT_>::fill_sizes(*static_cast<Halo<delta_, PLFace, MeshType_, WT_, ST_>*>(source), &target_sizes[0]);
            break;
          default:
            throw InternalError("Undefined polytope level!");
        }

        Geometry::MeshPart<Geometry::ConformalMesh<Shape::Hypercube<2> > > result(&target_sizes[0]);

        switch(source->get_level())
        {
          case pl_vertex:
            HaloControl<Dim2D::tag_value, WT_>::fill_target_set(*static_cast<Halo<delta_, PLVertex, MeshType_, WT_, ST_>*>(source), result);
            break;
          case pl_edge:
            HaloControl<Dim2D::tag_value, WT_>::fill_target_set(*static_cast<Halo<delta_, PLEdge, MeshType_, WT_, ST_>*>(source), result);
            break;
          case pl_face:
            HaloControl<Dim2D::tag_value, WT_>::fill_target_set(*static_cast<Halo<delta_, PLFace, MeshType_, WT_, ST_>*>(source), result);
            break;
          default:
            throw InternalError("Undefined polytope level!");
        }

        return result;
      }
    };

    template<unsigned delta_>
    struct HaloInterface<delta_, Dim3D>
    {
      ///for concrete halos
      template<typename Level_, typename MeshType_, typename WT_, template<typename, typename> class ST_>
      static Geometry::MeshPart<Geometry::ConformalMesh<Shape::Hypercube<3> > > convert(Halo<delta_, Level_, MeshType_, WT_, ST_> source)
      {
        typename MeshType_::index_type_ target_sizes[4];
        HaloControl<Dim3D::tag_value, WT_>::fill_sizes(source, &target_sizes[0]);
        Geometry::MeshPart<Geometry::ConformalMesh<Shape::Hypercube<3> > > result(&target_sizes[0]);
        HaloControl<Dim3D::tag_value, WT_>::fill_target_set(source, result);

        return result;
      }

      ///delegate-overload for anonymous halos
      template<typename MeshType_, typename WT_, template<typename, typename> class ST_>
      static Geometry::MeshPart<Geometry::ConformalMesh<Shape::Hypercube<3> > > convert(HaloBase<MeshType_, WT_, ST_>* source)
      {
        typename MeshType_::index_type_ target_sizes[4];

        switch(source->get_level())
        {
          case pl_vertex:
            HaloControl<Dim3D::tag_value, WT_>::fill_sizes(*static_cast<Halo<delta_, PLVertex, MeshType_, WT_, ST_>*>(source), &target_sizes[0]);
            break;
          case pl_edge:
            HaloControl<Dim3D::tag_value, WT_>::fill_sizes(*static_cast<Halo<delta_, PLEdge, MeshType_, WT_, ST_>*>(source), &target_sizes[0]);
            break;
          case pl_face:
            HaloControl<Dim3D::tag_value, WT_>::fill_sizes(*static_cast<Halo<delta_, PLFace, MeshType_, WT_, ST_>*>(source), &target_sizes[0]);
            break;
          case pl_polyhedron:
            HaloControl<Dim3D::tag_value, WT_>::fill_sizes(*static_cast<Halo<delta_, PLPolyhedron, MeshType_, WT_, ST_>*>(source), &target_sizes[0]);
            break;
          default:
            throw InternalError("Undefined polytope level!");
        }

        Geometry::MeshPart<Geometry::ConformalMesh<Shape::Hypercube<3> > > result(&target_sizes[0]);

        switch(source->get_level())
        {
          case pl_vertex:
            HaloControl<Dim3D::tag_value, WT_>::fill_target_set(*static_cast<Halo<delta_, PLVertex, MeshType_, WT_, ST_>*>(source), result);
            break;
          case pl_edge:
            HaloControl<Dim3D::tag_value, WT_>::fill_target_set(*static_cast<Halo<delta_, PLEdge, MeshType_, WT_, ST_>*>(source), result);
            break;
          case pl_face:
            HaloControl<Dim3D::tag_value, WT_>::fill_target_set(*static_cast<Halo<delta_, PLFace, MeshType_, WT_, ST_>*>(source), result);
            break;
          case pl_polyhedron:
            HaloControl<Dim3D::tag_value, WT_>::fill_target_set(*static_cast<Halo<delta_, PLPolyhedron, MeshType_, WT_, ST_>*>(source), result);
            break;
          default:
            throw InternalError("Undefined polytope level!");
        }

        return result;
      }
    };
  }
}

#endif
