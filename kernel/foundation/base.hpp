#pragma once
#ifndef KERNEL_FOUNDATION_BASE_HH
#define KERNEL_FOUNDATION_BASE_HH 1

#include<kernel/base_header.hpp>

namespace FEAST
{
  /**
   * \brief Foundation namespace
   */
  namespace Foundation
  {
    ///polytope level identifiers
    enum PolytopeLevels
    {
      pl_vertex = 0,
      pl_edge,
      pl_face,
      pl_polyhedron,
      pl_none
    };

    ///spatial dimensions
    enum Dimensions
    {
      dim_1D = 1,
      dim_2D = 2,
      dim_3D = 3
    };

    ///Required number of topologies for RDMeshes
    enum RequiredNumTopologies
    {
      rnt_1D = 2,
      rnt_2D = 4,
      rnt_3D = 6
    };

    ///Tag classes
    struct PLNone
    {
      typedef PLNone SubElementPolytopeType_;

      static const PolytopeLevels tag_value = pl_none;
    };

    struct PLVertex
    {
      typedef PLNone SubElementPolytopeType_;

      static const PolytopeLevels tag_value = pl_vertex;
    };

    struct PLEdge
    {
      typedef PLVertex SubElementPolytopeType_;

      static const PolytopeLevels tag_value = pl_edge;
    };

    struct PLFace
    {
      typedef PLEdge SubElementPolytopeType_;

      static const PolytopeLevels tag_value = pl_face;
    };

    struct PLPolyhedron
    {
      typedef PLFace SubElementPolytopeType_;

      static const PolytopeLevels tag_value = pl_polyhedron;
    };

    struct Dim0D
    {
      typedef Dim0D SubSpaceType_;
    };

    struct Dim1D
    {
      typedef PLEdge ElementPolytopeType_;
      typedef Dim0D SubSpaceType_;

      static const Dimensions tag_value = dim_1D;
      static const RequiredNumTopologies required_num_topologies = rnt_1D;
    };

    struct Dim2D
    {
      typedef PLFace ElementPolytopeType_;
      typedef Dim1D SubSpaceType_;

      static const Dimensions tag_value = dim_2D;
      static const RequiredNumTopologies required_num_topologies = rnt_2D;
    };

    struct Dim3D
    {
      typedef PLPolyhedron ElementPolytopeType_;
      typedef Dim2D SubSpaceType_;

      static const Dimensions tag_value = dim_3D;
      static const RequiredNumTopologies required_num_topologies = rnt_3D;
    };

#ifdef FEAST_HAVE_PARMETIS
    struct ParmetisModePartKway
    {
    };

    struct ParmetisModePartGeomKway
    {
    };
#endif

  }
}
#endif
