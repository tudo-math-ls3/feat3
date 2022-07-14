// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_MESH_DISTORTION_HPP
#define KERNEL_GEOMETRY_MESH_DISTORTION_HPP 1

// includes, FEAT
#include <kernel/geometry/factory.hpp>
#include <kernel/geometry/mesh_permutation.hpp>
#include <kernel/geometry/intern/facet_neighbors.hpp>
#include <kernel/geometry/intern/standard_index_refiner.hpp>
#include <kernel/geometry/intern/standard_vertex_refiner.hpp>
#include <kernel/geometry/index_calculator.hpp>

#include <kernel/util/random.hpp>

#include <kernel/util/tiny_algebra.hpp>


namespace FEAT
{
  namespace Geometry
  {
    /**
     * \brief Mesh distortion class template
     *
     * \tparam MeshType_
     * The type of the mesh that is to be distorted, i.e. ConformalMesh
     *
     * \author Chantal Jahner
     */
    template<typename MeshType_>
    class MeshDistortion
    {
    protected:
      /// the mesh we want to distort
      MeshType_& _mesh;
      typedef typename MeshType_::CoordType CoordType;
      static constexpr int _shape_dim = MeshType_::shape_dim;
      // random number generator for random distortion
      Random _rand;
      // vector to store which vertex is at the boundary
      std::vector<int> _boundary_vertices;
      //vector to story which facet is at the boundary
      std::vector<int> _boundary_facets;
      // vector to store the length of the shortest edge belonging to each vertex
      std::vector<CoordType> _shortest_edge;
      const CoordType _pi_val = Math::pi<CoordType>();

    public:
      /**
       * \brief Constructor
       *
       * \param[in] _mesh
       * A reference to the mesh that is to be distorted
       *
       * \param[in] _seed
       * The seed for the random number generator
       *
       */
      explicit MeshDistortion(MeshType_& mesh, Random::SeedType seed = Random::def_seed) :
        _mesh(mesh),
        _rand(seed),
        _boundary_vertices(_mesh.get_num_vertices(), 0),
        _boundary_facets(_mesh.get_num_entities(_shape_dim - 1), 0),
        _shortest_edge(_mesh.get_num_vertices(), Math::huge<CoordType>())
      {
      }

      virtual ~MeshDistortion()
      {
      }

      /**
       * \brief Distorts the mesh uniformly
       *
       * Moves each vertex on an arbitrary point on a circle (2D) or sphere (3D) around the original vertex
       *
       * \param[in] rad
       * Radius of the circle or sphere
       */
      void distort_uniform(CoordType rad)
      {
        _calc_boundary_vertices();

        auto& vtx = _mesh.get_vertex_set();

        if(_shape_dim == 2)
        {
          for (Index i(0); i < vtx.get_num_vertices(); ++i)
          {
            if(_boundary_vertices[i] == 0)
            {
              CoordType t(0);
              _rand >> t;
              t *= CoordType(2) * _pi_val;
              vtx[i][0] += rad * Math::cos(t);
              vtx[i][1] += rad * Math::sin(t);
            }
          }
        }
        else if (_shape_dim == 3)
        {
          for (Index i(0); i < vtx.get_num_vertices(); ++i)
          {
            if(_boundary_vertices[i] == 0)
            {
              CoordType theta(0), phi(0);
              _rand >> theta;
              _rand >> phi;
              theta *= CoordType(1) * _pi_val;
              phi *= CoordType(2) * _pi_val;
              vtx[i][0] += rad * Math::sin(theta) * Math::cos(phi);
              vtx[i][1] += rad * Math::sin(theta) * Math::sin(phi);
              vtx[i][2] += rad * Math::cos(theta);
            }
          }
        }
      }

      /**
       * \brief Distorts the mesh according to the shortest  edge
       *
       * Moves each vertex on an arbitrary point on a circle (2D) or sphere (3D) around the original vertex
       * The radius of the circle/sphere is determined by the shortest edge in the mesh:
       *
       * radius = scale * shortest_edge
       *
       * \param[in] scale
       * factor to scale the length of the shortest edge to determine the radius
       */
      void distort_shortest_edge_uniform(CoordType scale = CoordType(0.25))
      {
        _calc_boundary_vertices();
        _calc_shortest_edge();
        CoordType min = *std::min_element(_shortest_edge.begin(), _shortest_edge.end());
        distort_uniform(min * scale);
      }

      /**
       * \brief Distorts the mesh according to the shortest local edge
       *
       * Moves each vertex on an arbitrary point on a circle (2D) or sphere (3D) around the original vertex
       * The radius of the circle/sphere is determined by the shortest edge adjacent to the vertex:
       * radius = scale * length_shortest_edge
       *
       * \param[in] scale
       * factor to scale the length of the shortest edge to determine the radius
       */
      void distort_shortest_edge_local(CoordType scale = CoordType(0.25))
      {
        _calc_boundary_vertices();
        _calc_shortest_edge();

        CoordType rad;

        auto& vtx = _mesh.get_vertex_set();

        if(_shape_dim == 2)
        {
          for (Index i(0); i < vtx.get_num_vertices(); ++i)
          {
            if(_boundary_vertices[i] == 0)
            {
              rad = scale * _shortest_edge[i];
              CoordType t(0);
              _rand >> t;
              t *= CoordType(2) * _pi_val;
              vtx[i][0] += rad * Math::cos(t);
              vtx[i][1] += rad * Math::sin(t);
            }
          }
        }
        else if (_shape_dim == 3)
        {
          for (Index i(0); i < vtx.get_num_vertices(); ++i)
          {
            if(_boundary_vertices[i] == 0)
            {
              rad = scale * _shortest_edge[i];
              CoordType theta(0), phi(0);
              _rand >> theta;
              _rand >> phi;
              theta *= CoordType(1) * _pi_val;
              phi *= CoordType(2) * _pi_val;
              vtx[i][0] += rad * Math::sin(theta) * Math::cos(phi);
              vtx[i][1] += rad * Math::sin(theta) * Math::sin(phi);
              vtx[i][2] += rad * Math::cos(theta);
            }
          }
        }
      }

    protected:
      virtual void _calc_boundary_facets()
      {
        const auto& vtx_neighbors = _mesh.get_neighbors();

        const auto& faces_at_elem = _mesh.template get_index_set<_shape_dim, _shape_dim - 1>();

        int num_faces = faces_at_elem.num_indices;

        Index num_elements(_mesh.get_num_elements());

        // loop over all elements
        for (Index i = 0; i < num_elements; ++i)
        {
          // loop over faces of each element
          for (int j = 0; j < num_faces; ++j)
          {
            //check wether there is a neighbor at that edge
            if (vtx_neighbors[i][j] > num_elements)
            {
              // if not: facet is at the boundary
              _boundary_facets[faces_at_elem[i][j]] = 1;
            }
          }
        }
      }

      virtual void _calc_boundary_vertices()
      {
        _calc_boundary_facets();
        const auto& verts_at_face = _mesh.template get_index_set<_shape_dim - 1, 0>();

        int vert_per_face = verts_at_face.num_indices;

        // loop over all facets
        for (Index i = 0; i < _mesh.get_num_entities(_shape_dim - 1); ++i)
        {
          // check if facet is at boundary
          if(_boundary_facets[i])
          {
            // loop over adjacent vertices
            for (int j = 0; j < vert_per_face; ++j)
            {
              _boundary_vertices[verts_at_face[i][j]] = 1;
            }
          }
        }
      }

      // for each vertex: get the length of the shortest edge adjacent to the vertex
      virtual void _calc_shortest_edge()
      {
        auto& vtx = _mesh.get_vertex_set();

        const auto& verts_at_edg = _mesh.template get_index_set<1, 0>();

        // loop over all edges
        for (Index i = 0; i < verts_at_edg.get_num_entities(); ++i)
        {
          Index ind_vert1 = verts_at_edg[i][0];
          Index ind_vert2 = verts_at_edg[i][1];
          CoordType edge_length = (vtx[ind_vert1] - vtx[ind_vert2]).norm_euclid();

          // update the two nodes belonging to the edge
          _shortest_edge[ind_vert1] = Math::min(_shortest_edge[ind_vert1], edge_length);
          _shortest_edge[ind_vert2] = Math::min(_shortest_edge[ind_vert2], edge_length);
        }
      }

    }; // class MeshDistortion<...>
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_CONFORMAL_MESH_HPP
