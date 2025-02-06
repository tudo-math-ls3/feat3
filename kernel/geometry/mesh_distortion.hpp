// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/geometry/factory.hpp>
#include <kernel/geometry/mesh_permutation.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/intern/facet_neighbors.hpp>
#include <kernel/geometry/intern/standard_index_refiner.hpp>
#include <kernel/geometry/intern/standard_vertex_refiner.hpp>
#include <kernel/geometry/index_calculator.hpp>
#include <kernel/util/dist.hpp>

#include <kernel/util/random.hpp>

#include <kernel/util/tiny_algebra.hpp>


namespace FEAT
{
  namespace Geometry
  {
    /**
     * \brief (Unpartitioned) mesh distortion class template
     *
     * \tparam MeshType_
     * The type of the mesh that is to be distorted, e.g. ConformalMesh
     *
     * \attention
     * This class should only be used in serial applications, i.e. on unpartitioned meshes!
     * If you need to distort a partitioned mesh, use the derived DistributedMeshDistortion class
     * template instead!
     *
     * \author Chantal Jahner
     */
    template<typename MeshType_>
    class MeshDistortion
    {
    public:
      typedef typename MeshType_::CoordType CoordType;
      static constexpr int shape_dim = MeshType_::shape_dim;
      static constexpr int world_dim = MeshType_::world_dim;

    protected:
      /// the mesh we want to distort
      MeshType_& _mesh;
      /// random number generator for random distortion
      Random _rand;
      /// vector to store which vertex is at the boundary
      std::vector<int> _boundary_vertices;
      /// vector to story which facet is at the boundary
      std::vector<int> _boundary_facets;
      /// vector to store the length of the shortest edge belonging to each vertex
      std::vector<CoordType> _shortest_edge;
      /// have a guess...
      const CoordType _pi_val = Math::pi<CoordType>();

    public:
      /**
       * \brief Constructor
       *
       * \param[in] _mesh
       * A reference to the mesh that is to be distorted
       */
      explicit MeshDistortion(MeshType_& mesh) :
        _mesh(mesh),
        _rand(),
        _boundary_vertices(_mesh.get_num_vertices(), 0),
        _boundary_facets(_mesh.get_num_entities(shape_dim - 1), 0),
        _shortest_edge(_mesh.get_num_vertices(), Math::huge<CoordType>())
      {
      }
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
      explicit MeshDistortion(MeshType_& mesh, Random::SeedType seed) :
        _mesh(mesh),
        _rand(seed),
        _boundary_vertices(_mesh.get_num_vertices(), 0),
        _boundary_facets(_mesh.get_num_entities(shape_dim - 1), 0),
        _shortest_edge(_mesh.get_num_vertices(), Math::huge<CoordType>())
      {
      }

      virtual ~MeshDistortion()
      {
      }

      /**
       * \brief Distorts the mesh uniformly
       *
       * Moves each inner vertex on an arbitrary point on a circle (2D) or sphere (3D) around the original vertex
       *
       * \param[in] rad
       * Radius of the circle or sphere
       *
       * \note
       * This function only moves inner vertices and leaves boundary vertices untouched.
       */
      void distort_uniform(CoordType rad)
      {
        _calc_boundary_vertices();

        auto& vtx = _mesh.get_vertex_set();

        if(world_dim == 2)
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
        else if (world_dim == 3)
        {
          for (Index i(0); i < vtx.get_num_vertices(); ++i)
          {
            if(_boundary_vertices[i] == 0)
            {
              CoordType theta(0), phi(0);
              _rand >> theta;
              _rand >> phi;
              theta *= _pi_val;
              phi *= CoordType(2) * _pi_val;
              vtx[i][0] += rad * Math::sin(theta) * Math::cos(phi);
              vtx[i][1] += rad * Math::sin(theta) * Math::sin(phi);
              vtx[i][2] += rad * Math::cos(theta);
            }
          }
        }

        // synchronize over all processes, if necessary
        _synchronize();
      }

      /**
       * \brief Distorts the mesh according to the shortest  edge
       *
       * Moves each inner vertex on an arbitrary point on a circle (2D) or sphere (3D) around the original vertex
       * The radius of the circle/sphere is determined by the shortest edge in the mesh:
       *
       * radius = scale * shortest_edge
       *
       * \param[in] scale
       * factor to scale the length of the shortest edge to determine the radius
       *
       * \note
       * This function only moves inner vertices and leaves boundary vertices untouched.
       */
      void distort_shortest_edge_uniform(CoordType scale = CoordType(0.1))
      {
        _calc_boundary_vertices();
        _calc_shortest_edge();
        CoordType min = *std::min_element(_shortest_edge.begin(), _shortest_edge.end());
        distort_uniform(min * scale);
      }

      /**
       * \brief Distorts the mesh according to the shortest local edge
       *
       * Moves each inner vertex on an arbitrary point on a circle (2D) or sphere (3D) around the original vertex
       * The radius of the circle/sphere is determined by the shortest edge adjacent to the vertex:
       * radius = scale * length_shortest_edge
       *
       * \param[in] scale
       * factor to scale the length of the shortest edge to determine the radius
       *
       * \note
       * This function only moves inner vertices and leaves boundary vertices untouched.
       */
      void distort_shortest_edge_local(CoordType scale = CoordType(0.1))
      {
        _calc_boundary_vertices();
        _calc_shortest_edge();

        CoordType rad;

        auto& vtx = _mesh.get_vertex_set();

        if(world_dim == 2)
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
        else if (world_dim == 3)
        {
          for (Index i(0); i < vtx.get_num_vertices(); ++i)
          {
            if(_boundary_vertices[i] == 0)
            {
              rad = scale * _shortest_edge[i];
              CoordType theta(0), phi(0);
              _rand >> theta;
              _rand >> phi;
              theta *= _pi_val;
              phi *= CoordType(2) * _pi_val;
              vtx[i][0] += rad * Math::sin(theta) * Math::cos(phi);
              vtx[i][1] += rad * Math::sin(theta) * Math::sin(phi);
              vtx[i][2] += rad * Math::cos(theta);
            }
          }
        }

        // synchronize over all processes, if necessary
        _synchronize();
      }

    protected:
      /// determines for each facet whether it is a boundary facet or an inner facet
      virtual void _calc_boundary_facets()
      {
        const auto& vtx_neighbors = _mesh.get_neighbors();

        const auto& faces_at_elem = _mesh.template get_index_set<shape_dim, shape_dim - 1>();

        int num_faces = faces_at_elem.num_indices;

        Index num_elements(_mesh.get_num_elements());

        // loop over all elements
        for (Index i = 0; i < num_elements; ++i)
        {
          // loop over faces of each element
          for (int j = 0; j < num_faces; ++j)
          {
            // check whether there is a neighbor at that edge
            if (vtx_neighbors[i][j] > num_elements)
            {
              // if not: facet is at the boundary
              _boundary_facets[faces_at_elem[i][j]] = 1;
            }
          }
        }
      }

      /// determines for each vertex whether it is a boundary vertex or an inner vertex
      virtual void _calc_boundary_vertices()
      {
        _calc_boundary_facets();
        const auto& verts_at_face = _mesh.template get_index_set<shape_dim - 1, 0>();

        int vert_per_face = verts_at_face.num_indices;

        // loop over all facets
        for (Index i = 0; i < _mesh.get_num_entities(shape_dim - 1); ++i)
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

      /// determines for each vertex the length of the shortest adjacent edge
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

      /// synchronizes the distortion over all processes; is overridden in DistributedMeshDistortion
      virtual void _synchronize()
      {
        // nothing to do here
      }
    }; // class MeshDistortion<...>

    /**
     * \brief Distributed mesh distortion class template
     *
     * This class extends the functionality of the MeshDistortion class by addition MPI-based
     * synchronization to ensure that one obtains a consistently distorted mesh in a distributed
     * simulation. Note that this class needs to work with a RootMeshNode rather than a "naked"
     * mesh object, because it requires the halos for synchronization, which are stored in the
     * mesh node and not in the mesh object itself.
     *
     * See the documentation of the base-class MeshDistortion for more information a list of
     * available distortion algorithms.
     *
     * \tparam MeshType_
     * The type of the mesh that is to be distorted, e.g. ConformalMesh
     *
     * \author Peter Zajac
     */
    template<typename MeshType_>
    class DistributedMeshDistortion :
      public MeshDistortion<MeshType_>
    {
    public:
      typedef MeshDistortion<MeshType_> BaseClass;
      typedef MeshType_ MeshType;
      typedef MeshPart<MeshType> MeshPartType;
      typedef RootMeshNode<MeshType> RootMeshNodeType;

      static constexpr int world_dim = MeshType::world_dim;
      static constexpr int facet_dim = MeshType::shape_dim - 1;

      typedef typename BaseClass::CoordType CoordType;

    protected:
      /// a reference to our communicator
      const Dist::Comm& _comm;
      /// a reference to our root mesh node
      RootMeshNodeType& _root_mesh_node;
      /// a vector of vertex owner ranks
      std::vector<int> _vertex_owners;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] comm
       * A \resident reference to the communicator that was used for partitioning.
       *
       * \param[in] root_mesh_node
       * A \resident reference to the root mesh node containing the mesh that is to be distorted.
       */
      explicit DistributedMeshDistortion(const Dist::Comm& comm, RootMeshNodeType& root_mesh_node) :
        BaseClass(*root_mesh_node.get_mesh()),
        _comm(comm),
        _root_mesh_node(root_mesh_node),
        _vertex_owners()
      {
      }
      /**
       * \brief Constructor
       *
       * \param[in] comm
       * A \resident reference to the communicator that was used for partitioning.
       *
       * \param[in] root_mesh_node
       * A \resident reference to the root mesh node containing the mesh that is to be distorted.
       *
       * \param[in] seed
       * The seed for the internal random number generator.
       */
      explicit DistributedMeshDistortion(const Dist::Comm& comm, RootMeshNodeType& root_mesh_node, Random::SeedType seed) :
        BaseClass(*root_mesh_node.get_mesh(), seed),
        _comm(comm),
        _root_mesh_node(root_mesh_node),
        _vertex_owners()
      {
      }

      /// virtual destructor
      virtual ~DistributedMeshDistortion()
      {
      }

    protected:
      /// determines for each facet whether it is a boundary facet or an inner facet
      virtual void _calc_boundary_facets() override
      {
        // call base class to determine all boundary facets
        BaseClass::_calc_boundary_facets();

        // Each process has already calculated its local boundary facets, so now we need to loop
        // over all halos and remove all halo facets from the set of boundary facets, because
        // halo facets are by definition not a part of the actual domain boundary.

        // loop over all halos and kick out all halo facets
        const auto& halo_map = _root_mesh_node.get_halo_map();
        for(const auto& v : halo_map)
        {
          // loop over all facets of the halo and remove them from the boundary set
          const TargetSet& halo_facets = v.second->template get_target_set<facet_dim>();
          for(Index i(0); i < halo_facets.get_num_entities(); ++i)
            this->_boundary_facets.at(halo_facets[i]) = 0;
        }
      }

      /// determines for each vertex whether it is a boundary vertex or an inner vertex
      virtual void _calc_shortest_edge() override
      {
        // call base class to determine the shortest edges for this patch
        BaseClass::_calc_shortest_edge();

        // Each process has already calculated the shortest edge for each of its vertices, so now
        // we need to loop over all halos and exchange the shortest edge length for each halo
        // vertex with the corresponding neighbor process and determine for each vertex whether
        // that neighbor has found a shorter edge for it.

        // get our halo map
        const auto& halo_map = _root_mesh_node.get_halo_map();
        const std::size_t num_halos = halo_map.size();

        // nothing to do here?
        if(halo_map.empty())
          return;

        // create buffers
        std::vector<int> ranks;
        std::vector<std::vector<CoordType>> send_bufs, recv_bufs;
        Dist::RequestVector send_reqs, recv_reqs;
        ranks.reserve(num_halos);
        send_bufs.resize(num_halos);
        recv_bufs.resize(num_halos);
        send_reqs.reserve(num_halos);
        recv_reqs.reserve(num_halos);

        // loop over all our halos/neighbors
        for(auto it = halo_map.begin(); it != halo_map.end(); ++it)
        {
          // get halo and the number of halo vertices
          const int halo_rank = it->first;
          const TargetSet& halo_vtx = it->second->template get_target_set<0>();
          const Index halo_size = halo_vtx.get_num_entities();

          // allocate buffers
          auto& sbuf = send_bufs.at(ranks.size());
          auto& rbuf = recv_bufs.at(ranks.size());
          sbuf.resize(halo_size);
          rbuf.resize(halo_size);
          ranks.push_back(halo_rank);

          // post receive
          recv_reqs.push_back(_comm.irecv(rbuf.data(), rbuf.size(), halo_rank));

          // compute source buffer = shortest edge for each halo vertex
          for(Index k(0); k < halo_size; ++k)
            sbuf.at(k) = this->_shortest_edge.at(halo_vtx[k]);

          // post send
          send_reqs.push_back(_comm.isend(sbuf.data(), sbuf.size(), halo_rank));
        }

        // process all pending receives
        for(std::size_t idx(0u); recv_reqs.wait_any(idx); )
        {
          // get the halo rank
          const int halo_rank = ranks.at(idx);

          // find the halo mesh part
          auto it = halo_map.find(halo_rank);
          XASSERT(it != halo_map.end());

          // get the halo vertex target set
          const TargetSet& halo_vtx = it->second->template get_target_set<0>();
          const Index halo_size = halo_vtx.get_num_entities();
          auto& rbuf = recv_bufs.at(idx);

          // update the shortest edge length for each halo vertex
          for(Index k(0); k < halo_size; ++k)
            this->_shortest_edge.at(halo_vtx[k]) = Math::min(this->_shortest_edge.at(halo_vtx[k]), rbuf.at(k));
        }

        // wait for all sends to finish
        send_reqs.wait_all();
      }

      /// determines for each vertex the length of the shortest adjacent edge
      virtual void _calc_vertex_owners()
      {
        if(!_vertex_owners.empty())
          return;

        // initially assume that we own all our vertices
        _vertex_owners.resize(_root_mesh_node.get_mesh()->get_num_vertices(), _comm.rank());

        // loop over all our halos/neighbors
        const auto& halo_map = _root_mesh_node.get_halo_map();
        for(auto it = halo_map.begin(); it != halo_map.end(); ++it)
        {
          const int halo_rank = it->first;
          const TargetSet& halo_vtx = it->second->template get_target_set<0>();

          // update vertex owner ranks
          for(Index i(0); i <  halo_vtx.get_num_entities(); ++i)
          {
            _vertex_owners[halo_vtx[i]] = Math::min(_vertex_owners[halo_vtx[i]], halo_rank);
          }
        }
      }


      // synchronize over all processes, if necessary
      virtual void _synchronize() override
      {
        // Each process has already disturbed its local mesh, so now we need to ensure that the
        // coordinates of all vertices, which are shared by more than one process, have consistent
        // (= identical) coordinates on all sharing process. We do this by employing the convention
        // that the coordinates of a shared vertex are defined by the process with the smallest rank
        // that shares this particular vertex and we call that process the 'owner' of that vertex.
        // With this definition, our approach is as follows: Loop over all halos/neighbors and
        // determine whether we have small rank than that neighbor and if...
        // yes: then we collect the vertex coordinates for this halo and send these coordinates
        //      to our neighbor.
        // no:  then we post a receive request for the vertex coordinates and once that receive
        //      request is fulfilled, we loop over all vertices in that halo, check for each
        //      have vertex if that neighbor is the owner and if yes, then we overwrite our
        //      vertex coordinates with the ones we receives from our (owner) neighbor.

        // get our halo map
        const auto& halo_map = _root_mesh_node.get_halo_map();
        const std::size_t num_halos = halo_map.size();

        // nothing to do here?
        if(halo_map.empty())
          return;

        // compute the vertex owner ranks
        _calc_vertex_owners();

        auto& vertex_set = _root_mesh_node.get_mesh()->get_vertex_set();

        const int my_rank = _comm.rank();

        // create buffers
        std::vector<int> ranks;
        std::vector<std::vector<CoordType>> buffers;
        Dist::RequestVector requests;
        ranks.reserve(num_halos);
        buffers.resize(num_halos);
        requests.reserve(num_halos);

        // loop over all our halos/neighbors
        for(auto it = halo_map.begin(); it != halo_map.end(); ++it)
        {
          // get halo and buffer size
          const int halo_rank = it->first;
          const TargetSet& halo_vtx = it->second->template get_target_set<0>();
          const Index halo_size = halo_vtx.get_num_entities();
          const std::size_t buf_size = halo_size * Index(world_dim);

          // allocate buffers
          auto& buf = buffers.at(ranks.size());
          buf.resize(buf_size);
          ranks.push_back(halo_rank);

          // should we send or receive the vertex coordinates?
          if(my_rank < halo_rank)
          {
            // copy vertex coordinates into buffer
            for(Index i(0), j(0); i < halo_size; ++i)
            {
              const auto& vtx = vertex_set[halo_vtx[i]];
              for(int k(0); k < world_dim; ++j, ++k)
                buf[j] = vtx[k];
            }

            // post send
            requests.push_back(_comm.isend(buf.data(), buf.size(), halo_rank));
          }
          else
          {
            // post receive
            requests.push_back(_comm.irecv(buf.data(), buf.size(), halo_rank));
          }
        }

        // process all pending receives
        for(std::size_t idx(0u); requests.wait_any(idx); )
        {
          // get the halo rank
          const int halo_rank = ranks.at(idx);

          // was this a send-request?
          if(my_rank < halo_rank)
            continue; // nothing to do after a send

          // find the halo mesh part
          auto it = halo_map.find(halo_rank);
          XASSERT(it != halo_map.end());

          // get the halo vertex target set
          const TargetSet& halo_vtx = it->second->template get_target_set<0>();
          const Index halo_size = halo_vtx.get_num_entities();
          auto& buf = buffers.at(idx);

          // loop over all halo vertices
          for(Index i(0), j(0); i < halo_size; ++i, j += world_dim)
          {
            // is the sender process the actual owner of the vertex?
            if(_vertex_owners[halo_vtx[i]] == halo_rank)
            {
              // yup, so accept the owner's vertex coordinates
              auto& vtx = vertex_set[halo_vtx[i]];
              for(int k(0); k < world_dim; ++k)
                vtx[k] = buf[j+Index(k)];
            }
          }
        }
      }
    }; // class DistributedMeshDistortion
  } // namespace Geometry
} // namespace FEAT
