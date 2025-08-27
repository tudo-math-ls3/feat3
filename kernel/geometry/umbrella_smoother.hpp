// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <exception>
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
    * \brief Umbrella mesh smoother with both uniform and scale-dependent neighbor averaging
    * For more information see \cite{bray2004notes,nealen2006laplacian}
    *
    * \tparam MeshType_
    * The type of the mesh that is to be smoothed, e.g. ConformalMesh
    *
    * \attention
    * This class should only be used in serial applications, i.e. on unpartitioned meshes!
    * If you need to distort a partitioned mesh, use the derived DistributedUmbrellaSmoother class
    * template instead!
    *
    * \attention
    * Only use uniform = false if you know what you are doing
    *
    * \author Pia Ritter
    */
    template<typename MeshType_>
    class UmbrellaSmoother
    {
    public:
      typedef typename MeshType_::CoordType CoordType;
      typedef typename MeshType_::VertexSetType VertexSetType;
      static constexpr int shape_dim = MeshType_::shape_dim;
      static constexpr int world_dim = MeshType_::world_dim;

    protected:
      /// the mesh we want to smooth
      MeshType_& _mesh;
      /// vector to store which vertex is at the boundary
      std::vector<int> _boundary_vertices;
      /// vector to store which facet is at the boundary
      std::vector<int> _boundary_facets;
      // optional per-edge weights; empty in serial case
      std::vector<CoordType> _edge_weights;
      // vertex set of the mesh
      VertexSetType& _vtx;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] _mesh
       * A reference to the mesh that is to be smoothed
       *
       */
      explicit UmbrellaSmoother(MeshType_& mesh) :
        _mesh(mesh),
        _boundary_vertices(),
        _boundary_facets(),
        _edge_weights(),
        _vtx(_mesh.get_vertex_set())
      {
      }

      // virtual destructor
      virtual ~UmbrellaSmoother()
      {
      }

      void compile()
      {
        if(_boundary_vertices.empty())
          _calc_boundary_vertices();
      }

      /**
      * \brief Smooths the mesh with an umbrella operator
      *
      * Iteratively moves each inner vertex towards the (possibly weighted)
      * average of its adjacent vertices (Jacobi update).
      * The weights are either uniform or inversely proportional
      * to the edge lengths (w_ij = 1 / |e_ij|).
      *
      * \param[in] uniform
      * true for uniform operator and false for scale dependent operator
      *
      * \param[in] max_iter
      * Maximum number of smoothing iterations
      *
      * \param[in] tol
      * Iteration stops if the largest vertex displacement between
      * two steps falls below this threshold
      *
      * \note
      * This function only moves inner vertices and leaves boundary vertices untouched.
      */
      void smooth(Index max_iter, CoordType tol, bool uniform = true)
      {
        if(_boundary_vertices.empty())
          _calc_boundary_vertices();

        const auto& verts_at_edge = _mesh.template get_index_set<1, 0>();
        const Index num_edges = _mesh.get_num_entities(1);
        const Index num_vtx = _mesh.get_num_vertices();

        // make sure |e_ij| != 0
        const CoordType eps = Math::sqrt(Math::eps<CoordType>());

        // sum of neighbors for every vertex
        std::vector<Tiny::Vector<CoordType, world_dim>> sum_neigh(num_vtx);

        // number of neighbors for every vertex
        std::vector<CoordType> deg(num_vtx, CoordType(0));

        for (Index iter = 0; iter < max_iter; ++iter)
        {
          // format vectors
          for (Index i(0); i < num_vtx; ++i)
          {
            sum_neigh[i].format();
            deg[i] = CoordType(0);
          }

          // loop over all edges
          for (Index e(0); e < num_edges; ++e)
          {
            // find vertices on edge e
            const Index i = verts_at_edge(e,0);
            const Index j = verts_at_edge(e,1);

            CoordType weight = CoordType(1);

            // calculate edge length for scale dependent smoothing
            if(!uniform)
            {
              weight = CoordType(1) / std::max((_vtx[i] - _vtx[j]).norm_euclid(), eps);
            }

            // multiply by precomputed per-edge weight, if provided
            if (!_edge_weights.empty())
            {
              weight *= _edge_weights[e];
            }

            // add to counter
            if (_boundary_vertices[i] == 0)
            {
              sum_neigh[i].axpy(weight, _vtx[j]);
              deg[i] += weight;
            }
            if (_boundary_vertices[j] == 0)
            {
              sum_neigh[j].axpy(weight, _vtx[i]);
              deg[j] += weight;
            }
          }

          // synchronize over all processes, if necessary
          _synchronize(sum_neigh, deg);

          CoordType max_change = 0;

          // loop over vertices to update location
          for (Index i(0); i < num_vtx; ++i)
          {
            if ((_boundary_vertices[i] == 0) && (std::abs(deg[i]) > eps))
            {
              // compute new vertex coordinates
              auto x_old = _vtx[i];
              _vtx[i] = (CoordType(1) / deg[i]) * sum_neigh[i];

              // update maximum change
              Math::maxi(max_change, (_vtx[i] - x_old).norm_euclid_sqr());
            }
          }

          // calculate maximum over all processes if necessary
          _calc_max(max_change);

          // abort if we reached tolerance
          if (max_change < tol*tol)
          {
            break;
          }
        }
      }

    protected:
      /// determines for each facet whether it is a boundary facet or an inner facet
      virtual void _calc_boundary_facets()
      {
        _boundary_facets.resize(_mesh.get_num_entities(shape_dim - 1), 0);

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

        _boundary_vertices.resize(_mesh.get_num_vertices(), 0);

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

      /// synchronizes the smoother over all processes; is overridden in DistributedUmbrellaSmoother
      virtual void _synchronize(std::vector<Tiny::Vector<CoordType, world_dim>>& , std::vector<CoordType>&)
      {
        // nothing to do here
      }

      /// calculates global maximum over all processes; is overridden in DistributedUmbrellaSmoother
      virtual void _calc_max(CoordType&)
      {
        // nothing to do here
      }
    }; // class UmbrellaSmoother<...>

    /**
     * \brief Umbrella smoother class template
     *
     * This class extends the functionality of the UmbrellaSmoother class by addition MPI-based
     * synchronization to ensure that one obtains a consistently smoothed mesh in a distributed
     * simulation. Note that this class needs to work with a RootMeshNode rather than a "naked"
     * mesh object, because it requires the halos for synchronization, which are stored in the
     * mesh node and not in the mesh object itself.
     *
     * See the documentation of the base-class UmbrellaSmoother for more information
     *
     * \tparam MeshType_
     * The type of the mesh that is to be smoothed, e.g. ConformalMesh
     *
     * \author Pia Ritter
     */
    template<typename MeshType_>
    class DistributedUmbrellaSmoother :
      public UmbrellaSmoother<MeshType_>
    {
    public:
      typedef UmbrellaSmoother<MeshType_> BaseClass;
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

    public:
      /**
      * \brief Constructor
      *
      * \param[in] comm
      * A \resident reference to the communicator that was used for partitioning.
      *
      * \param[in] root_mesh_node
      * A \resident reference to the root mesh node containing the mesh that is to be smoothed.
      */
      explicit DistributedUmbrellaSmoother(const Dist::Comm& comm, RootMeshNodeType& root_mesh_node) :
        BaseClass(*root_mesh_node.get_mesh()),
        _comm(comm),
        _root_mesh_node(root_mesh_node)
      {
        _init_edge_multiplicity();
      }

      // virtual destructor
      virtual ~DistributedUmbrellaSmoother()
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
      /**
      * \brief Computes per-edge weights for distributed smoothing
      *
      * Counts how many partitions share each edge and sets the edge weight to 1/multiplicity.
      */
      void _init_edge_multiplicity()
      {
        if(_comm.size() == 1)
        {
          return;
        }
        const Index num_edges = this->_mesh.get_num_entities(1);

        // initialize edge weights and multiplicity with 1
        this->_edge_weights.assign(num_edges, CoordType(1.0));
        std::vector<int> multiplicity(num_edges, 1);

        // get our halo map
        const auto& halo_map = _root_mesh_node.get_halo_map();

        // loop over all halos
        for (const auto& hm : halo_map)
        {
          const int neighbor_rank = hm.first;
          const std::unique_ptr<MeshPartType>& part_ptr = hm.second;
          if (!part_ptr)
          {
            continue;
          }

          // skip if neighboring rank is our rank
          if (neighbor_rank == _comm.rank())
          {
            continue;
          }

          // get edges that are shared with other halos
          const auto& target_edges = part_ptr->template get_target_set<1>();

          for (Index k(0); k < target_edges.get_num_entities(); ++k)
          {
            // get edge number
            const Index e = target_edges[k];
            // increase multiplicity
            multiplicity[e] += 1;
          }
        }
        for (Index e(0); e < num_edges; ++e)
        {
          // make sure we don't divide by zero
          const int m = std::max(1, multiplicity[e]);
          this->_edge_weights[e] = CoordType(1.0) / CoordType(m);
        }
      }

      // synchronize over all processes, if necessary
      void _synchronize(std::vector<Tiny::Vector<CoordType, world_dim>>& sum_neigh, std::vector<CoordType>& deg) override
      {
        // get our halo map
        const auto& halo_map = _root_mesh_node.get_halo_map();

        // if our halo map is empty or it's only one process -> abort
        if (halo_map.empty() || _comm.size() == 1) return;

        // sum_neigh has dimension world_dim and deg has dimension 1
        const std::size_t stride = std::size_t(world_dim + 1);

        // create buffers
        std::vector<int> ranks;
        std::vector<std::vector<CoordType>> send_bufs, recv_bufs;
        Dist::RequestVector send_reqs, recv_reqs;
        ranks.reserve(halo_map.size());
        send_bufs.resize(halo_map.size());
        recv_bufs.resize(halo_map.size());
        send_reqs.reserve(halo_map.size());
        recv_reqs.reserve(halo_map.size());

        // post receives
        for (auto it = halo_map.begin(); it != halo_map.end(); ++it)
        {
          // get halo rank and mesh part
          const int halo_rank = it->first;
          const auto& part = it->second;
          if (!part)
            continue;
          if (halo_rank == _comm.rank())
            continue;

          // get halo vertices and size
          const auto& halo_vtx  = part->template get_target_set<0>();
          const Index halo_size = halo_vtx.get_num_entities();

          // calculate necessary receive buffer length
          const std::size_t len = std::size_t(halo_size) * stride;

          // allocate buffers
          const std::size_t idx = ranks.size();
          ranks.push_back(halo_rank);
          auto& rbuf = recv_bufs[idx];
          rbuf.resize(len);

          // post receive
          recv_reqs.push_back(_comm.irecv(rbuf.data(), rbuf.size(), halo_rank));
        }

        // post sends
        for (std::size_t idx = 0; idx < ranks.size(); ++idx)
        {
          // get halo rank
          const int halo_rank = ranks[idx];
          auto it = halo_map.find(halo_rank);
          XASSERT(it != halo_map.end());

          // get halo vertices and size
          const auto& halo_vtx  = it->second->template get_target_set<0>();
          const Index halo_size = halo_vtx.get_num_entities();

          // calculate necessary send buffer length
          const std::size_t len = std::size_t(halo_size) * stride;

          // allocate buffers
          auto& sbuf = send_bufs[idx];
          sbuf.resize(len);

          // loop over all halos
          for (Index k(0); k < halo_size; ++k)
          {
            const Index v = halo_vtx[k];
            const std::size_t base = std::size_t(k) * stride;
            // copy sum_neigh and deg in send buffer
            for (int d(0); d < world_dim; ++d)
            {
              sbuf[base + std::size_t(d)] = sum_neigh[v][d];
            }
            sbuf[base + std::size_t(world_dim)] = deg[v];
          }

          // post send
          send_reqs.push_back(_comm.isend(sbuf.data(), sbuf.size(), halo_rank));
        }

        // merge entries
        for (std::size_t idx(0); recv_reqs.wait_any(idx); )
        {
          // get halo rank
          const int halo_rank = ranks.at(idx);
          auto it = halo_map.find(halo_rank);
          XASSERT(it != halo_map.end());

          // get halo vertices and size
          const auto& halo_vtx  = it->second->template get_target_set<0>();
          const Index halo_size = halo_vtx.get_num_entities();
          const auto& rbuf = recv_bufs.at(idx);

          // check if the size is right
          XASSERT(rbuf.size() == std::size_t(halo_size) * stride);

          // loop over all halos
          for (Index k(0); k < halo_size; ++k)
          {
            const Index v = halo_vtx[k];
            const std::size_t base = std::size_t(k) * stride;

            // add up sum_neigh and deg
            for (int d(0); d < world_dim; ++d)
            {
              sum_neigh[v][d] += rbuf[base + std::size_t(d)];
            }
            deg[v] += rbuf[base + std::size_t(world_dim)];
          }
        }

        // wait for all sends
        send_reqs.wait_all();
      }

      // calculate global maximum, if necessary
      void _calc_max(CoordType& global_max) override
      {
         _comm.allreduce(&global_max, &global_max, std::size_t(1), Dist::op_max);
      }
    }; // class DistributedUmbrellaSmoother
  } // namespace Geometry
} // namespace FEAT
