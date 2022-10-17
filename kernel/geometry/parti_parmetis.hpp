// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2022 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_PARTI_PARMETIS_HPP
#define KERNEL_GEOMETRY_PARTI_PARMETIS_HPP 1

#include <kernel/base_header.hpp>

#if defined(FEAT_HAVE_PARMETIS) || defined(DOXYGEN)
#include <kernel/util/dist.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/adjacency/coloring.hpp>
#include <kernel/geometry/conformal_mesh.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /**
     * \brief ParMETIS mesh/graph partitioner backend
     *
     * \author Peter Zajac
     */
    class PartiParMETIS
    {
      // maximum number of processes to use by ParMETIS
      static constexpr int max_procs = 1000;
      // minimum number of elements per process for ParMETIS
      static constexpr int min_elems = 1000;

    private:
      /// our main communicator
      const Dist::Comm& _comm;
      /// a sub-communicator for the partitioner
      Dist::Comm _sub_comm;
      /// the total number of elements
      Index _num_elems;
      /// the desired number of partitions
      Index _num_parts;
      /// first element for this process
      Index _first_elem;
      /// number of elements for this process
      Index _num_local_elems;
      /// the adjacency sub-graph for this process
      Adjacency::Graph _subgraph;
      /// element midpoints
      std::vector<Real> _midpoints;
      /// element partitioning
      std::vector<Index> _parts;
      /// the element coloring
      Adjacency::Coloring _coloring;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] comm
       * A \resident reference to the communicator to be used by the partitioner.
       * In a recursive partitioning setting, this has to be the progeny communicator.
       *
       * Note that this class might create a sub-communicator of \p comm internally if there are
       * too many ranks in the communicator for the given number of mesh elements or desired number
       * of partitions to improve the computation vs communication trade-off.
       */
      explicit PartiParMETIS(const Dist::Comm& comm);

      /// virtual destructor
      virtual ~PartiParMETIS();

      /**
       * \brief Executes the ParMETIS graph partitioner
       *
       * \param[in] faces_at_elem
       * A \transient reference to the faces-at-element index set of the mesh that is to be partitioned.
       * The face-dimension of the index set can be chosen freely, but typically in the finite element
       * case this is the vertices-at-element index set of the mesh.
       *
       * \param[in] verts_at_elem
       * A \transient reference to the vertices-at-element index set of the mesh. This index set
       * is only used to compute the element midpoint coordinates for the partitioner, but it is
       * not used for the computation of the adjacency pattern.
       *
       * \param[in] vertices
       * A \transient reference to the vertex set of the mesh.
       *
       * \param[in] num_parts
       * The desired number of partitions to create. Must be > 0, but it does not necessarily need
       * to be equal to the number of processes in the communicator.
       *
       * \returns
       * \c true, if the partitioning was successful, or \c false, if some error occurred.
       */
      template<int nf_, int nv_, typename VertexSet_>
      bool execute(const IndexSet<nf_>& faces_at_elem, const IndexSet<nv_>& verts_at_elem,
        const VertexSet_& vertices, const Index num_parts)
      {
        Adjacency::Graph graph(Adjacency::RenderType::as_is, faces_at_elem);
        _create_subgraph(graph, num_parts);
        _compute_midpoints(verts_at_elem, vertices);
        return _execute();
      }

      /**
       * \brief Builds and returns the elements-at-rank graph representing the partitioning.
       */
      Adjacency::Graph build_elems_at_rank() const
      {
        return _coloring.create_partition_graph();
      }

    protected:
      /// auxiliary function: builds a hypergraph from the adjacency graph
      void _create_subgraph(const Adjacency::Graph& faces_at_elem, const Index num_parts);
      /// auxiliary function: executes the partitioner
      bool _execute();
      /// auxiliary function: applies the actual ParMETIS call
      bool _apply_parmetis();
      /// auxiliary function: gathers the partitioned coloring onto the root process
      bool _gather_coloring();
      /// auxiliary function: broadcasts the coloring from the root process
      bool _broadcast_coloring(bool metis_ok);
      /// auxiliary function: compute element midpoints from vertex coordinates
      template<int nv_, typename VertexSet_>
      void _compute_midpoints(const IndexSet<nv_>& verts_at_elem, const VertexSet_& vertices)
      {
        typename VertexSet_::VertexType vx;
        const Real sc = Real(1) / Real(verts_at_elem.num_indices);
        static constexpr int nc = VertexSet_::num_coords;
        _midpoints.reserve(_num_local_elems * Index(nc));
        for(Index i(0); i < _num_local_elems; ++i)
        {
          vx.format();
          for(int j(0); j < verts_at_elem.num_indices; ++j)
            vx += vertices[verts_at_elem(i,j)];
          for(int k(0); k < nc; ++k)
            _midpoints.push_back(Real(vx[k]) * sc);
        }
      }
    }; // class PartiParMETIS
  } // namespace Geometry
} // namespace FEAT

#endif // FEAT_HAVE_PARMETIS
#endif // KERNEL_GEOMETRY_PARTI_PARMETIS_HPP
