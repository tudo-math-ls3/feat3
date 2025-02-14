// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>

#if defined(FEAT_HAVE_ZOLTAN) || defined(DOXYGEN)
#include <kernel/util/dist.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/adjacency/coloring.hpp>
#include <kernel/geometry/conformal_mesh.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /**
     * \brief Zoltan hypergraph partitioner backend
     *
     * This class implements a backend for the Zoltan partitioner from the Sandia National Laboratories.
     * This class only makes use for Zoltan's homebrew Hypergraph partitioners and does not use any
     * of its other partitioning methods, although these could be integrated into this class later,
     * if desired/required. This class is used by the Control::Domain::PartiDomainControl class
     * and its functionality can be used by specifying the tag <c>zoltan</c> as a parameter to the
     * <c>--parti-type</c> command line option.
     *
     * Note that this class is only available if FEAT is configured and linked against the Zoltan
     * third-party library, of course.
     *
     * \todo figure out max_procs/min_elems
     * \todo allow the user to configure more Zoltan parameters
     *
     * \author Peter Zajac
     */
    class PartiZoltan
    {
    public:
      /// internal data class to be used by callback functions
      class Hypergraph
      {
      public:
        /// the first element of this process
        Index first_elem;
        /// the sub-graph for this process
        Adjacency::Graph sub_graph;
        /// the weights for the elements of this sub-graph
        std::vector<float> weights;

        Hypergraph() :
          first_elem(0u),
          sub_graph(),
          weights()
        {
        }
      }; // class Hypergraph

    private:
      /// our main communicator
      const Dist::Comm& _comm;
      /// a sub-communicator for the partitioner
      Dist::Comm _zoltan_comm;
      /// maximum number of MPI processes to use
      int _max_procs;
      /// minimum number of elements per MPI process
      Index _min_elems;
      /// the total number of elements
      Index _num_elems;
      /// the desired number of partitions
      Index _num_parts;
      /// our hypergraph structure
      Hypergraph _hypergraph;
      /// the element coloring
      Adjacency::Coloring _coloring;
      /// zoltan internal data
      void* _zz;

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
      explicit PartiZoltan(const Dist::Comm& comm, Index min_elems = 1000u, int max_procs = 1000);

      /// virtual destructor
      virtual ~PartiZoltan();

      /**
       * \brief Executes the Zoltan partitioner
       *
       * \param[in] faces_at_elem
       * A \transient reference to the faces-at-element index set of the mesh that is to be partitioned.
       * The face-dimension of the index set can be chosen freely, but typically in the finite element
       * case this is the vertices-at-element index set of the mesh.
       *
       * \param[in] num_parts
       * The desired number of partitions to create. Must be > 0, but it does not necessarily need
       * to be equal to the number of processes in the communicator.
       *
       * \param[in] weights
       * A vector containing the distribution weights for each element. May be empty if all elements are to be weighted equally.
       *
       * \returns
       * \c true, if the partitioning was successful, or \c false, if some error occurred.
       */
      template<int n_>
      bool execute(const IndexSet<n_>& faces_at_elem, const Index num_parts, const std::vector<Real>& weights)
      {
        Adjacency::Graph graph(Adjacency::RenderType::as_is, faces_at_elem);
        return execute(graph, num_parts, weights);
      }

      /**
      * \brief Executes the Zoltan partitioner
      *
      * \param[in] faces_at_elem
      * A \transient reference to the faces-at-element adjacency graph of the mesh that is to be partitioned.
      * The face-dimension of the graph can be chosen freely, but typically in the finite element
      * case this is the vertices-at-element index set of the mesh.
      *
      * \param[in] num_parts
      * The desired number of partitions to create. Must be > 0, but it does not necessarily need
      * to be equal to the number of processes in the communicator.
      *
      * \param[in] weights
      * A vector containing the distribution weights for each element. May be empty if all elements are to be weighted equally.
      *
      * \returns
      * \c true, if the partitioning was successful, or \c false, if some error occurred.
      */
      bool execute(const Adjacency::Graph& faces_at_elem, const Index num_parts, const std::vector<Real>& weights);

      /**
       * \brief Builds and returns the elements-at-rank graph representing the partitioning.
       */
      Adjacency::Graph build_elems_at_rank() const
      {
        return _coloring.create_partition_graph();
      }

      /**
      * \brief Returns the size of the internal sub-communicator.
      *
      * This only gives useful results after the execute() function has been called.
      */
      int get_sub_comm_size() const
      {
        return _zoltan_comm.size();
      }

    protected:
      /// auxiliary function: builds a hypergraph from the adjacency graph
      void _create_hypergraph(const Adjacency::Graph& faces_at_elem, const std::vector<Real>& weights);
      /// auxiliary function: applies Zoltan onto the hypergraph
      bool _apply_zoltan();
      /// auxiliary function: gathers the partitioned coloring onto the root process
      bool _gather_coloring(const int num_export, const int* export_parts);
      /// auxiliary function: broadcasts the coloring from the root process
      bool _broadcast_coloring(bool zoltan_ok);
    }; // class PartiZoltan
  } // namespace Geometry
} // namespace FEAT

#endif // FEAT_HAVE_ZOLTAN
