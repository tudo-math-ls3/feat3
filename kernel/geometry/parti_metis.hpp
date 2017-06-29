#pragma once
#ifndef KERNEL_GEOMETRY_PARTI_METIS_HPP
#define KERNEL_GEOMETRY_PARTI_METIS_HPP 1

#include <kernel/base_header.hpp>

#ifdef FEAT_HAVE_PARMETIS
#include <kernel/archs.hpp>
FEAT_DISABLE_WARNINGS
#include <parmetis.h>
FEAT_RESTORE_WARNINGS

#include <kernel/adjacency/graph.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/adjacency/colouring.hpp>

#include <map>
#include <vector>
#include <list>

namespace FEAT
{
  namespace Geometry
  {
    /** \brief Metis based Partitioner class template declaration */
    template<typename Mesh_>
    class PartiMetis;

    /**
     * \brief Metis based Partitioner class template specialisation for ConformalMesh
     *
     * The basic usage of this class is as follows:
     * -# Refine the mesh until it contains at least num_rank cells.
     * -# Create an object of this class and pass the to-be-partitioned mesh as well
     *    as the desired number of ranks/patches to the constructor.
     * -# Create the Elements-At-Rank graph using the #build_elems_at_rank() function.
     *
     * \author Dirk Ribbrock
     */
    template<typename Shape_, int num_coords_, int stride_, typename Coord_>
    class PartiMetis<ConformalMesh<Shape_, num_coords_, stride_, Coord_>>
    {
      private:
        /// number of elements in input mesh
        const Index _num_elems;
        /// number of desired ranks/patches
        const Index _num_ranks;
        /// the achieved optimisation value
        idx_t * _objval;
        /// node to partition mapping
        idx_t * _part;


      public:
      /// our mesh type
      typedef ConformalMesh<Shape_, num_coords_, stride_, Coord_> MeshType;

      /**
       * \brief Constructor
       *
       * \param[in] mesh
       * The mesh that is to be partitioned (on some refined level).
       *
       * \param[in] num_ranks
       * The number of ranks involved in the mesh
       */
      explicit PartiMetis(Adjacency::Graph& adj_graph, const Index num_ranks) :
        _num_elems(adj_graph.get_num_nodes_domain()),
        _num_ranks(num_ranks),
        _objval(nullptr),
        _part(nullptr)
      {
        idx_t* nvtxs = new idx_t[1];
        nvtxs[0] = (idx_t)_num_elems;
        idx_t* ncon = new idx_t[1];
        ncon[0] = 1;
        idx_t* nparts = new idx_t[1];
        nparts[0] = (idx_t)_num_ranks;

        idx_t* xadj = new idx_t[_num_elems + 1];
        idx_t* adjncy = new idx_t[adj_graph.get_num_indices() - _num_elems];

        // convert and remove self-adjacencies
        xadj[0] = 0;
        for(Index node(0); node < _num_elems ; ++node)
        {
          idx_t j = xadj[node];
          for(auto it = adj_graph.image_begin(node); it != adj_graph.image_end(node); ++it)
          {
            if(*it != node)
            {
              adjncy[std::size_t(j)] = idx_t(*it);
              ++j;
            }
          }
          xadj[node+1] = j;
        }

        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);
        options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM; //avoid randomisation effects
        options[METIS_OPTION_NUMBERING] = 0;
        options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;

        _objval = new idx_t[1];
        _part = new idx_t[_num_elems];
        //int result = METIS_PartGraphKway(nvtxs, ncon, xadj, adjncy, NULL, NULL, NULL, nparts, NULL, NULL, options, _objval, _part);
        int result = METIS_PartGraphRecursive(nvtxs, ncon, xadj, adjncy, NULL, NULL, NULL, nparts, NULL, NULL, options, _objval, _part);
        XASSERT(result == METIS_OK);

        delete[] nvtxs;
        delete[] ncon;
        delete[] nparts;
        delete[] xadj;
        delete[] adjncy;
      }

      /// Destructor
      ~PartiMetis()
      {
        delete[] _objval;
        delete[] _part;
      }

      /**
       * \brief Returns the Elements-at-Rank graph of the partitioning.
       *
       * \returns
       * The Elements-at-Rank graph of the partitioning.
       */
      Adjacency::Graph build_elems_at_rank() const
      {
        Adjacency::Graph graph(_num_ranks, _num_elems, _num_elems);
        Index* ptr = graph.get_domain_ptr();
        Index* idx = graph.get_image_idx();

        std::vector<std::list<Index>> temp;
        temp.resize(_num_ranks);
        for (Index i(0) ; i < _num_elems ; ++i)
        {
          temp.at(Index(_part[i])).push_back(i);
        }

        Index offset(0);
        ptr[0] = 0;
        for (Index i(0) ; i < _num_ranks ; ++i)
        {
          Index j(0);
          for (auto item : temp.at(i))
          {
            idx[offset + j] = item;
            ++j;
          }
          offset += Index(temp.at(i).size());
          ptr[i+1] = offset;
        }

        return graph;
      }
    }; // class PartiMetis
  } // namespace Geometry
} // namespace FEAT

#endif // FEAT_HAVE_PARMETIS
#endif // KERNEL_GEOMETRY_PARTI_METIS_HPP
