// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_PARTI_2LVL_HPP
#define KERNEL_GEOMETRY_PARTI_2LVL_HPP 1

#include <kernel/adjacency/graph.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/intern/standard_refinement_traits.hpp>
#include <kernel/util/assertion.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      template<typename Shape_>
      struct Parti2LvlHelper
      {
        static constexpr Index factor = StandardRefinementTraits<Shape_, Shape_::dimension>::count;
        static constexpr Index lvlinc = Index(1);
      };

      template<int dim_>
      struct Parti2LvlHelper<Shape::Hypercube<dim_>>
      {
        static constexpr Index factor = Index(2);
        static constexpr Index lvlinc = Index(dim_);
      };
    } // namespace Intern
    /// \endcond

    /** \brief 2-Level-Partitioner class template declaration */
    template<typename Mesh_>
    class Parti2Lvl;

    /**
     * \brief 2-Level-Partitioner class template specialization for ConformalMesh
     *
     * This class implements a <em>2-Level-Partitioner</em>, i.e. a simple partitioner
     * that exploits the 2-level refinement algorithm and numbering.
     *
     * The 2-Level-Partitioner can provide a valid partitioning of a (refined)
     * mesh, if the number of required patches matches the number of elements
     * in a refined mesh.
     *
     * For Hypercube meshes, this algorithm can do even better: Let \em n be the
     * number of elements in the mesh and \em p the desired number of patches, then
     * this algorithm can provide a valid partitioning if \f$p = n\cdot2^k\f$ for
     * some \f$k\in\mathbb{N}\f$.
     *
     * The basic usage of this class is as follows:
     * -# Create an object of this class and pass the to-be-partitioned mesh as well
     *    as the desired number of ranks/patches to the constructor.
     * -# If the success() function returns \c true, then a valid 2-level partitioning
     *    has been found. If the function returns \c false, then no 2-level partitioning
     *    exists and one has to try some other partitioning strategy.
     * -# Refine the mesh up to the required refinement level returned by the #parti_level()
     *    function.
     * -# Create the Elements-At-Rank graph using the #build_elems_at_rank() function.
     *
     * \author Peter Zajac
     */
    template<typename Shape_, int num_coords_, typename Coord_>
    class Parti2Lvl<ConformalMesh<Shape_, num_coords_, Coord_>>
    {
    public:
      /// our mesh type
      typedef ConformalMesh<Shape_, num_coords_, Coord_> MeshType;

    protected:
      /// number of elements in input mesh
      const Index _num_elems;
      /// number of desired ranks/patches
      const Index _num_ranks;
      /// specifies whether partitioning is possible
      bool _success;
      /// refinement level of partitioned mesh;
      Index _ref_lvl;
      /// number of elements in refined mesh
      Index _ref_elems;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] mesh
       * The mesh that is to be partitioned (on some refined level)
       *
       * \param[in] num_ranks
       * The desired number of ranks.
       */
      explicit Parti2Lvl(const MeshType& mesh, Index num_ranks) :
        _num_elems(mesh.get_num_elements()),
        _num_ranks(num_ranks),
        _success(false),
        _ref_lvl(0),
        _ref_elems(0)
      {
        XASSERT(num_ranks > Index(0));

        const Index factor = Intern::Parti2LvlHelper<Shape_>::factor;
        const Index lvlinc = Intern::Parti2LvlHelper<Shape_>::lvlinc;
        Index count(_num_elems);
        Index power(0);

        while(count < _num_ranks)
        {
          count *= factor;
          ++power;
        }

        if(count != _num_ranks)
        {
          // no 2-level partitioning possible
          return;
        }

        _success = true;

        // compute refinement level
        _ref_lvl = (power + lvlinc - Index(1)) / lvlinc;

        // compute refined element count
        const Index ref_fac(Intern::StandardRefinementTraits<Shape_, Shape_::dimension>::count);
        _ref_elems = _num_elems;
        for(Index i(0); i < _ref_lvl; ++i)
          _ref_elems *= ref_fac;
      }

      /**
       * \brief Specifies whether a 2-level partitioning exists.
       *
       * \returns
       * \c true, if a 2-level partitioning of the mesh passed to the constructor exists, otherwise \c false.
       */
      bool success() const
      {
        return _success;
      }

      /**
       * \brief Returns the required refinement level for the 2-level partitioning.
       */
      Index parti_level() const
      {
        return _ref_lvl;
      }

      /**
       * \brief Returns the Elements-at-Rank graph of the partitioning.
       *
       * \returns
       * The Elements-at-Rank graph of the partitioning.
       */
      Adjacency::Graph build_elems_at_rank() const
      {
        Adjacency::Graph graph(_num_ranks, _ref_elems, _ref_elems);
        Index* ptr = graph.get_domain_ptr();
        Index* idx = graph.get_image_idx();

        // build pointer array
        const Index elems_per_rank = _ref_elems / _num_ranks;
        for(Index i(0); i <= _num_ranks; ++i)
          ptr[i] = elems_per_rank * i;

        // build index array
        for(Index i(0); i < _ref_elems; ++i)
          idx[i] = i;

        return graph;
      }
    }; // class Parti2Lvl
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_PARTI_2LVL_HPP
