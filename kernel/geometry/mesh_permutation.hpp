// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
#pragma once
#ifndef KERNEL_GEOMETRY_MESH_PERMUTATION_HPP
#define KERNEL_GEOMETRY_MESH_PERMUTATION_HPP 1

// includes, FEAT
#include <kernel/util/assertion.hpp>
#include <kernel/adjacency/coloring.hpp>
#include <kernel/adjacency/cuthill_mckee.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/adjacency/permutation.hpp>
#include <kernel/geometry/index_set.hpp>
#include <kernel/geometry/target_set.hpp>
#include <kernel/geometry/vertex_set.hpp>

// includes, system
#include <algorithm>
#include <array>
#include <set>
#include <vector>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      template<typename Shape_, int face_dim_ = Shape_::dimension-1>
      class ISWRenderer
      {
      public:
        static void render(std::array<Adjacency::Graph, Shape_::dimension>& idx,
          const IndexSetWrapper<Shape_, face_dim_>& isw)
        {
          ISWRenderer<Shape_, face_dim_-1>::render(idx, isw);
          idx.at(face_dim_) = Adjacency::Graph(Adjacency::RenderType::as_is,
            isw.template get_index_set<face_dim_>());
        }
      };

      template<typename Shape_>
      class ISWRenderer<Shape_, 0>
      {
      public:
        static void render(std::array<Adjacency::Graph, Shape_::dimension>& idx,
          const IndexSetWrapper<Shape_, 0>& isw)
        {
          idx.at(0) = Adjacency::Graph(Adjacency::RenderType::as_is,
            isw.template get_index_set<0>());
        }
      };

      /// helper class for lexicographic ordering
      template<typename Coord_, int dim_>
      class LexiPoint
      {
      public:
        /// entity (vertex, edge, face, ...) index
        Index idx;
        /// entity barycenter point
        Tiny::Vector<Coord_, dim_> v;

        // tolerance for coordinate equality check
        static constexpr Coord_ tol_ = Coord_(1e-6);

        LexiPoint()
        {
        }

        explicit LexiPoint(Index _idx) :
          idx(_idx),
          v(Coord_(0))
        {
        }

        template<int s_>
        explicit LexiPoint(Index _idx, const Tiny::Vector<Coord_, dim_, s_>& _vtx) :
          idx(_idx),
          v(_vtx)
        {
        }

        void format()
        {
          v.format();
        }

        template<int s_>
        void add(const Tiny::Vector<Coord_, dim_, s_>& vtx)
        {
          v += vtx;
        }

        bool operator<(const LexiPoint& other) const
        {
          // loop over the higher dimension coordinates
          for(int i(dim_-1); i > 0; --i)
          {
            if(v[i] + tol_ < other.v[i])
              return true;
            if(other.v[i] + tol_ <= v[i])
              return false;
          }
          // all coordinates of dimension > 1 are equal, so check the X-coordinate
          return v[0] < other.v[0];
        }
      }; // class LexiPoint

      /// helper class to compute lexicographic permutation
      template<typename Shape_, int shape_dim_ = Shape_::dimension>
      class LexiPermuter
      {
      public:
        template<int nc_, typename Coord_>
        static void compute(
          std::array<Adjacency::Permutation, Shape_::dimension + 1>& perms,
          const IndexSetHolder<Shape_>& ish, const VertexSet<nc_, Coord_>& vtx)
        {
          // recurse down
          LexiPermuter<Shape_, shape_dim_-1>::compute(perms, ish, vtx);

          // get the vertices-at-shape index set
          const auto& idx = ish.template get_index_set<shape_dim_, 0>();
          const Index n = idx.get_num_entities();
          const int nidx = idx.get_num_indices();

          std::vector<LexiPoint<Coord_, nc_>> vset(n);

          // loop over all entities and accumulate all vertices adjacent to that entity
          for(Index i(0); i < n; ++i)
          {
            LexiPoint<Coord_, nc_> v(i);
            for(int j(0); j < nidx; ++j)
              v.add(vtx[idx(i,j)]);
            vset.at(i) = v;
          }

          // sort lexicographically
          std::sort(vset.begin(), vset.end());

          // create permutation
          Adjacency::Permutation p(n);
          Index* vp = p.get_perm_pos();
          for(auto it = vset.begin(); it != vset.end(); ++it, ++vp)
            *vp = it->idx;
          p.calc_swap_from_perm();
          perms.at(shape_dim_) = std::move(p);
        }
      }; // class LexiPermuter

      /// helper class to compute lexicographic permutation
      template<typename Shape_>
      class LexiPermuter<Shape_, 0>
      {
      public:
        template<int nc_, typename Coord_>
        static void compute(
          std::array<Adjacency::Permutation, Shape_::dimension + 1>& perms,
          const IndexSetHolder<Shape_>&, const VertexSet<nc_, Coord_>& vtx)
        {
          const Index n = vtx.get_num_vertices();
          std::vector<LexiPoint<Coord_, nc_>> vset(n);

          // loop over all entities and accumulate all vertices adjacent to that entity
          for(Index i(0); i < n; ++i)
          {
            vset.at(i) = LexiPoint<Coord_, nc_>(i, vtx[i]);
          }

          // sort lexicographically
          std::sort(vset.begin(), vset.end());

          // create permutation
          Adjacency::Permutation p(n);
          Index* vp = p.get_perm_pos();
          for(auto it = vset.begin(); it != vset.end(); ++it, ++vp)
            *vp = it->idx;
          p.calc_swap_from_perm();
          perms.front() = std::move(p);
        }
      }; // class LexiPermuter<...,0>
    } // namespace Intern
    /// \endcond

    /**
     * \brief Mesh permutation strategy enumeration
     */
    enum class PermutationStrategy
    {
      /**
       * \brief no permutation strategy
       *
       * This value indicates that the mesh has not been permuted.
       */
      none = 0,

      /**
       * \brief generic/other permutation strategy
       *
       * This value indicates that the mesh has been permuted by some strategy other
       * than any of the standard permutation strategies.
       *
       * This value can be used by user-defined sorting strategies.
       */
      other,

      /**
       * \brief random permutation strategy
       *
       * This value indicates that all mesh entities have been permuted by using a pseudo-random
       * number generator.
       *
       * This permutation strategy usually represents a worst-case scenario and is therefore
       * primarily only of interest for benchmarking.
       */
      random,

      /**
       * \brief lexicographic permutation strategy
       *
       * This value indicates that all mesh entities have been sorted in lexicographic order, i.e.
       * the entities are sorted by the partial ordering based on their X (1D), Y-X (2D) or
       * Z-Y-X (3D) coordinates of the entity barycenters.
       *
       * \note
       * The lexicographic order in 2D/3D is defined as:
       * - 2D: (x1,y1) < (x2,y2) :<==> (y1 < y2) or ((y1 == y2) and (x1 < x2))
       * - 3D: (x1,y1,z1) < (x2,y2,z2) :<==> (z1 < z2) or ((z1 == z2) and ((y1 < y2) or ((y1 == y2) and (x1 < x2))))
       */
      lexicographic,

      /**
       * \brief colored permutation strategy a.k.a. "red-black" strategy
       *
       * This value indicates that the mesh elements have been sorted based on a coloring of the
       * elements.
       *
       * This permutation strategy will automatically create an element coloring; see the
       * documentation of the MeshPermutation class for more information.
       *
       * \note
       * This permutation strategy is only applied to the mesh elements (cells of highest dimension),
       * i.e. all sub-dimensional entities are left unpermuted.
       */
      colored,

      /**
       * \brief (algebraic) Cuthill-McKee permutation strategy
       *
       * This value indicates that the mesh elements have been sorted based on an algebraic
       * Cuthill-McKee ordering of the elements-neighbors-via-common-vertex adjacency graph.
       *
       * This permutation strategy will automatically create an element layering; see the
       * documentation of the MeshPermutation class for more information.
       *
       * \note
       * This permutation strategy is only applied to the mesh elements (cells of highest dimension),
       * i.e. all sub-dimensional entities are left unpermuted.
       */
      cuthill_mckee,

      /**
       * \brief reversed (algebraic) Cuthill-McKee permutation strategy
       *
       * This value indicates that the mesh elements have been sorted based on a reversed algebraic
       * Cuthill-McKee ordering of the elements-neighbors-via-common-vertex adjacency graph.
       *
       * This permutation strategy will automatically create an element layering; see the
       * documentation of the MeshPermutation class for more information.
       *
       * \note
       * This permutation strategy is only applied to the mesh elements (cells of highest dimension),
       * i.e. all sub-dimensional entities are left unpermuted.
       */
      cuthill_mckee_reversed,

      /**
       * \brief geometric Cuthill-McKee permutation strategy
       *
       * This value indicates that all mesh entities have been sorted based on a geometric
       * Cuthill-McKee ordering.
       *
       * This permutation strategy will automatically create an element layering; see the
       * documentation of the MeshPermutation class for more information.
       */
      geometric_cuthill_mckee,

      /**
       * \brief reversed geometric Cuthill-McKee permutation strategy
       *
       * This value indicates that all mesh entities have been sorted based on a reversed geometric
       * Cuthill-McKee ordering.
       *
       * This permutation strategy will automatically create an element layering; see the
       * documentation of the MeshPermutation class for more information.
       */
      geometric_cuthill_mckee_reversed,
    }; // enum class PermutationStrategy

    /**
     * \brief Mesh permutation class template
     *
     * This class is responsible for computing and managing mesh permutations for conformal meshes.
     * An object of this class is stored as a member variable of the ConformalMesh class.
     *
     * \param[in] Shape_
     * The shape type of the mesh that is to be permuted.
     *
     * This class stores both the forward and inverse permutations for the entities of each
     * dimension. The forward permutation P defines the new positions within the old ordering,
     * i.e. x_new[k] = x_old[P[k]], whereas the corresponding inverse permutation Q defined the
     * old positions within the new array, i.e. x_new[Q[k]] = x_old[k], where x_old is the original
     * unpermuted array and x_new is the permuted array.
     *
     * <u><b>Details about Element Coloring strategies:</b></u>\n
     * <b>The fundamental property (definition) of a valid element coloring strategy is as follows:</b>\n
     * <em>Let E1 != E2 be two elements which share at least one common vertex, and let C(E1) and
     * C(E2) denote the colors of E1 and E2, respectively, then it must hold that C(E1) != C(E2).</em>\n
     * In other words: <em>Two different elements of the same color are not adjacent to each other.</em>
     *
     * If the mesh was permuted by a coloring strategy and the strategy required NC > 0 colors,
     * then (by definition of the coloring strategy) all elements are sorted by their corresponding
     * color, i.e. first come all elements of color 0, then all elements of color 1, etc.
     *
     * This class stores the offsets of the first element for each color in an element coloring
     * vector, whose total length is equal to NC+1, i.e. we have:
     * - _element_coloring[0] = index of first element of color 0 = 0
     * - _element_coloring[1] = index  of first element of color 1
     * - _element_coloring[2] = index  of first element of color 2
     * - ...
     * - _element_coloring[NC-1] =  index  of first element of color NC
     * - _element_coloring[NC] = total number of elements
     *
     * <u><b>Details about Element Layering strategies:</b></u>\n
     * <b>The fundamental property (definition) of a valid element layering strategy is as follows:</b>\n
     * <em>Let E1 != E2 be two elements which share at least one common vertex, and let L(E1) and
     * L(E2) denote the layers of E1 and E2, respectively, then it must hold that |L(E1)- L(E2)| <= 1.</em>\n
     * In other words: <em>An element in layer L can only be adjacent to elements in the (same or
     * adjacent) layers L or (L-1) or (L+1).</em>
     *
     * If the mesh was permuted by one of the layering strategies (e.g. one of the algebraic or
     * geometric Cuthill-McKee permutations) and the strategy produced NL > 0 layers, then (by
     * definition of the layering strategy) all elements are sorted by their corresponding layer,
     * i.e. first come all elements in layer 0, then all elements in layer 1, etc.
     *
     * This class stores the offsets of the first element for each layer in an element layering
     * vector, whose total length is equal to NL+1, i.e. we have:
     * - _element_layering[0] = index of first element of layer 0 = 0
     * - _element_layering[1] = index  of first element of layer 1
     * - _element_layering[2] = index  of first element of layer 2
     * - ...
     * - _element_layering[NL-1] =  index  of first element of layer NL
     * - _element_layering[NL] = total number of elements
     *
     * \author Peter Zajac
     */
    template<typename Shape_>
    class MeshPermutation
    {
    public:
      /// the underlying shape type
      typedef Shape_ ShapeType;

      /// the shape dimension
      static constexpr int shape_dim = ShapeType::dimension;

      /// permutation array typedef
      typedef std::array<Adjacency::Permutation, shape_dim+1> PermArray;

    protected:
      /// the permutation strategy
      PermutationStrategy _strategy;
      /// the actual permutations
      PermArray _perms;
      /// the inverse permutations
      PermArray _inv_perms;
      /// the element color array in case of colored strategy
      /// see the documentation of this class for details
      std::vector<Index> _element_coloring;
      /// the layer array in case of algebraic/geometric Cuthill-McKee strategy
      /// see the documentation of this class for details
      std::vector<Index> _element_layering;

    public:
      /// standard constructor
      MeshPermutation() :
        _strategy(PermutationStrategy::none)
      {
      }

      /**
       * \brief Creates a mesh permutation for a conformal mesh.
       *
       * \param[in] strategy
       * The mesh permutation strategy to be used.
       *
       * \param[in] ish
       * A reference to the index set holder of the conformal mesh.
       *
       * \param[in] vtx
       * A reference to the vertex set of the conformal mesh.
       */
      template<int num_coords_, typename Coord_>
      explicit MeshPermutation(
        PermutationStrategy strategy,
        const IndexSetHolder<Shape_>& ish,
        const VertexSet<num_coords_, Coord_>& vtx) :
        MeshPermutation()
      {
        this->create(strategy, ish, vtx);
      }

      /// move constructor
      MeshPermutation(MeshPermutation&& other) :
        _strategy(other._strategy),
        _perms(std::forward<PermArray>(other._perms)),
        _inv_perms(std::forward<PermArray>(other._inv_perms)),
        _element_coloring(std::forward<std::vector<Index>>(other._element_coloring)),
        _element_layering(std::forward<std::vector<Index>>(other._element_layering))
      {
      }

      /// move-assignment operator
      MeshPermutation& operator=(MeshPermutation&& other)
      {
        if(this == &other)
          return *this;
        _strategy = other._strategy;
        _perms = std::forward<PermArray>(other._perms);
        _inv_perms = std::forward<PermArray>(other._inv_perms);
        _element_coloring = std::forward<std::vector<Index>>(other._element_coloring);
        _element_layering = std::forward<std::vector<Index>>(other._element_layering);
        return *this;
      }

      /// delete copy constructor
      MeshPermutation(const MeshPermutation&) = delete;
      /// delete copy-assign operator
      MeshPermutation& operator=(const MeshPermutation&) = delete;

      /// virtual destructor
      virtual ~MeshPermutation()
      {
      }

      /**
       * \brief Clones another  mesh permutation object into \c this object.
       *
       * \param[in] other
       * A reference to the source object that is to be cloned into \c this object.
       */
      void clone(const MeshPermutation& other)
      {
        this->_strategy = other._strategy;
        for(std::size_t i(0); i <= std::size_t(shape_dim); ++i)
        {
          this->_perms[i] = other._perms[i].clone();
          this->_inv_perms[i] = other._inv_perms[i].clone();
        }
        this->_element_coloring = other._element_coloring;
        this->_element_layering = other._element_layering;
      }

      /// \returns An independent clone of \c this mesh permutation object.
      MeshPermutation clone() const
      {
        MeshPermutation mp;
        mp.clone(*this);
        return mp;
      }

      /// \returns The total size of this object in bytes.
      std::size_t bytes() const
      {
        std::size_t b = sizeof(Index) * (_element_coloring.size() + _element_layering.size());
        for(std::size_t i(0); i <= std::size_t(shape_dim); ++i)
          b +=  sizeof(Index) * 2u * (_perms.at(i).size() + _inv_perms.at(i).size());
        return b;
      }

      /**
       * \brief Checks whether this permutation is empty.
       *
       * A permutation is empty if the permutation strategy is PermutationStrategy::none.
       *
       * \returns \c true, if this permutation is empty, otherwise \c false.
       */
      bool empty() const
      {
        return this->_strategy == PermutationStrategy::none;
      }

      /// \returns The PermutationStrategy that has been used to create this permutation.
      PermutationStrategy get_strategy() const
      {
        return this->_strategy;
      }

      /// \returns A const reference to the internal (forward) permutations array.
      const PermArray& get_perms() const
      {
        return this->_perms;
      }

      /// \returns A const reference to the internal inverse permutations array.
      const PermArray& get_inv_perms() const
      {
        return this->_inv_perms;
      }

      /**
       * \brief Returns a const reference to a (forward) permutation.
       *
       * \param[in] dim
       * The dimension of the entities whose permutation is to be returned. Must be 0 <= \p dim <= shape_dim.
       *
       * \returns
       * A const reference to the (forward) permutation of the entities of dimension \p dim.
       */
      const Adjacency::Permutation& get_perm(int dim = shape_dim) const
      {
        XASSERT((0 <= dim) && (dim <= shape_dim));
        return this->_perms.at(std::size_t(dim));
      }

      /**
       * \brief Returns a const reference to an inverse permutation.
       *
       * \param[in] dim
       * The dimension of the entities whose permutation is to be returned. Must be 0 <= \p dim <= shape_dim.
       *
       * \returns
       * A const reference to the inverse permutation of the entities of dimension \p dim.
       */
      const Adjacency::Permutation& get_inv_perm(int dim = shape_dim) const
      {
        XASSERT((0 <= dim) && (dim <= shape_dim));
        return this->_inv_perms.at(std::size_t(dim));
      }

      /**
       * \brief Returns a const reference to the element coloring vector
       *
       * The returned vector may be empty, which indicates that the permutation is not a colored
       * permutation and therefore no element coloring exists.
       *
       * \returns A const reference to the element coloring vector.
       */
      const std::vector<Index>& get_element_coloring() const
      {
        return this->_element_coloring;
      }

      /**
       * \brief Returns a const reference to the element layering vector
       *
       * The returned vector may be empty, which indicates that the permutation is not a layered
       * permutation and therefore no element layering exists.
       *
       * \returns A const reference to the element layering vector.
       */
      const std::vector<Index>& get_element_layering() const
      {
        return this->_element_layering;
      }

      /**
       * \brief Creates a mesh permutation for a conformal mesh.
       *
       * This function is effectively just a switch-case wrapper for the other create_XXX()
       * member functions.
       *
       * \param[in] strategy
       * The mesh permutation strategy to be used. Must not be PermutationStrategy::other.
       *
       * \param[in] ish
       * A reference to the index set holder of the conformal mesh.
       *
       * \param[in] vtx
       * A reference to the vertex set of the conformal mesh.
       */
      template<int num_coords_, typename Coord_>
      void create(PermutationStrategy strategy, const IndexSetHolder<Shape_>& ish, const VertexSet<num_coords_, Coord_>& vtx)
      {
        switch(strategy)
        {
        case PermutationStrategy::none:
          break;

        case PermutationStrategy::other:
          XABORTM("use the MeshPermutation::create_other() function");
          break;

        case PermutationStrategy::random:
          create_random(ish);
          break;

        case PermutationStrategy::lexicographic:
          create_lexicographic(ish, vtx);
          break;

        case PermutationStrategy::colored:
          create_colored(ish);
          break;

        case PermutationStrategy::cuthill_mckee:
          create_cmk(ish, false);
          break;

        case PermutationStrategy::cuthill_mckee_reversed:
          create_cmk(ish, true);
          break;

        case PermutationStrategy::geometric_cuthill_mckee:
          create_gcmk(ish, vtx, false);
          break;

        case PermutationStrategy::geometric_cuthill_mckee_reversed:
          create_gcmk(ish, vtx, true);
          break;
        }
      }

      /**
       * \brief Initializes a custom permutation.
       *
       * This functions sets the internal strategy enum to PermutationStrategy::other and returns
       * a non-const reference to the internal forwards permutations array, which can then be
       * filled by the caller with the actual permutation data.
       *
       * \attention
       * Do not forget to call the create_inverse_permutations() function once you have filled the
       * forward permutations arrays!
       *
       * \returns
       * A non-const reference to the internal forward permutations array.
       */
      PermArray& create_other()
      {
        XASSERTM(this->_strategy == PermutationStrategy::none, "permutation already created!");
        this->_strategy = PermutationStrategy::other;
        return this->_perms;
      }

      /**
       * \brief Creates random mesh permutation for a conformal mesh.
       *
       * \param[in] ish
       * A reference to the index set holder of the conformal mesh.
       *
       * \note
       * This function creates permutations for the mesh entities for each dimension.
       */
      void create_random(const IndexSetHolder<Shape_>& ish)
      {
        Random rng;
        create_random(ish, rng);
      }

      /**
       * \brief Creates random mesh permutation for a conformal mesh.
       *
       * \param[in] ish
       * A reference to the index set holder of the conformal mesh.
       *
       * \param[in] rng
       * A reference to the pseudo-random number generator to be used.
       *
       * \note
       * This function creates permutations for the mesh entities for each dimension.
       */
      void create_random(const IndexSetHolder<Shape_>& ish, Random& rng)
      {
        XASSERTM(this->_strategy == PermutationStrategy::none, "permutation already created!");
        Index num_entities[shape_dim+1];
        NumEntitiesExtractor<shape_dim>::set_num_entities_with_numverts(ish, num_entities);
        for(std::size_t dim(0); dim <= std::size_t(shape_dim); ++dim)
        {
          _perms.at(dim) = Adjacency::Permutation(num_entities[dim], rng);
          _inv_perms.at(dim) = _perms.at(dim).inverse();
        }
        this->_strategy = PermutationStrategy::random;
      }

      /**
       * \brief Creates a lexicographic mesh permutation for a conformal mesh.
       *
       * \param[in] ish
       * A reference to the index set holder of the conformal mesh.
       *
       * \param[in] vtx
       * A reference to the vertex set of the conformal mesh.
       *
       * \note
       * This function creates permutations for the mesh entities for each dimension.
       */
      template<int num_coords_, typename Coord_>
      void create_lexicographic(const IndexSetHolder<Shape_>& ish, const VertexSet<num_coords_, Coord_>& vtx)
      {
        XASSERTM(this->_strategy == PermutationStrategy::none, "permutation already created!");

        // call the actual helper class
        Intern::LexiPermuter<ShapeType>::compute(this->_perms, ish, vtx);

        // compute inverse permutations
        for(std::size_t dim(0); dim <= std::size_t(shape_dim); ++dim)
        {
          _inv_perms.at(dim) = _perms.at(dim).inverse();
        }

        // save strategy
        this->_strategy = PermutationStrategy::lexicographic;
      }

      /**
       * \brief Creates a colored mesh permutation for a conformal mesh.
       *
       * \param[in] ish
       * A reference to the index set holder of the conformal mesh.
       *
       * \note
       * This function also creates an element coloring automatically.
       *
       * \note
       * This function only creates a permutation for the mesh elements (entities of highest dimension)
       * but not for the lower-dimensional entities.
       */
      void create_colored(const IndexSetHolder<Shape_>& ish)
      {
        XASSERTM(this->_strategy == PermutationStrategy::none, "permutation already created!");

        const auto& verts_at_elem = ish.template get_index_set<shape_dim, 0>();
        Adjacency::Graph elems_at_vert(Adjacency::RenderType::transpose, verts_at_elem);
        Adjacency::Graph elems_at_elem(Adjacency::RenderType::injectify, verts_at_elem, elems_at_vert);

        // create coloring
        Adjacency::Coloring col(elems_at_elem);

        // create coloring partition graph
        Adjacency::Graph cparti = col.create_partition_graph();

        // get the total number of colors
        const Index num_elems = elems_at_elem.get_num_nodes_domain();
        const Index num_colors = cparti.get_num_nodes_domain();

        XASSERT(num_elems == cparti.get_num_indices());

        // create element permutation
        _perms.back() = Adjacency::Permutation(num_elems, Adjacency::Permutation::type_perm, cparti.get_image_idx());
        _inv_perms.back() = _perms.back().inverse();

        // store element coloring
        const Index* dom_ptr = cparti.get_domain_ptr();
        this->_element_coloring.assign(dom_ptr, &dom_ptr[num_colors]);

        // save strategy
        this->_strategy = PermutationStrategy::colored;
      }

      /**
       * \brief Creates a (reversed) algebraic Cuthill-McKee mesh permutation for a conformal mesh.
       *
       * \param[in] ish
       * A reference to the index set holder of the conformal mesh.
       *
       * \param[in] reverse
       * Specifies whether the ordering is to be reversed.
       *
       * \note
       * This function also creates an element layering automatically.
       *
       * \note
       * This function only creates a permutation for the mesh elements (entities of highest dimension)
       * but not for the lower-dimensional entities.
       */
      void create_cmk(const IndexSetHolder<Shape_>& ish, bool reverse)
      {
        XASSERTM(this->_strategy == PermutationStrategy::none, "permutation already created!");

        const auto& verts_at_elem = ish.template get_index_set<shape_dim, 0>();
        Adjacency::Graph elems_at_vert(Adjacency::RenderType::transpose, verts_at_elem);
        Adjacency::Graph elems_at_elem(Adjacency::RenderType::injectify, verts_at_elem, elems_at_vert);

        // create Cuthill-McKee permutation
        _perms.back() = Adjacency::CuthillMcKee::compute(this->_element_layering, elems_at_elem, reverse,
          Adjacency::CuthillMcKee::root_minimum_degree, Adjacency::CuthillMcKee::sort_asc);
        _inv_perms.back() = _perms.back().inverse();

        // save strategy
        this->_strategy = (reverse ? PermutationStrategy::cuthill_mckee_reversed : PermutationStrategy::cuthill_mckee);
      }

      /**
       * \brief Creates a (reversed) geometric Cuthill-McKee mesh permutation for a conformal mesh.
       *
       * \param[in] ish
       * A reference to the index set holder of the conformal mesh.
       *
       * \param[in] vtx
       * A reference to the vertex set of the conformal mesh.
       *
       * \param[in] reverse
       * Specifies whether the ordering is to be reversed.
       *
       * \note
       * This function also creates an element layering automatically.
       *
       * \note
       * This function creates permutations for the mesh entities for each dimension.
       *
       * \remark
       * Currently, the \p vtx parameter is not used by this function, but it is still required
       * for upwards compatibility, because it is possible to improve this algorithm by also
       * using the vertex positioning information.
       */
      template<int num_coords_, typename Coord_>
      void create_gcmk(const IndexSetHolder<Shape_>& ish, const VertexSet<num_coords_, Coord_>& DOXY(vtx), bool reverse)
      {
        // render all <k-faces>-at-element index sets into graph array for convenience
        std::array<Adjacency::Graph, shape_dim> idx_set;
        Intern::ISWRenderer<Shape_>::render(idx_set, ish.template get_index_set_wrapper<shape_dim>());

        // transpose vertices-at-element graph
        Adjacency::Graph& verts_at_elem = idx_set.front();
        Adjacency::Graph elems_at_vert(Adjacency::RenderType::transpose, verts_at_elem);

        // get the number of entities in the mesh
        Index num_entities[shape_dim+1];
        NumEntitiesExtractor<shape_dim>::set_num_entities_with_numverts(ish, num_entities);

        // get number of vertices and elements
        const Index num_verts = num_entities[0];
        const Index num_elems = num_entities[shape_dim];

        // allocate temporary permutation vectors
        std::array<std::vector<Index>, shape_dim+1> perms;
        for(std::size_t i(0); i <= std::size_t(shape_dim); ++i)
          perms.at(i).reserve(num_entities[i]);

        // create mask vectors and initialize them to 0
        // a mask value != 0 indicates that the entity has already been processed
        std::array<std::vector<int>, shape_dim+1> masks;
        for(std::size_t i(0); i <= std::size_t(shape_dim); ++i)
          masks.at(i).resize(num_entities[i]);

        // initialize layer vectors
        std::vector<Index> layer_cur, layer_vert, layer_next;
        layer_cur.reserve(num_elems);
        layer_next.reserve(num_elems);
        layer_vert.reserve(num_verts);

        // initialize element layering
        this->_element_layering.clear();
        this->_element_layering.push_back(Index(0));

        // loop until all elements are processed; this loop performs 1 iteration per domain
        // connectivity region (so usually 1), but it may perform several iterations in case the
        // elements of this domain are not connected -- this may happen e.g. in the MPI parallel
        // case where our local domain is actually just a partition of the whole (connected) domain
        while(Index(perms.back().size()) < num_elems)
        {
          layer_cur.clear();

          // pick a root element of minimal degree (least amount of neighbors)
          // this chosen root element will usually be a boundary or even a corner element
          Index idx = _pick_gcmk_root(idx_set.front(), elems_at_vert, masks.back());
          layer_cur.push_back(idx);
          masks.back()[idx] = 1;

          // loop through the current element layer
          while(!layer_cur.empty())
          {
            // loop over all elements in current layer
            layer_next.clear();
            for(auto it = layer_cur.begin(); it != layer_cur.end(); ++it)
            {
              // get and add current element
              const Index cur_elem = *it;
              perms.back().push_back(cur_elem);

              // add all vertices adjacent to the current element to the next vertex layer
              layer_vert.clear();
              for(auto kt = verts_at_elem.image_begin(cur_elem); kt != verts_at_elem.image_end(cur_elem); ++kt)
              {
                const Index cur_vert = *kt;
                // already processed?
                if(masks.front()[cur_vert] != 0)
                  continue;
                // add layer vertex
                layer_vert.push_back(cur_vert);
                perms.front().push_back(cur_vert);
                masks.front()[cur_vert] = 1;
              }

              // add all k-faces (edges, faces, ...) adjacent to the current element to the next k-face layer
              for(std::size_t k(1); k < std::size_t(shape_dim); ++k)
              {
                for(auto kt = idx_set.at(k).image_begin(cur_elem); kt != idx_set.at(k).image_end(cur_elem); ++kt)
                {
                  const Index cur_face = *kt;
                  if(masks.at(k)[cur_face] != 0)
                    continue;
                  perms.at(k).push_back(cur_face);
                  masks.at(k)[cur_face] = 1;
                }
              }

              // add all elements adjacent via a vertex to the next element layer
              for(auto jt = layer_vert.begin(); jt != layer_vert.end(); ++jt)
              {
                const Index cur_vert = *jt;
                for(auto kt = elems_at_vert.image_begin(cur_vert); kt != elems_at_vert.image_end(cur_vert); ++kt)
                {
                  const Index next_elem = *kt;
                  if(masks.back()[next_elem] != 0)
                    continue;
                  layer_next.push_back(next_elem);
                  masks.back()[next_elem] = 1;
                }
              }
            }

            // continue with next element layer
            layer_cur.swap(layer_next);

            // store element layer size
            this->_element_layering.push_back(Index(perms.back().size()));

          } // while(!layer_cur.empty())
        } // while(perms.back().size() < num_elems)

        // revert if desired
        if(reverse)
        {
          // reverse permutations
          for(std::size_t dim(0); dim <= std::size_t(shape_dim); ++dim)
          {
            for(std::size_t i(0), j(perms[dim].size()-1); i < j; ++i, --j)
              std::swap(perms[dim][i], perms[dim][j]);
          }

          // reverse layers
          const Index lay_back = Index(this->_element_layering.size()) - 1u;
          for(Index k(lay_back); k > 0u; --k)
          {
            // compute layer sizes from offsets
            this->_element_layering.at(k) -= this->_element_layering.at(k-1u);
          }
          for(Index k(1u), l(lay_back); k < l; ++k, --l)
          {
            // reverse layer sizes
            std::swap(this->_element_layering.at(k), this->_element_layering.at(l));
          }
          for(Index k(1u); k <= lay_back; ++k)
          {
            // compute layering offsets from sizes
            this->_element_layering.at(k) += this->_element_layering.at(k-1u);
          }
        }

        // create the actual permutation objects
        for(std::size_t dim(0); dim <= std::size_t(shape_dim); ++dim)
        {
          // sanity check: all entities should have been processed
          XASSERT(Index(perms.at(dim).size()) == num_entities[dim]);

          // create actual permutation and its inverse
          _perms.at(dim) = Adjacency::Permutation(Index(perms.at(dim).size()),
            Adjacency::Permutation::type_perm, perms.at(dim).data());
          _inv_perms.at(dim) = _perms.at(dim).inverse();
        }

        // save strategy
        this->_strategy = (reverse ? PermutationStrategy::geometric_cuthill_mckee_reversed : PermutationStrategy::geometric_cuthill_mckee);
      }

      /**
       * \brief (Re)Creates the inverse permutations from the forward ones.
       */
      void create_inverse_permutations()
      {
        for(std::size_t dim(0); dim <= std::size_t(shape_dim); ++dim)
        {
          if(_perms.at(dim).empty())
            _inv_perms.at(dim) = Adjacency::Permutation();
          else
            _inv_perms.at(dim) = _perms.at(dim).inverse();
        }
      }

      /**
       * \brief Validates the permutation sizes.
       *
       * \param[in] num_entities
       * An array of size shape_dim+1 that contains the number of entities for each dimension.
       *
       * \returns
       * - 0, if all permutation sizes are valid, or
       * - j+1, if the permutations of dimension k have the wrong size
       */
      int validate_sizes(const Index* num_entities) const
      {
        XASSERT(num_entities != nullptr);

        // check permutations
        for(std::size_t i(0); i <= std::size_t(shape_dim); ++i)
        {
          const Adjacency::Permutation& perm = this->_perms.at(i);
          const Adjacency::Permutation& iperm = this->_inv_perms.at(i);

          if(!perm.empty() && (perm.size() != num_entities[i]))
            return 2*int(i)+1;
          if(!iperm.empty() && (iperm.size() != num_entities[i]))
            return 2*int(i)+2;
        }

        // ok
        return 0;
      }

    protected:
      /**
       * \brief Auxiliary helper function: pick an unmasked root element for GCMK
       */
      static Index _pick_gcmk_root(const Adjacency::Graph& v_at_e, const Adjacency::Graph& e_at_v,
        const std::vector<int>& mask)
      {
        std::set<Index> iset;

        // pick a root node of minimal degree
        Index deg = v_at_e.get_num_nodes_domain() + e_at_v.get_num_nodes_domain();
        Index idx(0);
        for(Index i(0); i < v_at_e.get_num_nodes_domain(); ++i)
        {
          // skip masked elements
          if(mask[i] != 0)
            continue;

          // compute degree
          iset.clear();
          for(auto it = v_at_e.image_begin(i); it != v_at_e.image_end(i); ++it)
            for(auto jt = e_at_v.image_begin(*it); jt != e_at_v.image_end(*it); ++jt)
              iset.insert(*jt);

          // is this the new minimum?
          Index d = Index(iset.size());
          if(d < deg)
          {
            idx = i;
            deg = d;
          }
        }

        return idx;
      }
    }; // class MeshPermutation<ConformalMesh<...>>
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_MESH_PERMUTATION_HPP
