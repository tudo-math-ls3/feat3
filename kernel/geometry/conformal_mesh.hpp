// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_CONFORMAL_MESH_HPP
#define KERNEL_GEOMETRY_CONFORMAL_MESH_HPP 1

// includes, FEAT
#include <kernel/geometry/factory.hpp>
#include <kernel/geometry/mesh_permutation.hpp>
#include <kernel/geometry/intern/facet_neighbors.hpp>
#include <kernel/geometry/intern/standard_index_refiner.hpp>
#include <kernel/geometry/intern/standard_vertex_refiner.hpp>
#include <kernel/geometry/index_calculator.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /**
     * \brief Conformal mesh class template
     *
     * \todo detailed documentation
     *
     * \tparam Shape_
     * The shape that is to be used for the mesh. Must be either Shape::Simplex<n> or Shape::Hypercube<n>
     * for some \c n > 0.
     *
     * For more details on meshes, see the related doxygen page \ref mesh_file_format.
     *
     * \author Peter Zajac
     */
    template<
      typename Shape_,
      int num_coords_ = Shape_::dimension,
      typename Coord_ = Real>
    class ConformalMesh
    {
      static_assert(num_coords_ >= Shape_::dimension, "invalid number of coordinates");

    public:
      /// Shape type
      typedef Shape_ ShapeType;

      /// Coordinate type
      typedef Coord_ CoordType;

      /// Vertex set type
      typedef VertexSet<num_coords_, Coord_> VertexSetType;

      /// Vertex type
      typedef typename VertexSetType::VertexType VertexType;

      /// index set holder type
      typedef IndexSetHolder<ShapeType> IndexSetHolderType;

      /// shape dimension
      static constexpr int shape_dim = ShapeType::dimension;
      /// world dimension
      static constexpr int world_dim = VertexSetType::num_coords;

      /// this mesh is not structured
      static constexpr bool is_structured = false;

      /// mesh permutation type
      typedef MeshPermutation<ShapeType> MeshPermutationType;

      /**
       * \brief Index set type class template
       *
       * This nested class template is used to define the return type of the ConformalMesh::get_index_set()
       * function template.
       *
       * \tparam cell_dim_, face_dim_
       * The cell and face dimension parameters as passed to the ConformalMesh::get_index_set() function template.
       */
      template<
        int cell_dim_,
        int face_dim_>
      struct IndexSet
      {
        static_assert(cell_dim_ <= shape_dim, "invalid cell dimension");
        static_assert(face_dim_ < cell_dim_, "invalid face/cell dimension");
        static_assert(face_dim_ >= 0, "invalid face dimension");

        /// index set type
        typedef FEAT::Geometry::IndexSet
        <
            Shape::FaceTraits
            <
              typename Shape::FaceTraits< ShapeType, cell_dim_> ::ShapeType,
              face_dim_
            >
            ::count
          > Type;
      }; // struct IndexSet<...>

      /// index set type for storing neighbor adjacency information
      typedef typename IndexSet<shape_dim, shape_dim-1>::Type NeighborSetType;

    protected:
      /// number of entities for each shape dimension
      Index _num_entities[shape_dim + 1];

      /// the vertex set of the mesh
      VertexSetType _vertex_set;

      /// the index sets of the mesh
      IndexSetHolderType _index_set_holder;

      /// Information about cells sharing a facet
      NeighborSetType _neighbors;

      /// mesh permutation (if permuted)
      MeshPermutationType _permutation;

    public:
      /**
       * \brief Constructor.
       *
       * \param[in] num_entities
       * An array of length at least #shape_dim + 1 holding the number of entities for each shape dimension.
       * Must not be \c nullptr.
       */
      explicit ConformalMesh(const Index num_entities[]) :
        _vertex_set(num_entities[0]),
        _index_set_holder(num_entities),
        _neighbors(num_entities[shape_dim]),
        _permutation()
      {
        for(int i(0); i <= shape_dim; ++i)
        {
          //XASSERTM(num_entities[i] > 0, "Number of entities must not be zero!");
          _num_entities[i] = num_entities[i];
        }
      }

      /**
       * \brief Factory constructor
       *
       * \param[in] factory
       * The factory that is to be used to create the mesh.
       */
      explicit ConformalMesh(Factory<ConformalMesh>& factory) :
        _vertex_set(factory.get_num_entities(0)),
        _index_set_holder(Intern::NumEntitiesWrapper<shape_dim>(factory).num_entities),
        _neighbors(Intern::NumEntitiesWrapper<shape_dim>(factory).num_entities[shape_dim]),
        _permutation()
      {
        // Compute entity counts
        Intern::NumEntitiesWrapper<shape_dim>::apply(factory, _num_entities);

        // fill vertex set
        factory.fill_vertex_set(_vertex_set);

        // fill index sets
        factory.fill_index_sets(_index_set_holder);

        // Recompute entity counts. This is mainly for the MeshStreamerFactory, as complete information about the
        // number of edges/faces might not be known until after fill_index_sets is called
        Intern::NumEntitiesWrapper<shape_dim>::apply(factory, _num_entities);

        // Fill neighbor information. This needs facet at cell information, so it needs to be called after
        // fill_index_sets() etc.
        fill_neighbors();
      }

      /// move constructor
      ConformalMesh(ConformalMesh&& other) :
        _vertex_set(std::forward<VertexSetType>(other._vertex_set)),
        _index_set_holder(std::forward<IndexSetHolderType>(other._index_set_holder)),
        _neighbors(std::forward<NeighborSetType>(other._neighbors)),
        _permutation(std::forward<MeshPermutationType>(other._permutation))
      {
        for(int i(0); i <= shape_dim; ++i)
        {
          _num_entities[i] = other.get_num_entities(i);
        }
      }

      /// move-assignment operator
      ConformalMesh& operator=(ConformalMesh&& other)
      {
        // avoid self-move
        if(this == &other)
          return *this;

        _vertex_set = std::forward<VertexSetType>(other._vertex_set);
        _index_set_holder = std::forward<IndexSetHolderType>(other._index_set_holder);
        _neighbors = std::forward<NeighborSetType>(other._neighbors);
        _permutation = std::forward<MeshPermutationType>(other._permutation);

        for(int i(0); i <= shape_dim; ++i)
        {
          _num_entities[i] = other.get_num_entities(i);
        }

        return *this;
      }

      /// delete copy constructor
      ConformalMesh(const ConformalMesh&) = delete;

      /// delete copy-assign operator
      ConformalMesh& operator=(const ConformalMesh&) = delete;

      /// virtual destructor
      virtual ~ConformalMesh()
      {
      }

      /**
       * \brief Clones another conformal mesh object into \c this object.
       *
       * \param[in] other
       * A reference to the source object that is to be cloned into \c this object.
       */
      void clone(const ConformalMesh& other)
      {
        for(int i(0); i <= shape_dim; ++i)
          this->_num_entities[i] = other._num_entities[i];
        this->_vertex_set = other._vertex_set.clone();
        this->_index_set_holder.clone(other._index_set_holder);
        this->_neighbors = other._neighbors.clone();
        this->_permutation.clone(other._permutation);
      }

      /// \returns An independent clone of  \c this mesh object.
      ConformalMesh clone() const
      {
        ConformalMesh mesh(this->_num_entities);
        mesh._vertex_set = this->_vertex_set.clone();
        mesh._index_set_holder.clone(this->_index_set_holder);
        mesh._neighbors = this->_neighbors.clone();
        mesh._permutation.clone(this->_permutation);
        return mesh;
      }

      /// \returns The size of dynamically allocated memory in bytes.
      std::size_t bytes() const
      {
        return _vertex_set.bytes() + _index_set_holder.bytes() + _neighbors.bytes() + _permutation.bytes();
      }

      /**
       * \brief Returns the number of entities.
       *
       * \param[in] dim
       * The dimension of the entity whose count is to be returned. Must be 0 <= \p dim <= #shape_dim.
       *
       * \returns
       * The number of entities of dimension \p dim.
       */
      Index get_num_entities(int dim) const
      {
        XASSERT(dim >= 0);
        XASSERT(dim <= shape_dim);
        return _num_entities[dim];
      }

      /**
       * \brief Returns the number of vertices.
       * \returns The number of vertices.
       */
      Index get_num_vertices() const
      {
        return _num_entities[0];
      }

      /**
       * \brief Returns the number of elements.
       * \returns The number of elements.
       */
      Index get_num_elements() const
      {
        return _num_entities[shape_dim];
      }

      /**
       * \brief Checks whether the mesh is permuted.
       *
       * \returns \c true if the mesh is permuted, otherwise \c false.
       */
      bool is_permuted() const
      {
        return !this->_permutation.empty();
      }

      /**
       * \brief Returns a reference to the underlying mesh permutation object.
       */
      const MeshPermutationType& get_mesh_permutation() const
      {
        return this->_permutation;
      }

      /**
       * \brief Creates a mesh permutation based on one of the standard permutation strategies.
       *
       * This function creates a new mesh permutation and also applies that permutation to the
       * vertex set and all the index sets stored in this mesh object.
       *
       * \param[in] strategy
       * The permutation strategy to use, see #MeshPermutation for more details.
       *
       * \note
       * If you want to use a custom permutation other than one of the standard permutation
       * strategies then use set_permutation() instead.
       *
       * \attention
       * A mesh can only be permuted once and therefore this function will fire an assertion if
       * the mesh is already permuted.
       */
      void create_permutation(PermutationStrategy strategy)
      {
        // make sure that we don't already have a permutation
        XASSERTM(this->_permutation.empty(), "mesh is already permuted!");

        // create the permutation
        this->_permutation.create(strategy, this->_index_set_holder, this->_vertex_set);

        // permute vertex set
        this->_vertex_set.permute(this->_permutation.get_perm(0));

        // permute index sets
        this->_index_set_holder.permute(this->_permutation.get_perms(), this->_permutation.get_inv_perms());
      }

      /**
       * \brief Sets a custom mesh permutation for this mesh.
       *
       * This function can be used to apply a mesh permutation that is created using some other
       * approach than the predefined standard permutation strategies.
       *
       * This function also applies that permutation to the  vertex set and all the index sets
       * stored in this mesh object.
       *
       * \param[in] mesh_perm
       * The mesh permutation to use.
       *
       * \attention
       * A mesh can only be permuted once and therefore this function will fire an assertion if
       * the mesh is already permuted.
       */
      void set_permutation(MeshPermutationType&& mesh_perm)
      {
        // make sure that we don't already have a permutation
        XASSERTM(this->_permutation.empty(), "mesh is already permuted!");

        // check the dimensions
        XASSERTM(mesh_perm.validate_sizes(this->_num_entities) == 0, "mesh permutation has invalid size!");

        // save the permutation
        this->_permutation = std::forward<MeshPermutationType>(mesh_perm);

        // permute vertex set
        this->_vertex_set.permute(this->_permutation.get_perm(0));

        // permute index sets
        this->_index_set_holder.permute(this->_permutation.get_perms(), this->_permutation.get_inv_perms());
      }

      /**
       * \brief Validates the element coloring.
       *
       * An element coloring is valid, if any pair of two different elements, which share at least
       * one common vertex, have different colors.
       *
       * \returns \c true, if the element coloring is either valid or empty, otherwise \c false.
       */
      bool validate_element_coloring() const
      {
        // no coloring?
        const std::vector<Index>& coloring = this->get_mesh_permutation().get_element_coloring();
        if(coloring.empty())
          return true;

        // get vertices-at-element index set
        const auto& verts_at_elem = this->template get_index_set<shape_dim, 0>();

        // render transpose
        Adjacency::Graph elems_at_vert(Adjacency::RenderType::transpose, verts_at_elem);

        // loop over all color blocks
        for(std::size_t icol(0); icol+1u < coloring.size(); ++icol)
        {
          // get the bounds of our current color block
          const Index iel_beg = coloring[icol];
          const Index iel_end = coloring[icol+1u];

          // loop over all elements in the current color block
          for(Index iel(iel_beg); iel < iel_end; ++iel)
          {
            // loop over all vertices adjacent to this element
            for(int ivt(0); ivt < verts_at_elem.num_indices; ++ivt)
            {
              // loop over all elements adjacent to this vertex
              const Index ivtx = verts_at_elem(iel, ivt);
              for(auto it = elems_at_vert.image_begin(ivtx); it != elems_at_vert.image_end(ivtx); ++it)
              {
                // two adjacent element must not be in the same color block
                if((iel_beg <= *it) && (*it < iel_end) && (*it != iel))
                  return false; // invalid coloring
              }
            }
          }
        } // next color block

        // ok, coloring is valid
        return true;
      }

      /**
       * \brief Validates the element layering.
       *
       * An element layering is valid, if any pair of two different elements, which share at least
       * one common vertex, have different colors.
       *
       * \returns \c true, if the element layering is either valid or empty, otherwise \c false.
       */
      bool validate_element_layering() const
      {
        // no layering?
        const std::vector<Index>& layering = this->get_mesh_permutation().get_element_layering();
        if(layering.empty())
          return true;

        // get vertices-at-element index set
        const auto& verts_at_elem = this->template get_index_set<shape_dim, 0>();

        // render transpose
        Adjacency::Graph elems_at_vert(Adjacency::RenderType::transpose, verts_at_elem);

        // loop over all layers
        for(std::size_t ilay(0); ilay+1u < layering.size(); ++ilay)
        {
          // get the bounds of our current layer
          const Index iel_beg = layering[ilay];
          const Index iel_end = layering[ilay+1u];

          // get the lower bound for valid neighbors of our current layer = beginning of previous layer
          const Index iel_lower = layering[Math::max(ilay, std::size_t(1)) - 1u];

          // get the upper bound for valid neighbors of our current layer = end of next layer
          const Index iel_upper = layering[Math::min(ilay+2u, layering.size()-1u)];

          // loop over all elements in the current layer
          for(Index iel(iel_beg); iel < iel_end; ++iel)
          {
            // loop over all vertices adjacent to this element
            for(int ivt(0); ivt < verts_at_elem.num_indices; ++ivt)
            {
              // loop over all elements adjacent to this vertex
              const Index ivtx = verts_at_elem(iel, ivt);
              for(auto it = elems_at_vert.image_begin(ivtx); it != elems_at_vert.image_end(ivtx); ++it)
              {
                // adjacent element outside of adjacent layers?
                if(!((iel_lower <= *it) && (*it < iel_upper)))
                  return false; // invalid layer
              }
            }
          }
        } // next color block

        // ok, layering is valid
        return true;
      }

      /// Fills the neighbor index set
      void fill_neighbors()
      {
        // Facet at cell index set
        auto& facet_idx = get_index_set<shape_dim, shape_dim -1>();

        XASSERTM(get_num_entities(shape_dim-1) == facet_idx.get_index_bound(), "mesh num_entities / index_set num_entities mismatch");

        if(_neighbors.get_num_indices() == Index(0))
          _neighbors = std::move(typename IndexSet<shape_dim, shape_dim-1>::Type(get_num_entities(shape_dim)));

        Intern::FacetNeighbors::compute(_neighbors, facet_idx);
      }

      /// \returns A reference to the facet neighbor relations
      typename IndexSet<shape_dim, shape_dim-1>::Type& get_neighbors()
      {
        return _neighbors;
      }

      /// \copydoc get_neighbors()
      const typename IndexSet<shape_dim, shape_dim-1>::Type& get_neighbors() const
      {
        return _neighbors;
      }

      /// Returns a reference to the vertex set of the mesh.
      VertexSetType& get_vertex_set()
      {
        return _vertex_set;
      }

      /** \copydoc get_vertex_set() */
      const VertexSetType& get_vertex_set() const
      {
        return _vertex_set;
      }

      /**
       * \brief Returns the reference to an index set.
       *
       * \tparam cell_dim_
       * The dimension of the entity whose index set is to be returned.
       *
       * \tparam face_dim_
       * The dimension of the face that the index set refers to.
       *
       * \returns
       * A reference to the index set.
       */
      template<
        int cell_dim_,
        int face_dim_>
      typename IndexSet<cell_dim_, face_dim_>::Type& get_index_set()
      {
        return _index_set_holder.template get_index_set_wrapper<cell_dim_>().template get_index_set<face_dim_>();
      }

      /** \copydoc get_index_set() */
      template<
        int cell_dim_,
        int face_dim_>
      const typename IndexSet<cell_dim_, face_dim_>::Type& get_index_set() const
      {
        return _index_set_holder.template get_index_set_wrapper<cell_dim_>().template get_index_set<face_dim_>();
      }

      /// \cond internal
      IndexSetHolderType& get_index_set_holder()
      {
        return _index_set_holder;
      }

      const IndexSetHolderType& get_index_set_holder() const
      {
        return _index_set_holder;
      }

      IndexSetHolderType* get_topology()
      {
        return &_index_set_holder;
      }

      const IndexSetHolderType* get_topology() const
      {
        return &_index_set_holder;
      }
      /// \endcond

      /**
       * \brief Deducts the topology from the Vertex-At-Shape index set.
       */
      void deduct_topology_from_top()
      {
        RedundantIndexSetBuilder<ShapeType>::compute(_index_set_holder);
        NumEntitiesExtractor<shape_dim>::set_num_entities(_index_set_holder, _num_entities);
        this->fill_neighbors();
      }

      /**
       * \brief Applies a "proper rigid" transformation onto the mesh.
       *
       * Let \e v denote the \p origin world point, \e w the \p offset world point and \e R
       * the rotation matrix corresponding to the \p angles, then this function applies the
       * following transformation for any vertex \e x of the vertex set:
       *
       *   \f[ x \mapsto w + R\cdot (x - v) \f]
       *
       * \param[in] origin
       * The origin of the transformation. This is subtracted from any vertex before applying the
       * rotation.
       *
       * \param[in] angles
       * The angles of the rotation matrix.
       * - 2D: the rotation angle in radians is stored as:
       *   - angles(0): rotation angle
       *   - angles(1): \e ignored
       * - 3D: the rotation angles in radians stored as:
       *   - angles(0): yaw angle
       *   - angles(1): pitch angle
       *   - angles(2): roll angle
       *
       * \param[in] offset
       * The offset of the transformation. This is added to any vertex after applying the rotation.
       */
      void transform(const VertexType& origin, const VertexType& angles, const VertexType& offset)
      {
        _vertex_set.transform(origin, angles, offset);
      }

      /**
       * \brief Returns the name of the class.
       * \returns
       * The name of the class as a String.
       */
      static String name()
      {
        return "ConformalMesh<"+Shape_::name()+","+stringify(num_coords_)+">";
      }
    }; // class ConformalMesh<...>

    /* ************************************************************************************************************* */

    /**
     * \brief Factory specialization for ConformalMesh class template
     *
     * \author Peter Zajac
     */
    template<
      typename Shape_,
      int num_coords_,
      typename CoordType_>
    class Factory< ConformalMesh<Shape_, num_coords_, CoordType_> >
    {
    public:
      /// mesh typedef
      typedef ConformalMesh<Shape_, num_coords_, CoordType_> MeshType;

      /// vertex set type
      typedef typename MeshType::VertexSetType VertexSetType;
      /// index holder type
      typedef typename MeshType::IndexSetHolderType IndexSetHolderType;

    public:
      /// virtual destructor
      virtual ~Factory()
      {
      }

      /**
       * \brief Returns the number of entities.
       *
       * \param[in] dim
       * The dimension of the entity whose count is to be returned. Must be 0 <= \p dim <= shape_dim.
       *
       * \returns
       * The number of entities of dimension \p dim.
       */
      virtual Index get_num_entities(int dim) = 0;

      /**
       * \brief Fills the vertex set.
       *
       * \param[in,out] vertex_set
       * The vertex set whose coordinates are to be filled.
       */
      virtual void fill_vertex_set(VertexSetType& vertex_set) = 0;

      /**
       * \brief Fills the index sets.
       *
       * \param[in,out] index_set_holder
       * The index set holder whose index sets are to be filled.
       */
      virtual void fill_index_sets(IndexSetHolderType& index_set_holder) = 0;

    }; // class Factory<ConformalMesh<...>>

    /* ************************************************************************************************************* */

    /**
     * \brief StandardRefinery implementation for ConformalMesh
     *
     * \author Peter Zajac
     */
    template<
      typename Shape_,
      int num_coords_,
      typename CoordType_>
    class StandardRefinery<ConformalMesh<Shape_, num_coords_, CoordType_> > :
      public Factory< ConformalMesh<Shape_, num_coords_, CoordType_> >
    {
    public:
      /// mesh type
      typedef ConformalMesh<Shape_, num_coords_, CoordType_> MeshType;
      /// shape type
      typedef typename MeshType::ShapeType ShapeType;
      /// vertex set type
      typedef typename MeshType::VertexSetType VertexSetType;
      /// index holder type
      typedef typename MeshType::IndexSetHolderType IndexSetHolderType;

      /// shape dimension
      static constexpr int shape_dim = ShapeType::dimension;

    protected:
      /// coarse mesh reference
      const MeshType& _coarse_mesh;
      /// number of entities for coarse mesh
      Index _num_entities_coarse[shape_dim + 1];
      /// number of entities for fine mesh
      Index _num_entities_fine[shape_dim + 1];

    public:
      /**
       * \brief Constructor.
       *
       * \param[in] coarse_mesh
       * A reference to the coarse mesh that is to be refined.
       */
      explicit StandardRefinery(const MeshType& coarse_mesh) :
        _coarse_mesh(coarse_mesh)
      {
        // get number of entities in coarse mesh
        for(int i(0); i <= shape_dim; ++i)
        {
          _num_entities_fine[i] = coarse_mesh.get_num_entities(i);
          _num_entities_coarse[i] = coarse_mesh.get_num_entities(i);
        }

        // calculate number of entities in fine mesh
        Intern::EntityCountWrapper<Intern::StandardRefinementTraits, ShapeType>::query(_num_entities_fine);
      }

      /// virtual destructor
      virtual ~StandardRefinery()
      {
      }

      /**
       * \brief Returns the number of entities of the refined mesh.
       *
       * \param[in] dim
       * The dimension of the entity whose count is to be returned. Must be 0 <= \p dim <= #shape_dim.
       *
       * \returns
       * The number of entities of dimension \p dim.
       */
      virtual Index get_num_entities(int dim) override
      {
        return _num_entities_fine[dim];
      }

      /**
       * \brief Fills the vertex set of the refined mesh.
       *
       * \param[in,out] vertex_set
       * The vertex set whose coordinates are to be filled.
       */
      virtual void fill_vertex_set(VertexSetType& vertex_set) override
      {
        // refine vertices
        Intern::StandardVertexRefineWrapper<ShapeType, VertexSetType>
          ::refine(vertex_set, _coarse_mesh.get_vertex_set(), _coarse_mesh.get_index_set_holder());
      }

      /**
       * \brief Fills the index sets of the refined mesh.
       *
       * \param[in,out] index_set_holder
       * The index set holder whose index sets are to be filled.
       */
      virtual void fill_index_sets(IndexSetHolderType& index_set_holder) override
      {
        // refine indices
        Intern::IndexRefineWrapper<ShapeType>
          ::refine(index_set_holder, _num_entities_coarse, _coarse_mesh.get_index_set_holder());
      }
    }; // class StandardRefinery<ConformalMesh<...>,Nil>

#ifdef FEAT_EICKT
    extern template class ConformalMesh<Shape::Simplex<2>, 2, Real>;
    extern template class ConformalMesh<Shape::Simplex<3>, 3, Real>;
    extern template class ConformalMesh<Shape::Hypercube<2>, 2, Real>;
    extern template class ConformalMesh<Shape::Hypercube<3>, 3, Real>;

    extern template class StandardRefinery<ConformalMesh<Shape::Simplex<2>, 2, Real>>;
    extern template class StandardRefinery<ConformalMesh<Shape::Simplex<3>, 3, Real>>;
    extern template class StandardRefinery<ConformalMesh<Shape::Hypercube<2>, 2, Real>>;
    extern template class StandardRefinery<ConformalMesh<Shape::Hypercube<3>, 3, Real>>;
#endif // FEAT_EICKT
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_CONFORMAL_MESH_HPP
