#pragma once
#ifndef KERNEL_GEOMETRY_CONFORMAL_MESH_HPP
#define KERNEL_GEOMETRY_CONFORMAL_MESH_HPP 1

// includes, FEAT
#include <kernel/geometry/factory.hpp>
#include <kernel/geometry/intern/facet_neighbours.hpp>
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


    protected:
      /// number of entities for each shape dimension
      Index _num_entities[shape_dim + 1];

      /// the vertex set of the mesh
      VertexSetType _vertex_set;

      /// the index sets of the mesh
      IndexSetHolderType _index_set_holder;

      /// Information about cells sharing a facet
      typename IndexSet<shape_dim, shape_dim-1>::Type _neighbours;

    private:
      /// \brief Copy assignment operator declared but not implemented
      ConformalMesh& operator=(const ConformalMesh&);

    public:
      /**
       * \brief Constructor.
       *
       * \param[in] num_entities
       * An array of length at least #shape_dim + 1 holding the number of entities for each shape dimension.
       * Must not be \c nullptr.
       *
       * Up until now, every application has just one mesh ("root"), but this might change.
       */
      explicit ConformalMesh(const Index num_entities[]) :
        _vertex_set(num_entities[0]),
        _index_set_holder(num_entities),
        _neighbours(num_entities[shape_dim])
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
       *
       * Up until now, every application has just one mesh ("root"), but this might change.
       */
      explicit ConformalMesh(Factory<ConformalMesh>& factory) :
        _vertex_set(factory.get_num_entities(0)),
        _index_set_holder(Intern::NumEntitiesWrapper<shape_dim>(factory).num_entities),
        _neighbours(Intern::NumEntitiesWrapper<shape_dim>(factory).num_entities[shape_dim])
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

        // Fill neighbour information. This needs facet at cell information, so it needs to be called after
        // fill_index_sets() etc.
        fill_neighbours();

      }

      /**
       * \brief Copy Constructor
       *
       * \param[in] other
       * The conformal mesh that is to be copied.
       */
      ConformalMesh(const ConformalMesh& other) :
        _vertex_set(other.get_vertex_set()),
        _index_set_holder(other.get_index_set_holder()),
        _neighbours(other.get_neighbours())
      {
        for(int i(0); i <= shape_dim; ++i)
        {
          _num_entities[i] = other.get_num_entities(i);
        }
      }

      /// virtual destructor
      virtual ~ConformalMesh()
      {
      }

      /// \returns The size of dynamically allocated memory in bytes.
      std::size_t bytes() const
      {
        return _vertex_set.bytes() + _index_set_holder.bytes() + _neighbours.bytes();
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

      void fill_neighbours()
      {
        // Facet at cell index set
        auto& facet_idx = get_index_set<shape_dim, shape_dim -1>();

        XASSERTM(get_num_entities(shape_dim-1) == facet_idx.get_index_bound(), "mesh num_entities / index_set num_entities mismatch");

        if(_neighbours.get_num_indices() == Index(0))
          _neighbours = std::move(typename IndexSet<shape_dim, shape_dim-1>::Type(get_num_entities(shape_dim)));

        Intern::FacetNeighbours::compute(_neighbours, facet_idx);

      }

      /// \returns A reference to the facet neighbour relations
      typename IndexSet<shape_dim, shape_dim-1>::Type& get_neighbours()
      {
        return _neighbours;
      }

      /// \copydoc get_neighbours()
      const typename IndexSet<shape_dim, shape_dim-1>::Type& get_neighbours() const
      {
        return _neighbours;
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
        this->fill_neighbours();
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
     * \brief Factory specialisation for ConformalMesh class template
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
