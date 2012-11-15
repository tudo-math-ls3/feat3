#pragma once
#ifndef KERNEL_GEOMETRY_CONFORMAL_SUB_MESH_HPP
#define KERNEL_GEOMETRY_CONFORMAL_SUB_MESH_HPP 1

// includes, FEAST
#include <kernel/geometry/factory.hpp>
#include <kernel/geometry/intern/standard_index_refiner.hpp>
#include <kernel/geometry/intern/standard_target_refiner.hpp>
#include <kernel/geometry/intern/standard_vertex_refiner.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /**
     * \brief Standard conformal sub-mesh policy
     *
     * This class defines a default policy for the ConformalSubMesh class template.
     *
     * \tparam Shape_
     * The shape that is to be used for the mesh. Must be either Shape::Simplex<n> or Shape::Hypercube<n>
     * for some \c n > 0.
     *
     * \tparam VertexSet_
     * The vertex set class to be used by the mesh. By default, VertexSetVariable is used.
     *
     * \author Peter Zajac
     */
    template<
      typename Shape_,
      typename VertexSet_ = VertexSetVariable<> >
    struct ConformalSubMeshPolicy
    {
      /// shape type
      typedef Shape_ ShapeType;

      /// Vertex set type
      typedef VertexSet_ VertexSetType;
    }; // struct ConformalSubMeshPolicy

    /**
     * \brief Conformal sub-mesh class template
     *
     * \author Peter Zajac
     */
    template<typename Policy_>
    class ConformalSubMesh
    {
      // friends
      template<typename Mesh_, typename Parent_>
      friend class StandardRefinery;

    public:
      /// policy type
      typedef Policy_ PolicyType;

      /// Shape type
      typedef typename PolicyType::ShapeType ShapeType;

      /// Vertex set type
      typedef typename PolicyType::VertexSetType VertexSetType;

      /// index set holder type
      typedef IndexSetHolder<ShapeType> IndexSetHolderType;

      /// target set holder type
      typedef TargetSetHolder<ShapeType> TargetSetHolderType;

      /// dummy enum
      enum
      {
        /// shape dimension
        shape_dim = ShapeType::dimension
      };

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
        typedef
          FEAST::Geometry::IndexSet<
            Shape::FaceTraits<
              typename Shape::FaceTraits<
                ShapeType,
                cell_dim_>
              ::ShapeType,
              face_dim_>
            ::count> Type;
      }; // struct IndexSet<...>

      /**
       * \brief Target set type class template.
       *
       * This nested class template is used to define the return type of the ConformalSubMesh::get_target_set()
       * function template.
       *
       * \tparam cell_dim_
       * The cell dimesnion parameter as passed to the ConformalSubMesh::get_target_set() function template.
       */
      template<int cell_dim_>
      struct TargetSet
      {
        /// target set type
        typedef FEAST::Geometry::TargetSet Type;
      }; // struct TargetSet<...>

    protected:
      /// number of entities for each shape dimension
      Index _num_entities[shape_dim + 1];

      /// the vertex set of the mesh
      VertexSetType _vertex_set;

      /// the index sets of the mesh
      IndexSetHolderType _index_set_holder;

      /// the target sets of the mesh.
      TargetSetHolderType _target_set_holder;

    private:
      ConformalSubMesh(const ConformalSubMesh&);
      ConformalSubMesh& operator=(const ConformalSubMesh&);

    public:
      /**
       * \brief Constructor
       *
       * \param[in] num_entities
       * An array of length at least #shape_dim + 1 holding the number of entities for each shape dimension.
       * Must not be \c nullptr.
       *
       * \param[in] num_coords
       * The number of coordinates per vertex. This parameter is passed to the constructor of the vertex set.
       *
       * \param[in] vertex_stride
       * The vertex stride. This parameter is passed to the constructor of the vertex set.
       */
      explicit ConformalSubMesh(const Index num_entities[],
        int num_coords = shape_dim,
        int vertex_stride = 0)
          :
        _vertex_set(num_entities[0], num_coords, vertex_stride),
        _index_set_holder(num_entities),
        _target_set_holder(num_entities)
      {
        CONTEXT(name() + "::ConformalSubMesh()");
        for(int i(0); i <= shape_dim; ++i)
        {
          _num_entities[i] = num_entities[i];
        }
      }

      /**
       * \brief Factory constructor
       *
       * \param[in] factory
       * The factory that is to be used to create the mesh.
       */
      explicit ConformalSubMesh(const Factory<ConformalSubMesh>& factory) :
        _vertex_set(factory.get_num_entities(0), factory.get_num_coords(), factory.get_vertex_stride()),
        _index_set_holder(Intern::NumEntitiesWrapper<shape_dim>(factory).num_entities),
        _target_set_holder(Intern::NumEntitiesWrapper<shape_dim>(factory).num_entities)
      {
        // compute entity counts
        Intern::NumEntitiesWrapper<shape_dim>::apply(factory, _num_entities);

        // fill vertex set
        factory.fill_vertex_set(_vertex_set);

        // fill index sets
        factory.fill_index_sets(_index_set_holder);

        // fill target sets
        factory.fill_target_sets(_target_set_holder);
      }

      /// virtual destructor
      virtual ~ConformalSubMesh()
      {
        CONTEXT(name() + "::~ConformalSubMesh()");
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
        CONTEXT(name() + "::get_num_entities()");
        ASSERT_(dim >= 0);
        ASSERT_(dim <= shape_dim);
        return _num_entities[dim];
      }

      /// Returns a reference to the vertex set of the mesh.
      VertexSetType& get_vertex_set()
      {
        CONTEXT(name() + "::get_vertex_set()");
        return _vertex_set;
      }

      /** \copydoc get_vertex_set() */
      const VertexSetType& get_vertex_set() const
      {
        CONTEXT(name() + "::get_vertex_set() [const]");
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
        CONTEXT(name() + "::get_index_set<" + stringify(cell_dim_) + "," + stringify(face_dim_) + ">()");
        return _index_set_holder.template get_index_set_wrapper<cell_dim_>().get_index_set<face_dim_>();
      }

      /** \copydoc get_index_set() */
      template<
        int cell_dim_,
        int face_dim_>
      const typename IndexSet<cell_dim_, face_dim_>::Type& get_index_set() const
      {
        CONTEXT(name() + "::get_index_set<" + stringify(cell_dim_) + "," + stringify(face_dim_) + ">() [const]");
        return _index_set_holder.template get_index_set_wrapper<cell_dim_>().get_index_set<face_dim_>();
      }

      /// \cond internal
      IndexSetHolderType& get_index_set_holder()
      {
        CONTEXT(name() + "::get_index_set_holder()");
        return _index_set_holder;
      }

      const IndexSetHolderType& get_index_set_holder() const
      {
        CONTEXT(name() + "::get_index_set_holder() [const]");
        return _index_set_holder;
      }
      /// \endcond

      /**
       * \brief Return the reference to a target set.
       *
       * \tparam cell_dim_
       * The dimension fo the entity whose target set is to be returned.
       *
       * \returns
       * A reference to the target set.
       */
      template<int cell_dim_>
      typename TargetSet<cell_dim_>::Type& get_target_set()
      {
        CONTEXT(name() + "::get_target_set<" + stringify(cell_dim_) + ">()");
        static_assert(cell_dim_ >= 0, "invalid cell dimension");
        static_assert(cell_dim_ <= shape_dim, "invalid cell dimension");
        return _target_set_holder.template get_target_set<cell_dim_>();
      }

      /** \copydoc get_target_set() */
      template<int cell_dim_>
      const typename TargetSet<cell_dim_>::Type& get_target_set() const
      {
        CONTEXT(name() + "::get_target_set<" + stringify(cell_dim_) + ">() [const]");
        static_assert(cell_dim_ >= 0, "invalid cell dimension");
        static_assert(cell_dim_ <= shape_dim, "invalid cell dimension");
        return _target_set_holder.template get_target_set<cell_dim_>();
      }

      /// \cond internal
      TargetSetHolderType& get_target_set_holder()
      {
        CONTEXT(name() + "::get_target_set_holder()");
        return _target_set_holder;
      }

      const TargetSetHolderType& get_target_set_holder() const
      {
        CONTEXT(name() + "::get_target_set_holder() [const]");
        return _target_set_holder;
      }
      /// \endcond

      /**
       * \brief Refines the submesh.
       *
       * This function applies the standard refinement algorithm onto the submesh and returns the refined submesh.
       *
       * \param[in] parent_mesh
       * A reference to the (coarse) parent mesh that this submesh refers to.
       *
       * \returns
       * A pointer to the refined submesh.
       */
      template<typename ParentMesh_>
      ConformalSubMesh* refine(const ParentMesh_& parent_mesh) const
      {
        CONTEXT(name() + "::refine()");
        return new ConformalSubMesh(StandardRefinery<ConformalSubMesh>(*this, parent_mesh));
      }

      /**
       * \brief Returns the name of the class.
       * \returns
       * The name of the class as a String.
       */
      static String name()
      {
        return "ConformalSubMesh<...>";
      }
    }; // class ConformalSubMesh<...>

    /* ************************************************************************************************************* */

    /**
     * \brief Factory specialisation for ConformalSubMesh class template
     *
     * \author Peter Zajac
     */
    template<typename MeshPolicy_>
    class Factory< ConformalSubMesh<MeshPolicy_> >
    {
    public:
      /// mesh typedef
      typedef ConformalSubMesh<MeshPolicy_> MeshType;

      /// vertex set type
      typedef typename MeshType::VertexSetType VertexSetType;
      /// index set holder type
      typedef typename MeshType::IndexSetHolderType IndexSetHolderType;
      /// target set holder type
      typedef typename MeshType::TargetSetHolderType TargetSetHolderType;

    public:
      /// virtual destructor
      virtual ~Factory()
      {
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
      virtual Index get_num_entities(int dim) const = 0;

      /**
       * \brief Returns the number of coorindates per vertex.
       */
      virtual int get_num_coords() const = 0;

      /**
       * \brief Returns the vertex stride.
       */
      virtual int get_vertex_stride() const = 0;

      /**
       * \brief Fills the vertex set.
       *
       * \param[in,out] vertex_set
       * The vertex set whose coordinates are to be filled.
       */
      virtual void fill_vertex_set(VertexSetType& vertex_set) const = 0;

      /**
       * \brief Fills the index sets.
       *
       * \param[in,out] index_set_holder
       * The index set holder whose index sets are to be filled.
       */
      virtual void fill_index_sets(IndexSetHolderType& index_set_holder) const = 0;

      /**
       * \brief Fills the target sets.
       *
       * \param[in,out] target_set_holder
       * The target set holder whose target sets are to be filled.
       */
      virtual void fill_target_sets(TargetSetHolderType& target_set_holder) const = 0;
    }; // class Factory<ConformalSubMesh<...>>

    /* ************************************************************************************************************* */

    /**
     * \brief StandardRefinery implementation for ConformalSubMesh
     *
     * \author Peter Zajac
     */
    template<typename MeshPolicy_, typename Parent_>
    class StandardRefinery<ConformalSubMesh<MeshPolicy_>, Parent_> :
      public Factory< ConformalSubMesh<MeshPolicy_> >
    {
    public:
      /// mesh type
      typedef ConformalSubMesh<MeshPolicy_> MeshType;
      /// parent type
      typedef Parent_ ParentType;
      /// shape type
      typedef typename MeshType::ShapeType ShapeType;
      /// vertex set type
      typedef typename MeshType::VertexSetType VertexSetType;
      /// index set holder type
      typedef typename MeshType::IndexSetHolderType IndexSetHolderType;
      /// target set holder type
      typedef typename MeshType::TargetSetHolderType TargetSetHolderType;

      /// dummy enum
      enum
      {
        /// shape dimension
        shape_dim = ShapeType::dimension
      };

    protected:
      /// coarse mesh reference
      const MeshType& _coarse_mesh;
      /// coarse parent reference
      const ParentType& _parent;
      /// number of entities for fine mesh
      Index _num_entities_fine[shape_dim + 1];
      /// number of entities in parent
      Index _num_entities_parent[shape_dim + 1];

    public:
      /**
       * \brief Constructor.
       *
       * \param[in] coarse_mesh
       * A reference to the coarse mesh that is to be refined.
       *
       * \param[in] parent
       * A reference to the coarse parent.
       */
      explicit StandardRefinery(const MeshType& coarse_mesh, const ParentType& parent) :
        _coarse_mesh(coarse_mesh),
        _parent(parent)
      {
        // get number of entities in coarse mesh
        for(int i(0); i <= shape_dim; ++i)
        {
          _num_entities_fine[i] = coarse_mesh.get_num_entities(i);
          _num_entities_parent[i] = parent.get_num_entities(i);
        }

        // calculate number of entities in fine mesh
        Intern::EntityCountWrapper<ShapeType>::query(_num_entities_fine);
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
      virtual Index get_num_entities(int dim) const
      {
        return _num_entities_fine[dim];
      }

      /**
       * \brief Returns the number of coorindates per vertex.
       */
      virtual int get_num_coords() const
      {
        return _coarse_mesh.get_vertex_set().get_num_coords();
      }

      /**
       * \brief Returns the vertex stride.
       */
      virtual int get_vertex_stride() const
      {
        return _coarse_mesh.get_vertex_set().get_stride();
      }

      /**
       * \brief Fills the vertex set of the refined mesh.
       *
       * \param[in,out] vertex_set
       * The vertex set whose coordinates are to be filled.
       */
      virtual void fill_vertex_set(VertexSetType& vertex_set) const
      {
        // refine vertices
        Intern::StandardVertexRefineWrapper<ShapeType, VertexSetType>
          ::refine(vertex_set, _coarse_mesh._vertex_set, _coarse_mesh._index_set_holder);
      }

      /**
       * \brief Fills the index sets.
       *
       * \param[in,out] index_set_holder
       * The index set holder whose index sets are to be filled.
       */
      virtual void fill_index_sets(IndexSetHolderType& index_set_holder) const
      {
        // refine indices
        Intern::IndexRefineWrapper<ShapeType>
          ::refine(index_set_holder, _coarse_mesh._num_entities, _coarse_mesh._index_set_holder);
      }

      /**
       * \brief Fills the target sets.
       *
       * \param[in,out] target_set_holder
       * The target set holder whose target sets are to be filled.
       */
      virtual void fill_target_sets(TargetSetHolderType& target_set_holder) const
      {
        // refine target indices
        Intern::TargetRefineWrapper<ShapeType>::refine(target_set_holder, _num_entities_parent,
          _coarse_mesh._target_set_holder, _coarse_mesh._index_set_holder, _parent.get_index_set_holder());
      }
    }; // class StandardRefinery<ConformalSubMesh<...>,...>
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_CONFORMAL_SUB_MESH_HPP
