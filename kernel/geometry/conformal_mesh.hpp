#pragma once
#ifndef KERNEL_GEOMETRY_CONFORMAL_MESH_HPP
#define KERNEL_GEOMETRY_CONFORMAL_MESH_HPP 1

// includes, FEAST
#include <kernel/geometry/target_set.hpp>
#include <kernel/geometry/vertex_set.hpp>
#include <kernel/geometry/conformal/index_set.hpp>
#include <kernel/geometry/conformal/standard_refinement/index_refine_wrappers.hpp>
#include <kernel/geometry/conformal/standard_refinement/target_refine_wrappers.hpp>
#include <kernel/geometry/conformal/standard_refinement/vertex_refine_wrappers.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /**
     * \brief Standard conformal mesh policy
     *
     * This class defines a default policy for the ConformalMesh class template.
     *
     * \tparam Shape_
     * The shape that is to be used for the mesh. Must be either Shape::Simplex<n> or Shape::Hypercube<n>
     * for some \c n > 0.
     *
     * \tparam VertexSet_
     * The vertex set class to be used by the mesh. By default, VertexSetFixed is used.
     *
     * \author Peter Zajac
     */
    template<
      typename Shape_,
      typename VertexSet_ = VertexSetFixed<Shape_::dimension> >
    struct ConformalMeshPolicy
    {
      /// shape type
      typedef Shape_ ShapeType;

      /// Vertex set type
      typedef VertexSet_ VertexSetType;
    }; // struct ConformalMeshPolicy

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
     * \brief Conformal mesh class template
     *
     * \todo detailed documentation
     *
     * \author Peter Zajac
     */
    template<typename Policy_>
    class ConformalMesh
    {
    public:
      /// policy type
      typedef Policy_ PolicyType;

      /// Shape type
      typedef typename PolicyType::ShapeType ShapeType;

      /// Vertex set type
      typedef typename PolicyType::VertexSetType VertexSetType;

      /// index set holder type
      typedef Conformal::IndexSetHolder<ShapeType> IndexSetHolderType;

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
          FEAST::Geometry::Conformal::IndexSet<
            Shape::FaceTraits<
              typename Shape::FaceTraits<
                ShapeType,
                cell_dim_>
              ::ShapeType,
              face_dim_>
            ::count> Type;
      }; // struct IndexSet<...>


    protected:
      /// number of entities for each shape dimension
      Index _num_entities[shape_dim + 1];

      /// the vertex set of the mesh
      VertexSetType _vertex_set;

      /// the index sets of the mesh
      IndexSetHolderType _index_set_holder;

    private:
      ConformalMesh(const ConformalMesh&);
      ConformalMesh& operator=(const ConformalMesh&);

    public:
      /**
       * \brief Constructor.
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
      explicit ConformalMesh(
        const Index num_entities[],
        int num_coords = shape_dim,
        int vertex_stride = 0)
          :
        _vertex_set(num_entities[0], num_coords, vertex_stride),
        _index_set_holder(num_entities)
      {
        CONTEXT(name() + "::ConformalMesh(const Index[])");
        for(int i(0); i <= shape_dim; ++i)
        {
          _num_entities[i] = num_entities[i];
        }
      }

      /// virtual destructor
      virtual ~ConformalMesh()
      {
        CONTEXT(name() + "::~ConformalMesh()");
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
       * \brief Refines the mesh.
       *
       * This function applies the standard refinement algorithm onto the mesh and returns the refined mesh.
       *
       * \returns
       * A pointer to the refined mesh.
       */
      ConformalMesh* refine() const
      {
        CONTEXT(name() + "::refine()");
        using namespace Conformal::StandardRefinement;

        // get number of coordinates and vertex stride
        int num_coords = _vertex_set.get_num_coords();
        int vertex_stride = _vertex_set.get_stride();

        // get number of entities in coarse mesh
        Index num_entities_fine[shape_dim + 1];
        for(int i(0); i <= shape_dim; ++i)
        {
          num_entities_fine[i] = _num_entities[i];
        }

        // calculate number of entities in fine mesh
        EntityCountWrapper<ShapeType>::query(num_entities_fine);

        // allocate a fine mesh
        ConformalMesh* fine_mesh = new ConformalMesh(num_entities_fine, num_coords, vertex_stride);

        // refine vertices
        VertexRefineWrapper<ShapeType, VertexSetType>::refine(fine_mesh->_vertex_set, _vertex_set, _index_set_holder);

        // refine indices
        IndexRefineWrapper<ShapeType>::refine(fine_mesh->_index_set_holder, _num_entities, _index_set_holder);

        // return fine mesh
        return fine_mesh;
      }

      /**
       * \brief Returns the name of the class.
       * \returns
       * The name of the class as a String.
       */
      static String name()
      {
        return "ConformalMesh<...>";
      }
    }; // class ConformalMesh<...>

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    /**
     * \brief Conformal sub-mesh class template
     *
     * \author Peter Zajac
     */
    template<typename Policy_>
    class ConformalSubMesh :
      public ConformalMesh<Policy_>
    {
    public:
      /// base class typedef
      typedef ConformalMesh<Policy_> BaseClass;

      /// target set holder type
      typedef TargetSetHolder<typename BaseClass::ShapeType> TargetSetHolderType;

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
        int num_coords /*= BaseClass::shape_dim*/,
        int vertex_stride = 0)
          :
        BaseClass(num_entities, num_coords, vertex_stride),
        _target_set_holder(num_entities)
      {
        CONTEXT(name() + "::ConformalSubMesh()");
      }

      /// virtual destructor
      virtual ~ConformalSubMesh()
      {
        CONTEXT(name() + "::~ConformalSubMesh()");
      }

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
        using namespace Conformal::StandardRefinement;

        typedef typename BaseClass::ShapeType ShapeType;
        typedef typename BaseClass::VertexSetType VertexSetType;

        // get number of coordinates and vertex stride
        int num_coords = BaseClass::_vertex_set.get_num_coords();
        int vertex_stride = BaseClass::_vertex_set.get_stride();

        // get number of entities in coarse mesh
        Index num_entities_fine[BaseClass::shape_dim + 1];
        Index num_entities_parent[BaseClass::shape_dim + 1];
        for(int i(0); i <= BaseClass::shape_dim; ++i)
        {
          num_entities_fine[i] = BaseClass::_num_entities[i];
          num_entities_parent[i] = parent_mesh.get_num_entities(i);
        }

        // calculate number of entities in fine mesh
        EntityCountWrapper<ShapeType>::query(num_entities_fine);

        // allocate a fine mesh
        ConformalSubMesh* fine_mesh = new ConformalSubMesh(num_entities_fine, num_coords, vertex_stride);

        // refine vertices
        VertexRefineWrapper<ShapeType, VertexSetType>::refine(fine_mesh->_vertex_set,
          this->_vertex_set, this->_index_set_holder);

        // refine indices
        IndexRefineWrapper<ShapeType>::refine(fine_mesh->_index_set_holder,
          this->_num_entities, this->_index_set_holder);

        // refine target indices
        TargetRefineWrapper<ShapeType>::refine(fine_mesh->_target_set_holder, num_entities_parent,
          this->_target_set_holder, this->_index_set_holder, parent_mesh.get_index_set_holder());

        // return fine mesh
        return fine_mesh;
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
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_CONFORMAL_MESH_HPP
