#pragma once
#ifndef KERNEL_GEOMETRY_MESH_PART_HPP
#define KERNEL_GEOMETRY_MESH_PART_HPP 1

// includes, FEAST
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/factory.hpp>
#include <kernel/geometry/intern/standard_index_refiner.hpp>
#include <kernel/geometry/intern/standard_subset_refiner.hpp>
#include <kernel/geometry/intern/standard_target_refiner.hpp>
#include <kernel/geometry/intern/standard_vertex_refiner.hpp>

namespace FEAST
{
  namespace Geometry
  {

    // Forward declarations
    namespace Intern
    {
      template<Index, Index>
      struct TargetSetComputer;
    }

    /// Alias for MeshAttribute
    template<typename A>
    using MeshAttribute = VertexSetVariable<A>;

    /**
     * \brief Class for holding MeshAttributes for all shape dimensions
     *
     * \tparam Shape_
     * Shape type of the highest dimensional cell type, i.e. 3 for Hypercube<3> meshes
     *
     * \tparam DataType_
     * DataType used in the attribute's values, i.e. double
     *
     * A MeshAttribute is some kind of data that refers to an entity type in a mesh, like vertices, edges, faces or
     * cells. The MeshAttributeHolder contains the set of attributes for its dimension (Shape_::dimension) and
     * inherits from its lower dimensional version. By this inheritance structure, it contains all MeshAttribute
     * information from dimension 0 (vertices) up to Shape_::dimension.
     *
     * A MeshAttribute's values themselves are vector valued in general. The most important use of a MeshAttribute
     * at the time of writing is the parametrisation variable in (boundary) submeshes. A d-dimensional boundary mesh
     * therefore can have a MeshAttribute of dimension 0 (meaning it refers to vertices) as parametrisation variable,
     * which in turn has values of dimension d.
     *
     * Up until now, only MeshAttributes of dimension 0 are implemented and used, although this might change in the
     * future.
     *
     * Note that a MeshAttribute needs to be refined if the mesh it refers to is refined. For dimension 0 attributes,
     * it is clear how to refine them: Linear interpolation. For attributes of the highest possible dimension
     * (meaning cell attributes), new cells inherit the value from their parent cells. For attributes of intermediate
     * dimension, this is not clear.
     *
     * \author Jordi Paul
     *
     */
    template<typename Shape_, typename DataType_>
    class MeshAttributeHolder :
      public MeshAttributeHolder<typename Shape::FaceTraits<Shape_, Shape_::dimension - 1>::ShapeType, DataType_>
    {
    public:
      /// Type of the highest dimensional entity
      typedef Shape_ ShapeType;
      /// Dimension of ShapeType
      static constexpr int shape_dim = ShapeType::dimension;
      /// Type the MeshAttribute uses
      typedef DataType_ DataType;
      /// Type of the Attribute
      typedef MeshAttribute<DataType_> AttributeType;
      /// Type for the set of all Attributes belonging to shape_dim
      typedef std::vector<MeshAttribute<DataType_>> AttributeSetType;

    protected:
      /// This type's base class, meaning a MeshAttributeHolder of one dimension less
      typedef MeshAttributeHolder
      <
        typename Shape::FaceTraits<ShapeType, shape_dim - 1>::ShapeType, DataType_
      > BaseClass;

      /// The set of all attributes belonging to shape_dim
      AttributeSetType _mesh_attributes;

    public:
      /// \brief Empty default constructor
      explicit MeshAttributeHolder()
      {
      }

      /// \brief Copy constructor
      MeshAttributeHolder(const MeshAttributeHolder& other) :
        BaseClass(other),
        _mesh_attributes(other._mesh_attributes)
      {
      }

      /// \brief Virtual destructor
      virtual ~MeshAttributeHolder()
      {
        CONTEXT(name() + "::~MeshAttributeHolder()");
      }

      /**
       * \brief Returns a reference to the set of attributes belonging to dimension dim_
       *
       * \tparam dim_
       * Shape dim of the mesh entity the attribute set retrieved refers to
       *
       */
      template<int dim_>
      AttributeSetType& get_mesh_attributes()
      {
        CONTEXT(name() + "::get_mesh_attributes()");
        static_assert(dim_ >= 0, "invalid dimension");
        static_assert(dim_ <= shape_dim, "invalid dimension");
        typedef typename Shape::FaceTraits<Shape_, dim_>::ShapeType CellType;
        return MeshAttributeHolder<CellType, DataType_>::_mesh_attributes;
      }

      /// \copydoc get_mesh_attributes()
      template<int dim_>
      const AttributeSetType& get_mesh_attributes() const
      {
        CONTEXT(name() + "::get_mesh_attributes() [const]");
        static_assert(dim_ >= 0, "invalid dimension");
        static_assert(dim_ <= shape_dim, "invalid dimension");
        typedef typename Shape::FaceTraits<Shape_, dim_>::ShapeType CellType;
        return MeshAttributeHolder<CellType, DataType_>::_mesh_attributes;
      }

      /**
       * \brief Returns the number of attributes in the attribute set of dimension dim
       *
       * \param[in] dim
       * Dimension the number of attributes is to be computed for
       *
       */
      virtual Index get_num_attributes(Index dim) const
      {
        CONTEXT(name() + "::get_num_attributes()");
        ASSERT(int(dim) <= shape_dim, "invalid dimension parameter");

        if(int(dim) == shape_dim)
          // If the requested dimension is mine...
          return Index(_mesh_attributes.size());

        // Otherwise recurse down
        return BaseClass::get_num_attributes(dim);
      }

      /// \brief Return the name of the class
      static String name()
      {
        return "MeshAttributeHolder<" + Shape_::name() + ">";
      }

      /**
       * \brief Adds one attribute of shape dimension dim to the corresponding set
       *
       * \param[in] attribute
       * The Attribute to be added.
       *
       * \param[in] dim
       * Shape dimension of attribute.
       *
       * \warning Checks whether attribute has the same number of entries as the corresponding mesh has entities of
       * dimension dim have to be performed by the caller!
       */
      virtual void add_attribute(const AttributeType& attribute, Index dim)
      {
        ASSERT(int(dim) <= shape_dim, "Attribute shape dim exceeds MeshPart shape dim!");

        if(int(dim) == shape_dim)
          // If the requested dimension is mine...
          _mesh_attributes.push_back(attribute);
        else
          // Otherwise recurse down
          BaseClass::add_attribute(attribute, dim);

      }
    }; // class MeshAttributeHolder<Shape_, DataType_>

    /// \cond internal
    /**
     * \brief Specialisation for dimension 0
     */
    template<typename DataType_>
    class MeshAttributeHolder<Shape::Vertex, DataType_>
    {
    public:
      typedef Shape::Vertex ShapeType;
      static constexpr int shape_dim = ShapeType::dimension;
      typedef std::vector<MeshAttribute<DataType_>> AttributeSetType;

    protected:
      AttributeSetType _mesh_attributes;

    public:
      explicit MeshAttributeHolder()
      {
        CONTEXT(name() + "::MeshAttributeHolder(const Index[])");
      }

      MeshAttributeHolder(const MeshAttributeHolder& other) :
        _mesh_attributes(other._mesh_attributes)
      {
      }

      virtual ~MeshAttributeHolder()
      {
        CONTEXT(name() + "::~MeshAttributeHolder()");
      }

      template<int dim_>
      AttributeSetType& get_mesh_attributes()
      {
        CONTEXT(name() + "::get_mesh_attributes()");
        static_assert(dim_ == 0, "invalid dimension");
        return _mesh_attributes;
      }

      template<int dim_>
      const AttributeSetType& get_mesh_attributes() const
      {
        CONTEXT(name() + "::get_mesh_attributes() [const]");
        static_assert(dim_ == 0, "invalid dimension");
        return _mesh_attributes;
      }

      virtual Index get_num_attributes(Index dim) const
      {
        CONTEXT(name() + "::get_num_attributes()");
#if defined DEBUG
        ASSERT(dim == 0, "invalid dimension parameter");
#else
        (void)dim;
#endif
        return Index(_mesh_attributes.size());
      }

      virtual void add_attribute(const MeshAttribute<DataType_>& attribute_, Index dim)
      {
#if defined DEBUG
        ASSERT(dim == 0, "Only attributes of shape dim 0 can be added to MeshParts of shape dim 0");
#else
        (void)dim;
#endif
        _mesh_attributes.push_back(attribute_);
      }

      //virtual void add_attribute(const MeshAttribute<DataType_>&& attribute_, Index dim)
      //{
      //  ASSERT(dim == 0, "Only attributes of shape dim 0 can be added to MeshParts of shape dim 0");
      //  _mesh_attributes.push_back(std::move(attribute_));
      //}

      static String name()
      {
        return "MeshAttributeHolder<Vertex>";
      }
    };
    /// \endcond

    /**
     * \brief Class template for partial meshes
     *
     * A MeshPart is a part of another mesh (called parent mesh) and is defined by mapping its own mesh entities
     * (like vertices, edges etc.) to the corresponding mesh entities of the parent mesh. For at least one arbitrary
     * shape dimension a mapping has to exist, but a MeshPart can have one mapping for each shape dimension of the
     * mesh it refers to.
     *
     * A MeshPart does not need to be connected or to have a topology, although it can implicitly use the parent's
     * topology. It can have sets of MeshAttributes for each shape dimension. If the parent mesh is refine, it is
     * possible to use that information to create a refined version of the MeshPart, although particular care has to
     * be taken when refining the MeshAttributes.
     *
     * A MeshPart can supply its own topology, which can be different from the parent's. One important example for
     * this is a MeshPart that defines the boundary of its parent mesh. If the boundary is parametrised using a
     * MeshAttribute and the boundary is closed, the parametrisation requires the doubling of some mesh entities to
     * correctly represent the parametrisation variable.
     *
     * \verbatim
     *
     *   Parent:      MeshPart representing the closed boundary:
     *   3----2       0--1--2--3--4
     *   |    |         map to
     *   |    |       0--1--2--3--0
     *   0----1
     *
     * \endverbatim
     *
     * Vertex 0 of the parent mesh has two different parametrisation variables, as it is both the first and the last
     * vertex of the closed boundary. This can be easily represented by having two vertices in the MeshPart point to
     * the same vertex of the parent mesh.
     *
     * At least one parent mapping has to be given, but it is possible to deduct the other parent mappings from this.
     * There are two possibilities:
     * 1. Bottom to top: Starting with the lowest dimensional parent mapping, an entity is considered to be in the
     *    MeshPart if all its sub shapes are in the MeshPart, i.e. for a given edge parent mapping in 3d, all faces
     *    consisting exclusively of edges in the MeshPart are considered to be in the MeshPart as well. This is then
     *    repeated for cells as well.
     * 2. Top to bottom: Starting with the highest dimensional parent mapping, all subshapes of these entities are
     *    considered to be in the MeshPart as well, i.e. for a given face parent mapping in 3d, all edges belonging
     *    to faces in the MeshPart are considered to be in the MeshPart as well.
     *
     * \warning Extra care has to be taken with this when it comes to refinement when adding higher dimensional
     * parent information. Consider the case of a the parent mesh partitioned into several patches and a MeshPart
     * consisting of two vertices where more than two patches meet. If the mesh is refined, there are still exactly
     * two vertices where more than two patches meet.
     * If by chance those two vertices share an edge in the coarse mesh and bottom to top parent deduction is used,
     * an edge is added to the MeshPart. Refinement of said edge then leads to the addition of a new vertex to the
     * MeshPart, which is correct in terms of mesh topology, but the MeshPart then no longer represents the set of
     * patch cross points.
     *
     */
    template<typename MeshType_>
    class MeshPart
    {
    };

    /**
     * \brief MeshPart Implementation for conformal meshes
     *
     * \tparam ShapeType_
     * Shape type of the ConformalMesh this Meshpart refers to.
     *
     * \tparam num_coords_
     * Number of coordinates per vertex
     *
     * \tparam stride_
     * Padded num_coords_
     *
     * \tparam Coord_
     * Data type for vertex coordinates, i.e. double
     *
     * \copydoc MeshPart
     *
     */
    template<typename ShapeType_, int num_coords_, int stride_, typename Coord_>
    class MeshPart<ConformalMesh<ShapeType_, num_coords_, stride_, Coord_>>
    {
      public:
        /// Shape type
        typedef ShapeType_ ShapeType;
        /// Mesh type
        typedef ConformalMesh<ShapeType> MeshType;
        /// Index set holder type
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;
        /// Target set holder type
        typedef TargetSetHolder<ShapeType> TargetSetHolderType;
        /// Data type for attributes
        typedef typename MeshType::VertexSetType::CoordType AttributeDataType;
        /// Type for mesh attributes
        typedef MeshAttribute<AttributeDataType> AttributeType;
        /// Mesh attribute holder type
        typedef MeshAttributeHolder<ShapeType, AttributeDataType> AttributeHolderType;
        /// Shape dimension
        static constexpr int shape_dim = ShapeType::dimension;

        /**
         * \brief Index set type class template
         *
         * This nested class template is used to define the return type of the ConformalMesh::get_index_set()
         * function template.
         *
         * \tparam cell_dim_, face_dim_
         * The cell and face dimension parameters as passed to the ConformalMesh::get_index_set() function template.
         */
        template<int cell_dim_, int face_dim_>
        struct IndexSet
        {
          static_assert(cell_dim_ <= shape_dim, "invalid cell dimension");
          static_assert(face_dim_ < cell_dim_, "invalid face/cell dimension");
          static_assert(face_dim_ >= 0, "invalid face dimension");

          /// Index set type
          typedef FEAST::Geometry::IndexSet
          <
            Shape::FaceTraits
            <
              typename Shape::FaceTraits<ShapeType, cell_dim_>::ShapeType,
              face_dim_
            >::count
          > Type;
        }; // struct IndexSet<...>

        /**
         * \brief Attribute set type class template
         *
         * This is used to define the return type of the get_attributes() function template.
         *
         * \tparam cell_dim_
         * Shape dim of the attribute set
         *
         * \tparam DataType_
         * Type used in members of the attribute set
         *
         */
        template<int cell_dim_, typename DataType_>
        struct AttributeSet
        {
          /// Attribute set type
          typedef std::vector<FEAST::Geometry::MeshAttribute<DataType_>> Type;
        };

        /**
         * \brief Target set type class template.
         *
         * This nested class template is used to define the return type of the MeshPart::get_target_set()
         * function template.
         *
         * \tparam cell_dim_
         * The cell dimension parameter as passed to the MeshPart::get_target_set() function template.
         */
        template<int cell_dim_>
        struct TargetSet
        {
          /// Target set type
          typedef FEAST::Geometry::TargetSet Type;
        }; // struct TargetSet<...>

      protected:
        /// Number of entities for each shape dimension
        Index _num_entities[shape_dim + 1];
        /// The vertex set of the mesh
        AttributeHolderType _attribute_holder;
        /// The index sets of the mesh
        IndexSetHolderType* _index_set_holder;
        /// The target sets of the mesh.
        TargetSetHolderType _target_set_holder;

      private:
        /// \brief Copy assignment operator
        MeshPart& operator=(const MeshPart&);

      public:
        /**
         * \brief Constructor
         *
         * \param[in] num_entities
         * An array of length at least #shape_dim + 1 holding the number of entities for each shape dimension.
         * Must not be \c nullptr.
         *
         * \param[in] create_topology
         * Determines if the MeshPart is to have a mesh topology.
         *
         */
        explicit MeshPart(const Index num_entities[], bool create_topology = false) :
          _attribute_holder(),
          _index_set_holder(nullptr),
          _target_set_holder(num_entities)
          {
            CONTEXT(name() + "::MeshPart()");
            for(int i(0); i <= shape_dim; ++i)
              _num_entities[i] = num_entities[i];

            if(create_topology)
              _index_set_holder = new IndexSetHolderType(num_entities);
          }

        /**
         * \brief Factory constructor
         *
         * \param[in] factory
         * The factory that is to be used to create the mesh.
         */
        explicit MeshPart(Factory<MeshPart>& factory) :
          _attribute_holder(),
          _index_set_holder(nullptr),
          _target_set_holder(Intern::NumEntitiesWrapper<shape_dim>(factory).num_entities)
          {
            CONTEXT(name() + "::MeshPart() [factory]");

            // Compute entity counts
            Intern::NumEntitiesWrapper<shape_dim>::apply(factory, _num_entities);

            // Fill target sets
            factory.fill_target_sets(_target_set_holder);

            // Fill index sets
            factory.fill_index_sets(_index_set_holder);

            // Fill attribute sets
            factory.fill_attribute_sets(_attribute_holder);
          }

        /**
         * \brief Copy Constructor
         *
         * \param[in] other
         * The MeshPart that is to be copied.
         */
        MeshPart(const MeshPart& other) :
          _attribute_holder(other.get_attribute_holder()),
          _index_set_holder(nullptr),
          _target_set_holder(other.get_target_set_holder())
          {
            CONTEXT(name() + "::MeshPart() [copy]");

            for(int i(0); i <= shape_dim; ++i)
              _num_entities[i] = other.get_num_entities(i);

            // Only create _index_set_holder if other has one as well
            if(other.has_topology())
              _index_set_holder = new IndexSetHolderType(*other._index_set_holder);
          }

        /// Virtual destructor
        virtual ~MeshPart()
        {
          CONTEXT(name() + "::~MeshPart()");

          if(_index_set_holder != nullptr)
            delete _index_set_holder;
        }

        /// \brief Checks if this MeshPart has a mesh topology
        bool has_topology() const
        {
          return (_index_set_holder != nullptr);
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

        /**
         * \brief Returns the AttributeSet belonging to dimension dim_
         *
         * \tparam dim_
         * Shape dimension we want to have the attribute set for
         */
        template<int dim_>
        typename AttributeSet<dim_, AttributeDataType>::Type& get_attributes()
        {
          CONTEXT(name() + "::get_attributes()");
          return _attribute_holder.template get_mesh_attributes<dim_>();
        }

        /// \copydoc get_attributes()
        template<int dim_>
        const typename AttributeSet<dim_, AttributeDataType>::Type& get_attributes() const
        {
          CONTEXT(name() + "::get_attributes( ) [const]");
          return _attribute_holder.template get_mesh_attributes<dim_>();
        }

        /**
         * \brief Computes the total number of Attributes in this MeshPart
         *
         * \returns The total number of Attributes in all attribute sets
         */
        Index get_num_attributes() const
        {
          Index num_attribute_sets(0);
          // Sum up the number of attributes over all attribute sets
          for(Index d(0); d < Index(shape_dim); ++d)
            num_attribute_sets += _attribute_holder.get_num_attributes(d);

          return num_attribute_sets;
        }

        /**
         * \brief Adds one attribute to this MeshPart
         *
         * \tparam dim_
         * Shape dimension of the Attribute to be added
         *
         * \param[in] attribute
         * Attribute to be added.
         */
        template<int dim_>
        void add_attribute(const AttributeType& attribute)
        {
          CONTEXT(name() + "::add_attribute()");
          // Check if the attribute to be added has the same number of entities as the MeshPart wrt. dim_
          ASSERT(attribute.get_num_vertices() == _target_set_holder.template get_target_set<dim_>().get_num_entities(), "Attribute/entity count mismatch!");
          _attribute_holder.template get_mesh_attributes<dim_>().push_back(attribute);
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
        template<int cell_dim_, int face_dim_>
        typename IndexSet<cell_dim_, face_dim_>::Type& get_index_set()
        {
          CONTEXT(name() + "::get_index_set<" + stringify(cell_dim_) + "," + stringify(face_dim_) + ">()");
          return _index_set_holder->template get_index_set_wrapper<cell_dim_>().template get_index_set<face_dim_>();
        }

        /** \copydoc get_index_set() */
        template<int cell_dim_, int face_dim_>
        const typename IndexSet<cell_dim_, face_dim_>::Type& get_index_set() const
        {
          CONTEXT(name() + "::get_index_set<" + stringify(cell_dim_) + "," + stringify(face_dim_) + ">() [const]");
          ASSERT(has_topology(), "Requested index_set of MeshPart without topology!");
          return _index_set_holder->template get_index_set_wrapper<cell_dim_>().template get_index_set<face_dim_>();
        }

        /// \cond internal

        /// Returns a reference to the attribute holder of the mesh.
        AttributeHolderType& get_attribute_holder()
        {
          CONTEXT(name() + "::get_attribute_holder()");
          return _attribute_holder;
        }

        /** \copydoc get_attribute_holder() */
        const AttributeHolderType& get_attribute_holder() const
        {
          CONTEXT(name() + "::get_attribute_holder() [const]");
          return _attribute_holder;
        }

        /// Returns a reference to the index set holder of the mesh.
        IndexSetHolderType* get_topology()
        {
          CONTEXT(name() + "::get_topology()");
          return _index_set_holder;
        }

        /// \copydoc get_topology()
        const IndexSetHolderType* get_topology() const
        {
          CONTEXT(name() + "::get_topology() [const]");
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

        /**
         * \cond internal
         *
         * Returns a reference to the target set holder of the mesh.
         */
        TargetSetHolderType& get_target_set_holder()
        {
          CONTEXT(name() + "::get_target_set_holder()");
          return _target_set_holder;
        }

        /// \copydoc get_target_set_holder()
        const TargetSetHolderType& get_target_set_holder() const
        {
          CONTEXT(name() + "::get_target_set_holder() [const]");
          return _target_set_holder;
        }
        /// \endcond

        /**
         * \brief Returns the name of the class.
         * \returns
         * The name of the class as a String.
         */
        static String name()
        {
          return "MeshPart<...>";
        }

        /**
         * \brief Fills the target sets from bottom to top
         *
         * \tparam end_dim_
         * Dimension to stop at (meaning the lowest dimension).
         *
         * \tparam current_dim_
         * Dimension to generate parent information for.
         *
         * \param[in,out] parent_mesh
         * Parent this MeshPart refers to.
         *
         */
        template<Index end_dim_, Index current_dim_ = ShapeType_::dimension>
        void compute_target_sets_from_bottom(const MeshType& parent_mesh)
        {
          Intern::template TargetSetComputer<end_dim_, current_dim_>::bottom_to_top(_target_set_holder, parent_mesh.get_index_set_holder());
        }

        /**
         * \brief Fills the target sets from top to bottom
         *
         * \tparam end_dim_
         * Dimension to stop at (meaning the highest dimension).
         *
         * \tparam current_dim_
         * Dimension to generate parent information for.
         *
         * \param[in,out] parent_mesh
         * Parent this MeshPart refers to.
         *
         */
        template<Index end_dim_, Index current_dim_ = 0>
        void compute_target_sets_from_top(const MeshType& parent_mesh)
        {
          Intern::template TargetSetComputer<end_dim_, current_dim_>::top_to_bottom(_target_set_holder, parent_mesh.get_index_set_holder());
        }

        /**
         * \brief Fills the mesh topology from parent information
         *
         * \todo: Not implemented yet
         */
        virtual void compute_index_set_holder()
        {
        }

    }; // class MeshPart<ConformalMesh<Shape_>>

    /**
     * \brief Factory specialisation for MeshPart class template
     *
     * \author Peter Zajac
     */
    template<typename MeshType_>
    class Factory< MeshPart<MeshType_> >
    {
      public:
        /// The shape type of the mesh
        typedef typename MeshType_::ShapeType ShapeType;
        /// Mesh typedef
        typedef MeshPart<MeshType_> MeshPartType;
        /// Data type for attributes
        typedef typename MeshType_::VertexSetType::CoordType AttributeDataType;
        /// Mesh attribute holder type
        typedef MeshAttributeHolder<ShapeType, AttributeDataType> AttributeHolderType;
        /// index set holder type
        typedef typename MeshPartType::IndexSetHolderType IndexSetHolderType;
        /// target set holder type
        typedef typename MeshPartType::TargetSetHolderType TargetSetHolderType;

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
         * \brief Fills the attribute sets.
         *
         * \param[in,out] attribute_set_holder
         * The attribute set holder whose attribute sets are to be filled.
         */
        virtual void fill_attribute_sets(AttributeHolderType& attribute_set_holder) = 0;

        /**
         * \brief Fills the index sets.
         *
         * \param[in,out] index_set_holder
         * The index set holder whose index sets are to be filled.
         */
        virtual void fill_index_sets(IndexSetHolderType*& index_set_holder) = 0;

        /**
         * \brief Fills the target sets.
         *
         * \param[in,out] target_set_holder
         * The target set holder whose target sets are to be filled.
         */
        virtual void fill_target_sets(TargetSetHolderType& target_set_holder) = 0;

    }; // class Factory<MeshPart<...>>

    /**
     * \brief StandardRefinery implementation for MeshPart
     *
     * \author Peter Zajac
     */
    template<typename Shape_, typename Parent_>
    class StandardRefinery<MeshPart<Shape_>, Parent_> :
    public Factory< MeshPart<Shape_> >
    {
      public:
        /// mesh type
        typedef MeshPart<Shape_> MeshType;
        /// parent type
        typedef Parent_ ParentType;
        /// shape type
        typedef typename MeshType::ShapeType ShapeType;
        /// Attribute set type
        typedef typename MeshType::AttributeType AttributeType;
        /// index set holder type
        typedef typename MeshType::AttributeHolderType AttributeHolderType;
        /// index set holder type
        typedef typename MeshType::IndexSetHolderType IndexSetHolderType;
        /// target set holder type
        typedef typename MeshType::TargetSetHolderType TargetSetHolderType;

        /// shape dimension
        static constexpr int shape_dim = ShapeType::dimension;

      protected:
        /// coarse mesh reference
        const MeshType& _coarse_mesh;
        /// coarse parent reference
        const ParentType& _parent;
        /// number of entities for coarse mesh
        Index _num_entities_coarse[shape_dim + 1];
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
              _num_entities_fine[i] = _num_entities_coarse[i] = coarse_mesh.get_num_entities(i);
              _num_entities_parent[i] = parent.get_num_entities(i);
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
         * The dimension of the entity whose count is to be returned. Must be 0 <= \p dim <= shape_dim.
         *
         * \returns
         * The number of entities of dimension \p dim.
         */
        virtual Index get_num_entities(int dim)
        {
          return _num_entities_fine[dim];
        }

        /**
         * \brief Fills attribute sets where applicable
         *
         * \param[in,out] attribute_set_holder
         * Container for the attribute sets being added to
         *
         */
        virtual void fill_attribute_sets(AttributeHolderType& attribute_set_holder)
        {
          // Attributes of shape dimension 0 only make sense if we have a mesh topology
          if(_coarse_mesh.has_topology())
          {
            // Iterate over the dim 0 attributes in the coarse mesh
            for(Index i(0); i < _coarse_mesh.get_attribute_holder().template get_mesh_attributes<0>().size(); ++i)
            {
              auto& coarse_attribute =_coarse_mesh.template get_attributes<0>()[i];

              // Create a new empty attribute of the desired size
              AttributeType refined_attribute(get_num_entities(0), coarse_attribute.get_num_coords());

              // Refine the attribute in the coarse mesh and write the result to the new attribute
              Intern::StandardVertexRefineWrapper<ShapeType, AttributeType>
                ::refine(refined_attribute, coarse_attribute, *_coarse_mesh.get_topology());

              // Add the attribute to the corresponding set
              attribute_set_holder.add_attribute(refined_attribute,0);
            }
          }
        }

        /**
         * \brief Fills the index sets.
         *
         * \param[in,out] index_set_holder
         * The index set holder whose index sets are to be filled.
         */
        virtual void fill_index_sets(IndexSetHolderType*& index_set_holder)
        {
          ASSERT(index_set_holder == nullptr, "fill_index_sets: index_set_holder != nullptr!");

          // Only create the topology for the refined mesh if the coarse mesh has a topology
          if(_coarse_mesh.has_topology())
          {
            index_set_holder = new IndexSetHolderType(_num_entities_fine);
            // refine indices
            Intern::IndexRefineWrapper<ShapeType>
              ::refine(*index_set_holder, _num_entities_coarse, *_coarse_mesh.get_topology());
          }
        }

        /**
         * \brief Fills the target sets.
         *
         * \param[in,out] target_set_holder
         * The target set holder whose target sets are to be filled.
         */
        virtual void fill_target_sets(TargetSetHolderType& target_set_holder)
        {
          // refine target indices
         const IndexSetHolderType* coarse_ish(_coarse_mesh.get_topology());

         // The refinement is different depending on the coarse mesh having a mesh topology
         if(_coarse_mesh.has_topology())
           Intern::TargetRefineWrapper<ShapeType>
             ::refine(target_set_holder, _num_entities_parent, _coarse_mesh.get_target_set_holder(),
             *coarse_ish, *_parent.get_topology());
         else
           Intern::SubSetRefineWrapper<ShapeType>
             ::refine(target_set_holder, _num_entities_parent, _coarse_mesh.get_target_set_holder());
        }

    }; // class StandardRefinery<MeshPart<...>,...>

    /// \cond internal
    namespace Intern
    {
      /**
       * \brief Wrapper class for filling TargetSetHolder objects
       *
       * \tparam end_dim_
       * Dimension to stop the template recursion at.
       *
       * \tparam current_dim_
       * Dimension for which a new container is filled.
       *
       * \author Jordi Paul
       *
       */
      template<Index end_dim_, Index current_dim_>
      struct TargetSetComputer
      {
        /**
         * \brief Fills a TargetSetHolder from bottom to top using information from an IndexSetHolder
         *
         * The IndexSetHolder contains the topology of the parent mesh that the TargetSetHolder refers to.
         *
         * \tparam TargetSetHolderType
         * Type of the TargetSetHolder to be filled.
         *
         * \tparam IndexSetHolderType
         * Type containing mesh topology information.
         *
         */
        template<typename TargetSetHolderType, typename IndexSetHolderType>
        static void bottom_to_top(TargetSetHolderType& _target_set_holder, const IndexSetHolderType& parent_ish)
        {
          // Template recursion: Call the lower dimensional version first
          TargetSetComputer<end_dim_, current_dim_-1>::template bottom_to_top(_target_set_holder, parent_ish);

          typedef typename IndexSetHolderType::template IndexSet<current_dim_, current_dim_-1>::Type ParentIndexSetType;

          // TargetSet to fill
          TargetSet& my_target_set(_target_set_holder.template get_target_set<current_dim_>());

          // Only do something if the target set is still empty
          if(my_target_set.get_num_entities() == 0)
          {
            // Lower dimensional target set known to exist
            TargetSet& ts_below(_target_set_holder.template get_target_set<current_dim_-1>());

            // Indexset for entities of dimension current_dim_-1 at entities of dimension current_dim_
            const ParentIndexSetType& is_parent_below(parent_ish.template get_index_set<current_dim_, current_dim_-1>());
            // IndexSet for storing the indices of the MeshPart entities lying at entities of the parent
            ParentIndexSetType lower_parent_to_mesh_part(Index(is_parent_below.get_num_indices()));

            // For every entity of the parent, this will save whether it is referenced by the TargetSet
            TargetSet lower_origin_set(is_parent_below.get_index_bound());
            // and this is the marker for that
            Index marker(is_parent_below.get_index_bound());
            for(Index i(0); i < lower_origin_set.get_num_entities(); ++i)
              lower_origin_set[i] = 0;
            // Set the values to something bogus for every index present
            for(Index i(0); i < ts_below.get_num_entities(); ++i)
              lower_origin_set[ts_below[i]] = marker;

            // Count the number of entities that get added
            Index num_entities_current_dim(0);

            // Temporary TargetSet initialised with the maximum possible size (= number of entities in the parent mesh,
            // which is the case if the MeshPart contains all entities of the parent mesh)
            TargetSet tmp_target_set(is_parent_below.get_num_entities());

            // Check all entities of dimension current_dim_ from the parent
            for(Index i(0); i < is_parent_below.get_num_entities(); ++i)
            {
              bool is_in_mesh_part(true);

              // Check if all entities of dimension current_dim_-1 at entity i are referenced by the TargetSet
              for(Index j(0); j < ParentIndexSetType::num_indices; ++j)
              {
                // This is the index in the parent
                Index parent_index(is_parent_below[i][j]);
                if(lower_origin_set[parent_index] != marker)
                  is_in_mesh_part = false;
              }

              // If all subshapes belonged to the MeshPart, create the new parent mapping
              if(is_in_mesh_part)
              {
                tmp_target_set[num_entities_current_dim] = i;
                num_entities_current_dim++;
              }
            }

            // tmp_target_set possibly has the wrong size, so manually create a correctly sized TargetSet
            TargetSet new_target_set(num_entities_current_dim);
            for(Index i(0); i < new_target_set.get_num_entities(); ++i)
              new_target_set[i] = tmp_target_set[i];

            // Update the information in the TargetSetHolder
            my_target_set = std::move(new_target_set);
          }
        } // void bottom_to_top<typename, typename>()

        /**
         * \brief Fills a TargetSetHolder from top to bottom using information from an IndexSetHolder
         *
         * The IndexSetHolder contains the topology of the parent mesh that the TargetSetHolder refers to.
         *
         * \tparam TargetSetHolderType
         * Type of the TargetSetHolder to be filled.
         *
         * \tparam IndexSetHolderType
         * Type containing mesh topology information.
         *
         */
        template<typename TargetSetHolderType, typename IndexSetHolderType>
        static void top_to_bottom(TargetSetHolderType& _target_set_holder, const IndexSetHolderType& parent_ish)
        {
          // Call higher dimensional version first
          TargetSetComputer<end_dim_, current_dim_+1>::template top_to_bottom(_target_set_holder, parent_ish);

          typedef typename IndexSetHolderType::template IndexSet<current_dim_+1, current_dim_>::Type ParentIndexSetType;

          // TargetSet to fill
          TargetSet& my_target_set(_target_set_holder.template get_target_set<current_dim_>());

          // Only do something if the target set is still empty
          if(my_target_set.get_num_entities() == 0)
          {
            // Higher dimensional target set known to exist
            TargetSet& ts_above(_target_set_holder.template get_target_set<current_dim_+1>());

            // Indexset for entities of dimension current_dim_ at entities of dimension current_dim_+1
            const ParentIndexSetType& is_parent_above(parent_ish.template get_index_set<current_dim_+1, current_dim_>());
            // IndexSet for storing the indices of the MeshPart entities lying at entities of the parent
            ParentIndexSetType upper_parent_to_mesh_part(Index(is_parent_above.get_num_indices()));

            // For every entity of current_dim_ of the parent, this will save whether it belongs to an entity of
            // dimensin current_dim_+1 that is referenced through the TargetSetHolder
            TargetSet upper_origin_set(is_parent_above.get_index_bound());
            // And this is the marker for that
            Index marker(is_parent_above.get_index_bound());
            for(Index i(0); i < upper_origin_set.get_num_entities(); ++i)
              upper_origin_set[i] = 0;

            for(Index i(0); i < ts_above.get_num_entities(); ++i)
            {
              // Index of the current entity in the parent
              Index parent_index(ts_above[i]);
              for(Index j(0); int(j) < is_parent_above.get_num_indices(); ++j)
              {
                // This is i.e. the global number of local edge j at face parent_index in the parent
                Index parent_sub_entity_index(is_parent_above[parent_index][j]);
                upper_origin_set[parent_sub_entity_index] = marker;
              }
            }

            // Temporary TargetSet initialised with the maximum possible size (= number of entities in the parent mesh,
            // which is the case if the MeshPart contains all entities of the parent mesh)
            TargetSet tmp_target_set(is_parent_above.get_num_entities());
            // Count the number of entities that get added
            Index num_entities_current_dim(0);
            for(Index i(0); i < upper_origin_set.get_num_entities(); ++i)
            {
              if(upper_origin_set[i] == marker)
              {
                tmp_target_set[num_entities_current_dim] = i;
                num_entities_current_dim++;
              }
            }

            // tmp_target_set possibly has the wrong size, so manually create a correctly sized TargetSet
            TargetSet new_target_set(num_entities_current_dim);
            for(Index i(0); i < new_target_set.get_num_entities(); ++i)
              new_target_set[i] = tmp_target_set[i];

            // Update the information in the TargetSetHolder
            my_target_set = std::move(new_target_set);
          }

        }
      }; // struct TargetSetComputer<Index, Index>

      /**
       * \brief Specialisation as end of template recursion
       */
      template<Index end_dim_>
      struct TargetSetComputer<end_dim_, end_dim_>
      {
        /// \brief End of template recursion: Nothing to do
        template<typename TargetSetHolderType, typename IndexSetHolderType>
        static void bottom_to_top(TargetSetHolderType& , const IndexSetHolderType&)
        {
        }

        /// \brief End of template recursion: Nothing to do
        template<typename TargetSetHolderType, typename IndexSetHolderType>
        static void top_to_bottom(TargetSetHolderType& , const IndexSetHolderType&)
        {
        }
      }; // struct TargetSetComputer<Index>
    } // namespace Intern
    /// \endcond

  } // namespace Geometry
} // namespace FEAST
#endif
