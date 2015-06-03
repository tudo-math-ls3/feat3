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

    template<typename A>
    using MeshAttribute = VertexSetVariable<A>;

    /// \cond internal
    template<typename Shape_, typename DataType_>
    class MeshAttributeHolder :
      public MeshAttributeHolder<typename Shape::FaceTraits<Shape_, Shape_::dimension - 1>::ShapeType, DataType_>
    {
    public:
      typedef Shape_ ShapeType;
      static constexpr int shape_dim = ShapeType::dimension;
      typedef DataType_ DataType;
      typedef MeshAttribute<DataType_> AttributeType;
      typedef std::vector<MeshAttribute<DataType_>> AttributeSetType;

    protected:
      typedef MeshAttributeHolder<typename Shape::FaceTraits<ShapeType, shape_dim - 1>::ShapeType, DataType_> BaseClass;

      AttributeSetType _mesh_attributes;

    public:
      explicit MeshAttributeHolder()
      {
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
        static_assert(dim_ >= 0, "invalid dimension");
        static_assert(dim_ <= shape_dim, "invalid dimension");
        typedef typename Shape::FaceTraits<Shape_, dim_>::ShapeType CellType;
        return MeshAttributeHolder<CellType, DataType_>::_mesh_attributes;
      }

      template<int dim_>
      const AttributeSetType& get_mesh_attributes() const
      {
        CONTEXT(name() + "::get_mesh_attributes() [const]");
        static_assert(dim_ >= 0, "invalid dimension");
        static_assert(dim_ <= shape_dim, "invalid dimension");
        typedef typename Shape::FaceTraits<Shape_, dim_>::ShapeType CellType;
        return MeshAttributeHolder<CellType, DataType_>::_mesh_attributes;
      }

      virtual Index get_num_attributes(Index dim) const
      {
        CONTEXT(name() + "::get_num_attributes()");
        ASSERT(int(dim) <= shape_dim, "invalid dimension parameter");
        if(dim == shape_dim)
        {
          return _mesh_attributes.size();
        }
        return BaseClass::get_num_attributes(dim);
      }

      static String name()
      {
        return "MeshAttributeHolder<" + Shape_::name() + ">";
      }

      virtual void add_attribute(const AttributeType& attribute, Index dim)
      {
        ASSERT(dim <= shape_dim, "Attribute shape dim exceeds MeshPart shape dim!");

        std::cout << "AttributeSetHolder<" << shape_dim << ">::add_attribute(...," << dim <<")" << std::endl;

        if(dim == shape_dim)
          _mesh_attributes.push_back(attribute);
        else
          BaseClass::add_attribute(attribute, dim);

      }
    };

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
        return _mesh_attributes.size();
      }

      virtual void add_attribute(const MeshAttribute<DataType_>& attribute_, Index dim)
      {
        ASSERT(dim == 0, "Only attributes of shape dim 0 can be added to MeshParts of shape dim 0");
        _mesh_attributes.push_back(attribute_);
      }

      virtual void add_attribute(const MeshAttribute<DataType_>&& attribute_, Index dim)
      {
        ASSERT(dim == 0, "Only attributes of shape dim 0 can be added to MeshParts of shape dim 0");
        _mesh_attributes.push_back(std::move(attribute_));
      }

      static String name()
      {
        return "MeshAttributeHolder<Vertex>";
      }
    };
    /// \endcond

    template<typename MeshType_>
    class MeshPart
    {

    };

    template<typename ShapeType_>
    class MeshPart<ConformalMesh<ShapeType_>>
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

        template<int cell_dim_, typename DataType_>
        struct AttributeSet
        {
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
        MeshPart& operator=(const MeshPart&);

      public:
        /**
         * \brief Constructor
         *
         * \param[in] num_entities
         * An array of length at least #shape_dim + 1 holding the number of entities for each shape dimension.
         * Must not be \c nullptr.
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

        template<int dim_>
        typename AttributeSet<dim_, AttributeDataType>::Type& get_attributes()
        {
          CONTEXT(name() + "::get_attributes()");
          return _attribute_holder.template get_mesh_attributes<dim_>();
        }

        template<int dim_>
        const typename AttributeSet<dim_, AttributeDataType>::Type& get_attributes() const
        {
          CONTEXT(name() + "::get_attributes( ) [const]");
          return _attribute_holder.template get_mesh_attributes<dim_>();
        }

        Index get_num_attributes() const
        {
          Index num_attribute_sets(0);
          for(Index d(0); d < Index(shape_dim); ++d)
            num_attribute_sets += _attribute_holder.get_num_attributes(d);
          return num_attribute_sets;
        }

        template<int dim_>
        void add_attribute(const AttributeType& attribute)
        {
          CONTEXT(name() + "::get_attributes()");
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

        IndexSetHolderType* get_topology()
        {
          CONTEXT(name() + "::get_topology()");
          return _index_set_holder;
        }

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
         * \brief Returns the name of the class.
         * \returns
         * The name of the class as a String.
         */
        static String name()
        {
          return "MeshPart<...>";
        }

        template<Index shape_dim_>
        void deduct_target_sets_from_bottom(const MeshType& parent_mesh)
        {
          if(shape_dim > 0)
            deduct_target_sets_from_bottom<shape_dim_-1>(parent_mesh);

          typename TargetSet<shape_dim_>::Type my_target_set(_target_set_holder.template get_target_set<shape_dim_>);

        }

        virtual void deduct_target_sets_from_top()
        {
        }

        virtual void deduct_index_set_holder()
        {
        }

    }; // class MeshPart

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

        virtual void fill_attribute_sets(AttributeHolderType& attribute_set_holder)
        {
          if(_coarse_mesh.has_topology())
          {
            for(Index i(0); i < _coarse_mesh.get_attribute_holder().template get_mesh_attributes<0>().size(); ++i)
            {
              auto& coarse_attribute =_coarse_mesh.template get_attributes<0>()[i];

              AttributeType refined_attribute(get_num_entities(0), coarse_attribute.get_num_coords());

              Intern::StandardVertexRefineWrapper<ShapeType, AttributeType>
                ::refine(refined_attribute, coarse_attribute, *_coarse_mesh.get_topology());

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

         if(_coarse_mesh.has_topology())
           Intern::TargetRefineWrapper<ShapeType>
             ::refine(target_set_holder, _num_entities_parent, _coarse_mesh.get_target_set_holder(),
             *coarse_ish, *_parent.get_topology());
         else
           Intern::SubSetRefineWrapper<ShapeType>
             ::refine(target_set_holder, _num_entities_parent, _coarse_mesh.get_target_set_holder());
        }

    }; // class StandardRefinery<MeshPart<...>,...>
  } // namespace Geometry
} // namespace FEAST
#endif
