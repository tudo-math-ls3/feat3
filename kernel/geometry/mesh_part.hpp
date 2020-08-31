// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_MESH_PART_HPP
#define KERNEL_GEOMETRY_MESH_PART_HPP 1

// includes, FEAT
#include <kernel/geometry/factory.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/structured_mesh.hpp>
#include <kernel/geometry/index_calculator.hpp>
#include <kernel/geometry/attribute_set.hpp>
#include <kernel/geometry/intern/simple_target_refiner.hpp>
#include <kernel/geometry/intern/standard_target_refiner.hpp>
#include <kernel/geometry/intern/structured_target_refiner.hpp>
#include <kernel/geometry/intern/standard_attrib_refiner.hpp>
#include <kernel/geometry/intern/standard_index_refiner.hpp>
#include <kernel/geometry/intern/standard_vertex_refiner.hpp>
#include <kernel/geometry/intern/index_set_filler.hpp>
#include <kernel/geometry/intern/target_set_computer.hpp>

// includes, system
#include<map>

namespace FEAT
{
  namespace Geometry
  {
    /**
     * \brief Class template for partial meshes
     *
     * A MeshPart is a part of another mesh (called parent mesh) and is defined by mapping its own mesh entities
     * (like vertices, edges etc.) to the corresponding mesh entities of the parent mesh. For at least one arbitrary
     * shape dimension a mapping has to exist, but a MeshPart can have one mapping for each shape dimension of the
     * mesh it refers to.
     *
     * A MeshPart does not need to be connected or to have a topology, although it can implicitly use the parent's
     * topology. It can have sets of AttributeSets for each shape dimension. If the parent mesh is refine, it is
     * possible to use that information to create a refined version of the MeshPart, although particular care has to
     * be taken when refining the AttributeSets.
     *
     * A MeshPart can supply its own topology, which can be different from the parent's. One important example for
     * this is a MeshPart that defines the boundary of its parent mesh. If the boundary is parametrised using a
     * AttributeSet and the boundary is closed, the parametrisation requires the doubling of some mesh entities to
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
     * \note This template works for both ConformalMesh as well as StructuredMesh meshes.
     *
     * \author Jordi Paul
     * \author Peter Zajac
     */
    template<typename MeshType_>
    class MeshPart
    {
      public:
        /// parent Mesh type
        typedef MeshType_ MeshType;
        /// Shape type
        typedef typename MeshType::ShapeType ShapeType;
        /// Topology (aka IndexSetHolder) of parent mesh type
        typedef typename MeshType::IndexSetHolderType ParentIndexSetHolderType;
        /// Index set holder type
        typedef IndexSetHolder<ShapeType> IndexSetHolderType;
        /// Target set holder type
        typedef TargetSetHolder<ShapeType> TargetSetHolderType;
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
          typedef FEAT::Geometry::IndexSet
          <
            Shape::FaceTraits
            <
              typename Shape::FaceTraits<ShapeType, cell_dim_>::ShapeType,
              face_dim_
            >::count
          > Type;
        }; // struct IndexSet<...>

        /// Data type for attributes
        typedef typename MeshType::VertexSetType::CoordType AttributeDataType;
        /// Type for mesh attributes
        typedef AttributeSet<AttributeDataType> AttributeSetType;

        /// submesh node bin container type
        typedef std::map<String, AttributeSetType*> AttributeSetContainer;
        /// submesh node iterator type
        typedef typename AttributeSetContainer::iterator AttributeSetIterator;
        /// submesh node const-iterator type
        typedef typename AttributeSetContainer::const_iterator AttributeSetConstIterator;
        /// submesh node reverse-iterator type
        typedef typename AttributeSetContainer::reverse_iterator AttributeSetReverseIterator;

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
          typedef FEAT::Geometry::TargetSet Type;
        }; // struct TargetSet<...>

      protected:
        /// Number of entities for each shape dimension
        Index _num_entities[shape_dim + 1];
        /// The index sets of the mesh
        IndexSetHolderType* _index_set_holder;
        /// The target sets of the mesh.
        TargetSetHolderType _target_set_holder;
        /// The attribute sets of the mesh
        AttributeSetContainer _mesh_attributes;

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
         */
        explicit MeshPart(const Index num_entities[], bool create_topology = false) :
          _index_set_holder(nullptr),
          _target_set_holder(num_entities),
          _mesh_attributes()
        {
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
          _index_set_holder(nullptr),
          _target_set_holder(Intern::NumEntitiesWrapper<shape_dim>(factory).num_entities),
          _mesh_attributes()
        {
          // Compute entity counts
          Intern::NumEntitiesWrapper<shape_dim>::apply(factory, _num_entities);

          // Fill target sets
          factory.fill_target_sets(_target_set_holder);

          // Fill index sets
          factory.fill_index_sets(_index_set_holder);

          // Fill attribute sets
          factory.fill_attribute_sets(_mesh_attributes);
        }

        /// move constructor
        MeshPart(MeshPart&& other) :
          _index_set_holder(other._index_set_holder),
          _target_set_holder(std::forward<TargetSetHolderType>(other._target_set_holder)),
          _mesh_attributes(std::forward<AttributeSetContainer>(other._mesh_attributes))
        {
          for(int i(0); i <= shape_dim; ++i)
          {
            _num_entities[i] = other._num_entities[i];
          }
          other._index_set_holder = nullptr;
        }

        /// move-assignment operator
        MeshPart& operator=(MeshPart&& other)
        {
          // avoid self-move
          if(this == &other)
            return *this;

          if(_index_set_holder != nullptr)
            delete _index_set_holder;

          _index_set_holder = other._index_set_holder;
          _target_set_holder = std::forward<TargetSetHolderType>(other._target_set_holder);
          _mesh_attributes = std::forward<AttributeSetContainer>(other._mesh_attributes);
          for(int i(0); i <= shape_dim; ++i)
          {
            _num_entities[i] = other._num_entities[i];
          }
          other._index_set_holder = nullptr;

          return *this;
        }

        /// Virtual destructor
        virtual ~MeshPart()
        {
          if(_index_set_holder != nullptr)
            delete _index_set_holder;

          // Loop over all mesh attributes in reverse order and delete them
          AttributeSetReverseIterator it(_mesh_attributes.rbegin());
          AttributeSetReverseIterator jt(_mesh_attributes.rend());
          for(; it != jt; ++it)
          {
            if(it->second != nullptr)
              delete it->second;
          }
        }

        /// \returns The size of dynamically allocated memory in bytes.
        std::size_t bytes() const
        {
          std::size_t my_bytes(0);

          for(const auto& it : _mesh_attributes)
          {
            if(it.second != nullptr)
              my_bytes += it.second->bytes();
          }

          my_bytes += _target_set_holder.bytes();
          my_bytes +=(_index_set_holder != nullptr ? _index_set_holder->bytes() : std::size_t(0));

          return my_bytes;
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
          XASSERT(dim >= 0);
          XASSERT(dim <= shape_dim);
          return _num_entities[dim];
        }

        /**
         * \brief Finds an Attribute by its identifier String and returns a pointer to it
         *
         * \param[in] identifier
         * Identifier to match.
         *
         * \returns A pointer to the Attribute of shape dimension dim if present.
         *
         * \warning Will return nullptr if no Attribute with the given identifier is found
         *
         */
        AttributeSetType* find_attribute(const String& identifier)
        {
          AttributeSetIterator it(_mesh_attributes.find(identifier));
          return (it != _mesh_attributes.end()) ? (*it).second: nullptr;
        }

        /** \copydoc find_attribute() */
        const AttributeSetType* find_attribute(const String& identifier) const
        {
          AttributeSetConstIterator it(_mesh_attributes.find(identifier));
          return (it != _mesh_attributes.end()) ? (*it).second: nullptr;
        }

        /**
         * \brief Returns a const reference to the mesh attributes container
         *
         * \returns A const reference to the mesh attributes container
         */
        const AttributeSetContainer& get_mesh_attributes() const
        {
          return _mesh_attributes;
        }

        /**
         * \brief Returns the total number of Attributes in this MeshPart
         *
         * \returns The total number of Attributes.
         */
        int get_num_attributes() const
        {
          return int(_mesh_attributes.size());
        }

        /**
         * \brief Copies one attribute to this MeshPart's AttributeHolder
         *
         * \param[in] attribute
         * Attribute to be added.
         *
         * \param[in] identifier
         * Identifier for the map.
         *
         * \returns
         * True if the attribute was successfully added, meaning no attribute with the appropriate identifier was
         * present, or false otherwise.
         *
         */
        virtual bool add_attribute(AttributeSetType* attribute, const String& identifier)
        {
          if(attribute != nullptr || (attribute->get_num_values() != get_num_entities(0)) )
          {
            return (_mesh_attributes.insert( std::make_pair(identifier, attribute))).second;
          }

          return false;
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
          return _index_set_holder->template get_index_set_wrapper<cell_dim_>().template get_index_set<face_dim_>();
        }

        /** \copydoc get_index_set() */
        template<int cell_dim_, int face_dim_>
        const typename IndexSet<cell_dim_, face_dim_>::Type& get_index_set() const
        {
          XASSERTM(has_topology(), "Requested index_set of MeshPart without topology!");
          return _index_set_holder->template get_index_set_wrapper<cell_dim_>().template get_index_set<face_dim_>();
        }

        /// \cond internal

        /// Returns a reference to the index set holder of the mesh.
        IndexSetHolderType* get_topology()
        {
          return _index_set_holder;
        }

        /// \copydoc get_topology()
        const IndexSetHolderType* get_topology() const
        {
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
          static_assert(cell_dim_ >= 0, "invalid cell dimension");
          static_assert(cell_dim_ <= shape_dim, "invalid cell dimension");
          return _target_set_holder.template get_target_set<cell_dim_>();
        }

        /** \copydoc get_target_set() */
        template<int cell_dim_>
        const typename TargetSet<cell_dim_>::Type& get_target_set() const
        {
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
          return _target_set_holder;
        }

        /// \copydoc get_target_set_holder()
        const TargetSetHolderType& get_target_set_holder() const
        {
          return _target_set_holder;
        }
        /// \endcond

        /**
         * \brief Returns the name of the class.
         *
         * \returns
         * The name of the class as a String.
         */
        static String name()
        {
          return "MeshPart<...>";
        }

        /**
         * \brief Deducts the target sets from bottom to top
         *
         * \tparam end_dim_
         * Dimension to stop at (meaning the lowest dimension).
         *
         * \tparam current_dim_
         * Dimension to generate parent information for.
         *
         * \param[in] parent_ish
         * Topology of the parent this MeshPart refers to.
         *
         * \warning This will in general change the implied topology of the MeshPart and meant to be used in mesh
         * preprocessing only.
         *
         */
        template<int end_dim_, int current_dim_ = ShapeType::dimension>
        void deduct_target_sets_from_bottom(const ParentIndexSetHolderType& parent_ish)
        {
          Intern::template TargetSetComputer<end_dim_, current_dim_>::bottom_to_top(_target_set_holder, parent_ish);

          // Update num_entities information from the possibly modified _target_set_holder
          for(int i(end_dim_); i <= current_dim_; ++i)
            _num_entities[i] = _target_set_holder.get_num_entities(i);
        }

        /**
         * \brief Deducts the target sets from top to bottom
         *
         * \tparam end_dim_
         * Dimension to stop at (meaning the highest dimension).
         *
         * \tparam current_dim_
         * Dimension to generate parent information for.
         *
         * \param[in] parent_ish
         * Topology of the parent this MeshPart refers to.
         *
         * \warning This will in general change the implied topology of the MeshPart and meant to be used in mesh
         * preprocessing only.
         *
         */
        template<int end_dim_, int current_dim_ = 0>
        void deduct_target_sets_from_top(const ParentIndexSetHolderType& parent_ish)
        {
          Intern::template TargetSetComputer<end_dim_, current_dim_>::top_to_bottom(_target_set_holder, parent_ish);

          // Update num_entities information from the possibly modified _target_set_holder
          for(int i(current_dim_); i <= end_dim_; ++i)
            _num_entities[i] = _target_set_holder.get_num_entities(i);
        }

        /**
         * \brief Fills the mesh topology from parent information
         *
         * \param[in] parent_ish
         * IndexSetHolder of the parent this MeshPart refers to.
         *
         * \warning This will in general change the implied topology of the MeshPart and meant to be used in mesh
         * preprocessing only.
         *
         */
        void deduct_topology(const ParentIndexSetHolderType& parent_ish)
        {
          if(_index_set_holder == nullptr)
            _index_set_holder = new IndexSetHolderType(_num_entities);

          Intern::IndexSetFiller<ShapeType::dimension>::fill_ish(
            *_index_set_holder, _target_set_holder, parent_ish);

          // build redundant index sets
          RedundantIndexSetBuilder<ShapeType>::compute(*_index_set_holder);

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
        typedef typename MeshPartType::AttributeSetContainer AttributeSetContainer;
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
        virtual void fill_attribute_sets(AttributeSetContainer& attribute_set_holder) = 0;

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

    /// \cond internal
    namespace Intern
    {
      template<typename ParentMesh_>
      struct TargetSetRefineParentWrapper;
    }
    /// \endcond

    /**
     * \brief StandardRefinery implementation for MeshPart
     *
     * \author Jordi Paul
     * \author Peter Zajac
     */
    template<typename ParentMesh_>
    class StandardRefinery< MeshPart<ParentMesh_> > :
      public Factory< MeshPart<ParentMesh_> >
    {
      public:
        /// mesh type
        typedef MeshPart<ParentMesh_> MeshType;
        /// shape type
        typedef typename MeshType::ShapeType ShapeType;
        /// Attribute set type
        typedef typename MeshType::AttributeSetType AttributeType;
        /// index set holder type
        typedef typename MeshType::AttributeSetContainer AttributeSetContainer;
        /// index set holder type
        typedef typename MeshType::IndexSetHolderType IndexSetHolderType;
        /// target set holder type
        typedef typename MeshType::TargetSetHolderType TargetSetHolderType;

        /// shape dimension
        static constexpr int shape_dim = ShapeType::dimension;

      protected:
        /// coarse mesh reference
        const MeshType& _coarse_meshpart;
        /// coarse parent mesh reference
        const ParentMesh_* _parent_mesh;
        /// coarse parent mesh-part reference
        const MeshType* _parent_meshpart;
        /// number of entities for coarse mesh
        Index _num_entities_coarse[shape_dim + 1];
        /// number of entities for fine mesh
        Index _num_entities_fine[shape_dim + 1];

      public:
        /**
         * \brief Constructor.
         *
         * Use this constructor if the meshpart, which is to be refined, is defined on the
         * root mesh (and not on another meshpart).
         *
         * \param[in] coarse_meshpart
         * A reference to the coarse mesh that is to be refined.
         *
         * \param[in] parent_mesh
         * A reference to the coarse parent mesh.
         */
        explicit StandardRefinery(const MeshType& coarse_meshpart, const ParentMesh_& parent_mesh) :
          _coarse_meshpart(coarse_meshpart),
          _parent_mesh(&parent_mesh),
          _parent_meshpart(nullptr)
        {
          // get number of entities in coarse mesh
          for(int i(0); i <= shape_dim; ++i)
          {
            _num_entities_fine[i] = _num_entities_coarse[i] = coarse_meshpart.get_num_entities(i);
          }

          // calculate number of entities in fine mesh
          Intern::EntityCountWrapper<Intern::StandardRefinementTraits, ShapeType>::query(_num_entities_fine);
        }

        /**
         * \brief Constructor.
         *
         * Use this constructor if the meshpart, which is to be refined, is defined on
         * another meshpart (and not on the root mesh).
         *
         * \param[in] coarse_meshpart
         * A reference to the coarse mesh that is to be refined.
         *
         * \param[in] parent_meshpart
         * A reference to the coarse parent mesh part.
         */
        explicit StandardRefinery(const MeshType& coarse_meshpart, const MeshPart<ParentMesh_>& parent_meshpart) :
          _coarse_meshpart(coarse_meshpart),
          _parent_mesh(nullptr),
          _parent_meshpart(&parent_meshpart)
        {
          // get number of entities in coarse mesh
          for(int i(0); i <= shape_dim; ++i)
          {
            _num_entities_fine[i] = _num_entities_coarse[i] = coarse_meshpart.get_num_entities(i);
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
        virtual Index get_num_entities(int dim) override
        {
          return _num_entities_fine[dim];
        }

        /**
         * \brief Fills attribute sets where applicable
         *
         * \param[in,out] attribute_container
         * Container for the attribute sets being added to
         *
         */
        virtual void fill_attribute_sets(AttributeSetContainer& attribute_container) override
        {
          // Attributes of shape dimension 0 only make sense if we have a mesh topology
          if(!_coarse_meshpart.has_topology())
            return;

          // Iterate over the attributes in the coarse mesh
          typename MeshType::AttributeSetConstIterator it(_coarse_meshpart.get_mesh_attributes().begin());
          typename MeshType::AttributeSetConstIterator jt(_coarse_meshpart.get_mesh_attributes().end());

          for(; it != jt; ++it)
          {
            // Create a new empty attribute of the desired size
            AttributeType* refined_attribute = new AttributeType(
              get_num_entities(0), it->second->get_dimension());

            // Refine the attribute in the coarse mesh and write the result to the new attribute
            Intern::StandardAttribRefineWrapper<ShapeType, AttributeType>
              ::refine(*refined_attribute, *(it->second), *_coarse_meshpart.get_topology());

            // Add the attribute to the corresponding set
            XASSERTM(attribute_container.insert(std::make_pair(it->first, refined_attribute)).second, "Error refining attribute " + it->first);
          }
        }

        /**
         * \brief Fills the index sets.
         *
         * \param[in,out] index_set_holder
         * The index set holder whose index sets are to be filled.
         */
        virtual void fill_index_sets(IndexSetHolderType*& index_set_holder) override
        {
          XASSERT(index_set_holder == nullptr);

          // Only create the topology for the refined mesh if the coarse mesh has a topology
          if(!_coarse_meshpart.has_topology())
            return;

          index_set_holder = new IndexSetHolderType(_num_entities_fine);

          // refine indices
          Intern::IndexRefineWrapper<ShapeType>::refine(
            *index_set_holder, _num_entities_coarse, *_coarse_meshpart.get_topology());
        }

        /**
         * \brief Fills the target sets.
         *
         * \param[in,out] target_set_holder
         * The target set holder whose target sets are to be filled.
         */
        virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
        {
          // We have to use a helper class here, because the target set refinement
          // depends heavily on the underlying mesh type. For this, we first need
          // to check whether our parent is the root mesh or another mesh-part.
          if(_parent_mesh != nullptr)
          {
            // the parent is the root mesh
            Intern::TargetSetRefineParentWrapper<ParentMesh_>::fill_target_sets(target_set_holder,
              _coarse_meshpart.get_target_set_holder(), _coarse_meshpart.get_topology(), *_parent_mesh);
          }
          else if(_parent_meshpart != nullptr)
          {
            // the parent is a mesh-part
            Intern::TargetSetRefineParentWrapper<ParentMesh_>::fill_target_sets(target_set_holder,
              _coarse_meshpart.get_target_set_holder(), _coarse_meshpart.get_topology(), *_parent_meshpart);
          }
          else
            XABORTM("no parent present");
        }

    }; // class StandardRefinery<MeshPart<...>,...>

    /// \cond internal
    namespace Intern
    {
      /**
       * \brief Specialisation for ConformalMesh
       */
      template<typename Shape_, int num_coords_, typename Coord_>
      struct TargetSetRefineParentWrapper<ConformalMesh<Shape_, num_coords_, Coord_>>
      {
        static constexpr int shape_dim = Shape_::dimension;

        // Note: This functions works for both mesh parents and mesh-part parents.
        template<typename TargetSetHolderType_, typename IndexSetHolder_, typename Parent_>
        static void fill_target_sets(
          TargetSetHolderType_& target_set_holder,
          const TargetSetHolderType_& coarse_target_set_holder,
          const IndexSetHolder_* coarse_ish,
          const Parent_& parent)
        {
          // get number of parent entities
          Index num_entities_parent[shape_dim+1];
          for(int i(0); i <= shape_dim; ++i)
            num_entities_parent[i] = parent.get_num_entities(i);

          if(coarse_ish != nullptr)
          {
            // The coarse mesh-part has a topology, so refine the target set to match.
            // The parent must have a topology, too, which is non-trivial if the
            // parent is a mesh-part itself.
            const auto* parent_topo = parent.get_topology();
            XASSERTM(parent_topo != nullptr, "mesh-part has topology, but parent doesn't");
            Intern::TargetRefineWrapper<Shape_>::refine(target_set_holder,
              num_entities_parent, coarse_target_set_holder, *coarse_ish, *parent_topo);
          }
          else
          {
            // The coarse mesh-part does not have a topology, stick to simple refinement
            Intern::SimpleTargetRefineWrapper<Shape_>::refine(
              target_set_holder, num_entities_parent, coarse_target_set_holder);
          }
        }
      }; // // struct TargetSetRefineParentWrapper<ConformalMesh<...>>

      /**
       * \brief Specialisation for StructuredMesh
       */
      template<int shape_dim_, int num_coords_, typename Coord_>
      struct TargetSetRefineParentWrapper<StructuredMesh<shape_dim_, num_coords_, Coord_>>
      {
        typedef Shape::Hypercube<shape_dim_> ShapeType;

        template<typename TargetSetHolderType_, typename IndexSetHolder_>
        static void fill_target_sets(
          TargetSetHolderType_& target_set_holder,
          const TargetSetHolderType_& coarse_target_set_holder,
          const IndexSetHolder_* coarse_ish,
          const StructuredMesh<shape_dim_, num_coords_, Coord_>& parent)
        {
          // get number of coarse/fine parent slices
          Index num_slices_c[shape_dim_], num_slices_f[shape_dim_];
          for(int i(0); i < shape_dim_; ++i)
          {
            num_slices_c[i] = parent.get_num_slices(i);
            num_slices_f[i] = Index(2) * num_slices_c[i];
          }

          if(coarse_ish != nullptr)
          {
            XABORTM("TargetSetRefineParentWrapper not implemented");
            //Intern::TargetRefineWrapper<ShapeType>::refine(
              //target_set_holder, num_entities_parent, coarse_target_set_holder, *coarse_ish, parent_ish);
          }
          else
          {
            Intern::StructuredTargetRefineWrapper<shape_dim_>::refine_simple(
              target_set_holder, coarse_target_set_holder, num_slices_c, num_slices_f);
          }
        }

        // Specialisation for MeshPart<StructuredMesh<...>>
        // Note: The topology of a MeshPart is always "conforming" (ie "unstructured"),
        // therefore the implementation is identical to the class template specialisation
        // for ConformalMesh above.
        template<typename TargetSetHolderType_, typename IndexSetHolder_>
        static void fill_target_sets(
          TargetSetHolderType_& target_set_holder,
          const TargetSetHolderType_& coarse_target_set_holder,
          const IndexSetHolder_* coarse_ish,
          const MeshPart<StructuredMesh<shape_dim_, num_coords_, Coord_>>& parent)
        {
          // get number of parent entities
          Index num_entities_parent[shape_dim_+1];
          for(int i(0); i <= shape_dim_; ++i)
            num_entities_parent[i] = parent.get_num_entities(i);

          if(coarse_ish != nullptr)
          {
            // The coarse mesh-part has a topology, so refine the target set to match.
            // The parent must have a topology, too, which is non-trivial if the
            // parent is a mesh-part itself.
            const auto* parent_topo = parent.get_topology();
            XASSERTM(parent_topo != nullptr, "mesh-part has topology, but parent doesn't");
            Intern::TargetRefineWrapper<ShapeType>::refine(target_set_holder,
              num_entities_parent, coarse_target_set_holder, *coarse_ish, *parent_topo);
          }
          else
          {
            // The coarse mesh-part does not have a topology, stick to simple refinement
            Intern::SimpleTargetRefineWrapper<ShapeType>::refine(
              target_set_holder, num_entities_parent, coarse_target_set_holder);
          }
        }
      }; // struct TargetSetRefineParentWrapper<StructuredMesh<...>>
    } // namespace Intern
    /// \endcond

#ifdef FEAT_EICKT
    extern template class MeshPart<ConformalMesh<Shape::Simplex<2>, 2, Real>>;
    extern template class MeshPart<ConformalMesh<Shape::Simplex<3>, 3, Real>>;
    extern template class MeshPart<ConformalMesh<Shape::Hypercube<2>, 2, Real>>;
    extern template class MeshPart<ConformalMesh<Shape::Hypercube<3>, 3, Real>>;

    extern template class StandardRefinery<MeshPart<ConformalMesh<Shape::Simplex<2>, 2, Real>>>;
    extern template class StandardRefinery<MeshPart<ConformalMesh<Shape::Simplex<3>, 3, Real>>>;
    extern template class StandardRefinery<MeshPart<ConformalMesh<Shape::Hypercube<2>, 2, Real>>>;
    extern template class StandardRefinery<MeshPart<ConformalMesh<Shape::Hypercube<3>, 3, Real>>>;
#endif // FEAT_EICKT
  } // namespace Geometry
} // namespace FEAT
#endif
