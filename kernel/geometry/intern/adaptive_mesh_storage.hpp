// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include "kernel/util/string.hpp"
#include "kernel/util/tiny_algebra.hpp"
#include <kernel/shape.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/slotmap.hpp>

#include <array>
#include <cstdint>
#include <utility>

namespace FEAT::Geometry::Intern
{

  /**
   * \brief Newtype wrapper for mesh layers
   *
   * This serves to keep callsites for AdaptiveMesh methods readable. The
   * AdaptiveMesh expands many of the ConformalMesh APIs with an additional
   * layer parameter. This newtype wrapper makes it easier to see at a glance
   * which parameter is the mesh layer and which is, for example, a dimension.
   *
   * \author Markus Muegge
   */
  struct Layer
  {
    /// Index of the layer
    Index idx;
  };

  /// Layer comparison operator
  inline bool operator==(Layer lhs, Layer rhs)
  {
    return lhs.idx == rhs.idx;
  }

  /// Layer comparison operator
  inline bool operator!=(Layer lhs, Layer rhs)
  {
    return !(lhs == rhs);
  }

  inline bool operator<(Layer lhs, Layer rhs)
  {
    return lhs.idx < rhs.idx;
  }

  inline bool operator>(Layer lhs, Layer rhs)
  {
    return lhs.idx > rhs.idx;
  }

  inline bool operator<=(Layer lhs, Layer rhs)
  {
    return lhs.idx <= rhs.idx;
  }

  inline bool operator>=(Layer lhs, Layer rhs)
  {
    return lhs.idx >= rhs.idx;
  }

  /**
   * \brief SlotMap key for use by the AdaptiveMeshStorage class
   *
   * See \c SlotMapKey for basic details. This key is additionally tagged with additional
   * information required by the mesh storage.
   *
   * \author Markus Muegge
   */
  template<int n>
  struct ElementKey
  {
    ElementKey() = default;
    ElementKey(std::uint32_t idx, std::uint32_t gen) : index(idx), generation(gen)
    {
    }

    ElementKey(const ElementKey& other) = default;
    ElementKey(ElementKey&& other) = default;

    ElementKey& operator=(const ElementKey& other) = default;
    ElementKey& operator=(ElementKey&& other) = default;

    /// Index into the SlotMaps' slots. See \c SlotMapKey for details.
    std::uint32_t index = 0;
    /// Generation of this key. See \c SlotMapKey for details.
    std::uint32_t generation = 0;

    /** \brief Element Layer
     *
     * This is required by the AdaptiveMeshStorage to store mesh entities per layer.
     * Because keys are not tagged with the SlotMap they belong to, the storage needs
     * to know what layer a key belongs to.
     * Unused by the SlotMap.
     */
    Layer layer = {0};
    /** \brief Element Kind
     *
     * Indicates whether the element this key points to is a permanent element,
     * which is part of all further layers, or a transient element.
     * Used by the AdaptiveMeshStorage to store permanent and transient elements
     * separately, which enables easier indexing and guarantees better caching behaviour.
     * Unused by the SlotMap.
     */
    bool is_permanent = false;
  };

  /// ElementKey comparison operators
  template<int dim_>
  inline bool operator==(const ElementKey<dim_>& lhs, const ElementKey<dim_>& rhs)
  {
    return lhs.index == rhs.index && lhs.generation == rhs.generation && lhs.layer == rhs.layer &&
           lhs.is_permanent == rhs.is_permanent;
  }

  /// ElementKey comparison operators
  template<int dim_>
  inline bool operator!=(const ElementKey<dim_>& lhs, const ElementKey<dim_>& rhs)
  {
    return !(lhs == rhs);
  }

  template<int dim_>
  inline std::ostream& operator<<(std::ostream& stream, const ElementKey<dim_> key)
  {
    stream
      << "ElementKey<" << stringify(dim_)
      << ">{ index = " << stringify(key.index)
      << ", generation = " << stringify(key.generation)
      << ", layer = Layer {" << stringify(key.layer.idx)
      << "}, is_permanent = " << stringify(key.is_permanent) << " }";
    return stream;
  }

  using VertexKey = ElementKey<0>;
  using EdgeKey = ElementKey<1>;
  using FaceKey = ElementKey<2>;
  using CellKey = ElementKey<3>;

  /**
   * \brief Orientation-aware reference to another mesh element
     *
   * \author Markus Muegge
   */
  template<int n>
  struct OrientedElement
  {
    OrientedElement() = default;
    OrientedElement(int o, ElementKey<n> k) : orientation(o), key(k)
    {
    }

    int orientation = 0;
    ElementKey<n> key = ElementKey<n>();
  };

  using OrientedEdge = OrientedElement<1>;
  using OrientedFace = OrientedElement<2>;

  template<int dim_>
  inline bool operator==(const OrientedElement<dim_>& lhs, const OrientedElement<dim_>& rhs)
  {
    return lhs.key == rhs.key && lhs.orientation == rhs.orientation;
  }

  template<int dim_>
  inline bool operator!=(const OrientedElement<dim_>& lhs, const OrientedElement<dim_>& rhs)
  {
    return !(lhs == rhs);
  }

  /**
   * \brief Children of a mesh element
   *
   * \tparam TemplateSet_ TemplateSet used for the mesh
   * \tparam Shape_ The shape of the parent entity
   * \tparam dim_ Current dimension of this data structure
   *
   * Note that the size of this struct is determined at compile time
   * and that this struct is stored inline with all other mesh element data.
   * Using a template set with a large amount of children for any template will
   * thus increase the memory usage of _all_ mesh elements.
   *
   * \author Markus Muegge
   */
  template<typename TemplateSet_, typename Shape_, int dim_ = Shape_::dimension>
  struct ElementChildren : public ElementChildren<TemplateSet_, Shape_, dim_ - 1>
  {
    /// Maximum number of children of dimension \c dim_
    static constexpr int max_children = TemplateSet_::template max_children<Shape_::dimension, dim_>();

    /// Children of dimension \c dim_
    std::array<Intern::ElementKey<dim_>, max_children> children;

    template<int query_dim_>
    std::array<ElementKey<query_dim_>, TemplateSet_::template max_children<Shape_::dimension, query_dim_>()>& by_dim()
    {
      static_assert(query_dim_ <= Shape_::dimension);

      return ElementChildren<TemplateSet_, Shape_, query_dim_>::children;
    }

    template<int query_dim_>
    const std::array<ElementKey<query_dim_>, TemplateSet_::template max_children<Shape_::dimension, query_dim_>()>& by_dim() const
    {
      static_assert(query_dim_ <= Shape_::dimension);

      return ElementChildren<TemplateSet_, Shape_, query_dim_>::children;
    }
  };

  /**
   * \brief Children of a mesh element
   *
   * \tparam TemplateSet_ TemplateSet used for the mesh
   * \tparam Shape_ The shape of the parent entity
   * \tparam dim_ Current dimension of this data structure
   *
   * Base case for recursive data structure
   *
   * Note that the size of this struct is determined at compile time
   * and that this struct is stored inline with all other mesh element data.
   * Using a template set with a large amount of children for any template will
   * thus increase the memory usage of _all_ mesh elements.
   *
   * \author Markus Muegge
   */
  template<typename TemplateSet_, typename Shape_>
  struct ElementChildren<TemplateSet_, Shape_, 0>
  {
    /// Maximum number of children of dimension 0
    static constexpr int max_children = TemplateSet_::template max_children<Shape_::dimension, 0>();

    /// Children of dimension 0
    std::array<Intern::ElementKey<0>, max_children> children;

    template<int query_dim_>
    std::array<ElementKey<query_dim_>, TemplateSet_::template max_children<Shape_::dimension, query_dim_>()>& by_dim()
    {
      static_assert(query_dim_ == 0);

      return children;
    }

    template<int query_dim_>
    const std::array<ElementKey<query_dim_>, TemplateSet_::template max_children<Shape_::dimension, query_dim_>()>& by_dim() const
    {
      static_assert(query_dim_ == 0);

      return children;
    }
  };

  /**
   * \brief Surrounding elements of a mesh element
   *
   * Stores references to the vertices/edges/faces that make up a mesh element.
   * References to edges and faces are stored together with their orientations.
   *
   * \author Markus Muegge
   */
  template<typename Shape_, int dim_ = Shape_::dimension - 1>
  struct ElementTopology : public ElementTopology<Shape_, dim_ - 1>
  {
    static_assert(dim_ >= 0);

    /// Number of surrounding elements at dimension \c dim_
    static constexpr int num_entities = Shape::FaceTraits<Shape_, dim_>::count;

    /// Surrounding entities at dimension \c dim_
    std::array<OrientedElement<dim_>, num_entities> entities;

    template<int key_dim_>
    OrientedElement<key_dim_> key_by_dim(Index idx) const
    {
      static_assert(key_dim_ < Shape_::dimension);
      ASSERT(idx < (ElementTopology<Shape_, key_dim_>::num_entities));

      return ElementTopology<Shape_, key_dim_>::entities[idx];
    }

    template<int query_dim_>
    std::array<OrientedElement<query_dim_>, Shape::FaceTraits<Shape_, query_dim_>::count>& by_dim()
    {
      static_assert(query_dim_ < Shape_::dimension);

      return ElementTopology<Shape_, query_dim_>::entities;
    }

    template<int query_dim_>
    const std::array<OrientedElement<query_dim_>, Shape::FaceTraits<Shape_, query_dim_>::count>& by_dim() const
    {
      static_assert(query_dim_ < Shape_::dimension);

      return ElementTopology<Shape_, query_dim_>::entities;
    }

    template<int query_dim_>
    constexpr int size()
    {
      static_assert(query_dim_ <= Shape_::dimension);

      return Shape::FaceTraits<Shape_, query_dim_>::count;
    }
  };

  template<typename Shape_>
  struct ElementTopology<Shape_, 0>
  {
    /// Number of surrounding elements at dimension \c dim_
    static constexpr int num_entities = Shape::FaceTraits<Shape_, 0>::count;

    /// Surrounding entities at dimension \c dim_
    std::array<OrientedElement<0>, num_entities> entities;

    template<int key_dim_>
    OrientedElement<key_dim_> key_by_dim(Index idx) const
    {
      static_assert(key_dim_ == 0);
      ASSERT(idx < num_entities);

      return entities[idx];
    }

    template<int query_dim_>
    std::array<OrientedElement<query_dim_>, Shape::FaceTraits<Shape_, query_dim_>::count>& by_dim()
    {
      static_assert(query_dim_ == 0);

      return entities;
    }

    template<int query_dim_>
    const std::array<OrientedElement<query_dim_>, Shape::FaceTraits<Shape_, query_dim_>::count>& by_dim() const
    {
      static_assert(query_dim_ == 0);

      return entities;
    }

    template<int query_dim_>
    constexpr int size()
    {
      static_assert(query_dim_ == 0);

      return Shape::FaceTraits<Shape_, query_dim_>::count;
    }
  };

  /**
   * \brief Mesh Element Tree Node
   *
   * Stores topology and children for a edge/face/cell.
   *
   * \author Markus Muegge
   */
  template<typename TemplateSet_, typename Shape_>
  struct AdaptiveElement
  {
    static_assert(Shape_::dimension > 0);

    using TemplateSet = TemplateSet_;

    using RefinementType = typename TemplateSet::template RefinementTypeByDim<Shape_::dimension>;

    AdaptiveElement(RefinementType t, Layer l, const ElementTopology<Shape_> topo) :
      type(t),
      layer(l),
      topology(std::move(topo))
    {
    }

    /// The elements refinement type
    RefinementType type;
    /// The elements mesh index
    Index index = 0;
    /// The elements layer
    Layer layer = {0};

    // OPTIMIZATION: We are allocating the space for the maximum number of
    // children of any template of the template set. We can instead store
    // children out-of-band in Util::SecondaryMaps with two advantages:
    // * we do not need to allocate space for children of zero type elements
    // * improved cache-behaviour when iterating over mesh elements, since most
    // algorithms do not need the parent-child relations
    /// Child element references
    ElementChildren<TemplateSet_, Shape_> children;
    /// Surrounding element references
    ElementTopology<Shape_> topology;

    static constexpr int num_neighbors = Shape::FaceTraits<Shape_, Shape_::dimension - 1>::count;
    std::array<Intern::ElementKey<Shape_::dimension>, num_neighbors> neighbors;
  };

  /**
   * \brief Vertex Info
   *
   * Stores information about a vertex
   *
   * \author Markus Muegge
   */
  template<typename VertexType>
  struct AdaptiveVertex
  {
    AdaptiveVertex(VertexType v, Layer l) : vertex(v), layer(l), last_changed(l)
    {
    }

    /// Vertex coordinates
    VertexType vertex;
    /// This elements index in the mesh layers it is a part of
    Index index = 0;
    /// This elements depth in the refinement tree
    Layer layer = {0};

    // NOTE: last_changed is (so far) exclusively used for the
    // BPXJacobiPreconditioner. The mesh should probably not directly know that this
    // property is needed. Maybe introduce something like the PMP property system
    // that allows dynamically adding data to any mesh entities.

    /// Depth the vertex was last used as a parent vertex
    Layer last_changed;
  };

  /**
   * \brief Storage class for a single layer of the adaptive mesh
   *
   * This class provides storage and utility methods to manage the elements of a single layer of the adaptive mesh.
   *
   * Uses self-inheritance to construct itself with support for all required dimensions.
   *
   * \tparam TemplateSet_ Template set for mesh refinement
   * \tparam VertexType_ Type of vertices to be stored
   * \tparam element_dim_ Dimension of this component of the mesh layer
   *
   * \author Markus Muegge
   */
  template<typename TemplateSet_, typename MeshShape_, typename VertexType_, int element_dim_ = MeshShape_::dimension>
  class MeshLayer : public MeshLayer<TemplateSet_, MeshShape_, VertexType_, element_dim_ - 1>
  {
  protected:
    using BaseClass = MeshLayer<TemplateSet_, MeshShape_, VertexType_, element_dim_>;

    /// Shape type of entities stored in this layer of dimension \c element_dim_
    using ShapeType = typename Shape::FaceTraits<MeshShape_, element_dim_>::ShapeType;

    /// Type of AdaptiveElement managed by this dimension of the MeshLayer
    using AdaptiveElementType = AdaptiveElement<TemplateSet_, ShapeType>;
    /// Key type for AdaptiveElementType
    using ElementKeyType = ElementKey<element_dim_>;

    /// Slotmap for permanent mesh elements
    Util::SlotMap<AdaptiveElementType, ElementKeyType> _permanent_entities;
    /// Slotmap for transient mesh elements
    Util::SlotMap<AdaptiveElementType, ElementKeyType> _transient_entities;

  public:
    /**
     * \brief Insert a permanent element into this layer
     *
     * \param element r-value reference to an adaptive element
     * \returns Key of the inserted element
     *
     * Note that the returned key is not yet tagged with the correct layer.
     * You need to do so yourself.
     */
    template<typename Shape_>
    ElementKey<Shape_::dimension> insert_permanent(const AdaptiveElement<TemplateSet_, Shape_>& element)
    {
      auto key = MeshLayer<TemplateSet_, MeshShape_, VertexType_, Shape_::dimension>::_permanent_entities.insert(element);
      key.is_permanent = true;
      return key;
    }

    /**
     * \brief Insert a transient element into this layer
     *
     * \param element r-value reference to an adaptive element
     * \returns Key of the inserted element
     *
     * Note that the returned key is not yet tagged with the correct layer.
     * You need to do so yourself.
     */
    template<typename Shape_>
    ElementKey<Shape_::dimension> insert_transient(const AdaptiveElement<TemplateSet_, Shape_>& element)
    {
      auto key = MeshLayer<TemplateSet_, MeshShape_, VertexType_, Shape_::dimension>::_transient_entities.insert(element);
      key.is_permanent = false;
      return key;
    }

    /**
     * \brief Erase an element from this layer
     *
     * \param key Key of the element to be erased
     */
    template<int dim_>
    void erase(ElementKey<dim_> key)
    {
      if constexpr(dim_ == 0)
      {
        MeshLayer<TemplateSet_, MeshShape_, VertexType_, 0>::_vertices.erase(key);
      }
      else
      {
        if(key.is_permanent)
        {
          MeshLayer<TemplateSet_, MeshShape_, VertexType_, dim_>::_permanent_entities.erase(key);
        }
        else
        {
          MeshLayer<TemplateSet_, MeshShape_, VertexType_, dim_>::_transient_entities.erase(key);
        }
      }
    }

    /**
     * \brief Get number of permanent entities of dimension dim_ in this layer
     *
     * \tparam dim_ Dimension of mesh elements to return number of
     * \returns Number of permanent mesh elements of dimension dim_
     */
    template<int dim_>
    Index num_permanent_entities() const
    {
      if constexpr(dim_ == 0)
      {
        return MeshLayer<TemplateSet_, MeshShape_, VertexType_, 0>::_vertices.size();
      }
      else
      {
        return MeshLayer<TemplateSet_, MeshShape_, VertexType_, dim_>::_permanent_entities.size();
      }
    }

    /**
     * \brief Get number of transient entities of dimension dim_ in this layer
     *
     * \tparam dim_ Dimension of mesh elements to return number of
     * \returns Number of transient mesh elements of dimension dim_
     */
    template<int dim_>
    Index num_transient_entities() const
    {
      if constexpr(dim_ == 0)
      {
        return 0;
      }
      else
      {
        return MeshLayer<TemplateSet_, MeshShape_, VertexType_, dim_>::_transient_entities.size();
      }
    }

    /// Permanent element SlotMap accessor
    template<int dim_>
    Util::SlotMap<AdaptiveElement<TemplateSet_, typename Shape::FaceTraits<MeshShape_, dim_>::ShapeType>, ElementKey<dim_>>&
    permanent_entities()
    {
      return MeshLayer<TemplateSet_, MeshShape_, VertexType_, dim_>::_permanent_entities;
    }

    /// Permanent element SlotMap accessor
    template<int dim_>
    const Util::SlotMap<AdaptiveElement<TemplateSet_, typename Shape::FaceTraits<MeshShape_, dim_>::ShapeType>, ElementKey<dim_>>&
    permanent_entities() const
    {
      return MeshLayer<TemplateSet_, MeshShape_, VertexType_, dim_>::_permanent_entities;
    }

    /// Transient element SlotMap accessor
    template<int dim_>
    Util::SlotMap<AdaptiveElement<TemplateSet_, typename Shape::FaceTraits<MeshShape_, dim_>::ShapeType>, ElementKey<dim_>>&
    transient_entities()
    {
      return MeshLayer<TemplateSet_, MeshShape_, VertexType_, dim_>::_transient_entities;
    }

    /// Transient element SlotMap accessor
    template<int dim_>
    const Util::SlotMap<AdaptiveElement<TemplateSet_, typename Shape::FaceTraits<MeshShape_, dim_>::ShapeType>, ElementKey<dim_>>&
    transient_entities() const
    {
      return MeshLayer<TemplateSet_, MeshShape_, VertexType_, dim_>::_transient_entities;
    }

    /**
     * \brief Returns number of bytes used by this layer
     */
    std::size_t bytes() const
    {
      return _permanent_entities.bytes() + _transient_entities.bytes() + BaseClass::bytes();
    }

    /**
     *  \brief Retrieve the AdaptiveElement belonging to the given key
     *
     *  \tparam dim_ Dimension of element to retrieve
     *  \param key Key of element to retrieve
     *  \returns A reference to the AdaptiveElement pointed to by the key
     */
    template<int dim_>
    AdaptiveElement<TemplateSet_, typename Shape::FaceTraits<MeshShape_, dim_>::ShapeType>& get(ElementKey<dim_> key)
    {
      static_assert(dim_ >= 1);
      if(key.is_permanent)
      {
        return MeshLayer<TemplateSet_, MeshShape_, VertexType_, dim_>::_permanent_entities[key];
      }
      else
      {
        return MeshLayer<TemplateSet_, MeshShape_, VertexType_, dim_>::_transient_entities[key];
      }
    }

    /**
     *  \brief Retrieve the AdaptiveElement belonging to the given key
     *
     *  \tparam dim_ Dimension of element to retrieve
     *  \param key Key of element to retrieve
     *  \returns A reference to the AdaptiveElement pointed to by the key
     */
    template<int dim_>
    const AdaptiveElement<TemplateSet_, typename Shape::FaceTraits<MeshShape_, dim_>::ShapeType>& get(ElementKey<dim_> key) const
    {
      static_assert(dim_ >= 1);
      if(key.is_permanent)
      {
        return MeshLayer<TemplateSet_, MeshShape_, VertexType_, dim_>::_permanent_entities[key];
      }
      else
      {
        return MeshLayer<TemplateSet_, MeshShape_, VertexType_, dim_>::_transient_entities[key];
      }
    }
  };

  /**
   * \brief Storage class for a single layer of the adaptive mesh
   *
   * This class provides storage and utility methods to manage the elements of a single layer of the adaptive mesh.
   *
   * This specialization is the base case of the self-inheritance.
   *
   * \tparam TemplateSet_ Template set for mesh refinement
   * \tparam VertexType_ Type of vertices to be stored
   *
   * \author Markus Muegge
   */
  template<typename TemplateSet_, typename MeshShape_, typename VertexType_>
  class MeshLayer<TemplateSet_, MeshShape_, VertexType_, 0>
  {
  protected:
    /// SlotMap for this layers vertices
    Util::SlotMap<AdaptiveVertex<VertexType_>, ElementKey<0>> _vertices;

  public:
    /**
     * \brief Insert a vertex into this layer
     *
     * \param vertex r-value reference to the AdaptiveVertex to be inserted
     */
    ElementKey<0> insert_vertex(const AdaptiveVertex<VertexType_>& vertex)
    {
      return _vertices.insert(vertex);
    }

    /**
     * \brief Erase a vertex from this layer
     *
     * \param key Key pointing to the vertex to be erased
     */
    void erase(ElementKey<0> key)
    {
      _vertices.erase(key);
    }

    /**
     * \brief Returns the number of vertices in this layer
     */
    Index num_vertices() const
    {
      return _vertices.size();
    }

    /// Vertex SlotMap accessor
    Util::SlotMap<AdaptiveVertex<VertexType_>, ElementKey<0>>& vertices()
    {
      return _vertices;
    }

    /// Vertex SlotMap accessor
    const Util::SlotMap<AdaptiveVertex<VertexType_>, ElementKey<0>>& vertices() const
    {
      return _vertices;
    }

    /**
     * \brief Returns the number of bytes used by this layer
     */
    std::size_t bytes() const
    {
      return _vertices.bytes();
    }

    /**
     * \brief Retrieve a reference to the AdaptiveVertex pointed to by \c key
     */
    AdaptiveVertex<VertexType_>& get_vertex(ElementKey<0> key)
    {
      return _vertices[key];
    }

    /**
     * \brief Retrieve a reference to the AdaptiveVertex pointed to by \c key
     */
    const AdaptiveVertex<VertexType_>& get_vertex(ElementKey<0> key) const
    {
      return _vertices[key];
    }
  };

  /** \brief Storage class for AdaptiveMeshes
   *
   * Handles storing refinement trees for the AdaptiveMesh. Mesh elements are
   * stored separated by layer.
   *
   * Allows accessing elements via SlotMap keys and by index. To allow access
   * via indices, the storage must first be indexed via
   * AdaptiveMeshStorage::reindex(). After creation and after any modification
   * of the tree structure the storage is unindexed. If the storage is
   * unindexed, retrieving an element via its index will abort the program.
   * Ensure you call AdaptiveMeshStorage::reindex after any modification before
   * retrieving elements
   *
   * \tparam MeshShape_ Shape type of the mesh to be stored \tparam
   * TemplateSet_ Template set for mesh refinement \tparam VertexType_ Type of
   * mesh vertices to be stored
   *
   * \author Markus Muegge
   */
  template<typename MeshShape_, typename TemplateSet_, typename VertexType_>
  class AdaptiveMeshStorage
  {
  public:
    /// Type of vertices in this mesh
    using VertexType = VertexType_;

    using CoordType = typename VertexType_::ValueType;

    /// Template set to be used for refinement
    using TemplateSet = TemplateSet_;

    /// Accessor for element reference types by dimension
    template<int dim_>
    using ElementRefByDim = ElementKey<dim_>;

    /// Accessor for oriented element types by dimension
    template<int dim_>
    using OrientedElementRefByDim = OrientedElement<dim_>;

    /// Accessor for element topologies by dimension
    template<int dim_>
    using ElementTopologyByDim = ElementTopology<typename Shape::FaceTraits<MeshShape_, dim_>::ShapeType>;

    /// Accessor for element children by dimension
    template<int dim_>
    using ElementChildrenByDim = ElementChildren<TemplateSet_, typename Shape::FaceTraits<MeshShape_, dim_>::ShapeType>;

    /// Accessor for adaptive elements by dimension
    template<int dim_>
    using AdaptiveElementByDim = AdaptiveElement<TemplateSet_, typename Shape::FaceTraits<MeshShape_, dim_>::ShapeType>;

    /// Adaptive vertex type
    using AdaptiveVertexType = AdaptiveVertex<VertexType>;

    /**
     * Traversal orders for iterating over refinement trees
     */
    enum class IterationOrder
    {
      PreOrder,
      PostOrder
    };

  private:
    /// Type of a mesh layer for this mesh
    using MeshLayerType = MeshLayer<TemplateSet, MeshShape_, VertexType>;

    static constexpr int num_coords = VertexType::n;

    /// Sentinel key for non-existing neighbors, e.g. for boundary elements.
    static const constexpr ElementKey<MeshShape_::dimension> neighbor_sentinel = ElementKey<MeshShape_::dimension>();

    /// List of MeshLayers for this mesh
    std::vector<MeshLayerType> _layers;

    /// Indicates whether this mesh has been indexed since the last modification
    bool _is_indexed = false;

    /// Indicates whether neighbors have been calculated for this mesh
    bool _has_neighbors = false;

  public:
    /**
     * \brief Get refinement type of element
     *
     * \param[in] key
     * The key of the element whose type should be retrieved
     *
     * \return The refinement type of the element
     */
    template<int dim>
    typename TemplateSet::template RefinementTypeByDim<dim> type(ElementKey<dim> key)
    {
      return (*this)[key].type;
    }

    /**
     * \brief Erases an element and all its children from the storage
     *
     * \param[in] key Key of the root element to be erased
     *
     * \returns An array indicating how many elemets of each dimension were erased
     *
     * \attention Unindexes the storage!
     */
    template<int dim_>
    std::array<Index, 4> erase(ElementKey<dim_> key)
    {
      ASSERTM(key.layer.idx < _layers.size(), "Trying to erase adaptive mesh element of non-existing layer");

      // Erasing an entity produces a gap in that layers indices. The mesh is
      // thus no longer correctly indexed.
      _is_indexed = false;
      // Erasing an entity (most-likely) removed a neighbor of an element. Neighbors need to be re-determined.
      _has_neighbors = false;

      std::array<Index, 4> num_erased = {};
      if constexpr(dim_ > 0)
      {
        // NOTE: Different mesh layers are stored separately. Deleting a child
        // does thus not invalidate the references to the current entity or its
        // children.
        auto& entity = (*this)[key];
        num_erased =
          _erase_children<typename Shape::FaceTraits<MeshShape_, dim_>::ShapeType, dim_>(entity.type, entity.children);
      }

      // Actually delete current entity
      _layers[key.layer.idx].erase(key);
      num_erased[dim_] += 1;

      return num_erased;
    }

    /**
     * \brief Insert a vertex into the storage
     *
     * \param[in] vertex r-value reference to the vertex to be inserted
     *
     * \attention Unindexes the storage!
     */
    VertexKey insert(const AdaptiveVertexType& vertex)
    {
      // Adding a new vertex to a layer shifts the indices of all vertices on
      // further layers. These vertices stored indices are thus wrong, and need
      // to be recomputed.
      _is_indexed = false;

      // Add additional mesh layers, if required
      Layer layer = vertex.layer;
      _add_mesh_layers(layer.idx);

      // Add vertex
      VertexKey key = _layers[layer.idx].insert_vertex(vertex);
      key.layer = layer;
      return key;
    }

    /**
     * \brief Insert a vertex into the storage
     *
     * \param[in] vertex Vertex coordinate vector
     * \param[in] layer Layer to insert the vertex into
     *
     * \attention Unindexes the storage!
     */
    VertexKey insert(VertexType vertex, Layer layer)
    {
      return insert(AdaptiveVertex(vertex, layer));
    }

    /**
     * \brief Insert an adaptive element into the storage
     *
     * \param[in] element r-value reference to the element to be inserted
     *
     * \attention Unindexes the storage!
     */
    template<typename Shape_>
    ElementKey<Shape_::dimension> insert(const AdaptiveElement<TemplateSet, Shape_>& element)
    {
      // Adding an entity might cause the indices of entities on further layers
      // to be shifted. These entites stored indices are thus wrong and need to
      // be recomputed.
      _is_indexed = false;
      // Neighbors for this new entity need to be determined.
      _has_neighbors = false;

      // Add additional mesh layers, if required
      Layer layer = element.layer;
      _add_mesh_layers(layer.idx);

      ElementKey<Shape_::dimension> key;
      bool is_permanent = element.type.is_zero_refinement();

      if(is_permanent)
      {
        key = _layers[layer.idx].insert_permanent(element);
      }
      else
      {
        key = _layers[layer.idx].insert_transient(element);
      }

      // Add metadata to key
      key.is_permanent = is_permanent;
      key.layer = layer;

      return key;
    }

    /**
     * \brief Return the number of layers in this mesh
     *
     * \note Layer 0 contains a duplicate of those parts of the underlying
     * mesh, that will be refined further. It is counted as a mesh layer by
     * this method, but you might not consider it a mesh layer for your use
     * case. In that case layer 1 is the first refined layer.
     */
    Index num_layers() const
    {
      return _layers.size();
    }

    /**
     * Traverses the sub-tree anchored at the element pointed to by \c root_key
     * in the given order.
     *
     * \tparam VisitorType_ Object with
     * * VisitorType_::operator()(VertexKey)
     * * VisitorType_::operator()(EdgeKey)
     * * VisitorType_::operator()(FaceKey)
     * * VisitorType_::operator()(CellKey)
     * methods.
     *
     * \param[in] root_key Key pointing to root of subtree
     * \param[in] visitor Visitor object
     * \param[in] max_depth Maximum depth (inclusive)
     * \param[in] order Traversal order
     */
    template<typename VisitorType_, int dim_>
    void walk_subtree(
      ElementKey<dim_> root_key,
      VisitorType_& visitor,
      Layer max_depth,
      IterationOrder order = IterationOrder::PreOrder) const
    {
      auto& root = (*this)[root_key];

      // Abort if max_depth is reached
      if(root.layer > max_depth)
      {
        return;
      }

      // Visit Root
      if(order == IterationOrder::PreOrder)
      {
        visitor(root_key);
      }

      if constexpr(dim_ > 0)
      {
        _walk_subtree_children<typename Shape::FaceTraits<MeshShape_, dim_>::ShapeType>(root, visitor, max_depth);
      }

      // Visit Root
      if(order == IterationOrder::PostOrder)
      {
        visitor(root_key);
      }
    }

    /**
     * Traverses all refinement trees in the given order
     *
     * \tparam VisitorType_ Object with
     * * VisitorType_::operator()(VertexKey)
     * * VisitorType_::operator()(EdgeKey)
     * * VisitorType_::operator()(FaceKey)
     * * VisitorType_::operator()(CellKey)
     * methods.
     *
     * \param[in] root_key Key pointing to root of subtree
     * \param[in] visitor Visitor object
     * \param[in] max_depth Maximum depth (inclusive)
     * \param[in] order Traversal order
     */
    template<typename VisitorType_>
    void walk_tree(VisitorType_& visitor, Layer max_depth, IterationOrder order = IterationOrder::PreOrder) const
    {
      // There is nothing to iterate over if the mesh has no layers
      if(_layers.size() == 0)
      {
        return;
      }

      const MeshLayerType& root_layer = _layers[0];

      if constexpr(MeshShape_::dimension >= 0)
      {
        for(auto entry : root_layer.vertices())
        {
          walk_subtree(entry.key, visitor, max_depth, order);
        }
      }

      if constexpr(MeshShape_::dimension >= 1)
      {
        for(auto entry : root_layer.template permanent_entities<1>())
        {
          walk_subtree(entry.key, visitor, max_depth, order);
        }
        for(auto entry : root_layer.template transient_entities<1>())
        {
          walk_subtree(entry.key, visitor, max_depth, order);
        }
      }

      if constexpr(MeshShape_::dimension >= 2)
      {
        for(auto entry : root_layer.template permanent_entities<2>())
        {
          walk_subtree(entry.key, visitor, max_depth, order);
        }
        for(auto entry : root_layer.template transient_entities<2>())
        {
          walk_subtree(entry.key, visitor, max_depth, order);
        }
      }

      if constexpr(MeshShape_::dimension >= 3)
      {
        for(auto entry : root_layer.template permanent_entities<3>())
        {
          walk_subtree(entry.key, visitor, max_depth, order);
        }
        for(auto entry : root_layer.template transient_entities<3>())
        {
          walk_subtree(entry.key, visitor, max_depth, order);
        }
      }
    }

    /**
     * \brief Returns the number of mesh elements of the given dimension on the given layer
     *
     * \param[in] layer Layer to count elements of
     * \param[in] dim_ Dimension to query
     *
     * \note O(num_layers) complexity
     */
    Index num_entities(Layer layer, int dim) const
    {
      if(dim == 0)
      {
        return _num_entities<0>(layer);
      }
      if(dim == 1)
      {
        return _num_entities<1>(layer);
      }
      if(dim == 2)
      {
        return _num_entities<2>(layer);
      }
      if(dim == 3)
      {
        return _num_entities<3>(layer);
      }
      XABORTM("unsupported dimension in MeshStorage::num_entities");
    }

    /**
     * \brief Returns the total number of all mesh elements in the storage
     *
     * \tparam dim_ Dimension to query
     *
     * \note O(num_layers) complexity
     */
    template<int dim_>
    Index num_total_entities() const
    {
      Index result = 0;
      for(auto& layer : _layers)
      {
        result += layer.template num_permanent_entities<dim_>();
        result += layer.template num_transient_entities<dim_>();
      }

      return result;
    }

    /**
     * \brief Determines neighboring elements for all mesh elements
     */
    void fill_neighbors()
    {
      static constexpr int shape_dim = MeshShape_::dimension;
      static constexpr int num_facets = Shape::FaceTraits<MeshShape_, shape_dim - 1>::count;

      // OPTIMIZATION: This is the ConformalMesh::fill_neighbors() algorithm
      // ported to adaptive mesh keys. This information can be retrieved more
      // cheaply from the refinement templates themselves. Add this capability to
      // the template query (maybe optionally, with a boolean indicating support
      // on the template set), support it in the AdaptiveMesh code, and don't
      // invalidate neighbors when adding elements here if the template provide
      // neighbor information.
      // This method can then be left as a fallback.

      if(_has_neighbors)
      {
        // Neighbors have already been determined since the last adaptation of
        // the mesh. No work to be done.
        return;
      }

      // In the following we are calling mesh entities of dimension shape_dim - 1 facets.
      // We are using the fact that each facet is shared by at most 2 elements.
      // We are goint to first determine which elements are adjacent to each
      // facet, by iterating over the elements and adding their keys to the
      // facet_neighbors of each of their facets. Then we can iterate over the
      // elements a second time and determine their neighbors as the element
      // registered for the facet, that is not itself.

      // Map of facet neighbors
      Util::SecondaryMap<std::array<ElementKey<shape_dim>, 2>, ElementKey<shape_dim - 1>>
        facet_neighbors;

      // Register element with their adjacent facets. Fills the facet_neighbors map.
      auto determine_facet_neighbors = [&](Layer layer, auto slotmap_entry)
      {
        // Iterate over all facets
        for(int facet_idx = 0; facet_idx < num_facets; facet_idx++)
        {
          // Retrieve key of facet
          auto key = slotmap_entry.value.topology.template key_by_dim<shape_dim - 1>(facet_idx).key;

          // Init secondary map, if this is the first time we see this facet
          if(!facet_neighbors.contains_key(key))
          {
            facet_neighbors.insert(key, {neighbor_sentinel, neighbor_sentinel});
          }

          auto& neighbors = facet_neighbors[key];

          // HACK: We have to recreate the "full" element key here.
          // The keys returned by the SlotMap iterators do not contain the
          // layer and permanency information we have added to the usual
          // slotmap data.
          // The proper way to this is to wrap the iterators in the MeshLayer
          // class and add this information there. As this is currently the
          // only place where we directly iterate over all elements of a layer
          // via the iterators, I am leaving that boilerplate out for now.
          ElementKey<shape_dim> element_key = slotmap_entry.key;
          element_key.layer = layer;
          element_key.is_permanent = slotmap_entry.value.type.is_zero_refinement();

          if(neighbors[0] == neighbor_sentinel)
          {
            neighbors[0] = element_key;
          }
          else if(neighbors[1] == neighbor_sentinel)
          {
            neighbors[1] = element_key;
          }
          else
          {
            XABORTM("Facet has more than two neighbors!");
          }
        }
      };

      // Determine neighbors by looking at registered facet neighbors.
      auto determine_neighbors = [&](Layer layer, auto slotmap_entry)
      {
        auto& element = slotmap_entry.value;

        // Iterate over all facets
        for(int facet_idx = 0; facet_idx < num_facets; facet_idx++)
        {
          // Retrieve key of facet
          auto key = slotmap_entry.value.topology.template key_by_dim<shape_dim - 1>(facet_idx).key;

          // HACK: We have to recreate the "full" element key here.
          // The keys returned by the SlotMap iterators do not contain the
          // layer and permanency information we have added to the usual
          // slotmap data.
          // The proper way to this is to wrap the iterators in the MeshLayer
          // class and add this information there. As this is currently the
          // only place where we directly iterate over all elements of a layer
          // via the iterators, I am leaving that boilerplate out for now.
          ElementKey<shape_dim> element_key = slotmap_entry.key;
          element_key.layer = layer;
          element_key.is_permanent = slotmap_entry.value.type.is_zero_refinement();

          // Neighbor is the entity that is not the current entity
          auto& candidates = facet_neighbors[key];
          if(candidates[0] == element_key)
          {
            element.neighbors[facet_idx] = candidates[1];
          }
          else if(candidates[1] == element_key)
          {
            element.neighbors[facet_idx] = candidates[0];
          }
        }
      };

      // Find elements adjacent to all facets
      for(Index i(0); i < _layers.size(); i++)
      {
        for(auto entry : _layers[i].template permanent_entities<shape_dim>())
        {
          determine_facet_neighbors(Layer{i}, entry);
        }
        for(auto entry : _layers[i].template transient_entities<shape_dim>())
        {
          determine_facet_neighbors(Layer{i}, entry);
        }
      }

      for(Index i(0); i < _layers.size(); i++)
      {
        for(auto entry : _layers[i].template permanent_entities<shape_dim>())
        {
          determine_neighbors(Layer{i}, entry);
        }
        for(auto entry : _layers[i].template transient_entities<shape_dim>())
        {
          determine_neighbors(Layer{i}, entry);
        }
      }

      _has_neighbors = true;
    }

    /**
     * \brief Applies a "proper rigid" transformation onto the mesh vertices.
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
     *
     * \warning Affects all layers of the mesh
     */
    void transform(const VertexType& origin, const VertexType& angles, const VertexType& offset)
    {
      // create rotation matrix
      Tiny::Matrix<CoordType, num_coords, num_coords> rot;
      _aux_rot_mat(rot, angles);

      // transform all vertices
      for(auto& layer : _layers)
      {
        for(auto entry : layer.vertices())
        {
          entry.value.set_mat_vec_mul(rot, entry.value - origin) += offset;
        }
      }
    }

    /**
     * \brief Re-indexes the storage
     *
     * Assigns indices to all mesh elements, such that each mesh layer can be
     * accessed via a continous range of numeric indices starting at 0.
     */
    void reindex()
    {
      // NOTE: The index-schema works as follows. On each layer of the mesh,
      // first the permanent entities get assigned indices, then the transient
      // entities get assigned indices. Indices are assigned in order,
      // starting on each layer with the index one beyond the last index of the
      // previous layers permanent indices. In other words, the transient
      // indices of each layer get reused by the next layer.  This is fine,
      // because transient entities are replaced with their children on the
      // next layer and thus no longer require an index.
      //
      // Visually the indeces are thus staggered like this:
      //
      // L0: |-- Perm. Indices 0 --|---- Trans. Indices 0 ----|
      // L1: |-- Perm. Indices 0 --|-- Perm. Indices 1 --|- Trans. Indices 1 -|
      // L2: |-- Perm. Indices 0 --|-- Perm. Indices 1 --|---- Perm. Indices 2 ----|-- Trans. Indices 2 --|

      Index next_vertex_index = 0;
      Index next_edge_index = 0;
      Index next_face_index = 0;
      Index next_cell_index = 0;

      for(auto& layer : _layers)
      {
        if constexpr(MeshShape_::dimension >= 0)
        {
          // There are only permanent vertices. We can just index them without
          // having to worry about transient indices.
          for(auto entry : layer.vertices())
          {
            entry.value.index = next_vertex_index++;
          }
        }

        if constexpr(MeshShape_::dimension >= 1)
        {
          _index_layer<1>(layer, next_edge_index);
        }

        if constexpr(MeshShape_::dimension >= 2)
        {
          _index_layer<2>(layer, next_face_index);
        }

        if constexpr(MeshShape_::dimension >= 3)
        {
          _index_layer<3>(layer, next_cell_index);
        }
      }

      // Mesh is freshly indexed, by definition
      _is_indexed = true;
    }

    ///////////////////////////////////
    // Index-based Accessors
    ///////////////////////////////////

    /**
     * \brief Retrieves a vertex of the mesh by its layer and index
     *
     * \param[in] layer The layer to retrieve a vertex of
     * \param[in] idx The numeric index of the vertex
     *
     * \returns A reference to the vertex
     *
     * \attention Requires the storage to be indexed!
     */
    const AdaptiveVertexType& get_vertex_by_index(Layer layer, Index idx) const
    {
      ASSERTM(layer.idx < _layers.size(), "Trying to retrieve vertex of non-existing mesh layer!");

      // NOTE: We require the storage to be indexed, despite not making use of
      // the index member of any mesh element. This is so that, given a
      // ElementKey k pointing to a vertex v with index i on layer l, it is
      // impossible that get_index(k) != i or get_vertex_by_index(l, i) != v
      ASSERT(_is_indexed);

      // OPTIMIZATION: This is a O(num_layers) search for the right index
      // range.  We can optimize this slightly by precomputing the index range
      // for each layer, eliminating the subtraction of the current idx.  We
      // can then optimize this to O(log num_layers) by building a search tree.
      // And finally we can optimize the average case by caching the
      // index-range of the last access in a thread-local variable, allowing
      // linear iteration in O(1).

      // Search through permanent edges of previous and current layer
      for(Index layer_idx = 0; layer_idx <= layer.idx; layer_idx++)
      {
        const MeshLayerType& current_layer = _layers[layer_idx];

        if(idx < current_layer.num_vertices())
        {
          return current_layer.vertices().data()[idx];
        }
        idx -= current_layer.num_vertices();
      }

      XABORTM("Tried retrieving non-existent vertex index from AdaptiveMeshStorage!");
    }

    /**
     * \brief Retrieves a vertex of the mesh by its layer and index
     *
     * \param[in] layer The layer to retrieve a vertex of
     * \param[in] idx The numeric index of the vertex
     *
     * \returns A reference to the vertex
     *
     * \attention Requires the storage to be indexed!
     */
    AdaptiveVertexType& get_vertex_by_index(Layer layer, Index idx)
    {
      // Reuse const accessor by casting this to const and then removing const from the return value
      // SAFE: Casting to const is always fine
      // SAFE: Removing const is fine because this method was called on a non const value to begin with
      return const_cast<AdaptiveVertexType&>(std::as_const(*this).get_vertex_by_index(layer, idx));
    }

    /**
     * \brief Retrieves an element of the mesh by its layer and index
     *
     * \tparam dim_ Dimension of the element to retrieve
     * \param[in] layer The layer to retrieve the element from
     * \param[in] idx The numeric index of the element
     *
     * \returns A reference to the element
     *
     * \attention Requires the storage to be indexed!
     */
    template<int dim_>
    const AdaptiveElementByDim<dim_>& get_by_index(Layer layer, Index idx) const
    {
      ASSERTM(layer.idx < _layers.size(), "Trying to retrieve entity of non-existing mesh layer!");

      // NOTE: We require the storage to be indexed, despite not making use of
      // the index member of any mesh element. This is so that, given a
      // ElementKey k pointing to a mesh element e with index i on layer l, it is
      // impossible that get_index(k) != i or get_by_index(l, i) != e
      ASSERT(_is_indexed);

      // Search through permanent edges of previous and current layer
      for(Index layer_idx = 0; layer_idx <= layer.idx; layer_idx++)
      {
        const MeshLayerType& current_layer = _layers[layer_idx];

        if(idx < current_layer.template num_permanent_entities<dim_>())
        {
          return current_layer.template permanent_entities<dim_>().data()[idx];
        }
        idx -= current_layer.template num_permanent_entities<dim_>();
      }

      const MeshLayerType& current_layer = _layers[layer.idx];
      if(idx < current_layer.template num_transient_entities<dim_>())
      {
        return current_layer.template transient_entities<dim_>().data()[idx];
      }

      XABORTM("Trying to retrieve non-exisiting element from AdaptiveMeshStorage");
    }

    /**
     * \brief Retrieves an element of the mesh by its layer and index
     *
     * \tparam dim_ Dimension of the element to retrieve
     * \param[in] layer The layer to retrieve the element from
     * \param[in] idx The numeric index of the element
     *
     * \returns A reference to the element
     *
     * \attention Requires the storage to be indexed!
     */
    template<int dim_>
    AdaptiveElementByDim<dim_>& get_by_index(Layer layer, Index idx)
    {
      // Reuse const accessor by casting this to const and then removing const from the return value
      // SAFE: Casting to const is always fine
      // SAFE: Removing const is fine because this method was called on a non const value to begin with
      return const_cast<AdaptiveElementByDim<dim_>&>(std::as_const(*this).template get_by_index<dim_>(layer, idx));
    }

    /**
     * \brief Returns the index of the n-th neighbor of the given element
     *
     * \returns The index of the neighbor, if one exists, ~Index(0) otherwise.
     */
    Index get_neighbor(Layer layer, Index element_idx, Index neighbor_idx) const
    {
      static constexpr int shape_dim = MeshShape_::dimension;
      ASSERT(_has_neighbors);
      ASSERT(_is_indexed);

      auto& element = get_by_index<shape_dim>(layer, element_idx);
      if(element.neighbors[neighbor_idx] != neighbor_sentinel)
      {
        return get_index(element.neighbors[neighbor_idx]);
      }
      return ~Index(0);
    }

    ///////////////////////////////////
    // Key-based Access Operators
    ///////////////////////////////////

    /**
     * \brief Get the numeric mesh index of the element pointed to by \c key
     *
     * \tparam dim_ Dimension to query
     * \param[in] key Key pointing to mesh element whose index should be retrieved
     *
     * \attention Requires the storage to be indexed!
     */
    template<int dim_>
    Index get_index(ElementKey<dim_> key) const
    {
      ASSERT(_is_indexed);
      return (*this)[key].index;
    }

    /**
     * \brief Get the numeric mesh index of the element pointed to by \c key
     *
     * \tparam dim_ Dimension to query
     * \param[in] key Key pointing to mesh element whose index should be retrieved
     *
     * \attention Requires the storage to be indexed!
     */
    template<int dim_>
    Index get_index(OrientedElement<dim_> oriented_key) const
    {
      ASSERT(_is_indexed);
      return (*this)[oriented_key.key].index;
    }

    /// Retrieve vertex by key
    AdaptiveVertexType& operator[](VertexKey key)
    {
      ASSERT(key.layer.idx < _layers.size());
      return _layers[key.layer.idx].get_vertex(key);
    }

    /// Retrieve vertex by key
    const AdaptiveVertexType& operator[](VertexKey key) const
    {
      ASSERT(key.layer.idx < _layers.size());
      return _layers[key.layer.idx].get_vertex(key);
    }

    /// Retrieve element by key
    template<int dim_>
    AdaptiveElementByDim<dim_>& operator[](ElementKey<dim_> key)
    {
      ASSERT(key.layer.idx < _layers.size());
      return _layers[key.layer.idx].get(key);
    }

    /// Retrieve element by key
    template<int dim_>
    const AdaptiveElementByDim<dim_>& operator[](ElementKey<dim_> key) const
    {
      ASSERT(key.layer.idx < _layers.size());
      return _layers[key.layer.idx].get(key);
    }

    AdaptiveVertexType& operator[](OrientedElement<0> ref)
    {
      return (*this)[ref.key];
    }

    const AdaptiveVertexType& operator[](OrientedElement<0> ref) const
    {
      return (*this)[ref.key];
    }

    /// Retrieve element by oriented reference
    template<int dim_>
    AdaptiveElementByDim<dim_>& operator[](OrientedElement<dim_> ref)
    {
      return (*this)[ref.key];
    }

    template<int dim_>
    const AdaptiveElementByDim<dim_>& operator[](OrientedElement<dim_> ref) const
    {
      return (*this)[ref.key];
    }

    /**
     * \brief Returns the number of bytes used by this class
     */
    std::size_t bytes() const
    {
      std::size_t total = 0;
      for(auto& layer : _layers)
      {
        total += layer.bytes();
      }
      return total;
    }

  private:
    /**
     * \brief Returns the number of mesh elements of the given dimension on the given layer
     *
     * \tparam dim_ Dimension to query
     * \param[in] layer Layer to count elements of
     *
     * \note O(num_layers) complexity
     */
    template<int dim_>
    Index _num_entities(Layer layer) const
    {
      if constexpr(dim_ > MeshShape_::dimension)
      {
        return 0;
      }
      else
      {
        ASSERT(layer.idx < _layers.size());

        Index result = 0;
        for(Index layer_idx = 0; layer_idx <= layer.idx; layer_idx++)
        {
          result += _layers[layer_idx].template num_permanent_entities<dim_>();
        }
        result += _layers[layer.idx].template num_transient_entities<dim_>();

        return result;
      }
    }

    template<typename Shape_, int dim_ = Shape_::dimension>
    std::array<Index, 4> _erase_children(
        typename TemplateSet::template RefinementTypeByDim<Shape_::dimension> type,
        ElementChildren<TemplateSet, Shape_>& children)
    {
      if constexpr(dim_ >= 0)
      {
        auto result = _erase_children<Shape_, dim_ - 1>(type, children);

        auto& array = children.template by_dim<dim_>();
        for(Index i(0); i < TemplateSet::template num_children<Shape_::dimension, dim_>(type); i++)
        {
          std::array<Index, 4> num_erased = erase(array[i]);

          result[0] += num_erased[0];
          result[1] += num_erased[1];
          result[2] += num_erased[2];
          result[3] += num_erased[3];
        }

        return result;
      }
      else
      {
        return {};
      }
    }

    template<typename Shape_, int dim_ = Shape_::dimension, typename VisitorType>
    void _walk_subtree_children(const AdaptiveElement<TemplateSet, Shape_>& root, VisitorType& visitor, Layer max_depth) const
    {
      if constexpr(dim_ >= 0)
      {
        _walk_subtree_children<Shape_, dim_ - 1>(root, visitor, max_depth);

        const auto& children = root.children.template by_dim<dim_>();
        for(Index i = 0; i < TemplateSet::template num_children<Shape_::dimension, dim_>(root.type); i++)
        {
          walk_subtree(children[i], visitor, max_depth);
        }
      }
    }

    /**
     * \brief Adds additional mesh layers, such that
     *
     * \code_layers[new_top_layer]
     *
     * exists.
     */
    void _add_mesh_layers(Index new_top_layer)
    {
      while(_layers.size() <= new_top_layer)
      {
        _layers.push_back(MeshLayerType{});
      }
    }

    /**
     * \brief Assigns indexes to all elemens of a layer
     *
     * \tparam dim_ Dimension of elements to assign indices to
     * \param[in] layer Layer whose elements will be assigned indices
     * \param[inout] Initial next permanent index, will contain first permanent
     * index for next layer after call finishes
     */
    template<int dim_>
    void _index_layer(MeshLayerType& layer, Index& next_permanent_index)
    {
      // First assign indices to all permanent mesh elements
      for(auto entry : layer.template permanent_entities<dim_>())
      {
        entry.value.index = next_permanent_index++;
      }

      // Creat a copy of the current next index
      Index next_transient_index = next_permanent_index;

      // Assign indices to all transient mesh elements, using the copied index
      // This causes these indices to be reused on the next layer, where the
      // transient elements have been replaced with their children and no longer
      // exist
      for(auto entry : layer.template transient_entities<dim_>())
      {
        entry.value.index = next_transient_index++;
      }
    }

    /**
     * \brief Helper for transform method
     */
    template<int sm_, int sn_, int sv_>
    static void _aux_rot_mat(
      Tiny::Matrix<CoordType, 1, 1, sm_, sn_>& r,
      const Tiny::Vector<CoordType, 1, sv_>& /*unsused*/)
    {
      r.set_identity();
    }

    /**
     * \brief Helper for transform method
     */
    template<int sm_, int sn_, int sv_>
    static void _aux_rot_mat(
      Tiny::Matrix<CoordType, 2, 2, sm_, sn_>& r,
      const Tiny::Vector<CoordType, 2, sv_>& a)
    {
      r.set_rotation_2d(a[0]);
    }

    /**
     * \brief Helper for transform method
     */
    template<int sm_, int sn_, int sv_>
    static void _aux_rot_mat(
      Tiny::Matrix<CoordType, 3, 3, sm_, sn_>& r,
      const Tiny::Vector<CoordType, 3, sv_>& a)
    {
      r.set_rotation_3d(a[0], a[1], a[2]);
    }
  };
} // namespace FEAT::Geometry::Intern
