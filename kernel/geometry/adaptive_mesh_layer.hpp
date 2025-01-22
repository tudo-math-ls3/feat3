// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include "kernel/geometry/mesh_permutation.hpp"
#include "kernel/shape.hpp"
#include "kernel/util/string.hpp"
#include <kernel/adjacency/permutation.hpp>
#include <kernel/base_header.hpp>
#include <kernel/geometry/adaptive_mesh.hpp>
#include <kernel/util/tiny_algebra.hpp>

#include <memory>
#include <optional>
#include <utility>

namespace FEAT::Geometry
{

  /**
   * \brief VertexSet interface for adaptive meshes
   *
   * This class template is a compatability layer between the AdaptiveMesh
   * class and the rest of FEAT. It presents a normal VertexSet interface for a
   * single layer of an AdaptiveMesh. This allows using the AdaptiveMeshLayer
   * like any other mesh class.
   *
   * \author Markus Muegge
   */
  template<typename AdaptiveMeshType_>
  class AdaptiveVertexSet
  {
  public:
    /// Type of underlying adaptive mesh
    using AdaptiveMeshType = AdaptiveMeshType_;

    /// Number of coordinates per vertex
    static constexpr int num_coords = AdaptiveMeshType::world_dim;

    /// Type of single vertex coordinate
    using CoordType = typename AdaptiveMeshType::CoordType;

    /// Vertex type
    using VertexType = typename AdaptiveMeshType::VertexType;

  private:
    /// Underlying adaptive mesh
    std::shared_ptr<AdaptiveMeshType> _mesh;

    // Mesh layer we refer to
    Layer _layer;

    // Vertex permutation
    // Stored, because we can not apply it to the vertex set directly, as the
    // underlying memory is shared between all AdaptiveVertexSets.
    // Stored as an optional to avoid allocating a permutation if no permutation
    // is required by the user.
    std::optional<Adjacency::Permutation> _perm;

  public:
    /**
     * \brief Constructor
     *
     * \param[in] mesh Shared pointer to underlying mesh
     * \param[in] layer Layer this vertex set is to refer to
     */
    AdaptiveVertexSet(std::shared_ptr<AdaptiveMeshType> mesh, Layer layer) : _mesh(std::move(mesh)), _layer(layer)
    {
    }

    /**
     * \brief Returns a clone of this vertex set
     *
     * \warning A AdaptiveVertexSet is a view into a underlying AdaptiveMesh.
     * The returned clone will thus not be independent of this
     * AdaptiveVertexSet.
     */
    AdaptiveVertexSet clone() const
    {
      return AdaptiveVertexSet(_mesh, _layer);
    }

    std::size_t bytes() const
    {
      // We are just a view into data stored elsewhere.
      return 0;
    }

    /**
     * \brief Returns the number of coordinates per vertex.
     */
    int get_num_coords() const
    {
      return num_coords;
    }

    /**
     * \brief Returns the number of vertices in the vertex set.
     */
    Index get_num_vertices() const
    {
      return _mesh->get_num_entities(_layer, 0);
    }

    /**
     * \brief Returns a reference to a vertex.
     *
     * \param[in] idx The index of the vertex to be returned.
     *
     * \returns A reference to the vertex.
     */
    VertexType& operator[](Index idx)
    {
      if(_perm)
      {
        return _mesh->vertex(_layer, _perm.value().map(idx));
      }
      return _mesh->vertex(_layer, idx);
    }

    /** \copydoc operator[]() */
    const VertexType& operator[](Index idx) const
    {
      if(_perm)
      {
        return _mesh->vertex(_layer, _perm.value().map(idx));
      }
      return _mesh->vertex(_layer, idx);
    }

    // NOTE: Unsupported on purpose. Transform underlying mesh if required.
    //void transform(const VertexType& origin, const VertexType& angles, const VertexType& offset)

    /**
     * \brief Permutes this vertex set
     *
     * \param[in] perm The permutation to apply
     * \param[in] invert Specified whether to apply the forward or inverse permutation
     *
     * \warning Does not permute the underlying memory. The permutation is applied by mapping indices one each future
     * access.
     *
     * \warning Does not apply to any other AdaptiveVertexSet, even if the other AdaptiveVertexSet points to the same
     * AdaptiveMesh.
     */
    void permute(const Adjacency::Permutation& perm, bool invert = false)
    {
      if(!_perm)
      {
        _perm = std::optional{invert ? perm.inverse() : perm.clone()};
      }
      else
      {
        if(invert)
        {
          _perm.value().concat(perm.inverse());
        }
        else
        {
          _perm.value().concat(perm);
        }
      }
    }

    static String name()
    {
      return "AdaptiveVertexSet<...>";
    }
  };

  /**
   * \brief IndexTuple interface for adaptive meshes
   *
   * This class template is a compatability layer between the AdaptiveMesh
   * class and the rest of FEAT. It presents a normal IndexTuple interface for a
   * single layer of an AdaptiveMesh. This allows using the AdaptiveMeshLayer
   * like any other mesh class.
   *
   * \author Markus Muegge
   */
  template<typename AdaptiveMeshType_, int cell_dim_, int face_dim_>
  class AdaptiveIndexTuple
  {
  public:
    /// Type of underlying adaptive mesh
    using AdaptiveMeshType = AdaptiveMeshType_;

  private:
    /// Pointer to underlying mesh
    std::shared_ptr<AdaptiveMeshType> _mesh;
    /// Layer this tuple refers to
    Layer _layer;
    /// Index of entity this tuple refers to
    Index _entity_idx;
    /// Permutation of this tuple
    std::optional<Adjacency::Permutation> _inv_perm;

  public:
    // Wrapper Code

    /// Number of indices per tuple
    static constexpr int num_indices = Shape::FaceTraits<Shape::Hypercube<cell_dim_>, face_dim_>::count;

    /**
     * \brief Constructor
     *
     * \param[in] mesh Underlying mesh
     * \param[in] layer Layer this tuple refers to
     * \param[in] entity_idx Entity this tuple refers to
     */
    AdaptiveIndexTuple(std::shared_ptr<AdaptiveMeshType> mesh, Layer layer, Index entity_idx) :
      _mesh(std::move(mesh)),
      _layer(layer),
      _entity_idx(entity_idx),
      _inv_perm(std::nullopt)
    {
    }

    AdaptiveIndexTuple(const AdaptiveIndexTuple& other) = default;
    AdaptiveIndexTuple(AdaptiveIndexTuple&& other) = default;

    AdaptiveIndexTuple& operator=(const AdaptiveIndexTuple& other) = default;
    AdaptiveIndexTuple& operator=(AdaptiveIndexTuple&& other) = default;

    /// access operator
    Index operator[](int idx) const
    {
      ASSERT(idx >= 0);
      ASSERT(idx < num_indices);
      Index result = _mesh->template get_face_index<cell_dim_, face_dim_>(_layer, _entity_idx, idx);
      if(_inv_perm)
      {
        return _inv_perm.value().map(result);
      }
      return result;
    }

    /**
     * \brief Permute this index tuple
     *
     * \warning Does not permute the underlying memory. The permutation is
     * applied by mapping indices one each future
     * access.
     *
     * \warning Does not apply to any other AdaptiveIndexTuple, even if the
     * other AdaptiveIndexTuple points to the same
     * mesh entity.
     */
    void permute_map(const Adjacency::Permutation& inv_perm)
    {
      if(_inv_perm)
      {
        _inv_perm.value().concat(inv_perm);
      }
      else
      {
        _inv_perm = inv_perm.clone();
      }
    }
  };

  /**
   * \brief IndexSet interface for adaptive meshes
   *
   * This class template is a compatability layer between the AdaptiveMesh
   * class and the rest of FEAT. It presents a normal IndexSet interface for a
   * single layer of an AdaptiveMesh. This allows using the AdaptiveMeshLayer
   * like any other mesh class.
   *
   * \author Markus Muegge
   */
  template<typename AdaptiveMeshType_, int cell_dim_, int face_dim_>
  class AdaptiveIndexSet
  {
  public:
    /// Type of underlying mesh
    using AdaptiveMeshType = AdaptiveMeshType_;

    /// Number of indices per entry
    static constexpr int num_indices = Shape::FaceTraits<Shape::Hypercube<cell_dim_>, face_dim_>::count;

    /// Adjactor for AdaptiveIndexSet
    class AdaptiveIndexSetAdjactor
    {
      // Adaptive mesh reference
      std::shared_ptr<AdaptiveMeshType> _mesh;
      Layer _layer;

      // Iterator state
      Index _domain_node;
      Index _next;

    public:
      AdaptiveIndexSetAdjactor(std::shared_ptr<AdaptiveMeshType> mesh, Layer layer, Index domain_node, Index next) :
        _mesh(std::move(mesh)),
        _layer(layer),
        _domain_node(domain_node),
        _next(next)
      {
      }

      bool operator==(const AdaptiveIndexSetAdjactor& other) const
      {
        return _domain_node == other._domain_node && _next == other._next && _layer == other._layer &&
               _mesh == other._mesh;
      }

      bool operator!=(const AdaptiveIndexSetAdjactor& other) const
      {
        return !(*this == other);
      }

      AdaptiveIndexSetAdjactor& operator++()
      {
        _next++;
        return *this;
      }

      Index operator*() const
      {
        return _mesh->template get_face_index<cell_dim_, face_dim_>(_layer, _domain_node, _next);
      }
    };

    /// Adjactor type
    using ImageIterator = AdaptiveIndexSetAdjactor;

    /// Tuple type
    using IndexTupleType = AdaptiveIndexTuple<AdaptiveMeshType, cell_dim_, face_dim_>;

  private:
    // Mesh reference
    std::shared_ptr<AdaptiveMeshType> _mesh;
    Layer _layer;

    std::optional<Adjacency::Permutation> _permutation;
    std::optional<Adjacency::Permutation> _inverse_face_permutation;

  public:
    /**
     * \brief Constructor
     *
     * \param[in] mesh Underlying mesh
     * \param[in] Layer Layer this index set refers to
     */
    AdaptiveIndexSet(std::shared_ptr<AdaptiveMeshType> mesh, Layer layer) : _mesh(std::move(mesh)), _layer(layer)
    {
    }

    AdaptiveIndexSet(const AdaptiveIndexSet& other) = default;
    AdaptiveIndexSet(AdaptiveIndexSet&& other) = default;

    AdaptiveIndexSet& operator=(const AdaptiveIndexSet& other) = default;
    AdaptiveIndexSet& operator=(AdaptiveIndexSet&& other) = default;

    /**
     * \brief Returns a clone of this index set
     *
     * \warning A AdaptiveIndexSet is a view into a underlying AdaptiveMesh.
     * The returned clone will thus not be independent of this
     * AdaptiveIndexSet.
     */
    AdaptiveIndexSet clone() const
    {
      return AdaptiveIndexSet(_mesh, _layer);
    }

    /**
     * \brief Returns the size of dynamically allocated memory in bytes
     */
    std::size_t bytes() const
    {
      return 0;
    }

    /// access operator
    Index operator()(Index i, int j)
    {
      if(_permutation)
      {
        return _mesh->template get_face_index<cell_dim_, face_dim_>(
          _layer,
          _permutation.value().map(i),
          _inverse_face_permutation.value().map(j)
        );
      }
      return _mesh->template get_face_index<cell_dim_, face_dim_>(_layer, i, j);
    }

    /// access operator
    Index operator()(Index i, int j) const
    {
      if(_permutation)
      {
        return _mesh->template get_face_index<cell_dim_, face_dim_>(
          _layer,
          _permutation.value().map(i),
          _inverse_face_permutation.value().map(j)
        );
      }
      return _mesh->template get_face_index<cell_dim_, face_dim_>(_layer, i, j);
    }

    /**
     * \brief Returns an adaptive index tuple
     *
     * \param[in] i The index of the tuples entity
     *
     * \warning This method need to copy the shared_ptr to the AdaptiveMesh on each call. Use operator() if possible.
     */
    IndexTupleType operator[](Index i)
    {
      ASSERT(i < get_num_entities());
      if(_permutation)
      {
        auto tuple = IndexTupleType(_mesh, _layer, _permutation.value().map(i));
        tuple.permute_map(_inverse_face_permutation.value());
        return tuple;
      }
      return IndexTupleType(_mesh, _layer, i);
    }

    /** \copydoc operator[]() */
    IndexTupleType operator[](Index i) const
    {
      ASSERT(i < get_num_entities());
      if(_permutation)
      {
        auto tuple = IndexTupleType(_mesh, _layer, _permutation.value().map(i));
        tuple.permute_map(_inverse_face_permutation.value());
        return tuple;
      }
      return IndexTupleType(_mesh, _layer, i);
    }

    /**
     * \brief Returns the number of entities.
     */
    Index get_num_entities() const
    {
      return _mesh->get_num_entities(_layer, cell_dim_);
    }

    /**
     * \brief Returns the number of indices per entity
     */
    int get_num_indices() const
    {
      return num_indices;
    }

    /**
     * \brief Returns the maximum index of any face returned by this index set
     */
    Index get_index_bound() const
    {
      return _mesh->get_num_entities(_layer, face_dim_);
    }

    // Unsupported.
    // IndexSet API depends on pointer arithmetic that we can not support.
    // Rewrite IndexSet API first, if this is required somewhere
    //AdaptiveIndexTuple* get_indices()
    //const AdaptiveIndexTuple* get_indices() const

    void permute(const Adjacency::Permutation& perm, const Adjacency::Permutation& inv_perm_face)
    {
      // Permute index tuples
      if(_permutation)
      {
        _permutation.value().concat(perm);
      }
      else
      {
        _permutation = perm.clone();
      }

      // Permute indices by the inverse face permutation
      if(_inverse_face_permutation)
      {
        _inverse_face_permutation.value().concat(inv_perm_face);
      }
      else
      {
        _inverse_face_permutation = inv_perm_face.clone();
      }
    }

    static String name()
    {
      return "AdaptiveIndexSet<" + stringify(cell_dim_) + ", " + stringify(face_dim_) + ">";
    }

    // Adjactor Interface

    Index get_num_nodes_domain() const
    {
      return get_num_entities();
    }

    Index get_num_nodes_image() const
    {
      return _mesh->get_num_entities(_layer, face_dim_);
    }

    AdaptiveIndexSetAdjactor image_begin(Index domain_node) const
    {
      return AdaptiveIndexSetAdjactor(_mesh, _layer, domain_node, 0);
    }

    AdaptiveIndexSetAdjactor image_end(Index domain_node) const
    {
      return AdaptiveIndexSetAdjactor(_mesh, _layer, domain_node, num_indices);
    }
  };

  template<
    typename AdaptiveMeshType_,
    typename Shape_,
    int face_dim_ = Shape_::dimension - 1>
  class AdaptiveIndexSetWrapper :
    public AdaptiveIndexSetWrapper<AdaptiveMeshType_, Shape_, face_dim_ - 1>
  {
    static_assert(face_dim_ < Shape_::dimension, "invalid face dimension");
    static_assert(face_dim_ > 0, "invalid face dimension");

  public:
    using AdaptiveMeshType = AdaptiveMeshType_;

    static constexpr int cell_dim = Shape_::dimension;

    using BaseClass =
      AdaptiveIndexSetWrapper<AdaptiveMeshType_, Shape_, face_dim_ - 1>;

    using IndexSetType =
      AdaptiveIndexSet<AdaptiveMeshType_, Shape_::dimension, face_dim_>;


  protected:
    template<int face_dim__>
    using AdaptiveIndexSetByFaceDim =
      AdaptiveIndexSet<AdaptiveMeshType, cell_dim, face_dim__>;

    std::shared_ptr<AdaptiveMeshType> _mesh;
    Layer _layer;

    IndexSetType _index_set;

  public:
    AdaptiveIndexSetWrapper(std::shared_ptr<AdaptiveMeshType> mesh, Layer layer) :
      BaseClass(mesh, layer),
      _mesh(mesh),
      _layer(layer),
      _index_set(mesh, layer)
    {
    }

    AdaptiveIndexSetWrapper(const AdaptiveIndexSetWrapper& other) = default;
    AdaptiveIndexSetWrapper(AdaptiveIndexSetWrapper&& other) = default;

    AdaptiveIndexSetWrapper& operator=(const AdaptiveIndexSetWrapper& other) = default;
    AdaptiveIndexSetWrapper& operator=(AdaptiveIndexSetWrapper&& other) = default;

    void clone(const AdaptiveIndexSetWrapper& other)
    {
      BaseClass::clone(other);
      _index_set = other._index_set;
    }

    AdaptiveIndexSetWrapper clone() const
    {
      return AdaptiveIndexSetWrapper(_mesh, _layer);
    }

    template<int face_dim__>
    AdaptiveIndexSetByFaceDim<face_dim__>& get_index_set()
    {
      static_assert(face_dim__ >= 0, "invalid face dimension");
      static_assert(face_dim__ < Shape_::dimension, "invalid face dimension");

      return AdaptiveIndexSetWrapper<AdaptiveMeshType, Shape_, face_dim__>::_index_set;
    }

    template<int face_dim__>
    const AdaptiveIndexSetByFaceDim<face_dim__>& get_index_set() const
    {
      static_assert(face_dim__ >= 0, "invalid face dimension");
      static_assert(face_dim__ < Shape_::dimension, "invalid face dimension");

      return AdaptiveIndexSetWrapper<AdaptiveMeshType, Shape_, face_dim__>::_index_set;
    }

    template<std::size_t np_>
    void permute(const Adjacency::Permutation& shape_perm,
      const std::array<Adjacency::Permutation, np_>& inv_perm)
    {
      BaseClass::permute(shape_perm, inv_perm);
      _index_set.permute(shape_perm, inv_perm.at(face_dim_));
    }

    static String name()
    {
      return "AdaptiveIndexSetWrapper<" + Shape_::name() + ", " + stringify(face_dim_) + ">";
    }

    std::size_t bytes()
    {
      return BaseClass::bytes() + _index_set.bytes();
    }
  };

  template<typename AdaptiveMeshType_, typename Shape_>
  class AdaptiveIndexSetWrapper<AdaptiveMeshType_, Shape_, 0>
  {
    static_assert(Shape_::dimension > 0, "invalid shape dimension");

  public:
    using AdaptiveMeshType = AdaptiveMeshType_;

    static constexpr int cell_dim = Shape_::dimension;

    using IndexSetType =
      AdaptiveIndexSet<AdaptiveMeshType_, Shape_::dimension, 0>;


  protected:
    template<int face_dim__>
    using AdaptiveIndexSetByFaceDim =
      AdaptiveIndexSet<AdaptiveMeshType, cell_dim, face_dim__>;

    std::shared_ptr<AdaptiveMeshType> _mesh;
    Layer _layer;

    IndexSetType _index_set;

  public:
    AdaptiveIndexSetWrapper(std::shared_ptr<AdaptiveMeshType> mesh, Layer layer) :
      _mesh(mesh),
      _layer(layer),
      _index_set(mesh, layer)
    {
    }

    AdaptiveIndexSetWrapper(const AdaptiveIndexSetWrapper& other) = default;
    AdaptiveIndexSetWrapper(AdaptiveIndexSetWrapper&& other) = default;

    AdaptiveIndexSetWrapper& operator=(const AdaptiveIndexSetWrapper& other) = default;
    AdaptiveIndexSetWrapper& operator=(AdaptiveIndexSetWrapper&& other) = default;

    void clone(const AdaptiveIndexSetWrapper& other)
    {
      _index_set = other._index_set;
    }

    AdaptiveIndexSetWrapper clone() const
    {
      return AdaptiveIndexSetWrapper(_mesh, _layer);
    }

    template<int face_dim__>
    AdaptiveIndexSetByFaceDim<face_dim__>& get_index_set()
    {
      static_assert(face_dim__ == 0, "invalid face dimension");

      return AdaptiveIndexSetWrapper<AdaptiveMeshType, Shape_, face_dim__>::_index_set;
    }

    template<int face_dim__>
    const AdaptiveIndexSetByFaceDim<face_dim__>& get_index_set() const
    {
      static_assert(face_dim__ == 0, "invalid face dimension");

      return AdaptiveIndexSetWrapper<AdaptiveMeshType, Shape_, face_dim__>::_index_set;
    }

    template<std::size_t np_>
    void permute(const Adjacency::Permutation& shape_perm,
      const std::array<Adjacency::Permutation, np_>& inv_perm)
    {
      _index_set.permute(shape_perm, inv_perm.at(0));
    }

    static String name()
    {
      return "AdaptiveIndexSetWrapper<" + Shape_::name() + ",0>";
    }

    std::size_t bytes() const
    {
      return _index_set.bytes();
    }
  };

  /**
   * \brief IndexSetHolder interface for adaptive meshes
   *
   * This class template is a compatability layer between the AdaptiveMesh
   * class and the rest of FEAT. It presents a normal IndexSetHolder interface for a
   * single layer of an AdaptiveMesh. This allows using the AdaptiveMeshLayer
   * like any other mesh class.
   *
   * \author Markus Muegge
   */
  template<typename AdaptiveMeshType_, typename Shape_>
  class AdaptiveIndexSetHolder :
    public AdaptiveIndexSetHolder<
      AdaptiveMeshType_,
      typename Shape::FaceTraits<Shape_, Shape_::dimension - 1>::ShapeType>
  {
    /// Type of underlying mesh
    using AdaptiveMeshType = AdaptiveMeshType_;

    using BaseClass = AdaptiveIndexSetHolder<
      AdaptiveMeshType_,
      typename Shape::FaceTraits<Shape_, Shape_::dimension - 1>::ShapeType>;

    using IndexSetWrapperType = AdaptiveIndexSetWrapper<AdaptiveMeshType, Shape_>;

  protected:
    template<int shape_dim_>
    using AdaptiveIndexSetWrapperByShapeDim =
      AdaptiveIndexSetWrapper<
        AdaptiveMeshType,
        typename Shape::FaceTraits<Shape_, shape_dim_>::ShapeType>;

    std::shared_ptr<AdaptiveMeshType> _mesh;
    Layer _layer;

    IndexSetWrapperType _index_set_wrapper;

  public:
    AdaptiveIndexSetHolder(std::shared_ptr<AdaptiveMeshType> mesh, Layer layer) :
      BaseClass(mesh, layer),
      _mesh(mesh),
      _layer(layer),
      _index_set_wrapper(mesh, layer)
    {
    }

    AdaptiveIndexSetHolder(const AdaptiveIndexSetHolder& other) = default;
    AdaptiveIndexSetHolder(AdaptiveIndexSetHolder&& other) = default;

    AdaptiveIndexSetHolder& operator=(const AdaptiveIndexSetHolder& other) = default;
    AdaptiveIndexSetHolder& operator=(AdaptiveIndexSetHolder&& other) = default;

    void clone(const AdaptiveIndexSetHolder& other)
    {
      BaseClass::clone(other);
      _index_set_wrapper.clone(other._index_set_wrapper);
    }

    AdaptiveIndexSetHolder clone() const
    {
      return AdaptiveIndexSetHolder(_mesh, _layer);
    }

    template<int shape_dim_>
    AdaptiveIndexSetWrapperByShapeDim<shape_dim_>& get_index_set_wrapper()
    {
      static_assert(shape_dim_ > 0, "invalid shape dimension");
      static_assert(shape_dim_ <= Shape_::dimension, "invalid shape dimension");
      typedef typename Shape::FaceTraits<Shape_, shape_dim_>::ShapeType ShapeType;
      return AdaptiveIndexSetHolder<AdaptiveMeshType, ShapeType>::_index_set_wrapper;
    }

    template<int shape_dim_>
    const AdaptiveIndexSetWrapperByShapeDim<shape_dim_>& get_index_set_wrapper() const
    {
      static_assert(shape_dim_ > 0, "invalid shape dimension");
      static_assert(shape_dim_ <= Shape_::dimension, "invalid shape dimension");
      typedef typename Shape::FaceTraits<Shape_, shape_dim_>::ShapeType ShapeType;
      return AdaptiveIndexSetHolder<AdaptiveMeshType, ShapeType>::_index_set_wrapper;
    }

    template<int cell_dim_, int face_dim_>
    AdaptiveIndexSet<AdaptiveMeshType, cell_dim_, face_dim_>& get_index_set()
    {
      static_assert(cell_dim_ <= Shape_::dimension, "invalid cell dimension");
      static_assert(face_dim_ < cell_dim_, "invalid face/cell dimension");
      static_assert(face_dim_ >= 0, "invalid face dimension");
      return get_index_set_wrapper<cell_dim_>().template get_index_set<face_dim_>();
    }

    template<int cell_dim_, int face_dim_>
    const AdaptiveIndexSet<AdaptiveMeshType, cell_dim_, face_dim_>& get_index_set() const
    {
      static_assert(cell_dim_ <= Shape_::dimension, "invalid cell dimension");
      static_assert(face_dim_ < cell_dim_, "invalid face/cell dimension");
      static_assert(face_dim_ >= 0, "invalid face dimension");
      return get_index_set_wrapper<cell_dim_>().template get_index_set<face_dim_>();
    }

    template<std::size_t np_>
    void permute(const std::array<Adjacency::Permutation, np_>& perms,
      const std::array<Adjacency::Permutation, np_>& inv_perms)
    {
      BaseClass::permute(perms, inv_perms);
      _index_set_wrapper.permute(perms.at(Shape_::dimension), inv_perms);
    }

    static String name()
    {
      return "AdaptiveIndexSetHolder<" + Shape_::name() + ">";
    }

    std::size_t bytes() const
    {
      return BaseClass::bytes() + _index_set_wrapper.bytes();
    }
  };

  template<typename AdaptiveMeshType_>
  class AdaptiveIndexSetHolder<AdaptiveMeshType_, Shape::Vertex>
  {
  public:
    AdaptiveIndexSetHolder() = default;

    AdaptiveIndexSetHolder(std::shared_ptr<AdaptiveMeshType_> /*mesh*/, Layer /*layer*/){
    }

    AdaptiveIndexSetHolder(const AdaptiveIndexSetHolder& other) = default;
    AdaptiveIndexSetHolder(AdaptiveIndexSetHolder&& other) = default;

    AdaptiveIndexSetHolder& operator=(const AdaptiveIndexSetHolder& other) = default;
    AdaptiveIndexSetHolder& operator=(AdaptiveIndexSetHolder&& other) = default;

    void clone(const AdaptiveIndexSetHolder& /*unused*/)
    {
    }

    AdaptiveIndexSetHolder clone() const
    {
      return AdaptiveIndexSetHolder();
    }

    template<std::size_t np_>
    void permute(const std::array<Adjacency::Permutation, np_>&,
      const std::array<Adjacency::Permutation, np_>&)
    {
    }

    static String name()
    {
      return "AdaptiveIndexSetHolder<Vertex>";
    }

    std::size_t bytes() const
    {
      return std::size_t(0);
    }
  };

  /**
   * \brief Gives access to neighbor information on a mesh layer
   */
  template<typename AdaptiveMeshType_>
  class AdaptiveNeighbors
  {
  protected:
    using AdaptiveMeshType = AdaptiveMeshType_;

    std::shared_ptr<AdaptiveMeshType> _mesh;
    Layer _layer;

  public:

    /// Number of indices per entry
    static constexpr int num_indices =
      Shape::FaceTraits<typename AdaptiveMeshType::ShapeType, AdaptiveMeshType_::shape_dim>::count;

    /// Adjactor for AdaptiveIndexSet
    class AdaptiveNeighborsAdjactor
    {
      // Adaptive mesh reference
      std::shared_ptr<AdaptiveMeshType> _mesh;
      Layer _layer;

      // Iterator state
      Index _domain_node;
      Index _next;

    public:
      AdaptiveNeighborsAdjactor(std::shared_ptr<AdaptiveMeshType> mesh, Layer layer, Index domain_node, Index next) :
        _mesh(std::move(mesh)),
        _layer(layer),
        _domain_node(domain_node),
        _next(next)
      {
      }

      bool operator==(const AdaptiveNeighborsAdjactor& other) const
      {
        return _domain_node == other._domain_node && _next == other._next && _layer == other._layer &&
               _mesh == other._mesh;
      }

      bool operator!=(const AdaptiveNeighborsAdjactor& other) const
      {
        return !(*this == other);
      }

      AdaptiveNeighborsAdjactor& operator++()
      {
        _next++;
        return *this;
      }

      Index operator*() const
      {
        return _mesh->get_neighbor(_layer, _domain_node, _next);
      }
    };

    /// Adjactor type
    using ImageIterator = AdaptiveNeighborsAdjactor;

  private:
    // Mesh reference
    std::optional<Adjacency::Permutation> _permutation;
    std::optional<Adjacency::Permutation> _inverse_face_permutation;

  public:
    /**
     * \brief Constructor
     *
     * \param[in] mesh Underlying mesh
     * \param[in] Layer Layer this index set refers to
     */
    AdaptiveNeighbors(std::shared_ptr<AdaptiveMeshType> mesh, Layer layer) : _mesh(std::move(mesh)), _layer(layer)
    {
    }

    AdaptiveNeighbors(const AdaptiveNeighbors& other) = default;
    AdaptiveNeighbors(AdaptiveNeighbors&& other) = default;

    AdaptiveNeighbors& operator=(const AdaptiveNeighbors& other) = default;
    AdaptiveNeighbors& operator=(AdaptiveNeighbors&& other) = default;

    /**
     * \brief Returns a clone of this index set
     *
     * \warning A AdaptiveIndexSet is a view into a underlying AdaptiveMesh.
     * The returned clone will thus not be independent of this
     * AdaptiveIndexSet.
     */
    AdaptiveNeighbors clone() const
    {
      return AdaptiveNeighbors(_mesh, _layer);
    }

    /**
     * \brief Returns the size of dynamically allocated memory in bytes
     */
    std::size_t bytes() const
    {
      return 0;
    }

    /// access operator
    Index operator()(Index i, int j)
    {
      if(_permutation)
      {
        return _mesh->get_neighbor(
          _layer,
          _permutation.value().map(i),
          _inverse_face_permutation.value().map(j)
        );
      }
      return _mesh->get_neighbor(_layer, i, j);
    }

    /// access operator
    Index operator()(Index i, int j) const
    {
      if(_permutation)
      {
        return _mesh->get_neighbor(
          _layer,
          _permutation.value().map(i),
          _inverse_face_permutation.value().map(j)
        );
      }
      return _mesh->get_neighbor(_layer, i, j);
    }

    /**
     * \brief Returns the number of entities.
     */
    Index get_num_entities() const
    {
      return _mesh->get_num_entities(_layer, AdaptiveMeshType::ShapeType::dimension);
    }

    /**
     * \brief Returns the number of indices per entity
     */
    int get_num_indices() const
    {
      return num_indices;
    }

    /**
     * \brief Returns the maximum index of any face returned by this index set
     */
    Index get_index_bound() const
    {
      return _mesh->get_num_entities(_layer, AdaptiveMeshType::ShapeType::dimension);
    }

    // Unsupported.
    // IndexSet API depends on pointer arithmetic that we can not support.
    // Rewrite IndexSet API first, if this is required somewhere
    //AdaptiveIndexTuple* get_indices()
    //const AdaptiveIndexTuple* get_indices() const

    void permute(const Adjacency::Permutation& perm, const Adjacency::Permutation& inv_perm_face)
    {
      // Permute index tuples
      if(_permutation)
      {
        _permutation.value().concat(perm);
      }
      else
      {
        _permutation = perm.clone();
      }

      // Permute indices by the inverse face permutation
      if(_inverse_face_permutation)
      {
        _inverse_face_permutation.value().concat(inv_perm_face);
      }
      else
      {
        _inverse_face_permutation = inv_perm_face.clone();
      }
    }

    static String name()
    {
      return "AdaptiveNeighbors";
    }

    // Adjactor Interface

    Index get_num_nodes_domain() const
    {
      return get_num_entities();
    }

    Index get_num_nodes_image() const
    {
      return _mesh->get_num_entities(_layer, AdaptiveMeshType::ShapeType::dimension);
    }

    AdaptiveNeighborsAdjactor image_begin(Index domain_node) const
    {
      return AdaptiveIndexSetAdjactor(_mesh, _layer, domain_node, 0);
    }

    AdaptiveNeighborsAdjactor image_end(Index domain_node) const
    {
      return AdaptiveIndexSetAdjactor(_mesh, _layer, domain_node, num_indices);
    }
  };

  /**
   * \brief ConformalMesh interface for adaptive meshes
   *
   * This class template is a compatability layer between the AdaptiveMesh
   * class and the rest of FEAT. It presents a normal ConformalMesh interface
   * for a single layer of an AdaptiveMesh. This allows using the
   * AdaptiveMeshLayer like any other mesh class.
   *
   * \author Markus Muegge
   */
  template<typename AdaptiveMeshType_>
  class AdaptiveMeshLayer
  {
  public:
    /// Type of underlying mesh
    using AdaptiveMeshType = AdaptiveMeshType_;

    /// Type of neighbors wrapper
    using AdaptiveNeighborsType = AdaptiveNeighbors<AdaptiveMeshType>;

    // Conformal Mesh Interface

    /// Shape type
    using ShapeType = typename AdaptiveMeshType::ShapeType;

    /// Vertex set type
    using VertexSetType = AdaptiveVertexSet<AdaptiveMeshType>;

    /// Vertex type
    using VertexType = typename AdaptiveMeshType::VertexType;

    /// IndexSetHolder type
    using IndexSetHolderType = AdaptiveIndexSetHolder<AdaptiveMeshType, ShapeType>;

    /// Shape dimension
    static constexpr int shape_dim = AdaptiveMeshType::shape_dim;

    /// World dimension
    static constexpr int world_dim = AdaptiveMeshType::world_dim;

    /// This is an unstructured mesh
    static constexpr bool is_structured = false;

    /// Mesh permutation type
    using MeshPermutationType = MeshPermutation<ShapeType>;

    /// Used for return type of get_index_set
    template<int cell_dim_, int face_dim_>
    struct IndexSet
    {
      using Type = AdaptiveIndexSet<AdaptiveMeshType, cell_dim_, face_dim_>;
    };

  private:
    std::shared_ptr<AdaptiveMeshType> _mesh;
    Layer _layer;

    VertexSetType _vertex_set;
    IndexSetHolderType _index_set_holder;

    MeshPermutation<ShapeType> _permutation;

  public:
    /**
     * \brief Constructor
     *
     * \param[in] mesh Underlying adaptive mesh
     * \param[in] Layer Mesh layer
     */
    AdaptiveMeshLayer(std::shared_ptr<AdaptiveMeshType> mesh, Layer layer) :
      _mesh(mesh),
      _layer(layer),
      _vertex_set(VertexSetType(_mesh, layer)),
      _index_set_holder(_mesh, layer),
      _permutation()
    {
    }

    AdaptiveMeshLayer(const AdaptiveMeshLayer& other) = default;
    AdaptiveMeshLayer(AdaptiveMeshLayer&& other) = default;

    AdaptiveMeshLayer& operator=(const AdaptiveMeshLayer& other) = default;
    AdaptiveMeshLayer& operator=(AdaptiveMeshLayer&& other) = default;

    void clone(const AdaptiveMeshLayer& other)
    {
      _mesh = other._mesh;
      _layer = other._layer;
      _vertex_set = other._vertex_set;
      _index_set_holder = other._index_set_holder;
      _permutation = other._permutation.clone();
    }

    AdaptiveMeshLayer clone() const
    {
      AdaptiveMeshLayer result(_mesh, _layer);
      result._vertex_set = _vertex_set;
      result._index_set_holder = _index_set_holder;
      result._permutation = _permutation;
      return result;
    }

    std::size_t bytes() const
    {
      return _vertex_set.bytes() + _index_set_holder.bytes() + _permutation.bytes();
    }

    /**
     * \brief Returns the number of entities of dimension \c dim in the mesh
     *
     * \param[in] dim Dimension to query
     */
    Index get_num_entities(int dim) const
    {
      return _mesh->get_num_entities(_layer, dim);
    }

    /**
     * \brief Returns the number of vertices in the mesh
     */
    Index get_num_vertices() const
    {
      return get_num_entities(0);
    }

    /**
     * \brief Returns the number of elements (highest dimension entities) in the mesh
     */
    Index get_num_elements() const
    {
      return get_num_entities(ShapeType::dimension);
    }

    /**
     * \brief Checks whether the mesh is permuted.
     *
     * \returns \c true if the mesh is permuted, otherwise \c false.
     */
    bool is_permuted() const
    {
      return !_permutation.empty();
    }

    /**
     * \brief Returns the permutation applied to the mesh
     */
    const MeshPermutation<ShapeType>& get_mesh_permutation() const
    {
      return _permutation;
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

      // rebuild neighbor index set
      static_assert(false);
      // We need to calculate neighbors on the adaptive mesh
      // We also need to mock the neighbors, to support permutations
      //fill_neighbors();
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
      const std::array<Index, 4> mesh_size {
        get_num_entities(0),
        get_num_entities(1),
        get_num_entities(2),
        get_num_entities(3)
      };
      XASSERTM(mesh_perm.validate_sizes(mesh_size.data()) == 0, "mesh permutation has invalid size!");

      // save the permutation
      this->_permutation = std::forward<MeshPermutationType>(mesh_perm);

      // permute vertex set
      this->_vertex_set.permute(this->_permutation.get_perm(0));

      // permute index sets
      this->_index_set_holder.permute(this->_permutation.get_perms(), this->_permutation.get_inv_perms());

      // rebuild neighbor index set
      static_assert(false);
      // We need to calculate neighbors on the adaptive mesh
      // We also need to mock the neighbors, to support permutations
      //fill_neighbors();
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

    void fill_neighbours()
    {
      _mesh->fill_neighbors();
    }

    AdaptiveNeighborsType get_neighbors()
    {
      return AdaptiveNeighbors(_mesh, _layer);
    }

    AdaptiveNeighborsType get_neighbors() const
    {
      return AdaptiveNeighbors(_mesh, _layer);
    }

    /**
     * \brief Returns the vertex set of this mesh
     */
    VertexSetType& get_vertex_set()
    {
      return _vertex_set;
    }

    /**
     * \brief Returns the vertex set of this mesh
     */
    const VertexSetType& get_vertex_set() const
    {
      return _vertex_set;
    }

    /**
     * \brief Returns a reference to an index set of this mesh
     *
     * \tparam cell_dim_ Cell dimension of the index set
     * \tparam face_dim_ Face dimension of the index set
     */
    template<int cell_dim_, int face_dim_>
    auto& get_index_set()
    {
      return _index_set_holder.template get_index_set<cell_dim_, face_dim_>();
    }

    /**
     * \brief Returns a reference to an index set of this mesh
     *
     * \tparam cell_dim_ Cell dimension of the index set
     * \tparam face_dim_ Face dimension of the index set
     */
    template<int cell_dim_, int face_dim_>
    const auto& get_index_set() const
    {
      return _index_set_holder.template get_index_set<cell_dim_, face_dim_>();
    }

    /**
     * \brief Returns a reference to the index set holder of this mesh
     */
    IndexSetHolderType& get_index_set_holder()
    {
      return _index_set_holder;
    }

    /**
     * \brief Returns a reference to the index set holder of this mesh
     */
    const IndexSetHolderType& get_index_set_holder() const
    {
      return _index_set_holder;
    }

    IndexSetHolderType& get_topology()
    {
      return _index_set_holder;
    }

    const IndexSetHolderType& get_topology() const
    {
      return _index_set_holder;
    }

    // Wrapper specific code

    /**
     * \brief Retrieve index of a child element
     *
     * \param[in] element Index of parent
     * \param[in] child Child to retrieve index of. Should fulfill 0 <= child <= get_num_children(layer, parent_idx).
     *
     * \returns The index of the parent on the next layer, if the parent has no
     * children, the requested child's index, if the parent has at least \c
     * child children, std::nullopt else.
     */
    std::optional<Index> get_child(Index elem, Index child)
    {
      return _mesh->template get_child<ShapeType::dimension>(_layer, elem, child);
    }

    /**
     * \brief Returns the number of children of mesh a mesh element
     *
     * \param[in] elem Element to return number of children of
     */
    Index get_num_children(Index elem)
    {
      return _mesh->template get_num_children<ShapeType::dimension>(_layer, elem);
    }

    static String name()
    {
      return "AdaptiveMeshLayer<...>";
    }

    /**
     * \brief Indicates whether any mesh element adjacent to the given vertex has changed on the given layer
     *
     * This is intended to be used with the BPX solver to determine which vertices participate  on the given layer.
     */
    bool has_vertex_changed(Index vertex_idx) const
    {
      return _mesh->has_vertex_changed(_layer, vertex_idx);
    }

    /// accessor for foundation mesh
    typename AdaptiveMeshType::FoundationMeshType& foundation_mesh()
    {
      return _mesh->foundation_mesh();
    }

    /// accessor for adaptive mesh
    AdaptiveMeshType& adaptive_mesh()
    {
      return *_mesh;
    }
  };

  /**
   * \brief Adjactor for iterating over children of adaptive mesh elements
   *
   * \author Markus Muegge
   */
  template<typename AdaptiveMeshLayerType>
  class AdaptiveChildMapping
  {
  public:
    class ChildCellIterator
    {
    public:
      ChildCellIterator(AdaptiveMeshLayerType& coarse, Index cell, Index child) :
        _coarse(coarse),
        _cell(cell),
        _child(child)
      {
      }

      Index operator*() const
      {
        return _coarse.get_child(_cell, _child).value_or(0);
      }

      ChildCellIterator& operator++()
      {
        _child++;
        return *this;
      }

      bool operator!=(const ChildCellIterator& other)
      {
        return _cell != other._cell || _child != other._child;
      }

    private:
      static constexpr int dim = AdaptiveMeshLayerType::ShapeType::dimension;

      AdaptiveMeshLayerType& _coarse;
      Index _cell;
      Index _child;
    };

    using ImageIterator = ChildCellIterator;

  private:
    // TODO: Should these be std::shared_ptrs?
    AdaptiveMeshLayerType& _coarse;
    AdaptiveMeshLayerType& _fine;

  public:
    AdaptiveChildMapping(AdaptiveMeshLayerType& coarse, AdaptiveMeshLayerType& fine) : _coarse(coarse), _fine(fine)
    {
    }

    Index get_num_nodes_domain() const
    {
      return _coarse.get_num_elements();
    }

    Index get_num_nodes_image() const
    {
      return _fine.get_num_elements();
    }

    ChildCellIterator image_begin(Index domain_node) const
    {
      return ChildCellIterator(_coarse, domain_node, 0);
    }

    ChildCellIterator image_end(Index domain_node) const
    {
      Index num_children = _coarse.get_num_children(domain_node);
      return ChildCellIterator(_coarse, domain_node, num_children);
    }
  };
} // namespace FEAT::Geometry
