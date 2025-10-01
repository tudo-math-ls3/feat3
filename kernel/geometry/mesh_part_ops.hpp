// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

/**
 * @file mesh_part_ops.hpp
 * \brief Contains boolean operations on meshparts (intersection, union, difference, symmetric difference)
 *
 * The core idea of this implementation is a sort of join over mesh entities, inspired by database joins.
 * The join contains an entry for each entity in either of the input meshparts.
 * That entry contains the following information:
 * - is the entitiy contained in the lefthand-side meshpart and what is its index there
 * - is the entitiy contained in the righthand-side meshpart and what is its index there
 * - is the entitiy contained in the output meshpart and what is its index there
 * After this join has been calculated all other operations are rather simple.
 *
 * We can create the target sets with a simple predicate.
 *
 * Selecting all entities that are contained in both input meshpart creates an intersection
 * Selecting all entities that are contained in either input meshpart creates an union.
 * Selecting all entities that are contained in the lefthand-side meshpart but not the righthand-side creates a
 * difference.
 *
 * After selecting the elements we can create target sets by simply copying the mesh-indices from the correct input
 * target sets. During this step we also record the indicices of the selected entities in the result meshpart.
 *
 * We can copy topology information by looking at the topologies of the input meshparts (if they have one), searching of
 * any facets in the join, and then recording the indices of the entities in the result meshpart in the result meshparts
 * topology.
 *
 * We can copy attributes by copying from the input attributes and merging values whenever an entitiy occurs in both
 * input meshparts.
 */

// includes, feat
#include <kernel/base_header.hpp>
#include <kernel/geometry/index_set.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/target_set.hpp>
#include <kernel/shape.hpp>
#include <kernel/util/assertion.hpp>

// includes, system
#include <algorithm>
#include <functional>
#include <memory>
#include <optional>

namespace FEAT::Geometry
{
  namespace Inner
  {
    /**
     * \brief Utility class for storing one vector per dimension
     *
     * \tparam Shape_ Shape type, determines number of dimensions
     * \tparam T_ Type of values to store
     *
     * Implemented in the usual dimension-recursive way.
     *
     * \author Markus Muegge
     */
    template<typename Shape_, typename T_>
    class DimVector : public DimVector<typename Shape::FaceTraits<Shape_, Shape_::dimension - 1>::ShapeType, T_>
    {
    protected:
      /// Vector for current dimension
      std::vector<T_> vec;

    public:
      /**
       * \brief accessor
       *
       * \tparam int dim_ Which dimension's vector to retrieve
       */
      template<int dim_>
      std::vector<T_>& get()
      {
        static_assert(dim_ <= Shape_::dimension);

        using DimShape = typename Shape::FaceTraits<Shape_, dim_>::ShapeType;
        return DimVector<DimShape, T_>::vec;
      }

      /**
       * \brief const accessor
       *
       * \tparam int dim_ Which dimension's vector to retrieve
       */
      template<int dim_>
      const std::vector<T_>& get() const
      {
        static_assert(dim_ <= Shape_::dimension);

        using DimShape = typename Shape::FaceTraits<Shape_, dim_>::ShapeType;
        return DimVector<DimShape, T_>::vec;
      }
    };

    /**
     * \brief Utility class for storing one vector per dimension
     *
     * \tparam Shape_ Shape type, determines number of dimensions
     * \tparam T_ Type of values to store
     *
     * Implemented in the usual dimension-recursive way.
     * This is the base case of the recursion.
     *
     * \author Markus Muegge
     */
    template<typename T_>
    class DimVector<Shape::Vertex, T_>
    {
    protected:
      std::vector<T_> vec;

    public:
      template<int dim_>
      std::vector<T_>& get()
      {
        static_assert(dim_ == 0);

        return vec;
      }

      template<int dim_>
      const std::vector<T_>& get() const
      {
        static_assert(dim_ == 0);

        return vec;
      }
    };

    /**
     * \brief Helper-struct for storing relationships between meshpart entities
     *
     * The MeshPartOps class uses this to keep track of which indices refer to the
     * same mesh entity in the left or right input meshparts of a meshpart op and the
     * resulting output meshpart.
     *
     * These relationships are determined for all entities of either of the input meshparts.
     * Thus an entity might not have an index in the left or right meshparts.
     * An entity is also not necessarily part of the output meshpart, and thus
     * might not have an index in the output.
     * Hence all members are optional.
     *
     * The name JoinEntry is derived from join operations in databases.
     *
     * \author Markus Muegge
     */
    struct JoinEntry
    {
      /// Index of this entry in the left mesh part of the join, it it exists
      std::optional<Index> left = std::nullopt;

      /// Index of this entry in the right mesh part of the join, if it exists
      std::optional<Index> right = std::nullopt;

      /// Index of this entry in the result mesh part, after target sets have been built
      std::optional<Index> index = std::nullopt;
    };
  } // namespace Inner

  /**
   * \brief Boolean operations on mesh parts
   *
   * \tparam MeshType_ The type of the mesh underlying the meshparts, i.e.
   * this class template can operate on meshparts of type \c MeshPart<MeshType_>.
   *
   * Available operations are:
   * - intersections
   * - union
   * - difference
   * - symmetric difference
   *
   * This class template serves as a common namespace for boolean
   * operations on meshparts. The operations are essentially free functions,
   * this template is just for convenience to store common type aliases.
   *
   * \author Markus Muegge
   */
  template<typename MeshType_>
  class MeshPartOperations
  {
  private:
    /// Mesh type underlying all mesh parts
    using MeshType = MeshType_;

    /// Shape type
    using ShapeType = typename MeshType::ShapeType;

    /// Meshpart type
    using MeshPartType = MeshPart<MeshType>;

    /// Data types for meshpart attributes
    using AttributeDataType = typename MeshPartType::AttributeDataType;

    /// AttributeSet type
    using AttributeSetType = typename MeshPartType::AttributeSetType;

    /// Type for tracking relationships between meshpart entities
    using JoinEntry = Inner::JoinEntry;

    /// Collection of JoinEntries
    using JoinVec = std::vector<JoinEntry>;

    /// JoinVecs for all dimensions
    using Join = Inner::DimVector<ShapeType, JoinEntry>;

    /// TargetSetHolder type
    using TargetSetHolderType = TargetSetHolder<ShapeType>;

    /// IndexSetHolder type
    using IndexSetHolderType = IndexSetHolder<ShapeType>;

  public:
    /**
     * \brief Compute the intersection of two meshparts
     *
     * \tparam AttributeFn_       Type of attribute merging function.
     * Must be a callable that accepts a signature of <tt>DT(const DT&, const DT&)</tt>,
     * where <tt>DT = typename MeshPart<MeshType_>::AttributeDataType</tt>
     *
     * \param[in] left            Lefthand-side of intersection
     * \param[in] right           Righthand-side of intersection
     * \param[in] attribute_merge Attribute merging function
     *
     * Computes a new meshpart that contains all entities that occur in both \c left and \c right.
     *
     * Assume the following 2D mesh with quads A, B, C.
     *
     * \verbatim
     *  +--------+--------+--------+
     *  |        |        |        |
     *  |   A    |   B    |   C    |
     *  |        |        |        |
     *  +--------+--------+--------+
     * \endverbatim
     *
     * Further assume that \c left is a meshpart consisting of the quads A, B and all their vertices and edges
     * and that \c right is a meshpart consisting of quads B, C and all their vertices and edges.
     * Then the intersection of \c left and \c right consists of the quad B and all its vertices and edges.
     *
     * This function supports meshparts with topology. If \c left or \c right contain a topology,
     * then the applicable parts of the topologies will be copied to the result meshpart.
     * Precedence is given to \c right, if an entitiy appears in both meshparts.
     *
     * This function supports meshparts with attributes. If both \c left and \c right contain an attribute
     * with the same name and dimension, then an attribute with the same name and dimension will be created on
     * the result meshpart. For each value the attribute_merge function gets called to decide how to merge the existing
     * attribute values.
     */
    template<typename AttributeFn_>
    static MeshPartType
    meshpart_intersection(const MeshPartType& left, const MeshPartType& right, const AttributeFn_& attribute_merge)
    {
      const auto join_pred = [](const JoinEntry& e) { return e.left.has_value() && e.right.has_value(); };

      // The result meshpart has a topology as long as any of the inputs have a topology
      const bool create_topology = left.has_topology() || right.has_topology();
      return meshpart_op(left, right, attribute_merge, join_pred, create_topology);
    }

    /**
     * \brief Compute the union of two meshparts
     *
     * \tparam AttributeFn_       Type of attribute merging function.
     * Must be a callable that accepts a signature of <tt>DT(const DT&, const DT&)</tt>,
     * where <tt>DT = typename MeshPart<MeshType_>::AttributeDataType</tt>
     *
     * \param[in] left            Lefthand-side of union
     * \param[in] right           Righthand-side of union
     * \param[in] attribute_merge Attribute merging function
     *
     * Computes a new meshpart that contains all entities that occur in either \c left or \c right.
     *
     * Assume the following 2D mesh with quads A, B, C.
     *
     * \verbatim
     *  +--------+--------+--------+
     *  |        |        |        |
     *  |   A    |   B    |   C    |
     *  |        |        |        |
     *  +--------+--------+--------+
     * \endverbatim
     *
     * Further assume that \c left is a meshpart consisting of the quad A and all its vertices and edges
     * and that \c right is a meshpart consisting of quad C and all its vertices and edges.
     * Then the union of \c left and \c right consists of the quads A, B, and all their vertices and edges.
     *
     * This function supports meshparts with topology. If \c left or \c right contain a topology,
     * then the applicable parts of the topologies will be copied to the result meshpart.
     * Precedence is given to \c right, if an entitiy appears in both meshparts.
     *
     * This function supports meshparts with attributes. If both \c left and \c right contain an attribute
     * with the same name and dimension, then an attribute with the same name and dimension will be created on
     * the result meshpart. If an entity occurs in both meshparts, then the attribute_merge function gets called
     * to decide how to merge the existing attribute values. Otherwise the attribute values get copied as is.
     */
    template<typename AttributeFn_>
    static MeshPartType
    meshpart_union(const MeshPartType& left, const MeshPartType& right, const AttributeFn_& attribute_merge)
    {
      const auto join_pred = [](const JoinEntry& e) { return e.left.has_value() || e.right.has_value(); };
      // The result meshpart has a topology as long as both of the inputs have a topology
      const bool create_topology = left.has_topology() && right.has_topology();
      return meshpart_op(left, right, attribute_merge, join_pred, create_topology);
    }

    /**
     * \brief Compute the difference of two meshparts
     *
     * \param[in] left            Lefthand-side of difference
     * \param[in] right           Righthand-side of difference
     *
     * Computes a new meshpart that contains all entities that occur in \c left but not in \c right.
     *
     * Assume the following 2D mesh with quads A, B, C.
     *
     * \verbatim
     *  +--------+--------+--------+
     *  |        |        |        |
     *  |   A    |   B    |   C    |
     *  |        |        |        |
     *  +--------+--------+--------+
     * \endverbatim
     *
     * Further assume that \c left is a meshpart consisting of the quads A, B and all their vertices and edges
     * and that \c right is a meshpart consisting of quad B and all its vertices and edges.
     * Then the difference of \c left and \c right consists of the quad A, the leftmost two vertices of the mesh,
     * the leftmost edge of the mesh, and the top and bottom edges of A. Note that the right edge of A is not part of
     * the result meshpart, as it is an entity of \c right and has as such been removed.
     *
     * This function does not support mesparts with topology.
     *
     * This function supports meshparts with attributes. If both \c left and \c right contain an attribute
     * with the same name and dimension, then an attribute with the same name and dimension will be created on
     * the result meshpart. Values will be copied over from \c left.
     */
    static MeshPartType meshpart_difference(const MeshPartType& left, const MeshPartType& right)
    {
      // NOTE(mmuegge): Differences between meshparts with topology are technically supported.
      // They currently behave exactly as defined, which means they do not include any boundary entities
      // along the common border of the left and right mesh parts.
      // The implementation puts ~Index(0) values into the meshpart topology for these entities.
      // That is the proper operation, but not very useful in the context of meshparts
      // and quite surprising.
      // Hence the assertion below.
      // We could include these boundary entities in the difference.
      // To do so we would need to either determine how many additional entities
      // we need to add, so that we can size the meshpart topology appropriately,
      // or create a dynamic data structure for topologies.
      // Both of these solution are quite a lot of work, so we leave out differences between
      // meshparts with topology for now.
      XASSERTM(
        !left.has_topology() && !right.has_topology(),
        "Difference between meshparts with topology are not supported.");

      const auto merge = [](AttributeDataType& a, AttributeDataType& /*b*/) { return a; };
      const auto join_pred = [](const JoinEntry& e) { return e.left.has_value() && !e.right.has_value(); };
      return meshpart_op(left, right, merge, join_pred, false);
    }

    /**
     * \brief Compute the symmetric difference of two meshparts
     *
     * \param[in] left            Lefthand-side of symmetric difference
     * \param[in] right           Righthand-side of symmetric difference
     *
     * Computes a new meshpart that contains all entities that occur either in \c left but not \c right, or \c right but
     * not \c left.
     *
     * Assume the following 2D mesh with quads A, B, C.
     *
     * \verbatim
     *  +--------+--------+--------+
     *  |        |        |        |
     *  |   A    |   B    |   C    |
     *  |        |        |        |
     *  +--------+--------+--------+
     * \endverbatim
     *
     * Further assume that \c left is a meshpart consisting of the quads A, B and all their vertices and edges
     * and that \c right is a meshpart consisting of quads B, C and all their vertices and edges.
     * Then the symmetric difference of \c left and \c right consists of the quads A, C, the two leftmost vertices of
     * the mesh the two rightmost vertices of the mesh, and the boundary edges of A and C. Note that the right edge of A
     * and the left edge of C are not part of the result meshpart. They have been removed as part of quad B.
     *
     * This function does not support mesparts with topology.
     *
     * This function supports meshparts with attributes. If both \c left and \c right contain an attribute
     * with the same name and dimension, then an attribute with the same name and dimension will be created on
     * the result meshpart. Values will be copied over from \c left and \c right.
     */
    static MeshPartType meshpart_symmetric_difference(const MeshPartType& left, const MeshPartType& right)
    {
      // Identity merge, because no common entities end up in result
      const auto merge = [](AttributeDataType& a, AttributeDataType& /*b*/) { return a; };
      MeshPartType a = meshpart_difference(left, right, merge);
      MeshPartType b = meshpart_difference(right, left, merge); // NOLINT
      return meshpart_union(a, b, merge, false);
    }

  private:
    /**
     * \brief Common logic for all meshpart operations
     *
     * \tparam AttributeFn_ Type of attribute merging function. See actual docs of operations for details
     * \tparam PredFn_      Type of join predicate. Callable with signature <tt>bool(const JoinEntry&)</tt>
     *
     * \param[in] left  Lefthand-side of operation
     * \param[in] right Righthand-side of operation
     * \param[in] attribute_merge Attribute merging function
     * \param[in] join_pred Join predicate
     *
     * This function forms the common basis for all current meshpart operations.
     * The \c left, \c right, and \c attribute_merge parameters are the same as the operations parameters.
     * The \c join_pred determines the actual operation, by selecting which entities should be included in the final
     * meshpart. A \c join_pred that only selects entities that occur in both \c left and \c right implements an
     * intersection, a \c join_pred that selects all entities implements a union, etc.
     */
    template<typename AttributeFn_, typename PredFn_>
    static MeshPartType meshpart_op(
      const MeshPartType& left,
      const MeshPartType& right,
      const AttributeFn_& attribute_merge,
      const PredFn_& join_pred,
      const bool create_topology)
    {
      // Determine join of mesh parts
      // The join tracks which indices belong to the same entity on the left, right, and result meshparts.
      Join join;
      join_elements(join, left.get_target_set_holder(), right.get_target_set_holder());

      // The join predicate determines which elements get added to the target set holder
      TargetSetHolderType tsh;
      build_target_set_holder(tsh, join, join_pred, left.get_target_set_holder(), right.get_target_set_holder());

      // Determine result size from target set
      std::array<Index, ShapeType::dimension + 1> size{};
      for(int i(0); i < ShapeType::dimension + 1; i++)
      {
        size[i] = tsh.get_num_entities(i);
      }


      // Create result mesh part
      MeshPartType result(size.data(), create_topology);
      result.get_target_set_holder() = std::move(tsh);

      if(create_topology)
      {
        // Empty dummy IndexSetHolder. Allows us to handle all configurations of has/has-no topology the same
        IndexSetHolderType empty;
        const IndexSetHolderType* left_ish = left.has_topology() ? left.get_topology() : &empty;
        const IndexSetHolderType* right_ish = right.has_topology() ? right.get_topology() : &empty;

        build_topology(
          *result.get_topology(),
          join,
          *left_ish,
          left.get_target_set_holder(),
          *right_ish,
          right.get_target_set_holder());
      }

      // Handle attributes
      copy_attributes(attribute_merge, join, result, left, right);

      // Return intersection mesh part
      return std::move(result);
    }

    /**
     * \brief Find common attributes of two meshparts
     *
     * \param[in] left   Left meshpart
     * \param[in] right  Right meshpart
     *
     * \returns Set of all attributes that occur in both meshparts with the same name and dimension
     */
    static std::set<String> common_attributes(const MeshPartType& left, const MeshPartType& right)
    {
      std::set<String> attributes;

      for(const auto& pair : right.get_mesh_attributes())
      {
        const AttributeSetType* attribute = left.find_attribute(pair.first);
        if(attribute != nullptr && pair.second->get_dimension() == attribute->get_dimension())
        {
          attributes.insert(pair.first);
        }
      }

      return attributes;
    }

    /**
     * \brief Attribute copy logic
     *
     * \tparam AttributeFn_ Type of attribute merge function
     *
     * \param[in]    fn      Attribute merge function instance
     * \param[in]    join    Entity join. Used to determine where to pull attribute values from
     * \param[inout] result  The result meshpart. Attributes will be added to this meshpart.
     * \param[in]    left    Left meshpart of operation
     * \param[in]    right   Right meshpart of operation
     *
     * Copies, and merges (if required), attributes that exist with the same dimension on both
     * \c left and \c right.
     */
    template<typename AttributeFn_>
    static void copy_attributes(
      const AttributeFn_& fn,
      const Join& join,
      MeshPartType& result,
      const MeshPartType& left,
      const MeshPartType& right)
    {
      // Determine attributes that exist on both the left and right meshpart
      std::set<String> attributes = common_attributes(left, right);

      // Copy common attributes to result meshpart
      for(const String& name : attributes)
      {
        const AttributeSetType* left_attribute = left.find_attribute(name);
        const AttributeSetType* right_attribute = right.find_attribute(name);

        const int dim = left_attribute->get_dimension();

        auto new_attribute = std::make_unique<AttributeSetType>(result.get_num_entities(0), dim);

        const JoinVec& join_vec = join.template get<0>();

        for(const JoinEntry& t : join_vec)
        {
          if(!t.index.has_value())
          {
            // Entity is not in result meshpart. Skip it.
            continue;
          }

          if(t.left.has_value() && t.right.has_value())
          {
            // Entity exists on both meshpart. Merge the attribute values
            for(int i(0); i < dim; i++)
            {
              AttributeDataType l = left_attribute->operator()(t.left.value(), i);
              AttributeDataType r = right_attribute->operator()(t.right.value(), i);
              new_attribute->operator()(t.index.value(), i) = fn(l, r);
            }
          }
          else if(t.left.has_value())
          {
            // Entity exists only on left meshpart. Copy the attribute value
            for(int i(0); i < dim; i++)
            {
              AttributeDataType l = left_attribute->operator()(t.left.value(), i);
              new_attribute->operator()(t.index.value(), i) = l;
            }
          }
          else if(t.right.has_value())
          {
            // Entity exists only on right meshpart. Copy the attribute value
            for(int i(0); i < dim; i++)
            {
              AttributeDataType r = right_attribute->operator()(t.right.value(), i);
              new_attribute->operator()(t.index.value(), i) = r;
            }
          }
          else
          {
            XABORTM("Entity in meshpart operation result without source index");
          }
        }
        result.add_attribute(std::move(new_attribute), name);
      }
    }

    /**
     * \brief Outer wrapper for topology construction
     *
     * \tparam dim_ Current dimension
     *
     * \param[inout] output_ish  Output index set holder
     * \param[in]   join         Entity join. Used to determine where to pull topology from
     * \param[in]   left_ish     Left index set holder
     * \param[in]   left_tsh     Left target set holder
     * \param[in]   right_ish    Right index set holder
     * \param[in]   right_tsh    Right target set holder
     *
     * In order to build the topology of the result meshpart, each a@b pair of dimension needs to be visited.
     * This wrapper recursively runs through the b dimensions.
     */
    template<int dim_ = ShapeType::dimension>
    static void build_topology(
      IndexSetHolderType& output_ish,
      Join& join,
      const IndexSetHolderType& left_ish,
      const TargetSetHolderType& left_tsh,
      const IndexSetHolderType& right_ish,
      const TargetSetHolderType& right_tsh)
    {
      if constexpr(dim_ == 0)
      {
        // Nothing to do for vertices
        return;
      }
      else
      {
        build_topology<dim_ - 1>(output_ish, join, left_ish, left_tsh, right_ish, right_tsh);
        build_topology_inner<dim_, dim_ - 1>(output_ish, join, left_ish, left_tsh, right_ish, right_tsh);
      }
    }

    /**
     * \brief Inner logc for topology construction
     *
     * \tparam diim_  Current "cell" dimension
     * \tparam codim_ Current "face" dimension
     *
     * \param[inout] output_ish  Output index set holder
     * \param[in]   join         Entity join. Used to determine where to pull topology from
     * \param[in]   left_ish     Left index set holder
     * \param[in]   left_tsh     Left target set holder
     * \param[in]   right_ish    Right index set holder
     * \param[in]   right_tsh    Right target set holder
     *
     * In order to build the topology of the result meshpart, each a@b pair of dimension needs to be visited.
     * This wrapper recursively runs through the a dimensions and constructs the actual topologies.
     */
    template<int dim_, int codim_>
    static void build_topology_inner(
      IndexSetHolderType& output_ish,
      Join& join,
      const IndexSetHolderType& left_ish,
      const TargetSetHolderType& left_tsh,
      const IndexSetHolderType& right_ish,
      const TargetSetHolderType& right_tsh)
    {
      if constexpr(codim_ < 0)
      {
        // There are no facets with dimension lower than 0. No more work to be done.
        return;
      }
      else
      {
        // Recursive call to sub entities
        build_topology_inner<dim_, codim_ - 1>(output_ish, join, left_ish, left_tsh, right_ish, right_tsh);

        using DimShape = typename Shape::FaceTraits<ShapeType, dim_>::ShapeType;
        constexpr int num_indices = Shape::FaceTraits<DimShape, codim_>::count;
        using IndexSetType = IndexSet<num_indices>;

        JoinVec& dim_entities = join.template get<dim_>();
        JoinVec& codim_entities = join.template get<codim_>();

        const IndexSetType& right_is = right_ish.template get_index_set<dim_, codim_>();
        const IndexSetType& left_is = left_ish.template get_index_set<dim_, codim_>();

        IndexSetType& out = output_ish.template get_index_set<dim_, codim_>();

        // NOTE(mmuegge): It is tempting to turn the JoinVecs into hashmaps,
        // but we first search for right-indices and then left-indicies.
        // This means we need to duplicate information two maps with different keys.
        // So this might actually be more elegant.

        //////////////////
        // Right Search
        /////////////////

        // In the following, we need to search for entities in the codim_entities vector.
        // To do so efficiently we sort it beforehand, to allow for binary searches later
        std::sort(
          codim_entities.begin(),
          codim_entities.end(),
          [](const JoinEntry& a, const JoinEntry& b) { return a.right < b.right; });

        // We first search for topology in the entities of the right meshpart.
        // If, for any entitiy in the result meshpart, we do not find the topology there,
        // we need additionally search in the left meshpart.
        // Because we need to re-sort the codim vector, we do these steps one after the other
        // and use this bool to indicate whether we need to do the left-search at all.
        bool left_entities = false;

        // Determine topology for all entities that are part of the result
        // mesh part and stem from the right mesh part
        for(const JoinEntry& t : dim_entities)
        {
          if(!t.index.has_value())
          {
            // This entry is not part of the result mesh part. Skip it.
            continue;
          }

          if(!t.right.has_value())
          {
            // This entry does not stem from the right mesh part. Skip it.

            // This also means we need to handle entities stemming from the left
            // mesh part after this
            left_entities = true;
            continue;
          }

          for(int j(0); j < num_indices; j++)
          {
            // Search for join entry of j-th sub-entity.

            // NOTE(mmuegge): We are searching among the "std::optional<Index> right"
            // members of the join entries. We thus wrap our target value itself into
            // an optional before starting the search

            const std::optional<Index> target(right_is[t.right.value()][j]);
            std::optional<Index> codim_index = search(
              codim_entities.begin(),
              codim_entities.end(),
              target,
              [](const JoinEntry& a, const std::optional<Index>& i) { return a.right < i; });
            XASSERT(codim_index.has_value());

            const JoinEntry& codim_entry = codim_entities[codim_index.value()];

            // NOTE(mmuegge): codim_entry.index is only empty when trying to produce differences on meshparts.
            // In that case we could be searching for a topology element that has been removed as part of the
            // shared boundary between left and right.
            // That case is forbidden via an assertion above, but is safely handled here anyway.
            out(t.index.value(), j) = codim_entry.index.value_or(~Index(0));
          }
        }

        if(!left_entities)
        {
          return;
        }

        //////////////////
        // Left Search
        /////////////////

        // Prepare join for searching entities of right mesh part
        std::sort(
          codim_entities.begin(),
          codim_entities.end(),
          [](const JoinEntry& a, const JoinEntry& b) { return a.left < b.left; });

        for(const JoinEntry& t : dim_entities)
        {
          if(!t.index.has_value() || t.right.has_value() || !t.left.has_value())
          {
            // We only handle entries that are part of the result,
            // part of the left meshpart,
            // and have not been handled by the right meshpart case above
            continue;
          }

          for(int j(0); j < num_indices; j++)
          {
            const std::optional<Index> target(left_is[t.left.value()][j]);
            std::optional<Index> codim_index = search(
              codim_entities.begin(),
              codim_entities.end(),
              target,
              [](const JoinEntry& a, const std::optional<Index>& i) { return a.left < i; });
            XASSERT(codim_index.has_value());

            const JoinEntry& codim_entry = codim_entities[codim_index.value()];
            out(t.index.value(), j) = codim_entry.index.value_or(~Index(0));
          }
        }
      }
    }

    /**
     * \brief Construct target sets out of the join
     *
     * \tparam dim_ Current dimension
     *
     * \param[inout] output_tsh  Output target set holder. Target sets are written to this.
     * \param[inout] join        Entity join. For entities that are part of the result meshpart, target set indices are added
     * to the join.
     * \param[in] pred      Join predicate
     * \param[in] left_tsh  Left target set holder
     * \param[in] right_tsh Right target set holder
     *
     * This function walks over the join and uses the join predicate to determine what entities should be added to the
     * result meshpart. If it finds an entity that should be added, it pulls the correct mesh-index from the left or
     * right target set.
     */
    template<int dim_ = ShapeType::dimension>
    static void build_target_set_holder(
      TargetSetHolderType& output_tsh,
      Join& join,
      const std::function<bool(const JoinEntry&)>& pred,
      const TargetSetHolderType& left_tsh,
      const TargetSetHolderType& right_tsh)
    {
      if constexpr(dim_ < 0)
      {
        // No more work to be done
        return;
      }
      else
      {
        // Recurse to lower dimension
        build_target_set_holder<dim_ - 1>(output_tsh, join, pred, left_tsh, right_tsh);

        JoinVec& entities = join.template get<dim_>();
        const TargetSet& right_ts = right_tsh.template get_target_set<dim_>();
        const TargetSet& left_ts = left_tsh.template get_target_set<dim_>();

        std::vector<Index> target_set_vec;
        for(JoinEntry& t : entities)
        {
          if(pred(t))
          {
            // Entity should be added; assign it an index
            t.index = target_set_vec.size();
            if(t.left)
            {
              // Pull target set information from left target set
              target_set_vec.push_back(left_ts[t.left.value()]);
            }
            else if(t.right)
            {
              // Pull target set information from right target set
              target_set_vec.push_back(right_ts[t.right.value()]);
            }
            else
            {
              // NOTE(mmuegge): This seems at first glance like the starting place
              // for a complement operation. But I think we need to special case that,
              // because we would need to grab any entitity not in one of the meshparts
              // here and keep track of which entities we have already added.
              XABORTM("No source index for entity!");
            }
          }
        }

        // Copy to actual target set, now that we know the final size.
        TargetSet result(target_set_vec.size());
        for(Index i(0); i < target_set_vec.size(); i++)
        {
          result[i] = target_set_vec[i];
        }

        // Move into output target set holder
        output_tsh.template get_target_set<dim_>() = std::move(result);
      }
    }

    /**
     * \brief Construct a join of entities
     *
     * \tparam dim_ Current dimension
     *
     * \param[inout] join    Join to write to
     * \param[in] left_tsh   Left target set holder
     * \param[in] right_tsh  Right target set holder
     *
     * This function construct the entity join.
     * The join tracks which entities of the left and right meshparts correspond to the same
     * mesh entity, (and later also to which entity of the result meshpart).
     * This is calculated for all entities, which allows implementing the different boolean operations via
     * a simple predicate later.
     */
    template<int dim_ = ShapeType::dimension>
    static void join_elements(Join& join, const TargetSetHolderType& left_tsh, const TargetSetHolderType& right_tsh)
    {
      if constexpr(dim_ < 0)
      {
        return;
      }
      else
      {
        // Recurse to lower dimension
        join_elements<dim_ - 1>(join, left_tsh, right_tsh);

        JoinVec& join_vec = join.template get<dim_>();

        const TargetSet& left = left_tsh.template get_target_set<dim_>();
        const TargetSet& right = right_tsh.template get_target_set<dim_>();

        //////////////////////////////////////////////
        // Prepare lookup vector for left mesh part
        //////////////////////////////////////////////

        // Create a copy of the target set we can safely modify
        TargetSet left_lookup = left.clone();

        Index* lbegin = left_lookup.get_indices();
        Index* lend = left_lookup.get_indices() + left_lookup.get_num_entities();

        // The target set will be sorted to enable binary searches.
        // We need a way to retrieve the original indices later.
        // We thus create this vector of indices and sort it the same way
        // we sort the target set.
        std::vector<Index> left_indices(left_lookup.get_num_entities());
        for(Index i(0); i < left_lookup.get_num_entities(); i++)
        {
          left_indices[i] = i;
        }

        std::sort(
          left_indices.begin(),
          left_indices.end(),
          [&](Index a, Index b) { return left_lookup[a] < left_lookup[b]; });
        std::sort(lbegin, lend);

        /////////////////////////////////////////////////////////////////////////////
        // Find entities exclusive to right mesh part and entities in intersection
        /////////////////////////////////////////////////////////////////////////////

        for(Index i(0); i < right.get_num_entities(); i++)
        {
          JoinEntry entry;
          entry.right = i;

          const std::optional<Index> sorted_left = search(lbegin, lend, right[i]);
          if(sorted_left)
          {
            entry.left = left_indices[sorted_left.value()];
          }

          join_vec.push_back(entry);
        }

        // OPTIMIZATION(mmuegge): At this point we have sufficient information to
        // determine the intersection of two mesh parts. We calculate the full join
        // only for other mesh part operations. If mesh part intersection ever
        // becomes a performance bottleneck, we can calculate these joins more fine grained.

        //////////////////////////////////////////////
        // Prepare lookup vector for left mesh part
        //////////////////////////////////////////////

        // Create a copy of the target set we can safely modify
        TargetSet right_lookup = right.clone();

        Index* rbegin = right_lookup.get_indices();
        Index* rend = right_lookup.get_indices() + right_lookup.get_num_entities();

        // The target set will be sorted to enable binary searches.
        // We need a way to retrieve the original indices later.
        // We thus create this vector of indices and sort it the same way
        // we sort the target set.
        std::vector<Index> right_indices(right_lookup.get_num_entities());
        for(Index i(0); i < right_lookup.get_num_entities(); i++)
        {
          right_indices[i] = i;
        }

        std::sort(
          right_indices.begin(),
          right_indices.end(),
          [&](Index a, Index b) { return right_lookup[a] < right_lookup[b]; });
        std::sort(rbegin, rend);

        //////////////////////////////////////////////
        // Find entities exclusive to left mesh part
        //////////////////////////////////////////////

        for(Index i(0); i < left.get_num_entities(); i++)
        {
          const std::optional<Index> sorted_right = search(rbegin, rend, left[i]);
          if(!sorted_right.has_value())
          {
            JoinEntry entry;
            entry.left = i;
            join_vec.push_back(entry);
          }
        }
      }
    }

    /**
     * \brief Binary search
     *
     * \tparam Iter_ Iterator type
     * \tparam T_    Target value type
     *
     * \param[in] begin  Start iterator
     * \param[in] end    End iterator
     * \param[in] value  Target value
     *
     * \returns Index of the target value in the range [begin, end), if the range contains the value. nullopt otherwise.
     */
    template<typename Iter_, typename T_>
    static std::optional<Index> search(Iter_ begin, Iter_ end, const T_& value)
    {
      Iter_ it = std::lower_bound(begin, end, value);
      if(it == end || *it != value)
      {
        return std::nullopt;
      }
      else
      {
        return std::distance(begin, it);
      }
    }

    /**
     * \brief Binary search
     *
     * \tparam Iter_   Iterator type
     * \tparam T_      Target value type
     * \tparam CompFn_ Custom comparison function
     *
     * \param[in] begin  Start iterator
     * \param[in] end    End iterator
     * \param[in] value  Target value
     * \param[in] fn     Comparison function
     *
     * \returns Index of the target value in the range [begin, end), if the range contains the value. nullopt otherwise.
     */
    template<typename Iter_, typename T_, typename CompFn_>
    static std::optional<Index> search(Iter_ begin, Iter_ end, const T_& value, const CompFn_& fn)
    {
      Iter_ it = std::lower_bound(begin, end, value, fn);
      if(it == end || fn(*it, value))
      {
        return std::nullopt;
      }
      else
      {
        return std::distance(begin, it);
      }
    }
  };
} // namespace FEAT::Geometry
