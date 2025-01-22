// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/geometry/raw_refinement_templates.hpp>
#include <kernel/geometry/refinement_types.hpp>
#include <kernel/shape.hpp>

#include <unordered_map>

namespace FEAT::Geometry
{
  class SunZhaoMaExpansionData
  {
  public:
    /// Highest-dimensional shape enabled by this data
    using MaxShape = Shape::Hypercube<3>;

    /// Raw edge template type
    using RawEdgeTemplate = RawTemplate<Shape::Hypercube<1>>;
    /// Raw face template type
    using RawFaceTemplate = RawTemplate<Shape::Hypercube<2>>;
    /// Raw cell template type
    using RawCellTemplate = RawTemplate<Shape::Hypercube<3>>;

    /// Accessor for refinement types used by these templates
    template<int dim_>
    using RefinementTypeByDim = IsolatedPointRefinementType<Shape::Hypercube<dim_>>;

    /// Accessor for refinement types used by these templates
    template<typename Shape__>
    using RefinementTypeByShape = IsolatedPointRefinementType<Shape__>;

    /// Type of map from refinement type to raw templates
    template<int dim_>
    using TemplateMapByDim = std::unordered_map<RefinementTypeByDim<dim_>, RawTemplate<Shape::Hypercube<dim_>>>;

    /// Refinement type for edges
    using EdgeRefinementType = RefinementTypeByDim<1>;
    /// Refinement type for faces
    using FaceRefinementType = RefinementTypeByDim<2>;
    /// Refinement type for cells
    using CellRefinementType = RefinementTypeByDim<3>;

    /// Type of map from edge refinement types to raw templates
    using EdgeMap = TemplateMapByDim<1>;
    /// Type of map from face refinement types to raw templates
    using FaceMap = TemplateMapByDim<2>;
    /// Type of map from cell refinement types to raw templates
    using CellMap = TemplateMapByDim<3>;

    /**
     * \brief Shape compatability test
     *
     * This template set is compatible with quadrilateral and hexahedral meshes.
     *
     * \tparam Shape_ Mesh shape this template set is to be applied to
     *
     * \returns True, if the template set is compatible with the mesh shape, false otherwise.
     */
    template<typename Shape_>
    static constexpr bool is_shape_compatible()
    {
      return std::is_same_v<Shape_, Shape::Quadrilateral> || std::is_same_v<Shape_, Shape::Hexahedron>;
    }

    /**
     * \brief Constexpr function for retrieving maximum number of children for any template of this template set
     *
     * \tparam template_dim_ Dimension of the template
     * \tparam child_dim_ Dimension of the children for which the maximum number should be returned
     *
     * \note This is used by the AdaptiveMesh to determine the required storage
     * for child entities. The results of this function must be at least the
     * maximum. Larger values are allowed, but inefficient.
     */
    template<int template_dim_, int child_dim_>
    static constexpr int max_children()
    {
      static_assert(template_dim_ <= 3);
      static_assert(child_dim_ <= 3);
      const constexpr std::array<std::array<int, 4>, 4> children = {{
        {0, 0, 0, 0},
        {2, 3, 0, 0},
        {8, 18, 11, 0},
        {62, 197, 216, 82},
      }};
      return children[template_dim_][child_dim_];
    }

    /**
     * \brief Accessor for raw template maps
     *
     * \tparam dim_ Dimension of template to access
     */
    template<int dim_>
    static TemplateMapByDim<dim_>& raw_templates()
    {
      if constexpr(dim_ == 1)
      {
        return raw_edges();
      }
      if constexpr(dim_ == 2)
      {
        return raw_faces();
      }
      if constexpr(dim_ == 3)
      {
        return raw_cells();
      }
      XABORTM("SchneidersData supplied no templates of dimension " + stringify(dim_));
    }

    /// Returns the raw templates for edges of this template set
    static EdgeMap& raw_edges()
    {
      static EdgeMap result = {};

      using V = typename RawEdgeTemplate::VertexType;

      if(result.empty())
      {
        // No refinement unless either both vertices are marked or one of them is isolated
        result.insert({RefinementTypeByDim<1>(0b00, false), RawEdgeTemplate()});
        result.insert({RefinementTypeByDim<1>(0b10, false), RawEdgeTemplate().add_entity(V{0.0}, V{1.0})});
        result.insert({RefinementTypeByDim<1>(0b01, false), RawEdgeTemplate().add_entity(V{0.0}, V{1.0})});

        // Split edge in half, if adjacent to isolated point
        result.insert(
          {RefinementTypeByDim<1>(0b01, true),
           RawEdgeTemplate().add_entity(V{0.0}, V{1.0 / 2.0}).add_entity(V{1.0 / 2.0}, V{1.0})});

        result.insert(
          {RefinementTypeByDim<1>(0b10, true),
           RawEdgeTemplate().add_entity(V{0.0}, V{1.0 / 2.0}).add_entity(V{1.0 / 2.0}, V{1.0})});

        // Split in thirds if adjacent to two marked vertices
        result.insert(
          {RefinementTypeByDim<1>(0b11, false),
           RawEdgeTemplate()
             .add_entity(V{0.0}, V{1.0 / 3.0})
             .add_entity(V{1.0 / 3.0}, V{2.0 / 3.0})
             .add_entity(V{2.0 / 3.0}, V{1.0})});
      }

      return result;
    }

    /// Returns the raw templates for faces of this template set
    static FaceMap& raw_faces()
    {
      static FaceMap result;

      using V = typename RawFaceTemplate::VertexType;

      if(result.empty())
      {
        // No refinement if no vertices are marked
        result.insert({RefinementTypeByDim<2>(0b0000, false), RawFaceTemplate()});

        // Templates for isolated point refinement, with 4 different orientations
        result.insert(
          {RefinementTypeByDim<2>(0b0001, true),
           RawFaceTemplate()
             .add_entity(V{0.0, 0.0}, V{1.0 / 2.0, 0.0}, V{0.0, 1.0 / 2.0}, V{1.0 / 2.0, 1.0 / 2.0})
             .add_entity(V{1.0 / 2.0, 0.0}, V{1.0, 0.0}, V{1.0 / 2.0, 1.0 / 2.0}, V{1.0, 1.0})
             .add_entity(V{0.0, 1.0 / 2.0}, V{1.0 / 2.0, 1.0 / 2.0}, V{0.0, 1.0}, V{1.0, 1.0})});

        XASSERT(create_2dtemplate_rotations(result, RefinementTypeByDim<2>(0b0001, true)) == 4);

        // No refinement if corner vertex is not isolated
        result.insert(
          {RefinementTypeByDim<2>(0b0001, false), RawFaceTemplate().axis_aligned(V{0.0, 0.0}, V{1.0, 1.0})});
        result.insert(
          {RefinementTypeByDim<2>(0b0010, false), RawFaceTemplate().axis_aligned(V{0.0, 0.0}, V{1.0, 1.0})});
        result.insert(
          {RefinementTypeByDim<2>(0b0100, false), RawFaceTemplate().axis_aligned(V{0.0, 0.0}, V{1.0, 1.0})});
        result.insert(
          {RefinementTypeByDim<2>(0b1000, false), RawFaceTemplate().axis_aligned(V{0.0, 0.0}, V{1.0, 1.0})});

        // Templates for edge refinement, with 4 different orientations
        result.insert(
          {RefinementTypeByDim<2>(0b0011, false),
           RawFaceTemplate()
             .add_entity(V{0.0, 0.0}, V{1.0 / 3.0, 0.0}, V{0.0, 1.0}, V{1.0 / 3.0, 1.0 / 2.0})
             .add_entity(V{1.0 / 3.0, 0.0}, V{2.0 / 3.0, 0.0}, V{1.0 / 3.0, 1.0 / 2.0}, V{2.0 / 3.0, 1.0 / 2.0})
             .add_entity(V{2.0 / 3.0, 0.0}, V{1.0, 0.0}, V{2.0 / 3.0, 1.0 / 2.0}, V{1.0, 1.0})
             .add_entity(V{1.0 / 3.0, 1.0 / 2.0}, V{2.0 / 3.0, 1.0 / 2.0}, V{0.0, 1.0}, V{1.0, 1.0})});

        XASSERT(create_2dtemplate_rotations(result, RefinementTypeByDim<2>(0b0011, false)) == 4);

        // Template for two-edge refinement, with 4 different orientations
        result.insert(
          {RefinementTypeByDim<2>(0b0111, false),
           RawFaceTemplate()
             .add_entity(V{0.4, 0.4}, V{0.6, 0.4}, V{0.4, 0.6}, V{0.6, 0.6})
             .add_entity(V{0.4, 0.6}, V{0.6, 0.6}, V{0.0, 1.0}, V{1.0, 1.0})
             .add_entity(V{0.6, 0.6}, V{0.6, 0.4}, V{1.0, 1.0}, V{1.0, 0.0})
             .add_entity(V{0.6, 0.4}, V{0.4, 0.4}, V{1.0, 0.0}, V{0.0, 0.0})
             .add_entity(V{0.4, 0.4}, V{0.4, 0.6}, V{0.0, 0.0}, V{0.0, 1.0})
             .recurse(4, result[RefinementTypeByDim<2>(0b1100, 0b0000)])
             .recurse(3, result[RefinementTypeByDim<2>(0b1100, 0b0000)])});

        XASSERT(create_2dtemplate_rotations(result, RefinementTypeByDim<2>(0b0111, false)) == 4);

        // Opposite Corner refinements
        // No refinement, but valid for refinement fields
        result.insert(
          {RefinementTypeByDim<2>(0b0110, false), RawFaceTemplate().axis_aligned(V{0.0, 0.0}, V{1.0, 1.0})});
        result.insert(
          {RefinementTypeByDim<2>(0b1001, false), RawFaceTemplate().axis_aligned(V{0.0, 0.0}, V{1.0, 1.0})});

        // Template for full refinement
        result.insert({RefinementTypeByDim<2>(0b1111, false), RawFaceTemplate().grid({3, 3}, V{1.0 / 3.0, 1.0 / 3.0})});
      }

      return result;
    }

    /// Returns the raw templates for cells of this template set
    static CellMap& raw_cells()
    {
      static CellMap result;

      using V = typename RawCellTemplate::VertexType;

      if(result.empty())
      {
        // No refinement
        result.insert({RefinementTypeByDim<3>(0b00000000, false), RawCellTemplate()});

        // Isolated point refinement, in eight orientations
        RefinementTypeByDim<3> isolated_base_type(0b00000100, true);
        result.insert(
          {isolated_base_type,
           RawCellTemplate()
             .add_entity(
               V{0.0, 0.0, 0.0},
               V{1.0, 0.0, 0.0},
               V{0.0, 1.0 / 2.0, 0.0},
               V{1.0 / 2.0, 1.0 / 2.0, 0.0},
               V{0.0, 0.0, 1.0},
               V{1.0, 0.0, 1.0},
               V{0.0, 1.0 / 2.0, 1.0 / 2.0},
               V{1.0 / 2.0, 1.0 / 2.0, 1.0 / 2.0})
             .add_entity(
               V{1.0 / 2.0, 1.0 / 2.0, 0.0},
               V{1.0, 0.0, 0.0},
               V{1.0 / 2.0, 1.0, 0.0},
               V{1.0, 1.0, 0.0},
               V{1.0 / 2.0, 1.0 / 2.0, 1.0 / 2.0},
               V{1.0, 0.0, 1.0},
               V{1.0 / 2.0, 1.0, 1.0 / 2.0},
               V{1.0, 1.0, 1.0})
             .add_entity(
               V{0.0, 1.0 / 2.0, 1.0 / 2.0},
               V{1.0 / 2.0, 1.0 / 2.0, 1.0 / 2.0},
               V{0.0, 1.0, 1.0 / 2.0},
               V{1.0 / 2.0, 1.0, 1.0 / 2.0},
               V{0.0, 0.0, 1.0},
               V{1.0, 0.0, 1.0},
               V{0.0, 1.0, 1.0},
               V{1.0, 1.0, 1.0})
             .axis_aligned(V{0.0, 1.0 / 2.0, 0.0}, V{1.0 / 2.0, 1.0, 1.0 / 2.0})});

        XASSERT(create_3dtemplate_rotations(result, isolated_base_type) == 8);

        // Non-Isolated point refinement leads to no refinement
        RefinementTypeByDim<3> non_isolated_point_base_type(0b00000100, false);
        result.insert(
          {non_isolated_point_base_type, RawCellTemplate().axis_aligned(V{0.0, 0.0, 0.0}, V{1.0, 1.0, 1.0})});

        XASSERT(create_3dtemplate_rotations(result, non_isolated_point_base_type) == 8);

        // Edge refinement, in twelve orientations
        RefinementTypeByDim<3> edge_refinement_base_type(0b00000101, false);
        result.insert(
          {edge_refinement_base_type,
           RawCellTemplate()
             .add_entity(
               V{0.0, 0.0, 0.0},
               V{1.0, 0.0, 0.0},
               V{0.0, 1.0 / 3.0, 0.0},
               V{1.0 / 2.0, 1.0 / 3.0, 0.0},
               V{0.0, 0.0, 1.0},
               V{1.0, 0.0, 1.0},
               V{0.0, 1.0 / 3.0, 1.0 / 2.0},
               V{1.0 / 2.0, 1.0 / 3.0, 1.0 / 2.0})
             .add_entity(
               V{1.0 / 2.0, 1.0 / 3.0, 0.0},
               V{1.0, 0.0, 0.0},
               V{1.0 / 2.0, 2.0 / 3.0, 0.0},
               V{1.0, 1.0, 0.0},
               V{1.0 / 2.0, 1.0 / 3.0, 1.0 / 2.0},
               V{1.0, 0.0, 1.0},
               V{1.0 / 2.0, 2.0 / 3.0, 1.0 / 2.0},
               V{1.0, 1.0, 1.0})
             .add_entity(
               V{0.0, 2.0 / 3.0, 0.0},
               V{1.0 / 2.0, 2.0 / 3.0, 0.0},
               V{0.0, 1.0, 0.0},
               V{1.0, 1.0, 0.0},
               V{0.0, 2.0 / 3.0, 1.0 / 2.0},
               V{1.0 / 2.0, 2.0 / 3.0, 1.0 / 2.0},
               V{0.0, 1.0, 1.0},
               V{1.0, 1.0, 1.0})
             .add_entity(
               V{0.0, 1.0 / 3.0, 1.0 / 2.0},
               V{1.0 / 2.0, 1.0 / 3.0, 1.0 / 2.0},
               V{0.0, 2.0 / 3.0, 1.0 / 2.0},
               V{1.0 / 2.0, 2.0 / 3.0, 1.0 / 2.0},
               V{0.0, 0.0, 1.0},
               V{1.0, 0.0, 1.0},
               V{0.0, 1.0, 1.0},
               V{1.0, 1.0, 1.0})
             .add_entity(
               V{0.0, 1.0 / 3.0, 0.0},
               V{1.0 / 2.0, 1.0 / 3.0, 0.0},
               V{0.0, 2.0 / 3.0, 0.0},
               V{1.0 / 2.0, 2.0 / 3.0, 0.0},
               V{0.0, 1.0 / 3.0, 1.0 / 2.0},
               V{1.0 / 2.0, 1.0 / 3.0, 1.0 / 2.0},
               V{0.0, 2.0 / 3.0, 1.0 / 2.0},
               V{1.0 / 2.0, 2.0 / 3.0, 1.0 / 2.0})});

        XASSERT(create_3dtemplate_rotations(result, edge_refinement_base_type) == 12);

        // Face refinement, in six orientations
        RefinementTypeByDim<3> face_base_refinement(0b00001111, false);
        result.insert(
          {face_base_refinement,
           RawCellTemplate()
             .add_entity(
               V{0.0, 0.0, 0.0},
               V{1.0 / 3.0, 0.0, 0.0},
               V{0.0, 1.0 / 3.0, 0.0},
               V{1.0 / 3.0, 1.0 / 3.0, 0.0},
               V{0.0, 0.0, 1.0},
               V{1.0 / 3.0, 0.0, 1.0 / 2.0},
               V{0.0, 1.0 / 3.0, 1.0 / 2.0},
               V{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0})
             .add_entity(
               V{1.0 / 3.0, 0.0, 0.0},
               V{2.0 / 3.0, 0.0, 0.0},
               V{1.0 / 3.0, 1.0 / 3.0, 0.0},
               V{2.0 / 3.0, 1.0 / 3.0, 0.0},
               V{1.0 / 3.0, 0.0, 1.0 / 2.0},
               V{2.0 / 3.0, 0.0, 1.0 / 2.0},
               V{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
               V{2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0})
             .add_entity(
               V{2.0 / 3.0, 0.0, 0.0},
               V{1.0, 0.0, 0.0},
               V{2.0 / 3.0, 1.0 / 3.0, 0.0},
               V{1.0, 1.0 / 3.0, 0.0},
               V{2.0 / 3.0, 0.0, 1.0 / 2.0},
               V{1.0, 0.0, 1.0},
               V{2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
               V{1.0, 1.0 / 3.0, 1.0 / 2.0})
             .add_entity(
               V{0.0, 1.0 / 3.0, 0.0},
               V{1.0 / 3.0, 1.0 / 3.0, 0.0},
               V{0.0, 2.0 / 3.0, 0.0},
               V{1.0 / 3.0, 2.0 / 3.0, 0.0},
               V{0.0, 1.0 / 3.0, 1.0 / 2.0},
               V{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
               V{0.0, 2.0 / 3.0, 1.0 / 2.0},
               V{1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0})
             .add_entity(
               V{1.0 / 3.0, 1.0 / 3.0, 0.0},
               V{2.0 / 3.0, 1.0 / 3.0, 0.0},
               V{1.0 / 3.0, 2.0 / 3.0, 0.0},
               V{2.0 / 3.0, 2.0 / 3.0, 0.0},
               V{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
               V{2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
               V{1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0},
               V{2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0})
             .add_entity(
               V{2.0 / 3.0, 1.0 / 3.0, 0.0},
               V{1.0, 1.0 / 3.0, 0.0},
               V{2.0 / 3.0, 2.0 / 3.0, 0.0},
               V{1.0, 2.0 / 3.0, 0.0},
               V{2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
               V{1.0, 1.0 / 3.0, 1.0 / 2.0},
               V{2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0},
               V{1.0, 2.0 / 3.0, 1.0 / 2.0})
             .add_entity(
               V{0.0, 2.0 / 3.0, 0.0},
               V{1.0 / 3.0, 2.0 / 3.0, 0.0},
               V{0.0, 1.0, 0.0},
               V{1.0 / 3.0, 1.0, 0.0},
               V{0.0, 2.0 / 3.0, 1.0 / 2.0},
               V{1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0},
               V{0.0, 1.0, 1.0},
               V{1.0 / 3.0, 1.0, 1.0 / 2.0})
             .add_entity(
               V{1.0 / 3.0, 2.0 / 3.0, 0.0},
               V{2.0 / 3.0, 2.0 / 3.0, 0.0},
               V{1.0 / 3.0, 1.0, 0.0},
               V{2.0 / 3.0, 1.0, 0.0},
               V{1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0},
               V{2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0},
               V{1.0 / 3.0, 1.0, 1.0 / 2.0},
               V{2.0 / 3.0, 1.0, 1.0 / 2.0})
             .add_entity(
               V{2.0 / 3.0, 2.0 / 3.0, 0.0},
               V{1.0, 2.0 / 3.0, 0.0},
               V{2.0 / 3.0, 1.0, 0.0},
               V{1.0, 1.0, 0.0},
               V{2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0},
               V{1.0, 2.0 / 3.0, 1.0 / 2.0},
               V{2.0 / 3.0, 1.0, 1.0 / 2.0},
               V{1.0, 1.0, 1.0})
             .add_entity(
               V{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
               V{2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
               V{1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0},
               V{2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0},
               V{1.0 / 3.0, 0.0, 1.0 / 2.0},
               V{2.0 / 3.0, 0.0, 1.0 / 2.0},
               V{1.0 / 3.0, 1.0, 1.0 / 2.0},
               V{2.0 / 3.0, 1.0, 1.0 / 2.0})
             .add_entity(
               V{0.0, 1.0 / 3.0, 1.0 / 2.0},
               V{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
               V{0.0, 2.0 / 3.0, 1.0 / 2.0},
               V{1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0},
               V{0.0, 0.0, 1.0},
               V{1.0 / 3.0, 0.0, 1.0 / 2.0},
               V{0.0, 1.0, 1.0},
               V{1.0 / 3.0, 1.0, 1.0 / 2.0})
             .add_entity(
               V{2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0},
               V{1.0, 1.0 / 3.0, 1.0 / 2.0},
               V{2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0},
               V{1.0, 2.0 / 3.0, 1.0 / 2.0},
               V{2.0 / 3.0, 0.0, 1.0 / 2.0},
               V{1.0, 0.0, 1.0},
               V{2.0 / 3.0, 1.0, 1.0 / 2.0},
               V{1.0, 1.0, 1.0})
             .add_entity(
               V{1.0 / 3.0, 0.0, 1.0 / 2.0},
               V{2.0 / 3.0, 0.0, 1.0 / 2.0},
               V{1.0 / 3.0, 1.0, 1.0 / 2.0},
               V{2.0 / 3.0, 1.0, 1.0 / 2.0},
               V{0.0, 0.0, 1.0},
               V{1.0, 0.0, 1.0},
               V{0.0, 1.0, 1.0},
               V{1.0, 1.0, 1.0})});

        XASSERT(create_3dtemplate_rotations(result, face_base_refinement) == 6);

        // Two Edge Refinement, with 24? orientations
        RefinementTypeByDim<3> two_edge_base_refinement(0b00001101, false);
        result.insert(
          {two_edge_base_refinement,
           RawCellTemplate()
             .add_entity(
               V{0.0, 0.0, 0.0},
               V{2.0 / 5.0, 2.0 / 5.0, 0.0},
               V{0.0, 1.0, 0.0},
               V{2.0 / 5.0, 3.0 / 5.0, 0.0},
               V{0.0, 0.0, 1.0},
               V{2.0 / 5.0, 2.0 / 5.0, 2.0 / 3.0},
               V{0.0, 1.0, 1.0},
               V{2.0 / 5.0, 3.0 / 5.0, 2.0 / 3.0})
             .add_entity(
               V{2.0 / 5.0, 3.0 / 5.0, 0.0},
               V{3.0 / 5.0, 3.0 / 5.0, 0.0},
               V{0.0, 1.0, 0.0},
               V{1.0, 1.0, 0.0},
               V{2.0 / 5.0, 3.0 / 5.0, 2.0 / 3.0},
               V{3.0 / 5.0, 3.0 / 5.0, 2.0 / 3.0},
               V{0.0, 1.0, 1.0},
               V{1.0, 1.0, 1.0})
             .add_entity(
               V{3.0 / 5.0, 2.0 / 5.0, 0.0},
               V{1.0, 0.0, 0.0},
               V{3.0 / 5.0, 3.0 / 5.0, 0.0},
               V{1.0, 1.0, 0.0},
               V{3.0 / 5.0, 2.0 / 5.0, 2.0 / 3.0},
               V{1.0, 0.0, 1.0},
               V{3.0 / 5.0, 3.0 / 5.0, 2.0 / 3.0},
               V{1.0, 1.0, 1.0})
             .add_entity(
               V{0.0, 0.0, 0.0},
               V{1.0, 0.0, 0.0},
               V{2.0 / 5.0, 2.0 / 5.0, 0.0},
               V{3.0 / 5.0, 2.0 / 5.0, 0.0},
               V{0.0, 0.0, 1.0},
               V{1.0, 0.0, 1.0},
               V{2.0 / 5.0, 2.0 / 5.0, 2.0 / 3.0},
               V{3.0 / 5.0, 2.0 / 5.0, 2.0 / 3.0})
             .add_entity(
               V{2.0 / 5.0, 2.0 / 5.0, 2.0 / 3.0},
               V{3.0 / 5.0, 2.0 / 5.0, 2.0 / 3.0},
               V{2.0 / 5.0, 3.0 / 5.0, 2.0 / 3.0},
               V{3.0 / 5.0, 3.0 / 5.0, 2.0 / 3.0},
               V{0.0, 0.0, 1.0},
               V{1.0, 0.0, 1.0},
               V{0.0, 1.0, 1.0},
               V{1.0, 1.0, 1.0})
             .add_entity(
               V{2.0 / 5.0, 2.0 / 5.0, 0.0},
               V{3.0 / 5.0, 2.0 / 5.0, 0.0},
               V{2.0 / 5.0, 3.0 / 5.0, 0.0},
               V{3.0 / 5.0, 3.0 / 5.0, 0.0},
               V{2.0 / 5.0, 2.0 / 5.0, 2.0 / 3.0},
               V{3.0 / 5.0, 2.0 / 5.0, 2.0 / 3.0},
               V{2.0 / 5.0, 3.0 / 5.0, 2.0 / 3.0},
               V{3.0 / 5.0, 3.0 / 5.0, 2.0 / 3.0})
             .recurse(0, result[RefinementTypeByDim<3>(0b00000101, false)])
             // Note: Face 0 was replaced in the previous recurse call
             // We recurse on face 0 again, but mean the face 1 of before the first recursion
             .recurse(0, result[RefinementTypeByDim<3>(0b00001100, false)])});

        XASSERT(create_3dtemplate_rotations(result, two_edge_base_refinement) == 24);

        // Two face Refinement, with 12? orientations
        RefinementTypeByDim<3> two_face_base_refinement(0b11011101, false);
        result.insert(
          {two_face_base_refinement,
           RawCellTemplate()
             .add_entity(
               V{0.0, 0.0, 0.0},
               V{2.0 / 5.0, 2.0 / 5.0, 0.0},
               V{0.0, 1.0, 0.0},
               V{2.0 / 5.0, 3.0 / 5.0, 0.0},
               V{0.0, 0.0, 1.0},
               V{2.0 / 5.0, 2.0 / 5.0, 1.0},
               V{0.0, 1.0, 1.0},
               V{2.0 / 5.0, 3.0 / 5.0, 1.0})
             .add_entity(
               V{2.0 / 5.0, 3.0 / 5.0, 0.0},
               V{3.0 / 5.0, 3.0 / 5.0, 0.0},
               V{0.0, 1.0, 0.0},
               V{1.0, 1.0, 0.0},
               V{2.0 / 5.0, 3.0 / 5.0, 1.0},
               V{3.0 / 5.0, 3.0 / 5.0, 1.0},
               V{0.0, 1.0, 1.0},
               V{1.0, 1.0, 1.0})
             .add_entity(
               V{3.0 / 5.0, 2.0 / 5.0, 0.0},
               V{1.0, 0.0, 0.0},
               V{3.0 / 5.0, 3.0 / 5.0, 0.0},
               V{1.0, 1.0, 0.0},
               V{3.0 / 5.0, 2.0 / 5.0, 1.0},
               V{1.0, 0.0, 1.0},
               V{3.0 / 5.0, 3.0 / 5.0, 1.0},
               V{1.0, 1.0, 1.0})
             .add_entity(
               V{0.0, 0.0, 0.0},
               V{1.0, 0.0, 0.0},
               V{2.0 / 5.0, 2.0 / 5.0, 0.0},
               V{3.0 / 5.0, 2.0 / 5.0, 0.0},
               V{0.0, 0.0, 1.0},
               V{1.0, 0.0, 1.0},
               V{2.0 / 5.0, 2.0 / 5.0, 1.0},
               V{3.0 / 5.0, 2.0 / 5.0, 1.0})
             .axis_aligned(V{2.0 / 5.0, 2.0 / 5.0, 0.0}, V{3.0 / 5.0, 3.0 / 5.0, 1.0})
             .recurse(0, result[RefinementTypeByDim<3>(0b01010101, false)])
             .recurse(0, result[RefinementTypeByDim<3>(0b11001100, false)])
             .recurse(0, result[RefinementTypeByDim<3>(0b10001000, false)])
             .recurse(0, result[RefinementTypeByDim<3>(0b00010001, false)])});

        XASSERT(create_3dtemplate_rotations(result, two_face_base_refinement) == 12);
        RefinementTypeByDim<3> three_face_base_refinement(0b11011111, false);
        result.insert(
          {three_face_base_refinement,
           RawCellTemplate()
             .add_entity( // Left Face
               V{0.0, 0.0, 0.0},
               V{2.0 / 5.0, 2.0 / 5.0, 2.0 / 5.0},
               V{0.0, 1.0, 0.0},
               V{2.0 / 5.0, 3.0 / 5.0, 2.0 / 5.0},
               V{0.0, 0.0, 1.0},
               V{2.0 / 5.0, 2.0 / 5.0, 3.0 / 5.0},
               V{0.0, 1.0, 1.0},
               V{2.0 / 5.0, 3.0 / 5.0, 3.0 / 5.0})
             .add_entity( // Back Face
               V{2.0 / 5.0, 3.0 / 5.0, 2.0 / 5.0},
               V{3.0 / 5.0, 3.0 / 5.0, 2.0 / 5.0},
               V{0.0, 1.0, 0.0},
               V{1.0, 1.0, 0.0},
               V{2.0 / 5.0, 3.0 / 5.0, 3.0 / 5.0},
               V{3.0 / 5.0, 3.0 / 5.0, 3.0 / 5.0},
               V{0.0, 1.0, 1.0},
               V{1.0, 1.0, 1.0})
             .add_entity( // Right Face
               V{3.0 / 5.0, 2.0 / 5.0, 2.0 / 5.0},
               V{1.0, 0.0, 0.0},
               V{3.0 / 5.0, 3.0 / 5.0, 2.0 / 5.0},
               V{1.0, 1.0, 0.0},
               V{3.0 / 5.0, 2.0 / 5.0, 3.0 / 5.0},
               V{1.0, 0.0, 1.0},
               V{3.0 / 5.0, 3.0 / 5.0, 3.0 / 5.0},
               V{1.0, 1.0, 1.0})
             .add_entity( // Front Face
               V{0.0, 0.0, 0.0},
               V{1.0, 0.0, 0.0},
               V{2.0 / 5.0, 2.0 / 5.0, 2.0 / 5.0},
               V{3.0 / 5.0, 2.0 / 5.0, 2.0 / 5.0},
               V{0.0, 0.0, 1.0},
               V{1.0, 0.0, 1.0},
               V{2.0 / 5.0, 2.0 / 5.0, 3.0 / 5.0},
               V{3.0 / 5.0, 2.0 / 5.0, 3.0 / 5.0})
             .add_entity( // Bottom Face
               V{0.0, 0.0, 0.0},
               V{1.0, 0.0, 0.0},
               V{0.0, 1.0, 0.0},
               V{1.0, 1.0, 0.0},
               V{2.0 / 5.0, 2.0 / 5.0, 2.0 / 5.0},
               V{3.0 / 5.0, 2.0 / 5.0, 2.0 / 5.0},
               V{2.0 / 5.0, 3.0 / 5.0, 2.0 / 5.0},
               V{3.0 / 5.0, 3.0 / 5.0, 2.0 / 5.0})
             .add_entity( // Top Face
               V{2.0 / 5.0, 2.0 / 5.0, 3.0 / 5.0},
               V{3.0 / 5.0, 2.0 / 5.0, 3.0 / 5.0},
               V{2.0 / 5.0, 3.0 / 5.0, 3.0 / 5.0},
               V{3.0 / 5.0, 3.0 / 5.0, 3.0 / 5.0},
               V{0.0, 0.0, 1.0},
               V{1.0, 0.0, 1.0},
               V{0.0, 1.0, 1.0},
               V{1.0, 1.0, 1.0})
             .axis_aligned(V{2.0 / 5.0, 2.0 / 5.0, 2.0 / 5.0}, V{3.0 / 5.0, 3.0 / 5.0, 3.0 / 5.0})
             .recurse(0, result[RefinementTypeByDim<3>(0b01010101, false)])
             .recurse(0, result[RefinementTypeByDim<3>(0b11001100, false)])
             .recurse(0, result[RefinementTypeByDim<3>(0b10001010, false)])
             .recurse(0, result[RefinementTypeByDim<3>(0b00010011, false)])
             .recurse(0, result[RefinementTypeByDim<3>(0b00001111, false)])
             .recurse(0, result[RefinementTypeByDim<3>(0b11010000, false)])});

        XASSERT(create_3dtemplate_rotations(result, three_face_base_refinement) == 8);
        // Full refinement
        result.insert(
          {RefinementTypeByDim<3>(0b11111111, false),
           RawCellTemplate().grid({3, 3, 3}, V{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0})});
      }

      return result;
    }
  };
} // namespace FEAT::Geometry
