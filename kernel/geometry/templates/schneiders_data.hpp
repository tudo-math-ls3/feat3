// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/geometry/raw_refinement_templates.hpp>
#include <kernel/geometry/refinement_types.hpp>
#include <kernel/shape.hpp>

namespace FEAT::Geometry
{
  class SchneidersData
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
    using RefinementTypeByDim = StandardRefinementType<Shape::Hypercube<dim_>>;

    /// Accessor for refinement types used by these templates
    template<typename Shape__>
    using RefinementTypeByShape = StandardRefinementType<Shape__>;

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
        {4, 12, 9, 0},
        {8, 36, 54, 27},
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
      using V = typename RawEdgeTemplate::VertexType;

      static EdgeMap result;

      if(result.empty())
      {
        result.insert({RefinementTypeByDim<1>(0b00), RawEdgeTemplate()});

        result.insert(
          {RefinementTypeByDim<1>(0b01),
           RawEdgeTemplate().add_entity(V{0.0}, V{1.0 / 3.0}).add_entity(V{1.0 / 3.0}, V{1.0})});

        result.insert(
          {RefinementTypeByDim<1>(0b10),
           RawEdgeTemplate().add_entity(V{0.0}, V{2.0 / 3.0}).add_entity(V{2.0 / 3.0}, V{1.0})});

        result.insert(
          {RefinementTypeByDim<1>(0b11),
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
      using V = typename RawFaceTemplate::VertexType;

      static FaceMap result;

      if(result.empty())
      {
        // Add base types
        result.insert({RefinementTypeByDim<2>(0b0000), RawFaceTemplate()});

        result.insert(
          {RefinementTypeByDim<2>(0b0001),
           RawFaceTemplate()
             .add_entity(V{0.0, 1.0 / 3.0}, V{1.0 / 3.0, 1.0 / 3.0}, V{0.0, 1.0}, V{1.0, 1.0})
             .add_entity(V{0.0, 0.0}, V{1.0 / 3.0, 0.0}, V{0.0, 1.0 / 3.0}, V{1.0 / 3.0, 1.0 / 3.0})
             .add_entity(V{1.0 / 3.0, 0.0}, V{1.0, 0.0}, V{1.0 / 3.0, 1.0 / 3.0}, V{1.0, 1.0})});

        result.insert(
          {RefinementTypeByDim<2>(0b0011),
           RawFaceTemplate()
             .add_entity(V{1.0 / 3.0, 2.0 / 3.0}, V{2.0 / 3.0, 2.0 / 3.0}, V{0.0, 1.0}, V{1.0, 1.0})
             .add_entity(V{0.0, 1.0 / 3.0}, V{1.0 / 3.0, 1.0 / 3.0}, V{0.0, 1.0}, V{1.0 / 3.0, 2.0 / 3.0})
             .add_entity(
               V{1.0 / 3.0, 1.0 / 3.0},
               V{2.0 / 3.0, 1.0 / 3.0},
               V{1.0 / 3.0, 2.0 / 3.0},
               V{2.0 / 3.0, 2.0 / 3.0})
             .add_entity(V{2.0 / 3.0, 1.0 / 3.0}, V{1.0, 1.0 / 3.0}, V{2.0 / 3.0, 2.0 / 3.0}, V{1.0, 1.0})
             .add_entity(V{0.0, 0.0}, V{1.0 / 3.0, 0.0}, V{0.0, 1.0 / 3.0}, V{1.0 / 3.0, 1.0 / 3.0})
             .add_entity(V{1.0 / 3.0, 0.0}, V{2.0 / 3.0, 0.0}, V{1.0 / 3.0, 1.0 / 3.0}, V{2.0 / 3.0, 1.0 / 3.0})
             .add_entity(V{2.0 / 3.0, 0.0}, V{1.0, 0.0}, V{2.0 / 3.0, 1.0 / 3.0}, V{1.0, 1.0 / 3.0})});

        result.insert(
          {RefinementTypeByDim<2>(0b1001),
           RawFaceTemplate()
             .add_entity(V{1.0 / 3.0, 2.0 / 3.0}, V{2.0 / 3.0, 2.0 / 3.0}, V{0.0, 1.0}, V{2.0 / 3.0, 1.0})
             .add_entity(V{2.0 / 3.0, 2.0 / 3.0}, V{1.0, 2.0 / 3.0}, V{2.0 / 3.0, 1.0}, V{1.0, 1.0})
             .add_entity(V{0.0, 1.0 / 3.0}, V{1.0 / 3.0, 1.0 / 3.0}, V{0.0, 1.0}, V{1.0 / 3.0, 2.0 / 3.0})
             .add_entity(
               V{1.0 / 3.0, 1.0 / 3.0},
               V{2.0 / 3.0, 1.0 / 3.0},
               V{1.0 / 3.0, 2.0 / 3.0},
               V{2.0 / 3.0, 2.0 / 3.0})
             .add_entity(V{2.0 / 3.0, 1.0 / 3.0}, V{1.0, 0.0}, V{2.0 / 3.0, 2.0 / 3.0}, V{1.0, 2.0 / 3.0})
             .add_entity(V{0.0, 0.0}, V{1.0 / 3.0, 0.0}, V{0.0, 1.0 / 3.0}, V{1.0 / 3.0, 1.0 / 3.0})
             .add_entity(V{1.0 / 3.0, 0.0}, V{1.0, 0.0}, V{1.0 / 3.0, 1.0 / 3.0}, V{2.0 / 3.0, 1.0 / 3.0})});
        result.insert(
          {RefinementTypeByDim<2>(0b1011),
           RawFaceTemplate()
             .add_entity(V{1.0 / 3.0, 2.0 / 3.0}, V{2.0 / 3.0, 2.0 / 3.0}, V{0.0, 1.0}, V{2.0 / 3.0, 1.0})
             .add_entity(V{2.0 / 3.0, 2.0 / 3.0}, V{1.0, 2.0 / 3.0}, V{2.0 / 3.0, 1.0}, V{1.0, 1.0})
             .add_entity(V{0.0, 1.0 / 3.0}, V{1.0 / 3.0, 1.0 / 3.0}, V{0.0, 1.0}, V{1.0 / 3.0, 2.0 / 3.0})
             .add_entity(
               V{1.0 / 3.0, 1.0 / 3.0},
               V{2.0 / 3.0, 1.0 / 3.0},
               V{1.0 / 3.0, 2.0 / 3.0},
               V{2.0 / 3.0, 2.0 / 3.0})
             .add_entity(V{2.0 / 3.0, 1.0 / 3.0}, V{1.0, 1.0 / 3.0}, V{2.0 / 3.0, 2.0 / 3.0}, V{1.0, 2.0 / 3.0})
             .add_entity(V{0.0, 0.0}, V{1.0 / 3.0, 0.0}, V{0.0, 1.0 / 3.0}, V{1.0 / 3.0, 1.0 / 3.0})
             .add_entity(V{1.0 / 3.0, 0.0}, V{2.0 / 3.0, 0.0}, V{1.0 / 3.0, 1.0 / 3.0}, V{2.0 / 3.0, 1.0 / 3.0})
             .add_entity(V{2.0 / 3.0, 0.0}, V{1.0, 0.0}, V{2.0 / 3.0, 1.0 / 3.0}, V{1.0, 1.0 / 3.0})});

        result.insert(
          {RefinementTypeByDim<2>(0b1111),
           RawFaceTemplate()
             .add_entity(V{0.0, 2.0 / 3.0}, V{1.0 / 3.0, 2.0 / 3.0}, V{0.0, 1.0}, V{1.0 / 3.0, 1.0})
             .add_entity(V{1.0 / 3.0, 2.0 / 3.0}, V{2.0 / 3.0, 2.0 / 3.0}, V{1.0 / 3.0, 1.0}, V{2.0 / 3.0, 1.0})
             .add_entity(V{2.0 / 3.0, 2.0 / 3.0}, V{1.0, 2.0 / 3.0}, V{2.0 / 3.0, 1.0}, V{1.0, 1.0})
             .add_entity(V{0.0, 1.0 / 3.0}, V{1.0 / 3.0, 1.0 / 3.0}, V{0.0, 2.0 / 3.0}, V{1.0 / 3.0, 2.0 / 3.0})
             .add_entity(
               V{1.0 / 3.0, 1.0 / 3.0},
               V{2.0 / 3.0, 1.0 / 3.0},
               V{1.0 / 3.0, 2.0 / 3.0},
               V{2.0 / 3.0, 2.0 / 3.0})
             .add_entity(V{2.0 / 3.0, 1.0 / 3.0}, V{1.0, 1.0 / 3.0}, V{2.0 / 3.0, 2.0 / 3.0}, V{1.0, 2.0 / 3.0})
             .add_entity(V{0.0, 0.0}, V{1.0 / 3.0, 0.0}, V{0.0, 1.0 / 3.0}, V{1.0 / 3.0, 1.0 / 3.0})
             .add_entity(V{1.0 / 3.0, 0.0}, V{2.0 / 3.0, 0.0}, V{1.0 / 3.0, 1.0 / 3.0}, V{2.0 / 3.0, 1.0 / 3.0})
             .add_entity(V{2.0 / 3.0, 0.0}, V{1.0, 0.0}, V{2.0 / 3.0, 1.0 / 3.0}, V{1.0, 1.0 / 3.0})});

        // Rotate base templates to find remaining types
        result.insert({RefinementTypeByDim<2>(2), rotate_template_2d(result[RefinementTypeByDim<2>(1)])});
        result.insert({RefinementTypeByDim<2>(8), rotate_template_2d(result[RefinementTypeByDim<2>(2)])});
        result.insert({RefinementTypeByDim<2>(4), rotate_template_2d(result[RefinementTypeByDim<2>(8)])});

        result.insert({RefinementTypeByDim<2>(10), rotate_template_2d(result[RefinementTypeByDim<2>(3)])});
        result.insert({RefinementTypeByDim<2>(12), rotate_template_2d(result[RefinementTypeByDim<2>(10)])});
        result.insert({RefinementTypeByDim<2>(5), rotate_template_2d(result[RefinementTypeByDim<2>(12)])});

        result.insert({RefinementTypeByDim<2>(6), rotate_template_2d(result[RefinementTypeByDim<2>(9)])});

        result.insert({RefinementTypeByDim<2>(14), rotate_template_2d(result[RefinementTypeByDim<2>(11)])});
        result.insert({RefinementTypeByDim<2>(13), rotate_template_2d(result[RefinementTypeByDim<2>(14)])});
        result.insert({RefinementTypeByDim<2>(7), rotate_template_2d(result[RefinementTypeByDim<2>(13)])});
      }

      return result;
    }

    /// Returns the raw templates for cells of this template set
    static CellMap& raw_cells()
    {
      using V = typename RawCellTemplate::VertexType;

      static CellMap result;

      if(result.empty())
      {
        // Add base types
        result.insert({RefinementTypeByDim<3>(0b00000000), RawCellTemplate()});

        // Vertex refinement
        result.insert(
          {RefinementTypeByDim<3>(0b00100000),
           RawCellTemplate()
             .add_entity(
               V{2.0 / 3.0, 0.0, 2.0 / 3.0},
               V{1.0, 0.0, 2.0 / 3.0},
               V{2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0},
               V{1.0, 1.0 / 3.0, 2.0 / 3.0},
               V{2.0 / 3.0, 0.0, 1.0},
               V{1.0, 0.0, 1.0},
               V{2.0 / 3.0, 1.0 / 3.0, 1.0},
               V{1.0, 1.0 / 3.0, 1.0})
             .add_entity(
               V{0.0, 0.0, 0.0},
               V{2.0 / 3.0, 0.0, 2.0 / 3.0},
               V{0.0, 1.0, 0.0},
               V{2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0},
               V{0.0, 0.0, 1.0},
               V{2.0 / 3.0, 0.0, 1.0},
               V{0.0, 1.0, 1.0},
               V{2.0 / 3.0, 1.0 / 3.0, 1.0})
             .add_entity(
               V{2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0},
               V{1.0, 1.0 / 3.0, 2.0 / 3.0},
               V{0.0, 1.0, 0.0},
               V{1.0, 1.0, 0.0},
               V{2.0 / 3.0, 1.0 / 3.0, 1.0},
               V{1.0, 1.0 / 3.0, 1.0},
               V{0.0, 1.0, 1.0},
               V{1.0, 1.0, 1.0})
             .add_entity(
               V{0.0, 0.0, 0.0},
               V{1.0, 0.0, 0.0},
               V{0.0, 1.0, 0.0},
               V{1.0, 1.0, 0.0},
               V{2.0 / 3.0, 0.0, 2.0 / 3.0},
               V{1.0, 0.0, 2.0 / 3.0},
               V{2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0},
               V{1.0, 1.0 / 3.0, 2.0 / 3.0})});

        // Edge refinement
        result.insert(
          {RefinementTypeByDim<3>(0b00110000),
           RawCellTemplate()
             .add_entity(
               V{0.0, 0.0, 2.0 / 3.0},
               V{1.0 / 3.0, 0.0, 2.0 / 3.0},
               V{0.0, 1.0 / 3.0, 2.0 / 3.0},
               V{1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0},
               V{0.0, 0.0, 1.0},
               V{1.0 / 3.0, 0.0, 1.0},
               V{0.0, 1.0 / 3.0, 1.0},
               V{1.0 / 3.0, 1.0 / 3.0, 1.0})
             .add_entity(
               V{1.0 / 3.0, 0.0, 2.0 / 3.0},
               V{2.0 / 3.0, 0.0, 2.0 / 3.0},
               V{1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0},
               V{2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0},
               V{1.0 / 3.0, 0.0, 1.0},
               V{2.0 / 3.0, 0.0, 1.0},
               V{1.0 / 3.0, 1.0 / 3.0, 1.0},
               V{2.0 / 3.0, 1.0 / 3.0, 1.0})
             .add_entity(
               V{2.0 / 3.0, 0.0, 2.0 / 3.0},
               V{1.0, 0.0, 2.0 / 3.0},
               V{2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0},
               V{1.0, 1.0 / 3.0, 2.0 / 3.0},
               V{2.0 / 3.0, 0.0, 1.0},
               V{1.0, 0.0, 1.0},
               V{2.0 / 3.0, 1.0 / 3.0, 1.0},
               V{1.0, 1.0 / 3.0, 1.0})
             .add_entity(
               V{0.0, 0.0, 0.0},
               V{1.0 / 3.0, 0.0, 1.0 / 3.0},
               V{0.0, 1.0, 0.0},
               V{1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0},
               V{0.0, 0.0, 2.0 / 3.0},
               V{1.0 / 3.0, 0.0, 2.0 / 3.0},
               V{0.0, 1.0 / 3.0, 2.0 / 3.0},
               V{1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0})
             .add_entity(
               V{1.0 / 3.0, 0.0, 1.0 / 3.0},
               V{2.0 / 3.0, 0.0, 1.0 / 3.0},
               V{1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0},
               V{2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0},
               V{1.0 / 3.0, 0.0, 2.0 / 3.0},
               V{2.0 / 3.0, 0.0, 2.0 / 3.0},
               V{1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0},
               V{2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0})
             .add_entity(
               V{2.0 / 3.0, 0.0, 1.0 / 3.0},
               V{1.0, 0.0, 0.0},
               V{2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0},
               V{1.0, 1.0, 0.0},
               V{2.0 / 3.0, 0.0, 2.0 / 3.0},
               V{1.0, 0.0, 2.0 / 3.0},
               V{2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0},
               V{1.0, 1.0 / 3.0, 2.0 / 3.0})
             .add_entity(
               V{0.0, 1.0 / 3.0, 2.0 / 3.0},
               V{1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0},
               V{0.0, 1.0, 0.0},
               V{1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0},
               V{0.0, 1.0 / 3.0, 1.0},
               V{1.0 / 3.0, 1.0 / 3.0, 1.0},
               V{0.0, 1.0, 1.0},
               V{1.0 / 3.0, 2.0 / 3.0, 1.0})
             .add_entity(
               V{1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0},
               V{2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0},
               V{1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0},
               V{2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0},
               V{1.0 / 3.0, 1.0 / 3.0, 1.0},
               V{2.0 / 3.0, 1.0 / 3.0, 1.0},
               V{1.0 / 3.0, 2.0 / 3.0, 1.0},
               V{2.0 / 3.0, 2.0 / 3.0, 1.0})
             .add_entity(
               V{2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0},
               V{1.0, 1.0 / 3.0, 2.0 / 3.0},
               V{2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0},
               V{1.0, 1.0, 0.0},
               V{2.0 / 3.0, 1.0 / 3.0, 1.0},
               V{1.0, 1.0 / 3.0, 1.0},
               V{2.0 / 3.0, 2.0 / 3.0, 1.0},
               V{1.0, 1.0, 1.0})
             .add_entity(
               V{0.0, 0.0, 0.0},
               V{1.0, 0.0, 0.0},
               V{0.0, 1.0, 0.0},
               V{1.0, 1.0, 0.0},
               V{1.0 / 3.0, 0.0, 1.0 / 3.0},
               V{2.0 / 3.0, 0.0, 1.0 / 3.0},
               V{1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0},
               V{2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0})
             .add_entity(
               V{1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0},
               V{2.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0},
               V{0.0, 1.0, 0.0},
               V{1.0, 1.0, 0.0},
               V{1.0 / 3.0, 2.0 / 3.0, 1.0},
               V{2.0 / 3.0, 2.0 / 3.0, 1.0},
               V{0.0, 1.0, 1.0},
               V{1.0, 1.0, 1.0})});

        // Face refinement
        result.insert(
          {RefinementTypeByDim<3>(0b11110000),
           RawCellTemplate()
             .add_entity(
               V{0.0, 0.0, 2.0 / 3.0},
               V{1.0 / 3.0, 0.0, 2.0 / 3.0},
               V{0.0, 1.0 / 3.0, 2.0 / 3.0},
               V{1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0},
               V{0.0, 0.0, 1.0},
               V{1.0 / 3.0, 0.0, 1.0},
               V{0.0, 1.0 / 3.0, 1.0},
               V{1.0 / 3.0, 1.0 / 3.0, 1.0})
             .add_entity(
               V{1.0 / 3.0, 0.0, 2.0 / 3.0},
               V{2.0 / 3.0, 0.0, 2.0 / 3.0},
               V{1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0},
               V{2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0},
               V{1.0 / 3.0, 0.0, 1.0},
               V{2.0 / 3.0, 0.0, 1.0},
               V{1.0 / 3.0, 1.0 / 3.0, 1.0},
               V{2.0 / 3.0, 1.0 / 3.0, 1.0})
             .add_entity(
               V{2.0 / 3.0, 0.0, 2.0 / 3.0},
               V{1.0, 0.0, 2.0 / 3.0},
               V{2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0},
               V{1.0, 1.0 / 3.0, 2.0 / 3.0},
               V{2.0 / 3.0, 0.0, 1.0},
               V{1.0, 0.0, 1.0},
               V{2.0 / 3.0, 1.0 / 3.0, 1.0},
               V{1.0, 1.0 / 3.0, 1.0})
             .add_entity(
               V{0.0, 1.0 / 3.0, 2.0 / 3.0},
               V{1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0},
               V{0.0, 2.0 / 3.0, 2.0 / 3.0},
               V{1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0},
               V{0.0, 1.0 / 3.0, 1.0},
               V{1.0 / 3.0, 1.0 / 3.0, 1.0},
               V{0.0, 2.0 / 3.0, 1.0},
               V{1.0 / 3.0, 2.0 / 3.0, 1.0})
             .add_entity(
               V{1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0},
               V{2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0},
               V{1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0},
               V{2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0},
               V{1.0 / 3.0, 1.0 / 3.0, 1.0},
               V{2.0 / 3.0, 1.0 / 3.0, 1.0},
               V{1.0 / 3.0, 2.0 / 3.0, 1.0},
               V{2.0 / 3.0, 2.0 / 3.0, 1.0})
             .add_entity(
               V{2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0},
               V{1.0, 1.0 / 3.0, 2.0 / 3.0},
               V{2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0},
               V{1.0, 2.0 / 3.0, 2.0 / 3.0},
               V{2.0 / 3.0, 1.0 / 3.0, 1.0},
               V{1.0, 1.0 / 3.0, 1.0},
               V{2.0 / 3.0, 2.0 / 3.0, 1.0},
               V{1.0, 2.0 / 3.0, 1.0})
             .add_entity(
               V{0.0, 2.0 / 3.0, 2.0 / 3.0},
               V{1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0},
               V{0.0, 1.0, 2.0 / 3.0},
               V{1.0 / 3.0, 1.0, 2.0 / 3.0},
               V{0.0, 2.0 / 3.0, 1.0},
               V{1.0 / 3.0, 2.0 / 3.0, 1.0},
               V{0.0, 1.0, 1.0},
               V{1.0 / 3.0, 1.0, 1.0})
             .add_entity(
               V{1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0},
               V{2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0},
               V{1.0 / 3.0, 1.0, 2.0 / 3.0},
               V{2.0 / 3.0, 1.0, 2.0 / 3.0},
               V{1.0 / 3.0, 2.0 / 3.0, 1.0},
               V{2.0 / 3.0, 2.0 / 3.0, 1.0},
               V{1.0 / 3.0, 1.0, 1.0},
               V{2.0 / 3.0, 1.0, 1.0})
             .add_entity(
               V{2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0},
               V{1.0, 2.0 / 3.0, 2.0 / 3.0},
               V{2.0 / 3.0, 1.0, 2.0 / 3.0},
               V{1.0, 1.0, 2.0 / 3.0},
               V{2.0 / 3.0, 2.0 / 3.0, 1.0},
               V{1.0, 2.0 / 3.0, 1.0},
               V{2.0 / 3.0, 1.0, 1.0},
               V{1.0, 1.0, 1.0})
             .add_entity(
               V{0.0, 0.0, 0.0},
               V{1.0 / 3.0, 0.0, 1.0 / 3.0},
               V{0.0, 1.0 / 3.0, 1.0 / 3.0},
               V{1.0 / 3.0, 1.0 / 3.0, 1.0 / 2.0},
               V{0.0, 0.0, 2.0 / 3.0},
               V{1.0 / 3.0, 0.0, 2.0 / 3.0},
               V{0.0, 1.0 / 3.0, 2.0 / 3.0},
               V{1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0})
             .add_entity(
               V{1.0 / 3.0, 0.0, 1.0 / 3.0},
               V{2.0 / 3.0, 0.0, 1.0 / 3.0},
               V{1.0 / 3.0, 1.0 / 3.0, 1.0 / 2.0},
               V{2.0 / 3.0, 1.0 / 3.0, 1.0 / 2.0},
               V{1.0 / 3.0, 0.0, 2.0 / 3.0},
               V{2.0 / 3.0, 0.0, 2.0 / 3.0},
               V{1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0},
               V{2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0})
             .add_entity(
               V{2.0 / 3.0, 0.0, 1.0 / 3.0},
               V{1.0, 0.0, 0.0},
               V{2.0 / 3.0, 1.0 / 3.0, 1.0 / 2.0},
               V{1.0, 1.0 / 3.0, 1.0 / 3.0},
               V{2.0 / 3.0, 0.0, 2.0 / 3.0},
               V{1.0, 0.0, 2.0 / 3.0},
               V{2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0},
               V{1.0, 1.0 / 3.0, 2.0 / 3.0})
             .add_entity(
               V{2.0 / 3.0, 1.0 / 3.0, 1.0 / 2.0},
               V{1.0, 1.0 / 3.0, 1.0 / 3.0},
               V{2.0 / 3.0, 2.0 / 3.0, 1.0 / 2.0},
               V{1.0, 2.0 / 3.0, 1.0 / 3.0},
               V{2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0},
               V{1.0, 1.0 / 3.0, 2.0 / 3.0},
               V{2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0},
               V{1.0, 2.0 / 3.0, 2.0 / 3.0})
             .add_entity(
               V{2.0 / 3.0, 2.0 / 3.0, 1.0 / 2.0},
               V{1.0, 2.0 / 3.0, 1.0 / 3.0},
               V{2.0 / 3.0, 1.0, 1.0 / 3.0},
               V{1.0, 1.0, 0.0},
               V{2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0},
               V{1.0, 2.0 / 3.0, 2.0 / 3.0},
               V{2.0 / 3.0, 1.0, 2.0 / 3.0},
               V{1.0, 1.0, 2.0 / 3.0})
             .add_entity(
               V{1.0 / 3.0, 2.0 / 3.0, 1.0 / 2.0},
               V{2.0 / 3.0, 2.0 / 3.0, 1.0 / 2.0},
               V{1.0 / 3.0, 1.0, 1.0 / 3.0},
               V{2.0 / 3.0, 1.0, 1.0 / 3.0},
               V{1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0},
               V{2.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0},
               V{1.0 / 3.0, 1.0, 2.0 / 3.0},
               V{2.0 / 3.0, 1.0, 2.0 / 3.0})
             .add_entity(
               V{0.0, 2.0 / 3.0, 1.0 / 3.0},
               V{1.0 / 3.0, 2.0 / 3.0, 1.0 / 2.0},
               V{0.0, 1.0, 0.0},
               V{1.0 / 3.0, 1.0, 1.0 / 3.0},
               V{0.0, 2.0 / 3.0, 2.0 / 3.0},
               V{1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0},
               V{0.0, 1.0, 2.0 / 3.0},
               V{1.0 / 3.0, 1.0, 2.0 / 3.0})
             .add_entity(
               V{0.0, 1.0 / 3.0, 1.0 / 3.0},
               V{1.0 / 3.0, 1.0 / 3.0, 1.0 / 2.0},
               V{0.0, 2.0 / 3.0, 1.0 / 3.0},
               V{1.0 / 3.0, 2.0 / 3.0, 1.0 / 2.0},
               V{0.0, 1.0 / 3.0, 2.0 / 3.0},
               V{1.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0},
               V{0.0, 2.0 / 3.0, 2.0 / 3.0},
               V{1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0})
             .add_entity(
               V{0.0, 0.0, 0.0},
               V{1.0, 0.0, 0.0},
               V{0.0, 1.0, 0.0},
               V{1.0, 1.0, 0.0},
               V{1.0 / 3.0, 0.0, 1.0 / 3.0},
               V{2.0 / 3.0, 0.0, 1.0 / 3.0},
               V{1.0 / 3.0, 1.0, 1.0 / 3.0},
               V{2.0 / 3.0, 1.0, 1.0 / 3.0})
             .add_entity(
               V{0.0, 0.0, 0.0},
               V{1.0 / 3.0, 0.0, 1.0 / 3.0},
               V{0.0, 1.0, 0.0},
               V{1.0 / 3.0, 1.0, 1.0 / 3.0},
               V{0.0, 1.0 / 3.0, 1.0 / 3.0},
               V{1.0 / 3.0, 1.0 / 3.0, 1.0 / 2.0},
               V{0.0, 2.0 / 3.0, 1.0 / 3.0},
               V{1.0 / 3.0, 2.0 / 3.0, 1.0 / 2.0})
             .add_entity(
               V{2.0 / 3.0, 0.0, 1.0 / 3.0},
               V{1.0, 0.0, 0.0},
               V{2.0 / 3.0, 1.0, 1.0 / 3.0},
               V{1.0, 1.0, 0.0},
               V{2.0 / 3.0, 1.0 / 3.0, 1.0 / 2.0},
               V{1.0, 1.0 / 3.0, 1.0 / 3.0},
               V{2.0 / 3.0, 2.0 / 3.0, 1.0 / 2.0},
               V{1.0, 2.0 / 3.0, 1.0 / 3.0})
             .add_entity(
               V{1.0 / 3.0, 0.0, 1.0 / 3.0},
               V{2.0 / 3.0, 0.0, 1.0 / 3.0},
               V{1.0 / 3.0, 1.0, 1.0 / 3.0},
               V{2.0 / 3.0, 1.0, 1.0 / 3.0},
               V{1.0 / 3.0, 1.0 / 3.0, 1.0 / 2.0},
               V{2.0 / 3.0, 1.0 / 3.0, 1.0 / 2.0},
               V{1.0 / 3.0, 2.0 / 3.0, 1.0 / 2.0},
               V{2.0 / 3.0, 2.0 / 3.0, 1.0 / 2.0})});

        // Cell refinement
        result.insert(
          {RefinementTypeByDim<3>(0b11111111), RawCellTemplate().grid({3, 3, 3}, V{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0})});

        result[RefinementTypeByDim<3>(128)] = rotate_template_zaxis(result[RefinementTypeByDim<3>(32)]);
        result[RefinementTypeByDim<3>(64)] = rotate_template_zaxis(result[RefinementTypeByDim<3>(128)]);
        result[RefinementTypeByDim<3>(16)] = rotate_template_zaxis(result[RefinementTypeByDim<3>(64)]);

        result[RefinementTypeByDim<3>(2)] = rotate_template_yaxis(result[RefinementTypeByDim<3>(32)]);
        result[RefinementTypeByDim<3>(8)] = rotate_template_zaxis(result[RefinementTypeByDim<3>(2)]);
        result[RefinementTypeByDim<3>(4)] = rotate_template_zaxis(result[RefinementTypeByDim<3>(8)]);
        result[RefinementTypeByDim<3>(1)] = rotate_template_zaxis(result[RefinementTypeByDim<3>(4)]);

        result[RefinementTypeByDim<3>(34)] = rotate_template_yaxis(result[RefinementTypeByDim<3>(48)]);
        result[RefinementTypeByDim<3>(3)] = rotate_template_yaxis(result[RefinementTypeByDim<3>(34)]);
        result[RefinementTypeByDim<3>(17)] = rotate_template_yaxis(result[RefinementTypeByDim<3>(3)]);
        result[RefinementTypeByDim<3>(160)] = rotate_template_zaxis(result[RefinementTypeByDim<3>(48)]);
        result[RefinementTypeByDim<3>(10)] = rotate_template_yaxis(result[RefinementTypeByDim<3>(160)]);
        result[RefinementTypeByDim<3>(5)] = rotate_template_yaxis(result[RefinementTypeByDim<3>(10)]);
        result[RefinementTypeByDim<3>(80)] = rotate_template_yaxis(result[RefinementTypeByDim<3>(5)]);
        auto tmp = rotate_template_zaxis(result[RefinementTypeByDim<3>(48)]);
        result[RefinementTypeByDim<3>(192)] = rotate_template_zaxis(tmp);
        result[RefinementTypeByDim<3>(136)] = rotate_template_yaxis(result[RefinementTypeByDim<3>(192)]);
        result[RefinementTypeByDim<3>(12)] = rotate_template_yaxis(result[RefinementTypeByDim<3>(136)]);
        result[RefinementTypeByDim<3>(68)] = rotate_template_yaxis(result[RefinementTypeByDim<3>(12)]);

        result[RefinementTypeByDim<3>(51)] = rotate_template_xaxis(result[RefinementTypeByDim<3>(240)]);
        result[RefinementTypeByDim<3>(15)] = rotate_template_xaxis(result[RefinementTypeByDim<3>(51)]);
        result[RefinementTypeByDim<3>(204)] = rotate_template_xaxis(result[RefinementTypeByDim<3>(15)]);
        result[RefinementTypeByDim<3>(85)] = rotate_template_yaxis(result[RefinementTypeByDim<3>(15)]);
        result[RefinementTypeByDim<3>(170)] = rotate_template_yaxis(result[RefinementTypeByDim<3>(240)]);
      }

      return result;
    }
  };
} // namespace FEAT::Geometry
