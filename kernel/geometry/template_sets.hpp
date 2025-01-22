// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/geometry/refinement_types.hpp>
#include <kernel/shape.hpp>
#include <kernel/geometry/template_builder.hpp>
#include <kernel/geometry/subdivision_levels.hpp>
#include <kernel/geometry/intern/refinement_field.hpp>

namespace FEAT::Geometry
{
  namespace Intern
  {
    template<typename Templates_, typename MeshType_>
    static RefinementField<std::uint64_t>
    make_standard_refinement_field(const MeshType_& mesh, Templates_& templates, const Geometry::SubdivisionLevels& sdls)
    {
      using ShapeType = typename MeshType_::ShapeType;
      static const constexpr int dim = MeshType_::shape_dim;
      static const constexpr int num_vertices = Shape::FaceTraits<ShapeType, 0>::count;
      //static const constexpr int num_neighbors = Shape::FaceTraits<ShapeType, dim - 1>::count;

      // Copy SDLs into refinement field
      RefinementField<std::uint64_t> result(sdls.size());
      for(Index i(0); i < sdls.size(); ++i)
      {
        result[i] = sdls[i];
      }

      // We need to ensure the refinement field contains subdivision levels
      // such that all cells (and their (grand-)children) have valid refinement
      // types, for which templates exist. We treat the markings for a cells
      // vertices as a series of refinement types and ensure that each type is
      // valid. This is equivalent to assuming that each of these types will occur
      // in some (grand-)child at least once. That is not necessarily the case, as
      // refinement levels get spread during mesh adaptation, but is a safe
      // assumption. As an example, consider a cell where two diagonally opposite
      // vertices are marked. These vertices will not be shared by any child cell,
      // and their levels could be adjusted individually after the first refinement.

      // Grab necessary mesh and template info
      const auto& v_at_c = mesh.template get_index_set<dim, 0>();
      //const auto& neighbors = mesh.get_neighbors();

      // The worklist indicates which cells still need to be checked. At the
      // beginning all cells need to be checked. If a cells refinement types are
      // all valid, it gets removed from the worklist. If a cells vertex markings
      // needed to be adjusted, its neighboring cells get added to the worklist, as
      // modifying the vertex markings changed those cells types as well.
      std::vector<bool> worklist(mesh.get_num_elements(), true);

      while(true)
      {
        // Find cell to check
        auto it = std::find(worklist.begin(), worklist.end(), true);

        // No remaining marked cells in worklist. We are done.
        if(it == worklist.end())
        {
          break;
        }

        Index cell = std::distance(worklist.begin(), it);

        RefinementFieldTuple<std::uint64_t, num_vertices> levels =
          result.get_tuple(v_at_c[cell]);

        StandardRefinementType<ShapeType> type(levels);

        while(true)
        {
          // Strip away levels until type is no longer valid or type is zero type.
          while(!type.is_zero_refinement() && templates.template has_template<dim>(type))
          {
            // Type is valid. Decrement markings to check next type.
            for(Index i(0); i < num_vertices; i++)
            {
              levels[i] = levels[i] > 0 ? levels[i] - 1 : 0;
            }
            type = StandardRefinementType<ShapeType>(levels);
          }

          if(type.is_zero_refinement())
          {
            // No markings left. All types of the cell are valid.
            worklist[cell] = false;
            break;
          }

          // Invalid type found. Adjust levels
          auto& adjustment = templates.template type_adjustment<dim>(type.to_number());
          for(int i(0); i < num_vertices; i++)
          {
            result[v_at_c[cell][i]] += adjustment[i];
            levels[i] += adjustment[i];
          }
          type = StandardRefinementType<ShapeType>(levels);

          // We changed the cells vertex markings. This also changed the markings
          // of its neighbors. We need to re-check those.
          /*for(Index i(0); i < num_neighbors; i++)
          {
            if(neighbors(cell, i) != ~Index(0))
            {
              worklist[neighbors(cell, i)] = true;
            }
          }*/
          for(Index idx = 0; idx < mesh.get_num_elements(); idx++)
          {
            worklist[idx] = true;
          }
          // Continue checking
        }
      }
      return result;
    }

    template<typename Templates_, typename MeshType_>
    static RefinementField<IsolatedPointVertexMarking>
    make_isolated_point_refinement_field(const MeshType_& mesh, Templates_& templates, const Geometry::SubdivisionLevels& sdls)
    {
      using ShapeType = typename MeshType_::ShapeType;
      static const constexpr int dim = MeshType_::shape_dim;
      static const constexpr int num_vertices = Shape::FaceTraits<ShapeType, 0>::count;
      //static const constexpr int num_neighbors = Shape::FaceTraits<ShapeType, dim - 1>::count;

      // Copy SDLs into refinement field. We assume that all marked vertices
      // are isolated.
      RefinementField<IsolatedPointVertexMarking> result(sdls.size());
      for(Index i(0); i < sdls.size(); ++i)
      {
        result[i] = IsolatedPointVertexMarking{sdls[i], sdls[i] > 0};
      }

      // We need to ensure the refinement field contains markings such that all
      // cells (and their (grand-)children) have valid refinement types, for which
      // templates exist. We treat the markings for a cells vertices as a series of
      // refinement types and ensure that each type is valid. This is equivalent to
      // assuming that each of these types will occur in some (grand-)child at
      // least once. That is not necessarily the case, as refinement levels get
      // spread during mesh adaptation, but is a safe assumption. As an example,
      // consider a cell where two diagonally opposite vertices are marked. These
      // vertices will not be shared by any child cell, and their levels could be
      // adjusted individually after the first refinement.

      // Grab necessary mesh and template info
      const auto& v_at_c = mesh.template get_index_set<dim, 0>();
      //const auto& neighbors = mesh.get_neighbors();

      // The worklist indicates which cells still need to be checked. At the
      // beginning all cells need to be checked. If a cells refinement types are
      // all valid, it gets removed from the worklist. If a cells vertex markings
      // needed to be adjusted, its neighboring cells get added to the worklist, as
      // modifying the vertex markings changed those cells types as well.
      std::vector<bool> worklist(mesh.get_num_elements(), true);

      // Build accurate information about isolated points
      auto determine_isolation = [&]()
      {
        // A vertex is not isolated, if it is marked together with another
        // vertex in at least one cell
        for(Index cell = 0; cell < mesh.get_num_elements(); cell++)
        {
          Index num_marked = 0;
          for(int vertex(0); vertex < v_at_c.num_indices; vertex++)
          {
            num_marked += result[v_at_c(cell, vertex)].level > 0 ? 1 : 0;
          }

          if(num_marked != 1)
          {
            for(int vertex(0); vertex < v_at_c.num_indices; vertex++)
            {
              result[v_at_c(cell, vertex)].is_isolated = false;
            }
          }
        }
      };


      while(true)
      {
        // (Re-)determine isolation of vertices. Adjustments might have caused
        // previously isolated vertices to no longer be isolated.
        determine_isolation();
        // Find cell to check
        auto it = std::find(worklist.begin(), worklist.end(), true);

        // No remaining marked cells in worklist. We are done.
        if(it == worklist.end())
        {
          break;
        }

        Index cell = std::distance(worklist.begin(), it);

        RefinementFieldTuple<IsolatedPointVertexMarking, num_vertices> levels =
          result.get_tuple(v_at_c[cell]);

        IsolatedPointRefinementType<ShapeType> type(levels);

        while(true)
        {
          // Strip away levels until type is no longer valid or type is zero type.
          while(!type.is_zero_refinement() && templates.template has_template<dim>(type))
          {
            // Type is valid. Decrement markings to check next type.
            for(Index i(0); i < num_vertices; i++)
            {
              levels[i].level = levels[i].level > 0 ? levels[i].level - 1 : 0;
            }
            type = IsolatedPointRefinementType<ShapeType>(levels);
          }

          if(type.is_zero_refinement())
          {
            // No markings left. All types of the cell are valid.
            worklist[cell] = false;
            break;
          }

          // Invalid type found. Adjust levels
          // NOTE(mmuegge): We do not need to consider isolation here.
          // levels contains either:
          // * no marked vertices, in which case we do not hit this code path
          // * a single isolated marked vertex, in which case all possible
          //   templates exist and we do not hit this code path
          // * multiple marked vertices, in which case they aren't isolated per
          //   definition, and adding additional markings will not cause a previously
          //   isolated vertex (of this tuple) to become unisolated.
          // Also note that this does not preclude some other vertex from no
          // longer being isolated. We handle this by re-determining isolation at
          // the start of checking each cell.
          auto& adjustment = templates.template type_adjustment<dim>(type.to_number());
          for(int i(0); i < num_vertices; i++)
          {
            result[v_at_c[cell][i]].level += adjustment[i];
            levels[i].level += adjustment[i];
          }
          type = IsolatedPointRefinementType<ShapeType>(levels);

          // We changed the cells vertex markings. This also changed the markings
          // of its neighbors. We need to re-check those.
          for(Index idx = 0; idx < mesh.get_num_elements(); idx++)
          {
            worklist[idx] = true;
          }
          /*for(Index i(0); i < num_neighbors; i++)
          {
            if(neighbors(cell, i) != ~Index(0))
            {
              worklist[neighbors(cell, i)] = true;
            }
          }*/
          // Continue checking
        }
      }
      return result;
    }

    template<typename Shape_>
    std::uint64_t
    spread_refinement_field(
        const EntityReference& ref,
        const Geometry::Intern::RefinementFieldTuple<std::uint64_t, Shape::FaceTraits<Shape_, 0>::count>& markings)
    {
      switch(ref.source)
      {
        case EntitySource::ParentTopology: return markings[ref.index] > 0 ? markings[ref.index] - 1 : 0;
        case EntitySource::Sibling:
        {
          std::uint64_t min = markings[0];
          for(Index i(1); i < markings.size; ++i)
          {
            min = std::min(min, markings[i]);
          }
          return min > 0 ? min - 1 : 0;
        }
        case EntitySource::BoundaryEdge:
        {
          if constexpr (Shape_::dimension >= 2)
          {
            using Mapping = Intern::FaceIndexMapping<Shape_, 1, 0>;

            int face = ref.entity;

            std::uint64_t min = markings[Mapping::map(face, 0)];
            for(int i(1); i < 2; ++i)
            {
              min = std::min(min, markings[Mapping::map(face, i)]);
            }
            return min > 0 ? min - 1 : 0;
          }
          return 0;
        }
        case EntitySource::BoundaryFace:
        {
          if constexpr (Shape_::dimension >= 3)
          {
            using Mapping = Intern::FaceIndexMapping<Shape_, 2, 0>;

            int face = ref.entity;

            std::uint64_t min = markings[Mapping::map(face, 0)];
            for(int i(1); i < 4; ++i)
            {
              min = std::min(min, markings[Mapping::map(face, i)]);
            }
            return min > 0 ? min - 1 : 0;
          }
          return 0;
        }
        default:
          return 0;
      }
    }

    template<typename Shape_>
    IsolatedPointVertexMarking
    spread_refinement_field(
        const EntityReference& ref,
        const Geometry::Intern::RefinementFieldTuple<IsolatedPointVertexMarking, Shape::FaceTraits<Shape_, 0>::count>& markings)
    {
      switch(ref.source)
      {
        case EntitySource::ParentTopology: return {
          markings[ref.index].level > 0 ? markings[ref.index].level - 1 : 0,
          markings[ref.index].is_isolated
        };
        case EntitySource::Sibling:
        {
          std::uint64_t min = markings[0].level;
          for(Index i(1); i < markings.size; ++i)
          {
            min = std::min(min, markings[i].level);
          }
          return {min > 0 ? min - 1 : 0, false};
        }
        case EntitySource::BoundaryEdge:
        {
          if constexpr (Shape_::dimension >= 2)
          {
            using Mapping = Intern::FaceIndexMapping<Shape_, 1, 0>;

            int face = ref.entity;

            std::uint64_t min = markings[Mapping::map(face, 0)].level;
            for(int i(1); i < 2; ++i)
            {
              min = std::min(min, markings[Mapping::map(face, i)].level);
            }
            return {min > 0 ? min - 1 : 0, false};
          }
          return {0, false};
        }
        case EntitySource::BoundaryFace:
        {
          if constexpr (Shape_::dimension >= 3)
          {
            using Mapping = Intern::FaceIndexMapping<Shape_, 2, 0>;

            int face = ref.entity;

            std::uint64_t min = markings[Mapping::map(face, 0)].level;
            for(int i(1); i < 4; ++i)
            {
              min = std::min(min, markings[Mapping::map(face, i)].level);
            }
            return {min > 0 ? min - 1 : 0, false};
          }
          return {0, false};
        }
        default:
          return {0, false};
      }
    }
  } // namespace Intern

  /**
   * \brief Schneiders template set for adaptive mesh refinement
   *
   * This class can be used as a TemplateSet for the AdaptiveMesh class. It supports Quadrilateral and Hexahedral
   * meshes. The templates are described in R. Schneiders, Refining quadrilateral and hexahedral element meshes, 1998,
   *
   * This template set defines a full set of templates for refinement of
   * quadrilateral meshes. For hexahedral meshes the templates are limited and
   * allow only refinement of convex areas.
   */
  template<typename RawData_>
  class StandardTemplateSet
  {
    inline static TemplateBuilder<RawData_> _templates = {};

    using MaxShape = typename RawData_::MaxShape;
  public:
    using VertexMarkerType = std::uint64_t;

    template<int dim_>
    using RefinementTypeByDim = StandardRefinementType<typename Shape::FaceTraits<MaxShape, dim_>::ShapeType>;

    template<typename Shape_>
    using RefinementTypeByShape = IsolatedPointRefinementType<Shape_>;

    /**
     * \brief Shape compatability test
     *
     * \tparam Shape_ Mesh shape this template set is to be applied to
     *
     * \returns True, if the template set is compatible with the mesh shape, false otherwise.
     */
    template<typename Shape_>
    static constexpr bool is_shape_compatible()
    {
      return RawData_::template is_shape_compatible<Shape_>();
    }

    static void stats()
    {
      _templates().stats();
    }

    /**
     * \brief Returns maximum number of children a template produces
     *
     * \tparam template_dim_ Dimension of template
     * \tparam child_dim_ Dimension of children
     */
    template<int template_dim_, int child_dim_>
    static constexpr int max_children()
    {
      return RawData_::template max_children<template_dim_, child_dim_>();
    }

    /**
     * \brief Adjusts SubdivisionLevels to match TemplateSet requirements
     *
     * Will only add subdivision levels to match requirements. No subdivision levels will be removed.
     * No-op for quadrilateral meshes, as all 2D templates exists and all markings are valid.
     *
     * \param[inout] sdls SubdivisionLevels to adjust
     */
    template<typename Mesh_>
    static Intern::RefinementField<std::uint64_t> make_refinement_field(const Mesh_& mesh, Geometry::SubdivisionLevels& sdls)
    {
      return Intern::make_standard_refinement_field(mesh, _templates, sdls);
    }

    template<typename Shape_>
    static VertexMarkerType spread_refinement_field(
        const EntityReference& ref,
        const Intern::RefinementFieldTuple<VertexMarkerType, Shape::FaceTraits<Shape_, 0>::count>& levels)
    {
      return Intern::spread_refinement_field<Shape_>(ref, levels);
    }

    template<typename Shape_>
    static const RefinementTemplate<Shape_>& get_template(StandardRefinementType<Shape_> type)
    {
      return _templates.template get_template<Shape_::dimension>(type);
    }

    template<typename Shape_>
    static bool has_template(StandardRefinementType<Shape_> type)
    {
      return _templates.template has_template<Shape_::dimension>(type);
    }

    template<typename Shape_, int dim_>
    static std::pair<Index, int> correct_for_orientation(StandardRefinementType<Shape_> type, int orientation, Index idx)
    {
      return _templates.template correct_for_orientation<Shape_::dimension, dim_>(type, orientation, idx);
    }

    /**
     * \brief Returns number of children a specific template produces
     *
     * This is intended as a utility function for iterating over children without creating full TemplateQuerys.
     *
     * \tparam template_dim_ Template dimension
     * \tparam child_dim_ Child dimension
     *
     * \param[in] type Type of template
     *
     * \returns Number of children of dimension \c child_dim_ of template with dimension \c template_dim_ and type \c
     * type
     */
    template<int template_dim_, int child_dim_>
    static Index num_children(StandardRefinementType<typename Shape::FaceTraits<MaxShape, template_dim_>::ShapeType> type)
    {
      return _templates.template get_template<template_dim_>(type).template num_entities<child_dim_>();
    }
  };

  /**
   * \brief Schneiders template set for adaptive mesh refinement
   *
   * This class can be used as a TemplateSet for the AdaptiveMesh class. It supports Quadrilateral and Hexahedral
   * meshes. The templates are described in R. Schneiders, Refining quadrilateral and hexahedral element meshes, 1998,
   *
   * This template set defines a full set of templates for refinement of
   * quadrilateral meshes. For hexahedral meshes the templates are limited and
   * allow only refinement of convex areas.
   */
  template<typename RawData_>
  class IsolatedPointTemplateSet
  {
    using MaxShape = typename RawData_::MaxShape;
  public:
    using VertexMarkerType = IsolatedPointVertexMarking;

    template<int dim_>
    using RefinementTypeByDim = IsolatedPointRefinementType<typename Shape::FaceTraits<MaxShape, dim_>::ShapeType>;

    template<typename Shape_>
    using RefinementTypeByShape = IsolatedPointRefinementType<Shape_>;

    /**
     * \brief Shape compatability test
     *
     * \tparam Shape_ Mesh shape this template set is to be applied to
     *
     * \returns True, if the template set is compatible with the mesh shape, false otherwise.
     */
    template<typename Shape_>
    static constexpr bool is_shape_compatible()
    {
      return RawData_::template is_shape_compatible<Shape_>();
    }

    static void stats()
    {
      _templates().stats();
    }

    /**
     * \brief Returns maximum number of children a template produces
     *
     * \tparam template_dim_ Dimension of template
     * \tparam child_dim_ Dimension of children
     */
    template<int template_dim_, int child_dim_>
    static constexpr int max_children()
    {
      return RawData_::template max_children<template_dim_, child_dim_>();
    }

    /**
     * \brief Adjusts SubdivisionLevels to match TemplateSet requirements
     *
     * Will only add subdivision levels to match requirements. No subdivision levels will be removed.
     * No-op for quadrilateral meshes, as all 2D templates exists and all markings are valid.
     *
     * \param[inout] sdls SubdivisionLevels to adjust
     */
    template<typename Mesh_>
    static Intern::RefinementField<IsolatedPointVertexMarking> make_refinement_field(const Mesh_& mesh, Geometry::SubdivisionLevels& sdls)
    {
      return Intern::make_isolated_point_refinement_field(mesh, _templates(), sdls);
    }

    template<typename Shape_>
    static VertexMarkerType spread_refinement_field(
        const EntityReference& ref,
        const Intern::RefinementFieldTuple<VertexMarkerType, Shape::FaceTraits<Shape_, 0>::count>& levels)
    {
      return Intern::spread_refinement_field<Shape_>(ref, levels);
    }

    template<typename Shape_>
    static const RefinementTemplate<Shape_>& get_template(IsolatedPointRefinementType<Shape_> type)
    {
      return _templates().template get_template<Shape_::dimension>(type);
    }

    template<typename Shape_>
    static bool has_template(StandardRefinementType<Shape_> type)
    {
      return _templates().template has_template<Shape_::dimension>(type);
    }

    template<typename Shape_, int dim_>
    static std::pair<Index, int> correct_for_orientation(IsolatedPointRefinementType<Shape_> type, int orientation, Index idx)
    {
      return _templates().template correct_for_orientation<Shape_::dimension, dim_>(type, orientation, idx);
    }

    /**
     * \brief Returns number of children a specific template produces
     *
     * This is intended as a utility function for iterating over children without creating full TemplateQuerys.
     *
     * \tparam template_dim_ Template dimension
     * \tparam child_dim_ Child dimension
     *
     * \param[in] type Type of template
     *
     * \returns Number of children of dimension \c child_dim_ of template with dimension \c template_dim_ and type \c
     * type
     */
    template<int template_dim_, int child_dim_>
    static Index num_children(IsolatedPointRefinementType<typename Shape::FaceTraits<MaxShape, template_dim_>::ShapeType> type)
    {
      return _templates().template get_template<template_dim_>(type).template num_entities<child_dim_>();
    }

  private:
    static TemplateBuilder<RawData_>& _templates()
    {
      static TemplateBuilder<RawData_> _static_templates = {};
      return _static_templates;
    }

  };
} // namespace FEAT::Geometry
