#pragma once
#ifndef KERNEL_GEOMETRY_STANDARD_REFINERY_HPP
#define KERNEL_GEOMETRY_STANDARD_REFINERY_HPP 1

// includes, FEAST
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/conformal/standard_refinement/index_refine_wrappers.hpp>
#include <kernel/geometry/conformal/standard_refinement/target_refine_wrappers.hpp>
#include <kernel/geometry/conformal/standard_refinement/vertex_refine_wrappers.hpp>
#include <kernel/geometry/structured_mesh.hpp>
#include <kernel/geometry/structured/vertex_refiner.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /**
     * \brief Standard Refinery class template
     *
     * \author Peter Zajac
     */
    template<typename Mesh_>
    class StandardRefinery;

    /* ************************************************************************************************************* */

    /**
     * \brief Standard Refinery implementation for ConformalMesh
     *
     * \author Peter Zajac
     */
    template<typename MeshPolicy_>
    class StandardRefinery< ConformalMesh<MeshPolicy_> >
    {
    public:
      typedef ConformalMesh<MeshPolicy_> MeshType;
      typedef typename MeshType::ShapeType ShapeType;

    protected:
      const MeshType& _coarse_mesh;
      MeshType* _fine_mesh;

    public:
      explicit StandardRefinery(const MeshType& coarse_mesh) :
        _coarse_mesh(coarse_mesh),
        _fine_mesh(nullptr)
      {
        CONTEXT(name() + "::StandardRefinery()");
      }

      virtual ~StandardRefinery()
      {
        CONTEXT(name() + "::~StandardRefinery()");
      }

      MeshType* refine()
      {
        CONTEXT(name() + "::refine()");
        using namespace Conformal::StandardRefinement;

        // get number of coordinates and vertex stride
        int num_coords = _coarse_mesh.get_vertex_set().get_num_coords();
        int vertex_stride = _coarse_mesh.get_vertex_set().get_stride();

        // get number of entities in coarse mesh
        Index num_entities_coarse[MeshType::shape_dim+1];
        Index num_entities_fine[MeshType::shape_dim+1];
        for(int i(0); i <= MeshType::shape_dim; ++i)
        {
          num_entities_fine[i] = num_entities_coarse[i] = _coarse_mesh.get_num_entities(i);
        }

        // calculate number of entities in fine mesh
        EntityCountWrapper<ShapeType>::query(num_entities_fine);

        // allocate a fine mesh
        _fine_mesh = new MeshType(num_entities_fine, num_coords, vertex_stride);

        // refine vertices
        VertexRefineWrapper<ShapeType, typename MeshType::VertexSetType>::refine(
          _fine_mesh->get_vertex_set(),
          _coarse_mesh.get_vertex_set(),
          _coarse_mesh.get_index_set_holder());

        // refine indices
        IndexRefineWrapper<ShapeType>::refine(
          _fine_mesh->get_index_set_holder(),
          num_entities_coarse,
          _coarse_mesh.get_index_set_holder());

        // return fine mesh
        return _fine_mesh;
      }

      static String name()
      {
        return "StandardRefinery<ConformalMesh<...>>";
      }
    }; // class StandardRefinery<ConformalMesh<...>>

    /* ************************************************************************************************************* */

    /**
     * \brief Standard Refinery implementation for ConformalSubMesh
     *
     * \author Peter Zajac
     */
    template<typename MeshPolicy_>
    class StandardRefinery< ConformalSubMesh<MeshPolicy_> >
    {
    public:
      typedef ConformalSubMesh<MeshPolicy_> MeshType;
      typedef typename MeshType::ShapeType ShapeType;

    protected:
      const MeshType& _coarse_mesh;
      MeshType* _fine_mesh;

    public:
      explicit StandardRefinery(const MeshType& coarse_mesh) :
        _coarse_mesh(coarse_mesh),
        _fine_mesh(nullptr)
      {
        CONTEXT(name() + "::StandardRefinery()");
      }

      virtual ~StandardRefinery()
      {
        CONTEXT(name() + "::~StandardRefinery()");
      }

      template<typename ParentMesh_>
      MeshType* refine(const ParentMesh_& parent_mesh)
      {
        CONTEXT(name() + "::refine()");
        using namespace Conformal::StandardRefinement;

        // get number of coordinates and vertex stride
        int num_coords = _coarse_mesh.get_vertex_set().get_num_coords();
        int vertex_stride = _coarse_mesh.get_vertex_set().get_stride();

        // get number of entities in coarse mesh
        Index num_entities_coarse[MeshType::shape_dim+1];
        Index num_entities_fine[MeshType::shape_dim+1];
        Index num_entities_parent[MeshType::shape_dim+1];
        for(int i(0); i <= MeshType::shape_dim; ++i)
        {
          num_entities_fine[i] = num_entities_coarse[i] = _coarse_mesh.get_num_entities(i);
          num_entities_parent[i] = parent_mesh.get_num_entities(i);
        }

        // calculate number of entities in fine mesh
        EntityCountWrapper<ShapeType>::query(num_entities_fine);

        // allocate a fine mesh
        _fine_mesh = new MeshType(num_entities_fine, num_coords, vertex_stride);

        // refine vertices
        VertexRefineWrapper<ShapeType, typename MeshType::VertexSetType>::refine(
          _fine_mesh->get_vertex_set(),
          _coarse_mesh.get_vertex_set(),
          _coarse_mesh.get_index_set_holder());

        // refine indices
        IndexRefineWrapper<ShapeType>::refine(
          _fine_mesh->get_index_set_holder(),
          num_entities_coarse,
          _coarse_mesh.get_index_set_holder());

        // refine target indices
        TargetRefineWrapper<ShapeType>::refine(
          _fine_mesh->get_target_set_holder(),
          num_entities_parent,
          _coarse_mesh.get_target_set_holder(),
          _coarse_mesh.get_index_set_holder(),
          parent_mesh.get_index_set_holder());

        // return fine mesh
        return _fine_mesh;
      }

      static String name()
      {
        return "StandardRefinery<ConformalSubMesh<...>>";
      }
    }; // class StandardRefinery<ConformalSubMesh<...>>

    /* ************************************************************************************************************* */

    /**
     * \brief Standard Refinery implementation for StructuredMesh
     *
     * \author Peter Zajac
     */
    template<typename MeshPolicy_>
    class StandardRefinery< StructuredMesh<MeshPolicy_> >
    {
    public:
      typedef StructuredMesh<MeshPolicy_> MeshType;
      typedef typename MeshType::ShapeType ShapeType;

    protected:
      const MeshType& _coarse_mesh;
      MeshType* _fine_mesh;

    public:
      explicit StandardRefinery(const MeshType& coarse_mesh) :
        _coarse_mesh(coarse_mesh),
        _fine_mesh(nullptr)
      {
        CONTEXT(name() + "::StandardRefinery()");
      }

      virtual ~StandardRefinery()
      {
        CONTEXT(name() + "::~StandardRefinery()");
      }

      MeshType* refine()
      {
        CONTEXT(name() + "::refine()");

        using namespace Structured;

        // get number of slices in coarse mesh
        Index num_slices_coarse[MeshType::shape_dim];
        Index num_slices_fine[MeshType::shape_dim];
        for(int i(0); i < MeshType::shape_dim; ++i)
        {
          num_slices_coarse[i] = _coarse_mesh.get_num_slices(i);
          num_slices_fine[i] = 2*num_slices_coarse[i];
        }

        // allocate a fine mesh
        _fine_mesh = new MeshType(num_slices_fine);

        // refine vertices
        VertexRefiner<ShapeType, typename MeshType::VertexSetType>
          ::refine(_fine_mesh->get_vertex_set(), _coarse_mesh.get_vertex_set(), num_slices_coarse);

        // return fine mesh
        return _fine_mesh;
      }

      static String name()
      {
        return "StandardRefinery<StructuredMesh<...>>";
      }
    }; // class StandardRefinery<StructuredMesh<...>>

  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_STANDARD_REFINERY_HPP
