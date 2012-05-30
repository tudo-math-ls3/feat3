#pragma once
#ifndef KERNEL_GEOMETRY_STANDARD_REFINERY_HPP
#define KERNEL_GEOMETRY_STANDARD_REFINERY_HPP 1

// includes, FEAST
#include <kernel/geometry/cell_sub_set.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/structured_mesh.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /**
     * \brief Standard Refinery class template
     *
     * \todo think about eliminating the specialisations...
     *
     * \author Peter Zajac
     */
    template<typename Mesh_>
#ifndef DOXYGEN
    class StandardRefinery;
#else
    class StandardRefinery
    {
    public:
      /// mesh type
      typedef Mesh_ MeshType;

      /**
       * \brief Constructor.
       *
       * \param[in] coarse_mesh
       * The coarse mesh that is to be refined.
       */
      explicit StandardRefinery(const MeshType& coarse_mesh);

      /**
       * \brief Refines the coarse mesh.
       *
       * This function applies the standard refinement algorithm to refine the coarse mesh that
       * has been passed to the constructor of this class.
       *
       * \note
       * The refinery does not delete the refined mesh returned by this function upon destruction,
       * therefore it is the responsibility of the caller to delete the refined mesh once it has
       * served its purpose.
       *
       * \returns
       * A pointer to the refined mesh.
       */
      MeshType* refine();
    };
#endif

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
        return (_fine_mesh = _coarse_mesh.refine());
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
        return (_fine_mesh = _coarse_mesh.refine(parent_mesh));
      }

      static String name()
      {
        return "StandardRefinery<ConformalSubMesh<...>>";
      }
    }; // class StandardRefinery<ConformalSubMesh<...>>

    /* ************************************************************************************************************* */

    /**
     * \brief Standard Refinery implementation for CellSubSet
     *
     * \author Peter Zajac
     */
    template<typename Shape_>
    class StandardRefinery< CellSubSet<Shape_> >
    {
    public:
      typedef CellSubSet<Shape_> MeshType;
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
        return (_fine_mesh = _coarse_mesh.refine(parent_mesh));
      }

      static String name()
      {
        return "StandardRefinery<CellSubSet<...>>";
      }
    }; // class StandardRefinery<CellSubSet<...>>

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
        return (_fine_mesh = _coarse_mesh.refine());
      }

      static String name()
      {
        return "StandardRefinery<StructuredMesh<...>>";
      }
    }; // class StandardRefinery<StructuredMesh<...>>

  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_STANDARD_REFINERY_HPP
