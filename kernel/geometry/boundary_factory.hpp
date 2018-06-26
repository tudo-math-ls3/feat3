#pragma once
#ifndef KERNEL_GEOMETRY_BOUNDARY_FACTORY_HPP
#define KERNEL_GEOMETRY_BOUNDARY_FACTORY_HPP 1

// includes, FEAT
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/intern/boundary_computer.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /**
     * \brief BoundaryFactory implementation
     *
     * The boundary factory is a MeshPart factory, which creates a MeshPart without
     * topology containing all boundary faces for a given conformal mesh.
     *
     * \author Peter Zajac
     */
    template<typename ParentMesh_>
    class BoundaryFactory :
      public Factory<MeshPart<ParentMesh_>>
    {
    public:
      /// Our base class
      typedef Factory<MeshPart<ParentMesh_>> BaseClass;
      /// the parent mesh type
      typedef ParentMesh_ ParentMeshType;
      /// The MeshPart type
      typedef MeshPart<ParentMeshType> MeshType;
      /// the shape type
      typedef typename MeshType::ShapeType ShapeType;
      /// target set holder type
      typedef typename MeshType::TargetSetHolderType TargetSetHolderType;

    private:
      /// a reference to the input mesh
      const ParentMeshType& _mesh_in;
      /// a boundary face computer object for the dirty work
      Intern::BoundaryFaceComputer<ShapeType> _face_computer;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] mesh_in
       * The mesh for which the boundary meshpart is to be computed.
       */
      explicit BoundaryFactory(const ParentMeshType& mesh_in) :
        _mesh_in(mesh_in),
        _face_computer(mesh_in.get_index_set_holder())
      {
      }

      /// Returns the number of entities.
      virtual Index get_num_entities(int dim) override
      {
        return _face_computer.get_num_entities(dim);
      }

      /// Fills the MeshPart's target set, except that it doesn't
      virtual void fill_attribute_sets(typename BaseClass::AttributeSetContainer&) override
      {
      }

      /// Fills the MeshPart's index_set_holder, except that it doesn't as there is no topology
      virtual void fill_index_sets(typename BaseClass::IndexSetHolderType*&) override
      {
      }

      /// Fills the MeshPart's target set.
      virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
      {
        _face_computer.fill_target_sets(target_set_holder);
      }

    }; // class BoundaryFactory<...>
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_BOUNDARY_FACTORY_HPP
