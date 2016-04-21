#pragma once
#ifndef KERNEL_GEOMETRY_BOUNDARY_FACTORY_HPP
#define KERNEL_GEOMETRY_BOUNDARY_FACTORY_HPP 1

// includes, FEAST
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/intern/boundary_computer.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /**
     * \brief BoundaryFactory class template
     */
    template<typename Mesh_>
    class BoundaryFactory DOXY({});

    /**
     * \brief BoundaryFactory implementation for ConformalMesh
     *
     * The boundary factory is a MeshPart factory, which creates a MeshPart without topology containing all boundary
     * faces for a given conformal mesh.
     *
     * \author Peter Zajac
     */
    template<
      typename Shape_,
      int num_coords_,
      int stride_,
      typename Coord_>
    class BoundaryFactory<ConformalMesh<Shape_, num_coords_, stride_, Coord_> > :
      public Factory<MeshPart<ConformalMesh<Shape_, num_coords_, stride_, Coord_>>>
    {
    public:
      /// Our base class
      typedef Factory<MeshPart<ConformalMesh<Shape_, num_coords_, stride_, Coord_>>> BaseClass;
      /// the input mesh type
      typedef ConformalMesh<Shape_, num_coords_, stride_, Coord_> InputMeshType;
      /// The MeshPart type
      typedef MeshPart<InputMeshType> MeshType;
      /// target set holder type
      typedef typename MeshType::TargetSetHolderType TargetSetHolderType;

    private:
      /// a reference to the input mesh
      const InputMeshType& _mesh_in;
      /// a boundary face computer object for the dirty work
      Intern::BoundaryFaceComputer<Shape_> _face_computer;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] mesh_in
       * The mesh for which the boundary cellset is to be computed.
       */
      explicit BoundaryFactory(const InputMeshType& mesh_in) :
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
      virtual void fill_attribute_sets(typename BaseClass::AttributeHolderType& DOXY(target_set_holder)) override
      {
      }
      /// Fills the MeshPart's index_set_holder, except that it doesn't as there is no topology
      virtual void fill_index_sets(typename BaseClass::IndexSetHolderType*& DOXY(index_set_holder)) override
      {
      }

      /// Fills the MeshPart's target set.
      virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
      {
        _face_computer.fill_target_sets(target_set_holder);
      }

    }; // BoundaryFactory<ConformalMesh<...>>
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_BOUNDARY_FACTORY_HPP
