#pragma once
#ifndef KERNEL_GEOMETRY_BOUNDARY_FACTORY_HPP
#define KERNEL_GEOMETRY_BOUNDARY_FACTORY_HPP 1

// includes, FEAST
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/cell_sub_set.hpp>
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
     * The boundary factory is a CellSubSet-factory, which creates a cell-sub-set containing all boundary
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
      public Factory<CellSubSet<Shape_> >
    {
    public:
      /// the input mesh type
      typedef ConformalMesh<Shape_, num_coords_, stride_, Coord_> InputMeshType;
      /// the cell set type
      typedef CellSubSet<Shape_> MeshType;
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
      virtual Index get_num_entities(int dim)
      {
        return _face_computer.get_num_entities(dim);
      }

      /// Fills the cellset's target set.
      virtual void fill_target_sets(TargetSetHolderType& target_set_holder)
      {
        _face_computer.fill_target_sets(target_set_holder);
      }
    }; // BoundaryFactory<ConformalMesh<...>>
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_BOUNDARY_FACTORY_HPP
