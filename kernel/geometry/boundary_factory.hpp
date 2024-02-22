// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

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
       * A \resident reference to mesh for which the boundary meshpart is to be computed.
       */
      explicit BoundaryFactory(const ParentMeshType& mesh_in) :
        _mesh_in(mesh_in)
      {
        _face_computer.compute_all(mesh_in.get_index_set_holder());
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
      virtual void fill_index_sets(std::unique_ptr<typename BaseClass::IndexSetHolderType>&) override
      {
      }

      /// Fills the MeshPart's target set.
      virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
      {
        _face_computer.fill_target_sets(target_set_holder);
      }

    }; // class BoundaryFactory<...>

    /**
     * \brief MaskedBoundaryFactory implementation
     *
     * This class is an extension of the BoundaryFactory class, which creates a mesh part containing
     * the set of all unmasked boundary faces for a given mesh. Facets of the mesh can be masked by
     * added them individually using the #add_mask_facet() function or by adding an entire mesh
     * part by using the #add_mask_mehspart() function. Once all masked facets/meshparts have been
     * added, one has to call the #compile() function before creating the mesh part using this factory.
     *
     * \author Peter Zajac
     */
    template<typename ParentMesh_>
    class MaskedBoundaryFactory :
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
      /// the facet dimension
      static constexpr int facet_dim = ShapeType::dimension - 1;

    private:
      /// a reference to the input mesh
      const ParentMeshType& _mesh_in;
      /// a boundary face computer object for the dirty work
      Intern::BoundaryFaceComputer<ShapeType> _face_computer;
      /// a mask vector for the facets
      std::vector<int> _facet_mask;

    public:
      /**
      * \brief Constructor
      *
      * \param[in] mesh_in
      * A \resident reference to mesh for which the boundary meshpart is to be computed.
      */
      explicit MaskedBoundaryFactory(const ParentMeshType& mesh_in) :
        _mesh_in(mesh_in),
        _facet_mask(mesh_in.get_num_entities(facet_dim), 0)
      {
      }

      // Adds a single facet to the mask.
      void add_mask_facet(Index facet_idx)
      {
        _facet_mask.at(facet_idx) = 1;
      }

      /// Adds all facets in a mesh-part to the mask.
      void add_mask_meshpart(const MeshPart<ParentMeshType>& part)
      {
        const TargetSet& trg_set = part.template get_target_set<facet_dim>();
        for(Index i(0); i < trg_set.get_num_entities(); ++i)
          _facet_mask.at(trg_set[i]) = 1;
      }

      /// Compiles the factory
      void compile()
      {
        _face_computer.compute_masked(_mesh_in.get_index_set_holder(), _facet_mask);
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
      virtual void fill_index_sets(std::unique_ptr<typename BaseClass::IndexSetHolderType>&) override
      {
      }

      /// Fills the MeshPart's target set.
      virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
      {
        _face_computer.fill_target_sets(target_set_holder);
      }

    }; // class MaskedBoundaryFactory<...>
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_BOUNDARY_FACTORY_HPP
