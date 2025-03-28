// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/intern/boundary_computer.hpp>
#include <kernel/util/dist.hpp>

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
      /// the facet dimension
      static constexpr int shape_dim = ShapeType::dimension;

    private:
      /// a reference to the input mesh
      const ParentMeshType& _mesh_in;
      /// a boundary face computer object for the dirty work
      Intern::BoundaryFaceComputer<ShapeType> _face_computer;
      /// the face masks for each face dimensions
      std::array<std::vector<int>, std::size_t(shape_dim)> _bnd_masks;
      /// the (cell_dim_-1)-dimensional boundary face indices for each face dimension
      std::array<std::vector<Index>, std::size_t(shape_dim)> _face_idx;

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
        _face_computer.compute_all(mesh_in.get_index_set_holder(), _bnd_masks, _face_idx);
      }

      /// Returns the number of entities.
      virtual Index get_num_entities(int dim) override
      {
        if(dim < shape_dim)
          return Index(_face_idx.at(std::size_t(dim)).size());
        else
          return Index(0);
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
        _face_computer.fill_target_sets(target_set_holder, _face_idx);
      }

    }; // class BoundaryFactory<...>

    /**
     * \brief Creates a new boundary mesh-part for a given mesh
     *
     * \param[in] mesh
     * A \transient reference to the mesh for which a boundary mesh-part is to be created.
     *
     * \returns
     * A mesh-part containing all entities for which lie on the mesh boundary
     */
    template<typename Mesh_>
    MeshPart<Mesh_> make_boundary_meshpart(const Mesh_& mesh)
    {
      BoundaryFactory<Mesh_> factory(mesh);
      return factory.make();
    }

    /**
     * \brief Creates a new boundary mesh-part for a given mesh
     *
     * \param[in] mesh
     * A \transient reference to the mesh for which a boundary mesh-part is to be created.
     *
     * \returns
     * A mesh-part containing all entities for which lie on the mesh boundary
     */
    template<typename Mesh_>
    std::unique_ptr<MeshPart<Mesh_>> make_unique_boundary_meshpart(const Mesh_& mesh)
    {
      BoundaryFactory<Mesh_> factory(mesh);
      return factory.make_unique();
    }


    /**
     * \brief MaskedBoundaryFactory implementation
     *
     * This class is an extension of the BoundaryFactory class, which creates a mesh part containing
     * the set of all unmasked boundary faces for a given mesh. Facets of the mesh can be masked by
     * adding them individually using the #add_mask_facet() function or by adding an entire mesh
     * part by using the #add_mask_meshpart() function. Once all masked facets/meshparts have been
     * added, one has to call the #compile() function before creating the mesh part using this factory.
     *
     * \attention
     * This class <b>must not</b> be used for partitioned meshes, because it may miss certain parts
     * of the boundary! In the case of partitioned meshes, use the GlobalMaskedBoundaryFactory
     * instead! Please consider the documentation of the GlobalMaskedBoundaryFactory for a simple
     * example where this MaskedBoundaryFactory class would fail to create the correct boundary
     * meshpart in a partitioned mesh case.
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
      /// the shape dimension
      static constexpr int shape_dim = ShapeType::dimension;
      /// the facet dimension
      static constexpr int facet_dim = shape_dim - 1;

    private:
      /// a reference to the input mesh
      const ParentMeshType& _mesh_in;
      /// a boundary face computer object for the dirty work
      Intern::BoundaryFaceComputer<ShapeType> _face_computer;
      /// a mask vector for the facets
      std::vector<int> _facet_mask;
      /// the face masks for each face dimensions
      std::array<std::vector<int>, std::size_t(shape_dim)> _bnd_masks;
      /// the (cell_dim_-1)-dimensional boundary face indices for each face dimension
      std::array<std::vector<Index>, std::size_t(shape_dim)> _face_idx;

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

      // Adds a single facet to the boundary facet mask.
      void add_mask_facet(Index facet_idx)
      {
        _facet_mask.at(facet_idx) = 1;
      }

      /**
       * \brief Adds all facets in a mesh-part to the boundary facet mask
       *
       * \param[in] part
       * A \transient reference to the mesh-part to be added
       */
      void add_mask_meshpart(const MeshPart<ParentMeshType>& part)
      {
        const TargetSet& trg_set = part.template get_target_set<facet_dim>();
        for(Index i(0); i < trg_set.get_num_entities(); ++i)
          _facet_mask.at(trg_set[i]) = 1;
      }

      /// Compiles the factory
      void compile()
      {
        _face_computer.compute_masks(_mesh_in.get_index_set_holder(), _bnd_masks, _face_idx, _facet_mask);
        _face_computer.compute_faces(_bnd_masks, _face_idx);
      }

      /// Returns the number of entities.
      virtual Index get_num_entities(int dim) override
      {
        if(dim < shape_dim)
          return Index(_face_idx.at(std::size_t(dim)).size());
        else
          return Index(0);
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
        _face_computer.fill_target_sets(target_set_holder, _face_idx);
      }
    }; // class MaskedBoundaryFactory<...>

    /**
     * \brief GlobalMaskedBoundaryFactory implementation
     *
     * This class is an extension of the MaskedBoundaryFactory class, which creates a mesh part
     * containing the set of all unmasked boundary faces for a given partitioned mesh.
     * Facets of the mesh can be masked by added them individually using the #add_mask_facet()
     * function or by adding an entire mesh part by using the #add_mask_meshpart() function.
     * Furthermore, one \b must add all halos and their corresponding ranks by calling the
     * #add_halo() function. Once all masked facets/meshparts and halos have been added, one has
     * to call the #compile() function before creating the mesh part using this factory.
     *
     * To illustrate why one cannot use the MaskedBoundaryFactory in the case of a partitioned mesh
     * consider the following L-domain example with 3 MPI ranks and 1 element per rank:
       \verbatim
       +---+
       | 2 |
       +---X---+
       | 0 | 1 |
       +---+---+
       \endverbatim
     * On ranks 1 and 2, the vertex 'X' is correctly considered to be a boundary vertex, because
     * it is adjacent to a boundary edge, however, on rank 0 this vertex is only adjacent to halo
     * edges, but not to any boundary edge, thus it would be incorrectly treated as an inner vertex.
     * This class solves this problem by explicitly synchronizing the detected boundary entities
     * via the supplied halos.
     *
     * \author Peter Zajac
     */
    template<typename ParentMesh_>
    class GlobalMaskedBoundaryFactory :
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
      /// the shape dimension
      static constexpr int shape_dim = ShapeType::dimension;
      /// the facet dimension
      static constexpr int facet_dim = shape_dim - 1;

    private:
      /// a reference to the input mesh
      const ParentMeshType& _mesh_in;
      /// a boundary face computer object for the dirty work
      Intern::BoundaryFaceComputer<ShapeType> _face_computer;
      /// a mask vector for the facets
      std::vector<int> _facet_mask;
      /// the face masks for each face dimensions
      std::array<std::vector<int>, std::size_t(shape_dim)> _bnd_masks;
      /// the (cell_dim_-1)-dimensional boundary face indices for each face dimension
      std::array<std::vector<Index>, std::size_t(shape_dim)> _face_idx;
      /// a map of all halos
      std::vector<std::pair<int, const MeshType*>> _halos;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] mesh_in
       * A \resident reference to mesh for which the boundary meshpart is to be computed.
       */
      explicit GlobalMaskedBoundaryFactory(const ParentMeshType& mesh_in) :
        _mesh_in(mesh_in),
        _facet_mask(mesh_in.get_num_entities(facet_dim), 0)
      {
      }

      // Adds a single facet to the mask.
      void add_mask_facet(Index facet_idx)
      {
        _facet_mask.at(facet_idx) = 1;
      }

      /**
       * \brief Adds all facets in a mesh-part to the boundary facet mask
       *
       * \param[in] part
       * A \transient reference to the mesh-part to be added
       */
      void add_mask_meshpart(const MeshPart<ParentMeshType>& part)
      {
        const TargetSet& trg_set = part.template get_target_set<facet_dim>();
        for(Index i(0); i < trg_set.get_num_entities(); ++i)
          _facet_mask.at(trg_set[i]) = 1;
      }

      /**
       * \brief Adds a halo to the mask and saves it for synchronization
       *
       * \param[in] rank
       * The rank of the halo to be added
       *
       * \param[in] halo
       * A \resident reference to the halo to be added.
       */
      void add_halo(int rank, const MeshPart<ParentMeshType>& halo)
      {
        // add the halo as a mask meshpart
        add_mask_meshpart(halo);
        // and save the halo
        //_halos.emplace(rank, &halo);
        _halos.push_back(std::make_pair(rank, &halo));
      }

      /**
       * \brief Compiles the factory
       *
       * This function performs the necessary synchronization over the communicator to ensure that
       * the boundary is correctly computes for all patches.
       *
       * \param[in] comm
       * A \transient reference to the communicator where the partitioning of the underlying mesh is defined on
       */
      void compile(const Dist::Comm& comm)
      {
        _face_computer.compute_masks(_mesh_in.get_index_set_holder(), _bnd_masks, _face_idx, _facet_mask);
        _sync_bnd_masks(comm);
        _face_computer.compute_faces(_bnd_masks, _face_idx);
      }

      /// Returns the number of entities.
      virtual Index get_num_entities(int dim) override
      {
        if(dim < shape_dim)
          return Index(_face_idx.at(std::size_t(dim)).size());
        else
          return Index(0);
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
        _face_computer.fill_target_sets(target_set_holder, _face_idx);
      }

    protected:
      /**
       * \brief Synchronizes the boundary masks over the communicators
       */
      void _sync_bnd_masks(const Dist::Comm& comm)
      {
        if(_halos.empty())
          return;

        // allocate halo buffers
        std::vector<std::vector<int>> halo_send_bufs(_halos.size()), halo_recv_bufs(_halos.size());

        // create send buffers and allocate receive buffers
        for(std::size_t i(0); i < _halos.size(); ++i)
        {
          // count total number of entries in halo
          Index count = 0u;
          for(int dim = 0; dim < shape_dim; ++dim) // we don't need to sync dim=shape_dim here
            count += _halos[i].second->get_num_entities(dim);

          halo_recv_bufs[i].resize(count);
          halo_send_bufs[i].resize(count); // build_halo_buffer uses push_back
          _face_computer.build_halo_buffer(halo_send_bufs[i], _bnd_masks, _halos[i].second->get_target_set_holder());
        }

        // post receive requests
        Dist::RequestVector recv_reqs(_halos.size());
        for(std::size_t i(0); i < _halos.size(); ++i)
          recv_reqs.push_back(comm.irecv(halo_recv_bufs[i].data(), halo_recv_bufs[i].size(), _halos[i].first));

        // post send requests
        Dist::RequestVector send_reqs(_halos.size());
        for(std::size_t i(0); i < _halos.size(); ++i)
          send_reqs.push_back(comm.isend(halo_send_bufs[i].data(), halo_send_bufs[i].size(), _halos[i].first));

        // wait for all receives to finish
        recv_reqs.wait_all();

        // process receives halo buffers
        for(std::size_t i(0); i < _halos.size(); ++i)
        {
          _face_computer.mask_halo_buffer(_bnd_masks, halo_recv_bufs[i], _halos[i].second->get_target_set_holder());
        }

        // wait for all sends to finish
        send_reqs.wait_all();
      }
    }; // class GlobalMaskedBoundaryFactory<...>
  } // namespace Geometry
} // namespace FEAT
