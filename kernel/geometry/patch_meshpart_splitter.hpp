// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/geometry/mesh_part.hpp>

// includes, system
#include <map>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    // Forward declarations
    class PatchPartMap;

    template<typename, int>
    class PatchPartMapHolder;
    /// \endcond

    /**
     * \brief Class for creating MeshParts referring to a patch BaseMesh
     *
     * \tparam Mesh_
     * The mesh type the MeshParts etc. refer to.
     *
     * Generic class template. The implementations are in specializations wrt. Mesh_.
     */
    template<typename Mesh_>
    class PatchMeshPartSplitter;

    /**
     * \brief Class for creating MeshParts referring to a patch BaseMesh, ConformalMesh specialization.
     *
     * \tparam Shape_
     * Shape type of the mesh cells
     *
     * \tparam num_coords_
     * World dimension
     *
     * \tparam Coord_
     * Floating point type for vertex coordinates
     *
     * Like the PatchPartMap and PatchPartMapHolder, this assembles information about the BaseMesh and PatchMeshPart
     * first and can later be used to create several MeshParts referring to the PatchMeshPart by calling build() and
     * then the Factory interface, one at a time.
     *
     * \author Jordi Paul
     */
    template<typename Shape_, int num_coords_, typename Coord_>
    class PatchMeshPartSplitter<Geometry::ConformalMesh<Shape_, num_coords_, Coord_>> :
      public Factory<MeshPart<Geometry::ConformalMesh<Shape_, num_coords_, Coord_>>>
    {
    public:
      /// Our shape type
      typedef Shape_ ShapeType;
      /// The shape dimension
      static constexpr int shape_dim = ShapeType::dimension;
      /// The mesh type
      typedef Geometry::ConformalMesh<Shape_, num_coords_, Coord_> MeshType;
      /// The MeshPart type
      typedef Geometry::MeshPart<MeshType> MeshPartType;
      /// Data type for attributes
      typedef typename MeshType::VertexSetType::CoordType AttributeDataType;
      /// Mesh attribute holder type
      typedef typename MeshPartType::AttributeSetContainer AttributeSetContainer;
      /// Index set holder type
      typedef typename MeshPartType::IndexSetHolderType IndexSetHolderType;
      /// Target set holder type
      typedef typename MeshPartType::TargetSetHolderType TargetSetHolderType;

    protected:
      /// The BaseMesh we refer to
      const MeshType& _base_mesh;
      /// The MeshPart identifying this patch
      const MeshPartType& _patch_mesh_part;

      /// The mesh attributes of the current BaseMesh MeshPart
      const AttributeSetContainer* _cur_part_attribute_container;
      /// The topology of the current BaseMesh MeshPart
      const IndexSetHolder<ShapeType>* _cur_part_topology;
      /// The BaseMesh MeshPart to PatchMeshPart MeshPart mapping information
      PatchPartMapHolder<ShapeType, ShapeType::dimension> _part_holder;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] base_mesh
       * A \resident reference to the (global) BaseMesh
       *
       * \param[in] patch_mesh_part
       * A \resident reference to the MeshPart identifying the Patch
       */
      explicit PatchMeshPartSplitter(const MeshType& base_mesh, const MeshPartType& patch_mesh_part) :
        _base_mesh(base_mesh),
        _patch_mesh_part(patch_mesh_part),
        _cur_part_attribute_container(nullptr),
        _cur_part_topology(nullptr),
        _part_holder(patch_mesh_part.get_target_set_holder())
        {
        }

      /**
       * \brief Virtual destructor
       */
      virtual ~PatchMeshPartSplitter()
      {
      }

      /**
       * \brief Assembles information from a MeshPart referring to the BaseMesh
       *
       * \param[in] mesh_part
       * A \transient reference to the basemesh patch meshpart
       *
       * \returns
       * true if there is a nonempty cut between mesh_part and the PatchMeshPart.
       *
       */
      bool build(const MeshPartType& mesh_part)
      {
        _cur_part_topology = mesh_part.get_topology();
        _cur_part_attribute_container = &(mesh_part.get_mesh_attributes());
        return _part_holder.build(mesh_part.get_target_set_holder(), mesh_part.get_topology());
      }

      /* *************************************************************************************** */
      /* F A C T O R Y   I N T E R F A C E   I M P L E M E N T A T I O N                         */
      /* *************************************************************************************** */

      /**
       * \brief Returns the number of entities.
       *
       * \param[in] dim
       * The dimension of the entity whose count is to be returned. Must be 0 <= \p dim <= shape_dim.
       *
       * \returns
       * The number of entities of dimension \p dim.
       */
      virtual Index get_num_entities(int dim) override
      {
        return _part_holder.get_num_entities(dim);
      }

      /**
       * \brief Fills the attribute sets.
       *
       * \param[in,out] attribute_container
       * A \transient reference to the attribute set holder whose attribute sets are to be filled.
       */
      virtual void fill_attribute_sets(AttributeSetContainer& attribute_container) override
      {
        _part_holder.fill_attribute_sets(attribute_container, *_cur_part_attribute_container);
      }

      /**
       * \brief Fills the index sets.
       *
       * \param[in,out] index_set_holder
       * A \transient reference to the index set holder whose index sets are to be filled.
       */
      virtual void fill_index_sets(std::unique_ptr<IndexSetHolderType>& index_set_holder) override
      {
        XASSERT(index_set_holder == nullptr);

        // If the base MeshPart has a topology, create one for the patch MeshPart, too
        if(_part_holder.has_topology())
        {
          // Get num_entities for all shape dimensions
          Index num_entities[shape_dim+1];
          for(int i(0); i < shape_dim+1; ++i)
            num_entities[i] = get_num_entities(i);

          // Create this MeshPart's topology
          index_set_holder.reset(new IndexSetHolderType(num_entities));
          // This creates vertex-at-shape information for all shapes
          _part_holder.fill_index_sets(*index_set_holder, *_cur_part_topology);
          // Build the redundant index sets
          RedundantIndexSetBuilder<ShapeType>::compute(*index_set_holder);
        }
      }

      /**
       * \brief Fills the target sets.
       *
       * \param[in,out] target_set_holder
       * A \transient reference to the target set holder whose target sets are to be filled.
       */
      virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
      {
        _part_holder.fill_target_sets(target_set_holder);
      }
    }; // class PatchMeshPartSplitter<ConformalMesh<...>>

    /**
     * \brief Mapping class for BaseMesh/PatchBaseMesh and MeshPart index conversions
     *
     * When it is being built, the base mesh for a patch is first represented as a MeshPart (let's call it
     * PatchMeshPart). For every shape dimension, the PatchMeshPart has a target set mapping the patch entities to
     * the BaseMesh entities , and this is used to construct the PatchPartMap. Then the PatchPartMap can be used
     * (via the build()-function) to create the mapping of a MeshPart referring to the BaseMesh to a MeshPart
     * referring to the PatchMeshPart.
     *
     * This information is later used to construct the PatchBaseMesh and the restriction of the BaseMesh MeshParts
     * to it.
     *
     * The PatchPartMap is always created from the PatchMeshPart TargetSet. Then, for any number of TargetSets of
     * MeshParts referring to the BaseMesh, the build() function can be called to assemble the required
     * information so that the Factory interface can then be used to create the real PatchMeshPart MeshPart. One
     * at a time, that is.
     *
     * \author Jordi Paul
     */
    class PatchPartMap
    {
    private:
      /// Index map from BaseMesh entity to PatchMeshPart entity
      std::map<Index,Index> _idx_map;
      /// Index map from base mesh MeshPart entities to patch boundary MeshPart entities
      std::map<Index,Index> _io_map;
      /// This is basically the TargetSet mapping PatchMeshPart MeshPart entities to the PatchMeshPart
      std::vector<Index> _indices;

    public:
      /**
       * \brief Standard constructor
       *
       * \param[in] target_set
       * A \transient reference to the TargetSet mapping Patch MeshPart entities to BaseMesh entities
       */
      explicit PatchPartMap(const TargetSet& target_set) :
        _idx_map(),
        _io_map(),
        _indices()
      {
        // loop over all indices
        for(Index i(0); i < target_set.get_num_entities(); ++i)
        {
          // This maps the BaseMesh entity j = target_set[i] to the patch MeshPart entity i
          XASSERTM(_idx_map.emplace(target_set[i], i).second, "internal error");
        }
      }

      /**
       * \brief Creates mapping information from a BaseMesh MeshPart to a Patch MeshPart boundary MeshPart
       *
       * \param[in] target_in
       * A \transient reference to the Target set of the BaseMesh MeshPart
       *
       * \returns
       * true if there is a nonempty cut between the entities referenced by target_in and the PatchMeshPart.
       */
      bool build(const TargetSet& target_in)
      {
        // Very important: Clear old information
        _indices.clear();
        _io_map.clear();

        // We need information for every entity in the BaseMesh MeshPart
        for(Index i(0); i < target_in.get_num_entities(); ++i)
        {
          // Check if the BaseMesh entity referenced by the BaseMesh MeshPart is in the PatchMeshPart
          auto it = _idx_map.find(target_in[i]);
          if(it != _idx_map.end())
          {
            // We create new entries in the PatchMeshPart MeshPart here, so the entity i in the BaseMesh MeshPart
            // refers to the newly created entity in the PatchMeshPart MeshPart
            _io_map.emplace(i, Index(_indices.size()));
            // Add the new entity
            _indices.push_back(it->second);
          }
        }
        return !_indices.empty();
      }

      /**
       * \brief The size of the PatchMeshPart MeshPart
       *
       * \returns The number of entities in the PatchMeshPart MeshPart
       */
      Index size() const
      {
        return Index(_indices.size());
      }

      /**
       * \brief Fills an attribute set
       *
       * \tparam Attribute_
       * Type of the attribute contained in the passed sets.
       *
       * \param[out] attribute_container_out
       * A \transient reference to the attribute set to fill
       *
       * \param[in] attribute_container_in
       * A \transient reference to the AttributeSet of the BaseMesh MeshPart to copy to the PatchMeshPart's MeshPart
       */
      template<typename Attribute_>
      void fill_attribute_set(
        std::map<String, std::unique_ptr<Attribute_>>& attribute_container_out,
        const std::map<String, std::unique_ptr<Attribute_>>& attribute_container_in) const
      {
        for(auto it(attribute_container_in.begin()); it != attribute_container_in.end(); ++it)
        {
          // Create new attribute for the PatchMeshPart MeshPart
          std::unique_ptr<Attribute_> new_attribute(new Attribute_(
            Index(_indices.size()), it->second->get_dimension()));

          const int dim = it->second->get_dimension();

          // Emplace all elements referred to by the _io_map into new_attribute
          for(auto jt: _io_map)
          {
            for(int k(0); k < dim; ++k)
              new_attribute->operator()(jt.second, k) = it->second->operator()(jt.first, k);
          }

          // Push the new attribute to set_out
          XASSERTM(attribute_container_out.insert(std::make_pair(it->first, std::move(new_attribute))).second, "Error inserting new AttributeSet");
        }
      }

      /**
       * \brief Fills the PatchMeshPart MeshPart's TargetSet
       *
       * \param[in] target_out
       * A \transient reference to the TargetSet mapping the PatchMeshPart MeshPart's entities to the PatchMeshPart
       */
      void fill_target_set(TargetSet& target_out) const
      {
        XASSERT(target_out.get_num_entities() == size());
        for(Index i(0); i < target_out.get_num_entities(); ++i)
        {
          target_out[i] = _indices[i];
        }
      }

      /**
       * \brief Fills the vertex-at-shape IndexSet of the PatchMeshPart MeshPart
       *
       * The whole purpose of MeshPart topologies is to be different from the BaseMesh's topology. Think of the
       * polygon representing the boundary of a 2d domain. Then the last vertex in the last edge is the same as
       * the first vertex in the first edge. Add a parametrization variable. This needs to be different at the
       * vertex from above depending from which edge it is referenced. This can be easily be achieved by adding
       * another vertex to the MeshPart and changing the topology accordingly, meaning the last edge's last vertex
       * gets mapped to the new vertex.
       *
       * So we really need to copy the BaseMesh MeshPart's topology to the PatchMeshPart MeshPart's topology.
       *
       * \tparam IndexSetType_
       * Type of the vertex-at-shape IndexSet to be filled.
       *
       * \param[out] index_set_out
       * A \transient reference to the vertex-at-shape IndexSet of the MeshPart referring to the PatchMeshPart.
       * Already has to have the correct size.
       *
       * \param[in] index_set_in
       * A \transient reference to the vertex-at-shape IndexSet of the MeshPart referring to the BaseMesh.
       *
       * \param[in] vertex_map
       * A \transient reference to the mapping of BaseMesh Meshpart vertices to PatchMeshPart MeshPart vertices.
       * This is the _io_map of the PatchPartMap to dimension-0 entries, but this class has no information about it,
       * so it need to be fetched and passed from the outside.
       */
      template<typename IndexSetType_>
      void fill_index_set(IndexSetType_& index_set_out, const IndexSetType_& index_set_in,
        const std::map<Index, Index>& vertex_map) const
      {
        // We need to check all shape entities in the MeshPart referring to the BaseMesh
        for(Index k(0); k < index_set_in.get_num_entities(); ++k)
        {
          // Check if the BaseMesh entity referred to by the BaseMesh MeshPart is also present in the PatchMeshPart
          auto kt = _io_map.find(k);
          if(kt != _io_map.end())
          {
            // For all vertices at that entity, we need to find them in the PatchMeshPart's MeshPart
            for(int j(0); j < index_set_in.get_num_indices(); ++j)
            {
              // Fetch the BaseMesh MeshPart's vertex number
              Index i_in(index_set_in(k,j));
              // That vertex is definitely in the PatchMeshPart's MeshPart, but we need to catch any errors here
              auto jt = vertex_map.find(i_in);
              if(jt != vertex_map.end())
                // This is the correct number
                index_set_out(kt->second, j) = jt->second;
              else
                XABORTM("Vertex "+stringify(i_in)+" missing in MeshPart topology!");
            }
          }
        }
      }

      /**
       * \brief Gets the mapping BaseMesh MeshPart entities to PatchMeshPart MeshPart entities
       *
       * \returns
       * _io_map
       */
      const std::map<Index, Index>& get_io_map() const
      {
        return _io_map;
      }
    }; // class PatchPartMap

    /**
     * \brief This is a container class for holding PatchPartMaps for all dimensions
     *
     * \tparam Shape_
     * The ShapeType this container refers to. This also defines the dimension etc.
     *
     * \tparam dim_
     * Shape dimension this container refers to
     *
     * \author Jordi Paul
     */
    template<typename Shape_, int dim_ = Shape_::dimension>
    class PatchPartMapHolder :
      public PatchPartMapHolder<Shape_, dim_ - 1>
    {
    public:
      /// Our base class
      typedef PatchPartMapHolder<Shape_, dim_ - 1> BaseClass;

    private:
      /// The PatchPartMap for our dimension
      PatchPartMap _patch_map;
      /// We need to know if the last BaseMesh MeshPart from build() had a topology
      bool _has_topology;

    public:
      /**
       * \brief Standard constructor
       *
       * \param[in] tsh
       * A \transient reference to the TargetSetHolder of the PatchMeshPart
       */
      explicit PatchPartMapHolder(const TargetSetHolder<Shape_>& tsh) :
        BaseClass(tsh),
        _patch_map(tsh.template get_target_set<dim_>()),
        _has_topology(false)
      {
      }

      /**
       * \brief Assembles information from a BaseMesh MeshPart
       *
       * \param[in] tsh
       * A \transient reference to the TargetSetHolder of the MeshPart referring to the BaseMesh.
       *
       * \param[in] ish
       * A \transient pointer to the IndexSetHolder of the MeshPart referring to the BaseMesh.
       * This is nullptr if it does not have a topology.
       *
       * \returns true if this class or its BaseClass added any information
       */
      bool build(const TargetSetHolder<Shape_>& tsh, const IndexSetHolder<Shape_>* ish)
      {
        ish == nullptr ? _has_topology = false : _has_topology = true;

        // Call the BaseClass version (meaning of lower shape dimension) first
        bool b1 = BaseClass::build(tsh, ish);
        bool b2 =  _patch_map.build(tsh.template get_target_set<dim_>());

        return  b1 || b2;
      }

      /**
       * \brief Gets the number of entities
       *
       * \param[in] dim
       * The dimension. Note that \f$ 0 \leq \mathrm{dim} \leq \mathrm{shape\_dim}.\f$
       *
       * \returns The number of entities of a certain dimension
       */
      Index get_num_entities(int dim) const
      {
        XASSERT(dim <= dim_);
        if(dim == dim_)
          return _patch_map.size();
        else
          return BaseClass::get_num_entities(dim);
      }

      /**
       * \brief Does the current BaseMesh MeshPart have a topology?
       *
       * \returns
       * true if the MeshPart does have a topology, false otherwise.
       */
      bool has_topology() const
      {
        return _has_topology;
      }

      /**
       * \brief MeshPart Factory interface: Fills the AttributeSetContainer
       *
       * \tparam AttributeSetContainer_
       * Type of the AttributeSetContainer_.
       *
       * \param[out] attribute_container_out
       * A \transient reference to the AttributeSetContainer of the PatchMeshPart's MeshPart to be filled
       *
       * \param[in] attribute_container_in
       * A \transient reference to the AttributeSetContainer of the BaseMesh's MeshPart to restrict to the PatchMeshPart's MeshPart.
       */
      template<typename AttributeSetContainer_>
      void fill_attribute_sets(
        AttributeSetContainer_& attribute_container_out,
        const AttributeSetContainer_& attribute_container_in) const
      {
        // Fill attribute sets of one shape dimension lower first
        BaseClass::fill_attribute_sets(attribute_container_out, attribute_container_in);
      }

      /**
       * \brief Factory interface: Fills the TargetSetHolder
       *
       * \param[in] tsh
       * A \transient reference to the PatchMeshPart MeshPart's TargetSetHolder
       */
      void fill_target_sets(TargetSetHolder<Shape_>& tsh) const
      {
        BaseClass::fill_target_sets(tsh);
        _patch_map.fill_target_set(tsh.template get_target_set<dim_>());
      }

      /**
       * \brief Factory interface: Fills the IndexSetHolder
       *
       * \tparam IndexSetHolder_
       * The type of the IndexSetHolder
       *
       * \param[in] ish
       * A \transient reference to the PatchMeshPart MeshPart's IndexSetHolder
       *
       * \param[in] ish_in
       * A \transient reference to the PatchMeshPart's IndexSetHolder
       *
       * Note that it has to be checked beforehand if the IndexSetHolder(s) are nullptr.
       */
      template<typename IndexSetHolder_>
      void fill_index_sets(IndexSetHolder_& ish, const IndexSetHolder_& ish_in) const
      {
        BaseClass::fill_index_sets(ish, ish_in);
        _patch_map.fill_index_set(ish.template get_index_set<dim_,0>(), ish_in.template get_index_set<dim_,0>(),
        get_vertex_map());
      }

      /**
       * \brief Gets the vertex map
       *
       * This is needed for the get_index_set function of the PatchPartMap.
       *
       * \returns The mapping of BaseMesh MeshPart vertices to PatchMeshPart MeshPart vertices
       */
      const std::map<Index, Index>& get_vertex_map() const
      {
        return BaseClass::get_vertex_map();
      }
    }; // class PatchPartMapHolder<...>

    /// \cond internal
    /**
     * \brief End of template recursion: Shape dimension 0
     */
    template<typename Shape_>
    class PatchPartMapHolder<Shape_, 0>
    {
    private:
      PatchPartMap _patch_map;

    public:
      explicit PatchPartMapHolder(const TargetSetHolder<Shape_>& tsh) :
        _patch_map(tsh.template get_target_set<0>())
        {
        }

      bool build(const TargetSetHolder<Shape_>& tsh, const IndexSetHolder<Shape_>* DOXY(ish))
      {
        return _patch_map.build(tsh.template get_target_set<0>());
      }

      Index get_num_entities(int dim) const
      {
        XASSERT(dim == 0);
        return _patch_map.size();
      }

      template<typename AttributeSetContainer_>
      void fill_attribute_sets(
        AttributeSetContainer_& attribute_container_out,
        const AttributeSetContainer_& attribute_container_in) const
      {
        // Fill attribute set to shape dimension 0
        _patch_map.fill_attribute_set( attribute_container_out, attribute_container_in);
      }

      template<typename IndexSetHolder_>
      void fill_index_sets(IndexSetHolder_& DOXY(ish), const IndexSetHolder_& DOXY(ish_in)) const
      {
      }

      void fill_target_sets(TargetSetHolder<Shape_>& tsh) const
      {
        _patch_map.fill_target_set(tsh.template get_target_set<0>());
      }


      const std::map<Index, Index>& get_vertex_map() const
      {
        return _patch_map.get_io_map();
      }
    }; // class PatchPartMapHolder<Shape_, 0>
    /// \endcond

  } // namespace Geometry
} // namespace FEAT
