// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/geometry/adaptive_mesh.hpp>
#include <kernel/geometry/adaptive_mesh_layer.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_permutation.hpp>
#include <kernel/geometry/subdivision_levels.hpp>
#include <kernel/util/dist.hpp>

#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <memory>
#include <optional>
#include <utility>

namespace FEAT::Control::Domain
{
  /**
   * \brief Wrapper class for domain levels for AdaptivePartiDomainControl
   *
   * Wrapper class that extends the given DomainLevel class by an additional mesh_layer member,
   * that gives access to the refinement tree of the adaptive mesh created by the AdaptivePartiDomainControl.
   *
   * \author Markus Muegge
   */
  template<typename DomainLevel_, typename TemplateSet_>
  class AdaptiveLevelWrapper : public DomainLevel_
  {
  public:
    /// our base-class; the actual domain level
    using BaseClass = DomainLevel_;
    /// our mesh type
    using typename BaseClass::MeshType;
    /// our mesh node type
    using typename BaseClass::MeshNodeType;
    /// our adaptive mesh type
    using AdaptiveMeshType = Geometry::AdaptiveMesh<TemplateSet_, typename MeshType::ShapeType>;
    /// our adaptive mesh layer type
    using AdaptiveMeshLayerType = Geometry::AdaptiveMeshLayer<AdaptiveMeshType>;

    /// View into adaptive mesh. If empty, this layer was created by a full refinement
    std::optional<AdaptiveMeshLayerType> mesh_layer;

    /// Contructor for fully refined levels
    explicit AdaptiveLevelWrapper(int _lvl_idx, std::unique_ptr<MeshNodeType> node) :
      BaseClass(_lvl_idx, std::move(node)),
      mesh_layer(std::nullopt)
    {
    }

    /// Constructor for partially refined levels
    explicit AdaptiveLevelWrapper(int _lvl_idx, std::unique_ptr<MeshNodeType> node, AdaptiveMeshLayerType&& layer) :
      BaseClass(_lvl_idx, std::move(node)),
      mesh_layer(std::make_optional<AdaptiveMeshLayerType>(std::move(layer)))
    {
    }

    /// Returns true if this level belongs a potentially partial refinement
    bool is_partial_level()
    {
      return bool(mesh_layer);
    }
  }; // class AdaptiveLevelWrapper<...>

  /**
   * \brief Adaptive Partitioned Domain Control
   *
   * \tparam DomainLevel_ Assumed to be a AdaptiveLevelWrapper around some other domain level type
   * \tparam TemplateSet_ The template set to use for partial refinement
   *
   * Extends the PartiDomainControl class with partial refinement of the mesh at the (former) finest level.
   *
   * The partial refinement is given by a lambda function, called the refinement strategy, passed to the domain
   * controller, which assigns subdivision levels to mesh vertices. The refinement strategy is also used to guide
   * partitioning, by creating weights for use by the partitioner. It may thus be called multiple times during domain
   * creation for different meshes. As a consequence the apriori 2 level partitioner is disabled by this class, because
   * it makes no use of these weights.
   *
   * The partial refinement is performed at the finest level produced by the normal PartiDomainControl class.
   * No partitioning of the partial refinement is performed.
   *
   * Supports the same command line/property map arguments as the normal PartiDomainControl
   *
   * \author Markus Muegge
   */
  template<typename DomainLevel_, typename TemplateSet_>
  class AdaptivePartiDomainControl : public PartiDomainControl<DomainLevel_>
  {
  public:
    /// Our base class
    using BaseClass = Control::Domain::PartiDomainControl<DomainLevel_>;
    /// our domain level type
    using typename BaseClass::LevelType;
    /// our domain layer type
    using typename BaseClass::LayerType;
    /// our mesh type
    using typename BaseClass::MeshType;
    /// our atlas type
    using typename BaseClass::AtlasType;
    /// our root mesh node type
    using typename BaseClass::MeshNodeType;
    /// our mesh-part type
    using typename BaseClass::MeshPartType;
    /// our ancestry type
    using typename BaseClass::Ancestor;
    /// our weight type
    using typename BaseClass::WeightType;

    using AdaptiveMeshType = Geometry::AdaptiveMesh<TemplateSet_, typename MeshType::ShapeType>;

  protected:
    /// the refinement strategy for the partially refined parts of the domain
    std::function<void(Geometry::SubdivisionLevels&, const MeshType&)> _refinement_strategy;

    std::shared_ptr<AdaptiveMeshType> _adaptive_mesh;

  public:
    /**
     * \brief Constructor
     *
     * \param[in] comm_
     * A \resident reference to the main communicator to be used.
     *
     * \param[in] support_multi_layered
     * Specifies whether the controller is allowed to create multi-layered hierarchies.
     */
    explicit AdaptivePartiDomainControl(
      const Dist::Comm& comm_,
      bool support_multi_layered,
      std::function<void(Geometry::SubdivisionLevels&, const MeshType&)> refinement_strategy) :
      BaseClass(comm_, support_multi_layered),
      _refinement_strategy(std::move(refinement_strategy))
    {
      this->_allow_parti_2level = false;
    }

    /// virtual destructor
    virtual ~AdaptivePartiDomainControl() = default;

    /// Accessor for underlying adaptive mesh
    const AdaptiveMeshType& get_adaptive_mesh() const
    {
      return *_adaptive_mesh;
    }

  protected:
    /// Compute weights for a-posteriori partitioners
    /**
     * \brief Compute weights for an a-posteriori partitioner
     *
     * Weights are computed by counting the number of vertex markings of each cell,
     * and multiplying that count by the average number of new elements created by each
     * marking in the template set.
     */
    std::vector<WeightType> _compute_weights(Ancestor& DOXY(ancestor), const MeshNodeType& base_mesh_node) override
    {
      static constexpr int dim = MeshType::ShapeType::dimension;
      const WeightType average_elements_per_marking =
        WeightType(TemplateSet_::template average_elements_per_marking<dim>());

      const MeshType& mesh = *base_mesh_node.get_mesh();
      const auto& v_at_c = mesh.template get_index_set<dim, 0>();

      Geometry::SubdivisionLevels sdls(mesh.get_num_vertices());
      std::vector<WeightType> weights(mesh.get_num_elements());

      // Get subdivision levels for mesh
      _refinement_strategy(sdls, mesh);

      // Determine weights from subdivision levels
      for(Index cell(0); cell < weights.size(); cell++)
      {
        std::uint64_t markings(0);
        for(auto i = v_at_c.image_begin(cell); i != v_at_c.image_end(cell); ++i)
        {
          markings += sdls[*i];
        }
        // Minimum weight is 1
        weights[cell] = Math::max(WeightType(1), WeightType(markings) * average_elements_per_marking);
      }

      return weights;
    }

    /**
     * \brief Creates a single-layered mesh hierarchy for a single process.
     *
     * \param[in] base_mesh_node
     * The base-mesh node from which the hierarchy is to be derived from.
     *
     * \note
     * This function does not keep the base-mesh levels explicitly even if
     * _keep_base_levels is set to true, because it is not required to use
     * the base splitter on a single process.
     */
    void _create_single_process(std::unique_ptr<MeshNodeType> base_mesh_node) override
    {
      // Perform regular refinement of conformal mesh
      BaseClass::_create_single_process(std::move(base_mesh_node));
      refine_partially();
    }

    // Note: all following member functions are only required for parallel builds,
    // so we enclose them in the following #if-block to reduce compile times.

#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)

    /**
     * \brief Creates a single-layered mesh hierarchy.
     *
     * \param[in] base_mesh_node
     * The base-mesh node from which the hierarchy is to be derived from.
     */
    void _create_single_layered(std::unique_ptr<MeshNodeType> base_mesh_node) override
    {
      BaseClass::_create_single_layered(std::move(base_mesh_node));
      refine_partially();
    }

    /**
     * \brief Creates a multi-layered mesh hierarchy.
     *
     * \param[in] base_mesh_node
     * The base-mesh node from which the hierarchy is to be derived from.
     */
    void _create_multi_layered(std::unique_ptr<MeshNodeType> base_mesh_node) override
    {
      BaseClass::_create_multi_layered(std::move(base_mesh_node));
      refine_partially();
    }

#endif // defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
    /**
     * \brief Create partially refined levels on top of the current finest level
     *
     * Creates an adaptive mesh from the currently finest level, performs a globally consistent
     * partial refinement, and pushes the resulting meshes as new level to the domain.
     * The final partial level is recorded as a new chosen level.
     */
    void refine_partially()
    {
      const MeshNodeType& mesh_node = *this->_layer_levels.at(0).front()->get_mesh_node();
      Geometry::SubdivisionLevels sdls(mesh_node.get_mesh()->get_num_vertices());

      // Get requested subdivision levels
      _refinement_strategy(sdls, *mesh_node.get_mesh());

      // Sync subdivision levels across processes
      sync_subdivision_levels(sdls, mesh_node);

      // Sync maximum subdivision level across ranks
      // This value will determine how many copy-layers we need to add
      const std::uint64_t sdl_max_local = sdls.maximum_level();
      std::uint64_t sdl_max;
      this->comm().allreduce(&sdl_max_local, &sdl_max, std::size_t(1), Dist::op_max);

      // Create adaptive mesh
      _adaptive_mesh = std::make_shared<AdaptiveMeshType>(*mesh_node.get_mesh());
      _adaptive_mesh->adapt(sdls, AdaptiveMeshType::ImportBehaviour::All);

      int lvl = this->_chosen_levels.front().first;
      for(Index i(1); i <= sdl_max; i++)
      {
        // Either push the next layer of the adaptive mesh or a copy of the last layer.
        // We push a copy to ensure we correctly participate in global operations on all levels.
        // TODO: Use a more efficient method, rather than copying the previous level
        Geometry::Layer layer_idx{Math::min(i, _adaptive_mesh->num_layers() - 1)};
        Geometry::AdaptiveMeshLayer<AdaptiveMeshType> layer(_adaptive_mesh, layer_idx);
        this->push_level_front(
          0,
          std::make_shared<LevelType>(
            lvl++,
            project_mesh_node(*_adaptive_mesh, layer_idx, mesh_node),
            std::move(layer)));
      }

      this->_chosen_levels.emplace_front(lvl, this->_comm.size());
    }

    /**
     * \brief Project a MeshNode from a regular mesh onto a layer of an adaptive mesh
     *
     * \param[in] mesh Adaptive mesh build on the mesh contained in root_node
     * \param[in] l Layer of the adaptive mesh to project onto
     * \param[in] root_node Mesh node to project
     * \param[in] adapt_mode
     *
     * This mesh takes a RootMeshNode containing some ConformalMesh and an adaptive mesh
     * build on top of the ConformalMesh. It then projects the mesh-parts, halos, and patches
     * of the MeshNode onto chosen layer of the adaptive mesh, by replacing all entries by
     * their (grand-)children in the adaptive mesh.
     *
     * \note
     * Adapting the resulting mesh node to a chart is supported.
     * Dual adapation is not supported.
     */
    std::unique_ptr<MeshNodeType> project_mesh_node(
      const AdaptiveMeshType& mesh,
      Geometry::Layer l,
      const MeshNodeType& root_node,
      Geometry::AdaptMode adapt_mode = Geometry::AdaptMode::chart)
    {

      auto cmesh = std::make_unique<MeshType>(mesh.to_conformal_mesh(l));
      auto result = std::make_unique<MeshNodeType>(std::move(cmesh));

      // Projects mesh parts
      for(const auto& name : root_node.get_mesh_part_names())
      {
        const auto& mesh_part_node = root_node.find_mesh_part_node(name);
        const auto* mesh_part = mesh_part_node->get_mesh();

        std::unique_ptr<MeshPartType> new_mesh_part;
        if(mesh_part)
        {
          new_mesh_part = std::make_unique<MeshPartType>(mesh.template project_meshpart<MeshType>(l, *mesh_part));
        }
        else
        {
          new_mesh_part = std::unique_ptr<MeshPartType>(nullptr);
        }
        auto new_mesh_part_node = std::make_unique<Geometry::MeshPartNode<MeshType>>(std::move(new_mesh_part));

        result->add_mesh_part_node(
          name,
          std::move(new_mesh_part_node),
          mesh_part_node->find_mesh_part_chart_name(name),
          mesh_part_node->find_mesh_part_chart(name));
      }

      // refine our halos
      for(const auto& [rank, part] : root_node.get_halo_map())
      {
        if(part)
        {
          XASSERT(bool(part));
          auto new_halo = std::make_unique<MeshPartType>(mesh.template project_meshpart<MeshType>(l, *part));
          result->add_halo(rank, std::move(new_halo));
        }
        else
        {
          result->add_halo(rank, nullptr);
        }
      }

      // refine our patch mesh-parts
      for(const auto& [rank, patch] : root_node.get_patch_map())
      {
        if(patch)
        {
          XASSERT(bool(patch));
          auto new_patch = std::make_unique<MeshPartType>(mesh.template project_meshpart<MeshType>(l, *patch));
          result->add_patch(rank, std::move(new_patch));
        }
        else
        {
          result->add_patch(rank, nullptr);
        }
      }

      // adapt by chart?
      if((adapt_mode & Geometry::AdaptMode::chart) != Geometry::AdaptMode::none)
      {
        result->adapt(true);
      }

      // TODO: Adapt dual

      return result;
    }

    /**
     * \brief Sync subdivision levels across processes by taking the maximum across processes
     *
     * During domain creation each rank created its own subidivison levels for its patch.
     * Because each rank only considers its own local patch of the mesh when creating its
     * subdivision levels, these markings might differ at halo vertices. This method
     * synchronizes the values at these shared vertices by computing the maximum for each
     * shared vertex. This way all ranks obtain the same markings for all shared vertices.
     *
     * This is a precondition for creating halos on the partially refined levels.
     * Halo creation there works by replacing the entries of the halos at the finest level
     * by all their (grand-)children on the partially refined levels.
     * Halos created this way are only globally consistent if all ranks perform the same
     * refinement along shared vertices.
     */
    void sync_subdivision_levels(Geometry::SubdivisionLevels& sdls, const MeshNodeType& mesh_node)
    {
      XASSERT(mesh_node.get_mesh()->get_num_vertices() == sdls.size());

      std::unordered_map<int, Dist::RequestVector> sends;
      std::unordered_map<int, Dist::RequestVector> receives;

      Geometry::SubdivisionLevels received_sdls(sdls.size(), 0);

      // Post send and receives
      for(const auto& [rank, halo] : mesh_node.get_halo_map())
      {
        const auto& vertices = halo->template get_target_set<0>();

        for(Index i(0); i < vertices.get_num_entities(); i++)
        {
          sends[rank].push_back(this->_comm.isend(&sdls[vertices[i]], 1, rank));
          receives[rank].push_back(this->_comm.irecv(&received_sdls[vertices[i]], 1, rank));
        }
      }

      // Wait for all receives and map sdls
      for(const auto& [rank, halo] : mesh_node.get_halo_map())
      {
        const auto& vertices = halo->template get_target_set<0>();
        for(Index idx(0); receives[rank].wait_any(idx);)
        {
          const Index vidx = vertices[idx];
          // Subdivision level is maximum of all processes values
          sdls[vidx] = Math::max(sdls[vidx], received_sdls[vidx]);
        }
      }

      // Wait for all sends to finish
      for(const auto& [rank, halo] : mesh_node.get_halo_map())
      {
        sends[rank].wait_all();
      }
    }
  };
} // namespace FEAT::Control::Domain
