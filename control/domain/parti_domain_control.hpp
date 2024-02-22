// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef CONTROL_DOMAIN_PARTI_DOMAIN_CONTROL_HPP
#define CONTROL_DOMAIN_PARTI_DOMAIN_CONTROL_HPP 1

#include <kernel/base_header.hpp>

#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/common_factories.hpp>

#include <control/domain/parti_domain_control_base.hpp>

namespace FEAT
{
  namespace Control
  {
    namespace Domain
    {
      /**
       * \brief Adds the supported arguments of the PartiDomainControl to an argument parser.
       *
       * This function adds all supported options, which can be parsed by calling
       * the PartiDomainControl::parse_args() function, to a SimpleArgParser object.
       *
       * \param[in] args
       * The argument parser.
       */
      inline void add_supported_pdc_args(SimpleArgParser& args)
      {
        args.support("parti-type", "<types...>\n"
          "Specifies which partitioner types are allowed to be used.\n"
          "May contain the following types:\n"
          "2level extern genetic metis naive zoltan"
        );
        args.support("parti-extern-name", "<names...>\n"
          "Specifies the names of the allowed extern partitions."
        );
        args.support("parti-rank-elems", "<count>\n"
          "Specifies the minimum number of elements per rank for a-posteriori partitioning."
        );
        args.support("parti-genetic-time", "<time-init> <time-mutate>\n"
          "Specifies the time for initial distribution and mutation for the genetic partitioner."
        );
      }

      /**
       * \brief Recursively Partitioned Domain Control
       *
       * \todo document this
       *
       * For more details on meshes, see the related doxygen page \ref mesh_file_format.
       *
       * \author Peter Zajac
       */
      template<typename DomainLevel_>
      class PartiDomainControl :
        public Control::Domain::PartiDomainControlBase<DomainLevel_>
      {
      public:
        /// Our base class
        typedef Control::Domain::PartiDomainControlBase<DomainLevel_> BaseClass;
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

      protected:
        /// the extern partition sets
        Geometry::PartitionSet _parti_set;

        /// allow extern partitioner?
        bool _allow_parti_extern;
        /// allow 2-level partitioner?
        bool _allow_parti_2level;

        /// required partition name for extern partitioning
        std::deque<String> _extern_parti_names;

        /// the permutation strategy for the mesh permutation
        Geometry::PermutationStrategy _permutation_strategy;

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
        explicit PartiDomainControl(const Dist::Comm& comm_, bool support_multi_layered) :
          BaseClass(comm_, support_multi_layered),
          _allow_parti_extern(true),
          _allow_parti_2level(true),
          _extern_parti_names(),
          _permutation_strategy(Geometry::PermutationStrategy::none)
        {
        }

        /// virtual destructor
        virtual ~PartiDomainControl()
        {
        }

        /**
         * \brief Parses the partitioner options from an argument parser.
         *
         * \see Control::Domain::add_supported_pdc_args()
         *
         * \param[in] args
         * The parser that is to be used.
         *
         * \returns
         * \c true, if the parsing was successful, or \c false,
         * if at least one option was invalid.
         */
        virtual bool parse_args(SimpleArgParser& args) override
        {
          if(!BaseClass::parse_args(args))
            return false;

          // parse --parti-extern-name <names...>
          {
            auto it = args.query("parti-extern-name");
            if(it != nullptr)
              _extern_parti_names = it->second;
          }

          // okay
          return true;
        }

        /**
         * \brief Parses the partitioner options from a PropertyMap
         *
         * \param[in] pmap
         * A \transient reference to the property map that contains the configuration
         *
         * \returns
         * \c true, if the parsing was successful, or \c false,
         * if at least one option was invalid.
         */
        virtual bool parse_property_map(PropertyMap& pmap) override
        {
          if(!BaseClass::parse_property_map(pmap))
            return false;

          auto parti_extern_name_p = pmap.query("parti-extern-name");
          if(parti_extern_name_p.second)
          {
            this->_extern_parti_names = parti_extern_name_p.first.split_by_whitespaces();
          }

          return true;
        }

        /**
         * \brief Creates a domain control from a list of filenames.
         *
         * \param[in] filenames
         * A list of mesh filenames from which the domain control is to be created.
         *
         * \param[in] dirpath
         * The path in which the mesh files are located.
         */
        virtual void create(const std::deque<String>& filenames, String dirpath = "")
        {
          // create a mesh file reader
          Geometry::MeshFileReader mesh_reader;

          // add the files
          mesh_reader.add_mesh_files(this->_comm, filenames, dirpath);

          // finally, create the control
          this->create(mesh_reader);
        }

        /**
         * \brief Creates a domain control from a MeshFileReader.
         *
         * \param[inout] mesh_reader
         * A \transient reference to the mesh reader from which the domain control is to be created.
         */
        virtual void create(Geometry::MeshFileReader& mesh_reader)
        {
          // ensure that the domain control is still empty
          XASSERTM(!this->_was_created, "domain has already been created");
          XASSERT(this->size_physical() == std::size_t(0));
          XASSERT(this->size_virtual() == std::size_t(0));

          TimeStamp stamp_create;

          // parse mesh and create control
          this->create(mesh_reader.parse(this->_atlas, &_parti_set));

          // clear partition set after creation
          _parti_set.clear();

          // collect partition statistics
          FEAT::Statistics::toe_partition = stamp_create.elapsed_now();
        }

        /**
         * \brief Creates a domain control from a base-mesh node
         *
         * \param[in] base_mesh_node
         * The base-mesh node that the domain control is to be created from
         */
        virtual void create(std::unique_ptr<MeshNodeType> base_mesh_node)
        {
          TimeStamp stamp_create;

          // call base-class create
          this->_create(std::move(base_mesh_node));

          // finally, apply mesh permutation
          this->_was_created = false;
          this->create_mesh_permutations();
          this->_was_created = true;

          // collect partition statistics
          FEAT::Statistics::toe_partition = stamp_create.elapsed_now();
        }

        /**
         * \brief Creates a domain control from a structured rectilinear base mesh
         *
         * \note This function also creates a mesh-part named "bnd" that represents the entire domain boundary.
         *
         * \param[in] num_elems_x, num_elems_y, num_elems_z
         * The number of elements in each dimension for the base-mesh
         */
        virtual void create_rectilinear(Index num_elems_x, Index num_elems_y = 0u, Index num_elems_z = 0u)
        {
          // ensure that the domain control is still empty
          XASSERTM(!this->_was_created, "domain has already been created");
          XASSERT(this->size_physical() == std::size_t(0));
          XASSERT(this->size_virtual() == std::size_t(0));

          TimeStamp stamp_create;

          // allocate a base-mesh node and create mesh
          std::unique_ptr<MeshNodeType> base_mesh_node = MeshNodeType::make_unique(
            Geometry::StructUnitCubeFactory<MeshType>::make_unique_from(num_elems_x, num_elems_y, num_elems_z));

          // get mesh and create boundary
          {
            Geometry::BoundaryFactory<MeshType> bnd_factory(*base_mesh_node->get_mesh());
            base_mesh_node->add_mesh_part("bnd", bnd_factory.make_unique());
          }

          // create control
          this->create(std::move(base_mesh_node));

          // collect partition statistics
          FEAT::Statistics::toe_partition = stamp_create.elapsed_now();
        }

        /**
         * \brief Sets the permutation strategy for mesh permutation.
         *
         * \param[in] strategy
         * The mesh permutation strategy to be used.
         */
        void set_permutation_strategy(Geometry::PermutationStrategy strategy)
        {
          XASSERTM(!this->_was_created, "This function has to be called before domain control creation!");
          _permutation_strategy = strategy;
        }

        /// \returns The permutation strategy used for mesh permutation.
        Geometry::PermutationStrategy get_permutation_strategy() const
        {
          return _permutation_strategy;
        }

        /**
         * \brief Creates the mesh permutations based on the chosen permutation strategy.
         */
        virtual void create_mesh_permutations()
        {
          // nothing to do?
          if(_permutation_strategy == Geometry::PermutationStrategy::none)
            return;

          // loop over all levels and permute
          for(std::size_t i(0); i < this->size_physical(); ++i)
          {
            this->at(i)->get_mesh_node()->create_permutation(this->_permutation_strategy);
          }
        }

      protected:
        /**
         * \brief Resets/disables all partitioner types
         *
         * This function is called by the parse_args and parse_property_map function when the user
         * has supplied a list of allowed partitioner types.
         */
        virtual void _reset_parti_types() override
        {
          _allow_parti_extern = _allow_parti_2level = false;
          BaseClass::_reset_parti_types();
        }

        /**
         * \brief Parses a partitioner type
         *
         * \param[in] type
         * The name of the partitioner type that is to be enabled.
         *
         * \returns \c true, if \p type represents a supported partitioner type, otherwise \c false.
         */
        virtual bool _parse_parti_type(const String& type) override
        {
          if(type.compare_no_case("extern") == 0)
            return _allow_parti_extern = true;
          else if(type.compare_no_case("2level") == 0)
            return _allow_parti_2level = true;
          else
            return BaseClass::_parse_parti_type(type);
        }

        /**
         * \brief Creates a domain control from a base-mesh node
         *
         * \param[in] base_mesh_node
         * The base-mesh node that the domain control is to be created from
         */
        virtual void _create(std::unique_ptr<MeshNodeType> base_mesh_node)
        {
          // ensure that the domain control is still empty
          XASSERTM(!this->_was_created, "domain has already been created");
          XASSERT(this->size_physical() == std::size_t(0));
          XASSERT(this->size_virtual() == std::size_t(0));

          TimeStamp stamp_create;

          // create the domain control
#ifdef FEAT_HAVE_MPI
          if(this->_comm.size() == 1)
          {
            // We've got just one process, so it's a simple choice:
            this->_create_single_process(std::move(base_mesh_node));
          }
          else if(this->_support_multi_layered && (this->_desired_levels.size() > std::size_t(2)))
          {
            // The application supports multi-layered domain controls and
            // the user wants that, so create a multi-layered one:
            this->_create_multi_layered(std::move(base_mesh_node));
          }
          else
          {
            // Create a single-layered domain control:
            this->_create_single_layered(std::move(base_mesh_node));
          }
#else // not FEAT_HAVE_MPI
          {
            // In the non-MPI case, we always have only one process:
            this->_create_single_process(std::move(base_mesh_node));
          }
#endif // FEAT_HAVE_MPI

          // compile virtual levels
          this->compile_virtual_levels();

          // collect partition statistics
          FEAT::Statistics::toe_partition = stamp_create.elapsed_now();

          // domain created
          this->_was_created = true;
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
        virtual void _create_single_process(std::unique_ptr<MeshNodeType> base_mesh_node)
        {
          // create and push single layer
          this->push_layer(std::make_shared<LayerType>(this->_comm.comm_dup(), 0));

          // create single ancestry
          this->_create_ancestry_single();
          Ancestor& ancestor = this->_ancestry.front();

          // use the whole base mesh here
          ancestor.parti_info = "Using base-mesh";

          // refine up to desired minimum level
          int lvl = 0;
          for(; lvl < ancestor.desired_level_min; ++lvl)
          {
            // refine the patch mesh
            base_mesh_node = base_mesh_node->refine_unique(this->_adapt_mode);
          }

          // save chosen minimum level if it is not equal to the desired maximum level
          if(lvl < ancestor.desired_level_max)
            this->_chosen_levels.push_front(std::make_pair(lvl, 0));

          // refine up to maximum level and push to control
          for(; lvl < ancestor.desired_level_max; ++lvl)
          {
            // refine the base mesh
            auto refined_node = base_mesh_node->refine_unique(this->_adapt_mode);

            // push this level
            this->push_level_front(0, std::make_shared<LevelType>(lvl, std::move(base_mesh_node)));

            // continue with refined node
            base_mesh_node = std::move(refined_node);
          }

          // save chosen maximum level
          this->_chosen_levels.push_front(std::make_pair(lvl, 1));

          // push finest level
          this->push_level_front(0, std::make_shared<LevelType>(lvl, std::move(base_mesh_node)));
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
        virtual void _create_single_layered(std::unique_ptr<MeshNodeType> base_mesh_node)
        {
          // create and push single layer
          std::shared_ptr<LayerType> layer = std::make_shared<LayerType>(this->_comm.comm_dup(), 0);
          this->push_layer(layer);

          // create single-layered ancestry
          this->_create_ancestry_single();
          Ancestor& ancestor = this->_ancestry.front();

          // do we have to keep the base-mesh levels on this process?
          const bool keep_base = this->_keep_base_levels && (this->_comm.rank() == 0);

          // choose a partitioning strategy
          this->_check_parti(ancestor, *base_mesh_node, true);

          // refine up to chosen partition level
          int lvl = 0;
          for(; lvl < ancestor.parti_level; ++lvl)
          {
            // refine the base mesh
            base_mesh_node = base_mesh_node->refine_unique(this->_adapt_mode);
          }

          // apply partitioner
          if(!this->_apply_parti(ancestor, *base_mesh_node))
          {
            XABORTM("Failed to find a suitable partitioning");
          }

          // extract our patch
          std::vector<int> neighbor_ranks;
          std::unique_ptr<MeshNodeType> patch_mesh_node(
            base_mesh_node->extract_patch(neighbor_ranks, ancestor.parti_graph, this->_comm.rank()));

          // set the neighbor ranks of our child layer
          layer->set_neighbor_ranks(neighbor_ranks);

          // if we have to keep the base-mesh levels, then we also need the patches
          if(keep_base)
          {
            // Note: patch mesh-part for rank = 0 was already created by 'extract_patch' call
            for(int i(1); i < ancestor.num_parts; ++i)
            {
              base_mesh_node->create_patch_meshpart(ancestor.parti_graph, i);
            }
          }

          // refine up to minimum level
          for(; lvl < ancestor.desired_level_min; ++lvl)
          {
            // refine the patch mesh
            patch_mesh_node = patch_mesh_node->refine_unique(this->_adapt_mode);

            // keep base mesh node?
            if(keep_base)
              base_mesh_node = base_mesh_node->refine_unique(this->_adapt_mode);
          }

          // save chosen minimum level
          if(lvl < ancestor.desired_level_max)
            this->_chosen_levels.push_front(std::make_pair(lvl, 0));

          // refine up to maximum level
          for(; lvl < ancestor.desired_level_max; ++lvl)
          {
            // refine the patch mesh
            auto refined_node = patch_mesh_node->refine_unique(this->_adapt_mode);

            // push this (unrefined) level
            this->push_level_front(0, std::make_shared<LevelType>(lvl, std::move(patch_mesh_node)));

            // continue with refined node
            patch_mesh_node = std::move(refined_node);

            // keep base mesh node?
            if(keep_base)
            {
              // refine base-mesh node
              auto refined_base_mesh_node = base_mesh_node->refine_unique(this->_adapt_mode);

              // push this (unrefined) level
              this->_base_levels.push_front(std::make_shared<LevelType>(lvl, std::move(base_mesh_node)));

              // continue with refined node
              base_mesh_node = std::move(refined_base_mesh_node);
            }
          }

          // save chosen maximum level
          this->_chosen_levels.push_front(std::make_pair(lvl, this->_comm.size()));

          // push finest level
          this->push_level_front(0, std::make_shared<LevelType>(lvl, std::move(patch_mesh_node)));

          // keep base mesh node?
          if(keep_base)
            this->_base_levels.push_front(std::make_shared<LevelType>(lvl, std::move(base_mesh_node)));
        }

        /**
         * \brief Creates a multi-layered mesh hierarchy.
         *
         * \param[in] base_mesh_node
         * The base-mesh node from which the hierarchy is to be derived from.
         */
        virtual void _create_multi_layered(std::unique_ptr<MeshNodeType> base_mesh_node)
        {
          // create layers
          this->_create_multi_layers_scattered();

          // create ancestry
          this->_create_ancestry_scattered();

          // do we have to keep the base-mesh levels on this process?
          const bool keep_base = this->_keep_base_levels && (this->_comm.rank() == 0);

          // make sure we have at most 2 layers if we want to keep the base-mesh levels
          XASSERTM(!keep_base || (this->_ancestry.size() <= std::size_t(1)), "cannot keep base levels for more than 2 domain layers (yet)");

          // we start counting at level 0
          int lvl = 0;

          // move the base-mesh pointer to a new parent-mesh pointer;
          // the base-mesh is always the one on the layer with only 1 process
          std::unique_ptr<MeshNodeType> parent_mesh_node = std::move(base_mesh_node);
          base_mesh_node.reset();

          // loop over all global layers in reverse order (coarse to fine)
          for(std::size_t slayer = this->_ancestry.size(); slayer > std::size_t(0); )
          {
            // is this the base-mesh layer aka the 1-process layer?
            const bool is_base_layer = (slayer == this->_ancestry.size());

            --slayer;

            // get the ancestor object
            Ancestor& ancestor = this->_ancestry.at(slayer);

            // determine the minimum desired level of our parent layer
            int parent_min_lvl = -1;
            if(!is_base_layer)
              parent_min_lvl = this->_chosen_levels.front().first;
            else if(this->_ancestry.size() + 1u < this->_desired_levels.size())
              parent_min_lvl = this->_desired_levels.back().first;

            // check available partitioning strategies
            this->_check_parti(ancestor, *parent_mesh_node, is_base_layer);

            // the check_parti function returns the partitioning level w.r.t. the current
            // level (which may be > 0), so we have to compensate that by adding our current level:
            ancestor.parti_level += lvl;

            // Note: each progeny group within the main communicator may have chosen a different
            // partitioning level at this point. We will compensate this by adjusting the minimum
            // refinement level of the child layer after the partitioning step below.

            // refine up to the chosen partitioning level
            for(; lvl < ancestor.parti_level; ++lvl)
            {
              // refine the base-mesh node
              auto refined_node = parent_mesh_node->refine_unique(this->_adapt_mode);

              // push the base mesh into our parent layer if desired
              if((ancestor.layer_p >= 0) && (parent_min_lvl >= 0) && (lvl >= parent_min_lvl))
              {
                this->push_level_front(ancestor.layer_p, std::make_shared<LevelType>(lvl, std::move(parent_mesh_node)));
              }

              // continue with refined node
              parent_mesh_node = std::move(refined_node);

              // refine base-mesh node if we have one
              if(base_mesh_node)
              {
                // refine base-mesh node
                auto refined_base_mesh_node = base_mesh_node->refine_unique(this->_adapt_mode);

                // push this (unrefined) level
                this->_base_levels.push_front(std::make_shared<LevelType>(lvl, std::move(base_mesh_node)));

                // continue with refined node
                base_mesh_node = std::move(refined_base_mesh_node);
              }
            }

            // we're now at the partitioning level, so apply the partitioner
            if(!this->_apply_parti(ancestor, *parent_mesh_node))
            {
              XABORTM("Failed to find a suitable partitioning");
            }

            // extract our patch
            std::vector<int> neighbor_ranks;
            std::unique_ptr<MeshNodeType> patch_mesh_node(
              parent_mesh_node->extract_patch(neighbor_ranks, ancestor.parti_graph, ancestor.progeny_child));

            // translate neighbor ranks by progeny group to obtain the neighbor ranks
            // w.r.t. this layer's communicator
            {
              std::map<int,int> halo_map;
              for(auto& i : neighbor_ranks)
              {
                int old_i(i);
                halo_map.emplace(old_i, i += ancestor.progeny_group);
              }
              patch_mesh_node->rename_halos(halo_map);
            }

            // does this process participate in the parent layer or do we need to keep the base-meshes anyways?
            if((ancestor.layer_p >= 0) || (keep_base && is_base_layer))
            {
              // Note: patch mesh-part for rank = 0 was already created by 'extract_patch' call
              for(int i(1); i < ancestor.num_parts; ++i)
              {
                parent_mesh_node->create_patch_meshpart(ancestor.parti_graph, i);
              }
            }

            // make sure we choose the same minimum level for all processes, because we may
            // have chosen different partitioning levels for each patch
            int global_level_min = Math::max(ancestor.desired_level_min, ancestor.parti_level);
            this->_comm.allreduce(&global_level_min, &global_level_min, std::size_t(1), Dist::op_max);

            // make sure our minimum level is greater than the minimum level of the previous layer,
            // because each layer must contain at least one non-ghost level
            if(!is_base_layer)
              global_level_min = Math::max(global_level_min, this->_chosen_levels.front().first+1);

            // refine up to desired minimum level of this layer
            for(; lvl < global_level_min; ++lvl)
            {
              // refine the base mesh node
              std::unique_ptr<MeshNodeType> coarse_mesh_node(std::move(parent_mesh_node));
              parent_mesh_node = coarse_mesh_node->refine_unique(this->_adapt_mode);

              // push base mesh to parent layer if desired
              if((ancestor.layer_p >= 0) && (parent_min_lvl >= 0) && (lvl >= parent_min_lvl))
              {
                // clear patches before pushing this node as they are redundant here
                coarse_mesh_node->clear_patches();
                this->push_level_front(ancestor.layer_p, std::make_shared<LevelType>(lvl, std::move(coarse_mesh_node)));
              }

              // refine the patch mesh
              patch_mesh_node = patch_mesh_node->refine_unique(this->_adapt_mode);

              // refine base-mesh node if we have one
              if(base_mesh_node)
              {
                // refine base-mesh node
                auto refined_base_mesh_node = base_mesh_node->refine_unique(this->_adapt_mode);

                // push this (unrefined) level
                this->_base_levels.push_front(std::make_shared<LevelType>(lvl, std::move(base_mesh_node)));

                // continue with refined node
                base_mesh_node = std::move(refined_base_mesh_node);
              }
            }

            // split the halos of our base-mesh and compute the halos of our patches from that
            this->_split_basemesh_halos(ancestor, *parent_mesh_node, *patch_mesh_node, neighbor_ranks);

            // does this process participate in the child layer?
            if(ancestor.layer >= 0)
            {
              // set the neighbor ranks in our child layer
              this->_layers.at(std::size_t(ancestor.layer))->set_neighbor_ranks(neighbor_ranks);
            }

            // set chosen minimum level for this layer
            if(!is_base_layer)
              this->_chosen_levels.push_front(std::make_pair(lvl, this->_ancestry.at(slayer+1u).num_procs));
            else if(parent_min_lvl < 0)
              this->_chosen_levels.push_front(std::make_pair(lvl, 0));
            else
            {
              this->_chosen_levels.push_front(std::make_pair(parent_min_lvl, 0));
              this->_chosen_levels.push_front(std::make_pair(lvl, 1));
            }

            // if we have to keep the base-mesh, then we have to create a clone here if we're still on the base layer;
            // we don't need to worry that this is a clone, because it will get refined in the next layer anyways
            if(keep_base && is_base_layer)
            {
              XASSERT(base_mesh_node.get() == nullptr);
              base_mesh_node = parent_mesh_node->clone_unique();
            }

            // push the finest base-mesh
            if(ancestor.layer_p >= 0)
            {
              this->push_level_front(ancestor.layer_p, std::make_shared<LevelType>(lvl, std::move(parent_mesh_node)));
            }

            // continue with the next layer
            parent_mesh_node = std::move(patch_mesh_node);
          }

          // get the desired maximum level
          // if we have more than one layer, make sure that the finest one contains at
          // least one level, as otherwise the finest global level would be a ghost level
          int desired_level_max = Math::max(this->_ancestry.front().desired_level_max, lvl+1);

          for(; lvl < desired_level_max; ++lvl)
          {
            // refine the patch mesh
            auto refined_node = parent_mesh_node->refine_unique(this->_adapt_mode);

            // push patch mesh to this level
            this->push_level_front(0, std::make_shared<LevelType>(lvl, std::move(parent_mesh_node)));

            // continue with refined node
            parent_mesh_node = std::move(refined_node);

            // refine base-mesh node if we have one
            if(base_mesh_node)
            {
              // refine base-mesh node
              auto refined_base_mesh_node = base_mesh_node->refine_unique(this->_adapt_mode);

              // push this (unrefined) level
              this->_base_levels.push_front(std::make_shared<LevelType>(lvl, std::move(base_mesh_node)));

              // continue with refined node
              base_mesh_node = std::move(refined_base_mesh_node);
            }
          }

          // push finest level
          this->push_level_front(0, std::make_shared<LevelType>(lvl, std::move(parent_mesh_node)));

          // set chosen maximum level for finest layer
          this->_chosen_levels.push_front(std::make_pair(lvl, this->_comm.size()));

          // refine base-mesh node if we have one
          if(base_mesh_node)
            this->_base_levels.push_front(std::make_shared<LevelType>(lvl, std::move(base_mesh_node)));
        }

        /**
         * \brief Checks for an appropriate partitioning strategy.
         *
         * This function examines the given base mesh and checks whether
         * one of the a-priori partitioning strategies (extern or 2-level)
         * yields a valid partitioning for the given number of partitions.
         * If so, the corresponding partitioning level and graph are stored
         * in ancestor.parti_level and ancestor.parti_graph.
         * If not, this function determines how often the base-mesh has to be
         * refined until its has enough cells so that one of the a-posteriori
         * partitioners can be applied.
         *
         * \param[inout] ancestor
         * The ancestor object for this layer.
         *
         * \param[in] mesh_node
         * The base-mesh node that is to be partitioned.
         *
         * \param[in] check_extern
         * Set to \c true, if #base_mesh_node is the mesh node that corresponds to the input file's mesh
         * and therefore extern partitioning is a valid option, or set to \c false in any other case.
         *
         * \returns
         * \c true, if an a-priori partitioning was found, otherwise \c false.
         */
        virtual bool _check_parti(Ancestor& ancestor, const MeshNodeType& mesh_node, bool is_base_layer) override
        {
          // Try to find an appropriate a-priori partitioning first:
          if(is_base_layer && this->_check_parti_extern(ancestor))
            return true;

          if(this->_check_parti_2level(ancestor, mesh_node))
            return true;

          // call base class
          return BaseClass::_check_parti(ancestor, mesh_node, is_base_layer);
        }

        /**
         * \brief Checks whether an extern partition is given
         *
         * \param[inout] ancestor
         * The ancestor object for this layer.
         *
         * \returns
         * \c true, if an extern partition is given, otherwise \c false.
         */
        bool _check_parti_extern(Ancestor& ancestor)
        {
          // is this even allowed?
          if(!this->_allow_parti_extern)
            return false;

          // check whether we have a suitable partition
          const Geometry::Partition* part = this->_parti_set.find_partition(ancestor.num_parts, _extern_parti_names);
          if(part == nullptr)
            return false;

          // set our ancestor partitioning
          ancestor.parti_apriori = true;
          ancestor.parti_info = String("Found extern partition '") + part->get_name() + "'";
          ancestor.parti_level = part->get_level();
          ancestor.parti_graph = part->get_patches().clone();
          return true;
        }

        /**
         * \brief Checks whether the 2-level partitioner can be applied.
         *
         * \param[inout] ancestor
         * The ancestor object for this layer.
         *
         * \param[in] base_mesh_node
         * The base-mesh node that is to be partitioned.
         *
         * \returns
         * \c true, if the 2-level partitioner can be applied, otherwise \c false.
         */
        bool _check_parti_2level(Ancestor& ancestor, const MeshNodeType& base_mesh_node)
        {
          // is this even allowed?
          if(!this->_allow_parti_2level)
            return false;

          // create a 2-level partitioner
          Geometry::Parti2Lvl<MeshType> partitioner(*base_mesh_node.get_mesh(), Index(ancestor.num_parts));

          // successful?
          if(!partitioner.success())
            return false;

          // found a valid 2-level partitioning
          ancestor.parti_apriori = true;
          ancestor.parti_info = String("Found 2-level partition");
          ancestor.parti_level = int(partitioner.parti_level());
          ancestor.parti_graph = partitioner.build_elems_at_rank();
          return true;
        }
#endif // defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
      }; // class PartiDomainControl<...>
    } // namespace Domain
  } // namespace Control
} // namespace FEAT

#endif // CONTROL_DOMAIN_PARTI_DOMAIN_CONTROL_HPP
