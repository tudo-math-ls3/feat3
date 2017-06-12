#pragma once
#ifndef CONTROL_DOMAIN_PARTI_DOMAIN_CONTROL_HPP
#define CONTROL_DOMAIN_PARTI_DOMAIN_CONTROL_HPP 1

#include <kernel/base_header.hpp>

#include <kernel/util/comm_base.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/dist_file_io.hpp>
#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/foundation/pexecutor.hpp>
#include <kernel/foundation/pgraph.hpp>
#include <kernel/foundation/psynch.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/partition_set.hpp>
#include <kernel/geometry/parti_2lvl.hpp>
#include <kernel/geometry/parti_iterative.hpp>

#include <control/domain/domain_control.hpp>

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
      void add_supported_pdc_args(SimpleArgParser& args)
      {
        args.support("parti-type", "<types...>\n"
          "Specifies which partitioner types are allowed to be used.\n"
          "Can contain the following types:\n"
          "2level"
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
       * \brief Partitioned Domain Control
       *
       * \todo document this
       *
       * \author Peter Zajac
       */
      template<typename DomainLevel_>
      class PartiDomainControl :
        public Control::Domain::DomainControl<DomainLevel_>
      {
      public:
        /// Our base class
        typedef Control::Domain::DomainControl<DomainLevel_> BaseClass;
        /// our domain level type
        typedef typename BaseClass::LevelType LevelType;
        /// our domain layer type
        typedef typename BaseClass::LayerType LayerType;
        /// our mesh type
        typedef typename BaseClass::MeshType MeshType;
        /// our atlas type
        typedef typename BaseClass::AtlasType AtlasType;
        /// our root mesh node type
        typedef Geometry::RootMeshNode<MeshType> MeshNodeType;

      protected:
        /// specifies whether the domain control was already created
        bool _was_created;
        /// the adapt mode for refinement
        Geometry::AdaptMode _adapt_mode;
        /// the extern partition sets
        Geometry::PartitionSet _parti_set;

        /// allow extern partitioner?
        bool _allow_parti_extern;
        /// allow 2-level partitioner?
        bool _allow_parti_2level;
        /// allow metis partitioner?
        bool _allow_parti_metis;
        /// allow genetic partitioner?
        bool _allow_parti_genetic;
        /// allow naive partitioner?
        bool _allow_parti_naive;

        /// support double layered hierarchy?
        bool _support_double_layered;
        /// desire double layered hierarchy?
        bool _desired_double_layered;

        /// desired maximum level
        int _desired_level_max;
        /// desired median level (double layered only)
        int _desired_level_med;
        /// desired minimum level
        int _desired_level_min;

        /// required partition name for extern partitioning
        std::deque<String> _extern_parti_names;
        /// required number of elements per rank for a-posteriori partitioning
        int _required_elems_per_rank;
        /// time for genetic partitioner initialisation
        double _genetic_time_init;
        /// time for genetic partitioner mutation
        double _genetic_time_mutate;

        /// information about the chosen partitioning
        String _chosen_parti_info;
        /// the chosen partitioning level
        int _chosen_parti_level;
        /// whether an a-priori partitioning was chosen
        bool _chosen_parti_apriori;
        /// the chosen partitioning graph (aka "elems-at-rank")
        Adjacency::Graph _chosen_parti_graph;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] comm
         * The communicator to be used.
         *
         * \param[in] support_double_layered
         * Specifies whether the controller is allowed to create double-layered hierarchies.
         */
        explicit PartiDomainControl(const Dist::Comm& comm_, bool support_double_layered = false) :
          BaseClass(comm_),
          _was_created(false),
          _adapt_mode(Geometry::AdaptMode::chart),
          _allow_parti_extern(true),
          _allow_parti_2level(true),
          _allow_parti_metis(false),   // this one sucks
          _allow_parti_genetic(false), // this one is exotic
          _allow_parti_naive(true),
          _support_double_layered(support_double_layered),
          _desired_double_layered(false),
          _desired_level_max(-1),
          _desired_level_med(-1),
          _desired_level_min(-1),
          _extern_parti_names(),
          _required_elems_per_rank(1),
          _genetic_time_init(1.0),
          _genetic_time_mutate(1.0),
          _chosen_parti_info(),
          _chosen_parti_level(-1),
          _chosen_parti_apriori(false),
          _chosen_parti_graph()
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
        bool parse_args(SimpleArgParser& args)
        {
          // Try to parse --parti-type <types...>
          {
            auto it = args.query("parti-type");
            if(it != nullptr)
            {
              // disallow all strategies by default
              _allow_parti_extern = _allow_parti_2level = false;
              _allow_parti_metis = _allow_parti_naive = false;
              _allow_parti_genetic = false;

              // loop over all allowed strategies
              for(const auto& t : it->second)
              {
                if     (t.compare_no_case("extern") == 0)
                  _allow_parti_extern = true;
                else if(t.compare_no_case("2level") == 0)
                  _allow_parti_2level = true;
                else if(t.compare_no_case("metis") == 0)
                  _allow_parti_metis = true;
                else if(t.compare_no_case("genetic") == 0)
                  _allow_parti_genetic = true;
                else if(t.compare_no_case("naive") == 0)
                  _allow_parti_naive = true;
                else
                {
                  this->_comm.print("ERROR: unknown partitioner type '" + t + "'");
                  return false;
                }
              }
            }
          }

          // parse --parti-extern-name <names...>
          {
            auto it = args.query("parti-extern-name");
            if(it != nullptr)
              _extern_parti_names = it->second;
          }

          // parse --parti-rank-elems <count>
          args.parse("parti-rank-elems", _required_elems_per_rank);

          // parse --parti-genetic-time <time-init> <time-mutate>
          args.parse("parti-genetic-time", _genetic_time_init, _genetic_time_mutate);

          // okay
          return true;
        }

        /**
         * \brief Parses the partitioner options from a PropertyMap
         *
         * \param[in] pmap
         * Contains the configuration, e.g. read from a file
         *
         * \returns
         * \c true, if the parsing was successful, or \c false,
         * if at least one option was invalid.
         */
        bool parse_property_map(PropertyMap& pmap)
        {
          auto parti_type_p = pmap.query("parti-type");
          if(parti_type_p.second)
          {
            // disallow all strategies by default
            _allow_parti_extern = _allow_parti_2level = false;
            _allow_parti_metis = _allow_parti_naive = false;
            _allow_parti_genetic = false;

            std::deque<String> allowed_partitioners;
            parti_type_p.first.split_by_charset(allowed_partitioners);

            for(const auto& t : allowed_partitioners)
            {
              if     (t == "extern")
                _allow_parti_extern = true;
              else if(t == "2level")
                _allow_parti_2level = true;
              else if(t == "genetic")
                _allow_parti_genetic = true;
              else if(t == "metis")
                _allow_parti_metis = true;
              else if(t == "naive")
                _allow_parti_naive = true;
              else
              {
                this->_comm.print("ERROR: unknown partitioner type '" + t + "'");
                return false;
              }
            }
          }

          auto parti_extern_name_p = pmap.query("parti-extern-name");
          if(parti_extern_name_p.second)
          {
            parti_extern_name_p.first.split_by_charset(_extern_parti_names);
          }

          auto parti_rank_elems_p = pmap.query("parti-rank-elems");
          if(parti_rank_elems_p.second)
          {
            if(!parti_rank_elems_p.first.parse(_required_elems_per_rank))
            {
              this->_comm.print("ERROR: Failed to parse 'parti-rank-elems'");
              return false;
            }
          }

          auto genetic_time_init_p = pmap.query("parti-genetic-time-init");
          if(genetic_time_init_p.second)
          {
            if(!genetic_time_init_p.first.parse(_genetic_time_init))
            {
              this->_comm.print("ERROR: Failed to parse 'parti-genetic-time-init'");
              return false;
            }
          }

          auto genetic_time_mutate_p = pmap.query("parti-genetic-time-mutate");
          if(genetic_time_mutate_p.second)
          {
            if(!genetic_time_mutate_p.first.parse(_genetic_time_mutate))
            {
              this->_comm.print("ERROR: Failed to parse 'parti-genetic-time-mutate'");
              return false;
            }
          }

          return true;
        }

        /**
         * \brief Sets the adapt-mode for refinement
         *
         * \param[in] adapt_mode
         * The adapt-mode that is to be used.
         */
        void set_adapt_mode(Geometry::AdaptMode adapt_mode)
        {
          _adapt_mode = adapt_mode;
        }

        /**
         * \brief Gets the adapt-mode for refinement
         *
         * \returns
         * The adapt-mode that is used.
         */
        Geometry::AdaptMode get_adapt_mode()
        {
          return _adapt_mode;
        }

        /**
         * \brief Sets the desired levels for the partitioned hierarchy.
         *
         * \param[in] slvls
         * A deque of strings containing the desired levels.
         */
        void set_desired_levels(const std::deque<String>& slvls)
        {
          // parse all strings into ints
          std::deque<int> ilvls(slvls.size(), 0);
          for(std::size_t i(0); i < slvls.size(); ++i)
          {
            if(!slvls.at(i).parse(ilvls.at(i)))
              throw InternalError(__func__, __FILE__, __LINE__, "Failed to parse level '" + slvls.at(i) + "'");
          }

          // okay
          set_desired_levels(ilvls);
        }

        /**
         * \brief Sets the desired levels for the partitioned hierarchy.
         *
         * \param[in] ilvls
         * A deque containing the desired levels.
         */
        void set_desired_levels(const std::deque<int>& lvls)
        {
          if(lvls.empty())
            throw InternalError(__func__, __FILE__, __LINE__, "Failed to set empty levels");
          else if(lvls.size() == std::size_t(1))
            set_desired_levels(lvls.front());
          else if(lvls.size() == std::size_t(2))
            set_desired_levels(lvls.front(), lvls.back());
          else if((lvls.size() == std::size_t(3)) && _support_double_layered)
            set_desired_levels(lvls.front(), lvls.at(1u), lvls.back());
          else
            throw InternalError(__func__, __FILE__, __LINE__, "Invalid number of levels specified");
        }

        /**
         * \brief Sets the desired refinement levels for a single-layered hierarchy.
         *
         * \param[in] lvl_max
         * The desired maximum refinement level. Must be >= 0.
         *
         * \param[in] lvl_min
         * The desired minimum refinement level. Must be <= lvl_max.
         * If < 0, then (lvl_max+lvl_min+1) is used as the minimum level.
         */
        void set_desired_levels(int lvl_max, int lvl_min = -1)
        {
          // single-layered hierarchy
          _desired_double_layered = false;
          _desired_level_max = lvl_max;
          _desired_level_min = (lvl_min >= 0 ? lvl_min : lvl_max + lvl_min + 1);

          if(_desired_level_max < 0)
            throw InternalError(__func__, __FILE__, __LINE__, "Invalid level-max");
          if(_desired_level_min < 0)
            throw InternalError(__func__, __FILE__, __LINE__, "Invalid level-min");
          if(_desired_level_max < _desired_level_min)
            throw InternalError(__func__, __FILE__, __LINE__, "Invalid level-min/max combination");
        }

        /**
         * \brief Sets the desired refinement levels for a double-layered hierarchy.
         *
         * \param[in] lvl_max
         * The desired maximum refinement level. Must be >= 0.
         *
         * \param[in] lvl_med
         * The desired median refinement level. Must be < lvl_max
         * If < 0, then (lvl_max+lvl_med+1) is used as the median level.
         *
         * \param[in] lvl_min
         * The desired minimum refinement level. Must be <= lvl_med.
         * If < 0, then (lvl_med+lvl_min+1) is used as the minimum level.
         */
        void set_desired_levels(int lvl_max, int lvl_med, int lvl_min)
        {
          // support double layered hierarchy?
          if(!_support_double_layered)
            throw InternalError(__func__, __FILE__, __LINE__, "double-layered hierarchy not supported");

          // double-layered hierarchy
          _desired_double_layered = true;
          _desired_level_max = lvl_max;
          _desired_level_med = (lvl_med >= 0 ? lvl_med : lvl_max + lvl_med + 1);
          _desired_level_min = (lvl_min >= 0 ? lvl_min : lvl_med + lvl_min + 1);

          if(_desired_level_max < 0)
            throw InternalError(__func__, __FILE__, __LINE__, "Invalid level-max");
          if(_desired_level_med < 0)
            throw InternalError(__func__, __FILE__, __LINE__, "Invalid level-med");
          if(_desired_level_min < 0)
            throw InternalError(__func__, __FILE__, __LINE__, "Invalid level-min");

          // level_max must be strictly greater than level_med
          if(_desired_level_max <= _desired_level_med)
            throw InternalError(__func__, __FILE__, __LINE__, "Invalid level-max/med combination");

          // level_med must be greater or equal level_min
          if(_desired_level_med < _desired_level_min)
            throw InternalError(__func__, __FILE__, __LINE__, "Invalid level-med/min combination");
        }

        /**
         * \brief Returns The desired maximum refinement level
         *
         * \returns The desired maximum refinement level
         */
        int get_desired_level_max() const
        {
          return _desired_level_max;
        }

        /**
         * \brief Returns The desired median refinement level
         *
         * \returns The desired median refinement level
         */
        int get_desired_level_med() const
        {
          return _desired_level_med;
        }

        /**
         * \brief Returns The desired minimum refinement level
         *
         * \returns The desired minimum refinement level
         */
        int get_desired_level_min() const
        {
          return _desired_level_min;
        }

        /**
         * \brief Returns an informative string about the chosen partitioning.
         *
         * \note
         * This function only returns something useful after the #create()
         * function has been called.
         */
        String get_chosen_parti_info() const
        {
          return _chosen_parti_info;
        }

        /**
         * \brief Returns the level on which the base mesh was partitioned.
         */
        int get_chosen_parti_level() const
        {
          return _chosen_parti_level;
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
        void create(const std::deque<String>& filenames, String dirpath = "")
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
         * The mesh reader from which the domain control is to be created.
         */
        void create(Geometry::MeshFileReader& mesh_reader)
        {
          // ensure that the domain control is still empty
          XASSERTM(!_was_created, "domain has already been created");
          XASSERT(this->size_physical() == std::size_t(0));
          XASSERT(this->size_virtual() == std::size_t(0));

          TimeStamp stamp_create;

          // create a new base-mesh node
          std::shared_ptr<MeshNodeType> base_mesh_node = std::make_shared<MeshNodeType>(nullptr, &this->_atlas);

          // try to the read the base-mesh
          mesh_reader.parse(*base_mesh_node, this->_atlas, &_parti_set);

          // create the domain control
#ifdef FEAT_HAVE_MPI
          if(this->_comm.size() == 1)
          {
            // We've got just one process, so it's a simple choice:
            this->_create_single_process(base_mesh_node);
          }
          else if(_support_double_layered && _desired_double_layered)
          {
            // The application supports double-layered domain controls and
            // the user wants that, so create a double-layered one:
            this->_create_double_layered(base_mesh_node);
          }
          else
          {
            // Create a single-layered domain control:
            this->_create_single_layered(base_mesh_node);
          }
#else // not FEAT_HAVE_MPI
          {
            // In the non-MPI case, we always have only one process:
            this->_create_single_process(base_mesh_node);
          }
#endif // FEAT_HAVE_MPI

          // cleanup
          this->_post_create_cleanup();

          // finally, compile virtual levels
          this->compile_virtual_levels();

          // collect partition statistics
          FEAT::Statistics::toe_partition = stamp_create.elapsed_now();

          // domain created
          _was_created = true;
        }

      protected:
        /**
         * \brief Creates a single-layered mesh hierarchy for a single process.
         *
         * \param[in] base_mesh_node
         * The base-mesh node from which the hierarchy is to be derived from.
         */
        void _create_single_process(std::shared_ptr<MeshNodeType> base_mesh_node)
        {
          // create and push single layer
          this->push_layer(std::make_shared<LayerType>(this->_comm.comm_dup(), 0));

          // refine up to desired minimum level
          int lvl = 0;
          for(; lvl < this->_desired_level_min; ++lvl)
          {
            // refine the patch mesh
            base_mesh_node = std::shared_ptr<MeshNodeType>(base_mesh_node->refine(this->_adapt_mode));
          }

          // refine up to maximum level and push to control
          for(; lvl < this->_desired_level_max; ++lvl)
          {
            // push this level
            this->push_level_front(0, std::make_shared<LevelType>(lvl, base_mesh_node));

            // refine the patch mesh
            base_mesh_node = std::shared_ptr<MeshNodeType>(base_mesh_node->refine(this->_adapt_mode));
          }

          // push finest level
          this->push_level_front(0, std::make_shared<LevelType>(lvl, base_mesh_node));
        }

        /**
         * \brief Creates a single-layered mesh hierarchy.
         *
         * \param[in] base_mesh_node
         * The base-mesh node from which the hierarchy is to be derived from.
         */
        void _create_single_layered(std::shared_ptr<MeshNodeType> base_mesh_node)
        {
          // get my rank and the number of procs
          const int rank = this->_comm.rank();
          //const int nprocs = this->_comm.size();

          // create and push single layer
          std::shared_ptr<LayerType> layer = std::make_shared<LayerType>(this->_comm.comm_dup(), 0);
          this->push_layer(layer);

          // choose a partitioning strategy
          this->_check_parti(*base_mesh_node);

          // refine up to chosen partition level
          int lvl = 0;
          for(; lvl < this->_chosen_parti_level; ++lvl)
          {
            // refine the base mesh
            base_mesh_node = std::shared_ptr<MeshNodeType>(base_mesh_node->refine(this->_adapt_mode));
          }

          // apply partitioner
          if(!this->_apply_parti(*base_mesh_node))
          {
            throw InternalError(__func__, __FILE__, __LINE__, "Failed to find a suitable partitioning");
          }

          // extract our patch
          std::vector<int> neighbour_ranks;
          std::shared_ptr<MeshNodeType> patch_mesh_node(
            base_mesh_node->extract_patch(neighbour_ranks, this->_chosen_parti_graph, rank));

          // set our neighbour ranks of our child layer
          layer->set_neighbour_ranks(neighbour_ranks);

          // refine up to minimum level
          for(; lvl < this->_desired_level_min; ++lvl)
          {
            // refine the patch mesh
            patch_mesh_node = std::shared_ptr<MeshNodeType>(patch_mesh_node->refine(this->_adapt_mode));
          }

          // refine up to maximum level
          for(; lvl < this->_desired_level_max; ++lvl)
          {
            // push this level
            this->push_level_front(0, std::make_shared<LevelType>(lvl, patch_mesh_node));

            // refine the patch mesh
            patch_mesh_node = std::shared_ptr<MeshNodeType>(patch_mesh_node->refine(this->_adapt_mode));
          }

          // push finest level
          this->push_level_front(0, std::make_shared<LevelType>(lvl, patch_mesh_node));
        }

        /**
         * \brief Creates a double-layered mesh hierarchy.
         *
         * \param[in] base_mesh_node
         * The base-mesh node from which the hierarchy is to be derived from.
         */
        void _create_double_layered(std::shared_ptr<MeshNodeType> base_mesh_node)
        {
          // get my rank and the number of procs
          const int rank = this->_comm.rank();
          const int nprocs = this->_comm.size();

          // is this the parent process?
          const bool is_parent = (rank == 0);

          // create first layer (child)
          std::shared_ptr<LayerType> layer_c = std::make_shared<LayerType>(this->_comm.comm_dup(), 0);

          // add parent rank to our child layer
          layer_c->push_parent(0);

          // push the first layer
          this->push_layer(layer_c);

          // is this the parent process?
          if(is_parent)
          {
            // create second layer (parent)
            std::shared_ptr<LayerType> layer_p = std::make_shared<LayerType>(Dist::Comm::self(), 1);

            // push all child process ranks to parent layer
            for(int i(0); i < nprocs; ++i)
              layer_p->push_child(i);

            // push the second layer (parent)
            this->push_layer(layer_p);
          }

          // choose a partitioning strategy
          this->_check_parti(*base_mesh_node);

          // refine up to chosen partitioning level
          int lvl = 0;
          for(; lvl < this->_chosen_parti_level; ++lvl)
          {
            // push base-mesh into second layer if desired
            if(is_parent && (lvl >= this->_desired_level_min))
            {
              this->push_level_front(1, std::make_shared<LevelType>(lvl, base_mesh_node));
            }

            // refine the base mesh
            base_mesh_node = std::shared_ptr<MeshNodeType>(base_mesh_node->refine(this->_adapt_mode));
          }

          // apply partitioner
          if(!this->_apply_parti(*base_mesh_node))
          {
            throw InternalError(__func__, __FILE__, __LINE__, "Failed to find a suitable partitioning");
          }

          // extract our patch
          std::vector<int> neighbour_ranks;
          std::shared_ptr<MeshNodeType> patch_mesh_node(
            base_mesh_node->extract_patch(neighbour_ranks, this->_chosen_parti_graph, rank));

          // create child patch mesh-parts on parent process
          if(is_parent)
          {
            // Note: patch mesh-part for rank = 0 was already created by 'extract_patch' call
            for(int i(1); i < nprocs; ++i)
            {
              base_mesh_node->create_patch_meshpart(this->_chosen_parti_graph, i);
            }
          }

          // set our neighbour ranks of our child layer
          layer_c->set_neighbour_ranks(neighbour_ranks);

          // refine up to desired median level (if necessary)
          for(; lvl < this->_desired_level_med; ++lvl)
          {
            // parent process?
            if(is_parent)
            {
              // push base-mesh into second layer if desired
              if(lvl >= this->_desired_level_min)
              {
                this->push_level_front(1, std::make_shared<LevelType>(lvl, base_mesh_node));
              }

              // refine the base mesh
              base_mesh_node = std::shared_ptr<MeshNodeType>(base_mesh_node->refine(this->_adapt_mode));
            }

            // refine the patch mesh
            patch_mesh_node = std::shared_ptr<MeshNodeType>(patch_mesh_node->refine(this->_adapt_mode));
          }

          // push the finest base-mesh
          if(is_parent)
          {
            this->push_level_front(1, std::make_shared<LevelType>(lvl, base_mesh_node));
          }

          // refine up to maximum level
          for(; lvl < this->_desired_level_max; ++lvl)
          {
            // push patch mesh to this level
            this->push_level_front(0, std::make_shared<LevelType>(lvl, patch_mesh_node));

            // refine the patch mesh
            patch_mesh_node = std::shared_ptr<MeshNodeType>(patch_mesh_node->refine(this->_adapt_mode));
          }

          // push finest level
          this->push_level_front(0, std::make_shared<LevelType>(lvl, patch_mesh_node));
        }

        /**
         * \brief Checks for an appropriate partitioning strategy.
         *
         * This function examines the given base mesh and checks whether
         * one of the a-priori partitioning strategies (extern or 2-level)
         * yields a valid partitioning for the given number of processes.
         * If so, the corresponding partitioning level and graph are stored
         * in _chosen_parti_level and _chosen_parti_graph.
         * If not, this function determines how often the base-mesh has to be
         * refined until its has enough cells so that one of the a-posteriori
         * partitioners can be applied.
         *
         * \param[in] base_mesh_node
         * The base-mesh node that is to be partitioned.
         *
         * \returns
         * \c true, if an a-priori partitioning was found, otherwise \c false.
         */
        bool _check_parti(const MeshNodeType& base_mesh_node)
        {
          // Try to find an appropriate a-priori partitioning first:
          if(this->_allow_parti_extern && this->_check_parti_extern())
            return true;
          if(this->_allow_parti_2level && this->_check_parti_2level(base_mesh_node))
            return true;

          // No a-priori partitioning found, so we need to determine the
          // required partitioning level for a-posteriori partitioning.
          // For this, first determine the factor by which the number of elements
          // increases for each refinement:
          typedef typename MeshType::ShapeType ShapeType;
          const Index factor = Index(Geometry::Intern::StandardRefinementTraits<ShapeType, ShapeType::dimension>::count);

          // compute minimum number of elements for a-posteriori partitioning
          const Index min_elems = Index(this->_comm.size()) * Index(_required_elems_per_rank);

          // Okay, get the number of elements on base-mesh level 0
          Index num_elems = base_mesh_node.get_mesh()->get_num_elements();

          // compute refinement level on which we have enough elements
          int level = 0;
          for(; num_elems < min_elems; ++level)
          {
            num_elems *= factor;
          }

          // finally, the user may have chosen a greater level for partitioning:
          this->_chosen_parti_level = Math::max(this->_chosen_parti_level, level);

          // No a-priori partitioning; try a-posteriori partitioner later
          return false;
        }

        /**
         * \brief Applies an a-posteriori partitioner.
         *
         * \param[in] base_mesh_node
         * The base-mesh node that is to be partitioned.
         *
         * \returns
         * \c true, if the base-mesh was partitioned, otherwise \c false.
         */
        bool _apply_parti(/*const*/ MeshNodeType& base_mesh_node)
        {
          // First of all, check whether an a-priori partitioning was selected;
          // if so, then we do not have to apply any a-posteriori partitioner
          if(this->_chosen_parti_apriori)
            return true;

          // ensure that the mesh has enough elements
          XASSERT(this->_comm.size() <= int(base_mesh_node.get_mesh()->get_num_elements()));

          // try the various a-posteriori partitioners
          if(this->_allow_parti_metis && this->_apply_parti_metis(base_mesh_node))
            return true;
          if(this->_allow_parti_genetic && this->_apply_parti_genetic(base_mesh_node))
            return true;
          if(this->_allow_parti_naive && this->_apply_parti_naive(base_mesh_node))
            return true;

          // we should not arrive here...
          return false;
        }

        /**
         * \brief Checks whether an extern partition is given
         *
         * \returns
         * \c true, if an extern partition is given, otherwise \c false.
         */
        bool _check_parti_extern()
        {
          const int num_procs = this->_comm.size();

          // check whether we have a suitable partition
          const Geometry::Partition* part = this->_parti_set.find_partition(num_procs, _extern_parti_names);
          if(part == nullptr)
            return false;

          // found a valid extern partitioning
          this->_chosen_parti_apriori = true;
          this->_chosen_parti_info = String("Found extern partition '") + part->get_name() +
            "' on level " + stringify(part->get_level());
          this->_chosen_parti_level = part->get_level();
          this->_chosen_parti_graph = part->get_patches().clone();
          return true;
        }

        /**
         * \brief Checks whether the 2-level partitioner can be applied.
         *
         * \param[in] base_mesh_node
         * The base-mesh node that is to be partitioned.
         *
         * \returns
         * \c true, if the 2-level partitioner can be applied, otherwise \c false.
         */
        bool _check_parti_2level(const MeshNodeType& base_mesh_node)
        {
          const Index num_procs = Index(this->_comm.size());

          // create a 2-level partitioner
          Geometry::Parti2Lvl<MeshType> partitioner(*base_mesh_node.get_mesh(), num_procs);

          // successful?
          if(!partitioner.success())
            return false;

          // found a valid 2-level partitioning
          this->_chosen_parti_apriori = true;
          this->_chosen_parti_info = String("Found 2-level partition on level ") + stringify(partitioner.parti_level());
          this->_chosen_parti_level = int(partitioner.parti_level());
          this->_chosen_parti_graph = partitioner.build_elems_at_rank();
          return true;
        }

        /**
         * \brief Applies the METIS partitioner onto the base-mesh.
         *
         * \param[in] base_mesh_node
         * The base-mesh node that is to be partitioned.
         *
         * \returns
         * \c true, if METIS was applied successfully, otherwise \c false.
         */
        bool _apply_parti_metis(const MeshNodeType& base_mesh_node)
        {
#ifdef FEAT_HAVE_PARMETIS
          // define partition executor
          typedef Foundation::PExecutorParmetis<Foundation::ParmetisModePartKway> PExeType;

          // get our base mesh
          const auto& base_root_mesh = *base_mesh_node.get_mesh();

          // get number of elements
          const Index num_global_elements(base_root_mesh.get_num_elements());

          // allocate graph
          typename PExeType::PGraphT global_dual(base_root_mesh, num_global_elements, this->_comm);

          // local input for k-way partitioning
          auto local_dual(global_dual.create_local());

          auto part(PExeType::part(*((typename PExeType::PGraphT*)local_dual.get())));

          auto synched_part(Foundation::PSynch<PExeType>::exec(part, typename PExeType::IndexType(num_global_elements)));

          PExeType::fill_comm_structs_global(synched_part, global_dual);

          // render elements-at-rank graph
          this->_chosen_parti_graph = Adjacency::Graph(Adjacency::rt_transpose, synched_part.rank_at_element());

          // verify that each process has at least one element
          {
            const Index* ptr = this->_chosen_parti_graph.get_domain_ptr();
            for(Index i(0); i < this->_chosen_parti_graph.get_num_nodes_domain(); ++i)
            {
              if(ptr[i] == ptr[i+1])
              {
                this->_chosen_parti_info = String("ERROR: Process ") + stringify(i) + " received empty patch from METIS";
                return false;
              }
            }
          }

          // set info string
          this->_chosen_parti_info = String("Applied METIS partitioner on level ") + stringify(this->_chosen_parti_level);

          // okay
          return true;
#else // not FEAT_HAVE_PARMETIS
          (void)base_mesh_node;
          return false;
#endif // FEAT_HAVE_PARMETIS
        }

        /**
         * \brief Applies the genetic partitioner onto the base-mesh.
         *
         * \param[in] base_mesh_node
         * The base-mesh node that is to be partitioned.
         *
         * \returns
         * \c true.
         */
        bool _apply_parti_genetic(/*const*/ MeshNodeType& base_mesh_node)
        {
          // create a genetic partitioner
          Geometry::PartiIterative<MeshType> partitioner(
            *base_mesh_node.get_mesh(),
            this->_comm,
            this->_genetic_time_init,
            this->_genetic_time_mutate);

          // create elems-at-rank graph
          this->_chosen_parti_graph = partitioner.build_elems_at_rank();

          // set info string
          this->_chosen_parti_info = String("Applied genetic partitioner on level ") + stringify(this->_chosen_parti_level);

          // okay
          return true;
        }

        /**
         * \brief Applies the naive partitioner onto the base-mesh.
         *
         * \param[in] base_mesh_node
         * The base-mesh node that is to be partitioned.
         *
         * \returns
         * \c true.
         */
        bool _apply_parti_naive(const MeshNodeType& base_mesh_node)
        {
          const Index num_parts = Index(this->_comm.size());
          const Index num_elems = base_mesh_node.get_mesh()->get_num_elements();
          XASSERTM(num_parts <= num_elems, "Base-Mesh does not have enough elements");

          // create elems-at-rank graph
          this->_chosen_parti_graph = Adjacency::Graph(num_parts, num_elems, num_elems);
          Index* ptr = this->_chosen_parti_graph.get_domain_ptr();
          Index* idx = this->_chosen_parti_graph.get_image_idx();
          ptr[0] = Index(0);
          for(Index i(1); i < num_parts; ++i)
            ptr[i] = (i*num_elems) / num_parts;
          ptr[num_parts] = num_elems;
          for(Index j(0); j < num_elems; ++j)
            idx[j] = j;

          // set info string
          this->_chosen_parti_info = String("Applied naive partitioner on level ") + stringify(this->_chosen_parti_level);

          // okay
          return true;
        }

        /**
         * \brief Cleans up after creation.
         *
         * This function releases some internal data which is not required
         * anymore once the domain control has been created.
         */
        void _post_create_cleanup()
        {
          // clear partition set
          _parti_set.clear();
        }
      }; // class PartiDomainControl<...>
    } // namespace Domain
  } // namespace Control
} // namespace FEAT

#endif // CONTROL_DOMAIN_PARTI_DOMAIN_CONTROL_HPP
