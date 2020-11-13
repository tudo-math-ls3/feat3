// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef CONTROL_DOMAIN_PARTI_DOMAIN_CONTROL_HPP
#define CONTROL_DOMAIN_PARTI_DOMAIN_CONTROL_HPP 1

#include <kernel/base_header.hpp>

#include <kernel/util/dist.hpp>
#include <kernel/util/dist_file_io.hpp>
#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/property_map.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/partition_set.hpp>
#include <kernel/geometry/parti_2lvl.hpp>
#include <kernel/geometry/parti_iterative.hpp>
#include <kernel/geometry/parti_metis.hpp>

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
      inline void add_supported_pdc_args(SimpleArgParser& args)
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
        /// our mesh-part type
        typedef typename LevelType::PartType MeshPartType;

      protected:
        /**
         * \brief Ancestor info class
         */
        class Ancestor
        {
        public:
          /// the index of the layer that this ancestor object belongs to
          int layer;
          /// the index of the parent layer or -1, if this process is not part of the parent layer
          int layer_p;
          /// the number of processors that participate in this layer
          int num_procs;
          /// the number of partitions for each patch of the parent layer
          int num_parts;
          /// the desired minimum and maximum refinement levels for this layer
          int desired_level_max, desired_level_min;

          int progeny_group, progeny_child;
          int progeny_first, progeny_count;
          Dist::Comm progeny_comm;

          /// a string containing some information about the chosen partitioning
          String parti_info;
          /// the refinement level on which the patch is to be partitioned
          int parti_level;
          /// specifies whether the chosen partitioning is an a-priori partitioning strategy
          bool parti_apriori;
          /// this is the actual elements-at-rank partitioning graph
          Adjacency::Graph parti_graph;

          Ancestor() :
            layer(0), layer_p(0), num_procs(0), num_parts(0),
            desired_level_max(0), desired_level_min(0),
            progeny_group(0), progeny_child(0),
            progeny_first(0), progeny_count(0),
            progeny_comm(),
            parti_info(),
            parti_level(0),
            parti_apriori(false),
            parti_graph()
          {
          }
        }; // class Ancestor

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

        /// support multi-layered hierarchy?
        bool _support_multi_layered;

        /// desired level deque
        std::deque<std::pair<int,int>> _desired_levels;

        /// chosen level deque
        std::deque<std::pair<int,int>> _chosen_levels;

        /// required partition name for extern partitioning
        std::deque<String> _extern_parti_names;
        /// required number of elements per rank for a-posteriori partitioning
        int _required_elems_per_rank;
        /// time for genetic partitioner initialization
        double _genetic_time_init;
        /// time for genetic partitioner mutation
        double _genetic_time_mutate;

        /// the partition ancestry deque
        std::deque<Ancestor> _ancestry;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] comm_
         * The main communicator to be used.
         *
         * \param[in] support_multi_layered
         * Specifies whether the controller is allowed to create multi-layered hierarchies.
         */
        explicit PartiDomainControl(const Dist::Comm& comm_, bool support_multi_layered) :
          BaseClass(comm_),
          _was_created(false),
          _adapt_mode(Geometry::AdaptMode::chart),
          _allow_parti_extern(true),
          _allow_parti_2level(true),
          _allow_parti_metis(false),
          _allow_parti_genetic(false), // this one is exotic
          _allow_parti_naive(true),
          _support_multi_layered(support_multi_layered),
          _desired_levels(),
          _extern_parti_names(),
          _required_elems_per_rank(1),
          _genetic_time_init(5),
          _genetic_time_mutate(5),
          _ancestry()
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

            std::deque<String> allowed_partitioners = parti_type_p.first.split_by_whitespaces();

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
            this->_extern_parti_names = parti_extern_name_p.first.split_by_whitespaces();
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
         * A string containing the list of desired levels.
         */
        void set_desired_levels(const String& slvls)
        {
          std::deque<String> sl = slvls.split_by_whitespaces();
          set_desired_levels(sl);
        }

        /**
         * \brief Sets the desired levels for the partitioned hierarchy.
         *
         * \param[in] slvls
         * A deque of strings containing the desired levels.
         *
         * Each entry of the deque is supposed to be in the format <c>level:nprocs</c>.
         */
        void set_desired_levels(const std::deque<String>& slvls)
        {
          const int nranks = this->_comm.size();
          std::deque<String> sv;

          for(std::size_t i(0); i < slvls.size(); ++i)
          {
            int ilvl = -1, nprocs = -1;
            sv = slvls.at(i).split_by_string(":");

            if((sv.size() < std::size_t(1)) || (sv.size() > std::size_t(2)))
              throw ParseError("Invalid input format", slvls.at(i), "an int-pair 'level:patches'");

            if(!sv.front().parse(ilvl))
              throw ParseError("Failed to parse level index", slvls.at(i), "an integer");
            if((sv.size() > std::size_t(1)) && !sv.back().parse(nprocs))
              throw ParseError("Failed to parse process count" , slvls.at(i), "an integer");

            // level must be non-negative
            if(ilvl < 0)
              throw ParseError("Invalid negative level index", slvls.at(i), "a non-negative level index");

            // first level index?
            if(i == std::size_t(0))
            {
              if((nprocs >= 0) && (nprocs != this->_comm.size()))
                throw ParseError("Invalid number of processes for global level: '" + slvls.at(i) +
                  "', expected " + stringify(this->_comm.size()) + " but got " + stringify(nprocs));
              _desired_levels.push_back(std::make_pair(ilvl, nranks));
              continue;
            }

            // make sure the level is non-ascending
            if(_desired_levels.back().first < ilvl)
              throw ParseError("Invalid non-descending level index: '" + slvls.at(i) +
                "', expected <= " + stringify(_desired_levels.back().first) + " but got " + stringify(ilvl));

            // make sure process count is valid
            if((i + 1) == slvls.size())
              nprocs = 0;
            else
            {
              // the process count must be descending
              if(_desired_levels.back().second <= nprocs)
                throw ParseError("Invalid non-descending process count: '" + slvls.at(i) +
                  "', expected < " + stringify(_desired_levels.back().second) + " but got " + stringify(nprocs));

              // the previous process count must be a multiple
              if(_desired_levels.back().second % nprocs != 0)
                throw ParseError("Invalid indivisible process count: '" + slvls.at(i) +
                  "', expected a divisor of " + stringify(_desired_levels.back().second) + " but got " + stringify(nprocs));
            }

            // push the level-proc pair
            _desired_levels.push_back(std::make_pair(ilvl, nprocs));
          }
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
          int _desired_level_max = lvl_max;
          int _desired_level_min = (lvl_min >= 0 ? lvl_min : lvl_max + lvl_min + 1);

          XASSERTM(_desired_level_max >= 0, "Invalid level-max");
          XASSERTM(_desired_level_min >= 0, "Invalid level-min");
          XASSERTM(_desired_level_max >= _desired_level_min, "Invalid level-min/max combination");

          _desired_levels.emplace_back(std::make_pair(_desired_level_max, this->_comm.size()));
          _desired_levels.emplace_back(std::make_pair(_desired_level_min, 0));
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
          // double-layered hierarchy
          int _desired_level_max = lvl_max;
          int _desired_level_med = (lvl_med >= 0 ? lvl_med : lvl_max + lvl_med + 1);
          int _desired_level_min = (lvl_min >= 0 ? lvl_min : lvl_med + lvl_min + 1);

          XASSERTM(_desired_level_max >= 0, "Invalid level-max");
          XASSERTM(_desired_level_med >= 0, "Invalid level-med");
          XASSERTM(_desired_level_min >= 0, "Invalid level-min");

          // level_max must be strictly greater than level_med
          XASSERTM(_desired_level_max > _desired_level_med, "Invalid level-max/med combination");

          // level_med must be greater or equal level_min
          XASSERTM(_desired_level_med >= _desired_level_min, "Invalid level-med/min combination");

          _desired_levels.emplace_back(std::make_pair(_desired_level_max, this->_comm.size()));
          _desired_levels.emplace_back(std::make_pair(_desired_level_med, 1));
          _desired_levels.emplace_back(std::make_pair(_desired_level_min, 0));
        }

        /**
         * \brief Returns the desired maximum refinement level
         *
         * \returns The desired maximum refinement level
         */
        int get_desired_level_max() const
        {
          return _desired_levels.front().first;
        }

        /**
         * \brief Returns the desired minimum refinement level
         *
         * \returns The desired minimum refinement level
         */
        int get_desired_level_min() const
        {
          return _desired_levels.back().first;
        }

        /**
         * \brief Returns the desired levels formatted as a parsable string.
         *
         * The string returned by this function can be parsed by the set_desired_levels() function.
         *
         * \returns The desired levels formatted as a parsable string.
         */
        String format_desired_levels() const
        {
          String s;
          for(std::size_t i(0); (i + 1) < _desired_levels.size(); ++i)
          {
            s += stringify(_desired_levels.at(i).first);
            s += ":";
            s += stringify(_desired_levels.at(i).second);
            s += "  ";
          }
          s += stringify(_desired_levels.back().first);
          return s;
        }

        /**
         * \brief Returns the chosen levels formatted as a parsable string.
         *
         * The string returned by this function can be parsed by the set_desired_levels() function.
         *
         * \returns The chosen levels formatted as a parsable string.
         */
        String format_chosen_levels() const
        {
          String s;
          for(std::size_t i(0); (i + 1) < _chosen_levels.size(); ++i)
          {
            s += stringify(_chosen_levels.at(i).first);
            s += ":";
            s += stringify(_chosen_levels.at(i).second);
            s += "  ";
          }
          s += stringify(_chosen_levels.back().first);
          return s;
        }

        /**
         * \brief Returns an informative string about the chosen partitioning.
         *
         * \note
         * This function only returns something useful after the #create()
         * function has been called.
         *
         * \note
         * In the case of a multi-layered recursive partitioning, the returned string is
         * a multi-line string containing the partitioning info for each recursive layer
         * in a separate line.
         */
        String get_chosen_parti_info() const
        {
          String s;

          for(auto it = this->_ancestry.rbegin(); it != this->_ancestry.rend(); ++it)
          {
            if(it != this->_ancestry.rbegin())
              s += "\n";
            s += it->parti_info;
            s += " on level ";
            s += stringify(it->parti_level);
            s += " for ";
            s += stringify(it->num_parts);
            s += (it->num_parts > 1 ? " patches on " : " patch on ");
            s += stringify(it->num_procs);
            s += (it->num_procs > 1 ? " processes" : " process");
          }

          return s;
        }

        /**
         * \brief Returns the deque of desired refinement levels.
         *
         * \returns The deque of desired refinement levels.
         */
        const std::deque<std::pair<int, int>> get_desired_levels() const
        {
          return this->_desired_levels;
        }

        /**
         * \brief Returns the deque of chosen refinement levels.
         *
         * \returns The deque of chosen refinement levels.
         */
        const std::deque<std::pair<int, int>> get_chosen_levels() const
        {
          return this->_chosen_levels;
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
          else if(_support_multi_layered && (_desired_levels.size() > std::size_t(2)))
          {
            // The application supports multi-layered domain controls and
            // the user wants that, so create a multi-layered one:
            this->_create_multi_layered(base_mesh_node);
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
          _parti_set.clear();

          // finally, compile virtual levels
          this->compile_virtual_levels();

          // collect partition statistics
          FEAT::Statistics::toe_partition = stamp_create.elapsed_now();

          // domain created
          _was_created = true;
        }

        /**
         * \brief Debugging function: Returns a string containing encoded ancestry information
         */
        String dump_ancestry() const
        {
          String s;
          for(auto it = this->_ancestry.begin(); it != this->_ancestry.end(); ++it)
          {
            if(it != this->_ancestry.begin())
              s += " | ";
            s += stringify((*it).num_procs);
            s += ":";
            s += stringify((*it).num_parts);
            s += "[";
            s += stringify((*it).progeny_group).pad_front(2);
            s += "+";
            s += stringify((*it).progeny_child);
            s += ":";
            s += stringify((*it).progeny_first).pad_front(2);
            s += "+";
            s += stringify((*it).progeny_count);
            s += "]";
            s += "(";
            s += stringify((*it).desired_level_max);
            s += ":";
            s += stringify((*it).desired_level_min);
            s += ")";
            s += ((*it).layer >= 0 ? "*" : " ");
            s += ((*it).layer_p >= 0 ? ">" : " ");
          }
          return s;
        }

        /**
         * \brief Serializes the partitioning information into a binary buffer.
         */
        std::vector<char> serialize_partitioning() const
        {
          BinaryStream bs;
          typedef std::uint64_t u64;

          // serialization magic number: "F3PaDoCo"
          const std::uint64_t magic = 0x6F436F4461503346;

          XASSERT(this->_num_global_layers >= this->_ancestry.size());

          if(this->_comm.rank() == 0)
          {
            // dump magic
            bs.write((const char*)&magic, 8u);

            // write finest level
            const u64 lev = u64(this->_virt_levels.front()->get_level_index());
            bs.write((const char*)&lev, 8u);

            // serialize ancestry/layers
            const u64 na = u64(this->_ancestry.size());
            bs.write((const char*)&na, 8u);

            // write number of ranks per layer in reverse order
            for(auto it = this->_ancestry.rbegin(); it != this->_ancestry.rend(); ++it)
            {
              // write ranks
              const u64 np = u64(it->num_procs);
              const u64 pl = u64(it->parti_level);
              bs.write((const char*)&np, 8u);
              bs.write((const char*)&pl, 8u);
            }
          }

          // loop over all ancestry graphs
          for(auto it = this->_ancestry.rbegin(); it != this->_ancestry.rend(); ++it)
          {
            if(it == this->_ancestry.rbegin())
            {
              std::vector<char> buf = it->parti_graph.serialize();
              bs.write(buf.data(), std::streamsize(buf.size()));
              continue;
            }

            // get layer index
            if(it->layer_p < 0)
              continue;

            // get parent layer communicator
            const std::size_t ilp = std::size_t(it->layer_p);
            XASSERT(ilp < this->_layers.size());
            const Dist::Comm& comm_p = this->_layers.at(ilp)->comm();

            // serialize graph
            std::vector<char> buf = it->parti_graph.serialize();

            // choose maximum size
            u64 buf_size = buf.size();

            // gather individual buffer sizes on rank 0
            std::vector<u64> all_sizes;
            if(comm_p.rank() == 0)
              all_sizes.resize(std::size_t(comm_p.size()));
            comm_p.gather(&buf_size, 1u, all_sizes.data(), 1u, 0);

            // allreduce maximum size
            comm_p.allreduce(&buf_size, &buf_size, std::size_t(1), Dist::op_max);

            // adjust buffer size
            if(buf_size > u64(buf.size()))
              buf.resize(buf_size);

            // on rank 0, allocate common buffer
            std::vector<char> com_buf;
            if(comm_p.rank() == 0)
              com_buf.resize(std::size_t(buf_size * u64(comm_p.size())));

            // gather all buffers on rank 0
            comm_p.gather(buf.data(), buf_size, com_buf.data(), buf_size, 0);

            // write each individual buffer
            if(comm_p.rank() == 0)
            {
              char* x = com_buf.data();
              for(u64 k(0); k < u64(comm_p.size()); ++k)
                bs.write(&x[k*buf_size], std::streamsize(all_sizes[k]));
            }
          }

          return bs.container();
        }

      protected:
        /**
         * \brief Creates the ancestry for a single layer (or a single process)
         */
        void _create_ancestry_single()
        {
          // for more than 2 desired levels, call _create_ancestry_scattered
          XASSERTM(this->_desired_levels.size() <= std::size_t(2), "multi-layered control is desired here");

          // allocate and create ancestry
          this->_ancestry.resize(std::size_t(1));
          Ancestor& ancestor = this->_ancestry.front();

          // set the layer index to 0
          ancestor.layer = 0;

          // set the parent layer index to -1
          ancestor.layer_p = -1;

          // set the total number of processes for this layer
          ancestor.num_procs = this->_comm.size();

          // set the total number of partitions per progeny group
          ancestor.num_parts = this->_comm.size();

          // set desired maximum level
          ancestor.desired_level_max = this->_desired_levels.front().first;

          // set desired minimum level (may be = level_max)
          ancestor.desired_level_min = this->_desired_levels.back().first;

          // set the progeny group
          ancestor.progeny_group = 0;
          ancestor.progeny_child = this->_comm.rank();

          // create the progeny communicator
          ancestor.progeny_count = this->_comm.size();
          ancestor.progeny_first = 0;
          ancestor.progeny_comm = this->_comm.comm_dup();
        }

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
            base_mesh_node = std::shared_ptr<MeshNodeType>(base_mesh_node->refine(this->_adapt_mode));
          }

          // save chosen minimum level if it is not equal to the desired maximum level
          if(lvl < ancestor.desired_level_max)
            this->_chosen_levels.push_front(std::make_pair(lvl, 0));

          // refine up to maximum level and push to control
          for(; lvl < ancestor.desired_level_max; ++lvl)
          {
            // push this level
            this->push_level_front(0, std::make_shared<LevelType>(lvl, base_mesh_node));

            // refine the patch mesh
            base_mesh_node = std::shared_ptr<MeshNodeType>(base_mesh_node->refine(this->_adapt_mode));
          }

          // save chosen maximum level
          this->_chosen_levels.push_front(std::make_pair(lvl, 1));

          // push finest level
          this->push_level_front(0, std::make_shared<LevelType>(lvl, base_mesh_node));
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
        void _create_single_layered(std::shared_ptr<MeshNodeType> base_mesh_node)
        {
          // create and push single layer
          std::shared_ptr<LayerType> layer = std::make_shared<LayerType>(this->_comm.comm_dup(), 0);
          this->push_layer(layer);

          // create single-layered ancestry
          this->_create_ancestry_single();
          Ancestor& ancestor = this->_ancestry.front();

          // choose a partitioning strategy
          this->_check_parti(ancestor, *base_mesh_node, true);

          // refine up to chosen partition level
          int lvl = 0;
          for(; lvl < ancestor.parti_level; ++lvl)
          {
            // refine the base mesh
            base_mesh_node = std::shared_ptr<MeshNodeType>(base_mesh_node->refine(this->_adapt_mode));
          }

          // apply partitioner
          if(!this->_apply_parti(ancestor, *base_mesh_node))
          {
            XABORTM("Failed to find a suitable partitioning");
          }

          // extract our patch
          std::vector<int> neighbor_ranks;
          std::shared_ptr<MeshNodeType> patch_mesh_node(
            base_mesh_node->extract_patch(neighbor_ranks, ancestor.parti_graph, this->_comm.rank()));

          // set the neighbor ranks of our child layer
          layer->set_neighbor_ranks(neighbor_ranks);

          // refine up to minimum level
          for(; lvl < ancestor.desired_level_min; ++lvl)
          {
            // refine the patch mesh
            patch_mesh_node = std::shared_ptr<MeshNodeType>(patch_mesh_node->refine(this->_adapt_mode));
          }

          // save chosen minimum level
          if(lvl < ancestor.desired_level_max)
            this->_chosen_levels.push_front(std::make_pair(lvl, 0));

          // refine up to maximum level
          for(; lvl < ancestor.desired_level_max; ++lvl)
          {
            // push this level
            this->push_level_front(0, std::make_shared<LevelType>(lvl, patch_mesh_node));

            // refine the patch mesh
            patch_mesh_node = std::shared_ptr<MeshNodeType>(patch_mesh_node->refine(this->_adapt_mode));
          }

          // save chosen maximum level
          this->_chosen_levels.push_front(std::make_pair(lvl, this->_comm.size()));

          // push finest level
          this->push_level_front(0, std::make_shared<LevelType>(lvl, patch_mesh_node));
        }

        /**
         * \brief Creates the layers for a multi-layered domain control in a scattered fashion.
         */
        void _create_multi_layers_scattered()
        {
          // create and push global layer
          this->push_layer(std::make_shared<LayerType>(this->_comm.comm_dup(), 0));

          // loop over all desired layers
          for(std::size_t ilay(1u); (ilay+1u) < this->_desired_levels.size(); ++ilay)
          {
            // get child layer
            std::shared_ptr<LayerType> layer_c = this->_layers.back();

            // get child layer communicator
            const Dist::Comm& comm_c = layer_c->comm();

            // get number of processes in child comm
            const int nprocs_c = comm_c.size();

            // get desired number of processes for parent comm
            const int nprocs_p = this->_desired_levels.at(ilay).second;
            XASSERT(nprocs_p < nprocs_c);
            XASSERT(nprocs_p > 0);
            XASSERT(nprocs_c % nprocs_p == 0); // same number of children for each parent

            // compute number of siblings = children per parent
            const int num_sibs = nprocs_c / nprocs_p;

            // get my rank in child comm
            const int my_rank_c = comm_c.rank();

            // compute my parent's rank in child comm
            const int parent_rank_c = (my_rank_c / num_sibs) * num_sibs;

            // create our sibling communicator and set the parent rank
            layer_c->set_parent(comm_c.comm_create_range_incl(num_sibs, parent_rank_c), 0);

            // next, create the actual parent communicator
            Dist::Comm comm_p = comm_c.comm_create_range_incl(nprocs_p, 0, num_sibs);

            // Are we the parent?
            if(parent_rank_c == my_rank_c)
            {
              // make sure we have a valid communicator
              XASSERT(!comm_p.is_null());

              // push the parent layer
              this->push_layer(std::make_shared<LayerType>(std::move(comm_p), int(ilay)));
            }
            else
            {
              // We are not a parent, so we must have received a null communicator
              XASSERT(comm_p.is_null());

              // Exit the loop, as we are not part of the party anymore...
              break;
            }
          }
        }

        /**
         * \brief Creates the layers for a multi-layered domain control in a scattered fashion.
         */
        void _create_ancestry_scattered()
        {
          // the layers must have been created already
          XASSERT(!this->_layers.empty());

          // create a deque with the number of processes for each layer
          std::deque<int> num_procs;
          for(auto it = this->_desired_levels.begin(); it != this->_desired_levels.end(); ++it)
          {
            if(it->second > 1)
              num_procs.push_back(it->second);
          }
          // manually add the base-layer
          num_procs.push_back(1);

          // allocate and create ancestry
          this->_ancestry.resize(num_procs.size()-std::size_t(1));

          // set up the ancestry
          const int main_rank = this->_comm.rank();
          const int main_size = this->_comm.size();
          for(std::size_t i(0); i < this->_ancestry.size(); ++i)
          {
            // get our ancestor info
            Ancestor& ancestor = this->_ancestry.at(i);

            // set the layer index (or -1, if this process is not part of that layer)
            ancestor.layer = (i < this->_layers.size() ? int(i) : -1);

            // set the parent layer index (or -1, if this process is not in the parent layer)
            ancestor.layer_p = ((i+1u) < this->_layers.size() ? int(i)+1 : -1);

            // set the total number of processes for this layer
            ancestor.num_procs = num_procs.at(i);

            // set the total number of partitions per progeny group
            ancestor.num_parts = num_procs.at(i) / num_procs.at(i+1u);

            // set desired maximum level
            XASSERT(i < this->_desired_levels.size());
            ancestor.desired_level_max = this->_desired_levels.at(i).first;

            // set desired minimum level
            if((i+1u) < this->_desired_levels.size())
              ancestor.desired_level_min = this->_desired_levels.at(i+1).first;
            else
              ancestor.desired_level_min = ancestor.desired_level_max;

            // set the progeny group
            ancestor.progeny_group = ((main_rank * num_procs.at(i+1u)) / main_size) * ancestor.num_parts;
            ancestor.progeny_child = ((main_rank * num_procs.at(i)) / main_size) % ancestor.num_parts;

            // create the progeny communicator
            ancestor.progeny_count = main_size / num_procs.at(i+1u);
            ancestor.progeny_first = (main_rank / ancestor.progeny_count) * ancestor.progeny_count;
            ancestor.progeny_comm = this->_comm.comm_create_range_incl(ancestor.progeny_count, ancestor.progeny_first);
          }
        }

        /**
         * \brief Creates a multi-layered mesh hierarchy.
         *
         * \param[in] base_mesh_node
         * The base-mesh node from which the hierarchy is to be derived from.
         */
        void _create_multi_layered(std::shared_ptr<MeshNodeType> base_mesh_node)
        {
          // create layers
          this->_create_multi_layers_scattered();

          // create ancestry
          this->_create_ancestry_scattered();

          // we start counting at level 0
          int lvl = 0;

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
            else if(_ancestry.size()+1u < _desired_levels.size())
              parent_min_lvl = this->_desired_levels.back().first;

            // check available partitioning strategies
            this->_check_parti(ancestor, *base_mesh_node, is_base_layer);

            // the check_parti function returns the partitioning level w.r.t. the current
            // level (which may be > 0), so we have to compensate that by adding our current level:
            ancestor.parti_level += lvl;

            // Note: each progeny group within the main communicator may have chosen a different
            // partitioning level at this point. We will compensate this by adjusting the minimum
            // refinement level of the child layer after the partitioning step below.

            // refine up to the chosen partitioning level
            for(; lvl < ancestor.parti_level; ++lvl)
            {
              // push the base mesh into our parent layer if desired
              if((ancestor.layer_p >= 0) && (parent_min_lvl >= 0) && (lvl >= parent_min_lvl))
              {
                this->push_level_front(ancestor.layer_p, std::make_shared<LevelType>(lvl, base_mesh_node));
              }

              // refine the base-mesh node
              base_mesh_node = std::shared_ptr<MeshNodeType>(base_mesh_node->refine(this->_adapt_mode));
            }

            // we're now at the partitioning level, so apply the partitioner
            if(!this->_apply_parti(ancestor, *base_mesh_node))
            {
              XABORTM("Failed to find a suitable partitioning");
            }

            // extract our patch
            std::vector<int> neighbor_ranks;
            std::shared_ptr<MeshNodeType> patch_mesh_node(
              base_mesh_node->extract_patch(neighbor_ranks, ancestor.parti_graph, ancestor.progeny_child));

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

            // does this process participate in the parent layer?
            if(ancestor.layer_p >= 0)
            {
              // Note: patch mesh-part for rank = 0 was already created by 'extract_patch' call
              for(int i(1); i < ancestor.num_parts; ++i)
              {
                base_mesh_node->create_patch_meshpart(ancestor.parti_graph, i);
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
              std::shared_ptr<MeshNodeType> coarse_mesh_node = base_mesh_node;
              base_mesh_node = std::shared_ptr<MeshNodeType>(coarse_mesh_node->refine(this->_adapt_mode));

              // push base mesh to parent layer if desired
              if((ancestor.layer_p >= 0) && (parent_min_lvl >= 0) && (lvl >= parent_min_lvl))
              {
                // clear patches before pushing this node as they are redundant here
                coarse_mesh_node->clear_patches();
                this->push_level_front(ancestor.layer_p, std::make_shared<LevelType>(lvl, coarse_mesh_node));
              }

              // refine the patch mesh
              patch_mesh_node = std::shared_ptr<MeshNodeType>(patch_mesh_node->refine(this->_adapt_mode));
            }

            // split the halos of our base-mesh and compute the halos of our patches from that
            this->_split_basemesh_halos(ancestor, *base_mesh_node, *patch_mesh_node, neighbor_ranks);

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

            // push the finest base-mesh
            if(ancestor.layer_p >= 0)
            {
              this->push_level_front(ancestor.layer_p, std::make_shared<LevelType>(lvl, base_mesh_node));
            }

            // continue with the next layer
            base_mesh_node = patch_mesh_node;
          }

          // get the desired maximum level
          // if we have more than one layer, make sure that the finest one contains at
          // least one level, as otherwise the finest global level would be a ghost level
          int desired_level_max = Math::max(this->_ancestry.front().desired_level_max, lvl+1);

          for(; lvl < desired_level_max; ++lvl)
          {
            // push patch mesh to this level
            this->push_level_front(0, std::make_shared<LevelType>(lvl, base_mesh_node));

            // refine the patch mesh
            base_mesh_node = std::shared_ptr<MeshNodeType>(base_mesh_node->refine(this->_adapt_mode));
          }

          // push finest level
          this->push_level_front(0, std::make_shared<LevelType>(lvl, base_mesh_node));

          // set chosen maximum level for finest layer
          this->_chosen_levels.push_front(std::make_pair(lvl, this->_comm.size()));
        }

        /**
         * \brief Splits the base-mesh halos and computes the inter-patch-mesh halos.
         *
         * This is where the magic happens in the case of the multi-layered hierarchy.
         *
         * \param[in] ancestor
         * The ancestor object for the current layer.
         *
         * \param[in] base_mesh_node
         * The base-mesh node that has been partitioned and whose halos are to be split.
         *
         * \param[inout] patch_mesh_node
         * The patch-mesh node that represents the partition whose halos are to be computed.
         *
         * \param[inout] neighbor_ranks
         * The vector that receives the ranks of the new neighbors that derive from halo splitting.
         */
        void _split_basemesh_halos(
          const Ancestor& ancestor,
          const MeshNodeType& base_mesh_node,
          MeshNodeType& patch_mesh_node,
          std::vector<int>& neighbor_ranks)
        {
          // get the map of the base-mesh halos
          const std::map<int, MeshPartType*>& base_halo_map = base_mesh_node.get_halo_map();

          // if the base mesh has no halos, then we can jump out of here
          if(base_halo_map.empty())
            return;

          // get number of halos
          const std::size_t num_halos = base_halo_map.size();

          // create a halo splitter
          Geometry::PatchHaloSplitter<MeshType> halo_splitter(*base_mesh_node.get_mesh(),
            *base_mesh_node.get_patch(ancestor.progeny_child));

          // add each base-mesh halo to our halo splitter and store the resulting split data size
          std::vector<int> halo_ranks;
          std::vector<std::size_t> halo_sizes;
          for(auto it = base_halo_map.begin(); it != base_halo_map.end(); ++it)
          {
            // store halo rank
            halo_ranks.push_back(it->first);

            // try to split the halo and store the resulting data size
            halo_sizes.push_back(halo_splitter.add_halo(it->first, *it->second));
          }

          // This vector will receive the split halo data from all our potential neighbor processes
          std::vector<Index> halo_recv_data;
          std::vector<std::vector<Index>> halo_send_data;

          /* ******************************************************************************************************* */
          /* ******************************************************************************************************* */
          // PHASE I: collect halo data from our siblings
          /* ******************************************************************************************************* */
          /* ******************************************************************************************************* */
          if(ancestor.layer >= 0)
          {
            XASSERT(this->_layers.size() > std::size_t(ancestor.layer));
            const LayerType& layer = *this->_layers.at(std::size_t(ancestor.layer));
            const Dist::Comm& sibling_comm = *layer.sibling_comm_ptr();
            XASSERT(!sibling_comm.is_null());
            const std::size_t num_sibls = std::size_t(sibling_comm.size());

            // get the rank of this process in various communicators
            const int sibl_rank = sibling_comm.rank();
            const int layer_rank = layer.comm().rank();

            // gather halo infos of all siblings at sibling rank 0
            std::vector<std::size_t> sibl_halo_sizes(sibl_rank == 0 ? num_halos*num_sibls : 0u);
            sibling_comm.gather(halo_sizes.data(), halo_sizes.size(), sibl_halo_sizes.data(), halo_sizes.size(), 0);

            if(sibl_rank > 0)
            {
              std::vector<std::vector<Index>> halo_split_data(num_halos);
              Dist::RequestVector send_reqs(num_halos);

              // serialize all split halos
              for(std::size_t i(0); i < num_halos; ++i)
              {
                // skip empty halos
                if(halo_sizes.at(i) == Index(0))
                  continue;

                // serialize split halo data
                halo_split_data.at(i) = halo_splitter.serialize_split_halo(halo_ranks[i], layer_rank);
                XASSERT(halo_split_data.at(i).size() == halo_sizes.at(i));

                // send split data over to our parent process
                send_reqs[i] = sibling_comm.isend(halo_split_data.at(i).data(), halo_sizes.at(i), 0);
              }

              // wait for all pending sends to finish
              send_reqs.wait_all();
            }
            else // if(sibl_rank == 0)
            {
              halo_send_data.resize(num_halos);

              // compute halo send data sizes
              for(std::size_t i(0); i < num_halos; ++i)
              {
                // determine the number of child processes for this halo
                // as well as the required send buffer size
                Index num_halo_childs = 0u;
                std::size_t buffer_size = 0u, offset = 0u;
                for(std::size_t j(0); j < num_sibls; ++j)
                {
                  // update required buffer size
                  buffer_size += sibl_halo_sizes.at(j*num_halos + i);

                  // this sibling is a child if its split halo is not empty
                  if(sibl_halo_sizes.at(j*num_halos + i) > std::size_t(0))
                    ++num_halo_childs;
                }

                // increase buffer size to store child data offsets
                buffer_size += num_halo_childs + Index(1);

                // allocate send buffer and get a pointer to its data array
                halo_send_data.at(i).resize(buffer_size);
                Index* halo_buffer = halo_send_data.at(i).data();

                // store child count as first entry of buffer
                halo_buffer[0u] = num_halo_childs;

                // initialize offset for first child
                offset = num_halo_childs + Index(1);
                Index coi = 0u; // child offset index

                // collect my own split halo
                if(sibl_halo_sizes.at(i) > Index(0))
                {
                  std::vector<Index> my_data(halo_splitter.serialize_split_halo(halo_ranks[i], layer_rank));
                  std::size_t data_size = sibl_halo_sizes.at(i);
                  halo_buffer[++coi] = Index(offset);
                  for(std::size_t k(0); k < data_size; ++k)
                    halo_buffer[offset+k] = my_data[k];
                  offset += my_data.size();
                }

                Dist::RequestVector sibl_recv_reqs(num_sibls);

                // collect the other siblings
                for(std::size_t j(1); j < num_sibls; ++j)
                {
                  // skips all sibling with empty halos
                  std::size_t data_size = sibl_halo_sizes.at(j*num_halos + i);
                  if(data_size == std::size_t(0))
                    continue;
                  XASSERT(offset+data_size <= buffer_size);

                  // store offset for this child
                  halo_buffer[++coi] = Index(offset);

                  // receive serialized data from this sibling
                  sibl_recv_reqs[j] = sibling_comm.irecv(&halo_buffer[offset], data_size, int(j));
                  offset += data_size;
                }
                XASSERT(offset == buffer_size);

                // wait for all pending receives to finish
                sibl_recv_reqs.wait_all();
              }
            }
          } // if(ancestor.layer >= 0)

          /* ******************************************************************************************************* */
          /* ******************************************************************************************************* */
          // PHASE II: exchange split halos over parent layer communicator
          /* ******************************************************************************************************* */
          /* ******************************************************************************************************* */

          if(ancestor.layer_p >= 0)
          {
            XASSERT(!halo_send_data.empty());

            // get the parent layer communicator
            const Dist::Comm& parent_comm = this->_layers.at(std::size_t(ancestor.layer_p))->comm();

            std::vector<std::size_t> halo_recv_sizes(num_halos), halo_send_sizes(num_halos);
            Dist::RequestVector halo_recv_reqs(num_halos), halo_send_reqs(num_halos);

            // exchange halo send data sizes
            for(std::size_t i(0); i < num_halos; ++i)
            {
              // get halo send size
              halo_send_sizes.at(i) = halo_send_data.at(i).size();

              // post receive requests
              halo_recv_reqs[i] = parent_comm.irecv(&halo_recv_sizes[i], std::size_t(1), halo_ranks[i]);

              // post send requests
              halo_send_reqs[i] = parent_comm.isend(&halo_send_sizes[i], std::size_t(1), halo_ranks[i]);
            }

            // wait for sends and receives to finish
            halo_recv_reqs.wait_all();
            halo_send_reqs.wait_all();

            // compute total receive buffer size
            std::size_t recv_buf_size = num_halos + std::size_t(1);
            for(std::size_t i(0); i < num_halos; ++i)
              recv_buf_size += halo_recv_sizes[i];

            // allocate receive buffer
            halo_recv_data.resize(recv_buf_size);

            // set up receive data pointers
            halo_recv_data[0] = Index(num_halos) + Index(1);
            for(std::size_t i(0); i < num_halos; ++i)
              halo_recv_data[i+1u] = Index(halo_recv_data[i] + halo_recv_sizes[i]);

            // allocate receive buffers and post receives
            for(std::size_t i(0); i < num_halos; ++i)
            {
              // resize buffer and post receive
              halo_recv_reqs[i] = parent_comm.irecv(&halo_recv_data[halo_recv_data[i]], halo_recv_sizes.at(i), halo_ranks[i]);

              // post send of actual halo buffer
              halo_send_reqs[i] = parent_comm.isend(halo_send_data.at(i).data(), halo_send_sizes.at(i), halo_ranks[i]);
            }

            // wait for sends and receives to finish
            halo_recv_reqs.wait_all();
            halo_send_reqs.wait_all();
          } // if(ancestor.layer_p >= 0)

          /* ******************************************************************************************************* */
          /* ******************************************************************************************************* */
          // PHASE III: broadcast halo receive data over progeny comm
          /* ******************************************************************************************************* */
          /* ******************************************************************************************************* */

          {
            // broadcast receive buffer size
            std::size_t recv_data_size = halo_recv_data.size();
            ancestor.progeny_comm.bcast(&recv_data_size, std::size_t(1), 0);

            // allocate buffer
            if(ancestor.progeny_comm.rank() != 0)
            {
              XASSERT(halo_recv_data.empty()); // should be empty until now
              halo_recv_data.resize(recv_data_size);
            }
            else
            {
              XASSERT(!halo_recv_data.empty()); // must not be empty
            }

            // broadcast buffer
            ancestor.progeny_comm.bcast(halo_recv_data.data(), recv_data_size, 0);
          }

          /*{
            String s;
            for(std::size_t i(0); i < num_halos; ++i)
            {
              s += stringify(halo_ranks[i]);
              s += " >";
              for(Index j(halo_recv_data[i]); j < halo_recv_data[i+1]; ++j)
                (s += " ") += stringify(halo_recv_data[j]);
              s += "\n";
            }
            this->_comm.allprint(s);
          }*/

          /* ******************************************************************************************************* */
          /* ******************************************************************************************************* */
          // PHASE IV: intersect received halo splits
          /* ******************************************************************************************************* */
          /* ******************************************************************************************************* */

          for(std::size_t i(0); i < num_halos; ++i)
          {
            // get the offset of the first date for this halo
            const Index offset = halo_recv_data.at(i);

            // get the number of child processes for this halo
            const Index num_childs = halo_recv_data.at(offset);

            // loop over all child processes
            for(Index j(0); j < num_childs; ++j)
            {
              // compute the offset of this child's data within the large buffer
              const Index buffer_offset = offset + halo_recv_data.at(offset+j+1u);

              // intersect with our other halo
              if(!halo_splitter.intersect_split_halo(halo_ranks[i], halo_recv_data, buffer_offset))
                continue; // no intersection

              // get the new neighbor's rank
              const int neighbor_rank = int(halo_recv_data.at(buffer_offset));

              // add the new neighbor to our list
              neighbor_ranks.push_back(neighbor_rank);

              // create new mesh-part
              patch_mesh_node.add_halo(neighbor_rank, new MeshPartType(halo_splitter));
            }
          }
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
        bool _check_parti(Ancestor& ancestor, const MeshNodeType& mesh_node, bool check_extern)
        {
          // Try to find an appropriate a-priori partitioning first:
          if(check_extern && this->_check_parti_extern(ancestor))
            return true;

          if(this->_check_parti_2level(ancestor, mesh_node))
            return true;

          // No a-priori partitioning found, so we need to determine the
          // required partitioning level for a-posteriori partitioning.
          // For this, first determine the factor by which the number of elements
          // increases for each refinement:
          const Index factor = Index(Geometry::Intern::StandardRefinementTraits<typename BaseClass::ShapeType, BaseClass::ShapeType::dimension>::count);

          // compute minimum number of elements for a-posteriori partitioning
          const Index min_elems = Index(ancestor.num_parts) * Index(_required_elems_per_rank);

          // Okay, get the number of elements on base-mesh level 0
          Index num_elems = mesh_node.get_mesh()->get_num_elements();

          // compute refinement level on which we have enough elements
          int level = 0;
          for(; num_elems < min_elems; ++level)
          {
            num_elems *= factor;
          }

          // finally, the user may have ancestor a greater level for partitioning:
          ancestor.parti_level = Math::max(ancestor.parti_level, level);

          // No a-priori partitioning; try a-posteriori partitioner later
          return false;
        }

        /**
         * \brief Applies an a-posteriori partitioner.
         *
         * \param[inout] ancestor
         * The ancestor object for this layer.
         *
         * \param[in] base_mesh_node
         * The base-mesh node that is to be partitioned.
         *
         * \returns
         * \c true, if the base-mesh was partitioned, otherwise \c false.
         */
        bool _apply_parti(Ancestor& ancestor, /*const*/ MeshNodeType& base_mesh_node)
        {
          // First of all, check whether an a-priori partitioning was selected;
          // if so, then we do not have to apply any a-posteriori partitioner
          if(ancestor.parti_apriori)
            return true;

          // ensure that the mesh has enough elements
          XASSERT(ancestor.num_parts <= int(base_mesh_node.get_mesh()->get_num_elements()));

          // try the various a-posteriori partitioners
          if(this->_apply_parti_metis(ancestor, base_mesh_node))
            return true;
          if(this->_apply_parti_genetic(ancestor, base_mesh_node))
            return true;
          if(this->_apply_parti_naive(ancestor, base_mesh_node))
            return true;

          // we should not arrive here...
          return false;
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

        /**
         * \brief Applies the METIS partitioner onto the base-mesh.
         *
         * \param[inout] ancestor
         * The ancestor object for this layer.
         *
         * \param[in] base_mesh_node
         * The base-mesh node that is to be partitioned.
         *
         * \returns
         * \c true, if METIS was applied successfully, otherwise \c false.
         */
        bool _apply_parti_metis(Ancestor& ancestor, const MeshNodeType& base_mesh_node)
        {
#if defined(FEAT_HAVE_METIS) || defined(FEAT_HAVE_PARMETIS)
          // is this even allowed?
          if(!this->_allow_parti_metis)
            return false;

          // build element adjacency graph
          // connectivity by facets
          /// \todo dirk: use centralized method for adj graph retrieval
          /// \todo dirk: does any partitioner need self-adjacencies? -> remove it in creation
          ///       peter: yes, other partitioners need self-adjacencies
          const auto dimension = MeshNodeType::MeshType::ShapeType::dimension;
          Adjacency::Graph facets_at_elem(Adjacency::RenderType::as_is, base_mesh_node.get_mesh()->template get_index_set<dimension, dimension-1>());
          Adjacency::Graph elems_at_facet(Adjacency::RenderType::transpose, facets_at_elem);
          Adjacency::Graph adj_graph = Adjacency::Graph(Adjacency::RenderType::injectify, facets_at_elem, elems_at_facet);

          // create a metis partitioner
          Geometry::PartiMetis<MeshType> partitioner(adj_graph, Index(ancestor.num_parts));

          // create elems-at-rank graph
          ancestor.parti_graph = partitioner.build_elems_at_rank();

          // set info string
          ancestor.parti_info = String("Applied METIS partitioner");

          // okay
          return true;
#else
          (void)ancestor;
          (void)base_mesh_node;
          return false;
#endif // defined(FEAT_HAVE_METIS) || defined(FEAT_HAVE_PARMETIS)
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
        bool _apply_parti_genetic(Ancestor& ancestor, /*const*/ MeshNodeType& base_mesh_node)
        {
          // is this even allowed?
          if(!this->_allow_parti_genetic)
            return false;

          // create a genetic partitioner
          Geometry::PartiIterative<MeshType> partitioner(
            *base_mesh_node.get_mesh(),
            ancestor.progeny_comm,
            (Index)ancestor.num_parts,
            this->_genetic_time_init,
            this->_genetic_time_mutate);

          // create elems-at-rank graph
          ancestor.parti_graph = partitioner.build_elems_at_rank();

          // set info string
          ancestor.parti_info = String("Applied genetic partitioner");

          // okay
          return true;
        }

        /**
         * \brief Applies the naive partitioner onto the base-mesh.
         *
         * \param[inout] ancestor
         * The ancestor object for this layer.
         *
         * \param[in] mesh_node
         * The base-mesh node that is to be partitioned.
         *
         * \returns
         * \c true, if a naive partition was created successfully, otherwise \c false.
         */
        bool _apply_parti_naive(Ancestor& ancestor, const MeshNodeType& mesh_node)
        {
          // is this even allowed?
          if(!this->_allow_parti_naive)
            return false;

          const Index num_parts = Index(ancestor.num_parts);
          const Index num_elems = mesh_node.get_mesh()->get_num_elements();
          XASSERTM(num_parts <= num_elems, "Base-Mesh does not have enough elements");

          // create elems-at-rank graph
          ancestor.parti_graph = Adjacency::Graph(num_parts, num_elems, num_elems);
          Index* ptr = ancestor.parti_graph.get_domain_ptr();
          Index* idx = ancestor.parti_graph.get_image_idx();
          ptr[0] = Index(0);
          for(Index i(1); i < num_parts; ++i)
            ptr[i] = (i*num_elems) / num_parts;
          ptr[num_parts] = num_elems;
          for(Index j(0); j < num_elems; ++j)
            idx[j] = j;

          // set info string
          ancestor.parti_info = String("Applied naive partitioner");

          // okay
          return true;
        }
#endif // defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
      }; // class PartiDomainControl<...>
    } // namespace Domain
  } // namespace Control
} // namespace FEAT

#endif // CONTROL_DOMAIN_PARTI_DOMAIN_CONTROL_HPP
