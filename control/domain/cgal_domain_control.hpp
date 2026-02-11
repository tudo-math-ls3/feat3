// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
#pragma once

#include <kernel/base_header.hpp>

#include <kernel/util/dist.hpp>
#include <kernel/util/dist_file_io.hpp>
#include <kernel/util/stop_watch.hpp>
#include <kernel/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/property_map.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/partition_set.hpp>
#include <kernel/geometry/parti_2lvl.hpp>
#include <kernel/geometry/parti_iterative.hpp>
#include <kernel/geometry/parti_parmetis.hpp>
#include <kernel/geometry/parti_zoltan.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/patch_mesh_factory.hpp>
#include <kernel/geometry/patch_meshpart_factory.hpp>
#include <kernel/geometry/patch_meshpart_splitter.hpp>
#include <kernel/geometry/cgal.hpp>
#include <kernel/trafo/standard/mapping.hpp>

#include <control/domain/parti_domain_control_base.hpp>

#if defined(FEAT_HAVE_CGAL) || defined(DOXYGEN)
namespace FEAT
{
  namespace Control
  {
    namespace Domain
    {
      /**
       * \brief Wrapper class for domain levels for VoxelDomainControl
       *
       * This class acts as a derived class for a domain level class and extends it by an additional "slag level"
       * object, which is required to assemble grid transfer operators (i.e. prolongation/restriction) for spaces
       * organized in the VoxelDomainControl class. Furthermore, this class also stores the optional element coloring
       * and element layering for parallel assembly.
       *
       * \author Maximilian Esser
       */
      template<typename DomainLevel_>
      class CGALDomainLevelWrapper :
        public DomainLevel_
      {
      public:
        /// our base-class; the actual domain level
        typedef DomainLevel_ BaseClass;
        /// our mesh type
        using typename BaseClass::MeshType;
        /// our mesh node type
        using typename BaseClass::MeshNodeType;

        /// attachable slag level
        std::shared_ptr<DomainLevel_> slag_level;

        /// element coloring
        Adjacency::Coloring element_coloring;

        explicit CGALDomainLevelWrapper(int _lvl_idx, std::unique_ptr<MeshNodeType> node,
          Adjacency::Coloring&& coloring) :
          BaseClass(_lvl_idx, std::move(node)),
          slag_level(),
          element_coloring(std::forward<Adjacency::Coloring>(coloring))
        {
        }

        explicit CGALDomainLevelWrapper(int _lvl_idx, std::unique_ptr<MeshNodeType> node,
          std::unique_ptr<MeshNodeType> slag_node, Adjacency::Coloring&& coloring) :
          BaseClass(_lvl_idx, std::move(node)),
          slag_level(),
          element_coloring(std::forward<Adjacency::Coloring>(coloring))
        {
          if(slag_node)
            slag_level = std::make_shared<DomainLevel_>(_lvl_idx, std::move(slag_node));
        }
      }; // class CGALDomainLevelWrapper<...>

      /**
       * \brief Hierarchical partitioned Voxel Domain Control
       *
       * \author Peter Zajac, Maximilian Esser
       */
      template<typename DomainLevel_>
      class CGALDomainControl :
        public Control::Domain::PartiDomainControlBase< DomainLevel_ >
      {
      public:
        /// our base class
        typedef Control::Domain::PartiDomainControlBase< DomainLevel_ > BaseClass;
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
        /// weight type for partitioners
        using typename BaseClass::WeightType;
        /// ancestor class
        using typename BaseClass::Ancestor;

        /// shape dimension
        static constexpr int shape_dim = MeshType::shape_dim;
        static_assert(shape_dim == 3);

        /// coordinate type
        typedef typename MeshType::CoordType CoordType;

      protected:
        /// creation stage to keep track what has already been initialized
        int _create_stage;

        /// the bounding box of the structured domain
        Tiny::Vector<CoordType, shape_dim> _bbox_min, _bbox_max;
        /// number of element slices on base mesh level 0
        std::array<Index, shape_dim> _base_slices;

        /// cgal wrapper
        std::unique_ptr<Geometry::CGALWrapper<CoordType>> _cgal_wrapper;

        /// the input base-mesh node on level 0
        std::unique_ptr<MeshNodeType> _base_mesh_node;

        /// a bunch of stop-watches
        mutable StopWatch _watch_create_base_mesh, _watch_create_cgal_wrapper, _watch_create_hierarchy, _watch_calc_weights;

        /// specifies whether to keep the voxel map after hierarchy creation
        bool _keep_cgal_wrapper;

        /// specifies whether we use a non voxel base mesh
        bool _unstructered_mesh = false;

        bool _cheap_weights = true;

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
        explicit CGALDomainControl(const Dist::Comm& comm_, bool support_multi_layered) :
          BaseClass(comm_, support_multi_layered),
          _create_stage(0),
          _keep_cgal_wrapper(false),
          _unstructered_mesh(false),
          _cheap_weights(true)
        {
        }

        /// virtual destructor
        virtual ~CGALDomainControl()
        {
        }

        /**
         * \brief Instructs the domain controller to keep the voxel map after hierarchy creation.
         *
         * Unless this function is called prior to the call of the create_hierarchy() function,
         * the voxel map will be released after the hierarchy creation is completed to free up the
         * memory used by the voxel map.
         *
         * It is recommended not to call this function unless it is really necessary; if you want
         * to use the voxel map for anything after it has been created, e.g. to write it to a file
         * or to print some statistics to cout, then you can use it before calling the create_hierarchy
         * function.
         */
        void keep_cgal_wrapper()
        {
          this->_keep_cgal_wrapper = true;
        }

        /// set whether we use an expensive or cheap method to calculate the element weights
        void cheap_weight_calculation(bool option)
        {
          this->_cheap_weights = option;
        }

        Tiny::Vector<CoordType, shape_dim> get_bounding_box_min() const
        {
          return _bbox_min;
        }

        Tiny::Vector<CoordType, shape_dim> get_bounding_box_max() const
        {
          return _bbox_max;
        }

        /**
         * \brief Creates a 3D structured cuboid base mesh and returns its base-mesh node
         *
         * This functions returns a reference to the base-mesh node that is created and managed internally
         * by this class. This base-mesh node reference can be used to create mesh-parts at the base-mesh level,
         * which will be refined and deslagged during the creation of the domain hierarchy.
         *
         * \attention The returned base-mesh node reference is only valid until the create() method is called and
         * is orphaned after the create() method returns, because the base-mesh node is then processed internally
         * to create the domain level hierarchy.
         *
         * \param[in] num_x, num_y, num_z
         * The number of element slices in X-, Y- and Z-direction.
         *
         * \param[in] x_min, x_max
         * The X-range of the cuboid domain.
         *
         * \param[in] y_min, y_max
         * The X-range of the cuboid domain.
         *
         * \param[in] z_min, z_max
         * The X-range of the cuboid domain.
         *
         * \returns A reference to the internal base-mesh node containing the base-mesh.
         */
        MeshNodeType& create_base_mesh_3d(Index num_x, Index num_y, Index num_z, CoordType x_min, CoordType x_max,
          CoordType y_min, CoordType y_max, CoordType z_min, CoordType z_max)
        {
          XASSERTM(_create_stage == 0, "base mesh already created!");

          XASSERT(num_x > Index(0));
          XASSERT(num_y > Index(0));
          XASSERT(num_z > Index(0));
          XASSERT(x_min < x_max);
          XASSERT(y_min < y_max);
          XASSERT(z_min < z_max);

          _watch_create_base_mesh.start();

          _base_slices[0] = num_x;
          _base_slices[1] = num_y;
          _base_slices[2] = num_z;
          _bbox_min[0] = x_min;
          _bbox_max[0] = x_max;
          _bbox_min[1] = y_min;
          _bbox_max[1] = y_max;
          _bbox_min[2] = z_min;
          _bbox_max[2] = z_max;

          // create a structured mesh from factory
          _base_mesh_node = MeshNodeType::make_unique(
            Geometry::StructUnitCubeFactory<MeshType>::make_unique_from(num_x, num_y, num_z), &this->_atlas);

          // get the node's vertex set
          auto& vtx = _base_mesh_node->get_mesh()->get_vertex_set();
          Index n = vtx.get_num_vertices();
          for(Index i(0); i < n; ++i)
          {
            vtx[i][0] = x_min + (x_max - x_min) * vtx[i][0];
            vtx[i][1] = y_min + (y_max - y_min) * vtx[i][1];
            vtx[i][2] = z_min + (z_max - z_min) * vtx[i][2];
          }

          _watch_create_base_mesh.stop();

          // update creation stage
          this->_create_stage = 1;

          // return the base-mesh node
          return *_base_mesh_node.get();
        }

        /**
         * \brief Set another mesh as base mesh
         *
         * \param[in] mesh_node Unique ptr to mesh node that is to be consumed by this domain
         */
        MeshNodeType& set_base_mesh(std::unique_ptr<MeshNodeType>&& mesh_node)
        {
          XASSERTM(_create_stage == 0, "base mesh already created!");

          _watch_create_base_mesh.start();
          _base_mesh_node = std::move(mesh_node);

          for(std::size_t k = 0; k < std::size_t(shape_dim); ++k)
          {
            _base_slices[k] = Math::Limits<Index>::max();
          }

          for(int k = 0; k < shape_dim; ++k)
          {
            _bbox_min[k] = Math::Limits<CoordType>::max();
            _bbox_max[k] = Math::Limits<CoordType>::lowest();
          }

          auto& vtx = _base_mesh_node->get_mesh()->get_vertex_set();
          Index n = vtx.get_num_vertices();
          for(Index i(0); i < n; ++i)
          {
            for(int k = 0; k < shape_dim; ++k)
            {
              _bbox_min[k] = Math::min(_bbox_min[k], vtx[i][k]);
              _bbox_max[k] = Math::max(_bbox_max[k], vtx[i][k]);
            }
          }

          _watch_create_base_mesh.stop();

          // update creation stage
          this->_create_stage = 1;

          // update unstructered mesh flag
          this->_unstructered_mesh = true;

          // return the base-mesh node
          return *_base_mesh_node.get();
        }

        /**
         * \brief Creates a voxel map based on a VoxelMasker object
         *
         * This function can only be called after the base-mesh has been created by calling
         * create_base_mesh_3d or set_base_mesh.
         *
         * \attention
         * This function is a "collective" function in distributed parallel (i.e. MPI) simulations and each process
         * must be given the same file content, else this will lead to incostinet inside out information
         *
         * \param[in] file
         * A \transient reference to the input stream
         *
         * \param[in] file_mode
         * The file mode of the cgal input, should be either off or obj
         *
         * \compilerhack avoid std::make_unique for clang
         */
        void create_cgal_wrapper(std::istream& file, Geometry::CGALFileMode file_mode)
        {
          XASSERTM(this->_create_stage <= 1, "invalid creation stage: hirarchy already created");

          this->_watch_create_cgal_wrapper.start();
#ifdef FEAT_COMPILER_CLANG
          this->_cgal_wrapper.reset(new Geometry::CGALWrapper<CoordType>(file, file_mode));
#else
          this->_cgal_wrapper = std::make_unique<Geometry::CGALWrapper<CoordType>>(file, file_mode);
#endif
          this->_watch_create_cgal_wrapper.stop();
        }

        /// has to be called from each process part of the domain hirarchy
        void create_cgal_wrapper(const FEAT::String& file, Geometry::CGALFileMode file_mode)
        {
          // read in istream in an mpi friendly matter
          this->_watch_create_cgal_wrapper.start();
          std::stringstream cgal_input;
          DistFileIO::read_common(cgal_input, file, this->_comm);
          this->_watch_create_cgal_wrapper.stop();
          create_cgal_wrapper(cgal_input, file_mode);
        }

        /**
         * \brief Creates the domain level hierarchy
         *
         * This function can only be called after the base-mesh has been created by calling create_base_mesh_2d or
         * create_base_mesh_3d and after the voxel map has been created by calling create_voxel_map().
         */
        void create_hierarchy()
        {
          XASSERTM(this->_create_stage >= 1, "invalid creation stage; create the base-mesh and set the cgal wrapper first");
          XASSERTM(this->_create_stage <= 1, "invalid creation stage: domain hierarchy already created");

          this->_watch_create_hierarchy.start();

          // create the domain control
#ifdef FEAT_HAVE_MPI
          if(this->_comm.size() == 1)
          {
            // We've got just one process, so it's a simple choice:
            this->_create_single_process(std::move(this->_base_mesh_node));
          }
          else if(this->_support_multi_layered && (this->_desired_levels.size() > std::size_t(2)))
          {
            // The application supports multi-layered domain controls and
            // the user wants that, so create a multi-layered one:
            this->_create_multi_layered(std::move(this->_base_mesh_node));
          }
          else
          {
            // Create a single-layered domain control:
            this->_create_single_layered(std::move(this->_base_mesh_node));
          }
#else // not FEAT_HAVE_MPI
          {
            // In the non-MPI case, we always have only one process:
            this->_create_single_process(std::move(this->_base_mesh_node));
          }
#endif // FEAT_HAVE_MPI

          // compile virtual levels
          this->compile_virtual_levels();

          this->_watch_create_hierarchy.stop();

          // collect statistics
          FEAT::Statistics::toe_partition = this->_watch_create_base_mesh.elapsed() +
            this->_watch_create_cgal_wrapper.elapsed() + this->_watch_create_hierarchy.elapsed();

          // update creation stage
          this->_create_stage = 2;
          this->_was_created = true;

          // delete the voxel map unless we're meant to keep it
          if(!this->_keep_cgal_wrapper)
            this->_cgal_wrapper.reset();
        }

        /**
         * \brief Releases the internal cgal wrapper
         *
         * \returns A pointer to the heap allocated cgal wrapper
         *
         * \attention After calling this function YOU are responsible to
         *            free the cgal wrapper afterwards
         */
        Geometry::CGALWrapper<CoordType>* release_cgal_wrapper()
        {
          return this->_cgal_wrapper.release();
        }

        /// \returns a const reference to the internal cgal wrapper
        Geometry::CGALWrapper<CoordType>* get_cgal_wrapper()
        {
          return this->_cgal_wrapper.get();
        }

        /// \returns a const reference to the internal cgal wrapper
        const Geometry::CGALWrapper<CoordType>* get_cgal_wrapper() const
        {
          return this->_cgal_wrapper.get();
        }

        void set_cgal_wrapper(std::unique_ptr<Geometry::CGALWrapper<CoordType>>&& cgal_wrapper)
        {
          this->_cgal_wrapper = std::move(cgal_wrapper);
        }

        bool cgal_wrapper_empty() const
        {
          return !this->_cgal_wrapper;
        }

        /// Returns a const reference to the StopWatch that measures the base-mesh creation phase.
        const StopWatch& get_watch_base_mesh() const
        {
          return this->_watch_create_base_mesh;
        }

        /// Returns a const reference to the StopWatch that measures the voxel-map creation phase.
        const StopWatch& get_watch_cgal_wrapper() const
        {
          return this->_watch_create_cgal_wrapper;
        }

        /// Returns a const reference to the StopWatch that measures the base-mesh creation phase.
        const StopWatch& get_watch_hierarchy() const
        {
          return this->_watch_create_hierarchy;
        }

        const StopWatch& get_watch_weights() const
        {
          return this->_watch_calc_weights;
        }

        /**
         * \brief Gathers the element map weights on a given domain level.
         *
         * The "slag weight" of an element is defined as the integral of the step function defined by the inside out test,
         * divided by the cell volume.
         *
         * \param[in] level
         * The physical domain level for which the element weights are to be gathered.
         *
         * \returns A vector of floats in range [0,1] representing the element slag weights on the desired level.
         * The length of the vector corresponds to the number of elements on the patch on the desired domain level.
         */
        std::vector<WeightType> gather_element_weights(Index level) const
        {
          XASSERT(level < this->size_physical());
          return this->_gather_element_weights(this->at(level)->get_mesh());
        }

        /**
         * \brief Debugging function: Returns a string containing encoded slag layer level information
         */
        String dump_slag_layer_levels() const
        {
          String msg;
          for(auto it = this->_layer_levels.begin(); it != this->_layer_levels.end(); ++it)
          {
            if(it != this->_layer_levels.begin())
              msg += " |";
            for(auto jt = it->begin(); jt != it->end(); ++jt)
            {
              msg += " " + stringify((*jt)->get_level_index());
              msg += "[";
              msg += stringify((*jt)->get_mesh().get_num_elements()).pad_front(4);
              msg += "]";
              if((*jt)->slag_level)
                ((msg += "<[") += stringify((*jt)->slag_level->get_mesh().get_num_elements()).pad_front(4)) += "]";
            }
          }
          return msg;
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
          //dummy implementation
          ASSERTM(false, "Thou shall not arrive here!");
          (void) filenames;
          (void) dirpath;
        }

        /**
         * \brief Creates a domain control from a MeshFileReader.
         *
         * \param[inout] mesh_reader
         * A \transient reference to the mesh reader from which the domain control is to be created.
         */
        virtual void create(Geometry::MeshFileReader&)
        {
          //dummy implementation
          ASSERTM(false, "Thou shall not arrive here!");
        }

        /**
         * \brief Creates a domain control from a base-mesh node
         *
         * \param[in] base_mesh_node
         * The base-mesh node that the domain control is to be created from
         */
        virtual void create(std::unique_ptr<MeshNodeType>&)
        {
          //dummy implementation
          ASSERTM(false, "Thou shall not arrive here!");
        }

      protected:
        /**
         * \brief Gather a vector of element indices that intersect the masked domain
         *
         * \attention
         * It is silently assumed that all elements in the mesh are axis-parallel rectangles/cuboids; this function
         * may yield incorrect results if the mesh violates this assumption!
         *
         * \param[in] mesh
         * A \transient reference to the mesh whose elements are to be tested against the cgal mesh
         *
         * \returns
         * A vector containing the indices of all mesh elements that intersect the domain that is represented by
         * the voxel map.
         */
        std::vector<Index> _gather_masked_elements(const MeshType& mesh) const
        {
          const Index num_elems = mesh.get_num_elements();

          // reserve the vector
          std::vector<Index> masked_elems;
          masked_elems.reserve(num_elems);

          // get the vertex set and the vertices-at-element index set
          const auto& vtx = mesh.get_vertex_set();
          const auto& verts_at_elem = mesh.template get_index_set<shape_dim, 0>();

          // loop over all elements
          FEAT_PRAGMA_OMP(parallel for)
          for(Index ielem = 0; ielem < num_elems; ++ielem)
          {
            std::array<Tiny::Vector<CoordType, shape_dim>, std::decay_t<decltype(verts_at_elem)>::num_indices> verts;
            for(int k = 0; k < int(verts.size()); ++k)
            {
              verts[k] = vtx[verts_at_elem(ielem, k)];
            }

            if(this->_cgal_wrapper->intersects_polygon(verts))
              masked_elems.push_back(ielem);
          }

          // return the list of all elements that we found
          return masked_elems;
        }

        /**
         * \brief Gather a vector of element indices that intersect the masked domain in a copperative manner
         *
         * \attention
         * This has to be called by all processes in the sibling communicator
         *
         * \param[in] mesh
         * A \transient reference to the mesh whose elements are to be tested against the cgal mesh
         *
         * \param[in] sibling_comm
         * The communicator accross the workload is to be spread.
         *
         * \returns
         * A vector containing the indices of all mesh elements that intersect the domain that is represented by
         * the voxel map.
         */
        std::vector<Index> _gather_masked_elements(const MeshType& mesh, const Dist::Comm& sibling_comm) const
        {
          const Index num_elems = mesh.get_num_elements();
          const Index num_procs = Index(sibling_comm.size());
          const Index cur_rank = Index(sibling_comm.rank());
          const Index elems_per_rank = num_elems/num_procs + Index((num_elems%num_procs)>0u);

          // reserve the vector
          std::vector<Index> masked_elems;
          masked_elems.reserve(elems_per_rank);

          // get the vertex set and the vertices-at-element index set
          const auto& vtx = mesh.get_vertex_set();
          const auto& verts_at_elem = mesh.template get_index_set<shape_dim, 0>();

          // loop over all elements
          for(Index ielem = cur_rank*elems_per_rank; ielem < Math::min((cur_rank+1u)*elems_per_rank, num_elems); ++ielem)
          {
            std::array<Tiny::Vector<CoordType, shape_dim>, std::decay_t<decltype(verts_at_elem)>::num_indices> verts;
            for(int k = 0; k < int(verts.size()); ++k)
            {
              verts[k] = vtx[verts_at_elem(ielem, k)];
            }

            if(this->_cgal_wrapper->intersects_polygon(verts))
              masked_elems.push_back(ielem);
          }

          // communicate sizes of buffers
          std::vector<int> recv_sizes(num_procs);
          recv_sizes.at(cur_rank) = int(masked_elems.size());

          sibling_comm.allgather(recv_sizes.data(), 1, recv_sizes.data(), 1);

          // and now, create big buffer
          std::vector<Index> gathered_mask(std::size_t(std::accumulate(recv_sizes.begin(), recv_sizes.end(), 0)), 0u);

          // and also, we have to caclulate our displacements into our big buffer
          std::vector<int> displacements;
          displacements.reserve(recv_sizes.size());
          std::exclusive_scan(recv_sizes.begin(), recv_sizes.end(), std::back_inserter(displacements), 0);

          // and now gather
          sibling_comm.allgatherv(masked_elems.data(), masked_elems.size(), gathered_mask.data(), recv_sizes.data(), displacements.data());

          // return the list of all elements that we found
          return gathered_mask;
        }

        /**
         * \brief Gathers the element slag weights for a given mesh
         *
         * The "slag weight" of an element is defined as the integral of the step function defined by the inside out test,
         * divided by the cell volume.
         *
         * \param[in] mesh
         * The mesh for which the element weights are to be gathered.
         *
         * \returns A vector of floats in range [0,1] representing the element slag weights on the desired level.
         * The length of the vector corresponds to the number of elements on the patch on the desired domain level.
         */
        std::vector<WeightType> _gather_element_weights(const MeshType& mesh) const
        {
          this->_watch_calc_weights.start();
          const Index num_elems = mesh.get_num_elements();

          // allocate the weight vector
          std::vector<WeightType> weights(num_elems, 0.0);

          if(!_cgal_wrapper)
            return weights;

          Cubature::DynamicFactory cub_fac("refine*4:gauss-legendre:3");
          Cubature::Rule<typename MeshType::ShapeType> rule;
          cub_fac.create(rule);

          // guarantee that cgal wrapper is initialized
          this->_cgal_wrapper->point_not_outside(CoordType(0), CoordType(0), CoordType(0));

          // create trafo
          typedef Trafo::Standard::Mapping<MeshType> TrafoType;
          typedef typename TrafoType::template Evaluator<typename MeshType::ShapeType, CoordType>::Type TrafoEvaluator;
          static constexpr TrafoTags trafo_config = TrafoTags::dom_point | TrafoTags::img_point | TrafoTags::jac_det;

          /// trafo evaluation data type
          typedef typename TrafoEvaluator::template ConfigTraits<trafo_config>::EvalDataType TrafoEvalData;

          const TrafoType trafo(const_cast<MeshType&>(mesh));

          FEAT_PRAGMA_OMP(parallel)
          {
            TrafoEvaluator trafo_eval(trafo);
            TrafoEvalData trafo_data;

            // loop over all elements
            FEAT_PRAGMA_OMP(for)
            for(Index ielem = 0; ielem < num_elems; ++ielem)
            {
              trafo_eval.prepare(ielem);
              CoordType loc_sum = CoordType(0);
              for(int k = 0; k < rule.get_num_points(); ++k)
              {
                const auto ref_point = rule.get_point(k);
                trafo_eval(trafo_data, ref_point);
                const auto& img_point = trafo_data.img_point;
                const CoordType weight = trafo_data.jac_det * rule.get_weight(k);
                loc_sum += weight * CoordType(this->_cgal_wrapper->point_not_outside(img_point[0], img_point[1], img_point[2]));
              }
              loc_sum /= trafo_eval.volume();
              ASSERTM(loc_sum >= CoordType(0) && loc_sum <= CoordType(1), "Local sum not between zero and one");

              // gather element weight based on the bounding box
              weights[ielem] = loc_sum;
            }
          }
          this->_watch_calc_weights.stop();

          // returns the element weights
          return weights;
        }

        std::vector<WeightType> _compute_weights(Ancestor& DOXY(ancestor), const MeshNodeType& base_mesh_node) override
        {
          std::vector<WeightType> weights;
          if(!this->_cheap_weights)
            weights = this->_gather_element_weights(*base_mesh_node.get_mesh());

          return weights;
        }

        /**
         * \brief Deslags an unpartitioned mesh node including its mesh-parts
         *
         * \attention
         * This function can only be used for unpartitioned mesh nodes. If the mesh node represents a patch, then
         * the other overload that requires and Ancestor object has to be used, because the halos also need to be
         * deslagged in this case!
         *
         * This function adds a patch with rank -1 to the mesh_node, which represents the deslagged patch.
         *
         * \param[inout] mesh_node
         * A \transient reference to the mesh node that is to be deslagged.
         *
         * \returns A new mesh node that contains the deslagged mesh as well as the deslagged mesh parts.
         */
        virtual std::unique_ptr<MeshNodeType> _deslag_mesh_node(MeshNodeType& mesh_node)
        {
          XASSERTM(mesh_node.get_halo_map().empty(), "This function must not be used for partitioned mesh nodes!");
          // it is always silenetly assumed that we call this function on the complete communicator
          return mesh_node.extract_patch(_gather_masked_elements(*mesh_node.get_mesh(), this->_comm), true, false, false);
        }

        /**
         * \brief Deslags a (potentially partitioned) mesh node including its mesh-parts and halos
         *
         * This function adds a patch with rank -1 to the mesh_node, which represents the deslagged patch.
         *
         * \param[inout] mesh_node
         * A \transient reference to the mesh node that is to be deslagged.
         *
         * \param[in] ancestor
         * Depending on whether is_child is \c false or \c true, a \transient reference to the ancestor object of either
         * the layer that the mesh node belongs to or the child layer of the layer that the mesh node belongs to, resp.
         *
         * \param[in] is_child
         * Specifies whether the ancestor object is actually the child ancestor object; this is only \c true for the
         * mesh nodes belonging to parent layers in the case of hierarchical partitioning.
         *
         * \returns A new mesh node that contains the deslagged mesh as well as the deslagged mesh parts and halos.
         */
        virtual std::unique_ptr<MeshNodeType> _deslag_mesh_node(MeshNodeType& mesh_node, const Ancestor& ancestor, bool is_child)
        {
          const auto* prod_comm = is_child ? &ancestor.progeny_comm : nullptr;
          std::unique_ptr<MeshNodeType> new_node = prod_comm ? mesh_node.extract_patch(_gather_masked_elements(*mesh_node.get_mesh(), *prod_comm), true, false, false):
                                                               mesh_node.extract_patch(_gather_masked_elements(*mesh_node.get_mesh()), true, false, false);
          this->_deslag_patch_halos(*new_node, mesh_node, ancestor, is_child);
          return new_node;
        }

        /**
         *  \brief Deslags the halos of a patch mesh node.
         *
         * \param[in] patch_mesh_node
         * A \transient reference to the deslagged patch mesh node whose halos are to be deslagged.
         *
         * \param[in] base_mesh_node
         * A \transient reference to the mesh node that the deslagged mesh node was created from.
         *
         * \param[in] ancestor
         * Depending on whether is_child is \c false or \c true, a \transient reference to the ancestor object of either
         * the layer that the mesh node belongs to or the child layer of the layer that the mesh node belongs to, resp.
         *
         * \param[in] is_child
         * Specifies whether the ancestor object is actually the child ancestor object; this is only \c true for the
         * mesh nodes belonging to parent layers in the case of hierarchical partitioning.
         */
        virtual void _deslag_patch_halos(
          MeshNodeType& patch_mesh_node,
          const MeshNodeType& base_mesh_node,
          const Ancestor& ancestor,
          bool is_child)
        {
          // get the map of the base-mesh halos
          const std::map<int, std::unique_ptr<MeshPartType>>& base_halo_map = base_mesh_node.get_halo_map();

          // if the base mesh has no halos, then we can jump out of here
          if(base_halo_map.empty())
            return;

          // get number of halos
          const std::size_t num_halos = base_halo_map.size();

          // create a halo splitter
          Geometry::PatchHaloSplitter<MeshType> halo_splitter(*base_mesh_node.get_mesh(), *base_mesh_node.get_patch(-1));

          // add each base-mesh halo to our halo splitter and store the resulting split data size
          std::vector<int> halo_ranks;
          std::vector<std::size_t> halo_send_sizes;
          for(auto it = base_halo_map.begin(); it != base_halo_map.end(); ++it)
          {
            // store halo rank
            halo_ranks.push_back(it->first);

            // add halo and compute send buffer size
            halo_send_sizes.push_back(halo_splitter.add_halo(it->first, *it->second));
          }

          // This vector will receive the split halo data from all our potential neighbor processes
          std::vector<std::size_t> halo_recv_sizes(num_halos);
          std::vector<std::vector<Index>> halo_send_data(num_halos), halo_recv_data(num_halos);
          Dist::RequestVector halo_recv_reqs(num_halos), halo_send_reqs(num_halos);

          // get the layer index of our comm layer; the ancestor layer unless we got the child ancestry during
          // multi-layered partitioning; this may be equal to -1 in either case if this process does not participate
          // in that corresponding layer
          const int layer_idx = is_child ? ancestor.layer_p : ancestor.layer;

          // split and serialize halos and store sizes
          for(std::size_t i(0); i < num_halos; ++i)
          {
            // serialize the split halo into the send data buffer
            if(halo_send_sizes.at(i) > std::size_t(0))
              halo_send_data.at(i) = halo_splitter.serialize_split_halo(halo_ranks[i], -1);
          }

          // exchange halo send data sizes over the corresponding layer communicator
          if(layer_idx >= 0)
          {
            // get the corresponding layer communicator
            const Dist::Comm& layer_comm = this->_layers.at(std::size_t(layer_idx))->comm();

            // exchange halos sizes
            for(std::size_t i(0); i < num_halos; ++i)
            {
              // post receive requests
              halo_recv_reqs[i] = layer_comm.irecv(&halo_recv_sizes[i], std::size_t(1), halo_ranks[i]);

              // post send requests
              halo_send_reqs[i] = layer_comm.isend(&halo_send_sizes[i], std::size_t(1), halo_ranks[i]);
            }

            // wait for sends and receives to finish
            halo_recv_reqs.wait_all();
            halo_send_reqs.wait_all();

            // exchange halo data
            for(std::size_t i(0); i < num_halos; ++i)
            {
              // do we receive any data from our neighbor?
              if(halo_recv_sizes[i] > Index(0))
              {
                // resize buffer and post receive
                halo_recv_data.at(i).resize(halo_recv_sizes[i]);
                halo_recv_reqs[i] = layer_comm.irecv(halo_recv_data.at(i).data(), halo_recv_sizes.at(i), halo_ranks[i]);
              }

              // do we send any data to our neighbor?
              if(halo_send_sizes.at(i) > Index(0))
              {
                // post send of actual halo buffer
                halo_send_reqs[i] = layer_comm.isend(halo_send_data.at(i).data(), halo_send_sizes.at(i), halo_ranks[i]);
              }
            }

            // wait for sends and receives to finish
            halo_recv_reqs.wait_all();
            halo_send_reqs.wait_all();
          }

          // broadcast halo receive data over progeny comm in case of multi-layered partitioning
          if(is_child)
          {
            // broadcast sizes
            ancestor.progeny_comm.bcast(halo_recv_sizes.data(), num_halos, 0);

            // allocate buffers
            if(ancestor.progeny_comm.rank() != 0)
            {
              for(std::size_t i(0); i < num_halos; ++i)
                halo_recv_data.at(i).resize(halo_recv_sizes[i]);
            }

            // broadcast all halos
            for(std::size_t i(0); i < num_halos; ++i)
            {
              if(halo_recv_sizes[i] > std::size_t(0))
                ancestor.progeny_comm.bcast(halo_recv_data[i].data(), halo_recv_sizes[i], 0);
            }
          }

          // create a vector of split halo sizes
          std::vector<std::array<Index, shape_dim+1>> halo_send_intsec_sizes(num_halos), halo_recv_intsec_sizes(num_halos);

          // process all halos
          for(std::size_t i(0); i < num_halos; ++i)
          {
            // did we receive any data from our neighbor?
            if(halo_recv_sizes.at(i) == Index(0))
              continue;

            // intersect with our other halo
            if(!halo_splitter.intersect_split_halo(halo_ranks[i], halo_recv_data.at(i), 0u))
              continue; // no intersection between halos

            // create the new halo
            std::unique_ptr<MeshPartType> split_halo = halo_splitter.make_unique();

            // store the entity counts of the halo in our split sizes vector
            for(int j(0); j <= shape_dim; ++j)
              halo_send_intsec_sizes[i][Index(j)] = split_halo->get_num_entities(j);

            // create new halo mesh-part
            patch_mesh_node.add_halo(halo_ranks[i], std::move(split_halo));
          }

          // exchange halo intersected halo sizes over progeny comm; we do this just to ensure consistency
          if((ancestor.layer_p >= 0) || (!is_child && (ancestor.layer >= 0)))
          {
            // get the corresponding layer communicator
            const LayerType& layer = *this->_layers.at(std::size_t(is_child ? ancestor.layer_p : ancestor.layer));

            // exchange intersected halo sizes
            for(std::size_t i(0); i < num_halos; ++i)
            {
              halo_recv_reqs[i] = layer.comm().irecv(halo_recv_intsec_sizes.at(i).data(), std::size_t(shape_dim+1), halo_ranks[i]);
              halo_send_reqs[i] = layer.comm().isend(halo_send_intsec_sizes.at(i).data(), std::size_t(shape_dim+1), halo_ranks[i]);
            }

            // wait for sends and receives to finish
            halo_recv_reqs.wait_all();
            halo_send_reqs.wait_all();

            // ensure that the halo sizes are consistent between neighboring processes
            for(std::size_t i(0); i < num_halos; ++i)
            {
              for(int j(0); j <= shape_dim; ++j)
              {
                if(halo_recv_intsec_sizes[i][Index(j)] != halo_send_intsec_sizes[i][Index(j)])
                {
                  String msg = "Inconsistent deslagged halo size between process ";
                  msg += stringify(this->_comm.rank());
                  msg += " and neighbor with layer rank ";
                  msg += stringify(halo_ranks[i]);

                  // abort execution
                  XABORTM(msg.c_str());
                }
              }
            }
          }
        }

        // virtual bool _parse_parti_type(const String& type) override
        // {
        //   if(type.compare_no_case("extern") == 0)
        //     return this->_allow_parti_extern = true;
        //   // else if(type.compare_no_case("2level") == 0)
        //   //   return _allow_parti_2level = true;
        //   else
        //     return BaseClass::_parse_parti_type(type);
        // }

        Adjacency::Coloring _compute_unstructered_coloring(const MeshNodeType& base_node) const
        {
          const auto& verts_at_elem = base_node.get_mesh()->template get_index_set<MeshType::shape_dim, 0>();
          Adjacency::Graph elems_at_vert(Adjacency::RenderType::transpose, verts_at_elem);
          Adjacency::Graph elems_at_elem(Adjacency::RenderType::injectify, verts_at_elem, elems_at_vert);
          return Adjacency::Coloring(elems_at_elem);
        }

        /**
         * \brief Computes the coloring for the unpartitioned base mesh
         *
         * \param[in] slag_node
         * The unpartitioned structured base-mesh node
         *
         * \returns
         * The coloring for the unpartitioned deslagged base-mesh node
         */
        Adjacency::Coloring _compute_base_mesh_coloring(const MeshNodeType& slag_node) const
        {
          // get element target set of slag patch
          const Geometry::MeshPart<MeshType>* slag_part = slag_node.get_patch(-1);
          XASSERT(slag_part != nullptr);
          const Geometry::TargetSet& slag_target = slag_part->template get_target_set<MeshType::shape_dim>();

          // number of colors = 2^shape_dim
          const Index num_colors = Index(1) << MeshType::shape_dim;
          const Index num_elems = slag_target.get_num_entities();

          // allocate coloring
          Adjacency::Coloring coloring(num_elems, num_colors);
          Index* colors = coloring.get_coloring();

          // unrefined base mesh: structured numbering
          if constexpr(MeshType::shape_dim == 1)
          {
            for(Index i(0); i < num_elems; ++i)
              colors[i] = slag_target[i] & 1;
          }
          else if constexpr(MeshType::shape_dim == 2)
          {
            for(Index i(0); i < num_elems; ++i)
            {
              Index iel = slag_target[i];
              Index ix = iel % this->_base_slices[0];
              Index iy = iel / this->_base_slices[0];
              colors[i] = (ix & 1) | ((iy & 1) << 1);
            }
          }
          else if constexpr(MeshType::shape_dim == 3)
          {
            for(Index i(0); i < num_elems; ++i)
            {
              Index iel = slag_target[i];
              Index ix = iel % this->_base_slices[0];
              Index iy = (iel / this->_base_slices[0]) % this->_base_slices[1];
              Index iz = iel / (this->_base_slices[0] * this->_base_slices[1]);
              colors[i] = (ix & 1) | ((iy & 1) << 1) | ((iz & 1) << (1 << 1));
            }
          }

          return coloring;
        }

        /**
         * \brief Computes the coloring for a refined mesh
         *
         * \param[in] slag_node
         * The refined mesh node
         *
         * \returns
         * The coloring for the refined and deslagged mesh node
         */
        Adjacency::Coloring _compute_refined_mesh_coloring(const MeshNodeType& slag_node) const
        {
          // get element target set of slag patch
          const Geometry::MeshPart<MeshType>* slag_part = slag_node.get_patch(-1);
          XASSERT(slag_part != nullptr);
          const Geometry::TargetSet& slag_target = slag_part->template get_target_set<MeshType::shape_dim>();

          // number of colors = 2^shape_dim
          const Index num_colors = Index(1) << MeshType::shape_dim;
          const Index num_elems = slag_target.get_num_entities();

          // allocate coloring
          Adjacency::Coloring coloring(num_elems, num_colors);
          Index* colors = coloring.get_coloring();

          // color of element i = i % num_colors
          for(Index i(0); i < num_elems; ++i)
            colors[i] = slag_target[i] % num_colors;

          return coloring;
        }

        /**
         * \brief Extracts a patch  coloring from the parent mesh coloring
         *
         * \param[in] parent_node
         * The parent node that the patch was extracted from
         *
         * \param[in] parent_coloring
         * The coloring of the parent node
         *
         * \param[in] child_rank
         * The rank of the extracted patch
         *
         * \returns
         * The coloring for the refined deslagged mesh node
         */
        Adjacency::Coloring _extract_patch_coloring(const MeshNodeType& slag_node, const Adjacency::Coloring& parent_coloring, int child_rank) const
        {
          // get element target set of patch
          const Geometry::MeshPart<MeshType>* patch_part = slag_node.get_patch(child_rank);
          XASSERT(patch_part != nullptr);
          const Geometry::TargetSet& patch_target = patch_part->template get_target_set<MeshType::shape_dim>();

          // number of colors = 2^shape_dim
          const Index num_colors = parent_coloring.get_num_colors();
          const Index num_elems = patch_target.get_num_entities();

          // allocate coloring
          Adjacency::Coloring coloring(num_elems, num_colors);
          Index* colors = coloring.get_coloring();

          // color of element i = i % num_colors
          for(Index i(0); i < num_elems; ++i)
            colors[i] = parent_coloring[patch_target[i]];

          return coloring;
        }

        /**
         * \brief Computes the element layering for the unpartitioned base mesh
         *
         * \param[in] slag_node
         * The unpartitioned structured base-mesh node
         *
         * \returns
         * The layering for the unpartitioned deslagged base-mesh node
         */
        Adjacency::Coloring _compute_base_mesh_layering(const MeshNodeType& slag_node) const
        {
          // get element target set of slag patch
          const Geometry::MeshPart<MeshType>* slag_part = slag_node.get_patch(-1);
          XASSERT(slag_part != nullptr);
          const Geometry::TargetSet& slag_target = slag_part->template get_target_set<MeshType::shape_dim>();

          // maximum number of layers = number of slices in highest dimensions; this may be less than
          // base_slices[shape_dim], since entire layers might have been removed during the deslagging process
          const Index max_layers = this->_base_slices[shape_dim-1];
          const Index num_elems = slag_target.get_num_entities();

          // compute denominator for layer index computation
          Index denom = 1u;
          for(int j(0); j+1 < shape_dim; ++j)
            denom *= this->_base_slices[Index(j)];

          // compute number of elements for each potential layer
          std::vector<Index> aux(max_layers, 0u);
          for(Index i(0); i < num_elems; ++i)
            ++aux[slag_target[i] / denom];

          // count number of actually used layers
          Index num_layers = 0u;
          for(Index i(0); i < max_layers; ++i)
          {
            Index k = num_layers;
            if(aux[i] > 0u)
              ++num_layers;
            aux[i] = k;
          }

          // compute layering
          Adjacency::Coloring layering(num_elems, num_layers);
          Index* layers = layering.get_coloring();
          for(Index i(0); i < num_elems; ++i)
            layers[i] = aux[slag_target[i] / denom];

          return layering;
        }

        /**
         * \brief Computes the layering for a refined mesh
         *
         * \param[in] slag_node
         * The refined mesh node
         *
         * \returns
         * The layering for the refined and deslagged mesh node
         */
        Adjacency::Coloring _compute_refined_mesh_layering(const MeshNodeType& slag_node, const Adjacency::Coloring& coarse_layering) const
        {
          if(_unstructered_mesh)
          {
            return Adjacency::Coloring();
          }
          // get element target set of slag patch
          const Geometry::MeshPart<MeshType>* slag_part = slag_node.get_patch(-1);
          XASSERT(slag_part != nullptr);
          const Geometry::TargetSet& slag_target = slag_part->template get_target_set<MeshType::shape_dim>();

          // maximum number of layers = twice the number of layers in coarse mesh; this may actually be less,
          // since entire layers might have been removed during the deslagging process
          const Index max_layers = 2 * coarse_layering.get_num_colors();
          const Index num_elems = slag_target.get_num_entities();

          // compute number of elements for each potential layer
          std::vector<Index> aux(max_layers, 0u);
          for(Index i(0); i < num_elems; ++i)
          {
            Index qux = slag_target[i] >> (shape_dim - 1);
            ++aux[(coarse_layering[qux >> 1] << 1) | (qux & 1)]; // all hail to bit-shift magic!
          }

          // count number of actually used layers
          Index num_layers = 0u;
          for(Index i(0); i < max_layers; ++i)
          {
            Index k = num_layers;
            if(aux[i] > 0u)
              ++num_layers;
            aux[i] = k;
          }

          // compute layering
          Adjacency::Coloring layering(num_elems, num_layers);
          Index* layers = layering.get_coloring();
          for(Index i(0); i < num_elems; ++i)
          {
            Index qux = slag_target[i] >> (shape_dim - 1);
            layers[i] = aux[(coarse_layering[qux >> 1] << 1) | (qux & 1)];
          }

          return layering;
        }

        /**
         * \brief Extracts a patch  layering from the parent mesh layering
         *
         * \param[in] parent_node
         * The parent node that the patch was extracted from
         *
         * \param[in] parent_layering
         * The layering of the parent node
         *
         * \param[in] child_rank
         * The rank of the extracted patch
         *
         * \returns
         * The layering for the refined deslagged mesh node
         */
        Adjacency::Coloring _extract_patch_layering(const MeshNodeType& slag_node, const Adjacency::Coloring& parent_layering, int child_rank) const
        {
          if(_unstructered_mesh)
            return Adjacency::Coloring();
          // get element target set of patch
          const Geometry::MeshPart<MeshType>* patch_part = slag_node.get_patch(child_rank);
          XASSERT(patch_part != nullptr);
          const Geometry::TargetSet& patch_target = patch_part->template get_target_set<MeshType::shape_dim>();

          // maximum number of layers = number of layers in parent layering; this may actually be less,
          // since entire layers might have been removed during the deslagging process
          const Index max_layers = parent_layering.get_num_colors();
          const Index num_elems = patch_target.get_num_entities();

          // compute number of elements for each potential layer
          std::vector<Index> aux(max_layers, 0u);
          for(Index i(0); i < num_elems; ++i)
            ++aux[parent_layering[patch_target[i]]];

          // count number of actually used layers
          Index num_layers = 0u;
          for(Index i(0); i < max_layers; ++i)
          {
            Index k = num_layers;
            if(aux[i] > 0u)
              ++num_layers;
            aux[i] = k;
          }

          // compute layering
          Adjacency::Coloring layering(num_elems, num_layers);
          Index* layers = layering.get_coloring();
          for(Index i(0); i < num_elems; ++i)
            layers[i] = aux[parent_layering[patch_target[i]]];

          return layering;
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

          // save slagged base mesh node
          std::unique_ptr<MeshNodeType> base_slag_node = std::move(base_mesh_node);
          // deslag base mesh node
          base_mesh_node = this->_deslag_mesh_node(*base_slag_node);

          // refine and deslag up to desired minimum level
          int lvl = 0;
          for(; lvl < ancestor.desired_level_min; ++lvl)
          {
            base_slag_node = base_mesh_node->refine_unique(this->_adapt_mode);
            base_mesh_node = this->_deslag_mesh_node(*base_slag_node);
          }

          // compute coloring for the coarse level
          Adjacency::Coloring base_coloring;
          if(_unstructered_mesh)
          {
            base_coloring = this->_compute_unstructered_coloring(*base_mesh_node);
          }
          else
          {
            if(lvl == 0)
              base_coloring = this->_compute_base_mesh_coloring(*base_slag_node);
            else
              base_coloring = this->_compute_refined_mesh_coloring(*base_slag_node);
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
            std::shared_ptr<LevelType> level_ptr = std::make_shared<LevelType>(lvl, std::move(base_mesh_node),
              std::move(base_slag_node), std::move(base_coloring));
            this->push_level_front(0, level_ptr);

            // continue with refined node
            base_slag_node = std::move(refined_node);
            base_mesh_node = this->_deslag_mesh_node(*base_slag_node);
            base_coloring = _unstructered_mesh ? this->_compute_unstructered_coloring(*base_mesh_node) : this->_compute_refined_mesh_coloring(*base_slag_node);
          }

          // save chosen maximum level
          this->_chosen_levels.push_front(std::make_pair(lvl, 1));

          // push finest level
          this->push_level_front(0, std::make_shared<LevelType>(lvl, std::move(base_mesh_node),
            std::move(base_slag_node), std::move(base_coloring)));
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

          XASSERTM(!this->_keep_base_levels, "VoxelDomainControl cannot keep base levels!");

          // save slagged base mesh node
          std::unique_ptr<MeshNodeType> base_slag_node = std::move(base_mesh_node);

          // deslag base mesh node
          base_mesh_node = this->_deslag_mesh_node(*base_slag_node);
          Adjacency::Coloring base_coloring;

          // refine up to chosen partition level
          int lvl = 0;
          for(; lvl <= ancestor.desired_level_max; ++lvl)
          {
            // can we apply a partitioner?
            if((lvl >= ancestor.parti_level) && this->_apply_parti(ancestor, *base_mesh_node))
              break;

            // no partitioning found?
            if(lvl >= ancestor.desired_level_max)
              break;

            // refine and deslag base mesh node
            base_slag_node = base_mesh_node->refine_unique(this->_adapt_mode);
            base_mesh_node = this->_deslag_mesh_node(*base_slag_node);
            base_coloring = _unstructered_mesh ? this->_compute_unstructered_coloring(*base_mesh_node) : this->_compute_refined_mesh_coloring(*base_slag_node);
          }

          // no valid partitioning found?
          XASSERTM(ancestor.parti_found, "VoxelDomainControl failed to find a valid partitioning");

          // set the selected partitioner level
          ancestor.parti_level = lvl;

          // compute coloring if we're still on level 0; otherwise the coloring has been computed in the loop above
          if(lvl == 0)
          {
            base_coloring = _unstructered_mesh ? this->_compute_unstructered_coloring(*base_mesh_node) : this->_compute_refined_mesh_coloring(*base_slag_node);
          }

          // extract our patch
          std::vector<int> neighbor_ranks;
          std::unique_ptr<MeshNodeType> patch_mesh_node(
            base_mesh_node->extract_patch(neighbor_ranks, ancestor.parti_graph, this->_comm.rank()));

          // extract our coloring (only if we need it on this level)
          Adjacency::Coloring patch_coloring;
          if(lvl == ancestor.desired_level_min)
            patch_coloring = this->_extract_patch_coloring(*base_mesh_node, base_coloring, this->_comm.rank());

          // set the neighbor ranks of our child layer
          layer->set_neighbor_ranks(neighbor_ranks);

          // create an empty patch slag node
          std::unique_ptr<MeshNodeType> patch_slag_node;

          // refine up to minimum level
          for(; lvl < ancestor.desired_level_min; ++lvl)
          {
            // refine the patch mesh
            patch_slag_node = patch_mesh_node->refine_unique(this->_adapt_mode);
            patch_mesh_node = this->_deslag_mesh_node(*patch_slag_node, ancestor, false);
          }

          // compute coloring if we don't have one yet
          if(patch_coloring.empty())
            patch_coloring = _unstructered_mesh ? this->_compute_unstructered_coloring(*patch_mesh_node) : this->_compute_refined_mesh_coloring(*patch_slag_node);

          // save chosen minimum level
          if(lvl < ancestor.desired_level_max)
            this->_chosen_levels.push_front(std::make_pair(lvl, 0));

          // refine up to maximum level
          for(; lvl < ancestor.desired_level_max; ++lvl)
          {
            // refine the patch mesh
            auto refined_node = patch_mesh_node->refine_unique(this->_adapt_mode);

            // create new level
            std::shared_ptr<LevelType> level_ptr = std::make_shared<LevelType>(lvl, std::move(patch_mesh_node),
              std::move(patch_slag_node), std::move(patch_coloring));

            // push this (unrefined) level
            this->push_level_front(0, level_ptr);

            // continue with refined node
            patch_slag_node = std::move(refined_node);
            patch_mesh_node = this->_deslag_mesh_node(*patch_slag_node, ancestor, false);
            patch_coloring = _unstructered_mesh ? this->_compute_unstructered_coloring(*patch_mesh_node) : this->_compute_refined_mesh_coloring(*patch_slag_node);
          }

          // save chosen maximum level
          this->_chosen_levels.push_front(std::make_pair(lvl, this->_comm.size()));

          // push finest level
          this->push_level_front(0, std::make_shared<LevelType>(lvl, std::move(patch_mesh_node),
            std::move(patch_slag_node), std::move(patch_coloring)));
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

          XASSERTM(!this->_keep_base_levels, "VoxelDomainControl cannot keep base levels!");

          // we start counting at level 0
          int lvl = 0;

          // deslag the base-mesh node and move pointer to a new parent-mesh pointer;
          // the base-mesh is always the one on the layer with only 1 process
          std::unique_ptr<MeshNodeType> parent_slag_node;
          std::unique_ptr<MeshNodeType> parent_mesh_node = this->_deslag_mesh_node(*base_mesh_node);

          // create the coloring and layering for the unrefined base mesh
          Adjacency::Coloring parent_coloring = _unstructered_mesh ? this->_compute_unstructered_coloring(*parent_mesh_node) : this->_compute_base_mesh_coloring(*base_mesh_node);

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
            // note: contrary to the standard parti domain control, we calculate our partition always to the level
            // we split to the subdomains, i.e. on the level we actually change our layer
            // due to this, we enforce, that we at least win 1 real level to the parent
            ancestor.parti_level += lvl;
            ancestor.parti_level = Math::max(ancestor.parti_level, ancestor.desired_level_min);
            ancestor.parti_level = Math::max(ancestor.parti_level,parent_min_lvl + int(!is_base_layer));

            // Note: each progeny group within the main communicator may have chosen a different
            // partitioning level at this point. We will compensate this by adjusting the minimum
            // refinement level of the child layer after the partitioning step below.

            // refine up to the chosen partitioning level
            //for(; lvl < ancestor.parti_level; ++lvl)
            for(; lvl <= ancestor.desired_level_max; ++lvl)
            {
              // can we apply a partitioner?
              if(lvl >= ancestor.parti_level)
              {
                // try to apply the partitioner
                this->_apply_parti(ancestor, *parent_mesh_node);

                // partitioning successful?
                int parti_ok = ancestor.parti_found ? 1 : 0;

                // check if all processes found a valid partitioning
                this->_comm.allreduce(&parti_ok, &parti_ok, std::size_t(1), Dist::op_min);
                if(parti_ok > 0)
                  break;

                // nope, at least one process did not receive a partition, so keep refining
                ancestor.parti_found = false;
              }

              // no partitioning found?
              if(lvl >= ancestor.desired_level_max)
                break;

              // refine and deslag the parent mesh node
              auto refined_node = parent_mesh_node->refine_unique(this->_adapt_mode);

              // a level pointer for the unrefined level; we need this to access the parent layering
              std::shared_ptr<LevelType> level_ptr;

              // push the base mesh into our parent layer if desired
              if((ancestor.layer_p >= 0) && (parent_min_lvl >= 0) && (lvl >= parent_min_lvl))
              {
                level_ptr = std::make_shared<LevelType>(lvl, std::move(parent_mesh_node),
                  std::move(parent_slag_node), std::move(parent_coloring));
                this->push_level_front(ancestor.layer_p, level_ptr);
              }

              // continue with refined node
              parent_slag_node = std::move(refined_node);
              parent_mesh_node = this->_deslag_mesh_node(*parent_slag_node, ancestor, !is_base_layer);
              parent_coloring = _unstructered_mesh ? this->_compute_unstructered_coloring(*parent_mesh_node) : this->_compute_refined_mesh_coloring(*parent_slag_node);
            }

            // no valid partitioning found?
            XASSERTM(ancestor.parti_found, "VoxelDomainControl failed to find a valid partitioning");

            // set the selected partitioner level
            ancestor.parti_level = lvl;

            // extract our patch
            std::vector<int> neighbor_ranks;
            std::unique_ptr<MeshNodeType> patch_mesh_node(
              parent_mesh_node->extract_patch(neighbor_ranks, ancestor.parti_graph, ancestor.progeny_child));

            // create an empty patch slag node (which we might need, if we refine further for this partition)
            std::unique_ptr<MeshNodeType> patch_slag_node;

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
            if(ancestor.layer_p >= 0)
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
            XASSERTM(lvl == global_level_min, "INTERNAL ERROR");
            // for(; lvl < global_level_min; ++lvl)
            // {
            //   // refine the base mesh node
            //   std::unique_ptr<MeshNodeType> coarse_slag_node(std::move(parent_slag_node));
            //   std::unique_ptr<MeshNodeType> coarse_mesh_node(std::move(parent_mesh_node));
            //   parent_slag_node = coarse_mesh_node->refine_unique(this->_adapt_mode);
            //   // parent_mesh_node = this->_deslag_mesh_node(*parent_slag_node);
            //   parent_mesh_node = this->_deslag_mesh_node(*parent_slag_node, ancestor, !is_base_layer);

            //   // push base mesh to parent layer if desired
            //   if((ancestor.layer_p >= 0) && (parent_min_lvl >= 0) && (lvl >= parent_min_lvl))
            //   {
            //     // clear patches before pushing this node as they are redundant here
            //     coarse_mesh_node->clear_patches();
            //     this->push_level_front(ancestor.layer_p, std::make_shared<LevelType>(lvl, std::move(coarse_mesh_node), std::move(coarse_slag_node), std::move(parent_coloring)));
            //   }
            //   parent_coloring = _unstructered_mesh ? this->_compute_unstructered_coloring(*parent_mesh_node) : this->_compute_refined_mesh_coloring(*parent_slag_node);

            //   // refine the patch mesh
            //   patch_slag_node = patch_mesh_node->refine_unique(this->_adapt_mode);
            //   patch_mesh_node = this->_deslag_mesh_node(*patch_slag_node, ancestor, !is_base_layer);
            //   // patch_mesh_node = this->_deslag_mesh_node(*patch_slag_node);
            //   // patch_mesh_node = this->_deslag_mesh_node(*patch_mesh_node->refine_unique(this->_adapt_mode));
            // }

            // extract our coloring
            Adjacency::Coloring patch_coloring = this->_extract_patch_coloring(*parent_mesh_node, parent_coloring, ancestor.progeny_child);

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

            // a level pointer for the unrefined level; we need this to access the parent layering
            std::shared_ptr<LevelType> level_ptr;

            // push the finest base-mesh
            if(ancestor.layer_p >= 0)
            {
              this->push_level_front(ancestor.layer_p, std::make_shared<LevelType>(lvl, std::move(parent_mesh_node),
                std::move(parent_slag_node), std::move(parent_coloring)));
            }

            // continue with the next layer
            parent_slag_node = std::move(patch_slag_node);
            parent_mesh_node = std::move(patch_mesh_node);
            parent_coloring = std::move(patch_coloring);
          }

          // get the desired maximum level
          // int desired_level_max = Math::max(this->_ancestry.front().desired_level_max, lvl+1);
          int desired_level_max = this->_ancestry.front().desired_level_max;
          XASSERTM(desired_level_max > lvl, "Trying to refine larger than provided level hirarchy");

          for(; lvl < desired_level_max; ++lvl)
          {
            // refine the patch mesh
            auto refined_node = parent_mesh_node->refine_unique(this->_adapt_mode);

            // a level pointer for the unrefined level; we need this to access the parent layering
            std::shared_ptr<LevelType> level_ptr = std::make_shared<LevelType>(lvl, std::move(parent_mesh_node),
              std::move(parent_slag_node), std::move(parent_coloring));


            // push patch mesh to this level
            this->push_level_front(0, level_ptr);

            // continue with refined node
            parent_slag_node = std::move(refined_node);
            parent_mesh_node = this->_deslag_mesh_node(*parent_slag_node, this->_ancestry.front(), false);
            parent_coloring = _unstructered_mesh ? this->_compute_unstructered_coloring(*parent_mesh_node) : this->_compute_refined_mesh_coloring(*parent_slag_node);
          }

          // set chosen maximum level for finest layer
          this->_chosen_levels.push_front(std::make_pair(lvl, this->_comm.size()));

          // push finest level
          this->push_level_front(0, std::make_shared<LevelType>(lvl, std::move(parent_mesh_node),
            std::move(parent_slag_node), std::move(parent_coloring)));
        }

#endif // defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
      }; // class CGALDomainControl<...>
    } // namespace Domain
  } // namespace Control
} // namespace FEAT
#endif
