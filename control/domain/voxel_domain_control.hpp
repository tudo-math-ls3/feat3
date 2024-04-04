// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2024 by Stefan Turek & the FEAT group
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

#include <control/domain/parti_domain_control_base.hpp>

namespace FEAT
{
  namespace Control
  {
    namespace Domain
    {
      /**
       * \brief Interface for slag masker for the VoxelDomainControl
       *
       * \author Peter Zajac
       */
      template<typename CoordType_, int dim_>
      class VoxelSlagMasker
      {
      public:
        /// the coordinate type
        typedef CoordType_ CoordType;
        /// the point type
        typedef Tiny::Vector<CoordType_, dim_> PointType;

        /// virtual destructor
        virtual ~VoxelSlagMasker()
        {
        }

        /**
         * \brief Computes the mask for an entire X-coordinate row
         *
         * \param[inout] mask
         * The mask vector for the current X-row. Is allocated to correct size, but its contents are undefined upon entry.
         * Assuming that the length of this vector is n, then the i-th entry of the mask vector is to be set to 1 if
         * the point with the X-coordinate equal to x_min + i/(n-1)*x_max is inside the mask, otherwise it is to be set to 0.
         *
         * \param[in] x_min, x_max
         * The minimal and maximal X-coordinate of the X-row for which the mask is to be computed.
         *
         * \param[in] other_coords
         * The point that contains the Y- and (in 3D) Z-coordinates of the X-row that for which the mask is to be computed.
         * The first coordinate, i.e. point[0], is undefined, so only point[1] and (in 3D) point[2] are set to the Y-/Z-coords.
         */
        virtual void mask_row(std::vector<int>& mask, const CoordType x_min, const CoordType x_max, const PointType& other_coords) = 0;
      }; // class VoxelSlagMasker<...>

#ifdef FEAT_HAVE_CGAL
      /**
       * \brief Voxel slag masker implementation for CGAL wrappers
       *
       * This class implements the voxel slag masker interface by using the I/O test functionality of a CGALWrapper object.
       *
       * \author Peter Zajac
       */
      template<typename CoordType_>
      class VoxelSlagCGALMasker :
        public VoxelSlagMasker<CoordType_, 3>
      {
      public:
        typedef VoxelSlagMasker<CoordType_, 3> BaseClass;
        using typename BaseClass::CoordType;
        using typename BaseClass::PointType;

      private:
        const Geometry::CGALWrapper<CoordType_>& _cgal_wrapper;
        const bool _invert;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] cgal_wrapper
         * A \resident reference to the CGAL wrapper object that is to be used for the inside/outside test.
         *
         * \param[in] invert
         * Specifies whether inside and outside are to be swapped.
         */
        explicit VoxelSlagCGALMasker(const Geometry::CGALWrapper<CoordType_>& cgal_wrapper, bool invert) :
          _cgal_wrapper(cgal_wrapper),
          _invert(invert)
        {
        }

        virtual void mask_row(std::vector<int>& mask, const CoordType x_min, const CoordType x_max, const PointType& other_coords) override
        {
          if(mask.empty())
            return;

          const std::size_t nv = mask.size();
          const CoordType s =  (x_max - x_min) / CoordType(nv - 1u);

          for(std::size_t i(0); i < nv; ++i)
          {
            mask[i] = (this->_cgal_wrapper.point_inside(x_min + CoordType(i) * s, other_coords[1], other_coords[2]) ? int(!_invert) : int(_invert));
          }
        }
      }; // class VoxelSlagCGALMasker<...>
#endif // FEAT_HAVE_CGAL

      /**
       * \brief Voxel slag masker implementation for chart classes
       *
       * This class implements the voxel slag masker interface by using the signed distance functionality of the Atlas Chart implementations
       *
       * \author Peter Zajac
       */
      template<typename MeshType>
      class VoxelSlagChartMasker :
        public VoxelSlagMasker<typename MeshType::CoordType, MeshType::shape_dim>
      {
      public:
        typedef VoxelSlagMasker<typename MeshType::CoordType, MeshType::shape_dim> BaseClass;
        using typename BaseClass::CoordType;
        using typename BaseClass::PointType;

      private:
        const Geometry::Atlas::ChartBase<MeshType>& _chart;
        const bool _invert;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] chart
         * A \resident reference to the chart object that is to be used for the inside/outside test.
         *
         * \param[in] invert
         * Specifies whether inside and outside are to be swapped.
         */
        explicit VoxelSlagChartMasker(const Geometry::Atlas::ChartBase<MeshType>& chart, bool invert) :
          _chart(chart),
          _invert(invert)
        {
        }

        virtual void mask_row(std::vector<int>& mask, const CoordType x_min, const CoordType x_max, const PointType& other_coords) override
        {
          if(mask.empty())
            return;

          const std::size_t nv = mask.size();
          const CoordType s =  (x_max - x_min) / CoordType(nv - 1u);
          PointType point(other_coords);

          for(std::size_t i(0); i < nv; ++i)
          {
            point[0] = x_min + CoordType(i) * s;
            mask[i] = (_chart.signed_dist(point) <= CoordType(0) ? int(!_invert) : int(_invert));
          }
        }
      }; // class VoxelSlagChartMasker<...>

      /**
       * \brief Voxel slag masker implementation for lambda expression tests
       *
       * This class implements the voxel slag masker interface by evaluating a user-supplied lambda expression, which
       * performs the inside/outside test for each point.
       *
       * \author Peter Zajac
       */
      template<typename Lambda_, typename CoordType_, int dim_>
      class VoxelSlagLambdaMasker :
        public VoxelSlagMasker<CoordType_, dim_>
      {
      public:
        typedef VoxelSlagMasker<CoordType_, dim_> BaseClass;
        using typename BaseClass::CoordType;
        using typename BaseClass::PointType;

      public:
        /// the lambda expression
        Lambda_ _lambda;

        /**
         * \brief Constructor
         *
         * \param[in] lambda
         * The lambda expression that is to be used to perform the inside/outside test.
         * The lambda expression is expected to take a const reference to a PointType object as a sole input
         * parameter and to return either true or false depending on whether the point is inside or outside.
         */
        explicit VoxelSlagLambdaMasker(Lambda_&& lambda) :
          _lambda(std::forward<Lambda_>(lambda))
        {
        }

        virtual void mask_row(std::vector<int>& mask, const CoordType x_min, const CoordType x_max, const PointType& other_coords) override
        {
          if(mask.empty())
            return;

          const std::size_t nv = mask.size();
          const CoordType s =  (x_max - x_min) / CoordType(nv - 1u);
          PointType point(other_coords);
          // define a const reference to ensure that the lambda expression does not modify the point
          const PointType& cpoint(point);

          for(std::size_t i(0); i < nv; ++i)
          {
            point[0] = x_min + CoordType(i) * s;
            mask[i] = (_lambda(cpoint) ? 1 : 0);
          }
        }
      }; // class VoxelSlagLambdaMasker<...>

      /**
       * \brief Wrapper class for domain levels for VoxelDomainControl
       *
       * This class acts as a derived class for a domain level class and extends it by an additional "slag level"
       * object, which is required to assemble grid transfer operators (i.e. prolongation/restriction) for spaces
       * organized in the VoxelDomainControl class. Furthermore, this class also stores the optional element coloring
       * and element layering for parallel assembly.
       *
       * \author Peter Zajac
       */
      template<typename DomainLevel_>
      class VoxelDomainLevelWrapper :
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

        /// element layering (stored as coloring)
        Adjacency::Coloring element_layering;

        explicit VoxelDomainLevelWrapper(int _lvl_idx, std::unique_ptr<MeshNodeType> node,
          Adjacency::Coloring&& coloring, Adjacency::Coloring&& layering) :
          BaseClass(_lvl_idx, std::move(node)),
          slag_level(),
          element_coloring(std::forward<Adjacency::Coloring>(coloring)),
          element_layering(std::forward<Adjacency::Coloring>(layering))
        {
        }

        explicit VoxelDomainLevelWrapper(int _lvl_idx, std::unique_ptr<MeshNodeType> node,
          std::unique_ptr<MeshNodeType> slag_node, Adjacency::Coloring&& coloring, Adjacency::Coloring&& layering) :
          BaseClass(_lvl_idx, std::move(node)),
          slag_level(),
          element_coloring(std::forward<Adjacency::Coloring>(coloring)),
          element_layering(std::forward<Adjacency::Coloring>(layering))
        {
          if(slag_node)
            slag_level = std::make_shared<DomainLevel_>(_lvl_idx, std::move(slag_node));
        }
      }; // class VoxelDomainLevelWrapper<...>

      /**
       * \brief Hierarchical partitioned Voxel Domain Control
       *
       * \author Peter Zajac
       */
      template<typename DomainLevel_>
      class VoxelDomainControl :
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

        /// coordinate type
        typedef typename MeshType::CoordType CoordType;

      protected:
        /// creation stage to keep track what has already been initialized
        int _create_stage;

        /// the bounding box of the structured domain
        Tiny::Vector<CoordType, shape_dim> _bbox_min, _bbox_max;
        /// number of element slices on base mesh level 0
        std::array<Index, shape_dim> _base_slices;
        /// number of element slices on slag mask level = fine mesh level
        std::array<Index, shape_dim> _mask_slices;

        /// length of a single X-row in slag mask = number of vertices in X-row divided by 64 and rounded up
        std::size_t _slag_mask_row_length;
        /// total number of X-rows in slag mask = product of number of vertices in Y-row and Z-row
        std::size_t _slag_mask_row_count;
        /// the actual slag mask vector
        std::vector<std::uint64_t> _slag_mask;

        /// the input base-mesh node on level 0
        std::unique_ptr<MeshNodeType> _base_mesh_node;

        /// a bunch of stop-watches
        StopWatch _watch_create_base_mesh, _watch_create_slag_mask, _watch_create_hierarchy;

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
        explicit VoxelDomainControl(const Dist::Comm& comm_, bool support_multi_layered) :
          BaseClass(comm_, support_multi_layered),
          _create_stage(0),
          _slag_mask_row_length(0u),
          _slag_mask_row_count(0u)
        {
        }

        /// virtual destructor
        virtual ~VoxelDomainControl()
        {
        }

        /**
         * \brief Creates a 2D structured rectangle base mesh and returns its base-mesh node
         *
         * This functions returns a reference to the base-mesh node that is created and managed internally
         * by this class. This base-mesh node reference can be used to create mesh-parts at the base-mesh level,
         * which will be refined and deslagged during the creation of the domain hierarchy.
         *
         * \attention The returned base-mesh node reference is only valid until the create() method is called and
         * is orphaned after the create() method returns, because the base-mesh node is then processed internally
         * to create the domain level hierarchy.
         *
         * \param[in] num_x, num_y
         * The number of element slices in X- and Y-direction.
         *
         * \param[in] x_min, x_max
         * The X-range of the cuboid domain.
         *
         * \param[in] y_min, y_max
         * The X-range of the cuboid domain.
         *
         * \returns A reference to the internal base-mesh node containing the base-mesh.
         */
        MeshNodeType& create_base_mesh_2d(Index num_x, Index num_y, CoordType x_min, CoordType x_max, CoordType y_min, CoordType y_max)
        {
          XASSERTM(shape_dim == 2, "invalid dimension: this function can only be called for 2D domains!");
          XASSERTM(_create_stage == 0, "base mesh already created!");

          XASSERT(num_x > Index(0));
          XASSERT(num_y > Index(0));
          XASSERT(x_min < x_max);
          XASSERT(y_min < y_max);

          _watch_create_base_mesh.start();

          _base_slices[0] = num_x;
          _base_slices[1] = num_y;
          _bbox_min[0] = x_min;
          _bbox_max[0] = x_max;
          _bbox_min[1] = y_min;
          _bbox_max[1] = y_max;

          // create a structured mesh from factory
          _base_mesh_node = MeshNodeType::make_unique(
            Geometry::StructUnitCubeFactory<MeshType>::make_unique_from(num_x, num_y), &this->_atlas);

          // get the node's vertex set
          auto& vtx = _base_mesh_node->get_mesh()->get_vertex_set();
          Index n = vtx.get_num_vertices();
          for(Index i(0); i < n; ++i)
          {
            vtx[i][0] = x_min + (x_max - x_min) * vtx[i][0];
            vtx[i][1] = y_min + (y_max - y_min) * vtx[i][1];
          }

          _watch_create_base_mesh.stop();

          // update creation stage
          this->_create_stage = 1;

          // return the base-mesh node
          return *_base_mesh_node.get();
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
          XASSERTM(shape_dim == 3, "invalid dimension: this function can only be called for 3D domains!");
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
         * \brief Creates a voxel slag mask based on a VoxelSlagMasker object
         *
         * This function can only be called after the base-mesh has been created by calling create_base_mesh_2d or
         * create_base_mesh_3d.
         *
         * \attention
         * This function is a "collective" function in distributed parallel (i.e. MPI) simulations and each process
         * must be given the same masker, because this function distributes the actual work over all processes before
         * gathering the computed slag mask on all processes afterwards!
         *
         * \param[in] masker
         * A \transient reference to the voxel slag masker.
         */
        void create_slag_mask(VoxelSlagMasker<CoordType, shape_dim>& masker)
        {
          XASSERTM(this->_create_stage >= 1, "invalid creation stage; create the base-mesh first");
          XASSERTM(this->_create_stage <= 1, "invalid creation stage: slag mask already set");

          this->_watch_create_slag_mask.start();
          this->_compute_slag_mask(masker);
          this->_watch_create_slag_mask.stop();

          // update creation stage
          this->_create_stage = 2;
        }

        /**
         * \brief Creates a voxel slag mask based on a lambda expression
         *
         * This function can only be called after the base-mesh has been created by calling create_base_mesh_2d or
         * create_base_mesh_3d.
         *
         * \attention
         * This function is a "collective" function in distributed parallel (i.e. MPI) simulations and each process
         * must be given the same lambda expression, because this function distributes the actual work over all
         * processes before gathering the computed slag mask on all processes afterwards!
         *
         * \param[in] lambda
         * A lambda expression that takes a const reference (or a copy) of a Tiny::Vector<CoordType_, shape_dim> that
         * represents the point to be tested and returns \c true, if the point is in the interior of the domain or
         * \c false if it is outside the domain.
         */
        template<typename Lambda_>
        void create_slag_mask_from_lambda(Lambda_&& lambda)
        {
          XASSERTM(this->_create_stage >= 1, "invalid creation stage; create the base-mesh first");
          XASSERTM(this->_create_stage <= 1, "invalid creation stage: slag mask already set");

          this->_watch_create_slag_mask.start();
          VoxelSlagLambdaMasker<Lambda_, CoordType, shape_dim> masker(std::forward<Lambda_>(lambda));
          this->_compute_slag_mask(masker);
          this->_watch_create_slag_mask.stop();

          // update creation stage
          this->_create_stage = 2;
        }

        /**
         * \brief Creates a voxel slag mask based on a chart
         *
         * This function can only be called after the base-mesh has been created by calling create_base_mesh_2d or
         * create_base_mesh_3d.
         *
         * \attention
         * This function is a "collective" function in distributed parallel (i.e. MPI) simulations and each process
         * must be given the same chart, because this function distributes the actual work over all
         * processes before gathering the computed slag mask on all processes afterwards!
         *
         * \param[in] chart
         * A \transient reference to the chart that is to be used for the inside outside test.
         *
         * \param[in] invert
         * Specifies, whether the inside/outside is to be inverted.
         */
        void create_slag_mask_from_chart(const Geometry::Atlas::ChartBase<MeshType>& chart, bool invert)
        {
          XASSERTM(this->_create_stage >= 1, "invalid creation stage; create the base-mesh first");
          XASSERTM(this->_create_stage <= 1, "invalid creation stage: slag mask already set");

          this->_watch_create_slag_mask.start();
          VoxelSlagChartMasker<MeshType> masker(chart, invert);
          this->_compute_slag_mask(masker);
          this->_watch_create_slag_mask.stop();

          // update creation stage
          this->_create_stage = 2;
        }

        /**
         * \brief Creates a voxel slag mask based on a surface triangulation stored in an OFF file
         *
         * This function can only be called after the base-mesh has been created by calling create_base_mesh_3d.
         *
         * \note
         * This function can only be used in 3D and if FEAT was compiled and linked against the CGAL library.
         *
         * \attention
         * This function is a "collective" function in distributed parallel (i.e. MPI) simulations, which will
         * read in the file on the root process and distribute the file's contents to all processes.
         *
         * \param[in] filename
         * The filename of the OFF file that is to be read in.
         *
         * \param[in] invert
         * Specifies, whether the inside/outside is to be inverted.
         */
        void create_slag_mask_from_cgal_off(const String& filename, bool invert)
        {
          std::stringstream sstr;
          DistFileIO::read_common(sstr, filename, this->_comm);
          this->create_slag_mask_from_cgal_off(sstr, invert);
        }

        /**
         * \brief Creates a voxel slag mask based on a surface triangulation stored in an OFF file
         *
         * This function can only be called after the base-mesh has been created by calling create_base_mesh_3d.
         *
         * \note
         * This function can only be used in 3D and if FEAT was compiled and linked against the CGAL library.
         *
         * \attention
         * This function is a "collective" function in distributed parallel (i.e. MPI) simulations and each process
         * must be given the same OFF file content stream, because this function distributes the actual work over all
         * processes before gathering the computed slag mask on all processes afterwards!
         *
         * \param[in] is
         * The input stream that contains the OFF data; this should ideally be a std::stringstring that was read in
         * by using the DistFileIO::read_common function to avoid a DDoS attack on the file system.
         *
         * \param[in] invert
         * Specifies, whether the inside/outside is to be inverted.
         */
        void create_slag_mask_from_cgal_off(std::istream& is, bool invert)
        {
          // this only makes sense in 3D...
          XASSERTM(shape_dim == 3, "CGAL OFF slag mask creation is only available in 3D");
#ifdef FEAT_HAVE_CGAL
          Geometry::CGALWrapper<CoordType> cgal_wrapper(is, Geometry::CGALFileMode::fm_off);
          VoxelSlagCGALMasker<CoordType> slag_masker(cgal_wrapper, invert);
          this->create_slag_mask(slag_masker);
#else // no FEAT_HAVE_CGAL
          XABORTM("FEAT is not build and linked against CGAL library!");
          (void)is;
          (void)invert;
#endif // FEAT_HAVE_CGAL
        }

        /**
         * \brief Creates the domain level hierarchy
         *
         * This function can only be called after the base-mesh has been created by calling create_base_mesh_2d or
         * create_base_mesh_3d and after the slag mask has been created by calling create_slag_mask().
         */
        void create_hierarchy()
        {
          XASSERTM(this->_create_stage >= 2, "invalid creation stage; create the base-mesh and slag mask first");
          XASSERTM(this->_create_stage <= 2, "invalid creation stage: domain hierarchy already created");

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
            this->_watch_create_slag_mask.elapsed() + this->_watch_create_hierarchy.elapsed();

          // update creation stage
          this->_create_stage = 3;
          this->_was_created = true;
        }

        /// Returns a const reference to the StopWatch that measures the base-mesh creation phase.
        const StopWatch& get_watch_base_mesh() const
        {
          return this->_watch_create_base_mesh;
        }

        /// Returns a const reference to the StopWatch that measures the slag-mask creation phase.
        const StopWatch& get_watch_slag_mask() const
        {
          return this->_watch_create_slag_mask;
        }

        /// Returns a const reference to the StopWatch that measures the base-mesh creation phase.
        const StopWatch& get_watch_hierarchy() const
        {
          return this->_watch_create_hierarchy;
        }

        /**
         * \brief Returns a const reference to the internal slag mask vector
         */
        const std::vector<std::uint64_t>& get_slag_mask_vector() const
        {
          return this->_slag_mask;
        }

        /**
         * \brief Gathers the vertex slag mask on a given domain level.
         *
         * \param[in] level
         * The physical domain level for which the vertex mask is to be computed.
         *
         * \returns A vector of ints representing the vertex slag mask on the desired level. The length of the vector
         * corresponds to the number of vertices on the patch on the desired domain level.
         */
        std::vector<int> gather_vertex_slag_mask(Index level) const
        {
          XASSERTM(this->_create_stage == 3, "domain level hierarchy has to be created before the vertex mask can be computed");
          XASSERT(level < this->size_physical());

          const auto& vtx = this->at(level)->get_mesh().get_vertex_set();
          const Index nv = vtx.get_num_vertices();
          std::vector<int> mask(nv, 0);

          for(Index i(0); i < nv; ++i)
            mask[i] = _check_slag_mask(vtx[i]);

          return mask;
        }

        /**
         * \brief Gathers the element slag weights on a given domain level.
         *
         * The "slag weight" of an element is defined as the number of slag mask points with a value of 1,
         * which lie inside the element cuboid, divided by the total number of slag mask points, which lie inside
         * the element cuboid.
         *
         * \param[in] level
         * The physical domain level for which the element weights are to be gathered.
         *
         * \returns A vector of floats in range [0,1] representing the element slag weights on the desired level.
         * The length of the vector corresponds to the number of elements on the patch on the desired domain level.
         */
        std::vector<WeightType> gather_element_slag_weights(Index level) const
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

      protected:
        /**
         * \brief Compresses a slag mask row into the slag mask map.
         *
         * \param[in] mask
         * The slag mask row as computed by a VoxelSlagMasker implementation
         *
         * \param[in] row
         * The index of the row that is to be compressed.
         */
        void _compress_slag_mask_row(const std::vector<int>& mask, const Index row)
        {
          ASSERTM(row < this->_slag_mask_row_count, "invalid drop mask row");
          const Index n = Index(mask.size());
          const Index off = row * this->_slag_mask_row_length;
          for(Index i(0); i < n; ++i)
          {
            _slag_mask[off + (i >> 6)] |= std::uint64_t(mask[i] != 0) << (i & 0x3F);
          }
        }

        /**
         * \brief Maps a coordinate to a X-/Y-/Z-index in the slag mask
         *
         * \param[in] xyz
         * The X-/Y-/Z-coordinate that is to be mapped
         *
         * \param[in] dim
         * Specifies which dimension is to be mapped: 0=X, 1=Y, 2=Z
         */
        Index _map_slag_idx(CoordType xyz, int dim) const
        {
          CoordType t = (xyz - _bbox_min[dim]) / (_bbox_max[dim] - _bbox_min[dim]);
          int idx = int(std::nearbyint(CoordType(_mask_slices[std::size_t(dim)]) * t));
          ASSERTM(idx >= 0, "tested vertex coordinate is not inside coarse bounding box");
          ASSERTM(idx <= int(this->_mask_slices[Index(dim)]), "tested vertex coordinate is not inside coarse bounding box");
          return Index(idx);
        }

        /**
         * \brief Checks the slag mask entry for a 3D vertex
         *
         * \param[in] vtx
         * The coordinates of the vertex whose entry in the slag mask is to be checked.
         *
         * \returns
         * \c true, if the vertex is inside the slag mask domain (i.e. the slag mask entry is 1),
         * or \c false, if the vertex is outside of the slag mask domain.
         */
        bool _check_slag_mask(const Tiny::Vector<CoordType, 3>& vtx) const
        {
          Index xidx = _map_slag_idx(vtx[0], 0);
          Index yidx = _map_slag_idx(vtx[1], 1);
          Index zidx = _map_slag_idx(vtx[2], 2);
          Index row = zidx * (this->_mask_slices[1] + 1u) + yidx;
          return (this->_slag_mask[row * this->_slag_mask_row_length + (xidx >> 6)] >> (xidx & 0x3F)) & 0x1;
        }

        /**
         * \brief Checks the slag mask entry for a 2D vertex
         *
         * \param[in] vtx
         * The coordinates of the vertex whose entry in the slag mask is to be checked.
         *
         * \returns
         * \c true, if the vertex is inside the slag mask domain (i.e. the slag mask entry is 1),
         * or \c false, if the vertex is outside of the slag mask domain.
         */
        bool _check_slag_mask(const Tiny::Vector<CoordType, 2>& vtx) const
        {
          Index xidx = _map_slag_idx(vtx[0], 0);
          Index yidx = _map_slag_idx(vtx[1], 1);
          return (this->_slag_mask[yidx * this->_slag_mask_row_length + (xidx >> 6)] >> (xidx & 0x3F)) & 0x1;
        }

        /**
         * \brief Checks the slag mask entry for a 1D vertex
         *
         * \param[in] vtx
         * The coordinates of the vertex whose entry in the slag mask is to be checked.
         *
         * \returns
         * \c true, if the vertex is inside the slag mask domain (i.e. the slag mask entry is 1),
         * or \c false, if the vertex is outside of the slag mask domain.
         */
        bool _check_slag_mask(const Tiny::Vector<CoordType, 1>& vtx) const
        {
          Index xidx = _map_slag_idx(vtx[0], 0);
          return (this->_slag_mask[(xidx >> 6)] >> (xidx & 0x3F)) & 0x1;
        }

        /**
         * \brief Computes a range of slag mask rows for a 3D domain
         *
         * \param[in] masker
         * A \transient reference to the masker that is to be used to compute the mask
         *
         * \param[in] slices
         * A reference to the _mask_slices member object; this is only used to choose the correct dimensional overload of this function.
         *
         * \param[in] beg
         * The index of the first slag mask row that is to be computed.
         *
         * \param[in] end
         * The index of the last slag mask row that is to be computed plus one.
         */
        void _compute_slag_mask_rows(VoxelSlagMasker<CoordType, shape_dim>& masker,
          const std::array<Index,3>& DOXY(slices), Index beg, Index end)
        {
          std::vector<int> row_mask(this->_mask_slices[0] + 1u, 0);
          const CoordType sy = (_bbox_max[1] - _bbox_min[1]) / CoordType(this->_mask_slices[1]);
          const CoordType sz = (_bbox_max[2] - _bbox_min[2]) / CoordType(this->_mask_slices[2]);
          Tiny::Vector<CoordType, 3> coords;

          for(Index i(beg); i < end; ++i)
          {
            // row = iz * (this->_mask_slices[1] + 1u) + iy
            coords[1] = _bbox_min[1] + CoordType(i % (this->_mask_slices[1] + 1u)) * sy;
            coords[2] = _bbox_min[2] + CoordType(i / (this->_mask_slices[1] + 1u)) * sz;
            masker.mask_row(row_mask, _bbox_min[0], _bbox_max[0], coords);
            this->_compress_slag_mask_row(row_mask, i);
          }
        }

        /**
         * \brief Computes a range of slag mask rows for a 2D domain
         *
         * \param[in] masker
         * A \transient reference to the masker that is to be used to compute the mask
         *
         * \param[in] slices
         * A reference to the _mask_slices member object; this is only used to choose the correct dimensional overload of this function.
         *
         * \param[in] beg
         * The index of the first slag mask row that is to be computed.
         *
         * \param[in] end
         * The index of the last slag mask row that is to be computed plus one.
         */
        void _compute_slag_mask_rows(VoxelSlagMasker<CoordType, shape_dim>& masker,
          const std::array<Index,2>& DOXY(slices), Index beg, Index end)
        {
          std::vector<int> row_mask(this->_mask_slices[0] + 1u, 0);
          const CoordType sy = (_bbox_max[1] - _bbox_min[1]) / CoordType(this->_mask_slices[1]);
          Tiny::Vector<CoordType, 2> coords;
          for(Index i(beg); i < end; ++i)
          {
            coords[1] = _bbox_min[1] + CoordType(i) * sy;
            masker.mask_row(row_mask, _bbox_min[0], _bbox_max[0], coords);
            this->_compress_slag_mask_row(row_mask, i);
          }
        }

        /**
         * \brief Computes a range of slag mask rows for a 1D domain
         *
         * \param[in] masker
         * A \transient reference to the masker that is to be used to compute the mask
         *
         * \param[in] slices
         * A reference to the _mask_slices member object; this is only used to choose the correct dimensional overload of this function.
         *
         * \param[in] beg
         * The index of the first slag mask row that is to be computed. Must be equal to 0 in 1D.
         *
         * \param[in] end
         * The index of the last slag mask row that is to be computed plus one. Must be equal to 1 in 1D.
         */
        void _compute_slag_mask_rows(VoxelSlagMasker<CoordType, shape_dim>& masker,
          const std::array<Index,1>& DOXY(slices), Index DOXY(beg), Index DOXY(end))
        {
          std::vector<int> row_mask(this->_mask_slices[0] + 1u, 0);
          Tiny::Vector<CoordType, 1> coords;
          masker.mask_row(row_mask, _bbox_min[0], _bbox_max[0], coords);
          this->_compress_slag_mask_row(row_mask, 0);
        }

        /**
         * \brief Gathers the slag mask weight for a 3D bounding box cuboid
         *
         * \param[in] box_min, box_max
         * The minimum and maximum X/Y/Z coordinates that make up the bounding box
         *
         * \returns The slag mask weight for the bounding box.
         */
        WeightType _gather_slag_mask_weight(const std::array<Index, 3>& box_min, const std::array<Index, 3>& box_max) const
        {
          Index count = 0u;
          for(Index zidx(box_min[2]); zidx <= box_max[2]; ++zidx)
          {
            for(Index yidx(box_min[1]); yidx <= box_max[1]; ++yidx)
            {
              Index row = zidx * (this->_mask_slices[1] + 1u) + yidx;
              for(Index xidx(box_min[0]); xidx <= box_max[0]; ++xidx)
                count += (this->_slag_mask[row * this->_slag_mask_row_length + (xidx >> 6)] >> (xidx & 0x3F)) & 0x1;
            }
          }
          return WeightType(count) / WeightType((box_max[2] - box_min[2] + 1u) * (box_max[1] - box_min[1] + 1u) * (box_max[0] - box_min[0] + 1u));
        }

        /**
         * \brief Gathers the slag mask weight for a 2D bounding box rectangle
         *
         * \param[in] box_min, box_max
         * The minimum and maximum X/Y coordinates that make up the bounding box
         *
         * \returns The slag mask weight for the bounding box.
         */
        WeightType _gather_slag_mask_weight(const std::array<Index, 2>& box_min, const std::array<Index, 2>& box_max) const
        {
          Index count = 0u;
          for(Index yidx(box_min[1]); yidx <= box_max[1]; ++yidx)
          {
            for(Index xidx(box_min[0]); xidx <= box_max[0]; ++xidx)
              count += (this->_slag_mask[yidx * this->_slag_mask_row_length + (xidx >> 6)] >> (xidx & 0x3F)) & 0x1;
          }
          return WeightType(count) / WeightType((box_max[1] - box_min[1] + 1u) * (box_max[0] - box_min[0] + 1u));
        }

        /**
         * \brief Gathers the slag mask weight for a 1D bounding box interval
         *
         * \param[in] box_min, box_max
         * The minimum and maximum X coordinates that make up the bounding box
         *
         * \returns The slag mask weight for the bounding box.
         */
        WeightType _gather_slag_mask_weight(const std::array<Index, 1>& box_min, const std::array<Index, 1>& box_max) const
        {
          Index count = 0u;
          for(Index xidx(box_min[0]); xidx <= box_max[0]; ++xidx)
            count += (this->_slag_mask[(xidx >> 6)] >> (xidx & 0x3F)) & 0x1;
          return WeightType(count) / WeightType(box_max[0] - box_min[0] + 1u);
        }

        /**
         * \brief Computes the slag mask for the domain
         *
         * In an MPI-parallel simulation, this function splits the computation of the mask across all processes in the
         * communicator and gathers the mask on all processes afterwards. In the borderline case where there are more
         * processes in the communicator than there are slag mask rows in the mask, this function creates a temporary
         * sub-communicator with as many processes as there are slag mask rows to distribute the work over these
         * processes and performs a broadcast on the original communicator after gathering the mask on rank 0.
         *
         * \param[in] masker
         * A \transient reference to the masker that is to be used to compute the mask. Must be the same on all
         * MPI processes
         */
        void _compute_slag_mask(VoxelSlagMasker<CoordType, shape_dim>& masker)
        {
          const int level_max = this->_desired_levels.front().first;

          // compute number of fine mesh bounding box slices in each dimension and total number of vertices
          for(Index i(0); i < Index(shape_dim); ++i)
           this->_mask_slices[i] = this->_base_slices[i] << level_max;

          // compute number of fine mesh bounding box vertices
          this->_slag_mask_row_length = this->_mask_slices[0] + 1u;
          this->_slag_mask_row_count = 1u;
          for(Index i(1); i < Index(shape_dim); ++i)
            this->_slag_mask_row_count *= (this->_mask_slices[i] + 1u);

          // translate row length from bits to 64-bit integers, i.e. divide by 64 and round up if necessary
          this->_slag_mask_row_length = (this->_slag_mask_row_length / 64u) + (this->_slag_mask_row_length % 64u != 0);

          // single process or MPI parallel?
          if(this->_comm.size() <= 1)
          {
            // allocate slag mask and format to zero
            this->_slag_mask.resize(this->_slag_mask_row_length * this->_slag_mask_row_count, 0u);

            // compute all slag mask rows
            this->_compute_slag_mask_rows(masker, this->_mask_slices, 0u, this->_slag_mask_row_count);
          }
          else if(Index(this->_comm.size()) <= this->_slag_mask_row_count)
          {
            // there are at least as many rows in the slag mask as there are MPI processes; this is the usual case
            // we can use the entire communicator (usually MPI_COMM_WORLD) for the slag mask computation
            Index comm_rank = Index(this->_comm.rank());
            Index comm_size = Index(this->_comm.size());

            // compute row count for a single process and round up; this may be bigger than _slag_mask_row_count,
            // so the last process may have less rows to chew through than all other processes
            Index row_count = (this->_slag_mask_row_count + comm_size - 1u) / comm_size;
            std::size_t block_len = this->_slag_mask_row_length * row_count;

            // allocate slag mask and format to zero
            this->_slag_mask.resize(block_len * Index(comm_size), 0u);

            // compute first and last+1 row of this process; make sure the last process doesn't go beyond the total count
            Index row_beg = comm_rank * row_count;
            Index row_end = Math::min((comm_rank + 1u) * row_count, this->_slag_mask_row_count);

            // debug output
            //this->_comm.allprint(String(">>>") + stringify(row_beg).pad_front(6) + ":" + stringify(row_count) + ">" + stringify(row_end).pad_front(6));

            // compute rows for this process
            this->_compute_slag_mask_rows(masker, this->_mask_slices, row_beg, row_end);

            // perform an in-place allgather to synchronize the slag mask over all processes
            this->_comm.allgather(this->_slag_mask.data(), block_len, this->_slag_mask.data(), block_len);
          }
          else if(this->_slag_mask_row_count < Index(this->_comm.size()))
          {
            // more processes than slag rows; this is a borderline scenario, but we have to support it anyways
            // create a sub-communicator and then compute mask on that sub-communicator
            Dist::Comm sub_comm = this->_comm.comm_create_range_incl(int(this->_slag_mask_row_count));

            // does this process participate in the sub-communicator?
            if(!sub_comm.is_null())
            {
              // get sub-communicator rank and size
              Index comm_rank = Index(sub_comm.rank());
              Index comm_size = Index(sub_comm.size());

              // compute row count for a single process and round up; this may be bigger than _slag_mask_row_count,
              // so the last process may have less rows to chew through than all other processes
              Index row_count = (this->_slag_mask_row_count + comm_size - 1u) / comm_size;
              std::size_t block_len = this->_slag_mask_row_length * row_count;

              // allocate slag mask and format to zero
              this->_slag_mask.resize(block_len * Index(comm_size), 0u);

              // compute first and last+1 row of this process; make sure the last process doesn't go beyond the total count
              Index row_beg = comm_rank * row_count;
              Index row_end = Math::min((comm_rank + 1u) * row_count, this->_slag_mask_row_count);

              // debug output
              //this->_comm.allprint(String(">>>") + stringify(row_beg).pad_front(6) + ":" + stringify(row_count) + ">" + stringify(row_end).pad_front(6));

              // compute rows for this process
              this->_compute_slag_mask_rows(masker, this->_mask_slices, row_beg, row_end);

              // perform an in-place gather to synchronize the slag mask on rank 0
              sub_comm.gather(this->_slag_mask.data(), block_len, this->_slag_mask.data(), block_len, 0);
            }
            else
            {
              // just allocate the slag mask vector for the following broadcast
              this->_slag_mask.resize(this->_slag_mask_row_length * this->_slag_mask_row_count, 0u);
            }

            // broadcast from rank 0 to all processes our entire domain communicator
            this->_comm.bcast(this->_slag_mask.data(), this->_slag_mask_row_length * this->_slag_mask_row_count, 0);
          }

          // that's it
        }

        /**
         * \brief Gather a vector of element indices that intersect the masked domain
         *
         * \attention
         * It is silently assumed that all elements in the mesh are axis-parallel rectangles/cuboids; this function
         * may yield incorrect results if the mesh violates this assumption!
         *
         * \param[in] mesh
         * A \transient reference to the mesh whose elements are to be tested against the slag mask
         *
         * \returns
         * A vector containing the indices of all mesh elements that intersect the domain that is represented by
         * the slag mask.
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

          // element bounding box coordinates
          Tiny::Vector<CoordType, shape_dim> elbox_min, elbox_max;

          // mapped element bounding box slag mask indices
          std::array<Index, shape_dim> slag_idx_min, slag_idx_max;

          // loop over all elements
          for(Index ielem(0); ielem < num_elems; ++ielem)
          {
            // compute bounding box of element
            // note: actually, we should always have
            //   elbox_min = vtx[verts_at_elem(ielem, 0)]
            //   elbox_max = vtx[verts_at_elem(ielem, 1 << shape_dim - 1)]
            // because of the element orientation of the base mesh, but we'll compute the bounding
            // box for each element manually anyways, just to be sure in case that a derived class
            // might give us a mesh object that is not directly refined from the base mesh using
            // the standard 2-level refinement...
            elbox_min = elbox_max = vtx[verts_at_elem(ielem, 0)];
            for(int j(1); j < verts_at_elem.num_indices; ++j)
            {
              const auto& v = vtx[verts_at_elem(ielem, j)];
              for(int k(0); k < shape_dim; ++k)
                Math::minimax(v[k], elbox_min[k], elbox_max[k]);
            }

            // map slag mask indices of bounding box in drop row indices
            for(int k(0); k < shape_dim; ++k)
            {
              slag_idx_min[Index(k)] = this->_map_slag_idx(elbox_min[k], k);
              slag_idx_max[Index(k)] = this->_map_slag_idx(elbox_max[k], k);
            }

            // gather element weight based on the bounding box; any non-negative weight means that
            // there was at least one masked domain point in the slag mask covered by the element
            if(this->_gather_slag_mask_weight(slag_idx_min, slag_idx_max) > 1E-12)
              masked_elems.push_back(ielem);
          }

          // return the list of all elements that we found
          return masked_elems;
        }

        /**
         * \brief Gathers the element slag weights for a given mesh
         *
         * The "slag weight" of an element is defined as the number of slag mask points with a value of 1,
         * which lie inside the element cuboid, divided by the total number of slag mask points, which lie inside
         * the element cuboid.
         *
         * \param[in] mesh
         * The mesh for which the element weights are to be gathered.
         *
         * \returns A vector of floats in range [0,1] representing the element slag weights on the desired level.
         * The length of the vector corresponds to the number of elements on the patch on the desired domain level.
         */
        std::vector<WeightType> _gather_element_weights(const MeshType& mesh) const
        {
          const Index num_elems = mesh.get_num_elements();

          // allocate the weight vector
          std::vector<WeightType> weights(num_elems, 0.0);

          // get the vertex set and the vertices-at-element index set
          const auto& vtx = mesh.get_vertex_set();
          const auto& verts_at_elem = mesh.template get_index_set<shape_dim, 0>();

          // element bounding box coordinates
          Tiny::Vector<CoordType, shape_dim> elbox_min, elbox_max;

          // mapped element bounding box slag mask indices
          std::array<Index, shape_dim> slag_idx_min, slag_idx_max;

          // loop over all elements
          for(Index ielem(0); ielem < num_elems; ++ielem)
          {
            // compute bounding box of element
            // note: actually, we should always have
            //   elbox_min = vtx[verts_at_elem(ielem, 0)]
            //   elbox_max = vtx[verts_at_elem(ielem, 1 << shape_dim - 1)]
            // because of the element orientation of the base mesh, but we'll compute the bounding
            // box for each element manually anyways, just to be sure in case that a derived class
            // might give us a mesh object that is not directly refined from the base mesh using
            // the standard 2-level refinement...
            elbox_min = elbox_max = vtx[verts_at_elem(ielem, 0)];
            for(int j(1); j < verts_at_elem.num_indices; ++j)
            {
              const auto& v = vtx[verts_at_elem(ielem, j)];
              for(int k(0); k < shape_dim; ++k)
                Math::minimax(v[k], elbox_min[k], elbox_max[k]);
            }

            // map slag mask indices of bounding box in drop row indices
            for(int k(0); k < shape_dim; ++k)
            {
              slag_idx_min[Index(k)] = this->_map_slag_idx(elbox_min[k], k);
              slag_idx_max[Index(k)] = this->_map_slag_idx(elbox_max[k], k);
            }

            // gather element weight based on the bounding box
            weights[ielem] = this->_gather_slag_mask_weight(slag_idx_min, slag_idx_max);
          }

          // returns the element weights
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
          return mesh_node.extract_patch(_gather_masked_elements(*mesh_node.get_mesh()), true, false, false);
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
          std::unique_ptr<MeshNodeType> new_node = mesh_node.extract_patch(_gather_masked_elements(*mesh_node.get_mesh()), true, false, false);
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
          const Index num_colors = Index(1) << MeshType::shape_dim;
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

          // compute base-mesh layering
          Adjacency::Coloring base_layering = this->_compute_base_mesh_layering(*base_slag_node);

          // refine and deslag up to desired minimum level
          int lvl = 0;
          for(; lvl < ancestor.desired_level_min; ++lvl)
          {
            base_slag_node = base_mesh_node->refine_unique(this->_adapt_mode);
            base_mesh_node = this->_deslag_mesh_node(*base_slag_node);
            base_layering = this->_compute_refined_mesh_layering(*base_slag_node, base_layering);
          }

          // compute coloring for the coarse level
          Adjacency::Coloring base_coloring;
          if(lvl == 0)
            base_coloring = this->_compute_base_mesh_coloring(*base_slag_node);
          else
            base_coloring = this->_compute_refined_mesh_coloring(*base_slag_node);

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
              std::move(base_slag_node), std::move(base_coloring), std::move(base_layering));
            this->push_level_front(0, level_ptr);

            // continue with refined node
            base_slag_node = std::move(refined_node);
            base_mesh_node = this->_deslag_mesh_node(*base_slag_node);
            base_coloring = this->_compute_refined_mesh_coloring(*base_slag_node);
            base_layering = this->_compute_refined_mesh_layering(*base_slag_node, level_ptr->element_layering);
          }

          // save chosen maximum level
          this->_chosen_levels.push_front(std::make_pair(lvl, 1));

          // push finest level
          this->push_level_front(0, std::make_shared<LevelType>(lvl, std::move(base_mesh_node),
            std::move(base_slag_node), std::move(base_coloring), std::move(base_layering)));
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
          Adjacency::Coloring base_layering = this->_compute_base_mesh_layering(*base_slag_node);

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
            base_coloring = this->_compute_refined_mesh_coloring(*base_slag_node);
            base_layering = this->_compute_refined_mesh_layering(*base_slag_node, base_layering);
          }

          // no valid partitioning found?
          XASSERTM(ancestor.parti_found, "VoxelDomainControl failed to find a valid partitioning");

          // set the selected partitioner level
          ancestor.parti_level = lvl;

          // compute coloring if we're still on level 0; otherwise the coloring has been computed in the loop above
          if(lvl == 0)
            base_coloring = this->_compute_base_mesh_coloring(*base_slag_node);

          // extract our patch
          std::vector<int> neighbor_ranks;
          std::unique_ptr<MeshNodeType> patch_mesh_node(
            base_mesh_node->extract_patch(neighbor_ranks, ancestor.parti_graph, this->_comm.rank()));

          // extract our coloring (only if we need it on this level)
          Adjacency::Coloring patch_coloring;
          if(lvl == ancestor.desired_level_min)
            patch_coloring = this->_extract_patch_coloring(*base_mesh_node, base_coloring, this->_comm.rank());

          // extract layering
          Adjacency::Coloring patch_layering = this->_extract_patch_layering(*base_mesh_node, base_layering, this->_comm.rank());

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
            patch_layering = this->_compute_refined_mesh_layering(*patch_slag_node, patch_layering);
          }

          // compute coloring if we don't have one yet
          if(patch_coloring.empty())
            patch_coloring = this->_compute_refined_mesh_coloring(*patch_slag_node);

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
              std::move(patch_slag_node), std::move(patch_coloring), std::move(patch_layering));

            // push this (unrefined) level
            this->push_level_front(0, level_ptr);

            // continue with refined node
            patch_slag_node = std::move(refined_node);
            patch_mesh_node = this->_deslag_mesh_node(*patch_slag_node, ancestor, false);
            patch_coloring = this->_compute_refined_mesh_coloring(*patch_slag_node);
            patch_layering = this->_compute_refined_mesh_layering(*patch_slag_node, level_ptr->element_layering);
          }

          // save chosen maximum level
          this->_chosen_levels.push_front(std::make_pair(lvl, this->_comm.size()));

          // push finest level
          this->push_level_front(0, std::make_shared<LevelType>(lvl, std::move(patch_mesh_node),
            std::move(patch_slag_node), std::move(patch_coloring), std::move(patch_layering)));
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
          Adjacency::Coloring parent_coloring = this->_compute_base_mesh_coloring(*base_mesh_node);
          Adjacency::Coloring parent_layering = this->_compute_base_mesh_layering(*base_mesh_node);

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
            ancestor.parti_level = Math::max(ancestor.parti_level, ancestor.desired_level_min);

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
                  std::move(parent_slag_node), std::move(parent_coloring), std::move(parent_layering));
                this->push_level_front(ancestor.layer_p, level_ptr);
              }

              // continue with refined node
              parent_slag_node = std::move(refined_node);
              parent_mesh_node = this->_deslag_mesh_node(*parent_slag_node, ancestor, !is_base_layer);
              parent_coloring = this->_compute_refined_mesh_coloring(*parent_slag_node);
              parent_layering = this->_compute_refined_mesh_layering(*parent_slag_node,
                level_ptr ? level_ptr->element_layering : parent_layering);
            }

            // no valid partitioning found?
            XASSERTM(ancestor.parti_found, "VoxelDomainControl failed to find a valid partitioning");

            // set the selected partitioner level
            ancestor.parti_level = lvl;

            // extract our patch
            std::vector<int> neighbor_ranks;
            std::unique_ptr<MeshNodeType> patch_mesh_node(
              parent_mesh_node->extract_patch(neighbor_ranks, ancestor.parti_graph, ancestor.progeny_child));

            // extract our coloring
            Adjacency::Coloring patch_coloring = this->_extract_patch_coloring(*parent_mesh_node, parent_coloring, ancestor.progeny_child);

            // extract our layering
            Adjacency::Coloring patch_layering = this->_extract_patch_layering(*parent_mesh_node, parent_layering, ancestor.progeny_child);

            // create an empty patch slag node
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
            /*for(; lvl < global_level_min; ++lvl)
            {
              // refine the base mesh node
              std::unique_ptr<MeshNodeType> coarse_slag_node(std::move(parent_slag_node));
              std::unique_ptr<MeshNodeType> coarse_mesh_node(std::move(parent_mesh_node));
              parent_slag_node = coarse_mesh_node->refine_unique(this->_adapt_mode);
              parent_mesh_node = this->_deslag_mesh_node(*parent_slag_node);

              // push base mesh to parent layer if desired
              if((ancestor.layer_p >= 0) && (parent_min_lvl >= 0) && (lvl >= parent_min_lvl))
              {
                // clear patches before pushing this node as they are redundant here
                //coarse_mesh_node->clear_patches();
                this->push_level_front(ancestor.layer_p, std::make_shared<LevelType>(lvl, std::move(coarse_mesh_node), std::move(coarse_slag_node)));
              }

              // refine the patch mesh
              //patch_slag_node = patch_mesh_node->refine_unique(this->_adapt_mode);
              //patch_mesh_node = this->_deslag_mesh_node(*patch_slag_node);
              patch_mesh_node = this->_deslag_mesh_node(*patch_mesh_node->refine_unique(this->_adapt_mode));
            }*/

            // split the halos of our base-mesh and compute the halos of our patches from that
            //this->_split_basemesh_halos(ancestor, *parent_slag_node, *patch_slag_node, neighbor_ranks);
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
                std::move(parent_slag_node), std::move(parent_coloring), std::move(parent_layering)));
            }

            // continue with the next layer
            parent_slag_node = std::move(patch_slag_node);
            parent_mesh_node = std::move(patch_mesh_node);
            parent_coloring = std::move(patch_coloring);
            parent_layering = std::move(patch_layering);
          }

          // get the desired maximum level
          // if we have more than one layer, make sure that the finest one contains at
          // least one level, as otherwise the finest global level would be a ghost level
          int desired_level_max = Math::max(this->_ancestry.front().desired_level_max, lvl+1);

          for(; lvl < desired_level_max; ++lvl)
          {
            // refine the patch mesh
            auto refined_node = parent_mesh_node->refine_unique(this->_adapt_mode);

            // a level pointer for the unrefined level; we need this to access the parent layering
            std::shared_ptr<LevelType> level_ptr = std::make_shared<LevelType>(lvl, std::move(parent_mesh_node),
              std::move(parent_slag_node), std::move(parent_coloring), std::move(parent_layering));


            // push patch mesh to this level
            this->push_level_front(0, level_ptr);

            // continue with refined node
            parent_slag_node = std::move(refined_node);
            parent_mesh_node = this->_deslag_mesh_node(*parent_slag_node, this->_ancestry.front(), false);
            parent_coloring = this->_compute_refined_mesh_coloring(*parent_slag_node);
            parent_layering = this->_compute_refined_mesh_layering(*parent_slag_node, level_ptr->element_layering);
          }

          // set chosen maximum level for finest layer
          this->_chosen_levels.push_front(std::make_pair(lvl, this->_comm.size()));

          // push finest level
          this->push_level_front(0, std::make_shared<LevelType>(lvl, std::move(parent_mesh_node),
            std::move(parent_slag_node), std::move(parent_coloring), std::move(parent_layering)));
        }

#endif // defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
      }; // class VoxelDomainControl<...>
    } // namespace Domain
  } // namespace Control
} // namespace FEAT
