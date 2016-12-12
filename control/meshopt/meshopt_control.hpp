#ifndef FEAT_CONTROL_MESHOPT_MESHOPT_CONTROL_HPP
#define FEAT_CONTROL_MESHOPT_MESHOPT_CONTROL_HPP 1
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/assembly/grid_transfer.hpp>
#include <kernel/assembly/mirror_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/slip_filter_assembler.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/sparse_matrix_bwrappedcsr.hpp>
#include <kernel/lafem/transfer.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/geometry/mesh_file_writer.hpp>
#include <kernel/global/filter.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/muxer.hpp>
#include <kernel/global/nonlinear_functional.hpp>
#include <kernel/global/transfer.hpp>
#include <kernel/global/vector.hpp>

#include <control/domain/domain_control.hpp>

namespace FEAT
{
  namespace Control
  {
    /**
     * \brief Control layer for mesh optimisation
     */
    namespace Meshopt
    {
      /// \cond internal
      /// Forward declarations
      template<typename>
      class MeshoptAssemblerLevel;

      template<typename, typename, typename, template<typename, typename, typename> class>
      class QuadraticSystemLevel;

      template<typename, typename, typename, template<typename, typename, typename> class, template<typename, typename, typename> class>
      class MeshoptSystemLevel;

      /// \endcond

      /**
       * \brief Base class for Meshopt control objects
       *
       * \tparam DomainControl_
       * Type of the underlying domain control
       *
       * \tparam Trafo_
       * The type of transformation defining all meshes, i.e. P1/Q1 transformation
       *
       * For every mesh quality functional, there is a control class which inherits from this. This class only knows
       * the base minimum.
       *
       */
      template<typename DomainControl_, typename Trafo_>
      class MeshoptControlBase
      {
        public:
          /// The type of the underlying domain control
          typedef DomainControl_ DomainControl;
          /// The world dimension, i.e. number of coordinates
          static constexpr int world_dim = DomainControl_::MeshType::world_dim;
          /// Floating point type for coordinates
          typedef typename DomainControl_::MeshType::CoordType CoordType;
          /// Type for buffer vectors exchanging information between mesh and some vector type on which the solvers
          /// operate
          typedef LAFEM::DenseVectorBlocked<Mem::Main, CoordType, Index, world_dim> LocalCoordsBuffer;
          /// corresponding vector mirror
          //typedef LAFEM::VectorMirrorBlocked<Mem::Main, CoordType, Index, world_dim> CoordsMirror;
          typedef LAFEM::VectorMirror<Mem::Main, CoordType, Index> CoordsMirror;
          /// The global version of LocalCoordsBuffer, needed for prepare setting the internal state variable
          typedef Global::Vector<LocalCoordsBuffer, CoordsMirror> GlobalCoordsBuffer;
          /// Type of the vtk exporter this (and derived classes) can write to
          typedef Geometry::ExportVTK<typename Trafo_::MeshType> VTKExporterType;

          /// Check if DomainControl and Trafo are compatible
          static_assert( std::is_same<typename DomainControl_::MeshType, typename Trafo_::MeshType>::value,
          "DomainControl/Trafo MeshType mismatch");

        protected:
          /// The domain control whose mesh objects can be modified
          DomainControl& _dom_ctrl;

        public:

          /**
           * \brief The simplest constructor possible
           */
          explicit MeshoptControlBase(DomainControl& dom_ctrl_)
            : _dom_ctrl(dom_ctrl_)
            {
            }

          /**
           * \brief Empty virtual destructor
           */
          virtual ~MeshoptControlBase()
          {
          }

          /**
           * \brief Computes the cell size quality indicator
           *
           * This can be different for every MeshQualityFunctional.
           *
           * \returns The cell size quality indicator.
           */
          virtual CoordType compute_cell_size_defect(CoordType& lambda_min, CoordType& lambda_max,
          CoordType& vol_min, CoordType& vol_max, CoordType& vol) const = 0;

          /**
           * \brief Optimises the current mesh
           */
          virtual void optimise() = 0;

          /**
           * \brief Sets the internal state variable
           *
           * \param[in] vec_state
           * Global buffer vector containing (possibly new) mesh vertex coordinates
           *
           */
          virtual void prepare(const GlobalCoordsBuffer& DOXY(vec_state)) = 0;

          /**
           * \brief Copies the mesh's vertex coordinates to the buffer vector
           */
          virtual void mesh_to_buffer() = 0;

          /**
           * \brief Copies the contents of the buffer vector to the mesh's vertex coordinates
           */
          virtual void buffer_to_mesh() = 0;

          /**
           * \brief Gets the coordinates buffer vector
           *
           * \returns A reference to the coordinates buffer
           */
          virtual GlobalCoordsBuffer& get_coords() = 0;

          /**
           * \brief Gets the names of all Dirichlet boundaries
           *
           * Note that each name refers to a boundary, but that boundary does not necessarily have to be present due
           * to partitioning etc.
           *
           * \returns A deque of Strings with all Dirichlet boundary names
           */
          virtual std::deque<String> get_dirichlet_boundaries() const = 0;

          /**
           * \brief Gets the names of all slip boundaries
           *
           * Note that each name refers to a boundary, but that boundary does not necessarily have to be present due
           * to partitioning etc.
           *
           * \returns A deque of Strings with all slip boundary names
           */
          virtual std::deque<String> get_slip_boundaries() const = 0;

          /**
           * \brief Get the number of levels in this object
           *
           * \returns The number of levels.
           */
          virtual size_t get_num_levels() const = 0;

          /**
           * \brief Returns a descriptive String
           *
           * \returns The class name as String
           */
          virtual String name() const = 0;

          /**
           * \brief Adds quantities of the underlying mesh quality functional to a given exporter object
           */
          virtual void add_to_vtk_exporter
            (Geometry::ExportVTK<typename Trafo_::MeshType>& DOXY(exporter), const int DOXY(lvl_index)) const
            {
            }

          /**
           * \brief Prints settings of the control object
           */
          virtual void print() const = 0;


      }; // class MeshoptControlBase

      /**
       * \brief (Non)linear system of equations on one mesh refinement level
       *
       * \tparam Mem_
       * Memory architecture for the system of equations wrt. the solver
       *
       * \tparam DT_
       * Floating point type
       *
       * \tparam IT_
       * Index type
       *
       * \tparam Op_
       * The (patch-) local (non-)linear mesh quality functional
       *
       * \tparam GlobalOp_
       * The global wrapper class around Op_
       *
       * Examples for Op_ are Meshopt::DuDvFunctional or MeshOpt::Hyperelasticityfunctional. If the mesh quality
       * functional is quadratic, its gradient gives the linear system of equations to solve, so GlobalOp_ is
       * Global::Matrix. If the mesh quality functional is not quadratic, the gradient is nonlinear and GlobalOp_
       * has to be Global::NonlinearFunctional.
       *
       */
      template
      <
        typename Mem_, typename DT_, typename IT_,
        template<typename, typename, typename> class Functional_
      >
      class QuadraticSystemLevel
      {
        public:
          /// Memory architecture for the solver
          typedef Mem_ MemType;
          /// Floating point precision for the solver
          typedef DT_ DataType;
          /// Index type for the solver
          typedef IT_ IndexType;

          /// (Patch-) Local mesh quality functional type
          typedef Functional_<Mem_, DT_, IT_> LocalFunctional;
          typedef typename LocalFunctional::MatrixType LocalMatrix;

          /// Local left-vectors (dual space)
          typedef typename LocalFunctional::VectorTypeL LocalSystemVectorL;
          /// Local right-vectors (primal space)
          typedef typename LocalFunctional::VectorTypeR LocalSystemVectorR;
          /// Local vectors of scalar quantities
          typedef typename LocalFunctional::ScalarVectorType LocalScalarVector;
          /// Local coordinates buffer type for passing information to or from the mesh
          typedef typename LocalFunctional::CoordsBufferType LocalCoordsBuffer;
          /// Local inter-level transfer matrix
          typedef LAFEM::SparseMatrixBWrappedCSR<Mem_, DT_, IT_, LocalFunctional::MeshType::world_dim>
            LocalSystemTransferMatrix;
          /// Local transfer operator
          typedef LAFEM::Transfer<LocalSystemTransferMatrix> LocalSystemTransfer;

          /// Filter for the local system
          typedef typename LocalFunctional::FilterType LocalSystemFilter;
          /// This is comprised of a sequence of SlipFilters...
          typedef typename LocalFunctional::SlipFilterSequence LocalSlipFilterSequence;
          /// ... and a sequence of UnitFilters
          typedef typename LocalFunctional::DirichletFilterSequence LocalDirichletFilterSequence;

          /// Mirrors for system vectors
          //typedef LAFEM::VectorMirrorBlocked<Mem_, DT_, IT_, LocalFunctional::BlockHeight> SystemMirror;
          typedef LAFEM::VectorMirror<Mem_, DT_, IT_> SystemMirror;
          /// Gates for the system
          typedef Global::Gate<LocalSystemVectorR, SystemMirror> SystemGate;
          /// Mirrors for scalar vectors
          typedef LAFEM::VectorMirror<Mem_, DT_, IT_> ScalarMirror;
          /// Gates for scalar vectors
          typedef Global::Gate<LocalScalarVector, ScalarMirror> ScalarGate;

          /// Global mesh quality functional type
          typedef Global::Matrix<LocalMatrix, SystemMirror, SystemMirror> GlobalSystemMatrix;
          /// Global system filter type
          typedef Global::Filter<LocalSystemFilter, SystemMirror> GlobalSystemFilter;
          /// Global scalar vector type
          typedef Global::Vector<LocalScalarVector, ScalarMirror> GlobalScalarVector;
          /// Global left-vectors
          typedef Global::Vector<LocalSystemVectorL, SystemMirror> GlobalSystemVectorL;
          /// Global right-vectors
          typedef Global::Vector<LocalSystemVectorR, SystemMirror> GlobalSystemVectorR;
          /// Global coordinates buffer
          typedef Global::Vector<LocalCoordsBuffer, SystemMirror> GlobalCoordsBuffer;
          /// Global system transfer operator
          typedef Global::Transfer<LocalSystemTransfer, SystemMirror> GlobalSystemTransfer;
          /// The global muxer for mapping data from one partitioning to the other on this level
          typedef Global::Muxer<LocalSystemVectorR, SystemMirror> GlobalSystemMuxer;

        public:
          LocalFunctional local_functional;
          /// The scalar gate
          ScalarGate gate_scalar;
          /// The system gate
          SystemGate gate_sys;
          /// The global filter
          GlobalSystemFilter filter_sys;
          /// The global system matrix
          GlobalSystemMatrix matrix_sys;
          /// This contains a shallow copy of the operator's coords_buffer
          GlobalCoordsBuffer coords_buffer;
          /// The global partition muxer on this level
          GlobalSystemMuxer coarse_muxer_sys;
          /// The global transfer operator from this level to the next coarser one
          GlobalSystemTransfer transfer_sys;

        private:
          /// This is the number of refines it took from the mesh at the file level to here
          const int _level_index;

        public:
          /**
           * \brief Variadic template constructor
           *
           * \param[in] dirichlet_list
           * List of the names of all Dirichlet boundaries
           *
           * \param[in] slip_list
           * List of the names of all slip boundaries
           */
          template<typename... Args_>
          explicit QuadraticSystemLevel(const int level_index,
          const std::deque<String>& dirichlet_list,
          const std::deque<String>& slip_list,
          Args_&&... args) :
            local_functional(std::forward<Args_>(args)...),
            gate_sys(),
            filter_sys(),
            matrix_sys(&gate_sys, &gate_sys, local_functional.matrix_sys.clone(LAFEM::CloneMode::Shallow)),
            coords_buffer(&gate_sys, local_functional.get_coords().clone(LAFEM::CloneMode::Shallow)),
            coarse_muxer_sys(),
            transfer_sys(&coarse_muxer_sys),
            _level_index(level_index)
            {

              LocalSlipFilterSequence slip_sequence(slip_list);
              LocalDirichletFilterSequence dirichlet_sequence(dirichlet_list);

              LocalSystemFilter local_filter (std::move(slip_sequence), std::move(dirichlet_sequence));

              filter_sys.local() = std::move(local_filter);

            }

          /**
           * \brief Empty virtual destructor
           */
          virtual ~QuadraticSystemLevel()
          {
          }

          /**
           * \brief Returns the level index
           *
           * The level index is the number of times the underlying mesh was refined since reading its construction.
           */
          int get_level_index()
          {
            return _level_index;
          }


      }; // class QuadraticSystemLevel<...>

      /**
       * \brief (Non)linear system of equations on one mesh refinement level
       *
       * \tparam Mem_
       * Memory architecture for the system of equations wrt. the solver
       *
       * \tparam DT_
       * Floating point type
       *
       * \tparam IT_
       * Index type
       *
       * \tparam Op_
       * The (patch-) local (non-)linear mesh quality functional
       *
       * \tparam GlobalOp_
       * The global wrapper class around Op_
       *
       * Examples for Op_ are Meshopt::DuDvFunctional or MeshOpt::Hyperelasticityfunctional. If the mesh quality
       * functional is quadratic, its gradient gives the linear system of equations to solve, so GlobalOp_ is
       * Global::Matrix. If the mesh quality functional is not quadratic, the gradient is nonlinear and GlobalOp_
       * has to be Global::NonlinearFunctional.
       *
       */
      template
      <
        typename Mem_, typename DT_, typename IT_,
        template<typename, typename, typename> class Functional_,
        template<typename, typename, typename> class GlobalOp_
        >
        class MeshoptSystemLevel
        {
          public:
            /// Memory architecture for the solver
            typedef Mem_ MemType;
            /// Floating point precision for the solver
            typedef DT_ DataType;
            /// Index type for the solver
            typedef IT_ IndexType;

            /// (Patch-) Local mesh quality functional type
            typedef Functional_<Mem_, DT_, IT_> LocalFunctional;

            /// Local left-vectors (dual space)
            typedef typename LocalFunctional::VectorTypeL LocalSystemVectorL;
            /// Local right-vectors (primal space)
            typedef typename LocalFunctional::VectorTypeR LocalSystemVectorR;
            /// Local vectors of scalar quantities
            typedef typename LocalFunctional::ScalarVectorType LocalScalarVector;
            /// Local coordinates buffer type for passing information to or from the mesh
            typedef typename LocalFunctional::CoordsBufferType LocalCoordsBuffer;
            /// Local inter-level transfer matrix
            typedef LAFEM::SparseMatrixBWrappedCSR<Mem_, DT_, IT_, LocalFunctional::MeshType::world_dim>
              LocalSystemTransferMatrix;
            /// Local transfer operator
            typedef LAFEM::Transfer<LocalSystemTransferMatrix> LocalSystemTransfer;

            /// Filter for the local system
            typedef typename LocalFunctional::FilterType LocalSystemFilter;
            /// This is comprised of a sequence of SlipFilters...
            typedef typename LocalFunctional::SlipFilterSequence LocalSlipFilterSequence;
            /// ... and a sequence of UnitFilters
            typedef typename LocalFunctional::DirichletFilterSequence LocalDirichletFilterSequence;

            /// Mirrors for system vectors
            //typedef LAFEM::VectorMirrorBlocked<Mem_, DT_, IT_, LocalFunctional::BlockHeight> SystemMirror;
            typedef LAFEM::VectorMirror<Mem_, DT_, IT_> SystemMirror;
            /// Gates for the system
            typedef Global::Gate<LocalSystemVectorR, SystemMirror> SystemGate;
            /// Mirrors for scalar vectors
            typedef LAFEM::VectorMirror<Mem_, DT_, IT_> ScalarMirror;
            /// Gates for scalar vectors
            typedef Global::Gate<LocalScalarVector, ScalarMirror> ScalarGate;

            /// Global mesh quality functional type
            typedef GlobalOp_<LocalFunctional, SystemMirror, SystemMirror> GlobalFunctional;
            /// Global system filter type
            typedef Global::Filter<LocalSystemFilter, SystemMirror> GlobalSystemFilter;
            /// Global scalar vector type
            typedef Global::Vector<LocalScalarVector, ScalarMirror> GlobalScalarVector;
            /// Global left-vectors
            typedef Global::Vector<LocalSystemVectorL, SystemMirror> GlobalSystemVectorL;
            /// Global right-vectors
            typedef Global::Vector<LocalSystemVectorR, SystemMirror> GlobalSystemVectorR;
            /// Global coordinates buffer
            typedef Global::Vector<LocalCoordsBuffer, SystemMirror> GlobalCoordsBuffer;
            /// Global system transfer operator
            typedef Global::Transfer<LocalSystemTransfer, SystemMirror> GlobalSystemTransfer;
            /// The global muxer for mapping data from one partitioning to the other on this level
            typedef Global::Muxer<LocalSystemVectorR, SystemMirror> GlobalSystemMuxer;

            typedef GlobalFunctional GlobalSystemMatrix;

          public:
            /// The scalar gate
            ScalarGate gate_scalar;
            /// The system gate
            SystemGate gate_sys;
            /// The global filter
            GlobalSystemFilter filter_sys;
            /// The global system matrix
            GlobalFunctional op_sys;
            GlobalSystemMatrix& matrix_sys;
            LocalFunctional& local_functional;
            /// This contains a shallow copy of the operator's coords_buffer
            GlobalCoordsBuffer coords_buffer;
            /// The global partition muxer on this level
            GlobalSystemMuxer coarse_muxer_sys;
            /// The global transfer operator from this level to the next coarser one
            GlobalSystemTransfer transfer_sys;

          private:
            /// This is the number of refines it took from the mesh at the file level to here
            const int _level_index;

          public:
            /**
             * \brief Variadic template constructor
             *
             * \param[in] dirichlet_list
             * List of the names of all Dirichlet boundaries
             *
             * \param[in] slip_list
             * List of the names of all slip boundaries
             */
            template<typename... Args_>
            explicit MeshoptSystemLevel(const int level_index,
            const std::deque<String>& dirichlet_list,
            const std::deque<String>& slip_list,
            Args_&&... args) :
              gate_sys(),
              filter_sys(),
              op_sys(&gate_sys, &gate_sys, std::forward<Args_>(args)...),
              matrix_sys(op_sys),
              local_functional(op_sys.local()),
              coords_buffer(&gate_sys, op_sys.local().get_coords().clone(LAFEM::CloneMode::Shallow)),
              coarse_muxer_sys(),
              transfer_sys(&coarse_muxer_sys),
              _level_index(level_index)
              {

                LocalSlipFilterSequence slip_sequence(slip_list);
                LocalDirichletFilterSequence dirichlet_sequence(dirichlet_list);

                LocalSystemFilter local_filter (std::move(slip_sequence), std::move(dirichlet_sequence));

                filter_sys.local() = std::move(local_filter);

              }

            /**
             * \brief Empty virtual destructor
             */
            virtual ~MeshoptSystemLevel()
            {
            }

            /**
             * \brief Returns the level index
             *
             * The level index is the number of times the underlying mesh was refined since reading its construction.
             */
            int get_level_index()
            {
              return _level_index;
            }


        }; // class MeshoptSystemLevel<...>

      /**
       * \brief Base class for assembler levels for mesh optimisation
       *
       * \tparam Space_
       * The finite element space the problem is solved on.
       *
       * Control objects inheriting from MeshControlBase may overwrite this, but it contains all basic
       * functionality.
       *
       */
      template<typename Space_>
      class MeshoptAssemblerLevel
      {
        public:
          /// The finite element space
          typedef Space_ SpaceType;
          /// The underlying transformation
          typedef typename SpaceType::TrafoType TrafoType;
          /// The type of mesh we use
          typedef typename TrafoType::MeshType MeshType;
          /// Type for one level of the DomainControl, needed i.e. for multigrid solvers
          typedef Control::Domain::DomainLevel<MeshType> DomainLevelType;
          /// Type for one layer of the DomainControl, needed i.e. for ScaRC
          typedef Control::Domain::DomainLayer<MeshType> DomainLayerType;

        public:
          /// The domain level holding the RootMeshNode this refers to
          DomainLevelType& domain_level;
          /// The mesh the transformation and our FE spaces live on
          MeshType& mesh;
          /// The transformation
          TrafoType trafo;
          /// The FE space we solve our problem on
          SpaceType trafo_space;
          /// Cubature factory for integration
          Cubature::DynamicFactory cubature;

          /// All UnitFilterAssemblers for mesh optimisation
          std::map<String, std::shared_ptr<Assembly::UnitFilterAssembler<MeshType>>> dirichlet_asm;
          /// All SlipFilterAssemblers for mesh optimisation
          std::map<String, std::shared_ptr<Assembly::SlipFilterAssembler<MeshType>>> slip_asm;

        public:
          /**
           * \brief
           */
          explicit MeshoptAssemblerLevel(DomainLevelType& dom_lvl, const std::deque<String>& dirichlet_list,
          const std::deque<String>& slip_list) :
            domain_level(dom_lvl),
            mesh(domain_level.get_mesh()),
            trafo(mesh),
            trafo_space(trafo),
            cubature("auto-degree:" + stringify(Math::sqr(SpaceType::local_degree)+2)),
            dirichlet_asm(),
            slip_asm()
            {
              // For every MeshPart specified in the list of Dirichlet boundaries, create a UnitFilterAssembler and
              // insert it into the std::map
              for(auto& it : dirichlet_list)
              {
                // Create empty assembler
                auto new_asm = std::make_shared<Assembly::UnitFilterAssembler<MeshType>>();

                // Add the MeshPart to the assembler if it is there. There are legimate reasons for it NOT to be
                // there, i.e. we are in parallel and our patch is not adjacent to that MeshPart
                auto* mpp = domain_level.get_mesh_node()->find_mesh_part(it);
                if(mpp != nullptr)
                {
                  new_asm->add_mesh_part(*mpp);
                }

                // Insert into the map
                String identifier(it);
                dirichlet_asm.emplace(identifier, new_asm);
              }

              // For every MeshPart specified in the list of slip boundaries, create a SlipFilterAssembler and
              // insert it into the std::map
              for(auto& it : slip_list)
              {
                // Create empty assembler
                auto new_asm = std::make_shared<Assembly::SlipFilterAssembler<MeshType>>(mesh);

                // Add the MeshPart to the assembler if it is there. There are legimate reasons for it NOT to be
                // there, i.e. we are in parallel and our patch is not adjacent to that MeshPart
                auto* mpp = domain_level.get_mesh_node()->find_mesh_part(it);
                if(mpp != nullptr)
                {
                  new_asm->add_mesh_part(*mpp);
                }

                // Insert into the map
                String identifier(it);
                slip_asm.emplace(identifier, new_asm);
              }

            }

          /// Explicitly delete default constructor
          MeshoptAssemblerLevel() = delete;
          /// Explicitly delete the copy constructor
          MeshoptAssemblerLevel(const MeshoptAssemblerLevel&) = delete;
          /// Explicitly delete the move constructor. This could be useful to have, though.
          MeshoptAssemblerLevel(MeshoptAssemblerLevel&&) = delete;

          /**
           * \brief Empty virtual destructor
           */
          virtual ~MeshoptAssemblerLevel()
          {
          }

          template<typename SystemLevel_>
          void assemble_system_filter(SystemLevel_& sys_level)
          {
            // get our global system filter
            typename SystemLevel_::GlobalSystemFilter& fil_glob = sys_level.filter_sys;
            // get our local system filter
            typename SystemLevel_::LocalSystemFilter& fil_loc = fil_glob.local();

            // Assemble homogeneous dirichlet filters
            auto& dirichlet_filters = fil_loc.template at<1>();
            for(auto& it : dirichlet_filters)
            {
              const auto& assembler = dirichlet_asm.find(it.first);
              if(assembler == dirichlet_asm.end())
              {
                throw InternalError(__func__,__FILE__,__LINE__,
                "Could not find dirichlet assembler for filter with key "+it.first);
              }

              assembler->second->assemble(it.second, trafo_space);
            }

            // The slip filter contains the outer unit normal, so reassemble it
            auto& slip_filters = fil_loc.template at<0>();
            for(auto& it : slip_filters)
            {
              const auto& assembler = slip_asm.find(it.first);
              if(assembler == slip_asm.end())
              {
                throw InternalError(__func__,__FILE__,__LINE__,
                "Could not find slip filter assembler for filter with key "+it.first);
              }

              assembler->second->assemble(it.second, trafo_space);
            }

            // Sync the filter vectors in the SlipFilters
            {
              for(auto& it : slip_filters)
              {
                // Get the filter vector
                auto& slip_filter_vector = it.second.get_filter_vector();

                if(slip_filter_vector.used_elements() > 0)
                {
                  // Temporary DenseVector for syncing
                  typename SystemLevel_::LocalSystemVectorL tmp(slip_filter_vector.size());
                  auto* tmp_elements = tmp.template elements<LAFEM::Perspective::native>();
                  auto* sfv_elements = slip_filter_vector.template elements<LAFEM::Perspective::native>();

                  // Copy sparse filter vector contents to DenseVector
                  for(Index isparse(0); isparse < slip_filter_vector.used_elements(); ++isparse)
                  {
                    Index idense(slip_filter_vector.indices()[isparse]);
                    tmp_elements[idense] = sfv_elements[isparse];
                  }

                  sys_level.gate_sys.sync_0(tmp);

                  // Copy sparse filter vector contents to DenseVector
                  for(Index isparse(0); isparse < slip_filter_vector.used_elements(); ++isparse)
                  {
                    Index idense(slip_filter_vector.indices()[isparse]);
                    sfv_elements[isparse] = tmp_elements[idense];
                  }
                }
                else
                {
                  // Temporary DenseVector for syncing
                  typename SystemLevel_::LocalSystemVectorL tmp(slip_filter_vector.size());
                  sys_level.gate_sys.sync_0(tmp);
                }
              }

            } // evil slip filter sync

          }

          template<typename SystemLevel_>
          void assemble_system_filter(SystemLevel_& sys_level,
          /* const */ typename SystemLevel_::GlobalSystemVectorR& vec)
          {
          // get our global system filter
          typename SystemLevel_::GlobalSystemFilter& fil_glob = sys_level.filter_sys;
          // get our local system filter
          typename SystemLevel_::LocalSystemFilter& fil_loc = fil_glob.local();

          typename SystemLevel_::LocalSystemVectorR& vec_loc = vec.local();

          // Assemble Dirichlet boundary conditions from sol, as this->_coords contains the new boundary
          // coordinates
          auto& dirichlet_filters = fil_loc.template at<1>();

          for(auto& it : dirichlet_filters)
          {
            const auto& assembler = dirichlet_asm.find(it.first);
            if(assembler == dirichlet_asm.end())
            {
              throw InternalError(__func__,__FILE__,__LINE__,
              "Could not find dirichlet assembler for filter with key "+it.first);
            }

            assembler->second->assemble(it.second, trafo_space, vec_loc);
          }

          // The slip filter contains the outer unit normal, so reassemble it
          auto& slip_filters = fil_loc.template at<0>();

          for(auto& it : slip_filters)
          {
            const auto& assembler = slip_asm.find(it.first);
            if(assembler == slip_asm.end())
            {
              throw InternalError(__func__,__FILE__,__LINE__,
              "Could not find slip filter assembler for filter with key "+it.first);
            }

            assembler->second->assemble(it.second, trafo_space);
          }

          // Sync the filter vectors in the SlipFilters
          {
            for(auto& it : slip_filters)
            {
              // Get the filter vector
              auto& slip_filter_vector = it.second.get_filter_vector();

              if(slip_filter_vector.used_elements() > 0)
              {
                // Temporary DenseVector for syncing
                typename SystemLevel_::LocalSystemVectorL tmp(slip_filter_vector.size());
                auto* tmp_elements = tmp.template elements<LAFEM::Perspective::native>();
                auto* sfv_elements = slip_filter_vector.template elements<LAFEM::Perspective::native>();

                // Copy sparse filter vector contents to DenseVector
                for(Index isparse(0); isparse < slip_filter_vector.used_elements(); ++isparse)
                {
                  Index idense(slip_filter_vector.indices()[isparse]);
                  tmp_elements[idense] = sfv_elements[isparse];
                }

                sys_level.gate_sys.sync_0(tmp);

                // Copy sparse filter vector contents to DenseVector
                for(Index isparse(0); isparse < slip_filter_vector.used_elements(); ++isparse)
                {
                  Index idense(slip_filter_vector.indices()[isparse]);
                  sfv_elements[isparse] = tmp_elements[idense];
                }
              }
              else
              {
                // Temporary DenseVector for syncing
                typename SystemLevel_::LocalSystemVectorL tmp;
                sys_level.gate_sys.sync_0(tmp);
              }
            }

          } // evil slip filter sync

          }

          template<typename SystemLevel_>
          typename SystemLevel_::GlobalSystemVectorL assemble_rhs_vector(SystemLevel_& sys_level)
          {
            //XASSERTM(!sys_level.local_functional.empty(), "assemble_rhs_vector for empty operator");
            // create new vector
            typename SystemLevel_::GlobalSystemVectorL vec_rhs(sys_level.matrix_sys.create_vector_l());

            vec_rhs.format();

            // get our local system filter
            typename SystemLevel_::LocalSystemFilter& fil_loc = sys_level.filter_sys.local();

            fil_loc.filter_rhs(vec_rhs.local());

            return vec_rhs;
          }

          template<typename SystemLevel_>
          typename SystemLevel_::GlobalSystemVectorR assemble_sol_vector(SystemLevel_& sys_level)
          {
            //XASSERTM(!(*sys_level.op_sys).empty(), "assemble_sol_vector for empty operator");
            typename SystemLevel_::GlobalSystemVectorR vec_sol(sys_level.matrix_sys.create_vector_r());
            vec_sol.clone(sys_level.coords_buffer, LAFEM::CloneMode::Deep);

            return vec_sol;
          }

          template<typename SystemLevel_>
          void assemble_gates(const DomainLayerType& dom_layer, SystemLevel_& sys_level)
          {
            // Create the system gate
            typename SystemLevel_::SystemGate& gate_sys = sys_level.gate_sys;
            // Create the scalar gate
            typename SystemLevel_::ScalarGate& gate_scalar = sys_level.gate_scalar;

            // set the gate comm
            gate_sys.set_comm(dom_layer.get_comm());
            gate_scalar.set_comm(dom_layer.get_comm());

            // Loop over all ranks
            for(Index i(0); i < dom_layer.size(); ++i)
            {
              Index rank = dom_layer.get_rank(i);

              // try to find our halo
              auto* halo = domain_level.find_halo_part(rank);
              XASSERTM(halo != nullptr, "Halo not found.");

              // assemble the system mirror
              typename SystemLevel_::SystemMirror sys_mirror;
              Assembly::MirrorAssembler::assemble_mirror(sys_mirror, this->trafo_space, *halo);

              // push mirror into gates
              gate_sys.push(int(rank), std::move(sys_mirror));

              // assemble the scalar mirror
              typename SystemLevel_::ScalarMirror scalar_mirror;
              Assembly::MirrorAssembler::assemble_mirror(scalar_mirror, this->trafo_space, *halo);

              // push mirror into gates
              gate_scalar.push(int(rank), std::move(scalar_mirror));
            }

            // create local template vectors
            typename SystemLevel_::LocalSystemVectorR tmp_sys(trafo_space.get_num_dofs());
            typename SystemLevel_::LocalScalarVector tmp_scalar(trafo_space.get_num_dofs());

            // compile gates
            gate_sys.compile(std::move(tmp_sys));
            gate_scalar.compile(std::move(tmp_scalar));
          }

          template<typename SystemLevel_>
          void assemble_system_transfer(SystemLevel_& sys_level_fine, MeshoptAssemblerLevel& level_coarse)
          {
            // get global transfer operator
            typename SystemLevel_::GlobalSystemTransfer& glob_trans = sys_level_fine.transfer_sys;

            // get local transfer operator
            typename SystemLevel_::LocalSystemTransfer& loc_trans = glob_trans.local();

            // get local transfer matrices
            typename SystemLevel_::LocalSystemTransferMatrix::BaseClass& loc_prol = loc_trans.get_mat_prol();
            typename SystemLevel_::LocalSystemTransferMatrix::BaseClass& loc_rest = loc_trans.get_mat_rest();

            // assemble structure?
            if (loc_prol.empty())
            {
              Assembly::SymbolicAssembler::assemble_matrix_2lvl(
                loc_prol, this->trafo_space, level_coarse.trafo_space);
            }

            // Create a global scalar weight vector
            typename SystemLevel_::GlobalScalarVector glob_vec_weight(
              &sys_level_fine.gate_scalar, loc_prol.create_vector_l());
            // Get local scalar weight vector
            auto& loc_vec_weight = glob_vec_weight.local();

            // Assemble the underlying scalar prolongation matrix
            {
              loc_prol.format();
              loc_vec_weight.format();

              // Assemble prolongation matrix
              Assembly::GridTransfer::assemble_prolongation(loc_prol, loc_vec_weight,
              this->trafo_space, level_coarse.trafo_space, this->cubature);

              // Synchronise weight vector
              glob_vec_weight.sync_0();

              // Invert components
              loc_vec_weight.component_invert(loc_vec_weight);

              // Scale prolongation matrix
              loc_prol.scale_rows(loc_prol, loc_vec_weight);

              // Copy and transpose
              loc_rest = loc_prol.transpose();
            }
          }

      }; // class MeshoptAssemblerLevel

    } // namespace Meshopt

  } // namespace Control
} //namespace FEAT
#endif // FEAT_CONTROL_MESHOPT_MESHOPT_CONTROL_HPP
