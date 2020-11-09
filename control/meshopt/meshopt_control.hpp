// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

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
     * \brief Control layer for mesh optimization
     */
    namespace Meshopt
    {
      /// \cond internal
      /// Forward declarations
      template<typename, typename, typename, template<typename, typename, typename> class>
      class MeshoptSystemLevel;

      template<typename, typename, typename, template<typename, typename, typename> class>
      class QuadraticSystemLevel;

      template<typename, typename, typename, template<typename, typename, typename> class>
      class NonlinearSystemLevel;

      /// \endcond

      /**
       * \brief Base class for Meshopt control objects
       *
       * \tparam DomainControl_
       * Type of the underlying domain control
       *
       * For every mesh quality functional, there is a control class which inherits from this. This class only knows
       * the bare minimum.
       *
       */
      template<typename DomainControl_>
      class MeshoptControlBase
      {
        public:
          /// The type of the underlying domain control
          typedef DomainControl_ DomainControl;
          /// The type of levels of DomainControl
          typedef typename DomainControl_::LevelType DomainLevel;
          /// The type of mesh the DomainControl uses
          typedef typename DomainControl_::MeshType MeshType;

          /// The world dimension, i.e. number of coordinates
          static constexpr int world_dim = DomainControl_::MeshType::world_dim;
          /// Floating point type for coordinates
          typedef typename DomainControl_::MeshType::CoordType CoordType;

          /// Type for buffer vectors exchanging information between mesh and some vector type on which the solvers
          /// operate
          typedef LAFEM::DenseVectorBlocked<Mem::Main, CoordType, Index, world_dim> LocalCoordsBuffer;
          /// corresponding vector mirror
          typedef LAFEM::VectorMirror<Mem::Main, CoordType, Index> CoordsMirror;
          /// The global version of LocalCoordsBuffer, needed for prepare setting the internal state variable
          typedef Global::Vector<LocalCoordsBuffer, CoordsMirror> GlobalCoordsBuffer;

          /// Type of the vtk exporter this (and derived classes) can write to
          typedef Geometry::ExportVTK<MeshType> VTKExporterType;

        protected:
          /// The domain control whose mesh objects can be modified
          DomainControl& _dom_ctrl;
          /// List of all meshparts with dirichlet boundary conditions
          std::deque<String> _dirichlet_boundaries;
          /// List of all meshparts with slipboundary conditions
          std::deque<String> _slip_boundaries;

        public:

          /**
           * \brief The simplest constructor possible
           */
          explicit MeshoptControlBase(DomainControl& dom_ctrl_, const std::deque<String>& dirichlet_boundaries,
          const std::deque<String>& slip_boundaries)
            : _dom_ctrl(dom_ctrl_),
            _dirichlet_boundaries(dirichlet_boundaries),
            _slip_boundaries(slip_boundaries)
            {
            }

          /**
           * \brief Empty virtual destructor
           */
          virtual ~MeshoptControlBase()
          {
          }

          ///**
          // * \brief Gets a const reference to the mesh at a certain level
          // */
          //const MeshType& get_mesh(Index lvl)
          //{
          //  return _dom_ctrl.at(lvl)->get_mesh();
          //}

          /**
           * \brief Computes a quality indicator concerning the cell sizes
           *
           * \param[out] lambda_min
           * Minimum of the optimal cell size lambda over all cells
           *
           * \param[out] lambda_max
           * Maximum of the optimal cell size lambda over all cells
           *
           * \param[out] vol_min
           * Minimum cell volume
           *
           * \param[out] vol_max
           * Maximum cell volume
           *
           * \param[out] vol
           * Total volume of the domain given by the mesh
           *
           * In a truly optimal mesh (consisting ONLY of reference cells of the right size), every cell's volume is
           * exactly lambda(cell). This is especially the goal for r-adaptivity.
           * So in an optimal mesh,
           * \f[
           *   \forall K \in \mathcal{T}_h: |K|/|\Omega| = \lambda(K)
           * \f]
           * so we compute the 1-norm of the vector
           * \f$(v)_i = \left| \frac{|K_i|}{\sum_j |K_j|} - \lambda(K_i) \right| \f$.
           *
           * \returns The relative cell size quality indicator.
           *
           * \note lambda_min, lambda_max, vol_min, and vol_max are all volume fractions.
           *
           * \note As these quantities are global, this function uses a pre_sync and post_sync part.
           *
           * \see Meshopt::HyperelasticityFunctional::compute_cell_size_defect()
           *
           */
          virtual CoordType compute_cell_size_defect(CoordType& lambda_min, CoordType& lambda_max,
          CoordType& vol_min, CoordType& vol_max, CoordType& vol) const = 0;

          /**
           * \brief Optimizes the current mesh
           */
          virtual void optimize() = 0;

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
          const std::deque<String>& get_dirichlet_boundaries() const
          {
            return _dirichlet_boundaries;
          }

          /**
           * \brief Gets the names of all slip boundaries
           *
           * Note that each name refers to a boundary, but that boundary does not necessarily have to be present due
           * to partitioning etc.
           *
           * \returns A deque of Strings with all slip boundary names
           */
          const std::deque<String>& get_slip_boundaries() const
          {
            return _slip_boundaries;
          }

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
           *
           * \param[in,out] exporter
           * The vtk exporter to add our data to
           *
           * \param[in] lvl_index
           * This level's data gets added.
           */
          virtual void add_to_vtk_exporter
            (Geometry::ExportVTK<MeshType>& DOXY(exporter), const int DOXY(lvl_index)) const
            {
            }

          /**
           * \brief Prints settings of the control object
           */
          virtual String info() const = 0;


        /**
         * \brief Computes mesh quality heuristics
         *
         * \param[out] edge_angle
         * The worst angle between two edges. Keep in mind that in 3d, this can be >0 even for deteriorated cells.
         *
         * \param[out] qi_min
         * The minimum quality indicator over all cells.
         *
         * \param[out] qi_mean
         * The mean quality indicator overa all cells.
         *
         * \param[out] edge_angle_cellwise
         * For debugging or visualization purposes, this can receive the worst edge angle for every cell.
         *
         * \param[out] qi_cellwise
         * For debugging or visualization purposes, this can receive the quality indicator for every cell.
         *
         * \param[in] lvl_index
         * Index of the level to compute everything for. Defaults to the maximum level.
         *
         */
        void compute_mesh_quality(CoordType& edge_angle, CoordType& qi_min, CoordType& qi_mean,
          CoordType* edge_angle_cellwise = nullptr, CoordType* qi_cellwise = nullptr,
          int lvl_index = -1) const
        {
          // max_level_index cannot be called for the default argument, so we do it here
          if(lvl_index == -1)
          {
            lvl_index = _dom_ctrl.max_level_index();
          }
          XASSERT(lvl_index >= _dom_ctrl.min_level_index());
          XASSERT(lvl_index <= _dom_ctrl.max_level_index());

          const Dist::Comm& comm = _dom_ctrl.comm();

          CoordType qi_sum(0);

          for(std::size_t i(0u); i < _dom_ctrl.size_physical(); ++i)
          {
            const auto& dom_lvl = *_dom_ctrl.at(i);
            if(dom_lvl.get_level_index() == lvl_index)
            {
              const auto& my_mesh = dom_lvl.get_mesh();

              Index ncells(my_mesh.get_num_entities(MeshType::shape_dim));
              comm.allreduce(&ncells, &ncells, std::size_t(1), Dist::op_sum);

              Geometry::MeshQualityHeuristic<typename MeshType::ShapeType>::compute(qi_min, qi_sum,
              my_mesh.template get_index_set<MeshType::shape_dim, 0>(), my_mesh.get_vertex_set(), qi_cellwise);

              comm.allreduce(&qi_min, &qi_min, std::size_t(1), Dist::op_min);
              comm.allreduce(&qi_sum, &qi_sum, std::size_t(1), Dist::op_sum);

              edge_angle = Geometry::MeshQualityHeuristic<typename MeshType::ShapeType>::angle(
              my_mesh.template get_index_set<MeshType::shape_dim, 0>(), my_mesh.get_vertex_set(), edge_angle_cellwise);
              comm.allreduce(&edge_angle, &edge_angle, std::size_t(1), Dist::op_min);

              qi_mean = qi_sum/CoordType(ncells);

              return;
            }
          }

          // We should never get to this point
          XABORTM("Could not find level with index "+stringify(lvl_index)+"!\n");

        }
      }; // class MeshoptControlBase


      /**
       * \brief SystemLevel base class
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
       * \tparam Functional_
       * The (patch-) local mesh quality functional
       *
       */
      template
      <
        typename Mem_, typename DT_, typename IT_,
        template<typename, typename, typename> class Functional_
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

          /// The (patch-)local mesh quality functional
          typedef Functional_<Mem_, DT_, IT_> LocalFunctional;

          /// The finite element space
          typedef typename LocalFunctional::SpaceType SpaceType;
          /// The underlying transformation
          typedef typename LocalFunctional::TrafoType TrafoType;
          /// The type of mesh we use
          typedef typename TrafoType::MeshType MeshType;

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
          /// The (patch-)local mesh quality functional
          LocalFunctional local_functional;
          /// The scalar gate
          ScalarGate gate_scalar;
          /// The system gate
          SystemGate gate_sys;
          /// The global filter
          GlobalSystemFilter filter_sys;
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
           * \tparam Args_
           * For passing to the constructor of local_functional
           *
           * \param[in] level_index
           * The refinement level since the mesh was constructed
           *
           * \param[in] rmn
           * The root mesh node on this level
           *
           * \param[in] trafo
           * The transformation on this level
           *
           * \param[in] dirichlet_list
           * List of the names of all Dirichlet boundaries
           *
           * \param[in] slip_list
           * List of the names of all slip boundaries
           *
           * \param[in] args
           * Additional arguments that get passed to the constructor of local_functional
           *
           */
          template<typename... Args_>
          explicit MeshoptSystemLevel(const int level_index,
          Geometry::RootMeshNode<MeshType>* rmn,
          TrafoType& trafo,
          const std::deque<String>& dirichlet_list,
          const std::deque<String>& slip_list,
          Args_&&... args) :
            local_functional(rmn, trafo, dirichlet_list, slip_list, std::forward<Args_>(args)...),
            gate_sys(),
            filter_sys(),
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

          /// Explicitly delete default constructor
          MeshoptSystemLevel() = delete;
          /// Explicitly delete the copy constructor
          MeshoptSystemLevel(const MeshoptSystemLevel&) = delete;
          /// Explicitly delete the move constructor. This could be useful to have, though.
          MeshoptSystemLevel(MeshoptSystemLevel&&) = delete;

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
           *
           * \returns This SystemLevel's level index.
           */
          int get_level_index() const
          {
            return _level_index;
          }

          /**
           * \brief Returns if the local functional is empty
           *
           * \returns True if the local functional is empty
           */
          bool empty() const
          {
            return local_functional.empty();
          }

          /**
           * \brief Creates a new (left) vector
           */
          GlobalSystemVectorL create_vector_l() const
          {
            return GlobalSystemVectorL(&gate_sys, local_functional.create_vector_l());
          }

          /**
           * \brief Creates a new (right) vector
           */
          GlobalSystemVectorR create_vector_r() const
          {
            return GlobalSystemVectorR(&gate_sys, local_functional.create_vector_r());
          }


          /**
           * \brief Synchronizes the system filters
           *
           * If there is a global functional (like Global::Functional<SomeClass>) that handles the synchronization
           * after each call to prepare(), is is not needed. For quadratic functionals, where there is no global
           * functional but only a Global::Matrix, this routine is provided.
           */
          void sync_system_filter()
          {
            // Sync the filter vectors in the SlipFilters
            auto& slip_filters = filter_sys.local().template at<0>();
            {
              for(auto& it : slip_filters)
              {
                // Get the filter vector
                auto& slip_filter_vector = it.second.get_filter_vector();

                if(slip_filter_vector.used_elements() > 0)
                {
                  // Temporary DenseVector for syncing
                  LocalSystemVectorL tmp(slip_filter_vector.size());
                  auto* tmp_elements = tmp.template elements<LAFEM::Perspective::native>();
                  auto* sfv_elements = slip_filter_vector.template elements<LAFEM::Perspective::native>();

                  // Copy sparse filter vector contents to DenseVector
                  for(Index isparse(0); isparse < slip_filter_vector.used_elements(); ++isparse)
                  {
                    Index idense(slip_filter_vector.indices()[isparse]);
                    tmp_elements[idense] = sfv_elements[isparse];
                  }

                  gate_sys.sync_0(tmp);

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
                  LocalSystemVectorL tmp(slip_filter_vector.size());
                  gate_sys.sync_0(tmp);
                }
              }

            } // evil slip filter sync

          }

          /**
           * \brief Assembles a right hand side vector
           *
           * Note that this is always homogeneous and filtered.
           *
           * \returns A right hand side vector for the mesh optimization problem on this level.
           */
          GlobalSystemVectorL assemble_rhs_vector() const
          {
            XASSERTM(!empty(), "Assemble_rhs_vector for empty functional");

            // create new vector
            GlobalSystemVectorL vec_rhs(create_vector_l());

            vec_rhs.format();
            filter_sys.filter_rhs(vec_rhs);

            return vec_rhs;
          }

          /**
           * \brief Assembles an intial guess vector
           *
           * \returns A vector to be passed to an iterative solver as initial guess.
           */
          GlobalSystemVectorR assemble_sol_vector() const
          {
            XASSERTM(!empty(), "Assemble_sol_vector for empty functional");
            GlobalSystemVectorR vec_sol(create_vector_r());

            vec_sol.local().copy(coords_buffer.local());
            filter_sys.local().template at<1>().filter_sol(vec_sol.local());

            return vec_sol;
          }

          /**
           * \brief Assembles the gates on this level on a given DomainLayer
           *
           * \tparam DomainLayer_
           * Type for the layer
           *
           * \param[in] dom_layer
           * The layer to assemble the gates on
           *
           */
          template<typename DomainLayer_>
          void assemble_gates(const DomainLayer_& dom_layer)
          {
            // set the gate comm
            gate_sys.set_comm(dom_layer.comm_ptr());
            gate_scalar.set_comm(dom_layer.comm_ptr());

            // Loop over all ranks
            for(Index i(0); i < dom_layer.neighbor_count(); ++i)
            {
              int rank(dom_layer.neighbor_rank(i));

              // try to find our halo
              auto* halo = local_functional.get_mesh_node()->get_halo(rank);
              XASSERTM(halo != nullptr, "Halo not found.");

              // assemble the system mirror
              SystemMirror sys_mirror;
              Assembly::MirrorAssembler::assemble_mirror(sys_mirror, local_functional.trafo_space, *halo);

              // push mirror into gates
              gate_sys.push(rank, std::move(sys_mirror));

              // assemble the scalar mirror
              ScalarMirror scalar_mirror;
              Assembly::MirrorAssembler::assemble_mirror(scalar_mirror, local_functional.trafo_space, *halo);

              // push mirror into gates
              gate_scalar.push(rank, std::move(scalar_mirror));
            }

            // create local template vectors
            LocalSystemVectorR tmp_sys(local_functional.trafo_space.get_num_dofs());
            LocalScalarVector tmp_scalar(local_functional.trafo_space.get_num_dofs());

            // compile gates
            gate_sys.compile(std::move(tmp_sys));
            gate_scalar.compile(std::move(tmp_scalar));
          }

          /**
           * \brief Assembles the transfer operators from this level to a coarser one
           *
           * \param[in] level_coarse
           * Target coarse level for the transfer operators.
           *
           */
          void assemble_system_transfer(const MeshoptSystemLevel& level_coarse)
          {
            // get global transfer operator
            GlobalSystemTransfer& glob_trans = this->transfer_sys;

            // get local transfer operator
            LocalSystemTransfer& loc_trans = glob_trans.local();

            // get local transfer matrices
            typename LocalSystemTransferMatrix::BaseClass& loc_prol = loc_trans.get_mat_prol();
            typename LocalSystemTransferMatrix::BaseClass& loc_rest = loc_trans.get_mat_rest();

            // assemble structure?
            if (loc_prol.empty())
            {
              Assembly::SymbolicAssembler::assemble_matrix_2lvl(
                loc_prol, local_functional.trafo_space, level_coarse.local_functional.trafo_space);
            }

            // Create a global scalar weight vector
            GlobalScalarVector glob_vec_weight( &this->gate_scalar, loc_prol.create_vector_l());
            // Get local scalar weight vector
            auto& loc_vec_weight = glob_vec_weight.local();

            // Assemble the underlying scalar prolongation matrix
            {
              Cubature::DynamicFactory cubature_factory("auto-degree:3");

              loc_prol.format();
              loc_vec_weight.format();

              // Assemble prolongation matrix
              Assembly::GridTransfer::assemble_prolongation(loc_prol, loc_vec_weight,
              local_functional.trafo_space, level_coarse.local_functional.trafo_space, cubature_factory);

              // Synchronize weight vector
              glob_vec_weight.sync_0();

              // Invert components
              loc_vec_weight.component_invert(loc_vec_weight);

              // Scale prolongation matrix
              loc_prol.scale_rows(loc_prol, loc_vec_weight);

              // Copy and transpose
              loc_rest = loc_prol.transpose();
            }
          }

      }; // class MeshoptSystemLevel

      /**
       * \brief SystemLevel for a quadratic mesh quality functional leading to a linear system
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
       * \tparam Functional_
       * The (patch-) local quadratic mesh quality functional
       *
       * Since the mesh quality functional is quadratic, it can be minimized by solving the the linear system of
       * equations given by its gradient, represented by the GlobalSystemMatrix.
       *
       */
      template
      <
        typename Mem_, typename DT_, typename IT_,
        template<typename, typename, typename> class Functional_
      >
      class QuadraticSystemLevel : public MeshoptSystemLevel<Mem_, DT_, IT_, Functional_>
      {
        public:
          /// Memory architecture for the solver
          typedef Mem_ MemType;
          /// Floating point precision for the solver
          typedef DT_ DataType;
          /// Index type for the solver
          typedef IT_ IndexType;

          /// Our base class
          typedef MeshoptSystemLevel<Mem_, DT_, IT_, Functional_> BaseClass;

          /// (Patch-) Local mesh quality functional type
          typedef Functional_<Mem_, DT_, IT_> LocalFunctional;

          /// The finite element space
          typedef typename LocalFunctional::SpaceType SpaceType;
          /// The underlying transformation
          typedef typename LocalFunctional::TrafoType TrafoType;
          /// The type of mesh we use
          typedef typename TrafoType::MeshType MeshType;

          /// The local system matrix type (for the gradient)
          typedef typename LocalFunctional::MatrixType LocalMatrix;
          /// Local left-vectors (dual space)
          typedef typename BaseClass::LocalSystemVectorL LocalSystemVectorL;
          /// Local right-vectors (primal space)
          typedef typename BaseClass::LocalSystemVectorR LocalSystemVectorR;
          /// Local vectors of scalar quantities
          typedef typename BaseClass::LocalScalarVector LocalScalarVector;
          /// Local coordinates buffer type for passing information to or from the mesh
          typedef typename BaseClass::LocalCoordsBuffer LocalCoordsBuffer;
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
          typedef typename BaseClass::SystemMirror SystemMirror;
          /// Gates for the system
          typedef typename BaseClass::SystemGate SystemGate;
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
          /// The global system matrix
          GlobalSystemMatrix matrix_sys;

        public:
          /**
           * \brief Variadic template constructor
           *
           * \tparam Args_
           * For passing to the constructor of local_functional
           *
           * \param[in] level_index
           * The refinement level since the mesh was constructed
           *
           * \param[in] rmn
           * The root mesh node on this level
           *
           * \param[in] trafo
           * The transformation on this level
           *
           * \param[in] dirichlet_list
           * List of the names of all Dirichlet boundaries
           *
           * \param[in] slip_list
           * List of the names of all slip boundaries
           *
           * \param[in] args
           * Additional arguments that get passed to the constructor of local_functional
           *
           */
          template<typename... Args_>
          explicit QuadraticSystemLevel(const int level_index,
          Geometry::RootMeshNode<MeshType>* rmn,
          TrafoType& trafo,
          const std::deque<String>& dirichlet_list,
          const std::deque<String>& slip_list,
          Args_&&... args) :
            BaseClass(level_index, rmn, trafo, dirichlet_list, slip_list, std::forward<Args_>(args)...),
            matrix_sys(&(BaseClass::gate_sys), &(BaseClass::gate_sys),
            BaseClass::local_functional.matrix_sys.clone(LAFEM::CloneMode::Shallow))
            {
            }

          /**
           * \brief Empty virtual destructor
           */
          virtual ~QuadraticSystemLevel()
          {
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
       * \tparam Functional_
       * The (patch-) local nonlinear, nonquadratic mesh quality functional
       *
       * Since the mesh quality functional is not quadratic, the gradient is a nonlinear mapping.
       *
       */
      template
      <
        typename Mem_, typename DT_, typename IT_,
        template<typename, typename, typename> class Functional_
      >
      class NonlinearSystemLevel : public MeshoptSystemLevel<Mem_, DT_, IT_, Functional_>
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

          /// Our base class
          typedef MeshoptSystemLevel<Mem_, DT_, IT_, Functional_> BaseClass;

          /// The finite element space
          typedef typename LocalFunctional::SpaceType SpaceType;
          /// The underlying transformation
          typedef typename LocalFunctional::TrafoType TrafoType;
          /// The type of mesh we use
          typedef typename TrafoType::MeshType MeshType;

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
          typedef LAFEM::VectorMirror<Mem_, DT_, IT_> SystemMirror;
          /// Gates for the system
          typedef Global::Gate<LocalSystemVectorR, SystemMirror> SystemGate;
          /// Mirrors for scalar vectors
          typedef LAFEM::VectorMirror<Mem_, DT_, IT_> ScalarMirror;
          /// Gates for scalar vectors
          typedef Global::Gate<LocalScalarVector, ScalarMirror> ScalarGate;

          /// Global mesh quality functional type
          typedef Global::NonlinearFunctional<LocalFunctional, SystemMirror, SystemMirror> GlobalFunctional;
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
          /// The global nonlinear functional
          GlobalFunctional global_functional;

        public:
          /**
           * \brief Variadic template constructor
           *
           * \tparam Args_
           * For passing to the constructor of local_functional
           *
           * \param[in] level_index
           * The refinement level since the mesh was constructed
           *
           * \param[in] rmn
           * The root mesh node on this level
           *
           * \param[in] trafo
           * The transformation on this level
           *
           * \param[in] dirichlet_list
           * List of the names of all Dirichlet boundaries
           *
           * \param[in] slip_list
           * List of the names of all slip boundaries
           *
           * \param[in] args
           * Additional arguments that get passed to the constructor of local_functional
           *
           */
          template<typename... Args_>
          explicit NonlinearSystemLevel(const int level_index,
          Geometry::RootMeshNode<MeshType>* rmn,
          TrafoType& trafo,
          const std::deque<String>& dirichlet_list,
          const std::deque<String>& slip_list,
          Args_&&... args) :
            BaseClass(level_index, rmn, trafo, dirichlet_list, slip_list, std::forward<Args_>(args)...),
            global_functional(&(BaseClass::gate_sys), &(BaseClass::gate_sys), BaseClass::local_functional)
            {
            }

          /**
           * \brief Empty virtual destructor
           */
          virtual ~NonlinearSystemLevel()
          {
          }

      }; // class NonlinearSystemLevel<...>

    } // namespace Meshopt

  } // namespace Control
} //namespace FEAT
#endif // FEAT_CONTROL_MESHOPT_MESHOPT_CONTROL_HPP
