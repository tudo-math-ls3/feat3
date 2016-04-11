#pragma once
#ifndef CONTROL_SCALAR_BASIC_HPP
#define CONTROL_SCALAR_BASIC_HPP 1

#include <kernel/foundation/comm_base.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/mean_filter.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/global/foundation_gate.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/filter.hpp>
#include <kernel/global/mean_filter.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/grid_transfer.hpp>
#include <kernel/assembly/mirror_assembler.hpp>
#include <kernel/solver/base.hpp>
#include <kernel/solver/iterative.hpp>
#include <kernel/util/property_map.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/bicgstab.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/scale_precond.hpp>
#include <kernel/solver/ilu_precond.hpp>
#include <kernel/solver/ssor_precond.hpp>
#include <kernel/solver/schwarz_precond.hpp>

#include <control/domain/domain_control.hpp>

#include <deque>

namespace FEAST
{
  namespace Control
  {
    template<
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_> >
    class ScalarBasicSystemLevel
    {
    public:
      static_assert(std::is_same<MemType_, typename ScalarMatrix_::MemType>::value, "MemType missmatch!");
      static_assert(std::is_same<DataType_, typename ScalarMatrix_::DataType>::value, "DataType missmatch!");
      static_assert(std::is_same<IndexType_, typename ScalarMatrix_::IndexType>::value, "IndexType missmatch!");

      // basic types
      typedef MemType_ MemType;
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;

      /// define local scalar matrix type
      typedef ScalarMatrix_ LocalScalarMatrix;

      /// define local system matrix type
      typedef ScalarMatrix_ LocalSystemMatrix;

      /// define local system vector type
      typedef typename LocalSystemMatrix::VectorTypeR LocalSystemVector;

      /// define system mirror type
      typedef LAFEM::VectorMirror<MemType_, DataType, IndexType> SystemMirror;

      /// define system gate
      typedef Global::FoundationGate<LocalSystemVector, SystemMirror> SystemGate;

      /// define global system vector type
      typedef Global::Vector<LocalSystemVector> GlobalSystemVector;

      /// define global system matrix type
      typedef Global::Matrix<LocalSystemMatrix> GlobalSystemMatrix;

      /* ***************************************************************************************** */

      /// our system gate
      SystemGate gate_sys;

      /// our global system matrix
      GlobalSystemMatrix matrix_sys;

      /// CTOR
      ScalarBasicSystemLevel() :
        matrix_sys(&gate_sys, &gate_sys)
      {
      }

      virtual ~ScalarBasicSystemLevel()
      {
      }

      template<typename M_, typename D_, typename I_, typename SM_>
      void convert(const ScalarBasicSystemLevel<M_, D_, I_, SM_> & other)
      {
        gate_sys.convert(other.gate_sys);
        matrix_sys.convert(&gate_sys, &gate_sys, other.matrix_sys);
      }
    }; // class ScalarBasicSystemLevel<...>

    template<
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_> >
    class ScalarUnitFilterSystemLevel :
      public ScalarBasicSystemLevel<MemType_, DataType_, IndexType_, ScalarMatrix_>
    {
    public:
      typedef ScalarBasicSystemLevel<MemType_, DataType_, IndexType_, ScalarMatrix_> BaseClass;

      /// define local filter type
      typedef LAFEM::UnitFilter<MemType_, DataType_, IndexType_> LocalSystemFilter;

      /// define global filter type
      typedef Global::Filter<LocalSystemFilter> GlobalSystemFilter;

      /// our global system filter
      GlobalSystemFilter filter_sys;

      /// CTOR
      ScalarUnitFilterSystemLevel() :
        BaseClass()
      {
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return this->matrix_sys.bytes () + filter_sys.bytes();
      }

      /**
       *
       * \brief Conversion method
       *
       * Use source ScalarUnitFilterSystemLevel content as content of current ScalarUnitFilterSystemLevel.
       *
       */
      template<typename M_, typename D_, typename I_, typename SM_>
      void convert(const ScalarUnitFilterSystemLevel<M_, D_, I_, SM_> & other)
      {
        BaseClass::convert(other);
        filter_sys.convert(other.filter_sys);
      }
    }; // class ScalarUnitFilterSystemLevel<...>

    template<
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_> >
    class ScalarMeanFilterSystemLevel :
      public ScalarBasicSystemLevel<MemType_, DataType_, IndexType_, ScalarMatrix_>
    {
    public:
      typedef ScalarBasicSystemLevel<MemType_, DataType_, IndexType_, ScalarMatrix_> BaseClass;

      /// define local filter type
      typedef Global::MeanFilter<MemType_, DataType_, IndexType_> LocalSystemFilter;

      /// define global filter type
      typedef Global::Filter<LocalSystemFilter> GlobalSystemFilter;

      /// our global system filter
      GlobalSystemFilter filter_sys;

      /// CTOR
      ScalarMeanFilterSystemLevel() :
        BaseClass()
      {
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return this->matrix_sys.bytes () + filter_sys.bytes();
      }

      /**
       *
       * \brief Conversion method
       *
       * Use source ScalarMeanFilterSystemLevel content as content of current ScalarMeanFilterSystemLevel.
       *
       */
      template<typename M_, typename D_, typename I_, typename SM_>
      void convert(const ScalarMeanFilterSystemLevel<M_, D_, I_, SM_> & other)
      {
        BaseClass::convert(other);
        filter_sys.convert(other.filter_sys);
      }
    }; // class ScalarMeanFilterSystemLevel<...>

    template<typename SystemLevel_, typename ScalarMatrix_ = typename SystemLevel_::LocalScalarMatrix>
    class ScalarBasicTransferLevel
    {
    public:
      /// our local transfer matrix type
      typedef ScalarMatrix_ LocalSystemTransferMatrix;

      /// our global transfer matrix type
      typedef Global::Matrix<LocalSystemTransferMatrix> GlobalSystemTransferMatrix;

      /// our global transfer matrices
      GlobalSystemTransferMatrix prol_sys;

      /// \copydoc ScalarBasicTransferLevel::prol_sys
      GlobalSystemTransferMatrix rest_sys;

      ScalarBasicTransferLevel()
      {
      }

      /// CTOR
      explicit ScalarBasicTransferLevel(SystemLevel_& lvl_coarse, SystemLevel_& lvl_fine) :
        prol_sys(&lvl_fine.gate_sys, &lvl_coarse.gate_sys),
        rest_sys(&lvl_coarse.gate_sys, &lvl_fine.gate_sys)
      {
      }

      virtual ~ScalarBasicTransferLevel()
      {
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return prol_sys.bytes() + rest_sys.bytes();
      }

      /**
       *
       * \brief Conversion method
       *
       * Use source ScalarBasicTransferLevel content as content of current ScalarBasicTransferLevel.
       *
       * \warning The provided SystemLevels must already be converted to the matching
       * configuration, as they contain the used gateways.
       *
       */
      template <typename SL_, typename SM_>
      void convert(SystemLevel_ & lvl_coarse , SystemLevel_ & lvl_fine, const ScalarBasicTransferLevel<SL_, SM_> & other)
      {
        prol_sys.convert(&lvl_fine.gate_sys, &lvl_coarse.gate_sys, other.prol_sys);
        rest_sys.convert(&lvl_coarse.gate_sys, &lvl_fine.gate_sys, other.rest_sys);
      }
    }; // class ScalarBasicTransferLevel<...>

    template<typename Space_>
    class ScalarBasicAssemblerLevel
    {
    public:
      typedef Space_ SpaceType;
      typedef typename SpaceType::TrafoType TrafoType;
      typedef typename TrafoType::MeshType MeshType;
      typedef Control::Domain::DomainLevel<MeshType> DomainLevelType;
      typedef Control::Domain::DomainLayer<MeshType> DomainLayerType;

      DomainLevelType& domain_level;
      MeshType& mesh;
      TrafoType trafo;
      SpaceType space;
      Cubature::DynamicFactory cubature;

      explicit ScalarBasicAssemblerLevel(DomainLevelType& dom_lvl) :
        domain_level(dom_lvl),
        mesh(domain_level.get_mesh()),
        trafo(mesh),
        space(trafo),
        cubature("auto-degree:" + stringify(Math::sqr(SpaceType::local_degree)+2))
      {
      }

      virtual ~ScalarBasicAssemblerLevel()
      {
      }

      //template<typename MemType_, typename DataType_, typename IndexType_, template<typename, typename, typename> class ScalarMatrix_>
      //void assemble_gates(const DomainLayerType& dom_layer, ScalarBasicSystemLevel<MemType_, DataType_, IndexType_, ScalarMatrix_>& sys_level)
      template<typename SystemLevel_>
      void assemble_gates(const DomainLayerType& dom_layer, SystemLevel_& sys_level)
      {
        // get our gate
        typename SystemLevel_::SystemGate& gate_sys = sys_level.gate_sys;

        // loop over all ranks
        for(Index i(0); i < dom_layer.size(); ++i)
        {
          Index rank = dom_layer.get_rank(i);
          Index ctag = dom_layer.get_ctag(i);

          // try to find our halo
          auto* halo = domain_level.find_halo_part(rank);
          if (halo == nullptr)
            throw InternalError("ERROR: Halo not found!");

          // assemble the mirror
          typename SystemLevel_::SystemMirror mirror_sys;
          Assembly::MirrorAssembler::assemble_mirror(mirror_sys, space, *halo);

          // push mirror into gate
          gate_sys.push(rank,ctag, std::move(mirror_sys));
        }

        // create local template vector
        typename SystemLevel_::LocalSystemVector tmpl_s(space.get_num_dofs());

        // compile gate
        gate_sys.compile(std::move(tmpl_s));
      }

      template<typename TransferLevel_>
      void assemble_system_transfer(TransferLevel_& trans_level, ScalarBasicAssemblerLevel& level_coarse)
      {
        // get global transfer matrices
        typename TransferLevel_::GlobalSystemTransferMatrix& glob_prol = trans_level.prol_sys;
        typename TransferLevel_::GlobalSystemTransferMatrix& glob_rest = trans_level.rest_sys;

        // get local transfer matrices
        typename TransferLevel_::LocalSystemTransferMatrix& loc_prol = (*glob_prol);
        typename TransferLevel_::LocalSystemTransferMatrix& loc_rest = (*glob_rest);

        // assemble structure?
        if (loc_prol.empty())
        {
          Assembly::SymbolicMatrixAssembler<Assembly::Stencil::StandardRefinement>::assemble(
            loc_prol, this->space, level_coarse.space);
        }

        // create a global pressure weight vector
        auto glob_vec_weight = glob_prol.create_vector_l();

        // get local pressure weight vector
        auto& loc_vec_weight = (*glob_vec_weight);

        // assemble prolongation matrix
        {
          loc_prol.format();
          loc_vec_weight.format();

          // assemble prolongation matrix
          Assembly::GridTransfer::assemble_prolongation(loc_prol, loc_vec_weight,
            this->space, level_coarse.space, this->cubature);

          // synchronise weight vector
          glob_vec_weight.sync_0();

          // invert components
          loc_vec_weight.component_invert(loc_vec_weight);

          // scale prolongation matrix
          loc_prol.scale_rows(loc_prol, loc_vec_weight);

          // copy and transpose
          loc_rest = loc_prol.transpose();
        }
      }
    }; // class ScalarBasicAssemblerLevel<...>

    struct SolverFactory
    {
      private:

      /// \todo make internal
      template <typename VectorType_>
      static void configure_iterative_solver(PropertyMap * section, std::shared_ptr<Solver::PreconditionedIterativeSolver<VectorType_> > solver)
      {
        Index rank = Foundation::Comm::rank();

        auto plot_p = section->get_entry("plot");
        if (plot_p.second)
        {
          Index plot(std::stoul(plot_p.first));
          if (plot == 0)
          {
            solver->set_plot(false);
          }
          else if (plot == 1)
          {
            solver->set_plot(rank == 0);
          }
          else
          {
            throw InternalError(__func__, __FILE__, __LINE__, "plot value " + stringify(plot) + " unknown!");
          }
        }

        auto tol_rel_p = section->get_entry("tol_rel");
        if (tol_rel_p.second)
          solver->set_tol_rel(std::stod(tol_rel_p.first));

        auto max_iter_p = section->get_entry("max_iter");
        if (max_iter_p.second)
          solver->set_max_iter(std::stoul(max_iter_p.first));

        auto min_iter_p = section->get_entry("min_iter");
        if (min_iter_p.second)
          solver->set_min_iter(std::stoul(min_iter_p.first));
      }

      /// \todo make internal
      static String get_section_path(PropertyMap * base, PropertyMap * section, String path, String name)
      {
        if (section->get_sub_section(name) != nullptr)
        {
          return path + "/" + name;
        }
        else if (base->get_sub_section(name) != nullptr)
        {
          return name;
        }
        else
          throw InternalError(__func__, __FILE__, __LINE__, "section " + name + " not found!");
      }

      /// returns the object, if T_ has a GateType, i.e. is a GlobalVector - SFINAE at its best
      template <typename Evaluator_, typename T_>
      static T_ & derefer(T_ & object, typename Evaluator_::GateType *)
      {
        return object;
      }

      /// returns the dereferenced object, if T_ has no GateType, i.e. is a LocalVector - SFINAE at its best
      template <typename Evaluator_, typename T_>
      static auto derefer(T_ & object, ...) -> decltype(*object)
      {
        return *object;
      }

      template <typename Evaluator_, typename SystemLevelType_, typename TransferLevelType_>
      static std::shared_ptr<Solver::SolverBase<typename SystemLevelType_::GlobalSystemVector> >
      create_schwarz_precon(std::deque<SystemLevelType_*> & system_levels, std::deque<TransferLevelType_*> & transfer_levels, PropertyMap * base, String solver_name, PropertyMap * section, typename Evaluator_::GateType *)
      {
        typedef typename SystemLevelType_::GlobalSystemVector SolverVectorType;
        std::shared_ptr<Solver::SolverBase<typename SolverVectorType::LocalVectorType> > precon_schwarz;
        auto schwarz_p = section->query("solver");
        if (schwarz_p.second)
        {
          precon_schwarz = create_scalar_solver<SystemLevelType_, TransferLevelType_, typename SolverVectorType::LocalVectorType>(system_levels, transfer_levels, base, get_section_path(base, section, solver_name, schwarz_p.first));
        }
        else
          throw InternalError(__func__, __FILE__, __LINE__, "Schwarz precon section without solver key is not allowed!");

        return Solver::new_schwarz_precond(precon_schwarz, system_levels.back()->filter_sys);
      }

      template <typename Evaluator_, typename SystemLevelType_, typename TransferLevelType_>
      static std::shared_ptr<Solver::SolverBase<typename SystemLevelType_::GlobalSystemVector::LocalVectorType> >
      create_schwarz_precon(std::deque<SystemLevelType_*> & /*system_levels*/, std::deque<TransferLevelType_*> & /*transfer_levels*/, PropertyMap * /*base*/, String /*solver_name*/,
        PropertyMap * /*section*/, ...)
      {
        throw InternalError(__func__, __FILE__, __LINE__, "Schwarz precon section is only allowed in global context! Maybe you have two in one solver branch?");
        return nullptr;
      }

      template <typename Evaluator_, typename SystemLevelType_>
      static std::shared_ptr<Solver::SolverBase<typename SystemLevelType_::GlobalSystemVector> >
      create_ilu_precon(std::deque<SystemLevelType_*> & /*system_levels*/, typename Evaluator_::GateType *)
      {
        throw InternalError(__func__, __FILE__, __LINE__, "ilu precon section is only allowed in local context!");
        return nullptr;
      }

      template <typename Evaluator_, typename SystemLevelType_>
      static std::shared_ptr<Solver::SolverBase<typename SystemLevelType_::GlobalSystemVector::LocalVectorType> >
      create_ilu_precon(std::deque<SystemLevelType_*> & system_levels, ...)
      {
        auto result = Solver::new_ilu_precond(*system_levels.back()->matrix_sys, *system_levels.back()->filter_sys, 0ul);
        return result;
      }

      template <typename Evaluator_, typename SystemLevelType_>
      static std::shared_ptr<Solver::SolverBase<typename SystemLevelType_::GlobalSystemVector> >
      create_ssor_precon(std::deque<SystemLevelType_*> & /*system_levels*/, PropertyMap * /*section*/, typename Evaluator_::GateType *)
      {
        throw InternalError(__func__, __FILE__, __LINE__, "ssor precon section is only allowed in local context!");
        return nullptr;
      }

      template <typename Evaluator_, typename SystemLevelType_>
      static std::shared_ptr<Solver::SolverBase<typename SystemLevelType_::GlobalSystemVector::LocalVectorType> >
      create_ssor_precon(std::deque<SystemLevelType_*> & system_levels, PropertyMap * section, ...)
      {
        auto omega_p = section->get_entry("omega");
        if (omega_p.second)
        {
          auto result = Solver::new_ssor_precond(*system_levels.back()->matrix_sys, *system_levels.back()->filter_sys, std::stod(omega_p.first));
          return result;
        }
        else
        {
          auto result = Solver::new_ssor_precond(*system_levels.back()->matrix_sys, *system_levels.back()->filter_sys);
          return result;
        }
      }

      public:

      /**
       * \brief Create solver tree based on PropertyMap
       *
       * \param[in] system_levels SystemLevel hierarchy
       * \param[in] transfer_levels TransferLevel hierarchy
       * \param[in] base A pointer to the PropertyMap that contains all solver related informations
       * \param[in] solver_name The name of the solver tree's root section
       */
      template <typename SystemLevelType_, typename TransferLevelType_, typename SolverVectorType_ = typename SystemLevelType_::GlobalSystemVector>
      static std::shared_ptr<Solver::SolverBase<SolverVectorType_> >
      create_scalar_solver(std::deque<SystemLevelType_*> & system_levels, std::deque<TransferLevelType_*> & transfer_levels, PropertyMap * base, String solver_name)
      {
        std::shared_ptr<Solver::SolverBase<SolverVectorType_> > result;

        auto section = base->query_section(solver_name);

        auto solver_p = section->query("type");
        if (!solver_p.second)
          throw InternalError(__func__, __FILE__, __LINE__, "no type key found in property map: " + solver_name + "!");
        String solver_type = solver_p.first;

        std::shared_ptr<Solver::SolverBase<SolverVectorType_> > precon = nullptr;
        auto precon_p = section->query("precon");
        if (precon_p.second)
        {
          if (precon_p.first == "none")
            precon = nullptr;
          else
            precon = create_scalar_solver<SystemLevelType_, TransferLevelType_, SolverVectorType_>(system_levels, transfer_levels, base, get_section_path(base, section, solver_name, precon_p.first));
        }

        if (solver_type == "pcg")
        {
          std::shared_ptr<Solver::PreconditionedIterativeSolver<SolverVectorType_> > solver;
          solver = Solver::new_pcg(derefer<SolverVectorType_>(system_levels.back()->matrix_sys, nullptr), derefer<SolverVectorType_>(system_levels.back()->filter_sys, nullptr), precon);
          configure_iterative_solver(section, solver);
          result = solver;
        }
        else if (solver_type == "bicgstab")
        {
          std::shared_ptr<Solver::PreconditionedIterativeSolver<SolverVectorType_> > solver;
          solver = Solver::new_bicgstab(derefer<SolverVectorType_>(system_levels.back()->matrix_sys, nullptr), derefer<SolverVectorType_>(system_levels.back()->filter_sys, nullptr), precon);
          configure_iterative_solver(section, solver);
          result = solver;
        }
        else if (solver_type == "fgmres")
        {
          auto krylov_dim_p = section->get_entry("krylov_dim");
          Index krylov_dim;
          if (krylov_dim_p.second)
            krylov_dim = std::stoul(krylov_dim_p.first);
          else
            throw InternalError(__func__, __FILE__, __LINE__, "no krylov_dim key found in fgmres section!");

          std::shared_ptr<Solver::PreconditionedIterativeSolver<SolverVectorType_> > solver;
          solver = Solver::new_fgmres(derefer<SolverVectorType_>(system_levels.back()->matrix_sys, nullptr), derefer<SolverVectorType_>(system_levels.back()->filter_sys, nullptr), krylov_dim, 0.0, precon);
          configure_iterative_solver(section, solver);
          result = solver;
        }
        else if (solver_type == "richardson")
        {
          auto omega_p = section->get_entry("omega");
          double omega;
          if (omega_p.second)
            omega = stod(omega_p.first);
          else
            omega = 1.0;

          std::shared_ptr<Solver::PreconditionedIterativeSolver<SolverVectorType_> > solver;
          solver = Solver::new_richardson(derefer<SolverVectorType_>(system_levels.back()->matrix_sys, nullptr), derefer<SolverVectorType_>(system_levels.back()->filter_sys, nullptr), omega, precon);
          configure_iterative_solver(section, solver);
          result = solver;
        }
        else if (solver_type == "jac")
        {
          auto omega_p = section->get_entry("omega");
          if (omega_p.second)
          {
            result = Solver::new_jacobi_precond(derefer<SolverVectorType_>(system_levels.back()->matrix_sys, nullptr), derefer<SolverVectorType_>(system_levels.back()->filter_sys, nullptr), std::stod(omega_p.first));
          }
          else
          {
            result = Solver::new_jacobi_precond(derefer<SolverVectorType_>(system_levels.back()->matrix_sys, nullptr), derefer<SolverVectorType_>(system_levels.back()->filter_sys, nullptr));
          }
        }
        else if (solver_type == "scale")
        {
          auto omega_p = section->get_entry("omega");
          if (omega_p.second)
          {
            result = Solver::new_scale_precond(derefer<SolverVectorType_>(system_levels.back()->filter_sys, nullptr), std::stod(omega_p.first));
          }
          else
            throw InternalError(__func__, __FILE__, __LINE__, "no omega key found in scale section!");
        }
        else if (solver_type == "ilu")
        {
          result = create_ilu_precon<SolverVectorType_>(system_levels, nullptr);
        }
        else if (solver_type == "ssor")
        {
          result = create_ssor_precon<SolverVectorType_>(system_levels, section, nullptr);
        }
        else if (solver_type == "mgv")
        {
          typename std::remove_reference<decltype(derefer<SolverVectorType_>(system_levels.front()->matrix_sys, nullptr))>::type dummy_sys;
          typename std::remove_reference<decltype(derefer<SolverVectorType_>(system_levels.front()->filter_sys, nullptr))>::type dummy_filter;
          typename std::remove_reference<decltype(derefer<SolverVectorType_>(transfer_levels.front()->prol_sys, nullptr))>::type dummy_prol;
          auto mgv = std::make_shared<
            Solver::BasicVCycle<
            decltype(dummy_sys),
            decltype(dummy_filter),
            decltype(dummy_prol)
          > >();

          std::shared_ptr<Solver::SolverBase<SolverVectorType_> > coarse_solver;
          auto coarse_solver_p = section->query("coarse");
          if (coarse_solver_p.second)
          {
            //artificial deque containing the coarsest level only
            typename std::remove_reference<decltype(system_levels)>::type coarse_system_level(system_levels.begin(), ++system_levels.begin());
            typename std::remove_reference<decltype(transfer_levels)>::type coarse_transfer_level(transfer_levels.begin(), ++transfer_levels.begin());
            coarse_solver = create_scalar_solver<SystemLevelType_, TransferLevelType_, SolverVectorType_>(coarse_system_level, coarse_transfer_level, base, get_section_path(base, section, solver_name, coarse_solver_p.first));
          }
          else
            throw InternalError(__func__, __FILE__, __LINE__, "mgv precon section without coarse key is not allowed!");

          mgv->set_coarse_level(derefer<SolverVectorType_>(system_levels.front()->matrix_sys, nullptr), derefer<SolverVectorType_>(system_levels.front()->filter_sys, nullptr), coarse_solver);

          auto jt = transfer_levels.begin();
          for (auto it = ++system_levels.begin(); it != system_levels.end(); ++it, ++jt)
          {
            std::shared_ptr<Solver::SolverBase<SolverVectorType_> > smoother;
            auto smoother_p = section->query("smoother");
            if (smoother_p.second)
            {
              //artificial deque containing all levels up to the current level, that shall be smoothed
              typename std::remove_reference<decltype(system_levels)>::type smoother_system_levels(system_levels.begin(), it+1);
              typename std::remove_reference<decltype(transfer_levels)>::type smoother_transfer_levels(transfer_levels.begin(), jt+1);
              smoother = create_scalar_solver<SystemLevelType_, TransferLevelType_, SolverVectorType_>(smoother_system_levels, smoother_transfer_levels, base, get_section_path(base, section, solver_name, smoother_p.first));
            }
            else
              throw InternalError(__func__, __FILE__, __LINE__, "mgv precon section without smoother key is not allowed!");

            mgv->push_level(derefer<SolverVectorType_>((*it)->matrix_sys, nullptr), derefer<SolverVectorType_>((*it)->filter_sys, nullptr), derefer<SolverVectorType_>((*jt)->prol_sys, nullptr),
            derefer<SolverVectorType_>((*jt)->rest_sys, nullptr), smoother, smoother);
          }

          result = mgv;
        }
        else if (solver_type == "schwarz")
        {
          result = create_schwarz_precon<SolverVectorType_>(system_levels, transfer_levels, base, solver_name, section, nullptr);
        }
        else
          throw InternalError(__func__, __FILE__, __LINE__, "solver with type " + solver_type + " unkown!");

        return result;
      }
    }; // SolverFactory
  } // namespace Control
} // namespace FEAST

#endif // CONTROL_SCALAR_BASIC_HPP
