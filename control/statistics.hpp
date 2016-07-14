#pragma once
#ifndef CONTROL_STATISTICS_HPP
#define CONTROL_STATISTICS_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/comm_base.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/util/memory_usage.hpp>

namespace FEAT
{
  namespace Control
  {
    class Statistics
    {
      public:

        /**
         * \brief Create a detailed report about the application execution
         *
         * \param solver_toe The execution time of the complete linear solver in seconds
         * \param statistics_check The result of args.check("statistics"), i.e. the detail level(0 = some details, 1 = many details and log file output)
         * \param shape_dimension The maximum dimension of the used shapes
         * \param system_levels The system level hierarchy
         * \param transfer_levels The transfer level hierarchy
         * \param solver The solver tree
         * \param domain The domain
         *
         */
        template <typename SystemLevelType_, typename TransferLevelType_, typename DomainType_, typename SolverType_>
        static void report(double solver_toe, int statistics_check, int shape_dimension,
        std::deque<SystemLevelType_*> & system_levels,
        std::deque<TransferLevelType_*> & transfer_levels,
        SolverType_ & solver,
        DomainType_ & domain)
        {
          Index rank(Util::Comm::rank());
          Index nranks(Util::Comm::size());

          FEAT::Statistics::toe_solve = solver_toe;

          std::size_t la_size(0);
          std::for_each(system_levels.begin(), system_levels.end(), [&] (SystemLevelType_ * n) { la_size += n->bytes(); });
          std::for_each(transfer_levels.begin(), transfer_levels.end(), [&] (TransferLevelType_ * n) { la_size += n->bytes(); });
          std::size_t mpi_size(0);
          std::for_each(system_levels.begin(), system_levels.end(), [&] (SystemLevelType_ * n) { mpi_size += n->gate_sys.bytes(); });

          String op_timings = FEAT::Statistics::get_formatted_times(solver_toe);

          Index cells_coarse_local = domain.get_levels().front()->get_mesh().get_num_entities(shape_dimension);
          Index cells_coarse_max;
          Index cells_coarse_min;
          Util::Comm::allreduce(&cells_coarse_local, &cells_coarse_max, 1, Util::CommOperationMax());
          Util::Comm::allreduce(&cells_coarse_local, &cells_coarse_min, 1, Util::CommOperationMin());
          Index cells_fine_local = domain.get_levels().back()->get_mesh().get_num_entities(shape_dimension);
          Index cells_fine_max;
          Index cells_fine_min;
          Util::Comm::allreduce(&cells_fine_local, &cells_fine_max, 1, Util::CommOperationMax());
          Util::Comm::allreduce(&cells_fine_local, &cells_fine_min, 1, Util::CommOperationMin());

          Index dofs_coarse_local = system_levels.front()->matrix_sys.columns();
          Index dofs_coarse_max;
          Index dofs_coarse_min;
          Util::Comm::allreduce(&dofs_coarse_local, &dofs_coarse_max, 1, Util::CommOperationMax());
          Util::Comm::allreduce(&dofs_coarse_local, &dofs_coarse_min, 1, Util::CommOperationMin());
          Index dofs_fine_local = system_levels.back()->matrix_sys.columns();
          Index dofs_fine_max;
          Index dofs_fine_min;
          Util::Comm::allreduce(&dofs_fine_local, &dofs_fine_max, 1, Util::CommOperationMax());
          Util::Comm::allreduce(&dofs_fine_local, &dofs_fine_min, 1, Util::CommOperationMin());

          Index nzes_coarse_local = system_levels.front()->matrix_sys.used_elements();
          Index nzes_coarse_max;
          Index nzes_coarse_min;
          Util::Comm::allreduce(&nzes_coarse_local, &nzes_coarse_max, 1, Util::CommOperationMax());
          Util::Comm::allreduce(&nzes_coarse_local, &nzes_coarse_min, 1, Util::CommOperationMin());
          Index nzes_fine_local = system_levels.back()->matrix_sys.used_elements();
          Index nzes_fine_max;
          Index nzes_fine_min;
          Util::Comm::allreduce(&nzes_fine_local, &nzes_fine_max, 1, Util::CommOperationMax());
          Util::Comm::allreduce(&nzes_fine_local, &nzes_fine_min, 1, Util::CommOperationMin());

          if (rank == 0 && statistics_check >= 0)
          {
            std::cout<< std::endl;
            std::cout << String("TOE partition:").pad_back(20) << FEAT::Statistics::toe_partition << std::endl;
            std::cout << String("TOE assembly:").pad_back(20) << FEAT::Statistics::toe_assembly << std::endl;
            std::cout << String("TOE solve:").pad_back(20) << FEAT::Statistics::toe_solve << std::endl;
            std::cout << std::endl << solver->get_formatted_solver_tree().trim() << std::endl;
            String flops = FEAT::Statistics::get_formatted_flops(solver_toe, nranks);
            std::cout<<flops<<std::endl<<std::endl;
            std::cout<<op_timings<<std::endl<<std::endl;
            std::cout<<String("Domain size:").pad_back(20) << double(domain.bytes())  / (1024. * 1024.)  << " MByte" << std::endl;
            std::cout<<String("MPI size:").pad_back(20) << double(mpi_size) / (1024. * 1024.) << " MByte" << std::endl;
            std::cout<<String("LA size:").pad_back(20) << double(la_size) / (1024. * 1024.) << " MByte" << std::endl << std::endl;
            std::cout<<Util::get_formatted_memory_usage()<<std::endl;
            std::cout<<String("#Mesh cells:").pad_back(20) << "coarse " << cells_coarse_max << "/" << cells_coarse_min << ", fine " << cells_fine_max << "/" << cells_fine_min << std::endl;
            std::cout<<String("#DOFs:").pad_back(20) << "coarse " << dofs_coarse_max << "/" << dofs_coarse_min << ", fine " << dofs_fine_max << "/" << dofs_fine_min << std::endl;
            std::cout<<String("#NZEs").pad_back(20) << "coarse " << nzes_coarse_max << "/" << nzes_coarse_min << ", fine " << nzes_fine_max << "/" << nzes_fine_min << std::endl;
            std::cout<<std::endl;
            /*if (statistics_check > 0) // provided parameter full or whatever
            {
              std::cout<<FEAT::Statistics::get_formatted_solvers();
            }*/
          }
          if (statistics_check > 0) // provided parameter full or whatever
          {
            FEAT::Statistics::write_out_solver_statistics_scheduled(rank, la_size, domain.bytes(), mpi_size,
            domain.get_levels().back()->get_mesh().get_num_entities(shape_dimension), system_levels.back()->matrix_sys.columns(),
            system_levels.back()->matrix_sys.used_elements());
          }
        }

    }; // StatisticsControl
  } // namespace Control
} // namespace FEAT

#endif // CONTROL_STATISTICS_HPP
