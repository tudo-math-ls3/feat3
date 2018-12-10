#pragma once
#ifndef CONTROL_STATISTICS_HPP
#define CONTROL_STATISTICS_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/dist.hpp>
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
         * \param domain The domain
         *
         */
        template <typename SystemLevelType_, typename DomainType_>
        static void report(double solver_toe, int statistics_check, int DOXY(shape_dimension),
          std::deque<std::shared_ptr<SystemLevelType_>> & system_levels,
          DomainType_ & domain)
        {
          Dist::Comm comm(Dist::Comm::world());
          int rank(comm.rank());
          int nranks(comm.size());

          FEAT::Statistics::toe_solve = solver_toe;

          std::size_t la_size(0);
          std::size_t mpi_size(0);
          for(const auto& sl : system_levels)
          {
            la_size += sl->bytes();
            mpi_size += sl->gate_sys.bytes() + sl->coarse_muxer_sys.bytes();
          }

          String op_timings = FEAT::Statistics::get_formatted_times(solver_toe);

          Index cells_coarse_local = domain.back()->get_mesh().get_num_elements();
          Index cells_coarse_max;
          Index cells_coarse_min;
          comm.allreduce(&cells_coarse_local, &cells_coarse_max, std::size_t(1), Dist::op_max);
          comm.allreduce(&cells_coarse_local, &cells_coarse_min, std::size_t(1), Dist::op_min);
          Index cells_fine_local = domain.front()->get_mesh().get_num_elements();
          Index cells_fine_max;
          Index cells_fine_min;
          comm.allreduce(&cells_fine_local, &cells_fine_max, std::size_t(1), Dist::op_max);
          comm.allreduce(&cells_fine_local, &cells_fine_min, std::size_t(1), Dist::op_min);

          Index dofs_coarse_local = system_levels.back()->matrix_sys.local().columns();
          Index dofs_coarse_max;
          Index dofs_coarse_min;
          comm.allreduce(&dofs_coarse_local, &dofs_coarse_max, std::size_t(1), Dist::op_max);
          comm.allreduce(&dofs_coarse_local, &dofs_coarse_min, std::size_t(1), Dist::op_min);
          Index dofs_fine_local = system_levels.front()->matrix_sys.local().columns();
          Index dofs_fine_max;
          Index dofs_fine_min;
          comm.allreduce(&dofs_fine_local, &dofs_fine_max, std::size_t(1), Dist::op_max);
          comm.allreduce(&dofs_fine_local, &dofs_fine_min, std::size_t(1), Dist::op_min);

          Index nzes_coarse_local = system_levels.back()->matrix_sys.local().used_elements();
          Index nzes_coarse_max;
          Index nzes_coarse_min;
          comm.allreduce(&nzes_coarse_local, &nzes_coarse_max, std::size_t(1), Dist::op_max);
          comm.allreduce(&nzes_coarse_local, &nzes_coarse_min, std::size_t(1), Dist::op_min);
          Index nzes_fine_local = system_levels.front()->matrix_sys.local().used_elements();
          Index nzes_fine_max;
          Index nzes_fine_min;
          comm.allreduce(&nzes_fine_local, &nzes_fine_max, std::size_t(1), Dist::op_max);
          comm.allreduce(&nzes_fine_local, &nzes_fine_min, std::size_t(1), Dist::op_min);

          if (rank == 0 && statistics_check >= 0)
          {
            std::cout<< std::endl;
            std::cout << String("TOE partition:").pad_back(20) << FEAT::Statistics::toe_partition << std::endl;
            std::cout << String("TOE assembly:").pad_back(20) << FEAT::Statistics::toe_assembly << std::endl;
            std::cout << String("TOE solve:").pad_back(20) << FEAT::Statistics::toe_solve << std::endl;
            std::cout << std::endl << FEAT::Statistics::get_formatted_solver_tree().trim() <<std::endl;

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
          }
          if (statistics_check > 0) // provided parameter full or whatever
          {
            /// \todo reimplement method based on expressions
            /*FEAT::Statistics::write_out_solver_statistics_scheduled(rank, la_size, domain.bytes(), mpi_size,
            domain.get_levels().back()->get_mesh().get_num_entities(shape_dimension), (*system_levels.back()->matrix_sys).columns(),
            (*system_levels.back()->matrix_sys).used_elements());*/
          }
        }

    }; // StatisticsControl
  } // namespace Control
} // namespace FEAT

#endif // CONTROL_STATISTICS_HPP
