// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/voxel_assembly/poisson_assembler.hpp>
#include <kernel/lafem/matrix_gather_scatter_helper.hpp>

namespace FEAT
{
  namespace VoxelAssembly
  {
    namespace Kernel
    {

      template<typename Space_, typename DT_, typename IT_, FEAT::Intern::MatrixGatherScatterPolicy pol_ = FEAT::Intern::MatrixGatherScatterPolicy::useLocalOps>
      void poisson_assembler_matrix1_csr_host(DT_*  matrix_data,
                const IT_*  matrix_row_ptr, const IT_*  matrix_col_idx, Index matrix_num_rows, Index matrix_num_cols,
                const Tiny::Vector<DT_, Space_::world_dim>*  cub_pt,
                const DT_*  cub_wg, int num_cubs, DT_ alpha,
                const IT_*  cell_to_dof, [[maybe_unused]] Index cell_num,
                const Tiny::Vector<DT_, Space_::world_dim>*  nodes, [[maybe_unused]] Index node_size,
                const int*  coloring_map, Index coloring_size,
                const IT_* cell_to_dof_sorter)
      {
        // define types
        typedef Space_ SpaceType;
        typedef DT_ DataType;
        typedef IT_ IndexType;

        typedef VoxelAssembly::SpaceHelper<SpaceType, DT_, IT_> SpaceHelp;

        constexpr int dim = SpaceHelp::dim;

        // define local sizes
        constexpr int num_loc_dofs = SpaceType::DofMappingType::dof_count;
        // get number of nodes per element
        constexpr int num_loc_verts = SpaceType::MeshType::template IndexSet<dim, 0>::Type::num_indices;
        // Our local matrix type
        typedef Tiny::Matrix<DataType, num_loc_dofs, num_loc_dofs> LocMatType;


        FEAT_PRAGMA_OMP(parallel for)
        for(Index idx = 0; idx < coloring_size; ++idx)
        {
          // define local coefficients
          Tiny::Matrix<DataType, dim, num_loc_verts> local_coeffs;
          // and local matrix
          Tiny::Matrix<DataType, num_loc_dofs, num_loc_dofs> loc_mat;
          // now do work for this cell
          IndexType cell = IndexType(coloring_map[idx]);
          // std::cout << "Starting with cell " << cell << "\n";
          const IndexSetWrapper<IndexType> local_dofs_w{cell_to_dof, IndexType(num_loc_dofs)};
          const IndexType* local_dofs = cell_to_dof + cell*num_loc_dofs;
          const IndexType* local_dof_sorter = cell_to_dof_sorter + cell*num_loc_dofs;

          SpaceHelp::set_coefficients(local_coeffs, local_dofs_w, nodes, cell);

          // To minimize code repetition, we call the same kernel from cuda and non cuda build...
          VoxelAssembly::Kernel::poisson_assembly_kernel<SpaceHelp, LocMatType>(loc_mat, local_coeffs, cub_pt, cub_wg, num_cubs);

          // scatter
          LAFEM::template MatrixGatherScatterHelper<SpaceType, DataType, IndexType, pol_>::scatter_matrix_csr(loc_mat, matrix_data, local_dofs, local_dofs, IndexType(matrix_num_rows), IndexType(matrix_num_cols), matrix_row_ptr, matrix_col_idx, alpha, local_dof_sorter);
        }


      }
    }

    namespace Arch
    {
      template<typename Space_, typename DT_, typename IT_>
      void assemble_poisson_csr_host([[maybe_unused]] const Space_& space,
              const CSRMatrixData<DT_, IT_>& matrix_data,
              const AssemblyCubatureData<DT_>& cubature,
              const AssemblyMappingData<DT_, IT_>& dof_mapping,
              const std::vector<int*>& coloring_maps,
              const std::vector<Index>& coloring_map_sizes,
              DT_ alpha)
      {
        for(Index col = 0; col < Index(coloring_maps.size()); ++col)
        {
          VoxelAssembly::Kernel::template poisson_assembler_matrix1_csr_host<Space_, DT_, IT_, FEAT::Intern::MatrixGatherScatterPolicy::useLocalSortHelper>(
            matrix_data.data, matrix_data.row_ptr, matrix_data.col_idx, matrix_data.num_rows, matrix_data.num_cols,
            (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
            cubature.cub_wg, cubature.num_cubs, alpha, dof_mapping.cell_to_dof, dof_mapping.cell_num,
            (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
            (const int*) coloring_maps[col], coloring_map_sizes[col], dof_mapping.cell_to_dof_sorter
          );
        }
      }
    }
  }
}

//instantiate the templates
using namespace FEAT;
using namespace FEAT::VoxelAssembly;


/*--------------------Poisson Assembler Q2Quad-------------------------------------------------*/
template void Arch::assemble_poisson_csr_host(const Q2StandardQuad&, const CSRMatrixData<double, std::uint32_t>&, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, double);
template void Arch::assemble_poisson_csr_host(const Q2StandardQuad&, const CSRMatrixData<float, std::uint32_t>&, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, float);
template void Arch::assemble_poisson_csr_host(const Q2StandardQuad&, const CSRMatrixData<double, std::uint64_t>&, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, double);
template void Arch::assemble_poisson_csr_host(const Q2StandardQuad&, const CSRMatrixData<float, std::uint64_t>&, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, float);


/*---------------------Poisson Assembler Q2Hexa------------------------------------------------------*/
template void Arch::assemble_poisson_csr_host(const Q2StandardHexa&, const CSRMatrixData<double, std::uint32_t>&, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, double);
template void Arch::assemble_poisson_csr_host(const Q2StandardHexa&, const CSRMatrixData<float, std::uint32_t>&, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, float);
template void Arch::assemble_poisson_csr_host(const Q2StandardHexa&, const CSRMatrixData<double, std::uint64_t>&, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, double);
template void Arch::assemble_poisson_csr_host(const Q2StandardHexa&, const CSRMatrixData<float, std::uint64_t>&, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, float);
