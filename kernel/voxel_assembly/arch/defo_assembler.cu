#include <kernel/base_header.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/voxel_assembly/defo_assembler.hpp>
#include <kernel/lafem/matrix_gather_scatter_helper.hpp>
#include <kernel/voxel_assembly/helper/space_helper.hpp>
#include <kernel/util/tiny_algebra.hpp>

namespace FEAT
{
  namespace VoxelAssembly
  {
    namespace Kernel
    {
      /**************************************************************************************************************/
      /*                                       CUDA Kernel                                                          */
      /**************************************************************************************************************/

      template<typename Space_, typename DT_, typename IT_, FEAT::Intern::MatrixGatherScatterPolicy pol_ = FEAT::Intern::MatrixGatherScatterPolicy::useLocalOps>
      __global__ void defo_assembler_matrix1_bcsr(DT_* __restrict__ matrix_data,
                const IT_* __restrict__  matrix_row_ptr, const IT_* __restrict__ matrix_col_idx, Index matrix_num_rows, Index matrix_num_cols,
                const Tiny::Vector<DT_, Space_::world_dim>* __restrict__ cub_pt,
                const DT_* __restrict__  cub_wg, int num_cubs, DT_ alpha,
                const IT_* __restrict__ cell_to_dof, [[maybe_unused]] Index cell_num,
                const Tiny::Vector<DT_, Space_::world_dim>* __restrict__ nodes, [[maybe_unused]] Index node_size,
                const int* __restrict__ coloring_map, Index coloring_size, const IT_* __restrict__ cell_to_dof_sorter, DT_ nu)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if(idx >= coloring_size)
          return;
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
        //our local datatypes
        // typedef Tiny::Vector<DataType, dim> VecValueType;
        typedef Tiny::Matrix<DataType, dim, dim> MatValueType;
        // typedef Tiny::Vector<VecValueType, num_loc_dofs> LocalVectorType;
        typedef Tiny::Matrix<MatValueType, num_loc_dofs, num_loc_dofs> LocalMatrixType;

        // define local coefficients
        DataType local_coeffs[dim][num_loc_verts];
        //define local array
        LocalMatrixType loc_mat;


        //now do work for this cell
        Index cell = Index(coloring_map[idx]);
        const IndexSetWrapper<IndexType> local_dofs_w{cell_to_dof, IndexType(num_loc_dofs)};
        const IndexType* local_dofs = cell_to_dof + cell*num_loc_dofs;
        const IndexType* local_dof_sorter = cell_to_dof_sorter + cell*num_loc_dofs;
        SpaceHelp::set_coefficients(local_coeffs, local_dofs_w, nodes, cell);

        // To minimize code repetition, we call the same kernel from cuda and non cuda build...
        VoxelAssembly::Kernel::defo_assembly_kernel<SpaceHelp, LocalMatrixType>(loc_mat, local_coeffs, cub_pt, cub_wg, num_cubs, nu);

        //scatter
        LAFEM::template MatrixGatherScatterHelper<SpaceType, DataType, IndexType, pol_>::scatter_matrix_csr(loc_mat, (MatValueType*)matrix_data, local_dofs, local_dofs, IndexType(matrix_num_rows), IndexType(matrix_num_cols), matrix_row_ptr, matrix_col_idx, alpha, local_dof_sorter);

      }


      /**************************************************************************************************************/
      /*                                       CUDA Host OMP Kernels                                                */
      /**************************************************************************************************************/

      template<typename Space_, typename DT_, typename IT_, FEAT::Intern::MatrixGatherScatterPolicy pol_ = FEAT::Intern::MatrixGatherScatterPolicy::useLocalOps>
      void defo_assembler_matrix1_csr_host(DT_*  matrix_data,
                const IT_*  matrix_row_ptr, const IT_*  matrix_col_idx, Index matrix_num_rows, Index matrix_num_cols,
                const Tiny::Vector<DT_, Space_::world_dim>*  cub_pt,
                const DT_*  cub_wg, int num_cubs, DT_ alpha,
                const IT_*  cell_to_dof, [[maybe_unused]] Index cell_num,
                const Tiny::Vector<DT_, Space_::world_dim>*  nodes, [[maybe_unused]] Index node_size,
                const int*  coloring_map, Index coloring_size,
                const IT_* cell_to_dof_sorter, DT_ nu)
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
        //our local datatypes
        // typedef Tiny::Vector<DataType, dim> VecValueType;
        typedef Tiny::Matrix<DataType, dim, dim> MatValueType;
        // typedef Tiny::Vector<VecValueType, num_loc_dofs> LocalVectorType;
        typedef Tiny::Matrix<MatValueType, num_loc_dofs, num_loc_dofs> LocalMatrixType;



        #pragma omp parallel for
        for(Index idx = 0; idx < coloring_size; ++idx)
        {
          // define local coefficients
          DataType local_coeffs[dim][num_loc_verts];
          // and local matrix
          LocalMatrixType loc_mat;
          // now do work for this cell
          Index cell = Index(coloring_map[idx]);
          // std::cout << "Starting with cell " << cell << std::endl;
          const IndexSetWrapper<IndexType> local_dofs_w{cell_to_dof, IndexType(num_loc_dofs)};
          const IndexType* local_dofs = cell_to_dof + cell*num_loc_dofs;
          const IndexType* local_dof_sorter = cell_to_dof_sorter + cell*num_loc_dofs;

          SpaceHelp::set_coefficients(local_coeffs, local_dofs_w, nodes, cell);

          // To minimize code repetition, we call the same kernel from cuda and non cuda build...
          VoxelAssembly::Kernel::defo_assembly_kernel<SpaceHelp, LocalMatrixType>(loc_mat, local_coeffs, cub_pt, cub_wg, num_cubs, nu);

          // scatter
          LAFEM::template MatrixGatherScatterHelper<SpaceType, DataType, IndexType, pol_>::scatter_matrix_csr(loc_mat, (MatValueType*)matrix_data, local_dofs, local_dofs, IndexType(matrix_num_rows), IndexType(matrix_num_cols), matrix_row_ptr, matrix_col_idx, alpha, local_dof_sorter);
        }
      }
    }


    /**************************************************************************************************************/
    /*                                       CUDA kernel wrapper                                                  */
    /**************************************************************************************************************/

    namespace Arch
    {
      template<typename Space_, typename DT_, typename IT_>
      void assemble_defo_csr([[maybe_unused]] const Space_& space,
              const CSRMatrixData<DT_, IT_>& matrix_data,
              const AssemblyCubatureData<DT_>& cubature,
              const AssemblyMappingData<DT_, IT_>& dof_mapping,
              const std::vector<std::vector<int>>& coloring_maps_host,
              const std::vector<void*>& coloring_maps_device,
              DT_ alpha, DT_ nu
              )
      {
        const Index blocksize = Util::cuda_blocksize_blocked_assembly;
        // const Index blocksize = 64;

        for(Index i = 0; i < coloring_maps_device.size(); ++i)
        {
          dim3 grid;
          dim3 block;
          block.x = (unsigned int)blocksize;
          grid.x = (unsigned int)ceil(double(coloring_maps_host[i].size())/double(block.x));

          //kernel call, since this uses the standard stream, sync before next call is enforced:
          VoxelAssembly::Kernel::template defo_assembler_matrix1_bcsr<Space_, DT_, IT_, FEAT::Intern::MatrixGatherScatterPolicy::useLocalSortHelper><<< grid, block >>>(
              matrix_data.data, matrix_data.row_ptr, matrix_data.col_idx, matrix_data.num_rows, matrix_data.num_cols,
              (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
              cubature.cub_wg, cubature.num_cubs, alpha,
              dof_mapping.cell_to_dof, dof_mapping.cell_num,
              (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
              (const int*) coloring_maps_device[i], coloring_maps_host[i].size(), dof_mapping.cell_to_dof_sorter, nu
          );
        }

        //check for cuda error in our kernel
        Util::cuda_check_last_error();
      }

      template<typename Space_, typename DT_, typename IT_>
      void assemble_defo_csr_host([[maybe_unused]] const Space_& space,
              const CSRMatrixData<DT_, IT_>& matrix_data,
              const AssemblyCubatureData<DT_>& cubature,
              const AssemblyMappingData<DT_, IT_>& dof_mapping,
              const std::vector<std::vector<int>>& coloring_maps_host,
              [[maybe_unused]] const std::vector<void*>& coloring_maps_device,
              DT_ alpha, DT_ nu
              )
      {
        for(Index col = 0; col < Index(coloring_maps_host.size()); ++col)
        {
          VoxelAssembly::Kernel::template defo_assembler_matrix1_csr_host<Space_, DT_, IT_, FEAT::Intern::MatrixGatherScatterPolicy::useLocalSortHelper>(
            matrix_data.data, matrix_data.row_ptr, matrix_data.col_idx, matrix_data.num_rows, matrix_data.num_cols,
            (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
            cubature.cub_wg, cubature.num_cubs, alpha, dof_mapping.cell_to_dof, dof_mapping.cell_num,
            (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
            (const int*) coloring_maps_host[col].data(), coloring_maps_host[col].size(), dof_mapping.cell_to_dof_sorter, nu
          );
        }
      }
    }
  }
}


using namespace FEAT;
using namespace FEAT::VoxelAssembly;

template void Arch::assemble_defo_csr(const Q2StandardQuad&, const CSRMatrixData<double, std::uint32_t>&, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double, double);
template void Arch::assemble_defo_csr(const Q2StandardQuad&, const CSRMatrixData<float, std::uint32_t>&, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float, float);
template void Arch::assemble_defo_csr(const Q2StandardQuad&, const CSRMatrixData<double, std::uint64_t>&, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double, double);
template void Arch::assemble_defo_csr(const Q2StandardQuad&, const CSRMatrixData<float, std::uint64_t>&, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float, float);
#ifdef FEAT_HAVE_HALFMATH
template void Arch::assemble_defo_csr(const Q2StandardQuad&, const CSRMatrixData<Half, std::uint32_t>&, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, Half, Half);
template void Arch::assemble_defo_csr(const Q2StandardQuad&, const CSRMatrixData<Half, std::uint64_t>&, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, Half, Half);
#endif

template void Arch::assemble_defo_csr_host(const Q2StandardQuad&, const CSRMatrixData<double, std::uint32_t>&, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double, double);
template void Arch::assemble_defo_csr_host(const Q2StandardQuad&, const CSRMatrixData<float, std::uint32_t>&, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float, float);
template void Arch::assemble_defo_csr_host(const Q2StandardQuad&, const CSRMatrixData<double, std::uint64_t>&, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double, double);
template void Arch::assemble_defo_csr_host(const Q2StandardQuad&, const CSRMatrixData<float, std::uint64_t>&, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float, float);


template void Arch::assemble_defo_csr(const Q2StandardHexa&, const CSRMatrixData<double, std::uint32_t>&, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double, double);
template void Arch::assemble_defo_csr(const Q2StandardHexa&, const CSRMatrixData<float, std::uint32_t>&, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float, float);
template void Arch::assemble_defo_csr(const Q2StandardHexa&, const CSRMatrixData<double, std::uint64_t>&, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double, double);
template void Arch::assemble_defo_csr(const Q2StandardHexa&, const CSRMatrixData<float, std::uint64_t>&, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float, float);
#ifdef FEAT_HAVE_HALFMATH
template void Arch::assemble_defo_csr(const Q2StandardHexa&, const CSRMatrixData<Half, std::uint32_t>&, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, Half, Half);
template void Arch::assemble_defo_csr(const Q2StandardHexa&, const CSRMatrixData<Half, std::uint64_t>&, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, Half, Half);
#endif

template void Arch::assemble_defo_csr_host(const Q2StandardHexa&, const CSRMatrixData<double, std::uint32_t>&, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double, double);
template void Arch::assemble_defo_csr_host(const Q2StandardHexa&, const CSRMatrixData<float, std::uint32_t>&, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float, float);
template void Arch::assemble_defo_csr_host(const Q2StandardHexa&, const CSRMatrixData<double, std::uint64_t>&, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double, double);
template void Arch::assemble_defo_csr_host(const Q2StandardHexa&, const CSRMatrixData<float, std::uint64_t>&, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float, float);