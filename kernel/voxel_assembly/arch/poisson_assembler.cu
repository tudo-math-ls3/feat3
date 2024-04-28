#include <kernel/base_header.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/voxel_assembly/poisson_assembler.hpp>
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

      /**
      * \brief Cuda kernel poisson csr matrix assembly
      *
      * \tparam Space_ The underlying spacetype.
      * \tparam DT_ The datatype used.
      * \tparam IT_ The indextype.
      * \tparam pol_ The scatter and gather policy. See MatrixGatherScatterPolicy for details.
      *
      * \warning All ptr used here have to be unified or device pointer.
      * \param[out] matrix_data The matrix data array initilized as a unified memory array. Is added upon, so format before!
      * \param[in] matrix_row_ptr Row ptr of underlying csr matrix.
      * \param[in] matrix_col_idx Column index array of underlying csr matrix.
      * \param[in] matrix_num_rows Number of rows of csr matrix.
      * \param[in] matrix_num_cols Number of columns of csr matrix.
      * \param[in] cub_pt Array of cubature coordinate points.
      * \param[in] cub_wg Array of cubature weights.
      * \param[in] num_cubs Number of cubature points.
      * \param[in] alpha Scaling parameter for assembly.
      * \param[in] cell_to_dof Mapping cell to dof index.
      * \param[in] cell_num Number of cells.
      * \param[in] nodes Array of node points. Needed for local trafo assembly. /TODO: Version with precalced trafos.
      * \param[in] node_size Number of nodes.
      * \param[in] coloring_map Array of cell mapping for this color.
      * \param[in] coloring_size Size of the coloring array.
      * \param[in] cell_to_dof_sorter Depending on the used policy, this should either be nullptr, a local sorted array (or a col_ptr)
      */
      template<typename Space_, typename DT_, typename IT_, FEAT::Intern::MatrixGatherScatterPolicy pol_ = FEAT::Intern::MatrixGatherScatterPolicy::useLocalOps>
      __global__ void poisson_assembler_matrix1_csr(DT_* __restrict__ matrix_data,
                const IT_* __restrict__  matrix_row_ptr, const IT_* __restrict__ matrix_col_idx, Index matrix_num_rows, Index matrix_num_cols,
                const Tiny::Vector<DT_, Space_::world_dim>* __restrict__ cub_pt,
                const DT_* __restrict__  cub_wg, Index num_cubs, DT_ alpha,
                const IT_* __restrict__ cell_to_dof, [[maybe_unused]] Index cell_num,
                const Tiny::Vector<DT_, Space_::world_dim>* __restrict__ nodes, [[maybe_unused]] Index node_size,
                const int* __restrict__ coloring_map, Index coloring_size, const IT_* __restrict__ cell_to_dof_sorter)
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
        // Our local matrix type
        typedef Tiny::Matrix<DataType, num_loc_dofs, num_loc_dofs> LocMatType;

        // define local coeffs
        DataType local_coeffs[dim][num_loc_verts];

        // define local arrays
        Tiny::Matrix<DataType, num_loc_dofs, num_loc_dofs> loc_mat;

        DataType jac_det = DataType(0);

        typename SpaceHelp::EvalData basis_data;

        // now do work for this cell
        const Index cell = Index(coloring_map[idx]);
        const IndexType* local_dofs = cell_to_dof + cell*num_loc_dofs;
        const IndexType* local_dof_sorter = cell_to_dof_sorter + cell*num_loc_dofs;

        // printf("On cell %i My local dofs: ", cell);
        // for(int i = 0; i < num_loc_dofs; ++i)
        // {
        //   printf("%i ", *(local_dofs + i));
        // }
        // printf("\n");
        SpaceHelp::set_coefficients(local_coeffs, local_dofs, nodes);

        // To minimize code repetition, we call the same kernel from cuda and non cuda build...
        VoxelAssembly::Kernel::poisson_assembly_kernel<SpaceHelp, LocMatType>(loc_mat, local_coeffs, cub_pt, cub_wg, num_cubs);

        //scatter
        LAFEM::template MatrixGatherScatterHelper<SpaceType, DataType, IndexType, pol_>::scatter_matrix_csr(loc_mat, matrix_data, local_dofs, local_dofs, IndexType(matrix_num_rows), IndexType(matrix_num_cols), matrix_row_ptr, matrix_col_idx, alpha, local_dof_sorter);


      }

      /**************************************************************************************************************/
      /*                                       CUDA Host OMP Kernels                                                */
      /**************************************************************************************************************/
      /** \copydoc(poisson_assembler_matrix1_csr())*/
      template<typename Space_, typename DT_, typename IT_, FEAT::Intern::MatrixGatherScatterPolicy pol_ = FEAT::Intern::MatrixGatherScatterPolicy::useLocalOps>
      void poisson_assembler_matrix1_csr_host(DT_*  matrix_data,
                const IT_*  matrix_row_ptr, const IT_*  matrix_col_idx, Index matrix_num_rows, Index matrix_num_cols,
                const Tiny::Vector<DT_, Space_::world_dim>*  cub_pt,
                const DT_*  cub_wg, Index num_cubs, DT_ alpha,
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


        #pragma omp parallel for
        for(Index idx = 0; idx < coloring_size; ++idx)
        {
          // define local coefficients
          DataType local_coeffs[dim][num_loc_verts];
          // and local matrix
          Tiny::Matrix<DataType, num_loc_dofs, num_loc_dofs> loc_mat;
          // now do work for this cell
          Index cell = Index(coloring_map[idx]);
          // std::cout << "Starting with cell " << cell << std::endl;
          const IndexType* local_dofs = cell_to_dof + cell*num_loc_dofs;
          const IndexType* local_dof_sorter = cell_to_dof_sorter + cell*num_loc_dofs;

          SpaceHelp::set_coefficients(local_coeffs, local_dofs, nodes);

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
      void assemble_poisson_csr([[maybe_unused]] const Space_& space,
              const CSRMatrixData<DT_, IT_>& matrix_data,
              const AssemblyCubatureData<DT_>& cubature,
              const AssemblyMappingData<DT_, IT_>& dof_mapping,
              const std::vector<std::vector<int>>& coloring_maps_host,
              const std::vector<void*>& coloring_maps_device, DT_ alpha)
      {
        const Index blocksize = Util::cuda_blocksize_spmv;
        // const Index blocksize = 64;

        for(Index i = 0; i < coloring_maps_device.size(); ++i)
        {
          dim3 grid;
          dim3 block;
          block.x = (unsigned int)blocksize;
          grid.x = (unsigned int)ceil(double(coloring_maps_host[i].size())/double(block.x));

          //kernel call, since this uses the standard stream, sync before next call is enforced:
          VoxelAssembly::Kernel::template poisson_assembler_matrix1_csr<Space_, DT_, IT_, FEAT::Intern::MatrixGatherScatterPolicy::useLocalSortHelper><<< grid, block >>>(
              matrix_data.data, matrix_data.row_ptr, matrix_data.col_idx, matrix_data.num_rows, matrix_data.num_cols,
              (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
              cubature.cub_wg, cubature.num_cubs, alpha,
              dof_mapping.cell_to_dof, dof_mapping.cell_num,
              (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
              (const int*) coloring_maps_device[i], coloring_maps_host[i].size(), dof_mapping.cell_to_dof_sorter
          );
        }

        //check for cuda error in our kernel
        Util::cuda_check_last_error();
      }

      template<typename Space_, typename DT_, typename IT_>
      void assemble_poisson_csr_host([[maybe_unused]] const Space_& space,
              const CSRMatrixData<DT_, IT_>& matrix_data,
              const AssemblyCubatureData<DT_>& cubature,
              const AssemblyMappingData<DT_, IT_>& dof_mapping,
              const std::vector<std::vector<int>>& coloring_maps_host,
              [[maybe_unused]] const std::vector<void*>& coloring_maps_device, DT_ alpha)
      {
        for(Index col = 0; col < Index(coloring_maps_host.size()); ++col)
        {
          VoxelAssembly::Kernel::template poisson_assembler_matrix1_csr_host<Space_, DT_, IT_, FEAT::Intern::MatrixGatherScatterPolicy::useLocalSortHelper>(
            matrix_data.data, matrix_data.row_ptr, matrix_data.col_idx, matrix_data.num_rows, matrix_data.num_cols,
            (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
            cubature.cub_wg, cubature.num_cubs, alpha, dof_mapping.cell_to_dof, dof_mapping.cell_num,
            (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
            (const int*) coloring_maps_host[col].data(), coloring_maps_host[col].size(), dof_mapping.cell_to_dof_sorter
          );
        }
      }
    }
  }
}

using namespace FEAT;
using namespace FEAT::VoxelAssembly;

/*--------------------Poisson Assembler Q2Quad-------------------------------------------------*/
template void Arch::assemble_poisson_csr(const Q2StandardQuad&, const CSRMatrixData<double, std::uint32_t>&, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double);
template void Arch::assemble_poisson_csr(const Q2StandardQuad&, const CSRMatrixData<float, std::uint32_t>&, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float);
template void Arch::assemble_poisson_csr(const Q2StandardQuad&, const CSRMatrixData<double, std::uint64_t>&, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double);
template void Arch::assemble_poisson_csr(const Q2StandardQuad&, const CSRMatrixData<float, std::uint64_t>&, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float);
#ifdef FEAT_HAVE_HALFMATH
template void Arch::assemble_poisson_csr(const Q2StandardQuad&, const CSRMatrixData<Half, std::uint32_t>&, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, Half);
template void Arch::assemble_poisson_csr(const Q2StandardQuad&, const CSRMatrixData<Half, std::uint64_t>&, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, Half);
#endif

template void Arch::assemble_poisson_csr_host(const Q2StandardQuad&, const CSRMatrixData<double, std::uint32_t>&, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double);
template void Arch::assemble_poisson_csr_host(const Q2StandardQuad&, const CSRMatrixData<float, std::uint32_t>&, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float);
template void Arch::assemble_poisson_csr_host(const Q2StandardQuad&, const CSRMatrixData<double, std::uint64_t>&, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double);
template void Arch::assemble_poisson_csr_host(const Q2StandardQuad&, const CSRMatrixData<float, std::uint64_t>&, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float);

/*---------------------Poisson Assembler Q2Hexa------------------------------------------------------*/
template void Arch::assemble_poisson_csr(const Q2StandardHexa&, const CSRMatrixData<double, std::uint32_t>&, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double);
template void Arch::assemble_poisson_csr(const Q2StandardHexa&, const CSRMatrixData<float, std::uint32_t>&, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float);
template void Arch::assemble_poisson_csr(const Q2StandardHexa&, const CSRMatrixData<double, std::uint64_t>&, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double);
template void Arch::assemble_poisson_csr(const Q2StandardHexa&, const CSRMatrixData<float, std::uint64_t>&, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float);
#ifdef FEAT_HAVE_HALFMATH
template void Arch::assemble_poisson_csr(const Q2StandardHexa&, const CSRMatrixData<Half, std::uint32_t>&, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, Half);
template void Arch::assemble_poisson_csr(const Q2StandardHexa&, const CSRMatrixData<Half, std::uint64_t>&, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, Half);
#endif

template void Arch::assemble_poisson_csr_host(const Q2StandardHexa&, const CSRMatrixData<double, std::uint32_t>&, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double);
template void Arch::assemble_poisson_csr_host(const Q2StandardHexa&, const CSRMatrixData<float, std::uint32_t>&, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float);
template void Arch::assemble_poisson_csr_host(const Q2StandardHexa&, const CSRMatrixData<double, std::uint64_t>&, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double);
template void Arch::assemble_poisson_csr_host(const Q2StandardHexa&, const CSRMatrixData<float, std::uint64_t>&, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float);