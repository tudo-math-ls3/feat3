#include <kernel/base_header.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/voxel_assembly/burgers_velo_material_assembler.hpp>
#include <kernel/lafem/matrix_gather_scatter_helper.hpp>
#include <kernel/lafem/vector_gather_scatter_helper.hpp>
#include <kernel/voxel_assembly/helper/space_helper.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <kernel/util/cuda_math.cuh>

#include "omp.h"
#include "assert.h"

#include <numeric>


namespace FEAT
{
  namespace VoxelAssembly
  {
    namespace Kernel
    {

      /**************************************************************************************************************/
      /*                                       CUDA Kernel                                                          */
      /**************************************************************************************************************/

      template<typename Space_, typename DT_, typename IT_, typename ViscoFunc_, typename ViscoDFunc_, bool stream_diff_, FEAT::Intern::MatrixGatherScatterPolicy pol_ = FEAT::Intern::MatrixGatherScatterPolicy::useLocalOps>
      __global__ void full_burgers_vm_assembler_matrix1_bcsr(DT_* __restrict__ matrix_data, const DT_* __restrict__ conv_data,
                const IT_* __restrict__  matrix_row_ptr, const IT_* __restrict__ matrix_col_idx, Index matrix_num_rows, Index matrix_num_cols,
                const Tiny::Vector<DT_, Space_::world_dim>* __restrict__ cub_pt,
                const DT_* __restrict__  cub_wg, int num_cubs, DT_ alpha,
                const IT_* __restrict__ cell_to_dof, [[maybe_unused]] Index cell_num,
                const Tiny::Vector<DT_, Space_::world_dim>* __restrict__ nodes, [[maybe_unused]] Index node_size,
                const int* __restrict__ coloring_map, Index coloring_size, const IT_* __restrict__ cell_to_dof_sorter,
                const VoxelAssembly::AssemblyBurgersData<DT_> burgers_params, const VoxelAssembly::AssemblyMaterialData<DT_> material_params,
                const ViscoFunc_ visco_func, const ViscoDFunc_ visco_d_func)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if(idx >= coloring_size)
          return;
        //define types
        typedef Space_ SpaceType;
        typedef DT_ DataType;
        typedef IT_ IndexType;

        const DataType& beta{burgers_params.beta};
        const DataType& sd_delta{burgers_params.sd_delta};
        const DataType& sd_v_norm{burgers_params.sd_v_norm};

        const DataType tol_eps = CudaMath::cuda_get_sqrt_eps<DataType>();
        const bool need_streamline = (CudaMath::cuda_abs(sd_delta) > DataType(0)) && (sd_v_norm > tol_eps);
        const bool need_convection = CudaMath::cuda_abs(beta) > DataType(0);

        typedef VoxelAssembly::SpaceHelper<SpaceType, DT_, IT_> SpaceHelp;

        constexpr int dim = SpaceHelp::dim;

        //define local sizes
        constexpr int num_loc_dofs = SpaceType::DofMappingType::dof_count;
        //get number of nodes per element
        constexpr int num_loc_verts = SpaceType::MeshType::template IndexSet<dim, 0>::Type::num_indices;

        //our local datatypes
        typedef Tiny::Vector<DataType, dim> VecValueType;
        typedef Tiny::Matrix<DataType, dim, dim> MatValueType;
        typedef Tiny::Vector<VecValueType, num_loc_dofs> LocalVectorType;
        typedef Tiny::Matrix<MatValueType, num_loc_dofs, num_loc_dofs> LocalMatrixType;


        // define local coefficients
          Tiny::Matrix<DataType, dim, num_loc_verts> local_coeffs;
        // and local matrix
        LocalMatrixType loc_mat;
        LocalVectorType local_conv_dofs(DataType(0));
        // now do work for this cell
        Index cell = Index(coloring_map[idx]);
        // std::cout << "Starting with cell " << cell << std::endl;
        const IndexSetWrapper<IndexType> local_dofs_w{cell_to_dof, IndexType(num_loc_dofs)};
        const IndexType* local_dofs = cell_to_dof + cell*num_loc_dofs;
        const IndexType* local_dof_sorter = cell_to_dof_sorter + cell*num_loc_dofs;

        SpaceHelp::set_coefficients(local_coeffs, local_dofs_w, nodes, cell);
        //if we need to, gather local convection vector
        LAFEM::template VectorGatherScatterHelper<SpaceType, DataType, IndexType>::gather_vector_dense(local_conv_dofs,
                  (const VecValueType*)conv_data, IndexType(matrix_num_rows), local_dofs,DataType(1));

            VoxelAssembly::Kernel::burgers_velo_material_mat_assembly_kernel<SpaceHelp, LocalMatrixType, LocalVectorType, dim, num_loc_verts, ViscoFunc_, ViscoDFunc_, stream_diff_>(
                            loc_mat, local_conv_dofs, local_coeffs, cub_pt, cub_wg, num_cubs, burgers_params, material_params, need_convection, tol_eps, visco_func, visco_d_func);

        //scatter
        LAFEM::template MatrixGatherScatterHelper<SpaceType, DataType, IndexType, pol_>::scatter_matrix_csr(loc_mat, (MatValueType*)matrix_data, local_dofs, local_dofs, IndexType(matrix_num_rows), IndexType(matrix_num_cols), matrix_row_ptr, matrix_col_idx, alpha, local_dof_sorter);

      }

      template<typename Space_, typename DT_, typename IT_, typename ViscoFunc_, typename ViscoDFunc_, bool stream_diff_, FEAT::Intern::MatrixGatherScatterPolicy pol_ = FEAT::Intern::MatrixGatherScatterPolicy::useLocalOps>
      __global__ void full_burgers_vm_assembler_matrix1_bcsr_warp_based(DT_* __restrict__ matrix_data, const DT_* __restrict__ conv_data,
                const IT_* __restrict__  matrix_row_ptr, const IT_* __restrict__ matrix_col_idx, Index matrix_num_rows, Index matrix_num_cols,
                const Tiny::Vector<DT_, Space_::world_dim>* __restrict__ cub_pt,
                const DT_* __restrict__  cub_wg, int num_cubs, DT_ alpha,
                const IT_* __restrict__ cell_to_dof, [[maybe_unused]] Index cell_num,
                const Tiny::Vector<DT_, Space_::world_dim>* __restrict__ nodes, [[maybe_unused]] Index node_size,
                const int* __restrict__ coloring_map, Index coloring_size, const IT_* __restrict__/* cell_to_dof_sorter*/,
                const VoxelAssembly::AssemblyBurgersData<DT_> burgers_params, const VoxelAssembly::AssemblyMaterialData<DT_> material_params,
                const ViscoFunc_ visco_func, const ViscoDFunc_ visco_d_func, const int loc_block_size)
      {
        // this is a warp based/block based approach
        // The number of threads calculating one element depend on the chosen block size!
        // Hence you should at least take a multiple of 32 as blocksize
        // the base idea is load node ptrs, local csr ptr and so on into shared arrays
        // but first of all, calculate thread idx, the block index and the local (warp intern) idx
        // get cooperative group
        cg::thread_block tb = cg::this_thread_block();
        const Index t_idx = tb.thread_rank();

        // for now, we will use 3 threads per local matrix row, to optimize for 3D kernels
        // define a few constexpr variables which decide on the size of our local arrays
        //define types
        typedef Space_ SpaceType;
        typedef DT_ DataType;
        typedef IT_ IndexType;

        typedef VoxelAssembly::SpaceHelper<SpaceType, DT_, IT_> SpaceHelp;

        constexpr int dim = SpaceHelp::dim;

        //define local sizes
        constexpr int num_loc_dofs = SpaceType::DofMappingType::dof_count;
        //get number of nodes per element
        constexpr int num_loc_verts = SpaceType::MeshType::template IndexSet<dim, 0>::Type::num_indices;

        // define our shared array, into which we define our local array by offsets
        // for now, this does not necessarily be allocated dynamiclly, since the local sizes are all static
        // so, its not necessary that the size of this has to be sufficiently initilized at kernel launch time
        // constexpr int shared_size = dim * num_loc_verts + num_loc_dofs*(num_loc_dofs+3);
        // __shared__ __align__(CudaMath::cuda_get_min_align<DataType, IndexType>()) unsigned char shared_array[shared_size];
        // we need a shared array to hold our vertices
        // define local coefficients as shared array
        // Our static size shared arrays <- this works only when we use whole thread group and not warp groups...
        // __shared__ DataType local_coeffs[dim*num_loc_verts];
        // __shared__ IndexType local_dofs[num_loc_dofs];
        // // __shared__ IndexType local_dofs_sorter[num_loc_dofs];
        // __shared__ DataType loc_mat[num_loc_dofs*dim*num_loc_dofs*dim];
        // __shared__ DataType local_conv_dofs[num_loc_dofs*dim];
        // with external array, we have to be extremly careful with alignment, for this reason, first index the DataTypes (except local mat, since this will require special handling)
        extern __shared__ unsigned char shared_array[];

        typedef BurgersSharedDataGlobalWrapper<SpaceHelp> BSDGWrapper;
        typedef BurgersMatSharedDataKernelWrapper<SpaceHelp, stream_diff_> BSDKWrapper;

        BSDGWrapper* shared_meta_wrapper =
            reinterpret_cast<BSDGWrapper*>(&shared_array[0]);
        // first define local_coeffs, since this will be cast to a tiny matrix, which is required to be aligned...
        DataType* local_coeffs = &(shared_meta_wrapper->local_coeffs[0]);
        DataType* local_conv_dofs = &(shared_meta_wrapper->local_conv_dofs[0]);
        IndexType* local_dofs = &(shared_meta_wrapper->local_dofs[0]);
        // dummy ptr, since we do not use this...
        IndexType* local_dofs_sorter = nullptr;

        // these can be defined in shared memory
        DataType& tol_eps = shared_meta_wrapper->tol_eps;
        bool& need_streamline = shared_meta_wrapper->need_streamline;
        bool& need_convection = shared_meta_wrapper->need_convection;

        BSDKWrapper* shared_kernel_wrapper = (BSDKWrapper*)(shared_meta_wrapper+1);
        // TODO: Also gather node coefficients into cell local arrays? <- Increased data size, but better access?
        // DataType* loc_nodes = (DataType*)&local_conv_dofs[std::size_t(num_loc_dofs)];
        // offset of shared data in kernel is...
        DataType* loc_mat = VoxelAssembly::get_aligned_address<DataType>((void*)(&shared_array[0]+sizeof(BSDGWrapper)+sizeof(BSDKWrapper)));

        //stride based for loop
        for(int b_idx = blockIdx.x; b_idx < int(coloring_size); b_idx += gridDim.x)
        {
          //get our actual cell index
          const Index cell = Index(coloring_map[b_idx]);

          // start the memcpy calls before formating to overlap dataloading and calculation
          cg::memcpy_async(tb, local_dofs, cell_to_dof+num_loc_dofs*cell, num_loc_dofs*sizeof(IndexType));
          // cg::memcpy_async(tb, local_dofs_sorter, cell_to_dof_sorter+num_loc_dofs*cell, num_loc_dofs*sizeof(IndexType));
          // cg::memcpy_async(tb, loc_nodes, (DT_*)(nodes+cell*num_loc_verts), num_loc_verts); // does not work yet, due to mesh nodes


          VoxelAssembly::coalesced_format(tb, local_coeffs, dim*num_loc_verts);
          VoxelAssembly::coalesced_format(tb, loc_mat, loc_block_size*dim*dim);
          VoxelAssembly::coalesced_format(tb, local_conv_dofs, num_loc_dofs*dim);



          const DataType& beta{burgers_params.beta};
          const DataType& sd_delta{burgers_params.sd_delta};
          const DataType& sd_v_norm{burgers_params.sd_v_norm};

          //let the first thread rank initialize these
          cg::invoke_one(tb, [&]()
          {
            tol_eps = CudaMath::cuda_get_sqrt_eps<DataType>();
            need_streamline = (CudaMath::cuda_abs(sd_delta) > DataType(0)) && (sd_v_norm > tol_eps);
            need_convection = CudaMath::cuda_abs(beta) > DataType(0);
          });

          //wait for everything to be finished
          cg::wait(tb);
          tb.sync();

          //our local datatypes
          typedef Tiny::Vector<DataType, dim> VecValueType;
          typedef Tiny::Matrix<DataType, dim, dim> MatValueType;
          typedef Tiny::Vector<VecValueType, num_loc_dofs> LocalVectorType;
          typedef Tiny::Matrix<MatValueType, num_loc_dofs, num_loc_dofs> LocalMatrixType;


          // std::cout << "Starting with cell " << cell << std::endl;
          const SharedIndexSetWrapper<IndexType> local_dofs_w{local_dofs};

          // gather local vector
          SpaceHelp::grouped_set_coefficients(tb, local_coeffs, local_dofs_w, nodes, cell);
          LAFEM::template VectorGatherScatterHelper<SpaceType, DataType, IndexType>::template grouped_gather_vector_dense<cg::thread_block, dim>(tb, local_conv_dofs,
                    conv_data, IndexType(matrix_num_rows), local_dofs, num_loc_dofs, DataType(1));
          tb.sync();

          for(int mat_offset = 0; mat_offset < num_loc_dofs*num_loc_dofs; mat_offset += loc_block_size)
          {
            // VoxelAssembly::coalesced_format(tb, local_coeffs, dim*num_loc_verts);
            VoxelAssembly::coalesced_format(tb, loc_mat, loc_block_size*dim*dim);
            // VoxelAssembly::coalesced_format(tb, local_conv_dofs, num_loc_dofs*dim);

            tb.sync();

            VoxelAssembly::Kernel::template grouped_burgers_velo_material_mat_assembly_kernel<cg::thread_block, SpaceHelp, ViscoFunc_, ViscoDFunc_, stream_diff_>(tb, shared_kernel_wrapper, CudaMath::cuda_min(loc_block_size, num_loc_dofs*num_loc_dofs-mat_offset), mat_offset, loc_mat, local_conv_dofs, *((Tiny::Matrix<DataType, dim, num_loc_verts>*)local_coeffs), cub_pt, cub_wg, num_cubs, burgers_params,
                                        material_params, need_convection, tol_eps, visco_func, visco_d_func);

            tb.sync();

            //scatter
            LAFEM::template MatrixGatherScatterHelper<SpaceType, DataType, IndexType, pol_>::template grouped_scatter_matrix_csr<cg::thread_block, dim, dim>(tb, CudaMath::cuda_min(loc_block_size, num_loc_dofs*num_loc_dofs-mat_offset), mat_offset, loc_mat, matrix_data, local_dofs, local_dofs, IndexType(matrix_num_rows), IndexType(matrix_num_cols), matrix_row_ptr, matrix_col_idx, num_loc_dofs, num_loc_dofs, alpha, local_dofs_sorter);

            tb.sync();
          }
        }

      }

      template<typename Space_, typename DT_, typename IT_, typename ViscoFunc_>
      __global__ void full_burgers_vm_assembler_vector_bd(DT_* __restrict__ vector_data,
                const DT_* conv_data, const DT_* primal_data, Index vec_size,
                const Tiny::Vector<DT_, Space_::world_dim>* __restrict__ cub_pt,
                const DT_* __restrict__  cub_wg, int num_cubs, DT_ alpha,
                const IT_* __restrict__ cell_to_dof, [[maybe_unused]] Index cell_num,
                const Tiny::Vector<DT_, Space_::world_dim>* __restrict__ nodes, [[maybe_unused]] Index node_size,
                const int* __restrict__ coloring_map, Index coloring_size,
                const VoxelAssembly::AssemblyBurgersData<DT_> burgers_params,
                const ViscoFunc_ visco_func)
      {
        Index idx = threadIdx.x + blockDim.x * blockIdx.x;
        if(idx >= coloring_size)
          return;
        //define types
        typedef Space_ SpaceType;
        typedef DT_ DataType;
        typedef IT_ IndexType;

        const DataType& beta{burgers_params.beta};

        const bool need_convection = CudaMath::cuda_abs(beta) > DataType(0);

        typedef VoxelAssembly::SpaceHelper<SpaceType, DT_, IT_> SpaceHelp;

        constexpr int dim = SpaceHelp::dim;

        //define local sizes
        constexpr int num_loc_dofs = SpaceType::DofMappingType::dof_count;
        //get number of nodes per element
        constexpr int num_loc_verts = SpaceType::MeshType::template IndexSet<dim, 0>::Type::num_indices;

        //our local datatypes
        typedef Tiny::Vector<DataType, dim> VecValueType;
        typedef Tiny::Matrix<DataType, dim, dim> MatValueType;
        typedef Tiny::Vector<VecValueType, num_loc_dofs> LocalVectorType;
        typedef Tiny::Matrix<MatValueType, num_loc_dofs, num_loc_dofs> LocalMatrixType;

        //define local array
        LocalVectorType loc_vec(DataType(0));
        LocalVectorType local_conv_dofs(DataType(0));
        LocalVectorType local_prim_dofs(DataType(0));
        // typename NewSpaceHelp::ImagePointType img_point;
          Tiny::Matrix<DataType, dim, num_loc_verts> local_coeffs;

        //now do work for this cell
        Index cell = Index(coloring_map[idx]);
        const IndexSetWrapper<IndexType> local_dofs_w{cell_to_dof, IndexType(num_loc_dofs)};
        const IndexType* local_dofs = cell_to_dof + cell*num_loc_dofs;
        SpaceHelp::set_coefficients(local_coeffs, local_dofs_w, nodes, cell);
        LAFEM::template VectorGatherScatterHelper<SpaceType, DataType, IndexType>::gather_vector_dense(local_prim_dofs,
                  (const VecValueType*)primal_data, IndexType(vec_size), local_dofs, DataType(1));

        LAFEM::template VectorGatherScatterHelper<SpaceType, DataType, IndexType>::gather_vector_dense(local_conv_dofs,
                  (const VecValueType*)conv_data, IndexType(vec_size), local_dofs, DataType(1));

        VoxelAssembly::Kernel::burgers_velo_material_defect_assembly_kernel<SpaceHelp>(loc_vec, local_prim_dofs, local_conv_dofs, local_coeffs, cub_pt, cub_wg, num_cubs,
                                                         burgers_params, need_convection, visco_func);
        //scatter
        LAFEM::template VectorGatherScatterHelper<SpaceType, DataType, IndexType>::scatter_vector_dense(loc_vec,
                  (VecValueType*)vector_data, IndexType(vec_size), local_dofs, alpha);

      }

      template<typename Space_, typename DT_, typename IT_, typename ViscoFunc_>
      __global__ void full_burgers_vm_assembler_vector_bd_warp_based(DT_* __restrict__ vector_data,
                const DT_* conv_data, const DT_* primal_data, Index vec_size,
                const Tiny::Vector<DT_, Space_::world_dim>* __restrict__ cub_pt,
                const DT_* __restrict__  cub_wg, int num_cubs, DT_ alpha,
                const IT_* __restrict__ cell_to_dof, [[maybe_unused]] Index cell_num,
                const Tiny::Vector<DT_, Space_::world_dim>* __restrict__ nodes, [[maybe_unused]] Index node_size,
                const int* __restrict__ coloring_map, Index coloring_size,
                const VoxelAssembly::AssemblyBurgersData<DT_> burgers_params,
                const ViscoFunc_ visco_func,
                const int loc_block_size)
      {
        cg::thread_block tb = cg::this_thread_block();
        const Index t_idx = tb.thread_rank();
        //define types
        typedef Space_ SpaceType;
        typedef DT_ DataType;
        typedef IT_ IndexType;

        const DataType& beta{burgers_params.beta};

        typedef VoxelAssembly::SpaceHelper<SpaceType, DT_, IT_> SpaceHelp;

        constexpr int dim = SpaceHelp::dim;

        //define local sizes
        constexpr int num_loc_dofs = SpaceType::DofMappingType::dof_count;
        //get number of nodes per element
        constexpr int num_loc_verts = SpaceType::MeshType::template IndexSet<dim, 0>::Type::num_indices;

        //our local datatypes
        typedef Tiny::Vector<DataType, dim> VecValueType;
        typedef Tiny::Matrix<DataType, dim, dim> MatValueType;
        typedef Tiny::Vector<VecValueType, num_loc_dofs> LocalVectorType;
        typedef Tiny::Matrix<MatValueType, num_loc_dofs, num_loc_dofs> LocalMatrixType;

        extern __shared__ unsigned char shared_array[];

        typedef BurgersDefectSharedDataGlobalWrapper<SpaceHelp> BSDGWrapper;
        typedef BurgersMatSharedDataKernelWrapper<SpaceHelp, false> BSDKWrapper;

        BSDGWrapper* shared_meta_wrapper =
            reinterpret_cast<BSDGWrapper*>(&shared_array[0]);
        // first define local_coeffs, since this will be cast to a tiny matrix, which is required to be aligned...
        DataType* local_coeffs = &(shared_meta_wrapper->local_coeffs[0]);
        DataType* local_conv_dofs = &(shared_meta_wrapper->local_conv_dofs[0]);
        DataType* local_prim_dofs = &(shared_meta_wrapper->local_prim_dofs[0]);
        IndexType* local_dofs = &(shared_meta_wrapper->local_dofs[0]);

        // these can be defined in shared memory
        bool& need_convection = shared_meta_wrapper->need_convection;
        cg::invoke_one(tb, [&](){
            need_convection = CudaMath::cuda_abs(beta) > DataType(0);
        });
        // shared array for kernel
        BSDKWrapper* shared_kernel_wrapper = (BSDKWrapper*)(shared_meta_wrapper+1);
        // the shared array for the local vector to be written out
        DataType* loc_vec = VoxelAssembly::get_aligned_address<DataType>((void*)(&shared_array[0]+sizeof(BSDGWrapper)+sizeof(BSDKWrapper)));
        // now do work for this cell
        //stride based for loop
        for(int b_idx = blockIdx.x; b_idx < int(coloring_size); b_idx += gridDim.x)
        {
          //get our actual cell index
          const Index cell = Index(coloring_map[b_idx]);

          cg::memcpy_async(tb, local_dofs, cell_to_dof+num_loc_dofs*cell, num_loc_dofs*sizeof(IndexType));
          VoxelAssembly::coalesced_format(tb, local_coeffs, dim*num_loc_verts);
          VoxelAssembly::coalesced_format(tb, loc_vec, loc_block_size*dim);
          VoxelAssembly::coalesced_format(tb, local_conv_dofs, num_loc_dofs*dim);
          VoxelAssembly::coalesced_format(tb, local_prim_dofs, num_loc_dofs*dim);
          cg::wait(tb);



          const SharedIndexSetWrapper<IndexType> local_dofs_w{local_dofs};
          SpaceHelp::grouped_set_coefficients(tb, local_coeffs, local_dofs_w, nodes, cell);
          LAFEM::template VectorGatherScatterHelper<SpaceType, DataType, IndexType>::template grouped_gather_vector_dense<cg::thread_block, dim>(tb, local_prim_dofs,
                    primal_data, IndexType(vec_size), local_dofs, num_loc_dofs, DataType(1));

          LAFEM::template VectorGatherScatterHelper<SpaceType, DataType, IndexType>::template grouped_gather_vector_dense<cg::thread_block, dim>(tb, local_conv_dofs,
                    conv_data, IndexType(vec_size), local_dofs, num_loc_dofs, DataType(1));
          tb.sync();
          for(int vec_offset = 0; vec_offset < num_loc_dofs; vec_offset += loc_block_size)
          {
            VoxelAssembly::coalesced_format(tb, loc_vec, loc_block_size*dim);

            tb.sync();

            VoxelAssembly::Kernel::template grouped_burgers_velo_material_defect_assembly_kernel<cg::thread_block, SpaceHelp, ViscoFunc_>(tb, shared_kernel_wrapper, CudaMath::cuda_min(loc_block_size, num_loc_dofs-vec_offset), vec_offset, loc_vec,
                                                          local_prim_dofs, local_conv_dofs, *((Tiny::Matrix<DataType, dim, num_loc_verts>*)local_coeffs), cub_pt, cub_wg, num_cubs,
                                                          burgers_params, need_convection, visco_func);

            tb.sync();
            //scatter
            LAFEM::template VectorGatherScatterHelper<SpaceType, DataType, IndexType>::template grouped_scatter_vector_dense<cg::thread_block, dim>(tb, CudaMath::cuda_min(loc_block_size, num_loc_dofs-vec_offset), vec_offset, loc_vec,
                      vector_data, IndexType(vec_size), local_dofs, num_loc_dofs, alpha);
            tb.sync();
          }
        }

      }

      /**************************************************************************************************************/
      /*                                       CUDA Host OMP Kernels                                                */
      /**************************************************************************************************************/

      template<typename Space_, typename DT_, typename IT_, typename ViscoFunc_, typename ViscoDFunc_,
              FEAT::Intern::MatrixGatherScatterPolicy pol_ = FEAT::Intern::MatrixGatherScatterPolicy::useLocalOps>
      void full_burgers_vm_assembler_matrix1_bcsr_host(DT_* matrix_data, const DT_* conv_data,
                const IT_*  matrix_row_ptr, const IT_* matrix_col_idx, Index matrix_num_rows, Index matrix_num_cols,
                const Tiny::Vector<DT_, Space_::world_dim>* cub_pt,
                const DT_*  cub_wg, int num_cubs, DT_ alpha,
                const IT_* cell_to_dof, [[maybe_unused]] Index cell_num,
                const Tiny::Vector<DT_, Space_::world_dim>* nodes, [[maybe_unused]] Index node_size,
                const int* coloring_map, Index coloring_size, const IT_* cell_to_dof_sorter,
                const VoxelAssembly::AssemblyBurgersData<DT_>& burgers_params,
                const VoxelAssembly::AssemblyMaterialData<DT_>& material_params,
                const ViscoFunc_& visco_func, const ViscoDFunc_& visco_d_func)
      {
        //define types
        typedef Space_ SpaceType;
        typedef DT_ DataType;
        typedef IT_ IndexType;

        const DataType& beta{burgers_params.beta};
        const DataType& sd_delta{burgers_params.sd_delta};
        const DataType& sd_v_norm{burgers_params.sd_v_norm};

        const DataType tol_eps = CudaMath::cuda_get_sqrt_eps<DataType>();
        const bool need_streamline = (CudaMath::cuda_abs(sd_delta) > DataType(0)) && (sd_v_norm > tol_eps);
        const bool need_convection = CudaMath::cuda_abs(beta) > DataType(0);

        typedef VoxelAssembly::SpaceHelper<SpaceType, DT_, IT_> SpaceHelp;

        constexpr int dim = SpaceHelp::dim;

        //define local sizes
        constexpr int num_loc_dofs = SpaceType::DofMappingType::dof_count;
        //get number of nodes per element
        constexpr int num_loc_verts = SpaceType::MeshType::template IndexSet<dim, 0>::Type::num_indices;

        //our local datatypes
        typedef Tiny::Vector<DataType, dim> VecValueType;
        typedef Tiny::Matrix<DataType, dim, dim> MatValueType;
        typedef Tiny::Vector<VecValueType, num_loc_dofs> LocalVectorType;
        typedef Tiny::Matrix<MatValueType, num_loc_dofs, num_loc_dofs> LocalMatrixType;

        #pragma omp parallel for
        for(Index idx = 0; idx < coloring_size; ++idx)
        {
          // define local coefficients
          Tiny::Matrix<DataType, dim, num_loc_verts> local_coeffs;
          // and local matrix
          LocalMatrixType loc_mat;
          LocalVectorType local_conv_dofs(DataType(0));
          // now do work for this cell
          int cell = coloring_map[idx];
          // std::cout << "Starting with cell " << cell << std::endl;
          const IndexSetWrapper<IndexType> local_dofs_w{cell_to_dof, IndexType(num_loc_dofs)};
          const IndexType* local_dofs = cell_to_dof + IndexType(cell)*num_loc_dofs;
          const IndexType* local_dof_sorter = cell_to_dof_sorter + IndexType(cell)*num_loc_dofs;

          SpaceHelp::set_coefficients(local_coeffs, local_dofs_w, nodes, cell);
          LAFEM::template VectorGatherScatterHelper<SpaceType, DataType, IndexType>::gather_vector_dense(local_conv_dofs,
                    (const VecValueType*)conv_data, IndexType(matrix_num_rows), local_dofs,DataType(1));

          if(need_streamline)
          {
            VoxelAssembly::Kernel::burgers_velo_material_mat_assembly_kernel<SpaceHelp, LocalMatrixType, LocalVectorType, dim, num_loc_verts, ViscoFunc_, ViscoDFunc_, true>(
                            loc_mat, local_conv_dofs, local_coeffs, cub_pt, cub_wg, num_cubs, burgers_params, material_params, need_convection, tol_eps, visco_func, visco_d_func);
          }
          else
          {
            VoxelAssembly::Kernel::burgers_velo_material_mat_assembly_kernel<SpaceHelp, LocalMatrixType, LocalVectorType, dim, num_loc_verts, ViscoFunc_, ViscoDFunc_, false>(
                            loc_mat, local_conv_dofs, local_coeffs, cub_pt, cub_wg, num_cubs, burgers_params, material_params, need_convection, tol_eps, visco_func, visco_d_func);
          }

          //scatter
          LAFEM::template MatrixGatherScatterHelper<SpaceType, DataType, IndexType, pol_>::scatter_matrix_csr(loc_mat, (MatValueType*)matrix_data, local_dofs, local_dofs, IndexType(matrix_num_rows), IndexType(matrix_num_cols), matrix_row_ptr, matrix_col_idx, alpha, local_dof_sorter);

        }
      }

      template<typename Space_, typename DT_, typename IT_, typename ViscoFunc_>
      void full_burgers_vm_assembler_vector_bd_host(DT_* vector_data,
                const DT_* conv_data, const DT_* primal_data, Index vec_size,
                const Tiny::Vector<DT_, Space_::world_dim>* cub_pt,
                const DT_*  cub_wg, int num_cubs, DT_ alpha,
                const IT_* cell_to_dof, [[maybe_unused]] Index cell_num,
                const Tiny::Vector<DT_, Space_::world_dim>* nodes, [[maybe_unused]] Index node_size,
                const int* coloring_map, Index coloring_size,
                const VoxelAssembly::AssemblyBurgersData<DT_>& burgers_params,
                const ViscoFunc_& visco_func)
      {
        //define types
        typedef Space_ SpaceType;
        typedef DT_ DataType;
        typedef IT_ IndexType;

        const DataType& beta{burgers_params.beta};

        const bool need_convection = CudaMath::cuda_abs(beta) > DataType(0);

        typedef VoxelAssembly::SpaceHelper<SpaceType, DT_, IT_> SpaceHelp;

        constexpr int dim = SpaceHelp::dim;

        //define local sizes
        constexpr int num_loc_dofs = SpaceType::DofMappingType::dof_count;
        //get number of nodes per element
        constexpr int num_loc_verts = SpaceType::MeshType::template IndexSet<dim, 0>::Type::num_indices;

        //our local datatypes
        typedef Tiny::Vector<DataType, dim> VecValueType;
        typedef Tiny::Matrix<DataType, dim, dim> MatValueType;
        typedef Tiny::Vector<VecValueType, num_loc_dofs> LocalVectorType;
        typedef Tiny::Matrix<MatValueType, num_loc_dofs, num_loc_dofs> LocalMatrixType;

        #pragma omp parallel for
        for(Index idx = 0; idx < coloring_size; ++idx)
        {
          //define local array
          LocalVectorType loc_vec(DataType(0));
          LocalVectorType local_conv_dofs(DataType(0));
          LocalVectorType local_prim_dofs(DataType(0));
          // typename NewSpaceHelp::ImagePointType img_point;
          Tiny::Matrix<DataType, dim, num_loc_verts> local_coeffs;

          //now do work for this cell
          int cell = coloring_map[idx];
          const IndexSetWrapper<IndexType> local_dofs_w{cell_to_dof, IndexType(num_loc_dofs)};
          const IndexType* local_dofs = cell_to_dof + IndexType(cell)*num_loc_dofs;
          SpaceHelp::set_coefficients(local_coeffs, local_dofs_w, nodes, cell);
          LAFEM::template VectorGatherScatterHelper<SpaceType, DataType, IndexType>::gather_vector_dense(local_prim_dofs,
                    (const VecValueType*)primal_data, IndexType(vec_size), local_dofs, DataType(1));

          //if we need to, gather local convection vector
          LAFEM::template VectorGatherScatterHelper<SpaceType, DataType, IndexType>::gather_vector_dense(local_conv_dofs,
                    (const VecValueType*)conv_data, IndexType(vec_size), local_dofs, DataType(1));

          VoxelAssembly::Kernel::burgers_velo_material_defect_assembly_kernel<SpaceHelp>(loc_vec, local_prim_dofs, local_conv_dofs, local_coeffs, cub_pt, cub_wg, num_cubs,
                                                         burgers_params, need_convection, visco_func);
          //scatter
          LAFEM::template VectorGatherScatterHelper<SpaceType, DataType, IndexType>::scatter_vector_dense(loc_vec,
                    (VecValueType*)vector_data, IndexType(vec_size), local_dofs, alpha);

        }
      }

    }

    namespace Arch
    {
      // template<typename Space_, typename DT_, typename IT_>
      // void assemble_burgers_velo_material_csr_old([[maybe_unused]] const Space_& space,
      //         const CSRMatrixData<DT_, IT_>& matrix_data,
      //         const DT_* conv_data,
      //         const AssemblyCubatureData<DT_>& cubature,
      //         const AssemblyMappingData<DT_, IT_>& dof_mapping,
      //         const std::vector<int*>& coloring_maps,
      //         const std::vector<Index>& coloring_map_sizes,
      //         DT_ alpha, const AssemblyBurgersData<DT_>& burgers_params,
      //         const AssemblyMaterialData<DT_>& material_params,
      //         MaterialType material_type)
      // {
      //   const Index blocksize = Util::cuda_blocksize_blocked_assembly;
      //   // const Index blocksize = 64;
      //   typedef DT_ DataType;
      //   const DataType tol_eps = CudaMath::cuda_get_sqrt_eps<DataType>();
      //   const bool need_streamline = (CudaMath::cuda_abs(burgers_params.sd_delta) > DataType(0)) && (burgers_params.sd_v_norm > tol_eps);

      //   for(Index col = 0; col < coloring_maps.size(); ++col)
      //   {
      //     dim3 grid;
      //     dim3 block;
      //     block.x = (unsigned int)blocksize;
      //     grid.x = (unsigned int)ceil(double(coloring_map_sizes[col])/double(block.x));

      //     if(need_streamline)
      //     {
      //       switch(material_type)
      //       {
      //         case FEAT::VoxelAssembly::MaterialType::carreau:
      //         {
      //           //kernel call, since this uses the standard stream, sync before next call is enforced:
      //           VoxelAssembly::Kernel::template full_burgers_vm_assembler_matrix1_bcsr<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreau>,
      //                           Intern::ViscoDFunctor<DT_, MaterialType::carreau>, true, FEAT::Intern::MatrixGatherScatterPolicy::useLocalSortHelper><<< grid, block >>>(
      //               matrix_data.data, conv_data, matrix_data.row_ptr, matrix_data.col_idx, matrix_data.num_rows, matrix_data.num_cols,
      //               (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
      //               cubature.cub_wg, cubature.num_cubs, alpha,
      //               dof_mapping.cell_to_dof, dof_mapping.cell_num,
      //               (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
      //               (const int*) coloring_maps[col], coloring_map_sizes[col], dof_mapping.cell_to_dof_sorter,
      //               burgers_params, material_params, Intern::ViscoFunctor<DT_, MaterialType::carreau>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T},
      //               Intern::ViscoDFunctor<DT_, MaterialType::carreau>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T}
      //           );
      //           break;
      //         }
      //         case FEAT::VoxelAssembly::MaterialType::carreauYasuda:
      //         {
      //           //kernel call, since this uses the standard stream, sync before next call is enforced:
      //           VoxelAssembly::Kernel::template full_burgers_vm_assembler_matrix1_bcsr<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreauYasuda>,
      //                           Intern::ViscoDFunctor<DT_, MaterialType::carreauYasuda>, true, FEAT::Intern::MatrixGatherScatterPolicy::useLocalSortHelper><<< grid, block >>>(
      //               matrix_data.data, conv_data, matrix_data.row_ptr, matrix_data.col_idx, matrix_data.num_rows, matrix_data.num_cols,
      //               (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
      //               cubature.cub_wg, cubature.num_cubs, alpha,
      //               dof_mapping.cell_to_dof, dof_mapping.cell_num,
      //               (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
      //               (const int*) coloring_maps[col], coloring_map_sizes[col], dof_mapping.cell_to_dof_sorter,
      //               burgers_params, material_params,
      //               Intern::ViscoFunctor<DT_, MaterialType::carreauYasuda>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T,material_params.yasuda_a, material_params.mu_inf},
      //               Intern::ViscoDFunctor<DT_, MaterialType::carreauYasuda>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T, material_params.yasuda_a, material_params.mu_inf}
      //           );
      //           break;
      //         }
      //       }
      //     }
      //     else
      //     {
      //       switch(material_type)
      //       {
      //         case FEAT::VoxelAssembly::MaterialType::carreau:
      //         {
      //           //kernel call, since this uses the standard stream, sync before next call is enforced:
      //           VoxelAssembly::Kernel::template full_burgers_vm_assembler_matrix1_bcsr<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreau>,
      //                           Intern::ViscoDFunctor<DT_, MaterialType::carreau>, false, FEAT::Intern::MatrixGatherScatterPolicy::useLocalSortHelper><<< grid, block >>>(
      //               matrix_data.data, conv_data, matrix_data.row_ptr, matrix_data.col_idx, matrix_data.num_rows, matrix_data.num_cols,
      //               (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
      //               cubature.cub_wg, cubature.num_cubs, alpha,
      //               dof_mapping.cell_to_dof, dof_mapping.cell_num,
      //               (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
      //               (const int*) coloring_maps[col], coloring_map_sizes[col], dof_mapping.cell_to_dof_sorter,
      //               burgers_params, material_params, Intern::ViscoFunctor<DT_, MaterialType::carreau>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T},
      //               Intern::ViscoDFunctor<DT_, MaterialType::carreau>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T}
      //           );
      //           break;
      //         }
      //         case FEAT::VoxelAssembly::MaterialType::carreauYasuda:
      //         {
      //           //kernel call, since this uses the standard stream, sync before next call is enforced:
      //           VoxelAssembly::Kernel::template full_burgers_vm_assembler_matrix1_bcsr<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreauYasuda>,
      //                           Intern::ViscoDFunctor<DT_, MaterialType::carreauYasuda>, false, FEAT::Intern::MatrixGatherScatterPolicy::useLocalSortHelper><<< grid, block >>>(
      //               matrix_data.data, conv_data, matrix_data.row_ptr, matrix_data.col_idx, matrix_data.num_rows, matrix_data.num_cols,
      //               (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
      //               cubature.cub_wg, cubature.num_cubs, alpha,
      //               dof_mapping.cell_to_dof, dof_mapping.cell_num,
      //               (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
      //               (const int*) coloring_maps[col], coloring_map_sizes[col], dof_mapping.cell_to_dof_sorter,
      //               burgers_params, material_params,
      //               Intern::ViscoFunctor<DT_, MaterialType::carreauYasuda>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T,material_params.yasuda_a, material_params.mu_inf},
      //               Intern::ViscoDFunctor<DT_, MaterialType::carreauYasuda>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T, material_params.yasuda_a, material_params.mu_inf}
      //           );
      //           break;
      //         }
      //       }
      //     }
      //   }
      //   //check for cuda error in our kernel
      //   Util::cuda_check_last_error();
      // }

      template<typename Space_, typename DT_, typename IT_>
      inline int get_num_block_warp(MaterialType mat_type, bool need_streamdiff, int blocksize, int shared_mem)
      {
        switch(mat_type)
        {
          case FEAT::VoxelAssembly::MaterialType::carreau:
          {
            if(need_streamdiff)
              return Util::cuda_get_occupancy(VoxelAssembly::Kernel::template full_burgers_vm_assembler_matrix1_bcsr_warp_based<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreau>, Intern::ViscoDFunctor<DT_, MaterialType::carreau>, true>, blocksize, shared_mem);
            else
              return Util::cuda_get_occupancy(VoxelAssembly::Kernel::template full_burgers_vm_assembler_matrix1_bcsr_warp_based<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreau>, Intern::ViscoDFunctor<DT_, MaterialType::carreau>, false>, blocksize, shared_mem);
            break;
          }
          case FEAT::VoxelAssembly::MaterialType::carreauYasuda:
          {
            if(need_streamdiff)
              return Util::cuda_get_occupancy(VoxelAssembly::Kernel::template full_burgers_vm_assembler_matrix1_bcsr_warp_based<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreauYasuda>, Intern::ViscoDFunctor<DT_, MaterialType::carreauYasuda>, true>, blocksize, shared_mem);
            else
              return Util::cuda_get_occupancy(VoxelAssembly::Kernel::template full_burgers_vm_assembler_matrix1_bcsr_warp_based<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreauYasuda>, Intern::ViscoDFunctor<DT_, MaterialType::carreauYasuda>, false>, blocksize, shared_mem);
            break;
          }
          default:
          {
            XABORTM("Not implemented!");
          }
        }
      }

      template<typename Space_, typename DT_, typename IT_>
      inline int get_num_block_warp_defect(MaterialType mat_type, int blocksize, int shared_mem)
      {
        switch(mat_type)
        {
          case FEAT::VoxelAssembly::MaterialType::carreau:
          {
            return Util::cuda_get_occupancy(VoxelAssembly::Kernel::template full_burgers_vm_assembler_vector_bd_warp_based<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreau>>, blocksize, shared_mem);
            break;
          }
          case FEAT::VoxelAssembly::MaterialType::carreauYasuda:
          {
            return Util::cuda_get_occupancy(VoxelAssembly::Kernel::template full_burgers_vm_assembler_vector_bd_warp_based<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreauYasuda>>, blocksize, shared_mem);
            break;
          }
          default:
          {
            XABORTM("Not implemented!");
          }
        }
      }

      template<typename Space_, typename DT_, typename IT_>
      inline bool set_shared_mem_attribute(MaterialType mat_type, bool need_streamdiff, int shared_mem)
      {
        switch(mat_type)
        {
          case FEAT::VoxelAssembly::MaterialType::carreau:
          {
            if(need_streamdiff)
              return cudaSuccess == cudaFuncSetAttribute(VoxelAssembly::Kernel::template full_burgers_vm_assembler_matrix1_bcsr_warp_based<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreau>, Intern::ViscoDFunctor<DT_, MaterialType::carreau>, true>, cudaFuncAttributeMaxDynamicSharedMemorySize, shared_mem);
            else
              return cudaSuccess == cudaFuncSetAttribute(VoxelAssembly::Kernel::template full_burgers_vm_assembler_matrix1_bcsr_warp_based<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreau>, Intern::ViscoDFunctor<DT_, MaterialType::carreau>, false>, cudaFuncAttributeMaxDynamicSharedMemorySize, shared_mem);
            break;
          }
          case FEAT::VoxelAssembly::MaterialType::carreauYasuda:
          {
            if(need_streamdiff)
              return cudaSuccess == cudaFuncSetAttribute(VoxelAssembly::Kernel::template full_burgers_vm_assembler_matrix1_bcsr_warp_based<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreauYasuda>, Intern::ViscoDFunctor<DT_, MaterialType::carreauYasuda>, true>, cudaFuncAttributeMaxDynamicSharedMemorySize, shared_mem);
            else
              return cudaSuccess == cudaFuncSetAttribute(VoxelAssembly::Kernel::template full_burgers_vm_assembler_matrix1_bcsr_warp_based<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreauYasuda>, Intern::ViscoDFunctor<DT_, MaterialType::carreauYasuda>, false>, cudaFuncAttributeMaxDynamicSharedMemorySize, shared_mem);
            break;
          }
          default:
          {
            XABORTM("Not implemented!");
          }
        }
      }

      template<typename Space_, typename DT_, typename IT_>
      inline bool set_shared_mem_attribute_defect(MaterialType mat_type, int shared_mem)
      {
        switch(mat_type)
        {
          case FEAT::VoxelAssembly::MaterialType::carreau:
          {
            return cudaSuccess == cudaFuncSetAttribute(VoxelAssembly::Kernel::template full_burgers_vm_assembler_vector_bd_warp_based<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreau>>,cudaFuncAttributeMaxDynamicSharedMemorySize, shared_mem);
            break;
          }
          case FEAT::VoxelAssembly::MaterialType::carreauYasuda:
          {
            return cudaSuccess == cudaFuncSetAttribute(VoxelAssembly::Kernel::template full_burgers_vm_assembler_vector_bd_warp_based<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreauYasuda>>,cudaFuncAttributeMaxDynamicSharedMemorySize, shared_mem);
            break;
          }
          default:
          {
            XABORTM("Not implemented!");
          }
        }
      }

      template<typename Space_, typename DT_, typename IT_>
      inline int get_min_shared_size(bool streamline)
      {
        return sizeof(VoxelAssembly::BurgersSharedDataGlobalWrapper<VoxelAssembly::SpaceHelper<Space_, DT_, IT_>>) +
          (streamline ? sizeof(VoxelAssembly::BurgersMatSharedDataKernelWrapper<VoxelAssembly::SpaceHelper<Space_, DT_, IT_>, true>) :
                        sizeof(VoxelAssembly::BurgersMatSharedDataKernelWrapper<VoxelAssembly::SpaceHelper<Space_, DT_, IT_>, false>));
      }

      template<typename Space_, typename DT_, typename IT_>
      void assemble_burgers_velo_material_csr([[maybe_unused]] const Space_& space,
              const CSRMatrixData<DT_, IT_>& matrix_data,
              const DT_* conv_data,
              const AssemblyCubatureData<DT_>& cubature,
              const AssemblyMappingData<DT_, IT_>& dof_mapping,
              const std::vector<int*>& coloring_maps,
              const std::vector<Index>& coloring_map_sizes,
              DT_ alpha, const AssemblyBurgersData<DT_>& burgers_params,
              const AssemblyMaterialData<DT_>& material_params,
              MaterialType material_type, int shared_mem, int blocksize, int gridsize, bool print_occupancy)
      {
        // get the size of all necessary shared memory
        typedef DT_ DataType;
        const DataType tol_eps = CudaMath::cuda_get_sqrt_eps<DataType>();
        constexpr int dim = Space_::ShapeType::dimension;
        constexpr int num_loc_dofs = Space_::DofMappingType::dof_count;
        const bool need_streamline = (CudaMath::cuda_abs(burgers_params.sd_delta) > DataType(0)) && (burgers_params.sd_v_norm > tol_eps);
        // overestimate necessary memory a bit... else this gets to complicated...?!
        const int shared_size_nec = get_min_shared_size<Space_, DT_, IT_>(need_streamline);
        const int loc_mat_shared_mem = shared_mem - shared_size_nec;
        XASSERTM(loc_mat_shared_mem > 0, String("Not enough assigned shared memory\n Minimum required: ") + stringify(shared_size_nec)
            + String("\nProvided: ") + stringify(shared_mem));

        const int loc_block_size = CudaMath::cuda_min(loc_mat_shared_mem / int(sizeof(DT_)*dim*dim), num_loc_dofs*num_loc_dofs);
        XASSERTM(loc_block_size > 0, "Not enough memory assigned to assemble a single matrix entry!");
        XASSERTM(loc_block_size >= num_loc_dofs, "Not enough memory assigned to assemble a single matrix row!");

        if(print_occupancy)
        {
          int numBlocksWarp = get_num_block_warp<Space_, DT_, IT_>(material_type, need_streamline, blocksize, shared_mem);
          const int max_blocks_per_sm = int(Util::cuda_get_max_blocks_per_sm());
          printf("Numblocks/Occupancy per SM for device number %i: %i, %f\n", Util::cuda_device_number, numBlocksWarp, double(numBlocksWarp*(blocksize/32))/double(max_blocks_per_sm));
        }

        if(shared_mem > 48000)
        {
          if(!set_shared_mem_attribute<Space_, DT_, IT_>(material_type, need_streamline, shared_mem))
          {
            XABORTM("cudaFuncSetAttribute failed.");
          }
        }
        // const Index blocksize = 64;

        for(Index col = 0; col < coloring_maps.size(); ++col)
        {
          dim3 grid;
          dim3 block;
          block.x = (unsigned int)blocksize;
          grid.x = (unsigned int)gridsize;

          if(need_streamline)
          {
            switch(material_type)
            {
              case FEAT::VoxelAssembly::MaterialType::carreau:
              {
                //kernel call, since this uses the standard stream, sync before next call is enforced:
                VoxelAssembly::Kernel::template full_burgers_vm_assembler_matrix1_bcsr_warp_based<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreau>,
                                Intern::ViscoDFunctor<DT_, MaterialType::carreau>, true, FEAT::Intern::MatrixGatherScatterPolicy::useLocalOps><<< grid, block, shared_mem >>>(
                    matrix_data.data, conv_data, matrix_data.row_ptr, matrix_data.col_idx, matrix_data.num_rows, matrix_data.num_cols,
                    (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
                    cubature.cub_wg, cubature.num_cubs, alpha,
                    dof_mapping.cell_to_dof, dof_mapping.cell_num,
                    (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
                    (const int*) coloring_maps[col], coloring_map_sizes[col], nullptr,
                    burgers_params, material_params, Intern::ViscoFunctor<DT_, MaterialType::carreau>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T},
                    Intern::ViscoDFunctor<DT_, MaterialType::carreau>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T},
                    loc_block_size
                );
                break;
              }
              case FEAT::VoxelAssembly::MaterialType::carreauYasuda:
              {
                //kernel call, since this uses the standard stream, sync before next call is enforced:
                VoxelAssembly::Kernel::template full_burgers_vm_assembler_matrix1_bcsr_warp_based<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreauYasuda>,
                                Intern::ViscoDFunctor<DT_, MaterialType::carreauYasuda>, true, FEAT::Intern::MatrixGatherScatterPolicy::useLocalOps><<< grid, block, shared_mem >>>(
                    matrix_data.data, conv_data, matrix_data.row_ptr, matrix_data.col_idx, matrix_data.num_rows, matrix_data.num_cols,
                    (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
                    cubature.cub_wg, cubature.num_cubs, alpha,
                    dof_mapping.cell_to_dof, dof_mapping.cell_num,
                    (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
                    (const int*) coloring_maps[col], coloring_map_sizes[col], nullptr,
                    burgers_params, material_params,
                    Intern::ViscoFunctor<DT_, MaterialType::carreauYasuda>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T,material_params.yasuda_a, material_params.mu_inf},
                    Intern::ViscoDFunctor<DT_, MaterialType::carreauYasuda>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T, material_params.yasuda_a, material_params.mu_inf},
                    loc_block_size
                );
                break;
              }
            }
          }
          else
          {
            switch(material_type)
            {
              case FEAT::VoxelAssembly::MaterialType::carreau:
              {
                //kernel call, since this uses the standard stream, sync before next call is enforced:
                VoxelAssembly::Kernel::template full_burgers_vm_assembler_matrix1_bcsr_warp_based<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreau>,
                                Intern::ViscoDFunctor<DT_, MaterialType::carreau>, false, FEAT::Intern::MatrixGatherScatterPolicy::useLocalOps><<< grid, block, shared_mem >>>(
                    matrix_data.data, conv_data, matrix_data.row_ptr, matrix_data.col_idx, matrix_data.num_rows, matrix_data.num_cols,
                    (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
                    cubature.cub_wg, cubature.num_cubs, alpha,
                    dof_mapping.cell_to_dof, dof_mapping.cell_num,
                    (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
                    (const int*) coloring_maps[col], coloring_map_sizes[col], nullptr,
                    burgers_params, material_params, Intern::ViscoFunctor<DT_, MaterialType::carreau>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T},
                    Intern::ViscoDFunctor<DT_, MaterialType::carreau>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T},
                    loc_block_size
                );
                break;
              }
              case FEAT::VoxelAssembly::MaterialType::carreauYasuda:
              {
                //kernel call, since this uses the standard stream, sync before next call is enforced:
                VoxelAssembly::Kernel::template full_burgers_vm_assembler_matrix1_bcsr_warp_based<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreauYasuda>,
                                Intern::ViscoDFunctor<DT_, MaterialType::carreauYasuda>, false, FEAT::Intern::MatrixGatherScatterPolicy::useLocalOps><<< grid, block, shared_mem >>>(
                    matrix_data.data, conv_data, matrix_data.row_ptr, matrix_data.col_idx, matrix_data.num_rows, matrix_data.num_cols,
                    (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
                    cubature.cub_wg, cubature.num_cubs, alpha,
                    dof_mapping.cell_to_dof, dof_mapping.cell_num,
                    (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
                    (const int*) coloring_maps[col], coloring_map_sizes[col], nullptr,
                    burgers_params, material_params,
                    Intern::ViscoFunctor<DT_, MaterialType::carreauYasuda>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T,material_params.yasuda_a, material_params.mu_inf},
                    Intern::ViscoDFunctor<DT_, MaterialType::carreauYasuda>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T, material_params.yasuda_a, material_params.mu_inf},
                    loc_block_size
                );
                break;
              }
            }
          }
        }
        //check for cuda error in our kernel
        Util::cuda_check_last_error();
      }


      template<typename Space_, typename DT_, typename IT_>
      void assemble_burgers_velo_material_defect_old([[maybe_unused]] const Space_& space,
              DT_* vector_data,
              const DT_* conv_data,
              const DT_* primal_data,
              const AssemblyCubatureData<DT_>& cubature,
              const AssemblyMappingData<DT_, IT_>& dof_mapping,
              const std::vector<int*>& coloring_maps,
              const std::vector<Index>& coloring_map_sizes,
              DT_ alpha, const AssemblyBurgersData<DT_>& burgers_params,
              const AssemblyMaterialData<DT_>& material_params,
              MaterialType material_type)
      {
        const Index blocksize = Util::cuda_blocksize_scalar_assembly;
        // const Index blocksize = 64;
        #ifdef DEBUG
        std::size_t pValue;
        cudaDeviceGetLimit(&pValue, cudaLimit::cudaLimitStackSize);
        std::cout << "Current maximum thread cache: " << pValue << std::endl;
        #endif
        for(Index col = 0; col < coloring_maps.size(); ++col)
        {
          dim3 grid;
          dim3 block;
          block.x = (unsigned int)blocksize;
          grid.x = (unsigned int)ceil(double(coloring_map_sizes[col])/double(block.x));
          switch(material_type)
          {
            case FEAT::VoxelAssembly::MaterialType::carreau:
            {
              VoxelAssembly::Kernel::template full_burgers_vm_assembler_vector_bd<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreau>><<< grid, block >>>(vector_data, conv_data, primal_data,
                  space.get_num_dofs(), (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
                  cubature.cub_wg, cubature.num_cubs, alpha,
                  dof_mapping.cell_to_dof, dof_mapping.cell_num,
                  (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
                  (const int*) coloring_maps[col], coloring_map_sizes[col], burgers_params,
                  Intern::ViscoFunctor<DT_, MaterialType::carreau>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T}
              );
              break;
            }
            case FEAT::VoxelAssembly::MaterialType::carreauYasuda:
            {
              VoxelAssembly::Kernel::template full_burgers_vm_assembler_vector_bd<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreauYasuda>><<< grid, block >>>(vector_data, conv_data, primal_data,
                  space.get_num_dofs(), (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
                  cubature.cub_wg, cubature.num_cubs, alpha,
                  dof_mapping.cell_to_dof, dof_mapping.cell_num,
                  (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
                  (const int*) coloring_maps[col], coloring_map_sizes[col], burgers_params,
                  Intern::ViscoFunctor<DT_, MaterialType::carreauYasuda>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T,material_params.yasuda_a, material_params.mu_inf}
              );
              break;
            }
          }
        }

        //check for cuda error (also synchronizes)
        Util::cuda_check_last_error();
      }

      template<typename Space_, typename DT_, typename IT_>
      void assemble_burgers_velo_material_defect([[maybe_unused]] const Space_& space,
              DT_* vector_data,
              const DT_* conv_data,
              const DT_* primal_data,
              const AssemblyCubatureData<DT_>& cubature,
              const AssemblyMappingData<DT_, IT_>& dof_mapping,
              const std::vector<int*>& coloring_maps,
              const std::vector<Index>& coloring_map_sizes,
              DT_ alpha, const AssemblyBurgersData<DT_>& burgers_params,
              const AssemblyMaterialData<DT_>& material_params,
              MaterialType material_type,
              int shared_mem, int blocksize, int gridsize, bool print_occupancy)
      {
        // get the size of all necessary shared memory
        constexpr int dim = Space_::ShapeType::dimension;
        constexpr int num_loc_dofs = Space_::DofMappingType::dof_count;
        const int shared_size_nec = get_min_shared_size<Space_, DT_, IT_>(false);
        const int loc_vec_shared_mem = shared_mem - shared_size_nec;
        XASSERTM(loc_vec_shared_mem > 0, String("Not enough assigned shared memory\n Minimum required: ") + stringify(shared_size_nec)
            + String("\nProvided: ") + stringify(shared_mem));

        const int loc_block_size = CudaMath::cuda_min(loc_vec_shared_mem / int(sizeof(DT_)*dim), num_loc_dofs);
        XASSERTM(loc_block_size > 0, "Not enough memory assigned to assemble a single defect entry!");

        if(print_occupancy)
        {
          int numBlocksWarp = get_num_block_warp_defect<Space_, DT_, IT_>(material_type, blocksize, shared_mem);
          const int max_blocks_per_sm = int(Util::cuda_get_max_blocks_per_sm());
          printf("Numblocks/Occupancy per SM for device number %i: %i, %f\n", Util::cuda_device_number, numBlocksWarp, double(numBlocksWarp*(blocksize/32))/double(max_blocks_per_sm));
        }

        if(shared_mem > 48000)
        {
          if(!set_shared_mem_attribute_defect<Space_, DT_, IT_>(material_type, shared_mem))
          {
            XABORTM("cudaFuncSetAttribute failed.");
          }
        }

        for(Index col = 0; col < coloring_maps.size(); ++col)
        {
          dim3 grid;
          dim3 block;
          block.x = (unsigned int)blocksize;
          grid.x = (unsigned int)gridsize;
          switch(material_type)
          {
            case FEAT::VoxelAssembly::MaterialType::carreau:
            {
              VoxelAssembly::Kernel::template full_burgers_vm_assembler_vector_bd_warp_based<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreau>><<< grid, block, shared_mem >>>(vector_data, conv_data, primal_data,
                  space.get_num_dofs(), (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
                  cubature.cub_wg, cubature.num_cubs, alpha,
                  dof_mapping.cell_to_dof, dof_mapping.cell_num,
                  (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
                  (const int*) coloring_maps[col], coloring_map_sizes[col], burgers_params,
                  Intern::ViscoFunctor<DT_, MaterialType::carreau>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T},
                  loc_block_size
              );
              break;
            }
            case FEAT::VoxelAssembly::MaterialType::carreauYasuda:
            {
              VoxelAssembly::Kernel::template full_burgers_vm_assembler_vector_bd_warp_based<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreauYasuda>><<< grid, block, shared_mem >>>(vector_data, conv_data, primal_data,
                  space.get_num_dofs(), (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
                  cubature.cub_wg, cubature.num_cubs, alpha,
                  dof_mapping.cell_to_dof, dof_mapping.cell_num,
                  (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
                  (const int*) coloring_maps[col], coloring_map_sizes[col], burgers_params,
                  Intern::ViscoFunctor<DT_, MaterialType::carreauYasuda>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T,material_params.yasuda_a, material_params.mu_inf},
                  loc_block_size
              );
              break;
            }
          }
        }

        //check for cuda error (also synchronizes)
        Util::cuda_check_last_error();
      }

      template<typename Space_, typename DT_, typename IT_>
      void assemble_burgers_velo_material_csr_host([[maybe_unused]] const Space_& space,
              const CSRMatrixData<DT_, IT_>& matrix_data,
              const DT_* conv_data,
              const AssemblyCubatureData<DT_>& cubature,
              const AssemblyMappingData<DT_, IT_>& dof_mapping,
              const std::vector<int*>& coloring_maps,
              const std::vector<Index>& coloring_map_sizes,
              DT_ alpha, const AssemblyBurgersData<DT_>& burgers_params,
              const AssemblyMaterialData<DT_>& material_params, MaterialType material_type)
      {
        switch(material_type)
        {
          case MaterialType::carreau:
          {
            for(Index col = 0; col < Index(coloring_maps.size()); ++col)
            {
              VoxelAssembly::Kernel::template full_burgers_vm_assembler_matrix1_bcsr_host<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreau>,
                                                                    Intern::ViscoDFunctor<DT_, MaterialType::carreau>, FEAT::Intern::MatrixGatherScatterPolicy::useLocalSortHelper>(
                  matrix_data.data, conv_data, matrix_data.row_ptr, matrix_data.col_idx, matrix_data.num_rows, matrix_data.num_cols,
                  (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
                  cubature.cub_wg, cubature.num_cubs, alpha,
                  dof_mapping.cell_to_dof, dof_mapping.cell_num,
                  (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
                  (const int*) coloring_maps[col], coloring_map_sizes[col], dof_mapping.cell_to_dof_sorter,
                  burgers_params, material_params, Intern::ViscoFunctor<DT_, MaterialType::carreau>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T},
                  Intern::ViscoDFunctor<DT_, MaterialType::carreau>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T}
              );
            }
            break;
          }
          case MaterialType::carreauYasuda:
          {
            for(Index col = 0; col < Index(coloring_maps.size()); ++col)
            {
              VoxelAssembly::Kernel::template full_burgers_vm_assembler_matrix1_bcsr_host<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreauYasuda>,
                                                                    Intern::ViscoDFunctor<DT_, MaterialType::carreauYasuda>, FEAT::Intern::MatrixGatherScatterPolicy::useLocalSortHelper>(
                  matrix_data.data, conv_data, matrix_data.row_ptr, matrix_data.col_idx, matrix_data.num_rows, matrix_data.num_cols,
                  (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
                  cubature.cub_wg, cubature.num_cubs, alpha,
                  dof_mapping.cell_to_dof, dof_mapping.cell_num,
                  (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
                  (const int*) coloring_maps[col], coloring_map_sizes[col], dof_mapping.cell_to_dof_sorter,
                  burgers_params, material_params,
                  Intern::ViscoFunctor<DT_, MaterialType::carreauYasuda>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T,material_params.yasuda_a, material_params.mu_inf},
                  Intern::ViscoDFunctor<DT_, MaterialType::carreauYasuda>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T, material_params.yasuda_a, material_params.mu_inf}
              );
            }
            break;
          }
        }
      }

      template<typename Space_, typename DT_, typename IT_>
      void assemble_burgers_velo_material_defect_host([[maybe_unused]] const Space_& space,
              DT_* vector_data,
              const DT_* conv_data,
              const DT_* primal_data,
              const AssemblyCubatureData<DT_>& cubature,
              const AssemblyMappingData<DT_, IT_>& dof_mapping,
              const std::vector<int*>& coloring_maps,
              const std::vector<Index>& coloring_map_sizes,
              DT_ alpha, const AssemblyBurgersData<DT_>& burgers_params,
              const AssemblyMaterialData<DT_>& material_params, MaterialType material_type)
      {
        #ifdef DEBUG
        std::size_t pValue;
        cudaDeviceGetLimit(&pValue, cudaLimit::cudaLimitStackSize);
        std::cout << "Current maximum thread cache: " << pValue << std::endl;
        #endif
        switch(material_type)
        {
          case MaterialType::carreau:
          {
            for(Index col = 0; col < Index(coloring_maps.size()); ++col)
            {
              VoxelAssembly::Kernel::template full_burgers_vm_assembler_vector_bd_host<Space_, DT_, IT_>(vector_data, conv_data, primal_data,
                  space.get_num_dofs(), (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
                  cubature.cub_wg, cubature.num_cubs, alpha,
                  dof_mapping.cell_to_dof, dof_mapping.cell_num,
                  (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
                  (const int*) coloring_maps[col], coloring_map_sizes[col], burgers_params,
                  Intern::ViscoFunctor<DT_, MaterialType::carreau>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T}
              );
            }
            break;
          }
          case MaterialType::carreauYasuda:
          {
            for(Index col = 0; col < Index(coloring_maps.size()); ++col)
            {
              VoxelAssembly::Kernel::template full_burgers_vm_assembler_vector_bd_host<Space_, DT_, IT_>(vector_data, conv_data, primal_data,
                  space.get_num_dofs(), (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
                  cubature.cub_wg, cubature.num_cubs, alpha,
                  dof_mapping.cell_to_dof, dof_mapping.cell_num,
                  (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
                  (const int*) coloring_maps[col], coloring_map_sizes[col], burgers_params,
                  Intern::ViscoFunctor<DT_, MaterialType::carreauYasuda>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T,material_params.yasuda_a, material_params.mu_inf}
              );
            }
            break;
          }
        }
      }
    }
  }
}

using namespace FEAT;
using namespace FEAT::VoxelAssembly;

/*******************************************************2D implementations***************************************************/
template void Arch::assemble_burgers_velo_material_csr(const Q2StandardQuad&, const CSRMatrixData<double, std::uint32_t>&, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, double, const AssemblyBurgersData<double>&, const AssemblyMaterialData<double>&, MaterialType, int, int, int, bool);
template void Arch::assemble_burgers_velo_material_csr(const Q2StandardQuad&, const CSRMatrixData<float, std::uint32_t>&, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, float, const AssemblyBurgersData<float>&, const AssemblyMaterialData<float>&, MaterialType, int, int, int, bool);
template void Arch::assemble_burgers_velo_material_csr(const Q2StandardQuad&, const CSRMatrixData<double, std::uint64_t>&, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, double, const AssemblyBurgersData<double>&, const AssemblyMaterialData<double>&, MaterialType, int, int, int, bool);
template void Arch::assemble_burgers_velo_material_csr(const Q2StandardQuad&, const CSRMatrixData<float, std::uint64_t>&, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, float, const AssemblyBurgersData<float>&, const AssemblyMaterialData<float>&, MaterialType, int, int, int, bool);
#ifdef FEAT_HAVE_HALFMATH
template void Arch::assemble_burgers_velo_material_csr(const Q2StandardQuad&, const CSRMatrixData<Half, std::uint32_t>&, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, Half, const AssemblyBurgersData<Half>&, const AssemblyMaterialData<Half>&, MaterialType, int, int, int, bool);
template void Arch::assemble_burgers_velo_material_csr(const Q2StandardQuad&, const CSRMatrixData<Half, std::uint64_t>&, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, Half, const AssemblyBurgersData<Half>&, const AssemblyMaterialData<Half>&, MaterialType, int, int, int, bool);
#endif

template void Arch::assemble_burgers_velo_material_csr_host(const Q2StandardQuad&, const CSRMatrixData<double, std::uint32_t>&, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, double, const AssemblyBurgersData<double>&, const AssemblyMaterialData<double>&, MaterialType);
template void Arch::assemble_burgers_velo_material_csr_host(const Q2StandardQuad&, const CSRMatrixData<float, std::uint32_t>&, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, float, const AssemblyBurgersData<float>&, const AssemblyMaterialData<float>&, MaterialType);
template void Arch::assemble_burgers_velo_material_csr_host(const Q2StandardQuad&, const CSRMatrixData<double, std::uint64_t>&, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, double, const AssemblyBurgersData<double>&, const AssemblyMaterialData<double>&, MaterialType);
template void Arch::assemble_burgers_velo_material_csr_host(const Q2StandardQuad&, const CSRMatrixData<float, std::uint64_t>&, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, float, const AssemblyBurgersData<float>&, const AssemblyMaterialData<float>&, MaterialType);
#ifdef FEAT_HAVE_HALFMATH
template void Arch::assemble_burgers_velo_material_csr_host(const Q2StandardQuad&, const CSRMatrixData<Half, std::uint32_t>&, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, Half, const AssemblyBurgersData<Half>&, const AssemblyMaterialData<Half>&, MaterialType);
template void Arch::assemble_burgers_velo_material_csr_host(const Q2StandardQuad&, const CSRMatrixData<Half, std::uint64_t>&, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, Half, const AssemblyBurgersData<Half>&, const AssemblyMaterialData<Half>&, MaterialType);
#endif

template void Arch::assemble_burgers_velo_material_defect(const Q2StandardQuad&, double*, const double*, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, double, const AssemblyBurgersData<double>&, const AssemblyMaterialData<double>&, MaterialType, int, int, int, bool);
template void Arch::assemble_burgers_velo_material_defect(const Q2StandardQuad&, float*, const float*, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, float, const AssemblyBurgersData<float>&, const AssemblyMaterialData<float>&, MaterialType, int, int, int, bool);
template void Arch::assemble_burgers_velo_material_defect(const Q2StandardQuad&, double*, const double*, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, double, const AssemblyBurgersData<double>&, const AssemblyMaterialData<double>&, MaterialType, int, int, int, bool);
template void Arch::assemble_burgers_velo_material_defect(const Q2StandardQuad&, float*, const float*, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, float, const AssemblyBurgersData<float>&, const AssemblyMaterialData<float>&, MaterialType, int, int, int, bool);
#ifdef FEAT_HAVE_HALFMATH
template void Arch::assemble_burgers_velo_material_defect(const Q2StandardQuad&, Half*, const Half*, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, Half, const AssemblyBurgersData<Half>&, const AssemblyMaterialData<Half>&, MaterialType, int, int, int, bool);
template void Arch::assemble_burgers_velo_material_defect(const Q2StandardQuad&, Half*, const Half*, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, Half, const AssemblyBurgersData<Half>&, const AssemblyMaterialData<Half>&, MaterialType, int, int, int, bool);
#endif

template void Arch::assemble_burgers_velo_material_defect_host(const Q2StandardQuad&, double*, const double*, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, double, const AssemblyBurgersData<double>&, const AssemblyMaterialData<double>&, MaterialType);
template void Arch::assemble_burgers_velo_material_defect_host(const Q2StandardQuad&, float*, const float*, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, float, const AssemblyBurgersData<float>&, const AssemblyMaterialData<float>&, MaterialType);
template void Arch::assemble_burgers_velo_material_defect_host(const Q2StandardQuad&, double*, const double*, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, double, const AssemblyBurgersData<double>&, const AssemblyMaterialData<double>&, MaterialType);
template void Arch::assemble_burgers_velo_material_defect_host(const Q2StandardQuad&, float*, const float*, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, float, const AssemblyBurgersData<float>&, const AssemblyMaterialData<float>&, MaterialType);
#ifdef FEAT_HAVE_HALFMATH
template void Arch::assemble_burgers_velo_material_defect_host(const Q2StandardQuad&, Half*, const Half*, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, Half, const AssemblyBurgersData<Half>&, const AssemblyMaterialData<Half>&, MaterialType);
template void Arch::assemble_burgers_velo_material_defect_host(const Q2StandardQuad&, Half*, const Half*, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, Half, const AssemblyBurgersData<Half>&, const AssemblyMaterialData<Half>&, MaterialType);
#endif

/*********************************************************3D implementations**************************************************************************************/

template void Arch::assemble_burgers_velo_material_csr(const Q2StandardHexa&, const CSRMatrixData<double, std::uint32_t>&, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, double, const AssemblyBurgersData<double>&, const AssemblyMaterialData<double>&, MaterialType, int, int, int, bool);
template void Arch::assemble_burgers_velo_material_csr(const Q2StandardHexa&, const CSRMatrixData<float, std::uint32_t>&, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, float, const AssemblyBurgersData<float>&, const AssemblyMaterialData<float>&, MaterialType, int, int, int, bool);
template void Arch::assemble_burgers_velo_material_csr(const Q2StandardHexa&, const CSRMatrixData<double, std::uint64_t>&, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, double, const AssemblyBurgersData<double>&, const AssemblyMaterialData<double>&, MaterialType, int, int, int, bool);
template void Arch::assemble_burgers_velo_material_csr(const Q2StandardHexa&, const CSRMatrixData<float, std::uint64_t>&, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, float, const AssemblyBurgersData<float>&, const AssemblyMaterialData<float>&, MaterialType, int, int, int, bool);
#ifdef FEAT_HAVE_HALFMATH
template void Arch::assemble_burgers_velo_material_csr(const Q2StandardHexa&, const CSRMatrixData<Half, std::uint32_t>&, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, Half, const AssemblyBurgersData<Half>&, const AssemblyMaterialData<Half>&, MaterialType, int, int, int, bool);
template void Arch::assemble_burgers_velo_material_csr(const Q2StandardHexa&, const CSRMatrixData<Half, std::uint64_t>&, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, Half, const AssemblyBurgersData<Half>&, const AssemblyMaterialData<Half>&, MaterialType, int, int, int, bool);
#endif

template void Arch::assemble_burgers_velo_material_csr_host(const Q2StandardHexa&, const CSRMatrixData<double, std::uint32_t>&, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, double, const AssemblyBurgersData<double>&, const AssemblyMaterialData<double>&, MaterialType);
template void Arch::assemble_burgers_velo_material_csr_host(const Q2StandardHexa&, const CSRMatrixData<float, std::uint32_t>&, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, float, const AssemblyBurgersData<float>&, const AssemblyMaterialData<float>&, MaterialType);
template void Arch::assemble_burgers_velo_material_csr_host(const Q2StandardHexa&, const CSRMatrixData<double, std::uint64_t>&, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, double, const AssemblyBurgersData<double>&, const AssemblyMaterialData<double>&, MaterialType);
template void Arch::assemble_burgers_velo_material_csr_host(const Q2StandardHexa&, const CSRMatrixData<float, std::uint64_t>&, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, float, const AssemblyBurgersData<float>&, const AssemblyMaterialData<float>&, MaterialType);
#ifdef FEAT_HAVE_HALFMATH
template void Arch::assemble_burgers_velo_material_csr_host(const Q2StandardHexa&, const CSRMatrixData<Half, std::uint32_t>&, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, Half, const AssemblyBurgersData<Half>&, const AssemblyMaterialData<Half>&, MaterialType);
template void Arch::assemble_burgers_velo_material_csr_host(const Q2StandardHexa&, const CSRMatrixData<Half, std::uint64_t>&, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, Half, const AssemblyBurgersData<Half>&, const AssemblyMaterialData<Half>&, MaterialType);
#endif

template void Arch::assemble_burgers_velo_material_defect(const Q2StandardHexa&, double*, const double*, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, double, const AssemblyBurgersData<double>&, const AssemblyMaterialData<double>&, MaterialType, int, int, int, bool);
template void Arch::assemble_burgers_velo_material_defect(const Q2StandardHexa&, float*, const float*, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, float, const AssemblyBurgersData<float>&, const AssemblyMaterialData<float>&, MaterialType, int, int, int, bool);
template void Arch::assemble_burgers_velo_material_defect(const Q2StandardHexa&, double*, const double*, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, double, const AssemblyBurgersData<double>&, const AssemblyMaterialData<double>&, MaterialType, int, int, int, bool);
template void Arch::assemble_burgers_velo_material_defect(const Q2StandardHexa&, float*, const float*, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, float, const AssemblyBurgersData<float>&, const AssemblyMaterialData<float>&, MaterialType, int, int, int, bool);
#ifdef FEAT_HAVE_HALFMATH
template void Arch::assemble_burgers_velo_material_defect(const Q2StandardHexa&, Half*, const Half*, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, Half, const AssemblyBurgersData<Half>&, const AssemblyMaterialData<Half>&, MaterialType, int, int, int, bool);
template void Arch::assemble_burgers_velo_material_defect(const Q2StandardHexa&, Half*, const Half*, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, Half, const AssemblyBurgersData<Half>&, const AssemblyMaterialData<Half>&, MaterialType, int, int, int, bool);
#endif

template void Arch::assemble_burgers_velo_material_defect_host(const Q2StandardHexa&, double*, const double*, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, double, const AssemblyBurgersData<double>&, const AssemblyMaterialData<double>&, MaterialType);
template void Arch::assemble_burgers_velo_material_defect_host(const Q2StandardHexa&, float*, const float*, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, float, const AssemblyBurgersData<float>&, const AssemblyMaterialData<float>&, MaterialType);
template void Arch::assemble_burgers_velo_material_defect_host(const Q2StandardHexa&, double*, const double*, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, double, const AssemblyBurgersData<double>&, const AssemblyMaterialData<double>&, MaterialType);
template void Arch::assemble_burgers_velo_material_defect_host(const Q2StandardHexa&, float*, const float*, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, float, const AssemblyBurgersData<float>&, const AssemblyMaterialData<float>&, MaterialType);
#ifdef FEAT_HAVE_HALFMATH
template void Arch::assemble_burgers_velo_material_defect_host(const Q2StandardHexa&, Half*, const Half*, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint32_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, Half, const AssemblyBurgersData<Half>&, const AssemblyMaterialData<Half>&, MaterialType);
template void Arch::assemble_burgers_velo_material_defect_host(const Q2StandardHexa&, Half*, const Half*, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint64_t>&,
                                          const std::vector<int*>&, const std::vector<Index>&, Half, const AssemblyBurgersData<Half>&, const AssemblyMaterialData<Half>&, MaterialType);
#endif
