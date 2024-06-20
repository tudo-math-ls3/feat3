#include <kernel/base_header.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/voxel_assembly/burgers_assembler.hpp>
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

      template<typename Space_, typename DT_, typename IT_, FEAT::Intern::MatrixGatherScatterPolicy pol_ = FEAT::Intern::MatrixGatherScatterPolicy::useLocalOps>
      __global__ void full_burgers_assembler_matrix1_bcsr(DT_* __restrict__ matrix_data, const DT_* __restrict__ conv_data,
                const IT_* __restrict__  matrix_row_ptr, const IT_* __restrict__ matrix_col_idx, Index matrix_num_rows, Index matrix_num_cols,
                const Tiny::Vector<DT_, Space_::world_dim>* __restrict__ cub_pt,
                const DT_* __restrict__  cub_wg, int num_cubs, DT_ alpha,
                const IT_* __restrict__ cell_to_dof, [[maybe_unused]] Index cell_num,
                const Tiny::Vector<DT_, Space_::world_dim>* __restrict__ nodes, [[maybe_unused]] Index node_size,
                const int* __restrict__ coloring_map, Index coloring_size, const IT_* __restrict__ cell_to_dof_sorter,
                const VoxelAssembly::AssemblyBurgersData<DT_> burgers_params)
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
        DataType local_coeffs[dim][num_loc_verts];
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
        if(need_convection || need_streamline) //need stream diff or convection?
        {
          LAFEM::template VectorGatherScatterHelper<SpaceType, DataType, IndexType>::gather_vector_dense(local_conv_dofs,
                    (const VecValueType*)conv_data, IndexType(matrix_num_rows), local_dofs,DataType(1));
        }

        VoxelAssembly::Kernel::burgers_mat_assembly_kernel<SpaceHelp>(loc_mat, local_conv_dofs, local_coeffs, cub_pt, cub_wg, num_cubs, burgers_params,
                                    need_streamline, need_convection, tol_eps);

        //scatter
        LAFEM::template MatrixGatherScatterHelper<SpaceType, DataType, IndexType, pol_>::scatter_matrix_csr(loc_mat, (MatValueType*)matrix_data, local_dofs, local_dofs, IndexType(matrix_num_rows), IndexType(matrix_num_cols), matrix_row_ptr, matrix_col_idx, alpha, local_dof_sorter);

      }

      template<typename Space_, typename DT_, typename IT_>
      __global__ void full_burgers_assembler_vector_bd(DT_* __restrict__ vector_data,
                const DT_* conv_data, const DT_* primal_data, Index vec_size,
                const Tiny::Vector<DT_, Space_::world_dim>* __restrict__ cub_pt,
                const DT_* __restrict__  cub_wg, int num_cubs, DT_ alpha,
                const IT_* __restrict__ cell_to_dof, [[maybe_unused]] Index cell_num,
                const Tiny::Vector<DT_, Space_::world_dim>* __restrict__ nodes, [[maybe_unused]] Index node_size,
                const int* __restrict__ coloring_map, Index coloring_size,
                const VoxelAssembly::AssemblyBurgersData<DT_> burgers_params)
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
        DataType local_coeffs[dim][num_loc_verts];

        //now do work for this cell
        Index cell = Index(coloring_map[idx]);
        const IndexSetWrapper<IndexType> local_dofs_w{cell_to_dof, IndexType(num_loc_dofs)};
        const IndexType* local_dofs = cell_to_dof + cell*num_loc_dofs;
        SpaceHelp::set_coefficients(local_coeffs, local_dofs_w, nodes, cell);
        LAFEM::template VectorGatherScatterHelper<SpaceType, DataType, IndexType>::gather_vector_dense(local_prim_dofs,
                  (const VecValueType*)primal_data, IndexType(vec_size), local_dofs, DataType(1));

        //if we need to, gather local convection vector
        if(need_convection) //need stream diff or convection?
        {
          LAFEM::template VectorGatherScatterHelper<SpaceType, DataType, IndexType>::gather_vector_dense(local_conv_dofs,
                    (const VecValueType*)conv_data, IndexType(vec_size), local_dofs, DataType(1));
        }

        VoxelAssembly::Kernel::burgers_defect_assembly_kernel<SpaceHelp>(loc_vec, local_prim_dofs, local_conv_dofs, local_coeffs, cub_pt, cub_wg, num_cubs,
                                                      burgers_params, need_convection);
        //scatter
        LAFEM::template VectorGatherScatterHelper<SpaceType, DataType, IndexType>::scatter_vector_dense(loc_vec,
                  (VecValueType*)vector_data, IndexType(vec_size), local_dofs, alpha);

      }

      /**
      * \brief Reduces the max local vector norm of a convetion vector.
      *
      * This uses multiple reductions with red black sumations on a single block and an atomic add at the end to reduce
      * the vector value.
      *
      * \attention Only one dimensional kernel and blockDim.x * num_kernels >= vec_size required.
      * \warning Will fail if if blockdim is not a multiple of 2 (which you should do in any case)
      *
      * \param[in] convect Device memory to the convection vector in native view.
      * \param[out] result Ptr to device global variable.
      * \param[in] vec_size The size of the convection vector in native view.
      */
      template<typename DT_, int dim_>
      __global__ void set_sd_v_norm(const Tiny::Vector<DT_, dim_>* __restrict__ convect, DT_* __restrict__ result, Index vec_size)
      {
        //are we a power of two?
        #ifdef DEBUG
        assert(int((blockDim.x & (blockDim.x-1)) == 0));
        #endif
        // get our block dim and thread id
        const int lidx = threadIdx.x;
        const int gidx = lidx + blockIdx.x * blockDim.x;
        // use a shared array, to extract the max values
        // size is determined depending on blocksize, hence extern... You have to provide correct shared memory amount in kernel
        // one caveat, this initilizes the extern parameter partial max for each template instantiated, which does not work
        // if these have different datatypes. For this reason, we have to declare a base array of correct size and then reintepret the cast
        extern __shared__ __align__(8) unsigned char partial_max_shared[];
        DT_* partial_max = reinterpret_cast<DT_*>(partial_max_shared);
        if(gidx < int(vec_size))
        {
          // extract local array
          partial_max[lidx] = convect[gidx].norm_euclid();
        }
        else
        {
          //extract zeros, so we do not have to branch in the next loop
          partial_max[lidx] = DT_(0);
        }

        //and now reduce over half the array size each time
        for(int stride = blockDim.x / 2; stride > 0; stride >>= 1)
        {
          //synchronize all threads in block, to guarentee that all values are written
          __syncthreads();
          if(lidx < stride)
          {
            partial_max[lidx] = CudaMath::cuda_max(partial_max[lidx], partial_max[lidx+stride]);
          }
        }

        //and now use atomic Operation to synronize value in 0 over the whole device
        if(lidx == 0)
          CudaMath::cuda_atomic_max(result, partial_max[0]);
        //and done
      }

      /**************************************************************************************************************/
      /*                                       CUDA Host OMP Kernels                                                */
      /**************************************************************************************************************/

      template<typename Space_, typename DT_, typename IT_, FEAT::Intern::MatrixGatherScatterPolicy pol_ = FEAT::Intern::MatrixGatherScatterPolicy::useLocalOps>
      void full_burgers_assembler_matrix1_bcsr_host(DT_* matrix_data, const DT_* conv_data,
                const IT_*  matrix_row_ptr, const IT_* matrix_col_idx, Index matrix_num_rows, Index matrix_num_cols,
                const Tiny::Vector<DT_, Space_::world_dim>* cub_pt,
                const DT_*  cub_wg, int num_cubs, DT_ alpha,
                const IT_* cell_to_dof, [[maybe_unused]] Index cell_num,
                const Tiny::Vector<DT_, Space_::world_dim>* nodes, [[maybe_unused]] Index node_size,
                const int* coloring_map, Index coloring_size, const IT_* cell_to_dof_sorter,
                const VoxelAssembly::AssemblyBurgersData<DT_>& burgers_params)
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
          DataType local_coeffs[dim][num_loc_verts];
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
          if(need_convection || need_streamline) //need stream diff or convection?
          {
            LAFEM::template VectorGatherScatterHelper<SpaceType, DataType, IndexType>::gather_vector_dense(local_conv_dofs,
                      (const VecValueType*)conv_data, IndexType(matrix_num_rows), local_dofs,DataType(1));
          }

          VoxelAssembly::Kernel::burgers_mat_assembly_kernel<SpaceHelp>(loc_mat, local_conv_dofs, local_coeffs, cub_pt, cub_wg, num_cubs, burgers_params,
                                      need_streamline, need_convection, tol_eps);

          //scatter
          LAFEM::template MatrixGatherScatterHelper<SpaceType, DataType, IndexType, pol_>::scatter_matrix_csr(loc_mat, (MatValueType*)matrix_data, local_dofs, local_dofs, IndexType(matrix_num_rows), IndexType(matrix_num_cols), matrix_row_ptr, matrix_col_idx, alpha, local_dof_sorter);

        }
      }

      template<typename Space_, typename DT_, typename IT_>
      void full_burgers_assembler_vector_bd_host(DT_* vector_data,
                const DT_* conv_data, const DT_* primal_data, Index vec_size,
                const Tiny::Vector<DT_, Space_::world_dim>* cub_pt,
                const DT_*  cub_wg, int num_cubs, DT_ alpha,
                const IT_* cell_to_dof, [[maybe_unused]] Index cell_num,
                const Tiny::Vector<DT_, Space_::world_dim>* nodes, [[maybe_unused]] Index node_size,
                const int* coloring_map, Index coloring_size,
                const VoxelAssembly::AssemblyBurgersData<DT_>& burgers_params)
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
          DataType local_coeffs[dim][num_loc_verts];

          //now do work for this cell
          Index cell = Index(coloring_map[idx]);
          const IndexSetWrapper<IndexType> local_dofs_w{cell_to_dof, IndexType(num_loc_dofs)};
          const IndexType* local_dofs = cell_to_dof + cell*num_loc_dofs;
          SpaceHelp::set_coefficients(local_coeffs, local_dofs_w, nodes, cell);
          LAFEM::template VectorGatherScatterHelper<SpaceType, DataType, IndexType>::gather_vector_dense(local_prim_dofs,
                    (const VecValueType*)primal_data, IndexType(vec_size), local_dofs, DataType(1));

          //if we need to, gather local convection vector
          if(need_convection) //need stream diff or convection?
          {
            LAFEM::template VectorGatherScatterHelper<SpaceType, DataType, IndexType>::gather_vector_dense(local_conv_dofs,
                      (const VecValueType*)conv_data, IndexType(vec_size), local_dofs, DataType(1));
          }

          VoxelAssembly::Kernel::burgers_defect_assembly_kernel<SpaceHelp>(loc_vec, local_prim_dofs, local_conv_dofs, local_coeffs, cub_pt, cub_wg, num_cubs,
                                                         burgers_params, need_convection);
          //scatter
          LAFEM::template VectorGatherScatterHelper<SpaceType, DataType, IndexType>::scatter_vector_dense(loc_vec,
                    (VecValueType*)vector_data, IndexType(vec_size), local_dofs, alpha);

        }
      }

      /**
      * \brief Reduces the max local vector norm of a convection vector.
      *
      * This uses multiple reductions with red black sumations on a single block and an atomic add at the end to reduce
      * the vector value.
      *
      * \attention Only one dimensional kernel and blockDim.x * num_kernels >= vec_size required.
      * \warning Will fail if if blockdim is not a multiple of 2 (which you should do in any case)
      *
      * \param[in] convect Device memory to the convection vector in native view.
      * \param[out] result Ptr to device global variable.
      * \param[in] vec_size The size of the convection vector in native view.
      */
      template<typename DT_, int dim_>
      void set_sd_v_norm_host(const Tiny::Vector<DT_, dim_>* convect, DT_* result, Index vec_size)
      {
        #ifdef FEAT_HAVE_OMP
        DT_ max_val(DT_(0));
        // since we potentially use Half values, we do the max reduction ourselfes
        //simply use a local array of size 128... if we have more available threads, we wont use them...
        // this is of course not future proof, maybe solve this with a compiletime variable?
        DT_ max_vals[128];
        // parallel region with at most 128 threads
        #pragma omp parallel num_threads(128)
        {
          const int num_threads = omp_get_num_threads();
          const int thread_id = omp_get_thread_num();
          max_vals[thread_id] = DT_(0);
          #pragma omp for
          for(int i = 0; i < int(vec_size); ++i)
          {
            //synchronize all threads in block, to guarentee that all values are written
            max_vals[thread_id] = CudaMath::cuda_max(convect[i].norm_euclid(), max_vals[thread_id]);
          }

          #pragma omp single
          max_val = std::reduce(max_vals, max_vals + num_threads, DT_(0), [](const DT_& a, const DT_& b){return CudaMath::cuda_max(a,b);});

        }
        //and overwrite
        *result = max_val;
        #else
        DT_ max_val(DT_(0));
        // since we potentially use Half values, we do the max reduction ourselfes
        for(int i = 0; i < int(vec_size); ++i)
        {
          //synchronize all threads in block, to guarentee that all values are written
          max_val = CudaMath::cuda_max(convect[i].norm_euclid(), max_val);
        }
        //and overwrite
        *result = max_val;
        #endif

      }


    }

    namespace Arch
    {
      template<typename Space_, typename DT_, typename IT_>
      void assemble_burgers_csr([[maybe_unused]] const Space_& space,
              const CSRMatrixData<DT_, IT_>& matrix_data,
              const DT_* conv_data,
              const AssemblyCubatureData<DT_>& cubature,
              const AssemblyMappingData<DT_, IT_>& dof_mapping,
              const std::vector<std::vector<int>>& coloring_maps_host,
              const std::vector<void*>& coloring_maps_device,
              DT_ alpha, const AssemblyBurgersData<DT_>& burgers_params)
      {
        const Index blocksize = Util::cuda_blocksize_blocked_assembly;
        // const Index blocksize = 64;

        for(Index col = 0; col < coloring_maps_device.size(); ++col)
        {
          dim3 grid;
          dim3 block;
          block.x = (unsigned int)blocksize;
          grid.x = (unsigned int)ceil(double(coloring_maps_host[col].size())/double(block.x));

          //kernel call, since this uses the standard stream, sync before next call is enforced:
          VoxelAssembly::Kernel::template full_burgers_assembler_matrix1_bcsr<Space_, DT_, IT_, FEAT::Intern::MatrixGatherScatterPolicy::useLocalSortHelper><<< grid, block >>>(
              matrix_data.data, conv_data, matrix_data.row_ptr, matrix_data.col_idx, matrix_data.num_rows, matrix_data.num_cols,
              (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
              cubature.cub_wg, cubature.num_cubs, alpha,
              dof_mapping.cell_to_dof, dof_mapping.cell_num,
              (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
              (const int*) coloring_maps_device[col], coloring_maps_host[col].size(), dof_mapping.cell_to_dof_sorter,
              burgers_params
          );
        }
        //check for cuda error in our kernel
        Util::cuda_check_last_error();
      }


      template<typename Space_, typename DT_, typename IT_>
      void assemble_burgers_defect([[maybe_unused]] const Space_& space,
              DT_* vector_data,
              const DT_* conv_data,
              const DT_* primal_data,
              const AssemblyCubatureData<DT_>& cubature,
              const AssemblyMappingData<DT_, IT_>& dof_mapping,
              const std::vector<std::vector<int>>& coloring_maps_host,
              const std::vector<void*>& coloring_maps_device,
              DT_ alpha, const AssemblyBurgersData<DT_>& burgers_params)
      {
        const Index blocksize = Util::cuda_blocksize_scalar_assembly;
        // const Index blocksize = 64;
        #ifdef DEBUG
        std::size_t pValue;
        cudaDeviceGetLimit(&pValue, cudaLimit::cudaLimitStackSize);
        std::cout << "Current maximum thread cache: " << pValue << std::endl;
        #endif

        for(Index col = 0; col < coloring_maps_device.size(); ++col)
        {
          dim3 grid;
          dim3 block;
          block.x = (unsigned int)blocksize;
          grid.x = (unsigned int)ceil(double(coloring_maps_host[col].size())/double(block.x));

          VoxelAssembly::Kernel::template full_burgers_assembler_vector_bd<Space_, DT_, IT_><<< grid, block >>>(vector_data, conv_data, primal_data,
              space.get_num_dofs(), (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
              cubature.cub_wg, cubature.num_cubs, alpha,
              dof_mapping.cell_to_dof, dof_mapping.cell_num,
              (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
              (const int*) coloring_maps_device[col], coloring_maps_host[col].size(), burgers_params
          );
        }

        //check for cuda error (also synchronizes)
        Util::cuda_check_last_error();
      }

      // attention: This requires at max only one device per mpi rank!
      template<typename DT_, typename IT_, int dim_>
      DT_ get_sd_v_norm(const LAFEM::DenseVectorBlocked<DT_, IT_, dim_>& convect)
      {
        const Index blocksize = Util::cuda_blocksize_reduction;
        dim3 grid;
        dim3 block;
        block.x = (unsigned int)blocksize;
        grid.x = (unsigned int)(ceil(convect.template size<LAFEM::Perspective::native>()/double(block.x)));

        //extract pointer
        const Tiny::Vector<DT_, dim_>* vec_data = (const Tiny::Vector<DT_, dim_>*)convect.template elements<LAFEM::Perspective::native>();
        //init result device global variable
        DT_* glob_res = (DT_*)Util::cuda_malloc(sizeof(DT_));
        //init value
        Util::cuda_set_memory(glob_res, DT_(0), 1);

        VoxelAssembly::Kernel::template set_sd_v_norm<DT_, dim_><<<grid, block, blocksize*sizeof(DT_)>>>(vec_data, glob_res, convect.template size<LAFEM::Perspective::native>());

        //check for cuda error in our kernel
        Util::cuda_check_last_error();

        // retrieve value from device
        DT_ ret_val;
        Util::cuda_copy((void*)&ret_val, (void*)glob_res, sizeof(DT_));
        Util::cuda_free((void*)glob_res);

        return ret_val;
      }

      template<typename DT_, typename IT_, int dim_>
      DT_ get_sd_v_norm(const Global::Vector<LAFEM::DenseVectorBlocked<DT_, IT_, dim_>, LAFEM::VectorMirror<DT_, IT_>>& convect)
      {
        auto local_norm = get_sd_v_norm(convect.local());
        const auto* gate = convect.get_gate();
        if(gate != nullptr)
        {
          local_norm = gate->max(local_norm);
        }
        return local_norm;
      }

      template<typename Space_, typename DT_, typename IT_>
      void assemble_burgers_csr_host([[maybe_unused]] const Space_& space,
              const CSRMatrixData<DT_, IT_>& matrix_data,
              const DT_* conv_data,
              const AssemblyCubatureData<DT_>& cubature,
              const AssemblyMappingData<DT_, IT_>& dof_mapping,
              const std::vector<std::vector<int>>& coloring_maps_host,
              [[maybe_unused]] const std::vector<void*>& coloring_maps_device,
              DT_ alpha, const AssemblyBurgersData<DT_>& burgers_params)
      {
        for(Index col = 0; col < Index(coloring_maps_host.size()); ++col)
        {
          VoxelAssembly::Kernel::template full_burgers_assembler_matrix1_bcsr_host<Space_, DT_, IT_, FEAT::Intern::MatrixGatherScatterPolicy::useLocalSortHelper>(
              matrix_data.data, conv_data, matrix_data.row_ptr, matrix_data.col_idx, matrix_data.num_rows, matrix_data.num_cols,
              (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
              cubature.cub_wg, cubature.num_cubs, alpha,
              dof_mapping.cell_to_dof, dof_mapping.cell_num,
              (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
              (const int*) coloring_maps_host[col].data(), coloring_maps_host[col].size(), dof_mapping.cell_to_dof_sorter,
              burgers_params
          );
        }
      }

      template<typename Space_, typename DT_, typename IT_>
      void assemble_burgers_defect_host([[maybe_unused]] const Space_& space,
              DT_* vector_data,
              const DT_* conv_data,
              const DT_* primal_data,
              const AssemblyCubatureData<DT_>& cubature,
              const AssemblyMappingData<DT_, IT_>& dof_mapping,
              const std::vector<std::vector<int>>& coloring_maps_host,
              [[maybe_unused]] const std::vector<void*>& coloring_maps_device,
              DT_ alpha, const AssemblyBurgersData<DT_>& burgers_params)
      {
        #ifdef DEBUG
        std::size_t pValue;
        cudaDeviceGetLimit(&pValue, cudaLimit::cudaLimitStackSize);
        std::cout << "Current maximum thread cache: " << pValue << std::endl;
        #endif
        for(Index col = 0; col < Index(coloring_maps_host.size()); ++col)
        {
          VoxelAssembly::Kernel::template full_burgers_assembler_vector_bd_host<Space_, DT_, IT_>(vector_data, conv_data, primal_data,
              space.get_num_dofs(), (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
              cubature.cub_wg, cubature.num_cubs, alpha,
              dof_mapping.cell_to_dof, dof_mapping.cell_num,
              (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
              (const int*) coloring_maps_host[col].data(), coloring_maps_host[col].size(), burgers_params
          );
        }
      }


      template<typename DT_, typename IT_, int dim_>
      DT_ get_sd_v_norm_host(const LAFEM::DenseVectorBlocked<DT_, IT_, dim_>& convect)
      {
        //extract pointer
        const Tiny::Vector<DT_, dim_>* vec_data = (const Tiny::Vector<DT_, dim_>*)convect.template elements<LAFEM::Perspective::native>();

        DT_ glob_res = DT_(0);

        VoxelAssembly::Kernel::template set_sd_v_norm_host<DT_, dim_>(vec_data, &glob_res, convect.template size<LAFEM::Perspective::native>());


        return glob_res;
      }

      template<typename DT_, typename IT_, int dim_>
      DT_ get_sd_v_norm_host(const Global::Vector<LAFEM::DenseVectorBlocked<DT_, IT_, dim_>, LAFEM::VectorMirror<DT_, IT_>>& convect)
      {
        auto local_norm = get_sd_v_norm_host(convect.local());
        const auto* gate = convect.get_gate();
        if(gate != nullptr)
        {
          local_norm = gate->max(local_norm);
        }
        return local_norm;
      }

    }
  }
}

using namespace FEAT;
using namespace FEAT::VoxelAssembly;

/*******************************************************2D implementations***************************************************/
template void Arch::assemble_burgers_csr(const Q2StandardQuad&, const CSRMatrixData<double, std::uint32_t>&, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double, const AssemblyBurgersData<double>&);
template void Arch::assemble_burgers_csr(const Q2StandardQuad&, const CSRMatrixData<float, std::uint32_t>&, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float, const AssemblyBurgersData<float>&);
template void Arch::assemble_burgers_csr(const Q2StandardQuad&, const CSRMatrixData<double, std::uint64_t>&, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double, const AssemblyBurgersData<double>&);
template void Arch::assemble_burgers_csr(const Q2StandardQuad&, const CSRMatrixData<float, std::uint64_t>&, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float, const AssemblyBurgersData<float>&);
#ifdef FEAT_HAVE_HALFMATH
template void Arch::assemble_burgers_csr(const Q2StandardQuad&, const CSRMatrixData<Half, std::uint32_t>&, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, Half, const AssemblyBurgersData<Half>&);
template void Arch::assemble_burgers_csr(const Q2StandardQuad&, const CSRMatrixData<Half, std::uint64_t>&, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, Half, const AssemblyBurgersData<Half>&);
#endif

template void Arch::assemble_burgers_csr_host(const Q2StandardQuad&, const CSRMatrixData<double, std::uint32_t>&, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double, const AssemblyBurgersData<double>&);
template void Arch::assemble_burgers_csr_host(const Q2StandardQuad&, const CSRMatrixData<float, std::uint32_t>&, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float, const AssemblyBurgersData<float>&);
template void Arch::assemble_burgers_csr_host(const Q2StandardQuad&, const CSRMatrixData<double, std::uint64_t>&, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double, const AssemblyBurgersData<double>&);
template void Arch::assemble_burgers_csr_host(const Q2StandardQuad&, const CSRMatrixData<float, std::uint64_t>&, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float, const AssemblyBurgersData<float>&);
#ifdef FEAT_HAVE_HALFMATH
template void Arch::assemble_burgers_csr_host(const Q2StandardQuad&, const CSRMatrixData<Half, std::uint32_t>&, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, Half, const AssemblyBurgersData<Half>&);
template void Arch::assemble_burgers_csr_host(const Q2StandardQuad&, const CSRMatrixData<Half, std::uint64_t>&, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, Half, const AssemblyBurgersData<Half>&);
#endif

template void Arch::assemble_burgers_defect(const Q2StandardQuad&, double*, const double*, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double, const AssemblyBurgersData<double>&);
template void Arch::assemble_burgers_defect(const Q2StandardQuad&, float*, const float*, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float, const AssemblyBurgersData<float>&);
template void Arch::assemble_burgers_defect(const Q2StandardQuad&, double*, const double*, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double, const AssemblyBurgersData<double>&);
template void Arch::assemble_burgers_defect(const Q2StandardQuad&, float*, const float*, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float, const AssemblyBurgersData<float>&);
#ifdef FEAT_HAVE_HALFMATH
template void Arch::assemble_burgers_defect(const Q2StandardQuad&, Half*, const Half*, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, Half, const AssemblyBurgersData<Half>&);
template void Arch::assemble_burgers_defect(const Q2StandardQuad&, Half*, const Half*, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, Half, const AssemblyBurgersData<Half>&);
#endif

template void Arch::assemble_burgers_defect_host(const Q2StandardQuad&, double*, const double*, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double, const AssemblyBurgersData<double>&);
template void Arch::assemble_burgers_defect_host(const Q2StandardQuad&, float*, const float*, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float, const AssemblyBurgersData<float>&);
template void Arch::assemble_burgers_defect_host(const Q2StandardQuad&, double*, const double*, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double, const AssemblyBurgersData<double>&);
template void Arch::assemble_burgers_defect_host(const Q2StandardQuad&, float*, const float*, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float, const AssemblyBurgersData<float>&);
#ifdef FEAT_HAVE_HALFMATH
template void Arch::assemble_burgers_defect_host(const Q2StandardQuad&, Half*, const Half*, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, Half, const AssemblyBurgersData<Half>&);
template void Arch::assemble_burgers_defect_host(const Q2StandardQuad&, Half*, const Half*, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, Half, const AssemblyBurgersData<Half>&);
#endif

template double Arch::get_sd_v_norm(const LAFEM::DenseVectorBlocked<double, std::uint32_t, 2>&);
template float Arch::get_sd_v_norm(const LAFEM::DenseVectorBlocked<float, std::uint32_t, 2>&);
template double Arch::get_sd_v_norm(const LAFEM::DenseVectorBlocked<double, std::uint64_t, 2>&);
template float Arch::get_sd_v_norm(const LAFEM::DenseVectorBlocked<float, std::uint64_t, 2>&);
#ifdef FEAT_HAVE_HALFMATH
template Half Arch::get_sd_v_norm(const LAFEM::DenseVectorBlocked<Half, std::uint32_t, 2>&);
template Half Arch::get_sd_v_norm(const LAFEM::DenseVectorBlocked<Half, std::uint64_t, 2>&);
#endif

template double Arch::get_sd_v_norm(const Global::Vector<LAFEM::DenseVectorBlocked<double, std::uint32_t, 2>, LAFEM::VectorMirror<double, std::uint32_t>>&);
template float Arch::get_sd_v_norm(const Global::Vector<LAFEM::DenseVectorBlocked<float, std::uint32_t, 2>, LAFEM::VectorMirror<float, std::uint32_t>>&);
template double Arch::get_sd_v_norm(const Global::Vector<LAFEM::DenseVectorBlocked<double, std::uint64_t, 2>, LAFEM::VectorMirror<double, std::uint64_t>>&);
template float Arch::get_sd_v_norm(const Global::Vector<LAFEM::DenseVectorBlocked<float, std::uint64_t, 2>, LAFEM::VectorMirror<float, std::uint64_t>>&);
#ifdef FEAT_HAVE_HALFMATH
template Half Arch::get_sd_v_norm(const Global::Vector<LAFEM::DenseVectorBlocked<Half, std::uint32_t, 2>, LAFEM::VectorMirror<Half, std::uint32_t>>&);
template Half Arch::get_sd_v_norm(const Global::Vector<LAFEM::DenseVectorBlocked<Half, std::uint64_t, 2>, LAFEM::VectorMirror<Half, std::uint64_t>>&);
#endif

template double Arch::get_sd_v_norm_host(const LAFEM::DenseVectorBlocked<double, std::uint32_t, 2>&);
template float Arch::get_sd_v_norm_host(const LAFEM::DenseVectorBlocked<float, std::uint32_t, 2>&);
template double Arch::get_sd_v_norm_host(const LAFEM::DenseVectorBlocked<double, std::uint64_t, 2>&);
template float Arch::get_sd_v_norm_host(const LAFEM::DenseVectorBlocked<float, std::uint64_t, 2>&);
#ifdef FEAT_HAVE_HALFMATH
template Half Arch::get_sd_v_norm_host(const LAFEM::DenseVectorBlocked<Half, std::uint32_t, 2>&);
template Half Arch::get_sd_v_norm_host(const LAFEM::DenseVectorBlocked<Half, std::uint64_t, 2>&);
#endif

template double Arch::get_sd_v_norm_host(const Global::Vector<LAFEM::DenseVectorBlocked<double, std::uint32_t, 2>, LAFEM::VectorMirror<double, std::uint32_t>>&);
template float Arch::get_sd_v_norm_host(const Global::Vector<LAFEM::DenseVectorBlocked<float, std::uint32_t, 2>, LAFEM::VectorMirror<float, std::uint32_t>>&);
template double Arch::get_sd_v_norm_host(const Global::Vector<LAFEM::DenseVectorBlocked<double, std::uint64_t, 2>, LAFEM::VectorMirror<double, std::uint64_t>>&);
template float Arch::get_sd_v_norm_host(const Global::Vector<LAFEM::DenseVectorBlocked<float, std::uint64_t, 2>, LAFEM::VectorMirror<float, std::uint64_t>>&);
#ifdef FEAT_HAVE_HALFMATH
template Half Arch::get_sd_v_norm_host(const Global::Vector<LAFEM::DenseVectorBlocked<Half, std::uint32_t, 2>, LAFEM::VectorMirror<Half, std::uint32_t>>&);
template Half Arch::get_sd_v_norm_host(const Global::Vector<LAFEM::DenseVectorBlocked<Half, std::uint64_t, 2>, LAFEM::VectorMirror<Half, std::uint64_t>>&);
#endif

/*********************************************************3D implementations**************************************************************************************/

template void Arch::assemble_burgers_csr(const Q2StandardHexa&, const CSRMatrixData<double, std::uint32_t>&, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double, const AssemblyBurgersData<double>&);
template void Arch::assemble_burgers_csr(const Q2StandardHexa&, const CSRMatrixData<float, std::uint32_t>&, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float, const AssemblyBurgersData<float>&);
template void Arch::assemble_burgers_csr(const Q2StandardHexa&, const CSRMatrixData<double, std::uint64_t>&, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double, const AssemblyBurgersData<double>&);
template void Arch::assemble_burgers_csr(const Q2StandardHexa&, const CSRMatrixData<float, std::uint64_t>&, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float, const AssemblyBurgersData<float>&);
#ifdef FEAT_HAVE_HALFMATH
template void Arch::assemble_burgers_csr(const Q2StandardHexa&, const CSRMatrixData<Half, std::uint32_t>&, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, Half, const AssemblyBurgersData<Half>&);
template void Arch::assemble_burgers_csr(const Q2StandardHexa&, const CSRMatrixData<Half, std::uint64_t>&, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, Half, const AssemblyBurgersData<Half>&);
#endif

template void Arch::assemble_burgers_csr_host(const Q2StandardHexa&, const CSRMatrixData<double, std::uint32_t>&, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double, const AssemblyBurgersData<double>&);
template void Arch::assemble_burgers_csr_host(const Q2StandardHexa&, const CSRMatrixData<float, std::uint32_t>&, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float, const AssemblyBurgersData<float>&);
template void Arch::assemble_burgers_csr_host(const Q2StandardHexa&, const CSRMatrixData<double, std::uint64_t>&, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double, const AssemblyBurgersData<double>&);
template void Arch::assemble_burgers_csr_host(const Q2StandardHexa&, const CSRMatrixData<float, std::uint64_t>&, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float, const AssemblyBurgersData<float>&);
#ifdef FEAT_HAVE_HALFMATH
template void Arch::assemble_burgers_csr_host(const Q2StandardHexa&, const CSRMatrixData<Half, std::uint32_t>&, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, Half, const AssemblyBurgersData<Half>&);
template void Arch::assemble_burgers_csr_host(const Q2StandardHexa&, const CSRMatrixData<Half, std::uint64_t>&, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, Half, const AssemblyBurgersData<Half>&);
#endif

template void Arch::assemble_burgers_defect(const Q2StandardHexa&, double*, const double*, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double, const AssemblyBurgersData<double>&);
template void Arch::assemble_burgers_defect(const Q2StandardHexa&, float*, const float*, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float, const AssemblyBurgersData<float>&);
template void Arch::assemble_burgers_defect(const Q2StandardHexa&, double*, const double*, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double, const AssemblyBurgersData<double>&);
template void Arch::assemble_burgers_defect(const Q2StandardHexa&, float*, const float*, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float, const AssemblyBurgersData<float>&);
#ifdef FEAT_HAVE_HALFMATH
template void Arch::assemble_burgers_defect(const Q2StandardHexa&, Half*, const Half*, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, Half, const AssemblyBurgersData<Half>&);
template void Arch::assemble_burgers_defect(const Q2StandardHexa&, Half*, const Half*, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, Half, const AssemblyBurgersData<Half>&);
#endif

template void Arch::assemble_burgers_defect_host(const Q2StandardHexa&, double*, const double*, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double, const AssemblyBurgersData<double>&);
template void Arch::assemble_burgers_defect_host(const Q2StandardHexa&, float*, const float*, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float, const AssemblyBurgersData<float>&);
template void Arch::assemble_burgers_defect_host(const Q2StandardHexa&, double*, const double*, const double*, const AssemblyCubatureData<double>&, const AssemblyMappingData<double, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, double, const AssemblyBurgersData<double>&);
template void Arch::assemble_burgers_defect_host(const Q2StandardHexa&, float*, const float*, const float*, const AssemblyCubatureData<float>&, const AssemblyMappingData<float, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, float, const AssemblyBurgersData<float>&);
#ifdef FEAT_HAVE_HALFMATH
template void Arch::assemble_burgers_defect_host(const Q2StandardHexa&, Half*, const Half*, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint32_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, Half, const AssemblyBurgersData<Half>&);
template void Arch::assemble_burgers_defect_host(const Q2StandardHexa&, Half*, const Half*, const Half*, const AssemblyCubatureData<Half>&, const AssemblyMappingData<Half, std::uint64_t>&,
                                          const std::vector<std::vector<int>>&, const std::vector<void*>&, Half, const AssemblyBurgersData<Half>&);
#endif

template double Arch::get_sd_v_norm(const LAFEM::DenseVectorBlocked<double, std::uint32_t, 3>&);
template float Arch::get_sd_v_norm(const LAFEM::DenseVectorBlocked<float, std::uint32_t, 3>&);
template double Arch::get_sd_v_norm(const LAFEM::DenseVectorBlocked<double, std::uint64_t, 3>&);
template float Arch::get_sd_v_norm(const LAFEM::DenseVectorBlocked<float, std::uint64_t, 3>&);
#ifdef FEAT_HAVE_HALFMATH
template Half Arch::get_sd_v_norm(const LAFEM::DenseVectorBlocked<Half, std::uint32_t, 3>&);
template Half Arch::get_sd_v_norm(const LAFEM::DenseVectorBlocked<Half, std::uint64_t, 3>&);
#endif

template double Arch::get_sd_v_norm(const Global::Vector<LAFEM::DenseVectorBlocked<double, std::uint32_t, 3>, LAFEM::VectorMirror<double, std::uint32_t>>&);
template float Arch::get_sd_v_norm(const Global::Vector<LAFEM::DenseVectorBlocked<float, std::uint32_t, 3>, LAFEM::VectorMirror<float, std::uint32_t>>&);
template double Arch::get_sd_v_norm(const Global::Vector<LAFEM::DenseVectorBlocked<double, std::uint64_t, 3>, LAFEM::VectorMirror<double, std::uint64_t>>&);
template float Arch::get_sd_v_norm(const Global::Vector<LAFEM::DenseVectorBlocked<float, std::uint64_t, 3>, LAFEM::VectorMirror<float, std::uint64_t>>&);
#ifdef FEAT_HAVE_HALFMATH
template Half Arch::get_sd_v_norm(const Global::Vector<LAFEM::DenseVectorBlocked<Half, std::uint32_t, 3>, LAFEM::VectorMirror<Half, std::uint32_t>>&);
template Half Arch::get_sd_v_norm(const Global::Vector<LAFEM::DenseVectorBlocked<Half, std::uint64_t, 3>, LAFEM::VectorMirror<Half, std::uint64_t>>&);
#endif

template double Arch::get_sd_v_norm_host(const LAFEM::DenseVectorBlocked<double, std::uint32_t, 3>&);
template float Arch::get_sd_v_norm_host(const LAFEM::DenseVectorBlocked<float, std::uint32_t, 3>&);
template double Arch::get_sd_v_norm_host(const LAFEM::DenseVectorBlocked<double, std::uint64_t, 3>&);
template float Arch::get_sd_v_norm_host(const LAFEM::DenseVectorBlocked<float, std::uint64_t, 3>&);
#ifdef FEAT_HAVE_HALFMATH
template Half Arch::get_sd_v_norm_host(const LAFEM::DenseVectorBlocked<Half, std::uint32_t, 3>&);
template Half Arch::get_sd_v_norm_host(const LAFEM::DenseVectorBlocked<Half, std::uint64_t, 3>&);
#endif

template double Arch::get_sd_v_norm_host(const Global::Vector<LAFEM::DenseVectorBlocked<double, std::uint32_t, 3>, LAFEM::VectorMirror<double, std::uint32_t>>&);
template float Arch::get_sd_v_norm_host(const Global::Vector<LAFEM::DenseVectorBlocked<float, std::uint32_t, 3>, LAFEM::VectorMirror<float, std::uint32_t>>&);
template double Arch::get_sd_v_norm_host(const Global::Vector<LAFEM::DenseVectorBlocked<double, std::uint64_t, 3>, LAFEM::VectorMirror<double, std::uint64_t>>&);
template float Arch::get_sd_v_norm_host(const Global::Vector<LAFEM::DenseVectorBlocked<float, std::uint64_t, 3>, LAFEM::VectorMirror<float, std::uint64_t>>&);
#ifdef FEAT_HAVE_HALFMATH
template Half Arch::get_sd_v_norm_host(const Global::Vector<LAFEM::DenseVectorBlocked<Half, std::uint32_t, 3>, LAFEM::VectorMirror<Half, std::uint32_t>>&);
template Half Arch::get_sd_v_norm_host(const Global::Vector<LAFEM::DenseVectorBlocked<Half, std::uint64_t, 3>, LAFEM::VectorMirror<Half, std::uint64_t>>&);
#endif