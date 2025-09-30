#include <kernel/voxel_assembly/burgers_velo_material_assembler.hpp>
#include <kernel/lafem/matrix_gather_scatter_helper.hpp>
#include <kernel/lafem/vector_gather_scatter_helper.hpp>
#include <kernel/util/math.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/global/vector.hpp>

namespace FEAT
{
  namespace VoxelAssembly
  {
    namespace Kernel
    {

      template<typename Space_, typename DT_, typename IT_, typename ViscFunc, typename ViscDFunc,
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
                ViscFunc visco_func, ViscDFunc visco_d_func)
      {
        //define types
        typedef Space_ SpaceType;
        typedef DT_ DataType;
        typedef IT_ IndexType;

        const DataType& beta{burgers_params.beta};
        const DataType& sd_delta{burgers_params.sd_delta};
        const DataType& sd_v_norm{burgers_params.sd_v_norm};

        const DataType tol_eps = Math::sqrt(Math::eps<DataType>());
        // const DataType tol_eps = CudaMath::cuda_sqrt(DBL_EPSILON);
        // const DataType tol_eps = DBL_EPSILON;
        const bool need_streamline = (Math::abs(sd_delta) > DataType(0)) && (sd_v_norm > tol_eps);
        const bool need_convection = Math::abs(beta) > DataType(0);

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

        FEAT_PRAGMA_OMP(parallel for)
        for(Index idx = 0; idx < coloring_size; ++idx)
        {
          // define local coefficients
          Tiny::Matrix<DataType, dim, num_loc_verts> local_coeffs;
          // and local matrix
          LocalMatrixType loc_mat;
          LocalVectorType local_conv_dofs(DataType(0));
          // now do work for this cell
          IndexType cell = IndexType(coloring_map[idx]);
          // std::cout << "Starting with cell " << cell << std::endl;
          const IndexSetWrapper<IndexType> local_dofs_w{cell_to_dof, IndexType(num_loc_dofs)};
          const IndexType* local_dofs = cell_to_dof + cell*num_loc_dofs;
          const IndexType* local_dof_sorter = cell_to_dof_sorter + cell*num_loc_dofs;

          SpaceHelp::set_coefficients(local_coeffs, local_dofs_w, nodes, cell);
          //always gather local conv data
          LAFEM::template VectorGatherScatterHelper<SpaceType, DataType, IndexType>::gather_vector_dense(local_conv_dofs,
                    (const VecValueType*)conv_data, IndexType(matrix_num_rows), local_dofs,DataType(1));
          if(need_streamline)
          {
            VoxelAssembly::Kernel::burgers_velo_material_mat_assembly_kernel<SpaceHelp, LocalMatrixType, LocalVectorType, dim, num_loc_verts, ViscFunc, ViscDFunc, true>(
                            loc_mat, local_conv_dofs, local_coeffs, cub_pt, cub_wg, num_cubs, burgers_params, material_params, need_convection, tol_eps, visco_func, visco_d_func);
          }
          else
          {
            VoxelAssembly::Kernel::burgers_velo_material_mat_assembly_kernel<SpaceHelp, LocalMatrixType, LocalVectorType, dim, num_loc_verts, ViscFunc, ViscDFunc, false>(
                            loc_mat, local_conv_dofs, local_coeffs, cub_pt, cub_wg, num_cubs, burgers_params, material_params, need_convection, tol_eps, visco_func, visco_d_func);
          }

          //scatter
          LAFEM::template MatrixGatherScatterHelper<SpaceType, DataType, IndexType, pol_>::scatter_matrix_csr(loc_mat, (MatValueType*)matrix_data, local_dofs, local_dofs, IndexType(matrix_num_rows), IndexType(matrix_num_cols), matrix_row_ptr, matrix_col_idx, alpha, local_dof_sorter);

        }
      }

      template<typename Space_, typename DT_, typename IT_, typename ViscFunc>
      void full_burgers_vm_assembler_vector_bd_host(DT_* vector_data,
                const DT_* conv_data, const DT_* primal_data, Index vec_size,
                const Tiny::Vector<DT_, Space_::world_dim>* cub_pt,
                const DT_*  cub_wg, int num_cubs, DT_ alpha,
                const IT_* cell_to_dof, [[maybe_unused]] Index cell_num,
                const Tiny::Vector<DT_, Space_::world_dim>* nodes, [[maybe_unused]] Index node_size,
                const int* coloring_map, Index coloring_size,
                const VoxelAssembly::AssemblyBurgersData<DT_>& burgers_params,
                const ViscFunc& visco_func)
      {
        //define types
        typedef Space_ SpaceType;
        typedef DT_ DataType;
        typedef IT_ IndexType;

        const DataType& beta{burgers_params.beta};

        const bool need_convection = Math::abs(beta) > DataType(0);

        typedef VoxelAssembly::SpaceHelper<SpaceType, DT_, IT_> SpaceHelp;

        constexpr int dim = SpaceHelp::dim;

        //define local sizes
        constexpr int num_loc_dofs = SpaceType::DofMappingType::dof_count;
        //get number of nodes per element
        constexpr int num_loc_verts = SpaceType::MeshType::template IndexSet<dim, 0>::Type::num_indices;

        //our local datatypes
        typedef Tiny::Vector<DataType, dim> VecValueType;
        // typedef Tiny::Matrix<DataType, dim, dim> MatValueType;
        typedef Tiny::Vector<VecValueType, num_loc_dofs> LocalVectorType;
        // typedef Tiny::Matrix<MatValueType, num_loc_dofs, num_loc_dofs> LocalMatrixType;

        FEAT_PRAGMA_OMP(parallel for)
        for(Index idx = 0; idx < coloring_size; ++idx)
        {
          //define local array
          LocalVectorType loc_vec(DataType(0));
          LocalVectorType local_conv_dofs(DataType(0));
          LocalVectorType local_prim_dofs(DataType(0));
          // typename NewSpaceHelp::ImagePointType img_point;
          Tiny::Matrix<DataType, dim, num_loc_verts> local_coeffs;

          //now do work for this cell
          IndexType cell = IndexType(coloring_map[idx]);
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
      }
    }

    namespace Arch
    {
      template<typename Space_, typename DT_, typename IT_>
      void assemble_burgers_velo_material_csr_host([[maybe_unused]]const Space_& space,
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
              VoxelAssembly::Kernel::template full_burgers_vm_assembler_matrix1_bcsr_host<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreau>, Intern::ViscoDFunctor<DT_, MaterialType::carreau>,
                                                                                          FEAT::Intern::MatrixGatherScatterPolicy::useLocalSortHelper>(
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
              VoxelAssembly::Kernel::template full_burgers_vm_assembler_matrix1_bcsr_host<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreauYasuda>, Intern::ViscoDFunctor<DT_, MaterialType::carreauYasuda>,
                                                                                          FEAT::Intern::MatrixGatherScatterPolicy::useLocalSortHelper>(
                  matrix_data.data, conv_data, matrix_data.row_ptr, matrix_data.col_idx, matrix_data.num_rows, matrix_data.num_cols,
                  (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
                  cubature.cub_wg, cubature.num_cubs, alpha,
                  dof_mapping.cell_to_dof, dof_mapping.cell_num,
                  (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
                  (const int*) coloring_maps[col], coloring_map_sizes[col], dof_mapping.cell_to_dof_sorter,
                  burgers_params, material_params, Intern::ViscoFunctor<DT_, MaterialType::carreauYasuda>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T,material_params.yasuda_a, material_params.mu_inf},
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
        switch(material_type)
        {
          case MaterialType::carreau:
          {
            for(Index col = 0; col < Index(coloring_maps.size()); ++col)
            {
              VoxelAssembly::Kernel::template full_burgers_vm_assembler_vector_bd_host<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreau>>(vector_data, conv_data, primal_data,
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
              VoxelAssembly::Kernel::template full_burgers_vm_assembler_vector_bd_host<Space_, DT_, IT_, Intern::ViscoFunctor<DT_, MaterialType::carreauYasuda>>(vector_data, conv_data, primal_data,
                  space.get_num_dofs(), (const typename Tiny::Vector<DT_, Space_::world_dim>*) cubature.cub_pt,
                  cubature.cub_wg, cubature.num_cubs, alpha,
                  dof_mapping.cell_to_dof, dof_mapping.cell_num,
                  (const typename Tiny::Vector<DT_, Space_::world_dim>*) dof_mapping.nodes, dof_mapping.node_size,
                  (const int*) coloring_maps[col], coloring_map_sizes[col], burgers_params,
                  Intern::ViscoFunctor<DT_, MaterialType::carreauYasuda>{material_params.mu_0, material_params.exp, material_params.lambda, material_params.a_T, material_params.yasuda_a, material_params.mu_inf}
              );
            }
            break;
          }
        }
      }

    }
  }
}

//instantiate the templates
using namespace FEAT;
using namespace FEAT::VoxelAssembly;


/*******************************************************2D implementations***************************************************/
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
