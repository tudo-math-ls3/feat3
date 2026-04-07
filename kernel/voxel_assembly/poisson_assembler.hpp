// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>
#include <kernel/voxel_assembly/voxel_assembly_common.hpp>
#include <kernel/voxel_assembly/helper/data_handler.hpp>
#include <kernel/voxel_assembly/helper/space_helper.hpp>
#include <kernel/cubature/rule.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/util/tiny_algebra.hpp>

#ifdef FEAT_HAVE_CUDA
#include <kernel/util/cuda_util.hpp>
#endif

namespace FEAT
{
  namespace VoxelAssembly
  {

    /**
     * \brief Poisson Voxel Assembly template.
     *
     * Has to be specialized for each specific FE element space.
     *
     * \tparam Space_ The FE space for which it should be assembled
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indextype to be used.
     */
    template<typename Space_, typename DT_, typename IT_>
    class VoxelPoissonAssembler DOXY({});


    /**
     * \brief Namespace for the actual assembly kernels.
     */
    namespace Kernel
    {
      template<typename SpaceHelp_, typename LocMatType_, int dim_, int num_verts_>
      CUDA_HOST_DEVICE void poisson_assembly_kernel(LocMatType_& loc_mat, const Tiny::Matrix<typename SpaceHelp_::DataType, dim_, num_verts_>& local_coeffs,
                                            const typename SpaceHelp_::DomainPointType* cub_pt, const typename SpaceHelp_::DataType* cub_wg, const int num_cubs)
      {
        typedef SpaceHelp_ SpaceHelp;
        // constexpr int dim = SpaceHelp::dim;
        typedef typename SpaceHelp::SpaceType SpaceType;
        typedef typename SpaceHelp::DataType DataType;

        //define local sizes
        constexpr int num_loc_dofs = SpaceType::DofMappingType::dof_count;
        //get number of nodes per element
        // constexpr int num_loc_verts = SpaceType::MeshType::template IndexSet<dim, 0>::Type::num_indices;
        // define local arrays
        typename SpaceHelp::JacobianMatrixType loc_jac, loc_jac_inv;
        typename SpaceHelp::DomainPointType dom_point;
        typename SpaceHelp::EvalData basis_data;
        DataType jac_det = DataType(0);

        loc_mat.format();
        for(int cub_ind = 0; cub_ind < num_cubs; ++cub_ind)
        {
          dom_point = cub_pt[cub_ind];
          // NewSpaceHelp::map_point(img_point, dom_point, local_coeffs);
          SpaceHelp::calc_jac_mat(loc_jac, dom_point, local_coeffs);
          loc_jac_inv.set_inverse(loc_jac);
          jac_det = loc_jac.det();

          SpaceHelp::eval_ref_gradients(basis_data, dom_point);
          SpaceHelp::trans_gradients(basis_data, loc_jac_inv);

          for(int i = 0; i < num_loc_dofs; ++i)
          {
            for(int j = 0; j < num_loc_dofs; ++j)
            {
              //simple assembly for poisson operator
              Tiny::axpy(loc_mat(i,j), Tiny::dot(basis_data.phi[j].grad, basis_data.phi[i].grad), jac_det * cub_wg[cub_ind]);
            }
          }

        }
      }
    }

    /**
     * \brief Namespace for the kernel wrapper.
     */
    namespace Arch
    {
      #ifdef FEAT_HAVE_CUDA
      template<typename Space_, typename DT_, typename IT_>
      void assemble_poisson_csr(const Space_& space,
              const CSRMatrixData<DT_, IT_>& matrix_data,
              const AssemblyCubatureData<DT_>& cubature,
              const AssemblyMappingData<DT_, IT_>& dof_mapping,
              const Adjacency::ColoringDataHandler& coloring_data,
              DT_ alpha
              );
      #endif

      template<typename Space_, typename DT_, typename IT_>
      void assemble_poisson_csr_host(const Space_& space,
              const CSRMatrixData<DT_, IT_>& matrix_data,
              const AssemblyCubatureData<DT_>& cubature,
              const AssemblyMappingData<DT_, IT_>& dof_mapping,
              const Adjacency::ColoringDataHandler& coloring_data,
              DT_ alpha
              );

    }

#ifndef __CUDACC__
    //specialize for Q2Lagrange with standard Trafo
    template<int dim_, typename DT_, typename IT_>
    class VoxelPoissonAssembler<Q2StandardHyperCube<dim_>, DT_, IT_>
    {
    public:
      /// typedefs
      typedef Q2StandardHyperCube<dim_> SpaceType;
      typedef LagrangeDataHandler<SpaceType, DT_, IT_> DataHandler;
      typedef SpaceHelper<SpaceType, DT_, IT_> SpaceHelp;
      typedef typename SpaceHelp::ShapeType ShapeType;
      typedef typename SpaceHelp::DataType DataType;
      typedef typename SpaceHelp::IndexType IndexType;

      /// constexpr
      static constexpr int dim = SpaceHelp::dim;

      typedef typename SpaceHelp::DomainPointType DomainPointType;
      typedef typename SpaceHelp::ImagePointType ImagePointType;
      typedef typename SpaceHelp::ValueType ValueType;
      typedef typename SpaceHelp::JacobianMatrixType JacobianMatrixType;

      /// The datahandler
      DataHandler mesh_data;


    public:
      explicit VoxelPoissonAssembler() = default;

      template<typename ColoringType_>
      explicit VoxelPoissonAssembler(const SpaceType& space, const ColoringType_& coloring, int hint = -1) :
      mesh_data(space, coloring, hint)
      {
        #ifdef FEAT_HAVE_CUDA
        #ifdef DEBUG
        const std::size_t stack_limit = Util::cuda_get_max_cache_thread();
        const std::size_t stack_limit_target = sizeof(DataType) * (dim == 3 ? 8096u : 1012u);
        if(stack_limit < stack_limit_target)
          Util::cuda_set_max_cache_thread(stack_limit_target);
        #endif
        #endif
      }

      // rule of 5
      VoxelPoissonAssembler(const VoxelPoissonAssembler&) = delete;

      VoxelPoissonAssembler& operator=(const VoxelPoissonAssembler&) = delete;

      VoxelPoissonAssembler(VoxelPoissonAssembler&&) = default;

      VoxelPoissonAssembler& operator=(VoxelPoissonAssembler&&) = delete;

      ~VoxelPoissonAssembler(){}

      template<typename CubatureFactory_>
      void assemble_matrix1(LAFEM::SparseMatrixCSR<DataType, IndexType>& matrix, const SpaceType& space, const CubatureFactory_& cubature_factory, DataType alpha = DataType(1)) const
      {
        XASSERTM(matrix.num_rows() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.num_cols() == space.get_num_dofs(), "invalid matrix dimensions");

        //define cubature
        typedef Cubature::Rule<ShapeType, DataType, DataType> CubatureRuleType;
        CubatureRuleType cubature(Cubature::ctor_factory, cubature_factory);

        //get cubature points and weights
        int num_cubs = cubature.get_num_points();
        typename CubatureRuleType::PointType* cub_pt = cubature.get_points();
        DataType* cub_wg = cubature.get_weights();

        if(Backend::get_preferred_backend() == PreferredBackend::cuda)
          assemble_matrix1_cuda(matrix, space, cub_pt, cub_wg, num_cubs, alpha);
        else
          assemble_matrix1_generic(matrix, space, cub_pt, cub_wg, num_cubs, alpha);
      }


      void assemble_matrix1_generic(LAFEM::SparseMatrixCSR<DataType, IndexType>& matrix, const SpaceType& space, typename Cubature::Rule<ShapeType, DataType, DataType>::PointType* cub_pt, const DataType* cub_wg, int num_cubs, DataType alpha) const
      {
        Memory::TypedView<DataType> val_view(matrix.val_view_rw());
        const Memory::TypedView<IndexType> row_ptr_view(matrix.row_ptr_view_r());
        const Memory::TypedView<IndexType> col_idx_view(matrix.col_idx_view_r());
        CSRMatrixData<DataType, IndexType> mat_data = {val_view.get_w(), row_ptr_view.get_r(), col_idx_view.get_r(), matrix.num_rows(), matrix.num_cols()};

        AssemblyCubatureData<DataType> cub_data = {(void*)cub_pt, cub_wg, num_cubs};
        AssemblyMappingData<DataType, IndexType> mapping_data = mesh_data.get_assembly_field();


        VoxelAssembly::Arch::assemble_poisson_csr_host(space, mat_data, cub_data, mapping_data, mesh_data.coloring_data, alpha);
      }

      #ifdef FEAT_HAVE_CUDA
      void assemble_matrix1_cuda(LAFEM::SparseMatrixCSR<DataType, IndexType>& matrix, const SpaceType& space, typename Cubature::Rule<ShapeType, DataType, DataType>::PointType* cub_pt, const DataType* cub_wg, int num_cubs, DataType alpha) const
      {
        Memory::TypedView<DataType> val_view(matrix.val_view_rw(Memory::Location::cuda));
        const Memory::TypedView<IndexType> row_ptr_view(matrix.row_ptr_view_r(Memory::Location::cuda));
        const Memory::TypedView<IndexType> col_idx_view(matrix.col_idx_view_r(Memory::Location::cuda));
        CSRMatrixData<DataType, IndexType> mat_data = {val_view.get_w(), row_ptr_view.get_r(), col_idx_view.get_r(), matrix.num_rows(), matrix.num_cols()};

        typedef typename Cubature::Rule<ShapeType, DataType, DataType>::PointType CubPointType;
        //initialize all necessary pointer arrays and values //maybe more sense to specify cubature rule and set this to a const mem location?
        void* cub_pt_device = Memory::alloc_cuda(std::size_t(num_cubs) * sizeof(CubPointType));
        Memory::memcopy_main_to_cuda(cub_pt_device, (const void*)cub_pt, std::size_t(num_cubs) * sizeof(CubPointType));

        void* cub_wg_device = Memory::alloc_cuda(std::size_t(num_cubs) * sizeof(DataType));
        Memory::memcopy_main_to_cuda(cub_wg_device, (const void*)cub_wg, std::size_t(num_cubs) * sizeof(DataType));

        VoxelAssembly::AssemblyCubatureData<DataType> d_cub_data = {cub_pt_device, (DataType*)cub_wg_device, num_cubs};
        VoxelAssembly::AssemblyMappingData<DataType, IndexType> d_mapping_data = mesh_data.get_assembly_field();


        VoxelAssembly::Arch::assemble_poisson_csr(space, mat_data, d_cub_data, d_mapping_data, mesh_data.coloring_data, alpha);
        //free resources
        Memory::free_cuda(cub_wg_device);
        Memory::free_cuda(cub_pt_device);
      }

      #else //FEAT_HAVE_CUDA
      void assemble_matrix1_cuda(LAFEM::SparseMatrixCSR<DataType, IndexType>& DOXY(matrix), const SpaceType& DOXY(space), typename Cubature::Rule<ShapeType, DataType, DataType>::PointType* DOXY(cub_pt), const DataType* DOXY(cub_wg), int DOXY(num_cubs), DataType DOXY(alpha)) const
      {
        XABORTM("What in the nine hells are you doing here?");
      }
      #endif

    }; //class GPUPoissonAssembler

    #endif // __CUDACC__
  }
}
