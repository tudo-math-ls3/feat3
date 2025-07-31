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
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
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
     * \brief Deformation Voxel Assembly template.
     *
     * Has to be specialized for each specific FE element space.
     *
     * \tparam Space_ The FE space for which it should be assembled
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indextype to be used.
     */
    template<typename Space_, typename DT_, typename IT_>
    class VoxelDefoAssembler DOXY({});

    namespace Kernel
    {
      template<typename SpaceHelp_, typename LocMatType_, int dim_, int num_verts_>
      CUDA_HOST_DEVICE void defo_assembly_kernel(LocMatType_& loc_mat, const Tiny::Matrix<typename SpaceHelp_::DataType, dim_, num_verts_>& local_coeffs,
                                          const typename SpaceHelp_::DomainPointType* cub_pt, const typename SpaceHelp_::DataType* cub_wg, const int num_cubs,
                                          const typename SpaceHelp_::DataType nu)
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

        loc_mat.format();
        for(int cub_ind = 0; cub_ind < num_cubs; ++cub_ind)
        {
          dom_point = cub_pt[cub_ind];
          // NewSpaceHelp::map_point(img_point, dom_point, local_coeffs);
          SpaceHelp::calc_jac_mat(loc_jac, dom_point, local_coeffs);
          loc_jac_inv.set_inverse(loc_jac);

          SpaceHelp::eval_ref_gradients(basis_data, dom_point);
          SpaceHelp::trans_gradients(basis_data, loc_jac_inv);

          const DataType weight = loc_jac.det() * cub_wg[cub_ind];
          for(int i = 0; i < num_loc_dofs; ++i)
          {
            for(int j = 0; j < num_loc_dofs; ++j)
            {
              //now assemble outer and inner products
              loc_mat[i][j].add_scalar_main_diag(nu * weight * Tiny::dot(basis_data.phi[i].grad, basis_data.phi[j].grad));
              loc_mat[i][j].add_outer_product(basis_data.phi[j].grad, basis_data.phi[i].grad, nu * weight);
            }
          }

        }

      }
    }

    namespace Arch
    {
      #ifdef FEAT_HAVE_CUDA
      template<typename Space_, typename DT_, typename IT_>
      void assemble_defo_csr(const Space_& space,
              const CSRMatrixData<DT_, IT_>& matrix_data,
              const AssemblyCubatureData<DT_>& cubature,
              const AssemblyMappingData<DT_, IT_>& dof_mapping,
              const std::vector<int*>& coloring_maps,
              const std::vector<Index>& coloring_map_sizes,
              DT_ alpha, DT_ nu
              );
      #endif

      template<typename Space_, typename DT_, typename IT_>
      void assemble_defo_csr_host(const Space_& space,
              const CSRMatrixData<DT_, IT_>& matrix_data,
              const AssemblyCubatureData<DT_>& cubature,
              const AssemblyMappingData<DT_, IT_>& dof_mapping,
              const std::vector<int*>& coloring_maps,
              const std::vector<Index>& coloring_map_sizes,
              DT_ alpha, DT_ nu
              );
    }

#ifndef __CUDACC__
    //specialize for Q2Lagrange with standard Trafo
    template<int dim_, typename DT_, typename IT_>
    class VoxelDefoAssembler<Q2StandardHyperCube<dim_>, DT_, IT_>
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

      /// Scaling parameter
      DataType nu;


    public:
      explicit VoxelDefoAssembler() = default;

      template<typename ColoringType_>
      explicit VoxelDefoAssembler(const SpaceType& space, const ColoringType_& coloring, int hint = -1) :
      mesh_data(space, coloring, hint), nu(DataType(0.))
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
      VoxelDefoAssembler(const VoxelDefoAssembler&) = delete;

      VoxelDefoAssembler& operator=(const VoxelDefoAssembler&) = delete;

      VoxelDefoAssembler(VoxelDefoAssembler&&) = default;

      VoxelDefoAssembler& operator=(VoxelDefoAssembler&&) = delete;

      ~VoxelDefoAssembler(){}


      template<typename CubatureFactory_>
      void assemble_matrix1(LAFEM::SparseMatrixBCSR<DataType, IndexType, dim, dim>& matrix, const SpaceType& space, const CubatureFactory_& cubature_factory, DataType alpha = DataType(1)) const
      {
        XASSERTM(matrix.rows() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.columns() == space.get_num_dofs(), "invalid matrix dimensions");

        //define cubature
        typedef Cubature::Rule<ShapeType, DataType, DataType> CubatureRuleType;
        CubatureRuleType cubature(Cubature::ctor_factory, cubature_factory);

        //get cubature points and weights
        int num_cubs = cubature.get_num_points();
        typename CubatureRuleType::PointType* cub_pt = cubature.get_points();
        DataType* cub_wg = cubature.get_weights();

        BACKEND_SKELETON_VOID(assemble_matrix1_cuda, assemble_matrix1_generic, assemble_matrix1_generic, matrix, space, cub_pt, cub_wg, num_cubs, alpha)
      }


      void assemble_matrix1_generic(LAFEM::SparseMatrixBCSR<DataType, IndexType, dim, dim>& matrix, const SpaceType& space, typename Cubature::Rule<ShapeType, DataType, DataType>::PointType* cub_pt, const DataType* cub_wg, int num_cubs, DataType alpha) const
      {
        CSRMatrixData<DataType, IndexType> mat_data = {matrix.template val<LAFEM::Perspective::pod>(), matrix.row_ptr(), matrix.col_ind(), matrix.rows(), matrix.columns()};

        AssemblyCubatureData<DataType> cub_data = {(void*)cub_pt, cub_wg, num_cubs};
        AssemblyMappingData<DataType, IndexType> mapping_data = mesh_data.get_assembly_field();


        VoxelAssembly::Arch::assemble_defo_csr_host(space, mat_data, cub_data, mapping_data, mesh_data.get_coloring_maps(), mesh_data.get_color_map_sizes(), alpha, nu);
      }

      #ifdef FEAT_HAVE_CUDA
      void assemble_matrix1_cuda(LAFEM::SparseMatrixBCSR<DataType, IndexType, dim, dim>& matrix, const SpaceType& space, typename Cubature::Rule<ShapeType, DataType, DataType>::PointType* cub_pt, const DataType* cub_wg, int num_cubs, DataType alpha) const
      {
        VoxelAssembly::CSRMatrixData<DataType, IndexType> mat_data = {matrix.template val<LAFEM::Perspective::pod>(), matrix.row_ptr(), matrix.col_ind(), matrix.rows(), matrix.columns()};

        typedef typename Cubature::Rule<ShapeType, DataType, DataType>::PointType CubPointType;
        //initialize all necessary pointer arrays and values //maybe more sense to specify cubature rule and set this to a const mem location?
        void* cub_pt_device = Util::cuda_malloc(num_cubs * sizeof(CubPointType));
        Util::cuda_copy_host_to_device(cub_pt_device, (void*)cub_pt, num_cubs * sizeof(CubPointType));

        void* cub_wg_device = Util::cuda_malloc(num_cubs * sizeof(DataType));
        Util::cuda_copy_host_to_device(cub_wg_device, (void*)cub_wg, num_cubs * sizeof(DataType));

        VoxelAssembly::AssemblyCubatureData<DataType> d_cub_data = {cub_pt_device, (DataType*)cub_wg_device, num_cubs};
        VoxelAssembly::AssemblyMappingData<DataType, IndexType> d_mapping_data = mesh_data.get_assembly_field();


        VoxelAssembly::Arch::assemble_defo_csr(space, mat_data, d_cub_data, d_mapping_data,  mesh_data.get_coloring_maps(), mesh_data.get_color_map_sizes(), alpha, nu);
        //free resources
        Util::cuda_free(cub_wg_device);
        Util::cuda_free(cub_pt_device);
      }

      #else //FEAT_HAVE_CUDA
      void assemble_matrix1_cuda(LAFEM::SparseMatrixBCSR<DataType, IndexType, dim, dim>& DOXY(matrix), const SpaceType& DOXY(space), typename Cubature::Rule<ShapeType, DataType, DataType>::PointType* DOXY(cub_pt), const DataType* DOXY(cub_wg), int DOXY(num_cubs), DataType DOXY(alpha)) const
      {
        XABORTM("What in the nine hells are you doing here?");
      }
      #endif

    }; //class GPUPoissonAssembler

    #endif // __CUDACC__


  }
}
