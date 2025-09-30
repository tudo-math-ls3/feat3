#pragma once
#ifndef FEAT_KERNEL_VOXEL_ASSEMBLY_BURGERS_VELO_MATERIAL_HPP
#define FEAT_KERNEL_VOXEL_ASSEMBLY_BURGERS_VELO_MATERIAL_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/voxel_assembly/voxel_assembly_common.hpp>
#include <kernel/voxel_assembly/helper/data_handler.hpp>
#include <kernel/voxel_assembly/helper/space_helper.hpp>
#include <kernel/cubature/rule.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/voxel_assembly/burgers_assembler.hpp>
#ifdef FEAT_HAVE_CUDA
#include <kernel/util/cuda_util.hpp>
#endif

namespace FEAT
{
  namespace VoxelAssembly
  {

    #if defined(FEAT_HAVE_CUDA) || defined(DOXYGEN)
    /// Helper struct wrapping the data required as shared data
    template<typename SpaceHelp_, bool need_streamdiff_ = true>
    struct BurgersMatSharedDataKernelWrapper
    {
      typedef SpaceHelp_ SpaceHelp;
      static constexpr int dim = SpaceHelp::dim;
      typedef typename SpaceHelp::SpaceType SpaceType;
      typedef typename SpaceHelp::DataType DataType;
      //define local sizes
      static constexpr int num_loc_dofs = SpaceType::DofMappingType::dof_count;
      // local vector and matrix defines
      typedef Tiny::Vector<DataType, dim> VecValueType;
      typedef Tiny::Matrix<DataType, dim, dim> MatValueType;

      typename SpaceHelp::EvalData basis_data;
      typename SpaceHelp::JacobianMatrixType loc_jac;
      typename SpaceHelp::JacobianMatrixType loc_jac_inv;
      MatValueType loc_grad_v;
      VecValueType loc_v;
      VecValueType mean_v;
      typename SpaceHelp::DomainPointType dom_point;
      Tiny::Vector<DataType, num_loc_dofs> streamdiff_coeffs;
      DataType local_delta;
      DataType nu_loc;
      DataType gamma_dot;
      DataType det;
      DataType weight;
      bool need_frechet;
    };

    template<typename SpaceHelp_>
    struct BurgersMatSharedDataKernelWrapper<SpaceHelp_, false>
    {
      typedef SpaceHelp_ SpaceHelp;
      static constexpr int dim = SpaceHelp::dim;
      typedef typename SpaceHelp::SpaceType SpaceType;
      typedef typename SpaceHelp::DataType DataType;
      //define local sizes
      static constexpr int num_loc_dofs = SpaceType::DofMappingType::dof_count;
      // local vector and matrix defines
      typedef Tiny::Vector<DataType, dim> VecValueType;
      typedef Tiny::Matrix<DataType, dim, dim> MatValueType;

      typename SpaceHelp::EvalData basis_data;
      typename SpaceHelp::JacobianMatrixType loc_jac;
      typename SpaceHelp::JacobianMatrixType loc_jac_inv;
      MatValueType loc_grad_v;
      VecValueType loc_v;
      typename SpaceHelp::DomainPointType dom_point;
      DataType local_delta;
      DataType nu_loc;
      DataType gamma_dot;
      DataType det;
      DataType weight;
      bool need_frechet;
    };

    #endif
    /**
     * \brief Burgers Velocity Material Voxel Assembly template.
     *
     * Has to be specialized for each specific FE element space.
     *
     * \tparam Space_ The FE space for which it should be assembled
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indextype to be used.
     */
    template<typename Space_, typename DT_, typename IT_>
    class VoxelBurgersVeloMaterialAssembler DOXY({});

    namespace Kernel
    {

      /**
       * \brief Burgers Velo Material Matrix assembly kernel
       *
       * Outsources the the most inner execution of the assembly of the full burgers term for a cell previously gathered.
       * Local always refers to the local dofs of the current cell.
       *
       * \tparam SpaceHelp_ The psacehelper to be used. Holds most information about the underlying space.
       * \tparam LocMatType_ The Matrixtype of the local matrix dofs.
       * \tparam LocVecType_ The Vectortype of the local convection dofs.
       * \tparam dim_ Dimension of our space. For ease of access.
       * \tparam num_verts_ The number of vertices of one cell.
       * \tparam ViscFunc_ Function mapping from R to R.
       * \tparam ViscDerFunc_ Function mapping from R to R.
       * \tparam need_stream_diff_ Streamdiff required?
       *
       * \param[out] loc_mat A reference to the local matrix to be assembled.
       * \param[in] local_conv_dofs A reference to the local convection dofs. If convection is not needed values are ignored.
       * \param[in] local_coeffs The local trafo coefficients.
       * \param[in] cub_pt A pointer to the cubature point coordinate vector.
       * \param[in] cub_wg A pointer to the cubature weights. Same order as cub_pt array.
       * \param[in] burgers_params A struct holding the burgers parameter configuration.
       * \param[in] material_params A struct holding the additional material parameters.
       * \param[in] need_convection Do we need convection?
       * \param[in] tol_eps Tolerance for local stream norm to be regarded as zero.
       * \param[in] visco_func Function mapping invariants of velocity jacobian to viscosity.
       * \param[in] visco_d_func Function representing the derivative of visc_func.
       */
      template<typename SpaceHelp_, typename LocMatType_, typename LocVecType_, int dim_, int num_verts_, typename ViscFunc_, typename ViscDerFunc_, bool need_stream_diff_>
      CUDA_HOST_DEVICE void burgers_velo_material_mat_assembly_kernel(LocMatType_& loc_mat, const LocVecType_& local_conv_dofs, const Tiny::Matrix<typename SpaceHelp_::DataType, dim_, num_verts_>& local_coeffs,
                                          const typename SpaceHelp_::DomainPointType* cub_pt, const typename SpaceHelp_::DataType* cub_wg, const int num_cubs,
                                          const VoxelAssembly::AssemblyBurgersData<typename SpaceHelp_::DataType>& burgers_params,
                                          const VoxelAssembly::AssemblyMaterialData<typename SpaceHelp_::DataType>& material_params,
                                          const bool need_convection, const typename SpaceHelp_::DataType tol_eps,
                                          ViscFunc_ visco_func, ViscDerFunc_ visco_d_func)
      {
        typedef SpaceHelp_ SpaceHelp;
        constexpr int dim = SpaceHelp::dim;
        typedef typename SpaceHelp::SpaceType SpaceType;
        typedef typename SpaceHelp::DataType DataType;

        constexpr bool need_streamline = need_stream_diff_;

        // burgers params
        // const DataType& nu{burgers_params.nu};
        const DataType& theta{burgers_params.theta};
        const DataType& beta{burgers_params.beta};
        const DataType& frechet_beta{burgers_params.frechet_beta};
        const DataType& sd_delta{burgers_params.sd_delta};
        const DataType& sd_nu{burgers_params.sd_nu};
        const DataType& sd_v_norm{burgers_params.sd_v_norm};
        // formulation only makes sense with deformation formulation
        // const bool& deformation{burgers_params.deformation};

        // additional material parameters
        const DataType& frechet_material{material_params.frechet_material};
        const DataType& reg_eps{material_params.reg_eps};
        const bool& need_frechet_material{material_params.need_frechet_material};

        // #ifdef __CUDACC__
        // const DataType tol_eps = CudaMath::cuda_get_sqrt_eps<DataType>();
        // #else
        // const DataType tol_eps = Math::sqrt(Math::eps<DataType>());
        // #endif

        // #ifdef __CUDACC__
        // const bool need_streamline = (CudaMath::cuda_abs(sd_delta) > DataType(0)) && (sd_v_norm > tol_eps);
        // const bool need_convection = CudaMath::cuda_abs(beta) > DataType(0);
        // #else
        // const bool need_streamline = (Math::abs(sd_delta) > DataType(0)) && (sd_v_norm > tol_eps);
        // const bool need_convection = Math::abs(beta) > DataType(0);
        // #endif


        //define local sizes
        constexpr int num_loc_dofs = SpaceType::DofMappingType::dof_count;
        //get number of nodes per element
        // constexpr int num_loc_verts = SpaceType::MeshType::template IndexSet<dim, 0>::Type::num_indices;
        // define local arrays
        typename SpaceHelp::JacobianMatrixType loc_jac, loc_jac_inv;
        typename SpaceHelp::EvalData basis_data;

        // local vector and matrix defines
        typedef Tiny::Vector<DataType, dim> VecValueType;
        typedef Tiny::Matrix<DataType, dim, dim> MatValueType;

        VecValueType loc_v(DataType(0));
        MatValueType loc_grad_v(DataType(0)), strain_rate_tensor_2(DataType(0));
        DataType local_delta(DataType(0));
        DataType nu_loc(DataType(0));
        DataType gamma_dot(DataType(0));

        loc_mat.format();


        if constexpr(need_streamline) //need streamdiff? constexpr since we need this relatively rarely
        {
          VecValueType mean_v(DataType(0));
          //set domain point to barycenter, which is zero for Hypercubes
          const VecValueType barycenter(DataType(0));
          // NewSpaceHelp::_map_point(img_point, barycenter, local_coeffs);
          //only reserve memory for reference values
          SpaceHelp::eval_ref_values(basis_data, barycenter);
          SpaceHelp::trans_values(basis_data);
          for(int i(0); i < num_loc_dofs; ++i)
          {
            mean_v.axpy(basis_data.phi[i].value, local_conv_dofs[i]);
          }
          const DataType local_norm_v = mean_v.norm_euclid();

          if(local_norm_v > tol_eps)
          {
            //we need this trafo data -> for what??
            // SpaceHelp::calc_jac_mat(loc_jac, barycenter, local_coeffs);
            // loc_jac_inv.set_inverse(loc_jac);

            const DataType local_h = SpaceHelp::width_directed(mean_v, local_coeffs) * local_norm_v;
            const DataType local_re = (local_norm_v * local_h) / sd_nu;
            local_delta = sd_delta * (local_h / sd_v_norm) * (DataType(2) * local_re) / (DataType(1) + local_re);
          }
        }

        for(int cub_ind = 0; cub_ind < num_cubs; ++cub_ind)
        {
          const typename SpaceHelp::DomainPointType& dom_point = cub_pt[cub_ind];
          SpaceHelp::eval_ref_values(basis_data, dom_point);
          SpaceHelp::trans_values(basis_data);
          // NewSpaceHelp::_map_point(img_point, dom_point, local_coeffs);//not needed?
          SpaceHelp::calc_jac_mat(loc_jac, dom_point, local_coeffs);
          loc_jac_inv.set_inverse(loc_jac);

          SpaceHelp::eval_ref_gradients(basis_data, dom_point);
          SpaceHelp::trans_gradients(basis_data, loc_jac_inv);

          const DataType weight = loc_jac.det() * cub_wg[cub_ind];
          if(need_convection || need_streamline) // need streamdiff or convection?
          {
            loc_v.format();
            for(int i = 0; i < num_loc_dofs; ++i)
            {
              loc_v.axpy(basis_data.phi[i].value, local_conv_dofs[i]);
            }
          }
          // assemble loc_grad_v, always required for this assembler
          {
            loc_grad_v.format();
            for(int i = 0; i < num_loc_dofs; ++i)
            {
              loc_grad_v.add_outer_product(local_conv_dofs[i], basis_data.phi[i].grad);
            }
          }
          // required terms for carreau assembly
          strain_rate_tensor_2.set_transpose(loc_grad_v);
          strain_rate_tensor_2.axpy(DataType(1.0), loc_grad_v);
          #ifdef __CUDACC__
          gamma_dot = CudaMath::cuda_sqrt(DataType(0.5))*strain_rate_tensor_2.norm_frobenius();
          #else
          gamma_dot = Math::sqrt(DataType(0.5))*strain_rate_tensor_2.norm_frobenius();
          #endif

          //default to nu ?
          nu_loc = visco_func(gamma_dot);


          // deformation assembly, always required
          for(int i(0); i < num_loc_dofs; ++i)
          {
            // trial function loop
            for(int j(0); j < num_loc_dofs; ++j)
            {
              // compute scalar value
              const DataType value = nu_loc * weight * Tiny::dot(basis_data.phi[i].grad, basis_data.phi[j].grad);

              // update local matrix
              loc_mat[i][j].add_scalar_main_diag(value);
              loc_mat[i][j].add_outer_product(basis_data.phi[j].grad, basis_data.phi[i].grad, nu_loc * weight);
            }
          }
          // frechet carreau term "only" makes sense in defo formulation
          if(need_frechet_material)
          {
            const DataType fac =  frechet_material * visco_d_func(gamma_dot)/(gamma_dot + reg_eps);
            //std::cout << fac ;
            // test function loop
            for(int i(0); i < num_loc_dofs; ++i)
            {
              Tiny::Vector<DataType, dim> du_grad_phi;
              du_grad_phi.set_mat_vec_mult(strain_rate_tensor_2, basis_data.phi[i].grad);

              // trial function loop
              for(int j(0); j < num_loc_dofs; ++j)
              {
                Tiny::Vector<DataType, dim> du_grad_psi;
                du_grad_psi.set_mat_vec_mult(strain_rate_tensor_2, basis_data.phi[j].grad);
                // add outer product of grad(phi) and grad(psi)
                loc_mat[i][j].add_outer_product(du_grad_phi, du_grad_psi, fac * weight);
              }
            }
          }

          if(need_convection) // assemble convection?
          {
            // test function loop
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // trial function loop
              for(int j(0); j < num_loc_dofs; ++j)
              {
                // compute scalar value
                const DataType value = beta * weight * basis_data.phi[i].value * Tiny::dot(loc_v, basis_data.phi[j].grad);

                // update local matrix
                loc_mat[i][j].add_scalar_main_diag(value);
              }
            }
          }

          // assemble convection Frechet?
          #ifdef __CUDACC__
          if(CudaMath::cuda_abs(frechet_beta) > DataType(0))
          #else
          if(Math::abs(frechet_beta) > DataType(0))
          #endif
          {
            // test function loop
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // trial function loop
              for(int j(0); j < num_loc_dofs; ++j)
              {
                // compute scalar value
                const DataType value = frechet_beta * weight * basis_data.phi[i].value * basis_data.phi[j].value;

                // update local matrix
                loc_mat[i][j].axpy(value, loc_grad_v);
              }
            }
          }

          // assemble reaction?
          #ifdef __CUDACC__
          if(CudaMath::cuda_abs(theta) > DataType(0))
          #else
          if(Math::abs(theta) > DataType(0))
          #endif
          {
            // test function loop
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // trial function loop
              for(int j(0); j < num_loc_dofs; ++j)
              {
                // compute scalar value
                const DataType value = theta * weight *  basis_data.phi[i].value * basis_data.phi[j].value;

                // update local matrix
                loc_mat[i][j].add_scalar_main_diag(value);
              }
            }
          }

          if constexpr(need_streamline) // need streamdiff?
          {
            Tiny::Vector<DataType, num_loc_dofs> streamdiff_coeffs(DataType(0));
            for(int i = 0; i < num_loc_dofs; ++i)
            {
              streamdiff_coeffs[i] = Tiny::dot(loc_v, basis_data.phi[i].grad);
            }


            // assemble streamline diffusion?
            if((local_delta > tol_eps))
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = local_delta * weight * streamdiff_coeffs[i] * streamdiff_coeffs[j];

                  // update local matrix
                  loc_mat[i][j].add_scalar_main_diag(value);
                }
              }
            }
          }
        // next cubature point
        }
      }

      #if defined(__CUDACC__)
      /**
       * \brief Grouped Burgers Material Matrix assembly kernel
       *
       * Outsources the the most inner execution of the assembly of the full burgers term for a cell previously gathered.
       * Local always refers to the local dofs of the current cell.
       *
       * \tparam ThreadGroup_ The thread group the assembly is synchronized on.
       * \tparam SpaceHelp_ The spacehelper to be used. Holds most information about the underlying space.
       * \tparam ViscFunc_ Function mapping from R to R.
       * \tparam ViscDerFunc_ Function mapping from R to R.
       * \tparam need_stream_diff_ Do we require streamline diffusion.
       *
       * \param[out] loc_mat A reference to the local matrix to be assembled.
       * \param[in] local_conv_dofs A reference to the local convection dofs. If convection is not needed values are ignored.
       * \param[in] local_coeffs The local trafo coefficients.
       * \param[in] cub_pt A pointer to the cubature point coordinate vector.
       * \param[in] cub_wg A pointer to the cubature weights. Same order as cub_pt array.
       * \param[in] burgers_params A struct holding the burgers parameter configuration.
       * \param[in] material_params A struct holding the additional material parameters.
       * \param[in] need_convection Do we need convection?
       * \param[in] tol_eps Tolerance for local stream norm to be regarded as zero.
       * \param[in] visco_func Function mapping invariants of velocity jacobian to viscosity.
       * \param[in] visco_d_func Function representing the derivative of visc_func.
       */
      template<typename ThreadGroup_, typename SpaceHelp_, typename ViscFunc_, typename ViscDerFunc_, bool need_stream_diff_>
      CUDA_DEVICE __forceinline__ void grouped_burgers_velo_material_mat_assembly_kernel(const ThreadGroup_& tg, BurgersMatSharedDataKernelWrapper<SpaceHelp_, need_stream_diff_>* shared_wrapper,
                                          int loc_assemble_size, int assemble_offset, typename SpaceHelp_::DataType* loc_mat, const typename SpaceHelp_::DataType* local_conv_dofs,
                                          const Tiny::Matrix<typename SpaceHelp_::DataType, SpaceHelp_::dim, SpaceHelp_::num_verts>& local_coeffs,
                                          const typename SpaceHelp_::DomainPointType* cub_pt, const typename SpaceHelp_::DataType* cub_wg, const int num_cubs,
                                          const VoxelAssembly::AssemblyBurgersData<typename SpaceHelp_::DataType>& burgers_params,
                                          const VoxelAssembly::AssemblyMaterialData<typename SpaceHelp_::DataType>& material_params,
                                          const bool need_convection, const typename SpaceHelp_::DataType tol_eps,
                                          ViscFunc_ visco_func, ViscDerFunc_ visco_d_func)
      {
        typedef SpaceHelp_ SpaceHelp;
        constexpr int dim = SpaceHelp::dim;
        typedef typename SpaceHelp::SpaceType SpaceType;
        typedef typename SpaceHelp::DataType DataType;
        constexpr bool need_streamline = need_stream_diff_;
        // constexpr int num_loc_verts = SpaceHelp::num_verts;
        const int t_idx = tg.thread_rank();
        const int g_size = tg.num_threads();

        // const DataType& nu{burgers_params.nu};
        const DataType& theta{burgers_params.theta};
        const DataType& beta{burgers_params.beta};
        const DataType& frechet_beta{burgers_params.frechet_beta};
        const DataType& sd_delta{burgers_params.sd_delta};
        const DataType& sd_nu{burgers_params.sd_nu};
        const DataType& sd_v_norm{burgers_params.sd_v_norm};

        const DataType& frechet_material{material_params.frechet_material};
        const DataType& reg_eps{material_params.reg_eps};
        const bool& need_frechet_material{material_params.need_frechet_material};

        // #ifdef __CUDACC__
        // const DataType tol_eps = CudaMath::cuda_get_sqrt_eps<DataType>();
        // #else
        // const DataType tol_eps = Math::sqrt(Math::eps<DataType>());
        // #endif

        // #ifdef __CUDACC__
        // const bool need_streamline = (CudaMath::cuda_abs(sd_delta) > DataType(0)) && (sd_v_norm > tol_eps);
        // const bool need_convection = CudaMath::cuda_abs(beta) > DataType(0);
        // #else
        // const bool need_streamline = (Math::abs(sd_delta) > DataType(0)) && (sd_v_norm > tol_eps);
        // const bool need_convection = Math::abs(beta) > DataType(0);
        // #endif


        //define local sizes
        constexpr int num_loc_dofs = SpaceType::DofMappingType::dof_count;
        // local vector and matrix defines
        typedef Tiny::Vector<DataType, dim> VecValueType;
        typedef Tiny::Matrix<DataType, dim, dim> MatValueType;
        //get number of nodes per element
        // constexpr int num_loc_verts = SpaceType::MeshType::template IndexSet<dim, 0>::Type::num_indices;
        // define local arrays
        // use extern shared for this?!
        typedef BurgersMatSharedDataKernelWrapper<SpaceHelp, need_stream_diff_> BSDKWrapper;
        typename SpaceHelp::JacobianMatrixType& loc_jac = shared_wrapper->loc_jac;
        typename SpaceHelp::JacobianMatrixType& loc_jac_inv = shared_wrapper->loc_jac_inv;
        MatValueType& loc_grad_v = shared_wrapper->loc_grad_v;
        VecValueType& loc_v = shared_wrapper->loc_v;
        VecValueType* mean_v_p = nullptr;
        if constexpr(need_streamline) mean_v_p = &(shared_wrapper->mean_v);
        typename SpaceHelp::DomainPointType& dom_point = shared_wrapper->dom_point;
        Tiny::Vector<DataType, num_loc_dofs>* streamdiff_coeffs_p = nullptr;
        if constexpr(need_streamline) streamdiff_coeffs_p = &(shared_wrapper->streamdiff_coeffs);
        DataType& local_delta = shared_wrapper->local_delta;
        DataType& nu_loc = shared_wrapper->nu_loc;
        DataType& gamma_dot = shared_wrapper->gamma_dot;
        DataType& det = shared_wrapper->det;
        DataType& weight = shared_wrapper->weight;
        typename SpaceHelp::EvalData& basis_data = shared_wrapper->basis_data;
        bool& need_frechet = shared_wrapper->need_frechet;

        // local tmp strainrate value
        MatValueType strain_rate_tensor_2(DataType(0));

        VoxelAssembly::coalesced_format(tg, (unsigned int*) shared_wrapper, sizeof(BSDKWrapper)/(sizeof(unsigned int)));
        tg.sync();

        cg::invoke_one(tg, [&](){need_frechet = CudaMath::cuda_abs(frechet_beta) > DataType(0);});
        tg.sync();


        if constexpr(need_streamline) //need streamdiff?
        {
          cg::invoke_one(tg, [&]() {*mean_v_p = DataType(0);});
          //set domain point to barycenter, which is zero for Hypercubes
          // NewSpaceHelp::_map_point(img_point, barycenter, local_coeffs);
          //only reserve memory for reference values
          cg::invoke_one(tg, [&](){ const VecValueType bc(0);
                                    SpaceHelp::eval_ref_values(basis_data, bc);
                                    SpaceHelp::trans_values(basis_data);});

          tg.sync();

          for(int idx = t_idx; idx < num_loc_dofs; idx += g_size)
          {
            VoxelAssembly::template grouped_axpy<DataType, cg::thread_group, dim>(tg, &(*mean_v_p)[0], basis_data.phi[idx].value, &local_conv_dofs[idx*dim]);
          }

          tg.sync();

          DataType local_norm_v = mean_v_p->norm_euclid();

          if(local_norm_v > tol_eps)
          {
            cg::invoke_one(tg, [&](){
              const DataType local_h = SpaceHelp::width_directed(*mean_v_p, local_coeffs) * local_norm_v;
              const DataType local_re = (local_norm_v * local_h) / sd_nu;
              local_delta = sd_delta * (local_h / sd_v_norm) * (DataType(2) * local_re) / (DataType(1) + local_re);
            });
          }
        }

        tg.sync();


        for(int cub_ind = 0; cub_ind < num_cubs; ++cub_ind)
        {
          cg::invoke_one(tg, [&](){dom_point = cub_pt[cub_ind];
                                  SpaceHelp::eval_ref_values(basis_data, dom_point);
                                  SpaceHelp::trans_values(basis_data);});
          // NewSpaceHelp::_map_point(img_point, dom_point, local_coeffs);//not needed?
          tg.sync();
          SpaceHelp::grouped_calc_jac_mat(tg, loc_jac, dom_point, &local_coeffs[0][0]);
          tg.sync();
          cg::invoke_one(tg, [&](){det = loc_jac.det();
                                  weight = det * cub_wg[cub_ind];});
          tg.sync();
          loc_jac_inv.grouped_set_inverse(tg, loc_jac, det);
          tg.sync();
          cg::invoke_one(tg, [&](){
                                  SpaceHelp::eval_ref_gradients(basis_data, dom_point);
                                  SpaceHelp::trans_gradients(basis_data, loc_jac_inv);});

          tg.sync();
          if(need_convection || need_streamline) // need streamdiff or convection?
          {
            VoxelAssembly::coalesced_format(tg, &loc_v[0], dim);
            tg.sync();

            for(int idx = t_idx; idx < num_loc_dofs; idx += g_size)
            {
              VoxelAssembly::template grouped_axpy<DataType, cg::thread_group, dim>(tg, &loc_v[0], basis_data.phi[idx].value, &local_conv_dofs[idx*dim]);
            }
            // tg.sync();
          }

          {
            VoxelAssembly::coalesced_format(tg, &loc_grad_v[0][0], dim*dim);
            tg.sync();
            for(int i = t_idx; i < num_loc_dofs; i += g_size)
            {
              VoxelAssembly::grouped_add_outer_product(tg, &loc_grad_v[0][0], &local_conv_dofs[i*dim], basis_data.phi[i].grad);
            }
            // tg.sync();
          }
          tg.sync();
          // assemble necessary values for material assembler
          {
              strain_rate_tensor_2.set_transpose(loc_grad_v);
              strain_rate_tensor_2.axpy(DataType(1.0), loc_grad_v);
            cg::invoke_one(tg, [&](){
              // strain_rate_tensor_2.set_transpose(loc_grad_v);
              // strain_rate_tensor_2.axpy(DataType(1.0), loc_grad_v);
              gamma_dot = CudaMath::cuda_sqrt(DataType(0.5))*strain_rate_tensor_2.norm_frobenius();
              nu_loc = visco_func(gamma_dot);
            });
          }


          if constexpr(need_streamline) // need streamdiff?
          {
            tg.sync();
            for(int i = t_idx; i < num_loc_dofs; i += g_size)
            {
              (*streamdiff_coeffs_p)[i] = Tiny::dot(loc_v, basis_data.phi[i].grad);
            }
            // tg.sync();
          }

          tg.sync();

          const DataType fac =  frechet_material * visco_d_func(gamma_dot)/(gamma_dot + reg_eps);

          for(int idx = t_idx; idx < loc_assemble_size; idx += g_size)
          {
            // the thread local test and trial function indices
            const int i = (idx+assemble_offset) / num_loc_dofs;
            const int j = (idx+assemble_offset) % num_loc_dofs;

            {
              // compute scalar value
              const DataType value = nu_loc * weight * Tiny::dot(basis_data.phi[i].grad, basis_data.phi[j].grad);

              ((Tiny::Matrix<DataType, dim, dim>*)(loc_mat + dim*dim*idx))->add_scalar_main_diag(value);
              ((Tiny::Matrix<DataType, dim, dim>*)(loc_mat + dim*dim*idx))->add_outer_product(basis_data.phi[j].grad, basis_data.phi[i].grad, nu_loc * weight);
            }

            if(need_frechet_material)
            {
              Tiny::Vector<DataType, dim> du_grad_phi;
              du_grad_phi.set_mat_vec_mult(strain_rate_tensor_2, basis_data.phi[i].grad);
              Tiny::Vector<DataType, dim> du_grad_psi;
              du_grad_psi.set_mat_vec_mult(strain_rate_tensor_2, basis_data.phi[j].grad);
              // printf("from thread %i du grad psi %f, %f \n", idx, du_grad_psi[0], du_grad_psi[1]);
              ((Tiny::Matrix<DataType, dim, dim>*)(loc_mat + dim*dim*idx))->add_outer_product(du_grad_phi, du_grad_psi, fac*weight);

            }

            if(need_convection) // assemble convection?
            {
              // compute scalar value
              const DataType value = beta * weight * basis_data.phi[i].value * VoxelAssembly::dot(basis_data.phi[j].grad, &loc_v[0]);

              // update local matrix
              ((Tiny::Matrix<DataType, dim, dim>*)(loc_mat + dim*dim*idx))->add_scalar_main_diag(value);
            }

            // assemble convection Frechet?
            if(need_frechet)
            {
              // compute scalar value
              const DataType value = frechet_beta * weight * basis_data.phi[i].value * basis_data.phi[j].value;

              // update local matrix
              ((Tiny::Matrix<DataType, dim, dim>*)(loc_mat + dim*dim*idx))->axpy(value, loc_grad_v);
            }

            // assemble reaction?
            #ifdef __CUDACC__
            if(CudaMath::cuda_abs(theta) > DataType(0))
            #else
            if(Math::abs(theta) > DataType(0))
            #endif
            {
              // compute scalar value
              const DataType value = theta * weight *  basis_data.phi[i].value * basis_data.phi[j].value;

              // update local matrix
              ((Tiny::Matrix<DataType, dim, dim>*)(loc_mat + dim*dim*idx))->add_scalar_main_diag(value);
            }

            // assemble streamline diffusion?
            if constexpr(need_streamline)
            {
              if(local_delta > tol_eps)
              {
                // compute scalar value
                const DataType value = local_delta * weight * (*streamdiff_coeffs_p)[i] * (*streamdiff_coeffs_p)[j];

                // update local matrix
                ((Tiny::Matrix<DataType, dim, dim>*)(loc_mat + dim*dim*idx))->add_scalar_main_diag(value);
              }
            }
          }
        // next cubature point
        }
      }
      #endif

      /**
       * \brief Burgers Vector assembly kernel
       *
       * Outsources the the most inner execution of the assembly of the full burgers term for a cell previously gathered.
       * Local always refers to the local dofs of the current cell.
       *
       * \tparam SpaceHelp_ The psacehelper to be used. Holds most information about the underlying space.
       * \tparam LocMatType_ The Matrixtype of the local matrix dofs.
       * \tparam LocVecType_ The Vectortype of the local convection dofs.
       * \tparam dim_ Dimension of our space. For ease of access.
       * \tparam num_verts_ The number of vertices of one cell.
       * \tparam ViscFunc_ Function mapping from R to R.
       *
       * \param[out] loc_vec A reference to the local defect to be assembled.
       * \param[in] local_prim_dofs A reference to the local primal values the matrix is implicitly applied to.
       * \param[in] local_conv_dofs A reference to the local convection dofs. If convection is not needed values are ignored.
       * \param[in] local_coeffs The local trafo coefficients.
       * \param[in] cub_pt A pointer to the cubature point coordinate vector.
       * \param[in] cub_wg A pointer to the cubature weights. Same order as cub_pt array.
       * \param[in] burgers_params A struct holding the burgers parameter configuration.
       * \param[in] need_convection Do we need convection?
       * \param[in] visco_func Function mapping invariants of velocity jacobian to viscosity.
       */
      template<typename SpaceHelp_, typename LocVecType_, int dim_, int num_verts_, typename ViscFunc_>
      CUDA_HOST_DEVICE void burgers_velo_material_defect_assembly_kernel(LocVecType_& loc_vec, const LocVecType_& local_prim_dofs, const LocVecType_& local_conv_dofs, const Tiny::Matrix<typename SpaceHelp_::DataType, dim_, num_verts_>& local_coeffs,
                                          const typename SpaceHelp_::DomainPointType* cub_pt, const typename SpaceHelp_::DataType* cub_wg, const int num_cubs,
                                          const VoxelAssembly::AssemblyBurgersData<typename SpaceHelp_::DataType>& burgers_params,
                                          const bool need_convection, ViscFunc_ visco_func)
      {
        typedef SpaceHelp_ SpaceHelp;
        constexpr int dim = SpaceHelp::dim;
        typedef typename SpaceHelp::SpaceType SpaceType;
        typedef typename SpaceHelp::DataType DataType;

        // const DataType& nu{burgers_params.nu};
        const DataType& theta{burgers_params.theta};
        const DataType& beta{burgers_params.beta};
        // const DataType& frechet_beta{burgers_params.frechet_beta};
        // const DataType& sd_delta{burgers_params.sd_delta};
        // const DataType& sd_nu{burgers_params.sd_nu};
        // const DataType& sd_v_norm{burgers_params.sd_v_norm};
        // const bool& deformation{burgers_params.deformation};

        // #ifdef __CUDACC__
        // const DataType tol_eps = CudaMath::cuda_get_sqrt_eps<DataType>();
        // #else
        // const DataType tol_eps = Math::sqrt(Math::eps<DataType>());
        // #endif

        // #ifdef __CUDACC__
        // const bool need_streamline = (CudaMath::cuda_abs(sd_delta) > DataType(0)) && (sd_v_norm > tol_eps);
        // const bool need_convection = CudaMath::cuda_abs(beta) > DataType(0);
        // #else
        // const bool need_streamline = (Math::abs(sd_delta) > DataType(0)) && (sd_v_norm > tol_eps);
        // const bool need_convection = Math::abs(beta) > DataType(0);
        // #endif


        //define local sizes
        constexpr int num_loc_dofs = SpaceType::DofMappingType::dof_count;
        //get number of nodes per element
        // constexpr int num_loc_verts = SpaceType::MeshType::template IndexSet<dim, 0>::Type::num_indices;
        // define local arrays
        typename SpaceHelp::JacobianMatrixType loc_jac, loc_jac_inv;
        typename SpaceHelp::EvalData basis_data;

        // local vector and matrix defines
        typedef Tiny::Vector<DataType, dim> VecValueType;
        typedef Tiny::Matrix<DataType, dim, dim> MatValueType;

        VecValueType loc_v(DataType(0));
        DataType nu_loc(DataType(0));

        loc_vec.format();

        Tiny::Vector<DataType, num_loc_dofs> streamdiff_coeffs(DataType(0));

        for(int cub_ind = 0; cub_ind < num_cubs; ++cub_ind)
        {
          const typename SpaceHelp::DomainPointType& dom_point = cub_pt[cub_ind];
          SpaceHelp::eval_ref_values(basis_data, dom_point);
          SpaceHelp::trans_values(basis_data);
          // NewSpaceHelp::_map_point(img_point, dom_point, local_coeffs);//not needed?
          SpaceHelp::calc_jac_mat(loc_jac, dom_point, local_coeffs);
          loc_jac_inv.set_inverse(loc_jac);

          SpaceHelp::eval_ref_gradients(basis_data, dom_point);
          SpaceHelp::trans_gradients(basis_data, loc_jac_inv);

          const DataType weight = loc_jac.det() * cub_wg[cub_ind];
          if(need_convection) // need streamdiff or convection?
          {
            loc_v.format();
            for(int i = 0; i < num_loc_dofs; ++i)
            {
              loc_v.axpy(basis_data.phi[i].value, local_conv_dofs[i]);
            }
          }

          // assemble loc_grad_v, always required for this assembler
          {
            MatValueType loc_grad_v(DataType(0)), strain_rate_tensor_2(DataType(0));
            loc_grad_v.format();
            for(int i = 0; i < num_loc_dofs; ++i)
            {
              loc_grad_v.add_outer_product(local_conv_dofs[i], basis_data.phi[i].grad);
            }
            // required terms for carreau assembly
            strain_rate_tensor_2.set_transpose(loc_grad_v);
            strain_rate_tensor_2.axpy(DataType(1.0), loc_grad_v);
            #ifdef __CUDACC__
            nu_loc = visco_func(CudaMath::cuda_sqrt(DataType(0.5))*strain_rate_tensor_2.norm_frobenius());
            #else
            nu_loc = visco_func(Math::sqrt(DataType(0.5))*strain_rate_tensor_2.norm_frobenius());
            #endif

          }

          // always deformation
          for(int i(0); i < num_loc_dofs; ++i)
          {
            // trial function loop
            for(int j(0); j < num_loc_dofs; ++j)
            {
              // compute scalar values
              const DataType value1 = nu_loc * weight * Tiny::dot(basis_data.phi[i].grad, basis_data.phi[j].grad);
              const DataType value2 = nu_loc * weight * Tiny::dot(local_prim_dofs[j], basis_data.phi[i].grad);
              // update local vector
              loc_vec[i].axpy(value1, local_prim_dofs[j]);
              loc_vec[i].axpy(value2, basis_data.phi[j].grad);
            }
          }

          if(need_convection) // assemble convection?
          {
            // test function loop
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // trial function loop
              for(int j(0); j < num_loc_dofs; ++j)
              {
                // compute scalar value
                const DataType value = beta * weight * basis_data.phi[i].value * Tiny::dot(loc_v, basis_data.phi[j].grad);

                // update local matrix
                loc_vec[i].axpy(value, local_prim_dofs[j]);
              }
            }
          }

          // assemble reaction?
          #ifdef __CUDACC__
          if(CudaMath::cuda_abs(theta) > DataType(0))
          #else
          if(Math::abs(theta) > DataType(0))
          #endif
          {
            // test function loop
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // trial function loop
              for(int j(0); j < num_loc_dofs; ++j)
              {
                // compute scalar value
                const DataType value = theta * weight *  basis_data.phi[i].value * basis_data.phi[j].value;

                // update local matrix
                loc_vec[i].axpy(value, local_prim_dofs[j]);
              }
            }
          }
        // next cubature point
        }
      }

     #ifdef __CUDACC__
      /**
       * \brief Goruped Burgers Vector assembly kernel
       *
       * Outsources the the most inner execution of the assembly of the full burgers term for a cell previously gathered.
       * Local always refers to the local dofs of the current cell.
       *
       * \tparam SpaceHelp_ The psacehelper to be used. Holds most information about the underlying space.
       * \tparam ViscFunc_ Function mapping from R to R.
       *
       */
      template<typename ThreadGroup_, typename SpaceHelp_, typename ViscFunc_>
      CUDA_DEVICE void grouped_burgers_velo_material_defect_assembly_kernel(const ThreadGroup_& tg, BurgersMatSharedDataKernelWrapper<SpaceHelp_, false>* shared_wrapper, int loc_assemble_size, int assemble_offset,
                                          typename SpaceHelp_::DataType* loc_vec, const typename SpaceHelp_::DataType* local_prim_dofs, const typename SpaceHelp_::DataType* local_conv_dofs,
                                          const Tiny::Matrix<typename SpaceHelp_::DataType, SpaceHelp_::dim, SpaceHelp_::num_verts>& local_coeffs,
                                          const typename SpaceHelp_::DomainPointType* cub_pt, const typename SpaceHelp_::DataType* cub_wg, const int num_cubs,
                                          const VoxelAssembly::AssemblyBurgersData<typename SpaceHelp_::DataType>& burgers_params,
                                          const bool need_convection, ViscFunc_ visco_func)
      {
        typedef SpaceHelp_ SpaceHelp;
        constexpr int dim = SpaceHelp::dim;
        typedef typename SpaceHelp::SpaceType SpaceType;
        typedef typename SpaceHelp::DataType DataType;
        const int t_idx = tg.thread_rank();
        const int g_size = tg.num_threads();

        const DataType& nu{burgers_params.nu};
        const DataType& theta{burgers_params.theta};
        const DataType& beta{burgers_params.beta};
        // const DataType& frechet_beta{burgers_params.frechet_beta};
        // const DataType& sd_delta{burgers_params.sd_delta};
        // const DataType& sd_nu{burgers_params.sd_nu};
        // const DataType& sd_v_norm{burgers_params.sd_v_norm};
        // const bool& deformation{burgers_params.deformation};

        // #ifdef __CUDACC__
        // const DataType tol_eps = CudaMath::cuda_get_sqrt_eps<DataType>();
        // #else
        // const DataType tol_eps = Math::sqrt(Math::eps<DataType>());
        // #endif

        // #ifdef __CUDACC__
        // const bool need_streamline = (CudaMath::cuda_abs(sd_delta) > DataType(0)) && (sd_v_norm > tol_eps);
        // const bool need_convection = CudaMath::cuda_abs(beta) > DataType(0);
        // #else
        // const bool need_streamline = (Math::abs(sd_delta) > DataType(0)) && (sd_v_norm > tol_eps);
        // const bool need_convection = Math::abs(beta) > DataType(0);
        // #endif


        //define local sizes
        constexpr int num_loc_dofs = SpaceType::DofMappingType::dof_count;
        //get number of nodes per element
        // constexpr int num_loc_verts = SpaceType::MeshType::template IndexSet<dim, 0>::Type::num_indices;
        // define local arrays
        typedef BurgersMatSharedDataKernelWrapper<SpaceHelp, false> BSDKWrapper;
        typename SpaceHelp::JacobianMatrixType& loc_jac = shared_wrapper->loc_jac;
        typename SpaceHelp::JacobianMatrixType& loc_jac_inv = shared_wrapper->loc_jac_inv;

        DataType& det = shared_wrapper->det;
        DataType& weight = shared_wrapper->weight;
        DataType& nu_loc = shared_wrapper->nu_loc;
        DataType& gamma_dot = shared_wrapper->gamma_dot;
        typename SpaceHelp::EvalData& basis_data = shared_wrapper->basis_data;
        // local vector and matrix defines
        typedef Tiny::Vector<DataType, dim> VecValueType;
        typedef Tiny::Matrix<DataType, dim, dim> MatValueType;

        MatValueType& loc_grad_v = shared_wrapper->loc_grad_v;
        VecValueType& loc_v = shared_wrapper->loc_v;

        typename SpaceHelp::DomainPointType& dom_point = shared_wrapper->dom_point;

        VoxelAssembly::coalesced_format(tg, (unsigned int*) shared_wrapper, sizeof(BSDKWrapper)/(sizeof(unsigned int)));
        tg.sync();
        // 4 threads should work on one entry, requires block size to be mutliple of 4
        cg::thread_block_tile<4> tile4 = cg::tiled_partition<4>(tg);

        for(int cub_ind = 0; cub_ind < num_cubs; ++cub_ind)
        {
          tg.sync();
          cg::invoke_one(tg, [&](){dom_point = cub_pt[cub_ind];
                                  SpaceHelp::eval_ref_values(basis_data, dom_point);
                                  SpaceHelp::trans_values(basis_data);});
          // NewSpaceHelp::_map_point(img_point, dom_point, local_coeffs);//not needed?
          tg.sync();
          SpaceHelp::grouped_calc_jac_mat(tg, loc_jac, dom_point, &local_coeffs[0][0]);
          tg.sync();
          cg::invoke_one(tg, [&](){det = loc_jac.det();
          weight = det * cub_wg[cub_ind];});
          tg.sync();
          loc_jac_inv.grouped_set_inverse(tg, loc_jac, det);
          tg.sync();
          cg::invoke_one(tg, [&](){
                                  SpaceHelp::eval_ref_gradients(basis_data, dom_point);
                                  SpaceHelp::trans_gradients(basis_data, loc_jac_inv);});

          tg.sync();
          if(need_convection) // need streamdiff or convection?
          {
            VoxelAssembly::coalesced_format(tg, &loc_v[0], dim);
            tg.sync();

            for(int idx = t_idx; idx < num_loc_dofs; idx += g_size)
            {
              VoxelAssembly::template grouped_axpy<DataType, cg::thread_group, dim>(tg, &loc_v[0], basis_data.phi[idx].value, &local_conv_dofs[idx*dim]);
            }
          }
          tg.sync();
          // assemble necessary values for material assembler
          {
            VoxelAssembly::coalesced_format(tg, &loc_grad_v[0][0], dim*dim);
            tg.sync();
            for(int i = t_idx; i < num_loc_dofs; i += g_size)
            {
              VoxelAssembly::grouped_add_outer_product(tg, &loc_grad_v[0][0], &local_conv_dofs[i*dim], basis_data.phi[i].grad);
            }
            tg.sync();
            cg::invoke_one(tg, [&](){
              MatValueType strain_rate_tensor_2(DataType(0));
              strain_rate_tensor_2.set_transpose(loc_grad_v);
              strain_rate_tensor_2.axpy(DataType(1.0), loc_grad_v);
              gamma_dot = CudaMath::cuda_sqrt(DataType(0.5))*strain_rate_tensor_2.norm_frobenius();
              nu_loc = visco_func(gamma_dot);
            });
          }

          tg.sync();

          for(int idx = tile4.meta_group_rank(); idx < loc_assemble_size*dim; idx += tile4.meta_group_size())
          {
            // compute scalar value
            DataType val = DataType(0);
            // the thread local test function index
            const int i = (idx/dim+assemble_offset);
            const int ii = idx%dim;
            {
              for(int j = tile4.thread_rank(); j < num_loc_dofs; j += tile4.size())
              {
                // compute scalar values
                val += nu_loc * weight * (Tiny::dot(basis_data.phi[i].grad, basis_data.phi[j].grad) * local_prim_dofs[j*dim + ii]
                      + Tiny::dot(((Tiny::Vector<DataType, dim>*)local_prim_dofs)[j], basis_data.phi[i].grad) * basis_data.phi[j].grad[ii]);
              }
            }

            if(need_convection) // assemble convection?
            {
              for(int j = tile4.thread_rank(); j < num_loc_dofs; j += tile4.size())
              {
                // compute scalar values
                val += beta * weight * basis_data.phi[i].value * Tiny::dot(loc_v, basis_data.phi[j].grad) * local_prim_dofs[j*dim + ii];
              }
            }

            // assemble reaction?
            #ifdef __CUDACC__
            if(CudaMath::cuda_abs(theta) > DataType(0))
            #else
            if(Math::abs(theta) > DataType(0))
            #endif
            {
              for(int j = tile4.thread_rank(); j < num_loc_dofs; j += tile4.size())
              {
                // compute scalar values
                val += theta * weight *  basis_data.phi[i].value * basis_data.phi[j].value * local_prim_dofs[j*dim + ii];
              }
            }
            // reduce result and store value with async call
            // tile4.sync();
            cuda::atomic_ref<DataType, cuda::thread_scope_block> a_ref(*(loc_vec+idx));
            cg::reduce_update_async(tile4, a_ref, val, cg::plus<DataType>());
          }
        // next cubature point
        }
      }
      #endif
    } // namespace Kernel

    namespace Arch
    {
      #if defined(FEAT_HAVE_CUDA) || defined(DOXYGEN)
      /**
       * \brief Device kernel wrapper for the full matrix burgers assembler.
       *
       * Assembles the standard burgers operator on the device side in an additive manner.
       *
       * \tparam Space_ The underlying spacetype.
       * \tparam DT_ The datatype.
       * \tparam IT_ The indextype.
       *
       * \param[in] space The actual space we use.
       * \param[in/out] matrix_data Holds the CSR typed matrix data. The const only refers to the inner pointer, not the data they refer to.
       * \param[in] conv_data The data array of the convection vector.
       * \param[in] cubature A wrapper around the cubature rule data. Datapointer have to be device pointer.
       * \param[in] dof_mapping A wrapper around the dof mapping for scattering gathering of local dofs. Datapointer have to be devicepointer.
       * \param[in] coloring_maps A vector of arrays of the dof indices per color.
       * \param[in] coloring_map_sizes A vector of the coloring map sizes.
       * \param[in] alpha Scaling parameter for gathering.
       * \param[in] burgers_params A struct holding all burgers parameters.
       * \param[in] material_params Additonal material data.
       * \param[in] material_type Function determining the material model used.
       */
      template<typename Space_, typename DT_, typename IT_>
      void assemble_burgers_velo_material_csr(const Space_& space,
              const CSRMatrixData<DT_, IT_>& matrix_data,
              const DT_* conv_data,
              const AssemblyCubatureData<DT_>& cubature,
              const AssemblyMappingData<DT_, IT_>& dof_mapping,
              const std::vector<int*>& coloring_maps,
              const std::vector<Index>& coloring_map_sizes,
              DT_ alpha, const AssemblyBurgersData<DT_>& burgers_params,
              const AssemblyMaterialData<DT_>& material_params,
              MaterialType mat_type, int shared_mem, int blocksize, int gridsize, bool print_occupancy);

      /**
       * \brief Device kernel wrapper for the defect burgers assembler.
       *
       * Assembles the defect for the standard burgers operator on the device side in an additive manner.
       *
       * \tparam Space_ The underlying spacetype.
       * \tparam DT_ The datatype.
       * \tparam IT_ The indextype.
       *
       * \param[in] space The actual space we use.
       * \param[in/out] vector_data Holds the (device side) data of the defect vector to be assembled. Data is added onto.
       * \param[in] conv_data The data array of the convection vector.
       * \param[in] primal_data The data array of the primal vector the matrix is applied to.
       * \param[in] cubature A wrapper around the cubature rule data. Datapointer have to be device pointer.
       * \param[in] dof_mapping A wrapper around the dof mapping for scattering gathering of local dofs. Datapointer have to be devicepointer.
       * \param[in] coloring_maps A vector of arrays of the dof indices per color.
       * \param[in] coloring_map_sizes A vector of the coloring map sizes.
       * \param[in] alpha Scaling parameter for gathering.
       * \param[in] burgers_params A struct holding all burgers parameters.
       * \param[in] material_params Additonal material data.
       * \param[in] material_type Function determining the material model used.
       */
      template<typename Space_, typename DT_, typename IT_>
      void assemble_burgers_velo_material_defect(const Space_& space,
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
              int shared_mem, int blocksize, int gridsize, bool print_occupancy);

      #endif

      /**
       * \brief Host kernel wrapper for the full matrix burgers assembler.
       *
       * Assembles the extended burgers operator on the host side in an additive manner.
       *
       * \tparam Space_ The underlying spacetype.
       * \tparam DT_ The datatype.
       * \tparam IT_ The indextype.
       *
       * \param[in] space The actual space we use.
       * \param[in/out] matrix_data Holds the CSR typed matrix data. The const only refers to the inner pointer, not the data they refer to.
       * \param[in] conv_data The data array of the convection vector.
       * \param[in] cubature A wrapper around the cubature rule data. Data has to be allocated on host side.
       * \param[in] dof_mapping A wrapper around the dof mapping for scattering gathering of local dofs.
       * \param[in] coloring_maps A vector of arrays of the dof indices per color.
       * \param[in] coloring_map_sizes A vector of the coloring map sizes.
       * \param[in] alpha Scaling parameter for gathering.
       * \param[in] burgers_params A struct holding all burgers parameters.
       * \param[in] material_params Additonal material data.
       * \param[in] material_type Function determining the material model used.
       */
      template<typename Space_, typename DT_, typename IT_>
      void assemble_burgers_velo_material_csr_host(const Space_& space,
              const CSRMatrixData<DT_, IT_>& matrix_data,
              const DT_* conv_data,
              const AssemblyCubatureData<DT_>& cubature,
              const AssemblyMappingData<DT_, IT_>& dof_mapping,
              const std::vector<int*>& coloring_maps,
              const std::vector<Index>& coloring_map_sizes,
              DT_ alpha, const AssemblyBurgersData<DT_>& burgers_params,
              const AssemblyMaterialData<DT_>& material_params,
              MaterialType material_type);

      /**
       * \brief Host kernel wrapper for the defect burgers assembler.
       *
       * Assembles the defect for the standard burgers operator on the host side in an additive manner.
       *
       * \tparam Space_ The underlying spacetype.
       * \tparam DT_ The datatype.
       * \tparam IT_ The indextype.
       *
       * \param[in] space The actual space we use.
       * \param[in/out] vector_data Holds the (host side) data of the defect vector to be assembled. Data is added onto.
       * \param[in] conv_data The data array of the convection vector.
       * \param[in] primal_data The data array of the primal vector the matrix is applied to.
       * \param[in] cubature A wrapper around the cubature rule data.
       * \param[in] dof_mapping A wrapper around the dof mapping for scattering gathering of local dofs.
       * \param[in] coloring_maps A vector of arrays of the dof indices per color.
       * \param[in] coloring_map_sizes A vector of the coloring map sizes.
       * \param[in] alpha Scaling parameter for gathering.
       * \param[in] burgers_params A struct holding all burgers parameters.
       * \param[in] material_params Additonal material data.
       * \param[in] material_type Function determining the material model used.
       */
      template<typename Space_, typename DT_, typename IT_>
      void assemble_burgers_velo_material_defect_host(const Space_& space,
              DT_* vector_data,
              const DT_* conv_data,
              const DT_* primal_data,
              const AssemblyCubatureData<DT_>& cubature,
              const AssemblyMappingData<DT_, IT_>& dof_mapping,
              const std::vector<int*>& coloring_maps,
              const std::vector<Index>& coloring_map_sizes,
              DT_ alpha, const AssemblyBurgersData<DT_>& burgers_params,
              const AssemblyMaterialData<DT_>& material_params,
              MaterialType material_type);


    }


#ifndef __CUDACC__
    /**
     * \brief Q2Lagrange thread parallel assembly class for the burgers velo material operator.
     *
     * This class provides host and device variants for a thread parallel assembly of the extended burgers operator.
     * For details on the burgers velocity material operator, see the burgers velocity material assembler in </kernel/assembly/>.
     * The threading strategy is based on a cell based coloring of the mesh, so in principle this version for Q2
     * space works as long as a valid coloring is provided.
     *
     * \tparam dim_ Dimension of the hypercubes.
     * \tparam DT_ Datatype to be used.
     * \tparam IT_ Indextype to be used.
     */
    template<int dim_, typename DT_, typename IT_>
    class VoxelBurgersVeloMaterialAssembler<Q2StandardHyperCube<dim_>, DT_, IT_> :
     public VoxelBurgersAssembler<Q2StandardHyperCube<dim_>, DT_, IT_>
    {
    public:
      typedef VoxelBurgersAssembler<Q2StandardHyperCube<dim_>, DT_, IT_> BaseClass;
      /// typedefs
      typedef typename BaseClass::SpaceType SpaceType;
      typedef typename BaseClass::DataHandler DataHandler;
      typedef typename BaseClass::SpaceHelp SpaceHelp;
      typedef typename SpaceHelp::ShapeType ShapeType;
      typedef typename SpaceHelp::DataType DataType;
      typedef typename SpaceHelp::IndexType IndexType;

      /// constexpr
      static constexpr int dim = SpaceHelp::dim;

      typedef typename SpaceHelp::DomainPointType DomainPointType;
      typedef typename SpaceHelp::ImagePointType ImagePointType;
      typedef typename SpaceHelp::ValueType ValueType;
      typedef typename SpaceHelp::JacobianMatrixType JacobianMatrixType;

      /// Scaling factor full material jacobian
      DataType frechet_material;

      /// regularity parameter
      DataType reg_eps;

      /// Material Parameter
      DataType mu_0;
      DataType exp;
      DataType lambda;
      DataType a_T;
      DataType yasuda_a;
      DataType mu_inf;

      /// material type
      MaterialType material_type;


    public:
      explicit VoxelBurgersVeloMaterialAssembler() = default;

      /**
       * \brief Constructor for burgers velo material assembler
       *
       * \tparam ColoringType_ The type of the coloring indexset.
       *
       * \param[in] space The underlying finite element space.
       * \param[in] coloring The coloring index set. Has to map each cell to its color. Obviously requires no two adjecent cells to have
       *                     the same color.
       * \param[in] hint Hint for the number of colors. If < 0 the mesh_data will determine this number itself. Undefined behavior if
       *                 number is *not* correct.
       *
       * //TODO: setting of base configuration does not fit completely, due to DataWrapperSize
       */
      template<typename ColoringType_>
      explicit VoxelBurgersVeloMaterialAssembler(const SpaceType& space, const ColoringType_& coloring, int hint = -1) :
      BaseClass(space, coloring, hint),
      frechet_material(DataType(0)), reg_eps(DataType(1E-100)), material_type(MaterialType::carreau)
      {
      }

      // rule of 5
      VoxelBurgersVeloMaterialAssembler(const VoxelBurgersVeloMaterialAssembler&) = delete;

      VoxelBurgersVeloMaterialAssembler& operator=(const VoxelBurgersVeloMaterialAssembler&) = delete;

      VoxelBurgersVeloMaterialAssembler(VoxelBurgersVeloMaterialAssembler&&) = default;

      VoxelBurgersVeloMaterialAssembler& operator=(VoxelBurgersVeloMaterialAssembler&&) = default;

      ~VoxelBurgersVeloMaterialAssembler(){}

      /**
       * \brief Wraps the set burgers parameter into convenient struct.
       *
       * \returns A wrapper for the burgers data.
       */
      VoxelAssembly::AssemblyMaterialData<DataType> wrap_material_params() const
      {
        return VoxelAssembly::AssemblyMaterialData<DataType>{frechet_material, reg_eps, mu_0, exp, lambda, a_T, yasuda_a, mu_inf, frechet_material>Math::pow(Math::eps<DataType>(), DataType(0.9))};
      }

      /**
       * \brief Assembles the burgers velo material operator for finitie element space with same test and trial space.
       *
       * \tparam CubatureFactory_ Factory type for cubature data.
       *
       * \param[in/out] matrix The matrix to be assembled. Adds onto residing data.
       * \param[in] convect The convection vector. Values are ignored if either beta or sd_delta is zero.
       * \param[in] space The underlying space.
       * \param[in] cubature_factory The cubature rule to be used.
       * \param[in] alpha Scaling parameter for assembly.
       */
      template<typename CubatureFactory_>
      void assemble_matrix1(LAFEM::SparseMatrixBCSR<DataType, IndexType, dim, dim>& matrix, const LAFEM::DenseVectorBlocked<DataType, IndexType, dim>& convect,
       const SpaceType& space, const CubatureFactory_& cubature_factory, DataType alpha = DataType(1)) const
      {
        XASSERTM(matrix.rows() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.columns() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");

        //define cubature
        typedef Cubature::Rule<ShapeType, DataType, DataType> CubatureRuleType;
        CubatureRuleType cubature(Cubature::ctor_factory, cubature_factory);

        //get cubature points and weights
        int num_cubs = cubature.get_num_points();
        typename CubatureRuleType::PointType* cub_pt = cubature.get_points();
        DataType* cub_wg = cubature.get_weights();

        BACKEND_SKELETON_VOID(assemble_matrix1_cuda, assemble_matrix1_generic, assemble_matrix1_generic, matrix, convect, space, cub_pt, cub_wg, num_cubs, alpha)
      }

      /**
       * \brief Assembles the application of the burgers operator to a given primal vector.
       *
       * \tparam CubatureFactory_ Factory type for cubature data.
       *
       * \param[in/out] vector The vector to be assembled. Adds onto residing data.
       * \param[in] convect The convection vector. Values are ignored if either beta or sd_delta is zero.
       * \param[in] primal The primal vector the operator is applied on.
       * \param[in] space The underlying space.
       * \param[in] cubature_factory The cubature rule to be used.
       * \param[in] alpha Scaling parameter for assembly.
       */
      template<typename CubatureFactory_>
      void assemble_vector(LAFEM::DenseVectorBlocked<DataType, IndexType, dim>& vector, const LAFEM::DenseVectorBlocked<DataType, IndexType, dim>& convect,
       const LAFEM::DenseVectorBlocked<DataType, IndexType, dim>& primal, const SpaceType& space, const CubatureFactory_& cubature_factory, DataType alpha = DataType(1)) const
      {
        XASSERTM(vector.size() == space.get_num_dofs(), "invalid vector size");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");

        //define cubature
        typedef Cubature::Rule<ShapeType, DataType, DataType> CubatureRuleType;
        CubatureRuleType cubature(Cubature::ctor_factory, cubature_factory);

        //get cubature points and weights
        int num_cubs = cubature.get_num_points();
        typename CubatureRuleType::PointType* cub_pt = cubature.get_points();
        DataType* cub_wg = cubature.get_weights();

        BACKEND_SKELETON_VOID(assemble_vector_cuda, assemble_vector_generic, assemble_vector_generic, vector, convect, primal, space, cub_pt, cub_wg, num_cubs, alpha)
      }

      void assemble_matrix1_generic(LAFEM::SparseMatrixBCSR<DataType, IndexType, dim, dim>& matrix, const LAFEM::DenseVectorBlocked<DataType, IndexType, dim>& convect,
       const SpaceType& space, typename Cubature::Rule<ShapeType, DataType, DataType>::PointType* cub_pt, const DataType* cub_wg, int num_cubs, DataType alpha) const
      {
        VoxelAssembly::CSRMatrixData<DataType, IndexType> mat_data = {matrix.template val<LAFEM::Perspective::pod>(), matrix.row_ptr(), matrix.col_ind(), matrix.rows(), matrix.columns()};
        const DataType* vec_data = convect.template elements<LAFEM::Perspective::pod>();

        VoxelAssembly::AssemblyCubatureData<DataType> cub_data = {(void*)cub_pt, cub_wg, num_cubs};
        VoxelAssembly::AssemblyMappingData<DataType, IndexType> mapping_data = this->mesh_data.get_assembly_field();
        auto burgers_params = this->wrap_burgers_params();
        auto material_params = this->wrap_material_params();


        VoxelAssembly::Arch::assemble_burgers_velo_material_csr_host(space, mat_data, vec_data, cub_data, mapping_data,
            this->mesh_data.get_coloring_maps(), this->mesh_data.get_color_map_sizes(), alpha,
                                  burgers_params, material_params, material_type);
      }

      void assemble_vector_generic(LAFEM::DenseVectorBlocked<DataType, IndexType, dim>& vector, const LAFEM::DenseVectorBlocked<DataType, IndexType, dim>& convect, const LAFEM::DenseVectorBlocked<DataType, IndexType, dim>& primal,
       const SpaceType& space, typename Cubature::Rule<ShapeType, DataType, DataType>::PointType* cub_pt, const DataType* cub_wg, int num_cubs, DataType alpha) const
      {
        DataType* vec_data = vector.template elements<LAFEM::Perspective::pod>();
        const DataType* conv_data = convect.template elements<LAFEM::Perspective::pod>();
        const DataType* primal_data = primal.template elements<LAFEM::Perspective::pod>();

        VoxelAssembly::AssemblyCubatureData<DataType> cub_data = {(void*)cub_pt, cub_wg, num_cubs};
        VoxelAssembly::AssemblyMappingData<DataType, IndexType> mapping_data = this->mesh_data.get_assembly_field();
        auto burgers_params = this->wrap_burgers_params();
        auto material_params = wrap_material_params();


        VoxelAssembly::Arch::assemble_burgers_velo_material_defect_host(space, vec_data, conv_data, primal_data, cub_data, mapping_data,
                                  this->mesh_data.get_coloring_maps(), this->mesh_data.get_color_map_sizes(), alpha,
                                  burgers_params, material_params, material_type);
        //free resources
      }

      #ifdef FEAT_HAVE_CUDA
      void assemble_matrix1_cuda(LAFEM::SparseMatrixBCSR<DataType, IndexType, dim, dim>& matrix, const LAFEM::DenseVectorBlocked<DataType, IndexType, dim>& convect,
       const SpaceType& space, typename Cubature::Rule<ShapeType, DataType, DataType>::PointType* cub_pt, const DataType* cub_wg, int num_cubs, DataType alpha) const
      {
        VoxelAssembly::CSRMatrixData<DataType, IndexType> mat_data = {matrix.template val<LAFEM::Perspective::pod>(), matrix.row_ptr(), matrix.col_ind(), matrix.rows(), matrix.columns()};
        const DataType* vec_data = convect.template elements<LAFEM::Perspective::pod>();

        typedef typename Cubature::Rule<ShapeType, DataType, DataType>::PointType CubPointType;
        //initialize all necessary pointer arrays and values //maybe more sense to specify cubature rule and set this to a const mem location?
        void* cub_pt_device = Util::cuda_malloc(Index(num_cubs) * sizeof(CubPointType));
        Util::cuda_copy_host_to_device(cub_pt_device, (void*)cub_pt, Index(num_cubs) * sizeof(CubPointType));

        void* cub_wg_device = Util::cuda_malloc(Index(num_cubs) * sizeof(DataType));
        Util::cuda_copy_host_to_device(cub_wg_device, (void*)cub_wg, Index(num_cubs) * sizeof(DataType));

        VoxelAssembly::AssemblyCubatureData<DataType> d_cub_data = {cub_pt_device, (DataType*)cub_wg_device, num_cubs};
        VoxelAssembly::AssemblyMappingData<DataType, IndexType> d_mapping_data = this->mesh_data.get_assembly_field();
        auto burgers_params = this->wrap_burgers_params();
        auto material_params = wrap_material_params();

        VoxelAssembly::Arch::assemble_burgers_velo_material_csr(space, mat_data, vec_data, d_cub_data, d_mapping_data,
                  this->mesh_data.get_coloring_maps(), this->mesh_data.get_color_map_sizes(), alpha, burgers_params, material_params, material_type,
                  this->shared_mem, this->blocksize, this->gridsize, this->print_occupancy);
        Util::cuda_free(cub_wg_device);
        Util::cuda_free(cub_pt_device);
      }


      void assemble_vector_cuda(LAFEM::DenseVectorBlocked<DataType, IndexType, dim>& vector, const LAFEM::DenseVectorBlocked<DataType, IndexType, dim>& convect, const LAFEM::DenseVectorBlocked<DataType, IndexType, dim>& primal,
       const SpaceType& space, typename Cubature::Rule<ShapeType, DataType, DataType>::PointType* cub_pt, const DataType* cub_wg, int num_cubs, DataType alpha) const
      {
        DataType* vec_data = vector.template elements<LAFEM::Perspective::pod>();
        const DataType* conv_data = convect.template elements<LAFEM::Perspective::pod>();
        const DataType* primal_data = primal.template elements<LAFEM::Perspective::pod>();

        typedef typename Cubature::Rule<ShapeType, DataType, DataType>::PointType CubPointType;
        //initialize all necessary pointer arrays and values //maybe more sense to specify cubature rule and set this to a const mem location?
        void* cub_pt_device = Util::cuda_malloc(Index(num_cubs) * sizeof(CubPointType));
        Util::cuda_copy_host_to_device(cub_pt_device, (void*)cub_pt, Index(num_cubs) * sizeof(CubPointType));

        void* cub_wg_device = Util::cuda_malloc(Index(num_cubs) * sizeof(DataType));
        Util::cuda_copy_host_to_device(cub_wg_device, (void*)cub_wg, Index(num_cubs) * sizeof(DataType));

        VoxelAssembly::AssemblyCubatureData<DataType> d_cub_data = {cub_pt_device, (DataType*)cub_wg_device, num_cubs};
        VoxelAssembly::AssemblyMappingData<DataType, IndexType> d_mapping_data = this->mesh_data.get_assembly_field();
        auto burgers_params = this->wrap_burgers_params();
        auto material_params = wrap_material_params();

        VoxelAssembly::Arch::assemble_burgers_velo_material_defect(space, vec_data, conv_data, primal_data, d_cub_data, d_mapping_data,
                  this->mesh_data.get_coloring_maps(), this->mesh_data.get_color_map_sizes(), alpha, burgers_params, material_params, material_type,
                  this->shared_mem, this->blocksize, this->gridsize, this->print_occupancy);
        Util::cuda_free(cub_wg_device);
        Util::cuda_free(cub_pt_device);
      }

      #else //FEAT_HAVE_CUDA
      void assemble_matrix1_cuda(LAFEM::SparseMatrixBCSR<DataType, IndexType, dim, dim>& DOXY(matrix), const SpaceType& DOXY(space), typename Cubature::Rule<ShapeType, DataType, DataType>::PointType* DOXY(cub_pt), const DataType* DOXY(cub_wg), int DOXY(num_cubs), DataType DOXY(alpha)) const
      {
        XABORTM("What in the nine hells are you doing here?");
      }

      void assemble_vector_cuda(LAFEM::DenseVectorBlocked<DataType, IndexType, dim>&, const LAFEM::DenseVectorBlocked<DataType, IndexType, dim>&, const LAFEM::DenseVectorBlocked<DataType, IndexType, dim>&,
       const SpaceType&, typename Cubature::Rule<ShapeType, DataType, DataType>::PointType*, const DataType*, int, DataType) const
      {
        XABORTM("This call was logged and your local FBI agent was informed about your transgression");
      }

      #endif

    }; // class VoxelBurgersAssembler<Q2StandardCube>
#endif
  }
}

#endif