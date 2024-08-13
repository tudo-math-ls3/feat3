#pragma once
#ifndef FEAT_KERNEL_VOXEL_ASSEMBLY_BURGERS_HPP
#define FEAT_KERNEL_VOXEL_ASSEMBLY_BURGERS_HPP 1

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

#ifdef FEAT_HAVE_CUDA
#include <kernel/util/cuda_util.hpp>
#endif

namespace FEAT
{
  namespace VoxelAssembly
  {
    /**
     * \brief Burgers Voxel Assembly template.
     *
     * Has to be specialized for each specific FE element space.
     *
     * \tparam Space_ The FE space for which it should be assembled
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indextype to be used.
     */
    template<typename Space_, typename DT_, typename IT_>
    class VoxelBurgersAssembler DOXY({});

    namespace Kernel
    {

      /**
       * \brief Burgers Matrix assembly kernel
       *
       * Outsources the the most inner execution of the assembly of the full burgers term for a cell previously gathered.
       * Local always refers to the local dofs of the current cell.
       *
       * \tparam SpaceHelp_ The psacehelper to be used. Holds most information about the underlying space.
       * \tparam LocMatType_ The Matrixtype of the local matrix dofs.
       * \tparam LocVecType_ The Vectortype of the local convection dofs.
       * \tparam dim_ Dimension of our space. For ease of access.
       * \tparam num_verts_ The number of vertices of one cell.
       *
       * \param[out] loc_mat A reference to the local matrix to be assembled.
       * \param[in] local_conv_dofs A reference to the local convection dofs. If convection is not needed values are ignored.
       * \param[in] local_coeffs The local trafo coefficients.
       * \param[in] cub_pt A pointer to the cubature point coordinate vector.
       * \param[in] cub_wg A pointer to the cubature weights. Same order as cub_pt array.
       * \param[in] burgers_params A struct holding the burgers parameter configuration.
       * \param[in] need_streamline Do we need streamline?
       * \param[in] need_convection Do we need convection?
       * \param[in] tol_eps Tolerance for local stream norm to be regarded as zero.
       */
      template<typename SpaceHelp_, typename LocMatType_, typename LocVecType_, int dim_, int num_verts_>
      CUDA_HOST_DEVICE void burgers_mat_assembly_kernel(LocMatType_& loc_mat, const LocVecType_& local_conv_dofs, const typename SpaceHelp_::DataType (&local_coeffs)[dim_][num_verts_],
                                          const typename SpaceHelp_::DomainPointType* cub_pt, const typename SpaceHelp_::DataType* cub_wg, const int num_cubs,
                                          const VoxelAssembly::AssemblyBurgersData<typename SpaceHelp_::DataType>& burgers_params,
                                          const bool need_streamline, const bool need_convection, const typename SpaceHelp_::DataType tol_eps)
      {
        typedef SpaceHelp_ SpaceHelp;
        constexpr int dim = SpaceHelp::dim;
        typedef typename SpaceHelp::SpaceType SpaceType;
        typedef typename SpaceHelp::DataType DataType;

        const DataType& nu{burgers_params.nu};
        const DataType& theta{burgers_params.theta};
        const DataType& beta{burgers_params.beta};
        const DataType& frechet_beta{burgers_params.frechet_beta};
        const DataType& sd_delta{burgers_params.sd_delta};
        const DataType& sd_nu{burgers_params.sd_nu};
        const DataType& sd_v_norm{burgers_params.sd_v_norm};
        const bool& deformation{burgers_params.deformation};

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

        VecValueType loc_v(DataType(0)), mean_v(DataType(0));
        MatValueType loc_grad_v(DataType(0));
        DataType local_delta(DataType(0));

        loc_mat.format();

        Tiny::Vector<DataType, num_loc_dofs> streamdiff_coeffs(DataType(0));

        if(need_streamline) //need streamdiff?
        {
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
            //we need this trafo data
            SpaceHelp::calc_jac_mat(loc_jac, barycenter, local_coeffs);
            loc_jac_inv.set_inverse(loc_jac);

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
          #ifdef __CUDACC__
          if(CudaMath::cuda_abs(frechet_beta) > DataType(0)) // need frechet beta?
          #else
          if(Math::abs(frechet_beta) > DataType(0)) // need frechet beta?
          #endif
          {
            loc_grad_v.format();
            for(int i = 0; i < num_loc_dofs; ++i)
            {
              loc_grad_v.add_outer_product(local_conv_dofs[i], basis_data.phi[i].grad);
            }
          }

          if(need_streamline) // need streamdiff?
          {
            for(int i = 0; i < num_loc_dofs; ++i)
            {
              streamdiff_coeffs[i] = Tiny::dot(loc_v, basis_data.phi[i].grad);
            }
          }

          if(deformation)
          {
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // trial function loop
              for(int j(0); j < num_loc_dofs; ++j)
              {
                // compute scalar value
                const DataType value = nu * weight * Tiny::dot(basis_data.phi[i].grad, basis_data.phi[j].grad);

                // update local matrix
                loc_mat[i][j].add_scalar_main_diag(value);
                loc_mat[i][j].add_outer_product(basis_data.phi[j].grad, basis_data.phi[i].grad, nu * weight);
              }
            }
          }
          else
          {
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // trial function loop
              for(int j(0); j < num_loc_dofs; ++j)
              {
                // compute scalar value
                const DataType value = nu * weight * Tiny::dot(basis_data.phi[i].grad, basis_data.phi[j].grad);

                // update local matrix
                loc_mat[i][j].add_scalar_main_diag(value);
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

          // assemble streamline diffusion?
          if((need_streamline) && (local_delta > tol_eps))
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
        // next cubature point
        }
      }

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
       *
       * \param[out] loc_vec A reference to the local defect to be assembled.
       * \param[in] local_prim_dofs A reference to the local primal values the matrix is implicitly applied to.
       * \param[in] local_conv_dofs A reference to the local convection dofs. If convection is not needed values are ignored.
       * \param[in] local_coeffs The local trafo coefficients.
       * \param[in] cub_pt A pointer to the cubature point coordinate vector.
       * \param[in] cub_wg A pointer to the cubature weights. Same order as cub_pt array.
       * \param[in] burgers_params A struct holding the burgers parameter configuration.
       * \param[in] need_convection Do we need convection?
       */
      template<typename SpaceHelp_, typename LocVecType_, int dim_, int num_verts_>
      CUDA_HOST_DEVICE void burgers_defect_assembly_kernel(LocVecType_& loc_vec, const LocVecType_& local_prim_dofs, const LocVecType_& local_conv_dofs, const typename SpaceHelp_::DataType (&local_coeffs)[dim_][num_verts_],
                                          const typename SpaceHelp_::DomainPointType* cub_pt, const typename SpaceHelp_::DataType* cub_wg, const int num_cubs,
                                          const VoxelAssembly::AssemblyBurgersData<typename SpaceHelp_::DataType>& burgers_params,
                                          const bool need_convection)
      {
        typedef SpaceHelp_ SpaceHelp;
        constexpr int dim = SpaceHelp::dim;
        typedef typename SpaceHelp::SpaceType SpaceType;
        typedef typename SpaceHelp::DataType DataType;

        const DataType& nu{burgers_params.nu};
        const DataType& theta{burgers_params.theta};
        const DataType& beta{burgers_params.beta};
        // const DataType& frechet_beta{burgers_params.frechet_beta};
        // const DataType& sd_delta{burgers_params.sd_delta};
        // const DataType& sd_nu{burgers_params.sd_nu};
        // const DataType& sd_v_norm{burgers_params.sd_v_norm};
        const bool& deformation{burgers_params.deformation};

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
        // typedef Tiny::Matrix<DataType, dim, dim> MatValueType;

        VecValueType loc_v(DataType(0));

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

          if(deformation)
          {
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // trial function loop
              for(int j(0); j < num_loc_dofs; ++j)
              {
                // compute scalar values
                const DataType value1 = nu * weight * Tiny::dot(basis_data.phi[i].grad, basis_data.phi[j].grad);
                const DataType value2 = nu * weight * Tiny::dot(local_prim_dofs[j], basis_data.phi[i].grad);
                // update local vector
                loc_vec[i].axpy(value1, local_prim_dofs[j]);
                loc_vec[i].axpy(value2, basis_data.phi[j].grad);
              }
            }
          }
          else
          {
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // trial function loop
              for(int j(0); j < num_loc_dofs; ++j)
              {
                // compute scalar value
                const DataType value = nu * weight * Tiny::dot(basis_data.phi[i].grad, basis_data.phi[j].grad);

                // update local matrix
                loc_vec[i].axpy(value, local_prim_dofs[j]);
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
       */
      template<typename Space_, typename DT_, typename IT_>
      void assemble_burgers_csr(const Space_& space,
              const CSRMatrixData<DT_, IT_>& matrix_data,
              const DT_* conv_data,
              const AssemblyCubatureData<DT_>& cubature,
              const AssemblyMappingData<DT_, IT_>& dof_mapping,
              const std::vector<int*>& coloring_maps,
              const std::vector<Index>& coloring_map_sizes,
              DT_ alpha, const AssemblyBurgersData<DT_>& burgers_params);

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
       */
      template<typename Space_, typename DT_, typename IT_>
      void assemble_burgers_defect(const Space_& space,
              DT_* vector_data,
              const DT_* conv_data,
              const DT_* primal_data,
              const AssemblyCubatureData<DT_>& cubature,
              const AssemblyMappingData<DT_, IT_>& dof_mapping,
              const std::vector<int*>& coloring_maps,
              const std::vector<Index>& coloring_map_sizes,
              DT_ alpha, const AssemblyBurgersData<DT_>& burgers_params);

      /**
       * \brief Device kernel wrapper for the local sd_v_norm calculation.
       *
       * Calculates the maximum euclidean norm of a given convection vector.
       *
       * \tparam DT_ The datatype.
       * \tparam IT_ The indextype.
       * \tparam dim_ Dimension of the convection vector.
       *
       * \param[in] convect The convection vector which max euclidean norm should be determined.
       *
       * \returns The maximum euclidean norm of all entries of convect.
       */
      template<typename DT_, typename IT_, int dim_>
      DT_ get_sd_v_norm(const LAFEM::DenseVectorBlocked<DT_, IT_, dim_>& convect);

      /**
       * \brief Device kernel wrapper for the global sd_v_norm calculation.
       *
       * Calculates the maximum euclidean norm of a given global convection vector.
       *
       * \tparam DT_ The datatype.
       * \tparam IT_ The indextype.
       * \tparam dim_ Dimension of the convection vector.
       * \param[in] convect The global convection vector which max euclidean norm should be determined.
       *
       * \returns The maximum euclidean norm of all entries of convect over all ranks.
       */
      template<typename DT_, typename IT_, int dim_>
      DT_ get_sd_v_norm(const Global::Vector<LAFEM::DenseVectorBlocked<DT_, IT_, dim_>, LAFEM::VectorMirror<DT_, IT_>>& convect);
      #endif

      /**
       * \brief Host kernel wrapper for the full matrix burgers assembler.
       *
       * Assembles the standard burgers operator on the host side in an additive manner.
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
       */
      template<typename Space_, typename DT_, typename IT_>
      void assemble_burgers_csr_host(const Space_& space,
              const CSRMatrixData<DT_, IT_>& matrix_data,
              const DT_* conv_data,
              const AssemblyCubatureData<DT_>& cubature,
              const AssemblyMappingData<DT_, IT_>& dof_mapping,
              const std::vector<int*>& coloring_maps,
              const std::vector<Index>& coloring_map_sizes,
              DT_ alpha, const AssemblyBurgersData<DT_>& burgers_params);

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
       */
      template<typename Space_, typename DT_, typename IT_>
      void assemble_burgers_defect_host(const Space_& space,
              DT_* vector_data,
              const DT_* conv_data,
              const DT_* primal_data,
              const AssemblyCubatureData<DT_>& cubature,
              const AssemblyMappingData<DT_, IT_>& dof_mapping,
              const std::vector<int*>& coloring_maps,
              const std::vector<Index>& coloring_map_sizes,
              DT_ alpha, const AssemblyBurgersData<DT_>& burgers_params);

      /**
       * \brief Host kernel wrapper for the local sd_v_norm calculation.
       *
       * Calculates the maximum euclidean norm of a given convection vector.
       *
       * \tparam DT_ The datatype.
       * \tparam IT_ The indextype.
       * \tparam dim_ Dimension of the convection vector.
       *
       * \param[in] convect The convection vector which max euclidean norm should be determined.
       *
       * \returns The maximum euclidean norm of all entries of convect.
       */
      template<typename DT_, typename IT_, int dim_>
      DT_ get_sd_v_norm_host(const LAFEM::DenseVectorBlocked<DT_, IT_, dim_>& convect);

      /**
       * \brief Host kernel wrapper for the global sd_v_norm calculation.
       *
       * Calculates the maximum euclidean norm of a given global convection vector.
       *
       * \tparam DT_ The datatype.
       * \tparam IT_ The indextype.
       * \tparam dim_ Dimension of the convection vector.
       *
       * \param[in] convect The global convection vector which max euclidean norm should be determined.
       *
       * \returns The maximum euclidean norm of all entries of convect over all ranks.
       */
      template<typename DT_, typename IT_, int dim_>
      DT_ get_sd_v_norm_host(const Global::Vector<LAFEM::DenseVectorBlocked<DT_, IT_, dim_>, LAFEM::VectorMirror<DT_, IT_>>& convect);

    }


#ifndef __CUDACC__
    /**
     * \brief Q2Lagrange thread parallel assembly class for the burgers operator.
     *
     * This class provides host and device variants for a thread parallel assembly of the burgers operator.
     * For details on the burgers operator, see the burgers assembler in </kernel/assembly/>.
     * The threading strategy is based on a cell based coloring of the mesh, so in principle this version for Q2
     * space works as long as a valid coloring is provided.
     *
     * \tparam dim_ Dimension of the hypercubes.
     * \tparam DT_ Datatype to be used.
     * \tparam IT_ Indextype to be used.
     */
    template<int dim_, typename DT_, typename IT_>
    class VoxelBurgersAssembler<Q2StandardHyperCube<dim_>, DT_, IT_>
    {
    public:
      /// typedefs
      typedef VoxelAssembly::Q2StandardHyperCube<dim_> SpaceType;
      typedef LagrangeDataHandler<SpaceType, DT_, IT_> DataHandler;
      typedef VoxelAssembly::SpaceHelper<SpaceType, DT_, IT_> SpaceHelp;
      typedef typename SpaceHelp::ShapeType ShapeType;
      typedef typename SpaceHelp::DataType DataType;
      typedef typename SpaceHelp::IndexType IndexType;

      /// constexpr
      static constexpr int dim = SpaceHelp::dim;

      typedef typename SpaceHelp::DomainPointType DomainPointType;
      typedef typename SpaceHelp::ImagePointType ImagePointType;
      typedef typename SpaceHelp::ValueType ValueType;
      typedef typename SpaceHelp::JacobianMatrixType JacobianMatrixType;

      /// The meshdata handler
      DataHandler mesh_data;

      /// The viscosity
      DataType nu;

      /// Reaction scaling
      DataType theta;

      /// Convection scaling
      DataType beta;

      /// Scaling parameter for full jacobian convection part
      DataType frechet_beta;

      /// Streamline diffusion parameter.
      DataType sd_delta;

      /// Streamline diffosuion viscosity. Generally equal to nu.
      DataType sd_nu;

      /// The currecnt maximum velocity norm.
      DataType sd_v_norm;

      /// Should we use full deformation formulation?
      bool deformation;


    public:
      explicit VoxelBurgersAssembler() = default;

      /**
       * \brief Constructor for burgers assembler
       *
       * \tparam ColoringType_ The type of the coloring indexset.
       *
       * \param[in] space The underlying finite element space.
       * \param[in] coloring The coloring index set. Has to map each cell to its color. Obviously requires no two adjecent cells to have
       *                     the same color.
       * \param[in] hint Hint for the number of colors. If < 0 the mesh_data will determine this number itself. Undefined behavior if
       *                 number is *not* correct.
       */
      template<typename ColoringType_>
      explicit VoxelBurgersAssembler(const SpaceType& space, const ColoringType_& coloring, int hint = -1) :
      mesh_data(space, coloring, hint),
      nu(DataType(0)), theta(DataType(0)), beta(DataType(0)), frechet_beta(DataType(0)), sd_delta(DataType(0)),
      sd_nu(DataType(0)), sd_v_norm(DataType(0)), deformation(true)
      {
        #ifdef FEAT_HAVE_CUDA
        #ifdef DEBUG
        std::cout << "Setting cache per thread up to meet debug demands!" << std::endl;
        if constexpr(dim == 3)
          Util::cuda_set_max_cache_thread(8096u * sizeof(DataType));
        else
          Util::cuda_set_max_cache_thread(1012u * sizeof(DataType));
        #endif
        #endif
      }

      /**
       * \brief Wraps the set burgers parameter into convenient struct.
       *
       * \returns A wrapper for the burgers data.
       */
      VoxelAssembly::AssemblyBurgersData<DataType> wrap_burgers_params() const
      {
        return VoxelAssembly::AssemblyBurgersData<DataType>{nu, theta, beta, frechet_beta, sd_delta, sd_nu, sd_v_norm, deformation};
      }

      /**
       * \brief Assembles the burgers operator for finitie element space with same test and trial space.
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

      /**
       * \brief Sets the internal maximum streamline diffusion velocity norm.
       *
       * \param[in] convect The local convection vector which is evaluated regarding its maximum euclidean norm.
       */
      void set_sd_v_norm(const LAFEM::DenseVectorBlocked<DataType, IndexType, dim>& convect)
      {
        BACKEND_SKELETON_VOID(set_sd_v_norm_cuda, set_sd_v_norm_generic, set_sd_v_norm_generic, convect)
      }

      /**
       * \brief Sets the internal maximum streamline diffusion velocity norm.
       *
       * \param[in] convect The global convection vector which is evaluated regarding its maximum euclidean norm.
       */
      void set_sd_v_norm(const Global::Vector<LAFEM::DenseVectorBlocked<DataType, IndexType, dim>, LAFEM::VectorMirror<DataType, IndexType>>& convect)
      {
        BACKEND_SKELETON_VOID(set_sd_v_norm_cuda, set_sd_v_norm_generic, set_sd_v_norm_generic, convect)
      }

      void assemble_matrix1_generic(LAFEM::SparseMatrixBCSR<DataType, IndexType, dim, dim>& matrix, const LAFEM::DenseVectorBlocked<DataType, IndexType, dim>& convect,
       const SpaceType& space, typename Cubature::Rule<ShapeType, DataType, DataType>::PointType* cub_pt, const DataType* cub_wg, int num_cubs, DataType alpha) const
      {
        VoxelAssembly::CSRMatrixData<DataType, IndexType> mat_data = {matrix.template val<LAFEM::Perspective::pod>(), matrix.row_ptr(), matrix.col_ind(), matrix.rows(), matrix.columns()};
        const DataType* vec_data = convect.template elements<LAFEM::Perspective::pod>();

        VoxelAssembly::AssemblyCubatureData<DataType> cub_data = {(void*)cub_pt, cub_wg, num_cubs};
        VoxelAssembly::AssemblyMappingData<DataType, IndexType> mapping_data = mesh_data.get_assembly_field();
        auto burgers_params = wrap_burgers_params();


        VoxelAssembly::Arch::assemble_burgers_csr_host(space, mat_data, vec_data, cub_data, mapping_data, mesh_data.get_coloring_maps(), mesh_data.get_color_map_sizes(), alpha,
                                  burgers_params);
      }

      void assemble_vector_generic(LAFEM::DenseVectorBlocked<DataType, IndexType, dim>& vector, const LAFEM::DenseVectorBlocked<DataType, IndexType, dim>& convect, const LAFEM::DenseVectorBlocked<DataType, IndexType, dim>& primal,
       const SpaceType& space, typename Cubature::Rule<ShapeType, DataType, DataType>::PointType* cub_pt, const DataType* cub_wg, int num_cubs, DataType alpha) const
      {
        DataType* vec_data = vector.template elements<LAFEM::Perspective::pod>();
        const DataType* conv_data = convect.template elements<LAFEM::Perspective::pod>();
        const DataType* primal_data = primal.template elements<LAFEM::Perspective::pod>();

        VoxelAssembly::AssemblyCubatureData<DataType> cub_data = {(void*)cub_pt, cub_wg, num_cubs};
        VoxelAssembly::AssemblyMappingData<DataType, IndexType> mapping_data = mesh_data.get_assembly_field();
        auto burgers_params = wrap_burgers_params();


        VoxelAssembly::Arch::assemble_burgers_defect_host(space, vec_data, conv_data, primal_data, cub_data, mapping_data, mesh_data.get_coloring_maps(), mesh_data.get_color_map_sizes(), alpha,
                                  burgers_params);
        //free resources
      }

      void set_sd_v_norm_generic(const LAFEM::DenseVectorBlocked<DataType, IndexType, dim>& convect)
      {
        this->sd_v_norm = VoxelAssembly::Arch::get_sd_v_norm_host(convect);
      }


      void set_sd_v_norm_generic(const Global::Vector<LAFEM::DenseVectorBlocked<DataType, IndexType, dim>, LAFEM::VectorMirror<DataType, IndexType>>& convect)
      {
        this->sd_v_norm = VoxelAssembly::Arch::get_sd_v_norm_host(convect);
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
        VoxelAssembly::AssemblyMappingData<DataType, IndexType> d_mapping_data = mesh_data.get_assembly_field();
        auto burgers_params = wrap_burgers_params();

        VoxelAssembly::Arch::template assemble_burgers_csr(space, mat_data, vec_data, d_cub_data, d_mapping_data, mesh_data.get_coloring_maps(), mesh_data.get_color_map_sizes(), alpha, burgers_params);
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
        VoxelAssembly::AssemblyMappingData<DataType, IndexType> d_mapping_data = mesh_data.get_assembly_field();
        auto burgers_params = wrap_burgers_params();

        VoxelAssembly::Arch::template assemble_burgers_defect(space, vec_data, conv_data, primal_data, d_cub_data, d_mapping_data, mesh_data.get_coloring_maps(), mesh_data.get_color_map_sizes(), alpha, burgers_params);
        Util::cuda_free(cub_wg_device);
        Util::cuda_free(cub_pt_device);
      }

      void set_sd_v_norm_cuda(const LAFEM::DenseVectorBlocked<DataType, IndexType, dim>& convect)
      {
        this->sd_v_norm = VoxelAssembly::Arch::get_sd_v_norm(convect);
      }


      void set_sd_v_norm_cuda(const Global::Vector<LAFEM::DenseVectorBlocked<DataType, IndexType, dim>, LAFEM::VectorMirror<DataType, IndexType>>& convect)
      {
        this->sd_v_norm = VoxelAssembly::Arch::get_sd_v_norm(convect);
      }


      #else //FEAT_HAVE_CUDA
      void assemble_matrix1_cuda(LAFEM::SparseMatrixBCSR<DataType, IndexType, dim, dim>& DOXY(matrix), const SpaceType& DOXY(space), typename Cubature::Rule<ShapeType, DataType, DataType>::PointType* DOXY(cub_pt), const DataType* DOXY(cub_wg), int DOXY(num_cubs), DataType DOXY(alpha)) const
      {
        XABORTM("What in the nine hells are you doing here?");
      }

      void assemble_vector_cuda(LAFEM::DenseVectorBlocked<DataType, IndexType, dim>&, const LAFEM::DenseVectorBlocked<DataType, IndexType, dim>&, const LAFEM::DenseVectorBlocked<DataType, IndexType, dim>&,
       const SpaceType&, typename Cubature::Rule<ShapeType, DataType, DataType>::PointType*, const DataType*, int, DataType) const
      {
        XABORTM("This call was logged and your local FBI agent was informed about your transgressions");
      }

      void set_sd_v_norm_cuda(const LAFEM::DenseVectorBlocked<DataType, IndexType, dim>&)
      {
        XABORTM("Should not have called that....");
      }


      void set_sd_v_norm_cuda(const Global::Vector<LAFEM::DenseVectorBlocked<DataType, IndexType, dim>, LAFEM::VectorMirror<DataType, IndexType>>&)
      {
        XABORTM("Should not have called that....");
      }

      #endif

    }; // class VoxelBurgersAssembler<Q2StandardCube>
#endif
  }
}

#endif