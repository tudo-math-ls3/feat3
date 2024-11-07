#pragma once

#include <kernel/base_header.hpp>
#include <kernel/shape.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange2/element.hpp>

#ifdef FEAT_HAVE_OMP
#include "omp.h"
#endif

#ifdef __CUDACC__
#include <cooperative_groups.h>
#include <cooperative_groups/memcpy_async.h>
#include <cooperative_groups/reduce.h>

namespace cg = cooperative_groups;
#endif

namespace FEAT
{
  /**
   * \brief Namespace for different voxel based assembly methods
   *
   * This namespace incoperates a number of function and classes handling
   * thread parallel assembly of finite element matrices.
   * Since the threading strategy is based on coloring of the elements,
   * a good (i.e. balanced) coloring is required, which is for example given
   * by a voxel based mesh.
   *
   * \author Maximilian Esser
   */
  namespace VoxelAssembly
  {
    /// Q2 quadliteral space
    typedef Space::Lagrange2::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2>>>> Q2StandardQuad;
    /// Q2 hexadron space
    typedef Space::Lagrange2::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<3>>>> Q2StandardHexa;

    /// Templated hypercube definition for standard mapping
    template<int dim_>
    using Q2StandardHyperCube = Space::Lagrange2::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<dim_>>>>;


    /**
     * \brief A data field for all necessary values that define the dof mapping for assembly.
     *
     * \tparam DT_ The datatype.
     * \tparam IT_ The indextype.
     */
    template<typename DT_, typename IT_>
    struct AssemblyMappingData
    {
      /// The cell to dof, where cell_to_dof[i],..., cell_to_dof[i+cell_dofs-1] are the dofs of one cell
      const IT_* cell_to_dof;
      /// Array of sortingindices of cell_to_dof
      const IT_* cell_to_dof_sorter;
      /// The number of cells
      Index cell_num;
      /// An array of the nodes fitting to the cell_to_dof mapping
      const void* nodes;
      /// The number of nodes
      Index node_size;
    };

    /**
     * \brief A data field for a cubature rule
     *
     * \tparam DT_ The dataype.
     */
    template<typename DT_>
    struct AssemblyCubatureData
    {
      /// The cubature point data array
      const void* cub_pt;
      /// The cubature weights
      const DT_* cub_wg;
      /// Number of cubtaure points
      int num_cubs;
    };

    /**
     * \brief CSR Matrix data
     *
     * \tparam DT_ The dataype.
     * \tparam IT_ The indextype.
     */
    template<typename DT_, typename IT_>
    struct CSRMatrixData
    {
      DT_* data;
      const IT_* row_ptr;
      const IT_* col_idx;
      Index num_rows;
      Index num_cols;
    };

    /**
     * \brief Data for burgers assembler
     *
     * \tparam DT_ The datatype.
     */
    template<typename DT_>
    struct AssemblyBurgersData
    {
      DT_ nu;
      DT_ theta;
      DT_ beta;
      DT_ frechet_beta;
      DT_ sd_delta;
      DT_ sd_nu;
      DT_ sd_v_norm;
      bool deformation;
    };

    template<typename SpaceHelp_>
    struct BurgersSharedDataGlobalWrapper
    {
      typedef SpaceHelp_ SpaceHelp;
      typedef typename SpaceHelp::SpaceType SpaceType;
      typedef typename SpaceHelp::DataType DataType;
      typedef typename SpaceHelp::IndexType IndexType;
      static constexpr int dim = SpaceHelp::dim;
      static constexpr int num_loc_dofs = SpaceType::DofMappingType::dof_count;
      static constexpr int num_loc_verts = SpaceHelp::num_verts;
      DataType local_coeffs[dim*num_loc_verts];
      DataType local_conv_dofs[dim*num_loc_dofs];
      IndexType local_dofs[num_loc_dofs];
      // these can be defined in shared memory
      DataType tol_eps;
      bool need_streamline;
      bool need_convection;
    };

    template<typename SpaceHelp_>
    struct BurgersDefectSharedDataGlobalWrapper
    {
      typedef SpaceHelp_ SpaceHelp;
      typedef typename SpaceHelp::SpaceType SpaceType;
      typedef typename SpaceHelp::DataType DataType;
      typedef typename SpaceHelp::IndexType IndexType;
      static constexpr int dim = SpaceHelp::dim;
      static constexpr int num_loc_dofs = SpaceType::DofMappingType::dof_count;
      static constexpr int num_loc_verts = SpaceHelp::num_verts;
      DataType local_coeffs[dim*num_loc_verts];
      DataType local_conv_dofs[dim*num_loc_dofs];
      DataType local_prim_dofs[dim*num_loc_dofs];
      IndexType local_dofs[num_loc_dofs];
      // these can be defined in shared memory
      bool need_convection;
    };

    #ifdef __CUDACC__
    template<typename DT_, int dim_>
    CUDA_DEVICE inline DT_ dot(const Tiny::Vector<DT_, dim_>& grad, const DT_* conv, DT_ alpha = DT_(1))
    {
      DT_ lv(DT_(0));
      for(int k = 0; k < dim_; ++k)
      {
        lv += alpha * grad[k] * conv[k];
      }
      return lv;
    }

    template<typename DT_, typename TG_>
    CUDA_DEVICE inline void coalesced_format(const TG_& thread_group, DT_* target, int size)
    {
      for(int i = thread_group.thread_rank(); i < size; i += thread_group.num_threads())
      {
        target[i] = DT_(0);
      }
    }

    template<typename DT_, typename TG_, int dim_>
    CUDA_DEVICE inline void grouped_special_dot(const TG_& tg, DT_& loc_val, const DT_ phi, const DT_* conv, DT_ alpha = DT_(1))
    {
      auto coal_g = cg::coalesced_threads();
      DT_ lv(DT_(0));
      for(int k = 0; k < dim_; ++k)
      {
        lv += alpha * phi * conv[k];
      }
      lv = cg::reduce(coal_g, lv, cg::plus<DT_>());
      cg::invoke_one(coal_g, [&](){atomicAdd(&loc_val, lv);});
    }

    template<typename DT_, typename TG_, int dim_>
    CUDA_DEVICE inline void grouped_dot(const TG_& tg, DT_& loc_val, const Tiny::Vector<DT_, dim_>& grad, const DT_* conv, DT_ alpha = DT_(1))
    {
      auto coal_g = cg::coalesced_threads();
      DT_ lv(DT_(0));
      for(int k = 0; k < dim_; ++k)
      {
        lv += alpha * grad[k] * conv[k];
      }
      lv = cg::reduce(coal_g, lv, cg::plus<DT_>());
      cg::invoke_one(coal_g, [&](){atomicAdd(&loc_val, lv);});
    }

    template<typename DT_, typename TG_, int dim_>
    CUDA_DEVICE inline void grouped_axpy(const TG_& tg, DT_* loc_v, const DT_ phi, const DT_* conv, DT_ alpha = DT_(1))
    {
      auto coal_g = cg::coalesced_threads();
      for(int k = 0; k < dim_; ++k)
      {
        //thos could also be done a bit nicer by using async reduce versions...
        DT_ lv(DT_(0));
        lv = cg::reduce(coal_g, alpha * phi * conv[k], cg::plus<DT_>());
        cg::invoke_one(coal_g, [&](){atomicAdd(loc_v+k, lv);});
        // coal_g.snyc();
      }
    }

    /// Note reuqired to be symmetric!
    template<typename DT_, typename TG_, int dim_>
    CUDA_DEVICE inline void grouped_add_outer_product(const TG_& tg, DT_* loc_grad_v, const DT_* conv, const Tiny::Vector<DT_, dim_>& grad, DT_ alpha = DT_(1))
    {
      auto coal_g = cg::coalesced_threads();
      for(int i(0); i < dim_; ++i)
      {
        for(int j(0); j < dim_; ++j)
        {
          DT_ lv = cg::reduce(coal_g, alpha * conv[i] * grad[j], cg::plus<DT_>());
          cg::invoke_one(coal_g, [&](){atomicAdd(&loc_grad_v[i*dim_ + j], lv);});
        }
      }
    }

    template<typename DT_>
    CUDA_DEVICE DT_* get_aligned_address(void* ptr)
    {
      // alignment is always a power of 2, so can use the following bit magic
      // first add the (alignment-1) of the datatype, which leads to the next aligned
      // address modulo bits that fall into the alignment.
      // The simply cut off everything from the address larger then our aligned address by using
      // AND with the negative alignment, which is simply all bits flipped of (alignment-1)
      constexpr std::uint64_t align = alignof(DT_);
      std::uint64_t ptr_a = reinterpret_cast<std::uint64_t>(ptr);
      return (DT_*)((ptr_a + align - 1u)&(-align));
    }

    #endif

  }
}
