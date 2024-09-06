#pragma once
#ifndef KERNEL_ASSEMBLY_SPACE_HELPER_HPP
#define KERNEL_ASSEMBLY_SPACE_HELPER_HPP 1

// base FEAT includes
#include <kernel/base_header.hpp>
#include <kernel/shape.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/geometry/conformal_mesh.hpp>

// space includes
#include <kernel/space/eval_data.hpp>
#include <kernel/space/details.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/lagrange2/details.hpp>

namespace FEAT
{
  namespace VoxelAssembly
  {
    template<typename IndexType_>
    class IndexSetWrapper
    {
      public:
      typedef IndexType_ IndexType;
      const IndexType* cell_to_dofs;
      IndexType num_loc_dofs;

      CUDA_HOST_DEVICE IndexType operator()(IndexType cell_index, IndexType vert_index) const
      {
        return *(cell_to_dofs + cell_index*num_loc_dofs + vert_index);
      }

      CUDA_HOST_DEVICE IndexSetWrapper(const IndexType* _cell_t, IndexType numi)
      :
      cell_to_dofs(_cell_t),
      num_loc_dofs(numi)
      {}

      CUDA_HOST_DEVICE IndexSetWrapper(const IndexSetWrapper&) = delete;

    };//class IndexSetWrapper

    template<typename IndexType_>
    class SharedIndexSetWrapper
    {
      public:
      typedef IndexType_ IndexType;
      const IndexType* loc_cell_to_dofs;

      CUDA_HOST_DEVICE IndexType operator()(IndexType DOXY(cell_index), IndexType vert_index) const
      {
        return *(loc_cell_to_dofs + vert_index);
      }

      CUDA_HOST_DEVICE SharedIndexSetWrapper(const IndexType* _cell_t)
      :
      loc_cell_to_dofs(_cell_t)
      {}

      CUDA_HOST_DEVICE SharedIndexSetWrapper(const SharedIndexSetWrapper&) = delete;

    };

    template<typename VertexType_>
    class VertexSetWrapper
    {
      public:
      typedef VertexType_ VertexType;
      const VertexType* vert_array;

      CUDA_HOST_DEVICE VertexType operator[](Index index) const
      {
        return vert_array[index];
      }

      CUDA_HOST_DEVICE VertexSetWrapper(const VertexType* _verts)
      :
      vert_array(_verts)
      {}

      CUDA_HOST_DEVICE VertexSetWrapper(const VertexSetWrapper&) = delete;

    };//class VertexSetWrapper

    /// we only use standard mappings for voxel assembly
    template<typename Shape_>
    using Q2StandardFE = Space::Lagrange2::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape_>>>;

    template<typename Space_, typename DT_, typename IT_>
    struct SpaceHelper;

    template<typename Shape_, typename DT_, typename IT_>
    struct SpaceHelper<Q2StandardFE<Shape_>, DT_, IT_>
    {
      typedef Q2StandardFE<Shape_> SpaceType;
      typedef Shape_ ShapeType;
      typedef IT_ IndexType;

      /// constexpr
      static constexpr int dim = SpaceType::world_dim;

      /// The value type
      typedef DT_ ValueType;

      /// The domainpoint type
      typedef Tiny::Vector<DT_, dim> DomainPointType;

      /// The image point type
      typedef Tiny::Vector<DT_, dim> ImagePointType;

      /// The vertex point type
      typedef Tiny::Vector<DT_, dim> VertexPointType;

      /// jacobian matrix type
      typedef Tiny::Matrix<DT_, dim, dim> JacobianMatrixType;

      /// jacobian inverse matrix type
      typedef Tiny::Matrix<DT_, dim, dim> JacobianInverseType;

      /// jacobian determinant type
      typedef DT_ JacobianDeterminantType;

      /// hessian tensor type
      typedef Tiny::Tensor3<DT_, dim, dim, dim> HessianTensorType;

      /// hessian inverse tensor type
      typedef Tiny::Tensor3<DT_, dim, dim, dim> HessianInverseType;

      /// Parametric eval helper
      typedef Space::ParametricEvalHelper ParamEvalHelp;

      /// The space helper detail class
      typedef Space::Lagrange2::EvalHelper<DomainPointType, DT_, ShapeType> SpaceEvalHelp;

      /// Define an EvalTraits struct
      struct EvalTraits
      {
        typedef DT_ DataType;
        typedef DataType BasisValueType;
        typedef DataType BasisReferenceValueType;
        typedef Tiny::Vector<DataType, dim> BasisGradientType;
        typedef Tiny::Vector<DataType, dim> BasisReferenceGradientType;
        typedef Tiny::Matrix<DataType, dim, dim> BasisHessianType;
        typedef Tiny::Matrix<DataType, dim, dim> BasisReferenceHessianType;
        static constexpr int max_local_dofs = SpaceEvalHelp::get_num_local_dofs();
      };

      typedef typename EvalTraits::DataType DataType;

      /// The trafo helper detail class
      typedef Trafo::Standard::EvalHelper<DataType, DomainPointType, ImagePointType, JacobianMatrixType,
                                        JacobianInverseType, JacobianDeterminantType, HessianTensorType,
                                        HessianInverseType, ShapeType, dim> TrafoEvalHelp;

      static constexpr int num_verts = TrafoEvalHelp::num_verts;

      /// A dummy space tags
      static constexpr SpaceTags config = SpaceTags::none;

      /// The eval data wrapping the basis functions
      typedef Space::EvalDataReduced<EvalTraits, config> EvalData;
      static constexpr bool reduced_data = true;

      /// Ref value calculation
      CUDA_HOST_DEVICE static inline void eval_ref_values(EvalData& data, const DomainPointType& point)
      {
        SpaceEvalHelp::eval_ref_values(data, point);
      }

      #ifdef __CUDACC__
      /// Group Ref value calculation
      template<typename ThreadGroup_>
      CUDA_DEVICE static __forceinline__ void group_eval_ref_values(const ThreadGroup_& tg, EvalData& data, const DomainPointType& point)
      {
        SpaceEvalHelp::group_eval_ref_values(tg, data, point);
      }
      #endif

      /// Ref value calculation
      CUDA_HOST_DEVICE static inline void eval_ref_gradients(EvalData& data, const DomainPointType& point)
      {
        SpaceEvalHelp::eval_ref_gradients(data, point);
      }

      /// ref value calculation
      CUDA_HOST_DEVICE static inline void eval_ref_hessians(EvalData& data, const DomainPointType& point)
      {
        SpaceEvalHelp::eval_ref_hessians(data, point);
      }

      CUDA_HOST_DEVICE static inline void set_coefficients(Tiny::Matrix<DataType, dim, num_verts>& coeffs, const IndexSetWrapper<IndexType>& local_dofs, const VertexPointType* vertex_set, IndexType cell_index)
      {
        TrafoEvalHelp::set_coefficients(coeffs, VertexSetWrapper(vertex_set), local_dofs, cell_index);
      }

      #ifdef __CUDACC__
      template<typename ThreadGroup_>
      CUDA_DEVICE static inline void grouped_set_coefficients(const ThreadGroup_& tg, DataType* coeffs, const SharedIndexSetWrapper<IndexType>& local_dofs, const VertexPointType* vertex_set, IndexType cell_index)
      {
        TrafoEvalHelp::grouped_set_coefficients(tg, coeffs, VertexSetWrapper(vertex_set), local_dofs, cell_index);
      }
      #endif

      CUDA_HOST_DEVICE static inline void map_point(typename TrafoEvalHelp::ImagePointType& img_point, const typename TrafoEvalHelp::DomainPointType& dom_point, const Tiny::Matrix<DataType, dim, num_verts>& coeffs)
      {
        TrafoEvalHelp::map_point(img_point, dom_point, coeffs);
      }

      CUDA_HOST_DEVICE static inline void calc_jac_mat(typename TrafoEvalHelp::JacobianMatrixType& jac_mat, const typename TrafoEvalHelp::DomainPointType& dom_point, const Tiny::Matrix<DataType, dim, num_verts>& coeffs)
      {
        TrafoEvalHelp::calc_jac_mat(jac_mat, dom_point, coeffs);
      }

      #ifdef __CUDACC__
      template<typename ThreadGroup_>
      CUDA_DEVICE static inline void grouped_calc_jac_mat(const ThreadGroup_& tg, typename TrafoEvalHelp::JacobianMatrixType& jac_mat, const typename TrafoEvalHelp::DomainPointType& dom_point, const DataType* coeffs)
      {
        TrafoEvalHelp::grouped_calc_jac_mat(tg, jac_mat, dom_point, coeffs);
      }
      #endif

      CUDA_HOST_DEVICE static inline void calc_hess_ten(typename TrafoEvalHelp::HessianTensorType& hess_ten, const typename TrafoEvalHelp::DomainPointType& dom_point, const Tiny::Matrix<DataType, dim, num_verts>& coeffs)
      {
        TrafoEvalHelp::calc_jac_mat(hess_ten, dom_point, coeffs);
      }

      CUDA_HOST_DEVICE static inline DataType volume(const Tiny::Matrix<DataType, dim, num_verts>& coeffs)
      {
        return TrafoEvalHelp::volume(coeffs);
      }

      CUDA_HOST_DEVICE static inline DataType width_directed(const typename TrafoEvalHelp::ImagePointType& ray, const Tiny::Matrix<DataType, dim, num_verts>& coeffs)
      {
        return TrafoEvalHelp::width_directed(ray, coeffs);
      }

      CUDA_HOST_DEVICE static inline void trans_values(EvalData& data)
      {
        if constexpr(!reduced_data)
          ParamEvalHelp::trans_values(data);
      }

      CUDA_HOST_DEVICE static inline void trans_gradients(EvalData& data, const typename TrafoEvalHelp::JacobianInverseType& jac_inv)
      {
        // if constexpr(!reduced_data)
          ParamEvalHelp::trans_gradients(data, jac_inv);
        // else
          // ParamEvalHelp::inplace_trans_gradients(data, jac_inv);
      }

      CUDA_HOST_DEVICE static inline void trans_hessian(EvalData& data, const typename TrafoEvalHelp::JacobianInverseType& jac_inv, const typename TrafoEvalHelp::HessianInverseType& hess_inv)
      {
        ParamEvalHelp::trans_hessians(data, jac_inv, hess_inv);
      }

    };

  }

}


#endif