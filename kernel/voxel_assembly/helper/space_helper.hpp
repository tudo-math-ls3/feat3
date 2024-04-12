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
      typedef Space::EvalData<EvalTraits, config> EvalData;

      /// Ref value calculation
      CUDA_FUNC static inline void eval_ref_values(EvalData& data, const DomainPointType& point)
      {
        SpaceEvalHelp::eval_ref_values(data, point);
      }

      /// Ref value calculation
      CUDA_FUNC static inline void eval_ref_gradients(EvalData& data, const DomainPointType& point)
      {
        SpaceEvalHelp::eval_ref_gradients(data, point);
      }

      /// ref value calculation
      CUDA_FUNC static inline void eval_ref_hessians(EvalData& data, const DomainPointType& point)
      {
        SpaceEvalHelp::eval_ref_hessians(data, point);
      }

      CUDA_FUNC static inline void set_coefficients(DataType (&coeffs)[dim][num_verts], const IndexType* local_dofs, const VertexPointType* vertex_set)
      {
        TrafoEvalHelp::set_coefficients(coeffs, local_dofs, vertex_set);
      }

      CUDA_FUNC static inline void map_point(typename TrafoEvalHelp::ImagePointType& img_point, const typename TrafoEvalHelp::DomainPointType& dom_point, const DataType (&coeffs)[dim][num_verts])
      {
        TrafoEvalHelp::map_point(img_point, dom_point, coeffs);
      }

      CUDA_FUNC static inline void calc_jac_mat(typename TrafoEvalHelp::JacobianMatrixType& jac_mat, const typename TrafoEvalHelp::DomainPointType& dom_point, const DataType (&coeffs)[dim][num_verts])
      {
        TrafoEvalHelp::calc_jac_mat(jac_mat, dom_point, coeffs);
      }

      CUDA_FUNC static inline void calc_hess_ten(typename TrafoEvalHelp::HessianTensorType& hess_ten, const typename TrafoEvalHelp::DomainPointType& dom_point, const DataType (&coeffs)[dim][num_verts])
      {
        TrafoEvalHelp::calc_jac_mat(hess_ten, dom_point, coeffs);
      }

      CUDA_FUNC static inline DataType volume(const DataType (&coeffs)[dim][num_verts])
      {
        return TrafoEvalHelp::volume(coeffs);
      }

      CUDA_FUNC static inline DataType width_directed(const typename TrafoEvalHelp::ImagePointType& ray, const DataType (&coeffs)[dim][num_verts])
      {
        return TrafoEvalHelp::width_directed(ray, coeffs);
      }

      CUDA_FUNC static inline void trans_values(EvalData& data)
      {
        ParamEvalHelp::trans_values(data);
      }

      CUDA_FUNC static inline void trans_gradients(EvalData& data, const typename TrafoEvalHelp::JacobianInverseType& jac_inv)
      {
        ParamEvalHelp::trans_gradients(data, jac_inv);
      }

      CUDA_FUNC static inline void trans_hessian(EvalData& data, const typename TrafoEvalHelp::JacobianInverseType& jac_inv, const typename TrafoEvalHelp::HessianInverseType& hess_inv)
      {
        ParamEvalHelp::trans_gradients(data, jac_inv, hess_inv);
      }

    };

  }

}


#endif