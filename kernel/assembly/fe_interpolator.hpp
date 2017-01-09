#pragma once
#ifndef KERNEL_ASSEMBLY_FE_INTERPOLATOR_HPP
#define KERNEL_ASSEMBLY_FE_INTERPOLATOR_HPP 1

// includes, FEAT
#include <kernel/analytic/function.hpp>
#include <kernel/assembly/asm_traits.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/util/math.hpp>

namespace FEAT
{
  namespace Assembly
  {

    template<typename Trafo_, int shape_dim_ = Trafo_::MeshType::shape_dim>
    struct Lagrange2DofAtEntity
    {
      typedef Space::Lagrange2::Element<Trafo_> ToSpace;
      typedef Space::Lagrange1::Element<Trafo_> FromSpace;

      template<typename Vector_>
      static void recurse(Vector_& to_coeffs, const Vector_& from_coeffs, const ToSpace& to_space)
      {
        typedef typename Vector_::DataType DataType;

        const auto& mesh = to_space.get_trafo().get_mesh();

        const auto& from_vertex_at_shape_idx = mesh.template get_index_set<shape_dim_, 0>();

        typename ToSpace::template DofAssignment<shape_dim_, DataType>::Type to_dof_assignment(to_space);

        for(Index k(0); k < mesh.get_num_entities(shape_dim_); ++k)
        {
          to_dof_assignment.prepare(k);

          for(int j(0); j < to_dof_assignment.get_num_assigned_dofs(); ++j)
          {

            typename Vector_::ValueType tmp(0);

            for(Index i(0); i < Index(from_vertex_at_shape_idx.get_num_indices()); ++i)
            {
              tmp += from_coeffs(from_vertex_at_shape_idx(k,i));
            }

            tmp *= ( DataType(1)/DataType(from_vertex_at_shape_idx.get_num_indices()) );

            to_coeffs(to_dof_assignment.get_index(j), tmp);
          }
        }

        Lagrange2DofAtEntity<Trafo_, shape_dim_-1>::recurse( to_coeffs, from_coeffs, to_space);
      }
    };

    template<typename Trafo_>
    struct Lagrange2DofAtEntity<Trafo_, 0>
    {

      typedef Space::Lagrange2::Element<Trafo_> ToSpace;
      typedef Space::Lagrange1::Element<Trafo_> FromSpace;

      template<typename Vector_>
      static void recurse(Vector_& to_coeffs, const Vector_& from_coeffs, const ToSpace& to_space)
      {

        const auto& mesh = to_space.get_trafo().get_mesh();

        for(Index i(0); i < mesh.get_num_entities(0); ++i)
        {
          to_coeffs(i, from_coeffs(i));
        }

      }
    };

    /**
    */
    template<typename ToSpace_, typename FromSpace_>
    struct FEInterpolator
    {
    };
    /**
      template
      <
        typename DT_,
        typename IT_,
      >
      static void interpolate(
        LAFEM::DenseVector<Mem::Main, DT_, IT_>& to_coeffs,
        LAFEM::DenseVector<Mem::Main, DT_, IT_>& from_coeffs,
        const ToSpace_& to_space,
        const FromSpace_& from_space
        );
        */


    template<typename Trafo_>
    struct FEInterpolator<Space::Lagrange2::Element<Trafo_>, Space::Lagrange1::Element<Trafo_>>
    {
      typedef Trafo_ TrafoType;
      typedef Space::Lagrange2::Element<Trafo_> ToSpace;
      typedef Space::Lagrange1::Element<Trafo_> FromSpace;

      template<typename DT_, typename IT_>
      static void interpolate(LAFEM::DenseVector<Mem::Main, DT_, IT_>& to_coeffs,
      LAFEM::DenseVector<Mem::Main, DT_, IT_>& from_coeffs,
      const ToSpace& to_space,
      const FromSpace& DOXY(from_space))
      {
        Lagrange2DofAtEntity<TrafoType>::recurse(to_coeffs, from_coeffs, to_space);
      }

      template<typename DT_, typename IT_, int blocksize_>
      static void interpolate(LAFEM::DenseVectorBlocked<Mem::Main, DT_, IT_, blocksize_>& to_coeffs,
      LAFEM::DenseVectorBlocked<Mem::Main, DT_, IT_, blocksize_>& from_coeffs,
      const ToSpace& to_space,
      const FromSpace& DOXY(from_space))
      {
        Lagrange2DofAtEntity<TrafoType>::recurse(to_coeffs, from_coeffs, to_space);
      }
    };

    //  template<
    //  typename Mem_,
    //  typename DT_,
    //  typename IT_,
    //  typename Function_,
    //  typename Space_,
    //  typename CubatureFactory_>
    //    static void interpolate(
    //      LAFEM::DenseVector<Mem_, DT_, IT_>& vector,
    //      const Function_& function,
    //      const Space_& space,
    //      const CubatureFactory_& cubature_factory,
    //      WeightType weight_type = wt_volume)
    //      {
    //        LAFEM::DenseVector<Mem::Main, DT_, IT_> vec;
    //        interpolate(vec, function, space, cubature_factory, weight_type);
    //        vector.convert(vec);
    //      }
    //}; // class FEInterpolator
} // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_FE_INTERPOLATOR_HPP
