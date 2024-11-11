// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/util/math.hpp>

namespace FEAT
{
  namespace Assembly
  {

    /**
     * \brief Interpolation operator between two finite element spaces
     *
     * \tparam ToSpace_
     * The FE space to interpolate to
     *
     * \tparam FromSpace
     * The FE space to interpolate from
     *
     * This is the (currently empty) generic inter-FE space interpolator. For this, the node functionals from
     * ToSpace_ need to be applied to a function from FromSpace_, which might not be trivial (e.g. the Lagrange1
     * node functionals map to the function values at the mesh's vertices, where a P1~ or Q1~ function is not
     * continuous). Therefore, these operators shall be implemented by specializations in ToSpace_ and FromSpace_.
     *
     * \author Jordi Paul
     */
    template<typename ToSpace_, typename FromSpace_>
    struct FEInterpolator
    {
    };

    /**
     * \brief Helper class that recurses through shape dimensions and evaluates Lagrange2 DoF based on Lagrange1 DoF
     *
     * \tparam Trafo_
     * The type of the transformation
     *
     * \tparam shape_dim_
     * Shape dim for mapping Lagrange2 DoF, e.g. 1 for edges
     *
     * This is needed for the template recursion and applies the Lagrange2 node functionals for DoF associated with
     * entities of shape_dim_ to a Lagrange1 function. (This means it does linear interpolation at the edge, face and
     * cell midpoints if required by the shape type used).
     *
     * \author Jordi Paul
     */
    template<typename Trafo_, int shape_dim_ = Trafo_::MeshType::shape_dim>
    struct Lagrange1To2DofAtEntity
    {
      /// We map to Lagrange2
      typedef Space::Lagrange2::Element<Trafo_> ToSpace;
      /// We map from Lagrange1
      typedef Space::Lagrange1::Element<Trafo_> FromSpace;

      /**
       * \brief Evaluates the node functionals from shape_dim_ on down to shape dimension 0
       *
       * \tparam Vector_
       * The DoF vector type
       *
       * \param[in, out] to_coeffs
       * The \transient Lagrange2 DoF vector to be filled
       *
       * \param[in] from_coeffs
       * The \transient Lagrange1 DoF vector representing the FE function we interpolate
       *
       * \param[in] to_space
       * The \transient Lagrange2 space, needed for DoF mappings
       *
       * \note This explicitly uses that the DoF for shape dim 0 are numbered like the mesh's vertices. Without this,
       * it is not possible to find out which DoF of shape dim 0 lie at an edge for the interpolation.
       *
       */
      template<typename Vector_>
      static void recurse(Vector_& to_coeffs, const Vector_& from_coeffs, const ToSpace& to_space)
      {
        typedef typename Vector_::DataType DataType;

        const auto& mesh = to_space.get_trafo().get_mesh();

        // We use this to identify the vertex DoF at our shape
        const auto& from_vertex_at_shape_idx = mesh.template get_index_set<shape_dim_, 0>();

        // This is for mapping to the correct Lagrange2 DoF
        typename ToSpace::template DofAssignment<shape_dim_, DataType>::Type to_dof_assignment(to_space);

        for(Index k(0); k < mesh.get_num_entities(shape_dim_); ++k)
        {
          to_dof_assignment.prepare(k);

          // There is only one Lagrange2 DoF at each entity, but this is easier to understand for future reference
          for(int j(0); j < to_dof_assignment.get_num_assigned_dofs(); ++j)
          {

            typename Vector_::ValueType tmp(0);

            for(int i(0); i < from_vertex_at_shape_idx.get_num_indices(); ++i)
            {
              tmp += from_coeffs(from_vertex_at_shape_idx(k,i));
            }

            // Linear interpolation
            tmp *= ( DataType(1)/DataType(from_vertex_at_shape_idx.get_num_indices()) );

            // Emplacement operator
            to_coeffs(to_dof_assignment.get_index(j), tmp);
          }
        }

        // recurse down
        Lagrange1To2DofAtEntity<Trafo_, shape_dim_-1>::recurse( to_coeffs, from_coeffs, to_space);
      }
    };

    /**
     * \brief End of helper class template recursion
     *
     * \tparam Trafo_
     * The type of the transformation
     *
     */
    template<typename Trafo_>
    struct Lagrange1To2DofAtEntity<Trafo_, 0>
    {
      typedef Space::Lagrange2::Element<Trafo_> ToSpace;
      typedef Space::Lagrange1::Element<Trafo_> FromSpace;

      /**
       * \brief Evaluates the node functionals for shape dimension 0
       *
       * \tparam Vector_
       * The DoF vector type
       *
       * \param[in, out] to_coeffs
       * The \transient Lagrange2 DoF vector to be filled
       *
       * \param[in] from_coeffs
       * The \transient Lagrange1 DoF vector representing the FE function we interpolate
       *
       * \param[in] to_space
       * The \transient Lagrange2 space, needed for DoF mappings
       *
       * \note This explicitly uses that the DoF for shape dim 0 are numbered like the mesh's vertices. Without this,
       * it is not possible to find out which DoF of shape dim 0 lie at an edge for the interpolation.
       *
       */
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
     * \brief Interpolator class from Lagrange1 to Lagrange2
     *
     * \tparam Trafo_
     * The transformation type.
     *
     * \author Jordi Paul
     */
    template<typename Trafo_>
    struct FEInterpolator<Space::Lagrange2::Element<Trafo_>, Space::Lagrange1::Element<Trafo_>>
    {
      /// Our trafo
      typedef Trafo_ TrafoType;
      /// We map to Lagrange2
      typedef Space::Lagrange2::Element<Trafo_> ToSpace;
      /// We map from Lagrange1
      typedef Space::Lagrange1::Element<Trafo_> FromSpace;

      /**
       * \brief Interpolates a scalar Lagrange1 to Lagrange2 FE function
       *
       * \tparam DT_
       * Floating point precision of the DoF vector
       *
       * \tparam IT_
       * Index type of the DoF vector
       *
       * \param[out] to_coeffs
       * The \transient coefficient vector of the Lagrange2 FE function
       *
       * \param[in] from_coeffs
       * The \transient coefficient vector of the Lagrange1 FE function
       *
       * \param[in] to_space
       * The \transient space we map to, needed for the DofMapping and DofAssignment
       *
       * \param[in] from_space
       * The \transient space we map from, unused here because we use the IndexSet of the
       * corresponding mesh instead of the DofMapping and DofAssignment
       */
      template<typename DT_, typename IT_>
      static void interpolate(LAFEM::DenseVector<DT_, IT_>& to_coeffs,
      const LAFEM::DenseVector<DT_, IT_>& from_coeffs,
      const ToSpace& to_space,
      const FromSpace& DOXY(from_space))
      {
        Lagrange1To2DofAtEntity<TrafoType>::recurse(to_coeffs, from_coeffs, to_space);
      }

      /**
       * \brief Interpolates a vector valued Lagrange1 to Lagrange2 FE function
       *
       * \tparam DT_
       * Floating point precision of the DoF vector
       *
       * \tparam IT_
       * Index type of the DoF vector
       *
       * \param[out] to_coeffs
       * The \transient coefficient vector of the Lagrange2 FE function
       *
       * \param[in] from_coeffs
       * The \transient coefficient vector of the Lagrange1 FE function
       *
       * \param[in] to_space
       * The \transient space we map to, needed for the DofMapping and DofAssignment
       *
       * \param[in] from_space
       * The \transient space we map from, unused here because we use the IndexSet of the
       * corresponding mesh instead of the DofMapping and DofAssignment
       */
      template<typename DT_, typename IT_, int blocksize_>
      static void interpolate(LAFEM::DenseVectorBlocked<DT_, IT_, blocksize_>& to_coeffs,
      const LAFEM::DenseVectorBlocked<DT_, IT_, blocksize_>& from_coeffs,
      const ToSpace& to_space,
      const FromSpace& DOXY(from_space))
      {
        Lagrange1To2DofAtEntity<TrafoType>::recurse(to_coeffs, from_coeffs, to_space);
      }
    }; //FEInterpolator Lagrange1 to Lagrange2

    /**
     * \brief Interpolator class from Lagrange2 to Lagrange1
     *
     * \tparam Trafo_
     * The transformation type.
     *
     * \author Maximilian Esser
     */
    template<typename Trafo_>
    struct FEInterpolator<Space::Lagrange1::Element<Trafo_>, Space::Lagrange2::Element<Trafo_>>
    {
      /// Our trafo
      typedef Trafo_ TrafoType;
      /// We map to Lagrange1
      typedef Space::Lagrange1::Element<Trafo_> ToSpace;
      /// We map from Lagrange2
      typedef Space::Lagrange2::Element<Trafo_> FromSpace;

      /**
       * \brief Interpolates a scalar Lagrange1 to Lagrange2 FE function
       *
       * \tparam DT_
       * Floating point precision of the DoF vector
       *
       * \tparam IT_
       * Index type of the DoF vector
       *
       * \param[out] to_coeffs
       * The coefficient vector of the Lagrange2 FE function
       *
       * \param[in] from_coeffs
       * The coefficient vector of the Lagrange1 FE function
       *
       * \param[in] to_space
       * The space we map to, needed for the DofMapping and DofAssignment
       *
       * \param[in] from_space
       * The space we map from, unused here because we use the IndexSet of the corresponding mesh instead of the
       * DofMapping and DofAssignment
       */
      template<typename DT_, typename IT_>
      static void interpolate(LAFEM::DenseVector<DT_, IT_>& to_coeffs,
      const LAFEM::DenseVector<DT_, IT_>& from_coeffs,
      const ToSpace& DOXY(to_space),
      const FromSpace& DOXY(from_space))
      {
        //sanity check for the vector sizes
        XASSERTM(to_coeffs.size() < from_coeffs.size(), "Coefficient vectors do not match!\n To_coeffs size is ");
        //due to the two level ordering of Lagrange1 and Lagrange2 elements it is always sufficient to just truncate our vector
        for(Index i(0); i < to_coeffs.size(); ++i)
        {
          to_coeffs(i, from_coeffs(i));
        }
        //And we are done
      }

      /**
       * \brief Interpolates a vector Lagrange1 to Lagrange2 FE function
       *
       * \tparam DT_
       * Floating point precision of the DoF vector
       *
       * \tparam IT_
       * Index type of the DoF vector
       *
       * \param[out] to_coeffs
       * The coefficient vector of the Lagrange2 FE function
       *
       * \param[in] from_coeffs
       * The coefficient vector of the Lagrange1 FE function
       *
       * \param[in] to_space
       * The space we map to, needed for the DofMapping and DofAssignment
       *
       * \param[in] from_space
       * The space we map from, unused here because we use the IndexSet of the corresponding mesh instead of the
       * DofMapping and DofAssignment
       */
      template<typename DT_, typename IT_, int blocksize_>
      static void interpolate(LAFEM::DenseVectorBlocked<DT_, IT_, blocksize_>& to_coeffs,
      const LAFEM::DenseVectorBlocked<DT_, IT_, blocksize_>& from_coeffs,
      const ToSpace& DOXY(to_space),
      const FromSpace& DOXY(from_space))
      {
        //sanity check for the vector sizes
        XASSERTM(to_coeffs.size() < from_coeffs.size(), "Coefficient vectors do not match!\n To_coeffs size is ");
        //due to the two level ordering of Lagrange1 and Lagrange2 elements it is always sufficient to only truncate our vector
        for(Index i(0); i < to_coeffs.size(); ++i)
        {
          to_coeffs(i, from_coeffs(i));
        }
        //And we are done
      }

    }; //FEInterpolator Lagrange2 to Lagrange1

    /**
     * \brief Interpolator class for trivial interpolator between same class
     *
     * \tparam Space_
     * The Space_ of both input and output.
     *
     * \author Maximilian Esser
     */
    template<typename Space_>
    struct FEInterpolator<Space_, Space_>
    {
      typedef Space_ ToSpace;

      typedef Space_ FromSpace;

      /**
       * \brief Interpolates a scalar/vector FE function on the same space by just copying the vector
       *
       * \warning This implementation is only due to consistency and should be avoided
       *
       * \tparam DT_
       * Floating point precision of the DoF vector
       *
       * \tparam IT_
       * Index type of the DoF vector
       *
       * \param[out] to_coeffs
       * The coefficient vector of the Lagrange2 FE function
       *
       * \param[in] from_coeffs
       * The coefficient vector of the Lagrange1 FE function
       *
       * \param[in] to_space
       * The space we map to, needed for the DofMapping and DofAssignment
       *
       * \param[in] from_space
       * The space we map from, unused here because we use the IndexSet of the corresponding mesh instead of the
       * DofMapping and DofAssignment
       */
      template<typename Vector_>
      static void interpolate(Vector_& to_coeffs,
      const Vector_& from_coeffs,
      const ToSpace& DOXY(to_space),
      const FromSpace& DOXY(from_space))
      {
        to_coeffs.copy(from_coeffs, true);
      }

    }; //FEInterpolator same space

  } // namespace Assembly
} // namespace FEAT
