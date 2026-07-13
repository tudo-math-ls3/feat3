// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/assembly/discrete_evaluator.hpp>
#include <kernel/base_header.hpp>
#include <kernel/global/remote_lambda.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/trafo/inverse_mapping.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <type_traits>
#include <utility>

namespace FEAT::Global
{
  namespace Intern
  {
    /// Get dimension of points in point range
    template<typename PointRange_>
    static constexpr int dim_from_point_range =
      std::remove_const_t<std::remove_reference_t<decltype(*std::declval<PointRange_>().begin())>>::n;

    /**
     * \brief Helper struct for determining types involved in a global evaluation of a finite-element function
     *
     * \tparam PointRange_ Type of point range
     * \tparam VectorType_ Type of vector that represents the finite-element function
     * \tparam Space_ Type of finite-element space
     * \tparam gradient_ Are we evaluating a gradient?
     *
     * This helper struct is mostly used to decide between ScalarDiscreteEvalData, VectorDiscreteEvalData, and MatrixDiscreteEvalData
     * based on the combination of VectorType_ and gradient_.
     */
    template<typename PointRange_, typename VectorType_, typename Space_, bool gradient_>
    struct DiscreteEvaluatorTypeHelper;

    // Specialization for non-gradient, non-blocked vector
    template<typename PointRange_, typename DataType_, typename IndexType_, typename Space_>
    struct DiscreteEvaluatorTypeHelper<PointRange_, LAFEM::DenseVector<DataType_, IndexType_>, Space_, false>
    {
      using DataType = DataType_;
      using IndexType = IndexType_;
      using PointType = std::remove_const_t<std::remove_reference_t<decltype(*std::declval<PointRange_>().begin())>>;
      static constexpr int point_dim = PointType::n;
      using SpaceType = Space_;
      using MeshType = typename SpaceType::MeshType;
      using EvalType = Assembly::ScalarDiscreteEvalData<DataType_, point_dim>;
      using EvalValueType = typename EvalType::ValueType;

      using BoundingBoxType = Tiny::Matrix<DataType, 2, point_dim>;
    };

    // Specialization for non-gradient, blocked vector
    template<typename PointRange_, typename DataType_, typename IndexType_, int dim_, typename Space_>
    struct DiscreteEvaluatorTypeHelper<
      PointRange_,
      LAFEM::DenseVectorBlocked<DataType_, IndexType_, dim_>,
      Space_,
      false>
    {
      using DataType = DataType_;
      using IndexType = IndexType_;
      using PointType = std::remove_const_t<std::remove_reference_t<decltype(*std::declval<PointRange_>().begin())>>;
      static constexpr int point_dim = PointType::n;
      using SpaceType = Space_;
      using MeshType = typename SpaceType::MeshType;
      using EvalType = Assembly::VectorDiscreteEvalData<DataType_, point_dim>;
      using EvalValueType = typename EvalType::ValueType;
      using BoundingBoxType = Tiny::Matrix<DataType, 2, point_dim>;

      static_assert(point_dim == dim_);
    };

    // Specialization for gradient, non-blocked vector
    template<typename PointRange_, typename DataType_, typename IndexType_, typename Space_>
    struct DiscreteEvaluatorTypeHelper<PointRange_, LAFEM::DenseVector<DataType_, IndexType_>, Space_, true>
    {
      using DataType = DataType_;
      using IndexType = IndexType_;
      using PointType = std::remove_const_t<std::remove_reference_t<decltype(*std::declval<PointRange_>().begin())>>;
      static constexpr int point_dim = PointType::n;
      using SpaceType = Space_;
      using MeshType = typename SpaceType::MeshType;
      using EvalType = Assembly::VectorDiscreteEvalData<DataType_, point_dim>;
      using EvalValueType = typename EvalType::ValueType;
      using BoundingBoxType = Tiny::Matrix<DataType, 2, point_dim>;
    };

    // Specialization for gradient, blocked vector
    template<typename PointRange_, typename DataType_, typename IndexType_, int dim_, typename Space_>
    struct DiscreteEvaluatorTypeHelper<
      PointRange_,
      LAFEM::DenseVectorBlocked<DataType_, IndexType_, dim_>,
      Space_,
      true>
    {
      using DataType = DataType_;
      using IndexType = IndexType_;
      using PointType = std::remove_const_t<std::remove_reference_t<decltype(*std::declval<PointRange_>().begin())>>;
      static constexpr int point_dim = PointType::n;
      using SpaceType = Space_;
      using MeshType = typename SpaceType::MeshType;
      using EvalType = Assembly::MatrixDiscreteEvalData<DataType_, point_dim, point_dim>;
      using EvalValueType = typename EvalType::ValueType;
      using BoundingBoxType = Tiny::Matrix<DataType, 2, point_dim>;

      static_assert(point_dim == dim_);
    };
  } // namespace Intern

  /**
   * \brief Global discrete evaluator
   *
   * This class implements various functions which can be used to evaluate
   * discrete (finite-element) functions defined on a mesh/trafo in arbitrary
   * points on a distributed domain.
   *
   * \sa Assembly::DiscreteEvaluator
   */
  class DiscreteEvaluator
  {
  public:
    /**
      * \brief Evaluates a scalar finite-element function in all given points.
      *
      * \tparam PointRange_ Type of point range
      * \tparam DT_ Data type of coefficient vector and evaluation results
      * \tparam IT_ Index type of coefficient vector
      * \tparam Space_ Space type of finite-element function
      *
      * \param[in] comm A \transient reference to the communicator that the domain is distributed over.
      *
      * \param[in] points A \transient reference to a range of points in real
      * world coordinates in which the finite element function is to be
      * evaluated.
      *
      * \param[in] vector The \transient coefficient vector of the
      * finite-element function.
      *
      * \param[in] space The \transient finite-element space corresponding to
      * the vector.
      *
      * \returns A vector of ScalarDiscreteEvalData objects that contain the
      * evaluation results. One ScalarDiscreteEvalData object is returned per
      * point in the point range. Each ScalarDiscreteEvalData object contains
      * all values determined by all ranks that managed to evaluate the function
      * in that point.
      */
    template<typename PointRange_, typename DT_, typename IT_, typename Space_>
    static std::vector<Assembly::ScalarDiscreteEvalData<DT_, Intern::dim_from_point_range<PointRange_>>>
    eval_fe_function(const Dist::Comm& comm, const PointRange_& points, const LAFEM::DenseVector<DT_, IT_>& vector, const Space_& space)
    {
      return _remote_eval<false>(comm, points, vector, space);
    }

    /**
      * \brief Evaluates a vector-valued finite-element function in all given points.
      *
      * \tparam PointRange_ Type of point range
      * \tparam DT_ Data type of coefficient vector and evaluation results
      * \tparam IT_ Index type of coefficient vector
      * \tparam Space_ Space type of finite-element function
      *
      * \param[in] comm A \transient reference to the communicator that the domain is distributed over.
      *
      * \param[in] points A \transient reference to a range of points in real
      * world coordinates in which the finite element function is to be
      * evaluated.
      *
      * \param[in] vector The \transient coefficient vector of the
      * finite-element function.
      *
      * \param[in] space The \transient finite-element space corresponding to
      * the vector.
      *
      * \returns A vector of VectorDiscreteEvalData objects that contain the
      * evaluation results. One VectorDiscreteEvalData object is returned per
      * point in the point range. Each VectorDiscreteEvalData object contains
      * all values determined by all ranks that managed to evaluate the function
      * in that point.
      */
    template<typename PointRange_, typename DT_, typename IT_, typename Space_>
    static std::vector<Assembly::VectorDiscreteEvalData<DT_, Intern::dim_from_point_range<PointRange_>>> eval_fe_function(
      const Dist::Comm& comm,
      const PointRange_& points,
      const LAFEM::DenseVectorBlocked<DT_, IT_, Intern::dim_from_point_range<PointRange_>>& vector,
      const Space_& space)
    {
      return _remote_eval<false>(comm, points, vector, space);
    }

    /**
      * \brief Evaluates the gradient of a scalar finite-element function in all given points.
      *
      * \tparam PointRange_ Type of point range
      * \tparam DT_ Data type of coefficient vector and evaluation results
      * \tparam IT_ Index type of coefficient vector
      * \tparam Space_ Space type of finite-element function
      *
      * \param[in] comm A \transient reference to the communicator that the domain is distributed over.
      *
      * \param[in] points A \transient reference to a range of points in real
      * world coordinates in which the finite element function is to be
      * evaluated.
      *
      * \param[in] vector The \transient coefficient vector of the
      * finite-element function.
      *
      * \param[in] space The \transient finite-element space corresponding to
      * the vector.
      *
      * \returns A vector of VectorDiscreteEvalData objects that contain the
      * evaluation results. One VectorDiscreteEvalData object is returned per
      * point in the point range. Each VectorDiscreteEvalData object contains
      * all values determined by all ranks that managed to evaluate the function
      * in that point.
      */
    template<typename PointRange_, typename DT_, typename IT_, typename Space_>
    static std::vector<Assembly::VectorDiscreteEvalData<DT_, Intern::dim_from_point_range<PointRange_>>>
    eval_fe_gradient(const Dist::Comm& comm, const PointRange_& points, const LAFEM::DenseVector<DT_, IT_>& vector, const Space_& space)
    {
      return _remote_eval<true>(comm, points, vector, space);
    }

    /**
      * \brief Evaluates the gradient of a vector-valued finite-element function in all given points.
      *
      * \tparam PointRange_ Type of point range
      * \tparam DT_ Data type of coefficient vector and evaluation results
      * \tparam IT_ Index type of coefficient vector
      * \tparam Space_ Space type of finite-element function
      *
      * \param[in] comm A \transient reference to the communicator that the domain is distributed over.
      *
      * \param[in] points A \transient reference to a range of points in real
      * world coordinates in which the finite element function is to be
      * evaluated.
      *
      * \param[in] vector The \transient coefficient vector of the
      * finite-element function.
      *
      * \param[in] space The \transient finite-element space corresponding to
      * the vector.
      *
      * \returns A vector of MatrixDiscreteEvalData objects that contain the
      * evaluation results. One MatrixDiscreteEvalData object is returned per
      * point in the point range. Each MatrixDiscreteEvalData object contains
      * all values determined by all ranks that managed to evaluate the function
      * in that point.
      */
    template<typename PointRange_, typename DT_, typename IT_, typename Space_>
    static std::vector<Assembly::MatrixDiscreteEvalData<
      DT_,
      Intern::dim_from_point_range<PointRange_>,
      Intern::dim_from_point_range<PointRange_>>>
    eval_fe_gradient(
      const Dist::Comm& comm,
      const PointRange_& points,
      const LAFEM::DenseVectorBlocked<DT_, IT_, Intern::dim_from_point_range<PointRange_>>& vector,
      const Space_& space)
    {
      return _remote_eval<true>(comm, points, vector, space);
    }

  private:
    /**
      * \brief Common function for evaluating a scalar or vector-valued finite-element function on a distributed domain.
      *
      * \tparam gradient_ Are we evaluating the function or its gradient?
      * \tparam PointRange_ Type of point range
      * \tparam VectorType_ Type of vector
      * \tparam Space_ Space type of finite-element function
      *
      * \param[in] comm A \transient reference to the communicator that the domain is distributed over.
      *
      * \param[in] points A \transient reference to a range of points in real
      * world coordinates in which the finite element function is to be
      * evaluated.
      *
      * \param[in] vector The \transient coefficient vector of the
      * finite-element function.
      *
      * \param[in] space The \transient finite-element space corresponding to
      * the vector.
      *
      * \returns A vector of *DiscreteE*valData objects that contain the
      * evaluation results. One *DiscreteEvalData object is returned per
      * point in the point range. Each *DiscreteEvalData object contains
      * all values determined by all ranks that managed to evaluate the function
      * in that point.
      */
    template<bool gradient_, typename PointRange_, typename VectorType_, typename Space_>
    static std::vector<typename Intern::DiscreteEvaluatorTypeHelper<PointRange_, VectorType_, Space_, gradient_>::EvalType>
    _remote_eval(const Dist::Comm& comm, const PointRange_& points, const VectorType_& vector, const Space_& space)
    {
      // We need to know the types of the values produced by evaluating the finite-element function
      // and the appropriate *DiscreteEvalData wrapper. The DiscreteEvaluatorTypeHelper
      // provides the correct types based on the VectorType_ and gradient_ params.
      using Helper = Intern::DiscreteEvaluatorTypeHelper<PointRange_, VectorType_, Space_, gradient_>;

      using DataType = typename Helper::DataType;
      using BoundingBoxType = typename Helper::BoundingBoxType;
      using PointType = typename Helper::PointType;
      using EvalType = typename Helper::EvalType;
      using EvalValueType = typename Helper::EvalValueType;

      const auto total_ranks = std::size_t(comm.size());
      const auto rank = std::size_t(comm.rank());
      const auto num_points = std::size_t(std::distance(points.begin(), points.end()));

      // Exchange bounding boxes with other ranks
      std::vector<BoundingBoxType> bbs(static_cast<std::size_t>(total_ranks));
      bbs.at(rank) = get_local_bb<BoundingBoxType>(space.get_mesh());
      comm.allgather(bbs.data(), 1, bbs.data(), 1);

      // Build structure for exchanging points
      std::vector<int> ranks;
      ranks.reserve(num_points);
      std::vector<Index> indices;
      ranks.reserve(num_points);

      Index point_idx(0);
      for(auto it = points.begin(); it != points.end(); ++it)
      {
        const PointType& p = *it;
        for(std::size_t i(0); i < total_ranks; i++)
        {
          const BoundingBoxType& bb = bbs[i];

          if(is_in_bb(bb, p))
          {
            ranks.push_back(int(i));
            indices.push_back(point_idx);
          }
        }
        point_idx++;
      }

      // Cache inverse mapping
      // The inverse mapping creates a list of bounding boxes for each cell.
      // This way we create this list only once per rank, rather than once per point
      Trafo::InverseMapping<typename Space_::TrafoType, DataType> inv_mapping(space.get_trafo());

      // Make remote calls. Each ranks returns a vector of values for each point, in case a rank manages to unmap the point on multiple cells
      ScalarRemoteLambda<EvalValueType(const PointType&)> rl(comm);
      std::vector<EvalValueType> evals = rl.call(
        ranks,
        indices,
        points,
        [&vector, &space, &inv_mapping](const PointType& p)
        {
          auto unmapping = inv_mapping.unmap_point(p);
          if constexpr(gradient_)
          {
            return Assembly::DiscreteEvaluator::eval_fe_gradient(unmapping, vector, space).mean_value();
          }
          else
          {
            return Assembly::DiscreteEvaluator::eval_fe_function(unmapping, vector, space).mean_value();
          }
        });

      // Assemble final results from individual calls
      std::vector<EvalType> result(std::size_t(std::distance(points.begin(), points.end())));

      for(std::size_t i(0); i < evals.size(); i++)
      {
        // To which point did this remote lambda call belong?
        const Index idx = indices[i];

        // Add evaluation result to total results for that point
        result[idx].values.push_back(evals[i]);
      }

      return result;
    }

    /// Check if a point is in a bounding box
    template<typename BoundingBoxType_, typename PointType_>
    static bool is_in_bb(const BoundingBoxType_& bb, const PointType_& p)
    {
      using DataType = typename BoundingBoxType_::ValueType;
      constexpr int dim = BoundingBoxType_::m;

      constexpr DataType bb_eps(1E-2);
      for(int i(0); i < dim; i++)
      {
        if(p[i] > (bb[1][i] + bb_eps) || p[i] < (bb[0][i] - bb_eps))
        {
          return false;
        }
      }
      return true;
    }

    /// Get bounding box of local part of domain
    template<typename BoundingBoxType_, typename MeshType_>
    static BoundingBoxType_ get_local_bb(const MeshType_& mesh)
    {
      const auto& vertex_set = mesh.get_vertex_set();

      if(vertex_set.get_num_vertices() == 0)
      {
        return BoundingBoxType_{0};
      }

      // Determine local bounding box
      BoundingBoxType_ local_bb;

      local_bb[0] = vertex_set[0];
      local_bb[1] = vertex_set[0];
      for(Index vtx_idx(0); vtx_idx < vertex_set.get_num_vertices(); vtx_idx++)
      {
        const auto& v = vertex_set[vtx_idx];

        for(int i(0); i < MeshType_::world_dim; i++)
        {
          Math::mini(local_bb[0][i], v[i]);
          Math::maxi(local_bb[1][i], v[i]);
        }
      }

      return local_bb;
    }
  }; // class InverseMapping
} // namespace FEAT::Global
