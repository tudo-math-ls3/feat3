#pragma once
#ifndef KERNEL_ASSEMBLY_DISCRETE_EVALUATOR_HPP
#define KERNEL_ASSEMBLY_DISCRETE_EVALUATOR_HPP

#include <kernel/assembly/asm_traits.hpp>
#include <kernel/trafo/inverse_mapping.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/util/dist.hpp>

#include <vector>

namespace FEAT
{
  namespace Assembly
  {
    /**
     * \brief Discrete evaluation data for scalar functions
     *
     * \tparam DT_
     * The datatype that is to be used.
     *
     * \tparam dim_
     * The dimension of the underlying mesh.
     *
     * \author Peter Zajac
     */
    template<typename DT_, int dim_>
    class ScalarDiscreteEvalData
    {
    public:
      /// the image dimension
      static constexpr int dim = dim_;
      /// the data type
      typedef DT_ DataType;
      /// the value type
      typedef DataType ValueType;

      /// the vector of values
      std::vector<ValueType> values;

      /**
       * \brief Checks whether the evaluation data is empty.
       */
      bool empty() const
      {
        return values.empty();
      }

      /**
       * \brief Computes the mean function value.
       *
       * \attention
       * This function does not work on distributed domains!\n
       * Use the #mean_value_dist() function in this case instead!
       *
       * \returns
       * The mean function value.
       */
      ValueType mean_value() const
      {
        if(values.empty())
          return DT_(0);

        // compute and return mean value
        DT_ r = DT_(0);
        for(const auto& v : values)
          r += v;
        return r / DT_(values.size());
      }

      /**
       * \brief Computes the mean function value on a distributed domain.
       *
       * \param[in] comm
       * The communicator that the domain was distributed on.
       *
       * \returns
       * The mean function value.
       */
      ValueType mean_value_dist(const Dist::Comm& comm) const
      {
        // In a distributed domain, we need to compute the mean value
        // over all patches, as the evaluation point may intersect
        // with several patches. Therefore, we also need to sum up
        // the number of values in addition to the values themselves,
        // so that we can compute the mean value over all patches.
        DT_ val[2] =
        {
          DT_(values.size()),
          DT_(0)
        };

        // add all values from this patch
        for(const auto& v : values)
          val[1] += v;

        // sum up over all patches
        comm.allreduce(val, val, std::size_t(2), Dist::op_sum);

        // compute global mean value
        return (val[0] > DT_(0)) ? (val[1] / val[0]) : DT_(0);
      }
    }; // class ScalarDiscreteEvalData

    /**
     * \brief Discrete evaluation data for vecotr-valued functions
     *
     * \tparam DT_
     * The datatype that is to be used.
     *
     * \tparam dim_
     * The dimension of the underlying mesh.
     *
     * \author Peter Zajac
     */
    template<typename DT_, int dim_>
    class VectorDiscreteEvalData
    {
    public:
      /// the image dimension
      static constexpr int dim = dim_;
      /// the data type
      typedef DT_ DataType;
      /// the value type
      typedef Tiny::Vector<DT_, dim_> ValueType;

      /// the vector of values
      std::vector<ValueType> values;

      /**
       * \brief Checks whether the evaluation data is empty.
       */
      bool empty() const
      {
        return values.empty();
      }

      /**
       * \brief Computes the mean function value.
       *
       * \attention
       * This function does not work on distributed domains!\n
       * Use the #mean_value_dist() function in this case instead!
       *
       * \returns
       * The mean function value.
       */
      ValueType mean_value() const
      {
        ValueType r;
        r.format();
        if(values.empty())
          return r;

        for(const auto& v : values)
          r += v;
        r  *= (DT_(1) / DT_(values.size()));

        return r;
      }

      /**
       * \brief Computes the mean function value on a distributed domain.
       *
       * \param[in] comm
       * The communicator that the domain was distributed on.
       *
       * \returns
       * The mean function value.
       */
      ValueType mean_value_dist(const Dist::Comm& comm) const
      {
        ValueType r;
        r.format();

        DT_ val[1+dim_] =
        {
          DT_(values.size()),
        };

        for(int i(0); i < dim_; ++i)
          val[i+1] = DT_(0);

        for(const auto& v : values)
        {
          for(int i(0); i < dim_; ++i)
            val[i+1] += v[i];
        }

        comm.allreduce(val, val, std::size_t(dim_+1), Dist::op_sum);

        for(int i(0); i < dim_; ++i)
          r[i] = val[i+1];

        if(val[0] > DT_(0))
          r *= (DT_(1) / val[0]);

        return r;
      }
    }; // class VectorDiscreteEvalData

    /**
     * \brief Discrete function evaluator
     *
     * This class implements various functions which can be used to evaluate
     * discrete (finite-element) functions defined on a mesh/trafo in arbitrary
     * points.
     *
     * <b>Remarks:</b>
     * A basic problem with discrete (finite-element) function evaluation is that
     * one single point may intersect several cells of the mesh and, depending on
     * the finite-element space, the discrete function may not be uniquely defined
     * in this case. On the other hand, the given input point may be outside of
     * the domain represented by the underlying mesh (or the patch in the parallel
     * domain decomposition case), in which case the discrete function is not
     * defined at all.
     *
     * In consequence, the functions provided by this class do not return a single
     * value, but a ScalarDiscreteEvalData or VectorDiscreteEvalData object instead,
     * depending on whether the discrete function is a scalar function or a vector
     * field. These objects contain a vector of all cells, which have been found
     * to intersect with the given input point, as well as the corresponding
     * values of the discrete function as evaluated on the cells. It is up to
     * you to decide what has to be done if a given input point was found to
     * intersect with more than one element -- or no element at all.
     * In most cases, you will simply want to compute the average of all these
     * values, which can be computed by the mean_value() functions of the
     * ScalarDiscreteEvalData/VectorDiscreteEvalData classes.
     *
     * \note
     * This class is strongly related to the Trafo::InverseMapping class template,
     * which is used for the unmapping of the input point, so you might want to
     * get familiar with that one, too.
     *
     * \author Peter Zajac
     */
    class DiscreteEvaluator
    {
    public:
      /**
       * \brief Evaluates a scalar finite-element function in a given point.
       *
       * \note
       * If you intend to evaluate several finite-element functions in the same point,
       * then you might consider unmapping the input point by yourself and using the
       * other overload of this function instead for performance reasons.
       *
       * \param[in] point
       * The point in real world coordinates in which the finite element function
       * is to be evaluated.
       *
       * \param[in] vector
       * The coefficient vector of the finite-element function.
       *
       * \param[in] space
       * The finite-element space corresponding to the vector.
       *
       * \returns
       * A ScalarDiscreteEvalData object that contains the evaluations results.
       */
      template<typename DT_, typename DTP_, typename IT_, int dim_, int s_, typename Space_>
      static ScalarDiscreteEvalData<DT_, dim_> eval_fe_function(
        const Tiny::Vector<DTP_, dim_, s_>& point,
        const LAFEM::DenseVector<Mem::Main, DT_, IT_>& vector,
        const Space_& space)
      {
        // create inverse trafo mapping
        Trafo::InverseMapping<typename Space_::TrafoType, DT_> inv_mapping(space.get_trafo());

        // unmap point
        auto inv_map_data = inv_mapping.unmap_point(point);

        // evaluate FE function
        return eval_fe_function(inv_map_data, vector, space);
      }

      /**
       * \brief Evaluates a scalar finite-element function in a given unmapped point.
       *
       * \param[in] inv_map_data
       * A Trafo::InverseMappingData object that represents the unmapped evaluation point,
       * as computed by the Trafo::InverseMapping class.
       *
       * \param[in] vector
       * The coefficient vector of the finite-element function.
       *
       * \param[in] space
       * The finite-element space.
       *
       * \returns
       * A ScalarDiscreteEvalData object that contains the evaluations results.
       */
      template<typename DT_, typename DTP_, typename IT_, int dim_, typename Space_>
      static ScalarDiscreteEvalData<DT_, dim_> eval_fe_function(
        const Trafo::InverseMappingData<DTP_, dim_, dim_>& inv_map_data,
        const LAFEM::DenseVector<Mem::Main, DT_, IT_>& vector,
        const Space_& space)
      {
        // vector type
        typedef LAFEM::DenseVector<Mem::Main, DT_, IT_> VectorType;
        // space type
        typedef Space_ SpaceType;
        // assembly traits
        typedef AsmTraits1<DT_, SpaceType, TrafoTags::none, SpaceTags::value> AsmTraits;
        // data type
        typedef typename AsmTraits::DataType DataType;

        // fetch the trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create local vector data
        typename AsmTraits::LocalVectorType local_vector;

        // create matrix scatter-axpy
        typename VectorType::GatherAxpy gather_axpy(vector);

        // create evaluation data
        ScalarDiscreteEvalData<DT_, dim_> eval_data;

        // loop over all cells of the mesh
        for(std::size_t i(0); i < inv_map_data.size(); ++i)
        {
          const Index& cell = inv_map_data.cells.at(i);
          const auto& dom_point = inv_map_data.dom_points.at(i);

          // format local vector
          local_vector.format();

          // initialise dof-mapping
          dof_mapping.prepare(cell);

          // gather local vector data
          gather_axpy(local_vector, dof_mapping);

          // finish dof-mapping
          dof_mapping.finish();

          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // compute trafo data
          trafo_eval(trafo_data, dom_point);

          // compute basis function data
          space_eval(space_data, trafo_data);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // compute function value
          DataType value = DataType(0);
          for(int j(0); j < num_loc_dofs; ++j)
            value += local_vector[j] * space_data.phi[j].value;

          // finally, add this result to the eval data
          eval_data.values.push_back(value);

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();

          // continue with next cell
        }

        // return evaluation data
        return eval_data;
      }

      /**
       * \brief Evaluates a vector-valued finite-element function in a given point.
       *
       * \note
       * If you intend to evaluate several finite-element functions in the same point,
       * then you might consider unmapping the input point by yourself and using the
       * other overload of this function instead for performance reasons.
       *
       * \param[in] point
       * The point in real world coordinates in which the finite element function
       * is to be evaluated.
       *
       * \param[in] vector
       * The coefficient vector of the finite-element function.
       *
       * \param[in] space
       * The finite-element space.
       *
       * \returns
       * A VectorDiscreteEvalData object that contains the evaluations results.
       */
      template<typename DT_, typename DTP_, typename IT_, int dim_, int s_, typename Space_>
      static VectorDiscreteEvalData<DT_, dim_> eval_fe_function(
        const Tiny::Vector<DTP_, dim_, s_>& point,
        const LAFEM::DenseVectorBlocked<Mem::Main, DT_, IT_, dim_>& vector,
        const Space_& space)
      {
        // create inverse trafo mapping
        Trafo::InverseMapping<typename Space_::TrafoType, DT_> inv_mapping(space.get_trafo());

        // unmap point
        auto inv_map_data = inv_mapping.unmap_point(point);

        // evaluate FE function
        return eval_fe_function(inv_map_data, vector, space);
      }

      /**
       * \brief Evaluates a vector-valued finite-element function in a given unmapped point.
       *
       * \param[in] inv_map_data
       * A Trafo::InverseMappingData object that represents the unmapped evaluation point,
       * as computed by the Trafo::InverseMapping class.
       *
       * \param[in] vector
       * The coefficient vector of the finite-element function.
       *
       * \param[in] space
       * The finite-element space.
       *
       * \returns
       * A VectorDiscreteEvalData object that contains the evaluations results.
       */
      template<typename DT_, typename DTP_, typename IT_, int dim_, typename Space_>
      static VectorDiscreteEvalData<DT_, dim_> eval_fe_function(
        const Trafo::InverseMappingData<DTP_, dim_, dim_>& inv_map_data,
        const LAFEM::DenseVectorBlocked<Mem::Main, DT_, IT_, dim_>& vector,
        const Space_& space)
      {
        // vector type
        typedef LAFEM::DenseVectorBlocked<Mem::Main, DT_, IT_, dim_> VectorType;
        // space type
        typedef Space_ SpaceType;
        // assembly traits
        typedef AsmTraits1<DT_, SpaceType, TrafoTags::none, SpaceTags::value> AsmTraits;
        // data type
        typedef typename AsmTraits::DataType DataType;

        // fetch the trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // get maximum number of local dofs
        static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> ValueType;
        typedef Tiny::Vector<ValueType, max_local_dofs> LocalVectorType;
        LocalVectorType local_vector;

        // create matrix scatter-axpy
        typename VectorType::GatherAxpy gather_axpy(vector);

        // create evaluation data
        VectorDiscreteEvalData<DT_, dim_> eval_data;

        ValueType value;

        // loop over all cells of the mesh
        for(std::size_t i(0); i < inv_map_data.size(); ++i)
        {
          const Index& cell = inv_map_data.cells.at(i);
          const auto& dom_point = inv_map_data.dom_points.at(i);

          // format local vector
          local_vector.format();

          // initialise dof-mapping
          dof_mapping.prepare(cell);

          // gather local vector data
          gather_axpy(local_vector, dof_mapping);

          // finish dof-mapping
          dof_mapping.finish();

          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // compute trafo data
          trafo_eval(trafo_data, dom_point);

          // compute basis function data
          space_eval(space_data, trafo_data);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // compute function value
          value.format();
          for(int j(0); j < num_loc_dofs; ++j)
            value.axpy(space_data.phi[j].value, local_vector[j]);

          // finally, add this result to the eval data
          eval_data.values.push_back(value);

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();

          // continue with next cell
        }

        // return evaluation data
        return eval_data;
      }
    }; // class DiscreteEvaluator<...>
  } // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_DISCRETE_EVALUATOR_HPP
