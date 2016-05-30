#pragma once
#ifndef KERNEL_SPACE_CRO_RAV_RAN_TUR_ELEMENT_HPP
#define KERNEL_SPACE_CRO_RAV_RAN_TUR_ELEMENT_HPP 1

// includes, FEAT
#include <kernel/space/element_base.hpp>
#include <kernel/space/dof_assignment_base.hpp>
#include <kernel/space/dof_mapping_common.hpp>
#include <kernel/space/cro_rav_ran_tur/dof_traits.hpp>
#include <kernel/space/cro_rav_ran_tur/evaluator.hpp>
#include <kernel/space/cro_rav_ran_tur/node_functional.hpp>

namespace FEAT
{
  namespace Space
  {
    /**
     * \brief Crouzeix-Raviart / Rannacher-Turek element namespace
     */
    namespace CroRavRanTur
    {
      /// \cond internal
      namespace Intern
      {
        template<typename Shape_>
        struct LocalDegree;

        template<int dim_>
        struct LocalDegree<Shape::Simplex<dim_>>
        {
          static constexpr int value = 1;
        };

        template<int dim_>
        struct LocalDegree<Shape::Hypercube<dim_>>
        {
          static constexpr int value = 2;
        };
      } // namespace Intern
      /// \endcond

      /**
       * \brief Crouzeix-Raviart / Rannacher-Turek element class template
       *
       * This class template implemements the second-order H1-non-conforming finite element spaces
       * known as the Crouzeix-Raviart (Simplex shapes) and Rannacher-Turek (Hypercube shapes).
       *
       * \note
       * In the case of Hypercube meshes, this class implements the non-parameteric integral-mean
       * based implementation, which is known to be the most stable of all variants.
       *
       * \see \cite CR73, \cite RT92
       *
       * \tparam Trafo_
       * The transformation that is to be used by this space.
       *
       * \author Peter Zajac
       */
      template<typename Trafo_>
      class Element :
        public ElementBase<Trafo_>
      {
      public:
        /// base-class typedef
        typedef ElementBase<Trafo_> BaseClass;
        /// transformation type
        typedef Trafo_ TrafoType;
        /// mesh type
        typedef typename TrafoType::MeshType MeshType;
        /// shape type
        typedef typename TrafoType::ShapeType ShapeType;

        /// node functionals available
        static constexpr bool have_node_func = true;

        /** \copydoc ElementBase::local_degree */
        static constexpr int local_degree = Intern::LocalDegree<ShapeType>::value;

        /** \copydoc ElementBase::Evaluator */
        template<
          typename TrafoEvaluator_,
          typename DataType_ = typename TrafoEvaluator_::DataType>
        class Evaluator
        {
        private:
          /// evaluation policy
          typedef typename TrafoEvaluator_::EvalPolicy EvalPolicy;

          /// number of local dofs := number of facets per cell
          static constexpr int num_loc_dofs = Shape::FaceTraits<ShapeType, ShapeType::dimension-1>::count;

          /// space evaluation traits
          typedef StandardScalarEvalTraits<EvalPolicy, num_loc_dofs, DataType_> Traits;

        public:
          /// space evaluator type
          typedef CroRavRanTur::Evaluator<Element, TrafoEvaluator_, Traits> Type;
        };

        /** \copydoc ElementBase::DofMappingType */
        typedef DofMappingSingleEntity<Element, 1> DofMappingType;

        /** \copydoc ElementBase::DofAssignment */
        template<
          int shape_dim_,
          typename DataType_ = Real>
        class DofAssignment
        {
        public:
          /// Dof-Assignment type
          typedef DofAssignmentUniform<Element, shape_dim_, DataType_, DofTraits, ShapeType> Type;
        };

        /** \copydoc ElementBase::NodeFunctional */
        template<
          int shape_dim_,
          typename DataType_ = Real>
        class NodeFunctional
        {
        private:
          /// co-dimension
          static constexpr int codim = ShapeType::dimension - shape_dim_;

        public:
          /// node functional type
          typedef CroRavRanTur::NodeFunctional<Element, ShapeType, codim, DataType_> Type;
        };

      public:
        /**
         * \brief Constructor
         *
         * \param[in] trafo
         * A reference to the transformation which is to be used by this space.
         */
        explicit Element(TrafoType& trafo)
          : BaseClass(trafo)
        {
        }

        /// virtual destructor
        virtual ~Element()
        {
        }

        /** \copydoc ElementBase::get_num_dofs() */
        Index get_num_dofs() const
        {
          // number of DOFs = number of facets in the mesh
          return this->get_mesh().get_num_entities(ShapeType::dimension - 1);
        }

        /** \copydoc ElementBase::name() */
        static String name()
        {
          return "CroRavRanTur";
        }
      }; // class Element
    } // namespace CroRavRanTur
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_CRO_RAV_RAN_TUR_ELEMENT_HPP
