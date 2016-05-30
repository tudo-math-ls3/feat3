#pragma once
#ifndef KERNEL_SPACE_DISCONTINUOUS_ELEMENT_HPP
#define KERNEL_SPACE_DISCONTINUOUS_ELEMENT_HPP 1

// includes, FEAT
#include <kernel/space/element_base.hpp>
#include <kernel/space/dof_assignment_base.hpp>
#include <kernel/space/dof_mapping_common.hpp>
#include <kernel/space/discontinuous/dof_traits.hpp>
#include <kernel/space/discontinuous/evaluator.hpp>
#include <kernel/space/discontinuous/node_functional.hpp>

namespace FEAT
{
  namespace Space
  {
    /**
     * \brief Discontinuous Element namespace
     */
    namespace Discontinuous
    {
      /**
       * \brief Discontinuous Finite-Element space class template
       *
       * \tparam Trafo_
       * The transformation that is to be used by this space.
       *
       * \author Peter Zajac
       */
      template<
        typename Trafo_,
        typename Variant_ = Variant::StdPolyP<0> >
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
        /// variant of the element
        typedef Variant_ VariantTag;

        /// node functionals available
        static constexpr bool have_node_func = true;

        /** \copydoc ElementBase::local_degree */
        static constexpr int local_degree = VariantTag::local_degree;

        /// number of local dofs
        static constexpr int num_local_dofs = DofTraits<DofTag<ShapeType, VariantTag>, ShapeType::dimension>::count;

        /** \copydoc ElementBase::Evaluator */
        template<
          typename TrafoEvaluator_,
          typename DataType_ = typename TrafoEvaluator_::DataType>
        class Evaluator
        {
        private:
          /// evaluation policy
          typedef typename TrafoEvaluator_::EvalPolicy EvalPolicy;

          /// space evaluation traits
          typedef StandardScalarEvalTraits<EvalPolicy, num_local_dofs, DataType_> Traits;

        public:
          /// space evaluator type
          typedef Discontinuous::Evaluator<Element, TrafoEvaluator_, Traits, VariantTag> Type;
        };

        /** \copydoc ElementBase::DofMappingType */
        typedef DofMappingSingleEntity<Element, 0, num_local_dofs> DofMappingType;

        /** \copydoc ElementBase::DofAssignment */
        template<
          int shape_dim_,
          typename DataType_ = Real>
        class DofAssignment
        {
        public:
          /// Dof-Assignment type
          //typedef DofAssignmentSingleEntity<Element, shape_dim_, DataType_, ShapeType::dimension> Type;
          typedef DofAssignmentUniform<Element, shape_dim_, DataType_, DofTraits, DofTag<ShapeType, VariantTag> > Type;
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
          typedef Discontinuous::NodeFunctional<Element, codim, VariantTag, DataType_> Type;
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
          // number of DOFs = number of cells in the mesh
          return this->get_mesh().get_num_entities(ShapeType::dimension) * Index(num_local_dofs);
        }

        /** \copydoc ElementBase::name() */
        static String name()
        {
          return String("Discontinuous<") + Variant_::name() + ">";
        }

      }; // class Element
    } // namespace Discontinuous
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_DISCONTINUOUS_ELEMENT_HPP
