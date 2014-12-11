#pragma once
#ifndef KERNEL_SPACE_DISCONTINUOUS_ELEMENT_HPP
#define KERNEL_SPACE_DISCONTINUOUS_ELEMENT_HPP 1

// includes, FEAST
#include <kernel/space/element_base.hpp>
#include <kernel/space/dof_assignment_base.hpp>
#include <kernel/space/dof_mapping_common.hpp>
#include <kernel/space/discontinuous/dof_traits.hpp>
#include <kernel/space/discontinuous/evaluator.hpp>
#include <kernel/space/discontinuous/node_functional.hpp>

namespace FEAST
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

        /** \copydoc ElementBase::ElementCapabilities */
        enum ElementCapabilities
        {
          /// node functionals available
          have_node_func = 1
        };

        /** \copydoc ElementBase::local_degree */
        static constexpr int local_degree = 0;

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
          typedef StandardScalarEvalTraits<EvalPolicy, 1, DataType_> Traits;

        public:
          /// space evaluator type
          typedef Discontinuous::Evaluator<Element, Traits, TrafoEvaluator_, VariantTag> Type;
        };

        /** \copydoc ElementBase::DofMappingType */
        typedef DofMappingSingleEntity<Element, 0> DofMappingType;

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
          typename Functor_,
          int shape_dim_,
          typename DataType_ = Real>
        class NodeFunctional
        {
        private:
          /// dummy enum
          enum
          {
            /// co-dimension
            codim = ShapeType::dimension - shape_dim_,
          };

        public:
          /// node functional type
          typedef Discontinuous::NodeFunctional<Element, Functor_, codim, VariantTag, DataType_> Type;
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
          return this->get_mesh().get_num_entities(ShapeType::dimension);
        }

        /** \copydoc ElementBase::name() */
        static String name()
        {
          return String("Discontinuous<") + Variant_::name() + ">";
        }

      }; // class Element
    } // namespace Discontinuous
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_DISCONTINUOUS_ELEMENT_HPP
