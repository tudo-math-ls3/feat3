#pragma once
#ifndef KERNEL_SPACE_RANNACHER_TUREK_ELEMENT_HPP
#define KERNEL_SPACE_RANNACHER_TUREK_ELEMENT_HPP 1

// includes, FEAST
#include <kernel/space/element_base.hpp>
#include <kernel/space/dof_assignment_base.hpp>
#include <kernel/space/dof_mapping_common.hpp>
#include <kernel/space/rannacher_turek/dof_traits.hpp>
#include <kernel/space/rannacher_turek/evaluator.hpp>
#include <kernel/space/rannacher_turek/node_functional.hpp>

namespace FEAST
{
  namespace Space
  {
    /**
     * \brief Rannacher-Turek element namespace
     *
     * This namespace encapsulates all classes related to the implementation of the H1-non-conforming
     * Rannacher-Turek Finite Element spaces, also known as Q1~, in several variants.
     *
     * \see <b>R. Rannacher, S. Turek</b>:<i>Simple Nonconforming Quadrilateral Stokes Elements</i>;\n
     * Numerical Methods for Partial Differential Equations, Volume 8 (1992), pp. 97-111
     * \see http://www.mathematik.tu-dortmund.de/lsiii/cms/papers/RannacherTurek1990.pdf
     */
    namespace RannacherTurek
    {
      /**
       * \brief Rannacher-Turek element class template
       *
       * \tparam Trafo_
       * The transformation that is to be used by this space.
       *
       * \tparam Variant_
       * The variant of the element to be used. See RannacherTurek::Variant for details.
       *
       * \author Peter Zajac
       */
      template<
        typename Trafo_,
        typename Variant_ = Variant::StdNonPar>
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
        static constexpr int local_degree = 2;

        /** \copydoc ElementBase::Evaluator */
        template<
          typename TrafoEvaluator_,
          typename DataType_ = typename TrafoEvaluator_::DataType>
        class Evaluator
        {
        private:
          /// evaluation policy
          typedef typename TrafoEvaluator_::EvalPolicy EvalPolicy;

          /// dummy enum
          enum
          {
            /// number of local dofs := number of facets per cell
            num_loc_dofs = Shape::FaceTraits<ShapeType, ShapeType::dimension-1>::count
          };

          /// space evaluation traits
          typedef StandardScalarEvalTraits<EvalPolicy, num_loc_dofs, DataType_> Traits;

        public:
          /// space evaluator type
          typedef RannacherTurek::Evaluator<Element, Traits, TrafoEvaluator_, VariantTag> Type;
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
          typedef DofAssignmentUniform<Element, shape_dim_, DataType_, DofTraits, DofTag<ShapeType, VariantTag> > Type;
        };

        /** \copydoc ElementBase::NodeFunctional */
        template<
          typename Function_,
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
          typedef RannacherTurek::NodeFunctional<Element, Function_, codim, VariantTag, DataType_> Type;
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
      }; // class Element
    } // namespace RannacherTurek
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_RANNACHER_TUREK_ELEMENT_HPP
