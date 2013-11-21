#pragma once
#ifndef KERNEL_SPACE_HERMITE3_ELEMENT_HPP
#define KERNEL_SPACE_HERMITE3_ELEMENT_HPP 1

// includes, FEAST
#include <kernel/space/element_base.hpp>
#include <kernel/space/dof_assignment_base.hpp>
#include <kernel/space/dof_mapping_common.hpp>
#include <kernel/space/hermite3/dof_traits.hpp>
#include <kernel/space/hermite3/evaluator.hpp>
#include <kernel/space/hermite3/node_functional.hpp>

namespace FEAST
{
  namespace Space
  {
    /**
     * \brief Hermite-3 element namespace
     *
     * This namespace encapsulates all classes related to the implementation of the "almost H2-conforming"
     * Hermite-3 Finite Element spaces.
     */
    namespace Hermite3
    {
      /**
       * \brief Hermite-3 element class template
       *
       * \attention
       * This element is experimental. Do not use this element unless you know exactly what you are doing!
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

        /// number of dofs per vertex
        static constexpr int num_vert_dofs = DofTraits<ShapeType, 0>::count;
        /// number of dofs per cell
        static constexpr int num_cell_dofs = DofTraits<ShapeType, ShapeType::dimension>::count;

        /// node functionals available
        static constexpr bool have_node_func = true;

        /** \copydoc ElementBase::local_degree */
        static constexpr int local_degree = 3;

        /** \copydoc ElementBase::Evaluator */
        template<
          typename TrafoEvaluator_,
          typename DataType_ = typename TrafoEvaluator_::DataType>
        class Evaluator
        {
        private:
          /// evaluation policy
          typedef typename TrafoEvaluator_::EvalPolicy EvalPolicy;

          /// number of local dofs
          static constexpr int num_loc_dofs = Shape::FaceTraits<ShapeType, 0>::count * num_vert_dofs + num_cell_dofs;

          /// space evaluation traits
          typedef StandardScalarEvalTraits<EvalPolicy, num_loc_dofs, DataType_> Traits;

        public:
          /// space evaluator type
          typedef Hermite3::Evaluator<Element, TrafoEvaluator_, Traits> Type;
        };

        /** \copydoc ElementBase::DofMappingType */
        typedef DofMappingUniform<Element, DofTraits, ShapeType> DofMappingType;

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
        public:
          typedef Hermite3::NodeFunctional<Element, ShapeType, shape_dim_, DataType_> Type;
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
          return
            // number of vertex dofs
            Index(num_vert_dofs) * this->get_mesh().get_num_entities(0) +
            // number of cell dofs
            Index(num_cell_dofs) * this->get_mesh().get_num_entities(ShapeType::dimension);
        }

        /** \copydoc ElementBase::name() */
        static String name()
        {
          return "Hermite3";
        }
      }; // class Element<...>
    } // namespace Hermite3
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_HERMITE3_ELEMENT_HPP
