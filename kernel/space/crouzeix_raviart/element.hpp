#pragma once
#ifndef KERNEL_SPACE_CROUZEIX_RAVIART_ELEMENT_HPP
#define KERNEL_SPACE_CROUZEIX_RAVIART_ELEMENT_HPP 1

// includes, FEAST
#include <kernel/space/element_base.hpp>
#include <kernel/space/dof_assignment_base.hpp>
#include <kernel/space/dof_mapping_common.hpp>
#include <kernel/space/crouzeix_raviart/dof_traits.hpp>
#include <kernel/space/crouzeix_raviart/evaluator.hpp>
#include <kernel/space/crouzeix_raviart/node_functional.hpp>

namespace FEAST
{
  namespace Space
  {
    /**
     * \brief Crouzeix-Raviart Element namespace
     *
     * This namespace encapsulates all classes related to the implementation of the second-order
     * H1-nonconforming Crouzeix-Raviart Finite Element space.
     *
     * \see \cite CR73
     */
    namespace CrouzeixRaviart
    {
      /**
       * \brief Crouzeix-Raviart Finite-Element space class template
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
        static constexpr int local_degree = 1;

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
          typedef CrouzeixRaviart::Evaluator<Element, TrafoEvaluator_, Traits> Type;
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
          typedef CrouzeixRaviart::NodeFunctional<Element, codim, DataType_> Type;
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
          return this->get_mesh().get_num_entities(ShapeType::dimension-1);
        }

        /** \copydoc ElementBase::name() */
        static String name()
        {
          return "CrouzeixRaviart";
        }
      }; // class Element
    } // namespace CrouzeixRaviart
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_CROUZEIX_RAVIART_ELEMENT_HPP
