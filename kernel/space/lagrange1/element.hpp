#pragma once
#ifndef KERNEL_SPACE_LAGRANGE1_ELEMENT_HPP
#define KERNEL_SPACE_LAGRANGE1_ELEMENT_HPP 1

// includes, FEAST
#include <kernel/space/element_base.hpp>
#include <kernel/space/dof_mapping_common.hpp>
#include <kernel/space/lagrange1/evaluator.hpp>

namespace FEAST
{
  namespace Space
  {
    /**
     * \brief Lagrange-1 Element namespace
     *
     * This namespace encapsulates all classes related to the implementation of the standard second-order
     * H1-conforming Finite Element spaces widely known as P1 and Q1, resp.
     */
    namespace Lagrange1
    {
      /**
       * \brief Standard Lagrange-1 Finite-Element space class template
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

        /** \copydoc ElementBase::TrafoConfig */
        template<typename SpaceConfig_>
        struct TrafoConfig :
          public Trafo::ConfigBase
        {
          /** \copydoc Trafo::ConfigBase::TrafoRequirements */
          enum
          {
            /// we always need domain coordinates
            need_dom_point = 1,
            /// we need jacobian inverses for gradient calculation
            need_jac_inv = SpaceConfig_::need_grad
          };
        };

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
            /// number of local dofs := number of vertices per cell
            num_loc_dofs = Shape::FaceTraits<ShapeType, 0>::count
          };

          /// space evaluation traits
          typedef StandardScalarEvalTraits<EvalPolicy, num_loc_dofs, DataType_> Traits;

        public:
          /// space evaluator type
          typedef Lagrange1::Evaluator<Element, TrafoEvaluator_, Traits> Type;
        };

        /** \copydoc ElementBase::DofMapping */
        template<int shape_dim_ = ShapeType::dimension>
        class DofMapping
        {
        public:
          /// Dof-Mapping type
          typedef DofMappingSingleEntity<Element, shape_dim_, 0> Type;
        }; // class DofMapper

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
          // number of DOFs = number of vertices in the mesh
          return this->get_mesh().get_num_entities(0);
        }
      }; // class Element
    } // namespace Lagrange1
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_LAGRANGE1_ELEMENT_HPP
