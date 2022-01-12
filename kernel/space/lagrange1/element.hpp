// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SPACE_LAGRANGE1_ELEMENT_HPP
#define KERNEL_SPACE_LAGRANGE1_ELEMENT_HPP 1

// includes, FEAT
#include <kernel/space/element_base.hpp>
#include <kernel/space/dof_assignment_base.hpp>
#include <kernel/space/dof_mapping_common.hpp>
#include <kernel/space/lagrange1/dof_traits.hpp>
#include <kernel/space/lagrange1/evaluator.hpp>
#include <kernel/space/lagrange1/node_functional.hpp>

namespace FEAT
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

          /// number of local dofs := number of vertices per cell
          static constexpr int num_loc_dofs = Shape::FaceTraits<ShapeType, 0>::count;

          /// space evaluation traits
          typedef StandardScalarEvalTraits<EvalPolicy, num_loc_dofs, DataType_> Traits;

        public:
          /// space evaluator type
          typedef Lagrange1::Evaluator<Element, TrafoEvaluator_, Traits> Type;
        };

        /** \copydoc ElementBase::DofMappingType */
        typedef DofMappingSingleEntity<Element, ShapeType::dimension> DofMappingType;

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
          /// node functional type
          typedef Lagrange1::NodeFunctional<Element, shape_dim_, DataType_> Type;
        };

      public:
        /**
         * \brief Constructor
         *
         * \param[in] trafo
         * A \resident reference to the transformation which is to be used by this space.
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

        /** \copydoc ElementBase::name() */
        static String name()
        {
          return "Lagrange1";
        }
      }; // class Element
    } // namespace Lagrange1
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_LAGRANGE1_ELEMENT_HPP
