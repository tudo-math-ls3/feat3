// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SPACE_LAGRANGE3_ELEMENT_HPP
#define KERNEL_SPACE_LAGRANGE3_ELEMENT_HPP 1

// includes, FEAT
#include <kernel/space/element_base.hpp>
#include <kernel/space/dof_assignment_base.hpp>
#include <kernel/space/dof_mapping_common.hpp>
#include <kernel/space/lagrange3/dof_traits.hpp>
#include <kernel/space/lagrange3/evaluator.hpp>
#include <kernel/space/lagrange3/node_functional.hpp>

namespace FEAT
{
  namespace Space
  {
    /**
     * \brief Lagrange-3 Element namespace
     *
     * This namespace encapsulates all classes related to the implementation of the standard fourth-order
     * H1-conforming Finite Element spaces widely known as P3 and Q3, resp.
     */
    namespace Lagrange3
    {
      /**
       * \brief Standard Lagrange-3 Finite-Element space class template
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
        static constexpr bool have_node_func = 1;

        /** \copydoc ElementBase::local_degree */
        static constexpr int local_degree = 3;

        /** \copydoc ElementBase::DofMappingType */
        typedef DofMappingUniform<Element, Lagrange3::DofTraits, ShapeType> DofMappingType;

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
          static constexpr int num_loc_dofs = DofMappingType::dof_count;

          /// space evaluation traits
          typedef StandardScalarEvalTraits<EvalPolicy, num_loc_dofs, DataType_> Traits;

        public:
          /// space evaluator type
          typedef Lagrange3::Evaluator<Element, TrafoEvaluator_, Traits> Type;
        };

        /** \copydoc ElementBase::DofAssignment */
        template<
          int shape_dim_,
          typename DataType_ = Real>
        class DofAssignment
        {
        public:
          /// Dof-Assignment type
          typedef DofAssignmentUniform<Element, shape_dim_, DataType_, Lagrange3::DofTraits, ShapeType> Type;
        };

        /** \copydoc ElementBase::NodeFunctional */
        template<
          int shape_dim_,
          typename DataType_ = Real>
        class NodeFunctional
        {
          typedef typename Shape::FaceTraits<ShapeType, shape_dim_>::ShapeType FaceType;
        public:
          /// node functional type
          typedef Lagrange3::NodeFunctional<Element, FaceType, DataType_> Type;
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
          DofMappingType dof_map(*this);
          return dof_map.get_num_global_dofs();
        }

        /** \copydoc ElementBase::name() */
        static String name()
        {
          return "Lagrange3";
        }
      }; // class Element
    } // namespace Lagrange3
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_LAGRANGE3_ELEMENT_HPP
