// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/space/element_base.hpp>
#include <kernel/space/dof_mapping_common.hpp>
#include <kernel/space/dof_assignment_common.hpp>
#include <kernel/space/argyris/dof_traits.hpp>
#include <kernel/space/argyris/evaluator.hpp>
#include <kernel/space/argyris/node_functional.hpp>

namespace FEAT
{
  namespace Space
  {
    /**
     * \brief Argyris Element namespace
     *
     * This namespace encapsulates all classes related to the implementation of the sixth-order
     * H2-conforming triangular Finite Element space known as the "Argyris" element.
     *
     * \see P.G. Ciarlet: The Finite Element Method for Elliptic Problems; pp. 71-73
     */
    namespace Argyris
    {
      /**
       * \brief Argyris element class template
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

        static constexpr int num_vert_dofs = 6;
        static constexpr int num_edge_dofs = 1;

        /// node functionals available
        static constexpr bool have_node_func = true;

        /** \copydoc ElementBase::local_degree */
        static constexpr int local_degree = 5;

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
          static constexpr int num_loc_dofs = 21;

          /// space evaluation traits
          typedef StandardScalarEvalTraits<EvalPolicy, num_loc_dofs, DataType_> Traits;

        public:
          /// space evaluator type
          typedef Argyris::Evaluator<Element, TrafoEvaluator_, Traits> Type;
        };

        /** \copydoc ElementBase::DofMappingType */
        typedef DofMappingUniform<Element, Argyris::DofTraits, ShapeType> DofMappingType;

        /** \copydoc ElementBase::DofAssignment */
        template<
          int shape_dim_,
          typename DataType_ = Real>
        class DofAssignment
        {
        public:
          /// Dof-Assignment type
          typedef DofAssignmentUniform<Element, shape_dim_, DataType_, Argyris::DofTraits, ShapeType> Type;
        };

        /** \copydoc ElementBase::NodeFunctional */
        template<
          int shape_dim_,
          typename DataType_ = Real>
        class NodeFunctional
        {
        public:
          typedef Argyris::NodeFunctional<Element, shape_dim_, DataType_> Type;
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
          return
            // number of vertex dofs
            Index(num_vert_dofs) * this->get_mesh().get_num_entities(0) +
            // number of edge dofs
            Index(num_edge_dofs) * this->get_mesh().get_num_entities(1);
        }

        /** \copydoc ElementBase::name() */
        static String name()
        {
          return "Argyris";
        }
      }; // class Element<...>
    } // namespace Argyris
  } // namespace Space
} // namespace FEAT
