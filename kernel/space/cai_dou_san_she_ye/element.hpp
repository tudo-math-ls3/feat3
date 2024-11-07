// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/space/element_base.hpp>
#include <kernel/space/dof_assignment_base.hpp>
#include <kernel/space/dof_mapping_common.hpp>
#include <kernel/space/cai_dou_san_she_ye/dof_traits.hpp>
#include <kernel/space/cai_dou_san_she_ye/evaluator.hpp>
#include <kernel/space/cai_dou_san_she_ye/node_functional.hpp>

namespace FEAT
{
  namespace Space
  {
    /**
     * \brief Cai-Douglas-Santos-Sheen-Ye element namespace
     */
    namespace CaiDouSanSheYe
    {
      /**
       * \brief Cai-Douglas-Santos-Sheen-Ye element class template
       *
       * This class template implements the second-order H1-non-conforming parametric finite element
       * space known as the Cai-Douglas-Santos-Sheen-Ye element, which is a modified and extended
       * version of the Rannacher-Turek element.
       * Currently, only the 2D quadrilateral element is implemented
       *
       * \see \cite DSSY98, \cite CDSSY00
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

          /// number of local dofs := number of facets per cell
          static constexpr int num_loc_dofs = Shape::FaceTraits<ShapeType, ShapeType::dimension-1>::count + 1;

          /// space evaluation traits
          typedef StandardScalarEvalTraits<EvalPolicy, num_loc_dofs, DataType_> Traits;

        public:
          /// space evaluator type
          typedef CaiDouSanSheYe::Evaluator<Element, TrafoEvaluator_, Traits> Type;
        };

        /** \copydoc ElementBase::DofMappingType */
        typedef DofMappingUniform<Element, CaiDouSanSheYe::DofTraits, ShapeType> DofMappingType;

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
          typedef CaiDouSanSheYe::NodeFunctional<Element, ShapeType, codim, DataType_> Type;
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
          // number of DOFs = number of facets in the mesh + number of elements
          return this->get_mesh().get_num_entities(ShapeType::dimension - 1) + this->get_mesh().get_num_elements();
        }

        /** \copydoc ElementBase::name() */
        static String name()
        {
          return "CaiDouSanSheYe";
        }
      }; // class Element
    } // namespace CaiDouSanSheYe
  } // namespace Space
} // namespace FEAT
