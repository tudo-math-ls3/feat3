// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SPACE_ELEMENT_BASE_HPP
#define KERNEL_SPACE_ELEMENT_BASE_HPP 1

// includes, FEAT
#include <kernel/space/dof_mapping_base.hpp>
#include <kernel/space/evaluator_base.hpp>

namespace FEAT
{
  namespace Space
  {
    /**
     * \brief Finite-Element base class
     *
     * This class acts as a base class for Finite-Element container implementations.
     *
     * \author Peter Zajac
     */
    template<typename Trafo_>
    class ElementBase
    {
    public:
      /// transformation type
      typedef Trafo_ TrafoType;
      /// mesh type
      typedef typename TrafoType::MeshType MeshType;
      /// shape type
      typedef typename TrafoType::ShapeType ShapeType;

      /// shape dimension
      static constexpr int shape_dim = ShapeType::dimension;
      /// our image/world dimension
      static constexpr int world_dim = MeshType::world_dim;

      /*
       * \brief Specifies whether the element has node functionals
       *
       * If this value is 1, the element implements node functionals, i.e. the Type member of the
       * nested NodeFunctional class template is defined.
       */
      static constexpr bool have_node_func = false;

      // Note:
      // The following block serves as an element interface documentation and is therefore only
      // visible to doxygen. The actual functionality has to be supplied by the implementation.
#ifdef DOXYGEN
      /**
       * \brief Space evaluator class wrapper template.
       *
       * \tparam TrafoEvaluator_
       * The trafo evaluator to be used by this space evaluator.
       *
       * \tparam DataType_
       * The data type that is to be used by the evaluator.\n
       * If not given, the data type of the trafo evaluator is used.
       */
      template<
        typename TrafoEvaluator_,
        typename DataType_ = typename TrafoEvaluator_::DataType>
      class Evaluator
      {
      public:
        /// space evaluator type
        typedef ... Type;
      };

      /**
       * \brief Local Polynomial degree count
       *
       * This constant describes the maximum local polynomial degree of the element.
       * This may be used, e.g. for automated determination of appropriate cubature rules.
       *
       * Please note that for Simplex elements, 'degree' refers to the P_k spaces, whereas for
       * Hypercube elements, 'degree' refers to the Q_k spaces.
       */
      static constexpr int local_degree = ...;

      /**
       * \brief Dof-Mapping class.
       */
      typedef ... DofMappingType;

      /**
       * \brief Dof-Assignment class wrapper template
       *
       * \tparam shape_dim_
       * The dimension of the shape for which the Dof-Assignment is to be defined.
       *
       * \tparam DataType_
       * The data-type that is to be used for evaluation.
       */
      template<
        int shape_dim_,
        typename DataType_ = Real>
      class DofAssignment
      {
      public:
        /// dof-assignment type
        typedef ... Type;
      };

      /**
       * \brief Node-Functional class wrapper template
       *
       * \tparam Function_
       * The class of the function that is to be evaluated by the node functional.
       * Must implement the Assembly::AnaylticFunction interface.
       *
       * \tparam shape_dim_
       * The dimension of the shape for which the node-functional is to be defined.
       *
       * \tparam DataType_
       * The data-type that is to be used for evaluation.
       */
      template<
        int shape_dim_,
        typename DataType_ = Real>
      class NodeFunctional
      {
      public:
        /// node-functional type
        typedef ... Type;
      };

      /**
       * \brief Returns the number of dofs.
       *
       * \returns
       * The total number of global degrees of freedom for the patch.
       */
      Index get_num_dofs() const;

      /**
       * \brief Returns the name of the element.
       */
      String name() const;
#endif // DOXYGEN

    private:
      /// deleted copy-constructor
      ElementBase(const ElementBase&);
      /// deleted assignment-operator
      ElementBase& operator=(const ElementBase&);

    protected:
      /// transformation reference
      TrafoType& _trafo;

      /**
       * \brief Constructor
       *
       * \param[in] trafo
       * A reference to the transformation which is to be used by this space.
       *
       * \note This constructor is protected so that it can only be called from a derived class.
       */
      explicit ElementBase(TrafoType& trafo)
        : _trafo(trafo)
      {
      }

    public:
      /**
       * \brief Returns a reference to the trafo.
       *
       * \returns
       * A (const) reference to the transformation used by this space.
       */
      TrafoType& get_trafo()
      {
        return _trafo;
      }

      /** \copydoc get_trafo() */
      const TrafoType& get_trafo() const
      {
        return _trafo;
      }

      /**
       * \brief Returns a reference to the underlying mesh.
       *
       * \returns
       * A (const) reference to the mesh used by this space.
       */
      MeshType& get_mesh()
      {
        return get_trafo().get_mesh();
      }

      /** \copydoc get_mesh() */
      const MeshType& get_mesh() const
      {
        return get_trafo().get_mesh();
      }

      /**
       * \brief Comparison operator
       *
       * This operator checks whether two element objects are equal.
       *
       * \returns
       * \c true, if this \c this and \p other name the same object, otherwise \c false.
       */
      bool operator==(const ElementBase& other) const
      {
        return this == &other;
      }
    }; // class ElementBase<...>
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_ELEMENT_BASE_HPP
