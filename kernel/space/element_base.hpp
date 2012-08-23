#pragma once
#ifndef KERNEL_SPACE_ELEMENT_BASE_HPP
#define KERNEL_SPACE_ELEMENT_BASE_HPP 1

// includes, FEAST
#include <kernel/space/dof_mapping_base.hpp>
#include <kernel/space/evaluator_base.hpp>

namespace FEAST
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

      // Note:
      // The following block serves as an element interface documentation and is therefore only
      // visible to doxygen. The actual functionality has to be supplied by the implementation.
#ifdef DOXYGEN
      /**
       * \brief Trafo Config class template
       *
       * This class template specifies the trafo-config of the element's evaluator for a given
       * space-config.
       *
       * \tparam SpaceConfig_
       * The space configuration for which the trafo configuration is to be set up.
       */
      template<typename SpaceConfig_>
      struct TrafoConfig :
        public Trafo::ConfigBase
      {
        /** \copydoc Trafo::ConfigBase::TrafoRequirements */
        enum
        {
          ...
        };
      };

      /**
       * \brief Space evaluator class template.
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
       * \brief Dof-Mapping class template.
       *
       * \tparam shape_dim_
       * The shape dimension for which the Dof-Mapping shall be defined.
       */
      template<int shape_dim_ = ShapeType::dimension>
      class DofMapping
      {
      public:
        /// Dof-Mapping type
        typedef ... Type;
      }; // class DofMapper

      /**
       * \brief Returns the number of dofs.
       *
       * \returns
       * The total number of global degrees of freedom for the patch.
       */
      Index get_num_dofs() const;
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
} // namespace FEAST

#endif // KERNEL_SPACE_ELEMENT_BASE_HPP
