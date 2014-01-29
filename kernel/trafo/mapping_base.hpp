#pragma once
#ifndef KERNEL_TRAFO_MAPPING_BASE_HPP
#define KERNEL_TRAFO_MAPPING_BASE_HPP 1

// includes, FEAST
#include <kernel/trafo/evaluator_base.hpp>

namespace FEAST
{
  namespace Trafo
  {
    /**
     * \brief Trafo-Mapping base class
     *
     * This class acts as a base class for transformation mapping implementations.
     *
     * \author Peter Zajac
     */
    template<typename Mesh_>
    class MappingBase
    {
    public:
      /// mesh type
      typedef Mesh_ MeshType;
      /// shape type
      typedef typename MeshType::ShapeType ShapeType;

      // Note:
      // The following block serves as an element interface documentation and is therefore only
      // visible to doxygen. The actual functionality has to be supplied by the implementation.
#ifdef DOXYGEN
      /**
       * \brief Trafo evaluator class template.
       *
       * \tparam Shape_
       * The shape for which the evaluator is to be defined. Must be a face of #ShapeType.
       *
       * \tparam DataType_
       * The data type that is to be used by the evaluator.
       */
      template<
        typename Shape_ = ShapeType,
        typename DataType_ = Real>
      class Evaluator
      {
      public:
        /// evaluator type
        typedef ... Type;
      };
#endif // DOXYGEN

    private:
      /// deleted copy-constructor
      MappingBase(const MappingBase&);
      /// deleted assignment-operator
      MappingBase& operator=(const MappingBase&);

    protected:
      /// mesh reference
      MeshType& _mesh;

      /**
       * \brief Constructor
       *
       * \param[in] mesh
       * A reference to the mesh that this trafo mapping is to be defined on.
       *
       * \note This constructor is protected so that it can only be called from a derived class.
       */
      explicit MappingBase(MeshType& mesh) :
        _mesh(mesh)
      {
      }

    public:
      /**
       * \brief Returns a reference to the underlying mesh.
       * \returns
       * A (const) reference to the underlying mesh.
       */
      MeshType& get_mesh()
      {
        return _mesh;
      }

      /** \copydoc get_mesh() */
      const MeshType& get_mesh() const
      {
        return _mesh;
      }

      /**
       * \brief Comparison operator
       *
       * This operator checks whether two trafo mapping objects are equal.
       *
       * \returns
       * \c true, if this \c this and \p other name the same object, otherwise \c false.
       */
      bool operator==(const MappingBase& other) const
      {
        return this == &other;
      }
    }; // class MappingBase<...>
  } // namespace Trafo
} // namespace FEAST

#endif // KERNEL_TRAFO_MAPPING_BASE_HPP
