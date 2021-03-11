// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_TRAFO_ISOSPHERE_MAPPING_HPP
#define KERNEL_TRAFO_ISOSPHERE_MAPPING_HPP 1

// includes, FEAT
#include <kernel/trafo/mapping_base.hpp>
#include <kernel/trafo/isosphere/evaluator.hpp>

namespace FEAT
{
  namespace Trafo
  {
    /**
     * \brief Iso-Parametric Unit-Sphere Transformation namespace
     */
    namespace IsoSphere
    {
      /**
       * \brief Iso-Sphere transformation mapping class template
       *
       * This transformation implements an iso-parametric transformation for the
       * 2D unit-circle or the 3D unit-sphere domain.
       *
       * \tparam Mesh_
       * The mesh class that this transformation is to be defined on.
       * This must be a mesh that represents a simplex surface mesh.
       *
       * \author Peter Zajac
       */
      template<typename Mesh_>
      class Mapping :
        public MappingBase<Mesh_>
      {
      public:
        /// base-class typedef
        typedef MappingBase<Mesh_> BaseClass;
        /// mesh type
        typedef Mesh_ MeshType;
        /// data type
        typedef typename MeshType::VertexSetType::CoordType CoordType;
        /// shape type
        typedef typename MeshType::ShapeType ShapeType;

        // make sure that the mesh is a sub-dimensional mesh
        static_assert(MeshType::shape_dim < MeshType::world_dim, "mesh must be sub-dimensional for IsoSphere trafo mapping");
        // only 1D or 2D Simplex meshes are allowed
        static_assert(MeshType::shape_dim <= 2, "only Simplex<1> and Simplex<2> shapes are allowed (yet)");

      public:
        /** \copydoc MappingBase::Evaluator */
        template<
          typename Shape_ = ShapeType,
          typename CoordType_ = Real>
        class Evaluator
        {
        private:
          /// evaluation policy
          typedef Trafo::StandardEvalPolicy<Shape_, CoordType_, MeshType::world_dim> EvalPolicy;

        public:
          /// evaluator type
          typedef Trafo::IsoSphere::Evaluator<Mapping, EvalPolicy> Type;
        };

      public:
        /**
         * \brief Constructor
         *
         * \param[in] mesh
         * A reference to the mesh that this trafo mapping is to be defined on.
         */
        explicit Mapping(MeshType& mesh) :
          BaseClass(mesh)
        {
        }
      }; // class Mapping<...>
    } // namespace IsoSphere
  } // namespace Trafo
} // namespace FEAT

#endif // KERNEL_TRAFO_ISOSPHERE_MAPPING_HPP
