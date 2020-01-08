// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_INTERNAL_ENTITY_COUNTER_HPP
#define KERNEL_GEOMETRY_INTERNAL_ENTITY_COUNTER_HPP 1

// includes, FEAT
#include <kernel/shape.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      /**
       * \brief Entity counter class template
       *
       * This is an auxiliary class templates used to compute the number of entities and refinement offsets
       * used by mesh refinement algorithms.
       *
       * \author Peter Zajac
       */
      template<
        template<typename,int> class RefineTraits_,
        typename Shape_,
        int face_dim_,
        int cell_dim_ = Shape_::dimension>
      struct EntityCounter
      {
        static_assert(face_dim_ >= 0, "invalid face dimension");
        static_assert(face_dim_ < cell_dim_, "invalid cell/face dimension");
        static_assert(cell_dim_ <= Shape_::dimension, "invalid cell dimension");

        static String name()
        {
          return "EntityCounter<" + Shape_::name() + "," + stringify(face_dim_) + "," + stringify(cell_dim_) + ">";
        }

        static Index count(const Index num_entities[])
        {
          typedef typename Shape::FaceTraits<Shape_, cell_dim_>::ShapeType CellType;
          return EntityCounter<RefineTraits_, Shape_, face_dim_, cell_dim_-1>::count(num_entities) +
            RefineTraits_<CellType, face_dim_>::count * num_entities[cell_dim_];
        }

        static void offset(Index offsets[], const Index num_entities[])
        {
          EntityCounter<RefineTraits_, Shape_, face_dim_, cell_dim_-1>::offset(offsets, num_entities);
          offsets[cell_dim_] = EntityCounter<RefineTraits_, Shape_, face_dim_, cell_dim_-1>::count(num_entities);
        }
      }; // struct EntityCounter<...>

      template<
        template<typename,int> class RefineTraits_,
        typename Shape_,
        int cell_dim_>
      struct EntityCounter<RefineTraits_, Shape_, cell_dim_, cell_dim_>
      {
        static_assert(cell_dim_ >= 0, "invalid cell/face dimension");
        static_assert(cell_dim_ <= Shape_::dimension, "invalid cell dimension");

        static String name()
        {
          return "EntityCounter<" + Shape_::name() + "," + stringify(cell_dim_) + ">";
        }

        static Index count(const Index num_entities[])
        {
          typedef typename Shape::FaceTraits<Shape_, cell_dim_>::ShapeType CellType;
          return RefineTraits_<CellType, cell_dim_>::count * num_entities[cell_dim_];
        }

        static void offset(Index offsets[], const Index* /*num_entities[]*/)
        {
          offsets[cell_dim_] = 0;
        }
      }; // struct EntityCounter<Shape_, cell_dim_, cell_dim_>

      /**
       * \brief Entity count wrapper class template
       *
       * This is an auxiliary class templates used to compute the number of entities and refinement offsets
       * used by mesh refinement algorithms.
       *
       * \author Peter Zajac
       */
      template<
        template<typename,int> class RefineTraits_,
        typename Shape_,
        int face_dim_ = Shape_::dimension>
      struct EntityCountWrapper
      {
        static_assert(face_dim_ > 0, "invalid face dimension");
        static_assert(face_dim_ <= Shape_::dimension, "invalid face dimension");

        /**
          * \brief Returns the name of the class.
          * \returns
          * The name of the class as a String.
          */
        static String name()
        {
          return "EntityCounterWrapper<" + Shape_::name() + "," + stringify(face_dim_) + ">";
        }

        /**
          * \brief Calculates the number of fine mesh entities.
          *
          * \param[in,out] num_entities
          * An array of length Shape_::dimension containing the number of entities in each dimension.
          * \li On entry, the number of entities of the coarse mesh.
          * \li On exit, the number of entities of the fine mesh.
          */
        static void query(Index num_entities[])
        {
          EntityCountWrapper<RefineTraits_, Shape_, face_dim_-1>::query(num_entities);
          num_entities[face_dim_] = EntityCounter<RefineTraits_, Shape_, face_dim_>::count(num_entities);
        }
      }; // struct EntityCounter<...>

      template<
        template<typename,int> class RefineTraits_,
        typename Shape_>
      struct EntityCountWrapper<RefineTraits_, Shape_, 0>
      {
        static String name()
        {
          return "EntityCounterWrapper<" + Shape_::name() + ",0>";
        }

        static void query(Index num_entities[])
        {
          num_entities[0] = EntityCounter<RefineTraits_, Shape_, 0>::count(num_entities);
        }
      }; // struct EntityCountWrapper<Shape_, 0>
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_INTERNAL_ENTITY_COUNTER_HPP
