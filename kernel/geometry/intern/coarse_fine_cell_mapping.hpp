// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_INTERN_COARSE_FINE_CELL_MAPPING_HPP
#define KERNEL_GEOMETRY_INTERN_COARSE_FINE_CELL_MAPPING_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/adjacency/adjactor.hpp>

namespace FEAT
{
  namespace Geometry
  {
    namespace Intern
    {
#ifdef DOXYGEN
      /**
       * \brief A coarse cell to fine cell mapping class.
       *
       * \tparam FineMesh_ The fine mesh type to be used.
       * \tparam CoarseMesh_ The coarse mesh type to be used.
       * \tparam is_structured_ Classification of the mesh types.
       *
       * This class provides the mapping of coarse cells to fine cells.
       *
       * \attention
       * This class silenty assumes that \p fine_mesh is the mesh that was obtained
       * by applying the standard 2-level refinement algorithm on the \p coarse_mesh.
       *
       * \author Heiko Poelstra
       */
      template<typename FineMesh_, typename CoarseMesh_=FineMesh_, bool is_structured_=FineMesh_::is_structured>
      class CoarseFineCellMapping
      {
      public:
        /**
         * \brief Constructor.
         *
         * \param[in] fine_mesh
         * The fine mesh.
         *
         * \param[in] coarse_mesh
         * The coarse mesh.
         */
        explicit CoarseFineCellMapping(const FineMesh_& fine_mesh, const CoarseMesh_& coarse_mesh);

        /**
         * \brief Calculates the fine mesh cell index.
         *
         * \param[in] ccell
         * The coarse mesh cell index.
         *
         * \param[in] child
         * The child cell index of the coarse mesh cell.
         *
         * \returns
         * The fine mesh cell index.
         */
        Index calc_fcell(Index ccell, Index child=0) const;

        /**
         * \brief Returns the number of fine mesh cells per coarse mesh cell.
         *
         * \returns
         * The number of fine mesh cells per coarse mesh cell.
         */
        Index get_num_children() const;

        /** \copydoc Adjacency::Adjactor::ImageIterator */
        typedef Adjacency::Adjactor::ImageIterator ImageIterator;

        /** \copydoc Adjacency::Adjactor::get_num_nodes_domain() */
        Index get_num_nodes_domain() const;

        /** \copydoc Adjacency::Adjactor::get_num_nodes_image() */
        Index get_num_nodes_image() const;

        /** \copydoc Adjacency::Adjactor::image_begin() */
        ImageIterator image_begin(Index domain_node) const;

        /** \copydoc Adjacency::Adjactor::image_end() */
        ImageIterator image_end(Index domain_node) const;
      }; // class CoarseFineCellMapping<FineMesh_, CoarseMesh_, false>
#else
      template<typename FineMesh_, typename CoarseMesh_=FineMesh_, bool is_structured=FineMesh_::is_structured>
      class CoarseFineCellMapping;
#endif // DOXYGEN
    } // namespace Intern

    /// \cond internal
    namespace Intern
    {
      // Mapping for structured meshes
      template<typename FineMesh_, typename CoarseMesh_>
      class CoarseFineCellMapping<FineMesh_, CoarseMesh_, true>
      {
        static_assert(FineMesh_::is_structured == true, "Do not provide the third template parameter manually!");

        Index _num_elements;
        Index _num_children;
        Index _dim[3]={1,1,1}; // fine mesh dimension

      public:
        explicit CoarseFineCellMapping(const FineMesh_& fine_mesh, const CoarseMesh_& coarse_mesh) :
          _num_elements(coarse_mesh.get_num_entities(CoarseMesh_::ShapeType::dimension)),
          _num_children(fine_mesh.get_num_entities(FineMesh_::ShapeType::dimension) / _num_elements)
        {
          for (int d(0); d < FineMesh_::ShapeType::dimension; ++d)
          {
            _dim[d] = fine_mesh.get_num_slices(d);
          }
        }

        Index calc_fcell(Index ccell, Index child=0) const
        {
          // extract cell position based on structured numbering
          Index xy_plane_size = _dim[0]/2;
          if (FineMesh_::shape_dim > 1)
          {
            xy_plane_size *= _dim[1]/2;
          }
          Index zpos = ccell / xy_plane_size; // yes, it is ok to point to a non-existent slice.
                                              // see also the implementation of the operator++ of ImageIterator
          Index ypos = (ccell - zpos * xy_plane_size) / (_dim[0]/2);
          Index xpos = ccell - zpos * xy_plane_size - ypos * (_dim[0]/2);

          // coarse position -> fine position
          zpos = 2*zpos + (child / 4);
          ypos = 2*ypos;
          if (child == 2 || child == 3 || child == 6 || child == 7) { ++ypos; }
          xpos = 2*xpos + (child % 2);

          return zpos * _dim[1] * _dim[0] + ypos * _dim[0] + xpos;
        }

        Index get_num_children() const
        {
          return _num_children;
        }

        class ImageIterator
        {
        protected:
          // the current index of the iterator
          Index _index;
          // the dimensions of the fine mesh
          Index _dim[3]={1,1,1};

        public:
          ImageIterator() :
            _index(0)
          {
          }

          explicit ImageIterator(Index index, const Index dim[]) :
            _index(index)
          {
            _dim[0] = dim[0];
            _dim[1] = dim[1];
            _dim[2] = dim[2];
          }

          ImageIterator(const ImageIterator& other) :
            _index(other._index),
            _dim{other._dim[0], other._dim[1], other._dim[2]}
          {
          }

          ImageIterator& operator=(const ImageIterator& other)
          {
            _index = other._index;
            _dim[0] = other._dim[0];
            _dim[1] = other._dim[1];
            _dim[2] = other._dim[2];
            return *this;
          }

          ImageIterator& operator++()
          {
            // extract cell position based on structured numbering
            Index zpos = _index / (_dim[0] * _dim[1]);
            Index ypos = (_index - zpos * _dim[0] * _dim[1]) / _dim[0];
            Index xpos = _index - zpos * _dim[0] * _dim[1] - ypos * _dim[0];

            // update cell position
            if (xpos % 2 == 0)
            {
              // in-element jump
              ++xpos;
            }
            else
            {
              --xpos;
              if (FineMesh_::ShapeType::dimension > 1 && ypos % 2 == 0)
              {
                // still in-element jump
                ++ypos;
              }
              else
              {
                if (FineMesh_::ShapeType::dimension > 1) --ypos;
                if (FineMesh_::ShapeType::dimension > 2 && zpos % 2 == 0)
                {
                  // still in-element jump
                  ++zpos;
                }
                else
                {
                  // jump to next element
                  if (FineMesh_::ShapeType::dimension > 2) --zpos; // now we point to the first fine element of the coarse element
                  if (xpos + 2 < _dim[0])
                  {
                    xpos += 2;
                  }
                  else if (ypos + 2 < _dim[1])
                  {
                    xpos = 0;
                    ypos += 2;
                  }
                  else // dont check for limit. this means we can also jump into a non-existent slice.
                       // that's needed to have the iterator at the correct endpoint.
                       // see also the implementation of calc_fcell
                  {
                    xpos = 0;
                    ypos = 0;
                    zpos += 2;
                  }
                }
              }
            }
            _index = zpos * _dim[0] * _dim[1] + ypos * _dim[0] + xpos;
            return *this;
          }

          Index operator*() const
          {
            return _index;
          }

          bool operator!=(const ImageIterator& other) const
          {
            return _index != other._index;
          }
        }; // class CoarseFineCellMapping<FineMesh_, CoarseMesh_, true>::ImageIterator

        Index get_num_nodes_domain() const
        {
          return _num_elements;
        }

        Index get_num_nodes_image() const
        {
          return _num_elements * _num_children;
        }

        ImageIterator image_begin(Index domain_node) const
        {
          return ImageIterator(calc_fcell(domain_node), _dim);
        }

        ImageIterator image_end(Index domain_node) const
        {
          return ImageIterator(calc_fcell(domain_node+1), _dim);
        }
      }; // class CoarseFineCellMapping<FineMesh_, CoarseMesh_, true>

      // Mapping for unstructured meshes
      template<typename FineMesh_, typename CoarseMesh_>
      class CoarseFineCellMapping<FineMesh_, CoarseMesh_, false>
      {
        static_assert(FineMesh_::is_structured == false, "Do not provide the third template parameter manually!");

        Index _num_elements;
        Index _num_children;

      public:
        explicit CoarseFineCellMapping(const FineMesh_& fine_mesh, const CoarseMesh_& coarse_mesh) :
          _num_elements(coarse_mesh.get_num_entities(CoarseMesh_::ShapeType::dimension)),
          _num_children(fine_mesh.get_num_entities(FineMesh_::ShapeType::dimension) / _num_elements)
        {
        }

        Index calc_fcell(Index ccell, Index child=0) const
        {
          return ccell*_num_children + child;
        }

        Index get_num_children() const
        {
          return _num_children;
        }

        typedef Adjacency::Adjactor::IndexImageIterator ImageIterator;

        Index get_num_nodes_domain() const
        {
          return _num_elements;
        }

        Index get_num_nodes_image() const
        {
          return _num_elements * _num_children;
        }

        ImageIterator image_begin(Index domain_node) const
        {
          return ImageIterator(calc_fcell(domain_node));
        }

        ImageIterator image_end(Index domain_node) const
        {
          return ImageIterator(calc_fcell(domain_node+1));
        }
      }; // class CoarseFineCellMapping<FineMesh_, CoarseMesh_, false>

    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_INTERN_COARSE_FINE_CELL_MAPPING_HPP
