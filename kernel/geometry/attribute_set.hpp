#pragma once
#ifndef KERNEL_GEOMETRY_ATTRIBUTE_SET_HPP
#define KERNEL_GEOMETRY_ATTRIBUTE_SET_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>

// includes, system
#include <vector>

namespace FEAT
{
  namespace Geometry
  {
    /**
     * \brief Container for saving data related to mesh entities
     *
     * This class template implements a container for saving attributes of vertices and is used in
     * MeshPart. The attributes do not need to be scalar.
     *
     * \tparam DataType_
     * The datatype for aan attribute value entry.
     *
     * \author Peter Zajac
     */
    template<typename DataType_>
    class AttributeSet
    {
      public:
        /// Type for the values in the attributes;
        /// This usually coincides with the CoordType of the mesh
        typedef DataType_ DataType;

      protected:
        /// Number of attribute values
        Index _num_values;
        /// Number of entries per attribute value
        int _dimension;
        /// Value array
        std::vector<DataType> _values;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] num_values
         * The number of values to be allocated.
         *
         * \param[in] dimension
         * The dimension of the attribute. Must be > 0.
         */
        explicit AttributeSet(Index num_values, int dimension = 1) :
          _num_values(num_values),
          _dimension(dimension)
        {
          XASSERT((dimension > 0) || (num_values == Index(0)));
          if(num_values > Index(0))
          {
            _values.resize(std::size_t(num_values) * std::size_t(dimension));
          }
        }

        /// virtual destructor
        virtual ~AttributeSet()
        {
        }

        /// \returns The size of dynamically allocated memory in bytes.
        std::size_t bytes() const
        {
          return _values.size() * std::size_t(_dimension) * sizeof(DataType);
        }

        /**
         * \brief Returns the number attribute dimension.
         *
         * \warning This function may return 0.
         */
        int get_dimension() const
        {
          return _dimension;
        }

        /// Returns the number of attribute values
        Index get_num_values() const
        {
          return _num_values;
        }

        /**
         * \brief Returns a reference to an attribute value entry.
         *
         * \param[in] i
         * The index of the attribute value to be returned.
         *
         * \param[in] j
         * The index of the attribute value component to be returned.
         *
         * \returns
         * A reference to j-th component of the i-th attribute value.
         */
        DataType& operator()(Index i, int j)
        {
          ASSERT(i < _num_values);
          ASSERT((j >= 0) && (j < _dimension));
          return _values[std::size_t(i)*std::size_t(_dimension) + std::size_t(j)];
        }

        /**
         * \brief Returns a const reference to an attribute value entry.
         *
         * \param[in] i
         * The index of the attribute value to be returned.
         *
         * \param[in] j
         * The index of the attribute value component to be returned.
         *
         * \returns
         * A const reference to j-th component of the i-th attribute value.
         */
        const DataType& operator()(Index i, int j) const
        {
          ASSERT(i < _num_values);
          ASSERT((j >= 0) && (j < _dimension));
          return _values[std::size_t(i)*std::size_t(_dimension) + std::size_t(j)];
        }

        /// \cond internal
        // interpret attribute i as raw array and return pointer
        DataType* raw_at(Index i)
        {
          ASSERT(i < _num_values);
          return &_values.data()[std::size_t(i)*std::size_t(_dimension)];
        }

        const DataType* raw_at(Index i) const
        {
          ASSERT(i < _num_values);
          return &_values.data()[std::size_t(i)*std::size_t(_dimension)];
        }
        /// \endcond

        /// Returns the name of the class.
        static String name()
        {
          return "AttributeSet<...>";
        }
    }; // class AttributeSet<...>
  } // namespace Geometry
} // namespace FEAT
#endif // KERNEL_GEOMETRY_ATTRIBUTE_SET_HPP
