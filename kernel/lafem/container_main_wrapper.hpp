// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_CONTAINER_MAIN_WRAPPER_HPP
#define KERNEL_LAFEM_CONTAINER_MAIN_WRAPPER_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/memory_pool.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/base.hpp>

namespace FEAT
{
  /**
   * \brief LAFEM namespace
   */
  namespace LAFEM
  {
    /**
     *  \brief This class simplifies operations in main memory on LAFEM containers, regardless of their actual memory type.
     *  Data will be transferred to main memory, if necessary, and copied back after the wrapper object is destroyed.
     *
     *  \note Modifications to the original container made during the lifetime of the wrapper object will result in undefined behaviour.
     */
    template <typename CT_>
    class ContainerMainWrapper
    {
      private:
      /// the corresponding container object in remove memory
      CT_ & _parent;

      /// copy mode in constructor (full copy or index arrays only and allocate new data array on host)
      bool _full_forward;

      /// copy mode in destructor
      bool _full_backward;

      using DT_ = typename CT_::DataType;
      using IT_ = typename CT_::IndexType;
      using MainType_ = typename CT_::template ContainerType<Mem::Main, DT_, IT_>;

      public:
      /// container representation in main memory
      MainType_ main;

      /**
       * Creates a wrapper of the given container object in main memory.
       *
       * \param[in] parent The container object, of which the wrapper shall provide access in main memory.
       * \param[in] full_forward Wether the initial copy to main memory operation shall use the full copy scheme or only copy the index arrays and allocate a new data array in main memory.
       * \param[in] full_backward Wether the final copy back to the device operation shall use the full copy scheme or only copy the data arrays of the container back to the original memory.
       */
      explicit ContainerMainWrapper(CT_ & parent, bool full_forward = false, bool full_backward = false) :
        _parent(parent),
        _full_forward(full_forward),
        _full_backward(full_backward)
      {
        /// \todo copy only what is needed with respect to full_forward flag
        main.convert(_parent);
      }

      ~ContainerMainWrapper()
      {
        _parent.copy(main, _full_backward);
      }

      /// retrieve container representation in main memory
      MainType_& operator*()
      {
        return main;
      }
    }; // class ContainerMainWrapper
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_CONTAINER_MAIN_WRAPPER_HPP
