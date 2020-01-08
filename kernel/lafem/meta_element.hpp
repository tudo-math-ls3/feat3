// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_META_ELEMENT_HPP
#define KERNEL_LAFEM_META_ELEMENT_HPP 1

#include <kernel/base_header.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Tuple container element helper class template
     *
     * This class template is a helper which is used for the implementation
     * of the 'at' function template of various Tuple meta-containers.
     *
     * \author Peter Zajac
     */
    template<
      int i_,
      typename First_,
      typename... Rest_>
    struct TupleElement
    {
      static_assert(i_ > 0, "invalid index");

      typedef typename TupleElement<i_-1, Rest_...>::Type Type;

      template<typename Meta_>
      static Type& get(Meta_& meta)
      {
        return TupleElement<i_-1, Rest_...>::get(meta.rest());
      }

      template<typename Meta_>
      static const Type& get(const Meta_& meta)
      {
        return TupleElement<i_-1, Rest_...>::get(meta.rest());
      }
    };

    /// \cond internal
    template<typename First_, typename... Rest_>
    struct TupleElement<0, First_, Rest_...>
    {
      typedef First_ Type;

      template<typename Meta_>
      static Type& get(Meta_& meta)
      {
        return meta.first();
      }

      template<typename Meta_>
      static const Type& get(const Meta_& meta)
      {
        return meta.first();
      }
    };
    /// \endcond

    /**
     * \brief Power container element helper class template
     *
     * This class template is a helper which is used for the implementation
     * of the 'at' function template of various Power meta-containers.
     *
     * \author Peter Zajac
     */
    template<int i_, typename SubType_>
    struct PowerElement
    {
      static_assert(i_ > 0, "invalid index");

      template<typename Meta_>
      static SubType_& get(Meta_& meta)
      {
        return PowerElement<i_-1, SubType_>::get(meta.rest());
      }

      template<typename Meta_>
      static const SubType_& get(const Meta_& meta)
      {
        return PowerElement<i_-1, SubType_>::get(meta.rest());
      }
    };

    /// \cond internal
    template<typename SubType_>
    struct PowerElement<0, SubType_>
    {
      template<typename Meta_>
      static SubType_& get(Meta_& meta)
      {
        return meta.first();
      }

      template<typename Meta_>
      static const SubType_& get(const Meta_& meta)
      {
        return meta.first();
      }
    };
    /// \endcond
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_META_ELEMENT_HPP
