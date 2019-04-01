// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef FEAT_CONTROL_TIME_BASE_HPP
#define FEAT_CONTROL_TIME_BASE_HPP 1
#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>
#include <iostream>

namespace FEAT
{
  namespace Control
  {
    /**
     * Time discretisations
     */
    namespace Time
    {
      /**
       * \brief For determining if a term is handled explicitly or implicitly
       */
      enum class TermHandling
      {
        off = 0,
        expl,
        impl
      };

      /// \cond internal
      inline std::ostream& operator<<(std::ostream& os, TermHandling term_handling)
      {
        switch(term_handling)
        {
          case TermHandling::off:
            return os << "off";
          case TermHandling::expl:
            return os << "expl";
          case TermHandling::impl:
            return os << "impl";
          default:
            return os << "-unknown-";
        }
      }

      inline std::istream& operator>>(std::istream& is, TermHandling& term_handling)
      {
        String term_handling_name;


        if( (is >> term_handling_name).fail() )
        {
          return is;
        }

        if(term_handling_name == "off")
        {
          term_handling = TermHandling::off;
        }
        else if(term_handling_name == "expl")
        {
          term_handling = TermHandling::expl;
        }
        else if(term_handling_name == "impl")
        {
          term_handling = TermHandling::impl;
        }
        else
        {
          is.setstate(std::ios_base::failbit);
        }

        return is;

      }

      /// \endcond
    } // namespace Time
  } // namespace Control
} // namespace FEAT
#endif // FEAT_CONTROL_TIME_BASE_HPP
