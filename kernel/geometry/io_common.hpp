// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/util/string.hpp>

namespace FEAT::Geometry
{
  /// mesh type enumeration
  enum class ParsedMeshType
  {
    /// unknown mesh type
    unknown = 0,
    /// conformal mesh type
    conformal
  };

  inline std::ostream& operator<<(std::ostream& s, ParsedMeshType meshtype)
  {
    switch(meshtype)
    {
      case ParsedMeshType::unknown: s << "unknown"; break;
      case ParsedMeshType::conformal: s << "conformal"; break;
    }
    return s;
  }

  inline std::istream& operator>>(std::istream& is, ParsedMeshType& meshtype)
  {
    String s;
    if((is >> s).fail())
    {
      return is;
    }

    if(s.compare_no_case("unknown") == 0)
    {
      meshtype = ParsedMeshType::unknown;
    }
    else if(s.compare_no_case("conformal") == 0)
    {
      meshtype = ParsedMeshType::conformal;
    }
    else
    {
      is.setstate(std::ios_base::failbit);
    }

    return is;
  }

  /// shape type enumeration
  enum class ParsedShapeType
  {
    /// unknown shape type
    unknown = 0,
    /// simplex shape type
    simplex,
    /// hypercube shape type
    hypercube
  };

  inline std::ostream& operator<<(std::ostream& s, ParsedShapeType shapetype)
  {
    switch(shapetype)
    {
      case ParsedShapeType::unknown: s << "unknown"; break;
      case ParsedShapeType::simplex: s << "simplex"; break;
      case ParsedShapeType::hypercube: s << "hypercube"; break;
    }
    return s;
  }

  inline std::istream& operator>>(std::istream& is, ParsedShapeType& shapetype)
  {
    String s;
    if((is >> s).fail())
    {
      return is;
    }

    if(s.compare_no_case("unknown") == 0)
    {
      shapetype = ParsedShapeType::unknown;
    }
    else if(s.compare_no_case("simplex") == 0)
    {
      shapetype = ParsedShapeType::simplex;
    }
    else if(s.compare_no_case("hypercube") == 0)
    {
      shapetype = ParsedShapeType::hypercube;
    }
    else
    {
      is.setstate(std::ios_base::failbit);
    }

    return is;
  }
} // namespace FEAT::Geometry
