// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#include <kernel/base_header.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/dist.hpp>

namespace FEAT::Control::Domain
{

  /**
   * \brief Helper class to calculate the global boundingbox of a domain
   */
  template<typename DomainControl_>
  struct DomainGlobalBoundingBoxHelper
  {
    typedef typename DomainControl_::LevelType::MeshType::CoordType DataType;
    static constexpr int dim = DomainControl_::LevelType::MeshType::shape_dim;
    typedef Tiny::Vector<DataType, dim> CoordVecType;

    /**
    * \brief Helper function to calculate global boundingbox on current layer
    *
    * \param[in] domain The domain control for which the global boundingbox should be calculated for
    * \param[in] level The level the boundingbox should be calculated on
    *
    * \returns A pair of Tiny::Vector representing the minimum coordinates and maximum coordinates of the boundingbox
    *
    * \note Has to be called from all processes
    */
    static std::pair<CoordVecType, CoordVecType> calculate_global_boundingbox(const DomainControl_& domain, const Index level = Index(0))
    {
      XASSERTM(domain.size_virtual() >= level, "Chosen level is larger than virtual size");
      // calculate local boundingbox for all relevant domain level
      CoordVecType bb_min(Math::Limits<DataType>::max()), bb_max(Math::Limits<DataType>::lowest());
      if(level < domain.size_physical())
      {
        // calculate local boundingbox
        const auto& mesh = domain.at(level)->get_mesh();
        const auto& vtx_set = mesh.get_vertex_set();

        const Index n = vtx_set.get_num_vertices();
        for(Index i(0); i < n; ++i)
        {
          for(int k = 0; k < dim; ++k)
          {
            bb_min[k] = Math::min(bb_min[k], vtx_set[i][k]);
            bb_max[k] = Math::max(bb_max[k], vtx_set[i][k]);
          }
        }
      }
      // and now do a communication over all processes
      const Dist::Comm& comm = domain.comm();
      comm.allreduce(static_cast<DataType*>(&bb_min[0]), static_cast<DataType*>(&bb_min[0]), std::size_t(dim), Dist::op_min);
      comm.allreduce(static_cast<DataType*>(&bb_max[0]), static_cast<DataType*>(&bb_max[0]), std::size_t(dim), Dist::op_max);

      return std::make_pair(bb_min, bb_max);
    }
  };

  /**
   * \brief Specilization for voxel domain control
   */
  template<typename DomainLevel_>
  struct DomainGlobalBoundingBoxHelper<Control::Domain::DeslaggedDomainControl<DomainLevel_>>
  {
    typedef typename DomainLevel_::MeshType::CoordType DataType;
    static constexpr int dim = DomainLevel_::MeshType::shape_dim;
    typedef Tiny::Vector<DataType, dim> CoordVecType;

    /**
    * \brief Helper function to calculate global boundingbox on current layer
    *
    * \param[in] domain The domain control for which the global boundingbox should be calculated for
    * \param[in] level The level the boundingbox should be calculated on
    *
    * \returns A pair of Tiny::Vector representing the minimum coordinates and maximum coordinates of the boundingbox
    *
    * \note Has to be called from all processes
    */
    static std::pair<CoordVecType, CoordVecType> calculate_global_boundingbox(const Control::Domain::DeslaggedDomainControl<DomainLevel_>& domain, [[maybe_unused]] const Index level = Index(0))
    {
      XASSERTM(domain.size_virtual() > level, "Chosen level is larger than virtual size");
      // calculate local boundingbox for all relevant domain level
      return std::make_pair(domain.get_bounding_box_min(), domain.get_bounding_box_max());
    }
  };

}
