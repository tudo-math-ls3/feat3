// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <functional>

#include <kernel/adjacency/adjactor.hpp>
#include <kernel/adjacency/base.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/base_header.hpp>
#include <kernel/global/remote_lambda.hpp>
#include <kernel/trafo/inverse_mapping.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/tiny_algebra.hpp>

namespace FEAT::Global
{
  /**
   * \brief Data structure for InverseMapping evaluations
   *
   * This class is used for storing the resulting data that arises
   * from an inverse mapping evaluation.
   * Similar to Trafo::InverseMappingData but extended by rank information.
   *
   * \tparam DataType_
   * The datatype that is used for coordinates.
   *
   * \tparam shape_dim_
   * The shape dimension of the underlying mesh.
   *
   * \tparam world_dim_
   * The world dimension of the underlying mesh.
   *
   * \author Markus Muegge
   */
  template<typename DataType_, int shape_dim_, int world_dim_ = shape_dim_>
  class GlobalInverseMappingData
  {
  public:
    /// the domain point type
    using DomainPointType = Tiny::Vector<DataType_, shape_dim_>;
    /// the image point type
    using ImagePointType = Tiny::Vector<DataType_, world_dim_>;

    /// the image point that was unmapped
    ImagePointType img_point;

    /// the ranks of the cells that intersect with the image point
    std::vector<int> ranks;

    /// the indices of the cells that intersect with the image point
    std::vector<Index> cells;

    /// the domain points on each cell that map onto the image point
    std::vector<DomainPointType> dom_points;

    /// default constructor
    GlobalInverseMappingData() = default;

    /// use default copy constructor
    GlobalInverseMappingData(const GlobalInverseMappingData& other) = default;

    /// default copy assignment operator
    GlobalInverseMappingData& operator=(const GlobalInverseMappingData& other) = default;

    /// move constructor
    GlobalInverseMappingData(GlobalInverseMappingData&& other) noexcept :
      img_point(other.img_point),
      ranks(std::move(other.ranks)),
      cells(std::move(other.cells)),
      dom_points(std::move(other.dom_points))
    {
    }

    /// move-assign operator
    GlobalInverseMappingData& operator=(GlobalInverseMappingData&& other)
    {
      if(this != &other)
      {
        img_point = other.img_point;
        ranks = std::move(other.ranks);
        cells = std::move(other.cells);
        dom_points = std::move(other.dom_points);
      }
      return *this;
    }

    /**
     * \brief Checks whether no cells were found.
     */
    bool empty() const
    {
      return cells.empty();
    }

    /// \returns The numbers of cells that were found.
    std::size_t size() const
    {
      return cells.size();
    }
  }; // class GlobalInverseMappingData

  /**
   * \brief Global inverse mapping
   *
   * This class template can be used to unmap points on a distributed domain.
   *
   * The unmapping works as follows:
   * - The ranks exchange the axis-aligned bounding boxes of their patches
   * - The ranks exchange points, such that each ranks has all points that
   *   are contained in the bounding box of its patch
   * - The ranks unmap their points and send resulting cells and domain points
   *   back to the ranks that originally asked to unmap that point
   */
  template<typename InverseMappingType_>
  class InverseMapping
  {
  public:
    /// Local inverse mapping type
    using InverseMappingType = InverseMappingType_;
    /// Trafo type
    using TrafoType = typename InverseMappingType::TrafoType;
    /// Datatype
    using DataType = typename InverseMappingType::DataType;
    /// Shape type
    using ShapeType = typename InverseMappingType::ShapeType;
    /// Shape dimension
    static constexpr int shape_dim = ShapeType::dimension;
    /// World dimension
    static constexpr int world_dim = TrafoType::world_dim;
    /// The type returned by the unmap method
    using InvMapDataType = GlobalInverseMappingData<DataType, shape_dim, world_dim>;
    /// Point type for points on the mesh
    using ImagePointType = typename InvMapDataType::ImagePointType;
    /// Point type for points on the reference element
    using DomainPointType = typename InvMapDataType::DomainPointType;
    /// Axis-aligned bounding box type
    using BoundingBoxType = Tiny::Matrix<DataType, 2, world_dim>;

  private:
    /// Communicator that domain is distributed across
    Dist::Comm _comm;
    /// Reference to local inverse mapping
    std::reference_wrapper<const InverseMappingType> _mapping;
    /// Bounding box of local patch
    BoundingBoxType _local_bb;

    /// Pair for communicating cells and associated domain points between ranks
    struct UnmappingResult
    {
      Index cell;
      DomainPointType dom_point;

      /// Make type MPI-compatible
      static Dist::Typemap typemap()
      {
        return Dist::Typemap()
          .add_entry<Index>(offsetof(UnmappingResult, cell))
          .add_entry<DomainPointType>(offsetof(UnmappingResult, dom_point));
      }
    };

  public:

    /// Constructor
    InverseMapping(const Dist::Comm& comm, const InverseMappingType& m) : _comm(comm.comm_dup()), _mapping(m), _local_bb(get_local_bb())
    {
    }

    /**
     * \brief Unmap a range of points
     *
     * \tparam PointRange_ Point range type
     *
     * \param[in] points Range of points to unmap
     *
     * \returns An InvMapDataType instance for each point in the range
     */
    template<typename PointRange_>
    std::vector<InvMapDataType> unmap_points(const PointRange_& points)
    {
      const int total_ranks = _comm.size();

      // Exchange bounding boxes with other ranks
      std::vector<BoundingBoxType> bbs(static_cast<std::size_t>(total_ranks));
      bbs.at(std::size_t(_comm.rank())) = _local_bb;
      _comm.allgather(bbs.data(), 1, bbs.data(), 1);

      // Build structure for exchanging points
      std::vector<int> ranks;
      std::vector<Index> indices;

      Index point_idx(0);
      for(auto it = points.begin(); it != points.end(); ++it)
      {
        const ImagePointType& p = *it;
        for(int i(0); i < total_ranks; i++)
        {
          const BoundingBoxType& bb = bbs[std::size_t(i)];

          if(is_in_bb(bb, p))
          {
            ranks.push_back(i);
            indices.push_back(point_idx);
          }
        }
        point_idx++;
      }

      VectorRemoteLambda<std::vector<UnmappingResult>(const ImagePointType&)> rl(_comm);
      std::vector<std::vector<UnmappingResult>> result_pairs = rl.call(ranks, indices, points, [this](const ImagePointType& p)
      {
        // NOTE(mmuegge): ignore_failure is always set to true, because the RemoteLambda call is a collective operations
        // and one process throwing an exception will cause a deadlock.
        auto local_result = _mapping.get().unmap_point(p, true);

        // Convert local result to result pairs for mpi sending
        std::vector<UnmappingResult> result;
        result.reserve(local_result.size());
        for(Index i(0); i < local_result.size(); i++)
        {
          result.push_back({local_result.cells[i], local_result.dom_points[i]});
        }

        return result;
      });


      // Assemble final result from ResultPairs
      std::vector<InvMapDataType> result(std::size_t(std::distance(points.begin(), points.end())));
      for(std::size_t i(0); i < ranks.size(); i++)
      {
        const int rank = ranks[i];
        const Index idx = indices[i];
        auto iter = points.begin();
        std::advance(iter, idx);
        const ImagePointType& image_point = *iter;

        std::vector<UnmappingResult>& pairs = result_pairs[i];

        result[idx].img_point = image_point;
        for(const UnmappingResult& pair : pairs)
        {
          result[idx].cells.push_back(pair.cell);
          result[idx].dom_points.push_back(pair.dom_point);
          result[idx].ranks.push_back(rank);
        }
      }

      return result;
    }

  private:
    /// Check if a point is in a bounding box
    static bool is_in_bb(const BoundingBoxType& bb, const ImagePointType& p)
    {
      constexpr DataType bb_eps(1E-2);
      for(int i(0); i < world_dim; i++)
      {
        if(p[i] > (bb[1][i] + bb_eps) || p[i] < (bb[0][i] - bb_eps))
        {
          return false;
        }
      }
      return true;
    }

    /// Get bounding box of local part of domain
    BoundingBoxType get_local_bb()
    {
      const auto& vertex_set = _mapping.get().get_trafo().get_mesh().get_vertex_set();

      if(vertex_set.get_num_vertices() == 0)
      {
        return BoundingBoxType{0};
      }

      // Determine local bounding box
      BoundingBoxType local_bb;

      local_bb[0] = vertex_set[0];
      local_bb[1] = vertex_set[0];
      for(Index vtx_idx(0); vtx_idx < vertex_set.get_num_vertices(); vtx_idx++)
      {
        const auto& v = vertex_set[vtx_idx];

        for(int i(0); i < world_dim; i++)
        {
          Math::mini(local_bb[0][i], v[i]);
          Math::maxi(local_bb[1][i], v[i]);
        }
      }

      return local_bb;
    }
  }; // class InverseMapping
} // namespace FEAT::Global
