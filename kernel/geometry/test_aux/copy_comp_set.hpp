// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_TEST_AUX_COPY_COMP_SET_HPP
#define KERNEL_GEOMETRY_TEST_AUX_COPY_COMP_SET_HPP 1

// includes, FEAT
#include <kernel/geometry/index_set.hpp>
#include <kernel/geometry/target_set.hpp>
#include <kernel/geometry/vertex_set.hpp>

// includes, CSL
#include <cmath>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace TestAux
    {
      template<int num_idx_>
      void copy_idx(IndexSet<num_idx_>& idx, const Index jdx[])
      {
        Index n = idx.get_num_entities();
        size_t k = 0;
        for(Index i(0); i < n; ++i)
        {
          for(int j(0); j < num_idx_; ++j)
          {
            idx[i][j] = jdx[k++];
          }
        }
      }

      template<int num_idx_>
      bool comp_idx(const IndexSet<num_idx_>& idx, const Index jdx[])
      {
        Index n = idx.get_num_entities();
        size_t k = 0;
        for(Index i(0); i < n; ++i)
        {
          for(int j(0); j < num_idx_; ++j)
          {
            if(idx[i][j] != jdx[k++])
            {
              return false;
            }
          }
        }
        return true;
      }

      template<typename VertexSet_>
      void copy_vtx(VertexSet_& vtx, const Real wtx[])
      {
        Index n = vtx.get_num_vertices();
        int num_coords = vtx.get_num_coords();
        size_t k = 0;
        for(Index i(0); i < n; ++i)
        {
          for(int j(0); j < num_coords; ++j)
          {
            vtx[i][j] = wtx[k++];
          }
        }
      }

      template<typename VertexSet_>
      bool comp_vtx(const VertexSet_& vtx, const Real wtx[], Real tol = 1e-8)
      {
        Index n = vtx.get_num_vertices();
        int num_coords = vtx.get_num_coords();
        size_t k = 0;
        for(Index i(0); i < n; ++i)
        {
          for(int j(0); j < num_coords; ++j)
          {
            if(std::fabs(vtx[i][j] - wtx[k++]) > tol)
            {
              return false;
            }
          }
        }
        return true;
      }

      inline void copy_trg(TargetSet& trg, const Index jdx[])
      {
        Index n = trg.get_num_entities();
        for(Index i(0); i < n; ++i)
        {
          trg[i] = jdx[i];
        }
      }

      inline bool comp_trg(const TargetSet& trg, const Index jdx[])
      {
        Index n = trg.get_num_entities();
        for(Index i(0); i < n; ++i)
        {
          if(trg[i] != jdx[i])
          {
            return false;
          }
        }
        return true;
      }

      template<typename AttribSet_>
      void copy_attr(AttribSet_& attr, const Real wtx[])
      {
        Index n = attr.get_num_values();
        int num_coords = attr.get_dimension();
        size_t k = 0;
        for(Index i(0); i < n; ++i)
        {
          for(int j(0); j < num_coords; ++j, ++k)
          {
            attr(i,j) = wtx[k];
          }
        }
      }

      template<typename AttribSet_>
      bool comp_attr(const AttribSet_& attr, const Real wtx[], Real tol = 1e-8)
      {
        Index n = attr.get_num_values();
        int num_coords = attr.get_dimension();
        size_t k = 0;
        for(Index i(0); i < n; ++i)
        {
          for(int j(0); j < num_coords; ++j, ++k)
          {
            if(std::fabs(attr(i,j) - wtx[k]) > tol)
            {
              return false;
            }
          }
        }
        return true;
      }
    } // namespace TestAux
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_TEST_AUX_COPY_COMP_SET_HPP
