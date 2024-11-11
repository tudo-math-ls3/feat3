// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_MIRROR_GENERIC_HPP
#define KERNEL_LAFEM_ARCH_MIRROR_GENERIC_HPP 1

#ifndef KERNEL_LAFEM_ARCH_MIRROR_HPP
#error "Do not include this implementation-only header file directly!"
#endif

/// \cond internal
namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      template<typename DT_, typename IT_>
      void Mirror::gather_dv_generic(const Index boff, const Index nidx, const IT_* idx, DT_* buf, const DT_* vec)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0; i < nidx; ++i)
        {
          buf[boff+i] = vec[idx[i]];
        }
      }

      template<typename DT_, typename IT_>
      void Mirror::scatter_dv_generic(const Index boff, const Index nidx, const IT_* idx, const DT_* buf, DT_* vec, const DT_ alpha)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0; i < nidx; ++i)
        {
          vec[idx[i]] += alpha*buf[boff+i];
        }
      }

      template<typename DT_, typename IT_>
      void Mirror::gather_dvb_generic(const Index bs, const Index boff, const Index nidx, const IT_* idx, DT_* buf, const DT_* vec)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0; i < nidx; ++i)
        {
          for(Index k(0); k < bs; ++k)
          {
            buf[boff+i*bs+k] = vec[idx[i]*bs+k];
          }
        }
      }

      template<typename DT_, typename IT_>
      void Mirror::scatter_dvb_generic(const Index bs, const Index boff, const Index nidx, const IT_* idx, const DT_* buf, DT_* vec, const DT_ alpha)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0; i < nidx; ++i)
        {
          for(Index k(0); k < bs; ++k)
          {
            vec[idx[i]*bs+k] += alpha*buf[boff+i*bs+k];
          }
        }
      }

      template<typename DT_, typename IT_>
      void Mirror::gather_sv_generic(const Index boff, const Index nidx, const IT_* idx, DT_* buf, const Index nvec, const DT_* vval, const IT_* vidx)
      {
        // loop over all mirror indices
        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0; i < nidx; ++i)
        {
          buf[boff+i] = DT_(0);

          // loop over all vector indices and try to find the match
          for(Index j = 0; j < nvec; ++j)
          {
            if(idx[i] == vidx[j])
            {
              buf[boff+i] = vval[j];
              break;
            }
          }
        }
      }

      template<typename DT_, typename IT_>
      void Mirror::scatter_sv_generic(const Index boff, const Index nidx, const IT_* idx, const DT_* buf, const Index nvec, DT_* vval, const IT_* vidx, const DT_ alpha)
      {
        // loop over all mirror indices
        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0; i < nidx; ++i)
        {
          // loop over all vector indices and try to find the match
          for(Index j(0); j < nvec; ++j)
          {
            if(idx[i] == vidx[j])
            {
              vval[j] += alpha*buf[boff+i];
              break;
            }
          }
        }
      }

      template<typename DT_, typename IT_>
      void Mirror::gather_svb_generic(const Index bs, const Index boff, const Index nidx, const IT_* idx, DT_* buf, const Index nvec, const DT_* vval, const IT_* vidx)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0; i < nidx; ++i)
        {
          for(Index k(0); k < bs; ++k)
            buf[boff+i*bs+k] = DT_(0);

          for(Index j(0); j < nvec; ++j)
          {
            if(idx[i] == vidx[j])
            {
              for(Index k(0); k < bs; ++k)
                buf[boff+i*bs+k] = vval[j*bs+k];
              break;
            }
          }
        }
      }

      template<typename DT_, typename IT_>
      void Mirror::scatter_svb_generic(const Index bs, const Index boff, const Index nidx, const IT_* idx, const DT_* buf, const Index nvec, DT_* vval, const IT_* vidx, const DT_ alpha)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for(Index i = 0; i < nidx; ++i)
        {
          for(Index j(0); j < nvec; ++j)
          {
            if(idx[i] == vidx[j])
            {
              for(Index k(0); k < bs; ++k)
                vval[j*bs+k] += alpha*buf[boff+i*bs+k];
              break;
            }
          }
        }
      }
    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT
/// \endcond

#endif // KERNEL_LAFEM_ARCH_MIRROR_GENERIC_HPP
