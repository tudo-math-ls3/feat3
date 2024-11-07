// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_ARCH_MIRROR_HPP
#define KERNEL_LAFEM_ARCH_MIRROR_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    namespace Arch
    {
      struct Mirror
      {
        template<typename DT_, typename IT_>
        static void gather_dv(const Index boff, const Index nidx, const IT_* idx, DT_* buf, const DT_* vec)
        {
          gather_dv_generic(boff, nidx, idx, buf, vec);
        }

        static void gather_dv(const Index boff, const Index nidx, const std::uint32_t * idx, float * buf, const float * vec)
        {
          BACKEND_SKELETON_VOID(gather_dv_cuda, gather_dv_generic, gather_dv_generic, boff, nidx, idx, buf, vec)
        }

        static void gather_dv(const Index boff, const Index nidx, const std::uint64_t * idx, float * buf, const float * vec)
        {
          BACKEND_SKELETON_VOID(gather_dv_cuda, gather_dv_generic, gather_dv_generic, boff, nidx, idx, buf, vec)
        }

        static void gather_dv(const Index boff, const Index nidx, const std::uint32_t * idx, double * buf, const double * vec)
        {
          BACKEND_SKELETON_VOID(gather_dv_cuda, gather_dv_generic, gather_dv_generic, boff, nidx, idx, buf, vec)
        }

        static void gather_dv(const Index boff, const Index nidx, const std::uint64_t * idx, double * buf, const double * vec)
        {
          BACKEND_SKELETON_VOID(gather_dv_cuda, gather_dv_generic, gather_dv_generic, boff, nidx, idx, buf, vec)
        }

        template<typename DT_, typename IT_>
        static void scatter_dv(const Index boff, const Index nidx, const IT_* idx, const DT_* buf, DT_* vec, const DT_ alpha)
        {
          scatter_dv_generic(boff, nidx, idx, buf, vec, alpha);
        }

        static void scatter_dv(const Index boff, const Index nidx, const std::uint32_t * idx, const float * buf, float * vec, const float alpha)
        {
          BACKEND_SKELETON_VOID(scatter_dv_cuda, scatter_dv_generic, scatter_dv_generic, boff, nidx, idx, buf, vec, alpha)
        }

        static void scatter_dv(const Index boff, const Index nidx, const std::uint64_t * idx, const float * buf, float * vec, const float alpha)
        {
          BACKEND_SKELETON_VOID(scatter_dv_cuda, scatter_dv_generic, scatter_dv_generic, boff, nidx, idx, buf, vec, alpha)
        }

        static void scatter_dv(const Index boff, const Index nidx, const std::uint32_t * idx, const double * buf, double * vec, const double alpha)
        {
          BACKEND_SKELETON_VOID(scatter_dv_cuda, scatter_dv_generic, scatter_dv_generic, boff, nidx, idx, buf, vec, alpha)
        }

        static void scatter_dv(const Index boff, const Index nidx, const std::uint64_t * idx, const double * buf, double * vec, const double alpha)
        {
          BACKEND_SKELETON_VOID(scatter_dv_cuda, scatter_dv_generic, scatter_dv_generic, boff, nidx, idx, buf, vec, alpha)
        }

        template<typename DT_, typename IT_>
        static void gather_dvb(const Index bs, const Index boff, const Index nidx, const IT_* idx, DT_* buf, const DT_* vec)
        {
          gather_dvb_generic(bs, boff, nidx, idx, buf, vec);
        }

        static void gather_dvb(const Index bs, const Index boff, const Index nidx, const std::uint32_t * idx, float * buf, const float * vec)
        {
          BACKEND_SKELETON_VOID(gather_dvb_cuda, gather_dvb_generic, gather_dvb_generic, bs, boff, nidx, idx, buf, vec)
        }

        static void gather_dvb(const Index bs, const Index boff, const Index nidx, const std::uint64_t * idx, float * buf, const float * vec)
        {
          BACKEND_SKELETON_VOID(gather_dvb_cuda, gather_dvb_generic, gather_dvb_generic, bs, boff, nidx, idx, buf, vec)
        }

        static void gather_dvb(const Index bs, const Index boff, const Index nidx, const std::uint32_t * idx, double * buf, const double * vec)
        {
          BACKEND_SKELETON_VOID(gather_dvb_cuda, gather_dvb_generic, gather_dvb_generic, bs, boff, nidx, idx, buf, vec)
        }

        static void gather_dvb(const Index bs, const Index boff, const Index nidx, const std::uint64_t * idx, double * buf, const double * vec)
        {
          BACKEND_SKELETON_VOID(gather_dvb_cuda, gather_dvb_generic, gather_dvb_generic, bs, boff, nidx, idx, buf, vec)
        }

        template<typename DT_, typename IT_>
        static void scatter_dvb(const Index bs, const Index boff, const Index nidx, const IT_* idx, const DT_* buf, DT_* vec, const DT_ alpha)
        {
          scatter_dvb_generic(bs, boff, nidx, idx, buf, vec, alpha);
        }

        static void scatter_dvb(const Index bs, const Index boff, const Index nidx, const std::uint32_t * idx, const float * buf, float * vec, const float alpha)
        {
          BACKEND_SKELETON_VOID(scatter_dvb_cuda, scatter_dvb_generic, scatter_dvb_generic, bs, boff, nidx, idx, buf, vec, alpha)
        }

        static void scatter_dvb(const Index bs, const Index boff, const Index nidx, const std::uint64_t * idx, const float * buf, float * vec, const float alpha)
        {
          BACKEND_SKELETON_VOID(scatter_dvb_cuda, scatter_dvb_generic, scatter_dvb_generic, bs, boff, nidx, idx, buf, vec, alpha)
        }

        static void scatter_dvb(const Index bs, const Index boff, const Index nidx, const std::uint32_t * idx, const double * buf, double * vec, const double alpha)
        {
          BACKEND_SKELETON_VOID(scatter_dvb_cuda, scatter_dvb_generic, scatter_dvb_generic, bs, boff, nidx, idx, buf, vec, alpha)
        }

        static void scatter_dvb(const Index bs, const Index boff, const Index nidx, const std::uint64_t * idx, const double * buf, double * vec, const double alpha)
        {
          BACKEND_SKELETON_VOID(scatter_dvb_cuda, scatter_dvb_generic, scatter_dvb_generic, bs, boff, nidx, idx, buf, vec, alpha)
        }

        template<typename DT_, typename IT_>
        static void gather_sv(const Index boff, const Index nidx, const IT_* idx, DT_* buf, const Index nvec, const DT_* vval, const IT_* vidx)
        {
          gather_sv_generic(boff, nidx, idx, buf, nvec, vval, vidx);
        }

        template<typename DT_, typename IT_>
        static void scatter_sv(const Index boff, const Index nidx, const IT_* idx, const DT_* buf, const Index nvec, DT_* vval, const IT_* vidx, const DT_ alpha)
        {
          scatter_sv_generic(boff, nidx, idx, buf, nvec, vval, vidx, alpha);
        }

        template<typename DT_, typename IT_>
        static void gather_svb(const Index bs, const Index boff, const Index nidx, const IT_* idx, DT_* buf, const Index nvec, const DT_* vval, const IT_* vidx)
        {
          gather_svb_generic(bs, boff, nidx, idx, buf, nvec, vval, vidx);
        }

        template<typename DT_, typename IT_>
        static void scatter_svb(const Index bs, const Index boff, const Index nidx, const IT_* idx, const DT_* buf, const Index nvec, DT_* vval, const IT_* vidx, const DT_ alpha)
        {
          scatter_svb_generic(bs, boff, nidx, idx, buf, nvec, vval, vidx, alpha);
        }

        // generic declarations

        template<typename DT_, typename IT_>
        static void gather_dv_generic(const Index boff, const Index nidx, const IT_* idx, DT_* buf, const DT_* vec);

        template<typename DT_, typename IT_>
        static void scatter_dv_generic(const Index boff, const Index nidx, const IT_* idx, const DT_* buf, DT_* vec, const DT_ alpha);

        template<typename DT_, typename IT_>
        static void gather_dvb_generic(const Index bs, const Index boff, const Index nidx, const IT_* idx, DT_* buf, const DT_* vec);

        template<typename DT_, typename IT_>
        static void scatter_dvb_generic(const Index bs, const Index boff, const Index nidx, const IT_* idx, const DT_* buf, DT_* vec, const DT_ alpha);

        template<typename DT_, typename IT_>
        static void gather_sv_generic(const Index boff, const Index nidx, const IT_* idx, DT_* buf, const Index nvec, const DT_* vval, const IT_* vidx);

        template<typename DT_, typename IT_>
        static void scatter_sv_generic(const Index boff, const Index nidx, const IT_* idx, const DT_* buf, const Index nvec, DT_* vval, const IT_* vidx, const DT_ alpha);

        template<typename DT_, typename IT_>
        static void gather_svb_generic(const Index bs, const Index boff, const Index nidx, const IT_* idx, DT_* buf, const Index nvec, const DT_* vval, const IT_* vidx);

        template<typename DT_, typename IT_>
        static void scatter_svb_generic(const Index bs, const Index boff, const Index nidx, const IT_* idx, const DT_* buf, const Index nvec, DT_* vval, const IT_* vidx, const DT_ alpha);

        template<typename DT_, typename IT_>
        static void gather_dv_cuda(const Index boff, const Index nidx, const IT_* idx, DT_* buf, const DT_* vec);

        template<typename DT_, typename IT_>
        static void scatter_dv_cuda(const Index boff, const Index nidx, const IT_* idx, const DT_* buf, DT_* vec, const DT_ alpha);

        template<typename DT_, typename IT_>
        static void gather_dvb_cuda(const Index bs, const Index boff, const Index nidx, const IT_* idx, DT_* buf, const DT_* vec);

        template<typename DT_, typename IT_>
        static void scatter_dvb_cuda(const Index bs, const Index boff, const Index nidx, const IT_* idx, const DT_* buf, DT_* vec, const DT_ alpha);

        /// \todo implement cuda operations for SparseVector(Blocked) if required

      }; // struct Mirror

    } // namespace Arch
  } // namespace LAFEM
} // namespace FEAT

#ifndef __CUDACC__
#include <kernel/lafem/arch/mirror_generic.hpp>
#endif
#endif // KERNEL_LAFEM_ARCH_MIRROR_HPP
