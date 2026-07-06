// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/util/memory_arbiter.hpp>
#include <kernel/util/omp_util.hpp>
#include <kernel/backend.hpp>



namespace FEAT::LAFEM::Arch
{
  struct BatchedScratchSpace
  {
    // free scratch space of batched kernel
    static void free();
  };


  // The local dot-product results are passed as a Memory::Arbiter (rather than a raw host pointer), so
  // the batched kernels acquire the view they need themselves and the arbiter transparently mediates
  // the host<->device transfers and the on-device accumulation across meta-vector components. The
  // arbiter must already hold the (zeroed or partially accumulated) values that are added upon.
  struct BatchedDotProduct
  {
    template<typename DT_, typename DTT_>
    static void exec_generic(Memory::Arbiter& local_dots, const Memory::Arbiter& target, Index dot_size, Memory::Arbiter const* const* others, Index num_dots)
    {
      // the generic kernel operates on main memory, so it acquires main-memory views itself
      Memory::TypedView<DT_> dots_view(local_dots.view(Memory::Location::main, Memory::Access::read_write));
      DT_* dots_ptr = dots_view.get_w();
      const Memory::TypedView<DTT_> target_view(target.view(Memory::Location::main, Memory::Access::read));
      const DTT_* target_ptr = target_view.get_r();

      for(Index l = 0; l < num_dots; ++l)
      {
        const Memory::TypedView<DTT_> other_view(others[l]->view(Memory::Location::main, Memory::Access::read));
        const DTT_* other_ptr = other_view.get_r();
        DTT_ loc_v = DTT_(0);
        FEAT_PRAGMA_OMP(parallel for reduction(+:loc_v))
        for(Index k = 0; k < dot_size; ++k)
        {
          loc_v += DT_(target_ptr[k] * other_ptr[k]);
        }
        dots_ptr[l] += loc_v;
      }
    }

    template<typename DT_, typename DTT_>
    static void exec_cuda(Memory::Arbiter& local_dots, const Memory::Arbiter& target, Index dot_size, Memory::Arbiter const * const * others, Index num_dots);

    template<typename DT_, typename DTT_>
    static void exec(Memory::Arbiter& local_dots, const Memory::Arbiter& target, Index dot_size, Memory::Arbiter const * const * others, Index num_dots)
    {
      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::cuda:
      #ifdef FEAT_HAVE_CUDA
        exec_cuda<DT_, DTT_>(local_dots, target, dot_size, others, num_dots);
      #else
        XABORTM("LAFEM::ARCH::BatchedDotProduct::exec: CUDA backend not available");
      #endif
        break;
      default:
        exec_generic<DT_, DTT_>(local_dots, target, dot_size, others, num_dots);
        break;
      }
    }

  };


  struct BatchedTripleDotProduct
  {

    template<typename DT_, typename DTT_>
    static inline void exec_generic(Memory::Arbiter& local_dots, const Memory::Arbiter& target, Index dot_size, Memory::Arbiter const * const * others, const Memory::Arbiter& freqs, Index num_dots)
    {
      // the generic kernel operates on main memory, so it acquires main-memory views itself
      Memory::TypedView<DT_> dots_view(local_dots.view(Memory::Location::main, Memory::Access::read_write));
      DT_* dots_ptr = dots_view.get_w();
      Memory::TypedView<DTT_> target_view(target.view(Memory::Location::main, Memory::Access::read));
      const DTT_* target_ptr = target_view.get_r();
      Memory::TypedView<DTT_> freqs_view(freqs.view(Memory::Location::main, Memory::Access::read));
      const DTT_* freqs_ptr = freqs_view.get_r();

      for(Index l = 0; l < num_dots; ++l)
      {
        const Memory::TypedView<DTT_> other_view(others[l]->view(Memory::Location::main, Memory::Access::read));
        const DTT_* other_ptr = other_view.get_r();
        DTT_ loc_v = DTT_(0);
        FEAT_PRAGMA_OMP(parallel for reduction(+:loc_v))
        for(Index k = 0; k < dot_size; ++k)
        {
          loc_v += DT_(target_ptr[k] * freqs_ptr[k] * other_ptr[k]);
        }
        dots_ptr[l] += loc_v;
      }
    }

    template<typename DT_, typename DTT_>
    static void exec_cuda(Memory::Arbiter& local_dots, const Memory::Arbiter& target, Index dot_size, Memory::Arbiter const * const * others, const Memory::Arbiter& freqs, Index num_dots);

    template<typename DT_, typename DTT_>
    static void exec(Memory::Arbiter& local_dots, const Memory::Arbiter& target, Index dot_size, Memory::Arbiter const * const * others, const Memory::Arbiter& freqs, Index num_dots)
    {
      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::cuda:
      #ifdef FEAT_HAVE_CUDA
        exec_cuda<DT_, DTT_>(local_dots, target, dot_size, others, freqs, num_dots);
      #else
        XABORTM("LAFEM::ARCH::BatchedTripleDotProduct::exec: CUDA backend not available");
      #endif
        break;
      default:
        exec_generic<DT_, DTT_>(local_dots, target, dot_size, others, freqs, num_dots);
        break;
      }
    }
  };



  struct CumulativeAxpy
  {
    template<typename DT_>
    static void exec_generic(Memory::Arbiter& target, Index vec_size, Memory::Arbiter const * const * others, DT_ const * const local_scales, Index num_axpy, DT_ alpha)
    {
      // the generic kernel operates on main memory, so it acquires its main-memory views itself
      Memory::TypedView<DT_> target_view(target.view(Memory::Location::main, Memory::Access::read_write));
      DT_* target_ptr = target_view.get_w();

      for(Index l = 0; l < num_axpy; ++l)
      {
        const Memory::TypedView<DT_> other_view(others[l]->view(Memory::Location::main, Memory::Access::read));
        const DT_* other_ptr = other_view.get_r();
        FEAT_PRAGMA_OMP(parallel for)
        for(Index k = 0; k < vec_size; ++k)
        {
          target_ptr[k] += DT_(alpha * local_scales[l]) * other_ptr[k];
        }
      }
    }

    template<typename DT_>
    static void exec_cuda(Memory::Arbiter& target, Index vec_size, Memory::Arbiter const * const * others, DT_ const * const local_scales, Index num_axpy, DT_ alpha);

    template<typename DT_>
    static void exec(Memory::Arbiter& target, Index vec_size, Memory::Arbiter const * const * others, DT_ const * const local_scales, Index num_axpy, DT_ alpha)
    {
      switch(Backend::get_preferred_backend())
      {
      case PreferredBackend::cuda:
      #ifdef FEAT_HAVE_CUDA
        exec_cuda<DT_>(target, vec_size, others, local_scales, num_axpy, alpha);
      #else
        XABORTM("LAFEM::ARCH::CumulativeAxpy::exec: CUDA backend not available");
      #endif
        break;
      default:
        exec_generic<DT_>(target, vec_size, others, local_scales, num_axpy, alpha);
        break;
      }
    }
  };

}