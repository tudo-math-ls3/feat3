// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/runtime.hpp>
#include <kernel/backend.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/random.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/memory_pool.hpp>
#include <kernel/util/stop_watch.hpp>
#include <kernel/lafem/dense_vector.hpp>

#include <iostream>

#ifdef FEAT_HAVE_OMP
#include <omp.h>
#endif

using namespace FEAT;

typedef std::uint64_t u64;

static constexpr u64 pad_len = 32;

// OpenMP-parallel copy; must be declared NOINLINE to prevent the compiler from reordering instructions,
// which could move the operation out of the time measurement block
NOINLINE void omp_copy(u64* x, u64* y, u64 num_entries, u64 repeat_count)
{
  FEAT_PRAGMA_OMP(parallel)
  for(u64 k = 0u; k < repeat_count; ++k)
  {
    FEAT_PRAGMA_OMP(for)
    for(u64 i = 0u; i < num_entries; ++i)
    {
      x[i] = y[i];
    }
  }
}

// memory write; must be declared NOINLINE to prevent the compiler from reordering instructions,
// which could move the operation out of the time measurement block
NOINLINE void write_mem(u64* x, u64 num_entries, u64 repeat_count)
{
  for(u64 k = 0u; k < repeat_count; ++k)
  {
    for(u64 i = 0u; i < num_entries; ++i)
    {
      // we have to perform some cheap operation here to prevent the compiler from removing
      // any of the two outer loops as part of overaggressive optimization attempts
      x[i] = i ^ k;
    }
  }
}

// memory read; must be declared NOINLINE to prevent the compiler from reordering instructions,
// which could move the operation out of the time measurement block
NOINLINE void read_mem(u64* x, u64 num_entries, u64 repeat_count)
{
  u64 t = u64(0);
  for(u64 k = 0u; k < repeat_count; ++k)
  {
    for(u64 i = 0u; i < num_entries; ++i)
    {
      // we have to perform some cheap operation here to prevent the compiler from removing
      // any of the two outer loops as part of overaggressive optimization attempts
      t ^= x[i];
    }
  }
  // prevent compiler from optimizing 't' away by storing its value in a volatile object
  volatile char foo = *((char*)&t);
  // prevent compiler from complaining about unused variables
  (void)foo;
}

// prints nicely padded output
void print_time(const Dist::Comm& comm, const String& name, const StopWatch& watch, u64 total_size)
{
  comm.print(name.pad_back(pad_len, '.') + ": " + stringify_bytes((1000000ull*total_size) / watch.elapsed_micros()) + "/sec");
}

int main(int argc, char** argv)
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  Dist::Comm comm = Dist::Comm::world();

  // default: 64 MB buffer and 100 repeat steps
  u64 buffer_size = 64*1024*1024; // 64 MB
  u64 repeat_count = 100;

  // check for help request
  for(int i = 1; i < argc; ++i)
  {
    bool need_help = false;
    need_help = need_help || (String(argv[i]).compare_no_case("?") == 0);
    need_help = need_help || (String(argv[i]).compare_no_case("-?") == 0);
    need_help = need_help || (String(argv[i]).compare_no_case("-help") == 0);
    need_help = need_help || (String(argv[i]).compare_no_case("--help") == 0);
    if(need_help)
    {
      comm.print("\nUSAGE: bandwidth-bench [buffer-size (in MB)] [repeat-count]\n");
      return 0;
    }
  }

  // parse buffer size if given
  if(argc > 1)
  {
    if(!String(argv[1]).parse(buffer_size))
    {
      comm.print(std::cerr, "\nERROR: Failed to parse '" + String(argv[1]) + "' as buffer size!\n");
      return 1;
    }
    buffer_size *= 1024*1024; // convert to MB
  }

  // parse repeat count if given
  if(argc > 2)
  {
    if(!String(argv[2]).parse(repeat_count))
    {
      comm.print(std::cerr, "\nERROR: Failed to parse '" + String(argv[2]) + "' as repeat count!\n");
      return 1;
    }
  }

  // print basic setup information
#ifdef FEAT_HAVE_MPI
  comm.print(String("MPI Process Count").pad_back(pad_len, '.') + ": " + stringify(comm.size()).pad_front(7));
#else
  comm.print(String("MPI Process Count").pad_back(pad_len, '.') + ":   -N/A-");
#endif
#ifdef FEAT_HAVE_OMP
  comm.print(String("OpenMP Thread Count").pad_back(pad_len, '.') + ": " + stringify(omp_get_max_threads()).pad_front(7));
#else
  comm.print(String("OpenMP Thread Count").pad_back(pad_len, '.') + ":   -N/A-");
#endif
  comm.print(String("Buffer Size").pad_back(pad_len, '.') + ": " + stringify_bytes(buffer_size));
  comm.print(String("Repeat Count").pad_back(pad_len, '.') + ": " + stringify(repeat_count).pad_front(7));

  const u64 num_entries = buffer_size / sizeof(u64);

  // allocate host memory via new
  u64* x = (u64*)malloc(buffer_size);
  u64* y = (u64*)malloc(buffer_size);

  // format buffers with random content
  Random rng;
  for(u64 i = 0u; i < num_entries; ++i)
    x[i] = y[i] = rng.next();

  // total number of bytes touched in each operation
  const u64 total_size = 2ull * u64(comm.size()) * repeat_count * buffer_size;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // measure host memory write bandwidth
  StopWatch watch_host_mem_write;
  comm.barrier();
  watch_host_mem_write.start();
  write_mem(x, num_entries, repeat_count);
  write_mem(y, num_entries, repeat_count);
  comm.barrier();
  watch_host_mem_write.stop();

  print_time(comm, "Host Memory Write", watch_host_mem_write, total_size);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // measure host memory read bandwidth
  StopWatch watch_host_mem_read;
  comm.barrier();
  watch_host_mem_read.start();
  read_mem(x, num_entries, repeat_count);
  read_mem(y, num_entries, repeat_count);
  comm.barrier();
  watch_host_mem_read.stop();

  print_time(comm, "Host Memory Read", watch_host_mem_read, total_size);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // measure host memcpy bandwidth
  StopWatch watch_host_memcpy;
  comm.barrier();
  watch_host_memcpy.start();
  for(u64 k = 0u; k < repeat_count; ++k)
  {
    memcpy(y, x, buffer_size);
  }
  comm.barrier();
  watch_host_memcpy.stop();

  print_time(comm, "Host Memcpy", watch_host_memcpy, total_size);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // measure host std::copy bandwidth
  StopWatch watch_host_stdcopy;
  watch_host_stdcopy.start();
  comm.barrier();
  for(u64 k = 0u; k < repeat_count; ++k)
  {
    std::copy(x, x + num_entries, y);
  }
  comm.barrier();
  watch_host_stdcopy.stop();

  print_time(comm, "Host std::copy", watch_host_stdcopy, total_size);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef FEAT_HAVE_OMP

  // measure host OpenMP copy bandwidth
  StopWatch watch_host_ompcpy;
  comm.barrier();
  watch_host_ompcpy.start();
  omp_copy(x, y, num_entries, repeat_count);
  comm.barrier();
  watch_host_ompcpy.stop();

  print_time(comm, "Host OpenMP Copy", watch_host_ompcpy, total_size);

#endif // FEAT_HAVE_OMP

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  typedef double DataType;
  typedef Index IndexType;

  Backend::set_preferred_backend(PreferredBackend::generic);

  const Index data_size = buffer_size / sizeof(DataType);
  LAFEM::DenseVector<DataType, IndexType> vec_x(data_size, DataType(0));
  LAFEM::DenseVector<DataType, IndexType> vec_y(data_size, DataType(0));

  // measure host DenseVector::format bandwidth

  StopWatch watch_host_dv_format;
  comm.barrier();
  watch_host_dv_format.start();
  for(u64 k = 0u; k < repeat_count; ++k)
  {
    vec_x.format(DataType(k));
    vec_y.format(DataType(k));
  }
  comm.barrier();
  watch_host_dv_format.stop();

  print_time(comm, "Host DenseVector::format", watch_host_dv_format, total_size);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // measure host DenseVector::copy bandwidth

  StopWatch watch_host_dv_copy;
  comm.barrier();
  watch_host_dv_copy.start();
  for(u64 k = 0u; k < repeat_count; ++k)
  {
    vec_y.copy(vec_x);
  }
  comm.barrier();
  watch_host_dv_copy.stop();

  print_time(comm, "Host DenseVector::copy", watch_host_dv_copy, total_size);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef FEAT_HAVE_CUDA

  u64* a = (u64*)Util::cuda_malloc(buffer_size);
  u64* b = (u64*)Util::cuda_malloc(buffer_size);

  // format memory to ensure that all pages are allocated
  Util::cuda_set_memory(a, u64(1), num_entries);
  Util::cuda_set_memory(b, u64(2), num_entries);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // measure CUDA host-to-device copy bandwidth

  StopWatch watch_cuda_h2dcpy;
  comm.barrier();
  watch_cuda_h2dcpy.start();
  for(u64 k = 0u; k < repeat_count; ++k)
  {
    Util::cuda_copy_host_to_device(a, x, buffer_size);
    Util::cuda_synchronize();
  }
  comm.barrier();
  watch_cuda_h2dcpy.stop();

  print_time(comm, "CUDA Host-2-Device Copy", watch_cuda_h2dcpy, total_size);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // measure CUDA device-to-device copy bandwidth

  StopWatch watch_cuda_d2dcpy;
  comm.barrier();
  watch_cuda_d2dcpy.start();
  for(u64 k = 0u; k < repeat_count; ++k)
  {
    Util::cuda_copy_device_to_device(b, a, buffer_size);
    Util::cuda_synchronize();
  }
  comm.barrier();
  watch_cuda_d2dcpy.stop();

  print_time(comm, "CUDA Device-2-Device Copy", watch_cuda_d2dcpy, total_size);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // measure CUDA device-to-host copy bandwidth

  StopWatch watch_cuda_d2hcpy;
  comm.barrier();
  watch_cuda_d2hcpy.start();
  for(u64 k = 0u; k < repeat_count; ++k)
  {
    Util::cuda_copy_device_to_host(y, b, buffer_size);
    Util::cuda_synchronize();
  }
  comm.barrier();
  watch_cuda_d2hcpy.stop();

  print_time(comm, "CUDA Device-2-Host Copy", watch_cuda_d2hcpy, total_size);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  Backend::set_preferred_backend(PreferredBackend::cuda);

  // format vectors to move memory to device
  vec_x.format(DataType(0));
  vec_y.format(DataType(0));

  // measure CUDA DenseVector::format bandwidth

  StopWatch watch_cuda_dv_format;
  comm.barrier();
  watch_cuda_dv_format.start();
  for(u64 k = 0u; k < repeat_count; ++k)
  {
    vec_x.format(DataType(k));
    vec_y.format(DataType(k));
  }
  comm.barrier();
  watch_cuda_dv_format.stop();

  print_time(comm, "CUDA DenseVector::format", watch_cuda_dv_format, total_size);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // measure CUDA DenseVector::copy bandwidth

  StopWatch watch_cuda_dv_copy;
  comm.barrier();
  watch_cuda_dv_copy.start();
  for(u64 k = 0u; k < repeat_count; ++k)
  {
    vec_y.copy(vec_x);
  }
  comm.barrier();
  watch_cuda_dv_copy.stop();

  print_time(comm, "CUDA DenseVector::copy", watch_cuda_dv_copy, total_size);

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  Util::cuda_free(b);
  Util::cuda_free(a);

#endif // FEAT_HAVE_CUDA

  free(y);
  free(x);

  return 0;
}
