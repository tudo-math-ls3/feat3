// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/runtime.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/kahan_accumulator.hpp>
#include <kernel/util/likwid_marker.hpp>
#include <kernel/util/memory_usage.hpp>
#include <kernel/util/stop_watch.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/memory_aux.hpp>
#include <memory>

using namespace FEAT;

void bench_send_receive(Dist::Comm& comm, SimpleArgParser& args)
{
  TimeStamp stamp_0;

  // approximate 4th root of num ranks
  const int num_ranks = comm.size();
  const int my_rank = comm.rank();
  int dir_rank = 2;
  while(dir_rank*dir_rank*dir_rank*dir_rank < num_ranks)
    ++dir_rank;

  Index max_size_mb = 1024; // in MB
  Index num_steps = 5;
  Index num_repeats = 3;

  if(args.parse("sendrecv", max_size_mb, num_steps, num_repeats) < 0)
  {
    comm.print("ERROR: invalid argument for --sendrecv");
    Runtime::abort();
  }

  // compute neighbor ranks (up to 8)
  std::vector<int> neighbors;
  auto add_rank = [&neighbors,my_rank, num_ranks](int i)
    {if((my_rank+i >= 0) && (my_rank+i < num_ranks)) neighbors.push_back(my_rank+i);};

  add_rank(-dir_rank*dir_rank*dir_rank);
  add_rank(-dir_rank*dir_rank);
  add_rank(-dir_rank);
  add_rank(-1);
  add_rank(1);
  add_rank(dir_rank);
  add_rank(dir_rank*dir_rank);
  add_rank(dir_rank*dir_rank*dir_rank);

  // number of neighbors
  const std::size_t nn = neighbors.size();

  TimeStamp stamp_1;

  // create and format buffers to perform first touch of memory
  const std::size_t buf_size = std::size_t(max_size_mb) * 1024u*64u; // 16 buffers in total
  std::vector<std::vector<char>> send_bufs, recv_bufs;
  send_bufs.resize(nn);
  recv_bufs.resize(nn);
  for(std::size_t i = 0; i < nn; ++i)
  {
    send_bufs.at(i).resize(buf_size);
    recv_bufs.at(i).resize(buf_size);
    Memory::memset_random_main(send_bufs.at(i).data(), buf_size, i+1u);
    Memory::memset_main(recv_bufs.at(i).data(), 0, buf_size);
  }

  Dist::RequestVector send_reqs(nn), recv_reqs(nn);

  TimeStamp stamp_2;

  if(my_rank == 0)
  {
    std::cout << String(120, '=') << "\n";
    std::cout << "Performing send-receive benchmark...\n";
    std::cout << "Total MPI ranks: " << num_ranks << "\n";
    std::cout << "Directional distance: " << dir_rank << "\n";
    std::cout << "Maximum memory buffer size: " << stringify_bytes(16u*buf_size) << "\n";
    std::cout << "Total memory usage: " << stringify_bytes(16u*buf_size*num_ranks) << "\n";
    std::cout << "Coarsening steps: " << num_steps << "\n";
    std::cout << "Number of repeats: " << num_repeats << "\n";
    std::cout << "\n";
    //std::cout << "Step    Buffer-Size    Min-Time    Max-Time    Avg-Time    Avg-Bandwidth        Latency\n";
    std::cout << "Step    Buffer-Size    Min-Time    Max-Time    Avg-Time    Avg-Bandwidth\n";
    std::cout.flush();
  }

  std::size_t step_buf_size = 1;

  //std::vector<double> time_min, time_max, time_avg, size_buf;

  for(Index step = 0; step < num_steps; ++step)
  {
    if(my_rank == 0)
    {
      //std::cout << "Step " << (step+1) << ", buffer size " << step_buf_size << "\n";
      std::cout << stringify(step).pad_front(4);
      //std::cout << stringify_bytes(step_buf_size).pad_front(12);
      std::cout << stringify(step_buf_size).pad_front(15);
    }

    //size_buf.push_back(double(step_buf_size));

    comm.barrier();
    TimeStamp stamp_t1;

    for(Index rep = 0; rep < num_repeats; ++rep)
    {
      for(std::size_t i(0); i < nn; ++i)
      {
        recv_reqs[i] = comm.irecv(recv_bufs[i].data(), step_buf_size, neighbors[i]);
      }
      for(std::size_t i(0); i < nn; ++i)
      {
        send_reqs[i] = comm.isend(send_bufs[i].data(), step_buf_size, neighbors[i]);
      }
      recv_reqs.wait_all();
      send_reqs.wait_all();
    }

    comm.barrier();
    TimeStamp stamp_t2;

    const double t = stamp_t2.elapsed(stamp_t1);

    // reduce runtimes
    double t_min = t;
    double t_max = t;
    double t_avg = t;
    comm.reduce(&t, &t_min, std::size_t(1), Dist::op_min, 0);
    comm.reduce(&t, &t_max, std::size_t(1), Dist::op_max, 0);
    comm.reduce(&t, &t_avg, std::size_t(1), Dist::op_sum, 0);
    t_min /= double(num_repeats);
    t_max /= double(num_repeats);
    t_avg /= double(num_repeats) * double(num_ranks);
    //time_min.push_back(t_min);
    //time_max.push_back(t_max);
    //time_avg.push_back(t_avg);

    if(my_rank == 0)
    {
      std::size_t bw = std::size_t(double(num_repeats)*double(step_buf_size) / t_avg);
      std::cout << stringify_fp_fix(t_min, 6, 12);
      std::cout << stringify_fp_fix(t_max, 6, 12);
      std::cout << stringify_fp_fix(t_avg, 6, 12);
      if(step == 0)
        std::cout << "              ---";
      else
        std::cout << stringify_bytes(bw).pad_front(15) << "/s";
      std::cout << "\n";
    }

    if(step == 0)
      step_buf_size = buf_size;
    else
      step_buf_size /= 2u;
  }

  TimeStamp stamp_3;
  if(my_rank == 0)
  {
    std::cout << "\nBenchmark Runtime: " << stamp_3.elapsed_string(stamp_0) << "\n";
    std::cout.flush();
  }
}

void bench_broadcast(Dist::Comm& comm, SimpleArgParser& args)
{
  TimeStamp stamp_0;

  const int num_ranks = comm.size();
  const int my_rank = comm.rank();

  Index max_size_mb = 128; // in MB
  Index num_steps = 5;
  Index num_repeats = 3;

  if(args.parse("broadcast", max_size_mb, num_steps, num_repeats) < 0)
  {
    comm.print("ERROR: invalid argument for --broadcast");
    Runtime::abort();
  }

  TimeStamp stamp_1;

  // create and format buffers to perform first touch of memory
  const std::size_t buf_size = std::size_t(max_size_mb) * 1024u*1024;
  std::vector<char> buffer(buf_size);
  Memory::memset_random_main(buffer.data(), buf_size);

  TimeStamp stamp_2;

  if(my_rank == 0)
  {
    std::cout << String(120, '=') << "\n";
    std::cout << "Performing broadcast benchmark...\n";
    std::cout << "Total MPI ranks: " << num_ranks << "\n";
    std::cout << "Maximum memory buffer size: " << stringify_bytes(buf_size) << "\n";
    std::cout << "Total memory usage: " << stringify_bytes(buf_size*num_ranks) << "\n";
    std::cout << "Coarsening steps: " << num_steps << "\n";
    std::cout << "Number of repeats: " << num_repeats << "\n";
    std::cout << "\n";
    //std::cout << "Step    Buffer-Size    Min-Time    Max-Time    Avg-Time    Avg-Bandwidth        Latency\n";
    std::cout << "Step    Buffer-Size    Min-Time    Max-Time    Avg-Time    Avg-Bandwidth\n";
    std::cout.flush();
  }

  std::size_t step_buf_size = 1;

  for(Index step = 0; step < num_steps; ++step)
  {
    if(my_rank == 0)
    {
      //std::cout << "Step " << (step+1) << ", buffer size " << step_buf_size << "\n";
      std::cout << stringify(step).pad_front(4);
      //std::cout << stringify_bytes(step_buf_size).pad_front(12);
      std::cout << stringify(step_buf_size).pad_front(15);
    }

    comm.barrier();
    TimeStamp stamp_t1;

    for(Index rep = 0; rep < num_repeats; ++rep)
    {
      comm.bcast(buffer.data(), step_buf_size, Dist::dt_byte, 0);
    }

    comm.barrier();
    TimeStamp stamp_t2;

    const double t = stamp_t2.elapsed(stamp_t1);

    double t_min = t;
    double t_max = t;
    double t_avg = t;
    comm.reduce(&t, &t_min, std::size_t(1), Dist::op_min, 0);
    comm.reduce(&t, &t_max, std::size_t(1), Dist::op_max, 0);
    comm.reduce(&t, &t_avg, std::size_t(1), Dist::op_sum, 0);
    t_min /= double(num_repeats);
    t_max /= double(num_repeats);
    t_avg /= double(num_repeats) * double(num_ranks);

    if(my_rank == 0)
    {
      std::size_t bw = std::size_t(double(num_repeats)*double(step_buf_size) / t_avg);
      std::cout << stringify_fp_fix(t_min, 6, 12);
      std::cout << stringify_fp_fix(t_max, 6, 12);
      std::cout << stringify_fp_fix(t_avg, 6, 12);
      if(step == 0)
        std::cout << "              ---";
      else
        std::cout << stringify_bytes(bw).pad_front(15) << "/s";
      std::cout << "\n";
      std::cout.flush();
    }

    if(step == 0)
      step_buf_size = buf_size;
    else
      step_buf_size /= 2u;
  }

  TimeStamp stamp_3;
  if(my_rank == 0)
  {
    std::cout << "\nBenchmark Runtime: " << stamp_3.elapsed_string(stamp_0) << "\n";
    std::cout.flush();
  }
}

void bench_reduce(Dist::Comm& comm, SimpleArgParser& args)
{
  TimeStamp stamp_0;

  const int num_ranks = comm.size();
  const int my_rank = comm.rank();

  Index max_size_mb = 128; // in MB
  Index num_steps = 5;
  Index num_repeats = 3;

  if(args.parse("reduce", max_size_mb, num_steps, num_repeats) < 0)
  {
    comm.print("ERROR: invalid argument for --reduce");
    Runtime::abort();
  }

  TimeStamp stamp_1;

  // create and format buffers to perform first touch of memory
  const std::size_t buf_size = std::size_t(max_size_mb) * 1024u*1024;
  std::vector<char> buffer(buf_size);
  Memory::memset_random_main(buffer.data(), buf_size);

  TimeStamp stamp_2;

  if(my_rank == 0)
  {
    std::cout << String(120, '=') << "\n";
    std::cout << "Performing reduce benchmark...\n";
    std::cout << "Total MPI ranks: " << num_ranks << "\n";
    std::cout << "Maximum memory buffer size: " << stringify_bytes(buf_size) << "\n";
    std::cout << "Total memory usage: " << stringify_bytes(buf_size*num_ranks) << "\n";
    std::cout << "Coarsening steps: " << num_steps << "\n";
    std::cout << "Number of repeats: " << num_repeats << "\n";
    std::cout << "\n";
    //std::cout << "Step    Buffer-Size    Min-Time    Max-Time    Avg-Time    Avg-Bandwidth        Latency\n";
    std::cout << "Step    Buffer-Size    Min-Time    Max-Time    Avg-Time    Avg-Bandwidth\n";
    std::cout.flush();
  }

  std::size_t step_buf_size = 1;

  Dist::Operation op_bxor(MPI_BXOR);

  for(Index step = 0; step < num_steps; ++step)
  {
    if(my_rank == 0)
    {
      //std::cout << "Step " << (step+1) << ", buffer size " << step_buf_size << "\n";
      std::cout << stringify(step).pad_front(4);
      //std::cout << stringify_bytes(step_buf_size).pad_front(12);
      std::cout << stringify(step_buf_size).pad_front(15);
    }

    comm.barrier();
    TimeStamp stamp_t1;

    for(Index rep = 0; rep < num_repeats; ++rep)
    {
      comm.reduce(buffer.data(), buffer.data(), step_buf_size, Dist::dt_byte, op_bxor, 0);
    }

    comm.barrier();
    TimeStamp stamp_t2;

    const double t = stamp_t2.elapsed(stamp_t1);

    double t_min = t;
    double t_max = t;
    double t_avg = t;
    comm.reduce(&t, &t_min, std::size_t(1), Dist::op_min, 0);
    comm.reduce(&t, &t_max, std::size_t(1), Dist::op_max, 0);
    comm.reduce(&t, &t_avg, std::size_t(1), Dist::op_sum, 0);
    t_min /= double(num_repeats);
    t_max /= double(num_repeats);
    t_avg /= double(num_repeats) * double(num_ranks);

    if(my_rank == 0)
    {
      std::size_t bw = std::size_t(double(num_repeats)*double(step_buf_size) / t_avg);
      std::cout << stringify_fp_fix(t_min, 6, 12);
      std::cout << stringify_fp_fix(t_max, 6, 12);
      std::cout << stringify_fp_fix(t_avg, 6, 12);
      if(step == 0)
        std::cout << "              ---";
      else
        std::cout << stringify_bytes(bw).pad_front(15) << "/s";
      std::cout << "\n";
      std::cout.flush();
    }

    if(step == 0)
      step_buf_size = buf_size;
    else
      step_buf_size /= 2u;
  }

  TimeStamp stamp_3;
  if(my_rank == 0)
  {
    std::cout << "\nBenchmark Runtime: " << stamp_3.elapsed_string(stamp_0) << "\n";
    std::cout.flush();
  }
}

void bench_gather(Dist::Comm& comm, SimpleArgParser& args)
{
  TimeStamp stamp_0;

  const int num_ranks = comm.size();
  const int my_rank = comm.rank();

  Index max_size_mb = 64; // in MB
  Index num_steps = 5;
  Index num_repeats = 3;

  if(args.parse("gather", max_size_mb, num_steps, num_repeats) < 0)
  {
    comm.print("ERROR: invalid argument for --gather");
    Runtime::abort();
  }

  TimeStamp stamp_1;

  // create and format buffers to perform first touch of memory
  const std::size_t send_buf_size = std::size_t(max_size_mb) * 1024u*1024;
  const std::size_t recv_buf_size = send_buf_size * std::size_t(num_ranks);
  std::vector<char> send_buffer(send_buf_size), recv_buffer;
  Memory::memset_random_main(send_buffer.data(), send_buf_size, my_rank+1);
  if(my_rank == 0)
  {
    recv_buffer.resize(recv_buf_size);
    Memory::memset_main(recv_buffer.data(), 0, recv_buffer.size());
  }

  TimeStamp stamp_2;

  if(my_rank == 0)
  {
    std::cout << String(120, '=') << "\n";
    std::cout << "Performing gather benchmark...\n";
    std::cout << "Total MPI ranks: " << num_ranks << "\n";
    std::cout << "Send buffer size: " << stringify_bytes(send_buf_size) << "\n";
    std::cout << "Receive buffer size: " << stringify_bytes(recv_buf_size) << "\n";
    std::cout << "Total memory usage: " << stringify_bytes(2*recv_buf_size) << "\n";
    std::cout << "Coarsening steps: " << num_steps << "\n";
    std::cout << "Number of repeats: " << num_repeats << "\n";
    std::cout << "\n";
    //std::cout << "Step    Buffer-Size    Min-Time    Max-Time    Avg-Time    Avg-Bandwidth        Latency\n";
    std::cout << "Step    Buffer-Size    Min-Time    Max-Time    Avg-Time    Avg-Bandwidth\n";
    std::cout.flush();
  }

  std::size_t step_buf_size = 1;

  for(Index step = 0; step < num_steps; ++step)
  {
    if(my_rank == 0)
    {
      //std::cout << "Step " << (step+1) << ", buffer size " << step_buf_size << "\n";
      std::cout << stringify(step).pad_front(4);
      //std::cout << stringify_bytes(step_buf_size).pad_front(12);
      std::cout << stringify(step_buf_size).pad_front(15);
    }

    comm.barrier();
    TimeStamp stamp_t1;

    for(Index rep = 0; rep < num_repeats; ++rep)
    {
      comm.gather(send_buffer.data(), step_buf_size, recv_buffer.data(), step_buf_size, 0);
    }

    comm.barrier();
    TimeStamp stamp_t2;

    const double t = stamp_t2.elapsed(stamp_t1);

    double t_min = t;
    double t_max = t;
    double t_avg = t;
    comm.reduce(&t, &t_min, std::size_t(1), Dist::op_min, 0);
    comm.reduce(&t, &t_max, std::size_t(1), Dist::op_max, 0);
    comm.reduce(&t, &t_avg, std::size_t(1), Dist::op_sum, 0);
    t_min /= double(num_repeats);
    t_max /= double(num_repeats);
    t_avg /= double(num_repeats) * double(num_ranks);

    if(my_rank == 0)
    {
      std::size_t bw = std::size_t(double(num_repeats)*double(step_buf_size) / t_avg);
      std::cout << stringify_fp_fix(t_min, 6, 12);
      std::cout << stringify_fp_fix(t_max, 6, 12);
      std::cout << stringify_fp_fix(t_avg, 6, 12);
      if(step == 0)
        std::cout << "              ---";
      else
        std::cout << stringify_bytes(bw).pad_front(15) << "/s";
      std::cout << "\n";
      std::cout.flush();
    }

    if(step == 0)
      step_buf_size = send_buf_size;
    else
      step_buf_size /= 2u;
  }

  TimeStamp stamp_3;
  if(my_rank == 0)
  {
    std::cout << "\nBenchmark Runtime: " << stamp_3.elapsed_string(stamp_0) << "\n";
    std::cout.flush();
  }
}

void bench_scatter(Dist::Comm& comm, SimpleArgParser& args)
{
  TimeStamp stamp_0;

  const int num_ranks = comm.size();
  const int my_rank = comm.rank();

  Index max_size_mb = 64; // in MB
  Index num_steps = 5;
  Index num_repeats = 3;

  if(args.parse("scatter", max_size_mb, num_steps, num_repeats) < 0)
  {
    comm.print("ERROR: invalid argument for --scatter");
    Runtime::abort();
  }

  TimeStamp stamp_1;

  // create and format buffers to perform first touch of memory
  const std::size_t recv_buf_size = std::size_t(max_size_mb) * 1024u*1024;
  const std::size_t send_buf_size = recv_buf_size * std::size_t(num_ranks);
  std::vector<char> recv_buffer(recv_buf_size), send_buffer;
  Memory::memset_random_main(recv_buffer.data(), recv_buf_size, my_rank+1);
  if(my_rank == 0)
  {
    send_buffer.resize(send_buf_size);
    Memory::memset_main(send_buffer.data(), 0, send_buffer.size());
  }

  TimeStamp stamp_2;

  if(my_rank == 0)
  {
    std::cout << String(120, '=') << "\n";
    std::cout << "Performing scatter benchmark...\n";
    std::cout << "Total MPI ranks: " << num_ranks << "\n";
    std::cout << "Receive buffer size: " << stringify_bytes(recv_buf_size) << "\n";
    std::cout << "Send buffer size: " << stringify_bytes(send_buf_size) << "\n";
    std::cout << "Total memory usage: " << stringify_bytes(2*send_buf_size) << "\n";
    std::cout << "Coarsening steps: " << num_steps << "\n";
    std::cout << "Number of repeats: " << num_repeats << "\n";
    std::cout << "\n";
    //std::cout << "Step    Buffer-Size    Min-Time    Max-Time    Avg-Time    Avg-Bandwidth        Latency\n";
    std::cout << "Step    Buffer-Size    Min-Time    Max-Time    Avg-Time    Avg-Bandwidth\n";
    std::cout.flush();
  }

  std::size_t step_buf_size = 1;

  for(Index step = 0; step < num_steps; ++step)
  {
    if(my_rank == 0)
    {
      //std::cout << "Step " << (step+1) << ", buffer size " << step_buf_size << "\n";
      std::cout << stringify(step).pad_front(4);
      //std::cout << stringify_bytes(step_buf_size).pad_front(12);
      std::cout << stringify(step_buf_size).pad_front(15);
    }

    comm.barrier();
    TimeStamp stamp_t1;

    for(Index rep = 0; rep < num_repeats; ++rep)
    {
      comm.scatter(send_buffer.data(), step_buf_size, recv_buffer.data(), step_buf_size, 0);
    }

    comm.barrier();
    TimeStamp stamp_t2;

    const double t = stamp_t2.elapsed(stamp_t1);

    double t_min = t;
    double t_max = t;
    double t_avg = t;
    comm.reduce(&t, &t_min, std::size_t(1), Dist::op_min, 0);
    comm.reduce(&t, &t_max, std::size_t(1), Dist::op_max, 0);
    comm.reduce(&t, &t_avg, std::size_t(1), Dist::op_sum, 0);
    t_min /= double(num_repeats);
    t_max /= double(num_repeats);
    t_avg /= double(num_repeats) * double(num_ranks);

    if(my_rank == 0)
    {
      std::size_t bw = std::size_t(double(num_repeats)*double(step_buf_size) / t_avg);
      std::cout << stringify_fp_fix(t_min, 6, 12);
      std::cout << stringify_fp_fix(t_max, 6, 12);
      std::cout << stringify_fp_fix(t_avg, 6, 12);
      if(step == 0)
        std::cout << "              ---";
      else
        std::cout << stringify_bytes(bw).pad_front(15) << "/s";
      std::cout << "\n";
      std::cout.flush();
    }

    if(step == 0)
      step_buf_size = recv_buf_size;
    else
      step_buf_size /= 2u;
  }

  TimeStamp stamp_3;
  if(my_rank == 0)
  {
    std::cout << "\nBenchmark Runtime: " << stamp_3.elapsed_string(stamp_0) << "\n";
    std::cout.flush();
  }
}

int main(int argc, char** argv)
{
  Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  Dist::Comm comm(Dist::Comm::world());

  SimpleArgParser args(argc, argv);
  args.support("sendrecv", "<max_mem> [<num_steps>] [<repeat>]\n");
  args.support("broadcast", "<max_mem> [<num_steps>] [<repeat>]\n");
  args.support("reduce", "<max_mem> [<num_steps>] [<repeat>]\n");
  args.support("gather", "<max_mem> [<num_steps>] [<repeat>]\n");
  args.support("scatter", "<max_mem> [<num_steps>] [<repeat>]\n");

  if(args.num_args() < 2)
  {
    comm.print(args.get_supported_help());
    return 0;
  }

  if(args.check("sendrecv") >= 0)
    bench_send_receive(comm, args);

  if(args.check("broadcast") >= 0)
    bench_broadcast(comm, args);

  if(args.check("reduce") >= 0)
    bench_reduce(comm, args);

  if(args.check("gather") >= 0)
    bench_gather(comm, args);

  if(args.check("scatter") >= 0)
    bench_scatter(comm, args);

  return 0;
}
