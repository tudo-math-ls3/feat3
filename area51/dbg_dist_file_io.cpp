// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/util/runtime.hpp>
#include <kernel/util/random.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/dist_file_io.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/string.hpp>

#include <cstring>

namespace DbgDistFileIO
{
  using namespace FEAT;

  char cx(int c)
  {
    return char(c < 10 ? '0' + c : 'A' + (c-10));
  }

  String dumpbin(const std::vector<char>& v)
  {
    String s;
    s.reserve(3*v.size());
    for(auto c : v)
    {
      s.push_back(' ');
      s.push_back(cx((c >> 4) & 0xF));
      s.push_back(cx(c & 0xF));
    }
    return s;
  }

  bool vcomp(const std::vector<char>& x, const std::vector<char>& y)
  {
    if(x.size() != y.size())
      return false;

    for(std::size_t i(0); i < x.size(); ++i)
    {
      if(x[i] != y[i])
        return false;
    }

    return true;
  }

  void test_combined(Dist::Comm& comm)
  {
    std::vector<char> shared1, buffer1, shared2, buffer2;

    //                         123456789-1234567
    const String str_shared = "This isn't a Test";
    const String filename = "dist_file_io_combined.tmp";
    const std::uint64_t my_seed = 0xDEADBEEFC0DEBABEull ^ ((unsigned)comm.rank() * 17ull);
    const std::size_t my_size = 4*(1 + (comm.rank() % 5));

    // initialize shared
    shared1.resize(str_shared.size());
    memcpy(shared1.data(), str_shared.c_str(), str_shared.size());

    // initialize buffer with random values
    {
      Random rng(my_seed);
      buffer1.resize(my_size);
      char* x = buffer1.data();
      for(std::size_t i(0); i < my_size; ++i)
        x[i] = (char)rng(0, 255);
    }

    // write via DistFileIO
    DistFileIO::write_combined(shared1, buffer1, filename, comm);

    // read via DistFileIO
    DistFileIO::read_combined(shared2, buffer2, filename, comm);

    // compare shared size
    if(shared1.size() != shared2.size())
    {
      std::cerr << "ERROR on rank " << comm.rank() << ": shared size mismatch: expected "
        << shared1.size() << " but got " << buffer2.size() << std::endl;
    }

    // compare buffer size
    if(buffer2.size() != my_size)
    {
      std::cerr << "ERROR on rank " << comm.rank() << ": buffer size mismatch: expected "
        << my_size << " but got " << buffer2.size() << std::endl;
      Runtime::abort();
    }

    // dump shared buffer
    comm.print("\nShared Buffer:");
    String sh;
    if(comm.rank() == 0)
      (sh += dumpbin(shared1)) += " vs\n";
    sh += dumpbin(shared2);
    sh += (vcomp(shared1, shared2) ? " OK" : " <-- FAILED");
    comm.allprint(sh);

    // dump individual buffers
    comm.print("\nRank Buffers:");
    String bh;
    bh += dumpbin(buffer1) + " vs\n";
    bh += dumpbin(buffer2);
    bh += (vcomp(buffer1, buffer2) ? " OK" : " <-- FAILED");
    comm.allprint(bh);

    // compare shared
    for(std::size_t i(0); i < shared1.size(); ++i)
    {
      if(shared1[i] != shared1[i])
      {
        std::cerr << "ERROR on rank " << comm.rank() << ": shared data mismatch at byte " << i
          << ": expected " << int(shared1[i]) << " but got " << int(shared2[i]) << std::endl;
        Runtime::abort();
      }
    }

    // compare buffer data
    for(std::size_t i(0); i < my_size; ++i)
    {
      if(buffer1[i] != buffer2[i])
      {
        std::cerr << "ERROR on rank " << comm.rank() << ": buffer data mismatch at byte " << i
          << ": expected " << int(buffer1[i]) << " but got " << int(buffer2[i]) << std::endl;
        Runtime::abort();
      }
    }

    // okay
    comm.barrier();
    comm.print("OK");
  }

  void main(/*int argc, char** argv*/)
  {
    Dist::Comm comm = Dist::Comm::world();

#ifdef FEAT_HAVE_MPI
    comm.print("Have MPI.: yes");
#else
    comm.print("Have MPI.: no");
#endif
    comm.print("Processes: " + stringify(comm.size()));

    test_combined(comm);
  }
}

int main(int argc, char** argv)
{
  FEAT::Runtime::initialize(argc, argv);
  DbgDistFileIO::main(/*argc, argv*/);
  return FEAT::Runtime::finalize();
}
