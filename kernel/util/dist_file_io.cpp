// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/util/assertion.hpp>
#include <kernel/util/dist_file_io.hpp>

#include <iostream>
#include <fstream>

namespace FEAT
{
  String DistFileIO::_rankname(const String& pattern, int rank)
  {
    // find the first asterisk in the filename
    const std::size_t p = pattern.find_first_of('*');
    XASSERTM(p != pattern.npos, "sequence filename template is missing rank wildcat pattern");

    // find the first character after the pattern
    const std::size_t q = pattern.find_first_not_of('*', p);

    // compute wildcat pattern length
    const std::size_t n = (q != pattern.npos ? q - p : pattern.size() - p);

    // build the filename of our rank
    String filename(pattern.substr(std::size_t(0), p));
    filename += stringify(rank).pad_front(n, '0');
    if(q != pattern.npos)
      filename += pattern.substr(q);

    // that's it
    return filename;
  }

  void DistFileIO::_read_file(std::stringstream& stream, const String& filename)
  {
    // open input file
    std::ifstream ifs(filename, std::ios_base::in);
    if(!ifs.is_open() || !ifs.good())
      throw FileNotFound(filename);

    // read into our stream buffer
    stream << ifs.rdbuf();

    // close file
    ifs.close();

    // sanity check
    XASSERT(stream.good());
  }

  void DistFileIO::_read_file(BinaryStream& stream, const String& filename)
  {
    // open input file
    std::ifstream ifs(filename, std::ios_base::in|std::ios_base::binary);
    if(!ifs.is_open() || !ifs.good())
      throw FileNotFound(filename);

    // read our input stream
    stream.read_stream(ifs);

    // close file
    ifs.close();

    // sanity check
    XASSERT(stream.good());
  }

  void DistFileIO::_write_file(std::stringstream& stream, const String& filename, bool truncate)
  {
    // determine output mode
    std::ios_base::openmode mode = std::ios_base::out;
    if(truncate)
      mode |= std::ios_base::trunc;

    // open output file
    std::ofstream ofs(filename, mode);
    if(!ofs.is_open() || !ofs.good())
      throw FileNotCreated(filename);

    // write stream
    ofs << stream.rdbuf();

    // close file
    ofs.close();
  }

  void DistFileIO::_write_file(BinaryStream& stream, const String& filename, bool truncate)
  {
    // determine output mode
    std::ios_base::openmode mode = std::ios_base::out|std::ios_base::binary;
    if(truncate)
      mode |= std::ios_base::trunc;

    // open output file
    std::ofstream ofs(filename, mode);
    if(!ofs.is_open() || !ofs.good())
      throw FileNotCreated(filename);

    // write buffer
    stream.write_stream(ofs);

    // close file
    ofs.close();
  }

#ifdef FEAT_HAVE_MPI

  void DistFileIO::read_common(std::stringstream& stream, const String& filename, const Dist::Comm& comm, int root_rank)
  {
    // root rank reads the file
    if(comm.rank() == root_rank)
    {
      _read_file(stream, filename);
    }

    // broadcast
    comm.bcast_stringstream(stream, root_rank);

    // sanity check
    XASSERT(stream.good());
  }

  void DistFileIO::read_common(BinaryStream& stream, const String& filename, const Dist::Comm& comm, int root_rank)
  {
    // root rank reads the file
    if(comm.rank() == root_rank)
    {
      _read_file(stream, filename);
    }

    // broadcast
    comm.bcast_binarystream(stream, root_rank);

    // seek to beginning
    stream.seekg(std::streamoff(0), std::ios_base::beg);

    // sanity check
    XASSERT(stream.good());
  }

  void DistFileIO::read_sequence(std::stringstream& stream, const String& pattern, const Dist::Comm& comm)
  {
    // retrieve our rank and nprocs
    int rank = comm.rank();
    int nprocs = comm.size();

    // build our filename
    String filename = _rankname(pattern, rank);

    // round-robin we go
    for(int i(0); i < nprocs; ++i)
    {
      // is it our turn?
      if(i == rank)
      {
        _read_file(stream, filename);
      }
      comm.barrier();
    }
  }

  void DistFileIO::read_sequence(BinaryStream& stream, const String& pattern, const Dist::Comm& comm)
  {
    // retrieve our rank and nprocs
    int rank = comm.rank();
    int nprocs = comm.size();

    // build our filename
    String filename = _rankname(pattern, rank);

    // round-robin we go
    for(int i(0); i < nprocs; ++i)
    {
      // is it our turn?
      if(i == rank)
      {
        _read_file(stream, filename);
      }
      comm.barrier();
    }
  }

  void DistFileIO::write_sequence(std::stringstream& stream, const String& pattern, const Dist::Comm& comm, bool truncate)
  {
    // retrieve our rank and nprocs
    int rank = comm.rank();
    int nprocs = comm.size();

    // build our filename
    String filename = _rankname(pattern, rank);

    // round-robin we go
    for(int i(0); i < nprocs; ++i)
    {
      // is it our turn?
      if(i == rank)
      {
        _write_file(stream, filename, truncate);
      }
      comm.barrier();
    }
  }

  void DistFileIO::write_sequence(BinaryStream& stream, const String& pattern, const Dist::Comm& comm, bool truncate)
  {
    // retrieve our rank and nprocs
    int rank = comm.rank();
    int nprocs = comm.size();

    // build our filename
    String filename = _rankname(pattern, rank);

    // round-robin we go
    for(int i(0); i < nprocs; ++i)
    {
      // is it our turn?
      if(i == rank)
      {
        _write_file(stream, filename, truncate);
      }
      comm.barrier();
    }
  }

  void DistFileIO::read_ordered(void* buffer, const std::size_t size, const String& filename, const Dist::Comm& comm)
  {
    XASSERT((buffer != nullptr) || (size == std::size_t(0)));

    // open file
    MPI_Status status;
    MPI_File file = MPI_FILE_NULL;
    MPI_File_open(comm.mpi_comm(), filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
    XASSERTM(file != MPI_FILE_NULL, "failed to open file via MPI_File_open");

    // read buffer via collective
    MPI_File_read_ordered(file, buffer, int(size), MPI_BYTE, &status);

    // close file
    MPI_File_close(&file);
  }

  void DistFileIO::write_ordered(const void* buffer, const std::size_t size, const String& filename, const Dist::Comm& comm, bool truncate)
  {
    XASSERT((buffer != nullptr) || (size == std::size_t(0)));

    // select file access mode
    const int amode = MPI_MODE_WRONLY | MPI_MODE_CREATE;

    // open file
    MPI_Status status;
    MPI_File file = MPI_FILE_NULL;
    MPI_File_open(comm.mpi_comm(), filename.c_str(), amode, MPI_INFO_NULL, &file);
    XASSERTM(file != MPI_FILE_NULL, "failed to open file via MPI_File_open");

    // truncate file?
    if(truncate)
    {
      MPI_File_set_size(file, MPI_Offset(0));
    }

    // write buffer via collective
    MPI_File_write_ordered(file, buffer, int(size), MPI_BYTE, &status);

    // close file
    MPI_File_close(&file);
  }

  void DistFileIO::read_combined(std::vector<char>& common, std::vector<char>& buffer,
    const String& filename, const Dist::Comm& comm, int root_rank, bool bcast_common)
  {
    typedef std::uint64_t u64;

    XASSERTM((root_rank >= 0) && (root_rank < comm.size()), "invalid root rank");

    std::vector<char> header(std::size_t(32), 0);//, padding(std::size_t(64), 0);
    u64* head_u64 = reinterpret_cast<u64*>(header.data());

    // open file
    MPI_Status status;
    MPI_File file = MPI_FILE_NULL;
    MPI_File_open(comm.mpi_comm(), filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
    XASSERTM(file != MPI_FILE_NULL, "failed to open file via MPI_File_open");

    // read header on root rank
    if(comm.rank() == root_rank)
    {
      // read header
      MPI_File_read_shared(file, header.data(), int(header.size()), MPI_BYTE, &status);

      // check file magic
      String msg1 = String("input file is not a valid FEAT3 file: ") + filename;
      XASSERTM(head_u64[0] == magic_combined, msg1.c_str());

      // check number of processes
      String msg2 = String("invalid number of processes: comm.size() is ") + stringify(comm.size()) +
        " but file contains data for " + stringify(head_u64[2]) + " processes";
      XASSERTM(head_u64[2] == u64(comm.size()), msg2.c_str());
    }

    // broadcast header from root
    comm.bcast(header.data(), header.size(), root_rank);

    // compute shared buffer padding size
    const u64 common_size = head_u64[3];

    // read buffer sizes
    u64 buffer_size(0u);
    MPI_File_read_ordered(file, &buffer_size, 8, MPI_BYTE, &status);

    // read shared data if given
    if(common_size > u64(0))
    {
      // allocate shared buffer
      if(bcast_common || (comm.rank() == root_rank))
        common.resize(common_size);

      // read shared data on root rank
      if(comm.rank() == root_rank)
      {
        // write shared data
        MPI_File_read_shared(file, common.data(), int(common_size), MPI_BYTE, &status);
      }

      // broadcast shared buffer if desired
      if(bcast_common)
        comm.bcast(common.data(), common.size(), root_rank);
    }

    // allocate output buffer
    buffer.resize(buffer_size, 0);

    // read buffers via collective
    MPI_File_read_ordered(file, buffer.data(), int(buffer_size), MPI_BYTE, &status);

    // close file
    MPI_File_close(&file);
  }

  void DistFileIO::write_combined(const std::vector<char>& common, const std::vector<char>& buffer,
    const String& filename, const Dist::Comm& comm, int root_rank)
  {
    typedef std::uint64_t u64;

    XASSERTM((root_rank >= 0) && (root_rank < comm.size()), "invalid root rank");

    // set up file header
    std::vector<char> header(std::size_t(32), 0);//, padding(std::size_t(64), 0);
    u64* head_u64 = reinterpret_cast<u64*>(header.data());

    // get buffer sizes
    const u64 common_size = common.size();
    const u64 buffer_size = buffer.size();

    // set magic
    head_u64[0] = magic_combined;

    // gather combined file size at rank 0
    u64 file_size = 0u;
    comm.allreduce(&buffer_size, &file_size, std::size_t(1), Dist::op_sum);
    file_size += common_size; // shared buffer size
    file_size += u64(comm.size())*u64(8); // buffer size for each rank
    file_size += u64(header.size()); // header size
    head_u64[1] = file_size;

    // number of processes
    head_u64[2] = u64(comm.size());

    // common buffer size without padding
    head_u64[3] = common_size;

    // select file access mode
    const int amode = MPI_MODE_WRONLY | MPI_MODE_CREATE;

    // open file
    MPI_Status status;
    MPI_File file = MPI_FILE_NULL;
    MPI_File_open(comm.mpi_comm(), filename.c_str(), amode, MPI_INFO_NULL, &file);
    XASSERTM(file != MPI_FILE_NULL, "failed to open file via MPI_File_open");

    // broadcast file size
    comm.bcast(&file_size, std::size_t(1), root_rank);

    // set final file size
    MPI_File_set_size(file, MPI_Offset(file_size));

    // write header on root rank
    if(comm.rank() == root_rank)
    {
      // write header
      MPI_File_write_shared(file, header.data(), int(header.size()), MPI_BYTE, &status);
    }

    // write buffer sizes via collective
    MPI_File_write_ordered(file, &buffer_size, 8, MPI_BYTE, &status);

    // write shared data if given
    if(!common.empty() && (comm.rank() == root_rank))
    {
      // write shared data
      MPI_File_write_shared(file, common.data(), int(common_size), MPI_BYTE, &status);
    }

    // write buffers via collective
    MPI_File_write_ordered(file, buffer.data(), int(buffer_size), MPI_BYTE, &status);

    // close file
    MPI_File_close(&file);
  }

#else // non-MPI implementation

  void DistFileIO::read_common(std::stringstream& stream, const String& filename, const Dist::Comm&, int)
  {
    _read_file(stream, filename);
  }

  void DistFileIO::read_common(BinaryStream& stream, const String& filename, const Dist::Comm&, int)
  {
    _read_file(stream, filename);
  }

  void DistFileIO::read_sequence(std::stringstream& stream, const String& pattern, const Dist::Comm&)
  {
    _read_file(stream, _rankname(pattern, Index(0)));
  }

  void DistFileIO::read_sequence(BinaryStream& stream, const String& pattern, const Dist::Comm&)
  {
    _read_file(stream, _rankname(pattern, Index(0)));
  }

  void DistFileIO::write_sequence(std::stringstream& stream, const String& pattern, const Dist::Comm&, bool truncate)
  {
    _write_file(stream, _rankname(pattern, Index(0)), truncate);
  }

  void DistFileIO::write_sequence(BinaryStream& stream, const String& pattern, const Dist::Comm&, bool truncate)
  {
    _write_file(stream, _rankname(pattern, Index(0)), truncate);
  }

  void DistFileIO::read_ordered(void* buffer, const std::size_t size, const String& filename, const Dist::Comm&)
  {
    XASSERT((buffer != nullptr) || (size == std::size_t(0)));

    // determine output mode
    std::ios_base::openmode mode = std::ios_base::in|std::ios_base::binary;

    // open output file
    std::ifstream ifs(filename, mode);
    if(!ifs.is_open() || !ifs.good())
      throw FileNotFound(filename);

    // read buffer
    if(size > std::size_t(0))
    {
      ifs.read(reinterpret_cast<char*>(buffer), std::streamsize(size));
    }

    // close stream
    ifs.close();
  }

  void DistFileIO::write_ordered(const void* buffer, const std::size_t size, const String& filename, const Dist::Comm&, bool truncate)
  {
    XASSERT((buffer != nullptr) || (size == std::size_t(0)));

    // determine output mode
    std::ios_base::openmode mode = std::ios_base::out|std::ios_base::binary;
    if(truncate)
      mode |= std::ios_base::trunc;

    // open output file
    std::ofstream ofs(filename, mode);
    if(!ofs.is_open() || !ofs.good())
      throw FileNotCreated(filename);

    // write buffer
    if(size > std::size_t(0))
    {
      ofs.write(reinterpret_cast<const char*>(buffer), std::streamsize(size));
    }

    // close stream
    ofs.close();
  }

  void DistFileIO::read_combined(std::vector<char>& shared, std::vector<char>& buffer,
    const String& filename, const Dist::Comm&, int, bool)
  {
    typedef std::uint64_t u64;

    std::vector<char> header(std::size_t(40), 0);
    u64* head_u64 = reinterpret_cast<u64*>(header.data());

    // determine output mode
    std::ios_base::openmode mode = std::ios_base::in|std::ios_base::binary;

    // open output file
    std::ifstream ifs(filename, mode);
    if(!ifs.is_open() || !ifs.good())
      throw FileNotFound(filename);

    // read header
    ifs.read(header.data(), std::streamsize(header.size()));

    // check file magic
    String msg1 = String("input file is not a valid FEAT3 file: ") + filename;
    XASSERTM(head_u64[0] == magic_combined, msg1.c_str());

    // check number of processes
    String msg2 = String("invalid number of processes: serial program but file contains data for ") +
      stringify(head_u64[2]) + " processes";
    XASSERTM(head_u64[2] == u64(1), msg2.c_str());

    // get shared buffer size
    const u64 shared_size = head_u64[3];
    const u64 buffer_size = head_u64[4];

    // read shared data if given
    if(shared_size > u64(0))
    {
      // allocate and read shared buffer
      shared.resize(shared_size);
      ifs.read(shared.data(), std::streamsize(shared.size()));
    }

    // read buffer data
    if(buffer_size > u64(0))
    {
      buffer.resize(buffer_size);
      ifs.read(buffer.data(), std::streamsize(buffer.size()));
    }

    // close stream
    ifs.close();
  }

  void DistFileIO::write_combined(const std::vector<char>& shared, const std::vector<char>& buffer,
    const String& filename, const Dist::Comm&, int)
  {
    typedef std::uint64_t u64;

    // determine output mode
    std::ios_base::openmode mode = std::ios_base::out|std::ios_base::binary|std::ios_base::trunc;

    // open output file
    std::ofstream ofs(filename, mode);
    if(!ofs.is_open() || !ofs.good())
      throw FileNotCreated(filename);

    // set up file header
    std::vector<char> header(std::size_t(40), 0);
    u64* head_u64 = reinterpret_cast<u64*>(header.data());

    // get buffer sizes
    const u64 shared_size = shared.size();
    const u64 buffer_size = buffer.size();

    // set magic
    head_u64[0] = magic_combined;

    // gather combined file size at rank 0
    head_u64[1] = u64(40) + buffer_size + shared_size;

    // number of processes
    head_u64[2] = u64(1);

    // shared buffer size without padding
    head_u64[3] = shared_size;

    // actually the buffer size for rank 0
    head_u64[4] = buffer_size;

    // write header
    ofs.write(header.data(), std::streamsize(header.size()));

    // write shared data
    if(!shared.empty())
      ofs.write(shared.data(), std::streamsize(shared.size()));

    // write process data
    if(!buffer.empty())
      ofs.write(buffer.data(), std::streamsize(buffer.size()));

    // close stream
    ofs.close();
  }
#endif // FEAT_HAVE_MPI

} // namespace FEAT
