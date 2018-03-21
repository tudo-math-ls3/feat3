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

  void DistFileIO::read_common(std::stringstream& stream, const String& filename, const Dist::Comm& comm)
  {
    // rank 0 reads the file
    if(comm.rank() == 0)
    {
      _read_file(stream, filename);
    }

    // synchronise
    comm.bcast_stringstream(stream);

    // sanity check
    XASSERT(stream.good());
  }

  void DistFileIO::read_common(BinaryStream& stream, const String& filename, const Dist::Comm& comm)
  {
    // rank 0 reads the file
    if(comm.rank() == 0)
    {
      _read_file(stream, filename);
    }

    // synchronsise
    comm.bcast_binarystream(stream);

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

  void DistFileIO::write_ordered(const void* buffer, const std::size_t size, const String& filename, const Dist::Comm& comm, bool truncate)
  {
    XASSERT((buffer != nullptr) || (size == std::size_t(0)));

    // select file access mode
    int amode = MPI_MODE_WRONLY | MPI_MODE_CREATE;

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

#else // non-MPI implementation

  void DistFileIO::read_common(std::stringstream& stream, const String& filename, const Dist::Comm&)
  {
    _read_file(stream, filename);
  }

  void DistFileIO::read_common(BinaryStream& stream, const String& filename, const Dist::Comm&)
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
    ofs.write(reinterpret_cast<const char*>(buffer), std::streamsize(size));

    // close stream
    ofs.close();
  }

#endif // FEAT_HAVE_MPI

} // namespace FEAT
