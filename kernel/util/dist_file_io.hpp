#pragma once
#ifndef KERNEL_UTIL_DIST_FILE_IO_HPP
#define KERNEL_UTIL_DIST_FILE_IO_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/binary_stream.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/dist.hpp>

#include <sstream>

namespace FEAT
{
  /**
   * \brief Distributed File Input/Output class
   *
   * This class implements various static functions for reading and writing files in
   * text or binary mode in a distributed (i.e. parallel) build.
   *
   * There are three sets of functions for reading files:
   * - #read_common: one process reads a single text or binary file and broadcasts its
   *   contents to all other processes
   * - #read_sequence: each process reads a single text or binary file, which is indexed
   *   by the process rank
   * - #read_ordered: all processes read from a single common binary file, ordered by ranks,
   *   and using the official MPI I/O routines
   *
   * Moreover, there are two sets of functions for writing files:
   * - #write_sequence: each process writes a single text or binary file, which is indexed
   *   by the process rank
   * - #write_ordered: all processes write into a single common binary file, ordered by ranks,
   *   and using the official MPI I/O routines
   *
   * Each function comes in two overloads: one for objects of type std::stringstream for text files
   * and another one for objects of type BinaryStream for binary files.
   *
   * For both the read_sequence and write_sequence functions, one specifies a filename pattern
   * which is used as a template for the generation of the filename for each rank.
   * The filename pattern has to contain a continuous block of asterisks (*), which serves as
   * a place-holder for the 0-padded rank number.
   *
   * <b>Example:</b>\n
   * If pattern is <c>./my_file_****.txt</c>, then
   * - rank 0 reads/writes the file <c>./my_file_0000.txt</c>
   * - rank 1 reads/writes the file <c>./my_file_0001.txt</c>
   * - rank 2 reads/writes the file <c>./my_file_0002.txt</c>
   * - etc.
   *
   * \author Peter Zajac
   */
  class DistFileIO
  {
  public:
    /**
     * \brief Reads a common text file for all ranks.
     *
     * This function reads a single file on rank 0 and broadcasts its contents to all other ranks.
     *
     * \param[out] stream
     * The string stream that receives the file contents.
     *
     * \param[in] filename
     * The name of the file to be read.
     *
     * \param[in] comm
     * The communicator to be used for synchronisation. Ignored if compiled without MPI.
     *
     * \throws FileNotFound if the file could not be opened.
     */
    static void read_common(std::stringstream& stream, const String& filename, const Dist::Comm& comm);

    static void read_common(std::stringstream& stream, const String& filename)
    {
      Dist::Comm comm(Dist::Comm::world());
      read_common(stream, filename, comm);
    }

    /**
     * \brief Reads a common binary file for all ranks.
     *
     * This function reads a single file on rank 0 and broadcasts its contents to all other ranks.
     *
     * \param[out] stream
     * The binary stream that receives the file contents.
     *
     * \param[in] filename
     * The name of the file to be read.
     *
     * \param[in] comm
     * The communicator to be used for synchronisation. Ignored if compiled without MPI.
     *
     * \throws FileNotFound if the file could not be opened.
     */
    static void read_common(BinaryStream& stream, const String& filename, const Dist::Comm& comm);

    static void read_common(BinaryStream& stream, const String& filename)
    {
      Dist::Comm comm(Dist::Comm::world());
      read_common(stream, filename, comm);
    }

    /**
     * \brief Reads a rank-indexed text file sequence.
     *
     * This function reads a sequence of files, where each rank reads a single indexed file.
     *
     * \param[out] stream
     * The string stream that receives the file contents.
     *
     * \param[in] pattern
     * The pattern that is to be used for filename generation.
     *
     * \param[in] comm
     * The communicator to be used for synchronisation. Ignored if compiled without MPI.
     *
     * \throws FileNotFound if the file could not be opened.
     */
    static void read_sequence(std::stringstream& stream, const String& pattern, const Dist::Comm& comm);

    static void read_sequence(std::stringstream& stream, const String& pattern)
    {
      Dist::Comm comm(Dist::Comm::world());
      read_sequence(stream, pattern, comm);
    }


    /**
     * \brief Reads a rank-indexed binary file sequence.
     *
     * This function reads a sequence of files, where each rank reads a single indexed file.
     *
     * \param[out] stream
     * The binary stream that receives the file contents.
     *
     * \param[in] pattern
     * The pattern that is to be used for filename generation.
     *
     * \param[in] comm
     * The communicator to be used for synchronisation. Ignored if compiled without MPI.
     *
     * \throws FileNotFound if the file could not be opened.
     */
    static void read_sequence(BinaryStream& stream, const String& pattern, const Dist::Comm& comm);

    static void read_sequence(BinaryStream& stream, const String& pattern)
    {
      Dist::Comm comm(Dist::Comm::world());
      read_sequence(stream, pattern, comm);
    }

    /**
     * \brief Writes a rank-indexed text file sequence.
     *
     * \param[out] stream
     * The string stream whose contents are to be written to the file.
     *
     * \param[in] pattern
     * The pattern that is to be used for filename generation.
     *
     * \param[in] comm
     * The communicator to be used for synchronisation. Ignored if compiled without MPI.
     *
     * \param[in] truncate
     * Specifies whether the output file(s) are to be truncated to the output size.
     *
     * \throws FileNotCreated if the file could not be opened.
     */
    static void write_sequence(std::stringstream& stream, const String& pattern, const Dist::Comm& comm, bool truncate = true);

    static void write_sequence(std::stringstream& stream, const String& pattern, bool truncate = true)
    {
      Dist::Comm comm(Dist::Comm::world());
      write_sequence(stream, pattern, comm, truncate);
    }

    /**
     * \brief Writes a rank-indexed binary file sequence.
     *
     * This function writes a sequence of files, where each rank writes a single indexed file.
     *
     * \param[out] stream
     * The binary stream whose contents are to be written to the file.
     *
     * \param[in] pattern
     * The pattern that is to be used for filename generation.
     *
     * \param[in] comm
     * The communicator to be used for synchronisation. Ignored if compiled without MPI.
     *
     * \param[in] truncate
     * Specifies whether the output file(s) are to be truncated to the output size.
     *
     * \throws FileNotCreated if the file could not be opened.
     */
    static void write_sequence(BinaryStream& stream, const String& pattern, const Dist::Comm& comm, bool truncate = true);

    static void write_sequence(BinaryStream& stream, const String& pattern, bool truncate = true)
    {
      Dist::Comm comm(Dist::Comm::world());
      write_sequence(stream, pattern, comm, truncate);
    }

    /**
     * \brief Reads a buffer from a common binary file in rank order.
     *
     * \note This function is effectively a wrapper around \b MPI_File_read_ordered.
     *
     * \param[out] buffer
     * A pointer to the binary buffer that is to read into. Must not be \c nullptr.
     *
     * \param[in] size
     * The size of this process's buffer in bytes. May differ on each process.
     *
     * \param[in] filename
     * The name of the common input file. Must be the same on all calling processes.
     *
     * \param[in] comm
     * The communicator to be used for synchronisation. Ignored if compiled without MPI.
     */
    static void read_ordered(void* buffer, const std::size_t size, const String& filename, const Dist::Comm& comm);

    static void read_ordered(void* buffer, const std::size_t size, const String& filename)
    {
      Dist::Comm comm(Dist::Comm::world());
      read_ordered(buffer, size, filename, comm);
    }

    /**
     * \brief Writes a buffer into a common binary file in rank order.
     *
     * This function write a single common binary file, where the individual processes
     *  write their outputs in ascending order by their rank.
     *
     * \note This function is effectively a wrapper around \b MPI_File_write_ordered.
     *
     * \param[in] buffer
     * A pointer to the binary buffer that is to be written. Must not be \c nullptr.
     *
     * \param[in] size
     * The size of this process's buffer in bytes. May differ on each process.
     *
     * \param[in] filename
     * The name of the common output file. Must be the same on all calling processes.
     *
     * \param[in] comm
     * The communicator to be used for synchronisation. Ignored if compiled without MPI.
     *
     * \param[in] truncate
     * Specifies whether the output file(s) are to be truncated to the output size.
     */
    static void write_ordered(const void* buffer, const std::size_t size, const String& filename, const Dist::Comm& comm, bool truncate = true);

    static void write_ordered(const void* buffer, const std::size_t size, const String& filename, bool truncate = true)
    {
      Dist::Comm comm(Dist::Comm::world());
      write_ordered(buffer, size, filename, comm, truncate);
    }

    /**
     * \brief Writes a binary stream into a common binary file in rank order.
     *
     * This function write a single common binary file, where the individual processes
     *  write their outputs in ascending order by their rank.
     *
     * \note This function is effectively a wrapper around \b MPI_File_write_ordered.
     *
     * \param[in] stream
     * The binary stream whose contents are to be written to the file.
     *
     * \param[in] filename
     * The name of the common output file. Must be the same on all calling processes.
     *
     * \param[in] comm
     * The communicator to be used for synchronisation. Ignored if compiled without MPI.
     *
     * \param[in] truncate
     * Specifies whether the output file(s) are to be truncated to the output size.
     */
    static void write_ordered(BinaryStream& stream, const String& filename, const Dist::Comm& comm, bool truncate = true)
    {
      write_ordered(stream.data(), std::size_t(stream.size()), filename, comm, truncate);
    }

    static void write_ordered(BinaryStream& stream, const String& filename, bool truncate = true)
    {
      Dist::Comm comm(Dist::Comm::world());
      write_ordered(stream.data(), std::size_t(stream.size()), filename, comm, truncate);
    }

  protected:
    /// auxiliary function: build a rank filename from a pattern
    static String _rankname(const String& pattern, int rank);
    /// auxiliary function: read a file into a string stream
    static void _read_file(std::stringstream& stream, const String& filename);
    /// auxiliary function: read a file into a binary stream
    static void _read_file(BinaryStream& stream, const String& filename);
    /// auxiliary function: write a string stream to a file
    static void _write_file(std::stringstream& stream, const String& filename, bool truncate);
    /// auxiliary function: write a binary stream to a file
    static void _write_file(BinaryStream& stream, const String& filename, bool truncate);
  }; // class DistFileIO
} // namespace FEAT

#endif // KERNEL_UTIL_DIST_FILE_IO_HPP
