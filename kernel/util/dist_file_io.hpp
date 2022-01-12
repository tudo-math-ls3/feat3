// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

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
   * There are three basic sets of functions for reading files as well as a combined version:
   * - #read_common: one process reads a single text or binary file and broadcasts its
   *   contents to all other processes
   * - #read_sequence: each process reads a single text or binary file, which is indexed
   *   by the process rank
   * - #read_ordered: all processes read from a single common binary file, ordered by ranks,
   *   and using the official MPI I/O routines
   * - #read_combined: all processes read from a single common binary file which contains
   *   a common data block, which is read by all processes, as well as a private data block,
   *   that is read by each process individually using a combination of the official
   *   MPI I/O routines. See the corresponding remarks section below for details.
   *
   * Moreover, there are two sets of functions for writing files:
   * - #write_sequence: each process writes a single text or binary file, which is indexed
   *   by the process rank
   * - #write_ordered: all processes write into a single common binary file, ordered by ranks,
   *   and using the official MPI I/O routines
   * - #write_combined: all processes write into a single common binary file, which contains
   *   a common data block as well as a private data block, that is written by each process
   *   ordered by ranks using a combination of the official MPI I/O routines.
   *   See the corresponding remarks section below for details.
   *
   * Each common/sequence function comes in two overloads: one for objects of type std::stringstream
   * for text files  and another one for objects of type BinaryStream for binary files.
   *
   * <b>Remarks regarding Sequence I/O functions:</b>\n
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
   * <b>Remarks regarding Combined I/O functions:</b>\n
   * In contrast to the  #write_ordered() functions also offered by this class, the
   * #write_combined() function does not simply write out the raw data arrays to the file,
   * but it also writes out a simple file header including all the required array size
   * informations for each data array. In consequence, when reading in a file using the
   * #read_combined() function, it is not necessary to allocate the output buffer for each
   * process beforehand (as opposed to the #read_ordered() function) because the required
   * buffer size is automatically read from file.
   *
   * In addition to that, the combined data file can also contain a \e common data array,
   * which is automatically broadcasted to all processes by the #read_combined() function
   * upon reading. This common data array can be used to store information that is
   * identical for all processes and therefore does not need to be stored by each process
   * individually to reduce redundancy within the file and therefore reduce the total
   * file size. When calling the #write_combined() function, the caller can decide which
   * process (called "root") is responsible for writing out the common data array.
   *
   * The binary output file is composed of four sections or chunks:
     \verbatim
       +---+---+---+---+
       | H | S | C | P |
       +---+---+---+---+
     \endverbatim
   * - 'H': a 32-byte file header which consists of four \c uint64 entries:
   *   - the file magic, which must be equal to #magic_combined = <c>"FEAT3CDF"</c>
   *   - the total filesize in bytes
   *   - the number of processes that were used to write the file
   *   - the size of the common buffer in bytes
   * - 'S': an \c uint64 array which contains the sizes of the process data arrays in bytes
   * - 'C': the raw common data array, may be empty
   * - 'P': the data arrays from each process, ordered by ranks
   *
   * Assume we have (n+1) processes and each process has an identical copy of a common data
   * array \c C of size <c>SC := sizeof(C)</c> bytes as well as a private data array <c>Pk</c>
   * of size <c>Sk := sizeof(Pk)</c> bytes:
     \verbatim
      +--------+   +--------+   +--------+       +--------+
      | Rank 0 |   | Rank 1 |   | Rank 2 |       | Rank n |
      | C | P0 |   | C | P1 |   | C | P2 |  ...  | C | Pn |
      +--------+   +--------+   +--------+       +--------+
       \_/ \__/     \_/ \__/     \_/ \__/         \_/ \__/
       SC   S0      SC   S1      SC   S2          SC   Sn
        bytes        bytes        bytes            bytes
     \endverbatim
   * The corresponding binary output file created by the #write_combined() function will have
   * a total fizesize of <c>32 + 8*(n+1) + SC + S0 + S1 + S2 + ... + Sn</c> bytes:
     \verbatim
      +--------+-----------------+----------+----------+----------+-----+----------+
      | Header | S0 S1 S2 ... Sn |    C     |    P0    |    P1    | ... |    Pn    |
      +--------+-----------------+----------+----------+----------+-----+----------+
       \______/ \_______________/ \________/ \________/ \________/       \________/
       32 bytes   8*(n+1) bytes    SC bytes   S0 bytes   S1 bytes   ...   Sn bytes
     \endverbatim
   *
   *
   * \author Peter Zajac
   */
  class DistFileIO
  {
  public:
    /**
     * \brief Magic number for combined distributed files as used by #read_combined() and #write_combined()
     *
     * This magic number represents the ASCII-string <c>"FEAT3CDF"</c> encoded as a little endian uint64.
     */
    static constexpr std::uint64_t magic_combined = 0x4644433354414546ull;

    /**
     * \brief Reads a common text file for all ranks.
     *
     * This function reads a single file on the root rank and broadcasts its contents to all other ranks.
     *
     * \param[out] stream
     * A \transient reference to the string stream that receives the file contents.
     *
     * \param[in] filename
     * The name of the file to be read.
     *
     * \param[in] comm
     * The communicator to be used for synchronization. Ignored if compiled without MPI.
     *
     * \param[in] root_rank
     * Specifies which rank should actually read the file. Ignored if compiled without MPI.
     *
     * \throws FileNotFound if the file could not be opened.
     */
    static void read_common(std::stringstream& stream, const String& filename, const Dist::Comm& comm, int root_rank = 0);

    /** \copydoc #read_common() */
    static void read_common(std::stringstream& stream, const String& filename)
    {
      Dist::Comm comm(Dist::Comm::world());
      read_common(stream, filename, comm);
    }

    /**
     * \brief Reads a common binary file for all ranks.
     *
     * This function reads a single file on the root rank and broadcasts its contents to all other ranks.
     *
     * \param[out] stream
     * A \transient reference to the binary stream that receives the file contents.
     *
     * \param[in] filename
     * The name of the file to be read.
     *
     * \param[in] comm
     * The communicator to be used for synchronization. Ignored if compiled without MPI.
     *
     * \param[in] root_rank
     * Specifies which rank should actually read the file. Ignored if compiled without MPI.
     *
     * \throws FileNotFound if the file could not be opened.
     */
    static void read_common(BinaryStream& stream, const String& filename, const Dist::Comm& comm, int root_rank = 0);

    /** \copydoc #read_common() */
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
     * A \transient reference to the string stream that receives the file contents.
     *
     * \param[in] pattern
     * The pattern that is to be used for filename generation.
     *
     * \param[in] comm
     * The communicator to be used for synchronization. Ignored if compiled without MPI.
     *
     * \throws FileNotFound if the file could not be opened.
     */
    static void read_sequence(std::stringstream& stream, const String& pattern, const Dist::Comm& comm);

    /** \copydoc #read_sequence() */
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
     * A \transient reference to the binary stream that receives the file contents.
     *
     * \param[in] pattern
     * The pattern that is to be used for filename generation.
     *
     * \param[in] comm
     * The communicator to be used for synchronization. Ignored if compiled without MPI.
     *
     * \throws FileNotFound if the file could not be opened.
     */
    static void read_sequence(BinaryStream& stream, const String& pattern, const Dist::Comm& comm);

    /** \copydoc #read_sequence() */
    static void read_sequence(BinaryStream& stream, const String& pattern)
    {
      Dist::Comm comm(Dist::Comm::world());
      read_sequence(stream, pattern, comm);
    }

    /**
     * \brief Writes a rank-indexed text file sequence.
     *
     * \param[out] stream
     * A \transient reference to the string stream whose contents are to be written to the file.
     *
     * \param[in] pattern
     * The pattern that is to be used for filename generation.
     *
     * \param[in] comm
     * The communicator to be used for synchronization. Ignored if compiled without MPI.
     *
     * \param[in] truncate
     * Specifies whether the output file(s) are to be truncated to the output size.
     *
     * \throws FileNotCreated if the file could not be opened.
     */
    static void write_sequence(std::stringstream& stream, const String& pattern, const Dist::Comm& comm, bool truncate = true);

    /** \copydoc #write_sequence() */
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
     * A \transient reference to the binary stream whose contents are to be written to the file.
     *
     * \param[in] pattern
     * The pattern that is to be used for filename generation.
     *
     * \param[in] comm
     * The communicator to be used for synchronization. Ignored if compiled without MPI.
     *
     * \param[in] truncate
     * Specifies whether the output file(s) are to be truncated to the output size.
     *
     * \throws FileNotCreated if the file could not be opened.
     */
    static void write_sequence(BinaryStream& stream, const String& pattern, const Dist::Comm& comm, bool truncate = true);

    /** \copydoc #write_sequence() */
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
     * A \transient pointer to the binary buffer that is to read into. Must not be \c nullptr.
     *
     * \param[in] size
     * The size of this process's buffer in bytes. May differ on each process.
     *
     * \param[in] filename
     * The name of the common input file. Must be the same on all calling processes.
     *
     * \param[in] comm
     * The communicator to be used for synchronization. Ignored if compiled without MPI.
     */
    static void read_ordered(void* buffer, const std::size_t size, const String& filename, const Dist::Comm& comm);

    /** \copydoc #read_ordered() */
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
     * A \transient pointer to the binary buffer that is to be written. Must not be \c nullptr.
     *
     * \param[in] size
     * The size of this process's buffer in bytes. May differ on each process.
     *
     * \param[in] filename
     * The name of the common output file. Must be the same on all calling processes.
     *
     * \param[in] comm
     * The communicator to be used for synchronization. Ignored if compiled without MPI.
     *
     * \param[in] truncate
     * Specifies whether the output file(s) are to be truncated to the output size.
     */
    static void write_ordered(const void* buffer, const std::size_t size, const String& filename, const Dist::Comm& comm, bool truncate = true);

    /** \copydoc #write_ordered() */
    static void write_ordered(const void* buffer, const std::size_t size, const String& filename, bool truncate = true)
    {
      Dist::Comm comm(Dist::Comm::world());
      write_ordered(buffer, size, filename, comm, truncate);
    }

    /**
     * \brief Writes a binary stream into a common binary file in rank order.
     *
     * This function write a single common binary file, where the individual processes
     * write their outputs in ascending order by their rank.
     *
     * \note This function is effectively a wrapper around \b MPI_File_write_ordered.
     *
     * \param[in] stream
     * A \transient reference to the binary stream whose contents are to be written to the file.
     *
     * \param[in] filename
     * The name of the common output file. Must be the same on all calling processes.
     *
     * \param[in] comm
     * The communicator to be used for synchronization. Ignored if compiled without MPI.
     *
     * \param[in] truncate
     * Specifies whether the output file(s) are to be truncated to the output size.
     */
    static void write_ordered(BinaryStream& stream, const String& filename, const Dist::Comm& comm, bool truncate = true)
    {
      write_ordered(stream.data(), std::size_t(stream.size()), filename, comm, truncate);
    }

    /** \copydoc #write_ordered() */
    static void write_ordered(BinaryStream& stream, const String& filename, bool truncate = true)
    {
      Dist::Comm comm(Dist::Comm::world());
      write_ordered(stream.data(), std::size_t(stream.size()), filename, comm, truncate);
    }

    /**
     * \brief Reads a combined shared/ordered binary file as written by the write_combined function.
     *
     * This function reads a single binary file, which contains both a common buffer, that
     * is shared amongst all processes, as well one private buffer for each individual process.
     *
     * \param[out] common
     * A byte-vector that receives the common buffer whose contents are to be read from the file.
     * This vector is allocated to its correct size by this function. Depending on the value
     * of the parameter \p bcast_common, this common buffer in only filled on the root process
     * (<c>bcast_common = false</c>) or it is filled (broadcasted) onto all processes
     * (<c>bcast_common = true</c>).
     *
     * \param[out] buffer
     * A byte-vector that receives the individual buffer which is read from file by each process.
     * This vector is allocated to its correct size on each process by this function.
     *
     * \param[in] filename
     * The name of the common output file. Must be the same on all calling processes.
     *
     * \param[in] comm
     * The communicator to be used for synchronization. Ignored if compiled without MPI.
     *
     * \param[in] root_rank
     * Specifies which process should actually read the file. Ignored if compiled without MPI.
     *
     * \param[in] bcast_common
     * Specifies whether the shared buffer is to be broadcasted to all processes or whether
     * it is only filled on the root process. Ignored if compiled without MPI.
     */
    static void read_combined(std::vector<char>& common, std::vector<char>& buffer,
      const String& filename, const Dist::Comm& comm, int root_rank = 0, bool bcast_common = true);

    /** \copydoc #read_combined() */
    static void read_combined(BinaryStream& common, BinaryStream& buffer,
      const String& filename, const Dist::Comm& comm, int root_rank = 0, bool bcast_common = true)
    {
      read_combined(common.container(), buffer.container(), filename, comm, root_rank, bcast_common);
    }

    /**
     * \brief Writes a combined shared/ordered binary file.
     *
     * This function writes a single binary file, which contains both a common buffer, that
     * is shared amongst all processes, as well one private buffer for each individual process.
     *
     * \param[in] common
     * The common buffer whose contents are to be written to the file.
     * This buffer is only required on the root process and is ignored on all other processes.
     *
     * \param[in] buffer
     * The individual buffer for each process whose contents are to be written to the file.
     * This buffer is required on each process and all process buffers will be written
     * to the file in rank-ascending order.
     *
     * \param[in] filename
     * The name of the common output file. Must be the same on all calling processes.
     *
     * \param[in] comm
     * The communicator to be used for synchronization. Ignored if compiled without MPI.
     *
     * \param[in] root_rank
     * Specifies which rank should actually write the file. Ignored if compiled without MPI.
     */
    static void write_combined(const std::vector<char>& common, const std::vector<char>& buffer,
      const String& filename, const Dist::Comm& comm, int root_rank = 0);

    /** \copydoc #write_combined() */
    static void write_combined(const BinaryStream& common, const BinaryStream& buffer,
      const String& filename, const Dist::Comm& comm, int root_rank = 0)
    {
      write_combined(common.container(), buffer.container(), filename, comm, root_rank);
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
