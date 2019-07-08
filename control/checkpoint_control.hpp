// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef CONTROL_CHECKPOINT_HPP
#define CONTROL_CHECKPOINT_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/dist_file_io.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/binary_stream.hpp>
#include <kernel/util/pack.hpp>

#include <cstdio>
#include <cstring>
#include <map>
#include <sys/stat.h>
#include <memory>

namespace FEAT
{
  namespace Control
  {
    /**
     * \brief Checkpoint interface class.
     *
     * This abstract class provides the necessary interface for all checkpointable objects.
     *
     */
    class Checkpointable
    {
    public:
      /**
         * \brief Calculate size
         *
         * Calculate size of complete object as it is stored in the checkpoint
         *
         * \return size of object
         *
         */
      virtual uint64_t get_checkpoint_size() = 0;

      /**
         * \brief Extract object from checkpoint
         *
         * Restores complete object with all its contents from checkpoint data
         *
         * \param[out] data object as bytestrem
         *
         */
      virtual void restore_from_checkpoint_data(std::vector<char> & data) = 0;

      /**
         * \brief Collect checkpoint data from object
         *
         * Adds the condensed complete object with all its contents to the end of the checkpoint buffer
         *
         * \param[out] data object as bytestream
         *
         */
      virtual void set_checkpoint_data(std::vector<char> & data) = 0;

    protected:
      virtual ~Checkpointable() {}
    }; // class Checkpointable

    /**
     * \brief Wrapper class for checkpointable objects
     *
     * Any class that implements the Checkpointable interface can be wrapped by CheckpointableWrapper
     * and then be used with CheckpointControl.
     */
    template <typename Object_>
    class CheckpointableWrapper : public Checkpointable
    {
    private:
      Object_ & _object;

    public:
      /// Contructor
      explicit CheckpointableWrapper(Object_ & obj) :
        _object(obj) {}

      /**
         * \brief Calculate size
         *
         * Calculate size of complete object as it is stored in the checkpoint
         *
         * \return size of object
         *
         */
      virtual uint64_t get_checkpoint_size() override
      {
        return _object.get_checkpoint_size();
      }

      /**
         * \brief Extract object from checkpoint
         *
         * Restores complete object with all its contents from checkpoint data
         *
         * \param[out] data object as bytestrem
         *
         */
      virtual void restore_from_checkpoint_data(std::vector<char> & data) override
      {
        _object.restore_from_checkpoint_data(data);
      }

      /**
         * \brief Collect checkpoint data from object
         *
         * Adds the condensed complete object with all its contents to the end of the checkpoint buffer
         *
         * \param[out] data object as bytestream
         *
         */
      virtual void set_checkpoint_data(std::vector<char> & data) override
      {
        _object.set_checkpoint_data(data);
      }
    }; // class CheckpointableWrapper

    /**
     * \brief Checkpoint/Restart infrastructure
     *
     * This class provides a general interface for checkpoint/restart usecases.
     * Compatible objects can be summed up to a single checkpoint.
     * This checkpoint can then be stored and restored via mpi io mechanics.
     * Afterwards, any object can be restored from the recreated checkpoint.
     *
     * \note The checkpoint control class object is used on a per rank basis. The contents of all
     * ranks objects will be written into one single file.
     * Later on, the loaded objects will, again, be local objects in every rank with distinct contents.
     *
     * \todo do we need extra handling for machines with no global acessible file servers?
     */
    class CheckpointControl
    {
    private:
      /// Mapping of identifier string to pointer to the checkpointable object
      std::map<String, std::shared_ptr<Checkpointable>> _checkpointable_by_identifier;
      /// Pointer to the char array read from the input file during restore
      char * _input_array;
      /// Mapping of identifier string to the offset in the input file
      std::map<String, uint64_t> _offset_by_identifier;
      /// the mpi communicator identifying our mpi context
      const Dist::Comm & _comm;

      /**
         * \brief Build checkpoint buffer
         *
         * Collect all data needed for a checkpoint from every object registered at the checkpoint control.
         *
         * \param[in] buffer buffer containing the collected data
         * \returns the size of the buffer
         *
         */
      uint64_t _collect_checkpoint_data(std::vector<char> & buffer)
      {
        uint64_t checkpoint_size(0);
        std::map<String, std::tuple<uint64_t, uint64_t>> sizes;

        for (auto const & it : _checkpointable_by_identifier)
        {
          uint64_t identifierlength = (uint64_t)it.first.length();
          uint64_t datalength = it.second->get_checkpoint_size();

          checkpoint_size += identifierlength + datalength + sizeof(uint64_t) + sizeof(uint64_t);
          sizes[it.first] = std::make_tuple(identifierlength, datalength);
        }

        buffer.reserve(checkpoint_size);

        for (auto const & it : _checkpointable_by_identifier)
        {
          char * cidentifierlength = reinterpret_cast<char *>(&std::get<0>(sizes[it.first]));
          buffer.insert(std::end(buffer), cidentifierlength, cidentifierlength + sizeof(uint64_t));
          buffer.insert(std::end(buffer), it.first.begin(), it.first.end());
          char * cdatalength = reinterpret_cast<char *>(&std::get<1>(sizes[it.first]));
          buffer.insert(std::end(buffer), cdatalength, cdatalength + sizeof(uint64_t));
          it.second->set_checkpoint_data(buffer);
        }

        return checkpoint_size;
      }

      /**
         * \brief Extract input buffer
         *
         * Extract the identifier string and offset for its data in the input array for every object
         * in the read checkpoint, and store this information for restoration in this checkpoint control object.
         *
         * \param[in] size current size of _input_array after reading the file / stream
         *
         */
      void _restore_checkpoint_data(const uint64_t size)
      {
        uint64_t stringsize;
        uint64_t datasize(0);
        size_t i = 0;

        // Loop over all objects stored in the checkpoint file
        while (i < size)
        {
          // Get size of identifier string
          ::memcpy(&stringsize, _input_array + i, sizeof(uint64_t));
          i += sizeof(uint64_t);

          // Read the identifier string and put it as key for the calculated offset to the map
          _offset_by_identifier[String(_input_array + i, stringsize)] = i + stringsize;
          i += stringsize;

          // Get the size of the Data holding the array, to get the start point of the next checkpointable object in the array
          ::memcpy(&datasize, _input_array + i, sizeof(uint64_t));
          i += sizeof(uint64_t) + datasize;
        }
      }

#ifdef FEAT_HAVE_ZLIB
      /**
         * \brief Write checkpoint files to disk
         *
         * Write two binary files:
         * [name].zcp: a compressed file holding the current status of all checkpointable objects added to the checkpoint control
         * [name].szcp: a file holding the size of the checkpoint data per rank, needed for loading the checkpoint
         *
         * \param[in] name String holding the name of the file (without the extension .zcp)
         */
      void _save_compressed(const String name)
      {
        std::vector<char> buffer;
        uint64_t sizes[2];
        sizes[0] = _collect_checkpoint_data(buffer);

        uLongf compressedLen = compressBound(static_cast<uLong>(sizes[0]));
        Bytef * compressed = new Bytef[compressedLen];
        Bytef * buffer_data = (Bytef *)buffer.data();

        int c_status = ::compress(compressed, &compressedLen, buffer_data, static_cast<uLongf>(sizes[0]));
        XASSERTM(c_status == Z_OK, "compression of checkpoint data failed");

        sizes[1] = static_cast<uint64_t>(compressedLen);

        DistFileIO::write_ordered(compressed, sizes[1], name + ".zcp", _comm);
        DistFileIO::write_ordered(sizes, 2 * sizeof(uint64_t), name + ".szcp", _comm);
      }

      /**
         * \brief Load checkpoint from disk
         *
         * Read the [name].szcp file, to restore the size of the checkpoint data per rank.
         * Read the compressed binary file [name].zcp holding a safed state of all checkpointable objects and store it to _input_array.
         *
         * Fill the _offset_by_identifier map with the data from the array.
         *
         * \param[in] name String holding the name of the file (without the extension .zcp)
         *
         * \warning Reading another input file / stream will overwrite values in _input_array and _offset_by_identifier
         */
      void _load_compressed(const String name)
      {
        struct stat stat_buf;
        stat((name + ".szcp").c_str(), &stat_buf);
        size_t filesize = (unsigned)stat_buf.st_size;
        size_t original_rank_count = (filesize / (2 * sizeof(uint64_t)));
        int world_size(_comm.size());
        XASSERTM(original_rank_count == unsigned(world_size), "number of ranks of checkpoint file and running program does not match");

        // read the file with the sizes first, so that every rank knows its checkpoint size (needed for read_ordered)
        uint64_t size_buffer[2];
        DistFileIO::read_ordered((char *)size_buffer, 2 * sizeof(uint64_t), name + ".szcp", _comm);
        uLongf size = static_cast<uLongf>(size_buffer[0]);
        uint64_t buffersize = size_buffer[1];

        auto buffer = new char[buffersize];
        DistFileIO::read_ordered(buffer, buffersize, name + ".zcp", _comm);

        _input_array = new char[size];

        Bytef * buffer_data = (Bytef *)buffer;
        int d_status = ::uncompress((Bytef *)_input_array, &size, buffer_data, static_cast<uLongf>(buffersize));
        XASSERTM(d_status == Z_OK, "decompression of checkpoint data failed");

        delete[] buffer;

        _restore_checkpoint_data(static_cast<uint64_t>(size));
      }
#endif //FEAT_HAVE_ZLIB

      /**
         * \brief Write checkpoint files to disk
         *
         * Write two binary files:
         * [name].cp: a uncompressed binary file holding the current status of all checkpointable objects added to the checkpoint control
         * [name].scp: a file holding the size of the checkpoint data per rank, needed for loading the checkpoint
         *
         * \param[in] name String holding the name of the file (without the extension .cp)
         */
      void _save(const String name)
      {
        std::vector<char> buffer;
        auto checkpoint_size = _collect_checkpoint_data(buffer);

        DistFileIO::write_ordered(buffer.data(), checkpoint_size, name + ".cp", _comm);
        DistFileIO::write_ordered(reinterpret_cast<char *>(&checkpoint_size), sizeof(uint64_t), name + ".scp", _comm);
      }

      /**
         * \brief Load checkpoint from disk
         *
         * Read the [name].scp file, to restore the size of the checkpoint data per rank.
         * Read binary file [name].cp holding a safed state of all checkpointable objects and store it to _input_array.
         *
         * Fill the _offset_by_identifier map with the data from the array.
         *
         * \param[in] name String holding the name of the file (without the extension .zcp)
         *
         * \warning Reading another input file / stream will overwrite values in _input_array and _offset_by_identifier
         */
      void _load(const String name)
      {
        struct stat stat_buf;
        stat((name + ".scp").c_str(), &stat_buf);
        size_t filesize = (unsigned)stat_buf.st_size;
        int world_size(_comm.size());
        XASSERTM((filesize / sizeof(uint64_t)) == unsigned(world_size), "number of ranks of checkpoint file and running program does not match");

        // read the file with the sizes first, so that every rank knows its checkpoint size (needed for read_ordered)
        uint64_t size_buffer[1];
        DistFileIO::read_ordered((char *)size_buffer, sizeof(uint64_t), name + ".scp", _comm);
        const uint64_t size = *size_buffer;

        // read the checkpoint file
        _input_array = new char[size];
        DistFileIO::read_ordered(_input_array, size, name + ".cp", _comm);

        _restore_checkpoint_data(size);
      }

    public:
      /**
         * \brief Constructor
         *
         * Initalise the input array as NULL pointer.
         *
         * \param[in] comm The communicator common to all stored objects
         */
      CheckpointControl(const Dist::Comm & comm) :
        _comm(comm)
      {
        _input_array = nullptr;
      }

      /**
         * \brief Destructor
         *
         * Destroy the checkpoint control instance and delete the input array.
         */
      ~CheckpointControl()
      {
        delete[] _input_array;
      }

      /**
         * \brief Delete all read input
         *
         * Clear the input read before from this checkpoint control object.
         * Does not touch the map holding the pointers to the checkpointable objects included into the following written checkpoints.
         */
      void clear_input()
      {
        delete[] _input_array;
        _offset_by_identifier.clear();
        _input_array = nullptr;
      }

      /**
         * \brief Register checkpointable object
         *
         * Adds a pointer to a checkpointable object to the checkpoint control.
         *
         * The object will be internally identified by its indentifier string (especially during restoration).
         *
         * \param[in] identifier String used as identifier
         * \param[in] object Pointer to the checkpointable object
         *
         */
      template <typename OT_>
      void add_object(String identifier, OT_ & object)
      {
        XASSERTM(_checkpointable_by_identifier.count(identifier) == 0, "a object with the given identifier is already registered; choose a unique identifier");
        auto checkpointable_object = std::make_shared<CheckpointableWrapper<OT_>>(object);
        _checkpointable_by_identifier[identifier] = checkpointable_object;
      }

      /**
         * \brief Remove checkpointable object
         *
         * Remove pointer to a checkpointable object from the checkpoint control.
         *
         * \param[in] identifier Identifier of the object to remove
         *
         */
      void remove_object(String identifier)
      {
        XASSERTM(_checkpointable_by_identifier.count(identifier) == 1, "a object with the given identifier is not registered at this checkpoint instance");
        _checkpointable_by_identifier.erase(identifier);
      }

      /**
         * \brief Restore checkpointable object
         *
         * Restore a checkpointable object from the checkpoint control.
         * The objects contents will be found by its identifier string and restored.
         *
         * \param[in] identifier Identifier of the object to restore
         * \param[out] object Reference to object
         * \param[in] add_to_checkpoint_control If true the restored object is added to the checkpoint control, default true
         *
         * \warning The user has to load an input file before and has to ensure proper size/data/whatsover type layout.
         */
      template <typename OT_>
      void restore_object(String identifier, OT_ & object, bool add_to_checkpoint_control = true)
      {
        XASSERTM(_input_array != nullptr, "no input file has been loaded before");
        XASSERTM(_offset_by_identifier.count(identifier) == 1, "input file doesn't contain data for the given identifier");

        auto checkpointable_object = std::make_shared<CheckpointableWrapper<OT_>>(object);

        uint64_t size(0);
        ::memcpy(&size, _input_array + _offset_by_identifier[identifier], sizeof(uint64_t));
        std::vector<char> buffer(&_input_array[_offset_by_identifier[identifier] + sizeof(uint64_t)], &_input_array[_offset_by_identifier[identifier] + sizeof(uint64_t) + size]);
        checkpointable_object->restore_from_checkpoint_data(buffer);

        if (add_to_checkpoint_control)
        {
          add_object(identifier, object);
        }
      }

      /**
         * \brief Write checkpoint files to disk
         *
         * Write two uncompressed binary files:
         * filename: [name.cp]: holding the current status of all checkpointable objects added to the checkpoint control
         *           [name.scp]: holding the size of the checkpoint data per rank, needed for loading the checkpoint
         * or write two compressed binary files:
         * filename: [name.zcp]: a compressed file holding the current status of all checkpointable objects added to the checkpoint control, needs zlib support
         *           [name.szcp]: holding the size of the checkpoint data per rank, needed for loading the checkpoint
         *
         * \param[in] filename String holding the complete name of the file with extension .cp, without compression, or .zcp, for compression
         */
      void save(const String filename)
      {
        size_t pos = filename.rfind('.');
        String extension = filename.substr(pos + 1);
        String name = filename.substr(0, pos);

        XASSERTM(name != "", "no complete filename consisting of name.extension given");

        if (extension == "cp")
        {
          _save(name);
        }
#ifdef FEAT_HAVE_ZLIB
        else if (extension == "zcp")
        {
          _save_compressed(name);
        }
#endif // FEAT_HAVE_ZLIB
        else
        {
          XASSERTM(extension != "zcp", "no zlib support, activate zlib or write uncompressed files with .cp extension instead");
          XASSERTM(extension == "zcp", "no valid checkpoint filename choosen");
        }
      }

      /**
         * \brief Load checkpoint from disk
         *
         * Read the [filename.sz] file, to restore the size of the checkpoint data per rank.
         *
         * Read binary file [filename] holding a safed state of all checkpointable objects and store it to _input_array.
         *
         * Fill the _offset_by_identifier map with the data from the array.
         *
         * \param[in] filename String holding the name of the file
         *
         * \warning Reading another input file / stream will overwrite values in _input_array and _offset_by_identifier
         */
      void load(const String filename)
      {
        XASSERTM(_input_array == nullptr, "another input file was read before");

        size_t pos = filename.rfind('.');
        String extension = filename.substr(pos + 1);
        String name = filename.substr(0, pos);

        XASSERTM(name != "", "no complete filename consisting of name.extension given");

        if (extension == "cp")
        {
          _load(name);
        }
#ifdef FEAT_HAVE_ZLIB
        else if (extension == "zcp")
        {
          _load_compressed(name);
        }
#endif // FEAT_HAVE_ZLIB
        else
        {
          XASSERTM(extension != "zcp", "no zlib support, activate zlip or read uncompressed files with .cp extension instead");
          XASSERTM(extension == "zcp", "no valid checkpoint file choosen");
        }
      }

#ifdef FEAT_HAVE_ZLIB
      /**
         * \brief Save checkpoint to a stream
         *
         * Save a checkpoint, holding the current status of all registered objects, in a BinaryStream
         *
         * \param[in] bs BinaryStream that shall be written to
         */
      void save(BinaryStream & bs)
      {
        std::vector<char> data;

        uLongf checkpoint_size = static_cast<uLongf>(_collect_checkpoint_data(data));

        uLongf compressedLen = static_cast<uLongf>(compressBound(checkpoint_size));

        Bytef * compressed = new Bytef[compressedLen];
        Bytef * uncompressed = (Bytef *)data.data();
        int c_status = ::compress(compressed, &compressedLen, uncompressed, checkpoint_size);
        XASSERTM(c_status == Z_OK, "compression of checkpoint data failed");

        std::uint64_t clen = compressedLen;
        std::uint64_t slen = checkpoint_size;

        bs.write(reinterpret_cast<char *>(&clen), sizeof(clen));
        bs.write(reinterpret_cast<char *>(&slen), sizeof(slen));
        bs.write(reinterpret_cast<char *>(compressed), static_cast<std::streamsize>(clen));

        delete[] compressed;
      }

      /**
         * \brief Load checkpoint from stream
         *
         * Read the BinaryStream, to restore the safed state of all objects registered to the checkpoint control
         * and store it to _input_array.
         *
         * Fill the _offset_by_identifier map with the data from the array.
         *
         * \param[in] bs BinaryStream that shall be read from
         *
         * \warning Reading another input file / stream will overwrite values in _input_array and _offset_by_identifier
         */
      void load(BinaryStream & bs)
      {
        XASSERTM(_input_array == nullptr, "another input file was read before");

        char * data = bs.data();
        uLongf compressed_size = static_cast<uLongf>(*(std::uint64_t *)(data));
        uLongf checkpoint_size = static_cast<uLongf>(*(std::uint64_t *)(data + sizeof(uint64_t)));

        Bytef * compressed = (Bytef *)(data + (2 * sizeof(uint64_t)));
        _input_array = new char[checkpoint_size];
        int d_status = ::uncompress((Bytef *)_input_array, &checkpoint_size, compressed, compressed_size);
        XASSERTM(d_status == Z_OK, "decompression of checkpoint data failed");

        _restore_checkpoint_data(static_cast<uint64_t>(checkpoint_size));
      }
#else //FEAT_HAVE_ZLIB
      /**
         * \brief Save checkpoint to a stream
         *
         * Save a checkpoint, holding the current status of all registered objects, in a BinaryStream
         *
         * \param[in] bs BinaryStream that shall be written to
         */
      void save(BinaryStream & bs)
      {
        std::vector<char> buffer;
        auto checkpoint_size = _collect_checkpoint_data(buffer);

        std::uint64_t slen = checkpoint_size;
        bs.write(reinterpret_cast<char *>(&slen), sizeof(slen));
        bs.write(buffer.data(), static_cast<std::streamsize>(buffer.size()));
      }

      /**
         * \brief Load checkpoint from stream
         *
         * Read the BinaryStream, to restore the safed state of all objects registered to the checkpoint control
         * and store it to _input_array.
         *
         * Fill the _offset_by_identifier map with the data from the array.
         *
         * \param[in] bs BinaryStream that shall be read from
         *
         * \warning Reading another input file / stream will overwrite values in _input_array and _offset_by_identifier
         */
      void load(BinaryStream & bs)
      {
        XASSERTM(_input_array == nullptr, "another input file was read before");

        char * buffer = bs.data();
        uint64_t size = *(std::uint64_t *)(buffer);

        _input_array = new char[size];
        std::copy(buffer + sizeof(uint64_t), buffer + sizeof(uint64_t) + size - 1, _input_array);

        _restore_checkpoint_data(size);
      }
#endif //FEAT_HAVE_ZLIB

    }; // class Checkpoint

  } // namespace Control
} // namespace FEAT

#endif // CONTROL_CHECKPOINT_HPP
