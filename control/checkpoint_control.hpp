// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
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
#include <kernel/lafem/container.hpp>

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
      virtual std::uint64_t get_checkpoint_size(LAFEM::SerialConfig & config) = 0;

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
      virtual std::uint64_t set_checkpoint_data(std::vector<char> & data, LAFEM::SerialConfig & config) = 0;

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
      virtual std::uint64_t get_checkpoint_size(LAFEM::SerialConfig & config) override
      {
        return _object.get_checkpoint_size(config);
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
      virtual std::uint64_t set_checkpoint_data(std::vector<char> & data, LAFEM::SerialConfig & config) override
      {
        return _object.set_checkpoint_data(data, config);
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
      /// Vector holding the array read from the input file during restore
      std::vector<char> _input_array;
      /// Mapping of identifier string to the offset in the input file
      std::map<String, std::uint64_t> _offset_by_identifier;
      /// the mpi communicator identifying our mpi context
      const Dist::Comm & _comm;
      /// The config class that controls the compression modes
      FEAT::LAFEM::SerialConfig _config;

      /**
         * \brief Build checkpoint buffer
         *
         * Collect all data needed for a checkpoint from every object registered at the checkpoint control.
         *
         * \param[in] buffer buffer containing the collected data
         * \returns the size of the buffer
         *
         */
      std::uint64_t _collect_checkpoint_data(std::vector<char> & buffer)
      {
        std::uint64_t checkpoint_size(0u);
        std::uint64_t real_size(0u);
        std::map<String, std::tuple<std::uint64_t, std::uint64_t>> sizes;

        for (auto const & it : _checkpointable_by_identifier)
        {
          std::uint64_t identifierlength = (std::uint64_t)it.first.length();
          std::uint64_t datalength = it.second->get_checkpoint_size(_config);

          checkpoint_size += identifierlength + datalength + sizeof(std::uint64_t) + sizeof(std::uint64_t);
          real_size += identifierlength + sizeof(std::uint64_t) + sizeof(std::uint64_t);
          sizes[it.first] = std::make_tuple(identifierlength, datalength);
        }

        buffer.reserve(checkpoint_size);

        for (auto const & it : _checkpointable_by_identifier)
        {
          char * cidentifierlength = reinterpret_cast<char *>(&std::get<0>(sizes[it.first]));
          buffer.insert(std::end(buffer), cidentifierlength, cidentifierlength + sizeof(std::uint64_t));
          buffer.insert(std::end(buffer), it.first.begin(), it.first.end());
          std::uint64_t old_size = buffer.size();
          buffer.insert(std::end(buffer), sizeof(std::uint64_t), 0); //set datalength to an arbitrary value
          std::uint64_t ireal_size = it.second->set_checkpoint_data(buffer, _config);
          char * csize = reinterpret_cast<char *>(&ireal_size);
          for(uint i(0) ; i < sizeof(std::uint64_t) ; ++i)  //overwrite the the datalength
          {
            buffer[old_size + i] = csize[i];
          }
          real_size += ireal_size;
        }

        buffer.resize(real_size);
        return real_size;
      }

      /**
         * \brief Extract input buffer
         *
         * Extract the identifier string and offset for its data in the input array for every object
         * in the read checkpoint, and store this information for restoration in this checkpoint control object.
         *
         */
      void _restore_checkpoint_data()
      {
        std::uint64_t stringsize;
        std::uint64_t datasize(0);
        std::size_t i = 0;
        std::size_t size = _input_array.size();

        // Loop over all objects stored in the checkpoint file
        while (i < size)
        {
          // Get size of identifier string
          ::memcpy(&stringsize, _input_array.data() + i, sizeof(std::uint64_t));
          i += sizeof(std::uint64_t);

          // Read the identifier string and put it as key for the calculated offset to the map
          _offset_by_identifier[String(_input_array.data() + i, stringsize)] = i + stringsize;
          i += stringsize;

          // Get the size of the Data holding the array, to get the start point of the next checkpointable object in the array
          ::memcpy(&datasize, _input_array.data() + i, sizeof(std::uint64_t));
          i += sizeof(std::uint64_t) + datasize;
        }
      }

      /**
         * \brief Write checkpoint files to disk
         *
         * Write one binary file:
         * [name].cp: a uncompressed binary file holding the current status of all checkpointable objects added to the checkpoint control
         *
         * \param[in] name String holding the name of the file (without the extension .cp)
         */
      void _save(const String name)
      {
        std::vector<char> buffer;
        _collect_checkpoint_data(buffer);
        std::vector<char> common(0);

        DistFileIO::write_combined(common, buffer, name + ".cp", _comm);
      }

      /**
         * \brief Load checkpoint from disk
         *
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

        //std::vector<char> temp(0);
        std::vector<char> common(0);
        // read the checkpoint file
        DistFileIO::read_combined(common, _input_array, name + ".cp", _comm);

        _restore_checkpoint_data();
      }

    public:
      /**
         * \brief Constructor
         *
         * Initalize the input array as NULL pointer.
         *
         * \param[in] comm The communicator common to all stored objects
         * \param[in] config A config class, controlling the compression of the written out data
         */
      CheckpointControl(const Dist::Comm & comm, const LAFEM::SerialConfig & config = LAFEM::SerialConfig()) :
        _comm(comm),
        _config(config)
      {
        _input_array = std::vector<char>(0);
      }

      /**
         * \brief Destructor
         *
         * Destroy the checkpoint control instance and delete the input array.
         */
      ~CheckpointControl()
      {
      }

      /**
       * \brief Set a new serialize configuration
       *
       * \param[in] conf LAFEM::SerialConfig, a config class holding information about the compression parameters.
       *
       * Set a new SerialConfig object.
       */
      void set_config(FEAT::LAFEM::SerialConfig& conf)
      {
        _config = conf;
      }
      /**
         * \brief Delete all read input
         *
         * Clear the input read before from this checkpoint control object.
         * Does not touch the map holding the pointers to the checkpointable objects included into the following written checkpoints.
         */
      void clear_input()
      {
        _offset_by_identifier.clear();
        _input_array.resize(0);
      }

      /// Retrieve a list of all items stored in the checkpoint
      String get_identifier_list()
      {
        String identifiers;
        for (auto item : _checkpointable_by_identifier)
        {
          if (identifiers.length() > 0)
          {
            identifiers += "\n";
          }
          identifiers += item.first;
        }
        return identifiers;
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
        XASSERTM(_input_array.size() > 0u, "no input file has been loaded before");
        XASSERTM(_offset_by_identifier.count(identifier) == 1, "input file doesn't contain data for the given identifier");

        auto checkpointable_object = std::make_shared<CheckpointableWrapper<OT_>>(object);

        std::uint64_t size(0);
        char * in_data = _input_array.data() + _offset_by_identifier[identifier];
        ::memcpy(&size, in_data, sizeof(std::uint64_t));
        std::vector<char> buffer(in_data + sizeof(std::uint64_t), in_data + sizeof(std::uint64_t) + size);
        checkpointable_object->restore_from_checkpoint_data(buffer);

        if (add_to_checkpoint_control)
        {
          add_object(identifier, object);
        }
      }

      /**
         * \brief Write checkpoint file to disk
         *
         * Write one uncompressed binary file:
         * filename: [name.cp]: holding the current status of all checkpointable objects added to the checkpoint control
         *
         * \param[in] filename String holding the complete name of the file with extension .cp
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
        else
        {
          XASSERTM(extension != "zcp", "no zlib support, activate zlib or write uncompressed files with .cp extension instead");
          XASSERTM(extension == "zcp", "no valid checkpoint filename choosen");
        }
      }

      //Should we ignore the extension and just give in filename without .cp?
      /**
         * \brief Load checkpoint from disk
         *
         * Read the [filename.cp] file, to restore the size of the checkpoint data per rank.
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
        XASSERTM(_input_array.size() == 0u, "another input file was read before");

        size_t pos = filename.rfind('.');
        String extension = filename.substr(pos + 1);
        String name = filename.substr(0, pos);

        XASSERTM(name != "", "no complete filename consisting of name.extension given");

        if (extension == "cp")
        {
          _load(name);
        }
        else
        {
          XABORTM("no valid checkpoint file choosen");
        }
      }

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
        std::uint64_t slen = _collect_checkpoint_data(buffer);

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
        XASSERTM(_input_array.size() == 0, "another input file was read before");

        char * buffer = bs.data();
        std::uint64_t size = *(std::uint64_t *)(buffer);

        _input_array.resize(size);
        std::copy(buffer + sizeof(std::uint64_t), buffer + sizeof(std::uint64_t) + size - 1, _input_array.data());

        _restore_checkpoint_data();
      }

    }; // class Checkpoint

  } // namespace Control
} // namespace FEAT

#endif // CONTROL_CHECKPOINT_HPP
