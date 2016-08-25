#pragma once
#ifndef KERNEL_UTIL_BINARY_STREAM_HPP
#define KERNEL_UTIL_BINARY_STREAM_HPP 1

// includes, system
#include <vector>
#include <iostream>

namespace FEAT
{
  /**
   * \brief Binary Stream class
   *
   * This class implements a dynamic binary memory input/output stream implementing the std::iostream interface.
   *
   * \warning This class is not yet fully matured and shall therefore be used with care!
   *
   * This class can be used to write and read binary data into/from memory to test functions which usually
   * use the std::[i|o]fstream classes for input and output from/into files.
   *
   * This class is the binary counterpart of the text-based std::stringstream class.
   *
   * \author Peter Zajac
   */
  class BinaryStream :
    public std::iostream
  {
  private:
    /// \cond internal
    /**
     * \brief Binary Stream Buffer class
     *
     * This class implements the std::basic_streambuf<char> interface using a std::vector as a data container.
     *
     * \author Peter Zajac
     */
    class BinaryBuffer :
      public std::basic_streambuf<char>
    {
    private:
      typedef std::vector<char>::size_type datasize;
      typedef std::streamsize streamsize;

      /// the data vector
      std::vector<char> _data;
      /// the current stream position
      datasize _data_pos;

    public:
      BinaryBuffer() :
        _data_pos(0)
      {
      }

      streamsize size() const
      {
        return streamsize(_data.size());
      }

      void clear()
      {
        _data.clear();
        _data_pos = datasize(0);
      }

      char* data()
      {
        return _data.data();
      }

      const char* data() const
      {
        return _data.data();
      }

      std::vector<char>& container()
      {
        return _data;
      }

      const std::vector<char>& container() const
      {
        return _data;
      }

    protected:
      virtual std::streamsize xsputn(const char* ptr, std::streamsize count) override
      {
        streamsize copied(0);

        // Ensure that the vector has sufficient length; if not, then resize the vector
        if(_data_pos + datasize(count) > _data.size())
          _data.resize(_data_pos + datasize(count));

        // copy the data
        for(; copied < count; ++_data_pos, ++copied, ++ptr)
          _data[_data_pos] = *ptr;

        return copied;
      }

      virtual std::streamsize xsgetn(char* ptr, std::streamsize count) override
      {
        streamsize copied(0);
        datasize len(_data.size());

        // Ensure that we do not read beyond the end of the buffer.
        if(_data_pos + datasize(count) > len)
          count = streamsize(len - _data_pos);

        // copy the data
        for(; copied < count; ++_data_pos, ++copied, ++ptr)
          *ptr = _data[_data_pos];

        return copied;
      }

      virtual pos_type seekoff(
        off_type off,
        std::ios_base::seekdir way,
        std::ios_base::openmode which = std::ios_base::in | std::ios_base::out) override
      {
        if((which & (std::ios_base::in|std::ios_base::out)) == 0)
          return pos_type(-1);

        switch(way)
        {
        case std::ios_base::beg:
          if((off < 0) || (off > off_type(_data.size())))
            return pos_type(-1);
          return pos_type(std::streamoff(_data_pos = datasize(off)));

        case std::ios_base::cur:
          if((off < 0) && (-off > off_type(_data_pos)))
            return pos_type(-1);
          if((off > 0) && (off_type(_data_pos) + off > off_type(_data.size())))
            return pos_type(-1);
          return pos_type(std::streamoff(_data_pos = datasize(off_type(_data_pos) + off)));

        case std::ios_base::end:
          if((off > 0) || (-off > off_type(_data.size())))
            return pos_type(-1);
          return pos_type(std::streamoff(_data_pos = datasize(off_type(_data.size()) + off)));

        default:
          return pos_type(-1);
        }
      }

      virtual pos_type seekpos(
        pos_type pos,
        std::ios_base::openmode which = std::ios_base::in | std::ios_base::out) override
      {
        if((which & (std::ios_base::in|std::ios_base::out)) == 0)
          return pos_type(-1);

        // Ensure the position is not out-of-bounds
        if (pos >= pos_type(std::streamoff(_data.size())))
          return pos_type(-1);

        // move stream position
        return pos_type(std::streamoff(_data_pos = datasize(pos)));
      }
    };
    /// \endcond

    // the stream's buffer
    BinaryBuffer _buffer;

  public:
    /// default CTOR
    BinaryStream() :
      std::iostream(&_buffer)
    {
    }

    /**
     * \brief Returns the current size of the stream.
     * \returns The current size of the stream in bytes.
     */
    std::streamsize size() const
    {
      return _buffer.size();
    }

    /**
     * \brief Clears the current stream.
     */
    void clear()
    {
      _buffer.clear();
    }

    /**
     * \brief Returns the data array of the stream.
     * \returns The data array of the stream.
     */
    char* data()
    {
      return _buffer.data();
    }

    /**
     * \brief Returns the data array of the stream.
     * \returns The data array of the stream.
     */
    const char* data() const
    {
      return _buffer.data();
    }

    /**
     * \brief Returns a reference to the internal vector container.
     * \returns A reference to the internal vector container.
     */
    std::vector<char>& container()
    {
      return _buffer.container();
    }

    /**
     * \brief Returns a reference to the internal vector container.
     * \returns A reference to the internal vector container.
     */
    const std::vector<char>& container() const
    {
      return _buffer.container();
    }
  }; // class BinaryStream
} // namespace FEAT
#endif // KERNEL_UTIL_BINARY_STREAM_HPP
