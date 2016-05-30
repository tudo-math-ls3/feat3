#pragma once
#ifndef KERNEL_ADJACENCY_EXPORT_TGA_HPP
#define KERNEL_ADJACENCY_EXPORT_TGA_HPP 1

// includes, FEAT
#include <kernel/adjacency/adjactor.hpp>

// includes, system
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdint>

namespace FEAT
{
  namespace Adjacency
  {
    /**
     * \brief Truevision TGA exporter
     *
     * This class can export the adjacency structure of any object implementing the
     * Adjactor interface into a monochrome TGA image file, which can then be admired
     * in a image viewing application.
     *
     * Important notes:
     * - The domain/image dimensions are limited to 32767. Trying to export a larger
     *   adjactor will result in an exception being thrown.
     * - Please keep in mind that uncompressed image files become very large very quickly,
     *   e.g. adjactor with an domain/image dimension of ~16000 will result in a TGA file
     *   size of ~250 MB!
     * - This exporter is implemented in an <em>endian-invariant</em> way, i.e. it should
     *   work correctly on both little-endian and big-endian platforms.
     *
     * \author Peter Zajac
     */
    class ExportTGA
    {
    public:
      /**
       * \brief Writes out an adjactor to a TGA image file.
       *
       * \param[in] filename
       * The filename for the TGA image file.
       *
       * \param[in] adj
       * The adjactor that is to be rendered into the image file.
       */
      template<typename Adjactor_>
      static void write(const String& filename, const Adjactor_& adj)
      {
        // try to open output file
        std::ofstream ofs(filename, std::ios_base::binary|std::ios_base::out);
        if(!(ofs.is_open() && ofs.good()))
          throw InternalError(String("Failed to open '") + filename + "'");

        // write
        write(ofs, adj);

        // close file
        ofs.close();
      }

      /**
       * \brief Writes out an adjactor to a TGA image file.
       *
       * \param[in] os
       * A binary output stream that receives the rendered adjactor.
       *
       * \param[in] adj
       * The adjactor that is to be rendered into the image file.
       */
      template<typename Adjactor_>
      static void write(std::ostream& os, const Adjactor_& adj)
      {
        typedef std::uint8_t u8;

        // set up header
        u8 header[18];
        for(int i(0); i < 18; ++i)
          header[i] = 0;

        // get dimensions
        const Index w = adj.get_num_nodes_image();
        const Index h = adj.get_num_nodes_domain();

        // make sure we do not exceed the 16-bit range
        if((w > Index(32767)) || (h > Index(32767)))
          throw InternalError("TGA dimensions are limited to 32767");

        // set dimensions
        header[12] = u8( w       & 0xFF);
        header[13] = u8((w >> 8) & 0x7F);
        header[14] = u8( h       & 0xFF);
        header[15] = u8((h >> 8) & 0x7F);

        // set basic stuff
        header[ 2] = u8( 3); // datatype code
        header[16] = u8( 8); // bits per pixel
        header[17] = u8(32); // image descriptor

        // write header to file
        os.write(reinterpret_cast<char*>(header), 18);

        // allocate a scan-line
        std::vector<char> scan(std::size_t(w), -1);

        // okay, let's loop over all domain nodes
        for(Index i(0); i < h; ++i)
        {
          // get the image iterators
          auto jt = adj.image_end(i);

          // loop over all image nodes and set black pixels
          for(auto it = adj.image_begin(i); it != jt; ++it)
          {
            // set pixel to black
            scan[*it] = 0;
          }

          // write scan-line
          os.write(scan.data(), std::streamsize(scan.size()));

          // clear pixels
          for(auto it = adj.image_begin(i); it != jt; ++it)
          {
            scan[*it] = -1;
          }
        }

        // finally, write the TGA v2 footer
        const char* footer = "\0\0\0\0\0\0\0\0TRUEVISION-XFILE.\0";
        os.write(footer, 26);
      }
    }; // class ExportTGA
  } // namespace Adjacency
} // namespace FEAT

#endif // KERNEL_ADJACENCY_EXPORT_TGA_HPP
