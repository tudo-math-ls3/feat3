// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/geometry/cgal.hpp>
#include <kernel/geometry/voxel_map.hpp>
#include <kernel/util/stop_watch.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdint>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace VoxelMapGenerator
{
  using namespace FEAT;

  int display_help(Dist::Comm& comm, SimpleArgParser& args)
  {
    comm.print(
    //          1         2         3         4         5         6         7         8         9         0
    // 123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
      "\nWelcome to FEAT's VoxelMap generation and analysis tool\n"
      "\n"
      "This tool can be used to create new or read in existing voxel map files to analyze, export,\n"
      "render or write out (compressed) voxel map files.\n"
      "To read in an existing voxel map file, use the '--in <file>' option.\n"
      "To create a new voxel map, you have to perform three basic steps:\n"
      "* First, specify the domain's bounding box by using the '--box ...' option.\n"
      "* Second, specify the voxel map resolution by either the '--num ...' or '--res ...' options.\n"
      "* Third, create the actual map from...\n"
      "  - a surface triangulation stored in an OFF file by using the '--off <file>' option.\n"
      "  - a formula in x,y,z coordinates by using the '--formula ...' option, in which case a voxel is\n"
      "    considered to be inside if the formula evaluates to a value > 0.\n"
      "\n"
      "The following options are supported by this tool:\n"
    );
    comm.print(args.get_supported_help());
    return 0;
  }

  int check_voxel_maps(Dist::Comm& comm, SimpleArgParser& args)
  {
    typedef std::uint64_t u64;

    // only rank 0 is working here
    if(comm.rank() > 0)
      return 0;

    // do all files look valid?
    bool seem_all_valid = true;

    // loop over all filenames
    for(const String& filename : args.query("check")->second)
    {
      std::cout << "Reading voxel map file '" << filename << "'...\n";

      std::ifstream ifs(filename, std::ios_base::in | std::ios_base::binary);
      if(!ifs.is_open() || !ifs.good())
      {
        std::cout << "ERROR: Failed to open '" << filename << "'!\n";
        return 1;
      }

      // get the file size
      ifs.seekg(0ll, std::ios_base::end);
      const u64 filesize = u64(ifs.tellg());
      ifs.seekg(0ll, std::ios_base::beg);

      // at the very least, we have to be able to read a header
      if(filesize < Geometry::VoxelMap::header_size)
      {
        std::cout << "ERROR: File is too small to even contain the voxel map header!\n";
        return 1;
      }

      // allocate a buffer of appropriate size
      std::vector<char> buffer(filesize);
      ifs.read(buffer.data(), std::streamsize(filesize));
      ifs.close();

      // try to interpret the header
      const Geometry::VoxelMap::FileHeader& header = *reinterpret_cast<Geometry::VoxelMap::FileHeader*>(buffer.data());

      // check magic number
      if(header.magic_number != Geometry::VoxelMap::magic_number)
      {
        std::cout << "ERROR: File has incorrect magic number; this does not seem to be a valid voxel map file!\n";
        return 1;
      }

      // check header size
      if(header.header_size != Geometry::VoxelMap::header_size)
      {
        std::cout << "ERROR: File has incorrect header size; maybe this file was created for a newer version?\n";
        std::cout << "Expected header size: " << Geometry::VoxelMap::header_size << " but found " << header.header_size << "\n";
        return 1;
      }

      // ok, it seems we have a valid voxel map!

      // compute unit scaling factor
      const double unit_scale = 1.0 / double(Geometry::VoxelMap::unit_size);

      // print dimensions
      std::cout << "Minimum X coord....: " << stringify_fp_fix(double(header.min_x) * unit_scale, 6, 12) << "\n";
      std::cout << "Maximum X coord....: " << stringify_fp_fix(double(header.max_x) * unit_scale, 6, 12) << "\n";
      std::cout << "Resolution X.......: " << stringify_fp_fix(header.num_x > 0u ? double(header.max_x - header.min_x) / double(header.num_x) * unit_scale : 0.0, 6, 12) << "\n";
      std::cout << "Number of Points X.: " << stringify(header.num_x).pad_front(12) << "\n";
      std::cout << "Minimum Y coord....: " << stringify_fp_fix(double(header.min_y) * unit_scale, 6, 12) << "\n";
      std::cout << "Maximum Y coord....: " << stringify_fp_fix(double(header.max_y) * unit_scale, 6, 12) << "\n";
      std::cout << "Resolution Y.......: " << stringify_fp_fix(header.num_y > 0u ? double(header.max_y - header.min_y) / double(header.num_y) * unit_scale : 0.0, 6, 12) << "\n";
      std::cout << "Number of Points Y.: " << stringify(header.num_y).pad_front(12) << "\n";
      std::cout << "Minimum Z coord....: " << stringify_fp_fix(double(header.min_z) * unit_scale, 6, 12) << "\n";
      std::cout << "Maximum Z coord....: " << stringify_fp_fix(double(header.max_z) * unit_scale, 6, 12) << "\n";
      std::cout << "Resolution Z.......: " << stringify_fp_fix(header.num_z > 0u ? double(header.max_z - header.min_z) / double(header.num_z) * 1E-6 : 0.0, 6, 12) << "\n";
      std::cout << "Number of Points Z.: " << stringify(header.num_z).pad_front(12) << "\n";
      std::cout << "Number of Voxels...: " << stringify(header.num_x * Math::max(header.num_y, u64(1)) * Math::max(header.num_z, u64(1))).pad_front(12) << "\n";
      std::cout << "Stride Line........: " << stringify(header.stride_line).pad_front(12) << " Bytes\n";
      std::cout << "Stride Plane.......: " << stringify(header.stride_plane).pad_front(12) << " Bytes\n";
      std::cout << "Stride Volume......: " << stringify(header.stride_volume).pad_front(12) << " Bytes\n";
      double box_vol = 1.0;
      if(header.num_x > 0u)
        box_vol *= double(header.max_x - header.min_x) * unit_scale;
      if(header.num_y > 0u)
        box_vol *= double(header.max_y - header.min_y) * unit_scale;
      if(header.num_z > 0u)
        box_vol *= double(header.max_z - header.min_z) * unit_scale;
      std::cout << "Bounding Box Volume: " << stringify_fp_sci(box_vol, 3, 12).pad_front(12) << "\n";
      std::cout << "Domain Volume......: " << stringify_fp_sci(box_vol * double(header.coverage) * unit_scale, 3, 12) << "\n";
      std::cout << "Domain Coverage....: " << stringify_fp_fix(double(header.coverage) * unit_scale, 6, 12) << "\n";

      //std::cout << "Out-Of-Bounds Value: " << (header.flags & Geometry::VoxelMap::file_header_flags_oob_value_true ? "true" : "false") << "\n";
      std::cout << "Planes per Block...: " << stringify(header.planes_per_block).pad_front(12) << "\n";
      std::cout << "Number of Blocks...: " << stringify(header.num_blocks).pad_front(12) << "\n";

      // let's assume that's a valid voxel map file
      bool seems_valid = true;

      // perform some sanity checks before we continue
      if((header.stride_line & 0xFu) != 0u)
      {
        std::cout << "WARNING: stride_line is not a multiple of 16!\n";
        seems_valid = false;
      }
      if((header.stride_plane & 0xFu) != 0u)
      {
        std::cout << "WARNING: stride_plane is not a multiple of 16!\n";
        seems_valid = false;
      }
      if((header.stride_volume & 0xFu) != 0u)
      {
        std::cout << "WARNING: stride_volume is not a multiple of 16!\n";
        seems_valid = false;
      }
      if(header.stride_line*8u < header.num_x)
      {
        std::cout << "ERROR: stride_line is too small to contain an entire voxel X-line!\n";
        seems_valid = false;
      }
      if(header.stride_plane < header.num_y * header.stride_line)
      {
        std::cout << "ERROR: stride_plane is too small to contain an entire voxel XY-plane!\n";
        seems_valid = false;
      }
      if(header.stride_volume < header.num_z * header.stride_plane)
      {
        std::cout << "ERROR: stride_volume is too small to contain an entire voxel XYZ-volume!\n";
        seems_valid = false;
      }
      if((header.num_x > 0u) && (header.max_x <= header.min_x))
      {
        std::cout << "ERROR: max_x is smaller than min_x!\n";
        seems_valid = false;
      }
      if((header.num_y > 0u) && (header.max_y <= header.min_y))
      {
        std::cout << "ERROR: max_y is smaller than min_y!\n";
        seems_valid = false;
      }
      if((header.num_z > 0u) && (header.max_z <= header.min_z))
      {
        std::cout << "ERROR: max_z is smaller than min_z!\n";
        seems_valid = false;
      }
      if((header.num_blocks > 0u) && (header.planes_per_block == 0u))
      {
        std::cout << "ERROR: num_blocks is > 0 but planes_per_block is 0!\n";
        seems_valid = false;
      }
      if((header.num_blocks == 0u) && (header.planes_per_block > 0u))
      {
        std::cout << "ERROR: planes_per_block is > 0 but num_compress is 0!\n";
        seems_valid = false;
      }
      if(header.num_z == 0u) // 2D case
      {
        if(header.num_blocks > 1u)
        {
          std::cout << "ERROR: num_blocks is > 1 but there is only 1 plane in 2D!\n";
          seems_valid = false;
        }
        /*if(header.planes_per_block > 1u)
        {
          std::cout << "ERROR: planes_per_block is > 1 but there is only 1 plane in 2D!\n";
          seems_valid = false;
        }*/
      }
      else if(header.num_blocks > 0u) // 3D case
      {
        if(header.num_blocks > header.num_z)
        {
          std::cout << "ERROR: num_blocks is bigger than num_z!\n";
          seems_valid = false;
        }
        /*if(header.planes_per_block > header.num_z)
        {
          std::cout << "ERROR: planes_per_block is bigger than num_z!\n";
          seems_valid = false;
        }*/
        if(header.num_blocks*header.planes_per_block < header.num_z)
        {
          std::cout << "ERROR: num_blocks*planes_per_block is less than num_z!\n";
          seems_valid = false;
        }
      }

      // check minimum file size
      u64 total_head_size = Geometry::VoxelMap::header_size + header.num_blocks * sizeof(u64);
      if(filesize < total_head_size)
      {
        std::cout << "ERROR: file is too small to contain header and compression blocks!\n";
        seems_valid = false;
      }

      // do we have compression blocks?
      u64 compress_size = 0u;
      if(header.planes_per_block > 0u)
      {
        const u64* blocks = reinterpret_cast<const u64*>(&buffer.data()[header.header_size]);
        for(u64 iblock = 0; iblock < header.num_blocks; ++iblock)
        {
          compress_size += blocks[iblock];
        }
      }

      // compute total map size
      u64 map_size = (header.planes_per_block > 0u ? compress_size : header.stride_volume);
      if(filesize < total_head_size + map_size)
      {
        std::cout << "ERROR: file is too small to contain voxel map data!\n";
        seems_valid = false;
      }

      // print the rest
      std::cout << "Voxel Map Size.....: " << stringify_bytes(header.stride_volume, 3, 12) << "\n";
      std::cout << "Total File Size....: " << stringify_bytes(filesize, 3, 12) << "\n";
      if(compress_size > 0u)
      {
        std::cout << "Compressed Map Size: " << stringify_bytes(compress_size, 3, 12) << "\n";
        std::cout << "Compression Rate...: " << stringify_fp_fix(100.0*double(compress_size) / double(header.stride_volume), 3, 12) << " %\n";
      }
      else
      {
        std::cout << "Compressed Map Size: -N/A-\n";
        std::cout << "Compression Rate...: -N/A-\n";
      }
      std::cout << "\n";

      if(!seems_valid)
        seem_all_valid = false;
    }

    return (seem_all_valid ? 0 : 1);
  }

  int main(int argc, char** argv)
  {
    Dist::Comm comm(Dist::Comm::world());

#ifndef FEAT_HAVE_ZLIB
    comm.print("WARNING: ZLIB library not found; you will not be able to read/write compressed voxel maps!");
#endif

#ifndef FEAT_HAVE_CGAL
    comm.print("WARNING: CGAL library not found; you will not be able to create voxel maps from OFF files!");
#endif

#ifndef FEAT_HAVE_FPARSER
    comm.print("WARNING: fparser library not found; you will not be able to create voxel maps from formulae!");
#endif

    SimpleArgParser args(argc, argv);
    args.support("help");
    args.support("check", "<voxel-map-files...>\n"
      "Reads in a voxel map and prints some basic information about its contents.");
    args.support("out", "<voxel-map>\n"
      "Specifies the output voxel map filename; the recommended extension is '.vxl'.");
    args.support("in", "<voxel-map>\n"
      "Specifies the input voxel map filename.");
    args.support("off", "<off-file>\n"
      "Specifies an OFF file to generate the voxel map from.\n"
      "Requires that FEAT is linked against the 'CGAL' third-party library.");
    args.support("formula", "<formula>\n"
      "Specifies a formula in x,y,z to generate the voxel map from; a point is considered to be inside\n"
      "the domain if the formula evaluates to a value > 0.\n"
      "Requires that FEAT is linked against the 'fparser' third-party library.");
    args.support("compress", "[<block-size> [<compression-level>]]\n"
      "Specifies that the voxel map should be compressed. The optional <block-size> argument specifies\n"
      "the maximum size of a single compression block in MB; this value defaults to 128. The second\n"
      "optional <compression-level> argument specifies the desired ZLIB compression level in range [1,9]\n"
      "where 1 is the fastest compression and 9 is the best compression; this value defaults to 9.\n"
      "Requires that FEAT is linked against the 'ZLIB' third-party library.");
    args.support("box", "<x-min> <x-max> <y-min> <y-max> [<z-min> <z-max>]\n"
      "Specifies the minimum and maximum X-, Y- and Z-coordinates.");
    args.support("num", "<num-x> <num-y> [<num-z>]\n"
      "Specifies the number of points in X-, Y- and Z-dimension. Mutually exclusive with --res.");
    args.support("res", "<delta>\n"
      "Specifies the maximum voxel resolution in all dimensions. Mutually exclusive with --num.");
    args.support("render", "<bmp-fileprefix> <x-res> <y-res> <z-res>\n"
      "Specifies the filename prefix and resolutions for BMP rendering.");
    args.support("export", "<bmp-fileprefix>\n"
      "Specifies the filename prefix for BMP export.");
    args.support("invert","\nIf specified, flips the definition of inside and outside.");
    args.support("threads", "<num-threads>\n"
      "Specifies how many threads are to be used (per MPI process)."
#ifndef _OPENMP
      "\nNote: This tool was compiled without OpenMP support, so multi-threading is not available."
#endif
      );

    int iarg = 0;

    // check for unsupported options
    auto unsupported = args.query_unsupported();
    if( !unsupported.empty() )
    {
      // print all unsupported options to cerr
      for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
      {
        comm.print(std::cerr, "ERROR: unsupported option '--" + (*it).second + "'");
      }

      display_help(comm, args);
      return 1;
    }

    // print help or info only?
    if((argc <= 1) || (args.check("help") >= 0))
      return display_help(comm, args);
    if(args.check("check") >= 0)
      return check_voxel_maps(comm, args);

    // create an empty voxel map object
    Geometry::VoxelMap voxel_map;

    // get input name
    String in_name;
    args.parse("in", in_name);

    // get output name
    String out_name;
    if((args.parse("out", out_name) > 0) && !out_name.ends_with(".vxl"))
      comm.print("INFO: It is recommended that Voxel Map files should use the extension .vxl");

    // get render name prefix
    String render_name;
    Index render_xres(0u), render_yres(0u), render_zres(0u);
    if((iarg = args.parse("render", render_name, render_xres, render_yres, render_zres)) < 0)
    {
      std::cerr << "ERROR: failed to parse --render parameter '" << args.get_arg(-iarg) << "' as resolution!\n";
      display_help(comm, args);
      return 1;
    }

    // get render name prefix
    String export_name;
    args.parse("export", export_name);

    // query OFF filename
    String off_name;
    args.parse("off", off_name);

    // query formula
    String formula;
    args.parse("formula", formula);

    // we'll figure this out later
    int dimension = 0;

    Geometry::VoxelMap::ReadResult read_result;

    if(!in_name.empty())
    {
      // read voxel map from file
      comm.print("\nReading voxel map from file '" + in_name + "'...");
      read_result = voxel_map.read(comm, in_name);
      dimension = (voxel_map.get_num_points(2) > 0u ? 3 : 2);
    }
    else // create voxel map somehow
    {
      if(off_name.empty() && formula.empty())
      {
        comm.print("ERROR: no input files specified via --in or --off and no formula specifies via --formula");
        display_help(comm, args);
        return 1;
      }

      if(args.check("box") < 4)
      {
        comm.print(std::cerr, "ERROR: mandatory option '--box <x-min> <x-max> <y-min> <y-max> [<z-min> <z-max>]' is missing!");
        display_help(comm, args);
        return 1;
      }

      // set bounding box
      Tiny::Vector<Real, 6> bbox;
      iarg = args.parse("box", bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5]);
      if(iarg == 6)
        voxel_map.set_bounding_box_3d(bbox[0], bbox[1], bbox[2], bbox[3], bbox[4], bbox[5]);
      else if(iarg == 4)
        voxel_map.set_bounding_box_2d(bbox[0], bbox[1], bbox[2], bbox[3]);
      else if(iarg >= 0)
      {
        std::cerr << "ERROR: invalid number of arguments for --box parameter: must be either 4 (2D) or 6 (3D)\n";
        return 1;
      }
      else
      {
        std::cerr << "ERROR: failed to parse --box parameter '" << args.get_arg(-iarg) << "' as coordinate!\n";
        return 1;
      }

      // remember what dimension we're dealing with
      dimension = iarg / 2; // either 2 or 3

      // set number of points
      if(args.check("num") > 0)
      {
        Index num[3];
        iarg = args.parse("num", num[0], num[1], num[2]);
        if(iarg < dimension)
        {
          std::cerr << "ERROR: failed to parse --num parameter '" << args.get_arg(-iarg) << "' as count!\n";
          return 1;
        }
        voxel_map.set_num_points(num[0], num[1], num[2]);
      }
      else if(args.check("res") > 0)
      {
        Real res(0.0);

        iarg = args.parse("res", res);
        if(iarg < 1)
        {
          std::cerr << "ERROR: failed to parse --res parameter '" << args.get_arg(-iarg) << "' as resolution!\n";
          return 1;
        }
        else if(res < 1e-6)
        {
          std::cerr << "ERROR: invalid resolution '" << res << "'!\n";
          return 1;
        }
        voxel_map.set_resolution(res);
      }
      else
      {
        std::cerr << "ERROR: either '--num <num-x> <num-y> <num-z>' or '--res <delta>' mus be given!\n";
        display_help(comm, args);
        return 1;
      }
    }

    // invert?
    bool invert = (args.check("invert") >= 0);

    // query compression settings
    Index compress_size = 128u;
    int compress_level = 9;
    bool compress_map = (args.check("compress") >= 0);
    if(args.parse("compress", compress_size, compress_level) < 0)
    {
      comm.print(std::cerr, "ERROR: failed to parse compression size/level");
      return 1;
    }
    // technically 0 is also a valid compression level
    if((compress_level < 1) || (compress_level > 9))
    {
      comm.print(std::cerr, "ERROR: invalid compression level " + stringify(compress_level));
      return 1;
    }
    if(!compress_map)
    {
      compress_size = 0u;
      compress_level = 0;
    }
#ifndef FEAT_HAVE_ZLIB
    if(compress_map)
    {
      comm.print(std::cerr, "ERROR: Cannot compress voxel map because FEAT has been compiled without ZLIB support");
      return 1;
    }
#endif

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // let's print a quick summary
    comm.print("Number of MPI Processes.: " + stringify(comm.size()));

    // specify number of threads?
#ifdef _OPENMP
    int num_threads(0);
    if(args.parse("threads", num_threads) > 0)
      omp_set_num_threads(num_threads);
    comm.print("Number of OpenMP Threads: " + stringify(omp_get_max_threads()));
#else
    comm.print("Number of OpenMP Threads: 1 (OpenMP not available)");
#endif
    comm.print("Out-Of-Bounds Value.....: " + String(voxel_map.get_out_of_bounds_value() ? "true" : "false"));
    comm.print("Invert Geometry.........: " + String(invert ? "yes" : "no"));
    comm.print("Create Map From.........: " + String(!off_name.empty() ? "OFF file" : (!formula.empty() ? "formula" : "???")));
    comm.print("Input Geometry File.....: " + (!off_name.empty() ? off_name : String("-N/A-")));
    comm.print("Formula.................: " + (!formula.empty() ? formula : String("-N/A-")));
    comm.print("Input Voxel Map File....: " + (!in_name.empty() ? in_name : String("-N/A-")));
    comm.print("Output Voxel Map File...: " + (!out_name.empty() ? out_name : String("-N/A-")));
    comm.print("Export BMP File Prefix..: " + (!export_name.empty() ? export_name : String("-N/A-")));
    comm.print("Render BMP File Prefix..: " + (!render_name.empty() ? render_name : String("-N/A-")));
    comm.print("Minimum X coord.........: " + stringify_fp_fix(voxel_map.get_bounding_box_min(0), 6, 12));
    comm.print("Maximum X coord.........: " + stringify_fp_fix(voxel_map.get_bounding_box_max(0), 6, 12));
    comm.print("Resolution X............: " + stringify_fp_fix(voxel_map.get_num_points(0) > 0u ?
      (voxel_map.get_bounding_box_max(0) - voxel_map.get_bounding_box_min(0)) / double(voxel_map.get_num_points(0)) : 0.0, 6, 12));
    comm.print("Number of Points X......: " + stringify(voxel_map.get_num_points(0)).pad_front(12));
    comm.print("Minimum Y coord.........: " + stringify_fp_fix(voxel_map.get_bounding_box_min(1), 6, 12));
    comm.print("Maximum Y coord.........: " + stringify_fp_fix(voxel_map.get_bounding_box_max(1), 6, 12));
    comm.print("Resolution Y............: " + stringify_fp_fix(voxel_map.get_num_points(1) > 0u ?
      (voxel_map.get_bounding_box_max(1) - voxel_map.get_bounding_box_min(0)) / double(voxel_map.get_num_points(1)) : 0.0, 6, 12));
    comm.print("Number of Points Y......: " + stringify(voxel_map.get_num_points(1)).pad_front(12));
    comm.print("Minimum Z coord.........: " + stringify_fp_fix(voxel_map.get_bounding_box_min(2), 6, 12));
    comm.print("Maximum Z coord.........: " + stringify_fp_fix(voxel_map.get_bounding_box_max(2), 6, 12));
    comm.print("Resolution Z............: " + stringify_fp_fix(voxel_map.get_num_points(2) ?
      (voxel_map.get_bounding_box_max(2) - voxel_map.get_bounding_box_min(2)) / double(voxel_map.get_num_points(2)): 0.0, 6, 12));
    comm.print("Number of Points Z......: " + stringify(voxel_map.get_num_points(2)).pad_front(12));
    comm.print("Number of Voxels........: " + stringify(voxel_map.get_num_voxels()).pad_front(12));
    comm.print("Bounding Box Volume.....: " + stringify_fp_sci(voxel_map.get_bounding_box_volume(), 3, 12));
    comm.print("Stride Line.............: " + stringify(voxel_map.get_stride_line()).pad_front(12) + " Bytes");
    comm.print("Stride Plane............: " + stringify(voxel_map.get_stride_plane()).pad_front(12) + " Bytes");
    comm.print("Stride Volume...........: " + stringify(voxel_map.get_stride_volume()).pad_front(12) + " Bytes");
    comm.print("Voxel Map Size..........: " + stringify_bytes(voxel_map.get_stride_volume(), 3, 12));

    if(read_result)
    {
      comm.print("Input Map File Size.....: " + stringify_bytes(read_result.filesize, 3, 12));
      comm.print("Compressed Map Size.....: " + stringify_bytes(read_result.compressed_size, 3, 12));
      comm.print("Compression Rate........: " + stringify_fp_fix(100.0*double(read_result.compressed_size) / double(voxel_map.get_stride_volume()), 3, 12) + " %");
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // create the actual map, unless it was read in
    if(!in_name.empty())
    {
      // voxel map was read in from file, so we don't have to create anything here
    }
    else if(!off_name.empty())
    {
      // ok, let's get to work
      comm.print("\nCreating voxel map from OFF file '" + off_name + "'... please have patience...");

      StopWatch watch_create;
      watch_create.start();
      try
      {
        voxel_map.compute_map_from_off_3d(comm, off_name, invert, false); // no gather to all necessary
      }
      catch(const std::exception& exc)
      {
        comm.print(std::cerr, exc.what());
        Runtime::abort();
      }
      watch_create.stop();

      comm.print("Done! Time for voxel map creation " + watch_create.elapsed_string() + " seconds");
    }
    else if(!formula.empty())
    {
      // ok, let's get to work
      comm.print("\nCreating voxel map from formula '" + formula + "'... please have patience...");

      StopWatch watch_create;
      watch_create.start();
      try
      {
        if(dimension == 2)
          voxel_map.compute_map_from_formula_2d(comm, formula, false);
        else
          voxel_map.compute_map_from_formula_3d(comm, formula, false);
      }
      catch(const std::exception& exc)
      {
        comm.print(std::cerr, exc.what());
        Runtime::abort();
      }
      watch_create.stop();

      comm.print("Done! Time for voxel map creation " + watch_create.elapsed_string() + " seconds");
    }

    comm.print("\nDomain Volume...........: " + stringify_fp_sci(voxel_map.get_domain_volume(), 3, 12));
    comm.print("Domain Coverage.........: " + stringify_fp_fix(voxel_map.get_domain_coverage(), 6, 12));

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // write the voxel map to file if desired
    if(!out_name.empty())
    {
      // now let's write the map to file
      comm.print("\nWriting voxel map to file '" + out_name + "'... please have patience...");

      StopWatch watch_write;
      watch_write.start();
      // only rank 0 writes
      Geometry::VoxelMap::WriteResult result;
      if(comm.rank() == 0)
        result = voxel_map.write(out_name, compress_size, compress_level);
      comm.barrier();
      watch_write.stop();
      comm.print("Done! Time for voxel map output " + watch_write.elapsed_string() + " seconds");

      comm.print("\nTotal File Size.........: " + stringify_bytes(result.filesize, 3, 12));
      if(result.compressed_size > 0u)
      {
        comm.print("Compressed Map Size.....: " + stringify_bytes(result.compressed_size, 3, 12));
        comm.print("Compression Rate........: " + stringify_fp_fix(100.0*double(result.compressed_size) / double(voxel_map.get_stride_volume()), 3, 12) + " %");
        comm.print("Compression Level.......: " + stringify(result.compress_level).pad_front(12));
        comm.print("Compression Block Size..: " + stringify(result.compress_block_size).pad_front(12) + " MiB");
        comm.print("Compression Block Count.: " + stringify(result.num_compress).pad_front(12));
      }
      else
      {
        comm.print("Compressed Map Size.....: -N/A-");
        comm.print("Compression Rate........: -N/A-");
        comm.print("Compression Level.......: -N/A-");
        comm.print("Compression Block Size..: -N/A-");
        comm.print("Compression Block Count.: -N/A-");
      }
    }

    // export voxel map to BMP file sequence
    if(!export_name.empty())
    {
      comm.print("\nExporting Voxel Map to BMP files '" + export_name + ".*****.bmp... please have patience...");
      if(comm.rank() == 0)
        voxel_map.export_to_bmp(export_name);
      comm.barrier();
    }

    // render voxel map to BMP file sequence
    if(!render_name.empty())
    {
      comm.print("\nRendering Voxel Map to BMP files '" + render_name + ".*****.bmp... please have patience...");
      if(comm.rank() == 0)
        voxel_map.render_to_bmp(render_name, render_xres, render_yres, render_zres);
      comm.barrier();
    }

    return 0;
  }
} // namespace VoxelMapGenerator

int main(int argc, char** argv)
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);
  return VoxelMapGenerator::main(argc, argv);
}
