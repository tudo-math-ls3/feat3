// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/geometry/voxel_map.hpp>

// includes, thirdparty
#ifdef FEAT_HAVE_ZLIB
#include <zlib.h>
#endif // FEAT_HAVE_ZLIB

namespace FEAT
{
  namespace Geometry
  {
    /// signed 64-bit integer type
    typedef std::int64_t i64;
    /// unsigned 64-bit integer type
    typedef std::uint64_t u64;

    VoxelMap::VoxelMap() :
      _create_stage(0),
      _bbox_min(),
      _bbox_max(),
      _num_points(),
      _stride_line(0u),
      _stride_plane(0u),
      _num_lines(0u),
      _num_planes(0u),
      _voxel_map(),
      _out_of_bounds_value(true),
      _coverage(0u)
    {
    }

    VoxelMap::VoxelMap(VoxelMap&& other) :
      _create_stage(other._create_stage),
      _bbox_min(other._bbox_min),
      _bbox_max(other._bbox_max),
      _num_points(other._num_points),
      _stride_line(other._stride_line),
      _stride_plane(other._stride_plane),
      _num_lines(other._num_lines),
      _num_planes(other._num_planes),
      _voxel_map(std::forward<std::vector<char>>(other._voxel_map)),
      _out_of_bounds_value(other._out_of_bounds_value),
      _coverage(other._coverage)
    {
    }

    /// move assign operator
    VoxelMap& VoxelMap::operator=(VoxelMap&& other)
    {
      if(this == &other)
        return *this;

      _create_stage = other._create_stage;
      _bbox_min = other._bbox_min;
      _bbox_max = other._bbox_max;
      _num_points = other._num_points;
      _stride_line = other._stride_line;
      _stride_plane = other._stride_plane;
      _num_lines = other._num_lines;
      _num_planes = other._num_planes;
      _voxel_map = std::forward<std::vector<char>>(other._voxel_map);
      _out_of_bounds_value = other._out_of_bounds_value;
      _coverage = other._coverage;

      return *this;
    }

    /// virtual destructor
    VoxelMap::~VoxelMap()
    {
    }

    void VoxelMap::set_bounding_box_2d(Real x_min, Real x_max, Real y_min, Real y_max)
    {
      XASSERTM(x_min < x_max, "invalid X dimensions for voxel map bounding box");
      XASSERTM(y_min < y_max, "invalid Y dimensions for voxel map bounding box");
      _bbox_min[0u] = i64(x_min * Real(unit_size));
      _bbox_max[0u] = i64(x_max * Real(unit_size));
      _bbox_min[1u] = i64(y_min * Real(unit_size));
      _bbox_max[1u] = i64(y_max * Real(unit_size));
      _bbox_min[2u] = _bbox_max[2u] = i64(0);
      XASSERTM(_bbox_min[0u] + 1 < _bbox_max[0u], "bounding box is too small in X dimension");
      XASSERTM(_bbox_min[1u] + 1 < _bbox_max[1u], "bounding box is too small in Y dimension");
    }

    void VoxelMap::set_bounding_box_3d(Real x_min, Real x_max, Real y_min, Real y_max, Real z_min, Real z_max)
    {
      XASSERTM(x_min < x_max, "invalid X dimensions for voxel map bounding box");
      XASSERTM(y_min < y_max, "invalid Y dimensions for voxel map bounding box");
      XASSERTM(z_min < z_max, "invalid Z dimensions for voxel map bounding box");
      _bbox_min[0u] = i64(x_min * Real(unit_size));
      _bbox_max[0u] = i64(x_max * Real(unit_size));
      _bbox_min[1u] = i64(y_min * Real(unit_size));
      _bbox_max[1u] = i64(y_max * Real(unit_size));
      _bbox_min[2u] = i64(z_min * Real(unit_size));
      _bbox_max[2u] = i64(z_max * Real(unit_size));
      XASSERTM(_bbox_min[0u] + 1 < _bbox_max[0u], "bounding box is too small in X dimension");
      XASSERTM(_bbox_min[1u] + 1 < _bbox_max[1u], "bounding box is too small in Y dimension");
      XASSERTM(_bbox_min[2u] + 1 < _bbox_max[2u], "bounding box is too small in Z dimension");
    }

    void VoxelMap::set_num_points(Index num_x, Index num_y, Index num_z)
    {
      XASSERT(num_x > 1ull);
      XASSERT(num_y > 1ull);
      XASSERT((num_z > 1ull) || (_bbox_min[2u] == _bbox_max[2u]));
      _num_points[0] = num_x;
      _num_points[1] = num_y;
      _num_points[2] = (num_z > 0u ? num_z : u64(1));

      // compute other quantities
      _stride_line = calc_line_stride(_num_points[0]);
      _stride_plane = _stride_line * _num_points[1];
      _num_planes = _num_points[2];
      _num_lines = _num_points[1] * _num_planes;

      // allocate voxel map
      _voxel_map.resize(_num_planes * _stride_plane, 0u);
    }

    void VoxelMap::set_resolution(Real max_res)
    {
      XASSERTM(max_res > 1E-12, "invalid resolution for voxel map");
      set_num_points(
        u64(_bbox_max[0u] - _bbox_min[0u]) / u64(max_res * Real(unit_size)) + 1u,
        u64(_bbox_max[1u] - _bbox_min[1u]) / u64(max_res * Real(unit_size)) + 1u,
        u64(_bbox_max[2u] - _bbox_min[2u]) / u64(max_res * Real(unit_size)) + 1u);
    }

    void VoxelMap::compute_map_from_formula_2d(const Dist::Comm& comm, const String& formula, bool gather_to_all)
    {
#ifdef FEAT_HAVE_FPARSER
      VoxelFormulaMasker<2> masker(formula);
      this->_compute_voxel_map(comm, masker, gather_to_all);
#else
      XABORTM("FEAT is not build and lined against FPARSER third-party library!");
      (void)comm;
      (void)formula;
      (void)gather_to_all;
#endif
    }

    void VoxelMap::compute_map_from_formula_3d(const Dist::Comm& comm, const String& formula, bool gather_to_all)
    {
#ifdef FEAT_HAVE_FPARSER
      VoxelFormulaMasker<3> masker(formula);
      this->_compute_voxel_map(comm, masker, gather_to_all);
#else
      XABORTM("FEAT is not build and lined against FPARSER third-party library!");
      (void)comm;
      (void)formula;
      (void)gather_to_all;
#endif
    }

    void VoxelMap::compute_map_from_off_3d(const Dist::Comm& comm, const String& filename, bool invert, bool gather_to_all)
    {
      XASSERTM(_bbox_min[2u] < _bbox_max[2u], "CGAL OFF voxel map creation is only available in 3D");
      std::stringstream sstr;
      DistFileIO::read_common(sstr, filename, comm);
      this->compute_map_from_off_3d(comm, sstr, invert, gather_to_all);
    }

    void VoxelMap::compute_map_from_off_3d(const Dist::Comm& comm, std::istream& is, bool invert, bool gather_to_all)
    {
      XASSERTM(_bbox_min[2u] < _bbox_max[2u], "CGAL OFF voxel map creation is only available in 3D");
#ifdef FEAT_HAVE_CGAL
      Geometry::CGALWrapper<double> cgal_wrapper(is, Geometry::CGALFileMode::fm_off);
      VoxelCGALMasker<double> cgal_masker(cgal_wrapper, invert);
      this->_compute_voxel_map(comm, cgal_masker, gather_to_all);
#else
      XABORTM("FEAT is not build and linked against CGAL third-party library!");
      (void)comm;
      (void)is;
      (void)invert;
      (void)gather_to_all;
#endif
    }

    VoxelMap::WriteResult VoxelMap::write(const String& filename, const u64 compress_block_size, const int compress_level) const
    {
      std::ofstream ofs(filename, std::ios_base::out | std::ios_base::binary);
      if(!ofs.is_open() || !ofs.good())
        throw FileNotCreated(filename);

      return write(ofs, compress_block_size, compress_level);
    }

    VoxelMap::WriteResult VoxelMap::write(std::ostream& os, const u64 compress_block_size, const int compress_level) const
    {
      std::size_t written = 0u;

      // set up file header
      FileHeader header;
      memset(&header, 0, sizeof(header));
      header.magic_number = magic_number;
      header.header_size = sizeof(header);
      header.min_x = _bbox_min[0];
      header.max_x = _bbox_max[0];
      header.num_x = _num_points[0];
      header.stride_line = _stride_line;
      header.min_y = _bbox_min[1];
      header.max_y = _bbox_max[1];
      header.num_y = _num_points[1];
      header.stride_plane = _stride_plane;
      header.min_z = _bbox_min[2];
      header.max_z = _bbox_max[2];
      header.num_z = _num_points[2];
      header.stride_volume = _stride_plane * _num_points[2];
      header.coverage = _coverage;

      // _voxel_map may contain trailing padding bytes
      XASSERT(header.stride_volume <= _voxel_map.size());

      // do we have to compress the map?
#ifdef FEAT_HAVE_ZLIB
      if(compress_block_size > 0u)
      {
        XASSERTM(compress_block_size < 1024u, "maximum compression block size is 1024 MB");

        // compute number of planes per compression block
        header.planes_per_block = (compress_block_size*1024ull*1024ull) / _stride_plane;
        XASSERTM(header.planes_per_block > 0ull, "compression block size is to small to hold a single plane!");

        header.num_blocks = (_num_planes + header.planes_per_block - 1u) / header.planes_per_block;
        XASSERTM(header.num_blocks > 0ull, "invalid number of compression blocks!");
      }
#endif // FEAT_HAVE_ZLIB

      // write header
      os.write(reinterpret_cast<char*>(&header), sizeof(header));
      written += sizeof(header);

      // write uncompressed voxel map data
      if(header.num_blocks == 0ull)
      {
        // write uncompressed voxel map data
        os.write(reinterpret_cast<const char*>(_voxel_map.data()), std::streamsize(header.stride_volume));
        return WriteResult(written + header.stride_volume);
      }

#ifdef FEAT_HAVE_ZLIB
      // ensure that a compression block is not larger than 1 GB
      XASSERTM(header.planes_per_block * _stride_plane < 1073741824ull, "voxel map plane compression block size is too big");

      // compression buffers
      std::vector<std::vector<char>> compress_buffers(header.num_blocks);

      // compression block entries
      std::vector<u64> compress_blocks(header.num_blocks, 0u);

      // get compression level
      const int cmp_lvl = Math::max(Math::min(compress_level, 9), 0);

      // accumulated ZLIB results
      int failures = 0;

      // loop over all compression blocks
      FEAT_PRAGMA_OMP(parallel for schedule(dynamic,1) reduction(+:failures))
      for(u64 iblock = 0; iblock < header.num_blocks; ++iblock)
      {
        // compute first and last planes
        u64 first_plane = iblock * header.planes_per_block;
        u64 block_bytes = Math::min(header.planes_per_block, _num_planes - first_plane) * _stride_plane;

        // estimate compression buffer size
        u64 buffer_size = u64(::compressBound(uLong(block_bytes)));

        // allocate compression buffer
        compress_buffers[iblock].resize(buffer_size);

        // compression buffer size
        uLongf dest_size = uLongf(buffer_size);

        // compress via ZLIB
        int rtn = ::compress2(reinterpret_cast<Bytef*>(compress_buffers[iblock].data()), &dest_size,
          reinterpret_cast<const Bytef*>(&_voxel_map[first_plane * _stride_plane]), uLongf(block_bytes), cmp_lvl);

        // failure?
        if(rtn != Z_OK)
        {
          ++failures;
          continue;
        }

        // save compression buffer size
        compress_blocks[iblock] = u64(dest_size);
      }

      // did any compression fail?
      XASSERTM(failures == 0, "ZLIB failed to compress at least one plane block!");

      // write the compression block headers
      os.write(reinterpret_cast<const char*>(compress_blocks.data()), std::streamsize(compress_blocks.size() * sizeof(u64)));
      written += compress_blocks.size() * sizeof(u64);

      // okay, loop over all compression blocks and write them out
      u64 compress_size = 0u;
      for(std::size_t i(0); i < compress_buffers.size(); ++i)
      {
        os.write(reinterpret_cast<const char*>(compress_buffers[i].data()), std::streamsize(compress_blocks[i]));
        written += compress_blocks[i];
        compress_size += compress_blocks[i];
      }

      // return result
      return WriteResult(written, compress_block_size, u64(compress_level), compress_buffers.size(), compress_size);
#else
      // we should never arrive in this #else case even if compiling without zlib
      (void)compress_block_size;
      (void)compress_level;
      XABORTM("INTERNAL ERROR");
#endif // FEAT_HAVE_ZLIB
    }

    VoxelMap::ReadResult VoxelMap::read(const Dist::Comm& comm, const String& filename)
    {
      BinaryStream stream;
      DistFileIO::read_common(stream, filename, comm);
      return read(stream, filename);
    }

    VoxelMap::ReadResult VoxelMap::read(std::istream& is, const String& filename)
    {
      // set up file header
      FileHeader header;
      memset(&header, 0, sizeof(header));

      // try to read header
      if(!is.read(reinterpret_cast<char*>(&header), sizeof(header)).good())
        throw VoxelMapFileError(filename, "Failed to read voxel map header");

      u64 bytes_read = sizeof(header);

      // check magic number
      if(header.magic_number != magic_number)
        throw VoxelMapFileError(filename, "File does not seems to be a voxel map file");

      // check header size
      if(header.header_size != header_size)
        throw VoxelMapFileError(filename, String("invalid header size; expected ") + stringify(header_size) + " but got " + stringify(header.header_size));

      // perform some sanity checks
      if((header.num_x > 0u) && (header.max_x <= header.min_x))
        throw VoxelMapFileError(filename, String("invalid X dimensions: x_max <= x_min"));
      if((header.num_y > 0u) && (header.max_y <= header.min_y))
        throw VoxelMapFileError(filename, String("invalid Y dimensions: y_max <= y_min"));
      if((header.num_z > 0u) && (header.max_z <= header.min_z))
        throw VoxelMapFileError(filename, String("invalid Z dimensions: z_max <= z_min"));
      if((header.num_x == 0u) && (header.max_x != header.min_x))
        throw VoxelMapFileError(filename, String("invalid X dimensions: x_max != x_min for num_x = 0"));
      if((header.num_y == 0u) && (header.max_y != header.min_y))
        throw VoxelMapFileError(filename, String("invalid Y dimensions: y_max != y_min for num_y = 0"));
      if((header.num_z == 0u) && (header.max_z != header.min_z))
        throw VoxelMapFileError(filename, String("invalid Z dimensions: z_max != z_min for num_z = 0"));
      if(header.stride_line*8u < header.num_x)
        throw VoxelMapFileError(filename, String("invalid line stride: too small for a single X-line"));
      if((header.stride_line & 0xF) != 0)
        throw VoxelMapFileError(filename, String("invalid line stride: not a multiple of 16"));
      if(header.stride_plane < header.num_y * header.stride_line)
        throw VoxelMapFileError(filename, String("invalid plane stride: too small for a single XY-plane"));
      if((header.stride_plane & 0xF) != 0)
        throw VoxelMapFileError(filename, String("invalid plane stride: not a multiple of 16"));
      if(header.stride_volume < header.num_z * header.stride_plane)
        throw VoxelMapFileError(filename, String("invalid volume stride: too small for a single XYZ-volume"));
      if((header.stride_volume & 0xF) != 0)
        throw VoxelMapFileError(filename, String("invalid plane stride: not a multiple of 16"));

      // ok, extract the vital data
      _bbox_min[0u] = header.min_x;
      _bbox_max[0u] = header.max_x;
      _bbox_min[1u] = header.min_y;
      _bbox_max[1u] = header.max_y;
      _bbox_min[2u] = header.min_z;
      _bbox_max[2u] = header.max_z;
      _num_points[0u] = header.num_x;
      _num_points[1u] = header.num_y;
      _num_points[2u] = header.num_z;
      _stride_line = header.stride_line;
      _stride_plane = header.stride_plane;
      _num_planes = _num_points[2];
      _num_lines = _num_points[1] * _num_planes;
      _coverage = header.coverage;

      // allocate voxel map
      _voxel_map.resize(header.stride_volume);

      // don't we have compression blocks?
      if(header.num_blocks == 0u)
      {
        // just read the raw map data
        if(!is.read(_voxel_map.data(), std::streamsize(header.stride_volume)).good())
          throw VoxelMapFileError(filename, "Failed to read uncompressed voxel map buffer");
        bytes_read += u64(is.gcount());
        return ReadResult(bytes_read);
      }

      // we're dealing with a ZLIB compressed voxel map

#ifdef FEAT_HAVE_ZLIB
      // read compression blocks
      std::vector<u64> blocks(header.num_blocks);
      if(!is.read(reinterpret_cast<char*>(blocks.data()), std::streamsize(header.num_blocks * sizeof(u64))).good())
        throw VoxelMapFileError(filename, "Failed to read voxel map compression blocks");
      bytes_read += u64(is.gcount());

      // compute buffer offsets and total buffer size
      std::vector<std::size_t> buffer_offset(header.num_blocks + 1u);
      buffer_offset[0u] = 0u;
      for(std::size_t i(0); i < blocks.size(); ++i)
        buffer_offset[i+1u] = buffer_offset[i] + blocks[i];

      // read compression buffer
      std::vector<char> compression_buffer(buffer_offset.back());
      if(!is.read(compression_buffer.data(), std::streamsize(compression_buffer.size())).good())
        throw VoxelMapFileError(filename, "Failed to read compressed voxel map buffer");
      bytes_read += u64(is.gcount());

      // accumulated ZLIB results
      int failures = 0;

      // loop over all compression blocks
      FEAT_PRAGMA_OMP(parallel for schedule(dynamic,1) reduction(+:failures))
      for(u64 iblock = 0; iblock < header.num_blocks; ++iblock)
      {
        u64 first_plane = iblock * header.planes_per_block;
        u64 block_bytes = Math::min(header.planes_per_block, _num_planes - first_plane) * _stride_plane;
        u64 compress_size = blocks[iblock];

        // decompress
        uLongf dest_size = uLongf(block_bytes);
        int rtn = ::uncompress(reinterpret_cast<Bytef*>(&_voxel_map[first_plane * _stride_plane]), &dest_size,
          reinterpret_cast<const Bytef*>(&compression_buffer[buffer_offset[iblock]]), uLong(compress_size));

        if(rtn != Z_OK)
        {
          ++failures;
          continue;
        }
      }

      // did any compression fail?
      XASSERTM(failures == 0, "ZLIB failed to decompress at least one plane block!");
      return ReadResult(bytes_read, header.num_blocks, compression_buffer.size());
#else
      XABORTM("Cannot read compressed voxel map file, because FEAT was compiled without the ZLIB third-party library");
#endif // FEAT_HAVE_ZLIB
    }

    void VoxelMap::export_to_bmp(const String& filename_prefix) const
    {
      XASSERTM(_num_points[2] <= 20000, "voxel map is too big for BMP export!");

      for(u64 iplane(0); iplane < _num_planes; ++iplane)
      {
        String filename = filename_prefix + "." + stringify(iplane).pad_front(5, '0') + ".bmp";
        std::ofstream ofs(filename, std::ios_base::binary);
        _export_plane_to_bmp(ofs, iplane);
        ofs.close();
      }
    }

    void VoxelMap::export_plane_to_bmp(const String& filename, Real z_coord) const
    {
      std::ofstream ofs(filename, std::ios_base::binary);
      export_plane_to_bmp(ofs, z_coord);
      ofs.close();
    }

    void VoxelMap::export_plane_to_bmp(std::ostream& os, Real z_coord) const
    {
      _export_plane_to_bmp(os, _map_coord_idx_nearest(i64(z_coord * Real(unit_size)), 2u));
    }

    void VoxelMap::render_to_bmp(const String& filename_prefix, Index width, Index height, Index depth) const
    {
      XASSERTM(depth <= 20000, "render width is too big for BMP export!");
      for(i64 iplane(0); iplane < i64(depth); ++iplane)
      {
        String filename = filename_prefix + "." + stringify(iplane).pad_front(5, '0') + ".bmp";
        std::ofstream ofs(filename, std::ios_base::binary);
        i64 z_min = _bbox_min[2] + ( iplane    * (_bbox_max[2] - _bbox_min[1])) / i64(depth);
        i64 z_max = _bbox_min[2] + ((iplane+1) * (_bbox_max[2] - _bbox_min[1])) / i64(depth);
        _render_plane_to_bmp(ofs, width, height, z_min, z_max);
        ofs.close();
      }
    }

    void VoxelMap::render_plane_to_bmp(const String& filename, Index width, Index height, Real z_min, Real z_max) const
    {
      std::ofstream ofs(filename, std::ios_base::binary);
      render_plane_to_bmp(ofs, width, height, z_min, z_max);
      ofs.close();
    }

    void VoxelMap::render_plane_to_bmp(std::ostream& os, Index width, Index height, Real z_min, Real z_max) const
    {
      _render_plane_to_bmp(os, width, height, i64(z_min * Real(unit_size)), i64(z_max * Real(unit_size)));
    }

    Real VoxelMap::get_bounding_box_volume() const
    {
      if(_num_points[2] > 0u)
        return Real(_bbox_max[2] - _bbox_min[2]) * Real(_bbox_max[1] - _bbox_min[1]) * Real(_bbox_max[0] - _bbox_min[0]) / Math::cub(Real(unit_size));
      if(_num_points[1] > 0u)
        return Real(_bbox_max[1] - _bbox_min[1]) * Real(_bbox_max[0] - _bbox_min[0]) / Math::sqr(Real(unit_size));
      if(_num_points[0] > 0u)
        return Real(_bbox_max[0] - _bbox_min[0]) / Real(unit_size);
      return Real(0);
    }

    u64 VoxelMap::_compute_domain_coverage() const
    {
      if(_voxel_map.empty())
        return u64(0);

      u64 count(0u), n(_voxel_map.size()), nv(get_num_voxels());
      FEAT_PRAGMA_OMP(parallel for reduction(+:count))
      for(u64 i = 0u; i < n; ++i)
      {
        for(int k = 0; k < 8; ++k)
          count += u64((_voxel_map[i] >> k) & 1);
      }

      // convert coverage relative to unit size
      return (count*unit_size) / nv;
    }

    u64 VoxelMap::_map_coord_idx_nearest(i64 xyz, std::size_t sdim) const
    {
      // The input coordinate X is expected to be in the range of the bounding box [BMIN,BMAX], which is discretized
      // by the NP voxel coordinates in the range {0, 1, ..., NP-1}, so applying standard linear transformation from
      // the interval [BMIN, BMAX] to [0, NP-1] yields the formula
      //
      //                           (NP - 1) * (X - BMIN)
      //                 X -> 0 + -----------------------
      //                               BMAX - BMIN
      //
      // However, we are working with integer arithmetic here, so the fraction will be handled by integer division,
      // which truncates all decimal places, so we have to add 1/2 to the result of the fraction to achieve a
      // "round-to-nearest" voxel behavior. We exploit that 1/2 = (BMAX - BMIN) / (2*BMAX - 2*MIN) and so we get
      //
      //                        (NP - 1) * (X - BMIN)         BMAX - BMIN
      //                 X ->  -----------------------  +  -----------------
      //                            BMAX - BMIN             2*BMAX - 2*BMIN
      //
      // which is equivalent to
      //
      //                       2 * (NP - 1) * (X - BMIN) + BMAX - BMIN
      //                 X -> -----------------------------------------
      //                               2 * BMAX - 2 * BMIN
      //
      i64 idx = (2*(i64(_num_points[sdim]) - 1) * (xyz - _bbox_min[sdim]) + _bbox_max[sdim] - _bbox_min[sdim])  / (2*_bbox_max[sdim] - 2*_bbox_min[sdim]);

      // also check whether the index is out of bounds
      return (0 <= idx) && (idx < i64(this->_num_points[sdim])) ? u64(idx) : ~u64(0);
    }

    i64 VoxelMap::_map_coord_idx_lower(i64 xyz, std::size_t sdim) const
    {
      // The input coordinate X is expected to be in the range of the bounding box [BMIN,BMAX], which is discretized
      // by the NP voxel coordinates in the range {0, 1, ..., NP-1}, so applying standard linear transformation from
      // the interval [BMIN, BMAX] to [0, NP-1] yields the formula
      //
      //                       (NP - 1) * (X - BMIN)
      //                 X -> ----------------------- =: IDX
      //                           BMAX - BMIN
      //
      // We are working with integer arithmetic here, so the fraction will be handled by integer division,
      // which truncates all decimal places, which is what we want in this case to obtain the lower voxel index.
      return ((i64(_num_points[sdim]) - 1) * (xyz - _bbox_min[sdim]))  / (_bbox_max[sdim] - _bbox_min[sdim]);
    }

    i64 VoxelMap::_map_coord_idx_upper(i64 xyz, std::size_t sdim) const
    {
      // The input coordinate X is expected to be in the range of the bounding box [BMIN,BMAX], which is discretized
      // by the NP voxel coordinates in the range {0, 1, ..., NP-1}, so applying standard linear transformation from
      // the interval [BMIN, BMAX] to [0, NP-1] yields the formula
      //
      //                       (NP - 1) * (X - BMIN)
      //                 X -> -----------------------
      //                           BMAX - BMIN
      //
      // We are working with integer arithmetic here, so the fraction will be handled by integer division,
      // which truncates all decimal places, so to round an integer fraction A/B up to the next integer, we have to
      // add (B-1) to the numerator to obtain (A+B-1)/B, which in our case results in the formula
      //
      //                       (NP - 1) * (X - BMIN) + BMAX - BMIN - 1
      //                 X -> ----------------------------------------- =: IDX
      //                                    BMAX - BMIN
      //
      return ((i64(_num_points[sdim]) - 1) * (xyz - _bbox_min[sdim]) + _bbox_max[sdim] - _bbox_min[sdim] - 1) / (_bbox_max[sdim] - _bbox_min[sdim]);
    }

    Real VoxelMap::_calc_sample_rate(i64 xyz, i64 idx, std::size_t sdim) const
    {
      // the sample rate is defined as   sr(xyz,idx) :=  1 - |xyz - V(idx)| / |V(idx+1) - V(idx)|  where
      //
      //        V(idx+1) - V(idx) = (V(n) - V(0)) / (n-1)
      // and
      //                   V(idx) = V(0) + idx*(V(idx+1) - V(idx)) = V(0) + idx*(V(n) - V(0)) / (n-1)
      //
      // so in total we have
      //
      // sr(xyz,idx) =  1 - |xyz - V(idx)| / |V(idx+1) - V(idx)|
      //             =  1 - |(xyz - V(idx)) / (V(idx+1) - V(idx))|
      //             =  1 - |(xyz - V(0) - idx*(V(n) - V(0)) / (n-1)) / ((V(n) - V(0)) / (n-1))|
      //             =  1 - |((xyz - V(0))*(n-1) - idx*(V(n) - V(0))) / (V(n) - V(0))|
      //             =  1 - |(xyz - V(0))*(n-1) / (V(n) - V(0)) - idx)|
      //
      return Real(1) - Math::abs(Real(_num_points[sdim] - 1) * Real(xyz - _bbox_min[sdim]) / Real(_bbox_max[sdim] - _bbox_min[sdim]) - Real(idx));
    }

    bool VoxelMap::_check_point(i64 xidx, i64 yidx, i64 zidx) const
    {
      if((xidx < i64(0)) || (i64(_num_points[0]) <= xidx))
        return Real(_out_of_bounds_value);
      if((yidx < i64(0)) || (i64(_num_points[1]) <= yidx))
        return Real(_out_of_bounds_value);
      if((zidx < i64(0)) || (i64(_num_points[2]) <= zidx))
        return Real(_out_of_bounds_value);
      return ((this->_voxel_map[zidx*this->_stride_plane + yidx*this->_stride_line + (xidx>>3)] >> (xidx & 7)) & 1u) != 0u;
    }

    bool VoxelMap::_check_point_nearest(const std::array<i64, 1>& p) const
    {
      u64 xidx = _map_coord_idx_nearest(p[0], 0u);
      if(xidx == ~u64(0))
        return _out_of_bounds_value;
      return (this->_voxel_map[(xidx >> 3)] >> (xidx & 0x7)) & 0x1;
    }

    bool VoxelMap::_check_point_nearest(const std::array<i64, 2>& p) const
    {
      u64 xidx = _map_coord_idx_nearest(p[0], 0u);
      u64 yidx = _map_coord_idx_nearest(p[1], 1u);
      if((xidx|yidx) == ~u64(0))
        return _out_of_bounds_value;
      return (this->_voxel_map[yidx * this->_stride_line + (xidx >> 3)] >> (xidx & 0x7)) & 0x1;
    }

    bool VoxelMap::_check_point_nearest(const std::array<i64, 3>& p) const
    {
      u64 xidx = _map_coord_idx_nearest(p[0], 0u);
      u64 yidx = _map_coord_idx_nearest(p[1], 1u);
      u64 zidx = _map_coord_idx_nearest(p[2], 2u);
      if((xidx|yidx|zidx) == ~u64(0))
        return _out_of_bounds_value;
      u64 line = zidx * this->_num_points[1] + yidx;
      return (this->_voxel_map[line * this->_stride_line + (xidx >> 3)] >> (xidx & 0x7)) & 0x1;
    }

    bool VoxelMap::_check_box(const std::array<u64, 1>& box_min, const std::array<u64, 1>& box_max) const
    {
      for(u64 xidx(box_min[0]); xidx <= box_max[0]; ++xidx)
      {
        if((this->_voxel_map[(xidx >> 3)] >> (xidx & 0x7)) & 0x1)
          return true;
      }
      return false;
    }

    bool VoxelMap::_check_box(const std::array<u64, 2>& box_min, const std::array<u64, 2>& box_max) const
    {
      for(u64 yidx(box_min[1]); yidx <= box_max[1]; ++yidx)
      {
        for(u64 xidx(box_min[0]); xidx <= box_max[0]; ++xidx)
        {
          if((this->_voxel_map[yidx * this->_stride_line + (xidx >> 3)] >> (xidx & 0x7)) & 0x1)
            return true;
        }
      }
      return false;
    }

    bool VoxelMap::_check_box(const std::array<u64, 3>& box_min, const std::array<u64, 3>& box_max) const
    {
      for(u64 zidx(box_min[2]); zidx <= box_max[2]; ++zidx)
      {
        for(u64 yidx(box_min[1]); yidx <= box_max[1]; ++yidx)
        {
          u64 line = zidx * this->_num_points[1] + yidx;
          for(u64 xidx(box_min[0]); xidx <= box_max[0]; ++xidx)
          {
            if((this->_voxel_map[line * this->_stride_line + (xidx >> 3)] >> (xidx & 0x7)) & 0x1)
              return true;
          }
        }
      }
      return false;
    }

    Real VoxelMap::_sample_point_1d_x(i64 px, i64 xidx, i64 yidx, i64 zidx) const
    {
      Real v0 = Real(check_point(xidx, yidx, zidx));
      Real v1 = Real(check_point(xidx+1, yidx, zidx));
      return  v0 + (v1 - v0) * _calc_sample_rate(px, xidx, 0u);
    }

    Real VoxelMap::_sample_point_1d_y(i64 py, i64 xidx, i64 yidx, i64 zidx) const
    {
      Real v0 = Real(check_point(xidx, yidx, zidx));
      Real v1 = Real(check_point(xidx, yidx+1, zidx));
      return  v0 + (v1 - v0) * _calc_sample_rate(py, yidx, 1u);
    }

    Real VoxelMap::_sample_point_1d_z(i64 pz, i64 xidx, i64 yidx, i64 zidx) const
    {
      Real v0 = Real(check_point(xidx, yidx, zidx));
      Real v1 = Real(check_point(xidx, yidx, zidx+1));
      return  v0 + (v1 - v0) * _calc_sample_rate(pz, zidx, 2u);
    }

    Real VoxelMap::_sample_point_2d(i64 px, i64 py, i64 xidx, i64 yidx, i64 zidx) const
    {
      Real v0 = _sample_point_1d_x(px, xidx, yidx, zidx);
      Real v1 = _sample_point_1d_x(px, xidx, yidx+1, zidx);
      return  v0 + (v1 - v0) * _calc_sample_rate(py, yidx, 1u);
    }

    Real VoxelMap::_sample_point_3d(i64 px, i64 py, i64 pz, i64 xidx, i64 yidx, i64 zidx) const
    {
      Real v0 = _sample_point_2d(px, py, xidx, yidx, zidx);
      Real v1 = _sample_point_2d(px, py, xidx, yidx, zidx+1);
      return  v0 + (v1 - v0) * _calc_sample_rate(pz, zidx, 2u);
    }

    Real VoxelMap::_sample_point(const std::array<i64, 1>& p) const
    {
      return _sample_point_1d_x(p[0], _map_coord_idx_lower(p[0], 0u));
    }

    Real VoxelMap::_sample_point(const std::array<i64, 2>& p) const
    {
      return _sample_point_2d(p[0], p[1], _map_coord_idx_lower(p[0], 0u), _map_coord_idx_lower(p[1], 1u));
    }

    Real VoxelMap::_sample_point(const std::array<i64, 3>& p) const
    {
      return _sample_point_3d(p[0], p[1], p[2], _map_coord_idx_lower(p[0], 0u),
        _map_coord_idx_lower(p[1], 1u), _map_coord_idx_lower(p[2], 2u));
    }

    Real VoxelMap::_sample_box_1d(i64 px0, i64 px1, i64 xidx0, i64 xidx1, i64 yidx, i64 zidx) const
    {
      ASSERT(xidx0 < xidx1);

      // is the box completely left or completely right of the voxel map domain?
      if(xidx1 < i64(0))
        return Real(_out_of_bounds_value);
      if(i64(_num_points[0]) <= xidx0)
        return Real(_out_of_bounds_value);

      // sample the two voxel map entries left and right of px0
      Real v00 = _check_point(xidx0, yidx, zidx);
      Real v01 = _check_point(xidx0+1, yidx, zidx);

      // sample the two voxel map entries left and right of px1
      Real v10 = _check_point(xidx1-1, yidx, zidx);
      Real v11 = _check_point(xidx1, yidx, zidx);

      // compute sample rate of px0 w.r.t. left voxel V(xidx0) and of px1 w.r.t. right voxel V(xidx1)
      Real s0 = _calc_sample_rate(px0, xidx0, 0);
      Real s1 = _calc_sample_rate(px1, xidx1, 0);

      // compute interpolated values v0 and v1 of px0 and px1, respectively, from their neighboring voxel values
      Real v0 = v01 + (v00 - v01) * s0;
      Real v1 = v10 + (v11 - v10) * s1;

      // if the entire box [px0, px1] is between two neighboring voxels, we can integrate the interpolated voxel
      // values v0 and v1 at px0 and px1 via the trapezoidal rule by simply computing the average of the two values:
      if(xidx1 == xidx0 + 1)
        return Real(0.5) * (v0 + v1);

      // in any other case, we have to integrate the voxel values over the interval [px0,px1], where v0 and v1
      // denote the previously interpolated voxel values at the interval endpoints px0 and px1, respectively:
      //
      //                     v0                                                      v1
      //            --V------+======V=============V=====...=====V=============V======+------V
      //            xidx0   px0  xidx0+1         ...           ...         xidx1-1  px1  xidx1
      //                     |--s0--|                                         |--s1--|
      //
      // we compute the integral in three steps and then compute the sum of those three sub-integrals:
      // 1. compute the integral over the left-most sub-interval [px0,V(xidx0+1)] via the trapezoidal rule
      // 2. compute the integral over the right-most sub-interval [V(xidx1-1),px1] via the trapezoidal rule
      // 3. compute the integral over the inner intervals [V(idx0+1),V(idx1-1)] via the summed trapezoidal rule
      //
      // for the outer two sub-interval we have that the length of...
      // ...the left-most sub-interval [px0,V(xidx0+1)] is given by the sample rate s0 of px0 w.r.t. V(idx0)
      // ...the right-most sub-interval [V(xidx1-1),px1] is given by the sample rate s1 of px1 w.r.t. V(idx1)
      //
      // therefore, the integral of
      /// ...the left-most interval equals to s0*(v0+V(idx0+1))/2
      /// ...the right-most interval equals to s1*(v1+V(idx1-1))/2
      //  ...the inner interval equals to V(idx0+1)/2 + V(idx1-1)/2 + sum_{xidx0+1<=k<=idx1-1} V(k)
      //
      // which results in
      //
      //   s0*(v0+V(idx0+1))/2 + V(idx0+1)/2 + sum_{xidx0+1 < k < idx1-1} V(k) + V(idx1-1)/2 + s1*(v1+V(idx1-1))/2
      //  \___________________/ \___________________________________________________________/ \___________________/
      //   left-most interval                          inner intervals                         right-most interval
      //
      // which in the case where idx0+1 < idx1-1 is equal to
      //
      //  (s0*v0)/2 + (s0+1)*V(idx0+1)/2 + sum_{xidx0+1<k<idx1-1} V(k) + (s1+1)*V(idx1-1)/2 + (s1*v1)/2
      //
      // there is one final special case where idx0+1 = idx1-1, i.e. there is just a single voxel between the
      // two points px0 and px1 -- or in other words: there are no inner intervals -- and therefore we just have
      // to add up the two integrals over the outer most sub-intervals, so for this case we get
      //
      //      s0*(v0+V(idx0+1))/2 + s1*(v1+V(idx1-1))/2 = (s0*v0 + s1*v1 + (s0+s1)*V(idx0+1))/2

      // add scaled endpoint values for trapezoidal rule of first and last intervals of length s0 and s1
      Real result = Real(0.5) * (s0*v0 + s1*v1);

      // add first and last inner voxel and scale properly to compensate for the following summed trapezoidal rule
      if(xidx1 == xidx0 + 2)
      {
        // no inner X-intervals
        result += Real(0.5) * (s0 + s1) * _check_point(xidx0+1, yidx, zidx);
      }
      else
      {
        // at least one inner X-interval exists
        result += Real(0.5) * (s0 + Real(1)) * _check_point(xidx0+1, yidx, zidx);
        result += Real(0.5) * (s1 + Real(1)) * _check_point(xidx1-1, yidx, zidx);
      }

      // now add all inner voxel values for summed trapezoidal rule
      for(i64 i = xidx0+2; i < xidx1-1; ++i)
        result += _check_point(i, yidx, zidx);

      // finally divide by the total integration interval width
      return result / (Real(xidx1 - xidx0 - 2) + s0 + s1);
    }

    Real VoxelMap::_sample_box_2d(i64 px0, i64 px1, i64 py0, i64 py1, i64 xidx0, i64 xidx1, i64 yidx0, i64 yidx1, i64 zidx) const
    {
      ASSERT(yidx0 < yidx1);

      // is the box completely below or completely above of the voxel map domain?
      if(yidx1 < i64(0))
        return Real(_out_of_bounds_value);
      if(i64(_num_points[1]) <= yidx0)
        return Real(_out_of_bounds_value);

      // sample the two X-parallel voxel map lines below and above of py0
      Real v00 = _sample_box_1d(px0, px1, xidx0, xidx1, yidx0, zidx);
      Real v01 = _sample_box_1d(px0, px1, xidx0, xidx1, yidx0+1, zidx);

      // sample the two X-parallel voxel map lines below and above of py1
      Real v10 = _sample_box_1d(px0, px1, xidx0, xidx1, yidx1-1, zidx);
      Real v11 = _sample_box_1d(px0, px1, xidx0, xidx1, yidx1, zidx);

      // compute sample rate of py0 w.r.t. lower voxel V(yidx0) and of py1 w.r.t. upper voxel V(yidx1)
      Real s0 = _calc_sample_rate(py0, yidx0, 1);
      Real s1 = _calc_sample_rate(py1, yidx1, 1);

      // compute interpolated values v0 and v1 of py0 and py1, respectively, from their neighboring voxel values
      Real v0 = v01 + (v00 - v01) * s0;
      Real v1 = v10 + (v11 - v10) * s1;

      // if the entire box [py0, py1] is between two neighboring voxel lines, we can integrate the interpolated voxel
      // line values v0 and v1 at py0 and py1 via the trapezoidal rule by simply computing the average of the two values:
      if(yidx1 == yidx0 + 1)
        return Real(0.5) * (v0 + v1);

      // see the comments inside _sample_box_1d() for more details about what is happening here

      // add scaled endpoint values for trapezoidal rule of first and last intervals of length s0 and s1
      Real result = Real(0.5) * (s0*v0 + s1*v1);

      // add first and last inner voxel rows and scale properly to compensate for the following summed trapezoidal rule
      if(yidx1 == yidx0 + 2)
      {
        // no inner Y-intervals
        result += Real(0.5) * (s0 + s1) * _sample_box_1d(px0, px1, xidx0, xidx1, yidx0+1, zidx);
      }
      else
      {
        // at least one inner Y-interval exists
        result += Real(0.5) * (s0 + Real(1)) * _sample_box_1d(px0, px1, xidx0, xidx1, yidx0+1, zidx);
        result += Real(0.5) * (s1 + Real(1)) * _sample_box_1d(px0, px1, xidx0, xidx1, yidx1-1, zidx);
      }

      // now add all inner X-parallel voxel line values for summed trapezoidal rule
      for(i64 i = yidx0+2; i < yidx1-1; ++i)
        result += _sample_box_1d(px0, px1, xidx0, xidx1, i, zidx);

      // finally divide by the total integration interval height
      return result / (Real(yidx1 - yidx0 - 2) + s0 +s1);
    }

    Real VoxelMap::_sample_box_3d(i64 px0, i64 px1, i64 py0, i64 py1, i64 pz0, i64 pz1, i64 xidx0, i64 xidx1, i64 yidx0, i64 yidx1, i64 zidx0, i64 zidx1) const
    {
      ASSERT(zidx0 < zidx1);

      // is the box completely below or completely above of the voxel map domain?
      if(zidx1 < i64(0))
        return Real(_out_of_bounds_value);
      if(i64(_num_points[2]) <= zidx0)
        return Real(_out_of_bounds_value);

      // sample the two XY-parallel voxel map planes in front and behind of pz0
      Real v00 = _sample_box_2d(px0, px1, py0, py1, xidx0, xidx1, yidx0, yidx1, zidx0);
      Real v01 = _sample_box_2d(px0, px1, py0, py1, xidx0, xidx1, yidx0, yidx1, zidx0+1);

      // sample the two XY-parallel voxel map planes in front and behind of pz1
      Real v10 = _sample_box_2d(px0, px1, py0, py1, xidx0, xidx1, yidx0, yidx1, zidx1-1);
      Real v11 = _sample_box_2d(px0, px1, py0, py1, xidx0, xidx1, yidx0, yidx1, zidx1);

      // compute sample rate of pz0 w.r.t. back voxel V(zidx0) and of pz1 w.r.t. front voxel V(zidx1)
      Real s0 = _calc_sample_rate(pz0, zidx0, 2);
      Real s1 = _calc_sample_rate(pz1, zidx1, 2);

      // compute interpolated values v0 and v1 of pz0 and pz1, respectively, from their neighboring voxel values
      Real v0 = v01 + (v00 - v01) * s0;
      Real v1 = v10 + (v11 - v10) * s1;

      // if the entire box [pz0, pz1] is between two neighboring voxel planes, we can integrate the interpolated voxel
      // plane values v0 and v1 at pz0 and pz1 via the trapezoidal rule by simply computing the average of the two values:
      if(zidx1 == zidx0 + 1)
        return Real(0.5) * (v0 + v1);

      // see the comments inside _sample_box_1d() for more details about what is happening here

      // add scaled endpoint values for trapezoidal rule of first and last intervals of length s0 and s1
      Real result = Real(0.5) * (s0*v0 + s1*v1);

      // add first and last inner voxel planes and scale properly to compensate for the following summed trapezoidal rule
      if(zidx1 == zidx0 + 2)
      {
        // no inner Z-intervals
        result += Real(0.5) * (s0 + s1) * _sample_box_2d(px0, px1, py0, py1, xidx0, xidx1, yidx0, yidx1, zidx0+1);
      }
      else
      {
        // at least one inner Z-interval exists
        result += Real(0.5) * (s0 + Real(1)) * _sample_box_2d(px0, px1, py0, py1, xidx0, xidx1, yidx0, yidx1, zidx0+1);
        result += Real(0.5) * (s1 + Real(1)) * _sample_box_2d(px0, px1, py0, py1, xidx0, xidx1, yidx0, yidx1, zidx1-1);
      }

      // now add all inner XY-parallel voxel plane values for summed trapezoidal rule
      for(i64 i = zidx0+2; i < zidx1-1; ++i)
        result += _sample_box_2d(px0, px1, py0, py1, xidx0, xidx1, yidx0, yidx1, i);

      // finally divide by the total integration interval depth
      return result / (Real(zidx1 - zidx0 - 2) + s0 +s1);
    }

    Real VoxelMap::_sample_box(const std::array<i64, 1>& box_min, const std::array<i64, 1>& box_max) const
    {
      if(_num_points[0] > 1u)
        return _sample_box_1d(box_min[0], box_max[0],
          _map_coord_idx_lower(box_min[0], 0u), _map_coord_idx_upper(box_max[0], 0u));
      else
        return Real(0);
    }

    Real VoxelMap::_sample_box(const std::array<i64, 2>& box_min, const std::array<i64, 2>& box_max) const
    {
      if(_num_points[1] > 1u)
        return _sample_box_2d(box_min[0], box_max[0], box_min[1], box_max[1],
        _map_coord_idx_lower(box_min[0], 0u), _map_coord_idx_upper(box_max[0], 0u),
        _map_coord_idx_lower(box_min[1], 1u), _map_coord_idx_upper(box_max[1], 1u));
      else if(_num_points[0] > 1u)
        return _sample_box_1d(box_min[0], box_max[0],
        _map_coord_idx_lower(box_min[0], 0u), _map_coord_idx_upper(box_max[0], 0u));
      else
        return Real(0);
    }

    Real VoxelMap::_sample_box(const std::array<i64, 3>& box_min, const std::array<i64, 3>& box_max) const
    {
      if(_num_points[2] > 1u)
        return _sample_box_3d(box_min[0], box_max[0], box_min[1], box_max[1], box_min[2], box_max[2],
          _map_coord_idx_lower(box_min[0], 0u), _map_coord_idx_upper(box_max[0], 0u),
          _map_coord_idx_lower(box_min[1], 1u), _map_coord_idx_upper(box_max[1], 1u),
          _map_coord_idx_lower(box_min[2], 2u), _map_coord_idx_upper(box_max[2], 2u));
      else if(_num_points[1] > 1u)
        return _sample_box_2d(box_min[0], box_max[0], box_min[1], box_max[1],
          _map_coord_idx_lower(box_min[0], 0u), _map_coord_idx_upper(box_max[0], 0u),
          _map_coord_idx_lower(box_min[1], 1u), _map_coord_idx_upper(box_max[1], 1u));
      else if(_num_points[0] > 1u)
        return _sample_box_1d(box_min[0], box_max[0],
          _map_coord_idx_lower(box_min[0], 0u), _map_coord_idx_upper(box_max[0], 0u));
      else
        return Real(0);
    }

    void VoxelMap::_compress_voxel_map_line(const std::vector<int>& mask, const u64 line)
    {
      ASSERTM(line < this->_num_lines, "invalid voxel map line");
      const u64 n = u64(mask.size());
      const u64 off = line * this->_stride_line;
      for(u64 i(0); i < n; ++i)
      {
        _voxel_map[off + (i >> 3)] |= char(mask[i] != 0) << (i & 0x7);
      }
    }

    void VoxelMap::_export_plane_to_bmp(std::ostream& os, u64 plane) const
    {
      // maximum BMP resolutions seems to be around 30K, but we won't export more than 20K,
      // which results in a file size of roughly 50 MB per image
      XASSERTM(_num_points[0] <= 20000, "voxel map is too big for BMP export!");
      XASSERTM(_num_points[1] <= 20000, "voxel map is too big for BMP export!");
      XASSERTM(plane < _num_points[2], "invalid plane index");

      typedef std::uint16_t u16;
      typedef std::uint32_t u32;

      // get dimensions
      u32 width = u32(this->_num_points[0]);
      u32 height = u32(this->_num_points[1]);

      // compute stride (multiple of 4 bytes)
      u32 stride = ((width + 31u) & ~31u) >> 3;

      // BITMAPFILEHEADER: 14 bytes
      // write magic
      u16 magic = 0x4D42;
      os.write((char*)&magic, 2u);
      // write file size
      u32 filesize = 62u + height * stride;
      os.write((char*)&filesize, 4u);
      // write reserved bytes
      u32 zeros = 0u;
      os.write((char*)&zeros, 4u);
      // wite bit offset
      u32 offbits = 54u + 2u*4u;
      os.write((char*)&offbits, 4u);

      // BITMAPINFOHEADER: 40 bytes
      // write header size
      u32 bisize = 40u;
      os.write((char*)&bisize, 4u);
      // write width and height
      os.write((char*)&width, 4u);
      os.write((char*)&height, 4u);
      // write plane count
      u16 planes = 1;
      os.write((char*)&planes, 2u);
      // write bit count
      u16 bitcount = 1;
      os.write((char*)&bitcount, 2u);
      // write compression
      os.write((char*)&zeros, 4u);
      // write image size
      u32 size_image = height * stride;
      os.write((char*)&size_image, 4u);
      // write pixels per meter (X and Y)
      u32 px_per_m = 3780;
      os.write((char*)&px_per_m, 4u);
      os.write((char*)&px_per_m, 4u);
      // number of colors
      u32 clr_used = 2u;
      os.write((char*)&clr_used, 4u);
      // number of important colors (WTF?)
      os.write((char*)&zeros, 4u);

      // write black and white colors
      u32 clr_black = 0u, clr_white = ~u32(0);
      os.write((char*)&clr_black, 4u);
      os.write((char*)&clr_white, 4u);

      std::vector<char> linebuffer(stride);
      char* oline = linebuffer.data();

      // okay, let's loop over all domain nodes
      for(u32 i(0); i < height; ++i)
      {
        const char* iline = &_voxel_map.data()[plane*_stride_plane + i*_stride_line];
        for(u32 j(0); j < stride; ++j)
          for(int k(0); k < 8; ++k)
            (oline[j] <<= 1) |= ((iline[j] >> k) & 1); // all hail to bitshift!
        os.write(oline, stride);
      }
    }

    void VoxelMap::_render_plane_to_bmp(std::ostream& os, Index width, Index height, i64 box_min_z, i64 box_max_z) const
    {
      // maximum BMP resolutions seems to be around 30K, but we won't export more than 20K,
      // which results in a file size of roughly 50 MB per image
      XASSERTM(width <= 20000, "render width is too big for BMP export!");
      XASSERTM(height <= 20000, "render width is too big for BMP export!");
      XASSERTM(box_min_z < box_max_z, "invalid depth for voxel map render!");
      XASSERTM(box_min_z >= _bbox_min[2], "invalid depth for voxel map render!");
      XASSERTM(box_max_z <= _bbox_max[2], "invalid depth for voxel map render!");

      typedef std::uint8_t u8;
      typedef std::uint16_t u16;
      typedef std::uint32_t u32;

      // get dimensions
      u32 width2 = u32(width);
      u32 height2 = u32(height);

      // compute stride (multiple of 4 bytes)
      u32 stride = ((width2*8 + 31u) & ~31u) >> 3;

      // BITMAPFILEHEADER: 14 bytes
      // write magic
      u16 magic = 0x4D42;
      os.write((char*)&magic, 2u);
      // write file size
      u32 filesize = 54u + height2 * stride + 256u*4u;
      os.write((char*)&filesize, 4u);
      // write reserved bytes
      u32 zeros = 0u;
      os.write((char*)&zeros, 4u);
      // wite bit offset
      u32 offbits = 54u + 256u*4u;
      os.write((char*)&offbits, 4u);

      // BITMAPINFOHEADER: 40 bytes
      // write header size
      u32 bisize = 40u;
      os.write((char*)&bisize, 4u);
      // write width and height
      os.write((char*)&width2, 4u);
      os.write((char*)&height2, 4u);
      // write plane count
      u16 planes = 1;
      os.write((char*)&planes, 2u);
      // write bit count
      u16 bitcount = 8;
      os.write((char*)&bitcount, 2u);
      // write compression
      os.write((char*)&zeros, 4u);
      // write image size
      u32 size_image = height2 * stride;
      os.write((char*)&size_image, 4u);
      // write pixels per meter (X and Y)
      u32 px_per_m = 3780;
      os.write((char*)&px_per_m, 4u);
      os.write((char*)&px_per_m, 4u);
      // number of colors
      u32 clr_used = 256u;
      os.write((char*)&clr_used, 4u);
      // number of important colors (WTF?)
      os.write((char*)&zeros, 4u);

      // write colors
      std::vector<u32> colors(256);
      for(u32 i = 0u; i < 256u; ++i)
      {
        colors[i] = i | (i << 8) | (i << 16);
      }
      os.write((char*)colors.data(), 256u*4u);

      std::vector<u8> linebuffer(stride);
      u8* oline = linebuffer.data();

      // set up bounding box
      std::array<i64, 3u> box_min, box_max;
      box_min[2u] = box_min_z;
      box_max[2u] = box_max_z;

      // okay, let's loop over all domain nodes
      for(Index i(0); i < height; ++i)
      {
        box_min[1u] = i64(_bbox_min[1] + (i * (_bbox_max[1] - _bbox_min[1])) / height);
        box_max[1u] = i64(_bbox_min[1] + ((i+1) * (_bbox_max[1] - _bbox_min[1])) / height);
        for(Index j(0); j < width; ++j)
        {
          box_min[0u] = i64(_bbox_min[0] + (j * (_bbox_max[0] - _bbox_min[0])) / width);
          box_max[0u] = i64(_bbox_min[0] + ((j+1) * (_bbox_max[0] - _bbox_min[0])) / width);
          Real w = _sample_box(box_min, box_max) * 255.0;
          oline[j] = (w <= 0.0 ? 0 : (w >= 255.0 ? 255 : u8(w)));
        }
        os.write((char*)oline, stride);
      }
    }
  } // namespace Geometry
} // namespace FEAT
