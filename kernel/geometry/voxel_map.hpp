// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/dist_file_io.hpp>
#include <kernel/geometry/cgal.hpp>
#include <kernel/geometry/atlas/chart.hpp>

// includes, system
#include <array>
#include <vector>
#include <iostream>
#include <cstdint>
#include <cstring>

// includes, thirdparty
#ifdef FEAT_HAVE_FPARSER
#include <fparser.hh>
#endif // FEAT_HAVE_FPARSER

namespace FEAT
{
  namespace Geometry
  {
    /**
     * \brief Error class for VoxelMap related file errors
     */
    class VoxelMapFileError :
      public FileError
    {
    public:
      explicit VoxelMapFileError(const String& filename, const String& msg) :
        FileError(filename.empty() ? msg : (filename + ": " + msg))
      {
      }
    }; // class VoxelMapFileError

    /**
     * \brief Error class for VoxelMap related formula parsing errors
     */
    class VoxelMapFormulaParseError :
      public ParseError
    {
    public:
      explicit VoxelMapFormulaParseError(const String& msg) :
        ParseError(msg)
      {
      }
    }; // class VoxelMapFormulaParseError

    /**
     * \brief Interface for voxel masker for the VoxelMap class
     *
     * A voxel masker is basically an object that performs an inside-outside-test for some sort of
     * geometrical object for a given point in X/Y/Z coordinates. To reduce the overhead introduced
     * by virtual function calls, this class provides an abstract method #mask_line that will
     * perform the I/O-test for an entire row/line of points with common Y/Z coordinates at once
     * rather than for a single point.
     *
     * \author Peter Zajac
     */
    template<typename CoordType_, int dim_>
    class VoxelMasker
    {
    public:
      /// the coordinate type
      typedef CoordType_ CoordType;
      /// the point type
      typedef Tiny::Vector<CoordType_, dim_> PointType;

      /// virtual destructor
      virtual ~VoxelMasker()
      {
      }

      /**
       * \brief Computes the X-coordinate for a given point
       *
       * This helper function can be used by derived classes to compute the X-coordinate for a
       * single point on the mask line.
       *
       * \param[in] x_min, x_max
       * The minimal and maximal X-coordinate of the X-line for which the mask is to be computed.
       *
       * \param[in] i
       * The index of the current point in the X-line that is being tested
       *
       * \param[in] n
       * The number of points in the line; corresponds to mask.size()
       *
       * \returns
       * The X-coordinate of the current point
       */
      static CoordType x_coord(const CoordType x_min, const CoordType x_max, const std::size_t i, const std::size_t n)
      {
        return x_min + (x_max - x_min) * CoordType(i) / CoordType(n - 1u);
      }

      /**
       * \brief Computes the mask for an entire X-coordinate row
       *
       * This function is called by the VoxelMap class to perform the inside/outside-test for a line of points with
       * common Y/Z coordinates, where the X coordinates of the points have to be computed from the X-coordinate range
       * and the number points by calling the #x_coord helper function.
       *
       * \attention
       * This function \b must be implemented in a thread-safe fashion, since the VoxelMap class makes use of OpenMP
       * parallelization to speed up the map computation time.
       *
       * \param[inout] mask
       * The mask vector for the current X-row. Is allocated to correct size, but its contents are undefined upon entry.
       * Assuming that the length of this vector is n, then the i-th entry of the mask vector is to be set to 1 if
       * the point with the X-coordinate equal to x_coord(x_min, x_max, i, mask.size()) is inside the mask, otherwise
       * it is to be set to 0.
       *
       * \param[in] x_min, x_max
       * The minimal and maximal X-coordinate of the X-line for which the mask is to be computed.
       *
       * \param[in] point
       * The point that contains the Y- and (in 3D) Z-coordinates of the X-line that for which the mask is to be computed.
       * The first coordinate, i.e. point[0], is undefined, so only point[1] and (in 3D) point[2] are set to the Y-/Z-coords.
       */
      virtual void mask_line(std::vector<int>& mask, const CoordType x_min, const CoordType x_max, const PointType& point) = 0;
    }; // class VoxelMasker<...>

#if defined(FEAT_HAVE_CGAL) || defined(DOXYGEN)
    /**
     * \brief Voxel masker implementation for CGAL wrapper
     *
     * This class implements the voxel masker interface by using the I/O test functionality of a CGALWrapper object.
     *
     * \author Peter Zajac
     */
    template<typename CoordType_>
    class VoxelCGALMasker :
      public VoxelMasker<CoordType_, 3>
    {
    public:
      typedef VoxelMasker<CoordType_, 3> BaseClass;
      using typename BaseClass::CoordType;
      using typename BaseClass::PointType;

    private:
      /// a reference to our CGAL wrapper
      const Geometry::CGALWrapper<CoordType_>& _cgal_wrapper;
      const bool _invert;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] cgal_wrapper
       * A \resident reference to the CGAL wrapper object that is to be used for the inside/outside test.
       *
       * \param[in] invert
       * Specifies whether inside and outside are to be swapped.
       */
      explicit VoxelCGALMasker(const Geometry::CGALWrapper<CoordType_>& cgal_wrapper, bool invert) :
        _cgal_wrapper(cgal_wrapper),
        _invert(invert)
      {
      }

      virtual void mask_line(std::vector<int>& mask, const CoordType x_min, const CoordType x_max, const PointType& point) override
      {
        if(mask.empty())
          return;

        const std::size_t nv = mask.size();
        for(std::size_t i(0); i < nv; ++i)
        {
          mask[i] = (this->_cgal_wrapper.point_inside(BaseClass::x_coord(x_min, x_max, i, nv),
            point[1], point[2]) ? int(!_invert) : int(_invert));
        }
      }
    }; // class VoxelCGALMasker<...>
#endif // FEAT_HAVE_CGAL || DOXYGEN

    /**
     * \brief Voxel masker implementation for chart classes
     *
     * This class implements the voxel masker interface by using the signed distance functionality of the Atlas Chart implementations
     *
     * \author Peter Zajac
     */
    template<typename MeshType>
    class VoxelChartMasker :
      public VoxelMasker<typename MeshType::CoordType, MeshType::shape_dim>
    {
    public:
      typedef VoxelMasker<typename MeshType::CoordType, MeshType::shape_dim> BaseClass;
      using typename BaseClass::CoordType;
      using typename BaseClass::PointType;

    private:
      const Geometry::Atlas::ChartBase<MeshType>& _chart;
      const bool _invert;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] chart
       * A \resident reference to the chart object that is to be used for the inside/outside test.
       *
       * \param[in] invert
       * Specifies whether inside and outside are to be swapped.
       */
      explicit VoxelChartMasker(const Geometry::Atlas::ChartBase<MeshType>& chart, bool invert) :
        _chart(chart),
        _invert(invert)
      {
      }

      virtual void mask_line(std::vector<int>& mask, const CoordType x_min, const CoordType x_max, const PointType& point) override
      {
        if(mask.empty())
          return;

        PointType pt(point);

        const std::size_t nv = mask.size();
        for(std::size_t i(0); i < nv; ++i)
        {
          pt[0] = BaseClass::x_coord(x_min, x_max, i, nv);
          mask[i] = (_chart.signed_dist(pt) <= CoordType(0) ? int(!_invert) : int(_invert));
        }
      }
    }; // class VoxelChartMasker<...>

    /**
     * \brief Voxel masker implementation for lambda expression tests
     *
     * This class implements the voxel masker interface by evaluating a user-supplied lambda expression, which
     * performs the inside/outside test for each point.
     * A point is considered to be inside if the lambda expression evaluates to true.
     *
     * \author Peter Zajac
     */
    template<typename Lambda_, typename CoordType_, int dim_>
    class VoxelLambdaMasker :
      public VoxelMasker<CoordType_, dim_>
    {
    public:
      typedef VoxelMasker<CoordType_, dim_> BaseClass;
      using typename BaseClass::CoordType;
      using typename BaseClass::PointType;

    public:
      /// the lambda expression
      Lambda_ _lambda;

      /**
       * \brief Constructor
       *
       * \param[in] lambda
       * The lambda expression that is to be used to perform the inside/outside test.
       * The lambda expression is expected to take a const reference to a PointType object as a sole input
       * parameter and to return either true or false depending on whether the point is inside or outside.
       */
      explicit VoxelLambdaMasker(Lambda_&& lambda) :
        _lambda(std::forward<Lambda_>(lambda))
      {
      }

      virtual void mask_line(std::vector<int>& mask, const CoordType x_min, const CoordType x_max, const PointType& point) override
      {
        if(mask.empty())
          return;

        PointType pt(point);
        // define a const reference to ensure that the lambda expression does not modify the point
        const PointType&  cpt(pt);

        const std::size_t nv = mask.size();
        for(std::size_t i(0); i < nv; ++i)
        {
          pt[0] = BaseClass::x_coord(x_min, x_max, i, nv);
          mask[i] = (_lambda(cpt) ? 1 : 0);
        }
      }
    }; // class VoxelLambdaMasker<...>

#if defined(FEAT_HAVE_FPARSER) || defined(DOXYGEN)
    /**
     * \brief Voxel masker implementation for function parser expression tests
     *
     * This class implements the voxel masker interface by evaluating a user-supplied formula
     * in x,y,z at runtime, which performs the inside/outside test for each point.
     * A point is considered to be inside if the formula evaluates to a value > 0.
     *
     * \author Peter Zajac
     */
    template<int dim_>
    class VoxelFormulaMasker :
      public VoxelMasker<double, dim_>
    {
    public:
      typedef VoxelMasker<double, dim_> BaseClass;
      using typename BaseClass::CoordType;
      using typename BaseClass::PointType;

    public:
      /// the formula
      String _formula;
      /// the function parser
      ::FunctionParser _parser;

      /**
       * \brief Constructor
       *
       * \param[in] formula
       * The formula in x,y,z that is to be used to perform the inside/outside test.
       * The formula should return a value > 0 for all points that are considered to be inside and
       * a value <= 0 for all points which are considered to be outside.
       */
      explicit VoxelFormulaMasker(const String& formula) :
        _formula(formula),
        _parser()
      {
        _parser.AddConstant("pi", Math::pi<double>());
        _parser.AddConstant("eps", Math::eps<double>());

        String vars("x");
        if(dim_ > 1) vars += ",y";
        if(dim_ > 2) vars += ",z";

        // try to parse the function
        const int ret = _parser.Parse(_formula.c_str(), vars.c_str());
        if(ret >= 0)
        {
          String msg(_parser.ErrorMsg());
          msg.append("\n>>> '");
          msg.append(_formula);
          msg.append("'");
          throw VoxelMapFormulaParseError(msg);
        }

        // optimize the parsed function
        _parser.Optimize();
      }

      virtual void mask_line(std::vector<int>& mask, const CoordType x_min, const CoordType x_max, const PointType& point) override
      {
        if(mask.empty())
          return;

        PointType pt(point);

        const std::size_t nv = mask.size();
        for(std::size_t i(0); i < nv; ++i)
        {
          pt[0] = BaseClass::x_coord(x_min, x_max, i, nv);
          mask[i] = (_parser.Eval(pt.v) > 0.0 ? 1 : 0);
        }
      }
    }; // class VoxelFormulaMasker<...>
#endif // FEAT_HAVE_FPARSER || DOXYGEN

    /**
     * \brief VoxelMap class
     *
     * This class is responsible for creating an managing voxel maps, which are used by the Control::Domain::VoxelDomainControl
     * class to determine which elements can be dropped during the mesh hierarchy creation process.
     *
     * A voxel map is basically a 2D or 3D monochrome raster image representing a domain bounding box, where each voxel
     * in the map is either 1 or 0, depending on whether the point at the voxel's position is considered to be inside
     * or outside of the domain that is supposed to be discretized.
     *
     * Voxel maps can be created from analytic expressions which are either given as compile-time lambda expressions
     * or as runtime string expressions to be parsed or from a surface triangulation given as an OFF file, assuming
     * that the required third-party libraries 'CGAL' and 'fparser' are present.
     *
     * This class also provides the functionality to write a voxel map into a binary output file, optionally compressed
     * by the ZLIB library, along with the functionality to read these binary files, of course. It is strongly recommended
     * to write voxel map in compressed format, because the compression often shrinks the resulting file sizes by more
     * than 95%.
     *
     * To create a voxel map object, one has to perform the following steps (in this particular order):
     * -# Specify the bounding box for the voxel map by calling one of the following functions:
     *    - #set_bounding_box_2d for a 2D voxel map
     *    - #set_bounding_box_3d for a 3D voxel map
     * -# Specify the number of points for the voxel map in each dimension by calling one of the following functions:
     *    - #set_num_points to explicitly specify the number of points for each dimension
     *    - #set_resolution to compute the number of points for a maximum voxel resolution
     * -# Compute the voxel map contents by calling one of the following functions:
     *    - #compute_map to compute the entries from an object implementing the VoxelMasker interface
     *    - #compute_map_from_lambda_2d to compute the entries from a 2D lambda expression
     *    - #compute_map_from_lambda_3d to compute the entries from a 3D lambda expression
     *    - #compute_map_from_formula_2d to compute the entries from a 2D formula string (requires 'fparser' third-party library)
     *    - #compute_map_from_formula_3d to compute the entries from a 3D formula string (requires 'fparser' third-party library)
     *    - #compute_map_from_off_3d to compute the entries from a 3D OFF model (requires 'CGAL' third-party library)
     * -# Optional: Specify the out-of-bounds value by calling #set_out_of_bounds_value
     *
     * \brief Peter Zajac
     */
    class VoxelMap
    {
    public:
      /// signed 64-bit integer type
      typedef std::int64_t i64;

      /// unsigned 64-bit integer type
      typedef std::uint64_t u64;

      /// voxel map magic number
      static constexpr u64 magic_number = 0x70614D6C65786F56ull; // "VoxelMap"

      /// size of voxel map header in bytes
      static constexpr u64 header_size = 128ull;

      /// size of a real unit (meter) in our internal units (nanometers)
      static constexpr i64 unit_size = 1'000'000'000;

      /// header structure of voxel map file header
      struct FileHeader
      {
        u64 magic_number;     //< "VoxelMap" = 0x70614D6C65786F56ll
        u64 header_size;      //< 128 bytes
        i64 min_x;            //< minimum X-coordinate in nanometer
        i64 max_x;            //< maximum X-coordinate in nanometer
        u64 num_x;            //< number of points in X dimension
        u64 stride_line;      //< number of bytes for a single uncompressed X-line; a multiple of 16; usually stride_line = (((num_x + 7ull)/8ull + 15ull) & ~15ull)
        i64 min_y;            //< minimum Y-coordinate in nanometer
        i64 max_y;            //< maximum Y-coordinate in nanometer
        u64 num_y;            //< number of points in Y dimension
        u64 stride_plane;     //< number of bytes for a single uncompressed XY-plane; a multiple of 16; usually stride_plane = stride_line * num_y
        i64 min_z;            //< minimum Z-coordinate in nanometer
        i64 max_z;            //< maximum Z-coordinate in nanometer
        u64 num_z;            //< number of points in Z dimension
        u64 stride_volume;    //< number of bytes for a single uncompressed XYZ-volume; a multiple of 16; usually stride_volume = stride_plane * num_z
        u64 planes_per_block; //< number of planes per compression block (0 if uncompressed)
        u64 num_blocks;       //< total number of compression blocks (0 if uncompressed)
      }; // struct FileHeader: 128 bytes

      /// ensure that we're not compiling on some freak architecture...
      static_assert(sizeof(FileHeader) == header_size, "invalid FileHeader size");

      /// auxiliary function: compute default line stride from number of points in X dimension
      static u64 calc_line_stride(u64 num_x)
      {
        // compute number of bytes by rounding to next multiple of 8 then dividing by 8,
        // afterwards round up to next multiple of 16
        return (((num_x + 7ull)/8ull + 15ull) & ~15ull);
      }

      /**
       * \brief helper class for write function result
       *
       * An instance of this class is returned by the VoxelMap::write() function containing details
       * about the written file.
       */
      class WriteResult
      {
      public:
        u64 filesize;            //< total size of written file in bytes
        u64 compress_block_size; //< chosen maximum compression block size in MB
        u64 compress_level;      //< chosen compression level
        u64 num_compress;        //< number of compression blocks
        u64 compressed_size;     //< total compressed buffer size in bytes

        WriteResult() :
          filesize(0u),
          compress_block_size(0u),
          compress_level(0u),
          num_compress(0u),
          compressed_size(0u)
        {
        }

        explicit WriteResult(u64 fs, u64 cbs = 0u, u64 cl = 0u, u64 nc = 0u, u64 cs = 0u) :
          filesize(fs),
          compress_block_size(cbs),
          compress_level(cl),
          num_compress(nc),
          compressed_size(cs)
        {
        }

        /// Returns true if the file was written successfully, otherwise false
        operator bool() const
        {
          return filesize > u64(0);
        }
      }; // class WriteResult

      /**
       * \brief helper class for read function result
       *
       * An instance of this class is returned by the VoxelMap::read() function containing details
       * about the read file.
       */
      class ReadResult
      {
      public:
        u64 filesize;            //< total size of read file in bytes
        u64 num_compress;        //< number of compression blocks
        u64 compressed_size;     //< total compressed buffer size in bytes

        ReadResult() :
          filesize(0u),
          num_compress(0u),
          compressed_size(0u)
        {
        }

        explicit ReadResult(u64 fs, u64 nc = 0u, u64 cs = 0u) :
          filesize(fs),
          num_compress(nc),
          compressed_size(cs)
        {
        }

        /// Returns true if the file was read successfully, otherwise false
        operator bool() const
        {
          return filesize > u64(0);
        }
      }; // class ReadResult

    protected:
      /// current creation stage
      int _create_stage;
      /// bounding box dimensions of domain
      std::array<i64, 3u> _bbox_min, _bbox_max;
      /// number of points in each dimension
      std::array<u64, 3u> _num_points;
      /// stride of a single voxel line in X dimension
      u64 _stride_line;
      /// stride of a single plane in XY dimensions
      u64 _stride_plane;
      /// number of lines; 2D = num_points[2]; 3D: = num_points[2] * num_points[3]
      u64 _num_lines;
      /// number of planes; 2D = 1; 3D = num_points[3]
      u64 _num_planes;
      /// the actual voxel map; its size may be larger than necessary due to padding
      std::vector<char> _voxel_map;
      /// the out-of-bounds-value for the voxel map
      bool _out_of_bounds_value;

    public:
      /// standard constructor
      VoxelMap();

      /// no copies, no problems
      VoxelMap(const VoxelMap&) = delete;
      VoxelMap& operator=(const VoxelMap&) = delete;

      /// move constructor
      VoxelMap(VoxelMap&& other);

      /// move assign operator
      VoxelMap& operator=(VoxelMap&& other);

      /// virtual destructor
      virtual ~VoxelMap();

      /**
       * \brief Sets the bounding box for a 2D voxel map
       *
       * \param[in] x_min, x_max
       * The range of the voxel map in X dimension
       *
       * \param[in] y_min, y_max
       * The range of the voxel map in Y dimension
       */
      void set_bounding_box_2d(Real x_min, Real x_max, Real y_min, Real y_max);

      /**
       * \brief Sets the bounding box for a 3D voxel map
       *
       * \param[in] x_min, x_max
       * The range of the voxel map in X dimension
       *
       * \param[in] y_min, y_max
       * The range of the voxel map in Y dimension
       *
       * \param[in] z_min, z_max
       * The range of the voxel map in Z dimension
       */
      void set_bounding_box_3d(Real x_min, Real x_max, Real y_min, Real y_max, Real z_min, Real z_max);

      /// Returns the minimum bounding box coordinate for a given dimensions
      Real get_bounding_box_min(int dim) const
      {
        XASSERT((dim >= 0) && (dim < 3));
        return Real(this->_bbox_min[std::size_t(dim)]) / Real(unit_size);
      }

      /// Returns the maximum bounding box coordinate for a given dimensions
      Real get_bounding_box_max(int dim) const
      {
        XASSERT((dim >= 0) && (dim < 3));
        return Real(this->_bbox_max[std::size_t(dim)]) / Real(unit_size);
      }

      /**
       * \brief Sets the out-of-bounds value for the voxel map
       *
       * The out-of-bounds-value specifies whether the outside of the voxel map bounding box is
       * being treated as 'inside' (true) or 'outside' (false), which is the value that is returned
       * by the check_point and check_box functions if the point/bounding box is not inside of
       * the voxel map bounding box.
       *
       * \param[in] value
       * The out-of-bounds-value for the voxel map
       */
      void set_out_of_bounds_value(bool value)
      {
        this->_out_of_bounds_value = value;
      }

      /// Returns the out-of-bounds value for the voxel map
      bool get_out_of_bounds_value() const
      {
        return this->_out_of_bounds_value;
      }

      /**
       * \brief Sets the number of points of the voxel map in each dimension
       *
       * \param[in] num_x, num_y, num_z
       * The number of points in each dimension
       */
      void set_num_points(Index num_x, Index num_y, Index num_z = 0u);

      /// Returns the number of points of the voxel map in each dimension
      Index get_num_points(int dim) const
      {
        return Index(this->_num_points.at(Index(dim)));
      }

      /**
       * \brief Sets the resolution for the voxel map, i.e. the maximum voxel size in each dimension
       *
       * This function computes the number of points in each dimension to ensure that the size
       * of a voxel is less or equal to max_res in any dimension.
       *
       * \param[in] max_res
       * The maximum voxel size; must be > 0.
       */
      void set_resolution(Real max_res);

      /// Returns the total number of voxels in the voxel map
      Index get_num_voxels() const
      {
        return Index(this->_num_points[0] * this->_num_points[1] * this->_num_points[2]);
      }

      /// Returns the line stride, i.e. the size of a single X-line in bytes
      u64 get_stride_line() const
      {
        return this->_stride_line;
      }

      /// Returns the plane stride, i.e. the size of a single XY-plane in bytes
      u64 get_stride_plane() const
      {
        return this->_stride_plane;
      }

      /// Returns the volume stride, i.e. the size of a single XYZ-volume in bytes
      u64 get_stride_volume() const
      {
        return u64(this->_voxel_map.size());
      }

      /// Returns a reference to the voxel map
      const std::vector<char>& get_map() const
      {
        return this->_voxel_map;
      }

      /**
       * \brief Creates the voxel map based on a given masker object by utilizing MPI parallelization
       *
       * In an MPI-parallel simulation, this function must be called collectively by all processes
       * in the communicator \p comm, since this function splits the workload across all these
       * processes.
       *
       * \attention
       * This function can only be called once the bounding box dimensions as well as the number
       * of points or the resolution have been set.
       *
       * \param[in] comm
       * A \transient reference to the communicator to use
       *
       * \param[in] masker
       * A \transient reference to the masker object
       *
       * \param[in] gather_to_all
       * Specifies whether the voxel map is to be gathered on all processes in the communicator.
       * The only reason to set this to \b false is when the map is only required on rank 0 for file output.
       */
      template<typename Coord_, int dim_>
      void compute_map(const Dist::Comm& comm, VoxelMasker<Coord_, dim_>& masker, bool gather_to_all = true)
      {
        XASSERT(!_voxel_map.empty());
        this->_compute_voxel_map(comm, masker, gather_to_all);
      }

      /**
       * \brief Creates the voxel map based on a chart by utilizing MPI parallelization
       *
       * In an MPI-parallel simulation, this function must be called collectively by all processes
       * in the communicator \p comm, since this function splits the workload across all these
       * processes.
       *
       * \attention
       * This function can only be called once the bounding box dimensions as well as the number
       * of points or the resolution have been set.
       *
       * \param[in] comm
       * A \transient reference to the communicator to use
       *
       * \param[in] chart
       * A \transient reference to the chart
       *
       * \param[in] invert
       * Specifies whether the definition of inside and outside of the chart is to be swapped.
       *
       * \param[in] gather_to_all
       * Specifies whether the voxel map is to be gathered on all processes in the communicator.
       * The only reason to set this to \b false is when the map is only required on rank 0 for file output.
       */
      template<typename MeshType_>
      void compute_map_from_chart(const Dist::Comm& comm, const Geometry::Atlas::ChartBase<MeshType_>& chart, bool invert, bool gather_to_all = true)
      {
        XASSERT(!_voxel_map.empty());
        VoxelChartMasker<MeshType_> masker(chart, invert);
        this->_compute_voxel_map(comm, masker, gather_to_all);
      }

      /**
       * \brief Creates the voxel map based on a 2D lambda expression by utilizing MPI parallelization
       *
       * In an MPI-parallel simulation, this function must be called collectively by all processes
       * in the communicator \p comm, since this function splits the workload across all these
       * processes.
       *
       * \attention
       * This function can only be called once the bounding box dimensions as well as the number
       * of points or the resolution have been set.
       *
       * \param[in] comm
       * A \transient reference to the communicator to use
       *
       * \param[in] lambda
       * A lambda expression that accepts a Tiny::Vector<Real,2> as the only input parameter and returns
       * \c true, if the point is inside the masked domain or \c false otherwise. This expression must be
       * identical on all MPI processes.
       *
       * \param[in] gather_to_all
       * Specifies whether the voxel map is to be gathered on all processes in the communicator.
       * The only reason to set this to \b false is when the map is only required on rank 0 for file output.
       */
      template<typename Lambda_>
      void compute_map_from_lambda_2d(const Dist::Comm& comm, Lambda_&& lambda, bool gather_to_all = true)
      {
        VoxelLambdaMasker<Lambda_, Real, 2> masker(std::forward<Lambda_>(lambda));
        this->_compute_voxel_map(comm, masker, gather_to_all);
      }

      /**
       * \brief Creates the voxel map based on a 3D lambda expression by utilizing MPI parallelization
       *
       * In an MPI-parallel simulation, this function must be called collectively by all processes
       * in the communicator \p comm, since this function splits the workload across all these
       * processes.
       *
       * \attention
       * This function can only be called once the bounding box dimensions as well as the number
       * of points or the resolution have been set.
       *
       * \param[in] comm
       * A \transient reference to the communicator to use
       *
       * \param[in] lambda
       * A lambda expression that accepts a Tiny::Vector<Real,3> as the only input parameter and returns
       * \c true, if the point is inside the masked domain or \c false otherwise. This expression must be
       * identical on all MPI processes.
       *
       * \param[in] gather_to_all
       * Specifies whether the voxel map is to be gathered on all processes in the communicator.
       * The only reason to set this to \b false is when the map is only required on rank 0 for file output.
       */
      template<typename Lambda_>
      void compute_map_from_lambda_3d(const Dist::Comm& comm, Lambda_&& lambda, bool gather_to_all = true)
      {
        VoxelLambdaMasker<Lambda_, Real, 3> masker(std::forward<Lambda_>(lambda));
        this->_compute_voxel_map(comm, masker, gather_to_all);
      }

      /**
       * \brief Creates the voxel map based on a 2D formula by utilizing MPI parallelization
       *
       * In an MPI-parallel simulation, this function must be called collectively by all processes
       * in the communicator \p comm, since this function splits the workload across all these
       * processes.
       *
       * \attention
       * This function can only be called once the bounding box dimensions as well as the number
       * of points or the resolution have been set.
       *
       * \param[in] comm
       * A \transient reference to the communicator to use
       *
       * \param[in] formula
       * A formula in x,y that evaluates to a value >= 0, if the points is inside the masked domain
       * or to a value <= 0 otherwise.
       *
       * \param[in] gather_to_all
       * Specifies whether the voxel map is to be gathered on all processes in the communicator.
       * The only reason to set this to \b false is when the map is only required on rank 0 for file output.
       */
      void compute_map_from_formula_2d(const Dist::Comm& comm, const String& formula, bool gather_to_all = true);

      /**
       * \brief Creates the voxel map based on a 3D formula by utilizing MPI parallelization
       *
       * In an MPI-parallel simulation, this function must be called collectively by all processes
       * in the communicator \p comm, since this function splits the workload across all these
       * processes.
       *
       * \attention
       * This function can only be called once the bounding box dimensions as well as the number
       * of points or the resolution have been set.
       *
       * \param[in] comm
       * A \transient reference to the communicator to use
       *
       * \param[in] formula
       * A formula in x,y,z that evaluates to a value >= 0, if the points is inside the masked domain
       * or to a value <= 0 otherwise.
       *
       * \param[in] gather_to_all
       * Specifies whether the voxel map is to be gathered on all processes in the communicator.
       * The only reason to set this to \b false is when the map is only required on rank 0 for file output.
       */
      void compute_map_from_formula_3d(const Dist::Comm& comm, const String& formula, bool gather_to_all = true);

      /**
       * \brief Creates the voxel map based on a 3D OFF model handled by CGAL by utilizing MPI parallelization
       *
       * In an MPI-parallel simulation, this function must be called collectively by all processes
       * in the communicator \p comm, since this function splits the workload across all these
       * processes.
       *
       * \attention
       * This function can only be called once the bounding box dimensions as well as the number
       * of points or the resolution have been set.
       *
       * \param[in] comm
       * A \transient reference to the communicator to use
       *
       * \param[in] filename
       * The filename of the OFF model file that is to be read in. This file will be read in using the
       * DistFileIO functionality and broadcasted across all processes in the communicator.
       *
       * \param[in] invert
       * Specifies whether the definition of inside and outside of the model is to be swapped.
       *
       * \param[in] gather_to_all
       * Specifies whether the voxel map is to be gathered on all processes in the communicator.
       * The only reason to set this to \b false is when the map is only required on rank 0 for file output.
       */
      void compute_map_from_off_3d(const Dist::Comm& comm, const String& filename, bool invert, bool gather_to_all = true);

      /**
       * \brief Creates the voxel map based on a 3D OFF model handled by CGAL by utilizing MPI parallelization
       *
       * In an MPI-parallel simulation, this function must be called collectively by all processes
       * in the communicator \p comm, since this function splits the workload across all these
       * processes.
       *
       * \attention
       * This function can only be called once the bounding box dimensions as well as the number
       * of points or the resolution have been set.
       *
       * \param[in] comm
       * A \transient reference to the communicator to use
       *
       * \param[in] is
       * The input stream that the OFF file contents are to be parsed from. Must have the identical content
       * on all MPI processes.
       *
       * \param[in]
       * Specifies whether the definition of inside and outside of the model is to be swapped.
       *
       * \param[in] gather_to_all
       * Specifies whether the voxel map is to be gathered on all processes in the communicator.
       * The only reason to set this to \b false is when the map is only required on rank 0 for file output.
       */
      void compute_map_from_off_3d(const Dist::Comm& comm, std::istream& is, bool invert, bool gather_to_all = true);

      /**
       * \brief Writes the voxel map into a file that can later be read in
       *
       * \attention
       * In an MPI-parallel simulation, only one process shall execute this function!
       *
       * \param[in] filename
       * The filename of the output file. The recommended extension is 'vxl'.
       *
       * \param[in] compress_block_size
       * The maximum size of a single ZLIB compression block in MB. Set to 0 to disable ZLIB compression.
       *
       * \param[in] compress_level
       * The ZLIB compression level in range [1,9], where 1 is fastest compression and 9 is best compression.
       *
       * \returns
       * A WriteResult object that contains more detailed information about the output file, including bytes
       * written and total compression buffer size. Can be cast to bool to check whether the write was successful.
       */
      WriteResult write(const String& filename, const u64 compress_block_size = 128u, const int compress_level = 9) const;

      /**
       * \brief Writes the voxel map into a binary output stream
       *
       * \attention
       * In an MPI-parallel simulation, only one process shall execute this function!
       *
       * \param[out] os
       * The binary output stream to write to.
       *
       * \param[in] compress_block_size
       * The maximum size of a single ZLIB compression block in MB. Set to 0 to disable ZLIB compression.
       *
       * \param[in] compress_level
       * The ZLIB compression level in range [1,9], where 1 is fastest compression and 9 is best compression.
       *
       * \returns
       * A WriteResult object that contains more detailed information about the output file, including bytes
       * written and total compression buffer size. Can be cast to bool to check whether the write was successful.
       */
      WriteResult write(std::ostream& os, const u64 compress_block_size = 128u, const int compress_level = 9) const;

      /**
       * \brief Reads a voxel map from a file
       *
       * \param[in] comm
       * A \transient reference to a communicator.
       *
       * \param[in] filename
       * The filename of the voxel map file to be read in.
       *
       * \returns
       * A ReadResult object that contains more detailed information about the input file, including bytes
       * read and total compression buffer size. Can be cast to bool to check whether the read was successful.
       */
      ReadResult read(const Dist::Comm& comm, const String& filename);

      /**
       * \brief Reads a voxel map from a binary input stream
       *
       * \param[in] is
       * A \transient reference to a binary input stream to read form
       *
       * \param[in] filename
       * The filename of the voxel map file to be read in. This is just used for error output.
       */
      ReadResult read(std::istream& is, const String& filename = "");

      /**
       * \brief Returns the value of a specific voxel point in the voxel map
       *
       * \param[in] ix, iy, iz
       * The X/Y/Z indices of the point in the voxel map; each must be in range {0, ..., num_points[dim]-1}.
       *
       * \returns
       * The value of the point at the desired voxel map position.
       */
      bool check_point(u64 ix, u64 iy = 0u, u64 iz = 0u) const
      {
        ASSERT((ix < this->_num_points[0]) || ((this->_num_points[0] == 0u) && (ix == 0u)));
        ASSERT((iy < this->_num_points[1]) || ((this->_num_points[1] == 0u) && (iy == 0u)));
        ASSERT((iz < this->_num_points[2]) || ((this->_num_points[2] == 0u) && (iz == 0u)));
        return this->_check_point(i64(ix), i64(iy), i64(iz));
      }

      /**
       * \brief Checks the voxel map entry for a given point
       *
       * \param[in] point
       * The point whose voxel map entry is to be checked.
       *
       * \returns
       * The value of the voxel map in the position represented by point if the point is inside the
       * voxel map bounding box, otherwise returns the out-of-bounds-value.
       */
      template<typename Coord_, int dim_>
      bool check_point_nearest(const Tiny::Vector<Coord_, dim_>& point) const
      {
        // convert point to internal units
        std::array<i64, dim_> ipt;
        for(int i(0); i < dim_; ++i)
          ipt[std::size_t(i)] = i64(point[i] * Coord_(unit_size));
        return this->_check_point_nearest(ipt);
      }

      /**
       * \brief Samples the voxel map entry for a given point by multi-linear interpolation
       *
       * \param[in] point
       * The point whose voxel map entry is to be samples.
       *
       * \returns
       * The multi-linearly interpolated value of the voxel map in the position represented by point.
       */
      template<typename Coord_, int dim_>
      Real sample_point(const Tiny::Vector<Coord_, dim_>& point) const
      {
        // convert point to internal units
        std::array<i64, dim_> ipt;
        for(int i(0); i < dim_; ++i)
          ipt[std::size_t(i)] = i64(point[i] * Coord_(unit_size));
        return this->_sample_point(ipt);
      }

      /**
       * \brief Checks the voxel map entries for a given box of voxels
       *
       * \param[in] box_min
       * The minimum coordinates for the voxel box that is to be checked
       *
       * \param[in] box_max
       * The maximum coordinates for the voxel box that is to be checked
       *
       * \note
       * If the box represented by \p box_min and \p box_max is clamped to the bounding box of the voxel map.
       *
       * \returns
       * \c true, if at least one voxel contained in the box given by \p box_min and \p box_max is
       * set to \c true, otherwise \c false.
       */
      template<typename Coord_, int dim_>
      bool check_box(const Tiny::Vector<Coord_, dim_>& box_min, const Tiny::Vector<Coord_, dim_>& box_max) const
      {
        // convert coordinates to map indices
        std::array<u64, dim_> ibox_min, ibox_max;
        for(std::size_t i(0); i < std::size_t(dim_); ++i)
        {
          // make sure the map indices do not exceed the bounding box of the voxel map
          i64 imin = _map_coord_idx_lower(i64(box_min[int(i)] * Coord_(unit_size)), i);
          if(imin < i64(0))
          {
            if(this->_out_of_bounds_value)
              return true;
            ibox_min[i] = u64(0);
          }
          else
            ibox_min[i] = u64(imin);

          i64 imax = _map_coord_idx_upper(i64(box_max[int(i)] * Coord_(unit_size)), i);
          if(i64(this->_num_points[i]) < imax + 1)
          {
            if(this->_out_of_bounds_value)
              return true;
            ibox_max[i] = i64(this->_num_points[i]) - 1;
          }
          else
            ibox_max[i] = u64(imax);
        }

        return this->_check_box(ibox_min, ibox_max);
      }

      /**
       * \brief Samples the voxel map entries for a given box of voxels
       *
       * \param[in] box_min
       * The minimum coordinates for the voxel box that is to be sampled
       *
       * \param[in] box_max
       * The maximum coordinates for the voxel box that is to be sampled
       *
       * \note
       * The box represented by \p box_min and \p box_max is clamped to the bounding box of the voxel map.
       *
       * \returns
       * The average number of voxel map entries in the box which are set to true.
       *
       * \note
       * The return value is 1 if all voxels in the box are set to 1 and the return value is 0 if
       * all voxels in the box are 0; in any other case the return value is > 0 and < 1, representing
       * the average value of all voxels in the sampled box.
       */
      template<typename Coord_, int dim_>
      Real sample_box(const Tiny::Vector<Coord_, dim_>& box_min, const Tiny::Vector<Coord_, dim_>& box_max) const
      {
        std::array<i64, dim_> ibox_min, ibox_max;
        for(std::size_t i(0); i < std::size_t(dim_); ++i)
        {
          ibox_min[i] = i64(box_min[int(i)] * Coord_(unit_size));
          ibox_max[i] = i64(box_max[int(i)] * Coord_(unit_size));
        }
        return this->_sample_box(ibox_min, ibox_max);
      }

      /**
       * \brief Exports the voxel map to a sequence of monochrome BMP image files
       *
       * This function will export the entire voxel map in full resolution as a sequence of
       * monochrome bitmap files. The maximum voxel map dimension is 20Kx20Kx20K, which would
       * result in 20000 BMP files with a resolution 400 MPixels each.
       *
       * \note
       * If you want to export the voxel map in a lower resolution, consider using the
       * render_to_bmp() function instead.
       *
       * \param[in] filename_prefix
       * The filename prefix for the image file sequence.
       */
      void export_to_bmp(const String& filename_prefix) const;

      /**
       * \brief Exports a single plane of the voxel map to a monochrome BMP image file
       *
       * This function will export an entire plane of the voxel map in full resolution as a
       * monochrome bitmap files.
       *
       * \note
       * If you want to export the voxel map plane in a lower resolution, consider using the
       * render_plane_to_bmp() function instead.
       *
       * \param[in] filename
       * The filename for the BMP image file
       *
       * \param[in] z_coord
       * The Z-coordinate of the plane that is to be exported.
       */
      void export_plane_to_bmp(const String& filename, Real z_coord) const;

      /**
       * \brief Exports a single plane of the voxel map to a monochrome BMP image file
       *
       * This function will export an entire plane of the voxel map in full resolution as a
       * monochrome bitmap files.
       *
       * \note
       * If you want to export the voxel map plane in a lower resolution, consider using the
       * render_plane_to_bmp() function instead.
       *
       * \param[out] os
       * The binary output stream to write to
       *
       * \param[in] z_coord
       * The Z-coordinate of the plane that is to be exported.
       */
      void export_plane_to_bmp(std::ostream& os, Real z_coord) const;

      /**
       * \brief Renders the voxel map to a sequence of gray-scale BMP image files
       *
       * This function will render the entire voxel map in reduced resolution as a sequence of
       * gray-scale bitmap files, where each pixel contains the average value of all voxels that
       * have been reduced to it.
       *
       * \param[in] filename_prefix
       * The filename prefix for the image file sequence
       *
       * \param[in] width
       * The desired bitmap width
       *
       * \param[in] height
       * The desired bitmap height
       *
       * \param[in] depth
       * The desired number of bitmap images in the sequence
       */
      void render_to_bmp(const String& filename_prefix, Index width, Index height, Index depth) const;

      /**
       * \brief Renders a plane range of voxel map to a single gray-scale BMP image file
       *
       * This function will render a range of plane of the voxel map in reduced resolution as a
       * single gray-scale bitmap file, where each pixel contains the average value of all voxels that
       * have been reduced to it.
       *
       * \param[in] filename
       * The filename for the BMP image
       *
       * \param[in] width
       * The desired bitmap width
       *
       * \param[in] height
       * The desired bitmap height
       *
       * \param[in] z_min, z_max
       * The Z-coordinates of the plane range that has to be sampled into a single output plane
       */
      void render_plane_to_bmp(const String& filename, Index width, Index height, Real z_min, Real z_max) const;

      /**
       * \brief Renders a plane range of voxel map to a single gray-scale BMP image file stream
       *
       * This function will render a range of plane of the voxel map in reduced resolution as a
       * single gray-scale bitmap file, where each pixel contains the average value of all voxels that
       * have been reduced to it.
       *
       * \param[out] os
       * The output stream to write to
       *
       * \param[in] width
       * The desired bitmap width
       *
       * \param[in] height
       * The desired bitmap height
       *
       * \param[in] z_min, z_max
       * The Z-coordinates of the plane range that has to be sampled into a single output plane
       */
      void render_plane_to_bmp(std::ostream& os, Index width, Index height, Real z_min, Real z_max) const;

    protected:
      /**
       * \brief Maps a unit coordinate to a X-/Y-/Z-index of the closest voxel in the voxel map
       *
       * \param[in] xyz
       * The X-/Y-/Z unit coordinate that is to be mapped
       *
       * \param[in] sdim
       * Specifies which dimension is to be mapped: 0=X, 1=Y, 2=Z
       *
       * \returns
       * The index of the closest voxel to the given unit coordinate or ~u64(0),
       * if the unit coordinate was outside of the voxel map bounding box.
       */
      u64 _map_coord_idx_nearest(i64 xyz, std::size_t sdim) const;

      /**
       * \brief Maps a unit coordinate to a X-/Y-/Z-index of the lower voxel in the voxel map
       *
       * \attention
       * This function intentionally does not check the point against the valid voxel map bounding box, i.e.
       * the returned index may be outside of the interval {0, ..., num_points[sdim]-1}.
       *
       * \param[in] xyz
       * The X-/Y-/Z unit coordinate that is to be mapped
       *
       * \param[in] sdim
       * Specifies which dimension is to be mapped: 0=X, 1=Y, 2=Z
       *
       * \returns
       * The index of the lower voxel to the given unit coordinate, i.e. the voxel whose unit coordinate
       * is less than or equal to the given point unit coordinate \p xyz.
       */
      i64 _map_coord_idx_lower(i64 xyz, std::size_t sdim) const;

      /**
       * \brief Maps a unit coordinate to a X-/Y-/Z-index of the upper voxel in the voxel map
       *
       * \attention
       * This function intentionally does not check the point against the valid voxel map bounding box, i.e.
       * the returned index may be outside of the interval {0, ..., num_points[sdim]-1}.
       *
       * \param[in] xyz
       * The X-/Y-/Z unit coordinate that is to be mapped
       *
       * \param[in] sdim
       * Specifies which dimension is to be mapped: 0=X, 1=Y, 2=Z
       *
       * \returns
       * The index of the upper voxel to the given unit coordinate, i.e. the voxel whose unit coordinate
       * is greater than or equal to the given point unit coordinate \p xyz.
       */
      i64 _map_coord_idx_upper(i64 xyz, std::size_t sdim) const;

      /**
       * \brief Helper function for _sample_voxel_map_Xd: compute the sample rate in a given dimension
       *
       * The sample rate of a point xyz with respect to the voxel v(idx) is defined as being equal to 1 if xyz = V(idx)
       * and to fall linearly to 0 for xyz = V(idx-1) and xyz = V(idx+1), which can be achieved by the formula
       *
       *                        sr(xyz,idx) :=  1 - |xyz - V(idx)| / |V(idx+1) - V(idx)|
       *
       * note that the denominator |V(idx+1) - V(idx)| := 1/h is identical for all idx
       *
       * \param[in] xyz
       * The X-/Y-/Z unit coordinate of the point for which the sample rate is to be computed
       *
       * \param[in] idx
       * The voxel index for the coordinate \p xyz as returned by _map_coord_idx_lower
       *
       * \returns The sample rate of the point \p xyz w.r.t. the voxel V(idx)
       */
      Real _calc_sample_rate(i64 xyz, i64 idx, std::size_t sdim) const;

      /// Helper function for _sample_point: checks a single voxel point value
      bool _check_point(i64 xidx, i64 yidx, i64 zidx) const;

      /**
       * \brief Checks the nearest voxel map entry for a given 1D point
       *
       * \param[in] p
       * The unit coordinates of the point whose closest voxel in the voxel map is to be checked.
       *
       * \returns
       * \c true, if the closest voxel to the point is inside the voxel map domain (i.e. the voxel map entry is 1),
       * or \c false, if the closest voxel is outside of the voxel map domain.
       */
      bool _check_point_nearest(const std::array<i64, 1>& p) const;

      /**
       * \brief Checks the nearest voxel map entry for a given 2D point
       *
       * \param[in] p
       * The unit coordinates of the point whose closest voxel in the voxel map is to be checked.
       *
       * \returns
       * \c true, if the closest voxel to the point is inside the voxel map domain (i.e. the voxel map entry is 1),
       * or \c false, if the closest voxel is outside of the voxel map domain.
       */
      bool _check_point_nearest(const std::array<i64, 2>& p) const;

      /**
       * \brief Checks the nearest voxel map entry for a given 3D point
       *
       * \param[in] p
       * The unit coordinates of the point whose closest voxel in the voxel map is to be checked.
       *
       * \returns
       * \c true, if the closest voxel to the point is inside the voxel map domain (i.e. the voxel map entry is 1),
       * or \c false, if the closest voxel is outside of the voxel map domain.
       */
      bool _check_point_nearest(const std::array<i64, 3>& p) const;

      /**
       * \brief Checks whether a 1D box contains at least one active voxel
       *
       * \param[in] box_min, box_max
       * The minimum and maximum X indices that make up the box to be checked
       *
       * \returns \c true, if at least one voxel in the given box is set to \c true, otherwise \c false.
       */
      bool _check_box(const std::array<u64, 1>& box_min, const std::array<u64, 1>& box_max) const;

      /**
       * \brief Checks whether a 2D box contains at least one active voxel
       *
       * \param[in] box_min, box_max
       * The minimum and maximum X/Y indices that make up the box to be checked
       *
       * \returns \c true, if at least one voxel in the given box is set to \c true, otherwise \c false.
       */
      bool _check_box(const std::array<u64, 2>& box_min, const std::array<u64, 2>& box_max) const;

      /**
       * \brief Checks whether a 3D box contains at least one active voxel
       *
       * \param[in] box_min, box_max
       * The minimum and maximum X/Y/Z indices that make up the box to be checked
       *
       * \returns \c true, if at least one voxel in the given box is set to \c true, otherwise \c false.
       */
      bool _check_box(const std::array<u64, 3>& box_min, const std::array<u64, 3>& box_max) const;

      /// Helper function for _sample_point: samples a 1D edge in X-parallel line
      Real _sample_point_1d_x(i64 px, i64 xidx, i64 yidx = i64(0), i64 zidx = i64(0)) const;

      /// Helper function for _sample_point: samples a 1D edge in Y-parallel line
      Real _sample_point_1d_y(i64 py, i64 xidx, i64 yidx, i64 zidx = i64(0)) const;

      /// Helper function for _sample_point: samples a 1D edge in Z-parallel line
      Real _sample_point_1d_z(i64 pz, i64 xidx, i64 yidx, i64 zidx) const;

      /// Helper function for _sample_point: samples a 2D square in XY-parallel plane
      Real _sample_point_2d(i64 px, i64 py, i64 xidx, i64 yidx, i64 zidx = i64(0)) const;

      /// Helper function for _sample_point: samples a 3D cube in XYZ volume
      Real _sample_point_3d(i64 px, i64 py, i64 pz, i64 xidx, i64 yidx, i64 zidx) const;

      /**
       * \brief Samples the voxel map for a given 1D point by multi-linear interpolation
       *
       * \param[in] p
       * The unit coordinates of the point whose voxel map value is to be sampled.
       *
       * \returns
       * The multi-linearly interpolated value of the eight voxels surrounding the point
       */
      Real _sample_point(const std::array<i64, 1>& p) const;

      /**
       * \brief Samples the voxel map for a given 2D point by multi-linear interpolation
       *
       * \param[in] p
       * The unit coordinates of the point whose voxel map value is to be sampled.
       *
       * \returns
       * The multi-linearly interpolated value of the eight voxels surrounding the point
       */
      Real _sample_point(const std::array<i64, 2>& p) const;

      /**
       * \brief Samples the voxel map for a given 3D point by multi-linear interpolation
       *
       * \param[in] p
       * The unit coordinates of the point whose voxel map value is to be sampled.
       *
       * \returns
       * The multi-linearly interpolated value of the eight voxels surrounding the point
       */
      Real _sample_point(const std::array<i64, 3>& p) const;

      /// Helper function for _sample_box: samples a 1D edge in X line
      Real _sample_box_1d(i64 px0, i64 px1, i64 xidx0, i64 xidx1, i64 yidx = i64(0), i64 zidx = i64(0)) const;

      /// Helper function for _sample_box: samples a 2D rectangle in XY plane
      Real _sample_box_2d(i64 px0, i64 px1, i64 py0, i64 py1, i64 xidx0, i64 xidx1, i64 yidx0, i64 yidx1, i64 zidx = i64(0)) const;

      /// Helper function for _sample_box: samples a 3D cuboid in XYZ volume
      Real _sample_box_3d(i64 px0, i64 px1, i64 py0, i64 py1, i64 pz0, i64 pz1, i64 xidx0, i64 xidx1, i64 yidx0, i64 yidx1, i64 zidx0, i64 zidx1) const;

      /**
       * \brief Samples a 1D box of the voxel map
       *
       * \param[in] box_min, box_max
       * The minimum and maximum X unit coordinates that make up the box to be sampled
       *
       * \returns The average voxel map value of the sample box
       */
      Real _sample_box(const std::array<i64, 1>& box_min, const std::array<i64, 1>& box_max) const;

      /**
       * \brief Samples a 2D box of the voxel map
       *
       * \param[in] box_min, box_max
       * The minimum and maximum X/Y unit coordinates that make up the box to be sampled
       *
       * \returns The average voxel map value of the sample box
       */
      Real _sample_box(const std::array<i64, 2>& box_min, const std::array<i64, 2>& box_max) const;

      /**
       * \brief Samples a 3D box of the voxel map
       *
       * \param[in] box_min, box_max
       * The minimum and maximum X/Y/Z unit coordinates that make up the box to be sampled
       *
       * \returns The average voxel map value of the sample box
       */
      Real _sample_box(const std::array<i64, 3>& box_min, const std::array<i64, 3>& box_max) const;

      /**
       * \brief Compresses a voxel map line into the voxel map.
       *
       * \param[in] mask
       * The voxel map line as computed by a VoxelMasker implementation
       *
       * \param[in] line
       * The index of the line that is to be compressed.
       */
      void _compress_voxel_map_line(const std::vector<int>& mask, const u64 line);

      /**
       * \brief Computes a range of voxel map lines for a 1D domain
       *
       * \param[in] masker
       * A \transient reference to the masker that is to be used to compute the mask
       *
       * \param[in] beg
       * The index of the first voxel map line that is to be computed. Must be equal to 0 in 1D.
       *
       * \param[in] end
       * The index of the last voxel map line that is to be computed plus one. Must be equal to 1 in 1D.
       *
       * \param[in] offset
       * The offset of the first voxel map line that is to be computed; Must be equal to 0 in 1D.
       */
      template<typename CoordType_>
      void _compute_voxel_map_lines(VoxelMasker<CoordType_, 1>& masker, u64 DOXY(beg), u64 DOXY(end), u64 DOXY(offset))
      {
        XASSERTM(_num_points[1] == 1ull, "cannot use a 1D VoxelMasker to create a 2D/3D voxel mask");
        XASSERTM(_num_points[2] == 1ull, "cannot use a 1D VoxelMasker to create a 2D/3D voxel mask");
        const CoordType_ x_min = CoordType_(_bbox_min[0]) / CoordType_(unit_size);
        const CoordType_ x_max = CoordType_(_bbox_max[0]) / CoordType_(unit_size);
        std::vector<int> line_mask(this->_num_points[0], 0);
        Tiny::Vector<CoordType_, 1> coords;
        masker.mask_line(line_mask, x_min, x_max, coords);
        this->_compress_voxel_map_line(line_mask, 0);
      }

      /**
       * \brief Computes a range of voxel map lines for a 2D domain
       *
       * \param[in] masker
       * A \transient reference to the masker that is to be used to compute the mask
       *
       * \param[in] beg
       * The index of the first voxel map line that is to be computed.
       *
       * \param[in] end
       * The index of the last voxel map line that is to be computed plus one.
       *
       * \param[in] offset
       * The offset of the first voxel map line that is to be computed; this is typically equal to 0,
       * unless we are only computing a sub-map that is going to be gathered on the root later on,
       * in which case it may be set equal to \p beg.
       */
      template<typename CoordType_>
      void _compute_voxel_map_lines(VoxelMasker<CoordType_, 2>& masker, u64 beg, u64 end, u64 offset)
      {
        XASSERTM(_num_points[2] == 1ull, "cannot use a 2D VoxelMasker to create a 3D voxel mask");
        const CoordType_ x_min = CoordType_(_bbox_min[0]) / CoordType_(unit_size);
        const CoordType_ x_max = CoordType_(_bbox_max[0]) / CoordType_(unit_size);
        const CoordType_ y_min = CoordType_(_bbox_min[1]) / CoordType_(unit_size);
        const CoordType_ y_max = CoordType_(_bbox_max[1]) / CoordType_(unit_size);
        std::vector<int> line_mask(this->_num_points[0], 0);
        Tiny::Vector<CoordType_, 2> coords;

        for(u64 i = beg; i < end; ++i)
        {
          coords[1] = y_min + (y_max - y_min) * CoordType_(i % this->_num_points[1]) / CoordType_(this->_num_points[1] - 1u);
          masker.mask_line(line_mask, x_min, x_max, coords);
          this->_compress_voxel_map_line(line_mask, i - offset);
        }
      }

      /**
       * \brief Computes a range of voxel map lines for a 3D domain
       *
       * \param[in] masker
       * A \transient reference to the masker that is to be used to compute the mask
       *
       * \param[in] beg
       * The index of the first voxel map line that is to be computed.
       *
       * \param[in] end
       * The index of the last voxel map line that is to be computed plus one.
       *
       * \param[in] offset
       * The offset of the first voxel map line that is to be computed; this is typically equal to 0,
       * unless we are only computing a sub-map that is going to be gathered on the root later on,
       * in which case it may be set equal to \p beg.
       */
      template<typename CoordType_>
      void _compute_voxel_map_lines(VoxelMasker<CoordType_, 3>& masker, u64 beg, u64 end, u64 offset)
      {
        const CoordType_ x_min = CoordType_(_bbox_min[0]) / CoordType_(unit_size);
        const CoordType_ x_max = CoordType_(_bbox_max[0]) / CoordType_(unit_size);
        const CoordType_ y_min = CoordType_(_bbox_min[1]) / CoordType_(unit_size);
        const CoordType_ y_max = CoordType_(_bbox_max[1]) / CoordType_(unit_size);
        const CoordType_ z_min = CoordType_(_bbox_min[2]) / CoordType_(unit_size);
        const CoordType_ z_max = CoordType_(_bbox_max[2]) / CoordType_(unit_size);

        FEAT_PRAGMA_OMP(parallel)
        {
          std::vector<int> line_mask(this->_num_points[0], 0);
          Tiny::Vector<CoordType_, 3> coords;

          FEAT_PRAGMA_OMP(for schedule(dynamic, 16))
          for(u64 i = beg; i < end; ++i)
          {
            // line = iz * this->_num_points[1] + iy
            coords[1] = y_min + (y_max - y_min) * CoordType_(i % this->_num_points[1]) / CoordType_(this->_num_points[1] - 1u);
            coords[2] = z_min + (z_max - z_min) * CoordType_(i / this->_num_points[1]) / CoordType_(this->_num_points[2] - 1u);
            masker.mask_line(line_mask, x_min, x_max, coords);
            this->_compress_voxel_map_line(line_mask, i - offset);
          }
        }
      }

      /**
       * \brief Computes the voxel map for the domain
       *
       * In an MPI-parallel simulation, this function splits the computation of the mask across all processes in the
       * communicator and gathers the mask on all processes afterwards. In the borderline case where there are more
       * processes in the communicator than there are voxel map lines in the mask, this function creates a temporary
       * sub-communicator with as many processes as there are voxel map lines to distribute the work over these
       * processes and performs a broadcast on the original communicator after gathering the mask on rank 0.
       *
       * \param[in] comm
       * A \transient reference to the communicator that is to be used to spread the workload onto.
       *
       * \param[in] masker
       * A \transient reference to the masker that is to be used to compute the mask. Must be the same on all
       * MPI processes
       *
       * \param[in] gather_to_all
       * Specifies whether the voxel map is to be gathered on all processes in the communicator.
       * The only reason to set this to \b false is when the map is only required on rank 0 for file output.
       */
      template<typename CoordType, int dim_>
      void _compute_voxel_map(const Dist::Comm& comm, VoxelMasker<CoordType, dim_>& masker, bool gather_to_all)
      {
        // single process or MPI parallel?
        if(comm.size() <= 1)
        {
          // compute all voxel map lines
          this->_compute_voxel_map_lines(masker, 0u, this->_num_lines, 0u);
        }
        else if(u64(comm.size()) <= this->_num_lines)
        {
          // there are at least as many lines in the voxel map as there are MPI processes; this is the usual case
          // we can use the entire communicator (usually MPI_COMM_WORLD) for the voxel map computation
          u64 comm_rank = u64(comm.rank());
          u64 comm_size = u64(comm.size());

          // compute line count for a single process and round up; this may be bigger than _num_lines,
          // so the last process may have less lines to chew through than all other processes
          u64 line_count = (this->_num_lines + comm_size - 1u) / comm_size;
          std::size_t block_len = this->_stride_line * line_count;

          // allocate voxel map and format to zero
          std::size_t vmap_size = block_len;
          // do we need the entire map on this process?
          if(gather_to_all || (comm_rank == 0))
            vmap_size *= u64(comm_size);
          this->_voxel_map.resize(vmap_size, 0u);

          // compute first and last+1 line of this process; make sure the last process doesn't go beyond the total count
          u64 line_beg = comm_rank * line_count;
          u64 line_end = Math::min((comm_rank + 1u) * line_count, this->_num_lines);

          // debug output
          //this->_comm.allprint(String(">>>") + stringify(line_beg).pad_front(6) + ":" + stringify(line_count) + ">" + stringify(line_end).pad_front(6));

          // compute lines for this process
          this->_compute_voxel_map_lines(masker, line_beg, line_end, gather_to_all ? u64(0) : line_beg);

          // perform an in-place allgather to synchronize the voxel map over all processes
          if(gather_to_all)
            comm.allgather(this->_voxel_map.data(), block_len, this->_voxel_map.data(), block_len);
          else
            comm.gather(this->_voxel_map.data(), block_len, this->_voxel_map.data(), block_len, 0);
        }
        else if(this->_num_lines < u64(comm.size()))
        {
          // more processes than slag lines; this is a borderline scenario, but we have to support it anyways
          // create a sub-communicator and then compute mask on that sub-communicator
          Dist::Comm sub_comm = comm.comm_create_range_incl(int(this->_num_lines));

          // does this process participate in the sub-communicator?
          if(!sub_comm.is_null())
          {
            // get sub-communicator rank and size
            u64 comm_rank = u64(sub_comm.rank());
            u64 comm_size = u64(sub_comm.size());

            // compute line count for a single process and round up; this may be bigger than _num_lines,
            // so the last process may have less lines to chew through than all other processes
            u64 line_count = (this->_num_lines + comm_size - 1u) / comm_size;
            std::size_t block_len = this->_stride_line * line_count;

            // allocate voxel map and format to zero
            std::size_t vmap_size = block_len;
            // do we need the entire map on this process?
            if(gather_to_all || (comm_rank == 0))
              vmap_size *= u64(comm_size);
            this->_voxel_map.resize(vmap_size, 0u);

            // compute first and last+1 line of this process; make sure the last process doesn't go beyond the total count
            u64 line_beg = comm_rank * line_count;
            u64 line_end = Math::min((comm_rank + 1u) * line_count, this->_num_lines);

            // debug output
            //this->_comm.allprint(String(">>>") + stringify(line_beg).pad_front(6) + ":" + stringify(line_count) + ">" + stringify(line_end).pad_front(6));

            // compute lines for this process
            this->_compute_voxel_map_lines(masker, line_beg, line_end, gather_to_all ? u64(0) : line_beg);

            // perform an in-place gather to synchronize the voxel map on rank 0
            sub_comm.gather(this->_voxel_map.data(), block_len, this->_voxel_map.data(), block_len, 0);
          }
          else
          {
            // just allocate the voxel map vector for the following broadcast
            if(gather_to_all)
              this->_voxel_map.resize(this->_stride_line * this->_num_lines, 0u);
          }

          // broadcast from rank 0 to all processes our entire domain communicator
          if(gather_to_all)
            comm.bcast(this->_voxel_map.data(), this->_stride_line * this->_num_lines, 0);
        }

        // that's it
      }

      /**
       * \brief Exports a single plane of the voxel map to a monochrome BMP image file
       *
       * This function will export an entire plane of the voxel map in full resolution as a
       * monochrome bitmap files.
       *
       * \param[in] filename
       * The filename for the BMP image file
       *
       * \param[in] plane
       * The index of the plane that is to be exported.
       */
      void _export_plane_to_bmp(std::ostream& os, u64 plane) const;

      /**
       * \brief Renders a plane range of voxel map to a single gray-scale BMP image file
       *
       * This function will render a range of plane of the voxel map in reduced resolution as a
       * single gray-scale bitmap file, where each pixel contains the average value of all voxels that
       * have been reduced to it.
       *
       * \param[in] os
       * The binary output stream for the BMP image
       *
       * \param[in] width
       * The desired bitmap width; must be less than or equal to the number of points in X-dimension.
       *
       * \param[in] height
       * The desired bitmap height; must be less than or equal to the number of points in Y-dimension.
       *
       * \param[in] box_min_z, box_max_z
       * The indices of the plane Z range that has to be sampled into a single output plane.
       */
      void _render_plane_to_bmp(std::ostream& os, Index width, Index height, i64 box_min_z, i64 box_max_z) const;
    }; // class VoxelMap
  } // namespace Geometry
} // namespace FEAT
