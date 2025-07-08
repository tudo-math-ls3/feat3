// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/stop_watch.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/memory_usage.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/trafo/isoparam/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/analytic/lambda_function.hpp>

/**
 * \brief CCND Simple namespace
 *
 * This namespaces encapsulates all classes and functions related to the simple CCND applications.
 */
namespace CCNDSimple
{
  /// we're using FEAT
  using namespace FEAT;

#ifndef FEAT_CCND_SIMPLE_DIM
#error Macro 'FEAT_CCND_SIMPLE_DIM' must be defined by build system!
#endif

  /// the dimension that this app is being compiled for
  static constexpr int dim = FEAT_CCND_SIMPLE_DIM;

  /// our one and only index type
  typedef FEAT::Index IndexType;
  /// our one and only data type
  typedef double DataType;

  /// the shape type of the mesh
  typedef Shape::Hypercube<dim> ShapeType;

  /// the mesh type
  typedef Geometry::ConformalMesh<ShapeType, dim, DataType> MeshType;

  /// the trafo type
#ifdef FEAT_CCND_SIMPLE_ISOPARAM
  typedef Trafo::Isoparam::Mapping<MeshType, 2> TrafoType;
#else
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
#endif

  /// the velocity space type
  typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;

  /// the pressure space type
  typedef Space::Discontinuous::ElementP1<TrafoType> SpacePresType;

  /// the Q1 space type
  typedef Space::Lagrange1::Element<TrafoType> SpaceTypeQ1;

  /// the Q2 space type
  typedef Space::Lagrange2::Element<TrafoType> SpaceTypeQ2;


  /// local scalar vector type
  typedef LAFEM::DenseVector<DataType, IndexType> LocalScalarVectorType;

  /// local field/blocked vector type
  typedef LAFEM::DenseVectorBlocked<DataType, IndexType, dim> LocalFieldVectorType;

  /**
   * \brief Auxiliary helper function for padded printing
   *
   * \param[in] comm
   * The communicator
   *
   * \param[in] key
   * The key to be printed
   *
   * \param[in] value
   * The value to be printed
   *
   * \param[in] pad_len
   * The padding length
   *
   * \param[in] pad_char
   * The character to pad with.
   */
  inline void print_pad(const Dist::Comm& comm, const String& key, const String& value, int pad_len = 40, char pad_char = '.')
  {
    comm.print(key.pad_back(std::size_t(pad_len), pad_char) + ": " + value);
  }

  /**
   * \brief Auxiliary helper function for padded runtime printing
   *
   * \param[in] comm
   * The communicator
   *
   * \param[in] key
   * The key to be printed
   *
   * \param[in] time
   * The runtime to be printed
   *
   * \param[in] total_runtime
   * The total runtime of the entire simulation
   *
   * \param[in] pad_len
   * The padding length
   *
   * \param[in] pad_char
   * The character to pad with.
   */
  inline void print_time(const Dist::Comm& comm, const String& key, double time, double total_runtime, int pad_len = 40, char pad_char = '.')
  {
    double min_time = time, max_time = time;
    comm.reduce(&min_time, &min_time, std::size_t(1), Dist::op_min, 0);
    comm.reduce(&max_time, &max_time, std::size_t(1), Dist::op_max, 0);
    if(comm.rank() == 0)
    {
      std::cout << key.pad_back(std::size_t(pad_len), pad_char);
      std::cout << ": ";
      std::cout << stringify_fp_fix(max_time, 3, 10);
      std::cout << " sec { ";
      std::cout << stringify_fp_fix(min_time > 1E-5 ? min_time/max_time : 0.0, 3, 5);
      std::cout << " } [";
      std::cout << stringify_fp_fix(100.0*max_time/total_runtime, 2, 6);
      std::cout << "% ]\n";
    }
  }

  /**
   * \brief Prints memory usage statistics
   *
   * \param[in] comm
   * The communicator
   */
  inline void print_memory_usage(const Dist::Comm& comm, int pad_len = 40, char pad_char = '.')
  {
    MemoryUsage mi;
    std::size_t mem_sum[2] = {mi.get_peak_physical(), mi.get_peak_virtual()};
    std::size_t mem_min[2] = {mi.get_peak_physical(), mi.get_peak_virtual()};
    std::size_t mem_max[2] = {mi.get_peak_physical(), mi.get_peak_virtual()};
    comm.reduce(mem_sum, mem_sum, std::size_t(2), Dist::op_sum, 0);
    comm.reduce(mem_min, mem_min, std::size_t(2), Dist::op_min, 0);
    comm.reduce(mem_max, mem_max, std::size_t(2), Dist::op_max, 0);
    print_pad(comm, "Peak Physical Memory Usage per Process", stringify_bytes(mem_max[0], 3, 10)
      + " { " + stringify_fp_fix(mem_min[0] > 0u ? double(mem_min[0])/double(mem_max[0]) : 0.0, 3, 5) + " }", pad_len, pad_char);
    print_pad(comm, "Peak Virtual Memory Usage per Process", stringify_bytes(mem_max[1], 3, 10)
      + " { " + stringify_fp_fix(mem_min[1] > 0u ? double(mem_min[1])/double(mem_max[1]) : 0.0, 3, 5) + " }", pad_len, pad_char);
    print_pad(comm, "Total Peak Physical Memory Usage", stringify_bytes(mem_sum[0], 3, 10), pad_len, pad_char);
    print_pad(comm, "Total Peak Virtual Memory Usage",  stringify_bytes(mem_sum[1], 3, 10), pad_len, pad_char);
  }
} // namespace CCNDSimple
