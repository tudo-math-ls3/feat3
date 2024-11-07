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

/**
 * \brief CCND Simple namespace
 *
 * This namespaces encapsulates all classes and functions related to the simple CCND applications.
 */
namespace CCNDSimple
{
  /// we're using FEAT
  using namespace FEAT;

  /// the dimension that this app is being compiled for
#ifdef FEAT_CCND_SIMPLE_3D
  static constexpr int dim = 3;
#else
  static constexpr int dim = 2;
#endif

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
  inline void print_pad(const Dist::Comm& comm, const String& key, const String& value, int pad_len = 30, char pad_char = '.')
  {
    comm.print(key.pad_back(std::size_t(pad_len), pad_char) + ": " + value);
  }

  /**
   * \brief Prints memory usage and timing statistics
   *
   * \param[in] comm
   * The communicator
   *
   * \param[in] watch_total
   * The stop watch for the total runtime
   */
  inline void print_final_stats(const Dist::Comm& comm, const StopWatch& watch_total)
  {
    // collect and print memory usage
    MemoryUsage mi;
    std::size_t mem_use[2] = {mi.get_peak_physical(), mi.get_peak_virtual()};
    comm.allreduce(mem_use, mem_use, std::size_t(2), Dist::op_sum);

    comm.print(String("\n") + String(100u, '=') + "\n");
    print_pad(comm, "Peak Physical Memory", stringify_fp_fix(double(mem_use[0]) * 9.313225746154785E-10, 3, 10) + " GiB");
    print_pad(comm, "Peak Virtual Memory",  stringify_fp_fix(double(mem_use[1]) * 9.313225746154785E-10, 3, 10) + " GiB");
    print_pad(comm, "Total Runtime", watch_total.elapsed_string().pad_front(10) + " seconds [" + watch_total.elapsed_string(TimeFormat::h_m_s_m) + "]");
  }
} // namespace CCNDSimple
