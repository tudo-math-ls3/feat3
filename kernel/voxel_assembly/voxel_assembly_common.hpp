#pragma once
#ifndef FEAT_KERNEL_VOXEL_ASSEMBLY_COMMON_HPP
#define FEAT_KERNEL_VOXEL_ASSEMBLY_COMMON_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/shape.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange2/element.hpp>

#ifdef FEAT_HAVE_OMP
#include "omp.h"
#endif

namespace FEAT
{
  /**
   * \brief Namespace for different voxel based assembly methods
   *
   * This namespace incoperates a number of function and classes handling
   * thread parallel assembly of finite element matrices.
   * Since the threading strategy is based on coloring of the elements,
   * a good (i.e. balanced) coloring is required, which is for example given
   * by a voxel based mesh.
   *
   * \author Maximilian Esser
   */
  namespace VoxelAssembly
  {
    /// Q2 quadliteral space
    typedef Space::Lagrange2::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2>>>> Q2StandardQuad;
    /// Q2 hexadron space
    typedef Space::Lagrange2::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<3>>>> Q2StandardHexa;

    /// Templated hypercube definition for standard mapping
    template<int dim_>
    using Q2StandardHyperCube = Space::Lagrange2::Element<Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<dim_>>>>;


    /**
     * \brief A data field for all necessary values that define the dof mapping for assembly.
     *
     * \tparam DT_ The datatype.
     * \tparam IT_ The indextype.
     */
    template<typename DT_, typename IT_>
    struct AssemblyMappingData
    {
      /// The cell to dof, where cell_to_dof[i],..., cell_to_dof[i+cell_dofs-1] are the dofs of one cell
      const IT_* cell_to_dof;
      /// Array of sortingindices of cell_to_dof
      const IT_* cell_to_dof_sorter;
      /// The number of cells
      Index cell_num;
      /// An array of the nodes fitting to the cell_to_dof mapping
      const void* nodes;
      /// The number of nodes
      Index node_size;
    };

    /**
     * \brief A data field for a cubature rule
     *
     * \tparam DT_ The dataype.
     */
    template<typename DT_>
    struct AssemblyCubatureData
    {
      /// The cubature point data array
      const void* cub_pt;
      /// The cubature weights
      const DT_* cub_wg;
      /// Number of cubtaure points
      int num_cubs;
    };

    /**
     * \brief CSR Matrix data
     *
     * \tparam DT_ The dataype.
     * \tparam IT_ The indextype.
     */
    template<typename DT_, typename IT_>
    struct CSRMatrixData
    {
      DT_* data;
      const IT_* row_ptr;
      const IT_* col_idx;
      Index num_rows;
      Index num_cols;
    };

    /**
     * \brief Data for burgers assembler
     *
     * \tparam DT_ The datatype.
     */
    template<typename DT_>
    struct AssemblyBurgersData
    {
      DT_ nu;
      DT_ theta;
      DT_ beta;
      DT_ frechet_beta;
      DT_ sd_delta;
      DT_ sd_nu;
      DT_ sd_v_norm;
      bool deformation;
    };
  }
}
#endif