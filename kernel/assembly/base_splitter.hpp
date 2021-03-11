// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ASSEMBLY_BASE_SPLITTER_HPP
#define KERNEL_ASSEMBLY_BASE_SPLITTER_HPP 1

#include <kernel/assembly/mirror_assembler.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/vector_mirror.hpp>

#include <vector>

namespace FEAT
{
  namespace Assembly
  {
    /**
     * \brief Base-Mesh Splitter class template.
     *
     * This class is used for splitting a solution vector on a refined mesh
     * into parts associated with the base-mesh cells, as well as for joining
     * those parts back into a single solution vector.
     * This class is required for dynamic load balancing.
     *
     * \author Peter Zajac
     */
    template<typename Space_, typename DT_ = Real, typename IT_ = Index>
    class BaseSplitter
    {
    public:
      typedef Space_ SpaceType;
      typedef typename SpaceType::MeshType MeshType;
      typedef Geometry::RootMeshNode<MeshType> RootMeshNodeType;

    protected:
      std::vector<LAFEM::VectorMirror<Mem::Main, DT_, IT_>> _mirrors;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] space
       * The finite-element space whose vectors are to be split or joined.
       *
       * \param[in] mesh_node
       * The root mesh node on which \p space is defined. The 'create_base_splitting'
       * function must have been called for the root mesh node representing the base
       * mesh from which this mesh node is refined.
       */
      explicit BaseSplitter(const SpaceType& space, const RootMeshNodeType& mesh_node)
      {
        // ensure that the mesh is valid
        XASSERTM(mesh_node.get_mesh() == &space.get_mesh(), "Space mesh and root mesh different");

        // loop over all base cell mesh-parts
        for(Index cell(0); ; ++cell)
        {
          // try to fetch base-cell mesh part
          const auto* mesh_part = mesh_node.find_mesh_part("_base:" + stringify(cell));
          if(mesh_part == nullptr)
            break;

          // assemble mirror
          LAFEM::VectorMirror<Mem::Main, DT_, IT_> mir;
          Assembly::MirrorAssembler::assemble_mirror(mir, space, *mesh_part);
          _mirrors.push_back(std::move(mir));
        }

        XASSERTM(!_mirrors.empty(), "No base cell splitting available in mesh node");
      }

      virtual ~BaseSplitter()
      {
      }

      /**
       * \brief Splits a FE vector into its base-cell parts.
       *
       * \param[out] splits
       * The base-cell parts of the FE vector.
       *
       * \param[in] vector
       * The vector that is to be split.
       */
      void split(
        std::vector<LAFEM::DenseVector<Mem::Main, DT_, IT_>>& splits,
        const LAFEM::DenseVector<Mem::Main, DT_, IT_>& vector) const
      {
        // clear splits
        splits.clear();
        splits.reserve(_mirrors.size());

        // loop over all mirrors and gather
        for(auto it = _mirrors.begin(); it != _mirrors.end(); ++it)
        {
          auto svec = (*it).create_buffer(vector);
          (*it).gather(svec, vector);
          splits.push_back(std::move(svec));
        }
      }

      /**
       * \brief Joins a FE vector from its base-cell parts.
       *
       * \param[in,out] vector
       * The joined vector. Is assumed to be allocated to correct length.
       *
       * \param[in] splits
       * The base-cell parts of the FE vector.
       */
      void join(
        LAFEM::DenseVector<Mem::Main, DT_, IT_>& vector,
        const std::vector<LAFEM::DenseVector<Mem::Main, DT_, IT_>>& splits) const
      {
        // check size
        XASSERTM(_mirrors.size() == splits.size(), "Invalid number of split vectors");

        // format output
        vector.format();

        // clone for a frequencies vector
        auto freqs = vector.clone();

        // loop over all mirrors and scatter splits
        for(std::size_t i(0); i < _mirrors.size(); ++i)
        {
          _mirrors.at(i).scatter_axpy(vector, splits.at(i));
          auto tmp = splits.at(i).clone();
          tmp.format(DT_(1));
          _mirrors.at(i).scatter_axpy(freqs, tmp);
        }

        // scale output
        freqs.component_invert(freqs);
        vector.component_product(vector, freqs);
      }
    }; // class BaseSplitter<...>
  } // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_BASE_SPLITTER_HPP
