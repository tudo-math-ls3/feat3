// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2024 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>

namespace FEAT
{
  namespace Control
  {
    namespace Asm
    {
      /**
       * \brief Assembles a scalar homogeneous unit filter for a set of mesh-parts given by their names
       *
       *       *
       * \param[in] dom_level
       * A \transient reference to the domain level on which is to be assembled
       *
       * \param[in] mesh_part_names
       * A string containing all mesh part names on which the unit filter is to be assembled.
       * If set to "*", all mesh-parts with be added to the filter.
       */
      template<typename DomainLevel_>
      Assembly::UnitFilterAssembler<typename DomainLevel_::MeshType> create_unit_filter_asm(const DomainLevel_& dom_level, const std::deque<String>& mesh_part_names)
      {

        Assembly::UnitFilterAssembler<typename DomainLevel_::MeshType> unit_asm;
        // loop over all mesh parts
        for(const auto& mp_name : mesh_part_names)
        {
          auto* mesh_part_node = dom_level.get_mesh_node()->find_mesh_part_node(mp_name);
          XASSERT(mesh_part_node != nullptr);

          // let's see if we have that mesh part
          // if it is nullptr, then our patch is not adjacent to that boundary part
          auto* mesh_part = mesh_part_node->get_mesh();
          if (mesh_part != nullptr)
          {
            // add to boundary assembler
            unit_asm.add_mesh_part(*mesh_part);
          }
        }
        return unit_asm;
      }


      /**
       * \brief Assembles a scalar homogeneous unit filter for a set of mesh-parts given by their names
       *
       *
       * \param[in] mesh_part_names
       * A string containing all mesh part names on which the unit filter is to be assembled.
       * If set to "*", all mesh-parts with be added to the filter.
       */
      template<typename DomainLevel_>
       Assembly::UnitFilterAssembler<typename DomainLevel_::MeshType> create_unit_filter_asm(const DomainLevel_& dom_level, const String& mesh_part_names)
      {
        std::deque<String> mp_names;
        if(mesh_part_names == "*")
          mp_names = dom_level.get_mesh_node()->get_mesh_part_names(true);
        else
          mp_names = mesh_part_names.split_by_whitespaces();

        return create_unit_filter_asm(dom_level, mp_names);
      }

      /**
       * \brief Assembles a scalar homogeneous unit filter for a set of mesh-parts given by their names
       *
       * \param[inout] filter
       * A \transient reference to the filter to be assembled
       *
       * \param[in] dom_level
       * A \transient reference to the domain level on which is to be assembled
       *
       * \param[in] space
       * A \transient reference to the space to be assemble on
       *
       * \param[in] mesh_part_names
       * A string containing all mesh part names on which the unit filter is to be assembled.
       * If set to "*", all mesh-parts with be added to the filter.
       */
      template<typename DataType_, typename IndexType_, typename DomainLevel_, typename Space_>
      void asm_unit_filter_scalar_homogeneous(LAFEM::UnitFilter<DataType_, IndexType_>& filter,
        const DomainLevel_& dom_level, const Space_& space, const String& mesh_part_names)
      {
        Assembly::UnitFilterAssembler<typename DomainLevel_::MeshType> unit_asm = create_unit_filter_asm(dom_level, mesh_part_names);
        unit_asm.assemble(filter, space);
      }

      /**
       * \brief Assembles a scalar unit filter for a set of mesh-parts given by their names
       *
       * \param[inout] filter
       * A \transient reference to the filter to be assembled
       *
       * \param[in] dom_level
       * A \transient reference to the domain level on which is to be assembled
       *
       * \param[in] space
       * A \transient reference to the space to be assemble on
       *
       * \param[in] mesh_part_names
       * A string containing all mesh part names on which the unit filter is to be assembled.
       * If set to "*", all mesh-parts with be added to the filter.
       *
       * \param[in] function
       * A \transient reference to the analytic function that defines the boundary values
       */
      template<typename DataType_, typename IndexType_, typename DomainLevel_, typename Space_, typename Function_>
      void asm_unit_filter_scalar(LAFEM::UnitFilter<DataType_, IndexType_>& filter,
        const DomainLevel_& dom_level, const Space_& space, const String& mesh_part_names, const Function_& function)
      {
        Assembly::UnitFilterAssembler<typename DomainLevel_::MeshType> unit_asm = create_unit_filter_asm(dom_level, mesh_part_names);
        unit_asm.assemble(filter, space, function);
      }

      /**
       * \brief Assembles a blocked homogeneous unit filter for a set of mesh-parts given by their names
       *
       * \param[inout] filter
       * A \transient reference to the filter to be assembled
       *
       * \param[in] dom_level
       * A \transient reference to the domain level on which is to be assembled
       *
       * \param[in] space
       * A \transient reference to the space to be assemble on
       *
       * \param[in] mesh_part_names
       * A string containing all mesh part names on which the unit filter is to be assembled.
       * If set to "*", all mesh-parts with be added to the filter.
       */
      template<typename DataType_, typename IndexType_, int block_size_, typename DomainLevel_, typename Space_>
      void asm_unit_filter_blocked_homogeneous(LAFEM::UnitFilterBlocked<DataType_, IndexType_, block_size_>& filter,
        const DomainLevel_& dom_level, const Space_& space, const String& mesh_part_names)
      {
        Assembly::UnitFilterAssembler<typename DomainLevel_::MeshType> unit_asm = create_unit_filter_asm(dom_level, mesh_part_names);
        unit_asm.assemble(filter, space);
      }

      /**
      * \brief Assembles a scalar unit filter for a set of mesh-parts given by their names
      *
      * \param[inout] filter
      * A \transient reference to the filter to be assembled
      *
      * \param[in] dom_level
      * A \transient reference to the domain level on which is to be assembled
      *
      * \param[in] space
      * A \transient reference to the space to be assemble on
      *
      * \param[in] mesh_part_names
      * A string containing all mesh part names on which the unit filter is to be assembled.
      * If set to "*", all mesh-parts with be added to the filter.
      *
      * \param[in] function
      * A \transient reference to the analytic function that defines the boundary values
      */
      template<typename DataType_, typename IndexType_, int block_size_, typename DomainLevel_, typename Space_, typename Function_>
      void asm_unit_filter_blocked(LAFEM::UnitFilterBlocked<DataType_, IndexType_, block_size_>& filter,
        const DomainLevel_& dom_level, const Space_& space, const String& mesh_part_names, const Function_& function)
      {
        Assembly::UnitFilterAssembler<typename DomainLevel_::MeshType> unit_asm = create_unit_filter_asm(dom_level, mesh_part_names);
        unit_asm.assemble(filter, space, function);
      }
    } // namespace Asm
  } // namespace Control
} // namespace FEAT
