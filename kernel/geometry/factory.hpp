// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>

// includes, system
#include <array>

namespace FEAT
{
  namespace Geometry
  {
    /**
     * \brief Mesh Factory class template
     *
     * See the specializations of this class template for the actual factory interfaces.
     *
     * \tparam Mesh_
     * The type of the mesh that is to be produced by the factory.
     *
     * \author Peter Zajac
     */
    template<typename Mesh_>
#ifndef DOXYGEN
    class Factory;
#else
    class Factory
    {
    public:
      virtual ~Factory()
      {
      }
    }; // class Factory<...>
#endif // DOXYGEN

    /**
     * \brief Standard Refinery class template
     *
     * See the specializations of this class template for the actual standard-refinery interfaces.
     *
     * \tparam Mesh_
     * The type of the mesh that is to be refined by the refinery. This may be either a Mesh or
     * a MeshPart instance.
     *
     * \author Peter Zajac
     */
    template<typename Mesh_>
#ifndef DOXYGEN
    class StandardRefinery;
#else
    class StandardRefinery :
      public Factory<Mesh_>
    {
    public:
      /// mesh type
      typedef Mesh_ MeshType;

      /**
       * \brief Constructor.
       *
       * \param[in] coarse_mesh
       * A \resident reference to the coarse mesh that is to be refined.
       */
      explicit StandardRefinery(const MeshType& coarse_mesh);
    };
#endif

    /// \cond internal
    namespace Intern
    {
      template<int n_>
      struct NumEntitiesWrapper
      {
        Index num_entities[n_+1];

        template<typename Factory_>
        explicit NumEntitiesWrapper(Factory_& factory)
        {
          for(int i(0); i <= n_; ++i)
          {
            num_entities[i] = factory.get_num_entities(i);
          }
        }

        template<typename Factory_>
        static void apply(Factory_& factory, Index* num_ent)
        {
          for(int i(0); i <= n_; ++i)
          {
            num_ent[i] = factory.get_num_entities(i);
          }
        }

        template<typename Factory_, std::size_t m_>
        static void apply(Factory_& factory, std::array<Index, m_>& num_ent)
        {
          static_assert(int(m_) >= n_+1, "invalid array size");
          for(int i(0); i <= n_; ++i)
          {
            num_ent[std::size_t(i)] = factory.get_num_entities(i);
          }
        }
      };

      template<int n_>
      struct NumSlicesWrapper
      {
        Index num_slices[n_];

        template<typename Factory_>
        explicit NumSlicesWrapper(Factory_& factory)
        {
          for(int i(0); i < n_; ++i)
          {
            num_slices[i] = factory.get_num_slices(i);
          }
        }

        template<typename Factory_>
        static void apply(Factory_& factory, Index* num_slic)
        {
          for(int i(0); i < n_; ++i)
          {
            num_slic[i] = factory.get_num_slices(i);
          }
        }

        template<typename Factory_, std::size_t m_>
        static void apply(Factory_& factory, std::array<Index, m_>& num_slic)
        {
          static_assert(int(m_) >= n_, "invalid array size");
          for(int i(0); i < n_; ++i)
          {
            num_slic[i] = factory.get_num_slices(i);
          }
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT
