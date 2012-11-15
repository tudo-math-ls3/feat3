#pragma once
#ifndef KERNEL_GEOMETRY_FACTORY_HPP
#define KERNEL_GEOMETRY_FACTORY_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /**
     * \brief Mesh Factory class template
     *
     * See the specialisations of this class template for the actual factory interfaces.
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
     * See the specialisations of this class template for the actual standard-refinery interfaces.
     *
     * \tparam Mesh_
     * The type of the mesh that is to be refined by the refinery.
     *
     * \tparam Parent_
     * The type of the parent mesh or cell subset for a submesh or cell subset.
     * Is set to \c Nil if the mesh to be refined is the root mesh.
     *
     * \author Peter Zajac
     */
    template<
      typename Mesh_,
      typename Parent_ = Nil>
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
       * The coarse mesh that is to be refined.
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
        NumEntitiesWrapper(const Factory_& factory)
        {
          for(int i(0); i <= n_; ++i)
          {
            num_entities[i] = factory.get_num_entities(i);
          }
        }

        template<typename Factory_>
        static void apply(const Factory_& factory, Index* num_ent)
        {
          for(int i(0); i <= n_; ++i)
          {
            num_ent[i] = factory.get_num_entities(i);
          }
        }
      };

      template<int n_>
      struct NumSlicesWrapper
      {
        Index num_slices[n_];

        template<typename Factory_>
        NumSlicesWrapper(const Factory_& factory)
        {
          for(int i(0); i < n_; ++i)
          {
            num_slices[i] = factory.get_num_slices(i);
          }
        }

        template<typename Factory_>
        static void apply(const Factory_& factory, Index* num_slic)
        {
          for(int i(0); i < n_; ++i)
          {
            num_slic[i] = factory.get_num_slices(i);
          }
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_FACTORY_HPP
