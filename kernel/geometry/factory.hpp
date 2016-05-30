#pragma once
#ifndef KERNEL_GEOMETRY_FACTORY_HPP
#define KERNEL_GEOMETRY_FACTORY_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>

namespace FEAT
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
      };
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_FACTORY_HPP
