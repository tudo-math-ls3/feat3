#pragma once
#ifndef KERNEL_LAFEM_CONTAINER_POOL_HPP
#define KERNEL_LAFEM_CONTAINER_POOL_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/instantiation_policy.hpp>
#include <kernel/archs.hpp>

#include <vector>


namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Container Pool.
     *
     * This class provides a pool of often used container instances.
     *
     * \author Dirk Ribbrock
     */
    template <typename CT_>
    class ContainerPool
        : public InstantiationPolicy<ContainerPool<CT_>, Singleton>
    {
      private:
        std::vector<CT_> _list;

      public:
        Index size() const
        {
          return _list.size();
        }

        CT_ & at(Index i)
        {
          return _list.at(i);
        }

        void push_back(CT_ & temp)
        {
          _list.push_back(temp);
        }
    };
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_CONTAINER_POOL_HPP
