#pragma once
#ifndef KERNEL_ADJACENCY_BASE_HPP
#define KERNEL_ADJACENCY_BASE_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>

namespace FEAST
{
  /**
   * \brief Adjacency namespace
   *
   * This namespace encapsulates classes and class templates related to handling adjacency information.
   */
  namespace Adjacency
  {
    /**
     * \brief Render type enumeration
     *
     * This enumeration specifies the different render modes available for the render constructors.
     */
    enum RenderType
    {
      /**
       * \brief Render-As-Is mode
       *
       * In this mode, the adjactor passed to the constructor will be rendered "as is", i.e. the
       * graph's adjacency information will be identical to the adjactor's, including adjacency
       * duplicates.
       */
      rt_as_is = 0,

      /**
       * \brief Render-Injectified mode
       *
       * In this mode, the adjactor passed to the constructor will be rendered "injective", i.e.
       * each domain node will contain at most one adjacency to a specific image node - in other
       * words: the graph will not contain adjacency duplicates.
       */
      rt_injectify = 1,

      /**
       * \brief Render-Transpose mode
       *
       * In this mode, the transpose of the adjactor passed to the constructor will be rendered
       * instead of the adjactor itself. The render process is performed "as is", i.e. if a domain
       * node \e D has \e k adjacencies with an image node \e I, then the graph's domain node \e I
       * will have \e k adjacencies with the image node \e D.
       */
      rt_transpose = 2,

      /**
       * \brief Render-Injectified-Transpose mode
       *
       * In this mode, the transpose of the adjactor passed to the constructor will be rendered
       * "injective".
       * \see rt_transpose, rt_injectify
       */
      rt_injectify_transpose = 3
    }; // enum RenderType
  } // namespace Adjacency
} // namespace FEAST

#endif // KERNEL_ADJACENCY_BASE_HPP
