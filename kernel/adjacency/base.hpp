#pragma once
#ifndef KERNEL_ADJACENCY_BASE_HPP
#define KERNEL_ADJACENCY_BASE_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>

namespace FEAT
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
    enum class RenderType
    {
      /**
       * \brief Render-As-Is mode
       *
       * In this mode, the adjactor passed to the constructor will be rendered "as is", i.e. the
       * graph's adjacency information will be identical to the adjactor's, including adjacency
       * duplicates.
       */
      as_is = 0,

      /**
       * \brief Render-Injectified mode
       *
       * In this mode, the adjactor passed to the constructor will be rendered "injective", i.e.
       * each domain node will contain at most one adjacency to a specific image node - in other
       * words: the graph will not contain adjacency duplicates.
       */
      injectify = 1,

      /**
       * \brief Render-Transpose mode
       *
       * In this mode, the transpose of the adjactor passed to the constructor will be rendered
       * instead of the adjactor itself. The render process is performed "as is", i.e. if a domain
       * node \e D has \e k adjacencies with an image node \e I, then the graph's domain node \e I
       * will have \e k adjacencies with the image node \e D.
       */
      transpose = 2,

      /**
       * \brief Render-Injectified-Transpose mode
       *
       * In this mode, the transpose of the adjactor passed to the constructor will be rendered
       * "injective".
       * \see RenderType::transpose, RenderType::injectify
       */
      injectify_transpose = 3
    }; // enum class RenderType
  } // namespace Adjacency
} // namespace FEAT

#endif // KERNEL_ADJACENCY_BASE_HPP
