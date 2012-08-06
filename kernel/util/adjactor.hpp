#pragma once
#ifndef KERNEL_UTIL_ADJACTOR_HPP
#define KERNEL_UTIL_ADJACTOR_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>

namespace FEAST
{
  // Note: The following class is just an abstract description of an adjactor for the documentation -- it is not
  // an actual class that should be compiled, therefore the class is inside a block reserved only for doxygen.
#ifdef DOXYGEN
  /**
   * \brief Adjactor interface
   *
   * \todo details and terminology
   *
   * \author Peter Zajac
   */
  class Adjactor
  {
  public:
    /**
     * \brief Adjactor image node iterator class.
     *
     * \todo details
     *
     * \author Peter Zajac
     */
    class ImageIterator
    {
    public:
      /**
       * \brief Standard constructor.
       */
      ImageIterator();

      /**
       * \brief Copy constructor.
       *
       * \param[in] other
       * The image node iterator that is to be copied.
       */
      ImageIterator(const ImageIterator& other);

      /**
       * \brief Assignment operator.
       *
       * \param[in] other
       * The image node iteration that is to be assigned.
       *
       * \returns <c>*this</c>
       */
      ImageIterator& operator=(const ImageIterator& other);

      /**
       * \brief Dereferentiation operator.
       *
       * \returns
       * The index of the image node the iterator currently refers to.
       */
      Index operator*() const;

      /**
       * \brief Pre-increment operator.
       *
       * Increments the iterator by one position.
       *
       * \returns
       * \c *this
       */
      ImageIterator& operator++();

      /**
       * \brief Checks two iterators for inequality.
       *
       * \param[in] other
       * The iterator to be checked against.
       *
       * \returns
       * \c false if the iterators are identical, otherwise \c true.
       */
      bool operator!=(const ImageIterator& other) const;
    }; // class ImageIterator

    /**
     * \brief Returns the number of domain nodes.
     *
     * \returns
     * The total number of domain nodes in the adjactor.
     */
    Index get_num_nodes_domain() const;

    /**
     * \brief Returns the number of image nodes.
     *
     * \returns
     * The total number of image nodes in the adjactor.
     */
    Index get_num_nodes_image() const;

    /**
     * \brief Returns an iterator for the first adjacent image node.
     *
     * \param[in] domain_node
     * The index of the domain node whose image node iterator is to be returned.
     *
     * \returns
     * An iterator representing the first image node adjacent to \p domain_node, or image_end(domain_node) if no
     * adjacent image nodes exist.
     */
    ImageIterator image_begin(Index domain_node) const;

    /**
     * \brief Returns an iterator for the first position past the last adjacent image node.
     *
     * \param[in] domain_node
     * The index of the domain node whose image node iterator is to be returned.
     *
     * \returns
     * An iterator representing the first position past the last image node adjacent to \p domain_node.
     */
    ImageIterator image_end(Index domain_node) const;
  }; // class Adjcuntor
#endif // DOXYGEN

  /**
   * \brief Composite Adjactor implementation
   *
   * \todo detailed description
   *
   * \author Peter Zajac
   */
  template<
    typename Adj1_,
    typename Adj2_>
  class CompositeAdjactor
  {
  private:
    /// first adjactor
    const Adj1_& _adj1;
    /// second adjactor
    const Adj2_& _adj2;

  public:
    /// typedef for first adjactor
    typedef Adj1_ Adjactor1;
    /// typedef for second adjactor
    typedef Adj2_ Adjactor2;

    /**
     * \brief Dual Iterator implementation for CompositeAdjactor
     *
     * \author Peter Zajac
     */
    class ImageIterator
    {
      // The CompositeAdjactor class needs to be a friend as it calls the protected constuctors.
      friend class CompositeAdjactor<Adj1_, Adj2_>;

    private:

      /// Dual iterator of first adjactor
      typedef typename Adj1_::ImageIterator Iterator1;
      /// Dual iterator of second adjactor
      typedef typename Adj2_::ImageIterator Iterator2;

      /// pointer to second adjactor
      const Adj2_* _adj2;

      /// current iterator for first adjactor
      mutable Iterator1 _cur1;
      /// end iterator for first adjactor
      mutable Iterator1 _end1;
      /// current iterator for second adjactor
      mutable Iterator2 _cur2;
      /// end iterator for second adjactor
      mutable Iterator2 _end2;

    protected:
      /// \cond internal
      // constructor for image_begin()
      ImageIterator(
        const Adj1_* adj1,
        const Adj2_* adj2,
        Index domain_node)
         :
        _adj2(adj2),
        _cur1(adj1->image_begin(domain_node)),
        _end1(adj1->image_end(domain_node)),
        _cur2(),
        _end2()
      {
        CONTEXT("CompositeAdjactor::ImageIterator:ImageIterator()");
        if(_cur1 != _end1)
        {
          _cur2 = adj2->image_begin(*_cur1);
          _end2 = adj2->image_end(*_cur1);
        }
      }

      // constructor for image_end()
      ImageIterator(
        const Adj1_* adj1,
        Index domain_node)
         :
        _adj2(nullptr),
        _cur1(adj1->image_end(domain_node)),
        _end1(_cur1),
        _cur2(),
        _end2()
      {
        CONTEXT("CompositeAdjactor::ImageIterator:ImageIterator()");
      }
      /// \endcond

    public:
      /// \brief Standard constructor
      inline ImageIterator() :
        _adj2(nullptr),
        _cur1(),
        _end1(),
        _cur2(),
        _end2()
      {
      }

      /** \copydoc Adjactor::ImageIterator::ImageIterator(const ImageIterator& other) */
      inline ImageIterator(const ImageIterator& other) :
        _adj2(other._adj2),
        _cur1(other._cur1),
        _cur2(other._cur2),
        _end1(other._end1),
        _end2(other._end2)
      {
        CONTEXT("CompositeAdjactor::ImageIterator:ImageIterator()");
      }

      /** \copydoc Adjactor::ImageIterator::operator*() */
      inline Index operator*() const
      {
        CONTEXT("CompositeAdjactor::ImageIterator::operator*()");
        return *_cur2;
      }

      /** \copydoc Adjactor::ImageIterator::operator!=() */
      inline bool operator!=(const ImageIterator& other) const
      {
        CONTEXT("CompositeAdjactor::ImageIterator::operator!=()");
        return (_cur1 != other._cur1) || (_cur2 != other._cur2);
      }

      /** \copydoc Adjactor::ImageIterator::operator++() */
      ImageIterator& operator++()
      {
        CONTEXT("CompositeAdjactor::ImageIterator::operator++()");

        // increment the second iterator; if it did not reach the end then we have a valid position
        if(++_cur2 != _end2)
        {
          return *this;
        }

        // second iterator reached its end, so keep incrementing the first iterator until we find the next
        // non-empty second adjacency list or until we reach its end
        while(++_cur1 != _end1)
        {
          // reset second iterator
          _cur2 = _adj2->image_begin(*_cur1);
          _end2 = _adj2->image_end(*_cur1);

          // check whether the second adjacency list is empty, if not then we have a valid position
          if(_cur2 != _end2)
          {
            return *this;
          }
        }

        // If we come out here, then there are no more adjacencies left, so we need to set
        // this iterator to the 'end' position - this setting must be compatible to the iterator
        // returned by the image_end() function of the CompositeAdjactor class!
        // Note: see the second protected constructor of this class for the following definition.
        _cur2 = Iterator2();
        _end2 = _cur2;

        return *this;
      }
    }; // class CompositeAdjactor::ImageIterator

    /**
     * \brief Constructor.
     *
     * \param[in] adjactor1
     * The first adjactor for the composition.
     *
     * \param[in] adjactor2
     * The second adjactor for the composition.
     */
    inline CompositeAdjactor(
      const Adj1_& adjactor1,
      const Adj2_& adjactor2)
       :
      _adj1(adjactor1),
      _adj2(adjactor2)
    {
      CONTEXT("CompositeAdjactor::CompositeAdjactor()");
      // Ensure that the number of image nodes in the first adjactor is not greater than the number of domain
      // nodes in the second adjactor; otherwise the composite adjactor is ill-formed.
      ASSERT(_adj1.get_num_nodes_image() <= _adj2.get_num_nodes_domain(), "Composite Adjuctor is ill-formed");
    }

    /** \copydoc Adjactor::num_nodes_domain() */
    inline Index get_num_nodes_domain() const
    {
      CONTEXT("CompositeAdjactor::num_nodes_domain()");
      // The number of domain nodes of the composite adjactor is given by the number of domain nodes of the
      // first adjactor in the composition.
      return _adj1.get_num_nodes_domain();
    }

    /** \copydoc Adjactor::num_nodes_image() */
    inline Index get_num_nodes_image() const
    {
      CONTEXT("CompositeAdjactor::num_nodes_image()");
      // The number of image nodes of the composite adjactor is given by the number of image nodes of the
      // second adjactor in the composition.
      return _adj2.get_num_nodes_image();
    }

    /** \copydoc Adjactor::image_begin() */
    inline ImageIterator image_begin(Index domain_node) const
    {
      CONTEXT("CompositeAdjuctor::image_begin()");
      ASSERT(domain_node < get_num_nodes_domain(), "domain_node out of valid range");
      return ImageIterator(&_adj1, &_adj2, domain_node);
    }

    /** \copydoc Adjactor::image_end() */
    inline ImageIterator image_end(Index domain_node) const
    {
      CONTEXT("CompositeAdjuctor::image_end()");
      ASSERT(domain_node < get_num_nodes_domain(), "domain_node out of valid range");
      return ImageIterator(&_adj1, domain_node);
    }

    /**
     * \brief Returns the first adjactor
     * \returns
     * A reference to the first adjactor in the composition.
     */
    inline const Adjactor1& get_adjactor1() const
    {
      return _adj1;
    }

    /**
     * \brief Returns the second adjactor.
     * \returns
     * A reference to the second adjactor in the composition.
     */
    inline const Adjactor2& get_adjactor2() const
    {
      return _adj2;
    }
  }; // class CompositeAdjactor
} // namespace FEAST

#endif // KERNEL_UTIL_ADJACTOR_HPP
