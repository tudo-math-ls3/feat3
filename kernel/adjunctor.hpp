#pragma once
#ifndef KERNEL_ADJUNCTOR_HPP
#define KERNEL_ADJUNCTOR_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>

namespace FEAST
{
  // Note: The following class is just an abstract description of an adjunctor for the documentation -- it is not
  // an actual class that should be compiled, therefore the class is inside a block reserved only for doxygen.
#ifdef DOXYGEN
  /**
   * \brief Adjunctor interface
   *
   * \todo details and terminology
   *
   * \author Peter Zajac
   */
  class Adjunctor
  {
  public:
    /**
     * \brief Adjunctor dual node iterator class.
     *
     * \todo details
     *
     * \warning Please note that an implementation of the DualIterator interface does not necessarily provide
     * a \em default constructor. Take this into account when you develop an algorithm based on adjunctors.
     *
     * \author Peter Zajac
     */
    class DualIterator
    {
    public:
      /**
       * \brief Default constructor.
       */
      DualIterator();

      /**
       * \brief Copy constructor.
       *
       * \param[in] other
       * The dual node iterator that is to be copied.
       */
      DualIterator(const DualIterator& other);

      /**
       * \brief Dereferentiation operator.
       *
       * \returns
       * The index of the dual node the iterator currently refers to.
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
      DualIterator& operator++();

      /**
       * \brief Checks two iterators for inequality.
       *
       * \param[in] other
       * The iterator to be checked against.
       *
       * \returns
       * \c false if the iterators are identical, otherwise \c true.
       */
      bool operator!=(const DualIterator& other) const;
    }; // class DualIterator

    /**
     * \brief Returns the number of primal nodes.
     *
     * \returns
     * The total number of primal nodes in the adjunctor.
     */
    Index num_nodes_primal() const;

    /**
     * \brief Returns the number of dual nodes.
     *
     * \returns
     * The total number of dual nodes in the adjunctor.
     */
    Index num_nodes_dual() const;

    /**
     * \brief Returns an iterator for the first adjacent dual node.
     *
     * \param[in] primal_node
     * The index of the primal node whose dual node iterator is to be returned.
     *
     * \returns
     * An iterator representing the first dual node adjacent to \p primal_node, or dual_end(primal_node) if no
     * adjacent dual nodes exist.
     */
    DualIterator dual_begin(Index primal_node) const;

    /**
     * \brief Returns an iterator for the first position past the last adjacent dual node.
     *
     * \param[in] primal_node
     * The index of the primal node whose dual node iterator is to be returned.
     *
     * \returns
     * An iterator representing the first position past the last dual node adjacent to \p primal_node.
     */
    DualIterator dual_end(Index primal_node) const;
  }; // class Adjcuntor
#endif // DOXYGEN

  /**
   * \brief Composite Adjunctor implementation
   *
   * \todo detailed description
   *
   * \author Peter Zajac
   */
  template<
    typename Adj1_,
    typename Adj2_>
  class CompositeAdjunctor
  {
  private:
    /// first adjunctor
    const Adj1_& _adj1;
    /// second adjunctor
    const Adj2_& _adj2;

  public:
    /// typedef for first adjunctor
    typedef Adj1_ Adjunctor1;
    /// typedef for second adjunctor
    typedef Adj2_ Adjunctor2;

    /**
     * \brief Dual Iterator implementation for CompositeAdjunctor
     *
     * \author Peter Zajac
     */
    class DualIterator
    {
      // The CompositeAdjunctor class needs to be a friend as it calls the protected constuctors.
      friend class CompositeAdjunctor<Adj1_, Adj2_>;

    private:

      /// Dual iterator of first adjunctor
      typedef typename Adj1_::DualIterator Iterator1;
      /// Dual iterator of second adjunctor
      typedef typename Adj2_::DualIterator Iterator2;

      /// pointer to second adjunctor
      const Adj2_* _adj2;

      /// current iterator for first adjunctor
      mutable Iterator1 _cur1;
      /// end iterator for first adjunctor
      mutable Iterator1 _end1;
      /// current iterator for second adjunctor
      mutable Iterator2 _cur2;
      /// end iterator for second adjunctor
      mutable Iterator2 _end2;

    protected:
      /// constructor for dual_begin()
      DualIterator(
        const Adj1_* adj1,
        const Adj2_* adj2,
        Index primal_node)
         :
        _adj2(adj2),
        _cur1(adj1->dual_begin(primal_node)),
        _end1(adj1->dual_end(primal_node)),
        _cur2(),
        _end2()
      {
        CONTEXT("CompositeAdjunctor::DualIterator:DualIterator()");
        if(_cur1 != _end1)
        {
          _cur2 = adj2->dual_begin(*_cur1);
          _end2 = adj2->dual_end(*_cur1);
        }
      }

      /// constructor for dual_end()
      DualIterator(
        const Adj1_* adj1,
        Index primal_node)
         :
        _adj2(nullptr),
        _cur1(adj1->dual_end(primal_node)),
        _end1(_cur1),
        _cur2(),
        _end2()
      {
        CONTEXT("CompositeAdjunctor::DualIterator:DualIterator()");
      }

    public:
      inline DualIterator() :
        _adj2(nullptr),
        _cur1(),
        _end1(),
        _cur2(),
        _end2()
      {
      }

      /** \copydoc Adjunctor::DualIterator::DualIterator() */
      inline DualIterator(const DualIterator& other) :
        _adj2(other._adj2),
        _cur1(other._cur1),
        _cur2(other._cur2),
        _end1(other._end1),
        _end2(other._end2)
      {
        CONTEXT("CompositeAdjunctor::DualIterator:DualIterator()");
      }

      /** \copydoc Adjunctor::DualIterator::operator*() */
      inline Index operator*() const
      {
        CONTEXT("CompositeAdjunctor::DualIterator::operator*()");
        return *_cur2;
      }

      /** \copydoc Adjunctor::DualIterator::operator!=() */
      inline bool operator!=(const DualIterator& other) const
      {
        CONTEXT("CompositeAdjunctor::DualIterator::operator!=()");
        return (_cur1 != other._cur1) || (_cur2 != other._cur2);
      }

      /** \copydoc Adjunctor::DualIterator::operator++() */
      DualIterator& operator++()
      {
        CONTEXT("CompositeAdjunctor::DualIterator::operator++()");

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
          _cur2 = _adj2->dual_begin(*_cur1);
          _end2 = _adj2->dual_end(*_cur1);

          // check whether the second adjacency list is empty, if not then we have a valid position
          if(_cur2 != _end2)
          {
            return *this;
          }
        }

        // If we come out here, then there are no more adjacencies left, so we need to set
        // this iterator to the 'end' position - this setting must be compatible to the iterator
        // returned by the dual_end() function of the CompositeAdjunctor class!
        // Note: see the second protected constructor of this class for the following definition.
        _cur2 = Iterator2();
        _end2 = _cur2;

        return *this;
      }
    }; // class CompositeAdjunctor::DualIterator

    /**
     * \brief Constructor.
     */
    inline CompositeAdjunctor(
      const Adj1_& adjunctor1,
      const Adj2_& adjunctor2)
       :
      _adj1(adjunctor1),
      _adj2(adjunctor2)
    {
      CONTEXT("CompositeAdjunctor::CompositeAdjunctor()");
      // Ensure that the number of dual nodes in the first adjunctor is not greater than the number of primal
      // nodes in the second adjunctor; otherwise the composite adjunctor is ill-formed.
      ASSERT(_adj1.num_nodes_dual() <= _adj2.num_nodes_primal(), "Composite Adjuctor is ill-formed");
    }

    /** \copydoc Adjunctor::num_nodes_primal() */
    inline Index num_nodes_primal() const
    {
      CONTEXT("CompositeAdjunctor::num_nodes_primal()");
      // The number of primal nodes of the composite adjunctor is given by the number of primal nodes of the
      // first adjunctor in the composition.
      return _adj1.num_nodes_primal();
    }

    /** \copydoc Adjunctor::num_nodes_dual() */
    inline Index num_nodes_dual() const
    {
      CONTEXT("CompositeAdjunctor::num_nodes_dual()");
      // The number of dual nodes of the composite adjunctor is given by the number of dual nodes of the
      // second adjunctor in the composition.
      return _adj2.num_nodes_dual();
    }

    /** \copydoc Adjunctor::dual_begin() */
    inline DualIterator dual_begin(Index primal_node) const
    {
      CONTEXT("CompositeAdjuctor::dual_begin()");
      ASSERT(primal_node < num_nodes_primal(), "primal_node out of valid range");
      return DualIterator(&_adj1, &_adj2, primal_node);
    }

    /** \copydoc Adjunctor::dual_end() */
    inline DualIterator dual_end(Index primal_node) const
    {
      CONTEXT("CompositeAdjuctor::dual_end()");
      ASSERT(primal_node < num_nodes_primal(), "primal_node out of valid range");
      return DualIterator(&_adj1, primal_node);
    }

    /**
     * \brief Returns the first adjunctor
     * \returns
     * A reference to the first adjunctor.
     */
    inline const Adjunctor1& get_adjunctor1() const
    {
      return _adj1;
    }

    /**
     * \brief Returns the second adjunctor.
     * \returns
     * A reference to the second adjunctor.
     */
    inline const Adjunctor2& get_adjunctor2() const
    {
      return _adj2;
    }
  }; // class CompositeAdjunctor
} // namespace FEAST

#endif // KERNEL_ADJUNCTOR_HPP
