#pragma once
#ifndef KERNEL_UTIL_PERMUTATION_HPP
#define KERNEL_UTIL_PERMUTATION_HPP 1

// includes, FEAST
#include <kernel/util/assertion.hpp>
#include <kernel/util/random.hpp>

namespace FEAST
{
  /**
   * \brief Permutation class
   *
   * \author Peter Zajac
   */
  class Permutation
  {
  public:
    /// Construction type enumeration
    enum ConstrType
    {
      /// create uninitialised permutation
      type_none,
      /// create identity permutation
      type_identity,
      /// create from permutation array
      type_perm,
      /// create from swap array
      type_swap,
      /// create from inverse permutation array
      type_inv_perm,
      /// create from inverse swap array
      type_inv_swap
    };

  private:
    /// the size of the permutation
    Index _num_entries;
    /// the permute-position array
    Index* _perm_pos;
    /// the swap-position array
    Index* _swap_pos;

  public:
    /**
     * \brief Constructor
     *
     * \param[in] num_entries
     * The size of the permutation to be created.
     *
     * \param[in] constr_type
     * Specifies the construction type:
     * - \p type_none\n Create an uninitialised permutation.\n
     *   The permutation array has to be set up after the object is created.\n
     *   The input array \p v is ignored.
     * - \p type_identity\n Create an identity permutation.\nThe input array \p v is ignored.
     * - \p type_perm\n Interpret the input array as a permute-position array.
     * - \p type_swap\n Interpret the input array as a swap-position array.
     * - \p type_inv_perm\n Interpret the input array as an inverse permute-position array.
     * - \p type_inv_swap\n Interpret the input array as an inverse swap-position array.
     *
     * \param[in] v
     * The input array for the initialisation. The interpretation of the array's content depends on the
     * \p constr_type parameter.
     */
    explicit Permutation(
      Index num_entries,
      ConstrType constr_type = type_none,
      const Index* v = nullptr)
       :
      _num_entries(num_entries),
      _perm_pos(new Index[num_entries]),
      _swap_pos(new Index[num_entries])
    {
      ASSERT_(num_entries > 0);
      switch(constr_type)
      {
      case type_none:
        // leave arrays uninitialised
        return;

      case type_identity:
        // initialise identity permutation
        for(Index i(0); i < _num_entries; ++i)
        {
          _perm_pos[i] = _swap_pos[i] = i;
        }
        return;

      default:
        // The GCC warns if this pointless default-case is missing...
        break;
      }

      // for any other construction type we need an input array
      ASSERT(v != nullptr, "invalid input array");

      switch(constr_type)
      {
      case type_perm:
        // construct from permutation array
        for(Index i(0); i < _num_entries; ++i)
        {
          _perm_pos[i] = v[i];
        }
        calc_swap_from_perm();
        break;

      case type_inv_perm:
        // construct from inverse permutation array
        for(Index i(0); i < _num_entries; ++i)
        {
          _perm_pos[v[i]] = i;
        }
        calc_swap_from_perm();
        break;

      case type_swap:
        // construct from swap array
        for(Index i(0); i < _num_entries; ++i)
        {
          _swap_pos[i] = v[i];
        }
        calc_perm_from_swap();
        break;

      case type_inv_swap:
        // construct from inverse swap array;
        // initialise identity permutation
        for(Index i(0); i < _num_entries; ++i)
        {
          _perm_pos[i] = i;
        }
        // apply swapping in reverse manner
        for(Index i(_num_entries-1); i > 0; --i)
        {
          Index t(_perm_pos[i-1]);
          _perm_pos[i-1] = _perm_pos[v[i-1]];
          _perm_pos[v[i-1]] = t;
        }
        calc_swap_from_perm();
        break;

      default:
        // The GCC warns if this pointless default-case is missing...
        break;
      }
    }

    /**
     * \brief Random constructor
     *
     * This constructor creates a random permutation using a random number generator.
     *
     * \param[in] num_entries
     * The size of the permutation to be created.
     *
     * \param[in] random
     * The random number generator to be used.
     */
    explicit Permutation(
      Index num_entries,
      Random& random)
       :
      _num_entries(num_entries),
      _perm_pos(new Index[num_entries]),
      _swap_pos(new Index[num_entries])
    {
      ASSERT_(num_entries > 0);
      for(Index i(0); i+1 < _num_entries; ++i)
      {
        _swap_pos[i] = random(i, _num_entries-1);
      }
      _swap_pos[_num_entries-1] = _num_entries-1;
      calc_perm_from_swap();
    }

    /**
     * \brief Copy/Inversion constructor
     *
     * This constructor either copies or inverts another permutation.
     *
     * \param[in] other
     * The permutation that is to be copied or inverted.
     *
     * \param[in] invert
     * If \c false, the other permutation will be copied, otherwise it will be inverted.
     */
    explicit Permutation(
      const Permutation& other,
      bool invert = false)
       :
      _num_entries(other.size()),
      _perm_pos(new Index[_num_entries]),
      _swap_pos(new Index[_num_entries])
    {
      if(!invert)
      {
        // copy both arrays
        for(Index i(0); i < _num_entries; ++i)
        {
          _perm_pos[i] = other._perm_pos[i];
          _swap_pos[i] = other._swap_pos[i];
        }
      }
      else
      {
        // invert permutation array
        for(Index i(0); i < _num_entries; ++i)
        {
          _perm_pos[other._perm_pos[i]] = i;
        }
        // and compute swap array
        calc_swap_from_perm();
      }
    }

    /// virtual destructor
    virtual ~Permutation()
    {
      if(_swap_pos != nullptr)
      {
        delete [] _swap_pos;
      }
      if(_perm_pos != nullptr)
      {
        delete [] _perm_pos;
      }
    }

    /// returns the size of the permutation
    Index size() const
    {
      return _num_entries;
    }

    /// returns the permute-position array
    Index* get_perm_pos()
    {
      return _perm_pos;
    }

    /** \copydoc get_perm_pos() */
    const Index* get_perm_pos() const
    {
      return _perm_pos;
    }

    /// returns the swap-position array
    Index* get_swap_pos()
    {
      return _swap_pos;
    }

    /** \copydoc get_swap_pos */
    const Index* get_swap_pos() const
    {
      return _swap_pos;
    }

    /// Computes the swap-position array from the permuation-position array.
    void calc_swap_from_perm()
    {
      for(Index i(0); i < _num_entries; ++i)
      {
        // fetch the permutation position and trace it through the swap array
        Index j(_perm_pos[i]);
        while(j < i)
        {
          j = _swap_pos[j];
        }
        _swap_pos[i] = j;
      }
    }

    /// Computes the permuation-position array from the swap-position array.
    void calc_perm_from_swap()
    {
      // initialise identity permuation
      for(Index i(0); i < _num_entries; ++i)
      {
        _perm_pos[i] = i;
      }

      // apply swapping to permutation
      operator()(_perm_pos);
    }

    /**
     * \brief In-Situ permutation operator
     *
     * This operator applies the permutation in-situ.
     *
     * \param[in] x
     * The array that is to be permuted.
     *
     * \param[in] invert
     * Specifies whether to apply the forward (\p false) or inverse (\p true) permutation.
     */
    template<typename Tx_>
    void operator()(Tx_* x, bool invert = false) const
    {
      ASSERT_(x != nullptr);

      if(!invert)
      {
        // apply in-situ swapping
        for(Index i(0); i+1 < _num_entries; ++i)
        {
          Index j(_swap_pos[i]);
          ASSERT((j >= i) && (j < _num_entries), "invalid swap position");
          if(j > i)
          {
            Tx_ t(x[i]);
            x[i] = x[j];
            x[j] = t;
          }
        }
      }
      else
      {
        // apply inverse in-situ swapping
        for(Index i(_num_entries-1); i > 0; --i)
        {
          Index j(_swap_pos[i-1]);
          ASSERT((j >= i-1) && (j < _num_entries), "invalid swap position");
          if(j > i-1)
          {
            Tx_ t(x[i-1]);
            x[i-1] = x[j];
            x[j] = t;
          }
        }
      }
    }

    /**
     * \brief Permutation operator
     *
     * This operator applies the permutation on an array.
     *
     * \param[out] y
     * The array that shall receive the permuted array.
     *
     * \param[in] x
     * The array that is to be permuted.
     *
     * \param[in] invert
     * Specifies whether to apply the forward (\p false) or inverse (\p true) permutation.
     */
    template<typename Ty_, typename Tx_>
    void operator()(Ty_* y, const Tx_* x, bool invert = false) const
    {
      ASSERT_(y != nullptr);
      ASSERT_(x != nullptr);
      if(!invert)
      {
        for(Index i(0); i < _num_entries; ++i)
        {
          y[i] = Ty_(x[_perm_pos[i]]);
        }
      }
      else
      {
        for(Index i(0); i < _num_entries; ++i)
        {
          y[_perm_pos[i]] = Ty_(x[i]);
        }
      }
    }
  }; // class Permutation
} // namespace FEAST

#endif // KERNEL_UTIL_PERMUTATION_HPP
