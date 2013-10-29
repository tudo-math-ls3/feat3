#pragma once
#ifndef KERNEL_ADJACENCY_PERMUTATION_HPP
#define KERNEL_ADJACENCY_PERMUTATION_HPP 1

// includes, FEAST
#include <kernel/util/assertion.hpp>

namespace FEAST
{
  class Random;

  namespace Adjacency
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
      /// default CTOR
      Permutation() :
        _num_entries(0),
        _perm_pos(nullptr),
        _swap_pos(nullptr)
      {
      }

      /**
       * \brief Constructor
       *
       * \param[in] num_entries
       * The size of the permutation to be created.
       *
       * \param[in] constr_type
       * Specifies the construction type:
       * - \c type_none \n Create an uninitialised permutation.\n
       *   The permutation array has to be set up after the object is created.\n
       *   The input array \p v is ignored.
       * - \c type_identity \n Create an identity permutation.\n The input array \p v is ignored.
       * - \c type_perm \n Interpret the input array as a permute-position array.
       * - \c type_swap \n Interpret the input array as a swap-position array.
       * - \c type_inv_perm \n Interpret the input array as an inverse permute-position array.
       * - \c type_inv_swap \n Interpret the input array as an inverse swap-position array.
       *
       * \param[in] v
       * The input array for the initialisation. The interpretation of the array's content depends on the
       * \p constr_type parameter.
       */
      explicit Permutation(
        Index num_entries,
        ConstrType constr_type = type_none,
        const Index* v = nullptr);

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
        Random& random);

      /// move ctor
      Permutation(Permutation&& other);
      /// move-assign operator
      Permutation& operator=(Permutation&&);

      /// virtual destructor
      virtual ~Permutation();

      /**
       * \brief Clones this permutation.
       *
       * \returns A deep-copy of this permutation
       */
      Permutation clone() const
      {
        if(_perm_pos != nullptr)
          return Permutation(_num_entries, type_perm, _perm_pos);
        else
          return Permutation();
      }

      /**
       * \brief Computes the inverse permutation.
       *
       * \returns The inverse permutiation.
       */
      Permutation inverse() const
      {
        if(_perm_pos != nullptr)
          return Permutation(_num_entries, type_inv_perm, _perm_pos);
        else
          return Permutation();
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
      void calc_swap_from_perm();

      /// Computes the permuation-position array from the swap-position array.
      void calc_perm_from_swap();

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
  } // namespace Adjacency
} // namespace FEAST

#endif // KERNEL_ADJACENCY_PERMUTATION_HPP
