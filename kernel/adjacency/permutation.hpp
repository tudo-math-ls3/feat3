// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/util/assertion.hpp>
#include <kernel/util/random.hpp>

//includes, system
#include<vector>

namespace FEAT
{
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
      using IndexVector = std::vector < Index>;

      /// Construction type enumeration
      enum class ConstrType
      {
        /// create uninitialized permutation
        none,
        /// create identity permutation
        identity,
        /// create from permutation array
        perm,
        /// create from swap array
        swap,
        /// create from inverse permutation array
        inv_perm,
        /// create from inverse swap array
        inv_swap
      };

    private:
      /// the permute-position vector
      IndexVector _perm_pos;
      /// the swap-position vector
      IndexVector _swap_pos;

    public:
      /// default CTOR
      Permutation() :
        _perm_pos(),
        _swap_pos()
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
       * - \c ConstrType::none \n Create an uninitialized permutation.\n
       *   The permutation array has to be set up after the object is created.\n
       *   The input array \p v is ignored.
       * - \c ConstrType::identity \n Create an identity permutation.\n The input array \p v is ignored.
       * - \c ConstrType::perm \n Interpret the input vector as a permute-position vector.
       * - \c ConstrType::swap \n Interpret the input vector as a swap-position vector.
       * - \c ConstrType::inv_perm \n Interpret the input array as an inverse permute-position array.
       * - \c ConstrType::inv_swap \n Interpret the input array as an inverse swap-position array.
       *
       * \param[in] v
       * The \transient input array for the initialization. The interpretation of the array's
       * content depends on the \p constr_type parameter.
       */
      explicit Permutation(
        Index num_entries,
        ConstrType constr_type = ConstrType::none,
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
       * The \transient random number generator to be used.
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
        if(!_perm_pos.empty())
          return Permutation(this->size(), ConstrType::perm, _perm_pos.data());
        else
          return Permutation();
      }

      /**
       * \brief Computes the inverse permutation.
       *
       * \returns The inverse permutation.
       */
      Permutation inverse() const
      {
        if(!_perm_pos.empty())
          return Permutation(this->size(), ConstrType::inv_perm, _perm_pos.data());
        else
          return Permutation();
      }

      /// returns the size of the permutation
      Index size() const
      {
        XASSERTM(_perm_pos.size()==_swap_pos.size(), "perm_pos and swap_pos have different sizes");
        return Index(_perm_pos.size());
      }

      /// Checks whether the permutation is empty.
      bool empty() const
      {
        return _perm_pos.empty();
      }

      /// returns the permute-position array
      Index* get_perm_pos()
      {
        return _perm_pos.data();
      }

      /** \copydoc get_perm_pos() */
      const Index* get_perm_pos() const
      {
        return _perm_pos.data();
      }

      /// returns the swap-position array
      Index* get_swap_pos()
      {
        return _swap_pos.data();
      }

      /** \copydoc get_swap_pos */
      const Index* get_swap_pos() const
      {
        return _swap_pos.data();
      }

      template<typename IT_>
      IT_ map(IT_ i) const
      {
        return IT_(_perm_pos.at(std::size_t(i)));
      }

      /// Computes the swap-position vector from the permutation-position vector.
      void calc_swap_from_perm();

      /// Computes the permutation-position vector from the swap-position vector.
      void calc_perm_from_swap();

      /**
       * \brief Applies In-Situ permutation
       *
       * This member function applies the permutation in-situ.
       *
       * \param[in,out] x
       * The \transient array that is to be permuted.
       *
       * \param[in] invert
       * Specifies whether to apply the forward (\p false) or inverse (\p true) permutation.
       *
       * \warning The function does not check, whether the referenced array x is large enough!
       */
      template<typename Tx_>
      void apply(Tx_* x, bool invert = false) const
      {
        XASSERT(x != nullptr);

        if(!invert)
        {
          // apply in-situ swapping
          for(Index i(0); i+1 < this->size(); ++i)
          {
            Index j(_swap_pos[i]);
            ASSERTM((j >= i) && (j < this->size()), "invalid swap position");
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
          for(Index i(this->size() -1); i > 0; --i)
          {
            Index j(_swap_pos[i-1]);
            ASSERTM((j >= i-1) && (j < this->size()), "invalid swap position");
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
       * \brief Applies permutation
       *
       * This function applies the permutation on an array.
       *
       * \param[out] y
       * The \transient array that shall receive the permuted array.
       *
       * \param[in] x
       * The \transient array that is to be permuted.
       *
       * \param[in] invert
       * Specifies whether to apply the forward (\p false) or inverse (\p true) permutation.
       *
       * \warning The function does not check, whether the referenced arrays x and y are large enough!
       */
      template<typename Ty_, typename Tx_>
     // void operator()(Ty_* y, const Tx_* x, bool invert = false) const
      void apply(Ty_* y, const Tx_* x, bool invert = false) const
      {
        XASSERT(y != nullptr);
        XASSERT(x != nullptr);
        if(!invert)
        {
          for(Index i(0); i < this->size(); ++i)
          {
            y[i] = Ty_(x[_perm_pos[i]]);
          }
        }
        else
        {
          for(Index i(0); i < this->size(); ++i)
          {
            y[_perm_pos[i]] = Ty_(x[i]);
          }
        }
      }

      /**
      * \brief concatenates two permutations
      *
      * Let P1 be the permutation denoted by this permutation and let P2 denote the permutation p,
      * then this function computes the contatenation P3 := P1 * P2, such that P3(x) = P1(P2(x)).
      *
      * \param[in, out] p
      * The concatenated permutation of the current permutation and p
      */
      void concat(const Permutation& p);

    }; // class Permutation
  } // namespace Adjacency
} // namespace FEAT
