// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/adjacency/permutation.hpp>
#include <kernel/util/random.hpp>

namespace FEAT
{
  namespace Adjacency
  {
    Permutation::Permutation(Index num_entries, ConstrType constr_type, const Index* v) :
      _perm_pos(num_entries),
      _swap_pos(num_entries)
    {
      XASSERTM(num_entries > 0, "cannot create empty permutation");
      switch(constr_type)
      {
      case type_none:
        // leave vectors uninitialized
        return;

      case type_identity:
        // initialize identity permutation
        for(Index i(0); i < num_entries; ++i)
        {
          _perm_pos[i] = _swap_pos[i] = i;
        }
        return;

      default:
        // The GCC warns if this pointless default-case is missing...
        break;
      }

      // for any other construction type we need an input array
      XASSERTM(v != nullptr, "invalid input array");

      switch(constr_type)
      {
      case type_perm:
        // construct from permutation array
        for(Index i(0); i < num_entries; ++i)
        {
          _perm_pos[i] = v[i];
        }
        calc_swap_from_perm();
        break;

      case type_inv_perm:
        // construct from inverse permutation array
        for(Index i(0); i < num_entries; ++i)
        {
          _perm_pos[v[i]] = i;
        }
        calc_swap_from_perm();
        break;

      case type_swap:
        // construct from swap array
        for(Index i(0); i < num_entries; ++i)
        {
          _swap_pos[i] = v[i];
        }
        calc_perm_from_swap();
        break;

      case type_inv_swap:
        // construct from inverse swap array;
        // initialize identity permutation
        for(Index i(0); i < num_entries; ++i)
        {
          _perm_pos[i] = i;
        }
        // apply swapping in reverse manner
        for(Index i(num_entries-1); i > 0; --i)
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

    Permutation::Permutation(Index num_entries, Random& random) :
      _perm_pos(num_entries),
      _swap_pos(num_entries)
    {
      XASSERTM(num_entries > 0, "cannot create empty random permutation");
      for(Index i(0); i+1 < num_entries; ++i)
      {
        _swap_pos[i] = random(i, num_entries-1);
      }
      _swap_pos[num_entries-1] = num_entries-1;
      calc_perm_from_swap();
    }

    Permutation::Permutation(Permutation&& other) :
      _perm_pos(std::forward<IndexVector>(other._perm_pos)),
      _swap_pos(std::forward<IndexVector>(other._swap_pos))
    {}

    Permutation& Permutation::operator=(Permutation&& other)
    {
      // avoid self-move
      if(this == &other)
        return *this;

      _perm_pos = std::forward<IndexVector>(other._perm_pos);
      _swap_pos = std::forward<IndexVector>(other._swap_pos);

      return *this;
    }

    Permutation::~Permutation()
    {}

    void Permutation::calc_swap_from_perm()
    {
      for(Index i(0); i < this->size(); ++i)
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

    void Permutation::calc_perm_from_swap()
    {
      // initialize identity permuation
      for(Index i(0); i < this->size(); ++i)
      {
        _perm_pos[i] = i;
      }

      // apply swapping to permutation
      this->apply(_perm_pos.data());
    }
  } // namespace Adjacency
} // namespace FEAT
