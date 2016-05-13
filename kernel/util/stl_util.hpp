#pragma once
#ifndef KERNEL_FOUNDATION_STL_UTIL_HPP
#define KERNEL_FOUNDATION_STL_UTIL_HPP 1

#include<kernel/base_header.hpp>
#include<vector>

namespace FEAST
{
  struct STLUtil
  {
    template<template<typename, typename> class ST_, typename DT_>
      static ST_<DT_, std::allocator<DT_> > set_union(ST_<DT_, std::allocator<DT_> >& a, ST_<DT_, std::allocator<DT_> >& b)
      {
        ST_<DT_, std::allocator<DT_> > r(a.size() + b.size());

        std::sort(a.begin(), a.end());
        std::sort(b.begin(), b.end());

        auto iter(std::set_union(a.begin(), a.end(), b.begin(), b.end(), r.begin()));
        r.resize(Index(iter - r.begin()));

        return r;
      }

    template<template<typename, typename> class ST_, typename DT_>
      static ST_<DT_, std::allocator<DT_> > set_intersection(ST_<DT_, std::allocator<DT_> >& a, ST_<DT_, std::allocator<DT_> >& b)
      {
        ST_<DT_, std::allocator<DT_> > r(std::max(a.size(), b.size()));

        std::sort(a.begin(), a.end());
        std::sort(b.begin(), b.end());

        auto iter(std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), r.begin()));
        r.resize(Index(iter - r.begin()));

        return r;
      }

    template<template<typename, typename> class ST_, typename DT_>
      static ST_<DT_, std::allocator<DT_> > set_difference(ST_<DT_, std::allocator<DT_> >& a, ST_<DT_, std::allocator<DT_> >& b)
      {
        ST_<DT_, std::allocator<DT_> > r(a.size());

        std::sort(a.begin(), a.end());
        std::sort(b.begin(), b.end());

        auto iter(std::set_difference(a.begin(), a.end(), b.begin(), b.end(), r.begin()));
        r.resize(Index(iter - r.begin()));

        return r;
      }

    template<template<typename, typename> class ST_, typename DT_>
      static ST_<DT_, std::allocator<DT_> > sorted_seq_union(const ST_<DT_, std::allocator<DT_> >& a, const ST_<DT_, std::allocator<DT_> >& b)
      {
        ST_<DT_, std::allocator<DT_> > r(a.size() + b.size());

        auto iter(std::set_union(a.begin(), a.end(), b.begin(), b.end(), r.begin()));
        r.resize(Index(iter - r.begin()));

        return r;
      }

    template<template<typename, typename> class ST_, typename DT_>
      static ST_<DT_, std::allocator<DT_> > sorted_seq_intersection(const ST_<DT_, std::allocator<DT_> >& a, const ST_<DT_, std::allocator<DT_> >& b)
      {
        ST_<DT_, std::allocator<DT_> > r(std::max(a.size(), b.size()));

        auto iter(std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), r.begin()));
        r.resize(Index(iter - r.begin()));

        return r;
      }

    template<template<typename, typename> class ST_, typename DT_>
      static ST_<DT_, std::allocator<DT_> > sorted_seq_difference(const ST_<DT_, std::allocator<DT_> >& a, const ST_<DT_, std::allocator<DT_> >& b)
      {
        ST_<DT_, std::allocator<DT_> > r(a.size());

        auto iter(std::set_difference(a.begin(), a.end(), b.begin(), b.end(), r.begin()));
        r.resize(Index(iter - r.begin()));

        return r;
      }
  };
}
#endif
