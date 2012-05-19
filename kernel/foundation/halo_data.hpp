#pragma once
#ifndef KERNEL_FOUNDATION_HALO_DATA_HH
#define KERNEL_FOUNDATION_HALO_DATA_HH 1

#include <kernel/foundation/halo.hpp>
#include <iostream>
#include <cmath>

namespace FEAST
{
  namespace Foundation
  {
    /**
     * \brief HaloData class template represents a static halo
     *
     * HaloData is used for storing static halos.
     *
     * \tparam HaloType_
     * type of halo
     *
     * \tparam VectorType_
     * type for storing the dense data
     *
     * \tparam IndexType_
     * type of the indices
     *
     * \author Markus Geveler
     */
    template<
      typename HaloType_,
      template<typename> class VectorType_,
      typename IndexType_ = Index>
    class HaloData
    {
      public:
        ///CTOR
        HaloData(HaloType_ & halo) :
          _halo(halo),
          _halo_elements(VectorType_<IndexType_>(halo.size())),
          _halo_element_counterparts(VectorType_<IndexType_>(halo.size()))
        {
          for(IndexType_ i(0) ; i < halo.size() ; ++i)
          {
            _halo_element_counterparts[i] = halo.get_element_counterpart(i);
            _halo_elements[i] = halo.get_element(i);
          }
        }

        ///DTOR
        ~HaloData()
        {
        }

        ///public access functions
        HaloType_ & get_halo()
        {
          return _halo;
        }

        IndexType_ get_element_counterpart(IndexType_ index)
        {
          return _halo_element_counterparts[index];
        }

        IndexType_ get_element(IndexType_ index)
        {
          return _halo_elements[index];
        }

        IndexType_ size()
        {
          return _halo_elements.size();
        }

        typename HaloType_::mesh_type_ & get_mesh()
        {
          return _halo.get_mesh();
        }

        IndexType_ & get_other()
        {
          return _halo.get_other();
        }

        unsigned get_overlap()
        {
          return _halo.get_overlap();
        }


      private:
        HaloType_ & _halo;

        VectorType_<IndexType_> _halo_elements;
        VectorType_<IndexType_> _halo_element_counterparts;
    };
  }
}
#endif
