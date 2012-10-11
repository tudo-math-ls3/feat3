#pragma once
#ifndef KERNEL_FOUNDATION_HALO_DATA_HH
#define KERNEL_FOUNDATION_HALO_DATA_HH 1

#include <kernel/foundation/mesh.hpp>
#include <kernel/foundation/communication.hpp>
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
     * \tparam Arch_
     * arch parameter of VectorType_
     *
     * \tparam IndexType_
     * type of the indices in VectorType_
     *
     * \author Markus Geveler
     */
    template<
      typename HaloType_,
      template<typename, typename> class VectorType_,
      typename Arch_,
      typename IndexType_ = Index>
    class HaloData
    {
      public:
        ///CTOR
        HaloData(HaloType_ & halo) :
          _halo(halo),
          _halo_elements(halo.size()),
          _halo_element_counterparts(halo.size()),
          _overlap(halo.get_overlap()),
          _level(halo.get_level())
        {
          for(IndexType_ i(0) ; i < halo.size() ; ++i)
          {
            _halo_element_counterparts(i, halo.get_element_counterpart(i));
            _halo_elements(i, halo.get_element(i));
          }
        }

        ///Copy-CTOR
        HaloData(const HaloData & other) :
          _halo(other._halo),
          _halo_elements(other._halo_elements),
          _halo_element_counterparts(other._halo_element_counterparts),
          _overlap(other._overlap),
          _level(other._level)
        {
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

        IndexType_ get_element_counterpart(IndexType_ index) const
        {
          return _halo_element_counterparts(index);
        }

        IndexType_ get_element(IndexType_ index) const
        {
          return _halo_elements(index);
        }

        IndexType_ size() const
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

        IndexType_ & get_other() const
        {
          return _halo.get_other();
        }

        unsigned* get_overlap()
        {
          return (unsigned *) &_overlap;
        }

        unsigned* get_level()
        {
          return (unsigned *) &_level;
        }

        HaloData& operator=(const HaloData& rhs)
        {
          if(this == &rhs)
            return *this;

          this->_halo = rhs._halo;
          this->_halo_elements = rhs._halo_elements;
          this->_halo_element_counterparts = rhs._halo_element_counterparts;
          this->_overlap = rhs._overlap;
          this->_level = rhs._level;
          return *this;
        }

        VectorType_<Arch_, IndexType_>& get_elements()
        {
          return _halo_elements;
        }

        const VectorType_<Arch_, IndexType_>& get_elements() const
        {
          return _halo_elements;
        }

        VectorType_<Arch_, IndexType_>& get_element_counterparts()
        {
          return _halo_element_counterparts;
        }

        const VectorType_<Arch_, IndexType_>& get_element_counterparts() const
        {
          return _halo_element_counterparts;
        }

      private:
        HaloType_ & _halo;

        VectorType_<Arch_, IndexType_> _halo_elements;
        VectorType_<Arch_, IndexType_> _halo_element_counterparts;

        unsigned _overlap;

        PolytopeLevels _level;
    };
  }
}
#endif
