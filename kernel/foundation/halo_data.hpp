#pragma once
#ifndef KERNEL_FOUNDATION_HALO_DATA_HH
#define KERNEL_FOUNDATION_HALO_DATA_HH 1

#include <kernel/foundation/halo.hpp>
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
    class HaloData : public Communicateable<HaloData<HaloType_, VectorType_, Arch_, IndexType_> >
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

        ///Implementation of Communicateable interface
        void send_recv(int otherrank,
                       HaloData<HaloType_, VectorType_, Arch_, IndexType_>& otherdata,
                       int myrank)
        {
          Comm<Archs::Parallel>::send_recv(_halo_elements.elements(), _halo_elements.size(), otherrank, otherdata.get_halo_elements().elements(), myrank);
          Comm<Archs::Parallel>::send_recv(_halo_element_counterparts.elements(), _halo_element_counterparts.size(), otherrank, otherdata.get_halo_element_counterparts().elements(), myrank);
          //Comm<Archs::Parallel>::send_recv((unsigned *)&_overlap, 1, otherrank, &otherdata.get_overlap(), myrank);
          //Comm<Archs::Parallel>::send_recv((unsigned *)&_level, 1, otherrank, &otherdata.get_level(), myrank);
        }

        VectorType_<Arch_, IndexType_>& get_halo_elements()
        {
          return _halo_elements;
        }

        const VectorType_<Arch_, IndexType_>& get_halo_elements() const
        {
          return _halo_elements;
        }

        const VectorType_<Arch_, IndexType_>& get_halo_element_counterparts() const
        {
          return _halo_element_counterparts;
        }

        VectorType_<Arch_, IndexType_>& get_halo_element_counterparts()
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
