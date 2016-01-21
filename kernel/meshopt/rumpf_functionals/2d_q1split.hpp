#pragma once
#ifndef KERNEL_MESHOPT_RUMPF_FUNCTIONALS_Q1SPLIT_HPP
#define KERNEL_MESHOPT_RUMPF_FUNCTIONALS_Q1SPLIT_HPP 1

#include <kernel/base_header.hpp>

namespace FEAST
{
  namespace Meshopt
  {
    /// \cond internal
    /**
     * \brief Class template for Rumpf functionals in 2d using the Q1 hack
     *
     * This variant is for Q1 transformations only and uses an idea by Steffen Basting: Split every hypercube
     * into simplices and evaluate the Rumpf functional for the P1 transformation of those. Do this for every
     * possible way of splitting the hypercube into simplices.
     *
     * \verbatim
        Hypercube:    Splitting 1:   Splitting 2
         2-----3       2-----3        2-----3
         |     |       | \   |        |   / |
         |     |       |   \ |        | /   |
         0-----1       0-----1        0-----1
     * \endverbatim
     * Splitting 1 gives permutations 1 and 2, splitting 2 gives permutations 3 and 4.
     *
     * This means that the local simplex reference cells have to be scaled differently, according to this:
     *
     * For given lambda = volume(cell), how do we need to chose h such that two Rumpf reference simplices of
     * scale h have the same volume?
     * They have volume = 2 * sqrt(3)/4 * h, so h = 2/sqrt(3) * volume =: hcoeff[0] * volume
     *
     * For given lambda = volume(cell), how do we need to chose h such that the det + 1/det part of the
     * functional has a local minimum for two Rumpf reference simplices of scale h?
     * Lengthy computation with maple says that hcoeff[1] = sqrt(hcoeff[0]) is the right constant.
     *
     * \author Jordi Paul
     */
    template<typename DataType_, template<typename, typename> class BaseClass_>
    class RumpfFunctionalQ1Split<DataType_, Shape::Hypercube<2>, BaseClass_> :
    public BaseClass_<DataType_, Shape::Simplex<2>>
    {
      public:
        /// Our data type
        typedef DataType_ DataType;
        /// Our shape type
        typedef Shape::Hypercube<2> ShapeType;
        /// Our base class
        typedef BaseClass_<DataType_, Shape::Simplex<2>> BaseClass;

        /// In 2d, each hypercube is split into 2 simplices and there are two possible permutations
        static constexpr int n_perms = 4;
        /// For each of the n_perms permutations, this has the indices of the ShapeType::dimension+1 vertices forming
        /// the corresponding simplex
        const int perm[4][3];
        /// Factors for scaling the simplex reference cells
        const DataType hcoeffs[ShapeType::dimension];

      private:
        /// Rescaled local optimal scale
        Tiny::Vector<DataType, ShapeType::dimension> _hsplit;
        /// Buffer in which to collect the vertex coordinates of the simplices
        Tiny::Matrix<DataType, ShapeType::dimension+1, ShapeType::dimension> _xsplit;

      public:
        /**
         * \brief Constructor
         */
        explicit RumpfFunctionalQ1Split(
          const DataType fac_norm_,
          const DataType fac_det_,
          const DataType fac_cof_,
          const DataType fac_reg_) : BaseClass(fac_norm_, fac_det_, fac_cof_, fac_reg_),
          perm{ {0, 1, 2}, {1, 3, 2}, {0, 3, 2}, {0, 1, 3} },
          hcoeffs{ DataType(2)/Math::sqrt(DataType(3)), Math::sqrt(DataType(2)/Math::sqrt(DataType(3)))},
          _hsplit(DataType(0)),
          _xsplit(DataType(0))
          {
          }

        /// Explicitly delete the default constructor
        RumpfFunctionalQ1Split() = delete;

        /**
         * \brief Computes value the Rumpf functional on one element.
         */
        template<typename Tx_, typename Th_>
        DataType compute_local_functional(const Tx_& x, const Th_& h)
        {

          DataType norm_A(0);
          DataType det_A(0);
          DataType rec_det_A(0);

          // Compute rescaled h for the split
          _hsplit = h;
          for(int i(0); i < ShapeType::dimension; ++i)
            _hsplit(i)*=hcoeffs[i];

          for(int p(0); p < n_perms; ++p)
          {
            for(int j(0); j < ShapeType::dimension+1; ++j)
              _xsplit[j] = x[perm[p][j]];

            norm_A += BaseClass::compute_norm_A(_xsplit, _hsplit);
            det_A += BaseClass::compute_det_A(_xsplit, _hsplit);
            rec_det_A += BaseClass::compute_rec_det_A(_xsplit, _hsplit);
          }

          return (this->_fac_norm*norm_A
              + this->_fac_det*det_A
              + this->_fac_rec_det*rec_det_A)/DataType(n_perms);

        }

        /**
         * \brief Computes value the Rumpf functional on one element.
         */
        template<typename Tx_, typename Th_>
        DataType compute_local_functional(const Tx_& x, const Th_& h,
        DataType& func_norm,
        DataType& func_det,
        DataType& func_rec_det)
        {
          DataType norm_A(0);
          DataType det_A(0);
          DataType rec_det_A(0);

          // Compute rescaled h for the split
          _hsplit = h;
          for(int i(0); i < ShapeType::dimension; ++i)
            _hsplit(i)*=hcoeffs[i];

          for(int p(0); p < n_perms; ++p)
          {
            for(int j(0); j < ShapeType::dimension+1; ++j)
              _xsplit[j] = x[perm[p][j]];

            norm_A += BaseClass::compute_norm_A(_xsplit, _hsplit);
            det_A += BaseClass::compute_det_A(_xsplit, _hsplit);
            rec_det_A += BaseClass::compute_rec_det_A(_xsplit, _hsplit);
          }

          func_norm = this->_fac_norm*norm_A/DataType(n_perms);
          func_det = this->_fac_det*det_A/DataType(n_perms);
          func_rec_det = this->_fac_rec_det*rec_det_A/DataType(n_perms);

          return (this->_fac_norm*norm_A
              + this->_fac_det*det_A
              + this->_fac_rec_det*rec_det_A)/DataType(n_perms);

        }
        /**
         * \brief Computes the functional gradient for one cell
         **/
        template<typename Tx_, typename Th_, typename Tgrad_>
        void NOINLINE compute_local_grad( const Tx_& x, const Th_& h, Tgrad_& grad)
        {
          grad.format(DataType(0));
          Tiny::Matrix<DataType, ShapeType::dimension+1, Tx_::n> grad_split(DataType(0));

          // Compute rescaled h for the split
          _hsplit = h;
          for(int i(0); i < ShapeType::dimension; ++i)
            _hsplit(i)*=hcoeffs[i];

          for(int p(0); p < n_perms; ++p)
          {
            for(int j(0); j < ShapeType::dimension+1; ++j)
              _xsplit[j] = x[perm[p][j]];

            BaseClass::compute_local_grad( _xsplit, _hsplit, grad_split);

            for(int j(0); j < ShapeType::dimension+1; ++j)
              grad[perm[p][j]] += grad_split[j];
          }

          grad*=(DataType(1)/DataType(n_perms));

        }
        /**
         * \brief The class name
         *
         * \returns String with the class name
         */
        static String name()
        {
          return "RumpfFunctionalQ1Split<"+ShapeType::name()+", "+BaseClass::name()+">";
        }

        /**
         * \brief Prints object parameters
         */
        void print()
        {
          std::cout << name() << std::endl;
          BaseClass::print();
        }

    }; // class RumpfFunctionalQ1Split

    /**
     * \brief This is a dummy implementation of the Q1Split functional for simplices
     *
     * There is nothing to split here, so this coincides with the functional for the Simplex<2> shape.
     * It still has to be here because of runtime switches for creating functionals need to be able to
     * chose this at compile time, even though it should never be instantiated.
     */
    template<typename DataType_, template<typename, typename> class BaseClass_>
    class RumpfFunctionalQ1Split<DataType_, Shape::Simplex<2>, BaseClass_> :
    public BaseClass_<DataType_, Shape::Simplex<2>>
    {
      public:
        /// Our data type
        typedef DataType_ DataType;
        /// Our shape type
        typedef Shape::Simplex<2> ShapeType;
        /// Our base class
        typedef BaseClass_<DataType_, Shape::Simplex<2>> BaseClass;

        using BaseClass::BaseClass;
        RumpfFunctionalQ1Split() = delete;
    };

    /// \endcond
  } // namespace Meshopt
} // namespace FEAST

#endif // KERNEL_MESHOPT_RUMPF_FUNCTIONALS_Q1SPLIT_HPP
