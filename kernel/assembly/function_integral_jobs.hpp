// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ASSEMBLY_FUNCTION_INTEGRAL_JOBS_HPP
#define KERNEL_ASSEMBLY_FUNCTION_INTEGRAL_JOBS_HPP 1

#include <kernel/assembly/base.hpp>
#include <kernel/analytic/function.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/util/dist.hpp>

namespace FEAT
{
  namespace Assembly
  {
    /**
     * \brief Function integral info class
     *
     * This class is used to store and manage the results of a AnalyticFunctionIntegralJob,
     * DiscreteFunctionIntegralJob or ErrorFunctionIntegralJob. It contains all the integral
     * values as well as a bunch of auxiliary helper functions which are used during the
     * assembly process.
     *
     * \tparam DataType_
     * The scalar datatype used during the assembly.
     *
     * \tparam ValueType
     * The type of a function value. This corresponds to DataType_ for scalar functions and to
     * Tiny::Vector<DataType_,...> for (blocked) vector fields.
     *
     * \tparam GradientType
     * The type of a function gradient. This corresponds to Tiny::Vector<DataType_,...> for scalar
     * functions and to Tiny::Matrix<DataType_, ...> for (blocked) vector fields.
     *
     * \tparam HessianType
     * The type of a function hessian. This corresponds to Tiny::Matrix<DataType_,...> for scalar
     * functions and to Tiny::Tensor3<DataType_, ...> for (blocked) vector fields.
     *
     * \author Peter Zajac
     */
    template<typename DataType_, typename ValueType_, typename GradientType_, typename HessianType_>
    class FunctionIntegralInfo
    {
    public:
      /// the underlying scalar data type
      typedef DataType_ DataType;
      /// the function value type
      typedef ValueType_ ValueType;
      /// the function gradient type
      typedef GradientType_ GradientType;
      /// the function hessian type
      typedef HessianType_ HessianType;

      /// the maximum derivative for which the integrals have been computed
      int max_der;

      /**
       * \brief Integral of the function value
       *
       * This entry stores the integral of the function value if max_der >= 0, i.e. it is
       * \f[v \in\mathbb{R};\quad v := \int_\Omega u(x)~dx \f]
       */
      ValueType_ value;

      /**
       * \brief Integral of the function gradient
       *
       * This entry stores the integral of the function gradient if max_der >= 1, i.e. it is
       * \f[g \in\mathbb{R}^n;\quad g_i := \int_\Omega \partial_i u(x)~dx \f]
       */
      GradientType_ grad;

      /**
       * \brief Integral of the function hessian
       *
       * This entry stores the integral of the function gradient if max_der >= 2, i.e. it is
       * \f[h \in\mathbb{R}^{n\times n};\quad h_{ij} := \int_\Omega \partial_i\partial_j u(x)~dx \f]
       */
      HessianType_ hess;

      /**
       * \brief squared H0-norm (aka L2-norm) of the function
       *
       * This entry stores the squared H0-norm of the function if max_der >= 0, i.e. it is
       * \f[|u|_{H^0} := \int_\Omega (u(x))^2~dx \f]
       *
       * For vector fields, this is simply the sum of all component-wise H0-norms.
       */
      DataType_ norm_h0_sqr;

      /**
       * \brief squared H1-semi-norm of the function
       *
       * This entry stores the squared H1-semi-norm of the function if max_der >= 1, i.e. it is
       * \f[|u|_{H^1} := \sum_{1\leq i\leq n}\int_\Omega (\partial_i u(x))^2~dx \f]
       *
       * For vector fields, this is simply the sum of all component-wise H1-semi-norms.
       */
      DataType_ norm_h1_sqr;

      /**
       * \brief squared H2-semi-norm of the function
       *
       * This entry stores the squared H2-semi-norm of the function if max_der >= 2, i.e. it is
       * \f[|u|_{H^2} := \sum_{1\leq i,j\leq n}\int_\Omega (\partial_i \partial_j u(x))^2~dx \f]
       *
       * For vector fields, this is simply the sum of all component-wise H2-semi-norms.
       */
      DataType_ norm_h2_sqr;

      /**
       * \brief L1-norm of the function
       *
       * This entry stores the L1-norm of the function if max_der >= 0, i.e. it is
       * \f[|u|_{L^1} := \int_\Omega |u(x)|~dx \f]
       *
       * For vector fields, this is simply the sum of all component-wise L1-norms.
       */
      DataType_ norm_l1;

      /**
       * \brief Lmax-norm of the function
       *
       * This entry stores the Lmax-norm of the function if max_der >= 0, i.e. it is
       * \f[|u|_{L^\infty} := \max_{x\in\Omega} |u(x)| \f]
       *
       * For vector fields, this is simply the sum of all component-wise Lmax-norms.
       */
      DataType_ norm_lmax;

      /// the component-wise squared H0-norm (only for vector fields)
      ValueType_ norm_h0_sqr_comp;
      /// the component-wise squared H1-semi-norm (only for vector fields)
      ValueType_ norm_h1_sqr_comp;
      /// the component-wise squared H2-semi-norm (only for vector fields)
      ValueType_ norm_h2_sqr_comp;
      /// the component-wise L1-norm (only for vector fields)
      ValueType_ norm_l1_comp;
      /// the component-wise Lmax-norm (only for vector fields)
      ValueType_ norm_lmax_comp;

    private:
      /// \cond internal
      static void _push_lmax(DataType_& a, const DataType_& b)
      {
        a = Math::max(a, b);
      }

      template<int n_, int sn_>
      static void _push_lmax(Tiny::Vector<DataType_, n_, sn_>& a, const Tiny::Vector<DataType_, n_, sn_>& b)
      {
        for(int i(0); i < n_; ++i)
          a[i] = Math::max(a[i], b[i]);
      }

      static void _write_to(DataType_* ptr, std::size_t k, const DataType_& x)
      {
        ptr[k] = x;
        ++k;
      }

      template<int n_, int sn_>
      static void _write_to(DataType_* ptr, std::size_t k, const Tiny::Vector<DataType_, n_, sn_>& x)
      {
        for(int i(0); i < n_; ++i, ++k)
          ptr[k] = x[i];
      }

      template<int m_, int n_, int sm_, int sn_>
      static void _write_to(DataType_* ptr, std::size_t k, const Tiny::Matrix<DataType_, m_, n_, sm_, sn_>& x)
      {
        for(int i(0); i < m_; ++i)
          _write_to(ptr, k, x[i]);
      }

      template<int l_, int m_, int n_, int sl_, int sm_, int sn_>
      static void _write_to(DataType_* ptr, std::size_t k, const Tiny::Tensor3<DataType_, l_, m_, n_, sl_, sm_, sn_>& x)
      {
        for(int i(0); i < l_; ++i)
          _write_to(ptr, k, x[i]);
      }

      static void _read_from(DataType_* ptr, std::size_t k, DataType_& x)
      {
        x = ptr[k];
        ++k;
      }

      template<int n_, int sn_>
      static void _read_from(DataType_* ptr, std::size_t k, Tiny::Vector<DataType_, n_, sn_>& x)
      {
        for(int i(0); i < n_; ++i, ++k)
          x[i] = ptr[k];
      }

      template<int m_, int n_, int sm_, int sn_>
      static void _read_from(DataType_* ptr, std::size_t k, Tiny::Matrix<DataType_, m_, n_, sm_, sn_>& x)
      {
        for(int i(0); i < m_; ++i)
          _read_from(ptr, k, x[i]);
      }

      template<int l_, int m_, int n_, int sl_, int sm_, int sn_>
      static void _read_from(DataType_* ptr, std::size_t k, Tiny::Tensor3<DataType_, l_, m_, n_, sl_, sm_, sn_>& x)
      {
        for(int i(0); i < l_; ++i)
          _read_from(ptr, k, x[i]);
      }

      static String _print_val(const DataType_& x, int prec)
      {
        return stringify_fp_sci(x, prec);
      }

      template<int n_, int sn_>
      static String _print_val(const Tiny::Vector<DataType_, n_, sn_>& x, int prec)
      {
        String s = "[";
        for(int i(0); i < n_; ++i)
          (s += " ") += stringify_fp_sci(x[i], prec);
        s += " ]";
        return s;
      }

      template<int n_, int sn_>
      static String _print_val_sqrt(const Tiny::Vector<DataType_, n_, sn_>& x, int prec)
      {
        String s = "[";
        for(int i(0); i < n_; ++i)
          (s += " ") += stringify_fp_sci(Math::sqrt(x[i]), prec);
        s += " ]";
        return s;
      }

      static String _print_norm(const DataType_& x, const DataType_&, int prec)
      {
        return _print_val(x, prec);
      }

      template<int n_, int sn_>
      static String _print_norm(const DataType_& x, const Tiny::Vector<DataType_, n_, sn_>& v, int prec)
      {
        return _print_val(x, prec) + " " + _print_val(v, prec);
      }

      static String _print_norm_sqrt(const DataType_& x, const DataType_&, int prec)
      {
        return _print_val(Math::sqrt(x), prec);
      }

      template<int n_, int sn_>
      static String _print_norm_sqrt(const DataType_& x, const Tiny::Vector<DataType_, n_, sn_>& v, int prec)
      {
        return _print_val(Math::sqrt(x), prec) + " " + _print_val_sqrt(v, prec);
      }
      /// \endcond

    public:
      /**
       * \brief Constructor
       */
      FunctionIntegralInfo() :
        max_der(0),
        value(DataType(0)),
        grad(DataType(0)),
        hess(DataType(0)),
        norm_h0_sqr(DataType(0)),
        norm_h1_sqr(DataType(0)),
        norm_h2_sqr(DataType(0)),
        norm_l1(DataType(0)),
        norm_lmax(DataType(0)),
        norm_h0_sqr_comp(DataType(0)),
        norm_h1_sqr_comp(DataType(0)),
        norm_h2_sqr_comp(DataType(0)),
        norm_l1_comp(DataType(0)),
        norm_lmax_comp(DataType(0))
      {
      }

      /**
       * \brief Assembly helper function: adds another info object
       *
       * This function is used by the assembly and is not meant to be used elsewhere.
       */
      void push(const FunctionIntegralInfo& other)
      {
        value += other.value;
        grad  += other.grad;
        hess  += other.hess;
        norm_h0_sqr += other.norm_h0_sqr;
        norm_h1_sqr += other.norm_h1_sqr;
        norm_h2_sqr += other.norm_h2_sqr;
        norm_l1 += other.norm_l1;
        _push_lmax(norm_lmax, other.norm_lmax);
        norm_h0_sqr_comp += other.norm_h0_sqr_comp;
        norm_h1_sqr_comp += other.norm_h1_sqr_comp;
        norm_h2_sqr_comp += other.norm_h2_sqr_comp;
        norm_l1_comp += other.norm_l1_comp;
        _push_lmax(norm_lmax_comp, other.norm_lmax_comp);
      }

      /**
       * \brief Assembly helper function: adds another function value
       *
       * This function is used by the assembly and is not meant to be used elsewhere.
       *
       * \param[in] omega
       * The cubature weight.
       *
       * \param[in] v
       * The function value in the cubature point.
       */
      void add_value(DataType omega, const DataType& v)
      {
        value += omega * v;
        const DataType av = Math::abs(v);
        const DataType sv = Math::sqr(v);
        norm_h0_sqr += omega * sv;
        norm_l1 += omega * av;
        norm_lmax = Math::max(norm_lmax, av);
      }

      /**
       * \brief Assembly helper function: adds another function value
       *
       * This function is used by the assembly and is not meant to be used elsewhere.
       *
       * \param[in] omega
       * The cubature weight.
       *
       * \param[in] v
       * The function value in the cubature point.
       */
      template<int n_, int sn_>
      void add_value(DataType omega, const Tiny::Vector<DataType, n_, sn_>& v)
      {
        value += omega * v;
        for(int i(0); i < n_; ++i)
        {
          const DataType av = Math::abs(v[i]);
          const DataType sv = Math::sqr(v[i]);
          norm_h0_sqr += omega * sv;
          norm_l1 += omega * av;
          norm_lmax = Math::max(norm_lmax, av);
          norm_h0_sqr_comp[i] += omega * sv;
          norm_l1_comp[i] += omega * av;
          norm_lmax_comp[i] = Math::max(norm_lmax_comp[i], av);
        }
      }

      /**
       * \brief Assembly helper function: adds another function gradient
       *
       * This function is used by the assembly and is not meant to be used elsewhere.
       *
       * \param[in] omega
       * The cubature weight.
       *
       * \param[in] g
       * The function gradient in the cubature point.
       */
      template<int n_, int sn_>
      void add_grad(DataType omega, const Tiny::Vector<DataType, n_, sn_>& g)
      {
        grad += omega * g;
        norm_h1_sqr += omega * g.norm_euclid_sqr();
      }

      /**
       * \brief Assembly helper function: adds another function gradient
       *
       * This function is used by the assembly and is not meant to be used elsewhere.
       *
       * \param[in] omega
       * The cubature weight.
       *
       * \param[in] g
       * The function gradient in the cubature point.
       */
      template<int m_, int n_, int sm_, int sn_>
      void add_grad(DataType omega, const Tiny::Matrix<DataType, m_, n_, sm_, sn_>& g)
      {
        grad += omega * g;
        for(int i(0); i < m_; ++i)
        {
          const DataType sv = g[i].norm_euclid_sqr();
          norm_h1_sqr += omega * sv;
          norm_h1_sqr_comp[i] += omega * sv;
        }
      }

      /**
       * \brief Assembly helper function: adds another function hessian
       *
       * This function is used by the assembly and is not meant to be used elsewhere.
       *
       * \param[in] omega
       * The cubature weight.
       *
       * \param[in] h
       * The function hessian in the cubature point.
       */
      template<int m_, int n_, int sm_, int sn_>
      void add_hess(DataType omega, const Tiny::Matrix<DataType, m_, n_, sm_, sn_>& h)
      {
        hess += omega * h;
        norm_h2_sqr += omega * h.norm_hessian_sqr();
      }

      /**
       * \brief Assembly helper function: adds another function hessian
       *
       * This function is used by the assembly and is not meant to be used elsewhere.
       *
       * \param[in] omega
       * The cubature weight.
       *
       * \param[in] h
       * The function hessian in the cubature point.
       */
      template<int l_, int m_, int n_, int sl_, int sm_, int sn_>
      void add_hess(DataType omega, const Tiny::Tensor3<DataType, l_, m_, n_, sl_, sm_, sn_>& h)
      {
        hess += omega * h;
        for(int i(0); i < l_; ++i)
        {
          const DataType sv = h[i].norm_hessian_sqr();
          norm_h2_sqr += omega * sv;
          norm_h2_sqr_comp[i] += omega * sv;
        }
      }

      /**
       * \brief Synchronizes the function integral information over a communicator
       *
       * This function sums up the function integral information of all patches in a
       * parallel simulation to obtain the information for the global mesh.
       *
       * \param[in] comm
       * The communication over which to synchronize.
       */
      void synchronize(const Dist::Comm& comm)
      {
        static constexpr std::size_t max_len = std::size_t(5) +
          (6u*sizeof(ValueType) + sizeof(GradientType) + sizeof(HessianType)) / sizeof(DataType);
        DataType v[max_len];

        // write all stuff to array
        std::size_t k = 0u;
        _write_to(v, k, value);
        _write_to(v, k, grad);
        _write_to(v, k, hess);
        _write_to(v, k, norm_h0_sqr);
        _write_to(v, k, norm_h1_sqr);
        _write_to(v, k, norm_h2_sqr);
        _write_to(v, k, norm_l1);
        _write_to(v, k, norm_h0_sqr_comp);
        _write_to(v, k, norm_h1_sqr_comp);
        _write_to(v, k, norm_h2_sqr_comp);
        _write_to(v, k, norm_l1_comp);
        XASSERT(k <= max_len);

        // synchronize
        comm.allreduce(v, v, k, Dist::op_sum);

        // read stuff form array
        k = std::size_t(0u);
        _read_from(v, k, value);
        _read_from(v, k, grad);
        _read_from(v, k, hess);
        _read_from(v, k, norm_h0_sqr);
        _read_from(v, k, norm_h1_sqr);
        _read_from(v, k, norm_h2_sqr);
        _read_from(v, k, norm_l1);
        _read_from(v, k, norm_h0_sqr_comp);
        _read_from(v, k, norm_h1_sqr_comp);
        _read_from(v, k, norm_h2_sqr_comp);
        _read_from(v, k, norm_l1_comp);

        // write lmax values
        k = std::size_t(0u);
        _write_to(v, k, norm_lmax);
        _write_to(v, k, norm_lmax_comp);

        // synchronize
        comm.allreduce(v, v, k, Dist::op_max);

        // read lmax values
        k = std::size_t(0u);
        _read_from(v, k, norm_lmax);
        _read_from(v, k, norm_lmax_comp);
      }

      /**
       * \brief Prints all (computed) norms to a formatted string and returns it.
       *
       * \param[in] precision
       * The number of digits to print. Default to 0, which prints a default number of digits.
       *
       * \param[in] pad_size
       * The number of characters for the padding. Should be >= 9 to enable a vertically aligned output.
       *
       * \param[in] pad_char
       * The character used for padding. Defaults to a dot.
       *
       * \returns
       * A formatted multi-line string that contains all the computed norms.
       */
      String print_norms(int precision = 0, std::size_t pad_size = 10u, char pad_char = '.') const
      {
        String s;
        s += String("H0-Norm").pad_back(pad_size, pad_char) + ": "
          + _print_norm_sqrt(norm_h0_sqr, norm_h0_sqr_comp, precision) + "\n";
        if(max_der >= 1)
          s += String("H1-Norm").pad_back(pad_size, pad_char) + ": "
            + _print_norm_sqrt(norm_h1_sqr, norm_h1_sqr_comp, precision) + "\n";
        if(max_der >= 2)
          s += String("H2-Norm").pad_back(pad_size, pad_char) + ": "
            + _print_norm_sqrt(norm_h2_sqr, norm_h2_sqr_comp, precision) + "\n";
        s += String("L1-Norm").pad_back(pad_size, pad_char) + ": "
          + _print_norm(norm_l1, norm_l1_comp, precision) + "\n";
        s += String("Lmax-Norm").pad_back(pad_size, pad_char) + ": "
          + _print_norm(norm_lmax, norm_lmax_comp, precision);
        return s;
      }
    }; // class FunctionIntegralInfo<...>

    /**
     * \brief Helper class to determine the FunctionIntegralInfo type for analytic functions
     *
     * The only purpose of this helper class is to provide the nested #Type typedef, which is an
     * instance of the FunctionIntegralInfo class with the appropriate types for the given function.
     *
     * \tparam DataType_
     * The scalar data type to be used
     *
     * \tparam Function_
     * The analytic function implementation
     */
    template<typename DataType_, typename Function_>
    struct AnalyticFunctionIntegral
    {
      /// the FunctionIntegralInfo type for this analytic function
      typedef FunctionIntegralInfo<DataType_,
        typename Analytic::EvalTraits<DataType_, Function_>::ValueType,
        typename Analytic::EvalTraits<DataType_, Function_>::GradientType,
        typename Analytic::EvalTraits<DataType_, Function_>::HessianType> Type;
    };

    /**
     * \brief Helper class to determine the FunctionIntegralInfo type for discrete functions
     *
     * The only purpose of this helper class is to provide the nested Type typedef, which is an
     * instance of the FunctionIntegralInfo class with the appropriate types for the given function.
     *
     * Note that this generic template does not contain an implementation, because this class is
     * explicitly specialized for all supported vector types.
     *
     * \tparam Vector_
     * The vector type of the coefficient vector. May be one of the following:
     * - LAFEM::DenseVector<...>
     * - LAFEM::DenseVectorBlocked<...>
     *
     * \tparam Space_
     * The finite element spaced used for the discretization.
     */
    template<typename Vector_, typename Space_>
    struct DiscreteFunctionIntegral;

    /// specialization for LAFEM::DenseVector coefficient vectors
    template<typename Mem_, typename DT_, typename IT_, typename Space_>
    struct DiscreteFunctionIntegral<LAFEM::DenseVector<Mem_, DT_, IT_>, Space_>
    {
      /// the FunctionIntegralInfo type for this scalar discrete function
      typedef FunctionIntegralInfo<DT_, DT_,
        Tiny::Vector<DT_, Space_::shape_dim>,
        Tiny::Matrix<DT_, Space_::shape_dim, Space_::shape_dim>> Type;
    };

    /// specialization for LAFEM::DenseVectorBlocked coefficient vectors
    template<typename Mem_, typename DT_, typename IT_, int bs_, typename Space_>
    struct DiscreteFunctionIntegral<LAFEM::DenseVectorBlocked<Mem_, DT_, IT_, bs_>, Space_>
    {
      /// the FunctionIntegralInfo type for this blocked discrete function
      typedef FunctionIntegralInfo<DT_, Tiny::Vector<DT_, bs_>,
        Tiny::Matrix<DT_, bs_, Space_::shape_dim>,
        Tiny::Tensor3<DT_, bs_, Space_::shape_dim, Space_::shape_dim>> Type;
    };

    /// \cond internal
    namespace Intern
    {
      template<int max_der_>
      struct AnaFunIntJobHelper;

      template<>
      struct AnaFunIntJobHelper<0>
      {
        template<typename FID_, typename DT_, typename FEV_, typename IPT_>
        static void work(FID_& fid, const DT_ omega, FEV_& fev, const IPT_& ipt)
        {
          fid.add_value(omega, fev.value(ipt));
        }
      };

      template<>
      struct AnaFunIntJobHelper<1>
      {
        template<typename FID_, typename DT_, typename FEV_, typename IPT_>
        static void work(FID_& fid, const DT_ omega, FEV_& fev, const IPT_& ipt)
        {
          fid.add_value(omega, fev.value(ipt));
          fid.add_grad(omega, fev.gradient(ipt));
        }
      };

      template<>
      struct AnaFunIntJobHelper<2>
      {
        template<typename FID_, typename DT_, typename FEV_, typename IPT_>
        static void work(FID_& fid, const DT_ omega, FEV_& fev, const IPT_& ipt)
        {
          fid.add_value(omega, fev.value(ipt));
          fid.add_grad(omega, fev.gradient(ipt));
          fid.add_hess(omega, fev.hessian(ipt));
        }
      };

      template<int max_der_>
      struct DiscFunIntJobHelper;

      template<>
      struct DiscFunIntJobHelper<0>
      {
        template<typename FID_, typename DT_, typename SED_, typename VT_, int n_, int sn_>
        static void work(FID_& fid, const DT_ omega, const SED_& sed,
          const Tiny::Vector<VT_, n_, sn_>& lvec, int n)
        {
          typename FID_::ValueType    v(DT_(0));
          for(int i(0); i < n; ++i)
          {
            v += lvec[i] * sed.phi[i].value;
          }
          fid.add_value(omega, v);
        }
      };

      template<>
      struct DiscFunIntJobHelper<1>
      {
        // version for scalar functions
        template<typename FID_, typename DT_, typename SED_, int n_, int sn_>
        static void work(FID_& fid, const DT_ omega, const SED_& sed,
          const Tiny::Vector<DT_, n_, sn_>& lvec, int n)
        {
          typename FID_::ValueType    v(DT_(0));
          typename FID_::GradientType g(DT_(0));
          for(int i(0); i < n; ++i)
          {
            v += lvec[i] * sed.phi[i].value;
            g += lvec[i] * sed.phi[i].grad;
          }
          fid.add_value(omega, v);
          fid.add_grad(omega, g);
        }
        // version for vector fields
        template<typename FID_, typename DT_, typename SED_, int d_, int sd_, int n_, int sn_>
        static void work(FID_& fid, const DT_ omega, const SED_& sed,
          const Tiny::Vector<Tiny::Vector<DT_, d_, sd_>, n_, sn_>& lvec, int n)
        {
          typename FID_::ValueType    v(DT_(0));
          typename FID_::GradientType g(DT_(0));
          for(int i(0); i < n; ++i)
          {
            v.axpy(sed.phi[i].value, lvec[i]);
            g.add_outer_product(lvec[i], sed.phi[i].grad);
          }
          fid.add_value(omega, v);
          fid.add_grad(omega, g);
        }
      };

      template<>
      struct DiscFunIntJobHelper<2>
      {
        // version for scalar functions
        template<typename FID_, typename DT_, typename SED_, int n_, int sn_>
        static void work(FID_& fid, const DT_ omega, const SED_& sed,
          const Tiny::Vector<DT_, n_, sn_>& lvec, int n)
        {
          typename FID_::ValueType    v(DT_(0));
          typename FID_::GradientType g(DT_(0));
          typename FID_::HessianType  h(DT_(0));
          for(int i(0); i < n; ++i)
          {
            v += lvec[i] * sed.phi[i].value;
            g += lvec[i] * sed.phi[i].grad;
            h += lvec[i] * sed.phi[i].hess;
          }

          fid.add_value(omega, v);
          fid.add_grad(omega, g);
          fid.add_hess(omega, h);
        }
        // version for vector fields
        template<typename FID_, typename DT_, typename SED_, int d_, int sd_, int n_, int sn_>
        static void work(FID_& fid, const DT_ omega, const SED_& sed,
          const Tiny::Vector<Tiny::Vector<DT_, d_, sd_>, n_, sn_>& lvec, int n)
        {
          typename FID_::ValueType v(DT_(0));
          typename FID_::GradientType g(DT_(0));
          typename FID_::HessianType h(DT_(0));
          for(int i(0); i < n; ++i)
          {
            v.axpy(sed.phi[i].value, lvec[i]);
            g.add_outer_product(lvec[i], sed.phi[i].grad);
            h.add_vec_mat_outer_product(lvec[i], sed.phi[i].hess);
          }

          fid.add_value(omega, v);
          fid.add_grad(omega, g);
          fid.add_hess(omega, h);
        }
      };

      template<int max_der_>
      struct ErrFunIntJobHelper;

      template<>
      struct ErrFunIntJobHelper<0>
      {
        // version for scalar functions
        template<typename FID_, typename DT_, typename FEV_, typename IPT_, typename SED_, int n_, int sn_>
        static void work(FID_& fid, const DT_ omega, FEV_& fev, const IPT_& ipt, const SED_& sed,
          const Tiny::Vector<DT_, n_, sn_>& lvec, int n)
        {
          typename FID_::ValueType    v = fev.value(ipt);
          for(int i(0); i < n; ++i)
          {
            v -= lvec[i] * sed.phi[i].value;
          }

          fid.add_value(omega, v);
        }

        // version for vector fields
        template<typename FID_, typename DT_, typename FEV_, typename IPT_, typename SED_, int d_, int sd_, int n_, int sn_>
        static void work(FID_& fid, const DT_ omega, FEV_& fev, const IPT_& ipt, const SED_& sed,
          const Tiny::Vector<Tiny::Vector<DT_, d_, sd_>, n_, sn_>& lvec, int n)
        {
          typename FID_::ValueType    v = fev.value(ipt);
          for(int i(0); i < n; ++i)
          {
            v.axpy(-sed.phi[i].value, lvec[i]);
          }

          fid.add_value(omega, v);
        }
      };

      template<>
      struct ErrFunIntJobHelper<1>
      {
        // version for scalar functions
        template<typename FID_, typename DT_, typename FEV_, typename IPT_, typename SED_, int n_, int sn_>
        static void work(FID_& fid, const DT_ omega, FEV_& fev, const IPT_& ipt, const SED_& sed,
          const Tiny::Vector<DT_, n_, sn_>& lvec, int n)
        {
          typename FID_::ValueType    v = fev.value(ipt);
          typename FID_::GradientType g = fev.gradient(ipt);
          for(int i(0); i < n; ++i)
          {
            v -= lvec[i] * sed.phi[i].value;
            g -= lvec[i] * sed.phi[i].grad;
          }

          fid.add_value(omega, v);
          fid.add_grad(omega, g);
        }

        // version for vector fields
        template<typename FID_, typename DT_, typename FEV_, typename IPT_, typename SED_, int d_, int sd_, int n_, int sn_>
        static void work(FID_& fid, const DT_ omega, FEV_& fev, const IPT_& ipt, const SED_& sed,
          const Tiny::Vector<Tiny::Vector<DT_, d_, sd_>, n_, sn_>& lvec, int n)
        {
          typename FID_::ValueType    v = fev.value(ipt);
          typename FID_::GradientType g = fev.gradient(ipt);
          for(int i(0); i < n; ++i)
          {
            v.axpy(-sed.phi[i].value, lvec[i]);
            g.add_outer_product(lvec[i], sed.phi[i].grad, -DT_(1));
          }

          fid.add_value(omega, v);
          fid.add_grad(omega, g);
        }
      };

      template<>
      struct ErrFunIntJobHelper<2>
      {
        // version for scalar functions
        template<typename FID_, typename DT_, typename FEV_, typename IPT_, typename SED_, int n_, int sn_>
        static void work(FID_& fid, const DT_ omega, FEV_& fev, const IPT_& ipt, const SED_& sed,
          const Tiny::Vector<DT_, n_, sn_>& lvec, int n)
        {
          typename FID_::ValueType    v = fev.value(ipt);
          typename FID_::GradientType g = fev.gradient(ipt);
          typename FID_::HessianType  h = fev.hessian(ipt);
          for(int i(0); i < n; ++i)
          {
            v -= lvec[i] * sed.phi[i].value;
            g -= lvec[i] * sed.phi[i].grad;
            h -= lvec[i] * sed.phi[i].hess;
          }

          fid.add_value(omega, v);
          fid.add_grad(omega, g);
          fid.add_hess(omega, h);
        }

        // version for vector fields
        template<typename FID_, typename DT_, typename FEV_, typename IPT_, typename SED_, int d_, int sd_, int n_, int sn_>
        static void work(FID_& fid, const DT_ omega, FEV_& fev, const IPT_& ipt, const SED_& sed,
          const Tiny::Vector<Tiny::Vector<DT_, d_, sd_>, n_, sn_>& lvec, int n)
        {
          typename FID_::ValueType    v = fev.value(ipt);
          typename FID_::GradientType g = fev.gradient(ipt);
          typename FID_::HessianType  h = fev.hessian(ipt);
          for(int i(0); i < n; ++i)
          {
            v.axpy(-sed.phi[i].value, lvec[i]);
            g.add_outer_product(lvec[i], sed.phi[i].grad, -DT_(1));
            h.add_vec_mat_outer_product(lvec[i], sed.phi[i].hess, -DT_(1));
          }

          fid.add_value(omega, v);
          fid.add_grad(omega, g);
          fid.add_hess(omega, h);
        }
      };

      /// helper: checks whether function and vector both represent scalar or blocked functions
      template<typename Fun_, typename Vec_>
      struct ErrCompatHelper
      {
        // invalid combination of function and vector
        static constexpr bool valid = false;
      };

      template<typename Fun_, typename Mem_, typename DT_, typename IT_>
      struct ErrCompatHelper<Fun_, LAFEM::DenseVector<Mem_, DT_, IT_>>
      {
        // scalar vector, so function must also be scalar
        static constexpr bool valid = Fun_::ImageType::is_scalar;
      };

      template<typename Fun_, typename Mem_, typename DT_, typename IT_, int dim_>
      struct ErrCompatHelper<Fun_, LAFEM::DenseVectorBlocked<Mem_, DT_, IT_, dim_>>
      {
        typedef typename Fun_::ImageType ImageType;
        // blocked vector, so function must be a vector field of same dimension
        static constexpr bool valid = ImageType::is_vector && (ImageType::image_dim == dim_);
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Assembly job for the integration of an analytic function
     *
     * \tparam DataType_
     * The (scalar) datatype in which the assembly is to be performed.
     *
     * \tparam Function_
     * The analytic function class that is to be integrated.
     *
     * \tparam Trafo_
     * The trafo on which to integrate.
     *
     * \tparam max_der_
     * The maximum derivative order of the integrals that are to be assembled.
     *
     * \author Peter Zajac
     */
    template<typename DataType_, typename Function_, typename Trafo_, int max_der_>
    class AnalyticFunctionIntegralJob
    {
    public:
      static_assert(Function_::can_value, "function cannot compute values");
      static_assert(Function_::can_grad || (max_der_ < 1), "function gradients are required but not available");
      static_assert(Function_::can_hess || (max_der_ < 2), "function hessians are required but not available");

      /// the datatype to be used by the assembly
      typedef DataType_ DataType;

      /// declare our analytic eval traits
      typedef Analytic::EvalTraits<DataType, Function_> AnalyticEvalTraits;

      /// our function integral type
      typedef typename AnalyticFunctionIntegral<DataType, Function_>::Type FunctionIntegralType;

    public:
      /**
       * \brief Analytic Function Integral Task implementation
       */
      class Task
      {
      public:
        /// this task doesn't need to scatter
        static constexpr bool need_scatter = false;
        /// this task needs to combine
        static constexpr bool need_combine = true;

      protected:
        static constexpr TrafoTags trafo_config = TrafoTags::img_point | TrafoTags::jac_det;
        /// trafo evaluator
        typedef typename Trafo_::template Evaluator<typename Trafo_::ShapeType, DataType>::Type TrafoEvaluator;
        /// trafo evaluator
        TrafoEvaluator trafo_eval;
        /// trafo eval data type
        typename TrafoEvaluator::template ConfigTraits<trafo_config>::EvalDataType trafo_data;
        /// the cubature rule
        typename Assembly::Intern::CubatureTraits<TrafoEvaluator>::RuleType cubature_rule;
        /// the function evaluator
        typename Function_::template Evaluator<AnalyticEvalTraits> func_eval;
        /// the local integrals
        FunctionIntegralType loc_integral;
        /// the accumulated job integral
        FunctionIntegralType& job_integral;

      public:
        explicit Task(AnalyticFunctionIntegralJob& job) :
          trafo_eval(job._trafo),
          trafo_data(),
          cubature_rule(Cubature::ctor_factory, job._cubature_factory),
          func_eval(job._function),
          loc_integral(),
          job_integral(job._integral)
        {
        }

        void prepare(Index cell)
        {
          trafo_eval.prepare(cell);
        }

        void assemble()
        {
          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(k));

            // call helper to integrate
            Intern::AnaFunIntJobHelper<max_der_>::work(loc_integral,
              cubature_rule.get_weight(k) * trafo_data.jac_det, func_eval, trafo_data.img_point);
          }
        }

        void scatter()
        {
          // nothing to do here
        }

        void finish()
        {
          trafo_eval.finish();
        }

        void combine()
        {
          job_integral.push(loc_integral);
        }
      }; // class Task

    protected:
      /// the trafo to use
      const Trafo_& _trafo;
      /// the function to be integrated
      const Function_& _function;
      /// the cubature factory
      Cubature::DynamicFactory _cubature_factory;
      /// the function integral
      FunctionIntegralType _integral;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] function
       * A reference to the analytic function that is to be integrated.
       *
       * \param[in] trafo
       * A reference to the trafo representing the domain over which to integrate.
       *
       * \param[in] cubature
       * The name of the cubature rule to use for integration.
       */
      explicit AnalyticFunctionIntegralJob(const Function_& function, const Trafo_& trafo, const String& cubature) :
        _trafo(trafo),
        _function(function),
        _cubature_factory(cubature),
        _integral()
      {
        _integral.max_der = max_der_;
      }

      /**
       * \brief Returns the assembled function integral info
       */
      const FunctionIntegralType& result() const
      {
        return _integral;
      }
    }; // class AnalyticFunctionIntegralJob<...>

    /**
     * \brief Assembly job for the integration of a discrete finite element function
     *
     * \tparam Vector_
     * The type of the coefficient vector of the discrete function. May be one of the following:
     * - LAFEM::DenseVector<...>
     * - LAFEM::DenseVectorBlocked<...>
     *
     * \tparam Space_
     * The finite element space used for the discretization.
     *
     * \tparam max_der_
     * The maximum derivative order of the integrals that are to be assembled.
     *
     * \author Peter Zajac
     */
    template<typename Vector_, typename Space_, int max_der_>
    class DiscreteFunctionIntegralJob
    {
    public:
      /// the data-type of the vector
      typedef typename Vector_::DataType DataType;
      /// the value-type of the vector
      typedef typename Vector_::ValueType ValueType;

      typedef typename DiscreteFunctionIntegral<Vector_, Space_>::Type FunctionIntegralType;

    public:
      class Task
      {
      public:
        /// this task doesn't need to scatter
        static constexpr bool need_scatter = false;
        /// this task needs to combine
        static constexpr bool need_combine = true;

      protected:
        static constexpr SpaceTags space_tags = SpaceTags::value |
          (max_der_ >= 1 ? SpaceTags::grad : SpaceTags::none) |
          (max_der_ >= 2 ? SpaceTags::hess : SpaceTags::none);

        /// our assembly traits
        typedef Assembly::AsmTraits1<DataType, Space_, TrafoTags::jac_det, space_tags> AsmTraits;
        /// the vector that is to be integrated
        const Vector_& vector;
        /// the finite element space to be used
        const Space_& space;
        /// the cubature factory used for integration
        const typename AsmTraits::TrafoType& trafo;
        /// the trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval;
        /// the space evaluator
        typename AsmTraits::SpaceEvaluator space_eval;
        /// the space dof-mapping
        typename AsmTraits::DofMapping dof_mapping;
        /// the cubature rule used for integration
        typename AsmTraits::CubatureRuleType cubature_rule;
        /// the trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;
        /// the space evaluation data
        typename AsmTraits::SpaceEvalData space_data;
        /// the local vector to be assembled
        typename AsmTraits::template TLocalVector<ValueType> local_vector;
        /// the vector gather object
        typename Vector_::GatherAxpy gather_axpy;
        /// the local integral
        FunctionIntegralType loc_integral;
        /// the function value
        FunctionIntegralType& job_integral;

      public:
        explicit Task(DiscreteFunctionIntegralJob& job) :
          vector(job._vector),
          space(job._space),
          trafo(space.get_trafo()),
          trafo_eval(trafo),
          space_eval(space),
          dof_mapping(space),
          cubature_rule(Cubature::ctor_factory, job._cubature_factory),
          trafo_data(),
          space_data(),
          local_vector(),
          gather_axpy(vector),
          loc_integral(),
          job_integral(job._integral)
        {
        }

        void prepare(Index cell)
        {
          // prepare dof mapping
          dof_mapping.prepare(cell);

          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);
        }

        void assemble()
        {
          // format local vector
          local_vector.format();

          // gather local vector data
          gather_axpy(local_vector, dof_mapping);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(k));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // do the dirty work
            Intern::DiscFunIntJobHelper<max_der_>::work(loc_integral,
              cubature_rule.get_weight(k) * trafo_data.jac_det, space_data, local_vector, num_loc_dofs);
          }
        }

        void scatter()
        {
          // nothing to do here
        }

        void finish()
        {
          trafo_eval.finish();
        }

        void combine()
        {
          job_integral.push(loc_integral);
        }
      }; // class Task

    protected:
      /// the vector that is to be integrated
      const Vector_& _vector;
      /// the finite element space to be used
      const Space_& _space;
      /// the cubature factory
      Cubature::DynamicFactory _cubature_factory;
      /// the function integral
      FunctionIntegralType _integral;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] vector
       * A reference to the coefficient vector of the discrete function.
       *
       * \param[in] space
       * A reference to the finite element space used for the discretization
       *
       * \param[in] cubature
       * The name of the cubature rule to use for integration.
       */
      explicit DiscreteFunctionIntegralJob(const Vector_& vector, const Space_& space, const String& cubature) :
        _vector(vector),
        _space(space),
        _cubature_factory(cubature),
        _integral()
      {
        _integral.max_der = max_der_;
      }

      // \returns The integral of the discrete function
      const FunctionIntegralType& result() const
      {
        return _integral;
      }
    }; // class DiscreteFunctionIntegralJob<...>

    /**
     * \brief Assembly job for the integration of a analytic vs discrete error function
     *
     * This class is used to integrate the error function (u - u_h), where u is an analytic function
     * and u_h is a discrete finite element function.
     *
     * \tparam Function_
     * The type of the analytic function u that is to be integrated.
     *
     * \tparam Vector_
     * The type of the coefficient vector of the discrete function u_h. May be one of the following:
     * - LAFEM::DenseVector<...>
     * - LAFEM::DenseVectorBlocked<...>
     *
     * \tparam Space_
     * The finite element space used for the discretization.
     *
     * \tparam max_der_
     * The maximum derivative order of the integrals that are to be assembled.
     *
     * \author Peter Zajac
     */
    template<typename Function_, typename Vector_, typename Space_, int max_der_>
    class ErrorFunctionIntegralJob
    {
    public:
      /// the data-type of the vector
      typedef typename Vector_::DataType DataType;
      /// the value-type of the vector
      typedef typename Vector_::ValueType ValueType;

      /// make sure that function and vector are compatible
      static_assert(Intern::ErrCompatHelper<Function_, Vector_>::valid, "function and vector are incompatible");

      /// our function integral type
      typedef typename DiscreteFunctionIntegral<Vector_, Space_>::Type FunctionIntegralType;

      /// declare our analytic eval traits
      typedef Analytic::EvalTraits<DataType, Function_> AnalyticEvalTraits;

    public:
      class Task
      {
      public:
        /// this task doesn't need to scatter
        static constexpr bool need_scatter = false;
        /// this task needs to combine
        static constexpr bool need_combine = true;

      protected:
        static constexpr TrafoTags trafo_config = TrafoTags::img_point | TrafoTags::jac_det;
        static constexpr SpaceTags space_tags = SpaceTags::value |
          (max_der_ >= 1 ? SpaceTags::grad : SpaceTags::none) |
          (max_der_ >= 2 ? SpaceTags::hess : SpaceTags::none);

        /// our assembly traits
        typedef Assembly::AsmTraits1<DataType, Space_, trafo_config, space_tags> AsmTraits;
        /// the vector that is to be integrated
        const Vector_& vector;
        /// the finite element space to be used
        const Space_& space;
        /// the cubature factory used for integration
        const typename AsmTraits::TrafoType& trafo;
        /// the trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval;
        /// the space evaluator
        typename AsmTraits::SpaceEvaluator space_eval;
        /// the space dof-mapping
        typename AsmTraits::DofMapping dof_mapping;
        /// the cubature rule used for integration
        typename AsmTraits::CubatureRuleType cubature_rule;
        /// the trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;
        /// the space evaluation data
        typename AsmTraits::SpaceEvalData space_data;
        /// the local vector to be assembled
        typename AsmTraits::template TLocalVector<ValueType> local_vector;
        /// the vector gather object
        typename Vector_::GatherAxpy gather_axpy;
        /// the function evaluator
        typename Function_::template Evaluator<AnalyticEvalTraits> func_eval;
        /// the local integral
        FunctionIntegralType loc_integral;
        /// the function value
        FunctionIntegralType& job_integral;

      public:
        explicit Task(ErrorFunctionIntegralJob& job) :
          vector(job._vector),
          space(job._space),
          trafo(space.get_trafo()),
          trafo_eval(trafo),
          space_eval(space),
          dof_mapping(space),
          cubature_rule(Cubature::ctor_factory, job._cubature_factory),
          trafo_data(),
          space_data(),
          local_vector(),
          gather_axpy(vector),
          func_eval(job._function),
          loc_integral(),
          job_integral(job._integral)
        {
        }

        void prepare(Index cell)
        {
          // prepare dof mapping
          dof_mapping.prepare(cell);

          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);
        }

        void assemble()
        {
          // format local vector
          local_vector.format();

          // gather local vector data
          gather_axpy(local_vector, dof_mapping);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(k));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // do the dirty work
            Intern::ErrFunIntJobHelper<max_der_>::work(loc_integral,
              cubature_rule.get_weight(k) * trafo_data.jac_det, func_eval, trafo_data.img_point,
              space_data, local_vector, num_loc_dofs);
          }
        }

        void scatter()
        {
          // nothing to do here
        }

        void finish()
        {
          trafo_eval.finish();
        }

        void combine()
        {
          job_integral.push(loc_integral);
        }
      }; // class Task

    protected:
      /// the function to be integrated
      const Function_& _function;
      /// the vector that is to be integrated
      const Vector_& _vector;
      /// the finite element space to be used
      const Space_& _space;
      /// the cubature factory
      Cubature::DynamicFactory _cubature_factory;
      /// the function integral
      FunctionIntegralType _integral;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] function
       * A reference to the analytic function u.
       *
       * \param[in] vector
       * A reference to the coefficient vector of the discrete function u_h.
       *
       * \param[in] space
       * A reference to the finite element space used for the discretization of u_h
       *
       * \param[in] cubature
       * The name of the cubature rule to use for integration.
       */
      explicit ErrorFunctionIntegralJob(const Function_& function, const Vector_& vector, const Space_& space, const String& cubature) :
        _function(function),
        _vector(vector),
        _space(space),
        _cubature_factory(cubature),
        _integral()
      {
        _integral.max_der = max_der_;
      }

      // \returns The integral of the discrete function
      const FunctionIntegralType& result() const
      {
        return _integral;
      }
    }; // class ErrorFunctionIntegralJob<...>
  } // namespace Assembly
} // namespace FEAT
#endif // KERNEL_ASSEMBLY_FUNCTION_INTEGRAL_JOBS_HPP
