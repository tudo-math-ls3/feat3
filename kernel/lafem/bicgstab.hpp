#pragma once
#ifndef KERNEL_LAFEM_BICGSTAB_HPP
#define KERNEL_LAFEM_BICGSTAB_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/math.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/preconditioner.hpp>


namespace FEAST
{
  namespace LAFEM
  {
    template <typename Algo_>
    struct BiCGStab
    {
      template <typename MT_, typename VT_>
      static void value(VT_ & x, const MT_ & A, const VT_ & b, Preconditioner<Algo_, MT_, VT_> & precon, Index max_iters, typename VT_::DataType eps_relative, bool verbose = true)
      {
        typedef typename VT_::DataType DT_;

        DT_ defnorm, defnorm_0, defnorm_00(DT_(1e14));
        Index iter = 0;
        DT_ rho_tilde, rho_tilde_old, alpha_tilde, omega_tilde, beta_tilde, gamma_tilde;
        //bool early_exit = 0;
        bool restarted = false;
        bool converged = 0;

        VT_ r(A.create_vector_l());
        VT_ r_tilde(A.create_vector_l());
        VT_ r_tilde_0(A.create_vector_l());
        VT_ p_tilde(A.create_vector_l());
        VT_ v(A.create_vector_l());
        VT_ v_tilde(A.create_vector_l());
        VT_ s(A.create_vector_l());
        VT_ s_tilde(A.create_vector_l());
        VT_ t(A.create_vector_l());
        VT_ t_tilde(A.create_vector_l());

        do
        {
          //r.template defect<Algo_>(b, A, x);
          A.template apply<Algo_>(r, x, b, -DT_(1));
          defnorm_0 = r.template norm2<Algo_>();
          defnorm = defnorm_0;
          precon.apply(r_tilde_0, r);

          if (restarted == false)
          {
            defnorm_00 = defnorm_0;
          }
          r_tilde.copy(r_tilde_0);
          p_tilde.copy(r_tilde_0);

          rho_tilde = r_tilde_0.template dot<Algo_>(r_tilde_0);

          // main BiCGStab loop
          do
          {
            iter = iter + 1;

            //v.template product_matvec<Algo_>(A, p_tilde);
            A.template apply<Algo_>(v, p_tilde);
            precon.apply(v_tilde, v);

            gamma_tilde = v_tilde.template dot<Algo_>(r_tilde_0);

            if (Math::abs(gamma_tilde) < Math::abs(rho_tilde)*1e-14)
            {
              restarted = true;
              //std::cout << "Breakpoint 1" << std::endl;
              break;
            }

            alpha_tilde = rho_tilde / gamma_tilde;

            if ((Math::abs(alpha_tilde) * v_tilde.template norm2<Algo_>()) / defnorm < 1e-5)
            {
              restarted = true;;
              //std::cout << "Breakpoint 2" << std::endl;
              // \TODO warum ist das break hier nicht aktiv?
              //break;
            }

            DT_ malpha_tilde(-alpha_tilde);
            s.template axpy<Algo_>(v, r, malpha_tilde);

            defnorm = s.template norm2<Algo_>();
            if (defnorm < eps_relative * defnorm_00)
            {
              x.template axpy<Algo_>(p_tilde, x, alpha_tilde);

              //early_exit = 1;
              converged = 1;
              //std::cout << "Breakpoint 3 (converged)" << std::endl;
              break;
            }
            s_tilde.template axpy<Algo_>(v_tilde, r_tilde, malpha_tilde);

            //t.template product_matvec<Algo_>(A, s_tilde);
            A.template apply<Algo_>(t, s_tilde);

            precon.apply(t_tilde, t);

            gamma_tilde = t_tilde.template dot<Algo_>(t_tilde);
            omega_tilde = t_tilde.template dot<Algo_>(s_tilde);

            if (Math::abs(gamma_tilde) < Math::abs(omega_tilde) * 1e-14)
            {
              restarted = true;
              //std::cout << "Breakpoint 4" << std::endl;
              break;
            }
            omega_tilde = omega_tilde / gamma_tilde;

            x.template axpy<Algo_>(s_tilde, x, omega_tilde);
            x.template axpy<Algo_>(p_tilde, x, alpha_tilde);

            DT_ momega_tilde(-omega_tilde);
            r.template axpy<Algo_>(t, s, momega_tilde);

            defnorm = r.template norm2<Algo_>();
            if (defnorm < eps_relative * defnorm_00)
            {
              converged = 1;
              //std::cout << "Breakpoint 5 (converged)" << std::endl;
              break;
            }

            r_tilde.template axpy<Algo_>(t_tilde, s_tilde, momega_tilde);

            rho_tilde_old = rho_tilde;
            rho_tilde = r_tilde.template dot<Algo_>(r_tilde_0);

            beta_tilde = (alpha_tilde / omega_tilde) * (rho_tilde / rho_tilde_old);

            p_tilde.template axpy<Algo_>(v_tilde, p_tilde, momega_tilde);
            p_tilde.template scale<Algo_>(p_tilde, beta_tilde);
            p_tilde.template axpy<Algo_>(p_tilde, r_tilde);

          } while (iter <= max_iters);

        } while (converged == 0 && iter < max_iters);


        Index used_iters = 0;
        used_iters = iter + 1;
        if (verbose)
          std::cout<<"Initial defect norm: " << defnorm_00 << ", Final defect norm: " << defnorm << ", used iters: " << used_iters << std::endl;
      }

      /**
       * \brief Variant with filter
       *
       * \author Jordi Paul
       **/
      template <typename MT_, typename VT_, typename FT_>
      static void value(VT_ & x, const MT_ & A, const VT_ & b, const FT_ & filter, Preconditioner<Algo_, MT_, VT_> & precon, Index max_iters, typename VT_::DataType eps_relative, bool verbose = true)
      {
        typedef typename VT_::DataType DT_;

        DT_ defnorm, defnorm_0, defnorm_00(DT_(1e14));
        Index iter = 0;
        DT_ rho_tilde, rho_tilde_old, alpha_tilde, omega_tilde, beta_tilde, gamma_tilde;
        //bool early_exit = 0;
        bool restarted = false;
        bool converged = 0;

        VT_ r(A.create_vector_l());
        VT_ r_tilde(A.create_vector_l());
        VT_ r_tilde_0(A.create_vector_l());
        VT_ p_tilde(A.create_vector_l());
        VT_ v(A.create_vector_l());
        VT_ v_tilde(A.create_vector_l());
        VT_ s(A.create_vector_l());
        VT_ s_tilde(A.create_vector_l());
        VT_ t(A.create_vector_l());
        VT_ t_tilde(A.create_vector_l());

        do
        {
          //r.template defect<Algo_>(b, A, x);
          A.template apply<Algo_>(r, x, b, -DT_(1));
          filter.template filter_def<Algo_>(r);

          defnorm_0 = r.template norm2<Algo_>();
          defnorm = defnorm_0;
          precon.apply(r_tilde_0, r);
          filter.template filter_cor<Algo_>(r_tilde_0);

          if (restarted == false)
          {
            defnorm_00 = defnorm_0;
          }
          r_tilde.copy(r_tilde_0);
          p_tilde.copy(r_tilde_0);

          rho_tilde = r_tilde_0.template dot<Algo_>(r_tilde_0);

          // main BiCGStab loop
          do
          {
            iter = iter + 1;

            //v.template product_matvec<Algo_>(A, p_tilde);
            A.template apply<Algo_>(v, p_tilde);
            filter.template filter_def<Algo_>(v);
            precon.apply(v_tilde, v);

            gamma_tilde = v_tilde.template dot<Algo_>(r_tilde_0);

            if (Math::abs(gamma_tilde) < Math::abs(rho_tilde)*1e-14)
            {
              restarted = true;
              //std::cout << "Breakpoint 1" << std::endl;
              break;
            }

            alpha_tilde = rho_tilde / gamma_tilde;

            if ((Math::abs(alpha_tilde) * v_tilde.template norm2<Algo_>()) / defnorm < 1e-5)
            {
              restarted = true;;
              //std::cout << "Breakpoint 2" << std::endl;
              // \TODO warum ist das break hier nicht aktiv?
              //break;
            }

            DT_ malpha_tilde(-alpha_tilde);
            s.template axpy<Algo_>(v, r, malpha_tilde);

            defnorm = s.template norm2<Algo_>();
            if (defnorm < eps_relative * defnorm_00)
            {
              x.template axpy<Algo_>(p_tilde, x, alpha_tilde);

              //early_exit = 1;
              converged = 1;
              //std::cout << "Breakpoint 3 (converged)" << std::endl;
              break;
            }
            s_tilde.template axpy<Algo_>(v_tilde, r_tilde, malpha_tilde);

            //t.template product_matvec<Algo_>(A, s_tilde);
            A.template apply<Algo_>(t, s_tilde);
            filter.template filter_def<Algo_>(t);

            precon.apply(t_tilde, t);
            filter.template filter_cor<Algo_>(t_tilde);

            gamma_tilde = t_tilde.template dot<Algo_>(t_tilde);
            omega_tilde = t_tilde.template dot<Algo_>(s_tilde);

            if (Math::abs(gamma_tilde) < Math::abs(omega_tilde) * 1e-14)
            {
              restarted = true;
              //std::cout << "Breakpoint 4" << std::endl;
              break;
            }
            omega_tilde = omega_tilde / gamma_tilde;

            x.template axpy<Algo_>(s_tilde, x, omega_tilde);
            x.template axpy<Algo_>(p_tilde, x, alpha_tilde);

            DT_ momega_tilde(-omega_tilde);
            r.template axpy<Algo_>(t, s, momega_tilde);

            defnorm = r.template norm2<Algo_>();
            if (defnorm < eps_relative * defnorm_00)
            {
              converged = 1;
              //std::cout << "Breakpoint 5 (converged)" << std::endl;
              break;
            }

            r_tilde.template axpy<Algo_>(t_tilde, s_tilde, momega_tilde);

            rho_tilde_old = rho_tilde;
            rho_tilde = r_tilde.template dot<Algo_>(r_tilde_0);

            beta_tilde = (alpha_tilde / omega_tilde) * (rho_tilde / rho_tilde_old);

            p_tilde.template axpy<Algo_>(v_tilde, p_tilde, momega_tilde);
            p_tilde.template scale<Algo_>(p_tilde, beta_tilde);
            p_tilde.template axpy<Algo_>(p_tilde, r_tilde);

          } while (iter <= max_iters);

        } while (converged == 0 && iter < max_iters);


        Index used_iters = 0;
        used_iters = iter + 1;
        if (verbose)
          std::cout<<"Initial defect norm: " << defnorm_00 << ", Final defect norm: " << defnorm << ", used iters: " << used_iters << std::endl;
      }
    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_BICGSTAB_HPP
