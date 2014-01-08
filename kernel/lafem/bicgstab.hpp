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
#include <kernel/lafem/algorithm.hpp>



namespace FEAST
{
  namespace LAFEM
  {
    template <typename Algo_>
    struct BiCGStab
    {
      template <typename MT_, typename VT_>
      static void value(VT_ & x, const MT_ & A, const VT_ & b, Preconditioner<Algo_, MT_, VT_> & precon, Index max_iters, typename VT_::DataType eps_relative)
      {
        typedef typename VT_::DataType DT_;
        typedef typename VT_::MemType Arch_;

        DT_ defnorm, defnorm_0, defnorm_00(1e14);
        Index iter = 0;
        DT_ rho_tilde, rho_tilde_old, alpha_tilde, omega_tilde, beta_tilde, gamma_tilde;
        //bool early_exit = 0;
        bool restarted = false;
        bool converged = 0;

        VT_ r(b.size());
        VT_ r_tilde(b.size());
        VT_ r_tilde_0(b.size());
        VT_ p_tilde(b.size());
        VT_ v(b.size());
        VT_ v_tilde(b.size());
        VT_ s(b.size());
        VT_ s_tilde(b.size());
        VT_ t(b.size());
        VT_ t_tilde(b.size());

        do
        {
          r.template defect<Algo_>(b, A, x);
          defnorm_0 = r.template norm2<Algo_>();
          defnorm = defnorm_0;
          precon.apply(r_tilde_0, r);

          if (restarted == false)
          {
            defnorm_00 = defnorm_0;
          }
          copy(r_tilde, r_tilde_0);
          copy(p_tilde, r_tilde_0);

          rho_tilde = r_tilde_0.template dot_product<Algo_>(r_tilde_0);

          // main BiCGStab loop
          do
          {
            iter = iter + 1;

            v.template product_matvec<Algo_>(A, p_tilde);
            precon.apply(v_tilde, v);

            gamma_tilde = v_tilde.template dot_product<Algo_>(r_tilde_0);

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
            s.template axpy<Algo_>(malpha_tilde, v, r);

            defnorm = s.template norm2<Algo_>();
            if (defnorm < eps_relative * defnorm_00)
            {
              x.template axpy<Algo_>(alpha_tilde, p_tilde, x);

              //early_exit = 1;
              converged = 1;
              //std::cout << "Breakpoint 3 (converged)" << std::endl;
              break;
            }
            s_tilde.template axpy<Algo_>(malpha_tilde, v_tilde, r_tilde);

            t.template product_matvec<Algo_>(A, s_tilde);

            precon.apply(t_tilde, t);

            gamma_tilde = t_tilde.template dot_product<Algo_>(t_tilde);
            omega_tilde = t_tilde.template dot_product<Algo_>(s_tilde);

            if (Math::abs(gamma_tilde) < Math::abs(omega_tilde) * 1e-14)
            {
              restarted = true;
              //std::cout << "Breakpoint 4" << std::endl;
              break;
            }
            omega_tilde = omega_tilde / gamma_tilde;

            x.template axpy<Algo_>(omega_tilde, s_tilde, x);
            x.template axpy<Algo_>(alpha_tilde, p_tilde, x);

            DT_ momega_tilde(-omega_tilde);
            r.template axpy<Algo_>(momega_tilde, t, s);

            defnorm = r.template norm2<Algo_>();
            if (defnorm < eps_relative * defnorm_00)
            {
              converged = 1;
              //std::cout << "Breakpoint 5 (converged)" << std::endl;
              break;
            }

            r_tilde.template axpy<Algo_>(momega_tilde, t_tilde, s_tilde);

            rho_tilde_old = rho_tilde;
            rho_tilde = r_tilde.template dot_product<Algo_>(r_tilde_0);

            beta_tilde = (alpha_tilde / omega_tilde) * (rho_tilde / rho_tilde_old);

            p_tilde.template axpy<Algo_>(momega_tilde, v_tilde, p_tilde);
            p_tilde.template scale<Algo_>(p_tilde, beta_tilde);
            p_tilde.template sum<Algo_>(p_tilde, r_tilde);

          } while (iter <= max_iters);

        } while (converged == 0 && iter < max_iters);


        Index used_iters = 0;
        used_iters = iter + 1;
        std::cout<<"Initial defect norm: " << defnorm_00 << ", Final defect norm: " << defnorm << ", used iters: " << used_iters << std::endl;
      }
    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_BICGSTAB_HPP
