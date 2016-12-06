#ifndef FEAT_KERNEL_GEOMETRY_MESH_QUALITY_HEURISTIC
#define FEAT_KERNEL_GEOMETRY_MESH_QUALITY_HEURISTIC 1

#include <kernel/base_header.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/tiny_algebra.hpp>

namespace FEAT
{
  namespace Geometry
  {

    /**
     * \brief Helper class for computing heuristic mesh quality
     *
     * \tparam Shape_
     * Shape of the cells to compute the quality for
     *
     * The exact type of quality indicator depends on this. This class is an empty generic template, all implementations
     * are specialisations in Shape_.
     *
     * \author Jordi Paul
     */
    template<typename Shape_>
    struct MeshQualityHeuristic
    {
      /**
       * \brief Returns a descriptive String
       *
       * \returns What quality indicator the class computes as String.
       */
      static String description()
      {
        return String("UNDEFINED");
      }

      /**
       * \brief Computes minimum cell quality
       *
       * \tparam IdxType_
       * The index set type
       *
       * \tparam VtxType_
       * The vertex set type
       *
       * \param[in] idx
       * The vertex at cell index set of the mesh.
       *
       * \param[in] idx
       * The vertex set of the mesh.
       *
       * \returns A value between 0 and 1, 0 meaning degenerated cells and 1 means only pretty cells.
       *
       */
      template<typename IdxType_, typename VtxType_>
      static typename VtxType_::CoordType compute(const IdxType_& idx, const VtxType_& vtx);
    };

    /// \cond internal
    template<>
    struct MeshQualityHeuristic<Shape::Simplex<2>>
    {

      static String description()
      {
        return String("max( h(T) / vol(T)^1/2):");
      }

      template<typename IdxType_, typename VtxType_>
      static typename VtxType_::CoordType angle(const IdxType_& idx, const VtxType_& vtx,
      typename VtxType_::CoordType* cell_worst_angle = nullptr)
      {

        typedef typename VtxType_::CoordType CoordType;

        CoordType min_angle(Math::huge<CoordType>());
        CoordType max_angle(0);
        for(Index cell(0); cell < idx.get_num_entities(); ++cell)
        {
          Index i0(idx(cell,Index(0)));
          Index i1(idx(cell,Index(1)));
          Index i2(idx(cell,Index(2)));

          auto v1 = vtx[i1] - vtx[i0];
          auto v2 = vtx[i2] - vtx[i0];
          auto v3 = vtx[i2] - vtx[i1];

          CoordType n1(v1.norm_euclid());
          CoordType n2(v2.norm_euclid());
          CoordType n3(v3.norm_euclid());

          CoordType a0 = Math::acos(Tiny::dot(v1,v2)/(n1*n2));
          CoordType a1 = Math::acos( - Tiny::dot(v3,v1)/(n3*n1));
          CoordType a2 = Math::acos( Tiny::dot(v3,v2)/(n3*n2));

          ASSERT(Math::abs(a0+a1+a2-Math::pi<CoordType>()) < Math::sqrt(Math::eps<CoordType>()));

          CoordType my_min_angle(Math::min(a2,Math::min(a0,a1)));
          CoordType my_max_angle(Math::max(a2,Math::max(a0,a1)));

          if(cell_worst_angle != nullptr)
          {
            cell_worst_angle[cell] = CoordType(360)/(CoordType(2)*Math::pi<CoordType>())
              *Math::min(my_min_angle, Math::abs(Math::pi<CoordType>() - my_max_angle));
          }

          min_angle = Math::min(my_min_angle, min_angle);
          max_angle = Math::max(my_max_angle, max_angle);
        }

        CoordType worst_angle(Math::min(min_angle, Math::abs(Math::pi<CoordType>() - max_angle)));

        return worst_angle*(CoordType(360))/(CoordType(2)*Math::pi<CoordType>());
      }

      template<typename IdxType_, typename VtxType_>
      static void compute(typename VtxType_::CoordType& qual_min, typename VtxType_::CoordType& qual_sum,
      const IdxType_& idx, const VtxType_& vtx, typename VtxType_::CoordType* cell_qual = nullptr)
      {
        typedef typename VtxType_::CoordType CoordType;

        CoordType rho_min(Math::huge<CoordType>());
        CoordType rho_sum(0);
        CoordType diam(0);
        CoordType vol(0);

        Tiny::Matrix<CoordType, VtxType_::num_coords, VtxType_::num_coords> A(CoordType(0));

        for(Index cell(0); cell < idx.get_num_entities(); ++cell)
        {
          diam = CoordType(0);

          Index i0(idx(cell,Index(0)));
          Index i1(idx(cell,Index(1)));
          Index i2(idx(cell,Index(2)));

          A[0] = vtx[i1] - vtx[i0];
          A[1] = vtx[i2] - vtx[i0];

          vol = A.vol()/CoordType(2);

          CoordType h0((vtx[i2] - vtx[i0]).norm_euclid());
          diam = Math::max(diam, h0);

          CoordType h1(A[1].norm_euclid());
          diam = Math::max(diam, h1);

          CoordType h2(A[0].norm_euclid());
          diam = Math::max(diam, h2);

          CoordType my_rho(Math::sqrt(vol)/diam);

          if(cell_qual != nullptr)
          {
            cell_qual[cell] = Math::sqrt(CoordType(4)/Math::sqrt(CoordType(3)))*my_rho;
          }

          rho_sum += my_rho;
          rho_min = Math::min(rho_min, my_rho);
        }
        // We want the quality to go to zero if rho_max goes to infinity, so we take 1/rho_max. The absolute minimum
        // of rho is for the triangle with all angles = 60 degrees, so we scale with its rho for normalisation.
        qual_min = Math::sqrt(CoordType(4)/Math::sqrt(CoordType(3)))*rho_min;
        qual_sum = Math::sqrt(CoordType(4)/Math::sqrt(CoordType(3)))*rho_sum;
      }

    };

    template<>
    struct MeshQualityHeuristic<Shape::Simplex<3>>
    {

      static String description()
      {
        return String("max( h(T) / vol(T)^1/2):");
      }

      template<typename IdxType_, typename VtxType_>
      static typename VtxType_::CoordType angle(const IdxType_& idx, const VtxType_& vtx,
      typename VtxType_::CoordType* cell_worst_angle = nullptr)
      {

        typedef typename VtxType_::CoordType CoordType;

        Tiny::Matrix<CoordType, 6, VtxType_::num_coords> E(CoordType(0));
        CoordType edgelengths[6];
        CoordType angle[6];

        CoordType min_angle(Math::huge<CoordType>());
        CoordType max_angle(0);
        for(Index cell(0); cell < idx.get_num_entities(); ++cell)
        {

          CoordType my_min_angle(Math::huge<CoordType>());
          CoordType my_max_angle(0);

          Index i0(idx(cell,Index(0)));
          Index i1(idx(cell,Index(1)));
          Index i2(idx(cell,Index(2)));
          Index i3(idx(cell,Index(3)));

          E[0] = vtx[i1] - vtx[i0];
          E[1] = vtx[i2] - vtx[i0];
          E[2] = vtx[i3] - vtx[i0];
          E[3] = vtx[i2] - vtx[i1];
          E[4] = vtx[i3] - vtx[i1];
          E[5] = vtx[i3] - vtx[i2];

          for(int i(0); i < 6; ++i)
          {
            edgelengths[i] = E[i].norm_euclid();
          }

          angle[0] = Math::acos(Tiny::dot(E[0],E[1])/(edgelengths[0]*edgelengths[1]));
          angle[1] = Math::acos(Tiny::dot(E[0],E[2])/(edgelengths[0]*edgelengths[2]));
          angle[2] = Math::acos(Tiny::dot(E[1],E[2])/(edgelengths[1]*edgelengths[2]));
          angle[3] = Math::acos(Tiny::dot(E[3],E[4])/(edgelengths[3]*edgelengths[4]));
          angle[4] = Math::acos(Tiny::dot(E[4],E[5])/(edgelengths[4]*edgelengths[5]));
          angle[5] = Math::acos(Tiny::dot(E[0],E[4])/(edgelengths[0]*edgelengths[4]));

          for(int i(0); i < 6; ++i)
          {
            my_min_angle = Math::min(my_min_angle, angle[i]);
            my_max_angle = Math::max(my_max_angle, angle[i]);
          }

          min_angle = Math::min(my_min_angle, min_angle);
          max_angle = Math::max(my_max_angle, max_angle);

          if(cell_worst_angle != nullptr)
          {
            cell_worst_angle[cell] = CoordType(360)/(CoordType(2)*Math::pi<CoordType>())
              *Math::min(my_min_angle, Math::abs(Math::pi<CoordType>() - my_max_angle));
          }

        }

        CoordType worst_angle(Math::min(min_angle, Math::abs(Math::pi<CoordType>() - max_angle)));

        return worst_angle*(CoordType(360))/(CoordType(2)*Math::pi<CoordType>());
      }

      template<typename IdxType_, typename VtxType_>
      static void compute(typename VtxType_::CoordType& qual_min, typename VtxType_::CoordType& qual_sum,
      const IdxType_& idx, const VtxType_& vtx, typename VtxType_::CoordType* cell_qual = nullptr)
      {
        typedef typename VtxType_::CoordType CoordType;
        const CoordType vol_ref(Math::sqrt(CoordType(3))/CoordType(4)*Math::sqrt(CoordType(6))/CoordType(3*3));
        const CoordType fac(Math::pow(vol_ref, -CoordType(1)/CoordType(3)));

        CoordType rho_min(Math::huge<CoordType>());
        CoordType rho_sum(0);
        CoordType diam(0);
        CoordType vol(0);

        Tiny::Matrix<CoordType, VtxType_::num_coords, VtxType_::num_coords> A(CoordType(0));
        Tiny::Matrix<CoordType, 6, VtxType_::num_coords> E(CoordType(0));
        CoordType edgelengths[6];

        for(Index cell(0); cell < idx.get_num_entities(); ++cell)
        {
          diam = CoordType(0);

          Index i0(idx(cell,Index(0)));
          Index i1(idx(cell,Index(1)));
          Index i2(idx(cell,Index(2)));
          Index i3(idx(cell,Index(3)));

          A[0] = vtx[i1] - vtx[i0];
          A[1] = vtx[i2] - vtx[i0];
          A[2] = vtx[i3] - vtx[i0];

          vol = A.det()/CoordType(6);

          E[0] = vtx[i1] - vtx[i0];
          E[1] = vtx[i2] - vtx[i0];
          E[2] = vtx[i3] - vtx[i0];
          E[3] = vtx[i2] - vtx[i1];
          E[4] = vtx[i3] - vtx[i1];
          E[5] = vtx[i3] - vtx[i2];

          for(int i(0); i < 6; ++i)
          {
            edgelengths[i] = E[i].norm_euclid();
            diam = Math::max(diam, edgelengths[i]);
          }

          CoordType my_rho(Math::pow(vol, CoordType(1)/CoordType(3))/diam);

          if(cell_qual != nullptr)
          {
            cell_qual[cell] = fac*my_rho;
          }

          rho_sum += my_rho;
          rho_min = Math::min(rho_min, my_rho);
        }
        // We want the quality to go to zero if rho_max goes to infinity, so we take 1/rho_max. The absolute minimum
        // of rho is for the triangle with all angles = 60 degrees, so we scale with its rho for normalisation.
        qual_min = fac*rho_min;
        qual_sum = fac*rho_sum;
      }

    };

    template<>
    struct MeshQualityHeuristic<Shape::Hypercube<2>>
    {

      static String description()
      {
        return String("h_min/h_max*sqrt(1 - cos(alpha_min)):");
      }

      template<typename IdxType_, typename VtxType_>
      static typename VtxType_::CoordType angle(const IdxType_& idx, const VtxType_& vtx,
      typename VtxType_::CoordType* cell_worst_angle = nullptr)
      {

        typedef typename VtxType_::CoordType CoordType;

        CoordType min_angle(Math::huge<CoordType>());
        CoordType max_angle(0);
        for(Index cell(0); cell < idx.get_num_entities(); ++cell)
        {

          CoordType my_min_angle(Math::huge<CoordType>());
          CoordType my_max_angle(0);

          Index i0(idx(cell,Index(0)));
          Index i1(idx(cell,Index(1)));
          Index i2(idx(cell,Index(2)));
          Index i3(idx(cell,Index(3)));

          auto v0 = vtx[i1] - vtx[i0];
          auto v1 = vtx[i3] - vtx[i2];
          auto v2 = vtx[i2] - vtx[i0];
          auto v3 = vtx[i3] - vtx[i1];

          CoordType h0(v0.norm_euclid());
          CoordType h1(v1.norm_euclid());
          CoordType h2(v2.norm_euclid());
          CoordType h3(v3.norm_euclid());

          CoordType a0 = Math::acos(Tiny::dot(v2,v0)/(h2*h0));
          my_min_angle = Math::min(my_min_angle, Math::abs(a0));
          my_max_angle = Math::max(my_max_angle, Math::abs(a0));

          CoordType a1 = Math::acos( - Tiny::dot(v0,v3)/(h0*h3));
          my_min_angle = Math::min(my_min_angle, Math::abs(a1));
          my_max_angle = Math::max(my_max_angle, Math::abs(a1));

          CoordType a2 = Math::acos( - Tiny::dot(v1,v2)/(h1*h2));
          my_min_angle = Math::min(my_min_angle, Math::abs(a2));
          my_max_angle = Math::max(my_max_angle, Math::abs(a2));

          CoordType a3 = Math::acos( Tiny::dot(v3,v1)/(h3*h1));
          my_min_angle = Math::min(my_min_angle, Math::abs(a3));
          my_max_angle = Math::max(my_max_angle, Math::abs(a3));

          min_angle = Math::min(my_min_angle, min_angle);
          max_angle = Math::max(my_max_angle, max_angle);

          if(cell_worst_angle != nullptr)
          {
            cell_worst_angle[cell] = CoordType(360)/(CoordType(2)*Math::pi<CoordType>())
              *Math::min(my_min_angle, Math::abs(Math::pi<CoordType>() - my_max_angle));
          }

          ASSERT(Math::abs(a0+a1+a2+a3-CoordType(2)*Math::pi<CoordType>()) < Math::sqrt(Math::eps<CoordType>()));
        }

        CoordType worst_angle(Math::min(min_angle, Math::abs(Math::pi<CoordType>() - max_angle)));

        return worst_angle*CoordType(360)/(CoordType(2)*Math::pi<CoordType>());
      }

      template<typename IdxType_, typename VtxType_>
      static void compute(typename VtxType_::CoordType& qual_min, typename VtxType_::CoordType& qual_sum,
      const IdxType_& idx, const VtxType_& vtx, typename VtxType_::CoordType* cell_qual = nullptr)
      {

        typedef typename VtxType_::CoordType CoordType;

        qual_min = Math::huge<CoordType>();
        qual_sum = CoordType(0);

        for(Index cell(0); cell < idx.get_num_entities(); ++cell)
        {
          CoordType gamma(0);
          CoordType angle(Math::huge<CoordType>());

          Index i0(idx(cell,Index(0)));
          Index i1(idx(cell,Index(1)));
          Index i2(idx(cell,Index(2)));
          Index i3(idx(cell,Index(3)));

          auto v0 = vtx[i1] - vtx[i0];
          auto v1 = vtx[i3] - vtx[i2];
          auto v2 = vtx[i2] - vtx[i0];
          auto v3 = vtx[i3] - vtx[i1];

          CoordType h_min(Math::huge<CoordType>());
          CoordType h_max(Math::eps<CoordType>());

          CoordType h0(v0.norm_euclid());
          h_min = Math::min(h_min, h0);
          h_max = Math::max(h_max, h0);

          CoordType h1(v1.norm_euclid());
          h_min = Math::min(h_min, h1);
          h_max = Math::max(h_max, h1);

          CoordType h2(v2.norm_euclid());
          h_min = Math::min(h_min, h2);
          h_max = Math::max(h_max, h2);

          CoordType h3(v3.norm_euclid());
          h_min = Math::min(h_min, h3);
          h_max = Math::max(h_max, h3);

          CoordType a0 = Math::acos(Tiny::dot(v2,v0)/(h2*h0));
          angle = Math::min(angle, Math::abs(a0));
          gamma = Math::max(gamma, Tiny::dot(v2,v0)/(h2*h0));

          CoordType a1 = Math::acos( - Tiny::dot(v0,v3)/(h0*h3));
          angle = Math::min(angle, Math::abs(a1));
          gamma = Math::max(gamma, - Tiny::dot(v0,v3)/(h0*h3));

          CoordType a2 = Math::acos( - Tiny::dot(v1,v2)/(h1*h2));
          angle = Math::min(angle, Math::abs(a2));
          gamma = Math::max(gamma, - Tiny::dot(v1,v2)/(h1*h2));

          CoordType a3 = Math::acos( Tiny::dot(v3,v1)/(h3*h1));
          angle = Math::min(angle, Math::abs(a3));
          gamma = Math::max(gamma,  Tiny::dot(v3,v1)/(h3*h1));

          ASSERT(Math::abs(a0+a1+a2+a3-CoordType(2)*Math::pi<CoordType>()) < Math::sqrt(Math::eps<CoordType>()));

          CoordType my_qual(h_min/h_max*Math::sqrt(CoordType(1) - gamma));

          if(cell_qual != nullptr)
          {
            cell_qual[cell] = my_qual;
          }

          qual_min = Math::min(qual_min, my_qual);
          qual_sum += my_qual;

        }
      }
    };
    template<>
    struct MeshQualityHeuristic<Shape::Hypercube<3>>
    {


      static String description()
      {
        return String("h_min*min(det( nabla Phi))/(h_max*max(det nabla Phi)):");
      }

      template<typename IdxType_, typename VtxType_>
      static typename VtxType_::CoordType angle(const IdxType_& idx, const VtxType_& vtx,
      typename VtxType_::CoordType* cell_worst_angle = nullptr)
      {

        typedef typename VtxType_::CoordType CoordType;

        Tiny::Matrix<CoordType, 12, VtxType_::num_coords> E(CoordType(0));
        CoordType edgelengths[12];
        CoordType angle[24];

        CoordType min_angle(Math::huge<CoordType>());
        CoordType max_angle(0);
        for(Index cell(0); cell < idx.get_num_entities(); ++cell)
        {

          CoordType my_min_angle(Math::huge<CoordType>());
          CoordType my_max_angle(0);

          E[0] = vtx[idx(cell,Index(1))] - vtx[idx(cell,Index(0))];
          E[1] = vtx[idx(cell,Index(3))] - vtx[idx(cell,Index(2))];
          E[2] = vtx[idx(cell,Index(5))] - vtx[idx(cell,Index(4))];
          E[3] = vtx[idx(cell,Index(7))] - vtx[idx(cell,Index(6))];
          E[4] = vtx[idx(cell,Index(2))] - vtx[idx(cell,Index(0))];
          E[5] = vtx[idx(cell,Index(3))] - vtx[idx(cell,Index(1))];
          E[6] = vtx[idx(cell,Index(6))] - vtx[idx(cell,Index(4))];
          E[7] = vtx[idx(cell,Index(7))] - vtx[idx(cell,Index(5))];
          E[8] = vtx[idx(cell,Index(4))] - vtx[idx(cell,Index(0))];
          E[9] = vtx[idx(cell,Index(5))] - vtx[idx(cell,Index(1))];
          E[10] = vtx[idx(cell,Index(6))] - vtx[idx(cell,Index(2))];
          E[11] = vtx[idx(cell,Index(7))] - vtx[idx(cell,Index(3))];

          for(int i(0); i < 12; ++i)
          {
            edgelengths[i] = E[i].norm_euclid();
          }

          // Angles at vertex 0
          angle[0] = Math::acos(Tiny::dot(E[0],E[4])/(edgelengths[0]*edgelengths[4]));
          angle[1] = Math::acos(Tiny::dot(E[0],E[8])/(edgelengths[0]*edgelengths[8]));
          angle[2] = Math::acos(Tiny::dot(E[4],E[8])/(edgelengths[4]*edgelengths[8]));

          // Angles at vertex 1
          angle[3] = Math::acos(-Tiny::dot(E[0],E[5])/(edgelengths[0]*edgelengths[5]));
          angle[4] = Math::acos(-Tiny::dot(E[0],E[9])/(edgelengths[9]*edgelengths[9]));
          angle[5] = Math::acos( Tiny::dot(E[5],E[9])/(edgelengths[5]*edgelengths[9]));

          // Angles at vertex 2
          angle[6] = Math::acos(-Tiny::dot(E[1],E[4])/(edgelengths[1]*edgelengths[ 4]));
          angle[7] = Math::acos( Tiny::dot(E[1],E[10])/(edgelengths[1]*edgelengths[10]));
          angle[8] = Math::acos(-Tiny::dot(E[4],E[10])/(edgelengths[4]*edgelengths[10]));

          // Angles at vertex 3
          angle[9]  = Math::acos(-Tiny::dot(E[1],E[ 5])/(edgelengths[1]*edgelengths[ 5]));
          angle[10] = Math::acos(-Tiny::dot(E[1],E[11])/(edgelengths[1]*edgelengths[11]));
          angle[11] = Math::acos(-Tiny::dot(E[5],E[11])/(edgelengths[5]*edgelengths[11]));

          // Angles at vertex 4
          angle[12] = Math::acos( Tiny::dot(E[2],E[6])/(edgelengths[2]*edgelengths[6]));
          angle[13] = Math::acos(-Tiny::dot(E[2],E[8])/(edgelengths[2]*edgelengths[8]));
          angle[14] = Math::acos( Tiny::dot(E[6],E[8])/(edgelengths[6]*edgelengths[8]));

          // Angles at vertex 5
          angle[15] = Math::acos(-Tiny::dot(E[2],E[7])/(edgelengths[2]*edgelengths[7]));
          angle[16] = Math::acos(-Tiny::dot(E[2],E[9])/(edgelengths[2]*edgelengths[9]));
          angle[17] = Math::acos(-Tiny::dot(E[7],E[9])/(edgelengths[7]*edgelengths[9]));

          // Angles at vertex 6
          angle[18] = Math::acos(-Tiny::dot(E[3],E[ 6])/(edgelengths[3]*edgelengths[ 6]));
          angle[19] = Math::acos(-Tiny::dot(E[3],E[10])/(edgelengths[3]*edgelengths[10]));
          angle[20] = Math::acos(-Tiny::dot(E[6],E[10])/(edgelengths[6]*edgelengths[10]));

          // Angles at vertex 7
          angle[21] = Math::acos(Tiny::dot(E[3],E[ 7])/(edgelengths[3]*edgelengths[ 7]));
          angle[22] = Math::acos(Tiny::dot(E[3],E[11])/(edgelengths[3]*edgelengths[11]));
          angle[23] = Math::acos(Tiny::dot(E[7],E[11])/(edgelengths[7]*edgelengths[11]));

          for(int i(0); i < 24; ++i)
          {
            my_min_angle = Math::min(my_min_angle,angle[i]);
            my_max_angle = Math::max(my_max_angle,angle[i]);
          }

          min_angle = Math::min(my_min_angle, min_angle);
          max_angle = Math::max(my_max_angle, max_angle);

          if(cell_worst_angle != nullptr)
          {
            cell_worst_angle[cell] = CoordType(360)/(CoordType(2)*Math::pi<CoordType>())
              *Math::min(my_min_angle, Math::abs(Math::pi<CoordType>() - my_max_angle));
          }

        }

        CoordType worst_angle(Math::min(min_angle, Math::abs(Math::pi<CoordType>() - max_angle)));

        return worst_angle*(CoordType(360))/(CoordType(2)*Math::pi<CoordType>());
      }

      template<typename IdxType_, typename VtxType_>
      static void compute(typename VtxType_::CoordType& qual_min, typename VtxType_::CoordType& qual_sum,
      const IdxType_& idx, const VtxType_& vtx, typename VtxType_::CoordType* cell_qual = nullptr)
      {

        typedef typename VtxType_::CoordType CoordType;
        typedef typename VtxType_::ConstVertexReference ConstVertexReference;
        typedef typename VtxType_::VertexType DomainPointType;

        static constexpr int num_coords = VtxType_::num_coords;

        typedef Tiny::Matrix<CoordType, 3, num_coords> JacobiMatrixType;

        CoordType coeff[num_coords][8];
        DomainPointType xq(0);
        JacobiMatrixType jac_mat;
        Tiny::Matrix<CoordType, 12, VtxType_::num_coords> E(CoordType(0));
        CoordType edgelengths[12];

        qual_min = Math::huge<CoordType>();
        qual_sum = CoordType(0);

        for(Index cell(0); cell < idx.get_num_entities(); ++cell)
        {
          CoordType h_min(Math::huge<CoordType>());
          CoordType h_max(0);

          CoordType jac_det_min(Math::huge<CoordType>());
          CoordType jac_det_max(0);

          ConstVertexReference v0 = vtx[idx(cell, 0)];
          ConstVertexReference v1 = vtx[idx(cell, 1)];
          ConstVertexReference v2 = vtx[idx(cell, 2)];
          ConstVertexReference v3 = vtx[idx(cell, 3)];
          ConstVertexReference v4 = vtx[idx(cell, 4)];
          ConstVertexReference v5 = vtx[idx(cell, 5)];
          ConstVertexReference v6 = vtx[idx(cell, 6)];
          ConstVertexReference v7 = vtx[idx(cell, 7)];

          E[0] = v1 - v0;
          E[1] = v3 - v2;
          E[2] = v5 - v4;
          E[3] = v7 - v6;
          E[4] = v2 - v0;
          E[5] = v3 - v1;
          E[6] = v6 - v4;
          E[7] = v7 - v5;
          E[8] = v4 - v0;
          E[9] = v5 - v1;
          E[10] = v6 - v2;
          E[11] = v7 - v3;

          for(int i(0); i < 12; ++i)
          {
            edgelengths[i] = E[i].norm_euclid();
            h_min = Math::min(edgelengths[i], h_min);
            h_max = Math::max(edgelengths[i], h_max);
          }

          // calculate transformation coefficients
          // j = _coeff[i][j] for all j = 0....7
          // v = 0 + 1*x + 2*y + 3*z + x*y*4 + x*z*5 + y*z*6 + x*y*z*7
          for(int i(0); i < num_coords; ++i)
          {
            coeff[i][0] = CoordType(0.125) * CoordType( + v0[i] + v1[i] + v2[i] + v3[i] + v4[i] + v5[i] + v6[i] + v7[i]);
            coeff[i][1] = CoordType(0.125) * CoordType( - v0[i] + v1[i] - v2[i] + v3[i] - v4[i] + v5[i] - v6[i] + v7[i]);
            coeff[i][2] = CoordType(0.125) * CoordType( - v0[i] - v1[i] + v2[i] + v3[i] - v4[i] - v5[i] + v6[i] + v7[i]);
            coeff[i][3] = CoordType(0.125) * CoordType( - v0[i] - v1[i] - v2[i] - v3[i] + v4[i] + v5[i] + v6[i] + v7[i]);
            coeff[i][4] = CoordType(0.125) * CoordType( + v0[i] - v1[i] - v2[i] + v3[i] + v4[i] - v5[i] - v6[i] + v7[i]);
            coeff[i][5] = CoordType(0.125) * CoordType( + v0[i] - v1[i] + v2[i] - v3[i] - v4[i] + v5[i] - v6[i] + v7[i]);
            coeff[i][6] = CoordType(0.125) * CoordType( + v0[i] + v1[i] - v2[i] - v3[i] - v4[i] - v5[i] + v6[i] + v7[i]);
            coeff[i][7] = CoordType(0.125) * CoordType( - v0[i] + v1[i] + v2[i] - v3[i] + v4[i] - v5[i] - v6[i] + v7[i]);
          }

          for(int i(0); i < 8; ++i)
          {
            // Set point coords using bloody bitshifts
            for(int j(0); j < 3; ++j)
            {
              xq(j) = (CoordType(((i >> j) & 1) << 1) - CoordType(1)) * CoordType(1);
            }

            for(int j(0); j < jac_mat.m; ++j)
            {
              jac_mat(j,0) = coeff[j][1] + xq[1] * coeff[j][4] + xq[2] * (coeff[j][5] + xq[1] * coeff[j][7]);
              jac_mat(j,1) = coeff[j][2] + xq[0] * coeff[j][4] + xq[2] * (coeff[j][6] + xq[0] * coeff[j][7]);
              jac_mat(j,2) = coeff[j][3] + xq[0] * coeff[j][5] + xq[1] * (coeff[j][6] + xq[0] * coeff[j][7]);
            }

            CoordType jac_det(jac_mat.vol());
            jac_det_min = Math::min(jac_det, jac_det_min);
            jac_det_max = Math::max(jac_det, jac_det_max);

          }

          CoordType my_qual(h_min*jac_det_min/(h_max*jac_det_max));

          if(cell_qual != nullptr)
          {
            cell_qual[cell] = my_qual;
          }

          qual_min = Math::min(qual_min, my_qual);
          qual_sum += my_qual;

        }
      }
    };
    /// \endcond
  } // namespace Geometry
} // namespace FEAT
#endif // FEAT_KERNEL_GEOMETRY_MESH_QUALITY_HEURISTIC
