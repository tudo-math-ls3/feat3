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
      static typename VtxType_::CoordType angle(const IdxType_& idx, const VtxType_& vtx)
      {

        typedef typename VtxType_::CoordType CoordType;

        CoordType angle(Math::huge<CoordType>());
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

          angle = Math::min(angle,Math::min(a2,Math::min(a0,a1)));
        }

        return angle*(CoordType(360))/(CoordType(2)*Math::pi<CoordType>());
      }

      template<typename IdxType_, typename VtxType_>
      static void compute(typename VtxType_::CoordType& qual_min, typename VtxType_::CoordType& qual_sum,
      const IdxType_& idx, const VtxType_& vtx)
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

          vol = Math::sqrt(A.vol());

          CoordType h0((vtx[i2] - vtx[i0]).norm_euclid());
          diam = Math::max(diam, h0);

          CoordType h1(A[1].norm_euclid());
          diam = Math::max(diam, h1);

          CoordType h2(A[0].norm_euclid());
          diam = Math::max(diam, h2);

          CoordType my_rho(vol/diam);

          rho_sum += my_rho;
          rho_min = Math::min(rho_min, my_rho);
        }
        // We want the quality to go to zero if rho_max goes to infinity, so we take 1/rho_max. The absolute minimum
        // of rho is for the triangle with all angles = 60 degrees, so we scale with its rho for normalisation.
        qual_min = Math::sqrt(CoordType(2)/Math::sqrt(CoordType(3)))*rho_min;
        qual_sum = Math::sqrt(CoordType(2)/Math::sqrt(CoordType(3)))*rho_sum;
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
      static typename VtxType_::CoordType angle(const IdxType_& idx, const VtxType_& vtx)
      {

        typedef typename VtxType_::CoordType CoordType;

        CoordType worst_angle(Math::huge<CoordType>());
        for(Index cell(0); cell < idx.get_num_entities(); ++cell)
        {

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
          worst_angle = Math::min(worst_angle, Math::abs(a0));

          CoordType a1 = Math::acos( - Tiny::dot(v0,v3)/(h0*h3));
          worst_angle = Math::min(worst_angle, Math::abs(a1));

          CoordType a2 = Math::acos( - Tiny::dot(v1,v2)/(h1*h2));
          worst_angle = Math::min(worst_angle, Math::abs(a2));

          CoordType a3 = Math::acos( Tiny::dot(v3,v1)/(h3*h1));
          worst_angle = Math::min(worst_angle, Math::abs(a3));

          ASSERT(Math::abs(a0+a1+a2+a3-CoordType(2)*Math::pi<CoordType>()) < Math::sqrt(Math::eps<CoordType>()));
        }
        return worst_angle*(CoordType(360))/(CoordType(2)*Math::pi<CoordType>());
      }

      template<typename IdxType_, typename VtxType_>
      static void compute(typename VtxType_::CoordType& qual_min, typename VtxType_::CoordType& qual_sum,
      const IdxType_& idx, const VtxType_& vtx)
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

          qual_min = Math::min(qual_min, my_qual);
          qual_sum += my_qual;

        }
      }
    };
    /// \endcond
  } // namespace Geometry
} // namespace FEAT
#endif // FEAT_KERNEL_GEOMETRY_MESH_QUALITY_HEURISTIC
