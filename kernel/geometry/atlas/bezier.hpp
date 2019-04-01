// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_ATLAS_BEZIER_HPP
#define KERNEL_GEOMETRY_ATLAS_BEZIER_HPP 1

#include <kernel/geometry/atlas/chart.hpp>

namespace FEAT
{
  namespace Geometry
  {
    namespace Atlas
    {
      /// Bezier chart traits
      struct BezierTraits
      {
        /// we support explicit map
        static constexpr bool is_explicit = true;
        /// we don't support implicit project
        static constexpr bool is_implicit = true;
        /// this is a 2D object
        static constexpr int world_dim = 2;
        /// we have 1D parameters
        static constexpr int param_dim = 1;
      };

      /**
       * \brief Bezier chart class template
       *
       * This chart represents a 2D composite Bezier curve of varying degree, which is represented
       * by a vector of vertex points, Bezier control points and parameters.
       *
       * This chart implements both the implicit and explicit chart interfaces.
       *
       * \tparam Mesh_
       * The type of the mesh to be parameterised by this chart.
       *
       * \author Peter Zajac
       */
      template<typename Mesh_>
      class Bezier :
        public ChartCRTP<Bezier<Mesh_>, Mesh_, BezierTraits>
      {
      public:
        /// The CRTP base class
        typedef ChartCRTP<Bezier<Mesh_>, Mesh_, BezierTraits> BaseClass;
        /// Floating point type for coordinates
        typedef typename BaseClass::CoordType DataType;
        /// Vector type for world points aka. image points
        typedef typename BaseClass::WorldPoint WorldPoint;
        /// Vector type for parameter points aka. domain points
        typedef typename BaseClass::ParamPoint ParamPoint;

        /// maximum allowed spline degree
        static constexpr int max_degree = 5;

      protected:
        /// The vertex point pointer array
        std::deque<std::size_t> _vtx_ptr;
        /// The vertex and control points for the spline
        std::deque<WorldPoint> _world;
        /// The parameter values for these world points
        std::deque<ParamPoint> _param;
        /// Specifies whether the spline is closed
        bool _closed;
        /// Specifies if the chart is to be oriented negatively or not
        DataType _orientation;

      public:
        /// default CTOR
        explicit Bezier(bool closed = false, DataType orientation = DataType(1)) :
          BaseClass(),
          _closed(closed),
          _orientation(orientation)
        {
        }

        /**
         * \brief Constructor.
         *
         * \param[in] vtx_ptr
         * A deque containing the indices of the vertex points in the world point array.
         *
         * \param[in] world
         * A deque of vertex and control points representing the spline.
         *
         * \param[in] param
         * A deque of parameter points representing the parameter space.
         *
         * \param[in] closed
         * Specifies whether the spline is closed.
         */
        explicit Bezier(
          const std::deque<std::size_t>& vtx_ptr,
          const std::deque<WorldPoint>& world,
          const std::deque<ParamPoint>& param,
          bool closed, DataType orientation = DataType(1)) :
          BaseClass(),
          _vtx_ptr(vtx_ptr),
          _world(world),
          _param(param),
          _closed(closed),
          _orientation(orientation)
        {
          // we need at least 2 points
          XASSERTM(_vtx_ptr.size() > std::size_t(1), "Bezier needs at least 2 points");
          XASSERTM(_param.empty() || (_param.size() == _vtx_ptr.size()), "Bezier world/param size mismatch");

          // check last vertex pointer
          XASSERTM(_vtx_ptr.front() == std::size_t(0), "invalid vertex pointer array");
          XASSERTM(_vtx_ptr.back() == _world.size(), "invalid vertex pointer array");

          // validate inner vertex pointers
          for(std::size_t i(0); (i+1) < _vtx_ptr.size(); ++i)
          {
            XASSERTM(_vtx_ptr[i] < _vtx_ptr[i+1], "invalid vertex pointer array");
            XASSERTM(int(_vtx_ptr[i+1]-_vtx_ptr[i]) <= max_degree, "invalid spline degree");
          }

          XASSERT(_orientation == DataType(1) || _orientation == -DataType(1));
        }

        // use default copy ctor; this one is required by the MeshExtruder !
        Bezier(const Bezier&) = default;

        /** \copydoc ChartBase::get_type() */
        virtual String get_type() const override
        {
          return "bezier";
        }

        /**
         * \brief Pushes a new vertex point to the spline
         *
         * \param[in] point
         * The vertex point to be pushed
         */
        void push_vertex(const WorldPoint& point)
        {
          // update vertex pointer
          _vtx_ptr.push_back(_world.size());
          _world.push_back(point);
        }

        /**
         * \brief Pushes a new control point to the spline
         *
         * \param[in] point
         * The control point to be pushed
         */
        void push_control(const WorldPoint& point)
        {
          _world.push_back(point);
        }

        /**
         * \brief Pushes a new parameter to the spline
         *
         * \param[in] param
         * The parameter to be pushed
         */
        void push_param(const ParamPoint& param)
        {
          this->_param.push_back(param);
        }

        virtual bool can_explicit() const override
        {
          return !_param.empty();
        }

        std::deque<WorldPoint>& get_world_points()
        {
          return _world;
        }

        const std::deque<WorldPoint>& get_world_points() const
        {
          return _world;
        }

        /// \copydoc ChartBase:transform()
        virtual void transform(const WorldPoint& origin, const WorldPoint& angles, const WorldPoint& offset) override
        {
          // create rotation matrix
          Tiny::Matrix<DataType, 2, 2> rot;
          rot.set_rotation_2d(angles(0));

          // transform all world points
          WorldPoint tmp;
          for(auto& pt : _world)
          {
            tmp = pt - origin;
            pt.set_mat_vec_mult(rot, tmp) += offset;
          }
        }

        /**
         * \brief Maps a local segment parameter point
         *
         * \param[in] i
         * The segment index onto which to map
         *
         * \param[in] t
         * The local segment parameter
         */
        WorldPoint map_on_segment(const Index i, const DataType t) const
        {
          ASSERTM((i+1) < Index(this->_vtx_ptr.size()), "invalid segment index");

          const DataType t1 = DataType(1) - t;

          // get number of points on segment
          const std::size_t n = _vtx_ptr[i+1] - _vtx_ptr[i];

          // check for simple line segment (1st degree spline)
          if(n == std::size_t(1))
          {
            // perform linear interpolation
            std::size_t k = _vtx_ptr[i];
            return (t1 * _world[k]) + (t * _world[k+1]);
          }

          // temporary work array
          WorldPoint v[max_degree+1];

          // fetch all points on current segment
          for(std::size_t k(0); k <= n; ++k)
          {
            v[k] = _world[_vtx_ptr[i] + k];
          }

          // apply recursive Bezier interpolation scheme
          for(std::size_t j(0); j < n; ++j)
          {
            for(std::size_t k(0); (k+j) < n; ++k)
            {
              // v_new[k] <- (1-t)*v[k] + t*v[k+1]  <==>
              (v[k] *= t1) += (t * v[k+1]);
            }
          }

          // the interpolated result is stored in v[0]
          return v[0];
        }

        /**
         * \brief Computes the outer unit normal in a single point
         *
         * \param[in] i
         * The segment index onto which to compute the normal.
         *
         * \param[in] t
         * The local segment parameter.
         *
         * \returns The outer unit normal in the point given by t.
         */
        WorldPoint get_normal_on_segment(const Index i, const DataType t) const
        {
          // compute first order derivative
          ASSERTM((i+1) < Index(this->_vtx_ptr.size()), "invalid segment index");

          // our normal vector
          WorldPoint nu;

          // get number of points on segment
          const std::size_t n = _vtx_ptr[i+1] - _vtx_ptr[i];

          // check for simple line segment (1st degree spline)
          if(n == std::size_t(1))
          {
            std::size_t k = _vtx_ptr[i];
            nu[0] = _world[k+1][1] - _world[k][1];
            nu[1] = _world[k][0] - _world[k+1][0];
            return nu;
          }

          // we have a spline of degree > 1

          WorldPoint vals[max_degree+1];
          WorldPoint der1[max_degree+1];

          // fetch all points on current segment
          for(std::size_t k(0); k <= n; ++k)
          {
            vals[k] = _world[_vtx_ptr[i] + k];
            der1[k] = DataType(0);
          }

          const DataType t1 = DataType(1) - t;

          // apply recursive Bezier interpolation scheme
          for(std::size_t j(0); j < n; ++j)
          {
            for(std::size_t k(0); (k+j) < n; ++k)
            {
              // update 1st order derivative
              // v_new'[k] <- (1-t)*v'[k] + t*v'[k+1] - v[k] + v[k+1]
              (((der1[k] *= t1) += (t * der1[k+1])) -= vals[k]) += vals[k+1];

              // update function value
              // v_new[k] <- (1-t)*v[k] + t*v[k+1]
              (vals[k] *= t1) += (t * vals[k+1]);
            }
          }

          // rotate to obtain normal
          nu[0] =  _orientation*der1[0][1];
          nu[1] = -_orientation*der1[0][0];
          return nu;
        }

        /**
         * \brief Projects a point onto one segment of the spline.
         *
         * \note
         * For the description of the Selimovic elimination test used in this function, see \cite Selimovic2006
         *
         * \param[in] i
         * The index of the spline segment onto which to project
         *
         * \param[in] point
         * The world point that is to be projected.
         *
         * \returns
         * The projected world point parameter.
         */
        DataType project_on_segment(const Index i, const WorldPoint& point) const
        {
          ASSERTM((i+1) < Index(this->_vtx_ptr.size()), "invalid segment index");

          // get number of points on segment
          const std::size_t n = _vtx_ptr[i+1] - _vtx_ptr[i];

          // check for simple line segment (1st degree spline)
          if(n == std::size_t(1))
          {
            // get the ends of the line segment
            const std::size_t k = _vtx_ptr[i];
            const WorldPoint& x0 = this->_world.at(k);
            const WorldPoint& x1 = this->_world.at(k+1);

            // compute xe := x1-x0 and xp := p-x0
            WorldPoint xe = x1 - x0;
            WorldPoint xp = point - x0;

            // compute t = <xp,xe>/<xe,xe> and clamp it to [0,1]
            return Math::clamp(Tiny::dot(xp,xe) / Tiny::dot(xe,xe), DataType(0), DataType(1));
          }

          // apply Selimovic elimination test
          {
            // get curve start- and end-points
            const std::size_t ks = _vtx_ptr[i];
            const std::size_t ke = _vtx_ptr[i+1];
            const WorldPoint& xs = this->_world.at(ks);
            const WorldPoint& xe = this->_world.at(ke);
            const WorldPoint ps = xs - point;
            const WorldPoint pe = xe - point;

            // which is closer to our query point: the start or the end-point?
            if(ps.norm_euclid_sqr() < pe.norm_euclid_sqr())
            {
              // test whether the start-point is the projection
              bool prj_xs = true;
              for(std::size_t k(ks + 1); k < ke; ++k)
              {
                prj_xs = prj_xs && (Tiny::dot(this->_world.at(k) - xs, ps) > DataType(0));
              }

              if(prj_xs)
                return DataType(0);
            }
            else
            {
              // test whether the end-point is the projection
              bool prj_xe = true;
              for(std::size_t k(ks + 1); k < ke; ++k)
              {
                prj_xe = prj_xe && (Tiny::dot(this->_world.at(k) - xe, pe) > DataType(0));
              }
              if(prj_xe)
                return DataType(1);
            }
          }

          //
          // Algorithm description
          // ---------------------
          // This function is meant to find the parameter value 't' of the point
          // B(t) closest to the given input point 'X' (named 'point'), i.e.
          //
          //                 argmin     { (B(t) - X)^2 }
          //             {0 <= t <= 1}
          //
          // As Bezier curves a polynomials (and therefore smooth), the parameter
          // 't' we are looking for fulfills the orthogonal projection property
          //
          //              f(t) :=   < B'(t), B(t) - X > = 0
          //
          // The algorithm implemented below is a simple Newton iteration applied
          // onto the non-linear equation above.
          // For Newton, we require the derivative f'(t) of f(t), which is
          //
          //           f'(t) = < B"(t), B(t) - X > + < B'(t), B'(t) >
          //
          // and the Newton iteration is defined as
          //
          //               t_{k+1} := t_k - f(t_k) / f'(t_k)
          //
          // One problem is to find an appropriate initial guess t_0 for our
          // algorithm. In many cases, the initial guess t_0 = 1/2 will do the
          // job. However, if the Bezier curve is S-shaped, this algorithm may
          // fail, therefore we have to choose another initial guess.
          // To circumvent this problem, we explicitly test the distance
          // (B(t) - X)^2 for t in {1/5, 1/2, 4/5} and choose the t with
          // the minimal distance as an initial guess.
          //
          // Unfortunately, there is no guarantee that this algorithm will
          // converge.

          // temporary work arrays
          WorldPoint vals[max_degree+1]; // B (t)
          WorldPoint der1[max_degree+1]; // B'(t)
          WorldPoint der2[max_degree+1]; // B"(t)

          // choose an appropriate initial parameter value
          DataType t = DataType(0.5);
          if(n <= std::size_t(3))
          {
            // choose t in {0.2, 0.5, 0.8}
            DataType d1 = (point - map_on_segment(i, DataType(0.2))).norm_euclid_sqr();
            DataType d2 = (point - map_on_segment(i, DataType(0.5))).norm_euclid_sqr();
            DataType d3 = (point - map_on_segment(i, DataType(0.8))).norm_euclid_sqr();
            if(d1 < Math::min(d2,d3)) t = DataType(0.2);
            if(d2 < Math::min(d3,d1)) t = DataType(0.5);
            if(d3 < Math::min(d1,d2)) t = DataType(0.8);
          }
          else // n > 3
          {
            // choose t in {0.0, 0.25, 0.5, 0.75, 1.0}
            DataType d1 = (point - map_on_segment(i, DataType(0.00))).norm_euclid_sqr();
            DataType d2 = (point - map_on_segment(i, DataType(0.25))).norm_euclid_sqr();
            DataType d3 = (point - map_on_segment(i, DataType(0.50))).norm_euclid_sqr();
            DataType d4 = (point - map_on_segment(i, DataType(0.75))).norm_euclid_sqr();
            DataType d5 = (point - map_on_segment(i, DataType(1.00))).norm_euclid_sqr();
            if(d1 < Math::min(Math::min(d2,d3), Math::min(d4,d5))) t = DataType(0.00);
            if(d2 < Math::min(Math::min(d3,d4), Math::min(d5,d1))) t = DataType(0.25);
            if(d3 < Math::min(Math::min(d4,d5), Math::min(d1,d2))) t = DataType(0.50);
            if(d4 < Math::min(Math::min(d5,d1), Math::min(d2,d3))) t = DataType(0.75);
            if(d5 < Math::min(Math::min(d1,d2), Math::min(d3,d4))) t = DataType(1.00);
          }

          // Newton-Iteration
          for(int iter(0); iter < 10; ++iter)
          {
            // pre-compute (1-t)
            const DataType t1 = DataType(1) - t;

            // fetch all points on current segment
            for(std::size_t k(0); k <= n; ++k)
            {
              vals[k] = _world[_vtx_ptr[i] + k];
              der1[k] = der2[k] = DataType(0);
            }

            // apply recursive Bezier interpolation
            for(std::size_t j(0); j < n; ++j)
            {
              for(std::size_t k(0); (k+j) < n; ++k)
              {
                // update 2nd order derivative
                // v_new"[k] <- (1-t)*v"[k] + t*v"[k+1] - 2*v'[k] + 2*v'[k+1]
                (((der2[k] *= t1) += (t * der2[k+1])) -= DataType(2)*der1[k]) += DataType(2)*der1[k+1];

                // update 1st order derivative
                // v_new'[k] <- (1-t)*v'[k] + t*v'[k+1] - v[k] + v[k+1]
                (((der1[k] *= t1) += (t * der1[k+1])) -= vals[k]) += vals[k+1];

                // update function value
                // v_new[k] <- (1-t)*v[k] + t*v[k+1]
                (vals[k] *= t1) += (t * vals[k+1]);
              }
            }

            // map local point and subtract world point: B(t) - X
            DataType vpx = vals[0][0] - point[0];
            DataType vpy = vals[0][1] - point[1];

            // compute first derivatives: B'(t)
            DataType d1x = der1[0][0];
            DataType d1y = der1[0][1];

            // compute second derivatives: B"(t)
            DataType d2x = der2[0][0];
            DataType d2y = der2[0][1];

            // compute function value: f(t) := < B'(t), B(t) - X >
            DataType fv = (d1x*vpx + d1y*vpy);
            if(Math::abs(fv) < DataType(1E-8))
              return t; // finished!

            // compute function derivative:  f'(t) := < B"(t), B(t) - X > + < B'(t), B'(t) >
            DataType fd = (d2x*vpx + d2y*vpy + d1x*d1x + d1y*d1y);
            if(Math::abs(fd) < DataType(1E-8))
            {
              // This should not happen...
              break;
            }

            // compute t-update:
            DataType tu = -fv / fd;

            // ensure that we do not try to run beyond our segment;
            // in this case the projected point is one of our segment ends
            if((t <= DataType(0)) && (tu <= DataType(0))) break;
            if((t >= DataType(1)) && (tu >= DataType(0))) break;

            // compute new t and clamp to [0,1]
            t = Math::clamp(t + tu, DataType(0), DataType(1));
          }

          // Maximum number of iterations reached...
          return t;
        }

        /**
         * \brief Maps a single parameter point
         *
         * \param[out] point
         * The image of the parameter point under the chart mapping
         *
         * \param[in] param
         * The parameter point to be mapped
         */
        void map_param(WorldPoint& point, const ParamPoint& param) const
        {
          XASSERTM(!this->_param.empty(), "Bezier has no parameters");

          // find enclosing segment
          const Index i = Index(this->find_param(param[0]));

          // compute local segment parameter
          const DataType t = (param[0] - this->_param[i][0]) / (this->_param[i+1][0] - this->_param[i][0]);

          // map on segment
          point = map_on_segment(i, t);
        }

        /**
         * \brief Projects a single world point
         *
         * \param[in,out] point
         * The world point to be projected
         */
        void project_point(WorldPoint& point) const
        {
          // create a const copy of our input point
          const WorldPoint inpoint(point);

          // find the vertex that is closest to our input point and remember the distance
          point =  this->_world.front();
          DataType min_dist = (inpoint - point).norm_euclid_sqr();
          for(std::size_t i(1); i < this->_vtx_ptr.size(); ++i)
          {
            DataType distance = (inpoint - this->_world.at(this->_vtx_ptr[i])).norm_euclid_sqr();
            if(distance < min_dist)
            {
              point = this->_world.at(this->_vtx_ptr[i]);
              min_dist = distance;
            }
          }

          // compute square root to obtain real distance (and not squared one)
          min_dist = Math::sqrt(min_dist);

          // loop over all curve segments
          for(Index i(0); (i+1) < Index(this->_vtx_ptr.size()); ++i)
          {
            // we now check the bounding box of the bezier curve segment against a
            // "test box" of the current minimum distance around our input point

            // can we skip this curve segment?
            bool gt_x0 = false; // not left of our test box
            bool lt_x1 = false; // not right of our test box
            bool gt_y0 = false; // not below our test box
            bool lt_y1 = false; // not above our test box

            // loop over all points on this curve segment
            for(std::size_t k(this->_vtx_ptr[i]); k <= this->_vtx_ptr[i+1]; ++k)
            {
              // get the vertex/control point - input point
              const WorldPoint p = this->_world.at(k) - inpoint;

              // check the position of this point
              gt_x0 = gt_x0 || (p[0] > -min_dist);
              lt_x1 = lt_x1 || (p[0] < +min_dist);
              gt_y0 = gt_y0 || (p[1] > -min_dist);
              lt_y1 = lt_y1 || (p[1] < +min_dist);
            }

            // we can skip the curve if all points of the curve are on one side of our test box:
            if(!(gt_x0 && lt_x1 && gt_y0 && lt_y1)) // <==> (!gt_x0 || !lt_x1 || !gt_y0 || !lt_y1)
              continue;

            // project on current segment
            const DataType t = project_on_segment(i, inpoint);

            // map point on current segment
            const WorldPoint x = map_on_segment(i, t);

            // compute distance to original point
            DataType distance = (x - inpoint).norm_euclid();

            // is that a new projection candidate?
            if(distance < min_dist)
            {
              point = x;
              min_dist = distance;
            }
          }
        }

        /**
         * \brief Projects all mesh points identified by a meshpart
         *
         * \param[in,out] mesh
         * The mesh whose points will be projected
         *
         * \param[in] meshpart
         * The MeshPart identifying the point to be projected
         *
         * \todo implement fancy variant here
         */
        void project_meshpart(Mesh_& mesh, const MeshPart<Mesh_>& meshpart) const
        {
          auto& vtx = mesh.get_vertex_set();
          const auto& target_vtx = meshpart.template get_target_set<0>();

          for(Index i(0); i < meshpart.get_num_entities(0); ++i)
          {
            project_point(reinterpret_cast<WorldPoint&>(vtx[target_vtx[i]]));
          }
        }

        /// \copydoc ChartBase::dist()
        DataType compute_dist(const WorldPoint& point) const
        {
          WorldPoint projected(point);
          project_point(projected);
          return (projected - point).norm_euclid();
        }

        /// \copydoc ChartBase::dist()
        DataType compute_dist(const WorldPoint& point, WorldPoint& grad_dist) const
        {
          WorldPoint projected(point);
          project_point(projected);

          grad_dist = (projected - point);

          return grad_dist.norm_euclid();
        }

        /// \copydoc ChartBase::signed_dist()
        DataType compute_signed_dist(const WorldPoint& point) const
        {
          WorldPoint grad_dist(DataType(0));

          return compute_signed_dist(point, grad_dist);
        }

        /// \copydoc ChartBase::signed_dist()
        DataType compute_signed_dist(const WorldPoint& point, WorldPoint& grad_dist) const
        {
          DataType best_distance_sqr(Math::huge<DataType>());
          DataType best_sign(_orientation*Math::huge<DataType>());

          WorldPoint projected(DataType(0));
          WorldPoint best_nu(DataType(0));

          // loop over all line segments
          for(Index i(0); (i+1) < Index(this->_vtx_ptr.size()); ++i)
          {
            // project on current segment
            const DataType t = project_on_segment(i, point);

            // map point on current segment
            projected = map_on_segment(i, t);

            WorldPoint difference(projected - point);

            // compute squared distance to original point
            DataType my_distance_sqr(difference.norm_euclid_sqr());
            WorldPoint nu(get_normal_on_segment(i, t));

            // If we have a point inside this segment we are happy and can leave
            if(my_distance_sqr <= Math::sqr(Math::eps<DataType>()))
            {
              best_distance_sqr = DataType(0);
              best_nu = nu;
              best_sign = -DataType(1);
              break;
            }

            DataType my_sign(1);
            if(_closed)
            {
              // Compute normal so we can compute the sign
              my_sign = Math::signum(Tiny::dot(nu, difference));

              // If the sign is 0, the projected point lies in the continuation of the current segment, meaning
              // t got clamped to 0 or 1. So we compute the sign according to the normal wrt. to the other segment.
              // If the sign is still 0 then, the projected point lies in the continuation of the next segment,
              // so the same situation will arise in the next iteration anyway.
              if(my_sign == DataType(0))
              {
                if(Math::abs(t) < Math::sqrt(Math::eps<DataType>()))
                {
                  Index other_segment;
                  (i == Index(0)) ? other_segment = Index(this->_vtx_ptr.size()-2) : other_segment = i-Index(1);
                  WorldPoint other_nu(get_normal_on_segment(other_segment, DataType(1)));
                  my_sign = Math::signum(Tiny::dot(difference, other_nu));
                }
                else if(Math::abs(t - DataType(1)) < Math::sqrt(Math::eps<DataType>()))
                {
                  Index other_segment;
                  (i == Index(this->_vtx_ptr.size()-2)) ? other_segment = Index(0) : other_segment = i+Index(1);
                  WorldPoint other_nu(get_normal_on_segment(other_segment, DataType(0)));
                  my_sign = Math::signum(Tiny::dot(difference, other_nu));
                }
                //else
                //  throw InternalError(__func__,__FILE__,__LINE__,
                //  "signum == "+stringify_fp_sci(my_sign)+" but t = "+stringify_fp_sci(t));
              }
            }

            // Update the projection candidate iff it is the first segment, or the distance is lower
            if(i == Index(0) || (my_distance_sqr < best_distance_sqr))
            {
              best_distance_sqr = my_distance_sqr;
              best_nu = nu;
              // best_sign can be 0 if this happens for i==0
              best_sign = my_sign;
              grad_dist = difference;
            }
          }

          // If the point was far enough away from the interface, we can normalise the difference vector
          if(best_distance_sqr > Math::sqr(Math::eps<DataType>()) && best_sign != DataType(0))
          {
            grad_dist.normalise();
            grad_dist *= -best_sign;
          }
          // If we are too close, we take the normal in the projected point. We do not ALWAYS do this because the
          // projected point might be where the normal is discontinuous, but in this case it is ok.
          else
          {
            grad_dist = DataType(-1)*best_nu;
            grad_dist.normalise();
            XASSERT(Math::abs(grad_dist.norm_euclid()-DataType(1)) < Math::sqrt(Math::eps<DataType>()));
          }

          return best_sign*Math::sqrt(best_distance_sqr);
        }

        /** \copydoc ChartBase::write */
        virtual void write(std::ostream& os, const String& sindent) const override
        {
          String sind(sindent), sind2(sindent);
          if(!sind.empty())
          {
            sind.append("  ");
            sind2.append("    ");
          }

          os << sindent << "<Bezier dim=\"2\" size=\"" << this->_vtx_ptr.size() << "\"";
          os << " type=\"" << (this->_closed ? "closed" : "open") << "\"";
          os << (_orientation == -DataType(1) ? " orientation=\"-1\"" : "" )<<">" << std::endl;

          // write points
          os << sind << "<Points>" << std::endl;
          // write first vertex point
          os << sind2 << 0;
            for(int j(0); j < BaseClass::world_dim; ++j)
              os << " " << _world[_vtx_ptr.front()][j];
          // write remaining points
          for(std::size_t i(1); i < _vtx_ptr.size(); ++i)
          {
            // write number of control points
            os << sind2 << ( _vtx_ptr[i] - _vtx_ptr[i-1] - std::size_t(1));
            // write point coordinates
            for(std::size_t k(_vtx_ptr[i-1]+1); k <= _vtx_ptr[i]; ++k)
              for(int j(0); j < BaseClass::world_dim; ++j)
                os << " " << _world[k][j];
            os << std::endl;
          }
          os << sind << "</Points>" << std::endl;

          // write parameters
          if(!_param.empty())
          {
            os << sind << "<Params>" << std::endl;
            for(std::size_t i(0); i < _param.size(); ++i)
              os << sind2 << _param[i][0] << std::endl;;
            os << sind << "</Params>" << std::endl;
          }

          // write terminator
          os << sindent << "</Bezier>" << std::endl;
        }

      protected:
        std::size_t find_param(const DataType x) const
        {
          XASSERTM(!_param.empty(),"Bezier has no parameters");

          // check for boundary
          if(x <= _param.front()[0])
            return std::size_t(0);
          if(x >= _param.back()[0])
            return _param.size()-std::size_t(2);

          // apply binary search
          std::size_t il(0), ir(_param.size()-1);
          while(il+1 < ir)
          {
            // test median
            std::size_t im = (il+ir)/2;
            DataType xm = _param.at(im)[0];
            if(x < xm)
            {
              ir = im;
            }
            else
            {
              il = im;
            }
          }

          // return interval index
          return (il+1 < _param.size() ? il : (_param.size() - std::size_t(2)));
        }
      }; // class Bezier<...>

      template<typename Mesh_>
      class BezierPointsParser :
        public Xml::MarkupParser
      {
        typedef Bezier<Mesh_> ChartType;

      private:
        Bezier<Mesh_>& _bezier;
        Index _size, _read;

      public:
        explicit BezierPointsParser(Bezier<Mesh_>& bezier, Index size) :
          _bezier(bezier), _size(size), _read(0)
        {
        }

        virtual bool attribs(std::map<String,bool>&) const override
        {
          return true;
        }

        virtual void create(int iline, const String& sline, const String&, const std::map<String, String>&, bool closed) override
        {
          if(closed)
            throw Xml::GrammarError(iline, sline, "Invalid closed markup");
        }

        virtual void close(int iline, const String& sline) override
        {
          // ensure that we have read all vertices
          if(_read < _size)
            throw Xml::GrammarError(iline, sline, "Invalid terminator; expected point");
        }

        virtual std::shared_ptr<MarkupParser> markup(int, const String&, const String&) override
        {
          // no children allowed
          return nullptr;
        }

        virtual bool content(int iline, const String& sline) override
        {
          // make sure that we do not read more points than expected
          if(_read >= _size)
            throw Xml::ContentError(iline, sline, "Invalid content; exprected terminator");

          // split line by whitespaces
          std::deque<String> scoords = sline.split_by_whitespaces();

          // try to parse the control point count
          std::size_t num_ctrl(0);
          if(!scoords.front().parse(num_ctrl))
            throw Xml::ContentError(iline, sline, "Failed to parse control point count");

          // the first point must not have control points
          if((_read == Index(0)) && (num_ctrl > std::size_t(0)))
            throw Xml::ContentError(iline, sline, "First point must be a vertex point");

          // check size: 1 + (#ctrl+1) * #coords
          if(scoords.size() != ((num_ctrl+1) * std::size_t(ChartType::world_dim) + std::size_t(1)))
            throw Xml::ContentError(iline, sline, "Invalid number of coordinates");

          typename ChartType::WorldPoint point;

          // try to parse the control points
          for(int k(0); k < int(num_ctrl); ++k)
          {
            for(int i(0); i < ChartType::world_dim; ++i)
            {
              if(!scoords.at(std::size_t(k*ChartType::world_dim+i+1)).parse(point[i]))
                throw Xml::ContentError(iline, sline, "Failed to parse control point coordinate");
            }
            // push the control point
            _bezier.push_control(point);
          }

          // try to parse all coords of the vertex point
          for(int i(0); i < ChartType::world_dim; ++i)
          {
            if(!scoords.at(std::size_t(int(num_ctrl)*ChartType::world_dim+i+1)).parse(point[i]))
              throw Xml::ContentError(iline, sline, "Failed to parse vertex point coordinate");
          }

          // push the vertex point
          _bezier.push_vertex(point);

          // okay, another point done
          ++_read;

          return true;
        }
      }; // BezierParamsParser

      template<typename Mesh_>
      class BezierParamsParser :
        public Xml::MarkupParser
      {
        typedef Bezier<Mesh_> ChartType;

      private:
        Bezier<Mesh_>& _bezier;
        Index _size, _read;

      public:
        explicit BezierParamsParser(Bezier<Mesh_>& bezier, Index size) :
          _bezier(bezier), _size(size), _read(0)
        {
        }

        virtual bool attribs(std::map<String,bool>&) const override
        {
          return true;
        }

        virtual void create(int iline, const String& sline, const String&, const std::map<String, String>&, bool closed) override
        {
          if(closed)
            throw Xml::GrammarError(iline, sline, "Invalid closed markup");
        }

        virtual void close(int iline, const String& sline) override
        {
          // ensure that we have read all vertices
          if(_read < _size)
            throw Xml::GrammarError(iline, sline, "Invalid terminator; expected point");
        }

        virtual std::shared_ptr<MarkupParser> markup(int, const String&, const String&) override
        {
          // no children allowed
          return nullptr;
        }

        virtual bool content(int iline, const String& sline) override
        {
          // make sure that we do not read more points than expected
          if(_read >= _size)
            throw Xml::ContentError(iline, sline, "Invalid content; exprected terminator");

          typename ChartType::ParamPoint point;

          // try to parse all coords
          if(!sline.parse(point[0]))
            throw Xml::ContentError(iline, sline, "Failed to parse param coordinate");

          // push
          _bezier.push_param(point);

          // okay, another point done
          ++_read;

          return true;
        }
      }; // class BezierParamsParser

      template<typename Mesh_, typename ChartReturn_ = ChartBase<Mesh_>>
      class BezierChartParser :
        public Xml::MarkupParser
      {
      private:
        typedef Bezier<Mesh_> ChartType;
        typedef typename ChartType::DataType DataType;
        ChartReturn_*& _chart;
        Bezier<Mesh_>* _bezier;
        Index _size;

      public:
        explicit BezierChartParser(ChartReturn_*& chart) :
          _chart(chart),
          _bezier(nullptr),
          _size(0)
        {
        }

        virtual bool attribs(std::map<String,bool>& attrs) const override
        {
          attrs.emplace("dim", true);
          attrs.emplace("size", true);
          attrs.emplace("type", false);
          attrs.emplace("orientation", false);
          return true;
        }

        virtual void create(
          int iline,
          const String& sline,
          const String&,
          const std::map<String, String>& attrs,
          bool closed) override
        {
          // make sure this one isn't closed
          if(closed)
            throw Xml::GrammarError(iline, sline, "Invalid closed Bezier markup");

          Index dim(0);

          // try to parse the dimension
          if(!attrs.find("dim")->second.parse(dim))
            throw Xml::GrammarError(iline, sline, "Failed to parse Bezier dimension");
          if(dim != Index(2))
            throw Xml::GrammarError(iline, sline, "Invalid Bezier dimension");

          // try to parse the size
          if(!attrs.find("size")->second.parse(_size))
            throw Xml::GrammarError(iline, sline, "Failed to parse Bezier size");
          if(_size < Index(2))
            throw Xml::GrammarError(iline, sline, "Invalid Bezier size");

          // try to check type
          bool poly_closed(false);
          auto it = attrs.find("type");
          if(it != attrs.end())
          {
            String stype = it->second;
            if(it->second == "closed")
              poly_closed = true;
            else if (it->second != "open")
              throw Xml::ContentError(iline, sline, "Invalid Bezier type; must be either 'closed' or 'open'");
          }

          DataType orientation(1);
          it = attrs.find("orientation");
          if(it != attrs.end())
            it->second.parse(orientation);

          // up to now, everything's fine
          _bezier = new ChartType(poly_closed, orientation);
        }

        virtual void close(int, const String&) override
        {
          // okay
          _chart = _bezier;
        }

        virtual bool content(int, const String&) override
        {
          return false;
        }

        virtual std::shared_ptr<Xml::MarkupParser> markup(int, const String&, const String& name) override
        {
          if(name == "Points") return std::make_shared<BezierPointsParser<Mesh_>>(*_bezier, _size);
          if(name == "Params") return std::make_shared<BezierParamsParser<Mesh_>>(*_bezier, _size);
          return nullptr;
        }
      }; // class BezierChartParser<...>
    } // namespace Atlas
  } // namespace Geometry
} // namespace FEAT
#endif // KERNEL_GEOMETRY_ATLAS_BEZIER_HPP
