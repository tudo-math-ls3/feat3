#pragma once
#ifndef KERNEL_GEOMETRY_ATLAS_SPLINE_BASE_HPP
#define KERNEL_GEOMETRY_ATLAS_SPLINE_BASE_HPP 1

#include <kernel/geometry/atlas/chart.hpp>

#include <deque>

namespace FEAT
{
  namespace Geometry
  {
    namespace Atlas
    {
      /// Spline chart traits
      struct SplineBaseTraits
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

      template<typename Derived_, typename Mesh_, typename Traits_ = SplineBaseTraits>
      class SplineBaseCRTP :
        public ChartCRTP<Derived_, Mesh_, Traits_>
      {
      public:
        /// The CRTP base class
        typedef ChartCRTP<Derived_, Mesh_, Traits_> BaseClass;
        /// Floating point type for coordinates
        typedef typename BaseClass::CoordType DataType;
        /// Vector type for world points aka. image points
        typedef typename BaseClass::WorldPoint WorldPoint;
        /// Vector type for parameter points aka. domain points
        typedef typename BaseClass::ParamPoint ParamPoint;

      protected:
        /// The world points for the spline
        std::deque<WorldPoint> _world;
        /// The parameter values for these world points
        std::deque<ParamPoint> _param;
        /// Specifies whether the spline is closed
        bool _closed;

        Derived_& cast() {return static_cast<Derived_&>(*this);}
        const Derived_& cast() const {return static_cast<const Derived_&>(*this);}

      public:
        /// default CTOR
        explicit SplineBaseCRTP(bool closed = false) :
          _world(),
          _param(),
          _closed(closed)
        {
        }

        void push_world(const WorldPoint& world)
        {
          _world.push_back(world);
        }

        void push_param(const ParamPoint& param)
        {
          _param.push_back(param);
        }

        bool can_explicit() const
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

        /// \copydoc ChartBase::move_by()
        virtual void move_by(const WorldPoint& translation) override
        {
          for(auto& it:_world)
            it+=translation;
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
        void map(WorldPoint& point, const ParamPoint& param) const
        {
          XASSERTM(!this->_param.empty(), "Parameter is missing");

          // find enclosing segment
          const Index i = Index(this->find_param(param[0]));

          // compute local segment parameter
          const DataType t = (param[0] - this->_param[i][0]) / (this->_param[i+1][0] - this->_param[i][0]);

          // map on segment
          point = cast().map_on_segment(i, t);
        }

        WorldPoint get_normal_on_segment(const Index i, const DataType t) const
        {
          return cast().get_normal_on_segment(i, t);
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
          DataType mindist(0.0);

          // loop over all line segments
          for(Index i(0); (i+1) < Index(this->_world.size()); ++i)
          {
            // project on current segment
            const DataType t = cast().project_on_segment(i, inpoint);

            // map point on current segment
            const WorldPoint x = cast().map_on_segment(i, t);

            // compute squared distance to original point
            DataType distance = (x - inpoint).norm_euclid_sqr();

            // is that a new projection candidate?
            if((i == Index(0)) || (distance < mindist))
            {
              point = x;
              mindist = distance;
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
            cast().project_point(reinterpret_cast<WorldPoint&>(vtx[target_vtx[i]]));
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
          DataType signed_distance(0.0);

          Index best_segment(0);
          WorldPoint projected(DataType(0));

          // loop over all line segments
          for(Index i(0); (i+1) < Index(this->_world.size()); ++i)
          {
            // project on current segment
            const DataType t = cast().project_on_segment(i, point);

            // map point on current segment
            projected = cast().map_on_segment(i, t);

            // compute squared distance to original point
            DataType my_distance = (projected - point).norm_euclid_sqr();

            // is that a new projection candidate?
            if((i == Index(0)) || (my_distance < signed_distance))
            {
              signed_distance = my_distance;
              best_segment = i;
            }
          }

          // Project on best segment
          const DataType t = cast().project_on_segment(best_segment, point);
          // Map point on best segment
          projected = cast().map_on_segment(best_segment, t);

          grad_dist = (projected - point);
          // This has no sign yet
          signed_distance = grad_dist.norm_euclid();

          // If the distance is too small, we set the gradient vector to zero
          if(signed_distance <= Math::eps<DataType>())
            grad_dist.format(DataType(0));
          else
          {
            grad_dist.normalise();
            WorldPoint nu(get_normal_on_segment(best_segment, t));
            signed_distance *= Math::signum(Tiny::dot(nu, grad_dist));
          }

          grad_dist *= Math::signum(signed_distance);

          return signed_distance;
        }


      protected:
        std::size_t find_param(const DataType x) const
        {
          XASSERTM(!_param.empty(),"Spline has no parameters");

          // check for boundary
          if(x <= _param.front()[0])
            return std::size_t(0);
          if(x >= _param.back()[0])
            return _param.size()-std::size_t(2);

          // apply binary search
          DataType xl = _param.front()[0];
          DataType xr = _param.back()[0];
          std::size_t il(0), ir(_param.size()-1);
          while(il+1 < ir)
          {
            // test median
            std::size_t im = (il+ir)/2;
            DataType xm = _param.at(im)[0];
            if(x < xm)
            {
              xr = xm;
              ir = im;
            }
            else
            {
              xl = xm;
              il = im;
            }
          }

          // return interval index
          return (il+1 < _param.size() ? il : (_param.size() - std::size_t(2)));
        }

        void write_points(std::ostream& os, const String& sindent) const
        {
          String sind(sindent);
          if(!sind.empty())
            sind.append("  ");

          os << sindent << "<Points>" << std::endl;
          for(std::size_t i(0); i < _world.size(); ++i)
          {
            os << sind << _world[i][0];
            for(int j(1); j < BaseClass::world_dim; ++j)
              os << " " << _world[i][j];
            os << std::endl;
          }
          os << sindent << "</Points>" << std::endl;
        }

        void write_params(std::ostream& os, const String& sindent) const
        {
          String sind(sindent);
          if(!sind.empty())
            sind.append("  ");

          if(!_param.empty())
          {
            os << sindent << "<Params>" << std::endl;
            for(std::size_t i(0); i < _param.size(); ++i)
              os << sind << _param[i][0] << std::endl;;
            os << sindent << "</Params>" << std::endl;
          }
        }
      }; // class SplineBaseCRTP<...>


      template<typename Chart_>
      class SplineBasePointsParser :
        public Xml::MarkupParser
      {
        typedef Chart_ ChartType;

      private:
        ChartType& _spline;
        Index _size, _read;

      public:
        explicit SplineBasePointsParser(ChartType& spline, Index size) :
          _spline(spline), _size(size), _read(0)
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
          std::deque<String> scoords;
          sline.split_by_charset(scoords);

          // check size
          if(scoords.size() != std::size_t(ChartType::world_dim))
            throw Xml::ContentError(iline, sline, "Invalid number of coordinates");

          typename ChartType::WorldPoint point;

          // try to parse all coords
          for(int i(0); i < ChartType::world_dim; ++i)
          {
            if(!scoords.at(std::size_t(i)).parse(point[i]))
              throw Xml::ContentError(iline, sline, "Failed to parse point coordinate");
          }

          // push
          _spline.push_world(point);

          // okay, another point done
          ++_read;

          return true;
        }
      }; // SplineBasePointsParser


      template<typename Chart_>
      class SplineBaseParamsParser :
        public Xml::MarkupParser
      {
        typedef Chart_ ChartType;

      private:
        ChartType& _spline;
        Index _size, _read;

      public:
        explicit SplineBaseParamsParser(ChartType& spline, Index size) :
          _spline(spline), _size(size), _read(0)
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
          _spline.push_param(point);

          // okay, another point done
          ++_read;

          return true;
        }
      }; // class SplineBaseParamsParser

    } // namespace Atlas
  } // namespace Geometry
} // namespace FEAT
#endif // KERNEL_GEOMETRY_ATLAS_SPLINE_BASE_HPP
