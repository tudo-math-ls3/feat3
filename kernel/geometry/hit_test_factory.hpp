// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/geometry/target_set.hpp>
#include <kernel/geometry/index_set.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/atlas/chart.hpp>
#include <kernel/util/tiny_algebra.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      template<typename VertexSet_, typename IndexTuple_>
      inline Tiny::Vector<typename VertexSet_::CoordType, VertexSet_::num_coords>
      get_midpoint(const VertexSet_& vertex_set, const IndexTuple_& index_tuple)
      {
        using DataType = typename VertexSet_::CoordType;
        constexpr int num_vertex = IndexTuple_::num_indices;

        Tiny::Vector<DataType, VertexSet_::num_coords> mid_point(DataType(0));
        for(int i(0); i < num_vertex; ++i)
        {
          for(int j(0); j < VertexSet_::num_coords; ++j)
          {
            mid_point[j] += vertex_set[index_tuple[i]][j];
          }
        }
        return mid_point * (DataType(1) / DataType(num_vertex));
      }

      template<typename VertexSet_>
      static Tiny::Vector<typename VertexSet_::CoordType, VertexSet_::num_coords>
      get_midpoint(const VertexSet_& vertex_set, const Index vertex)
      {
        using DataType = typename VertexSet_::CoordType;

        Tiny::Vector<DataType, VertexSet_::num_coords> mid_point;
        for(int j(0); j < VertexSet_::num_coords; ++j)
        {
          mid_point[j] = vertex_set[vertex][j];
        }
        return mid_point;
      }

      template<typename HitFunc_, typename Mesh_, int dim_ = Mesh_::shape_dim>
      inline void hittest_compute_target_data(
        std::vector<std::vector<Index>>& target_data,
        const Mesh_& mesh,
        const HitFunc_& hit_test)
      {
        if constexpr (dim_ == 0)
        {
          const Index num_cells(mesh.get_num_entities(0));

          target_data[0].reserve(num_cells);
          for(Index i(0); i < num_cells; ++i)
          {
            if (hit_test(get_midpoint(mesh.get_vertex_set(), i)))
            {
              target_data[0].push_back(i);
            }
          }
        }
        else
        {
          // Recurse
          hittest_compute_target_data<HitFunc_, Mesh_, dim_ - 1>(target_data, mesh, hit_test);

          const Index num_cells(mesh.get_num_entities(dim_));
          const auto& index_set(mesh.template get_index_set<dim_,0>());

          target_data[dim_].reserve(num_cells);
          for(Index i(0); i < num_cells; ++i)
          {
            if (hit_test(get_midpoint(mesh.get_vertex_set(), index_set[i])))
            {
              target_data[dim_].push_back(i);
            }
          }
        }
      }

      template<typename HitFunc_, typename Mesh_, int dim_ = Mesh_::shape_dim>
      inline void hittest_compute_filtered_target_data(
        std::vector<std::vector<Index>>& target_data,
        const Mesh_& mesh,
        const MeshPart<Mesh_>& part,
        const HitFunc_& hit_test)
      {
        if constexpr (dim_ == 0)
        {
          const Index num_cells(part.get_num_entities(0));
          const TargetSet& vertex_target_set = part.template get_target_set<0>();

          target_data[0].reserve(num_cells);
          for(Index i(0); i < num_cells; ++i)
          {
            const Index vertex = vertex_target_set[i];
            if (hit_test(get_midpoint(mesh.get_vertex_set(), vertex)))
            {
              target_data[0].push_back(vertex);
            }
          }
        }
        else
        {
          // Recurse
          hittest_compute_filtered_target_data<HitFunc_, Mesh_, dim_ - 1>(target_data, mesh, part, hit_test);

          const Index num_cells(part.get_num_entities(dim_));
          const TargetSet& filter_target_set = part.template get_target_set<dim_>();
          const auto& index_set(mesh.template get_index_set<dim_,0>());

          target_data[dim_].reserve(num_cells);
          for(Index i(0); i < num_cells; ++i)
          {
            const Index entity = filter_target_set[i];
            if (hit_test(get_midpoint(mesh.get_vertex_set(), index_set[entity])))
            {
              target_data[dim_].push_back(entity);
            }
          }
        }
      }

      template<typename Shape_, int dim_ = Shape_::dimension>
      inline void write_to_target_set(TargetSetHolder<Shape_>& tsh, std::vector<std::vector<Index>>& source)
      {
        if constexpr (dim_ < 0)
        {
          return;
        }
        else
        {
          // Recurse
          write_to_target_set<Shape_, dim_ - 1>(tsh, source);

          const auto num_cells(Index(source[dim_].size()));
          TargetSet& trg = tsh.template get_target_set<dim_>();
          for(Index i(0); i < num_cells; ++i)
          {
            trg[i] = source[dim_][i];
          }

        }
      }
    } // namespace Intern
    /// \endcond

    /**
     * \brief Hit-Test Factory class template
     *
     * This class template can be used to create a MeshPart for a particular mesh,
     * which consists of all entities that are inside the region characterized by a
     * hit-test function.
     *
     * \tparam HitFunc_
     * A class implementing the HitTest-Function interface. See SphereHitTestFunction
     * for an example.
     *
     * \tparam Mesh_
     * The type of the mesh for which the cell sub-set is to be computed.
     *
     * \author Peter Zajac, Stefan Wahlers
     */
    template<
      typename HitFunc_,
      typename Mesh_>
    class HitTestFactory :
      public Factory< MeshPart<Mesh_> >
    {
    public:
      /// The shape type of the mesh
      typedef typename Mesh_::ShapeType ShapeType;
      /// mesh part type
      typedef MeshPart<Mesh_> MeshType;
      /// target set holder type
      typedef typename MeshType::TargetSetHolderType TargetSetHolderType;

    protected:
      /// reference to the hit-test function
      const HitFunc_& _hit_func;
      /// reference to the input mesh
      const Mesh_& _mesh;
      /// internal data storing the indices
      std::vector<std::vector<Index>> _target_data;

    public:
      /**
       * \brief Creates the factory.
       *
       * \param[in] hit_func
       * A \resident reference to the hit-test function characterizing the region.
       *
       * \param[in] mesh
       * A \resident reference to the mesh for which the cell sub-set is to be computed.
       */
      explicit HitTestFactory(const HitFunc_& hit_func, const Mesh_& mesh) :
        _hit_func(hit_func),
        _mesh(mesh),
        _target_data(std::size_t(_mesh.shape_dim+1))
      {
        Intern::hittest_compute_target_data<HitFunc_, Mesh_>(_target_data, _mesh, _hit_func);
      }

      /**
       * \brief Creates the factory.
       *
       * \param[in] hit_func
       * A \resident reference to the hit-test function characterizing the region.
       *
       * \param[in] mesh
       * A \resident reference to the mesh for which the cell sub-set is to be computed.
       *
       * \param[in] filter
       * A \resident reference to a filter meshpart. Only mesh entities of the filter
       * meshpart are tested. All other entities are considered as outside the region.
       */
      explicit HitTestFactory(const HitFunc_& hit_func, const Mesh_& mesh, const MeshType& filter) :
        _hit_func(hit_func),
        _mesh(mesh),
        _target_data(std::size_t(_mesh.shape_dim+1))
      {
        Intern::hittest_compute_filtered_target_data<HitFunc_, Mesh_>(_target_data, _mesh, filter, _hit_func);
      }

      /// \copydoc Factory::get_num_entities()
      virtual Index get_num_entities(int dim) override
      {
        return Index(_target_data.at(std::size_t(dim)).size());
      }

      virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
      {
        Intern::write_to_target_set(target_set_holder, _target_data);
      }

      virtual void fill_attribute_sets(typename MeshType::AttributeSetContainer&) override
      {
        // do nothing as the object has no attribute sets
      }

      virtual void fill_index_sets(std::unique_ptr<typename MeshType::IndexSetHolderType>&) override
      {
        // do nothing as the object has no index sets
      }

    }; // class HitTestFactory

    /**
     * \brief Creates a new mesh-part from a hit-test function
     *
     * \param[in] mesh
     * A \transient reference to the mesh for which a mesh-part is to be created.
     *
     * \param[in] hit_func
     * An object implementing the hit-test function interface.
     *
     * \returns
     * A mesh-part containing all entities for which the hit-test function evaluated to \c true.
     */
    template<typename Mesh_, typename HitFunc_>
    MeshPart<Mesh_> make_meshpart_by_hit_test(const Mesh_& mesh, const HitFunc_& hit_func)
    {
      HitTestFactory<HitFunc_, Mesh_> factory(mesh, hit_func);
      return factory.make();
    }

    /**
     * \brief Creates a new mesh-part from a hit-test function
     *
     * \param[in] mesh
     * A \transient reference to the mesh for which a mesh-part is to be created.
     *
     * \param[in] filter
     * A \transient reference to a filter meshpart. Only mesh entities of the filter
     * meshpart are tested. All other entities are considered as outside the region.
     *
     * \param[in] hit_func
     * An object implementing the hit-test function interface.
     *
     * \returns
     * A mesh-part containing all entities for which the hit-test function evaluated to \c true.
     */
    template<typename Mesh_, typename HitFunc_>
    MeshPart<Mesh_> make_meshpart_by_filtered_hit_test(const Mesh_& mesh, const MeshPart<Mesh_>& filter, const HitFunc_& hit_func)
    {
      HitTestFactory<HitFunc_, Mesh_> factory(mesh, filter, hit_func);
      return factory.make();
    }

    /**
     * \brief Creates a new mesh-part from a hit-test function
     *
     * \param[in] mesh
     * A \transient reference to the mesh for which a mesh-part is to be created.
     *
     * \param[in] hit_func
     * An object implementing the hit-test function interface.
     *
     * \returns
     * A unique-pointer to a mesh-part containing all entities for which the hit-test function evaluated to \c true.
     */
    template<typename Mesh_, typename HitFunc_>
    std::unique_ptr<MeshPart<Mesh_>> make_unique_meshpart_by_hit_test(const Mesh_& mesh, const HitFunc_& hit_func)
    {
      HitTestFactory<HitFunc_, Mesh_> factory(mesh, hit_func);
      return factory.make_unique();
    }

    /**
     * \brief Creates a new mesh-part from a hit-test function
     *
     * \param[in] mesh
     * A \transient reference to the mesh for which a mesh-part is to be created.
     *
     * \param[in] filter
     * A \transient reference to a filter meshpart. Only mesh entities of the filter
     * meshpart are tested. All other entities are considered as outside the region.
     *
     * \param[in] hit_func
     * An object implementing the hit-test function interface.
     *
     * \returns
     * A unique-pointer to a mesh-part containing all entities for which the hit-test function evaluated to \c true.
     */
    template<typename Mesh_, typename HitFunc_>
    std::unique_ptr<MeshPart<Mesh_>> make_unique_meshpart_by_filtered_hit_test(const Mesh_& mesh, const MeshPart<Mesh_>& filter, const HitFunc_& hit_func)
    {
      HitTestFactory<HitFunc_, Mesh_> factory(mesh, filter, hit_func);
      return factory.make_unique();
    }

    /**
     * \brief Sphere hit-test function
     *
     * This function implements a spherical hit-test, i.e. it checks whether a point
     * is inside or outside of a sphere.
     *
     * \tparam DataType_
     * The data-type to be used by the function.
     *
     * \tparam dim_
     * The dimension of the sphere.
     *
     * \author Peter Zajac
     */
    template<typename DataType_, int dim_>
    class SphereHitTestFunction
    {
    public:
      /// The point type
      typedef Tiny::Vector<DataType_, dim_> PointType;

    private:
      /// the sphere's midpoint
      PointType _midpoint;
      /// the sphere's radius
      DataType_ _radius;

    public:
      /**
       * \brief Creates a sphere around the origin.
       *
       * \param[in] radius
       * The desired radius of the sphere. Must be > 0.
       */
      explicit SphereHitTestFunction(DataType_ radius) :
        _midpoint(DataType_(0)),
        _radius(radius)
      {
        XASSERTM(radius > DataType_(0), "sphere radius must be > 0");
      }

      /**
       * \brief Creates a sphere.
       *
       * \param[in] midpoint
       * The midpoint of the sphere.
       *
       * \param[in] radius
       * The desired radius of the sphere. Must be > 0.
       */
      SphereHitTestFunction(PointType midpoint, DataType_ radius) :
        _midpoint(midpoint),
        _radius(radius)
      {
        XASSERTM(radius > DataType_(0), "sphere radius must be > 0");
      }

      /**
       * \brief Performs the hit-test.
       *
       * \param[in] point
       * The point to be tested.
       *
       * \returns
       * \c true, if \p point is inside the sphere, otherwise \c false.
       */
      bool operator()(PointType point) const
      {
        return (point - _midpoint).norm_euclid() <= _radius;
      }
    }; // class SphereHitTestFunction<...>


    /**
     * \brief Chart-based Hit-Test Factory class template
     *
     * This class template can be used to create a MeshPart for a particular mesh, which consists of
     * all entities that are inside or outside the region characterized by a given chart.
     *
     * \tparam Mesh_
     * The type of the mesh for which the cell sub-set is to be computed.
     *
     * \tparam Chart_
     * The type of the chart that is to be tested; typically Atlas::ChartBase.
     *
     * \author Peter Zajac, Stefan Wahlers
     */
    template<typename Mesh_, typename Chart_ = Atlas::ChartBase<Mesh_>>
    class ChartHitTestFactory :
      public Factory< MeshPart<Mesh_> >
    {
    public:
      /// The shape type of the mesh
      typedef typename Mesh_::ShapeType ShapeType;
      /// mesh part type
      typedef MeshPart<Mesh_> MeshType;
      /// target set holder type
      typedef typename MeshType::TargetSetHolderType TargetSetHolderType;

      class ChartHitFunction
      {
        typedef typename Mesh_::CoordType CoordType;
        typedef typename Mesh_::VertexType WorldPointType;

      protected:
        const Geometry::Atlas::ChartBase<Mesh_>& _chart;
        const CoordType _inv;

      public:
        explicit ChartHitFunction(const Geometry::Atlas::ChartBase<Mesh_>& chart, bool invert) :
          _chart(chart),
          _inv(CoordType(invert ? -1 : 1))
        {
        }

        bool operator()(const WorldPointType& point) const
        {
          return _inv * _chart.signed_dist(point) <= CoordType(0);
        }
      };

    protected:
      /// reference to the input mesh
      const Mesh_& _mesh;
      /// our hit function
      ChartHitFunction _hit_func;
      /// internal data storing the indices
      std::vector<std::vector<Index>> _target_data;

    public:
      /**
       * \brief Creates the factory.
       *
       * \param[in] mesh
       * A \resident reference to the mesh for which the cell sub-set is to be computed.
       *
       * \param[in] chart
       * A \resident reference to the chart for which the mesh part is to be created.
       *
       * \param[in] invert
       * If \c true, then the mesh part will contain all mesh entities, which are outside of the chart,
       * i.e. have a signed distance >= 0, otherwise it will contain all entities which are inside of
       * the chart, i.e. have a signed distance <= 0.
       */
      explicit ChartHitTestFactory(const Mesh_& mesh, const Chart_& chart, bool invert = false) :
        _mesh(mesh),
        _hit_func(chart, invert),
        _target_data(std::size_t(_mesh.shape_dim+1))
      {
        Intern::hittest_compute_target_data<ChartHitFunction, Mesh_>(_target_data, _mesh, _hit_func);
      }

      /**
       * \brief Creates the factory.
       *
       * \param[in] mesh
       * A \resident reference to the mesh for which the cell sub-set is to be computed.
       *
       * \param[in] filter
       * A \resident reference to a filter meshpart. Only mesh entities of the filter
       * meshpart are tested. All other entities are considered as outside the region.
       *
       * \param[in] chart
       * A \resident reference to the chart for which the mesh part is to be created.
       *
       * \param[in] invert
       * If \c true, then the mesh part will contain all mesh entities, which are outside of the chart,
       * i.e. have a signed distance >= 0, otherwise it will contain all entities which are inside of
       * the chart, i.e. have a signed distance <= 0.
       */
      explicit ChartHitTestFactory(const Mesh_& mesh, const MeshType& filter, const Chart_& chart, bool invert = false) :
        _mesh(mesh),
        _hit_func(chart, invert),
        _target_data(std::size_t(_mesh.shape_dim+1))
      {
        Intern::hittest_compute_filtered_target_data<ChartHitFunction, Mesh_>(_target_data, _mesh, filter, _hit_func);
      }

      /// \copydoc Factory::get_num_entities()
      virtual Index get_num_entities(int dim) override
      {
        return Index(_target_data.at(std::size_t(dim)).size());
      }

      virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
      {
        Intern::write_to_target_set<ShapeType>(target_set_holder, _target_data);
      }

      virtual void fill_attribute_sets(typename MeshType::AttributeSetContainer&) override
      {
        // do nothing as the object has no attribute sets
      }

      virtual void fill_index_sets(std::unique_ptr<typename MeshType::IndexSetHolderType>&) override
      {
        // do nothing as the object has no index sets
      }
    }; // class ChartHitTestFactory

    /**
     * \brief Creates a new mesh-part from a chart hit-test function
     *
     * \param[in] mesh
     * A \transient reference to the mesh for which a mesh-part is to be created.
     *
     * \param[in] chart
     * A chart object whose signed-distance function is to be used for hit-testing.
     *
     * \returns
     * A mesh-part containing all entities for which the chart signed distance is <= 0.
     */
    template<typename Mesh_, typename Chart_>
    MeshPart<Mesh_> make_meshpart_by_chart_hit_test(const Mesh_& mesh, const Chart_& chart, bool invert = false)
    {
      ChartHitTestFactory<Mesh_, Chart_> factory(mesh, chart, invert);
      return factory.make();
    }

    /**
     * \brief Creates a new mesh-part from a chart hit-test function
     *
     * \param[in] mesh
     * A \transient reference to the mesh for which a mesh-part is to be created.
     *
     * \param[in] filter
     * A \resident reference to a filter meshpart. Only mesh entities of the filter
     * meshpart are tested. All other entities are considered as outside the region.
     *
     * \param[in] chart
     * A chart object whose signed-distance function is to be used for hit-testing.
     *
     * \returns
     * A mesh-part containing all entities for which the chart signed distance is <= 0.
     */
    template<typename Mesh_, typename Chart_>
    MeshPart<Mesh_> make_meshpart_by_filtered_chart_hit_test(const Mesh_& mesh, const MeshPart<Mesh_>& filter, const Chart_& chart, bool invert = false)
    {
      ChartHitTestFactory<Mesh_, Chart_> factory(mesh, filter, chart, invert);
      return factory.make();
    }

    /**
     * \brief Creates a new mesh-part from a chart hit-test function
     *
     * \param[in] mesh
     * A \transient reference to the mesh for which a mesh-part is to be created.
     *
     * \param[in] chart
     * A chart object whose signed-distance function is to be used for hit-testing.
     *
     * \returns
     * A unique-pointer to a mesh-part containing all entities for which the chart signed distance is <= 0.
     */
    template<typename Mesh_, typename Chart_>
    std::unique_ptr<MeshPart<Mesh_>> make_unique_meshpart_by_chart_hit_test(const Mesh_& mesh, const Chart_& chart, bool invert = false)
    {
      ChartHitTestFactory<Mesh_, Chart_> factory(mesh, chart, invert);
      return factory.make_unique();
    }

    /**
     * \brief Creates a new mesh-part from a chart hit-test function
     *
     * \param[in] mesh
     * A \transient reference to the mesh for which a mesh-part is to be created.
     *
     * \param[in] filter
     * A \resident reference to a filter meshpart. Only mesh entities of the filter
     * meshpart are tested. All other entities are considered as outside the region.
     *
     * \param[in] chart
     * A chart object whose signed-distance function is to be used for hit-testing.
     *
     * \returns
     * A unique-pointer to a mesh-part containing all entities for which the chart signed distance is <= 0.
     */
    template<typename Mesh_, typename Chart_>
    std::unique_ptr<MeshPart<Mesh_>> make_unique_meshpart_by_filtered_chart_hit_test(const Mesh_& mesh, const MeshPart<Mesh_>& filter, const Chart_& chart, bool invert = false)
    {
      ChartHitTestFactory<Mesh_, Chart_> factory(mesh, filter, chart, invert);
      return factory.make_unique();
    }

    /**
     * \brief Lambda Hit-Test Factory class template
     *
     * This class template can be used to create a MeshPart for a particular mesh which consists of all entities that
     * are inside the region characterized by a lambda expression hit-test function.
     *
     * \tparam Mesh_
     * The type of the mesh for which the cell sub-set is to be computed.
     *
     * \tparam Lambda_
     * A lambda expression in the X/Y/Z coordinates returning a bool that specifies whether the tested point is
     * inside the mesh-part region or not.
     *
     * \author Peter Zajac
     */
    template<typename Mesh_, typename Lambda_>
    class LambdaHitTestFactory :
      public Factory< MeshPart<Mesh_> >
    {
    public:
      /// The shape type of the mesh
      typedef typename Mesh_::ShapeType ShapeType;
      /// mesh part type
      typedef MeshPart<Mesh_> MeshType;
      /// target set holder type
      typedef typename MeshType::TargetSetHolderType TargetSetHolderType;

    protected:
      /// helper class for hit test
      class HitTestLambda
      {
      public:
        Lambda_ lambda;

        explicit HitTestLambda(Lambda_&& lb) :
          lambda(std::forward<Lambda_>(lb))
        {
        }

        /// 1D point overload
        template<typename T_, int s_>
        bool operator()(const Tiny::Vector<T_, 1, s_>& p) const
        {
          return lambda(p[0]);
        }

        /// 2D point overload
        template<typename T_, int s_>
        bool operator()(const Tiny::Vector<T_, 2, s_>& p) const
        {
          return lambda(p[0], p[1]);
        }

        /// 3D point overload
        template<typename T_, int s_>
        bool operator()(const Tiny::Vector<T_, 3, s_>& p) const
        {
          return lambda(p[0], p[1], p[3]);
        }
      } _hit_func;
      /// reference to the input mesh
      const Mesh_& _mesh;
      /// internal data storing the indices
      std::vector<std::vector<Index>> _target_data;

    public:
      /**
       * \brief Creates the factory.
       *
       * \param[in] mesh
       * A \resident reference to the mesh for which the cell sub-set is to be computed.
       *
       * \param[in] lambda
       * A lambda expression in the point coordinates for the hit-test
       */
      explicit LambdaHitTestFactory(const Mesh_& mesh,  Lambda_&& lambda) :
        _hit_func(std::forward<Lambda_>(lambda)),
        _mesh(mesh),
        _target_data(std::size_t(_mesh.shape_dim+1))
      {
        Intern::hittest_compute_target_data<HitTestLambda, Mesh_>(_target_data, _mesh, _hit_func);
      }

      /**
       * \brief Creates the factory.
       *
       * \param[in] mesh
       * A \resident reference to the mesh for which the cell sub-set is to be computed.
       *
       * \param[in] filter
       * A \resident reference to a filter meshpart. Only mesh entities of the filter
       * meshpart are tested. All other entities are considered as outside the region.
       *
       * \param[in] lambda
       * A lambda expression in the point coordinates for the hit-test
       */
      explicit LambdaHitTestFactory(const Mesh_& mesh, const MeshType& filter, Lambda_&& lambda) :
        _hit_func(std::forward<Lambda_>(lambda)),
        _mesh(mesh),
        _target_data(std::size_t(_mesh.shape_dim+1))
      {
        Intern::hittest_compute_filtered_target_data<HitTestLambda, Mesh_>(_target_data, _mesh, filter, _hit_func);
      }

      /// \copydoc Factory::get_num_entities()
      virtual Index get_num_entities(int dim) override
      {
        return Index(_target_data.at(std::size_t(dim)).size());
      }

      virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
      {
        Intern::write_to_target_set<ShapeType>(target_set_holder, _target_data);
      }

      virtual void fill_attribute_sets(typename MeshType::AttributeSetContainer&) override
      {
        // do nothing as the object has no attribute sets
      }

      virtual void fill_index_sets(std::unique_ptr<typename MeshType::IndexSetHolderType>&) override
      {
        // do nothing as the object has no index sets
      }
    }; // class LambdaHitTestFactory

    /**
     * \brief Creates a new mesh-part from a lambda hit-test function
     *
     * \param[in] mesh
     * A \transient reference to the mesh for which a mesh-part is to be created.
     *
     * \param[in] lambda
     * A lambda expression in X/Y/Z coordinates that returns \c true, if a point is inside, otherwise \c false.
     *
     * \returns
     * A mesh-part containing all entities for which the lambda evaluated to \c true.
     */
    template<typename Mesh_, typename Lambda_>
    MeshPart<Mesh_> make_meshpart_by_lambda_hit_test(const Mesh_& mesh, Lambda_&& lambda)
    {
      LambdaHitTestFactory<Mesh_, Lambda_> factory(mesh, std::forward<Lambda_>(lambda));
      return factory.make();
    }

    /**
     * \brief Creates a new mesh-part from a lambda hit-test function
     *
     * \param[in] mesh
     * A \transient reference to the mesh for which a mesh-part is to be created.
     *
     * \param[in] filter
     * A \resident reference to a filter meshpart. Only mesh entities of the filter
     * meshpart are tested. All other entities are considered as outside the region.
     *
     * \param[in] lambda
     * A lambda expression in X/Y/Z coordinates that returns \c true, if a point is inside, otherwise \c false.
     *
     * \returns
     * A mesh-part containing all entities for which the lambda evaluated to \c true.
     */
    template<typename Mesh_, typename Lambda_>
    MeshPart<Mesh_> make_meshpart_by_filtered_lambda_hit_test(const Mesh_& mesh, const MeshPart<Mesh_>& filter, Lambda_&& lambda)
    {
      LambdaHitTestFactory<Mesh_, Lambda_> factory(mesh, filter, std::forward<Lambda_>(lambda));
      return factory.make();
    }

    /**
     * \brief Creates a new mesh-part from a lambda hit-test function
     *
     * \param[in] mesh
     * A \transient reference to the mesh for which a mesh-part is to be created.
     *
     * \param[in] lambda
     * A lambda expression in X/Y/Z coordinates that returns \c true, if a point is inside, otherwise \c false.
     *
     * \returns
     * A unique-pointer to a mesh-part containing all entities for which the lambda evaluated to \c true.
     */
    template<typename Mesh_, typename Lambda_>
    std::unique_ptr<MeshPart<Mesh_>> make_unique_meshpart_by_lambda_hit_test(const Mesh_& mesh, Lambda_&& lambda)
    {
      LambdaHitTestFactory<Mesh_, Lambda_> factory(mesh, std::forward<Lambda_>(lambda));
      return factory.make_unique();
    }

    /**
     * \brief Creates a new mesh-part from a lambda hit-test function
     *
     * \param[in] mesh
     * A \transient reference to the mesh for which a mesh-part is to be created.
     *
     * \param[in] filter
     * A \resident reference to a filter meshpart. Only mesh entities of the filter
     * meshpart are tested. All other entities are considered as outside the region.
     *
     * \param[in] lambda
     * A lambda expression in X/Y/Z coordinates that returns \c true, if a point is inside, otherwise \c false.
     *
     * \returns
     * A unique-pointer to a mesh-part containing all entities for which the lambda evaluated to \c true.
     */
    template<typename Mesh_, typename Lambda_>
    std::unique_ptr<MeshPart<Mesh_>> make_unique_meshpart_by_filtered_lambda_hit_test(const Mesh_& mesh, const MeshPart<Mesh_>& filter, Lambda_&& lambda)
    {
      LambdaHitTestFactory<Mesh_, Lambda_> factory(mesh, filter, std::forward<Lambda_>(lambda));
      return factory.make_unique();
    }
  } // namespace Geometry
} // namespace FEAT
