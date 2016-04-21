#pragma once
#ifndef KERNEL_GEOMETRY_HIT_TEST_FACTORY_HPP
#define KERNEL_GEOMETRY_HIT_TEST_FACTORY_HPP 1

// includes, FEAST
#include <kernel/geometry/index_set.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/util/tiny_algebra.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      template<
        typename HitFunc_,
        typename Mesh_,
        typename Shape_>
      struct HitTestCompute;

      template<typename Shape_>
      struct HitTestTargeter;
    } // namespace Intern
    /// \endcond

    /**
     * \brief Hit-Test Factory class template
     *
     * This class template can be used to create a MeshPart for a particular mesh,
     * which consists of all entities that are inside the region characterised by a
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
    class HitTestFactory
      : public Factory< MeshPart<Mesh_> >
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
       * A reference to the hit-test function characterising the region.
       *
       * \param[in] mesh
       * A reference to the mesh for which the cell sub-set is to be computed.
       */
      explicit HitTestFactory(const HitFunc_& hit_func, const Mesh_& mesh) :
        _hit_func(hit_func),
        _mesh(mesh),
        _target_data(std::size_t(_mesh.shape_dim+1))
      {
        // call wrapper
        Intern::HitTestCompute<HitFunc_, Mesh_, ShapeType>::wrap(_target_data, _mesh, _hit_func);
      }

      /// \copydoc Factory::get_num_entities()
      virtual Index get_num_entities(int dim) override
      {
        return Index(_target_data.at(std::size_t(dim)).size());
      }

      virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
      {
        // call wrapper
        Intern::HitTestTargeter<ShapeType>::wrap(target_set_holder, _target_data);
      }

      virtual void fill_attribute_sets(typename MeshType::MeshAttributeContainer& DOXY(attribute_container)) override
      {
        // do nothing as the object has no attribute sets
      }

      virtual void fill_index_sets(typename MeshType::IndexSetHolderType*& DOXY(index_set_holder)) override
      {
        // do nothing as the object has no index sets
      }

    }; // class HitTestFactory

    /// \cond internal
    namespace Intern
    {
      template<
        typename HitFunc_,
        typename Mesh_,
        typename Shape_>
      struct HitTestCompute
      {
        typedef std::vector<std::vector<Index>> TargetData;
        static void wrap(TargetData& target_data, const Mesh_& mesh, const HitFunc_& hit_test)
        {
          typedef typename Shape::FaceTraits<Shape_, Shape_::dimension-1>::ShapeType FacetType;
          HitTestCompute<HitFunc_, Mesh_, FacetType>::wrap(target_data, mesh, hit_test);
          apply(target_data, mesh, hit_test);
        }

        static void apply(TargetData& trg, const Mesh_& mesh, const HitFunc_& hit_test)
        {
          static constexpr int shape_dim(Shape_::dimension);
          const Index num_cells(mesh.get_num_entities(shape_dim));
          const IndexSet<Shape::FaceTraits<Shape_, 0>::count> &index_set(mesh.template get_index_set<shape_dim,0>());
          trg[shape_dim].reserve(num_cells);
          for(Index i(0); i < num_cells; ++i)
          {
            if (hit_test(get_midpoint(mesh.get_vertex_set(), index_set[i])))
            {
              trg[shape_dim].push_back(i);
            }
          }
        }

        template<
          typename VertexSet_,
          typename IndexTuple_>
        static Tiny::Vector<
          typename VertexSet_::CoordType,
          VertexSet_::num_coords>
        get_midpoint(
          const VertexSet_& vertex_set,
          const IndexTuple_& index_tuple)
        {
          typedef typename VertexSet_::CoordType DataType;
          int num_vertex(Shape::FaceTraits<Shape_,0>::count);
          Tiny::Vector<DataType, VertexSet_::num_coords> mid_point;
          mid_point.format();
          for(int i(0); i < num_vertex; ++i)
          {
            for(int j(0); j < VertexSet_::num_coords; ++j)
            {
              mid_point[j] += vertex_set[index_tuple[i]][j];
            }
          }
          return mid_point * (DataType(1) / DataType(num_vertex));
        }
      };

      template<
        typename HitFunc_,
        typename Mesh_>
      struct HitTestCompute<HitFunc_, Mesh_, Shape::Vertex>
      {
        typedef std::vector<std::vector<Index>> TargetData;
        static void wrap(TargetData& target_data, const Mesh_& mesh, const HitFunc_& hit_test)
        {
          apply(target_data, mesh, hit_test);
        }

        static void apply(TargetData& trg, const Mesh_& mesh, const HitFunc_& hit_test)
        {
          const Index num_cells(mesh.get_num_entities(0));
          trg[0].reserve(num_cells);
          for(Index i(0); i < num_cells; ++i)
          {
            if (hit_test(get_midpoint(mesh.get_vertex_set(), i)))
            {
              trg[0].push_back(i);
            }
          }
        }

        template<typename VertexSet_>
        static Tiny::Vector<
          typename VertexSet_::CoordType,
          VertexSet_::num_coords>
        get_midpoint(const VertexSet_& vertex_set, const Index vertex)
        {
          typedef typename VertexSet_::CoordType DataType;
          Tiny::Vector<DataType, VertexSet_::num_coords> mid_point;
          for(int j(0); j < VertexSet_::num_coords; ++j)
          {
            mid_point[j] = vertex_set[vertex][j];
          }
          return mid_point;
        }
      };

      template<typename Shape_>
      struct HitTestTargeter
      {
        typedef std::vector<std::vector<Index>> SourceData;
        static void wrap(TargetSetHolder<Shape_>& tsh, SourceData& sd)
        {
          typedef typename Shape::FaceTraits<Shape_, Shape_::dimension-1>::ShapeType FacetType;
          HitTestTargeter<FacetType>::wrap(tsh, sd);
          apply(tsh.template get_target_set<Shape_::dimension>(), sd[Shape_::dimension]);
        }

        static void apply(TargetSet& trg, std::vector<Index>& sd)
        {
          const Index num_cells(Index(sd.size()));
          for(Index i(0); i < num_cells; ++i)
          {
            trg[i] = sd[i];
          }
        }
      };

      template<>
      struct HitTestTargeter<Shape::Vertex>
      {
        typedef std::vector<std::vector<Index>> SourceData;
        static void wrap(TargetSetHolder<Shape::Vertex>& tsh, SourceData& sd)
        {
          apply(tsh.get_target_set<0>(), sd[0]);
        }

        static void apply(TargetSet& trg, std::vector<Index>& sd)
        {
          const Index num_cells(Index(sd.size()));
          for(Index i(0); i < num_cells; ++i)
          {
            trg[i] = sd[i];
          }
        }
      };
    } // namespace Intern
    /// \endcond

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
    template<typename DataType_, Index dim_>
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
        ASSERT_(radius > DataType_(0));
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
        ASSERT_(radius > DataType_(0));
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
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_HIT_TEST_FACTORY_HPP
