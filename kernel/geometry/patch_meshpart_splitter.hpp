#pragma once
#ifndef KERNEL_GEOMETRY_PATCH_MESHPART_SPLITTER_HPP
#define KERNEL_GEOMETRY_PATCH_MESHPART_SPLITTER_HPP 1

// includes, FEAST
#include <kernel/geometry/mesh_part.hpp>

// includes, system
#include <map>


namespace FEAST
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      //template<typename Shape_>
      class PatchPartMap
      {
      private:
        std::map<Index,Index> _idx_map;
        std::vector<Index> _indices;

      public:
        explicit PatchPartMap(const TargetSet& target_set) :
          _idx_map(),
          _indices()
        {
          // loop over all indices
          for(Index i(0); i < target_set.get_num_entities(); ++i)
          {
            if(!_idx_map.emplace(target_set[i], i).second)
            {
              throw InternalError("PatchMeshPartSplitter: internal error");
            }
          }
        }

        bool build(const TargetSet& target_in)
        {
          _indices.clear();
          for(Index i(0); i < target_in.get_num_entities(); ++i)
          {
            auto it = _idx_map.find(target_in[i]);
            if(it != _idx_map.end())
            {
              _indices.push_back(it->second);
            }
          }
          return !_indices.empty();
        }

        Index size() const
        {
          return _indices.size();
        }

        void fill(TargetSet& target_out) const
        {
          ASSERT_(target_out.get_num_entities() == size());
          for(Index i(0); i < target_out.get_num_entities(); ++i)
          {
            target_out[i] = _indices[i];
          }
        }
      };

      template<typename Shape_, int dim_ = Shape_::dimension>
      class PatchPartMapper :
        public PatchPartMapper<Shape_, dim_ - 1>
      {
      public:
        typedef PatchPartMapper<Shape_, dim_ - 1> BaseClass;

      private:
        PatchPartMap _patch_map;

      public:
        explicit PatchPartMapper(const TargetSetHolder<Shape_>& tsh) :
          BaseClass(tsh),
          _patch_map(tsh.template get_target_set<dim_>())
        {
        }

        bool build(const TargetSetHolder<Shape_>& tsh)
        {
          bool b1 = BaseClass::build(tsh);
          bool b2 =  _patch_map.build(tsh.template get_target_set<dim_>());
          return  b1 || b2;
        }

        Index get_num_entities(int dim) const
        {
          ASSERT_(dim <= dim_);
          if(dim == dim_)
            return Index(_patch_map.size());
          else
            return BaseClass::get_num_entities(dim);
        }

        void fill(TargetSetHolder<Shape_>& tsh) const
        {
          BaseClass::fill(tsh);
          _patch_map.fill(tsh.template get_target_set<dim_>());
        }
      };

      template<typename Shape_>
      class PatchPartMapper<Shape_, 0>
      {
      private:
        PatchPartMap _patch_map;

      public:
        explicit PatchPartMapper(const TargetSetHolder<Shape_>& tsh) :
          _patch_map(tsh.template get_target_set<0>())
        {
        }

        bool build(const TargetSetHolder<Shape_>& tsh)
        {
          return _patch_map.build(tsh.template get_target_set<0>());
        }

        Index get_num_entities(int dim) const
        {
#ifdef DEBUG
          ASSERT_(dim == 0);
#else
          (void)dim;
#endif
          return Index(_patch_map.size());
        }

        void fill(TargetSetHolder<Shape_>& tsh) const
        {
          _patch_map.fill(tsh.template get_target_set<0>());
        }
      };
    } // namespace Intern
    /// \endcond


    template<typename Mesh_>
    class PatchMeshPartSplitter;

    template<typename Shape_, int num_coords_, int stride_, typename Coord_>
    class PatchMeshPartSplitter<Geometry::ConformalMesh<Shape_, num_coords_, stride_, Coord_>> :
      public Factory<MeshPart<Geometry::ConformalMesh<Shape_, num_coords_, stride_, Coord_>>>
    {
    public:
      typedef Shape_ ShapeType;
      static constexpr int shape_dim = ShapeType::dimension;

      typedef Geometry::ConformalMesh<Shape_, num_coords_, stride_, Coord_> MeshType;
      typedef Geometry::MeshPart<MeshType> MeshPartType;

      /// Data type for attributes
      typedef typename MeshType::VertexSetType::CoordType AttributeDataType;
      /// Mesh attribute holder type
      typedef MeshAttributeHolder<ShapeType, AttributeDataType> AttributeHolderType;
      /// index set holder type
      typedef typename MeshPartType::IndexSetHolderType IndexSetHolderType;
      /// target set holder type
      typedef typename MeshPartType::TargetSetHolderType TargetSetHolderType;

    protected:
      const MeshType& _base_mesh;
      const MeshPartType& _patch_mesh_part;
      String _cur_part_name;
      Intern::PatchPartMapper<ShapeType> _part_mapper;

    public:
      explicit PatchMeshPartSplitter(const MeshType& base_mesh, const MeshPartType& patch_mesh_part) :
        _base_mesh(base_mesh),
        _patch_mesh_part(patch_mesh_part),
        _cur_part_name(),
        _part_mapper(patch_mesh_part.get_target_set_holder())
      {
      }

      virtual ~PatchMeshPartSplitter()
      {
      }

      bool build(const MeshPartType& mesh_part)
      {
        _cur_part_name = mesh_part.get_identifier();
        return _part_mapper.build(mesh_part.get_target_set_holder());
      }

      /* *************************************************************************************** */
      /* F A C T O R Y   I N T E R F A C E   I M P L E M E N T A T I O N                         */
      /* *************************************************************************************** */

      virtual String get_identifier() const override
      {
        return _cur_part_name;
      }

      virtual String get_parent_identifier() const override
      {
        return "root";
      }

      virtual Index get_num_entities(int dim) override
      {
        return _part_mapper.get_num_entities(dim);
      }

      virtual void fill_attribute_sets(AttributeHolderType&) override
      {
        // nothing to do
      }

      virtual void fill_index_sets(IndexSetHolderType*&) override
      {
        // nothing to do
      }

      virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
      {
        _part_mapper.fill(target_set_holder);
      }
    };

  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_PATCH_MESHPART_SPLITTER_HPP
