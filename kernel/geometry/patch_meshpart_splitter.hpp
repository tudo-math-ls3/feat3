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
        std::map<Index,Index> _io_map;
        std::vector<Index> _indices;

      public:
        explicit PatchPartMap(const TargetSet& target_set) :
          _idx_map(),
          _io_map(),
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
          _io_map.clear();
          for(Index i(0); i < target_in.get_num_entities(); ++i)
          {
            auto it = _idx_map.find(target_in[i]);
            if(it != _idx_map.end())
            {
              _io_map.emplace(i, Index(_indices.size()));
              _indices.push_back(it->second);
            }
          }
          return !_indices.empty();
        }

        Index size() const
        {
          return Index(_indices.size());
        }


        void fill_target_set(TargetSet& target_out) const
        {
          ASSERT_(target_out.get_num_entities() == size());
          for(Index i(0); i < target_out.get_num_entities(); ++i)
          {
            target_out[i] = _indices[i];
          }
        }

        template<typename IndexSetType_>
        void fill_index_set(IndexSetType_& index_set_out, const IndexSetType_& index_set_in, const std::map<Index, Index>& vertex_map) const
        {
          for(Index k(0); k < index_set_in.get_num_entities(); ++k)
          {
            auto kt = _io_map.find(k);
            if(kt != _io_map.end())
            {
              for(Index j(0); j < Index(index_set_in.get_num_indices()); ++j)
              {
                Index i_in(index_set_in(k,j));
                auto jt = vertex_map.find(i_in);
                if(jt != vertex_map.end())
                  index_set_out(kt->second, j) = jt->second;
                else
                  throw InternalError("Vertex "+stringify(i_in)+" missing in MeshPart topology!");
              }
            }
          }
        }

        const std::map<Index, Index>& get_io_map() const
        {
          return _io_map;
        }

      };

      template<typename Shape_, int dim_ = Shape_::dimension>
      class PatchPartMapHolder :
        public PatchPartMapHolder<Shape_, dim_ - 1>
      {
      public:
        typedef PatchPartMapHolder<Shape_, dim_ - 1> BaseClass;

      private:
        PatchPartMap _patch_map;

        bool _has_topology;

      public:
        explicit PatchPartMapHolder(const TargetSetHolder<Shape_>& tsh) :
          BaseClass(tsh),
          _patch_map(tsh.template get_target_set<dim_>()),
          _has_topology(false)
        {
        }

        bool build(const TargetSetHolder<Shape_>& tsh, const IndexSetHolder<Shape_>* ish)
        {
          ish == nullptr ? _has_topology = false : _has_topology = true;

          bool b1 = BaseClass::build(tsh, ish);
          bool b2 =  _patch_map.build(tsh.template get_target_set<dim_>());

          return  b1 || b2;
        }

        Index get_num_entities(int dim) const
        {
          ASSERT_(dim <= dim_);
          if(dim == dim_)
            return _patch_map.size();
          else
            return BaseClass::get_num_entities(dim);
        }

        bool has_topology() const
        {
          return _has_topology;
        }

        void fill_target_sets(TargetSetHolder<Shape_>& tsh) const
        {
          BaseClass::fill_target_sets(tsh);
          _patch_map.fill_target_set(tsh.template get_target_set<dim_>());
        }

        template<typename IndexSetHolder_>
        void fill_index_sets(IndexSetHolder_& ish, const IndexSetHolder_& ish_in) const
        {
          BaseClass::fill_index_sets(ish, ish_in);
          _patch_map.fill_index_set(ish.template get_index_set<dim_,0>(), ish_in.template get_index_set<dim_,0>(),
          get_vertex_map());
        }

        const std::map<Index, Index>& get_vertex_map() const
        {
          return BaseClass::get_vertex_map();
        }
      };

      template<typename Shape_>
      class PatchPartMapHolder<Shape_, 0>
      {
      private:
        PatchPartMap _patch_map;

      public:
        explicit PatchPartMapHolder(const TargetSetHolder<Shape_>& tsh) :
          _patch_map(tsh.template get_target_set<0>())
        {
        }

        bool build(const TargetSetHolder<Shape_>& tsh, const IndexSetHolder<Shape_>* DOXY(ish))
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
          return _patch_map.size();
        }

        void fill_target_sets(TargetSetHolder<Shape_>& tsh) const
        {
          _patch_map.fill_target_set(tsh.template get_target_set<0>());
        }

        template<typename IndexSetHolder_>
        void fill_index_sets(IndexSetHolder_& DOXY(ish), const IndexSetHolder_& DOXY(ish_in)) const
        {
        }

        const std::map<Index, Index>& get_vertex_map() const
        {
          return _patch_map.get_io_map();
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
      const IndexSetHolder<ShapeType>* _cur_part_topology;
      Intern::PatchPartMapHolder<ShapeType> _part_holder;

    public:
      explicit PatchMeshPartSplitter(const MeshType& base_mesh, const MeshPartType& patch_mesh_part) :
        _base_mesh(base_mesh),
        _patch_mesh_part(patch_mesh_part),
        _cur_part_name(),
        _cur_part_topology(nullptr),
        _part_holder(patch_mesh_part.get_target_set_holder())
      {
      }

      virtual ~PatchMeshPartSplitter()
      {
      }

      bool build(const MeshPartType& mesh_part)
      {
        _cur_part_name = mesh_part.get_identifier();
        _cur_part_topology = mesh_part.get_topology();
        return _part_holder.build(mesh_part.get_target_set_holder(), mesh_part.get_topology());
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
        return _part_holder.get_num_entities(dim);
      }

      virtual void fill_attribute_sets(AttributeHolderType&) override
      {
        // nothing to do
      }

      virtual void fill_index_sets(IndexSetHolderType*& index_set_holder) override
      {
        ASSERT(index_set_holder == nullptr, "fill_index_sets: index_set_holder != nullptr!");

        // If the base MeshPart has a topology, create one for the patch MeshPart, too
        if(_part_holder.has_topology())
        {
          Index num_entities[shape_dim+1];
          for(int i(0); i < shape_dim+1; ++i)
            num_entities[i] = get_num_entities(i);

          index_set_holder = new IndexSetHolderType(num_entities);
          _part_holder.fill_index_sets(*index_set_holder, *_cur_part_topology);

          //RedundantIndexSetBuilder<ShapeType>::compute(*index_set_holder);
        }

      }

      virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
      {
        _part_holder.fill_target_sets(target_set_holder);
      }
    };

  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_PATCH_MESHPART_SPLITTER_HPP
