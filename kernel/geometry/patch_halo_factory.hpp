// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/adjacency/graph.hpp>

#include <vector>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      template<typename Shape_, int codim_>
      class PatchHaloBuild
      {
        static constexpr int face_dim = Shape_::dimension - codim_;

      public:
        const TargetSet& _target_face;
        Adjacency::Graph _elem_at_face;
        std::vector<Index> _indices;

      public:
        explicit PatchHaloBuild(const TargetSetHolder<Shape_>& tsh, const IndexSetHolder<Shape_>& ish) :
          _target_face(tsh.template get_target_set<face_dim>()),
          _elem_at_face(Adjacency::RenderType::transpose, ish.template get_index_set<Shape_::dimension, face_dim>())
        {
        }

        void build(const Index halo_rank, const Adjacency::Graph& ranks_at_elem)
        {
          _indices.clear();

          // loop over all faces in the patch mesh part target set
          for(Index face(0); face < _target_face.get_num_entities(); ++face)
          {
            // get index of base-mesh face
            const Index base_face = _target_face[face];

            // check whether one of them has our desired rank
            if(_has_face_rank(base_face, halo_rank, ranks_at_elem))
            {
              _indices.push_back(face);
            }
          }
        }

        Index size() const
        {
          return Index(_indices.size());
        }

        void fill(TargetSet& target_set) const
        {
          XASSERT(target_set.get_num_entities() == size());
          for(Index i(0); i < size(); ++i)
          {
            target_set[i] = _indices[i];
          }
        }

      private:
        bool _has_face_rank(const Index base_face, const Index halo_rank, const Adjacency::Graph& ranks_at_elem)
        {
          // loop over all elements adjacent to that face and
          auto ib = _elem_at_face.image_begin(base_face);
          auto ie = _elem_at_face.image_end(base_face);
          for(auto ii(ib); ii != ie; ++ii)
          {
            // now loop over all ranks of this element
            auto jb = ranks_at_elem.image_begin(*ii);
            auto je = ranks_at_elem.image_end(*ii);
            for(auto jj(jb); jj != je; ++jj)
            {
              // check rank
              if(*jj == halo_rank)
                return true;
            }
          }
          // nope
          return false;
        }
      };

      template<typename Shape_>
      class PatchHaloBuild<Shape_, 0>
      {
      public:
        const TargetSet& _target_face;
        std::vector<Index> _indices;

      public:
        explicit PatchHaloBuild(const TargetSetHolder<Shape_>& tsh, const IndexSetHolder<Shape_>&) :
          _target_face(tsh.template get_target_set<Shape_::dimension>())
        {
        }

        void build(const Index halo_rank, const Adjacency::Graph& ranks_at_elem)
        {
          _indices.clear();

          // loop over all faces in the patch mesh part target set
          for(Index elem(0); elem < _target_face.get_num_entities(); ++elem)
          {
            // get index of base-mesh face
            const Index base_elem = _target_face[elem];

            // now loop over all ranks of this element
            auto jb = ranks_at_elem.image_begin(base_elem);
            auto je = ranks_at_elem.image_end(base_elem);
            for(auto jj(jb); jj != je; ++jj)
            {
              // check rank
              if(*jj == halo_rank)
              {
                _indices.push_back(elem);
                break;
              }
            }
          }
        }

        Index size() const
        {
          return Index(_indices.size());
        }

        void fill(TargetSet& target_set) const
        {
          XASSERT(target_set.get_num_entities() == size());
          for(Index i(0); i < size(); ++i)
          {
            target_set[i] = _indices[i];
          }
        }
      };

      template<typename Shape_, int face_dim_ = Shape_::dimension>
      class PatchHaloBuildWrapper :
        public PatchHaloBuildWrapper<Shape_, face_dim_ - 1>
      {
        typedef PatchHaloBuildWrapper<Shape_, face_dim_ - 1> BaseClass;

        static constexpr int face_codim = Shape_::dimension - face_dim_;

        PatchHaloBuild<Shape_, face_codim> _hbuild;

      public:
        explicit PatchHaloBuildWrapper(const TargetSetHolder<Shape_>& tsh, const IndexSetHolder<Shape_>& ish) :
          BaseClass(tsh, ish),
          _hbuild(tsh, ish)
        {
        }

        void build(const Index halo_rank, const Adjacency::Graph& ranks_at_elem)
        {
          BaseClass::build(halo_rank, ranks_at_elem);
          _hbuild.build(halo_rank, ranks_at_elem);
        }

        Index get_num_entities(int dim) const
        {
          XASSERT(dim <= face_dim_);
          if(dim == face_dim_)
            return Index(_hbuild.size());
          else
            return BaseClass::get_num_entities(dim);
        }

        void fill(TargetSetHolder<Shape_>& tsh) const
        {
          BaseClass::fill(tsh);
          _hbuild.fill(tsh.template get_target_set<face_dim_>());
        }
      };

      template<typename Shape_>
      class PatchHaloBuildWrapper<Shape_, 0>
      {
        static constexpr int face_codim = Shape_::dimension;

        PatchHaloBuild<Shape_, face_codim> _hbuild;

      public:
        explicit PatchHaloBuildWrapper(const TargetSetHolder<Shape_>& tsh, const IndexSetHolder<Shape_>& ish) :
          _hbuild(tsh, ish)
        {
        }

        void build(const Index halo_rank, const Adjacency::Graph& ranks_at_elem)
        {
          _hbuild.build(halo_rank, ranks_at_elem);
        }

        Index get_num_entities(int dim) const
        {
          XASSERT(dim == 0);
          return Index(_hbuild.size());
        }

        void fill(TargetSetHolder<Shape_>& tsh) const
        {
          _hbuild.fill(tsh.template get_target_set<0>());
        }
      };
    } // namespace Intern
    /// \endcond

    template<typename Mesh_>
    class PatchHaloFactory;

    /**
     * \brief Factory for creating halo mesh-parts between neighbor patches on a base-mesh.
     *
     * \attention
     * Note to self:
     * Do *NOT* switch this class from "ranks-at-elem" version back to "elems-at-rank" version,
     * because then the complexity of the PatchHaloBuild::build() function will have quadratic
     * runtime in the "many elements, few processes" case then.
     *
     * \author Peter Zajac
     */
    template<typename Shape_, int num_coords_, typename Coord_>
    class PatchHaloFactory<Geometry::ConformalMesh<Shape_, num_coords_, Coord_>> :
      public Factory<MeshPart<Geometry::ConformalMesh<Shape_, num_coords_, Coord_>>>
    {
    public:
      typedef Shape_ ShapeType;
      static constexpr int shape_dim = ShapeType::dimension;

      typedef Geometry::ConformalMesh<Shape_, num_coords_, Coord_> MeshType;
      typedef Geometry::MeshPart<MeshType> MeshPartType;

      /// Data type for attributes
      typedef typename MeshType::VertexSetType::CoordType AttributeDataType;
      /// Mesh attribute holder type
      typedef typename MeshPartType::AttributeSetContainer AttributeSetContainer;
      /// index set holder type
      typedef typename MeshPartType::IndexSetHolderType IndexSetHolderType;
      /// target set holder type
      typedef typename MeshPartType::TargetSetHolderType TargetSetHolderType;

    protected:
      const Adjacency::Graph& _ranks_at_elem;
      const MeshType& _base_mesh;
      const MeshPartType& _patch_mesh_part;
      Index _cur_halo_rank;
      Intern::PatchHaloBuildWrapper<ShapeType> _halo_wrapper;

    public:
      explicit PatchHaloFactory(
        const Adjacency::Graph& ranks_at_elem, const MeshType& base_mesh,
        const MeshPartType& patch_mesh_part) :
        _ranks_at_elem(ranks_at_elem),
        _base_mesh(base_mesh),
        _patch_mesh_part(patch_mesh_part),
        _cur_halo_rank(0),
        _halo_wrapper(patch_mesh_part.get_target_set_holder(), base_mesh.get_index_set_holder())
      {
      }

      virtual ~PatchHaloFactory()
      {
      }

      void build(Index halo_rank)
      {
        _halo_wrapper.build(halo_rank, _ranks_at_elem);
      }

      /* *************************************************************************************** */
      /* F A C T O R Y   I N T E R F A C E   I M P L E M E N T A T I O N                         */
      /* *************************************************************************************** */

      virtual Index get_num_entities(int dim) override
      {
        return _halo_wrapper.get_num_entities(dim);
      }

      virtual void fill_attribute_sets(AttributeSetContainer&) override
      {
        // nothing to do
      }

      virtual void fill_index_sets(std::unique_ptr<IndexSetHolderType>&) override
      {
        // nothing to do
      }

      virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
      {
        _halo_wrapper.fill(target_set_holder);
      }
    }; // class PatchHaloFactory<ConformalMesh<...>>
  } // namespace Geometry
} // namespace FEAT
