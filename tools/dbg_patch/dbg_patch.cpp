#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/macro_factory.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/patch_factory.hpp>
#include <kernel/geometry/patch_halo_factory.hpp>
#include <kernel/geometry/test_aux/tetris_factory.hpp>

using namespace FEAST;

typedef Geometry::ConformalMesh<Shape::Quadrilateral> QuadMesh;
typedef Geometry::MeshPart<QuadMesh> QuadCellSet;
typedef Geometry::TestAux::TetrisFactory<QuadMesh> TetrisFactory;
typedef Geometry::StandardRefinery<QuadMesh> Refinery;
typedef Geometry::PatchFactory<QuadMesh> PatchFactory;
typedef Geometry::MacroFactory<QuadMesh> MacroFactory;
typedef Geometry::PatchHaloFactory<QuadCellSet> PatchMapFactory;

//*
//template<typename T_>
class PatchSetFactory :
  public Geometry::Factory<QuadCellSet>
{
public:
  enum Region
  {
    lower_patch,
    upper_patch,
    lower_interface,
    upper_interface,
    base_interface
  };

protected:
  Region _region;

public:
  explicit PatchSetFactory(Region region) :
    _region(region)
  {
  }

  virtual Index get_num_entities(int dim)
  {
    switch(_region)
    {
    case lower_patch:
    case upper_patch:
      switch(dim)
      {
      case 0:
        return 6;
      case 1:
        return 7;
      case 2:
        return 2;
      }

    case lower_interface:
    case upper_interface:
    case base_interface:
      switch(dim)
      {
      case 0:
        return 2;
      case 1:
        return 1;
      }
    }
    return 0;
  }

  virtual void fill_target_sets(TargetSetHolderType& target_set_holder)
  {
    switch(_region)
    {
    case lower_patch:
      _fill_lower_patch(target_set_holder);
      break;
    case upper_patch:
      _fill_upper_patch(target_set_holder);
      break;
    case lower_interface:
      _fill_lower_interface(target_set_holder);
      break;
    case upper_interface:
      _fill_upper_interface(target_set_holder);
      break;
    case base_interface:
      _fill_base_interface(target_set_holder);
      break;
    }
  }

  virtual void fill_attribute_sets(typename QuadCellSet::AttributeHolderType& DOXY(attribute_holder))
  {
    // do nothing as the class has no attribute sets
  }

  virtual void fill_index_sets(typename QuadCellSet::IndexSetHolderType*& DOXY(index_set_holder))
  {
    // do nothing as the class has no index sets
  }

protected:
  void _fill_lower_patch(TargetSetHolderType& target_set_holder)
  {
    Geometry::TargetSet& vi(target_set_holder.get_target_set<0>());
    vi[0] = 0;
    vi[1] = 1;
    vi[2] = 2;
    vi[3] = 4;
    vi[4] = 5;
    vi[5] = 6;

    Geometry::TargetSet& ei(target_set_holder.get_target_set<1>());
    ei[0] = 0;
    ei[1] = 1;
    ei[2] = 2;
    ei[3] = 3;
    ei[4] = 4;
    ei[5] = 6;
    ei[6] = 7;

    Geometry::TargetSet& qi(target_set_holder.get_target_set<2>());
    qi[0] = 0;
    qi[1] = 1;
  }

  void _fill_upper_patch(TargetSetHolderType& target_set_holder)
  {
    Geometry::TargetSet& vi(target_set_holder.get_target_set<0>());
    vi[0] = 3;
    vi[1] = 4;
    vi[2] = 5;
    vi[3] = 7;
    vi[4] = 8;
    vi[5] = 9;

    Geometry::TargetSet& ei(target_set_holder.get_target_set<1>());
    ei[0] = 5;
    ei[1] = 6;
    ei[2] = 8;
    ei[3] = 9;
    ei[4] = 10;
    ei[5] = 11;
    ei[6] = 12;

    Geometry::TargetSet& qi(target_set_holder.get_target_set<2>());
    qi[0] = 2;
    qi[1] = 3;
  }

  void _fill_lower_interface(TargetSetHolderType& target_set_holder)
  {
    Geometry::TargetSet& vi(target_set_holder.get_target_set<0>());
    vi[0] = 3;
    vi[1] = 4;

    Geometry::TargetSet& ei(target_set_holder.get_target_set<1>());
    ei[0] = 5;
  }

  void _fill_upper_interface(TargetSetHolderType& target_set_holder)
  {
    Geometry::TargetSet& vi(target_set_holder.get_target_set<0>());
    vi[0] = 1;
    vi[1] = 2;

    Geometry::TargetSet& ei(target_set_holder.get_target_set<1>());
    ei[0] = 1;
  }

  void _fill_base_interface(TargetSetHolderType& target_set_holder)
  {
    Geometry::TargetSet& vi(target_set_holder.get_target_set<0>());
    vi[0] = 4;
    vi[1] = 5;

    Geometry::TargetSet& ei(target_set_holder.get_target_set<1>());
    ei[0] = 6;
  }
};
//*/

int main(int /*argc*/, char** /*argv*/)
{
  // create base-mesh
  TetrisFactory factory;
  QuadMesh base_mesh(factory);

  // build patch cellsets
  PatchSetFactory patch_set_factory_l(PatchSetFactory::lower_patch);
  PatchSetFactory patch_set_factory_u(PatchSetFactory::upper_patch);
  QuadCellSet patchset_l(patch_set_factory_l);
  QuadCellSet patchset_u(patch_set_factory_u);

  // create patch meshes
  PatchFactory patch_factory_l(base_mesh, patchset_l);
  PatchFactory patch_factory_u(base_mesh, patchset_u);
  QuadMesh patch_l(patch_factory_l);
  QuadMesh patch_u(patch_factory_u);

  // create base-interface cellset
  PatchSetFactory interf_set_factory_b(PatchSetFactory::base_interface);
  QuadCellSet interf_set_b(interf_set_factory_b);

  // create patch-interface cellsets
  PatchSetFactory interf_set_factory_l(PatchSetFactory::lower_interface);
  PatchSetFactory interf_set_factory_u(PatchSetFactory::upper_interface);
  QuadCellSet interf_set_l(interf_set_factory_l);
  QuadCellSet interf_set_u(interf_set_factory_u);

  // create interface cellsets by base-interface cellset and patchsets
  PatchMapFactory interf_map_factory_l(patchset_l, interf_set_b);
  PatchMapFactory interf_map_factory_u(patchset_u, interf_set_b);
  QuadCellSet interf_map_l(interf_map_factory_l);
  QuadCellSet interf_map_u(interf_map_factory_u);

  // TODO: validate interface cellsets
}
