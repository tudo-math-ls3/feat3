

#pragma once

#include <kernel/geometry/cgal.hpp>
#include <kernel/geometry/voxel_map.hpp>
#include <kernel/geometry/hit_test_factory.hpp>
#include <kernel/analytic/lambda_function.hpp>
#include <kernel/analytic/distance_function.hpp>
#include <kernel/util/dist_file_io.hpp>


namespace Gendie
{
  template<typename DT_>
  class CGALHitTestFunction
  {
  public:
    typedef DT_ DataType;
    const FEAT::Geometry::CGALWrapper<DataType>& cgal;

    explicit CGALHitTestFunction(const FEAT::Geometry::CGALWrapper<DataType>& cg) : cgal(cg) {}

    template<typename DT2_, int n_, int sn_>
    bool operator()(const FEAT::Tiny::Vector<DT2_, n_, sn_>& point) const
    {
      if constexpr(n_ == 3)
        return !cgal.point_inside(DataType(point[0]), DataType(point[1]), DataType(point[2]));
      else
      {
        XABORTM("Called CGAL hit test in 2D");
        return false;
      }
    }
  }; // class CGALHitTestFunction<...>

  class VoxelMapHitTestFunction
  {
  public:
    const FEAT::Geometry::VoxelMap& vxl_map;
    const Real tol;

    explicit VoxelMapHitTestFunction(const FEAT::Geometry::VoxelMap& vxl) : vxl_map(vxl), tol(Math::pow(Math::Limits<Real>::epsilon(), Real(0.7))) {}

    template<typename DT2_, int n_, int sn_>
    bool operator()(const FEAT::Tiny::Vector<DT2_, n_, sn_>& point) const
    {
      // return (vxl_map.sample_point(point) < tol);
      return !vxl_map.check_point_nearest(point);
    }
  }; // class VoxelMapHitTestFunction<...>

  /// Collective call
  template<typename DT_ = FEAT::Real>
  inline std::unique_ptr<FEAT::Geometry::CGALWrapper<DT_>> create_cgal_wrapper(const FEAT::String& filename, const FEAT::Dist::Comm& comm)
  {
    std::stringstream iss;
    FEAT::DistFileIO::read_common(iss, filename, comm);
    return std::make_unique<FEAT::Geometry::CGALWrapper<DT_>>(iss, FEAT::Geometry::CGALFileMode::fm_off);
  }

  template<typename Mesh_, typename DT_>
  inline std::unique_ptr<FEAT::Geometry::MeshPart<Mesh_>> make_cgal_mesh_part(const Mesh_& mesh, const FEAT::Geometry::CGALWrapper<DT_>& cgal)
  {
    CGALHitTestFunction cghtf(cgal);
    FEAT::Geometry::HitTestFactory factory(cghtf, mesh);
    return factory.make_unique();
  }

  template<typename Mesh_>
  inline std::unique_ptr<FEAT::Geometry::MeshPart<Mesh_>> make_voxel_mesh_part(const Mesh_& mesh, const FEAT::Geometry::VoxelMap& vxl_map)
  {
    VoxelMapHitTestFunction vxlhtf(vxl_map);
    FEAT::Geometry::HitTestFactory factory(vxlhtf, mesh);
    return factory.make_unique();
  }

  template<typename MeshNode_, typename DT_>
  inline typename MeshNode_::MeshPartNodeType* add_cgal_mesh_part(MeshNode_* mesh_node, const FEAT::Geometry::CGALWrapper<DT_>* cgal_wrapper,
                                                            const FEAT::String& part_name = "gendie_fbm")
  {
    if(!mesh_node)
      return nullptr;
    if(!cgal_wrapper)
      return nullptr;

    return mesh_node->add_mesh_part(part_name, make_cgal_mesh_part(*mesh_node->get_mesh(), *cgal_wrapper));
  }

  template<typename MeshNode_>
  inline typename MeshNode_::MeshPartNodeType* add_voxel_mesh_part(MeshNode_* mesh_node, const FEAT::Geometry::VoxelMap& vxl_map,
                                                            const FEAT::String& part_name)
  {
    if(!mesh_node)
      return nullptr;

    return mesh_node->add_mesh_part(part_name, make_voxel_mesh_part(*mesh_node->get_mesh(), vxl_map));
  }
}