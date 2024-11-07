// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <iostream>
#include <fstream>

#ifdef FEAT_HAVE_CGAL

#define CGAL_HEADER_ONLY 1

#ifdef FEAT_COMPILER_MICROSOFT
#pragma warning(disable: 4244)
#endif // FEAT_COMPILER_MICROSOFT

FEAT_DISABLE_WARNINGS
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/algorithm.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/boost/graph/IO/OBJ.h>
#include <CGAL/boost/graph/IO/OFF.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
FEAT_RESTORE_WARNINGS

// check if cgal threading support is enabled
#ifndef CGAL_HAS_THREADS
static_assert(false, "No cgal multithreading support");
#endif

#ifndef BOOST_HAS_THREADS
static_assert(false, "Boost has no threads");
#endif


#ifdef FEAT_COMPILER_INTEL
#pragma warning disable 3280
#pragma warning disable 1224
#endif //FEAT_COMPILER_INTEL

#ifdef FEAT_COMPILER_MICROSOFT
#pragma warning(default: 4244)
#endif // FEAT_COMPILER_MICROSOFT

#include <kernel/geometry/cgal.hpp>

#ifdef FEAT_HAVE_HALFMATH
#include <kernel/util/half.hpp>
#endif

template<typename DT_>
struct CGALTypeWrapper
{
  typedef DT_ CGALDT_;
  typedef typename CGAL::Simple_cartesian<DT_> K;
  typedef typename K::Point_3 Point_;
  typedef typename CGAL::Polyhedron_3<K> Polyhedron_;
  typedef typename CGAL::AABB_face_graph_triangle_primitive<Polyhedron_> Primitive_;
  typedef typename CGAL::AABB_traits<K, Primitive_> Traits_;
  typedef typename CGAL::AABB_tree<Traits_> Tree_;
  typedef typename CGAL::Side_of_triangle_mesh<Polyhedron_, K> Point_inside_;
  typedef typename CGAL::Aff_transformation_3<K> Transformation_;
  typedef typename CGAL::AABB_traits<K, Primitive_>::Point_and_primitive_id Point_and_primitive_id_;
};

template<typename DT_>
struct CGALWrapperData
{
  typedef CGALTypeWrapper<DT_> TW_;
  typename TW_::Polyhedron_ * _polyhedron;
  typename TW_::Tree_ * _tree;
  typename TW_::Point_inside_ * _inside_tester;

  CGALWrapperData() :
    _polyhedron(nullptr),
    _tree(nullptr),
    _inside_tester(nullptr)
  {
  }

  CGALWrapperData(const CGALWrapperData&) = delete;
  CGALWrapperData& operator=(const CGALWrapperData&) = delete;

  void swap(CGALWrapperData& other)
  {
    std::swap(this->_polyhedron, other._polyhedron);
    std::swap(this->_tree, other._tree);
    std::swap(this->_inside_tester, other._inside_tester);
  }

  CGALWrapperData(CGALWrapperData&& other) noexcept
  : _polyhedron(nullptr), _tree(nullptr), _inside_tester(nullptr)
  {
    this->swap(other);
  }

  CGALWrapperData& operator=(CGALWrapperData&& other) noexcept
  {
    this->swap(other);
    return *this;
  }


  ~CGALWrapperData()
  {
    if(_inside_tester)
      delete _inside_tester;
    if(_tree)
      delete _tree;
    if(_polyhedron)
      delete _polyhedron;
  }
};

#ifdef FEAT_HAVE_QUADMATH
template<>
struct CGALTypeWrapper<__float128>
{
  typedef double CGALDT_;
  typedef typename CGAL::Simple_cartesian<double> K;
  typedef typename K::Point_3 Point_;
  typedef typename CGAL::Polyhedron_3<K> Polyhedron_;
  typedef typename CGAL::AABB_face_graph_triangle_primitive<Polyhedron_> Primitive_;
  typedef typename CGAL::AABB_traits<K, Primitive_> Traits_;
  typedef typename CGAL::AABB_tree<Traits_> Tree_;
  typedef typename CGAL::Side_of_triangle_mesh<Polyhedron_, K> Point_inside_;
  typedef typename CGAL::Aff_transformation_3<K> Transformation_;
  typedef typename CGAL::AABB_traits<K, Primitive_>::Point_and_primitive_id Point_and_primitive_id_;
};
#endif

template<typename DT_>
using PointTypeAlias = typename FEAT::Geometry::template CGALWrapper<DT_>::PointType;

template<typename DT_>
using TransformMatrixAlias = typename FEAT::Geometry::template CGALWrapper<DT_>::TransformMatrix;

template<typename DT_>
FEAT::Geometry::CGALWrapper<DT_>::CGALWrapper(std::istream & file, CGALFileMode file_mode) :
  _cgal_data(nullptr)
{
  _parse_mesh(file, file_mode);
}

template<typename DT_>
FEAT::Geometry::CGALWrapper<DT_>::CGALWrapper(const String & filename, CGALFileMode file_mode) :
  _cgal_data(nullptr)
{
  std::ifstream file(filename.c_str());
  XASSERTM(file.is_open(), "CGALWrapper: file read error: " + filename + " !");
  _parse_mesh(file, file_mode);
  file.close();
}

template<typename DT_>
FEAT::Geometry::CGALWrapper<DT_>::CGALWrapper::~CGALWrapper()
{
  if(_cgal_data)
    delete (CGALWrapperData<DT_>*)_cgal_data;
}

template<typename DT_>
bool FEAT::Geometry::CGALWrapper<DT_>::point_inside(DT_ x, DT_ y, DT_ z) const
{
  typedef typename CGALTypeWrapper<DT_>::CGALDT_ IDT_;
  typename CGALTypeWrapper<DT_>::Point_ query{IDT_(x), IDT_(y), IDT_(z)};

  // Determine the side and return true if inside!
  CGALWrapperData<DT_> * cd = (CGALWrapperData<DT_>*)_cgal_data;
  return (*(cd->_inside_tester))(query) == CGAL::ON_BOUNDED_SIDE;
}

template<typename DT_>
DT_ FEAT::Geometry::CGALWrapper<DT_>::squared_distance(DT_ x, DT_ y, DT_ z) const
{
  typedef typename CGALTypeWrapper<DT_>::CGALDT_ IDT_;
  typename CGALTypeWrapper<DT_>::Point_ query{IDT_(x), IDT_(y), IDT_(z)};

  CGALWrapperData<DT_> * cd = (CGALWrapperData<DT_>*)_cgal_data;
  return DT_(cd->_tree->squared_distance(query));
}

template<typename DT_>
typename FEAT::Geometry::CGALWrapper<DT_>::PointType FEAT::Geometry::CGALWrapper<DT_>::closest_point(const PointType& point) const
{
  return closest_point(point[0], point[1], point[2]);
}

template<typename DT_>
typename FEAT::Geometry::CGALWrapper<DT_>::PointType FEAT::Geometry::CGALWrapper<DT_>::closest_point(DT_ x, DT_ y, DT_ z) const
{
  typedef typename CGALTypeWrapper<DT_>::CGALDT_ IDT_;
  typename CGALTypeWrapper<DT_>::Point_ query{IDT_(x), IDT_(y), IDT_(z)};
  CGALWrapperData<DT_> * cd = (CGALWrapperData<DT_>*)_cgal_data;
  typename CGALTypeWrapper<DT_>::Point_ n_point = cd->_tree->closest_point(query);
  return PointType{{DT_(n_point[0]), DT_(n_point[1]), DT_(n_point[2])}};
}

template<typename DT_>
typename FEAT::Geometry::CGALWrapper<DT_>::PointType FEAT::Geometry::CGALWrapper<DT_>::closest_point(const PointType& point, PointType& primitive_grad) const
{
  typedef typename CGALTypeWrapper<DT_>::CGALDT_ IDT_;
  typename CGALTypeWrapper<DT_>::Point_ query{IDT_(point[0]), IDT_(point[1]), IDT_(point[2])};
  CGALWrapperData<DT_> * cd = (CGALWrapperData<DT_>*)_cgal_data;
  const typename CGALTypeWrapper<DT_>::Point_and_primitive_id_& cl_p_query = cd->_tree->closest_point_and_primitive(query);
  typename CGALTypeWrapper<DT_>::Polyhedron_::Face_handle f = cl_p_query.second;
  auto v = CGAL::Polygon_mesh_processing::compute_face_normal(f,*(cd->_polyhedron));
  for(int i = 0; i < 3; ++i)
    primitive_grad[i] = DT_(v[i]);
  return PointType{{DT_(cl_p_query.first[0]), DT_(cl_p_query.first[1]), DT_(cl_p_query.first[2])}};
}

template<typename DT_>
void FEAT::Geometry::CGALWrapper<DT_>::_parse_mesh(std::istream & file, CGALFileMode file_mode)
{
  delete (CGALWrapperData<DT_>*)_cgal_data;

  _cgal_data = new CGALWrapperData<DT_>;
  CGALWrapperData<DT_> * cd = (CGALWrapperData<DT_>*)_cgal_data;

  cd->_polyhedron = new typename CGALTypeWrapper<DT_>::Polyhedron_;

  bool status(false);
  switch(file_mode)
  {
    case CGALFileMode::fm_off:
      status = CGAL::IO::read_OFF(file, *(cd->_polyhedron));
      XASSERTM(status == true, "CGAL::IO::read_OFF: read/parse error !");
      break;
    case CGALFileMode::fm_obj:
      status = CGAL::IO::read_OBJ(file, *(cd->_polyhedron));
      XASSERTM(status == true, "CGAL::IO::read_OBJ: read/parse error !");
      break;
    default:
      XASSERTM(false, "CGAL FileMode not supported!");
  }

  _init_wrapper();

}

template<typename DT_>
void FEAT::Geometry::CGALWrapper<DT_>::transform(const TransformMatrix& scale_rot, const PointType& translation, DT_ scale)
{
  typedef typename CGALTypeWrapper<DT_>::CGALDT_ IDT_;
  CGALWrapperData<DT_> * cd = (CGALWrapperData<DT_>*)_cgal_data;
  //delete tree and _inside_tester
  _delete_tree();
  //build affine transformation
  typename CGALTypeWrapper<DT_>::Transformation_ trafo_mat{(IDT_)(scale_rot[0][0]), (IDT_)(scale_rot[0][1]), (IDT_)(scale_rot[0][2]), (IDT_)(translation[0]),
                            (IDT_)(scale_rot[1][0]), (IDT_)(scale_rot[1][1]), (IDT_)(scale_rot[1][2]), (IDT_)(translation[1]),
                            (IDT_)(scale_rot[2][0]), (IDT_)(scale_rot[2][1]), (IDT_)(scale_rot[2][2]), (IDT_)(translation[2]),
                            IDT_(scale)};
  //apply affine transformation on polyhedron
  CGAL::Polygon_mesh_processing::transform(trafo_mat, *(cd->_polyhedron));
  //create new trees
  _init_wrapper();


}

/// TODO: Implement me
template<typename DT_>
std::size_t FEAT::Geometry::CGALWrapper<DT_>::bytes() const
{
  return 0;
}

template<typename DT_>
void FEAT::Geometry::CGALWrapper<DT_>::_delete_tree()
{
  CGALWrapperData<DT_> * cd = (CGALWrapperData<DT_>*)_cgal_data;
  XASSERTM(cd->_polyhedron != nullptr, "ERROR: Polyhedron is not initialized!");
  delete cd->_inside_tester;
  cd->_inside_tester = nullptr;
  delete cd->_tree;
  cd->_tree = nullptr;
}

template<typename DT_>
void FEAT::Geometry::CGALWrapper<DT_>::_init_wrapper()
{
  CGALWrapperData<DT_> * cd = (CGALWrapperData<DT_>*)_cgal_data;
  XASSERTM(cd->_polyhedron != nullptr, "ERROR: Polyhedron is not initialized!");
  XASSERTM(cd->_tree == nullptr && cd->_inside_tester == nullptr, "ERROR: Tree or Inside Tester are already initialized");

    // Construct AABB tree with a KdTree
  cd->_tree = new typename CGALTypeWrapper<DT_>::Tree_(faces(*(cd->_polyhedron)).first, faces(*(cd->_polyhedron)).second, *(cd->_polyhedron));
  cd->_tree->accelerate_distance_queries();
  // Initialize the point-in-polyhedron tester
  cd->_inside_tester = new typename CGALTypeWrapper<DT_>::Point_inside_(*(cd->_tree));
}

template<typename DT_>
FEAT::Geometry::CGALWrapper<DT_>::CGALWrapper(FEAT::Geometry::CGALWrapper<DT_>&& other) noexcept
 : _cgal_data(nullptr)
{
  std::swap(this->_cgal_data, other._cgal_data);
}

template<typename DT_>
FEAT::Geometry::CGALWrapper<DT_>& FEAT::Geometry::CGALWrapper<DT_>::operator=(FEAT::Geometry::CGALWrapper<DT_>&& other) noexcept
{
  std::swap(this->_cgal_data, other._cgal_data);
  return *this;
}


// explicitly instantiate templates for all sensible datatype

template class FEAT::Geometry::CGALWrapper<double>;
template class FEAT::Geometry::CGALWrapper<float>;
#ifdef FEAT_HAVE_HALFMATH
template class FEAT::Geometry::CGALWrapper<Half>;
#endif
#ifdef FEAT_HAVE_QUADMATH
template class FEAT::Geometry::CGALWrapper<__float128>;
#endif


#elif !defined(DOXYGEN)

// Dummy class instance to silence ipo linker optimization warnings about empty libgeometry
class ipo_foobar_geometry_cgal
{
public:
  int i;
  ipo_foobar_geometry_cgal() :
    i(0)
  {
    (void)i;
  }
} ipo_barfoo_geometry_cgal;

#endif // FEAT_HAVE_CGAL
