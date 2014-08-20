#pragma once
#ifndef KERNEL_FOUNDATION_MESH_UTIL_HPP
#define KERNEL_FOUNDATION_MESH_UTIL_HPP 1

#include<kernel/foundation/mesh.hpp>
#include<kernel/foundation/stl_util.hpp>
#include<cmath>

namespace FEAST
{
  namespace Foundation
  {
    enum EdgeTypes
    {
      et_iz_x = 0,
      et_c_x,
      et_iz_y,
      et_c_y
    };

    struct MeshUtil
    {
      ///for QUAD 2D meshes
      template<typename TopologyType_,
               template <typename, typename> class OuterStorageType_,
               typename AT_>
      static bool iz_property(const Mesh<Dim2D, TopologyType_, OuterStorageType_>& m, const AT_& x, const AT_& y)
      {
        for(Index i(0) ; i < m.num_polytopes(pl_face) ; ++i)
        {
          auto v_fi(m.get_adjacent_polytopes(pl_face, pl_vertex, i));

          if(v_fi.size() != 4)
          {
            std::cout << "WARNING: not a pure quad mesh!" << std::endl;
          }

          typename AT_::data_type_ e0_x(x.at(v_fi.at(1)) - x.at(v_fi.at(0)));
          typename AT_::data_type_ e0_y(y.at(v_fi.at(1)) - y.at(v_fi.at(0)));
          typename AT_::data_type_ ez_x(x.at(v_fi.at(2)) - x.at(v_fi.at(1)));
          typename AT_::data_type_ ez_y(y.at(v_fi.at(2)) - y.at(v_fi.at(1)));

          ///check sanity of vertex-based iz-curve (cross-edges)
          if(e0_x > 0. && e0_y > 0)
          {
            if(!(ez_x < 0.))
            {
              std::cout << "WARNING: malformed cross-edge in iz-curve! e0_x > 0 AND e0_y > 0 => ez_x < 0, but ez_x is " << ez_x <<  "!" << std::endl;
              return false;
            }
          }
          else if(e0_x > 0. && e0_y <= 0)
          {
            if(!(ez_y > 0.))
            {
              std::cout << "WARNING: malformed cross-edge in iz-curve! e0_x > 0 AND e0_y <= 0 => ez_y > 0, but ez_y is " << ez_y << "!" << std::endl;
              return false;
            }
          }
          else if(e0_x <= 0. && e0_y > 0.)
          {
            if(!(ez_y < 0.))
            {
              std::cout << "WARNING: malformed cross-edge in iz-curve! e0_x <= 0 AND e0_y > 0 => ez_y < 0, but ez_y is " << ez_y << "!" << std::endl;
              return false;
            }
          }
          else if(e0_x <= 0. && e0_y <= 0)
          {
            if(!(ez_x > 0.))
            {
              std::cout << "WARNING: malformed cross-edge in iz-curve! e0_x <= 0 AND e0_y <= 0 => ez_x > 0, but ez_x is " << ez_x << "!" << std::endl;
              return false;
            }
          }

          ///check existence of iz-edges
          ///check existence of completion edges
          auto e_fi(m.get_adjacent_polytopes(pl_face, pl_edge, i));
          bool found_e0(false);
          bool found_e1(false);
          bool found_c0(false);
          bool found_c1(false);
          for(auto e_fi_j : e_fi)
          {
            auto v_e_fi_j(m.get_adjacent_polytopes(pl_edge, pl_vertex, e_fi_j));
            found_e0 = v_e_fi_j.at(0) == v_fi.at(0) && v_e_fi_j.at(1) == v_fi.at(1) ? true : found_e0;
            found_e1 = v_e_fi_j.at(0) == v_fi.at(2) && v_e_fi_j.at(1) == v_fi.at(3) ? true : found_e1;
            found_c0 = v_e_fi_j.at(0) == v_fi.at(0) && v_e_fi_j.at(1) == v_fi.at(2) ? true : found_c0;
            found_c1 = v_e_fi_j.at(0) == v_fi.at(1) && v_e_fi_j.at(1) == v_fi.at(3) ? true : found_c1;
          }
          if(!(found_e0 && found_e1 && found_c0 && found_c1))
          {
            if(!found_e0)
              std::cout << "WARNING: no matching iz-edge to iz-curve (e0) at face " << i << "!" << std::endl;
            if(!found_e1)
              std::cout << "WARNING: no matching iz-edge to iz-curve (e1) at face " << i << "!" << std::endl;
            if(!found_c0)
              std::cout << "WARNING: no matching completion-edge to iz-curve (c0) at face " << i << "!" << std::endl;
            if(!found_c1)
              std::cout << "WARNING: no matching completion-edge to iz-curve (c1) at face " << i << "!" << std::endl;

            return false;
          }
        }
        return true;
      }

      template<typename TopologyType_,
               template <typename, typename> class OuterStorageType_,
               typename AT_>
      static void establish_iz_property(Mesh<Dim2D, TopologyType_, OuterStorageType_>& m, const AT_& x, const AT_& y)
      {
        if(iz_property(m, x, y))
          return;

        std::cout << "Establishing iz-property..." << std::endl;

        OuterStorageType_<Index, std::allocator<Index> > faces_processed;
        OuterStorageType_<Index, std::allocator<Index> > edges_processed;
        OuterStorageType_<EdgeTypes, std::allocator<EdgeTypes> > edge_types;

        ///start with face 0
        auto E_f0(m.get_adjacent_polytopes(pl_face, pl_edge, 0));
        ///pick edge 0 and find e1 = E_f0 - E_V_e0
        auto e0 = E_f0.at(0);
        std::sort(E_f0.begin(), E_f0.end());
        auto V_e0(m.get_adjacent_polytopes(pl_edge, pl_vertex, e0));
        decltype(V_e0) E_V_e0;
        for(auto V_e0_j : V_e0)
        {
          auto E_V_e0_j(m.get_adjacent_polytopes(pl_vertex, pl_edge, V_e0_j));
          std::sort(E_V_e0_j.begin(), E_V_e0_j.end());

          decltype(E_V_e0) tmp(E_V_e0.size() + E_V_e0_j.size());
          auto iter(std::set_union(E_V_e0.begin(), E_V_e0.end(), E_V_e0_j.begin(), E_V_e0_j.end(), tmp.begin()));
          tmp.resize(Index(iter - tmp.begin()));
          E_V_e0 = tmp;
        }
        decltype(V_e0) E_f0_minus_E_V_e0(E_f0.size());
        std::set_difference(E_f0.begin(), E_f0.end(), E_V_e0.begin(), E_V_e0.end(), E_f0_minus_E_V_e0.begin());
        auto e1(E_f0_minus_E_V_e0.at(0));
        auto V_e1(m.get_adjacent_polytopes(pl_edge, pl_vertex, e1));


        ///direct e0, e1 to positive x direction (if possible), heuristics: take smallest y-coord-sum-edge for real ez0
        auto x_diff_ez0(x.at(V_e0.at(1)) - x.at(V_e0.at(0)));
        auto x_diff_ez1(x.at(V_e1.at(1)) - x.at(V_e1.at(0)));
        auto y_diff_ez0(y.at(V_e0.at(1)) - y.at(V_e0.at(0)));
        auto y_diff_ez1(y.at(V_e1.at(1)) - y.at(V_e1.at(0)));
        if(x_diff_ez0 != 0) //x pos mode
        {
          auto y_sum_e0(y.at(V_e0.at(0)) + y.at(V_e0.at(1)));
          auto y_sum_e1(y.at(V_e1.at(0)) + y.at(V_e1.at(1)));
          auto ez0(y_sum_e0 < y_sum_e1 ? e0 : e1);
          edges_processed.push_back(ez0);
          edge_types.push_back(et_iz_x);

          auto ez1(y_sum_e0 < y_sum_e1 ? e1 : e0);
          edges_processed.push_back(ez1);
          edge_types.push_back(et_iz_x);

          if(x_diff_ez0 < 0)
          {
            m.get_topologies().at(ipi_edge_vertex).at(ez0).at(0) = V_e0.at(1);
            m.get_topologies().at(ipi_edge_vertex).at(ez0).at(1) = V_e0.at(0);
          }

          if(x_diff_ez1 < 0)
          {
            m.get_topologies().at(ipi_edge_vertex).at(ez1).at(0) = V_e1.at(1);
            m.get_topologies().at(ipi_edge_vertex).at(ez1).at(1) = V_e1.at(0);
          }

          ///find completion edges
          auto v_ec0_0(ez0 == e0 ? V_e0.at(0) : V_e1.at(0));
          auto v_ec0_1(ez0 == e0 ? V_e1.at(0) : V_e0.at(0));
          auto v_ec1_0(ez0 == e0 ? V_e0.at(1) : V_e1.at(1));
          auto v_ec1_1(ez0 == e0 ? V_e1.at(1) : V_e0.at(1));

          for(auto E_f0_j : E_f0)
          {
            auto V_E_f0_j(m.get_adjacent_polytopes(pl_edge, pl_vertex, E_f0_j));

            auto iter00(std::find(V_E_f0_j.begin(), V_E_f0_j.end(), v_ec0_0));
            auto iter01(std::find(V_E_f0_j.begin(), V_E_f0_j.end(), v_ec0_1));

            auto iter10(std::find(V_E_f0_j.begin(), V_E_f0_j.end(), v_ec1_0));
            auto iter11(std::find(V_E_f0_j.begin(), V_E_f0_j.end(), v_ec1_1));

            if(iter00 != V_E_f0_j.end() && iter01 != V_E_f0_j.end())
            {
              if(x.at(V_E_f0_j.at(0)) > x.at(V_E_f0_j.at(1)))
              {
                m.get_topologies().at(ipi_edge_vertex).at(E_f0_j).at(0) = V_E_f0_j.at(1);
                m.get_topologies().at(ipi_edge_vertex).at(E_f0_j).at(1) = V_E_f0_j.at(0);
              }

              edges_processed.push_back(E_f0_j);
              edge_types.push_back(et_c_y);
            }

            if(iter10 != V_E_f0_j.end() && iter11 != V_E_f0_j.end())
            {
              if(x.at(V_E_f0_j.at(0)) > x.at(V_E_f0_j.at(1)))
              {
                m.get_topologies().at(ipi_edge_vertex).at(E_f0_j).at(0) = V_E_f0_j.at(1);
                m.get_topologies().at(ipi_edge_vertex).at(E_f0_j).at(1) = V_E_f0_j.at(0);
              }
              edges_processed.push_back(E_f0_j);
              edge_types.push_back(et_c_y);
            }
          }

          ///set iz-curve
          auto V_ez0(m.get_adjacent_polytopes(pl_edge, pl_vertex, ez0));
          auto V_ez1(m.get_adjacent_polytopes(pl_edge, pl_vertex, ez1));

          m.get_topologies().at(ipi_face_vertex).at(0).at(0) = V_ez0.at(0);
          m.get_topologies().at(ipi_face_vertex).at(0).at(1) = V_ez0.at(1);
          m.get_topologies().at(ipi_face_vertex).at(0).at(2) = V_ez1.at(0);
          m.get_topologies().at(ipi_face_vertex).at(0).at(3) = V_ez1.at(1);
        }
        else //y pos mode
        {
          auto x_sum_e0(x.at(V_e0.at(0)) + x.at(V_e0.at(1)));
          auto x_sum_e1(x.at(V_e1.at(0)) + x.at(V_e1.at(1)));
          auto ez0(x_sum_e0 > x_sum_e1 ? e0 : e1);
          edges_processed.push_back(ez0);
          edge_types.push_back(et_iz_y);

          auto ez1(x_sum_e0 > x_sum_e1 ? e1 : e0);
          edges_processed.push_back(ez1);
          edge_types.push_back(et_iz_y);

          if(y_diff_ez0 < 0)
          {
            m.get_topologies().at(ipi_edge_vertex).at(ez0).at(0) = V_e0.at(1);
            m.get_topologies().at(ipi_edge_vertex).at(ez0).at(1) = V_e0.at(0);
          }

          if(y_diff_ez1 < 0)
          {
            m.get_topologies().at(ipi_edge_vertex).at(ez1).at(0) = V_e1.at(1);
            m.get_topologies().at(ipi_edge_vertex).at(ez1).at(1) = V_e1.at(0);
          }

          ///find completion edges
          auto v_ec0_0(ez0 == e0 ? V_e0.at(0) : V_e1.at(0));
          auto v_ec0_1(ez0 == e0 ? V_e1.at(0) : V_e0.at(0));
          auto v_ec1_0(ez0 == e0 ? V_e0.at(1) : V_e1.at(1));
          auto v_ec1_1(ez0 == e0 ? V_e1.at(1) : V_e0.at(1));

          for(auto E_f0_j : E_f0)
          {
            auto V_E_f0_j(m.get_adjacent_polytopes(pl_edge, pl_vertex, E_f0_j));

            auto iter00(std::find(V_E_f0_j.begin(), V_E_f0_j.end(), v_ec0_0));
            auto iter01(std::find(V_E_f0_j.begin(), V_E_f0_j.end(), v_ec0_1));

            auto iter10(std::find(V_E_f0_j.begin(), V_E_f0_j.end(), v_ec1_0));
            auto iter11(std::find(V_E_f0_j.begin(), V_E_f0_j.end(), v_ec1_1));

            if(iter00 != V_E_f0_j.end() && iter01 != V_E_f0_j.end())
            {
              if(y.at(V_E_f0_j.at(0)) < y.at(V_E_f0_j.at(1)))
              {
                m.get_topologies().at(ipi_edge_vertex).at(E_f0_j).at(0) = V_E_f0_j.at(1);
                m.get_topologies().at(ipi_edge_vertex).at(E_f0_j).at(1) = V_E_f0_j.at(0);
              }
              edges_processed.push_back(E_f0_j);
              edge_types.push_back(et_c_x);
            }

            if(iter10 != V_E_f0_j.end() && iter11 != V_E_f0_j.end())
            {
              if(y.at(V_E_f0_j.at(0)) < y.at(V_E_f0_j.at(1)))
              {
                m.get_topologies().at(ipi_edge_vertex).at(E_f0_j).at(0) = V_E_f0_j.at(1);
                m.get_topologies().at(ipi_edge_vertex).at(E_f0_j).at(1) = V_E_f0_j.at(0);
              }
              edges_processed.push_back(E_f0_j);
              edge_types.push_back(et_c_x);
            }
          }
          ///set iz-curve
          auto V_ez0(m.get_adjacent_polytopes(pl_edge, pl_vertex, ez0));
          auto V_ez1(m.get_adjacent_polytopes(pl_edge, pl_vertex, ez1));

          m.get_topologies().at(ipi_face_vertex).at(0).at(0) = V_ez0.at(0);
          m.get_topologies().at(ipi_face_vertex).at(0).at(1) = V_ez0.at(1);
          m.get_topologies().at(ipi_face_vertex).at(0).at(2) = V_ez1.at(0);
          m.get_topologies().at(ipi_face_vertex).at(0).at(3) = V_ez1.at(1);
        }

        faces_processed.push_back(0);

        ///start on all adjacent faces
        auto E_fi(m.get_adjacent_polytopes(pl_face, pl_edge, 0));
        decltype(E_fi) F_E_fi;
        for(auto E_fi_j : E_fi)
        {
          auto F_E_fi_j(m.get_adjacent_polytopes(pl_edge, pl_face, E_fi_j));
          std::sort(F_E_fi_j.begin(), F_E_fi_j.end());

          decltype(F_E_fi_j) tmp(F_E_fi.size() + F_E_fi_j.size());
          auto iter(std::set_union(F_E_fi.begin(), F_E_fi.end(), F_E_fi_j.begin(), F_E_fi_j.end(), tmp.begin()));
          tmp.resize(Index(iter - tmp.begin()));
          F_E_fi = tmp;
        }

        for(auto F_E_fi_j : F_E_fi)
          _establish_iz_property(m, x, y, faces_processed, edges_processed, edge_types, 0, F_E_fi_j);
      }


      ///for QUAD 3D meshes
      template<typename TopologyType_,
               template <typename, typename> class OuterStorageType_,
               typename AT_>
      static bool iz_property(const Mesh<Dim3D, TopologyType_, OuterStorageType_>& m, const AT_& x, const AT_& y,
                              const AT_& z)
      {
        for(Index i(0) ; i < m.num_polytopes(pl_face) ; ++i)
        {
          auto v_fi(m.get_adjacent_polytopes(pl_face, pl_vertex, i));

          if(v_fi.size() != 4)
          {
            std::cout << "WARNING: not a pure quad mesh!" << std::endl;
          }

          typename AT_::data_type_ e0_x(x.at(v_fi.at(1)) - x.at(v_fi.at(0)));
          typename AT_::data_type_ e0_y(y.at(v_fi.at(1)) - y.at(v_fi.at(0)));
          typename AT_::data_type_ e0_z(z.at(v_fi.at(1)) - z.at(v_fi.at(0)));
          typename AT_::data_type_ ez_x(x.at(v_fi.at(2)) - x.at(v_fi.at(1)));
          typename AT_::data_type_ ez_y(y.at(v_fi.at(2)) - y.at(v_fi.at(1)));
          typename AT_::data_type_ ez_z(z.at(v_fi.at(2)) - z.at(v_fi.at(1)));

          ///check sanity of vertex-based iz-curve (cross-edges)
          if(e0_x > 0. && e0_y > 0.)
          {
            if(!(e0_z > 0. && ez_x < 0.))
            {
              std::cout << "WARNING: malformed cross-edge in iz-curve! e0_x >= 0 AND e0_y >= 0 => ez_x < 0 AND e0_z > 0, but ez_x is " << ez_x << " and e0_z is " << e0_z << "!" << std::endl;
              return false;
            }
          }
          else if(e0_x < 0. && e0_y < 0. )
          {
            if(!(e0_z < 0. && ez_x > 0.))
            {
              std::cout << "WARNING: malformed cross-edge in iz-curve! e0_x < 0 AND e0_y < 0 => ez_x > 0 AND e0_z < 0, but ez_x is " << ez_x << " and e0_z is " << e0_z << "!" << std::endl;
              return false;
            }
          }
          else if(e0_x <= 0. && e0_y > 0. )
          {
            if(!(ez_y < 0. && ez_z < 0.))
            {
              std::cout << "WARNING: malformed cross-edge in iz-curve! e0_x <= 0 AND e0_y > 0 => ez_y < 0 AND ez_z < 0, but ez_x is " << ez_x << " and ez_z is " << ez_z << "!" << std::endl;
              return false;
            }
          }
          else if(e0_x > 0. && e0_y <= 0. )
          {
            if(!(ez_y >= 0. && ez_z >= 0.))
            {
              std::cout << "WARNING: malformed cross-edge in iz-curve! e0_x > 0 AND e0_y <= 0 => ez_y >= 0 AND ez_z >= 0, but ez_y is " << ez_y << " and ez_z is " << ez_z << "!" << std::endl;
              return false;
            }
          }


          ///check existence of iz-edges
          ///check existence of completion edges
          auto e_fi(m.get_adjacent_polytopes(pl_face, pl_edge, i));
          bool found_e0(false);
          bool found_e1(false);
          bool found_c0(false);
          bool found_c1(false);
          for(auto e_fi_j : e_fi)
          {
            auto v_e_fi_j(m.get_adjacent_polytopes(pl_edge, pl_vertex, e_fi_j));
            found_e0 = v_e_fi_j.at(0) == v_fi.at(0) && v_e_fi_j.at(1) == v_fi.at(1) ? true : found_e0;
            found_e1 = v_e_fi_j.at(0) == v_fi.at(2) && v_e_fi_j.at(1) == v_fi.at(3) ? true : found_e1;
            found_c0 = v_e_fi_j.at(0) == v_fi.at(0) && v_e_fi_j.at(1) == v_fi.at(2) ? true : found_c0;
            found_c1 = v_e_fi_j.at(0) == v_fi.at(1) && v_e_fi_j.at(1) == v_fi.at(3) ? true : found_c1;
          }
          if(!(found_e0 && found_e1 && found_c0 && found_c1))
          {
            if(!found_e0)
              std::cout << "WARNING: no matching iz-edge to iz-curve (e0) at face " << i << "!" << std::endl;
            if(!found_e1)
              std::cout << "WARNING: no matching iz-edge to iz-curve (e1) at face " << i << "!" << std::endl;
            if(!found_c0)
              std::cout << "WARNING: no matching completion-edge to iz-curve (c0) at face " << i << "!" << std::endl;
            if(!found_c1)
              std::cout << "WARNING: no matching completion-edge to iz-curve (c1) at face " << i << "!" << std::endl;

            return false;
          }
        }
        return true;
      }

      template<typename TopologyType_,
               template <typename, typename> class OuterStorageType_,
               typename AT_>
      static void establish_iz_property(Mesh<Dim3D, TopologyType_, OuterStorageType_>& m, const AT_& x, const AT_& y,
                                        const AT_& z)
      {
        if(iz_property(m, x, y, z))
          return;

        OuterStorageType_<Index, std::allocator<Index> > faces_processed;
        OuterStorageType_<Index, std::allocator<Index> > edges_processed;
        OuterStorageType_<EdgeTypes, std::allocator<EdgeTypes> > edge_types;

        ///start with face 0
        auto E_f0(m.get_adjacent_polytopes(pl_face, pl_edge, 0));
        ///pick edge 0 and find e1 = E_f0 - E_V_e0
        auto e0 = E_f0.at(0);
        std::sort(E_f0.begin(), E_f0.end());
        auto V_e0(m.get_adjacent_polytopes(pl_edge, pl_vertex, e0));
        decltype(V_e0) E_V_e0;
        for(auto V_e0_j : V_e0)
        {
          auto E_V_e0_j(m.get_adjacent_polytopes(pl_vertex, pl_edge, V_e0_j));
          std::sort(E_V_e0_j.begin(), E_V_e0_j.end());

          decltype(E_V_e0) tmp(E_V_e0.size() + E_V_e0_j.size());
          auto iter(std::set_union(E_V_e0.begin(), E_V_e0.end(), E_V_e0_j.begin(), E_V_e0_j.end(), tmp.begin()));
          tmp.resize(Index(iter - tmp.begin()));
          E_V_e0 = tmp;
        }
        decltype(V_e0) E_f0_minus_E_V_e0(E_f0.size());
        std::set_difference(E_f0.begin(), E_f0.end(), E_V_e0.begin(), E_V_e0.end(), E_f0_minus_E_V_e0.begin());
        auto e1(E_f0_minus_E_V_e0.at(0));
        auto V_e1(m.get_adjacent_polytopes(pl_edge, pl_vertex, e1));


        ///direct e0, e1 to positive x direction (if possible), heuristics: take smallest y-coord-sum-edge for real ez0
        auto x_diff_ez0(x.at(V_e0.at(1)) - x.at(V_e0.at(0)));
        auto x_diff_ez1(x.at(V_e1.at(1)) - x.at(V_e1.at(0)));
        auto y_diff_ez0(y.at(V_e0.at(1)) - y.at(V_e0.at(0)));
        auto y_diff_ez1(y.at(V_e1.at(1)) - y.at(V_e1.at(0)));
        if(x_diff_ez0 != 0 && y_diff_ez0 != y_diff_ez1) //x pos mode
        {
          auto y_sum_e0(y.at(V_e0.at(0)) + y.at(V_e0.at(1)));
          auto y_sum_e1(y.at(V_e1.at(0)) + y.at(V_e1.at(1)));

          auto ez0(e0);
          auto V_ez0(V_e0);
          auto ez1(e1);
          auto V_ez1(V_e1);
          if(y_sum_e0 > y_sum_e1)
          {
            ez0=e1;
            V_ez0=V_e1;
            ez1=e0;
            V_ez1=V_e0;
          }
          edges_processed.push_back(ez0);
          edge_types.push_back(et_iz_x);

          //auto ez1(x_sum_e0 < x_sum_e1 ? e1 : e0);
          edges_processed.push_back(ez1);
          edge_types.push_back(et_iz_x);

          x_diff_ez0 = x.at(V_ez0.at(1)) - x.at(V_ez0.at(0));
          x_diff_ez1 = x.at(V_ez1.at(1)) - x.at(V_ez1.at(0));

          if(x_diff_ez0 < 0)
          {
            m.get_topologies().at(ipi_edge_vertex).at(ez0).at(0) = V_e0.at(1);
            m.get_topologies().at(ipi_edge_vertex).at(ez0).at(1) = V_e0.at(0);
          }

          if(x_diff_ez1 < 0)
          {
            m.get_topologies().at(ipi_edge_vertex).at(ez1).at(0) = V_e1.at(1);
            m.get_topologies().at(ipi_edge_vertex).at(ez1).at(1) = V_e1.at(0);
          }

          ///find completion edges
          auto v_ec0_0(ez0 == e0 ? V_e0.at(0) : V_e1.at(0));
          auto v_ec0_1(ez0 == e0 ? V_e1.at(0) : V_e0.at(0));
          auto v_ec1_0(ez0 == e0 ? V_e0.at(1) : V_e1.at(1));
          auto v_ec1_1(ez0 == e0 ? V_e1.at(1) : V_e0.at(1));

          for(auto E_f0_j : E_f0)
          {
            auto V_E_f0_j(m.get_adjacent_polytopes(pl_edge, pl_vertex, E_f0_j));

            auto iter00(std::find(V_E_f0_j.begin(), V_E_f0_j.end(), v_ec0_0));
            auto iter01(std::find(V_E_f0_j.begin(), V_E_f0_j.end(), v_ec0_1));

            auto iter10(std::find(V_E_f0_j.begin(), V_E_f0_j.end(), v_ec1_0));
            auto iter11(std::find(V_E_f0_j.begin(), V_E_f0_j.end(), v_ec1_1));

            if(iter00 != V_E_f0_j.end() && iter01 != V_E_f0_j.end())
            {
              if(x.at(V_E_f0_j.at(0)) < x.at(V_E_f0_j.at(1)))//hier ggf. z
              {
                m.get_topologies().at(ipi_edge_vertex).at(E_f0_j).at(0) = V_E_f0_j.at(1);
                m.get_topologies().at(ipi_edge_vertex).at(E_f0_j).at(1) = V_E_f0_j.at(0);
              }

              edges_processed.push_back(E_f0_j);
              edge_types.push_back(et_c_y);
            }

            if(iter10 != V_E_f0_j.end() && iter11 != V_E_f0_j.end())
            {
              if(x.at(V_E_f0_j.at(0)) < x.at(V_E_f0_j.at(1)))//hier ggf. z
              {
                m.get_topologies().at(ipi_edge_vertex).at(E_f0_j).at(0) = V_E_f0_j.at(1);
                m.get_topologies().at(ipi_edge_vertex).at(E_f0_j).at(1) = V_E_f0_j.at(0);
              }
              edges_processed.push_back(E_f0_j);
              edge_types.push_back(et_c_y);
            }
          }

          ///set iz-curve
          //auto V_ez0(m.get_adjacent_polytopes(pl_edge, pl_vertex, ez0));
          //auto V_ez1(m.get_adjacent_polytopes(pl_edge, pl_vertex, ez1));

          m.get_topologies().at(ipi_face_vertex).at(0).at(0) = V_ez0.at(0);
          m.get_topologies().at(ipi_face_vertex).at(0).at(1) = V_ez0.at(1);
          m.get_topologies().at(ipi_face_vertex).at(0).at(2) = V_ez1.at(0);
          m.get_topologies().at(ipi_face_vertex).at(0).at(3) = V_ez1.at(1);
        }
        else if(y_diff_ez0 != 0 )//y pos mode
        {
          auto x_sum_e0(x.at(V_e0.at(0)) + x.at(V_e0.at(1)));
          auto x_sum_e1(x.at(V_e1.at(0)) + x.at(V_e1.at(1)));
          //auto ez0(x_sum_e0 < x_sum_e1 ? e0 : e1);

          auto ez0(e0);
          auto V_ez0(V_e0);
          auto ez1(e1);
          auto V_ez1(V_e1);
          if(x_sum_e0 > x_sum_e1)
          {
            ez0=e1;
            V_ez0=V_e1;
            ez1=e0;
            V_ez1=V_e0;
          }
          edges_processed.push_back(ez0);
          edge_types.push_back(et_iz_y);

          //auto ez1(x_sum_e0 < x_sum_e1 ? e1 : e0);
          edges_processed.push_back(ez1);
          edge_types.push_back(et_iz_y);

          y_diff_ez0 = y.at(V_ez0.at(1)) - y.at(V_ez0.at(0));
          y_diff_ez1 = y.at(V_ez1.at(1)) - y.at(V_ez1.at(0));

          if(y_diff_ez0 > 0)
          {
            m.get_topologies().at(ipi_edge_vertex).at(ez0).at(0) = V_e0.at(1);
            m.get_topologies().at(ipi_edge_vertex).at(ez0).at(1) = V_e0.at(0);
          }

          if(y_diff_ez1 > 0)
          {
            m.get_topologies().at(ipi_edge_vertex).at(ez1).at(0) = V_e1.at(1);
            m.get_topologies().at(ipi_edge_vertex).at(ez1).at(1) = V_e1.at(0);
          }

          ///find completion edges
          auto v_ec0_0(ez0 == e0 ? V_e0.at(0) : V_e1.at(0));
          auto v_ec0_1(ez0 == e0 ? V_e1.at(0) : V_e0.at(0));
          auto v_ec1_0(ez0 == e0 ? V_e0.at(1) : V_e1.at(1));
          auto v_ec1_1(ez0 == e0 ? V_e1.at(1) : V_e0.at(1));

          for(auto E_f0_j : E_f0)
          {
            auto V_E_f0_j(m.get_adjacent_polytopes(pl_edge, pl_vertex, E_f0_j));

            auto iter00(std::find(V_E_f0_j.begin(), V_E_f0_j.end(), v_ec0_0));
            auto iter01(std::find(V_E_f0_j.begin(), V_E_f0_j.end(), v_ec0_1));

            auto iter10(std::find(V_E_f0_j.begin(), V_E_f0_j.end(), v_ec1_0));
            auto iter11(std::find(V_E_f0_j.begin(), V_E_f0_j.end(), v_ec1_1));

            if(iter00 != V_E_f0_j.end() && iter01 != V_E_f0_j.end())
            {
              if(y.at(V_E_f0_j.at(0)) > y.at(V_E_f0_j.at(1)))
              {
                m.get_topologies().at(ipi_edge_vertex).at(E_f0_j).at(0) = V_E_f0_j.at(1);
                m.get_topologies().at(ipi_edge_vertex).at(E_f0_j).at(1) = V_E_f0_j.at(0);
              }
              edges_processed.push_back(E_f0_j);
              edge_types.push_back(et_c_x);
            }

            if(iter10 != V_E_f0_j.end() && iter11 != V_E_f0_j.end())
            {
              if(y.at(V_E_f0_j.at(0)) > y.at(V_E_f0_j.at(1)))
              {
                m.get_topologies().at(ipi_edge_vertex).at(E_f0_j).at(0) = V_E_f0_j.at(1);
                m.get_topologies().at(ipi_edge_vertex).at(E_f0_j).at(1) = V_E_f0_j.at(0);
              }
              edges_processed.push_back(E_f0_j);
              edge_types.push_back(et_c_x);
            }
          }
          ///set iz-curve
          //auto V_ez0(m.get_adjacent_polytopes(pl_edge, pl_vertex, ez0));
          //auto V_ez1(m.get_adjacent_polytopes(pl_edge, pl_vertex, ez1));

          m.get_topologies().at(ipi_face_vertex).at(0).at(0) = V_ez0.at(0);
          m.get_topologies().at(ipi_face_vertex).at(0).at(1) = V_ez0.at(1);
          m.get_topologies().at(ipi_face_vertex).at(0).at(2) = V_ez1.at(0);
          m.get_topologies().at(ipi_face_vertex).at(0).at(3) = V_ez1.at(1);
        }
        else // z-pos mode
        {
          auto y_sum_e0(y.at(V_e0.at(0)) + y.at(V_e0.at(1)));
          auto y_sum_e1(y.at(V_e1.at(0)) + y.at(V_e1.at(1)));

          auto ez0(e0);
          auto V_ez0(V_e0);
          auto ez1(e1);
          auto V_ez1(V_e1);
          if(y_sum_e0 > y_sum_e1)
          {
            ez0=e1;
            V_ez0=V_e1;
            ez1=e0;
            V_ez1=V_e0;
          }
          edges_processed.push_back(ez0);
          edge_types.push_back(et_iz_x);

          //auto ez1(x_sum_e0 < x_sum_e1 ? e1 : e0);
          edges_processed.push_back(ez1);
          edge_types.push_back(et_iz_x);

          x_diff_ez0 = x.at(V_ez0.at(1)) - x.at(V_ez0.at(0));
          x_diff_ez1 = x.at(V_ez1.at(1)) - x.at(V_ez1.at(0));

          if(x_diff_ez0 < 0)
          {
            m.get_topologies().at(ipi_edge_vertex).at(ez0).at(0) = V_e0.at(1);
            m.get_topologies().at(ipi_edge_vertex).at(ez0).at(1) = V_e0.at(0);
          }

          if(x_diff_ez1 < 0)
          {
            m.get_topologies().at(ipi_edge_vertex).at(ez1).at(0) = V_e1.at(1);
            m.get_topologies().at(ipi_edge_vertex).at(ez1).at(1) = V_e1.at(0);
          }

          ///find completion edges
          auto v_ec0_0(ez0 == e0 ? V_e0.at(0) : V_e1.at(0));
          auto v_ec0_1(ez0 == e0 ? V_e1.at(0) : V_e0.at(0));
          auto v_ec1_0(ez0 == e0 ? V_e0.at(1) : V_e1.at(1));
          auto v_ec1_1(ez0 == e0 ? V_e1.at(1) : V_e0.at(1));

          for(auto E_f0_j : E_f0)
          {
            auto V_E_f0_j(m.get_adjacent_polytopes(pl_edge, pl_vertex, E_f0_j));

            auto iter00(std::find(V_E_f0_j.begin(), V_E_f0_j.end(), v_ec0_0));
            auto iter01(std::find(V_E_f0_j.begin(), V_E_f0_j.end(), v_ec0_1));

            auto iter10(std::find(V_E_f0_j.begin(), V_E_f0_j.end(), v_ec1_0));
            auto iter11(std::find(V_E_f0_j.begin(), V_E_f0_j.end(), v_ec1_1));

            if(iter00 != V_E_f0_j.end() && iter01 != V_E_f0_j.end())
            {
              if(x.at(V_E_f0_j.at(0)) < x.at(V_E_f0_j.at(1)))//hier ggf. z
              {
                m.get_topologies().at(ipi_edge_vertex).at(E_f0_j).at(0) = V_E_f0_j.at(1);
                m.get_topologies().at(ipi_edge_vertex).at(E_f0_j).at(1) = V_E_f0_j.at(0);
              }

              edges_processed.push_back(E_f0_j);
              edge_types.push_back(et_c_y);
            }

            if(iter10 != V_E_f0_j.end() && iter11 != V_E_f0_j.end())
            {
              if(x.at(V_E_f0_j.at(0)) < x.at(V_E_f0_j.at(1)))//hier ggf. z
              {
                m.get_topologies().at(ipi_edge_vertex).at(E_f0_j).at(0) = V_E_f0_j.at(1);
                m.get_topologies().at(ipi_edge_vertex).at(E_f0_j).at(1) = V_E_f0_j.at(0);
              }
              edges_processed.push_back(E_f0_j);
              edge_types.push_back(et_c_y);
            }
          }

          ///set iz-curve
          //auto V_ez0(m.get_adjacent_polytopes(pl_edge, pl_vertex, ez0));
          //auto V_ez1(m.get_adjacent_polytopes(pl_edge, pl_vertex, ez1));

          m.get_topologies().at(ipi_face_vertex).at(0).at(0) = V_ez0.at(0);
          m.get_topologies().at(ipi_face_vertex).at(0).at(1) = V_ez0.at(1);
          m.get_topologies().at(ipi_face_vertex).at(0).at(2) = V_ez1.at(0);
          m.get_topologies().at(ipi_face_vertex).at(0).at(3) = V_ez1.at(1);
        }

        faces_processed.push_back(0);

        ///start on all adjacent faces
        auto E_fi(m.get_adjacent_polytopes(pl_face, pl_edge, 0));
        decltype(E_fi) F_E_fi;
        for(auto E_fi_j : E_fi)
        {
          auto F_E_fi_j(m.get_adjacent_polytopes(pl_edge, pl_face, E_fi_j));
          std::sort(F_E_fi_j.begin(), F_E_fi_j.end());

          decltype(F_E_fi_j) tmp(F_E_fi.size() + F_E_fi_j.size());
          auto iter(std::set_union(F_E_fi.begin(), F_E_fi.end(), F_E_fi_j.begin(), F_E_fi_j.end(), tmp.begin()));
          tmp.resize(Index(iter - tmp.begin()));
          F_E_fi = tmp;
        }

        for(auto F_E_fi_j : F_E_fi)
          _establish_iz_property(m, x, y, z, faces_processed, edges_processed, edge_types, 0, F_E_fi_j);auto test(m.get_adjacent_polytopes(pl_face,pl_vertex,0));
      }


      private:
        template<typename TopologyType_,
               template <typename, typename> class OuterStorageType_,
               typename AT_>
        static void _establish_iz_property(Mesh<Dim2D, TopologyType_, OuterStorageType_>& m, const AT_& x, const AT_& y,
                                          OuterStorageType_<Index, std::allocator<Index> >& fp,
                                          OuterStorageType_<Index, std::allocator<Index> >& ep,
                                          OuterStorageType_<EdgeTypes, std::allocator<EdgeTypes> >& et,
                                          Index face_from,
                                          Index face_num)
        {
          ///face already processed -> recursion end
          if(std::find(fp.begin(), fp.end(), face_num) != fp.end())
            return;

          ///retrieve information about already processed edges
          auto E_fi(m.get_adjacent_polytopes(pl_face, pl_edge, face_num));
          OuterStorageType_<Index, std::allocator<Index> > local_ep;
          OuterStorageType_<EdgeTypes, std::allocator<EdgeTypes> > local_et;

          for(auto E_fi_j : E_fi)
          {
            auto iter(std::find(ep.begin(), ep.end(), E_fi_j));

            if(iter != ep.end())
            {
              local_et.push_back(et.at(Index(iter - ep.begin())));
              local_ep.push_back(ep.at(Index(iter - ep.begin())));
            }
          }

          ///all edges processed -> recursion end
          if(local_et.size() >= 4)
          {
            fp.push_back(face_num);
            return;
          }

          ///pull directions from origin face
          auto E0(m.get_comm_intersection(pl_face, pl_edge, face_from, face_num));

          if(E0.size() == 0) ///coming from diagonal face
            return;

          auto e0(E0.at(0));

          auto V_e0(m.get_adjacent_polytopes(pl_edge, pl_vertex, e0));
          decltype(V_e0) E_V_e0;
          for(auto V_e0_i : V_e0)
          {
            auto E_V_e0_i(m.get_adjacent_polytopes(pl_vertex, pl_edge, V_e0_i));
            E_V_e0 = STLUtil::set_union(E_V_e0, E_V_e0_i);
          }
          auto E_f1(m.get_adjacent_polytopes(pl_face, pl_edge, face_num));
          auto E1(STLUtil::set_difference(E_f1, E_V_e0));
          auto e1(E1.at(0));

          auto iter0(std::find(ep.begin(), ep.end(), e0));
          Index access0(Index(iter0 - ep.begin()));
          _direct(m, et.at(access0) == et_iz_x || et.at(access0) == et_c_x ? x : y, e0, e1); //capture x_diff == 0 case
          ep.push_back(e1);
          et.push_back(et.at(access0));

          //et1 = E_f0 - E_V_e0
          auto Ef0(m.get_adjacent_polytopes(pl_face, pl_edge, face_from));
          auto et1(STLUtil::set_difference(Ef0, E_V_e0).at(0));

          //et2 = any of E_f0 - {e0, et1}
          decltype(Ef0) e0et1;
          e0et1.push_back(e0);
          e0et1.push_back(et1);
          auto et2(STLUtil::set_difference(Ef0, e0et1).at(0));

          //{e2, e3} = E_f1 - {e0, e1}
          auto E23(STLUtil::set_difference(E_f1, E0));
          E23 = STLUtil::set_difference(E23, E1);
          auto e2(E23.at(0));
          auto e3(E23.at(1));

          //_direct(et2, {e2, e3})
          auto iter1(std::find(ep.begin(), ep.end(), et2));
          Index access1(Index(iter1 - ep.begin()));
          _direct(m, et.at(access1) == et_iz_x || et.at(access1) == et_c_x ? x : y, et2, e2); //capture x_diff == 0 case
          ep.push_back(e2);
          et.push_back(et.at(access1));
          _direct(m, et.at(access1) == et_iz_x || et.at(access1) == et_c_x ? x : y, et2, e3); //capture x_diff == 0 case
          ep.push_back(e3);
          et.push_back(et.at(access1));

          if(et.at(access0) == et_iz_x || et.at(access0) == et_iz_y) //e0, e1 is iz curve
          {
            auto V_e0_new(m.get_adjacent_polytopes(pl_edge, pl_vertex, e0));
            auto V_e1_new(m.get_adjacent_polytopes(pl_edge, pl_vertex, e1));
            if(et.at(access0) == et_iz_x)
            {
              auto y_sum_e0(y.at(V_e0_new.at(0)) + y.at(V_e0_new.at(1)));
              auto y_sum_e1(y.at(V_e1_new.at(0)) + y.at(V_e1_new.at(1)));

              if(y_sum_e0 < y_sum_e1)
              {
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(0) = V_e0_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(1) = V_e0_new.at(1);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(2) = V_e1_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(3) = V_e1_new.at(1);
              }
              else
              {
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(0) = V_e1_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(1) = V_e1_new.at(1);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(2) = V_e0_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(3) = V_e0_new.at(1);
              }
            }
            else
            {
              auto x_sum_e0(x.at(V_e0_new.at(0)) + x.at(V_e0_new.at(1)));
              auto x_sum_e1(x.at(V_e1_new.at(0)) + x.at(V_e1_new.at(1)));

              if(x_sum_e0 > x_sum_e1)
              {
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(0) = V_e0_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(1) = V_e0_new.at(1);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(2) = V_e1_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(3) = V_e1_new.at(1);
              }
              else
              {
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(0) = V_e1_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(1) = V_e1_new.at(1);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(2) = V_e0_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(3) = V_e0_new.at(1);
              }
            }
          }
          else if(et.at(access1) == et_iz_x || et.at(access1) == et_iz_y) //e2, e3 is iz curve
          {
            auto V_e2_new(m.get_adjacent_polytopes(pl_edge, pl_vertex, e2));
            auto V_e3_new(m.get_adjacent_polytopes(pl_edge, pl_vertex, e3));
            if(et.at(access1) == et_iz_x)
            {
              auto y_sum_e2(y.at(V_e2_new.at(0)) + y.at(V_e2_new.at(1)));
              auto y_sum_e3(y.at(V_e3_new.at(0)) + y.at(V_e3_new.at(1)));

              if(y_sum_e2 < y_sum_e3)
              {
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(0) = V_e2_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(1) = V_e2_new.at(1);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(2) = V_e3_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(3) = V_e3_new.at(1);
              }
              else
              {
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(0) = V_e3_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(1) = V_e3_new.at(1);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(2) = V_e2_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(3) = V_e2_new.at(1);
              }
            }
            else
            {
              auto x_sum_e2(x.at(V_e2_new.at(0)) + x.at(V_e2_new.at(1)));
              auto x_sum_e3(x.at(V_e3_new.at(0)) + x.at(V_e3_new.at(1)));

              if(x_sum_e2 > x_sum_e3)
              {
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(0) = V_e2_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(1) = V_e2_new.at(1);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(2) = V_e3_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(3) = V_e3_new.at(1);
              }
              else
              {
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(0) = V_e3_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(1) = V_e3_new.at(1);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(2) = V_e2_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(3) = V_e2_new.at(1);
              }
            }
          }

          //end and recursion
          decltype(E_fi) F_E_fi;
          for(auto E_fi_j : E_fi)
          {
            auto F_E_fi_j(m.get_adjacent_polytopes(pl_edge, pl_face, E_fi_j));
            std::sort(F_E_fi_j.begin(), F_E_fi_j.end());

            decltype(F_E_fi_j) tmp(F_E_fi.size() + F_E_fi_j.size());
            auto iter3(std::set_union(F_E_fi.begin(), F_E_fi.end(), F_E_fi_j.begin(), F_E_fi_j.end(), tmp.begin()));
            tmp.resize(Index(iter3 - tmp.begin()));
            F_E_fi = tmp;
          }
          fp.push_back(face_num);
          ///do recursion
          for(auto F_E_fi_j : F_E_fi)
            _establish_iz_property(m, x, y, fp, ep, et, face_num, F_E_fi_j);
        }

        template<typename TopologyType_,
                 template <typename, typename> class OuterStorageType_,
                 typename AT_>
        static void _direct(Mesh<Dim2D, TopologyType_, OuterStorageType_>& m, const AT_& xy, Index e0, Index e1)
        {
          auto V_e0(m.get_adjacent_polytopes(pl_edge, pl_vertex, e0));
          auto xy_diff_e0(xy.at(V_e0.at(1)) - xy.at(V_e0.at(0)));

          if(xy_diff_e0 == typename AT_::data_type_(0))
            throw MeshError("Edge cannot be directed like this!");

          auto V_e1(m.get_adjacent_polytopes(pl_edge, pl_vertex, e1));
          auto xy_diff_e1(xy.at(V_e1.at(1)) - xy.at(V_e1.at(0)));

          if((xy_diff_e0 > 0 && xy_diff_e1 > 0) || (xy_diff_e0 < 0 && xy_diff_e1 < 0))
            return;
          else
          {
            m.get_topologies().at(ipi_edge_vertex).at(e1).at(0) = V_e1.at(1);
            m.get_topologies().at(ipi_edge_vertex).at(e1).at(1) = V_e1.at(0);
          }
        }

        template<typename TopologyType_,
               template <typename, typename> class OuterStorageType_,
               typename AT_>
        static void _establish_iz_property(Mesh<Dim3D, TopologyType_, OuterStorageType_>& m, const AT_& x, const AT_& y,
                                          const AT_& z,
                                          OuterStorageType_<Index, std::allocator<Index> >& fp,
                                          OuterStorageType_<Index, std::allocator<Index> >& ep,
                                          OuterStorageType_<EdgeTypes, std::allocator<EdgeTypes> >& et,
                                          Index face_from,
                                          Index face_num)
        {
          ///face already processed -> recursion end
          if(std::find(fp.begin(), fp.end(), face_num) != fp.end())
            return;

          ///retrieve information about already processed edges
          auto E_fi(m.get_adjacent_polytopes(pl_face, pl_edge, face_num));
          OuterStorageType_<Index, std::allocator<Index> > local_ep;
          OuterStorageType_<EdgeTypes, std::allocator<EdgeTypes> > local_et;

          for(auto E_fi_j : E_fi)
          {
            auto iter(std::find(ep.begin(), ep.end(), E_fi_j));

            if(iter != ep.end())
            {
              local_et.push_back(et.at(Index(iter - ep.begin())));
              local_ep.push_back(ep.at(Index(iter - ep.begin())));
            }
          }

          ///all edges processed -> recursion end
          if(local_et.size() >= 4)
          {
            fp.push_back(face_num);
            return;
          }

          ///pull directions from origin face
          auto E0(m.get_comm_intersection(pl_face, pl_edge, face_from, face_num));

          if(E0.size() == 0) ///coming from diagonal face
            return;

          auto e0(E0.at(0));

          auto V_e0(m.get_adjacent_polytopes(pl_edge, pl_vertex, e0));
          decltype(V_e0) E_V_e0;
          for(auto V_e0_i : V_e0)
          {
            auto E_V_e0_i(m.get_adjacent_polytopes(pl_vertex, pl_edge, V_e0_i));
            E_V_e0 = STLUtil::set_union(E_V_e0, E_V_e0_i);
          }
          auto E_f1(m.get_adjacent_polytopes(pl_face, pl_edge, face_num));
          auto E1(STLUtil::set_difference(E_f1, E_V_e0));
          auto e1(E1.at(0));

          auto iter0(std::find(ep.begin(), ep.end(), e0));
          Index access0(Index(iter0 - ep.begin()));
          _direct(m, et.at(access0) == et_iz_x || et.at(access0) == et_c_x ? x : y, e0, e1); //capture x_diff == 0 case
          ep.push_back(e1);
          et.push_back(et.at(access0));

          //et1 = E_f0 - E_V_e0
          auto Ef0(m.get_adjacent_polytopes(pl_face, pl_edge, face_from));
          auto et1(STLUtil::set_difference(Ef0, E_V_e0).at(0));

          //et2 = any of E_f0 - {e0, et1}
          decltype(Ef0) e0et1;
          e0et1.push_back(e0);
          e0et1.push_back(et1);
          auto et2(STLUtil::set_difference(Ef0, e0et1).at(0));

          //{e2, e3} = E_f1 - {e0, e1}
          auto E23(STLUtil::set_difference(E_f1, E0));
          E23 = STLUtil::set_difference(E23, E1);
          auto e2(E23.at(0));
          auto e3(E23.at(1));

          //_direct(et2, {e2, e3})
          auto iter1(std::find(ep.begin(), ep.end(), et2));
          Index access1(Index(iter1 - ep.begin()));
          _direct(m, et.at(access1) == et_iz_x || et.at(access1) == et_c_x ? x : y, et2, e2); //capture x_diff == 0 case
          ep.push_back(e2);
          et.push_back(et.at(access1));
          _direct(m, et.at(access1) == et_iz_x || et.at(access1) == et_c_x ? x : y, et2, e3); //capture x_diff == 0 case
          ep.push_back(e3);
          et.push_back(et.at(access1));

          if(et.at(access0) == et_iz_x || et.at(access0) == et_iz_y) //e0, e1 is iz curve
          {
            auto V_e0_new(m.get_adjacent_polytopes(pl_edge, pl_vertex, e0));
            auto V_e1_new(m.get_adjacent_polytopes(pl_edge, pl_vertex, e1));
            if(et.at(access0) == et_iz_x)
            {
              auto y_sum_e0(y.at(V_e0_new.at(0)) + y.at(V_e0_new.at(1)));
              auto y_sum_e1(y.at(V_e1_new.at(0)) + y.at(V_e1_new.at(1)));

              if(y_sum_e0 < y_sum_e1)
              {
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(0) = V_e0_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(1) = V_e0_new.at(1);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(2) = V_e1_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(3) = V_e1_new.at(1);
              }
              else
              {
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(0) = V_e1_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(1) = V_e1_new.at(1);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(2) = V_e0_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(3) = V_e0_new.at(1);
              }
            }
            else
            {
              auto x_sum_e0(x.at(V_e0_new.at(0)) + x.at(V_e0_new.at(1)));
              auto x_sum_e1(x.at(V_e1_new.at(0)) + x.at(V_e1_new.at(1)));

              if(x_sum_e0 > x_sum_e1)
              {
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(0) = V_e0_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(1) = V_e0_new.at(1);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(2) = V_e1_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(3) = V_e1_new.at(1);
              }
              else
              {
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(0) = V_e1_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(1) = V_e1_new.at(1);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(2) = V_e0_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(3) = V_e0_new.at(1);
              }
            }
          }
          else if(et.at(access1) == et_iz_x || et.at(access1) == et_iz_y) //e2, e3 is iz curve
          {
            auto V_e2_new(m.get_adjacent_polytopes(pl_edge, pl_vertex, e2));
            auto V_e3_new(m.get_adjacent_polytopes(pl_edge, pl_vertex, e3));
            if(et.at(access1) == et_iz_x)
            {
              auto y_sum_e2(y.at(V_e2_new.at(0)) + y.at(V_e2_new.at(1)));
              auto y_sum_e3(y.at(V_e3_new.at(0)) + y.at(V_e3_new.at(1)));

              if(y_sum_e2 < y_sum_e3)
              {
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(0) = V_e2_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(1) = V_e2_new.at(1);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(2) = V_e3_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(3) = V_e3_new.at(1);
              }
              else
              {
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(0) = V_e3_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(1) = V_e3_new.at(1);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(2) = V_e2_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(3) = V_e2_new.at(1);
              }
            }
            else
            {
              auto x_sum_e2(x.at(V_e2_new.at(0)) + x.at(V_e2_new.at(1)));
              auto x_sum_e3(x.at(V_e3_new.at(0)) + x.at(V_e3_new.at(1)));

              if(x_sum_e2 > x_sum_e3)
              {
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(0) = V_e2_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(1) = V_e2_new.at(1);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(2) = V_e3_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(3) = V_e3_new.at(1);
              }
              else
              {
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(0) = V_e3_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(1) = V_e3_new.at(1);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(2) = V_e2_new.at(0);
                m.get_topologies().at(ipi_face_vertex).at(face_num).at(3) = V_e2_new.at(1);
              }
            }
          }

          //end and recursion
          decltype(E_fi) F_E_fi;
          for(auto E_fi_j : E_fi)
          {
            auto F_E_fi_j(m.get_adjacent_polytopes(pl_edge, pl_face, E_fi_j));
            std::sort(F_E_fi_j.begin(), F_E_fi_j.end());

            decltype(F_E_fi_j) tmp(F_E_fi.size() + F_E_fi_j.size());
            auto iter3(std::set_union(F_E_fi.begin(), F_E_fi.end(), F_E_fi_j.begin(), F_E_fi_j.end(), tmp.begin()));
            tmp.resize(Index(iter3 - tmp.begin()));
            F_E_fi = tmp;
          }
          fp.push_back(face_num);
          ///do recursion
          for(auto F_E_fi_j : F_E_fi)
            _establish_iz_property(m, x, y, z, fp, ep, et, face_num, F_E_fi_j);
        }

        template<typename TopologyType_,
                 template <typename, typename> class OuterStorageType_,
                 typename AT_>
        static void _direct(Mesh<Dim3D, TopologyType_, OuterStorageType_>& m, const AT_& xy, Index e0, Index e1)
        {
          auto V_e0(m.get_adjacent_polytopes(pl_edge, pl_vertex, e0));
          auto xy_diff_e0(xy.at(V_e0.at(1)) - xy.at(V_e0.at(0)));

          if(xy_diff_e0 == typename AT_::data_type_(0))
            throw MeshError("Edge cannot be directed like this!");

          auto V_e1(m.get_adjacent_polytopes(pl_edge, pl_vertex, e1));
          auto xy_diff_e1(xy.at(V_e1.at(1)) - xy.at(V_e1.at(0)));

          if((xy_diff_e0 > 0 && xy_diff_e1 > 0) || (xy_diff_e0 < 0 && xy_diff_e1 < 0))
            return;
          else
          {
            m.get_topologies().at(ipi_edge_vertex).at(e1).at(0) = V_e1.at(1);
            m.get_topologies().at(ipi_edge_vertex).at(e1).at(1) = V_e1.at(0);
          }
        }


    };

  }
}
#endif
