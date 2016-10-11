#pragma once
#ifndef KERNEL_GEOMETRY_MESH_ATLAS_HPP
#define KERNEL_GEOMETRY_MESH_ATLAS_HPP 1

// includes, FEAT
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/atlas/chart.hpp>

// includes, system
#include <map>
#include <deque>

namespace FEAT
{
  namespace Geometry
  {
    /**
     * \brief Mesh Atlas class template
     *
     * \tparam Mesh_
     * The type of the mesh that is to be parameterised by this atlas, i.e. ConformalMesh
     *
     * \author Peter Zajac
     */
    template<typename Mesh_>
    class MeshAtlas
    {
    public:
      /// our mesh type
      typedef Mesh_ MeshType;
      /// our chart base-class type
      typedef Atlas::ChartBase<MeshType> MeshChartType;
      /// our vertex type (aka world point type)
      typedef typename MeshType::VertexType VertexType;

    protected:
      /// our chart map type
      typedef std::map<String, MeshChartType*> MeshChartMap;

      /// our chart map
      MeshChartMap _chart_map;

    public:
      /// default CTOR
      MeshAtlas()
      {
      }

      /// move CTOR
      MeshAtlas(MeshAtlas&& other) :
        _chart_map(std::move(other._chart_map))
      {
      }

      /// virtual DTOR
      virtual ~MeshAtlas()
      {
        // delete all charts
        for(auto it = _chart_map.begin(); it != _chart_map.end(); ++it)
        {
          delete (*it).second;
        }
        _chart_map.clear();
      }

      /// \returns The size of dynamically allocated memory in bytes.
      std::size_t bytes() const
      {
        std::size_t s(0);
        for(auto it = _chart_map.begin(); it != _chart_map.end(); ++it)
          s += it->second->bytes();
        return s;
      }

      /**
       * \brief Inserts a new chart into the map.
       *
       * \param[in] name
       * The name of the chart to be inserted.
       *
       * \param[in] chart
       * The mesh chart to be inserted.
       *
       * \param[in] replace
       * Specifies whether to replace an old chart with the same name if it already exists.
       *
       * \returns
       * True if the MeshChart was successfully inserted, meaning a chart of that name was either not present
       * yet, or it was present and replace was set to true.
       */
      bool add_mesh_chart(const String& name, MeshChartType* chart, bool replace = false)
      {
        // try to insert
        auto it = _chart_map.insert(std::make_pair(name, chart));

        // insertion successful?
        if(it.second)
          return true;

        // a chart with that name already exists. replace?
        if(!replace)
          return false;

        // delete old chart
        if((it.first->second) != nullptr)
          delete (it.first->second);

        // assign new chart
        it.first->second = chart;
        return true;
      }

      /**
       * \brief Returns the names of all charts of this atlas.
       */
      std::deque<String> get_chart_names() const
      {
        std::deque<String> names;
        for(auto it = _chart_map.begin(); it != _chart_map.end(); ++it)
        {
          names.push_back((*it).first);
        }
        return names;
      }

      /**
       * \brief Returns a const reference to the chart map
       *
       * \returns A const reference to the chart map
       */
      const MeshChartMap& get_mesh_chart_map() const
      {
        return _chart_map;
      }

      /**
       * \brief Searches for a mesh chart.
       *
       * \param[in] name
       * The name of the mesh chart to the found.
       *
       * \returns
       * A pointer to the mesh chart associated with \p name, or \c nullptr if no chart
       * with that name was found.
       */
      MeshChartType* find_mesh_chart(const String& name)
      {
        auto it = _chart_map.find(name);
        if(it == _chart_map.end())
          return nullptr;
        return (*it).second;
      }

      /** \copydoc find_mesh_chart() */
      const MeshChartType* find_mesh_chart(const String& name) const
      {
        auto it = _chart_map.find(name);
        if(it == _chart_map.end())
          return nullptr;
        return (*it).second;
      }

      /**
       * \brief Applies a "proper rigid" transformation onto the atlas.
       *
       * Let \e v denote the \p origin world point, \e w the \p offset world point and \e R
       * the rotation matrix corresponding to the \p angles, then this function applies the
       * following transformation for any chart point \e x of the atlas:
       *
       *   \f[ x \mapsto w + R\cdot (x - v) \f]
       *
       * \param[in] origin
       * The origin of the transformation. This is subtracted from any vertex before applying the
       * rotation.
       *
       * \param[in] angles
       * The angles of the rotation matrix.
       * - 2D: the rotation angle in radians is stored as:
       *   - angles(0): rotation angle
       *   - angles(1): \e ignored
       * - 3D: the rotation angles in radians stored as:
       *   - angles(0): yaw angle
       *   - angles(1): pitch angle
       *   - angles(2): roll angle
       *
       * \param[in] offset
       * The offset of the transformation. This is added to any vertex after applying the rotation.
       */
      void transform(const VertexType& origin, const VertexType& angles, const VertexType& offset)
      {
        // transform all charts
        for(auto it = _chart_map.begin(); it != _chart_map.end(); ++it)
          it->second->transform(origin, angles, offset);
      }
    }; // class MeshAltas<...>
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_MESH_ATLAS_HPP
