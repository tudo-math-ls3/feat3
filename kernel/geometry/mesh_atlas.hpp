#pragma once
#ifndef KERNEL_GEOMETRY_MESH_ATLAS_HPP
#define KERNEL_GEOMETRY_MESH_ATLAS_HPP 1

// includes, FEAST
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/atlas/chart.hpp>

// includes, system
#include <map>
#include <deque>

namespace FEAST
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

    }; // class MeshAltas<...>
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_MESH_ATLAS_HPP
