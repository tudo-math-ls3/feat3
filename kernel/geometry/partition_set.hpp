// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_PARTITION_SET_HPP
#define KERNEL_GEOMETRY_PARTITION_SET_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/adjacency/dynamic_graph.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/util/string.hpp>

#include <deque>

namespace FEAT
{
  namespace Geometry
  {
    /**
     * \brief Partition class
     *
     * This class represents a (manual) partitioning of a mesh, which consists of an
     * adjacency graph storing the elements for each rank as well as a few additional
     * entities as the name, the priority and the overlap.
     *
     * \author Peter Zajac
     */
    class Partition
    {
    protected:
      /// the name of this partitioning
      String _name;
      /// the priority for automatic choice
      int _prio;
      /// the level that this partitioning refers to
      int _level;

      /// our actual partitioning
      Adjacency::Graph _patches;

    public:
      /// default construcotr
      Partition() :
        _prio(0),
        _level(0)
      {
      }

      /**
       * \brief Constructor
       *
       * \param[in] graph_
       * The elements-at-rank graph defining the partition.
       *
       * \param[in] name_
       * The name of the partition. May be an empty string
       *
       * \param[in] prio_
       * The priority of the partition.
       *
       * \param[in] level_
       * The refinement level that the partition is defined on.
       */
      explicit Partition(Adjacency::Graph&& graph_, const String& name_, int prio_ = 0, int level_ = 0) :
        _name(name_),
        _prio(prio_),
        _level(level_),
        _patches(std::forward<Adjacency::Graph>(graph_))
      {
      }


      /**
       * \brief Constructor
       *
       * \param[in] graph_
       * The elements-at-rank graph defining the partition.
       *
       * \param[in] name_
       * The name of the partition. May be an empty string
       *
       * \param[in] prio_
       * The priority of the partition.
       *
       * \param[in] level_
       * The refinement level that the partition is defined on.
       */
      explicit Partition(const Adjacency::DynamicGraph& graph_, const String& name_, int prio_ = 0, int level_ = 0) :
        _name(name_),
        _prio(prio_),
        _level(level_),
        _patches(Adjacency::RenderType::as_is, graph_)
      {
      }

      /// move constructor
      Partition(Partition&& other) :
        _name(other._name),
        _prio(other._prio),
        _level(other._level),
        _patches(std::forward<Adjacency::Graph>(other._patches))
      {
      }

      /// move-assign operator
      Partition& operator=(Partition&& other)
      {
        if(this == &other)
          return *this;
        _name = other._name;
        _prio = other._prio;
        _level = other._level;
        _patches = std::forward<Adjacency::Graph>(other._patches);
        return *this;
      }

      /// destructor
      virtual ~Partition()
      {
      }

      /// \returns the number of ranks/patches in the partition.
      Index size() const
      {
        return Index(_patches.get_num_nodes_domain());
      }

      /// \returns the number of ranks/patches in the partition.
      Index get_num_patches() const
      {
        return Index(_patches.get_num_nodes_domain());
      }

      /// \returns the number of elements in the partition.
      Index get_num_elements() const
      {
        return Index(_patches.get_num_nodes_image());
      }

      /// \returns the name of the partition.
      String get_name() const
      {
        return _name;
      }

      /// \returns the priority of the partition.
      int get_priority() const
      {
        return _prio;
      }

      /// \returns the refinement level of the partition.
      int get_level() const
      {
        return _level;
      }

      /// \returns the elements-at-rank adjacency graph of the partition.
      const Adjacency::Graph& get_patches() const
      {
        return _patches;
      }
    }; // class Partition

    /**
     * \brief Partition set class
     *
     * This class acts as a container, which manages a set of partitions.
     * This class is effectively just a wrapper around a deque of Partition objects.
     *
     * \author Peter Zajac
     */
    class PartitionSet
    {
    protected:
      /// the actual partitions in the set
      std::deque<Partition> _parts;

    public:
      PartitionSet()
      {
      }

      PartitionSet(PartitionSet&& other) :
        _parts(std::move(other._parts))
      {
      }

      virtual ~PartitionSet()
      {
      }

      void clear()
      {
        _parts.clear();
      }

      std::deque<Partition>& get_partitions()
      {
        return _parts;
      }

      const std::deque<Partition>& get_partitions() const
      {
        return _parts;
      }

      void add_partition(Partition&& part)
      {
        _parts.emplace_back(std::forward<Partition>(part));
      }

      /**
       * \brief Tries to find a suitable partition in the set.
       *
       * \param[in] size
       * The required size (number of ranks) of the partition. Must be > 0.
       *
       * \param[in] name
       * The required name of the partition. If \p name is an empty string, the name
       * of a partition is considered irrelevant.
       *
       * \param[in] prio
       * The required minimal priority of the partition. Only partitions with a priority greater or
       * equal to \p prio are considered valid candidates.
       *
       * \returns
       * A pointer to a Partition object that meets the requirements or \c nullptr,
       * if no suitable partition was found in the set.
       */
      const Partition* find_partition(int size, const String& name = "", int prio = 0) const
      {
        std::deque<String> names;
        if(!name.empty())
          names.push_back(name);
        return find_partition(size, names, prio);
      }

      const Partition* find_partition(int size, const std::deque<String>& names, int prio = 0) const
      {
        const Partition* part(nullptr);

        // loop over all partitions
        for(auto it = _parts.begin(); it != _parts.end(); ++it)
        {
          // check size
          if(int(it->size()) != size)
            continue;

          // check name
          if(!names.empty())
          {
            bool match = false;
            for(const auto& name : names)
            {
              match = match || (name == it->get_name());
            }
            if(!match)
              continue;
          }

          // check priority
          if((it->get_priority() <= 0) || (it->get_priority() < prio))
            continue;

          // okay, that's a candidate
          // does it have a greater or equal priority to the previously found one?
          if((part != nullptr) && (part->get_priority() > it->get_priority()))
            continue;

          // we have found a matching partition
          part = &(*it);
        }

        return part;
      }
    }; // class PartitionSet
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_PARTITION_SET_HPP
