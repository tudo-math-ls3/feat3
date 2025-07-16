// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/assembly/asm_traits.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/intern/face_ref_trafo.hpp>
#include <kernel/geometry/intern/congruency_sampler.hpp>
#include <kernel/geometry/intern/congruency_trafo.hpp>
#include <kernel/geometry/intern/index_representative.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>

#include <vector>

namespace FEAT
{
  namespace Assembly
  {
    /// \cond internal
    namespace Intern
    {
      template<typename Outer_, typename Shape_, int cell_dim_>
      class CompIndexMap
      {
      protected:
        const Outer_& _outer;
        int _idx;

      public:
        explicit CompIndexMap(const Outer_& outer, int idx) : _outer(outer), _idx(idx) {}

        Index operator[](int i) const
        {
          return _outer[Geometry::Intern::FaceIndexMapping<Shape_, cell_dim_, 0>::map(_idx,i)];
        }
      }; // class CompIndexMap
    } // namespace Intern
    /// \endcond

#ifdef DOXYGEN
    /**
     * \brief Interface description of a trace assembly job
     *
     * \attention
     * This class is only defined for the doxygen documentation, because it only serves as a
     * duck-type class interface documentation.
     *
     * The only mandatory member of a job class is the nested Task class, which is instantiated by
     * each worker thread in the trace assembler. This outer job class can be used to store
     * coefficients and references to data that has to be used by each task during the assembly
     * process. Each newly created instance of the nested Task class receives a reference to the
     * outer job object, so it can retrieve all the necessary data from it inside the constructor.
     *
     * \author Peter Zajac
     */
    class TraceAssemblyJob
    {
    public:
      /**
       * \brief Trace assembly task class
       *
       * \attention
       * This class is only defined for the doxygen documentation, because it only serves as a
       * duck-type class interface documentation.
       *
       * An instance of this Task class is created by each worker thread at the beginning of the
       * assembly process and is used to perform the assembly on the set of the elements that
       * are assigned to that particular worker thread.
       *
       * At the very least, this class contains a mandatory constructor, the five member functions
       * prepare(), assemble(), scatter(), finish() and combine() as well as three constexpr
       * members #assemble_pairwise, #need_scatter and #need_combine, which determine whether the
       * assembly tasks supports cell-pairwise assembly (see below for details) and whether the
       * scatter and combine functions have to be called by the trace assembler at all. All member
       * functions except for scatter() and combine() are silently assumed to be thread-safe.
       *
       * Because the trace assembler assembles surface integrals over facets which may be shared by
       * two neighboring cells (unless the facet is a boundary or halo facet), there are two ways
       * how the assembler may proceed in this case: it may either assemble the trace integral on
       * both adjacent cells separately or it may assemble both cells at the same time.
       * If the task's #assemble_pairwise member is set to true and the assembler processes a facet
       * that is shared by two neighboring cells, then the assembler will call the task's prepare
       * function overload with 7 function parameters indicating that the task shall assemble on
       * that pair of neighboring cells at the same time, otherwise the assembler will process the
       * two neighboring cells successively and it will call the task's prepare function overload
       * with 4 function parameters for each of the two cells.
       *
       * In the case of boundary (or halo) facets, the assembler will always call the task's prepare
       * function overload with 4 parameters, independently of whether #assemble_pairwise is set to
       * true or false.
       *
       * \note Currently, the TraceAssembler does not perform any sort of multi-threading and so it
       * does not employ any worker threads, but this may change in future.
       *
       * <b><u>Basic Job Assembly work flow:</u></b>\n
       * Each worker thread employed by the trace assembler will perform the following
       * (simplified) steps (pseudo-code):
       \code{.cpp}
         // create a Task object for this worker thread
         Job::Task task(job);

         // loop over all facets of this worker thread
         for(auto facet : this_workers_facet_list)
         {
           // is this an inner facet and do we assemble pairwise?
           if(facet.is_inner && task.assemble_pairwise)
           {
             // prepare task for assembly on a pair of neighboring cells sharing a facet
             task.prepare(facet, cell_a, cell_b, ...);

             // assemble on current cell pair
             task.assemble();

             // does the task need to scatter?
             if(task.need_scatter)
             {
               // scatter local data
               neighbor_mutex.lock();
               task.scatter();
               neighbor_mutex.unlock();
             }

             // finish assembly on current cell pair
             task.finish();

             // continue with next facet
             continue;
           }

           // no pairwise assembly
           for(auto cell : facet.adjacent_cells)
           {
             // prepare task for assembly on one cell adjacent to the facet
             task.prepare(facet, cell, ...);

             // assemble on current cell pair
             task.assemble();

             // does the task need to scatter?
             if(task.need_scatter)
             {
               // scatter local data
               neighbor_mutex.lock();
               task.scatter();
               neighbor_mutex.unlock();
             }

             // finish assembly on current cell pair
             task.finish();
           }
         }

         // does the task need to combine?
         if(task.need_combine)
         {
           // combine local data
           domain_mutex.lock();
           task.combine();
           domain_mutex.unlock();
         }
       \endcode
       * In the pseudo-code above, the code lines between <c>neighbor_mutex.lock()</c> and
       * <c>neighbor_mutex.unlock()</c> are executed in a neighbor-mutually exclusive fashion,
       * i.e. it is guaranteed that this code (<c>task.scatter()</c>) is never executed by two
       * threads, which are active on two neighbor elements, at the same time. However, the code
       * may be executed by two threads simultaneously if the cells, that the two threads are
       * currently active on, are not adjacent to each other.\n
       * Furthermore, the code lines between <c>domain_mutex.lock()</c> and
       * <c>domain_mutex.unlock()</c> are executed in a (domain-) mutually exclusive fashion, i.e.
       * it is guaranteed that this code (<c>task.combine()</c>) is never executed by two threads
       * at the same time.
       *
       * \author Peter Zajac
       */
      class Task
      {
      public:
        /**
         * \brief Specifies whether this task assembles element pairs on inner facets simultaenously.
         */
        static constexpr bool assemble_pairwise = true or false;

        /**
         * \brief Specifies whether this task has a #scatter() function, which is required to be
         * called from within a critical section after the local assembly on each facet.
         */
        static constexpr bool need_scatter = true or false;

        /**
         * \brief Specifies whether this task fas a #combine() function, which is required to be
         * called from within a critical section after the overall assembly is finished.
         */
        static constexpr bool need_combine = true or false;

        /**
         * \brief Mandatory Constructor
         *
         * This constructor is used by each of the trace assemblers worker threads to create its
         * own task object.
         *
         * \param[in] job
         * A \resident reference to the encapsulating job class object.
         *
         * \attention
         * This constructor is silently assumed to be thread-safe, i.e. no thread may try to write
         * to a common shared resource inside this function (without manual mutexing).
         */
        explicit Task(TraceAssemblyJob& job);

        /**
         * \brief Prepares the task for assembly for a facet on a single one of its adjacent elements
         *
         * This function is called by the assembly worker thread each time the assembly on a new
         * element begins. This function is called for boundary facets as well as for each inner
         * facet and one of its adjacent elements if assemble_pairwise is false, otherwise the
         * 7-parameter overload of this function is called.
         * This function is typically used to gather various local data from the
         * mesh and/or linear algebra containers prior to the actual local assembly.
         *
         * \param[in] facet
         * The index of the mesh facet that the following assembly is going to take place on.
         *
         * \param[in] cell
         * The index of the mesh element that the following assembly is going to take place on.
         *
         * \param[in] local_facet
         * The index of the local facet on the element that the global facet corresponds to.
         *
         * \param[in] facet_ori
         * The facet orientation code that encodes the orientation between the local reference facet
         * and the global facet.
         *
         * \attention
         * This function is silently assumed to be thread-safe, i.e. no thread may try to write
         * to a common shared resource inside this function (without manual mutexing).
         */
        void prepare(Index facet, Index cell, int local_facet, int facet_ori);

        /**
         * \brief Prepares the task for assembly for a facet on both of its adjacent elements
         *
         * This function is called by the assembly worker thread each time the assembly on a new
         * inner facet begins, but only if #assemble_pairwise is true, otherwise the assembler will
         * call the 4-parameter overload of this function when assembling both adjacent elements
         * successively. This function is typically used to gather various local data from the
         * mesh and/or linear algebra containers prior to the actual local assembly.
         *
         * \param[in] facet
         * The index of the mesh facet that the following assembly is going to take place on.
         *
         * \param[in] cell_a, cell_b
         * The indices of the two mesh elements that share the common facet that the following assembly is going to take place on.
         *
         * \param[in] local_facet_a, local_facet_b
         * The indices of the local facets on the two elements that the global facet corresponds to.
         *
         * \param[in] facet_ori_a, facet_ori_b
         * The facet orientation codes that encodes the orientation between the local reference facets of both elements
         * and the global facet.
         *
         * \attention
         * This function is silently assumed to be thread-safe, i.e. no thread may try to write
         * to a common shared resource inside this function (without manual mutexing).
         */
        void prepare(Index facet, Index cell_a, Index cell_b, int local_facet_a, int local_facet_b, int facet_ori_a, int facet_ori_b);

        /**
         * \brief Performs the local assembly on the current facet
         *
         * This function is called by the assembly worker thread directly after the #prepare()
         * function and is typically used to perform the actual local assembly, i.e. this is where
         * most of the magic happens.
         *
         * \attention
         * This function is silently assumed to be thread-safe, i.e. no thread may try to write
         * to a common shared resource inside this function (without manual mutexing).
         */
        void assemble();

        /**
         * \brief Scatters the local assembly into the global system
         *
         * This function is called by the assembly worker thread after the #assemble() function
         * and is typically used to scatter the local (i.e. element-wise) assembled data into a
         * global (i.e. mesh/patch-wise) container.
         *
         * \note
         * This function is called by the worker threads in a neighbor-mutually exclusive fashion,
         * i.e. this function is never called by two (or more) threads at the same time, which are
         * active on vertex-adjacent elements. This ensures that no race condition can appear
         * when scattering a local matrix/vector into a global matrix/vector.
         *
         * \attention
         * This function is only called by the assembler if #need_scatter is true, however, it
         * must always exist even if #need_scatter is false to prevent linker errors.
         */
        void scatter();

        /**
         * \brief Finishes the task on the current cell.
         *
         * This function is called by the assembly worker thread after the #scatter() function
         * and is typically only used to call the finish member functions of all corresponding
         * member functions like e.g. evaluators.
         *
         * \attention
         * This function is silently assumed to be thread-safe, i.e. no thread may try to write
         * to a common shared resource inside this function (without manual mutexing).
         */
        void finish();

        /**
         * \brief Finishes the overall assembly and combines all local results
         *
         * This function is called once per assembly run after the worker thread has finished
         * assembling all its cells. This function can be used to collect (or reduce) the
         * information assembled in each task into a combined result in the job object.
         *
         * \note
         * This function is called by the worker threads in a mutually exclusive fashion,
         * i.e. this function is never called by two threads at the same time, so accessing
         * a common resource shared by all threads is safe inside this function.
         *
         * \attention
         * This function is only called by the assembler if #need_combine is true, however, it
         * must always exist even if #need_combine is false to prevent linker errors.
         */
        void combine();
      }; // class Task
    }; // class TraceAssemblyJob
#endif // DOXYGEN

    /**
     * \brief Trace Integral Assembler class template
     *
     * This class can be used to assemble operators, functionals and various other quantities on
     * (a part of) the boundary of a mesh or any other set of mesh facets.
     *
     * After you have constructed the trace assembler object, you still have to tell the assembler
     * on which set of facets the assembly should take place before you can actually assemble anything.
     * There are 2 ways to do this:
     * - You can tell the assembler to assemble on all inner and/or outer facets of the mesh by
     *   calling the #compile_all_facets() function.
     * - You can add individual facets or whole mesh parts to the assembler by calling the
     *   #add_facet() and #add_mesh_part() functions. Afterwards, you need to compile the assembler
     *   by calling the #compile() function.
     *
     * \note
     * Each trace assembler object is tied to the trafo object that is has been constructed on, i.e.
     * if you have several different trafo objects (such as in a multigrid scenario) then you also
     * need several trace assembler objects -- one for each trafo object. Also, a trace assembler
     * always assembles over the elements of its trafo, so all FE spaces used in the assembly have
     * to be defined on the same trafo object as the trace assembler.
     *
     * \tparam Trafo_
     * The transformation on whose underlying mesh the assembly should take place
     *
     * \note
     * This class currently does not employ any sort of multi-threading, however, the documentation
     * of the classes and functions related to this class may already be written as if this class
     * already supported assembly with multiple threads -- which it actually may do sometime in the
     * future.
     *
     * \author Peter Zajac
     */
    template<typename Trafo_>
    class TraceAssembler
    {
    public:
      typedef Trafo_ TrafoType;
      typedef typename TrafoType::MeshType MeshType;
      typedef typename TrafoType::ShapeType ShapeType;

      static constexpr int shape_dim = ShapeType::dimension;
      static constexpr int facet_dim = shape_dim-1;

    protected:
      /// a reference to the trafo
      const TrafoType& _trafo;
      /// the facet masks, local cell facet indices and facet orientation codes
      std::vector<int> _facet_mask, _cell_facet, _facet_ori;
      /// the indices of all cells and facets to loop over during assembly
      std::vector<Index> _cells, _facets;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] trafo
       * A \resident reference to the trafo on which to assemble
       */
      explicit TraceAssembler(const TrafoType& trafo) :
        _trafo(trafo),
        _facet_mask(trafo.get_mesh().get_num_entities(facet_dim), 0)
      {
      }

      /// \returns a reference to the trafo
      const TrafoType& get_trafo() const
      {
        return _trafo;
      }

      /**
       * \brief Clears the assembler
       */
      void clear()
      {
        _cell_facet.clear();
        _facet_ori.clear();
        _cells.clear();
        _facets.clear();
        for(auto& f : _facets)
          f = Index(0);
      }

      /**
       * \brief Adds a single facet to the assembler.
       *
       * \param[in] facet
       * The index of the facet to be added to the assembler.
       */
      void add_facet(Index ifacet)
      {
        ASSERTM(ifacet < Index(_facet_mask.size()), "invalid facet index");
        _facet_mask.at(ifacet) = 1;
      }

      /**
       * \brief Adds all facets of a mesh part to the assembler.
       *
       * \param[in] mesh_part
       * The mesh part whose facets are to be added.
       */
      void add_mesh_part(const Geometry::MeshPart<MeshType>& mesh_part)
      {
        const auto& trg = mesh_part.template get_target_set<facet_dim>();
        for(Index i(0); i < trg.get_num_entities(); ++i)
          _facet_mask.at(trg[i]) = 1;
      }

      /**
       * \brief Compiles the assembler for all facets that have been added manually.
       *
       * This function compiles the assembler for all facets that have been added manually by
       * previous calls of the #add_facet or #add_mesh_part functions.
       *
       * \note
       * If you want to compile the assembler for all facets of the trafo's underlying mesh,
       * consider using the #compile_all_facets() function instead.
       */
      void compile()
      {
        _cell_facet.clear();
        _facet_ori.clear();
        _cells.clear();
        _facets.clear();

        // build elements-at-facet graph
        Adjacency::Graph elem_at_facet(Adjacency::RenderType::injectify_transpose,
          _trafo.get_mesh().template get_index_set<shape_dim, facet_dim>());

        // loop over all facets
        for(Index iface(0); iface < Index(_facet_mask.size()); ++iface)
        {
          // use this facet?
          if(_facet_mask[iface] == 0)
            continue;

          // ensure that this is a boundary facet if required
          //if(only_boundary && (elem_at_facet.degree(iface) != Index(1)))
            //XABORTM("facet is adjacent to more than 1 element");

          // add all elements
          for(auto it = elem_at_facet.image_begin(iface); it != elem_at_facet.image_end(iface); ++it)
          {
            Index icell = *it;

            // try to compute local facet index and orientation
            int loc_face(0), face_ori(0);
            if(!_find_local_facet(iface, icell, loc_face, face_ori))
              XABORTM("failed to find local facet");

            // alright, add this facet to our list
            _facets.push_back(iface);
            _cells.push_back(icell);
            _cell_facet.push_back(loc_face);
            _facet_ori.push_back(face_ori);
          }
        }
      }

      /**
       * \brief Compiles the assembler for all inner and/our outer facets of the underlying mesh
       *
       * \param[in] inner
       * Specifies whether the mesh's inner facets are to be added to the assembler or not.
       *
       * \param[in] outer
       * Specifies whether the mesh's outer (i.e. boundary) facets are to be added to the
       * assembler or not.
       */
      void compile_all_facets(bool inner, bool outer)
      {
        _cell_facet.clear();
        _facet_ori.clear();
        _cells.clear();
        _facets.clear();

        // build elements-at-facet graph
        Adjacency::Graph elem_at_facet(Adjacency::RenderType::injectify_transpose,
          _trafo.get_mesh().template get_index_set<shape_dim, facet_dim>());

        // loop over all facets
        for(Index iface(0); iface < Index(_facet_mask.size()); ++iface)
        {
          // how many elements are adjacent to this facet?
          const int degree = int(elem_at_facet.degree(iface));
          XASSERT(degree > 0);
          XASSERT(degree < 3);

          // outer facet with only 1 adjacent element?
          if((degree == 1) && !outer)
            continue;

          // inner facet with 2 adjacent elements?
          if((degree == 2) && !inner)
            continue;

          // add all elements
          for(auto it = elem_at_facet.image_begin(iface); it != elem_at_facet.image_end(iface); ++it)
          {
            Index icell = *it;

            // try to compute local facet index and orientation
            int loc_face(0), face_ori(0);
            if(!_find_local_facet(iface, icell, loc_face, face_ori))
              XABORTM("failed to find local facet");

            // alright, add this facet to our list
            _facets.push_back(iface);
            _cells.push_back(icell);
            _cell_facet.push_back(loc_face);
            _facet_ori.push_back(face_ori);
          }
        }
      }

      /**
       * \brief Executes a trace assembly job
       *
       * \param[inout] job
       * A \transient reference to the job that is to be assembled.
       * See Assembly::TraceAssemblyJob for details.
       */
      template<typename Job_>
      void assemble(Job_& job)
      {
        typedef typename Job_::Task TaskType;

        //XASSERTM(_compiled, "assembler has not been compiled yet");

        // no facets to assemble on?
        if(this->_facets.empty())
          return;

        // create a task for this thread
        std::unique_ptr<TaskType> task(new TaskType(job));

        // note: in the case of inner facets, a single inner facet is stored twice in the _facets
        // vector in two successive elements, i.e. one has _facets[i] == _facets[i+1], and the
        // corresponding two adjacent cells are given by _cells[i] and _cells[i+1]

        // loop over all facets that have been added to the assembler
        for(Index fi(0); fi < Index(_facets.size()); ++fi)
        {
          // prepare task -- either pairwise or not
          if constexpr(TaskType::assemble_pairwise)
          {
            // task supports pairwise assembly, so let's check if the next facet is an inner facet
            Index fj = fi+1u;
            if((fj < Index(_facets.size())) && (_facets[fi] == _facets[fj]))
            {
              // that's an inner facet, so assemble pairwise here
              task->prepare(_facets[fi], _cells[fi], _cells[fj], _cell_facet[fi], _cell_facet[fj], _facet_ori[fi], _facet_ori[fj]);

              // make sure we don't try to assemble on fi+1 again
              ++fi;
            }
            else
            {
              // not an inner facet, so call single-cell prepare overload
              task->prepare(_facets[fi], _cells[fi], _cell_facet[fi], _facet_ori[fi]);
            }
          }
          else
          {
            // task does not support pairwise assembly, so always call the single-cell prepare overload
            task->prepare(_facets[fi], _cells[fi], _cell_facet[fi], _facet_ori[fi]);
          }

          // assemble task
          task->assemble();

          // scatter
          if(task->need_scatter)
          {
            task->scatter();
          }

          // finish
          task->finish();
        }

        // do we have to combine the assembly?
        if(task->need_combine)
        {
          // combine the assembly
          task->combine();
        }

        // delete task object
        task.reset();
      }

    protected:
      /**
       * \brief Helper function: tries to find the local facet index for a given facet/cell pair
       *
       * \param[in] face
       * The global index of the facet that is adjacent to \p cell
       *
       * \param[in] cell
       * The global index of the cell/element that is adjacent to \p face
       *
       * \param[out] facet
       * The local facet index of \p face with respect to \p cell
       *
       * \param[out] ori
       * The local facet orientation code of \p face with respect to the corresponding reference
       * element facet of \p cell
       *
       * \returns
       * \c true, if the local facet was identified, or \c false, if \p face does not seem to be
       * a facet of \p cell
       */
      bool _find_local_facet(Index face, Index cell, int& facet, int& ori)
      {
        typedef typename Shape::FaceTraits<ShapeType, facet_dim>::ShapeType FacetType;
        static constexpr int num_facets = Shape::FaceTraits<ShapeType, shape_dim-1>::count;
        static constexpr int num_vaf = Shape::FaceTraits<FacetType, 0>::count;

        typedef Geometry::Intern::IndexRepresentative<FacetType> FacetRepr;

        const auto& vert_at_elem = _trafo.get_mesh().template get_index_set<shape_dim, 0>();
        const auto& vert_at_face = _trafo.get_mesh().template get_index_set<facet_dim, 0>();

        Index face_verts[num_vaf];
        FacetRepr::compute(face_verts, vert_at_face[face]);

        for(int li(0); li < num_facets; ++li)
        {
          Intern::CompIndexMap<decltype(vert_at_elem[0]), ShapeType, shape_dim-1> cim(vert_at_elem[cell], li);
          Index lf_verts[num_vaf];
          FacetRepr::compute(lf_verts, cim);
          bool match = true;
          for(int k(0); k < num_vaf; ++k)
          {
            if(lf_verts[k] != face_verts[k])
            {
              match = false;
              break;
            }
          }
          if(match)
          {
            ori = Geometry::Intern::CongruencySampler<FacetType>::compare(vert_at_face[face], cim);
            facet = li;
            return true;
          }
        }
        return false;
      }
    }; // class TraceAssembler<...>
  } // namespace Assembly
} // namespace FEAT
