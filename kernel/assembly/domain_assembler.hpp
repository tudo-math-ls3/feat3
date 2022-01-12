// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ASSEMBLY_DOMAIN_ASSEMBLER_HPP
#define KERNEL_ASSEMBLY_DOMAIN_ASSEMBLER_HPP 1

// includes, FEAT
#include <kernel/assembly/asm_traits.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/adjacency/coloring.hpp>
#include <kernel/util/thread.hpp>

// includes, system
#include <algorithm>
#include <memory>
#include <vector>

namespace FEAT
{
  namespace Assembly
  {
    /**
     * \brief Threading Strategy for multi-threaded assembler
     *
     * This enumeration defines the various threading strategies supported by
     * the DomainAssembler class template.
     */
    enum class ThreadingStrategy
    {
      /**
       * \brief Automatic threading strategy
       *
       * This is the recommended default threading strategy.
       *
       * This value specifies that the domain assembler should choose an appropriate threading
       * strategy for the given number of maximum worker threads however it sees fit.
       *
       * In its current implementation, the domain assembler will automatically choose the
       * strategy based on the following criterions:
       * - If only 1 thread is requested, ThreadingStrategy::single is chosen.
       * - If the underlying mesh is permuted using the Geometry::PermutationStrategy::colored
       *   permutation strategy, then ThreadingStrategy::colored is chosen.
       * - In any other case, ThreadingStrategy::layered is chosen.
       */
      automatic = 0,

      /**
       * \brief Single-Threaded strategy
       *
       * This strategy corresponds to the single-threading approach, i.e. no multi-threading
       * is performed.
       */
      single,

      /**
       * \brief Layered threading strategy
       *
       * This strategy performs a (reverse) Cuthill-McKee ordering of the mesh
       * elements to obtain a sequence of element layers which can be processed in
       * parallel. By construction of the Cuthill-McKee layers, race conditions can
       * only appear within two adjacent element layers, so it is sufficient to
       * ensure that no two threads process two adjacent layers at the same time.
       *
       * <b>Example:</b>\n
       * Consider the following 5x5 square mesh, then the Cuthill-McKee algorithm
       * could produce the following 5 element layers A, B, C, D and E.
       *
         \verbatim
         +---+---+---+---+---+
         | E | E | E | E | E |
         +---+---+---+---+---+
         | D | D | D | D | E |
         +---+---+---+---+---+
         | C | C | C | D | E |
         +---+---+---+---+---+
         | B | B | C | D | E |
         +---+---+---+---+---+
         | A | B | C | D | E |
         +---+---+---+---+---+
         \endverbatim
       *
       * In a two-threaded assembly, the first thread may process the layers
       * E and D, where as the second thread may process the layers C, B, and A.
       * The only critical section are the two adjacent layers D and C, which
       * may not be processed at the same time by the two threads.
       * To avoid race conditions here, the second thread will notify the first
       * thread once it has finished processing its first layer C, so that the
       * first thread may proceed to process its last layer D.
       *
       * In the general case, each thread is assigned 2 or more adjacent
       * Cuthill-McKee layers, so that the only required interaction between
       * threads takes place at the adjacent layers of two adjacent threads.
       * Just as in the simple example above, each thread notifies its preceding
       * thread once it has finished processing its first layer and each thread
       * has to wait for notification from its succeeding thread before its starts
       * processing its last layer to avoid race conditions on adjacent layers.
       */
      layered,

      /**
       * \brief Layered + sorted threading strategy
       *
       * This threading strategy is mostly identical to the (simple) layered strategy, where the
       * only difference is that all entities in each layer are sorted by ascending degree, i.e.
       * by the number of neighbor entities. This represents the 'classical' Cuthill-McKee
       * ordering algorithms for simple graphs.
       *
       * In practical experiments -- especially when combined with mesh permutations -- this
       * layered+sorted strategy is often LESS efficient than the simple (unsorted) layered
       * strategy, so therefore it is recommended to use the unsorted strategy.
       * This sorted strategy is implemented primarily for benchmarking reasons.
       */
      layered_sorted,

      /**
       * \brief Colored threading strategy
       *
       * This strategy performs a coloring of the mesh elements, so that
       * any two adjacent elements have different colors. By construction,
       * no race conditions can appear when elements of a single color
       * are processed in parallel. Therefore, the assembly processes all
       * colors sequentially, where all elements of a single color are
       * processed in parallel by multiple threads.
       *
       * <b>Example:</b>\n
       * Consider the following 5x5 square mesh, then the elements could
       * be partitioned into the four colors R, G, B and W.
       *
         \verbatim
         +---+---+---+---+---+
         | R | G | R | G | R |
         +---+---+---+---+---+
         | B | W | B | W | B |
         +---+---+---+---+---+
         | R | G | R | G | R |
         +---+---+---+---+---+
         | B | W | B | W | B |
         +---+---+---+---+---+
         | R | G | R | G | R |
         +---+---+---+---+---+
         \endverbatim
       *
       * Although this strategy is the more prominent one in the literature,
       * it is usually less efficient than the layered strategy due to higher
       * synchronization costs and less favorable memory access patterns.
       */
      colored
    }; // enum class ThreadingStrategy

#ifdef DOXYGEN
    /**
     * \brief Interface description of a domain assembly job
     *
     * \attention
     * This class is only defined for the doxygen documentation, because it only serves as a
     * duck-type class interface documentation.
     *
     * The only mandatory member of a job class is the nested Task class, which is instantiated by
     * each worker thread in the domain assembler. This outer job class can be used to store
     * coefficients and references to data that has to be used by each task during the assembly
     * process. Each newly created instance of the nested Task class receives a reference to the
     * outer job object, so it can retrieve all the necessary data from it inside the constructor.
     *
     * \author Peter Zajac
     */
    class DomainAssemblyJob
    {
    public:
      /**
       * \brief Domain assembly task class
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
       * prepare(), assemble(), scatter(), finish() and combine() as well as two constexpr
       * members #need_scatter and #need_combine, which determine whether the scatter and combine
       * functions have to be called by the domain assembler at all. All member functions except
       * for scatter() and combine() are silently assumed to be thread-safe.
       *
       * <b><u>Basic Job Assembly workflow:</u></b>\n
       * Each worker thread employed by the domain assembler will perform the following
       * (simplified) steps:
         \code{.cpp}
         // create a Task object for this worker thread
         Job::Task task(job);

         // loop over all cells of this worker thread
         for(auto cell : this_workers_cell_list)
         {
           // prepare task for assembly on current cell
           task.prepare(cell);

           // assemble on current cell
           task.assemble();

           // does the task need to scatter?
           if(task.need_scatter)
           {
             // scatter local data
             neighbor_mutex.lock();
             task.scatter();
             neighbor_mutex.unlock();
           }

           // finish assembly on current cell
           task.finish();
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
         * \brief Specifies whether this task has a #scatter() function, which is required to be
         * called from within a critical section after the local assembly on each cell.
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
         * This constructor is used by each of the domain assemblers worker threads to create its
         * own task object.
         *
         * \param[in] job
         * A \resident reference to the encapsulating job class object.
         *
         * \attention
         * This constructor is silently assumed to be thread-safe, i.e. no thread may try to write
         * to a common shared resource inside this function (without manual mutexing).
         */
        explicit Task(DomainAssemblyJob& job);

        /**
         * \brief Prepares the task for assembly on a element/cell
         *
         * This function is called by the assembly worker thread each time the assembly on a new
         * element begins. This function is typically used to gather various local data from the
         * mesh and/or linear algebra containers prior to the actual local assembly.
         *
         * \param[in] cell
         * The index of the mesh element/cell that the following assembly is going to take place on.
         *
         * \attention
         * This function is silently assumed to be thread-safe, i.e. no thread may try to write
         * to a common shared resource inside this function (without manual mutexing).
         */
        void prepare(Index cell);

        /**
         * \brief Performs the local assembly on the current cell
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
    }; // class DomainAssemblyJob
#endif // DOXYGEN

    /**
     * \brief Domain Integral Assembler class template
     *
     * This class implements a multi-threaded domain assembler, which is capable of executing
     * domain assembly jobs.
     *
     * After you have constructed the domain assembler object, there are still a few steps left
     * that have to be performed before you can assemble a assembly job by the domain assembler:
     * - Optional: You can choose a specific threading strategy by calling the #set_threading_strategy
     *   member function or you can just stick with the default automatic strategy that is preselected.
     * - Optional: You can choose the maximum number of desired worker threads for the multi-threaded
     *   assembly by calling the #set_max_worker_threads member function. Skipping this step always
     *   results in single-threaded assembly.
     * - If you want to assemble on the whole domain, i.e. on all elements of the mesh (which
     *   is the usual case), all you have to do as a final mandatory step is to call the
     *   #compile_all_elements() member function.
     * - If you want to assemble on a subset of the domain, you have to add all elements to the
     *   domain assembler by either adding them individually using the #add_element() function
     *   or by adding them as mesh-parts using the #add_mesh_part() member function. Once all
     *   elements have been added, you have to call the #compile() member function.
     *
     * \note
     * Each domain assembler object is tied to the trafo object that is has been constructed on, i.e.
     * if you have several different trafo objects (such as in a multigrid scenario) then you also
     * need several domain assembler objects -- one for each trafo object. Also, a domain assembler
     * always assembles over the elements of its trafo, so all FE spaces used in the assembly have
     * to be defined on the same trafo object as the domain assembler.
     *
     * \attention
     * Please be aware that the setup of a DomainAssembler object is a rather expensive operation,
     * because the assembler has to perform a lot of internal initialization including an analysis
     * of the element adjacency pattern as well as applying graph coloring/layering algorithms,
     * which (in general) may have super-linear runtime. Therefore you should always hang on to a
     * DomainAssembler object once it has been created and avoid re-creation of several
     * DomainAssembler objects that could have been avoided by saving a previously created
     * DomainAssembler object.
     *
     * \author Peter Zajac
     */
    template<typename Trafo_>
    class DomainAssembler
    {
    public:
      /// the underlying trafo type
      typedef Trafo_ TrafoType;
      /// the underlying mesh type
      typedef typename TrafoType::MeshType MeshType;
      /// the shape dimension
      static constexpr int shape_dim = MeshType::shape_dim;

      /**
       * \brief Thread statistics helper class
       *
       * This class collects some assembly statistics for a given worker thread, most notably
       * the assembly and waiting-for-mutex-lock timings.
       */
      class ThreadStats
      {
      public:
        /// microseconds assembling in total
        long long micros_total;
        /// microseconds spend in actual assembly
        long long micros_assemble;
        /// microseconds spend waiting for mutex locks
        long long micros_wait;

      public:
        ThreadStats() :
          micros_total(0ll),
          micros_assemble(0ll),
          micros_wait(0ll)
        {
        }

        void reset()
        {
          micros_total = 0ll;
          micros_assemble = 0ll;
          micros_wait = 0ll;
        }

        ThreadStats& operator+=(const ThreadStats& other)
        {
          micros_total    += other.micros_total;
          micros_assemble += other.micros_assemble;
          micros_wait     += other.micros_wait;
          return *this;
        }
      }; // class ThreatStats

    protected:
      /// \cond internal
      /// helper class for std::sort; sort domain nodes by their degree
      class DegreeCompare
      {
      public:
        const Adjacency::Graph& _graph;
        explicit DegreeCompare(const Adjacency::Graph& graph) : _graph(graph) {}
        bool operator()(Index i, Index j) const {return _graph.degree(i) < _graph.degree(j);}
      }; // class DegreeCompare
      /// \endcond

      /**
       * \brief Worker thread data class
       *
       * This class encapsulates all data that is required by a worker thread
       * to execute an assembly job.
       */
      template<typename Job_>
      class Worker
      {
      private:
        /// a typedef for the task
        typedef typename Job_::Task TaskType;
        /// a reference to the assembly job
        Job_& _job;
        /// id of this worker thread and total number of worker threads
        const std::size_t _my_id, _num_workers;
        /// the chosen threading strategy
        const ThreadingStrategy _strategy;
        /// the free-to-use thread mutex
        std::mutex& _thread_mutex;
        /// the thread statistics
        ThreadStats& _thread_stats;
        /// the thread fences vector
        std::vector<ThreadFence>& _thread_fences;
        /// the element indices vector
        const std::vector<Index>& _element_indices;
        /// the color elements vector
        const std::vector<Index>& _color_elements;
        /// the layer elements vector
        const std::vector<Index>& _layer_elements;
        /// the thread layers vector
        const std::vector<Index>& _thread_layers;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] job
         * A \resident reference to the assembly job.
         *
         * \param[in] id
         * The id of this worker thread.
         * Is 0 if the master thread is assembling, otherwise 1 <= id <= num_workers.
         *
         * \param[in] num_workers
         * The total number of worker threads.
         * May be 0 if the master thread is assembling.
         *
         * \param[in] strategy
         * The chosen threading strategy
         *
         * \param[inout] thread_stats
         * A \resident reference to the thread's statistics object.
         *
         * \param[inout] thread_mutex
         * A \resident reference to the common thread mutex.
         *
         * \param[inout] thread_fences
         * A \resident reference to the thread fence vector.
         *
         * \param[in] element_indices
         * A \resident reference to the element indices vector.
         *
         * \param[in] color_elements
         * A \resident reference to the color elements vector.
         *
         * \param[in] layer_elements
         * A \resident reference to the layer elements vector.
         *
         * \param[in] thread_layers
         * A \resident reference to the thread layers vector.
         */
        explicit Worker(Job_& job, std::size_t id, std::size_t num_workers,
          ThreadingStrategy strategy,
          ThreadStats& thread_stats,
          std::mutex& thread_mutex,
          std::vector<ThreadFence>& thread_fences,
          const std::vector<Index>& element_indices,
          const std::vector<Index>& color_elements,
          const std::vector<Index>& layer_elements,
          const std::vector<Index>& thread_layers) :
          _job(job),
          _my_id(id),
          _num_workers(num_workers),
          _strategy(strategy),
          _thread_mutex(thread_mutex),
          _thread_stats(thread_stats),
          _thread_fences(thread_fences),
          _element_indices(element_indices),
          _color_elements(color_elements),
          _layer_elements(layer_elements),
          _thread_layers(thread_layers)
        {
        }

        /// virtual destructor
        virtual ~Worker()
        {
        }

        /**
         * \brief Evaluation operator
         *
         * This operator is called by each thread to execute the assembly.
         */
        void operator()()
        {
          bool okay = false;

          TimeStamp stamp_total;

          // put everything in a try-catch block
          try
          {
            // create the task
            std::unique_ptr<TaskType> task(new TaskType(_job));

            // choose the appropriate work function for this task
            if(this->_num_workers <= std::size_t(1))
            {
              // we only have 1 thread
              okay = this->_work_single(std::move(task));
            }
            else if(!task->need_scatter)
            {
              // we have multiple threads, but the task does not need to scatter
              okay = this->_work_no_scatter(std::move(task));
            }
            else if(this->_strategy == ThreadingStrategy::colored)
            {
              // multiple threads and colored threading strategy
              okay = this->_work_colored(std::move(task));
            }
            else
            {
              // multiple threads and layered threading strategy
              okay = this->_work_layered(std::move(task));
            }
          }
          catch(...)
          {
            okay = false;
          }

          // if something went wrong, then we'll notify the master
          // thread by setting our worker thread fence's status to false.
          if(!okay)
          {
            this->_thread_fences.at(this->_my_id).open(false);
          }

          // update timing statistics
          this->_thread_stats.micros_total += stamp_total.elapsed_micros_now();
        }

      protected:
        /**
         * \brief Assembly worker implementation for single-threaded strategy
         *
         * This member function implements the actual assembly for the single-threaded strategy.
         *
         * \param[in] task
         * A unique pointer to the task that is to be executed.
         */
        bool _work_single(std::unique_ptr<TaskType> task)
        {
          XASSERTM(this->_my_id == std::size_t(0), "invalid threading strategy");
          XASSERTM(this->_num_workers <= std::size_t(1), "invalid threading strategy");

          Index elem_beg = Index(0);
          Index elem_end = Index(this->_element_indices.size());

          // create assembly time stamp
          TimeStamp stamp_asm;

          // loop over all elements in this thread's layers
          for(Index elem(elem_beg); elem < elem_end; ++elem)
          {
            // prepare task
            task->prepare(this->_element_indices.at(elem));

            // assemble task
            task->assemble();

            // scatter
            if(task->need_scatter)
              task->scatter();

            // finish
            task->finish();
          }

          // finalize the assembly
          if(task->need_combine)
            task->combine();

          // save elapsed time
          this->_thread_stats.micros_assemble += stamp_asm.elapsed_micros_now();

          // delete task object
          task.reset();

          // okay, we're done here
          return true;
        }

        /**
         * \brief Assembly worker implementation for no-scatter assembly
         *
         * This member function implements the actual no-scatter assembly, which does not perform
         * any synchronization to avoid race conditions except for the optional mutexed combine.
         *
         * \param[in] task
         * A unique pointer to the task that is to be executed.
         */
        bool _work_no_scatter(std::unique_ptr<TaskType> task)
        {
          XASSERTM(this->_num_workers > std::size_t(1), "invalid threading strategy");

          Index elem_beg = Index(((this->_my_id-1u) * this->_element_indices.size()) / this->_num_workers);
          Index elem_end = Index(((this->_my_id   ) * this->_element_indices.size()) / this->_num_workers);

          // create assembly time stamp
          TimeStamp stamp_asm, stamp_wait;

          // loop over all elements in this thread's layers
          for(Index elem(elem_beg); elem < elem_end; ++elem)
          {
            // prepare task
            task->prepare(this->_element_indices.at(elem));

            // assemble task
            task->assemble();

            // finish
            task->finish();
          }

          // do we have to combine the assembly?
          if(task->need_combine)
          {
            // start waiting stamp and update assembly time before opening the fence
            this->_thread_stats.micros_assemble += stamp_wait.stamp().elapsed_micros(stamp_asm);

            // acquire lock for the thread mutex
            std::unique_lock<std::mutex> lock(this->_thread_mutex);

            // start assembly stamp and update waiting time
            this->_thread_stats.micros_wait += stamp_asm.stamp().elapsed_micros(stamp_wait);

            // combine the assembly
            task->combine();
          }

          // save elapsed time
          this->_thread_stats.micros_assemble += stamp_asm.elapsed_micros_now();

          // delete task object
          task.reset();

          // okay, we're done here
          return true;
        }

        /**
         * \brief Assembly worker implementation for layered (+sorted) strategy
         *
         * This member function implements the actual assembly for the layered threading strategy.
         *
         * \param[in] task
         * A unique pointer to the task that is to be executed.
         */
        bool _work_layered(std::unique_ptr<TaskType> task)
        {
          Index elem_beg = Index(0);
          Index elem_end = Index(this->_element_indices.size());
          Index elem_fence_open = ~Index(0);
          Index elem_fence_wait = ~Index(0);

          TimeStamp stamp_asm, stamp_wait;

          // multi-threaded assembly ?
          if(this->_my_id > 0u)
          {
            // first and last element for this thread
            elem_beg = this->_layer_elements.at(this->_thread_layers.at(_my_id-1));
            elem_end = this->_layer_elements.at(this->_thread_layers.at(_my_id));

            // last element of first layer: open fence of previous thread
            if(this->_my_id > 1u)
              elem_fence_open = this->_layer_elements.at(this->_thread_layers.at(_my_id-1) + 1u) - 1u;

            // first element of last layer: wait for fence of next thread
            if(this->_my_id < this->_num_workers)
              elem_fence_wait = this->_layer_elements.at(this->_thread_layers.at(_my_id) - 1u);

            // make sure that we do not enter a deadlock
            if((elem_fence_open != ~Index(0)) && (elem_fence_wait != ~Index(0)))
            {
              // note that this case should never happen, because it should be prevented
              // by the DomainAssembler::_build_thread_layer() function
              XASSERTM(elem_fence_open < elem_fence_wait, "potential deadlock detected");
            }
          }

          stamp_wait.stamp();

          // wait for start fence to open
          if(!this->_thread_fences.front().wait())
            return false;

          this->_thread_stats.micros_wait += stamp_wait.elapsed_micros_now();

          // start assembly stamp
          stamp_asm.stamp();

          // loop over all elements in this thread's layers
          for(Index elem(elem_beg); elem < elem_end; ++elem)
          {
            // prepare task
            task->prepare(this->_element_indices.at(elem));

            // assemble task
            task->assemble();

            // do we have to scatter?
            if(task->need_scatter)
            {
              // first element of last layer?
              if(elem == elem_fence_wait)
              {
                // start waiting stamp and update assembly time before waiting for the fence
                this->_thread_stats.micros_assemble += stamp_wait.stamp().elapsed_micros(stamp_asm);

                // wait for next thread to finish its first layer
                if(!this->_thread_fences.at(this->_my_id+1).wait())
                  return false;

                // start assembly stamp and update waiting time
                this->_thread_stats.micros_wait += stamp_asm.stamp().elapsed_micros(stamp_wait);
              }

              // perform the scatter
              task->scatter();

              // last element of first layer?
              if(elem == elem_fence_open)
              {
                // start waiting stamp and update assembly time before opening the fence
                this->_thread_stats.micros_assemble += stamp_wait.stamp().elapsed_micros(stamp_asm);

                // signal the previous thread that we have finished our first layer
                this->_thread_fences.at(this->_my_id).open(true);

                // start assembly stamp and update waiting time
                this->_thread_stats.micros_wait += stamp_asm.stamp().elapsed_micros(stamp_wait);
              }
            }

            // finish
            task->finish();
          }

          // do we have to combine the assembly?
          if(task->need_combine)
          {
            // start waiting stamp and update assembly time before opening the fence
            this->_thread_stats.micros_assemble += stamp_wait.stamp().elapsed_micros(stamp_asm);

            // acquire lock for the thread mutex
            std::unique_lock<std::mutex> lock(this->_thread_mutex);

            // start assembly stamp and update waiting time
            this->_thread_stats.micros_wait += stamp_asm.stamp().elapsed_micros(stamp_wait);

            // combine the assembly
            task->combine();
          }

          // update assembly time
          this->_thread_stats.micros_assemble += stamp_asm.elapsed_micros_now();

          // delete task object
          task.reset();

          // okay, we're done here
          return true;
        }

        /**
         * \brief Assembly worker implementation for colored strategy
         *
         * This member function implements the actual assembly for the colored threading strategy.
         *
         * \param[in] task
         * A unique pointer to the task that is to be executed.
         */
        bool _work_colored(std::unique_ptr<TaskType> task)
        {
          TimeStamp stamp_asm, stamp_wait;

          // loop over all layers/colors
          for(std::size_t icol(0); icol+1u < this->_color_elements.size(); ++icol)
          {
            // get number of elements for this color
            Index color_offs = this->_color_elements.at(icol);
            Index color_size = this->_color_elements.at(icol+1u) - this->_color_elements.at(icol);

            // first/last element to assemble on for this thread
            // note: it is perfectly valid that elem_beg = elem_end if the current color has less
            // elements than there are worker threads available; in this case the corresponding
            // worker threads simply do not assemble anything for this color
            Index elem_beg = Index(0);
            Index elem_end = color_size;
            if(this->_my_id > 0u)
            {
              elem_beg = (color_size * Index(this->_my_id-1)) / Index(this->_num_workers);
              elem_end = (color_size * Index(this->_my_id  )) / Index(this->_num_workers);
            }

            // start waiting stamp and update assembly time before waiting for the fence
            this->_thread_stats.micros_assemble += stamp_wait.stamp().elapsed_micros(stamp_asm);

            // wait for start signal
            if(!this->_thread_fences.front().wait())
              return false;

            // start assembly stamp and update waiting time
            this->_thread_stats.micros_wait += stamp_asm.stamp().elapsed_micros(stamp_wait);

            // loop over all elements
            for(Index elem(elem_beg); elem < elem_end; ++elem)
            {
              // prepare task
              task->prepare(this->_element_indices.at(color_offs + elem));

              // assemble task
              task->assemble();

              // scatter
              task->scatter();

              // finish
              task->finish();
            }

            // start waiting stamp and update assembly time before waiting for the fence
            this->_thread_stats.micros_assemble += stamp_wait.stamp().elapsed_micros(stamp_asm);

            // notify master that we're ready
            this->_thread_fences.at(this->_my_id).open(true);

            // wait for end signal
            if(!this->_thread_fences.back().wait())
              return false;

            // notify master that we're ready
            this->_thread_fences.at(this->_my_id).open(true);

            // start assembly stamp and update waiting time
            this->_thread_stats.micros_wait += stamp_asm.stamp().elapsed_micros(stamp_wait);
          } // next color layer

          // do we have to combine the assembly?
          if(task->need_combine)
          {
            // start waiting stamp and update assembly time before opening the fence
            this->_thread_stats.micros_assemble += stamp_wait.stamp().elapsed_micros(stamp_asm);

            // acquire lock for the thread mutex
            std::unique_lock<std::mutex> lock(this->_thread_mutex);

            // start assembly stamp and update waiting time
            this->_thread_stats.micros_wait += stamp_asm.stamp().elapsed_micros(stamp_wait);

            // combine the assembly
            task->combine();
          }

          // update assembly time
          this->_thread_stats.micros_assemble += stamp_asm.elapsed_micros_now();

          // delete task object
          task.reset();

          // okay
          return true;
        }
      }; // template class Worker<Job_>

    protected:
      /// a reference to the underlying trafo
      const TrafoType& _trafo;
      /// adjacency graph for vertices-at-element
      Adjacency::Graph _verts_at_elem;
      /// adjacency graph for elements-at-vertex
      Adjacency::Graph _elems_at_vert;
      /// adjacency graph for element neighbors
      Adjacency::Graph _elem_neighbors;
      /// an element mask vector
      std::vector<char> _element_mask;
      /// a vector of all elements to assemble on
      std::vector<Index> _element_indices;
      /// a vector of element color offsets
      std::vector<Index> _color_elements;
      /// a vector of element layer offsets
      std::vector<Index> _layer_elements;
      /// a vector of thread layer blocks
      std::vector<Index> _thread_layers;
      /// a vector of thread fences
      std::vector<ThreadFence> _thread_fences;
      /// a vector of thread statistics
      std::vector<ThreadStats> _thread_stats;
      /// a vector of worker threads
      std::vector<std::thread> _threads;
      /// a mutex for free use by the worker threads
      std::mutex _thread_mutex;
      /// specifies the chosen threading strategy
      ThreadingStrategy _strategy;
      /// specifies the maximum number of worker threads to use
      std::size_t _max_worker_threads;
      /// specifies the actual number of worker threads to use
      std::size_t _num_worker_threads;
      /// specifies whether the assembler has already been compiled
      bool _compiled;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] trafo
       * A \resident reference to the transformation that defines the assembly domain.
       */
      explicit DomainAssembler(const TrafoType& trafo) :
        _trafo(trafo),
        _verts_at_elem(),
        _elems_at_vert(),
        _elem_neighbors(),
        _element_mask(trafo.get_mesh().get_num_elements(), 0),
        _element_indices(),
        _color_elements(),
        _layer_elements(),
        _thread_layers(),
        _thread_fences(),
        _thread_stats(),
        _threads(),
        _strategy(ThreadingStrategy::automatic),
        _max_worker_threads(0),
        _num_worker_threads(0),
        _compiled(false)
      {
      }

      /// delete copy constructor
      DomainAssembler(const DomainAssembler&) = delete;
      /// delete copy assignment operator
      DomainAssembler& operator=(const DomainAssembler&) = delete;

      /// virtual destructor
      virtual ~DomainAssembler()
      {
        // go on, there's nothing to see here
      }

      /**
       * \brief Clears the assembler.
       */
      void clear()
      {
        XASSERTM(_threads.empty(), "currently executing a job");
        _verts_at_elem.clear();
        _elems_at_vert.clear();
        _elem_neighbors.clear();
        _element_indices.clear();
        _element_mask.clear();
        _color_elements.clear();
        _layer_elements.clear();
        _thread_layers.clear();
        _thread_fences.clear();
        _thread_stats.clear();
        _threads.clear();
        _compiled = false;
      }

      /**
       * \brief Adds a single element to the assembler.
       *
       * \param[in] ielem
       * The index of the mesh element that is to be added.
       *
       * \note
       * If the element \p ielem has already been added to the assembler by calling this function
       * or the #add_mesh_part function, then it will not be added a second time, so every element
       * is processed at most once, no matter how often it has been added to the assembler.
       */
      void add_element(Index ielem)
      {
        XASSERTM(!_compiled, "assembler has already been compiled!");
        this->_element_mask.at(ielem) = 1;
      }

      /**
       * \brief Adds all elements of a mesh-part to the assembler.
       *
       * \param[in] mesh_part
       * A \transient reference to the mesh part whose elements are to be added.
       *
       * \note
       * If an element in the mesh-part has already been added to the assembler by calling this
       * function or the #add_element function, then it will not be added a second time, so every
       * element is processed at most once, no matter how often it has been added to the assembler.
       */
      void add_mesh_part(const Geometry::MeshPart<MeshType>& mesh_part)
      {
        XASSERTM(!_compiled, "assembler has already been compiled!");
        const auto& trg = mesh_part.template get_target_set<shape_dim>();
        for(Index i(0); i < trg.get_num_entities(); ++i)
          this->_element_mask.at(trg[i]) = 1;
      }

      /**
       * \brief Returns a reference to the domain assembler's trafo.
       */
      const Trafo_& get_trafo() const
      {
        return this->_trafo;
      }

      /**
       * \brief Returns the element indices vector.
       */
      const std::vector<Index>& get_element_indices() const
      {
        return this->_element_indices;
      }

      /**
       * \brief Sets the maximum number of worker threads.
       *
       * \attention
       * The number of worker threads has to be set before the assembler is compiled.
       */
      void set_max_worker_threads(std::size_t max_worker_threads)
      {
        XASSERTM(!_compiled, "assembler has already been compiled!");
        this->_max_worker_threads = max_worker_threads;
      }

      /**
       * \brief Returns the maximum number of worker threads.
       */
      std::size_t get_max_worker_threads() const
      {
        return this->_max_worker_threads;
      }

      /**
       * \brief Returns the actual number of worker threads.
       *
       * This returned count corresponds to the actual number of worker threads used for assembling,
       * which may be lower than the maximum number of worker requested by the user.
       *
       * \note The assembler has to be compiled before calling this function, because the actual
       * number of worker threads to use is determined in the compilation process.
       */
      std::size_t get_num_worker_threads() const
      {
        return this->_num_worker_threads;
      }

      /**
       * \brief Sets the desired threading strategy
       *
       * See Assembly::ThreadingStrategy for details.
       *
       * \attention
       * The threading strategy has to be set before the assembler is compiled.
       */
      void set_threading_strategy(ThreadingStrategy strategy)
      {
        XASSERTM(!_compiled, "assembler has already been compiled!");
        this->_strategy = strategy;
      }

      /**
       * \brief Returns the threading strategy.
       */
      ThreadingStrategy get_threading_strategy() const
      {
        return this->_strategy;
      }

      /**
       * \brief Compiles the assembler for all elements that have been added manually.
       *
       * This function compiles the assembler for all cells/elements that have been added manually
       * by previous calls of the #add_element or #add_mesh_part functions according to the chosen
       * threading strategy and desired maximum worker thread count.
       *
       * \note
       * If you want to compile the assembler for all cells/elements of the trafo's
       * underlying mesh, consider using the #compile_all_elements() instead.
       */
      void compile()
      {
        XASSERTM(!_compiled, "assembler has already been compiled!");
        XASSERT(_element_indices.empty());

        // build element indices from element mask
        Index num_elems = this->_trafo.get_mesh().get_num_elements();
        this->_element_indices.reserve(num_elems);
        for(std::size_t i(0); i < this->_element_mask.size(); ++i)
        {
          if(this->_element_mask[i] != 0)
          {
            _element_indices.push_back(Index(i));
          }
        }

        this->_compile();
      }

      /**
       * \brief Compiles the assembler for all elements of the underlying mesh.
       *
       * This function compiles the assembler for all elements of the trafo's underlying mesh
       * according to the chosen threading strategy and desired maximum worker thread count.
       */
      void compile_all_elements()
      {
        XASSERTM(!_compiled, "assembler has already been compiled!");
        XASSERT(_element_indices.empty());

        // assemble on all elements
        Index num_elems = this->_trafo.get_mesh().get_num_elements();
        _element_indices.resize(num_elems);
        for(Index i(0); i < num_elems; ++i)
          _element_indices[i] = i;

        this->_compile();
      }

      /**
       * \brief Executes a domain assembly job (in parallel) by (multiple) worker threads.
       *
       * \param[inout] job
       * A \transient reference to the job that is to be assembled.
       * See Assembly::DomainAssemblyJob for details.
       */
      template<typename Job_>
      void assemble(Job_& job)
      {
        XASSERTM(_compiled, "assembler has not been compiled yet");
        XASSERTM(_threads.empty(), "already executing a job");

        // no worker threads?
        if(this->_num_worker_threads <= 0)
        {
          // assemble on master thread instead
          assemble_master(job);
          return;
        }

        // reset all fences
        for(auto& s : this->_thread_fences)
          s.close();

        // create worker threads
        for(std::size_t i(0); i < this->_num_worker_threads; ++i)
        {
          // create worker thread
          _threads.emplace_back(std::thread(Worker<Job_>(
            job, i+1, this->_num_worker_threads, this->_strategy,
            this->_thread_stats.at(i),
            this->_thread_mutex,
            this->_thread_fences,
            this->_element_indices,
            this->_color_elements,
            this->_layer_elements,
            this->_thread_layers
          )));
        }

        // assemble based on the chosen strategy
        switch(this->_strategy)
        {
        case ThreadingStrategy::single:
        case ThreadingStrategy::layered:
        case ThreadingStrategy::layered_sorted:
          // single/layered assembly is straight forward:
          // start all threads and wait for them to finish
          this->_thread_fences.front().open(true);
          for(std::size_t i(0); i < this->_threads.size(); ++i)
            this->_threads.at(i).join();
          break;

        case ThreadingStrategy::colored:
          // colored assembly is significantly more complex:
          // each layer represents a single color and all threads have
          // to traverse the layers simultaneously to avoid race conditions
          for(std::size_t icol(0); icol+1u < this->_color_elements.size(); ++icol)
          {
            // start threads by opening the front fence
            this->_thread_fences.front().open(true);

            bool all_okay = true;

            // wait for all threads to finish
            for(std::size_t i(0); i < this->_threads.size(); ++i)
            {
              all_okay = (this->_thread_fences.at(i+1u).wait() && all_okay);
              this->_thread_fences.at(i+1u).close();
            }

            // reset start signal
            this->_thread_fences.front().close();
            this->_thread_fences.back().open(all_okay);

            // break out in case of an error
            if(!all_okay)
              break;

            // wait for all threads to finish
            for(std::size_t i(0); i < this->_threads.size(); ++i)
            {
              this->_thread_fences.at(i+1u).wait();
              this->_thread_fences.at(i+1u).close();
            }
            this->_thread_fences.back().close();
          }

          // wait for all threads to finish
          for(std::size_t i(0); i < this->_threads.size(); ++i)
            this->_threads.at(i).join();
          break;

        default:
          XABORTM("invalid threading strategy!");
          break;
        }

        // clear thread vector
        this->_threads.clear();
      }

      /**
       * \brief Executes a domain assembly job directly on the calling thread.
       *
       * This function executes the assembly job directly on the calling thread without spawning
       * dedicated worker threads, which makes debugging of assembly jobs and tasks much easier.
       * Note that as a consequence, the assembly performed by this function is always
       * single-threaded even if the assembler has been configured to work with multiple threads.
       *
       * \param[inout] job
       * A \transient reference to the job that is to be assembled.
       * See Assembly::DomainAssemblyJob for details.
       */
      template<typename Job_>
      void assemble_master(Job_& job)
      {
        XASSERTM(_compiled, "assembler has not been compiled yet");
        XASSERTM(_threads.empty(), "already executing a job");

        // reset all fences
        for(auto& s : this->_thread_fences)
          s.close();

        // create worker object
        Worker<Job_> worker(job, 0, 0, this->_strategy, this->_thread_stats.front(),
          this->_thread_mutex, this->_thread_fences, this->_element_indices,
          this->_color_elements, this->_layer_elements, this->_thread_layers);

        // signal begin
        this->_thread_fences.front().open(true);
        this->_thread_fences.back().open(true);

        // perform work
        worker();
      }

      /**
       * \brief Returns a string dump of various debugging information
       *
       * \returns A string dump of various debugging information
       */
      String dump() const
      {
        std::ostringstream oss;

        oss << "Elements: " << stringify(this->_element_indices.size()) << " of " <<
          stringify(this->_trafo.get_mesh().get_num_elements()) << std::endl;

        oss << "Strategy: " << (this->_strategy == ThreadingStrategy::layered ? "layered" : "colored") << std::endl;

        if(_strategy == ThreadingStrategy::layered)
        {
          oss << std::endl << "Layers:" << std::endl;
          for(std::size_t i(0); i+1u < this->_layer_elements.size(); ++i)
          {
            auto jb = this->_layer_elements.at(i);
            auto je = this->_layer_elements.at(i+1);
            oss << stringify(i).pad_front(3) << ": "
              << stringify(je - jb).pad_front(6)
              << std::endl;
          }
           oss << std::endl;

          double desired = double(this->_element_indices.size()) / Math::max(double(this->_num_worker_threads), 1.0);
          oss << "Desired : " << stringify(desired) << std::endl;
          oss << "Threads:" << std::endl;
          for(std::size_t i(0); i+1 < this->_thread_layers.size(); ++i)
          {
            auto jb = this->_thread_layers.at(i);
            auto je = this->_thread_layers.at(i+1);
            auto nel = this->_layer_elements.at(je) - this->_layer_elements.at(jb);
            oss << stringify(i).pad_front(3) << ": " << stringify(je-jb).pad_front(4) << " : "
              << stringify(nel).pad_front(6) << " :" << stringify_fp_fix(double(nel) / desired, 3, 6)
              << std::endl;
          }
        }
        else
        {
          oss << std::endl << "Colors:" << std::endl;
          for(std::size_t i(0); i+1u < this->_color_elements.size(); ++i)
            oss << stringify(i).pad_front(3) << ": "
              << stringify(this->_color_elements.at(i+1) - this->_color_elements.at(i)).pad_front(6)
              << std::endl;
        }

        return oss.str();
      }

      /**
       * \brief Resets the thread statistics.
       */
      void reset_thread_stats()
      {
        for(auto& x : this->_thread_stats)
          x.reset();
      }

      /**
       * \brief Reduces the thread statistics to a single object.
       *
       * \returns A combined thread statistics object.
       */
      ThreadStats reduce_thread_stats() const
      {
        ThreadStats r;
        for(auto& x : this->_thread_stats)
          r += x;
        return r;
      }

    protected:
      /**
       * \brief Compiles the domain assembler.
       *
       * This function performs all the required initializations based on the chosen threading
       * strategy and desired worker thread count.
       */
      void _compile()
      {
        // nothing to do?
        if(this->_element_indices.empty())
        {
          _compiled = true;
          return;
        }

        // choose automatic strategy?
        if(this->_strategy == ThreadingStrategy::automatic)
        {
          // only 1 thread? => use single-threaded strategy
          if(this->_max_worker_threads <= std::size_t(1))
            this->_strategy = ThreadingStrategy::single;
          // multi-threaded: is the mesh permuted using colored strategy? => use colored threading strategy
          else if(this->_trafo.get_mesh().get_mesh_permutation().get_strategy() == Geometry::PermutationStrategy::colored)
            this->_strategy = ThreadingStrategy::colored;
          // multi-threaded: use layered threading strategy
          else
            this->_strategy = ThreadingStrategy::layered;
        }

        // do we want multi-threading?
        this->_num_worker_threads = 0;
        if(this->_max_worker_threads > 0)
        {
          // build the graphs for our element list
          this->_build_graphs();

          // note that one of the following function calls may set the
          // _num_worker_threads member variable to zero to disable
          // multi-threading if there are too few elements available and
          // therefore multi-threading is pointless or even impossible

          // build layers/colors for the threading strategy
          switch(this->_strategy)
          {
          case ThreadingStrategy::layered:
            this->_build_layers(false, false);
            this->_build_thread_layers();
            break;

          case ThreadingStrategy::layered_sorted:
            this->_build_layers(true, true);
            this->_build_thread_layers();
            break;

          case ThreadingStrategy::colored:
            this->_build_colors();
            break;

          default:
            // go no, there's nothing to see here...
            break;
          }
        }

        // we need one fence per worker thread and two additional fences
        // (the first and the last) for the master thread
        _thread_fences = std::vector<ThreadFence>(this->_num_worker_threads + 2u);
        _threads.reserve(std::size_t(this->_num_worker_threads));

        // stats are reserved for the maximum desired number of threads
        _thread_stats = std::vector<ThreadStats>(this->_max_worker_threads + 1u);

        _compiled = true;
      }

      /**
       * \brief Builds the element adjacencies graphs
       */
      void _build_graphs()
      {
        // get vertices-at-element index set
        const auto& idx_set = this->_trafo.get_mesh().template get_index_set<shape_dim,0>();

        // query dimensions
        const Index nel = Index(this->_element_indices.size());
        const Index nvt = idx_set.get_index_bound();
        const Index nix = Index(idx_set.get_num_indices());

        // allocate graph
        this->_verts_at_elem = Adjacency::Graph(nel, nvt, nel*nix);
        Index* dom_ptr = this->_verts_at_elem.get_domain_ptr();
        Index* img_idx = this->_verts_at_elem.get_image_idx();

        // build domain pointer array
        for(Index i(0); i <= nel; ++i)
          dom_ptr[i] = nix * i;

        // build image index array
        for(Index i(0), k(0); i < nel; ++i)
        {
          const auto& idx = idx_set[this->_element_indices.at(i)];
          for(Index j(0); j < nix; ++j, ++k)
            img_idx[k] = idx[int(j)];
        }

        // sort indices
        this->_verts_at_elem.sort_indices();

        // transpose graph
        this->_elems_at_vert = Adjacency::Graph(Adjacency::RenderType::transpose, this->_verts_at_elem);

        // build neighbors graph
        this->_elem_neighbors = Adjacency::Graph(Adjacency::RenderType::injectify_sorted, this->_verts_at_elem, this->_elems_at_vert);
      }

      /**
       * \brief Builds the Cuthill-McKee layer graphs
       *
       * \param[in] reverse
       * Specifies whether to reverse the ordering.
       *
       * \param[in] sorted
       * Specifies whether to sort entities in each layer by their degree.
       *
       * \todo utilize mesh permutation layering if available
       */
      void _build_layers(bool reverse, bool sorted)
      {
        // get number of elements and vertices
        const Index num_elems = this->_elem_neighbors.get_num_nodes_domain();

        // get neighbor arrays
        const Index* neigh_ptr = this->_elem_neighbors.get_domain_ptr();
        const Index* neigh_idx = this->_elem_neighbors.get_image_idx();

        // a compare function for std::stable_sort
        DegreeCompare degree_compare(this->_elem_neighbors);

        // allocate an element mask vector and initialize to 0
        std::vector<int> elem_mask(num_elems, 0);

        // allocate a new element vector
        std::vector<Index> elements, layers;
        elements.reserve(num_elems);
        layers.reserve(num_elems);

        // push beginning of first layer
        layers.push_back(0u);

        // main Cuthill-McKee loop
        while(Index(elements.size()) < num_elems)
        {
          // pick a new root element
          Index root = num_elems + 1u;

          // choose the node of minimum degree
          {
            Index min = num_elems + 1u;
            for(Index j(0); j < num_elems; ++j)
            {
              Index deg = this->_elem_neighbors.degree(j);
              if((deg < min) && (elem_mask[j] == 0))
              {
                root = j;
                min = deg;
              }
            }
          }

          XASSERTM(root < num_elems, "no valid root element found");

          // push next layer
          elements.push_back(root);
          elem_mask[root] = int(layers.size());

          // loop over the adjacency levels of the root element
          while(Index(elements.size()) < num_elems)
          {
            // get layer start
            const Index layer_beg = layers.back();
            const Index layer_end = Index(elements.size());
            layers.push_back(layer_end);

            // loop over all elements in the current layer
            for(Index i(layer_beg); i < layer_end; ++i)
            {
              // get the element's index
              const Index elem_idx = elements.at(i);

              // loop over all element neighbors
              for(Index j(neigh_ptr[elem_idx]); j < neigh_ptr[elem_idx+1]; ++j)
              {
                // get the neighbor's element index
                const Index elem_jdx = neigh_idx[j];

                // did we already process this element?
                if(elem_mask[elem_jdx] == 0)
                {
                  // add element
                  elements.push_back(elem_jdx);
                  elem_mask[elem_jdx] = int(layers.size());
                }
              }
            }

            // no new elements?
            const Index layer_nxt = Index(elements.size());
            if(layer_nxt <= layer_end)
              break;

            // sort elements in layer by neighbor degree
            if(sorted)
              std::stable_sort(elements.begin() + std::ptrdiff_t(layer_end), elements.begin() + std::ptrdiff_t(layer_nxt), degree_compare);

            // continue with next layer
          }
          // continue with next root
        }

        // push final layer end
        const Index num_layers = Index(layers.size());
        layers.push_back(num_elems);

        // translate and sort element indices
        for(Index ilay(0); ilay < num_layers; ++ilay)
        {
          Index ibeg = layers.at(ilay);
          Index iend = layers.at(ilay+1u);

          for(Index j(ibeg); j < iend; ++j)
            elements[j] = this->_element_indices[elements[j]];

          // sort elements in layer
          if(sorted)
            std::sort(elements.begin() + std::ptrdiff_t(ibeg), elements.begin() + std::ptrdiff_t(iend));
        }

        // reverse layers?
        if(reverse)
        {
          this->_element_indices.clear();
          this->_layer_elements.clear();
          this->_layer_elements.reserve(layers.size());
          this->_layer_elements.push_back(0u);

          // reverse layers
          for(Index ilay(0); ilay < num_layers; ++ilay)
          {
            Index ibeg = layers.at(num_layers-ilay-1u);
            Index iend = layers.at(num_layers-ilay);

            for(Index j(ibeg); j < iend; ++j)
            {
              this->_element_indices.push_back(elements[j]);
            }

            this->_layer_elements.push_back(this->_layer_elements.back() + iend - ibeg);
          }
        }
        else
        {
          this->_element_indices = std::move(elements);
          this->_layer_elements = std::move(layers);
        }
      }

      /**
       * \brief Build the actual thread layers for the layered strategy
       */
      bool _build_thread_layers()
      {
        // no threading?
        if(this->_max_worker_threads < std::size_t(1))
          return false;

        // set the number of actual worker threads; we want at least 3 layers per thread on average
        this->_num_worker_threads = Math::min(this->_max_worker_threads,
          this->_layer_elements.size() / std::size_t(3));

        const Index num_elems = Index(this->_element_indices.size());
        const Index num_layers = Index(this->_layer_elements.size() - 1u);

        this->_thread_layers.reserve(this->_num_worker_threads+1u);
        this->_thread_layers.push_back(0u);

        // In the following, we need to assign a set of consecutive layers to each thread.
        // We want to distribute the layers so that each thread is assigned roughly the
        // same number elements to avoid imbalance, but at the same time we have to make
        // sure that each thread is assigned at least two layers to rule out race conditions.
        // This is performed in a 3-step approach: first we distribute the layers to avoid
        // imbalance and afterwards then we enforce the minimum of two layers per thread
        // by a backward and a forward sweep. Note that in general this approach will not
        // yield an optimal partitioning if there are only few layers, but this shouldn't
        // really matter in practice.

        // Example:
        // Assume we have 60 elements, 8 layers of different sizes and 3 worker threads
        // In a perfect world, each worker thread should process the same number of
        // elements (20), but in our case, the threads can only be assigned whole layers.
        //
        // Step 1: for each thread, we choose the starting layer whose first element
        //         is greater or equal to the desired first element for that thread
        // Step 2: ensure that each thread has at least 2 layer by a backward sweep
        // Step 3: ensure that each thread has at least 2 layer by a forward sweep
        //
        // elements:     .123456789.123456789.123456789.123456789.123456789.123456789
        // layers:       |0-|1--|2---|3----|4-----|5-------|6---------|7------------|
        // desired:      |0------------------|1------------------|2-----------------|
        // threads #1:   |0-----------------------|1------------------|2------------|
        // threads #2:   |0----------------|1--------------|2-----------------------|
        // threads #3:   |0----------------|1--------------|2-----------------------|

        // Step 1: build thread layers by desired elements
        for(std::size_t i(1); i < this->_num_worker_threads; ++i)
        {
          // compute the desired first element of the i-th thread
          const Index desired_first = (num_elems * Index(i)) / Index(this->_num_worker_threads);

          // choose the first layer whose first element is greater or equal to our desired element
          Index j(this->_thread_layers.back()+1);
          while((j < num_layers) && (this->_layer_elements.at(j) < desired_first))
            ++j;
          this->_thread_layers.push_back(j);
        }
        this->_thread_layers.push_back(num_layers);

        // Step 2: make sure each thread has at least two layers by backward sweep
        for(std::size_t i(this->_num_worker_threads-1); i > 0; --i)
        {
          // make sure there are at least two layers for this thread
          // if not, then decrease the preceding thread's starting layer index
          if(this->_thread_layers.at(i+1) < this->_thread_layers.at(i) + Index(2))
            this->_thread_layers.at(i) = this->_thread_layers.at(i+1) - Index(2);

          // bail out if we have less that two layers left;
          // the next loop takes care of this case
          if(this->_thread_layers.at(i) < Index(2))
            break;
        }

        // Step 3: make sure each thread has at least two layers by forward sweep
        for(std::size_t i(0); i < this->_num_worker_threads; ++i)
        {
          // make sure there are at least two layers for this thread
          // if not, then increase the succeeding thread's starting layer index
          if(this->_thread_layers.at(i+1) < this->_thread_layers.at(i) + Index(2))
            this->_thread_layers.at(i+1) = this->_thread_layers.at(i) + Index(2);
        }

        // make sure that we didn't change the first and last entries
        XASSERT(this->_thread_layers.front() == Index(0));
        XASSERT(this->_thread_layers.back() == num_layers);

        // okay
        return true;
      }

      /**
       * \brief Builds the color element vectors for the colored threading strategy.
       *
       * \todo utilize mesh permutation coloring if available
       */
      void _build_colors()
      {
        // create coloring from our neighbors graph
        Adjacency::Coloring coloring(this->_elem_neighbors);

        // create the partitioning graph from our coloring
        Adjacency::Graph color_parti = coloring.create_partition_graph();

        const Index num_elems = color_parti.get_num_nodes_image();
        const Index num_colors = color_parti.get_num_nodes_domain();
        const Index* dom_ptr = color_parti.get_domain_ptr();
        const Index* img_idx = color_parti.get_image_idx();

        // sanity check
        XASSERT(num_elems == Index(this->_element_indices.size()));

        // set number of worker threads to the minimum of the desired maximum number of threads
        // and the maximum number of elements per color; note that it is perfectly legal if some
        // colors contain less elements than we have worker threads because in this case some of
        // the threads will simply twiddle their thumbs during the assembly
        this->_num_worker_threads = Math::min(this->_max_worker_threads, std::size_t(color_parti.degree()));

        // backup element indices
        std::vector<Index> elems(this->_element_indices);

        // translate image indices
        for(Index i(0); i < num_elems; ++i)
          this->_element_indices.at(i) = elems.at(img_idx[i]);

        // store color layers
        for(Index i(0); i <= num_colors; ++i)
          this->_color_elements.push_back(dom_ptr[i]);
      }
    }; // DomainAssembler
  } // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_DOMAIN_ASSEMBLER_HPP
