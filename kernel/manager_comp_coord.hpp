/* GENERAL_REMARK_BY_HILMAR:
 * See GENERAL_REMARK_BY_HILMAR in manager.hpp.
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
#pragma once
#ifndef KERNEL_MANAGER_COMP_COORD_HPP
#define KERNEL_MANAGER_COMP_COORD_HPP 1

// includes, Feast
#include <kernel/base_header.hpp>
#ifdef PARALLEL
#include <kernel/util/exception.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/error_handler.hpp>
#include <kernel/manager_comp.hpp>
#include <kernel/load_balancer.hpp>
#include <kernel/base_mesh/file_parser.hpp>
#include <kernel/base_mesh/bm.hpp>

// includes, system
#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <vector>

namespace FEAST
{

  /// forward declaration of class Manager (needed for friend declarations)
  template<
    unsigned char space_dim_,
    unsigned char world_dim_>
  class Manager;

  /**
  * \brief defines manager of a compute process group on the coordinator process
  *
  * See the description of class Manager.
  *
  * \author Hilmar Wobker
  */
  template<
    unsigned char space_dim_,
    unsigned char world_dim_>
  class ManagerCompCoord
    : public ManagerComp<space_dim_, world_dim_>
  {

    /// declare Manager function as friend to grant access to function _create_work_groups()
    friend void Manager<space_dim_, world_dim_>::create_work_groups();
    friend unsigned int Manager<space_dim_, world_dim_>::test_communication();

  private:

    /* *****************
    * member variables *
    *******************/
    /**
    * \brief vector of graph structures representing the process topology within the work groups
    *
    * These are only known to the coordinator process. The worker processes only get 'their' portions of the global
    * graph.
    */
    std::vector<Graph*> _graphs;

    /// base mesh the manager works with
    BaseMesh::BM<space_dim_, world_dim_>* _base_mesh;

    /// pointer to the load balancer of the process group
    LoadBalancer<space_dim_, world_dim_>* _load_balancer;


    /* **********
    * shortcuts *
    ************/
    /**
    * \brief getter for the vector of extended work groups
    *
    * Just a shortcut to save typing of 'ManagerComp<space_dim_, world_dim_>::'
    *
    * \return reference to the vector of extended work groups
    */
    inline const std::vector<WorkGroupExt*>& _work_groups() const
    {
      CONTEXT("ManagerCompCoord::_work_groups()");
      return ManagerComp<space_dim_, world_dim_>::_work_groups;
    }


    /**
    * \brief getter for the number of work groups
    *
    * Just a shortcut to save typing of 'ManagerComp<space_dim_, world_dim_>::'
    *
    * \return number of work groups
    */
    inline unsigned int _num_work_groups() const
    {
      CONTEXT("ManagerCompCoord::_num_work_groups()");
      return ManagerComp<space_dim_, world_dim_>::_num_work_groups;
    }


    /**
    * \brief getter for the array indicating to which work groups this process belongs
    *
    * Just a shortcut to save typing of 'ManagerComp<space_dim_, world_dim_>::'
    *
    * \return array indicating to which work groups this process belongs
    */
    inline bool* _belongs_to_work_group() const
    {
      CONTEXT("ManagerCompCoord::_belongs_to_work_group()");
      return ManagerComp<space_dim_, world_dim_>::_belongs_to_work_group;
    }


    /**
    * \brief getter for the number of interlevel groups
    *
    * Just a shortcut to save typing of 'ManagerComp<space_dim_, world_dim_>::'
    *
    * \return number of interlevel groups
    */
    inline unsigned int _num_interlevel_groups() const
    {
      CONTEXT("ManagerCompCoord::_num_interlevel_groups()");
      return ManagerComp<space_dim_, world_dim_>::_num_interlevel_groups;
    }


    /**
    * \brief getter for the process group
    *
    * Just a shortcut to save typing of 'ManagerComp<space_dim_, world_dim_>::'
    *
    * \return ProcessGroup pointer #_process_group
    */
    inline ProcessGroup* _process_group() const
    {
      CONTEXT("ManagerCompCoord::_process_group()");
      return ManagerComp<space_dim_, world_dim_>::_process_group;
    }


    /* *****************
    * member functions *
    *******************/
    /**
    * \brief function that sets up (extended) work groups and interlevel groups basing on the provided information
    *
    * This function is called on the coordinator of the manager's compute process group.
    *
    * \author Hilmar Wobker
    */
    void _create_work_groups()
    {
      CONTEXT("ManagerCompCoord::_create_work_groups()");

      // set data/pointers on the coordinator process (where the data is already available)
      ManagerComp<space_dim_, world_dim_>::_num_work_groups = _load_balancer->num_work_groups();
      ManagerComp<space_dim_, world_dim_>::_num_proc_in_ext_work_group = _load_balancer->num_proc_in_ext_work_group();
      ManagerComp<space_dim_, world_dim_>::_group_contains_extra_coord = _load_balancer->group_contains_extra_coord();
      ManagerComp<space_dim_, world_dim_>::_ext_work_group_ranks = _load_balancer->ext_work_group_ranks();

      // now the coordinator broadcasts the relevant data to the other processes, that is:
      //   - _num_work_groups
      //   - _num_proc_in_ext_work_group
      //   - _group_contains_extra_coord

      unsigned int temp(_num_work_groups());
      int mpi_error_code = MPI_Bcast(&temp, 1, MPI_UNSIGNED, _process_group()->rank_coord(), _process_group()->comm());
      validate_error_code_mpi(mpi_error_code, "MPI_Bcast");

      mpi_error_code = MPI_Bcast(ManagerComp<space_dim_, world_dim_>::_num_proc_in_ext_work_group, _num_work_groups(),
                                 MPI_UNSIGNED, _process_group()->rank_coord(), _process_group()->comm());
      validate_error_code_mpi(mpi_error_code, "MPI_Bcast");

      mpi_error_code = MPI_Bcast(ManagerComp<space_dim_, world_dim_>::_group_contains_extra_coord, _num_work_groups(),
                                 MPI_UNSIGNED_CHAR, _process_group()->rank_coord(), _process_group()->comm());
      validate_error_code_mpi(mpi_error_code, "MPI_Bcast");

      // set data/pointers on the coordinator process (where the data is already available)
      ManagerComp<space_dim_, world_dim_>::_num_interlevel_groups = _load_balancer->num_interlevel_groups();
      ManagerComp<space_dim_, world_dim_>::_num_proc_in_interlevel_group
        = _load_balancer->num_proc_in_interlevel_group();
      ManagerComp<space_dim_, world_dim_>::_interlevel_group_ranks = _load_balancer->interlevel_group_ranks();
      ManagerComp<space_dim_, world_dim_>::_pos_of_root_in_interlevel_group
        = _load_balancer->pos_of_root_in_interlevel_group();

      // now the coordinator broadcasts the relevant data to the other processes, that is:
      //   - _num_interlevel_groups
      //   - _num_proc_in_interlevel_group
      //   - _pos_of_root_in_interlevel_group

      temp = _num_interlevel_groups();
      mpi_error_code = MPI_Bcast(&temp, 1, MPI_UNSIGNED, _process_group()->rank_coord(), _process_group()->comm());
      validate_error_code_mpi(mpi_error_code, "MPI_Bcast");

      if(_num_interlevel_groups() > 0)
      {
        mpi_error_code = MPI_Bcast(ManagerComp<space_dim_, world_dim_>::_num_proc_in_interlevel_group,
                                   _num_interlevel_groups(), MPI_UNSIGNED, _process_group()->rank_coord(),
                                   _process_group()->comm());
        validate_error_code_mpi(mpi_error_code, "MPI_Bcast");

        mpi_error_code = MPI_Bcast(ManagerComp<space_dim_, world_dim_>::_pos_of_root_in_interlevel_group,
                                   _num_interlevel_groups(), MPI_UNSIGNED, _process_group()->rank_coord(),
                                   _process_group()->comm());
        validate_error_code_mpi(mpi_error_code, "MPI_Bcast");
      }

      // call routine for creating work groups (within this routine the arrays _ext_work_group_ranks and
      // _interlevel_group_ranks are transferred)
      ManagerComp<space_dim_, world_dim_>::_create_work_groups_common();
    } // _create_work_groups()


    /**
    * \brief sends the relevant parts of the global graph to the corresponding workers of each work group
    *
    * The local graphs tell the workers with which workers of their work group they have to communicate.
    * This function must be called on the coordinator process of the process group, at the same time all
    * non-coordinator processes must call the function receive_and_set_graphs(). Before this function can be used, the
    * function create_work_groups() must have been called.
    *
    * \param[in] graphs
    * array of Graph pointers representing the connectivity of work group processes
    *
    * \sa receive_and_set_graphs()
    *
    * \author Hilmar Wobker
    */
    void _transfer_graphs_to_workers(Graph** graphs)
    {
      CONTEXT("ManagerCompCoord::_transfer_graphs_to_workers()");

      for(unsigned int igroup(0) ; igroup < _num_work_groups() ; ++igroup)
      {
        if(_belongs_to_work_group()[igroup])
        {
          ASSERT(_work_groups()[igroup]->is_coordinator(), "Routine must be called on the coordinator process.");
          // number of neighbours of the graph node corresponding to this process (use unsigned int datatype here
          // instead of index_glob_t since MPI routines expect it)
          unsigned int num_neighbours_local;
          index_glob_t* neighbours_local(nullptr);
          int rank_coord = _work_groups()[igroup]->rank_coord();

          // set the graph pointer for the current work group
          _graphs.push_back(graphs[igroup]);

          // since the MPI routines used below expect integer arrays, we have to copy two index_glob_t arrays
          // within the graph structures to corresponding int arrays
// COMMENT_HILMAR: Is there a possibility to avoid this? A reinterpret_cast<int*>(unsigned long) does not work!
          unsigned int* num_neighbours_aux;
          unsigned int* index_aux;

          index_glob_t num_nodes;
          if(_work_groups()[igroup]->contains_extra_coordinator())
          {
            // In case there is an extra coordinator process, we have to add one pseudo node to the graph and the
            // index array must be modified correspondingingly. This pseudo node corresponds to the coordinator
            // process itself which has to be included in the call of MPI_Scatterv(...). It has to appear at the
            // position in the arrays num_neighbours_aux[] and index_aux[] corresponding to the rank of the
            // coordinator. Although we know, that this is rank 0, we do not explicitly exploit this information here
            // since this might be changed in future.
            num_nodes = _graphs[igroup]->num_nodes() + 1;
            index_aux = new unsigned int[num_nodes + 1];
            // copy the first part of the graphs's index array to the aux array, performing implicit cast from
            // index_glob_t to unsigned int
            for(index_glob_t i(0) ; i < (index_glob_t)rank_coord+1 ; ++i)
            {
              index_aux[i] = _graphs[igroup]->index()[i];
            }
            // insert the pseudo node
            index_aux[rank_coord+1] = index_aux[rank_coord];
            // copy the remaining part of the graphs's index array to the aux array, performing implicit cast from
            // index_glob_t to unsigned int
            for(index_glob_t i(rank_coord+1) ; i < num_nodes ; ++i)
            {
              index_aux[i+1] = _graphs[igroup]->index()[i];
            }
          }
          else
          {
            // in case there is no extra coordinator process, the number of neighbours equals the number of nodes
            // in the graph and the index array does not have to be modified
            num_nodes = _graphs[igroup]->num_nodes();
            index_aux = new unsigned int[num_nodes + 1];
            // copy the graphs's index array to the aux array, performing implicit cast from index_glob_t to
            // unsigned int
            for(index_glob_t i(0) ; i < num_nodes+1 ; ++i)
            {
              index_aux[i] = _graphs[igroup]->index()[i];
            }
          }

          // now determine the number of neighbours per node (eventually including the pseudo node for the extra
          // coordinator process)
          num_neighbours_aux = new unsigned int[num_nodes];
          for(index_glob_t i(0) ; i < num_nodes ; ++i)
          {
            num_neighbours_aux[i] = index_aux[i+1] - index_aux[i];
          }

          if(_work_groups()[igroup]->contains_extra_coordinator())
          {
            // send the number of neighbours to the non-coordinator processes (use MPI_IN_PLACE to indicate that the
            // coordinator does not receive/store any data)
            MPI_Scatter(num_neighbours_aux, 1, MPI_UNSIGNED, MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                        rank_coord, _work_groups()[igroup]->comm());
            // send the neighbours to the non-coordinator processes
            MPI_Scatterv(_graphs[igroup]->neighbours(), reinterpret_cast<int*>(num_neighbours_aux),
                         reinterpret_cast<int*>(index_aux), MPIType<index_glob_t>::value(), MPI_IN_PLACE, 0,
                         MPI_DATATYPE_NULL, rank_coord, _work_groups()[igroup]->comm());
          }
          else
          {
            // When there is no extra coordinator process, then the coordinator is part of the compute work group and
            // also sends data to itself.

            // scatter the number of neighbours to the non-coordinator processes and to the coordinator process itself
            MPI_Scatter(reinterpret_cast<int*>(num_neighbours_aux), 1, MPI_UNSIGNED, &num_neighbours_local, 1,
                        MPI_INTEGER, rank_coord, _work_groups()[igroup]->comm());
            neighbours_local = new index_glob_t[num_neighbours_local];
            // scatter the neighbours to the non-coordinator processes and to the coordinator process itself
            MPI_Scatterv(_graphs[igroup]->neighbours(), reinterpret_cast<int*>(num_neighbours_aux),
                         reinterpret_cast<int*>(index_aux), MPIType<index_glob_t>::value(), neighbours_local,
                         num_neighbours_local, MPIType<index_glob_t>::value(), rank_coord,
                         _work_groups()[igroup]->comm());
          }
          // delete aux. arrays again
          delete [] num_neighbours_aux;
          delete [] index_aux;

          // now create distributed graph structure within the compute work groups. The coordinater only performs this
          // task when it is not an extra coordinator process, i.e., only if it actually is a worker process.
          if (!_work_groups()[igroup]->contains_extra_coordinator())
          {
            _work_groups()[igroup]->work_group()->set_graph_distributed(num_neighbours_local, neighbours_local);
            // The array neighbours_local is copied inside the constructor of the distributed graph object, hence it
            // can be deallocated again.
            delete [] neighbours_local;
          }
        } // if(_belongs_to_work_group()[igroup])
      } // for(unsigned int igroup(0) ; igroup < _num_work_groups() ; ++igroup)
    } // _transfer_graphs_to_workers()


    /**
    * \brief dummy function for testing work group and interlevel group communication
    *
    * The function simply performs an exchange of integer arrays between the corresponding processes.
    * Note that this function and the corresponding version in ManagerCompNonCoord are slightly different, which is
    * why they are not implemented in ManagerComp.
    *
    * \return flag whether test was succesful (0: succesful, >0: failed)
    */
    unsigned int _test_communication()
    {
      CONTEXT("ManagerCompCoord::_test_communication()");

      unsigned int work_group_test_ok = 0;
      unsigned int interlevel_group_test_ok = 0;
      // test local neighbourhood communication within work groups
      for(unsigned int igroup(0) ; igroup < _num_work_groups() ; ++igroup)
      {
        if (_belongs_to_work_group()[igroup] && !_work_groups()[igroup]->contains_extra_coordinator())
        {
          work_group_test_ok = work_group_test_ok + _work_groups()[igroup]->work_group()->test_communication();
        }
      }

      // test communication within interlevel groups
      for(unsigned int igroup(0) ; igroup < _num_interlevel_groups() ; ++igroup)
      {
        if (ManagerComp<space_dim_, world_dim_>::_belongs_to_interlevel_group[igroup])
        {
          interlevel_group_test_ok = interlevel_group_test_ok
            + ManagerComp<space_dim_, world_dim_>::_interlevel_groups[igroup]->test_communication();
        }
      }
      return work_group_test_ok + interlevel_group_test_ok;
    } // _test_communication



  public:


    /* *************************
    * constructor & destructor *
    ***************************/
    /// CTOR
    ManagerCompCoord(ProcessGroup* process_group)
      : ManagerComp<space_dim_, world_dim_>(process_group),
        _base_mesh(nullptr),
        _load_balancer(nullptr)
    {
      CONTEXT("ManagerCompCoord::ManagerCompCoord()");
    }

    /// DTOR
    ~ManagerCompCoord()
    {
      CONTEXT("ManagerCompCoord::~ManagerCompCoord()");
      if (_base_mesh != nullptr)
      {
        delete _base_mesh;
      }

      // assume that the _graphs vector only holds copies of the graph object pointers and that the graph objects
      // themselves are destroyed where they were created
      _graphs.clear();

    }


    /* ******************
    * getters & setters *
    ********************/
    /**
    * \brief getter for the base mesh
    *
    * \return BaseMesh::BM<space_dim_, world_dim_> pointer #_base_mesh
    */
    inline BaseMesh::BM<space_dim_, world_dim_>* base_mesh() const
    {
      CONTEXT("ManagerCompCoord::base_mesh()");
      return _base_mesh;
    }


    /**
    * \brief setter for the pointer to the load balancer
    *
    * \param[in] pointer to load balancer object
    */
    inline void set_load_balancer(LoadBalancer<space_dim_, world_dim_>* load_bal)
    {
      CONTEXT("ManagerCompCoord::set_load_balancer()");
      _load_balancer = load_bal;
    }


    /* *****************
    * member functions *
    *******************/
    /// read in a mesh file and set up base mesh
    void read_mesh(std::string const & mesh_file)
    {
      CONTEXT("ManagerCompCoord::read_mesh()");
      // the mesh is read by the process group coordinator

      _base_mesh = new BaseMesh::BM<space_dim_, world_dim_>();
      BaseMesh::FileParser<space_dim_, world_dim_> parser;
      Logger::log("Reading mesh file " + mesh_file + "...\n", Logger::master);
      try
      {
        parser.parse(mesh_file, _base_mesh);
      }
      catch(Exception& e)
      {
        // abort the program
        ErrorHandler::exception_occured(e);
      }
      // set cell numbers (equal to indices since all cells are active)
      _base_mesh->set_cell_numbers();
      // create base mesh's graph structure
      _base_mesh->create_graph();
      // print base mesh
      std::string s = _base_mesh->print();
      Logger::log(s, Logger::master_standard);
      Logger::log(s);
      // validate base mesh
      std::ostringstream strstream;
      _base_mesh->validate(strstream);
      Logger::log(strstream.str());
    }
  }; // class ManagerCompCoord
} // namespace FEAST

#endif // PARALLEL
#endif // KERNEL_MANAGER_COMP_COORD_HPP
