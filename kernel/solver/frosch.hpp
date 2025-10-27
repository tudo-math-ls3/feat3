#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/solver/adp_solver_base.hpp>

#if defined(FEAT_HAVE_TRILINOS) || defined(DOXYGEN)

#ifndef FEAT_HAVE_MPI
#  error FROSch without MPI does not make sense.
#endif

// includes, system
#include <vector>
#include <regex>

namespace FEAT
{
  namespace Solver
  {
    /// \cond internal
    namespace Trilinos
    {
      /**
       * \brief Interface-Class for a Teuchos::ParameterList object
       *
       * The "core wrapper object", which is returned by this function, Trilinos Teuchos::ParameterList object,
       * which is required for the configuration of the FROSch preconditioner
       *
       * \param[in] comm_
       * A Dist::Comm reference.  MUST be logically the same as in the matrix which will used for the preconditioner
       *
       * \param[in] problem_
       * FROSchParameterList::ProblemType [PRESSUREPOISSON | SADDLEPOINT].  Depending on the given argument slightly differnt
       * parameter lists will be constructed.  __ATTENTION:__ For the pressure Poisson:  The configuration is still done by the
       * pressure parameters.
       *
       * \param[in] nlevels_
       * Number of levels (Domain Decomposition levels) before the coarse problem is solved
       *     For example:
       *         The two-level preconditioner has _nlevels = 1
       *         The three-level preconditioner has _nlevels = 2
       *
       * \param[in] name
       * Name of the constructed parameterlist.
       *
       * \author Stephan Köhler
       */
      class FROSchParameterList
      {
        public:
          /* Use MUMPS or UMFPACK or ILU */
          enum DirectSolver {
            KLU,          /* Trilinos internal solver:  works always */
            MUMPS,        /* Trilinos needs to be compiled with Amesos2 and Mumps */
            UMFPACK,      /* Trilinos needs to be compiled with Amesos2 and Umfpack */
            ILU           /* Trilinos needs to be compiled with Ifpack2 */
          };
          /* Normally, use TWOLEVELBLOCK */
          enum Preconditioner {
            ONELEVEL,          /*  one-level overlapping schwarz; just for testing */
            TWOLEVEL,          /* works only for pressure Poisson and two-level neither for multi-level pressure Poisson nor for the saddle point problem */
            TWOLEVELBLOCK      /* works for pressure Poisson (two-level and multi-level) and the saddle point problem (two-level and multi-level) */
          };

          enum ProblemType {
            PRESSUREPOISSON,
            SADDLEPOINT
          };

          enum PartitionType {
            PHG,       /* Parallel hyper-graph, Zoltan */
            PARMETIS,  /* Trilinos needs to be compiled with Parmetis */
            BLOCK,     /* Should work always (not good) */
            SCOTCH,    /* Scotch */
            PTSCOTCH   /* Scotch */
          };

          enum PartitionApproach {
            PARTITION,
            REPARTITION    /* probably the best, but I'm not sure */
          };

          enum CombineOverlap {
            FULL,          /* produces a symmetric preconditioner */
            RESTRICTED,    /* __ATTENTION:__ produces a non-symmetric preconditioner.  As default: Restricted with (F)GMRES should work better than FULL with PCG. */
            AVERAGING      /* __ATTENTION:__ produces a non-symmetric preconditioner.  Should work better than FULL with PCG. */
          };

          enum CombineLevel {
            ADDITIVE,         /* Combine first and second level in an additive way */
            MULTIPLICATIVE    /* Combine first and second level in a multiplicative way: __ATTENTION:__ this leads to an unsymmetric preconditioner */
          };

          /* interface partition of unity
           *
           * Normally:
           *     velocity:
           *         In 2D use GDSW.
           *         In 3D use GDSWSTAR.
           *     pressure:
           *         Use GDSW, but for discont. pressure it does not matter which is choosen
           *
           */
          enum IPOU {
            GDSW,
            GDSWSTAR,   /* Only for three dimensional problems */
            RGDSW
          };

          /**
           * \brief Interface-Class for a Teuchos::ParameterList object
           *
           * The "core wrapper object", which is returned by this function, Trilinos Teuchos::ParameterList object,
           * which is required for the configuration of the FROSch preconditioner
           *
           * \param[in] comm_
           * A Dist::Comm reference.  MUST be logically the same as in the matrix which will used for the preconditioner
           *
           * \param[in] problem_
           * FROSchParameterList::ProblemType [PRESSUREPOISSON | SADDLEPOINT].  Depending on the given argument slightly differnt
           * parameter lists will be constructed.  __ATTENTION:__ For the pressure Poisson:  The configuration is still done by the
           * pressure parameters.
           *
           * \param[in] nlevels_
           * Number of levels (Domain Decomposition levels) before the coarse problem is solved
           *     For example:
           *         The two-level preconditioner has _nlevels = 1
           *         The three-level preconditioner has _nlevels = 2
           *
           * \param[in] name
           * Name of the constructed parameterlist.
           *
           * \author Stephan Köhler
           */
          explicit FROSchParameterList(const Dist::Comm& comm_, const int dim_, const FROSchParameterList::ProblemType problem_, const int nlevels_ = 1, const std::string& name_ = "FEAT-FROSch Parameterlist");

          /**
           * \brief Interface-Class for a Teuchos::ParameterList object
           *
           * The "core wrapper object", which is returned by this function, Trilinos Teuchos::ParameterList object,
           * which is required for the configuration of the FROSch preconditioner
           *
           * \param[in] comm_
           * A Dist::Comm reference.  MUST be logically the same as in the matrix which will used for the preconditioner
           *
           * \param[in] xml_file_
           * A file path to a Trilinos parameter list XML file.  The file will be read with DistFileIO
           *     For example:
           *         The two-level preconditioner has _nlevels = 1
           *         The three-level preconditioner has _nlevels = 2
           *
           * \param[in] name
           * Name of the constructed parameterlist.
           *
           * \author Stephan Köhler
           */
          explicit FROSchParameterList(const Dist::Comm& comm_, const int dim_, const std::string& xml_file_, const std::string& name_ = "Parameter List");
          ~FROSchParameterList();

          void read_from_xml_str(const std::string &xml_str_);
          void read_from_xml_file(const std::string &xml_file_);

          /* Print parameter list to rank zero */
          void print() const;

          /* return a pointer to the core */
          void* get_core() const;

          void create_core();
          void delete_core();

          int get_dim() const;
          int get_dpe_velo() const;
          int get_dpe_pres() const;
          int get_nlevels() const;
          bool get_use_timer() const;
          bool get_print() const;

          void set_dim(const int dim);
          void set_nlevels(const int nlevels);
          void set_problem_type(const ProblemType problem_);
          void set_use_timer(const bool use_timer_);
          void set_print(const bool print_);

          void set_precond(const Preconditioner precond_);
          void set_coarse_solver(const DirectSolver solver_);

          void set_gsteps(const int gsteps_);
          void set_gsteps(const std::vector<int>& gsteps_);

          void set_overlaps(const int overlap_);
          void set_overlaps(const std::vector<int>& overlaps_);

          void set_combine_overlap(const CombineOverlap combine_mode_);
          void set_combine_overlap(const std::vector<CombineOverlap>& combine_mode_);

          void set_combine_lvl(const CombineLevel combine_mode_);
          void set_combine_lvl(const std::vector<CombineLevel>& combine_lvl_);

          void set_solvers(const DirectSolver solver_ao_, const DirectSolver solver_ext_);
          void set_solvers(const std::vector<DirectSolver>& solvers_ao_, const std::vector<DirectSolver>& solvers_ext_);
          void set_solvers_ao(const DirectSolver solver_);
          void set_solvers_ao(const std::vector<DirectSolver>& solvers_);
          void set_solvers_ext(const DirectSolver solver_);
          void set_solvers_ext(const std::vector<DirectSolver>& solvers_);

          void set_paritioner(const std::vector<int>& subregions_, const PartitionType parti_type_, const PartitionApproach parti_approach_, const double parti_imbl_);
          void set_paritioner(const std::vector<int>& subregions_, const std::vector<PartitionType>& parti_types_, const std::vector<PartitionApproach>& parti_approach_, const std::vector<double>& parti_imbl_);

          void set_subregions(const std::vector<int>& subregions_);

          void set_parti_types(const PartitionType parti_type_);
          void set_parti_types(const std::vector<PartitionType>& parti_types_);

          void set_parti_approach(const PartitionApproach parti_approach_);
          void set_parti_approach(const std::vector<PartitionApproach>& parti_approach_);

          void set_parti_imbl(const double parti_imbl_);
          void set_parti_imbl(const std::vector<double>& parti_imbl_);

          void set_excludes(const bool velo_pres_, const bool pres_velo_);
          void set_excludes(const std::vector<bool>& velo_pres_, const std::vector<bool>& pres_velo_);
          void set_excludes_velo_pres(const bool velo_pres_);
          void set_excludes_velo_pres(const std::vector<bool>& velo_pres_);
          void set_excludes_pres_velo(const bool pres_velo_);
          void set_excludes_pres_velo(const std::vector<bool>& pres_velo_);

          void set_ipous(const IPOU ipou_velo_, const IPOU ipou_pres_);
          void set_ipous(const std::vector<IPOU>& ipous_velo_, const std::vector<IPOU>& ipous_pres_);
          void set_ipous_velo(const IPOU ipou_);
          void set_ipous_velo(const std::vector<IPOU>& ipous_);
          void set_ipous_pres(const IPOU ipou_);
          void set_ipous_pres(const std::vector<IPOU>& ipous_);

          void set_cspace(const bool cspace_velo_, const bool cspace_pres_);
          void set_cspace(const std::vector<bool>& cspace_velo_, const std::vector<bool>& cspace_pres_);
          void set_cspace_velo(const bool cspace_);
          void set_cspace_velo(const std::vector<bool>& cspace_);
          void set_cspace_pres(const bool cspace_);
          void set_cspace_pres(const std::vector<bool>& cspace_);

          void set_reuse_sf(const bool rsf_ao_, const bool rsf_cm_);
          void set_reuse_sf(const std::vector<bool> &rsf_ao_, const std::vector<bool> &rsf_cm_);
          void set_reuse_sf_ao(const bool rsf_ao_);
          void set_reuse_sf_ao(const std::vector<bool> &rsf_ao_);
          void set_reuse_sf_cm(const bool rsf_cm_);
          void set_reuse_sf_cm(const std::vector<bool> &rsf_cm_);

          void set_reuse_coarse(const bool reuse_cm_, const bool reuse_cb_);
          void set_reuse_coarse(const std::vector<bool> &reuse_cm_, const std::vector<bool> &reuse_cb_);
          void set_reuse_coarse_matrix(const bool reuse_cm_);
          void set_reuse_coarse_matrix(const std::vector<bool> &reuse_cm_);
          void set_reuse_coarse_basis(const bool reuse_cb_);
          void set_reuse_coarse_basis(const std::vector<bool> &reuse_cb_);

          void set_phi_dropping_threshold(const double phi_dt_);
          void set_phi_dropping_threshold(const std::vector<double> &phi_dt_);

          void set_verbose(const bool verbose_);

          bool parse_args(SimpleArgParser& args);

        protected:
          /// a pointer to our opaque core wrapper object
          /// Teuchos::ParameterList* pointer
          void* _core;

          /* for parsing of arguments */
          const Dist::Comm& _comm;

          /* spatial dimensios */
          int _dimension;
          /* dofs per entity (velocity) */
          int _dpe_velo;
          /* dofs per pressure (velocity) */
          int _dpe_pres;

          /* Number of levels (Domain Decomposition levels) before the coarse problem is solved
           *
           * For example:
           *     The two-level preconditioner has _nlevels = 1
           *     The three-level preconditioner has _nlevels = 2
           * */
          int _nlevels;

          ProblemType _problem;
          Preconditioner _preconditioner;
          DirectSolver _solver_coarse;

          std::vector<int> _gsteps;

          std::vector<int> _overlaps;
          std::vector<DirectSolver> _solvers_ao;
          std::vector<DirectSolver> _solvers_ext;

          /* __ATTENTION:__ the last entry has to be one */
          std::vector<int> _subregions;
          std::vector<PartitionType> _partition_types;
          std::vector<PartitionApproach> _partition_approach;
          std::vector<double> _partition_imbl_tol;

          /* exclude blocks for the coarse space:
           *
           *     _exclude_velo_pres:  exclude the velocity-pressure coupling
           *     _exclude_pres_velo:  exclude the pressure-velocity coupling
           */
          std::vector<bool> _exclude_pres_velo;
          std::vector<bool> _exclude_velo_pres;
          /*
           * use block for coarse space.  Normally, set this both to true.
           */
          std::vector<bool> _cspace_velo;
          std::vector<bool> _cspace_pres;
          /* interface partition of unity for blocks */
          std::vector<IPOU> _ipous_velo;
          std::vector<IPOU> _ipous_pres;

          /* combine mode of overlap */
          std::vector<CombineOverlap> _combine_ov;

          /* combine mode of lvl */
          std::vector<CombineLevel> _combine_lvl;

          /* reuse symbolic factorization: algebraic overlapping Operator */
          std::vector<bool> _ao_rsf;
          /* reuse symbolic factorization: coarse matrix */
          std::vector<bool> _cm_rsf;

          /* reuse the complete coarse matrix */
          std::vector<bool> _cm_r;
          /* reuse the complete coarse basis */
          std::vector<bool> _cb_r;
          /* Phi: Dropping Threshold ==> default 1.e-8.
           * The higher the dropping threshold the more Krylov iterations
           * but the faster the construction of the coarse basis */
          std::vector<double> _phi_dt;

          std::string _name;

          /* print the parameter list after construction */
          bool _print_list;
          /* use Trilinos timers for FROSch solver */
          bool _use_timer;

          /* verbose flag for the FROSch/Trilinos class */
          bool _verbose;

          void init_defaults();

      }; // FROSchParameterList

      inline void add_supported_fpl_args(SimpleArgParser& args)
      {
        args.support("frosch-plist", "<XML file>\n"
                     "The frosch parameter list as a Trilinos parameter list.\n"
                     "__ATTENTION:__ If this option is given, then all the other options are ignored.\n"
                     "May contain the following types: XML file path."
        );
        args.support("frosch-print-list", "<bool>\n"
                     "The Trilinos parameter list is print after the construction.\n"
                     "May contain the following types: true, false."
        );
        args.support("frosch-verbose", "<bool>\n"
                     "The Trilinos/FROSch verbose.  Print output of FROSch.\n"
                     "May contain the following types: true, false."
        );
        args.support("frosch-use-timer", "<bool>\n"
                     "Use the Trilinos timers for FROSch solver\n"
                     "May contain the following types: true, false."
        );
        args.support("frosch-nlevels", "<int>\n"
                     "Specifies the (additional) number of levels.\n"
                     "__ATTENTION:__ The used number of levels are one more, since the coarsest level is not counted.\n"
                     "               Therefore, two-level overlapping Schwarz has nlevels == 1.\n"
                     "May contain the following types: int > 0"
        );
        args.support("frosch-dim", "<int>\n"
                     "Specifies the (spatial) dimension.\n"
                     "May contain the following types: int > 0"
        );
        args.support("frosch-overlap", "<int | int...>\n"
                     "Specifies which number of overlaps are used.\n"
                     "Number of arguments has be one or same as the value of the argument: frosch-nlevels (int).\n"
                     "May contain the following types: int > 0\n"
                     "Default is 1."
        );
        args.support("frosch-subreg", "<int | int...>\n"
                     "Specifies number of subregions for multi-level approach.\n"
                     "Number of arguments has same as the value of the argument: frosch-nlevels (int).\n"
                     "The last entry has to be one.\n"
                     "May contain the following types: int > 0"
        );
        args.support("frosch-gstep", "<int | int...>\n"
                     "Specifies number of gathering steps for the coarse problem.\n"
                     "Number of arguments has same as the value of the argument: frosch-nlevels (int).\n"
                     "May contain the following types: int > 0"
        );
        args.support("frosch-parti-type", "<types | types...>\n"
                     "Specifies the paritioner for multi-level approach.\n"
                     "Number of arguments has be one or same as the value of the argument: frosch-nlevels (int).\n"
                     "May contain the following types: block, parmetis, phg, zoltan, scotch, ptscotch"
        );
        args.support("frosch-parti-approach", "<types | types...>\n"
                     "Specifies the paritioner approach for multi-level approach.\n"
                     "Only used when parmetis partitioning is used\n."
                     "Number of arguments has be one or same as the value of the argument: frosch-nlevels (int).\n"
                     "May contain the following types: partition, repartition."
        );
        args.support("frosch-parti-imbl", "<double | double...>\n"
                     "Specifies the paritioner imblance for multi-level approach.\n"
                     "Number of arguments has be one or same as the value of the argument: frosch-nlevels (int).\n"
                     "May contain the following types: double > 1."
        );
        args.support("frosch-ipou-pres", "<types | types...>\n"
                     "Specifies the interface partition of unity for the pressure space.\n"
                     "Number of arguments has be one or same as the value of the argument: frosch-nlevels (int).\n"
                     "May contain the following types: gdsw, gdswstar, rgdsw."
        );
        args.support("frosch-ipou-velo", "<types | types...>\n"
                     "Specifies the interface partition of unity for the velocity space.\n"
                     "Number of arguments has be one or same as the value of the argument: frosch-nlevels (int).\n"
                     "May contain the following types: gdsw, gdswstar, rgdsw."
        );
        args.support("frosch-solver-ao", "<types | types...>\n"
                     "Specifies the overlapping direct solver.\n"
                     "Number of arguments has be one or same as the value of the argument: frosch-nlevels (int).\n"
                     "May contain the following types: klu, mumps, umfpack."
        );
        args.support("frosch-solver-ext", "<types | types...>\n"
                     "Specifies the extension direct solver.\n"
                     "Number of arguments has be one or same as the value of the argument: frosch-nlevels (int).\n"
                     "May contain the following types: klu, mumps, umfpack."
        );
        args.support("frosch-solver-coarse", "<types>\n"
                     "Specifies the coarse direct solver.\n"
                     "May contain the following types: klu, mumps, umfpack."
        );
        args.support("frosch-precond-type", "<types>\n"
                     "Specifies FROSch preconditioner type.\n"
                     "May contain the following types: onelevel, twolevel, twolevelblock."
        );
        args.support("frosch-cspace-pres", "<bool | bool...>\n"
                     "Specifies the use of the pressure for the coarse space (monolithic only).\n"
                     "Number of arguments has be one or same as the value of the argument: frosch-nlevels (int).\n"
                     "May contain the following types: true, false."
        );
        args.support("frosch-cspace-velo", "<bool | bool...>\n"
                     "Specifies the use of the velocity for the coarse space (monolithic only).\n"
                     "Number of arguments has be one or same as the value of the argument: frosch-nlevels (int).\n"
                     "May contain the following types: true, false."
        );
        args.support("frosch-exclude-velo-pres", "<bool | bool...>\n"
                     "Exclude the velocity-pressure coupling in coarse space (monolithic only).\n"
                     "Number of arguments has be one or same as the value of the argument: frosch-nlevels (int).\n"
                     "May contain the following types: true, false."
        );
        args.support("frosch-exclude-pres-velo", "<bool | bool...>\n"
                     "Exclude the pressure-velocity coupling in coarse space (monolithic only).\n"
                     "Number of arguments has be one or same as the value of the argument: frosch-nlevels (int).\n"
                     "May contain the following types: true, false."
        );
        args.support("frosch-combine-overlap", "<types | types...>\n"
                     "Combine mode of the overlap.\n"
                     "Number of arguments has be one or same as the value of the argument: frosch-nlevels (int).\n"
                     "May contain the following types: full, restricted, averaging."
        );
        args.support("frosch-combine-lvl", "<types | types...>\n"
                     "Combine mode of the different levels.\n"
                     "Number of arguments has be one or same as the value of the argument: frosch-nlevels (int).\n"
                     "May contain the following types: additive, multiplicative"
        );
        args.support("frosch-reuse-sf-ao", "<bool | bool...>\n"
                     "Reuse the symbolic factorization of the algebraic overlapping operator.\n"
                     "Number of arguments has be one or same as the value of the argument: frosch-nlevels (int).\n"
                     "May contain the following types: true, false"
        );
        args.support("frosch-reuse-sf-cm", "<bool | bool...>\n"
                     "Reuse the symbolic factorization of the coarse matrix.\n"
                     "Number of arguments has be one or same as the value of the argument: frosch-nlevels (int).\n"
                     "May contain the following types: true, false"
        );
        args.support("frosch-reuse-coarse-matrix", "<bool | bool...>\n"
                     "Reuse the complete coarse matrix.\n"
                     "Number of arguments has be one or same as the value of the argument: frosch-nlevels (int).\n"
                     "May contain the following types: true, false"
        );
        args.support("frosch-reuse-coarse-basis", "<bool | bool...>\n"
                     "Reuse the coarse basis functions.\n"
                     "Number of arguments has be one or same as the value of the argument: frosch-nlevels (int).\n"
                     "May contain the following types: true, false"
        );
        args.support("frosch-phi-dropping-threshold", "<double | double...>\n"
                     "Dropping threshold for the coarse matrix.\n"
                     "Number of arguments has be one or same as the value of the argument: frosch-nlevels (int).\n"
                     "May contain the following types: double > 0."
        );
      }

      inline std::shared_ptr<FROSchParameterList> init_params_two_level_pressure_poisson(const Dist::Comm& comm_,
                                                                                         const int dim_,
                                                                                         const int overlap_,
                                                                                         const FROSchParameterList::IPOU ipou_pres_ = FROSchParameterList::GDSW,
                                                                                         const FROSchParameterList::DirectSolver solver_ao_ = FROSchParameterList::KLU,
                                                                                         const FROSchParameterList::DirectSolver solver_ext_ = FROSchParameterList::KLU,
                                                                                         const FROSchParameterList::DirectSolver solver_coarse_ = FROSchParameterList::KLU,
                                                                                         const FROSchParameterList::Preconditioner precond_ = FROSchParameterList::TWOLEVELBLOCK,
                                                                                         const FROSchParameterList::CombineOverlap mode_ = FROSchParameterList::RESTRICTED,
                                                                                         const bool use_timer_ = false,
                                                                                         const bool print_list_ = false)
      {
        std::shared_ptr<FROSchParameterList> params = std::make_shared<FROSchParameterList>(comm_, dim_, FROSchParameterList::PRESSUREPOISSON);
        params->set_overlaps(overlap_);
        /* IPOU for velocity will not be used, set it anyway */
        params->set_ipous(FROSchParameterList::GDSW, ipou_pres_);
        params->set_solvers(solver_ao_, solver_ext_);
        params->set_coarse_solver(solver_coarse_);
        params->set_precond(precond_);
        params->set_combine_overlap(mode_);
        params->set_use_timer(use_timer_);
        params->set_print(print_list_);

        params->set_cspace(true, true);
        /* Will not be used, set it anyway */
        params->set_excludes(false, false);
        params->set_paritioner(std::vector<int>(1, 1),
                               Solver::Trilinos::FROSchParameterList::PARMETIS,     /* DUMMY NOT USED, since  std::vector<int>(1, 1) */
                               Solver::Trilinos::FROSchParameterList::REPARTITION,  /* DUMMY NOT USED, since  std::vector<int>(1, 1) */
                               1.1                                            /* DUMMY NOT USED, since  std::vector<int>(1, 1) */);
        return params;
      }

      inline std::shared_ptr<FROSchParameterList> init_params_three_level_pressure_poisson(const Dist::Comm& comm_,
                                                                                           const int dim_,
                                                                                           const int overlap_,
                                                                                           const int nsubreg_,
                                                                                           const FROSchParameterList::PartitionType parti_type_ = FROSchParameterList::PARMETIS,
                                                                                           const FROSchParameterList::PartitionApproach parti_approach_ = FROSchParameterList::REPARTITION,
                                                                                           const double parti_imbl_ = 1.1,
                                                                                           const FROSchParameterList::IPOU ipou_pres_ = FROSchParameterList::GDSW,
                                                                                           const FROSchParameterList::DirectSolver solver_ao_ = FROSchParameterList::KLU,
                                                                                           const FROSchParameterList::DirectSolver solver_ext_ = FROSchParameterList::KLU,
                                                                                           const FROSchParameterList::DirectSolver solver_coarse_ = FROSchParameterList::KLU,
                                                                                           const FROSchParameterList::Preconditioner precond_ = FROSchParameterList::TWOLEVELBLOCK,
                                                                                           const FROSchParameterList::CombineOverlap mode_ = FROSchParameterList::RESTRICTED,
                                                                                           const bool use_timer_ = false,
                                                                                           const bool print_list_ = false)
      {
        std::shared_ptr<FROSchParameterList> params = std::make_shared<FROSchParameterList>(comm_, dim_, FROSchParameterList::PRESSUREPOISSON, 2);
        params->set_overlaps(overlap_);
        params->set_paritioner(std::vector<int>({nsubreg_, 1}), parti_type_, parti_approach_, parti_imbl_);
        /* IPOU for velocity will not be used, set it anyway */
        params->set_ipous(FROSchParameterList::GDSW, ipou_pres_);
        params->set_solvers(solver_ao_, solver_ext_);
        params->set_coarse_solver(solver_coarse_);
        params->set_precond(precond_);
        params->set_combine_overlap(mode_);
        params->set_use_timer(use_timer_);
        params->set_print(print_list_);

        params->set_cspace(true, true);
        /* Will not be used, set it anyway */
        params->set_excludes(false, false);
        return params;
      }

      inline std::shared_ptr<FROSchParameterList> init_params_two_level_monolithic(const Dist::Comm& comm_,
                                                                                   const int dim_,
                                                                                   const int overlap_,
                                                                                   const FROSchParameterList::IPOU ipou_velo_ = FROSchParameterList::GDSW,
                                                                                   const FROSchParameterList::IPOU ipou_pres_ = FROSchParameterList::GDSW,
                                                                                   const FROSchParameterList::DirectSolver solver_ao_ = FROSchParameterList::KLU,
                                                                                   const FROSchParameterList::DirectSolver solver_ext_ = FROSchParameterList::KLU,
                                                                                   const FROSchParameterList::DirectSolver solver_coarse_ = FROSchParameterList::KLU,
                                                                                   const bool excl_velo_pres_ = false,
                                                                                   const bool excl_pres_velo_ = false,
                                                                                   const FROSchParameterList::CombineOverlap mode_ = FROSchParameterList::RESTRICTED,
                                                                                   const bool use_timer_ = false,
                                                                                   const bool print_list_ = false)
      {
        std::shared_ptr<FROSchParameterList> params = std::make_shared<FROSchParameterList>(comm_, dim_, FROSchParameterList::PRESSUREPOISSON);
        params->set_overlaps(overlap_);
        params->set_ipous(ipou_velo_, ipou_pres_);
        params->set_solvers(solver_ao_, solver_ext_);
        params->set_coarse_solver(solver_coarse_);
        params->set_excludes(excl_velo_pres_, excl_pres_velo_);
        params->set_combine_overlap(mode_);
        params->set_use_timer(use_timer_);
        params->set_print(print_list_);

        params->set_precond(FROSchParameterList::TWOLEVELBLOCK);
        params->set_cspace(true, true);
        params->set_paritioner(std::vector<int>(1, 1),
                               Solver::Trilinos::FROSchParameterList::PARMETIS,     /* DUMMY NOT USED, since  std::vector<int>(1, 1) */
                               Solver::Trilinos::FROSchParameterList::REPARTITION,  /* DUMMY NOT USED, since  std::vector<int>(1, 1) */
                               1.1                                            /* DUMMY NOT USED, since  std::vector<int>(1, 1) */);
        return params;
      }

      inline std::shared_ptr<FROSchParameterList> init_params_three_level_monolithic(const Dist::Comm& comm_,
                                                                                     const int dim_,
                                                                                     const int overlap_,
                                                                                     const int nsubreg_,
                                                                                     const FROSchParameterList::PartitionType parti_type_ = FROSchParameterList::PARMETIS,
                                                                                     const FROSchParameterList::PartitionApproach parti_approach_ = FROSchParameterList::REPARTITION,
                                                                                     const double parti_imbl_ = 1.1,
                                                                                     const FROSchParameterList::IPOU ipou_velo_ = FROSchParameterList::GDSW,
                                                                                     const FROSchParameterList::IPOU ipou_pres_ = FROSchParameterList::GDSW,
                                                                                     const FROSchParameterList::DirectSolver solver_ao_ = FROSchParameterList::KLU,
                                                                                     const FROSchParameterList::DirectSolver solver_ext_ = FROSchParameterList::KLU,
                                                                                     const FROSchParameterList::DirectSolver solver_coarse_ = FROSchParameterList::KLU,
                                                                                     const bool excl_velo_pres_ = false,
                                                                                     const bool excl_pres_velo_ = false,
                                                                                     const FROSchParameterList::CombineOverlap mode_ = FROSchParameterList::RESTRICTED,
                                                                                     const bool use_timer_ = false,
                                                                                     const bool print_list_ = false)
      {
        std::shared_ptr<FROSchParameterList> params = std::make_shared<FROSchParameterList>(comm_, dim_, FROSchParameterList::PRESSUREPOISSON, 2);
        params->set_overlaps(overlap_);
        params->set_paritioner(std::vector<int>({nsubreg_, 1}), parti_type_, parti_approach_, parti_imbl_);
        params->set_ipous(ipou_velo_, ipou_pres_);
        params->set_solvers(solver_ao_, solver_ext_);
        params->set_coarse_solver(solver_coarse_);
        params->set_excludes(excl_velo_pres_, excl_pres_velo_);
        params->set_combine_overlap(mode_);
        params->set_use_timer(use_timer_);
        params->set_print(print_list_);

        params->set_precond(FROSchParameterList::TWOLEVELBLOCK);
        params->set_cspace(true, true);
        return params;
      }

      /**
       * \brief Creates a core wrapper object for Tpetra.
       *
       * The "core wrapper object", which is returned by this function, is basically a collection
       * of all Tpetra data structures, which are required to represent a linear system A*x=b, i.e.
       * a partitioned matrix as well as two vectors.
       *
       * \param[in] comm
       * A pointer to an MPI_Comm object that represents the communicator.
       * This argument is ignored in non-MPI builds.
       *
       * \param[in] num_global_dofs
       * The number of global DOFs across all processes.
       *
       * \param[in] my_dof_offset
       * The global DOF offset for this process.
       *
       * \param[in] num_owned_dofs
       * The number of global DOFs owned by this process.
       *
       * \param[in] num_nonzeros
       * The number of non-zero entries for this process.
       *
       * \param[in] params
       * The FROSchParameterList object for Trilinos.
       *
       * \returns
       * A pointer to a newly allocated core wrapper object.
       *
       * \author Stephan Köhler
       */
      void* create_core_scalar(const void* comm, Index num_global_dofs, Index my_dof_offset,
        Index num_owned_dofs, Index num_nonzeros, const FROSchParameterList & params);

      // Trilinos FROSch Stokes
      void* create_core_stokes(const void* comm, Index num_global_dofs, Index my_dof_offset,
        Index num_owned_dofs, Index num_nonzeros, Index num_owned_velo_dofs, Index num_owned_pres_dofs,
        Index first_owned_velo_dof, Index first_owned_pres_dof,
        Index num_global_velo_dofs, Index num_global_pres_dofs,
        const Trilinos::FROSchParameterList & params);

      void destroy_core(void* core);

      void set_parameter_list(void* core, const FROSchParameterList & params);

      void init_symbolic(void* core);
      void init_numeric(void* core);

      std::int64_t* get_row_ptr(void* core);
      std::int32_t* get_col_idx(void* core);
      double* get_mat_val(void* core);
      double* get_vec_def(void* core);
      double* get_vec_cor(void* core);
      void format_vec_cor(void* core);

      // FROSch Wrappers
      void* create_frosch(void* core);
      void reinit_frosch(void *core);
      void compute_frosch(void *solver);
      void destroy_frosch(void* solver);
      void solve_frosch(void* core, void* solver);

    } // namespace Trilinos
    /// \endcond

    /**
     * \brief Base-Class for solvers/preconditioners borrowed from Trilinos/FROSch library
     *
     * This class acts as a base-class for all solvers that we borrow from the FROSch library.
     * We use a CSR format, which corresponds to the "algebraic DOF partitioning", this class
     * derives from the ADPSolverBase class, which
     * takes care of the translation between the system matrix and the ADP data structures.
     *
     * This base-class takes care of allocating, initialising and updating the required
     * Tpetra matrix and vector objects, which are outsourced into an opaque core wrapper object.
     *
     * \author Stephan Köhler
     */
    template<typename Matrix_, typename Filter_, typename SolverBase_ = Solver::SolverBase<typename Matrix_::VectorTypeL>>
    class TpetraSolverBase :
      public ADPSolverBase<Matrix_, Filter_, SolverBase_>
    {
    public:
      /// our base-class
      typedef ADPSolverBase<Matrix_, Filter_, SolverBase_> BaseClass;
      // the vector type
      typedef typename BaseClass::VectorType VectorType;

    protected:
      /// a pointer to our opaque core wrapper object
      void* _core;
      const Trilinos::FROSchParameterList& _params;

      explicit TpetraSolverBase(const Matrix_& matrix, const Filter_& filter, const Trilinos::FROSchParameterList & params) :
        BaseClass(matrix, filter),
        _core(nullptr),
        _params(params)
      {
      }

      /**
       * \brief Uploads the Tpetra defect vector
       *
       * This function first uploads the given defect vector into its ADP defect vector
       * counterpart and afterwards uploads that into the Tpetra vector counterpart.
       *
       * \param[in] vec_def
       * The defect vector to be uploaded from.
       */
      void _upload_def(const VectorType& vec_def)
      {
        this->_upload_vector(Trilinos::get_vec_def(this->_core), vec_def.local());
      }

      /**
       * \brief Format the Tpetra correction vector
       */
      void _format_cor()
      {
        // format Tpetra correction vector
        Trilinos::format_vec_cor(this->_core);
      }

      /**
       * \brief Downloads the Tpetra correction vector
       *
       * This function first downloads the Tpetra vector into its ADP correction vector
       * counterpart and afterwards downloads that into the given correction vector.
       *
       * \param[out] vec_cor
       * The correction vector to download to.
       */
      void _download_cor(VectorType& vec_cor)
      {
        this->_download_vector(vec_cor.local(), Trilinos::get_vec_cor(this->_core));

        // apply correction filter
        this->_system_filter.filter_cor(vec_cor);
      }

    public:
      /**
       * \brief Numeric Initialization
       *
       * This function uploads the numerical values of the ADP matrix, which
       * is managed by the base-class, to the Tpetra matrix.
       */
      virtual void init_numeric() override
      {
        BaseClass::init_numeric();

        XASSERT(this->_core != nullptr);

        // update matrix values of Trilinos matrix
        this->_upload_numeric(Trilinos::get_mat_val(this->_core),
          Trilinos::get_row_ptr(this->_core), Trilinos::get_col_idx(this->_core));

        Trilinos::init_numeric(this->_core);
      }

      /**
       * \brief Symbolic Finalization
       *
       * This function destroys all Tpetra objects managed by this class
       * and resets all auxiliary vectors and pointers.
       */
      virtual void done_symbolic() override
      {
        XASSERT(this->_core != nullptr);

        Trilinos::destroy_core(this->_core);
        this->_core = nullptr;

        BaseClass::done_symbolic();
      }
    }; // class TpetraSolverBase<...>

    /**
     * \brief Trilinos FROSchPreconditioner class template
     *
     * This class acts as a wrapper around the FROSchPreconditioner from the Trilinos library.
     *
     * \todo support setting of solver parameters
     *
     * \author Stephan Köhler
     */
    template<typename Matrix_, typename Filter_>
    class FROSchPreconditioner :
      public TpetraSolverBase<Matrix_, Filter_>
    {
    public:
      /// our base-class
      typedef TpetraSolverBase<Matrix_, Filter_> BaseClass;
      /// the vector type
      typedef typename BaseClass::VectorType VectorType;

    protected:
      /// the FROSchPreconditioner solver object
      void* _solver;

    public:
      explicit FROSchPreconditioner(const Matrix_& matrix, const Filter_& filter, const Trilinos::FROSchParameterList & params) :
        BaseClass(matrix, filter, params),
        _solver(nullptr)
      {
      }

      explicit FROSchPreconditioner(const String& DOXY(section_name), PropertyMap* DOXY(section),
        const Matrix_& matrix, const Filter_& filter)
         :
        FROSchPreconditioner(matrix, filter)
      {
        /// \todo parse parameters
      }

      virtual String name() const override
      {
        return "FROSchPreconditioner";
      }

      /**
       * \brief Symbolic Initialization
       *
       * This function creates the Tpetra matrix and vector objects and initializes their
       * structure/layout by using the algebraic DOF partitioning (ADP) that is managed by the
       * base-class. This function also performs an initial upload of the matrix and vector
       * values from the ADP structures (because Tpetra requires this), although these values
       * may be undefined (but existent) at this point.
       */
      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();

        XASSERT(this->_core == nullptr);

        // create our Tpetra core wrapper object
        this->_core = Trilinos::create_core_scalar(
          &this->_get_comm()->mpi_comm(),
          this->_get_num_global_dofs(),
          this->_get_global_dof_offset(),
          this->_get_num_owned_dofs(),
          this->_get_adp_matrix_num_nzes(),
          this->_params);

        BaseClass::_upload_symbolic(Trilinos::get_row_ptr(this->_core), Trilinos::get_col_idx(this->_core));

        Trilinos::init_symbolic(this->_core);
      }

      virtual void init_numeric() override
      {
        BaseClass::init_numeric();

        // create FROSchPreconditioner
        if(_solver == nullptr)
        {
          this->_solver = Trilinos::create_frosch(this->_core);
        }
        else
        {
          Trilinos::reinit_frosch(this->_solver);
        }
      }

      virtual void done_numeric() override
      {
        // if(_solver != nullptr)
        // {
        //   Trilinos::destroy_frosch(_solver);
        //   _solver = nullptr;
        // }
      }

      virtual ~FROSchPreconditioner() override
      {
        if (_solver != nullptr)
        {
          Trilinos::destroy_frosch(_solver);
          _solver = nullptr;
        }
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        // upload defect vector and format correction
        this->_upload_def(vec_def);
        this->_format_cor();

        // apply FROSchPreconditioner
        Trilinos::solve_frosch(this->_core, this->_solver);

        // download correction
        this->_download_cor(vec_cor);

        // okay
        return Status::success;
      }
    }; // class FROSchPreconditioner

    /**
     * \brief Creates a new FROSchPreconditioner solver object
     *
     * \param[in] matrix
     * The global system matrix.
     *
     * \param[in] filter
     * The global system filter.
     *
     * \returns
     * A shared pointer to a new FROSchPreconditioner object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<FROSchPreconditioner<Matrix_, Filter_>> new_frosch(
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<FROSchPreconditioner<Matrix_, Filter_>>(matrix, filter);
    }

    /**
     * \brief Creates a new FROSchPreconditioner solver object based on a PropertyMap
     *
     * \param[in] section_name
     * The name of the config section, which it does not know by itself
     *
     * \param[in] section
     * A pointer to the PropertyMap section configuring this solver
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \returns
     * A shared pointer to a new FROSchPreconditioner object.
     **/
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<FROSchPreconditioner<Matrix_, Filter_>> new_frosch(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<FROSchPreconditioner<Matrix_, Filter_>>(section_name, section, matrix, filter);
    }

    /**
     * \brief Creates a new FROSchPreconditioner solver object based on a parameter list
     *
     * \param[in] params
     * A pointer to the parameter list section configuring this solver
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \returns
     * A shared pointer to a new FROSchPreconditioner object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<FROSchPreconditioner<Matrix_, Filter_>> new_frosch(
      const Matrix_& matrix, const Filter_& filter, const Trilinos::FROSchParameterList & params)
    {
      return std::make_shared<FROSchPreconditioner<Matrix_, Filter_>>(matrix, filter, params);
    }

    /**
     * \brief Trilinos StokesFROSchPreconditioner class template
     *
     * This class acts as a wrapper around the FROSch Preconditioner (block variant) from the Trilinos library.
     *
     * \todo support setting of solver parameters
     *
     * \author Stephan Köhler
     */
    template<typename Matrix_, typename Filter_>
    class StokesFROSchPreconditioner :
      public TpetraSolverBase<Matrix_, Filter_>
    {
    public:
      /// our base-class
      typedef TpetraSolverBase<Matrix_, Filter_> BaseClass;
      /// the vector type
      typedef typename BaseClass::VectorType VectorType;

    protected:
      /// the FROSchPreconditioner solver object
      void* _solver;
      // number of velocity/pressure dofs owned by this process
      Index _num_owned_velo_dofs;
      Index _num_owned_pres_dofs;
      // index of first velocity/pressure dof owned by this process
      Index _first_owned_velo_dof;
      Index _first_owned_pres_dof;
      // total number of velocity/pressure dofs
      Index _num_global_velo_dofs;
      Index _num_global_pres_dofs;

    public:
      explicit StokesFROSchPreconditioner(const Matrix_& matrix, const Filter_& filter, const Trilinos::FROSchParameterList & params) :
        BaseClass(matrix, filter, params),
        _solver(nullptr),
        _num_owned_velo_dofs(0),
        _num_owned_pres_dofs(0),
        _first_owned_velo_dof(0),
        _first_owned_pres_dof(0),
        _num_global_velo_dofs(0),
        _num_global_pres_dofs(0)
      {
      }

      explicit StokesFROSchPreconditioner(const String& DOXY(section_name), PropertyMap* DOXY(section),
                                          const Matrix_& matrix, const Filter_& filter) :
        StokesFROSchPreconditioner(matrix, filter)
      {
        /// \todo parse parameters
      }

      virtual ~StokesFROSchPreconditioner() override
      {
        if (_solver != nullptr)
        {
          Trilinos::destroy_frosch(_solver);
          _solver = nullptr;
        }
      }

      virtual String name() const override
      {
        return "StokesFROSchPreconditioner";
      }

      /**
       * \brief Symbolic Initialization
       *
       * This function creates the Tpetra matrix and vector objects and initializes their
       * structure/layout by using the algebraic DOF partitioning (ADP) that is managed by the
       * base-class. This function also performs an initial upload of the matrix and vector
       * values from the ADP structures (because Tpetra requires this), although these values
       * may be undefined (but existent) at this point.
       *
       * Furthermore, the maps for FROSch (block variant) are built and written into the
       * parameterlist.
       */
      virtual void init_symbolic() override
      {
        // This is the init_symbolic of FROSchPreconditioner,
        // but the FROSchPreconditioner class does not override
        // the method of the TpetraSolverClass

        BaseClass::init_symbolic();

        this->_deduct_sizes();

        XASSERT(this->_num_global_velo_dofs + this->_num_global_pres_dofs == this->_get_num_global_dofs());
        XASSERT(this->_first_owned_velo_dof + this->_first_owned_pres_dof == this->_get_global_dof_offset());
        XASSERT(this->_core == nullptr);

        /*this->_get_comm()->allprint(
          "this->_get_num_global_dofs() = "     + stringify(this->_get_num_global_dofs()) + "\n" +
          "this->_get_global_dof_offset() = "   + stringify(this->_get_global_dof_offset()) + "\n" +
          "this->_get_num_owned_dofs() = "      + stringify(this->_get_num_owned_dofs()) + "\n" +
          "this->_get_adp_matrix_num_nzes() = " + stringify(this->_get_adp_matrix_num_nzes()) + "\n" +
          "this->_num_owned_velo_dofs = "       + stringify(this->_num_owned_velo_dofs) + "\n" +
          "this->_num_owned_pres_dofs = "       + stringify(this->_num_owned_pres_dofs) + "\n" +
          "this->_first_owned_velo_dof = "      + stringify(this->_first_owned_velo_dof) + "\n" +
          "this->_first_owned_pres_dof = "      + stringify(this->_first_owned_pres_dof) + "\n" +
          "this->_num_global_velo_dofs = "      + stringify(this->_num_global_velo_dofs) + "\n" +
          "this->_num_global_pres_dofs = "      + stringify(this->_num_global_pres_dofs));*/

        // create our Tpetra Stokes core wrapper object
        this->_core = Trilinos::create_core_stokes(
          &this->_get_comm()->mpi_comm(),
          this->_get_num_global_dofs(),
          this->_get_global_dof_offset(),
          this->_get_num_owned_dofs(),
          this->_get_adp_matrix_num_nzes(),
          this->_num_owned_velo_dofs,
          this->_num_owned_pres_dofs,
          this->_first_owned_velo_dof,
          this->_first_owned_pres_dof,
          this->_num_global_velo_dofs,
          this->_num_global_pres_dofs,
          this->_params);

        BaseClass::_upload_symbolic(Trilinos::get_row_ptr(this->_core), Trilinos::get_col_idx(this->_core));

        Trilinos::init_symbolic(this->_core);
       }

      virtual void init_numeric() override
      {
        BaseClass::init_numeric();

        // create FROSchPreconditioner
        if(_solver == nullptr)
        {
          this->_solver = Trilinos::create_frosch(this->_core);
        }
        else
        {
          Trilinos::reinit_frosch(this->_solver);
        }
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        // upload defect vector and format correction
        this->_upload_def(vec_def);
        this->_format_cor();

        // apply FROSchPreconditioner
        Trilinos::solve_frosch(this->_core, this->_solver);

        // download correction
        this->_download_cor(vec_cor);

        // okay
        return Status::success;
      }

    private:
      void _deduct_sizes()
      {
        // get and split the block information
        auto vbi = this->_get_adp_block_information().split_by_charset("\n");

        // block information must consist of 4 entries
        XASSERTM(vbi.size() == std::size_t(4), "invalid block information for StokesFROSchPreconditioner");

        // last block must terminate the tuple block
        XASSERTM(vbi[3] == "</Tuple>", "invalid block information for StokesFROSchPreconditioner");

        // try to parse the block information for the entire tuple
        std::regex rext("<Tuple gc=\"(\\d+)\" gf=\"(\\d+)\" go=\"(\\d+)\" lc=\"(\\d+)\" lo=\"0\">");
        std::smatch rmt;
        if(!std::regex_match(vbi.at(0), rmt, rext))
          throw InternalError("invalid Tuple block information for StokesFROSchPreconditioner:\n" + vbi.at(0));

        // try to parse the block information for the velocity component
        std::regex rexv("<Blocked bs=\"(\\d+)\" gc=\"(\\d+)\" gf=\"(\\d+)\" go=\"(\\d+)\" lc=\"(\\d+)\" lo=\"0\"/>");
        std::smatch rmv;
        if(!std::regex_match(vbi.at(1), rmv, rexv))
          throw InternalError("invalid velocity block information for StokesFROSchPreconditioner:\n" + vbi.at(1));

        if(!String(rmv[2]).parse(this->_num_global_velo_dofs))
          throw InternalError("invalid velocity block information for StokesFROSchPreconditioner: failed to parse global velocity dof count");
        if(!String(rmv[3]).parse(this->_first_owned_velo_dof))
          throw InternalError("invalid velocity block information for StokesFROSchPreconditioner: failed to parse first owned global velocity dof");
        if(!String(rmv[5]).parse(this->_num_owned_velo_dofs))
          throw InternalError("invalid velocity block information for StokesFROSchPreconditioner: failed to parse local velocity dof count");

        // try to parse the block information for the pressure component
        std::regex rexp("<Scalar gc=\"(\\d+)\" gf=\"(\\d+)\" go=\"(\\d+)\" lc=\"(\\d+)\" lo=\"(\\d+)\"/>");
        std::smatch rmp;
        if(!std::regex_match(vbi.at(2), rmp, rexp))
          throw InternalError("invalid velocity block information for StokesFROSchPreconditioner:\n" + vbi.at(2));

        if(!String(rmp[1]).parse(this->_num_global_pres_dofs))
          throw InternalError("invalid velocity block information for StokesFROSchPreconditioner: failed to parse global velocity dof count");
        if(!String(rmp[2]).parse(this->_first_owned_pres_dof))
          throw InternalError("invalid velocity block information for StokesFROSchPreconditioner: failed to parse first owned global velocity dof");
        if(!String(rmp[4]).parse(this->_num_owned_pres_dofs))
          throw InternalError("invalid velocity block information for StokesFROSchPreconditioner: failed to parse local velocity dof count");
      }

      /*void _deduct_sizes_old()
      {
        this->_num_owned_velo_dofs = this->_get_num_owned_dofs() - this->_num_owned_pres_dofs;

        this->_num_global_velo_dofs = this->_num_owned_velo_dofs;
        this->_num_global_pres_dofs = this->_num_owned_pres_dofs;
        this->_get_comm()->allreduce(&this->_num_global_velo_dofs, &this->_num_global_velo_dofs, 1u, Dist::op_sum);
        this->_get_comm()->allreduce(&this->_num_global_pres_dofs, &this->_num_global_pres_dofs, 1u, Dist::op_sum);

        this->_first_owned_velo_dof = 0;
        this->_first_owned_pres_dof = 0;
        this->_get_comm()->exscan(&this->_num_owned_velo_dofs, &this->_first_owned_velo_dof, 1u, Dist::op_sum);
        this->_get_comm()->exscan(&this->_num_owned_pres_dofs, &this->_first_owned_pres_dof, 1u, Dist::op_sum);
      }*/
    }; // class StokesFROSchPreconditioner

    /**
     * \brief Creates a new StokesFROSchPreconditioner solver object
     *
     * \param[in] matrix
     * The global system matrix.
     *
     * \param[in] filter
     * The global system filter.
     *
     * \param[in] num_owned_pres_dofs
     * The number of owned pressure dofs.
     *
     * \returns
     * A shared pointer to a new StokesFROSchPreconditioner object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<FROSchPreconditioner<Matrix_, Filter_>> new_stokes_frosch(
      const Matrix_& matrix, const Filter_& filter, const Index num_owned_pres_dofs)
    {
      return std::make_shared<StokesFROSchPreconditioner<Matrix_, Filter_>>(matrix, filter, num_owned_pres_dofs);
    }

    /**
     * \brief Creates a new StokesFROSchPreconditioner solver object based on a PropertyMap
     *
     * \param[in] section_name
     * The name of the config section, which it does not know by itself
     *
     * \param[in] section
     * A pointer to the PropertyMap section configuring this solver
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \returns
     * A shared pointer to a new FROSchPreconditioner object.
     **/
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<StokesFROSchPreconditioner<Matrix_, Filter_>> new_stokes_frosch(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<StokesFROSchPreconditioner<Matrix_, Filter_>>(section_name, section, matrix, filter);
    }

    /**
     * \brief Creates a new StokesFROSchPreconditioner solver object based on a parameter list
     *
     * \param[in] params
     * A pointer to the parameter list section configuring this solver
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \returns
     * A shared pointer to a new StokesFROSchPreconditioner object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<StokesFROSchPreconditioner<Matrix_, Filter_>> new_stokes_frosch(
      const Matrix_& matrix, const Filter_& filter, const Trilinos::FROSchParameterList & params)
    {
      return std::make_shared<StokesFROSchPreconditioner<Matrix_, Filter_>>(matrix, filter, params);
    }
  } // namespace Solver
} // namespace FEAT

#endif // defined(FEAT_HAVE_TRILINOS) || defined(DOXYGEN)
