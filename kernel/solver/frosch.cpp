#include <kernel/base_header.hpp>

#if defined(FEAT_HAVE_TRILINOS)

// preprocessor macro sanity checks
#ifdef FROSCH_SEQUENTIAL
#   error FROSCH_SEQUENTIAL does not make sense
#endif
#ifdef FEAT_HAVE_MPI
#  ifndef FROSCH_HAVE_MPI
#    define FROSCH_HAVE_MPI 1
#  endif
#else
#  error FROSch without MPI does not make sense.
#endif // FEAT_HAVE_MPI

#define FROSCH_TIMER_DETAILS 0

#include <ShyLU_DDFROSch_config.h>

#include <mpi.h>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#include <Teuchos_StackedTimer.hpp>

// Thyra includes
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>
#include <Thyra_SolveSupportTypes.hpp>
#include <Thyra_TpetraLinearOp.hpp>
#include <Thyra_TpetraMultiVector.hpp>
#include <Thyra_TpetraThyraWrappers.hpp>
#include <Thyra_TpetraVector.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_VectorStdOps.hpp>
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
#include <Thyra_EpetraLinearOp.hpp>
#endif
#include <Thyra_VectorSpaceBase_decl.hpp>
#include <Thyra_VectorSpaceBase_def.hpp>

// Stratimikos includes
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Stratimikos_FROSch_def.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

// FROSCH thyra includes
#include <Thyra_FROSchLinearOp_def.hpp>
#include <Thyra_FROSchFactory_def.hpp>
#include <FROSch_Tools_def.hpp>

// undefine assert to define it in feat, which was previously defined in trilinos
#undef ASSERT

// includes, FEAT
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/dist_file_io.hpp>

// includes, system
#include <vector>

using DataType = double;
using IndexType = FEAT::Index;

namespace FEAT
{
  namespace Solver
  {
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
          BLOCK      /* Should work always (not good) */
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
          ADDITIVE,
          MULTIPLICATIVE
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

          void set_reuse_sf(const bool rsf_ao_, const bool rsf_cm_, const bool rsf_ext_);
          void set_reuse_sf(const std::vector<bool> &rsf_ao_, const std::vector<bool> &rsf_cm_, const std::vector<bool> &rsf_ext_);
          void set_reuse_sf_ao(const bool rsf_ao_);
          void set_reuse_sf_ao(const std::vector<bool> &rsf_ao_);
          void set_reuse_sf_cm(const bool rsf_cm_);
          void set_reuse_sf_cm(const std::vector<bool> &rsf_cm_);
          void set_reuse_sf_ext(const bool rsf_ext_);
          void set_reuse_sf_ext(const std::vector<bool> &rsf_ext_);

          void set_reuse_coarse(const bool reuse_cm_, const bool reuse_cb_);
          void set_reuse_coarse(const std::vector<bool> &reuse_cm_, const std::vector<bool> &reuse_cb_);
          void set_reuse_coarse_matrix(const bool reuse_cm_);
          void set_reuse_coarse_matrix(const std::vector<bool> &reuse_cm_);
          void set_reuse_coarse_basis(const bool reuse_cb_);
          void set_reuse_coarse_basis(const std::vector<bool> &reuse_cb_);

          void set_phi_dropping_threshold(const double phi_dt_);
          void set_phi_dropping_threshold(const std::vector<double> &phi_dt_);

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
          /* reuse symbolic factorization: extension solver */
          std::vector<bool> _ext_rsf;

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

          void init_defaults();

      }; // FROSchParameterList

      inline void parameterlist_set_level_id(const int levelid, Teuchos::ParameterList &params)
      {
        params.set("Level ID", (unsigned int)levelid);
      }

      inline void parameterlist_set_default_block_id(const int defblockid, Teuchos::ParameterList &params)
      {
        params.set("Default BlockId Subdomain Graph", defblockid);
      }

      inline void parameterlist_set_dim(const int dimension, Teuchos::ParameterList &params)
      {
        params.set("Dimension", dimension);
      }

      inline void parameterlist_set_overlap(const int overlap, Teuchos::ParameterList &params)
      {
        params.set("Overlap", overlap);
      }

      inline void parameterlist_set_precond(const FROSchParameterList::Preconditioner precond, const bool reuse, Teuchos::ParameterList &params)
      {
        if(precond == FROSchParameterList::TWOLEVEL)
        {
          params.set("FROSch Preconditioner Type", "TwoLevelPreconditioner");
        } else if(precond == FROSchParameterList::TWOLEVELBLOCK)
        {
          params.set("FROSch Preconditioner Type", "TwoLevelBlockPreconditioner");
        } else if(precond == FROSchParameterList::ONELEVEL)
        {
          params.set("FROSch Preconditioner Type", "OneLevelPreconditioner");
        } else
        {
          XASSERTM(false, "Unkown FROSchParameterList::Preconditioner type");
        }
        params.set("Recycling", reuse);
      }

      inline void parameterlist_set_combine_lvl(const FROSchParameterList::CombineLevel combine_lvl, Teuchos::ParameterList &params)
      {
        if(combine_lvl == FROSchParameterList::ADDITIVE)
          params.set("Level Combination", "Additive");
        else if(combine_lvl == FROSchParameterList::MULTIPLICATIVE)
          params.set("Level Combination", "Multiplicative");
        else
          XASSERTM(false, "Unkown FROSchParameterList::CombineLevel type");
      }

      inline void parameterlist_header(const int dim, const int overlap, const int levelid, const int defblockid, const FROSchParameterList::Preconditioner precond, const std::string& nullspace, const FROSchParameterList::CombineLevel combine_lvl, const bool reuse, Teuchos::ParameterList &params)
      {
        parameterlist_set_dim(dim, params);
        parameterlist_set_overlap(overlap, params);
        parameterlist_set_default_block_id(defblockid, params);
        parameterlist_set_level_id(levelid, params);
        parameterlist_set_precond(precond, reuse, params);
        parameterlist_set_combine_lvl(combine_lvl, params);
        params.set("Null Space Type", nullspace);
        params.set("OverlappingOperator Type", "AlgebraicOverlappingOperator");
        params.set("CoarseOperator Type", "IPOUHarmonicCoarseOperator");
      }

      inline void parameterlist_set_solver(const FROSchParameterList::DirectSolver solver, Teuchos::ParameterList &params)
      {
        if(solver == FROSchParameterList::ILU)
        {
          params.set("SolverType", "Ifpack2");
          params.set("Solver", "RILUK");
          params.set("Ifpack2", Teuchos::ParameterList("Ifpack2"));
          params.sublist("Ifpack2").set("fact: iluk level-of-fill", 0);
        }else
        {
          params.set("SolverType", "Amesos2");
          if(solver == FROSchParameterList::KLU)
            params.set("Solver", "Klu");
          else if(solver == FROSchParameterList::MUMPS)
            params.set("Solver", "Mumps");
          else if(solver == FROSchParameterList::UMFPACK)
            params.set("Solver", "Umfpack");
          else
            XASSERTM(false, "Unkown FROSchParameterList::DirectSolver type");
        }
      }

      inline void parameterlist_set_aoo(const FROSchParameterList::DirectSolver solver, Teuchos::ParameterList &params)
      {
        params.set("AlgebraicOverlappingOperator", Teuchos::ParameterList("AlgebraicOverlappingOperator"));
        if(solver == FROSchParameterList::MUMPS)
          params.sublist("AlgebraicOverlappingOperator").set("Reuse: Symbolic Factorization", false);
        params.sublist("AlgebraicOverlappingOperator").set("Solver", Teuchos::ParameterList("Solver"));
        parameterlist_set_solver(solver, params.sublist("AlgebraicOverlappingOperator").sublist("Solver"));
      }

      inline void parameterlist_set_extsolv(const FROSchParameterList::DirectSolver solver, Teuchos::ParameterList &params)
      {
        params.set("ExtensionSolver", Teuchos::ParameterList("ExtensionSolver"));
        parameterlist_set_solver(solver, params.sublist("ExtensionSolver"));
      }

      inline void parameterlist_set_coarsesolv(const FROSchParameterList::DirectSolver solver, Teuchos::ParameterList &params)
      {
        params.set("CoarseSolver", Teuchos::ParameterList("CoarseSolver"));
        parameterlist_set_solver(solver, params.sublist("CoarseSolver"));
      }

      inline void parameterlist_combine_overlap(const FROSchParameterList::CombineOverlap mode, Teuchos::ParameterList &params)
      {
        if(mode == FROSchParameterList::FULL)
          params.set("Combine Values in Overlap", "Full");
        else if(mode == FROSchParameterList::RESTRICTED)
          params.set("Combine Values in Overlap", "Restricted");
        else if(mode == FROSchParameterList::AVERAGING)
          params.set("Combine Values in Overlap", "Averaging");
        else
          XASSERTM(false, "Unkown FROSchParameterList::DirectSolver type");
      }

      inline void parameterlist_set_zoltan2(const FROSchParameterList::PartitionType parti_type, const FROSchParameterList::PartitionApproach parti_app, const double parti_imbl, Teuchos::ParameterList &params)
      {
        params.set("Zoltan2 Parameter", Teuchos::ParameterList("Zoltan2 Parameter"));
        if(parti_type == FROSchParameterList::PARMETIS)
          params.sublist("Zoltan2 Parameter").set("algorithm", "parmetis");
        else if(parti_type == FROSchParameterList::BLOCK)
          params.sublist("Zoltan2 Parameter").set("algorithm", "block");
        else if(parti_type == FROSchParameterList::PHG)
          params.sublist("Zoltan2 Parameter").set("algorithm", "phg");
        else
          XASSERTM(false, "Unkown FROSchParameterList::PartitionType type");

        /* __TODO:__ Funktioniert das auch für "block" oder nur für "parmetis" */
        params.sublist("Zoltan2 Parameter").set("imbalance_tolerance", parti_imbl);
        if(parti_app == FROSchParameterList::PARTITION)
          params.sublist("Zoltan2 Parameter").set("partitioning_approach", "partition");
        else if(parti_app == FROSchParameterList::REPARTITION)
          params.sublist("Zoltan2 Parameter").set("partitioning_approach", "repartition");
        else
          XASSERTM(false, "Unkown FROSchParameterList::PartitionApproach type");
      }

      inline void parameterlist_set_distribution(const int subreg, const FROSchParameterList::PartitionType parti_type, const FROSchParameterList::PartitionApproach parti_app, const double parti_imbl, const int gsteps, Teuchos::ParameterList &params)
      {
        params.set("Distribution", Teuchos::ParameterList("Distribution"));
        params.sublist("Distribution").set("NumProcs", subreg);
        if(subreg == (int)1)
        {
          params.sublist("Distribution").set("Type", "linear");
        }else
        {
          params.sublist("Distribution").set("Type", "ZoltanDual");
          parameterlist_set_zoltan2(parti_type, parti_app, parti_imbl, params.sublist("Distribution"));
        }
        params.sublist("Distribution").set("Factor", (double)1.0);
        params.sublist("Distribution").set("GatheringSteps", gsteps);
        params.sublist("Distribution").set("Gathering Communication", Teuchos::ParameterList("Gathering Communication"));
        params.sublist("Distribution").sublist("Gathering Communication").set("Send type", "Send");
      }

      inline void parameterlist_set_ipou(const FROSchParameterList::IPOU ipou, Teuchos::ParameterList &params)
      {
        std::string ipou_str;
        if(ipou == FROSchParameterList::GDSW)
          ipou_str = "GDSW";
        else if(ipou == FROSchParameterList::GDSWSTAR)
          ipou_str = "GDSWStar";
        else if(ipou == FROSchParameterList::RGDSW)
          ipou_str = "RGDSW";
        else
          XASSERTM(false, "Unkown FROSchParameterList::IPOU type");
        params.set("InterfacePartitionOfUnity", Teuchos::ParameterList("InterfacePartitionOfUnity"));
        params.sublist("InterfacePartitionOfUnity").set("Type", ipou_str);
        params.sublist("InterfacePartitionOfUnity").set(ipou_str, Teuchos::ParameterList(ipou_str));
        params.sublist("InterfacePartitionOfUnity").sublist(ipou_str).set("Type", "Full");
        params.sublist("InterfacePartitionOfUnity").sublist(ipou_str).set("Distance Function", "Constant");
        params.sublist("InterfacePartitionOfUnity").sublist(ipou_str).set("Custom", Teuchos::ParameterList("Custom"));
        params.sublist("InterfacePartitionOfUnity").sublist(ipou_str).sublist("Custom").set("Roots", true);
      }

      inline void parameterlist_set_block(const int blockid, const bool use_cspace, const int exclude, const FROSchParameterList::IPOU ipou, Teuchos::ParameterList &params)
      {
        std::string blockid_str = std::to_string(blockid);
        std::string excl_str = std::to_string(exclude);

        params.set(blockid_str, Teuchos::ParameterList(blockid_str));
        params.sublist(blockid_str).set("Use For Coarse Space", use_cspace);
        if(exclude > 0)
          params.sublist(blockid_str).set("Exclude", excl_str);
        params.sublist(blockid_str).set("Rotations", false);
        parameterlist_set_ipou(ipou, params.sublist(blockid_str));
      }

      inline void parameterlist_set_reuse_coarse(const bool rsf_cm, const bool rsf_ext, const bool reuse_cm, const bool reuse_cb, Teuchos::ParameterList &params)
      {
        params.set("Reuse: Coarse Matrix Symbolic Factorization", rsf_cm);
        params.set("Reuse: Extension Symbolic Factorization", rsf_ext);
        params.set("Reuse: Coarse Matrix", reuse_cm);
        params.set("Reuse: Coarse Basis", reuse_cb);
      }

      FROSchParameterList::FROSchParameterList(const Dist::Comm& comm_, const int dim_, const FROSchParameterList::ProblemType problem_, const int nlevels_, const std::string& name_) :
        _core(nullptr),
        _comm(comm_),
        _dimension(dim_),
        _dpe_velo(dim_),
        _dpe_pres(dim_+1),
        _nlevels(nlevels_),
        _problem(problem_),
        _name(name_),
        _print_list(false),
        _use_timer(false)
      {
        XASSERTM(dim_ == 2 || dim_ == 3, "Only implemented for dim == 2 or dim == 3");
        init_defaults();
      }

      FROSchParameterList::FROSchParameterList(const Dist::Comm& comm_, const int dim_, const std::string& xml_file_, const std::string& name_) :
        _core(nullptr),
        _comm(comm_),
        _dimension(dim_),
        _dpe_velo(dim_),
        _dpe_pres(dim_+1),
        _nlevels(1),
        _name(name_),
        _print_list(false),
        _use_timer(false)
      {
        read_from_xml_file(xml_file_);
      }

      FROSchParameterList::~FROSchParameterList()
      {
        delete_core();
      }

      void FROSchParameterList::read_from_xml_str(const std::string& xml_str_)
      {
        delete_core();

        Teuchos::ParameterList *param = new Teuchos::ParameterList(_name);
        Teuchos::Ptr<Teuchos::ParameterList> ptr_param(param);
        Teuchos::updateParametersFromXmlString(xml_str_, ptr_param);

        _core = (void*)ptr_param.get();
      }

      void FROSchParameterList::read_from_xml_file(const std::string& xml_file_)
      {
        std::stringstream xml_str;
        DistFileIO::read_common(xml_str, xml_file_, _comm);
        read_from_xml_str(xml_str.str());
      }

      void FROSchParameterList::init_defaults()
      {
        set_overlaps(1);
        set_gsteps(1);

        /* Es fehlen noch die Partitionsparameter für nlevel > 1 */
        if(_nlevels == 1)
          set_subregions(std::vector<int>(1, 1));

        set_parti_types(FROSchParameterList::PARMETIS);
        set_parti_approach(FROSchParameterList::REPARTITION);
        set_parti_imbl(1.1);
        if(_dimension == 2)
          set_ipous(GDSW, GDSW);
        else
          set_ipous(GDSWSTAR, GDSW);

        set_solvers(FROSchParameterList::KLU, FROSchParameterList::KLU);
        _solver_coarse = FROSchParameterList::KLU;
        _preconditioner = FROSchParameterList::TWOLEVELBLOCK;
        set_cspace(true, true);
        set_excludes(false, false);
        /* __ATTENTION:__ Is produces a non-symmetric preconditioner */
        set_combine_overlap(FROSchParameterList::RESTRICTED);
        set_combine_lvl(FROSchParameterList::ADDITIVE);

        set_reuse_sf(false, false, false);
        set_reuse_coarse(false, false);
        set_phi_dropping_threshold(1.e-8);
      }

      void FROSchParameterList::print() const
      {
        if(_comm.rank() == 0)
        {
          XASSERTM(_core != nullptr, "_core == nullptr");
          (reinterpret_cast<Teuchos::ParameterList*>(_core))->print(std::cout, 4);
        }
      }

      void* FROSchParameterList::get_core() const
      {
        return _core;
      }

      void FROSchParameterList::create_core()
      {
        if(_core != nullptr)
        {
          /* core is already created  */
          return;
        }

        XASSERTM(_nlevels == int(_subregions.size()), "number of subregions is not consistent with number of levels.");
        XASSERTM(_subregions.back() == 1, "last entry of subregions is not 1.");

        std::string precond_str;
        if(_preconditioner == FROSchParameterList::TWOLEVEL)
          precond_str = "TwoLevelPreconditioner";
        else if(_preconditioner == FROSchParameterList::TWOLEVELBLOCK)
          precond_str = "TwoLevelBlockPreconditioner";
        else if(_preconditioner == FROSchParameterList::ONELEVEL)
          precond_str = "OneLevelPreconditioner";
        else
          XASSERTM(false, "Unkown FROSchParameterList::Preconditioner type");

        Teuchos::ParameterList *params = new Teuchos::ParameterList(_name);
        params->set("Preconditioner Type", "FROSch");
        params->set("Preconditioner Types", Teuchos::ParameterList("Preconditioner Types"));
        params->sublist("Preconditioner Types").set("FROSch", Teuchos::ParameterList("FROSch"));

        int overlap, gsteps;
        DirectSolver solver_ao, solver_ext;
        CombineOverlap mode;
        IPOU ipou_velo, ipou_pres;
        bool cspace_velo, cspace_pres, excl_velo_pres, excl_pres_velo;
        bool rsf_ao, rsf_ext, rsf_cm, reuse_cm, reuse_cb;
        double phi_dt;
        PartitionType partitype;
        PartitionApproach partiapp;
        double partiimbl;
        CombineLevel combine_lvl;
        bool reuse;

        Teuchos::ParameterList *tmp_params = &(params->sublist("Preconditioner Types").sublist("FROSch"));
        for(int i = 0; i < _nlevels; ++i)
        {
          overlap = _overlaps.size() == 1 ? _overlaps.at(0) : _overlaps.at(i);
          solver_ao = _solvers_ao.size() == 1 ? _solvers_ao.at(0) : _solvers_ao.at(i);
          solver_ext = _solvers_ext.size() == 1 ? _solvers_ext.at(0) : _solvers_ext.at(i);
          mode = _combine_ov.size() == 1 ? _combine_ov.at(0) : _combine_ov.at(i);
          combine_lvl = _combine_ov.size() == 1 ? _combine_lvl.at(0) : _combine_lvl.at(i);
          ipou_velo = _ipous_velo.size() == 1 ? _ipous_velo.at(0) : _ipous_velo.at(i);
          ipou_pres = _ipous_pres.size() == 1 ? _ipous_pres.at(0) : _ipous_pres.at(i);
          cspace_velo = _cspace_velo.size() == 1 ? _cspace_velo.at(0) : _cspace_velo.at(i);
          cspace_pres = _cspace_pres.size() == 1 ? _cspace_pres.at(0) : _cspace_pres.at(i);
          excl_velo_pres = _exclude_velo_pres.size() == 1 ? _exclude_velo_pres.at(0) : _exclude_velo_pres.at(i);
          excl_pres_velo = _exclude_pres_velo.size() == 1 ? _exclude_pres_velo.at(0) : _exclude_pres_velo.at(i);
          partitype = _partition_types.size() == 1 ? _partition_types.at(0) : _partition_types.at(i);
          partiapp = _partition_approach.size() == 1 ? _partition_approach.at(0) : _partition_approach.at(i);
          partiimbl = _partition_imbl_tol.size() == 1 ? _partition_imbl_tol.at(0) : _partition_imbl_tol.at(i);
          gsteps = _gsteps.size() == 1 ? _gsteps.at(0) : _gsteps.at(i);
          rsf_ao = _ao_rsf.size() == 1 ? _ao_rsf.at(0) : _ao_rsf.at(i);
          rsf_cm = _cm_rsf.size() == 1 ? _cm_rsf.at(0) : _cm_rsf.at(i);
          rsf_ext = _ext_rsf.size() == 1 ? _ext_rsf.at(0) : _ext_rsf.at(i);
          reuse_cm = _cm_r.size() == 1 ? _cm_r.at(0) : _cm_r.at(i);
          reuse_cb = _cb_r.size() == 1 ? _cb_r.at(0) : _cb_r.at(i);
          phi_dt = _phi_dt.size() == 1 ? _phi_dt.at(0) : _phi_dt.at(i);

          reuse = rsf_ao || rsf_cm || rsf_ext || reuse_cb || reuse_cm;

          parameterlist_header(_dimension, overlap, i+1, 1, _preconditioner, (i == 0 ? "Input" : "Laplace"), combine_lvl, reuse, *tmp_params);
          parameterlist_set_aoo(solver_ao, *tmp_params);
          parameterlist_combine_overlap(mode, tmp_params->sublist("AlgebraicOverlappingOperator"));
          tmp_params->sublist("AlgebraicOverlappingOperator").set("Reuse: Symbolic Factorization", rsf_ao);

          tmp_params->set("IPOUHarmonicCoarseOperator", Teuchos::ParameterList("IPOUHarmonicCoarseOperator"));
          tmp_params->sublist("IPOUHarmonicCoarseOperator").set("Blocks", Teuchos::ParameterList("Blocks"));
          parameterlist_set_reuse_coarse(rsf_cm, rsf_ext, reuse_cm, reuse_cb, tmp_params->sublist("IPOUHarmonicCoarseOperator"));
          tmp_params->sublist("IPOUHarmonicCoarseOperator").set("Phi: Dropping Threshold", phi_dt);

          if(_problem == FROSchParameterList::PRESSUREPOISSON)
          {
            /* for pressure poisson:  we always use the pressure in the coarse space and do not exclude the velocity (not existent) */
            /* pressure block */
            parameterlist_set_block(1, true, -1, ipou_pres, tmp_params->sublist("IPOUHarmonicCoarseOperator").sublist("Blocks"));
          }else if(_problem == FROSchParameterList::SADDLEPOINT)
          {
            /* velocity block */
            parameterlist_set_block(1, cspace_velo, (excl_velo_pres ? 2 : -1), ipou_velo, tmp_params->sublist("IPOUHarmonicCoarseOperator").sublist("Blocks"));
            /* pressure block */
            parameterlist_set_block(2, cspace_pres, (excl_pres_velo ? 1 : -1), ipou_pres, tmp_params->sublist("IPOUHarmonicCoarseOperator").sublist("Blocks"));

          }else
            XASSERTM(false, "Unkown FROSchParameterList::ProblemType type");

          parameterlist_set_extsolv(solver_ext, tmp_params->sublist("IPOUHarmonicCoarseOperator"));
          parameterlist_set_distribution(_subregions.at(i), partitype, partiapp, partiimbl, gsteps, tmp_params->sublist("IPOUHarmonicCoarseOperator"));

          /* next level */
          tmp_params->sublist("IPOUHarmonicCoarseOperator").set("CoarseSolver", Teuchos::ParameterList("CoarseSolver"));

          if(i+1 < _nlevels)
          {
            /* there is another level */
            tmp_params->sublist("IPOUHarmonicCoarseOperator").sublist("CoarseSolver").set("SolverType", "FROSchPreconditioner");
            tmp_params->sublist("IPOUHarmonicCoarseOperator").sublist("CoarseSolver").set("Solver", precond_str);
            tmp_params->sublist("IPOUHarmonicCoarseOperator").sublist("CoarseSolver").set("FROSchPreconditioner", Teuchos::ParameterList("FROSchPreconditioner"));
            tmp_params->sublist("IPOUHarmonicCoarseOperator").sublist("CoarseSolver").sublist("FROSchPreconditioner").set(precond_str, Teuchos::ParameterList(precond_str));

            tmp_params = &(tmp_params->sublist("IPOUHarmonicCoarseOperator").sublist("CoarseSolver").sublist("FROSchPreconditioner").sublist(precond_str));
          }else
          {
            /* add coarse solver */
            parameterlist_set_solver(_solver_coarse, tmp_params->sublist("IPOUHarmonicCoarseOperator").sublist("CoarseSolver"));
          }
        }

        _core = (void*) params;

        if(_print_list)
        {
          print();
        }

      }

      void FROSchParameterList::delete_core()
      {
        if(_core != nullptr)
        {
          // /* DEBUG */
          // if(_print_list)
          //   print();
          delete reinterpret_cast<Teuchos::ParameterList*>(_core);
        }
      }

      int FROSchParameterList::get_dim() const
      {
        return _dimension;
      }

      int FROSchParameterList::get_dpe_velo() const
      {
        return _dpe_velo;
      }

      int FROSchParameterList::get_dpe_pres() const
      {
        return _dpe_pres;
      }

      int FROSchParameterList::get_nlevels() const
      {
        return _nlevels;
      }

      bool FROSchParameterList::get_use_timer() const
      {
        return _use_timer;
      }

      bool FROSchParameterList::get_print() const
      {
        return _print_list;
      }

      void FROSchParameterList::set_print(const bool print_)
      {
        _print_list = print_;
      }

      void FROSchParameterList::set_use_timer(const bool use_timer_)
      {
        _use_timer = use_timer_;
      }

      void FROSchParameterList::set_nlevels(const int nlevels_)
      {
        _nlevels = nlevels_;
      }

      void FROSchParameterList::set_dim(const int dim_)
      {
        XASSERTM(dim_ == 2 || dim_ == 3, "Only implemented for dim_ == 2 or dim_ == 3");
        _dimension = dim_;
        _dpe_velo = dim_;
        _dpe_pres = dim_ + 1;
      }

      void FROSchParameterList::set_problem_type(const ProblemType problem_)
      {
        _problem = problem_;
      }

      void FROSchParameterList::set_precond(const Preconditioner precond_)
      {
        _preconditioner = precond_;
      }

      void FROSchParameterList::set_coarse_solver(const DirectSolver solver_)
      {
        _solver_coarse = solver_;
      }

      void FROSchParameterList::set_gsteps(const int gsteps_)
      {
        _gsteps.assign(_nlevels, gsteps_);
      }

      void FROSchParameterList::set_gsteps(const std::vector<int>& gsteps_)
      {
        XASSERT(int(gsteps_.size()) == _nlevels);
        _gsteps.assign(gsteps_.begin(), gsteps_.end());
      }

      void FROSchParameterList::set_overlaps(const int overlap_)
      {
        _overlaps.assign(_nlevels, overlap_);
      }

      void FROSchParameterList::set_overlaps(const std::vector<int>& overlaps_)
      {
        XASSERT(int(overlaps_.size()) == _nlevels);
        _overlaps.assign(overlaps_.begin(), overlaps_.end());
      }

      void FROSchParameterList::set_combine_overlap(const CombineOverlap combine_mode_)
      {
        _combine_ov.assign(_nlevels, combine_mode_);
      }

      void FROSchParameterList::set_combine_overlap(const std::vector<CombineOverlap>& combine_mode_)
      {
        XASSERT(int(combine_mode_.size()) == _nlevels);
        _combine_ov.assign(combine_mode_.begin(), combine_mode_.end());
      }

      void FROSchParameterList::set_combine_lvl(const CombineLevel combine_lvl_)
      {
        _combine_lvl.assign(_nlevels, combine_lvl_);
      }

      void FROSchParameterList::set_combine_lvl(const std::vector<CombineLevel>& combine_lvl_)
      {
        XASSERT(int(_combine_lvl.size()) == _nlevels);
        _combine_lvl.assign(combine_lvl_.begin(), combine_lvl_.end());
      }

      void FROSchParameterList::set_solvers(const DirectSolver solver_ao_, const DirectSolver solver_ext_)
      {
        set_solvers_ao(solver_ao_);
        set_solvers_ext(solver_ext_);
      }

      void FROSchParameterList::set_solvers(const std::vector<DirectSolver>& solver_ao_, const std::vector<DirectSolver>& solver_ext_)
      {
        set_solvers_ao(solver_ao_);
        set_solvers_ext(solver_ext_);
      }

      void FROSchParameterList::set_solvers_ao(const DirectSolver solver_)
      {
        _solvers_ao.assign(_nlevels, solver_);
      }

      void FROSchParameterList::set_solvers_ao(const std::vector<DirectSolver>& solvers_)
      {
        XASSERT(int(solvers_.size()) == _nlevels);
        _solvers_ao.assign(solvers_.begin(), solvers_.end());
      }

      void FROSchParameterList::set_solvers_ext(const DirectSolver solver_)
      {
        _solvers_ext.assign(_nlevels, solver_);
      }

      void FROSchParameterList::set_solvers_ext(const std::vector<DirectSolver>& solvers_)
      {
        XASSERT(int(solvers_.size()) == _nlevels);
        _solvers_ext.assign(solvers_.begin(), solvers_.end());
      }

      void FROSchParameterList::set_paritioner(const std::vector<int>& subregions_, const PartitionType parti_type_, const PartitionApproach parti_approach_, const double parti_imbl_)
      {
        set_subregions(subregions_);
        set_parti_types(parti_type_);
        set_parti_approach(parti_approach_);
        set_parti_imbl(parti_imbl_);
      }

      void FROSchParameterList::set_paritioner(const std::vector<int>& subregions_, const std::vector<PartitionType>& parti_types_, const std::vector<PartitionApproach>& parti_approach_, const std::vector<double>& parti_imbl_)
      {
        set_subregions(subregions_);
        set_parti_types(parti_types_);
        set_parti_approach(parti_approach_);
        set_parti_imbl(parti_imbl_);
      }

      void FROSchParameterList::set_subregions(const std::vector<int>& subregions_)
      {
        XASSERT(int(subregions_.size()) == _nlevels);
        XASSERTM(int(subregions_.back()) == int(1), "Last entry has to be one.");
        _subregions.assign(subregions_.begin(), subregions_.end());
      }

      void FROSchParameterList::set_parti_types(const PartitionType parti_type_)
      {
        _partition_types.assign(_nlevels, parti_type_);
      }

      void FROSchParameterList::set_parti_types(const std::vector<PartitionType>& parti_types_)
      {
        XASSERT(int(parti_types_.size()) == _nlevels);
        _partition_types.assign(parti_types_.begin(), parti_types_.end());
      }

      void FROSchParameterList::set_parti_approach(const PartitionApproach parti_approach_)
      {
        _partition_approach.assign(_nlevels, parti_approach_);
      }

      void FROSchParameterList::set_parti_approach(const std::vector<PartitionApproach>& parti_approach_)
      {
        XASSERT(int(parti_approach_.size()) == _nlevels);
        _partition_approach.assign(parti_approach_.begin(), parti_approach_.end());
      }

      void FROSchParameterList::set_parti_imbl(const double parti_imbl_)
      {
        _partition_imbl_tol.assign(_nlevels, parti_imbl_);
      }

      void FROSchParameterList::set_parti_imbl(const std::vector<double>& parti_imbl_)
      {
        XASSERT(int(parti_imbl_.size()) == _nlevels);
        _partition_imbl_tol.assign(parti_imbl_.begin(), parti_imbl_.end());
      }

      void FROSchParameterList::set_excludes(const bool velo_pres_, const bool pres_velo_)
      {
        _exclude_velo_pres.assign(_nlevels, velo_pres_);
        _exclude_pres_velo.assign(_nlevels, pres_velo_);
      }

      void FROSchParameterList::set_excludes(const std::vector<bool>& velo_pres_, const std::vector<bool>& pres_velo_)
      {
        XASSERT(int(velo_pres_.size()) == _nlevels);
        XASSERT(int(pres_velo_.size()) == _nlevels);
        _exclude_velo_pres.assign(velo_pres_.begin(), velo_pres_.end());
        _exclude_pres_velo.assign(pres_velo_.begin(), pres_velo_.end());
      }

      void FROSchParameterList::set_excludes_velo_pres(const bool velo_pres_)
      {
        _exclude_velo_pres.assign(_nlevels, velo_pres_);
      }

      void FROSchParameterList::set_excludes_velo_pres(const std::vector<bool>& velo_pres_)
      {
        XASSERT(int(velo_pres_.size()) == _nlevels);
        _exclude_velo_pres.assign(velo_pres_.begin(), velo_pres_.end());
      }

      void FROSchParameterList::set_excludes_pres_velo(const bool pres_velo_)
      {
        _exclude_pres_velo.assign(_nlevels, pres_velo_);
      }

      void FROSchParameterList::set_excludes_pres_velo(const std::vector<bool>& pres_velo_)
      {
        XASSERT(int(pres_velo_.size()) == _nlevels);
        _exclude_pres_velo.assign(pres_velo_.begin(), pres_velo_.end());
      }

      void FROSchParameterList::set_ipous(const IPOU ipou_velo_, const IPOU ipou_pres_)
      {
        set_ipous_velo(ipou_velo_);
        set_ipous_pres(ipou_pres_);
      }

      void FROSchParameterList::set_ipous(const std::vector<IPOU>& ipous_velo_, const std::vector<IPOU>& ipous_pres_)
      {
        set_ipous_velo(ipous_velo_);
        set_ipous_pres(ipous_pres_);
      }

      void FROSchParameterList::set_ipous_velo(const IPOU ipou_)
      {
        _ipous_velo.assign(_nlevels, ipou_);
      }

      void FROSchParameterList::set_ipous_velo(const std::vector<IPOU>& ipous_)
      {
        XASSERT(int(ipous_.size()) == _nlevels);
        _ipous_velo.assign(ipous_.begin(), ipous_.end());
      }

      void FROSchParameterList::set_ipous_pres(const IPOU ipou_)
      {
        _ipous_pres.assign(_nlevels, ipou_);
      }

      void FROSchParameterList::set_ipous_pres(const std::vector<IPOU>& ipous_)
      {
        XASSERT(int(ipous_.size()) == _nlevels);
        _ipous_pres.assign(ipous_.begin(), ipous_.end());
      }

      void FROSchParameterList::set_cspace(const bool cspace_velo_, const bool cspace_pres_)
      {
        set_cspace_velo(cspace_velo_);
        set_cspace_pres(cspace_pres_);
      }

      void FROSchParameterList::set_cspace(const std::vector<bool>& cspace_velo_, const std::vector<bool>& cspace_pres_)
      {
        set_cspace_velo(cspace_velo_);
        set_cspace_pres(cspace_pres_);
      }

      void FROSchParameterList::set_cspace_velo(const bool cspace_)
      {
        _cspace_velo.assign(_nlevels, cspace_);
      }

      void FROSchParameterList::set_cspace_velo(const std::vector<bool>& cspace_)
      {
        XASSERT(int(cspace_.size()) == _nlevels);
        _cspace_velo.assign(cspace_.begin(), cspace_.end());
      }

      void FROSchParameterList::set_cspace_pres(const bool cspace_)
      {
        _cspace_pres.assign(_nlevels, cspace_);
      }

      void FROSchParameterList::set_cspace_pres(const std::vector<bool>& cspace_)
      {
        XASSERT(int(cspace_.size()) == _nlevels);
        _cspace_pres.assign(cspace_.begin(), cspace_.end());
      }

      void FROSchParameterList::set_reuse_sf(const bool rsf_ao_, const bool rsf_cm_, const bool rsf_ext_)
      {
        set_reuse_sf_ao(rsf_ao_);
        set_reuse_sf_cm(rsf_cm_);
        set_reuse_sf_ext(rsf_ext_);
      }

      void FROSchParameterList::set_reuse_sf(const std::vector<bool> &rsf_ao_, const std::vector<bool> &rsf_cm_, const std::vector<bool> &rsf_ext_)
      {
        set_reuse_sf_ao(rsf_ao_);
        set_reuse_sf_cm(rsf_cm_);
        set_reuse_sf_ext(rsf_ext_);
      }

      void FROSchParameterList::set_reuse_sf_ao(const bool rsf_ao_)
      {
        _ao_rsf.assign(_nlevels, rsf_ao_);
      }

      void FROSchParameterList::set_reuse_sf_ao(const std::vector<bool> &rsf_ao_)
      {
        XASSERT(int(rsf_ao_.size()) == _nlevels);
        _ao_rsf.assign(rsf_ao_.begin(), rsf_ao_.end());
      }

      void FROSchParameterList::set_reuse_sf_cm(const bool rsf_cm_)
      {
        _cm_rsf.assign(_nlevels, rsf_cm_);
      }

      void FROSchParameterList::set_reuse_sf_cm(const std::vector<bool> &rsf_cm_)
      {
        XASSERT(int(rsf_cm_.size()) == _nlevels);
        _cm_rsf.assign(rsf_cm_.begin(), rsf_cm_.end());
      }

      void FROSchParameterList::set_reuse_sf_ext(const bool rsf_ext_)
      {
        _ext_rsf.assign(_nlevels, rsf_ext_);
      }

      void FROSchParameterList::set_reuse_sf_ext(const std::vector<bool> &rsf_ext_)
      {
        XASSERT(int(rsf_ext_.size()) == _nlevels);
        _ext_rsf.assign(rsf_ext_.begin(), rsf_ext_.end());
      }

      void FROSchParameterList::set_reuse_coarse(const bool reuse_cm_, const bool reuse_cb_)
      {
        set_reuse_coarse_matrix(reuse_cm_);
        set_reuse_coarse_basis(reuse_cb_);
      }

      void FROSchParameterList::set_reuse_coarse(const std::vector<bool> &reuse_cm_, const std::vector<bool> &reuse_cb_)
      {
        set_reuse_coarse_matrix(reuse_cm_);
        set_reuse_coarse_basis(reuse_cb_);
      }

      void FROSchParameterList::set_reuse_coarse_matrix(const bool reuse_cm_)
      {
        _cm_r.assign(_nlevels, reuse_cm_);
      }

      void FROSchParameterList::set_reuse_coarse_matrix(const std::vector<bool> &reuse_cm_)
      {
        XASSERT(int(reuse_cm_.size()) == _nlevels);
        _cm_r.assign(reuse_cm_.begin(), reuse_cm_.end());
      }

      void FROSchParameterList::set_reuse_coarse_basis(const bool rsf_cm_)
      {
        _cb_r.assign(_nlevels, rsf_cm_);
      }

      void FROSchParameterList::set_reuse_coarse_basis(const std::vector<bool> &reuse_cb_)
      {
        XASSERT(int(reuse_cb_.size()) == _nlevels);
        _cb_r.assign(reuse_cb_.begin(), reuse_cb_.end());
      }

      void FROSchParameterList::set_phi_dropping_threshold(const double phi_dt_)
      {
        _phi_dt.assign(_nlevels, phi_dt_);
      }

      void FROSchParameterList::set_phi_dropping_threshold(const std::vector<double> &phi_dt_)
      {
        XASSERT(int(phi_dt_.size()) == _nlevels);
        _phi_dt.assign(phi_dt_.begin(), phi_dt_.end());
      }

      inline bool parse_fpl_levels(const Dist::Comm &comm, SimpleArgParser &args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-nlevels");
        if(it != nullptr)
        {
          if(it->second.size() == 1)
          {
            int nlevels = std::stoi(it->second.front());
            if(nlevels < 1)
            {
              comm.print("ERROR: given number of levels is smaller then 1");
            }
            params.set_nlevels(nlevels);
          }else
          {
            comm.print("ERROR: multiple levels given");
            return false;
          }
        }
        return true;
      }

      inline bool parse_fpl_dim(const Dist::Comm& comm, SimpleArgParser& args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-dim");
        if(it != nullptr)
        {
          if(it->second.size() == 1)
          {
            int dim = std::stoi(it->second.front());
            if(dim != 2 && dim != 3)
            {
              comm.print("ERROR: given dimension has to be 2 or 3");
            }
            params.set_dim(dim);
          }else
          {
            comm.print("ERROR: multiple dimensions given");
            return false;
          }
        }

        return true;
      }

      inline bool parse_fpl_overlaps(const Dist::Comm& comm, SimpleArgParser& args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-overlap");
        if(it != nullptr)
        {
          if(it->second.size() == 1u)
          {
            auto ov = std::stoi(it->second.front());
            if(ov < 0 )
            {
              comm.print("ERROR: overlap is smaller than 0");
              return false;
            }
            params.set_overlaps((int)ov);
          }else if(int(it->second.size()) == params.get_nlevels())
          {
            std::vector<int> overlaps;
            overlaps.reserve(params.get_nlevels());
            for(const auto& t : it->second)
            {
              auto ov = std::stoi(t);
              if(ov < 0)
              {
                comm.print("ERROR: overlap is smaller than 0");
                return false;
              }
              overlaps.push_back(ov);
            }
            params.set_overlaps(overlaps);
          }else
          {
            comm.print("ERROR: not matching number of overlaps: nlevels:  " + std::to_string(params.get_nlevels()) + " number of given overlaps: " + std::to_string(it->second.size()) + "'");
            return false;
          }
        }
        return true;
      }

      inline bool parse_fpl_subreg(const Dist::Comm& comm, SimpleArgParser& args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-subreg");
        if(it != nullptr)
        {
          if(int(it->second.size()) == params.get_nlevels())
          {
            std::vector<int> subreg;
            subreg.reserve(params.get_nlevels());
            for(const auto& t : it->second)
            {
              auto sr = std::stoi(t);
              if(sr < 0)
              {
                /* __TODO:__ Eigentlich müsste man hier testen auf "< 1", aber der letzte Eintrag muss 1 sein */
                comm.print("ERROR: subreg is smaller than 0");
                return false;
              }
              subreg.push_back(sr);
            }
            if(subreg.back() != 1)
            {
              comm.print("ERROR: last entry of subreg is not 1: " + std::to_string(subreg.back()));
            }
            params.set_subregions(subreg);
          }else
          {
            comm.print("ERROR: not matching number of subregions: nlevels:  " + std::to_string(params.get_nlevels()) + " number of given subregions: " + std::to_string(it->second.size()) + "'");
            return false;
          }
        }
        return true;
      }

      inline bool parse_fpl_gsteps(const Dist::Comm& comm, SimpleArgParser& args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-gstep");
        if(it != nullptr)
        {
          if(it->second.size() == 1u)
          {
            auto gsteps = std::stoi(it->second.front());
            if(gsteps < 1 )
            {
              comm.print("ERROR: gstep is smaller than 1");
              return false;
            }
            params.set_gsteps((int)gsteps);
          }else if(int(it->second.size()) == params.get_nlevels())
          {
            std::vector<int> gsteps;
            gsteps.reserve(params.get_nlevels());
            for(const auto& t : it->second)
            {
              auto gs = std::stoi(t);
              if(gs < 1)
              {
                comm.print("ERROR: gstep is smaller than 0");
                return false;
              }
              gsteps.push_back(gs);
            }
            params.set_gsteps(gsteps);
          }else
          {
            comm.print("ERROR: not matching number of gsteps: nlevels:  " + std::to_string(params.get_nlevels()) + " number of given gsteps: " + std::to_string(it->second.size()) + "'");
            return false;
          }
        }
        return true;
      }

      inline bool parse_fpl_parti_type(const Dist::Comm& comm, SimpleArgParser& args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-parti-type");
        if(it != nullptr)
        {
          if(it->second.size() == 1u)
          {
            FROSchParameterList::PartitionType parti_type;
            std::string type = it->second.front();
            std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

            if(!type.compare("parmetis"))
            {
              parti_type = Solver::Trilinos::FROSchParameterList::PARMETIS;
            }else if(!type.compare("block"))
            {
              parti_type = Solver::Trilinos::FROSchParameterList::BLOCK;
            }else if(!type.compare("phg"))
            {
              parti_type = Solver::Trilinos::FROSchParameterList::PHG;
            }else
            {
              comm.print("ERROR: unknown partitioner type '" + it->second.front() + "'");
              return false;
            }
            params.set_parti_types(parti_type);
          }else if(int(it->second.size()) == params.get_nlevels())
          {
            std::vector<FROSchParameterList::PartitionType> parti_types;
            parti_types.reserve(params.get_nlevels());
            for(const auto& t : it->second)
            {
              std::string type = t;
              std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

              if(!type.compare("parmetis"))
              {
                parti_types.push_back(FROSchParameterList::PARMETIS);
              }else if(!type.compare("block"))
              {
                parti_types.push_back(FROSchParameterList::BLOCK);
              }else if(!type.compare("phg"))
              {
                parti_types.push_back(Solver::Trilinos::FROSchParameterList::PHG);
              }else
              {
                comm.print("ERROR: unknown partitioner type '" + t + "'");
                return false;
              }
            }
            params.set_parti_types(parti_types);
          }else
          {
            comm.print("ERROR: not matching number of parti-type: nlevels:  " + std::to_string(params.get_nlevels()) + " number of given parti-type: " + std::to_string(it->second.size()) + "'");
            return false;
          }
        }
        return true;
      }

      inline bool parse_fpl_parti_approach(const Dist::Comm& comm, SimpleArgParser& args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-parti-approach");
        if(it != nullptr)
        {
          if(it->second.size() == 1u)
          {
            FROSchParameterList::PartitionApproach parti_app;
            std::string type = it->second.front();
            std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

            if(!type.compare("partition"))
            {
              parti_app = FROSchParameterList::PARTITION;
            }else if(!type.compare("repartition"))
            {
              parti_app = FROSchParameterList::REPARTITION;
            }else
            {
              comm.print("ERROR: unknown partitioner approach '" + it->second.front() + "'");
              return false;
            }
            params.set_parti_approach(parti_app);
          }else if(int(it->second.size()) == params.get_nlevels())
          {
            std::vector<FROSchParameterList::PartitionApproach> parti_apps;
            parti_apps.reserve(params.get_nlevels());
            for(const auto& t : it->second)
            {
              std::string type = t;
              std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

              if(!type.compare("partition"))
              {
                parti_apps.push_back(FROSchParameterList::PARTITION);
              }else if(!type.compare("repartitionblock"))
              {
                parti_apps.push_back(FROSchParameterList::REPARTITION);
              }else
              {
                comm.print("ERROR: unknown partitioner approach '" + t + "'");
                return false;
              }
            }
            params.set_parti_approach(parti_apps);
          }else
          {
            comm.print("ERROR: not matching number of parti-approach: nlevels:  " + std::to_string(params.get_nlevels()) + " number of given parti-approach: " + std::to_string(it->second.size()) + "'");
            return false;
          }
        }
        return true;
      }

      inline bool parse_fpl_parti_imbl(const Dist::Comm& comm, SimpleArgParser& args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-parti-imbl");
        if(it != nullptr)
        {
          if(it->second.size() == 1u)
          {
            double parti_imbl = std::stod(it->second.front());
            if(parti_imbl < 1.)
            {
              /* __TODO:__ Stimmt das? */
              comm.print("ERROR: unknown partition imblance smaller than 1.:  '" + it->second.front() + "'");
              return false;
            }
            params.set_parti_imbl(parti_imbl);
          }else if(int(it->second.size()) == params.get_nlevels())
          {
            std::vector<double> parti_imbl;
            parti_imbl.reserve(params.get_nlevels());
            for(const auto& t : it->second)
            {
              double imbl = std::stod(t);
              if(imbl < 1.)
              {
                /* __TODO:__ Stimmt das? */
                comm.print("ERROR: unknown partition imblance smaller than 1.:  '" + t + "'");
                return false;
              }
              parti_imbl.push_back(imbl);
            }
            params.set_parti_imbl(parti_imbl);
          }else
          {
            comm.print("ERROR: not matching number of parti-imbl: nlevels:  " + std::to_string(params.get_nlevels()) + " number of given parti-imbl: " + std::to_string(it->second.size()) + "'");
            return false;
          }
        }
        return true;
      }

      inline bool parse_fpl_ipou_pres(const Dist::Comm& comm, SimpleArgParser& args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-ipou-pres");
        if(it != nullptr)
        {
          if(it->second.size() == 1u)
          {
            FROSchParameterList::IPOU ipou;
            std::string type = it->second.front();
            std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

            if(!type.compare("gdsw"))
            {
              ipou = FROSchParameterList::GDSW;
            }else if(!type.compare("gdswstar"))
            {
              ipou = FROSchParameterList::GDSWSTAR;
            }else if(!type.compare("rgdsw"))
            {
              ipou = FROSchParameterList::RGDSW;
            }else
            {
              comm.print("ERROR: unknown ipou pressure type '" + it->second.front() + "'");
              return false;
            }
            params.set_ipous_pres(ipou);
          }else if(int(it->second.size()) == params.get_nlevels())
          {
            std::vector<FROSchParameterList::IPOU> ipous;
            ipous.reserve(params.get_nlevels());
            for(const auto& t : it->second)
            {
              std::string type = t;
              std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

              if(!type.compare("gdsw"))
              {
                ipous.push_back(FROSchParameterList::GDSW);
              }else if(!type.compare("gdswstar"))
              {
                ipous.push_back(FROSchParameterList::GDSWSTAR);
              }else if(!type.compare("rgdsw"))
              {
                ipous.push_back(FROSchParameterList::RGDSW);
              }else
              {
                comm.print("ERROR: unknown ipou pressure type '" + t + "'");
                return false;
              }
            }
            params.set_ipous_pres(ipous);
          }else
          {
            comm.print("ERROR: not matching number of ipou-pres: nlevels:  " + std::to_string(params.get_nlevels()) + " number of given ipou-pres: " + std::to_string(it->second.size()) + "'");
            return false;
          }
        }
        return true;
      }

      inline bool parse_fpl_ipou_velo(const Dist::Comm& comm, SimpleArgParser& args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-ipou-velo");
        if(it != nullptr)
        {
          if(it->second.size() == 1u)
          {
            FROSchParameterList::IPOU ipou;
            std::string type = it->second.front();
            std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

            if(!type.compare("gdsw"))
            {
              ipou = FROSchParameterList::GDSW;
            }else if(!type.compare("gdswstar"))
            {
              ipou = FROSchParameterList::GDSWSTAR;
            }else if(!type.compare("rgdsw"))
            {
              ipou = FROSchParameterList::RGDSW;
            }else
            {
              comm.print("ERROR: unknown ipou velocity type '" + it->second.front() + "'");
              return false;
            }
            params.set_ipous_velo(ipou);
          }else if(int(it->second.size()) == params.get_nlevels())
          {
            std::vector<FROSchParameterList::IPOU> ipous;
            ipous.reserve(params.get_nlevels());
            for(const auto& t : it->second)
            {
              std::string type = t;
              std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

              if(!type.compare("gdsw"))
              {
                ipous.push_back(FROSchParameterList::GDSW);
              }else if(!type.compare("gdswstar"))
              {
                ipous.push_back(FROSchParameterList::GDSWSTAR);
              }else if(!type.compare("rgdsw"))
              {
                ipous.push_back(FROSchParameterList::RGDSW);
              }else
              {
                comm.print("ERROR: unknown ipou velo type '" + t + "'");
                return false;
              }
            }
            params.set_ipous_pres(ipous);
          }else
          {
            comm.print("ERROR: not matching number of ipou-velo: nlevels:  " + std::to_string(params.get_nlevels()) + " number of given ipou-velo: " + std::to_string(it->second.size()) + "'");
            return false;
          }
        }
        return true;
      }

      inline bool parse_fpl_solver_ao(const Dist::Comm& comm, SimpleArgParser& args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-solver-ao");
        if(it != nullptr)
        {
          if(it->second.size() == 1u)
          {
            FROSchParameterList::DirectSolver solver;
            std::string type = it->second.front();
            std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

            if(!type.compare("klu"))
            {
              solver = FROSchParameterList::KLU;
            }else if(!type.compare("mumps"))
            {
              solver = FROSchParameterList::MUMPS;
            }else if(!type.compare("umfpack"))
            {
              solver = FROSchParameterList::UMFPACK;
            }else
            {
              comm.print("ERROR: unknown solver-ao type '" + it->second.front() + "'");
              return false;
            }
            params.set_solvers_ao(solver);
          }else if(int(it->second.size()) == params.get_nlevels())
          {
            std::vector<FROSchParameterList::DirectSolver> solvers;
            solvers.reserve(params.get_nlevels());
            for(const auto& t : it->second)
            {
              std::string type = t;
              std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

              if(!type.compare("klu"))
              {
                solvers.push_back(FROSchParameterList::KLU);
              }else if(!type.compare("mumps"))
              {
                solvers.push_back(FROSchParameterList::MUMPS);
              }else if(!type.compare("umfpack"))
              {
                solvers.push_back(FROSchParameterList::UMFPACK);
              }else
              {
                comm.print("ERROR: unknown solver-ao type '" + t + "'");
                return false;
              }
            }
            params.set_solvers_ao(solvers);
          }else
          {
            comm.print("ERROR: not matching number of solver-ao: nlevels:  " + std::to_string(params.get_nlevels()) + " number of given solver-ao: " + std::to_string(it->second.size()) + "'");
            return false;
          }
        }
        return true;
      }

      inline bool parse_fpl_solver_ext(const Dist::Comm& comm, SimpleArgParser& args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-solver-ext");
        if(it != nullptr)
        {
          if(it->second.size() == 1u)
          {
            FROSchParameterList::DirectSolver solver;
            std::string type = it->second.front();
            std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

            if(!type.compare("klu"))
            {
              solver = FROSchParameterList::KLU;
            }else if(!type.compare("mumps"))
            {
              solver = FROSchParameterList::MUMPS;
            }else if(!type.compare("umfpack"))
            {
              solver = FROSchParameterList::UMFPACK;
            }else
            {
              comm.print("ERROR: unknown solver-ext type '" + it->second.front() + "'");
              return false;
            }
            params.set_solvers_ext(solver);
          }else if(int(it->second.size()) == params.get_nlevels())
          {
            std::vector<FROSchParameterList::DirectSolver> solvers;
            solvers.reserve(params.get_nlevels());
            for(const auto& t : it->second)
            {
              std::string type = t;
              std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

              if(!type.compare("klu"))
              {
                solvers.push_back(FROSchParameterList::KLU);
              }else if(!type.compare("mumps"))
              {
                solvers.push_back(FROSchParameterList::MUMPS);
              }else if(!type.compare("umfpack"))
              {
                solvers.push_back(FROSchParameterList::UMFPACK);
              }else
              {
                comm.print("ERROR: unknown solver-ext type '" + t + "'");
                return false;
              }
            }
            params.set_solvers_ext(solvers);
          }else
          {
            comm.print("ERROR: not matching number of solver-ext: nlevels:  " + std::to_string(params.get_nlevels()) + " number of given solver-ext: " + std::to_string(it->second.size()) + "'");
            return false;
          }
        }
        return true;
      }

      inline bool parse_fpl_solver_coarse(const Dist::Comm& comm, SimpleArgParser& args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-solver-coarse");
        if(it != nullptr)
        {
          if(it->second.size() == 1u)
          {
            FROSchParameterList::DirectSolver solver;
            std::string type = it->second.front();
            std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

            if(!type.compare("klu"))
            {
              solver = FROSchParameterList::KLU;
            }else if(!type.compare("mumps"))
            {
              solver = FROSchParameterList::MUMPS;
            }else if(!type.compare("umfpack"))
            {
              solver = FROSchParameterList::UMFPACK;
            }else
            {
              comm.print("ERROR: unknown solver-coarse type '" + it->second.front() + "'");
              return false;
            }
            params.set_coarse_solver(solver);
          }else
          {
            comm.print("ERROR: multiple coarse solvers are given");
            return false;
          }
        }
        return true;
      }

      inline bool parse_fpl_precond(const Dist::Comm& comm, SimpleArgParser& args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-precond-type");
        if(it != nullptr)
        {
          if(it->second.size() == 1u)
          {
            FROSchParameterList::Preconditioner precond;
            std::string type = it->second.front();
            std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

            if(!type.compare("twolevel"))
            {
              precond = FROSchParameterList::TWOLEVEL;
            }else if(!type.compare("twolevelblock"))
            {
              precond = FROSchParameterList::TWOLEVELBLOCK;
            }else if(!type.compare("onelevel"))
            {
              precond = FROSchParameterList::ONELEVEL;
            }else
            {
              comm.print("ERROR: unknown precond-type '" + it->second.front() + "'");
              return false;
            }
            params.set_precond(precond);
          }else
          {
            comm.print("ERROR: multiple preconditioner types are given");
            return false;
          }
        }
        return true;
      }

      inline bool parse_fpl_cspace_pres(const Dist::Comm& comm, SimpleArgParser& args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-cspace-pres");
        if(it != nullptr)
        {
          if(it->second.size() == 1u)
          {
            bool cspace;
            std::string type = it->second.front();
            std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

            if(!type.compare("true"))
            {
              cspace = true;
            }else if(!type.compare("false"))
            {
              cspace = false;
            }else
            {
              comm.print("ERROR: unknown cspace-pres type '" + it->second.front() + "'");
              return false;
            }
            params.set_cspace_pres(cspace);
          }else if(int(it->second.size()) == params.get_nlevels())
          {
            std::vector<bool> cspaces;
            cspaces.reserve(params.get_nlevels());
            for(const auto& t : it->second)
            {
              std::string type = t;
              std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

              if(!type.compare("true"))
              {
                cspaces.push_back(true);
              }else if(!type.compare("false"))
              {
                cspaces.push_back(false);
              }else
              {
                comm.print("ERROR: unknown cspace-pres type '" + t + "'");
                return false;
              }
            }
            params.set_cspace_pres(cspaces);
          }else
          {
            comm.print("ERROR: not matching number of cspace-pres: nlevels:  " + std::to_string(params.get_nlevels()) + " number of given cspace-pres: " + std::to_string(it->second.size()) + "'");
            return false;
          }
        }
        return true;
      }

      inline bool parse_fpl_cspace_velo(const Dist::Comm& comm, SimpleArgParser& args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-cspace-velo");
        if(it != nullptr)
        {
          if(it->second.size() == 1u)
          {
            bool cspace;
            std::string type = it->second.front();
            std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

            if(!type.compare("true"))
            {
              cspace = true;
            }else if(!type.compare("false"))
            {
              cspace = false;
            }else
            {
              comm.print("ERROR: unknown cspace-velo type '" + it->second.front() + "'");
              return false;
            }
            params.set_cspace_velo(cspace);
          }else if(int(it->second.size()) == params.get_nlevels())
          {
            std::vector<bool> cspaces;
            cspaces.reserve(params.get_nlevels());
            for(const auto& t : it->second)
            {
              std::string type = t;
              std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

              if(!type.compare("true"))
              {
                cspaces.push_back(true);
              }else if(!type.compare("false"))
              {
                cspaces.push_back(false);
              }else
              {
                comm.print("ERROR: unknown cspace-velo type '" + t + "'");
                return false;
              }
            }
            params.set_cspace_velo(cspaces);
          }else
          {
            comm.print("ERROR: not matching number of cspace-velo: nlevels:  " + std::to_string(params.get_nlevels()) + " number of given cspace-velo: " + std::to_string(it->second.size()) + "'");
            return false;
          }
        }
        return true;
      }

      inline bool parse_fpl_exclude_velo_pres(const Dist::Comm& comm, SimpleArgParser& args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-exclude-velo-pres");
        if(it != nullptr)
        {
          if(it->second.size() == 1u)
          {
            bool excl;
            std::string type = it->second.front();
            std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

            if(!type.compare("true"))
            {
              excl = true;
            }else if(!type.compare("false"))
            {
              excl = false;
            }else
            {
              comm.print("ERROR: unknown exclude-velo-pres type '" + it->second.front() + "'");
              return false;
            }
            params.set_excludes_velo_pres(excl);
          }else if(int(it->second.size()) == params.get_nlevels())
          {
            std::vector<bool> excl;
            excl.reserve(params.get_nlevels());
            for(const auto& t : it->second)
            {
              std::string type = t;
              std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

              if(!type.compare("true"))
              {
                excl.push_back(true);
              }else if(!type.compare("false"))
              {
                excl.push_back(false);
              }else
              {
                comm.print("ERROR: unknown exclude-velo-pres type '" + t + "'");
                return false;
              }
            }
            params.set_excludes_velo_pres(excl);
          }else
          {
            comm.print("ERROR: not matching number of exclude-velo-pres: nlevels:  " + std::to_string(params.get_nlevels()) + " number of given exclude-velo-pres: " + std::to_string(it->second.size()) + "'");
            return false;
          }
        }
        return true;
      }

      inline bool parse_fpl_exclude_pres_velo(const Dist::Comm& comm, SimpleArgParser& args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-exclude-pres-velo");
        if(it != nullptr)
        {
          if(it->second.size() == 1u)
          {
            bool excl;
            std::string type = it->second.front();
            std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

            if(!type.compare("true"))
            {
              excl = true;
            }else if(!type.compare("false"))
            {
              excl = false;
            }else
            {
              comm.print("ERROR: unknown exclude-velo-pres type '" + it->second.front() + "'");
              return false;
            }
            params.set_excludes_pres_velo(excl);
          }else if(int(it->second.size()) == params.get_nlevels())
          {
            std::vector<bool> excl;
            excl.reserve(params.get_nlevels());
            for(const auto& t : it->second)
            {
              std::string type = t;
              std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

              if(!type.compare("true"))
              {
                excl.push_back(true);
              }else if(!type.compare("false"))
              {
                excl.push_back(false);
              }else
              {
                comm.print("ERROR: unknown exclude-pres-velo type '" + t + "'");
                return false;
              }
            }
            params.set_excludes_pres_velo(excl);
          }else
          {
            comm.print("ERROR: not matching number of exclude-pres-velo: nlevels:  " + std::to_string(params.get_nlevels()) + " number of given exclude-pres-velo: " + std::to_string(it->second.size()) + "'");
            return false;
          }
        }
        return true;
      }

      inline bool parse_fpl_print(const Dist::Comm& comm, SimpleArgParser& args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-print-list");
        if(it != nullptr)
        {
          if(it->second.size() == 1u)
          {
            bool print;
            std::string type = it->second.front();
            std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });
            if(!type.compare("true"))
            {
              print = true;
            }else if(!type.compare("false"))
            {
              print = false;
            }else
            {
              comm.print("ERROR: unknown print-list '" + it->second.front() + "'");
              return false;
            }
            params.set_print(print);
          }else
          {
            comm.print("ERROR: more than one print-list arg given: number of args: " + std::to_string(it->second.size()) + "'");
            return false;
          }
        }
        return true;
      }

      inline bool parse_fpl_timer(const Dist::Comm& comm, SimpleArgParser& args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-use-timer");
        if(it != nullptr)
        {
          if(it->second.size() == 1u)
          {
            bool timer;
            std::string type = it->second.front();
            std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

            if(!type.compare("true"))
            {
              timer = true;
            }else if(!type.compare("false"))
            {
              timer = false;
            }else
            {
              comm.print("ERROR: unknown use-timer '" + it->second.front() + "'");
              return false;
            }
            params.set_use_timer(timer);
          }else
          {
            comm.print("ERROR: more than one use-timer arg given: number of args: " + std::to_string(it->second.size()) + "'");
            return false;
          }
        }
        return true;
      }

      inline bool parse_fpl_combine(const Dist::Comm& comm, SimpleArgParser& args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-combine-overlap");
        if(it != nullptr)
        {
          if(it->second.size() == 1u)
          {
            FROSchParameterList::CombineOverlap mode;
            std::string type = it->second.front();
            std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

            if(!type.compare("full"))
            {
              mode = FROSchParameterList::FULL;
            }else if(!type.compare("restricted"))
            {
              mode = FROSchParameterList::RESTRICTED;
            }else if(!type.compare("averaging"))
            {
              mode = FROSchParameterList::AVERAGING;
            }else
            {
              comm.print("ERROR: unknown combine-overlap type '" + it->second.front() + "'");
              return false;
            }
            params.set_combine_overlap(mode);
          }else if(int(it->second.size()) == params.get_nlevels())
          {
            std::vector<FROSchParameterList::CombineOverlap> modes;
            modes.reserve(params.get_nlevels());
            for(const auto& t : it->second)
            {
              std::string type = t;
              std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

              if(!type.compare("full"))
              {
                modes.push_back(FROSchParameterList::FULL);
              }else if(!type.compare("restricted"))
              {
                modes.push_back(FROSchParameterList::RESTRICTED);
              }else if(!type.compare("umfpack"))
              {
                modes.push_back(FROSchParameterList::AVERAGING);
              }else
              {
                comm.print("ERROR: unknown combine-overlap type '" + t + "'");
                return false;
              }
            }
            params.set_combine_overlap(modes);
          }else
          {
            comm.print("ERROR: not matching number of combine-overlap: nlevels:  " + std::to_string(params.get_nlevels()) + " number of given combine-overlap: " + std::to_string(it->second.size()) + "'");
            return false;
          }
        }
        return true;
      }

      inline bool parse_fpl_combine_lvl(const Dist::Comm& comm, SimpleArgParser& args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-combine-lvl");
        if(it != nullptr)
        {
          if(it->second.size() == 1u)
          {
            FROSchParameterList::CombineLevel mode;
            std::string type = it->second.front();
            std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

            if(!type.compare("additive"))
            {
              mode = FROSchParameterList::ADDITIVE;
            }else if(!type.compare("multiplicative"))
            {
              mode = FROSchParameterList::MULTIPLICATIVE;
            }else
            {
              comm.print("ERROR: unknown combine-lvl type '" + it->second.front() + "'");
              return false;
            }
            params.set_combine_lvl(mode);
          }else if(int(it->second.size()) == params.get_nlevels())
          {
            std::vector<FROSchParameterList::CombineLevel> modes;
            modes.reserve(params.get_nlevels());
            for(const auto& t : it->second)
            {
              std::string type = t;
              std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

              if(!type.compare("additive"))
              {
                modes.push_back(FROSchParameterList::ADDITIVE);
              }else if(!type.compare("multiplicative"))
              {
                modes.push_back(FROSchParameterList::MULTIPLICATIVE);
              }else
              {
                comm.print("ERROR: unknown combine-lvl type '" + t + "'");
                return false;
              }
            }
            params.set_combine_lvl(modes);
          }else
          {
            comm.print("ERROR: not matching number of combine-lvl: nlevels:  " + std::to_string(params.get_nlevels()) + " number of given combine-lvl: " + std::to_string(it->second.size()) + "'");
            return false;
          }
        }
        return true;
      }

      inline bool parse_fpl_reuse_sf_ao(const Dist::Comm &comm, SimpleArgParser &args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-reuse-sf-ao");
        if(it != nullptr)
        {
          if(it->second.size() == 1u)
          {
            bool rsf_ao;
            std::string type = it->second.front();
            std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

            if(!type.compare("true"))
            {
              rsf_ao = true;
            }else if(!type.compare("false"))
            {
              rsf_ao = false;
            }else
            {
              comm.print("ERROR: unknown reuse-sf-ao '" + it->second.front() + "'");
              return false;
            }
            params.set_reuse_sf_ao(rsf_ao);
          }else if(int(it->second.size()) == params.get_nlevels())
          {
            std::vector<bool> rsf_ao;
            rsf_ao.reserve(params.get_nlevels());
            for(const auto& t : it->second)
            {
              std::string type = t;
              std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

              if(!type.compare("true"))
              {
                rsf_ao.push_back(true);
              }else if(!type.compare("false"))
              {
                rsf_ao.push_back(false);
              }else
              {
                comm.print("ERROR: unknown reuse-sf-ao type '" + t + "'");
                return false;
              }
            }
            params.set_reuse_sf_ao(rsf_ao);
          }else
          {
            comm.print("ERROR: not matching number of reuse-sf-ao: nlevels:  " + std::to_string(params.get_nlevels()) + " number of given reuse-sf-ao: " + std::to_string(it->second.size()) + "'");
            return false;
          }
        }
        return true;
      }

      inline bool parse_fpl_reuse_sf_cm(const Dist::Comm &comm, SimpleArgParser &args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-reuse-sf-cm");
        if(it != nullptr)
        {
          if(it->second.size() == 1u)
          {
            bool rsf_cm;
            std::string type = it->second.front();
            std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

            if(!type.compare("true"))
            {
              rsf_cm = true;
            }else if(!type.compare("false"))
            {
              rsf_cm = false;
            }else
            {
              comm.print("ERROR: unknown reuse-sf-cm '" + it->second.front() + "'");
              return false;
            }
            params.set_reuse_sf_cm(rsf_cm);
          }else if(int(it->second.size()) == params.get_nlevels())
          {
            std::vector<bool> rsf_cm;
            rsf_cm.reserve(params.get_nlevels());
            for(const auto& t : it->second)
            {
              std::string type = t;
              std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

              if(!type.compare("true"))
              {
                rsf_cm.push_back(true);
              }else if(!type.compare("false"))
              {
                rsf_cm.push_back(false);
              }else
              {
                comm.print("ERROR: unknown reuse-sf-cm type '" + t + "'");
                return false;
              }
            }
            params.set_reuse_sf_cm(rsf_cm);
          }else
          {
            comm.print("ERROR: not matching number of reuse-sf-cm: nlevels:  " + std::to_string(params.get_nlevels()) + " number of given reuse-sf-cm: " + std::to_string(it->second.size()) + "'");
            return false;
          }
        }
        return true;
      }

      inline bool parse_fpl_reuse_sf_ext(const Dist::Comm &comm, SimpleArgParser &args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-reuse-sf-ext");
        if(it != nullptr)
        {
          if(it->second.size() == 1u)
          {
            bool rsf_ext;
            std::string type = it->second.front();
            std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

            if(!type.compare("true"))
            {
              rsf_ext = true;
            }else if(!type.compare("false"))
            {
              rsf_ext = false;
            }else
            {
              comm.print("ERROR: unknown reuse-sf-ext '" + it->second.front() + "'");
              return false;
            }
            params.set_reuse_sf_ext(rsf_ext);
          }else if(int(it->second.size()) == params.get_nlevels())
          {
            std::vector<bool> rsf_ext;
            rsf_ext.reserve(params.get_nlevels());
            for(const auto& t : it->second)
            {
              std::string type = t;
              std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

              if(!type.compare("true"))
              {
                rsf_ext.push_back(true);
              }else if(!type.compare("false"))
              {
                rsf_ext.push_back(false);
              }else
              {
                comm.print("ERROR: unknown reuse-sf-cm type '" + t + "'");
                return false;
              }
            }
            params.set_reuse_sf_ext(rsf_ext);
          }else
          {
            comm.print("ERROR: not matching number of reuse-sf-ext: nlevels:  " + std::to_string(params.get_nlevels()) + " number of given reuse-sf-ext: " + std::to_string(it->second.size()) + "'");
            return false;
          }
        }
        return true;
      }

      inline bool parse_fpl_reuse_coarse_matrix(const Dist::Comm &comm, SimpleArgParser &args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-reuse-coarse-matrix");
        if(it != nullptr)
        {
          if(it->second.size() == 1u)
          {
            bool reuse_cm;
            std::string type = it->second.front();
            std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

            if(!type.compare("true"))
            {
              reuse_cm = true;
            }else if(!type.compare("false"))
            {
              reuse_cm = false;
            }else
            {
              comm.print("ERROR: unknown reuse-coarse-matrix '" + it->second.front() + "'");
              return false;
            }
            params.set_reuse_coarse_matrix(reuse_cm);
          }else if(int(it->second.size()) == params.get_nlevels())
          {
            std::vector<bool> reuse_cm;
            reuse_cm.reserve(params.get_nlevels());
            for(const auto& t : it->second)
            {
              std::string type = t;
              std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

              if(!type.compare("true"))
              {
                reuse_cm.push_back(true);
              }else if(!type.compare("false"))
              {
                reuse_cm.push_back(false);
              }else
              {
                comm.print("ERROR: unknown reuse-coarse-matrix type '" + t + "'");
                return false;
              }
            }
            params.set_reuse_coarse_matrix(reuse_cm);
          }else
          {
            comm.print("ERROR: not matching number of reuse-coarse-matrix: nlevels:  " + std::to_string(params.get_nlevels()) + " number of given reuse-coarse-matrix: " + std::to_string(it->second.size()) + "'");
            return false;
          }
        }
        return true;
      }

      inline bool parse_fpl_reuse_coarse_basis(const Dist::Comm &comm, SimpleArgParser &args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-reuse-coarse-basis");
        if(it != nullptr)
        {
          if(it->second.size() == 1u)
          {
            bool reuse_cb;
            std::string type = it->second.front();
            std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

            if(!type.compare("true"))
            {
              reuse_cb = true;
            }else if(!type.compare("false"))
            {
              reuse_cb = false;
            }else
            {
              comm.print("ERROR: unknown reuse-coarse-basis '" + it->second.front() + "'");
              return false;
            }
            params.set_reuse_coarse_basis(reuse_cb);
          }else if(int(it->second.size()) == params.get_nlevels())
          {
            std::vector<bool> reuse_cb;
            reuse_cb.reserve(params.get_nlevels());
            for(const auto& t : it->second)
            {
              std::string type = t;
              std::for_each(type.begin(), type.end(), [](char& c) { c = std::tolower(c); });

              if(!type.compare("true"))
              {
                reuse_cb.push_back(true);
              }else if(!type.compare("false"))
              {
                reuse_cb.push_back(false);
              }else
              {
                comm.print("ERROR: unknown reuse-coarse-basis type '" + t + "'");
                return false;
              }
            }
            params.set_reuse_coarse_basis(reuse_cb);
          }else
          {
            comm.print("ERROR: not matching number of reuse-coarse-basis: nlevels:  " + std::to_string(params.get_nlevels()) + " number of given reuse-coarse-basis: " + std::to_string(it->second.size()) + "'");
            return false;
          }
        }
        return true;
      }

      inline bool parse_fpl_phi_dt(const Dist::Comm &comm, SimpleArgParser &args, FROSchParameterList &params)
      {
        auto it = args.query("frosch-phi-dropping-threshold");
        if(it != nullptr)
        {
          if(it->second.size() == 1u)
          {
            double phi_dt = std::stod(it->second.front());
            if (phi_dt < 0.)
            {
              /* __TODO:__ Stimmt das? */
              comm.print("ERROR: unknown dropping threshold smaller than 0.:  '" + it->second.front() + "'");
              return false;
            }
            params.set_phi_dropping_threshold(phi_dt);
          }else if(int(it->second.size()) == params.get_nlevels())
          {
            std::vector<double> phi_dt;
            phi_dt.reserve(params.get_nlevels());
            for(const auto& t : it->second)
            {
              double dt = std::stod(t);
              if (dt < 0.)
              {
                /* __TODO:__ Stimmt das? */
                comm.print("ERROR: unknown partition imblance smaller than 0.:  '" + t + "'");
                return false;
              }
              phi_dt.push_back(dt);
            }
            params.set_phi_dropping_threshold(phi_dt);
          }else
          {
            comm.print("ERROR: not matching number of phi-dropping-threshold: nlevels:  " + std::to_string(params.get_nlevels()) + " number of given phi-dropping-threshold: " + std::to_string(it->second.size()) + "'");
            return false;
          }
        }
        return true;
      }

      bool FROSchParameterList::parse_args(SimpleArgParser& args)
      {
        {
          auto it = args.query("frosch-plist");
          if(it != nullptr)
          {
            if(it->second.size() == 1u)
            {
              read_from_xml_file(it->second.front());
              /* given parameter list wins */
              return true;
            }else
            {
              _comm.print("ERROR: multiple XML files are given");
              return false;
            }
          }
        }

        if(!parse_fpl_print(_comm, args, *this))
        {
          return false;
        }
        if(!parse_fpl_timer(_comm, args, *this))
        {
          return false;
        }
        if(!parse_fpl_levels(_comm, args, *this))
        {
          return false;
        }
        if(!parse_fpl_dim(_comm, args, *this))
        {
          return false;
        }
        if(!parse_fpl_overlaps(_comm, args, *this))
        {
          return false;
        }
        if(!parse_fpl_subreg(_comm, args, *this))
        {
          return false;
        }
        if(!parse_fpl_gsteps(_comm, args, *this))
        {
          return false;
        }
        if(!parse_fpl_parti_type(_comm, args, *this))
        {
          return false;
        }
        if(!parse_fpl_parti_approach(_comm, args, *this))
        {
          return false;
        }
        if(!parse_fpl_parti_imbl(_comm, args, *this))
        {
          return false;
        }
        if(!parse_fpl_ipou_pres(_comm, args, *this))
        {
          return false;
        }
        if(!parse_fpl_ipou_velo(_comm, args, *this))
        {
          return false;
        }
        if(!parse_fpl_solver_ao(_comm, args, *this))
        {
          return false;
        }
        if(!parse_fpl_solver_ext(_comm, args, *this))
        {
          return false;
        }
        if(!parse_fpl_solver_coarse(_comm, args, *this))
        {
          return false;
        }
        if(!parse_fpl_precond(_comm, args, *this))
        {
          return false;
        }
        if(!parse_fpl_cspace_pres(_comm, args, *this))
        {
          return false;
        }
        if(!parse_fpl_cspace_velo(_comm, args, *this))
        {
          return false;
        }
        if(!parse_fpl_exclude_velo_pres(_comm, args, *this))
        {
          return false;
        }
        if(!parse_fpl_exclude_pres_velo(_comm, args, *this))
        {
          return false;
        }
        if(!parse_fpl_combine(_comm, args, *this))
        {
          return false;
        }
        if(!parse_fpl_combine_lvl(_comm, args, *this))
        {
          return false;
        }
        if (!parse_fpl_reuse_sf_ao(_comm, args, *this))
        {
          return false;
        }
        if (!parse_fpl_reuse_sf_cm(_comm, args, *this))
        {
          return false;
        }
        if (!parse_fpl_reuse_sf_ext(_comm, args, *this))
        {
          return false;
        }
        if (!parse_fpl_reuse_coarse_matrix(_comm, args, *this))
        {
          return false;
        }
        if (!parse_fpl_reuse_coarse_basis(_comm, args, *this))
        {
          return false;
        }
        if (!parse_fpl_phi_dt(_comm, args, *this))
        {
          return false;
        }

        return true;
      }

      /**
       * \brief Tpetra wrapper core data class
       */
      class TpetraCore
      {
        public:
          // some typedefs
          using UN = unsigned;
          using SC = double;
          using LO = int;
          using GO = FROSch::DefaultGlobalOrdinal;
          using NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType;

          using AGS   = Teuchos::Array<GO>::size_type;
          using ARCPS = Teuchos::ArrayRCP<std::size_t>::size_type;
          using TGS   = Tpetra::global_size_t;
          using VGS   = std::vector<GO>::size_type;

          //----------------------------------------
          Teuchos::RCP<const Teuchos::Comm<int>> _comm;

          const LO _num_owned_dofs;
          const GO _my_dof_offset;
          const int _dpe_pres;
          /// Tpetra number of non-zero entries per row array
          Teuchos::ArrayRCP<TGS> _num_nze;
          /// Tpetra column indices array
          Teuchos::ArrayRCP<GO> _col_idx;

          /// row and column maps are the same
          Teuchos::RCP<const Tpetra::Map<LO, GO, NO>> _map;

          /// handles right hand side and solution
          Teuchos::RCP<Tpetra::MultiVector<SC, LO,GO,NO>> _vec_def, _vec_cor;
          Teuchos::RCP<Tpetra::CrsMatrix<SC,LO,GO,NO>> _matrix;

          /// parameter list
          Teuchos::RCP<Teuchos::ParameterList> _params;

          bool _use_timer;
          Teuchos::RCP<Teuchos::StackedTimer> _stacked_timer;

          //template<typename IT_>
          explicit TpetraCore(const void* comm, IndexType my_dof_offset, IndexType num_owned_dofs, const Trilinos::FROSchParameterList & params) :
            _comm(Teuchos::rcp(new Teuchos::MpiComm<int>(*reinterpret_cast<const MPI_Comm*>(comm)))),
            _num_owned_dofs(static_cast<LO>(num_owned_dofs)),
            _my_dof_offset(static_cast<GO>(my_dof_offset)),
            _dpe_pres(params.get_dpe_pres()),
            _num_nze(static_cast<ARCPS>(_num_owned_dofs), static_cast<TGS>(0)),
            _col_idx(static_cast<ARCPS>(_num_owned_dofs), static_cast<GO>(0)),
            _map(nullptr),
            _vec_def(nullptr),
            _vec_cor(nullptr),
            _matrix(nullptr),
            _params(Teuchos::null),
            _use_timer(params.get_use_timer()),
            _stacked_timer(Teuchos::rcp(new Teuchos::StackedTimer("TpetraCore")))
          {
            void* pcore = params.get_core();
            if(pcore != nullptr)
            {
              _params = Teuchos::rcp(reinterpret_cast<Teuchos::ParameterList*>(pcore), false);
            }

            if(_use_timer)
            {
              _comm->barrier();
              Teuchos::TimeMonitor::setStackedTimer(_stacked_timer);
            }
          }

          virtual ~TpetraCore()
          {
            if(_use_timer)
            {
              _comm->barrier();
              _stacked_timer->stop("TpetraCore");
              Teuchos::StackedTimer::OutputOptions options;
              options.output_fraction = options.output_histogram = options.output_minmax = true;
              Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
              _stacked_timer->report(*out, _comm, options);

            }
          }

          template<typename IT_>
          void init_core(const IT_* row_ptr, const IT_* col_idx)
          {
            if(_use_timer)
            {
              _stacked_timer->start("TpetraCore::init_core");
              _stacked_timer->start("TpetraCore::init_core::build_map");
            }
            /// create Tpetra index vectors
            _num_nze = Teuchos::ArrayRCP<TGS>(static_cast<ARCPS>(_num_owned_dofs), static_cast<TGS>(0));

            /// loop over all owned matrix rows
            for(LO i(0); i < _num_owned_dofs; ++i)
              _num_nze[i] = static_cast<TGS>(row_ptr[i+1] - row_ptr[i]);

            /// create Tpetra map
            /**
             *  The default constructor of a contiguous distribution map is the same as
             *  FEAT layout:
             *
             * Nevertheless, we use an slightly more complex constructor, just be sure
             **/
            std::vector<GO> element_list_vec(static_cast<VGS>(_num_owned_dofs));
            for(IndexType i = 0; i < static_cast<IndexType>(_num_owned_dofs); ++i)
              element_list_vec.at(i) = static_cast<GO>(_my_dof_offset + i);

            GO num_total_dofs = 0;
            GO my_dof_offset_tmp = static_cast<GO>(_num_owned_dofs);
            Teuchos::reduceAll(*_comm, Teuchos::REDUCE_SUM, 1, &my_dof_offset_tmp, &num_total_dofs);

            const GO _index_base = 0;
            const Teuchos::ArrayView<const GO> element_list = Teuchos::arrayViewFromVector(element_list_vec);
            _map = Teuchos::rcp(new Tpetra::Map(static_cast<TGS>(num_total_dofs), element_list, _index_base, _comm));

            if(_use_timer)
            {
              _stacked_timer->stop("TpetraCore::init_core::build_map");
              _stacked_timer->start("TpetraCore::init_core::build_system_matrix");
            }

            /// create matrix
            auto crsgraph = Teuchos::rcp(new Tpetra::CrsGraph<LO, GO, NO >(_map, _num_nze()));
            const TGS max_nnz_per_row = *std::max_element(_num_nze.begin(), _num_nze.end());
            std::vector<GO> col_buffer(max_nnz_per_row);  // Reusable buffer

            /// initialize crs graph matrix
            /**
             *  Good practice would be to initialize a crs graph, but this seems not to
             *  be very practically for the Xpetra interfaces to maps and crs graphs
             **/
            AGS entries_start = 0;
            for(GO i = 0; i < static_cast<GO>(_num_owned_dofs); ++i)
            {
              const auto num_entries = static_cast<AGS>(_num_nze[i]);

              for(Teuchos::Array<GO>::size_type j = 0; j < num_entries; ++j)
                col_buffer[j] = static_cast<GO>(col_idx[entries_start + j]);

              crsgraph->insertGlobalIndices(static_cast<GO>(_my_dof_offset + i), num_entries, col_buffer.data());
              entries_start += num_entries;
            }
            crsgraph->fillComplete();
            _matrix = Teuchos::rcp(new Tpetra::CrsMatrix<SC, LO, GO, NO>(crsgraph));
            if(_use_timer)
            {
              _stacked_timer->stop("TpetraCore::init_core::build_system_matrix");
            }

            // create vectors
            const LO num_vectors = 1;
            _vec_def = Teuchos::rcp(new Tpetra::MultiVector<SC, LO, GO, NO>(_map, num_vectors));
            _vec_cor = Teuchos::rcp(new Tpetra::MultiVector<SC, LO, GO, NO>(_map, num_vectors));

            GO _col_idx_size = static_cast<GO>(row_ptr[_num_owned_dofs]);
            _col_idx.resize(static_cast<AGS>(_col_idx_size));
            for(GO i = 0; i < _col_idx_size; ++i)
              _col_idx[(static_cast<ARCPS>(i))] = static_cast<GO>(col_idx[i]);

            // set entries for TwoLevelBlockPreconditionier
            if (! sublist(sublist(_params, "Preconditioner Types"), "FROSch")->get("FROSch Preconditioner Type","TwoLevelPreconditioner").compare("TwoLevelBlockPreconditioner"))
            {
              const int NumberOfBlocks = 1;
              Teuchos::ArrayRCP<Teuchos::RCP<Xpetra::Map<LO, GO, NO>>> XRepeatedMaps(NumberOfBlocks);
              auto xcrsgraph = Teuchos::rcp(new const Xpetra::TpetraCrsGraph<LO, GO, NO>(Teuchos::rcp_const_cast<Tpetra::CrsGraph<LO, GO, NO >>(_matrix->getCrsGraph ())));
              XRepeatedMaps[0] = FROSch::BuildRepeatedMapNonConst<LO, GO, NO>(xcrsgraph);
              sublist(sublist(_params, "Preconditioner Types"), "FROSch")->set("Repeated Map Vector", XRepeatedMaps);

              Teuchos::ArrayRCP<FROSch::DofOrdering> dofOrderings(NumberOfBlocks);
              dofOrderings[0] = FROSch::NodeWise;
              sublist(sublist(_params, "Preconditioner Types"), "FROSch")->set("DofOrdering Vector", dofOrderings);

              Teuchos::ArrayRCP<unsigned> dofsPerNodeVector(NumberOfBlocks);
              dofsPerNodeVector[0] = _dpe_pres;
              sublist(sublist(_params, "Preconditioner Types"), "FROSch")->set("DofsPerNode Vector", dofsPerNodeVector);

              /* NULL SPACE */
              Teuchos::ArrayRCP<Teuchos::RCP<Xpetra::MultiVector<SC,LO,GO,NO>>> nullSpaceBasis(NumberOfBlocks);
              nullSpaceBasis[0] = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(XRepeatedMaps[0].getConst(), 1);
              unsigned i = 0;
              for (unsigned j = 0; j < XRepeatedMaps[i]->getLocalNumElements(); j += dofsPerNodeVector[i]) {
                nullSpaceBasis[0]->getDataNonConst(i)[XRepeatedMaps[i]->getLocalElement(XRepeatedMaps[i]->getGlobalElement(j))] = Teuchos::ScalarTraits<SC>::one();
              }
              sublist(sublist(_params, "Preconditioner Types"), "FROSch")->set("Null Space Vector", nullSpaceBasis);
            }else if(! sublist(sublist(_params, "Preconditioner Types"), "FROSch")->get("FROSch Preconditioner Type","TwoLevelPreconditioner").compare("TwoLevelPreconditioner"))
            {
              /* NULL SPACE */
              auto xcrsgraph = Teuchos::rcp(new const Xpetra::TpetraCrsGraph<LO, GO, NO>(Teuchos::rcp_const_cast<Tpetra::CrsGraph<LO, GO, NO >>(_matrix->getCrsGraph ())));
              Teuchos::RCP<Xpetra::Map<LO, GO, NO>> XRepeatedMap = FROSch::BuildRepeatedMapNonConst<LO, GO, NO>(xcrsgraph);
              sublist(sublist(_params, "Preconditioner Types"), "FROSch")->set("Repeated Map", XRepeatedMap);

              Teuchos::RCP<Xpetra::MultiVector<SC,LO,GO,NO>> nullSpaceBasis = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(XRepeatedMap.getConst(), 1);
              for (unsigned j = 0; j < XRepeatedMap->getLocalNumElements(); j += _dpe_pres) {
                nullSpaceBasis->getDataNonConst(0)[XRepeatedMap->getLocalElement(XRepeatedMap->getGlobalElement(j))] = Teuchos::ScalarTraits<SC>::one();
              }
              sublist(sublist(_params, "Preconditioner Types"), "FROSch")->set("Null Space", nullSpaceBasis);
            }

            if(_use_timer)
            {
              _stacked_timer->stop("TpetraCore::init_core");
            }
          }

          virtual void set_matrix_values(const double* vals)
          {
            if(this->_use_timer)
            {
              this->_stacked_timer->start("TpetraCore::SetValues");
            }

            Teuchos::ArrayView<const LO> Indices;

            _matrix->resumeFill();
            AGS entries_start = 0;

            for(LO i(0); i < _num_owned_dofs; ++i)
            {
              const auto num_entries = static_cast<Teuchos::Array<SC>::size_type>(_num_nze[i]);
              // copy vals to Teuchos::Array
              Teuchos::Array<SC> vals_ar(num_entries);
              for(Teuchos::Array<SC>::size_type j = 0; j < num_entries; ++j)
                vals_ar[j] = static_cast<SC>(vals[static_cast<VGS>(entries_start + j)]);

              Teuchos::Array<GO> cols(num_entries);
              for(Teuchos::Array<GO>::size_type j = 0; j < num_entries; ++j)
                cols[j] = static_cast<GO>(_col_idx[entries_start + j]);
              _matrix->replaceGlobalValues(_my_dof_offset + static_cast<GO>(i),
                                           cols.view(0, static_cast<AGS>(num_entries)),
                                           vals_ar.view(0, static_cast<Teuchos::Array<SC>::size_type>(num_entries)));

              entries_start += num_entries;
            }

            _matrix->fillComplete();

            if(_use_timer)
            {
              _stacked_timer->stop("TpetraCore::SetValues");
            }
          }

          void set_vec_cor_values(const double* vals)
          {
            if(vals != nullptr)
              for(LO i(0); i < _num_owned_dofs; ++i)
                _vec_cor->replaceLocalValue(i, 0, static_cast<SC>(vals[i]));
            else
              _vec_cor->putScalar(static_cast<SC>(0.));
          }

          void set_vec_def_values(const double* vals)
          {
            if(vals != nullptr)
              for(LO i(0); i < _num_owned_dofs; ++i)
                _vec_def->replaceLocalValue(i, 0, static_cast<SC>(vals[i]));
            else
              _vec_def->putScalar(static_cast<SC>(0.));
          }

          void get_vec_cor_values(double* vals) const
          {
            auto data_vec = _vec_cor->getData(0);
            for(LO i(0); i < _num_owned_dofs; ++i)
              vals[i] = data_vec[i];
          }

          void get_vec_def_values(double* vals) const
          {
            auto data_vec = _vec_def->getData(0);
            for(LO i(0); i < _num_owned_dofs; ++i)
              vals[i] = data_vec[i];
          }
      }; // class TpetraCore

      void set_parameter_list(void* core, const FROSchParameterList& params)
      {
        void* pcore = params.get_core();
        if(pcore == nullptr)
        {
          reinterpret_cast<TpetraCore*>(core)->_params = Teuchos::null;
        }else
        {
          reinterpret_cast<TpetraCore*>(core)->_params
            = Teuchos::rcp(reinterpret_cast<Teuchos::ParameterList*>(pcore), false);
        }
      }

      void* create_core(const void* comm, IndexType dof_offset, IndexType num_owned_dofs,
                        const unsigned int* row_ptr, const unsigned int* col_idx, const FROSchParameterList& params)
      {
        TpetraCore *core = new TpetraCore(comm, dof_offset, num_owned_dofs, params);
        core->init_core(row_ptr, col_idx);
        return (void*)core;
      }

      void* create_core(const void* comm, IndexType dof_offset, IndexType num_owned_dofs,
                        const unsigned long* row_ptr, const unsigned long* col_idx, const FROSchParameterList& params)
      {
        TpetraCore *core = new TpetraCore(comm, dof_offset, num_owned_dofs, params);
        core->init_core(row_ptr, col_idx);
        return (void*)core;
      }

      void* create_core(const void* comm, IndexType dof_offset, IndexType num_owned_dofs,
                        const unsigned long long* row_ptr, const unsigned long long* col_idx, const FROSchParameterList& params)
      {
        TpetraCore *core = new TpetraCore(comm, dof_offset, num_owned_dofs, params);
        core->init_core(row_ptr, col_idx);
        return (void*)core;
      }

      void destroy_core(void* core)
      {
        delete reinterpret_cast<TpetraCore*>(core);
      }

      void set_matrix_values(void* core, const double* vals)
      {
        reinterpret_cast<TpetraCore*>(core)->set_matrix_values(vals);
      }

      void set_vec_cor_values(void* core, const double* vals)
      {
        reinterpret_cast<TpetraCore*>(core)->set_vec_cor_values(vals);
      }

      void set_vec_def_values(void* core, const double* vals)
      {
        reinterpret_cast<TpetraCore*>(core)->set_vec_def_values(vals);
      }

      void get_vec_cor_values(const void* core, double* vals)
      {
        reinterpret_cast<const TpetraCore*>(core)->get_vec_cor_values(vals);
      }

      void get_vec_def_values(const void* core, double* vals)
      {
        reinterpret_cast<const TpetraCore*>(core)->get_vec_def_values(vals);
      }

      /**
       * \brief Tpetra wrapper core Stokes data class
       *
       * __ATTENTION:__ Only for discontinuous pressure
       */
      class TpetraCoreStokes : public TpetraCore
      {
        public:
          // some typedefs
          typedef TpetraCore BaseClass;
          using UN = unsigned;
          using SC = double;
          using LO = int;
          using GO = FROSch::DefaultGlobalOrdinal;
          using NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType;

          using AGS   = Teuchos::Array<GO>::size_type;
          using ARCPS = Teuchos::ArrayRCP<std::size_t>::size_type;
          using TGS   = Tpetra::global_size_t;
          using VGS   = std::vector<GO>::size_type;

          //----------------------------------------
          const int _dpe_velo;
          // number of velocity/pressure dofs owned by this process
          const IndexType _num_owned_velo_dofs;
          const IndexType _num_owned_pres_dofs; // <- given to constructor
          // index of first velocity/pressure dof owned by this process
          const IndexType _first_owned_velo_dof;
          const IndexType _first_owned_pres_dof;
          // total number of velocity/pressure dofs
          const IndexType _num_global_velo_dofs;
          const IndexType _num_global_pres_dofs;

          Teuchos::ArrayRCP<GO> _all_global_dof_offset;
          Teuchos::ArrayRCP<IndexType> _all_num_owned_velo_dofs, _all_num_owned_pres_dofs;
          Teuchos::ArrayRCP<IndexType> _all_first_owned_velo_dof, _all_first_owned_pres_dof;

          explicit TpetraCoreStokes(const void* comm,
                                    IndexType my_dof_offset, IndexType num_owned_dofs,
                                    IndexType num_owned_velo_dofs, IndexType num_owned_pres_dofs,
                                    IndexType first_owned_velo_dof, IndexType first_owned_pres_dof,
                                    IndexType num_global_velo_dofs, IndexType num_global_pres_dofs,
                                    const Trilinos::FROSchParameterList & params) :
            BaseClass(comm, my_dof_offset, num_owned_dofs, params),
            _dpe_velo(params.get_dpe_velo()),
            _num_owned_velo_dofs(num_owned_velo_dofs),
            _num_owned_pres_dofs(num_owned_pres_dofs),
            _first_owned_velo_dof(first_owned_velo_dof),
            _first_owned_pres_dof(first_owned_pres_dof),
            _num_global_velo_dofs(num_global_velo_dofs),
            _num_global_pres_dofs(num_global_pres_dofs),
            _all_global_dof_offset(static_cast<ARCPS>(this->_comm->getSize())),
            _all_num_owned_velo_dofs(static_cast<ARCPS>(this->_comm->getSize())),
            _all_num_owned_pres_dofs(static_cast<ARCPS>(this->_comm->getSize())),
            _all_first_owned_velo_dof(static_cast<ARCPS>(this->_comm->getSize())),
            _all_first_owned_pres_dof(static_cast<ARCPS>(this->_comm->getSize()))
          {
          }

          template<typename IT_>
          void init_core(const IT_* row_ptr, const IT_* col_idx)
          {
            if(_use_timer)
            {
              _comm->barrier();
              Teuchos::TimeMonitor::setStackedTimer(this->_stacked_timer);
              _stacked_timer->start("TpetraCoreStokes::init_core");
            }

            if(_use_timer)
            {
              _stacked_timer->start("TpetraCoreStokes::init_core::gather");
            }
            /// gather global information
            Teuchos::gatherAll(*(_comm), static_cast<int>(1), &_my_dof_offset, static_cast<int>(1*_comm->getSize()), _all_global_dof_offset.getRawPtr());
            Teuchos::gatherAll(*(_comm), static_cast<int>(1), &_num_owned_velo_dofs, static_cast<int>(1*_comm->getSize()), _all_num_owned_velo_dofs.getRawPtr());
            Teuchos::gatherAll(*(_comm), static_cast<int>(1), &_num_owned_pres_dofs, static_cast<int>(1*_comm->getSize()), _all_num_owned_pres_dofs.getRawPtr());
            Teuchos::gatherAll(*(_comm), static_cast<int>(1), &_first_owned_velo_dof, static_cast<int>(1*_comm->getSize()), _all_first_owned_velo_dof.getRawPtr());
            Teuchos::gatherAll(*(_comm), static_cast<int>(1), &_first_owned_pres_dof, static_cast<int>(1*_comm->getSize()), _all_first_owned_pres_dof.getRawPtr());
            if(_use_timer)
            {
              _stacked_timer->stop("TpetraCoreStokes::init_core::gather");
              _stacked_timer->start("TpetraCoreStokes::init_core::build_map");
            }

            /// create Tpetra index vectors
            _num_nze = Teuchos::ArrayRCP<TGS>(static_cast<ARCPS>(_num_owned_dofs), static_cast<TGS>(0));

            /// loop over all owned matrix rows
            for(LO i(0); i < _num_owned_dofs; ++i)
              _num_nze[i] = static_cast<TGS>(row_ptr[i+1] - row_ptr[i]);

            /// construct map.  we need to change the global order:
            /// The order of the scalarized matrix is, e.g.:
            ///
            ///       170  48  138  48  144  48
            /// 170 : A11 B11  A12 B12  A13 B13
            ///  48 : D11  0   D12 0    D13 0
            /// -------------------------
            /// 138 : A21 B21  A22 B22  A23 B23
            ///  48 : D21  0   D22  0   D23  0
            /// -------------------------
            /// 144 : A31 B31  A32 B32  A33 B33
            ///  48 : D31  0   D32  0   D33  0
            ///
            const GO _index_base = 0;
            std::vector<GO> element_list_vec(static_cast<VGS>(_num_owned_dofs));
            for(IndexType i = 0; i < static_cast<IndexType>(_num_owned_velo_dofs); ++i)
              element_list_vec.at(i) = static_cast<GO>(_first_owned_velo_dof + i);
            for(IndexType i = 0; i < static_cast<IndexType>(_num_owned_pres_dofs); ++i)
              element_list_vec.at(_num_owned_velo_dofs + i) = static_cast<GO>(_num_global_velo_dofs + _first_owned_pres_dof + i);
            const Teuchos::ArrayView<const GO> element_list = Teuchos::arrayViewFromVector(element_list_vec);

            _map = Teuchos::rcp(new Tpetra::Map(static_cast<TGS>(_num_global_velo_dofs + _num_global_pres_dofs), element_list, _index_base, _comm));

            if(_use_timer)
            {
              _stacked_timer->stop("TpetraCoreStokes::init_core::build_map");
              _stacked_timer->start("TpetraCoreStokes::init_core::build_repeated_maps");
            }
            const int NumberOfBlocks = 2;

            // Pre-compute rank boundaries for O(1) access
            std::vector<GO> rank_start(_comm->getSize());
            std::vector<GO> rank_velo_end(_comm->getSize());
            std::vector<GO> rank_total_end(_comm->getSize());
            for (int rank = 0; rank < _comm->getSize(); ++rank)
            {
              rank_start[rank] = _all_global_dof_offset[rank];
              rank_velo_end[rank] = _all_global_dof_offset[rank] + _all_num_owned_velo_dofs[rank];
              rank_total_end[rank] = _all_global_dof_offset[rank] + _all_num_owned_velo_dofs[rank] + _all_num_owned_pres_dofs[rank];
            }

            /// build repeated maps and set
            {
              if(_use_timer)
              {
                _stacked_timer->start("TpetraCoreStokes::init_core::build_repeated_maps::velocity_graph");
              }

              std::vector<GO> element_list_vec_velo(static_cast<VGS>(_num_owned_velo_dofs));
              for(IndexType i = 0; i < static_cast<IndexType>(_num_owned_velo_dofs); ++i)
                element_list_vec_velo.at(i) = static_cast<GO>(_first_owned_velo_dof + i);
              const Teuchos::ArrayView<const GO> element_list_velo = Teuchos::arrayViewFromVector(element_list_vec_velo);
              Teuchos::RCP<const Tpetra::Map<LO, GO, NO>> _map_velo = Teuchos::rcp(new Tpetra::Map(static_cast<TGS>(_num_global_velo_dofs), element_list_velo, _index_base, _comm));

              auto crsgraph = Teuchos::rcp(new Tpetra::CrsGraph<LO, GO, NO >(_map_velo, _num_nze.view(0, static_cast<AGS>(_num_owned_velo_dofs))));
              const TGS max_nnz_per_row = *std::max_element(_num_nze.begin(), _num_nze.end());
              std::vector<GO> col_buffer(max_nnz_per_row);  // Reusable buffer

              /// insert entries for velocity
              for(IndexType row(0); row < _num_owned_velo_dofs; ++row)
              {
                Teuchos::Array<GO> cols(static_cast<GO>(_num_nze[static_cast<ARCPS>(row)]));
                GO counter = 0;

                // Start scanning from rank 0 for each row
                int current_rank = 0;

                for(auto colidx = row_ptr[row]; colidx < row_ptr[row+1]; ++colidx)
                {
                  auto col = col_idx[colidx];

                  // Advance rank until we find one that might contain this column
                  while(current_rank < int(_comm->getSize()) && GO(col) >= rank_total_end[current_rank])
                  {
                    current_rank++;
                  }

                  // Check if current rank contains this column
                  if(current_rank < int(_comm->getSize()) && GO(col) >= rank_start[current_rank])
                  {
                    if(GO(col) < rank_velo_end[current_rank])
                    {
                      // Velocity DOF
                      GO mapped_dof = _all_first_owned_velo_dof[current_rank] + (col - rank_start[current_rank]);
                      cols[counter++] = mapped_dof;
                    }
                  }
                  // Note: current_rank is NOT reset - it continues from where it left off
                }

                crsgraph->insertGlobalIndices(static_cast<GO>(_first_owned_velo_dof + row),
                                              cols.view(0, static_cast<AGS>(counter)));
              }
              crsgraph->fillComplete();

              if(_use_timer)
              {
                _stacked_timer->stop("TpetraCoreStokes::init_core::build_repeated_maps::velocity_graph");
              }

              Teuchos::ArrayRCP<Teuchos::RCP<Tpetra::Map<LO, GO, NO>>> RepeatedMaps(NumberOfBlocks);
              Teuchos::ArrayRCP<Teuchos::RCP<Xpetra::Map<LO, GO, NO>>> XRepeatedMaps(NumberOfBlocks);

              const GO _index_base_tmp = 0;
              /* pressure map */
              std::vector<GO> element_list_vec_tmp(static_cast<VGS>(_num_owned_pres_dofs));
              for(IndexType i = 0; i < static_cast<IndexType>(_num_owned_pres_dofs); ++i)
                element_list_vec_tmp.at(i) = static_cast<GO>(_first_owned_pres_dof + i);
              const Teuchos::ArrayView<const GO> element_list_tmp = Teuchos::arrayViewFromVector(element_list_vec_tmp);
              RepeatedMaps[1] = Teuchos::rcp(new Tpetra::Map(static_cast<TGS>(_num_owned_pres_dofs), element_list_tmp, _index_base_tmp, _comm));

              auto xcrsgraph = Teuchos::rcp(new Xpetra::TpetraCrsGraph<LO, GO, NO>(crsgraph));
              XRepeatedMaps[0] = FROSch::BuildRepeatedMapNonConst<LO, GO, NO>(xcrsgraph);
              // __ATTENTION:__ This works only for a discontinuous pressure
              XRepeatedMaps[1] = Teuchos::rcp(new Xpetra::TpetraMap<LO, GO, NO>(RepeatedMaps[1]));

              sublist(sublist(_params, "Preconditioner Types"), "FROSch")->set("Repeated Map Vector", XRepeatedMaps);

              /// set dof ordering and dofs per node
              Teuchos::ArrayRCP<FROSch::DofOrdering> dofOrderings(NumberOfBlocks);
              dofOrderings[0] = FROSch::NodeWise;
              dofOrderings[1] = FROSch::NodeWise;
              sublist(sublist(_params, "Preconditioner Types"), "FROSch")->set("DofOrdering Vector", dofOrderings);

              Teuchos::ArrayRCP<unsigned> dofsPerNodeVector(NumberOfBlocks);
              dofsPerNodeVector[0] = _dpe_velo;
              dofsPerNodeVector[1] = _dpe_pres;
              sublist(sublist(_params, "Preconditioner Types"), "FROSch")->set("DofsPerNode Vector", dofsPerNodeVector);

              /* Null Space */
              Teuchos::ArrayRCP<Teuchos::RCP<Xpetra::MultiVector<SC,LO,GO,NO>>> nullSpaceBasis(NumberOfBlocks);
              Teuchos::ArrayRCP<Teuchos::RCP<const Xpetra::Map<LO, GO, NO>>> XRepeatedMapsTmp(NumberOfBlocks);
              XRepeatedMapsTmp[0] = XRepeatedMaps[0].getConst();
              XRepeatedMapsTmp[1] = XRepeatedMaps[1].getConst();
              Teuchos::RCP<const Xpetra::Map<LO, GO, NO>> repeatedMap = FROSch::MergeMaps(XRepeatedMapsTmp);

              nullSpaceBasis[0] = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(repeatedMap, dofsPerNodeVector[0]);
              nullSpaceBasis[1] = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(repeatedMap, 1);
              for(unsigned i = 0; i < dofsPerNodeVector[0]; i++)
                for (unsigned j = i; j < XRepeatedMaps[0]->getLocalNumElements(); j += dofsPerNodeVector[0]) {
                  nullSpaceBasis[0]->getDataNonConst(i)[XRepeatedMaps[0]->getLocalElement(XRepeatedMaps[0]->getGlobalElement(j))] = Teuchos::ScalarTraits<SC>::one();
                }
              LO offset = XRepeatedMaps[0]->getLocalNumElements();
              for (unsigned j = 0; j < XRepeatedMaps[1]->getLocalNumElements(); j += dofsPerNodeVector[1]) {
                nullSpaceBasis[1]->getDataNonConst(0)[offset + XRepeatedMaps[1]->getLocalElement(XRepeatedMaps[1]->getGlobalElement(j))] = Teuchos::ScalarTraits<SC>::one();
              }

              sublist(sublist(_params, "Preconditioner Types"), "FROSch")->set("Null Space Vector", nullSpaceBasis);
            }
            if(_use_timer)
            {
              _stacked_timer->stop("TpetraCoreStokes::init_core::build_repeated_maps");
              _stacked_timer->start("TpetraCoreStokes::init_core::build_crsgraph_system");
            }

            /// _col_idx is for the method set_matrix_values
            GO _col_idx_size = static_cast<GO>(row_ptr[_num_owned_dofs]);
            _col_idx.resize(static_cast<AGS>(_col_idx_size));
            GO counter_col_idx = 0;

            /// build crsgraph
            auto crsgraph = Teuchos::rcp(new Tpetra::CrsGraph<LO, GO, NO >(_map, _num_nze.view(0, static_cast<AGS>(_num_owned_dofs))));

            /// insert entries for velocity
            for(IndexType row(0); row < _num_owned_velo_dofs; ++row)
            {
              Teuchos::Array<GO> cols(static_cast<GO>(_num_nze[static_cast<ARCPS>(row)]));
              GO counter = 0;

              // Start scanning from rank 0 for each row
              int current_rank = 0;

              for(auto colidx = row_ptr[row]; colidx < row_ptr[row+1]; ++colidx)
              {
                auto col = col_idx[colidx];

                // Advance rank until we find one that might contain this column
                while(current_rank < int(_comm->getSize()) && GO(col) >= rank_total_end[current_rank])
                {
                  current_rank++;
                }

                // Check if current rank contains this column
                if(current_rank < int(_comm->getSize()) && GO(col) >= rank_start[current_rank])
                {
                  if(GO(col) < rank_velo_end[current_rank])
                  {
                    // Velocity DOF
                    GO mapped_dof = _all_first_owned_velo_dof[current_rank] + (col - rank_start[current_rank]);
                    cols[counter++] = mapped_dof;
                    _col_idx[counter_col_idx++] = mapped_dof;
                  } else if (GO(col) < rank_total_end[current_rank]) {
                    // Pressure DOF
                    GO pres_offset = col - rank_velo_end[current_rank];
                    GO mapped_dof = _num_global_velo_dofs + _all_first_owned_pres_dof[current_rank] + pres_offset;
                    cols[counter++] = mapped_dof;
                    _col_idx[counter_col_idx++] = mapped_dof;
                  }
                }
                // Note: current_rank is NOT reset - it continues from where it left off
              }

              crsgraph->insertGlobalIndices(static_cast<GO>(_first_owned_velo_dof + row), cols.view(0, static_cast<AGS>(counter)));
            }
            /// insert entries for pressure
            for(IndexType row(_num_owned_velo_dofs); row < (_num_owned_velo_dofs + _num_owned_pres_dofs); ++row)
            {
              Teuchos::Array<GO> cols(static_cast<GO>(_num_nze[static_cast<ARCPS>(row)]));
              GO counter = 0;

              // Start scanning from rank 0 for each row
              int current_rank = 0;

              for(auto colidx = row_ptr[row]; colidx < row_ptr[row+1]; ++colidx)
              {
                auto col = col_idx[colidx];

                // Advance rank until we find one that might contain this column
                while(current_rank < int(_comm->getSize()) && GO(col) >= rank_total_end[current_rank])
                {
                  current_rank++;
                }

                // Check if current rank contains this column
                if(current_rank < int(_comm->getSize()) && GO(col) >= rank_start[current_rank])
                {

                  if(GO(col) < rank_velo_end[current_rank])
                  {
                    // Velocity DOF
                    GO mapped_dof = _all_first_owned_velo_dof[current_rank] + (col - rank_start[current_rank]);
                    cols[counter++] = mapped_dof;
                    _col_idx[counter_col_idx++] = mapped_dof;
                  } else if(GO(col) < rank_total_end[current_rank])
                  {
                    // Pressure DOF
                    GO pres_offset = col - rank_velo_end[current_rank];
                    GO mapped_dof = _num_global_velo_dofs + _all_first_owned_pres_dof[current_rank] + pres_offset;
                    cols[counter++] = mapped_dof;
                    _col_idx[counter_col_idx++] = mapped_dof;
                  }
                }
              }
              crsgraph->insertGlobalIndices(static_cast<GO>(_num_global_velo_dofs + _first_owned_pres_dof + (row - _num_owned_velo_dofs)),
                                            cols.view(0, static_cast<AGS>(counter)));
            }
            crsgraph->fillComplete();

            if(_use_timer)
            {
              _stacked_timer->stop("TpetraCoreStokes::init_core::build_crsgraph_system");
              _stacked_timer->start("TpetraCoreStokes::init_core::build_system_matrix");
            }

            /// create matrix
            _matrix = Teuchos::rcp(new Tpetra::CrsMatrix<SC, LO, GO, NO>(crsgraph, _params));

            // create vectors
            const LO num_vectors = 1;
            _vec_def = Teuchos::rcp(new Tpetra::MultiVector<SC, LO, GO, NO>(_map, num_vectors));
            _vec_cor = Teuchos::rcp(new Tpetra::MultiVector<SC, LO, GO, NO>(_map, num_vectors));

            if(_use_timer)
            {
              _stacked_timer->stop("TpetraCoreStokes::init_core::build_system_matrix");
              _stacked_timer->stop("TpetraCoreStokes::init_core");
            }

          }

          virtual void set_matrix_values(const double* vals) override
          {
            if(_use_timer)
              _stacked_timer->start("TpetraCoreStokes::SetValues");

            Teuchos::ArrayView<const LO> Indices;

            _matrix->resumeFill();

            AGS entries_start = 0;
            /// velocity
            for(IndexType row(0); row < _num_owned_velo_dofs; ++row)
            {
              const auto num_entries = static_cast<Teuchos::Array<SC>::size_type>(_num_nze[static_cast<AGS>(row)]);

              // copy vals to Teuchos::Array
              Teuchos::Array<SC> vals_ar(num_entries);
              for(Teuchos::Array<SC>::size_type j = 0; j < num_entries; ++j)
                vals_ar[j] = static_cast<SC>(vals[static_cast<VGS>(entries_start + j)]);

              _matrix->replaceGlobalValues(static_cast<GO>(_first_owned_velo_dof + row),
                                           _col_idx.view(entries_start, static_cast<AGS>(num_entries)),
                                           vals_ar.view(0, static_cast<Teuchos::Array<SC>::size_type>(num_entries)));

              entries_start += num_entries;
            }

            /// pressure
            for(IndexType row(0); row < _num_owned_pres_dofs; ++row)
            {
              const auto num_entries = static_cast<Teuchos::Array<SC>::size_type>(_num_nze[static_cast<AGS>(_num_owned_velo_dofs+row)]);

              // copy vals to Teuchos::Array
              Teuchos::Array<SC> vals_ar(num_entries);
              for(Teuchos::Array<SC>::size_type j = 0; j < num_entries; ++j)
                vals_ar[j] = static_cast<SC>(vals[static_cast<VGS>(entries_start + j)]);

              _matrix->replaceGlobalValues(static_cast<GO>(_num_global_velo_dofs + _first_owned_pres_dof + row),
                                           _col_idx.view(entries_start, static_cast<AGS>(num_entries)),
                                           vals_ar.view(0, static_cast<Teuchos::Array<SC>::size_type>(num_entries)));

              entries_start += num_entries;
            }

            _matrix->fillComplete();

            if(this->_use_timer)
            {
              this->_stacked_timer->stop("TpetraCoreStokes::SetValues");
            }
          }
      };

      void* create_core_stokes(const void* comm, IndexType dof_offset, IndexType num_owned_dofs, IndexType num_owned_velo_dofs, IndexType num_owned_pres_dofs, IndexType first_owned_velo_dof, IndexType first_owned_pres_dof, IndexType num_global_velo_dofs, IndexType num_global_pres_dofs, const unsigned int* row_ptr, const unsigned int* col_idx, const Trilinos::FROSchParameterList & params)
      {
        TpetraCoreStokes *core = new TpetraCoreStokes(comm, dof_offset, num_owned_dofs, num_owned_velo_dofs, num_owned_pres_dofs, first_owned_velo_dof, first_owned_pres_dof, num_global_velo_dofs, num_global_pres_dofs, params);
        core->init_core(row_ptr, col_idx);
        return core;
      }

      void* create_core_stokes(const void* comm, IndexType dof_offset, IndexType num_owned_dofs, IndexType num_owned_velo_dofs, IndexType num_owned_pres_dofs, IndexType first_owned_velo_dof, IndexType first_owned_pres_dof, IndexType num_global_velo_dofs, IndexType num_global_pres_dofs, const unsigned long* row_ptr, const unsigned long* col_idx, const Trilinos::FROSchParameterList & params)
      {
        TpetraCoreStokes *core = new TpetraCoreStokes(comm, dof_offset, num_owned_dofs, num_owned_velo_dofs, num_owned_pres_dofs, first_owned_velo_dof, first_owned_pres_dof, num_global_velo_dofs, num_global_pres_dofs, params);
        core->init_core(row_ptr, col_idx);
        return core;
      }

      void* create_core_stokes(const void* comm, IndexType dof_offset, IndexType num_owned_dofs, IndexType num_owned_velo_dofs, IndexType num_owned_pres_dofs, IndexType first_owned_velo_dof, IndexType first_owned_pres_dof, IndexType num_global_velo_dofs, IndexType num_global_pres_dofs, const unsigned long long* row_ptr, const unsigned long long* col_idx, const Trilinos::FROSchParameterList & params)
      {
        TpetraCoreStokes *core = new TpetraCoreStokes(comm, dof_offset, num_owned_dofs, num_owned_velo_dofs, num_owned_pres_dofs, first_owned_velo_dof, first_owned_pres_dof, num_global_velo_dofs, num_global_pres_dofs, params);
        core->init_core(row_ptr, col_idx);
        return core;
      }

      void destroy_core_stokes(void* core)
      {
        delete reinterpret_cast<TpetraCoreStokes*>(core);
      }

      /* *********************************************************************************************************** */
      /* *********************************************************************************************************** */
      /* *********************************************************************************************************** */

      typedef struct FROSchPrecond
      {
        using SC = TpetraCore::SC;
        Teuchos::RCP<Thyra::PreconditionerFactoryBase<SC>> pfbFactory;
        Teuchos::RCP<const Thyra::LinearOpBase<SC>> K_thyra;
        Teuchos::RCP<Thyra::PreconditionerBase<SC>> ThyraPrec;
      } FROSchPrecond;

      //----------------------------------------
      void* create_frosch(void* cr)
      {
        using SC = TpetraCore::SC;
        using LO = TpetraCore::LO;
        using GO = TpetraCore::GO;
        using NO = TpetraCore::NO;

        TpetraCore *core = reinterpret_cast<TpetraCore *>(cr);

        Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
        Stratimikos::enableFROSch<SC, LO, GO, NO>(linearSolverBuilder);
        // set parameter list
        linearSolverBuilder.setParameterList(core->_params);

        Teuchos::RCP<Thyra::PreconditionerFactoryBase<SC>> pfbFactory = linearSolverBuilder.createPreconditioningStrategy("");

        if(core->_use_timer)
        {
          Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
          pfbFactory->setOStream(out);
          pfbFactory->setVerbLevel(Teuchos::VERB_HIGH);
        }

        // system matrix and right-hand side to Thyra interface
        Teuchos::RCP<const Tpetra::RowMatrix<SC, LO, GO, NO>> tpRowMat  = Teuchos::rcp_dynamic_cast<const Tpetra::RowMatrix<SC, LO, GO, NO>>(core->_matrix, true);
        Teuchos::RCP<const Tpetra::Operator<SC, LO, GO, NO>> tpOperator = Teuchos::rcp_dynamic_cast<const Tpetra::Operator<SC, LO, GO, NO>> (tpRowMat, true);
        // thyraOp = Thyra::createConstLinearOp(tpOperator);

        // // Xpetra::CrsMatrixWrap<SC, LO, GO, NO> K(core->_matrix);
        // // Teuchos::RCP<const Thyra::LinearOpBase<SC>> K_thyra      = FROSch::ThyraUtils<SC,LO,GO,NO>::toThyra(K.getCrsMatrix());
        // Teuchos::RCP<const Thyra::LinearOpBase<SC>> K_thyra      = Thyra::createConstLinearOp(tpOperator);
        // Teuchos::RCP<Thyra::PreconditionerBase<SC>> ThyraPrec    = Thyra::prec(*pfbFactory, K_thyra);
        // Teuchos::RCP<Thyra::LinearOpBase<SC>>       LinearPrecOp = ThyraPrec->getNonconstUnspecifiedPrecOp();

        // /**
        //  *  It is a bit weird to use a raw pointer to a smart pointer.
        //  *  It would be nicer to use a pointer to Thyra::LinearOpWithSolveBase<SC>, but
        //  *  such things are not intended by the Thyra interface
        //  *  or at least I'm to stupid to achieve this with a little effort.
        //  **/
        // Teuchos::RCP<Thyra::LinearOpBase<SC>> *precond = new Teuchos::RCP<Thyra::LinearOpBase<SC>>();
        // *precond = LinearPrecOp;

        FROSchPrecond *precond = new FROSchPrecond;
        precond->pfbFactory = pfbFactory;
        precond->K_thyra = Thyra::createConstLinearOp(tpOperator);
        precond->ThyraPrec = Thyra::prec(*pfbFactory, precond->K_thyra);

        return precond;
      }

      void reinit_frosch(void *pre)
      {
        using SC = TpetraCore::SC;
        FROSchPrecond *precond = reinterpret_cast<FROSchPrecond *>(pre);
        Thyra::initializePrec<SC>(*(precond->pfbFactory), precond->K_thyra, precond->ThyraPrec.ptr());
      }

      void destroy_frosch(void *pre)
      {
        // using SC = TpetraCore::SC;
        // Teuchos::RCP<Thyra::LinearOpBase<SC>> *precond = reinterpret_cast<Teuchos::RCP<Thyra::LinearOpBase<SC>> *>(pre);
        // delete precond;
        FROSchPrecond *precond = reinterpret_cast<FROSchPrecond*>(pre);
        delete precond;
      }

      void solve_frosch(void* cr, void* pre)
      {
        using SC = TpetraCore::SC;
        using LO = TpetraCore::LO;
        using GO = TpetraCore::GO;
        using NO = TpetraCore::NO;

        TpetraCore *core = reinterpret_cast<TpetraCore *>(cr);
        FROSchPrecond *fprecond = reinterpret_cast<FROSchPrecond *>(pre);
        // Teuchos::RCP<Thyra::LinearOpBase<SC>> *precond = reinterpret_cast<Teuchos::RCP<Thyra::LinearOpBase<SC>> *>(pre);

        auto thyTpMap  = Thyra::tpetraVectorSpace<SC, LO, GO, NO>(core->_vec_cor->getMap());
        auto thyDomMap = Thyra::tpetraVectorSpace<SC, LO, GO, NO>(Tpetra::createLocalMapWithNode<LO, GO, NO>(core->_vec_cor->getNumVectors(), core->_vec_cor->getMap()->getComm()) );
        auto thyra_X   = Teuchos::rcp(new Thyra::TpetraMultiVector<SC, LO, GO, NO>());
        thyra_X->initialize(thyTpMap, thyDomMap, core->_vec_cor);

        thyTpMap  = Thyra::tpetraVectorSpace<SC, LO, GO, NO>(core->_vec_def->getMap());
        thyDomMap = Thyra::tpetraVectorSpace<SC, LO, GO, NO>(Tpetra::createLocalMapWithNode<LO, GO, NO>(core->_vec_def->getNumVectors(), core->_vec_def->getMap()->getComm()) );
        auto thyra_B   = Teuchos::rcp(new Thyra::TpetraMultiVector<SC, LO, GO, NO>());
        thyra_B->initialize(thyTpMap, thyDomMap, core->_vec_def);

        //Teuchos::RCP<Thyra::MultiVectorBase<SC>>       thyra_X = Teuchos::rcp_const_cast<Thyra::MultiVectorBase<SC> >(FROSch::ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector(core->_vec_cor));
        //Teuchos::RCP<const Thyra::MultiVectorBase<SC>> thyra_B = FROSch::ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector((core->_vec_def));

        // Thyra::apply<SC>(*(*precond), Thyra::NOTRANS, *(thyra_B.getConst()), thyra_X.ptr());
        Teuchos::RCP<Thyra::LinearOpBase<SC>> LinearPrecOp = fprecond->ThyraPrec->getNonconstUnspecifiedPrecOp();
        Thyra::apply<SC>(*LinearPrecOp, Thyra::NOTRANS, *(thyra_B.getConst()), thyra_X.ptr());
      }
    } // namespace Trilinos
  } // namespace Solver
} // namespace FEAT
#else
// insert dummy function to suppress linker warnings
void dummy_frosch_function() {}
#endif // FEAT_HAVE_TRILINOS
