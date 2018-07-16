#include <kernel/util/runtime.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/geometry/common_factories.hpp>


#if defined(FEAT_HAVE_HYPRE) && defined(FEAT_HAVE_MPI)
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"

void run(int /*argc*/, char** /*argv*/)
{
  // taken from example #5 of hypre
  FEAT::Dist::Comm comm = FEAT::Dist::Comm::world();
  int myid = comm.rank();
  int num_procs = comm.size();

  HYPRE_Int i;
  int N, n;

  HYPRE_Int ilower, iupper;
  HYPRE_Int local_size, extra;

  double h, h2;

  HYPRE_IJMatrix A;
  HYPRE_ParCSRMatrix parcsr_A;
  HYPRE_IJVector b;
  HYPRE_ParVector par_b;
  HYPRE_IJVector x;
  HYPRE_ParVector par_x;
  HYPRE_Solver solver;

  /* Default problem parameters */
  n = 33;

  /* Preliminaries: want at least one processor per row */
  if (n*n < num_procs) n = int(sqrt(num_procs)) + 1;
  N = n*n; /* global number of rows */
  h = 1.0/(n+1); /* mesh size*/
  h2 = h*h;

  /* Each processor knows only of its own rows - the range is denoted by ilower
    and upper.  Here we partition the rows. We account for the fact that
    N may not divide evenly by the number of processors. */
  local_size = N/num_procs;
  extra = N - local_size*num_procs;

  ilower = local_size*myid;
  ilower += hypre_min(myid, extra);

  iupper = local_size*(myid+1);
  iupper += hypre_min(myid+1, extra);
  iupper = iupper - 1;

  /* How many rows do I have? */
  local_size = iupper - ilower + 1;

  /* Create the matrix.
    Note that this is a square matrix, so we indicate the row partition
    size twice (since number of rows = number of cols) */
  HYPRE_IJMatrixCreate(comm.mpi_comm(), ilower, iupper, ilower, iupper, &A);

  /* Choose a parallel csr format storage (see the User's Manual) */
  HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);

  /* Initialize before setting coefficients */
  HYPRE_IJMatrixInitialize(A);

  /* Now go through my local rows and set the matrix entries.
    Each row has at most 5 entries. For example, if n=3:

    A = [M -I 0; -I M -I; 0 -I M]
    M = [4 -1 0; -1 4 -1; 0 -1 4]

    Note that here we are setting one row at a time, though
    one could set all the rows together (see the User's Manual).
  */
  {
    HYPRE_Int nnz;
    double values[5];
    HYPRE_Int cols[5];

    for (i = ilower; i <= iupper; i++)
    {
        nnz = 0;

        /* The left identity block:position i-n */
        if ((i-n)>=0)
        {
          cols[nnz] = i-n;
          values[nnz] = -1.0;
          nnz++;
        }

        /* The left -1: position i-1 */
        if (i%n)
        {
          cols[nnz] = i-1;
          values[nnz] = -1.0;
          nnz++;
        }

        /* Set the diagonal: position i */
        cols[nnz] = i;
        values[nnz] = 4.0;
        nnz++;

        /* The right -1: position i+1 */
        if ((i+1)%n)
        {
          cols[nnz] = i+1;
          values[nnz] = -1.0;
          nnz++;
        }

        /* The right identity block:position i+n */
        if ((i+n)< N)
        {
          cols[nnz] = i+n;
          values[nnz] = -1.0;
          nnz++;
        }

        /* Set the values for row i */
        HYPRE_IJMatrixSetValues(A, 1, &nnz, &i, cols, values);
    }
  }

  /* Assemble after setting the coefficients */
  HYPRE_IJMatrixAssemble(A);

  /* Note: for the testing of small problems, one may wish to read
    in a matrix in IJ format (for the format, see the output files
    from the -print_system option).
    In this case, one would use the following routine:
    HYPRE_IJMatrixRead( <filename>, MPI_COMM_WORLD,
                        HYPRE_PARCSR, &A );
    <filename>  = IJ.A.out to read in what has been printed out
    by -print_system (processor numbers are omitted).
    A call to HYPRE_IJMatrixRead is an *alternative* to the
    following sequence of HYPRE_IJMatrix calls:
    Create, SetObjectType, Initialize, SetValues, and Assemble
  */


  /* Get the parcsr matrix object to use */
  HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);


  /* Create the rhs and solution */
  HYPRE_IJVectorCreate(comm.mpi_comm(), ilower, iupper,&b);
  HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(b);

  HYPRE_IJVectorCreate(comm.mpi_comm(), ilower, iupper,&x);
  HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
  HYPRE_IJVectorInitialize(x);

  /* Set the rhs values to h^2 and the solution to zero */
  {
    double *rhs_values, *x_values;
    HYPRE_Int *rows;

    rhs_values = (double*) calloc(std::size_t(local_size), sizeof(double));
    x_values = (double*) calloc(std::size_t(local_size), sizeof(double));
    rows = (HYPRE_Int*) calloc(std::size_t(local_size), sizeof(HYPRE_Int));

    for (i=0; i<local_size; i++)
    {
        rhs_values[i] = h2;
        x_values[i] = 0.0;
        rows[i] = ilower + i;
    }

    HYPRE_IJVectorSetValues(b, local_size, rows, rhs_values);
    HYPRE_IJVectorSetValues(x, local_size, rows, x_values);

    free(x_values);
    free(rhs_values);
    free(rows);
  }


  HYPRE_IJVectorAssemble(b);
  /*  As with the matrix, for testing purposes, one may wish to read in a rhs:
      HYPRE_IJVectorRead( <filename>, MPI_COMM_WORLD,
                                HYPRE_PARCSR, &b );
      as an alternative to the
      following sequence of HYPRE_IJVectors calls:
      Create, SetObjectType, Initialize, SetValues, and Assemble
  */
  HYPRE_IJVectorGetObject(b, (void **) &par_b);

  HYPRE_IJVectorAssemble(x);
  HYPRE_IJVectorGetObject(x, (void **) &par_x);


  /* Choose a solver and solve the system */

  /* AMG */
  {
    HYPRE_Int num_iterations;
    double final_res_norm;

    /* Create solver */
    HYPRE_BoomerAMGCreate(&solver);

    /* Set some parameters (See Reference Manual for more parameters) */
    HYPRE_BoomerAMGSetPrintLevel(solver, 3);  /* print solve info + parameters */
    HYPRE_BoomerAMGSetOldDefault(solver); /* Falgout coarsening with modified classical interpolation */
    HYPRE_BoomerAMGSetRelaxType(solver, 3);   /* G-S/Jacobi hybrid relaxation */
    HYPRE_BoomerAMGSetRelaxOrder(solver, 1);   /* Uses C/F relaxation */
    HYPRE_BoomerAMGSetNumSweeps(solver, 1);   /* Sweeeps on each level */
    HYPRE_BoomerAMGSetMaxLevels(solver, 20);  /* maximum number of levels */
    HYPRE_BoomerAMGSetTol(solver, 1e-7);      /* conv. tolerance */

    /* Now setup and solve! */
    HYPRE_BoomerAMGSetup(solver, parcsr_A, par_b, par_x);
    HYPRE_BoomerAMGSolve(solver, parcsr_A, par_b, par_x);

    /* Run info - needed logging turned on */
    HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
    HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
    if (myid == 0)
    {
        printf("\n");
        printf("Iterations = %lld\n", num_iterations);
        printf("Final Relative Residual Norm = %e\n", final_res_norm);
        printf("\n");
    }

    /* Destroy solver */
    HYPRE_BoomerAMGDestroy(solver);
  }

  /* Clean up */
  HYPRE_IJMatrixDestroy(A);
  HYPRE_IJVectorDestroy(b);
  HYPRE_IJVectorDestroy(x);
}

int main(int argc, char** argv)
{
  FEAT::Runtime::initialise(argc, argv);
  run(argc, argv);
  return FEAT::Runtime::finalise();
}

#else

int main()
{
  return 0;
}

#endif // FEAT_HAVE_HYPRE && FEAT_HAVE_MPI
