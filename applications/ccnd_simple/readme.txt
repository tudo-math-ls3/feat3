===================================================================================================
Simple Monolithic Steady-State/Unsteady Navier-Stokes Solver Applications
===================================================================================================

This directory houses a set of both steady-state and unsteady Navier-Stokes solvers, which are
based upon a monolithic solver, i.e. a solver which solves both velocity and pressure components
simultaneously rather than by a splitting approach. The nonlinear systems are usually solved by
a Newton iteration (or a Picard iteration, if the user wishes so), which uses a Multigrid-Vanka
solver for the solution of the linearized Jacobian system.

To keep things simple, all applications perform a spatial discretization by using the Q2/P1dc
finite element space pair. Furthermore, all unsteady applications use the BDF(2) time stepping
scheme, unless explicitly stated otherwise in the application.

The actual code of the applications is split across several (non-template) classes, each of which
resides in its own source file. Typically, each application can be compiled as a 2D or a 3D
application, so that all class templates can instantiated for a compilation, which should make
development and debugging of these applications simpler than in an fully templatized code.

This document offers a rough description of the files that can be found in this directory and which
are shared by the various applications. For details about a particular application, please check
the comment block at the beginning of the application.

Important: All the apps in this directory require at least the 'fparser' third-party library to
be built; additionally, all apps use UMFPACK as their coarse grid solver by default, so it is
also recommended (but not mandatory) to include the UMFPACK third-party library as well.

---------------------------------------------------------------------------------------------------

The following files are common to (almost) all applications:

------------
* base.hpp *
------------
This is the base header for all applications in this directory. This header defines the most basic
classes used throughout the entire code base, such as mesh types, FE spaces and linear algebra
containers, as well as some auxiliary helper functions. This header also defines whether the
application is compiled as a 2D or a 3D application by checking for the definition of the
FEAT_CCND_SIMPLE_3D preprocessor macro. This header is included by any other header and source
file in this directory, either directly or indirectly.


-------------------
* domain.hpp/.cpp *
-------------------
This header/source file pair defines the DomainLevel class, which is a required class for the use
of the Control::Domain::PartiDomainControl class, as well as a custom DomainControl class, which
derives from the PartiDomainControl class. The DomainLevel class merely defines the basic mandatory
contents of a domain level, i.e. the trafo and the FE spaces. The DomainControl class extends the
PartiDomainControl with some additional output.
The DomainControl class is used by the application to read in the mesh from a mesh file and to
partition and distribute it to the individual processes in an MPI-parallel simulation.

-------------------------
* stokes_level.hpp/.cpp *
-------------------------
This header/source file pair defines the StokesLevel class, which is derived from an instance of
the Control::StokesBlockedCombinedSystemLevel class template. This class contains all linear
algebra containers which are required on each level of the multigrid hierarchy, i.e. all gates,
muxers, matrices, filters and grid transfers. This class also provides a set of member functions
for the assembly of matrices and filters.

--------------------------
* steady_solver.hpp/.cpp *
--------------------------
This header/source file pair defines the SteadySolver class, which implements the non-linear
Newton-/Picard-solver along with the setup of its linear Multigrid-AmaVanka preconditioner/solver.
This class is also responsible for assembling the non-linear Burgers matrices as well as the
non-linear defects. This class also contains some helper functions for the assembly of boundary
condition filters, the assembly of right-hand-side vectors and the interpolation of initial
solution vectors. Finally, this class also offers functionality for the I/O of joined solution
vectors.
This class is also the base-class for the UnsteadySolver class.

----------------------------
* unsteady_solver.hpp/.cpp *
----------------------------
This header/source file pair defines the UnsteadySolver class, which derives from the SteadySolver
class and extends it by adding an outer BDF(2) time stepping scheme. This class overrides the
matrix and vector assemblies of the base-class to incorporate the terms from the time-stepping
scheme and it provides the implementation of the basic time stepping scheme. This class also
implements a checkpoint/restart system, which allows the user to restart (or better: to continue)
a previously aborted/killed/crashed/stopped simulation by reading in a checkpoint from a file --
assuming that the previous simulation did save checkpoints periodically.

-----------------------
* vtk_writer.hpp/.cpp *
-----------------------
This header/source file pair defines the VtkWriter class, which is an extension of the very simple
Geometry::ExportVTK class template. The VtkWriter class allows to export of the solutions on a
once refined mesh, which is generally useful since the spatial discretization uses a Q2/P1dc
FE space pair, which contains more DOFs than could be exported into a VTK file on the same level.

---------------------------------------------------------------------------------------------------

\author Peter Zajac
