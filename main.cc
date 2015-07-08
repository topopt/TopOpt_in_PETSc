#include <petsc.h>
#include <TopOpt.h>
#include <LinearElasticity.h>
#include <MMA.h>
#include <Filter.h>
#include <MPIIO.h>
#include <mpi.h>
/*
 Authors: Niels Aage, Erik Andreassen, Boyan Lazarov, August 2013

 Disclaimer:
 The authors reserves all rights but does not guaranty that the code is
 free from errors. Furthermore, we shall not be liable in any event
 caused by the use of the program.
*/

static char help[] = "3D TopOpt using KSP-MG on PETSc's DMDA (structured grids) \n";

int main(int argc, char *argv[]){

  // Error code for debugging
  PetscErrorCode ierr;

  // Initialize PETSc / MPI and pass input arguments to PETSc
  PetscInitialize(&argc,&argv,PETSC_NULL,help);

// STEP 1: THE OPTIMIZATION PARAMETERS, DATA AND MESH (!!! THE DMDA !!!)
  TopOpt *opt = new TopOpt();

// STEP 2: THE PHYSICS
  LinearElasticity *physics = new LinearElasticity(opt);

// STEP 3: THE FILTERING
  Filter *filter = new Filter(opt);

// STEP 4: VISUALIZATION USING VTK
  MPIIO *output = new MPIIO(opt->da_nodes,3,"ux, uy, uz",2,"x, xPhys");

// STEP 5: THE OPTIMIZER MMA
  MMA *mma;
  PetscInt itr=0;
  opt->AllocateMMAwithRestart(&itr, &mma); // allow for restart !

// STEP 6: FILTER THE INITIAL DESIGN/RESTARTED DESIGN
  ierr = filter->FilterProject(opt); CHKERRQ(ierr);

// STEP 7: OPTIMIZATION LOOP
  PetscScalar ch = 1.0;
  double t1,t2;
  while (itr < opt->maxItr && ch > 0.01){
	// Update iteration counter
	itr++;

	// start timer
	t1 = MPI_Wtime();

	// Compute (a) obj+const, (b) sens, (c) obj+const+sens
	ierr = physics->ComputeObjectiveConstraintsSensitivities(opt); CHKERRQ(ierr);

	// Compute objective scale
	if (itr==1){
		opt->fscale = 10.0/opt->fx;
	}
	// Scale objectie and sens
	opt->fx = opt->fx*opt->fscale;
	VecScale(opt->dfdx,opt->fscale);

	// Filter sensitivities (chainrule)
	ierr = filter->Gradients(opt); CHKERRQ(ierr);

	// Sets outer movelimits on design variables
	ierr = mma->SetOuterMovelimit(opt->Xmin,opt->Xmax,opt->movlim,opt->x,opt->xmin,opt->xmax); CHKERRQ(ierr);

	// Update design by MMA
	ierr = mma->Update(opt->x,opt->dfdx,opt->gx,opt->dgdx,opt->xmin,opt->xmax); CHKERRQ(ierr);

	// Inf norm on the design change
	ch = mma->DesignChange(opt->x,opt->xold);

	// Filter design field
	ierr = filter->FilterProject(opt); CHKERRQ(ierr);

	// stop timer
	t2 = MPI_Wtime();

	// Print to screen
	PetscPrintf(PETSC_COMM_WORLD,"It.: %i, obj.: %f, g[0]: %f, ch.: %f, time: %f\n",
				itr,opt->fx,opt->gx[0], ch,t2-t1);

	// Write field data: first 10 iterations and then every 20th
	if (itr<11 || itr%20==0){
		output->WriteVTK(opt->da_nodes,physics->GetStateField(),opt, itr);
	}

	// Dump data needed for restarting code at termination
	if (itr%3==0)	{
		opt->WriteRestartFiles(&itr, mma);
		physics->WriteRestartFiles();
	}
  }
// Write restart WriteRestartFiles
  opt->WriteRestartFiles(&itr, mma);
  physics->WriteRestartFiles();

// Dump final design
  output->WriteVTK(opt->da_nodes,physics->GetStateField(),opt, itr+1);

// STEP 7: CLEAN UP AFTER YOURSELF
  delete mma;
  delete output;
  delete filter;
  delete opt;
  delete physics;

  // Finalize PETSc / MPI
  PetscFinalize();
  return 0;
}
