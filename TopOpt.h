#ifndef TOPOPT_H
#define TOPOPT_H

#include <petsc.h>
//#include <petsc-private/dmdaimpl.h>
#include <petsc/private/dmdaimpl.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <MMA.h>

/*
 Authors: Niels Aage, Erik Andreassen, Boyan Lazarov, August 2013

 Disclaimer:
 The authors reserves all rights but does not guaranty that the code is
 free from errors. Furthermore, we shall not be liable in any event
 caused by the use of the program.
*/


/*
 *
 * Parameter container for the topology optimization problem
 *
 * min_x fx
 * s.t. gx_j <= 0, j=1..m
 *      xmin_i <= x_i <= xmax_i, i=1..n
 *
 * with filtering and a volume constraint
 *
 */

class TopOpt {

public:

	// Constructor/Destructor
	TopOpt();
	~TopOpt();

	// Method to allocate MMA with/without restarting
	PetscErrorCode AllocateMMAwithRestart(PetscInt *itr, MMA **mma);
	PetscErrorCode WriteRestartFiles(PetscInt *itr, MMA *mma);

	// Physical domain variables
	PetscScalar xc[6]; // Domain coordinates
	PetscScalar dx,dy,dz; // Element size
	PetscInt nxyz[3]; // Number of nodes in each direction
	PetscInt nlvls; // Number of multigrid levels
	PetscScalar nu; // Poisson's ratio
	// Nodal mesh (contains physics)
	DM da_nodes;
	// element mesh (contains design)
	DM da_elem;

	// Optimization parameters
	PetscInt n; // Total number of design variables
	PetscInt nloc	; // Local number of local nodes?
	PetscInt m; // Number of constraints
	PetscScalar fx; // Objective value
	PetscScalar fscale; // Scaling factor for objective
	PetscScalar *gx; // Array with constraint values
	PetscScalar Xmin; // Min. value of design variables
	PetscScalar Xmax; // Max. value of design variables

	PetscScalar movlim; // Max. change of design variables
	PetscScalar volfrac; // Volume fraction
	PetscScalar penal; // Penalization parameter
	PetscScalar Emin, Emax; // Modified SIMP, max and min E

	PetscScalar rmin; // filter radius

	PetscInt maxItr; // Max iterations
	PetscInt filter; // Filter type

	Vec x; // Design variables
	Vec xPhys; // Physical variables (filtered x)
	Vec dfdx; // Sensitivities of objective
	Vec xmin, xmax; // Vectors with max and min values of x
	Vec xold; // x from previous iteration
	Vec *dgdx; // Sensitivities of constraints (vector array)

	// Restart data for MMA:
	PetscBool restart, flip;
	std::string restdens_1,restdens_2;
	Vec xo1, xo2, U, L;

private:
	// Allocate and set default values
	PetscErrorCode SetUp();

	PetscErrorCode SetUpMESH();
	PetscErrorCode SetUpOPT();

	// Restart filenames
	std::string filename00, filename00Itr, filename01, filename01Itr;

	// File existence
	inline PetscBool fexists(const std::string& filename) {
	      std::ifstream ifile(filename.c_str());
	      if (ifile) {
		return PETSC_TRUE;
	      }
	      return PETSC_FALSE;
	}

};

#endif
