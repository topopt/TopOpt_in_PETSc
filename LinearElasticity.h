#ifndef __LINEARELASTICITY__
#define __LINEARELASTICITY__

#include <petsc.h>
//#include <petsc-private/dmdaimpl.h>
#include <petsc/private/dmdaimpl.h>
#include <iostream>
#include <math.h>
#include <TopOpt.h>

/*
 Authors: Niels Aage, Erik Andreassen, Boyan Lazarov, August 2013

 Disclaimer:                                                              
 The authors reserves all rights but does not guaranty that the code is   
 free from errors. Furthermore, we shall not be liable in any event     
 caused by the use of the program.                                     
*/


class LinearElasticity{
  
public:
    // Constructor
    LinearElasticity(TopOpt *opt);
    
    // Destructor
    ~LinearElasticity();

    // Compute objective and constraints for the optimiation
    PetscErrorCode ComputeObjectiveConstraints(TopOpt *opt); 
    
    // Compute sensitivities
    PetscErrorCode ComputeSensitivities(TopOpt *opt); // needs ....
    
    //  Compute objective and constraints and sensitivities at once: GOOD FOR SELF_ADJOINT PROBLEMS
    PetscErrorCode ComputeObjectiveConstraintsSensitivities(TopOpt *opt);
    
    // Restart writer
    PetscErrorCode WriteRestartFiles();
    
    // Get pointer to the FE solution
    Vec GetStateField(){ return(U); };
    
private:
  
    PetscScalar KE[24*24]; // Element stiffness matrix
    Mat K; // Global stiffness matrix
    Vec U; // Displacement vector
    Vec RHS; // Load vector 
    Vec N; // Dirichlet vector (used when imposing BCs) 
    // Solver 
    KSP ksp;	// Pointer to the KSP object i.e. the linear solver+prec
  
    // Set up the FE mesh and data structures
    PetscErrorCode SetUpLoadAndBC(TopOpt *opt);
    
    // Solve the FE problem
    PetscErrorCode SolveState(TopOpt *opt);
    
    // Assemble the stiffness matrix	
    PetscErrorCode AssembleStiffnessMatrix(TopOpt *opt);
  
    // Start the solver
    PetscErrorCode SetUpSolver(TopOpt *opt);
    
    // Routine that doesn't change the element type upon repeated calls
    PetscErrorCode DMDAGetElements_3D(DM dm,PetscInt *nel,PetscInt *nen,const PetscInt *e[]);
	
    // Methods used to assemble the element stiffness matrix
    PetscInt Hex8Isoparametric(PetscScalar *X, PetscScalar *Y, PetscScalar *Z, PetscScalar nu, PetscInt redInt, PetscScalar *ke);
    PetscScalar Dot(PetscScalar *v1, PetscScalar *v2, PetscInt l);
    void DifferentiatedShapeFunctions(PetscScalar xi, PetscScalar eta, PetscScalar zeta, PetscScalar *dNdxi, PetscScalar *dNdeta, PetscScalar *dNdzeta);
    PetscScalar Inverse3M(PetscScalar J[][3], PetscScalar invJ[][3]);

    // Restart
    PetscBool restart, flip;
    std::string filename00,filename01;	
    
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
