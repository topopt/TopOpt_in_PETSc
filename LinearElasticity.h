#ifndef __LINEARELASTICITY__
#define __LINEARELASTICITY__

#include <fstream>
#include <iostream>
#include <math.h>
#include <petsc.h>
#include <petsc/private/dmdaimpl.h>

/*
 Authors: Niels Aage, Erik Andreassen, Boyan Lazarov, August 2013
 Updated: June 2019, Niels Aage
 Copyright (C) 2013-2019,

 Disclaimer:
 The authors reserves all rights but does not guaranty that the code is
 free from errors. Furthermore, we shall not be liable in any event
 caused by the use of the program.
*/

class LinearElasticity {

  public:
    // Constructor
    LinearElasticity(DM da_nodes);

    // Destructor
    ~LinearElasticity();

    //  Compute objective and constraints and sensitivities at once: GOOD FOR
    //  SELF_ADJOINT PROBLEMS
    PetscErrorCode ComputeObjectiveConstraintsSensitivities(PetscScalar* fx, PetscScalar* gx, Vec dfdx, Vec dgdx,
                                                            Vec xPhys, PetscScalar Emin, PetscScalar Emax,
                                                            PetscScalar penal, PetscScalar volfrac);

    // Compute objective and constraints for the optimiation
    PetscErrorCode ComputeObjectiveConstraints(PetscScalar* fx, PetscScalar* gx, Vec xPhys, PetscScalar Emin,
                                               PetscScalar Emax, PetscScalar penal, PetscScalar volfrac);

    // Compute sensitivities
    PetscErrorCode ComputeSensitivities(Vec dfdx, Vec dgdx, Vec xPhys, PetscScalar Emin, PetscScalar Emax,
                                        PetscScalar penal,
                                        PetscScalar volfrac); // needs ....

    // Restart writer
    PetscErrorCode WriteRestartFiles();

    // Get pointer to the FE solution
    Vec GetStateField() { return (U); };

    // Get pointer to DMDA
    DM GetDM() { return (da_nodal); };

    // Logical mesh
    DM da_nodal; // Nodal mesh

  private:
    // Logical mesh
    PetscInt    nn[3]; // Number of nodes in each direction
    PetscInt    ne[3]; // Number of elements in each direction
    PetscScalar xc[6]; // Domain coordinates

    // Linear algebra
    Mat         K;           // Global stiffness matrix
    Vec         U;           // Displacement vector
    Vec         RHS;         // Load vector
    Vec         N;           // Dirichlet vector (used when imposing BCs)
    PetscScalar KE[24 * 24]; // Element stiffness matrix
    // Solver
    KSP         ksp; // Pointer to the KSP object i.e. the linear solver+prec
    PetscInt    nlvls;
    PetscScalar nu; // Possions ratio

    // Set up the FE mesh and data structures
    PetscErrorCode SetUpLoadAndBC(DM da_nodes);

    // Solve the FE problem
    PetscErrorCode SolveState(Vec xPhys, PetscScalar Emin, PetscScalar Emax, PetscScalar penal);

    // Assemble the stiffness matrix
    PetscErrorCode AssembleStiffnessMatrix(Vec xPhys, PetscScalar Emin, PetscScalar Emax, PetscScalar penal);

    // Start the solver
    PetscErrorCode SetUpSolver();

    // Routine that doesn't change the element type upon repeated calls
    PetscErrorCode DMDAGetElements_3D(DM dm, PetscInt* nel, PetscInt* nen, const PetscInt* e[]);

    // Methods used to assemble the element stiffness matrix
    PetscInt    Hex8Isoparametric(PetscScalar* X, PetscScalar* Y, PetscScalar* Z, PetscScalar nu, PetscInt redInt,
                                  PetscScalar* ke);
    PetscScalar Dot(PetscScalar* v1, PetscScalar* v2, PetscInt l);
    void        DifferentiatedShapeFunctions(PetscScalar xi, PetscScalar eta, PetscScalar zeta, PetscScalar* dNdxi,
                                             PetscScalar* dNdeta, PetscScalar* dNdzeta);
    PetscScalar Inverse3M(PetscScalar J[][3], PetscScalar invJ[][3]);

    // Restart
    PetscBool   restart, flip;
    std::string filename00, filename01;

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
