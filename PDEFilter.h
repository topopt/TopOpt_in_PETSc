#ifndef PDE_FILTER_H
#define PDE_FILTER_H
#include "TopOpt.h"
#include <petsc.h>

/* -----------------------------------------------------------------------------
Authors: Niels Aage, Erik Andreassen, Boyan Lazarov, August 2013
Updated: June 2019, Niels Aage
Copyright (C) 2013-2019,

This PDEFilter implementation is licensed under Version 2.1 of the GNU
Lesser General Public License.

This MMA implementation is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This Module is distributed in the hope that it will be useful,implementation
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this Module; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
-------------------------------------------------------------------------- */

class PDEFilt {

  public:
    PDEFilt(DM da_nodes, PetscScalar rmin);
    ~PDEFilt();

    PetscErrorCode FilterProject(Vec XX, Vec F);
    PetscErrorCode Gradients(Vec OS, Vec FS);

  private:
    PetscInt    nn[3];   // Number of nodes in each direction
    PetscInt    ne[3];   // Number of elements in each direction
    PetscScalar xc[6];   // Domain coordinates
    PetscScalar elemVol; // element volume

    PetscScalar R; // filter parameter

    PetscScalar KF[8 * 8]; // PDE filter stiffness matrix
    PetscScalar TF[8];     // PDE filter transformation matrix

    PetscInt nloc; // Number of local nodes?

    PetscInt nlvls; // Number of multigrid levels for the filter

    DM da_nodal;
    DM da_element;

    Mat K;   // Global stiffness matrix
    Mat T;   // Transformation matrix   RHS=T*X
    Vec RHS; // Load vector - nodal
    Vec U;
    Vec X; // filtered filed - element

    KSP ksp; // linear solver

    void PDEFilterMatrix(PetscScalar dx, PetscScalar dy, PetscScalar dz, PetscScalar R, PetscScalar* KK,
                         PetscScalar* T);

    PetscErrorCode DMDAGetElements_3D(DM dm, PetscInt* nel, PetscInt* nen, const PetscInt* e[]);

    void MatAssemble(); // assemble K and T
                        // RHS = T*elvol*RHO

    PetscErrorCode SetUpSolver();
    PetscErrorCode Free();
};

#endif
