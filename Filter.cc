#include "Filter.h"

/* -----------------------------------------------------------------------------
Authors: Niels Aage, Erik Andreassen, Boyan Lazarov, August 2013
Copyright (C) 2013-2014,

This Filter implementation is licensed under Version 2.1 of the GNU
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

Filter::Filter(DM da_nodes, Vec x, PetscInt filterT, PetscScalar Rin) {
    // Set all pointers to NULL
    H       = NULL;
    Hs      = NULL;
    da_elem = NULL;
    pdef    = NULL;

    // Get parameters
    R          = Rin;
    filterType = filterT;

    // Call the setup method
    SetUp(da_nodes, x);
}

Filter::~Filter() {
    // Deallocate data
    if (Hs != NULL) {
        VecDestroy(&Hs);
    }
    if (H != NULL) {
        MatDestroy(&H);
    }
    if (da_elem != NULL) {
        DMDestroy(&da_elem);
    }
    if (pdef != NULL) {
        delete pdef;
    }
    if (dx != NULL) {
        VecDestroy(&dx);
    }
}

// Filter design variables
PetscErrorCode Filter::FilterProject(Vec x, Vec xTilde, Vec xPhys, PetscBool projectionFilter, PetscScalar beta,
                                     PetscScalar eta) {
    PetscErrorCode ierr;

    // Filter the design variables or copy to xPhys
    // STANDARD FILTER
    if (filterType == 1) {
        // Filter the densitities
        ierr = MatMult(H, x, xTilde);
        CHKERRQ(ierr);
        VecPointwiseDivide(xTilde, xTilde, Hs);
    }
    // PDE FILTER
    else if (filterType == 2) {
        ierr = pdef->FilterProject(x, xTilde);
        CHKERRQ(ierr);
        // Check for bound violation: simple, but cheap check!
        PetscScalar* xp;
        PetscInt     locsiz;
        VecGetArray(xTilde, &xp);
        VecGetLocalSize(xTilde, &locsiz);
        for (PetscInt i = 0; i < locsiz; i++) {
            if (xp[i] < 0.0) {
                if (PetscAbsReal(xp[i]) > 1.0e-4) {
                    PetscPrintf(PETSC_COMM_WORLD,
                                "BOUND VIOLATION IN PDEFILTER - INCREASE RMIN OR MESH "
                                "RESOLUTION: xPhys = %f\n",
                                xp[i]);
                }
                xp[i] = 0.0;
            }
            if (xp[i] > 1.0) {
                if (PetscAbsReal(xp[i] - 1.0) > 1.0e-4) {
                    PetscPrintf(PETSC_COMM_WORLD,
                                "BOUND VIOLATION IN PDEFILTER - INCREASE RMIN OR MESH "
                                "RESOLUTION: xPhys = %f\n",
                                xp[i]);
                }
                xp[i] = 1.0;
            }
        }
        VecRestoreArray(xTilde, &xp);
    }
    // COPY IN CASE OF SENSITIVITY FILTER
    else {
        ierr = VecCopy(x, xTilde);
        CHKERRQ(ierr);
    }

    // Check for projection
    if (projectionFilter) {
        HeavisideFilter(xPhys, xTilde, beta, eta);
    } else {
        VecCopy(xTilde, xPhys);
    }

    return ierr;
}

// Filter the sensitivities
PetscErrorCode Filter::Gradients(Vec x, Vec xTilde, Vec dfdx, PetscInt m, Vec* dgdx, PetscBool projectionFilter,
                                 PetscScalar beta, PetscScalar eta) {

    PetscErrorCode ierr;
    // Cheinrule for projection filtering
    if (projectionFilter) {

        // Get correction
        ChainruleHeavisideFilter(dx, xTilde, beta, eta);

        PetscScalar *xt, *dg, *df, *dxp;
        PetscInt     locsiz;

        ierr = VecGetLocalSize(xTilde, &locsiz);
        CHKERRQ(ierr);
        ierr = VecGetArray(xTilde, &xt);
        CHKERRQ(ierr);
        ierr = VecGetArray(dx, &dxp);
        CHKERRQ(ierr);
        // Objective function
        ierr = VecGetArray(dfdx, &df);
        CHKERRQ(ierr);
        for (PetscInt j = 0; j < locsiz; j++) {
            df[j] = df[j] * dxp[j];
        }
        ierr = VecRestoreArray(dfdx, &df);
        CHKERRQ(ierr);
        // Run through all constraints
        for (PetscInt i = 0; i < m; i++) {
            ierr = VecGetArray(dgdx[i], &dg);
            CHKERRQ(ierr);
            // The eta item corresponding to the correct realization
            for (PetscInt j = 0; j < locsiz; j++) {
                dg[j] = dg[j] * dxp[j];
            }
            ierr = VecRestoreArray(dgdx[i], &dg);
            CHKERRQ(ierr);
        }
        ierr = VecRestoreArray(dx, &dxp);
        CHKERRQ(ierr);
        ierr = VecRestoreArray(dgdx[0], &dg);
        CHKERRQ(ierr);
        ierr = VecRestoreArray(xTilde, &xt);
        CHKERRQ(ierr);
    }

    // Chainrule/Filter for the sensitivities
    if (filterType == 0)
    // Filter the sensitivities, df,dg
    {
        Vec xtmp;
        ierr = VecDuplicate(xTilde, &xtmp);
        CHKERRQ(ierr);
        VecPointwiseMult(xtmp, dfdx, x);
        MatMult(H, xtmp, dfdx);
        VecPointwiseDivide(xtmp, dfdx, Hs);
        VecPointwiseDivide(dfdx, xtmp, x);
        VecDestroy(&xtmp);
    } else if (filterType == 1) {
        // Filter the densities, df,dg: STANDARD FILTER
        Vec xtmp;
        ierr = VecDuplicate(x, &xtmp);
        CHKERRQ(ierr);
        // dfdx
        VecPointwiseDivide(xtmp, dfdx, Hs);
        MatMult(H, xtmp, dfdx);
        // dgdx
        for (PetscInt i = 0; i < m; i++) {
            VecPointwiseDivide(xtmp, dgdx[i], Hs);
            MatMult(H, xtmp, dgdx[i]);
        }
        // tidy up
        VecDestroy(&xtmp);
    } else if (filterType == 2) {
        // Filter the densities, df,dg: PDE FILTER
        ierr = pdef->Gradients(dfdx, dfdx);
        CHKERRQ(ierr);
        for (PetscInt i = 0; i < m; i++) {
            ierr = pdef->Gradients(dgdx[i], dgdx[i]);
            CHKERRQ(ierr);
        }
    }

    return ierr;
}

PetscScalar Filter::GetMND(Vec x) {

    PetscScalar mnd, mndloc = 0.0;

    PetscScalar* xv;
    PetscInt     nelloc, nelglob;
    VecGetLocalSize(x, &nelloc);
    VecGetSize(x, &nelglob);

    // Compute power sum
    VecGetArray(x, &xv);
    for (PetscInt i = 0; i < nelloc; i++) {
        mndloc += 4 * xv[i] * (1.0 - xv[i]);
    }
    // Collect from procs
    MPI_Allreduce(&mndloc, &mnd, 1, MPIU_SCALAR, MPI_SUM, PETSC_COMM_WORLD);
    mnd = mnd / ((PetscScalar)nelglob);

    return mnd;
}

PetscErrorCode Filter::HeavisideFilter(Vec y, Vec x, PetscReal beta, PetscReal eta) {
    PetscErrorCode ierr;

    PetscScalar *yp, *xp;
    PetscInt     nelloc;
    VecGetLocalSize(x, &nelloc);
    ierr = VecGetArray(x, &xp);
    CHKERRQ(ierr);
    ierr = VecGetArray(y, &yp);
    CHKERRQ(ierr);

    for (PetscInt i = 0; i < nelloc; i++) {
        yp[i] = SmoothProjection(xp[i], beta, eta);
    }
    ierr = VecRestoreArray(x, &xp);
    CHKERRQ(ierr);
    ierr = VecRestoreArray(y, &yp);
    CHKERRQ(ierr);
}

PetscErrorCode Filter::ChainruleHeavisideFilter(Vec y, Vec x, PetscReal beta, PetscReal eta) {
    PetscErrorCode ierr;

    PetscScalar *yp, *xp;
    PetscInt     nelloc;
    VecGetLocalSize(x, &nelloc);
    ierr = VecGetArray(x, &xp);
    CHKERRQ(ierr);
    ierr = VecGetArray(y, &yp);
    CHKERRQ(ierr);

    for (PetscInt i = 0; i < nelloc; i++) {
        yp[i] = ChainruleSmoothProjection(xp[i], beta, eta);
    }
    ierr = VecRestoreArray(x, &xp);
    CHKERRQ(ierr);
    ierr = VecRestoreArray(y, &yp);
    CHKERRQ(ierr);
}

// Continuation function
PetscBool Filter::IncreaseBeta(PetscReal* beta, PetscReal betaFinal, PetscScalar gx, PetscInt itr, PetscReal ch) {

    PetscBool changeBeta = PETSC_FALSE;

    // Increase beta when fitting
    if ((ch < 0.01 || itr % 10 == 0) && beta[0] < betaFinal && gx < 0.000001) {
        changeBeta = PETSC_TRUE;
        if (beta[0] < 7) {
            beta[0] = beta[0] + 1;
        } else {
            beta[0] = beta[0] * 1.2;
        }
        if (beta[0] > betaFinal) {
            beta[0]    = betaFinal;
            changeBeta = PETSC_FALSE;
        }
        PetscPrintf(PETSC_COMM_WORLD, "Beta has been increased to: %f\n", beta[0]);
    }

    return changeBeta;
}

PetscErrorCode Filter::SetUp(DM da_nodes, Vec x) {

    PetscErrorCode ierr;

    VecDuplicate(x, &dx);
    VecSet(dx, 1.0);

    if (filterType == 0 || filterType == 1) {

        // Extract information from the nodal mesh
        PetscInt        M, N, P, md, nd, pd;
        DMBoundaryType  bx, by, bz;
        DMDAStencilType stype;
        ierr = DMDAGetInfo(da_nodes, NULL, &M, &N, &P, &md, &nd, &pd, NULL, NULL, &bx, &by, &bz, &stype);
        CHKERRQ(ierr);

        // Find the element size
        Vec lcoor;
        DMGetCoordinatesLocal(da_nodes, &lcoor);
        PetscScalar* lcoorp;
        VecGetArray(lcoor, &lcoorp);

        PetscInt        nel, nen;
        const PetscInt* necon;
        DMDAGetElements_3D(da_nodes, &nel, &nen, &necon);

        PetscScalar dx, dy, dz;
        // Use the first element to compute the dx, dy, dz
        dx = lcoorp[3 * necon[0 * nen + 1] + 0] - lcoorp[3 * necon[0 * nen + 0] + 0];
        dy = lcoorp[3 * necon[0 * nen + 2] + 1] - lcoorp[3 * necon[0 * nen + 1] + 1];
        dz = lcoorp[3 * necon[0 * nen + 4] + 2] - lcoorp[3 * necon[0 * nen + 0] + 2];
        VecRestoreArray(lcoor, &lcoorp);

        // Create the minimum element connectivity shit
        PetscInt ElemConn;
        // Check dx,dy,dz and find max conn for a given rmin
        ElemConn = (PetscInt)PetscMax(ceil(R / dx) - 1, PetscMax(ceil(R / dy) - 1, ceil(R / dz) - 1));
        ElemConn = PetscMin(ElemConn, PetscMin((M - 1) / 2, PetscMin((N - 1) / 2, (P - 1) / 2)));

        // The following is needed due to roundoff errors
        PetscInt tmp;
        MPI_Allreduce(&ElemConn, &tmp, 1, MPIU_INT, MPI_MAX, PETSC_COMM_WORLD);
        ElemConn = tmp;

        // Print to screen: mesh overlap!
        PetscPrintf(PETSC_COMM_WORLD, "# Filter radius rmin = %f results in a stencil of %i elements \n", R, ElemConn);

        // Find the geometric partitioning of the nodal mesh, so the element mesh
        // will coincide
        PetscInt* Lx = new PetscInt[md];
        PetscInt* Ly = new PetscInt[nd];
        PetscInt* Lz = new PetscInt[pd];

        // get number of nodes for each partition
        const PetscInt *LxCorrect, *LyCorrect, *LzCorrect;
        DMDAGetOwnershipRanges(da_nodes, &LxCorrect, &LyCorrect, &LzCorrect);

        // subtract one from the lower left corner.
        for (int i = 0; i < md; i++) {
            Lx[i] = LxCorrect[i];
            if (i == 0) {
                Lx[i] = Lx[i] - 1;
            }
        }
        for (int i = 0; i < nd; i++) {
            Ly[i] = LyCorrect[i];
            if (i == 0) {
                Ly[i] = Ly[i] - 1;
            }
        }
        for (int i = 0; i < pd; i++) {
            Lz[i] = LzCorrect[i];
            if (i == 0) {
                Lz[i] = Lz[i] - 1;
            }
        }

        // Create the element grid:
        DMDACreate3d(PETSC_COMM_WORLD, bx, by, bz, stype, M - 1, N - 1, P - 1, md, nd, pd, 1, ElemConn, Lx, Ly, Lz,
                     &da_elem);
        // Initialize
        DMSetFromOptions(da_elem);
        DMSetUp(da_elem);

        // Set the coordinates: from 0+dx/2 to xmax-dx/2 and so on
        PetscScalar xmax = (M - 1) * dx;
        PetscScalar ymax = (N - 1) * dy;
        PetscScalar zmax = (P - 1) * dz;
        DMDASetUniformCoordinates(da_elem, dx / 2.0, xmax - dx / 2.0, dy / 2.0, ymax - dy / 2.0, dz / 2.0,
                                  zmax - dz / 2.0);

        // Allocate and assemble
        DMCreateMatrix(da_elem, &H);
        DMCreateGlobalVector(da_elem, &Hs);

        // Set the filter matrix and vector
        DMGetCoordinatesLocal(da_elem, &lcoor);
        VecGetArray(lcoor, &lcoorp);
        DMDALocalInfo info;
        DMDAGetLocalInfo(da_elem, &info);
        // The variables from info that are used are described below:
        // -------------------------------------------------------------------------
        // sw = Stencil width
        // mx, my, mz = Global number of "elements" in each direction
        // xs, ys, zs = Starting point of this processor, excluding ghosts
        // xm, ym, zm = Number of grid points on this processor, excluding ghosts
        // gxs, gys, gzs = Starting point of this processor, including ghosts
        // gxm, gym, gzm = Number of grid points on this processor, including ghosts
        // -------------------------------------------------------------------------

        // Outer loop is local part = find row
        // What is done here, is:
        //
        // 1. Run through all elements in the mesh - should not include ghosts
        for (PetscInt k = info.zs; k < info.zs + info.zm; k++) {
            for (PetscInt j = info.ys; j < info.ys + info.ym; j++) {
                for (PetscInt i = info.xs; i < info.xs + info.xm; i++) {
                    // The row number of the element we are considering:
                    PetscInt row =
                        (i - info.gxs) + (j - info.gys) * (info.gxm) + (k - info.gzs) * (info.gxm) * (info.gym);
                    //
                    // 2. Loop over nodes (including ghosts) within a cubic domain with
                    // center at (i,j,k)
                    //    For each element, run through all elements in a box of size
                    //    stencilWidth * stencilWidth * stencilWidth Remark, we want to
                    //    make sure we are not running "out of the domain", therefore k2
                    //    etc. are limited to the max global index (info.mz-1 etc.)
                    for (PetscInt k2 = PetscMax(k - info.sw, 0); k2 <= PetscMin(k + info.sw, info.mz - 1); k2++) {
                        for (PetscInt j2 = PetscMax(j - info.sw, 0); j2 <= PetscMin(j + info.sw, info.my - 1); j2++) {
                            for (PetscInt i2 = PetscMax(i - info.sw, 0); i2 <= PetscMin(i + info.sw, info.mx - 1);
                                 i2++) {
                                PetscInt col = (i2 - info.gxs) + (j2 - info.gys) * (info.gxm) +
                                               (k2 - info.gzs) * (info.gxm) * (info.gym);
                                PetscScalar dist = 0.0;
                                // Compute the distance from the "col"-element to the
                                // "row"-element
                                for (PetscInt kk = 0; kk < 3; kk++) {
                                    dist = dist + PetscPowScalar(lcoorp[3 * row + kk] - lcoorp[3 * col + kk], 2.0);
                                }
                                dist = PetscSqrtScalar(dist);
                                if (dist < R) {
                                    // Longer distances should have less weight
                                    dist = R - dist;
                                    MatSetValuesLocal(H, 1, &row, 1, &col, &dist, INSERT_VALUES);
                                }
                            }
                        }
                    }
                }
            }
        }
        // Assemble H:
        MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(H, MAT_FINAL_ASSEMBLY);
        // Compute the Hs, i.e. sum the rows
        Vec dummy;
        VecDuplicate(Hs, &dummy);
        VecSet(dummy, 1.0);
        MatMult(H, dummy, Hs);

        // Clean up
        VecRestoreArray(lcoor, &lcoorp);
        VecDestroy(&dummy);
        delete[] Lx;
        delete[] Ly;
        delete[] Lz;

    } else if (filterType == 2) {
        // ALLOCATE AND SETUP THE PDE FILTER CLASS
        pdef = new PDEFilt(da_nodes, R);
    }

    return ierr;
}

PetscErrorCode Filter::DMDAGetElements_3D(DM dm, PetscInt* nel, PetscInt* nen, const PetscInt* e[]) {
    PetscErrorCode ierr;
    DM_DA*         da = (DM_DA*)dm->data;
    PetscInt       i, xs, xe, Xs, Xe;
    PetscInt       j, ys, ye, Ys, Ye;
    PetscInt       k, zs, ze, Zs, Ze;
    PetscInt       cnt = 0, cell[8], ns = 1, nn = 8;
    PetscInt       c;
    if (!da->e) {
        if (da->elementtype == DMDA_ELEMENT_Q1) {
            ns = 1;
            nn = 8;
        }
        ierr = DMDAGetCorners(dm, &xs, &ys, &zs, &xe, &ye, &ze);
        CHKERRQ(ierr);
        ierr = DMDAGetGhostCorners(dm, &Xs, &Ys, &Zs, &Xe, &Ye, &Ze);
        CHKERRQ(ierr);
        xe += xs;
        Xe += Xs;
        if (xs != Xs)
            xs -= 1;
        ye += ys;
        Ye += Ys;
        if (ys != Ys)
            ys -= 1;
        ze += zs;
        Ze += Zs;
        if (zs != Zs)
            zs -= 1;
        da->ne = ns * (xe - xs - 1) * (ye - ys - 1) * (ze - zs - 1);
        PetscMalloc((1 + nn * da->ne) * sizeof(PetscInt), &da->e);
        for (k = zs; k < ze - 1; k++) {
            for (j = ys; j < ye - 1; j++) {
                for (i = xs; i < xe - 1; i++) {
                    cell[0] = (i - Xs) + (j - Ys) * (Xe - Xs) + (k - Zs) * (Xe - Xs) * (Ye - Ys);
                    cell[1] = (i - Xs + 1) + (j - Ys) * (Xe - Xs) + (k - Zs) * (Xe - Xs) * (Ye - Ys);
                    cell[2] = (i - Xs + 1) + (j - Ys + 1) * (Xe - Xs) + (k - Zs) * (Xe - Xs) * (Ye - Ys);
                    cell[3] = (i - Xs) + (j - Ys + 1) * (Xe - Xs) + (k - Zs) * (Xe - Xs) * (Ye - Ys);
                    cell[4] = (i - Xs) + (j - Ys) * (Xe - Xs) + (k - Zs + 1) * (Xe - Xs) * (Ye - Ys);
                    cell[5] = (i - Xs + 1) + (j - Ys) * (Xe - Xs) + (k - Zs + 1) * (Xe - Xs) * (Ye - Ys);
                    cell[6] = (i - Xs + 1) + (j - Ys + 1) * (Xe - Xs) + (k - Zs + 1) * (Xe - Xs) * (Ye - Ys);
                    cell[7] = (i - Xs) + (j - Ys + 1) * (Xe - Xs) + (k - Zs + 1) * (Xe - Xs) * (Ye - Ys);
                    if (da->elementtype == DMDA_ELEMENT_Q1) {
                        for (c = 0; c < ns * nn; c++)
                            da->e[cnt++] = cell[c];
                    }
                }
            }
        }
    }
    *nel = da->ne;
    *nen = nn;
    *e   = da->e;
    return (0);
}
