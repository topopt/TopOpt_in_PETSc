#include <Filter.h>


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


Filter::Filter(TopOpt *opt){
  // Set all pointers to NULL
  H=NULL;  
  Hs=NULL;
  da_elem=NULL;
  pdef=NULL;
  
  // Call the setup method
  SetUp(opt);
}

Filter::~Filter(){
    // Deallocate data
    if (Hs!=NULL){ VecDestroy(&Hs); }
    if (H!=NULL){ MatDestroy(&H); }
    if (da_elem!=NULL){ DMDestroy(&da_elem); }
    if (pdef!=NULL){delete pdef; }
    
}

// Filter design variables
PetscErrorCode Filter::FilterProject(TopOpt *opt){
	PetscErrorCode ierr;
  	
	// Filter the design variables or copy to xPhys
	// STANDARD FILTER
	if (opt->filter == 1){
		// Filter the densitities
		ierr = MatMult(H,opt->x,opt->xPhys); CHKERRQ(ierr);
		VecPointwiseDivide(opt->xPhys,opt->xPhys,Hs);
	}
	// PDE FILTER
	else if (opt->filter == 2 ){
		ierr = pdef->FilterProject(opt->x, opt->xPhys); CHKERRQ(ierr);
		// Check for bound violation: simple, but cheap check!
                PetscScalar *xp;
                PetscInt locsiz;
                VecGetArray(opt->xPhys,&xp);
                VecGetLocalSize(opt->xPhys,&locsiz);
                for (PetscInt i=0;i<locsiz;i++){
                        if (xp[i] < 0.0){
                                if (PetscAbsReal(xp[i]) > 1.0e-4){
                                        PetscPrintf(PETSC_COMM_WORLD,"BOUND VIOLATION IN PDEFILTER - INCREASE RMIN OR MESH RESOLUTION: xPhys = %f\n",xp[i]);
                                }
                                xp[i]= 0.0;
                        }
                        if (xp[i] > 1.0){
                                if (PetscAbsReal(xp[i]-1.0) > 1.0e-4){
                                        PetscPrintf(PETSC_COMM_WORLD,"BOUND VIOLATION IN PDEFILTER - INCREASE RMIN OR MESH RESOLUTION: xPhys = %f\n",xp[i]);
                                }
                                xp[i]=1.0;
                        }

                }
                VecRestoreArray(opt->xPhys,&xp);
	}
	// COPY IN CASE OF SENSITIVITY FILTER
	else {	ierr = VecCopy(opt->x,opt->xPhys); CHKERRQ(ierr); }
		
	return ierr;
}
    
// Filter the sensitivities
PetscErrorCode Filter::Gradients(TopOpt *opt){
	
	PetscErrorCode ierr;
        // Chainrule/Filter for the sensitivities
	if (opt->filter == 0)
	// Filter the sensitivities, df,dg
        {
	        Vec xtmp;
                ierr = VecDuplicate(opt->x,&xtmp);  CHKERRQ(ierr);
	        VecPointwiseMult(xtmp,opt->dfdx,opt->x);
                MatMult(H,xtmp,opt->dfdx);
	        VecPointwiseDivide(xtmp,opt->dfdx,Hs);
	        VecPointwiseDivide(opt->dfdx,xtmp,opt->x);
                VecDestroy(&xtmp);
        }
        else if (opt->filter == 1) {
		// Filter the densities, df,dg: STANDARD FILTER
	        Vec xtmp;
                ierr = VecDuplicate(opt->x,&xtmp);  CHKERRQ(ierr);
		// dfdx
                VecPointwiseDivide(xtmp,opt->dfdx,Hs);
		MatMult(H,xtmp,opt->dfdx);
	        // dgdx
		VecPointwiseDivide(xtmp,opt->dgdx[0],Hs);
		MatMult(H,xtmp,opt->dgdx[0]);
		// tidy up
                VecDestroy(&xtmp);
	}
	else if (opt->filter == 2){
		// Filter the densities, df,dg: PDE FILTER  
		ierr = pdef->Gradients(opt->dfdx,opt->dfdx); CHKERRQ(ierr);
		ierr = pdef->Gradients(opt->dgdx[0],opt->dgdx[0]); CHKERRQ(ierr);
	}
  
	return ierr;
}


PetscErrorCode Filter::SetUp(TopOpt *opt){
  
	PetscErrorCode ierr;
	
  	// Get rank from MPI
	PetscInt myrank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &myrank);
  
	if (opt->filter==0 || opt->filter==1){
	
// Extract information from the nodal mesh
	PetscInt M,N,P,md,nd,pd; 
	DMBoundaryType bx, by, bz;
	DMDAStencilType stype;
	ierr = DMDAGetInfo(opt->da_nodes,NULL,&M,&N,&P,&md,&nd,&pd,NULL,NULL,&bx,&by,&bz,&stype); CHKERRQ(ierr);
	
// Find the element size
	Vec lcoor;
	DMGetCoordinatesLocal(opt->da_nodes,&lcoor);
	PetscScalar *lcoorp;
	VecGetArray(lcoor,&lcoorp);
	
	PetscInt nel, nen;
	const PetscInt *necon;
	DMDAGetElements_3D(opt->da_nodes,&nel,&nen,&necon);
	
	PetscScalar dx,dy,dz;
	// Use the first element to compute the dx, dy, dz
	dx = lcoorp[3*necon[0*nen + 1]+0]-lcoorp[3*necon[0*nen + 0]+0];
	dy = lcoorp[3*necon[0*nen + 2]+1]-lcoorp[3*necon[0*nen + 1]+1];
	dz = lcoorp[3*necon[0*nen + 4]+2]-lcoorp[3*necon[0*nen + 0]+2];
	VecRestoreArray(lcoor,&lcoorp);
	
// Create the minimum element connectivity shit
	PetscInt ElemConn;
	// Check dx,dy,dz and find max conn for a given rmin
	ElemConn = (PetscInt)PetscMax(ceil(opt->rmin/dx)-1,PetscMax(ceil(opt->rmin/dy)-1,ceil(opt->rmin/dz)-1));
	ElemConn = PetscMin(ElemConn,PetscMin((M-1)/2,PetscMin((N-1)/2,(P-1)/2)));
	
	// The following is needed due to roundoff errors 
	PetscInt tmp;
        MPI_Allreduce(&ElemConn, &tmp, 1,MPIU_INT, MPI_MAX,PETSC_COMM_WORLD );
        ElemConn = tmp;
	
	// Print to screen: mesh overlap!
	PetscPrintf(PETSC_COMM_WORLD,"# Filter radius rmin = %f results in a stencil of %i elements \n",opt->rmin,ElemConn);
	
// Find the geometric partitioning of the nodal mesh, so the element mesh will coincide 
	PetscInt *Lx=new PetscInt[md];
	PetscInt *Ly=new PetscInt[nd];
	PetscInt *Lz=new PetscInt[pd];
	
	// get number of nodes for each partition
	const PetscInt *LxCorrect, *LyCorrect, *LzCorrect;
	DMDAGetOwnershipRanges(opt->da_nodes, &LxCorrect, &LyCorrect, &LzCorrect); 

	// subtract one from the lower left corner.
	for (int i=0; i<md; i++){
		Lx[i] = LxCorrect[i];
		if (i==0){Lx[i] = Lx[i]-1;}
	}
	for (int i=0; i<nd; i++){
		Ly[i] = LyCorrect[i];
		if (i==0){Ly[i] = Ly[i]-1;}
	}
	for (int i=0; i<pd; i++){
		Lz[i] = LzCorrect[i];
		if (i==0){Lz[i] = Lz[i]-1;}
	}

	// Create the element grid:
	DMDACreate3d(PETSC_COMM_WORLD,bx,by,bz,stype,M-1,N-1,P-1,md,nd,pd,
			1,ElemConn,Lx,Ly,Lz,&da_elem);
	
	
	// Set the coordinates: from 0+dx/2 to xmax-dx/2 and so on
	PetscScalar xmax = (M-1)*dx;
	PetscScalar ymax = (N-1)*dy;
	PetscScalar zmax = (P-1)*dz;
	DMDASetUniformCoordinates(da_elem , dx/2.0,xmax-dx/2.0, dy/2.0,ymax-dy/2.0, dz/2.0,zmax-dz/2.0);

	// Allocate and assemble
	DMCreateMatrix(da_elem,&H);
	DMCreateGlobalVector(da_elem,&Hs);
	
	//Set the filter matrix and vector
	// 1. assemble H
	DMGetCoordinatesLocal(da_elem,&lcoor);
	VecGetArray(lcoor,&lcoorp);
	DMDALocalInfo info;
	DMDAGetLocalInfo(da_elem,&info);
	// Outer loop is local part = find row
	for (PetscInt k=info.zs; k<info.zs+info.zm; k++) {
		for (PetscInt j=info.ys; j<info.ys+info.ym; j++) {
			for (PetscInt i=info.xs; i<info.xs+info.xm; i++) {
				PetscInt col = (i-info.gxs) + (j-info.gys)*(info.gxm) + (k-info.gzs  )*(info.gxm)*(info.gym);
				// Loop over nodes (including ghosts) within a cubic domain with center at (i,j,k)
				for (PetscInt k2=std::max(k-info.sw,0);k2<=std::min(k+info.sw,info.mz);k2++){
					for (PetscInt j2=std::max(j-info.sw,0);j2<=std::min(j+info.sw,info.my);j2++){
						for (PetscInt i2=std::max(i-info.sw,0);i2<=std::min(i+info.sw,info.mx);i2++){
							PetscInt row = (i2-info.gxs) + (j2-info.gys)*(info.gxm) + (k2-info.gzs  )*(info.gxm)*(info.gym);
							PetscScalar dist = 0.0;
							for(PetscInt kk=0; kk<3; kk++){
								dist = dist + PetscPowScalar(lcoorp[3*col+kk]-lcoorp[3*row+kk],2.0);
							}
							dist = PetscSqrtScalar(dist);
							if (dist<opt->rmin){
								dist = opt->rmin-dist;
								MatSetValuesLocal(H, 1, &col, 1, &row, &dist, INSERT_VALUES); 
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
	// 2. compute the Hs, i.e. sum the rows
	Vec dummy;
	VecDuplicate(Hs,&dummy);
	VecSet(dummy,1.0);
	MatMult(H,dummy,Hs);

	// Clean up
	VecRestoreArray(lcoor,&lcoorp);
	VecDestroy(&dummy);
	delete [] Lx;
	delete [] Ly;
	delete [] Lz;
	
	} 
	else if (opt->filter==2){
		// ALLOCATE AND SETUP THE PDE FILTER CLASS
		pdef = new PDEFilt(opt);
	}
	  
	return ierr;
  
}


PetscErrorCode Filter::DMDAGetElements_3D(DM dm,PetscInt *nel,PetscInt *nen,const PetscInt *e[]) {
	PetscErrorCode ierr;
	DM_DA          *da = (DM_DA*)dm->data;
	PetscInt       i,xs,xe,Xs,Xe;
	PetscInt       j,ys,ye,Ys,Ye;
	PetscInt       k,zs,ze,Zs,Ze;
	PetscInt       cnt=0, cell[8], ns=1, nn=8;
	PetscInt       c; 
	if (!da->e) {
		if (da->elementtype == DMDA_ELEMENT_Q1) {ns=1; nn=8;}
		ierr = DMDAGetCorners(dm,&xs,&ys,&zs,&xe,&ye,&ze);
		CHKERRQ(ierr);
		ierr = DMDAGetGhostCorners(dm,&Xs,&Ys,&Zs,&Xe,&Ye,&Ze);
		CHKERRQ(ierr);
		xe    += xs; Xe += Xs; if (xs != Xs) xs -= 1;
		ye    += ys; Ye += Ys; if (ys != Ys) ys -= 1;
		ze    += zs; Ze += Zs; if (zs != Zs) zs -= 1;
		da->ne = ns*(xe - xs - 1)*(ye - ys - 1)*(ze - zs - 1);
		PetscMalloc((1 + nn*da->ne)*sizeof(PetscInt),&da->e);
		for (k=zs; k<ze-1; k++) {
			for (j=ys; j<ye-1; j++) {
				for (i=xs; i<xe-1; i++) {
					cell[0] = (i-Xs  ) + (j-Ys  )*(Xe-Xs) + (k-Zs  )*(Xe-Xs)*(Ye-Ys);
					cell[1] = (i-Xs+1) + (j-Ys  )*(Xe-Xs) + (k-Zs  )*(Xe-Xs)*(Ye-Ys);
					cell[2] = (i-Xs+1) + (j-Ys+1)*(Xe-Xs) + (k-Zs  )*(Xe-Xs)*(Ye-Ys);
					cell[3] = (i-Xs  ) + (j-Ys+1)*(Xe-Xs) + (k-Zs  )*(Xe-Xs)*(Ye-Ys);
					cell[4] = (i-Xs  ) + (j-Ys  )*(Xe-Xs) + (k-Zs+1)*(Xe-Xs)*(Ye-Ys);
					cell[5] = (i-Xs+1) + (j-Ys  )*(Xe-Xs) + (k-Zs+1)*(Xe-Xs)*(Ye-Ys);
					cell[6] = (i-Xs+1) + (j-Ys+1)*(Xe-Xs) + (k-Zs+1)*(Xe-Xs)*(Ye-Ys);
					cell[7] = (i-Xs  ) + (j-Ys+1)*(Xe-Xs) + (k-Zs+1)*(Xe-Xs)*(Ye-Ys);
					if (da->elementtype == DMDA_ELEMENT_Q1) {
						for (c=0; c<ns*nn; c++) da->e[cnt++] = cell[c];
					}
				}
			}
		}
	}
	*nel = da->ne;
	*nen = nn;
	*e   = da->e;
	return(0);
}
