#include <PDEFilter.h>
#include <TopOpt.h>
//#include <petsc-private/dmdaimpl.h>
#include <petsc/private/dmdaimpl.h>
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


PDEFilt::PDEFilt(DM da_nodes, PetscScalar rmin)
{

	R=rmin/2.0/sqrt(3); // conversion factor for the PDEfilter

	nlvls=3; // MG levels

	// number of nodal dofs
	PetscInt numnodaldof = 1;

	// Stencil width: each node connects to a box around it - linear elements
	PetscInt stencilwidth = 1;

	PetscScalar dx,dy,dz;
	DMBoundaryType bx, by, bz;
	DMDAStencilType stype;
	{
		// Extract information from the nodal mesh
		PetscInt M,N,P,md,nd,pd; 
		DMDAGetInfo(da_nodes,NULL,&M,&N,&P,&md,&nd,&pd,NULL,NULL,&bx,&by,&bz,&stype); 

		// Find the element size
		Vec lcoor;
		DMGetCoordinatesLocal(da_nodes,&lcoor);
		PetscScalar *lcoorp;
		VecGetArray(lcoor,&lcoorp);

		PetscInt nel, nen;
		const PetscInt *necon;
		DMDAGetElements_3D(da_nodes,&nel,&nen,&necon);

		// Use the first element to compute the dx, dy, dz
		dx = lcoorp[3*necon[0*nen + 1]+0]-lcoorp[3*necon[0*nen + 0]+0];
		dy = lcoorp[3*necon[0*nen + 2]+1]-lcoorp[3*necon[0*nen + 1]+1];
		dz = lcoorp[3*necon[0*nen + 4]+2]-lcoorp[3*necon[0*nen + 0]+2];
		VecRestoreArray(lcoor,&lcoorp);

		// ELement volume
		elemVol = dx*dy*dz;

		nn[0]=M;
		nn[1]=N;
		nn[2]=P;

		ne[0]=nn[0]-1; 
		ne[1]=nn[1]-1; 
		ne[2]=nn[2]-1; 


		xc[0]=0.0;
		xc[1]=ne[0]*M;
		xc[2]=0.0;
		xc[3]=ne[1]*N;
		xc[4]=0.0;
		xc[5]=ne[2]*P;

	}

	// Create the nodal mesh
	DMDACreate3d(PETSC_COMM_WORLD,bx,by,bz,stype,nn[0],nn[1],nn[2],PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
			numnodaldof,stencilwidth,0,0,0,&(da_nodal));
	// Initialize
  	DMSetFromOptions(da_nodal);
  	DMSetUp(da_nodal);

	
	// Set the coordinates
	DMDASetUniformCoordinates(da_nodal, xc[0],xc[1], xc[2],xc[3], xc[4],xc[5]);
	// Set the element type to Q1: Otherwise calls to GetElements will change to P1 !
	// STILL DOESN*T WORK !!!!
	DMDASetElementType(da_nodal, DMDA_ELEMENT_Q1);

	//Create the element mesh

	// find the geometric partitioning of the nodal mesh, so the element mesh will coincide
	PetscInt md,nd,pd;
	DMDAGetInfo(da_nodal,NULL,NULL,NULL,NULL,&md,&nd,&pd,NULL,NULL,NULL,NULL,NULL,NULL);
	PetscInt *Lx=new PetscInt[md];
	PetscInt *Ly=new PetscInt[nd];
	PetscInt *Lz=new PetscInt[pd];
	// get number of nodes for each partition
	const PetscInt *LxCorrect, *LyCorrect, *LzCorrect;
	DMDAGetOwnershipRanges(da_nodal, &LxCorrect, &LyCorrect, &LzCorrect);
	// subtract one from the lower left corner
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

	PetscInt overlap=0; 
	// Create the element grid:
	DMDACreate3d(PETSC_COMM_WORLD,bx,by,bz,stype,nn[0]-1,nn[1]-1,nn[2]-1,md,nd,pd,
			1,overlap,Lx,Ly,Lz,&(da_element));
	// Initialize
	DMSetFromOptions(da_element);
        DMSetUp(da_element);


	delete [] Lx;
	delete [] Ly;
	delete [] Lz;

	PDEFilterMatrix(dx,dy,dz, R, KF, TF);

	//create the stiffness matrix
	DMCreateMatrix(da_nodal,&(K));
	//create RHS
	DMCreateGlobalVector(da_nodal,&(RHS));
	DMCreateGlobalVector(da_element,&(X));
	VecDuplicate(RHS, &U);


	//Create T matrix
	{
		PetscInt m;
		PetscInt n;
		//PetscInt M;
		//PetscInt N;

		//m,M  extract it from RHS
		//n,N  extract it from X
		VecGetLocalSize(RHS,&m);
		VecGetLocalSize(X,&n);

		MatCreateAIJ(PETSC_COMM_WORLD , m, n, PETSC_DETERMINE, PETSC_DETERMINE, 8, NULL,7,NULL, &T);

		ISLocalToGlobalMapping rmapping;
		ISLocalToGlobalMapping cmapping;	

		DMGetLocalToGlobalMapping(da_nodal, &rmapping);	
		DMGetLocalToGlobalMapping(da_element, &cmapping);


		MatSetLocalToGlobalMapping(T,rmapping,cmapping);

	}


	MatAssemble();
	SetUpSolver();

	//test
	PetscRandom rctx;
	PetscRandomCreate(PETSC_COMM_WORLD,&rctx);
	PetscRandomSetType(rctx,PETSCRAND48);
	VecSetRandom(X,rctx);
	PetscRandomDestroy(&rctx);

	FilterProject(X,X);
	Gradients(X,X);

	//
	PetscPrintf(PETSC_COMM_WORLD,"Done setting up the PDEFilter\n");
}

PetscErrorCode PDEFilt::FilterProject(Vec OX, Vec FX)
{

	PetscErrorCode ierr;

	double t1,t2;
	PetscScalar rnorm;
	PetscInt niter;

	t1 = MPI_Wtime();
	ierr = MatMult(T,OX,RHS); CHKERRQ(ierr);
	ierr = VecCopy(RHS,U); CHKERRQ(ierr);
	ierr = VecScale(RHS,elemVol);CHKERRQ(ierr);
	ierr = KSPSolve(ksp,RHS,U); CHKERRQ(ierr);
	ierr = KSPGetIterationNumber(ksp,&niter); CHKERRQ(ierr);
	ierr = KSPGetResidualNorm(ksp,&rnorm);  CHKERRQ(ierr);
	ierr = MatMultTranspose(T,U,FX); CHKERRQ(ierr);

	t2 = MPI_Wtime();
	PetscPrintf(PETSC_COMM_WORLD,"PDEFilter solver:  iter: %i, rerr.: %e, time: %f\n",niter,rnorm,t2-t1);
	return ierr;
}

PetscErrorCode PDEFilt::Gradients(Vec OS, Vec FS)
{
	return FilterProject(OS,FS);
}

PDEFilt::~PDEFilt()
{
	Free();
}

PetscErrorCode PDEFilt::Free()
{

	PetscErrorCode ierr;

	KSPDestroy(&ksp);

	VecDestroy(&RHS);
	VecDestroy(&X);
	VecDestroy(&U);

	MatDestroy(&T);
	MatDestroy(&K);

	ierr=DMDestroy(&da_nodal);    CHKERRQ(ierr);
	ierr=DMDestroy(&da_element);  CHKERRQ(ierr);

	return ierr;

}

void PDEFilt::MatAssemble()
{
	// Get the FE mesh structure (from the nodal mesh)
	PetscInt nel, nen;
	const PetscInt *necon;
	DMDAGetElements_3D(da_nodal,&nel,&nen,&necon);
	MatZeroEntries(K);
	MatZeroEntries(T);
	PetscInt *edof = new PetscInt[8];
	for (PetscInt i=0;i<nel;i++)
	{
		// loop over element nodes
		for (PetscInt j=0;j<nen;j++)
		{
			edof[j]=necon[i*nen+j];
		}

		MatSetValuesLocal(K,8,edof,8,edof,KF,ADD_VALUES);
		//assemble the T matrix 
		MatSetValuesLocal(T,8,edof,1,&i,TF,ADD_VALUES);	


	}
	MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(T, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(T, MAT_FINAL_ASSEMBLY);

	delete [] edof;

}


PetscErrorCode PDEFilt::SetUpSolver()
{
	//make sure ksp is not allocated before 
	PetscErrorCode ierr;
	PC pc;

	// The fine grid Krylov method
	KSPCreate(PETSC_COMM_WORLD,&ksp);
	ierr = KSPSetType(ksp,KSPFGMRES); // KSPCG, KSPGMRES
	PetscInt restart = 20;
	ierr = KSPGMRESSetRestart(ksp,restart);

	PetscScalar rtol = 1.0e-8;
	PetscScalar atol = 1.0e-50;
	PetscScalar dtol = 1.0e3;
	PetscInt maxitsGlobal = 60;
	ierr = KSPSetTolerances(ksp,rtol,atol,dtol,maxitsGlobal);
	ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
	KSPSetOperators(ksp,K,K); // ,SAME_PRECONDITIONER is now set in the prec

	//preconditioner
	KSPGetPC(ksp,&pc);
	PCSetType(pc,PCMG);
	// Set solver from options
	KSPSetFromOptions(ksp);
	// Get the prec again - check if it has changed
	KSPGetPC(ksp,&pc);
	ierr = PCSetReusePreconditioner(pc,PETSC_TRUE); CHKERRQ(ierr);
	// Flag for pcmg pc
	PetscBool pcmg_flag = PETSC_TRUE;
	PetscObjectTypeCompare((PetscObject)pc,PCMG,&pcmg_flag);
	// Only if PCMG is used
	if (pcmg_flag){
		// DMs for grid hierachy
		DM  *da_list,*daclist;
		Mat R;
		PetscMalloc(sizeof(DM)*nlvls,&da_list);
		for (PetscInt k=0; k<nlvls; k++) da_list[k] = NULL;
		PetscMalloc(sizeof(DM)*nlvls,&daclist);
		for (PetscInt k=0; k<nlvls; k++) daclist[k] = NULL;
		// Set 0 to the finest level
		daclist[0] = da_nodal;

		// Coordinates
		PetscReal xmin=xc[0], xmax=xc[1], ymin=xc[2], ymax=xc[3], zmin=xc[4], zmax=xc[5];

		// Set up the coarse meshes
		ierr=DMCoarsenHierarchy(da_nodal, nlvls-1,&daclist[1]); CHKERRQ(ierr);
		for (PetscInt k=0; k<nlvls; k++) {
			// NOTE: finest grid is nlevels - 1: PCMG MUST USE THIS ORDER ???
			da_list[k] = daclist[nlvls-1-k];
			DMDASetUniformCoordinates(da_list[k],xmin,xmax,ymin,ymax,zmin,zmax);
		}

		// the PCMG specific options
		PCMGSetLevels(pc,nlvls,NULL);
		PCMGSetType(pc,PC_MG_MULTIPLICATIVE); // Default
		PCMGSetCycleType(pc,PC_MG_CYCLE_V);
		PCMGSetGalerkin(pc,PC_MG_GALERKIN_BOTH);
		for (PetscInt k=1; k<nlvls; k++) {
			DMCreateInterpolation(da_list[k-1],da_list[k],&R,NULL);
			PCMGSetInterpolation(pc,k,R);
			MatDestroy(&R);
		}

		for (PetscInt k=1; k<nlvls; k++) { //0 level should be dealocated in the destructor
			DMDestroy(&daclist[k]);
		}
		PetscFree(da_list);
		PetscFree(daclist);

		// AVOID THE DEFAULT FOR THE MG PART
		{
			// SET the coarse grid solver:
			// i.e. get a pointer to the ksp and change its settings
			KSP cksp;
			PCMGGetCoarseSolve(pc,&cksp);
			// The solver
			ierr = KSPSetType(cksp,KSPGMRES); // KSPCG, KSPFGMRES
			// PetscInt restarts[nlvls] = {10, 1 , 1}; // coarse .... fine
			restart = 10;
			ierr = KSPGMRESSetRestart(cksp,restart);
			rtol = 1.0e-8;
			atol = 1.0e-50;
			dtol = 1e3;
			PetscInt maxits = 10;
			ierr = KSPSetTolerances(cksp,rtol,atol,dtol,maxits);
			// The preconditioner
			PC cpc;
			KSPGetPC(cksp,&cpc);
			//PCSetType(cpc,PCSOR); // PCSOR, PCSPAI (NEEDS TO BE COMPILED), PCJACOBI
			PCSetType(cpc,PCJACOBI);

			// Set smoothers on all levels (except for coarse grid):
			for (PetscInt k=1;k<nlvls;k++){
				KSP dksp;
				PCMGGetSmoother(pc,k,&dksp);
				PC dpc;
				KSPGetPC(dksp,&dpc);
				ierr = KSPSetType(dksp,KSPGMRES); // KSPCG, KSPGMRES, KSPCHEBYSHEV (VERY GOOD FOR SPD)
				restart = 1;
				ierr = KSPGMRESSetRestart(dksp,restart);
				ierr = KSPSetTolerances(dksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,restart); // NOTE maxitr=restart;
				PCSetType(dpc,PCJACOBI);// PCJACOBI, PCSOR for KSPCHEBYSHEV very good
			}
		}


	}


	// 	// Write check to screen:
	//         // Check the overall Krylov solver
	//         KSPType ksptype;
	//         KSPGetType(ksp,&ksptype);
	//         PCType pctype;
	//         PCGetType(pc,&pctype);
	//         PetscInt mmax;
	//         KSPGetTolerances(ksp,NULL,NULL,NULL,&mmax);
	//         PetscPrintf(PETSC_COMM_WORLD,"##############################################################\n");
	//         PetscPrintf(PETSC_COMM_WORLD,"################# Linear solver settings #####################\n");
	//         PetscPrintf(PETSC_COMM_WORLD,"# Main solver: %s, prec.: %s, maxiter.: %i \n",ksptype,pctype,mmax);
	// 
	//         // Only if pcmg is used
	//         if (pcmg_flag){
	//                 // Check the smoothers and coarse grid solver:
	//                 for (PetscInt k=0;k<nlvls;k++){
	//                         KSP dksp;
	//                         PC dpc;
	//                         KSPType dksptype;
	//                         PCMGGetSmoother(pc,k,&dksp);
	//                         KSPGetType(dksp,&dksptype);
	//                         KSPGetPC(dksp,&dpc);
	//                         PCType dpctype;
	//                         PCGetType(dpc,&dpctype);
	//                         PetscInt mmax;
	//                         KSPGetTolerances(dksp,NULL,NULL,NULL,&mmax);
	//                         PetscPrintf(PETSC_COMM_WORLD,"# Level %i smoother: %s, prec.: %s, sweep: %i \n",k,dksptype,dpctype,mmax);
	//                 }
	//         }
	//         PetscPrintf(PETSC_COMM_WORLD,"##############################################################\n");


	return 0;
}


PetscErrorCode PDEFilt::DMDAGetElements_3D(DM dm,PetscInt *nel,PetscInt *nen,const PetscInt *e[]) 
{
	DM_DA          *da = (DM_DA*)dm->data;
	PetscInt       i,xs,xe,Xs,Xe;
	PetscInt       j,ys,ye,Ys,Ye;
	PetscInt       k,zs,ze,Zs,Ze;
	PetscInt       cnt=0, cell[8], ns=1, nn=8;
	PetscInt       c;
	if (!da->e) {
		if (da->elementtype == DMDA_ELEMENT_Q1) {ns=1; nn=8;}
		DMDAGetCorners(dm,&xs,&ys,&zs,&xe,&ye,&ze);
		DMDAGetGhostCorners(dm,&Xs,&Ys,&Zs,&Xe,&Ye,&Ze);
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



void PDEFilt::PDEFilterMatrix(PetscScalar dx, PetscScalar dy, PetscScalar dz,
		PetscScalar RR,
		PetscScalar *KK, PetscScalar *T)
{
	PetscScalar t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t15,t16,t18,t22,t23,t27,t28,t32,t36,t37,t41,t45,t49,t53;
	t3 = 1.0/dx/dy;
	t4 = 1/dz;
	t5 = RR*RR;
	t6 = dx*dx;
	t7 = t5*t6;
	t8 = dy*dy;
	t9 = t7*t8;
	t10 = 3.0*t9;
	t11 = dz*dz;
	t12 = t7*t11;
	t13 = 3.0*t12;
	t15 = t5*t8*t11;
	t16 = 3.0*t15;
	t18 = t6*t8*t11;
	t22 = t3*t4*(t10+t13+t16+t18)/27.0;
	t23 = 6.0*t15;
	t27 = t3*t4*(t10+t13-t23+t18)/54.0;
	t28 = 6.0*t12;
	t32 = t3*t4*(t10-t28-t23+t18)/108.0;
	t36 = t3*t4*(t10-t28+t16+t18)/54.0;
	t37 = 6.0*t9;
	t41 = t3*t4*(t37-t13-t16-t18)/54.0;
	t45 = t3*t4*(t37-t13+t23-t18)/108.0;
	t49 = t3*t4*(t37+t28+t23-t18)/216.0;
	t53 = t3*t4*(t37+t28-t16-t18)/108.0;

	KK[0] = t22;
	KK[1] = t27;
	KK[2] = t32;
	KK[3] = t36;
	KK[4] = -t41;
	KK[5] = -t45;
	KK[6] = -t49;
	KK[7] = -t53;
	KK[8] = t27;
	KK[9] = t22;
	KK[10] = t36;
	KK[11] = t32;
	KK[12] = -t45;
	KK[13] = -t41;
	KK[14] = -t53;
	KK[15] = -t49;
	KK[16] = t32;
	KK[17] = t36;
	KK[18] = t22;
	KK[19] = t27;
	KK[20] = -t49;
	KK[21] = -t53;
	KK[22] = -t41;
	KK[23] = -t45;
	KK[24] = t36;
	KK[25] = t32;
	KK[26] = t27;
	KK[27] = t22;
	KK[28] = -t53;
	KK[29] = -t49;
	KK[30] = -t45;
	KK[31] = -t41;
	KK[32] = -t41;
	KK[33] = -t45;
	KK[34] = -t49;
	KK[35] = -t53;
	KK[36] = t22;
	KK[37] = t27;
	KK[38] = t32;
	KK[39] = t36;
	KK[40] = -t45;
	KK[41] = -t41;
	KK[42] = -t53;
	KK[43] = -t49;
	KK[44] = t27;
	KK[45] = t22;
	KK[46] = t36;
	KK[47] = t32;
	KK[48] = -t49;
	KK[49] = -t53;
	KK[50] = -t41;
	KK[51] = -t45;
	KK[52] = t32;
	KK[53] = t36;
	KK[54] = t22;
	KK[55] = t27;
	KK[56] = -t53;
	KK[57] = -t49;
	KK[58] = -t45;
	KK[59] = -t41;
	KK[60] = t36;
	KK[61] = t32;
	KK[62] = t27;
	KK[63] = t22;

	PetscScalar vol=1.0;
	T[0]=0.125*vol;
	T[1]=0.125*vol;
	T[2]=0.125*vol;
	T[3]=0.125*vol;
	T[4]=0.125*vol;
	T[5]=0.125*vol;
	T[6]=0.125*vol;
	T[7]=0.125*vol;

}
