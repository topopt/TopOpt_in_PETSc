#include <LinearElasticity.h>

/*
 Authors: Niels Aage, Erik Andreassen, Boyan Lazarov, August 2013

 Disclaimer:                                                              
 The authors reserves all rights but does not guaranty that the code is   
 free from errors. Furthermore, we shall not be liable in any event     
 caused by the use of the program.                                     
*/


LinearElasticity::LinearElasticity(TopOpt *opt){
	// Set pointers to null
	K = NULL;
	U = NULL;
	RHS = NULL;
	N = NULL;
	ksp = NULL;

	// Setup sitffness matrix, load vector and bcs (Dirichlet) for the design problem
	SetUpLoadAndBC(opt);

}

LinearElasticity::~LinearElasticity(){
	// Deallocate
	VecDestroy(&(U));
	VecDestroy(&(RHS));
	VecDestroy(&(N));
	MatDestroy(&(K));
	KSPDestroy(&(ksp));

}

PetscErrorCode LinearElasticity::SetUpLoadAndBC(TopOpt *opt){

	PetscErrorCode ierr;

	// Allocate matrix and the RHS and Solution vector and Dirichlet vector
	ierr = DMCreateMatrix(opt->da_nodes,&(K)); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(opt->da_nodes,&(U)); CHKERRQ(ierr);
	VecDuplicate(U,&(RHS));
	VecDuplicate(U,&(N));

	// Set the local stiffness matrix
	PetscScalar a = opt->dx; // x-side length
	PetscScalar b = opt->dy; // y-side length
	PetscScalar c = opt->dz; // z-side length 
	PetscScalar X[8] = {0.0, a, a, 0.0, 0.0, a, a, 0.0};
	PetscScalar Y[8] = {0.0, 0.0, b, b, 0.0, 0.0, b, b};
	PetscScalar Z[8] = {0.0, 0.0, 0.0, 0.0, c, c, c, c};

	// Compute the element stiffnes matrix - constant due to structured grid
	Hex8Isoparametric(X, Y, Z, 0.3, false, KE); 

	// Set the RHS and Dirichlet vector
	VecSet(N,1.0);
	VecSet(RHS,0.0);

	// Global coordinates and a pointer
	Vec lcoor; // borrowed ref - do not destroy!
	PetscScalar *lcoorp;

	// Get local coordinates in local node numbering including ghosts       
	ierr = DMGetCoordinatesLocal(opt->da_nodes,&lcoor); CHKERRQ(ierr);
	VecGetArray(lcoor,&lcoorp);

	// Get local dof number
	PetscInt nn;
	VecGetSize(lcoor,&nn); 
	
	// Compute epsilon parameter for finding points in space:
	PetscScalar epsi = PetscMin(a*0.05,PetscMin(b*0.05,c*0.05));

	// Set the values:
	// In this case: N = the wall at x=xmin is fully clamped
	//               RHS(z) = sin(pi*y/Ly) at x=xmax,z=zmin;
	// OR
	//               RHS(z) = -0.1 at x=xmax,z=zmin;
	for (PetscInt i=0;i<nn;i++){
		// Make a wall with all dofs clamped
		if (i % 3 == 0 && PetscAbsScalar(lcoorp[i]-opt->xc[0]) < epsi){
			VecSetValueLocal(N,i,0.0,INSERT_VALUES);
			VecSetValueLocal(N,++i,0.0,INSERT_VALUES);
			VecSetValueLocal(N,++i,0.0,INSERT_VALUES);
		}
		// Line load
		if (i % 3 == 0 && PetscAbsScalar(lcoorp[i]-opt->xc[1]) < epsi && 
				  PetscAbsScalar(lcoorp[i+2]-opt->xc[4]) < epsi){
			VecSetValueLocal(RHS,i+2,-0.1,INSERT_VALUES);
		}
		// Adjust the corners
		if (i % 3 == 0 && PetscAbsScalar(lcoorp[i]-opt->xc[1]) < epsi && 
				  PetscAbsScalar(lcoorp[i+1]-opt->xc[2]) < epsi && 
				  PetscAbsScalar(lcoorp[i+2]-opt->xc[4]) < epsi ){
			VecSetValueLocal(RHS,i+2,-0.05,INSERT_VALUES);
		}
		if (i % 3 == 0 && PetscAbsScalar(lcoorp[i]-opt->xc[1]) < epsi && 
				  PetscAbsScalar(lcoorp[i+1]-opt->xc[3]) < epsi && 
   				  PetscAbsScalar(lcoorp[i+2]-opt->xc[4]) < epsi){
			VecSetValueLocal(RHS,i+2,-0.05,INSERT_VALUES);
		}
	}

	VecAssemblyBegin(N);
	VecAssemblyBegin(RHS);
	VecAssemblyEnd(N);
	VecAssemblyEnd(RHS);
	VecRestoreArray(lcoor,&lcoorp);

	return ierr;

}

PetscErrorCode LinearElasticity::SolveState(TopOpt *opt){

	PetscErrorCode ierr;

	double t1,t2;
	t1 = MPI_Wtime();

	// Assemble the stiffness matrix
	ierr = AssembleStiffnessMatrix(opt);
	CHKERRQ(ierr);

	// Setup the solver
	if (ksp==NULL){
		ierr = SetUpSolver(opt);
		CHKERRQ(ierr);
	}
	else {
		ierr = KSPSetOperators(ksp,K,K);
		CHKERRQ(ierr);
		KSPSetUp(ksp); 
	}


	// Solve
	ierr = KSPSolve(ksp,RHS,U); CHKERRQ(ierr);
	CHKERRQ(ierr);

	// DEBUG
	// Get iteration number and residual from KSP
	PetscInt niter;
	PetscScalar rnorm;
	KSPGetIterationNumber(ksp,&niter);
	KSPGetResidualNorm(ksp,&rnorm); 

	t2 = MPI_Wtime();

	PetscPrintf(PETSC_COMM_WORLD,"State solver:  iter: %i, rerr.: %e, time: %f\n",niter,rnorm,t2-t1);

	return ierr;
}

PetscErrorCode LinearElasticity::ComputeObjectiveConstraints(TopOpt *opt) {

	// Error code
	PetscErrorCode ierr;

	// Solve state eqs 
	ierr = SolveState(opt); CHKERRQ(ierr); 

	// Get the FE mesh structure (from the nodal mesh)
	PetscInt nel, nen;
	const PetscInt *necon;
	ierr = DMDAGetElements_3D(opt->da_nodes,&nel,&nen,&necon); CHKERRQ(ierr);

	// Get pointer to the densities
	PetscScalar *xp;
	VecGetArray(opt->xPhys,&xp);

	// Get Solution
	Vec Uloc;
	DMCreateLocalVector(opt->da_nodes,&Uloc);
	DMGlobalToLocalBegin(opt->da_nodes,U,INSERT_VALUES,Uloc);
	DMGlobalToLocalEnd(opt->da_nodes,U,INSERT_VALUES,Uloc);

	// get pointer to local vector
	PetscScalar *up;
	VecGetArray(Uloc,&up);

	// Edof array
	PetscInt *edof = new PetscInt[24];

	opt->fx = 0.0;
	// Loop over elements
	for (PetscInt i=0;i<nel;i++){
		// loop over element nodes
		for (PetscInt j=0;j<nen;j++){
			// Get local dofs
			for (PetscInt k=0;k<3;k++){
				edof[j*3+k] = 3*necon[i*nen+j]+k;
			}
		}
		// Use SIMP for stiffness interpolation
		PetscScalar uKu=0.0;
		for (PetscInt k=0;k<24;k++){
			for (PetscInt h=0;h<24;h++){	
				uKu += up[edof[k]]*KE[k*24+h]*up[edof[h]];
			}
		}
		// Add to objective
		opt->fx += (opt->Emin + PetscPowScalar(xp[i],opt->penal)*(opt->Emax - opt->Emin))*uKu;
	}

	// Allreduce fx[0]
	PetscScalar tmp=opt->fx;
	opt->fx=0.0;
	MPI_Allreduce(&tmp,&(opt->fx),1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);		

	// Compute volume constraint gx[0]
	opt->gx[0]=0;
	VecSum(opt->xPhys, &(opt->gx[0]));
	opt->gx[0]=opt->gx[0]/(((PetscScalar)opt->n)*opt->volfrac)-1.0;


	VecRestoreArray(opt->xPhys,&xp);
	VecRestoreArray(Uloc,&up);
	VecDestroy(&Uloc);

	return(ierr);  

}

PetscErrorCode LinearElasticity::ComputeSensitivities(TopOpt *opt) {

	PetscErrorCode ierr;

	// Get the FE mesh structure (from the nodal mesh)
	PetscInt nel, nen;
	const PetscInt *necon;
	ierr = DMDAGetElements_3D(opt->da_nodes,&nel,&nen,&necon); CHKERRQ(ierr);
	
	// Get pointer to the densities
	PetscScalar *xp;
	VecGetArray(opt->xPhys,&xp);

	// Get Solution
	Vec Uloc;
	DMCreateLocalVector(opt->da_nodes,&Uloc);
	DMGlobalToLocalBegin(opt->da_nodes,U,INSERT_VALUES,Uloc);
	DMGlobalToLocalEnd(opt->da_nodes,U,INSERT_VALUES,Uloc);

	// get pointer to local vector
	PetscScalar *up;
	VecGetArray(Uloc,&up);

	// Get dfdx
	PetscScalar *df;
	VecGetArray(opt->dfdx,&df);

	// Edof array
	PetscInt *edof = new PetscInt[24];

	// Loop over elements
	for (PetscInt i=0;i<nel;i++){
		// loop over element nodes
		for (PetscInt j=0;j<nen;j++){
			// Get local dofs
			for (PetscInt k=0;k<3;k++){
				edof[j*3+k] = 3*necon[i*nen+j]+k;
			}
		}
		// Use SIMP for stiffness interpolation
		PetscScalar uKu=0.0;
		for (PetscInt k=0;k<24;k++){
			for (PetscInt h=0;h<24;h++){	
				uKu += up[edof[k]]*KE[k*24+h]*up[edof[h]];
			}
		}
		// Set the Senstivity
		df[i]= -1.0 * opt->penal*PetscPowScalar(xp[i],opt->penal-1)*(opt->Emax - opt->Emin)*uKu;
	}
	// Compute volume constraint gx[0]
	VecSet(opt->dgdx[0],1.0/(((PetscScalar)opt->n)*opt->volfrac));

	VecRestoreArray(opt->xPhys,&xp);
	VecRestoreArray(Uloc,&up);
	VecRestoreArray(opt->dfdx,&df);
	VecDestroy(&Uloc);

	return(ierr);  

}

PetscErrorCode LinearElasticity::ComputeObjectiveConstraintsSensitivities(TopOpt *opt) {
	// Errorcode
	PetscErrorCode ierr;

	// Solve state eqs 
	ierr = SolveState(opt); CHKERRQ(ierr); 

	// Get the FE mesh structure (from the nodal mesh)
	PetscInt nel, nen;
	const PetscInt *necon;
	ierr = DMDAGetElements_3D(opt->da_nodes,&nel,&nen,&necon); CHKERRQ(ierr);
	//DMDAGetElements(da_nodes,&nel,&nen,&necon); // Still issue with elemtype change !

	// Get pointer to the densities
	PetscScalar *xp;
	VecGetArray(opt->xPhys,&xp);

	// Get Solution
	Vec Uloc;
	DMCreateLocalVector(opt->da_nodes,&Uloc);
	DMGlobalToLocalBegin(opt->da_nodes,U,INSERT_VALUES,Uloc);
	DMGlobalToLocalEnd(opt->da_nodes,U,INSERT_VALUES,Uloc);

	// get pointer to local vector
	PetscScalar *up;
	VecGetArray(Uloc,&up);

	// Get dfdx
	PetscScalar *df;
	VecGetArray(opt->dfdx,&df);

	// Edof array
	PetscInt *edof = new PetscInt[24];

	opt->fx = 0.0;
	// Loop over elements
	for (PetscInt i=0;i<nel;i++){
		// loop over element nodes
		for (PetscInt j=0;j<nen;j++){
			// Get local dofs
			for (PetscInt k=0;k<3;k++){
				edof[j*3+k] = 3*necon[i*nen+j]+k;
			}
		}
		// Use SIMP for stiffness interpolation
		PetscScalar uKu=0.0;
		for (PetscInt k=0;k<24;k++){
			for (PetscInt h=0;h<24;h++){	
				uKu += up[edof[k]]*KE[k*24+h]*up[edof[h]];
			}
		}
		// Add to objective
		opt->fx += (opt->Emin + PetscPowScalar(xp[i],opt->penal)*(opt->Emax - opt->Emin))*uKu;
		// Set the Senstivity
		df[i]= -1.0 * opt->penal*PetscPowScalar(xp[i],opt->penal-1)*(opt->Emax - opt->Emin)*uKu;
	}

	// Allreduce fx[0]
	PetscScalar tmp=opt->fx;
	opt->fx=0.0;
	MPI_Allreduce(&tmp,&(opt->fx),1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD);		

	// Compute volume constraint gx[0]
	opt->gx[0]=0;
	VecSum(opt->xPhys, &(opt->gx[0]));
	opt->gx[0]=opt->gx[0]/(((PetscScalar)opt->n)*opt->volfrac)-1.0;
	VecSet(opt->dgdx[0],1.0/(((PetscScalar)opt->n)*opt->volfrac));

	VecRestoreArray(opt->xPhys,&xp);
	VecRestoreArray(Uloc,&up);
	VecRestoreArray(opt->dfdx,&df);
	VecDestroy(&Uloc);

	return(ierr);  

}

//##################################################################
//##################################################################
//##################################################################
// ######################## PRIVATE ################################
//##################################################################
//##################################################################

PetscErrorCode LinearElasticity::AssembleStiffnessMatrix(TopOpt *opt){

	PetscErrorCode ierr;

	// Get the FE mesh structure (from the nodal mesh)
	PetscInt nel, nen;
	const PetscInt *necon;
	ierr = DMDAGetElements_3D(opt->da_nodes,&nel,&nen,&necon);
	CHKERRQ(ierr);

	// Get pointer to the densities
	PetscScalar *xp;
	VecGetArray(opt->xPhys,&xp);

	// Zero the matrix
	MatZeroEntries(K);	

	// Edof array
	PetscInt *edof = new PetscInt[24];
	PetscScalar ke[24*24];

	// Loop over elements
	for (PetscInt i=0;i<nel;i++){
		// loop over element nodes
		for (PetscInt j=0;j<nen;j++){
			// Get local dofs
			for (PetscInt k=0;k<3;k++){
				edof[j*3+k] = 3*necon[i*nen+j]+k;
			}
		} 
		// Use SIMP for stiffness interpolation
		PetscScalar dens = opt->Emin + PetscPowScalar(xp[i],opt->penal)*(opt->Emax-opt->Emin);
		for (PetscInt k=0;k<24*24;k++){
			ke[k]=KE[k]*dens;
		}
		// Add values to the sparse matrix
		ierr = MatSetValuesLocal(K,24,edof,24,edof,ke,ADD_VALUES);
		CHKERRQ(ierr);
	}
	MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);

	// Impose the dirichlet conditions, i.e. K = N'*K*N - (N-I)
	// 1.: K = N'*K*N
	MatDiagonalScale(K,N,N);
	// 2. Add ones, i.e. K = K + NI, NI = I - N
	Vec NI;
	VecDuplicate(N,&NI);
	VecSet(NI,1.0);
	VecAXPY(NI,-1.0,N);
	MatDiagonalSet(K,NI,ADD_VALUES);

	// Zero out possible loads in the RHS that coincide
	// with Dirichlet conditions
	VecPointwiseMult(RHS,RHS,N);

	delete [] edof;
	VecDestroy(&NI);
	VecRestoreArray(opt->xPhys,&xp);
	DMDARestoreElements(opt->da_nodes,&nel,&nen,&necon);

	return ierr;
}

PetscErrorCode LinearElasticity::SetUpSolver(TopOpt *opt){

	PetscErrorCode ierr;
	PC pc;

	// The fine grid Krylov method
	KSPCreate(PETSC_COMM_WORLD,&(ksp));

	// SET THE DEFAULT SOLVER PARAMETERS
	// The fine grid solver settings
	PetscScalar rtol = 1.0e-5;
	PetscScalar atol = 1.0e-50;
	PetscScalar dtol = 1.0e3;
	PetscInt restart = 100;
	PetscInt maxitsGlobal = 200;

	// Coarsegrid solver
	PetscScalar coarse_rtol = 1.0e-8;
	PetscScalar coarse_atol = 1.0e-50;
	PetscScalar coarse_dtol = 1e3;
	PetscInt coarse_maxits = 30;
	PetscInt coarse_restart = 30;

	// Number of smoothening iterations per up/down smooth_sweeps
	PetscInt smooth_sweeps = 4;

	// Set up the solver
	ierr = KSPSetType(ksp,KSPFGMRES); // KSPCG, KSPGMRES
	CHKERRQ(ierr);

	ierr = KSPGMRESSetRestart(ksp,restart);
	CHKERRQ(ierr);

	ierr = KSPSetTolerances(ksp,rtol,atol,dtol,maxitsGlobal);
	CHKERRQ(ierr);

	ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
	CHKERRQ(ierr);

	ierr = KSPSetOperators(ksp,K,K);
	CHKERRQ(ierr);

	// The preconditinoer
	KSPGetPC(ksp,&pc);
	// Make PCMG the default solver
	PCSetType(pc,PCMG);

	// Set solver from options
	KSPSetFromOptions(ksp);

	// Get the prec again - check if it has changed
	KSPGetPC(ksp,&pc);

	// Flag for pcmg pc
	PetscBool pcmg_flag = PETSC_TRUE;
	PetscObjectTypeCompare((PetscObject)pc,PCMG,&pcmg_flag);

	// Only if PCMG is used
	if (pcmg_flag){ 

		// DMs for grid hierachy
		DM  *da_list,*daclist;
		Mat R;

		PetscMalloc(sizeof(DM)*opt->nlvls,&da_list);
		for (PetscInt k=0; k<opt->nlvls; k++) da_list[k] = NULL;
		PetscMalloc(sizeof(DM)*opt->nlvls,&daclist);
		for (PetscInt k=0; k<opt->nlvls; k++) daclist[k] = NULL;

		// Set 0 to the finest level
		daclist[0] = opt->da_nodes;

		// Coordinates
		PetscReal xmin=opt->xc[0], xmax=opt->xc[1], ymin=opt->xc[2], ymax=opt->xc[3], zmin=opt->xc[4], zmax=opt->xc[5];

		// Set up the coarse meshes
		DMCoarsenHierarchy(opt->da_nodes,opt->nlvls-1,&daclist[1]);
		for (PetscInt k=0; k<opt->nlvls; k++) {
			// NOTE: finest grid is nlevels - 1: PCMG MUST USE THIS ORDER ??? 
			da_list[k] = daclist[opt->nlvls-1-k];
			// THIS SHOULD NOT BE NECESSARY
			DMDASetUniformCoordinates(da_list[k],xmin,xmax,ymin,ymax,zmin,zmax);
		}

		// the PCMG specific options
		PCMGSetLevels(pc,opt->nlvls,NULL);
		PCMGSetType(pc,PC_MG_MULTIPLICATIVE); // Default
		PCMGSetCycleType(pc,PC_MG_CYCLE_V);
		PCMGSetGalerkin(pc,PETSC_TRUE);
		for (PetscInt k=1; k<opt->nlvls; k++) {
			DMCreateInterpolation(da_list[k-1],da_list[k],&R,NULL);
			PCMGSetInterpolation(pc,k,R);
			MatDestroy(&R);
		}

		// tidy up 
		for (PetscInt k=1; k<opt->nlvls; k++) { // DO NOT DESTROY LEVEL 0
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

			ierr = KSPGMRESSetRestart(cksp,coarse_restart);

			ierr = KSPSetTolerances(cksp,coarse_rtol,coarse_atol,coarse_dtol,coarse_maxits);
			// The preconditioner
			PC cpc;
			KSPGetPC(cksp,&cpc);
			PCSetType(cpc,PCSOR); // PCSOR, PCSPAI (NEEDS TO BE COMPILED), PCJACOBI     

			// Set smoothers on all levels (except for coarse grid):
			for (PetscInt k=1;k<opt->nlvls;k++){
				KSP dksp;
				PCMGGetSmoother(pc,k,&dksp);
				PC dpc;
				KSPGetPC(dksp,&dpc);
				ierr = KSPSetType(dksp,KSPGMRES); // KSPCG, KSPGMRES, KSPCHEBYSHEV (VERY GOOD FOR SPD)

				ierr = KSPGMRESSetRestart(dksp,smooth_sweeps);
				ierr = KSPSetTolerances(dksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,smooth_sweeps); // NOTE in the above maxitr=restart;
				PCSetType(dpc,PCSOR);// PCJACOBI, PCSOR for KSPCHEBYSHEV very good   
			}
		}
	}

	// Write check to screen:
	// Check the overall Krylov solver
	KSPType ksptype;
	KSPGetType(ksp,&ksptype);
	PCType pctype;
	PCGetType(pc,&pctype);
	PetscInt mmax;
	KSPGetTolerances(ksp,NULL,NULL,NULL,&mmax);
	PetscPrintf(PETSC_COMM_WORLD,"##############################################################\n");
	PetscPrintf(PETSC_COMM_WORLD,"################# Linear solver settings #####################\n");
	PetscPrintf(PETSC_COMM_WORLD,"# Main solver: %s, prec.: %s, maxiter.: %i \n",ksptype,pctype,mmax);

	// Only if pcmg is used
	if (pcmg_flag){
		// Check the smoothers and coarse grid solver:
		for (PetscInt k=0;k<opt->nlvls;k++){
			KSP dksp;
			PC dpc;
			KSPType dksptype;
			PCMGGetSmoother(pc,k,&dksp);
			KSPGetType(dksp,&dksptype);
			KSPGetPC(dksp,&dpc);
			PCType dpctype;
			PCGetType(dpc,&dpctype);
			PetscInt mmax;
			KSPGetTolerances(dksp,NULL,NULL,NULL,&mmax);
			PetscPrintf(PETSC_COMM_WORLD,"# Level %i smoother: %s, prec.: %s, sweep: %i \n",k,dksptype,dpctype,mmax);
		}
	}
	PetscPrintf(PETSC_COMM_WORLD,"##############################################################\n");


	return(ierr);
}

PetscErrorCode LinearElasticity::DMDAGetElements_3D(DM dm,PetscInt *nel,PetscInt *nen,const PetscInt *e[]) {
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

PetscInt LinearElasticity::Hex8Isoparametric(PetscScalar *X, PetscScalar *Y, PetscScalar *Z, PetscScalar nu, PetscInt redInt, PetscScalar *ke){
	// HEX8_ISOPARAMETRIC - Computes HEX8 isoparametric element matrices
	// The element stiffness matrix is computed as:
	//
	//       ke = int(int(int(B^T*C*B,x),y),z)
	//
	// For an isoparameteric element this integral becomes:
	//
	//       ke = int(int(int(B^T*C*B*det(J),xi=-1..1),eta=-1..1),zeta=-1..1)
	//
	// where B is the more complicated expression:
	// B = [dx*alpha1 + dy*alpha2 + dz*alpha3]*N
	// where
	// dx = [invJ11 invJ12 invJ13]*[dxi deta dzeta]
	// dy = [invJ21 invJ22 invJ23]*[dxi deta dzeta]
	// dy = [invJ31 invJ32 invJ33]*[dxi deta dzeta]
	//
	// Remark: The elasticity modulus is left out in the below
	// computations, because we multiply with them afterwards (the aim is
	// topology optimization).
	// Furthermore, this is not the most efficient code, but it is readible.
	//
	/////////////////////////////////////////////////////////////////////////////////
	//////// INPUT:
	// X, Y, Z  = Vectors containing the coordinates of the eight nodes
	//               (x1,y1,z1,x2,y2,z2,...,x8,y8,z8). Where node 1 is in the lower
	//               left corner, and node 2 is the next node counterclockwise
	//               (looking in the negative z-dir).
	//               Finish the x-y-plane and then move in the positive z-dir.
	// redInt   = Reduced integration option boolean (here an integer).
	//           	redInt == 0 (false): Full integration
	//           	redInt == 1 (true): Reduced integration
	// nu 		= Poisson's ratio.
	//
	//////// OUTPUT:
	// ke  = Element stiffness matrix. Needs to be multiplied with elasticity modulus
	//
	//   Written 2013 at
	//   Department of Mechanical Engineering
	//   Technical University of Denmark (DTU).
	/////////////////////////////////////////////////////////////////////////////////

	//// COMPUTE ELEMENT STIFFNESS MATRIX
	// Lame's parameters (with E=1.0):
	PetscScalar lambda = nu/((1.0+nu)*(1.0-2.0*nu));
	PetscScalar mu = 1.0/(2.0*(1.0+nu));
	// Constitutive matrix
	PetscScalar C[6][6] = {{lambda+2.0*mu, lambda, lambda, 0.0, 0.0, 0.0},
		{lambda, lambda+2.0*mu, lambda, 0.0, 0.0, 0.0},
		{lambda, lambda, lambda+2.0*mu, 0.0, 0.0, 0.0},
		{0.0,    0.0,    0.0,           mu,  0.0, 0.0},
		{0.0, 	0.0, 	0.0, 		   0.0, mu,  0.0}, 
		{0.0, 	0.0,	0.0, 		   0.0, 0.0, mu}};
	// Gauss points (GP) and weigths
	// Two Gauss points in all directions (total of eight)
	PetscScalar GP[2] = {-0.577350269189626, 0.577350269189626}; 
	// Corresponding weights
	PetscScalar W[2] = {1.0, 1.0};
	// If reduced integration only use one GP
	if (redInt){
		GP[1] = 0.0;
		W[1] = 2.0;
	}
	// Matrices that help when we gather the strain-displacement matrix:
	PetscScalar alpha1[6][3]; PetscScalar alpha2[6][3]; PetscScalar alpha3[6][3];
	memset(alpha1, 0, sizeof(alpha1[0][0])*6*3); // zero out
	memset(alpha2, 0, sizeof(alpha2[0][0])*6*3); // zero out
	memset(alpha3, 0, sizeof(alpha3[0][0])*6*3); // zero out
	alpha1[0][0] = 1.0; alpha1[3][1] = 1.0; alpha1[5][2] = 1.0;
	alpha2[1][1] = 1.0; alpha2[3][0] = 1.0; alpha2[4][2] = 1.0;
	alpha3[2][2] = 1.0; alpha3[4][1] = 1.0; alpha3[5][0] = 1.0;
	PetscScalar dNdxi[8]; PetscScalar dNdeta[8]; PetscScalar dNdzeta[8];
	PetscScalar J[3][3];
	PetscScalar invJ[3][3];
	PetscScalar beta[6][3];
	PetscScalar B[6][24]; // Note: Small enough to be allocated on stack
	PetscScalar *dN;
	// Make sure the stiffness matrix is zeroed out:
	memset(ke, 0, sizeof(ke[0])*24*24);
	// Perform the numerical integration
	for (PetscInt ii=0; ii<2-redInt; ii++){
		for (PetscInt jj=0; jj<2-redInt; jj++){
			for (PetscInt kk=0; kk<2-redInt; kk++){
				// Integration point
				PetscScalar xi = GP[ii]; 
				PetscScalar eta = GP[jj]; 
				PetscScalar zeta = GP[kk];
				// Differentiated shape functions
				DifferentiatedShapeFunctions(xi, eta, zeta, dNdxi, dNdeta, dNdzeta);
				// Jacobian
				J[0][0] = Dot(dNdxi,X,8); J[0][1] = Dot(dNdxi,Y,8); J[0][2] = Dot(dNdxi,Z,8);
				J[1][0] = Dot(dNdeta,X,8); J[1][1] = Dot(dNdeta,Y,8); J[1][2] = Dot(dNdeta,Z,8);
				J[2][0] = Dot(dNdzeta,X,8); J[2][1] = Dot(dNdzeta,Y,8); J[2][2] = Dot(dNdzeta,Z,8);
				// Inverse and determinant
				PetscScalar detJ = Inverse3M(J, invJ);
				// Weight factor at this point
				PetscScalar weight = W[ii]*W[jj]*W[kk]*detJ;
				// Strain-displacement matrix
				memset(B, 0, sizeof(B[0][0])*6*24); // zero out
				for (PetscInt ll=0; ll<3; ll++){
					// Add contributions from the different derivatives
					if (ll==0) {dN = dNdxi;}
					if (ll==1) {dN = dNdeta;}
					if (ll==2) {dN = dNdzeta;}
					// Assemble strain operator
					for (PetscInt i=0; i<6; i++){
						for (PetscInt j=0; j<3; j++){
							beta[i][j] = invJ[0][ll]*alpha1[i][j]
								+invJ[1][ll]*alpha2[i][j]
								+invJ[2][ll]*alpha3[i][j];
						}
					}
					// Add contributions to strain-displacement matrix
					for (PetscInt i=0; i<6; i++){
						for (PetscInt j=0; j<24; j++){
							B[i][j] = B[i][j] + beta[i][j%3]*dN[j/3];
						}
					}
				}
				// Finally, add to the element matrix
				for (PetscInt i=0; i<24; i++){
					for (PetscInt j=0; j<24; j++){
						for (PetscInt k=0; k<6; k++){
							for (PetscInt l=0; l<6; l++){
								
								ke[j+24*i] = ke[j+24*i] + weight*(B[k][i] * C[k][l] * B[l][j]);
							}
						}
					}
				}
			}
		}
	}
	return 0;
}
PetscScalar LinearElasticity::Dot(PetscScalar *v1, PetscScalar *v2, PetscInt l){
	// Function that returns the dot product of v1 and v2,
	// which must have the same length l
	PetscScalar result = 0.0;
	for (PetscInt i=0; i<l; i++){
		result = result + v1[i]*v2[i];
	}
	return result;
}

void LinearElasticity::DifferentiatedShapeFunctions(PetscScalar xi, PetscScalar eta, PetscScalar zeta, PetscScalar *dNdxi, PetscScalar *dNdeta, PetscScalar *dNdzeta){
	//differentiatedShapeFunctions - Computes differentiated shape functions
	// At the point given by (xi, eta, zeta).
	// With respect to xi:
	dNdxi[0]  = -0.125*(1.0-eta)*(1.0-zeta);
	dNdxi[1]  =  0.125*(1.0-eta)*(1.0-zeta);
	dNdxi[2]  =  0.125*(1.0+eta)*(1.0-zeta);
	dNdxi[3]  = -0.125*(1.0+eta)*(1.0-zeta);
	dNdxi[4]  = -0.125*(1.0-eta)*(1.0+zeta);
	dNdxi[5]  =  0.125*(1.0-eta)*(1.0+zeta);
	dNdxi[6]  =  0.125*(1.0+eta)*(1.0+zeta);
	dNdxi[7]  = -0.125*(1.0+eta)*(1.0+zeta);
	// With respect to eta:
	dNdeta[0] = -0.125*(1.0-xi)*(1.0-zeta);
	dNdeta[1] = -0.125*(1.0+xi)*(1.0-zeta);
	dNdeta[2] =  0.125*(1.0+xi)*(1.0-zeta);
	dNdeta[3] =  0.125*(1.0-xi)*(1.0-zeta);
	dNdeta[4] = -0.125*(1.0-xi)*(1.0+zeta);
	dNdeta[5] = -0.125*(1.0+xi)*(1.0+zeta);
	dNdeta[6] =  0.125*(1.0+xi)*(1.0+zeta);
	dNdeta[7] =  0.125*(1.0-xi)*(1.0+zeta);
	// With respect to zeta:
	dNdzeta[0]= -0.125*(1.0-xi)*(1.0-eta);
	dNdzeta[1]= -0.125*(1.0+xi)*(1.0-eta);
	dNdzeta[2]= -0.125*(1.0+xi)*(1.0+eta);
	dNdzeta[3]= -0.125*(1.0-xi)*(1.0+eta);
	dNdzeta[4]=  0.125*(1.0-xi)*(1.0-eta);
	dNdzeta[5]=  0.125*(1.0+xi)*(1.0-eta);
	dNdzeta[6]=  0.125*(1.0+xi)*(1.0+eta);
	dNdzeta[7]=  0.125*(1.0-xi)*(1.0+eta);
}

PetscScalar LinearElasticity::Inverse3M(PetscScalar J[][3], PetscScalar invJ[][3]){
	//inverse3M - Computes the inverse of a 3x3 matrix
	PetscScalar detJ = J[0][0]*(J[1][1]*J[2][2]-J[2][1]*J[1][2])-J[0][1]*(J[1][0]*J[2][2]-J[2][0]*J[1][2])+J[0][2]*(J[1][0]*J[2][1]-J[2][0]*J[1][1]);
	invJ[0][0] = (J[1][1]*J[2][2]-J[2][1]*J[1][2])/detJ;
	invJ[0][1] = -(J[0][1]*J[2][2]-J[0][2]*J[2][1])/detJ;
	invJ[0][2] = (J[0][1]*J[1][2]-J[0][2]*J[1][1])/detJ;
	invJ[1][0] = -(J[1][0]*J[2][2]-J[1][2]*J[2][0])/detJ;
	invJ[1][1] = (J[0][0]*J[2][2]-J[0][2]*J[2][0])/detJ;
	invJ[1][2] = -(J[0][0]*J[1][2]-J[0][2]*J[1][0])/detJ;
	invJ[2][0] = (J[1][0]*J[2][1]-J[1][1]*J[2][0])/detJ;
	invJ[2][1] = -(J[0][0]*J[2][1]-J[0][1]*J[2][0])/detJ;
	invJ[2][2] = (J[0][0]*J[1][1]-J[1][0]*J[0][1])/detJ;
	return detJ;
}

