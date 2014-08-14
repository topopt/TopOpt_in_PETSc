#include <TopOpt.h>
#include <cmath>

/*
 Authors: Niels Aage, Erik Andreassen, Boyan Lazarov, August 2013

 Disclaimer:                                                              
 The authors reserves all rights but does not guaranty that the code is   
 free from errors. Furthermore, we shall not be liable in any event     
 caused by the use of the program.                                     
*/


TopOpt::TopOpt(){
      
  x=NULL;
  xPhys=NULL;
  dfdx=NULL;
  dgdx=NULL;
  gx=NULL;
  da_nodes=NULL;
  da_elem=NULL;
  
  xo1=NULL;
  xo2=NULL;
  U=NULL;
  L=NULL;
  
  SetUp();
  
}

TopOpt::~TopOpt(){
  
  
    // Delete vectors
    if (x!=NULL){ VecDestroy(&x); }
    if (dfdx!=NULL){ VecDestroy(&dfdx); }
    if (dgdx!=NULL){ VecDestroyVecs(m,&dgdx); }
    // Densities
    if (xPhys!=NULL){ VecDestroy(&xPhys); }
    if (xold!=NULL){ VecDestroy(&xold); }
    if (xmin!=NULL){ VecDestroy(&xmin); }
    if (xmax!=NULL){ VecDestroy(&xmax); }
    
    if (da_nodes!=NULL){ DMDestroy(&(da_nodes)); }
    if (da_elem!=NULL){ DMDestroy(&(da_elem)); }
    
    // Delete constraints
    if (gx!=NULL){ delete [] gx; }
    
    // mma restart method    		
    if (xo1!=NULL){ VecDestroy(&xo1); }
    if (xo2!=NULL){ VecDestroy(&xo2); }
    if (L!=NULL){ VecDestroy(&L); }
    if (U!=NULL){ VecDestroy(&U);  }
}


// NO METHODS !
//PetscErrorCode TopOpt::SetUp(Vec CRAPPY_VEC){
PetscErrorCode TopOpt::SetUp(){
	PetscErrorCode ierr;
  
	// SET DEFAULTS for FE mesh and levels for MG solver
        nxyz[0] = 65;//129;
        nxyz[1] = 33;//65;
        nxyz[2] = 33;//65;
        xc[0] = 0.0;
        xc[1] = 2.0;
        xc[2] = 0.0;
        xc[3] = 1.0;
        xc[4] = 0.0;
        xc[5] = 1.0;
	nu=0.3;
        nlvls = 4;
	
	// SET DEFAULTS for optimization problems
	volfrac = 0.12;
	maxItr = 400;
        rmin = 0.08;
        penal = 3.0;
        Emin = 1.0e-9;
        Emax = 1.0;
        filter = 0; // 0=sens,1=dens,2=PDE - other val == no filtering
        m = 1; // volume constraint
        Xmin = 0.0;
        Xmax = 1.0;
        movlim = 0.2;
	
	ierr = SetUpMESH(); CHKERRQ(ierr);

	ierr = SetUpOPT(); CHKERRQ(ierr);

	return(ierr);
}


PetscErrorCode TopOpt::SetUpMESH(){
	
	PetscErrorCode ierr;
	
	// Read input from arguments
	PetscBool flg;
	
	// Physics parameters
	PetscOptionsGetInt(NULL,"-nx",&(nxyz[0]),&flg);
	PetscOptionsGetInt(NULL,"-ny",&(nxyz[1]),&flg);
	PetscOptionsGetInt(NULL,"-nz",&(nxyz[2]),&flg);
	PetscOptionsGetReal(NULL,"-xcmin",&(xc[0]),&flg);	
	PetscOptionsGetReal(NULL,"-xcmax",&(xc[1]),&flg);
	PetscOptionsGetReal(NULL,"-ycmin",&(xc[2]),&flg);
	PetscOptionsGetReal(NULL,"-ycmax",&(xc[3]),&flg);
	PetscOptionsGetReal(NULL,"-zcmin",&(xc[4]),&flg);
	PetscOptionsGetReal(NULL,"-zcmax",&(xc[5]),&flg);
        PetscOptionsGetReal(NULL,"-penal",&penal,&flg);
// 	PetscOptionsGetReal(NULL,"-nu",&(nu),&flg);
	PetscOptionsGetInt(NULL,"-nlvls",&nlvls,&flg);

	
	// Write parameters for the physics _ OWNED BY TOPOPT
	PetscPrintf(PETSC_COMM_WORLD,"########################################################################\n");
	PetscPrintf(PETSC_COMM_WORLD,"############################ FEM settings ##############################\n");
	PetscPrintf(PETSC_COMM_WORLD,"# Number of nodes: (-nx,-ny,-nz):        (%i,%i,%i) \n",nxyz[0],nxyz[1],nxyz[2]);
        PetscPrintf(PETSC_COMM_WORLD,"# Number of degree of freedom:           %i \n",3*nxyz[0]*nxyz[1]*nxyz[2]);
	PetscPrintf(PETSC_COMM_WORLD,"# Number of elements:                    (%i,%i,%i) \n",nxyz[0]-1,nxyz[1]-1,nxyz[2]-1);
	PetscPrintf(PETSC_COMM_WORLD,"# Dimensions: (-xcmin,-xcmax,..,-zcmax): (%f,%f,%f)\n",xc[1]-xc[0],xc[3]-xc[2],xc[5]-xc[4]);
	PetscPrintf(PETSC_COMM_WORLD,"# -nlvls: %i\n",nlvls);
       	PetscPrintf(PETSC_COMM_WORLD,"########################################################################\n");

	// Check if the mesh supports the chosen number of MG levels
	PetscScalar divisor = PetscPowScalar(2.0,(PetscScalar)nlvls-1.0);
	// x - dir
	if ( std::floor((PetscScalar)(nxyz[0]-1)/divisor) != (nxyz[0]-1.0)/((PetscInt)divisor) ) {
		PetscPrintf(PETSC_COMM_WORLD,"MESH DIMENSION NOT COMPATIBLE WITH NUMBER OF MULTIGRID LEVELS!\n");
                PetscPrintf(PETSC_COMM_WORLD,"X - number of nodes %i is cannot be halfened %i times\n",nxyz[0],nlvls-1);
		exit(0);
	}	
	// y - dir
        if ( std::floor((PetscScalar)(nxyz[1]-1)/divisor) != (nxyz[1]-1.0)/((PetscInt)divisor) ) {
                PetscPrintf(PETSC_COMM_WORLD,"MESH DIMENSION NOT COMPATIBLE WITH NUMBER OF MULTIGRID LEVELS!\n");
                PetscPrintf(PETSC_COMM_WORLD,"Y - number of nodes %i is cannot be halfened %i times\n",nxyz[1],nlvls-1);
		exit(0);
        }
	// z - dir
        if ( std::floor((PetscScalar)(nxyz[2]-1)/divisor) != (nxyz[2]-1.0)/((PetscInt)divisor) ) {
                PetscPrintf(PETSC_COMM_WORLD,"MESH DIMENSION NOT COMPATIBLE WITH NUMBER OF MULTIGRID LEVELS!\n");
                PetscPrintf(PETSC_COMM_WORLD,"Z - number of nodes %i is cannot be halfened %i times\n",nxyz[2],nlvls-1);
		exit(0);
        }


	// Start setting up the FE problem
	// Boundary types: DMDA_BOUNDARY_NONE, DMDA_BOUNDARY_GHOSTED, DMDA_BOUNDARY_PERIODIC
	DMBoundaryType bx = DM_BOUNDARY_NONE;
	DMBoundaryType by = DM_BOUNDARY_NONE;
	DMBoundaryType bz = DM_BOUNDARY_NONE;

	// Stencil type - box since this is closest to FEM (i.e. STAR is FV/FD)
	DMDAStencilType  stype = DMDA_STENCIL_BOX;

	// Discretization: nodes:
	// For standard FE - number must be odd
	// FOr periodic: Number must be even
	PetscInt nx = nxyz[0];
	PetscInt ny = nxyz[1];
	PetscInt nz = nxyz[2];

	// number of nodal dofs
	PetscInt numnodaldof = 3;

	// Stencil width: each node connects to a box around it - linear elements
	PetscInt stencilwidth = 1;

	// Coordinates and element sizes: note that dx,dy,dz are half the element size
	PetscReal xmin=xc[0], xmax=xc[1], ymin=xc[2], ymax=xc[3], zmin=xc[4], zmax=xc[5];
	dx = (xc[1]-xc[0])/(PetscScalar(nxyz[0]-1));
	dy = (xc[3]-xc[2])/(PetscScalar(nxyz[1]-1));
	dz = (xc[5]-xc[4])/(PetscScalar(nxyz[2]-1));

	// Create the nodal mesh
	ierr = DMDACreate3d(PETSC_COMM_WORLD,bx,by,bz,stype,nx,ny,nz,PETSC_DECIDE,PETSC_DECIDE,PETSC_DECIDE,
			numnodaldof,stencilwidth,0,0,0,&(da_nodes));
	CHKERRQ(ierr);

	// Set the coordinates
	ierr = DMDASetUniformCoordinates(da_nodes, xmin,xmax, ymin,ymax, zmin,zmax);
	CHKERRQ(ierr);
	
	// Set the element type to Q1: Otherwise calls to GetElements will change to P1 !
	// STILL DOESN*T WORK !!!!
	ierr = DMDASetElementType(da_nodes, DMDA_ELEMENT_Q1);
	CHKERRQ(ierr);
	
	// Create the element mesh: NOTE THIS DOES NOT INCLUDE THE FILTER !!!
	// find the geometric partitioning of the nodal mesh, so the element mesh will coincide 
	// with the nodal mesh
	PetscInt md,nd,pd; 
	ierr = DMDAGetInfo(da_nodes,NULL,NULL,NULL,NULL,&md,&nd,&pd,NULL,NULL,NULL,NULL,NULL,NULL);
	CHKERRQ(ierr);
	
	// vectors with partitioning information
	PetscInt *Lx=new PetscInt[md];
	PetscInt *Ly=new PetscInt[nd];
	PetscInt *Lz=new PetscInt[pd];

	// get number of nodes for each partition
	const PetscInt *LxCorrect, *LyCorrect, *LzCorrect;
	ierr = DMDAGetOwnershipRanges(da_nodes, &LxCorrect, &LyCorrect, &LzCorrect); 
	CHKERRQ(ierr);
	
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

	// Create the element grid: NOTE CONNECTIVITY IS 0
	PetscInt conn = 0;
	ierr = DMDACreate3d(PETSC_COMM_WORLD,bx,by,bz,stype,nx-1,ny-1,nz-1,md,nd,pd,
			1,conn,Lx,Ly,Lz,&(da_elem));
	CHKERRQ(ierr);
	
	// Set element center coordinates
	ierr = DMDASetUniformCoordinates(da_elem , xmin+dx/2.0,xmax-dx/2.0, ymin+dy/2.0,ymax-dy/2.0, zmin+dz/2.0,zmax-dz/2.0);
	CHKERRQ(ierr);

	// Clean up
	delete [] Lx;
	delete [] Ly;
	delete [] Lz;
  
  
  	return(ierr);
}

PetscErrorCode TopOpt::SetUpOPT(){
  
	PetscErrorCode ierr;
  
	//ierr = VecDuplicate(CRAPPY_VEC,&xPhys); CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(da_elem,&xPhys);  CHKERRQ(ierr);
	// Total number of design variables
	VecGetSize(xPhys,&n);

	PetscBool flg;
	
	// Optimization paramteres
	PetscOptionsGetReal(NULL,"-Emin",&Emin,&flg);
	PetscOptionsGetReal(NULL,"-Emax",&Emax,&flg);
	PetscOptionsGetReal(NULL,"-volfrac",&volfrac,&flg);
        PetscOptionsGetReal(NULL,"-penal",&penal,&flg);
	PetscOptionsGetReal(NULL,"-rmin",&rmin,&flg);
	PetscOptionsGetInt(NULL,"-maxItr",&maxItr,&flg);
	PetscOptionsGetInt(NULL,"-filter",&filter,&flg);
	PetscOptionsGetReal(NULL,"-Xmin",&Xmin,&flg);
        PetscOptionsGetReal(NULL,"-Xmax",&Xmax,&flg);
	PetscOptionsGetReal(NULL,"-movlim",&movlim,&flg);
        
	PetscPrintf(PETSC_COMM_WORLD,"################### Optimization settings ####################\n");
	PetscPrintf(PETSC_COMM_WORLD,"# Problem size: n= %i, m= %i\n",n,m);
	PetscPrintf(PETSC_COMM_WORLD,"# -filter: %i  (0=sens., 1=dens, 2=PDE)\n",filter);
	PetscPrintf(PETSC_COMM_WORLD,"# -rmin: %f\n",rmin);
	PetscPrintf(PETSC_COMM_WORLD,"# -volfrac: %f\n",volfrac);
        PetscPrintf(PETSC_COMM_WORLD,"# -penal: %f\n",penal);
	PetscPrintf(PETSC_COMM_WORLD,"# -Emin/-Emax: %e - %e \n",Emin,Emax);
	PetscPrintf(PETSC_COMM_WORLD,"# -maxItr: %i\n",maxItr);
	PetscPrintf(PETSC_COMM_WORLD,"# -movlim: %f\n",movlim);
       	PetscPrintf(PETSC_COMM_WORLD,"##############################################################\n");

        // Allocate after input
        gx = new PetscScalar[m];
	if (filter==0){
		Xmin = 0.001; // Prevent division by zero in filter
	}
	
	// Allocate the optimization vectors
	ierr = VecDuplicate(xPhys,&x); CHKERRQ(ierr);
	VecSet(x,volfrac); // Initialize to volfrac !
	VecSet(xPhys,volfrac); // Initialize to volfrac !
  
	// Sensitivity vectors
	ierr = VecDuplicate(x,&dfdx); CHKERRQ(ierr);
	ierr = VecDuplicateVecs(x,m, &dgdx); CHKERRQ(ierr);

	// Bounds and 
	VecDuplicate(x,&xmin);
	VecDuplicate(x,&xmax);
	VecDuplicate(x,&xold);	
	VecSet(xold,volfrac);
  
  	return(ierr);
}
  
void TopOpt::AllocMMAwithRestart(int *itr, MMA **mma)  {

  	// Check if restart is desired
	restart = PETSC_FALSE; // DEFAULT DOES NOT USE RESTART
	flip = PETSC_TRUE;     // BOOL to ensure that two dump streams are kept
	
	PetscBool flg;
	PetscOptionsGetBool(NULL,"-restart",&restart,&flg);
	// default name of the restart dir
        std::string filename = "./";

        // Check PETSc input for a data directory
        char filenameChar[PETSC_MAX_PATH_LEN];
        PetscOptionsGetString(NULL,"-restartDir",filenameChar,sizeof(filenameChar),&flg);
	// If input, change path of the file in filename
        if (flg){
                filename="";
                filename.append(filenameChar);
        }
	
	// Which solution to use for restarting
	PetscInt restartNumber;
	PetscOptionsGetInt(NULL,"-restartNumber",&restartNumber,&flg);
	if (!flg){
		restartNumber=1;
	}

	PetscPrintf(PETSC_COMM_WORLD,"##############################################################\n");
	PetscPrintf(PETSC_COMM_WORLD,"# Continue from previous iteration (-restart): %i \n",restart);
	PetscPrintf(PETSC_COMM_WORLD,"# Restart files are located in folder (-restartDir): %s \n",filename.c_str());
	
	// Append the dummyname for restart files	
	filename.append("/restore_V");
	
	PetscPrintf(PETSC_COMM_WORLD,"# The restart point is restore_V%i****.dat  (where %i is the -restartNumber) \n",restartNumber,restartNumber);
	
	// RESTORE FROM BREAKDOWN
	PetscInt myrank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &myrank);
        std::ifstream indes;
	restdens_1=filename;
	restdens_1.append("1");
	restdens_2=filename;
	restdens_2.append("2");

	std::stringstream ss;
	if(myrank<10){
                ss<<restdens_1<<"_0000"<<myrank<<".dat";
		ss<<" ";
		ss<<restdens_2<<"_0000"<<myrank<<".dat";
        }
	else
        if(myrank<100){
		ss<<restdens_1<<"_000"<<myrank<<".dat";
		ss<<" ";
		ss<<restdens_2<<"_000"<<myrank<<".dat";
        }
	else
        if(myrank<1000){
                ss<<restdens_1<<"_00"<<myrank<<".dat";
		ss<<" ";
		ss<<restdens_2<<"_00"<<myrank<<".dat";
        }
	else 
	if(myrank<10000){
                ss<<restdens_1<<"_0"<<myrank<<".dat";
                ss<<" ";
                ss<<restdens_2<<"_0"<<myrank<<".dat";
        }

        ss>>restdens_1;
	ss>>restdens_2;
	
	// Allocate the data needed for a MMA restart
	VecDuplicate(x,&xo1);
	VecDuplicate(x,&xo2);
	VecDuplicate(x,&U);
	VecDuplicate(x,&L);
	
	// Read from restart point
	if (restartNumber==1){
		indes.open(restdens_1.c_str(),std::ios::in);
	}
	else
	if (restartNumber==2){
                indes.open(restdens_2.c_str(),std::ios::in);
	}

	if(indes && restart)
        {
                PetscInt nlocsiz;
		PetscScalar *xp, *xpp, *xo1p, *xo2p, *Up, *Lp;
                	
		VecGetArray(x,&xp);
		VecGetArray(xPhys,&xpp);
		
		VecGetArray(xo1,&xo1p);
		VecGetArray(xo2,&xo2p);
		VecGetArray(U,&Up);
		VecGetArray(L,&Lp);
		
                indes.read((char*)&nlocsiz,sizeof(PetscInt));
                indes.read((char*)xp,sizeof(PetscScalar)*nlocsiz);
                indes.read((char*)xpp,sizeof(PetscScalar)*nlocsiz);
                indes.read((char*)xo1p,sizeof(PetscScalar)*nlocsiz);
                indes.read((char*)xo2p,sizeof(PetscScalar)*nlocsiz);
                indes.read((char*)Up,sizeof(PetscScalar)*nlocsiz);
                indes.read((char*)Lp,sizeof(PetscScalar)*nlocsiz);
                indes.read((char*)itr,sizeof(PetscInt));
		indes.read((char*)&fscale,sizeof(PetscScalar));
                indes.close();

		VecRestoreArray(x,&xp);
		VecRestoreArray(xPhys,&xpp);
		VecRestoreArray(xo1,&xo1p);
		VecRestoreArray(xo2,&xo2p);
		VecRestoreArray(U,&Up);
		VecRestoreArray(L,&Lp);
		
		*mma = new MMA(n,m,*itr,xo1,xo2,U,L);
	
		if (restartNumber==1){
			PetscPrintf(PETSC_COMM_WORLD,"# Succesfull restart from from file (starting from): %s \n",restdens_1.c_str());
		}
		else
        	if (restartNumber==2){
			PetscPrintf(PETSC_COMM_WORLD,"# Succesfull restart from from file (starting from): %s \n",restdens_2.c_str());	
		}
	

	}
	else {
		*mma = new MMA(n,m,x);
	}  
  
	indes.close();
  
}  

void TopOpt::WriteRestartFiles(int *itr, MMA *mma) {

	// Always dump data if correct allocater has been used
	if (xo1!=NULL){
	  	// Get data from MMA
		mma->Restart(xo1,xo2,U,L);

		// Open the stream (and make sure there always is one working copy)
		std::string dens_iter;
		std::stringstream ss_iter;
		if (flip) {
		    ss_iter << restdens_1;
		    ss_iter >> dens_iter;
		    flip = PETSC_FALSE;
		} else {
		    ss_iter << restdens_2;
		    ss_iter >> dens_iter;
		    flip = PETSC_TRUE;
		}

		// Open stream
		std::ofstream out(dens_iter.c_str(),std::ios::out);

		// poniters to data
		PetscInt nlocsiz;
		VecGetLocalSize(x,&nlocsiz);
		PetscScalar *xp, *xpp, *xo1p, *xo2p, *Up, *Lp;
                	
		VecGetArray(x,&xp);
		VecGetArray(xPhys,&xpp);
		VecGetArray(xo1,&xo1p);
		VecGetArray(xo2,&xo2p);
		VecGetArray(U,&Up);
		VecGetArray(L,&Lp);
	
		// Write to file
		out.write((char*)&nlocsiz,sizeof(PetscInt));
		out.write((char*)xp,sizeof(PetscScalar)*nlocsiz);
                out.write((char*)xpp,sizeof(PetscScalar)*nlocsiz);
                out.write((char*)xo1p,sizeof(PetscScalar)*nlocsiz);
                out.write((char*)xo2p,sizeof(PetscScalar)*nlocsiz);
                out.write((char*)Up,sizeof(PetscScalar)*nlocsiz);
                out.write((char*)Lp,sizeof(PetscScalar)*nlocsiz);
                out.write((char*)itr,sizeof(PetscInt));
		out.write((char*)&fscale,sizeof(PetscScalar));
		out.close();
		
		// Tidy up
		VecRestoreArray(x,&xp);
		VecRestoreArray(xPhys,&xpp);
		VecRestoreArray(xo1,&xo1p);
		VecRestoreArray(xo2,&xo2p);
		VecRestoreArray(U,&Up);
		VecRestoreArray(L,&Lp);
	}

  
}
