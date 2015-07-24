// Import necessary stuff
#include <iostream>
#include <MPIIO.h>
#include <cstdlib> // To get the exit function

/* -----------------------------------------------------------------------------
Authors: Niels Aage, Erik Andreassen, Boyan Lazarov, August 2013 
Copyright (C) 2013-2014,

This MPIIO implementation is licensed under Version 2.1 of the GNU
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

// Constructor
#define _NO_SUCH_FILE 35

MPIIO::MPIIO(DM da_nodes, int nPf, std::string pnames,int nCf, std::string cnames){
	
	// User defined string
	std::string infoString = "TopOpt result version 1.1";
	// Maximum number of points per element
	int nPEl = 8;
	// Number of domains (where domain refers to different regions in the optimization)
	const int nDom = 1;
	
	// Number of point fields per domain:
	int nPFields[nDom]= {nPf};
	// The names of the point fields
	std::string pFieldNames = pnames;
	// Number of cell (element) fields per domain:
	int nCFields[nDom]= {nCf};
	// The names of the cell fields
	std::string cFieldNames = cnames;

	
	// Points/cells in each of the domains:
	nPointsMyrank = new unsigned long int[nDom];
	nCellsMyrank = new unsigned long int[nDom];

	// ----------------- Get the mesh: -------------------
	// Coordinates and a pointer
	Vec coordinates;
	PetscScalar *coordinatesPointer;

	// Get coordinates in local node numbering (including ghosts) 
	DMGetCoordinatesLocal(da_nodes,&coordinates); 
	VecGetArray(coordinates,&coordinatesPointer); 

	// Get the global dof number (and ghosts)
	PetscInt nn;
	VecGetSize(coordinates,&nn); 
		
	// Get the FE mesh structure (from the nodal mesh)
	PetscInt nel, nen;
	const PetscInt *necon;
	DMDAGetElements_3D(da_nodes,&nel,&nen,&necon);

	// Number of points/cells in each domain for this rank
	PetscInt numnodaldof = 3;
	nPointsMyrank[0] = nn/numnodaldof;
	nCellsMyrank[0] = nel; // We have this number from when we called DMDAGetElements_3D

	// --------- Allocate the output object: -------
	Allocate(infoString, nDom, nPFields, nCFields, nPointsMyrank, nCellsMyrank, nPEl, pFieldNames, cFieldNames);
	// --------- Allocate the output object: -------
	
	// Write the points (or coordinates of the points)
	float *pointsDomain0 = new float[3*nPointsMyrank[0]];
	for (unsigned long int i=0; i<3*nPointsMyrank[0]; i++){
		// Convert to single precission to save space
		pointsDomain0[i] = float(coordinatesPointer[i]);
	}
	writePoints(0, pointsDomain0);
	// Restore coordinates array
	VecRestoreArray(coordinates,&coordinatesPointer); 

	// Run through the elements of the domain (we have already called DMDAGetElements)
	unsigned long int *cellsDomain0 = new unsigned long int[nPEl*nCellsMyrank[0]];
	unsigned long int *cellsOffset0 = new unsigned long int[nCellsMyrank[0]]; 
	unsigned long int *cellsTypes0 = new unsigned long int[nCellsMyrank[0]];
	unsigned long int CellOffset = 0;
	for (unsigned long int i=0;i<nCellsMyrank[0];i++){
		// Element type is the first number outputted:
		if (nen == 8){ // Hex element
			cellsTypes0[i] = 12; // (in vtk hex element type is 12)
		}
		// Then run through the nodes
		for (int j=0;j<nen;j++){
			cellsDomain0[i*nPEl+j] = necon[i*nen+j];
		}
		// Create the offset
		if (nen == 8){ // Hex element
			CellOffset+=nen;
                        cellsOffset0[i] = CellOffset; // (in vtk hex element type is 12)
                } 
		// Finally, in case we have varying elements size, make sure the extra nodes
		// are put to zero
		// REMARK: This is only an example, and is never used in this implementation
		for(int j=8; j<nPEl; j++){
			cellsDomain0[i*(nPEl+1)+j+1] = 0;
		}
	}
	writeCells(0, cellsDomain0,cellsOffset0,cellsTypes0); // First domain

	// Allocate working arrays for outputting fields from timesteps:
	workPointField = new float[nPointsMyrank[0]*nPFields[0]]; // For first domain
	workCellField  = new float[nCellsMyrank[0]*nCFields[0]];  // For first domain

	delete [] pointsDomain0;
	delete [] cellsDomain0;
	delete [] cellsOffset0;
	delete [] cellsTypes0;

  
}

// Destructor
MPIIO::~MPIIO(){
	// Delete the allocated arrays
	delete [] workPointField;
	delete [] workCellField;
	delete [] nPointsMyrank;
	delete [] nCellsMyrank;
	
	delete [] nPoints;
	delete [] nCells;
	delete [] nPointsT;
	delete [] nCellsT;
	delete [] nPFields;
	delete [] nCFields;
}


PetscErrorCode MPIIO::WriteVTK(DM da_nodes, Vec U, TopOpt *opt, PetscInt itr){
  
  // Here we only have one "timestep" (no optimization)
	unsigned long int timestep = itr;

	PetscErrorCode ierr;

	// POINT FIELD(S)
	// Displacement
	Vec Ulocal;
	DMCreateLocalVector(da_nodes,&Ulocal);
	ierr = VecSet(Ulocal, 0.0); CHKERRQ(ierr);
	// Update the local vector from global solution
	ierr = DMGlobalToLocalBegin(da_nodes,U,INSERT_VALUES,Ulocal); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(da_nodes,U,INSERT_VALUES,Ulocal); CHKERRQ(ierr);
	// We need a pointer to the local vector
	PetscScalar *UlocalPointer;
	ierr = VecGetArray(Ulocal, &UlocalPointer); CHKERRQ(ierr);
	for (unsigned long int i=0; i<nPointsMyrank[0]; i++){
		// Ux
		workPointField[i] = float(UlocalPointer[3*i]);
		// Uy
		workPointField[i+nPointsMyrank[0]] = float(UlocalPointer[3*i+1]);
		// Uz
		workPointField[i+2*nPointsMyrank[0]] = float(UlocalPointer[3*i+2]);
	}
	writePointFields(timestep, 0, workPointField);
	// Restore Ulocal array
	ierr = VecRestoreArray(Ulocal, &UlocalPointer); CHKERRQ(ierr);

	// CELL FIELD(S)
	PetscScalar *xpp, *xp;
	VecGetArray(opt->x,&xp);
	VecGetArray(opt->xPhys,&xpp);
	for (unsigned long int i=0; i<nCellsMyrank[0]; i++){
		// Density
		workCellField[i] = float(xp[i]);
		workCellField[i+nCellsMyrank[0]] = float(xpp[i]);
		
	}
	writeCellFields(0, workCellField);
	// Restore arrays
	VecRestoreArray(opt->x,&xp);
	VecRestoreArray(opt->xPhys,&xpp);
	

	// clean up
	ierr = VecDestroy(&Ulocal); CHKERRQ(ierr);
  
	return ierr;
}

void MPIIO::Allocate(std::string info, const int nDom, const int nPFields[], 
					 const int nCFields[], unsigned long int nPointsMyrank[], 
					unsigned long int nCellsMyrank[], unsigned long int nodesPerElement,
					std::string pFNames, std::string cFNames) //, std::string filename)
/* 	info = string with user defined info
	nDom = number of domains
   	nPFields = array with number of point fields in each domain (the number we want to write)
    pFNames = string with point field names
	nCFields = array with number of cell fields in each domain (the number we want to write)
    cFNames = string with cell field names
	nPointsMyrank = array with number of points in each domain in thread=Myrank
	nCellsMyrank = array with number of cells in each domain in thread=Myrank
    nodesPerElement = (max) number of nodes per element. 
	filename = name of output file (default = "output.dat")
*/ 
{
	// default name
	std::string filename = "output.dat";
	
	// Check PETSc input for a work directory
	char filenameChar[PETSC_MAX_PATH_LEN];
	PetscBool flg = PETSC_FALSE;	
	PetscOptionsGetString(NULL,"-workdir",filenameChar,sizeof(filenameChar),&flg);

	// If input, change path of the file in filename
	if (flg){
		filename="";
		filename.append(filenameChar);
		filename.append("/output.dat");
	}

	PetscPrintf(PETSC_COMM_WORLD,"########################################################################\n");
	PetscPrintf(MPI_COMM_WORLD,"Outputfile is written to: %s \n",filename.c_str());
	PetscPrintf(MPI_COMM_WORLD,"To change the working directory, specify '-workdir' at runtime\n");

	// Continue to allocate
	firstFieldOutputDone = false; // Initialize
	this->filename = filename;
	int ierror;
	int headerLen;
	MPI_File fh; // File handle
	// Find how many bytes are used to store an unsigned long integer
	MPI_Type_size(MPI_UNSIGNED_LONG, &MPI_IS);
	// Bytes used to store a float
	MPI_Type_size(MPI_FLOAT, &MPI_FS);
	// Bytes used to store a char
	MPI_Type_size(MPI_CHAR, &MPI_CS);
	// Find out how many cpus we have and their ranks
	ierror = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (ierror) {abort("Problems getting rank", "MPIIO:MPIIO");}
	ierror = MPI_Comm_size(MPI_COMM_WORLD, &ncpu);
	if (ierror) {abort("Problems getting number of cpus", "MPIIO:MPIIO");}

	// Communicate number of points and cells to all processors
	// First allocate space for the arrays
	// Check how many domains your trying to allocate (to avoid crash)
	if (nDom > 1000){abort("ERROR: More than 1000 domains!", "MPIIO:MPIIO");}
	nPoints = new unsigned long int[nDom*ncpu];
	nCells = new unsigned long int[nDom*ncpu];
	for (int i=0; i < nDom; i++){
		ierror = MPI_Allgather(&nPointsMyrank[i], 1, MPI_UNSIGNED_LONG, &nPoints[i*ncpu], 
				1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
		if (ierror) {abort("Problems exchanging number of points", "MPIIO:MPIIO");}
		ierror = MPI_Allgather(&nCellsMyrank[i], 1, MPI_UNSIGNED_LONG, &nCells[i*ncpu], 
				1, MPI_UNSIGNED_LONG, MPI_COMM_WORLD);
		if (ierror) {abort("Problems exchanging number of cells", "MPIIO:MPIIO");}
	}
	// Also allocate and save the other data provided
	this->nDom = nDom;
	this->nodesPerElement = nodesPerElement;
	this->nPFields = new int[nDom];
	this->nCFields = new int[nDom];
	this->nPointsT = new unsigned long int[nDom]; // Total number of points
	this->nCellsT  = new unsigned long int[nDom]; // Total number of cells
	// All processors position in the file is moved below the outputted data
	headerLen = 2+4*nDom; 
	unsigned long int* header = new unsigned long int[headerLen];
	// Put the number of domains into the buffer
	header[0] = nDom;
	for (int i=0; i<nDom ; i++){
		// Sum up total number of points and cells in the domain
		header[1+i] = 0;
		header[1+nDom+i] = 0;
		for (int j=0; j<ncpu; j++){
			header[1+i] += nPoints[i*ncpu + j];
			header[1+nDom+i] += nCells[i*ncpu + j]; 
		}
		nPointsT[i] = header[1+i];
		nCellsT[i]  = header[1+nDom+i];
		// And then the number of point/cell fields in each domain
		header[1+2*nDom+i] = nPFields[i];
		header[1+3*nDom+i] = nCFields[i];
		this->nPFields[i] = nPFields[i];
		this->nCFields[i] = nCFields[i];
	}	
	// Last entry is the nodes per element
	header[headerLen-1] = nodesPerElement;
	// Save the number of characters to output
	int numberOfCharacters = info.size() + pFNames.size() + cFNames.size() + 4;
	// The first processor outputs total number of points, cells, and fields
	if(rank == 0){ // Can this part be done as standard C++ binary output? 
		// If there is an old file, delete this
		ierror = MPI_File_delete(&filename[0], MPI_INFO_NULL);
		// The below line is commented since it caused errors with some MPI implementations
		//if (ierror != _NO_SUCH_FILE && ierror) {abort("Problems deleting old file", "MPIIO:MPIIO");}
		// Then open file
		ierror = MPI_File_open(MPI_COMM_SELF, &filename[0], MPI_MODE_CREATE | MPI_MODE_WRONLY,
				MPI_INFO_NULL, &fh);
		if (ierror) {abort("Problems opening file", "MPIIO::MPIIO");}
		// No need to create filetype here, its just a buffer with integers
		offset = 0; // Start at the beginning of the file 
		// Set view
		ierror = MPI_File_set_view(fh, offset, MPI_CHAR, MPI_CHAR, (char *)"native", MPI_INFO_NULL);
		if (ierror) {abort("Problems setting view", "MPIIO::MPIIO");}
		// Write to the file
		info.append("\n\x01"); // Make sure the string ends with an endline
		ierror = MPI_File_write(fh, (char *)info.c_str(), info.size(), MPI_CHAR, MPI_STATUS_IGNORE);
		if (ierror) {abort("Problems writing to file", "MPIIO::MPIIO");}
		// Set view
		offset += MPI_CS*info.size(); // Adjust offset
		ierror = MPI_File_set_view(fh, offset, MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG, (char *)"native", MPI_INFO_NULL);
		if (ierror) {abort("Problems setting view", "MPIIO::MPIIO");}
		// Write to the file
		ierror = MPI_File_write(fh, header, headerLen, MPI_UNSIGNED_LONG, MPI_STATUS_IGNORE);
		if (ierror) {abort("Problems writing to file", "MPIIO::MPIIO");}
		// Set view
		offset += MPI_IS*headerLen; // Adjust offset
		ierror = MPI_File_set_view(fh, offset, MPI_CHAR, MPI_CHAR, (char *)"native", MPI_INFO_NULL);
		if (ierror) {abort("Problems setting view", "MPIIO::MPIIO");}
		// Write to the file
		pFNames.append("\x01"); // Make sure the string ends with an endline
		cFNames.append("\x01"); // Make sure the string ends with an endline
		pFNames.append(cFNames); // Output both strings at once
		// Write to the file
		ierror = MPI_File_write(fh, (char *)pFNames.c_str(), pFNames.size(), MPI_CHAR, MPI_STATUS_IGNORE);
		// Close the file (I don't think we need a barrier here)
		ierror = MPI_File_close(&fh);
		if (ierror) {abort("Problems closing file", "MPIIO::MPIIO");}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	// All processors position in the file is moved below the outputted data
	offset = MPI_IS*headerLen + MPI_CS*numberOfCharacters;
	// ALWAYS remember to deallocate:
	delete [] header;
}



// Output coordinates - only done once
void MPIIO::writePoints(int domain, float coordinates[])
/*
domain = The domain number (not the name)
coordinates = An array with the coordinates 
*/
{
	int ierror;
	// Open file
	ierror = MPI_File_open(MPI_COMM_WORLD, &filename[0], MPI_MODE_CREATE | 
			MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
	if (ierror) {abort("Problems opening file", "MPIIO::writePoints");}
	// Compute offset for the different processors
	if (domain == 0){ // For the first domain sum up to the current position
		offset += 3*sum(nPoints, rank)*MPI_FS;
	}
	else{ // For the rest of the domains add the points written since last time 
		offset += 3*sum(&nPoints[ncpu*(domain-1)+rank], ncpu)*MPI_FS;
	}
	// Set view
	ierror = MPI_File_set_view(fh, offset, MPI_FLOAT, MPI_FLOAT, (char*)"native", MPI_INFO_NULL);
	if (ierror) {abort("Problems setting view", "MPIIO::writePoints");}
	// Write to the file
	int len = 3*nPoints[domain*ncpu + rank]; // Number of floats to write
	ierror = MPI_File_write_all(fh, coordinates, len, MPI_FLOAT, MPI_STATUS_IGNORE);
	if (ierror) {abort("Problems writing to file", "MPIIO::writePoints");}
	// Close the file (I don't think we need a barrier here)
	ierror = MPI_File_close(&fh);
	if (ierror) {abort("Problems closing file", "MPIIO::writePoints");}

}


// Output cells - only done once
// You provide the elements with "local numbering" - 
// the method will convert to "global numbering".
void MPIIO::writeCells(int domain, unsigned long int elements[],unsigned long int cellsOffset0[], unsigned long int cellsTypes0[])
/*
domain = The domain number
elements = Array with node connectivities. It is assumed that
           	the node numbering is local in both domain and thread.
			And that we want a an overall global numbering for all
			domains.
			Furthermore, every nodesPerElement number should be 
			an element type number.
*/
{
	int ierror;
	unsigned long int shift;
	// Compute the shift number (from local to global node number)
	// This is done by summing up all points outputted before the points in this
	// domain from this rank
	shift = sum(nPoints, ncpu*domain+rank);
	// Run through all "elements" and make them global by adding "shift":
	for (unsigned long int i=0; i<(nodesPerElement)*nCells[ncpu*domain+rank];i++){
		// but shift all the nodes to global numbering
		elements[i] += shift;
	}
	// Open file
	ierror = MPI_File_open(MPI_COMM_WORLD, &filename[0], MPI_MODE_CREATE | 
							MPI_MODE_WRONLY,	MPI_INFO_NULL, &fh);
	if (ierror) {abort("Problems opening file", "MPIIO::writeCells");}

// Compute offset FOR ELEMENT CONN for the different processors
	if (domain == 0){ // For the first domain sum up to the current position
		// First, add for the remaining points written
		offset += 3*sum(&nPoints[ncpu*(nDom-1)+rank], ncpu-rank)*MPI_FS;
		// Then, add the cells
		offset += (nodesPerElement)*sum(nCells, rank)*MPI_IS;
	}
	else{ // For the rest of the domains add the elements written since last time 
		offset += (nodesPerElement)*sum(&nCells[ncpu*(domain-1)+rank], ncpu)*MPI_IS;
	}
	// Set view
	ierror = MPI_File_set_view(fh, offset, MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG,(char*)"native",MPI_INFO_NULL);
	if (ierror) {abort("Problems setting view", "MPIIO::writeCells");}
	// Length of data stream to write
	unsigned long int len = (nodesPerElement)*nCells[domain*ncpu + rank]; // Number of integers
	// Write ELEMENT Conn to file
	ierror = MPI_File_write_all(fh, elements, len, MPI_UNSIGNED_LONG, MPI_STATUS_IGNORE);
	if (ierror) {abort("Problems writing ELEMENTS to file", "MPIIO::writeCells");}

// Write the VTK OFFSET
	// Update the write offset	
	// First jump past ALL the connectivity
	offset += (nodesPerElement)*sum(&nCells[ncpu*(nDom-1)+rank], ncpu-rank)*MPI_IS;
	// Next jump past the previous ranks offset list
	offset += sum(nCells, rank)*MPI_IS;
	unsigned long int addToOffsetList = nodesPerElement*sum(nCells, rank);
	for (int i=0;i<(int)nCells[ncpu*domain+rank];i++){
		cellsOffset0[i]+=addToOffsetList;
	}
	// Length of the offset to write
	len = nCells[domain*ncpu + rank]; // Number of integers
	// Move the offset in the file
	ierror = MPI_File_set_view(fh, offset, MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG,(char*)"native",MPI_INFO_NULL);
	if (ierror) {abort("Problems setting view OFFSET", "MPIIO::writeCells");}
	// write the offset list
        ierror = MPI_File_write_all(fh, cellsOffset0, len, MPI_UNSIGNED_LONG, MPI_STATUS_IGNORE);

// Write the VTK ELEMENT TYPE
	// First jump past ALL the offsets
	offset += sum(&nCells[ncpu*(nDom-1)+rank], ncpu-rank)*MPI_IS;
	// Nextjump past the previous ranks type list
	offset += sum(nCells, rank)*MPI_IS;
	// Length of type list for this rank
	len = nCells[domain*ncpu + rank]; // Number of integers
	// Move the offset in the file
        ierror = MPI_File_set_view(fh, offset, MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG,(char*)"native",MPI_INFO_NULL);
	// write the type list to file
        ierror = MPI_File_write_all(fh, cellsTypes0, len, MPI_UNSIGNED_LONG, MPI_STATUS_IGNORE);

	// Close the file (I don't think we need a barrier here)
	ierror = MPI_File_close(&fh);
	if (ierror) {abort("Problems closing file", "MPIIO::writeCells");}

}

// Output point fields
void MPIIO::writePointFields(unsigned long int timeStep, int domain, float fields[], 
								std::string newFilename)
// timeStep = integer specifying time step
// domain = domain number (integer, should always start at zero and increase by one)
// fields = array containing field values (for all fields, the fields come after each other),
//			which means the user has to put the field values in an array before passing it to
//			this method. This is done because, the entries have to be converted from double
//			to single precision anyway, and then it does not matter that the user
//			has to allocate an extra array and store the single precision values in this.
//			Furthermore, it simplifies the MPI_IO commands.
// newFilename = optional argument with filename to the file you want to write
// NB: It is assumed that all domains are written to the same file
//     and that they are always written in the same chronological order (this could be changed)
{
	int ierror;
	if (newFilename != "notDefined" && newFilename != filename){
		if (domain != 0) {abort("Only new filename when first domain!", "MPIIO::writePointFields");}
		filename = newFilename;
		// Reset positions
	    offset = 0;	
	} 
	// Compute the offset
	else if (domain == 0) {
		if (!firstFieldOutputDone){ 
			// Add for the remaining elements written
			offset += sum(&nCells[ncpu*(nDom-1)+rank], ncpu-rank)*MPI_IS;
		}
		else { // We have outputted fields earlier
			// Add for the remaining cell field values written
			offset += sum(&nCells[ncpu*(nDom-1)+rank], ncpu-rank)*MPI_FS;
		}
	}
	if (domain == 0){ // For the first domain sum up to the current position
		// Add the point field values written by other ranks
		offset += sum(nPoints, rank)*MPI_FS;
	}
	else{// For the rest of the domains add the point field values written since last time 
		// From the previous domain
		offset += sum(&nPoints[ncpu*(domain-1)+rank], ncpu-rank)*MPI_FS;
		// From the current domain (if rank != 0)
		offset += sum(&nPoints[ncpu*domain], rank)*MPI_FS;
	}
	// If it is the first domain (domain=0), the timeStep should be written by rank 0 and all 
	// offsets should be increased by one
	if (domain==0){
		if (rank==0){	
		//abort(" --------- HERE---------------", "MPIIO::writePointFields");
			ierror = MPI_File_open(MPI_COMM_SELF, &filename[0], MPI_MODE_CREATE | 
							MPI_MODE_WRONLY,	MPI_INFO_NULL, &fh);
			if (ierror) {abort("Problems opening file", "MPIIO::writePointFields");}
			// Set view
			ierror = MPI_File_set_view(fh, offset, MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG, (char *)"native", MPI_INFO_NULL);
			if (ierror) {abort("Problems setting view", "MPIIO::writePointFields");}
			// Write to the file
			ierror = MPI_File_write(fh, &timeStep, 1, MPI_UNSIGNED_LONG, MPI_STATUS_IGNORE);
			if (ierror) {abort("Problems writing to file", "MPIIO::writePointFields");}
			// Close the file
			ierror = MPI_File_close(&fh);
			if (ierror) {abort("Problems closing file", "MPIIO::writePointFields");}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		// Increase all offsets by one integer:
		offset += 1*MPI_IS;
	}
	//MPI_Barrier(MPI_COMM_WORLD);
	// Open the file
	ierror = MPI_File_open(MPI_COMM_WORLD, &filename[0], MPI_MODE_CREATE | 
							MPI_MODE_WRONLY,	MPI_INFO_NULL, &fh);
	if (ierror) {abort("Problems opening file", "MPIIO::writePointFields");}
	// Set filetype
	MPI_Datatype filetype;
	// The number of points for the specified rank in this domain
	int blocklength = nPoints[ncpu*domain + rank];
	// The total number of points in the domain
	int stride = nPointsT[domain];
	// Number of blocks to write
	int count = nPFields[domain];
	ierror = MPI_Type_vector(count, blocklength, stride, MPI_FLOAT, &filetype);
	if (ierror) {abort("Problems creating MPI vector", "MPIIO::writePointFields");}
	ierror = MPI_Type_commit(&filetype);
	if (ierror) {abort("Problems creating filetype", "MPIIO::writePointFields");}
	// Set view
	ierror = MPI_File_set_view(fh, offset, MPI_FLOAT, filetype, (char *)"native", MPI_INFO_NULL);
	if (ierror) {abort("Problems setting view", "MPIIO::writePointFields");}
	
	// Since the array fields contain all fields, we can output them all at once
	// Remember: datatype specifies the layout in memory, while filetype specifies 
	// the layout in the file. They are both of the MPI_Datatype kind.
	// Set datatype such that the access in memory is correct. In this case the 
	// field values are in the right order already, so no need to make a datatype.
	// Write to the file (one datatype is written)
	ierror = MPI_File_write_all(fh, fields, count*blocklength, MPI_FLOAT, MPI_STATUS_IGNORE);  
	if (ierror) {abort("Problems writing to file", "MPIIO::writePointFields");}
	// Close the file
	ierror = MPI_File_close(&fh);
	if (ierror) {abort("Problems closing file", "MPIIO::writePointFields");}
	// Check if it was the first time this function has been called 
	if (!firstFieldOutputDone){firstFieldOutputDone = true;}
	// Free the memory used for filetype
	ierror = MPI_Type_free(&filetype);
	if (ierror) {abort("Problems freeing datatype", "MPIIO::writePointFields");}
	// Finally, update the offset to the beginning of the last field we wrote
	offset += stride*(count-1)*MPI_FS;
}

// Output cell fields
void MPIIO::writeCellFields(int domain, float fields[])
// timeStep = integer specifying time step
// domain = domain number
// fields = field values at cell points. Should contain all fields!
// NB: It is always assumed that cell fields are written to the same file as point fields
//     and that all domains are written to the same file, and that point fields are called first
{
	int ierror;
	if (domain == 0){ 
		// First, add for the remaining points written
		offset += sum(&nPoints[ncpu*(nDom-1)+rank], ncpu-rank)*MPI_FS;
		// Add the cell field values written by other ranks
		offset += sum(nCells, rank)*MPI_FS;
	}
	else{// For the rest of the domains add the cell field values written since last time 
		// From the previous domain
		offset += sum(&nCells[ncpu*(domain-1)+rank], ncpu-rank)*MPI_FS;
		// From the current domain (if rank != 0)
		offset += sum(&nCells[ncpu*domain], rank)*MPI_FS;
	}
	// Open the file
	ierror = MPI_File_open(MPI_COMM_WORLD, &filename[0], MPI_MODE_CREATE | 
							MPI_MODE_WRONLY,	MPI_INFO_NULL, &fh);
	if (ierror) {abort("Problems opening file", "MPIIO::writeCellFields");}
	
	// Set filetype
	MPI_Datatype filetype;
	// The number of points for the specified rank in this domain
	int blocklength = nCells[ncpu*domain + rank];
	// The total number of points in the domain
	unsigned long int stride = nCellsT[domain];
	// Number of blocks to write
	int count = nCFields[domain];
	ierror = MPI_Type_vector(count, blocklength, stride, MPI_FLOAT, &filetype);
	if (ierror) {abort("Problems creating MPI vector", "MPIIO::writeCellFields");}
	ierror = MPI_Type_commit(&filetype);
	if (ierror) {abort("Problems creating filetype", "MPIIO::writeCellFields");}
	// Set view
	ierror = MPI_File_set_view(fh, offset, MPI_FLOAT, filetype, (char *)"native", MPI_INFO_NULL);
	if (ierror) {abort("Problems setting view", "MPIIO::writeCellFields");}
	
	// Since the array fields contain all fields, we can output them all at once
	// Remember: datatype specifies the layout in memory, while filetype specifies 
	// the layout in the file. They are both of the MPI_Datatype kind.
	// Set datatype such that the access in memory is correct. In this case the 
	// field values are in the right order already, so no need to make a datatype.
	// Write to the file (one datatype is written)
	ierror = MPI_File_write_all(fh, fields, count*blocklength, MPI_FLOAT, MPI_STATUS_IGNORE);  
	if (ierror) {abort("Problems writing to file", "MPIIO::writeCellFields");}
	// Close the file
	ierror = MPI_File_close(&fh);
	if (ierror) {abort("Problems closing file", "MPIIO::writeCellFields");}
	// Free the memory used for filetype
	ierror = MPI_Type_free(&filetype);
	if (ierror) {abort("Problems freeing datatype", "MPIIO::writePointFields");}
	// Finally, update the offset to the beginning of the last field we wrote
	offset += stride*(count-1)*MPI_FS;
}

// Method to do MPI errors
void MPIIO::abort(std::string errorMsg, std::string position)
// errorMsg = the error message the programmer has written
// position = the methode in which the error occured
{
	std::cerr << errorMsg << " in " << position << std::endl;
	// Stop the execution of the program
	MPI_Barrier(MPI_COMM_WORLD);
	std::cerr << "rank = " << rank << std::endl;
	MPI_Abort(MPI_COMM_WORLD, -1);
	exit(0);
	
}

unsigned long int MPIIO::sum(unsigned long int *startPos, unsigned long int nel)
{
	unsigned long int total=0; // The number of points
	for (unsigned long int i = 0; i < nel; i++){
		total += startPos[i]; 
	}	
	return total;
}

PetscErrorCode MPIIO::DMDAGetElements_3D(DM dm,PetscInt *nel,PetscInt *nen,const PetscInt *e[]) {
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
