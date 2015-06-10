#ifndef MPIIO_H
#define MPIIO_H

// Include necessary libraries
#include <petsc.h>
//#include <petsc-private/dmdaimpl.h>
#include <petscdmda.h>
#include <string>
#include <mpi.h>
#include <TopOpt.h>

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


class MPIIO
{
	public:
		// ------------- METHODS ------------------------------------------
		
		MPIIO(DM da_nodes,int nPfields, std::string pnames,int nCfields, std::string cnames);
		~MPIIO();

		// NOT CLEAN INTERFACE: REPLACE BY STD::PAIR OR SUCH !!!!!!
		PetscErrorCode WriteVTK(DM da_nodes, Vec U, TopOpt *opt, PetscInt itr);
		
	private:
		// -------------- METHODS -----------------------------------------
		
		
		void abort(std::string errorMsg, std::string position);
		
		unsigned long int sum(unsigned long int *startPos, unsigned long int nel);
		// --------------- MEMBERS ----------------------------------------
		int MPI_IS;			//!< The size of an MPI unsigned long integer in bytes 
		int MPI_FS;			//!< The size of an MPI float in bytes
		int MPI_CS;			//!< The size of an MPI char in bytes
		int nDom;			//!< Number of domains
		int *nPFields;		//!< Number of point fields in each domain
		int *nCFields;		//!< Number of cell/element fields in each domain
		MPI_Offset offset;	//!< The offset of each thread in the file
		int rank;			//!< The processor rank
		int ncpu;			//!< The number of cpus
		int nodesPerElement;//!< Number of nodes per element
		bool firstFieldOutputDone; //!< Will be set to true after first field output
		std::string filename;//!< Output filename
		MPI_File fh; 		 //!< Filehandle
		
		void Allocate(std::string info, const int nDom, const int nPFields[],
				const int nCFields[], unsigned long int nPointsMyrank[], 
				unsigned long int nCellsMyrank[], unsigned long int nodesPerElement,
				std::string pFNames, std::string cFNames);
				//std::string filename = "/home/naage/PETSc/output2.dat");
		
		void writePoints(int domain, float coordinates[]);
		
		void writeCells(int domain, unsigned long int elements[],unsigned long int cellsOffset0[], unsigned long int cellsTypes0[]);

		void writePointFields(unsigned long int timeStep, int domain, float fields[],
				std::string newFilename = "notDefined" );	

		void writeCellFields(int domain, float fields[]);	
		// ------------ MEMBERS  - maybe they can be private too -----------
		unsigned long int *nPoints;	//!< The number of points in each domain in each thread
		unsigned long int *nCells;	//!< The number of elements/cells in each domain in each thread
		unsigned long int *nPointsT;//!< The total number of points in each domain
		unsigned long int *nCellsT; //!< The total number of cells in each domain
		
		// Converters needed for PETSc adaptation
		unsigned long int *nPointsMyrank, *nCellsMyrank;
		float *workPointField, *workCellField;
		PetscErrorCode DMDAGetElements_3D(DM dm,PetscInt *nel,PetscInt *nen,const PetscInt *e[]);
		
		
};
/** @example 
	An illustrative example explaining how to use the class
	@code 
  string userDefined = "something"
  int nDom 		 = 2;
  int nPFields[nDom]= {2, 1};
  string pFieldNames = "density, pdensity, density";
  int nCFields[nDom]= {1, 1};
  string cFieldNames = "density, density";
  int nPointsMyrank[nDom] = {10, 20}; // Arbitrary numbers 
  int nCellsMyrank[nDom] = {9, 19};  
  
  // Initialize
  MPIIO outputObject = DFEMMMPIIO(userDefined, nDom, nPFields, pFieldNames, 
   								nCFields, cFieldNames, rank, nPoints, nCells);
  // Both the total number of points/cells and the corresponding
  // processor dependent numbers are determined in the constructor!
  // This does of course require MPI communication, but is only done once.
  
   
  // Output coordinates/points/vertices:
  // First domain:
  float points1[3*nPoints[0]]; 
  // ... put coordinates into the array points1
  outputObject.writePoints(0, points1);
  // Second domain:
  float points1[3*nPoints[1]]; 
  // ... put coordinates into the array points2
  outputObject.writePoints(1, points2);
 
  // Output elements:	
  // First domain:
  int elements1[USER_SPECIFIED_SIZE]; 
  // ... put elements into the array elements1
  outputObject.writeCells(0, elements1);
  // Second domain
  int elements2[USER_SPECIFIED_SIZE];
  // ... put elements into the array
  outputObject.writeCells(1,  elements2);

  // The thought is that elements should be a long list with integers, where one integer
  // specifies the element type and the following integers specify point connectivity.
  
  // ...  Perform your computations
  
  // Output fields at a time point (timeStep)
  // First point fields are outputted, then cell fields
  
  // Output point field:
  // First domain:
  // Put all field variables into an array pFieldsD1, should of course 
  // be in the same order as your point list, and the fields should come
  // after each other.
  outputObject.writePointField(timeStep, 0, pFieldsD1); 
  // Second domain:
  // Put all field variables into an array pFieldsD2.
  outputObject.writePointField(1, pFieldsD2); 
  
  // Output cell/element fields:
  // First domain:
  // Put all field variables into an array cFieldsD1, should of course 
  // be in the same order as your element list, and the fields should come after each other.
  outputObject.writePointField(0, cFieldsD1); 
  // Second domain:
  // Put all field variables into an array cFieldsD2.
  outputObject.writePointField(1, cFieldsD2); 
   
  // Repeat as many times as you want.
	@endcode
*/

   
#endif

