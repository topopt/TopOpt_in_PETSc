#ifndef __FILTER__
#define __FILTER__

#include <petsc.h>
//#include <petsc-private/dmdaimpl.h>
#include <petsc/private/dmdaimpl.h>
#include <iostream>
#include <math.h>
#include <PDEFilter.h>
#include <TopOpt.h>

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


class Filter{
  
public:
    // Constructor
    Filter(TopOpt *opt); 
    
    // Destructor
    ~Filter(); 
    
    // Filter design variables
    PetscErrorCode FilterProject(TopOpt *opt);
    
    // Filter the sensitivities
    PetscErrorCode Gradients(TopOpt *opt);
    
private:
  
    // Standard density/sensitivity filter matrix
    Mat H; 		// Filter matrix
    Vec Hs; 		// Filter "sum weight" (normalization factor) vector   
    
    // Mesh used for standard filtering
    DM da_elem;  	// da for image-filter field mesh
    
    // PDE filtering
    PDEFilt *pdef;	// PDE filter class
    
    // Setup datastructures for the filter
    PetscErrorCode SetUp(TopOpt *opt);
    
    // Routine that doesn't change the element type upon repeated calls
    PetscErrorCode DMDAGetElements_3D(DM dm,PetscInt *nel,PetscInt *nen,const PetscInt *e[]);
    
};

#endif
