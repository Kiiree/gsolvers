/*---------------------------------------------------------------------------*\
	This file is part of gSolvers
	
	gSolvers is based on cusp-library(http://code.google.com/p/cusp-library/) 
	and provides GPU solver interface for OpenFOAM(http://www.openfoam.com/).

    gSolvers is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    gSolvers is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with gSolvers.  If not, see <http://www.gnu.org/licenses/>.
    
\*---------------------------------------------------------------------------*/

#ifndef GSOLVERSMATRIX_H
#define GSOLVERSMATRIX_H

#include "gSolvers.h"
#include <cusp/csr_matrix.h>
#include <iostream>

template<typename RealType, typename MemorySpace>
void ldu2csr(
		Index const rows,
		Index const nonZeroNumber,
		Index const faces,//matrix.upper().size()
		Index const* upperAddr,//matrix.lduAddr().upperAddr()
		Index const* lowerAddr,//matrix.lduAddr().lowerAddr()
		RealType const* upper,//matrix.upper()
		RealType const* lower,//matrix.lower()
		RealType const* diag,//matrix.diag()
		cusp::csr_matrix<Index, RealType, MemorySpace> &csrmatrix)
{

	cusp::array1d<Index, cusp::host_memory> row_offsets(rows+1,1);
	
	cusp::array1d<RealType, cusp::host_memory> v(nonZeroNumber);
	
	cusp::array1d<Index, cusp::host_memory> c_idx(nonZeroNumber);
	
	csrmatrix.resize(rows,rows,nonZeroNumber);

	for(Index i=0; i < faces; i++)
	{
		Index lr = lowerAddr[i];
		Index ur = upperAddr[i];
		row_offsets[lr+1]++;
		row_offsets[ur+1]++;
	}
	row_offsets[0] = 0;
	
	for(Index i = 0;i<rows; i++)
	{
		row_offsets[i+1]+= row_offsets[i];	
	}

	csrmatrix.row_offsets = row_offsets;

	//L
	for(Index i = 0; i < faces; i++)
	{
		Index row = upperAddr[i];
		Index col = lowerAddr[i];

		Index offsets = row_offsets[row]++;
		v[offsets] = lower[i];
	    c_idx[offsets] = col;	
	}

	//D
	for(Index i = 0; i < rows; i++)
	{
		Index offsets = row_offsets[i]++;
		v[offsets] = diag[i];
		c_idx[offsets] = i;	
	}

	//U
	for(Index i = 0; i < faces; i++)
	{
		Index row = lowerAddr[i];
		Index col = upperAddr[i];

		Index offsets = row_offsets[row]++;
		v[offsets] = upper[i];
		c_idx[offsets] = col;
	}
	
	csrmatrix.values = v;
	csrmatrix.column_indices = c_idx;

};    


#endif
