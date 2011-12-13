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

#include <cusp/krylov/bicgstab.h>

#include <cusp/precond/diagonal.h>

#include <cusp/csr_matrix.h>

#include <cusp/dia_matrix.h>

#include "gmonitor.h"
#include "gSolversMatrix.h"
#include "gSolvers.h"

//#include <cusp/io/matrix_market.h>

void pbicgsolve(Index const rows, Index const nonZeroNumber, 
			RealType *hx,
			RealType const *hb,
			Index const faces,//matrix.upper().size()
			Index const* upperAddr,//matrix.lduAddr().upperAddr()
			Index const* lowerAddr,//matrix.lduAddr().lowerAddr()
			RealType const* upper,//matrix.upper()
			RealType const* lower,//matrix.lower()
			RealType const* diag,//matrix.diag()
			Index const maxIter,
			RealType const relTol,
			RealType const tolerance,
			RealType &initialResidual,
			RealType &finalResidual,
			Index &nIterations,
			bool &converged
			)
{
	cusp::array1d<RealType,cusp::device_memory> x(hx, hx+rows);
	
	cusp::array1d<RealType,cusp::device_memory> b(hb, hb+rows);
	
	cusp::csr_matrix<Index, RealType, cusp::device_memory> A;
	 
	ldu2csr(rows, nonZeroNumber, faces, upperAddr, lowerAddr, upper, lower, diag, A);
    
    //More
    cusp::precond::diagonal<RealType, cusp::device_memory> M(A);

	gSolvers::GMonitor<RealType> monitor(b, maxIter,  relTol, tolerance);
    
    cusp::krylov::bicgstab(A, x, b, monitor, M);

    thrust::copy(x.begin(), x.end(), hx);
    
    initialResidual= monitor.initialResidual();
    
    finalResidual = monitor.finalResidual();
    
    nIterations= monitor.iteration_count();
    
    converged = monitor.converged();
}
