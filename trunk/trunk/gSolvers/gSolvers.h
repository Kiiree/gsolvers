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

#ifndef GSOLVERS_H
#define GSOLVERS_H

#ifdef WM_SP
typedef float RealType;
#else
typedef double RealType;
#endif
typedef int Index;
		
void pcgsolve(Index const rows, Index const nonZeroNumber, 
			RealType *hx,
			RealType const *hb,
			Index const faces,//matrix.upper().size()
			Index const* lduUppAdd,//matrix.lduAddr().upperAddr()
			Index const* lduLowAdd,//matrix.lduAddr().lowerAddr()
			RealType const* upperv,//matrix.upper()
			RealType const* lowerv,//matrix.lower()
			RealType const* diagv,//matrix.diag()
			Index const maxIter,
			RealType const relTol,
			RealType const tolerance,
			RealType &initialResidual,
			RealType &finalResidual,
			Index &nIterations,
			bool &converged
			);
		
void pbicgsolve(Index const rows, Index const nonZeroNumber, 
			RealType *hx,
			RealType const *hb,
			Index const faces,//matrix.upper().size()
			Index const* lduUppAdd,//matrix.lduAddr().upperAddr()
			Index const* lduLowAdd,//matrix.lduAddr().lowerAddr()
			RealType const* upperv,//matrix.upper()
			RealType const* lowerv,//matrix.lower()
			RealType const* diagv,//matrix.diag()
			Index const maxIter,
			RealType const relTol,
			RealType const tolerance,
			RealType &initialResidual,
			RealType &finalResidual,
			Index &nIterations,
			bool &converged
			);
#endif
