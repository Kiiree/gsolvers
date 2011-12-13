/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "gPBiCG.H"

#include "gSolvers.h"

namespace Foam
{
	defineTypeNameAndDebug(gPBiCG,0);

	lduMatrix::solver::addasymMatrixConstructorToTable<gPBiCG>
		addgPBiCGAsymMatrixConstructorToTable_;

}

Foam::gPBiCG::gPBiCG
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::solver
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    )
{}

Foam::lduMatrix::solverPerformance Foam::gPBiCG::solve
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt
) const
{
    scalar initialResidual;
    scalar finalResidual;
    int nIterations = 0;
    bool converged = false;

    label rows = matrix_.diag().size();
    label nonZeroNumber= matrix_.lower().size() + matrix_.upper().size() + matrix_.diag().size();
		
	pbicgsolve(rows, nonZeroNumber, 
		psi.data(),
		source.cdata(),
		matrix_.upper().size(),
		matrix_.lduAddr().upperAddr().cdata(),
		matrix_.lduAddr().lowerAddr().cdata(),
		matrix_.upper().cdata(),
		matrix_.lower().cdata(),
		matrix_.diag().cdata(),
		maxIter_,
		relTol_,
		tolerance_,
		initialResidual,
		finalResidual,
		nIterations,
		converged
		);

	// --- Setup class containing solver performance data
    lduMatrix::solverPerformance solverPerf
    (
        lduMatrix::preconditioner::getName(controlDict_) + typeName,
        fieldName_,
        initialResidual,
        finalResidual,
        nIterations,
        converged
    );

	return solverPerf;

}
