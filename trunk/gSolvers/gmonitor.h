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

#ifndef GMONITOR_H
#define GMONITOR_H

#include <cusp/monitor.h>

#include "gSolvers.h"

namespace gSolvers
{

template <typename Real>
class GMonitor : public cusp::default_monitor<Real>
{
	typedef cusp::default_monitor<Real> super;

public:
	template <typename Vector>
	GMonitor(const Vector &b, size_t maxIterations = 500, Real relTol = 1e-5, Real absTol = 0)
		: super(b, maxIterations, relTol, absTol)
	{}
	
	template <typename Vector>
    bool finished(const Vector& r)
    {   
    	super::r_norm = cusp::blas::nrm2(r);
    	
		if (0 == super::iteration_count()) {
			initialResidual_ = super::r_norm;
		}

		return super::converged() || super::iteration_count() >= super::iteration_limit();
	}
	

	Real initialResidual() const{ return initialResidual_; }
    
    Real finalResidual() const{ return super::r_norm; }

private:

	Real initialResidual_;
	
	Real finalResidual_;
	
};

}

#endif
