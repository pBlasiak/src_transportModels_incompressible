/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Class
    Foam::singlePhaseHeliumTransportModel

Description
    A simple single-phase transport model based on HeliumModel.

    Used by the incompressible single-phase solver buoyantBoussinesqSuperFluidPimpleFoam.

SourceFiles
    singlePhaseHeliumTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef singlePhaseHeliumTransportModel_H
#define singlePhaseHeliumTransportModel_H

//#include "incompressible/transportModel/transportModel.H"
#include "transportModel.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "HeliumModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//class HeliumModel;

/*---------------------------------------------------------------------------*\
                Class singlePhaseHeliumTransportModel Declaration
\*---------------------------------------------------------------------------*/

class singlePhaseHeliumTransportModel
:
    public IOdictionary,
    public transportModel
{

    // Private Member Functions

        //- Disallow copy construct
        singlePhaseHeliumTransportModel(const singlePhaseHeliumTransportModel&);

        //- Disallow default bitwise assignment
        void operator=(const singlePhaseHeliumTransportModel&);

protected:		
    // Protected Data

        autoPtr<HeliumModel> HeliumModelPtr_;



public:

    // Constructors

        //- Construct from components
        singlePhaseHeliumTransportModel
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~singlePhaseHeliumTransportModel();


    // Member Functions

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const
		{
		    return HeliumModelPtr_->nu();
		}

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const
		{
		    return HeliumModelPtr_->nu(patchi);
		}

		const volScalarField& betaHe() const
		{
			return HeliumModelPtr_->beta();
		}

		const volScalarField& onebyf() const
		{
			return HeliumModelPtr_->onebyf();
		}

		const volScalarField& cpHe() const
		{
			return HeliumModelPtr_->cp();
		}

		const volScalarField& sHe() const
		{
			return HeliumModelPtr_->s();
		}

		const volScalarField& etaHe() const
		{
			return HeliumModelPtr_->eta();
		}

		const volScalarField& rhoHe() const
		{
			return HeliumModelPtr_->rho();
		}

		const volScalarField& rhon() const
		{
			return HeliumModelPtr_->rhon();
		}

		const volScalarField& rhos() const
		{
			return HeliumModelPtr_->rhos();
		}

		const volScalarField& AGMHe() const
		{
			return HeliumModelPtr_->AGM();
		}

        //- Correct the laminar viscosity
        virtual void correct()
		{
		    HeliumModelPtr_->correct();
		}

        //- Read transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
