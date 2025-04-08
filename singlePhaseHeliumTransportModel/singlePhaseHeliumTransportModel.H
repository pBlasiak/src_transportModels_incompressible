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
#include "HeliumThermalConductivityModel.H"

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

		//- Pointer to helium model
        autoPtr<HeliumModel> HeliumModelPtr_;

		//- Pointer to helium thermal conductivity model
        autoPtr<HeliumThermalConductivityModel> kHeModelPtr_; 

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

		void calckHe() 
		{
			kHeModelPtr_->calckHe(*HeliumModelPtr_);
		}

		const volScalarField& kHe() const 
		{
			return kHeModelPtr_->kHe();
		}

		//- Calculates ds/dT - needed in two-fluid model
		volScalarField& dSdT(const word& p) 
		{
			return HeliumModelPtr_->dSdT(p);
		}

		//- Calculates drho/dT - needed in two-fluid model
		volScalarField& drhodT(const word& p) 
		{
			return HeliumModelPtr_->drhodT(p);
		}

		//- Calculates thermal property of HeII at specific 
		//  temperature and pressure
		dimensionedScalar heThermProp
		(
			const dimensionedScalar& temp, 
			const word& wtp,
			const word& wp = "SVP"
		) const
		{
			return HeliumModelPtr_->thermProp(temp,wtp,wp);
		}

        //- Read transportProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
