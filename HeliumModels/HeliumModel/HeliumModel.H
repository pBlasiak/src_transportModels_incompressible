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

Namespace
    Foam::HeliumModels

Description
    A namespace for various incompressible HeliumModel implementations.

Class
    Foam::HeliumModel

Description
    An abstract base class for incompressible HeliumModels.

SourceFiles
    HeliumModel.C
    HeliumModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef HeliumModel_H
#define HeliumModel_H

#include "dictionary.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "runTimeSelectionTables.H"
#include "typeInfo.H"
#include "volFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class HeliumModel Declaration
\*---------------------------------------------------------------------------*/

class HeliumModel
{
    // Private data

protected:

    // Protected data
	
		//- Maximum index number for thermophysical tables i = 667
		static const label indexMax_;

		//- Step of temperature in the thermophysical tables dT = 0.001 K
		static const dimensionedScalar dT_;
	
		//- lambda Temperature 2.1711132461
		static const dimensionedScalar Tlambda_;

		//- Minimum temperature TMin = 1.5 K
		static const dimensionedScalar TMin_;

		//- Maximum temperature TMax = 2.167 K
		static const dimensionedScalar TMax_;

		//- Minimum index number for thermophysical tables i = 0
		static const label indexMin_;

		//- Table for betaHe
		static const List<scalar> betaHeTable_;

		//- Table for AGMHe
		static const List<scalar> AGMHeTable_;

		//- Table for sHe
		static const List<scalar> sHeTable_;

		//- Table for etaHe
		static const List<scalar> etaHeTable_;

		//- Table for cpHe
		static const List<scalar> cpHeTable_;

		//- Table for onebyf
		static const List<scalar> onebyfTable_;

		//- Table for rhoHe
		static const List<scalar> rhoHeTable_;

        word name_;
        dictionary HeliumProperties_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;

		volScalarField TMinField_;
		volScalarField TMaxField_;

        //- Thermal expansion coefficient
        volScalarField betaHe_;
		
        //- Goerter-Mellink parameter
        volScalarField AGMHe_;

        //- Entropy
        volScalarField sHe_;
		
        //- Dynamic Helium
        volScalarField etaHe_;
		
        //- Specific heat at constant pressure
        volScalarField cpHe_;
		
        //- Thermal conductivity function
        volScalarField onebyf_;

        //- Density of the mixture 
        volScalarField rhoHe_;

        //- Kinematic Helium
        volScalarField nu_;

		//- Calculates helium thermal property based on temperature
        void calcHeProp
		(
		    volScalarField& vsf, 
			const List<scalar>& vsfTable,
			const volScalarField& T,
			const label maxIndex = indexMax_,
			const dimensionedScalar dt = dT_
		);

    // Private Member Functions

        //- Disallow copy construct
        HeliumModel(const HeliumModel&);

        //- Disallow default bitwise assignment
        void operator=(const HeliumModel&);


public:

    //- Runtime type information
    TypeName("HeliumModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            HeliumModel,
            dictionary,
            (
                const word& name,
                const dictionary& HeliumProperties,
                const volVectorField& U,
                const surfaceScalarField& phi
            ),
            (name, HeliumProperties, U, phi)
        );


    // Selectors

        //- Return a reference to the selected Helium model
        static autoPtr<HeliumModel> New
        (
            const word& name,
            const dictionary& HeliumProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Constructors

        //- Construct from components
        HeliumModel
        (
            const word& name,
            const dictionary& HeliumProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~HeliumModel()
    {}


    // Member Functions

        //- Return the helium properties dictionary
        const dictionary& HeliumProperties() const
        {
            return HeliumProperties_;
        }

        //- Return the laminar kinematic viscosity of Helium
        virtual tmp<volScalarField> nu() const
        {
            return nu_;
        }

        //- Return the laminar kinematic viscosity of Helium for patch
        virtual tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

		//- Returns thermal expansion coefficient 
        virtual const volScalarField& beta() const
		{
			return betaHe_;
		}

		//- Returns Goerter-Mellink coefficient 
        virtual const volScalarField& AGM() const
		{
			return AGMHe_;
		}

		//- Returns entropy 
        virtual const volScalarField& s() const
		{
			return sHe_;
		}

		//- Returns dynamic Helium 
        virtual const volScalarField& eta() const
		{
			return etaHe_;
		}
	
		//- Returns specific heat capacity 
        virtual const volScalarField& cp() const
		{
			return cpHe_;
		}

		//- Returns reciprocal of thermal conductivity function for helium based on temperature
        virtual const volScalarField& kappa() const
		{
			return onebyf_;
		}
		
        //- Returns the density
        virtual const volScalarField& rho() const
		{
			return rhoHe_;
		}

		static const dimensionedScalar Tlambda() 
		{
			return Tlambda_;
		}
		
        //- Correct the laminar Helium
        virtual void correct() = 0;

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& HeliumProperties) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
