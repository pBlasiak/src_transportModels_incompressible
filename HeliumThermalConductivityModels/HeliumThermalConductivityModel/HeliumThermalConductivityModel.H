/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::HeliumThermalConductivityModel

Description

SourceFiles
    HeliumThermalConductivityModel.C
    newThermalProperty.C

\*---------------------------------------------------------------------------*/

#ifndef HeliumThermalConductivityModel_H
#define HeliumThermalConductivityModel_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "HeliumModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class thermalIncompressibleTwoPhaseMixture;

/*---------------------------------------------------------------------------*\
                Class HeliumThermalConductivityModel Declaration
\*---------------------------------------------------------------------------*/

class HeliumThermalConductivityModel
: 
	IOdictionary
{

protected:

    // Protected data

		//- Access to dictionary
        dictionary HeliumPropertiesDict_;

		//- Reference to velocity field
        const volVectorField& U_;

        const surfaceScalarField& phi_;

		//- Thermal conductivity field
		volScalarField kHe_;


    // Private Member Functions

        //- Disallow copy construct
        HeliumThermalConductivityModel(const HeliumThermalConductivityModel&);

        //- Disallow default bitwise assignment
        void operator=(const HeliumThermalConductivityModel&);

public:

    //- Runtime type information
    TypeName("HeliumThermalConductivityModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            HeliumThermalConductivityModel,
            components,
            (
                const volVectorField& U,
                const surfaceScalarField& phi
            ),
            (U, phi)
        );


    // Selectors

        //- Return a reference to the selected thermal conductivity model
        static autoPtr<HeliumThermalConductivityModel> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    // Constructors

        //- Construct from components
        HeliumThermalConductivityModel
        (
            const word& type,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~HeliumThermalConductivityModel()
    {}


    // Member Functions
	
		//- Sets thermal conductivity field
		void setkHe(const volScalarField& k)
		{
			kHe_ = k;	
		}

        //- Return the helium properties dictionary
        const dictionary& HeliumPropertiesDict() const
        {
            return HeliumPropertiesDict_;
        }

		//- Returns reference to squared magnitude 
		//  of gradient temperature
		const volScalarField& magGradT2() const
		{
			return U_.db().lookupObject<volScalarField>("magGradT2");
		}

        //- Return the helium thermal conductivity 
        const volScalarField& kHe() const
        {
            return kHe_;
        }

        //- Correct the thermal conductivity of He II
        virtual void calckHe(const HeliumModel& hm) = 0;

        //- Read transportProperties dictionary
        //virtual bool read(const dictionary& HeliumProperties) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
