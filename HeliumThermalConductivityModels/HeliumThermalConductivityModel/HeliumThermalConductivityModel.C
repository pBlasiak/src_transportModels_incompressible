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

\*---------------------------------------------------------------------------*/

#include "HeliumThermalConductivityModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(HeliumThermalConductivityModel, 0);
    defineRunTimeSelectionTable(HeliumThermalConductivityModel, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HeliumThermalConductivityModel::HeliumThermalConductivityModel
(
    const word& type,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
	IOdictionary
	(
	    IOobject
        (
            "transportProperties", // dictionary name
            U.time().constant(),     // dict is found in "constant"
            U.db(),                  // registry for the dict
            IOobject::MUST_READ,     // must exist, otherwise failure
            IOobject::NO_WRITE       // dict is only read by the solver
        )
	),
    HeliumPropertiesDict_
	(
	    IOdictionary
	    (
	        IOobject
            (
                "transportProperties", // dictionary name
                U.time().constant(),     // dict is found in "constant"
                U.db(),                  // registry for the dict
                IOobject::MUST_READ,     // must exist, otherwise failure
                IOobject::NO_WRITE       // dict is only read by the solver
            )
	    )
	),
    U_(U),
    phi_(phi),
	kHe_
    {
        IOobject
        {
            "kHe",
            U_.mesh().time().timeName(),
            U_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
		},
        U_.mesh(),
		dimensionedScalar{dimPower/dimLength/dimTemperature, 0.0}
	} 
{ }


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//bool Foam::HeliumThermalConductivityModel::read(const dictionary& HeliumProperties)
//{
//    HeliumPropertiesDict_ = HeliumProperties;
//
//    return true;
//}



// ************************************************************************* //
