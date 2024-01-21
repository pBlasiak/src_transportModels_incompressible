/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::HeliumThermalConductivityModel>
Foam::HeliumThermalConductivityModel::New
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
{
    IOdictionary heliumThermalConductivityDict
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const word HeliumThermalConductivityModelTypeName
    (
        heliumThermalConductivityDict.lookup("kHeModel")
    );

    Info<< "Selecting helium thermal conductivity model "
        << HeliumThermalConductivityModelTypeName << endl;

	auto* cstrPtr = componentsConstructorTable(HeliumThermalConductivityModelTypeName);
	
    //componentsConstructorTable::iterator cstrIter =
    //    componentsConstructorTablePtr_
    //        ->find(HeliumThermalConductivityModelTypeName);

   // if (cstrIter == componentsConstructorTablePtr_->end())
    if (!cstrPtr)
    {
        FatalErrorInFunction
            << "Unknown HeliumThermalConductivityModel type "
            << HeliumThermalConductivityModelTypeName << endl << endl
            << "Valid  HeliumThermalConductivityModels are : " << endl
            << componentsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<HeliumThermalConductivityModel>(cstrPtr(U, phi));
}


// ************************************************************************* //
