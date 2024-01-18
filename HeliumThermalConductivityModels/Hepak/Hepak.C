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

\*---------------------------------------------------------------------------*/

#include "Hepak.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace HeliumThermalConductivityModels
{
    defineTypeNameAndDebug(Hepak, 0);
    addToRunTimeSelectionTable(HeliumThermalConductivityModel, Hepak, components);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HeliumThermalConductivityModels::Hepak::Hepak
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    HeliumThermalConductivityModel(typeName, U, phi)

    //Cc_(HeliumThermalConductivityModelCoeffs_.subDict(type() + "Coeffs").lookup("Cc")),
    //Cv_(HeliumThermalConductivityModelCoeffs_.subDict(type() + "Coeffs").lookup("Cv")),
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::HeliumThermalConductivityModels::Hepak::calckHe(const HeliumModel& hm)
{
	Info<< "Hepak thermal conductivity model  updates kHe..." << endl;
	setkHe(pow(hm.onebyf()/magGradT2(), 1./3));
}

//bool Foam::HeliumThermalConductivityModels::Hepak::read
//(
//    const dictionary& HeliumProperties
//)
//{
//    HeliumThermalConductivityModel::read(HeliumProperties);
//
//    //HepakCoeffs_ = HeliumProperties.subDict(typeName + "Coeffs");
//
//    //HepakCoeffs_.lookup("TMean") >> TMean_;
//
//    return true;
//}


// ************************************************************************* //
