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

#include "HeliumModel.H"
//#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(HeliumModel, 0);
    defineRunTimeSelectionTable(HeliumModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HeliumModel::HeliumModel
(
    const word& name,
    const dictionary& HeliumProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    name_(name),
    HeliumProperties_(HeliumProperties),
    U_(U),
    phi_(phi),
	hl_(U, phi),
	p_(HeliumProperties_.lookup("heliumPressure")),
//    TMinField_
//    (
//        IOobject
//        (
//            "TMin",
//            U.mesh().time().timeName(),
//            U.mesh(),
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//        U.mesh(),
//		TMin_
//    ),
//    TMaxField_
//    (
//        IOobject
//        (
//            "TMax",
//            U.mesh().time().timeName(),
//            U.mesh(),
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//        U.mesh(),
//		TMax_
//    ),
    betaHe_
    (
        IOobject
        (
            "betaHe",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
		dimensionedScalar("betaHe", dimless/dimTemperature, 0.0)//,
		//"zeroGradient"
    ),

    AGMHe_
    (
        IOobject
        (
            "AGMHe",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
		dimensionedScalar("AGM", dimensionSet(-1,1,1,0,0,0,0), 0.0)//,
		//"zeroGradient"
    ),

    sHe_
    (
        IOobject
        (
            "sHe",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
		dimensionedScalar("sHe", dimensionSet(0,2,-2,-1,0,0,0), 0.0)//,
		//"zeroGradient"
    ),

    etaHe_
    (
        IOobject
        (
            "etaHe",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
		dimensionedScalar("etaHe", dimensionSet(1,-1,-1,0,0,0,0), 0.0)//,
		//"zeroGradient"
    ),

    cpHe_
    (
        IOobject
        (
            "cpHe",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
		dimensionedScalar("cpHe", dimensionSet(0,2,-2,-1,0,0,0), 0.0)//,
		//"zeroGradient"
    ),
	
    onebyf_
    (
        IOobject
        (
            "onebyf",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
		dimensionedScalar("onebyf", dimensionSet(3,1,-9,-1,0,0,0), 0.0)//,
		//"zeroGradient"
    ),

    rhoHe_
    (
        IOobject
        (
            "rhoHe",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
		dimensionedScalar("rhoHe", dimDensity, 0.0)//,
		//"zeroGradient"
    ),

    rhon_
    (
        IOobject
        (
            "rhon",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
		dimensionedScalar("rhon", dimDensity, 0.0)//,
		//"zeroGradient"
    ),

    rhos_
    (
        IOobject
        (
            "rhos",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
		dimensionedScalar("rhos", dimDensity, 0.0)//,
		//"zeroGradient"
    ),

    nu_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        //calcNu()
        U.mesh(),
		dimensionedScalar("nuHe", dimViscosity, 0.0)//,
		//"zeroGradient"
    )

{

	Info<< "\nHelium operating pressure is " << p_ << endl << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


bool Foam::HeliumModel::read(const dictionary& HeliumProperties)
{
    HeliumProperties_ = HeliumProperties;

    return true;
}


// ************************************************************************* //
