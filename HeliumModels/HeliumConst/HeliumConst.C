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

#include "HeliumConst.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace HeliumModels
{
    defineTypeNameAndDebug(HeliumConst, 0);

    addToRunTimeSelectionTable
    (
        HeliumModel,
        HeliumConst,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::HeliumModels::HeliumConst::calcNu() 
{
	calcHeProp(etaHe_, etaHeTable_, TMean_);
	calcHeProp(rhoHe_, rhoHeTable_, TMean_);

	return tmp<volScalarField>
	(
	    new volScalarField
		(
		    "nu",
			etaHe_/rhoHe_
		)
	);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HeliumModels::HeliumConst::HeliumConst
(
    const word& name,
    const dictionary& HeliumProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    HeliumModel(name, HeliumProperties, U, phi),
    HeliumConstCoeffs_(HeliumProperties.subDict(typeName + "Coeffs")),
    TMean0_("TMean", dimTemperature, HeliumConstCoeffs_),
	TMean_
    (
        IOobject
        (
            "TMean",
            U.mesh().time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
		TMean0_
    )
{
	Info<< "HeliumConst calculates thermal properties for TMean = " << TMean_ << endl;
	nu_ = calcNu();
	calcHeProp(betaHe_, betaHeTable_, TMean_);
	calcHeProp(AGMHe_, AGMHeTable_, TMean_);
	calcHeProp(sHe_, sHeTable_, TMean_);
	calcHeProp(cpHe_, cpHeTable_, TMean_);
	calcHeProp(onebyf_, onebyfTable_, TMean_);

	Info<< "TMean = " << TMean_ << endl;
	Info<< "betaHe = " << betaHe_ << endl;
	Info<< "AGMHe = " << AGMHe_ << endl;
	Info<< "sHe = " << sHe_ << endl;
	Info<< "cpHe = " << cpHe_ << endl;
	Info<< "onebyf = " << onebyf_ << endl;
	Info<< "rhoHe = " << rhoHe_ << endl;
	Info<< "nuHe = " << nu_ << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::HeliumModels::HeliumConst::correct()
{
	const volScalarField& T = U_.db().lookupObject<volScalarField>("T");
	Info<< "HeliumComst updates rhon and rhos..." << endl;
	rhon_ = rhoHe_*pow(max(T/Tlambda_, dimensionedScalar("small", dimless, SMALL)), scalar(5.6));
	rhos_ = rhoHe_ - rhon_; 
}

bool Foam::HeliumModels::HeliumConst::read
(
    const dictionary& HeliumProperties
)
{
    HeliumModel::read(HeliumProperties);

    HeliumConstCoeffs_ = HeliumProperties.subDict(typeName + "Coeffs");

    HeliumConstCoeffs_.lookup("TMean") >> TMean_;

	Info<< "HeliumConst updates thermal properties..." << endl;
	nu_ = calcNu();
	calcHeProp(betaHe_, betaHeTable_, TMean_);
	calcHeProp(AGMHe_, AGMHeTable_, TMean_);
	calcHeProp(sHe_, sHeTable_, TMean_);
	calcHeProp(cpHe_, cpHeTable_, TMean_);
	calcHeProp(onebyf_, onebyfTable_, TMean_);

	Info<< "TMean = " << TMean_ << endl;
	Info<< "betaHe = " << betaHe_ << endl;
	Info<< "AGMHe = " << AGMHe_ << endl;
	Info<< "sHe = " << sHe_ << endl;
	Info<< "cpHe = " << cpHe_ << endl;
	Info<< "onebyf = " << onebyf_ << endl;
	Info<< "rhoHe = " << rhoHe_ << endl;
	Info<< "nuHe = " << nu_ << endl;

    return true;
}


// ************************************************************************* //
