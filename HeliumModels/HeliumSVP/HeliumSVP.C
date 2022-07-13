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

#include "HeliumSVP.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace HeliumModels
{
    defineTypeNameAndDebug(HeliumSVP, 0);

    addToRunTimeSelectionTable
    (
        HeliumModel,
        HeliumSVP,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::HeliumModels::HeliumSVP::calcNu() 
{
	Info<< "Jestem w calcNu() w HeliumSVP. " << endl;
	calcHeProp(etaHe_, etaHeTable_, T_);
	calcHeProp(rhoHe_, rhoHeTable_, T_);

    volScalarField nu
    (
        IOobject
        (
            "nuLocal",
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
		etaHe_/rhoHe_
    );

    return max
    (
        nuMin_,
        min
        (
            nuMax_,
			nu
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HeliumModels::HeliumSVP::HeliumSVP
(
    const word& name,
    const dictionary& HeliumSVPProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    HeliumModel(name, HeliumSVPProperties, U, phi),
    HeliumSVPCoeffs_(HeliumSVPProperties.subDict(typeName + "Coeffs")),
	T_(U.db().lookupObject<volScalarField>("T")),
    nuMin_("nuMin", dimViscosity, etaHeTable_[indexMin_]/rhoHeTable_[indexMin_]),
    nuMax_("nuMax", dimViscosity, etaHeTable_[indexMax_]/rhoHeTable_[indexMax_])
{

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::HeliumModels::HeliumSVP::correct()
{
	Info<< "HeliumSVP updates thermal properties..." << endl;
	nu_ = calcNu();
	// TODO: dodano update rhon i rhos w pozostalych modelach HeII
	// ale nie sprawdzono czy to jest ok zrobione
	rhon_ = rhoHe_*pow(max(T_/Tlambda_, dimensionedScalar("small", dimless, SMALL)), scalar(5.6));
	rhos_ = rhoHe_ - rhon_;
	calcHeProp(betaHe_, betaHeTable_, T_);
	calcHeProp(AGMHe_, AGMHeTable_, T_);
	calcHeProp(sHe_, sHeTable_, T_);
	calcHeProp(cpHe_, cpHeTable_, T_);
	calcHeProp(onebyf_, onebyfSVPTable_, T_);
}

bool Foam::HeliumModels::HeliumSVP::read
(
    const dictionary& HeliumSVPProperties
)
{
    HeliumModel::read(HeliumSVPProperties);

    //HeliumSVPCoeffs_ = HeliumSVPProperties.subDict(typeName + "Coeffs");

    //HeliumSVPCoeffs_.lookup("rhoHe") >> rhoHe_;

    return true;
}


// ************************************************************************* //
