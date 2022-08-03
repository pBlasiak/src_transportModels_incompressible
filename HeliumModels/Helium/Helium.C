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

#include "Helium.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace HeliumModels
{
    defineTypeNameAndDebug(Helium, 0);

    addToRunTimeSelectionTable
    (
        HeliumModel,
        Helium,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::HeliumModels::Helium::calcNu() 
{
	Info<< "Jestem w calcNu() w Helium. " << endl;
	hl_.calcHeProp(etaHe_, T_, "eta", p_);
	hl_.calcHeProp(rhoHe_, T_, "rho", p_);

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

Foam::HeliumModels::Helium::Helium
(
    const word& name,
    const dictionary& HeliumProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    HeliumModel(name, HeliumProperties, U, phi),
    HeliumCoeffs_(HeliumProperties.subDict(typeName + "Coeffs")),
	T_(U.db().lookupObject<volScalarField>("T")),
    //nuMin_("nuMin", dimViscosity, etaHeTable_[indexMin_]/rhoHeTable_[indexMin_]),
    //nuMax_("nuMax", dimViscosity, etaHeTable_[indexMax_]/rhoHeTable_[indexMax_])
    nuMin_
	(
		"nuMin", 
		dimViscosity, 
		hl_.getThermProp("eta", p_)[hl_.indexMin()]/
		hl_.getThermProp("rho", p_)[hl_.indexMin()]
	),
    nuMax_
	(
		"nuMax", 
		dimViscosity,
		hl_.getThermProp("eta", p_)[hl_.indexMax()]/
		hl_.getThermProp("rho", p_)[hl_.indexMax()]
	)
{

}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::HeliumModels::Helium::correct()
{
	Info<< "Helium updates thermal properties..." << endl;
	nu_ = calcNu();
	// TODO: dodano update rhon i rhos w pozostalych modelach HeII
	// ale nie sprawdzono czy to jest ok zrobione
	const dimensionedScalar Tlambda{hl_.Tlambda()};
	rhon_ = rhoHe_*pow(max(T_/Tlambda, dimensionedScalar("small", dimless, SMALL)), scalar(5.6));
	rhos_ = rhoHe_ - rhon_;
	hl_.calcHeProp(betaHe_, T_, "beta", p_);
	hl_.calcHeProp(AGMHe_, T_, "AGM", p_);
	hl_.calcHeProp(sHe_, T_, "s", p_);
	hl_.calcHeProp(cpHe_, T_, "cp", p_);
	hl_.calcHeProp(onebyf_, T_, "oneByf", p_);
}

bool Foam::HeliumModels::Helium::read
(
    const dictionary& HeliumProperties
)
{
    HeliumModel::read(HeliumProperties);

    //HeliumCoeffs_ = HeliumProperties.subDict(typeName + "Coeffs");

    //HeliumCoeffs_.lookup("rhoHe") >> rhoHe_;

    return true;
}


// ************************************************************************* //
