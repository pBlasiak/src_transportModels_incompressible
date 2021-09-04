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
    Foam::HeliumModels::HeliumConstRho

Description
     Helium II model with thermophysical properties dependent on temperature
	 except density (rho is constant).

SourceFiles
    HeliumConstRho.C

\*---------------------------------------------------------------------------*/

#ifndef HeliumConstRho_H
#define HeliumConstRho_H

#include "HeliumModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "HeliumConst.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace HeliumModels
{

/*---------------------------------------------------------------------------*\
                           Class HeliumConstRho Declaration
\*---------------------------------------------------------------------------*/

class HeliumConstRho
:
    public HeliumConst
{
    // Protected data
protected:

//        dictionary HeliumConstRhoCoeffs_;

        const volScalarField& T_;

        dimensionedScalar nuMin_;
        dimensionedScalar nuMax_;


    // Private Member Functions

        //- Calculate and return the laminar Helium
        tmp<volScalarField> calcNu() ;


public:

    //- Runtime type information
    TypeName("HeliumConstRho");


    // Constructors

        //- Construct from components
        HeliumConstRho
        (
            const word& name,
            const dictionary& HeliumProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~HeliumConstRho()
    {}


    // Member Functions

        //- Correct the laminar Helium and thermal properties of He
        virtual void correct() override;

        //- Read transportProperties dictionary
//        virtual bool read(const dictionary& HeliumProperties) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace HeliumModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
