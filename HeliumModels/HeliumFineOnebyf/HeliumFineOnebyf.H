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
    Foam::HeliumModels::HeliumFineOnebyf

Description
     Helium II model with thermophysical properties dependent on temperature.
	 Onebyf function is much more accurate than in the Helium model (dT = 0.0001 K).

SourceFiles
    HeliumFineOnebyf.C

\*---------------------------------------------------------------------------*/

#ifndef HeliumFineOnebyf_H
#define HeliumFineOnebyf_H

#include "HeliumModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"
#include "Helium.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace HeliumModels
{

/*---------------------------------------------------------------------------*\
                           Class HeliumFineOnebyf Declaration
\*---------------------------------------------------------------------------*/

class HeliumFineOnebyf
:
    public Helium
{
    // Private data

		static const label largeIndex_;

		//- Step of temperature in the thermophysical tables dT = 0.0001 K
		static const scalar dTFineOnebyf_;

		//- Table for betaHe
		static const List<scalar> fineOnebyfTable_;


    // Private Member Functions


public:

    //- Runtime type information
    TypeName("HeliumFineOnebyf");


    // Constructors

        //- Construct from components
        HeliumFineOnebyf
        (
            const word& name,
            const dictionary& HeliumProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~HeliumFineOnebyf()
    {}


    // Member Functions

        //- Correct the laminar Helium and thermal properties of He
        virtual void correct() override;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace HeliumModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
