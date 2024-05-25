/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::HeliumThermalConductivityModelModels::Hepak

Description
    Hepak thermal property model for a mixture. 
	thermProp = alpha1*thermProp1 + alpha2*thermProp2

SourceFiles
    Hepak.C

\*--------------------------------------------------------------------*/

#ifndef Hepak_H
#define Hepak_H

#include "HeliumThermalConductivityModel.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace HeliumThermalConductivityModels
{

/*--------------------------------------------------------------------*\
                              Class Hepak
\*--------------------------------------------------------------------*/

class Hepak
:
    public HeliumThermalConductivityModel
{
    // Flag to change exponent in Gorter-Mellink relation
	// that defines gradT and q relationship
	bool isExponentDifferentFrom3_{false};

	// Exponent in Gorter-Mellink relation
	scalar exponentGM_{3.0};
		

public:

    //- Runtime type information
    TypeName("Hepak");


    // Constructors

        //- construct from components
        Hepak
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~Hepak()
    {}


    // Member Functions
	
        //- Correct the Helium thermal conductivity 
        virtual void calckHe(const HeliumModel& hm) override;

        //- Read transportProperties dictionary
        //virtual bool read(const dictionary& HeliumProperties) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace HeliumThermalConductivityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
