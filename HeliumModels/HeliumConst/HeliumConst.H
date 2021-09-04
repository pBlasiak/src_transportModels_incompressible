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
    Foam::HeliumModels::HeliumConst

Description
     Helium II model with constant thermophysical properties.
	 The properties are calculated at a mean temperature.

SourceFiles
    HeliumConst.C

\*---------------------------------------------------------------------------*/

#ifndef HeliumConst_H
#define HeliumConst_H

#include "HeliumModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace HeliumModels
{

/*---------------------------------------------------------------------------*\
                           Class HeliumConst Declaration
\*---------------------------------------------------------------------------*/


class Helium;

class HeliumConst
:
    public HeliumModel
{
    // Protected data
protected:

        dictionary HeliumConstCoeffs_;

		dimensionedScalar TMean0_;

		volScalarField TMean_;


    // Private Member Functions

        //- Calculate and return the laminar Helium
        tmp<volScalarField> calcNu() ;

		////- Calculates helium thermal property based on mean temperature
        //void calcHeProp(volScalarField& vsf, const List<scalar>& vsfTable);


		//friend class Helium;
		//friend class HeliumConstRho;
public:

    //- Runtime type information
    TypeName("HeliumConst");


    // Constructors

        //- Construct from components
        HeliumConst
        (
            const word& name,
            const dictionary& HeliumProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~HeliumConst()
    {}


    // Member Functions


        //- Correct the laminar Helium and thermal properties of He
        virtual void correct() override;

        //- Read transportProperties dictionary
        virtual bool read(const dictionary& HeliumProperties) override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace HeliumModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
