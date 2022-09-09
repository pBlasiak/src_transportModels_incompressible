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

Namespace
    Foam::

Description
    A namespace for various incompressible HeliumLibrary implementations.

Class
    Foam::HeliumLibrary

Description
	He II thermo-physical properties library.

SourceFiles
    HeliumLibrary.C

\*---------------------------------------------------------------------------*/

#ifndef HeliumLibrary_H
#define HeliumLibrary_H

#include "dictionary.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "dimensionedScalar.H"
#include "runTimeSelectionTables.H"
#include "typeInfo.H"
#include "volFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class HeliumLibrary Declaration
\*---------------------------------------------------------------------------*/

class HeliumLibrary
{
    // Public data
public:

        //- Enumeration listing the possible helium thermal properties 
        enum HeliumThermalPropertyType
        {
            thermalExpansion = 0,
			AGMCoeff,
			entropy,
			dynamicViscosity,
			specificHeatCapacity,
			oneByf,
			//TODO: dodac fineonebyf???
			density,
			thermPropsSize_
        };
		
        //- Helium thermal properties types
        static const Enum<HeliumThermalPropertyType> HeliumThermalPropertyTypeNames_;

        //- Enumeration listing the possible helium pressures
        enum HeliumPressure
        {
			SVP = 0,
			onebar,
			pressuresSize_
        };

        //- Helium pressures names
        static const Enum<HeliumPressure> HeliumPressureNames_;

private:

    // Private data
	
	//TODO zamien const na constexpr????
	
		//- Minimum index number for thermophysical tables i = 0
		static const label indexMin_;

		//- Maximum index number for thermophysical tables i = 667
		static const label indexMax_;

		//- Step of temperature in the thermophysical tables dT = 0.001 K
		static const dimensionedScalar dT_;
	
		//- lambda Temperature 2.16795 K from Hepak at 1 atm
		//  previous value was 2.1711132461 (probably taken from Wiki) 
		//  Van Sciver str. 163 Tlambda = 2.177 K
		static const dimensionedScalar Tlambda_;

		//- Minimum temperature TMin = 1.5 K
		static const dimensionedScalar TMin_;

		//- Maximum temperature TMax = 2.167 K
		static const dimensionedScalar TMax_;

		//- Static tables defined in staticTables.H file
		//- Table for thermal expansion coefficient, 1/K
		static const List<scalar> betaHeTableSVP_;
		static const List<scalar> betaHeTable1bar_;

		//- Table for Gorter-Mellink coefficient, (m s)/kg
		static const List<scalar> AGMHeTableSVP_;
		static const List<scalar> AGMHeTable1bar_;

		//- Table for entropy, J/(kg K)
		static const List<scalar> sHeTableSVP_;
		static const List<scalar> sHeTable1bar_;

		//- Table for dynamic viscosity, Pa s
		static const List<scalar> etaHeTableSVP_;
		static const List<scalar> etaHeTable1bar_;

		//- Table for specific heat capacity at constant pressure, J/(kg K)
		static const List<scalar> cpHeTableSVP_;
		static const List<scalar> cpHeTable1bar_;

		//- Table for thermalc conductivity function f^(-1), W^3/(m^5 K)
		static const List<scalar> onebyfTableSVP_;
		static const List<scalar> onebyfTable1bar_;

		//- Table for density, kg/m^3
		static const List<scalar> rhoHeTableSVP_;
		static const List<scalar> rhoHeTable1bar_;


        //word name_;
        //dictionary HeliumProperties_;

        const volVectorField& U_;
        const surfaceScalarField& phi_;

		//TODO: add fineOneByf ??
		//- Pointers to table data of helium
		PtrList<const List<scalar>> betaHeTables_;
		PtrList<const List<scalar>> AGMHeTables_;
		PtrList<const List<scalar>> sHeTables_;
		PtrList<const List<scalar>> etaHeTables_;
		PtrList<const List<scalar>> cpHeTables_;
		PtrList<const List<scalar>> onebyfTables_;
		PtrList<const List<scalar>> rhoHeTables_;

		PtrList<PtrList<const List<scalar>>> HeThermPropsTables_;


private:

    // Private Member Functions

        //- Disallow copy construct
        HeliumLibrary(const HeliumLibrary&);

        //- Disallow default bitwise assignment
        void operator=(const HeliumLibrary&);


public:


    // Constructors

        //- Construct from components
        HeliumLibrary
        (
            //const word& name,
            //const dictionary& HeliumProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~HeliumLibrary()
    {}


    // Member Functions

		static const label& indexMin() 
		{
			return indexMin_;
		}

		static const label& indexMax() 
		{
			return indexMax_;
		}

		static const dimensionedScalar& Tlambda() 
		{
			return Tlambda_;
		}

		//- Calculates helium thermal property based on temperature
        void calcHeProp
		(
		    volScalarField& vsf, 
			//const List<scalar>& vsfTable,
			const volScalarField& T,
			const word wthermalProperty,
			const word wpressure = "SVP",
			//HeliumThermalPropertyType tp,
			//HeliumPressure p = SVP,
			const label maxIndex = indexMax_,
			const dimensionedScalar dt = dT_
		);
		
		const List<scalar>& getThermProp
		(
			const word wthermalProperty, 
			const word wpressure
		) const;
        ////- Correct the laminar Helium
        //virtual void correct() = 0;

        ////- Read transportProperties dictionary
        //virtual bool read(const dictionary& HeliumProperties) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
