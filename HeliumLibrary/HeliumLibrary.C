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

#include "HeliumLibrary.H"
//#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	const Enum<HeliumLibrary::HeliumThermalPropertyType>
	HeliumLibrary::HeliumThermalPropertyTypeNames_
	({
        { HeliumThermalPropertyType::thermalExpansion, "beta" },
        { HeliumThermalPropertyType::AGMCoeff, "AGM" },
        { HeliumThermalPropertyType::entropy, "s" },
        { HeliumThermalPropertyType::dynamicViscosity, "eta" },
        { HeliumThermalPropertyType::specificHeatCapacity, "cp" },
        { HeliumThermalPropertyType::oneByf, "oneByf" },
        { HeliumThermalPropertyType::density, "rho" },
        { HeliumThermalPropertyType::thermPropsSize_, "thermPropsSize" }
	 }
	);

	const Enum<HeliumLibrary::HeliumPressure>
	HeliumLibrary::HeliumPressureNames_
	({
        { HeliumPressure::SVP,    "SVP"  },
        { HeliumPressure::onebar, "1bar" },
        { HeliumPressure::pressuresSize_,  "pressuresSize" }
	 }
	);

	const Foam::label
	Foam::HeliumLibrary::indexMin_(0);
	
	const Foam::label
	Foam::HeliumLibrary::indexMax_(667);

	const Foam::dimensionedScalar
	Foam::HeliumLibrary::dT_("dT", dimTemperature, 0.001);

	const Foam::dimensionedScalar
	Foam::HeliumLibrary::Tlambda_("Tlambda", dimTemperature, 2.16795);

	const Foam::dimensionedScalar
	Foam::HeliumLibrary::TMin_("TMin", dimTemperature, 1.5);
	
	const Foam::dimensionedScalar
	Foam::HeliumLibrary::TMax_("TMax", dimTemperature, 2.167);
	

	#include "staticTables.H"
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::HeliumLibrary::HeliumLibrary
(
    //const word& name,
    //const dictionary& HeliumProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    //name_(name),
    //HeliumProperties_(HeliumProperties),
	//heliumPressure_(HeliumPressures::onebar),
    U_(U),
    phi_(phi),
    //TMinField_
    //(
    //    IOobject
    //    (
    //        "TMin",
    //        U.mesh().time().timeName(),
    //        U.mesh(),
    //        IOobject::NO_READ,
    //        IOobject::NO_WRITE
    //    ),
    //    U.mesh(),
	//	TMin_
    //),
    //TMaxField_
    //(
    //    IOobject
    //    (
    //        "TMax",
    //        U.mesh().time().timeName(),
    //        U.mesh(),
    //        IOobject::NO_READ,
    //        IOobject::NO_WRITE
    //    ),
    //    U.mesh(),
	//	TMax_
    //),
	betaHeTables_(pressuresSize_),
	AGMHeTables_(pressuresSize_),
	sHeTables_(pressuresSize_),
	etaHeTables_(pressuresSize_),
	cpHeTables_(pressuresSize_),
	onebyfTables_(pressuresSize_),
	rhoHeTables_(pressuresSize_),

	HeThermPropsTables_(thermPropsSize_)
	//HeThermProps_(7),
	//HeThermPropsTables_(7)
{ 
	//TODO: sprawdz czy dobrze pointers sa przypisane
	betaHeTables_.set(SVP,    &betaHeTableSVP_);
	betaHeTables_.set(onebar, &betaHeTable1bar_);
	AGMHeTables_.set(SVP,     &AGMHeTableSVP_);
	AGMHeTables_.set(onebar,  &AGMHeTable1bar_);
	sHeTables_.set(SVP,       &sHeTableSVP_);
	sHeTables_.set(onebar,    &sHeTable1bar_);
	etaHeTables_.set(SVP,     &etaHeTableSVP_);
	etaHeTables_.set(onebar,  &etaHeTable1bar_);
	cpHeTables_.set(SVP,      &cpHeTableSVP_);
	cpHeTables_.set(onebar,   &cpHeTable1bar_);
	onebyfTables_.set(SVP,    &onebyfTableSVP_);
	onebyfTables_.set(onebar, &onebyfTable1bar_);
	rhoHeTables_.set(SVP,     &rhoHeTableSVP_);
	rhoHeTables_.set(onebar,  &rhoHeTable1bar_); 

	HeThermPropsTables_.set(thermalExpansion, &betaHeTables_);
	HeThermPropsTables_.set(AGMCoeff, &AGMHeTables_);
	HeThermPropsTables_.set(entropy, &sHeTables_);
	HeThermPropsTables_.set(dynamicViscosity, &etaHeTables_);
	HeThermPropsTables_.set(specificHeatCapacity, &cpHeTables_);
	HeThermPropsTables_.set(oneByf, &onebyfTables_);
	HeThermPropsTables_.set(density, &rhoHeTables_);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::HeliumLibrary::calcHeProp
(
    volScalarField& vsf, 
	//const List<scalar>& vsfTable,
	//TODO: sprawdz potem czy sprawdza jak dasz inna wartosc np. 10bar
	const volScalarField& T,
	const word wtp,
	const word wp,
	//HeliumThermalPropertiesType tp,
	//HeliumPressures p,
	const label maxIndex, 
	const dimensionedScalar dt
)
{
	const scalar TMin(TMin_.value());
	const scalar TMax(TMax_.value());
	const scalar dT(dt.value());

	// Solution with iterators
	//forAll(T, celli)
	//{
	//	if (T[celli] < TMin)
	//	{
	//		PtrList<const List<scalar> >::const_iterator iterTable = HeThermPropsTables_.begin();
	//		forAllIters(HeThermProps_, iter)
	//		{
	//			iter()[celli] = iterTable()[indexMin_];
	//			iterTable++;
	//		}
	//	}
	//	else if (T[celli] > TMax)
	//	{
	//		PtrList<const List<scalar> >::const_iterator iterTable = HeThermPropsTables_.begin();
	//		forAllIters(HeThermProps_, iter)
	//		{
	//			iter()[celli] = iterTable()[maxIndex];
	//			iterTable++;
	//		}
	//	}
	//	else
	//	{
	//	    label index = (T[celli] - TMin)/dT;
	//	    if (index == maxIndex)
	//	    {
	//	    	PtrList<const List<scalar> >::const_iterator iterTable = HeThermPropsTables_.begin();
	//	    	forAllIters(HeThermProps_, iter)
	//	    	{
	//	    		iter()[celli] = iterTable()[maxIndex];
	//				iterTable++;
	//	    	}
	//	    }
	//	    else
	//	    {
	//	    	scalar Ti1 = TMin + index*dT;
	//	    	scalar Ti2 = Ti1 + dT;
	//	    	PtrList<const List<scalar> >::const_iterator iterTable = HeThermPropsTables_.begin();
	//	    	forAllIters(HeThermProps_, iter)
	//	    	{
	//	    		scalar a = (iterTable()[index + 1] - iterTable()[index])/(Ti2 - Ti1);
	//	    		scalar b = iterTable()[index] - a*Ti1;
	//	    		iter()[celli] = a*T[celli] + b;
	//				iterTable++;
	//	    	}
	//	    }
    //    }
	//}

	//forAllIters(HeThermProps_, iter)
	//{
	//	iter->correctBoundaryConditions();
	//}

	// Old solution with forAll loops 
	
	const List<scalar>& vsfTable = getThermProp(wtp, wp); 
	forAll(vsf, celli)
	{
		if (T[celli] < TMin)
		{
			vsf[celli] = vsfTable[indexMin_];
		}
		else if (T[celli] > TMax)
		{
			vsf[celli] = vsfTable[maxIndex];
		}
		else
		{
			label index = (T[celli] - TMin)/dT;
			if (index == maxIndex)
			{
				vsf[celli] = vsfTable[maxIndex];
			}
			else
			{
				scalar Ti1 = TMin + index*dT;
				scalar Ti2 = Ti1 + dT;
				scalar a = (vsfTable[index + 1] - vsfTable[index])/(Ti2 - Ti1);
				scalar b = vsfTable[index] - a*Ti1;
				vsf[celli] = a*T[celli] + b;
			}
		}
	}

	// if we have this line probably we do not need 
	// the second forAll loop over boundaries but it has to be checked
	//vsf.correctBoundaryConditions();

	forAll(vsf.boundaryField(), patchi)
	{
		forAll(vsf.boundaryField()[patchi], facei)
		{
			if (T[facei] < TMin)
			{
				vsf.boundaryFieldRef()[patchi][facei] = vsfTable[indexMin_];
			}
			else if (T[facei] > TMax)
			{
				vsf.boundaryFieldRef()[patchi][facei] = vsfTable[maxIndex];
			}
			else
			{
				label index = (T[facei] - TMin)/dT;
				if (index == maxIndex)
				{
					vsf.boundaryFieldRef()[patchi][facei] = vsfTable[maxIndex];
				}
				else
				{
					scalar Ti1 = TMin + index*dT;
					scalar Ti2 = Ti1 + dT;
					scalar a = (vsfTable[index + 1] - vsfTable[index])/(Ti2 - Ti1);
					scalar b = vsfTable[index] - a*Ti1;
					vsf.boundaryFieldRef()[patchi][facei] = a*T[facei] + b;
				}
			}
		}
	}
}

const Foam::List<Foam::scalar>&
Foam::HeliumLibrary::getThermProp(const word wtp, const word wp) const
{
	const label tp(HeliumThermalPropertyTypeNames_.get(wtp));
	const label p(HeliumPressureNames_.get(wp));
	return HeThermPropsTables_[tp][p];
}
//
//
//bool Foam::HeliumLibrary::read(const dictionary& HeliumProperties)
//{
//    HeliumProperties_ = HeliumProperties;
//
//    return true;
//}


// ************************************************************************* //
