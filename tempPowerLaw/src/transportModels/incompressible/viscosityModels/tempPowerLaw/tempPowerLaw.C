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

#include "tempPowerLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(tempPowerLaw, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        tempPowerLaw,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::tempPowerLaw::calcNu() const
{
	
//inicio da edição
	const volScalarField& T=U_.mesh().lookupObject<volScalarField>("T");
//final da edição

    return max
    (
        nuMin_,
        min
        (
            nuMax_,
//inicio da edição
			(k_-kSlope_*(T-TBase_))*pow
//final da edição
            (
                max
                (
                    dimensionedScalar("one", dimTime, 1.0)*strainRate(),
                    dimensionedScalar("VSMALL", dimless, VSMALL)
                ),
                n_.value() - scalar(1.0)
            )
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::tempPowerLaw::tempPowerLaw
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    tempPowerLawCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
    k_("k", dimViscosity, tempPowerLawCoeffs_),
    n_("n", dimless, tempPowerLawCoeffs_),
//inicio da edição
	kSlope_(tempPowerLawCoeffs_.lookup("kSlope")),
	TBase_(tempPowerLawCoeffs_.lookup("TBase")),
//final da edição
    nuMin_("nuMin", dimViscosity, tempPowerLawCoeffs_),
    nuMax_("nuMax", dimViscosity, tempPowerLawCoeffs_),
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
        calcNu()
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::tempPowerLaw::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    tempPowerLawCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

    tempPowerLawCoeffs_.lookup("k") >> k_;
    tempPowerLawCoeffs_.lookup("n") >> n_;
//inicio da edição
	tempPowerLawCoeffs_.lookup("kSlope") >> kSlope_;
	tempPowerLawCoeffs_.lookup("TBase") >> TBase_;
//final da edição
    tempPowerLawCoeffs_.lookup("nuMin") >> nuMin_;
    tempPowerLawCoeffs_.lookup("nuMax") >> nuMax_;

    return true;
}


// ************************************************************************* //
