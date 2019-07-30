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

#include "ArrheniusPowerLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(ArrheniusPowerLaw, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        ArrheniusPowerLaw,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::ArrheniusPowerLaw::calcNu() const
{

	dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);

    tmp<volScalarField> sr(strainRate());
	
//inicio da edição
	const volScalarField& T=U_.mesh().lookupObject<volScalarField>("T");

	return max
    (
        nuMin_,
        min
        (
            nuMax_,
            kRef_*(
					exp(
						Sk_*(
								dimensionedScalar("one", dimless, 1.0)/ max(
																			T, dimensionedScalar ("VSMALL", dimTemperature, VSMALL)
																			)
																		
								-dimensionedScalar("one", dimless, 1.0)/TRef_
								)
						)
				)
			*pow
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

//final da edição
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::ArrheniusPowerLaw::ArrheniusPowerLaw
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    ArrheniusPowerLawCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
//inicio da edição    
	kRef_("kRef", dimViscosity, ArrheniusPowerLawCoeffs_),	
	Sk_("Sk", dimTemperature, ArrheniusPowerLawCoeffs_),
	TRef_("TRef", dimTemperature, ArrheniusPowerLawCoeffs_),
	n_("n", dimless, ArrheniusPowerLawCoeffs_),
    nuMin_("nuMin", dimViscosity, ArrheniusPowerLawCoeffs_),
    nuMax_("nuMax", dimViscosity, ArrheniusPowerLawCoeffs_),
//final da edição
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

bool Foam::viscosityModels::ArrheniusPowerLaw::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    ArrheniusPowerLawCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

//inicio da edição
	ArrheniusPowerLawCoeffs_.lookup("kRef") >> kRef_;
	ArrheniusPowerLawCoeffs_.lookup("Sk") >> Sk_;
	ArrheniusPowerLawCoeffs_.lookup("TRef") >> TRef_;
	ArrheniusPowerLawCoeffs_.lookup("n") >> n_;
	ArrheniusPowerLawCoeffs_.lookup("nuMin") >> nuMin_;
	ArrheniusPowerLawCoeffs_.lookup("nuMax") >> nuMax_;
//final da edição

    return true;
}


// ************************************************************************* //
