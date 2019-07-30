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

#include "tempHerschelBukley.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(tempHerschelBukley, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        tempHerschelBukley,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::tempHerschelBukley::calcNu() const
{

	dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);

    tmp<volScalarField> sr(strainRate());
	
//inicio da edição
	const volScalarField& T=U_.mesh().lookupObject<volScalarField>("T");

	return
	(
        (tauTIAC_*(
					exp(
						STau_*(
								dimensionedScalar("one", dimless, 1.0)/(
																		max(
																			T, dimensionedScalar ("VSMALL", dimTemperature, VSMALL)
																			)
																		)
								-dimensionedScalar("one", dimless, 1.0)/TIAC_
								)
						)
					-dimensionedScalar("one", dimless, 1.0)
					) 
		+ kTIAC_*(
					exp(
						kTau_*(
								dimensionedScalar("one", dimless, 1.0)/(
																		max(
																			T, dimensionedScalar ("VSMALL", dimTemperature, VSMALL)
																			)
																		)
								-dimensionedScalar("one", dimless, 1.0)/TIAC_
								)
						)
					-dimensionedScalar("one", dimless, 1.0)
					)
				*rtone*pow(tone*sr(), n_)
		+nuRef_*(
					exp(
						SNu_*(
								dimensionedScalar("one", dimless, 1.0)/(
																		max(
																			T, dimensionedScalar ("VSMALL", dimTemperature, VSMALL)
																			)
																		)
								-dimensionedScalar("one", dimless, 1.0)/TRef_
								)
						)
					)
				*sr()
		)
        /(
			max(
				sr(), dimensionedScalar ("VSMALL", dimless/dimTime, VSMALL)
				)
		 )
    );
//final da edição
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::tempHerschelBukley::tempHerschelBukley
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    tempHerschelBukleyCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
//inicio da edição    
	n_("n", dimless, tempHerschelBukleyCoeffs_),
	tauTIAC_("tauTIAC", dimViscosity/dimTime, tempHerschelBukleyCoeffs_),
	STau_("STau", dimTemperature, tempHerschelBukleyCoeffs_),
	TIAC_("TIAC", dimTemperature, tempHerschelBukleyCoeffs_),
	kTIAC_("kTIAC", dimViscosity, tempHerschelBukleyCoeffs_),
	kTau_("kTau", dimTemperature, tempHerschelBukleyCoeffs_),
	nuRef_("nuRef", dimViscosity, tempHerschelBukleyCoeffs_),	
	SNu_("SNu", dimTemperature, tempHerschelBukleyCoeffs_),
	TRef_("TRef", dimTemperature, tempHerschelBukleyCoeffs_),
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

bool Foam::viscosityModels::tempHerschelBukley::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    tempHerschelBukleyCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

//inicio da edição
	tempHerschelBukleyCoeffs_.lookup("n") >> n_;
	tempHerschelBukleyCoeffs_.lookup("tauTIAC") >> tauTIAC_;
	tempHerschelBukleyCoeffs_.lookup("STau") >> STau_;
	tempHerschelBukleyCoeffs_.lookup("TIAC") >> TIAC_;
	tempHerschelBukleyCoeffs_.lookup("kTIAC") >> kTIAC_;
	tempHerschelBukleyCoeffs_.lookup("kTau") >> kTau_;
	tempHerschelBukleyCoeffs_.lookup("nuRef") >> nuRef_;
	tempHerschelBukleyCoeffs_.lookup("SNu") >> SNu_;
	tempHerschelBukleyCoeffs_.lookup("TRef") >> TRef_;
//final da edição

    return true;
}


// ************************************************************************* //
