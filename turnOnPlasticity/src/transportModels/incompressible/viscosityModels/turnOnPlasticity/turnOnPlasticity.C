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

#include "turnOnPlasticity.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(turnOnPlasticity, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        turnOnPlasticity,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::turnOnPlasticity::calcNu() const
{

	dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);
	dimensionedScalar TVSMALL("VSMALL", dimTemperature, VSMALL);

    tmp<volScalarField> sr(strainRate());
	
//inicio da edição
	const volScalarField& T=U_.mesh().lookupObject<volScalarField>("T");

	return
	(
		min
        	(
			nuMax_,
			(tauY_*(
						max(
							scalar(0.0),
							exp(
								STau_*(
										scalar(1.0)/max(T, TVSMALL)-scalar(1.0)/TYS_
										)
								)
							-scalar(1.0)
							)
						) 
			+ k_*(
						max(
							scalar(0.0),
							exp(
								Sk_*(
										scalar(1.0)/max(T, TVSMALL)-scalar(1.0)/TPL_
										)
								)
							-scalar(1.0)
							)
						)
					*rtone*pow(tone*sr(), n_)
			+nuRef_*(
						exp(
							SNu_*(
									scalar(1.0)/max(T, TVSMALL)-scalar(1.0)/TRef_
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
			)
    );
//final da edição
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::turnOnPlasticity::turnOnPlasticity
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    turnOnPlasticityCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
//inicio da edição    
	n_("n", dimless, turnOnPlasticityCoeffs_),
	tauY_("tauY", dimViscosity/dimTime, turnOnPlasticityCoeffs_),
	STau_("STau", dimTemperature, turnOnPlasticityCoeffs_),
	TYS_("TYS", dimTemperature, turnOnPlasticityCoeffs_),
	TPL_("TPL", dimTemperature, turnOnPlasticityCoeffs_),
	k_("k", dimViscosity, turnOnPlasticityCoeffs_),
	Sk_("Sk", dimTemperature, turnOnPlasticityCoeffs_),
	nuRef_("nuRef", dimViscosity, turnOnPlasticityCoeffs_),	
	SNu_("SNu", dimTemperature, turnOnPlasticityCoeffs_),
	TRef_("TRef", dimTemperature, turnOnPlasticityCoeffs_),
	nuMax_("nuMax", dimViscosity, turnOnPlasticityCoeffs_),
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

bool Foam::viscosityModels::turnOnPlasticity::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    turnOnPlasticityCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

//inicio da edição
	turnOnPlasticityCoeffs_.lookup("n") >> n_;
	turnOnPlasticityCoeffs_.lookup("tauY") >> tauY_;
	turnOnPlasticityCoeffs_.lookup("STau") >> STau_;
	turnOnPlasticityCoeffs_.lookup("TYS") >> TYS_;
	turnOnPlasticityCoeffs_.lookup("TPL") >> TPL_;
	turnOnPlasticityCoeffs_.lookup("k") >> k_;
	turnOnPlasticityCoeffs_.lookup("Sk") >> Sk_;
	turnOnPlasticityCoeffs_.lookup("nuRef") >> nuRef_;
	turnOnPlasticityCoeffs_.lookup("SNu") >> SNu_;
	turnOnPlasticityCoeffs_.lookup("TRef") >> TRef_;
	turnOnPlasticityCoeffs_.lookup("nuMax") >> nuMax_;
//final da edição

    return true;
}


// ************************************************************************* //
