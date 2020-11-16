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

#include "turnOnPlasticity2.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(turnOnPlasticity2, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        turnOnPlasticity2,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::turnOnPlasticity2::calcNu() const
{

	dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);
	dimensionedScalar TVSMALL("VSMALL", dimTemperature, VSMALL);

    tmp<volScalarField> sr(strainRate());
	
	const volScalarField& T=U_.mesh().lookupObject<volScalarField>("T");
	
	return
	(
		min
        	(
			nuMax_,
			(tauYRef_*(
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
			+ kRef_*(
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
			+nu0Ref_*(
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
}


Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::turnOnPlasticity2::calcTauY() const
{
	dimensionedScalar TVSMALL("VSMALL", dimTemperature, VSMALL);
	const volScalarField& T=U_.mesh().lookupObject<volScalarField>("T");

	return
	(     	
		tauYRef_*(
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
    );
}

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::turnOnPlasticity2::calcK() const
{
	dimensionedScalar TVSMALL("VSMALL", dimTemperature, VSMALL);
	const volScalarField& T=U_.mesh().lookupObject<volScalarField>("T");

	return
	(
		kRef_*(
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
    );
}
Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::turnOnPlasticity2::calcNu0() const
{
	dimensionedScalar TVSMALL("VSMALL", dimTemperature, VSMALL);
	const volScalarField& T=U_.mesh().lookupObject<volScalarField>("T");

	return
	(
		nu0Ref_*(
					exp(
						SNu_*(
								scalar(1.0)/max(T, TVSMALL)-scalar(1.0)/TRef_
							)
						)
				)
			
	);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityModels::turnOnPlasticity2::turnOnPlasticity2
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    turnOnPlasticity2Coeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
//inicio da edição    
	n_("n", dimless, turnOnPlasticity2Coeffs_),
	tauYRef_("tauYRef", dimViscosity/dimTime, turnOnPlasticity2Coeffs_),
	STau_("STau", dimTemperature, turnOnPlasticity2Coeffs_),
	TYS_("TYS", dimTemperature, turnOnPlasticity2Coeffs_),
	TPL_("TPL", dimTemperature, turnOnPlasticity2Coeffs_),
	kRef_("kRef", dimViscosity, turnOnPlasticity2Coeffs_),
	Sk_("Sk", dimTemperature, turnOnPlasticity2Coeffs_),
	nu0Ref_("nu0Ref", dimViscosity, turnOnPlasticity2Coeffs_),	
	SNu_("SNu", dimTemperature, turnOnPlasticity2Coeffs_),
	TRef_("TRef", dimTemperature, turnOnPlasticity2Coeffs_),
	nuMax_("nuMax", dimViscosity, turnOnPlasticity2Coeffs_),
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
    ),
    tauY_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcTauY()
    ),
    k_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcK()
    ),
    nu0_
    (
        IOobject
        (
            name,
            U_.time().timeName(),
            U_.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        calcNu0()
    )   
   
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::viscosityModels::turnOnPlasticity2::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    turnOnPlasticity2Coeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

//inicio da edição
	turnOnPlasticity2Coeffs_.lookup("n") >> n_;
	turnOnPlasticity2Coeffs_.lookup("tauYRef") >> tauYRef_;
	turnOnPlasticity2Coeffs_.lookup("STau") >> STau_;
	turnOnPlasticity2Coeffs_.lookup("TYS") >> TYS_;
	turnOnPlasticity2Coeffs_.lookup("TPL") >> TPL_;
	turnOnPlasticity2Coeffs_.lookup("kRef") >> kRef_;
	turnOnPlasticity2Coeffs_.lookup("Sk") >> Sk_;
	turnOnPlasticity2Coeffs_.lookup("nu0Ref") >> nu0Ref_;
	turnOnPlasticity2Coeffs_.lookup("SNu") >> SNu_;
	turnOnPlasticity2Coeffs_.lookup("TRef") >> TRef_;
	turnOnPlasticity2Coeffs_.lookup("nuMax") >> nuMax_;
//final da edição

    return true;
}


// ************************************************************************* //
