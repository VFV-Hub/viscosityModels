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

#include "ArrheniusHB.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{
    defineTypeNameAndDebug(ArrheniusHB, 0);

    addToRunTimeSelectionTable
    (
        viscosityModel,
        ArrheniusHB,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::viscosityModels::ArrheniusHB::calcNu() const
{

	dimensionedScalar tone("tone", dimTime, 1.0);
    dimensionedScalar rtone("rtone", dimless/dimTime, 1.0);

    tmp<volScalarField> sr(strainRate());
	
//inicio da edição
	const volScalarField& T=U_.mesh().lookupObject<volScalarField>("T");
	return
    (
        min
        (
            nu0_,
            (
			 tauRef_ *(
						exp(
							SNu_*(
									dimensionedScalar("one", dimless, 1.0)/ max(
																				T, dimensionedScalar ("VSMALL", dimTemperature, VSMALL)
																				)
									-dimensionedScalar("one", dimless, 1.0)/TRef_
								)
							)
					)
			
			+ kRef_*(
						exp(
							Sk_*(
									dimensionedScalar("one", dimless, 1.0)/ max(
																				T, dimensionedScalar ("VSMALL", dimTemperature, VSMALL)
																				)
									-dimensionedScalar("one", dimless, 1.0)/TRef_
								)
							)
					)
			*rtone*pow(tone*sr(), n_)
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

Foam::viscosityModels::ArrheniusHB::ArrheniusHB
(
    const word& name,
    const dictionary& viscosityProperties,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    viscosityModel(name, viscosityProperties, U, phi),
    ArrheniusHBCoeffs_(viscosityProperties.subDict(typeName + "Coeffs")),
//inicio da edição    
	nu0_("nu0", dimViscosity, ArrheniusHBCoeffs_),
	tauRef_("tauRef", dimViscosity/dimTime, ArrheniusHBCoeffs_),	
	SNu_("SNu", dimTemperature, ArrheniusHBCoeffs_),
	TRef_("TRef", dimTemperature, ArrheniusHBCoeffs_),
	kRef_("kRef", dimViscosity, ArrheniusHBCoeffs_),	
	Sk_("Sk", dimTemperature, ArrheniusHBCoeffs_),
	n_("n", dimless, ArrheniusHBCoeffs_),
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

bool Foam::viscosityModels::ArrheniusHB::read
(
    const dictionary& viscosityProperties
)
{
    viscosityModel::read(viscosityProperties);

    ArrheniusHBCoeffs_ = viscosityProperties.subDict(typeName + "Coeffs");

//inicio da edição
	ArrheniusHBCoeffs_.lookup("nu0") >> nu0_;
	ArrheniusHBCoeffs_.lookup("tauRef") >> tauRef_;
	ArrheniusHBCoeffs_.lookup("SNu") >> SNu_;
	ArrheniusHBCoeffs_.lookup("TRef") >> TRef_;
	ArrheniusHBCoeffs_.lookup("kRef") >> kRef_;
	ArrheniusHBCoeffs_.lookup("Sk") >> Sk_;
	ArrheniusHBCoeffs_.lookup("n") >> n_;
//final da edição

    return true;
}


// ************************************************************************* //
