/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "mySinglePhaseTransportModel.H"
#include "myViscosityModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mySinglePhaseTransportModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mySinglePhaseTransportModel::mySinglePhaseTransportModel
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    myViscosityModelPtr_(myViscosityModel::New("nu","tauY", *this, U, phi))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mySinglePhaseTransportModel::~mySinglePhaseTransportModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::mySinglePhaseTransportModel::nu() const
{
    return myViscosityModelPtr_->nu();
}

Foam::tmp<Foam::volScalarField>
Foam::mySinglePhaseTransportModel::tauY() const
{
    return myViscosityModelPtr_->tauY();
}


Foam::tmp<Foam::scalarField>
Foam::mySinglePhaseTransportModel::nu(const label patchi) const
{
    return myViscosityModelPtr_->nu(patchi);
}

Foam::tmp<Foam::scalarField>
Foam::mySinglePhaseTransportModel::tauY(const label patchi) const
{
    return myViscosityModelPtr_->tauY(patchi);
}


void Foam::mySinglePhaseTransportModel::correct()
{
    myViscosityModelPtr_->correct();
}


bool Foam::mySinglePhaseTransportModel::read()
{
    if (regIOobject::read())
    {
        return myViscosityModelPtr_->read(*this);
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
