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
    Foam::viscosityModels::turnOnPlasticity2

Description
     Standard power-law non-Newtonian viscosity model.

SourceFiles
    turnOnPlasticity2.C

\*---------------------------------------------------------------------------*/

#ifndef turnOnPlasticity2_H
#define turnOnPlasticity2_H

#include "viscosityModel.H"
#include "dimensionedScalar.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace viscosityModels
{

/*---------------------------------------------------------------------------*\
                           Class turnOnPlasticity2 Declaration
\*---------------------------------------------------------------------------*/

class turnOnPlasticity2
:
    public viscosityModel
{
    // Private data

        dictionary turnOnPlasticity2Coeffs_;

//inicio da edição
		dimensionedScalar n_;
		dimensionedScalar tauYRef_;
		dimensionedScalar STau_;
		dimensionedScalar TYS_;
		dimensionedScalar TPL_;
		dimensionedScalar kRef_;
		dimensionedScalar Sk_;
		dimensionedScalar nu0Ref_;
		dimensionedScalar SNu_;
		dimensionedScalar TRef_;
		dimensionedScalar nuMax_;		


        volScalarField nu_;
        volScalarField tauY_;
        volScalarField k_;
        volScalarField nu0_;

//final da edição


    // Private Member Functions

        //- Calculate and return the laminar viscosity
        tmp<volScalarField> calcNu() const;
        tmp<volScalarField> calcTauY() const;
        tmp<volScalarField> calcK() const;
        tmp<volScalarField> calcNu0() const;

public:

    //- Runtime type information
    TypeName("turnOnPlasticity2");


    // Constructors

        //- Construct from components
        turnOnPlasticity2
        (
            const word& name,
            const dictionary& viscosityProperties,
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    ~turnOnPlasticity2()
    {}


    // Member Functions

        //- Return the laminar viscosity
        tmp<volScalarField> nu() const
        {
            return nu_;
        }
        tmp<volScalarField> tauY() const
        {
            return tauY_;
        }
        
        //- Return the laminar viscosity for patch
        tmp<scalarField> nu(const label patchi) const
        {
            return nu_.boundaryField()[patchi];
        }

        //- Correct the laminar viscosity
        void correct()
        {
            nu_ = calcNu();
            tauY_=calcTauY();
            k_=calcK();
            nu0_=calcNu0();
        }

        //- Read transportProperties dictionary
        bool read(const dictionary& viscosityProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace viscosityModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
