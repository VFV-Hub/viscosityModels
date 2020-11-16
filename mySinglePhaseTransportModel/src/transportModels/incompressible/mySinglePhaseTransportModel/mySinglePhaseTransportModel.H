/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::mySinglePhaseTransportModel

Description
    A simple single-phase transport model based on myViscosityModel.

    Used by the incompressible single-phase solvers like simpleFoam,
    turbFoam etc.

SourceFiles
    mySinglePhaseTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef mySinglePhaseTransportModel_H
#define mySinglePhaseTransportModel_H

#include "incompressible/transportModel/transportModel.H"
#include "IOdictionary.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class myViscosityModel;

/*---------------------------------------------------------------------------*\
                Class mySinglePhaseTransportModel Declaration
\*---------------------------------------------------------------------------*/

class mySinglePhaseTransportModel
:
    public IOdictionary,
    public transportModel
{
    // Private Data

        autoPtr<myViscosityModel> myViscosityModelPtr_;


public:

    //- Runtime type information
    TypeName("mySinglePhaseTransportModel");


    // Constructors

        //- Construct from components
        mySinglePhaseTransportModel
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );

        //- Disallow default bitwise copy construction
        mySinglePhaseTransportModel(const mySinglePhaseTransportModel&);


    //- Destructor
    virtual ~mySinglePhaseTransportModel();


    // Member Functions

        //- Return the laminar viscosity
        virtual tmp<volScalarField> nu() const;
        
        //- Return the laminar viscosity
        virtual tmp<volScalarField> tauY() const;

        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> nu(const label patchi) const;
        
        //- Return the laminar viscosity for patch
        virtual tmp<scalarField> tauY(const label patchi) const;

        //- Correct the laminar viscosity
        virtual void correct();

        //- Read transportProperties dictionary
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const mySinglePhaseTransportModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
