/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    Foam::fv::limitP

Group
    grpFvOptionsCorrections

Description
    Corrects pressure field (i.e. \c p) within a specified
    region by applying limits between a given minimum and maximum.

    Corrections applied to:
    \verbatim
      p    | Temperature                   [Pa]
    \endverbatim

    Required fields:
    \verbatim
      p    | Temperature                   [Pa]
    \endverbatim

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    limitP1
    {
        // Mandatory entries (unmodifiable)
        type            limitP;

        // Mandatory entries (runtime modifiable)
        min             <minValue>;
        max             <maxValue>;

        // Optional entries (runtime modifiable)
        phase           <phaseName>;

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                   | Type   | Reqd  | Dflt
      type      | Type name: limitP   | word   | yes   | -
      min       | Minimum pressure limit [Pa] | scalar | yes   | -
      max       | Maximum pressure limit [Pa] | scalar | yes   | -
      phase     | Name of operand phase field   | word   | no    | ""
    \endtable

    The inherited entries are elaborated in:
      - \link fvOption.H \endlink
      - \link cellSetOption.H \endlink
      - \link writeFile.H \endlink

See also
  - Foam::fv::fixedTemperatureConstraint

SourceFiles
    limitP.C

\*---------------------------------------------------------------------------*/

#ifndef fv_limitP_H
#define fv_limitP_H

#include "cellSetOption.H"
#include "writeFile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {

        /*---------------------------------------------------------------------------*\
                            Class limitP Declaration
        \*---------------------------------------------------------------------------*/

        class limitP
            :
            public fv::cellSetOption,
            public functionObjects::writeFile
        {
        protected:

            // Protected Data

            //- Minimum pressure limit [K]
            scalar Pmin_;

            //- Maximum pressure limit [K]
            scalar Pmax_;

            //- Optional phase name
            word phase_;


            // Protected Member Functions

            //- Write file header information
            void writeFileHeader( Ostream& os );


        public:

            //- Runtime type information
            TypeName( "limitP" );


            // Constructors

            //- Construct from components
            limitP
            (
                const word& name,
                const word& modelType,
                const dictionary& dict,
                const fvMesh& mesh
            );

            //- No copy construct
            limitP( const limitP& ) = delete;

            //- No copy assignment
            void operator=( const limitP& ) = delete;


            //- Destructor
            virtual ~limitP() = default;


            // Member Functions

            //- Read dictionary
            virtual bool read( const dictionary& dict );

            //- Correct the energy field
            virtual void correct( volScalarField& p1 );
        };


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    } // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
