/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "limitP.H"
#include "fvMesh.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug( limitP, 0 );
        addToRunTimeSelectionTable( option, limitP, dictionary );
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::limitP::writeFileHeader( Ostream& os )
{
    writeHeaderValue( os, "Pmin", Foam::name( Pmin_ ) );
    writeHeaderValue( os, "Pmax", Foam::name( Pmax_ ) );
    writeCommented( os, "Time" );
    writeTabbed( os, "nDampedCellsMin_[count]" );
    writeTabbed( os, "nDampedCellsMin_[%]" );
    writeTabbed( os, "nDampedCellsMax_[count]" );
    writeTabbed( os, "nDampedCellsMax_[%]" );

    os  << endl;

    writtenHeader_ = true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::limitP::limitP
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
    :
    fv::cellSetOption( name, modelType, dict, mesh ),
    writeFile( mesh, name, typeName, dict, false ),
    Pmin_( 0 ),
    Pmax_( 0 ),
    phase_()
{
    if( isActive() )
    {
        read( dict );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::limitP::read( const dictionary& dict )
{
    if( !( fv::cellSetOption::read( dict ) && writeFile::read( dict ) ) )
    {
        return false;
    }

    coeffs_.readEntry( "min", Pmin_ );
    coeffs_.readEntry( "max", Pmax_ );
    coeffs_.readIfPresent( "phase", phase_ );

    if( Pmax_ < Pmin_ )
    {
        FatalIOErrorInFunction( dict )
                << "Minimum pressure limit cannot exceed maximum limit" << nl
                << "min = " << Pmin_ << nl
                << "max = " << Pmax_
                << exit( FatalIOError );
    }

    if( Pmin_ < 0 )
    {
        FatalIOErrorInFunction( dict )
                << "Minimum pressure limit cannot be negative" << nl
                << "min = " << Pmin_
                << exit( FatalIOError );
    }

    // Set the field name to that of the energy
    // field from which the pressure is obtained
     auto& thermo =
        mesh_.lookupObjectRef<basicThermo>
        (
            IOobject::groupName( basicThermo::dictName, phase_ )
        );

    fieldNames_.resize( 1, thermo.p().name() );

    fv::option::resetApplied();


    if( canResetFile() )
    {
        resetFile( typeName );
    }

    if( canWriteHeader() )
    {
        writeFileHeader( file() );
    }


    return true;
}


void Foam::fv::limitP::correct( volScalarField& p1 )
{
    auto& thermo =
        mesh_.lookupObjectRef<basicThermo>
        (
            IOobject::groupName( basicThermo::dictName, phase_ )
        );

    volScalarField& p = thermo.p();

    scalar Pmin0 = min( p ).value();
    scalar Pmax0 = max( p ).value();

    // Count nTotCells ourselves
    // (maybe only applying on a subset)
    label nBelowMin( 0 );
    label nAboveMax( 0 );
    const label nTotCells( returnReduce( cells_.size(), sumOp<label>() ) );

    forAll( cells_, i )
    {
        const label celli = cells_[i];

        if( p[celli] < Pmin_ )
        {
            p[celli] = Pmin_;
            ++nBelowMin;
        }
        else if( p[celli] > Pmax_ )
        {
            p[celli] = Pmax_;
            ++nAboveMax;
        }
    }

    reduce( nBelowMin, sumOp<label>() );
    reduce( nAboveMax, sumOp<label>() );

    reduce( Pmin0, minOp<scalar>() );
    reduce( Pmax0, maxOp<scalar>() );

    // Percent, max 2 decimal places
    const auto percent = []( scalar num, label denom ) -> scalar
    {
        return ( denom ? 1e-2*round( 1e4*num/denom ) : 0 );
    };

    const scalar nBelowMinPercent = percent( nBelowMin, nTotCells );
    const scalar nAboveMaxPercent = percent( nAboveMax, nTotCells );

    Info<< type() << ' ' << name_ << " Lower limited " << nBelowMin << " ("
        << nBelowMinPercent
        << "%) of cells, with min limit " << Pmin_ << endl;

    Info<< type() << ' ' << name_ << " Upper limited " << nAboveMax << " ("
        << nAboveMaxPercent
        << "%) of cells, with max limit " << Pmax_ << endl;

    Info<< type() << ' ' << name_ << " Unlimited Pmin " << Pmin0 << endl;
    Info<< type() << ' ' << name_ << " Unlimited Pmax " << Pmax0 << endl;


    if( canWriteToFile() )
    {
        file()
                << mesh_.time().timeOutputValue() << token::TAB
                << nBelowMin << token::TAB
                << nBelowMinPercent << token::TAB
                << nAboveMax << token::TAB
                << nAboveMaxPercent
                << endl;
    }


    // Handle boundaries in the case of 'all'
    bool changedValues = ( nBelowMin || nAboveMax );

    if( !cellSetOption::useSubMesh() )
    {
        volScalarField::Boundary& bf = p.boundaryFieldRef();

        forAll( bf, patchi )
        {
            fvPatchScalarField& pp = bf[patchi];

            if( !pp.fixesValue() )
            {          
                forAll( pp, facei )
                {
                    if( pp[facei] < Pmin_ )
                    {
                        pp[facei] = Pmin_;
                        changedValues = true;
                    }
                    else if( pp[facei] > Pmax_ )
                    {
                        pp[facei] = Pmax_;
                        changedValues = true;
                    }
                }
            }
        }
    }


    if( returnReduceOr( changedValues ) )
    {
        // We've changed internal values so give
        // boundary conditions opportunity to correct
        p.correctBoundaryConditions();
    }
}


// ************************************************************************* //
