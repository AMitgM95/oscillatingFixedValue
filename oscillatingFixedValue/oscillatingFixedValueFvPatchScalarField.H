/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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
    Foam::oscillatingFixedValueFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides an oscillating condition in terms of
    amplitude and frequency.

        \f[
            x_p = (1 + a sin(\pi f t))x_{ref} + x_o
        \f]

    where

    \vartable
        x_p     | patch values
        x_{ref} | patch reference values
        x_o     | patch offset values
        a       | amplitude
        f       | frequency [1/s]
        t       | time [s]
    \endvartable

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        refValue     | reference value         | yes         |
        offset       | offset value            | no          | 0.0
        amplitude    | oscillation amplitude   | yes         |
        frequency    | oscillation frequency   | yes         |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            oscillatingFixedValue;
        refValue        uniform 5.0;
        offset          0.0;
        amplitude       constant 0.5;
        frequency       constant 10;
    }
    \endverbatim

Note
    The amplitude and frequency entries are DataEntry types, able to describe
    time varying functions.  The example above gives the usage for supplying
    constant values.

SeeAlso
    Foam::DataEntry

SourceFiles
    oscillatingFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef BC_H
#define BC_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class oscillatingFixedValueFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class oscillatingFixedValueFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private Data
    
        //- Amplitude
        autoPtr<Function1<scalar>> amplitude_;

        //- Frequency
        autoPtr<Function1<scalar>> frequency_;
        
        //- Offset
        scalar offset_;

        //- Reference value
        scalarField refValue_;

        //- Label, e.g. patch index, current time index
        label curTimeIndex_;


    // Private Member Functions

        //- Return current time
        scalar t() const;


public:

    //- Runtime type information
    TypeName("oscillatingFixedValue");


    // Constructors

        //- Construct from patch and internal field
        oscillatingFixedValueFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        oscillatingFixedValueFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedValueTypeFvPatchField
        //  onto a new patch
        oscillatingFixedValueFvPatchScalarField
        (
            const oscillatingFixedValueFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        oscillatingFixedValueFvPatchScalarField
        (
            const oscillatingFixedValueFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        oscillatingFixedValueFvPatchScalarField
        (
            const oscillatingFixedValueFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new oscillatingFixedValueFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchScalarField&, const labelList&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
