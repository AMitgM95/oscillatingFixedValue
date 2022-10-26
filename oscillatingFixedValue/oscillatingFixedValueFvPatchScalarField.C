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

\*---------------------------------------------------------------------------*/

#include "oscillatingFixedValueFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::oscillatingFixedValueFvPatchScalarField::t() const
{
    return db().time().timeOutputValue();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oscillatingFixedValueFvPatchScalarField::
oscillatingFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    amplitude_(),
    frequency_(),
    offset_(0.0),
    refValue_(p.size(), Zero),
    curTimeIndex_(-1)
{
}


Foam::oscillatingFixedValueFvPatchScalarField::
oscillatingFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    amplitude_(Function1<scalar>::New("amplitude", dict)),
    frequency_(Function1<scalar>::New("frequency", dict)),
    offset_(dict.lookup<scalar>("offset")),
    refValue_("refValue", dict, p.size()),
    curTimeIndex_(-1)
{


    fixedValueFvPatchScalarField::evaluate();

    /*
    // Initialise with the value entry if evaluation is not possible
    fvPatchScalarField::operator=
    (
        scalarField("value", dict, p.size())
    );
    */
}


Foam::oscillatingFixedValueFvPatchScalarField::
oscillatingFixedValueFvPatchScalarField
(
    const oscillatingFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    amplitude_(ptf.amplitude_, false),
    frequency_(ptf.frequency_, false),
    offset_(ptf.offset_),
    refValue_(mapper(ptf.refValue_)),
    curTimeIndex_(-1)
{}


Foam::oscillatingFixedValueFvPatchScalarField::
oscillatingFixedValueFvPatchScalarField
(
    const oscillatingFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    amplitude_(ptf.amplitude_, false),
    frequency_(ptf.frequency_, false),
    offset_(ptf.offset_),
    refValue_(ptf.refValue_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::oscillatingFixedValueFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    m(refValue_, refValue_);
}


void Foam::oscillatingFixedValueFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const oscillatingFixedValueFvPatchScalarField& tiptf =
        refCast<const oscillatingFixedValueFvPatchScalarField>(ptf);

    refValue_.rmap(tiptf.refValue_, addr);
}


void Foam::oscillatingFixedValueFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
    fixedValueFvPatchScalarField::operator==
    (
        (1.0
      + amplitude_->value(t())*sin(constant::mathematical::twoPi*frequency_->value(t())*t()))*refValue_
      + offset_
    );
    
        curTimeIndex_ = this->db().time().timeIndex();
    }


    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::oscillatingFixedValueFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, amplitude_());
    writeEntry(os, frequency_());
    writeEntry(os, "offset", offset_);
    writeEntry(os, "refValue", refValue_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        oscillatingFixedValueFvPatchScalarField
    );
}

// ************************************************************************* //
